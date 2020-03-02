#include "macros.h"
#:include "macros.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module gasci
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
        near_zero
    use sort_mod, only : sort
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles, excit_gen_store_type
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        last_tgt_unknown, set_last_tgt
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    implicit none

    private
    public :: is_valid, is_connected, GAS_specification, GAS_specification_t, &
        init_GAS, clear_GAS, get_nGAS, &
        generate_nGAS_excitation, contains_nI

    public :: iGAS_from_spatorb, iGAS_from_spinorb

    !> Speficies the GAS spaces.
    !> It is assumed, that the GAS spaces are contigous.
    !> The indices are:
    !>  n_orbs_per_GAS(1:nGAS), n_min(1:nGAS), n_max(1:nGAS)
    !> n_orbs_per_GAS(iGAS) specifies how many orbitals are in
    !> the `iGAS` GAS space in the `iRep` Irrep.
    !> n_min(iGAS) specifies the cumulated! minimum particle number per GAS space.
    !> n_max(iGAS) specifies the cumulated! maximum particle number per GAS space.
    type :: GAS_specification_t
        integer, allocatable :: n_orbs(:), n_min(:), n_max(:)
    end type

    type(GAS_specification_t) :: GAS_specification

contains

    pure integer function get_nGAS(GAS_spec)
        type(GAS_specification_t), intent(in) :: GAS_spec
        get_nGAS = size(GAS_spec%n_orbs)
    end function

    logical pure function is_valid(GAS_spec, n_particles)
        type(GAS_specification_t), intent(in) :: GAS_spec
        integer, intent(in), optional :: n_particles

        logical :: shapes_match, nEl_correct
        integer :: nGAS

        nGAS = get_nGAS(GAS_spec)
        shapes_match = size(GAS_spec%n_min) == nGAS &
                       .and. size(GAS_spec%n_max) == nGAS &
                       .and. size(GAS_spec%n_orbs) == nGAS

        if (present(n_particles)) then
            nEl_correct = GAS_spec%n_min(nGAS)  == n_particles &
                          .and. GAS_spec%n_max(nGAS) == n_particles
        else
            nEl_correct = .true.
        end if

        ! TODO(Oskar): Switch back on
        ! nOrbs_correct = 2 * GAS_specification%n_orbs(nGAS) == nBasis

        is_valid = shapes_match .and. nEl_correct
    end function

    pure function is_connected(GAS_spec) result(res)
        type(GAS_specification_t), intent(in) :: GAS_spec
        logical :: res
        res = all(GAS_spec%n_min(:) == GAS_spec%n_max(:))
    end function

    subroutine init_GAS()
        character(*), parameter :: this_routine = 'init_GAS'
        write(*, *) this_routine, 'called'
    end subroutine

    subroutine clear_GAS()
        character(*), parameter :: this_routine = 'clear_GAS'
        write(*, *) this_routine, 'called'
    end subroutine

    !> Could be changed into lookup variable
    elemental function iGAS_from_spatorb(GAS_spec, orbidx) result(res)
        type(GAS_specification_t), intent(in) :: GAS_spec
        integer, intent(in) :: orbidx
        integer :: res

        integer :: iGAS

        ! We assume that orbidx <= GAS_specification%n_orbs(nGAS)
        do iGAS = 1, get_nGAS(GAS_spec)
            if (orbidx <= GAS_spec%n_orbs(iGAS)) then
                res = iGAS
                return
            end if
        end do
        res = -1
    end function

    elemental function iGAS_from_spinorb(GAS_spec, orbidx) result(res)
        type(GAS_specification_t), intent(in) :: GAS_spec
        integer, intent(in) :: orbidx
        integer :: res

        res = iGAS_from_spatorb(GAS_spec, (orbidx + 1) .div. 2)
    end function


    function contains_nI(GAS_spec, nI) result(res)
        type(GAS_specification_t), intent(in) :: GAS_spec
        integer, intent(in) :: nI(:)
        logical :: res

        !> Returns iGAS for the i-th electron.
        integer :: GAS_table(size(nI))
        integer :: n_particle, iGAS

        res = .true.
        GAS_table = iGAS_from_spinorb(GAS_spec, nI)
        n_particle = 0
        do iGAS = 1, get_nGAS(GAS_spec)
            ! Calculate cumulated particles in iGAS
            n_particle = n_particle + count(GAS_table == iGAS)
            if (n_particle < GAS_spec%n_min(iGAS) &
                .or. GAS_spec%n_max(iGAS) < n_particle) then
                res = .false.
                return
            end if
        end do
    end function


    subroutine generate_nGAS_excitation(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex, tParity, pGen, hel, store, part_type)
        ! particle number conserving GAS excitation generator:
        ! we create only excitations, that preserver the number of electrons within
        ! each active space

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        real(dp) :: r

        @:unused_var(exFlag, part_type, store)
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        ! ! single or double excitation?
        ! r = genrand_real2_dSFMT()
        ! if (r < pDoubles) then
        !     call generate_nGAS_double(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
        !     pgen = pgen * pDoubles
        !     ic = 2
        ! else
        !     call generate_nGAS_single(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
        !     pgen = pgen * (1.0_dp - pDoubles)
        !     ic = 1
        ! end if
    end subroutine generate_nGAS_excitation


end module gasci
