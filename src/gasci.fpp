#include "macros.h"
#:include "macros.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']
#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']
#:set excitation_t = 'ASDF'

module gasci
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
        near_zero, cumsum
    use sort_mod, only : sort
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles, excit_gen_store_type
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        last_tgt_unknown, set_last_tgt
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    implicit none

    private
    public :: is_valid, is_connected, GAS_specification, GASSpec_t, &
        init_GAS, clear_GAS, get_nGAS, &
!         generate_nGAS_excitation, &
        contains_det, particles_per_GAS, &
        get_possible_spaces

    public :: get_iGAS, operator(.contains.)


    integer, parameter :: idx_alpha = 1, idx_beta = 2

    !> Speficies the GAS spaces.
    !> It is assumed, that the GAS spaces are contigous.
    !> The indices are:
    !>  n_orbs_per_GAS(1:nGAS), n_min(1:nGAS), n_max(1:nGAS)
    !> n_orbs_per_GAS(iGAS) specifies how many orbitals are in
    !> the `iGAS` GAS space in the `iRep` Irrep.
    !> n_min(iGAS) specifies the cumulated! minimum particle number per GAS space.
    !> n_max(iGAS) specifies the cumulated! maximum particle number per GAS space.
    type :: GASSpec_t
        integer, allocatable :: n_orbs(:), n_min(:), n_max(:)
    end type

    interface GASSpec_t
        module procedure construction_from_array_GAS_specification_t
    end interface

    type(GASSpec_t) :: GAS_specification

    !>  @brief
    !>      Return the GAS spaces, where one particle can be created.
    !>
    !>  @details
    !>  It can be proven, that the spaces where a particle can
    !>  be created are contigous, so a two element integer array
    !>  [lower_bound, upper_bound] is returned.
    !>  As long as lower_bound <= iGAS <= upper_bound is true,
    !>  the created particle will lead to a valid Slater-Determinant.
    !>
    !>  It is **assumed** that the input determinant is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  Checks are only performed in DEBUG compilation mode and
    !>  the return value is undefined, if this is not the case!
    !>
    !>  It is possible to delete additional particles with the
    !>  optional argument `additional_holes` before checking the
    !>  validity of particle creation.
    !>  It is assumed, that they are occupied in det_I!
    !>  Checks are only performed in DEBUG compilation mode and
    !>  the return value is undefined, if this is not the case!
    !>
    !>  If more than one particle should be created, the optional argument
    !>  n_particles (default 1) should be used.
    !>  Note, that this function returns possible spaces where the
    !>  first of the n_particles can be created.
    !>  After modifying `det_I`, or `additional_holes` the function
    !>  has to be called again with `n_particles - 1`.
    !>
    !>  If no creation is allowed by the GAS constraints,
    !>  the bounds will be returned as zero.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] det_I, An index of occupied spatial
    !>      or spin orbitals (SpinOrbIdx_t, SpatOrbIdx_t).
    !>  @param[in] additional_holes, optional, An index of orbitals
    !>      where particles should be deleted.
    !>      It is assumed, that they are occupied in det_I.
    !>      (SpinOrbIdx_t, SpatOrbIdx_t).
    !>  @param[in] n_particles, optional, The total number of particles
    !>      that will be created. Defaults to one (integer).
    interface get_possible_spaces
    #:for orb_idx_type in OrbIdxTypes
        module procedure get_possible_spaces_${orb_idx_type}$
    #:endfor
    end interface

    !>  @brief
    !>      Query wether a determinant is contained in the GAS space.
    !>
    !>  @details
    !>  It is **assumed** that the determinant is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  The return value is not defined, if that is not the case!
    !>
    !>  An operator for this function is also implemented that allows
    !>  to write::
    !>
    !>  GAS_spec .contains. det_I
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] occupied, An index of occupied spatial
    !>      or spin orbitals (SpinOrbIdx_t, SpatOrbIdx_t).
    interface contains_det
        #:for orb_idx_type in OrbIdxTypes
            module procedure contains_det_${orb_idx_type}$
        #:endfor
    end interface


    !>  @brief
    !>      Get the number of particles per GAS space.
    !>
    !>  @details
    !>  This generic function works for spin- and spatial orbitals
    !>  and requires a GAS specification.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] occupied, An index of occupied spatial
    !>      or spin orbitals (SpinOrbIdx_t, SpatOrbIdx_t).
    interface particles_per_GAS
        #:for orb_idx_type in OrbIdxTypes
            module procedure particles_per_GAS_${orb_idx_type}$
        #:endfor
    end interface


    !>  @brief
    !>      Get the GAS space for each orbital in idx.
    !>
    !>  @details
    !>  This generic function works for spin- and spatial orbitals
    !>  and requires a GAS specification.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] idx, An index of either spatial, or a spin orbitals.
    interface get_iGAS
        module procedure get_iGAS_SpatOrbIdx_t, get_iGAS_SpinOrbIdx_t
    end interface

    !>  @brief
    !>      Operator version of contains_det
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] occupied, An index of either spatial
    !>      or spin orbitals (SpinOrbIdx_t, SpatOrbIdx_t).
    interface operator(.contains.)
    #:for orb_idx_type in OrbIdxTypes
        module procedure contains_det_${orb_idx_type}$
    #:endfor
    end interface

contains

    pure function construction_from_array_GAS_specification_t(n_orbs, n_min, n_max) result(GAS_spec)
        integer, intent(in) :: n_orbs(:), n_min(:), n_max(:)
        type(GASSpec_t) :: GAS_spec
        GAS_spec%n_orbs = n_orbs(:)
        GAS_spec%n_min = n_min(:)
        GAS_spec%n_max = n_max(:)
    end function

    pure integer function get_nGAS(GAS_spec)
        type(GASSpec_t), intent(in) :: GAS_spec
        get_nGAS = size(GAS_spec%n_orbs)
    end function

    logical pure function is_valid(GAS_spec, n_particles, n_basis)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in), optional :: n_particles, n_basis

        logical :: shapes_match, nEl_correct, pauli_principle, monotonic, &
            n_orbs_correct
        integer :: nGAS, iGAS

        nGAS = get_nGAS(GAS_spec)

        associate(n_orbs => GAS_spec%n_orbs, n_min => GAS_spec%n_min, &
                  n_max => GAS_spec%n_max)

            shapes_match = &
                all([size(n_orbs), size(n_min), size(n_max)] == nGAS)

            if (present(n_particles)) then
                nEl_correct = all([n_min(nGAS), n_max(nGAS)] == n_particles)
            else
                nEl_correct = n_min(nGAS) == n_max(nGAS)
            end if

            pauli_principle = all(n_min(:) <= n_orbs(:) * 2)

            if (nGAS >= 2) then
                monotonic = all([all(n_orbs(2:) >= n_orbs(: nGAS - 1)), &
                                 all(n_max(2:) >= n_max(: nGAS - 1)), &
                                 all(n_min(2:) >= n_min(: nGAS - 1))])
            else
                monotonic = .true.
            end if

            if (present(n_basis)) then
                n_orbs_correct = 2 * n_orbs(nGAS) == n_basis
            else
                n_orbs_correct = .true.
            end if
        end associate

        is_valid = all([shapes_match, nEl_correct, pauli_principle, &
                        monotonic, n_orbs_correct])
    end function

    pure function is_connected(GAS_spec) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
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
    pure function get_iGAS_SpatOrbIdx_t(GAS_spec, idx) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpatOrbIdx_t), intent(in) :: idx
        integer :: res(size(idx%idx))
        res = f(idx%idx)
        contains
            elemental function f(orbidx) result(iGAS)
                integer, intent(in) :: orbidx
                integer :: iGAS
                ! We assume that orbidx <= GAS_specification%n_orbs(nGAS)
                do iGAS = 1, get_nGAS(GAS_spec)
                    if (orbidx <= GAS_spec%n_orbs(iGAS)) return
                end do
                iGAS =  -1
            end function
    end function

    pure function get_iGAS_SpinOrbIdx_t(GAS_spec, idx) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: idx
        integer :: res(size(idx%idx))
        res = f(idx%idx)
        contains
            elemental function f(orbidx) result(iGAS)
                integer, intent(in) :: orbidx
                integer :: iGAS
                ! We assume that orbidx <= GAS_specification%n_orbs(nGAS)
                do iGAS = 1, get_nGAS(GAS_spec)
                    if (((orbidx + 1) .div. 2) <= GAS_spec%n_orbs(iGAS)) return
                end do
                iGAS =  -1
            end function
    end function


#:for orb_idx_type in OrbIdxTypes
    pure function particles_per_GAS_${orb_idx_type}$(GAS_spec, occupied) result(n_particle)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        integer :: n_particle(get_nGAS(GAS_spec))
        !> Returns iGAS for the i-th electron.
        integer :: GAS_table(size(occupied%idx)), iGAS
        GAS_table(:) = get_iGAS(GAS_spec, occupied)
        n_particle(:) = [(count(GAS_table(:) == iGAS), iGAS = 1, get_nGAS(GAS_spec))]
    end function
#:endfor


#:for orb_idx_type in OrbIdxTypes
    pure function contains_det_${orb_idx_type}$(GAS_spec, occupied) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        logical :: res

        !> Cumulated number of particles per iGAS
        integer :: cum_n_particle(get_nGAS(GAS_spec))
        cum_n_particle(:) = cumsum(particles_per_GAS(GAS_spec, occupied))
        res = all(GAS_spec%n_min(:) <= cum_n_particle(:) &
            .and. cum_n_particle(:) <= GAS_spec%n_max(:))
    end function
#:endfor


#:for orb_idx_type in OrbIdxTypes
    DEBUG_IMPURE function get_possible_spaces_${orb_idx_type}$(GAS_spec, det_I, additional_holes, n_particles) result(spaces)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: det_I
        type(${orb_idx_type}$), intent(in), optional :: additional_holes
        integer, intent(in), optional :: n_particles
        character(*), parameter :: this_routine = 'get_possible_spaces_${orb_idx_type}$'

        !> Lower and upper bound for spaces where a particle can be created.
        !> If no particle can be created, then spaces == 0 .
        integer :: spaces(2), n_particles_, i, iGAS

        integer :: &
        !> Cumulated number of particles per iGAS
            cum_n_particle(get_nGAS(GAS_spec)), &
        !> Cumulated deficit per iGAS
            deficit(get_nGAS(GAS_spec)), &
        !> Cumulated vacant orbitals per iGAS
            vacant(get_nGAS(GAS_spec))

        @:def_default(n_particles_, n_particles, 1)

        associate(lower_bound => spaces(1), upper_bound => spaces(2))

            if (present(additional_holes)) then
                ! Check if all additional_holes are occupied in det_I
                do i = 1, size(additional_holes%idx)
                    ASSERT(any(additional_holes%idx(i) == det_I%idx(:)))
                end do
                cum_n_particle(:) = cumsum( &
                    particles_per_GAS(GAS_spec, det_I) &
                    - particles_per_GAS(GAS_spec, additional_holes))
            else
                cum_n_particle(:) = cumsum(particles_per_GAS(GAS_spec, det_I))
            end if

            deficit = GAS_spec%n_min(:) - cum_n_particle(:)
            vacant = GAS_spec%n_max(:) - cum_n_particle(:)
            if (any(n_particles_ < deficit) .or. all(vacant < n_particles_)) then
                spaces(:) = 0
                return
            end if

            ! Find the first index, where a particle has to be created.
            do iGAS = 1, get_nGAS(GAS_spec)
                if (deficit(iGAS) == n_particles_) exit
            end do
            upper_bound = iGAS

            ! We assume that at it is possible to create a particle at least in
            ! the last GAS space.
            ! Search from behind the first occurence where it is not possible
            ! anymore to create a particle.
            ! The lower bound is one GAS index above.
            do iGAS = get_nGAS(GAS_spec), 1, -1
                if (vacant(iGAS) <= 0) exit
            end do
            lower_bound = iGAS + 1

            if (lower_bound > upper_bound) spaces(:) = 0
        end associate
    end function
#:endfor


!     subroutine generate_nGAS_excitation(nI, ilutI, nJ, ilutJ, exFlag, ic, &
!                                         ex, tParity, pGen, hel, store, part_type)
!         ! particle number conserving GAS excitation generator:
!         ! we create only excitations, that preserver the number of electrons within
!         ! each active space
!
!         integer, intent(in) :: nI(nel), exFlag
!         integer(n_int), intent(in) :: ilutI(0:NIfTot)
!         integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
!         integer(n_int), intent(out) :: ilutJ(0:NifTot)
!         real(dp), intent(out) :: pGen
!         logical, intent(out) :: tParity
!         HElement_t(dp), intent(out) :: hel
!         type(excit_gen_store_type), intent(inout), target :: store
!         integer, intent(in), optional :: part_type
!
!         real(dp) :: r
!         type(SpinOrbIdx_t) :: det_J
!         class(Excitation_t), allocatable :: exc
!
! ! #ifdef WARNING_WORKAROUND_
! !         hel = 0.0_dp
! ! #endif
! !         ! single or double excitation?
! !         r = genrand_real2_dSFMT()
! !         if (r >= pDoubles) then
! !             exc = SingleExc_t()
! !         else
! !             exc = DoubleExc_t()
! !         end if
! !
! !         select type(exc)
! !         type is(SingleExc_t)
! !             call gen_exc(SpinOrbIdx_t(nI), exc, det_J, par, pgen)
! !             pgen = pgen * (1.0_dp - pDoubles)
! !             ic = 1
! !         type is(DoubleExc_t)
! !             call gen_exc(SpinOrbIdx_t(nI), exc, det_J, par, pgen)
! !             pgen = pgen * pDoubles
! !             ic = 2
! !         end select
!     end subroutine generate_nGAS_excitation

!     subroutine gen_exc_SingleExc_t(det_I, exc, det_J, par, pgen)
!         type(SpinOrbIdx_t), intent(in) :: det_I
!         type(SingleExc_t), intent(out) :: exc
!         type(SpinOrbIdx_t), intent(out) :: det_J
!         logical, intent(out) :: par
!         real(dp), intent(out) :: pgen
!
!         integer :: elec
!
!         ! adjust pgen
!         pgen = 1.0_dp / real(nEl, kind=dp)
!
!         ! Pick any random electron
!         elec = int(genrand_real2_dSFMT() * nEl) + 1
!         src = det_I%idx(elec)
!         exc%val(1) = src
!
!         ! Get a hole with the same spin projection
!         ! while fullfilling GAS-constraints.
!         spin_idx = get_spin(src)
!         tgt = pick_weighted_hole(det_I, exc, spin_idx, pgen)
!     end subroutine
!
!
!     function pick_weighted_hole_SingleExc_t(det_I, exc, spin_idx, pgen) result(tgt)
!         ! pick a hole of nI with spin ms from the active space with index
!         ! srcGASInd the random number is to be supplied as r
!         ! nI is the source determinant, nJBase the one from which we obtain
!         ! the ket of the matrix element by single excitation
!         type(SpinOrbIdx_t), intent(in) :: det_I
!         type(SingleExc_t), intent(in) :: exc
!         integer, intent(in) :: spin_idx
!         real(dp), intent(inout) :: pgen
!         character(*), parameter :: this_routine = 'pick_weighted_hole_${excitation_t}$'
!
!         integer :: tgt, nOrbs, GAS_list(GAS_size(iGAS))
!         real(dp) :: r, cSum(GAS_size(iGAS))
!
!
!         ASSERT(last_tgt_unknown(exc))
!         ! initialize auxiliary variables
! !         nOrbs = GAS_size(iGAS)
! !         GAS_list = GAS_spin_orb_list(1:nOrbs, iGAS, spin_idx)
!         possible_holes = get_possible_holes(GAS_spec, det_I)
!
!         ASSERT(last_tgt_unknown(exc))
!         ! build the cumulative list of matrix elements <src|H|tgt>
!         cSum = get_cumulative_list(det_I, exc, possible_holes)
!
!         ! now, pick with the weight from the cumulative list
!         r = genrand_real2_dSFMT() * cSum(nOrbs)
!
!         ! there might not be such an excitation
!         if (cSum(nOrbs) > 0) then
!             ! find the index of the target orbital in the gasList
!             tgt = binary_search_first_ge(cSum, r)
!
!             ! adjust pgen with the probability for picking tgt from the cumulative list
!             if (tgt == 1) then
!                 pgen = pgen * cSum(1)
!             else
!                 pgen = pgen * (cSum(tgt) - cSum(tgt - 1))
!             end if
!
!             ! convert to global orbital index
!             tgt = GAS_list(tgt)
!         else
!             tgt = 0
!         end if
!     end function pick_weighted_hole_SingleExc_t
!
!     pure function get_possible_holes(GAS_spec, det_I) result(
!         type(GAS_specification_t), intent(in) :: GAS_spec
!         type(SpatOrbIdx_t), intent(in) :: det_I
!
!         integer ::
!
!
!
!     end function
!
!
!     function get_cumulative_list_SingleExc_t(GAS_list, nI, incomplete_exc) result(cSum)
!         integer, intent(in) :: GAS_list(:), nI(nel)
!         type(${excitation_t}$), intent(in) :: incomplete_exc
!         real(dp) :: cSum(size(GAS_list))
!         character(*), parameter :: this_routine = 'get_cumulative_list_${excitation_t}$'
!
!         real(dp) :: previous
!         type(${excitation_t}$) :: exc
!         integer :: i, nOrbs
!
!         nOrbs = size(GAS_list)
!
!         exc = incomplete_exc
!         ASSERT(last_tgt_unknown(exc))
!
!         ! build the cumulative list of matrix elements <src|H|tgt>
!         previous = 0.0_dp
!         do i = 1, nOrbs
!             call set_last_tgt(exc, GAS_list(i))
!             cSum(i) = get_mat_element(nI, exc) + previous
!             previous = cSum(i)
!         end do
!
!         ! Normalize
!         if (near_zero(cSum(nOrbs))) then
!             cSum(:) = 0.0_dp
!         else
!             cSum(:) = cSum(:) / cSum(nOrbs)
!         end if
!     end function get_cumulative_list_${excitation_t}$
!
!
!     function get_mat_element_SingleExc_t(det_I, exc) result(res)
!         type(SpinOrbIdx_t), intent(in) :: det_I
!         type(SingleExc_t), intent(in) :: exc
!         real(dp) :: res
!
!         associate (tgt => exc%val(2))
!             if (all(det_I%idx /= tgt)) then
!                 res = abs(sltcnd_excit(det_I, exc, .false.))
!             else
!                 res = 0.0_dp
!             end if
!         end associate
!     end function get_mat_element_SingleExc_t

end module gasci
