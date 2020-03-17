#include "macros.h"
#:include "macros.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']
#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']
#:set excitation_t = 'ASDF'

module gasci
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
        near_zero, cumsum, operator(.isclose.)
    use sort_mod, only: sort
    use sets_mod, only: is_sorted, complement, union, disjoint
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles, excit_gen_store_type
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        last_tgt_unknown, set_last_tgt, excite, defined
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, Spin_t, size, calc_spin, &
        calc_spin_raw, spin_values
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    implicit none

    private
    public :: is_valid, is_connected, GAS_specification, GASSpec_t, &
        init_GAS, clear_GAS, get_nGAS, &
        generate_nGAS_excitation, &
        contains_det, particles_per_GAS, &
        get_possible_spaces, get_possible_holes, split_per_GAS

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

    integer, parameter :: EMPTY_BOUNDS = -1

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
    !>  optional argument `add_holes` before checking the
    !>  validity of particle creation.
    !>  It is assumed, that they are occupied in det_I!
    !>  Checks are only performed in DEBUG compilation mode and
    !>  the return value is undefined, if this is not the case!
    !>
    !>  If more than one particle should be created, the optional argument
    !>  n_particles (default 1) should be used.
    !>  Note, that this function returns possible spaces where the
    !>  first of the n_particles can be created.
    !>  After modifying `det_I`, or `add_holes` the function
    !>  has to be called again with `n_particles - 1`.
    !>
    !>  If no creation is allowed by the GAS constraints,
    !>  the bounds will be returned as integer constant EMPTY_BOUNDS.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] det_I, An index of occupied spatial
    !>      or spin orbitals (SpinOrbIdx_t, SpatOrbIdx_t).
    !>  @param[in] add_holes, optional, An index of orbitals
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


    interface get_cumulative_list
        #:for excitation_t in ExcitationTypes
            module procedure get_cumulative_list_${excitation_t}$
        #:endfor
    end interface

    interface split_per_GAS
        #:for orb_idx_type in OrbIdxTypes
            module procedure split_per_GAS_${orb_idx_type}$
        #:endfor
    end interface

    interface
        subroutine stop_all(sub_name, error_msg)
            character(*), intent(in) :: sub_name, error_msg
        end subroutine
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
        integer :: res(size(idx))
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
        integer :: res(size(idx))
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
    pure function particles_per_GAS_${orb_idx_type}$(splitted_det_I) result(n_particle)
        type(${orb_idx_type}$), intent(in) :: splitted_det_I(:)
        integer :: n_particle(size(splitted_det_I))
        integer :: i
        n_particle = [(size(splitted_det_I(i)), i = 1, size(splitted_det_I))]
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    function split_per_GAS_${orb_idx_type}$(GAS_spec, occupied) result(splitted_orbitals)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        type(${orb_idx_type}$), allocatable :: splitted_orbitals(:)

        integer :: iGAS, prev, i, GAS_table(size(occupied))

        allocate(splitted_orbitals(get_nGAS(GAS_spec)))
        GAS_table = get_iGAS(GAS_spec, occupied)

        ! We assume that GAS_table is sorted and looks like:
        ! [1, 1, ..., 2, 2, ..., nGAS, nGAS]
        i = 0
        do iGAS = 1, get_nGAS(GAS_spec)
            prev = i
            if (i < size(GAS_table)) then
                do while (GAS_table(i + 1) == iGAS)
                    i = i + 1
                    if (i == size(GAS_table)) exit
                end do
            end if
            if (prev + 1 <= size(occupied%idx) .and. prev + 1 <= i) then
                splitted_orbitals(iGAS)%idx = occupied%idx(prev + 1 : i)
            else
                splitted_orbitals(iGAS)%idx = [integer::]
            end if
        end do
    end function
#:endfor


#:for orb_idx_type in OrbIdxTypes
    function contains_det_${orb_idx_type}$(GAS_spec, occupied) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        logical :: res

        !> Cumulated number of particles per iGAS
        integer :: cum_n_particle(get_nGAS(GAS_spec)), i

        associate(splitted => split_per_GAS(GAS_spec, occupied))
            cum_n_particle = cumsum([(size(splitted(i)), i = 1, get_nGAS(GAS_spec))])
        end associate

        res = all(GAS_spec%n_min(:) <= cum_n_particle(:) &
            .and. cum_n_particle(:) <= GAS_spec%n_max(:))
    end function
#:endfor


#:for orb_idx_type in OrbIdxTypes
    function get_possible_spaces_${orb_idx_type}$(GAS_spec, splitted_det_I, add_holes, add_particles, n_total) result(spaces)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: splitted_det_I(:)
        type(${orb_idx_type}$), intent(in), optional :: add_holes, add_particles
        integer, intent(in), optional :: n_total
        character(*), parameter :: this_routine = 'get_possible_spaces_${orb_idx_type}$'
        integer, allocatable :: spaces(:)

        !> Lower and upper bound for spaces where a particle can be created.
        !> If no particle can be created, then spaces == 0 .
        integer :: n_total_, i, iGAS, lower_bound, upper_bound

        integer :: &
        !> Cumulated number of particles per iGAS
            cum_n_particle(get_nGAS(GAS_spec)), &
        !> Cumulated deficit per iGAS
            deficit(get_nGAS(GAS_spec)), &
        !> Cumulated vacant orbitals per iGAS
            vacant(get_nGAS(GAS_spec))

        @:def_default(n_total_, n_total, 1)

        associate(A => particles_per_GAS(splitted_det_I))
            if (present(add_holes) .and. present(add_particles)) then
                associate(B => particles_per_GAS(split_per_GAS(GAS_spec, add_holes)), &
                          C => particles_per_GAS(split_per_GAS(GAS_spec, add_particles)))
                    cum_n_particle = cumsum(A - B + C)
                end associate
            else if (present(add_holes)) then
                associate(B => particles_per_GAS(split_per_GAS(GAS_spec, add_holes)))
                    cum_n_particle = cumsum(A - B)
                end associate
            else if (present(add_particles)) then
                associate(C => particles_per_GAS(split_per_GAS(GAS_spec, add_particles)))
                    cum_n_particle = cumsum(A + C)
                end associate
            else
                cum_n_particle = cumsum(A)
            end if
        end associate

        deficit = GAS_spec%n_min(:) - cum_n_particle(:)
        vacant = GAS_spec%n_max(:) - cum_n_particle(:)
        if (any(n_total_ < deficit) .or. all(vacant < n_total_)) then
            spaces = [integer::]
            return
        end if

        ! Find the first index, where a particle has to be created.
        do iGAS = 1, get_nGAS(GAS_spec)
            if (deficit(iGAS) == n_total_) exit
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

        if (lower_bound > upper_bound .or. lower_bound > get_nGAS(GAS_spec)) then
            spaces = [integer::]
        else
            spaces = [lower_bound, upper_bound]
        end if
    end function
#:endfor


    function get_possible_holes(GAS_spec, det_I, add_holes, add_particles, n_total, m_s) result(possible_holes)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SpinOrbIdx_t), intent(in), optional :: add_holes
        type(SpinOrbIdx_t), intent(in), optional :: add_particles
        integer, intent(in), optional :: n_total
        type(Spin_t), intent(in), optional :: m_s

        type(SpinOrbIdx_t) :: possible_holes

        type(SpinOrbIdx_t), allocatable :: splitted_det_I(:)
        integer, allocatable :: spaces(:)

        allocate(splitted_det_I(get_nGAS(GAS_spec)))
        splitted_det_I = split_per_GAS(GAS_spec, det_I)

        ! Note, that non-present optional arguments can be passed
        ! into optional arguments without checking!
        spaces = get_possible_spaces(&
             GAS_spec, splitted_det_I, add_holes=add_holes, &
             add_particles=add_particles, n_total=n_total)

        if (size(spaces) == 0) then
            possible_holes%idx = [integer::]
            return
        end if

        block
            integer :: i, lower_bound, upper_bound
            type(SpinOrbIdx_t) :: possible_values, occupied

            if (spaces(1) == 1) then
                lower_bound = 1
            else
                lower_bound = GAS_spec%n_orbs(spaces(1) - 1) + 1
            end if
            upper_bound = GAS_spec%n_orbs(spaces(2))

            occupied = SpinOrbIdx_t([(splitted_det_I(i)%idx, i = spaces(1), spaces(2))], m_s)
            possible_values = SpinOrbIdx_t(SpatOrbIdx_t([(i, i = lower_bound, upper_bound)]), m_s)

            if (present(add_particles)) then
                possible_holes%idx = complement(possible_values%idx, union(occupied%idx, add_particles%idx))
            else
                possible_holes%idx = complement(possible_values%idx, occupied%idx)
            end if
        end block
    end function



    subroutine generate_nGAS_excitation(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex_mat, tParity, pGen, hel, store, part_type)
        ! particle number conserving GAS excitation generator:
        ! we create only excitations, that preserver the number of electrons within
        ! each active space

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'generate_nGAS_excitation'

        type(SpinOrbIdx_t) :: det_J
        class(Excitation_t), allocatable :: exc

        logical :: par

#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        ! single or double excitation?
        @:ASSERT(0.0_dp <= pDoubles .and. pDoubles <= 1.0_dp, pDoubles)
        if (genrand_real2_dSFMT() >= pDoubles) then
            ic = 1
            associate(r => [genrand_real2_dSFMT(), genrand_real2_dSFMT()])
                call gen_exc_single(GAS_specification, SpinOrbIdx_t(nI), ilutI, r, &
                                    nJ, ilutJ, ex_mat, par, pgen)
            end associate
            pgen = pgen * (1.0_dp - pDoubles)
        else
            ic = 2
            write(*, *) 'call gen_exc_double'
            associate(r => [genrand_real2_dSFMT(), genrand_real2_dSFMT()])
                call gen_exc_double(GAS_specification, SpinOrbIdx_t(nI), ilutI, r, &
                                    nJ, ilutJ, ex_mat, par, pgen)
            end associate
            pgen = pgen * pDoubles
        end if
        @:ASSERT(0.0_dp <= pgen .and. pgen <= 1.0_dp, pgen)

    end subroutine generate_nGAS_excitation

    subroutine gen_exc_single(GAS_spec, det_I, ilutI, r, nJ, ilutJ, ex_mat, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp), intent(in) :: r(2)
        integer, intent(out) :: nJ(nel), ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_exc_SingleExc_t'

        type(SingleExc_t) :: exc
        type(SpinOrbIdx_t) :: possible_holes
        real(dp) :: pgen_particle, pgen_hole
        real(dp), allocatable :: c_sum(:)
        integer :: i, elec

        ! two random numbers \in [0, 1]
        @:ASSERT(all(0.0_dp <= r .and. r <= 1.0))

        ! Pick any random electron
        elec = int(r(1) * nEl) + 1
        exc%val(1) = det_I%idx(elec)
        pgen_particle = 1.0_dp / real(nEl, kind=dp)

        ! Get a hole with the same spin projection
        ! while fullfilling GAS-constraints.
        associate(deleted => SpinOrbIdx_t([exc%val(1)]))
            possible_holes = get_possible_holes(&
                GAS_spec, det_I, add_holes=deleted, m_s=calc_spin(deleted))
        end associate


        if (size(possible_holes) == 0) then
            exc%val(2) = 0
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt \in possible_holes
        c_sum = get_cumulative_list(det_I, exc, possible_holes)
        call draw_from_cum_list(c_sum, r(2), i, pgen_hole)
        exc%val(2) = possible_holes%idx(i)


        if (i /= 0) then
            call make_single(det_I%idx, nJ, elec, exc%val(2), ex_mat, par)
            ilutJ = excite(ilutI, exc)
        else
            ilutJ = 0
        end if

        pgen = pgen_particle * pgen_hole
        write(*, *) 'succesful single', exc%val
    end subroutine

    subroutine gen_exc_double(GAS_spec, det_I, ilutI, r, nJ, ilutJ, ex_mat, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp), intent(in) :: r(2)
        integer, intent(out) :: nJ(nel), ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_exc_double'

        type(DoubleExc_t) :: exc, reverted_exc
        type(SpinOrbIdx_t) :: possible_holes
        ! Spin of second electron
        type(Spin_t) :: m_s_1, m_s_2
        real(dp) :: pgen_particles, &
            ! These are arrays, because the pgen might be different if we
            ! pick AB or BA.
            ! pgen_first_pick == p(A) == p(B)
            ! pgen_second_pick == [p(B | A), p(A | B)]
            pgen_first_pick, pgen_second_pick(2)
        real(dp), allocatable :: c_sum(:)
        integer :: i, elec

        integer :: elecs(2), sym_product, ispn, sum_ml, tgt(2)
        integer :: ms, nJBase(nel)
        logical :: tExchange

        ! two random numbers \in [0, 1]
        @:ASSERT(all(0.0_dp <= r .and. r <= 1.0))


        call pick_biased_elecs(det_I%idx, elecs, exc%val(1, :), &
                               sym_product, ispn, sum_ml, pgen_particles)
        @:ASSERT(exc%val(1, 1) /= exc%val(2, 1), exc%val)

        associate(deleted => SpinOrbIdx_t(exc%val(1, :)), &
                  deleted_spin => Spin_t(sum(calc_spin_raw(exc%val(1, :)))))
        ! Get possible holes for the first particle, while fullfilling GAS-constraints.
        ! and knowing that a second particle will be created afterwards!
        possible_holes = get_possible_holes(&
            GAS_spec, det_I, add_holes=deleted, n_total=2)
        @:ASSERT(disjoint(possible_holes%idx, det_I%idx))

        write(*, *) 'PEW', 1
        if (size(possible_holes) == 0) then
            write(*, *) 'ASDF', 1
            pgen = pgen_particles
            call zeroResult()
            return
        end if

        ! Pick randomly one hole with arbitrary spin
        exc%val(2, 1) = possible_holes%idx(int(r(1) * real(size(possible_holes), kind=dp)) + 1)
        pgen_first_pick = 1.0_dp / real(size(possible_holes), dp)
        m_s_1 = Spin_t(calc_spin_raw(exc%val(2, 1)))
        @:ASSERT(any(m_s_1%m_s == [spin_values%alpha, spin_values%beta]))


        ! Pick second hole.
        ! The total spin projection of the created particles has to add up
        ! to the total spin projection of the deleted particles.
        m_s_2 = Spin_t(deleted_spin%m_s - m_s_1%m_s)
        @:ASSERT(any(m_s_2%m_s == [spin_values%alpha, spin_values%beta]), m_s_2%m_s)

        ! Get possible holes for the second particle,
        ! while fullfilling GAS- and Spin-constraints.
        possible_holes = get_possible_holes(&
                GAS_spec, det_I, add_holes=deleted, &
                add_particles=SpinOrbIdx_t([exc%val(2, 1)]), &
                n_total=1, m_s=m_s_2)
        @:ASSERT(disjoint(possible_holes%idx, [exc%val(2, 1)]))
        @:ASSERT(disjoint(possible_holes%idx, det_I%idx), GAS_spec%n_orbs, GAS_spec%n_min, GAS_spec%n_max, det_I%idx, deleted%idx, exc%val(2, 1), m_s_2%m_s)

        write(*, *) 'PEW', 2
        if (size(possible_holes) == 0) then
            write(*, *) 'ASDF', 2
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt2 in possible_holes
        c_sum = get_cumulative_list(det_I, exc, possible_holes)
        call draw_from_cum_list(c_sum, r(2), i, pgen_second_pick(1))

        if (i /= 0) then
            exc%val(2, 2) = possible_holes%idx(i)
            write(*, *) 'found tgt2', exc%val
        else
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if

        ! We could have picked the holes the other way round and have to
        ! determine the probability of picking tgt1 with spin m_s_1 upon picking tgt2 first.
        associate (src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                   src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
            @:ASSERT(src1 /= tgt2 .and. src2 /= tgt2, src1, tgt2, src2)
            possible_holes = get_possible_holes(&
                    GAS_spec, det_I, add_holes=deleted, &
                    add_particles=SpinOrbIdx_t([tgt2]), &
                    n_total=1, m_s=m_s_1)
            @:ASSERT(disjoint(possible_holes%idx, [tgt2]))
            @:ASSERT(disjoint(possible_holes%idx, det_I%idx))

            if (size(possible_holes) == 0) then
                pgen = pgen_particles * pgen_first_pick * pgen_second_pick(1)
                call zeroResult()
                return
            end if
            ! Possible_holes has to contain tgt1.
            ! we look up its index with binary search
            i = binary_search_first_ge(possible_holes%idx, tgt1)
            @:ASSERT(i /= -1, tgt1, possible_holes%idx)

            reverted_exc = DoubleExc_t(src1=src1, tgt1=tgt2, src2=src2)
            c_sum = get_cumulative_list(det_I, reverted_exc, possible_holes)
            if (i == 1) then
                pgen_second_pick(2) = c_sum(1)
            else
                pgen_second_pick(2) = (c_sum(i) - c_sum(i - 1))
            end if

            if (i /= 0) then
                call make_double(det_I%idx, nJ, elecs(1), elecs(2), tgt1, tgt2, ex_mat, par)
                ilutJ = excite(ilutI, exc)
            else
                ilutJ = 0
            end if
        end associate
        end associate

        pgen = pgen_particles * pgen_first_pick * sum(pgen_second_pick)
        write(*, *) 'succesfull double', exc%val
        contains

            subroutine zeroResult()
                integer :: src_copy(2)

                src_copy(:) = exc%val(1, :)
                call sort(src_copy)
                ex_mat(1, :2) = src_copy
                ex_mat(2, :2) = exc%val(2, :)
                nJ(1) = 0
                ilutJ = 0_n_int
            end subroutine zeroResult
    end subroutine gen_exc_double

    DEBUG_IMPURE subroutine draw_from_cum_list(c_sum, r, idx, pgen)
        real(dp), intent(in) :: c_sum(:), r
        integer, intent(out) :: idx
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'draw_from_cum_list'

        @:ASSERT((c_sum(size(c_sum)) .isclose. 0.0_dp) &
          .or. (c_sum(size(c_sum)) .isclose. 1.0_dp))
        @:ASSERT(is_sorted(c_sum))

        ! there might not be such an excitation
        if (c_sum(size(c_sum)) > 0) then
            ! find the index of the target hole in the cumulative list
            idx = binary_search_first_ge(c_sum, r)

            ! adjust pgen with the probability for picking tgt from the cumulative list
            if (idx == 1) then
                pgen = c_sum(1)
            else
                pgen = c_sum(idx) - c_sum(idx - 1)
            end if
        else
            idx = 0
            pgen = 0.0_dp
        end if
    end subroutine


#:for excitation_t in ExcitationTypes
    function get_cumulative_list_${excitation_t}$(det_I, incomplete_exc, possible_holes) result(cSum)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(${excitation_t}$), intent(in) :: incomplete_exc
        type(SpinOrbIdx_t), intent(in) :: possible_holes
        real(dp) :: cSum(size(possible_holes))
        character(*), parameter :: this_routine = 'get_cumulative_list_${excitation_t}$'

        real(dp) :: previous
        type(${excitation_t}$) :: exc
        integer :: i

        exc = incomplete_exc
        @:ASSERT(last_tgt_unknown(exc))

        ! build the cumulative list of matrix elements <src|H|tgt>
        previous = 0.0_dp
        do i = 1, size(possible_holes)
            call set_last_tgt(exc, possible_holes%idx(i))

            @:ASSERT(defined(exc), exc%val)
            cSum(i) = abs(sltcnd_excit(det_I, exc, .false.)) + previous
            previous = cSum(i)
        end do

        ! Normalize
        if (near_zero(cSum(size(cSum)))) then
            cSum(:) = 0.0_dp
        else
            cSum(:) = cSum(:) / cSum(size(cSum))
        end if
    end function get_cumulative_list_${excitation_t}$
#:endfor
end module gasci
