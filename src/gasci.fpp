#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']
#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']

module gasci
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
                        near_zero, cumsum, operator(.isclose.)
    use sort_mod, only: sort
    use DetBitOps, only: ilut_lt, ilut_gt
    use sets_mod, only: is_sorted, complement, union, disjoint
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pSingles, pDoubles, excit_gen_store_type
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
                                get_last_tgt, set_last_tgt, excite, defined, UNKNOWN
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, SpinProj_t, size, calc_spin, &
                           calc_spin_raw, alpha, beta, sum, operator(==), operator(/=), &
                           operator(+), operator(-), lex_leq, to_ilut, write_det
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    implicit none

    private
    public :: &
        possible_GAS_exc_gen, GAS_exc_gen, &
        is_valid, is_connected, GAS_specification, GASSpec_t, &
        get_nGAS, user_input_GAS_exc_gen, &
        generate_nGAS_excitation, &
        contains_det, count_per_GAS, &
        get_possible_spaces, get_possible_holes, split_per_GAS, &
        get_available_singles, get_available_doubles, operator(==), &
        operator(/=), gen_all_excits, &
        get_iGAS, operator(.contains.)

    !> Speficies the GAS spaces.
    !> The indices are:
    !>  n_orbs_per_GAS(1:nGAS), n_min(1:nGAS), n_max(1:nGAS), GAS_table(1 : nBasis)
    !> n_orbs_per_GAS(iGAS) specifies how many orbitals are in
    !> the `iGAS` GAS space in the `iRep` Irrep.
    !> n_min(iGAS) specifies the cumulated! minimum particle number per GAS space.
    !> n_max(iGAS) specifies the cumulated! maximum particle number per GAS space.
    !> GAS_table(i) returns the GAS space for the i-th spin orbital
    type :: GASSpec_t
        integer, allocatable :: n_orbs(:), n_min(:), n_max(:), GAS_table(:)

        !> splitted_orbitals orbitals is the preimage of GAS_specification%GAS_table.
        !> An array that contains the spin orbitals per GAS space.
        !> splitted_orbitals(1 : nGAS)

        ! NOTE: Should be allocatable, but has to be pointer for constraints in GFortran.
        ! (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87568)
        type(SpinOrbIdx_t), private, pointer :: splitted_orbitals(:) => null()
    contains
        procedure :: init => init_GAS_spec
        ! NOTE: Should be a final procedure, but not possible because of GFortran 5
        procedure :: destroy => destroy_GAS_spec
    end type

    type(GASSpec_t) :: GAS_specification

    type :: GAS_exc_gen_t
        integer :: val
    end type

    type :: possible_GAS_exc_gen_t
        type(GAS_exc_gen_t) :: &
            DISCONNECTED = GAS_exc_gen_t(1), &
            GENERAL = GAS_exc_gen_t(2)
    end type

    type(possible_GAS_exc_gen_t), parameter :: possible_GAS_exc_gen = possible_GAS_exc_gen_t()

    type(GAS_exc_gen_t) :: GAS_exc_gen = possible_GAS_exc_gen%GENERAL
    type(GAS_exc_gen_t), allocatable :: user_input_GAS_exc_gen

    interface operator(==)
        module procedure eq_GAS_exc_gen_t
    end interface

    interface operator(/=)
        module procedure neq_GAS_exc_gen_t
    end interface

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
    interface count_per_GAS
        #:for orb_idx_type in OrbIdxTypes
        module procedure particles_per_GAS_splitted_${orb_idx_type}$
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
        #:for Excitation_t in ExcitationTypes
        module procedure get_cumulative_list_${Excitation_t}$
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

    pure integer function get_nGAS(GAS_spec)
        type(GASSpec_t), intent(in) :: GAS_spec
        get_nGAS = size(GAS_spec%n_orbs)
    end function

    logical pure function is_valid(GAS_spec, n_particles, n_basis)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in), optional :: n_particles, n_basis

        logical :: shapes_match, nEl_correct, pauli_principle, monotonic, &
                   n_orbs_correct, GAS_sizes_match
        integer :: nGAS, iGAS, i

        associate (n_orbs => GAS_spec%n_orbs, n_min => GAS_spec%n_min, &
                   n_max => GAS_spec%n_max)

            nGAS = get_nGAS(GAS_spec)

            shapes_match = &
                all([size(n_orbs), size(n_min), size(n_max)] == nGAS) &
                .and. maxval(GAS_spec%GAS_table) == size(n_orbs)

            if (present(n_particles)) then
                nEl_correct = all([n_min(nGAS), n_max(nGAS)] == n_particles)
            else
                nEl_correct = n_min(nGAS) == n_max(nGAS)
            end if

            block
                integer :: GAS_sizes(nGAS), iGAS
                GAS_sizes = 0
                do i = 1, nBasis - 1, 2
                    iGAS = GAS_spec%GAS_table(i)
                    GAS_sizes(iGAS) = GAS_sizes(iGAS) + 1
                end do
                GAS_sizes_match = all(cumsum(GAS_sizes) == n_orbs)
            end block

            pauli_principle = all(n_min(:) <= n_orbs(:) * 2)

            if (nGAS >= 2) then
                monotonic = all([all(n_orbs(2:) >= n_orbs(:nGAS - 1)), &
                                 all(n_max(2:) >= n_max(:nGAS - 1)), &
                                 all(n_min(2:) >= n_min(:nGAS - 1))])
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
        res = any(GAS_spec%n_min(:) /= GAS_spec%n_max(:))
    end function

    subroutine init_GAS_spec(GAS_spec)
        class(GASSpec_t), intent(inout) :: GAS_spec
        character(*), parameter :: this_routine = 'init_GAS_spec'

        integer :: i
        type(SpinOrbIdx_t) :: all_orbitals

        associate(n_spin_orbs => 2 * GAS_spec%n_orbs(size(GAS_spec%n_orbs)))
            all_orbitals = SpinOrbIdx_t([(i, i=1, n_spin_orbs)])
        end associate

        allocate(GAS_spec%splitted_orbitals(get_nGAS(GAS_spec)))
        GAS_spec%splitted_orbitals = split_per_GAS(GAS_spec, all_orbitals)

        @:ASSERT(is_valid(GAS_spec))
    end subroutine

    subroutine destroy_GAS_spec(GAS_spec)
        class(GASSpec_t), intent(inout) :: GAS_spec
        character(*), parameter :: this_routine = 'destroy_GAS_spec'

        integer :: iGAS

        if (associated(GAS_spec%splitted_orbitals)) then
            do iGAS = 1, size(GAS_spec%splitted_orbitals)
                deallocate(GAS_spec%splitted_orbitals(iGAS)%idx)
            end do
            deallocate(GAS_spec%splitted_orbitals)
            nullify (GAS_spec%splitted_orbitals)
        end if
        ! The other components are allocatables and not pointers
        ! and are automatically cleaned up.
    end subroutine

    pure function get_iGAS_SpatOrbIdx_t(GAS_spec, idx) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpatOrbIdx_t), intent(in) :: idx
        integer :: res(size(idx))
        res = GAS_spec%GAS_table(idx%idx * 2 - 1)
    end function

    pure function get_iGAS_SpinOrbIdx_t(GAS_spec, idx) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: idx
        integer :: res(size(idx))
        res = GAS_spec%GAS_table(idx%idx)
    end function

#:for orb_idx_type in OrbIdxTypes
    pure function particles_per_GAS_splitted_${orb_idx_type}$ (splitted_det_I) result(n_particle)
        type(${orb_idx_type}$), intent(in) :: splitted_det_I(:)
        integer :: n_particle(size(splitted_det_I))
        integer :: i
        n_particle = [(size(splitted_det_I(i)), i=1, size(splitted_det_I))]
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    function split_per_GAS_${orb_idx_type}$ (GAS_spec, occupied) result(splitted_orbitals)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        type(${orb_idx_type}$), allocatable :: splitted_orbitals(:)

        integer :: iel, iGAS, GAS_table(size(occupied)), n_particle(get_nGAS(GAS_spec)), counter(get_nGAS(GAS_spec))

        allocate(splitted_orbitals(get_nGAS(GAS_spec)))
        GAS_table = get_iGAS(GAS_spec, occupied)

        n_particle = 0
        do iel = 1, size(GAS_table)
            iGAS = GAS_table(iel)
            n_particle(iGAS) = n_particle(iGAS) + 1
        end do

        do iGAS = 1, size(n_particle)
            allocate(splitted_orbitals(iGAS)%idx(n_particle(iGAS)))
        end do

        counter = 1
        do iel = 1, size(occupied)
            iGAS = GAS_table(iel)
            splitted_orbitals(iGAS)%idx(counter(iGAS)) = occupied%idx(iel)
            counter(iGAS) = counter(iGAS) + 1
        end do
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    function contains_det_${orb_idx_type}$ (GAS_spec, occupied) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        logical :: res

        !> Cumulated number of particles per iGAS
        integer :: cum_n_particle(get_nGAS(GAS_spec)), i

        associate(splitted => split_per_GAS(GAS_spec, occupied))
            cum_n_particle = cumsum([(size(splitted(i)), i=1, get_nGAS(GAS_spec))])
        end associate

        res = all(GAS_spec%n_min(:) <= cum_n_particle(:) &
                  .and. cum_n_particle(:) <= GAS_spec%n_max(:))
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    function get_possible_spaces_${orb_idx_type}$ (GAS_spec, splitted_det_I, add_holes, add_particles, n_total) result(spaces)
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

        associate(A => count_per_GAS(splitted_det_I))
            if (present(add_holes) .and. present(add_particles)) then
                associate(B => count_per_GAS(split_per_GAS(GAS_spec, add_holes)), &
                           C => count_per_GAS(split_per_GAS(GAS_spec, add_particles)))
                    cum_n_particle = cumsum(A - B + C)
                end associate
            else if (present(add_holes)) then
                associate(B => count_per_GAS(split_per_GAS(GAS_spec, add_holes)))
                    cum_n_particle = cumsum(A - B)
                end associate
            else if (present(add_particles)) then
                associate(C => count_per_GAS(split_per_GAS(GAS_spec, add_particles)))
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

        ! We assume that it is possible to create a particle at least in
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

    function get_possible_holes(GAS_spec, det_I, add_holes, add_particles, n_total, excess) result(possible_holes)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SpinOrbIdx_t), intent(in), optional :: add_holes
        type(SpinOrbIdx_t), intent(in), optional :: add_particles
        integer, intent(in), optional :: n_total
        type(SpinProj_t), intent(in), optional :: excess
        character(*), parameter :: this_routine = 'get_possible_holes'

        type(SpinOrbIdx_t) :: possible_holes

        type(SpinOrbIdx_t), allocatable :: splitted_det_I(:)
        integer, allocatable :: spaces(:)
        integer :: n_total_

        @:def_default(n_total_, n_total, 1)
        @:ASSERT(1 <= n_total_)
        if (present(excess)) then
            @:ASSERT(abs(excess%val) <= n_total_, excess, n_total_)
        end if

        allocate(splitted_det_I(get_nGAS(GAS_spec)))
        splitted_det_I = split_per_GAS(GAS_spec, det_I)

        ! Note, that non-present optional arguments can be passed
        ! into optional arguments without checking!
        spaces = get_possible_spaces( &
                 GAS_spec, splitted_det_I, add_holes=add_holes, &
                 add_particles=add_particles, n_total=n_total_)

        if (size(spaces) == 0) then
            possible_holes%idx = [integer::]
            return
        end if

        associate(splitted_orbitals => GAS_spec%splitted_orbitals)
            block
                integer :: i, k, iGAS, incr, curr_value, iGAS_min_val
                integer :: L(spaces(1):spaces(2)), counter(spaces(1):spaces(2))
                type(SpinOrbIdx_t) :: possible_values
                type(SpinProj_t) :: m_s

                m_s = SpinProj_t(0)
                if (present(excess)) then
                    if (abs(excess%val) == n_total_) m_s = -SpinProj_t(sign(1, excess%val))
                end if

                L = count_per_GAS(splitted_orbitals(spaces(1):spaces(2)))
                if (m_s == beta) then
                    allocate(possible_values%idx(sum(L) .div.2))
                    counter = 1
                    incr = 2
                else if (m_s == alpha) then
                    allocate(possible_values%idx(sum(L) .div.2))
                    counter = 2
                    incr = 2
                else
                    allocate(possible_values%idx(sum(L)))
                    counter = 1
                    incr = 1
                end if

                ! Here we merge the values from splitted_orbitals sortedly into
                ! possible_values
                i = 1
                do while (any(counter <= L))
                    curr_value = huge(curr_value)
                    do iGAS = spaces(1), spaces(2)
                        if (counter(iGAS) <= size(splitted_orbitals(iGAS))) then
                            if (splitted_orbitals(iGAS)%idx(counter(iGAS)) < curr_value) then
                                curr_value = splitted_orbitals(iGAS)%idx(counter(iGAS))
                                iGAS_min_val = iGAS
                            end if
                        end if
                    end do

                    counter(iGAS_min_val) = counter(iGAS_min_val) + incr
                    possible_values%idx(i) = curr_value
                    i = i + 1
                end do

                if (present(add_particles)) then
                    possible_holes%idx = complement(possible_values%idx, union(det_I%idx, add_particles%idx))
                else
                    possible_holes%idx = complement(possible_values%idx, det_I%idx)
                end if
            end block
        end associate
    end function

    subroutine generate_nGAS_excitation(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex_mat, tParity, pGen, hel, store, part_type)
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

        @:unused_var(exFlag, part_type, store)

#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        if (genrand_real2_dSFMT() >= pDoubles) then
            ic = 1
            call gen_exc_single(GAS_specification, SpinOrbIdx_t(nI), ilutI, &
                                nJ, ilutJ, ex_mat, tParity, pgen)
            pgen = pgen * (1.0_dp - pDoubles)
        else
            ic = 2
            call gen_exc_double(GAS_specification, SpinOrbIdx_t(nI), ilutI, &
                                nJ, ilutJ, ex_mat, tParity, pgen)
            pgen = pgen * pDoubles
        end if

        @:ASSERT(0.0_dp <= pgen .and. pgen <= 1.0_dp, pgen)
        @:ASSERT(all(ex_mat(:, :ic) /= UNKNOWN) .or. nJ(1) == 0)
    end subroutine generate_nGAS_excitation

    subroutine gen_exc_single(GAS_spec, det_I, ilutI, nJ, ilutJ, ex_mat, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
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

        ! Pick any random electron
        elec = int(genrand_real2_dSFMT() * nEl) + 1
        exc%val(1) = det_I%idx(elec)
        pgen_particle = 1.0_dp / real(nEl, kind=dp)

        ! Get a hole with the same spin projection
        ! while fullfilling GAS-constraints.
        block
            type(SpinOrbIdx_t) :: deleted
            deleted = SpinOrbIdx_t([exc%val(1)])
            possible_holes = get_possible_holes( &
                             GAS_spec, det_I, add_holes=deleted, excess=-calc_spin_raw(exc%val(1)))
        end block

        if (size(possible_holes) == 0) then
            exc%val(2) = 0
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt \in possible_holes
        c_sum = get_cumulative_list(det_I, exc, possible_holes)
        call draw_from_cum_list(c_sum, i, pgen_hole)
        @:ASSERT(i == 0 .neqv. (0.0_dp < pgen_hole .and. pgen_hole <= 1.0_dp))
        if (i /= 0) then
            exc%val(2) = possible_holes%idx(i)
            call make_single(det_I%idx, nJ, elec, exc%val(2), ex_mat, par)
            ilutJ = excite(ilutI, exc)
        else
            nJ = 0
            ilutJ = 0
        end if

        pgen = pgen_particle * pgen_hole
        @:ASSERT(all(nJ == 0) .neqv. 0.0_dp < pgen .and. pgen <= 1.0_dp)
    end subroutine

    subroutine gen_exc_double(GAS_spec, det_I, ilutI, nJ, ilutJ, ex_mat, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_exc_double'

        type(DoubleExc_t) :: exc, reverted_exc
        type(SpinOrbIdx_t) :: possible_holes, deleted
        ! Spin of second electron
        type(SpinProj_t) :: m_s_1, m_s_2
        real(dp) :: pgen_particles, &
                    ! These are arrays, because the pgen might be different if we
                    ! pick AB or BA. Let i, j denote the particles and a, b the holes.
                    ! pgen_particles == p (i j)
                    ! pgen_first_pick == p(a | i j) == p(b | i j)
                    ! pgen_second_pick == [p(b | a i j), p(a | b i j)]
                    ! pgen == p_double * p (i j) * p(a | i j) * sum([p(b | a i j), p(a | b i j)])
                    !      == p_double * p (i j) * p(b | i j) * sum([p(b | a i j), p(a | b i j)])
                    pgen_first_pick, pgen_second_pick(2)
        real(dp) :: r
        real(dp), allocatable :: c_sum(:)
        integer :: i, elec

        integer :: elecs(2), sym_product, ispn, sum_ml, tgt(2)
        integer :: ms, nJBase(nel)
        logical :: tExchange

        call pick_biased_elecs(det_I%idx, elecs, exc%val(1, :), &
                               sym_product, ispn, sum_ml, pgen_particles)
        @:ASSERT(exc%val(1, 1) /= exc%val(2, 1), exc%val)

        deleted = SpinOrbIdx_t(exc%val(1, :))
        ! Get possible holes for the first particle, while fullfilling and spin- and GAS-constraints.
        ! and knowing that a second particle will be created afterwards!
        possible_holes = get_possible_holes( &
                         GAS_spec, det_I, add_holes=deleted, n_total=2, excess=-sum(calc_spin(deleted)))
        @:ASSERT(disjoint(possible_holes%idx, det_I%idx))

        if (size(possible_holes) == 0) then
            pgen = pgen_particles
            call zeroResult()
            return
        end if

        r = genrand_real2_dSFMT()
        ! Pick randomly one hole with arbitrary spin
        exc%val(2, 1) = possible_holes%idx(int(r * real(size(possible_holes), kind=dp)) + 1)
        pgen_first_pick = 1.0_dp / real(size(possible_holes), dp)
        m_s_1 = calc_spin_raw(exc%val(2, 1))
        @:ASSERT(any(m_s_1 == [alpha, beta]))

        ! Pick second hole.
        ! The total spin projection of the created particles has to add up
        ! to the total spin projection of the deleted particles.
        ! Get possible holes for the second particle,
        ! while fullfilling GAS- and Spin-constraints.
        block
            type(SpinProj_t) :: excess

            excess = m_s_1 - sum(calc_spin(deleted))

            @:ASSERT(abs(excess%val) <= 1)

            possible_holes = get_possible_holes( &
                             GAS_spec, det_I, add_holes=deleted, &
                             add_particles=SpinOrbIdx_t([exc%val(2, 1)]), &
                             n_total=1, excess=excess)
        end block

        @:ASSERT(disjoint(possible_holes%idx, [exc%val(2, 1)]))
        @:ASSERT(disjoint(possible_holes%idx, det_I%idx))

        if (size(possible_holes) == 0) then
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt2 in possible_holes
        c_sum = get_cumulative_list(det_I, exc, possible_holes)
        call draw_from_cum_list(c_sum, i, pgen_second_pick(1))

        if (i /= 0) then
            exc%val(2, 2) = possible_holes%idx(i)
        else
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if
        @:ASSERT(defined(exc))

        ! We could have picked the holes the other way round and have to
        ! determine the probability of picking tgt1 with spin m_s_1 upon picking tgt2 first.
        associate(src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                   src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
            @:ASSERT(src1 /= tgt2 .and. src2 /= tgt2, src1, tgt2, src2)
            possible_holes = get_possible_holes( &
                             GAS_spec, det_I, add_holes=deleted, &
                             add_particles=SpinOrbIdx_t([tgt2]), &
                             n_total=1, excess=calc_spin_raw(tgt2) - sum(calc_spin(deleted)))
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

        pgen = pgen_particles * pgen_first_pick * sum(pgen_second_pick)

        @:ASSERT(0.0_dp < pgen .and. pgen <= 1.0_dp)
        @:ASSERT(all(ex_mat(:, :2) /= UNKNOWN))

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

    DEBUG_IMPURE subroutine draw_from_cum_list(c_sum, idx, pgen)
        real(dp), intent(in) :: c_sum(:)
        integer, intent(out) :: idx
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'draw_from_cum_list'

        @:ASSERT((c_sum(size(c_sum)) .isclose. 0.0_dp) &
            .or. (c_sum(size(c_sum)) .isclose. 1.0_dp))
        @:ASSERT(is_sorted(c_sum))

        ! there might not be such an excitation
        if (c_sum(size(c_sum)) > 0) then
            ! find the index of the target hole in the cumulative list
            idx = binary_search_first_ge(c_sum, genrand_real2_dSFMT())
            @:ASSERT(1 <= idx .and. idx <= size(c_sum))

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

    DEBUG_IMPURE function get_available_singles(GAS_spec, det_I) result(singles_exc_list)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SpinOrbIdx_t), allocatable :: singles_exc_list(:)
        character(*), parameter :: this_routine = 'get_available_singles'

        type(SpinOrbIdx_t) :: possible_holes
        type(SpinOrbIdx_t), allocatable :: tmp_buffer(:)
        integer :: i, j, i_buffer, src1, tgt1
        type(SpinProj_t) :: m_s_1

        @:ASSERT(GAS_spec .contains. det_I)

        i_buffer = 0
        do i = 1, size(det_I%idx)
            src1 = det_I%idx(i)
            m_s_1 = calc_spin_raw(src1)
            possible_holes = get_possible_holes( &
                             GAS_spec, det_I, add_holes=SpinOrbIdx_t([src1]), &
                             excess=-m_s_1)
            do j = 1, size(possible_holes)
                tgt1 = possible_holes%idx(j)
                i_buffer = i_buffer + 1
                call grow_assign(tmp_buffer, i_buffer, &
                                 excite(det_I, SingleExc_t(src1, tgt1)))
            end do
        end do
        singles_exc_list = tmp_buffer(:i_buffer)

        @:sort(SpinOrbIdx_t, singles_exc_list, lex_leq)
    end function

    DEBUG_IMPURE function get_available_doubles(GAS_spec, det_I) result(doubles_exc_list)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SpinOrbIdx_t), allocatable :: doubles_exc_list(:)
        character(*), parameter :: this_routine = 'get_available_doubles'

        type(SpinOrbIdx_t) :: possible_holes(2), deleted
        type(SpinOrbIdx_t), allocatable :: tmp_buffer(:)
        integer :: i, j, k, l, i_buffer, src1, src2, tgt1, tgt2
        type(SpinProj_t) :: m_s_1

        @:ASSERT(GAS_spec .contains. det_I)

        i_buffer = 0
        do i = 1, size(det_I)
            do j = i + 1, size(det_I)
                src1 = det_I%idx(i)
                src2 = det_I%idx(j)
                deleted = SpinOrbIdx_t([src1, src2])
                possible_holes(1) = get_possible_holes(GAS_spec, det_I, &
                                                       add_holes=deleted, excess=-sum(calc_spin(deleted)), &
                                                       n_total=2)
                @:ASSERT(disjoint(possible_holes(1)%idx, det_I%idx))
                do k = 1, size(possible_holes(1))
                    tgt1 = possible_holes(1)%idx(k)
                    m_s_1 = calc_spin_raw(tgt1)
                    @:ASSERT(any(m_s_1 == [alpha, beta]))

                    possible_holes(2) = get_possible_holes( &
                                        GAS_spec, det_I, add_holes=deleted, &
                                        add_particles=SpinOrbIdx_t([tgt1]), &
                                        n_total=1, excess=m_s_1 - sum(calc_spin(deleted)))

                    @:ASSERT(disjoint(possible_holes(2)%idx, [tgt1]))
                    @:ASSERT(disjoint(possible_holes(2)%idx, det_I%idx))

                    do l = 1, size(possible_holes(2))
                        tgt2 = possible_holes(2)%idx(l)
                        i_buffer = i_buffer + 1
                        call grow_assign( &
                            tmp_buffer, i_buffer, &
                            excite(det_I, DoubleExc_t(src1, tgt1, src2, tgt2)))
                    end do
                end do
            end do
        end do

        doubles_exc_list = tmp_buffer(:i_buffer)

        @:sort(SpinOrbIdx_t, doubles_exc_list, lex_leq)

        ! Remove double appearances
        j = 1
        tmp_buffer(j) = doubles_exc_list(1)
        do i = 2, size(doubles_exc_list)
            if (any(doubles_exc_list(i - 1) /= doubles_exc_list(i))) then
                j = j + 1
                tmp_buffer(j) = doubles_exc_list(i)
            end if
        end do
        doubles_exc_list = tmp_buffer(:j)
    end function

    subroutine grow_assign(lhs, i, rhs)
        type(SpinOrbIdx_t), intent(inout), allocatable :: lhs(:)
        integer, intent(in) :: i
        type(SpinOrbIdx_t), intent(in) :: rhs

        type(SpinOrbIdx_t), allocatable :: buffer(:)
        integer :: n
        real(dp), parameter :: grow_factor = 2.0_dp
        integer, parameter :: start_n = 10

        if (.not. allocated(lhs)) allocate(lhs(start_n))

        if (i > size(lhs)) then
            buffer = lhs(:)
            deallocate(lhs)
            allocate(lhs(int(i * grow_factor)))
            lhs(:size(buffer)) = buffer
        end if
        lhs(i) = rhs
    end subroutine

#:for excitation_t in ExcitationTypes
    function get_cumulative_list_${excitation_t}$ (det_I, incomplete_exc, possible_holes) result(cSum)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(${excitation_t}$), intent(in) :: incomplete_exc
        type(SpinOrbIdx_t), intent(in) :: possible_holes
        real(dp) :: cSum(size(possible_holes))
        character(*), parameter :: this_routine = 'get_cumulative_list_${excitation_t}$'

        real(dp) :: previous
        type(${excitation_t}$) :: exc
        integer :: i

        @:ASSERT(get_last_tgt(exc) == UNKNOWN)
        exc = incomplete_exc

        ! build the cumulative list of matrix elements <src|H|tgt>
        previous = 0.0_dp
        do i = 1, size(possible_holes)
            call set_last_tgt(exc, possible_holes%idx(i))
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

    logical pure function eq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        eq_GAS_exc_gen_t = lhs%val == rhs%val
    end function

    logical pure function neq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        neq_GAS_exc_gen_t = lhs%val /= rhs%val
    end function

    subroutine gen_all_excits(nI, n_excits, det_list)
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:, :)

        type(SpinOrbIdx_t) :: det_I
        type(SpinOrbIdx_t), allocatable :: singles(:), doubles(:)
        integer :: i, j, k

        det_I = SpinOrbIdx_t(nI)

        singles = get_available_singles(GAS_specification, det_I)
        doubles = get_available_doubles(GAS_specification, det_I)

        n_excits = size(singles) + size(doubles)
        allocate(det_list(0:niftot, n_excits))
        j = 1
        do i = 1, size(singles)
            det_list(:, j) = to_ilut(singles(i))
            j = j + 1
        end do

        do i = 1, size(doubles)
            det_list(:, j) = to_ilut(doubles(i))
            j = j + 1
        end do

        call sort(det_list, ilut_lt, ilut_gt)
    end subroutine gen_all_excits
end module gasci
