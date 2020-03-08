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
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, Spin_t, size
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    implicit none(type, external)

    private
    public :: is_valid, is_connected, GAS_specification, GASSpec_t, &
        init_GAS, clear_GAS, get_nGAS, &
!         generate_nGAS_excitation, &
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


    interface gen_exc
        #:for excitation_t in ExcitationTypes
            module procedure gen_exc_${excitation_t}$
        #:endfor
    end interface

    interface get_cumulative_list
        module procedure get_cumulative_list_SingleExc_t
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
    pure function split_per_GAS_${orb_idx_type}$(GAS_spec, occupied) result(splitted_orbitals)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: occupied

        type(${orb_idx_type}$) :: splitted_orbitals(get_nGAS(GAS_spec))

        integer :: iGAS, prev, i, GAS_table(size(occupied))

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
            splitted_orbitals(iGAS) = ${orb_idx_type}$(occupied%idx(prev + 1 : i))
        end do
    end function
#:endfor


#:for orb_idx_type in OrbIdxTypes
    pure function contains_det_${orb_idx_type}$(GAS_spec, occupied) result(res)
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
    DEBUG_IMPURE function get_possible_spaces_${orb_idx_type}$(GAS_spec, splitted_det_I, add_holes, n_particles) result(spaces)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(${orb_idx_type}$), intent(in) :: splitted_det_I(:)
        type(${orb_idx_type}$), intent(in), optional :: add_holes
        integer, intent(in), optional :: n_particles
        character(*), parameter :: this_routine = 'get_possible_spaces_${orb_idx_type}$'
        integer, allocatable :: spaces(:)

        !> Lower and upper bound for spaces where a particle can be created.
        !> If no particle can be created, then spaces == 0 .
        integer :: n_particles_, i, iGAS, lower_bound, upper_bound

        integer :: &
        !> Cumulated number of particles per iGAS
            cum_n_particle(get_nGAS(GAS_spec)), &
        !> Cumulated deficit per iGAS
            deficit(get_nGAS(GAS_spec)), &
        !> Cumulated vacant orbitals per iGAS
            vacant(get_nGAS(GAS_spec))

        @:def_default(n_particles_, n_particles, 1)

        associate(A => particles_per_GAS(splitted_det_I))
            if (present(add_holes)) then
                associate(splitted_holes => split_per_GAS(GAS_spec, add_holes))
                    associate(B => particles_per_GAS(splitted_holes))
                        cum_n_particle = cumsum(A - B)
                    end associate
                end associate
! TODO(Kai, Oskar):
! If one uncomments the following statement, the gasci unit test crashes.
! I don't understand why. :-(
                ! cum_n_particle = cumsum(particles_per_GAS(splitted_det_I) &
                !                         - particles_per_GAS(split_per_GAS(GAS_spec, add_holes)))
            else
                cum_n_particle = cumsum(A)
            end if
        end associate

        deficit = GAS_spec%n_min(:) - cum_n_particle(:)
        vacant = GAS_spec%n_max(:) - cum_n_particle(:)
        if (any(n_particles_ < deficit) .or. all(vacant < n_particles_)) then
            spaces = [integer::]
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

        if (lower_bound > upper_bound .or. lower_bound > get_nGAS(GAS_spec)) then
            spaces = [integer::]
        else
            spaces = [lower_bound, upper_bound]
        end if
    end function
#:endfor


    pure function get_possible_holes(GAS_spec, det_I, add_holes, n_particles, m_s) result(holes)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        ! Note, that non-present optional arguments can be passed
        ! into optional arguments without checking!
        type(SpinOrbIdx_t), intent(in), optional :: add_holes
        integer, intent(in), optional :: n_particles
        type(Spin_t), intent(in), optional :: m_s

        type(SpinOrbIdx_t) :: holes, splitted_det_I(get_nGAS(GAS_spec))

        integer, allocatable :: spaces(:)

        splitted_det_I = split_per_GAS(GAS_spec, det_I)
        spaces = get_possible_spaces(&
             GAS_spec, splitted_det_I, add_holes=add_holes, n_particles=n_particles)

        if (size(spaces) == 0) then
            holes = SpinOrbIdx_t([integer::])
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

            associate(splitted_occ => splitted_det_I(spaces(1) : spaces(2)))
                occupied = SpinOrbIdx_t([(splitted_occ(i)%idx, i = 1, size(splitted_occ))], m_s)
            end associate

            possible_values = SpinOrbIdx_t(SpatOrbIdx_t([(i, i = lower_bound, upper_bound)]), m_s)
            holes%idx = if_not_in(possible_values%idx, occupied%idx)
        end block

        contains
            !> Return all values of A that are not in B
            !> Assume:
            !>      1. All values of B appear in A.
            !>      2. A and B are sorted.
            pure function if_not_in(A, B) result(D)
                integer, intent(in) :: A(:), B(:)
                integer, allocatable :: D(:)

                integer :: i, j, l

                allocate(D(size(A) - size(B)))

                i = 1; j = 1; l = 1
                do while (l <= size(D))
                    if (j > size(B)) then
                        D(l) = A(i)
                        i = i + 1
                        l = l + 1
                    else if (A(i) /= B(j)) then
                        D(l) = A(i)
                        i = i + 1
                        l = l + 1
                    else if (A(i) == B(j)) then
                        i = i + 1
                        j = j + 1
                    end if
                end do
            end function
    end function



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
!         logical :: par
!
! #ifdef WARNING_WORKAROUND_
!         hel = 0.0_dp
! #endif
!         ! single or double excitation?
!         r = genrand_real2_dSFMT()
!         if (r >= pDoubles) then
!             exc = SingleExc_t()
!         else
!             exc = DoubleExc_t()
!         end if
!
!         select type(exc)
!         type is(SingleExc_t)
!             call gen_exc(GAS_specification, SpinOrbIdx_t(nI), exc, det_J, par, pgen)
!             pgen = pgen * (1.0_dp - pDoubles)
!             ic = 1
!         type is(DoubleExc_t)
!             call gen_exc(GAS_specification, SpinOrbIdx_t(nI), exc, det_J, par, pgen)
!             pgen = pgen * pDoubles
!             ic = 2
!         end select
!     end subroutine generate_nGAS_excitation

    subroutine gen_exc_SingleExc_t(GAS_spec, det_I, exc, det_J, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SingleExc_t), intent(out) :: exc
        type(SpinOrbIdx_t), intent(out) :: det_J
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        ! adjust pgen
        pgen = 1.0_dp / real(nEl, kind=dp)
        associate(src => exc%val(1))
            ! Pick any random electron
            src = det_I%idx(int(genrand_real2_dSFMT() * nEl) + 1)
        end associate
        ! Get a hole with the same spin projection
        ! while fullfilling GAS-constraints.
        ! call pick_weighted_hole_SingleExc_t(GAS_spec, det_I, exc, pgen)
    end subroutine

    subroutine gen_exc_DoubleExc_t(GAS_spec, det_I, exc, det_J, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(DoubleExc_t), intent(out) :: exc
        type(SpinOrbIdx_t), intent(out) :: det_J
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_exc_DoubleExc_t'
        call stop_all(this_routine, 'Not implemented')
    end subroutine

    subroutine pick_weighted_hole_SingleExc_t(GAS_spec, det_I, exc, pgen)
        ! pick a hole of nI with spin ms from the active space with index
        ! srcGASInd the random number is to be supplied as r
        ! nI is the source determinant, nJBase the one from which we obtain
        ! the ket of the matrix element by single excitation
        type(GASSpec_t) :: GAS_spec
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SingleExc_t), intent(inout) :: exc
        real(dp), intent(inout) :: pgen
        character(*), parameter :: this_routine = 'pick_weighted_hole_${excitation_t}$'

        real(dp) :: r
        type(SpinOrbIdx_t) :: possible_holes
        real(dp), allocatable :: cSum(:)

        ASSERT(last_tgt_unknown(exc))

        associate(src => exc%val(1), tgt => exc%val(2))

            possible_holes = get_possible_holes(&
                    GAS_spec, det_I, add_holes=SpinOrbIdx_t([src]))
            if (size(possible_holes) == 0) then
                tgt = 0
                return
            end if

            ! build the cumulative list of matrix elements <src|H|tgt>
            cSum = get_cumulative_list(det_I, exc, possible_holes)

            ! ! now, pick with the weight from the cumulative list
            ! r = genrand_real2_dSFMT() * cSum(nOrbs)
            !
            ! ! there might not be such an excitation
            ! if (cSum(nOrbs) > 0) then
            !     ! find the index of the target orbital in the gasList
            !     tgt = binary_search_first_ge(cSum, r)
            !
            !     ! adjust pgen with the probability for picking tgt from the cumulative list
            !     if (tgt == 1) then
            !         pgen = pgen * cSum(1)
            !     else
            !         pgen = pgen * (cSum(tgt) - cSum(tgt - 1))
            !     end if
            !
            !     ! convert to global orbital index
            !     tgt = GAS_list(tgt)
            ! else
            !     tgt = 0
            ! end if
        end associate
    end subroutine pick_weighted_hole_SingleExc_t


    function get_cumulative_list_SingleExc_t(det_I, incomplete_exc, possible_holes) result(cSum)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SingleExc_t), intent(in) :: incomplete_exc
        type(SpinOrbIdx_t), intent(in) :: possible_holes
        real(dp) :: cSum(size(possible_holes))
        character(*), parameter :: this_routine = 'get_cumulative_list_${excitation_t}$'

        real(dp) :: previous
        type(SingleExc_t) :: exc
        integer :: i

        exc = incomplete_exc
        ASSERT(last_tgt_unknown(exc))

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
    end function get_cumulative_list_SingleExc_t
end module gasci
