#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

! A list (order matters!) of k non-negative integers zero that sum to the positive integer
! n is an integer composition of n.
! If we apply this thinking to GAS wavefunctions, then k is the number of GAS spaces and
!   n is the number of particles.
! All possible compositions of 3 with 3 summands are for example
!      [3, 0, 0]
!      [2, 1, 0]
!      [2, 0, 1]
!      [1, 2, 0]
!      [1, 1, 1]
!      [1, 0, 2]
!      [0, 3, 0]
!      [0, 2, 1]
!      [0, 1, 2]
!      [0, 0, 3]
! The number of all possible compositions is defined as p(k, n) and given by p(k, n) = nCr(n + k - 1, k - 1)
!
! We assume in the following lexicographical decreasing order and assign a composition index
!   based on this order.
! The composition index of [1, 0, 2] is for example 6.
! Due to the lexicographical order it is possible to calculate the index for a given composition
!   by "jumping" over the leading terms.
! For example the index of [1, 0, 2] is given by "jumping"
!    over all terms with a leading 3 ([3, ?, ?]) and 2 ([2, ?, ?]) and then applying the same logic to the next index.
! There are p(k - 1, n - leading term) compositions to jump over.
!                  [3, ?, ?]     [2, ?, ?]   [1, 2, ?]    [1, 1, ?]
! idx([1, 0, 2]) = p(2, 0)    +    p(2, 1) +   p(1, 0) +    p(1, 1) +         1
!                =      1     +         2  +        1  +          1 +         1
!                = 6
!
! A supergroup is a composition that is allowed under constraints
!   to the cumulative minimum and maximum particle number.
! If the cumulative minimum is [0, 1, 3]
! and the cumulative maximum is [2, 2, 3]
! then the supergroups are:
! 3    [2, 0, 1]
! 5    [1, 1, 1]
! 6    [1, 0, 2]
! 8    [0, 2, 1]
! 9    [0, 1, 2]
! In the first column the composition index was written.
! We again assume lexicographical decreasing order and assign a supergroup index
!   based on this order.
! The supergroup index of [1, 0, 2] is for example 3, while the composition index was 6.
!
! There is a nice and elegant recursive solution to calculate the supergroup index,
! but in practice it is faster to calculate the composition index for a given supergroup
! and search the corresponding supergroup index in a precomputed list of composition indices
! for all allowed super groups.
module gasci_supergroup_index
    use constants, only: int64, n_int
    use util_mod, only: choose_i64, cumsum, binary_search_first_ge, custom_findloc, stop_all
    use bit_rep_data, only: nIfD
    use gasci, only: GASSpec_t, LocalGASSpec_t, CumulGASSpec_t, FlexibleGASSpec_t, CAS_spec
    use hash, only: hash_table_lookup
    use guga_bitRepOps, only: CSF_Info_t
    use growing_buffers, only: buffer_int_2D_t, buffer_int64_1D_t
    use composition_utils, only: n_compositions, get_compositions, composition_idx, &
        composition_from_idx, next_composition
    use excitation_types, only: Excite_1_t, Excite_2_t
    use guga_bitRepOps, only: CSF_Info_t
    use guga_data, only: DistinctDouble_t
    use generalized_supergroup, only: AlsoGUGA_SuperGroupIndexer_t

    better_implicit_none

    private

    public :: GAS_SuperGroupIndexer_t, get_CAS_indexer
    public :: get_first_supergroup, get_last_supergroup, get_supergroups, next_supergroup

    type, extends(AlsoGUGA_SuperGroupIndexer_t) :: GAS_SuperGroupIndexer_t
        class(GASSpec_t), allocatable :: GAS_spec
    contains
        private
        procedure, public :: get_supergroups => indexer_get_supergroups
            !! We overwrite the `get_supergroups` method of the generalized indexer
            !! because this allows us a much larger number of GAS spaces
            !! for disconnected spaces.

        procedure :: to_supergroup_nI
        procedure :: to_supergroup_csf

        procedure :: is_allowed_single
        procedure :: is_allowed_double
        procedure :: is_allowed_distinct_double

        procedure, public :: gen_all_excits => idxer_gen_all_excits
    end type

    interface GAS_SuperGroupIndexer_t
        module procedure construct_SuperGroupIndexer_t
    end interface

contains

    pure function get_supergroups(GAS_spec) result(res)
        !! Get the ordered compositions of n into k summands
        !!  constrained by cumulative minima and maxima.
        !!
        !! GAS allowed compositions are called supergroups.
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, allocatable :: res(:, :)
        character(*), parameter :: this_routine = 'get_supergroups'

        integer :: start_comp(GAS_spec%nGAS()), &
            end_comp(GAS_spec%nGAS()), &
            sg(GAS_spec%nGAS())
        integer(int64) :: start_idx, end_idx
        type(buffer_int_2D_t) :: supergroups

        start_comp = get_first_supergroup(GAS_spec)
        if (.not. GAS_spec%is_connected()) then
            res = reshape(start_comp, [size(start_comp), 1])
        else
            end_comp = get_last_supergroup(GAS_spec)
            start_idx = composition_idx(start_comp)
            end_idx = composition_idx(end_comp)

            call supergroups%init(rows=GAS_spec%nGAS())
            sg = start_comp
            do while (sg(1) /= -1)
                call supergroups%push_back(sg)
                sg = next_supergroup(GAS_spec, end_idx, sg)
            end do
            call supergroups%dump_reset(res)
            @:pure_ASSERT(all(res(:, 1) == start_comp))
            @:pure_ASSERT(all(res(:, size(res, 2)) == end_comp))
        end if
    end function


    pure function get_first_supergroup(GAS_spec) result(res)
        !! Return the first supergroup
        !!
        !! "First" as defined by lexicographically decreasing order.
        class(GASSpec_t), intent(in) :: GAS_spec
            !! GAS constraints
        integer :: res(GAS_spec%nGAS())

        character(*), parameter :: this_routine = 'get_lowest_supergroup'
        integer :: remaining, iGAS, to_add

        select type(GAS_spec)
        type is(LocalGASSpec_t)
            res = GAS_spec%get_min()
            remaining = GAS_spec%nEl() - sum(res)

            iGAS = 1
            do while (remaining > 0)
                to_add = min(GAS_spec%get_max(iGAS) - res(iGAS), &
                             GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS + 1
            end do
            @:pure_ASSERT(remaining == 0)
        type is(CumulGASSpec_t)
            res = GAS_spec%get_cmin() &
                    - eoshift(GAS_spec%get_cmax(), shift=-1)
            remaining = GAS_spec%nEl() - sum(res)
            iGAS = 1
            do while (remaining > 0)
                to_add = min(GAS_spec%get_cmax(iGAS) - sum(res(:iGAS)), &
                             GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS + 1
            end do
            @:pure_ASSERT(remaining == 0)
        type is(FlexibleGASSpec_t)
            res = GAS_spec%supergroups(:, 1)
        class default
            call stop_all(this_routine, 'Wrong class.')
        end select
    end function


    pure function get_last_supergroup(GAS_spec) result(res)
        !! Return the last supergroup
        !!
        !! "Last" as defined by lexicographically decreasing order.
        class(GASSpec_t), intent(in) :: GAS_spec
            !! GAS constraints
        integer :: res(GAS_spec%nGAS())

        character(*), parameter :: this_routine = 'get_highest_supergroup'
        integer :: remaining, iGAS, to_add

        select type(GAS_spec)
        type is(LocalGASSpec_t)
            res = GAS_spec%get_min()
            remaining = GAS_spec%nEl() - sum(res)

            iGAS = GAS_spec%nGAS()
            do while (remaining > 0)
                to_add = min(GAS_spec%get_max(iGAS) - res(iGAS), &
                             GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS - 1
            end do
            @:pure_ASSERT(remaining == 0)
        type is(CumulGASSpec_t)
            res = min(GAS_spec%get_cmin() - eoshift(GAS_spec%get_cmin(), shift=-1), &
                      GAS_spec%GAS_size())
            remaining = GAS_spec%nEl() - sum(res)

            iGAS = GAS_spec%nGAS()
            do while (remaining > 0)
                to_add = min(GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS - 1
            end do
            @:pure_ASSERT(remaining == 0)
        type is(FlexibleGASSpec_t)
            res = GAS_spec%supergroups(:, size(GAS_spec%supergroups, 2))
        class default
            call stop_all(this_routine, 'Wrong class.')
        end select
    end function


    pure function next_supergroup(GAS_spec, comp_idx_last, previous) result(res)
        !! Return the next supergoup
        class(GASSpec_t), intent(in) :: GAS_spec
        integer(int64), intent(in) :: comp_idx_last
            !! The composition index of the last supergroup that can be generated.
            !! "Last" as defined by lexicographically decreasing order.
        integer, intent(in) :: previous(:)
            !! The previous supergroup.
        routine_name("next_supergroup")
        integer :: res(size(previous))
        integer :: k, n, src, tgt
        integer(int64) :: comp_idx

        k = size(previous)
        n = sum(previous)
        if (k == 0) then
            return
        else if (previous(1) == -1) then
            res = -1
            return
        end if
        comp_idx = composition_idx(previous)
        if (comp_idx >= comp_idx_last) then
            res(:) = -1
        else
            @:pure_ASSERT(GAS_spec%contains_supergroup(previous))
            select type(GAS_spec)
            type is(LocalGASSpec_t)
                call find_flip_local(GAS_spec, previous, src, tgt)
                res = move_particles_local(GAS_spec, previous, src, tgt)
            type is(CumulGASSpec_t)
                call find_flip_cumul(GAS_spec, previous, src, tgt)
                res = move_particles_cumul(GAS_spec, previous, src, tgt)
            type is (FlexibleGASSpec_t)
            block
                integer :: i_sg_previous
                integer(int64) :: new_comp_idx
                i_sg_previous = custom_findloc(GAS_spec%composition_indices, composition_idx(previous))
                new_comp_idx = GAS_spec%composition_indices(i_sg_previous + 1)
                res = composition_from_idx(GAS_spec%nGAS(), GAS_spec%nEl(), new_comp_idx)
            end block
            class default
                call stop_all(this_routine, 'Wrong class.')
            end select
        end if

        contains
    end function


    pure subroutine find_flip_local(GAS_spec, sg, src, tgt)
        !! Find source and target to transfer particle
        !!
        !! Find the largest `src`, and the smallest `tgt`,
        !!  with `src <= tgt`,
        !! to transfer one particle from `src -> tgt`
        !! while adherring to GAS constraints.
        !! This is **not yet** the next supergroup, because
        !! one still has to move as many particles as possible
        !! from the right of `tgt` closer to `tgt`.
        type(LocalGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: sg(:)
            !! A given supergroup
        integer, intent(out) :: src, tgt
            !! Source and target GAS space.
            !! If no possible movement with `src <= tgt` exists,
            !! then `tgt == 0`.

        integer, allocatable :: sources(:), targets(:)
        integer :: i_src, i_tgt, i, idx(size(sg))

        idx(:) = [(i, i = 1, size(sg))]

        sources = pack(idx, mask=sg - GAS_spec%get_min() > 0)
        targets = pack(idx, mask=GAS_spec%get_max() - sg > 0)
        tgt = 0
        do i_src = size(sources), 1, -1
            src = sources(i_src)
            i_tgt = custom_findloc(targets - src > 0, .true.)
            if (i_tgt /= 0) then
                tgt = targets(i_tgt)
                return
            end if
        end do
    end subroutine


    pure function move_particles_local(GAS_spec, previous, src, tgt) result(res)
        !! Move particles from `src -> tgt` and as many particles
        !!  as possible from the right of `tgt`.
        type(LocalGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: previous(:)
            !! The previous supergroup
        integer, intent(in) :: src, tgt
            !! Source and target GAS space
        integer :: res(size(previous))
        debug_function_name("move_particles_local")

        if (tgt == 0) then
            res(:) = -1
            return
        end if
        @:pure_ASSERT(src < tgt)
        res = previous
        ! Transfer 1 from src to target
        res(src) = res(src) - 1
        res(tgt) = res(tgt) + 1

        ! Transfer everything from all right neighbours to tgt, if possible
        from_right_to_left: block
            integer :: n_open(size(res)), n_available(size(res))
            integer :: n_move, from, to, idx(size(previous)), i


            idx(:) = [(i, i = 1, size(previous))]
            n_open = GAS_spec%get_max() - res(:)
            n_available = res(:) - GAS_spec%get_min()
            to = tgt
            from = size(res)
            do while (to < from)
                if (n_open(to) == 0) then
                    to = to + 1
                else if (n_available(from) == 0) then
                    from = from - 1
                else
                    n_move = min(n_open(to), n_available(from))
                    res(from) = res(from) - n_move
                    res(to) = res(to) + n_move
                    n_open(to) = n_open(to) - n_move
                    n_available(from) = n_available(from) - n_move
                end if
            end do
        end block from_right_to_left
    end function


    pure subroutine find_flip_cumul(GAS_spec, sg, src, tgt)
        !! Find source and target to transfer particle
        !!
        !! Find the largest `src`, and the smallest `tgt`,
        !!  with `src <= tgt`,
        !! to transfer one particle from `src -> tgt`
        !! while adherring to GAS constraints.
        !! This is **not yet** the next supergroup, because
        !! one still has to move as many particles as possible
        !! from the right of `tgt` closer to `tgt`.
        type(CumulGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: sg(:)
            !! A given supergroup
        integer, intent(out) :: src, tgt
            !! Source and target GAS space.
            !! If no possible movement with `src <= tgt` exists,
            !! then `tgt == 0`.

        integer :: k, i, idx(size(sg))

        idx(:) = [(i, i = 1, size(sg))]
        k = size(sg)

        src = custom_findloc(&
            cumsum(sg(: k - 1)) - GAS_spec%get_cmin(idx(: k - 1)) > 0 &
            .and. sg(: k - 1) > 0, &
            .true., back=.true.)
        ! Also account for Pauli
        tgt = custom_findloc(GAS_spec%GAS_size(idx(src + 1 :)) - sg(src + 1 :) > 0, .true.)
        if (tgt /= 0) tgt = tgt + src
    end subroutine


    pure function move_particles_cumul(GAS_spec, previous, src, tgt) result(res)
        !! Move particles from `src -> tgt` and as many particles
        !!  as possible from the right of `tgt`.
        type(CumulGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: previous(:)
            !! The previous supergroup
        integer, intent(in) :: src, tgt
            !! Source and target GAS space
        integer :: res(size(previous))
        debug_function_name("move_particles_cumul")
        @:pure_ASSERT(GAS_spec%is_valid())

        if (tgt == 0) then
            res(:) = -1
            return
        end if
        @:pure_ASSERT(src < tgt)
        res = previous
        ! Transfer 1 from src to target
        res(src) = res(src) - 1
        res(tgt) = res(tgt) + 1

        ! Transfer everything from all right neighbours to tgt, if possible
        from_right_to_left: block
            integer :: n_open(tgt : size(res)), n_available(tgt + 1 : size(res)), &
                n_move, from, to, idx(size(previous)), i, c_res(size(res))


            idx(:) = [(i, i = 1, size(previous))]
            c_res = cumsum(res)
            n_open = min(GAS_spec%get_cmax(idx(tgt : )) - c_res(tgt :), &
                         GAS_spec%GAS_size(idx(tgt : )) - res(tgt:))
            n_available = res(tgt + 1 : )
            to = tgt
            from = size(res)
            do while (to < from)
                if (n_open(to) == 0) then
                    to = to + 1
                else if (n_available(from) == 0) then
                    from = from - 1
                else
                    n_move = min(n_open(to), n_available(from))
                    res(from) = res(from) - n_move
                    res(to) = res(to) + n_move
                    n_open(to : from - 1) = n_open(to : from - 1) - n_move
                    n_available(from) = n_available(from) + n_move
                end if
            end do
        end block from_right_to_left
    end function


    pure function get_allowed_composition_indices(GAS_spec) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer(int64), allocatable :: res(:)
        debug_function_name("get_allowed_composition_indices")

        integer :: sg(GAS_spec%nGAS())
        integer(int64) :: end_idx
        type(buffer_int64_1D_t) :: buf_index
        @:pure_ASSERT(GAS_spec%is_valid())

        sg = get_first_supergroup(GAS_spec)
        end_idx = composition_idx(get_last_supergroup(GAS_spec))

        call buf_index%init()
        do while (sg(1) /= -1)
            call buf_index%push_back(composition_idx(sg))
            sg = next_supergroup(GAS_spec, end_idx, sg)
        end do
        call buf_index%dump_reset(res)
    end function


    pure function construct_SuperGroupIndexer_t(GAS_spec) result(idxer)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_SuperGroupIndexer_t) :: idxer
        integer(int64), allocatable :: allowed_indices(:)
        character(*), parameter :: this_routine = 'construct_SuperGroupIndexer_t'
        @:pure_ASSERT(GAS_spec%is_valid())

        if (GAS_spec%is_connected()) then
            allowed_indices = get_allowed_composition_indices(GAS_spec)
        else
            allowed_indices = [-1_int64]
        end if


        allocate(idxer%GAS_spec, source=GAS_spec)
        allocate(idxer%allowed_composition_indices, source=allowed_indices)
        idxer%n_El_ = GAS_spec%nEl()
        idxer%n_spin_orbs_ = GAS_spec%n_spin_orbs()
        idxer%n_sg_dim_ = GAS_spec%nGAS()
    end function

    pure function indexer_get_supergroups(this) result(res)
        !! Get the ordered compositions of n into k summands
        !!  constrained by cumulative minima and maxima.
        !!
        !! GAS allowed compositions are called supergroups.
        !!
        !! We overwrite the `get_supergroups` method of the generalized indexer
        !! because this allows us a much larger number of GAS spaces
        !! for disconnected spaces.
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        integer, allocatable :: res(:, :)
        integer(int64) :: i
        routine_name("indexer_get_supergroups")

        allocate(res(this%GAS_spec%nGAS(), this%n_supergroups()))
        if (this%GAS_spec%is_connected()) then
            do i = 1_int64, this%n_supergroups()
                res(:, i) = composition_from_idx(&
                    this%GAS_spec%nGAS(), this%GAS_spec%nEl(), this%allowed_composition_indices(i))
            end do
        else
            associate(GAS_spec => this%GAS_spec)
            select type(GAS_spec)
            type is(LocalGASSpec_t)
                res(:, 1) = GAS_spec%get_min()
            type is(CumulGASSpec_t)
                res(:, 1) = GAS_spec%get_cmin() - eoshift(GAS_spec%get_cmin(), shift=-1)
            type is(FlexibleGASSpec_t)
                res(:, 1) = GAS_spec%supergroups(:, 1)
            class default
                call stop_all(this_routine, 'Wrong class.')
            end select
            end associate
        end if
    end function


    logical pure function is_allowed_single(this, exc, supergroup)
        !! Check if a single excitation is allowed.
        !!
        !! Is called once at initialization, so it does not have to be super fast.
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        type(Excite_1_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)

        is_allowed_single = this%GAS_spec%is_allowed(exc, supergroup)
    end function


    logical pure function is_allowed_double(this, exc, supergroup)
        !! Check if a double excitation is allowed.
        !!
        !! Is called once at initialization, so it does not have to be super fast.
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        type(Excite_2_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)

        is_allowed_double = this%GAS_spec%is_allowed(exc, supergroup)
    end function

    logical pure function is_allowed_distinct_double(this, exc, supergroup)
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        type(DistinctDouble_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)

        is_allowed_distinct_double = this%GAS_spec%is_allowed(exc, supergroup)
    end function


    subroutine idxer_gen_all_excits(this, nI, n_excits, det_list, ic)
        !!  Get all excitated determinants from
        !!  det_I that are allowed under GAS constraints.
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
            !!  Starting determinant
        integer, intent(out) :: n_excits
            !!  Number of determinants
        integer(n_int), intent(out), allocatable :: det_list(:,:)
            !!  Allocatable array of determinants in ilut format
        integer, intent(in), optional :: ic
            !!  Optional input for excitation level (ic=1 => singles, ic=2 => doubles)
            !!      If ommited generate all.

        call this%GAS_spec%gen_all_excits(nI, n_excits, det_list, ic)
    end subroutine

    pure function to_supergroup_nI(this, nI) result(supergroup)
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        integer :: supergroup(this%dim_sg())
        supergroup = this%GAS_spec%count_per_GAS(nI)
    end function

    pure function to_supergroup_csf(this, csf_i) result(supergroup)
        class(GAS_SuperGroupIndexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer :: supergroup(this%dim_sg())
        supergroup = this%GAS_spec%count_per_GAS(csf_i)
    end function

    pure function get_CAS_indexer(n_el, n_spat_orbs) result(res)
        integer, intent(in) :: n_el, n_spat_orbs
        type(GAS_SuperGroupIndexer_t) :: res
        res = GAS_SuperGroupIndexer_t(CAS_spec(n_el, n_spat_orbs))
    end function


end module gasci_supergroup_index
