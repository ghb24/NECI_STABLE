#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

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
    use util_mod, only: choose, cumsum, binary_search_first_ge, custom_findloc
    use bit_rep_data, only: nIfD
    use gasci, only: GASSpec_t, LocalGASSpec_t, CumulGASSpec_t
    use hash, only: hash_table_lookup
    use growing_buffers, only: buffer_int_2D_t
    use util_mod, only: stop_all, cumsum

    better_implicit_none

    private

    public :: SuperGroupIndexer_t, lookup_supergroup_indexer

    public :: n_compositions, get_compositions, composition_idx
    public :: get_first_supergroup, get_last_supergroup, get_supergroups, next_supergroup

    type :: SuperGroupIndexer_t
        private
        class(GASSpec_t), allocatable :: GASspec
        integer(int64), allocatable :: allowed_composition_indices(:)
        !> The particle number.
        integer :: N
    contains
        private
        procedure, public :: nEl => get_nEl
        procedure, public :: idx_supergroup => get_supergroup_idx
        procedure, public :: idx_nI => get_supergroup_idx_det
        procedure, public :: lookup_supergroup_idx
        procedure, public :: n_supergroups => get_n_supergroups
        procedure, public :: get_supergroups => indexer_get_supergroups
    end type

    interface SuperGroupIndexer_t
        module procedure construct_SuperGroupIndexer_t
    end interface

    type(SuperGroupIndexer_t), pointer :: lookup_supergroup_indexer => null()

contains

    pure function construct_SuperGroupIndexer_t(GASspec, N) result(idxer)
        class(GASSpec_t), intent(in) :: GASspec
        integer, intent(in) :: N
        type(SuperGroupIndexer_t) :: idxer

        idxer%GASspec = GASspec
        idxer%N = N
        idxer%allowed_composition_indices = get_allowed_composition_indices(GASspec, N)
    end function


    elemental function n_compositions(k, n) result(res)
        integer, intent(in) :: k, n
        integer(int64) :: res
        res = choose(n + k - 1, k - 1)
    end function


    !> @brief
    !> Get the ordered compositions of n into k summands.
    !>
    !> @details
    !> Get all possible solutions for the k dimensional hypersurface.
    !> \f[x_1 + ... + x_k = n  \f]
    !> by taking into account the order.
    !> \f[ 1 + 0 = 1 \f] is different from
    !> \f[ 0 + 1 = 1 \f].
    !> The German wikipedia has a nice article
    !> https://de.wikipedia.org/wiki/Partitionsfunktion#Geordnete_Zahlcompositionen
    !>
    !> The compositions are returned in lexicographically decreasing order.
    pure function get_compositions(k, n) result(res)
        integer, intent(in) :: k, n
        integer, allocatable :: res(:, :)
        integer :: i

        allocate(res(k, n_compositions(k, n)))

        res(:, 1) = 0
        res(1, 1) = n

        do i = 2, size(res, 2)
            res(:, i) = next_composition(res(:, i - 1))
        end do
    end function


    pure function next_composition(previous) result(res)
        integer, intent(in) :: previous(:)
        integer :: res(size(previous))
        integer :: k, n, i
        k = size(previous)
        n = sum(previous)
        if (k == 0) then
            continue
        else if (n == previous(k) .or. previous(1) == -1) then
            res(:) = -1
        else
            i = custom_findloc(previous(: k - 1) > 0, .true., back=.true.) + 1
            ! Transfer 1 from left neighbour and everything from all right neighbours to res(i)
            res(: i - 2) = previous(: i - 2)
            res(i - 1) = previous(i - 1) - 1
            res(i) = previous(i) + 1 + sum(previous(i + 1 :))
            res(i + 1 :) = 0
        end if
    end function

    pure function next_supergroup(GAS_spec, comp_idx_last, previous) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer(int64), intent(in) :: comp_idx_last
        integer, intent(in) :: previous(:)
        integer :: res(size(previous))
        integer :: k, n, src, tgt
        integer(int64) :: comp_idx

        debug_function_name("next_supergroup")
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
                ! res = next_composition(previous)
                ! do while (.not. GAS_spec%contains_supergroup(res))
                !     res = next_composition(res)
                ! end do
                call find_flip_cumul(GAS_spec, previous, src, tgt)
                res = move_particles_cumul(GAS_spec, previous, src, tgt)
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
        debug_function_name("move_particles")

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
        debug_function_name("move_particles")

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

    pure function composition_idx(composition) result(idx)
        integer, intent(in) :: composition(:)
        integer(int64) :: idx

        integer :: reminder, i_summand, leading_term

        idx = 1_int64
        i_summand = 1
        reminder = sum(composition)
        do while (reminder /= 0)
            do leading_term = reminder, composition(i_summand) + 1, -1
                idx = idx + n_compositions(size(composition) - i_summand, reminder - leading_term)
            end do
            reminder = reminder - composition(i_summand)
            i_summand = i_summand + 1
        end do
    end function


    pure function get_allowed_composition_indices(GAS_spec, N) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: N

        integer(int64), allocatable :: res(:)
        integer, allocatable :: supergroups(:, :)
        integer :: i

        supergroups = get_supergroups(GAS_spec, N)
        allocate(res(size(supergroups, 2)))
        do i = 1, size(supergroups, 2)
            res(i) = composition_idx(supergroups(:, i))
        end do
    end function


    !> @brief
    !> Get the number of possible supergroups.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    pure function get_n_supergroups(this) result(res)
        class(SuperGroupIndexer_t), intent(in) :: this
        integer(int64) :: res
        res = size(this%allowed_composition_indices, kind=int64)
    end function

    integer elemental function get_nEl(this)
        class(SuperGroupIndexer_t), intent(in) :: this
        get_nEl = this%N
    end function


    !> @brief
    !> Get the ordered compositions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    pure function indexer_get_supergroups(this) result(res)
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, allocatable :: res(:, :)
        res = get_supergroups(this%GASspec, this%N)
    end function

    !> @brief
    !> Get the ordered compositions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    pure function get_supergroups(GAS_spec, N) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: N
        integer, allocatable :: res(:, :)
        character(*), parameter :: this_routine = 'indexer_get_supergroups'

        integer :: start_comp(GAS_spec%nGAS()), &
            end_comp(GAS_spec%nGAS()), &
            sg(GAS_spec%nGAS())
        integer(int64) :: start_idx, end_idx
        type(buffer_int_2D_t) :: supergroups

        start_comp = get_first_supergroup(GAS_spec, N)
        end_comp = get_last_supergroup(GAS_spec, N)
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
    end function

    pure function get_first_supergroup(GAS_spec, N) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
            !! GAS constraints
        integer, intent(in) :: N
            !! Particle number
        integer :: res(GAS_spec%nGAS())

        character(*), parameter :: this_routine = 'get_lowest_supergroup'
        integer :: remaining, iGAS, to_add

        @:pure_ASSERT(sum(GAS_spec%GAS_size()) >= N)

        select type(GAS_spec)
        type is(LocalGASSpec_t)
            res = GAS_spec%get_min()
            remaining = N - sum(res)

            iGAS = 1
            do while (remaining > 0)
                to_add = min(GAS_spec%get_max(iGAS) - res(iGAS), &
                             GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS + 1
            end do
        type is(CumulGASSpec_t)
            res = GAS_spec%get_cmin() &
                    - eoshift(GAS_spec%get_cmax(), shift=-1)
            remaining = N - sum(res)
            iGAS = 1
            do while (remaining > 0)
                to_add = min(GAS_spec%get_cmax(iGAS) - sum(res(:iGAS)), &
                             GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS + 1
            end do
        end select
        @:pure_ASSERT(remaining == 0)
    end function

    pure function get_last_supergroup(GAS_spec, N) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
            !! GAS constraints
        integer, intent(in) :: N
            !! Particle number
        integer :: res(GAS_spec%nGAS())

        character(*), parameter :: this_routine = 'get_highest_supergroup'
        integer :: remaining, iGAS, to_add

        select type(GAS_spec)
        type is(LocalGASSpec_t)
            res = GAS_spec%get_min()
            remaining = N - sum(res)

            iGAS = GAS_spec%nGAS()
            do while (remaining > 0)
                to_add = min(GAS_spec%get_max(iGAS) - res(iGAS), &
                             GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS - 1
            end do
        type is(CumulGASSpec_t)
            res = min(GAS_spec%get_cmin() - eoshift(GAS_spec%get_cmin(), shift=-1), &
                      GAS_spec%GAS_size())
            remaining = N - sum(res)

            iGAS = GAS_spec%nGAS()
            do while (remaining > 0)
                to_add = min(GAS_spec%GAS_size(iGAS) - res(iGAS), remaining)
                res(iGAS) = res(iGAS) + to_add
                remaining = remaining - to_add
                iGAS = iGAS - 1
            end do
        end select
        @:pure_ASSERT(remaining == 0)
    end function

    pure function get_supergroup_idx(this, supergroup) result(idx)
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: supergroup(:)
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx'

        if (this%GASspec%is_connected()) then
            idx = int(binary_search_first_ge(this%allowed_composition_indices, composition_idx(supergroup)))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function


    !>  @brief
    !>  Calculate the supergroup index for a determinant nI
    !>
    !>  @param[in] nI The determinant for which the supergroup index should be calculated.
    pure function get_supergroup_idx_det(this, nI) result(idx)
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx_det'

        @:pure_ASSERT(this%GASspec%contains_conf(nI))
        if (this%GASspec%is_connected()) then
            idx = this%idx_supergroup(this%GASspec%count_per_GAS(nI))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function


    !>  @brief
    !>  Use a precomputed supergroup index from global_det_data.
    !>
    !>  @details
    !>  This function heavily relies on correctly initialized global data
    !>  outside the control of this class.
    !>  Carefully make sure, that global_det_data is correctly initialized.
    !>
    !>  @param[in] idet The index of nI in the FciMCData::CurrentDets array.
    !>  @param[in] nI The determinant for which the supergroup index should be calculated.
    function lookup_supergroup_idx(this, idet, nI) result(idx)
        use global_det_data, only: global_lookup => get_supergroup_idx
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: idet
        integer, intent(in) :: nI(:)
        integer :: idx
        debug_function_name('lookup_supergroup_idx')

        if (this%GASspec%is_connected()) then
            idx = global_lookup(idet)
            ! Assert that looked up and computed index agree.
            @:pure_ASSERT(idx == this%idx_nI(nI))
        else
            idx = 1
        end if
    end function

end module gasci_supergroup_index
