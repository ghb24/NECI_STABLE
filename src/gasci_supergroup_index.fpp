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
    use constants, only: int64
    use util_mod, only: choose, cumsum, binary_search_first_ge
    use gasci, only: GASSpec_t

    implicit none

    private

    public :: SuperGroupIndexer_t

    public :: n_compositions, get_compositions, composition_idx
    public :: n_supergroups, get_supergroups, supergroup_idx
    public :: supergroup_idx_precomputed, get_allowed_composition_indices

    type :: SuperGroupIndexer_t
        private
        type(GASSpec_t) :: GASspec
        integer(int64), allocatable :: allowed_composition_indices(:)
    contains
        private
        procedure, public :: idx_supergroup => get_supergroup_idx
        procedure, public :: idx_nI => get_supergroup_idx_det
        procedure, public :: n_supergroups => indexer_n_supergroups
        procedure, public :: get_supergroups => indexer_get_supergroups
    end type

    interface SuperGroupIndexer_t
        module procedure construct_SuperGroupIndexer_t
    end interface

contains

    pure function get_supergroup_idx(self, supergroup) result(idx)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer, intent(in) :: supergroup(:)
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx'

        if (self%GASspec%is_connected()) then
            idx = int(binary_search_first_ge(self%allowed_composition_indices, composition_idx(supergroup)))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function


    pure function get_supergroup_idx_det(self, nI) result(idx)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer, intent(in) :: nI(:)
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx_det'

        @:pure_ASSERT(self%GASspec%contains_det(nI))
        if (self%GASspec%is_connected()) then
            idx = self%idx_supergroup(self%GASspec%count_per_GAS(nI))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function


    pure function construct_SuperGroupIndexer_t(GASspec) result(idxer)
        type(GASSpec_t), intent(in) :: GASspec
        type(SuperGroupIndexer_t) :: idxer

        integer :: i

        idxer%GASspec = GASspec
        idxer%allowed_composition_indices = get_allowed_composition_indices(&
                GASspec%cumulated_min([(i, i = 1, GASspec%nGAS())]) , &
                GASspec%cumulated_max([(i, i = 1, GASspec%nGAS())]))
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
    pure function get_compositions(k, n) result(res)
        integer, intent(in) :: k, n
        integer :: res(k, n_compositions(k, n))
        integer :: idx_part, i

        idx_part = 1
        res(:, idx_part) = 0
        res(1, idx_part) = n

        if (k == 1) return

        do idx_part = 2, size(res, 2) - 1
            res(:, idx_part) = res(:, idx_part - 1)

            do i = size(res, 1), 2, -1
                if (res(i - 1, idx_part) > 0) exit
            end do

            ! Transfer 1 from left neighbour and everything from all right neighbours to res(j)
            res(i, idx_part) = res(i, idx_part) + 1 + sum(res(i + 1 :, idx_part))
            res(i + 1 :, idx_part) = 0
            res(i - 1, idx_part) = res(i - 1, idx_part) - 1
        end do

        res(:, idx_part) = 0
        res(size(res, 1), idx_part) = n
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


    !> @brief
    !> Get the ordered compositions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    pure function get_supergroups(cn_min, cn_max) result(res)
        integer, intent(in) :: cn_min(:), cn_max(:)
        integer, allocatable :: res(:, :)

        integer :: i, j, k, n
        integer, allocatable :: all_compositions(:, :)

        k = size(cn_min)
        n = cn_min(k)

        all_compositions = get_compositions(k, n)


        allocate(res(size(cn_min), n_supergroups(cn_min, cn_max)))
        j = 1
        do i = 1, size(res, 2)
            do while (any(cn_min > cumsum(all_compositions(:, j)) .or. cumsum(all_compositions(:, j)) > cn_max))
                j = j + 1
            end do
            res(:, i) = all_compositions(:, j)
            j = j + 1
        end do
    end function


    !> @brief
    !> Get the idx of a given supergroup.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    !> Assume lexical decreasing sortedness.
    pure function supergroup_idx(composition, in_cn_min, in_cn_max) result(idx)
        integer, intent(in) :: composition(:), in_cn_min(:), in_cn_max(:)
        integer(int64) :: idx
        character(*), parameter :: this_routine = 'supergroup_idx'

        integer :: reminder
        integer :: i_summand, leading_term
        integer :: cn_min(size(in_cn_min)), cn_max(size(in_cn_max))

        @:pure_ASSERT(all(in_cn_min <= cumsum(composition) .and. cumsum(composition) <= in_cn_max))

        cn_min = in_cn_min; cn_max = in_cn_max

        idx = 1_int64
        i_summand = 1
        reminder = sum(composition)

        do while (reminder /= 0)
            do leading_term = min(reminder, cn_max(i_summand)), composition(i_summand) + 1, -1
                idx = idx + n_supergroups(cn_min(i_summand + 1 : ) - leading_term, cn_max( i_summand + 1 :) - leading_term)
            end do

            reminder = reminder - composition(i_summand)
            cn_min = cn_min - composition(i_summand)
            cn_max = cn_max - composition(i_summand)
            i_summand = i_summand + 1
        end do
    end function


    pure function supergroup_idx_precomputed(composition, supergroup_indices) result(idx)
        integer, intent(in) :: composition(:)
        integer(int64), intent(in) :: supergroup_indices(:)
        integer :: idx

        idx = binary_search_first_ge(supergroup_indices, composition_idx(composition))
    end function


    pure function get_allowed_composition_indices(cn_min, cn_max) result(res)
        integer, intent(in) :: cn_min(:), cn_max(:)

        integer(int64), allocatable :: res(:)
        integer, allocatable :: supergroups(:, :)
        integer :: i

        supergroups = get_supergroups(cn_min, cn_max)
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
    recursive pure function n_supergroups(cn_min, cn_max) result(n_part)
        integer, intent(in) :: cn_min(:), cn_max(:)
        integer(int64) :: n_part

        integer :: k, n, i
        character(*), parameter :: this_routine = 'n_supergroups'

        @:pure_ASSERT(size(cn_min) == size(cn_max))
        k = size(cn_min)
        @:pure_ASSERT(0 <= cn_max(1) .and. cn_min(k) == cn_max(k))
        n = cn_min(k)
        @:pure_ASSERT(all(cn_min(2:) >= cn_min(: k - 1)) .and. all(cn_max(2:) >= cn_max(: k - 1)))

        if (k == 1 .or. n == 0) then
            n_part = merge(1_int64, 0_int64, cn_min(1) <= n .and. n <= cn_max(1))
        else
            n_part = 0_int64
            do i = max(0, cn_min(1)), min(n, cn_max(1))
                n_part = n_part + n_supergroups(cn_min(2:) - i, cn_max(2:) - i)
            end do
        end if
    end function


    !> @brief
    !> Get the number of possible supergroups.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    pure function indexer_n_supergroups(self) result(n_part)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer(int64) :: n_part

        integer :: i, idx(self%GASspec%nGAS())
        idx = [(i, i = 1, size(idx))]
        n_part = n_supergroups(self%GASspec%cumulated_min(idx), &
                               self%GASspec%cumulated_max(idx))
    end function


    !> @brief
    !> Get the ordered compositions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    !> GAS allowed compositions are called supergroups.
    pure function indexer_get_supergroups(self) result(res)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer :: res(self%GASspec%nGAS(), self%n_supergroups())

        integer :: i, idx(self%GASspec%nGAS())

        idx = [(i, i = 1, size(idx))]
        res = get_supergroups(self%GASspec%cumulated_min(idx), &
                              self%GASspec%cumulated_max(idx))
    end function

end module gasci_supergroup_index
