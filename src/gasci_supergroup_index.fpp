#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_supergroup_index
    use constants, only: int64
    use util_mod, only: choose, cumsum, binary_search_first_ge
    use gasci, only: GASSpec_t

    implicit none

    private

    public :: n_partitions, get_partitions, partition_idx
    public :: n_supergroups, get_supergroups, supergroup_idx

    public :: supergroup_idx_precomputed, get_supergroup_indices
    public :: SuperGroupIndexer_t

    type :: SuperGroupIndexer_t
        private
        type(GASSpec_t) :: GASspec
        integer(int64), allocatable :: supergroup_idx(:)
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
        integer(int64) :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx'

        idx = binary_search_first_ge(self%supergroup_idx, partition_idx(supergroup))
        @:pure_ASSERT(idx /= -1)
    end function

    pure function get_supergroup_idx_det(self, nI) result(idx)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer, intent(in) :: nI(:)
        integer(int64) :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx_det'

        @:pure_ASSERT(self%GASspec%contains_det(nI))

        idx = binary_search_first_ge(&
                    self%supergroup_idx, &
                    partition_idx(self%GASspec%count_per_GAS(nI)))
        @:pure_ASSERT(idx /= -1)
    end function

    pure function construct_SuperGroupIndexer_t(GASspec) result(idxer)
        type(GASSpec_t), intent(in) :: GASspec
        type(SuperGroupIndexer_t) :: idxer

        integer :: i

        idxer%GASspec = GASspec

        idxer%supergroup_idx = get_supergroup_indices(&
                GASspec%cumulated_min([(i, i = 1, GASspec%nGAS())]) , &
                GASspec%cumulated_max([(i, i = 1, GASspec%nGAS())]))
    end function

    elemental function n_partitions(k, n) result(res)
        integer, intent(in) :: k, n
        integer(int64) :: res
        res = choose(n + k - 1, k - 1)
    end function

    !> @brief
    !> Get the ordered partitions of n into k summands.
    !>
    !> @details
    !> Get all possible solutions for the k dimensional hypersurface.
    !> \f[x_1 + ... + x_k = n  \f]
    !> by taking into account the order.
    !> \f[ 1 + 0 = 1 \f] is different from
    !> \f[ 0 + 1 = 1 \f].
    !> The German wikipedia has a nice article
    !> https://de.wikipedia.org/wiki/Partitionsfunktion#Geordnete_Zahlpartitionen
    pure function get_partitions(k, n) result(res)
        integer, intent(in) :: k, n
        integer :: res(k, n_partitions(k, n))
        integer :: idx_part, j, i

        idx_part = 1
        res(:, idx_part) = 0
        res(1, idx_part) = n

        if (k == 1) return

        do idx_part = 2, size(res, 2) - 1
            res(:, idx_part) = res(:, idx_part - 1)

            do j = size(res, 1), 2, -1
                if (res(j - 1, idx_part) > 0) exit
            end do

            ! Transfer 1 from left neighbour and everything from all right neighbours to res(j)
            res(j, idx_part) = res(j, idx_part) + 1 + sum(res(j + 1 :, idx_part))
            res(j + 1 :, idx_part) = 0
            res(j - 1, idx_part) = res(j - 1, idx_part) - 1
        end do

        res(:, idx_part) = 0
        res(size(res, 1), idx_part) = n
    end function

    pure function partition_idx(partition) result(idx)
        integer, intent(in) :: partition(:)
        integer(int64) :: idx
        character(*), parameter :: this_routine = 'get_partition_idx'

        integer :: reminder, i_summand, leading_term

        idx = 1_int64
        i_summand = 1
        reminder = sum(partition)
        do while (reminder /= 0)
            do leading_term = reminder, partition(i_summand) + 1, -1
                idx = idx + n_partitions(size(partition) - i_summand, reminder - leading_term)
            end do
            reminder = reminder - partition(i_summand)
            i_summand = i_summand + 1
        end do
    end function




    !> @brief
    !> Get the ordered partitions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    !> GAS allowed partitions are called supergroups.
    pure function get_supergroups(cn_min, cn_max) result(res)
        integer, intent(in) :: cn_min(:), cn_max(:)
        integer, allocatable :: res(:, :)
        integer :: k, n

        integer :: i, j
        integer, allocatable :: all_partitions(:, :)

        k = size(cn_min)
        n = cn_min(k)

        all_partitions = get_partitions(k, n)


        allocate(res(size(cn_min), n_supergroups(cn_min, cn_max)))
        j = 1
        do i = 1, size(res, 2)
            do while (any(cn_min > cumsum(all_partitions(:, j)) .or. cumsum(all_partitions(:, j)) > cn_max))
                j = j + 1
            end do
            res(:, i) = all_partitions(:, j)
            j = j + 1
        end do
    end function

    !> @brief
    !> Get the idx of a given supergroup.
    !>
    !> @details
    !> GAS allowed partitions are called supergroups.
    !> Assume lexical decreasing sortedness.
    pure function supergroup_idx(partition, in_cn_min, in_cn_max) result(idx)
        integer, intent(in) :: partition(:), in_cn_min(:), in_cn_max(:)
        integer(int64) :: idx
        character(*), parameter :: this_routine = 'supergroup_idx'

        integer :: reminder, rhs
        integer :: i_summand, leading_term
        integer :: cn_min(size(in_cn_min)), cn_max(size(in_cn_max))

        @:pure_ASSERT(all(in_cn_min <= cumsum(partition) .and. cumsum(partition) <= in_cn_max))

        cn_min = in_cn_min; cn_max = in_cn_max

        idx = 1_int64
        i_summand = 1
        reminder = sum(partition)

        do while (reminder /= 0)
            do leading_term = min(reminder, cn_max(i_summand)), partition(i_summand) + 1, -1
                idx = idx + n_supergroups(cn_min(i_summand + 1 : ) - leading_term, cn_max( i_summand + 1 :) - leading_term)
            end do

            reminder = reminder - partition(i_summand)
            cn_min = cn_min - partition(i_summand)
            cn_max = cn_max - partition(i_summand)
            i_summand = i_summand + 1
        end do
    end function

    pure function supergroup_idx_precomputed(partition, supergroup_indices) result(idx)
        integer, intent(in) :: partition(:)
        integer(int64), intent(in) :: supergroup_indices(:)
        integer :: idx

        idx = binary_search_first_ge(supergroup_indices, partition_idx(partition))
    end function

    pure function get_supergroup_indices(cn_min, cn_max) result(res)
        integer, intent(in) :: cn_min(:), cn_max(:)

        integer(int64), allocatable :: res(:)
        integer, allocatable :: supergroups(:, :)
        integer :: i

        supergroups = get_supergroups(cn_min, cn_max)
        allocate(res(size(supergroups, 2)))
        do i = 1, size(supergroups, 2)
            res(i) = partition_idx(supergroups(:, i))
        end do
    end function

    !> @brief
    !> Get the number of possible supergroups.
    !>
    !> @details
    !> GAS allowed partitions are called supergroups.
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
    !> GAS allowed partitions are called supergroups.
    pure function indexer_n_supergroups(self) result(n_part)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer(int64) :: n_part

        integer :: i, idx(self%GASspec%nGAS())
        idx = [(i, i = 1, size(idx))]
        n_part = n_supergroups(self%GASspec%cumulated_min(idx), &
                               self%GASspec%cumulated_max(idx))
    end function


    !> @brief
    !> Get the ordered partitions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    !> GAS allowed partitions are called supergroups.
    pure function indexer_get_supergroups(self) result(res)
        class(SuperGroupIndexer_t), intent(in) :: self
        integer :: res(self%GASspec%nGAS(), self%n_supergroups())

        integer :: i, idx(self%GASspec%nGAS())
        idx = [(i, i = 1, size(idx))]

        res = get_supergroups(self%GASspec%cumulated_min(idx), &
                              self%GASspec%cumulated_max(idx))
    end function

end module gasci_supergroup_index
