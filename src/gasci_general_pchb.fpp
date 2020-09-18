#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_general_pchb
    use constants, only: n_int, dp, int64, maxExcit
    use util_mod, only: choose, factrl, stop_all, cumsum
    use FciMCData, only: excit_gen_store_type
    use SystemData, only: nEl
    use bit_rep_data, only: NIfTot

    use gasci, only: GAS_specification, GASSpec_t

    implicit none

    private

    public :: gen_general_GASCI_pchb, get_partitions

    public :: get_partition_index, new_get_n_partitions, new_get_partition_index, new_get_partitions, pure_new_get_n_partitions

!     type :: SuperGroupIndexer_t
!         private
!         type(GASSpec_t) :: GAS_spec
!         integer :: n_el
!
!         integer, allocatable :: partitions(:, :)
!         integer :: n_supergroups, start_supergroups
!     contains
!         private
!         procedure :: get_partition_index
!         procedure :: idx => get_supergroup_index
!     end type

contains

    !>  @brief
    !>  The disconnected_GAS_PCHB excitation generator subroutine.
    !>
    !>  @details
    !>  This is a wrapper around `disconnected_GAS_PCHB%gen_excit`
    !>  to match the function pointer interface.
    !>  The interface is common to all excitation generators, see proc_ptrs.F90
    subroutine gen_general_GASCI_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
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
        character(*), parameter :: this_routine = 'gen_GASCI_pchb'

        @:unused_var(exFlag, part_type, store)
        @:ASSERT(GAS_specification%contains(nI))


        @:unused_var(nI, ilutI, nJ, ilutJ, exFlag, ic, ex_mat, tParity, pGen, hel, store, part_type)

    end subroutine gen_general_GASCI_pchb

    pure function get_partition_index(partition) result(idx)
        integer, intent(in) :: partition(:)
        integer :: idx
        character(*), parameter :: this_routine = 'get_partition_index'

        integer :: reminder
        integer :: i_summand, leading_term

        idx = 1
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


    !> @brief
    !> Get the ordered partitions of n into k summands
    !>  constrained by cumulative minima and maxima.
    !>
    !> @details
    pure function new_get_partitions(k, n, cn_min, cn_max) result(res)
        integer, intent(in) :: k, n, cn_min(:), cn_max(:)
        integer :: res(k, pure_new_get_n_partitions(k, n, cn_min, cn_max))

        integer :: i, j, all_partitions(k, n_partitions(k, n))

        all_partitions(:, :) = get_partitions(k, n)

        j = 1
        do i = 1, size(res, 2)
            do while (any(cn_min > cumsum(all_partitions(:, j)) .or. cumsum(all_partitions(:, j)) > cn_max))
                j = j + 1
            end do
            res(:, i) = all_partitions(:, j)
            j = j + 1
        end do
    end function

    elemental function n_partitions(k, n) result(res)
        integer, intent(in) :: k, n
        integer :: res
        res = choose(n + k - 1, k - 1)
    end function

    function new_get_partition_index(partition, in_cn_min, in_cn_max) result(idx)
        integer, intent(in) :: partition(:), in_cn_min(:), in_cn_max(:)
        integer :: idx
        character(*), parameter :: this_routine = 'new_get_partition_index'

        integer :: reminder, rhs
        integer :: i_summand, leading_term
        integer :: cn_min(size(in_cn_min)), cn_max(size(in_cn_max))

        @:pure_ASSERT(all(in_cn_min <= cumsum(partition) .and. cumsum(partition) <= in_cn_max))

        cn_min = in_cn_min; cn_max = in_cn_max

        idx = 1
        i_summand = 1
        reminder = sum(partition)

        do while (reminder /= 0)
            do leading_term = min(reminder, cn_max(i_summand)), partition(i_summand) + 1, -1
                idx = idx + new_get_n_partitions(size(partition) - i_summand, reminder - leading_term, cn_min(i_summand + 1 : ) - leading_term, cn_max( i_summand + 1 :) - leading_term)
            end do

            reminder = reminder - partition(i_summand)
            cn_min = cn_min - partition(i_summand)
            cn_max = cn_max - partition(i_summand)
            i_summand = i_summand + 1
        end do
    end function

    recursive function new_get_n_partitions(k, n, cn_min, cn_max) result(n_part)
        integer, intent(in) :: k, n
        integer, intent(in) :: cn_min(:), cn_max(:)
        integer :: n_part
        integer :: i
        character(*), parameter :: this_routine = 'new_get_n_partitions'

        @:ASSERT(k == size(cn_min) .and. k == size(cn_max), k, size(cn_min), size(cn_max))

        if (k == 1) then
            n_part = merge(1, 0, cn_min(1) <= n .and. n <= cn_max(1))
            return
        else
            n_part = 0
            do i = max(0, cn_min(1)), min(n, cn_max(1))
                n_part = n_part + new_get_n_partitions(k - 1, n - i, cn_min(2:) - i, cn_max(2:) - i)
            end do
            return
        end if
    end function


    recursive pure function pure_new_get_n_partitions(k, n, cn_min, cn_max) result(n_part)
        integer, intent(in) :: k, n
        integer, intent(in) :: cn_min(:), cn_max(:)
        integer :: n_part
        integer :: i
        character(*), parameter :: this_routine = 'pure_new_get_n_partitions'

        @:pure_ASSERT(k == size(cn_min) .and. k == size(cn_max), k, size(cn_min), size(cn_max))

        if (k == 0) then
            n_part = 0
        else if (k == 1) then
            n_part = merge(1, 0, cn_min(1) <= n .and. n <= cn_max(1))
            return
        else
            n_part = 0
            do i = max(0, cn_min(1)), min(n, cn_max(1))
                n_part = n_part + pure_new_get_n_partitions(k - 1, n - i, cn_min(2:) - i, cn_max(2:) - i)
            end do
            return
        end if
    end function
end module gasci_general_pchb
