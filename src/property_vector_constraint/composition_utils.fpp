#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module composition_utils

    use constants, only: int64
    use util_mod, only: choose_i64, custom_findloc
    use growing_buffers, only: buffer_int_2D_t, buffer_int64_1D_t
    use util_mod, only: stop_all, cumsum

    better_implicit_none

    private
    public :: n_compositions, get_compositions, composition_idx, composition_from_idx, &
        next_composition

contains

    elemental function n_compositions(k, n) result(res)
        !! Return the number of compositions for `k` summands and a sum of `n`
        !!
        !! A composition is a solution to the integer equation
        !! \( n = x_1 + ... + x_k \)
        !! with \( x_i, n \in \mathbb{N}^0, k \in \mathbb{N}\).
        integer, intent(in) :: k, n
        integer(int64) :: res
        res = choose_i64(n + k - 1, k - 1)
    end function


    pure function next_composition(previous) result(res)
        !! Return the next composition.
        !!
        !! If there is no next composition or the first element is -1,
        !! then the result is -1 everywhere.
        !! This means that the "iterator" is exhausted.
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


    pure function get_compositions(k, n) result(res)
        !! Get the ordered compositions of n into k summands.
        !!
        !! Get all possible solutions for the k dimensional hypersurface.
        !! \( x_1 + ... + x_k = n  \)
        !! by taking into account the order.
        !! \( 1 + 0 = 1 \) is different from
        !! \( 0 + 1 = 1 \).
        !! The German wikipedia has a nice article
        !! https://de.wikipedia.org/wiki/Partitionsfunktion#Geordnete_Zahlcompositionen
        !!
        !! The compositions are returned in lexicographically decreasing order.
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


    pure function composition_idx(composition) result(idx)
        !! Return the composition index for a given composition.
        !!
        !! The index is assigned by **lexicographically decreasing** order.
        integer, intent(in) :: composition(:)
        integer(int64) :: idx

        integer :: remaining, i_summand, leading_term

        idx = 1_int64
        i_summand = 1
        remaining = sum(composition)
        do while (remaining /= 0)
            do leading_term = remaining, composition(i_summand) + 1, -1
                idx = idx + n_compositions(size(composition) - i_summand, remaining - leading_term)
            end do
            remaining = remaining - composition(i_summand)
            i_summand = i_summand + 1
        end do
    end function


    pure function composition_from_idx(k, N, idx) result(composition)
        !! Return the composition for a given composition index
        !!
        !! The index is assigned by **lexicographically decreasing** order.
        !! This function is the inverse of `composition_idx`.
        integer, intent(in) :: k
            !! `k` is the number of summands (`k == size(composition)`).
        integer, intent(in) :: N
            !! `N` is the sum over the composition (`N == sum(composition)`).
        integer(int64), intent(in) :: idx
            !! The composition index.
        integer :: composition(k)

        integer(int64) :: new_idx, prev_idx
        integer :: remaining, i_summand, leading_term

        composition(:) = -1
        remaining = N
        new_idx = 1_int64
        do i_summand = 1, k - 1
            loop_leading_term: do leading_term = remaining, 0, -1
                prev_idx = new_idx
                new_idx = new_idx + n_compositions(k - i_summand, remaining - leading_term)
                if (new_idx > idx) then
                    new_idx = prev_idx
                    composition(i_summand) = leading_term
                    remaining = remaining - leading_term
                    exit loop_leading_term
                end if
            end do loop_leading_term
        end do
        composition(k) = remaining
    end function

end module composition_utils
