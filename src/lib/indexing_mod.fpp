#include "macros.h"
#:include "../algorithms.fpph"
#:include "../macros.fpph"

#:set kinds = ['int32', 'int64']
module indexing_mod
    use constants, only: dp, int32, int64
    use util_mod, only: operator(.implies.), operator(.div.)
    better_implicit_none
    private
    public :: fuse_symm_idx, inv_fuse_symm_idx, fuse_idx, inv_fuse_idx, &
        new_fuse_symm_idx, triangle_number, inv_triangle_number

    interface fuse_symm_idx
#:for kind in kinds
        module procedure fuse_symm_idx_${kind}$
#:endfor
    end interface

    interface inv_fuse_symm_idx
#:for kind in kinds
        module procedure inv_fuse_symm_idx_${kind}$
#:endfor
    end interface

    interface fuse_idx
#:for kind in kinds
        module procedure fuse_idx_${kind}$
#:endfor
    end interface

    interface inv_fuse_idx
#:for kind in kinds
        module procedure inv_fuse_idx_${kind}$
#:endfor
    end interface

    interface new_fuse_symm_idx
#:for kind in kinds
        module procedure new_fuse_symm_idx_${kind}$
#:endfor
    end interface

    interface triangle_number
#:for kind in kinds
        module procedure triangle_number_${kind}$
#:endfor
    end interface

    interface inv_triangle_number
#:for kind in kinds
        module procedure inv_triangle_number_${kind}$
#:endfor
    end interface
contains

#:for kind in kinds
    elemental function triangle_number_${kind}$(n) result(res)
        integer(${kind}$), intent(in) :: n
        integer(${kind}$) :: res
        res = n * (n + 1) / 2
    end function

    elemental function inv_triangle_number_${kind}$(n) result(res)
        integer(${kind}$), intent(in) :: n
        integer(${kind}$) :: res
        routine_name("inv_triangle_number_${kind}$")
        res = nint(-0.5_dp + sqrt(0.25_dp + 2._dp * n), ${kind}$)
        @:pure_ASSERT(triangle_number(res) == n)
    end function

    elemental function fuse_symm_idx_${kind}$(p, q) result(pq)
        !! Fuse an index
        !!
        !! create a composite index out of two indices, assuming they are unordered
        !! i.e. their ordering does not matter.
        !!
        !! @note
        !!  Note that this function does not preserve order if we assume lexicographical ordering.
        !!  If `(p_1, q_1) <= (p_2, q_2)`, then we don't have necessarily
        !!  `fuse_sym_idx(p_1, q_1) <= fuse_symm_idx(p_2, q_2)`.
        !! @endote
        !!
        !! The ordering looks like this:
        !!
        !!   1   2   4   7
        !!       3   5   8
        !!           6   9
        !!              10
        integer(${kind}$), intent(in) :: p, q
            !! 2d-array indices
        integer(${kind}$) :: pq
            !! 1d-array index assuming the array is symmetric w.r. p<->q

        if (p <= q) then
            pq = p + q * (q - 1) / 2
        else
            pq = q + p * (p - 1) / 2
        endif
    end function

    elemental function new_fuse_symm_idx_${kind}$(p, q, n_dim) result(pq)
        !! Fuse an index
        !!
        !! create a composite index out of two indices, assuming they are unordered
        !! i.e. their ordering does not matter.
        !!
        !! @note
        !!  Note that this function does not preserve order if we assume lexicographical ordering.
        !!  If `(p_1, q_1) <= (p_2, q_2)`, then we don't have necessarily
        !!  `fuse_sym_idx(p_1, q_1) <= fuse_symm_idx(p_2, q_2)`.
        !! @endote
        !!
        !! The ordering looks like this:
        !!
        !!   1   2   3   4
        !!       5   6   7
        !!           8   9
        !!              10
        integer(${kind}$), intent(in) :: p, q, n_dim
            !! 2d-array indices
        integer(${kind}$) :: pq
            !! 1d-array index assuming the array is symmetric w.r. p<->q

        if (p <= q) then
            pq = (q - (p - 1)) + triangle_number(n_dim) - triangle_number(n_dim - p + 1)
        else
            pq = (p - (q - 1)) + triangle_number(n_dim) - triangle_number(n_dim - q + 1)
        endif
    end function

    elemental subroutine inv_fuse_symm_idx_${kind}$(pq, p, q)
        !! Invert the `fuse_symm_idx` routine
        !!
        !! Returns index pair x <= y for the fused index xy.
        !! This assumes that the underlying matrix is symmetric,
        !! i.e. (x, y) and (y, x) map to the same xy.
        integer(${kind}$), intent(in) :: pq
            !! The fused index
        integer(${kind}$), intent(out) :: p, q
            !! The index pair x <= y

        q = ceiling(-0.5_dp + sqrt(0.25_dp + 2._dp * pq), kind=${kind}$)
        p = pq - fuse_symm_idx(1_${kind}$, q) + 1_${kind}$
    end subroutine

    elemental function fuse_idx_${kind}$(p, q, n_dim, fortran) result(pq)
        !! Fuse `p, q` into one contiguous index
        !!
        !! If `fortran` is true, then this assumes columnwise Fortran order, i.e. the index looks like this:
        !!
        !!    1     5
        !!    2     6
        !!    3     7
        !!    4     8
        !!
        !! Otherwise it is `C`-order, i.e. row major (but still one-indexed).
        !!
        !! @note
        !!  Note that this function does preserve lexicographical ordering, if `fortran == .false.`.
        !!  If `(p_1, q_1) <= (p_2, q_2)`
        !!  `fuse_idx(p_1, q_1, n_dim, fortran=.false.) <= fuse_idx(p_2, q_2, fortran=.false.)`.
        !! @endote
        integer(${kind}$), intent(in) :: p, q
            !! 2d-array indices
        integer(${kind}$), intent(in) :: n_dim
            !! number of rows for Fortran (column-major), number of columns for
            !!  C (row-major).
        logical, intent(in) :: fortran
            !! Use column major order
        integer(${kind}$) :: pq
            !! contiguous 1d-array index
        debug_function_name("fuse_idx_${kind}$")

        @:pure_ASSERT(all([p, q, n_dim] > 0))
        @:pure_ASSERT(fortran .implies. p <= n_dim)
        @:pure_ASSERT((.not. fortran) .implies. q <= n_dim)
        if (fortran) then
            pq = p + (q - 1_${kind}$) * n_dim
        else
            pq = q + (p - 1_${kind}$) * n_dim
        end if
    end function

    elemental subroutine inv_fuse_idx_${kind}$(pq, n_dim, p, q, fortran)
        !! Fuse `p, q` into one contiguous index
        !!
        !! This assumes columnwise Fortran order, i.e. the index looks like this:
        !!
        !!    1     5
        !!    2     6
        !!    3     7
        !!    4     8
        integer(${kind}$), intent(in) :: n_dim
            !! number of rows
        integer(${kind}$), intent(in) :: pq
            !! contiguous 1d-array index
        logical, intent(in) :: fortran
            !! Use column major order
        integer(${kind}$), intent(out) :: p, q
            !! 2d-array indices

        if (fortran) then
            q = (pq - 1) / n_dim + 1
            p = pq - (q - 1) * n_dim
        else
            p = (pq - 1) / n_dim + 1
            q = pq - (p - 1) * n_dim
        end if
    end subroutine


#:endfor
end module indexing_mod
