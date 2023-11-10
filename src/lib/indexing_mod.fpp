#include "macros.h"
#:include "../algorithms.fpph"
#:include "../macros.fpph"

module indexing_mod
    use constants, only: dp, int32, int64
    better_implicit_none
    private
    public :: fuse_symm_idx, inv_fuse_symm_idx, fuse_idx, inv_fuse_idx

    interface fuse_symm_idx
        module procedure fuse_symm_idx_int32, fuse_symm_idx_int64
    end interface

    interface inv_fuse_symm_idx
        module procedure inv_fuse_symm_idx_int32, inv_fuse_symm_idx_int64
    end interface

    interface fuse_idx
        module procedure fuse_idx_int32, fuse_idx_int64
    end interface

    interface inv_fuse_idx
        module procedure inv_fuse_idx_int32, inv_fuse_idx_int64
    end interface

contains

    elemental function fuse_symm_idx_int32(q, p) result(pq)
        !! Fuse an index
        !!
        !! create a composite index out of two indices, assuming they are unordered
        !! i.e. their ordering does not matter.
        !!
        !! The ordering looks like this:
        !!
        !!   1   2   4   7
        !!       3   5   8
        !!           6   9
        !!              10
        integer(int32), intent(in) :: p, q
            !! 2d-array indices
        integer(int32) :: pq
            !! 1d-array index assuming the array is symmetric w.r. p<->q

        pq = int(fuse_symm_idx_int64(int(q, int64), int(p, int64)), kind(pq))
    end function


    elemental function fuse_symm_idx_int64(p, q) result(pq)
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
        integer(int64), intent(in) :: p, q
            !! 2d-array indices
        integer(int64) :: pq
            !! 1d-array index assuming the array is symmetric w.r. p<->q

        if (p <= q) then
            pq = p + q * (q - 1) / 2
        else
            pq = q + p * (p - 1) / 2
        endif
    end function


    elemental subroutine inv_fuse_symm_idx_int64(pq, p, q)
        !! Invert the `fuse_symm_idx` routine
        !!
        !! Returns index pair x <= y for the fused index xy.
        !! This assumes that the underlying matrix is symmetric,
        !! i.e. (x, y) and (y, x) map to the same xy.
        integer(int64), intent(in) :: pq
            !! The fused index
        integer(int64), intent(out) :: p, q
            !! The index pair x <= y

        q = ceiling(-0.5_dp + sqrt(0.25_dp + real(2_int64 * pq, kind=dp)), kind=int64)
        p = pq - fuse_symm_idx(1_int64, q) + 1_int64
    end subroutine


    elemental subroutine inv_fuse_symm_idx_int32(pq, p, q)
        !! Invert the `fuse_symm_idx` routine
        !!
        !! Returns index pair x <= y for the fused index xy.
        !! This assumes that the underlying matrix is symmetric,
        !! i.e. (x, y) and (y, x) map to the same xy.
        integer(int32), intent(in) :: pq
            !! The fused index
        integer(int32), intent(out) :: p, q
            !! The index pair x <= y

        integer(int64) :: tmp_p, tmp_q

        call inv_fuse_symm_idx(int(pq, int64), tmp_p, tmp_q)
        p = int(tmp_p, kind(p))
        q = int(tmp_q, kind(q))
    end subroutine

    elemental function fuse_idx_int64(p, q, n_rows, fortran) result(pq)
        !! Fuse `p, q` into one contiguous index
        !!
        !! If `fortran` is true, then this assumes columnwise Fortran order, i.e. the index looks like this:
        !!
        !!    1     5
        !!    2     6
        !!    3     7
        !!    4     8
        !!
        !! Otherwise it is `C`-order, i.e. row mayor (but still one-indexed).
        !!
        !! @note
        !!  Note that this function does preserve lexicographical ordering, if `fortran == .false.`.
        !!  If `(p_1, q_1) <= (p_2, q_2)`
        !!  `fuse_idx(p_1, q_1, n_rows, fortran=.false.) <= fuse_idx(p_2, q_2, fortran=.false.)`.
        !! @endote
        integer(int64), intent(in) :: p, q
            !! 2d-array indices
        integer(int64), intent(in) :: n_rows
            !! number of rows
        logical, intent(in) :: fortran
            !! Use column mayor order
        integer(int64) :: pq
            !! contiguous 1d-array index
        debug_function_name("fuse_idx_int64")

        @:pure_ASSERT(all([p, q, n_rows] > 0))
        @:pure_ASSERT(p <= n_rows)
        if (fortran) then
            pq = p + (q - 1_int64) * n_rows
        else
            pq = q + (p - 1_int64) * n_rows
        end if
    end function

    elemental function fuse_idx_int32(p, q, n_rows, fortran) result(pq)
        !! Fuse `p, q` into one contiguous index
        !!
        !! This assumes columnwise Fortran order, i.e. the index looks like this:
        !!
        !!    1     5
        !!    2     6
        !!    3     7
        !!    4     8
        !!
        !! Otherwise it is `C`-order, i.e. row mayor.
        !!
        !! @note
        !!  Note that this function does preserve lexicographical ordering, if `fortran == .false.`.
        !!  If `(p_1, q_1) <= (p_2, q_2)`
        !!  `fuse_idx(p_1, q_1, n_rows, fortran=.false.) <= fuse_idx(p_2, q_2, fortran=.false.)`.
        !! @endote
        integer(int32), intent(in) :: p, q
            !! 2d-array indices
        integer(int32), intent(in) :: n_rows
            !! number of rows
        logical, intent(in) :: fortran
            !! Use column mayor order
        integer(int32) :: pq
            !! contiguous 1d-array index

        pq = int(fuse_idx(int(p, int64), int(q, int64), int(n_rows, int64), fortran), kind(pq))
    end function

    elemental subroutine inv_fuse_idx_int64(pq, n_rows, p, q, fortran)
        !! Fuse `p, q` into one contiguous index
        !!
        !! This assumes columnwise Fortran order, i.e. the index looks like this:
        !!
        !!    1     5
        !!    2     6
        !!    3     7
        !!    4     8
        integer(int64), intent(in) :: n_rows
            !! number of rows
        integer(int64), intent(in) :: pq
            !! contiguous 1d-array index
        logical, intent(in) :: fortran
            !! Use column mayor order
        integer(int64), intent(out) :: p, q
            !! 2d-array indices

        if (fortran) then
            q = (pq - 1_int64) / n_rows + 1_int64
            p = pq - (q - 1_int64) * n_rows
        else
            p = (pq - 1_int64) / n_rows + 1_int64
            q = pq - (p - 1_int64) * n_rows
        end if
    end subroutine

    elemental subroutine inv_fuse_idx_int32(pq, n_rows, p, q, fortran)
        !! Fuse `p, q` into one contiguous index
        !!
        !! This assumes columnwise Fortran order, i.e. the index looks like this:
        !!
        !!    1     5
        !!    2     6
        !!    3     7
        !!    4     8
        integer(int32), intent(in) :: n_rows
            !! number of rows
        integer(int32), intent(in) :: pq
            !! contiguous 1d-array index
        logical, intent(in) :: fortran
            !! Use column mayor order
        integer(int32), intent(out) :: p, q
            !! 2d-array indices

        integer(int64) :: large_p, large_q

        call inv_fuse_idx(int(pq, int64), int(n_rows, int64), large_p, large_q, fortran)
        p = int(large_p, kind(p))
        q = int(large_q, kind(q))
    end subroutine

end module indexing_mod
