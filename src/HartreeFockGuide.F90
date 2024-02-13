#include "macros.h"
module HartreeFockGuide

    !This module contains the new stuff about hartree fock guide function
    !and the stuff about the hf state approx

    use OneEInts, only: GetTMatEl
    use constants, only: dp, EPS
    use SystemData, only: nel, nBasis, FixedNodeMinThresh, FixedNodeLowestVal
    use CalcData, only: tControlZeroDenom, GutzwillerParam
    use lattice_mod, only: lat
    use error_handling_neci, only: stop_all

    better_implicit_none
    private
    public :: calc_det_ratio, U_matrix

    real(dp), allocatable :: U_matrix(:, :)

contains

    pure subroutine calc_det_ratio(nJ, nI, det_ratio)
        !! Returns the ratio of weights of `nJ / nI` in the HF basis.
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(out) :: det_ratio
        real(dp) :: determinantJ, determinantI

        determinantJ = determinant(U_matrix(nJ(:), 1 : nEl) )
        determinantI = determinant(U_matrix(nI(:), 1 : nEl) )

        ! Fix a minimum value to avoid the 0.0
        if (.not. abs(determinantJ) > 0.0_dp) determinantJ = 10.0_dp**(-15.0_dp)
        if (.not. abs(determinantI) > 0.0_dp) determinantI = 10.0_dp**(-15.0_dp)
        det_ratio = determinantJ / determinantI
    end subroutine calc_det_ratio

    pure function determinant(mat) result(det)
        !! Evaluates the determinant of a Matrix `M`.
        real(dp), intent(in) :: mat(:, :)
        real(dp):: det

        real(dp) :: work_mat(size(mat, 1), size(mat, 2))
        real(dp) :: tmp, c(size(work_mat, dim=1), size(work_mat, dim=2))
        real(dp) :: maxi
        integer:: i, k, l, m, num(size(work_mat, dim=1)), n

        work_mat(:, :) = mat(:, :)

        n = size(work_mat, dim=1)
        det = 1.0_dp
        do k = 1, n
            maxi = work_mat(k, k)
            num(k) = k
            do i = k + 1, n
                if (abs(maxi) < abs(work_mat(i, k))) then
                    maxi = work_mat(i, k)
                    num(k) = i
                end if
            end do
            if (num(k) /= k) then
                do l = k, n
                    tmp = work_mat(k, l)
                    work_mat(k, l) = work_mat(num(k), l)
                    work_mat(num(k), l) = tmp
                end do
                det = -1.0_dp * det
            end if
            do m = k + 1, n
                c(m, k) = work_mat(m, k) / work_mat(k, k)
                do l = k, n
                    work_mat(m, l) = work_mat(m, l) - c(m, k) * work_mat(k, l)
                end do
            end do !There we made matrix triangular!
        end do

        do i = 1, n
            det = det * work_mat(i, i)
        end do
    end function determinant

end module HartreeFockGuide
