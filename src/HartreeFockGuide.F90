#include "macros.h"
module HartreeFockGuide

    ! This module contains the routines to calculate the HF guiding wavefunction
    ! for the real-space Hubbard model

    use OneEInts, only: GetTMatEl
    use constants, only: dp, EPS
    use SystemData, only: nel, nBasis, FixedNodeMinThresh, FixedNodeLowestVal
    use CalcData, only: tControlZeroDenom, GutzwillerParam
    use lattice_mod, only: lat

    implicit none

    real(dp), allocatable :: U_matrix(:, :)

    character(len=1024) :: orb_coeff_name

    private
    public :: U_matrix, calc_det_ratio, calc_abinitio_det_ratio

contains
    subroutine calc_det_ratio(nJ, nI, det_ratio)
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(out) :: det_ratio
        real(dp), allocatable :: UI(:, :), UJ(:, :)
        real(dp) :: determinantJ, determinantI

        ! Here we compute the U Matrix in order to access the det ratio
        UI = get_U(nI)
        UJ = get_U(nJ)

        determinantJ = determinant(UJ)
        determinantI = determinant(UI)
        ! Fix a minimum value to avoid the 0.0
        if (.not. abs(determinantJ) > 0.0_dp) determinantJ = 10.0_dp**(-15.0_dp)
        if (.not. abs(determinantI) > 0.0_dp) determinantI = 10.0_dp**(-15.0_dp)
        det_ratio = determinantJ / determinantI
    end subroutine calc_det_ratio

    subroutine calc_abinitio_det_ratio(nJ, nI, orb_coeff, det_ratio)
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(in) :: orb_coeff(nBasis / 2, nBasis / 2)
        real(dp), intent(out) :: det_ratio
        real(dp) :: det_I, det_J

        det_J = abinitio_det(nJ, orb_coeff)
        det_I = abinitio_det(nI, orb_coeff)
        if (tControlZeroDenom .and. (.not. abs(det_I) > 0.0_dp) .and. (.not. abs(det_J) > 0.0_dp)) then
            det_ratio = 1.0_dp
        else if (tControlZeroDenom .and. (.not. abs(det_J) > 0.0_dp)) then
            det_ratio = 10.0_dp**(-15) / det_I
        else
            det_ratio = det_J / det_I
        end if
    end subroutine calc_abinitio_det_ratio

    function get_U(nD) result(U_final)
        integer :: nD(nel)
        real(dp) :: U_final(nel, nel)
        integer :: k, l

        do k = 1, nel
            do l = 1, nel
                U_final(k, l) = U_matrix(nD(k), l)
            end do
        end do
    end function get_U

    function determinant(aa) result(det)
        real(dp), intent(in) :: aa(:, :)
        real(dp) :: tmp, c(size(aa, dim=1), size(aa, dim=2))
        real(dp) :: maxi
        integer:: i, k, l, m, num(size(aa, dim=1)), n
        real(dp) :: det

        n = size(aa, dim=1)
        det = 1.0_dp
        do k = 1, n
            maxi = aa(k, k)
            num(k) = k
            do i = k + 1, n
                if (abs(maxi) < abs(aa(i, k))) then
                    maxi = aa(i, k)
                    num(k) = i
                end if
            end do
            if (num(k) /= k) then
                do l = k, n
                    tmp = aa(k, l)
                    aa(k, l) = aa(num(k), l)
                    aa(num(k), l) = tmp
                end do
                det = -1.0_dp * det
            end if
            do m = k + 1, n
                c(m, k) = aa(m, k) / aa(k, k)
                do l = k, n
                    aa(m, l) = aa(m, l) - c(m, k) * aa(k, l)
                end do
            end do ! There we made matrix triangular!
        end do

        do i = 1, n
            det = det * aa(i, i)
        end do
    end function determinant

end module HartreeFockGuide
