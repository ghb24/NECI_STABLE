#include "macros.h"
module HartreeFockGuide

    !This module contains the new stuff about hartree fock guide function
    !and the stuff about the hf state approx

    use OneEInts, only: GetTMatEl
    use constants, only: dp, EPS
    use SystemData, only: nel, nBasis, FixedNodeMinThresh, FixedNodeLowestVal
    use CalcData, only: tControlZeroDenom, GutzwillerParam
    use lattice_mod, only: lat

    implicit none

    real(dp), allocatable :: U_matrix(:, :)

    character(len=1024) :: orb_coeff_name

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

    subroutine calc_hf_weight(nI, hf_weight)
        integer, intent(in) :: nI(nel)
        real(dp), intent(out) :: hf_weight
        real(dp), allocatable:: U(:, :)

        U = get_U(nI)
        hf_weight = determinant(U)
        if (abs(hf_weight) < FixedNodeMinThresh) hf_weight = FixedNodeLowestVal
    end subroutine calc_hf_weight

    subroutine calc_abinitio_hf_weight(nI, orb_coeff, hf_weight)
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: orb_coeff(nBasis / 2, nBasis / 2)
        real(dp), intent(out) :: hf_weight

        hf_weight = abinitio_det(nI, orb_coeff)
    end subroutine calc_abinitio_hf_weight

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

    subroutine read_orb_coeff(orb_coeff)
        ! Read an orbital coefficient matrix output and
        ! construct an overlap matrix
        ! Currently only supports Molpro coeff file

        ! avoid using dimension
        real(dp), intent(out) :: orb_coeff(nBasis / 2, nBasis / 2)
        real(dp) :: line(nBasis / 2)
        integer :: iunit, i

        open(newunit=iunit, FILE=orb_coeff_name)

        do i = 1, nBasis / 2
            read(iunit, *) line
            orb_coeff(:, i) = line
        end do

        close(iunit)
    end subroutine read_orb_coeff

    subroutine split_spin_orbs(nI, a_arr, b_arr, num_swaps)
        ! Split nI into alpha (odd) and beta (even) arrays
        ! Also count the number of row swaps that necessary to
        ! construct a block diagonal matrix
        integer, intent(in)  :: nI(nel)
        integer, intent(out), allocatable :: a_arr(:), b_arr(:)
        integer, intent(out) :: num_swaps
        integer :: i, b_counter

        allocate(a_arr(0))
        allocate(b_arr(0))

        num_swaps = 0
        b_counter = 0

        do i = 1, size(nI)
            if (mod(nI(i), 2) == 0) then
                b_arr = [b_arr, nI(i)]
                b_counter = b_counter + 1
            else
                a_arr = [a_arr, nI(i)]
                num_swaps = num_swaps + b_counter
            end if
        end do
    end subroutine split_spin_orbs

    function abinitio_det(nI, orb_coeff) result(det)
        ! return abinitio overlap determinant value
        integer :: nI(nel), num_row_swaps, total_num_swaps
        integer, allocatable :: a_arr(:), b_arr(:)
        real(dp) :: det, det_a, det_b, orb_coeff(nBasis / 2, nBasis / 2)

        call split_spin_orbs(nI, a_arr, b_arr, num_row_swaps)

        ! final det = sign_factor * det_alpha * det_beta

        ! Selection rules (valid only for the HF guiding wf)
        ! Even nel: #alpha = #beta
        ! Odd nel: #alpha = #beta + 1
        ! Otherwise det = 0
        if (mod(nel, 2) == 0) then
            if (size(a_arr) == size(b_arr)) then
                det_a = determinant(orb_coeff(1:nel / 2,(a_arr + 1) / 2))
                det_b = determinant(orb_coeff(1:nel / 2, b_arr / 2))
                total_num_swaps = num_row_swaps + nel * (nel - 2) / 8
                det = (1 - 2 * mod(total_num_swaps, 2)) * det_a * det_b
            else
                det = 0.0_dp
            end if
        else
            if (size(a_arr) == size(b_arr) + 1) then
                det_a = determinant(orb_coeff(1:(nel + 1) / 2,(a_arr + 1) / 2))
                det_b = determinant(orb_coeff(1:nel / 2, b_arr / 2))
                total_num_swaps = num_row_swaps + (nel * nel - 1) / 8
                det = (1 - 2 * mod(total_num_swaps, 2)) * det_a * det_b
            else
                det = 0.0_dp
            end if
        end if
    end function abinitio_det

    function determinant(aa) result(det)
        real(dp) :: aa(:, :)
        real(dp) :: tmp, c(size(aa, dim=1), size(aa, dim=2))
        real(dp) :: maxi
        integer:: i, k, l, m, num(size(aa, dim=1)), n
        real(dp):: det

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
            end do !There we made matrix triangular!
        end do

        do i = 1, n
            det = det * aa(i, i)
        end do
    end function determinant

end module HartreeFockGuide
