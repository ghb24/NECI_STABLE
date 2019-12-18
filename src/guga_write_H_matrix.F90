#include "macros.h"

module guga_write_H_matrix
    use constants, only: n_int, dp
    use display_matrices, only: write_matrix
    use unit_test_helpers, only: print_matrix

    use guga_data, only: ExcitationInformation_t
    use guga_excitations, only: calc_guga_matrix_element
    use guga_bitRepOps, onlY: write_det_guga
    use bit_reps, only: nifguga
    implicit none
    private
    public :: write_H_mat

contains

    function get_H_mat(ilutG) result(H_mat)
        integer(n_int), intent(in) :: ilutG(:, :)
        HElement_t(dp), allocatable :: H_mat(:, :)
        integer :: i, j
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_H_mat"
#endif

        ASSERT(lbound(ilutG, 1) == 0 .and. ubound(ilutG, 1) == nIfGUGA)

        allocate(H_mat(size(ilutG, 2), size(ilutG, 2)))
        do j = 1, size(H_mat, 2)
            do i = j, size(H_mat, 1)
                H_mat(i, j) = calc_mat_ele(ilutG(:, i), ilutG(:, j))
                H_mat(j, i) = H_mat(i, j)
            end do
        end do

    contains
        !> Calculates < ilutG_i | H | ilutG_j >
        function calc_mat_ele(ilutG_i, ilutG_j) result(res)
            integer(n_int), intent(in) :: ilutG_i(:), ilutG_j(:)
            HElement_t(dp) :: res
            type(ExcitationInformation_t) :: excit_info
            call calc_guga_matrix_element( &
                ilutG_i, ilutG_j, excit_info, res, &
                t_hamil=.true., calc_type=2)
        end function
    end function

    subroutine write_H_mat(ilutG, path)
        integer(n_int), intent(in) :: ilutG(:, :)
        character(*), intent(in) :: path

        HElement_t(dp), allocatable :: H_mat(:, :)
        integer :: file_id
#ifdef DEBUG_
        character(*), parameter :: this_routine = "write_H_mat"
#endif

        ASSERT(lbound(ilutG, 1) == 0 .and. ubound(ilutG, 2) == nIfGUGA)

        open(file_id, file=path)
            call write_header(unit_id=file_id)
            call write_CSF_repr(ilutG, unit_id=file_id)
            call print_matrix(get_H_mat(ilutG), file_id)
        close(file_id)

    contains

        subroutine write_header(unit_id)
            integer, intent(in) :: unit_id
            write(unit_id, '(A)') '# This file contains the CSF strings &
                &and the FCI-matrix H = < i | H | j >.'
            write(unit_id, '(A)') '# The H-matrix is sorted according to &
                &the order of the CSF strings.'
        end subroutine

        subroutine write_CSF_repr(ilutG, unit_id)
            integer(n_int), intent(in) :: ilutG(:, :)
            integer, intent(in) :: unit_id

            integer :: i
#ifdef DEBUG_
            character(*), parameter :: this_routine = "write_CSF_repr"
#endif

            ASSERT(lbound(ilutG, 1) == 0 .and. ubound(ilutG, 2) == nIfGUGA)

            do i = 1, size(ilutG, 2)
                write(unit_id, '(A2, I15, A2)', advance='no') '# ', i, '. '
                call write_det_guga(ilut=ilutG(:, i), nunit=unit_id)
            end do
        end subroutine
    end subroutine

end module guga_write_H_matrix
