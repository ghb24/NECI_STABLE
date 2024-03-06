#include "macros.h"
module test_util_mod
    use fruit, only: run_test_case, assert_equals, assert_true
    use util_mod, only: stop_all
    use constants, only: int64, dp
    use util_mod, only: choose_i64, ByteSize_t, defined_ByteSize_t, defined_sizes, &
        operator(*), operator(/)
    use indexing_mod, only: fuse_symm_idx, inv_fuse_symm_idx, &
        fuse_idx, inv_fuse_idx, &
        triangle_number, inv_triangle_number
    better_implicit_none
    private
    public :: test_driver

contains

    subroutine test_driver()

        call run_test_case(test_binomial, "test_binomial")

        call run_test_case(test_fuse, "test_fuse")

        call run_test_case(test_byte_sizes, "test_byte_sizes")

    end subroutine

    subroutine test_binomial()
        call assert_equals(2, choose_i64(2, 1))
        call assert_equals(3, choose_i64(3, 2))
        call assert_equals(15, choose_i64(6, 4))
        call assert_equals(720720, choose_i64(66, 4))
        call assert_equals(7219428434016265740_int64, choose_i64(66, 33))

        call assert_equals(9657648, choose_i64(67, 5))

        call assert_equals(4950, choose_i64(100, 2))

        call assert_equals(161700, choose_i64(100, 3))

        call assert_equals(-1, choose_i64(100, 50, signal_overflow=.true.))

    end subroutine


    subroutine test_fuse()
        test_triangle_number: block
            integer(int64) :: i

            do i = 1_int64, 100_int64
                call assert_equals(i, inv_triangle_number(triangle_number(i)))
            end do

        end block test_triangle_number

        symmetric_matrix: block
            integer, parameter :: n_rows = 4, n_cols = n_rows

            idx: block
                integer :: i, j, ij, out_i, out_j

                do i = 1, n_rows
                    do j = 1, n_cols
                        ij = fuse_idx(i, j, n_rows, fortran=.true.)
                        call assert_true(ij == fuse_idx(j, i, n_rows, fortran=.false.))
                        call inv_fuse_idx(ij, n_rows, out_i, out_j, fortran= .true.)
                        call assert_true(i == out_i .and. j == out_j)

                        ij = fuse_idx(i, j, n_rows, fortran=.false.)
                        call inv_fuse_idx(ij, n_rows, out_i, out_j, fortran= .false.)
                        call assert_true(i == out_i .and. j == out_j)
                    end do
                end do
            end block idx


            symm_idx: block
                integer :: i, j, ij, out_i, out_j
                do i = 1, n_rows
                    do j = 1, n_cols
                        ij = fuse_symm_idx(i, j)
                        call inv_fuse_symm_idx(ij, out_i, out_j)
                        if (i <= j) then
                            call assert_true(i == out_i .and. j == out_j)
                        else
                            call assert_true(j == out_i .and. i == out_j)
                        end if
                    end do
                end do
            end block symm_idx
        end block symmetric_matrix

        non_symmetric_matrix: block
            integer, parameter :: n_rows = 10, n_cols = 100

            integer :: i, j

            do i = 1, n_rows
                do j = 1, n_cols
                    fortran: block
                        integer :: ij, out_i, out_j
                        ij = fuse_idx(i, j, n_rows, fortran=.true.)
                        call inv_fuse_idx(ij, n_rows, out_i, out_j, fortran=.true.)
                        call assert_true(i == out_i .and. j == out_j)
                    end block fortran
                    C: block
                        integer :: ij, out_i, out_j
                        ij = fuse_idx(i, j, n_cols, fortran=.false.)
                        call inv_fuse_idx(ij, n_cols, out_i, out_j, fortran=.false.)
                        call assert_true(i == out_i .and. j == out_j)
                    end block C
                end do
            end do
        end block non_symmetric_matrix
    end subroutine


    subroutine test_byte_sizes()

        call assert_true(2 * defined_sizes%GB == 2_int64 * defined_sizes%GB)

        call assert_true(2 * defined_sizes%GiB >= 2_int64 * defined_sizes%GB)

        call assert_true(2 * defined_sizes%GiB == defined_sizes%GiB + defined_sizes%GiB)

        call assert_true(2.1_dp * defined_sizes%GiB == defined_sizes%GiB + 1.1_dp * defined_sizes%GiB)

        call assert_true(defined_sizes%GiB == (2**30_int64 * defined_sizes%GB) / 1e9_dp)

    end subroutine


end module

program test_prog_util_mod
    use test_util_mod, only: test_driver
    use Parallel_neci, only: MPIInit, MPIEnd
    use util_mod, only: stop_all
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count

    implicit none

    integer :: failed_count

    call MPIInit(.false.)
    call init_fruit()
    call test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) call stop_all('test_util_mod', 'failed_tests')

    call MPIEnd(.false.)
end program test_prog_util_mod
