module test_matrix_util_mod
    use fruit, only: assert_true, assert_false, run_test_case
    use constants, only: dp
    use matrix_util, only: gram_schmidt, find_degeneracies, canonicalize
    use util_mod, only: operator(.isclose.)
    use display_matrices, only: write_matrix
    implicit none
    private
    public :: test_driver



contains

    subroutine test_find_degeneracies()
        integer, allocatable :: ind(:, :), pairs(:, :)
        real(dp), allocatable :: e_values(:)

        e_values = [1._dp, 1._dp, 3._dp, 5._dp]

        call find_degeneracies(e_values, ind, pairs)

        ! call write_matrix(e_values)
        ! call write_matrix(ind)
        ! call write_matrix(pairs)

    end subroutine

    subroutine test_canonicalize()
        real(dp), allocatable :: V(:, :), lambda(:)
        real(dp), parameter :: &
            expected_V(3, 3) = transpose(reshape(&
                              [ 1, 0, 0, &
                                0, 1, 0, &
                                0, 0, 1  ], [3, 3])), &
            expected_lambda(3) = [1, 2, 2]

        V = transpose(reshape([ 0, 2, 0, &
                                3, 0,-7, &
                                1, 0, 1  ], [3, 3]))
        lambda = [1.99999999999_dp, 1._dp, 2.000000000001_dp]

        call canonicalize(V, lambda)

        call assert_true(all(lambda .isclose. expected_lambda))
        call assert_true(all(V .isclose. expected_V))

    end subroutine


    subroutine test_conversion()

    end subroutine

    subroutine test_driver()
        call run_test_case(test_canonicalize, "test_if_numbers")
        call run_test_case(test_find_degeneracies, "test_split")
        call run_test_case(test_conversion, "test_conversion")
    end subroutine

end module

program test_matrix_util_program

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use test_matrix_util_mod, only: test_driver

    implicit none
    integer :: failed_count

    block

        call init_fruit()

        call test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_matrix_util_program', 'failed_tests')

    end block

contains

end program
