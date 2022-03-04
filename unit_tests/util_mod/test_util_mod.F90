#include "macros.h"
module test_util_mod
    use fruit
    use constants, only: int64
    use util_mod, only: choose_i64
    better_implicit_none
    private
    public :: test_driver

contains

    subroutine test_driver()

        call run_test_case(test_binomial, "test_binomial")

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


end module

program test_prog_util_mod
    use test_util_mod, only: test_driver

    use fruit

    implicit none

    integer :: failed_count

    call init_fruit()
    call test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) call stop_all('test_util_mod', 'failed_tests')

end program test_prog_util_mod

