#include "macros.h"

module test_SD_spin_purification_mod
    use fruit
    use constants, only: dp, n_int
    use SD_spin_purification_mod, only: S2_expval, spin_momentum, get_open_shell, spin_q_num, S2_expval_exc
    use excitation_types, only: Excitation_t, NoExc_t, SingleExc_t, DoubleExc_t, TripleExc_t
    use util_mod, only: operator(.isclose.)
    implicit none
    private
    public :: test_S2_expval, test_spin_momentum, test_spin_q_num, test_get_open_shell, test_S2_expval_exc



contains

    subroutine test_spin_momentum()
        call assert_equals(0._dp, spin_momentum(s=0._dp))

        call assert_true(0.75_dp .isclose. spin_momentum(s=0.5_dp)**2)
    end subroutine

    subroutine test_get_open_shell()
        integer, allocatable :: input(:), expected(:), calculated(:)

        expected = [integer::]
        input = [integer::]
        calculated = get_open_shell(input)
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(expected))

        expected = [integer::]
        input = [1, 2, 3, 4]
        calculated = get_open_shell(input)
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(expected))

        expected = [1]
        input = [1, 3, 4]
        calculated = get_open_shell(input)
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(expected))

        expected = [1, 3, 6]
        input = [1, 3, 6, 7, 8]
        calculated = get_open_shell(input)
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(expected))
    end subroutine

    subroutine test_spin_q_num()
        call assert_true(spin_q_num(0._dp) .isclose. 0._dp)
        call assert_true(spin_q_num(sqrt(2._dp)) .isclose. 1._dp)
        call assert_true(spin_q_num(sqrt(1._dp)) .isclose. -0.5_dp + sqrt(1.25_dp))
    end subroutine

    subroutine test_S2_expval()
        call assert_equals(0._dp, S2_expval([integer::], [integer::]))

        call assert_equals(0._dp, S2_expval([1, 2, 3, 4], [1, 2, 3, 4]))

        call assert_equals(0._dp, S2_expval([1, 3], [5, 7]))

        call assert_true(spin_momentum(0.5_dp)**2 .isclose. S2_expval([1], [1]))

        call assert_true(spin_momentum(0.5_dp)**2 .isclose. S2_expval([2], [2]))

        call assert_true(2.0_dp .isclose. S2_expval([1, 3], [1, 3]))

        call assert_true(spin_momentum(1.0_dp)**2 .isclose. S2_expval([2, 4], [2, 4]))

        call assert_true(1.0_dp .isclose. S2_expval([1, 4], [1, 4]))

        call assert_true(1.0_dp .isclose. S2_expval([1, 4, 5, 6], [2, 3, 5, 6]))
    end subroutine


    subroutine test_S2_expval_exc()
        call assert_equals(0._dp, S2_expval_exc([integer::], NoExc_t()))

        call assert_equals(0._dp, S2_expval_exc([integer::], DoubleExc_t(1, 2, 3, 4)))

        call assert_equals(0._dp, S2_expval_exc([1, 2, 3, 4], NoExc_t()))

        call assert_equals(0._dp, S2_expval_exc([1, 3], DoubleExc_t(1, 3, 5, 7)))

        call assert_true(spin_momentum(0.5_dp)**2 .isclose. S2_expval_exc([1], NoExc_t()))

        call assert_true(spin_momentum(0.5_dp)**2 .isclose. S2_expval_exc([2], NoExc_t()))

        call assert_equals(2._dp, S2_expval_exc([1, 3], NoExc_t()))

        call assert_true(spin_momentum(1.0_dp)**2 .isclose. S2_expval_exc([2, 4], NoExc_t()))

        call assert_true(1.0_dp .isclose. S2_expval_exc([1, 4], NoExc_t()))

        call assert_true(1.0_dp .isclose. S2_expval_exc([1, 4, 5, 6], DoubleExc_t(1, 2, 4, 3)))
    end subroutine


end module test_SD_spin_purification_mod

program test_SD_spin_purification_prog

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd

    use test_SD_spin_purification_mod

    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_SD_spin_purification_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_SD_spin_purification_prog', 'failed_tests')

        call MPIEnd(.false.)
    end block

contains

    subroutine test_SD_spin_purification_driver()
        call run_test_case(test_spin_momentum, "test_spin_momentum")
        call run_test_case(test_spin_q_num, "test_spin_q_num")
        call run_test_case(test_S2_expval, "test_S2_expval")
        call run_test_case(test_S2_expval_exc, "test_S2_expval_exc")
        call run_test_case(test_get_open_shell, "test_get_open_shell")
    end subroutine
end program test_SD_spin_purification_prog
