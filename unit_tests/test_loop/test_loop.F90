program test_loop_program

    use Parallel_neci, only: MPIInit, MPIEnd
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use test_loop_testcases, only: test_loop_4ind_wghtd_2, test_loop_pchb, &
                                   test_loop_guga

    implicit none
    integer :: failed_count
    logical :: err

    call MPIInit(err)

    call init_fruit()

    call test_loop_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_loop_program', 'failed_tests')

    call MPIEnd(err)

contains

    subroutine test_loop_driver()
         call run_test_case(test_loop_4ind_wghtd_2, "test_loop")
         call run_test_case(test_loop_pchb, "test_loop_pchb")
         call run_test_case(test_loop_guga, "test_loop_guga")
    end subroutine test_loop_driver
end program test_loop_program
