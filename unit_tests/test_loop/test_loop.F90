program test_loop_program

    use mpi
    use fruit
    use test_loop_testcases, only: test_loop_4ind_wghtd_2, test_loop_pchb, &
                                   test_loop_guga

    implicit none
    integer :: failed_count, err

    call mpi_init(err)

    call init_fruit()

    call test_loop_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_loop_program', 'failed_tests')

    call mpi_finalize(err)

contains

    subroutine test_loop_driver()
         call run_test_case(test_loop_4ind_wghtd_2, "test_loop")
         call run_test_case(test_loop_pchb, "test_loop_pchb")
         call run_test_case(test_loop_guga, "test_loop_guga")
    end subroutine test_loop_driver
end program test_loop_program
