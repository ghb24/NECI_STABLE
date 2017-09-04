! replace _template by whatever one needs 
program test_cc_amplitudes

    use fruit
    use cc_amplitudes

    implicit none 

    integer :: failed_count 

    call init_fruit() 
    call cc_amplitudes_test_driver
    call fruit_summary() 
    call fruit_finalize()

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine cc_amplitudes_test_driver

        call run_test_case(dongxia_amplitudes_test, "dongxia_amplitudes_test")

    end subroutine cc_amplitudes_test_driver

    subroutine dongxia_amplitudes_test

        print *, ""
        print *, "test dongxia_amplitudes "

        call assert_true(.false.)

    end subroutine dongxia_amplitudes_test

end program test_cc_amplitudes
