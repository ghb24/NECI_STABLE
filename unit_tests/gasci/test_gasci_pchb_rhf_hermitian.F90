program test_gasci_program

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case, assert_true
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd
    use gasci_pchb_test_helper, only: test_pgen_RHF_hermitian

    implicit none
    integer :: failed_count

    block

        call MPIInit(.false.)

        call init_fruit()

        call test_gasci_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

contains

    subroutine test_gasci_driver()
        call run_test_case(test_pgen_RHF_hermitian, "test_pgen_RHF_hermitian")
        call run_test_case(test_string_conversion, "test_string_conversion")
    end subroutine test_gasci_driver

    subroutine test_string_conversion()
        use gasci_pchb_doubles_main, only: PCHB_DoublesOptions_t, doubles_options_vals

        type(PCHB_DoublesOptions_t) :: options

        options = doubles_options_vals%from_str("UNIF-FULL:FULL-FULL")
        call assert_true(options%particle_selection == doubles_options_vals%particle_selection%UNIF_FULL)
        call assert_true(options%hole_selection == doubles_options_vals%hole_selection%FULL_FULL)

        options = doubles_options_vals%from_str("FULL-FULL:FULL-FULL")
        call assert_true(options%particle_selection == doubles_options_vals%particle_selection%FULL_FULL)
        call assert_true(options%hole_selection == doubles_options_vals%hole_selection%FULL_FULL)

        options = doubles_options_vals%from_str("FULL-FULL:FAST-FAST")
        call assert_true(options%particle_selection == doubles_options_vals%particle_selection%FULL_FULL)
        call assert_true(options%hole_selection == doubles_options_vals%hole_selection%FAST_FAST)

    end subroutine
end program test_gasci_program
