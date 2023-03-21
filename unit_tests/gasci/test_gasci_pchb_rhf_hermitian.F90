program test_gasci_program

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
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
    end subroutine test_gasci_driver
end program test_gasci_program
