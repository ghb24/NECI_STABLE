program test_pchb_excitgen
    use fruit, only: get_failed_count, init_fruit, fruit_summary, fruit_finalize, &
                     run_test_case
    use Parallel_neci, only: MPIInit, MPIEnd
    use pchb_excitgen_test_helper, only: test_pgen_uhf_nonhermitian
    use util_mod, only: stop_all
    implicit none

    call MPIInit(.false.)
    call init_fruit()
    call test_driver()
    call fruit_summary()
    call fruit_finalize()
    block
        integer :: failed_count
        call get_failed_count(failed_count)
        if (failed_count /= 0) call stop_all('test_pchb_excitgen', 'failed_tests')
    end block
    call MPIEnd(.false.)

contains

    subroutine test_driver()
        call run_test_case(test_pgen_uhf_nonhermitian, "test_pgen_uhf_nonhermitian")
    end subroutine test_driver

end program test_pchb_excitgen
