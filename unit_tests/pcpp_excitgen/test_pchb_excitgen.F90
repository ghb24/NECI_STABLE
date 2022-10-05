module test_pchb_excitgen_mod
    use fruit, only: assert_true
    use constants, only: dp, maxExcit
    use gasci_pchb, only: possible_PCHB_particle_selection
    use pchb_excitgen, only: PCHB_FCI_excit_generator_t
    use FciMCData, only: pSingles, pDoubles, pParallel
    use SystemData, only: nEl
    use sltcnd_mod, only: dyn_sltcnd_excit_old
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use util_mod, only: near_zero
    use unit_test_helper_excitgen, only: test_excitation_generator, &
                                         init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
                                         FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester
    implicit none
    private
    public :: pchb_test_driver

contains

    subroutine pchb_test_driver()
        integer, parameter :: n_iters = 5 * 10**7
        type(PCHB_FCI_excit_generator_t) :: exc_generator
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10], n_spat_orb = 10
        logical :: successful

        pParallel = 0.05_dp
        pSingles = 0.3_dp
        pDoubles = 1.0_dp - pSingles

        call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))

        call exc_generator%init(possible_PCHB_particle_selection%UNIFORM)

        call run_excit_gen_tester( &
            exc_generator, 'PCHB FCI', &
            opt_nI=det_I, opt_n_dets=n_iters, &
            problem_filter=is_problematic, &
            successful=successful)

        call assert_true(successful)
        call exc_generator%finalize()
        call finalize_excitgen_test()

    contains

        subroutine random_fcidump(iunit)
            integer, intent(in) :: iunit

            call generate_random_integrals( &
                iunit, n_el=size(det_I), n_spat_orb=n_spat_orb, &
                sparse=0.7_dp, sparseT=0.7_dp, total_ms=sum(calc_spin_raw(det_I)))
        end subroutine

        logical function is_problematic(nI, exc, ic, pgen_diagnostic)
            integer, intent(in) :: nI(nEl), exc(2, maxExcit), ic
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = &
                (abs(1._dp - pgen_diagnostic) >= 0.1_dp) &
                .and. .not. near_zero(dyn_sltcnd_excit_old(nI, ic, exc, .true.))
        end function

    end subroutine pchb_test_driver

end module

program test_pchb_excitgen
    use fruit, only: get_failed_count, init_fruit, fruit_summary, fruit_finalize
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_pchb_excitgen_mod, only: pchb_test_driver
    use util_mod, only: stop_all
    implicit none

    call MPIInit(.false.)
    call init_fruit()
    call pchb_test_driver()
    call fruit_summary()
    call fruit_finalize()
    block
        integer :: failed_count
        call get_failed_count(failed_count)
        if (failed_count /= 0) call stop_all('test_pchb_excitgen', 'failed_tests')
    end block
    call MPIEnd(.false.)

end program test_pchb_excitgen
