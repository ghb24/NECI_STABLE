module test_gasci_disconnected_mod
    use fruit, only: assert_true
    use constants, only: dp, n_int, maxExcit
    use SystemData, only: nEl
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use excitation_types, only: Excitation_t

    use excitation_generators, only: ExcitationGenerator_t
    use gasci, only: LocalGASSpec_t
    use gasci_disconnected, only: GAS_disc_ExcGenerator_t

    use sltcnd_mod, only: dyn_sltcnd_excit_old
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        RandomFcidumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester
    implicit none
    private
    public :: test_pgen



contains


    subroutine test_pgen()
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(LocalGASSpec_t) :: GAS_spec
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10]

        logical :: successful
        class(ExcitationGenerator_t), allocatable :: exc_generator
        integer, parameter :: n_iters=10 * 10**5

        pParallel = 0.05_dp
        pSingles = 0.6_dp
        pDoubles = 1.0_dp - pSingles

        GAS_spec = LocalGASSpec_t(n_min=[3, 3], n_max=[3, 3], &
                             spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
        call assert_true(GAS_spec%is_valid())
        call assert_true(GAS_spec%contains_conf(det_I))

        call init_excitgen_test(det_I, &
            RandomFcidumpWriter_t(&
                GAS_spec, det_I, sparse=1.0_dp, sparseT=1.0_dp) &
        )

        exc_generator = GAS_disc_ExcGenerator_t(GAS_spec)
        call run_excit_gen_tester( &
            exc_generator, 'only disconnected implementation, random_fcidump', &
            opt_nI=det_I, opt_n_dets=n_iters, &
            problem_filter=is_problematic,&
            successful=successful)
        call assert_true(successful)
        call finalize_excitgen_test()

    contains

        logical function is_problematic(nI, exc, ic, pgen_diagnostic)
            integer, intent(in) :: nI(nEl), exc(2, maxExcit), ic
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = &
                (abs(1._dp - pgen_diagnostic) >= 0.15_dp) &
                .and. .not. near_zero(dyn_sltcnd_excit_old(nI, ic, exc, .true.))
        end function
    end subroutine test_pgen
end module test_gasci_disconnected_mod

program test_gasci_program
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_disconnected_mod, only: test_pgen


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
        call run_test_case(test_pgen, "test_pgen")
    end subroutine
end program test_gasci_program
