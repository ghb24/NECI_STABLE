module test_gasci_discarding_mod
    use fruit
    use constants, only: dp, n_int
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use excitation_types, only: Excitation_t

    use gasci, only: GASSpec_t
    use gasci_discarding, only: gen_GASCI_discarding, init_GASCI_discarding, finalize_GASCI_discarding
    use gasci_general, only: gen_all_excits

    use sltcnd_mod, only: dyn_sltcnd_excit
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester
    implicit none
    private
    public :: test_pgen



contains


    subroutine test_pgen()
        use gasci, only: global_GAS_spec => GAS_specification
        use SystemData, only: tGASSpinRecoupling
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GASSpec_t) :: GAS_spec
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10]

        logical :: successful
        integer, parameter :: n_iters=5 * 10**6

        pParallel = 0.05_dp
        pSingles = 0.3_dp
        pDoubles = 1.0_dp - pSingles

        call assert_true(tGASSpinRecoupling)

        GAS_spec = GASSpec_t(n_min=[3, size(det_I)], n_max=[3, size(det_I)], &
                             spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
        call assert_true(GAS_spec%is_valid())
        call assert_true(GAS_spec%contains(det_I))
        global_GAS_spec = GAS_spec

        call init_excitgen_test(size(det_I), FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
        call init_GASCI_discarding()
        call run_excit_gen_tester( &
            gen_GASCI_discarding, 'discarding GASCI implementation, random fcidump', &
            opt_nI=det_I, opt_n_iters=n_iters, &
            gen_all_excits=gen_all_excits, &
            problem_filter=is_problematic,&
            successful=successful)
        call assert_true(successful)
        call finalize_GASCI_discarding()
        call finalize_excitgen_test()

    contains

        subroutine random_fcidump(iunit)
            integer, intent(in) :: iunit
            integer :: n_spat_orb, iGAS

            n_spat_orb = sum([(GAS_spec%GAS_size(iGAS), iGAS = 1, GAS_spec%nGAS())]) .div. 2

            call generate_random_integrals(&
                iunit, n_el=size(det_I), n_spat_orb=n_spat_orb, &
                sparse=1.0_dp, sparseT=1.0_dp, total_ms=sum(calc_spin_raw(det_I)))
        end subroutine

        logical function is_problematic(det_I, exc, pgen_diagnostic)
            type(SpinOrbIdx_t), intent(in) :: det_I
            class(Excitation_t), intent(in) :: exc
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = (abs(1.0_dp - pgen_diagnostic) >= 0.15_dp) &
                              .and. .not. near_zero(dyn_sltcnd_excit(det_I%idx, exc, .true.))
        end function

    end subroutine test_pgen
end module test_gasci_discarding_mod

program test_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_discarding_mod, only: test_pgen


    implicit none
    integer :: failed_count, err

    integer :: n
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
