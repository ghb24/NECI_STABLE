module pchb_excitgen_test_helper
    use fruit, only: assert_true, run_test_case
    use constants, only: dp, maxExcit
    use pchb_excitgen, only: PCHB_FCI_excit_generator_t, &
        FCI_PCHB_options_t, FCI_PCHB_options_vals, &
        FCI_PCHB_SinglesOptions_t, PCHB_DoublesOptions_t

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
    public :: test_pgen_rhf_hermitian, test_pgen_rhf_nonhermitian, test_pgen_uhf_hermitian, test_pgen_uhf_nonhermitian

contains

    subroutine test_pgen_rhf_hermitian()
        call pchb_test_general(.false., .true.)
    end subroutine test_pgen_rhf_hermitian

    subroutine test_pgen_uhf_hermitian()
        call pchb_test_general(.true., .true.)
    end subroutine test_pgen_uhf_hermitian

    subroutine test_pgen_rhf_nonhermitian()
        call pchb_test_general(.false., .false.)
    end subroutine test_pgen_rhf_nonhermitian

    subroutine test_pgen_uhf_nonhermitian()
        call pchb_test_general(.true., .false.)
    end subroutine test_pgen_uhf_nonhermitian

    subroutine pchb_test_general(uhf, hermitian)
        use System, only: SetSysDefaults
        use Calc, only: SetCalcDefaults
        use SystemData, only: t_non_hermitian_1_body, t_non_hermitian_2_body, tUHF, tMolpro
        logical, intent(in) :: uhf, hermitian
        integer, parameter :: n_iters = 5 * 10**7
        type(PCHB_FCI_excit_generator_t) :: exc_generator
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10], n_spat_orb = 10
        logical :: successful
        type(FCI_PCHB_options_t) :: options
        character(len=128) :: message

        call SetCalcDefaults()
        call SetSysDefaults()
        t_non_hermitian_1_body = .not. hermitian
        t_non_hermitian_2_body = .not. hermitian
        tUHF = UHF
        tMolpro = UHF

        write(message, *) "Failed for uhf=", uhf, ", hermitian=", hermitian

        pParallel = 0.05_dp
        pSingles = 0.3_dp
        pDoubles = 1.0_dp - pSingles

        call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'), setdefaults=.false.)

        if (uhf) then
            options = FCI_PCHB_options_t(&
                            FCI_PCHB_SinglesOptions_t(&
                                FCI_PCHB_options_vals%singles%algorithm%UNIFORM &
                            ), &
                            PCHB_DoublesOptions_t( &
                                FCI_PCHB_options_vals%doubles%particle_selection%UNIFORM, &
                                FCI_PCHB_options_vals%doubles%hole_selection%SPINORB_FAST_WEIGHTED &
                            ) &
                        )
        else
            options = FCI_PCHB_options_t(&
                            FCI_PCHB_SinglesOptions_t(&
                                FCI_PCHB_options_vals%singles%algorithm%UNIFORM &
                            ), &
                            PCHB_DoublesOptions_t( &
                                FCI_PCHB_options_vals%doubles%particle_selection%UNIFORM, &
                                FCI_PCHB_options_vals%doubles%hole_selection%SPATORB_FAST_WEIGHTED &
                            ) &
                        )
        end if
        call exc_generator%init(options)

        call run_excit_gen_tester( &
            exc_generator, 'PCHB FCI', &
            opt_nI=det_I, opt_n_dets=n_iters, &
            problem_filter=is_problematic, &
            successful=successful)

        call assert_true(successful, trim(message))
        call exc_generator%finalize()
        call finalize_excitgen_test()

    contains

        subroutine random_fcidump(iunit)
            integer, intent(in) :: iunit
            call random_fcidump_general(iunit, uhf, hermitian, det_I, n_spat_orb)
        end subroutine random_fcidump

        logical function is_problematic(nI, exc, ic, pgen_diagnostic)
            integer, intent(in) :: nI(nEl), exc(2, maxExcit), ic
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = &
                (abs(1._dp - pgen_diagnostic) >= 0.15_dp) &
                .and. .not. near_zero(dyn_sltcnd_excit_old(nI, ic, exc, .true.))
        end function is_problematic

    end subroutine pchb_test_general

    subroutine random_fcidump_general(iunit, is_uhf, is_hermitian, det_I, n_spat_orb)
        integer, intent(in) :: iunit
        logical, intent(in) :: is_uhf, is_hermitian
        integer, intent(in) :: det_I(:), n_spat_orb

        call generate_random_integrals( &
            iunit, n_el=size(det_I), n_spat_orb=n_spat_orb, &
            sparse=0.7_dp, sparseT=0.7_dp, total_ms=sum(calc_spin_raw(det_I)), &
            uhf=is_uhf, hermitian=is_hermitian)
    end subroutine random_fcidump_general

end module pchb_excitgen_test_helper
