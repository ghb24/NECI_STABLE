! NOTE
! For meaningful pgen tests `n_iters` has to be higher
! and the `is_problematic` function has to be tighter!

module gasci_pchb_test_helper
    use fruit, only: assert_equals, assert_true, assert_false
    use fortran_strings, only: str
    use constants, only: dp, int64, n_int, maxExcit
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero
    use util_mod, only: factrl, swap, cumsum
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use SystemData, only: nEl
    use excitation_types, only: Excitation_t

    use gasci, only: LocalGASSpec_t
    use gasci_supergroup_index, only: GAS_SuperGroupIndexer_t
    use property_vector_pchb_main, only: GAS_PCHB_ExcGenerator_t, GAS_PCHB_options_t, &
        GAS_PCHB_options_vals
    use property_vector_pchb_doubles_main, only: PCHB_DoublesOptions_t
    use gasci_singles_main, only: GAS_PCHB_SinglesOptions_t, PC_WeightedSinglesOptions_t

    use excitation_generators, only: ExcitationGenerator_t

    use sltcnd_mod, only: dyn_sltcnd_excit_old
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test
    use unit_test_helper_fcidumps, only: RandomFciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester

    use SystemData, only: nEl
    implicit none
    private
    public :: test_pgen_RHF_hermitian, &
              test_pgen_RHF_nonhermitian, &
              test_pgen_UHF_hermitian, &
              test_pgen_UHF_nonhermitian

contains

    subroutine test_pgen_RHF_hermitian()
        call test_pgen_general(.false., .true.)
    end subroutine test_pgen_RHF_hermitian

    subroutine test_pgen_UHF_hermitian()
        call test_pgen_general(.true., .true.)
    end subroutine test_pgen_UHF_hermitian

    subroutine test_pgen_RHF_nonhermitian()
        call test_pgen_general(.false., .false.)
    end subroutine test_pgen_RHF_nonhermitian

    subroutine test_pgen_UHF_nonhermitian()
        call test_pgen_general(.true., .false.)
    end subroutine test_pgen_UHF_nonhermitian

    subroutine test_pgen_general(UHF, hermitian)
        ! TODO make this also accept an options object and parallelise over that
        ! (make this subroutine public)
        use FciMCData, only: pSingles, pDoubles, pParallel
        use SystemData, only: t_non_hermitian_1_body, t_non_hermitian_2_body, &
                              tUHF, tZeroDelimitedFCIDUMP
        use System, only: SetSysDefaults
        use Calc, only: SetCalcDefaults
        logical, intent(in) :: UHF, hermitian
        type(GAS_PCHB_ExcGenerator_t) :: exc_generator
        type(LocalGASSpec_t) :: GAS_spec
        type(GAS_SuperGroupIndexer_t) :: indexer
        type(GAS_PCHB_options_t), allocatable :: settings(:)
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10]

        logical :: successful
        integer :: n_interspace_exc, i
        integer, parameter :: n_iters=10**6
            ! NOTE:
            ! For meaningful pgen tests `n_iters` has to be higher
            ! and the `is_problematic` function has to be tighter.
        call SetCalcDefaults()
        call SetSysDefaults()
        t_non_hermitian_1_body = .not. hermitian
        t_non_hermitian_2_body = .not. hermitian
        tUHF = UHF
        tZeroDelimitedFCIDUMP = UHF

        pParallel = 0.05_dp
        pSingles = 0.2_dp
        pDoubles = 1.0_dp - pSingles

        if (UHF) then ! UHF cannot use spatorbs
            settings = [&
                GAS_PCHB_options_t(&
                        GAS_PCHB_SinglesOptions_t(&
                            GAS_PCHB_options_vals%singles%algorithm%BITMASK_UNIFORM &
                        ), &
                        PCHB_DoublesOptions_t( &
                            GAS_PCHB_options_vals%doubles%particle_selection%UNIF_UNIF, &
                            GAS_PCHB_options_vals%doubles%hole_selection%FAST_FAST, &
                            spin_orb_resolved=.true. &
                        ), &
                        use_lookup=.false. &
                ), &
                GAS_PCHB_options_t(&
                        GAS_PCHB_SinglesOptions_t(&
                            GAS_PCHB_options_vals%singles%algorithm%PC_WEIGHTED, &
                            PC_WeightedSinglesOptions_t(&
                                GAS_PCHB_options_vals%singles%PC_weighted%drawing%UNIF_FULL &
                            ) &
                        ), &
                        PCHB_DoublesOptions_t( &
                            GAS_PCHB_options_vals%doubles%particle_selection%UNIF_FULL, &
                            GAS_PCHB_options_vals%doubles%hole_selection%FULL_FULL, &
                            spin_orb_resolved=.true. &
                        ), &
                        use_lookup=.false. &
                ), &
                GAS_PCHB_options_t(&
                        GAS_PCHB_SinglesOptions_t(&
                            GAS_PCHB_options_vals%singles%algorithm%PC_WEIGHTED, &
                            PC_WeightedSinglesOptions_t(&
                                GAS_PCHB_options_vals%singles%PC_weighted%drawing%FULL_FULL &
                            ) &
                        ), &
                        PCHB_DoublesOptions_t( &
                            GAS_PCHB_options_vals%doubles%particle_selection%FULL_FULL, &
                            GAS_PCHB_options_vals%doubles%hole_selection%FULL_FULL, &
                            spin_orb_resolved=.true. &
                        ), &
                        use_lookup=.false. &
                ) &
            ]
        else
            settings = [&
                GAS_PCHB_options_t(&
                        GAS_PCHB_SinglesOptions_t(&
                            GAS_PCHB_options_vals%singles%algorithm%BITMASK_UNIFORM &
                        ), &
                        PCHB_DoublesOptions_t( &
                            GAS_PCHB_options_vals%doubles%particle_selection%UNIF_UNIF, &
                            GAS_PCHB_options_vals%doubles%hole_selection%FAST_FAST, &
                            spin_orb_resolved=.false. &
                        ), &
                        use_lookup=.false. &
                ), &
                GAS_PCHB_options_t(&
                        GAS_PCHB_SinglesOptions_t(&
                            GAS_PCHB_options_vals%singles%algorithm%PC_WEIGHTED, &
                            PC_WeightedSinglesOptions_t(&
                                GAS_PCHB_options_vals%singles%PC_weighted%drawing%UNIF_FULL &
                            ) &
                        ), &
                        PCHB_DoublesOptions_t( &
                            GAS_PCHB_options_vals%doubles%particle_selection%UNIF_FULL, &
                            GAS_PCHB_options_vals%doubles%hole_selection%FULL_FULL, &
                            spin_orb_resolved=.true. &
                        ), &
                        use_lookup=.false. &
                ), &
                GAS_PCHB_options_t(&
                        GAS_PCHB_SinglesOptions_t(&
                            GAS_PCHB_options_vals%singles%algorithm%PC_WEIGHTED, &
                            PC_WeightedSinglesOptions_t(&
                                GAS_PCHB_options_vals%singles%PC_weighted%drawing%FULL_FULL &
                            ) &
                        ), &
                        PCHB_DoublesOptions_t( &
                            GAS_PCHB_options_vals%doubles%particle_selection%FULL_FULL, &
                            GAS_PCHB_options_vals%doubles%hole_selection%FULL_FULL, &
                            spin_orb_resolved=.true. &
                        ), &
                        use_lookup=.false. &
                ) &
            ]
        end if

        over_settings: do i = 1, size(settings)
            over_interspace_excitations: do n_interspace_exc = 0, 1
                GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                        spat_GAS_orbs=[1, 1, 1, 2, 2, 2], nEl=6)
                call assert_true(GAS_spec%is_valid())
                call assert_true(GAS_spec%contains_conf(det_I))
                indexer = GAS_SuperGroupIndexer_t(GAS_spec)

                call init_excitgen_test(det_I, &
                    RandomFcidumpWriter_t(&
                        GAS_spec, det_I, sparse=1.0_dp, sparseT=1.0_dp, &
                        uhf=uhf, hermitian=hermitian), &
                    setdefaults=.false. &
                )

                call exc_generator%init(&
                    indexer, options=settings(i) &
                )

                call run_excit_gen_tester( &
                    exc_generator, 'general implementation, Li2 like system', &
                    opt_nI=det_I, &
                    opt_n_dets=n_iters, &
                    problem_filter=is_problematic,&
                    successful=successful)
                call exc_generator%finalize()
                call assert_true(successful, "Failed for UHF="//str(UHF)// &
                        ", hermitian="//str(hermitian)//", settings index="//str(i))
                call finalize_excitgen_test()

            end do over_interspace_excitations
        end do over_settings

    contains

        ! NOTE:
        ! For meaningful pgen tests `n_iters` has to be higher
        ! and the `is_problematic` function has to be tighter.
        logical function is_problematic(nI, exc, ic, pgen_diagnostic)
            integer, intent(in) :: nI(nEl), exc(2, maxExcit), ic
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = &
                (abs(1._dp - pgen_diagnostic) >= 1.6_dp) &
                .and. .not. near_zero(dyn_sltcnd_excit_old(nI, ic, exc, .true.))
        end function
    end subroutine test_pgen_general

end module gasci_pchb_test_helper
