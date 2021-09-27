#include "macros.h"
module FciMCParMod

    ! This module contains the main loop for FCIMC calculations, and the
    ! main per-iteration processing loop.
    use SystemData, only: nel, tUEG2, tGen_4ind_2, &
                          tGen_4ind_weighted, t_test_excit_gen, tGUGA, &
                          t_new_real_space_hubbard, t_tJ_model, t_heisenberg_model, &
                          t_k_space_hubbard, max_ex_level, t_uniform_excits, &
                          tGen_guga_mixed, t_guga_mixed_init, t_guga_mixed_semi, &
                          tReal, t_mixed_excits, &
                          t_crude_exchange_noninits, t_approx_exchange_noninits, &
                          is_init_guga, tGen_sym_guga_ueg, t_analyze_pchb

    use CalcData, only: tFTLM, tSpecLanc, tExactSpec, tDetermProj, tMaxBloom, &
                        tUseRealCoeffs, tWritePopsNorm, tExactDiagAllSym, &
                        AvMCExcits, pops_norm_unit, iExitWalkers, tAdaptiveShift, &
                        iFullSpaceIter, semistoch_shift_iter, tEN2, tOutputInitsRDM, &
                        tOrthogonaliseReplicas, orthogonalise_iter, tNonInitsForRDMs, &
                        tDetermHFSpawning, use_spawn_hash_table, tInitsRDMRef, &
                        ss_space_in, s_global_start, tContTimeFCIMC, tInitsRDM, &
                        trial_shift_iter, tStartTrialLater, &
                        tTrialWavefunction, tSemiStochastic, ntrial_ex_calc, &
                        t_hist_tau_search_option, t_back_spawn, back_spawn_delay, &
                        t_back_spawn_flex, t_back_spawn_flex_option, tSimpleInit, &
                        t_back_spawn_option, tDynamicCoreSpace, coreSpaceUpdateCycle, &
                        DiagSft, tDynamicTrial, trialSpaceUpdateCycle, semistochStartIter, &
                        tSkipRef, tTrialShift, tSpinProject, t_activate_decay, &
                        t_trunc_guga_pgen_noninits, &
                        tLogAverageSpawns, tActivateLAS, eq_cyc, &
                        t_guga_back_spawn, tEN2Init, tEN2Rigorous, tDeathBeforeComms, &
                        tDetermProjApproxHamil, tCoreAdaptiveShift, &
                        tScaleBlooms, max_allowed_spawn

    use adi_data, only: tReadRefs, tDelayGetRefs, allDoubsInitsDelay, &
                        tDelayAllDoubsInits, tReferenceChanged, &
                        SIUpdateInterval, tSuppressSIOutput, nRefUpdateInterval, &
                        SIUpdateOffset

    use LoggingData, only: tJustBlocking, tCompareTrialAmps, tChangeVarsRDM, &
                           tWriteCoreEnd, tNoNewRDMContrib, tPrintPopsDefault, &
                           compare_amps_period, PopsFileTimer, &
                           write_end_core_size, t_calc_double_occ, t_calc_double_occ_av, &
                           equi_iter_double_occ, t_print_frq_histograms, ref_filename, &
                           t_spin_measurements, &
                           tCoupleCycleOutput, StepsPrint, &
                           tFValEnergyHist, tFValPopHist, &
                           FvalEnergyHist_EnergyBins, FvalEnergyHist_FValBins, &
                           FvalPopHist_PopBins, FvalPopHist_FValBins, &
                           tHDF5TruncPopsWrite, iHDF5TruncPopsEx, tAccumPops, &
                           tAccumPopsActive, iAccumPopsIter, iAccumPopsExpireIters, &
                           tPopsProjE, iHDF5TruncPopsIter, iAccumPopsCounter, &
                           AccumPopsExpirePercent, t_print_core_vec

    use rdm_data, only: print_2rdm_est, ThisRDMIter, inits_one_rdms, two_rdm_inits_spawn, &
                        two_rdm_inits, rdm_inits_defs, RDMCorrectionFactor, inits_estimates, tSetupInitsEst, &
                        tApplyLC

    use rdm_finalising, only: finalise_rdms
    use rdm_general, only: init_rdms, SumCorrectionContrib, UpdateRDMCorrectionTerm
    use rdm_filling, only: fill_rdm_offdiag_deterministic, fill_rdm_diag_wrapper
    use rdm_explicit, only: fill_explicitrdm_this_iter, fill_hist_explicitrdm_this_iter
    use procedure_pointers, only: attempt_die_t, generate_excitation_t, &
                                  get_spawn_helement_t
    use semi_stoch_gen, only: write_most_pop_core_at_end, init_semi_stochastic, &
                              refresh_semistochastic_space
    use semi_stoch_procs, only: is_core_state, check_determ_flag, &
                                determ_projection, average_determ_vector, &
                                determ_projection_no_death, core_space_pos, &
                                print_determ_vec_av, print_determ_vec
    use trial_wf_gen, only: update_compare_trial_file, init_trial_wf, refresh_trial_wf
    use hist, only: write_zero_hist_excit_tofrom
    use orthogonalise, only: orthogonalise_replicas, calc_replica_overlaps, &
                             orthogonalise_replica_pairs
    use bit_reps, only: set_flag, clr_flag, any_run_is_initiator, &
                        all_runs_are_initiator
    use load_balance, only: tLoadBalanceBlocks, adjust_load_balance, RemoveHashDet, &
                            need_load_balancing, loadBalanceInterval
    use exact_diag, only: perform_exact_diag_all_symmetry
    use spectral_lanczos, only: perform_spectral_lanczos
    use bit_rep_data, only: IlutBits, flag_determ_parent, test_flag, flag_prone
    use errors, only: standalone_errors, error_analysis
    use PopsFileMod, only: WriteToPopsFileParOneArr
    use AnnihilationMod, only: DirectAnnihilation, communicate_and_merge_spawns, &
                               rm_non_inits_from_spawnedparts
    use exact_spectrum, only: get_exact_spectrum
    use determ_proj, only: perform_determ_proj, perform_determ_proj_approx_ham
    use cont_time, only: iterate_cont_time
    use global_det_data, only: det_diagH, reset_tau_int, get_all_spawn_pops, &
                               reset_shift_int, update_shift_int, &
                               update_tau_int, set_spawn_pop, &
                               get_tot_spawns, get_acc_spawns, &
                               update_pops_sum_all, get_pops_iter, &
                               replica_est_len, get_max_ratio, update_max_ratio
    use DetBitOps, only: tAccumEmptyDet

    use RotateOrbsMod, only: RotateOrbs
    use NatOrbsMod, only: PrintOrbOccs
    use ftlm_neci, only: perform_ftlm
    use hash, only: clear_hash_table
    use soft_exit, only: ChangeVars
    use adi_references, only: setup_reference_space, enable_adi, adjust_nRefs, &
                              update_reference_space
    use fcimc_initialisation
    use fcimc_iter_utils
    use replica_estimates
    use neci_signals
    use fcimc_helper, only: create_particle, create_particle_with_hash_table
    use fcimc_output
    use FciMCData
    use constants

    use excit_gen_5, only: gen_excit_4ind_weighted2

    use guga_excitations, only: deallocate_projE_list, generate_excitation_guga, &
                                global_excitInfo
    use guga_bitrepops, only: fill_csf_i, current_csf_i
    use tJ_model, only: init_guga_heisenberg_model, init_guga_tj_model

    use real_time_data, only: t_prepare_real_time, n_real_time_copies, &
                              cnt_real_time_copies

    use real_time_init, only: init_overlap_buffers

    use bit_reps, only: decode_bit_det

    use util_mod, only: operator(.div.), toggle_lprof

    use hdiag_from_excit, only: get_hdiag_from_excit, get_hdiag_bare_hphf

    use double_occ_mod, only: get_double_occupancy, inst_double_occ, &
                              rezero_double_occ_stats, write_double_occ_stats, &
                              sum_double_occ, sum_norm_psi_squared, finalize_double_occ_and_spin_diff, &
                              measure_double_occ_and_spin_diff, rezero_spin_diff, &
                              write_spin_diff_stats, write_spat_doub_occ_stats, &
                              all_sum_double_occ, calc_double_occ_from_rdm
    use tau_search_hist, only: print_frequency_histograms, deallocate_histograms
    use back_spawn, only: init_back_spawn
    use real_space_hubbard, only: init_real_space_hubbard, gen_excit_rs_hubbard
    use tJ_model, only: init_tJ_model, init_heisenberg_model

    use k_space_hubbard, only: init_k_space_hubbard, gen_excit_k_space_hub_transcorr, &
                               gen_excit_uniform_k_space_hub_transcorr, &
                               gen_excit_mixed_k_space_hub_transcorr, &
                               gen_excit_k_space_hub
    use cc_amplitudes, only: t_cc_amplitudes, init_cc_amplitudes, cc_delay, &
                             t_plot_cc_amplitudes, print_cc_amplitudes

    use analyse_wf_symmetry, only: analyze_wavefunction_symmetry, t_symmetry_analysis

    use sltcnd_mod, only: sltcnd_excit
    use hdf5_popsfile, only: write_popsfile_hdf5
    use local_spin, only: measure_local_spin, write_local_spin_stats, &
                          finalize_local_spin_measurement

    implicit none

    !array for timings of the main compute loop
    real(dp), dimension(100) :: lt_arr

    integer, parameter :: lb_measure_cycle = 100

contains

    subroutine FciMCPar(energy_final_output)

        use rdm_data, only: rdm_estimates, two_rdm_main, two_rdm_recv, two_rdm_recv_2
        use rdm_data, only: two_rdm_spawn, one_rdms, rdm_definitions, en_pert_main
        use rdm_estimators, only: calc_2rdm_estimates_wrapper, write_rdm_estimates
        use rdm_data_utils, only: clear_rdm_list_t, add_rdm_1_to_rdm_2
        USE MolproPlugin, only: MolproPluginResult

        real(dp), intent(out), allocatable :: energy_final_output(:)

        integer :: iroot, isymh
        real(dp) :: Weight, Energyxw, BestEnergy
        INTEGER :: error, irdm
        LOGICAL :: TIncrement, tWritePopsFound, tSingBiasChange, tPrintWarn
        REAL(dp) :: s_start, s_end, tstart(2), tend(2), totaltime
        real(dp) :: TotalTime8
        real(dp) :: Inpair(2), Outpair(2)
        integer, dimension(lenof_sign) :: tmp_sgn
        integer :: tmp_int(lenof_sign), i, istart, iRDMSamplingIter
        real(dp) :: grow_rate, EnergyDiff, Norm_2RDM
        TYPE(BasisFn) RefSym
        real(dp) :: mean_ProjE_re, mean_ProjE_im, mean_Shift
        real(dp) :: ProjE_Err_re, ProjE_Err_im, Shift_Err
        logical :: tNoProjEValue, tNoShiftValue
        real(dp) :: BestErr
        real(dp) :: start_time, stop_time
        logical :: tStartedFromCoreGround
        real(dp), dimension(100) :: lt_sum, lt_max

        real(dp):: lt_imb, lt_imb_cycle
        integer:: rest, err, allErr

        real(dp) :: CurrentSign(lenof_sign)
        integer :: j
        integer :: pops_iter, nJ(nel)

        HElement_t(dp):: InstE(inum_runs)
        HElement_t(dp):: AccumE(inum_runs)
        integer :: ExcitLevel

        logical :: t_comm_done, tScheduledLoadBalance
        integer :: run

        ! Procedure pointer temporaries
        procedure(generate_excitation_t), pointer :: ge_tmp
        procedure(get_spawn_helement_t), pointer :: gs_tmp
        procedure(attempt_die_t), pointer :: ad_tmp

        character(*), parameter :: this_routine = 'FciMCPar'
        character(6), parameter :: excit_descriptor(0:3) = &
                                   (/"IC0   ", "single", "double", "triple"/)

        integer :: tmp_det(nel)

        if (tJustBlocking) then
            ! Just reblock the current data, and do not perform an fcimc calculation.
            write(stdout, "(A)") "Skipping FCIQMC calculation and simply reblocking previous output"
            call Standalone_Errors()
            return
        end if
        err = 0
        ProjE_Err_re = 1.0_dp
        ProjE_Err_im = 1.0_dp
        shift_err = 1.0_dp

        TDebug = .false.  ! Set debugging flag

        ! This is set here not in SetupParameters, as otherwise it would be
        ! wiped just when we need it!
        tPopsAlreadyRead = .false.

        call SetupParameters()
        call init_fcimc_fn_pointers()
        call InitFCIMCCalcPar()

        if (tLogGreensfunction .and. .not. t_real_time_fciqmc) then
            call init_overlap_buffers()
        end if

        if (t_new_real_space_hubbard) then
            call init_real_space_hubbard()
        end if
        if (t_tJ_model) then
            if (tGUGA) then
                call init_guga_tj_model()
            else
                call init_tJ_model()
            end if
        end if
        if (t_heisenberg_model) then
            if (tGUGA) then
                call init_guga_heisenberg_model()
            else
                call init_heisenberg_model()
            end if
        end if
        ! try to call this earlier..
        ! just do it twice for now..
        if (t_k_space_hubbard) then
            call init_k_space_hubbard()
        end if

#ifdef DEBUG_
        call decode_bit_det(tmp_det, ilutHF)
        write(stdout, *) "HF: ", tmp_det
        call decode_bit_det(tmp_det, ilutHF_true)
        write(stdout, *) "HF_true: ", tmp_det
        call decode_bit_det(tmp_det, ilutRef(:, 1))
        write(stdout, *) "Ref: ", tmp_det
        write(stdout, *) "ProjEDet: ", ProjEDet
#endif

        ! Attach signal handlers to give a more graceful death-mannerism
        call init_signals()

        ! We want to do some population checking before we run any iterations.
        ! In the normal case this is run between iterations, but it is
        ! helpful to do it here.
        call population_check()

        if (n_int /= int64) then
            call stop_all('setup parameters', 'Use of realcoefficients requires 64 bit integers.')
        end if

        if (tDetermProj) then
            ! If performing a deterministic projection instead of an FCIQMC calc:
            if (tDetermProjApproxHamil) then
                call perform_determ_proj_approx_ham()
            else
                call perform_determ_proj()
            end if
            return
        else if (tFTLM) then
            ! If performing a finite-temperature Lanczos method job instead of FCIQMC:
            call perform_ftlm()
            return
        else if (tSpecLanc) then
            call perform_spectral_lanczos()
            return
        else if (tExactSpec) then
            call get_exact_spectrum()
            return
        else if (tExactDiagAllSym) then
            call perform_exact_diag_all_symmetry()
            return
        end if

        ! Initial output
        if (tFCIMCStats2) then
            call write_fcimcstats2(iter_data_fciqmc, initial=.true.)
            call write_fcimcstats2(iter_data_fciqmc)
        else
            call WriteFciMCStatsHeader()
            ! Prepend a # to the initial status line so analysis doesn't pick
            ! up repetitions in the FCIMCStats or INITIATORStats files from
            ! restarts.
            !write(stdout,'("#")', advance='no')
            if (iProcIndex == root) then
                write(fcimcstats_unit, '("#")', advance='no')
                if (inum_runs == 2) &
                    write(fcimcstats_unit2, '("#")', advance='no')
                write(initiatorstats_unit, '("#")', advance='no')
                if (tLogEXLEVELStats) &
                    write(EXLEVELStats_unit, '("#")', advance='no')
            end if
            call WriteFCIMCStats()
        end if

        if (t_measure_local_spin) then
            call write_local_spin_stats(initial = .true.)
            call write_local_spin_stats()
        end if
        ! double occupancy:
        if (t_calc_double_occ) then
            call write_double_occ_stats(initial=.true.)
            call write_double_occ_stats()

            if (t_spin_measurements) then

                call write_spin_diff_stats(initial=.true.)
                call write_spin_diff_stats()

                call write_spat_doub_occ_stats(initial=.true.)
                call write_spat_doub_occ_stats()
            end if
        end if
        ! Put a barrier here so all processes synchronise before we begin.
        call MPIBarrier(error)

        ! Start MC simulation...
        !
        ! If tIncrument is true, it means that when it comes out of the loop,
        ! it wants to subtract 1 from the iteration count to get the true
        ! number of iterations
        tIncrement = .true.
        Iter = 1
        iRDMSamplingIter = 1    !For how many iterations have we accumulated the RDM

        SumSigns = 0.0_dp
        SumSpawns = 0.0_dp

        ! In we go - start the timer for scaling curve!
        start_time = neci_etime(tstart)
        lt_imb = 0.
        lt_max = 0.
        lt_sum = 0.
        ! For calculations with only few iterations
        lt_arr = 0.
        lt_imb_cycle = 0.


        call toggle_lprof()
        main_iteration_loop: do while (.true.)
            if (TestMCExit(Iter, iRDMSamplingIter)) then
                ! The popsfile requires the right total walker number, so
                ! update it (TotParts is updated in the annihilation step)
                call MPISumAll(TotParts, AllTotParts)
                exit
            end if

            ! start logging average spawns once we enter the variable shift mode
            if (.not. any(tSinglePartPhase) .and. tActivateLAS) tLogAverageSpawns = .true.

            IFDEBUG(FCIMCDebug, 2) write(stdout, *) 'Iter', iter

            if (iProcIndex == root) s_start = neci_etime(tstart)

            ! Update the semistochastic space if requested
            if (tSemiStochastic .and. tDynamicCoreSpace .and. &
                mod(iter - semistochStartIter, &
                    coreSpaceUpdateCycle) == 0) then
                call refresh_semistochastic_space()
                write(stdout, *) "Refereshing semistochastic space at iteration ", iter
            end if

            ! Is this an iteration where semi-stochastic is turned on?
            if (semistoch_shift_iter /= 0 .and. .not. any(tSinglePartPhase)) then
                if ((Iter - maxval(VaryShiftIter)) == semistoch_shift_iter + 1) then
                    tSemiStochastic = .true.
                    call init_semi_stochastic(ss_space_in, tStartedFromCoreGround)
                    if (tStartedFromCoreGround) call set_initial_run_references()
                    ! Count iterations for corespace updates from here
                    semistochStartIter = iter
                    ! and switch how iterations for SI updates are counted
                    SIUpdateOffset = semistochStartIter
                end if
            end if

            ! Update the trial wavefunction if requested
            if (tTrialWavefunction .and. tDynamicTrial .and. &
                mod(iter - trial_shift_iter, trialSpaceUpdateCycle) == 0) then
                if (tPairedReplicas) then
                    call refresh_trial_wf(trial_space_in, ntrial_ex_calc, &
                                          inum_runs.div.2, .true.)
                else
                    call refresh_trial_wf(trial_space_in, ntrial_ex_calc, &
                                          inum_runs, .false.)
                end if
                write(stdout, *) "Refreshing trial wavefunction at iteration ", iter
            end if

            if (((Iter - maxval(VaryShiftIter)) == allDoubsInitsDelay + 1 &
                 .and. all(.not. tSinglePartPhase))) then
                ! Start the all-doubs-initiator procedure
                if (tDelayAllDoubsInits) call enable_adi()
                ! If desired, we now set up the references for the purpose of the
                ! all-doubs-initiators
                if (tDelayGetRefs) then
                    ! Re-initialize the reference space
                    call setup_reference_space(.true.)
                end if
            end if

            ! turn on double occ measurement after equilibration
            if (equi_iter_double_occ /= 0 .and. all(.not. tSinglePartPhase)) then
                if ((iter - maxval(VaryShiftIter)) == equi_iter_double_occ + 1) then
                    t_calc_double_occ_av = .true.
                end if
            end if

            ! [W.D]
            ! option to enable back-spawning after a certain number of iterations
            ! but this should be independent if the shift is varied already (or?)
            if (back_spawn_delay /= 0 .and. iter == back_spawn_delay + 1) then
                if (t_back_spawn_flex_option) then
                    t_back_spawn_flex = .true.
                else if (t_back_spawn_option) then
                    t_back_spawn = .true.
                end if
                call init_back_spawn()
            end if

            if (t_cc_amplitudes .and. cc_delay /= 0 .and. all(.not. tSinglePartPhase)) then
                if ((iter - maxval(VaryShiftIter)) == cc_delay + 1) then
                    ! for now just test if it works
                    call init_cc_amplitudes()
                end if
            end if

            ! Is this an iteration where trial-wavefunction estimators are
            ! turned on?
            if (tStartTrialLater .and. all(.not. tSinglePartPhase)) then
                if ((Iter - maxval(VaryShiftIter)) == trial_shift_iter + 1) then
                    tTrialWavefunction = .true.

                    if (tPairedReplicas) then
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs.div.2, .true.)
                    else
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs, .false.)
                    end if
                end if
            end if

            if (equi_iter_double_occ /= 0 .and. all(.not. tSinglePartPhase)) then
                if ((iter - maxval(VaryShiftIter)) == equi_iter_double_occ + 1) then
                    t_calc_double_occ_av = .true.
                end if
            end if

            if (tRDMonFly .and. (.not. tFillingExplicRDMonFly) &
                & .and. (.not. tFillingStochRDMonFly)) call check_start_rdm()

            if (tContTimeFCIMC) then
                call iterate_cont_time(iter_data_fciqmc)
            else
                call PerformFciMCycPar(iter_data_fciqmc, err)
                if (err /= 0) call stop_all(this_routine, 'Failure from PerformFciMCycPar')
            end if

            if (tAccumPops .and. iter + PreviousCycles >= iAccumPopsIter) then
                if (.not. tAccumPopsActive) then
                    tAccumPopsActive = .true.
                    write(stdout, *) "Starting to accumulate populations ..."
                end if

                call update_pops_sum_all(TotWalkers, iter + PreviousCycles)
                iAccumPopsCounter = iAccumPopsCounter + 1

                ! The currentdets is almost full, we should start removing
                ! dets which have been empty long enough
                if (iAccumPopsExpireIters > 0 .and. TotWalkers > AccumPopsExpirePercent * real(MaxWalkersPart, dp)) then
                    do j = 1, int(TotWalkers, sizeof_int)
                        ! The loop is over empty dets only
                        call extract_sign(CurrentDets(:, j), CurrentSign)
                        if (.not. IsUnoccDet(CurrentSign)) cycle

                        ! Keep semi-stochastic dets
                        if (check_determ_flag(CurrentDets(:, j))) cycle

                        ! Keep up to double excitations
                        ExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:, j))
                        if (ExcitLevel <= 2) cycle

                        ! Check if the det has already been removed
                        if (test_flag(CurrentDets(:, j), flag_removed)) cycle

                        ! Now if it has been empty for a long time, remove it
                        pops_iter = INT(get_pops_iter(j))
                        if (iter + PreviousCycles - pops_iter > iAccumPopsExpireIters) then
                            call decode_bit_det(nJ, CurrentDets(:, j))
                            call RemoveHashDet(HashIndex, nJ, j)
                        end if
                    end do
                end if
            end if

            ! if the iteration failed, stop the calculation now
            if (err /= 0) then
                call stop_all(this_routine, 'if the iteration failed, stop the calculation now')
            end if

            if (iProcIndex == root) then
                s_end = neci_etime(tend)
                IterTime = real(IterTime + (s_end - s_start), kind=sp)
            end if

            if (loadBalanceInterval > 0) then
                tScheduledLoadBalance = mod(iter, loadBalanceInterval) == 0
            else
                tScheduledLoadBalance = .false.
            end if

            ! Add some load balancing magic!
            if (tLoadBalanceBlocks .and. (tScheduledLoadBalance .or. &
                                          (mod(iter, lb_measure_cycle) == 1 .and. loadBalanceInterval == 0)) .and. &
                .not. tSemiStochastic .and. .not. tFillingStochRDMOnFly) then
                ! Use the ratio of time lost due to load imbalance as an estimtor
                ! whether load balancing should be used
                if (iter > lb_measure_cycle) then
                    lt_imb_cycle = lt_imb_cycle / sum(lt_sum)
                else
                    lt_imb_cycle = 0.0
                end if
                if (need_load_balancing(lt_imb_cycle) .or. tScheduledLoadBalance) then
                    call adjust_load_balance(iter_data_fciqmc)
                end if
            end if

            if (SIUpdateInterval > 0) then
                ! Regular update of the superinitiators. Use with care as it
                ! is still rather expensive if secondary superinitiators are used
                if (mod(iter - SIUpdateOffset, SIUpdateInterval) == 0) then
                    ! the reference population needs some time to equilibrate
                    ! hence, nRefs cannot be updated that often
                    if (mod(iter, nRefUpdateInterval) == 0) call adjust_nRefs()
                    call update_reference_space(tReadPops .or. all(.not. tSinglePartPhase))
                end if
            end if

            ! as alternative to the update cycle bound output, a specified output interval
            ! may be given that is not correlated with the update cycle
            ! -> do the output independently, only requires communication + write
            ! An output frequency below 1 means no output
            t_comm_done = .false.
            if (StepsPrint > 0 .and. .not. tCoupleCycleOutput) then
                if (mod(Iter, StepsPrint) == 0) then
                    ! just perform a communication + output, without update
                    call iteration_output_wrapper(iter_data_fciqmc, TotParts, tPairedReplicas)
                    ! mark that the communication has been done
                    t_comm_done = .true.
                end if
            end if

            if (mod(Iter, StepsSft) == 0) then

                ! Has there been a particle bloom this update cycle? Loop
                ! through the spawned particle types, and output details.
                ! If outputting only the biggest blooms to date, keep
                ! track of that.
                if (iProcIndex == Root) then
                    istart = 1
                    ! if (tSpinProjDets) istart = 0
                    do i = istart, max_ex_level
                        if (bloom_count(i) /= 0) then
                            if (.not. tMaxBloom .or. &
                                bloom_sizes(i) > bloom_max(i)) then
                                bloom_max(i) = bloom_sizes(i)

                                write(stderr, bloom_warn_string) &
                                    trim(excit_descriptor(i)), bloom_sizes(i), &
                                    bloom_count(i)

                                if (tUEG2) &
                                    call stop_all(this_routine, "Should never &
                                                 &bloom in UEG2")

                            end if
                        end if
                    end do
                end if
                ! Zero the accumulators
                bloom_sizes = 0
                bloom_count = 0
                ! Calculate the a new value for the shift (amongst other
                ! things). Generally, collate information from all processors,
                ! update statistics and output them to the user.
                call set_timer(Stats_Comms_Time)
                call calculate_new_shift_wrapper(iter_data_fciqmc, TotParts, &
                                                 tPairedReplicas, t_comm_req=.not. t_comm_done)
                ! If the output is in sync, do a non-communicated output
                if (tCoupleCycleOutput) call iteration_output_wrapper(iter_data_fciqmc, TotParts, &
                                                                      tPairedReplicas, t_comm_req=.false.)
                call halt_timer(Stats_Comms_Time)

                if (t_measure_local_spin) then
                    call write_local_spin_stats()
                end if
                ! in calculate_new_shift_wrapper output is plotted too!
                ! so for now do it here for double occupancy
                if (t_calc_double_occ) then
                    call write_double_occ_stats()
                    if (t_spin_measurements) then
                        call write_spin_diff_stats()
                        call write_spat_doub_occ_stats()
                    end if
                end if

                if (tRestart) cycle

                IF ((tTruncCAS .or. tTruncSpace .or. tTruncInitiator) .and. (Iter > iFullSpaceIter) &
                    .and. (iFullSpaceIter /= 0)) THEN
!Test if we want to expand to the full space if an EXPANDSPACE variable has been set
                    IF (tHistSpawn .or. tCalcFCIMCPsi) THEN
                        IF (iProcIndex == 0) write(stdout, *) "Unable to expand space since histgramming the wavefunction..."
                    ELSE
                        ICILevel = 0
                        tTruncSpace = .false.
                        tTruncCAS = .false.
                        IF (tTruncInitiator) tTruncInitiator = .false.
                        IF (iProcIndex == 0) THEN
                            write(stdout, *) "Expanding to the full space on iteration ", Iter + PreviousCycles
                        end if
                    end if
                end if

                if (iProcIndex == root) TotalTime8 = real(s_end - s_global_start, dp)

                call MPIBCast(TotalTime8)    !TotalTime is local - broadcast to all procs

!This routine will check for a CHANGEVARS file and change the parameters of the calculation accordingly.
                CALL ChangeVars(tSingBiasChange, tWritePopsFound)
                IF (tSoftExitFound) THEN
                    !Now changed such that we do one more block of iterations and then exit, to allow for proper
                    !inclusion of all the RDM elements
                    NMCyc = Iter + StepsSft
                    tSoftExitFound = .false.

                end if
                IF (tTimeExit .and. (TotalTime8 >= MaxTimeExit)) THEN
                    !Is it time to exit yet?
                    write(stdout, "(A,F8.2,A)") "Time limit reached for simulation of: ", MaxTimeExit / 60.0_dp, " minutes - exiting..."
                    NMCyc = Iter + StepsSft
                    ! Set this to false so that this if statement won't be entered next time.
                    tTimeExit = .false.
                end if
                IF (iExitWalkers /= -1_int64 .and. sum(AllTotParts) > iExitWalkers) THEN
                    !Exit criterion based on total walker number met.
                    write(stdout, "(A,I15)") "Total walker population exceeds that given by &
                        &EXITWALKERS criteria - exiting...", sum(AllTotParts)
                    tIncrement = .false.
                    exit
                end if
                IF (tWritePopsFound) THEN
!We have explicitly asked to write out the POPSFILE from the CHANGEVARS file.
                    CALL WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
                end if
                IF (tSingBiasChange) THEN
                    CALL CalcApproxpDoubles()
                end if

                if ((PopsfileTimer > 0.0_dp) .and. ((iPopsTimers * PopsfileTimer) < (TotalTime8 / 3600.0_dp))) then
                    !Write out a POPSFILE every PopsfileTimer hours
                    if (iProcIndex == Root) then
                        CALL RENAME('popsfile.h5', 'popsfile.h5.bk')
                        CALL RENAME('POPSFILEBIN', 'POPSFILEBIN.bk')
                        CALL RENAME('POPSFILEHEAD', 'POPSFILEHEAD.bk')
                        write(stdout, "(A,F7.3,A)") "Writing out a popsfile after ", iPopsTimers * PopsfileTimer, " hours..."
                    end if
                    call WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
                    iPopsTimers = iPopsTimers + 1
                    if (iProcIndex == Root) then
                        s_end = neci_etime(tend)
                        write(stdout, "(A,F7.3,A)") "Time taken to write out POPSFILE: ", real(s_end, dp) - TotalTime8, " seconds."
                    end if
                end if
                if (tHistExcitToFrom) &
                    call write_zero_hist_excit_tofrom()

            end if   !Endif end of update cycle

            IF (TPopsFile .and. (.not. tPrintPopsDefault) .and. (mod(Iter, iWritePopsEvery) == 0)) THEN
                ! differentiate between normal routine and the real-time
                ! preperation
                if (t_prepare_real_time) then
                    if (cnt_real_time_copies < n_real_time_copies) then
                        CALL WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
                        cnt_real_time_copies = cnt_real_time_copies + 1
                    else
                        ! if the number of wanted copies is reached exit the
                        ! calculation. REMINDER: at the end of the calc. an
                        ! additional popsfile is printed -> so start count at 1
                        ! todo: cleanly exit calc. here!
                        NMCyc = Iter + StepsSft
                        t_prepare_real_time = .false.
                        ! do want to finish all the rdm stuff to go on, but
                        ! do not want anymore popsfile to be printed except
                        ! at the very end..
                        tPrintPopsDefault = .true.
                    end if
                else
                    !This will write out the POPSFILE if wanted
                    CALL WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
                end if
            end if

            if (TPopsFile .and. tHDF5TruncPopsWrite .and. iHDF5TruncPopsIter > 0) then
                if (mod(Iter, iHDF5TruncPopsIter) == 0) then
                    call write_popsfile_hdf5(iHDF5TruncPopsEx, .true.)
                end if
            end if
            IF (tHistSpawn .and. (mod(Iter, iWriteHistEvery) == 0) .and. (.not. tRDMonFly)) THEN
                CALL WriteHistogram()
            end if

            ! accumulate the rdm correction due to adaptive shift
            if (tAdaptiveShift .and. all(.not. tSinglePartPhase)) call UpdateRDMCorrectionTerm()

            if (tRDMonFly .and. all(.not. tSinglePartPhase)) then
                ! If we wish to calculate the energy, have started accumulating the RDMs,
                ! and this is an iteration where the energy should be calculated, do so.
                if (print_2rdm_est .and. ((Iter - maxval(VaryShiftIter)) > IterRDMonFly) &
                    .and. (mod(Iter + PreviousCycles - IterRDMStart + 1, RDMEnergyIter) == 0)) then

                    ! rezero the count of how many iterations we have been averaging over
                    ThisRDMIter = 0.0_dp
                    if (tOutputInitsRDM) then
                        call calc_2rdm_estimates_wrapper(rdm_inits_defs, inits_estimates, &
                                                         two_rdm_inits, en_pert_main)
                    end if
                    if (tInitsRDMRef .and. tNonInitsForRDMs) then
                        ! add the initiator-only rdm of this cycle to the main rdm (rescaled
                        ! with the correction factor)
                        call add_rdm_1_to_rdm_2(two_rdm_inits, two_rdm_main, RDMCorrectionFactor)
                        ! get the inits-only energy
                        call calc_2rdm_estimates_wrapper(rdm_inits_defs, inits_estimates, &
                                                         two_rdm_inits, en_pert_main)
                        tSetupInitsEst = .true.
                        ! and reset the initiator-only rdm
                        call clear_rdm_list_t(two_rdm_inits)
                    end if
                    call calc_2rdm_estimates_wrapper(rdm_definitions, rdm_estimates, two_rdm_main, en_pert_main)

                    if (iProcIndex == 0) then
                        if (.not. tInitsRDMRef .or. tSetupInitsEst) &
                            call write_rdm_estimates(rdm_definitions, rdm_estimates, .false., print_2rdm_est, .false.)
                        if (tOutputInitsRDM) call write_rdm_estimates(rdm_inits_defs, &
                                                                      inits_estimates, .false., print_2rdm_est, .true.)
                    end if

                    if (tEN2) then
                        ! If calculating the Epstein-Nesbet perturbation, reset the
                        ! array and hash table where contributions are accumulated.
                        en_pert_main%ndets = 0
                        call clear_hash_table(en_pert_main%hash_table)
                    end if
                end if

            end if

            if (tChangeVarsRDM) then
                ! Decided during the CHANGEVARS that the RDMs should be calculated.
                call init_rdms(rdm_definitions%nrdms_standard, rdm_definitions%nrdms_transition)
                tRDMonFly = .true.
                tChangeVarsRDM = .false.
            end if

            if (tDiagWalkerSubspace .and. (mod(Iter, iDiagSubspaceIter) == 0)) then
                ! Diagonalise a subspace consisting of the occupied determinants
                call DiagWalkerSubspace()
            end if

            ! If requested and on a correct iteration, update the COMPARETRIAL file.
            if (tCompareTrialAmps) then
                ASSERT(compare_amps_period /= 0)
                if (mod(Iter, compare_amps_period) == 0) then
                    call update_compare_trial_file(.false.)
                end if
            end if

            ! Compute the time lost due to load imbalance - aggregation done for 100 iterations
            ! at a time to avoid unnecessary synchronisation points
            if (mod(Iter, lb_measure_cycle) == 0) then
                call MPIAllReduce(lt_arr, MPI_SUM, lt_sum)
                call MPIAllReduce(lt_arr, MPI_MAX, lt_max)
                lt_imb_cycle = sum(lt_max - lt_sum / nProcessors)
                lt_imb = lt_imb + lt_imb_cycle
            end if

            Iter = Iter + 1
            if (tFillingStochRDMonFly) iRDMSamplingIter = iRDMSamplingIter + 1

        end do main_iteration_loop
        call toggle_lprof()

        ! Final output is always enabled
        tSuppressSIOutput = .false.

        ! We are at the end - get the stop-time. Output the timing details
        stop_time = neci_etime(tend)
        write(stdout, *) '- - - - - - - - - - - - - - - - - - - - - - - -'
        write(stdout, *) 'Total loop-time: ', stop_time - start_time

        !add load imbalance from remaining iterations (if any)
        rest = mod(Iter - 1, lb_measure_cycle)
        if (rest > 0) then
            call MPIReduce(lt_arr, MPI_SUM, lt_sum)
            call MPIReduce(lt_arr, MPI_MAX, lt_max)
            lt_imb = lt_imb + sum(lt_max(1:rest) - lt_sum(1:rest) / nProcessors)
        end if
        if (iProcIndex == 0) write(stdout, *) 'Time lost due to load imbalance: ', lt_imb
        write(stdout, *) '- - - - - - - - - - - - - - - - - - - - - - - -'

        ! [Werner Dobrautz 4.4.2017]
        ! for now always print out the frequency histograms for the
        ! tau-search.. maybe change that later to be an option
        ! to be turned off
        if (t_print_frq_histograms .and. t_hist_tau_search_option) then
            call print_frequency_histograms()

            ! also deallocate here after no use of the histograms anymore
            call deallocate_histograms()
        end if


        if (t_cc_amplitudes .and. t_plot_cc_amplitudes) then
            call print_cc_amplitudes()
        end if

        if (tFValEnergyHist) call print_fval_energy_hist(FvalEnergyHist_EnergyBins, FvalEnergyHist_FValBins)
        if (tFValPopHist) call print_fval_pop_hist(FvalPopHist_PopBins, FvalPopHist_FValBins)

        ! Remove the signal handlers now that there is no way for the
        ! soft-exit part to work
        call clear_signals()

        ! Reduce the iteration count fro the POPSFILE since it is incremented
        ! upon leaving the loop (If done naturally).
        IF (TIncrement) Iter = Iter - 1
        IF (TPopsFile) THEN
            CALL WriteToPopsfileParOneArr(CurrentDets, TotWalkers)

            if (tHDF5TruncPopsWrite) then
                ! If we have already written a file in the last iteration,
                ! we should not write it again
                if (iHDF5TruncPopsIter == 0 .or. (mod(Iter, iHDF5TruncPopsIter) /= 0)) then
                    call write_popsfile_hdf5(iHDF5TruncPopsEx)
                else
                    write(stdout, *)
                    write(stdout, *) "============== Writing Truncated HDF5 popsfile =============="
                    write(stdout, *) "Unnecessary duplication of truncated popsfile is avoided."
                    write(stdout, *) "It has already been written in the last iteration."
                end if
            end if

            if (tPopsProjE) then
                call calc_proje(InstE, AccumE)
                write(stdout, *)
                write(stdout, *) 'Instantaneous projected energy of popsfile:', InstE + Hii
                if (tAccumPopsActive) &
                    write(stdout, *) 'Accumulated projected energy of popsfile:', AccumE + Hii
            end if
        end if

        IF (tCalcFCIMCPsi) THEN
!This routine will actually only print the matrix if tPrintFCIMCPsi is on
            CALL PrintFCIMCPsi()

            IF (tFindCINatOrbs) THEN
                ! This routine takes the wavefunction Psi, calculates the one
                ! electron density matrix, and rotates the HF orbitals to
                ! produce a new ROFCIDUMP file.
                CALL RotateOrbs()
                CALL MPIBarrier(error)
            end if
        end if

        ! If requested, write the most populated states in CurrentDets to a
        ! CORESPACE file, for use in future semi-stochastic calculations.
        if (tWriteCoreEnd) call write_most_pop_core_at_end(write_end_core_size)

        IF (tHistSpawn) CALL WriteHistogram()

        IF (tHistEnergies) CALL WriteHistogramEnergies()

        IF (tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
        end if

        if (t_measure_local_spin) then
            call finalize_local_spin_measurement()
        end if

        if (t_calc_double_occ) then
            ! also output the final estimates from the summed up
            ! variable:
            if (iProcIndex == root) then
                print *, " ===== "
                print *, " Double occupancy from direct measurement: ", &
                    sum_double_occ / (sum_norm_psi_squared * real(StepsSft, dp))
                print *, " ===== "
            end if
            if (t_spin_measurements) then
                call finalize_double_occ_and_spin_diff()
            end if
        end if


        if (tFillingStochRDMonFly .or. tFillingExplicRDMonFly) then
            call finalise_rdms(rdm_definitions, one_rdms, two_rdm_main, two_rdm_recv, &
                               two_rdm_recv_2, en_pert_main, two_rdm_spawn, rdm_estimates, &
                               .false.)
            ! if available, also output the initiator-rdms
            if (tInitsRDM) call finalise_rdms(rdm_inits_defs, inits_one_rdms, two_rdm_inits, &
                                              two_rdm_recv, two_rdm_recv_2, en_pert_main, two_rdm_inits_spawn, inits_estimates, &
                                              .true.)
        end if
        if (tAdaptiveShift) then
            write(stdout, *) "Prefactor of RDM correction due to adaptive shift", RDMCorrectionFactor
        end if

        call PrintHighPops()

        if (tSemiStochastic .and. t_print_core_vec) then
            call print_determ_vec()
            if (tFillingStochRDMonFly) then
                call print_determ_vec_av()
            end if
        end if

        if (t_symmetry_analysis) then
            call analyze_wavefunction_symmetry()
        end if

        !Close open files.
        IF (iProcIndex == Root) THEN
            close(fcimcstats_unit)
            if (inum_runs == 2) close(fcimcstats_unit2)
            IF (tTruncInitiator) close(initiatorstats_unit)
            IF (tLogComplexPops) close(complexstats_unit)
            if (tWritePopsNorm) close(pops_norm_unit)
            if (tLogEXLEVELStats) close(EXLEVELStats_unit)
        end if
        IF (TDebug) close(11)

        if (tHistSpawn) then
            close(Tot_Unique_Dets_Unit)
        end if
        if (tDiagWalkerSubspace) then
            close(unitWalkerDiag)
        end if

        if (tReplicaEstimates .and. iProcIndex == 0) then
            close(replica_est_unit)
        end if

        ! Print out some load balancing stats nicely to end.
        ! n.B. TotWalkers is the number of determinants. (Horrible naming indeed).
        !       TotParts is the number of walkers.
        if (iProcIndex == Root) write(stdout, *) ! linebreak
        if (iProcIndex == Root) then
            write(stdout, '(/,1X,a55)') 'Load balancing information based on the last iteration:'
        end if
        if (iProcIndex == Root) write(stdout, *) ! linebreak
        block
            INTEGER(int64) :: MaxWalkers, MinWalkers

            CALL MPIReduce(TotWalkers, MPI_MAX, MaxWalkers)
            CALL MPIReduce(TotWalkers, MPI_MIN, MinWalkers)
            CALL MPIAllReduce(TotWalkers, MPI_SUM, AllTotWalkers)
            if (iProcIndex == Root) then
                write(stdout, '(1X,a35,1X,f18.10)') 'Mean number of determinants/process:', &
                        real(AllTotWalkers, dp) / real(nNodes, dp)
                write(stdout, '(1X,a34,1X,i18)') 'Min number of determinants/process:', MinWalkers
                write(stdout, '(1X,a34,1X,i18,/)') 'Max number of determinants/process:', MaxWalkers
            end if
        end block
        if (iProcIndex == Root) write(stdout, *) ! linebreak
        block
            real(dp) :: total_n_walkers, min_n_walkers, max_n_walkers

            CALL MPIReduce(sum(abs(TotParts)), MPI_MAX, max_n_walkers)
            CALL MPIReduce(sum(abs(TotParts)), MPI_MIN, min_n_walkers)
            CALL MPIReduce(sum(abs(TotParts)), MPI_SUM, total_n_walkers)
            if (iProcIndex == Root) then
                write(stdout, '(/,1X,a55)') 'Load balancing information based on the last iteration:'
                write(stdout, '(1X,a35,1X,f18.10)') 'Mean number of walkers/process:', &
                        total_n_walkers / real(nNodes, dp)
                write(stdout, '(1X,a34,1X,f18.5)') 'Min number of walkers/process:', min_n_walkers
                write(stdout, '(1X,a34,1X,f18.5,/)') 'Max number of walkers/process:', max_n_walkers
            end if
        end block
        if (iProcIndex == Root) write(stdout, *) ! linebreak


        ! Automatic error analysis.
        call error_analysis(tSinglePartPhase(1), iBlockingIter(1), mean_ProjE_re, ProjE_Err_re, &
                            mean_ProjE_im, ProjE_Err_im, mean_Shift, Shift_Err, tNoProjEValue, tNoShiftValue)

        call MPIBCast(ProjectionE)
        call MPIBCast(mean_ProjE_re)
        call MPIBCast(ProjE_Err_re)
        call MPIBCast(mean_ProjE_im)
        call MPIBCast(ProjE_Err_im)
        call MPIBCast(mean_Shift)
        call MPIBCast(Shift_Err)
        call MPIBCast(tNoProjEValue)
        call MPIBCast(tNoShiftValue)

        if (tTrialWavefunction) then
            allocate(energy_final_output(size(tot_trial_denom)))
            energy_final_output = tot_trial_numerator / tot_trial_denom
        else
            allocate(energy_final_output(size(ProjectionE)))
            energy_final_output = ProjectionE + Hii
        end if

        ! get the value of OutputHii - the offset used for the projected energy
        ! (not necessarily the one used for the Shift!)
        call getProjEOffset()
        iroot = 1
        CALL GetSym(ProjEDet(:, 1), NEl, G1, NBasisMax, RefSym)
        isymh = int(RefSym%Sym%S, sizeof_int) + 1
        write(stdout, '('' Current reference energy'',T52,F19.12)') OutputHii
        if (tNoProjEValue) then
            write(stdout, '('' Projected correlation energy'',T52,F19.12)') real(ProjectionE(1), dp)
            write(stdout, "(A)") " No automatic errorbar obtained for projected energy"
        else
            write(stdout, '('' Projected correlation energy'',T52,F19.12)') mean_ProjE_re
            write(stdout, '('' Estimated error in Projected correlation energy'',T52,F19.12)') ProjE_Err_re
            if (lenof_sign == 2 .and. inum_runs == 1) then
                write(stdout, '('' Projected imaginary energy'',T52,F19.12)') mean_ProjE_im
                write(stdout, '('' Estimated error in Projected imaginary energy'',T52,F19.12)') ProjE_Err_im
            end if
        end if
        if (.not. tNoShiftValue) then
            write(stdout, '('' Shift correlation energy'',T52,F19.12)') mean_Shift
            write(stdout, '('' Estimated error in shift correlation energy'',T52,F19.12)') shift_err
        else
            write(stdout, "(A)") " No reliable averaged shift correlation energy could be obtained automatically"
        end if
        if ((.not. tNoProjEValue) .and. (.not. tNoShiftValue)) then
            !Do shift and projected energy agree?
            write(stdout, "(A)")
            EnergyDiff = abs(mean_Shift - mean_ProjE_re)
            if (EnergyDiff <= sqrt(shift_err**2 + ProjE_Err_re**2)) then
                write(stdout, "(A,F15.8)") " Projected and shift energy estimates agree " &
                    & //"within errorbars: EDiff = ", EnergyDiff
            else if (EnergyDiff <= sqrt((max(shift_err, ProjE_Err_re) * 2)**2 + min(shift_err, ProjE_Err_re)**2)) then
                write(stdout, "(A,F15.8)") " Projected and shift energy estimates agree to within " &
                    & //"two sigma of largest error: EDiff = ", EnergyDiff
            else
                write(stdout, "(A,F15.8)") " Projected and shift energy estimates do not agree to " &
                    & //"within approximate errorbars: EDiff = ", EnergyDiff
            end if
            if (ProjE_Err_re < shift_err) then
                BestEnergy = mean_ProjE_re + OutputHii
                BestErr = ProjE_Err_re
            else
                BestEnergy = mean_shift + Hii
                BestErr = shift_err
            end if
        else if (tNoShiftValue .and. (.not. tNoProjEValue)) then
            BestEnergy = mean_ProjE_re + OutputHii
            BestErr = ProjE_Err_re
        else if (tNoProjEValue .and. (.not. tNoShiftValue)) then
            BestEnergy = mean_shift + Hii
            BestErr = shift_err
        else
            BestEnergy = ProjectionE(1) + OutputHii
            BestErr = 0.0_dp
        end if
        write(stdout, "(A)")
        if (tNoProjEValue) then
            write(stdout, "(A,F20.8)") " Total projected energy ", real(ProjectionE(1), dp) + OutputHii
        else
            write(stdout, "(A,F20.8,A,G15.6)") " Total projected energy ", &
                mean_ProjE_re + OutputHii, " +/- ", ProjE_Err_re
        end if
        if (.not. tNoShiftValue) then
            write(stdout, "(A,F20.8,A,G15.6)") " Total shift energy     ", &
                mean_shift + Hii, " +/- ", shift_err
        end if

        ! Output replica_estimates data for the test suite
        if (tReplicaEstimates) then
            write(stdout, '(/,1x,"THE FOLLOWING IS FOR TEST SUITE USE ONLY!",/)')
            do i = 1, replica_est_len
                write(stdout, '(1x,"REPLICA ESTIMATES FOR STATE",1x,'//int_fmt(i)//',":",)') i

                write(stdout, '(1x,"Variational energy from replica_estimates:",1x,es20.13)') &
                    var_e_num_all(i) / rep_est_overlap_all(i)
                write(stdout, '(1x,"Energy squared from replica_estimates:",1x,es20.13)') &
                    e_squared_num_all(i) / rep_est_overlap_all(i)
                if (tEN2Init) then
                    write(stdout, '(1x,"EN2 estimate from replica_estimates:",1x,es20.13)') &
                        en2_pert_all(i) / rep_est_overlap_all(i)
                end if
                if (tEN2Rigorous) then
                    write(stdout, '(1x,"EN2 New estimate from replica_estimates:",1x,es20.13)') &
                        en2_new_all(i) / rep_est_overlap_all(i)
                end if
                write(stdout, '(1x,"Preconditioned energy from replica_estimates:",1x,es20.13,/)') &
                    precond_e_num_all(i) / precond_denom_all(i)
            end do
        end if

        CALL MolproPluginResult('ENERGY', [BestEnergy])
        CALL MolproPluginResult('FCIQMC_ERR', [min(ProjE_Err_re, shift_err)])
        write(stdout, "(/)")

        ! Deallocate memory
        call DeallocFCIMCMemPar()

    end subroutine FciMCPar

    subroutine PerformFCIMCycPar(iter_data, err)
        use mpi
        use global_det_data, only: get_iter_occ_tot, get_av_sgn_tot
        use global_det_data, only: set_av_sgn_tot, set_iter_occ_tot
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: two_rdm_spawn, two_rdm_recv, two_rdm_main, one_rdms
        use rdm_data, only: rdm_definitions
        use rdm_data_utils, only: communicate_rdm_spawn_t, add_rdm_1_to_rdm_2, clear_rdm_list_t
        use symrandexcit_Ex_Mag, only: test_sym_excit_ExMag
        ! Iteration specific data
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(out) :: err

        ! Now the local, iteration specific, variables
        integer :: j, p, error, proc_temp, i, HFPartInd, isym
        integer :: DetCurr(nel), nJ(nel), FlagsCurr, parent_flags
        real(dp), dimension(lenof_sign) :: SignCurr, child, child_for_stats, SpawnSign
        integer(kind=n_int) :: iLutnJ(0:niftot)
        integer :: IC, walkExcitLevel, walkExcitLevel_toHF, ex(2, 3), TotWalkersNew, part_type, run
        integer(int64) :: tot_parts_tmp(lenof_sign)
        logical :: tParity, tSuccess, tCoreDet(inum_runs), tGlobalCoreDet
        real(dp) :: prob, HDiagCurr, EnergyCurr, hdiag_bare, TempTotParts, Di_Sign_Temp
        real(dp) :: RDMBiasFacCurr
        real(dp) :: lstart
        real(dp) :: AvSignCurr(len_av_sgn_tot), IterRDMStartCurr(len_iter_occ_tot)
        real(dp) :: av_sign(len_av_sgn_tot), iter_occ(len_iter_occ_tot)
        HElement_t(dp) :: HDiagTemp, HElGen
        character(*), parameter :: this_routine = 'PerformFCIMCycPar'
        HElement_t(dp), dimension(inum_runs) :: delta
        integer :: proc, pos, determ_index, irdm
        real(dp) :: r, sgn(lenof_sign), prob_extra_walker
        integer :: DetHash, FinalVal, clash, PartInd, k, y, MaxIndex
        type(ll_node), pointer :: TempNode

        integer :: ms, allErr
        real(dp) :: precond_fac
        ! average number of excitations per walker for a given determinant
        real(dp) :: AvMCExcitsLoc, scale, max_spawn
        HElement_t(dp) :: hdiag_spawn, h_diag_correct

        logical :: signChanged, newlyOccupied, flag_mixed

        real(dp) :: currArg, spawnArg
        integer :: scaleFactor
        ! how many entries were added to (the end of) CurrentDets in the last iteration
        integer, save :: detGrowth = 0

        real(dp) :: inst_rdm_occ

        call set_timer(Walker_Time, 30)
        err = 0
        allErr = 0
        MaxInitPopPos = 0.0_dp
        MaxInitPopNeg = 0.0_dp
        HighPopNeg = 1
        HighPopPos = 1
        FlagsCurr = 0

        ! Reset iteration variables
        ! Next free position in newly spawned list.
        ValidSpawnedList = InitialSpawnedSlots
        FreeSlot(1:iEndFreeSlot) = 0  !Does this cover enough?
        iStartFreeSlot = 1
        iEndFreeSlot = 0

        ! Clear the hash table for the spawning array.
        if (use_spawn_hash_table) call clear_hash_table(spawn_ht)

        ! Index for counting deterministic states.
        determ_index = 0

        call rezero_iter_stats_each_iter(iter_data, rdm_definitions)

        ! The processor with the HF determinant on it will have to check
        ! through each determinant until it's found. Once found, tHFFound is
        ! true, and it no longer needs to be checked.

        precond_fac = 1.0_dp

        IFDEBUGTHEN(FCIMCDebug, stdout)
            write(stdout, "(A)") "Hash Table: "
            write(stdout, "(A)") "Position in hash table, Position in CurrentDets"
            do j = 1, nWalkerHashes
                TempNode => HashIndex(j)
                if (TempNode%Ind /= 0) then
                    write(stdout, '(i9)', advance='no') j
                    do while (associated(TempNode))
                        write(stdout, '(i9)', advance='no') TempNode%Ind
                        TempNode => TempNode%Next
                    end do
                    write(stdout, '()', advance='yes')
                end if
            end do
        ENDIFDEBUG

        IFDEBUG(FCIMCDebug, 3) write(stdout, "(A,I12)") "Walker list length: ", TotWalkers
        IFDEBUG(FCIMCDebug, 3) write(stdout, "(A)") "TW: Walker  Det"

        ! This block decides whether or not to calculate the contribution to the RDMs from
        ! the diagonal elements (and explicit connections to the HF) for each occupied determinant.
        ! For efficiency, this is only done on the final iteration, or one where the RDM energy is
        ! being printed.
        tFill_RDM = .false.
        if (tFillingStochRDMonFly) then
            if (mod((Iter + PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0) then
                ! RDM energy is being printed, calculate the diagonal elements for
                ! the last RDMEnergyIter iterations.
                tFill_RDM = .true.
                IterLastRDMFill = RDMEnergyIter
            else if(Iter == NMCyc .or. Iter - maxval(VaryShiftIter) == eq_cyc) then
                ! Last iteration, calculate the diagonal element for the iterations
                ! since the last time they were included.
                tFill_RDM = .true.
                IterLastRDMFill = mod((Iter + PreviousCycles - IterRDMStart + 1), RDMEnergyIter)
            end if
        end if
        lstart = mpi_wtime()
        loop_over_determinants: do j = 1, int(TotWalkers, sizeof_int)

            ! N.B. j indicates the number of determinants, not the number
            !      of walkers.

            ! Indicate that the scratch storage used for excitation generation
            ! from the same walker has not been filled (it is filled when we
            ! excite from the first particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            ! let the interested excitation generator know about the index in the CurrentDets array.
            fcimc_excit_gen_store%idx_curr_dets = j

            ! Make sure that the parent flags from the last walker don't through.
            parent_flags = 0

            ! If we're not calculating the RDM (or we're calculating some HFSD combination of the
            ! RDM) this just extracts info from the bit representation like normal.
            ! IterRDMStartCurr and AvSignCurr just come out as 1.0_dp.
            ! Otherwise, it extracts the Curr info, and calculates the iteration this determinant
            ! became occupied (IterRDMStartCurr) and the average population during that time
            ! (AvSignCurr).

            ! Is this state is in the deterministic space?
            do run = 1, inum_runs
                tCoreDet(run) = check_determ_flag(CurrentDets(:, j), run)
            end do
            tGlobalCoreDet = all(tCoreDet)

            call extract_bit_rep_avsign(rdm_definitions, CurrentDets(:, j), j, DetCurr, SignCurr, FlagsCurr, &
                                        IterRDMStartCurr, AvSignCurr, fcimc_excit_gen_store)

            ! if we CurrentDets is almost full, make some space
            if (t_activate_decay) then
                if (test_flag(CurrentDets(:, j), flag_prone)) then
                    ! kill a prone determinant with probability given by the
                    ! ratio of how many space we probably need to how many prone dets exist
                    ! (prone dets always have only a single scaled walker)
                    r = genrand_real2_dSFMT()
                    if (n_prone_dets > 0) then
                        if (r < detGrowth / n_prone_dets) then
                            ! log the removal
                            iter_data%nremoved = iter_data%nremoved + abs(SignCurr)
                            ! remove all walkers here (this will be counted as unocc. later
                            ! and thus become an empty slot)
                            SignCurr = 0.0_dp
                            ! kill all walkers on the determinant
                            call nullify_ilut(CurrentDets(:, j))
                            ! and remove it from the hashtable
                            if (.not. tAccumEmptyDet(CurrentDets(:, j))) &
                                call RemoveHashDet(HashIndex, DetCurr, j)
                            cycle
                        end if
                    end if
                end if
            end if

            ! We only need to find out if determinant is connected to the
            ! reference (so no ex. level above 2 required,
            ! truncated etc.)
            walkExcitLevel = FindBitExcitLevel(iLutRef(:, 1), CurrentDets(:, j), &
                                               max_calc_ex_level, .true.)
            if (tRef_Not_HF) then
                walkExcitLevel_toHF = FindBitExcitLevel(iLutHF_true, CurrentDets(:, j), &
                                                        max_calc_ex_level, .true.)
            else
                walkExcitLevel_toHF = walkExcitLevel
            end if

            if (tFillingStochRDMonFly) then
                ! Set the average sign and occupation iteration which were
                ! found in extract_bit_rep_avsign.
                call set_av_sgn_tot(j, AvSignCurr)
                call set_iter_occ_tot(j, IterRDMStartCurr)
                ! If this is an iteration where we print out the RDM energy,
                ! add in the diagonal contribution to the RDM for this
                ! determinant, for each rdm.
                if (tFill_RDM .and. (.not. tNoNewRDMContrib)) then

                    av_sign = get_av_sgn_tot(j)
                    iter_occ = get_iter_occ_tot(j)
                    if (tInitsRDM .and. all_runs_are_initiator(CurrentDets(:, j))) &
                        call fill_rdm_diag_currdet(two_rdm_inits_spawn, inits_one_rdms, &
                                                   CurrentDets(:, j), DetCurr, walkExcitLevel_toHF, av_sign, iter_occ, &
                                                   tGlobalCoreDet, .false.)
                    call fill_rdm_diag_currdet(two_rdm_spawn, one_rdms, CurrentDets(:, j), &
                                               DetCurr, walkExcitLevel_toHF, av_sign, iter_occ, tGlobalCoreDet, tApplyLC)
                end if
            end if

            ! This if-statement is only entered when using semi-stochastic and
            ! only if this determinant is in the core space.
            ! Potential optimization: This list does not change between iterations, only set it up once
            if (allocated(cs_replicas)) then
                do run = 1, size(cs_replicas)
                    associate(rep => cs_replicas(run))
                        if (tCoreDet(run)) then
                            ! Store the index of this state, for use in annihilation later.

                            ! A global core-space is ordered in the same fashion
                            ! in currentdets as internally in the core_space
                            if (t_global_core_space) then
                                determ_index = determ_index + 1
                            else
                                ! In general, that is not true however, the
                                ! core-space index has to be looked up
                                determ_index = core_space_pos(CurrentDets(:, j), DetCurr, &
                                                              run) - rep%determ_displs(iProcIndex)
                            end if

                            rep%indices_of_determ_states(determ_index) = j

                            ! Add this amplitude to the deterministic vector.
                            rep%partial_determ_vecs(:, determ_index) = SignCurr(rep%min_part():rep%max_part())
                        end if
                    end associate
                end do
            end if

            ! As the main list (which is storing a hash table) no longer needs
            ! to be contiguous, we need to skip sites that are empty.
            if (IsUnoccDet(SignCurr)) then
                ! The deterministic states are always kept in CurrentDets, even when
                ! the amplitude is zero. Hence we must check if the amplitude is zero,
                ! and if so, skip the state.
                if (any(tCoreDet) .or. tAccumEmptyDet(CurrentDets(:, j))) cycle

                !It has been removed from the hash table already
                !Just add to the "freeslot" list
                iEndFreeSlot = iEndFreeSlot + 1
                FreeSlot(iEndFreeSlot) = j
                cycle

            end if

            ! sum in (fmu-1)*cmu^2 for the purpose of RDMs
            if (tAdaptiveShift .and. all(.not. tSinglePartPhase) .and. tFillingStochRDMOnFly) then
                ! Only add the contribution from the corespace if it is explicitly demanded
                if ((.not. tGlobalCoreDet) .or. tCoreAdaptiveShift) &
                    call SumCorrectionContrib(SignCurr, j)
            end if

            ! The current diagonal matrix element is stored persistently.
            HDiagCurr = det_diagH(j)
            EnergyCurr = det_diagH(j) + Hii

            if (tSeniorInitiators) then
                SpawnSign = get_all_spawn_pops(j)
                do run = 1, inum_runs
                    if (.not. is_run_unnocc(SignCurr, run) .and. .not. is_run_unnocc(SpawnSign, run)) then
#ifdef CMPLX_
                        !For complex walkers, we consider the sign changed when the argument of the complex
                        !number changes more than pi/2.

                        CurrArg = DATAN2(SignCurr(max_part_type(run)), SignCurr(min_part_type(run)))
                        SpawnArg = DATAN2(SpawnSign(max_part_type(run)), SpawnSign(min_part_type(run)))
                        signChanged = mod(abs(CurrArg - SpawnArg), PI) > PI / 2.0_dp
#else
                        signChanged = SpawnSign(min_part_type(run)) * SignCurr(min_part_type(run)) < 0.0_dp
#endif
                    else
                        signChanged = .false.
                    end if
                    newlyOccupied = is_run_unnocc(SignCurr, run) .and. .not. is_run_unnocc(SpawnSign, run)
                    if (signChanged .or. newlyOccupied) then
                        call reset_tau_int(j, run)
                        call reset_shift_int(j, run)
                        call set_spawn_pop(j, min_part_type(run), SignCurr(min_part_type(run)))
#ifdef CMPLX_
                        call set_spawn_pop(j, max_part_type(run), SignCurr(max_part_type(run)))
#endif
                    else
                        call update_tau_int(j, run, tau)
                        call update_shift_int(j, run, DiagSft(run) * tau)
                    end if
                end do
            end if
            if (tTruncInitiator) then
                call CalcParentFlag(j, DetCurr, WalkExcitLevel, parent_flags)
            end if

            !Debug output.
            IFDEBUGTHEN(FCIMCDebug, 3)
                write(stdout, "(A,I10,a)", advance='no') 'TW:', j, '['
                do part_type = 1, lenof_sign
                    write(stdout, "(f10.5)", advance='no') SignCurr(part_type)
                end do
                write(stdout, '(a,i7)', advance='no') '] ', FlagsCurr
                call WriteBitDet(stdout, CurrentDets(:, j), .true.)
                call neci_flush(stdout)
            ENDIFDEBUG

            if (walkExcitLevel_toHF == 0) HFInd = j

            IFDEBUGTHEN(FCIMCDebug, 1)
                if (j > 1) then
                    if (DetBitEQ(CurrentDets(:, j - 1), CurrentDets(:, j), nifd)) then
                        call stop_all(this_routine, "Shouldn't have the same determinant twice")
                    end if
                end if
            ENDIFDEBUG

            ! in the guga simulation it is probably better to initialize
            ! the csf_irmation out here, so it can be used in the
            ! new way to calculate the reference energy and then the flag
            ! does not have to be checked each time we loop over the walkers..
            if (tGUGA) call fill_csf_i(CurrentDets(0:nifd, j), current_csf_i)

            ! Sum in any energy contribution from the determinant, including
            ! other parameters, such as excitlevel info.
            ! This is where the projected energy is calculated.
            call SumEContrib(DetCurr, WalkExcitLevel, SignCurr, CurrentDets(:, j), HDiagCurr, 1.0_dp, tPairedReplicas, j)

            if (t_calc_double_occ) then
                inst_double_occ = inst_double_occ + &
                                  get_double_occupancy(CurrentDets(:, j), SignCurr)
            end if

            if (t_measure_local_spin) then
                if (tGUGA) then
                    call stop_all(this_routine, "measure_local_spin works only for GUGA")
                else
                    call measure_local_spin(SignCurr, current_csf_i)
                end if
            end if

            if (t_spin_measurements) then
                call measure_double_occ_and_spin_diff(CurrentDets(:, j), &
                                                      DetCurr, SignCurr)
            end if

            ! If we're on the Hartree-Fock, and all singles and doubles are in
            ! the core space, then there will be no stochastic spawning from
            ! this determinant, so we can the rest of this loop.
            if (tSemiStochastic .and. ss_space_in%tDoubles .and. walkExcitLevel_toHF == 0 .and. tDetermHFSpawning) cycle

            ! For HPHFs, get the diagonal Hamiltonian element without the
            ! cross-term correction. This will make calculating the diagonal
            ! elements for new spawnins more efficient
            if (tPreCond .or. tReplicaEstimates) then
                if (tHPHF) then
                    hdiag_bare = get_hdiag_bare_hphf(DetCurr, CurrentDets(:, j), EnergyCurr)
                else
                    hdiag_bare = EnergyCurr
                end if
            end if

            ! Loop over the 'type' of particle.
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
            !                 --> OR double run
            !                 --> part_type == 1, 2; population sets 1 and 2, both real

            ! alis additional idea to skip the number of attempted excitations
            ! for noninititators in the back-spawning approach
            ! remove that for now

            ! try a mixed excitation scheme for guga, where we only do a full
            ! excitation for initiators and the crude one for non-inits
            if (tGen_guga_mixed) then
                if (t_guga_mixed_init) then
                    flag_mixed = any_run_is_initiator(CurrentDets(:, j))

                else if (t_guga_mixed_semi) then
                    if (tSemiStochastic) then
                        flag_mixed = tGlobalCoreDet
                    else
                        flag_mixed = any_run_is_initiator(CurrentDets(:, j))
                    end if
                end if

                if (flag_mixed) then
                    generate_excitation => generate_excitation_guga
                else
                    if (tGen_sym_guga_mol) then
                        generate_excitation => gen_excit_4ind_weighted2

                    else if (tGen_sym_guga_ueg) then
                        if (t_new_real_space_hubbard) then
                            generate_excitation => gen_excit_rs_hubbard
                        else if (t_k_space_hubbard) then
                            generate_excitation => gen_excit_k_space_hub
                        end if
                    end if
                end if
            end if

            if (t_crude_exchange_noninits .or. t_approx_exchange_noninits) then
                is_init_guga = any_run_is_initiator(CurrentDets(:, j))
            end if

            if (t_trunc_guga_pgen_noninits) then
                is_init_guga = any_run_is_initiator(CurrentDets(:, j))
            end if

            if (t_guga_back_spawn) then
                is_init_guga = any_run_is_initiator(CurrentDets(:, j))
            end if

            loop_over_type : do part_type = 1, lenof_sign

                run = part_type_to_run(part_type)
                TempSpawnedPartsInd = 0

                ! Loop over all the particles of a given type on the
                ! determinant. CurrentSign gives number of walkers. Multiply
                ! up by AvMCExcits if attempting multiple excitations from
                ! each walker (default 1.0_dp).

                AvMCExcitsLoc = AvMCExcits
                ! optional: Adjust the number of spawns to the expected maximum
                ! Hij/pgen ratio of this determinant -> prevent blooms
                ! Only done while not updating tau (to prevent interdependencies)
                ! or, for hist-tau-search, in vairable shift mode
                ! Usually, this means: done in variable shift mode
                if (tScaleBlooms .and. .not. tSearchTau &
                    .and. .not. (t_hist_tau_search .and. tSinglePartPhase( &
                                 part_type_to_run(part_type)))) then
                    max_spawn = tau * get_max_ratio(j)
                    if (max_spawn > max_allowed_spawn) then
                        scale = max_spawn / max_allowed_spawn
                        AvMCExcitsLoc = AvMCExcitsLoc * scale
                    end if
                end if
                call decide_num_to_spawn(SignCurr(part_type), AvMCExcitsLoc, WalkersToSpawn)
                loop_over_particles : do p = 1, WalkersToSpawn

                    ! Zero the bit representation, to ensure no extraneous
                    ! data gets through.
                    ilutnJ = 0_n_int
                    child = 0.0_dp

                    ! for the 3-body excitations i really do not want to change
                    ! all the interfaces to the other excitation generators,
                    ! which all just assume ex(2,2) as size.. so use a
                    ! if here..
                    if (t_k_space_hubbard .and. t_3_body_excits) then
                        if (t_uniform_excits) then
                            call gen_excit_uniform_k_space_hub_transcorr(DetCurr, CurrentDets(:, j), &
                                                                         nJ, ilutnJ, exFlag, ic, ex, tParity, prob, &
                                                                         HElGen, fcimc_excit_gen_store, part_type)

                        else if (t_mixed_excits) then
                            call gen_excit_mixed_k_space_hub_transcorr(DetCurr, CurrentDets(:, j), &
                                                                       nJ, ilutnJ, exFlag, ic, ex, tParity, prob, &
                                                                       HElGen, fcimc_excit_gen_store, part_type)

                        else
                            call gen_excit_k_space_hub_transcorr(DetCurr, CurrentDets(:, j), &
                                                                 nJ, ilutnJ, exFlag, ic, ex, tParity, prob, &
                                                                 HElGen, fcimc_excit_gen_store, part_type)
                        end if
                    else
                        ! Generate a (random) excitation
                        call generate_excitation(DetCurr, CurrentDets(:,j), nJ, &
                            ilutnJ, exFlag, IC, ex, tParity, prob, &
                            HElGen, fcimc_excit_gen_store, part_type)

                    end if

                    !If we are fixing the population of reference det, skip spawing into it.
                    if (tSkipRef(run) .and. all(nJ == projEdet(:, run))) then
                        !Set nJ to null
                        nJ(1) = 0
                    end if

                    ! If a valid excitation, see if we should spawn children.
                    if (IsNullDet(nJ)) then
                        nInvalidExcits = nInvalidExcits + 1
                    else
                        nValidExcits = nValidExcits + 1

                        if (tSemiStochastic) then
                            call encode_child(CurrentDets(:, j), iLutnJ, ic, ex)

                            ! Temporary fix: FindExcitBitDet copies the flags of the parent onto the
                            ! child, which causes semi-stochastic simulations to crash. Should it copy
                            ! these flags? There are comments questioning this in create_particle, too.
                            iLutnJ(IlutBits%ind_flag) = 0_n_int

                            ! If the parent state in the core space.
                            if (check_determ_flag(CurrentDets(:, j), run)) then
                                ! Is the spawned state in the core space?
                                tInDetermSpace = is_core_state(iLutnJ, nJ, run)
                                ! If spawning is from and to the core space, cancel it.
                                if (tInDetermSpace) cycle
                                ! Set the flag to specify that the spawning is occuring
                                ! from the core space.
                                call set_flag(iLutnJ, flag_determ_parent)
                            end if

                        end if

                        if (tPreCond .or. tReplicaEstimates) then
                            hdiag_spawn = get_hdiag_from_excit(DetCurr, nJ, &
                                                               iLutnJ, ic, ex, hdiag_bare)

                            if (tPreCond) then
                                precond_fac = hdiag_spawn - proj_e_for_precond(part_type) - &
                                              proje_ref_energy_offsets(part_type) - Hii
                            end if
                        end if

                        child = attempt_create(DetCurr, &
                                               CurrentDets(:, j), SignCurr, &
                                               nJ, iLutnJ, Prob, HElGen, IC, ex, &
                                               tParity, walkExcitLevel, part_type, &
                                               AvSignCurr, AvMCExcitsLoc, RDMBiasFacCurr, precond_fac)
                        ! Note these last two, AvSignCurr and
                        ! RDMBiasFacCurr are not used unless we're
                        ! doing an RDM calculation.
                    end if

                    IFDEBUG(FCIMCDebug, 3) then
                        write(stdout, '(a)', advance='no') 'SP: ['
                        do y = 1, lenof_sign
                            write(stdout, '(f12.5)', advance='no') &
                                child(y)
                        end do
                        write(stdout, '("] ")', advance='no')
                        call write_det(6, nJ, .true.)
                        call neci_flush(stdout)
                    end if

                    ! Children have been chosen to be spawned.
                    is_child_created : if (.not. all(near_zero(child))) then

                        ! Encode child if not done already.
                        if (.not. tSemiStochastic) then
                            call encode_child(CurrentDets(:, j), iLutnJ, ic, ex)
                        end if
                        ! FindExcitBitDet copies the parent flags so that unwanted flags must be unset.
                        ! Should it really do this?
                        if (tTrialWavefunction) then
                            call clr_flag(iLutnJ, flag_trial)
                            call clr_flag(iLutnJ, flag_connected)
                        end if

                        ! If using a preconditioner, update the child weight for statistics
                        ! (mainly for blooms and hence updating the time step).
                        ! Note, the final preconiditoner is applied in annihilation, once
                        ! the exact projected energy is know.
                        child_for_stats = child

                        if (tPreCond) child_for_stats = child_for_stats / precond_fac

                        if (tScaleBlooms) call update_max_ratio(abs(HElGen) / prob, j)

                        call new_child_stats(iter_data, CurrentDets(:, j), &
                                             nJ, iLutnJ, ic, walkExcitLevel, &
                                             child_for_stats, parent_flags, part_type)

                        if (use_spawn_hash_table) then
                            call create_particle_with_hash_table(nJ, ilutnJ, child, part_type, &
                                                                 CurrentDets(:, j), iter_data, err)
                            if (err /= 0) call stop_all(this_routine, 'failure from create_particle_with_hash_table')
                        else
                            call create_particle(nJ, iLutnJ, child, part_type, hdiag_spawn, err, &
                                                 CurrentDets(:, j), SignCurr, p, &
                                                 RDMBiasFacCurr, WalkersToSpawn, abs(HElGen), j)
                            if (err /= 0) call stop_all(this_routine, 'failure from create_particle')
                        end if
                    end if is_child_created ! (child /= 0), Child created.

                end do loop_over_particles ! Cycling over mulitple particles on same determinant.

            end do loop_over_type  ! Cycling over 'type' of particle on a given determinant.

            ! If we are performing a semi-stochastic simulation and this state
            ! is in the deterministic space, then the death step is performed
            ! deterministically later. Otherwise, perform the death step now.
            ! If using a preconditioner, then death is done in the annihilation
            ! routine, after the energy has been calculated.
            if (tDeathBeforeComms) then
                call walker_death(iter_data, DetCurr, CurrentDets(:, j), &
                                  HDiagCurr, SignCurr, j, WalkExcitLevel, &
                                  t_core_die_=.false.)
            end if

        end do loop_over_determinants

        !loop timing for this iteration on this MPI rank
        lt_arr(mod(Iter - 1, 100) + 1) = mpi_wtime() - lstart

        ! if any proc ran out of memory, terminate
        call MPISumAll(err, allErr)
        if (allErr /= 0) then
            err = allErr
            return
        end if

        IFDEBUGTHEN(FCIMCDebug, 2)
            write(stdout, *) 'Finished loop over determinants'
            write(stdout, *) "Holes in list: ", iEndFreeSlot
        ENDIFDEBUG

        if (tSemiStochastic) then
            ! For semi-stochastic calculations only: Gather together the parts
            ! of the deterministic vector stored on each processor, and then
            ! perform the multiplication of the exact projector on this vector.
            if (tDeathBeforeComms) then
                call determ_projection()
            else
                call determ_projection_no_death()
            end if

            if (tFillingStochRDMonFly) then
                ! For RDM calculations, add the current core amplitudes into the
                ! running average.
                call average_determ_vector()
                ! If this is an iteration where the RDM energy is printed then
                ! add the off-diagonal contributions from the core determinants
                ! (the diagonal contributions are done in the same place for
                ! all determinants, regardless of whether they are core or not,
                ! so are not added in here).
                if (tFill_RDM) then
                    call fill_RDM_offdiag_deterministic(rdm_definitions, two_rdm_spawn, one_rdms)
                    ! deterministic space is always only initiators, so it fully counts towards
                    ! the initiator-only RDMs
                    if (tInitsRDM) call fill_RDM_offdiag_deterministic(rdm_inits_defs, &
                                                                       two_rdm_inits_spawn, inits_one_rdms)
                end if
            end if
        end if

        ! With this algorithm, the determinants do not move, and therefore
        ! TotWalkersNew is simply equal to TotWalkers
        TotWalkersNew = int(TotWalkers, sizeof_int)

        ! Update the statistics for the end of an iteration.
        ! Why is this done here - before annihilation!
        call end_iter_stats(TotWalkersNew)

        ! Print bloom/memory warnings
        call end_iteration_print_warn(totWalkersNew)
        call halt_timer(walker_time)

        ! For the direct annihilation algorithm. The newly spawned
        ! walkers should be in a seperate array (SpawnedParts) and the other
        ! list should be ordered.
        call set_timer(annihil_time, 30)
        !HolesInList is returned from direct annihilation with the number of unoccupied determinants in the list
        !They have already been removed from the hash table though.

        call communicate_and_merge_spawns(MaxIndex, iter_data, .false.)

        if (tSimpleInit) then
            call get_ests_from_spawns_simple(MaxIndex, proj_e_for_precond)
        else
            call get_ests_from_spawns(MaxIndex, proj_e_for_precond)
        end if

        if (tSimpleInit) call rm_non_inits_from_spawnedparts(MaxIndex, iter_data)

        ! If performing FCIQMC with preconditioning, then apply the
        ! the preconditioner to the spawnings, and perform the death step.
        if (tPreCond) then
            call set_timer(rescale_time, 30)
            call rescale_spawns(MaxIndex, proj_e_for_precond, iter_data)
            call halt_timer(rescale_time)
        end if

        ! If we haven't performed the death step yet, then do it now.
        if (.not. tDeathBeforeComms) then
            call set_timer(death_time, 30)
            call perform_death_all_walkers(iter_data)
            call halt_timer(death_time)
        end if

        call DirectAnnihilation(TotWalkersNew, MaxIndex, iter_data, err)

        ! if any proc ran out of memory, terminate
        call MPISumAll(err, allErr)
        if (allErr /= 0) then
            err = allErr
            return
        end if
        ! The growth in the size of the occupied part of CurrentDets
        ! this is important for the purpose of prone_walkers
        detGrowth = int(TotWalkersNew - TotWalkers)

        ! This indicates the number of determinants in the list + the number
        ! of holes that have been introduced due to annihilation.
        TotWalkers = TotWalkersNew

        ! if we still have plenty of empty slots in the list, deactivate the decay
        if (t_activate_decay .and. TotWalkers < 0.95_dp * real(MaxWalkersPart, dp)) then
            t_activate_decay = .false.
        end if

        ! The superinitiators are now the same as they will be at the beginning of
        ! the next cycle (this flag is reset if they change)
        tReferenceChanged = .false.

        CALL halt_timer(Annihil_Time)
        IFDEBUG(FCIMCDebug, 2) write(stdout, *) "Finished Annihilation step"

        ! If we are orthogonalising the replica wavefunctions, to generate
        ! excited states, then do that here.
        if (tOrthogonaliseReplicas .and. iter > orthogonalise_iter) then
            call orthogonalise_replicas(iter_data)
        else if (tPrintReplicaOverlaps .and. inum_runs > 1) then
            call calc_replica_overlaps()
        end if

        if (tFillingStochRDMonFly) then
            ! if we use the initiator-only rdms as gamma_0, get them in their own entity
            if (tInitsRDM) call fill_rdm_diag_wrapper(rdm_inits_defs, two_rdm_inits_spawn, &
                                                      inits_one_rdms, CurrentDets, int(TotWalkers, sizeof_int), .false., .false.)
            call fill_rdm_diag_wrapper(rdm_definitions, two_rdm_spawn, one_rdms, &
                                       CurrentDets, int(TotWalkers, sizeof_int), .true., tApplyLC)
        end if

        if (tTrialWavefunction .and. tTrialShift) then
            call fix_trial_overlap(iter_data)
        end if
        call update_iter_data(iter_data)

        ! This routine will take the CurrentDets and search the array to find all single and double
        ! connections - adding them into the RDM's.
        ! This explicit way of doing this is very expensive, but o.k for very small systems.
        if (tFillingExplicRDMonFly) then
            if (tHistSpawn) THEN
                call Fill_Hist_ExplicitRDM_this_Iter()
            else
                call Fill_ExplicitRDM_this_Iter(TotWalkers)
            end if
        end if

        if (tFillingStochRDMonFly .or. tFillingExplicRDMonFly) then
            ! Fill the receiving RDM list from the beginning.
            two_rdm_recv%nelements = 0
            call communicate_rdm_spawn_t(two_rdm_spawn, two_rdm_recv)
            call add_rdm_1_to_rdm_2(two_rdm_recv, two_rdm_main)
            two_rdm_recv%nelements = 0
            if (tInitsRDM) then
                call communicate_rdm_spawn_t(two_rdm_inits_spawn, two_rdm_recv)
                call add_rdm_1_to_rdm_2(two_rdm_recv, two_rdm_inits)
            end if
        end if

    end subroutine PerformFCIMCycPar

END MODULE FciMCParMod
