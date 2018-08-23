#include "macros.h"
module FciMCParMod

    ! This module contains the main loop for FCIMC calculations, and the
    ! main per-iteration processing loop.
    use SystemData, only: nel, tUEG2, hist_spin_dist_iter, tReltvy, tHub, & 
                          t_new_real_space_hubbard, t_tJ_model, t_heisenberg_model, & 
                          t_k_space_hubbard, max_ex_level, t_uniform_excits, &
                          t_mixed_excits

    use CalcData, only: tFTLM, tSpecLanc, tExactSpec, tDetermProj, tMaxBloom, &
                        tUseRealCoeffs, tWritePopsNorm, tExactDiagAllSym, &
                        AvMCExcits, pops_norm_unit, iExitWalkers, &
                        iFullSpaceIter, semistoch_shift_iter, tEN2, &
                        tOrthogonaliseReplicas, orthogonalise_iter, &
                        tDetermHFSpawning, use_spawn_hash_table, &
                        ss_space_in, s_global_start, tContTimeFCIMC, &
                        trial_shift_iter, tStartTrialLater, tAVReps, &
                        tTrialWavefunction, tSemiStochastic, ntrial_ex_calc, &
                        t_hist_tau_search_option, t_back_spawn, back_spawn_delay, &
                        t_back_spawn_flex, t_back_spawn_flex_option, &
                        t_back_spawn_option, tDynamicCoreSpace, coreSpaceUpdateCycle, &
                        DiagSft, tDynamicTrial, trialSpaceUpdateCycle, semistochStartIter, &
                        tSkipRef, tFixTrial, tTrialShift, tSpinProject, tRCCheck

    use adi_data, only: tReadRefs, tDelayGetRefs, allDoubsInitsDelay, tDelayAllSingsInits, &
                        tDelayAllDoubsInits, tDelayAllSingsInits, tReferenceChanged, &
                        SIUpdateInterval, tSuppressSIOutput, nRefUpdateInterval, &
                        SIUpdateOffset

    use LoggingData, only: tJustBlocking, tCompareTrialAmps, tChangeVarsRDM, &
                           tWriteCoreEnd, tNoNewRDMContrib, tPrintPopsDefault,&
                           compare_amps_period, PopsFileTimer, tOldRDMs, &
                           write_end_core_size, t_calc_double_occ, t_calc_double_occ_av, &
                           equi_iter_double_occ, t_print_frq_histograms, &
                           t_spin_measurements

    use spin_project, only: spin_proj_interval, disable_spin_proj_varyshift, &
                            spin_proj_iter_count, generate_excit_spin_proj, &
                            get_spawn_helement_spin_proj, iter_data_spin_proj,&
                            attempt_die_spin_proj

    use rdm_data, only: print_2rdm_est

    use rdm_data_old, only: rdms, one_rdms_old, rdm_estimates_old

    use rdm_finalising, only: finalise_rdms

    use rdm_general, only: init_rdms

    use rdm_general_old, only: InitRDMs_old, FinaliseRDMs_old
    use rdm_filling_old, only: fill_rdm_offdiag_deterministic_old, fill_rdm_diag_wrapper_old
    use rdm_filling, only: fill_rdm_offdiag_deterministic, fill_rdm_diag_wrapper
    use rdm_explicit, only: fill_explicitrdm_this_iter, fill_hist_explicitrdm_this_iter
    use procedure_pointers, only: attempt_die_t, generate_excitation_t, &
                                  get_spawn_helement_t
    use semi_stoch_gen, only: write_most_pop_core_at_end, init_semi_stochastic, &
         refresh_semistochastic_space
    use semi_stoch_procs, only: is_core_state, check_determ_flag, &
                                determ_projection, average_determ_vector
    use trial_wf_gen, only: update_compare_trial_file, init_trial_wf, refresh_trial_wf
    use hist, only: write_zero_hist_excit_tofrom, write_clear_hist_spin_dist
    use orthogonalise, only: orthogonalise_replicas, calc_replica_overlaps, &
                             orthogonalise_replica_pairs
    use load_balance, only: tLoadBalanceBlocks, adjust_load_balance
    use bit_reps, only: set_flag, clr_flag, add_ilut_lists, get_initiator_flag
    use exact_diag, only: perform_exact_diag_all_symmetry
    use spectral_lanczos, only: perform_spectral_lanczos
    use bit_rep_data, only: nOffFlag, flag_determ_parent, test_flag
    use errors, only: standalone_errors, error_analysis
    use PopsFileMod, only: WriteToPopsFileParOneArr
    use AnnihilationMod, only: DirectAnnihilation
    use exact_spectrum, only: get_exact_spectrum
    use determ_proj, only: perform_determ_proj
    use cont_time, only: iterate_cont_time
    use global_det_data, only: det_diagH, reset_tau_int, get_all_spawn_pops, &
                               reset_shift_int, update_shift_int, &
                               update_tau_int, set_spawn_pop
    use RotateOrbsMod, only: RotateOrbs
    use NatOrbsMod, only: PrintOrbOccs
    use ftlm_neci, only: perform_ftlm
    use hash, only: clear_hash_table
    use soft_exit, only: ChangeVars
    use adi_references, only: setup_reference_space, enable_adi, adjust_nRefs, &
         update_reference_space
    use fcimc_initialisation
    use fcimc_iter_utils
    use neci_signals
    use fcimc_helper
    use fcimc_output
    use FciMCData
    use constants

    use real_time_data, only: t_prepare_real_time, n_real_time_copies, &    
                              cnt_real_time_copies    

    use real_time_init, only: init_overlap_buffers

    use bit_reps, only: decode_bit_det    

    use double_occ_mod, only: get_double_occupancy, inst_double_occ, &
                        rezero_double_occ_stats, write_double_occ_stats, & 
                        sum_double_occ, sum_norm_psi_squared, finalize_double_occ_and_spin_diff, &
                        measure_double_occ_and_spin_diff, rezero_spin_diff, &
                        write_spin_diff_stats, write_spat_doub_occ_stats, &
                        all_sum_double_occ, calc_double_occ_from_rdm

    use tau_search_hist, only: print_frequency_histograms, deallocate_histograms
    use back_spawn, only: init_back_spawn
    use real_space_hubbard, only: init_real_space_hubbard
    use tJ_model, only: init_tJ_model, init_heisenberg_model
    use k_space_hubbard, only: init_k_space_hubbard, gen_excit_k_space_hub_transcorr, & 
                               gen_excit_uniform_k_space_hub_transcorr, &
                               gen_excit_mixed_k_space_hub_transcorr
    use cc_amplitudes, only: t_cc_amplitudes, init_cc_amplitudes, cc_delay, &
                            t_plot_cc_amplitudes, print_cc_amplitudes

    use analyse_wf_symmetry, only: analyze_wavefunction_symmetry, t_symmetry_analysis

#ifdef MOLPRO
    use outputResult
#endif

    implicit none

    contains

    subroutine FciMCPar(energy_final_output)

        use rdm_data, only: rdm_estimates, two_rdm_main, two_rdm_recv, two_rdm_recv_2
        use rdm_data, only: two_rdm_spawn, one_rdms, rdm_definitions, en_pert_main
        use rdm_estimators, only: calc_2rdm_estimates_wrapper, write_rdm_estimates
        use rdm_estimators_old, only: rdm_output_wrapper_old, write_rdm_estimates_old

        real(dp), intent(out), allocatable :: energy_final_output(:)

#ifdef MOLPRO
        integer :: nv, ityp(1)
#endif
        integer :: iroot, isymh
        real(dp) :: Weight, Energyxw, BestEnergy
        INTEGER :: error, irdm
        LOGICAL :: TIncrement, tWritePopsFound, tSingBiasChange, tPrintWarn
        REAL(dp) :: s_start, s_end, tstart(2), tend(2), totaltime
        real(dp) :: TotalTime8
        INTEGER(int64) :: MaxWalkers, MinWalkers
        real(dp) :: AllTotWalkers,MeanWalkers, Inpair(2), Outpair(2)
        integer, dimension(lenof_sign) :: tmp_sgn
        integer :: tmp_int(lenof_sign), i, istart, iRDMSamplingIter
        real(dp) :: grow_rate, EnergyDiff, Norm_2RDM
        TYPE(BasisFn) RefSym
        real(dp) :: mean_ProjE_re, mean_ProjE_im, mean_Shift
        real(dp) :: ProjE_Err_re, ProjE_Err_im, Shift_Err
        logical :: tNoProjEValue, tNoShiftValue
        real(dp) :: BestErr
        real(dp) :: start_time, stop_time
#ifdef MOLPRO
        real(dp) :: get_scalar
        include "common/molen"
#endif

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
            write(6,"(A)") "Skipping FCIQMC calculation and simply reblocking previous output"
            call Standalone_Errors()
            return
        endif

        TDebug = .false.  ! Set debugging flag
                    
!OpenMPI does not currently support MPI_Comm_set_errhandler - a bug in its F90 interface code.
!Ask Nick McLaren if we need to change the err handler - he has a fix/bypass.
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,error)
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL,error)
        
        ! This is set here not in SetupParameters, as otherwise it would be
        ! wiped just when we need it!
        tPopsAlreadyRead = .false.

        call SetupParameters()
        call InitFCIMCCalcPar()

        if(tLogGreensfunction .and. .not. t_real_time_fciqmc) then
           call init_overlap_buffers()
        endif

        call init_fcimc_fn_pointers() 

        if (t_new_real_space_hubbard) then 
            call init_real_space_hubbard()
        end if
        if (t_tJ_model) then 
            call init_tJ_model()
        end if
        if (t_heisenberg_model) then 
            call init_heisenberg_model()
        end if
        ! try to call this earlier..
        ! just do it twice for now.. 
        if (t_k_space_hubbard) then 
            call init_k_space_hubbard()
        end if

#ifdef __DEBUG
        call decode_bit_det(tmp_det, ilutHF)
        write(iout, *) "HF: ", tmp_det
        call decode_bit_det(tmp_det, ilutHF_true)
        write(iout, *) "HF_true: ", tmp_det
        call decode_bit_det(tmp_det, ilutRef(:,1))
        write(iout, *) "Ref: ", tmp_det 
        write(iout, *) "ProjEDet: ", ProjEDet
#endif

        ! Attach signal handlers to give a more graceful death-mannerism
        call init_signals()

        ! We want to do some population checking before we run any iterations.
        ! In the normal case this is run between iterations, but it is
        ! helpful to do it here.
        call population_check()

        if(n_int.eq.4) CALL Stop_All('Setup Parameters', &
                'Use of RealCoefficients does not work with 32 bit integers due to the use &
                &of the transfer operation from dp reals to 64 bit integers.')

        if (tDetermProj) then
            ! If performing a deterministic projection instead of an FCIQMC calc:
            call perform_determ_proj()
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
            !write (iout,'("#")', advance='no')
            if (iProcIndex == root) then
                write (fcimcstats_unit,'("#")', advance='no')
                if (inum_runs == 2) &
                    write(fcimcstats_unit2, '("#")', advance='no')
                write (initiatorstats_unit,'("#")', advance='no')
                if (tLogEXLEVELStats) &
                      write(EXLEVELStats_unit,'("#")', advance='no')
            end if
            call WriteFCIMCStats()
        end if

        ! double occupancy: 
        if (t_calc_double_occ) then 
            call write_double_occ_stats(iter_data_fciqmc, initial = .true.)
            call write_double_occ_stats(iter_data_fciqmc)
            if (t_spin_measurements) then

                call write_spin_diff_stats(iter_data_fciqmc, initial = .true.)
                call write_spin_diff_stats(iter_data_fciqmc)

                call write_spat_doub_occ_stats(iter_data_fciqmc, initial = .true.)
                call write_spat_doub_occ_stats(iter_data_fciqmc)

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
        Iter=1
        iRDMSamplingIter = 1    !For how many iterations have we accumulated the RDM

        SumSigns = 0.0_dp
        SumSpawns = 0.0_dp

        ! In we go - start the timer for scaling curve!
        start_time = neci_etime(tstart)

        do while (.true.)
!Main iteration loop...
            if(TestMCExit(Iter,iRDMSamplingIter)) then
               ! The popsfile requires the right total walker number, so 
               ! update it (TotParts is updated in the annihilation step)
               call MPISumAll(TotParts,AllTotParts)
               exit
            endif

            IFDEBUG(FCIMCDebug, 2) write(iout,*) 'Iter', iter

            if(iProcIndex.eq.root) s_start=neci_etime(tstart)

            ! Update the semistochastic space if requested
            if(tSemiStochastic .and. tDynamicCoreSpace .and. &
                 mod(iter-semistochStartIter, &
                 coreSpaceUpdateCycle) == 0) then
               call refresh_semistochastic_space()
               write(6,*) "Refereshing semistochastic space at iteration ", iter
            end if

            ! Is this an iteration where semi-stochastic is turned on?
            if (semistoch_shift_iter /= 0 .and. all(.not. tSinglePartPhase)) then
                if ((Iter - maxval(VaryShiftIter)) == semistoch_shift_iter + 1) then
                    tSemiStochastic = .true.
                    call init_semi_stochastic(ss_space_in)
                    ! Count iterations for corespace updates from here
                    semistochStartIter = iter
                    ! and switch how iterations for SI updates are counted
                    SIUpdateOffset = semistochStartIter
                end if
            end if

            ! Update the trial wavefunction if requested
            if(tTrialWavefunction .and. tDynamicTrial .and. &
                 mod(iter - trial_shift_iter, trialSpaceUpdateCycle) == 0) then
               if(tPairedReplicas) then
                  call refresh_trial_wf(trial_space_in,ntrial_ex_calc, &
                       inum_runs/2,.true.)
               else
                  call refresh_trial_wf(trial_space_in,ntrial_ex_calc, &
                       inum_runs,.false.)
               endif
               write(6,*) "Refreshing trial wavefunction at iteration ", iter
            endif
       
            if(((Iter - maxval(VaryShiftIter)) == allDoubsInitsDelay + 1 &
                 .and. all(.not. tSinglePartPhase))) then
               ! Start the all-doubs-initiator procedure
               if(tDelayAllDoubsInits) call enable_adi()
               ! And/or the all-sings-initiator procedure
               if(tDelayAllSingsInits) then 
                  tAllSingsInitiators = .true.
                  tAdiActive = .true.
               endif
               ! If desired, we now set up the references for the purpose of the
               ! all-doubs-initiators
               if(tDelayGetRefs) then 
                  ! Re-initialize the reference space
                  call setup_reference_space(.true.)
               endif
            endif

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
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs/2, .true.)
                    else
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs, .false.)
                    end if
                end if
            end if
            
            if(tRDMonFly .and. (.not. tFillingExplicRDMonFly) &
                & .and. (.not.tFillingStochRDMonFly)) call check_start_rdm()

            if (tContTimeFCIMC) then
                call iterate_cont_time(iter_data_fciqmc)
            else
                if (.not. (tSpinProject .and. spin_proj_interval == -1)) &
                    call PerformFciMCycPar(iter_data_fciqmc)
            end if

            ! Are we projecting the spin out between iterations?
            if (tSpinProject .and. (mod(Iter, spin_proj_interval) == 0 .or. &
                                    spin_proj_interval == -1) .and. &
                (tSinglePartPhase(1) .or. .not. disable_spin_proj_varyshift))then

                ! Set this up for a different type of iteration
                ge_tmp => generate_excitation
                gs_tmp => get_spawn_helement
                ad_tmp => attempt_die
                generate_excitation => generate_excit_spin_proj
                get_spawn_helement => get_spawn_helement_spin_proj
                attempt_die => attempt_die_spin_proj

                do i = 1, max(spin_proj_iter_count, 1)
                    call PerformFciMCycPar (iter_data_spin_proj)
                enddo

                ! Return to prior config
                generate_excitation => ge_tmp
                get_spawn_helement => gs_tmp
                attempt_die => ad_tmp
            endif

            if(iProcIndex.eq.root) then
                s_end=neci_etime(tend)
                IterTime=IterTime+(s_end-s_start)
            endif

            ! Add some load balancing magic!
            if (tLoadBalanceBlocks .and. mod(iter, 1000) == 0 .and. &
                .not. tSemiStochastic .and. .not. tFillingStochRDMOnFly) then
                call adjust_load_balance(iter_data_fciqmc)
            end if
            
            if(SIUpdateInterval > 0) then
               ! Regular update of the superinitiators. Use with care as it 
               ! is still rather expensive if secondary superinitiators are used
               if(mod(iter-SIUpdateOffset,SIUpdateInterval) == 0) then
                  ! the reference population needs some time to equilibrate
                  ! hence, nRefs cannot be updated that often
                  if(mod(iter,nRefUpdateInterval) == 0) call adjust_nRefs()
                  call update_reference_space(tReadPops .or. all(.not. tSinglePartPhase))
                  endif
            endif

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

                                write(iout, bloom_warn_string) &
                                    trim(excit_descriptor(i)), bloom_sizes(i),&
                                    bloom_count(i)

                                if (tUEG2) &
                                    call stop_all(this_routine, "Should never &
                                                 &bloom in UEG2")

                            end if
                        end if
                    end do
                endif
                ! Zero the accumulators
                bloom_sizes = 0
                bloom_count = 0
                ! Calculate the a new value for the shift (amongst other
                ! things). Generally, collate information from all processors,
                ! update statistics and output them to the user.
                call set_timer(Stats_Comms_Time)
                call calculate_new_shift_wrapper (iter_data_fciqmc, TotParts, tPairedReplicas)
                call halt_timer(Stats_Comms_Time)

                ! in calculate_new_shift_wrapper output is plotted too! 
                ! so for now do it here for double occupancy
                if (t_calc_double_occ) then 
                    call write_double_occ_stats(iter_data_fciqmc)
                    if (t_spin_measurements) then
                        call write_spin_diff_stats(iter_data_fciqmc)
                        call write_spat_doub_occ_stats(iter_data_fciqmc)
                    end if
                end if
                if(tRestart) cycle

                IF((tTruncCAS.or.tTruncSpace.or.tTruncInitiator).and.(Iter.gt.iFullSpaceIter)&
                            .and.(iFullSpaceIter.ne.0)) THEN
!Test if we want to expand to the full space if an EXPANDSPACE variable has been set
                    IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
                        IF(iProcIndex.eq.0) WRITE(iout,*) "Unable to expand space since histgramming the wavefunction..."
                    ELSE
                        ICILevel=0
                        tTruncSpace=.false.
                        tTruncCAS=.false.
                        IF(tTruncInitiator) tTruncInitiator=.false.
                        IF(iProcIndex.eq.0) THEN
                            WRITE(iout,*) "Expanding to the full space on iteration ",Iter + PreviousCycles
                        ENDIF
                    ENDIF
                ENDIF

                if(iProcIndex.eq.root) &
                    TotalTime8 = real(s_end - s_global_start, dp)
                call MPIBCast(TotalTime8)    !TotalTime is local - broadcast to all procs

!This routine will check for a CHANGEVARS file and change the parameters of the calculation accordingly.
                CALL ChangeVars(tSingBiasChange, tWritePopsFound)
                IF(tSoftExitFound) THEN
                    !Now changed such that we do one more block of iterations and then exit, to allow for proper
                    !inclusion of all the RDM elements
                    NMCyc=Iter+StepsSft  
                    tSoftExitFound = .false.
                    
                    !TIncrement=.false.
                    !! The diagonal elements of the RDM will not have been calculated (as we didn't know 
                    !! it was the last iteration), so this must be done now.
                    !It's problematic trying to det the core space off-diags done here too, hence my change
                    !to the way softexit is handled.
                    !if(tFillingStochRDMonFly) call fill_rdm_softexit(TotWalkers)
                    !EXIT
                ENDIF
                IF(tTimeExit.and.(TotalTime8.ge.MaxTimeExit)) THEN
                    !Is it time to exit yet?
                    write(iout,"(A,F8.2,A)") "Time limit reached for simulation of: ",MaxTimeExit/60.0_dp," minutes - exiting..."
                    NMCyc=Iter+StepsSft  
                    ! Set this to false so that this if statement won't be entered next time.
                    tTimeExit = .false.
                    !tIncrement=.false.
                    !if(tFillingStochRDMonFly) call fill_rdm_softexit(TotWalkers)
                    !EXIT
                ENDIF
                IF((iExitWalkers.ne.-1.0_dp).and.(sum(AllTotParts).gt.iExitWalkers)) THEN
                    !Exit criterion based on total walker number met.
                    write(iout,"(A,I15)") "Total walker population exceeds that given by &
                        &EXITWALKERS criteria - exiting...",sum(AllTotParts)
                    tIncrement=.false.
                    exit
                ENDIF
                IF(tWritePopsFound) THEN
!We have explicitly asked to write out the POPSFILE from the CHANGEVARS file.
                    CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                ENDIF
                IF(tSingBiasChange) THEN
                    CALL CalcApproxpDoubles()
                ENDIF

                if((PopsfileTimer.gt.0.0_dp).and.((iPopsTimers*PopsfileTimer).lt.(TotalTime8/3600.0_dp))) then
                    !Write out a POPSFILE every PopsfileTimer hours
                    if(iProcIndex.eq.Root) then
                        CALL RENAME('popsfile.h5','popsfile.h5.bk')
                        CALL RENAME('POPSFILEBIN','POPSFILEBIN.bk')
                        CALL RENAME('POPSFILEHEAD','POPSFILEHEAD.bk')
                        write(iout,"(A,F7.3,A)") "Writing out a popsfile after ",iPopsTimers*PopsfileTimer, " hours..."
                    endif
                    call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                    iPopsTimers=iPopsTimers+1
                    if(iProcIndex.eq.Root) then
                        s_end=neci_etime(tend)
                        write(iout,"(A,F7.3,A)") "Time taken to write out POPSFILE: ",real(s_end,dp)-TotalTime8," seconds."
                    endif
                endif
                if (tHistExcitToFrom) &
                    call write_zero_hist_excit_tofrom()

            ENDIF   !Endif end of update cycle

            if (tHistSpinDist .and. (mod(iter, hist_spin_dist_iter) == 0)) then
                if (inum_runs.eq.2) then
                    !COMPLEX
                    call stop_all(this_routine, "Not set up to combine HistSpinDist with double runs. &
                                    & Level of changes required to get this working: unknown.")
                endif
                call write_clear_hist_spin_dist (iter, hist_spin_dist_iter)
            endif

            IF(TPopsFile.and.(.not.tPrintPopsDefault).and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
                ! differentiate between normal routine and the real-time 
                ! preperation 
                if (t_prepare_real_time) then
                    if (cnt_real_time_copies < n_real_time_copies) then
                        CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                        cnt_real_time_copies = cnt_real_time_copies + 1
                    else 
                        ! if the number of wanted copies is reached exit the 
                        ! calculation. REMINDER: at the end of the calc. an
                        ! additional popsfile is printed -> so start count at 1
                        ! todo: cleanly exit calc. here! 
                        NMCyc=Iter+StepsSft  
                        t_prepare_real_time = .false.
                        ! do want to finish all the rdm stuff to go on, but 
                        ! do not want anymore popsfile to be printed except
                        ! at the very end.. 
                        tPrintPopsDefault = .true.
                    end if
                else
                    !This will write out the POPSFILE if wanted
                    CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                end if
            ENDIF
!            IF(TAutoCorr) CALL WriteHistogrammedDets()

            IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0).and.(.not.tRDMonFly)) THEN
                CALL WriteHistogram()
            ENDIF

            if (tRDMonFly .and. all(.not. tSinglePartPhase)) then
                ! If we wish to calculate the energy, have started accumulating the RDMs, 
                ! and this is an iteration where the energy should be calculated, do so.
                if (print_2rdm_est .and. ((Iter - maxval(VaryShiftIter)) > IterRDMonFly) &
                    .and. (mod(Iter+PreviousCycles-IterRDMStart+1, RDMEnergyIter) == 0) ) then
           
                    call calc_2rdm_estimates_wrapper(rdm_definitions, rdm_estimates, two_rdm_main, en_pert_main)
                    if (tOldRDMs) then
                        do irdm = 1, rdm_definitions%nrdms
                            call rdm_output_wrapper_old(rdms(irdm), one_rdms_old(irdm), irdm, rdm_estimates_old(irdm), .false.)
                        end do
                    end if

                    if (iProcIndex == 0) then
                        call write_rdm_estimates(rdm_definitions, rdm_estimates, .false., print_2rdm_est)
                        if (tOldRDMs) call write_rdm_estimates_old(rdm_estimates_old, .false.)
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
                if (tOldRDMs) call InitRDMs_old(rdm_definitions%nrdms)
                tRDMonFly = .true.
                tChangeVarsRDM = .false.
            endif

            if (tDiagWalkerSubspace.and.(mod(Iter,iDiagSubspaceIter).eq.0)) then
                ! Diagonalise a subspace consisting of the occupied determinants
                call DiagWalkerSubspace()
            endif

            ! If requested and on a correct iteration, update the COMPARETRIAL file.
            if (tCompareTrialAmps) then
                ASSERT(compare_amps_period /= 0)
                if (mod(Iter, compare_amps_period) == 0) then
                    call update_compare_trial_file(.false.)
                endif
            endif

            Iter=Iter+1
            if(tFillingStochRDMonFly) iRDMSamplingIter = iRDMSamplingIter + 1 

        ! End of MC cycle
        end do

        ! Final output is always enabled
        tSuppressSIOutput = .false.

        ! We are at the end - get the stop-time. Output the timing details
        stop_time = neci_etime(tend)
        write(iout,*) '- - - - - - - - - - - - - - - - - - - - - - - -'
        write(iout,*) 'Total loop-time: ', stop_time - start_time
        write(iout,*) '- - - - - - - - - - - - - - - - - - - - - - - -'

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

        ! Remove the signal handlers now that there is no way for the
        ! soft-exit part to work
        call clear_signals()

        ! Reduce the iteration count fro the POPSFILE since it is incremented
        ! upon leaving the loop (If done naturally).
        IF(TIncrement) Iter=Iter-1
        IF(TPopsFile) THEN
            CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
        ENDIF
        IF(tCalcFCIMCPsi) THEN
!This routine will actually only print the matrix if tPrintFCIMCPsi is on
            CALL PrintFCIMCPsi()

            IF(tFindCINatOrbs) THEN
               ! This routine takes the wavefunction Psi, calculates the one 
               ! electron density matrix, and rotates the HF orbitals to 
               ! produce a new ROFCIDUMP file.
                CALL RotateOrbs() 
                CALL MPIBarrier(error)
            ENDIF
        ENDIF

        ! mswalkercounts
!        if (tReltvy) then
!            call writeMsWalkerCountsAndCloseUnit()
!        endif

        ! If requested, write the most populated states in CurrentDets to a
        ! CORESPACE file, for use in future semi-stochastic calculations.
        if (tWriteCoreEnd) call write_most_pop_core_at_end(write_end_core_size)

        IF(tHistSpawn) CALL WriteHistogram()


        IF(tHistEnergies) CALL WriteHistogramEnergies()

        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
        ENDIF

        if (t_calc_double_occ) then 
            ! also output the final estimates from the summed up 
            ! variable: 
            if (iProcIndex == root) then
                print *, " ===== " 
                print *, " Double occupancy from direct measurement: ", & 
                    sum_double_occ / (sum_norm_psi_squared * real(StepsSft,dp))
                print *, " ===== "
            end if
            if (t_spin_measurements) then
                call finalize_double_occ_and_spin_diff()
            end if
        end if

        if (tFillingStochRDMonFly .or. tFillingExplicRDMonFly) then
            call finalise_rdms(rdm_definitions, one_rdms, two_rdm_main, two_rdm_recv, &
                               two_rdm_recv_2, en_pert_main, two_rdm_spawn, rdm_estimates)
            if (tOldRDMs) call FinaliseRDMs_old(rdms, one_rdms_old, rdm_estimates_old)
        end if

        call PrintHighPops()

        if (t_symmetry_analysis) then
            call analyze_wavefunction_symmetry()
        end if

        !Close open files.
        IF(iProcIndex.eq.Root) THEN
            CLOSE(fcimcstats_unit)
            if (inum_runs.eq.2) CLOSE(fcimcstats_unit2)
            IF(tTruncInitiator) CLOSE(initiatorstats_unit)
            IF(tLogComplexPops) CLOSE(complexstats_unit)
            if (tWritePopsNorm) close(pops_norm_unit)
            if (tLogEXLEVELStats) close(EXLEVELStats_unit)
        ENDIF
        IF(TDebug) CLOSE(11)

        if(tHistSpawn) then 
            close(Tot_Unique_Dets_Unit)
        endif
        if(tDiagWalkerSubspace) then
            close(unitWalkerDiag)
        endif

        ! Print out some load balancing stats nicely to end.
        CALL MPIReduce(TotWalkers,MPI_MAX,MaxWalkers)
        CALL MPIReduce(TotWalkers,MPI_MIN,MinWalkers)
        CALL MPIAllReduce(Real(TotWalkers,dp),MPI_SUM,AllTotWalkers)
        if (iProcIndex.eq.Root) then
            MeanWalkers=AllTotWalkers/nNodes
            write (iout,'(/,1X,a55)') 'Load balancing information based on the last iteration:'
            write (iout,'(1X,a35,1X,f18.10)') 'Mean number of determinants/process:',MeanWalkers
            write (iout,'(1X,a34,1X,i18)') 'Min number of determinants/process:',MinWalkers
            write (iout,'(1X,a34,1X,i18,/)') 'Max number of determinants/process:',MaxWalkers
        end if
       
        ! Automatic error analysis.
        call error_analysis(tSinglePartPhase(1),iBlockingIter(1),mean_ProjE_re,ProjE_Err_re,  &
            mean_ProjE_im,ProjE_Err_im,mean_Shift,Shift_Err,tNoProjEValue,tNoShiftValue)

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
            energy_final_output = tot_trial_numerator/tot_trial_denom
        else
            allocate(energy_final_output(size(ProjectionE)))
            energy_final_output = ProjectionE + Hii
        end if
        
        iroot=1
        CALL GetSym(ProjEDet(:,1),NEl,G1,NBasisMax,RefSym)
        isymh=int(RefSym%Sym%S,sizeof_int)+1
        write (iout,'('' Current reference energy'',T52,F19.12)') Hii 
        if(tNoProjEValue) then
            write (iout,'('' Projected correlation energy'',T52,F19.12)') real(ProjectionE(1),dp)
            write (iout,"(A)") " No automatic errorbar obtained for projected energy"
        else
            write (iout,'('' Projected correlation energy'',T52,F19.12)') mean_ProjE_re
            write (iout,'('' Estimated error in Projected correlation energy'',T52,F19.12)') ProjE_Err_re
            if(lenof_sign.eq.2) then
                write (iout,'('' Projected imaginary energy'',T52,F19.12)') mean_ProjE_im
                write (iout,'('' Estimated error in Projected imaginary energy'',T52,F19.12)') ProjE_Err_im
            endif
        endif
        if(.not.tNoShiftValue) then
            write (iout,'('' Shift correlation energy'',T52,F19.12)') mean_Shift
            write (iout,'('' Estimated error in shift correlation energy'',T52,F19.12)') shift_err
        else
            write(6,"(A)") " No reliable averaged shift correlation energy could be obtained automatically"
        endif
        if((.not.tNoProjEValue).and.(.not.tNoShiftValue)) then
           !Do shift and projected energy agree?
            write(iout,"(A)")
            EnergyDiff = abs(mean_Shift-mean_ProjE_re)
            if(EnergyDiff.le.sqrt(shift_err**2+ProjE_Err_re**2)) then
                write(iout,"(A,F15.8)") " Projected and shift energy estimates agree " &
                    & //"within errorbars: EDiff = ",EnergyDiff
            elseif(EnergyDiff.le.sqrt((max(shift_err,ProjE_Err_re)*2)**2+min(shift_err,ProjE_Err_re)**2)) then
                write(iout,"(A,F15.8)") " Projected and shift energy estimates agree to within " &
                    & //"two sigma of largest error: EDiff = ",EnergyDiff
            else
                write(iout,"(A,F15.8)") " Projected and shift energy estimates do not agree to " &
                    & //"within approximate errorbars: EDiff = ",EnergyDiff
            endif
            if(ProjE_Err_re.lt.shift_err) then
                BestEnergy = mean_ProjE_re + Hii
                BestErr = ProjE_Err_re
            else
                BestEnergy = mean_shift + Hii
                BestErr = shift_err
            endif
        elseif(tNoShiftValue.and.(.not.tNoProjEValue)) then
            BestEnergy = mean_ProjE_re + Hii
            BestErr = ProjE_Err_re
        elseif(tNoProjEValue.and.(.not.tNoShiftValue)) then
            BestEnergy = mean_shift + Hii
            BestErr = shift_err 
        else
            BestEnergy = ProjectionE(1) + Hii
            BestErr = 0.0_dp
        endif
        write(iout,"(A)")
        if(tNoProjEValue) then
            write(iout,"(A,F20.8)") " Total projected energy ",real(ProjectionE(1),dp) + Hii
        else
            write(iout,"(A,F20.8,A,G15.6)") " Total projected energy ", &
                mean_ProjE_re+Hii," +/- ",ProjE_Err_re
        endif
        if(.not.tNoShiftValue) then
            write(iout,"(A,F20.8,A,G15.6)") " Total shift energy     ", &
                mean_shift+Hii," +/- ",shift_err
        endif
#ifdef MOLPRO
        call output_result('FCIQMC','ENERGY',BestEnergy,iroot,isymh)
        if (iroot.eq.1) call clearvar('ENERGY')
        ityp(1)=1
        call setvar('ENERGY',BestEnergy,'AU',ityp,1,nv,iroot)
        do i=10,2,-1
            gesnam(i)=gesnam(i-1)
            energ(i)=energ(i-1)
        enddo
        gesnam(i) = 'FCIQMC'
        energ(i) = get_scalar("ENERGY")
        if(.not.(tNoShiftValue.and.tNoProjEValue)) then
            call output_result('FCIQMC','FCIQMC_ERR',BestErr,iroot,isymh)
            if (iroot.eq.1) call clearvar('FCIQMC_ERR')
            call setvar('FCIQMC_ERR',BestErr,'AU',ityp,1,nv,iroot)
!            do i=10,2,-1
!                gesnam(i)=gesnam(i-1)
!                energ(i)=energ(i-1)
!            enddo
!            gesnam(i) = 'FCIQMC_ERR'
!            energ(i) = get_scalar("FCIQMC_ERR")
        endif
#endif
        write(iout,"(/)")

        ! Deallocate memory
        call DeallocFCIMCMemPar()

    end subroutine FciMCPar

    subroutine PerformFCIMCycPar(iter_data)

        use global_det_data, only: get_iter_occ_tot, get_av_sgn_tot
        use global_det_data, only: set_av_sgn_tot, set_iter_occ_tot
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: two_rdm_spawn, two_rdm_recv, two_rdm_main, one_rdms
        use rdm_data, only: rdm_definitions, rdm_estimates
        use rdm_data_utils, only: communicate_rdm_spawn_t, add_rdm_1_to_rdm_2
        use symrandexcit_Ex_Mag, only: test_sym_excit_ExMag 
        ! Iteration specific data
        type(fcimc_iter_data), intent(inout) :: iter_data

        ! Now the local, iteration specific, variables
        integer :: j, p, error, proc_temp, i, HFPartInd,isym
        integer :: DetCurr(nel), nJ(nel), FlagsCurr, parent_flags
        real(dp), dimension(lenof_sign) :: SignCurr, child, SpawnSign
        integer(kind=n_int) :: iLutnJ(0:niftot)
        integer :: IC, walkExcitLevel, walkExcitLevel_toHF, ex(2,3), TotWalkersNew, part_type, run
        integer(int64) :: tot_parts_tmp(lenof_sign)
        logical :: tParity, tSuccess, tCoreDet
        real(dp) :: prob, HDiagCurr, TempTotParts, Di_Sign_Temp
        real(dp) :: RDMBiasFacCurr
        real(dp) :: AvSignCurr(len_av_sgn_tot), IterRDMStartCurr(len_iter_occ_tot)
        real(dp) :: av_sign(len_av_sgn_tot), iter_occ(len_iter_occ_tot)
        HElement_t(dp) :: HDiagTemp,HElGen
        character(*), parameter :: this_routine = 'PerformFCIMCycPar' 
        HElement_t(dp), dimension(inum_runs) :: delta
        integer :: proc, pos, determ_index, irdm
        real(dp) :: r, sgn(lenof_sign), prob_extra_walker
        integer :: DetHash, FinalVal, clash, PartInd, k, y
        type(ll_node), pointer :: TempNode

        integer :: ms
        logical :: signChanged, newlyOccupied
        real(dp) :: currArg, spawnArg

        real(dp) :: inst_rdm_occ

        call set_timer(Walker_Time,30)

        MaxInitPopPos=0.0_dp
        MaxInitPopNeg=0.0_dp
        HighPopNeg=1
        HighPopPos=1
        FlagsCurr=0
        ! Synchronise processors
!        CALL MPIBarrier(error)

        ! Reset iteration variables
        ! Next free position in newly spawned list.
        ValidSpawnedList = InitialSpawnedSlots
        FreeSlot(1:iEndFreeSlot)=0  !Does this cover enough?
        iStartFreeSlot=1
        iEndFreeSlot=0
       
        ! Clear the hash table for the spawning array.
        if (use_spawn_hash_table) call clear_hash_table(spawn_ht)

        ! Index for counting deterministic states.
        determ_index = 1
        
        call rezero_iter_stats_each_iter(iter_data, rdm_definitions)

        ! [W.D] i should not rezero in here or?
        ! otherwise i waste calculated stuff..
        ! quick and dirty double occupancy measurement: 
!         if (t_calc_double_occ) then 
!             call rezero_double_occ_stats()
!         end if

!         if (t_inst_spin_diff) then 
!             call rezero_spin_diff()
!         end if

        ! The processor with the HF determinant on it will have to check 
        ! through each determinant until it's found. Once found, tHFFound is
        ! true, and it no longer needs to be checked.

        ! This is a bit of a hack based on the fact that we mean something 
        ! different by exFlag for CSFs than in normal determinential code.
        ! It would be nice to fix this properly
        if (tCSF) exFlag = 7

        IFDEBUGTHEN(FCIMCDebug,iout)
            write(iout,"(A)") "Hash Table: "
            write(iout,"(A)") "Position in hash table, Position in CurrentDets"
            do j=1,nWalkerHashes
                TempNode => HashIndex(j)
                if (TempNode%Ind /= 0) then
                    write(iout,'(i9)',advance='no') j
                    do while (associated(TempNode))
                        write(iout,'(i9)',advance='no') TempNode%Ind
                        TempNode => TempNode%Next
                    end do
                    write(iout,'()',advance='yes')
                end if
            end do
        ENDIFDEBUG

        IFDEBUG(FCIMCDebug,3) write(iout,"(A,I12)") "Walker list length: ",TotWalkers
        IFDEBUG(FCIMCDebug,3) write(iout,"(A)") "TW: Walker  Det"

        ! This block decides whether or not to calculate the contribution to the RDMs from 
        ! the diagonal elements (and explicit connections to the HF) for each occupied determinant. 
        ! For efficiency, this is only done on the final iteration, or one where the RDM energy is 
        ! being printed.
        tFill_RDM = .false.
        if(tFillingStochRDMonFly) then
            if(mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter).eq.0) then 
                ! RDM energy is being printed, calculate the diagonal elements for 
                ! the last RDMEnergyIter iterations.
                tFill_RDM = .true.
                IterLastRDMFill = RDMEnergyIter
            elseif(Iter.eq.NMCyc) then
                ! Last iteration, calculate the diagonal element for the iterations 
                ! since the last time they were included.
                tFill_RDM = .true.
                IterLastRDMFill = mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter)
            endif
        endif

        do j = 1, int(TotWalkers,sizeof_int)

            ! N.B. j indicates the number of determinants, not the number
            !      of walkers.

            ! reset this flag for each det:
            ! W.D. remove this option for now..

            ! Indicate that the scratch storage used for excitation generation
            ! from the same walker has not been filled (it is filled when we
            ! excite from the first particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            ! Make sure that the parent flags from the last walker don't through.
            parent_flags=0

            ! If we're not calculating the RDM (or we're calculating some HFSD combination of the 
            ! RDM) this just extracts info from the bit representation like normal.
            ! IterRDMStartCurr and AvSignCurr just come out as 1.0_dp.  
            ! Otherwise, it extracts the Curr info, and calculates the iteration this determinant 
            ! became occupied (IterRDMStartCurr) and the average population during that time 
            ! (AvSignCurr).

            ! Is this state is in the deterministic space?
            tCoreDet = check_determ_flag(CurrentDets(:,j))

            call extract_bit_rep_avsign(rdm_definitions, CurrentDets(:,j), j, DetCurr, SignCurr, FlagsCurr, &
                                        IterRDMStartCurr, AvSignCurr, fcimc_excit_gen_store)

            !call test_sym_excit_ExMag(DetCurr,100000000)
            !call stop_all(this_routine, "Test complete")

            ! We only need to find out if determinant is connected to the
            ! reference (so no ex. level above 2 required, 
            ! truncated etc.)
            walkExcitLevel = FindBitExcitLevel (iLutRef(:,1), CurrentDets(:,j), &
                                                max_calc_ex_level, .true.)
            
            if(tRef_Not_HF) then
                walkExcitLevel_toHF = FindBitExcitLevel (iLutHF_true, CurrentDets(:,j), &
                                                max_calc_ex_level, .true.)
            else
                walkExcitLevel_toHF = walkExcitLevel
            endif

            ! if requested, average the sign over replicas if not coherent
            if(inum_runs > 1 .and. tWriteConflictLvls) call replica_coherence_check(&
                 CurrentDets(:,j), SignCurr, walkExcitLevel)

            if (tFillingStochRDMonFly) then
                ! Set the average sign and occupation iteration which were
                ! found in extract_bit_rep_avsign.
                call set_av_sgn_tot(j, AvSignCurr)
                call set_iter_occ_tot(j, IterRDMStartCurr)
                ! If this is an iteration where we print out the RDM energy,
                ! add in the diagonal contribution to the RDM for this
                ! determinant, for each rdm.
                if (tFill_RDM .and. (.not. tNoNewRDMContrib)) then
                    if (tOldRDMs) then
                        do irdm = 1, rdm_definitions%nrdms
                            call fill_rdm_diag_currdet_old(rdms(irdm), one_rdms_old(irdm), irdm, CurrentDets(:,j), &
                                                        DetCurr, j, walkExcitLevel_toHF, tCoreDet)
                        end do
                    end if

                    av_sign = get_av_sgn_tot(j)
                    iter_occ = get_iter_occ_tot(j)
                    call fill_rdm_diag_currdet(two_rdm_spawn, one_rdms, CurrentDets(:,j), DetCurr, &
                                                walkExcitLevel_toHF, av_sign, iter_occ, tCoreDet)
                endif
            endif

            ! This if-statement is only entered when using semi-stochastic and
            ! only if this determinant is in the core space.
            if (tCoreDet) then
                ! Store the index of this state, for use in annihilation later.
                indices_of_determ_states(determ_index) = j

                ! Add this amplitude to the deterministic vector.
                partial_determ_vecs(:,determ_index) = SignCurr

                determ_index = determ_index + 1

                ! The deterministic states are always kept in CurrentDets, even when
                ! the amplitude is zero. Hence we must check if the amplitude is zero,
                ! and if so, skip the state.
                if (IsUnoccDet(SignCurr)) cycle
            end if

            ! As the main list (which is storing a hash table) no longer needs
            ! to be contiguous, we need to skip sites that are empty.
            if(IsUnoccDet(SignCurr)) then
                !It has been removed from the hash table already
                !Add to the "freeslot" list
                iEndFreeSlot=iEndFreeSlot+1
                FreeSlot(iEndFreeSlot)=j
                cycle
            endif

            ! The current diagonal matrix element is stored persistently.
            HDiagCurr = det_diagH(j)

            if (tSeniorInitiators) then
                SpawnSign = get_all_spawn_pops(j)
                do run = 1, inum_runs
                    if(.not. is_run_unnocc(SignCurr, run) .and. .not. is_run_unnocc(SpawnSign, run)) then
#ifdef __CMPLX
                        !For complex walkers, we consider the sign changed when the argument of the complex 
                        !number changes more than pi/2.

                        CurrArg = DATAN2(SignCurr(max_part_type(run)), SignCurr(min_part_type(run)))
                        SpawnArg = DATAN2(SpawnSign(max_part_type(run)), SpawnSign(min_part_type(run)))
                        signChanged = mod(abs(CurrArg-SpawnArg), PI) > PI/2.0_dp
#else
                        signChanged = SpawnSign(min_part_type(run))*SignCurr(min_part_type(run)) < 0.0_dp
#endif
                    else
                        signChanged = .false.
                    end if
                    newlyOccupied = is_run_unnocc(SignCurr, run) .and. .not. is_run_unnocc(SpawnSign, run)
                    if ( signChanged .or. newlyOccupied) then
                        call reset_tau_int(j, run)
                        call reset_shift_int(j, run)
                        call set_spawn_pop(j, min_part_type(run), SignCurr(min_part_type(run)))
#ifdef __CMPLX
                        call set_spawn_pop(j, max_part_type(run), SignCurr(max_part_type(run)))
#endif
                    else
                        call update_tau_int(j, run, tau)
                        call update_shift_int(j, run, DiagSft(run)*tau)
                    end if
                end do
            end if

            if (tTruncInitiator) then
                call CalcParentFlag (j, parent_flags)
            end if

            !Debug output.
            IFDEBUGTHEN(FCIMCDebug,3)
                write(iout, "(A,I10,a)", advance='no') 'TW:', j, '['
                do part_type = 1, lenof_sign
                    write(iout, "(f10.5)", advance='no') SignCurr(part_type)
                end do
                write(iout, '(a,i7)', advance='no') '] ', FlagsCurr
                call WriteBitDet(iout,CurrentDets(:,j),.true.)
                call neci_flush(iout) 
            ENDIFDEBUG

!            call test_sym_excit3 (DetCurr, 1000000, pDoubles, 3)

            if(walkExcitLevel_toHF.eq.0) HFInd = j
            
            IFDEBUGTHEN(FCIMCDebug,1)
                if(j.gt.1) then
                    if(DetBitEQ(CurrentDets(:,j-1),CurrentDets(:,j),NIfDBO)) then
                        call stop_all(this_routine,"Shouldn't have the same determinant twice")
                    endif
                endif
            ENDIFDEBUG


            ! double occupancy measurement: quick and dirty for now
            if (t_calc_double_occ) then
                inst_double_occ = inst_double_occ + &
                    get_double_occupancy(CurrentDets(:,j), SignCurr)

                if (t_spin_measurements) then
                    call measure_double_occ_and_spin_diff(CurrentDets(:,j), &
                        DetCurr, SignCurr)
                end if
            end if

            call SumEContrib (DetCurr, WalkExcitLevel,SignCurr, CurrentDets(:,j), HDiagCurr, 1.0_dp, tPairedReplicas, j)

            ! If we're on the Hartree-Fock, and all singles and doubles are in
            ! the core space, then there will be no stochastic spawning from
            ! this determinant, so we can the rest of this loop.
            if (tSemiStochastic .and. ss_space_in%tDoubles .and. walkExcitLevel_toHF == 0 .and. tDetermHFSpawning) cycle

            ! Loop over the 'type' of particle. 
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
            !                 --> OR double run
            !                 --> part_type == 1, 2; population sets 1 and 2, both real

            ! alis additional idea to skip the number of attempted excitations
            ! for noninititators in the back-spawning approach
            ! remove that for now
            do part_type = 1, lenof_sign

                run = part_type_to_run(part_type)
                TempSpawnedPartsInd = 0

                ! Loop over all the particles of a given type on the 
                ! determinant. CurrentSign gives number of walkers. Multiply 
                ! up by AvMCExcits if attempting multiple excitations from 
                ! each walker (default 1.0_dp).
!                 call decide_num_to_spawn(SignCurr(part_type), AvMCExcits, WalkersToSpawn)
                call decide_num_to_spawn(SignCurr(part_type), HDiagCurr, AvMCExcits, WalkersToSpawn)

                do p = 1, WalkersToSpawn

                    ! Zero the bit representation, to ensure no extraneous
                    ! data gets through.
                    ilutnJ = 0_n_int
                    child = 0.0_dp

                    ! for the 3-body excitations i really do not want to change 
                    ! all the interfaces to the other excitation generators, 
                    ! which all just assume ex(2,2) as size.. so use a 
                    ! if here.. 
                    if (t_3_body_excits) then 
                        if (t_uniform_excits) then 
                            call gen_excit_uniform_k_space_hub_transcorr(DetCurr, CurrentDets(:,j), &
                                nJ, ilutnJ, exFlag, ic, ex, tParity, prob, & 
                                HElGen, fcimc_excit_gen_store, part_type) 

                        else if (t_mixed_excits) then 
                            call gen_excit_mixed_k_space_hub_transcorr(DetCurr, CurrentDets(:,j), &
                                nJ, ilutnJ, exFlag, ic, ex, tParity, prob, & 
                                HElGen, fcimc_excit_gen_store, part_type) 

                        else
                            call gen_excit_k_space_hub_transcorr(DetCurr, CurrentDets(:,j), &
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
                    if(tSkipRef(run) .and. all(nJ==projEdet(:,run))) then
                        !Set nJ to null
                        nJ(1) = 0
                    end if

                    ! If a valid excitation, see if we should spawn children.
                    if (.not. IsNullDet(nJ)) then

                        if (tSemiStochastic) then
                            call encode_child (CurrentDets(:,j), iLutnJ, ic, ex)

                            ! Temporary fix: FindExcitBitDet copies the flags of the parent onto the
                            ! child, which causes semi-stochastic simulations to crash. Should it copy
                            ! these flags? There are comments questioning this in create_particle, too.
                            iLutnJ(nOffFlag) = 0_n_int
                            
                            ! If the parent state in the core space.
                            if (test_flag(CurrentDets(:,j), flag_deterministic)) then
                                ! Is the spawned state in the core space?
                                tInDetermSpace = is_core_state(iLutnJ, nJ)
                                ! If spawning is from and to the core space, cancel it.
                                if (tInDetermSpace) cycle
                                ! Set the flag to specify that the spawning is occuring
                                ! from the core space.
                                call set_flag(iLutnJ, flag_determ_parent)
                            end if

                        end if

                        child = attempt_create (DetCurr, &
                                            CurrentDets(:,j), SignCurr, &
                                            nJ,iLutnJ, Prob, HElGen, IC, ex, &
                                            tParity, walkExcitLevel,part_type, &
                                            AvSignCurr,RDMBiasFacCurr)     
                                            ! Note these last two, AvSignCurr and 
                                            ! RDMBiasFacCurr are not used unless we're 
                                            ! doing an RDM calculation.
                    else
                                            
                        ! and rescale in the back-spawning algorithm.
                        ! this should always be a factor of 1 in the other 
                        ! cases so it is safe to rescale all i guess
                        ! (remove this option for now!)
                        ! [W.D. 20.3.2018:]
                        ! i dont quite understand this.. ask ghaldoon
                        ! because above this is set to 0.. and not 
                        ! changed ... 
                        child = child 
!                         child = 0.0_dp
                    endif

                    IFDEBUG(FCIMCDebug, 3) then
                        write(iout, '(a)', advance='no') 'SP: ['
                        do y = 1, lenof_sign
                            write(iout, '(f12.5)', advance='no') &
                                child(y)
                        end do
                        write(iout, '("] ")', advance='no')
                        call write_det(6, nJ, .true.)
                        call neci_flush(iout) 
                    endif

                    ! Children have been chosen to be spawned.
                    if (any(child /= 0)) then

                        ! Encode child if not done already.
                        if(.not. (tSemiStochastic)) call encode_child (CurrentDets(:,j), iLutnJ, ic, ex)
                        ! FindExcitBitDet copies the parent flags so that unwanted flags must be unset.
                        ! Should it really do this?
                        if (tTrialWavefunction) then
                            call clr_flag(iLutnJ, flag_trial)
                            call clr_flag(iLutnJ, flag_connected)
                        end if

                        call new_child_stats (iter_data, CurrentDets(:,j), &
                                              nJ, iLutnJ, ic, walkExcitLevel, &
                                              child, parent_flags, part_type)

                        if (use_spawn_hash_table) then
                            call create_particle_with_hash_table (nJ, ilutnJ, child, part_type, &
                                                                   CurrentDets(:,j), iter_data)
                        else
                            call create_particle (nJ, iLutnJ, child, part_type, & 
                                                  CurrentDets(:,j), SignCurr, p, &
                                                  RDMBiasFacCurr, WalkersToSpawn)
                        end if

                    endif ! (child /= 0), Child created.

                enddo ! Cycling over mulitple particles on same determinant.

            enddo   ! Cycling over 'type' of particle on a given determinant.

                ! If we are performing a semi-stochastic simulation and this state
                ! is in the deterministic space, then the death step is performed
                ! deterministically later. Otherwise, perform the death step now.
                if (.not. tCoreDet) call walker_death (iter_data, DetCurr, CurrentDets(:,j), &
                                                       HDiagCurr, SignCurr, j, WalkExcitLevel)

        enddo ! Loop over determinants.
        IFDEBUGTHEN(FCIMCDebug,2) 
            write(iout,*) 'Finished loop over determinants'
            write(iout,*) "Holes in list: ", iEndFreeSlot
        ENDIFDEBUG
        if (tSemiStochastic) then
            ! For semi-stochastic calculations only: Gather together the parts
            ! of the deterministic vector stored on each processor, and then
            ! perform the multiplication of the exact projector on this vector.
            call determ_projection()

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
                    if (tOldRDMs) call fill_RDM_offdiag_deterministic_old(rdms, one_rdms_old)
                    call fill_RDM_offdiag_deterministic(rdm_definitions, two_rdm_spawn, one_rdms)
                end if
            end if
        end if

        ! With this algorithm, the determinants do not move, and therefore
        ! TotWalkersNew is simply equal to TotWalkers
        TotWalkersNew = int(TotWalkers,sizeof_int)

        ! Update the statistics for the end of an iteration.
        ! Why is this done here - before annihilation!
        call end_iter_stats(TotWalkersNew)
        
        ! Print bloom/memory warnings
        call end_iteration_print_warn (totWalkersNew)
        call halt_timer (walker_time)

        ! For the direct annihilation algorithm. The newly spawned 
        ! walkers should be in a seperate array (SpawnedParts) and the other 
        ! list should be ordered.
        call set_timer (annihil_time, 30)
        !HolesInList is returned from direct annihilation with the number of unoccupied determinants in the list
        !They have already been removed from the hash table though.

        call DirectAnnihilation (totWalkersNew, iter_data, .false.) !.false. for not single processor

        ! This indicates the number of determinants in the list + the number
        ! of holes that have been introduced due to annihilation.
        TotWalkers = TotWalkersNew

        ! The superinitiators are now the same as they will be at the beginning of
        ! the next cycle (this flag is reset if they change)
        tReferenceChanged = .false.

        CALL halt_timer(Annihil_Time)
        IFDEBUG(FCIMCDebug,2) WRITE(iout,*) "Finished Annihilation step"
        

        ! If we are orthogonalising the replica wavefunctions, to generate
        ! excited states, then do that here.
        if (tOrthogonaliseReplicas .and. iter > orthogonalise_iter) then
            call orthogonalise_replicas(iter_data)
        else if (tPrintReplicaOverlaps .and. inum_runs > 1) then
            call calc_replica_overlaps()
        end if

        if (tFillingStochRDMonFly) then
            if (tOldRDMs) call fill_rdm_diag_wrapper_old(rdms, one_rdms_old, CurrentDets, int(TotWalkers, sizeof_int))
            call fill_rdm_diag_wrapper(rdm_definitions, two_rdm_spawn, one_rdms, CurrentDets, int(TotWalkers, sizeof_int))
        end if

        if(tTrialWavefunction .and. tTrialShift)then
            call fix_trial_overlap(iter_data)
        end if

        call update_iter_data(iter_data)

        ! This routine will take the CurrentDets and search the array to find all single and double 
        ! connections - adding them into the RDM's. 
        ! This explicit way of doing this is very expensive, but o.k for very small systems.
        if (tFillingExplicRDMonFly) then
            if (tHistSpawn) THEN
                call Fill_Hist_ExplicitRDM_this_Iter(TotWalkers)
            else
                call Fill_ExplicitRDM_this_Iter(TotWalkers)
            end if
        end if

        if (tFillingStochRDMonFly .or. tFillingExplicRDMonFly) then
            ! Fill the receiving RDM list from the beginning.
            two_rdm_recv%nelements = 0
            call communicate_rdm_spawn_t(two_rdm_spawn, two_rdm_recv)
            call add_rdm_1_to_rdm_2(two_rdm_recv, two_rdm_main)

!             if (t_calc_double_occ) then 
!                 call calc_double_occ_from_rdm(two_rdm_main, rdm_estimates%norm, &
!                     inst_rdm_occ)
!             end if

        end if
    end subroutine PerformFCIMCycPar

    subroutine test_routine()

        integer(n_int) :: list_1(0:NIfTot, 10)
        integer(n_int) :: list_2(0:NIfTot, 12)
        integer(n_int) :: list_out(0:NIfTot, 11)

        integer :: i
        integer :: ndets_out
        integer(n_int) :: dets_1(6), dets_2(6)
        real(dp) :: signs_1(6), signs_2(6)
        real(dp) :: real_sign(lenof_sign)

        dets_1 =  (/28,   47,  1,   23,   32,  57/)
        signs_1 = (/-1.0_dp, 2.0_dp, 1.3_dp, 4.0_dp, 1.0_dp, 7.0_dp/)

        dets_2 =  (/47,  23,  32,  38,  57,  63/)
        signs_2 = (/7.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 1.2_dp, 2.4_dp/)

        do i = 1, 6
            list_1(0,i) = dets_1(i) 
            real_sign = signs_1(i)
            call encode_sign(list_1(:,i), real_sign)
            list_1(3,i) = 0
            if (i == 1 .or. i == 2) call set_flag(list_1(:,i), flag_deterministic)
        end do

        do i = 1, 6
            list_2(0,i) = dets_2(i) 
            real_sign = signs_2(i)
            call encode_sign(list_2(:,i), real_sign)
            list_2(3,i) = 0
            if (i == 1) call set_flag(list_2(:,i), flag_deterministic)
        end do

        write(6,*) "List 1:"
        do i = 1, 6
            call extract_sign(list_1(:,i), real_sign)
            write(6,*) i, list_1(0,i), real_sign, list_1(3,i)
        end do

        write(6,*) "List 2:"
        do i = 1, 6
            call extract_sign(list_2(:,i), real_sign)
            write(6,*) i, list_2(0,i), real_sign, list_2(3,i)
        end do

        call add_ilut_lists(6, 6, .false., list_1, list_2, list_out, ndets_out)

        write(6,*) "Summed list:"
        do i = 1, ndets_out
            call extract_sign(list_out(:,i), real_sign)
            write(6,*) i, list_out(0,i), real_sign, list_out(3,i)
        end do

    end subroutine test_routine

END MODULE FciMCParMod

