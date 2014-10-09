#include "macros.h"
!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

!   The module now has the same structure with and without PARALLEL being defined.
!   Some routines require MPI and are enclosed in the #ifdef PARALLEL section.  These
!   should have dummy replacements in the #else of this if required.
!   At the end are functions which do not require parallel directives, and are accessible
!   for both parallel and non-parallel.
MODULE FciMCParMod
    use SystemData, only: nel, Brr, nBasis, nBasisMax, LMS, tHPHF, tHub, &
                          tReal, tRotatedOrbs, tFindCINatOrbs, tFixLz, &
                          LzTot, tUEG, tLatticeGens, tCSF, G1, Arr, &
                          tNoBrillouin, tKPntSym, tPickVirtUniform, &
                          tMolpro, csf_trunc_level, tMolproMimic, &
                          tTruncateCSF, tRef_Not_HF, &
                          MolproID, tGenHelWeighted, &
                          tGen_4ind_weighted, tGen_4ind_reverse
    use bit_rep_data, only: extract_sign, flag_trial, flag_connected, tUseFlags
    use bit_reps, only: NIfD, NIfTot, NIfDBO, NOffY, decode_bit_det, &
                        encode_bit_rep, encode_det, extract_bit_rep, &
                        test_flag, set_flag, extract_flags, &
                        flag_is_initiator, clear_all_flags, &
                        extract_sign, nOffSgn, flag_make_initiator, &
                        flag_parent_initiator, encode_sign, &
                        clr_flag, flag_trial, flag_connected, nOffFlag, &
                        flag_deterministic, flag_determ_parent, clr_flag, &
                        extract_part_sign, encode_part_sign, encode_first_iter
    use CalcData, only: InitWalkers, NMCyc, DiagSft, Tau, SftDamp, StepsSft, &
                        OccCASorbs, VirtCASorbs, NEquilSteps,&
                        tReadPops, iFullSpaceIter, MaxNoAtHF,&
                        tStartSinglePart, tCCMC, &
                        HFPopThresh, tTruncCAS, AvMCExcits, &
                        tTruncInitiator, &
                        NShiftEquilSteps, tWalkContGrow, &
                        tAddToInitiator, InitiatorWalkNo, &
                        tReadPopsRestart, tCheckHighestPopOnce, &
                        iRestartWalkNum, tRestartHighPop, FracLargerDet, &
                        tChangeProjEDet, tCheckHighestPop, tSpawnSpatialInit,&
                        MemoryFacInit, tMaxBloom, tTruncNOpen, tFCIMC, &
                        trunc_nopen_max, RealSpawnCutoff, &
                        TargetGrowRate, TargetGrowRateWalk, tShiftonHFPop, &
                        iExitWalkers,MemoryFacPart, &
                        tAllRealCoeff, tRealCoeffByExcitLevel, &
                        RealCoeffExcitThresh, &
                        tRealSpawnCutoff, RealSpawnCutoff, tDetermProj, &
                        tJumpShift, tUseRealCoeffs, tSpatialOnlyHash, &
                        tSemiStochastic, tTrialWavefunction, &
                        InitiatorCutoffEnergy, InitiatorCutoffWalkNo, &
                        tLetInitialPopDie, tFTLM, tWritePopsNorm, pops_norm_unit, &
                        tSpecLanc, tExactSpec, tExactDiagAllSym, tKP_FCIQMC
    use spatial_initiator, only: add_initiator_list, rm_initiator_list
    use HPHFRandExcitMod, only: FindExcitBitDetSym, gen_hphf_excit
    use Determinants, only: FDet, get_helement, write_det, &
                            get_helement_det_only, lexicographic_store, &
                            get_lexicographic_dets, DefDet
    use DetCalcData, only: ICILevel, nDet, Det, FCIDetIndex, FCIDets, lab, &
                           nRow, Hamil, ReIndex
    use GenRandSymExcitNUMod, only: gen_rand_excit, GenRandSymExcitNU, &
                                    ScratchSize, TestGenRandSymExcitNU, &
                                    ScratchSize1, ScratchSize2, ScratchSize3,&
                                    init_excit_gen_store,clean_excit_gen_store
    use GenRandSymExcitCSF, only: gen_csf_excit
    use IntegralsData, only: tPartFreezeCore, NPartFrozen, &
                             NHolesFrozen, tPartFreezeVirt, NVirtPartFrozen, &
                             NElVirtFrozen
    use LoggingData, only: iWritePopsEvery, TPopsFile, iPopsPartEvery, tBinPops, &
                           iWriteHistEvery, tHistEnergies, FCIMCDebug, &
                           AllHistInitPops, &
                           OffDiagBinRange, OffDiagMax, AllHistInitPopsTag, &
                           tLogComplexPops, tPrintFCIMCPsi, tCalcFCIMCPsi, &
                           NHistEquilSteps, tPrintOrbOcc, &
                           HistInitPopsTag, OrbOccs, OrbOccsTag, &
                           tPrintPopsDefault, tHistInitPops, HistInitPops, &
                           tCalcInstantS2, instant_s2_multiplier, tMCOutput, &
                           tDiagWalkerSubspace,iDiagSubspaceIter, &
                           tRDMonFly, IterRDMonFly,RDMExcitLevel, RDMEnergyIter, &
                           tChangeVarsRDM, tExplicitAllRDM, &
                           tDiagWalkerSubspace, iDiagSubspaceIter, &
                           tCalcInstantS2Init, instant_s2_multiplier_init, &
                           tJustBlocking, iBlockEquilShift, iBlockEquilProjE, &
                           tDiagAllSpaceEver, tCalcVariationalEnergy, tCompareTrialAmps, &
                           compare_amps_period, tNoNewRDMContrib, &
                           tFCIMCStats2, tHistExcitToFrom, &
                           tSpawnGhostChild, GhostThresh, &
                           tWriteCoreEnd, write_end_core_size
    use hist, only: init_hist_spin_dist, clean_hist_spin_dist, &
                    hist_spin_dist, ilut_spindist, tHistSpinDist, &
                    write_clear_hist_spin_dist, hist_spin_dist_iter, &
                    tHistSpawn, AllHistogramEnergy, &
                    AllHistogram, HistogramEnergy, Histogram, AllInstHist, &
                    InstHist, HistMinInd, project_spins, calc_s_squared, &
                    project_spin_csfs, calc_s_squared_multi, &
                    calc_s_squared_star, init_hist_excit_tofrom, &
                    add_hist_excit_tofrom, write_zero_hist_excit_tofrom, &
                    clean_hist_excit_tofrom
    use hist_data, only: beforenormhist, HistMinInd2, BinRange, iNoBins
    USE SymData , only : nSymLabels, Sym_Psi
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE Parallel_neci
    USE FciMCData
    USE AnnihilationMod, only: DirectAnnihilation
    use PopsfileMod, only: ReadFromPopsfilePar, FindPopsfileVersion, &
                           WriteToPopsFileParOneArr, open_pops_head, &
                           readpopsheadv3, readpopsheadv4, CheckPopsParams, &
                           ReadFromPopsFile, InitFCIMC_pops
    use sort_mod
    use DetBitops, only: EncodeBitDet, DetBitEQ, DetBitLT, FindExcitBitDet, &
                         FindBitExcitLevel, countbits, TestClosedShellDet, &
                         FindSpatialBitExcitLevel, IsAllowedHPHF, count_open_orbs, &
                         ilut_gt, get_bit_excitmat
    use hash, only: DetermineDetNode, FindWalkerHash, remove_hash_table_entry
    use csf, only: get_csf_bit_yama, iscsf, csf_orbital_mask, get_csf_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_spawn_sign, hphf_off_diag_helement_spawn
    use util_mod
    use constants
    use soft_exit, only: ChangeVars 
    use RotateOrbsMod, only: RotateOrbs
    use NatOrbsMod, only: PrintOrbOccs
    use spin_project, only: tSpinProject, spin_proj_interval, &
                            spin_proj_gamma, get_spawn_helement_spin_proj, &
                            generate_excit_spin_proj, attempt_die_spin_proj, &
                            iter_data_spin_proj, test_spin_proj, &
                            spin_proj_shift, spin_proj_iter_count, &
                            init_yama_store, clean_yama_store, &
                            disable_spin_proj_varyshift
    use symrandexcit3, only: gen_rand_excit3, test_sym_excit3
    USE SymExcit3 , only : GenExcitations3
    use errors, only: error_analysis
    use nElRDMMod, only: FinaliseRDM,Fill_ExplicitRDM_this_Iter,calc_energy_from_rdm, &
                         fill_hist_explicitrdm_this_iter, tCalc_RDMEnergy, &
                         extract_bit_rep_avsign_norm, &
                         extract_bit_rep_avsign_no_rdm, &
                         zero_rdms, store_parent_with_spawned, &
                         fill_rdm_diag_currdet_norm, calc_rdmbiasfac, &
                         det_removed_fill_diag_rdm, fill_rdm_offdiag_deterministic
    use determ_proj, only: perform_determ_proj
    use semi_stoch_gen, only: init_semi_stochastic, enumerate_sing_doub_kpnt, &
                              write_most_pop_core_at_end
    use semi_stoch_procs, only: determ_projection, return_most_populated_states, &
                                end_semistoch, is_core_state, return_mp1_amp_and_mp2_energy, &
                                recalc_core_hamil_diag, check_determ_flag, &
                                average_determ_vector
    use trial_wf_gen, only: init_trial_wf, update_compare_trial_file, end_trial_wf
    use ftlm_neci, only: perform_ftlm
    use spectral_lanczos, only: perform_spectral_lanczos
    use exact_spectrum, only: get_exact_spectrum
    use exact_diag, only: perform_exact_diag_all_symmetry
    use gndts_mod, only: gndts
    use sort_mod
    use get_excit, only: make_double
    use sltcnd_mod, only: sltcnd_excit
    use excit_gens_int_weighted, only: gen_excit_hel_weighted, &
                                       gen_excit_4ind_weighted, &
                                       test_excit_gen_4ind, &
                                       gen_excit_4ind_reverse
    use procedure_pointers
    use fcimc_helper
    use tau_search, only: init_tau_search, log_spawn_magnitude, update_tau, &
                          log_death_magnitude
    use global_det_data, only: init_global_det_data, clean_global_det_data, &
                               global_determinant_data, set_det_diagH, &
                               det_diagH, get_av_sgn, set_av_sgn, set_iter_occ

    implicit none

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)
        use LoggingData, only: PopsfileTimer
        use nElRDMMod, only: InitRDM
        use sym_mod, only: getsym
        use SystemData, only: tUEG2, SymRestrict
#ifdef MOLPRO
        use outputResult
        integer :: nv,ityp(1)
#endif
        integer :: iroot,isymh
        real(dp) :: Weight, Energyxw,BestEnergy
        INTEGER :: error
        LOGICAL :: TIncrement,tWritePopsFound,tSoftExitFound,tSingBiasChange,tPrintWarn
        REAL(sp) :: s_start,s_end,tstart(2),tend(2),totaltime
        real(dp) :: TotalTime8
        INTEGER(int64) :: MaxWalkers,MinWalkers
        real(dp) :: AllTotWalkers,MeanWalkers,Inpair(2),Outpair(2)
        integer, dimension(lenof_sign) :: tmp_sgn
        integer :: tmp_int(lenof_sign), i, istart
        real(dp) :: grow_rate,EnergyDiff,Norm_2RDM
        TYPE(BasisFn) RefSym
        real(dp) :: mean_ProjE_re,mean_ProjE_im,mean_Shift
        real(dp) :: ProjE_Err_re,ProjE_Err_im,Shift_Err
        logical :: tNoProjEValue,tNoShiftValue
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
        character(6), parameter :: excit_descriptor(0:2) = &
                                        (/"IC0   ", "single", "double"/)

        if(tJustBlocking) then
            !Just reblock the current data, and do not perform an fcimc calculation
            write(6,"(A)") "Skipping FCIQMC calculation and simply reblocking previous output"
            call Standalone_Errors()
            return
        endif

        TDebug=.false.  !Set debugging flag
                    
!OpenMPI does not currently support MPI_Comm_set_errhandler - a bug in its F90 interface code.
!Ask Nick McLaren if we need to change the err handler - he has a fix/bypass.
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,error)
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL,error)
        
        ! This is set here not in SetupParameters, as otherwise it would be
        ! wiped just when we need it!
        tPopsAlreadyRead = .false.

        call SetupParameters()
        call InitFCIMCCalcPar()
        call init_fcimc_fn_pointers() 

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
            end if
            call WriteFCIMCStats()
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

        SumSigns = 0.0_dp
        SumSpawns = 0.0_dp

        ! In we go - start the timer for scaling curve!
        start_time = neci_etime(tstart)

        do while (Iter <= NMCyc .or. NMCyc == -1)
!Main iteration loop...
            IFDEBUG(FCIMCDebug, 2) write(iout,*) 'Iter', iter

            if(iProcIndex.eq.root) s_start=neci_etime(tstart)
            
            if(tRDMonFly .and. (.not. tFillingExplicRDMonFly) &
                & .and. (.not.tFillingStochRDMonFly)) call check_start_rdm()

            if(tRDMonFly .and. tSpawnGhostChild) &
                    call stop_all("FciMCPar", "tSpawnGhostChild is not yet working correctly with &
                    & the RDMs.  I need to introduce a ghost flag for the ghost progeny so that &
                    & we know which population the spawning event came from once we get to annihilation. &
                    & See approx line 449 in Annihilation.F90 where we assign &
                    & Spawned_Parents NIfDBO+2,Parent_Array_Ind to say which pop the spawning event &
                    & came from")


            if (tCCMC) then
                if (tUseRealCoeffs) &
                    call stop_all("FciMCPar", "Continuous spawning not yet set up to work &
                    & with CCMC.  Attention will be required in a number of routines, &
                    & including AttemptCreatePar and PerformCCMCCycParInt")

                if(inum_runs.eq.2) &
                    call stop_all("FciMCPar", "Double runs not set up to work with CCMC")

                CALL PerformCCMCCycPar()
            else
                if (.not. (tSpinProject .and. spin_proj_interval == -1)) &
                    call PerformFciMCycPar(iter_data_fciqmc)
            endif

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

            if (mod(Iter, StepsSft) == 0) then

                ! Has there been a particle bloom this update cycle? Loop
                ! through the spawned particle types, and output details.
                ! If outputting only the biggest blooms to date, keep
                ! track of that.
                if (iProcIndex == Root) then
                    istart = 1
                   ! if (tSpinProjDets) istart = 0
                    do i = istart, 2
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
                if (tCCMC) then
                    call calculate_new_shift_wrapper (iter_data_ccmc, &
                                                      TotParts)
                else
                    call calculate_new_shift_wrapper (iter_data_fciqmc, &
                                                      TotParts)
                endif

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

                if(iProcIndex.eq.root) TotalTime8=real(s_end,dp)
                call MPIBCast(TotalTime8)    !TotalTime is local - broadcast to all procs

!This routine will check for a CHANGEVARS file and change the parameters of the calculation accordingly.
                CALL ChangeVars(tSingBiasChange,tSoftExitFound,tWritePopsFound)
                IF(tSoftExitFound) THEN
                    !Now changed such that we do one more block of iterations and then exit, to allow for proper
                    !inclusion of all the RDM elements
                    NMCyc=Iter+StepsSft  
                    
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
                    write(iout,"(A,F8.2,A)") "Time limit reached for simulation of: ",MaxTimeExit/60.0," minutes - exiting..."
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
!This will write out the POPSFILE if wanted
                CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
            ENDIF
!            IF(TAutoCorr) CALL WriteHistogrammedDets()

            IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0).and.(.not.tRDMonFly)) THEN
                CALL WriteHistogram()
            ENDIF
            IF(tRDMonFly.and.(.not.tSinglePartPhase(1)).and. &
                        (.not.(tSinglePartPhase(inum_runs)))) THEN
                ! If we wish to calculate the energy, have started accumulating the RDMs, 
                ! and this is an iteration where the energy should be calculated, do so.
                if(tCalc_RDMEnergy .and. ((Iter - maxval(VaryShiftIter)).gt.IterRDMonFly) &
                    .and. (mod((Iter+PreviousCycles - IterRDMStart)+1,RDMEnergyIter).eq.0) ) &
                        CALL Calc_Energy_from_RDM(Norm_2RDM)  
            ENDIF
            if(tChangeVarsRDM) then
                ! Decided during the CHANGEVARS that the RDMs should be calculated.
                call InitRDM() 
                tRDMonFly = .true.
                tChangeVarsRDM = .false.
            endif

            if(tDiagWalkerSubspace.and.(mod(Iter,iDiagSubspaceIter).eq.0)) then
                !Diagonalise a subspace consisting of the occupied determinants
                call DiagWalkerSubspace()
            endif

            ! If requested and on a correct iteration, update the COMPARETRIAL file.
            if (tCompareTrialAmps .and. mod(Iter, compare_amps_period) == 0) &
                call update_compare_trial_file(.false.)

            Iter=Iter+1

!End of MC cycle
        enddo

        ! We are at the end - get the stop-time. Output the timing details
        stop_time = neci_etime(tend)
        write(iout,*) '- - - - - - - - - - - - - - - - - - - - - - - -'
        write(iout,*) 'Total loop-time: ', stop_time - start_time
        write(iout,*) '- - - - - - - - - - - - - - - - - - - - - - - -'

        ! Reduce the iteration count fro the POPSFILE since it is incremented
        ! upon leaving the loop (If done naturally).
        IF(TIncrement) Iter=Iter-1
        IF(TPopsFile) THEN
            WRITE(6,*) "Totwalkers", TotWalkers
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

        ! If requested, write the most populated states in CurrentDets to a
        ! CORESPACE file, for use in future semi-stochastic calculations.
        if (tWriteCoreEnd) call write_most_pop_core_at_end(write_end_core_size)

        IF(tHistSpawn) CALL WriteHistogram()


        IF(tHistEnergies) CALL WriteHistogramEnergies()

        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
        ENDIF

        IF(tFillingStochRDMonFly.or.&
            tFillingExplicRDMonFly) CALL FinaliseRDM()
            !tFillingExplicRDMonFly.or.tHF_Ref_Explicit) CALL FinaliseRDM()

        call PrintHighPops()

        !Close open files.
        IF(iProcIndex.eq.Root) THEN
            CLOSE(fcimcstats_unit)
            if (inum_runs.eq.2) CLOSE(fcimcstats_unit2)
            IF(tTruncInitiator) CLOSE(initiatorstats_unit)
            IF(tLogComplexPops) CLOSE(complexstats_unit)
            if (tWritePopsNorm) close(pops_norm_unit)
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
       
        !Automatic error analysis
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
        Weight=(0.0_dp)
        if (tTrialWavefunction) then
            Energyxw = tot_trial_numerator(1)/tot_trial_denom(1) + Hii
        else
            Energyxw=(ProjectionE(1)+Hii)
        end if
        
        iroot=1
        CALL GetSym(ProjEDet,NEl,G1,NBasisMax,RefSym)
        isymh=int(RefSym%Sym%S,sizeof_int)+1
        write (iout,10101) iroot,isymh
10101   format(//'RESULTS FOR STATE',i2,'.',i1/'====================='/)
        write (iout,'('' Current reference energy'',T52,F19.12)') Hii 
        if(tNoProjEValue) then
            write (iout,'('' Projected correlation energy'',T52,F19.12)') ProjectionE
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
            BestEnergy = ProjectionE(1)+Hii
            BestErr = 0.0_dp
        endif
        write(iout,"(A)")
        if(tNoProjEValue) then
            write(iout,"(A,F20.8)") " Total projected energy ",ProjectionE+Hii
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
 
        !Deallocate memory
        CALL DeallocFCIMCMemPar()

    END SUBROUTINE FciMCPar


    subroutine PerformFCIMCycPar(iter_data)
        
        ! Iteration specific data
        type(fcimc_iter_data), intent(inout) :: iter_data

        ! Now the local, iteration specific, variables
        integer :: VecSlot, j, p, error, proc_temp, i, HFPartInd,isym
        integer :: DetCurr(nel), nJ(nel), FlagsCurr, parent_flags
        real(dp), dimension(lenof_sign) :: SignCurr, child
        integer(kind=n_int) :: iLutnJ(0:niftot)
        integer :: IC, walkExcitLevel, walkExcitLevel_toHF, ex(2,2), TotWalkersNew, part_type
        integer(int64) :: tot_parts_tmp(lenof_sign)
        logical :: tParity, tSuccess, tCoreDet
        real(dp) :: prob, HDiagCurr, TempTotParts, Di_Sign_Temp
        real(dp) :: RDMBiasFacCurr
        real(dp), dimension(lenof_sign) :: AvSignCurr, IterRDMStartCurr
        HElement_t :: HDiagTemp,HElGen
        character(*), parameter :: this_routine = 'PerformFCIMCycPar' 
        HElement_t, dimension(inum_runs) :: delta
        integer :: proc, pos
        real(dp) :: r, sgn(lenof_sign), prob_extra_walker
        integer :: determ_index, gen_ind
        integer :: DetHash, FinalVal, clash, PartInd, k
        type(ll_node), pointer :: TempNode

        call set_timer(Walker_Time,30)

        MaxInitPopPos=0.0
        MaxInitPopNeg=0.0
        HighPopNeg=1
        HighPopPos=1
        FlagsCurr=0
        ! Synchronise processors
!        CALL MPIBarrier(error)

        ! Reset iteration variables
        VecSlot = 1    ! Next position to write into CurrentDets
        ! Next free position in newly spawned list.
        ValidSpawnedList = InitialSpawnedSlots
        if(tHashWalkerList) then
            FreeSlot(1:iEndFreeSlot)=0  !Does this cover enough?
            iStartFreeSlot=1
            iEndFreeSlot=0
        endif

        ! Index for counting deterministic states.
        determ_index = 1
        
        call rezero_iter_stats_each_iter (iter_data)

        ! The processor with the HF determinant on it will have to check 
        ! through each determinant until it's found. Once found, tHFFound is
        ! true, and it no longer needs to be checked.

        ! This is a bit of a hack based on the fact that we mean something 
        ! different by exFlag for CSFs than in normal determinential code.
        ! It would be nice to fix this properly
        if (tCSF) exFlag = 7

        IFDEBUGTHEN(FCIMCDebug,iout)
            if(tHashWalkerList) then
                write(iout,"(A)") "Hash Table: "
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
            end if
        ENDIFDEBUG

        IFDEBUG(FCIMCDebug,3) write(iout,"(A,I12)") "Walker list length: ",TotWalkers
        IFDEBUG(FCIMCDebug,3) write(iout,"(A)") "TW: Walker  Det"

        ! This block decides whether or not to calculate the contribution to the RDMs from 
        ! the diagonal elements (and explicit connections to the HF) for each occupied determinant. 
        ! For efficiency, this is only done on the final iteration, or one where the RDM energy is 
        ! being printed.
        tFill_RDM = .false.
        if(tFillingStochRDMonFly) then
            if(mod((Iter+PreviousCycles - IterRDMStart + 1),RDMEnergyIter).eq.0) then 
                ! RDM energy is being printed, calculate the diagonal elements for 
                ! the last RDMEnergyIter iterations.
                tFill_RDM = .true.
                IterLastRDMFill = RDMEnergyIter
            elseif(Iter.eq.NMCyc) then
                ! Last iteration, calculate the diagonal element for the iterations 
                ! since the last time they were included.
                tFill_RDM = .true.
                IterLastRDMFill = mod((Iter+PreviousCycles - IterRDMStart + 1),RDMEnergyIter)
            endif
        endif
        
        do j=1,int(TotWalkers,sizeof_int)
            ! N.B. j indicates the number of determinants, not the number
            !      of walkers.

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

            call extract_bit_rep_avsign (CurrentDets(:,j), j, &
                                        DetCurr, SignCurr, FlagsCurr, IterRDMStartCurr, &
                                        AvSignCurr, fcimc_excit_gen_store)

            ! We only need to find out if determinant is connected to the
            ! reference (so no ex. level above 2 required, 
            ! truncated etc.)
            walkExcitLevel = FindBitExcitLevel (iLutRef, CurrentDets(:,j), &
                                                max_calc_ex_level)
            
            if(tRef_Not_HF) then
                walkExcitLevel_toHF = FindBitExcitLevel (iLutHF_true, CurrentDets(:,j), &
                                                max_calc_ex_level)
            else
                walkExcitLevel_toHF = walkExcitLevel
            endif

            if (tFillingStochRDMonFly) then
                ! Set the average sign and occupation iteration which were
                ! found in extract_bit_rep_avsign.
                call set_av_sgn(j, AvSignCurr)
                call set_iter_occ(j, IterRDMStartCurr)
                ! If this is an iteration where print out the RDM energy,
                ! calculate the diagonal contribution to the RDM for this
                ! determinant.
                if(tFill_RDM .and. (.not. tNoNewRDMContrib)) then
                    call fill_rdm_diag_currdet(CurrentDets(:,j), DetCurr, j, &
                                                walkExcitLevel_toHF, tCoreDet)
                endif
            endif

            ! A general index whose value depends on whether the following option is used.
            if (tHashWalkerList) then
                gen_ind = j
            else
                gen_ind = VecSlot
            end if

            ! This if-statement is only entered when using semi-stochastic and
            ! only if this determinant is in the core space.
            if (tCoreDet) then
                ! Store the index of this state, for use in annihilation later.
                indices_of_determ_states(determ_index) = gen_ind

                ! Add this amplitude to the deterministic vector.
                partial_determ_vector(:,determ_index) = SignCurr

                determ_index = determ_index + 1

                ! The deterministic states are always kept in CurrentDets, even when
                ! the amplitude is zero. Hence we must check if the amplitude is zero,
                ! and if so, skip the state.
                if (IsUnoccDet(SignCurr)) then
                    CurrentDets(:,gen_ind) = CurrentDets(:,j)
                    if (tFillingStochRDMonFly) then
                        call set_av_sgn(gen_ind, AvSignCurr)
                        call set_iter_occ(gen_ind, IterRDMStartCurr)
                    endif
                    VecSlot = VecSlot + 1
                    cycle
                end if
            end if

            ! The current diagonal matrix element is stored persistently.
            HDiagCurr = det_diagH(j)

            if (tTruncInitiator) &
                call CalcParentFlag (j, VecSlot, parent_flags, HDiagCurr)

            if(tHashWalkerList) then
                !Test here as to whether this is a "hole" or not...
                !Unfortunately, the main list no longer needs to be contiguous
                if(IsUnoccDet(SignCurr)) then
                    !It has been removed from the hash table already
                    !Add to the "freeslot" list
                    iEndFreeSlot=iEndFreeSlot+1
                    FreeSlot(iEndFreeSlot)=j
                    cycle
                endif
            endif

            !Debug output.
            IFDEBUGTHEN(FCIMCDebug,3)
                if(lenof_sign.eq.2) then
                    write(iout,"(A,I10,2f12.5,I5)",advance='no') "TW:", j,SignCurr,FlagsCurr
                else
                    write(iout,"(A,I10,f12.5,I5)",advance='no') "TW:", j,SignCurr,FlagsCurr
                endif
                call WriteBitDet(iout,CurrentDets(:,j),.true.)
                call neci_flush(iout) 
            ENDIFDEBUG

!            call test_sym_excit3 (DetCurr, 1000000, pDoubles, 3)

            if(walkExcitLevel_toHF.eq.0) HFInd = VecSlot
            
            IFDEBUGTHEN(FCIMCDebug,1)
                if(j.gt.1) then
                    if(DetBitEQ(CurrentDets(:,j-1),CurrentDets(:,j),NIfDBO)) then
                        call stop_all(this_routine,"Shouldn't have the same determinant twice")
                    endif
                endif
            ENDIFDEBUG

            ! Sum in any energy contribution from the determinant, including 
            ! other parameters, such as excitlevel info.
            ! This is where the projected energy is calculated.
            call SumEContrib (DetCurr, WalkExcitLevel,SignCurr, CurrentDets(:,j), HDiagCurr, 1.0_dp, j)

            ! Loop over the 'type' of particle. 
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
            !                 --> OR double run
            !                 --> part_type == 1, 2; population sets 1 and 2, both real
            do part_type=1,lenof_sign
            
                TempSpawnedPartsInd = 0

                ! Loop over all the particles of a given type on the 
                ! determinant. CurrentSign gives number of walkers. Multiply 
                ! up by AvMCExcits if attempting multiple excitations from 
                ! each walker (default 1.0_dp).
                call decide_num_to_spawn(SignCurr(part_type), AvMCExcits, WalkersToSpawn)

                do p = 1, WalkersToSpawn
                    ! Zero the bit representation, to ensure no extraneous
                    ! data gets through.
                    ilutnJ = 0_n_int

                    ! Generate a (random) excitation
                    call generate_excitation (DetCurr, CurrentDets(:,j), nJ, &
                                        ilutnJ, exFlag, IC, ex, tParity, prob, &
                                        HElGen, fcimc_excit_gen_store)
                    
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
                        child = 0.0_dp
                    endif

                    IFDEBUG(FCIMCDebug, 3) then
#ifdef __CMPLX
                        write(iout,"(a,2f12.5)",advance='no') &
#else
                        write(iout,"(a,f12.5)",advance='no') &
#endif
                            "SP:", child
                        call write_det(6, nJ, .true.)
                        call neci_flush(iout) 
                    endif

                    ! Children have been chosen to be spawned.
                    if (any(child /= 0) .or. tGhostChild ) then

                        !Encode child if not done already
                        if(.not. (tSemiStochastic)) call encode_child (CurrentDets(:,j), iLutnJ, ic, ex)
                        ! FindExcitBitDet copies the parent flags so that unwanted flags must be unset.
                        ! Should it really do this?
                        if (tTrialWavefunction) then
                            call clr_flag(iLutnJ, flag_trial)
                            call clr_flag(iLutnJ, flag_connected)
                        end if

                        call new_child_stats (iter_data, CurrentDets(:,j), &
                                              nJ, iLutnJ, ic, walkExcitLevel,&
                                              child, parent_flags, part_type)
                        call create_particle (nJ, iLutnJ, child, &
                                              parent_flags, part_type,& 
                                              CurrentDets(:,j),SignCurr,p,&
                                              RDMBiasFacCurr, WalkersToSpawn)
                                              ! RDMBiasFacCurr is only used if we're 
                                              ! doing an RDM calculation.

                    endif ! (child /= 0). Child created

                enddo ! Cycling over mulitple particles on same determinant.

            enddo   ! Cycling over 'type' of particle on a given determinant.
            
            ! DEBUG
            ! if (VecSlot > j) call stop_all (this_routine, 'vecslot > j')

            if (tSemiStochastic) then
                ! If we are performing a semi-stochastic simulation and this state is in the
                ! deterministic space, then the death step is performed deterministically later.
                if (.not. tCoreDet) then
                    call walker_death (iter_data, DetCurr, &
                                       CurrentDets(:,j), HDiagCurr, SignCurr, &
                                       AvSignCurr, IterRDMStartCurr, VecSlot, j, WalkExcitLevel)
                else
                    CurrentDets(:,gen_ind) = CurrentDets(:,j)
                    if (tFillingStochRDMonFly) then
                        call set_av_sgn(gen_ind, AvSignCurr)
                        call set_iter_occ(gen_ind, IterRDMStartCurr)
                    endif
                    ! We never overwrite the deterministic states, so move the
                    ! next slot in CurrentDets to the next state.
                    VecSlot = VecSlot + 1
                end if
            else
                call walker_death (iter_data, DetCurr, &
                                   CurrentDets(:,j), HDiagCurr, SignCurr, &
                                   AvSignCurr, IterRDMStartCurr, VecSlot, j, WalkExcitLevel)
            end if

        enddo ! Loop over determinants.
        IFDEBUGTHEN(FCIMCDebug,2) 
            write(iout,*) 'Finished loop over determinants'
            if(tHashWalkerList) then
                write(iout,*) "Holes in list: ",iEndFreeSlot
            endif
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
                if(tFill_RDM) call fill_RDM_offdiag_deterministic()
            end if
        end if

        if(tHashWalkerList) then
            ! With this algorithm, the determinants do not move, and therefore
            ! TotWalkersNew is simply equal to TotWalkers
            TotWalkersNew=int(TotWalkers,sizeof_int)
        else
            ! Since VecSlot holds the next vacant slot in the array, TotWalkers
            ! should be one less than this. TotWalkersNew is now the number of particles
            ! in the main array, before annihilation
            TotWalkersNew = VecSlot - 1
        endif

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
        call DirectAnnihilation (totWalkersNew, iter_data,.false.) !.false. for not single processor

        TotWalkers=TotWalkersNew    !with tHashWalkerList this indicates the number of determinants 
                                    !in list + holes from annihilation
        CALL halt_timer(Annihil_Time)
        IFDEBUG(FCIMCDebug,2) WRITE(iout,*) "Finished Annihilation step"

        call update_iter_data(iter_data)

        ! This routine will take the CurrentDets and search the array to find all single and double 
        ! connections - adding them into the RDM's. 
        ! This explicit way of doing this is very expensive, but o.k for very small systems.
        IF(tFillingExplicRDMonFly) THEN
            IF(tHistSpawn) THEN
                CALL Fill_Hist_ExplicitRDM_this_Iter(TotWalkers)
            ELSE
                CALL Fill_ExplicitRDM_this_Iter(TotWalkers)
            ENDIF
        ENDIF

    end subroutine

    subroutine update_iter_data(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data

        iter_data%update_growth = iter_data%update_growth + iter_data%nborn &
                                - iter_data%ndied - iter_data%nannihil &
                                - iter_data%naborted - iter_data%nremoved
        iter_data%update_iters = iter_data%update_iters + 1

    end subroutine update_iter_data

    subroutine end_iteration_print_warn (totWalkersNew)
        
        ! Worker function for PerformFciMCycPar. Prints warnings about 
        ! particle blooms and memory usage.
        integer, intent(in) :: totWalkersNew
        integer :: i
        real(dp) :: rat

        ! Too many particles?
        rat = real(TotWalkersNew,dp) / real(MaxWalkersPart,dp)
        if (rat > 0.95_dp) then
            if(tMolpro) then
                write (iout, '(a)') '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory. &
                                 &Errors imminent. Increase MEMORYFACWALKERS, or reduce rate of growth.'
            else
                write (iout, '(a)') '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory. &
                                 &Errors imminent. Increase MEMORYFACPART, or reduce rate of growth.'
            endif
            call neci_flush(iout)
        end if

        ! Are ony of the sublists near the end of their alloted space?
        if (nNodes > 1) then
            do i = 0, nNodes-1
                rat = real(ValidSpawnedList(i) - InitialSpawnedSlots(i),dp) /&
                             real(InitialSpawnedSlots(1), dp)
                if (rat > 0.95_dp) then
                    if(tMolpro) then
                        write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory.&
                                         &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                    else
                        write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory.&
                                         &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                    endif
                    call neci_flush(iout)
                endif
            enddo
        else
            rat = real(ValidSpawnedList(0), dp) / real(MaxSpawned, dp)
            if (rat > 0.95_dp) then
                if(tMolpro) then
                    write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory.&
                                     &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                else
                    write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory.&
                                     &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                endif
                call neci_flush(iout)
            endif
        endif

        ! Are we near the end of the spatial initiator list
        if (tSpawnSpatialInit) then
            rat = real(no_spatial_init_dets,dp) / real(max_inits,dp)
            if (rat > 0.95_dp) then
                write(iout, '(a)') '*WARNING* - Number of spatial initiators has&
                                & reached over 95% f max_inits.'
                call neci_flush(iout)
            endif
        endif

    end subroutine end_iteration_print_warn 

    subroutine decide_num_to_spawn(parent_pop, av_spawns_per_walker, nspawn)

        real(dp), intent(in) :: parent_pop
        real(dp), intent(in) :: av_spawns_per_walker
        integer, intent(out) :: nspawn
        real(dp) :: prob_extra_walker, r

        nspawn = abs(int(parent_pop*av_spawns_per_walker))
        if (abs(parent_pop*av_spawns_per_walker) - real(nspawn,dp) > 0) then
            prob_extra_walker = abs(parent_pop*av_spawns_per_walker) - real(nspawn,dp)
            r = genrand_real2_dSFMT()
            if (prob_extra_walker > r) nspawn = nspawn + 1
        end if

    end subroutine decide_num_to_spawn



    SUBROUTINE CheckOrdering(DetArray,SignArray,NoDets,tCheckSignCoher)
        INTEGER :: NoDets,i,j
        INTEGER(KIND=n_int) :: DetArray(0:NIfTot,1:NoDets)
        INTEGER :: SignArray(1:NoDets),Comp
        LOGICAL :: tCheckSignCoher

        IF(NoDets.gt.0) THEN
            IF(SignArray(1).eq.0) THEN
                WRITE(iout,*) "Iter: ",Iter,1
                CALL Stop_All("CheckOrdering","Array has annihilated particles in it...")
            ENDIF
        ENDIF
        do i=2,NoDets
            IF(SignArray(i).eq.0) THEN
                WRITE(iout,*) "Iter: ",Iter,i
                CALL Stop_All("CheckOrdering","Array has annihilated particles in it...")
            ENDIF
            Comp=DetBitLT(DetArray(:,i-1),DetArray(:,i),NIfDBO)
            IF(Comp.eq.-1) THEN
!Array is in reverse numerical order for these particles
                do j=max(i-5,1),min(i+5,NoDets)
                    WRITE(iout,*) Iter,j,DetArray(:,j),SignArray(j)
                enddo
                CALL Stop_All("CheckOrdering","Array not ordered correctly")
            ELSEIF(Comp.eq.0) THEN
!Dets are the same - see if we want to check sign-coherence
                IF(tCheckSignCoher) THEN
!!This bit checks that there is only one copy of the determinants in the list
                    do j=max(i-5,1),min(i+5,NoDets)
                        WRITE(iout,*) Iter,j,DetArray(:,j),SignArray(j)
                    enddo
                    CALL Stop_All("CheckOrdering","Determinant same as previous one...")
                ENDIF
                IF(tCheckSignCoher.and.(SignArray(i-1).ne.SignArray(i))) THEN
!This checks that any multple copies in the list are sign-coherent...
                    do j=i-5,i+5
                        WRITE(iout,*) Iter,j,DetArray(:,j),SignArray(j)
                    enddo
                    CALL Stop_All("CheckOrdering","Array not sign-coherent")
                ENDIF
            ENDIF
        enddo

    END SUBROUTINE CheckOrdering
        


    subroutine walker_death (iter_data, DetCurr, iLutCurr, Kii, &
                             RealwSign, wAvSign, IterRDMStartCurr, VecSlot, &
                             DetPosition, walkExcitLevel)

        integer, intent(in) :: DetCurr(nel) 
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        integer, intent(inout) :: VecSlot
        real(dp), intent(in) :: Kii
        real(dp), intent(in), dimension(lenof_sign) :: wAvSign, IterRDMStartCurr
        integer, intent(in) :: DetPosition
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign) :: iDie
        real(dp), dimension(lenof_sign) :: CopySign
        integer, intent(in) :: walkExcitLevel
        integer :: i
        character(len=*), parameter :: t_r="walker_death"

        ! Do particles on determinant die? iDie can be both +ve (deaths), or
        ! -ve (births, if shift > 0)
        iDie = attempt_die (DetCurr, Kii, realwSign, WalkExcitLevel)

        IFDEBUG(FCIMCDebug,3) then 
            if(sum(abs(iDie)).ne.0) write(iout,"(A,2f10.5)") "Death: ",iDie(:)
        endif

        ! Update death counter
        iter_data%ndied = iter_data%ndied + min(iDie, abs(RealwSign))
#ifdef __CMPLX
        NoDied = NoDied + sum(min(iDie, abs(RealwSign)))
#else
        NoDied = NoDied + min(iDie, abs(RealwSign))
#endif

        ! Count any antiparticles
        iter_data%nborn = iter_data%nborn + max(iDie - abs(RealwSign), 0.0_dp)
#ifdef __CMPLX
        NoBorn = NoBorn + sum(max(iDie - abs(RealwSign), 0.0_dp))
#else
        NoBorn = NoBorn + max(iDie - abs(RealwSign), 0.0_dp)
#endif

        ! Calculate new number of signed particles on the det.
        CopySign = RealwSign - (iDie * sign(1.0_dp, RealwSign))

        ! In the initiator approximation, abort any anti-particles.
        if (tTruncInitiator .and. any(CopySign /= 0)) then
            do i = 1, lenof_sign
                if (sign(1.0_dp, CopySign(i)) /= &
                        sign(1.0_dp, RealwSign(i))) then
                    NoAborted = NoAborted + abs(CopySign(i))
                    iter_data%naborted(i) = iter_data%naborted(i) &
                                          + abs(CopySign(i))
                    if (test_flag(ilutCurr, flag_is_initiator(i))) &
                        NoAddedInitiators = NoAddedInitiators - 1
                    CopySign(i) = 0
                end if
            end do
        end if

        if (any(CopySign /= 0)) then
            ! Normal method slots particles back in main array at position
            ! vecslot. The list is shuffled up if a particle is destroyed.
            ! For HashWalkerList, the particles don't move, so just adjust
            ! the weight.
            if (tHashWalkerList) then
                call encode_sign (CurrentDets(:,DetPosition), CopySign)
                if (tFillingStochRDMonFly) then
                    call set_av_sgn(DetPosition, wAvSign)
                    call set_iter_occ(DetPosition, IterRDMStartCurr)
                endif
            else
                call encode_bit_rep(CurrentDets(:,VecSlot),iLutCurr,CopySign,extract_flags(iLutCurr))
                call set_det_diagH(VecSlot, Kii)
                if (tFillingStochRDMonFly) then
                    call set_av_sgn(VecSlot, wAvSign)
                    call set_iter_occ(VecSlot, IterRDMStartCurr)
                endif
                VecSlot=VecSlot+1
            endif
        else
            ! All walkers died.
            if(tFillingStochRDMonFly) then
                call det_removed_fill_diag_rdm(CurrentDets(:,DetPosition), DetPosition)
                if (tHashWalkerList) then
                    ! Set the average sign and occupation iteration to zero, so
                    ! that the same contribution will not be added in in
                    ! CalcHashTableStats, if this determinant is not overwritten
                    ! before then
                    global_determinant_data(:,DetPosition) = 0.0_dp
                end if
            endif
            if(tTruncInitiator) then
                ! All particles on this determinant have gone. If the determinant was an initiator, update the stats
                if(test_flag(iLutCurr,flag_is_initiator(1))) then
                    NoAddedInitiators=NoAddedInitiators-1
                    if (tSpawnSpatialInit) call rm_initiator_list (ilutCurr)
                elseif(test_flag(iLutCurr,flag_is_initiator(lenof_sign))) then
                    NoAddedInitiators(inum_runs)=NoAddedInitiators(inum_runs)-1
                    if (tSpawnSpatialInit) call rm_initiator_list (ilutCurr)
                endif
            endif
            if(tHashWalkerList) then
                !Remove the determinant from the indexing list
                call remove_hash_table_entry(HashIndex, DetCurr, DetPosition)
                !Add to the "freeslot" list
                iEndFreeSlot=iEndFreeSlot+1
                FreeSlot(iEndFreeSlot)=DetPosition
                !Encode a null det to be picked up
                call encode_sign(CurrentDets(:,DetPosition),null_part)
            endif
        endif

        !Test - testsuite, RDM still work, both still work with Linscalealgo (all in debug)
        !Null particle not kept if antiparticles aborted.
        !When are the null particles removed?

    end subroutine

    

    ! This routine will find the largest weighted MP1 determinants, from
    ! which we can construct energy level splitting dependant on the sign.
    
    ! This routine will take the particle and sign lists, sort them and then 
    ! compress them so each determinant is only specified once. This requires 
    ! sign-coherent lists.
    ! The 'length' will be returned as the length of the new list.
    SUBROUTINE SortCompressLists(Length,PartList,SignList)
        INTEGER :: Length
        real(dp) :: SignList(Length)
        INTEGER(KIND=n_int) :: PartList(0:NIfTot,Length)
        INTEGER :: i,VecInd,DetsMerged

        call sort (PartList, SignList)

!Now compress the list.
        VecInd=1
        DetsMerged=0
        TotParts=0.0
        IF(Length.gt.0) THEN
            TotParts(1)=TotParts(1)+abs(SignList(1))
        ENDIF
        do i=2,Length
            TotParts(1)=TotParts(1)+abs(SignList(i))
            IF(.not.DetBitEQ(PartList(0:NIfTot,i),PartList(0:NIfTot,VecInd),NIfDBO)) THEN
                VecInd=VecInd+1
                PartList(:,VecInd)=PartList(:,i)
                SignList(VecInd)=SignList(i)
            ELSE
                SignList(VecInd)=SignList(VecInd)+SignList(i)
                DetsMerged=DetsMerged+1
            ENDIF
        enddo

        Length=Length-DetsMerged

    END SUBROUTINE SortCompressLists

    ! This routine will take the particle and sign lists, sort them and then 
    ! compress them so each determinant is only specified once. This requires
    ! sign-coherent lists.

    ! The 'length' will be returned as the length of the new list.
    ! In this version, the hamiltonian matrix elements will be fed through
    ! with the rest of the list and taken with the particles.
    SUBROUTINE SortCompressListswH(Length,PartList,SignList,HList)
        INTEGER :: Length
        INTEGER(KIND=n_int) :: PartList(0:NIfTot,Length)
        real(dp) :: HList(Length), SignList(Length)
        INTEGER :: i,DetsMerged,VecInd

        call sort (PartList, SignList, HList)
!        CALL CheckOrdering(PartList(:,1:Length),SignList(1:Length),Length,.false.)
  
!Now compress...
        VecInd=1
        DetsMerged=0
        TotParts=0
        IF(Length.gt.0) THEN
            TotParts(1)=TotParts(1)+abs(SignList(1))
        ENDIF
        do i=2,Length
            TotParts(1)=TotParts(1)+abs(SignList(i))
            IF(.not.DetBitEQ(PartList(0:NIfTot,i),PartList(0:NIfTot,VecInd),NIfDBO)) THEN
                VecInd=VecInd+1
                PartList(:,VecInd)=PartList(:,i)
                SignList(VecInd)=SignList(i)
                HList(VecInd)=HList(i)
            ELSE
                SignList(VecInd)=SignList(VecInd)+SignList(i)
                DetsMerged=DetsMerged+1
            ENDIF
        enddo

        Length=Length-DetsMerged


!        CALL CheckOrdering(PartList(:,1:Length),CurrentSign(1:Length),Length,.true.)

    END SUBROUTINE SortCompressListswH




    ! TODO: Move to hist.F90
    SUBROUTINE WriteHistogramEnergies()
        INTEGER :: i, io(8)
        real(dp) :: Norm,EnergyBin

        IF(iProcIndex.eq.Root) THEN
            AllHistogramEnergy(:)=0.0_dp
            AllAttemptHist(:)=0.0_dp
            AllSpawnHist(:)=0.0_dp
            AllDoublesHist(:)=0.0_dp
            AllDoublesAttemptHist(:)=0.0_dp
            AllSinglesHist(:)=0.0_dp
            AllSinglesAttemptHist(:)=0.0_dp
            AllSinglesHistOccOcc(:)=0.0_dp
            AllSinglesHistOccVirt(:)=0.0_dp
            AllSinglesHistVirtOcc(:)=0.0_dp
            AllSinglesHistVirtVirt(:)=0.0_dp
        ENDIF
        CALL MPIReduce(HistogramEnergy,MPI_SUM,AllHistogramEnergy)
        CALL MPIReduce(AttemptHist,MPI_SUM,AllAttemptHist)
        CALL MPIReduce(SpawnHist,MPI_SUM,AllSpawnHist)
        CALL MPIReduce(SinglesHist,MPI_SUM,AllSinglesHist)
        CALL MPIReduce(SinglesAttemptHist,MPI_SUM,AllSinglesAttemptHist)
        CALL MPIReduce(DoublesHist,MPI_SUM,AllDoublesHist)
        CALL MPIReduce(DoublesAttemptHist,MPI_SUM,AllDoublesAttemptHist)
        CALL MPIReduce(SinglesHistOccOcc,MPI_SUM,AllSinglesHistOccOcc)
        CALL MPIReduce(SinglesHistOccVirt,MPI_SUM,AllSinglesHistOccVirt)
        CALL MPIReduce(SinglesHistVirtOcc,MPI_SUM,AllSinglesHistVirtOcc)
        CALL MPIReduce(SinglesHistVirtVirt,MPI_SUM,AllSinglesHistVirtVirt)

        IF(iProcIndex.eq.Root) THEN
            AllHistogramEnergy=AllHistogramEnergy/sum(AllHistogramEnergy)
            AllAttemptHist=AllAttemptHist/sum(AllAttemptHist)
            AllSpawnHist=AllSpawnHist/sum(AllSpawnHist)
            AllSinglesAttemptHist=AllSinglesAttemptHist/sum(AllSinglesAttemptHist)
            AllDoublesHist=AllDoublesHist/sum(AllDoublesHist)
            Norm=sum(AllDoublesAttemptHist)
            AllDoublesAttemptHist=AllDoublesAttemptHist/Norm
            do i=1,iOffDiagNoBins
                AllDoublesAttemptHist(i)=AllDoublesAttemptHist(i)/Norm
            enddo
            Norm=0.0_dp
            do i=1,iOffDiagNoBins
                Norm=Norm+AllSinglesHist(i)
            enddo
!            WRITE(iout,*) "AllSinglesHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllSinglesHist(i)=AllSinglesHist(i)/Norm
            enddo
 
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccOcc(i)=AllSinglesHistOccOcc(i)/Norm
            enddo
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccVirt(i)=AllSinglesHistOccVirt(i)/Norm
            enddo
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtOcc(i)=AllSinglesHistVirtOcc(i)/Norm
            enddo
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtVirt(i)=AllSinglesHistVirtVirt(i)/Norm
            enddo

 
            io(1) = get_free_unit()
            OPEN(io(1),FILE='EVERYENERGYHIST',STATUS='UNKNOWN')
            io(2) = get_free_unit()
            OPEN(io(2),FILE='ATTEMPTENERGYHIST',STATUS='UNKNOWN')
            io(3) = get_free_unit()
            OPEN(io(3),FILE='SPAWNENERGYHIST',STATUS='UNKNOWN')

            EnergyBin=BinRange/2.0_dp
            do i=1,iNoBins
                IF(AllHistogramEnergy(i).gt.0.0_dp) WRITE(io(1),*) EnergyBin, AllHistogramEnergy(i)
                IF(AllAttemptHist(i).gt.0.0_dp) WRITE(io(2),*) EnergyBin, AllAttemptHist(i)
                IF(AllSpawnHist(i).gt.0.0_dp) WRITE(io(3),*) EnergyBin, AllSpawnHist(i)
                EnergyBin=EnergyBin+BinRange
            enddo
            CLOSE(io(1))
            CLOSE(io(2))
            CLOSE(io(3))
            OPEN(io(1),FILE='SINGLESHIST',STATUS='UNKNOWN')
            OPEN(io(2),FILE='ATTEMPTSINGLESHIST',STATUS='UNKNOWN')
            OPEN(io(3),FILE='DOUBLESHIST',STATUS='UNKNOWN')
            io(4) = get_free_unit()
            OPEN(io(4),FILE='ATTEMPTDOUBLESHIST',STATUS='UNKNOWN')
            io(5) = get_free_unit()
            OPEN(io(5),FILE='SINGLESHISTOCCOCC',STATUS='UNKNOWN')
            io(6) = get_free_unit()
            OPEN(io(6),FILE='SINGLESHISTOCCVIRT',STATUS='UNKNOWN')
            io(7) = get_free_unit()
            OPEN(io(7),FILE='SINGLESHISTVIRTOCC',STATUS='UNKNOWN')
            io(8) = get_free_unit()
            OPEN(io(8),FILE='SINGLESHISTVIRTVIRT',STATUS='UNKNOWN')

            EnergyBin=-OffDiagMax+OffDiagBinRange/2.0_dp
            do i=1,iOffDiagNoBins
                IF(AllSinglesHist(i).gt.0.0_dp) WRITE(io(1),*) EnergyBin, AllSinglesHist(i)
                IF(AllSinglesAttemptHist(i).gt.0.0_dp) WRITE(io(2),*) EnergyBin, AllSinglesAttemptHist(i)
                IF(AllDoublesHist(i).gt.0.0_dp) WRITE(io(3),*) EnergyBin, AllDoublesHist(i)
                IF(AllDoublesAttemptHist(i).gt.0.0_dp) WRITE(io(4),*) EnergyBin, AllDoublesAttemptHist(i)
                IF(AllSinglesHistOccOcc(i).gt.0.0_dp) WRITE(io(5),*) EnergyBin, AllSinglesHistOccOcc(i)
                IF(AllSinglesHistOccVirt(i).gt.0.0_dp) WRITE(io(6),*) EnergyBin, AllSinglesHistOccVirt(i)
                IF(AllSinglesHistVirtOcc(i).gt.0.0_dp) WRITE(io(7),*) EnergyBin, AllSinglesHistVirtOcc(i)
                IF(AllSinglesHistVirtVirt(i).gt.0.0_dp) WRITE(io(8),*) EnergyBin, AllSinglesHistVirtVirt(i)
                EnergyBin=EnergyBin+OffDiagBinRange
!                WRITE(6,*) i
            enddo

            CLOSE(io(1))
            CLOSE(io(2))
            CLOSE(io(3))
            CLOSE(io(4))
            CLOSE(io(5))
            CLOSE(io(6))
            CLOSE(io(7))
            CLOSE(io(8))
        ENDIF

    END SUBROUTINE WriteHistogramEnergies


    subroutine gennct(n,t,icomb)
       integer n,t
       integer c(t+2)
       integer j
       integer icomb
       integer iunit
       iunit = get_free_unit()
       icomb=0
       open(iunit,file='COMBINATIONS',STATUS='UNKNOWN')
!      WRITE(iout,*) ' Writing combinations to file COMBINATIONS'
       do j=1,t
           c(j)=j-1
       enddo
       c(t+1)=n
       c(t+2)=0
20     continue
!      visit c(t)
       
       icomb=icomb+1
!      write(iunit,'(20i3)' ) (c(j)+1,j=1,t)
       
       write(iunit,'(20i3)' ) (c(j),j=1,t)
       
       do j=1,n
           if(c(j+1).ne.(c(j)+1)) goto 30
           c(j)=j-1
       enddo
30     continue
       if(j.gt.t) then 
!           write(iout,*) ' Generated combinations:',ICOMB
           CLOSE(iunit)
           RETURN
       endif
       c(j)=c(j)+1
       goto 20
       
   end subroutine gennct
       

!Similar to WriteHistogram, but will only print out in order of maximum component, and only the averaged wavefunction
    SUBROUTINE PrintFCIMCPsi()
        use DetCalcData , only : FCIDets
        INTEGER :: i,nI(NEl),ExcitLevel,j, iunit
        real(dp) :: norm2
        real(dp), dimension(lenof_sign) :: norm1, norm

        CALL MPISumAll(Histogram,AllHistogram)
        norm1=0.0_dp
        do i=1,Det
            do j=1,lenof_sign
                norm1(j)=norm1(j)+AllHistogram(j,i)**2
            enddo
        enddo
#ifdef __CMPLX
        norm2=SQRT(sum(norm1))
#else
        norm1=SQRT(norm1)
#endif
        WRITE(iout,*) "Total FCIMC Wavefuction normalisation:",norm1
        do i=1,Det
            do j=1,lenof_sign
#ifdef __CMPLX
                AllHistogram(j,i)=AllHistogram(j,i)/norm2
#else
                AllHistogram(j,i)=AllHistogram(j,i)/norm1(j)
#endif
            enddo
        enddo

        iunit = 0
        IF(tPrintFCIMCPsi) THEN
!Order and print wavefunction

            
            IF(iProcIndex.eq.0) THEN

                ! We now want to order AllHistogram, taking the corresponding
                ! element(s) of FCIDets with it...
                call sort (AllHistogram, FCIDets)
                
                OPEN(iunit,FILE='FCIMCPsi',STATUS='UNKNOWN')

                norm=0.0_dp
                do i=1,Det
                    do j=1,lenof_sign
                        norm(j)=norm(j)+AllHistogram(j,i)**2
                    enddo
!write out FCIMC Component weight (normalised), current normalisation, excitation level
                    ExcitLevel = FindBitExcitLevel(iLutHF, FCIDets(:,i), nel)
                    CALL decode_bit_det(nI,FCIDets(0:NIfTot,i))
#ifdef __CMPLX
                    WRITE(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,sum(norm)
#else
                    WRITE(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,norm(1)
#endif
                    do j=1,NEl-1
                        WRITE(iunit,"(I5)",advance='no') nI(j)
                    enddo
                    WRITE(iunit,"(I5)") nI(NEl)
                enddo

                CLOSE(iunit)

            ENDIF
        ENDIF

    END SUBROUTINE PrintFCIMCPsi

!This routine will write out the average wavevector from the spawning run up until now.
    SUBROUTINE WriteHistogram()
        INTEGER :: i,IterRead, io1, io2, io3, Tot_No_Unique_Dets,ierr,val,j
        integer :: posn, FinalPop
        real(dp) :: ShiftRead,AllERead,NumParts,DDOT
        real(dp),dimension(lenof_sign) :: norm,norm1,norm2,norm3
        real(dp) :: norm_c,norm1_c,norm2_c,norm3_c
        real(dp) :: AvVarEnergy, VarEnergy, GroundE_Ever, ProjGroundE
        real(dp) , allocatable :: CKN(:),HOrderedHist(:),HOrderedInstHist(:)
        integer , allocatable :: ExpandedWalkerDets(:,:)
        CHARACTER(len=22) :: abstr,abstr2
        LOGICAL :: exists
        character(len=*), parameter :: t_r='WriteHistogram'

!This will open a file called SpawnHist-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='SpawnHist-'//adjustl(abstr)
        IF(iProcIndex.eq.0) THEN
            WRITE(iout,*) "Writing out the average wavevector up to iteration number: ", Iter
            CALL neci_flush(iout)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistogram(:,:)=0.0_dp
            AllInstHist(:,:)=0.0_dp
            AllAvAnnihil(:,:)=0.0_dp
            AllInstAnnihil(:,:)=0.0_dp
        ENDIF

        CALL MPIReduce(Histogram,MPI_SUM,AllHistogram)
        CALL MPIReduce(InstHist,MPI_SUM,AllInstHist)
        CALL MPIReduce(InstAnnihil,MPI_SUM,AllInstAnnihil)
        CALL MPIReduce(AvAnnihil,MPI_SUM,AllAvAnnihil)

        BeforeNormHist(:) = AllHistogram(1,:)
        
        IF(iProcIndex.eq.0) THEN

            norm=0.0_dp
            norm1=0.0_dp
            norm2=0.0_dp
            norm3=0.0_dp
            do i=1,Det
                do j=1,lenof_sign
                    norm(j)=norm(j)+AllHistogram(j,i)**2
                    norm1(j)=norm1(j)+AllInstHist(j,i)**2
                    norm2(j)=norm2(j)+AllInstAnnihil(j,i)**2
                    norm3(j)=norm3(j)+AllAvAnnihil(j,i)**2
                enddo
            enddo
#ifdef __CMPLX
            norm_c=SQRT(sum(norm))
            norm1_c=SQRT(sum(norm1))
            norm2_c=SQRT(sum(norm2))
            norm3_c=SQRT(sum(norm3))
#else
            norm=SQRT(norm)
            norm1=SQRT(norm1)
            norm2=SQRT(norm2)
            norm3=SQRT(norm3)
#endif
            do i=1,Det
                do j=1,lenof_sign
#ifdef __CMPLX
                    AllHistogram(j,i)=AllHistogram(j,i)/norm_c
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1_c
                    IF(norm2_c.ne.0.0_dp) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(j,i)/norm2_c
                    ENDIF
                    IF(norm3_c.ne.0.0_dp) THEN
                        AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3_c
                    ENDIF
#else
                    AllHistogram(j,i)=AllHistogram(j,i)/norm(j)
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1(j)
                    IF(norm2(j).ne.0.0_dp) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(1,i)/norm2(j)
                    ENDIF
                    IF(norm3(j).ne.0.0_dp) THEN
                    AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3(j)
                    ENDIF
#endif
                enddo
            enddo
            
            io1 = get_free_unit()
            OPEN(io1,FILE=abstr,STATUS='UNKNOWN')

            abstr=''
            write(abstr,'(I12)') Iter-iWriteHistEvery
            abstr='Energies-'//adjustl(abstr)

            abstr2=''
            write(abstr2,'(I12)') Iter
            abstr2='Energies-'//adjustl(abstr2)
            io2 = get_free_unit()
            OPEN(io2,FILE=abstr2,STATUS='UNKNOWN')

            INQUIRE(FILE=abstr,EXIST=exists)
            IF(exists) THEN
                io3 = get_free_unit()
                OPEN(io3,FILE=abstr,STATUS='OLD',POSITION='REWIND',ACTION='READ')
                do while(.true.)
                    READ(io3,"(I13,3G25.16)",END=99) IterRead,ShiftRead,AllERead,NumParts
                    WRITE(io2,"(I13,3G25.16)") IterRead,ShiftRead,AllERead,NumParts
                enddo
99              CONTINUE
#ifdef __CMPLX
                IF(AllHFCyc(1).eq.0.0_dp) THEN
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllERead,SUM(AllTotPartsOld)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,SUM(AllTotPartsOld)
                ENDIF
#else
                IF(AllHFCyc(1).eq.0.0_dp) THEN
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllERead,AllTotPartsOld(1)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumCyc(1)/AllHFCyc(1),AllTotPartsOld(1)
                ENDIF
#endif
                CLOSE(io2)
                CLOSE(io3)

            ELSE
                OPEN(io2,FILE=abstr2,STATUS='UNKNOWN')
#ifdef __CMPLX
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,SUM(AllTotPartsOld)
#else
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumCyc(1)/AllHFCyc(1),AllTotPartsOld(1)
#endif
                CLOSE(io2)
            ENDIF


            norm=0.0_dp
            norm1=0.0_dp
            Tot_No_Unique_Dets = 0
            do i=1,Det

                posn=binary_search(CurrentDets(:,1:TotWalkers), FCIDETS(:,i), NifD+1)

                if (posn.lt.0) then
                    FinalPop = 0
                else
                    FinalPop = int(CurrentDets(NifD+1,posn))
                endif

                do j=1,lenof_sign
                    norm(j)=norm(j)+AllHistogram(j,i)**2
                    norm1(j)=norm1(j)+AllAvAnnihil(j,i)**2
                enddo
                IF(lenof_sign.eq.1) THEN
                    WRITE(io1,"(I13,6G25.16,I13,G25.16)") i, &
                          AllHistogram(1,i), norm, AllInstHist(1,i), &
                          AllInstAnnihil(1,i), AllAvAnnihil(1,i), norm1, &
                          FinalPop, BeforeNormHist(i)
                ELSE
#ifdef __CMPLX
                    WRITE(io1,"(I13,6G25.16)") i, AllHistogram(1,i), sum(norm), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), sum(norm1)
#else
                    WRITE(io1,"(I13,6G25.16)") i, AllHistogram(1,i), norm(1), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), norm1(1)
#endif
                ENDIF
                IF(AllHistogram(1,i).ne.0.0_dp) Tot_No_Unique_Dets = Tot_No_Unique_Dets + 1
            enddo
            if(tCalcVariationalEnergy) then
                !Calculate the variational FCIMC energy

                !Check that Det eq NDet
                if(Det.ne.NDet) call stop_all(t_r,'Error')
                allocate(CKN(1:Det))
                CKN = 0.0_dp
                !We need to reorder the vectors so that they are in the same order as the hamiltonian
                allocate(HOrderedHist(1:Det))
                HOrderedHist(:) = AllHistogram(1,:)
                allocate(HOrderedInstHist(1:Det))
                HOrderedInstHist(:) = AllInstHist(1,:)
                call sort(ReIndex(1:Det),HOrderedHist(1:Det),HOrderedInstHist(1:Det))

!                write(6,*) ReIndex(1:10)    !Must be in increasing order
!                write(6,*) HOrderedHist(:)

                call my_hpsi(Det,1,NROW,LAB,HAMIL,HOrderedHist,CKN,.true.)
                AvVarEnergy = DDOT(Det,HOrderedHist,1,CKN,1)
                
                CKN = 0.0_dp
                call my_hpsi(Det,1,NROW,LAB,HAMIL,HOrderedInstHist,CKN,.true.)
                VarEnergy = DDOT(Det,HOrderedInstHist,1,CKN,1)

                deallocate(CKN)
                deallocate(HOrderedHist,HOrderedInstHist)
                if(tDiagAllSpaceEver) then

                    allocate(ExpandedWalkerDets(NEl,Tot_No_Unique_Dets),stat=ierr)
                    val=1
                    do i=1,Det
                        if(AllHistogram(1,i).ne.0.0_dp) then
                            call decode_bit_det(ExpandedWalkerDets(:,val),FCIDets(0:NIfTot,i))
                            val=val+1
                        endif
                    enddo
                    if((val-1).ne.Tot_No_Unique_Dets) call stop_all(t_r,'Wrong counting')

                    call LanczosFindGroundE(ExpandedWalkerDets,Tot_No_Unique_Dets,GroundE_Ever,ProjGroundE,.true.)
                    if(abs(GroundE_Ever-ProjGroundE).gt.1.0e-7_dp) call stop_all(t_r,'Why do these not agree?!')
                    write(Tot_Unique_Dets_Unit,"(2I14,3G25.15)") Iter,Tot_No_Unique_Dets,AvVarEnergy,VarEnergy,GroundE_Ever
                else
                    write(Tot_Unique_Dets_Unit,"(2I14,2G25.15)") Iter,Tot_No_Unique_Dets,AvVarEnergy,VarEnergy
                endif
            else
                if(tDiagAllSpaceEver) then

                    allocate(ExpandedWalkerDets(NEl,Tot_No_Unique_Dets),stat=ierr)
                    val=1
                    do i=1,Det
                        if(AllHistogram(1,i).ne.0.0_dp) then
                            call decode_bit_det(ExpandedWalkerDets(:,val),FCIDets(0:NIfTot,i))
                            val=val+1
                        endif
                    enddo
                    if((val-1).ne.Tot_No_Unique_Dets) call stop_all(t_r,'Wrong counting')

                    call LanczosFindGroundE(ExpandedWalkerDets,Tot_No_Unique_Dets,GroundE_Ever,ProjGroundE,.true.)
                    if(abs(GroundE_Ever-ProjGroundE).gt.1.0e-7_dp) call stop_all(t_r,'Why do these not agree?!')
                    write(Tot_Unique_Dets_Unit,"(2I14,3G25.15)") Iter,Tot_No_Unique_Dets,GroundE_Ever
                else
                    write(Tot_Unique_Dets_Unit,"(2I14)") Iter, Tot_No_Unique_Dets
                endif
            endif

            CLOSE(io1)
        ENDIF
        InstHist(:,:)=0.0_dp
        InstAnnihil(:,:)=0.0_dp

    END SUBROUTINE WriteHistogram

!This routine will write out the average hamiltonian from the spawning run up until now.
    SUBROUTINE WriteHamilHistogram()
        INTEGER :: i,j
        integer :: iunit
        CHARACTER(len=22) :: abstr

!This will open a file called HamilHist-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='HamilHist-'//adjustl(abstr)
        IF(iProcIndex.eq.0) THEN
            WRITE(iout,*) "Writing out the average hamiltonian up to iteration number: ", Iter
            CALL neci_flush(iout)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistHamil(:,:)=0.0_dp
            AllAvHistHamil(:,:)=0.0_dp
        ENDIF

        CALL MPIReduce(HistHamil,MPI_SUM,AllHistHamil)
        CALL MPIReduce(AvHistHamil,MPI_SUM,AllAvHistHamil)
        
        IF(iProcIndex.eq.0) THEN
!How do we normalise this!
            iunit = get_free_unit()
            OPEN(iunit,FILE=abstr,STATUS='UNKNOWN')
            do i=1,Det
                do j=1,Det
                    WRITE(iunit,*) j,i,AllAvHistHamil(j,i),AllHistHamil(j,i)
                enddo
                WRITE(iunit,*) ""
            enddo
            CLOSE(iunit)
        ENDIF
        HistHamil(:,:)=0.0_dp

    END SUBROUTINE WriteHamilHistogram



    ! TODO: COMMENTING
    subroutine iter_diagnostics ()

        character(*), parameter :: this_routine = 'iter_diagnostics'
        character(*), parameter :: t_r = this_routine
        real(dp) :: mean_walkers
        integer :: part_type

        ! Update the total imaginary time passed
        TotImagTime = TotImagTime + StepsSft * Tau

        ! Set Iter time to equal the average time per iteration in the
        ! previous update cycle.
        IterTime = IterTime / real(StepsSft,sp)

        ! Calculate the acceptance ratio
        AccRat = real(Acceptances, dp) / SumWalkersCyc


#ifndef __CMPLX
        if (.not. tKP_FCIQMC) then
            do part_type = 1, lenof_sign
                if ((.not.tFillingStochRDMonFly).or.(inum_runs.eq.1)) then
                    if (AllNoAtHF(part_type) < 0.0_dp) then
                        root_print 'No. at HF < 0 - flipping sign of entire ensemble &
                                   &of particles in simulation: ', part_type
                        root_print AllNoAtHF(part_type)

                        ! And do the flipping
                        call FlipSign(part_type)
                        AllNoatHF(part_type) = -AllNoatHF(part_type)
                        NoatHF(part_type) = -NoatHF(part_type)

                        if (tFillingStochRDMonFly) then
                            ! Want to flip all the averaged signs.
                            AvNoatHF = -AVNoatHF
                            InstNoatHF(part_type) = -InstNoatHF(part_type)
                        end if
                    endif
                end if
            end do
        end if
#endif

        if (iProcIndex == Root) then
            ! Have all of the particles died?
#ifdef __CMPLX
            if (AllTotwalkers == 0)  then
                write(iout,"(A)") "All particles have died. Restarting."
                tRestart=.true.
            else
                tRestart=.false.
            endif
#else
            if ((AllTotParts(1).eq.0).or.(AllTotParts(inum_runs).eq.0))  then
                write(iout,"(A)") "All particles have died. Restarting."
                tRestart=.true.
            else
                tRestart=.false.
            endif
            !TODO CMO: Work out how to wipe the walkers on the second population if double run
#endif
        endif
        call MPIBCast(tRestart)
        if(tRestart) then
!Initialise variables for calculation on each node
            Iter=1
            CALL DeallocFCIMCMemPar()
            IF(iProcIndex.eq.Root) THEN
                CLOSE(fcimcstats_unit)
                if (inum_runs.eq.2) CLOSE(fcimcstats_unit2)
                IF(tTruncInitiator) CLOSE(initiatorstats_unit)
                IF(tLogComplexPops) CLOSE(complexstats_unit)
            ENDIF
            IF(TDebug) CLOSE(11)
            CALL SetupParameters()
            CALL InitFCIMCCalcPar()
            if (tFCIMCStats2) then
                call write_fcimcstats2(iter_data_fciqmc, initial=.true.)
                call write_fcimcstats2(iter_data_fciqmc)
            else
                call WriteFciMCStatsHeader()
                ! Prepend a # to the initial status line so analysis doesn't pick up
                ! repetitions in the FCIMCStats or INITIATORStats files from restarts.
                if (iProcIndex == root) then
                    write (fcimcstats_unit,'("#")', advance='no')
                    if (inum_runs == 2) &
                        write(fcimcstats_unit2, '("#")', advance='no')
                    write (initiatorstats_unit,'("#")', advance='no')
                end if
                call WriteFCIMCStats()
            end if
            return
        endif

        if(iProcIndex.eq.Root) then
! AJWT dislikes doing this type of if based on a (seeminly unrelated) input option, but can't see another easy way.
!  TODO:  Something to make it better
            if(.not.tCCMC) then
               ! Check how balanced the load on each processor is (even though
               ! we cannot load balance with direct annihilation).
               WalkersDiffProc = int(MaxWalkersProc - MinWalkersProc,sizeof_int)
               ! Do the same for number of particles
               PartsDiffProc = int(MaxPartsProc - MinPartsProc, sizeof_int)

               mean_walkers = AllTotWalkers / real(nNodes,dp)
               if (WalkersDiffProc > nint(mean_walkers / 10.0_dp) .and. &
                   sum(AllTotParts) > real(nNodes * 500, dp)) then
                   root_write (iout, '(a, i13,a,2i11)') &
                       'Potential load-imbalance on iter ',iter + PreviousCycles,' Min/Max determinants on node: ', &
                       MinWalkersProc,MaxWalkersProc
               endif
            endif
        endif

    end subroutine iter_diagnostics

    subroutine population_check ()
        use HPHFRandExcitMod, only: ReturnAlphaOpenDet
        integer :: pop_highest, proc_highest
        real(dp) :: pop_change, old_Hii
        integer :: det(nel), i, error
        integer(int32) :: int_tmp(2)
        logical :: tSwapped
        HElement_t :: h_tmp

        if (tCheckHighestPop) then

            ! Obtain the determinant (and its processor) with the highest
            ! population. To keep this simple, do it only for set 1 if using double run,
            ! as we need to keep a consistent HF det for the two runs.
            call MPIAllReduceDatatype ((/int(iHighestPop,int32), int(iProcIndex,int32)/), 1, &
                                       MPI_MAXLOC, MPI_2INTEGER, int_tmp)
            pop_highest = int_tmp(1)
            proc_highest = int_tmp(2)

            ! NB: the use if int(iHighestPop) obviously introduces a small amount of error
            ! by ignoring the fractional population here

            ! How many walkers do we need to switch dets?
            
            ! If doing a double run, we only test population 1. abs_sign considers element 1
            ! unless we're running the complex code.
            if((lenof_sign.eq.2).and.(inum_runs.eq.1)) then 
                pop_change = FracLargerDet * abs_sign(AllNoAtHF)
            else
                pop_change = FracLargerDet * abs(AllNoAtHF(1))
            endif
!            write(iout,*) "***",AllNoAtHF,FracLargerDet,pop_change, pop_highest,proc_highest
            if (pop_change < pop_highest .and. pop_highest > 50) then

                ! Write out info!
                    root_print 'Highest weighted determinant not reference &
                               &det: ', pop_highest, abs_sign(AllNoAtHF)
                    

                ! Are we changing the reference determinant?
                if (tChangeProjEDet) then
                    ! Communicate the change to all dets and print out.
                    call MPIBcast (HighestPopDet(0:NIfTot), NIfTot+1, proc_highest)
                    iLutRef = 0
                    iLutRef(0:NIfDBO) = HighestPopDet(0:NIfDBO)
                    call decode_bit_det (ProjEDet, iLutRef)
                    write (iout, '(a)', advance='no') 'Changing projected &
                          &energy reference determinant for the next update cycle to: '
                    call write_det (iout, ProjEDet, .true.)
                    tRef_Not_HF = .true.

                    if(tHPHF) then
                        if(.not.TestClosedShellDet(iLutRef)) then
                            !Complications. We are now effectively projecting onto a LC of two dets.
                            !Ensure this is done correctly.
                            if(.not.Allocated(RefDetFlip)) then
                                allocate(RefDetFlip(NEl))
                                allocate(iLutRefFlip(0:NIfTot))
                                RefDetFlip = 0
                                iLutRefFlip = 0
                            endif
                            call ReturnAlphaOpenDet(ProjEDet,RefDetFlip,iLutRef,iLutRefFlip,.true.,.true.,tSwapped)
                            if(tSwapped) then
                                !The iLutRef should already be the correct one, since it was obtained by the normal calculation!
                                call stop_all("population_check","Error in changing reference determinant to open shell HPHF")
                            endif
                            write(iout,"(A)") "Now projecting onto open-shell HPHF as a linear combo of two determinants..."
                            tSpinCoupProjE=.true.
                        endif
                    else
                        tSpinCoupProjE=.false.  !In case it was already on, and is now projecting onto a CS HPHF.
                    endif

                    ! We can't use Brillouin's theorem if not a converged,
                    ! closed shell, ground state HF det.
                    tNoBrillouin = .true.
                    root_print "Ensuring that Brillouin's theorem is no &
                               &longer used."

                    ! Update the reference energy
                    old_Hii = Hii
                    if (tHPHF) then
                        h_tmp = hphf_diag_helement (ProjEDet, iLutRef)
                    else
                        h_tmp = get_helement (ProjEDet, ProjEDet, 0)
                    endif
                    Hii = real(h_tmp, dp)
                    write (iout, '(a, g25.15)') 'Reference energy now set to: ',&
                                             Hii

                    ! Reset averages
                    SumENum(:)=0
                    sum_proje_denominator(:) = 0
                    cyc_proje_denominator(:) = 0
                    SumNoatHF(:) = 0.0_dp
                    VaryShiftCycles(:) = 0
                    SumDiagSft(:) = 0
                    root_print 'Zeroing all energy estimators.'

                    !Since we have a new reference, we must block only from after this point
                    iBlockingIter = Iter + PreviousCycles

                    ! Regenerate all the diagonal elements relative to the
                    ! new reference det.
                    write (iout,*) 'Regenerating the stored diagonal HElements &
                                &for all walkers.'
                    do i = 1, int(Totwalkers,sizeof_int)
                        call decode_bit_det (det, CurrentDets(:,i))
                        if (tHPHF) then
                            h_tmp = hphf_diag_helement (det, CurrentDets(:,i))
                        else
                            h_tmp = get_helement (det, det, 0)
                        endif
                        call set_det_diagH(i, real(h_tmp, dp) - Hii)
                    enddo
                    if (tSemiStochastic) call recalc_core_hamil_diag(old_Hii, Hii)

                    ! Reset values introduced in soft_exit (CHANGEVARS)
                    if (tCHeckHighestPopOnce) then
                        tChangeProjEDet = .false.
                        tCheckHighestPop = .false.
                        tCheckHighestPopOnce = .false.
                    endif

                ! Or are we restarting the calculation with the reference 
                ! det switched?
#ifdef __CMPLX
                elseif (tRestartHighPop .and. &
                        iRestartWalkNum < sum(AllTotParts)) then
#else
                elseif (tRestartHighPop .and. &
                        iRestartWalkNum < AllTotParts(1)) then
#endif
                    
                    ! Broadcast the changed det to all processors
                    call MPIBcast (HighestPopDet, NIfTot+1, proc_highest)
                    iLutRef = 0
                    iLutRef(0:NIfDBO) = HighestPopDet(0:NIfDBO)
                    tRef_Not_HF = .true.

                    call decode_bit_det (ProjEDet, iLutRef)
                    write (iout, '(a)', advance='no') 'Changing projected &
                             &energy reference determinant to: '
                    call write_det (iout, ProjEDet, .true.)

                    ! We can't use Brillouin's theorem if not a converged,
                    ! closed shell, ground state HF det.
                    tNoBrillouin = .true.
                    root_print "Ensuring that Brillouin's theorem is no &
                               &longer used."
                    
                    ! Update the reference energy
                    if (tHPHF) then
                        h_tmp = hphf_diag_helement (ProjEDet, iLutRef)
                    else
                        h_tmp = get_helement (ProjEDet, ProjEDet, 0)
                    endif
                    Hii = real(h_tmp, dp)
                    write (iout, '(a, g25.15)') 'Reference energy now set to: ',&
                                             Hii

                    ! Reset values introduced in soft_exit (CHANGEVARS)
                    if (tCHeckHighestPopOnce) then
                        tChangeProjEDet = .false.
                        tCheckHighestPop = .false.
                        tCheckHighestPopOnce = .false.
                    endif

                    call ChangeRefDet (ProjEDet)
                endif

            endif
        endif
                    
    end subroutine

    subroutine collate_iter_data (iter_data, tot_parts_new, tot_parts_new_all)
        integer :: int_tmp(5+2*lenof_sign), proc, pos, i
        real(dp) :: sgn(lenof_sign)
        HElement_t :: helem_tmp(3*inum_runs)
        HElement_t :: real_tmp(2*inum_runs) !*lenof_sign
        integer(int64) :: int64_tmp(8),TotWalkersTemp
        type(fcimc_iter_data) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new
        real(dp), dimension(lenof_sign), intent(out) :: tot_parts_new_all
        character(len=*), parameter :: this_routine='collate_iter_data'
        real(dp), dimension(max(lenof_sign,inum_runs)) :: RealAllHFCyc
        real(dp), dimension(inum_runs) :: all_norm_psi_squared, all_norm_semistoch_squared
        real(dp) :: bloom_sz_tmp(0:2)
        integer :: run
    
        ! Communicate the integers needing summation

        call MPIReduce(SpawnFromSing, MPI_SUM, AllSpawnFromSing)
        call MPIReduce(iter_data%update_growth, MPI_SUM, iter_data%update_growth_tot)
        call MPIReduce(NoBorn, MPI_SUM, AllNoBorn)
        call MPIReduce(NoDied, MPI_SUM, AllNoDied)
        call MPIReduce(HFCyc, MPI_SUM, RealAllHFCyc)
        call MPIReduce(NoAtDoubs, MPI_SUM, AllNoAtDoubs)
        call MPIReduce(Annihilated, MPI_SUM, AllAnnihilated)
        
        do run=1,inum_runs
            AllHFCyc(run)=ARR_RE_OR_CPLX(RealAllHFCyc,run)
        enddo
        
        ! Integer summations required for the initiator method
        if (tTruncInitiator) then
            call MPISum ((/NoAddedInitiators(1), NoInitDets(1), &
                           NoNonInitDets(1), NoExtraInitdoubs(1), InitRemoved(1)/),&
                          int64_tmp(1:5))
            AllNoAddedInitiators(1) = int64_tmp(1)
            AllNoInitDets(1) = int64_tmp(2)
            AllNoNonInitDets(1) = int64_tmp(3)
            AllNoExtraInitDoubs(1) = int64_tmp(4)
            AllInitRemoved(1) = int64_tmp(5)

            call MPIReduce(NoAborted, MPI_SUM, AllNoAborted)
            call MPIReduce(NoRemoved, MPI_SUM, AllNoRemoved)
            call MPIReduce(NoNonInitWalk, MPI_SUM, AllNoNonInitWalk)
            call MPIReduce(NoInitWalk, MPI_SUM, AllNoInitWalk)
        endif

        ! 64bit integers
        !Remove the holes in the main list when wanting the number of uniquely occupied determinants
        !this should only change the number for tHashWalkerList
        TotWalkersTemp=TotWalkers-HolesInList
        call MPIReduce(TotwalkersTemp, MPI_SUM, AllTotWalkers)
        call MPIReduce(norm_psi_squared,MPI_SUM,all_norm_psi_squared)
        call MPIReduce(norm_semistoch_squared,MPI_SUM,all_norm_semistoch_squared)
        call MPIReduce(Totparts,MPI_SUM,AllTotParts)
        call MPIReduce(tot_parts_new,MPI_SUM,tot_parts_new_all)
#ifdef __CMPLX
        norm_psi = sqrt(sum(all_norm_psi_squared))
        norm_semistoch = sqrt(sum(all_norm_semistoch_squared))
#else
        norm_psi = sqrt(all_norm_psi_squared)
        norm_semistoch = sqrt(all_norm_semistoch_squared)
#endif
        
        call MPIReduce(SumNoatHF, MPI_SUM, AllSumNoAtHF)
        ! HElement_t values (Calculates the energy by summing all on HF and 
        ! doubles)

        call MPISum ((/ENumCyc, SumENum, ENumCycAbs/), helem_tmp)
        AllENumCyc(:) = helem_tmp(1:inum_runs)
        AllSumENum(:) = helem_tmp(1+inum_runs:2*inum_runs)
        AllENumCycAbs(:) = helem_tmp(1+2*inum_runs:3*inum_runs)
        
        ! Deal with particle blooms
!        if (tSpinProjDets) then
!            call MPISum(bloom_count(0:2), all_bloom_count(0:2))
!            call MPIReduce_inplace(bloom_sizes(0:2), MPI_MAX)
!        else
            call MPISum(bloom_count(1:2), all_bloom_count(1:2))
            call MPIReduce(bloom_sizes(1:2), MPI_MAX, bloom_sz_tmp(1:2))
            bloom_sizes(1:2) = bloom_sz_tmp(1:2)
!        end if

        ! real(dp) values
        call MPISum((/cyc_proje_denominator, sum_proje_denominator/),real_tmp)
        all_cyc_proje_denominator = real_tmp(1:inum_runs)!(1:lenof_sign)
        all_sum_proje_denominator = real_tmp(1+inum_runs:2*inum_runs)!(lenof_sign+1:2*lenof_sign)

        ! Max/Min values (check load balancing)
        call MPIReduce (TotWalkersTemp, MPI_MAX, MaxWalkersProc)
        call MPIReduce (TotWalkersTemp, MPI_MIN, MinWalkersProc)
        call MPIReduce (max_cyc_spawn, MPI_MAX, all_max_cyc_spawn)
        !call MPIReduce (sum(TotParts), MPI_MAX, MaxPartsProc)
        !call MPIReduce (sum(TotParts), MPI_MIN, MinPartsProc)

        ! We need the total number on the HF and SumWalkersCyc to be valid on
        ! ALL processors (Both double precision reals)
        call MPISumAll (NoatHF, AllNoatHF)
        call MPISumAll (SumWalkersCyc, AllSumWalkersCyc)

        !        WRITE(iout,*) "***",iter_data%update_growth_tot,AllTotParts-AllTotPartsOld

        if (tSearchTau .and. (.not. tFillingStochRDMonFly)) &
            call update_tau()

        !TODO CMO:Make sure these are length 2 as well
        if (tTrialWavefunction) then
            call MPIAllReduce(trial_numerator, MPI_SUM, tot_trial_numerator)
            call MPIAllReduce(trial_denom, MPI_SUM, tot_trial_denom)
        end if
        
#ifdef __DEBUG
        !Write this 'ASSERTROOT' out explicitly to avoid line lengths problems
        if ((iProcIndex == root) .and. .not. tSpinProject .and. &
         all(abs(iter_data%update_growth_tot-(AllTotParts-AllTotPartsOld)) > 1.0e-5)) then
            write(iout,*) "update_growth: ",iter_data%update_growth_tot
            write(iout,*) "AllTotParts: ",AllTotParts
            write(iout,*) "AllTotPartsOld: ", AllTotPartsOld
            call stop_all (this_routine, &
                "Assertation failed: all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld)")
        endif
#endif
    
    end subroutine collate_iter_data

    subroutine update_shift (iter_data)
        use CalcData, only : tInstGrowthRate
     
        type(fcimc_iter_data), intent(in) :: iter_data
        integer(int64) :: tot_walkers
        logical, dimension(inum_runs) :: tReZeroShift
        real(dp) :: AllGrowRateRe, AllGrowRateIm
        real(dp), dimension(inum_runs)  :: AllHFGrowRate
        real(dp), dimension(lenof_sign) :: denominator, all_denominator
        integer :: error, i, proc, pos, run
        logical, dimension(inum_runs) :: defer_update
        logical :: start_varying_shift

        ! Normally we allow the shift to vary depending on the conditions
        ! tested. Sometimes we want to defer this to the next cycle...
        defer_update(:) = .false.

!        call neci_flush(iout)
!        CALL MPIBarrier(error)

        ! collate_iter_data --> The values used are only valid on Root
        if (iProcIndex == Root) then
            ! Calculate the growth rate
!            WRITE(iout,*) "iter_data%nborn: ",iter_data%nborn(:)
!            WRITE(iout,*) "iter_data%ndied: ",iter_data%ndied(:)
!            WRITE(iout,*) "iter_data%nannihil: ",iter_data%nannihil(:)
!            WRITE(iout,*) "iter_data%naborted: ",iter_data%naborted(:)
!            WRITE(iout,*) "iter_data%update_growth: ",iter_data%update_growth(:)
!            WRITE(iout,*) "iter_data%update_growth_tot: ",iter_data%update_growth_tot(:)
!            WRITE(iout,*) "iter_data%tot_parts_old: ",iter_data%tot_parts_old(:)
!            WRITE(iout,*) "iter_data%update_iters: ",iter_data%update_iters
!            CALL neci_flush(iout)


            if(tInstGrowthRate) then
!Calculate the growth rate simply using the two points at the beginning and the
!end of the update cycle. 
                if ((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                    !COMPLEX
                    AllGrowRate = (sum(iter_data%update_growth_tot &
                               + iter_data%tot_parts_old)) &
                              / real(sum(iter_data%tot_parts_old), dp)
                else
                    do run=1,inum_runs
                        AllGrowRate(run) = (iter_data%update_growth_tot(run) &
                                   + iter_data%tot_parts_old(run)) &
                                  / real(iter_data%tot_parts_old(run), dp)
                    enddo
                endif
            else
!Instead attempt to calculate the average growth over every iteration
!over the update cycle
                if ((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                    !COMPLEX
                    AllGrowRate = (sum(AllSumWalkersCyc)/real(StepsSft,dp)) &
                                    /sum(OldAllAvWalkersCyc)
                else

                    do run=1,inum_runs
                        AllGrowRate(run) = (AllSumWalkersCyc(run)/real(StepsSft,dp)) &
                                        /OldAllAvWalkersCyc(run)
                    enddo
                endif
            endif

            ! For complex case, obtain both Re and Im parts
            if ((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                IF(iter_data%tot_parts_old(1).gt.0) THEN
                    AllGrowRateRe = (iter_data%update_growth_tot(1) + &
                                     iter_data%tot_parts_old(1)) / &
                                     iter_data%tot_parts_old(1)
                ENDIF
                IF(iter_data%tot_parts_old(lenof_sign).gt.0) THEN
                    AllGrowRateIm = (iter_data%update_growth_tot(lenof_sign) + &
                                         iter_data%tot_parts_old(lenof_sign)) / &
                                         iter_data%tot_parts_old(lenof_sign)
                ENDIF
            endif

!AJWT commented this out as DMC says it's not being used, and it gave a divide by zero
            ! Initiator abort growth rate
!            if (tTruncInitiator) then
!                AllGrowRateAbort = (sum(iter_data%update_growth_tot + &
!                                    iter_data%tot_parts_old) + AllNoAborted) &
!                                    / (sum(iter_data%tot_parts_old) &
!                                       + AllNoAbortedOld)
!            endif

            ! Exit the single particle phase if the number of walkers exceeds
            ! the value in the input file. If particle no has fallen, re-enter
            ! it.
            tReZeroShift = .false.
            do run=1,inum_runs
                if (TSinglePartPhase(run)) then
    ! AJWT dislikes doing this type of if based on a (seeminly unrelated) input option, but can't see another easy way.
    !  TODO:  Something to make it better
                    if(.not.tCCMC) then
                        tot_walkers = InitWalkers * int(nNodes,int64)
                    else
                        tot_walkers = InitWalkers
                    endif

#ifdef __CMPLX
                    if ((sum(AllTotParts) > tot_walkers) .or. &
                         (abs_sign(AllNoatHF) > MaxNoatHF)) then
    !                     WRITE(iout,*) "AllTotParts: ",AllTotParts(1),AllTotParts(2),tot_walkers
                        write (iout, '(a,i13,a)') 'Exiting the single particle growth phase on iteration: ',iter + PreviousCycles, &
                                     ' - Shift can now change'
                        VaryShiftIter = Iter
                        iBlockingIter = Iter + PreviousCycles
                        tSinglePartPhase = .false.
                        if(TargetGrowRate(1).ne.0.0_dp) then
                            write(iout,"(A)") "Setting target growth rate to 1."
                            TargetGrowRate=0.0_dp
                        endif

                        ! If enabled, jump the shift to the value preducted by the
                        ! projected energy!
                        if (tJumpShift) then
                            DiagSft = real(proje_iter,dp)
                            defer_update = .true.
                        end if
                    elseif (abs_sign(AllNoatHF) < (MaxNoatHF - HFPopThresh)) then
                        write (iout, '(a,i13,a)') 'No at HF has fallen too low - reentering the &
                                     &single particle growth phase on iteration',iter + PreviousCycles,' - particle number &
                                     &may grow again.'
                        tSinglePartPhase = .true.
                        tReZeroShift = .true.
                    endif
#else
                    start_varying_shift = .false.
                    if (tLetInitialPopDie) then
                        if (AllTotParts(run) < tot_walkers) start_varying_shift = .true.
                    else
                        if ((AllTotParts(run) > tot_walkers) .or. &
                             (abs(AllNoatHF(run)) > MaxNoatHF)) start_varying_shift = .true.
                    end if

                    if (start_varying_shift) then
    !                     WRITE(iout,*) "AllTotParts: ",AllTotParts(1),AllTotParts(2),tot_walkers
                        write (iout, '(a,i13,a,i1)') 'Exiting the single particle growth phase on iteration: ' &
                                     ,iter + PreviousCycles, ' - Shift can now change for population', run
                        VaryShiftIter(run) = Iter
                        iBlockingIter(run) = Iter + PreviousCycles
                        tSinglePartPhase(run) = .false.
                        if(TargetGrowRate(run).ne.0.0_dp) then
                            write(iout,"(A)") "Setting target growth rate to 1."
                            TargetGrowRate(run)=0.0_dp
                        endif

                        ! If enabled, jump the shift to the value preducted by the
                        ! projected energy!
                        if (tJumpShift) then
                            DiagSft(run) = real(proje_iter(run),dp)
                            defer_update(run) = .true.
                        end if
                    elseif (abs(AllNoatHF(run)) < (MaxNoatHF - HFPopThresh)) then
                        write (iout, '(a,i13,a)') 'No at HF has fallen too low - reentering the &
                                     &single particle growth phase on iteration',iter + PreviousCycles,' - particle number &
                                     &may grow again.'
                        tSinglePartPhase(run) = .true.
                        tReZeroShift(run) = .true.
                    endif
#endif
                endif
                ! How should the shift change for the entire ensemble of walkers 
                ! over all processors.
                if (((.not. tSinglePartPhase(run)).or.(TargetGrowRate(run).ne.0.0_dp)) .and.&
                    .not. defer_update(run)) then

                    !In case we want to continue growing, TargetGrowRate > 0.0_dp
                    ! New shift value
                    if(TargetGrowRate(run).ne.0.0_dp) then
                        if((lenof_sign.eq.2).and.(inum_runs.eq.1))then
                        
                            if(sum(AllTotParts).gt.TargetGrowRateWalk(1)) then
                                !Only allow targetgrowrate to kick in once we have > TargetGrowRateWalk walkers.
                                DiagSft = DiagSft - (log(AllGrowRate-TargetGrowRate) * SftDamp) / &
                                                    (Tau * StepsSft)
                            endif
                        else
                            if(AllTotParts(run).gt.TargetGrowRateWalk(run)) then
                                !Only allow targetgrowrate to kick in once we have > TargetGrowRateWalk walkers.
                                DiagSft(run) = DiagSft(run) - (log(AllGrowRate(run)-TargetGrowRate(run)) * SftDamp) / &
                                                    (Tau * StepsSft)
                            endif
                        endif
                    else
                        if(tShiftonHFPop) then
                            !Calculate the shift required to keep the HF population constant

                            AllHFGrowRate(run) = abs(AllHFCyc(run)/real(StepsSft,dp)) / abs(OldAllHFCyc(run))

                            DiagSft(run) = DiagSft(run) - (log(AllHFGrowRate(run)) * SftDamp) / &
                                                (Tau * StepsSft)
                        else
                            !"WRITE(6,*) "AllGrowRate, TargetGrowRate", AllGrowRate, TargetGrowRate
                            DiagSft(run) = DiagSft(run) - (log(AllGrowRate(run)) * SftDamp) / &
                                                (Tau * StepsSft)
                        endif
                    endif

                    if ((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                        !COMPLEX
                        DiagSftRe = DiagSftRe - (log(AllGrowRateRe-TargetGrowRate(1)) * SftDamp) / &
                                                (Tau * StepsSft)
                        DiagSftIm = DiagSftIm - (log(AllGrowRateIm-TargetGrowRate(1)) * SftDamp) / &
                                                (Tau * StepsSft)
                    endif

                    ! Update the shift averages
                    if ((iter - VaryShiftIter(run)) >= nShiftEquilSteps) then
                        if ((iter-VaryShiftIter(run)-nShiftEquilSteps) < StepsSft) &
                            write (iout, '(a,i14)') 'Beginning to average shift value on iteration: ',iter + PreviousCycles
                        VaryShiftCycles(run) = VaryShiftCycles(run) + 1
                        SumDiagSft(run) = SumDiagSft(run) + DiagSft(run)
                        AvDiagSft(run) = SumDiagSft(run) / real(VaryShiftCycles(run), dp)
                    endif

    !                ! Update DiagSftAbort for initiator algorithm
    !                if (tTruncInitiator) then
    !                    DiagSftAbort = DiagSftAbort - &
    !                              (log(real(AllGrowRateAbort-TargetGrowRate, dp)) * SftDamp) / &
    !                              (Tau * StepsSft)
    !
    !                    if (iter - VaryShiftIter >= nShiftEquilSteps) then
    !                        SumDiagSftAbort = SumDiagSftAbort + DiagSftAbort
    !                        AvDiagSftAbort = SumDiagSftAbort / &
    !                                         real(VaryShiftCycles, dp)
    !                    endif
    !                endif
                endif
                if((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                    ! Calculate the instantaneous 'shift' from the HF population
                    HFShift(run) = -1.0_dp / abs_sign(AllNoatHF) * &
                                        (abs_sign(AllNoatHF) - abs_sign(OldAllNoatHF) / &
                                      (Tau * real(StepsSft, dp)))
                    InstShift(run) = -1.0_dp / sum(AllTotParts) * &
                                ((sum(AllTotParts) - sum(AllTotPartsOld)) / &
                                 (Tau * real(StepsSft, dp)))
                 else
                    ! Calculate the instantaneous 'shift' from the HF population
                    HFShift(run) = -1.0_dp / abs(AllNoatHF(run)) * &
                                        (abs(AllNoatHF(run)) - abs(OldAllNoatHF(run)) / &
                                      (Tau * real(StepsSft, dp)))
                    InstShift(run) = -1.0_dp / AllTotParts(run) * &
                                ((AllTotParts(run) - AllTotPartsOld(run)) / &
                                 (Tau * real(StepsSft, dp)))
                 endif

                 ! When using a linear combination, the denominator is summed
                 ! directly.
                 all_sum_proje_denominator(run) = ARR_RE_OR_CPLX(AllSumNoatHF,run)
                 all_cyc_proje_denominator(run) = AllHFCyc(run)

                 ! Calculate the projected energy.
                 if((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                     if (any(AllSumNoatHF /= 0.0)) then
                         ProjectionE = (AllSumENum) / (all_sum_proje_denominator) 
                         proje_iter = (AllENumCyc) / (all_cyc_proje_denominator) 
                        AbsProjE = (AllENumCycAbs) / (all_cyc_proje_denominator)
                    endif
                 else
                     if ((AllSumNoatHF(run) /= 0.0)) then
                         ProjectionE(run) = (AllSumENum(run)) / (all_sum_proje_denominator(run)) 
                         proje_iter(run) = (AllENumCyc(run)) / (all_cyc_proje_denominator(run)) 
                        AbsProjE(run) = (AllENumCycAbs(run)) / (all_cyc_proje_denominator(run))
                    endif
                endif
                ! If we are re-zeroing the shift
                if (tReZeroShift(run)) then
                    DiagSft(run) = 0
                    VaryShiftCycles(run) = 0
                    SumDiagSft(run) = 0
                    AvDiagSft(run) = 0
                endif
            enddo

        endif ! iProcIndex == root

        ! Broadcast the shift from root to all the other processors
        call MPIBcast (tSinglePartPhase)
        call MPIBcast (VaryShiftIter)
        call MPIBcast (DiagSft)
        
        do run=1,inum_runs
            if(.not.tSinglePartPhase(run)) then
                TargetGrowRate(run)=0.0_dp
                tSearchTau=.false.
            endif
        enddo

    end subroutine update_shift 



    subroutine rezero_iter_stats_update_cycle (iter_data, tot_parts_new_all)
        
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new_all
        
        ! Zero all of the variables which accumulate for each iteration.

        IterTime = 0.0
        SumWalkersCyc(:)=0.0_dp
        Annihilated = 0
        Acceptances = 0
        NoBorn = 0
        SpawnFromSing = 0
        NoDied = 0
        ENumCyc = 0
        ENumCycAbs = 0
        HFCyc = 0.0_dp
        cyc_proje_denominator=0
        trial_numerator = 0.0_dp
        trial_denom = 0.0_dp

        ! Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld = TotWalkers
        TotPartsOld = TotParts

        ! Save the number at HF to use in the HFShift
        OldAllNoatHF = AllNoatHF
        !OldAllHFCyc is the average HF value for this update cycle
        OldAllHFCyc = AllHFCyc/real(StepsSft,dp)
        !OldAllAvWalkersCyc gives the average number of walkers per iteration in the last update cycle
      !TODO CMO: are these summed across real/complex? 
        OldAllAvWalkersCyc = AllSumWalkersCyc/real(StepsSft,dp)

        ! Also the cumulative global variables
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        AllNoAbortedOld = AllNoAborted


        ! Reset the counters
        iter_data%update_growth = 0.0_dp
        iter_data%update_iters = 0
        iter_data%tot_parts_old = tot_parts_new_all

        max_cyc_spawn = 0

    end subroutine

    subroutine calculate_new_shift_wrapper (iter_data, tot_parts_new)

        type(fcimc_iter_data) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new
        real(dp), dimension(lenof_sign) :: tot_parts_new_all

        call collate_iter_data (iter_data, tot_parts_new, tot_parts_new_all)
        call iter_diagnostics ()
        if(tRestart) return
        call population_check ()
        call update_shift (iter_data)
        if (tFCIMCStats2) then
            call write_fcimcstats2(iter_data_fciqmc)
        else
            call WriteFCIMCStats ()
        end if
        
        call rezero_iter_stats_update_cycle (iter_data, tot_parts_new_all)

    end subroutine calculate_new_shift_wrapper

    subroutine FlipSign(part_type)

        ! This routine flips the sign of all particles on the node with a
        ! given particle type. This allows us to keep multiple parallel
        ! simulations sign coherent.

        integer, intent(in) :: part_type
        character(*), parameter :: t_r = 'FlipSign'
        integer :: i
        real(dp) :: sgn
    
        do i = 1, int(TotWalkers, sizeof_int)

            sgn = extract_part_sign(CurrentDets(:,i), part_type)
            sgn = -sgn
            call encode_part_sign(CurrentDets(:,i), sgn, part_type)

            ! Flip average signs too.
            if (tFillingStochRDMOnFly) then
                if (lenof_sign /= 1) &
                    call stop_all(t_r, 'Not yet implemented')
                call set_av_sgn(i, -get_av_sgn(i))
            end if

        enddo

        ! Reverse the flag for whether the sign of the particles has been
        ! flipped so the ACF can be correctly calculated
        tFlippedSign = .not. tFlippedSign
    
    end subroutine FlipSign



    SUBROUTINE WriteFciMCStatsHeader()

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(iout,*) ""
            IF(tTruncInitiator) THEN
                WRITE(initiatorstats_unit,"(A2,A10,11A20)") "# ","1.Step","2.TotWalk","3.Annihil","4.Died", &
                & "5.Born","6.TotUniqDets",&
&               "7.InitDets","8.NonInitDets","9.InitWalks","10.NonInitWalks","11.AbortedWalks"
            ENDIF
            IF(tLogComplexPops) THEN
                WRITE(complexstats_unit,"(A)") '#   1.Step  2.Shift     3.RealShift     4.ImShift   5.TotParts      " &
                & //"6.RealTotParts      7.ImTotParts'
            ENDIF

#ifdef __CMPLX
            if(tMCOutput) then
                write(iout, '(a)') "       Step     Shift      WalkerCng(Re)  &
                       &WalkerCng(Im)    TotWalkers(Re)   TotWalkers(Im)    &
                       &Proj.E(Re)   ProjE(Im)     Proj.E.ThisCyc(Re)  &
                       &Proj.E.ThisCyc(Im)   NoatHF(Re)   NoatHF(Im)   &
                       &NoatDoubs      AccRat     UniqueDets     IterTime"
            endif
            write(fcimcstats_unit, "(a,i4,a,l1,a,l1,a,l1)") &
                   "# FCIMCStats VERSION 2 - COMPLEX : NEl=", nel, &
                   " HPHF=", tHPHF, ' Lz=', tFixLz, &
                   ' Initiator=', tTruncInitiator
            write(fcimcstats_unit, "(a)") &
                   "#     1.Step   2.Shift    3.WalkerCng(Re)  &
                   &4.WalkerCng(Im)   5.TotWalkers(Re)  6.TotWalkers(Im)  &
                   &7.Proj.E(Re)   8.Proj.E(Im)   9.Proj.E.ThisCyc(Re)  &
                   &10.Proj.E.ThisCyc(Im)  11.NoatHF(Re)   12.NoatHF(Im)  &
                   &13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime &
                   &17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime &
                   &  20.HFInstShift  21.TotInstShift  &
                   &22.HFContribtoE(Both)  &
                   &23.NumContribtoE(Re)  &
                   &24.NumContribtoE(Im)  25.HF weight   26.|Psi|    &
                   &27.Inst S^2  28.SpawnedParts  29.MergedParts  &
                   &30.Zero elems   31.PartsDiffProc   32.MaxCycSpawn"
#elif __DOUBLERUN
            write(fcimcstats_unit2, "(a,i4,a,l1,a,l1,a,l1)") &
                  "# FCIMCStats VERSION 2 - REAL : NEl=", nel, &
                  " HPHF=", tHPHF, ' Lz=', tFixLz, &
                  ' Initiator=', tTruncInitiator
            write(fcimcstats_unit2, "(A)", advance = 'no') &
                  "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     &
                  &5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  &
                  &9.Proj.E       10.Av.Shift 11.Proj.E.ThisCyc  12.NoatHF &
                  &13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime &
                  &17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  &
                  &20.ProjE.ThisIter  21.HFInstShift  22.TotInstShift  &
                  &23.Tot-Proj.E.ThisCyc   24.HFContribtoE  25.NumContribtoE &
                  &26.HF weight    27.|Psi|     28.Inst S^2 29.Inst S^2 30.AbsProjE &
                  &31.|Semistoch|/|Psi|   32.PartsDiffProc"
           if (tTrialWavefunction) then 
                  write(fcimcstats_unit2, "(A)", advance = 'no') &
                  "  33.TrialNumerator  34.TrialDenom  35.TrialOverlap"
           end if

           write(fcimcstats_unit2, "()", advance = 'yes')
#endif
#ifndef __CMPLX
            if(tMCOutput) then
                write(iout, "(A)", advance = 'no') "        Step    Shift           &
                      &WalkerCng       GrowRate        TotWalkers      Annihil         &
                      &NoDied          NoBorn          Proj.E          Av.Shift        &
                      &Proj.E.Cyc"
                if (tTrialWavefunction) write(iout, "(A)", advance = 'no') &
                      "    Trial.E.Cyc "
                write(iout, "(A)", advance = 'yes') "      NoatHF          NoatDoubs       &
                &AccRat        UniqueDets   IterTime"
            endif
            write(fcimcstats_unit, "(a,i4,a,l1,a,l1,a,l1)") &
                  "# FCIMCStats VERSION 2 - REAL : NEl=", nel, &
                  " HPHF=", tHPHF, ' Lz=', tFixLz, &
                  ' Initiator=', tTruncInitiator
            write(fcimcstats_unit, "(A)", advance = 'no') &
                  "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     &
                  &5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  &
                  &9.Proj.E       10.Av.Shift 11.Proj.E.ThisCyc  12.NoatHF &
                  &13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime &
                  &17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  &
                  &20.ProjE.ThisIter  21.HFInstShift  22.TotInstShift  &
                  &23.Tot-Proj.E.ThisCyc   24.HFContribtoE  25.NumContribtoE &
                  &26.HF weight    27.|Psi|     28.Inst S^2 &
                  &29.Inst S^2   30.AbsProjE   31.PartsDiffProc &
                  &32.|Semistoch|/|Psi|  33.MaxCycSpawn"
           if (tTrialWavefunction) then 
                  write(fcimcstats_unit, "(A)", advance = 'no') &
                  "  34.TrialNumerator  35.TrialDenom  36.TrialOverlap"
           end if

           write(fcimcstats_unit, "()", advance = 'yes')

#endif
            
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

    subroutine WriteFCIMCStats()
        INTEGER :: i, run
        real(dp),dimension(inum_runs) :: FracFromSing

        ! What is the current value of S2
        if (tCalcInstantS2) then
            if (mod(iter / StepsSft, instant_s2_multiplier) == 0) then
                if (tSpatialOnlyhash) then
                    curr_S2 = calc_s_squared (.false.)
                else
                    curr_S2 = calc_s_squared_star (.false.)
                end if
            end if
        else
            curr_S2 = -1
        end if

        ! What is the current value of S2 considering only initiators
        if (tCalcInstantS2Init) then
            if (mod(iter / StepsSft, instant_s2_multiplier_init) == 0) then
                if (tSpatialOnlyhash) then
                    curr_S2_init = calc_s_squared (.true.)
                else
                    curr_S2_init = calc_s_squared_star (.true.)
                end if
            end if
        else
            curr_S2_init = -1
        endif

        !To prevent /0 problems
        do run=1,inum_runs
            if(AllNoBorn(run).ne.0) then
                FracFromSing(run)=real(AllSpawnFromSing(run),dp) / real(AllNoBorn(run),dp)
            else
                FracFromSing(run)=0.0_dp
            endif
        enddo

        if (iProcIndex == root) then

        ! Becuase tot_trial_numerator/tot_trial_denom is the energy relative to the the trial
        ! energy, add on this contribution to make it relative to the HF energy.
        if (tTrialWavefunction) then
            tot_trial_numerator = tot_trial_numerator + (tot_trial_denom*(trial_energy-Hii))
        end if

#ifdef __CMPLX
            write(fcimcstats_unit,"(I12,5G16.7,7G17.9,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,8G17.9,I13,&

                                  &g16.7)") &
                Iter + PreviousCycles, &                !1.
                DiagSft, &                              !2.
                AllTotParts(1) - AllTotPartsOld(1), &   !3.
                AllTotParts(2) - AllTotPartsOld(2), &   !4.
                AllTotParts(1), &                       !5.
                AllTotParts(2), &                       !6.
                real(ProjectionE, dp), &                !7.     real \sum[ nj H0j / n0 ]
                aimag(projectionE), &                   !8.     Im   \sum[ nj H0j / n0 ]
                real(proje_iter, dp), &                 !9.     
                aimag(proje_iter), &                    !10.
                AllNoatHF(1), &                         !11.
                AllNoatHF(2), &                         !12.
                AllNoatDoubs, &                         !13.
                AccRat, &                               !14.
                AllTotWalkers, &                        !15.
                IterTime, &                             !16.
                FracFromSing(1), &                      !17.
                WalkersDiffProc, &                           !18.
                TotImagTime, &                               !19.
                HFShift, &                                   !20.
                InstShift, &                                 !21.
                real((AllHFCyc*conjg(AllHFCyc)),dp), &     !22 |n0|^2  denominator for both calcs
                real((AllENumCyc*conjg(AllHFCyc)),dp), &   !23. Re[\sum njH0j]xRe[n0]+Im[\sum njH0j]xIm[n0]   No div by StepsSft
                aimag(AllENumCyc*conjg(AllHFCyc)), &       !24.Im[\sum njH0j]xRe[n0]-Re[\sum njH0j]xIm[n0]   since no physicality
                sqrt(sum(AllNoatHF**2)) / norm_psi, & !25
                norm_psi, &                           !26
                curr_S2, &                            !27
                PartsDiffProc, &                      !28
                all_max_cyc_spawn                     !29.

            if(tMCOutput) then
                write (iout, "(I12,13G16.7,I12,G13.5)") &
                    Iter + PreviousCycles, &
                    DiagSft, &
                    AllTotParts(1) - AllTotPartsOld(1), &
                    AllTotParts(2) - AllTotPartsOld(2), &
                    AllTotParts(1), &
                    AllTotParts(2), &
                    real(ProjectionE, dp), &
                    aimag(ProjectionE), &
                    real(proje_iter, dp), &
                    aimag(proje_iter), &
                    AllNoatHF(1), &
                    AllNoatHF(2), &
                    AllNoatDoubs, &
                    AccRat, &
                    AllTotWalkers, &
                    IterTime
            endif
            if (tTruncInitiator) then
               write(initiatorstats_unit,"(I12,4G16.7,3I20,4G16.7)")&
                   Iter + PreviousCycles, sum(AllTotParts), &
                   AllAnnihilated(1), AllNoDied(1), AllNoBorn(1), AllTotWalkers,&
                   AllNoInitDets(1), AllNoNonInitDets(1), AllNoInitWalk(1), &
                   AllNoNonInitWalk(1),AllNoAborted(1), AllNoRemoved(1)
            endif
            if (tLogComplexPops) then
                write (complexstats_unit,"(I12,6G16.7)") &
                    Iter + PreviousCycles, DiagSft, DiagSftRe, DiagSftIm, &
                    sum(AllTotParts), AllTotParts(1), AllTotParts(lenof_sign)
            endif
#elif __DOUBLERUN
            write(fcimcstats_unit2,"(i12,7g16.7,5g17.9,g13.5,i12,g13.5,g17.5,&
                                   &i13,g13.5,11g17.9,i13,2g16.7)",advance = 'no') &
                Iter + PreviousCycles, &                   ! 1.
                DiagSft(2), &                              ! 2.
                AllTotParts(2) - AllTotPartsOld(2), &      ! 3.
                AllGrowRate(2), &                          ! 4.
                AllTotParts(2), &                          ! 5.
                AllAnnihilated(2), &                       ! 6.
                AllNoDied(2), &                            ! 7.
                AllNoBorn(2), &                            ! 8.
                ProjectionE(2), &                          ! 9.
                AvDiagSft(2), &                            ! 10.
                proje_iter(2), &                           ! 11.
                AllNoatHF(2), &                            ! 12.
                AllNoatDoubs(2), &                         ! 13.
                AccRat(2), &                               ! 14.
                AllTotWalkers, &                           ! 15.
                IterTime, &                                ! 16.
                FracFromSing(2), &                         ! 17.
                WalkersDiffProc, &                         ! 18.
                TotImagTime, &                             ! 19.
                0.0_dp, &                                  ! 20.
                HFShift(2), &                              ! 21.
                InstShift(2), &                            ! 22.
                proje_iter(2) + Hii, &                     ! 23.
                (AllHFCyc(2) / StepsSft), &                ! 24.
                (AllENumCyc(2) / StepsSft), &              ! 25.
                AllNoatHF(2) / norm_psi(2), &              ! 26.
                norm_psi(2), &                             ! 27.
                curr_S2(2), curr_S2_init(2), &             ! 28, 29.
                AbsProjE(2), &                             ! 30.
                PartsDiffProc, &                           ! 31.
                norm_semistoch(2)/norm_psi(2), &           ! 32.
                all_max_cyc_spawn                          ! 33.
                if (tTrialWavefunction) then
                    write(fcimcstats_unit2, "(3G16.7)", advance = 'no') &
                    (tot_trial_numerator(2) / StepsSft), &
                    (tot_trial_denom(2) / StepsSft), &
                    abs(tot_trial_denom(2) / (norm_psi(2)*StepsSft))
                end if
                
                write(fcimcstats_unit2, "()", advance = 'yes')
#endif
#ifndef __CMPLX

            write(fcimcstats_unit,"(i12,7g16.7,5g17.9,g13.5,i12,g13.5,g17.5,&
                                  &i13,g13.5,11g17.9,i13,2g16.7)",advance = 'no') &
                Iter + PreviousCycles, &                   ! 1.
                DiagSft(1), &                              ! 2.
                AllTotParts(1) - AllTotPartsOld(1), &      ! 3.
                AllGrowRate(1), &                          ! 4.
                AllTotParts(1), &                          ! 5.
                AllAnnihilated(1), &                       ! 6.
                AllNoDied(1), &                            ! 7.
                AllNoBorn(1), &                            ! 8.
                ProjectionE(1), &                          ! 9.
                AvDiagSft(1), &                            ! 10.
                proje_iter(1), &                           ! 11.
                AllNoatHF(1), &                            ! 12.
                AllNoatDoubs(1), &                         ! 13.
                AccRat(1), &                               ! 14.
                AllTotWalkers, &                           ! 15.
                IterTime, &                                ! 16.
                FracFromSing(1), &                         ! 17.
                WalkersDiffProc, &                         ! 18.
                TotImagTime, &                             ! 19.
                0.0_dp, &                                  ! 20.
                HFShift(1), &                              ! 21.
                InstShift(1), &                            ! 22.
                proje_iter(1) + Hii, &                     ! 23.
                (AllHFCyc(1) / StepsSft), &                ! 24.
                (AllENumCyc(1) / StepsSft), &              ! 25.
                AllNoatHF(1) / norm_psi(1), &              ! 26.
                norm_psi(1), &                             ! 27.
                curr_S2(1), curr_S2_init(1), &             ! 28, 29.
                AbsProjE(1), &                             ! 30.
                PartsDiffProc, &                           ! 31.
                norm_semistoch(1)/norm_psi(1), &           ! 32.
                all_max_cyc_spawn                          ! 33.
                if (tTrialWavefunction) then
                    write(fcimcstats_unit, "(3g16.7)", advance = 'no') &
                    (tot_trial_numerator / StepsSft), &             ! 34.
                    (tot_trial_denom / StepsSft), &                 ! 35.
                    abs((tot_trial_denom / (norm_psi*StepsSft)))    ! 36.
                end if
                write(fcimcstats_unit, "()", advance = 'yes')

            if(tMCOutput) then
                write (iout, "(I12,10G16.7)", advance = 'no') &
                    Iter + PreviousCycles, &
                    DiagSft(1), &
                    AllTotParts(1) - AllTotPartsOld(1), &
                    AllGrowRate(1), &
                    AllTotParts(1), &
                    AllAnnihilated(1), &
                    AllNoDied(1), &
                    AllNoBorn(1), &
                    ProjectionE(1), &
                    AvDiagSft(1), &
                    proje_iter(1)
                if (tTrialWavefunction) write(iout, "(G16.7)", advance = 'no') &
                    (tot_trial_numerator(1)/tot_trial_denom(1))
                write (iout, "(3G16.7,I12,G13.5)", advance = 'yes') &
                    AllNoatHF(1), &
                    AllNoatDoubs(1), &
                    AccRat(1), &
                    AllTotWalkers, &
                    IterTime
            endif
            if (tTruncInitiator) then
               write(initiatorstats_unit,"(I12,4G16.7,3I20,4G16.7)")&
                   Iter + PreviousCycles, AllTotParts(1), &
                   AllAnnihilated(1), AllNoDied(1), AllNoBorn(1), AllTotWalkers,&
                   AllNoInitDets(1), AllNoNonInitDets(1), AllNoInitWalk(1), &
                   AllNoNonInitWalk(1),AllNoAborted(1), AllNoRemoved(1)
            endif
#endif


            if(tMCOutput) then
                call neci_flush(iout)
            endif
            call neci_flush(fcimcstats_unit)
            if (inum_runs.eq.2) call neci_flush(fcimcstats_unit2)
            
        endif

    end subroutine WriteFCIMCStats


    subroutine open_create_fciqmc_stats(funit)

        integer, intent(in) :: funit
        character(*), parameter :: t_r = 'open_create_fciqmc_stats'

        character(30) :: filename
        character(43) :: filename2
        character(12) :: num
        integer :: i, ierr
        logical :: exists

        ! If we are using Molpro, then append the molpro ID to uniquely
        ! identify the output
        if (tMolpro .and. .not. tMolproMimic) then
            filename = 'fciqmc_stats_' // adjustl(MolproID)
        else
            filename = 'fciqmc_stats'
        end if
        

        if (tReadPops) then

            ! If we are reading from a POPSFILE, then we want to continue an
            ! existing fciqmc_stats file if it exists.
            open(funit, file=filename, status='unknown', position='append')

        else

            ! If we are doing a normal calculation, move existing fciqmc_stats
            ! files so that they are not overwritten, and then create a new one
            inquire(file=filename, exist=exists)
            if (exists) then

                ! Loop until we find an available spot to move the existing
                ! file to.
                i = 1
                do while(exists)
                    write(num, '(i12)') i
                    filename2 = trim(adjustl(filename)) // "." // &
                                trim(adjustl(num))
                    inquire(file=filename2, exist=exists)
                    if (i > 10000) &
                        call stop_all(t_r, 'Error finding free fciqmc_stats.*')
                    i = i + 1
                end do

                ! Move the file
                call rename(filename, filename2)

            end if

            ! And finally open the file
            open(funit, file=filename, status='unknown', iostat=ierr)

        end if

    end subroutine



    subroutine write_fcimcstats2(iter_data, initial)

        ! Write output to our FCIMCStats file.
        ! This should be done in a nice general way, so that we make merges
        ! somewhat easier.
        !
        ! --> The column numbers are no longer as well fixed (oh well...)

        type(fcimc_iter_data), intent(in) :: iter_data
        logical, intent(in), optional :: initial

        ! Use a state type to keep things compact and tidy below.
        type(write_state_t), save :: state
        logical, save :: inited = .false.
        character(5) :: tmpc
        integer :: p
        logical :: init

        ! Provide default 'initial' option
        if (present(initial)) then
            state%init = initial
        else
            state%init = .false.
        end if

        ! If the output file hasn't been opened yet, then create it.
        if (iProcIndex == Root .and. .not. inited) then
            state%funit = get_free_unit()
            call open_create_fciqmc_stats(state%funit)
            inited = .true.
        end if

        ! ------------------------------------------------
        ! This is where any calculation that needs multiple nodes should go
        ! ------------------------------------------------
        ! ------------------------------------------------
    
        if (iProcIndex == root) then

            ! Only do the actual outputting on the head node.

            ! Don't treat the header line as data. Add padding to align the
            ! other columns. We also add a # to the first line of data, so
            ! that there aren't repeats if starting from POPSFILES
            if (state%init .or. state%prepend) then
                write(state%funit, '("#")', advance='no')
                if (tMCOutput) write(iout, '("#")', advance='no')
                state%prepend = state%init
            else if (.not. state%prepend) then
                write(state%funit, '(" ")', advance='no')
                if (tMCOutput) write(iout, '(" ")', advance='no')
            end if

            ! And output the actual data!
            state%cols = 0
            state%cols_mc = 0
            state%mc_out = tMCOutput
            call stats_out(state,.true., iter, 'Iter.')
            call stats_out(state,.true., sum(abs(AllTotParts)), 'Tot. parts')
            call stats_out(state,.true., sum(abs(AllNoatHF)), 'Tot. ref')
#ifdef __CMPLX
            call stats_out(state,.true., real(proje_iter(1)), 'Re Proj. E')
            call stats_out(state,.true., aimag(proje_iter(1)), 'Im Proj. E')
#else
            call stats_out(state,.true., proje_iter(1), 'Proj. E (cyc)')
#endif
            call stats_out(state,.true., DiagSft(1), 'Shift. (cyc)')
            call stats_out(state,.true., IterTime, 'Iter. time')
            call stats_out(state,.false., AllNoBorn(1), 'No. born')
            call stats_out(state,.false., AllNoDied(1), 'No. died')
            call stats_out(state,.false., AllAnnihilated(1), 'No. annihil')
            call stats_out(state,.false., AllGrowRate(1), 'Growth fac.')
            call stats_out(state,.false., AccRat(1), 'Acc. rate')
            call stats_out(state,.false., TotImagTime, 'Im. time')
#ifdef __CMPLX
            call stats_out(state,.true., real(proje_iter(1)) + Hii, &
                           'Tot. Proj. E')
            call stats_out(state,.true., aimag(proje_iter(1)) + Hii, &
                           'Tot. Proj. E')
#else
            call stats_out(state,.true., proje_iter(1) + Hii, 'Tot. Proj. E')
#endif

            ! If we are running multiple (replica) simulations, then we
            ! want to record the details of each of these
#ifdef __PROG_LENOFSIGN
            do p = 1, lenof_sign
                write(tmpc, '(i5)') p
                call stats_out (state, .false., AllTotParts(p), &
                                'Parts (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., AllNoatHF(p), &
                                'Ref (' // trim(adjustl(tmpc)) // ")")
            end do
#endif

            ! Put the conditional columns at the end, so that the column
            ! numbers of the data are as stable as reasonably possible (for
            ! people who want to use gnuplot/not analyse column headers too
            ! frequently).
            ! This also makes column contiguity on resumes as likely as
            ! possible.
            if (tTruncInitiator) &
                call stats_out(state,.false., AllNoAborted(1), 'No. aborted')

            ! And we are done
            write(state%funit, *)
            if (tMCOutput) write(iout, *)
            call neci_flush(state%funit)
            call neci_flush(iout)

        end if

    end subroutine write_fcimcstats2



    !Ensure that the new FCIMCStats file which is about to be opened does not overwrite any other FCIMCStats
    !files. If there is already an FCIMCStats file present, then move it to FCIMCStats.x, where x is a largest
    !free filename.
    subroutine MoveFCIMCStatsFiles()
#ifdef NAGF95
        USe f90_unix_dir, only: rename
#endif
        integer :: extension,stat
        logical :: exists
        character(len=22) :: abstr
!        character(len=36) :: command
        character(len=*), parameter :: t_r='MoveFCIMCStatsFiles'

        if(tMolpro) then
            inquire(file='FCIQMCStats',exist=exists)
        else
            inquire(file='FCIMCStats',exist=exists)
        endif
        if(exists) then
            !We already have an FCIMCStats file - move it to the end of the list of FCIMCStats files.

            extension=1
            do while(.true.)
                abstr=''
                write(abstr,'(I12)') extension
                if(tMolpro) then
                    abstr='FCIQMCStats.'//adjustl(abstr)
                else
                    abstr='FCIMCStats.'//adjustl(abstr)
                endif
                inquire(file=abstr,exist=exists)
                if(.not.exists) exit
                extension=extension+1
                if(extension.gt.10000) then
                    call stop_all(t_r,"Error finding free FCIMCStats name")
                endif
            enddo
            
            !We have got a unique filename
            !Do not use system call
!            command = 'mv' // ' FCIMCStats ' // abstr
!            stat = neci_system(trim(command))

            if(tMolpro) then
                call rename('FCIQMCStats',abstr)
            else
                call rename('FCIMCStats',abstr)
            endif
            !Doesn't like the stat argument
!            if(stat.ne.0) then
!                call stop_all(t_r,"Error with renaming FCIMCStats file")
!            endif
        endif

    end subroutine MoveFCIMCStatsFiles


    subroutine check_start_rdm()
! This routine checks if we should start filling the RDMs - and does so if we should.        
        use nElRDMMod , only : DeAlloc_Alloc_SpawnedParts
        use LoggingData, only: tReadRDMs
        implicit none
        logical :: tFullVaryshift
        integer :: iunit_4

        tFullVaryShift=.false.

        if (.not. tSinglePartPhase(1).and.(.not.tSinglePartPhase(inum_runs))) tFullVaryShift=.true.

        !If we're reading in the RDMs we've already started accumulating them in a previous calculation
        ! We don't want to put in an arbitrary break now!
        if(tReadRDMs)   IterRDMonFly=0

        IF(tFullVaryShift .and. ((Iter - maxval(VaryShiftIter)).eq.(IterRDMonFly+1))) THEN
        ! IterRDMonFly is the number of iterations after the shift has changed that we want 
        ! to fill the RDMs.  If this many iterations have passed, start accumulating the RDMs! 
        
            IterRDMStart = Iter+PreviousCycles
            IterRDM_HF = Iter+PreviousCycles

            !if(tReadRDMs .and. tReadRDMAvPop) then
                !We need to read in the values of IterRDMStart and IterRDM_HF
            !    iunit_4=get_free_unit()
            !    OPEN(iunit_4,FILE='ITERRDMSTART',status='old')
            !    read(iunit_4, *) IterRDMStart, IterRDM_HF, AvNoAtHF

            !endif

            !We have reached the iteration where we want to start filling the RDM.
            if(tExplicitAllRDM) then
                ! Explicitly calculating all connections - expensive...
                if(inum_runs.eq.2) call stop_all('check_start_rdm',"Cannot yet do replica RDM sampling with explicit RDMs. &
                    & e.g Hacky bit in Gen_Hist_ExcDjs to make it compile")
                
                tFillingExplicRDMonFly = .true.
                if(tHistSpawn) NHistEquilSteps = Iter
            else
                extract_bit_rep_avsign => extract_bit_rep_avsign_norm
                !By default - we will do a stochastic calculation of the RDM.
                tFillingStochRDMonFly = .true.
                !if(.not.tHF_Ref_Explicit) call DeAlloc_Alloc_SpawnedParts()
                call DeAlloc_Alloc_SpawnedParts()
                !The SpawnedParts array now needs to carry both the spawned parts Dj, and also it's 
                !parent Di (and it's sign, Ci). - We deallocate it and reallocate it with the larger size.
                !Don't need any of this if we're just doing HF_Ref_Explicit calculation.
                !This is all done in the add_rdm_hfconnections routine.
            endif
            if(RDMExcitLevel.eq.1) then
                WRITE(6,'(A)') 'Calculating the 1 electron density matrix on the fly.'
            else
                WRITE(6,'(A)') 'Calculating the 2 electron density matrix on the fly.'
            endif
            WRITE(6,'(A,I10)') 'Beginning to fill the RDMs during iteration',Iter
        ENDIF

    end subroutine check_start_rdm


    LOGICAL FUNCTION TestifDETinCAS(CASDet)
        INTEGER :: k,z,CASDet(NEl), orb
        LOGICAL :: tElecInVirt, bIsCsf

        ! CASmax is the max spin orbital number (when ordered energetically) 
        ! within the chosen active space. Spin orbitals with energies larger 
        ! than this maximum value must be unoccupied for the determinant to 
        ! be in the active space.
!        CASmax=NEl+VirtCASorbs

        ! CASmin is the max spin orbital number below the active space.  As 
        ! well as the above criteria, spin orbitals with energies equal to,
        ! or below that of the CASmin orbital must be completely occupied for
        ! the determinant to be in the active space.
        !
        ! (These have been moved to the InitCalc subroutine so they're not
        !  calculated each time).
!        CASmin=NEl-OccCASorbs

        bIsCsf = iscsf(CASDet)

        z=0
        tElecInVirt=.false.
        do k=1,NEl      ! running over all electrons
            ! TODO: is it reasonable to just apply the orbital mask anyway?
            !       it is probably faster than running iscsf...
            if (bIsCsf) then
                orb = iand(CASDet(k), csf_orbital_mask)
            else
                orb = CASDet(k)
            endif

            if (SpinInvBRR(orb).gt.CASmax) THEN
                tElecInVirt=.true.
                EXIT            
                ! if at any stage an electron has an energy greater than the 
                ! CASmax value, the determinant can be ruled out of the active
                ! space.  Upon identifying this, it is not necessary to check
                ! the remaining electrons.
            else
                if (SpinInvBRR(orb).le.CASmin) THEN
                    z=z+1
                endif
                ! while running over all electrons, the number that occupy 
                ! orbitals equal to or below the CASmin cutoff are counted.
            endif
        enddo

        if(tElecInVirt.or.(z.ne.CASmin)) THEN
            ! if an electron is in an orbital above the active space, or the 
            ! inactive orbitals are not full, the determinant is automatically
            ! ruled out.
            TestifDETinCAS=.false.
        else
            ! no orbital in virtual and the inactive orbitals are completely 
            ! full - det in active space.        
            TestifDETinCAS=.true.
        endif
        
        RETURN

    END FUNCTION TestifDETinCAS




    ! This is the same as BinSearchParts1, but this time, the list to search 
    ! is passed in as an argument. The list goes from 1 to Length, but only 
    ! between MinInd and MaxInd is actually searched.
    SUBROUTINE BinSearchParts3(iLut,List,Length,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER :: Length
        INTEGER(KIND=n_int) :: iLut(0:NIfTot), List(0:NIfTot,Length)
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

!        WRITE(iout,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL neci_flush(iout)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(List(:,MaxInd),iLut(:),NIfDBO)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(iout,*) i,j,n

            ! Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it 
            ! is more or 0 if they are the same
            Comp=DetBitLT(List(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
                ! The value of the determinant at N is LESS than the 
                ! determinant we're looking for. Therefore, move the lower 
                ! bound of the search up to N. However, if the lower bound is
                ! already equal to N then the two bounds are consecutive and 
                ! we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
                    ! This deals with the case where we are interested in the
                    ! final/first entry in the list. Check the final entry of
                    ! the list and leave We need to check the last index.
                    Comp=DetBitLT(List(:,i+1),iLut(:),NIfDBO)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
                ! The value of the determinant at N is MORE than the 
                ! determinant we're looking for. Move the upper bound of the 
                ! search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

        ! If we have failed, then we want to find the index that is one less 
        ! than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts3
    

    SUBROUTINE CreateSpinInvBRR()
    ! Create an SpinInvBRR containing spin orbitals, 
    ! unlike 'createInvBRR' which only has spatial orbitals.
    ! This is used for the FixCASshift option in establishing whether or not
    ! a determinant is in the complete active space.
    ! In:
    !    BRR(i)=j: orbital i is the j-th lowest in energy.
    !    nBasis: size of basis
    ! SpinInvBRR is the inverse of BRR.  SpinInvBRR(j)=i: the j-th lowest energy
    ! orbital corresponds to the i-th orbital in the original basis.
    ! i.e the position in SpinInvBRR now corresponds to the orbital number and 
    ! the value to the relative energy of this orbital. 
    
        IMPLICIT NONE
        INTEGER :: I,t,ierr
        CHARACTER(len=*), PARAMETER :: this_routine='CreateSpinInvBrr'

        IF(ALLOCATED(SpinInvBRR)) RETURN
            
        ALLOCATE(SpinInvBRR(NBASIS),STAT=ierr)
        CALL LogMemAlloc('SpinInvBRR',NBASIS,4,this_routine,SpinInvBRRTag,ierr)
            

!        IF(iProcIndex.eq.root) THEN
!            WRITE(iout,*) "================================"
!            WRITE(iout,*) "BRR is "
!            WRITE(iout,*) BRR(:)
!        ENDIF
        
        SpinInvBRR(:)=0
        
        t=0
        DO I=1,NBASIS
            t=t+1
            SpinInvBRR(BRR(I))=t
        ENDDO

!        IF(iProcIndex.eq.root) THEN
!            WRITE(iout,*) "================================"
!            WRITE(iout,*) "SpinInvBRR is "
!            WRITE(iout,*) SpinInvBRR(:)
!        ENDIF
        
        RETURN
        
    END SUBROUTINE CreateSpinInvBRR

!This will store all the double excitations.
    SUBROUTINE StoreDoubs()
        use SystemData , only : tUseBrillouin
        use SymExcit3 , only : CountExcitations3,GenExcitations3
        INTEGER :: nJ(NEl),ierr,VecSlot,nSingles,ExcitMat3(2,2)
        LOGICAL :: tAllExcitFound,tParity

        IF(tUseBrillouin) THEN
            CALL Stop_All("StoreDoubs","Cannot have Brillouin theorem as now storing singles too...")
        ENDIF
        
!NoDoubs here is actually the singles + doubles of HF
        exflag=3
        CALL CountExcitations3(HFDet,exflag,nSingles,NoDoubs)
        NoDoubs=nSingles+NoDoubs

        ALLOCATE(DoublesDets(NEl,NoDoubs),stat=ierr)
        CALL LogMemAlloc('DoublesDets',NoDoubs*NEl,4,"StoreDoubs",DoublesDetsTag,ierr)
        DoublesDets(1:NEl,1:NoDoubs)=0
        
        VecSlot=1           !This is the next free slot in the DoublesDets array

        tAllExcitFound=.false.
        ExcitMat3(:,:)=0
!An exflag of anything but 1 or 2 indicates both the single and double excitations should be found.            
        exflag=3

        do while (.not.tAllExcitFound)
            CALL GenExcitations3(HFDet,iLutHF,nJ,exflag,ExcitMat3,tParity,tAllExcitFound,.false.)
            IF(tAllExcitFound) EXIT
            DoublesDets(1:NEl,VecSlot)=nJ(:)
            VecSlot=VecSlot+1
        enddo

!This means that now NoDoubs is double excitations AND singles
!        NoDoubs=VecSlot-1

        IF(VecSlot.ne.(NoDoubs+1)) THEN
            WRITE(iout,*) VecSlot,NoDoubs
            CALL Stop_All("StoreDoubs","Problem enumerating all double excitations")
        ENDIF

    END SUBROUTINE StoreDoubs



    subroutine test_routine()

        use bit_reps, only: add_ilut_lists

        integer(n_int) :: list_1(0:NIfTot, 10)
        integer(n_int) :: list_2(0:NIfTot, 12)
        integer(n_int) :: list_out(0:NIfTot, 11)

        integer :: i
        integer :: ndets_out
        integer(n_int) :: dets_1(6), dets_2(6)
        real(dp) :: signs_1(6), signs_2(6)
        real(dp) :: real_sign(lenof_sign)

        dets_1 =  (/28,   47,  1,   23,   32,  57/)
        signs_1 = (/-1.0, 2.0, 1.3, 4.0, 1.0, 7.0/)

        dets_2 =  (/47,  23,  32,  38,  57,  63/)
        signs_2 = (/7.0, 4.0, 1.0, 2.0, 1.2, 2.4/)

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



    subroutine CalcUEGMP2()
        use SymExcitDataMod, only: kPointToBasisFn
        use SystemData, only: ElecPairs,NMAXX,NMAXY,NMAXZ,OrbECutOff,tGCutoff,GCutoff, &
                                tMP2UEGRestrict,kiRestrict,kiMsRestrict,kjRestrict,kjMsRestrict, &
                                Madelung,tMadelung,tUEGFreeze,FreezeCutoff, kvec, tUEG2
        use Determinants, only: GetH0Element4, get_helement_excit
        integer :: Ki(3),Kj(3),Ka(3),LowLoop,HighLoop,X,i,Elec1Ind,Elec2Ind,K,Orbi,Orbj
        integer :: iSpn,FirstA,nJ(NEl),a_loc,Ex(2,2),kx,ky,kz,OrbB,FirstB
        integer :: ki2,kj2
        logical :: tParity,tMom
        real(dp) :: Ranger,mp2,mp2all,length,length_g,length_g_2
        HElement_t :: hel,H0tmp

        !Divvy up the ij pairs
        Ranger=real(ElecPairs,dp)/real(nProcessors,dp)
        LowLoop=int(iProcIndex*Ranger)+1
        Highloop=int((iProcIndex+1)*Ranger)

        if((iProcIndex+1).eq.nProcessors) Highloop=ElecPairs
        if(iProcIndex.eq.0) then
            if(lowLoop.ne.1) write(iout,*) "Error here!"
        endif
        write(iout,*) "Total ij pairs: ",ElecPairs
        write(iout,*) "Considering ij pairs from: ",LowLoop," to ",HighLoop  
!        write(iout,*) "HFDet: ",HFDet(:)

        do i=LowLoop,HighLoop   !Looping over electron pairs on this processor

            X=ElecPairs-i
            K=INT((SQRT(8.0_dp*REAL(X,dp)+1.0_dp)-1.0_dp)/2.0_dp)
            Elec1Ind=NEl-1-K
            Elec2Ind=NEl-X+((K*(K+1))/2)
            Orbi=HFDet(Elec1Ind)
            Orbj=HFDet(Elec2Ind)
            Ki=G1(Orbi)%k
            Kj=G1(Orbj)%k
            !=======================================
            if (tUEG2) then
                Ki=kvec(Orbi, 1:3)
                Kj=kvec(Orbj, 1:3)
            end if
            !=======================================
            if (tUEGFreeze) then
                ki2=ki(1)**2+ki(2)**2+ki(3)**2
                kj2=kj(1)**2+kj(2)**2+kj(3)**2
                if (.not.(ki2.gt.FreezeCutoff.and.kj2.gt.FreezeCutoff)) cycle
            endif
            if (tMP2UEGRestrict) then
                if (.not. ( &
                ( kiRestrict(1).eq.ki(1).and.kiRestrict(2).eq.ki(2).and.kiRestrict(3).eq.ki(3) .and. & 
                kjRestrict(1).eq.kj(1).and.kjRestrict(2).eq.kj(2).and.kjRestrict(3).eq.kj(3) .and. & 
                kjMsRestrict.eq.G1(Orbi)%Ms.and.kiMsRestrict.eq.G1(Orbj)%Ms ) .or. &
                ! the other way round
                ( kiRestrict(1).eq.kj(1).and.kiRestrict(2).eq.kj(2).and.kiRestrict(3).eq.kj(3) .and. & 
                kjRestrict(1).eq.ki(1).and.kjRestrict(2).eq.ki(2).and.kjRestrict(3).eq.ki(3) .and. & 
                kiMsRestrict.eq.G1(Orbi)%Ms.and.kjMsRestrict.eq.G1(Orbj)%Ms ) ) &
                ) cycle
                write(iout,*) "Restricting calculation to i,j pair: ",Orbi,Orbj
            endif

            IF((G1(Orbi)%Ms)*(G1(Orbj)%Ms).eq.-1) THEN
!We have an alpha beta pair of electrons.
                iSpn=2
            ELSE
                IF(G1(Orbi)%Ms.eq.1) THEN
!We have an alpha alpha pair of electrons.
                    iSpn=3
                ELSE
!We have a beta beta pair of electrons.
                    iSpn=1
                ENDIF
            ENDIF

!            write(iout,*) "ijpair: ",Orbi,Orbj

            if((iSpn.eq.3).or.(iSpn.eq.1)) then
                if(iSpn.eq.3) then
                    FirstA=2    !Loop over alpha
                else
                    FirstA=1    !Loop over beta
                endif

                do a_loc=FirstA,nBasis,2
                    !Loop over all a

                    !Reject if a is occupied
                    if(IsOcc(iLutHF,a_loc)) cycle

                    Ka=G1(a_loc)%k
                    !=======================================
                    if (tUEG2) then
                        Ka=kvec(a_loc, 1:3)
                    end if
                    !=======================================

                    !Find k labels of b
                    kx=Ki(1)+Kj(1)-Ka(1)
                    if(abs(kx).gt.NMAXX) cycle
                    ky=Ki(2)+Kj(2)-Ka(2)
                    if(abs(ky).gt.NMAXY) cycle
                    kz=Ki(3)+Kj(3)-Ka(3)
                    if(abs(kz).gt.NMAXZ) cycle
                    !if(tGCutoff) then
                    !    length_g=real((kx-kj(1))**2+(ky-kj(2))**2+(kz-kj(3))**2)
                    !    length_g_2=real((kx-ki(1))**2+(ky-ki(2))**2+(kz-ki(3))**2)
                    !    if(length_g.gt.gCutoff.and.length_g_2.gt.gCutoff) cycle
                    !endif
                    length=real((kx**2)+(ky**2)+(kz**2),dp)
                    if(length.gt.OrbECutoff) cycle

                    !Find the actual k orbital
                    if(iSpn.eq.3) then
                        !want alpha
                        OrbB=kPointToBasisFn(kx,ky,kz,2)
                    else
                        !want beta
                        OrbB=kPointToBasisFn(kx,ky,kz,1)
                    endif

                    !Reject k orbital if it is occupied or gt a
                    if(IsOcc(iLutHF,OrbB)) cycle
                    if(OrbB.ge.a_loc) cycle

                    !Find det
!                    write(iout,*) "OrbB: ",OrbB
                    call make_double (HFDet, nJ, elec1ind, elec2ind, a_loc, &
                                      orbB, ex, tParity)
                    !Sum in mp2 contrib
                    hel=get_helement_excit(HFDet,nJ,2,Ex,tParity)

                    H0tmp=getH0Element4(nJ,HFDet)
                    H0tmp=Fii-H0tmp
                    if(tMadelung) then
                        H0tmp=H0tmp+2.0_dp*Madelung
                    endif
                    mp2=mp2+(hel**2)/H0tmp
!                    write(iout,*) (hel**2),H0tmp
                enddo

            elseif(iSpn.eq.2) then
                do a_loc=1,nBasis
                    !Loop over all a_loc
!                    write(iout,*) "a_loc: ",a_loc

                    !Reject if a is occupied
                    if(IsOcc(iLutHF,a_loc)) cycle

                    Ka=G1(a_loc)%k
                    !=======================================
                    if (tUEG2) then
                        Ka=kvec(a_loc, 1:3)
                    end if
                    !=======================================

                    !Find k labels of b
                    kx=Ki(1)+Kj(1)-Ka(1)
                    if(abs(kx).gt.NMAXX) cycle
                    ky=Ki(2)+Kj(2)-Ka(2)
                    if(abs(ky).gt.NMAXY) cycle
                    kz=Ki(3)+Kj(3)-Ka(3)
                    if(abs(kz).gt.NMAXZ) cycle
                    !if(tGCutoff) then
                    !    length_g=real((kx-kj(1))**2+(ky-kj(2))**2+(kz-kj(3))**2)
                    !    length_g_2=real((kx-ki(1))**2+(ky-ki(2))**2+(kz-ki(3))**2)
                    !    if(length_g.gt.gCutoff.and.length_g_2.gt.gCutoff) cycle
                    !endif
                    length=real((kx**2)+(ky**2)+(kz**2),dp)
                    if(length.gt.OrbECutoff) cycle

                    !Find the actual k orbital
                    if(is_beta(a_loc)) then
                        !want alpha b orbital
                        OrbB=kPointToBasisFn(kx,ky,kz,2)
                    else
                        !want beta
                        OrbB=kPointToBasisFn(kx,ky,kz,1)
                    endif

                    !Reject k orbital if it is occupied or gt a
                    if(IsOcc(iLutHF,OrbB)) cycle
                    if(OrbB.ge.a_loc) cycle

!                    write(iout,*) "OrbB: ",OrbB
                    !Find det
                    call make_double (HFDet, nJ, elec1ind, elec2ind, a_loc, &
                                      orbB, ex, tParity)
                    !Sum in mp2 contrib
                    hel=get_helement_excit(HFDet,nJ,2,Ex,tParity)
                    H0tmp=getH0Element4(nJ,HFDet)
                    H0tmp=Fii-H0tmp
                    if(tMadelung) then
                        H0tmp=H0tmp+2.0_dp*Madelung
                    endif
                    mp2=mp2+(hel**2)/H0tmp
!                    write(iout,*) (hel**2),H0tmp
                enddo
            endif

        enddo

!        write(iout,*) "mp2: ",mp2
        mp2all=0.0_dp
        
        !Sum contributions across nodes.
        call MPISumAll(mp2,mp2all)
        write(iout,"(A,2G25.15)") "MP2 energy calculated: ",MP2All,MP2All+Hii
        call neci_flush(iout)

    end subroutine CalcUEGMP2
            

   subroutine SetupValidSpawned(WalkerListSize)
      use CalcData, only: MemoryFacSpawn
      implicit none
      integer(int64), intent(in) :: WalkerListSize
      integer ierr,i,j
      real(dp) Gap

      !When running normally, WalkerListSize will be equal to initwalkers
      !However, when reading in (and not continuing to grow) it should be equal to the number of dets in the popsfile
      MaxSpawned=NINT(MemoryFacSpawn*WalkerListSize*inum_runs)
!            WRITE(iout,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned
            
!      WRITE(iout,"(A)") "*Direct Annihilation* in use...Explicit load-balancing disabled."
      ALLOCATE(ValidSpawnedList(0:nNodes-1),stat=ierr)
      ! InitialSpawnedSlots is now filled later, once the number of particles
      ! wanted is known
      !(it can change according to the POPSFILE).
      ALLOCATE(InitialSpawnedSlots(0:nNodes-1),stat=ierr)
      ! InitialSpawnedSlots now holds the first free position in the 
      ! newly-spawned list for each processor, so it does not need to be 
      ! reevaluated each iteration.
!      MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
      Gap=REAL(MaxSpawned,dp)/REAL(nNodes,dp)
      do j=0,nNodes-1
          InitialSpawnedSlots(j)=NINT(Gap*j)+1
      enddo
      ! ValidSpawndList now holds the next free position in the newly-spawned
      ! list, but for each processor.
      ValidSpawnedList(:)=InitialSpawnedSlots(:)

   end subroutine

!This routine takes the walkers from all subspaces, constructs the hamiltonian, and diagonalises it.
!Currently, this only works in serial.
    subroutine DiagWalkerSubspace()
        implicit none
        integer :: i,iSubspaceSize,ierr,iSubspaceSizeFull
        real(dp) :: CurrentSign(lenof_sign)
        integer, allocatable :: ExpandedWalkerDets(:,:)
        integer(TagIntType) :: ExpandedWalkTag=0
        real(dp) :: GroundEFull,GroundEInit,CreateNan,ProjGroundEFull,ProjGroundEInit
        character(*), parameter :: t_r='DiagWalkerSubspace'

        if(nProcessors.gt.1) call stop_all(t_r,"Walker subspace diagonalisation only works in serial")
        if(lenof_sign.ne.1) call stop_all(t_r,'Cannot do Lanczos on complex orbitals.')

        if(tTruncInitiator) then
            !First, diagonalise initiator subspace
            write(iout,'(A)') 'Diagonalising initator subspace...'

            iSubspaceSize = 0
            do i=1,int(TotWalkers,sizeof_int)
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if((abs(CurrentSign(1)) > InitiatorWalkNo) .or. &
                        (DetBitEQ(CurrentDets(:,i),iLutHF,NIfDBO))) then
                    !Is allowed initiator. Add to subspace.
                    iSubspaceSize = iSubspaceSize + 1
                endif
            enddo

            write(iout,'(A,I12)') "Number of initiators found to diagonalise: ",iSubspaceSize
            allocate(ExpandedWalkerDets(NEl,iSubspaceSize),stat=ierr)
            call LogMemAlloc('ExpandedWalkerDets',NEl*iSubspaceSize,4,t_r,ExpandedWalkTag,ierr)

            iSubspaceSize = 0
            do i=1,int(TotWalkers,sizeof_int)
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if((abs(CurrentSign(1)) > InitiatorWalkNo) .or. &
                        (DetBitEQ(CurrentDets(:,i),iLutHF,NIfDBO))) then
                    !Is allowed initiator. Add to subspace.
                    iSubspaceSize = iSubspaceSize + 1
                    call decode_bit_det(ExpandedWalkerDets(:,iSubspaceSize),CurrentDets(:,i))
                endif
            enddo

            if(iSubspaceSize.gt.0) then
!Routine to diagonalise a set of determinants, and return the ground state energy
                call LanczosFindGroundE(ExpandedWalkerDets,iSubspaceSize,GroundEInit,ProjGroundEInit,.true.)
                write(iout,'(A,G25.10)') 'Ground state energy of initiator walker subspace = ',GroundEInit
            else
                CreateNan=-1.0_dp
                write(iout,'(A,G25.10)') 'Ground state energy of initiator walker subspace = ',sqrt(CreateNan)
            endif


            deallocate(ExpandedWalkerDets)
            call LogMemDealloc(t_r,ExpandedWalkTag)

        endif

        iSubspaceSizeFull = int(TotWalkers,sizeof_int)

        !Allocate memory for walker list.
        write(iout,'(A)') "Allocating memory for diagonalisation of full walker subspace"
        write(iout,'(A,I12,A)') "Size = ",iSubspaceSizeFull," walkers."

        allocate(ExpandedWalkerDets(NEl,iSubspaceSizeFull),stat=ierr)
        call LogMemAlloc('ExpandedWalkerDets',NEl*iSubspaceSizeFull,4,t_r,ExpandedWalkTag,ierr)
        do i=1,iSubspaceSizeFull
            call decode_bit_det(ExpandedWalkerDets(:,i),CurrentDets(:,i))
        enddo

!Routine to diagonalise a set of determinants, and return the ground state energy
        call LanczosFindGroundE(ExpandedWalkerDets,iSubspaceSizeFull,GroundEFull,ProjGroundEFull,.false.)

        write(iout,'(A,G25.10)') 'Ground state energy of full walker subspace = ',GroundEFull

        if(tTruncInitiator) then
            write(unitWalkerDiag,'(3I14,4G25.15)') Iter,iSubspaceSize,iSubspaceSizeFull,GroundEInit-Hii,    &
                    GroundEFull-Hii,ProjGroundEInit-Hii,ProjGroundEFull-Hii
        else
            write(unitWalkerDiag,'(2I14,2G25.15)') Iter,iSubspaceSizeFull,GroundEFull-Hii,ProjGroundEFull-Hii
        endif
        
        deallocate(ExpandedWalkerDets)
        call LogMemDealloc(t_r,ExpandedWalkTag)

    end subroutine DiagWalkerSubspace

!Routine which takes a set of determinants and returns the ground state energy
    subroutine LanczosFindGroundE(Dets,DetLen,GroundE,ProjGroundE,tInit)
        use DetCalcData, only : NKRY,NBLK,B2L,nCycle
        implicit none
        real(dp) , intent(out) :: GroundE,ProjGroundE
        logical , intent(in) :: tInit
        integer, intent(in) :: DetLen
        integer, intent(in) :: Dets(NEl,DetLen)
        integer :: gc,ierr,icmax,lenhamil,nKry1,nBlock,LScr,LIScr
        integer(n_int) :: ilut(0:NIfTot)
        logical :: tmc,tSuccess
        real(dp) , allocatable :: A_Arr(:,:),V(:),AM(:),BM(:),T(:),WT(:),SCR(:),WH(:),WORK2(:),V2(:,:),W(:)
        HElement_t, allocatable :: Work(:)
        HElement_t, allocatable :: CkN(:,:),Hamil(:),TruncWavefunc(:)
        HElement_t, allocatable :: Ck(:,:)  !This holds eigenvectors in the end
        HElement_t :: Num,Denom,HDiagTemp
        integer , allocatable :: LAB(:),NROW(:),INDEX(:),ISCR(:)
        integer(TagIntType) :: LabTag=0,NRowTag=0,ATag=0,VTag=0,AMTag=0,BMTag=0,TTag=0,tagCKN=0,tagCK=0
        integer(TagIntType) :: WTTag=0,SCRTag=0,ISCRTag=0,INDEXTag=0,WHTag=0,Work2Tag=0,V2Tag=0,tagW=0
        integer(TagIntType) :: tagHamil=0,WorkTag=0
        integer :: nEval,nBlocks,nBlockStarts(2),ExcitLev,iGetExcitLevel,i,nDoubles
        character(*), parameter :: t_r='LanczosFindGroundE'
        real(dp) :: norm
        integer :: PartInd,ioTrunc
        character(24) :: abstr
        
        if(DetLen.gt.1300) then
            nEval = 3
        else
            nEval = DetLen
        endif

        write(iout,'(A,I10)') 'DetLen = ',DetLen
!C..
        Allocate(Ck(DetLen,nEval), stat=ierr)
        call LogMemAlloc('CK',DetLen*nEval, 8,t_r, tagCK,ierr)
        CK=(0.0_dp)
!C..
        allocate(W(nEval), stat=ierr)
        call LogMemAlloc('W', nEval,8,t_r,tagW,ierr)
        W=0.0_dp
!        write(iout,*) "Calculating H matrix"
!C..We need to measure HAMIL and LAB first 
        allocate(NROW(DetLen),stat=ierr)
        call LogMemAlloc('NROW',DetLen,4,t_r,NROWTag,ierr)
        NROW(:)=0
        ICMAX=1
        TMC=.FALSE.
        call DETHAM(DetLen,NEL,Dets,HAMIL,LAB,NROW,.TRUE.,ICMAX,GC,TMC)
!        WRITE(iout,*) ' FINISHED COUNTING '
        WRITE(iout,*) "Allocating memory for hamiltonian: ",GC*2
        CALL neci_flush(iout)
!C..Now we know size, allocate memory to HAMIL and LAB
        LENHAMIL=GC
        Allocate(Hamil(LenHamil), stat=ierr)
        call LogMemAlloc('HAMIL', LenHamil, 8, t_r,tagHamil,ierr)
        HAMIL=(0.0_dp)
!C..
        ALLOCATE(LAB(LENHAMIL),stat=ierr)
        CALL LogMemAlloc('LAB',LenHamil,4,t_r,LabTag,ierr)

        LAB(1:LENHAMIL)=0
!C..Now we store HAMIL and LAB 
        CALL DETHAM(DetLen,NEL,Dets,HAMIL,LAB,NROW,.FALSE.,ICMAX,GC,TMC)

        if(DetLen.gt.1300) then
            !Lanczos diag

            Allocate(CkN(DetLen,nEval), stat=ierr)
            call LogMemAlloc('CKN',DetLen*nEval, 8,t_r, tagCKN,ierr)
            CKN=(0.0_dp)

            NKRY1=NKRY+1
            NBLOCK=MIN(NEVAL,NBLK)
            LSCR=MAX(DetLen*NEVAL,8*NBLOCK*NKRY)
            LISCR=6*NBLOCK*NKRY
    !C..
    !       write (iout,'(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
    !C..Set up memory for FRSBLKH

            ALLOCATE(A_Arr(NEVAL,NEVAL),stat=ierr)
            CALL LogMemAlloc('A_Arr',NEVAL**2,8,t_r,ATag,ierr)
            A_Arr=0.0_dp
    !C..
    !C,, W is now allocated with CK
    !C..
            ALLOCATE(V(DetLen*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('V',DetLen*NBLOCK*NKRY1,8,t_r,VTag,ierr)
            V=0.0_dp
    !C..   
            ALLOCATE(AM(NBLOCK*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('AM',NBLOCK*NBLOCK*NKRY1,8,t_r,AMTag,ierr)
            AM=0.0_dp
    !C..
            ALLOCATE(BM(NBLOCK*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('BM',NBLOCK*NBLOCK*NKRY,8,t_r,BMTag,ierr)
            BM=0.0_dp
    !C..
            ALLOCATE(T(3*NBLOCK*NKRY*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('T',3*NBLOCK*NKRY*NBLOCK*NKRY,8,t_r,TTag,ierr)
            T=0.0_dp
    !C..
            ALLOCATE(WT(NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('WT',NBLOCK*NKRY,8,t_r,WTTag,ierr)
            WT=0.0_dp
    !C..
            ALLOCATE(SCR(LScr),stat=ierr)
            CALL LogMemAlloc('SCR',LScr,8,t_r,SCRTag,ierr)
            SCR=0.0_dp
            ALLOCATE(ISCR(LIScr),stat=ierr)
            CALL LogMemAlloc('IScr',LIScr,4,t_r,IScrTag,ierr)
            ISCR(1:LISCR)=0
            ALLOCATE(INDEX(NEVAL),stat=ierr)
            CALL LogMemAlloc('INDEX',NEVAL,4,t_r,INDEXTag,ierr)
            INDEX(1:NEVAL)=0
    !C..
            ALLOCATE(WH(DetLen),stat=ierr)
            CALL LogMemAlloc('WH',DetLen,8,t_r,WHTag,ierr)
            WH=0.0_dp
            ALLOCATE(WORK2(3*DetLen),stat=ierr)
            CALL LogMemAlloc('WORK2',3*DetLen,8,t_r,WORK2Tag,ierr)
            WORK2=0.0_dp
            ALLOCATE(V2(DetLen,NEVAL),stat=ierr)
            CALL LogMemAlloc('V2',DetLen*NEVAL,8,t_r,V2Tag,ierr)
            V2=0.0_dp
    !C..Lanczos iterative diagonalising routine
            CALL NECI_FRSBLKH(DetLen,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A_Arr,W,V,AM,BM,T,WT, &
             &  SCR,ISCR,INDEX,NCYCLE,B2L,.true.,.false.,.false.,.false.)

            !Eigenvalues may come out wrong sign - multiply by -1
            if(W(1).gt.0.0_dp) then
                GroundE = -W(1)
            else
                GroundE = W(1)
            endif

            deallocate(A_Arr,V,AM,BM,T,WT,SCR,WH,V2,CkN,Index,IScr)
            call LogMemDealloc(t_r,ATag)
            call LogMemDealloc(t_r,VTag)
            call LogMemDealloc(t_r,AMTag)
            call LogMemDealloc(t_r,BMTag)
            call LogMemDealloc(t_r,TTag)
            call LogMemDealloc(t_r,tagCKN)
            call LogMemDealloc(t_r,WTTag)
            call LogMemDealloc(t_r,SCRTag)
            call LogMemDealloc(t_r,ISCRTag)
            call LogMemDealloc(t_r,IndexTag)
            call LogMemDealloc(t_r,WHTag)
            call LogMemDealloc(t_r,V2Tag)
        else
            !Full diag
            allocate(Work(4*DetLen),stat=ierr)
            call LogMemAlloc('Work',4*DetLen,8,t_r,WorkTag,ierr)
            allocate(Work2(3*DetLen),stat=ierr)
            call logMemAlloc('Work2',3*DetLen,8,t_r,Work2Tag,ierr)
            nBlockStarts(1) = 1
            nBlockStarts(2) = DetLen+1
            nBlocks = 1
            call HDIAG_neci(DetLen,Hamil,Lab,nRow,CK,W,Work2,Work,nBlockStarts,nBlocks)
            GroundE = W(1)
            deallocate(Work)
            call LogMemDealloc(t_r,WorkTag)
        endif

        !Calculate proje.
        nDoubles = 0
        Num = 0.0_dp
        do i=1,DetLen
            ExcitLev = iGetExcitLevel(HFDet,Dets(:,i),NEl)
            if((ExcitLev.eq.1).or.(ExcitLev.eq.2)) then
                HDiagTemp = get_helement(HFDet,Dets(:,i),ExcitLev)
                Num = Num + (HDiagTemp * CK(i,1))
                nDoubles = nDoubles + 1
            elseif(ExcitLev.eq.0) then
                Denom = CK(i,1)
            endif
        enddo
        ProjGroundE = (Num / Denom) + Hii
!        write(iout,*) "***", nDoubles

        if(tHistSpawn.and..not.tInit) then
            !Write out the ground state wavefunction for this truncated calculation.
            !This needs to be sorted, into the order given in the original FCI in DetCalc.
            allocate(TruncWavefunc(Det),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc error")
            TruncWavefunc(:) = 0.0_dp

            Norm = 0.0_dp
            do i=1,DetLen
            !Loop over iluts in instantaneous list.

                call EncodeBitDet(Dets(:,i),iLut)
                !Calculate excitation level (even if hf isn't in list)
                ExcitLev = FindBitExcitLevel(iLutHF,iLut,NEl)

                if(ExcitLev.eq.nel) then
                    call BinSearchParts2(iLut,FCIDetIndex(ExcitLev),Det,PartInd,tSuccess)
                elseif(ExcitLev.eq.0) then
                    PartInd = 1
                    tSuccess = .true.
                else
                    call BinSearchParts2(iLut,FCIDetIndex(ExcitLev),FCIDetIndex(ExcitLev+1)-1,PartInd,tSuccess)
                endif
                if(.not.tSuccess) then
                    call stop_all(t_r,"Error here as cannot find corresponding determinant in FCI expansion")
                endif
                TruncWavefunc(PartInd) = CK(i,1)

                !Check normalisation
                Norm = Norm + CK(i,1)**2.0_dp
            enddo
            Norm = sqrt(Norm)
            if(abs(Norm-1.0_dp).gt.1.0e-7_dp) then
                write(iout,*) "***",norm
                call warning_neci(t_r,"Normalisation not correct for diagonalised wavefunction!")
            endif

            !Now write out...
            abstr=''
            write(abstr,'(I12)') Iter
            abstr='TruncWavefunc-'//adjustl(abstr)
            ioTrunc = get_free_unit()
            open(ioTrunc,file=abstr,status='unknown')
            do i=1,Det
                write(ioTrunc,"(I13,F25.16)") i,TruncWavefunc(i)/norm
            enddo
            close(ioTrunc)

            deallocate(TruncWavefunc)

        endif


        deallocate(Work2,W,Hamil,Lab,nRow,CK)
        call LogMemDealloc(t_r,Work2Tag)
        call LogMemDealloc(t_r,tagW)
        call LogMemDealloc(t_r,tagHamil)
        call LogMemDealloc(t_r,LabTag)
        call LogMemDealloc(t_r,NRowTag)
        call LogMemDealloc(t_r,tagCK)
    
    end subroutine LanczosFindGroundE
    
    !Routine to print the highest populated determinants at the end of a run
    SUBROUTINE PrintHighPops()
        use DetBitOps, only : sign_lt,sign_gt
        use LoggingData, only: iHighPopWrite
        real(dp), dimension(lenof_sign) :: SignCurr, LowSign
        integer :: ierr,i,j,counter,ExcitLev,SmallestPos,HighPos,nopen
        real(dp) :: HighSign,reduce_in(1:2),reduce_out(1:2),Norm,AllNorm
        integer(n_int) , allocatable :: LargestWalkers(:,:)
        integer(n_int) , allocatable :: GlobalLargestWalkers(:,:)
        integer(n_int) :: HighestDet(0:NIfTot)
        integer, allocatable :: GlobalProc(:)
        character(len=*), parameter :: t_r='PrintHighPops'

        !Allocate memory to hold highest iHighPopWrite determinants
        allocate(LargestWalkers(0:NIfTot,iHighPopWrite),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"error allocating here")

        ! Return the most populated states in CurrentDets on *this* processor only.
        call return_most_populated_states(iHighPopWrite, LargestWalkers, norm)

        call MpiSum(norm,allnorm)
        norm=sqrt(allnorm)

!        write(iout,*) "Highest weighted dets on this process:"
!        do i=1,iHighPopWrite
!            write(iout,*) LargestWalkers(:,i)
!        enddo

        !Now have sorted list of the iHighPopWrite largest weighted determinants on the process
        if(iProcIndex.eq.Root) then
            allocate(GlobalLargestWalkers(0:NIfTot,iHighPopWrite),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"error allocating here")
            GlobalLargestWalkers(:,:)=0
            allocate(GlobalProc(iHighPopWrite),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"error allocating here")
            GlobalProc(:)=0
        endif

        do i=1,iHighPopWrite

            ! Find highest sign on each processor. Since all lists are
            ! sorted, this is just teh first nonzero value.
            do j=iHighPopWrite,1,-1
                call extract_sign (LargestWalkers(:,j), SignCurr)
                if (any(LargestWalkers(:,j) /= 0)) then
                    
#ifdef __CMPLX
                    HighSign = sqrt(SignCurr(1)**2 + SignCurr(lenof_sign)**2)
#else
                    HighSign = real(abs(SignCurr(1)),dp)
#endif

                    ! We have the largest sign
                    HighPos = j
                    exit
                end if
            end do
            reduce_in=(/ HighSign,real(iProcIndex,dp) /)
            call MPIAllReduceDatatype(reduce_in,1,MPI_MAXLOC,MPI_2DOUBLE_PRECISION,reduce_out)
            !Now, reduce_out(2) has the process of the largest weighted determinant - broadcast it!
            if(iProcIndex.eq.nint(reduce_out(2))) then
                HighestDet(0:NIfTot) = LargestWalkers(:,HighPos)
            else
                HighestDet(0:NIfTot) = 0
            endif
            call MPIBCast(HighestDet(:),NIfTot+1,nint(reduce_out(2)))
            if(iProcIndex.eq.Root) then
                GlobalLargestWalkers(0:NIfTot,i) = HighestDet(:)
                GlobalProc(i) = nint(reduce_out(2))
            endif

            !Now delete this highest determinant from the list on the corresponding processor
            if(iProcIndex.eq.nint(reduce_out(2))) then
                LargestWalkers(:,HighPos) = 0
                !No need to resort any more
!                call sort(LargestWalkers(:,1:iHighPopWrite), sign_lt, sign_gt)
            endif
        enddo

        if(iProcIndex.eq.Root) then
            !Now print out the info contained in GlobalLargestWalkers and GlobalProc

            counter=0
            do i=1,iHighPopWrite
                !How many non-zero determinants do we actually have?
                call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
#ifdef __CMPLX
                HighSign=sqrt(real(SignCurr(1),dp)**2+real(SignCurr(lenof_sign),dp)**2)
#else
                HighSign=real(abs(SignCurr(1)),dp)
#endif
                if (HighSign > 1.0e-7_dp) counter = counter + 1
            enddo


            write(iout,*) ""
            write(iout,'(A)') "Current reference: "
            call write_det (iout, ProjEDet, .true.)
            call writeDetBit(iout,iLutRef,.true.)
            write(iout,*) ""
            write(iout,"(A,I10,A)") "Most occupied ",counter," determinants as excitations from reference: "
            write(iout,*) 
            if(lenof_sign.eq.1) then
                if(tHPHF) then
                    write(iout,"(A)") " Excitation   ExcitLevel   Seniority    Walkers    Weight    Init?   Proc  Spin-Coup?"    
                else
                    write(iout,"(A)") " Excitation   ExcitLevel   Seniority    Walkers    Weight    Init?   Proc"    
                endif
            else
                if(tHPHF) then
                    write(iout,"(A)") " Excitation   ExcitLevel Seniority  Walkers(Re)   Walkers(Im)  Weight   &
                                        &Init?(Re)   Init?(Im)   Proc  Spin-Coup?"
                else
                    write(iout,"(A)") " Excitation   ExcitLevel Seniority   Walkers(Re)   Walkers(Im)  Weight   &
                                        &Init?(Re)   Init?(Im)   Proc"
                endif
            endif
            do i=1,counter
!                call WriteBitEx(iout,iLutRef,GlobalLargestWalkers(:,i),.false.)
                call WriteDetBit(iout,GlobalLargestWalkers(:,i),.false.)
                Excitlev=FindBitExcitLevel(iLutRef,GlobalLargestWalkers(:,i),nEl)
                write(iout,"(I5)",advance='no') Excitlev
                nopen=count_open_orbs(GlobalLargestWalkers(:,i))
                write(iout,"(I5)",advance='no') nopen
                call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
                do j=1,lenof_sign
                    write(iout,"(G16.7)",advance='no') SignCurr(j)
                enddo
#ifdef __CMPLX
                HighSign=sqrt(real(SignCurr(1),dp)**2+real(SignCurr(lenof_sign),dp)**2)
#else
                HighSign=real(abs(SignCurr(1)),dp)
#endif
                if(tHPHF.and.(.not.TestClosedShellDet(GlobalLargestWalkers(:,i)))) then 
                    !Weight is proportional to nw/sqrt(2)
                    write(iout,"(F9.5)",advance='no') (HighSign/sqrt(2.0_dp))/norm 
                else
                    write(iout,"(F9.5)",advance='no') HighSign/norm 
                endif
                do j=1,lenof_sign
                    if(.not.tTruncInitiator) then
                        write(iout,"(A3)",advance='no') 'Y'
                    else
                        if(test_flag(GlobalLargestWalkers(:,i),flag_is_initiator(j))) then
                            write(iout,"(A3)",advance='no') 'Y'
                        else
                            write(iout,"(A3)",advance='no') 'N'
                        endif
                    endif
                enddo
                if(tHPHF.and.(.not.TestClosedShellDet(GlobalLargestWalkers(:,i)))) then 
                    write(iout,"(I7)",advance='no') GlobalProc(i)
                    write(iout,"(A3)") "*"
                else
                    write(iout,"(I7)") GlobalProc(i)
                endif
            enddo

            if(tHPHF) then
                write(iout,"(A)") " * = Spin-coupled function implicitly has time-reversed determinant with same weight."
            endif

            deallocate(GlobalLargestWalkers,GlobalProc)
            write(iout,*) ""
        endif

        deallocate(LargestWalkers)

    END SUBROUTINE PrintHighPops
            
    !Routine to just calculate errors from FCIMCStats file
    subroutine Standalone_Errors()
        use sym_mod, only: getsym
#ifdef MOLPRO
        use outputResult
        integer :: nv,ityp(1)
#endif
        real(dp) :: mean_ProjE_re,mean_ProjE_im,mean_Shift
        real(dp) :: ProjE_Err_re,ProjE_Err_im,Shift_Err
        logical :: tNoProjEValue,tNoShiftValue
        integer :: iroot,isymh,i
        TYPE(BasisFn) RefSym
        HElement_t :: h_tmp
        real(dp) :: Hii,BestEnergy,EnergyDiff
#ifdef MOLPRO
        real(dp) :: get_scalar
        include "common/molen"
#endif

        !Automatic error analysis
        call error_analysis(tSinglePartPhase(1),iBlockingIter(1),mean_ProjE_re,ProjE_Err_re,  &
            mean_ProjE_im,ProjE_Err_im,mean_Shift,Shift_Err,tNoProjEValue,tNoShiftValue, &
            equilshift=iBlockEquilShift,equilproje=iBlockEquilProjE)
        call MPIBCast(ProjectionE)
        call MPIBCast(mean_ProjE_re)
        call MPIBCast(ProjE_Err_re)
        call MPIBCast(mean_ProjE_im)
        call MPIBCast(ProjE_Err_im)
        call MPIBCast(mean_Shift)
        call MPIBCast(Shift_Err)
        call MPIBCast(tNoProjEValue)
        call MPIBCast(tNoShiftValue)

        h_tmp = get_helement (FDet, FDet, 0)
        Hii = real(h_tmp, dp)

        iroot=1
        CALL GetSym(FDet,NEl,G1,NBasisMax,RefSym)
        isymh=int(RefSym%Sym%S,sizeof_int)+1
        write (iout,10101) iroot,isymh
10101   format(//'RESULTS FOR STATE',i2,'.',i1/'====================='/)
        write (iout,'('' Current reference energy'',T52,F19.12)') Hii 
        if(tNoProjEValue) then
            write(iout,'('' No projected energy value could be obtained'',T52)')
        else
            write (iout,'('' Projected correlation energy'',T52,F19.12)') mean_ProjE_re
            write (iout,'('' Estimated error in Projected correlation energy'',T52,F19.12)') ProjE_Err_re
        endif
        if(lenof_sign.eq.2) then
            write (iout,'('' Projected imaginary energy'',T52,F19.12)') mean_ProjE_im
            write (iout,'('' Estimated error in Projected imaginary energy'',T52,F19.12)') ProjE_Err_im
        endif
        if(tNoShiftValue) then
            write(iout,'('' No shift energy value could be obtained'',T52)')
        else
            write (iout,'('' Shift correlation energy'',T52,F19.12)') mean_Shift
            write (iout,'('' Estimated error in shift correlation energy'',T52,F19.12)') shift_err
        endif

        !Do shift and projected energy agree?
        write(iout,"(A)")
        if(tNoProjEValue.and.tNoShiftValue) return 
        EnergyDiff = abs(mean_Shift-mean_ProjE_re)
        if(EnergyDiff.le.sqrt(shift_err**2+ProjE_Err_re**2)) then
            write(iout,"(A,F15.8)") " Projected and shift energy estimates agree " &
                & //"within errorbars: EDiff = ",EnergyDiff
        elseif(EnergyDiff.le.sqrt((max(shift_err,ProjE_Err_re)*2)**2+min(shift_err,ProjE_Err_re)**2)) then
            write(iout,"(A,F15.8)") " Projected and shift energy estimates agree to within two sigma " &
                & //"of largest error: EDiff = ",EnergyDiff
        else
            write(iout,"(A,F15.8)") " Projected and shift energy estimates do not agree to " &
                & //"within approximate errorbars: EDiff = ",EnergyDiff
        endif
        if(ProjE_Err_re.lt.shift_err) then
            BestEnergy = mean_ProjE_re + Hii
        else
            BestEnergy = mean_shift + Hii
        endif
        write(iout,"(A)")
        write(iout,"(A,F20.8,A,G15.6)") " Total projected energy ", &
            mean_ProjE_re+Hii," +/- ",ProjE_Err_re
        write(iout,"(A,F20.8,A,G15.6)") " Total shift energy     ", &
            mean_shift+Hii," +/- ",shift_err

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
        call output_result('FCIQMC','FCIQMC_ERR',min(ProjE_Err_re,shift_err),iroot,isymh)
        if (iroot.eq.1) call clearvar('FCIQMC_ERR')
        call setvar('FCIQMC_ERR',min(ProjE_Err_re,shift_err),'AU',ityp,1,nv,iroot)
#endif
        write(iout,"(/)")

    end subroutine Standalone_Errors

END MODULE FciMCParMod

        

! This is the same as BinSearchParts1, but this time, it searches though the 
! full list of determinants created by the full diagonalizer when the 
! histogramming option is on. 
!
! This is outside the module so it is accessible to AnnihilateMod
SUBROUTINE BinSearchParts2(iLut,MinInd,MaxInd,PartInd,tSuccess)
    use DetCalcData , only : FCIDets
    use DetBitOps, only: DetBitLT
    use constants, only: n_int
    use bit_reps, only: NIfTot,NIfDBO
    IMPLICIT NONE
    INTEGER :: MinInd,MaxInd,PartInd
    INTEGER(KIND=n_int) :: iLut(0:NIfTot)
    INTEGER :: i,j,N,Comp
    LOGICAL :: tSuccess

!    WRITE(iout,*) "Binary searching between ",MinInd, " and ",MaxInd
!    CALL neci_flush(iout)
    i=MinInd
    j=MaxInd
    IF(i-j.eq.0) THEN
        Comp=DetBitLT(FCIDets(:,MaxInd),iLut(:),NIfDBO)
        IF(Comp.eq.0) THEN
            tSuccess=.true.
            PartInd=MaxInd
            RETURN
        ELSE
            tSuccess=.false.
            PartInd=MinInd
        ENDIF
    ENDIF
    do while(j-i.gt.0)  !End when the upper and lower bound are the same.
        N=(i+j)/2       !Find the midpoint of the two indices
!        WRITE(iout,*) i,j,n

        ! Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is 
        ! more or 0 if they are the same
        Comp=DetBitLT(FCIDets(:,N),iLut(:),NIfDBO)

        IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
            tSuccess=.true.
            PartInd=N
            RETURN
        ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
            ! The value of the determinant at N is LESS than the determinant 
            ! we're looking for. Therefore, move the lower bound of the 
            ! search up to N. However, if the lower bound is already equal to
            ! N then the two bounds are consecutive and we have failed...
            i=N
        ELSEIF(i.eq.N) THEN


            IF(i.eq.MaxInd-1) THEN
                ! This deals with the case where we are interested in the 
                ! final/first entry in the list. Check the final entry of the
                ! list and leave. We need to check the last index.
                Comp=DetBitLT(FCIDets(:,i+1),iLut(:),NIfDBO)
                IF(Comp.eq.0) THEN
                    tSuccess=.true.
                    PartInd=i+1
                    RETURN
                ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                    tSuccess=.false.
                    PartInd=i+1
                    RETURN
                ELSE
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ENDIF

            ELSEIF(i.eq.MinInd) THEN
                tSuccess=.false.
                PartInd=i
                RETURN
            ELSE
                i=j
            ENDIF


        ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
            j=N
        ELSE
!We have failed - exit loop
            i=j
        ENDIF

    enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
    tSuccess=.false.
    PartInd=MAX(MinInd,i-1)

END SUBROUTINE BinSearchParts2
   
