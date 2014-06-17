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
                          tAntiSym_MI, MolproID, tGenHelWeighted, &
                          tGen_4ind_weighted, tMomInv, tGen_4ind_reverse
    use bit_rep_data, only: extract_sign, flag_trial, flag_connected
    use bit_reps, only: NIfD, NIfTot, NIfDBO, NOffY, decode_bit_det, &
                        encode_bit_rep, encode_det, extract_bit_rep, &
                        test_flag, set_flag, extract_flags, &
                        flag_is_initiator, clear_all_flags, &
                        extract_sign, nOffSgn, flag_make_initiator, &
                        flag_parent_initiator, encode_sign, &
                        decode_bit_det_chunks, &
                        clr_flag, flag_trial, flag_connected, nOffFlag, &
                        flag_deterministic, flag_determ_parent, clr_flag, &
                        extract_part_sign, encode_part_sign
    use CalcData, only: InitWalkers, NMCyc, DiagSft, Tau, SftDamp, StepsSft, &
                        OccCASorbs, VirtCASorbs, tFindGroundDet, NEquilSteps,&
                        tReadPops, tRegenDiagHEls, iFullSpaceIter, MaxNoAtHF,&
                        GrowMaxFactor, CullFactor, tStartSinglePart, tCCMC, &
                        ScaleWalkers, HFPopThresh, tTruncCAS, AvMCExcits, &
                        tTruncInitiator, tDelayTruncInit, IterTruncInit, &
                        NShiftEquilSteps, tWalkContGrow, &
                        tAddToInitiator, InitiatorWalkNo, tInitIncDoubs, &
                        tRetestAddtoInit, tReadPopsChangeRef, &
                        tReadPopsRestart, tCheckHighestPopOnce, &
                        iRestartWalkNum, tRestartHighPop, FracLargerDet, &
                        tChangeProjEDet, tCheckHighestPop, tSpawnSpatialInit,&
                        MemoryFacInit, tMaxBloom, tTruncNOpen, tFCIMC, &
                        trunc_nopen_max, tSpawn_Only_Init, RealSpawnCutoff, &
                        TargetGrowRate, TargetGrowRateWalk, tShiftonHFPop, &
                        tContinueAfterMP2,iExitWalkers,MemoryFacPart, &
                        tAllRealCoeff, tRealCoeffByExcitLevel, tPopsMapping, &
                        tSpawn_Only_Init_Grow, RealCoeffExcitThresh, &
                        tRealSpawnCutoff, RealSpawnCutoff, tDetermProj, &
                        tJumpShift, tUseRealCoeffs, tSpatialOnlyHash, &
                        tSemiStochastic, tTrialWavefunction
    use spatial_initiator, only: add_initiator_list, rm_initiator_list
    use HPHFRandExcitMod, only: FindExcitBitDetSym, gen_hphf_excit
    use MomInvRandExcit, only: gen_MI_excit
    use MomInv, only: IsMomSelfInv, CalcMomAllowedBitDet
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
    use IntegralsData, only: fck, NMax, UMat, tPartFreezeCore, NPartFrozen, &
                             NHolesFrozen, tPartFreezeVirt, NVirtPartFrozen, &
                             NElVirtFrozen
    use LoggingData, only: iWritePopsEvery, TPopsFile, iPopsPartEvery, tBinPops, &
                           iWriteHistEvery, tHistEnergies, FCIMCDebug, &
                           IterShiftBlock, AllHistInitPops, &
                           OffDiagBinRange, OffDiagMax, AllHistInitPopsTag, &
                           tLogComplexPops, tPrintFCIMCPsi, tCalcFCIMCPsi, &
                           NHistEquilSteps, tPrintOrbOcc, StartPrintOrbOcc, &
                           tPrintOrbOccInit, tHFPopStartBlock, tIterStartBlock, &
                           IterStartBlocking, HFPopStartBlocking, &
                           tInitShiftBlocking, tHistHamil, iWriteHamilEvery, &
                           HistInitPopsTag, OrbOccs, OrbOccsTag, &
                           tPrintPopsDefault, iWriteBlockingEvery, &
                           tBlockEveryIteration, tHistInitPops, HistInitPopsIter,&
                           HistInitPops, DoubsUEG, DoubsUEGLookup, DoubsUEGStore,&
                           tPrintDoubsUEG, StartPrintDoubsUEG, tCalcInstantS2, &
                           instant_s2_multiplier, tMCOutput, tSplitProjEHist, &
                           tSplitProjEHistG, tSplitProjEHistK3, iProjEBins, &
                           tDiagWalkerSubspace,iDiagSubspaceIter, &
                           tRDMonFly, IterRDMonFly,RDMExcitLevel, RDMEnergyIter, &
                           tChangeVarsRDM, tExplicitAllRDM, tHF_Ref_Explicit, &
                           tHF_S_D_Ref, tHF_S_D, &
                           tDiagWalkerSubspace, iDiagSubspaceIter, &
                           tCalcInstantS2Init, instant_s2_multiplier_init, &
                           tJustBlocking, iBlockEquilShift, iBlockEquilProjE, &
                           tDiagAllSpaceEver, tCalcVariationalEnergy, tCompareTrialAmps, &
                           compare_amps_period, tNoNewRDMContrib, &
                           log_cont_time_survivals, tNoWarnIC0Bloom, &
                           tLogPopsMaxTau, tFCIMCStats2, tHistExcitToFrom
    use hist, only: init_hist_spin_dist, clean_hist_spin_dist, &
                    hist_spin_dist, ilut_spindist, tHistSpinDist, &
                    write_clear_hist_spin_dist, hist_spin_dist_iter, &
                    test_add_hist_spin_dist_det, add_hist_energies, &
                    add_hist_spawn, tHistSpawn, AllHistogramEnergy, &
                    AllHistogram, HistogramEnergy, Histogram, AllInstHist, &
                    InstHist, HistMinInd, project_spins, calc_s_squared, &
                    project_spin_csfs, calc_s_squared_multi, &
                    calc_s_squared_star, init_hist_excit_tofrom, &
                    add_hist_excit_tofrom, write_zero_hist_excit_tofrom, &
                    clean_hist_excit_tofrom
    use hist_data, only: beforenormhist, HistMinInd2, HistMinInd2, BinRange, &
                         iNoBins
    USE SymData , only : nSymLabels, Sym_Psi
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE Parallel_neci
    USE FciMCData
    USE AnnihilationMod, only: DirectAnnihilation, RemoveDetHashIndex
    use PopsfileMod, only: ReadFromPopsfilePar, FindPopsfileVersion, &
                           WriteToPopsFileParOneArr, open_pops_head, &
                           readpopsheadv3, readpopsheadv4, CheckPopsParams, &
                           ReadFromPopsFile
    use sort_mod
    use DetBitops, only: EncodeBitDet, DetBitEQ, DetBitLT, FindExcitBitDet, &
                         FindBitExcitLevel, countbits, TestClosedShellDet, &
                         FindSpatialBitExcitLevel, IsAllowedHPHF, count_open_orbs, &
                         ilut_gt, get_bit_excitmat
    use hash , only : DetermineDetNode, FindWalkerHash
    use csf, only: get_csf_bit_yama, iscsf, csf_orbital_mask, get_csf_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_spawn_sign, hphf_off_diag_helement_spawn
    use MI_integrals, only: MI_diag_helement, MI_spawn_sign, &
                            MI_off_diag_helement_spawn, MI_off_diag_helement
    use util_mod
    use constants
    use soft_exit, only: ChangeVars 
    use FciMCLoggingMod, only: FinaliseBlocking, FinaliseShiftBlocking, &
                               PrintShiftBlocking, PrintBlocking, &
                               SumInErrorContrib, WriteInitPops, &
                               InitErrorBlocking, InitShiftErrorBlocking, &
                               SumInShiftErrorContrib
    use RotateOrbsMod, only: RotateOrbs
    use NatOrbsMod, only: PrintOrbOccs,PrintDoubUEGOccs
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
                         zero_rdms, fill_rdm_softexit, store_parent_with_spawned, &
                         fill_rdm_diag_currdet_norm,fill_rdm_diag_currdet_hfsd, calc_rdmbiasfac, tFinalRDMIter
    use determ_proj, only: perform_determ_proj
    use semi_stoch_gen, only: init_semi_stochastic, enumerate_sing_doub_kpnt
    use semi_stoch_procs, only: deterministic_projection, return_most_populated_states, &
                                end_semistoch, is_core_state, return_mp1_amp_and_mp2_energy, &
                                recalc_core_hamil_diag
    use trial_wf_gen, only: init_trial_wf, update_compare_trial_file, end_trial_wf
    use gndts_mod, only: gndts
    use sort_mod
    use get_excit, only: make_double
    use sltcnd_mod, only: sltcnd_excit
    use excit_gens_int_weighted, only: gen_excit_hel_weighted, &
                                       gen_excit_4ind_weighted, &
                                       test_excit_gen_4ind, &
                                       gen_excit_4ind_reverse
    use procedure_pointers
    use tau_search, only: init_tau_search, log_spawn_magnitude, update_tau, &
                          log_death_magnitude

    implicit none

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)
        use LoggingData, only: PopsfileTimer
        use nElRDMMod, only: InitRDM
        use sym_mod, only: getsym
        use SystemData, only: tUEG2
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
        call init_fcimc_fn_pointers () 

        if(n_int.eq.4) CALL Stop_All('Setup Parameters', &
                'Use of RealCoefficients does not work with 32 bit integers due to the use &
                & of the transfer operation from dp reals to 64 bit integers.')

        ! If performing a deterministic projection instead of an FCIQMC calc:
        if (tDetermProj) then
            call perform_determ_proj()
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
!            WRITE(iout,*) 'Iter',Iter

            if(iProcIndex.eq.root) s_start=neci_etime(tstart)
            
            if(tRDMonFly .and. (.not. tFillingExplicRDMonFly) &
                & .and. (.not.tFillingStochRDMonFly)) call check_start_rdm()

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

                IF((tTruncCAS.or.tTruncSpace.or.tTruncInitiator).and.(Iter.gt.iFullSpaceIter).and.(iFullSpaceIter.ne.0)) THEN
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
                IF (proje_update_comb) CALL update_linear_comb_coeffs()
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
            IF(tHistHamil.and.(mod(Iter,iWriteHamilEvery).eq.0)) THEN
                CALL WriteHamilHistogram()
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

        IF(tHistSpawn) CALL WriteHistogram()

        IF(tHistHamil) CALL WriteHamilHistogram()

        IF(tHistEnergies) CALL WriteHistogramEnergies()

        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
        ENDIF

        IF(tFillingStochRDMonFly.or.&
            tFillingExplicRDMonFly.or.tHF_Ref_Explicit) CALL FinaliseRDM()

        IF(tPrintDoubsUEG) THEN
            CALL PrintDoubUEGOccs(DoubsUEG)
        ENDIF

        call PrintHighPops()

        !Close open files.
        IF(iProcIndex.eq.Root) THEN
            CLOSE(fcimcstats_unit)
            if (inum_runs.eq.2) CLOSE(fcimcstats_unit2)
            IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(initiatorstats_unit)
            IF(tLogComplexPops) CLOSE(complexstats_unit)
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

    ! 
    ! This is a null routine for encoding spawned sites
    ! --> DOES NOTHING!!!
    subroutine null_encode_child (ilutI, ilutJ, ic, ex)
        use SystemData, only: nel
        use bit_reps, only: niftot
        use constants, only: n_int
        implicit none
        integer(kind=n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: ic, ex(2,2)
        integer(kind=n_int), intent(inout) :: ilutj(0:niftot)

        ! Avoid compiler warnings
        integer :: iUnused
        integer(n_int) :: iUnused2
        iLutJ(0) = iLutJ(0); iUnused = IC; iUnused = ex(2,2)
        iUnused2 = iLutI(0)
    end subroutine

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
        logical :: tParity, tSuccess
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

        IF(tDelayTruncInit.and.(Iter.ge.IterTruncInit)) THEN 
            IF(Iter.eq.IterTruncInit) THEN
                ! Why do this? Why not just get all procs to do division?
                IF(iProcIndex.eq.root) THEN
                    Tau=Tau/10.0_dp
                    WRITE(iout,'(A,F10.5)') 'Beginning truncated initiator calculation and reducing timestep by " &
                        &//"a factor of 10. New tau is : ',Tau
                ENDIF
                CALL MPIBCast(Tau)
            ENDIF
            tTruncInitiator=.true.
        ENDIF
        MaxInitPopPos=0.0
        MaxInitPopNeg=0.0
        HighPopNeg=1
        HighPopPos=1
        parent_flags=0
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

        ! Initialise histograms if necessary
        call InitHistMin()

        ! This is a bit of a hack based on the fact that we mean something 
        ! different by exFlag for CSFs than in normal determinential code.
        ! It would be nice to fix this properly
        if (tCSF) exFlag = 7

        ! Update here the sum for the projected-energy denominator if projecting
        ! onto a linear combination of determinants.
        if (proje_linear_comb .and. (.not. proje_spatial) .and. nproje_sum > 1) then
           sum_proje_denominator=0
           do i = 1, nproje_sum
              proc = DetermineDetNode (proje_ref_dets(:,i), 0)
              if (iProcIndex == proc) then
                    pos = binary_search (CurrentDets(:,1:TotWalkers), &
                         proje_ref_iluts(:,i), NIfD+1)
                    if (pos > 0) then
                        call extract_sign (CurrentDets(:,pos), sgn)
                        do k=1, inum_runs
                            delta(k) = ARR_RE_OR_CPLX(sgn,k) * proje_ref_coeffs(i)
                            cyc_proje_denominator(k) = cyc_proje_denominator(k) + delta(k)
                            sum_proje_denominator(k) = sum_proje_denominator(k) + delta(k)
                        enddo
                    endif
                endif
            enddo
        endif

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
            if(mod((Iter - IterRDMStart + 1),RDMEnergyIter).eq.0) then 
                ! RDM energy is being printed, calculate the diagonal elements for 
                ! the last RDMEnergyIter iterations.
                tFill_RDM = .true.
                IterLastRDMFill = RDMEnergyIter
            elseif(Iter.eq.NMCyc) then
                ! Last iteration, calculate the diagonal element for the iterations 
                ! since the last time they were included.
                tFill_RDM = .true.
                IterLastRDMFill = mod((Iter - IterRDMStart + 1),RDMEnergyIter)
            endif
        endif

        do j=1,int(TotWalkers,sizeof_int)
            ! N.B. j indicates the number of determinants, not the number
            !      of walkers.

            ! Indicate that the scratch storage used for excitation generation
            ! from the same walker has not been filled (it is filled when we
            ! excite from the first particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            ! If we're not calculating the RDM (or we're calculating some HFSD combination of the 
            ! RDM) this just extracts info from the bit representation like normal.
            ! IterRDMStartCurr and AvSignCurr just come out as 1.0_dp.  
            ! Otherwise, it extracts the Curr info, and calculates the iteration this determinant 
            ! became occupied (IterRDMStartCurr) and the average population during that time 
            ! (AvSignCurr).

            call extract_bit_rep_avsign (CurrentDets(:,j), CurrentH(1:NCurrH,j), &
                                        DetCurr, SignCurr, FlagsCurr, IterRDMStartCurr, &
                                        AvSignCurr, fcimc_excit_gen_store)
            
            ! We only need to find out if determinant is connected to the
            ! reference (so no ex. level above 2 required, 
            ! truncated etc.)
            walkExcitLevel = FindBitExcitLevel (iLutRef, CurrentDets(:,j), &
                                                max_calc_ex_level)
           
             !   if(WalkExcitLevel.eq.0) WRITE(6,*) "AvSignCurr: updating", CurrentH(2,j), AvSignCurr, SignCurr

            if(tRef_Not_HF) then
                walkExcitLevel_toHF = FindBitExcitLevel (iLutHF_true, CurrentDets(:,j), &
                                                max_calc_ex_level)
            else
                walkExcitLevel_toHF = walkExcitLevel
            endif

            ! A general index whose value depends on whether the following option is used.
            if (tHashWalkerList) then
                gen_ind = j
            else
                gen_ind = VecSlot
            end if
            
            ! If this state is in the deterministic space.
            if (tSemiStochastic) then
                if (test_flag(CurrentDets(:,j), flag_deterministic)) then

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
                            CurrentH(2:1+lenof_sign,gen_ind) = AvSignCurr
                            CurrentH(2+lenof_sign:1+2*lenof_sign,gen_ind) = IterRDMStartCurr
                        endif
                        VecSlot = VecSlot + 1
                        cycle
                    end if
                 
                end if
            end if

            if (tTruncInitiator) call CalcParentFlag (j, VecSlot, parent_flags)

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

            ! Should be able to make this function pointer-able
            if (tRegenDiagHEls) then
                ! We are not storing the diagonal hamiltonian elements for 
                ! each particle. Therefore, we need to regenerate them.
                if (DetBitEQ(CurrentDets(:,j), iLutRef, NIfDBO) .and. &
                    (.not.(tHub .and. tReal))) then
                    HDiagCurr = 0
                else
                    if (tHPHF) then
                        HDiagtemp = hphf_diag_helement (DetCurr,CurrentDets(:,j))
                    elseif(tMomInv) then
                        HDiagtemp = MI_diag_helement(DetCurr,CurrentDets(:,j))
                    else
                        HDiagTemp = get_helement (DetCurr, DetCurr, 0)
                    endif
                    HDiagCurr = real(HDiagTemp, dp)-Hii
                endif
            else
                ! HDiags are stored.
                HDiagCurr = CurrentH(1,j)
            endif

            ! Test if we have found a determinant which is lower in E than
            ! the 'root' determinant. Should not happen in an (unrotated)
            ! HF basis.
            if (tFindGroundDet .and. HDiagCurr < 0) then
                call ChangeRefDet (DetCurr)
                exit
            endif

            ! Sum in any energy contribution from the determinant, including 
            ! other parameters, such as excitlevel info.
            ! This is where the projected energy is calculated.
            call SumEContrib (DetCurr, WalkExcitLevel,SignCurr, CurrentDets(:,j), HDiagCurr, 1.0_dp, j)

!            ! If we're filling the RDM, this calculates the explicitly connected singles and doubles.
!            ! Or in the case of HFSD (or some combination of this), it calculates the 
!            ! diagonal elements too - as we need the walkExcitlevel for the diagonal elements as well.
            TempSpawnedPartsInd = 0

            ! Loop over the 'type' of particle. 
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
            !                 --> OR double run
            !                 --> part_type == 1, 2; population sets 1 and 2, both real
            do part_type=1,lenof_sign

                !The logic here is a little messy, but relates to the tSpawn_only_init option.
                !With this, only initiators are allowed to spawn, therefore we are testing whether 
                !the 'type' of walker we are currently considering is an initiator or not.
                !This is determined uniquely for real walkers, where there is only one 'type' of walker,
                !but otherwise we need to test each type independently.
                
                if((lenof_sign.gt.1).and.tSpawn_Only_Init) then
                    if(.not.test_flag (CurrentDets(:,j), flag_parent_initiator(part_type))) cycle
                elseif(tSpawn_Only_Init) then
                    if(.not.tCurr_initiator) cycle
                endif

                ! Loop over all the particles of a given type on the 
                ! determinant. CurrentSign gives number of walkers. Multiply 
                ! up by AvMCExcits if attempting multiple excitations from 
                ! each walker (default 1.0_dp).

                WalkersToSpawn=abs(int(SignCurr(part_type)*AvMCExcits))
                if ((abs(SignCurr(part_type)*AvMCExcits)-real(WalkersToSpawn,dp)).gt.0) then
                    prob_extra_walker=abs(SignCurr(part_type)*AvMCExcits) - real(WalkersToSpawn,dp)
                    r = genrand_real2_dSFMT ()
                    if (prob_extra_walker > r) WalkersToSpawn=WalkersToSpawn+1
                endif
                    

                do p = 1, WalkersToSpawn
                    ! Zero the bit representation, to ensure no extraneous
                    ! data gets through.
                    ilutnJ = 0

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
                            iLutnJ(nOffFlag) = 0
                            
                            ! Is the spawned state in the core space?
                            tInDetermSpace = is_core_state(iLutnJ)

                            ! Is the parent state in the core space?
                            if (test_flag(CurrentDets(:,j), flag_deterministic)) then
                                ! If spawning is from and to the core space, cancel it.
                                if (tInDetermSpace) cycle
                            else
                                if (tInDetermSpace) call set_flag(iLutnJ, flag_deterministic)
                            end if

                            ! If the walker being spawned is spawned from the deterministic space,
                            ! then set the corresponding flag to specify this.
                            if (test_flag(CurrentDets(:,j), flag_deterministic)) &
                                call set_flag(iLutnJ, flag_determ_parent)

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
                        write(iout,"(a,f12.5)",advance='no') &
#else
                        write(iout,"(a,f12.5)",advance='no') &
#endif
                            "SP:", child
                        call write_det(6, nJ, .true.)
                        call neci_flush(iout) 
                    endif

                    ! Children have been chosen to be spawned.
                    if (any(child /= 0)) then

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
                if (.not. test_flag(CurrentDets(:,j), flag_deterministic)) then
                    call walker_death (iter_data, DetCurr, &
                                       CurrentDets(:,j), HDiagCurr, SignCurr, &
                                       AvSignCurr, IterRDMStartCurr, VecSlot, j, WalkExcitLevel)
                else
                    CurrentDets(:,gen_ind) = CurrentDets(:,j)
                    if (tFillingStochRDMonFly) then
                        CurrentH(2:1+lenof_sign,gen_ind) = AvSignCurr
                        CurrentH(2+lenof_sign:1+2*lenof_sign,gen_ind) = IterRDMStartCurr
                        ! NB We never overwrite the deterministic states, so move the next spawning slot
                        ! in CurrentDets to the next state.
                    endif
                    VecSlot = VecSlot + 1
                end if
            else
                call walker_death (iter_data, DetCurr, &
                                   CurrentDets(:,j), HDiagCurr, SignCurr, &
                                   AvSignCurr, IterRDMStartCurr, VecSlot, j, WalkExcitLevel)
            end if

            if(tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) &
                    & call fill_rdm_diag_currdet(CurrentDets(:,gen_ind), DetCurr, SignCurr, &
                    & walkExcitLevel_toHF, test_flag(CurrentDets(:,j), flag_deterministic))  

        enddo ! Loop over determinants.
        IFDEBUGTHEN(FCIMCDebug,2) 
            write(iout,*) 'Finished loop over determinants'
            if(tHashWalkerList) then
                write(iout,*) "Holes in list: ",iEndFreeSlot
            endif
        ENDIFDEBUG

        if(tHashWalkerList) then
            !With this algorithm, the determinants do not move, and therefore TotWalkersNew is simply equal
            !to TotWalkers
            TotWalkersNew=int(TotWalkers,sizeof_int)
        else
            ! Since VecSlot holds the next vacant slot in the array, TotWalkers
            ! should be one less than this. TotWalkersNew is now the number of particles
            ! in the main array, before annihilation
            TotWalkersNew = VecSlot - 1
        endif

        ! For semi-stochastic calculations only: Gather together the parts of the deterministic vector stored
        ! on each processor, and then perform the multiplication of the exact projector on this vector.
        if (tSemiStochastic) call deterministic_projection()

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

        ! Update iteration data
        iter_data%update_growth = iter_data%update_growth + iter_data%nborn &
                                - iter_data%ndied - iter_data%nannihil &
                                - iter_data%naborted - iter_data%nremoved
        iter_data%update_iters = iter_data%update_iters + 1

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

    subroutine end_iter_stats (TotWalkersNew)

        integer, intent(in) :: TotWalkersNew
        HElement_t :: delta(inum_runs)
        integer :: proc, pos, i, k
        real(dp) :: sgn(lenof_sign)

        ! SumWalkersCyc calculates the total number of walkers over an update
        ! cycle on each process.
#ifdef __CMPLX
        SumWalkersCyc = SumWalkersCyc + sum(TotParts)
#else
        SumWalkersCyc = SumWalkersCyc + TotParts
#endif
        ! Write initiator histograms if on the correct iteration.
        ! Why is this done here - before annihilation!
        if ((tHistInitPops .and. mod(iter, HistInitPopsIter) == 0) &
            .or. tPrintHighPop) then
            call FindHighPopDet (TotWalkersNew)
            root_write(iout,'(a)') 'Writing out the spread of the initiator &
                                &determinant populations.'
            call WriteInitPops (iter + PreviousCycles)
        endif

        ! Update the sum for the projected-energy denominatior if projecting
        ! onto a linear combination of determinants.
        ! Why is this done here - before annihilation!
        if (proje_spatial .and. proje_linear_comb .and. nproje_sum > 1) then
            do i = 1, nproje_sum
                proc = DetermineDetNode (proje_ref_dets(:,i), 0)
                if (iProcIndex == proc) then
                    pos = binary_search (CurrentDets(:,1:TotWalkers), &
                                         proje_ref_iluts(:,i), NIfD+1)
                    if (pos > 0) then
                        call extract_sign (CurrentDets(:,pos), sgn)
                        do k=1, inum_runs
                            delta(k) = ARR_RE_OR_CPLX(sgn,k) * proje_ref_coeffs(i)
                            cyc_proje_denominator(k) = cyc_proje_denominator(k) + delta(k)
                            sum_proje_denominator(k) = sum_proje_denominator(k) + delta(k)
                        enddo
                    endif
                endif
            enddo
        endif

    end subroutine end_iter_stats

    subroutine new_child_stats_normal (iter_data, iLutI, nJ, iLutJ, ic, &
                                       walkExLevel, child, parent_flags, &
                                       part_type)

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        real(dp), dimension(lenof_sign), intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(n_int) :: iUnused

        ! Write out some debugging information if asked
        IFDEBUG(FCIMCDebug,3) then
            if(lenof_sign.eq.2) then
                write(iout,"(A,2f10.5,A)", advance='no') &
                               "Creating ", child(1:lenof_sign), " particles: "
            else
                write(iout,"(A,f10.5,A)",advance='no') &
                                         "Creating ", child(1), " particles: "
            endif
            write(iout,"(A,2I4,A)",advance='no') &
                                      "Parent flag: ", parent_flags, part_type
            call writebitdet (iout, ilutJ, .true.)
            call neci_flush(iout)
        endif

        ! Count the number of children born
        NoBorn = NoBorn + sum(abs(child))
        iter_data%nborn = iter_data%nborn + abs(child)

        if (ic == 1) SpawnFromSing = SpawnFromSing + int(sum(abs(child)))

        ! Count particle blooms, and their sources
        if (sum(abs(child)) > INitiatorWalkNo) then
            bloom_count(ic) = bloom_count(ic) + 1
            bloom_sizes(ic) = max(real(sum(abs(child)), dp), bloom_sizes(ic))
        end if

        ! Histogram the excitation levels as required
        if (tHistExcitToFrom) &
            call add_hist_excit_tofrom(ilutI, ilutJ, child)

        ! Avoid compiler warnings
        iUnused = iLutI(0); iUnused = iLutJ(0)

    end subroutine

    subroutine create_particle (nJ, iLutJ, child, parent_flags, part_type, iLutI, SignI, &
                                WalkerNumber, RDMBiasFacCurr, WalkersToSpawn)

        ! Create a child in the spawned particles arrays. We spawn particles
        ! into a separate array, but non-contiguously. The processor that the
        ! newly-spawned particle is going to be sent to has to be determined,
        ! and then it will be put into the appropriate element determined by
        ! ValidSpawnedList

        use bit_rep_data, only: flag_bit_offset

        integer, intent(in) :: nJ(nel)
        integer(kind=n_int), intent(in) :: iLutJ(0:niftot)
        real(dp), dimension(lenof_sign), intent(in) :: child
        integer(kind=n_int), intent(in), optional :: iLutI(0:niftot)
        real(dp), dimension(lenof_sign), intent(in), optional :: SignI
        integer, intent(in) :: parent_flags
        integer, intent(in), optional :: WalkerNumber
        integer, intent(in), optional :: WalkersToSpawn
        ! 'type' of the particle - i.e. real/imag
        integer, intent(in) :: part_type
        real(dp) , intent(in) , optional :: RDMBiasFacCurr
        integer :: proc, flags, j
        logical :: parent_init

        proc = DetermineDetNode(nJ,0)    ! 0 -> nNodes-1)

        ! We need to include any flags set both from the parent and from the
        ! spawning steps. No we don't! - ghb
        ! This is highly highly yucky and needs cleaning up.
        ! Potentially, ilutJ can be given the flag of its parent in 
        ! FindExcitBitDet routine. I don't think it should be.
        ! To make things more confusing, this only happens for non-HPHF/CSF 
        ! runs.
        ! Things are even more confusing given the fact that CCMC is using this routine
        ! Who know whether they require the flag there or not...
        !TODO: CLEAN THIS UP. Make it clear, and transparent, with one way to change the
        ! flag. Otherwise, this will trip many people up in the future.
        flags = ior(parent_flags, extract_flags(ilutJ))

!        WRITE(6,*) 'Encoding',iLutJ
!        WRITE(6,*) 'To position',ValidSpawnedList(proc)
!        WRITE(6,*) 'Parent',iLutI
        call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), iLutJ, &
                            child, flags)

        ! For double run, we're going to try only including off-diagonal elements from
        ! successful spawning attempts within part_type=1.  Later, we should include both,
        ! but for now, this may be simplest.
        IF(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)) then
            call store_parent_with_spawned(RDMBiasFacCurr, WalkerNumber, iLutI, WalkersToSpawn, iLutJ, proc, part_type)
        endif
        !We are spawning from iLutI to SpawnedParts(:,ValidSpawnedList(proc)).
        !We want to store the parent (D_i) with the spawned child (D_j) so that we can
        !add in Ci.Cj to the RDM later on.
        !The parent is NIfDBO integers long, and stored in the second part of the SpawnedParts array 
        !from NIfTot+1 -> NIfTot+1 + NIfDBO.

        IF(lenof_sign.eq.2) THEN
            !With complex walkers, things are a little more tricky.
            !We want to transfer the flag for all particles created (both real and imag)
            !from the specific type of parent particle.
            !This can mean real walker flags being transfered to imaginary children and
            !vice versa.
            !This is unneccesary for real walkers.
            !Test the specific flag corresponding to the parent, of type 'part_type'
            parent_init=test_flag(SpawnedParts(:,ValidSpawnedList(proc)),flag_parent_initiator(part_type))
            !Assign this flag to all spawned children
            do j=1,lenof_sign
                if(child(j).ne.0) then
                    call set_flag(SpawnedParts(:,ValidSpawnedList(proc)),flag_parent_initiator(j),parent_init)
                endif
            enddo
        ENDIF

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        
        ! Sum the number of created children to use in acceptance ratio.
        acceptances = acceptances + int(sum(abs(child)), kind(acceptances))

    end subroutine

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

    subroutine CalcParentFlag(j, VecSlot, parent_flags)
!In the CurrentDets array, the flag at NIfTot refers to whether that determinant *itself* is an initiator or not.    
!We need to decide if this willchange due to the determinant acquiring a certain population, or its population dropping
!below the threshold.
!The CurrentDets(:,j) is the determinant we are currently spawning from, so this determines the ParentInitiator flag
!which is passed to the SpawnedDets array and refers to whether or not the walkers *parent* is an initiator or not.
!A flag of 0 means the determinant is an initiator, and 1 it is a non-initiator.
        integer, intent(in) :: j, VecSlot
        integer, intent(out) :: parent_flags
        real(dp), dimension(lenof_sign) :: CurrentSign
        integer :: part_type
        real(dp) :: init_thresh
        logical :: tDetinCAS, parent_init

        call extract_sign (CurrentDets(:,j), CurrentSign)

        init_thresh = InitiatorWalkNo

        tcurr_initiator = .false.
        do part_type=1,lenof_sign
            ! By default, the parent_flags are the flags of the parent...
            parent_init = test_flag (CurrentDets(:,j), flag_parent_initiator(part_type))

            ! The default path through this section makes no changes, leaving
            ! the initiator status of each parent unchanged.  If 
            ! tAddToInitiator is set, then the state of the parent may change.
            if (tAddToInitiator) then

                if (.not. parent_init) then
                    ! Determinant wasn't previously initiator 
                    ! - want to test if it has now got a large enough 
                    !   population to become an initiator.
                    if (abs(CurrentSign(part_type)) > init_thresh) then
                        parent_init = .true.
                        NoAddedInitiators = NoAddedInitiators + 1
                        if (tSpawnSpatialInit) &
                            call add_initiator_list (CurrentDets(:,j))
                    endif
                elseif (tRetestAddToInit) then
                    ! The source determinant is already an initiator.            
                    ! If tRetestAddToInit is on, the determinants become 
                    ! non-initiators again if their population falls below 
                    ! n_add (this is on by default).
                    tDetInCas = .false.
                    if (tTruncCAS) &
                        tDetInCas = TestIfDetInCASBit (CurrentDets(0:NIfD,j))
                   
                    ! If det. in fixed initiator space, or is the HF det, or it
                    ! is in the deterministic space, then it must remain an initiator.
                    if (.not. tDetInCas .and. &
                        .not. (DetBitEQ(CurrentDets(:,j), iLutHF, NIfDBO)) &
                        .and. .not. test_flag(CurrentDets(:,j), flag_deterministic) &
                        .and. abs(CurrentSign(part_type)) <= init_thresh &
                        .and. .not. test_flag(CurrentDets(:,j), &
                        flag_make_initiator(part_type))) then
                        ! Population has fallen too low. Initiator status 
                        ! removed.
                        parent_init = .false.
                        NoAddedInitiators = NoAddedInitiators - 1
                        if (tSpawnSpatialInit) &
                            call rm_initiator_list (CurrentDets(:,j))
                    endif
                endif
            endif

            ! Update counters as required.
            if (parent_init) then
                NoInitDets = NoInitDets + 1
                NoInitWalk = NoInitWalk + abs(CurrentSign(part_type))
            else
                NoNonInitDets = NoNonInitDets + 1
                NoNonInitWalk = NoNonInitWalk + abs(CurrentSign(part_type))
            endif

            ! Update the parent flag as required.
            call set_flag (CurrentDets(:,j), flag_parent_initiator(part_type), &
                           parent_init)
            
            if(parent_init) then                           
                tcurr_initiator = .true.
            endif

        enddo

        ! Store this flag for use in the spawning routines...
        parent_flags = extract_flags (CurrentDets(:,j))

        ! We don't want the deterministic flag to be set in parent_flags, as this
        ! would set the same flag in the child in create_particle, which we don't
        ! want in general.
        parent_flags = ibclr(parent_flags, flag_deterministic)

        ! We don't want the child to have trial or connected flags.
        parent_flags = ibclr(parent_flags, flag_trial)
        parent_flags = ibclr(parent_flags, flag_connected)

        if ((tHistInitPops .and. mod(iter, histInitPopsIter) == 0) &
            .or. tPrintHighPop) then
             call HistInitPopulations (CurrentSign(1), VecSlot)
        endif

    end subroutine CalcParentFlag


    SUBROUTINE HistInitPopulations(RealSignCurr,VecSlot)
        USE FciMCLoggingMOD, only : InitBinMin,InitBinIter
        INTEGER , INTENT(IN) :: VecSlot
        REAL(dp), INTENT(IN) :: RealSignCurr
        INTEGER :: InitBinNo

        IF(ABS(RealSignCurr).gt.InitiatorWalkNo) THEN
!Just summing in those determinants which are initiators. 

!Need to figure out which bin to put them in though.
            IF(RealSignCurr.lt.0.0_dp) THEN
                InitBinNo=(FLOOR(((log(ABS(RealSignCurr)))-InitBinMin)/InitBinIter))+1
                IF((InitBinNo.ge.1).and.(InitBinNo.le.25000)) THEN
                    HistInitPops(1,InitBinNo)=HistInitPops(1,InitBinNo)+1
                ELSE
                    CALL Stop_All('HistInitPopulations','Trying to histogram outside the range of the bins.')
                ENDIF

            ELSE
                InitBinNo=(FLOOR(((log(RealSignCurr))-InitBinMin)/InitBinIter))+1

                IF((InitBinNo.ge.1).and.(InitBinNo.le.25000)) THEN
                    HistInitPops(2,InitBinNo)=HistInitPops(2,InitBinNo)+1
                ELSE
                    CALL Stop_All('HistInitPopulations','Trying to histogram outside the range of the bins.')
                ENDIF
            ENDIF
        ENDIF

        IF(RealSignCurr.lt.MaxInitPopNeg) THEN
            MaxInitPopNeg=RealSignCurr
            HighPopNeg=VecSlot
        ENDIF
        IF(RealSignCurr.gt.MaxInitPopPos) THEN
            MaxInitPopPos=RealSignCurr
            HighPopPos=VecSlot
        ENDIF


    END SUBROUTINE HistInitPopulations


    SUBROUTINE FindHighPopDet(TotWalkersNew)
!Found the highest population on each processor, need to find out which of these has the highest of all.
        INTEGER(KIND=n_int) :: DetPos(0:NIfTot),DetNeg(0:NIfTot)
        INTEGER :: TotWalkersNew,ProcBCastNeg,ProcBCastPos
        !integer(int32) :: HighPopInNeg(2),HighPopInPos(2),HighPopoutNeg(2),HighPopoutPos(2)
        real(dp) :: HighPopInNeg(2),HighPopInPos(2),HighPopoutNeg(2),HighPopoutPos(2)
        real(dp) :: TempSign(lenof_sign)

!        WRITE(iout,*) 'HighPopPos',HighPopPos
!        WRITE(iout,*) 'CurrentSign(HighPopPos)',CurrentSign(HighPopPos)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopNeg),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF

        HighPopInNeg(1)=TempSign(1)
        HighPopInNeg(2)=real(iProcIndex,dp)

        CALL MPIAllReduceDatatype(HighPopinNeg,1,MPI_MINLOC,MPI_2DOUBLE_PRECISION,HighPopoutNeg)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopPos),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF
        
        HighPopInPos(1)=TempSign(1)
        HighPopInPos(2)=real(iProcIndex,dp)

        CALL MPIAllReduceDatatype(HighPopinPos,1,MPI_MAXLOC,MPI_2DOUBLE_PRECISION,HighPopoutPos)

        ! Now, on all processors, HighPopoutPos(1) is the highest positive 
        ! population, and HighPopoutNeg(1) is the highest negative population.
        ! HighPopoutPos(2) is the processor the highest population came from.

        IF(real(iProcIndex,dp).eq.HighPopOutNeg(2)) DetNeg(:)=CurrentDets(:,HighPopNeg)
        IF(real(iProcIndex,dp).eq.HighPopOutPos(2)) DetPos(:)=CurrentDets(:,HighPopPos)

        ! This is a horrible hack, because the process argument should be of 
        ! type 'integer' - whatever that is, but the highpopoutneg is 
        ! explicitly an int(4), so that it works with MPI_2INTEGER. Because
        ! of the explicit interfaces, we need to do this.
        CALL MPIBcast(DetNeg,NIfTot+1,int(HighPopOutNeg(2)))
        CALL MPIBcast(DetPos,NIfTot+1,int(HighPopOutPos(2)))

        if (iProcIndex == 0) then
            write (iout, '(a, f12.5, a)') 'The most highly populated determinant &
                                  & with the opposite sign to the HF has ', &
                                  HighPopoutNeg(1), ' walkers.'
            call WriteBitDet (iout, DetNeg, .true.)

            write (iout, '(a, f12.5, a9)') 'The most highly populated determinant &
                                  & with the same sign as the HF has ', &
                                  HighPopoutPos(1), ' walkers.'
            call WriteBitDet (iout, DetPos, .true.)
        endif

        tPrintHighPop=.false.


    END SUBROUTINE FindHighPopDet

    

    subroutine init_fcimc_fn_pointers ()

        ! Select the excitation generator
        if (tHPHF) then
            generate_excitation => gen_hphf_excit
        elseif (tCSF) then
            generate_excitation => gen_csf_excit
        elseif (tPickVirtUniform) then
            generate_excitation => gen_rand_excit3
        elseif (tMomInv) then
            generate_excitation => gen_MI_excit
        elseif (tGenHelWeighted) then
            generate_excitation => gen_excit_hel_weighted
        elseif (tGen_4ind_weighted) then
            generate_excitation => gen_excit_4ind_weighted
        elseif (tGen_4ind_reverse) then
            generate_excitation => gen_excit_4ind_reverse
        else
            generate_excitation => gen_rand_excit
        endif

        ! In the main loop, we only need to find out if a determinant is
        ! connected to the reference det or not (so no ex. level above 2 is
        ! required). Except in some cases where we need to know the maximum
        ! excitation level
        if (tTruncSpace .or. tHistSpawn .or. tCalcFCIMCPsi .or. &
            tHistHamil) then
            max_calc_ex_level = nel
        else
            max_calc_ex_level = 2
        endif

        ! How many children should we spawn given an excitation?
        if ((tTruncCas .and. (.not. tTruncInitiator)) .or. tTruncSpace .or. &
            tPartFreezeCore .or. tPartFreezeVirt .or. tFixLz .or. &
            (tUEG .and. .not. tLatticeGens) .or. tTruncNOpen) then
            if (tHPHF .or. tCSF .or. tMomInv .or. tSemiStochastic) then
                attempt_create => attempt_create_trunc_spawn
            else
                attempt_create => att_create_trunc_spawn_enc
            endif
        else
            attempt_create => attempt_create_normal
        endif

        ! In attempt create, how should we evaluate the off diagonal matrix
        ! elements between a parent and its (potentially) spawned offspring?
        if (tCSF) then
            get_spawn_helement => get_csf_helement
        elseif (tHPHF) then
            if (tGenMatHEL) then
                get_spawn_helement => hphf_spawn_sign
            else
                get_spawn_helement => hphf_off_diag_helement_spawn
            endif
        elseif (tMomInv) then
            if (tGenMatHEl) then
                get_spawn_helement => MI_spawn_sign
            else
                get_spawn_helement => MI_off_diag_helement_spawn
            endif
        else
            get_spawn_helement => get_helement_det_only
        endif

        ! Once we have generated the children, do we need to encode them?
        if (.not. (tCSF .or. tHPHF .or. tMomInv .or. tGen_4ind_weighted)) then
            encode_child => FindExcitBitDet
        else
            encode_child => null_encode_child
        endif

        ! What message should we display for a particle bloom?
        if (tAddToInitiator) then
            bloom_warn_string = '("Bloom of more than n_add on ", a, " excit: &
                                &A max of ", f10.2, " particles created. ", &
                                &i8, " blooms occurred.")'
        else
            ! Use this variable to store the bloom cutoff level.
            InitiatorWalkNo = 25.0_dp
            bloom_warn_string = '("Bloom of more than 25 on ", a, " excit: &
                                &A max of ", f10.2, " particles created. ", &
                                &i8, " blooms occurred.")'
        endif
        bloom_max = 0

        ! Perform the correct statistics on new child particles
        if (tHistHamil) then
            new_child_stats => new_child_stats_hist_hamil
        else
            new_child_stats => new_child_stats_normal
        endif

        attempt_die => attempt_die_normal

        extract_bit_rep_avsign => extract_bit_rep_avsign_no_rdm

        if(tHF_Ref_Explicit.or.tHF_S_D.or.tHF_S_D_Ref) then
            fill_rdm_diag_currdet => fill_rdm_diag_currdet_hfsd
        else
            fill_rdm_diag_currdet => fill_rdm_diag_currdet_norm
        endif

    end subroutine


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
        
    function attempt_create_trunc_spawn (DetCurr,&
                                         iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type, &
                                         AvSignCurr, RDMBiasFacCurr) result(child)
        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t, intent(in) :: HElGen

        if (CheckAllowedTruncSpawn (walkExcitLevel, nJ, iLutnJ, IC)) then
            child = attempt_create_normal (DetCurr, &
                               iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, ic, ex, &
                               tParity, walkExcitLevel, part_type, AvSignCurr, RDMBiasFacCurr)
        else
            child = 0
        endif
    end function

!Decide whether to spawn a particle at nJ from DetCurr. (bit strings iLutnJ and iLutCurr respectively).  
!  ic and ex specify the excitation of nJ from DetCurr, along with the sign change tParity.
!  part_type:           Is the parent real (1) or imaginary (2)
!  wSign:               wSign gives the sign of the particle we are trying to spawn from
!                          if part_type is 1, then it will only use wsign(1)
!                                          2,                       wsign(2)
!                       Only the sign, not magnitude is used.
!  prob:                prob is the generation probability of the excitation in order to unbias.
!                       The probability of spawning is divided by prob to do this.
!  HElGen:              If the matrix element has already been calculated, it is sent in here.
!  get_spawn_helement:  A function pointer for looking up or calculating the relevant matrix element.
!  walkExcitLevel:      Is Unused
! 
!  child:      A lenof_sign array containing the particles spawned.
    function att_create_trunc_spawn_enc (DetCurr,&
                                         iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type, &
                                         AvSignCurr,RDMBiasFacCurr) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t , intent(in) :: HElGen

        call EncodeBitDet (nJ, iLutnJ)
        if (CheckAllowedTruncSpawn (walkExcitLevel, nJ, iLutnJ, IC)) then
            child = attempt_create_normal (DetCurr, &
                               iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, ic, ex, &
                               tParity, walkExcitLevel, part_type, AvSignCurr, RDMBiasFacCurr)
        else
            child = 0
        endif
    end function

    function attempt_create_normal (DetCurr, iLutCurr, &
                                    RealwSign, nJ, iLutnJ, prob, HElGen, ic, ex, tParity,&
                                    walkExcitLevel, part_type, AvSignCurr, RDMBiasFacCurr) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! 1 = Real parent particle, 2 = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t , intent(in) :: HElGen
        character(*), parameter :: this_routine = 'attempt_create_normal'

        real(dp) :: rat, r, walkerweight, pSpawn, nSpawn, MatEl, p_spawn_rdmfac
        integer :: extracreate, tgt_cpt, component, i, iUnused
        integer :: TargetExcitLevel
        logical :: tRealSpawning
        HElement_t :: rh

        ! Just in case
        child = 0

        ! If each walker does not have exactly one spawning attempt
        ! (if AvMCExcits /= 1.0_dp) then the probability of an excitation
        ! having been chosen, prob, must be altered accordingly.
        prob = prob * AvMCExcits

        ! In the case of using HPHF, and when tGenMatHEl is on, the matrix
        ! element is calculated at the time of the excitation generation, 
        ! and returned in HElGen. In this case, get_spawn_helement simply
        ! returns HElGen, rather than recomputing the matrix element.
        rh = get_spawn_helement (DetCurr, nJ, iLutCurr, iLutnJ, ic, ex, &
                                 tParity, HElGen)

        !write(6,*) 'p,rh', prob, rh

        ! The following is useful for debugging the contributions of single
        ! excitations, and double excitations of spin-paired/opposite
        ! electron pairs to the value of tau.
!        if (ic == 2) then
!            if (G1(ex(1,1))%Ms /= G1(ex(1,2))%Ms) then
!                write(6,*) 'OPP', rh, prob
!            else
!                write(6,*) 'SAM', rh, prob
!            end if
!        else
!            write(6,*) 'IC1', rh, prob
!        end if

        ! Are we doing real spawning?
        
        tRealSpawning = .false.
        if (tAllRealCoeff) then
            tRealSpawning = .true.
        elseif (tRealCoeffByExcitLevel) then
            TargetExcitLevel = FindBitExcitLevel (iLutRef, iLutnJ)
            if (TargetExcitLevel <= RealCoeffExcitThresh) &
                tRealSpawning = .true.
        endif

        ! We actually want to calculate Hji - take the complex conjugate, 
        ! rather than swap around DetCurr and nJ.
#ifdef __CMPLX
        rh = conjg(rh)
#endif
        
        ! Spawn to real and imaginary particles. Note that spawning from
        ! imaginary parent particles has slightly different rules:
        !       - Attempt to spawn REAL walkers with prob +AIMAG(Hij)/P
        !       - Attempt to spawn IMAG walkers with prob -REAL(Hij)/P
        do tgt_cpt = 1, (lenof_sign/inum_runs)

            ! Real, single run:    inum_runs=1, lenof_sign=1 --> 1 loop
            ! Complex, single run: inum_runs=1, lenof_sign=2 --> 2 loops
            ! Real, double run:    inum_runs=2, lenof_sign=2 --> 1 loop

            ! For spawning from imaginary particles, we cross-match the 
            ! real/imaginary matrix-elements/target-particles.
            component = tgt_cpt

            if ((part_type.eq.2).and.(inum_runs.eq.1)) component = 3 - tgt_cpt

            ! Get the correct part of the matrix element
            walkerweight = sign(1.0_dp, RealwSign(part_type))
            if (component == 1) then
                MatEl = real(rh, dp)
            else
#ifdef __CMPLX
                MatEl = real(aimag(rh), dp)
                ! n.b. In this case, spawning is of opposite sign.
                if (part_type == 2) walkerweight = -walkerweight
#else
                call stop_all('attempt_create_normal', &
                        & "ERROR: We shouldn't reach this part of the code unless complex calc")
#endif
            end if

            nSpawn = - tau * MatEl * walkerweight / prob
            
            ! n.b. if we ever end up with |walkerweight| /= 1, then this
            !      will need to ffed further through.
            if (tSearchTau) &
                call log_spawn_magnitude (ic, ex, matel, prob)

            ! Keep track of the biggest spawn this cycle
            max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)
            
            if (tRealSpawning) then
                ! Continuous spawning. Add in acceptance probabilities.
                
                if (tRealSpawnCutoff .and. &
                    abs(nSpawn) < RealSpawnCutoff) then
                    p_spawn_rdmfac=abs(nSpawn)/RealSpawnCutoff
                    nSpawn = RealSpawnCutoff &
                           * stochastic_round (nSpawn / RealSpawnCutoff)
               else
                    p_spawn_rdmfac=1.0_dp !The acceptance probability of some kind of child was equal to 1
               endif
            else
                if(abs(nSpawn).ge.1) then
                    p_spawn_rdmfac=1.0_dp !We were certain to create a child here.
                    ! This is the special case whereby if P_spawn(j | i) > 1, 
                    ! then we will definitely spawn from i->j.
                    ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                    ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                    ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                else
                    p_spawn_rdmfac=abs(nSpawn)
                endif
                
                ! How many children should we spawn?

                ! And round this to an integer in the usual way
                ! HACK: To use the same number of random numbers for the tests.
                if (nspawn - real(int(nspawn),dp) == 0.0_dp) r = genrand_real2_dSFMT()
                nSpawn = real(stochastic_round (nSpawn), dp)
            endif
#ifdef __CMPLX
            ! And create the parcticles
            child(tgt_cpt) = nSpawn
#else
            ! And create the parcticles
            child(part_type) = nSpawn
#endif
        enddo
       
        if(tFillingStochRDMonFly) then
            if((child(part_type).ne.0).and.(.not.tHF_Ref_Explicit).and.(part_type.eq.1)) then
                !Only add in contributions for spawning events within population 1
                !(Otherwise it becomes tricky in annihilation as spawnedparents doesn't tell you which population the event came from at present)
                call calc_rdmbiasfac(p_spawn_rdmfac, prob, AvSignCurr, realwSign, RDMBiasFacCurr) 
            else
                RDMBiasFacCurr = 0.0_dp
            endif
        else
            ! Not filling the RDM stochastically, bias is zero.
            RDMBiasFacCurr = 0.0_dp
        endif

        ! Avoid compiler warnings
        iUnused = walkExcitLevel

    end function

    ! Depreciated function - only used for CCMC
    !
    ! This function tells us whether we should create a child particle on nJ,
    ! from a parent particle on DetCurr with sign WSign, created with 
    ! probability Prob. It returns zero if we are not going to create a child,
    ! or -1/+1 if we are to create a child, giving the sign of the new
    ! particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,iLutCurr,RealWSign,nJ,iLutnJ,Prob,IC,Ex,tParity)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        use LoggingData, only : CCMCDebug
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,ExtraCreate,Ex(2,2),Bin
        INTEGER(KIND=n_int) :: iLutCurr(0:NIfTot),iLutnJ(0:NIfTot)
        logical :: tParity
        real(dp) :: Prob, r, rat, RealAttemptCreatePar
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        HElement_t :: rh

        ! If each walker does not have exactly one spawning attempt
        ! (if AvMCExcits /= 1.0_dp) then the probability of an excitation
        ! having been chosen, prob, must be altered accordingly.
        prob = prob * AvMCExcits
            

!Calculate off diagonal hamiltonian matrix element between determinants
!        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        IF(tHPHF) THEN
            IF(tGenMatHEl) THEN
!The prob is actually prob/HEl, since the matrix element was generated at the same time as the excitation

                write(iout,*) "AttemptCreatePar is depreciated, and not compatible with HPHF - use attempt_create_normal"
                call stop_all("AttemptCreatePar","This is a depreciated routine.")

!                rat=Tau/abs(Prob)

!                rh=Prob ! to get the signs right for later on.
!                WRITE(iout,*) Prob, DetCurr(:),"***",nJ(:)
!                WRITE(iout,*) "******"
!                CALL HPHFGetOffDiagHElement(DetCurr,nJ,iLutCurr,iLutnJ,rh)

            ELSE
                ! The IC given doesn't really matter. It just needs to know 
                ! whether it is a diagonal or off-diagonal matrix element.
                ! However, the excitation generator can generate the same 
                ! HPHF again. If this is done, the routine will send the 
                ! matrix element back as zero.
                rh = hphf_off_diag_helement (DetCurr, nJ, iLutCurr, iLutnJ)

                ! Divide by the probability of creating the excitation to 
                ! negate the fact that we are only creating a few determinants
                rat=Tau*abs(rh)/Prob
!                WRITE(iout,*) Prob/rh, DetCurr(:),"***",nJ(:)
!                WRITE(iout,*) "******"

            ENDIF
        elseif(tMomInv) then

            write(iout,*) "AttemptCreatePar is a depreciated routine, and is &
                          &not compatible with MomInv - use &
                          &attempt_create_normal"
            call stop_all("AttemptCreatePar","This is a depreciated routine.")
        ELSE
!Normal determinant spawn

            rh = get_helement (DetCurr, nJ, IC, Ex, tParity)
            !WRITE(iout,*) rh

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
            rat=Tau*abs(rh)/Prob
        ENDIF
        IF(CCMCDebug.gt.5) WRITE(iout,*) "Connection H-element to spawnee:",rh
!        CALL IsSymAllowedExcit(DetCurr,nJ,IC,Ex,SymAllowed) 
!        IF((.not.SymAllowed).and.(abs(rh).gt.0.0_dp)) THEN
!            WRITE(17,*) rh
!        ENDIF

!        rhcheck=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!        IF(rh.ne.rhcheck) THEN
!            WRITE(iout,*) "DetCurr: ",DetCurr(:)
!            WRITE(iout,*) "nJ: ",nJ(:)
!            WRITE(iout,*) "EX: ",Ex(1,:),Ex(2,:)
!            WRITE(iout,*) "tParity: ",tParity
!            STOP
!        ENDIF

!        IF(abs(rh).le.HEpsilon) THEN
!            AttemptCreatePar=0
!            RETURN
!        ENDIF


!If probability is > 1, then we can just create multiple children at the chosen determinant
        ExtraCreate=INT(rat)
        rat=rat-REAL(ExtraCreate,dp)

!Stochastically choose whether to create or not according to ranlux 
        r = genrand_real2_dSFMT() 
        IF(rat.gt.r) THEN
!            IF(Iter.eq.18925) THEN
!                WRITE(iout,*) "Created",rh,rat
!            ENDIF

!Child is created - what sign is it?
            IF(RealwSign(1).gt.0) THEN
!Parent particle is positive
                IF(real(rh,dp).gt.0.0_dp) THEN
                    AttemptCreatePar=-1     !-ve walker created
                ELSE
                    AttemptCreatePar=1      !+ve walker created
                ENDIF

            ELSE
!Parent particle is negative
                IF(real(rh,dp).gt.0.0_dp) THEN
                    AttemptCreatePar=1      !+ve walker created
                ELSE
                    AttemptCreatePar=-1     !-ve walker created
                ENDIF
            ENDIF

        ELSE
!No child particle created
!            IF(Iter.eq.18925) THEN
!                WRITE(iout,*) "Not Created",rh,rat
!            ENDIF
            AttemptCreatePar=0
        ENDIF

        IF(ExtraCreate.ne.0) THEN
!Need to include the definitely create additional particles from a initial probability > 1

            IF(AttemptCreatePar.lt.0) THEN
!In this case particles are negative
                AttemptCreatePar=AttemptCreatePar-ExtraCreate
            ELSEIF(AttemptCreatePar.gt.0) THEN
!Include extra positive particles
                AttemptCreatePar=AttemptCreatePar+ExtraCreate
            ELSEIF(AttemptCreatePar.eq.0) THEN
!No particles were stochastically created, but some particles are still definatly created - we need to determinant their sign...
                if (RealwSign(1) > 0) then
                    if (real(rh,dp) > 0.0_dp) then
                        AttemptCreatePar=-ExtraCreate    !Additional particles are negative
                    ELSE
                        AttemptCreatePar=ExtraCreate       !Additional particles are positive
                    ENDIF
                ELSE
                    IF(real(rh,dp).gt.0.0_dp) THEN
                        AttemptCreatePar=ExtraCreate
                    ELSE
                        AttemptCreatePar=-ExtraCreate
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        RealAttemptCreatePar=real(AttemptCreatePar,dp)
        AttemptCreatePar=transfer(RealAttemptCreatePar, AttemptCreatePar)

        
!We know we want to create a particle. Return the bit-representation of this particle (if we have not already got it)
        IF(.not.tHPHF.and.AttemptCreatePar.ne.0) THEN
            if (tCSF) then
                ! This makes sure that the Yamanouchi symbol is correct. It
                ! also makes it work if we have tTruncateCSF on, and ex would
                ! therefore leave all singles as beta, when we switch to dets.
                call EncodeBitDet (nJ, iLutnJ)
            else
                call FindExcitBitDet(iLutCurr,iLutnJ,IC,Ex)
            endif
        ENDIF

!        IF(AttemptCreatePar.ne.0) THEN
!            WRITE(iout,"(A,F15.5,I5,G25.15,I8,G25.15)") "Outwards ", rat,ExtraCreate,real(rh),Prob
!        ENDIF

        IF(tHistEnergies) THEN
!First histogram off-diagonal matrix elements.
            Bin=INT((real(rh,dp)+OffDiagMax)/OffDiagBinRange)+1
            IF(Bin.le.0.or.Bin.gt.iOffDiagNoBins) THEN
                CALL Stop_All("AttemptCreatePar","Trying to histogram off-diagonal matrix elements, "&
                & //"but outside histogram array bounds.")
            ENDIF
            IF(IC.eq.1) THEN
                SinglesAttemptHist(Bin)=SinglesAttemptHist(Bin)+(Tau/Prob)
                IF(RealAttemptCreatePar.ne.0) THEN
                    SinglesHist(Bin)=SinglesHist(Bin)+abs(RealAttemptCreatePar)
                    IF(BRR(Ex(1,1)).le.NEl) THEN
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistOccOcc(Bin)=SinglesHistOccOcc(Bin)+abs(RealAttemptCreatePar)
                        ELSE
                            SinglesHistOccVirt(Bin)=SinglesHistOccVirt(Bin)+abs(RealAttemptCreatePar)
                        ENDIF
                    ELSE
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistVirtOcc(Bin)=SinglesHistVirtOcc(Bin)+abs(RealAttemptCreatePar)
                        ELSE
                            SinglesHistVirtVirt(Bin)=SinglesHistVirtVirt(Bin)+abs(RealAttemptCreatePar)
                        ENDIF
                    ENDIF
                ENDIF
            ELSEIF(IC.eq.2) THEN
                DoublesAttemptHist(Bin)=DoublesAttemptHist(Bin)+(Tau/Prob)
                IF(RealAttemptCreatePar.ne.0) THEN
                    DoublesHist(Bin)=DoublesHist(Bin)+abs(RealAttemptCreatePar)
                ENDIF
            ENDIF

            IF(tHPHF) THEN
                rh = hphf_diag_helement (nJ, iLutnJ)
            ELSE
                rh = get_helement (nJ, nJ, 0)
            ENDIF
            Bin=INT((real(rh,dp)-Hii)/BinRange)+1
            IF(Bin.gt.iNoBins) THEN
                CALL Stop_All("AttemptCreatePar","Histogramming energies higher than the arrays can cope with. "&
                & //"Increase iNoBins or BinRange")
            ENDIF
            IF(RealAttemptCreatePar.ne.0) THEN
                SpawnHist(Bin)=SpawnHist(Bin)+abs(RealAttemptCreatePar)
!                WRITE(iout,*) "Get Here!", abs(RealAttemptCreatePar),Bin
            ENDIF
            AttemptHist(Bin)=AttemptHist(Bin)+(Tau/Prob)
        ENDIF

    END FUNCTION AttemptCreatePar

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
        NoDied = NoDied + sum(min(iDie, abs(RealwSign)))

        ! Count any antiparticles
        iter_data%nborn = iter_data%nborn + max(iDie - abs(RealwSign), 0.0_dp)
        NoBorn = NoBorn + sum(max(iDie - abs(RealwSign), 0.0_dp))

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
                    CurrentH(2:1+lenof_sign,DetPosition) = wAvSign
                    CurrentH(2+lenof_sign:1+2*lenof_sign,DetPosition) = IterRDMStartCurr
                endif
            else
                call encode_bit_rep(CurrentDets(:,VecSlot),iLutCurr,CopySign,extract_flags(iLutCurr))
                if(.not.tRegenDiagHEls) CurrentH(1,VecSlot) = Kii
                if (tFillingStochRDMonFly) then
                    CurrentH(2:1+lenof_sign,VecSlot) = wAvSign
                    CurrentH(2+lenof_sign:1+2*lenof_sign,VecSlot) = IterRDMStartCurr
                endif
                VecSlot=VecSlot+1
            endif
        else
            !All walkers died
            if(tFillingStochRDMonFly) then
                ! If we're stochastically filling the RDMs, we want to keep determinants even if 
                ! their walkers have all died.
                ! This is because we're using the sign of each determinant before death.
                ! But if the walker is removed from the list altogether, this would be lost too.
                if(tHashWalkerList) then
                    call encode_sign(CurrentDets(:,DetPosition),CopySign)
                    CurrentH(2:1+lenof_sign,DetPosition) = wAvSign
                    CurrentH(2+lenof_sign:1+2*lenof_sign,DetPosition) = IterRDMStartCurr
                else
                    call encode_bit_rep(CurrentDets(:,VecSlot),iLutCurr,CopySign,extract_flags(iLutCurr))
                    if (.not.tRegenDiagHEls) CurrentH(1,VecSlot) = Kii
                    CurrentH(2:1+lenof_sign,VecSlot) = wAvSign
                    CurrentH(2+lenof_sign:1+2*lenof_sign,VecSlot) = IterRDMStartCurr
                    VecSlot=VecSlot+1
                endif
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
            if(tHashWalkerList.and.(.not.tFillingStochRDMonFly)) then
                !Remove the determinant from the indexing list
                call RemoveDetHashIndex(DetCurr,DetPosition)
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

    function attempt_die_normal (DetCurr, Kii, realwSign, WalkExcitLevel) result(ndie)
        
        ! Should we kill the particle at determinant DetCurr. 
        ! The function allows multiple births (if +ve shift), or deaths from
        ! the same particle. The returned number is the number of deaths if
        ! positive, and the
        !
        ! In:  DetCurr - The determinant to consider
        !      Kii     - The diagonal matrix element of DetCurr (-Ecore)
        !      wSign   - The sign of the determinant being considered. If
        !                |wSign| > 1, attempt to die multiple particles at
        !                once (multiply probability of death by |wSign|)
        ! Ret: ndie    - The number of deaths (if +ve), or births (If -ve).

        integer, intent(in) :: DetCurr(nel)
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        real(dp), intent(in) :: Kii
        real(dp), dimension(lenof_sign) :: ndie
        integer, intent(in) :: WalkExcitLevel

        real(dp) :: r, rat, probsign
        real(dp), dimension(inum_runs) :: fac
        integer :: i, iUnused

        do i=1, inum_runs
            fac(i)=tau*(Kii-DiagSft(i))

            ! And for tau searching purposes
            call log_death_magnitude (Kii - DiagSft(i))
        enddo

        if(fac(1).gt.1.0_dp) then
            if(fac(1).gt.2.0_dp) then
                call stop_all("attempt_die_normal","Death probability > 2: Algorithm unstable. Reduce timestep.")
            else
                write(iout,"(A,F20.10)") "** WARNING ** Death probability > 1: Creating Antiparticles. "&
                    & //"Timestep errors possible: ",fac
            endif
        endif


        if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel .le. RealCoeffExcitThresh)) &
            .or. tAllRealCoeff ) then
            do i=1, lenof_sign
                if(inum_runs.eq.2) then
                    ndie(i)=fac(i)*abs(realwSign(i))
                else
                    ndie(i)=fac(1)*abs(realwSign(i))
                endif
            enddo
        else
            do i=1,lenof_sign
                
                ! Subtract the current value of the shift, and multiply by tau.
                ! If there are multiple particles, scale the probability.
                if(inum_runs.eq.2) then
                    rat = fac(i) * abs(realwSign(i))
                else
                    rat = fac(1) * abs(realwSign(i))
                endif

                ndie(i) = real(int(rat), dp)
                rat = rat - ndie(i)

                ! Choose to die or not stochastically
                r = genrand_real2_dSFMT() 
                if (abs(rat) > r) ndie(i) = ndie(i) + real(nint(sign(1.0_dp, rat)), dp)
            enddo
        endif

        ! Avoid compiler warnings
        iUnused = DetCurr(1)

    end function

    

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

    subroutine new_child_stats_hist_hamil (iter_data, iLutI, nJ, iLutJ, ic, &
                                           walkExLevel, child, parent_flags, &
                                           part_type)
        ! Based on old AddHistHamilEl. Histograms the hamiltonian matrix, and 
        ! then calls the normal statistics routine.

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        real(dp), dimension(lenof_sign) , intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'new_child_stats_hist_hamil'
        integer :: partInd, partIndChild, childExLevel
        logical :: tSuccess

        if (walkExLevel == nel) then
            call BinSearchParts2 (iLutI, FCIDetIndex(walkExLevel), Det, &
                                  PartInd, tSuccess)
        else
            call BinSearchParts2 (iLutI, FCIDetIndex(walkExLevel), &
                                  FciDetIndex(walkExLevel+1)-1, partInd, &
                                  tSuccess)
        endif

        if (.not. tSuccess) &
            call stop_all (this_routine, 'Cannot find determinant nI in list')

        childExLevel = FindBitExcitLevel (iLutHF, iLutJ, nel)
        if (childExLevel == nel) then
            call BinSearchParts2 (iLutJ, FCIDetIndex(childExLevel), Det, &
                                  partIndChild, tSuccess)
        elseif (childExLevel == 0) then
            partIndChild = 1
            tSuccess = .true.
        else
            call BinSearchParts2 (iLutJ, FCIDetIndex(childExLevel), &
                                  FciDetIndex(childExLevel+1)-1, &
                                  partIndChild, tSuccess)
        endif

        histHamil (partIndChild, partInd) = &
                histHamil (partIndChild, partInd) + (1.0_dp * child(1))
        histHamil (partInd, partIndChild) = &
                histHamil (partInd, partIndChild) + (1.0_dp * child(1))
        avHistHamil (partIndChild, partInd) = &
                avHistHamil (partIndChild, partInd) + (1.0_dp * child(1))
        avHistHamil (partInd, partIndChild) = &
                avHistHamil (partInd, partIndChild) + (1.0_dp * child(1))

        ! Call the normal stats routine
        call new_child_stats_normal (iter_data, iLutI, nJ, iLutJ, ic, &
                                     walkExLevel, child, parent_flags, &
                                     part_type)

    end subroutine


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
        do part_type = 1, lenof_sign
            if (AllNoAtHF(part_type) < 0.0_dp) then
                root_print 'No. at HF < 0 - flipping sign of entire ensemble &
                           &of particles in simulation: ', part_type
                root_print AllNoAtHF(part_type)

                ! And do the flipping
                call FlipSign(part_type)
                AllNoatHF(part_type) = -AllNoatHF(part_type)
                NoatHF(part_type) = -NoatHF(part_type)

                if (tFillingStochRDMonFly) then
                    if (lenof_sign /= 1) &
                        call stop_all(t_r, 'Not yet implmented')
                    ! Want to flip all the averaged signs.
                    AvNoatHF = -AVNoatHF
                    InstNoatHF(part_type) = -InstNoatHF(part_type)
                end if
            end if
        end do
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
                IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(initiatorstats_unit)
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
                    elseif(tMomInv) then
                        if(.not.IsMomSelfInv(ProjEDet,iLutRef)) then
                            !Complications. We are now effectively projecting onto a LC of two dets.
                            !Ensure this is done correctly.
                            if(.not.Allocated(RefDetFlip)) then
                                allocate(RefDetFlip(NEl))
                                allocate(iLutRefFlip(0:NIfTot))
                                RefDetFlip = 0
                                iLutRefFlip = 0
                            endif
                            call CalcMomAllowedBitDet(ProjEDet,RefDetFlip,iLutRef,iLutRefFlip,.true.,.true.,tSwapped)
                            if(tSwapped) then
                                !The iLutRef should already be the correct one, since it was obtained by the normal calculation!
                                call stop_all("population_check","Error in changing reference det to momentum-coupled function")
                            endif
                            write(iout,"(A)") "Now projecting onto a momentum-coupled function as a linear combo of two dets..."
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
                    elseif(tMomInv) then
                        h_tmp = MI_diag_helement(ProjEDet,iLutRef)
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
                        elseif(tMomInv) then
                            h_tmp = MI_diag_helement(det,CurrentDets(:,i))
                        else
                            h_tmp = get_helement (det, det, 0)
                        endif
                        CurrentH(1,i) = real(h_tmp, dp) - Hii
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
                    elseif(tMomInv) then
                        h_tmp = MI_diag_helement(ProjEDet,iLutRef)
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
        
        if(tSplitProjEHist) then
            if(tSplitProjEHistG) then
                AllENumCycHistG=0.0_dp
                call MPISum(ENumCycHistG,AllENumCycHistG)
            endif
            if(tSplitProjEHistK3) then
                AllENumCycHistK3=0.0_dp
                call MPISum(ENumCycHistK3,AllENumCycHistK3)
            endif
        endif

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

        if (tSearchTau) &
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
                        if(tSpawn_Only_Init.and.tSpawn_Only_Init_Grow) then
                            !Remove the restriction that only initiators can spawn.
                            write(iout,*) "All determinants now with the ability to spawn new walkers."
                            tSpawn_Only_Init=.false.
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
                    if ((AllTotParts(run) > tot_walkers) .or. &
                         (abs(AllNoatHF(run)) > MaxNoatHF)) then
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
                        if(tSpawn_Only_Init.and.tSpawn_Only_Init_Grow) then
                            !Remove the restriction that only initiators can spawn.
                            write(iout,*) "All determinants now with the ability to spawn new walkers."
                            tSpawn_Only_Init=.false.
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
                 if (.not. (proje_linear_comb .and. nproje_sum > 1)) then
                     all_sum_proje_denominator(run) = ARR_RE_OR_CPLX(AllSumNoatHF,run)
                     all_cyc_proje_denominator(run) = AllHFCyc(run)
                 endif

                 ! Calculate the projected energy.
                 if((lenof_sign.eq.2).and.(inum_runs.eq.1)) then
                     if (any(AllSumNoatHF /= 0.0) .or. &
                         (proje_linear_comb .and. nproje_sum > 1)) then
                         ProjectionE = (AllSumENum) / (all_sum_proje_denominator) 
                         proje_iter = (AllENumCyc) / (all_cyc_proje_denominator) 
                        AbsProjE = (AllENumCycAbs) / (all_cyc_proje_denominator)
                    endif
                 else
                     if ((AllSumNoatHF(run) /= 0.0) .or. &
                         (proje_linear_comb .and. nproje_sum > 1)) then
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
        if(tSpawn_Only_Init_Grow) then
            !if tSpawn_Only_Init_Grow is on, the the tSpawn_Only_Init variable can change,
            !and thus needs to be broadcast in case it does.
            call MPIBcast(tSpawn_Only_Init)
        endif
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


    subroutine rezero_iter_stats_each_iter (iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp) :: TempTotParts
        real(dp), dimension(lenof_sign) :: AllInstNoatHF
        real(dp), dimension(lenof_sign) :: Prev_AvNoatHF

        NoInitDets = 0
        NoNonInitDets = 0
        NoInitWalk = 0.0_dp
        NoNonInitWalk = 0.0_dp
        InitRemoved = 0

        NoAborted = 0.0_dp
        NoRemoved = 0.0_dp
        NoatHF = 0.0_dp     ! Number at HF and doubles for stats
        NoatDoubs = 0.0_dp

        iter_data%nborn = 0.0
        iter_data%ndied = 0.0
        iter_data%nannihil = 0.0
        iter_data%naborted = 0.0
        iter_data%nremoved = 0.0

        if(tFillingStochRDMonFly) then
            call MPISumAll(InstNoatHF, AllInstNoAtHF)
            if(all(InstNoatHF(1:lenof_sign) == 0) &
                .and. (.not. tSemiStochastic)) then
                !The HF determinant won't be in currentdets, so the CurrentH averages will have been wiped.
                !NB - there will be a small issue here if the HF determinant isn't in the core space
                IterRDM_HF(1) = Iter+PreviousCycles + 1 
                AvNoatHF(1) = 0.0_dp
                if(inum_runs.eq.2) then
                    IterRDM_HF(inum_runs) = Iter+PreviousCycles + 1 
                    AvNoatHF(inum_runs) = 0.0_dp
                endif
               ! WRITE(6,*) "zeroed AvNoAtHF", Prev_AvNoAtHF, AvNoAtHF, InstNoAtHF
            else
                Prev_AvNoatHF(1) = AvNoatHF(1)
                AvNoatHF(1) = ( (real((Iter+PreviousCycles - IterRDM_HF(1)),dp) * Prev_AvNoatHF(1)) &
                    + InstNoatHF(1) ) / real((Iter+PreviousCycles - IterRDM_HF(1)) + 1,dp)
                if(inum_runs.eq.2) then
                    Prev_AvNoatHF(inum_runs) = AvNoatHF(inum_runs)
                    AvNoatHF(inum_runs) = ( (real((Iter+PreviousCycles - IterRDM_HF(inum_runs)),dp) * Prev_AvNoatHF(inum_runs)) &
                        + InstNoatHF(inum_runs) ) / real((Iter+PreviousCycles - IterRDM_HF(inum_runs)) + 1,dp)
                endif
!                WRITE(6,*) "AvNoAtHF", Prev_AvNoAtHF, AvNoAtHF, InstNoAtHF
            endif
        endif
        HFInd = 0            
        
        trial_ind = 0
        con_ind = 0
        min_trial_ind = 1
        min_conn_ind = 1

    end subroutine


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

        ! Reset the linear combination coefficients
        ! TODO: Need to rethink how/when this is done. This is just for tests
        if (proje_spatial .and. proje_linear_comb) &
             call update_linear_comb_coeffs()


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
                CurrentH(2:1+lenof_sign, i) = -CurrentH(2:1+lenof_sign, i)
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
            IF(tTruncInitiator.or.tDelayTruncInit) THEN
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
            if (tTruncInitiator .or. tDelayTruncInit) then
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
                                   &i13,g13.5,11g17.9,i13,g16.7)",advance = 'no') &
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
            if (tTruncInitiator .or. tDelayTruncInit) then
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
            if (tTruncInitiator .or. tDelayTruncInit) &
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

    SUBROUTINE SetupParameters()
        use SystemData, only : tUseBrillouin,iRanLuxLev,tSpn,tHPHFInts,tRotateOrbs,tROHF,tFindCINatOrbs,nOccBeta,nOccAlpha,tUHF
        use SystemData, only : tBrillouinsDefault,ECore,tNoSingExcits,tOddS_HPHF
        USE dSFMT_interface , only : dSFMT_init
        use CalcData, only: G_VMC_Seed, &
                            MemoryFacPart, MemoryFacAnnihil, TauFactor, &
                            StepsSftImag, tCheckHighestPop, tSpatialOnlyHash,tStartCAS, &
                            MaxWalkerBloom
        use Determinants , only : GetH0Element3,GetH0Element4
        use SymData , only : SymLabelList,SymLabelCounts,TwoCycleSymGens
        use DeterminantData , only : write_det
        use LoggingData , only : tTruncRODump,tCalcVariationalEnergy,tDiagAllSpaceEver
        use DetCalcData, only : NMRKS,tagNMRKS,FCIDets
        use SymExcit3, only : CountExcitations3 
        use constants, only: bits_n_int
        use HPHFRandExcitMod, only: ReturnAlphaOpenDet
        use sym_mod
        use HElem
        INTEGER :: ierr,i,j,HFDetTest(NEl),Seed,alpha,beta,symalpha,symbeta,endsymstate
        INTEGER :: HFConn,LargestOrb,nBits,HighEDet(NEl),orb
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        HElement_t :: TempHii
        TYPE(BasisFn) HFSym
        real(dp) :: TotDets,SymFactor,r,Gap,UpperTau
        CHARACTER(len=*), PARAMETER :: t_r='SetupParameters'
        CHARACTER(len=12) :: abstr
        character(len=24) :: filename, filename2
        LOGICAL :: tSuccess,tFoundOrbs(nBasis),FoundPair,tSwapped,tAlreadyOcc
        INTEGER :: HFLz,ChosenOrb,KPnt(3), step,FindProjEBins,SymFinal
        integer(int64) :: ExcitLevPop,SymHF

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors
        WRITE(iout,*) 
        if(nProcessors.gt.1) then
            WRITE(iout,*) "Performing Parallel FCIQMC...."
        else
            write(iout,*) "Performing FCIQMC...."
        endif
        WRITE(iout,*) 
        
!Set timed routine names
        Walker_Time%timer_name='WalkerTime'
        Annihil_Time%timer_name='AnnihilTime'
        Sort_Time%timer_name='SortTime'
        Comms_Time%timer_name='CommsTime'
        ACF_Time%timer_name='ACFTime'
        AnnSpawned_time%timer_name='AnnSpawnedTime'
        AnnMain_time%timer_name='AnnMainTime'
        BinSearch_time%timer_name='BinSearchTime'
        Sync_Time%timer_name='SyncTime'
        SemiStoch_Comms_Time%timer_name='SemiStochCommsTime'
        SemiStoch_Multiply_Time%timer_name='SemiStochMultiplyTime'
        Trial_Search_Time%timer_name='TrialSearchTime'
        SemiStoch_Init_Time%timer_name='SemiStochInitTime'
        Trial_Init_Time%timer_name='TrialInitTime'

        IF(TDebug) THEN
!This will open a file called LOCALPOPS-"iprocindex" on unit number 11 on every node.
            abstr=''
            write(abstr,'(I2)') iProcIndex
            abstr='LOCALPOPS-'//adjustl(abstr)
            OPEN(11,FILE=abstr,STATUS='UNKNOWN')
        ENDIF

        IF(iProcIndex.eq.Root) THEN
            if (.not. tFCIMCStats2) then
                fcimcstats_unit = get_free_unit()
                if (tReadPops) then
                    ! Restart calculation.  Append to stats file (if it exists).
                    if(tMolpro .and. .not. tMolproMimic) then
                        filename = 'FCIQMCStats_' // adjustl(MolproID)
                        OPEN(fcimcstats_unit,file=filename,status='unknown',position='append')
                    else
                        OPEN(fcimcstats_unit,file='FCIMCStats',status='unknown',position='append')
                    endif
                else
                    call MoveFCIMCStatsFiles()          !This ensures that FCIMCStats files are not overwritten
                    if(tMolpro .and. .not. tMolproMimic) then
                        filename = 'FCIQMCStats_' // adjustl(MolproID)
                        OPEN(fcimcstats_unit,file=filename,status='unknown')
                    else
                        OPEN(fcimcstats_unit,file='FCIMCStats',status='unknown')
                    endif
                end if
            end if
            if(inum_runs.eq.2) then
                fcimcstats_unit2 = get_free_unit()
                if (tReadPops) then
                    ! Restart calculation.  Append to stats file (if it exists).
                    if(tMolpro) then
                        filename2 = 'FCIQMCStats2_' // adjustl(MolproID)
                        OPEN(fcimcstats_unit2,file=filename2,status='unknown',position='append')
                    else
                        OPEN(fcimcstats_unit2,file='FCIMCStats2',status='unknown',position='append')
                    endif
                else
                    if(tMolpro) then
                        filename = 'FCIQMCStats_' // adjustl(MolproID)
                        OPEN(fcimcstats_unit2,file=filename2,status='unknown')
                    else
                        OPEN(fcimcstats_unit2,file='FCIMCStats2',status='unknown')
                    endif
                end if
            endif


            IF(tTruncInitiator.or.tDelayTruncInit) THEN
                initiatorstats_unit = get_free_unit()
                if (tReadPops) then
! Restart calculation.  Append to stats file (if it exists)
                    OPEN(initiatorstats_unit,file='INITIATORStats',status='unknown',form='formatted',position='append')
                else
                    OPEN(initiatorstats_unit,file='INITIATORStats',status='unknown',form='formatted')
                endif
            ENDIF
            IF(tLogComplexPops) THEN
                ComplexStats_unit = get_free_unit()
                OPEN(ComplexStats_unit,file='COMPLEXStats',status='unknown')
            ENDIF
        ENDIF

!Store information specifically for the HF determinant
        ALLOCATE(HFDet(NEl),stat=ierr)
        CALL LogMemAlloc('HFDet',NEl,4,t_r,HFDetTag)
        ALLOCATE(iLutHF(0:NIfTot),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(t_r,"Cannot allocate memory for iLutHF")

        do i=1,NEl
            HFDet(i)=FDet(i)
        enddo
        CALL EncodeBitDet(HFDet,iLutHF)

        if(tHPHF.and.TestClosedShellDet(iLutHF).and.tOddS_HPHF.and.TwoCycleSymGens) then
            !This is not a compatible reference function.
            !Create single excitation of the correct symmetry
            !Use this as the reference.
            write(6,"(A)") "Converging to ODD S eigenstate"

            SymFinal = int((G1(HFDet(nel))%Sym%S),sizeof_int)+1

            tAlreadyOcc=.false.
            do i=SymLabelCounts(1,SymFinal),SymLabelCounts(1,SymFinal)+SymLabelCounts(2,SymFinal)-1
                if(G1(HFDet(nel))%Ms.eq.-1) then
                    !Choose beta ones
                    orb=(2*SymLabelList(i))-1
                else
                    orb=(2*SymLabelList(i))
                endif
                tAlreadyOcc=.false.
                do j=1,nel
                    if(orb.eq.HFDet(j)) then
                        tAlreadyOcc=.true.
                        exit
                    endif
                enddo
                if(.not.tAlreadyOcc) then
                    HFDet(nel) = orb
                    call sort(HFDet)
                    exit
                endif
            enddo
            if(tAlreadyOcc)     &
                call stop_all(t_r,"Cannot automatically detect open-shell determinant for reference to use with odd S")
            call EncodeBitDet(HFDet,iLutHF)
            if(TestClosedShellDet(iLutHF))  &
                call stop_all(t_r,"Fail to create open-shell determinant for reference to use with odd S")
            write(6,*) "Reference determinant changed to the open-shell:"
            call write_det(iout,HFDet,.true.)
        endif

        !iLutRef is the reference determinant for the projected energy.
        !Initially, it is chosen to be the same as the inputted reference determinant
        ALLOCATE(iLutRef(0:NIfTot),stat=ierr)
        ALLOCATE(ProjEDet(NEl),stat=ierr)

        IF(ierr.ne.0) CALL Stop_All(t_r,"Cannot allocate memory for iLutRef")
        
        ! The reference / projected energy determinants are the same as the
        ! HF determinant.
        ! TODO: Make these pointers rather than copies?
        iLutRef = iLutHF
        ProjEDet = HFDet

        ALLOCATE(iLutHF_True(0:NIfTot),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(t_r,"Cannot allocate memory for iLutHF_True")
        ALLOCATE(HFDet_True(NEl),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(t_r,"Cannot allocate memory for HFDet_True")

        if(tRef_Not_HF) then
            do i = 1, NEl
                HFDet_True(i) = BRR(i)
            enddo
            call sort(HFDet_True(1:NEl))
            CALL EncodeBitDet(HFDet_True,iLutHF_True)
        else
            iLutHF_True = iLutHF
            HFDet_True = HFDet
        endif

        if(tHPHF) then
            if(.not.TestClosedShellDet(iLutRef)) then
                !We test here whether the reference determinant actually corresponds to an open shell HPHF.
                !If so, we need to ensure that we are specifying the correct determinant of the pair, and also
                !indicate that the projected energy needs to be calculated as a projection onto both of these determinants.
                ALLOCATE(RefDetFlip(NEl))
                ALLOCATE(iLutRefFlip(0:NIfTot))
                !We need to ensure that the correct pair of the HPHF det is used to project onto/start from.
                call ReturnAlphaOpenDet(ProjEDet,RefDetFlip,iLutRef,iLutRefFlip,.true.,.true.,tSwapped)
                if(tSwapped) then
                    write(iout,*) "HPHF used, and open shell determinant spin-flipped for consistency."
                endif
                write(iout,*) "Two *different* determinants contained in initial HPHF"
                write(iout,*) "Projected energy will be calculated as projection onto both of these"
                tSpinCoupProjE=.true.
            else
                tSpinCoupProjE=.false.
            endif
            if(tMomInv) then
                call stop_all(t_r,"Cannot currently have MomInv and HPHF functions. If this is important, bug ghb")
            endif
        elseif(tMomInv) then
            if(tAntisym_MI) then
                write(iout,*) "Using hilbert space of antisymmetric momentum-coupled determinants..."
            else
                write(iout,*) "Using hilbert space of symmetric momentum-coupled determinants..."
            endif
            if(.not.IsMomSelfInv(ProjEDet,iLutRef)) then
                !We test here whether the reference determinant actually corresponds to a momentum-coupled function.
                !If so, we need to ensure that we are specifying the correct determinant of the pair, and also
                !indicate that the projected energy needs to be calculated as a projection onto both of these determinants.
                ALLOCATE(RefDetFlip(NEl))
                ALLOCATE(iLutRefFlip(0:NIfTot))
                !We need to ensure that the correct pair of the reference det is used to project onto/start from.
                call CalcMomAllowedBitDet(ProjEDet,RefDetFlip,iLutRef,iLutRefFlip,.true.,.true.,tSwapped)
                if(tSwapped) then
                    write(iout,*) "Momentum-coupled function used, and reference determinant momentum-flipped for consistency."
                endif
                write(iout,*) "Two *different* determinants contained in initial reference function"
                write(iout,*) "Projected energy will be calculated as projection onto both of these"
                tSpinCoupProjE=.true.
            endif
            tSpinCoupProjE=.false.
        else
            tSpinCoupProjE=.false.
        endif

!Init hash shifting data
        hash_iter=0

        IF(tFixLz) THEN
            CALL GetLz(HFDet,NEl,HFLz)
            WRITE(iout,"(A,I5)") "Ml value of reference determinant is: ",HFLz
            IF(HFLz.ne.LzTot) THEN
                CALL Stop_All("SetupParameters","Chosen reference determinant does not have the " &
                & //"same Lz value as indicated in the input.")
            ENDIF
        ENDIF

!Do a whole lot of tests to see if we can use Brillouins theorem or not.
        IF(tBrillouinsDefault) CALL CheckforBrillouins() 
        
!test the encoding of the HFdet to bit representation.
        ! Test that the bit operations are working correctly...
        ! TODO: Move this to using the extract_bit_det routines to test those
        !       too...
        call decode_bit_det (HFDetTest, iLutHF)
        do i=1,NEl
            IF(HFDetTest(i).ne.HFDet(i)) THEN
                WRITE(iout,*) "HFDet: ",HFDet(:)
                WRITE(iout,*) "HFDetTest: ",HFDetTest(:)
                CALL Stop_All(t_r,"HF Determinant incorrectly decoded.")
            ENDIF
        enddo
        CALL LargestBitSet(iLutHF,NIfD,LargestOrb)
        IF(LargestOrb .ne. iand(HFDet(NEl), csf_orbital_mask)) THEN
            write(6,*) 'ilut HF', ilutHF
            write(6,*) 'largest orb', largestorb
            write(6,*) 'HFDet', iand(hfdet, csf_orbital_mask)
            CALL Stop_All(t_r,"LargestBitSet FAIL")
        ENDIF
        nBits = CountBits(iLutHF, NIfD, NEl)
        IF(nBits.ne.NEl) THEN
            CALL Stop_All(t_r,"CountBits FAIL")
        ENDIF

        ALLOCATE(HighestPopDet(0:NIfTot),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(t_r,"Cannot allocate memory for HighestPopDet")
        HighestPopDet(:)=0

!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess=.true.
        tFoundOrbs(:)=.false.

        IF((.not.tHub).and.(.not.tUEG).and.TwoCycleSymGens) THEN
            do i=1,nSymLabels
!                WRITE(iout,*) "NSymLabels: ",NSymLabels,i-1
                EndSymState=SymLabelCounts(1,i)+SymLabelCounts(2,i)-1
!                WRITE(iout,*) "Number of states: ",SymLabelCounts(2,i)
                do j=SymLabelCounts(1,i),EndSymState

                    Beta=(2*SymLabelList(j))-1
                    Alpha=(2*SymLabelList(j))
                    SymAlpha=INT((G1(Alpha)%Sym%S),sizeof_int)
                    SymBeta=INT((G1(Beta)%Sym%S),sizeof_int)
!                    WRITE(iout,*) "***",Alpha,Beta

                    IF(.not.tFoundOrbs(Beta)) THEN
                        tFoundOrbs(Beta)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF
                    IF(.not.tFoundOrbs(Alpha)) THEN
                        tFoundOrbs(Alpha)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF

                    IF(G1(Beta)%Ms.ne.-1) THEN
                        tSuccess=.false.
                    ELSEIF(G1(Alpha)%Ms.ne.1) THEN
                        tSuccess=.false.
                    ELSEIF((SymAlpha.ne.(i-1)).or.(SymBeta.ne.(i-1))) THEN
                        tSuccess=.false.
                    ENDIF
                enddo
            enddo
            do i=1,nBasis
                IF(.not.tFoundOrbs(i)) THEN
                    WRITE(iout,*) "Orbital: ",i, " not found."
                    CALL Stop_All("SetupParameters","Orbital not found")
                ENDIF
            enddo
        ENDIF
        IF(.not.tSuccess) THEN
            WRITE(iout,*) "************************************************"
            WRITE(iout,*) "**                 WARNING!!!                 **"
            WRITE(iout,*) "************************************************"
            WRITE(iout,*) "Symmetry information of orbitals not the same in alpha and beta pairs."
            WRITE(iout,*) "Symmetry now set up in terms of spin orbitals"
            WRITE(iout,*) "I strongly suggest you check that the reference energy is correct."
        ENDIF
        ! From now on, the orbitals are contained in symlabellist2 and 
        ! symlabelcounts2 rather than the original arrays. These are stored 
        ! using spin orbitals. Assume that if we want to use the non-uniform 
        ! random excitation generator, we also want to use the NoSpinSym full
        ! excitation generators if they are needed. 

        CALL GetSym(iand(HFDet, csf_orbital_mask),NEl,G1,NBasisMax,HFSym)
        
        Sym_Psi=INT(HFSym%Sym%S,sizeof_int)  !Store the symmetry of the wavefunction for later
        WRITE(iout,"(A,I10)") "Symmetry of reference determinant is: ",INT(HFSym%Sym%S,sizeof_int)
        SymHF=0
        do i=1,NEl
            SymHF=IEOR(SymHF,G1(iand(HFDet(i), csf_orbital_mask))%Sym%S)
        enddo
        WRITE(iout,"(A,I10)") "Symmetry of reference determinant from spin orbital symmetry info is: ",SymHF
        if(SymHF.ne.HFSym%Sym%S) then
            call stop_all(t_r,"Inconsistency in the symmetry arrays.")
        endif
        IF(tKPntSym) THEN
            CALL DecomposeAbelianSym(HFSym%Sym%S,KPnt)
            WRITE(iout,"(A,3I5)") "Crystal momentum of reference determinant is: ",KPnt(1),KPnt(2),KPnt(3)
        ENDIF

!If using a CAS space truncation, write out this CAS space
        IF(tTruncCAS) THEN
            IF(tTruncInitiator) THEN
                WRITE(iout,'(A)') " *********** INITIATOR METHOD IN USE ***********"
                WRITE(iout,'(A)') " Fixed initiator space defined using the CAS method."
            ELSE
                WRITE(iout,*) "Truncated CAS space detected. Writing out CAS space..."
            ENDIF
            WRITE(iout,'(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (", &
                (OccCASOrbs+VirtCASOrbs)/2,",",OccCASOrbs,")"
            DO I=NEl-OccCASorbs+1,NEl
                WRITE(iout,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(iout,G1(BRR(I))%SYM,.FALSE.)
                WRITE(iout,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(iout,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
            WRITE(iout,'(A)') " ================================================================================================="
            DO I=NEl+1,NEl+VirtCASOrbs
                WRITE(iout,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(iout,G1(BRR(I))%SYM,.FALSE.)
                WRITE(iout,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(iout,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
        ELSEIF(tTruncInitiator) THEN
            WRITE(iout,'(A)') "*********** INITIATOR METHOD IN USE ***********"
            WRITE(iout,'(A)') "Starting with only the reference determinant in the fixed initiator space."
        ENDIF

        if(tSpawn_Only_Init.and.tSpawn_Only_Init_Grow) then
            if(inum_runs.eq.2) call stop_all("SetupParameters", "Double run not set up to use with tSpawn_Only_Init_Grow.& 
                                & May need a separate logical for each population")
            write(iout,"(A)") "Only allowing initiator determinants to spawn during growth phase, but allowing all to "&
                & //"spawn after fixed shift."
        elseif(tSpawn_Only_Init) then
            if(inum_runs.eq.2) call stop_all("SetupParameters", "Double run not set up to use with tSpawn_Only_Init.& 
                                & May need a separate logical for each population")
            write(iout,"(A)") "Only allowing initiator determinants to spawn"
        endif
 
        ! Setup excitation generator for the HF determinant. If we are using 
        ! assumed sized excitgens, this will also be assumed size.
        IF(tUEG.or.tHub.or.tNoSingExcits) THEN
            exflag=2
        ELSE
            exflag=3
        ENDIF
        IF(.not.tKPntSym) THEN
!Count all possible excitations - put into HFConn
!TODO: Get CountExcitations3 working with tKPntSym
            CALL CountExcitations3(iand(HFDet, csf_orbital_mask),exflag,nSingles,nDoubles)
        ELSE
            ! Use Alex's old excitation generators to enumerate all excitations.
            call enumerate_sing_doub_kpnt(exflag,nSingles,nDoubles)
        ENDIF
        HFConn=nSingles+nDoubles

        ! Initialise random number seed - since the seeds need to be different
        ! on different processors, subract processor rank from random number
        if(.not.tRestart) then
            Seed=abs(G_VMC_Seed-iProcIndex)
            WRITE(iout,"(A,I12)") "Value for seed is: ",Seed
            !Initialise...
            CALL dSFMT_init(Seed)
            if(tMolpro) then
                if((NMCyc.eq.-1).and.(.not.tTimeExit)) then
                    !No iteration number, or TIME option has been specified.
                    call warning_neci(t_r,          &
                    "No iteration number specified. Only running for 100 iterations initially. Change with ITERATIONS option.")
                    NMCyc=100   !Only run for 100 iterations.
                elseif(tTimeExit.and.(NMCyc.eq.-1)) then
                    write(iout,"(A,F10.3,A)") "Running FCIQMC for ",MaxTimeExit/60.0_dp," minutes."
                elseif(tTimeExit.and.(NMCyc.ne.-1)) then
                    write(iout,"(A,F10.3,A,I15,A)") "Running FCIQMC for ",MaxTimeExit/60.0_dp," minutes OR ",NMCyc," iterations."
                elseif((.not.tTimeExit).and.(NMCyc.gt.0)) then
                    write(iout,"(A,I15,A)") "Running FCIQMC for ",NMCyc," iterations."
                else
                    call stop_all(t_r,"Iteration number/Time unknown for simulation - contact ghb")
                endif
            endif
        else
            !Reset the DiagSft to its original value
            DiagSft = InputDiagSft
        endif
        
        ! Option tRandomiseHashOrbs has now been removed.
        ! Its behaviour is now considered default
        ! --> Create a random mapping for the orbitals 
        ALLOCATE(RandomHash(nBasis),stat=ierr)
        IF(ierr.ne.0) THEN
            CALL Stop_All(t_r,"Error in allocating RandomHash")
        ENDIF
        RandomHash(:)=0
        if(tHashWalkerList .or. tSemiStochastic .or. (tTrialWavefunction .and. tTrialHash)) then
            !We want another independent randomizing array for the hash table, so we do not introduce
            !correlations between the two
            ALLOCATE(RandomHash2(nBasis),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All(t_r,"Error in allocating RandomHash2")
            ENDIF
            RandomHash2(:)=0
        endif
        IF(iProcIndex.eq.root) THEN
            do i=1,nBasis
                ! If we want to hash only by spatial orbitals, then the
                ! spin paired orbitals must be set equal
                if (tSpatialOnlyHash) then
                    if (.not. btest(i, 0)) then
                        RandomHash(i) = RandomHash(i - 1)
                        cycle
                    endif
                endif

                ! Ensure that we don't set two values to be equal accidentally
                FoundPair=.false.
                do while(.not.FoundPair)
                    r = genrand_real2_dSFMT()
                    ChosenOrb=INT(nBasis*r*1000)+1

                    ! Check all values which have already been set.
                    do j=1,nBasis
                        IF(RandomHash(j).eq.ChosenOrb) EXIT
                    enddo

                    ! If not already used, then we can move on
                    if (j == nBasis+1) FoundPair = .true.
                    RandomHash(i) = ChosenOrb
                enddo
            enddo
            if(tHashWalkerList .or. tSemiStochastic .or. (tTrialWavefunction .and. tTrialHash)) then
                !Do again for RandomHash2
                do i=1,nBasis
                    ! If we want to hash only by spatial orbitals, then the
                    ! spin paired orbitals must be set equal
                    if (tSpatialOnlyHash) then
                        if (.not. btest(i, 0)) then
                            RandomHash2(i) = RandomHash2(i - 1)
                            cycle
                        endif
                    endif
                    ! Ensure that we don't set two values to be equal accidentally
                    FoundPair=.false.
                    do while(.not.FoundPair)
                        r = genrand_real2_dSFMT()
                        ChosenOrb=INT(nBasis*r*1000)+1

                        ! Check all values which have already been set.
                        do j=1,nBasis
                            IF(RandomHash2(j).eq.ChosenOrb) EXIT
                        enddo

                        ! If not already used, then we can move on
                        if (j == nBasis+1) FoundPair = .true.
                        RandomHash2(i) = ChosenOrb
                    enddo
                enddo
            endif

!            WRITE(iout,*) "Random Orbital Indexing for hash:"
!            WRITE(iout,*) RandomHash2(:)
            if (tSpatialOnlyHash) then
                step = 2
            else
                step = 1
            endif
            do i=1,nBasis
                IF((RandomHash(i).eq.0).or.(RandomHash(i).gt.nBasis*1000)) THEN
                    CALL Stop_All(t_r,"Random Hash incorrectly calculated")
                ENDIF
                if(tHashWalkerList .or. tSemiStochastic .or. (tTrialWavefunction .and. tTrialHash)) then
                    IF((RandomHash2(i).eq.0).or.(RandomHash2(i).gt.nBasis*1000)) THEN
                        CALL Stop_All(t_r,"Random Hash 2 incorrectly calculated")
                    ENDIF
                endif
                do j = i+step, nBasis, step
                    IF(RandomHash(i).eq.RandomHash(j)) THEN
                        CALL Stop_All(t_r,"Random Hash incorrectly calculated")
                    ENDIF
                    if(tHashWalkerList .or. tSemiStochastic .or. (tTrialWavefunction .and. tTrialHash)) then
                        IF(RandomHash2(i).eq.RandomHash2(j)) THEN
                            CALL Stop_All(t_r,"Random Hash 2 incorrectly calculated")
                        ENDIF
                    endif
                enddo
            enddo
        ENDIF
        !Now broadcast to all processors
        CALL MPIBCast(RandomHash,nBasis)
        if(tHashWalkerList .or. tSemiStochastic .or. (tTrialWavefunction .and. tTrialHash)) call MPIBCast(RandomHash2,nBasis)

        IF(tHPHF) THEN
            !IF(tLatticeGens) CALL Stop_All("SetupParameters","Cannot use HPHF with model systems currently.")
            IF(tROHF.or.(LMS.ne.0)) CALL Stop_All("SetupParameters","Cannot use HPHF with high-spin systems.")
            tHPHFInts=.true.
        ENDIF

        if(tMomInv) then
            if(.not.tFixLz) then
                call stop_all("SetupParameters","Cannot use MI functions without Lz conservation")
            endif
            if(LzTot.ne.0) then
                call stop_all("SetupParameters","Cannot use MI functions if Lz is not zero")
            endif
        endif

!Calculate Hii
        IF(tHPHF) THEN
            TempHii = hphf_diag_helement (HFDet, iLutHF)
        elseif(tMomInv) then
            TempHii = MI_diag_helement(HFDet,iLutHF)
        ELSE
            TempHii = get_helement (HFDet, HFDet, 0)
        ENDIF
        Hii=REAL(TempHii,dp)
        WRITE(iout,"(A,F20.10)") "Reference Energy set to: ",Hii
        if(tUEG) then
            !We require calculation of the sum of fock eigenvalues,
            !without knowing them - calculate from the full 1e matrix elements
            !of full hamiltonian removing two electron terms.
            TempHii=GetH0Element4(iand(HFDet, csf_orbital_mask),iand(HFDet, csf_orbital_mask))
        else
            !Know fock eigenvalues, so just use these.
            TempHii=GetH0Element3(iand(HFDet, csf_orbital_mask))
        endif
        Fii=REAL(TempHii,dp)

!Find the highest energy determinant...
        IF(LMS.eq.0) THEN
            do i=1,NEl
                HighEDet(i)=Brr(nBasis-(i-1))
            enddo
            IF(tHPHF) THEN
                call EncodeBitDet (HighEDet, iLutTemp)
                TempHii = hphf_diag_helement (HighEDet, iLutTemp)
            elseif(tMomInv) then
                call EncodeBitDet (HighEDet, iLutTemp)
                TempHii = MI_diag_helement (HighEDet, iLutTemp)
            ELSE
                TempHii = get_helement (HighEDet, HighEDet, 0)
            ENDIF
            UpperTau = 1.0_dp/REAL(TempHii-Hii,dp)
!            WRITE(iout,"(A,G25.15)") "Highest energy determinant is (approximately): ",REAL(TempHii,dp)
!            WRITE(iout,"(A,F25.15)") "This means tau should be no more than about ",UpperTau
!            WRITE(iout,*) "Highest energy determinant is: ", HighEDet(:)
        else
            UpperTau=0.0_dp
        ENDIF

        IF(tHub) THEN
            IF(tReal) THEN
!We also know that in real-space hubbard calculations, there are only single excitations.
                exFlag=1
            ELSE
!We are doing a momentum space hubbard calculation - set exFlag to 2 since only doubles are connected for momentum conservation.
                exFlag=2
            ENDIF
        ENDIF


        IF(LMS.ne.0) THEN
            IF(tNoBrillouin.or.(tHub.and.tReal).or.tRotatedOrbs) THEN
                WRITE(iout,*) "No brillouin theorem assumed. Single excitations also used to calculate projected energy."
            ELSEIF(tUHF) THEN
                WRITE(iout,*) "High spin calculation - but single excitations will *NOT* be used to calculate energy as "&
                & //"this is an unrestricted calculation."
            ELSE
                CALL Stop_All("SetupParameters","High-spin, restricted calculation detected, but single excitations are "&
                & //"not being used to calculate the energy.  &
                  & Either use the UHF keyword, or turn off brillouins theorem using NOBRILLOUINS, ROHF or ROTATEDORBS.")
            ENDIF
!            tRotatedOrbs=.true.
!        ELSEIF(LMS.ne.0) THEN
!            CALL Stop_All(t_r,"Ms not equal to zero, but tSpn is false. Error here")
        ENDIF

!Initialise variables for calculation on each node
        iter=0          !This is set so that calls to CalcParentFlag in the initialisation are ok with the logging.
        iPopsTimers=1   !Number of timed popsfiles written out
        iBlockingIter=0
        IterTime=0.0
        ProjectionE(:)=0.0_dp
        AvSign=0.0_dp
        AvSignHFD=0.0_dp
        SumENum(:)=0.0_dp
        SumNoatHF(:)=0.0_dp
        NoatHF(:)=0.0_dp
        InstNoatHF(:)=0.0_dp
        AvNoatHF(:) = 0.0_dp
        Annihilated(:)=0.0_dp
        Acceptances(:)=0.0_dp
        PreviousCycles=0
        NoBorn=0.0_dp
        SpawnFromSing=0
        NoDied=0
        HFCyc=0.0_dp
        ENumCyc=0.0_dp
        VaryShiftCycles=0
        AvDiagSft(:)=0.0_dp
        SumDiagSft(:)=0.0_dp
        SumWalkersCyc(:)=0.0_dp
!        SumDiagSftAbort=0.0_dp
!        AvDiagSftAbort=0.0_dp
        NoAborted(:)=0.0_dp
        NoRemoved(:)=0.0_dp
        NoAddedInitiators=0
        NoInitDets=0
        NoNonInitDets=0
        NoInitWalk(:)=0.0_dp
        NoNonInitWalk(:)=0.0_dp
        NoExtraInitDoubs=0
        InitRemoved=0
        TotImagTime=0.0_dp
        DiagSftRe=0.0_dp
        DiagSftIm=0.0_dp
        sum_proje_denominator = 0
        cyc_proje_denominator = 0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum(:)=0.0_dp
        AllNoatHF(:)=0.0_dp
        AllNoatDoubs(:)=0.0_dp
        AllSumNoatHF(:)=0.0_dp
        AllGrowRate(:)=0.0_dp
        AllGrowRateAbort(:)=0
!        AllMeanExcitLevel=0.0_dp
        AllSumWalkersCyc(:)=0
        AllAvSign=0.0_dp
        AllAvSignHFD=0.0_dp
        AllNoBorn(:)=0
        AllSpawnFromSing(:)=0
        AllNoDied(:)=0
        AllAnnihilated(:)=0
        AllENumCyc(:)=0.0_dp
        AllHFCyc(:)=0.0_dp
!        AllDetsNorm=0.0_dp
        AllNoAborted=0
        AllNoRemoved=0
        AllNoAddedInitiators=0
        AllNoInitDets=0
        AllNoNonInitDets=0
        AllNoInitWalk=0.0_dp
        AllNoNonInitWalk=0.0_dp
        AllNoExtraInitDoubs=0
        AllInitRemoved=0

        ! Initialise the fciqmc counters
        iter_data_fciqmc%update_growth = 0.0_dp
        iter_data_fciqmc%update_iters = 0
 
        IF(tHistSpawn.or.(tCalcFCIMCPsi.and.tFCIMC).or.tHistHamil) THEN
            ALLOCATE(HistMinInd(NEl))
            ALLOCATE(HistMinInd2(NEl))
            maxdet=0
            do i=1,nel
                maxdet=maxdet+2**(nbasis-i)
            enddo

            IF(.not.allocated(FCIDets)) THEN
                CALL Stop_All(t_r,"A Full Diagonalization is required before histogramming can occur.")
            ENDIF

            IF(tHistHamil) THEN
                WRITE(iout,*) "Histogramming total Hamiltonian, with Dets=", Det
                ALLOCATE(HistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays") 
                HistHamil(:,:)=0.0_dp
                ALLOCATE(AvHistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                AvHistHamil(:,:)=0.0_dp
                IF(iProcIndex.eq.0) THEN
                    ALLOCATE(AllHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                    AllHistHamil(:,:)=0.0_dp
                    ALLOCATE(AllAvHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                    AllAvHistHamil(:,:)=0.0_dp
                ENDIF
            ELSE
                WRITE(iout,*) "Histogramming spawning wavevector, with Dets=", Det
                ALLOCATE(Histogram(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays ")
                ENDIF
                Histogram(:,:)=0.0_dp
                ALLOCATE(AllHistogram(1:lenof_sign,1:det),stat=ierr)
                ALLOCATE(BeforeNormHist(1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
            ENDIF
            IF(tHistSpawn) THEN
                ALLOCATE(InstHist(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
                InstHist(:,:)=0.0_dp
                ALLOCATE(AvAnnihil(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
                AvAnnihil(:,:)=0.0_dp
                ALLOCATE(InstAnnihil(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
                InstAnnihil(:,:)=0.0_dp
            ENDIF

            IF(iProcIndex.eq.0) THEN
                IF(tHistSpawn) THEN
                    Tot_Unique_Dets_Unit = get_free_unit()
                    OPEN(Tot_Unique_Dets_Unit,FILE='TOTUNIQUEDETS',STATUS='UNKNOWN')
                    if(tCalcVariationalEnergy.and.tDiagAllSpaceEver) then
                        write(Tot_Unique_Dets_Unit,"(A)") "# Iter  UniqueDetsEver  AvVarEnergy   InstVarEnergy   GroundE_Ever"
                    elseif(tCalcVariationalEnergy) then
                        write(Tot_Unique_Dets_Unit,"(A)") "# Iter  UniqueDetsEver  AvVarEnergy   InstVarEnergy"
                    elseif(tDiagAllSpaceEver) then
                        write(Tot_Unique_Dets_Unit,"(A)") "# Iter  UniqueDetsEver  GroundE_Ever"
                    else
                        write(Tot_Unique_Dets_Unit,"(A)") "# Iter  UniqueDetsEver"
                    endif
                    ALLOCATE(AllInstHist(1:lenof_sign,1:det),stat=ierr)
                    ALLOCATE(AllInstAnnihil(1:lenof_sign,1:det),stat=ierr)
                    ALLOCATE(AllAvAnnihil(1:lenof_sign,1:det),stat=ierr)
                ENDIF
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
            ENDIF
        ELSEIF(tHistEnergies) THEN
            WRITE(iout,*) "Histogramming the energies of the particles, with iNoBins=",iNoBins, " and BinRange=", BinRange
            WRITE(iout,*) "Histogramming spawning events from ",-OffDiagMax, " with BinRange = ", OffDiagBinRange
            iOffDiagNoBins=INT((2.0_dp*OffDiagMax)/OffDiagBinRange)+1
            WRITE(iout,*) "This gives ",iOffDiagNoBins," bins to histogram the off-diagonal matrix elements."
            ALLOCATE(HistogramEnergy(1:iNoBins))
            ALLOCATE(AttemptHist(1:iNoBins))
            ALLOCATE(SpawnHist(1:iNoBins))
            ALLOCATE(SinglesHist(1:iOffDiagNoBins))
            ALLOCATE(SinglesAttemptHist(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistOccOcc(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistOccVirt(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistVirtOcc(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistVirtVirt(1:iOffDiagNoBins))
            ALLOCATE(DoublesHist(1:iOffDiagNoBins))
            ALLOCATE(DoublesAttemptHist(1:iOffDiagNoBins))
            HistogramEnergy(:)=0.0_dp
            AttemptHist(:)=0.0_dp
            SpawnHist(:)=0.0_dp
            SinglesHist(:)=0.0_dp
            SinglesAttemptHist(:)=0.0_dp
            SinglesHistOccOcc(:)=0.0_dp
            SinglesHistOccVirt(:)=0.0_dp
            SinglesHistVirtOcc(:)=0.0_dp
            SinglesHistVirtVirt(:)=0.0_dp
            DoublesHist(:)=0.0_dp
            DoublesAttemptHist(:)=0.0_dp
            IF(iProcIndex.eq.Root) THEN
                ALLOCATE(AllHistogramEnergy(1:iNoBins))
                ALLOCATE(AllAttemptHist(1:iNoBins))
                ALLOCATE(AllSpawnHist(1:iNoBins))
                ALLOCATE(AllSinglesHist(1:iOffDiagNoBins))
                ALLOCATE(AllDoublesHist(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesAttemptHist(1:iOffDiagNoBins))
                ALLOCATE(AllDoublesAttemptHist(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistOccOcc(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistOccVirt(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistVirtOcc(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistVirtVirt(1:iOffDiagNoBins))
            ENDIF
        ENDIF

        if((iProcIndex.eq.Root).and.tDiagWalkerSubspace) then
            write(iout,'(A,I9,A)') "Diagonalising walker subspace every ",iDiagSubspaceIter," iterations"
            unitWalkerDiag = get_free_unit()
            open(unitWalkerDiag,file='WalkerSubspaceDiag',status='unknown')
            if(tTruncInitiator) then
                write(unitWalkerDiag,'(A)') "# Iter    NoInitDets   NoOccDets  InitiatorSubspaceEnergy   &
                        & FullSubspaceEnergy   ProjInitEnergy   ProjFullEnergy"
            else
                write(unitWalkerDiag,'(A)') "# Iter    NoOccDets    InitiatorSubspaceEnergy         &
                        & FullSubspaceEnergy    ProjFullEnergy"
            endif
        endif

        ! Initialise the spin distribution histogramming
        if (tHistSpinDist) then
            call init_hist_spin_dist ()
        endif

        if (tHistExcitToFrom) &
            call init_hist_excit_tofrom()

!Need to declare a new MPI type to deal with the long integers we use in the hashing, and when reading in from POPSFILEs
!        CALL MPI_Type_create_f90_integer(18,mpilongintegertype,error)
!        CALL MPI_Type_commit(mpilongintegertype,error)
        IF(tUseBrillouin) THEN
            WRITE(iout,"(A)") "Brillouin theorem in use for calculation of projected energy." 
        ENDIF
!        WRITE(iout,*) "Non-uniform excitation generators in use."
        CALL CalcApproxpDoubles()
        IF(TauFactor.ne.0.0_dp) THEN
            WRITE(iout,*) "TauFactor detected. Resetting Tau based on connectivity of: ",HFConn
            Tau=TauFactor/REAL(HFConn,dp)
            WRITE(iout,*) "Timestep set to: ",Tau
        ENDIF

        if(tSearchTau) &
            call init_tau_search()

        IF(StepsSftImag.ne.0.0_dp) THEN
            WRITE(iout,*) "StepsShiftImag detected. Resetting StepsShift."
            StepsSft=NINT(StepsSftImag/Tau)
            IF(StepsSft.eq.0) StepsSft=1
            WRITE(iout,*) "StepsShift set to: ",StepsSft
        ENDIF
        ! Once Alex's CCMC/FCIMC unification project is finished, we can
        ! remove this. The problem is that ValidSpawnedList is now setup in
        ! InitFCIMCCalcPar - this is for compatibility with POPSFILE v.3,
        ! where MaxSpawned is calculated from the number in the POPSFILE.
        !
        ! Set it pu here for CCMC.
        IF(tCCMC) then
            if(inum_runs.eq.2) call stop_all('SetupParameters',"CCMC not set up to work with double run")
            Call SetupValidSpawned(int(InitWalkers, int64))
        endif

        IF(TPopsFile) THEN
            IF(mod(iWritePopsEvery,StepsSft).ne.0) then
                CALL Warning_neci(t_r,"POPSFILE writeout should be a multiple of the update cycle length.")
            endif
        ENDIF

        if (TReadPops) then
            if (tStartSinglePart .and. .not. tReadPopsRestart) then
                call warning_neci(t_r, &
                               "ReadPOPS cannot work with StartSinglePart: ignoring StartSinglePart")
                tStartSinglePart = .false.
            end if
        endif

        IF(.not.TReadPops) THEN
            WRITE(iout,"(A,F20.10)") "Initial Diagonal Shift: ", DiagSft(1)
        ENDIF
!        WRITE(iout,*) "Damping parameter for Diag Shift set to: ", SftDamp
        if(tShiftonHFPop) then
            write(iout,*) "Shift will be varied in order to keep the population on the reference determinant fixed"
        endif
        WRITE(iout,*) "Connectivity of HF determinant is: ",HFConn
        IF(TStartSinglePart) THEN
            TSinglePartPhase(:)=.true.
        ELSE
            TSinglePartPhase(:)=.false.
        ENDIF
        
        IF(ICILevel.ne.0) THEN
!We are truncating the excitations at a certain value
            TTruncSpace=.true.
            WRITE(iout,'(A,I4)') "Truncating the S.D. space at determinants will an excitation level w.r.t. HF of: ",ICILevel
        ENDIF
        IF(tTruncCAS.or.tStartCAS) THEN
            ! We are truncating the FCI space by only allowing excitations 
            ! in a predetermined CAS space.
            ! The following line has already been written out if we are doing
            ! a CAS calculation.

!            WRITE(iout,'(A,I4,A,I5)') "Truncating the S.D. space as &
!                                   &determinants must be within a CAS of ", &
!                                   OccCASOrbs, " , ", VirtCASOrbs
            ! The SpinInvBRR array is required for the tTruncCAS option. Its 
            ! properties are explained more fully in the subroutine. 

            CALL CreateSpinInvBRR()

            ! CASmax is the max spin orbital number (when ordered 
            ! energetically) within the chosen active space.
            ! Spin orbitals with energies larger than this maximum value must 
            ! be unoccupied for the determinant to be in the active space.
            CASmax=NEl+VirtCASorbs

            ! CASmin is the max spin orbital number below the active space.  
            ! As well as the above criteria, spin orbitals with energies
            ! equal to, or below that of the CASmin orbital must be completely
            ! occupied for the determinant to be in the active space.
            CASmin=NEl-OccCASorbs

            IF(OccCASOrbs.gt.NEl) then
                CALL Stop_All("SetupParameters","Occupied orbitals in CAS space specified is greater than number of electrons")
            endif
            IF(VirtCASOrbs.gt.(nBasis-NEl)) then
                CALL Stop_All("SetupParameters","Virtuals in CAS space specified greater than number of unoccupied orbitals")
            endif

!Create the bit masks for the bit calculation of these properties.
            ALLOCATE(CASMask(0:NIfD))
            ALLOCATE(CoreMask(0:NIfD))
            CASMask(:)=0
            CoreMask(:)=0
            do i=1,nBasis
                IF(SpinInvBRR(i).gt.CASmax) THEN
                    !Orbital is in external space
                    CASMask((SpinInvBRR(i)-1)/bits_n_int) = ibset(CASMask((i-1)/bits_n_int),MOD((i-1),bits_n_int))

                ELSEIF(SpinInvBRR(i).le.CASmin) THEN
                    !Orbital is in core space
                    CoreMask((SpinInvBRR(i)-1)/bits_n_int) = ibset(CoreMask((i-1)/bits_n_int),MOD((i-1),bits_n_int))
                    CASMask((SpinInvBRR(i)-1)/bits_n_int) = ibset(CASMask((i-1)/bits_n_int),MOD((i-1),bits_n_int))
                ENDIF
            enddo

        ENDIF
        IF(tPartFreezeCore) THEN
            WRITE(iout,'(A,I4,A,I5)') 'Partially freezing the lowest ',NPartFrozen,' spin orbitals so that no more than ', &
                NHolesFrozen,' holes exist within this core.'
            CALL CreateSpinInvBRR()
        ENDIF
        IF(tPartFreezeVirt) THEN
            WRITE(iout,'(A,I4,A,I5)') 'Partially freezing the highest ',NVirtPartFrozen, &
                ' virtual spin orbitals so that no more than ',NElVirtFrozen,' electrons occupy these orbitals.'
            CALL CreateSpinInvBRR()
        ENDIF

        if (tTruncNOpen) then
            write(iout, '("Truncating determinant space at a maximum of ",i3," &
                    &unpaired electrons.")') trunc_nopen_max
        endif

!        SymFactor=(Choose(NEl,2)*Choose(nBasis-NEl,2))/(HFConn+0.0_dp)
!        TotDets=1.0_dp
!        do i=1,NEl
!            ExcitLevPop=int((Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor,int64)
!            if(ExcitLevPop.lt.0) cycle  !int overflow
!            WRITE(iout,"(A,I5,I20)") "Approximate excitation level population: ",i,ExcitLevPop
!            TotDets=TotDets+(Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor
!        enddo
!        if(TotDets.gt.0) then
!            WRITE(iout,"(A,I20)") "Approximate size of determinant space is: ",NINT(TotDets)
!        endif

    END SUBROUTINE SetupParameters

    subroutine check_start_rdm()
! This routine checks if we should start filling the RDMs - and does so if we should.        
        use nElRDMMod , only : DeAlloc_Alloc_SpawnedParts
        use LoggingData, only: tReadRDMs, tReadRDMAvPop
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

            if(tReadRDMs .and. tReadRDMAvPop) then
                !We need to read in the values of IterRDMStart and IterRDM_HF
                iunit_4=get_free_unit()
                OPEN(iunit_4,FILE='ITERRDMSTART',status='old')
                read(iunit_4, *) IterRDMStart, IterRDM_HF, AvNoAtHF

            endif

            !We have reached the iteration where we want to start filling the RDM.
            if(tExplicitAllRDM) then
                ! Explicitly calculating all connections - expensive...
                if(inum_runs.eq.2) call stop_all('check_start_rdm',"Cannot yet do replica RDM sampling with explicit RDMs. &
                    & e.g Hacky bit in Gen_Hist_ExcDjs to make it compile")
                
                tFillingExplicRDMonFly = .true.
                if(tHistSpawn) NHistEquilSteps = Iter
            else
                !extract_bit_rep_avsign => extract_bit_rep_avsign_norm
                !By default - we will do a stochastic calculation of the RDM.
                tFillingStochRDMonFly = .true.
                if(.not.tHF_Ref_Explicit) call DeAlloc_Alloc_SpawnedParts()
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

    SUBROUTINE CheckforBrillouins()
        use SystemData, only : tUseBrillouin,tNoBrillouin,tUHF
        use Determinants, only : tDefineDet
        INTEGER :: i,j
        LOGICAL :: tSpinPair
        
!Standard cases.
        IF((tHub.and.tReal).or.(tRotatedOrbs).or.((LMS.ne.0).and.(.not.tUHF))) THEN
!Open shell, restricted.            
            tNoBrillouin=.true.
        ELSE
!Closed shell restricted, or open shell unrestricted are o.k.            
            tNoBrillouin=.false.
            tUseBrillouin=.true.
        ENDIF

!Special case of complex orbitals.        
        IF(tFixLz.and.(.not.tNoBrillouin)) THEN
            WRITE(iout,*) "Turning Brillouins theorem off since we are using non-canonical complex orbitals"
            tNoBrillouin=.true.
        ENDIF

        ! Special case of defining a det with LMS=0, but which is open shell.
        ! No Brillouins if it's a restricted HF calc.
        tSpinPair = .false.
        IF(tDefineDet.and.(LMS.eq.0).and.(.not.tUHF)) THEN
            ! If we are defining our own reference determinant, we want to 
            ! find out if it is open shell or closed to know whether or not 
            ! brillouins theorem holds.
            !
            ! If LMS/=0, then it is easy and must be open shell, otherwise 
            ! we need to consider the occupied orbitals.
            do i=1,(NEl-1),2
                ! Assuming things will probably go alpha beta alpha beta, 
                ! run through each alpha and see if there's a corresponding 
                ! beta.
                tSpinPair=.false.
                IF(MOD(BRR(FDet(i)),2).ne.0) THEN
!Odd energy, alpha orbital.                    
                    IF(BRR(FDet(i+1)).ne.(BRR(FDet(i))+1)) THEN
                        ! Check the next orbital to see if it's the beta (will
                        ! be alpha+1 when ordered by energy). If not, check 
                        ! the other orbitals for the beta, as it's possible 
                        ! the orbitals are ordered weird (?).
                        do j=1,NEl
                            IF(BRR(FDet(j)).eq.(BRR(FDet(i))+1)) tSpinPair=.true. 
                        enddo
                    ELSE
                        tSpinPair=.true.
                    ENDIF
                ELSE
!Even energy, beta orbital. The corresponding alpha will be beta-1.                    
                    IF(BRR(FDet(i+1)).ne.(BRR(FDet(i))-1)) THEN
                        do j=1,NEl
                            IF(BRR(FDet(j)).eq.(BRR(FDet(i))-1)) tSpinPair=.true. 
                        enddo
                    ELSE
                        tSpinPair=.true.
                    ENDIF
                ENDIF
                IF(.not.tSpinPair) EXIT
            enddo
            IF(.not.tSpinPair) THEN
!Open shell LMS=0 determinant.
!If restricted HF orbitals are being used, brillouins theorem does not hold.
                tNoBrillouin=.true.
                tUseBrillouin=.false.
                WRITE(iout,'(A)') " Using an open shell reference determinant in a basis of restricted HF orbitals; " &
                & //"Brillouins theorem is being turned off. "
            ENDIF
        ENDIF

    ENDSUBROUTINE CheckforBrillouins

    LOGICAL FUNCTION TestifDETinCASBit(iLutnI)
        ! In:
        !    iLutNI: bit string representation of a determinant.
        ! Returns:
        !    true if the determinant is in the complete active space.
        INTEGER(KIND=n_int), INTENT(IN) :: iLutnI(0:NIfD)

        ! A determinant is in the CAS iff
        !  a) all orbitals in the core space are occupied;
        !  b) no orbitals in the external space are occupied;
        ! Thus ANDing the determinant with CASMask (containing set bits for the
        ! core and external orbitals) will give precisely the core orbitals
        ! if the determinant is in the CAS.
        TestifDETinCASBit = all(iand(iLutNI,CASMask) == CoreMask)

    END FUNCTION TestifDETinCASBit

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




    function CheckAllowedTruncSpawn (WalkExcitLevel, nJ, ilutnJ, IC) &
                                    result(bAllowed)

        ! Under any currently applied truncation schemes, is an excitation to
        ! this determinant allowed?
        !
        ! In:  WalkExcitLevel - Current excitation level relative to HF
        !      nJ             - Natural integer representation of det
        !                       (not Needed for HPHF/tTruncNOpen/MomInv)
        !      ilutnJ         - Bit representation of det
        !      IC             - Excitation level relative to parent
        ! Ret: bAllowed       - .true. if excitation is allowed

        integer, intent(in) :: nJ(nel), WalkExcitLevel, IC
        integer(n_int), intent(in) :: ilutnJ(0:NIfTot)
        logical :: bAllowed

        integer :: NoInFrozenCore, MinVirt, ExcitLevel, nopen, i
        ! For UEG
        integer :: k(3)

        bAllowed = .true.

        ! Truncate space by excitation level
        if (tTruncSpace) then
! If parent walker is one below excitation cutoff, could be
! disallowed if double. If higher, then all excits could
! be disallowed. If HPHF, excit could be single or double,
! and IC not returned --> Always test.
            if (tMomInv .or. tHPHF .or. WalkExcitLevel >= ICILevel .or. &
                (WalkExcitLevel == (ICILevel-1) .and. IC == 2)) then
                ExcitLevel = FindBitExcitLevel (iLutHF, ilutnJ, ICILevel)
                if (ExcitLevel > ICILevel) &
                    bAllowed = .false.
            endif
        endif

        ! Is the number of unpaired electrons too high?
        if (tTruncNOpen .and. bAllowed) then
            if (count_open_orbs(ilutnJ) > trunc_nopen_max) &
                bAllowed = .false.
        endif


        ! If the FCI space is restricted by a predetermined CAS space
        if (tTruncCAS .and. .not. tTruncInitiator .and. bAllowed) then
            if (.not. TestIfdetinCASBit(ilutnJ(0:NIfD))) &
                bAllowed = .false.
        endif


        ! Does the spawned determinant have more than the restricted number
        ! of holes in the partially frozen core?
        !
        ! --> Run through the e- in nJ, count the number in the partially
        !     frozen core (i.e. with energy, from BRR, less than the frozen
        !     core limit). If too few, then forbidden.
        if (tPartFreezeCore .and. bAllowed) then
            NoInFrozenCore = 0
            bAllowed = .false.
            do i = 1, nel
                if (SpinInvBRR(nJ(i)) <= NPartFrozen) &
                    NoInFrozenCore = NoInFrozenCore + 1
                if (NoInFrozenCore == (NPartFrozen - NHolesFrozen)) then
                    bAllowed = .true.
                    exit
                endif
            enddo
        endif


        ! Does the spawned determinant have more than the restricted number
        ! of e- in the partially frozen virtual orbitals?
        !
        ! --> Run through the e- in nJ, count the number in the partially
        !     frozen orbitals (i.e. with energy, from BRR, greater than
        !     minumum unfrozen virtual). If too many, then forbidden
        if (tPartFreezeVirt .and. bAllowed) then
            NoInFrozenCore = 0
            MinVirt = nBasis - NVirtPartFrozen
            ! BRR(i) = j: orbital i is the j-th lowest in energy
            do i = 1, nel
                if (SpinInvBRR(nJ(i)) > MinVirt) &
                    NoInFrozenCore = NoInFrozenCore + 1
                if (NoInFrozenCore > NElVirtFrozen) then
                    ! Too many e- in part-frozen orbs
                    bAllowed = .false.
                    exit
                endif
            enddo
        endif


        ! Check to see if UEG excitation is allowed, by summing kx, ky, kz
        ! over all the electrons
        if (tUEG .and. .not. tLatticeGens .and. bAllowed) then
            k = 0
            do i = 1, nel
                k = k + G1(nJ(i))%k
            enddo
            if (.not. all(k == 0)) &
                bAllowed = .false.
        endif


    end function CheckAllowedTruncSpawn



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
    
    SUBROUTINE CalcApproxpDoubles()
        use SystemData , only : tAssumeSizeExcitgen,tUseBrillouin,tNoSingExcits
        use CalcData , only : SinglesBias
        use SymData , only : SymClassSize
        use SymExcit3 , only : CountExcitations3
        INTEGER :: iTotal
        integer :: nSing, nDoub, ncsf, ierr
        integer :: hfdet_loc(nel)

        ! A quick hack. Count excitations as though we were a determinant.
        ! We could fix this later...
        hfdet_loc = iand(hfdet, csf_orbital_mask)

        ! TODO: A better approximation for ncsf.
        if (tCSF) then
            ncsf = 10
        else
            ncsf = 0
        endif
        nSing=0
        nDoub=0

        IF(tHub.or.tUEG) THEN
            IF(tReal) THEN
                WRITE(iout,*) "Since we are using a real-space hubbard model, only single excitations are connected &
                &   and will be generated."
                pDoubles=0.0_dp
                RETURN
            ELSE
                WRITE(iout,*) "Since we are using a momentum-space hubbard model/UEG, only double excitaitons &
     &                          are connected and will be generated."
                pDoubles=1.0_dp
                RETURN
            ENDIF
        elseif(tNoSingExcits) then
            write(iout,*) "Only double excitations will be generated"
            return
        ENDIF

!NSing=Number singles from HF, nDoub=No Doubles from HF

        WRITE(iout,"(A)") " Calculating approximate pDoubles for use with &
                       &excitation generator by looking a excitations from &
                       &reference."
        exflag=3
        IF(tKPntSym) THEN
            call enumerate_sing_doub_kpnt(exFlag, nSing, nDoub) 
        ELSE
            CALL CountExcitations3(HFDet_loc,exflag,nSing,nDoub)
        ENDIF
        iTotal=nSing + nDoub + ncsf

        WRITE(iout,"(I7,A,I7,A)") NDoub, " double excitations, and ",NSing, &
            " single excitations found from reference. This will be used to calculate pDoubles."

        IF(SinglesBias.ne.1.0_dp) THEN
            WRITE(iout,*) "Singles Bias detected. Multiplying single excitation connectivity of HF determinant by ", &
                SinglesBias," to determine pDoubles."
        ENDIF

        IF((NSing+nDoub+ncsf).ne.iTotal) THEN
            CALL Stop_All("CalcApproxpDoubles","Sum of number of singles and number of doubles does " &
            & //"not equal total number of excitations")
        ENDIF
        IF((NSing.eq.0).or.(NDoub.eq.0)) THEN
            WRITE(iout,*) "Number of singles or doubles found equals zero. pDoubles will be set to 0.95. Is this correct?"
            pDoubles = 0.95_dp
            pSingles = 0.05_dp
            RETURN
        elseif ((NSing < 0) .or. (NDoub < 0) .or. (ncsf < 0)) then
            call stop_all("CalcApproxpDoubles", &
                          "Number of singles, doubles or Yamanouchi symbols &
                          &found to be a negative number. Error here.")
        endif

        ! Set pDoubles to be the fraction of double excitations.
        ! If using CSFs, also consider only changing Yamanouchi Symbol
        if (tCSF) then
            pDoubles = real(nDoub,dp) / &
                   ((real(nSing,dp)*SinglesBias)+real(nDoub,dp)+real(ncsf,dp))
            pSingles = real(nSing,dp) / &
                   ((real(nSing,dp)*SinglesBias)+real(nDoub,dp)+real(ncsf,dp))

        else
            pDoubles = real(nDoub,dp) / &
                   ((real(NSing,dp)*SinglesBias) + real(NDoub,dp))
            pSingles = real(nSing,dp) * SinglesBias/ &
                   ((real(nSing,dp)*SinglesBias) + real(nDoub,dp))
        endif

        IF(SinglesBias.ne.1.0_dp) THEN
            write (iout, '("pDoubles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pDoubles, real(nDoub,dp) / real(iTotal,dp)
            write (iout, '("pSingles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pSingles, real(nSing,dp) / real(iTotal,dp)

!            WRITE(iout,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ", &
!                & real(nDoub,dp)/real(iTotal,dp)
        ELSE
            write (iout,'(A,F14.6)') " pDoubles set to: ", pDoubles
            write (iout,'(A,F14.6)') " pSingles set to: ", pSingles
        ENDIF

    END SUBROUTINE CalcApproxpDoubles

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

!Initialize the Histogramming searching arrays if necessary
    SUBROUTINE InitHistMin()
        IF(tHistSpawn.or.tCalcFCIMCPsi.and.(Iter.ge.NHistEquilSteps)) THEN
            IF(Iter.eq.NHistEquilSteps) THEN
                IF(iProcIndex.eq.Root) WRITE(iout,*) 'The iteration is equal to HISTEQUILSTEPS.  Beginning to histogram.'
            ENDIF
            HistMinInd(1:NEl)=FCIDetIndex(1:NEl)    !This is for the binary search when histogramming
        ENDIF
    END SUBROUTINE InitHistMin


    subroutine setup_linear_comb ()

        ! Take the specified spatial orbital, and configure to calculate the
        ! projected energy on a linear combination of these orbitals.
        !
        ! For now, weight the linear combination as per the current weightings
        ! --> works for from a pops file.

        ! --> Need current dets

        type(lexicographic_store) :: store
        integer(n_int) :: ilut_init(0:NIfTot), ilut_tmp(0:NIfTot)
        integer :: nopen, nup, pos, nfound
        real(dp) :: norm, sgn(lenof_sign)
        character(*), parameter :: t_r = 'setup_linear_comb'

        write(iout,*) 'Initialising projection onto linear combination of &
                   &determinants to calculate projected energy.'

        ! This currently only works with real walkers
        if (lenof_sign > 1) then
            call stop_all (t_r, 'Currently only available for real walkers')
            ! To enable for complex walkers, need to deal with complex
            ! amplitudes more sensibly in determining coeffs
        endif

        ! Get the initial ilut
        call EncodeBitDet (proje_ref_det_init, ilut_init)

        ! How many dets do we need?
        nopen = count_open_orbs(ilut_init)
        nup = (nopen + LMS) / 2
        nproje_sum = int(choose(nopen, nup))

        ! We only need a linear combination if there is more than one det...
        if (nproje_sum > 1) then

            ! TODO: log these.
            allocate(proje_ref_dets(nel, nproje_sum))
            allocate(proje_ref_iluts(0:NIfTot, nproje_sum))
            allocate(proje_ref_coeffs(nproje_sum))
            allocate(All_proje_ref_coeffs(nproje_sum))
            proje_ref_coeffs = 0

            ! Get all the dets
            nfound = 0
            call get_lexicographic_dets (ilut_init, store, ilut_tmp)
            do while (.not. all(ilut_tmp == 0))
                
                ! Store the ilut/det for later usage
                nfound = nfound + 1
                proje_ref_iluts(:,nfound) = ilut_tmp
                call decode_bit_det (proje_ref_dets(:,nfound), ilut_tmp)

                ! Find the ilut in CurrentDets, and use it to get coeffs
                ! TODO: 
                pos = binary_search(CurrentDets(:,1:TotWalkers), ilut_tmp, &
                                    NIfD+1)
                if (pos > 0) then
                    call extract_sign(CurrentDets(:,pos), sgn)
                    proje_ref_coeffs(nfound) = sgn(1)
                endif

                call get_lexicographic_dets (ilut_init, store, ilut_tmp)
            enddo

            if (nfound /= nproje_sum) &
                call stop_all (t_r, 'Incorrect number of determinants found')

            ! Get the total number of walkers on each site
            call MPISumAll(proje_ref_coeffs, All_proje_ref_coeffs)
            proje_ref_coeffs=All_proje_ref_coeffs
            
            norm = sqrt(sum(proje_ref_coeffs**2))
            if (norm == 0) norm = 1
            proje_ref_coeffs = proje_ref_coeffs / norm

        endif

    end subroutine setup_linear_comb

    subroutine update_linear_comb_coeffs ()

        integer :: i, pos, nfound, ierr
        real(dp) :: norm,reduce_in(1:2),reduce_out(1:2), sgn(lenof_sign)
                
        allocate(All_proje_ref_coeffs(nproje_sum))

        if(nproje_sum>1) then
           if(proje_spatial) then
              proje_ref_coeffs = 0
              do i = 1, nproje_sum
                 
                 pos = binary_search (CurrentDets(:,1:TotWalkers), &
                      proje_ref_iluts(:,i), NIfD+1)
                 if (pos > 0) then
                    call extract_sign (CurrentDets(:,pos), sgn)
                    proje_ref_coeffs(i) = sgn(1)
                 endif
              enddo
              
              call MPISumAll(proje_ref_coeffs, All_proje_ref_coeffs)
              proje_ref_coeffs=All_proje_ref_coeffs
              
              norm = sqrt(sum(proje_ref_coeffs**2))
              if (norm == 0) norm = 1
              proje_ref_coeffs = proje_ref_coeffs / norm
           else
              ! Finds the nproje_sum highest weight determinants in CurrentDets and
              ! forms a linear combination reference using the corresponding instantaneous weights.
              proje_linear_comb=.true.
              call clean_linear_comb()             
              allocate(proje_ref_dets(1:nel, 1:nproje_sum))
              allocate(proje_ref_iluts(0:NIfTot, 1:nproje_sum))
              allocate(proje_ref_coeffs(1:nproje_sum))
              proje_ref_coeffs=0
              proje_ref_dets=0
              proje_ref_iluts=0
              
              do pos=1,int(TotWalkers,sizeof_int)
                 call extract_sign(CurrentDets(:,pos),sgn)
                 if(abs(sgn(1))>abs(proje_ref_coeffs(nproje_sum))) then
                    inner: do nfound=1,nproje_sum
                       if(abs(sgn(1))>abs(proje_ref_coeffs(nfound))) then
                          proje_ref_coeffs((nfound+1):nproje_sum)=proje_ref_coeffs((nfound):(nproje_sum-1))
                          proje_ref_iluts(0:NIfTot,(nfound+1):nproje_sum)=proje_ref_iluts(0:NIfTot,(nfound):(nproje_sum-1))
                          proje_ref_coeffs(nfound)=sgn(1)
                          proje_ref_iluts(0:NIfTot,nfound)=CurrentDets(:,pos)
                          exit inner
                       endif
                    enddo inner
                 endif
              enddo
              
              do nfound=1,nproje_sum
                 reduce_in=(/abs(proje_ref_coeffs(nfound)),real(iProcIndex,dp)/)
                 CALL MPIAllReduceDatatype(reduce_in,1,MPI_MAXLOC,MPI_2DOUBLE_PRECISION,reduce_out)
                 IF(nint(reduce_out(2))/=iProcIndex) THEN 
                    proje_ref_coeffs((nfound+1):nproje_sum)=proje_ref_coeffs(nfound:(nproje_sum-1))
                    proje_ref_iluts(0:NIfTot,(nfound+1):nproje_sum)=proje_ref_iluts(0:NIfTot,(nfound):(nproje_sum-1))
                 ENDIF
                 CALL MPIBCast(proje_ref_coeffs(nfound),1,nint(reduce_out(2)))
                 CALL MPIBCast(proje_ref_iluts(0:NIfTot,nfound),NIfTot+1,nint(reduce_out(2)))
              enddo
              do nfound=1,nproje_sum
                 call decode_bit_det(proje_ref_dets(:,nfound),proje_ref_iluts(0:NIfTot,nfound))
              enddo
              norm = sqrt(sum(proje_ref_coeffs**2))
              if (norm == 0) norm = 1
              proje_ref_coeffs = proje_ref_coeffs / norm
              proje_update_comb=.false.
              write(iout,*) ' Linear combination coeffs successfully updated. '
              write(iout,*) ' Number of linear combination coeffs is now',nproje_sum
              write(iout,*) ' Coeff | Determinant '
              do nfound=1,nproje_sum
                 write(iout,*) nfound, proje_ref_coeffs(nfound),'DET:',proje_ref_dets(:,nfound)
              enddo
              write(iout,*) ' end of linear combination list'
           endif
        endif
    end subroutine

    subroutine clean_linear_comb ()

        if (allocated(proje_ref_dets)) &
            deallocate(proje_ref_dets)
        if (allocated(proje_ref_iluts)) &
            deallocate(proje_ref_iluts)
        if (allocated(proje_ref_coeffs)) &
            deallocate(proje_ref_coeffs)
        if (allocated(All_proje_ref_coeffs)) &
            deallocate(All_proje_ref_coeffs)

    end subroutine clean_linear_comb

    ! This routine sums in the energy contribution from a given walker and 
    ! updates stats such as mean excit level AJWT added optional argument 
    ! dProbFin which is a probability that whatever gave this contribution 
    ! was generated. It defaults to 1, and weights the contribution of this 
    ! det (only in the projected energy) by dividing its contribution by this
    ! number 
    subroutine SumEContrib (nI, ExcitLevel, RealWSign, ilut, HDiagCurr, dProbFin, ind)

        integer, intent(in) :: nI(nel), ExcitLevel
        real(dp), intent(in) :: RealwSign(lenof_sign)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: HDiagCurr, dProbFin
        integer, intent(in), optional :: ind

        integer :: i, bin, pos, ExcitLevel_local, ExcitLevelSpinCoup, run
        integer :: PartInd, OpenOrbs, spatial_ic
        integer(n_int) :: iLutSym(0:NIfTot)
        logical tSuccess
        integer :: iUEG1, iUEG2, ProjEBin
        HElement_t :: HOffDiag
        HElement_t :: HDoubDiag
        integer :: DoubEx(2,2),DoubEx2(2,2),kDoub(3) ! For histogramming UEG doubles
        integer :: ExMat(2,2),FindSplitProjEBinG,FindSplitProjEBinK3
        logical :: tDoubParity,tDoubParity2,tSign ! As above

        ! Are we performing a linear sum over various determinants?
        ! TODO: If we use this, function pointer it.

        HOffDiag = 0

        ! Add in the contributions to the numerator and denominator of the trial
        ! estimator, if it is being used.
        if (tTrialWavefunction .and. present(ind)) then
            if (tHashWalkerlist) then
                if (test_flag(ilut, flag_trial)) then
                    do run=1,inum_runs
                        trial_denom(run) = trial_denom(run) + current_trial_amps(ind)*RealwSign(run)
                    enddo
                else if (test_flag(ilut, flag_connected)) then
                    do run=1,inum_runs
                        trial_numerator(run) = trial_numerator(run) + current_trial_amps(ind)*RealwSign(run)
                    enddo
                end if
            else
                if (test_flag(ilut, flag_trial)) then
                    ! Take the next element in the occupied trial vector.
                    trial_ind = trial_ind + 1
                    do run=1,inum_runs
                        trial_denom(run) = trial_denom(run) + occ_trial_amps(trial_ind)*RealwSign(run)
                    enddo
                else if (test_flag(ilut, flag_connected)) then
                    ! Take the next element in the occupied connected vector.
                    con_ind = con_ind + 1
                    do run=1,inum_runs
                        trial_numerator(run) = trial_numerator(run) + occ_con_amps(con_ind)*RealwSign(run)
                    enddo
                end if
            end if
        end if

        if (proje_linear_comb .and. nproje_sum > 1) then
           if(proje_spatial) then
              spatial_ic = FindSpatialBitExcitLevel (ilut, proje_ref_iluts(:,1))
              if (spatial_ic <= 2) then
                 do i = 1, nproje_sum
                    if (proje_ref_coeffs(i) /= 0) then
                       HOffDiag = HOffDiag + proje_ref_coeffs(i) &
                            * get_helement (proje_ref_dets(:,i), nI, &
                            proje_ref_iluts(:,i), ilut)
                    endif
                 enddo
              endif
           else
              do i=1, nproje_sum
                 !this is thought for a spin-polarized system.
                 spatial_ic = FindBitExcitLevel (ilut, proje_ref_iluts(:,i))
                 if (spatial_ic <= 2) then
                    if (proje_ref_coeffs(i) /= 0._dp) then
                       IF(spatial_ic==0) THEN
                          HOffDiag = HOffDiag + proje_ref_coeffs(i) &
                               * (get_helement (proje_ref_dets(:,i),nI,spatial_ic, &
                               proje_ref_iluts(:,i), ilut)-Hii)
                       ELSE
                          HOffDiag = HOffDiag + proje_ref_coeffs(i) &
                               * get_helement (proje_ref_dets(:,i),nI,spatial_ic, &
                               proje_ref_iluts(:,i), ilut)
                       ENDIF
                    endif
                 endif
              enddo  
           endif
        else
            ! ExcitLevel indicates the excitation level between the det and
            ! *one* of the determinants in an HPHF/MomInv function. If needed,
            ! calculate the connection between it and the other one. If either
            ! is connected, then it has to be counted. Since the excitation
            ! level is the same to either det, we don't need to consider the
            ! spin-coupled det of both reference and current HPHFs.
            !
            ! For determinants, set ExcitLevel_local == ExcitLevel.
            ExcitLevel_local = ExcitLevel
            if (tSpinCoupProjE .and. (ExcitLevel /= 0)) then
                ExcitLevelSpinCoup = FindBitExcitLevel (iLutRefFlip, &
                                                        ilut, 2)
                if (ExcitLevelSpinCoup <= 2 .or. ExcitLevel <= 2) &
                    ExcitLevel_local = 2
            endif

            ! Perform normal projection onto reference determinant
            if (ExcitLevel_local == 0) then
                
                if (iter > NEquilSteps) SumNoatHF = SumNoatHF + RealwSign
                NoatHF = NoatHF + RealwSign
                ! Number at HF * sign over course of update cycle
                HFCyc = HFCyc + RealwSign

            elseif (ExcitLevel_local == 2 .or. &
                    (ExcitLevel_local == 1 .and. tNoBrillouin)) then

                ! For the real-space Hubbard model, determinants are only
                ! connected to excitations one level away, and Brillouins
                ! theorem cannot hold.
                !
                ! For Rotated orbitals, Brillouins theorem also cannot hold,
                ! and energy contributions from walkers on singly excited
                ! determinants must also be included in the energy values
                ! along with the doubles
                
                if (ExcitLevel == 2) then
#ifdef __CMPLX
                    NoatDoubs(1) = NoatDoubs(1) + sum(abs(RealwSign))
#else           
                    do run=1, inum_runs
                        NoatDoubs(run) = NoatDoubs(run) + abs(RealwSign(run))
                    enddo
#endif
                endif
                ! Obtain off-diagonal element
                if (tHPHF) then
                    HOffDiag = hphf_off_diag_helement (ProjEDet, nI, iLutRef,&
                                                       ilut)
                elseif(tMomInv) then
                    HOffDiag = MI_off_diag_helement (ProjEDet, nI, iLutRef, ilut)
                else
                    HOffDiag = get_helement (ProjEDet, nI, ExcitLevel, &
                                             ilutRef, ilut)
                endif

            endif ! ExcitLevel_local == 1, 2, 3

        endif ! sume_linear_contrib

        ! Sum in energy contribution
        do run=1, inum_runs
            if (iter > NEquilSteps) &
                SumENum(run) = SumENum(run) + (HOffDiag * ARR_RE_OR_CPLX(RealwSign,run)) / dProbFin
       
            ENumCyc(run) = ENumCyc(run) + (HOffDiag * ARR_RE_OR_CPLX(RealwSign,run)) / dProbFin
            ENumCycAbs(run) = ENumCycAbs(run) + abs(HOffDiag * ARR_RE_OR_CPLX(RealwSign,run)) / dProbFin
        enddo
        

        ! -----------------------------------
        ! HISTOGRAMMING
        ! -----------------------------------

        if ((tHistSpawn .or. (tCalcFCIMCPsi .and. tFCIMC)) .and. &
            (iter >= NHistEquilSteps)) then
            ! Histogram particles by determinant
            call add_hist_spawn (ilut, RealwSign, ExcitLevel_local, dProbFin)
        elseif (tHistEnergies) then
            ! Histogram particles by energy
            call add_hist_energies (ilut, RealwSign, HDiagCurr, ExcitLevel)
        endif

        ! Are we doing a spin-projection histogram?
        if (tHistSpinDist) then
            if (tUseRealCoeffs) call stop_all("SumEContrib", "Not set up to use real coeffs with tHistSpinDist")
            call test_add_hist_spin_dist_det (ilut, RealwSign)
        endif

        ! Maintain a list of the degree of occupation of each orbital
        if (tPrintOrbOcc .and. (iter >= StartPrintOrbOcc)) then
            if((Iter.eq.StartPrintOrbOcc).and.(DetBitEQ(ilut,ilutHF,NIfDBO))) then
                write(6,*) 'Beginning to fill the HF orbital occupation list &
                                &during iteration',Iter
                if(tPrintOrbOccInit) &
                    write(6,*) 'Only doing so for initiator determinants.'
            endif
            if ((tPrintOrbOccInit .and. test_flag(ilut,flag_is_initiator(1)))&
                .or. .not. tPrintOrbOccInit) then
                forall (i = 1:nel) OrbOccs(nI(i)) = OrbOccs(nI(i)) &
                                          + (RealwSign(1) * RealwSign(1))
            endif
        endif
        
        if (tPrintDoubsUEG) then
            if (Iter.ge.StartPrintDoubsUEG) then
                if (ExcitLevel.eq.2) then
                    DoubEx2=0
                    DoubEx2(1,1)=2
                    call GetExcitation (ProjEDet,nI,NEl,DoubEx2,tDoubParity2)
                    DoubEx=0
                    DoubEx(1,1)=2
                    call GetBitExcitation(iLutRef,ilut,DoubEx,tDoubParity)
                    if (DoubEx2(1,1).ne.DoubEx(1,1) &
                        .or. DoubEx2(1,2).ne.DoubEx(1,2) &
                        .or. DoubEx2(2,2).ne.DoubEx(2,2) &
                        .or. DoubEx2(2,1).ne.DoubEx(2,1) &
                        .or. tDoubParity.neqv.tDoubParity2) then
                        call stop_all("SumEContrib","GetBitExcitation doesn't agree with GetExcitation")
                    endif
                    iUEG1=0
                    iUEG2=0
                    iUEG1=int(DoubsUEGLookup(DoubEx(1,1)),sizeof_int)
                    iUEG2=int(DoubsUEGLookup(DoubEx(1,2)),sizeof_int)
                    if (iUEG1.eq.0.or.iUEG2.eq.0) call stop_all("SumEContrib","Array bounds issue")
                    DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),1)=DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),1)+RealWSign(1)
    ! Test against natural orbital generation. For a two electron system, this should just be the same 
    ! as the nat orbs if WSign is squared
    !                    DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),1)=DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),1)+(REAL(WSign(1))*REAL(WSign(1)))
                    if (DoubsUEGStore(iUEG1,iUEG2,DoubEx(2,1))) then
                        DoubsUEGStore(iUEG1,iUEG2,DoubEx(2,1))=.false.
                        DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),2)=HOffDiag
                        if(tHPHF) then
                            HDoubDiag = hphf_diag_helement (nI,ilut)
                        elseif(tMomInv) then
                            HDoubDiag = MI_diag_helement(nI,ilut)
                        else
                            HDoubDiag = get_helement (nI,nI,0) !, iLutCurr, &
                                                    ! iLutCurr)
                        endif
                        DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),3)=HDoubDiag
                        kDoub=0
                        kDoub=G1(DoubEx(2,1))%k
                        DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),4)=REAL(kDoub(1)*kDoub(1),dp)+ &
                            REAL(kDoub(2)*kDoub(2),dp)+REAL(kDoub(3)*kDoub(3),dp)
                    endif
                endif
            endif
        endif

    end subroutine SumEContrib


!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        use FciMCLoggingMOD , only : InitHistInitPops
        use SystemData , only : tRotateOrbs
        use CalcData , only : InitialPart,tstartmp1,tStartCAS, InitialPartVec
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,iReadWalkersRoot
        use constants , only : size_n_int
        use DeterminantData , only : write_det
        use nElRDMMod, only: InitRDM
        use LoggingData , only : tReadRDMs
        INTEGER :: ierr,iunithead,DetHash,Slot,MemTemp,run
        LOGICAL :: formpops,binpops
        INTEGER :: error,MemoryAlloc,PopsVersion,j,iLookup,WalkerListSize
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'
        integer :: ReadBatch    !This parameter determines the length of the array to batch read in walkers from a popsfile
        integer :: PopBlockingIter
        real(dp) :: Gap,ExpectedMemWalk,read_tau, read_psingles, read_par_bias
        integer(int64) :: read_walkers_on_nodes(0:nProcessors-1)
        integer :: read_nnodes
        !Variables from popsfile header...
        logical :: tPop64Bit,tPopHPHF,tPopLz
        integer :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,Popinum_runs
        integer(int64) :: iPopAllTotWalkers
        integer :: i
        real(dp) :: PopDiagSft, PopDiagSft2
        real(dp) , dimension(lenof_sign) :: InitialSign
        real(dp) , dimension(lenof_sign/inum_runs) :: PopSumNoatHF
        HElement_t :: PopAllSumENum
        
        !default
        Popinum_runs=1

        if(tReadPops.and..not.tPopsAlreadyRead) then
            call open_pops_head(iunithead,formpops,binpops)
            PopsVersion=FindPopsfileVersion(iunithead)
            if(iProcIndex.eq.root) close(iunithead)
            write(iout,*) "POPSFILE VERSION ",PopsVersion," detected."
        endif

        if(tPopsMapping.and.(PopsVersion.lt.3)) then
            write(iout,*) "Popsfile mapping cannot work with old POPSFILEs"
            call stop_all("InitFCIMCCalcPar","Popsfile mapping cannot work with old POPSFILEs")
        endif

        ! Initialise measurement of norm, to avoid divide by zero
        norm_psi = 1.0_dp

        if (tReadPops .and. (PopsVersion.lt.3) .and..not.tPopsAlreadyRead) then
!Read in particles from multiple POPSFILES for each processor
            !Ugh - need to set up ValidSpawnedList here too...
            call SetupValidSpawned(int(InitWalkers,int64))
            WRITE(iout,*) "Reading in initial particle configuration from *OLD* POPSFILES..."
            CALL ReadFromPopsFilePar()
        ELSE
!initialise the particle positions - start at HF with positive sign
!Set the maximum number of walkers allowed
            if(tReadPops.and..not.tPopsAlreadyRead) then
                !Read header.
                call open_pops_head(iunithead,formpops,binpops)
                if(PopsVersion.eq.3) then
                    call ReadPopsHeadv3(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                            iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                            PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)

                    ! The following values were not read in...
                    read_tau = 0.0_dp
                    read_nnodes = 0
                elseif(PopsVersion.eq.4) then
                    call ReadPopsHeadv4(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                            iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                            PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot, &
                            read_tau,PopBlockingIter, read_psingles, &
                            read_par_bias, read_nnodes, read_walkers_on_nodes)
                    ! The only difference between 3 & 4 is just that 4 reads 
                    ! in via a namelist, so that we can add more details 
                    ! whenever we want.
                else
                    call stop_all(this_routine,"Popsfile version invalid")
                endif

                call CheckPopsParams(tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                        iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                        PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot, &
                        WalkerListSize,read_tau,PopBlockingIter, &
                        read_psingles, read_par_bias)

                if(iProcIndex.eq.root) close(iunithead)
            else
                WalkerListSize=int(InitWalkers,sizeof_int)
            endif

            MaxWalkersPart=NINT(MemoryFacPart*WalkerListSize)
            ExpectedMemWalk=real((NIfTot+1)*MaxWalkersPart*size_n_int+8*MaxWalkersPart,dp)/1048576.0_dp
            if(ExpectedMemWalk.lt.20.0) then
                !Increase memory allowance for small runs to a min of 20mb
                MaxWalkersPart=int(20.0*1048576.0/real((NIfTot+1)*size_n_int+8,dp),sizeof_int)
                write(iout,"(A)") "Low memory requested for walkers, so increasing memory to 20Mb to avoid memory errors"
            endif
            WRITE(iout,"(A,I14)") "Memory allocated for a maximum particle number per node of: ",MaxWalkersPart
            !Here is where MaxSpawned is set up - do we want to set up a minimum allocation here too?
            Call SetupValidSpawned(int(InitWalkers, int64))

!Put a barrier here so all processes synchronise
            CALL MPIBarrier(error)
!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
            MemoryAlloc=(NIfTot+1)*MaxWalkersPart*size_n_int    !Memory Allocated in bytes

            IF(.not.tRegenDiagHEls) THEN
                if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
                    ! If calculating the RDMs stochastically, need to include the average sign and the 
                        ! iteration it became occupied in the CurrentH array (with the Hii elements).
                        ! Need to store this for each set of walkers if we're doing a double run.
                        ALLOCATE(WalkVecH(1+2*lenof_sign,MaxWalkersPart),stat=ierr)
                        CALL LogMemAlloc('WalkVecH',(1+2*lenof_sign)*MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                        WalkVecH(:,:)=0.0_dp
                        MemoryAlloc=MemoryAlloc+8*MaxWalkersPart*(1+2*lenof_sign)
                        WRITE(6,"(A)") " The average current signs before death will be stored for use in the RDMs."
                        WRITE(6,"(A,F14.6,A)") " This requires ", &
                            real(MaxWalkersPart*8*2*lenof_sign,dp)/1048576.0_dp, &
                            " Mb/Processor"
                    NCurrH = 1+2*lenof_sign
                else
                    ALLOCATE(WalkVecH(1,MaxWalkersPart),stat=ierr)
                    CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                    WalkVecH(:,:)=0.0_dp
                    MemoryAlloc=MemoryAlloc+8*MaxWalkersPart
                    NCurrH = 1
                endif
            ELSE
                WRITE(iout,"(A,F14.6,A)") " Diagonal H-Elements will not be stored. This will *save* ", &
                    & REAL(MaxWalkersPart*8,dp)/1048576.0_dp," Mb/Processor"
            ENDIF

            if(tRDMonFly.and.(.not.tExplicitAllRDM).and.(.not.tHF_Ref_Explicit)) then
                !Allocate memory to hold walkers spawned from one determinant at a time.
                !Walkers are temporarily stored here, so we can check if we're spawning onto the same Dj multiple times.
                !If only using connections to the HF (tHF_Ref_Explicit), no stochastic RDM construction is done, and this 
                !is not necessary.
                ALLOCATE(TempSpawnedParts(0:NIfDBO,20000),stat=ierr)
                CALL LogMemAlloc('TempSpawnedParts',20000*(NIfDBO+1),size_n_int,this_routine,TempSpawnedPartsTag,ierr)
                TempSpawnedParts(0:NIfDBO,1:20000)=0
                MemoryAlloc=MemoryAlloc + (NIfDBO+1)*20000*size_n_int    !Memory Allocated in bytes
                WRITE(6,"(A)") " Allocating temporary array for walkers spawned from a particular Di."
                WRITE(6,"(A,F14.6,A)") " This requires ", REAL(((NIfDBO+1)*20000*size_n_int),dp)/1048576.0_dp," Mb/Processor"
            endif
            
            WRITE(iout,"(A,I12,A)") "Spawning vectors allowing for a total of ",MaxSpawned, &
                    " particles to be spawned in any one iteration per core."
            ALLOCATE(SpawnVec(0:NIftot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVecTag,ierr)
            ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)

            SpawnVec(:,:)=0
            SpawnVec2(:,:)=0

!Point at correct spawning arrays
            SpawnedParts=>SpawnVec
            SpawnedParts2=>SpawnVec2

            MemoryAlloc=MemoryAlloc+(NIfTot+1)*MaxSpawned*2*size_n_int

            if(tHashWalkerList) then
                write(iout,"(A)") "Storing walkers in hash-table. Algorithm is now formally linear scaling with walker number"
                write(iout,"(A,I15)") "Length of hash-table: ",nWalkerHashes
                write(iout,"(A,F20.5)") "Length of hash-table as a fraction of targetwalkers: ",HashLengthFrac
                ! TODO: Correct the memory usage.
                !MemTemp=2*(8*(nClashMax+1)*nWalkerHashes)+8*MaxWalkersPart
                !write(iout,"(A,F10.3,A)") "This will use ",real(MemTemp,dp)/1048576.0_dp,&
                !    " Mb of memory per process, although this is likely to increase as it expands"
                allocate(HashIndex(nWalkerHashes),stat=ierr)
                if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")
                do i = 1, nWalkerHashes
                    HashIndex(i)%ind = 0
                end do
                !Also need to allocate memory for the freeslot array
                allocate(FreeSlot(MaxWalkersPart),stat=ierr)
                if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")
                freeslot(:)=0
                !MemoryAlloc=MemoryAlloc+MemTemp
            endif

!Allocate pointers to the correct walker arrays
            CurrentDets=>WalkVecDets
            IF(.not.tRegenDiagHEls) THEN
                CurrentH=>WalkVecH
            ENDIF

            ! Get the (0-based) processor index for the HF det.
            iHFProc = DetermineDetNode(HFDet,0)
            WRITE(iout,"(A,I8)") "Reference processor is: ",iHFProc
            write(iout,"(A)",advance='no') "Initial reference is: "
            call write_det(iout,HFDet,.true.)

            TotParts(:)=0.0
            TotPartsOld(:)=0.0
            NoatHF(:)=0.0_dp
            InstNoatHF(:)=0.0_dp
            AvNoatHF(:) = 0.0_dp

            if (tSpawnSpatialInit) then
                max_inits = int(MemoryFacInit * INitWalkers)
                no_spatial_init_dets = 0
                allocate(CurrentInits(0:nIfTot, max_inits), stat=ierr)
                call LogMemAlloc('CurrentInits', max_inits * (NIfTot+1), &
                                 size_n_int, this_routine, CurrentInitTag, &
                                 ierr)
            endif

            !
            ! If we have a popsfile, read the walkers in now.
            if(tReadPops.and..not.tPopsAlreadyRead) then

                if(iReadWalkersRoot.eq.0) then

                    ! ReadBatch is the number of walkers to read in from the 
                    ! popsfile at one time. The larger it is, the fewer
                    ! communictions will be needed to scatter the particles.
                    !
                    ! By default, the new array (which is only created on the 
                    ! root processors) is the same length as the spawning 
                    ! arrays.
                    ReadBatch=MaxSpawned
                else
                    ReadBatch = iReadWalkersRoot
                endif

                ! TotWalkers and TotParts are returned as the dets and parts 
                ! on each processor.
                call ReadFromPopsfile(iPopAllTotWalkers, ReadBatch, &
                                      TotWalkers ,TotParts, NoatHF, &
                                      CurrentDets, MaxWalkersPart, &
                                      read_nnodes, read_walkers_on_nodes, &
                                      PopNIfSgn)

                !Setup global variables
                TotWalkersOld=TotWalkers
                TotPartsOld = TotParts
                call MPISumAll(TotWalkers,AllTotWalkers)
                AllTotWalkersOld = AllTotWalkers
                call MPISumAll(TotParts,AllTotParts)
                AllTotPartsOld=AllTotParts
                call MPISumAll(NoatHF,AllNoatHF)
                OldAllNoatHF=AllNoatHF
#ifdef __CMPLX
                OldAllAvWalkersCyc=sum(AllTotParts)
#else
                OldAllAvWalkersCyc=AllTotParts
#endif
                
                do run=1,inum_runs
                    OldAllHFCyc(run) = ARR_RE_OR_CPLX(AllNoatHF,run)
                enddo
                
                AllNoAbortedOld(:)=0.0_dp
                iter_data_fciqmc%tot_parts_old = AllTotParts
                
                ! Calculate the projected energy for this iteration.
                do run=1,inum_runs
                    if (ARR_RE_OR_CPLX(AllSumNoAtHF,run)/=0) &
                        ProjectionE(run) = AllSumENum(run) / ARR_RE_OR_CPLX(AllSumNoatHF,run)
                enddo 

                if(iProcIndex.eq.iHFProc) then
                    !Need to store SumENum and SumNoatHF, since the global variable All... gets wiped each iteration. 
                    !Rather than POPSFILE v2, where the average values were scattered, just store the previous
                    !energy contributions on the root node.
                    SumNoatHF(:)=AllSumNoatHF(:)
                    SumENum(:)=AllSumENum(:)
                    InstNoatHF(:) = NoatHF(:)

                    if((AllNoatHF(1).ne.NoatHF(1)).or.(AllNoatHF(lenof_sign).ne.NoatHF(lenof_sign))) then
                        call stop_all(this_routine,"HF particles spread across different processors.")
                    endif
                endif
            
            else

                if(tStartMP1) then
                    !Initialise walkers according to mp1 amplitude.
                    call InitFCIMC_MP1()

                elseif(tStartCAS) then
                    !Initialise walkers according to a CAS diagonalisation.
                    call InitFCIMC_CAS()

                else !Set up walkers on HF det

                    if(tStartSinglePart) then
                        WRITE(iout,"(A,I10,A,F9.3,A,I15)") "Initial number of particles set to ",int(InitialPart), &
                            " and shift will be held at ",DiagSft(1)," until particle number gets to ", int(InitWalkers*nNodes)
                    else
                        write(iout,"(A,I16)") "Initial number of walkers per processor chosen to be: ", InitWalkers
                    endif
                   
                    do run=1,inum_runs
                        InitialPartVec(run)=InitialPart
                    enddo

                    !Setup initial walker local variables for HF walkers start
                    IF(iProcIndex.eq.iHFProc) THEN

                        ! Encode the reference determinant identification.
                        call encode_det(CurrentDets(:,1), iLutHF)
                        if(tHashWalkerList) then
                            !Point at the correct position for the first walker
                            DetHash=FindWalkerHash(HFDet, nWalkerHashes)    !Find det hash position
                            HashIndex(DetHash)%Ind = 1
                        endif

                        ! Clear the flags
                        call clear_all_flags (CurrentDets(:,1))

                        ! Set reference determinant as an initiator if
                        ! tTruncInitiator is set, for both imaginary and real flags
                        if (tTruncInitiator) then
                            call set_flag (CurrentDets(:,1), flag_is_initiator(1))
                            call set_flag (CurrentDets(:,1), flag_is_initiator(2))
                            if (tSpawnSpatialInit) &
                                call add_initiator_list (CurrentDets(:,1))
                        endif

                        ! If running a semi-stochastic simulation, set flag to specify the Hartree-Fock is in the
                        ! deterministic space.
                        if (tSemiStochastic) call set_flag (CurrentDets(:,1), flag_deterministic)

                        ! HF energy is equal to 0 (by definition)
                        if (.not. tRegenDiagHEls) CurrentH(1,1) = 0
                        HFInd = 1

                        ! Obtain the initial sign
                        InitialSign = 0
                        if (tStartSinglePart) then
                            InitialSign(:) = InitialPartVec(:)
                            TotParts(:) = InitialPartVec(:)
                            TotPartsOld(:) = InitialPartVec(:)
                        else
                            do run=1, inum_runs
                                InitialSign(run) = InitWalkers
                                TotParts(run) = real(InitWalkers,dp)
                                TotPartsOld(run) = real(InitWalkers,dp)
                            enddo
                        endif

                        ! set initial values for global control variables.
                        
                        TotWalkers = 1
                        TotWalkersOld = 1
                        NoatHF(:) = InitialSign(:)
                        call encode_sign (CurrentDets(:,1), InitialSign)
                    ELSE
                        NoatHF(:) = 0.0_dp
                        TotWalkers = 0.0_dp
                        TotWalkersOld = 0.0_dp
                    ENDIF

                    OldAllNoatHF(:)=0.0_dp
                    AllNoatHF(:)=0.0_dp
                    IF(TStartSinglePart) THEN
        !Initialise global variables for calculation on the root node
                        IF(iProcIndex.eq.root) THEN
                            OldAllNoatHF=InitialPartVec
                            do run=1,inum_runs
                                OldAllAvWalkersCyc(run) = InitialPartVec(run)
                            enddo
                            AllNoatHF=InitialPartVec
                            InstNoatHF = InitialPartVec
                            AllTotParts=InitialPartVec
                            AllTotPartsOld=InitialPartVec
                            AllNoAbortedOld(:)=0.0_dp
                            iter_data_fciqmc%tot_parts_old = InitialPartVec
                            AllTotWalkers = 1
                            AllTotWalkersOld = 1
                            do run=1,inum_runs
                                OldAllHFCyc(run) = ARR_RE_OR_CPLX(InitialPartVec,run)
                            enddo
                        ENDIF
                    ELSE
        !In this, only one processor has initial particles.
                        IF(iProcIndex.eq.Root) THEN
                            AllTotWalkers = 1
                            AllTotWalkersOld = 1
                            do run=1,inum_runs
                                iter_data_fciqmc%tot_parts_old(run) = real(InitWalkers,dp)
                                AllTotParts(run)=InitWalkers
                                AllTotPartsOld(run)=InitWalkers
                                AllNoAbortedOld(run)=0.0_dp
                            enddo
                        ENDIF
                    ENDIF

                endif   !tStartmp1
            endif  
        
            WRITE(iout,"(A,F14.6,A)") " Initial memory (without excitgens + temp arrays) consists of : ", &
                & REAL(MemoryAlloc,dp)/1048576.0_dp," Mb/Processor"
            WRITE(iout,*) "Only one array of memory to store main particle list allocated..."
            WRITE(iout,*) "Initial memory allocation sucessful..."
            CALL neci_flush(iout)

        ENDIF   !End if initial walkers method

        ! If we are projecting onto a linear combination to calculate projE,
        ! then do the setup
        if (proje_spatial .and. proje_linear_comb) &
             call setup_linear_comb ()

            
!Put a barrier here so all processes synchronise
        CALL MPIBarrier(error)

        IF(tTruncInitiator.or.tDelayTruncInit) THEN
            IF(tDelayTruncInit) tTruncInitiator=.false.
        ENDIF

        IF(tPrintOrbOcc) THEN
            ALLOCATE(OrbOccs(nBasis),stat=ierr)
            CALL LogMemAlloc('OrbOccs',nBasis,8,this_routine,OrbOccsTag,ierr)
            OrbOccs(:)=0.0_dp
        ENDIF

        IF(tPrintDoubsUEG) THEN
            ALLOCATE(DoubsUEG(NEl,NEl,nBasis,4),stat=ierr)
            DoubsUEG(:,:,:,:)=0.0_dp
            ALLOCATE(DoubsUEGLookup(nBasis),stat=ierr)
            DoubsUEGLookup(:)=0
            ALLOCATE(DoubsUEGStore(NEl,NEl,nBasis),stat=ierr)
            DoubsUEGStore(:,:,:)=.true.
!Add LogMemAllocs
            do iLookup=1,NEl                    
                DoubsUEGLookup(HFDet(iLookup))=iLookup
            enddo
        ENDIF

        IF(tHistInitPops) THEN
            CALL InitHistInitPops()
        ENDIF
        tPrintHighPop=.false.
        MaxInitPopPos=0.0
        MaxInitPopNeg=0.0

        IF(MaxNoatHF.eq.0) THEN
            MaxNoatHF=InitWalkers*nNodes
            HFPopThresh=int(MaxNoatHF,int64)
        ENDIF

        ! Initialise excitation generation storage
        call init_excit_gen_store (fcimc_excit_gen_store)

        IF((NMCyc.ne.0).and.(tRotateOrbs.and.(.not.tFindCINatOrbs))) then 
            CALL Stop_All(this_routine,"Currently not set up to rotate and then go straight into a spawning &
            & calculation.  Ordering of orbitals is incorrect.  This may be fixed if needed.")
        endif
        
        if(tSpawn_Only_Init .and. (.not.tTruncInitiator)) then
            CALL Stop_All(this_routine,"Cannot use the SPAWNONLYINIT option without the TRUNCINITIATOR option.")
        endif

        if (tSpinProject) then
            if (inum_runs.eq.2) call stop_all(this_routine,"Code not yet set up to do a double run &
                    & with tSpinProject. E.g. when calling the main loop, tSinglePartPhase is now length 2")
            call init_yama_store ()
        endif

        IF(tRDMonFly) CALL InitRDM()
        !This keyword (tRDMonFly) is on from the beginning if we eventually plan to calculate the RDM's.
        !Initialises RDM stuff for both explicit and stochastic calculations of RDM.

        tFillingStochRDMonFly = .false.      
        tFillingExplicRDMonFly = .false.      
        !One of these becomes true when we have reached the relevant iteration to begin filling the RDM.

        !If the iteration specified to start filling the RDM has already been, want to 
        !start filling as soon as possible.
        if(tRDMonFly) then
            do run=1,inum_runs
                if(.not.tSinglePartPhase(run)) VaryShiftIter(run) = 0
            enddo
        endif

        ! Perform all semi-stochastic initialisation. This includes generating all the states in the
        ! deterministic space, finding their processors, ordering them, inserting them into
        ! CurrentDets, calculating and storing all Hamiltonian matrix elements and initalising all
        ! arrays required to store and distribute the vectors in the deterministic space later.
        if (tSemiStochastic) then
            if(inum_runs.eq.2 .and. tRDMonFly) call stop_all('InitFCIMCCalcPar','Cannot yet do SS and double &
                    &run simultaneously with RDMs. In addition to other things, need to deal with &
                    &fill_RDM_offdiag_deterministic().')
            call init_semi_stochastic()
        endif


       if (tSpawnSpatialInit .and. (inum_runs.eq.2)) call stop_all('InitFCIMCCalcPar', &
                    & "Double run not set up to use with tSpawnSpatialInit.  e.g. likely problem with rm_initiator_list")

        ! Initialise the trial wavefunction information which can be used for the energy estimator.
        ! This includes generating the trial space, generating the space connected to the trial space,
        ! diagonalising the trial space to find the trial wavefunction and calculating the vector
        ! in the connected space, required for the energy estimator.
        if (tTrialWavefunction) call init_trial_wf()
        
    end subroutine InitFCIMCCalcPar

!Routine to initialise the particle distribution according to a CAS diagonalisation. 
!This hopefully will help with close-lying excited states of the same sym.
    subroutine InitFCIMC_CAS()
        use SystemData, only : tSpn,tHPHFInts
        use CalcData, only: InitialPart
        use DeterminantData, only : write_det,write_det_len 
        use DetCalcData, only : NKRY,NBLK,B2L,nCycle
        use DetBitOps, only: FindBitExcitLevel
        use sym_mod , only : Getsym, writesym
        use MomInv, only: IsAllowedMI 
        type(BasisFN) :: CASSym
        integer :: i, j, ierr, nEval, NKRY1, NBLOCK, LSCR, LISCR, DetIndex
        integer :: iNode, nBlocks, nBlockStarts(2), DetHash, Slot
        integer :: CASSpinBasisSize, elec, nCASDet, ICMax, GC, LenHamil, iInit
        integer :: nHPHFCAS, iCasDet, ExcitLevel
        real(dp) :: NoWalkers
        integer , allocatable :: CASBrr(:),CASDet(:),CASFullDets(:,:),nRow(:),Lab(:),ISCR(:),INDEX(:)
        integer , pointer :: CASDetList(:,:) => null()
        integer(n_int) :: iLutnJ(0:NIfTot)
        logical :: tMC
        HElement_t :: HDiagTemp
        real(dp) , allocatable :: CK(:,:),W(:),CKN(:,:),Hamil(:),A_Arr(:,:),V(:),BM(:),T(:),WT(:)
        real(dp) , allocatable :: SCR(:),WH(:),Work2(:),V2(:,:),AM(:)
        real(dp) , allocatable :: Work(:)
        integer(TagIntType) :: ATag=0,VTag=0,BMTag=0,TTag=0,WTTag=0,SCRTag=0,WHTag=0,Work2Tag=0,V2Tag=0
        integer(TagIntType) :: ISCRTag=0,IndexTag=0,AMTag=0
        integer(TagIntType) :: WorkTag=0
        real(dp) :: CASRefEnergy,TotWeight,PartFac,amp,rat,r,GetHElement
        real(dp), dimension(lenof_sign) :: temp_sign
        real(dp) :: energytmp(nel), max_wt
        integer  :: tmp_det(nel), det_max, run
        type(ll_node), pointer :: TempNode
        character(len=*) , parameter :: this_routine='InitFCIMC_CAS'
#ifdef __CMPLX
        call stop_all(this_routine,"StartCAS currently does not work with complex walkers")
#endif
        if(tReadPops) call stop_all(this_routine,"StartCAS cannot work with with ReadPops")
        if(tStartSinglePart) call stop_all(this_routine,"StartCAS cannot work with StartSinglePart")
        if(tRestartHighPop) call stop_all(this_routine,"StartCAS cannot with with dynamically restarting calculations")


        write(iout,*) "Initialising walkers proportional to a CAS diagonalisation..."
        write(iout,'(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (" &
            ,(OccCASOrbs+VirtCASOrbs)/2,",",OccCASOrbs,")"
        DO I=NEl-OccCASorbs+1,NEl
            WRITE(iout,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(iout,G1(BRR(I))%SYM,.FALSE.)
            WRITE(iout,'(I4)',advance='no') G1(BRR(I))%Ml
            WRITE(iout,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
        ENDDO
        WRITE(iout,'(A)') " ================================================================================================="
        DO I=NEl+1,NEl+VirtCASOrbs
            WRITE(iout,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(iout,G1(BRR(I))%SYM,.FALSE.)
            WRITE(iout,'(I4)',advance='no') G1(BRR(I))%Ml
            WRITE(iout,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
        ENDDO

        CASSpinBasisSize=OccCASorbs+VirtCASorbs
        allocate(CASBrr(1:CASSpinBasisSize))
        allocate(CASDet(1:OccCasOrbs))
        do i=1,CASSpinBasisSize
            !Run through the cas space, and create an array which will map these orbtials to the 
            !orbitals they actually represent.
            CASBrr(i)=BRR(i+(NEl-OccCasorbs))
        enddo

        !Calculate symmetry of CAS determinants, and check that this will be the same as the reference determinant
        !for the rest of the FCIMC calculations.
!        do i=1,OccCASOrbs
!            CASDet(i)=CASBrr(i)
!        enddo
        elec=1
        do i=NEl-OccCasOrbs+1,NEl
            CASDet(elec)=ProjEDet(i)
            elec=elec+1
        enddo

        write(iout,*) "CAS Det is: "
        call write_det_len(iout,CASDet,OccCASOrbs,.true.)
        call GetSym(CASDet,OccCASOrbs,G1,nBasisMax,CASSym)
        write(iout,*) "Spatial symmetry of CAS determinants: ",CASSym%Sym%S
        write(iout,*) "Ms of CAS determinants: ",CASSym%Ms
        if(tFixLz) then
            write(iout,*) "Ml of CAS determinants: ",CASSym%Ml
        endif
        call neci_flush(iout)

        if(CASSym%Ml.ne.LzTot) call stop_all(this_routine,"Ml of CAS ref det does not match Ml of full reference det")
        if(CASSym%Ms.ne.0) call stop_all(this_routine,"CAS diagonalisation can only work with closed shell CAS spaces initially")
        if(CASSym%Sym%S.ne.HFSym%Sym%S) then
            call stop_all(this_routine,"Sym of CAS ref det does not match Sym of fulll reference det")
        endif

        !First, we need to generate all the excitations.
        call gndts(OccCASorbs,CASSpinBasisSize,CASBrr,nBasisMax,CASDetList,.true.,G1,tSpn,LMS,.true.,CASSym,nCASDet,iCASDet)

        if(nCASDet.eq.0) call stop_all(this_routine,"No CAS determinants found.")
        write(iout,*) "Number of symmetry allowed CAS determinants found to be: ",nCASDet
        Allocate(CASDetList(OccCASorbs,nCASDet),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating CASDetList")
        CASDetList(:,:)=0

        !Now fill up CASDetList...
        call gndts (OccCASorbs, CASSpinBasisSize, CASBrr, nBasisMax, &
                    CASDetList, .false., G1, tSpn, LMS, .true., CASSym, &
                    nCASDet, iCASDet)

        !We have a complication here. If we calculate the hamiltonian from these CAS determinants, then we are not
        !including the mean-field generated from the other occupied orbitals. We need to either 'freeze' the occupied
        !orbitals and modify the 1 & two electron integrals, or add the other electrons back into the list. We do the latter.
        allocate(CASFullDets(NEl,nCASDet),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating CASFullDets")
        CASFullDets(:,:)=0

        ! Get the first part of a determinant with the lowest energy, rather
        ! than lowest index number orbitals
        energytmp = ARR(ProjEDet, 2)
        tmp_det = ProjEDet
        call sort(energytmp, tmp_det)

        ! Construct the determinants resulting from the CAS expansion.
        do i=1,nCASDet
            CASFullDets(1:nel-OccCASorbs,i) = tmp_det(1:nel-OccCASOrbs)
            CASFullDets(nel-OccCASorbs+1:nel,i) = CASDetList(1:OccCASorbs, i)
            call sort(CASFullDets(:,i))
        enddo
        deallocate(CASDetList)

        write(iout,*) "First CAS determinant in list is: "
        call write_det(iout,CASFullDets(:,1),.true.)

        if(nCASDet.gt.1300) then
            !Do lanczos
            nEval=4
        else
            nEval = nCASDet
        endif
        write(iout,"(A,I4,A)") "Calculating lowest ",nEval," eigenstates of CAS Hamiltonian..."
        Allocate(Ck(nCASDet,nEval),stat=ierr)
        Ck=0.0_dp
        Allocate(W(nEval),stat=ierr)    !Eigenvalues
        W=0.0_dp
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating")
        
        write(iout,*) "Calculating hamiltonian..."
        allocate(nRow(nCASDet),stat=ierr)
        nRow=0
        ICMax=1
        tMC=.false.

        !HACK ALERT!! Need to fill up array in space of determinants, not HPHF functions.
        !Turn of tHPHFInts and turn back on when hamiltonian constructed.
        tHPHFInts=.false.

        CALL Detham(nCASDet,NEl,CASFullDets,Hamil,Lab,nRow,.true.,ICMax,GC,tMC)
        LenHamil=GC
        write(iout,*) "Allocating memory for hamiltonian: ",LenHamil*2
        Allocate(Hamil(LenHamil),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating Hamil")
        Hamil=0.0_dp
        Allocate(Lab(LenHamil),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating Lab")
        Lab=0
        call Detham(nCASDet,NEl,CASFullDets,Hamil,Lab,nRow,.false.,ICMax,GC,tMC)

        CASRefEnergy=GETHELEMENT(1,1,HAMIL,LAB,NROW,NCASDET)
        write(iout,*) "Energy of first CAS det is: ",CASRefEnergy

        !Turn back on HPHF integrals if needed.
        if(tHPHF) tHPHFInts=.true.

!        if(abs(CASRefEnergy-Hii).gt.1.0e-7_dp) then
!            call stop_all(this_routine,"CAS reference energy does not match reference energy of full space")
!        endif

        if(nCASDet.gt.1300) then
            !Lanczos
            NKRY1=NKRY+1
            NBLOCK=MIN(NEVAL,NBLK)
            LSCR=MAX(nCASDet*NEVAL,8*NBLOCK*NKRY)
            LISCR=6*NBLOCK*NKRY
            ALLOCATE(A_Arr(NEVAL,NEVAL),stat=ierr)
            CALL LogMemAlloc('A_Arr',NEVAL**2,8,this_routine,ATag,ierr)
            A_Arr=0.0_dp
            ALLOCATE(V(nCASDet*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('V',nCASDet*NBLOCK*NKRY1,8,this_routine,VTag,ierr)
            V=0.0_dp
            ALLOCATE(AM(NBLOCK*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('AM',NBLOCK*NBLOCK*NKRY1,8,this_routine,AMTag,ierr)
            AM=0.0_dp
            ALLOCATE(BM(NBLOCK*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('BM',NBLOCK*NBLOCK*NKRY,8,this_routine,BMTag,ierr)
            BM=0.0_dp
            ALLOCATE(T(3*NBLOCK*NKRY*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('T',3*NBLOCK*NKRY*NBLOCK*NKRY,8,this_routine,TTag,ierr)
            T=0.0_dp
            ALLOCATE(WT(NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('WT',NBLOCK*NKRY,8,this_routine,WTTag,ierr)
            WT=0.0_dp
            ALLOCATE(SCR(LScr),stat=ierr)
            CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag,ierr)
            SCR=0.0_dp
            ALLOCATE(ISCR(LIScr),stat=ierr)
            CALL LogMemAlloc('IScr',LIScr,4,this_routine,IScrTag,ierr)
            ISCR(1:LISCR)=0
            ALLOCATE(INDEX(NEVAL),stat=ierr)
            CALL LogMemAlloc('INDEX',NEVAL,4,this_routine,INDEXTag,ierr)
            INDEX(1:NEVAL)=0
            ALLOCATE(WH(nCASDet),stat=ierr)
            CALL LogMemAlloc('WH',nCASDet,8,this_routine,WHTag,ierr)
            WH=0.0_dp
            ALLOCATE(WORK2(3*nCASDet),stat=ierr)
            CALL LogMemAlloc('WORK2',3*nCASDet,8,this_routine,WORK2Tag,ierr)
            WORK2=0.0_dp
            ALLOCATE(V2(nCASDet,NEVAL),stat=ierr)
            CALL LogMemAlloc('V2',nCASDet*NEVAL,8,this_routine,V2Tag,ierr)
            V2=0.0_dp
            Allocate(CkN(nCASDet,nEval), stat=ierr)
            CkN=0.0_dp
    !C..Lanczos iterative diagonalising routine
            CALL NECI_FRSBLKH(nCASDet,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A_Arr,W,V,AM,BM,T,WT, &
         &  SCR,ISCR,INDEX,NCYCLE,B2L,.false.,.false.,.false.,.true.)
    !Multiply all eigenvalues by -1.
            CALL DSCAL(NEVAL,-1.0_dp,W,1)
            if(CK(1,1).lt.0.0_dp) then
                do i=1,nCASDet
                    CK(i,1)=-CK(i,1)
                enddo
            endif

            deallocate(CKN,A_Arr,V,BM,T,WT,SCR,WH,V2,iscr,index,AM)
            call logmemdealloc(this_routine,ATag)
            call logmemdealloc(this_routine,VTag)
            call logmemdealloc(this_routine,BMTag)
            call logmemdealloc(this_routine,TTag)
            call logmemdealloc(this_routine,WTTag)
            call logmemdealloc(this_routine,SCRTag)
            call logmemdealloc(this_routine,WHTag)
            call logmemdealloc(this_routine,V2Tag)
            call logmemdealloc(this_routine,iscrTag)
            call logmemdealloc(this_routine,indexTag)
            call logmemdealloc(this_routine,AMTag)
        else
            !complete diagonalisation
            allocate(Work(4*nCASDet),stat=ierr)
            call LogMemAlloc('Work',4*nCASDet,8,this_routine,WorkTag,ierr)
            allocate(Work2(3*nCASDet),stat=ierr)
            call logMemAlloc('Work2',3*nCASDet,8,this_routine,Work2Tag,ierr)
            nBlockStarts(1) = 1
            nBlockStarts(2) = nCASDet+1
            nBlocks = 1
            call HDIAG_neci(nCASDet,Hamil,Lab,nRow,CK,W,Work2,Work,nBlockStarts,nBlocks)
            deallocate(Work)
            call LogMemDealloc(this_routine,WorkTag)
        endif
        !Deallocate all the lanczos arrays now.
        deallocate(nrow,lab,work2)
        call logmemdealloc(this_routine,Work2Tag)


        write(iout,*) "Diagonalisation complete. Lowest energy CAS eigenvalues/corr E are: "
        do i=1,min(NEval,10)
            write(iout,*) i,W(i),W(i)-CASRefEnergy
        enddo

        TotWeight=0.0_dp
        nHPHFCAS=0
        max_wt = 0
        det_max = 0
        do i=1,nCASDet
            if(tHPHF) then
                !Only allow valid HPHF functions
                call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                if(IsAllowedHPHF(iLutnJ)) then
                    nHPHFCAS=nHPHFCAS+1
                    if(.not.TestClosedShellDet(iLutnJ)) then
                        !Open shell. Weight is sqrt(2) of det weight.
                        TotWeight=TotWeight+(abs(CK(i,1))*sqrt(2.0_dp))
                        !Return this new weight to the CK array, so that we do not need to do this a second time.
                        CK(i,1)=CK(i,1)*sqrt(2.0_dp)
                    else
                        !Closed Shell
                        TotWeight=TotWeight+abs(CK(i,1))
                    endif
                endif
            elseif(tMomInv) then
                !Only allow valid HPHF functions
                call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                if(IsAllowedMI(CASFullDets(:,i),iLutnJ)) then
                    nHPHFCAS=nHPHFCAS+1
                    if(.not.IsAllowedMI(CASFullDets(:,i),iLutnJ)) then
                        !Momentum-coupled. Weight is sqrt(2) of det weight.
                        TotWeight=TotWeight+(abs(CK(i,1))*sqrt(2.0_dp))
                        !Return this new weight to the CK array, so that we do not need to do this a second time.
                        CK(i,1)=CK(i,1)*sqrt(2.0_dp)
                    else
                        !Closed Shell
                        TotWeight=TotWeight+abs(CK(i,1))
                    endif
                endif
            else
                TotWeight=TotWeight+abs(CK(i,1))
            endif

            ! Find the maximum weighted determinant.
            if (abs(ck(i, 1)) > max_wt) then
                max_wt = abs(ck(i, 1))
                det_max = i
            end if
        enddo

        ! Output details
        write(iout,*) "Total weight of lowest eigenfunction: ", TotWeight
        write(iout,*) "Maximum weighted det: ", det_max, max_wt

         !If the reference det is not the maximum weighted det, suggest that
         !we change it!
        if (.not. all(CASFullDets(:, det_max) == ProjEDet)) then
            write(iout,*) 'The specified reference determinant is not the &
                       &maximum weighted determinant in the CAS expansion'
            write(iout,*) 'Use following det as reference:'
            call write_det(6, CASFullDets(:, det_max), .true.)
            call warning_neci(this_routine, "Poor reference chosen")
        end if

        if(tMomInv) write(iout,*) "Converting into momentum-coupled space. Total MI functions: ",nHPHFCAS
        if(tHPHF) write(iout,*) "Converting into HPHF space. Total HPHF CAS functions: ",nHPHFCAS

        if((InitialPart.eq.1).or.(InitialPart.ge.(InitWalkers*nNodes)-50)) then
            !Here, all the walkers will be assigned to the CAS wavefunction.
            !InitialPart = 1 by default
            write(iout,"(A)") "All walkers specified in input will be distributed according to the CAS wavefunction."
            write(iout,"(A)") "Shift will be allowed to vary from the beginning"
            write(iout,"(A,F20.9)") "Setting initial shift to equal CAS correlation energy",W(1)-CASRefEnergy
            DiagSft=W(1)-CASRefEnergy
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac=(real(InitWalkers,dp)* real(nNodes,dp))/TotWeight
        else
            !Here, not all walkers allowed will be initialised to the CAS wavefunction.
            write(iout,"(A,I15,A)") "Initialising ",int(InitialPart), " walkers according to the CAS distribution."
            write(iout,"(A,I15)") "Shift will remain fixed until the walker population reaches ",int(InitWalkers*nNodes)
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac=real(InitialPart,dp)/TotWeight
            tSinglePartPhase(:)=.true.
        endif

        !Now generate all excitations again, creating the required number of walkers on each one.
        DetIndex=1
        NoatHF(:)=0.0
        TotParts(:)=0.0
        do i=1,nCASDet
            if(tHPHF) then
                call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                if(.not.IsAllowedHPHF(iLutnJ)) cycle
            elseif(tMomInv) then
                call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                if(.not.IsAllowedMI(CASFullDets(:,i),iLutnJ)) cycle
            endif
            iNode=DetermineDetNode(CASFullDets(:,i),0)
            if(iProcIndex.eq.iNode) then
                !Number parts on this det = PartFac*Amplitude
                amp=CK(i,1)*PartFac
                
                if (tRealCoeffByExcitLevel) ExcitLevel=FindBitExcitLevel(iLutnJ, iLutRef, nEl)
                if (tAllRealCoeff .or. &
                    & (tRealCoeffByExcitLevel.and.(ExcitLevel.le.RealCoeffExcitThresh))) then
                    NoWalkers=amp
                else
                    NoWalkers=int(amp)
                    rat=amp-real(NoWalkers,dp)
                    r=genrand_real2_dSFMT()
                    if(abs(rat).gt.r) then
                        if(amp < 0.0_dp) then
                            NoWalkers=NoWalkers-1
                        else
                            NoWalkers=NoWalkers+1
                        endif
                    endif
                endif

                if(NoWalkers.ne.0.0) then
                    call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                    if(DetBitEQ(iLutnJ,iLutRef,NIfDBO)) then
                        !Check if this determinant is reference determinant, so we can count number on hf.
                        do run=1,inum_runs
                            NoatHF(run) = NoWalkers
                        enddo
                    endif
                    call encode_det(CurrentDets(:,DetIndex),iLutnJ)
                    call clear_all_flags(CurrentDets(:,DetIndex))
                    do run=1,inum_runs
                        temp_sign(run) = NoWalkers
                    enddo
                    call encode_sign(CurrentDets(:,DetIndex),temp_sign)
                    if(tTruncInitiator) then
                        !Set initiator flag if needed (always for HF)
                        call CalcParentFlag(DetIndex,1,iInit)
                    endif
                    if(.not.tRegenDiagHEls) then
                        if(tHPHF) then
                            HDiagTemp = hphf_diag_helement(CASFullDets(:,i),iLutnJ)
                        elseif(tMomInv) then
                            HDiagTemp = MI_diag_helement(CASFullDets(:,i),iLutnJ)
                        else
                            HDiagTemp = get_helement(CASFullDets(:,i),CASFullDets(:,i),0)
                        endif
                        CurrentH(1,DetIndex)=real(HDiagTemp,dp)-Hii
                    endif

                    if (tHashWalkerList) then
                        DetHash = FindWalkerHash(CASFullDets(:,i), nWalkerHashes)
                        TempNode => HashIndex(DetHash)
                        ! If the first element in the list has not been used.
                        if (TempNode%Ind == 0) then
                            TempNode%Ind = DetIndex
                        else
                            do while (associated(TempNode%Next))
                                TempNode => TempNode%Next
                            end do
                            allocate(TempNode%Next)
                            nullify(TempNode%Next%Next)
                            TempNode%Next%Ind = DetIndex
                        end if
                        nullify(TempNode)
                    end if

                    DetIndex=DetIndex+1
                    do run=1,inum_runs
                        TotParts(run)=TotParts(run)+abs(NoWalkers)
                    enddo
                endif
            endif   !End if desired node
        enddo

        TotWalkers=DetIndex-1   !This is the number of occupied determinants on each node
        TotWalkersOld=TotWalkers
        if(.not.tHashWalkerList) then
            call sort(CurrentDets(:,1:TotWalkers),CurrentH(:,1:TotWalkers))
        endif

        !Set local&global variables
        TotPartsOld=TotParts
        call mpisumall(TotParts,AllTotParts)
        call mpisumall(NoatHF,AllNoatHF)
        call mpisumall(TotWalkers,AllTotWalkers)
        OldAllNoatHF=AllNoatHF
        do run=1,inum_runs
            OldAllHFCyc(run)=AllNoatHF(run)
            OldAllAvWalkersCyc(run)=AllTotParts(run)
        enddo
        AllTotWalkersOld=AllTotWalkers
        AllTotPartsOld=AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld=0.0_dp

        deallocate(CK,W,Hamil,CASBrr,CASDet,CASFullDets)

    end subroutine InitFCIMC_CAS 

!Routine to initialise the particle distribution according to the MP1 wavefunction.
!This hopefully will help with close-lying excited states of the same sym.
    subroutine InitFCIMC_MP1()
        use MomInv, only: IsAllowedMI
        use Determinants, only: GetH0Element3,GetH0Element4
        use SymExcit3 , only : GenExcitations3
        use CalcData , only : InitialPart
        real(dp) :: TotMP1Weight,amp,MP2Energy,PartFac,H0tmp,rat,r,energy_contrib
        HElement_t :: hel,HDiagtemp
        integer :: iExcits, exflag, Ex(2,2), nJ(NEl), ic, DetIndex, iNode
        integer :: iInit, Slot, DetHash, ExcitLevel, run, i
        integer(n_int) :: iLutnJ(0:NIfTot)
        real(dp) :: NoWalkers, temp_sign(lenof_sign)
        logical :: tAllExcitsFound, tParity
        type(ll_node), pointer :: TempNode
        character(len=*), parameter :: this_routine="InitFCIMC_MP1"

#ifdef __CMPLX
        call stop_all(this_routine,"StartMP1 currently does not work with complex walkers")
#endif
        if(tReadPops) call stop_all(this_routine,"StartMP1 cannot work with with ReadPops")
        if(tStartSinglePart) call stop_all(this_routine,"StartMP1 cannot work with StartSinglePart")
        if(tRestartHighPop) call stop_all(this_routine,"StartMP1 cannot with with dynamically restarting calculations")

        write(iout,*) "Initialising walkers proportional to the MP1 amplitudes..."

        if(tHPHF) then
            if(.not.TestClosedShellDet(iLutHF)) then
                call stop_all(this_routine,"Cannot use HPHF with StartMP1 if your reference is open-shell")
            endif
        endif

        if(tUEG) then
            !Parallel N^2M implementation of MP2 for UEG
            call CalcUEGMP2()
        endif

        !First, calculate the total weight - TotMP1Weight
        mp2energy=0.0_dp
        TotMP1Weight=1.0_dp
        iExcits=0
        tAllExcitsFound=.false.
        if(tUEG) then
            exflag=2
        else
            exflag=3
        endif
        Ex(:,:)=0
        do while(.true.)
            call GenExcitations3(HFDet,iLutHF,nJ,exflag,Ex,tParity,tAllExcitsFound,.false.)
            if(tAllExcitsFound) exit !All excits found

            call EncodeBitDet(nJ,iLutnJ)
            if(tHPHF) then
                !Working in HPHF Space. Check whether determinant generated is an 'HPHF'
                if(.not.IsAllowedHPHF(iLutnJ)) cycle
            elseif(tMomInv) then
                !Working in MI space, Check whether determinant generated is allowed
                if(.not.IsAllowedMI(nJ,iLutnJ)) cycle
            endif
            iExcits = iExcits + 1

            call return_mp1_amp_and_mp2_energy(nJ,iLutnJ,Ex,tParity,amp,energy_contrib)
            TotMP1Weight=TotMP1Weight+abs(amp)
            MP2Energy=MP2Energy+energy_contrib
        enddo

        if((.not.tHPHF).and.(.not.tMomInv).and.(iExcits.ne.(nDoubles+nSingles))) then
            write(iout,*) nDoubles,nSingles,iExcits
            call stop_all(this_routine,"Not all excitations accounted for in StartMP1")
        endif

        write(iout,"(A,2G25.15)") "MP2 energy calculated: ",MP2Energy,MP2Energy+Hii

        if((InitialPart.eq.1).or.(InitialPart.ge.(InitWalkers*nNodes)-50)) then
            !Here, all the walkers will be assigned to the MP1 wavefunction.
            !InitialPart = 1 by default
            write(iout,"(A)") "All walkers specified in input will be distributed according to the MP1 wavefunction."
            write(iout,"(A)") "Shift will be allowed to vary from the beginning"
            write(iout,"(A)") "Setting initial shift to equal MP2 correlation energy"
            DiagSft=MP2Energy
            !PartFac is the number of walkers that should reside on the HF determinant
            !in an intermediate normalised MP1 wavefunction. 
            PartFac=(real(InitWalkers,dp)* real(nNodes,dp))/TotMP1Weight
        else
            !Here, not all walkers allowed will be initialised to the MP1 wavefunction.
            write(iout,"(A,G15.5,A)") "Initialising ",InitialPart, " walkers according to the MP1 distribution."
            write(iout,"(A,G15.5)") "Shift will remain fixed until the walker population reaches ",InitWalkers*nNodes
            !PartFac is the number of walkers that should reside on the HF determinant
            !in an intermediate normalised MP1 wavefunction. 
            PartFac=real(InitialPart,dp)/TotMP1Weight
            tSinglePartPhase(:)=.true.
        endif


        !Now generate all excitations again, creating the required number of walkers on each one.
        DetIndex=1
        TotParts=0.0
        tAllExcitsFound=.false.
        if(tUEG) then
            exflag=2
        else
            exflag=3
        endif
        Ex(:,:)=0
        do while(.true.)
            call GenExcitations3(HFDet,iLutHF,nJ,exflag,Ex,tParity,tAllExcitsFound,.false.)
            if(tAllExcitsFound) exit !All excits found
            call EncodeBitDet(nJ,iLutnJ)
            if(tHPHF) then
                !Working in HPHF Space. Check whether determinant generated is an 'HPHF'
                if(.not.IsAllowedHPHF(iLutnJ)) cycle
            elseif(tMomInv) then
                !Working in MI space, Check whether determinant generated is allowed
                if(.not.IsAllowedMI(nJ,iLutnJ)) cycle
            endif

            iNode=DetermineDetNode(nJ,0)
            if(iProcIndex.eq.iNode) then
                call return_mp1_amp_and_mp2_energy(nJ,iLutnJ,Ex,tParity,amp,energy_contrib)
                amp = amp*PartFac

                if (tRealCoeffByExcitLevel) ExcitLevel=FindBitExcitLevel(iLutnJ, iLutRef, nEl)
                if (tAllRealCoeff .or. &
                    & (tRealCoeffByExcitLevel.and.(ExcitLevel.le.RealCoeffExcitThresh))) then
                    NoWalkers=amp
                else
                    NoWalkers=int(amp)
                    rat=amp-real(NoWalkers,dp)
                    r=genrand_real2_dSFMT()
                    if(abs(rat).gt.r) then
                        if(amp < 0.0_dp) then
                            NoWalkers=NoWalkers-1
                        else
                            NoWalkers=NoWalkers+1
                        endif
                    end if
                end if
                
                if(NoWalkers.ne.0.0) then
                    call encode_det(CurrentDets(:,DetIndex),iLutnJ)
                    call clear_all_flags(CurrentDets(:,DetIndex))
                    do run=1, inum_runs
                        temp_sign(run) = NoWalkers
                    enddo
                    call encode_sign(CurrentDets(:,DetIndex),temp_sign)
                    if(tTruncInitiator) then
                        !Set initiator flag if needed (always for HF)
                        call CalcParentFlag(DetIndex,1,iInit)
                    endif
                    if(.not.tRegenDiagHEls) then
                        if(tHPHF) then
                            HDiagTemp = hphf_diag_helement(nJ,iLutnJ) 
                        elseif(tMomInv) then
                            HDiagTemp = MI_diag_helement(nJ,iLutnJ)
                        else
                            HDiagTemp = get_helement(nJ,nJ,0)
                        endif
                        CurrentH(1,DetIndex)=real(HDiagTemp,dp)-Hii
                    endif

                    if (tHashWalkerList) then
                        DetHash = FindWalkerHash(nJ, nWalkerHashes)
                        TempNode => HashIndex(DetHash)
                        ! If the first element in the list has not been used.
                        if (TempNode%Ind == 0) then
                            TempNode%Ind = DetIndex
                        else
                            do while (associated(TempNode%Next))
                                TempNode => TempNode%Next
                            end do
                            allocate(TempNode%Next)
                            nullify(TempNode%Next%Next)
                            TempNode%Next%Ind = DetIndex
                        end if
                        nullify(TempNode)
                    end if

                    DetIndex=DetIndex+1
                    do run=1,inum_runs
                        TotParts(run)=TotParts(run)+abs(NoWalkers)
                    enddo
                endif
            endif   !End if desired node

            
        enddo

        !Now for the walkers on the HF det
        if(iHFProc.eq.iProcIndex) then
            if (tAllRealCoeff .or. tRealCoeffByExcitLevel) then
                    NoWalkers=PartFac
            else
                NoWalkers=int(PartFac)
                rat=PartFac-real(NoWalkers,dp)
                if (rat < 0.0_dp) &
                    call stop_all(this_routine,"Should not have negative weight on HF")
                r=genrand_real2_dSFMT()
                if(abs(rat).gt.r) NoWalkers=NoWalkers+1
            endif
            if(NoWalkers.ne.0) then
                call encode_det(CurrentDets(:,DetIndex),iLutHF)
                call clear_all_flags(CurrentDets(:,DetIndex))
                do run=1,inum_runs
                    temp_sign(run) = NoWalkers
                enddo
                call encode_sign(CurrentDets(:,DetIndex),temp_sign)
                if(tTruncInitiator) then
                    !Set initiator flag (always for HF)
                    call set_flag(CurrentDets(:,DetIndex),flag_is_initiator(1))
                    call set_flag(CurrentDets(:,DetIndex),flag_is_initiator(2))
                endif
                if(.not.tRegenDiagHEls) CurrentH(1,DetIndex)=0.0_dp
                if (tHashWalkerList) then
                    ! Now add the Hartree-Fock determinant (not with index 1).
                    DetHash = FindWalkerHash(HFDet, nWalkerHashes)
                    TempNode => HashIndex(DetHash)
                    ! If the first element in the list has not been used.
                    if (TempNode%Ind == 0) then
                        TempNode%Ind = DetIndex
                    else
                        do while (associated(TempNode%Next))
                            TempNode => TempNode%Next
                        end do
                        allocate(TempNode%Next)
                        nullify(TempNode%Next%Next)
                        TempNode%Next%Ind = DetIndex
                    end if
                    nullify(TempNode)
                end if
                DetIndex=DetIndex+1
                do run=1,inum_runs
                    TotParts(run)=TotParts(run)+abs(NoWalkers)
                    NoatHF(run) = NoWalkers
                enddo
            else
                call stop_all(this_routine,"No walkers initialised on the HF det with StartMP1")
            endif
        else
            NoatHF(:)=0.0_dp
        endif
            
        TotWalkers=DetIndex-1   !This is the number of occupied determinants on each node
        TotWalkersOld=TotWalkers


        if(.not.tHashWalkerList) then
            call sort(CurrentDets(:,1:TotWalkers),CurrentH(:,1:TotWalkers))
        endif

        !Set local&global variables
        TotPartsOld=TotParts
        call mpisumall(TotParts,AllTotParts)
        call mpisumall(NoatHF,AllNoatHF)
        call mpisumall(TotWalkers,AllTotWalkers)
        OldAllNoatHF=AllNoatHF
        do run=1,inum_runs
            OldAllHFCyc(run)=AllNoatHF(run)
            OldAllAvWalkersCyc(run)=AllTotParts(run)
        enddo
        AllTotWalkersOld=AllTotWalkers
        AllTotPartsOld=AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld=0.0_dp

    end subroutine InitFCIMC_MP1

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
        if (.not.tContinueAfterMP2) then
            call stop_all("CalcUEGMP2","Dying after calculation of MP2 energy...")
        endif

    end subroutine CalcUEGMP2
            
    SUBROUTINE DeallocFCIMCMemPar()
        use nElRDMMod, only: DeallocateRDM
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocFciMCMemPar'
        type(ll_node), pointer :: Curr, Prev
        integer :: i, ierr

        if(tHashWalkerList .or. tSemiStochastic .or. (tTrialWavefunction .and. tTrialHash)) then
            deallocate(RandomHash2,stat=ierr)
            if(ierr.ne.0) call stop_all(this_routine,"Err deallocating")
        end if
        if (tHashWalkerList) then
            ! Deallocate the linked list
            do i = 1, nWalkerHashes
                Curr => HashIndex(i)%Next
                Prev => HashIndex(i)
                nullify(Prev%Next)
                do while (associated(Curr))
                    Prev => Curr
                    Curr => Curr%Next
                    deallocate(Prev)
                    if(ierr.ne.0) call stop_all(this_routine,"Err deallocating")
                end do
            end do
            deallocate(HashIndex,stat=ierr)
            if(ierr.ne.0) call stop_all(this_routine,"Err deallocating")
            nullify(Curr)
            nullify(Prev)

            deallocate(FreeSlot,stat=ierr)
            if(ierr.ne.0) call stop_all(this_routine,"Err deallocating")
        endif

        IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
            DEALLOCATE(Histogram)
            DEALLOCATE(AllHistogram)
            IF(tHistSpawn) THEN
                DEALLOCATE(InstHist)
                DEALLOCATE(InstAnnihil)
                DEALLOCATE(AvAnnihil)
            ENDIF
            IF(iProcIndex.eq.0) THEN
                IF(tHistSpawn) THEN
                    DEALLOCATE(AllInstHist)
                    DEALLOCATE(AllAvAnnihil)
                    DEALLOCATE(AllInstAnnihil)
                ENDIF
            ENDIF
        ELSEIF(tHistEnergies) THEN
            DEALLOCATE(HistogramEnergy)
            DEALLOCATE(AttemptHist)
            DEALLOCATE(SpawnHist)
            DEALLOCATE(SinglesHist)
            DEALLOCATE(DoublesHist)
            DEALLOCATE(DoublesAttemptHist)
            DEALLOCATE(SinglesAttemptHist)
            DEALLOCATE(SinglesHistOccOcc)
            DEALLOCATE(SinglesHistVirtOcc)
            DEALLOCATE(SinglesHistOccVirt)
            DEALLOCATE(SinglesHistVirtVirt)
            IF(iProcIndex.eq.Root) THEN
                DEALLOCATE(AllHistogramEnergy)
                DEALLOCATE(AllAttemptHist)
                DEALLOCATE(AllSpawnHist)
                DEALLOCATE(AllSinglesAttemptHist)
                DEALLOCATE(AllSinglesHist)
                DEALLOCATE(AllDoublesAttemptHist)
                DEALLOCATE(AllDoublesHist)
                DEALLOCATE(AllSinglesHistOccOcc)
                DEALLOCATE(AllSinglesHistVirtOcc)
                DEALLOCATE(AllSinglesHistOccVirt)
                DEALLOCATE(AllSinglesHistVirtVirt)
            ENDIF
        ENDIF
        IF(tHistHamil) THEN
            DEALLOCATE(HistHamil)
            DEALLOCATE(AvHistHamil)
            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(AllHistHamil)
                DEALLOCATE(AllAvHistHamil)
            ENDIF
        ENDIF
        if (tHistSpinDist) call clean_hist_spin_dist()
        if (tHistExcitToFrom) call clean_hist_excit_tofrom()
        DEALLOCATE(WalkVecDets)
        CALL LogMemDealloc(this_routine,WalkVecDetsTag)
        IF(.not.tRegenDiagHEls) THEN
            DEALLOCATE(WalkVecH)
            CALL LogMemDealloc(this_routine,WalkVecHTag)
        ENDIF
        DEALLOCATE(SpawnVec)
        CALL LogMemDealloc(this_routine,SpawnVecTag)
        DEALLOCATE(SpawnVec2)
        CALL LogMemDealloc(this_routine,SpawnVec2Tag)
        if(tRDMonFly.and.(.not.tExplicitAllRDM).and.(.not.tHF_Ref_Explicit)) then
            DEALLOCATE(TempSpawnedParts)
            CALL LogMemDealloc(this_routine,TempSpawnedPartsTag)
        endif
        DEALLOCATE(HFDet)
        CALL LogMemDealloc(this_routine,HFDetTag)
        DEALLOCATE(iLutHF)
        DEALLOCATE(iLutRef)
        DEALLOCATE(ProjEDet)
        DEALLOCATE(iLutHF_True)
        DEALLOCATE(HFDet_True)
        IF(ALLOCATED(HighestPopDet)) DEALLOCATE(HighestPopDet)
        IF(ALLOCATED(RandomHash)) DEALLOCATE(RandomHash)

        IF(ALLOCATED(SpinInvBrr)) THEN
            CALL LogMemDealloc(this_routine,SpinInvBRRTag)
            DEALLOCATE(SpinInvBRR)
        ENDIF
        IF(ALLOCATED(CoreMask)) THEN
            DEALLOCATE(CoreMask)
            DEALLOCATE(CASMask)
        ENDIF
        IF(tPrintOrbOcc) THEN
            DEALLOCATE(OrbOccs)
            CALL LogMemDeAlloc(this_routine,OrbOccsTag)
        ENDIF
        IF(tPrintDoubsUEG) THEN
            DEALLOCATE(DoubsUEG)
            DEALLOCATE(DoubsUEGLookup)
        ENDIF

        IF(tHistInitPops) THEN
            if (allocated(HistInitPops)) then
                deallocate (HistInitPops)
                call LogMemDeAlloc (this_routine, HistInitPopsTag)
            endif
            IF(iProcIndex.eq.0) THEN
                if (allocated(AllHistInitPops)) then
                    deallocate (AllHistInitPops)
                    call LogMemDeAlloc(this_routine,AllHistInitPopsTag)
                endif
            ENDIF
        ENDIF

        IF(tRDMonFly) CALL DeallocateRDM()
        if (allocated(refdetflip)) deallocate(refdetflip)
        if (allocated(ilutrefflip)) deallocate(ilutrefflip)

        ! Cleanup excitation generation storage
        call clean_excit_gen_store (fcimc_excit_gen_store)

        ! Cleanup linear combination projected energy
        call clean_linear_comb ()

        ! Cleanup storage for spin projection
        call clean_yama_store ()

        if (tSemiStochastic) call end_semistoch()

        if (tTrialWavefunction) call end_trial_wf()


!There seems to be some problems freeing the derived mpi type.
!        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!Free the mpi derived type that we have created for the hashes.
!            CALL MPI_Type_free(mpilongintegertype,error)
!            IF(error.ne.MPI_SUCCESS) THEN
!                CALL MPI_Error_string(error,message,length,temp)
!                IF(temp.ne.MPI_SUCCESS) THEN
!                    WRITE(iout,*) "REALLY SERIOUS PROBLEMS HERE!",temp
!                    CALL neci_flush(iout)
!                ENDIF
!                WRITE(iout,*) message(1:length)
!                WRITE(iout,*) "ERROR FOUND"
!                CALL neci_flush(iout)
!            ENDIF
!        ENDIF

    END SUBROUTINE DeallocFCIMCMemPar

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

! This routine will change the reference determinant to DetCurr. It will 
! also re-zero all the energy estimators, since they now correspond to
! projection onto a different determinant.
SUBROUTINE ChangeRefDet(DetCurr)
    use constants, only: inum_runs
    use FciMCParMod, only: DeallocFCIMCMemPar, SetupParameters, &
                           InitFCIMCCalcPar
    use DeterminantData, only: FDet
    use FciMCData, only: initiatorstats_unit, tDebug, iter, fcimcstats_unit, &
                         complexstats_unit, fcimcstats_unit2
    use CalcData, only: tTruncInitiator, tDelayTruncInit
    use SystemData, only: NEl
    use LoggingData, only: tLogComplexPops
    use Parallel_neci, only: iProcIndex, root
    use constants
    IMPLICIT NONE
    INTEGER :: DetCurr(NEl),i

    do i=1,NEl
        FDet(i)=DetCurr(i)
    enddo

    WRITE(iout,"(A)") "*** Changing the reference determinant ***"
    WRITE(iout,"(A)") "Switching reference and zeroing energy counters - restarting simulation"
!        
!Initialise variables for calculation on each node
    Iter=1
    CALL DeallocFCIMCMemPar()
    IF(iProcIndex.eq.Root) THEN
        CLOSE(fcimcstats_unit)
        if(inum_runs.eq.2) CLOSE(fcimcstats_unit2)
        IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(initiatorstats_unit)
        IF(tLogComplexPops) CLOSE(complexstats_unit)
    ENDIF
    IF(TDebug) CLOSE(11)
    CALL SetupParameters()
    CALL InitFCIMCCalcPar()

END SUBROUTINE ChangeRefDet
        

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
    
