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
    use fcimc_initialisation
    use fcimc_iter_utils
    use fcimc_helper
    use fcimc_output
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

        

   
