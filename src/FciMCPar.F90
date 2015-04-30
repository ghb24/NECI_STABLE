#include "macros.h"
module FciMCParMod

    ! This module contains the main loop for FCIMC calculations, and the
    ! main per-iteration processing loop.
    use SystemData, only: nel, tUEG2, hist_spin_dist_iter
    use CalcData, only: tFTLM, tSpecLanc, tExactSpec, tDetermProj, tMaxBloom, &
                        tUseRealCoeffs, tWritePopsNorm, tExactDiagAllSym, &
                        AvMCExcits, pops_norm_unit, iExitWalkers, &
                        iFullSpaceIter, semistoch_shift_iter, &
                        tOrthogonaliseReplicas, orthogonalise_iter, &
                        tDetermHFSpawning, use_spawn_hash_table, &
                        semistoch_shift_iter, ss_space_in, s_global_start, &
                        tContTimeFCIMC, trial_shift_iter, tStartTrialLater, &
                        tTrialWavefunction, tSemiStochastic, ntrial_ex_calc
    use LoggingData, only: tJustBlocking, tCompareTrialAmps, tChangeVarsRDM, &
                           tWriteCoreEnd, tNoNewRDMContrib, tPrintPopsDefault,&
                           compare_amps_period, PopsFileTimer, &
                           write_end_core_size
    use spin_project, only: spin_proj_interval, disable_spin_proj_varyshift, &
                            spin_proj_iter_count, generate_excit_spin_proj, &
                            get_spawn_helement_spin_proj, iter_data_spin_proj,&
                            attempt_die_spin_proj
    use nElRDMMod, only: tCalc_RDMEnergy, FinaliseRDM, calc_energy_from_rdm, &
                         fill_explicitrdm_this_iter, &
                         fill_rdm_offdiag_deterministic, &
                         fill_hist_explicitrdm_this_iter
    use procedure_pointers, only: attempt_die_t, generate_excitation_t, &
                                  get_spawn_helement_t
    use semi_stoch_gen, only: write_most_pop_core_at_end, init_semi_stochastic
    use semi_stoch_procs, only: is_core_state, check_determ_flag, &
                                determ_projection, average_determ_vector
    use trial_wf_gen, only: update_compare_trial_file, &
                            update_compare_trial_file, init_trial_wf
    use hash, only: clear_hash_table
    use hist, only: write_zero_hist_excit_tofrom, write_clear_hist_spin_dist
    use bit_reps, only: set_flag, clr_flag, add_ilut_lists
    use exact_diag, only: perform_exact_diag_all_symmetry
    use spectral_lanczos, only: perform_spectral_lanczos
    use bit_rep_data, only: nOffFlag, flag_determ_parent
    use errors, only: standalone_errors, error_analysis
    use orthogonalise, only: orthogonalise_replicas
    use PopsFileMod, only: WriteToPopsFileParOneArr
    use AnnihilationMod, only: DirectAnnihilation
    use exact_spectrum, only: get_exact_spectrum
    use determ_proj, only: perform_determ_proj
    use cont_time, only: iterate_cont_time
    use global_det_data, only: det_diagH
    use RotateOrbsMod, only: RotateOrbs
    use NatOrbsMod, only: PrintOrbOccs
    use ftlm_neci, only: perform_ftlm
    use soft_exit, only: ChangeVars
    use fcimc_initialisation
    use fcimc_iter_utils
    use fcimc_helper
    use fcimc_output
    use FciMCData
    use constants

#ifdef MOLPRO
    use outputResult
#endif

    implicit none

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)

#ifdef MOLPRO
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
        integer :: tmp_int(lenof_sign), i, istart, iRDMSamplingIter
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

        ! Set the output energy value to zero, as it is not used by all
        ! of the calculations, and is therefore giving spurious uninitialised
        ! values for the testcode to pick up...
        Energyxw = 0.0_dp

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
        iRDMSamplingIter = 1    !For how many iterations have we accumulated the RDM

        SumSigns = 0.0_dp
        SumSpawns = 0.0_dp

        ! In we go - start the timer for scaling curve!
        start_time = neci_etime(tstart)

        do while (.true.)
!Main iteration loop...
            if(TestMCExit(Iter,iRDMSamplingIter)) exit

            IFDEBUG(FCIMCDebug, 2) write(iout,*) 'Iter', iter

            if(iProcIndex.eq.root) s_start=neci_etime(tstart)

            ! Is this an iteration where semi-stochastic is turned on?
            if (semistoch_shift_iter /= 0 .and. all(.not. tSinglePartPhase)) then
                if ((Iter - maxval(VaryShiftIter)) == semistoch_shift_iter + 1) then
                    tSemiStochastic = .true.
                    call init_semi_stochastic(ss_space_in)
                end if
            end if

            ! Is this an iteration where trial-wavefunction estimators are
            ! turned on?
            if (tStartTrialLater .and. all(.not. tSinglePartPhase)) then
                if ((Iter - maxval(VaryShiftIter)) == trial_shift_iter + 1) then
                    tTrialWavefunction = .true.

                    if (tOrthogonaliseReplicas .or. (tExcitedStateKP .and. .not. tPairedKPReplicas)) then
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs)
                    else if (tExcitedStateKP .and. tPairedKPReplicas) then
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs/2)
                    else
                        call init_trial_wf(trial_space_in, ntrial_ex_calc, 1)
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
                call set_timer(Stats_Comms_Time)
                call calculate_new_shift_wrapper (iter_data_fciqmc, TotParts, .false.)
                call halt_timer(Stats_Comms_Time)

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
            if(tFillingStochRDMonFly) iRDMSamplingIter = iRDMSamplingIter + 1 

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
            Energyxw = tot_trial_numerator(1)/tot_trial_denom(1)
        else
            Energyxw=(ProjectionE(1)+Hii)
        end if
        
        iroot=1
        CALL GetSym(ProjEDet(:,1),NEl,G1,NBasisMax,RefSym)
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
        integer :: j, p, error, proc_temp, i, HFPartInd,isym
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
        integer :: proc, pos, determ_index
        real(dp) :: r, sgn(lenof_sign), prob_extra_walker
        integer :: DetHash, FinalVal, clash, PartInd, k, y
        type(ll_node), pointer :: TempNode

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
        
        call rezero_iter_stats_each_iter (iter_data)

        ! The processor with the HF determinant on it will have to check 
        ! through each determinant until it's found. Once found, tHFFound is
        ! true, and it no longer needs to be checked.

        ! This is a bit of a hack based on the fact that we mean something 
        ! different by exFlag for CSFs than in normal determinential code.
        ! It would be nice to fix this properly
        if (tCSF) exFlag = 7

        IFDEBUGTHEN(FCIMCDebug,iout)
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
                if (IsUnoccDet(SignCurr)) then
                    if (tFillingStochRDMonFly) then
                        call set_av_sgn(j, AvSignCurr)
                        call set_iter_occ(j, IterRDMStartCurr)
                    endif
                    cycle
                end if
            end if

            ! The current diagonal matrix element is stored persistently.
            HDiagCurr = det_diagH(j)

            if (tTruncInitiator) &
                call CalcParentFlag (j, parent_flags, HDiagCurr)

            ! As the main list (which is storing a hash table) no longer needs
            ! to be contiguous, we need to skip sites that are empty.
            if(IsUnoccDet(SignCurr)) then
                !It has been removed from the hash table already
                !Add to the "freeslot" list
                iEndFreeSlot=iEndFreeSlot+1
                FreeSlot(iEndFreeSlot)=j
                cycle
            endif

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

            ! Sum in any energy contribution from the determinant, including 
            ! other parameters, such as excitlevel info.
            ! This is where the projected energy is calculated.
            call SumEContrib (DetCurr, WalkExcitLevel,SignCurr, CurrentDets(:,j), HDiagCurr, 1.0_dp, .false., j)

            ! If we're on the Hartree-Fock, and all singles and doubles are in
            ! the core space, then there will be no stochastic spawning from
            ! this determinant, so we can the rest of this loop.
            if (ss_space_in%tDoubles .and. walkExcitLevel_toHF == 0 .and. tDetermHFSpawning) then
                if (tFillingStochRDMonFly) then
                    call set_av_sgn(j, AvSignCurr)
                    call set_iter_occ(j, IterRDMStartCurr)
                endif
                cycle
            end if

            ! Loop over the 'type' of particle. 
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
            !                 --> OR double run
            !                 --> part_type == 1, 2; population sets 1 and 2, both real
            do part_type = 1, lenof_sign
            
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
                            call create_particle_with_hash_table (nJ, ilutnJ, child, &
                                                                  part_type, CurrentDets(:,j))
                        else
                            call create_particle (nJ, iLutnJ, child, part_type, & 
                                                  CurrentDets(:,j),SignCurr,p, &
                                                  RDMBiasFacCurr, WalkersToSpawn)
                        end if

                    endif ! (child /= 0). Child created

                enddo ! Cycling over mulitple particles on same determinant.

            enddo   ! Cycling over 'type' of particle on a given determinant.

            if (tSemiStochastic) then
                ! If we are performing a semi-stochastic simulation and this state is in the
                ! deterministic space, then the death step is performed deterministically later.
                if (.not. tCoreDet) then
                    call walker_death (iter_data, DetCurr, &
                                       CurrentDets(:,j), HDiagCurr, SignCurr, &
                                       AvSignCurr, IterRDMStartCurr, j, WalkExcitLevel)
                else
                    if (tFillingStochRDMonFly) then
                        call set_av_sgn(j, AvSignCurr)
                        call set_iter_occ(j, IterRDMStartCurr)
                    endif
                end if
            else
                call walker_death (iter_data, DetCurr, &
                                   CurrentDets(:,j), HDiagCurr, SignCurr, &
                                   AvSignCurr, IterRDMStartCurr, j, WalkExcitLevel)
            end if

        enddo ! Loop over determinants.
        IFDEBUGTHEN(FCIMCDebug,2) 
            write(iout,*) 'Finished loop over determinants'
            write(iout,*) "Holes in list: ",iEndFreeSlot
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

        ! With this algorithm, the determinants do not move, and therefore
        ! TotWalkersNew is simply equal to TotWalkers
        TotWalkersNew=int(TotWalkers,sizeof_int)

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

        ! This indicates the number of determinants in the list + the number
        ! of holes that have been introduced due to annihilation.
        TotWalkers=TotWalkersNew

        CALL halt_timer(Annihil_Time)
        IFDEBUG(FCIMCDebug,2) WRITE(iout,*) "Finished Annihilation step"
        
        ! If we are orthogonalising the replica wavefunctions, to generate
        ! excited states, then do that here.
        if (tOrthogonaliseReplicas .and. iter > orthogonalise_iter) &
            call orthogonalise_replicas(iter_data)

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

