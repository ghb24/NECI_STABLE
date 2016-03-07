#include "macros.h"

module fcimc_iter_utils

    use SystemData, only: nel, tHPHF, tNoBrillouin, tRef_Not_HF
    use CalcData, only: tSemiStochastic, tChangeProjEDet, tTrialWavefunction, &
                        tCheckHighestPopOnce, tRestartHighPop, StepsSft, tau, &
                        tTruncInitiator, tJumpShift, TargetGrowRate, &
                        tLetInitialPopDie, InitWalkers, tCheckHighestPop, &
                        HFPopThresh, DiagSft, tShiftOnHFPop, iRestartWalkNum, &
                        FracLargerDet, tKP_FCIQMC, MaxNoatHF, SftDamp, &
                        nShiftEquilSteps, TargetGrowRateWalk, tContTimeFCIMC, &
                        tContTimeFull, pop_change_min, tPositiveHFSign, &
                        qmc_trial_wf
    use cont_time_rates, only: cont_spawn_success, cont_spawn_attempts
    use LoggingData, only: tFCIMCStats2, tPrintDataTables
    use semi_stoch_procs, only: recalc_core_hamil_diag
    use fcimc_helper, only: update_run_reference
    use bit_rep_data, only: NIfD, NIfTot, NIfDBO
    use hphf_integrals, only: hphf_diag_helement
    use global_det_data, only: set_det_diagH
    use Determinants, only: get_helement
    use LoggingData, only: tFCIMCStats2
    use tau_search, only: update_tau
    use Parallel_neci
    use fcimc_initialisation
    use fcimc_output
    use fcimc_helper
    use FciMCData
    use constants
    use util_mod
#ifdef __REALTIME
    use real_time_data, only: SpawnFromSing_1, AllSpawnFromSing_1, &
        NoBorn_1, NoDied_1, AllNoBorn_1, AllNoDied_1, NoAtDoubs_1, AllNoAtDoubs_1, &
        Annihilated_1, AllAnnihilated_1, AllNoAddedInitiators_1, AllNoInitDets_1, &
        AllNoNonInitDets_1, AllInitRemoved_1, bloom_count_1, bloom_sizes_1, &
        AllNoAborted_1, AllNoInitWalk_1, AllNoNonInitWalk_1, AllNoRemoved_1, &
        all_bloom_count_1, NoAddedInitiators_1, AccRat_1, SumWalkersCyc_1, &
        nspawned_1, nspawned_tot_1, second_spawn_iter_data, TotParts_1, &
        AllTotParts_1, AllTotPartsOld_1, TotWalkers_1, AllTotWalkers_1, &
        AllTotWalkersOld_1, AllSumWalkersCyc_1, OldAllAvWalkersCyc_1
#endif 

    implicit none

contains

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
        if (tContTimeFCIMC .and. .not. tContTimeFull) then
            AccRat = real(cont_spawn_success) / real(cont_spawn_attempts)
        else
            ! in the real-time fciqmc keep track of both distinct RK steps
            ! -> also need to keep track of the SumWalkersCyc... todo..
#ifdef __REALTIME
            AccRat_1 = real(Acceptances_1, dp) / SumWalkersCyc_1
#endif
            AccRat = real(Acceptances, dp) / SumWalkersCyc
        end if


#ifndef __CMPLX
        if (tPositiveHFSign) then
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
            ! a restart not wanted in the real-time fciqmc.. 
#ifdef __REALTIME 
            call stop_all(this_routine, &
                "a restart due to all died walkers not wanted in the real-time fciqmc!")
#endif
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

    end subroutine iter_diagnostics

    subroutine population_check ()

        use HPHFRandExcitMod, only: ReturnAlphaOpenDet

        integer :: pop_highest(inum_runs), proc_highest(inum_runs)
        real(dp) :: pop_change, old_Hii
        integer :: det(nel), i, error, ierr, run
        integer(int32) :: int_tmp(2)
        logical :: tSwapped, allocate_temp_parts, changed_any
        HElement_t(dp) :: h_tmp
        character(*), parameter :: this_routine = 'population_check'
        character(*), parameter :: t_r = this_routine

        ! If we aren't doing this, then bail out...
        if (.not. tCheckHighestPop) return

        ! If we are accumulating RDMs, then a temporary spawning array is
        ! required of <~ the size of the largest occupied det.
        !
        ! This memory holds walkers spawned from one determinant. This
        ! allows us to test if we are spawning onto the same Dj multiple
        ! times. If only using connections to the HF (tHF_Ref_Explicit)
        ! no stochastic RDM construction is done, and this is not
        ! necessary.
        if (tRDMOnFly .and. .not. tExplicitAllRDM) then

            ! Test if we need to allocate or re-allocate the temporary
            ! spawned parts array
            allocate_temp_parts = .false.
            if (.not. allocated(TempSpawnedParts)) then
                allocate_temp_parts = .true.
                TempSpawnedPartsSize = 1000
            end if
            if (1.5 * maxval(iHighestPop) > TempSpawnedPartsSize) then
                ! This testing routine is only called once every update
                ! cycle. The 1.5 gives us a buffer to cope with particle
                ! growth
                TempSpawnedPartsSize = maxval(iHighestPop) * 1.5
                allocate_temp_parts = .true.
                write(6,*) 1.5 * maxval(iHighestPop), TempSpawnedPartsSize
            end if

            ! If we need to allocate this array, then do so.
            if (allocate_temp_parts) then
                if (allocated(TempSpawnedParts)) then
                    deallocate(TempSpawnedParts)
                    log_dealloc(TempSpawnedPartsTag)
                end if
                allocate(TempSpawnedParts(0:NIfDBO, TempSpawnedPartsSize), &
                         stat=ierr)
                log_alloc(TempSpawnedParts,TempSpawnedPartsTag,ierr)
                TempSpawnedParts = 0
                write(6,"(' Allocating temporary array for walkers spawned &
                           &from a particular Di.')")
                write(6,"(a,f14.6,a)") " This requires ", &
                    real(((NIfDBO+1) * TempSpawnedPartsSize * size_n_int), dp)&
                        /1048576.0_dp, " Mb/Processor"
            end if

        end if ! Allocating memory for RDMs

        ! Obtain the determinant (and its processor) with the highest pop
        ! in each of the runs.
        ! n.b. the use of int(iHighestPop) obviously introduces a small amount
        !      of error here, by ignoring the fractional part...
        if (tReplicaReferencesDiffer) then

            do run = 1, inum_runs
                call MPIAllReduceDatatype (&
                    (/int(iHighestPop(run),int32), int(iProcIndex,int32)/), 1, &
                                           MPI_MAXLOC, MPI_2INTEGER, int_tmp)
                pop_highest(run) = int_tmp(1)
                proc_highest(run) = int_tmp(2)
            end do
        else

            call MPIAllReduceDatatype (&
                (/int(iHighestPop(1),int32), int(iProcIndex,int32)/), 1, &
                                       MPI_MAXLOC, MPI_2INTEGER, int_tmp)
            pop_highest = int_tmp(1)
            proc_highest = int_tmp(2)

        end if


        changed_any = .false.
        do run = 1, inum_runs

            ! If using the same reference for all, then we don't consider the
            ! populations seperately...
            if (run /= 1 .and. .not. tReplicaReferencesDiffer) &
                exit
            
            ! What are the change conditions?
            if((lenof_sign == 2) .and. (inum_runs == 1)) then 
                pop_change = FracLargerDet * abs_sign(AllNoAtHF)
            else if (lenof_sign /= inum_runs) then
                call stop_all(this_routine, "Complex not yet supported in multi-run mode")
            else if (tReplicaReferencesDiffer) then
                pop_change = FracLargerDet * abs(AllNoAtHF(run))
            else
                pop_change = FracLargerDet * abs(AllNoATHF(1))
            endif
!            write(iout,*) "***",AllNoAtHF,FracLargerDet,pop_change, pop_highest,proc_highest
            ! Do we need to do a change?
            if (pop_change < pop_highest(run) .and. pop_highest(run) > pop_change_min) then

                if (tChangeProjEDet) then

                    ! Write out info!
                    changed_any = .true.
                    root_print 'Highest weighted determinant on run', run, &
                               'not reference det: ', pop_highest, abs_sign(AllNoAtHF)

                    !
                    ! Here we are changing the reference det on the fly.
                    ! --> else block for restarting simulation.
                    !

                    ! Communicate the change to all dets and print out.
                    call MPIBcast (HighestPopDet(0:NIfTot, run), NIfTot+1, &
                                   proc_highest(run))

                    call update_run_reference(HighestPopDet(:, run), run)

                    ! Reset averages
                    SumENum = 0
                    sum_proje_denominator = 0
                    cyc_proje_denominator = 0
                    SumNoatHF = 0.0_dp
                    VaryShiftCycles = 0
                    SumDiagSft = 0
                    root_print 'Zeroing all energy estimators.'

                    !Since we have a new reference, we must block only from after this point
                    iBlockingIter = Iter + PreviousCycles

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
                    !
                    ! Here we are restarting the simulation with a new
                    ! reference. See above block for doing it on the fly.
                    !
                    
                    ! Broadcast the changed det to all processors
                    call MPIBcast (HighestPopDet(:,run), NIfTot+1, &
                                   proc_highest(run))

                    call update_run_reference(HighestPopDet(:, run), run)
                    
                    ! Only update the global reference energies if they
                    ! correspond to run 1 (which is used for those)
                    if (run == 1) then
                        call ChangeRefDet (ProjEDet(:, 1))
                    end if

                    ! Reset values introduced in soft_exit (CHANGEVARS)
                    if (tCHeckHighestPopOnce) then
                        tChangeProjEDet = .false.
                        tCheckHighestPop = .false.
                        tCheckHighestPopOnce = .false.
                    endif

                endif

            endif
        end do

    end subroutine population_check

    subroutine collate_iter_data (iter_data, tot_parts_new, tot_parts_new_all, replica_pairs)

        type(fcimc_iter_data) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new
        real(dp), dimension(lenof_sign), intent(out) :: tot_parts_new_all
        logical, intent(in) :: replica_pairs

        integer :: int_tmp(5+2*lenof_sign), proc, pos, i
        real(dp) :: sgn(lenof_sign)
        HElement_t(dp) :: helem_tmp(3*inum_runs)
        HElement_t(dp) :: real_tmp(2*inum_runs) !*lenof_sign
        integer(int64) :: int64_tmp(8),TotWalkersTemp
        character(len=*), parameter :: this_routine='collate_iter_data'
        real(dp), dimension(max(lenof_sign,inum_runs)) :: RealAllHFCyc
        real(dp), dimension(inum_runs) :: all_norm_psi_squared, all_norm_semistoch_squared
        real(dp) :: bloom_sz_tmp(0:2)
        logical :: ltmp
        integer :: run
    
        ! Communicate the integers needing summation

        call MPIReduce(SpawnFromSing, MPI_SUM, AllSpawnFromSing)
        call MPIReduce(iter_data%update_growth, MPI_SUM, iter_data%update_growth_tot)
        call MPIReduce(NoBorn, MPI_SUM, AllNoBorn)
        call MPIReduce(NoDied, MPI_SUM, AllNoDied)
        call MPIReduce(HFCyc, MPI_SUM, RealAllHFCyc)
        call MPIReduce(NoAtDoubs, MPI_SUM, AllNoAtDoubs)
        call MPIReduce(Annihilated, MPI_SUM, AllAnnihilated)
        
        ! also do the same for the stats of the intermediate RK step in the 
        ! real-time fciqmc
#ifdef __REALTIME
        call MPIReduce(SpawnFromSing_1, MPI_SUM, AllSpawnFromSing_1)
        call MPIReduce(NoBorn_1, MPI_SUM, AllNoBorn_1)
        call MPIReduce(NoDied_1, MPI_SUM, AllNoDied_1)
!         call MPIReduce(HFCyc, MPI_SUM, RealAllHFCyc)
        call MPIReduce(NoAtDoubs_1, MPI_SUM, AllNoAtDoubs_1)
        call MPIReduce(Annihilated_1, MPI_SUM, AllAnnihilated_1)
        call MPIReduce(iter_data_fciqmc%update_growth, MPI_SUM, &
            iter_data_fciqmc%update_growth_tot)
        ! NOTE: i kind of mix up where what is stored... in the second_spawn 
        ! the most necessary info of the actual second step is stored.. 
        ! and in the NoDied vars. above, but the info of the first step is 
        ! in the general iter_data_fciqmc, and in the specific NoDied_1 ..
        ! dont mix that up!
#endif

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
#ifdef __REALTIME 
            call MPISum((/NoAddedInitiators_1(1), NoInitDets_1(1), &
                NoNonInitDets_1(1), NoExtraInitdoubs(1), InitRemoved_1(1)/),&
                int64_tmp(1:5))

            AllNoAddedInitiators_1(1) = int64_tmp(1)
            AllNoInitDets_1(1) = int64_tmp(2)
            AllNoNonInitDets_1(1) = int64_tmp(3)
            AllInitRemoved_1(1) = int64_tmp(5)

            call MPIReduce(NoAborted_1, MPI_SUM, AllNoAborted_1)
            call MPIReduce(NoRemoved_1, MPI_SUM, AllNoRemoved_1)
            call MPIReduce(NoNonInitWalk_1, MPI_SUM, AllNoNonInitWalk_1)
            call MPIReduce(NoInitWalk_1, MPI_SUM, AllNoInitWalk_1)
#endif
        endif

        ! 64bit integers
        !Remove the holes in the main list when wanting the number of uniquely occupied determinants
        TotWalkersTemp=TotWalkers-HolesInList
        call MPIReduce(TotwalkersTemp, MPI_SUM, AllTotWalkers)
        call MPIReduce(norm_psi_squared,MPI_SUM,all_norm_psi_squared)
        call MPIReduce(norm_semistoch_squared,MPI_SUM,all_norm_semistoch_squared)
        call MPIReduce(Totparts,MPI_SUM,AllTotParts)
        call MPIReduce(tot_parts_new,MPI_SUM,tot_parts_new_all)

#ifdef __REALTIME 
        ! hm.. how do i remove the number of old holes from the 1st RK step..
        ! damn.. this gets messy.. 
        TotWalkersTemp = TotWalkers_1 - HolesInList
        call MPIReduce(TotWalkersTemp, MPI_SUM, AllTotWalkers_1)
#endif
        ! also keep track of the real-time fciqmc specific runge-kutta steps
        ! quantities 
#ifdef __REALTIME 
        call MPIReduce(TotParts_1, MPI_SUM, AllTotParts_1) 
        ! do i need that:
!         call MPIReduce(tot_parts_new_1, MPI_SUM, tot_parts_new_all_1)
#endif

#ifdef __CMPLX
        norm_psi = sqrt(sum(all_norm_psi_squared))
        norm_semistoch = sqrt(sum(all_norm_semistoch_squared))
#else
        norm_psi = sqrt(all_norm_psi_squared)
        norm_semistoch = sqrt(all_norm_semistoch_squared)
#endif
        
        call MPIReduce(SumNoatHF, MPI_SUM, AllSumNoAtHF)
        ! HElement_t(dp) values (Calculates the energy by summing all on HF and 
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

#ifdef __REALTIME
            call MPISum(bloom_count_1(1:2), all_bloom_count_1(1:2))
            call MPIReduce(bloom_sizes_1(1:2), MPI_MAX, bloom_sz_tmp(1:2))
            bloom_sizes_1(1:2) = bloom_sz_tmp(1:2)
#endif
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

        ! hm.. do i need the grow rate.. and how to calc. it correctly.. 
#ifdef __REALTIME 
        call MPISum(SumWalkersCyc_1, AllSumWalkersCyc_1)
#endif

        ! The total number of spawned determinants.
        call MPISum (nspawned, nspawned_tot)

        ! seperate again the 2 RK steps
#ifdef __REALTIME 
        call MPISum(nspawned_1, nspawned_tot_1)
#endif

        !        WRITE(iout,*) "***",iter_data%update_growth_tot,AllTotParts-AllTotPartsOld

        ! We should update tau searching if it is enabled, or if it has been
        ! enabled, and now tau is outside the range acceptable for tau
        ! searching
        if (.not. tSearchTau) then
            call MPIAllLORLogical(tSearchTauDeath, ltmp)
           tSearchTauDeath = ltmp
        end if
        if ((tSearchTau .or. (tSearchTauOption .and. tSearchTauDeath)) .and. &
                            .not. tFillingStochRDMOnFly) then   
            call update_tau()
        end if

        if (tTrialWavefunction) then
            call MPIAllReduce(trial_numerator, MPI_SUM, tot_trial_numerator)
            call MPIAllReduce(trial_denom, MPI_SUM, tot_trial_denom)

            if (.not. qmc_trial_wf) then
                ! Becuase tot_trial_numerator/tot_trial_denom is the energy
                ! relative to the the trial energy, add on this contribution to
                ! make it relative to the HF energy.
                if (ntrial_excits == 1) then
                    tot_trial_numerator = tot_trial_numerator + (tot_trial_denom*trial_energies(1))
                else
                    if (replica_pairs) then
                        do run = 2, inum_runs, 2
                            tot_trial_numerator(run-1:run) = tot_trial_numerator(run-1:run) + &
                                tot_trial_denom(run-1:run)*trial_energies(run/2)
                        end do
                    else
                        tot_trial_numerator = tot_trial_numerator + (tot_trial_denom*trial_energies)
                    end if
                end if
            end if
        end if
        
#ifdef __DEBUG
#ifdef __REALTIME
        ! change that here in the real-time to handle it correctly for the 2 
        ! distinct RK steps..
        if ((iProcIndex == root) .and. .not. tSpinProject .and. &
            ! iter_data is the first step -> so i neet info about the first totparts
            ! i "just" could change the input so that second_spawn_iter_data 
            ! is the iter_data input -> so the correct info is stored and 
            ! handled in the same way..
            all(abs(iter_data_fciqmc%update_growth_tot - &
                (AllTotParts_1 - AllTotPartsOld_1)) > 1.0e-5)) then
            write(iout,*) " wrong info in the FIRST Runge-Kutta step: "
            write(iout,*) "update_growth: ", iter_data_fciqmc%update_growth_tot
            write(iout,*) "AllTotParts: ", AllTotParts_1
            write(iout,*) "AllTotPartsOld: ", AllTotPartsOld_1
            call stop_all(this_routine, &
                "Assertation failed: all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld)")
        end if

        if ((iProcIndex == root) .and. .not. tSpinProject .and. &
            all(abs(iter_data%update_growth_tot - &
            (AllTotParts - AllTotPartsOld)) > 1.0e-5)) then
            write(iout,*) " wrong info in the SECOND Runge-Kutta step: "
            write(iout,*) "update_growth: ", iter_data%update_growth_tot
            write(iout,*) "AllTotParts: ", AllTotParts
            write(iout,*) "AllTotPartsOld: ", AllTotPartsOld
            call stop_all(this_routine, &
                "Assertation failed: all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld)")
        end if

#else

        ! Write this 'ASSERTROOT' out explicitly to avoid line lengths problems
        if ((iProcIndex == root) .and. .not. tSpinProject .and. &
         all(abs(iter_data%update_growth_tot-(AllTotParts-AllTotPartsOld)) > 1.0e-5)) then
            write(iout,*) "update_growth: ",iter_data%update_growth_tot
            write(iout,*) "AllTotParts: ",AllTotParts
            write(iout,*) "AllTotPartsOld: ", AllTotPartsOld
            call stop_all (this_routine, &
                "Assertation failed: all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld)")
        endif
#endif
#endif
    
    end subroutine collate_iter_data

    subroutine update_shift (iter_data)

        use CalcData, only: tInstGrowthRate, tShiftProjectGrowth
     
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

        ! collate_iter_data --> The values used are only valid on Root
        if (iProcIndex == Root) then

            if(tInstGrowthRate) then

                ! Calculate the growth rate simply using the two points at
                ! the beginning and the end of the update cycle. 
                if (lenof_sign == 2 .and. inum_runs == 1) then
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

            else if (tShiftProjectGrowth) then

                ! Extrapolate the expected number of walkers at the end of the
                ! _next_ update cycle for calculating the shift. i.e. use
                !
                ! log((N_t + (N_t - N_(t-1))) / N_t)
                if (lenof_sign == 2 .and. inum_runs == 1) then
                    !COMPLEX
                        AllGrowRate(run) = &
                            (2 * sum(AllSumWalkersCyc) - sum(OldAllAvWalkersCyc)) &
                            / sum(AllSumWalkersCyc)
                else
                    do run = 1, inum_runs
                        AllGrowRate(run) = &
                            (2 * AllSumWalkersCyc(run) - OldAllAvWalkersCyc(run)) &
                            / AllSumWalkersCyc(run)
                    end do
                endif

            else

                ! Instead attempt to calculate the average growth over every
                ! iteration over the update cycle
                if (lenof_sign == 2 .and. inum_runs == 1) then
                    !COMPLEX
                    AllGrowRate = (sum(AllSumWalkersCyc)/real(StepsSft,dp)) &
                                    /sum(OldAllAvWalkersCyc)
#ifdef __REALTIME
                    print *, "toto: am i here?"
                    AllGrowRate_1 = (sum(AllSumWalkersCyc_1)/real(StepsSft,dp))  & 
                                    / sum(OldAllAvWalkersCyc_1)
#endif
                else
                    do run=1,inum_runs
                        AllGrowRate(run) = (AllSumWalkersCyc(run)/real(StepsSft,dp)) &
                                        /OldAllAvWalkersCyc(run)
                    enddo
                endif

            end if

            ! For complex case, obtain both Re and Im parts
            if (lenof_sign == 2 .and. inum_runs == 1) then
                if (iter_data%tot_parts_old(1) > 0) then
                    AllGrowRateRe = (iter_data%update_growth_tot(1) + &
                                     iter_data%tot_parts_old(1)) / &
                                     iter_data%tot_parts_old(1)
                end if
                if (iter_data%tot_parts_old(lenof_sign) > 0) then
                    AllGrowRateIm = (iter_data%update_growth_tot(lenof_sign) + &
                                         iter_data%tot_parts_old(lenof_sign)) / &
                                         iter_data%tot_parts_old(lenof_sign)
                end if
            endif

            ! Exit the single particle phase if the number of walkers exceeds
            ! the value in the input file. If particle no has fallen, re-enter
            ! it.
            tReZeroShift = .false.
            do run=1,inum_runs
                if (TSinglePartPhase(run)) then
                    tot_walkers = InitWalkers * int(nNodes,int64)

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
                         ProjectionE = (AllSumENum) / (all_sum_proje_denominator) &
                                     + proje_ref_energy_offsets
                         proje_iter = (AllENumCyc) / (all_cyc_proje_denominator) &
                                    + proje_ref_energy_offsets
                        AbsProjE = (AllENumCycAbs) / (all_cyc_proje_denominator) &
                                 + proje_ref_energy_offsets
                    endif
                 else
                     if ((AllSumNoatHF(run) /= 0.0)) then
                         ProjectionE(run) = (AllSumENum(run)) / (all_sum_proje_denominator(run)) &
                                          + proje_ref_energy_offsets(run)
                         proje_iter(run) = (AllENumCyc(run)) / (all_cyc_proje_denominator(run)) &
                                         + proje_ref_energy_offsets(run)
                        AbsProjE(run) = (AllENumCycAbs(run)) / (all_cyc_proje_denominator(run)) &
                                      + proje_ref_energy_offsets(run)
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

            ! Get some totalled values
#ifdef __CMPLX
            projectionE_tot = ProjectionE(1)
            proje_iter_tot = proje_iter(1)
#else
            projectionE_tot = sum(AllSumENum(1:inum_runs)) &
                            / sum(all_sum_proje_denominator(1:inum_runs))
            proje_iter_tot = sum(AllENumCyc(1:inum_runs)) &
                           / sum(all_cyc_proje_denominator(1:inum_runs))
#endif

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

        ! also reset the real-time specific quantities: 
        ! and maybe have to call this routine twice to rezero also the 
        ! inputted iter_data for both RK steps..
#ifdef __REALTIME 
        SumWalkersCyc_1(:) = 0.0_dp
        Annihilated_1 = 0.0_dp
        Acceptances_1 = 0.0_dp
        NoBorn_1 = 0.0_dp
        SpawnFromSing_1 = 0.0_dp
        NoDied = 0.0_dp
#endif

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

#ifdef __REALTIME
        OldAllAvWalkersCyc_1 = AllSumWalkersCyc_1 / real(StepsSft,dp)
#endif

        ! Also the cumulative global variables
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        AllNoAbortedOld = AllNoAborted

#ifdef __REALTIME 
        AllTotPartsOld_1 = AllTotParts_1
        AllTotWalkersOld_1 = AllTotWalkers_1
        ! do i need old det numner and aborted number? 
#endif

        ! Reset the counters
        iter_data%update_growth = 0.0_dp
        iter_data%update_iters = 0
        iter_data%tot_parts_old = tot_parts_new_all

        max_cyc_spawn = 0

        cont_spawn_attempts = 0
        cont_spawn_success = 0

    end subroutine rezero_iter_stats_update_cycle

    subroutine calculate_new_shift_wrapper (iter_data, tot_parts_new, replica_pairs)

        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new
        real(dp), dimension(lenof_sign) :: tot_parts_new_all
        logical, intent(in) :: replica_pairs

        call collate_iter_data (iter_data, tot_parts_new, tot_parts_new_all, replica_pairs)
        call iter_diagnostics ()
        if(tRestart) return
        call population_check ()
        call update_shift (iter_data)
        if (tPrintDataTables) then
            if (tFCIMCStats2) then
                call write_fcimcstats2(iter_data_fciqmc)
            else
                call WriteFCIMCStats ()
            end if
        end if
        
        call rezero_iter_stats_update_cycle (iter_data, tot_parts_new_all)

    end subroutine calculate_new_shift_wrapper

    subroutine update_iter_data(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data

        iter_data%update_growth = iter_data%update_growth + iter_data%nborn &
                                - iter_data%ndied - iter_data%nannihil &
                                - iter_data%naborted - iter_data%nremoved
        iter_data%update_iters = iter_data%update_iters + 1

    end subroutine update_iter_data

end module fcimc_iter_utils
