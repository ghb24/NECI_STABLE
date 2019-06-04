#include "macros.h"

module fcimc_iter_utils

    use SystemData, only: nel, tHPHF, tNoBrillouin, tRef_Not_HF, max_ex_level
    use CalcData, only: tSemiStochastic, tChangeProjEDet, tTrialWavefunction, &
                        tCheckHighestPopOnce, tRestartHighPop, StepsSft, tau, &
                        tTruncInitiator, tJumpShift, TargetGrowRate, &
                        tLetInitialPopDie, InitWalkers, tCheckHighestPop, &
                        HFPopThresh, DiagSft, tShiftOnHFPop, iRestartWalkNum, &
                        FracLargerDet, tKP_FCIQMC, MaxNoatHF, SftDamp, &
                        nShiftEquilSteps, TargetGrowRateWalk, tContTimeFCIMC, &
                        tContTimeFull, pop_change_min, tPositiveHFSign, &
                        qmc_trial_wf, nEquilSteps, t_hist_tau_search, &
                        t_hist_tau_search_option, corespaceWalkers, &
                        allCorespaceWalkers, tSpinProject, &
                        tFixedN0, tSkipRef, N0_Target, &
                        tTrialShift, tFixTrial, TrialTarget, tEN2

    use cont_time_rates, only: cont_spawn_success, cont_spawn_attempts

    use LoggingData, only: tFCIMCStats2, tPrintDataTables, tLogEXLEVELStats, &
                           t_spin_measurements

    use semi_stoch_procs, only: recalc_core_hamil_diag

    use bit_rep_data, only: NIfD, NIfTot, NIfDBO
    use hphf_integrals, only: hphf_diag_helement
    use Determinants, only: get_helement
    use LoggingData, only: tFCIMCStats2, t_calc_double_occ, t_calc_double_occ_av, tWriteUnocc, &
         AllInitsPerExLvl, initsPerExLvl
    use tau_search, only: update_tau
    use rdm_data, only: en_pert_main, InstRDMCorrectionFactor
    use Parallel_neci
    use fcimc_initialisation
    use fcimc_output
    use fcimc_helper
    use FciMCData
    use constants
    use util_mod
    use double_occ_mod, only: inst_double_occ, all_inst_double_occ, sum_double_occ, &
                              sum_norm_psi_squared, inst_spin_diff, all_inst_spin_diff, &
                              inst_spatial_doub_occ, all_inst_spatial_doub_occ, &
                              sum_double_occ_vec, sum_spin_diff, rezero_spin_diff, &
                              rezero_double_occ_stats

    use tau_search_hist, only: update_tau_hist

    implicit none

contains

    ! TODO: COMMENTING
    subroutine iter_diagnostics ()

        character(*), parameter :: this_routine = 'iter_diagnostics'
        character(*), parameter :: t_r = this_routine
        integer :: run, part_type

        ! Update the total imaginary time passed
        TotImagTime = TotImagTime + StepsSft * Tau

        ! Set Iter time to equal the average time per iteration in the
        ! previous update cycle.
        IterTime = IterTime / real(StepsSft,sp)

        ! Calculate the acceptance ratio
        if (tContTimeFCIMC .and. .not. tContTimeFull) then
           if(abs(real(cont_spawn_attempts)) > eps) then 
              AccRat = real(cont_spawn_success) / real(cont_spawn_attempts)
           else
              AccRat = 0.0_dp
           endif
        else
            ! in the real-time fciqmc keep track of both distinct RK steps
            ! -> also need to keep track of the SumWalkersCyc... todo..
           if(all(abs(SumWalkersCyc) > eps)) then
              AccRat = real(Acceptances, dp) / SumWalkersCyc
           else
              AccRat = 0.0_dp
           endif
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
            tRestart = .false.
            do run = 1, inum_runs
                if (sum(AllTotParts(min_part_type(run):max_part_type(run))) ==0 )  then
                    write(iout,"(A)") "All particles have died. Restarting."
                    tRestart=.true.
                    exit
                endif
            enddo
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
!Initialise variables for calculation on each node
            CALL DeallocFCIMCMemPar()
            IF(iProcIndex.eq.Root) THEN
                CLOSE(fcimcstats_unit)
                if (inum_runs.eq.2) CLOSE(fcimcstats_unit2)
                IF(tTruncInitiator) CLOSE(initiatorstats_unit)
                IF(tLogComplexPops) CLOSE(complexstats_unit)
                if (tLogEXLEVELStats) close(EXLEVELStats_unit)
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
            Iter=1
            if (iProcIndex==root .and. tLogEXLEVELStats) &
                  write(EXLEVELStats_unit,'("#")', advance='no')
            return
        endif

    end subroutine iter_diagnostics

    subroutine population_check ()

        use HPHFRandExcitMod, only: ReturnAlphaOpenDet

        integer(int32) :: pop_highest(inum_runs), proc_highest(inum_runs)
        real(dp) :: pop_change
        integer :: ierr, run
        integer(int32) :: int_tmp(2)
        logical :: tSwapped, allocate_temp_parts, changed_any
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
                !write(6,*) 1.5 * maxval(iHighestPop), TempSpawnedPartsSize
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
        ! [Werner Dobrautz 15.5.2017:]
        ! maybe this samll error here is the cause of the failed test_suite
        ! runs.. 
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
#ifdef __CMPLX
            if (tReplicaReferencesDiffer) then
                pop_change = FracLargerDet * abs_sign(AllNoAtHF(min_part_type(run):max_part_type(run)))
            else
                pop_change = FracLargerDet * abs_sign(AllNoAtHF(1:(lenof_sign/inum_runs)))
            endif
#else
            if (tReplicaReferencesDiffer) then
                pop_change = FracLargerDet * abs(AllNoAtHF(run))
            else
                pop_change = FracLargerDet * abs(AllNoAtHF(1))
            endif
#endif
!            write(iout,*) "***",AllNoAtHF,FracLargerDet,pop_change, pop_highest,proc_highest
            ! Do we need to do a change?
            ! is this a valid comparison?? we ware comparing a real(dp) pop_change 
            ! with a (now) 32 bit integer..
            if (pop_change < real(pop_highest(run),dp) .and. & 
                real(pop_highest(run),dp) > pop_change_min) then

                if (tChangeProjEDet) then

                    ! Write out info!
                    changed_any = .true.
                    root_print 'Highest weighted determinant on run', run, &
                         'not reference det: ', pop_highest, abs_sign(AllNoAtHF( &
                         min_part_type(run):max_part_type(run)))

                    !
                    ! Here we are changing the reference det on the fly.
                    ! --> else block for restarting simulation.
                    !

                    ! Communicate the change to all dets and print out.
                    ! [W.D. 15.5.2017:]
                    ! we are typecasting here too.. 
                    ! we are casting a 32 bit int to a 64 bit ... 
                    ! that could cause troubles! 
!                     call MPIBcast (HighestPopDet(0:NIfTot, run), NIfTot+1, &
!                                    int(proc_highest(run),n_int))
                    call MPIBcast (HighestPopDet(0:NIfTot, run), NIfTot+1, &
                                   int(proc_highest(run),sizeof_int))

                    call update_run_reference(HighestPopDet(:, run), run)

                    ! Reset averages
                    SumENum = 0.0_dp
                    sum_proje_denominator = 0.0_dp
                    cyc_proje_denominator = 0.0_dp
                    SumNoatHF = 0.0_dp
                    VaryShiftCycles = 0
                    SumDiagSft = 0.0_dp
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
                        iRestartWalkNum < sum(AllTotParts(1:2))) then
#else
                elseif (tRestartHighPop .and. &
                        iRestartWalkNum < AllTotParts(1)) then
#endif
                    !
                    ! Here we are restarting the simulation with a new
                    ! reference. See above block for doing it on the fly.
                    !
                    
                    ! Broadcast the changed det to all processors
!                     call MPIBcast (HighestPopDet(:,run), NIfTot+1, &
!                                    int(proc_highest(run),n_int))
                    call MPIBcast (HighestPopDet(:,run), NIfTot+1, &
                                   int(proc_highest(run),sizeof_int))

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

    subroutine communicate_estimates(iter_data, tot_parts_new, tot_parts_new_all)
      
        ! This routine sums all estimators and stats over all processes.

        ! We want this to be done in as few MPI calls as possible. Therefore, all
        ! quantities are first placed into one of two arrays. There is one array
        ! for HElement_t(dp) estimates, and another array (real(dp)) for all other
        ! kinds and types. A single MPISumAll is then performed for both of these
        ! combined arrays. After communication, the summed results are copied to
        ! the appropriate final arrays, with the type and kind corrected when
        ! necessary.

        ! There are also a few separate MPI calls which reduce using MPI_MAX and
        ! MPI_MIN at the end.

        ! -- *IMPORTANT FOR DEVELOPERS* ---------------------------------------
        ! To add a new quantity to be communicated, you must give it a new entry
        ! in send_arr or send_arr_helem array. It is hopefully clear how to do
        ! this by analogy. You should also update the indices in the appropriate
        ! stop_all, so that it can be checked if enough memory has been assigned.
      use adi_data, only: nCoherentDoubles, nIncoherentDets, nConnection, &
           AllCoherentDoubles, AllIncoherentDets, AllConnection
        type(fcimc_iter_data) :: iter_data
        real(dp), intent(in) :: tot_parts_new(lenof_sign)
        real(dp), intent(out) :: tot_parts_new_all(lenof_sign)
        ! RT_M_Merge: Added real-time statistics for the newer communication scheme
        integer, parameter :: real_arr_size = 1000
        integer, parameter :: hel_arr_size = 100
        integer, parameter :: NoArrs = 35
        integer, parameter :: size_arr_size = 100
        ! RT_M_Merge: Doubled all array sizes since there are now two
        ! copies of most of the variables (necessary?)
        
        ! Allow room to send up to 1000 (2000 for rt) elements.
        real(dp) :: send_arr(real_arr_size)
        ! Allow room to receive up to 1000 (2000 for rt) elements.
        real(dp) :: recv_arr(real_arr_size)
        ! Equivalent arrays for HElement_t variables.
        HElement_t(dp) :: send_arr_helem(hel_arr_size)
        HElement_t(dp) :: recv_arr_helem(hel_arr_size)
        ! Equivalent arrays for EXLEVELStats (of exactly required size).
        real(dp) :: send_arr_WNorm(3*(NEl+1)*inum_runs), &
                    recv_arr_WNorm(3*(NEl+1)*inum_runs)
        ! Allow room for 100 different arrays to be communicated.
        integer :: sizes(size_arr_size)
        integer :: low, upp, run

        integer(int64) :: TotWalkersTemp
        ! [W.D.12.12.2017]
        ! allow for triples now: 
        ! Todo: make that more flexible in the future! 
        real(dp) :: bloom_sz_tmp(0:3)
        real(dp) :: RealAllHFCyc(max(lenof_sign,inum_runs))
!         real(dp) :: all_norm_psi_squared(inum_runs)
        real(dp) :: all_norm_semistoch_squared(inum_runs)
        character(len=*), parameter :: t_r = 'communicate_estimates'

        ! Remove the holes in the main list when wanting the number of uniquely
        ! occupied determinants.
        TotWalkersTemp = TotWalkers - HolesInList

        sizes = 0

        ! low will represent the lower bound of an array slice.
        low = 0
        ! upp will represent the upper bound of an array slice.
        upp = 0

        sizes(1 ) = size(SpawnFromSing)
        sizes(2 ) = size(iter_data%update_growth)
        sizes(3 ) = size(NoBorn)
        sizes(4 ) = size(NoDied)
        sizes(5 ) = size(HFCyc)
        sizes(6 ) = size(NoAtDoubs)
        sizes(7 ) = size(Annihilated)
        if (tTruncInitiator) then
            sizes(8 ) = size(NoAddedInitiators)
            sizes(9 ) = size(NoInitDets)
            sizes(10) = size(NoNonInitDets)
            sizes(11) = size(NoExtraInitDoubs)
            sizes(12) = size(InitRemoved)
            sizes(13) = size(NoAborted)
            sizes(14) = size(NoRemoved)
            sizes(15) = size(NoNonInitWalk)
            sizes(16) = size(NoInitWalk)
        end if
        sizes(17) = 1 ! TotWalkersTemp (single int, not an array)
        sizes(18) = size(norm_psi_squared)
        sizes(19) = size(norm_semistoch_squared)
        sizes(20) = size(TotParts)
        sizes(21) = size(tot_parts_new)
        sizes(22) = size(SumNoAtHF)
        sizes(23) = size(bloom_count(0:max_ex_level))
        sizes(24) = size(NoAtHF)
        sizes(25) = size(SumWalkersCyc)
        sizes(26) = 1 ! nspawned (single int, not an array)

        sizes(27) = 1 ! inst_double_occ
        if(tTruncInitiator) sizes(28) = 1 ! doubleSpawns
        ! communicate the coherence numbers for SI
        sizes(29) = 1
        sizes(30) = 1
        ! Perturbation correction
        sizes(31) = 1
        ! communicate the instant spin diff.. although i am not sure if this 
        ! gets too big..
        if (t_spin_measurements) then 
            sizes(32) = nBasis/2
            sizes(33) = nBasis/2
        end if
        ! truncated weight
        sizes(34) = 1
        ! inits per ex lvl
        sizes(35) = size(initsPerExLvl)

        if (sum(sizes(1:NoArrs)) > real_arr_size) call stop_all(t_r, &
             "No space left in arrays for communication of estimates. Please increase &
             & the size of the send_arr and recv_arr arrays in the source code.")


        low = upp + 1; upp = low + sizes(1 ) - 1; send_arr(low:upp) = SpawnFromSing;
        low = upp + 1; upp = low + sizes(2 ) - 1; send_arr(low:upp) = iter_data%update_growth;
        low = upp + 1; upp = low + sizes(3 ) - 1; send_arr(low:upp) = NoBorn;
        low = upp + 1; upp = low + sizes(4 ) - 1; send_arr(low:upp) = NoDied;
        low = upp + 1; upp = low + sizes(5 ) - 1; send_arr(low:upp) = HFCyc;
        low = upp + 1; upp = low + sizes(6 ) - 1; send_arr(low:upp) = NoAtDoubs;
        low = upp + 1; upp = low + sizes(7 ) - 1; send_arr(low:upp) = Annihilated;
        if (tTruncInitiator) then
            low = upp + 1; upp = low + sizes(8 ) - 1; send_arr(low:upp) = NoAddedInitiators;
            low = upp + 1; upp = low + sizes(9 ) - 1; send_arr(low:upp) = NoInitDets;
            low = upp + 1; upp = low + sizes(10) - 1; send_arr(low:upp) = NoNonInitDets;
            low = upp + 1; upp = low + sizes(11) - 1; send_arr(low:upp) = NoExtraInitDoubs;
            low = upp + 1; upp = low + sizes(12) - 1; send_arr(low:upp) = InitRemoved;
            low = upp + 1; upp = low + sizes(13) - 1; send_arr(low:upp) = NoAborted;
            low = upp + 1; upp = low + sizes(14) - 1; send_arr(low:upp) = NoRemoved;
            low = upp + 1; upp = low + sizes(15) - 1; send_arr(low:upp) = NoNonInitWalk;
            low = upp + 1; upp = low + sizes(16) - 1; send_arr(low:upp) = NoInitWalk;
        end if
        low = upp + 1; upp = low + sizes(17) - 1; send_arr(low:upp) = TotWalkersTemp;
        low = upp + 1; upp = low + sizes(18) - 1; send_arr(low:upp) = norm_psi_squared;
        low = upp + 1; upp = low + sizes(19) - 1; send_arr(low:upp) = norm_semistoch_squared;
        low = upp + 1; upp = low + sizes(20) - 1; send_arr(low:upp) = TotParts;
        low = upp + 1; upp = low + sizes(21) - 1; send_arr(low:upp) = tot_parts_new;

        low = upp + 1; upp = low + sizes(22) - 1; send_arr(low:upp) = SumNoAtHf;
        low = upp + 1; upp = low + sizes(23) - 1; send_arr(low:upp) = bloom_count(0:max_ex_level);
        low = upp + 1; upp = low + sizes(24) - 1; send_arr(low:upp) = NoAtHF;
        low = upp + 1; upp = low + sizes(25) - 1; send_arr(low:upp) = SumWalkersCyc;
        low = upp + 1; upp = low + sizes(26) - 1; send_arr(low:upp) = nspawned;
        ! double occ change:
        low = upp + 1; upp = low + sizes(27) - 1; send_arr(low:upp) = inst_double_occ

        if(tTruncInitiator) then
           low = upp + 1; upp = low + sizes(28) -1; send_arr(low:upp) = doubleSpawns;
        endif
        low = upp + 1; upp = low + sizes(29) - 1; send_arr(low:upp) = nCoherentDoubles
        low = upp + 1; upp = low + sizes(30) - 1; send_arr(low:upp) = nIncoherentDets
        low = upp + 1; upp = low + sizes(31) - 1; send_arr(low:upp) = nConnection

        if (t_spin_measurements) then
            low = upp + 1; upp = low + sizes(32) -1; send_arr(low:upp) = inst_spin_diff
            low = upp + 1; upp = low + sizes(33) - 1; send_arr(low:upp) = inst_spatial_doub_occ
        end if
        ! truncated weight
        low = upp + 1; upp = low + sizes(34) - 1; send_arr(low:upp) = truncatedWeight;        
        ! initiators per excitation level
        low = upp + 1; upp = low + sizes(35) - 1; send_arr(low:upp) = initsPerExLvl;     


        ! Perform the communication.
        call MPISumAll (send_arr(1:upp), recv_arr(1:upp))

        ! Now we just need each result to be extracted to the correct array, with
        ! the correct type.

        low = 0; upp = 0

        low = upp + 1; upp = low + sizes(1 ) - 1; AllSpawnFromSing = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(2 ) - 1; iter_data%update_growth_tot = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(3 ) - 1; AllNoBorn = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(4 ) - 1; AllNoDied = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(5 ) - 1; RealAllHFCyc = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(6 ) - 1; AllNoAtDoubs = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(7 ) - 1; AllAnnihilated = recv_arr(low:upp);
        if (tTruncInitiator) then
            low = upp + 1; upp = low + sizes(8 ) - 1; AllNoAddedInitiators = nint(recv_arr(low:upp), int64);
            low = upp + 1; upp = low + sizes(9 ) - 1; AllNoInitDets = nint(recv_arr(low:upp), int64);
            low = upp + 1; upp = low + sizes(10) - 1; AllNoNonInitDets = nint(recv_arr(low:upp), int64);
            low = upp + 1; upp = low + sizes(11) - 1; AllNoExtraInitDoubs = nint(recv_arr(low:upp), int64);
            low = upp + 1; upp = low + sizes(12) - 1; AllInitRemoved = nint(recv_arr(low:upp), int64);
            low = upp + 1; upp = low + sizes(13) - 1; AllNoAborted = recv_arr(low:upp);
            low = upp + 1; upp = low + sizes(14) - 1; AllNoRemoved = recv_arr(low:upp);
            low = upp + 1; upp = low + sizes(15) - 1; AllNoNonInitWalk = recv_arr(low:upp);
            low = upp + 1; upp = low + sizes(16) - 1; AllNoInitWalk = recv_arr(low:upp);
        end if
        low = upp + 1; upp = low + sizes(17) - 1; AllTotWalkers = nint(recv_arr(low), int64);
        low = upp + 1; upp = low + sizes(18) - 1; all_norm_psi_squared = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(19) - 1; all_norm_semistoch_squared = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(20) - 1; AllTotParts = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(21) - 1; tot_parts_new_all = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(22) - 1; AllSumNoAtHF = recv_arr(low:upp);

        low = upp + 1; upp = low + sizes(23) - 1; all_bloom_count(0:max_ex_level) = nint(recv_arr(low:upp));
        low = upp + 1; upp = low + sizes(24) - 1; AllNoAtHf = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(25) - 1; AllSumWalkersCyc = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(26) - 1; nspawned_tot = nint(recv_arr(low));
        ! double occ: 
        low = upp + 1; upp = low + sizes(27) - 1; all_inst_double_occ = recv_arr(low);
        if(tTruncInitiator) then
           low = upp + 1; upp = low + sizes(28) - 1; allDoubleSpawns = nint(recv_arr(low));
           doubleSpawns = 0
        endif
        low = upp + 1; upp = low + sizes(29) - 1; AllCoherentDoubles = recv_arr(low);
        low = upp + 1; upp = low + sizes(30) - 1; AllIncoherentDets = recv_arr(low);
        low = upp + 1; upp = low + sizes(31) - 1; AllConnection = recv_arr(low);

        if (t_spin_measurements) then
            low = upp + 1; upp = low + sizes(32) - 1; all_inst_spin_diff = recv_arr(low:upp)
            low = upp + 1; upp = low + sizes(33) - 1; all_inst_spatial_doub_occ = recv_arr(low:upp)
        end if

        ! truncated weight
        low = upp + 1; upp = low + sizes(34) - 1; AllTruncatedWeight = recv_arr(low);
        ! initiators per excitation level
        low = upp + 1; upp = low + sizes(35) - 1; AllInitsPerExLvl = recv_arr(low:upp);

        ! Communicate HElement_t variables:

        low = 0; upp = 0;

        sizes(1) = size(ENumCyc)
        sizes(2) = size(SumENum)
        sizes(3) = size(ENumCycAbs)
        sizes(4) = size(cyc_proje_denominator)
        sizes(5) = size(sum_proje_denominator)
        if (tTrialWavefunction) then
            sizes(6) = size(trial_numerator)
            sizes(7) = size(trial_denom)
            sizes(8) = size(trial_num_inst)
            sizes(9) = size(trial_denom_inst)
            sizes(10) = size(init_trial_numerator)
            sizes(11) = size(init_trial_denom)
        end if
        if (tEN2) sizes(12) = 1
        sizes(13) = size(InitsEnumCyc)

        if (sum(sizes(1:11)) > 100) call stop_all(t_r, "No space left in arrays for communication of estimates. Please &
                                                        & increase the size of the send_arr_helem and recv_arr_helem &
                                                        & arrays in the source code.")

        low = upp + 1; upp = low + sizes(1) - 1; send_arr_helem(low:upp) = ENumCyc;
        low = upp + 1; upp = low + sizes(2) - 1; send_arr_helem(low:upp) = SumENum;
        low = upp + 1; upp = low + sizes(3) - 1; send_arr_helem(low:upp) = ENumCycAbs;
        low = upp + 1; upp = low + sizes(4) - 1; send_arr_helem(low:upp) = cyc_proje_denominator;
        low = upp + 1; upp = low + sizes(5) - 1; send_arr_helem(low:upp) = sum_proje_denominator;
        if (tTrialWavefunction) then
            low = upp + 1; upp = low + sizes(6) - 1; send_arr_helem(low:upp) = trial_numerator;
            low = upp + 1; upp = low + sizes(7) - 1; send_arr_helem(low:upp) = trial_denom;
            low = upp + 1; upp = low + sizes(8) - 1; send_arr_helem(low:upp) = trial_num_inst;
            low = upp + 1; upp = low + sizes(9) - 1; send_arr_helem(low:upp) = trial_denom_inst;
            low = upp + 1; upp = low + sizes(10) - 1; send_arr_helem(low:upp) = init_trial_numerator;
            low = upp + 1; upp = low + sizes(11) - 1; send_arr_helem(low:upp) = init_trial_denom;
        end if
        if (tEN2) then
           low = upp + 1; upp = low + sizes(12) - 1; send_arr_helem(low) = en_pert_main%ndets;
        endif
        low = upp + 1; upp = low + sizes(13) - 1; send_arr_helem(low:upp) = InitsENumCyc;

        call MPISumAll (send_arr_helem(1:upp), recv_arr_helem(1:upp))

        low = 0; upp = 0;

        low = upp + 1; upp = low + sizes(1) - 1; AllENumCyc = recv_arr_helem(low:upp);
        low = upp + 1; upp = low + sizes(2) - 1; AllSumENum = recv_arr_helem(low:upp);
        low = upp + 1; upp = low + sizes(3) - 1; AllENumCycAbs = recv_arr_helem(low:upp);
        low = upp + 1; upp = low + sizes(4) - 1; all_cyc_proje_denominator = recv_arr_helem(low:upp);
        low = upp + 1; upp = low + sizes(5) - 1; all_sum_proje_denominator = recv_arr_helem(low:upp);
        if (tTrialWavefunction) then
            low = upp + 1; upp = low + sizes(6) - 1; tot_trial_numerator = recv_arr_helem(low:upp);
            low = upp + 1; upp = low + sizes(7) - 1; tot_trial_denom = recv_arr_helem(low:upp);
            low = upp + 1; upp = low + sizes(8) - 1; tot_trial_num_inst = recv_arr_helem(low:upp);
            low = upp + 1; upp = low + sizes(9) - 1; tot_trial_denom_inst = recv_arr_helem(low:upp);
            low = upp + 1; upp = low + sizes(10) - 1; tot_init_trial_numerator = recv_arr_helem(low:upp);
            low = upp + 1; upp = low + sizes(11) - 1; tot_init_trial_denom = recv_arr_helem(low:upp);
        end if
        if (tEN2) then
           low = upp + 1; upp = low + sizes(12) - 1; en_pert_main%ndets_all = recv_arr_helem(low);
        endif
        low = upp + 1; upp = low + sizes(13) - 1; AllInitsENumCyc = recv_arr_helem(low:upp);

        ! Optionally communicate EXLEVEL_WNorm.
        if (tLogEXLEVELStats) then
            upp = size(EXLEVEL_WNorm)
            send_arr_WNorm(1:upp) = reshape(EXLEVEL_WNorm, (/upp/))
            call MPISumAll (send_arr_WNorm(1:upp), recv_arr_WNorm(1:upp))
            AllEXLEVEL_WNorm = reshape(recv_arr_WNorm(1:upp), &
                                       shape(AllEXLEVEL_WNorm))
            ! Apply square root for L2 norm.
            AllEXLEVEL_WNorm(2,:,:) = sqrt(AllEXLEVEL_WNorm(2,:,:))
        end if

        ! Do some processing of the received data.

        ! Convert real array into HElement_t one.
        do run = 1, inum_runs
            AllHFCyc(run) = ARR_RE_OR_CPLX(RealAllHFCyc, run)
        end do

#ifdef __CMPLX
        norm_psi = sqrt(sum(all_norm_psi_squared))
        norm_semistoch = sqrt(sum(all_norm_semistoch_squared))
#else
        norm_psi = sqrt(all_norm_psi_squared)
        norm_semistoch = sqrt(all_norm_semistoch_squared)
#endif

        ! These require a different type of reduce operation, so are communicated
        ! separately to the above communication.
        call MPIAllReduce(bloom_sizes(1:max_ex_level), MPI_MAX, bloom_sz_tmp(1:max_ex_level))
        bloom_sizes(1:max_ex_level) = bloom_sz_tmp(1:max_ex_level)

        ! Arrays for checking load balancing.
        call MPIReduce(TotWalkersTemp, MPI_MAX, MaxWalkersProc)
        call MPIReduce(TotWalkersTemp, MPI_MIN, MinWalkersProc)
        call MPIReduce(max_cyc_spawn, MPI_MAX, all_max_cyc_spawn)

    end subroutine communicate_estimates

    subroutine collate_iter_data(iter_data, replica_pairs)

        type(fcimc_iter_data) :: iter_data
        logical, intent(in) :: replica_pairs

        integer :: run
        logical :: ltmp
        character(len=*), parameter :: this_routine = 'collate_iter_data'

        ! We should update tau searching if it is enabled, or if it has been
        ! enabled, and now tau is outside the range acceptable for tau
        ! searching
        if (.not. tSearchTau) then
            call MPIAllLORLogical(tSearchTauDeath, ltmp)
           tSearchTauDeath = ltmp
        end if

        if ((tSearchTau .or. (tSearchTauOption .and. tSearchTauDeath)) .and. .not. tFillingStochRDMOnFly) then   
            call update_tau()

        ! [Werner Dobrautz 4.4.2017:]
        else if (((t_hist_tau_search .or. (t_hist_tau_search_option .and. tSearchTauDeath)) &
            .and. (.not. tFillingStochRDMonFly))) then
            call update_tau_hist()
        end if

        if (tTrialWavefunction) then
            if (.not. qmc_trial_wf) then
                ! Becuase tot_trial_numerator/tot_trial_denom is the energy
                ! relative to the the trial energy, add on this contribution to
                ! make it relative to the HF energy.
                if (ntrial_excits == 1) then
                    tot_trial_numerator = tot_trial_numerator + (tot_trial_denom*trial_energies(1))
                    tot_init_trial_numerator = tot_init_trial_numerator + (tot_init_trial_denom*&
                         trial_energies(1))
                else
                    if (replica_pairs) then
                        do run = 2, inum_runs, 2
                            tot_trial_numerator(run-1:run) = tot_trial_numerator(run-1:run) + &
                                tot_trial_denom(run-1:run)*trial_energies(run/2)
                            tot_init_trial_numerator(run-1:run) = tot_init_trial_numerator(run-1:run) + &
                                tot_init_trial_denom(run-1:run)*trial_energies(run/2)
                        end do
                    else
                        tot_trial_numerator = tot_trial_numerator + (tot_trial_denom*trial_energies)
                        tot_init_trial_numerator = tot_init_trial_numerator + (tot_init_trial_denom*trial_energies)
                    end if
                end if
            end if
        end if
        
        ! [W.D]
        ! quick fix for the double occupancy: 
        if (t_calc_double_occ_av) then 
            ! sum up the squared norm after shift has set in TODO
            ! and use the mean value if multiple runs are used
            ! still thinking about if i only want to calc it after 
            ! equilibration
            sum_norm_psi_squared = sum_norm_psi_squared + & 
                sum(all_norm_psi_squared)/real(inum_runs,dp)

            ! and also sum up the double occupancy: 
            sum_double_occ = sum_double_occ + all_inst_double_occ
            ! the averaging is also controlled by the t_calc_double_occ_av
            ! logical.. maybe change that in the future to be more clear
            if (t_spin_measurements) then 
                sum_double_occ_vec = sum_double_occ_vec + all_inst_spatial_doub_occ
                sum_spin_diff = sum_spin_diff + all_inst_spin_diff
            end if
        end if

#ifdef __DEBUG
        if(.not. tfirst_cycle) then
           ! realtime case is handled seperately with the check_update_growth function
           ! as each RK step has to be monitored separately

           ! Write this 'ASSERTROOT' out explicitly to avoid line lengths problems
           if ((iProcIndex == root) .and. .not. tSpinProject .and. &
                all(abs(iter_data%update_growth_tot-(AllTotParts-AllTotPartsOld)) > 1.0e-5)) then
              write(iout,*) "update_growth: ",iter_data%update_growth_tot
              write(iout,*) "AllTotParts: ",AllTotParts
              write(iout,*) "AllTotPartsOld: ", AllTotPartsOld
              call stop_all (this_routine, &
                   "Assertation failed: all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld)")
           endif
        end if
#endif
    
    end subroutine collate_iter_data

    subroutine update_shift (iter_data)

        use CalcData, only: tInstGrowthRate, tL2GrowRate
     
        type(fcimc_iter_data), intent(in) :: iter_data
        integer(int64) :: tot_walkers
        logical, dimension(inum_runs) :: tReZeroShift
        real(dp), dimension(inum_runs) :: AllGrowRateRe, AllGrowRateIm
        real(dp), dimension(inum_runs)  :: AllHFGrowRate
        integer :: i, run, lb, ub
        logical, dimension(inum_runs) :: defer_update
        logical :: start_varying_shift

        ! Normally we allow the shift to vary depending on the conditions
        ! tested. Sometimes we want to defer this to the next cycle...
        defer_update(:) = .false.

        ! collate_iter_data --> The values used are only valid on Root
        if (iProcIndex == Root) then

           if(tL2GrowRate) then 
              ! use the L2 norm to determine the growrate
              do run = 1, inum_runs
                 AllGrowRate(run) = norm_psi(run) / old_norm_psi(run)
              end do

           else if(tInstGrowthRate) then

              ! Calculate the growth rate simply using the two points at
              ! the beginning and the end of the update cycle. 
              do run = 1, inum_runs
                 lb = min_part_type(run)
                 ub = max_part_type(run)
                 AllGrowRate(run) = (sum(iter_data%update_growth_tot(lb:ub) &
                      + iter_data%tot_parts_old(lb:ub))) &
                      / real(sum(iter_data%tot_parts_old(lb:ub)), dp)
              enddo

           else

              ! Instead attempt to calculate the average growth over every
              ! iteration over the update cycle
              do run = 1, inum_runs
                 AllGrowRate(run) = AllSumWalkersCyc(run)/real(StepsSft,dp) &
                      /OldAllAvWalkersCyc(run)
              enddo

           end if
            ! For complex case, obtain both Re and Im parts
#ifdef __CMPLX
            do run = 1, inum_runs
                lb = min_part_type(run)
                ub = max_part_type(run)
                if (iter_data%tot_parts_old(lb) > 0) then
                    AllGrowRateRe(run) = (iter_data%update_growth_tot(lb) + &
                                     iter_data%tot_parts_old(lb)) / &
                                     iter_data%tot_parts_old(lb)
                end if
                if (iter_data%tot_parts_old(ub) > 0) then
                    AllGrowRateIm(run) = (iter_data%update_growth_tot(ub) + &
                                         iter_data%tot_parts_old(ub)) / &
                                         iter_data%tot_parts_old(ub)
                end if
            enddo
#endif

            ! Exit the single particle phase if the number of walkers exceeds
            ! the value in the input file. If particle no has fallen, re-enter
            ! it.
            tReZeroShift = .false.
            do run=1,inum_runs
                lb = min_part_type(run)
                ub = max_part_type(run)


                if (tTrialShift .and. .not. tFixTrial(run) .and. tTrialWavefunction .and. abs(tot_trial_denom(run))>=TrialTarget) then
                    !When reaching target overlap with trial wavefunction, set flag to keep it fixed.
                    tFixTrial(run) = .True.

                    write (iout, '(a,i13,a,i1)') 'Exiting the varaible shift phase on iteration: ' &
                                 ,iter + PreviousCycles, ' - overlap with trial wavefunction of the following run is now fixed: ', run
                end if


                if(tFixedN0)then
                    if (.not. tSkipRef(run) .and. abs(AllHFCyc(run))>=N0_Target) then
                        !When reaching target N0, set flag to keep the population of reference det fixed.
                        tSkipRef(run) = .True.

                        write (iout, '(a,i13,a,i1)') 'Exiting the fixed shift phase on iteration: ' &
                                     ,iter + PreviousCycles, ' - reference population of the following run is now fixed: ', run
                        !Set these parameters because other parts of the code depends on them
                        VaryShiftIter(run) = Iter
                        iBlockingIter(run) = Iter + PreviousCycles
                        tSinglePartPhase(run) = .false.
                    end if

                    if(tSkipRef(run))then
                        !Use the projected energy as the shift to fix the
                        !population of the reference det and thus reduce the
                        !fluctuations of the projected energy.

                        !ToDo: Make DiafSft complex
                        DiagSft(run) = (AllENumCyc(run)) / (AllHFCyc(run))+proje_ref_energy_offsets(run)

                        ! Update the shift averages
                        if ((iter - VaryShiftIter(run)) >= nShiftEquilSteps) then
                            if ((iter-VaryShiftIter(run)-nShiftEquilSteps) < StepsSft) &
                                write (iout, '(a,i14)') 'Beginning to average shift value on iteration: ',iter + PreviousCycles
                            VaryShiftCycles(run) = VaryShiftCycles(run) + 1
                            SumDiagSft(run) = SumDiagSft(run) + DiagSft(run)
                            AvDiagSft(run) = SumDiagSft(run) / real(VaryShiftCycles(run), dp)                            
                        endif
                    else
                        !Keep shift equal to input till target reference population is reached.
                        DiagSft(run) = InputDiagSft(run)
                    end if

                else if(tFixTrial(run))then
                        !Use the trial energy as the shift to fix the
                        !overlap with trial wavefunction and thus reduce the
                        !fluctuations of the trial energy.

                        !ToDo: Make DiafSft complex
                        DiagSft(run) = (tot_trial_numerator(run) / tot_trial_denom(run))-Hii

                        ! Update the shift averages
                        if ((iter - VaryShiftIter(run)) >= nShiftEquilSteps) then
                            if ((iter-VaryShiftIter(run)-nShiftEquilSteps) < StepsSft) &
                                write (iout, '(a,i14)') 'Beginning to average shift value on iteration: ',iter + PreviousCycles
                            VaryShiftCycles(run) = VaryShiftCycles(run) + 1
                            SumDiagSft(run) = SumDiagSft(run) + DiagSft(run)
                            AvDiagSft(run) = SumDiagSft(run) / real(VaryShiftCycles(run), dp)
                        endif

                else !not Fixed-N0 and not Trial-Shift
                    if (TSinglePartPhase(run)) then
                        tot_walkers = InitWalkers * int(nNodes,int64)

#ifdef __CMPLX
                        if ((sum(AllTotParts(lb:ub)) > tot_walkers) .or. &
                             (abs_sign(AllNoatHF(lb:ub)) > MaxNoatHF)) then
                           !                     WRITE(iout,*) "AllTotParts: ",AllTotParts(1),AllTotParts(2),tot_walkers
                           write (iout, '(a,i13,a)') 'Exiting the single particle growth phase on iteration: ',iter + PreviousCycles, &
                                ' - Shift can now change'
                           VaryShiftIter(run) = Iter
                           iBlockingIter(run) = Iter + PreviousCycles
                           tSinglePartPhase(run) = .false.
                           ! [W.D.15.5.2017:]
                           ! we should remove these equal 0 comparisons..
                           !                         if(TargetGrowRate(run).ne.0.0_dp) then
                           if(abs(TargetGrowRate(run)) > EPS) then
                              write(iout,"(A)") "Setting target growth rate to 1."
                              TargetGrowRate=0.0_dp
                           endif

                           ! If enabled, jump the shift to the value preducted by the
                           ! projected energy!
                           if (tJumpShift) then
                              if (tJumpShift .and. & 
                                   (.not. (isnan(real(proje_iter(run),dp))) .or. & 
                                   .not. (is_inf(real(proje_iter(run),dp))))) then
                                 DiagSft(run) = real(proje_iter(run),dp)
                                 defer_update(run) = .true.
                              end if
                           endif
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
                            ! [W.D. 15.5.2017]
                            ! change equal 0 comps
    !                         if(TargetGrowRate(run).ne.0.0_dp) then
                            if(abs(TargetGrowRate(run)) > EPS) then
                                write(iout,"(A)") "Setting target growth rate to 1."
                                TargetGrowRate(run)=0.0_dp
                            endif

                            ! If enabled, jump the shift to the value preducted by the
                            ! projected energy!
                            if (tJumpShift) then
                                DiagSft(run) = real(proje_iter(run),dp)
                                defer_update(run) = .true.
                            end if
                        endif
#endif
                    else ! .not.tSinglePartPhase(run)

#ifdef __CMPLX
                        if (abs_sign(AllNoatHF(lb:ub)) < MaxNoatHF-HFPopThresh) then
#else
                        if (abs(AllNoatHF(run)) < MaxNoatHF-HFPopThresh) then
#endif
                            write (iout, '(a,i13,a)') 'No at HF has fallen too low - reentering the &
                                         &single particle growth phase on iteration',iter + PreviousCycles,' - particle number &
                                         &may grow again.'
                            tSinglePartPhase(run) = .true.
                            tReZeroShift(run) = .true.
                        endif

                    endif ! tSinglePartPhase(run) or not

                    ! How should the shift change for the entire ensemble of walkers 
                    ! over all processors.
                    if (((.not. tSinglePartPhase(run)).or.(TargetGrowRate(run).ne.0.0_dp)) .and.&
                        .not. defer_update(run)) then

                        !In case we want to continue growing, TargetGrowRate > 0.0_dp
                        ! New shift value
    !                     if(TargetGrowRate(run).ne.0.0_dp) then
                        ! [W.D. 15.5.2017]
                        if(abs(TargetGrowRate(run)) > EPS) then
#ifdef __CMPLX
                            if(sum(AllTotParts(lb:ub)).gt.TargetGrowRateWalk(run)) then
#else
                            if(AllTotParts(run).gt.TargetGrowRateWalk(run)) then
#endif
                                !Only allow targetgrowrate to kick in once we have > TargetGrowRateWalk walkers.
                                DiagSft(run) = DiagSft(run) - (log(AllGrowRate(run)-TargetGrowRate(run)) * SftDamp) / &
                                                    (Tau * StepsSft)
                                ! Same for the info shifts for complex walkers
#ifdef __CMPLX
                        DiagSftRe(run) = DiagSftRe(run) - (log(AllGrowRateRe(run)-TargetGrowRate(run)) * SftDamp) / &
                                                    (Tau * StepsSft)
                        DiagSftIm(run) = DiagSftIm(run) - (log(AllGrowRateIm(run)-TargetGrowRate(run)) * SftDamp) / &
                                                    (Tau * StepsSft)
#endif
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
            end if !tFixedN0 or not
                ! only update the shift this way if possible
                if(abs_sign(AllNoatHF(lb:ub)) > EPS) then
#ifdef __CMPLX
                    ! Calculate the instantaneous 'shift' from the HF population
                    HFShift(run) = -1.0_dp / abs_sign(AllNoatHF(lb:ub)) * &
                                        (abs_sign(AllNoatHF(lb:ub)) - abs_sign(OldAllNoatHF(lb:ub)) / &
                                      (Tau * real(StepsSft, dp)))
                    InstShift(run) = -1.0_dp / sum(AllTotParts(lb:ub)) * &
                                ((sum(AllTotParts(lb:ub)) - sum(AllTotPartsOld(lb:ub))) / &
                                 (Tau * real(StepsSft, dp)))
#else
                    ! Calculate the instantaneous 'shift' from the HF population
                    HFShift(run) = -1.0_dp / abs(AllNoatHF(run)) * &
                                        (abs(AllNoatHF(run)) - abs(OldAllNoatHF(run)) / &
                                      (Tau * real(StepsSft, dp)))
                    InstShift(run) = -1.0_dp / AllTotParts(run) * &
                                ((AllTotParts(run) - AllTotPartsOld(run)) / &
                                 (Tau * real(StepsSft, dp)))
#endif
             endif

                 ! When using a linear combination, the denominator is summed
                 ! directly.
                 all_sum_proje_denominator(run) = ARR_RE_OR_CPLX(AllSumNoatHF,run)
                 all_cyc_proje_denominator(run) = AllHFCyc(run)

                 ! Calculate the projected energy.

!<<<<<<< HEAD
                ! did khaldoon move this part?
!                  if(tFixedN0)then
!                      !When reaching target N0, set flag to keep the population of reference det fixed.
!                      if(.not. tSkipRef(run) .and. AllHFCyc(run)>=N0_Target) tSkipRef(run) = .True.
!  
!                      if(tSkipRef(run))then
!                          !Use the projected energy as the shift to fix the
!                          !population of the reference det and thus reduce the
!                          !fluctuations of the projected energy.
!                          DiagSft(run) = (AllENumCyc(run)) / (AllHFCyc(run))+proje_ref_energy_offsets(run)
!                      else
!                          !Keep shift equal to input till target reference population is reached.
!                          DiagSft(run) = InputDiagSft(run)
!                      end if
!  
!                  end if
!>>>>>>>>>>>>> check that!

                if (abs(AllSumNoAtHF(run)) > EPS) then 
                    ProjectionE(run) = (AllSumENum(run)) / (all_sum_proje_denominator(run)) &
                         + proje_ref_energy_offsets(run)
                 endif
                 if (abs(AllHFCyc(run)) > EPS) then
                    proje_iter(run) = (AllENumCyc(run)) / (all_cyc_proje_denominator(run)) &
                         + proje_ref_energy_offsets(run)
                    AbsProjE(run) = (AllENumCycAbs(run)) / (all_cyc_proje_denominator(run)) &
                         + proje_ref_energy_offsets(run)
                    inits_proje_iter(run) = (AllInitsENumCyc(run)) / (all_cyc_proje_denominator(run)) &
                         + proje_ref_energy_offsets(run)
                 endif
                ! If we are re-zeroing the shift
                if (tReZeroShift(run)) then
                    DiagSft(run) = 0.0_dp
                    VaryShiftCycles(run) = 0
                    SumDiagSft(run) = 0.0_dp
                    AvDiagSft(run) = 0.0_dp
                endif
            enddo

            ! Get some totalled values
            if(abs(sum(all_sum_proje_denominator(1:inum_runs))) > EPS) then
               projectionE_tot = sum(AllSumENum(1:inum_runs)) &
                    / sum(all_sum_proje_denominator(1:inum_runs))
            endif
            if(abs(sum(all_cyc_proje_denominator(1:inum_runs))) > EPS) then
               proje_iter_tot = sum(AllENumCyc(1:inum_runs)) &
                    / sum(all_cyc_proje_denominator(1:inum_runs))
               inits_proje_iter_tot = sum(AllInitsENumCyc(1:inum_runs)) &
                    / sum(all_cyc_proje_denominator(1:inum_runs))
            endif

        endif ! iProcIndex == root

        ! Broadcast the shift from root to all the other processors
        call MPIBcast (tSinglePartPhase)
        call MPIBcast (VaryShiftIter)
        call MPIBcast (DiagSft)
        call MPIBcast (tSkipRef)
        call MPIBcast (tFixTrial)
        call MPIBcast (VaryShiftCycles)
        call MPIBcast (SumDiagSft)
        call MPIBcast (AvDiagSft)

        do run = 1, inum_runs
            if (.not. tSinglePartPhase(run)) then
                TargetGrowRate(run)=0.0_dp

                if (tPreCond) then
                    if (iter > 80) tSearchTau = .false.
                else
                    tSearchTau = .false.
                end if
            endif
        enddo
       
    end subroutine update_shift 

    subroutine rezero_iter_stats_update_cycle (iter_data, tot_parts_new_all)
        
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new_all
        
        ! Zero all of the variables which accumulate for each iteration.

        IterTime = 0.0_sp
        SumWalkersCyc(:)=0.0_dp
        Annihilated = 0.0_dp
        Acceptances = 0.0_dp
        NoBorn = 0.0_dp
        SpawnFromSing = 0.0_dp
        NoDied = 0.0_dp
        ENumCyc = 0.0_dp
        InitsENumCyc = 0.0_dp
        ENumCycAbs = 0.0_dp
        HFCyc = 0.0_dp
        cyc_proje_denominator=0.0_dp
        trial_numerator = 0.0_dp
        trial_denom = 0.0_dp

        ! also reset the real-time specific quantities: 
        ! and maybe have to call this routine twice to rezero also the 
        ! inputted iter_data for both RK steps..
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

        ! and the norm
        old_norm_psi = norm_psi

        ! Reset the counters
        iter_data%update_growth = 0.0_dp
        iter_data%update_iters = 0
        iter_data%tot_parts_old = tot_parts_new_all

        max_cyc_spawn = 0.0_dp

        cont_spawn_attempts = 0
        cont_spawn_success = 0

        tfirst_cycle = .false.
        if (t_calc_double_occ) then
            call rezero_double_occ_stats()
            if (t_spin_measurements) then
                call rezero_spin_diff()
            end if
        end if
        ! reset the truncated weight
        truncatedWeight = 0.0_dp

        ! reset the logged number of initiators
        initsPerExLvl = 0

    end subroutine rezero_iter_stats_update_cycle

    subroutine calculate_new_shift_wrapper (iter_data, tot_parts_new, replica_pairs)

        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign), intent(in) :: tot_parts_new
        real(dp), dimension(lenof_sign) :: tot_parts_new_all
        logical, intent(in) :: replica_pairs

        call communicate_estimates(iter_data, tot_parts_new, tot_parts_new_all)
        call collate_iter_data (iter_data, replica_pairs)
        call iter_diagnostics ()
        if(tRestart) return
        call population_check ()
        call update_shift (iter_data)
        if(tSemiStochastic) call getCoreSpaceWalkers()
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
        
!        write(6,*) '===================================='
!        write(6,*) 'Nborn', iter_data%nborn, NoBorn
!        write(6,*) 'Ndied', iter_data%ndied, NoDied
!        write(6,*) 'Nannihil', iter_data%nannihil, Annihilated
!        write(6,*) 'Nabrt', iter_data%naborted, NoAborted
!        write(6,*) 'Nremvd', iter_data%nremoved, NoRemoved
!        write(6,*) '===================================='

        iter_data%update_growth = iter_data%update_growth + iter_data%nborn &
                                - iter_data%ndied - iter_data%nannihil &
                                - iter_data%naborted - iter_data%nremoved
        iter_data%update_iters = iter_data%update_iters + 1

    end subroutine update_iter_data

    function get_occ_dets() result(nOccDets)
      implicit none
      integer :: nOccDets
      integer :: i
      real(dp) :: check_sign(lenof_sign)

      nOccDets = 0
      do i = 1, TotWalkers
         call extract_sign(CurrentDets(:,i),check_sign)
         if(.not. IsUnoccDet(check_sign)) nOccDets = nOccDets + 1
      enddo
      
    end function get_occ_dets

    subroutine getCoreSpaceWalkers
      use semi_stoch_procs, only: check_determ_flag

      implicit none
      integer :: i
      real(dp) :: sgn(lenof_sign)

      corespaceWalkers = 0.0_dp
      do i = 1, TotWalkers
         if(check_determ_flag(CurrentDets(:,i))) then
            call extract_sign(CurrentDets(:,i),sgn)
            ! Just sum up all walkers
            corespaceWalkers = corespaceWalkers + sum(abs(sgn))
         endif
      enddo
            
    end subroutine getCoreSpaceWalkers

    !Fix the overlap with trial wavefunction by enforcing the value of a random determinant of the trial space
    !As long as the shift equals the trial energy, this should still give the right dynamics.
    subroutine fix_trial_overlap(iter_data)
        use util_mod, only: binary_search_first_ge
        type(fcimc_iter_data), intent(inout) :: iter_data

        HElement_t(dp), dimension(inum_runs) :: new_trial_denom, new_tot_trial_denom
        real(dp), dimension(lenof_sign) :: trial_delta, SignCurr, newSignCurr
        integer :: j, rand, det_idx, proc_idx, run, part_type, lbnd, ubnd,err
        integer :: trial_count, trial_indices(tot_trial_space_size)
        real(dp) :: amps(tot_trial_space_size), total_amp, total_amps(nProcessors)
        logical :: tIsStateDeterm

#ifdef __CMPLX
        call stop_all("fix_trial_overlap", "Complex wavefunction is not supported yet!")
#else

        !Calculate the new overlap
        new_trial_denom = 0.0
        new_tot_trial_denom = 0.0

        trial_count = 0
        total_amp = 0.0
        do j = 1, int(TotWalkers,sizeof_int)
            call extract_sign (CurrentDets(:,j), SignCurr)
            if (.not. IsUnoccDet(SignCurr) .and. test_flag(CurrentDets(:,j), flag_trial)) then
                trial_count = trial_count + 1
                trial_indices(trial_count) = j 
                amps(trial_count) = abs(current_trial_amps(1,j))
                total_amp = total_amp + amps(trial_count)
                !Update the overlap
                if (ntrial_excits == 1) then
                    new_trial_denom = new_trial_denom + current_trial_amps(1,j)*SignCurr
                else if (tReplicaReferencesDiffer.and. tPairedReplicas) then
                    do run = 2, inum_runs, 2
                        new_trial_denom(run-1:run) = new_trial_denom(run-1:run) + current_trial_amps(run/2,j)*SignCurr(run-1:run)
                    end do
                else if (ntrial_excits == lenof_sign) then
                    new_trial_denom = new_trial_denom + current_trial_amps(:,j)*SignCurr
                end if
            end if
        end do

        !Collecte overlaps from call processors
        call MPIAllReduce(new_trial_denom,MPI_SUM,new_tot_trial_denom)

        !Choose a random processor propotioanl to the sum of amplitudes of its trial space
        call MPIGather(total_amp, total_amps, err)
        if(iProcIndex .eq. root) then
            !write(6,*) "total_amps: ", total_amps
            do j=2, nProcessors
                total_amps(j) = total_amps(j)+total_amps(j-1)
            end do
            proc_idx = binary_search_first_ge(total_amps, genrand_real2_dSFMT() * total_amps(nProcessors))-1
        end if
        call MPIBCast(proc_idx)


        !write(6,*) "proc_idx", proc_idx
        !write(6,*) "total_count: ", trial_count
        !write(6,*) "amps: ", amps(1:trial_count)
        !Enforcing an update of the random determinant of the random processor
        if(iProcIndex .eq. proc_idx) then
            !Choose a random determinant
            do j=2, trial_count 
                amps(j) = amps(j)+amps(j-1)
            end do
            det_idx = trial_indices(binary_search_first_ge(amps(1:trial_count), genrand_real2_dSFMT() * amps(trial_count)))
            do part_type = 1, lenof_sign
                run = part_type_to_run(part_type)
                if(tFixTrial(run))then
                    trial_delta(part_type) = (tot_trial_denom(run)-new_tot_trial_denom(run))/current_trial_amps(part_type,det_idx)
                else
                    trial_delta(part_type) = 0.0
                end if
            end do

            call extract_sign (CurrentDets(:,det_idx), SignCurr)
            newSignCurr = SignCurr+trial_delta
            call encode_sign (CurrentDets(:,det_idx), newSignCurr)

            !Correct statistics filled by CalcHashTableStats
            iter_data%ndied = iter_data%ndied + abs(SignCurr)
            iter_data%nborn = iter_data%nborn + abs(newSignCurr)
            TotParts = TotParts + abs(newSignCurr) - abs(SignCurr)

            tIsStateDeterm = .False.
            if (tSemiStochastic) tIsStateDeterm = test_flag(CurrentDets(:,det_idx), flag_deterministic)

            norm_psi_squared = norm_psi_squared + (newSignCurr)**2 - SignCurr**2
            if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + (newSignCurr)**2 - SignCurr**2

            if (tCheckHighestPop) then
                do run = 1, inum_runs
                    lbnd = min_part_type(run)
                    ubnd = max_part_type(run)
                    if (abs_sign(newSignCurr(lbnd:ubnd)) > iHighestPop(run)) then
                        iHighestPop(run) = int(abs_sign(newSignCurr(lbnd:ubnd)))
                        HighestPopDet(:,run)=CurrentDets(:,det_idx)
                    end if
                end do
            end if
            if (tFillingStochRDMonFly) then
                if (IsUnoccDet(newSignCurr) .and. (.not. tIsStateDeterm)) then
                    if (DetBitEQ(CurrentDets(:,det_idx), iLutHF_True, NIfDBO)) then
                        AvNoAtHF = 0.0_dp
                        IterRDM_HF = Iter + 1
                    end if
                end if
            end if

            if (DetBitEQ(CurrentDets(:,det_idx), iLutHF_True, NIfDBO)) then
                InstNoAtHF=newSignCurr
            end if
        end if
#endif
    end subroutine fix_trial_overlap

end module fcimc_iter_utils
