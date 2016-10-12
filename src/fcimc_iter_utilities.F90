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

    implicit none

contains

    ! TODO: COMMENTING
    subroutine iter_diagnostics ()

        character(*), parameter :: this_routine = 'iter_diagnostics'
        character(*), parameter :: t_r = this_routine
        real(dp) :: mean_walkers
        integer :: part_type, run

        ! Update the total imaginary time passed
        TotImagTime = TotImagTime + StepsSft * Tau

        ! Set Iter time to equal the average time per iteration in the
        ! previous update cycle.
        IterTime = IterTime / real(StepsSft,sp)

        ! Calculate the acceptance ratio
        if (tContTimeFCIMC .and. .not. tContTimeFull) then
            AccRat = real(cont_spawn_success) / real(cont_spawn_attempts)
        else
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
            if (pop_change < pop_highest(run) .and. pop_highest(run) > pop_change_min) then

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

        type(fcimc_iter_data) :: iter_data
        real(dp), intent(in) :: tot_parts_new(lenof_sign)
        real(dp), intent(out) :: tot_parts_new_all(lenof_sign)

        ! Allow room to send up to 1000 elements.
        real(dp) :: send_arr(1000)
        ! Allow room to receive up to 1000 elements.
        real(dp) :: recv_arr(1000)
        ! Equivalent arrays for HElement_t variables.
        HElement_t(dp) :: send_arr_helem(100)
        HElement_t(dp) :: recv_arr_helem(100)
        ! Allow room for 100 different arrays to be communicated.
        integer :: sizes(100)
        integer :: low, upp, run

        integer(int64) :: TotWalkersTemp
        real(dp) :: bloom_sz_tmp(0:2)
        real(dp) :: RealAllHFCyc(max(lenof_sign,inum_runs))
        real(dp) :: all_norm_psi_squared(inum_runs), all_norm_semistoch_squared(inum_runs)
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
        sizes(23) = size(bloom_count)
        sizes(24) = size(NoAtHF)
        sizes(25) = size(SumWalkersCyc)
        sizes(26) = 1 ! nspawned (single int, not an array)

        if (sum(sizes(1:26)) > 1000) call stop_all(t_r, "No space left in arrays for communication of estimates. Please increase &
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
        low = upp + 1; upp = low + sizes(23) - 1; send_arr(low:upp) = bloom_count;
        low = upp + 1; upp = low + sizes(24) - 1; send_arr(low:upp) = NoAtHF;
        low = upp + 1; upp = low + sizes(25) - 1; send_arr(low:upp) = SumWalkersCyc;
        low = upp + 1; upp = low + sizes(26) - 1; send_arr(low:upp) = nspawned;

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

        low = upp + 1; upp = low + sizes(23) - 1; all_bloom_count = nint(recv_arr(low:upp));
        low = upp + 1; upp = low + sizes(24) - 1; AllNoAtHf = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(25) - 1; AllSumWalkersCyc = recv_arr(low:upp);
        low = upp + 1; upp = low + sizes(26) - 1; nspawned_tot = nint(recv_arr(low));

        ! Communicate HElement_t variables:

        low = 0; upp = 0;

        sizes(1) = size(ENumCyc)
        sizes(2) = size(SumENum)
        sizes(3) = size(ENumCycAbs)
        sizes(4) = size(cyc_proje_denominator)
        sizes(5) = size(sum_proje_denominator)

        if (sum(sizes(1:5)) > 100) call stop_all(t_r, "No space left in arrays for communication of estimates. Please &
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
        end if

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
        call MPIReduce(bloom_sizes(1:2), MPI_MAX, bloom_sz_tmp(1:2))
        bloom_sizes(1:2) = bloom_sz_tmp(1:2)

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
        end if

        if (tTrialWavefunction) then
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
    
    end subroutine collate_iter_data

    subroutine update_shift (iter_data)

        use CalcData, only: tInstGrowthRate
     
        type(fcimc_iter_data), intent(in) :: iter_data
        integer(int64) :: tot_walkers
        logical, dimension(inum_runs) :: tReZeroShift
        real(dp), dimension(inum_runs) :: AllGrowRateRe, AllGrowRateIm
        real(dp), dimension(inum_runs)  :: AllHFGrowRate
        real(dp), dimension(lenof_sign) :: denominator, all_denominator
        integer :: error, i, proc, pos, run, lb, ub
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
                        if(TargetGrowRate(run).ne.0.0_dp) then
                            write(iout,"(A)") "Setting target growth rate to 1."
                            TargetGrowRate=0.0_dp
                        endif

                        ! If enabled, jump the shift to the value preducted by the
                        ! projected energy!
                        if (tJumpShift) then
                            DiagSft(run) = real(proje_iter(run),dp)
                            defer_update(run) = .true.
                        end if
                    elseif (abs_sign(AllNoatHF(lb:ub)) < (MaxNoatHF - HFPopThresh)) then
                        write (iout, '(a,i13,a)') 'No at HF has fallen too low - reentering the &
                                     &single particle growth phase on iteration',iter + PreviousCycles,' - particle number &
                                     &may grow again.'
                        tSinglePartPhase(run) = .true.
                        tReZeroShift(run) = .true.
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
#ifdef __CMPLX
                        if(sum(AllTotParts(lb:ub)).gt.TargetGrowRateWalk(run)) then
#else
                        if(AllTotParts(run).gt.TargetGrowRateWalk(run)) then
#endif
                            !Only allow targetgrowrate to kick in once we have > TargetGrowRateWalk walkers.
                            DiagSft(run) = DiagSft(run) - (log(AllGrowRate(run)-TargetGrowRate(run)) * SftDamp) / &
                                                (Tau * StepsSft)
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

#ifdef __CMPLX
                    DiagSftRe(run) = DiagSftRe(run) - (log(AllGrowRateRe(run)-TargetGrowRate(run)) * SftDamp) / &
                                                (Tau * StepsSft)
                    DiagSftIm(run) = DiagSftIm(run) - (log(AllGrowRateIm(run)-TargetGrowRate(run)) * SftDamp) / &
                                                (Tau * StepsSft)
#endif

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

                 ! When using a linear combination, the denominator is summed
                 ! directly.
                 all_sum_proje_denominator(run) = ARR_RE_OR_CPLX(AllSumNoatHF,run)
                 all_cyc_proje_denominator(run) = AllHFCyc(run)

                 ! Calculate the projected energy.
#ifdef __CMPLX
                 if (any(AllSumNoatHF(lb:ub) /= 0.0)) then
#else
                 if ((AllSumNoatHF(run) /= 0.0)) then
#endif
                         ProjectionE(run) = (AllSumENum(run)) / (all_sum_proje_denominator(run)) &
                                          + proje_ref_energy_offsets(run)
                         proje_iter(run) = (AllENumCyc(run)) / (all_cyc_proje_denominator(run)) &
                                         + proje_ref_energy_offsets(run)
                        AbsProjE(run) = (AllENumCycAbs(run)) / (all_cyc_proje_denominator(run)) &
                                      + proje_ref_energy_offsets(run)
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
            projectionE_tot = sum(AllSumENum(1:inum_runs)) &
                            / sum(all_sum_proje_denominator(1:inum_runs))
            proje_iter_tot = sum(AllENumCyc(1:inum_runs)) &
                           / sum(all_cyc_proje_denominator(1:inum_runs))

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

        cont_spawn_attempts = 0
        cont_spawn_success = 0

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
