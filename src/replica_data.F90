#include "macros.h"

module replica_data

    use constants
    use FciMCData
    use CalcData
    use util_mod
    use kp_fciqmc_data_mod
    use real_time_data
    implicit none

contains

    subroutine init_replica_arrays()

        ! We run the specified number of replicas in parallel.
        !
        ! This is initialisation of _global_ data that depends on the number
        ! of 

        character(*), parameter :: this_routine = 'init_replica_arrays'
        integer :: ierr

        if (inum_runs > inum_runs_max .or. lenof_sign > lenof_sign_max) then
            write(6,*) "System has too many replicas"
            write(6,*) "This breaks the use of initiator flags, and &
                       &potentially other things."
            call stop_all(this_routine, "Too many replicas requested")
        end if

        ! Global simulation properties, which apply per particle stream (i.e.
        ! depend on lenof_sign). There will be tow of these per simulation if
        ! complex walkers are used.
        allocate(InstNoatHF(lenof_sign), &
                 SumNoatHF(lenof_sign), AllSumNoatHF(lenof_sign), &
                 NoatHF(lenof_sign), AllNoAtHF(lenof_sign), &
                                     OldAllNoatHF(lenof_sign), &


                 TotParts(lenof_sign), AllTotParts(lenof_sign), &
                 TotPartsOld(lenof_sign), AllTotPartsOld(lenof_sign), &
                 TotPartsStorage(lenof_sign), TotPartsLastAlpha(lenof_sign), &
                 ! n.b. AllHFCyc is in inum_runs, with different type
                 HFCyc(lenof_sign), &
                 proje_denominator_cyc(lenof_sign), &
                 proje_denominator_sum(lenof_sign), &
                 InitialPartVec(lenof_sign), &

                 NoAborted(lenof_sign), AllNoAborted(lenof_sign), &
                                        AllNoAbortedOld(lenof_sign), &
                 NoRemoved(lenof_sign), AllNoRemoved(lenof_sign), &
                                        AllNoRemovedOld(lenof_sign), &

                 ! Initiator related
                 NoAddedInitiators(lenof_sign), &
                 NoInitDets(lenof_sign), NoNonInitDets(lenof_sign), &
                 NoInitWalk(lenof_sign), NoNonInitWalk(lenof_sign), &
                 NoExtraInitDoubs(lenof_sign), InitRemoved(lenof_sign), &

                 ! And initiator accumulator related
                 AllNoAddedInitiators(lenof_sign), &
                 AllNoInitDets(lenof_sign), AllNoNonInitDets(lenof_sign), &
                 AllNoInitWalk(lenof_sign), AllNoNonInitWalk(lenof_sign), &
                 AllNoExtraInitDoubs(lenof_sign), AllInitRemoved(lenof_sign), &
                 AllGrowRateAbort(lenof_sign), &

                 stat=ierr)


        ! Global simulation properties that apply per simulation, i.e. will
        ! not be duplicated if complex walkers are used. --> inum_runs
        allocate(&
                 ! Track the dynamics
                 NoBorn(inum_runs), AllNoBorn(inum_runs), &
                 NoDied(inum_runs), AllNoDied(inum_runs), &
                 Annihilated(inum_runs), AllAnnihilated(inum_runs), &
                 Acceptances(inum_runs), &
                 SpawnFromSing(inum_runs), AllSpawnFromSing(inum_runs), &
                 iRefProc(inum_runs), proje_ref_energy_offsets(inum_runs), &
                 iHighestPop(inum_runs), &
                 replica_overlaps_real(inum_runs, inum_runs), &
#ifdef __CMPLX
                 replica_overlaps_imag(inum_runs, inum_runs), &
#endif
                 tSpinCoupProjE(inum_runs), &

                 NoatDoubs(inum_runs), AllNoatDoubs(inum_runs), &
                 AccRat(inum_runs), &
                 AllHFCyc(inum_runs), OldAllHFCyc(inum_runs), &
                 ENumCyc(inum_runs), AllENumCyc(inum_runs), &
                 ENumCycAbs(inum_runs), AllENumCycAbs(inum_runs), &
                 ProjECyc(inum_runs), &
                 AllGrowRate(inum_runs), &
                 SumWalkersCyc(inum_runs), AllSumWalkersCyc(Inum_runs), &
                 OldAllAvWalkersCyc(inum_runs), &

                 ! Overall wavefunction properties
                 norm_psi(inum_runs), norm_psi_squared(inum_runs), &
                 all_norm_psi_squared(inum_runs), &
                 norm_semistoch(inum_runs), norm_semistoch_squared(inum_runs),&
                 curr_S2(inum_runs), curr_S2_init(inum_runs), &

                 ! Wavefunction output values
                 SumENum(inum_runs), AllSumENum(inum_runs), &
                 ProjectionE(inum_runs), &
                 proje_iter(inum_runs), &
                 AbsProjE(inum_runs), &
                 trial_numerator(inum_runs), tot_trial_numerator(inum_runs), &
                 trial_denom(inum_runs), tot_trial_denom(inum_runs), &
                 sum_proje_denominator(inum_runs), &
                 all_sum_proje_denominator(inum_runs), &
                 cyc_proje_denominator(inum_runs), &
                 all_cyc_proje_denominator(inum_runs), &

                 ! Control variables
                 iBlockingIter(inum_runs), &
                 TargetGrowRate(inum_runs), TargetGrowRateWalk(inum_runs), &

                 ! Variables to do with the shift
                 VaryShiftCycles(inum_runs), VaryShiftIter(inum_runs), &
                 HFShift(inum_runs), &
                 InstShift(inum_runs), &
                 AvDiagSft(inum_runs), SumDiagSft(inum_runs), &
                 DiagSft(inum_runs), &
                 DiagSftRe(inum_runs), &
                 DiagSftIm(inum_runs), &
                 tSinglePartPhase(inum_runs), stat=ierr)

        ! Iteration data
        call allocate_iter_data(iter_data_fciqmc)

        ! real-time FCIQMC: keep track of first and second Runge-Kutta step 
        ! seperately, think of which stats i need for it!
        ! maybe move that to real-time init module..
#ifdef __REALTIME
        allocate(NoAborted_1(lenof_sign), AllNoAborted_1(lenof_sign), &
                 AllNoAbortedOld_1(lenof_sign), &
                 NoRemoved_1(lenof_sign), AllNoRemoved_1(lenof_sign), &
                 AllNoRemovedOld_1(lenof_sign), &

                 ! initiator related: 
                 NoAddedInitiators_1(lenof_sign), InitRemoved_1(lenof_sign), &
                 AllNoAddedInitiators_1(lenof_sign), AllInitRemoved_1(lenof_sign), &

                 ! dynamics: 
                 NoBorn_1(inum_runs), AllNoBorn_1(inum_runs), &
                 NoDied_1(inum_runs), AllNoDied_1(inum_runs), &
                 Annihilated_1(inum_runs), AllAnnihilated_1(inum_runs), &
                 Acceptances_1(inum_runs), &
                 SpawnFromSing_1(inum_runs), AllSpawnFromSing_1(inum_runs), &
                 NoatDoubs_1(inum_runs), AllNoatDoubs_1(inum_runs), &
                 AccRat_1(inum_runs), &
                 AllGrowRateAbort_1(inum_runs), AllGrowRate_1(inum_runs), &

                 ! additional, maybe unused stat. vars. 
                 NoInitDets_1(lenof_sign), NoNonInitDets_1(lenof_sign), &
                 NoInitWalk_1(lenof_sign), NoNonInitWalk_1(lenof_sign), &
                 NoatHF_1(lenof_sign), AllNoInitDets_1(lenof_sign), &
                 AllNoNonInitDets_1(lenof_sign), AllNoInitWalk_1(lenof_sign), &
                 AllNoNonInitWalk_1(lenof_sign), SumWalkersCyc_1(inum_runs), &
                 TotParts_1(lenof_sign), AllTotParts_1(lenof_sign), &
                 AllTotPartsOld_1(lenof_sign), AllNoatHF_1(lenof_sign), &
                 AllSumWalkersCyc_1(inum_runs), OldAllAvWalkersCyc_1(inum_runs), &
                 stat = ierr)
#endif
        ! KPFCIQMC
        allocate(TotPartsInit(lenof_sign), &
                 AllTotPartsInit(lenof_sign), &
                 tSinglePartPhaseKPInit(inum_runs), stat=ierr)

        ! Hacky bugfixes, for variables that aren't clearly set elsewhere.
        VaryShiftIter = 0

    end subroutine

    subroutine clean_replica_arrays()

        ! The reverse of the above routine...

        deallocate(InstNoatHF, &
                   SumNoatHF, AllSumNoatHF, &
                   NoatHF, AllNoatHF, OldAllNoatHF, &
                   iRefProc, proje_ref_energy_offsets, &
                   iHighestPop, &
                   replica_overlaps_real, &
#ifdef __CMPLX
                   replica_overlaps_imag, &
#endif
                   tSpinCoupProjE, &

                   TotParts, AllTotParts, &
                   TotPartsOld, AllTotPartsOld, &
                   HFCyc, &
                   proje_denominator_cyc, proje_denominator_sum, &
                   InitialPartVec, &

                   NoAborted, AllNoAborted, AllNoAbortedOld, &
                   NoRemoved, AllNoRemoved, AllNoRemovedOld, &

                   ! Initiator related
                   NoAddedInitiators, &
                   NoInitDets, NoNonInitDets, &
                   NoInitWalk, NoNonInitWalk, &
                   NoExtraInitDoubs, InitRemoved, &

                   ! And initiator accumulator related
                   AllNoAddedInitiators, &
                   AllNoInitDets, AllNoNonInitDets, &
                   AllNoInitWalk, AllNoNonInitWalk, &
                   AllNoExtraInitDoubs, AllInitRemoved, &

                   ! Track the dynamics
                   NoBorn, AllNoBorn, &
                   NoDied, AllNoDied, &
                   Annihilated, AllAnnihilated, &
                   Acceptances, &
                   SpawnFromSing, AllSpawnFromSing, &
                   AllGrowRateAbort, &

                   NoatDoubs, AllNoatDoubs, &
                   AccRat, &
                   AllHFCyc, OldAllHFCyc, &
                   ENumCyc, AllENumCyc, ENumCycAbs, AllENumCycAbs, &
                   ProjECyc, &
                   AllGrowRate, &
                   SumWalkersCyc, AllSumWalkersCyc, &
                   OldAllAvWalkersCyc, &

                   norm_psi, norm_psi_squared, &
                   all_norm_psi_squared, &
                   norm_semistoch, norm_semistoch_squared, &
                   curr_S2, curr_S2_init, &
                   
                   SumENum, AllSumENum, ProjectionE, &
                   proje_iter, AbsProjE, &
                   trial_numerator, tot_trial_numerator, &
                   trial_denom, tot_trial_denom, &
                   sum_proje_denominator, all_sum_proje_denominator, &
                   cyc_proje_denominator, all_cyc_proje_denominator, &

                   ! Control variables
                   iBlockingIter, &
                   TargetGrowRate, TargetGrowRateWalk, &

                   ! Variables to do with the shift
                   VaryShiftCycles, VaryShiftIter, &
                   HFShift, &
                   InstShift, &
                   AvDiagSft, SumDiagSft, &
                   DiagSft, &
                   tSinglePartPhase, &

                   ! KPFCIQMC
                   TotPartsInit, &
                   AllTotPartsInit, &
                   tSinglePartPhaseKPInit)

                   ! real-time fciqmc
#ifdef __REALTIME 
                   deallocate(NoAborted_1, AllNoAborted_1, AllNoAbortedOld_1, &
                       NoRemoved_1, AllNoRemoved_1, AllNoRemovedOld_1, &
                       NoAddedInitiators_1, AllNoAddedInitiators_1, &
                       InitRemoved_1, AllInitRemoved_1, NoBorn_1, AllNoBorn_1, &
                       NoDied_1, AllNoDied_1, Annihilated_1, AllAnnihilated_1, &
                       Acceptances_1, SpawnFromSing_1, AllSpawnFromSing_1, &
                       NoatDoubs_1, AllNoatDoubs_1, AccRat_1, AllGrowRate_1, &
                       AllGrowRateAbort_1)
#endif

        call clean_iter_data(iter_data_fciqmc)

    end subroutine

    subroutine allocate_iter_data(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: ierr

        allocate(iter_data%nborn(lenof_sign), &
                 iter_data%ndied(lenof_sign), &
                 iter_data%nannihil(lenof_sign), &
                 iter_data%naborted(lenof_sign), &
                 iter_data%nremoved(lenof_sign), &
                 iter_data%update_growth(lenof_sign), &
                 iter_data%update_growth_tot(lenof_sign), &
                 iter_data%tot_parts_old(lenof_sign), stat=ierr)

    end subroutine

    subroutine clean_iter_data(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data

        deallocate(iter_data%nborn, &
                   iter_data%ndied, &
                   iter_data%nannihil, &
                   iter_data%naborted, &
                   iter_data%nremoved, &
                   iter_data%update_growth, &
                   iter_data%update_growth_tot, &
                   iter_data%tot_parts_old)

    end subroutine

    subroutine set_initial_global_data(ndets, ilut_list)

        use bit_rep_data, only: NIfTot, NIfDBO, extract_sign
        use Parallel_neci, only: iProcIndex, MPISumAll

        ! Take in a list of determinants and calculate and set all of the
        ! global data needed for the start of a FCIQMC calculation.

        integer(int64), intent(in) :: ndets
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot,ndets)

        integer :: i, run
        real(dp) :: real_sign(lenof_sign)
        character(*), parameter :: t_r = 'set_initial_global_data'

        TotParts = 0.0_dp
        NoAtHF = 0.0_dp
        iHighestPop = 0

        ! First, find the population of the walkers in walker_list.
        do i = 1, ndets
            call extract_sign(ilut_list(:,i), real_sign)
            TotParts = TotParts + abs(real_sign)
            if ( all(ilut_list(0:NIfDBO,i) == iLutRef(0:NIfDBO, 1)) ) NoAtHF = real_sign

            do run = 1, inum_runs
                if (abs_sign(real_sign(min_part_type(run):max_part_type(run))) > iHighestPop(run)) then
                    iHighestPop(run) = int(abs_sign(real_sign(min_part_type(run):max_part_type(run))))
                    HighestPopDet(:,run) = ilut_list(:, i)
                end if
            end do
        end do

        TotWalkers = ndets
        TotWalkersOld = TotWalkers
        call MPISumAll(TotWalkers,AllTotWalkers)
        AllTotWalkersOld = AllTotWalkers

        TotPartsOld = TotParts
        call MPISumAll(TotParts, AllTotParts)
        AllTotPartsOld = AllTotParts

        call MPISumAll(NoatHF, AllNoatHF)
        OldAllNoatHF = AllNoatHF

#ifdef __PROG_NUMRUNS
        do run = 1, inum_runs
            OldAllAvWalkersCyc(run) = sum(AllTotParts(min_part_type(run):max_part_type(run)))
        enddo
#else
        OldAllAvWalkersCyc = AllTotParts
#endif

        do run = 1, inum_runs
            OldAllHFCyc(run) = ARR_RE_OR_CPLX(AllNoatHF, run)
        end do

        AllNoAbortedOld(:) = 0.0_dp

        iter_data_fciqmc%tot_parts_old = AllTotParts
        
        do run = 1, inum_runs

            ! Calculate the projected energy for this iteration.
            if (ARR_RE_OR_CPLX(AllSumNoAtHF,run) /= 0) &
                ProjectionE(run) = AllSumENum(run) / ARR_RE_OR_CPLX(AllSumNoatHF,run)
            
            ! Keep track of where the particles are
            if (iProcIndex == iRefProc(run)) then
                SumNoatHF(run) = AllSumNoatHF(run)
                SumENum(run) = AllSumENum(run)
                InstNoatHF(run) = NoatHF(run)
            end if

        enddo 

    end subroutine set_initial_global_data


end module
