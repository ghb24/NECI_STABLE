#include "macros.h"
!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

! AJWT
! Bringing you a better FciMCPar.  A vision for the future...
!
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
                          tMomInv
    use bit_reps, only: NIfD, NIfTot, NIfDBO, NIfY, decode_bit_det, &
                        encode_bit_rep, encode_det, extract_bit_rep, &
                        test_flag, set_flag, extract_flags, &
                        flag_is_initiator, clear_all_flags, set_flag_general,&
                        extract_sign, nOffSgn
    use CalcData, only: InitWalkers, NMCyc, DiagSft, Tau, SftDamp, StepsSft, &
                        OccCASorbs, VirtCASorbs, tFindGroundDet, NEquilSteps,&
                        tReadPops, tRegenDiagHEls, iFullSpaceIter, MaxNoAtHF,&
                        GrowMaxFactor, CullFactor, tStartSinglePart, tCCMC, &
                        ScaleWalkers, HFPopThresh, tTruncCAS, NoMCExcits, &
                        tTruncInitiator, tDelayTruncInit, IterTruncInit, &
                        NShiftEquilSteps, tWalkContGrow, tMCExcits, &
                        tAddToInitiator, InitiatorWalkNo, tInitIncDoubs, &
                        tRetestAddtoInit, tReadPopsChangeRef, &
                        tReadPopsRestart, tCheckHighestPopOnce, &
                        iRestartWalkNum, tRestartHighPop, FracLargerDet, &
                        tChangeProjEDet, tCheckHighestPop, tSpawnSpatialInit,&
                        MemoryFacInit, tMaxBloom, tTruncNOpen, tFCIMC, &
                        trunc_nopen_max, tSpawn_Only_Init, tSpawn_Only_Init_Grow, &
                        TargetGrowRate, TargetGrowRateWalk
    use HPHFRandExcitMod, only: FindExcitBitDetSym, gen_hphf_excit
    use MomInvRandExcit, only: gen_MI_excit
    use Determinants, only: FDet, get_helement, write_det, &
                            get_helement_det_only, lexicographic_store, &
                            get_lexicographic_dets, DefDet
    USE DetCalcData , only : ICILevel,nDet,Det,FCIDetIndex
    use GenRandSymExcitNUMod, only: gen_rand_excit, GenRandSymExcitNU, &
                                    ScratchSize, TestGenRandSymExcitNU, &
                                    ScratchSize1, ScratchSize2, ScratchSize3,&
                                    init_excit_gen_store,clean_excit_gen_store
    use GenRandSymExcitCSF, only: gen_csf_excit
    use IntegralsData , only : fck,NMax,UMat,tPartFreezeCore,NPartFrozen,NHolesFrozen,tPartFreezeVirt,NVirtPartFrozen,NElVirtFrozen
    use Logging, only: iWritePopsEvery, TPopsFile, iPopsPartEvery, tBinPops, &
                       iWriteHistEvery, tHistEnergies, FCIMCDebug, &
                       IterShiftBlock, AllHistInitPops, BinRange, iNoBins, &
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
                       instant_s2_multiplier, tMCOutput, &
                       tRDMonFly, IterRDMonFly, tHF_S_D_Ref, &
                       RDMExcitLevel, RDMEnergyIter
    use hist, only: init_hist_spin_dist, clean_hist_spin_dist, &
                    hist_spin_dist, ilut_spindist, tHistSpinDist, &
                    write_clear_hist_spin_dist, hist_spin_dist_iter, &
                    test_add_hist_spin_dist_det, add_hist_energies, &
                    add_hist_spawn, tHistSpawn, AllHistogramEnergy, &
                    AllHistogram, HistogramEnergy, Histogram, AllInstHist, &
                    InstHist, HistMinInd, project_spins, calc_s_squared, &
                    project_spin_csfs, calc_s_squared_multi, &
                    calc_s_squared_star
    USE SymData , only : nSymLabels
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE Parallel
    USE FciMCData
    USE AnnihilationMod
    use PopsfileMod
    use DetBitops, only: EncodeBitDet, DetBitEQ, DetBitLT, FindExcitBitDet, &
                         FindBitExcitLevel, countbits, TestClosedShellDet, &
                         FindSpatialBitExcitLevel, IsAllowedHPHF
    use csf, only: get_csf_bit_yama, iscsf, csf_orbital_mask, get_csf_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_spawn_sign, hphf_off_diag_helement_spawn
    use MI_integrals
    use util_mod, only: choose, abs_int_sign, abs_int8_sign, binary_search
    use constants, only: dp, int64, n_int, lenof_sign
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
    use nElRDMMod, only: FinaliseRDM,Fill_ExplicitRDM_this_Iter,calc_energy_from_rdm, &
                         fill_diag_rdm, fill_sings_rdm, fill_doubs_rdm, &
                         Add_RDM_From_IJ_Pair, tCalc_RDMEnergy, Add_StochRDM_Diag, &
                         DeAlloc_Alloc_SpawnedParts


#ifdef __DEBUG                            
    use DeterminantData, only: write_det
#endif

    implicit none

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)
        use Logging, only: PopsfileTimer
        use util_mod, only: get_free_unit
        real(dp) :: Weight, Energyxw
        INTEGER :: error
        LOGICAL :: TIncrement,tWritePopsFound,tSoftExitFound,tSingBiasChange,tPrintWarn
        REAL(4) :: s_start,s_end,etime,tstart(2),tend(2),totaltime
        real(dp) :: TotalTime8
        INTEGER(int64) :: MaxWalkers,MinWalkers
        real(dp) :: AllTotWalkers,MeanWalkers,Inpair(2),Outpair(2)
        integer, dimension(lenof_sign) :: tmp_sgn
        integer :: tmp_int(lenof_sign), i
        real(dp) :: grow_rate

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

        if(tHistSpawn) then 
            Tot_Unique_Dets_Unit = get_free_unit()
            OPEN(Tot_Unique_Dets_Unit,FILE='TOTUNIQUEDETS',STATUS='UNKNOWN')
        endif

        ! Initial output
        call WriteFciMCStatsHeader()
        ! Prepend a # to the initial status line so analysis doesn't pick up
        ! repetitions in the FCIMCStats or INITIATORStats files from restarts.
!        write (6,'("#")', advance='no')
        write (fcimcstats_unit,'("#")', advance='no')
        write (initiatorstats_unit,'("#")', advance='no')
        call WriteFCIMCStats()

        ! Put a barrier here so all processes synchronise before we begin.
        call MPIBarrier(error)

!Start MC simulation...
        TIncrement=.true.   !If TIncrement is true, it means that when it comes out of the loop, it wants to subtract 1 from the Iteration count to get the true number of iterations
        Iter=1

        SumSigns = 0.D0
        SumSpawns = 0.D0

        do while (Iter <= NMCyc .or. NMCyc == -1)
!Main iteration loop...
!            WRITE(6,*) 'Iter',Iter

            if(iProcIndex.eq.root) s_start=etime(tstart)

            if (tCCMC) then
                CALL PerformCCMCCycPar()
            else
                if (.not. (tSpinProject .and. spin_proj_interval == -1)) then
                    call sub_dispatcher_7 (PerformFciMCycPar, &
                                           ptr_excit_generator, &
                                           ptr_attempt_create, &
                                           ptr_get_spawn_helement, &
                                           ptr_encode_child, &
                                           ptr_new_child_stats, &
                                           ptr_attempt_die, &
                                           ptr_iter_data)
                endif
            endif

            ! Are we projecting the spin out between iterations?
            if (tSpinProject .and. (mod(Iter, spin_proj_interval) == 0 .or. &
                                    spin_proj_interval == -1) .and. &
                (tSinglePartPhase .or. .not. disable_spin_proj_varyshift))then
                do i = 1, max(spin_proj_iter_count, 1)
                    call PerformFciMCycPar (generate_excit_spin_proj, &
                                           attempt_create_normal, &
                                           get_spawn_helement_spin_proj, &
                                           null_encode_child, &
                                           new_child_stats_normal, &
                                           attempt_die_spin_proj, &
                                           iter_data_spin_proj)
                enddo
            endif

            if(iProcIndex.eq.root) then
                s_end=etime(tend)
                IterTime=IterTime+(s_end-s_start)
            endif
            
            if (mod(Iter, StepsSft) == 0) then

                ! Has there been a particle bloom this update cycle?
                if(iProcIndex.eq.Root) then
                    if(tMaxBloom) then
                        ! Only print out warnings if it is a larger bloom that has been experienced before.
                        if(abs(iPartBloom) > iMaxBloom) then
                            if(abs(iPartBloom) > InitiatorWalkNo) tPrintWarn=.true.
                            iMaxBloom=abs(iPartBloom)
                        else
                            tPrintWarn=.false.
                        endif
                    else
                        if(abs(iPartBloom) > InitiatorWalkNo) then
                            tPrintWarn=.true.
                        else
                            tPrintWarn=.false.
                        endif
                    endif
                    if (tPrintWarn) then
                        write (6, bloom_warn_string, advance='no') abs(iPartBloom)
                        if (iPartBloom > 0) then
                            write (6, '("double excit.")')
                        else
                            write (6, '("single excit.")')
                        endif
                    endif
                    iPartBloom = 0 ! Max number spawned from an excitation
                endif

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
                !call project_spin_csfs()

                if(tRestart) cycle

                !! Quick hack for spin projection
                !if (tSpinProject .and. (mod(iter/stepssft, spin_proj_interval) == 0 .or. &
                !    spin_proj_interval == 0)) then


                !    call MPIReduce(iter_data_spin_proj%update_growth, &
                !                   MPI_SUM, tmp_int)
                !    if (iProcIndex == Root) then
                !        grow_rate = (sum(tmp_int + &
                !                     iter_data_spin_proj%tot_parts_old)) &
                !                    / real(sum(iter_data_spin_proj%tot_parts_old), dp)
                !        iter_data_spin_proj%tot_parts_old = AllTotPartsOld

                !        if (.not. tSinglePartPhase) then
                !            spin_proj_shift = spin_proj_shift - ((log(grow_rate) * SftDamp) / &
                !                              (spin_proj_gamma * iter_data_spin_proj%update_iters))
                !        endif

                !    endif

                !    call MPIBCast (spin_proj_shift)

                !    iter_data_spin_proj%update_growth = 0
                !    iter_data_spin_proj%update_iters = 0
                !endif


                IF((tTruncCAS.or.tTruncSpace.or.tTruncInitiator).and.(Iter.gt.iFullSpaceIter).and.(iFullSpaceIter.ne.0)) THEN
!Test if we want to expand to the full space if an EXPANDSPACE variable has been set
                    IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
                        IF(iProcIndex.eq.0) WRITE(6,*) "Unable to expand space since histgramming the wavefunction..."
                    ELSE
                        ICILevel=0
                        tTruncSpace=.false.
                        tTruncCAS=.false.
                        IF(tTruncInitiator) tTruncInitiator=.false.
                        IF(iProcIndex.eq.0) THEN
                            WRITE(6,*) "Expanding to the full space on iteration ",Iter
                        ENDIF
                    ENDIF
                ENDIF

                TotalTime8=real(s_end,dp)
                call MPIBCast(TotalTime8)    !TotalTime is local - broadcast to all procs

!This routine will check for a CHANGEVARS file and change the parameters of the calculation accordingly.
                CALL ChangeVars(tSingBiasChange,tSoftExitFound,tWritePopsFound)
                IF(tSoftExitFound) THEN
                    TIncrement=.false.
                    EXIT
                ENDIF
                IF(tTimeExit.and.(TotalTime8.ge.MaxTimeExit)) THEN
                    !Is it time to exit yet?
                    write(6,"(A,F8.2,A)") "Time limit reached for simulation of: ",MaxTimeExit/60.0," minutes - exiting..."
                    tIncrement=.false.
                    EXIT
                ENDIF
                IF(tWritePopsFound) THEN
!We have explicitly asked to write out the POPSFILE from the CHANGEVARS file.
                    CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                ENDIF
                IF(tSingBiasChange) THEN
                    CALL CalcApproxpDoubles()
                ENDIF

                if((PopsfileTimer.gt.0.D0).and.((iPopsTimers*PopsfileTimer).lt.(TotalTime8/3600.D0))) then
                    !Write out a POPSFILE every PopsfileTimer hours
                    if(iProcIndex.eq.Root) then
                        write(6,"(A,F7.3,A)") "Writing out a popsfile after ",iPopsTimers*PopsfileTimer, " hours..."
                    endif
                    call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                    iPopsTimers=iPopsTimers+1
                    if(iProcIndex.eq.Root) then
                        s_end=etime(tend)
                        write(6,"(A,F7.3,A)") "Time taken to write out POPSFILE: ",real(s_end,dp)-TotalTime8," seconds."
                    endif
                endif
            
            ENDIF   !Endif end of update cycle

!            IF(mod(Iter,iWriteBlockingEvery).eq.0) THEN
!                !Every 100 update cycles, write out a new blocking file.
!                IF(tErrorBlocking.and.(Iter.gt.IterStartBlocking)) CALL PrintBlocking(Iter) 
!                IF(tShiftBlocking.and.(Iter.gt.(VaryShiftIter+IterShiftBlock))) CALL PrintShiftBlocking(Iter)
!            ENDIF
            
            if (tHistSpinDist .and. (mod(iter, hist_spin_dist_iter) == 0)) &
                call write_clear_hist_spin_dist (iter, hist_spin_dist_iter)

            IF(TPopsFile.and.(.not.tPrintPopsDefault).and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
!This will write out the POPSFILE if wanted
                CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
            ENDIF
!            IF(TAutoCorr) CALL WriteHistogrammedDets()

            IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
                CALL WriteHistogram()
            ENDIF
            IF(tHistHamil.and.(mod(Iter,iWriteHamilEvery).eq.0)) THEN
                CALL WriteHamilHistogram()
            ENDIF
            IF(tRDMonFly.and.tCalc_RDMEnergy) THEN
                IF((Iter.gt.IterRDMonFly).and.(mod(Iter-IterRDMonFly,RDMEnergyIter).eq.0)) &
                                CALL Calc_Energy_from_RDM()  
            ENDIF

            Iter=Iter+1
!End of MC cycle
        enddo

        IF(TIncrement) Iter=Iter-1     !Reduce the iteration count for the POPSFILE since it is incremented upon leaving the loop (if done naturally)
        IF(TPopsFile) THEN
            CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
        ENDIF
        IF(tCalcFCIMCPsi) THEN
!This routine will actually only print the matrix if tPrintFCIMCPsi is on
            CALL PrintFCIMCPsi()

            IF(tFindCINatOrbs) THEN
!This routine takes the wavefunction Psi, calculates the one electron density matrix, and rotates the HF orbitals to produce a new ROFCIDUMP file.
                CALL RotateOrbs() 
                CALL MPIBarrier(error)
            ENDIF
        ENDIF

!        IF(tErrorBlocking) CALL FinaliseBlocking(Iter)
!        IF(tShiftBlocking) CALL FinaliseShiftBlocking(Iter)

        IF(tHistSpawn) CALL WriteHistogram()

        IF(tHistHamil) CALL WriteHamilHistogram()

        Weight=(0.D0)
        Energyxw=(ProjectionE+Hii)

        IF(tHistEnergies) CALL WriteHistogramEnergies()

        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
        ENDIF

        IF(tRDMonFly) CALL FinaliseRDM()

        IF(tPrintDoubsUEG) THEN
            CALL PrintDoubUEGOccs(DoubsUEG)
        ENDIF

! Print out some load balancing stats nicely to end.
        CALL MPIReduce(TotWalkers,MPI_MAX,MaxWalkers)
        CALL MPIReduce(TotWalkers,MPI_MIN,MinWalkers)
        CALL MPIAllReduce(Real(TotWalkers,dp),MPI_SUM,AllTotWalkers)
        if (iProcIndex.eq.Root) then
            MeanWalkers=AllTotWalkers/nNodes
            write (6,'(/,1X,a55)') 'Load balancing information based on the last iteration:'
            write (6,'(1X,a33,1X,f18.10)') 'Mean number of walkers/processor:',MeanWalkers
            write (6,'(1X,a32,1X,i18)') 'Min number of walkers/processor:',MinWalkers
            write (6,'(1X,a32,1X,i18,/)') 'Max number of walkers/processor:',MaxWalkers
        end if

!Deallocate memory
        CALL DeallocFCIMCMemPar()

        IF(iProcIndex.eq.Root) THEN
            CLOSE(fcimcstats_unit)
            IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(initiatorstats_unit)
            IF(tLogComplexPops) CLOSE(complexstats_unit)
        ENDIF
        IF(TDebug) CLOSE(11)

        if(tHistSpawn) then 
            close(Tot_Unique_Dets_Unit)
        endif

        RETURN

    END SUBROUTINE FciMCPar

    ! **********************************************************
    ! ************************* NOTE ***************************
    ! ANY changes to the following interfaces MUST be replicated in the
    ! interface declarations for the function pointers passed into
    ! FciMCycPar
    ! --> Otherwise BAD things (may) happen at runtime, in a
    !     non-deterministic (but probably segfault) manner.
    ! **********************************************************

    ! These wrapper functions exist only to enforce interfaces at compile
    ! time. And to contain the access to the (slightly hackish) assign_proc
    subroutine set_excit_generator (gen)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            subroutine gen (nI, iLutI, nJ, iLutJ, exFlag, IC, ex, tParity, &
                            pGen, HEl, store)

                use SystemData, only: nel
                use bit_reps, only: niftot
                use GenRandSymExcitNUMod, only: scratchsize
                use constants, only: n_int,dp
                use FciMCData, only: excit_gen_store_type
                implicit none

                integer, intent(in) :: nI(nel) 
                integer(kind=n_int), intent(in) :: iLutI(0:niftot)
                integer, intent(in) :: exFlag
                integer, intent(out) :: nJ(nel) 
                integer(kind=n_int), intent(out) :: iLutJ(0:niftot)
                integer, intent(out) :: ic, ex(2,2)
                real(dp), intent(out) :: pGen
                logical, intent(out) :: tParity
                HElement_t, intent(out) :: HEl
                type(excit_gen_store_type), intent(inout), target :: store
            end subroutine
        end interface

        call assign_proc (ptr_excit_generator, gen)
    end subroutine

    subroutine set_attempt_create (attempt_create)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            function attempt_create (get_spawn_helement, nI, iLutI, wSign, &
                                     nJ, iLutJ, prob, HElGen, ic, ex, tPar, exLevel, &
                                     part_type, wSignDied) result(child)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int, dp, lenof_sign
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel), part_type 
                integer(kind=n_int), intent(in) :: iLutI(0:nIfTot)
                integer(kind=n_int), intent(inout) :: iLutJ(0:nIfTot)
                integer, intent(in) :: ic, ex(2,2), exLevel
                integer, dimension(lenof_sign), intent(in) :: wSign
                integer, dimension(lenof_sign), intent(in), optional :: wSignDied
                logical, intent(in) :: tPar
                real(dp), intent(inout) :: prob
                integer , dimension(lenof_sign) :: child      
                HElement_t , intent(in) :: HElGen

                interface
                    function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, &
                                                 ex, tParity, HElGen) &
                                                 result (hel)
                        use SystemData, only: nel
                        use bit_reps, only: niftot
                        use constants, only: n_int,dp
                        implicit none
                        integer, intent(in) :: nI(nel), nJ(nel)
                        integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
                        integer, intent(in) :: ic, ex(2,2)
                        logical, intent(in) :: tParity
                        HElement_t, intent(in) :: HElGen
                        HElement_t :: hel
                    end function
                end interface
            end function
        end interface
    
        call assign_proc (ptr_attempt_create, attempt_create)
    end subroutine

    subroutine set_get_spawn_helement (get_spawn_helement)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, &
                                         ex, tParity, HElGen) result (hel)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int,dp
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
                integer, intent(in) :: ic, ex(2,2)
                logical, intent(in) :: tParity
                HElement_t, intent(in) :: HElGen
                HElement_t :: hel
            end function
        end interface

        call assign_proc (ptr_get_spawn_helement, get_spawn_helement)
    end subroutine

    subroutine set_encode_child (encode_child)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            subroutine encode_child (ilutI, ilutJ, ic, ex)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int
                implicit none
                integer(kind=n_int), intent(in) :: iLutI(0:nifTot)
                integer, intent(in) :: ic, ex(2,2)
                integer(kind=n_int), intent(inout) :: iLutJ(0:nIfTot)
            end subroutine
        end interface
    
        call assign_proc (ptr_encode_child, encode_child)
    end subroutine

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

    subroutine set_new_child_stats (new_child_stats)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            subroutine new_child_stats (iter_data, iLutI, nJ, iLutJ, ic, &
                                        walkExLevel, child, parent_flags, &
                                        part_type)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int, lenof_sign
                use FciMCData, only: fcimc_iter_data
                implicit none
                integer(kind=n_int), intent(in) :: ilutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
                integer, intent(in) :: part_type
                integer, dimension(lenof_sign) , intent(in) :: child
                type(fcimc_iter_data), intent(inout) :: iter_data
            end subroutine
        end interface

        call assign_proc (ptr_new_child_stats, new_child_stats)
    end subroutine

    subroutine set_attempt_die (attempt_die)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            function attempt_die (nI, Kii, wSign) result(ndie)
                use SystemData, only: nel
                use constants, only: lenof_sign, dp
                implicit none
                integer, intent(in) :: nI(nel)
                integer, dimension(lenof_sign), intent(in) :: wSign
                real(dp), intent(in) :: Kii
                integer, dimension(lenof_sign) :: ndie
            end function
        end interface

        call assign_proc (ptr_attempt_die, attempt_die)
    end subroutine

    subroutine set_fcimc_iter_data (data_struct)
        use, intrinsic :: iso_c_binding
        implicit none
        type(fcimc_iter_data), target :: data_struct

        call assign_proc (ptr_iter_data, data_struct)
    end subroutine
    ! This is the heart of FCIMC, where the MC Cycles are performed.
    !
    ! Note: This should only be called indirectly:
    !       call fn_dispatcher_5 (PerformFCIMCycPar, ptr_...)
    subroutine PerformFCIMCycPar(generate_excitation, attempt_create, &
                                 get_spawn_helement, encode_child, &
                                 new_child_stats, attempt_die, &
                                 iter_data)

        ! **********************************************************
        ! ************************* NOTE ***************************
        ! ANY changes to the following interfaces MUST be replicated in the
        ! interface declarations for the function pointer setting functions
        ! at the top of the module.
        ! --> Otherwise BAD things (may) happen at runtime, in a
        !     non-deterministic (but probably segfault) manner.
        ! **********************************************************
        interface
            subroutine generate_excitation (nI, iLutI, nJ, iLutJ, &
                             exFlag, IC, ex, tParity, pGen, HEl, store)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use GenRandSymExcitNUMod, only: scratchsize
                use constants, only: dp,n_int
                use FciMCData, only: excit_gen_store_type
                implicit none
                integer, intent(in) :: nI(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot)
                integer, intent(in) :: exFlag
                integer, intent(out) :: nJ(nel) 
                integer(kind=n_int), intent(out) :: iLutJ(0:niftot)
                integer, intent(out) :: ic, ex(2,2)
                real(dp), intent(out) :: pGen
                logical, intent(out) :: tParity
                HElement_t, intent(out) :: HEl
                type(excit_gen_store_type), intent(inout), target :: store
            end subroutine
            function attempt_create (get_spawn_helement, nI, iLutI, wSign, &
                                     nJ, iLutJ, prob, HElGen, ic, ex, tPar, exLevel, &
                                     part_type, wSignDied) result(child)
                use systemdata, only: nel
                use bit_reps, only: niftot
                use constants, only: dp, n_int, lenof_sign
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel), part_type
                integer(kind=n_int), intent(in) :: iLutI(0:nifTot)
                integer(kind=n_int), intent(inout) :: iLutJ(0:nIfTot)
                integer, intent(in) :: ic, ex(2,2), exLevel
                integer, dimension(lenof_sign), intent(in) :: wSign
                integer, dimension(lenof_sign), intent(in), optional :: wSignDied
                logical, intent(in) :: tPar
                real(dp), intent(inout) :: prob
                integer, dimension(lenof_sign) :: child
                HElement_t , intent(in) :: HElGen

                interface
                    function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, &
                                                 ex, tParity, HElGen) &
                                                 result (hel)
                        use systemdata, only: nel
                        use bit_reps, only: niftot
                        use constants, only: dp,n_int
                        implicit none
                        integer, intent(in) :: nI(nel), nJ(nel)
                        integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
                        integer, intent(in) :: ic, ex(2,2)
                        logical, intent(in) :: tParity
                        HElement_t, intent(in) :: HElGen
                        HElement_t :: hel
                    end function
                end interface
            end function
            function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, &
                                         ex, tParity, HElGen) &
                                         result (hel)
                use systemdata, only: nel
                use bit_reps, only: niftot
                use constants, only: dp,n_int
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
                integer, intent(in) :: ic, ex(2,2)
                logical, intent(in) :: tParity
                HElement_t, intent(in) :: HElGen
                HElement_t :: hel
            end function
            subroutine encode_child (ilutI, ilutJ, ic, ex)
                use systemdata, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int
                implicit none
                integer(kind=n_int), intent(in) :: ilutI(0:niftot)
                integer, intent(in) :: ic, ex(2,2)
                integer(kind=n_int), intent(inout) :: iLutJ(0:nIfTot)
            end subroutine
            subroutine new_child_stats (iter_data, iLutI, nJ, iLutJ, ic, &
                                        walkExLevel, child, parent_flags, &
                                        part_type)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int, lenof_sign
                use FciMCData, only: fcimc_iter_data
                implicit none
                integer(kind=n_int), intent(in) :: ilutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
                integer, intent(in) :: part_type
                integer, dimension(lenof_sign), intent(in) :: child
                type(fcimc_iter_data), intent(inout) :: iter_data
            end subroutine
            function attempt_die (nI, Kii, wSign) result(ndie)
                use SystemData, only: nel
                use constants, only: lenof_sign, dp
                implicit none
                integer, intent(in) :: nI(nel)
                integer, dimension(lenof_sign), intent(in) :: wSign
                real(dp), intent(in) :: Kii
                integer, dimension(lenof_sign) :: ndie
            end function
        end interface

        ! Iteration specific data
        type(fcimc_iter_data), intent(inout) :: iter_data

        ! Now the local, iteration specific, variables
        integer :: VecSlot, j, p, error, proc_temp, i
        integer :: DetCurr(nel), nJ(nel), FlagsCurr, parent_flags
        integer, dimension(lenof_sign) :: SignCurr, child, DiedSignCurr
        integer(kind=n_int) :: iLutnJ(0:niftot)
        integer :: IC, walkExcitLevel, ex(2,2), TotWalkersNew, part_type, k
        integer(int64) :: tot_parts_tmp(lenof_sign)
        logical :: tParity
        real(dp) :: prob, HDiagCurr, TempTotParts, Di_Sign_Temp
        HElement_t :: HDiagTemp,HElGen
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'PerformFCIMCycPar' 
#endif

        call set_timer(Walker_Time,30)

        IF(tDelayTruncInit.and.(Iter.ge.IterTruncInit)) THEN 
            IF(Iter.eq.IterTruncInit) THEN
                ! Why do this? Why not just get all procs to do division?
                IF(iProcIndex.eq.root) THEN
                    Tau=Tau/10.D0
                    WRITE(6,'(A,F10.5)') 'Beginning truncated initiator calculation and reducing tau by " &
                        &//"a factor of 10. New tau is : ',Tau
                ENDIF
                CALL MPIBCast(Tau)
            ENDIF
            tTruncInitiator=.true.
        ENDIF

        IF(tRDMonFly.and.(Iter.eq.IterRDMonFly)) THEN
            !We have reached the iteration where we want to start filling the RDM.
            if(tExplicitAllRDM) then
                tFillingExplicRDMonFly = .true.
            else
                !By default - we will do a stochastic calculation of the RDM.
                tFillingStochRDMonFly = .true.
            endif
            !RDMExcitLevel of 3 means we calculate both the 1 and 2 RDM's - otherwise we 
            !calculated only the RDMExcitLevel-RDM.
            IF(RDMExcitLevel.eq.3) THEN
                WRITE(6,'(A)') 'Beginning to calculate both the 1 and 2 electron density matrices on the fly.'
            ELSE
                WRITE(6,'(A28,I1,A36)') 'Beginning to calculate the ',RDMExcitLevel,' electron density matrix on the fly.'
            ENDIF
            call DeAlloc_Alloc_SpawnedParts()
        ENDIF

        MaxInitPopPos=0
        MaxInitPopNeg=0
        HighPopNeg=1
        HighPopPos=1

        ! Synchronise processors
        CALL MPIBarrier(error)

        ! Reset iteration variables
        VecSlot = 1    ! Next position to write into CurrentDets
        NoatHF = 0     ! Number at HF and doubles for stats
        NoatDoubs = 0
        ! Next free position in newly spawned list.
        ValidSpawnedList = InitialSpawnedSlots
        
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

        IFDEBUG(FCIMCDebug,3) write(6,"(A)") "TW: Walker  Det"
        do j=1,TotWalkers
            ! N.B. j indicates the number of determinants, not the number
            !      of walkers.

            ! Indicate that the scratch storage used for excitation generation
            ! from the same walker has not been filled (it is filled when we
            ! excite from the first particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.


            if (tSpawn_Only_Init) then
                call extract_sign (CurrentDets(:,j), SignCurr)

                ! TODO: Ensure that the HF determinant has its flags setup
                !       correctly at the start of a run.
                call CalcParentFlag (j, VecSlot, parent_flags)

                ! If we are only spawning from initiators, we don't need to do this decoding
                ! if CurrentDet is a non-initiator otherwise it is necessary either way.

                if (tcurr_initiator)                                         & 
                    ! tcurr_initiator is a global variable which indicates that at least one of the 'types' of
                    ! walker on this determinant is an initiator.
                    ! Decode determinant from (stored) bit-representation.
                    call extract_bit_rep (CurrentDets(:,j), DetCurr, SignCurr, &
                                      FlagsCurr, fcimc_excit_gen_store)
            else                                      
                call extract_bit_rep (CurrentDets(:,j), DetCurr, SignCurr, &
                                  FlagsCurr, fcimc_excit_gen_store)

                if (tTruncInitiator) call CalcParentFlag (j, VecSlot, &
                                                          parent_flags)
            endif                                                          

            !Debug output.
            IFDEBUGTHEN(FCIMCDebug,3)
                if(lenof_sign.eq.2) then
                    write(6,"(A,I10,2I7,I5)",advance='no') "TW:", j,SignCurr,FlagsCurr
                else
                    write(6,"(A,I10,I7,I5)",advance='no') "TW:", j,SignCurr,FlagsCurr
                endif
                call WriteBitDet(6,CurrentDets(:,j),.true.)
                call Flush(6) 
            ENDIFDEBUG

!            call test_sym_excit3 (DetCurr, 1000000, pDoubles, 3)

            ! TODO: The next couple of bits could be done automatically

            ! We only need to find out if determinant is connected to the
            ! reference (so no ex. level above 2 required, 
            ! truncated etc.)
            walkExcitLevel = FindBitExcitLevel (iLutRef, CurrentDets(:,j), &
                                                max_calc_ex_level)

            ! Should be able to make this function pointer-able
            if (tRegenDiagHEls) then
                ! We are not storing the diagonal hamiltonian elements for 
                ! each particle. Therefore, we need to regenerate them.
                if (DetBitEQ(CurrentDets(:,j), iLutHF, NIfDBO) .and. &
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
                HDiagCurr = CurrentH(j)
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
            call SumEContrib (DetCurr, WalkExcitLevel, SignCurr, &
                              CurrentDets(:,j), HDiagCurr, 1.d0)

            call calc_walker_death (attempt_die, iter_data, DetCurr, &
                                        HDiagCurr, SignCurr, DiedSignCurr)

            ! Add in any reduced density matrix info we want while running 
            ! over CurrentDets - i.e. diagonal elements and connections to HF.
            if(tFillingStochRDMonFly) &
                call Add_StochRDM_Diag(CurrentDets(:,j),DetCurr,DiedSignCurr,walkExcitLevel)

            ! Loop over the 'type' of particle. 
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
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
                ! up by noMCExcits if attempting multiple excitations from 
                ! each walker (default 1)
                do p = 1, abs(SignCurr(part_type)) * noMCExcits
                    ! Zero the bit representation, to ensure no extraneous
                    ! data gets through.
                    ilutnJ = 0

                    ! Generate a (random) excitation
                    call generate_excitation (DetCurr, CurrentDets(:,j), nJ, &
                                   ilutnJ, exFlag, IC, ex, tParity, prob, &
                                   HElGen, fcimc_excit_gen_store)

                    ! If a valid excitation, see if we should spawn children.
                    if (.not. IsNullDet(nJ)) then
                        child = attempt_create (get_spawn_helement, DetCurr, &
                                            CurrentDets(:,j), SignCurr, &
                                            nJ,iLutnJ, Prob, HElGen, IC, ex, &
                                            tParity, walkExcitLevel,part_type, &
                                            DiedSignCurr)
                    else
                        child = 0
                    endif

                    ! Children have been chosen to be spawned.
                    if (any(child /= 0)) then
                        ! We know we want to create a particle of this type.
                        ! Encode the bit representation if it isn't already.

                        call encode_child (CurrentDets(:,j), iLutnJ, ic, ex)

                        call new_child_stats (iter_data, CurrentDets(:,j), &
                                              nJ, iLutnJ, ic, walkExcitLevel,&
                                              child, parent_flags, part_type)

                        call create_particle (nJ, iLutnJ, child, &
                                              parent_flags, part_type,& 
                                              CurrentDets(:,j),SignCurr)

                    endif ! (child /= 0). Child created

                enddo ! Cycling over mulitple particles on same determinant.

            enddo   ! Cycling over 'type' of particle on a given determinant.

            ! DEBUG
            ! if (VecSlot > j) call stop_all (this_routine, 'vecslot > j')
            call walker_death (iter_data, CurrentDets(:,j), HDiagCurr, &
                                SignCurr, DiedSignCurr, VecSlot)

        enddo ! Loop over determinants.
        IFDEBUG(FCIMCDebug,2) write(6,*) 'Finished loop over determinants'
        
        ! Since VecSlot holds the next vacant slot in the array, TotWalkers
        ! should be one less than this. TotWalkersNew is now the number of particles
        ! in the main array, before annihilation
        TotWalkersNew = VecSlot - 1

        ! Update the statistics for the end of an iteration.
        call end_iter_stats (TotWalkersNew)
        
        ! Print bloom/memory warnings
        call end_iteration_print_warn (totWalkersNew)
        call halt_timer (walker_time)
        
        ! For the direct annihilation algorithm. The newly spawned 
        ! walkers should be in a seperate array (SpawnedParts) and the other 
        ! list should be ordered.
        call set_timer (annihil_time, 30)

!        WRITE(6,*) 'Spawned before annihilation'
!        do j = 1, ValidSpawnedList(iProcIndex)-1
!            WRITE(6,*) SpawnedParts(:,j)
!        enddo

        call DirectAnnihilation (totWalkersNew, iter_data,.false.) !.false. for not single processor
        TotWalkers=TotWalkersNew
        CALL halt_timer(Annihil_Time)
        IFDEBUG(FCIMCDebug,2) WRITE(6,*) "Finished Annihilation step"

        ! Update iteration data
        iter_data%update_growth = iter_data%update_growth + iter_data%nborn &
                                - iter_data%ndied - iter_data%nannihil &
                                - iter_data%naborted
        iter_data%update_iters = iter_data%update_iters + 1

        ! This routine will take the CurrentDets and search the array to find all single and double 
        ! connections - adding them into the RDM's. 
        ! This explicit way of doing this is very expensive, but o.k for very small systems.
        IF(tFillingExplicRDMonFly) CALL Fill_ExplicitRDM_this_Iter(TotWalkers)

    end subroutine

    subroutine end_iter_stats (TotWalkersNew)

        integer, intent(in) :: TotWalkersNew
        HElement_t :: delta
        integer :: proc, pos, sgn(lenof_sign), i

        ! SumWalkersCyc calculates the total number of walkers over an update
        ! cycle on each process.
        SumWalkersCyc = SumWalkersCyc + int(sum(TotParts), int64)

        ! Write initiator histograms if on the correct iteration.
        if ((tHistInitPops .and. mod(iter, HistInitPopsIter) == 0) &
            .or. tPrintHighPop) then
            call FindHighPopDet (TotWalkersNew)
            root_write(6,'(a)') 'Writing out the spread of the initiator &
                                &determinant populations.'
            call WriteInitPops (iter + PreviousCycles)
        endif

        ! Update the sum for the projected-energy denominatior if projecting
        ! onto a linear combination of determinants.
        if (proje_linear_comb .and. nproje_sum > 1) then
            do i = 1, nproje_sum
                proc = DetermineDetNode (proje_ref_dets(:,i), 0)
                if (iProcIndex == proc) then
                    pos = binary_search (CurrentDets(:,1:TotWalkers), &
                                         proje_ref_iluts(:,i), NIfD+1)
                    if (pos > 0) then
                        call extract_sign (CurrentDets(:,pos), sgn)
                        delta = ARR_RE_OR_CPLX(sgn) * proje_ref_coeffs(i)
                        cyc_proje_denominator = cyc_proje_denominator + delta
                        sum_proje_denominator = sum_proje_denominator + delta
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
        integer, dimension(lenof_sign), intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(n_int) :: iUnused

        ! Write out some debugging information if asked
        IFDEBUG(FCIMCDebug,3) then
#ifdef __CMPLX
            write(6,"(A,2I4,A)", advance='no') &
                               "Creating ", child(1), child(2), " particles: "
#else
            write(6,"(A,I4,A)",advance='no') &
                                         "Creating ", child(1), " particles: "
#endif
            write(6,"(A,2I4,A)",advance='no') &
                                      "Parent flag: ", parent_flags, part_type
            call writebitdet (6, ilutJ, .true.)
            call flush(6)
        endif
       
        ! Count the number of children born
        NoBorn = NoBorn + sum(abs(child))
        iter_data%nborn = iter_data%nborn + abs(child)

        if (ic == 1) SpawnFromSing = SpawnFromSing + sum(abs(child))

        if(iProcIndex.eq.Root) then
            if (sum(abs(child)) > abs(iPartBloom)) then
                iPartBloom = sign(sum(abs(child)), 2*ic - 3)
            endif
        endif

        ! Avoid compiler warnings
        iUnused = iLutI(0); iUnused = iLutJ(0); iUnused = walkExLevel

    end subroutine
                        
    subroutine create_particle (nJ, iLutJ, child, parent_flags, part_type, iLutI, SignI)

        ! Create a child in the spawned particles arrays. We spawn particles
        ! into a separate array, but non-contiguously. The processor that the
        ! newly-spawned particle is going to be sent to has to be determined,
        ! and then it will be put into the appropriate element determined by
        ! ValidSpawnedList

        integer, intent(in) :: nJ(nel)
        integer(kind=n_int), intent(in) :: iLutJ(0:niftot)
        integer(kind=n_int), intent(in), optional :: iLutI(0:niftot)
        integer, dimension(lenof_sign), intent(in), optional :: SignI
        integer, dimension(lenof_sign), intent(in) :: child
        integer, intent(in) :: parent_flags
        ! 'type' of the particle - i.e. real/imag
        integer, intent(in) :: part_type
        integer :: proc, flags, j, BiasFac
        logical :: parent_init

        proc = DetermineDetNode(nJ,0)    ! 0 -> nNodes-1)
        ! We need to include any flags set both from the parent and from the
        ! spawning steps
        flags = ior(parent_flags, extract_flags(ilutJ))

!        WRITE(6,*) 'Encoding',iLutJ
!        WRITE(6,*) 'To position',ValidSpawnedList(proc)
!        WRITE(6,*) 'Parent',iLutI

        call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), iLutJ, &
                            child, flags)

        IF(tFillingStochRDMonFly) THEN                            
            !We are spawning from iLutI to SpawnedParts(:,ValidSpawnedList(proc)).
            !We want to store the parent (D_i) with the spawned child (D_j) so that we can
            !add in Di.Dj to the RDM later on.
            !The parent is NIfDBO integers long, and stored in the second part of the SpawnedParts array 
            !from NIfTot+1 -> NIfTot+1 + NIfDBO.
            SpawnedParts(niftot+1:niftot+nifdbo+1, ValidSpawnedList(proc)) = iLutI(0:nifdbo) 

            !We also need to carry with the child (and the parent), the sign of the parent.
            !In actual fact this is the sign of the parent divided by the probability of generating that pair Di and Dj, to account for the 
            !fact that Di and Dj are not always added to the RDM, but only when Di spawns on Dj.
            !This turns the real RDMBiasFacI into an integer to pass around to the relevant processors.
            SpawnedParts(niftot+nifdbo+2, ValidSpawnedList(proc)) = transfer(RDMBiasFacI,SpawnedParts(niftot+nifdbo+2, ValidSpawnedList(proc)))

        ENDIF

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
                    call set_flag_general(SpawnedParts(:,ValidSpawnedList(proc)),flag_parent_initiator(j),parent_init)
                endif
            enddo
        ENDIF

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        

        ! Sum the number of created children to use in acceptance ratio.
        acceptances = acceptances + sum(abs(child))
    end subroutine

    subroutine end_iteration_print_warn (totWalkersNew)
        
        ! Worker function for PerformFciMCycPar. Prints warnings about 
        ! particle blooms and memory usage.
        integer, intent(in) :: totWalkersNew
        integer :: i
        real(dp) :: rat

        ! Too many particles?
        rat = real(TotWalkersNew,dp) / real(MaxWalkersPart,dp)
        if (rat > 0.95) then
            write (6, '(a)') '*WARNING* - Number of particles/determinants &
                             &has increased to over 95% of MaxWakersPart.'
            call flush(6)
        end if

        ! Are ony of the sublists near the end of their alloted space?
        if (nNodes > 1) then
            do i = 0, nNodes-1
                rat = real(ValidSpawnedList(i) - InitialSpawnedSlots(i),dp) /&
                             real(InitialSpawnedSlots(1), dp)
                if (rat > 0.95) then
                    write (6, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of &
                                     &MaxSpawned.'
                    call flush(6)
                endif
            enddo
        else
            rat = real(ValidSpawnedList(0), dp) / real(MaxSpawned, dp)
            if (rat > 0.95) then
                write (6, '(a)') '*WARNING* - Number of spawned particles has&
                                 & reached over 95% of MaxSpawned.'
                call flush(6)
            endif
        endif

        ! Are we near the end of the spatial initiator list
        if (tSpawnSpatialInit) then
            rat = real(no_spatial_init_dets) / real(max_inits)
            if (rat > 0.95) then
                write(6, '(a)') '*WARNING* - Number of spatial initiators has&
                                & reached over 95% f max_inits.'
                call flush(6)
            endif
        endif

    end subroutine


    subroutine CalcParentFlag(j, VecSlot, parent_flags)
!In the CurrentDets array, the flag at NIfTot refers to whether that determinant *itself* is an initiator or not.    
!We need to decide if this willchange due to the determinant acquiring a certain population, or its population dropping
!below the threshold.
!The CurrentDets(:,j) is the determinant we are currently spawning from, so this determines the ParentInitiator flag
!which is passed to the SpawnedDets array and refers to whether or not the walkers *parent* is an initiator or not.
!A flag of 0 means the determinant is an initiator, and 1 it is a non-initiator.
        integer, intent(in) :: j, VecSlot
        integer, intent(out) :: parent_flags
        integer, dimension(lenof_sign) :: CurrentSign
        integer :: part_type
        logical :: tDetinCAS, parent_init

        call extract_sign (CurrentDets(:,j), CurrentSign)

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
                    if (abs(CurrentSign(part_type)) > InitiatorWalkNo) then
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
                   
                    ! If det. in fixed initiator space, or is the HF det, it 
                    ! must remain an initiator.
                    if (.not. tDetInCas .and. &
                        .not. (DetBitEQ(CurrentDets(:,j), iLutHF, NIfDBO)) &
                        .and. abs(CurrentSign(part_type)) <= InitiatorWalkNo &
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

        if ((tHistInitPops .and. mod(iter, histInitPopsIter) == 0) &
            .or. tPrintHighPop) then
             call HistInitPopulations (CurrentSign(1), VecSlot)
        endif

    end subroutine CalcParentFlag


    SUBROUTINE HistInitPopulations(SignCurr,VecSlot)
        USE FciMCLoggingMOD, only : InitBinMin,InitBinIter
        INTEGER , INTENT(IN) :: VecSlot,SignCurr
        INTEGER :: InitBinNo

        IF(ABS(SignCurr).gt.InitiatorWalkNo) THEN
!Just summing in those determinants which are initiators. 

!Need to figure out which bin to put them in though.
            IF(SignCurr.lt.0) THEN
                InitBinNo=(FLOOR(((log(REAL(ABS(SignCurr))))-InitBinMin)/InitBinIter))+1
                IF((InitBinNo.ge.1).and.(InitBinNo.le.25000)) THEN
                    HistInitPops(1,InitBinNo)=HistInitPops(1,InitBinNo)+1
                ELSE
                    CALL Stop_All('HistInitPopulations','Trying to histogram outside the range of the bins.')
                ENDIF

            ELSE
                InitBinNo=(FLOOR(((log(REAL(SignCurr)))-InitBinMin)/InitBinIter))+1

                IF((InitBinNo.ge.1).and.(InitBinNo.le.25000)) THEN
                    HistInitPops(2,InitBinNo)=HistInitPops(2,InitBinNo)+1
                ELSE
                    CALL Stop_All('HistInitPopulations','Trying to histogram outside the range of the bins.')
                ENDIF
            ENDIF
        ENDIF

        IF(SignCurr.lt.MaxInitPopNeg) THEN
            MaxInitPopNeg=SignCurr
            HighPopNeg=VecSlot
        ENDIF
        IF(SignCurr.gt.MaxInitPopPos) THEN
            MaxInitPopPos=SignCurr
            HighPopPos=VecSlot
        ENDIF


    END SUBROUTINE HistInitPopulations


    SUBROUTINE FindHighPopDet(TotWalkersNew)
        USE constants, only : MpiDetInt
!Found the highest population on each processor, need to find out which of these has the highest of all.
        INTEGER(KIND=n_int) :: DetPos(0:NIfTot),DetNeg(0:NIfTot)
        INTEGER :: TotWalkersNew,ProcBCastNeg,ProcBCastPos
        integer(int32) :: HighPopInNeg(2),HighPopInPos(2),HighPopoutNeg(2),HighPopoutPos(2)
        INTEGER, DIMENSION(lenof_sign) :: TempSign

!        WRITE(6,*) 'HighPopPos',HighPopPos
!        WRITE(6,*) 'CurrentSign(HighPopPos)',CurrentSign(HighPopPos)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopNeg),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF
        HighPopInNeg(1)=int(TempSign(1),int32)
        HighPopInNeg(2)=int(iProcIndex,int32)

        CALL MPIAllReduceDatatype(HighPopinNeg,1,MPI_MINLOC,MPI_2INTEGER,HighPopoutNeg)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopPos),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF
        HighPopInPos(1)=int(TempSign(1),int32)
        HighPopInPos(2)=int(iProcIndex,int32)

        CALL MPIAllReduceDatatype(HighPopinPos,1,MPI_MAXLOC,MPI_2INTEGER,HighPopoutPos)

!Now, on all processors, HighPopoutPos(1) is the highest positive population, and HighPopoutNeg(1) is the highest negative population.
!HighPopoutPos(2) is the processor the highest population came from.

        IF(iProcIndex.eq.HighPopOutNeg(2)) DetNeg(:)=CurrentDets(:,HighPopNeg)
        IF(iProcIndex.eq.HighPopOutPos(2)) DetPos(:)=CurrentDets(:,HighPopPos)

        !This is a horrible hack, because the process argument should be of type 'integer' - whatever that is,
        !but the highpopoutneg is explicitly an int(4), so that it works with MPI_2INTEGER. Because
        !of the explicit interfaces, we need to do this.
        ProcBCastNeg=HighPopOutNeg(2)
        ProcBCastPos=HighPopOutPos(2)
        CALL MPIBcast(DetNeg,NIfTot+1,ProcBCastNeg)
        CALL MPIBcast(DetPos,NIfTot+1,ProcBCastPos)

        if (iProcIndex == 0) then
            write (6, '(a, i10, a)') 'The most highly populated determinant &
                                  & with the opposite sign to the HF has ', &
                                  HighPopoutNeg(1), ' walkers.'
            call WriteBitDet (6, DetNeg, .true.)

            write (6, '(a, i10, a9)') 'The most highly populated determinant &
                                  & with the same sign as the HF has ', &
                                  HighPopoutPos(1), ' walkers.'
            call WriteBitDet (6, DetPos, .true.)
        endif

        tPrintHighPop=.false.


    END SUBROUTINE FindHighPopDet

    

    subroutine init_fcimc_fn_pointers ()

        ! Call wrapper functions in C to assign a collection of function 
        ! pointers. These will be passed as arguments to FciMCycPar, allowing
        ! it to directly call the correct subroutines for various actions.
        !
        ! The main advantage of this is that it avoids testing all of the
        ! conditionals for every single particle, during every iteration.

        ! Select the excitation generator
        if (tHPHF) then
            call set_excit_generator (gen_hphf_excit)
        elseif (tCSF) then
            call set_excit_generator (gen_csf_excit)
        elseif (tPickVirtUniform) then
            call set_excit_generator (gen_rand_excit3)
        elseif (tMomInv) then
            call set_excit_generator (gen_MI_excit)
        else
            call set_excit_generator (gen_rand_excit)
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
            if (tHPHF .or. tCSF .or. tMomInv) then
                call set_attempt_create (attempt_create_trunc_spawn)
            else
                call set_attempt_create (attempt_create_trunc_spawn_encode)
            endif
        else
            call set_attempt_create (attempt_create_normal)
        endif

        ! In attempt create, how should we evaluate the off diagonal matrix
        ! elements between a parent and its (potentially) spawned offspring?
        if (tCSF) then
            call set_get_spawn_helement (get_csf_helement)
        elseif (tHPHF) then
            if (tGenMatHEL) then
                call set_get_spawn_helement (hphf_spawn_sign)
            else
                call set_get_spawn_helement (hphf_off_diag_helement_spawn)
            endif
        elseif (tMomInv) then
            if (tGenMatHEl) then
                call set_get_spawn_helement (MI_spawn_sign)
            else
                call set_get_spawn_helement (MI_off_diag_helement_spawn)
            endif
        else
            call set_get_spawn_helement (get_helement_det_only)
        endif

        ! Once we have generated the children, do we need to encode them?
        if (.not. (tCSF .or. tHPHF .or. tMomInv)) then
            call set_encode_child (FindExcitBitDet)
        else
            call set_encode_child (null_encode_child)
        endif

        ! What message should we display for a particle bloom?
        if (tAddToInitiator) then
            bloom_warn_string = '("Bloom of more than n_add: &
                                &A max of ", i8, &
                                &" particles created from ")'
        else
            ! Use this variable to store the bloom cutoff level.
            InitiatorWalkNo = 25
            bloom_warn_string = '("Particle blooms of more than 25 in &
                                &iteration ", i14, ": A max of ", i8, &
                                &" particles created in one attempt from ")'
        endif
        iMaxBloom=0 !The largest bloom to date.

        ! Perform the correct statistics on new child particles
        if (tHistHamil) then
            call set_new_child_stats (new_child_stats_hist_hamil)
        else
            call set_new_child_stats (new_child_stats_normal)
        endif

        call set_attempt_die (attempt_die_normal)

        call set_fcimc_iter_data (iter_data_fciqmc)

    end subroutine


    SUBROUTINE CheckOrdering(DetArray,SignArray,NoDets,tCheckSignCoher)
        INTEGER :: NoDets,i,j
        INTEGER(KIND=n_int) :: DetArray(0:NIfTot,1:NoDets)
        INTEGER :: SignArray(1:NoDets),Comp
        LOGICAL :: tCheckSignCoher

        IF(NoDets.gt.0) THEN
            IF(SignArray(1).eq.0) THEN
                WRITE(6,*) "Iter: ",Iter,1
                CALL Stop_All("CheckOrdering","Array has annihilated particles in it...")
            ENDIF
        ENDIF
        do i=2,NoDets
            IF(SignArray(i).eq.0) THEN
                WRITE(6,*) "Iter: ",Iter,i
                CALL Stop_All("CheckOrdering","Array has annihilated particles in it...")
            ENDIF
            Comp=DetBitLT(DetArray(:,i-1),DetArray(:,i),NIfDBO)
            IF(Comp.eq.-1) THEN
!Array is in reverse numerical order for these particles
                do j=max(i-5,1),min(i+5,NoDets)
                    WRITE(6,*) Iter,j,DetArray(:,j),SignArray(j)
                enddo
                CALL Stop_All("CheckOrdering","Array not ordered correctly")
            ELSEIF(Comp.eq.0) THEN
!Dets are the same - see if we want to check sign-coherence
                IF(tCheckSignCoher) THEN
!!This bit checks that there is only one copy of the determinants in the list
                    do j=max(i-5,1),min(i+5,NoDets)
                        WRITE(6,*) Iter,j,DetArray(:,j),SignArray(j)
                    enddo
                    CALL Stop_All("CheckOrdering","Determinant same as previous one...")
                ENDIF
                IF(tCheckSignCoher.and.(SignArray(i-1).ne.SignArray(i))) THEN
!This checks that any multple copies in the list are sign-coherent...
                    do j=i-5,i+5
                        WRITE(6,*) Iter,j,DetArray(:,j),SignArray(j)
                    enddo
                    CALL Stop_All("CheckOrdering","Array not sign-coherent")
                ENDIF
            ENDIF
        enddo

    END SUBROUTINE CheckOrdering
        
    function attempt_create_trunc_spawn (get_spawn_helement, DetCurr,&
                                         iLutCurr, wSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type, &
                                         wSignDied) result(child)
        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        integer, dimension(lenof_sign), intent(in) :: wSign
        integer, dimension(lenof_sign), intent(in), optional :: wSignDied
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        integer, dimension(lenof_sign) :: child
        HElement_t, intent(in) :: HElGen

        interface
            function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, ex, &
                                         tParity, HElGen) result (hel)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: dp,n_int
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, ex(2,2)
                logical, intent(in) :: tParity
                HElement_t :: hel
                HElement_t , intent(in) :: HElGen 
            end function
        end interface

        if (CheckAllowedTruncSpawn (walkExcitLevel, nJ, iLutnJ, IC)) then
            child = attempt_create_normal (get_spawn_helement, DetCurr, &
                               iLutCurr, wSign, nJ, iLutnJ, prob, HElGen, ic, ex, &
                               tParity, walkExcitLevel, part_type, wSignDied)
        else
            child = 0
        endif
    end function

    function attempt_create_trunc_spawn_encode (get_spawn_helement, DetCurr,&
                                         iLutCurr, wSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type, &
                                         wSignDied) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        integer, dimension(lenof_sign), intent(in) :: wSign
        integer, dimension(lenof_sign), intent(in), optional :: wSignDied
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        integer, dimension(lenof_sign) :: child
        HElement_t , intent(in) :: HElGen

        interface
            function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, ex, &
                                         tParity, HElGen) result (hel)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: dp,n_int
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, ex(2,2)
                logical, intent(in) :: tParity
                HElement_t, intent(in) :: HElGen
                HElement_t :: hel
            end function
        end interface

        call EncodeBitDet (nJ, iLutnJ)
        if (CheckAllowedTruncSpawn (walkExcitLevel, nJ, iLutnJ, IC)) then
            child = attempt_create_normal (get_spawn_helement, DetCurr, &
                               iLutCurr, wSign, nJ, iLutnJ, prob, HElGen, ic, ex, &
                               tParity, walkExcitLevel, part_type, wSignDied)
        else
            child = 0
        endif
    end function

    function attempt_create_normal (get_spawn_helement, DetCurr, iLutCurr, &
                                    wSign, nJ, iLutnJ, prob, HElGen, ic, ex, tparity,&
                                    walkExcitLevel, part_type, wSignDied) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! 1 = Real parent particle, 2 = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        integer, dimension(lenof_sign), intent(in) :: wSign
        integer, dimension(lenof_sign), intent(in), optional :: wSignDied
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        integer, dimension(lenof_sign) :: child
        HElement_t , intent(in) :: HElGen

        interface
            function get_spawn_helement (nI, nJ, ilutI, ilutJ, ic, ex, &
                                         tParity, HElGen) result (hel)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: dp,n_int
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, ex(2,2)
                logical, intent(in) :: tParity
                HElement_t, intent(in) :: HElGen
                HElement_t :: hel
            end function
        end interface
        
        real(dp) :: rat, r, p_spawn_rdmfac, p_notlist_rdmfac
        integer :: extracreate, iUnused,j
        logical :: tpspawn_gt1
        HElement_t :: rh
#ifdef __CMPLX
        ! Avoid compiler warnings when compiling the real version.
        integer :: i
        real(dp) :: MatEl
#endif

        ! If we are generating multiple excitotions, then the probability of
        ! spawning on them must be reduced by the number of excitations
        ! generated (i.e. the excitation is likely to arise a factor of
        ! NoMCExcits more often)
        prob = prob * real(NoMCExcits, dp)

        ! In the case of using HPHF, and when tGenMatHEl is on, the matrix
        ! element is calculated at the time of the excitation generation, 
        ! and returned in HElGen. In this case, get_spawn_helement simply
        ! returns HElGen, rather than recomputing the matrix element.
        rh = get_spawn_helement (DetCurr, nJ, iLutCurr, iLutnJ, ic, ex, &
                                 tParity, HElGen)

        !print*, 'p,rh', prob, rh
#ifdef __CMPLX

!We actually want to calculate Hji - take the complex conjugate, rather than swap around DetCurr and nJ.
            rh = CONJG(rh)

            !We are dealing with spawning from real and imaginary elements, and assume
            !that rh is complex
            IF(part_type.eq.1) THEN
                !Real parent particle

                do i=1,lenof_sign
                    !Run over spawnings from both the real and imaginary part of the matrix element

                    IF(i.eq.1) THEN
                        !We want to use the real part of the matrix element to create real walkers
                        MatEl=REAL(rh,dp)
                    ELSE
                        !We want to use the imaginary part of the matrix element to create imaginary walkers
                        MatEl=AIMAG(rh)
                    ENDIF

                    !Attempt spawning
                    rat = tau * abs(MatEl / prob)

                    ! If probability > 1, then we just create multiple children at the
                    ! chosen determinant.
                    extraCreate = int(rat)
                    rat = rat - real(extraCreate, dp)

                    ! Stochastically choose whether to create or not.
                    r = genrand_real2_dSFMT ()
                    if (rat > r) then
                        !Create child
                        child(i) = -nint(sign(1.0_dp, wsign(part_type)*MatEl))  !Will return +- one depending on the desired sign of the stochastically created child.
                        child(i) = child(i) + sign(extraCreate, child(i))
                    else
                        !Just return if any extra particles created
                        child(i) = -extraCreate*nint(sign(1.0_dp, wsign(part_type)*MatEl))
                    endif
                enddo

            ELSE
                !Imaginary parent particle - rules are slightly different...
                !Attempt to spawn REAL walkers with prob +AIMAG(Hij)/P
                !Attempt to spawn IMAG walkers with prob -REAL(Hij)/P
                do i=1,lenof_sign
                    !Run over spawnings from both the real and imaginary part of the matrix element

                    IF(i.eq.1) THEN
                        !We want to use the imaginary part of the matrix element to create real walkers
                        MatEl=AIMAG(rh)
                    ELSE
                        !We want to use the real part of the matrix element to create imaginary walkers
                        MatEl=REAL(rh,dp)
                    ENDIF

                    !Attempt spawning
                    rat = tau * abs(MatEl / prob)

                    ! If probability > 1, then we just create multiple children at the
                    ! chosen determinant.
                    extraCreate = int(rat)
                    rat = rat - real(extraCreate, dp)

                    ! Stochastically choose whether to create or not.
                    r = genrand_real2_dSFMT ()
                    IF(i.eq.1) THEN
                        !Prob = +AIMAG(Hij)/P to create real children
                        if (rat > r) then
                            !Create child
                            child(i) = nint(sign(1.0_dp, wsign(part_type)*MatEl))  !Will return +- one depending on the desired sign of the stochastically created child.
                            child(i) = child(i) + sign(extraCreate, child(i))
                        else
                            !Just return if any extra particles created
                            child(i) = extraCreate*nint(sign(1.0_dp, wsign(part_type)*MatEl))
                        endif
                    ELSE
                        !Prob = -REAL(Hij)/P to create imaginary children
                        if (rat > r) then
                            !Create child
                            child(i) = -nint(sign(1.0_dp, wsign(part_type)*MatEl))  !Will return +- one depending on the desired sign of the stochastically created child.
                            child(i) = child(i) + sign(extraCreate, child(i))
                        else
                            !Just return if any extra particles created
                            child(i) = -extraCreate*nint(sign(1.0_dp, wsign(part_type)*MatEl))
                        endif
                    ENDIF
                enddo

            ENDIF   ! Type of parent

#else
            !We are dealing with real particles always here.

            rat = tau * abs(rh / prob)

            ! If probability > 1, then we just create multiple children at the
            ! chosen determinant.
            extraCreate = int(rat)
            rat = rat - real(extraCreate, dp)

            ! Stochastically choose whether to create or not.
            r = genrand_real2_dSFMT ()
            if (rat > r) then
                !Create child
                child(1) = -nint(sign(1.0_dp, wsign(1)*real(rh,dp)))
                child(1) = child(1) + sign(extraCreate, child(1))
            else
                !Just return if any extra particles created
                child(1) = -extraCreate*nint(sign(1.0_dp, wsign(1)*real(rh,dp)))
            endif

            ! Avoid compiler warnings
            iUnused = part_type

            if(tFillingStochRDMonFly.and.(child(1).ne.0)) then

                ! We eventually turn this real bias factor into an integer to be passed around 
                ! with the spawned children and their parents - this only works with 64 bit at the mo.
                if(n_int.eq.4) CALL Stop_All('attempt_create_normal', &
                                'the bias factor currently does not work with 32 bit integers.')

                ! Otherwise calculate the 'sign' of Di we are eventually going to add in as Di.Dj.                                
                ! Because we only add in Di.Dj when we successfully spawn from Di.Dj, we need to unbias (scale up) 
                ! Di by the probability of this happening.
                ! We need the probability that the determinant i, with population n_i, will spawn on j.
                ! We only consider one instance of a pair Di,Dj, so just want the probability of any of the n_i 
                ! walkers spawning on j.
                ! P_successful_spawn(j | i)[n_i] =  1 - P_not_spawn(j | i)[n_i]
                ! P_not_spawn(j | i )[n_i] is the probability of none of the n_i walkers spawning on j from i.
                ! This requires either not generating j, or generating j and not succesfully spawning, n_i times.
                ! P_not_spawn(j | i )[n_i] = [(1 - P_gen(j | i)) + ( P_gen( j | i ) * (1 - P_spawn(j | i))]^n_i

                if(extraCreate.ne.0) then
                    ! This is the special case whereby if P_spawn(j | i) > 1, then we will definitely spawn from i -> j.
                    ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                    ! We don't care about multiple spawns - if it's in the list, it gets added in regardless of 
                    ! the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                    p_spawn_rdmfac = 1.D0
                else
!                    p_spawn_rdmfac = tau * abs( real(rh,dp) / prob )
                    p_spawn_rdmfac = rat
                endif
                p_notlist_rdmfac = ( 1.D0 - prob ) + ( prob * (1.D0 - p_spawn_rdmfac) )

                ! The bias fac is now n_i / P_successful_spawn(j | i)[n_i]
                ! However when we add in Di -> Dj, we also add in Dj -> Di, so the probability of generating this pair 
                ! is twice that of just generating Di -> Dj.
                RDMBiasFacI = real(wSignDied(1),dp) / abs( ( 1.D0 - ( p_notlist_rdmfac ** (abs(real(wSign(1),dp)))) ) * 2.D0 )
!                RDMBiasFacI = real(wSign(1),dp) / abs( 1.D0 - ( p_notlist_rdmfac ** (abs(real(wSign(1),dp)))) )
                    
            endif

#endif
        ! Avoid compiler warnings
        iUnused = walkExcitLevel

    end function

!Depreciated function - only used for CCMC
!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,iLutCurr,WSign,nJ,iLutnJ,Prob,IC,Ex,tParity)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        use Logging, only : CCMCDebug
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,ExtraCreate,Ex(2,2),Bin
        INTEGER(KIND=n_int) :: iLutCurr(0:NIfTot),iLutnJ(0:NIfTot)
        LOGICAL :: tParity
        real(dp) :: Prob,r,rat
        integer, dimension(lenof_sign), intent(in) :: wSign
        HElement_t :: rh

        IF(tMCExcits) THEN
!If we are generating multiple excitations, then the probability of spawning on them must be reduced by the number of excitations generated.
!This is equivalent to saying that the excitation is likely to arise a factor of NoMCExcits more often.
            Prob=Prob*REAL(NoMCExcits,dp)
        ENDIF
            

!Calculate off diagonal hamiltonian matrix element between determinants
!        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        IF(tHPHF) THEN
            IF(tGenMatHEl) THEN
!The prob is actually prob/HEl, since the matrix element was generated at the same time as the excitation

                write(6,*) "AttemptCreatePar is a depreciated routine, and is not compatible with HPHF - use attempt_create_normal"
                call stop_all("AttemptCreatePar","This is a depreciated routine.")

!                rat=Tau/abs(Prob)

!                rh=Prob ! to get the signs right for later on.
!                WRITE(6,*) Prob, DetCurr(:),"***",nJ(:)
!                WRITE(6,*) "******"
!                CALL HPHFGetOffDiagHElement(DetCurr,nJ,iLutCurr,iLutnJ,rh)

            ELSE
!The IC given doesn't really matter. It just needs to know whether it is a diagonal or off-diagonal matrix element.
!However, the excitation generator can generate the same HPHF again. If this is done, the routine will send the matrix element back as zero.
                rh = hphf_off_diag_helement (DetCurr, nJ, iLutCurr, iLutnJ)
!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
                rat=Tau*abs(rh)/Prob
!                WRITE(6,*) Prob/rh, DetCurr(:),"***",nJ(:)
!                WRITE(6,*) "******"

            ENDIF
        elseif(tMomInv) then

            write(6,*) "AttemptCreatePar is a depreciated routine, and is not compatible with MomInv - use attempt_create_normal"
            call stop_all("AttemptCreatePar","This is a depreciated routine.")
        ELSE
!Normal determinant spawn

            rh = get_helement (DetCurr, nJ, IC, Ex, tParity)
            !WRITE(6,*) rh

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
            rat=Tau*abs(rh)/Prob
        ENDIF
        IF(CCMCDebug.gt.5) WRITE(6,*) "Connection H-element to spawnee:",rh
!        CALL IsSymAllowedExcit(DetCurr,nJ,IC,Ex,SymAllowed) 
!        IF((.not.SymAllowed).and.(abs(rh).gt.0.D0)) THEN
!            WRITE(17,*) rh
!        ENDIF

!        rhcheck=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!        IF(rh.ne.rhcheck) THEN
!            WRITE(6,*) "DetCurr: ",DetCurr(:)
!            WRITE(6,*) "nJ: ",nJ(:)
!            WRITE(6,*) "EX: ",Ex(1,:),Ex(2,:)
!            WRITE(6,*) "tParity: ",tParity
!            STOP
!        ENDIF

!        IF(abs(rh).le.HEpsilon) THEN
!            AttemptCreatePar=0
!            RETURN
!        ENDIF


!If probability is > 1, then we can just create multiple children at the chosen determinant
        ExtraCreate=INT(rat)
        rat=rat-REAL(ExtraCreate)

!Stochastically choose whether to create or not according to ranlux 
        r = genrand_real2_dSFMT() 
        IF(rat.gt.r) THEN
!            IF(Iter.eq.18925) THEN
!                WRITE(6,*) "Created",rh,rat
!            ENDIF

!Child is created - what sign is it?
            IF(WSign(1).gt.0) THEN
!Parent particle is positive
                IF(real(rh).gt.0.D0) THEN
                    AttemptCreatePar=-1     !-ve walker created
                ELSE
                    AttemptCreatePar=1      !+ve walker created
                ENDIF

            ELSE
!Parent particle is negative
                IF(real(rh).gt.0.D0) THEN
                    AttemptCreatePar=1      !+ve walker created
                ELSE
                    AttemptCreatePar=-1     !-ve walker created
                ENDIF
            ENDIF

        ELSE
!No child particle created
!            IF(Iter.eq.18925) THEN
!                WRITE(6,*) "Not Created",rh,rat
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
                IF(WSign(1).gt.0) THEN
                    IF(real(rh).gt.0.D0) THEN
                        AttemptCreatePar=-ExtraCreate    !Additional particles are negative
                    ELSE
                        AttemptCreatePar=ExtraCreate       !Additional particles are positive
                    ENDIF
                ELSE
                    IF(real(rh).gt.0.D0) THEN
                        AttemptCreatePar=ExtraCreate
                    ELSE
                        AttemptCreatePar=-ExtraCreate
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        
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
!            WRITE(6,"(A,F15.5,I5,G25.15,I8,G25.15)") "Outwards ", rat,ExtraCreate,real(rh),Prob
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
                IF(AttemptCreatePar.ne.0) THEN
                    SinglesHist(Bin)=SinglesHist(Bin)+real(abs(AttemptCreatePar),dp)
                    IF(BRR(Ex(1,1)).le.NEl) THEN
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistOccOcc(Bin)=SinglesHistOccOcc(Bin)+real(abs(AttemptCreatePar),dp)
                        ELSE
                            SinglesHistOccVirt(Bin)=SinglesHistOccVirt(Bin)+real(abs(AttemptCreatePar),dp)
                        ENDIF
                    ELSE
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistVirtOcc(Bin)=SinglesHistVirtOcc(Bin)+real(abs(AttemptCreatePar),dp)
                        ELSE
                            SinglesHistVirtVirt(Bin)=SinglesHistVirtVirt(Bin)+real(abs(AttemptCreatePar),dp)
                        ENDIF
                    ENDIF
                ENDIF
            ELSEIF(IC.eq.2) THEN
                DoublesAttemptHist(Bin)=DoublesAttemptHist(Bin)+(Tau/Prob)
                IF(AttemptCreatePar.ne.0) THEN
                    DoublesHist(Bin)=DoublesHist(Bin)+real(abs(AttemptCreatePar),dp)
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
            IF(AttemptCreatePar.ne.0) THEN
                SpawnHist(Bin)=SpawnHist(Bin)+real(abs(AttemptCreatePar),dp)
!                WRITE(6,*) "Get Here!", real(abs(AttemptCreatePar),dp),Bin
            ENDIF
            AttemptHist(Bin)=AttemptHist(Bin)+(Tau/Prob)
        ENDIF

    END FUNCTION AttemptCreatePar

    subroutine calc_walker_death (attempt_die, iter_data, DetCurr, Kii, &
                                        wSign, CopySign)
! This routine calculates the number of walkers that will die in this iteration, 
! and the consequent new sign for this current determinant, but does not actually 
! change the current sign.
! The new sign is CopySign, and becomes DiedSignCurr in the main routine.
! This is kept separate until the end of the loop over n_i (wSign), at which point 
! walker_death is called to encode the new sign (after death).
        interface
            function attempt_die (nI, Kii, wSign) result(ndie)
                use SystemData, only: nel
                use constants, only: lenof_sign, dp
                implicit none
                integer, intent(in) :: nI(nel)
                integer, dimension(lenof_sign), intent(in) :: wSign
                real(dp), intent(in) :: Kii
                integer, dimension(lenof_sign) :: ndie
            end function
        end interface

        integer, intent(in) :: DetCurr(nel) 
        integer, dimension(lenof_sign), intent(inout) :: wSign
        real(dp), intent(in) :: Kii
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, dimension(lenof_sign), intent(out) :: CopySign
        integer, dimension(lenof_sign) :: iDie

        ! Do particles on determinant die? iDie can be both +ve (deaths), or
        ! -ve (births, if shift > 0)
        iDie = attempt_die (DetCurr, Kii, wSign)

!        iDie(1) = 0     !dmc

!        IF(iDie.ne.0) WRITE(6,*) "Death: ",iDie
        
        IFDEBUG(FCIMCDebug,3) then 
            if(sum(abs(iDie)).ne.0) write(6,"(A,2I4)") "Death: ",iDie(:)
        endif

        ! Update death counter
        iter_data%ndied = iter_data%ndied + min(iDie, abs(wSign))
        NoDied = NoDied + sum(min(iDie, abs(wSign)))

        ! Count any antiparticles
        iter_data%nborn = iter_data%nborn + max(iDie - abs(wSign), 0)
        NoBorn = NoBorn + sum(max(iDie - abs(wSign), 0))

        ! Calculate new number of signed particles on the det.
        CopySign = wSign - (iDie * sign(1, wSign))

    end subroutine

    subroutine walker_death (iter_data, iLutCurr, Kii, wSign, CopySign, VecSlot)
! This is the routine when the walkers are actually killed off, i.e. when the new sign 
! (after walker death) is encoded into the CurrentDets array, and determinants with no 
! walkers left after death are removed from the list.
        integer, dimension(lenof_sign), intent(in) :: wSign, CopySign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        integer, intent(inout) :: VecSlot
        real(dp), intent(in) :: Kii
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, dimension(lenof_sign) :: iDie

        ! Normally slot particles back into main array at position vecslot.
        ! This will normally increment with j, except when a particle dies
        ! completely (so VecSlot <= j, and we can't overwrite a walker we
        ! haven't got to yet).
#ifndef __CMPLX            
        if (CopySign(1).ne.0) then
            IF(tTruncInitiator.and.(sign(1,CopySign(1)).ne.sign(1,wSign(1)))) THEN
                !Abort creation of antiparticles if using initiator
!                WRITE(6,*) "Creating Antiparticles"
                NoAborted=NoAborted+abs(CopySign(1)) 
                iter_data%naborted(1) = iter_data%naborted(1) + abs(CopySign(1))
                if(test_flag(iLutCurr,flag_is_initiator(1))) then
                    NoAddedInitiators=NoAddedInitiators-1.D0
                    if (tSpawnSpatialInit) &
                        call rm_initiator_list (ilutCurr)
                endif

            ELSE
                !Normally we will go in this block
                call encode_bit_rep(CurrentDets(:,VecSlot),iLutCurr,CopySign,extract_flags(iLutCurr))
                if (.not.tRegenDiagHEls) CurrentH(VecSlot) = Kii
                VecSlot = VecSlot + 1
            ENDIF
        elseif(tTruncInitiator) then
            ! All particles on this determinant have gone. If the determinant was an initiator, update the stats
            if(test_flag(iLutCurr,flag_is_initiator(1))) then
                NoAddedInitiators=NoAddedInitiators-1.D0
                if (tSpawnSpatialInit) &
                    call rm_initiator_list (ilutCurr)
            endif
        endif
#else
        !In complex case, fill slot if either real or imaginary particle still there.
        IF((CopySign(1).ne.0).or.(CopySign(2).ne.0)) THEN
            call encode_bit_rep(CurrentDets(:,VecSlot),iLutCurr,CopySign,extract_flags(iLutCurr))
            if (.not.tRegenDiagHEls) CurrentH(VecSlot) = Kii
            VecSlot=VecSlot+1
        ENDIF
#endif            

    end subroutine

    function attempt_die_normal (DetCurr, Kii, wSign) result(ndie)
        
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
        integer, dimension(lenof_sign), intent(in) :: wSign
        real(dp), intent(in) :: Kii
        integer, dimension(lenof_sign) :: ndie

        real(dp) :: r, rat, fac, probsign
        integer :: i, iUnused

        fac = tau * (Kii-DiagSft)

        if(fac.gt.1.D0) then
            if(fac.gt.2.D0) then
                call stop_all("attempt_die_normal","Death probability > 2: Algorithm unstable. Reduce timestep.")
            else
                write(6,"(A,F20.10)") "** WARNING ** Death probability > 1: Creating Antiparticles. Timestep errors possible: ",fac
            endif
        endif

        do i=1,lenof_sign
            ! Subtract the current value of the shift, and multiply by tau.
            ! If there are multiple particles, scale the probability.
            rat = fac * abs(wSign(i))

            ndie(i) = int(rat)

            rat = rat - real(ndie(i), dp)

            ! Choose to die or not stochastically
            r = genrand_real2_dSFMT() 
            if (abs(rat) > r) ndie(i) = ndie(i) + nint(sign(1.0_dp, rat))

        enddo

        ! Avoid compiler warnings
        iUnused = DetCurr(1)

    end function

    

!This routine will find the largest weighted MP1 determinants, from which we can construct energy level splitting dependant on the sign.
    

!This routine will take the particle and sign lists, sort them and then compress them so each determinant is only specified once. This requires sign-coherent lists.
!The 'length' will be returned as the length of the new list.
    SUBROUTINE SortCompressLists(Length,PartList,SignList)
        INTEGER :: Length,SignList(Length)
        INTEGER(KIND=n_int) :: PartList(0:NIfTot,Length)
        INTEGER :: i,VecInd,DetsMerged

        call sort (PartList, SignList)

!Now compress the list.
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
            ELSE
                SignList(VecInd)=SignList(VecInd)+SignList(i)
                DetsMerged=DetsMerged+1
            ENDIF
        enddo

        Length=Length-DetsMerged

    END SUBROUTINE SortCompressLists

!This routine will take the particle and sign lists, sort them and then compress them so each determinant is only specified once. This requires sign-coherent lists.
!The 'length' will be returned as the length of the new list.
!In this version, the hamiltonian matrix elements will be fed through with the rest of the list and taken with the particles.
    SUBROUTINE SortCompressListswH(Length,PartList,SignList,HList)
        INTEGER :: Length,SignList(Length)
        INTEGER(KIND=n_int) :: PartList(0:NIfTot,Length)
        real(dp) :: HList(Length)
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
        use util_mod, only: get_free_unit
        INTEGER :: i, io(8)
        real(dp) :: Norm,EnergyBin

        IF(iProcIndex.eq.Root) THEN
            AllHistogramEnergy(:)=0.D0
            AllAttemptHist(:)=0.D0
            AllSpawnHist(:)=0.D0
            AllDoublesHist(:)=0.D0
            AllDoublesAttemptHist(:)=0.D0
            AllSinglesHist(:)=0.D0
            AllSinglesAttemptHist(:)=0.D0
            AllSinglesHistOccOcc(:)=0.D0
            AllSinglesHistOccVirt(:)=0.D0
            AllSinglesHistVirtOcc(:)=0.D0
            AllSinglesHistVirtVirt(:)=0.D0
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
            Norm=0.D0
            do i=1,iOffDiagNoBins
                Norm=Norm+AllSinglesHist(i)
            enddo
!            WRITE(6,*) "AllSinglesHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllSinglesHist(i)=AllSinglesHist(i)/Norm
            enddo
 
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccOcc(i)=AllSinglesHistOccOcc(i)/Norm
            enddo
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccVirt(i)=AllSinglesHistOccVirt(i)/Norm
            enddo
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtOcc(i)=AllSinglesHistVirtOcc(i)/Norm
            enddo
!            Norm=0.D0
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

            EnergyBin=BinRange/2.D0
            do i=1,iNoBins
                IF(AllHistogramEnergy(i).gt.0.D0) WRITE(io(1),*) EnergyBin, AllHistogramEnergy(i)
                IF(AllAttemptHist(i).gt.0.D0) WRITE(io(2),*) EnergyBin, AllAttemptHist(i)
                IF(AllSpawnHist(i).gt.0.D0) WRITE(io(3),*) EnergyBin, AllSpawnHist(i)
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

            EnergyBin=-OffDiagMax+OffDiagBinRange/2.D0
            do i=1,iOffDiagNoBins
                IF(AllSinglesHist(i).gt.0.D0) WRITE(io(1),*) EnergyBin, AllSinglesHist(i)
                IF(AllSinglesAttemptHist(i).gt.0.D0) WRITE(io(2),*) EnergyBin, AllSinglesAttemptHist(i)
                IF(AllDoublesHist(i).gt.0.D0) WRITE(io(3),*) EnergyBin, AllDoublesHist(i)
                IF(AllDoublesAttemptHist(i).gt.0.D0) WRITE(io(4),*) EnergyBin, AllDoublesAttemptHist(i)
                IF(AllSinglesHistOccOcc(i).gt.0.D0) WRITE(io(5),*) EnergyBin, AllSinglesHistOccOcc(i)
                IF(AllSinglesHistOccVirt(i).gt.0.D0) WRITE(io(6),*) EnergyBin, AllSinglesHistOccVirt(i)
                IF(AllSinglesHistVirtOcc(i).gt.0.D0) WRITE(io(7),*) EnergyBin, AllSinglesHistVirtOcc(i)
                IF(AllSinglesHistVirtVirt(i).gt.0.D0) WRITE(io(8),*) EnergyBin, AllSinglesHistVirtVirt(i)
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
       use util_mod, only: get_free_unit
       integer n,t
       integer c(t+2)
       integer j
       integer icomb
       integer iunit
       iunit = get_free_unit()
       icomb=0
       open(iunit,file='COMBINATIONS',STATUS='UNKNOWN')
!      WRITE(6,*) ' Writing combinations to file COMBINATIONS'
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
!           write(6,*) ' Generated combinations:',ICOMB
           CLOSE(iunit)
           RETURN
       endif
       c(j)=c(j)+1
       goto 20
       
   end subroutine gennct
       

!Similar to WriteHistogram, but will only print out in order of maximum component, and only the averaged wavefunction
    SUBROUTINE PrintFCIMCPsi()
        use DetCalcData , only : FCIDets
        use util_mod, only: get_free_unit
        INTEGER :: i,nI(NEl),ExcitLevel,j, iunit
        real(dp) :: norm,norm1

        CALL MPISumAll(Histogram,AllHistogram)
        norm1=0.D0
        do i=1,Det
            IF(lenof_sign.eq.1) THEN
                norm1=norm1+AllHistogram(1,i)**2
            ELSE
                norm1=norm1+(AllHistogram(1,i)**2)+(AllHistogram(lenof_sign,i)**2)
            ENDIF
        enddo
        norm1=SQRT(norm1)
        WRITE(6,*) "Total FCIMC Wavefuction normalisation:",norm1
        do i=1,Det
            do j=1,lenof_sign
                AllHistogram(j,i)=AllHistogram(j,i)/norm1
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

                norm=0.D0
                do i=1,Det
                    IF(lenof_sign.eq.1) THEN
                        norm=norm+AllHistogram(1,i)**2
                    ELSE
                        norm=norm+(AllHistogram(1,i)**2)+(AllHistogram(lenof_sign,i)**2)
                    ENDIF
!write out FCIMC Component weight (normalised), current normalisation, excitation level
                    ExcitLevel = FindBitExcitLevel(iLutHF, FCIDets(:,i), nel)
                    CALL decode_bit_det(nI,FCIDets(0:NIfTot,i))
                    WRITE(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,norm
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
        use SystemData , only : BasisFN
        use util_mod, only: get_free_unit
        INTEGER :: i,IterRead, io1, io2, io3, Tot_No_Unique_Dets
        real(dp) :: norm,norm1,norm2,norm3,ShiftRead,AllERead,NumParts
        CHARACTER(len=22) :: abstr,abstr2
        LOGICAL :: exists

!This will open a file called SpawnHist-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='SpawnHist-'//adjustl(abstr)
        IF(iProcIndex.eq.0) THEN
            WRITE(6,*) "Writing out the average wavevector up to iteration number: ", Iter
            CALL FLUSH(6)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistogram(:,:)=0.D0
            AllInstHist(:,:)=0.D0
            AllAvAnnihil(:,:)=0.D0
            AllInstAnnihil(:,:)=0.D0
        ENDIF

        CALL MPIReduce(Histogram,MPI_SUM,AllHistogram)
        CALL MPIReduce(InstHist,MPI_SUM,AllInstHist)
        CALL MPIReduce(InstAnnihil,MPI_SUM,AllInstAnnihil)
        CALL MPIReduce(AvAnnihil,MPI_SUM,AllAvAnnihil)
        
        IF(iProcIndex.eq.0) THEN

!            IF(.not.associated(NMRKS)) THEN
!                CALL Stop_All("WriteHistogram","A Full Diagonalization is required in the same calculation before histogramming can occur.")
!            ENDIF

            norm=0.D0
            norm1=0.D0
            norm2=0.D0
            norm3=0.D0
            do i=1,Det
                IF(lenof_sign.eq.1) THEN
                    norm=norm+AllHistogram(1,i)**2
                    norm1=norm1+AllInstHist(1,i)**2
                    norm2=norm2+AllInstAnnihil(1,i)**2
                    norm3=norm3+AllAvAnnihil(1,i)**2
                ELSE
                    norm=norm+(AllHistogram(1,i)**2)+(AllHistogram(lenof_sign,i)**2)
                    norm1=norm1+(AllInstHist(1,i)**2)+(AllInstHist(lenof_sign,i)**2)
                    norm2=norm2+(AllInstAnnihil(1,i)**2)+(AllInstAnnihil(lenof_sign,i)**2)
                    norm3=norm3+(AllAvAnnihil(1,i)**2)+(AllAvAnnihil(lenof_sign,i)**2)
                ENDIF
            enddo
            norm=SQRT(norm)
            norm1=SQRT(norm1)
            norm2=SQRT(norm2)
            norm3=SQRT(norm3)
!            WRITE(6,*) "NORM",norm
            do i=1,Det
                IF(lenof_sign.eq.1) THEN
                    AllHistogram(1,i)=AllHistogram(1,i)/norm
                    AllInstHist(1,i)=AllInstHist(1,i)/norm1
                    IF(norm2.ne.0.D0) THEN
                        AllInstAnnihil(1,i)=AllInstAnnihil(1,i)/norm2
                    ENDIF
                    IF(norm3.ne.0.D0) THEN
                        AllAvAnnihil(1,i)=AllAvAnnihil(1,i)/norm3
                    ENDIF
                ELSE
                    AllHistogram(1,i)=AllHistogram(1,i)/norm
                    AllHistogram(lenof_sign,i)=AllHistogram(lenof_sign,i)/norm
                    AllInstHist(1,i)=AllInstHist(1,i)/norm1
                    AllInstHist(lenof_sign,i)=AllInstHist(lenof_sign,i)/norm1
                    IF(norm2.ne.0.D0) THEN
                        AllInstAnnihil(1,i)=AllInstAnnihil(1,i)/norm2
                        AllInstAnnihil(lenof_sign,i)=AllInstAnnihil(lenof_sign,i)/norm2
                    ENDIF
                    IF(norm3.ne.0.D0) THEN
                        AllAvAnnihil(1,i)=AllAvAnnihil(1,i)/norm3
                        AllAvAnnihil(lenof_sign,i)=AllAvAnnihil(lenof_sign,i)/norm3
                    ENDIF
                ENDIF
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
            OPEN(io2,FILE=abstr2,STATUS='NEW')

            INQUIRE(FILE=abstr,EXIST=exists)
            IF(exists) THEN
                io3 = get_free_unit()
                OPEN(io3,FILE=abstr,STATUS='OLD',POSITION='REWIND',ACTION='READ')
                do while(.true.)
                    READ(io3,"(I13,3G25.16)",END=99) IterRead,ShiftRead,AllERead,NumParts
                    WRITE(io2,"(I13,3G25.16)") IterRead,ShiftRead,AllERead,NumParts
                enddo
99              CONTINUE
                IF(AllHFCyc.eq.0.D0) THEN
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllERead,SUM(AllTotPartsOld)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,SUM(AllTotPartsOld)
                ENDIF
                CLOSE(io2)
                CLOSE(io3)

            ELSE
                OPEN(io2,FILE=abstr2,STATUS='UNKNOWN')
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,SUM(AllTotPartsOld)
                CLOSE(io2)
            ENDIF


            norm=0.D0
            norm1=0.D0
            Tot_No_Unique_Dets = 0
            do i=1,Det
                IF(lenof_sign.eq.1) THEN
                    norm=norm+AllHistogram(1,i)**2
                    norm1=norm1+(AllAvAnnihil(1,i))**2
                ELSE
                    norm=norm+(AllHistogram(1,i)**2)+(AllHistogram(lenof_sign,i)**2)
                    norm1=norm1+(AllAvAnnihil(1,i)**2)+(AllAvAnnihil(lenof_sign,i)**2)
                ENDIF
                IF(lenof_sign.eq.1) THEN
                    WRITE(io1,"(I13,6G25.16)") i,AllHistogram(1,i),norm,AllInstHist(1,i),AllInstAnnihil(1,i),AllAvAnnihil(1,i),norm1
                ELSE
                    WRITE(io1,"(I13,6G25.16)") i,AllHistogram(1,i),norm,AllInstHist(1,i),AllInstAnnihil(1,i),AllAvAnnihil(1,i),norm1
                ENDIF
                IF(AllHistogram(1,i).ne.0.D0) Tot_No_Unique_Dets = Tot_No_Unique_Dets + 1
            enddo
            write(Tot_Unique_Dets_Unit,"(2A20)") Iter, Tot_No_Unique_Dets


!            do i=1,Maxdet
!                bits=0
!                do j=0,nbasis-1
!                    IF(BTEST(i,j)) THEN
!                        Bits=Bits+1
!                    ENDIF
!                enddo
!                IF(Bits.eq.NEl) THEN
!
!                    do j=1,ndet
!                        CALL EncodeBitDet(NMRKS(:,j),iLut)
!                        IF(iLut(0).eq.i) THEN
!                            CALL GETSYM(NMRKS(1,j),NEL,G1,NBASISMAX,ISYM)
!                            IF(ISym%Sym%S.eq.0) THEN
!                                Det=Det+1
!                                WRITE(io1,"(3I12)",advance='no') Det,iLut(0)
!                                HEL=GetHElement3(NMRKS(:,j),NMRKS(:,j),0)
!                                norm=norm+(AllHistogram(i))**2
!                                norm1=norm1+(AllInstHist(i))**2
!                                WRITE(io1,"(5G25.16)") REAL(HEL,8),AllHistogram(i),norm,AllInstHist(i),norm1
!                            ENDIF
!                            EXIT
!                        ENDIF
!                    ENDDO
!                ENDIF
!
!            ENDDO
            CLOSE(io1)
        ENDIF
        InstHist(:,:)=0.D0
        InstAnnihil(:,:)=0.D0

    END SUBROUTINE WriteHistogram

!This routine will write out the average hamiltonian from the spawning run up until now.
    SUBROUTINE WriteHamilHistogram()
        use util_mod, only: get_free_unit
        INTEGER :: i,j
        integer :: iunit
        CHARACTER(len=22) :: abstr

!This will open a file called HamilHist-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='HamilHist-'//adjustl(abstr)
        IF(iProcIndex.eq.0) THEN
            WRITE(6,*) "Writing out the average hamiltonian up to iteration number: ", Iter
            CALL FLUSH(6)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistHamil(:,:)=0.D0
            AllAvHistHamil(:,:)=0.D0
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
        HistHamil(:,:)=0.D0

    END SUBROUTINE WriteHamilHistogram

    subroutine new_child_stats_hist_hamil (iter_data, iLutI, nJ, iLutJ, ic, &
                                           walkExLevel, child, parent_flags, &
                                           part_type)
        ! Based on old AddHistHamilEl. Histograms the hamiltonian matrix, and 
        ! then calls the normal statistics routine.

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        integer, dimension(lenof_sign) , intent(in) :: child
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
        real(dp) :: mean_walkers

        ! Update the total imaginary time passed
        TotImagTime = TotImagTime + StepsSft * Tau

        ! Set Iter time to equal the average time per iteration in the
        ! previous update cycle.
        IterTime = IterTime / real(StepsSft)

        ! Calculate the acceptance ratio
        AccRat = real(Acceptances, dp) / SumWalkersCyc


        IF(lenof_sign.eq.1) THEN
            ! Flip sign of entire ensemble if negative pop on HF
            IF(AllNoatHF(1).lt.0) THEN
                root_print "No. at HF < 0 - flipping sign of entire ensemble &
                           &of particles..."
                root_print AllNoatHF(1)

                call FlipSign ()
                AllNoatHF(1) = -AllNoatHF(1)
                NoatHF(1) = -NoatHF(1)
            endif
        ENDIF

        if (iProcIndex == Root) then
            ! Have all of the particles died?
            if (AllTotwalkers == 0) then
!                write(6,*) AllTotWalkers, TotWalkers
!                call stop_all (this_routine, 'All particles have died. &
!                              &Consider choosing new seed, or raising shift &
!                              &value.')
                write(6,"(A)") "All particles have died. Restarting."
                tRestart=.true.
            else
                tRestart=.false.
            endif
        endif
        call MPIBCast(tRestart)
        if(tRestart) then
!Initialise variables for calculation on each node
            Iter=1
            CALL DeallocFCIMCMemPar()
            IF(iProcIndex.eq.Root) THEN
                CLOSE(fcimcstats_unit)
                IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(initiatorstats_unit)
                IF(tLogComplexPops) CLOSE(complexstats_unit)
            ENDIF
            IF(TDebug) CLOSE(11)
            CALL SetupParameters()
            CALL InitFCIMCCalcPar()
            call WriteFciMCStatsHeader()
            ! Prepend a # to the initial status line so analysis doesn't pick up
            ! repetitions in the FCIMCStats or INITIATORStats files from restarts.
    !        write (6,'("#")', advance='no')
            write (fcimcstats_unit,'("#")', advance='no')
            write (initiatorstats_unit,'("#")', advance='no')
            call WriteFCIMCStats()
            return
        endif

        if(iProcIndex.eq.Root) then
! AJWT dislikes doing this type of if based on a (seeminly unrelated) input option, but can't see another easy way.
!  TODO:  Something to make it better

            if(.not.tCCMC) then
               ! Check how balanced the load on each processor is (even though
               ! we cannot load balance with direct annihilation).
               WalkersDiffProc = MaxWalkersProc - MinWalkersProc
               mean_walkers = AllTotWalkers / real(nNodes,dp)
               if (WalkersDiffProc > nint(mean_walkers / 10.d0) .and. &
                   sum(AllTotParts) > real(nNodes * 500, dp)) then
                   root_write (6, '(a, i13,a,2i11)') &
                       'Potential load-imbalance on iter ',iter,' Min/Max determinants on node: ', &
                       MinWalkersProc,MaxWalkersProc
               endif
            endif

!            ! Deal with blocking analysis - this has been temporarily disabled
!            !
!            ! If we are waiting to start error blocking, check if we should 
!            ! enable it. The conditional causes this test to be skipped once
!            ! blocking has been enabled.            
!            if (tIterStartBlock) then
!                ! If IterStartBlocking is positive, then start blocking when
!                ! we are at that iteration. Otherwise, wait until out of
!                ! fixed shift.
!                if ( (IterStartBlocking > 0 .and. iter > IterStartBlocking) &
!                     .or. (IterStartBlocking <= 0 .and. &
!                           .not. tSinglePartPhase)) then
!                    call InitErrorBlocking (iter)
!                    tIterStartBlock = .false.
!                    tErrorBlocking = .true.
!                endif
!            elseif (tHFPopStartBlock) then
!                if (abs(AllHFCyc) / StepsSft >= HFPopStartBlocking) then
!                    call InitErrorBlocking (iter)
!                    tHFPopStartBlock = .false.
!                    tErrorBlocking = .true.
!                endif
!            endif
!
!            ! If we are waiting to start shift blocking, check if we should
!            ! enable it. The test is skipped once blocknig is enabled.
!            if ((.not. tSinglePartPhase) .and. tInitShiftBlocking .and. &
!                (iter == VaryShiftIter + IterShiftBlock)) then
!                call InitShiftErrorBlocking (iter)
!                tInitShiftBlocking = .false.
!                tShiftBlocking = .true.
!            endif
!
!            ! Perform the blocking at the end of each update
!            if (tErrorBlocking .and. .not. tBlockEveryIteration) &
!                call SumInErrorContrib (iter, AllENumCyc, AllHFCyc)
!            if (tShiftBlocking .and. iter >= VaryShiftIter + IterShiftBlock) &
!                call SumInShiftErrorContrib (iter, DiagSft)
        endif

    end subroutine iter_diagnostics

    subroutine population_check ()
        use HPHFRandExcitMod, only: ReturnAlphaOpenDet
        integer :: pop_highest, proc_highest, pop_change
        integer :: det(nel), i, error
        integer(int32) :: int_tmp(2)
        logical :: tSwapped
        HElement_t :: h_tmp

        if (tCheckHighestPop) then

            ! Obtain the determinant (and its processor) with the highest
            ! population.
            call MPIAllReduceDatatype ((/int(iHighestPop,int32), int(iProcIndex,int32)/), 1, &
                                       MPI_MAXLOC, MPI_2INTEGER, int_tmp)
            pop_highest = int_tmp(1)
            proc_highest = int_tmp(2)

            ! How many walkers do we need to switch dets?
            pop_change = int(FracLargerDet * real(abs_int_sign(AllNoAtHF), dp))
!            write(6,*) "***",AllNoAtHF,FracLargerDet,pop_change, pop_highest,proc_highest
            if (pop_change < pop_highest .and. pop_highest > 50) then

                ! Write out info!
                    root_print 'Highest weighted determinant not reference &
                               &det: ', pop_highest, abs_int_sign(AllNoAtHF)

                ! Are we changing the reference determinant?
                if (tChangeProjEDet) then
                    ! Communicate the change to all dets and print out.
                    call MPIBcast (HighestPopDet(0:NIfTot), NIfTot+1, proc_highest)
                    iLutRef = 0
                    iLutRef(0:NIfDBO) = HighestPopDet(0:NIfDBO)
                    call decode_bit_det (ProjEDet, iLutRef)
                    write (6, '(a)', advance='no') 'Changing projected &
                          &energy reference determinant for the next update cycle to: '
                    call write_det (6, ProjEDet, .true.)

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
                            write(6,"(A)") "Now projecting onto open-shell HPHF as a linear combo of two determinants..."
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
                            write(6,"(A)") "Now projecting onto a momentum-coupled function as a linear combo of two dets..."
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
                    if (tHPHF) then
                        h_tmp = hphf_diag_helement (ProjEDet, iLutRef)
                    elseif(tMomInv) then
                        h_tmp = MI_diag_helement(ProjEDet,iLutRef)
                    else
                        h_tmp = get_helement (ProjEDet, ProjEDet, 0)
                    endif
                    Hii = real(h_tmp, dp)
                    write (6, '(a, g25.15)') 'Reference energy now set to: ',&
                                             Hii

                    ! Reset averages
                    SumENum = 0
                    sum_proje_denominator = 0
                    cyc_proje_denominator = 0
                    SumNoatHF = 0
                    VaryShiftCycles = 0
                    SumDiagSft = 0
                    root_print 'Zeroing all energy estimators.'

                    ! Regenerate all the diagonal elements relative to the
                    ! new reference det.
                    write (6,*) 'Regenerating the stored diagonal HElements &
                                &for all walkers.'
                    do i = 1, Totwalkers
                        call decode_bit_det (det, CurrentDets(:,i))
                        if (tHPHF) then
                            h_tmp = hphf_diag_helement (det, CurrentDets(:,i))
                        elseif(tMomInv) then
                            h_tmp = MI_diag_helement(det,CurrentDets(:,i))
                        else
                            h_tmp = get_helement (det, det, 0)
                        endif
                        CurrentH(i) = real(h_tmp, dp) - Hii
                    enddo

                    ! Reset values introduced in soft_exit (CHANGEVARS)
                    if (tCHeckHighestPopOnce) then
                        tChangeProjEDet = .false.
                        tCheckHighestPop = .false.
                        tCheckHighestPopOnce = .false.
                    endif

                ! Or are we restarting the calculation with the reference 
                ! det switched?
                elseif (tRestartHighPop .and. &
                        iRestartWalkNum < sum(AllTotParts)) then
                    
                    ! Broadcast the changed det to all processors
                    call MPIBcast (HighestPopDet, NIfTot+1, proc_highest)
                    iLutRef = 0
                    iLutRef(0:NIfDBO) = HighestPopDet(0:NIfDBO)

                    call decode_bit_det (ProjEDet, iLutRef)
                    write (6, '(a)', advance='no') 'Changing projected &
                             &energy reference determinant to: '
                    call write_det (6, ProjEDet, .true.)

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
                    write (6, '(a, g25.15)') 'Reference energy now set to: ',&
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
        integer :: int_tmp(5+2*lenof_sign), proc, sgn(lenof_sign), pos, i
        HElement_t :: helem_tmp(2)
        HElement_t :: real_tmp(2)!*lenof_sign)
        integer(int64) :: int64_tmp(9)
        type(fcimc_iter_data) :: iter_data
        integer(int64), dimension(lenof_sign), intent(in) :: tot_parts_new
        integer(int64), dimension(lenof_sign), intent(out) :: tot_parts_new_all
        character(len=*), parameter :: this_routine='collate_iter_data'
    
        ! Communicate the integers needing summation
        call MPIReduce ((/Annihilated, NoAtDoubs, NoBorn, NoDied, HFCyc, &
                          SpawnFromSing, iter_data%update_growth/), &
                          MPI_SUM, int_tmp)
        AllAnnihilated = int_tmp(1)
        AllNoAtDoubs = int_tmp(2)
        AllNoBorn = int_tmp(3)
        AllNoDied = int_tmp(4)
        !AllHFCyc = ARR_RE_OR_CPLX(int_tmp(5:4+lenof_sign))
        AllSpawnFromSing = int_tmp(5+lenof_sign)
        iter_data%update_growth_tot = int_tmp(6+lenof_sign:5+2*lenof_sign)
        if (lenof_sign == 1) then
            AllHFCyc = real(int_tmp(5), dp)
        else
            AllHFCyc = cmplx(int_tmp(5),int_tmp(6), dp)
        endif

        ! Integer summations required for the initiator method
        if (tTruncInitiator) then
            call MPISum ((/NoAborted, NoAddedInitiators, NoInitDets, &
                           NoNonInitDets, NoInitWalk, NoNonInitWalk, &
                           NoExtraInitdoubs, InitRemoved/),&
                          int64_tmp(1:8))
            AllNoAborted = int64_tmp(1)
            AllNoAddedInitiators = int64_tmp(2)
            AllNoInitDets = int64_tmp(3)
            AllNoNonInitDets = int64_tmp(4)
            AllNoInitWalk = int64_tmp(5)
            AllNoNonInitWalk = int64_tmp(6)
            AllNoExtraInitDoubs = int64_tmp(7)
            AllInitRemoved = int64_tmp(8)
        endif

        ! 64bit integers
        call MPISum ((/TotWalkers, norm_psi_squared, TotParts, SumNoatHF, &
                       tot_parts_new/), int64_tmp(1:2+3*lenof_sign))
        AllTotWalkers = int64_tmp(1)
        norm_psi_squared = int64_tmp(2)
        norm_psi = sqrt(real(norm_psi_squared, dp))
        AllTotParts = int64_tmp(3:2+lenof_sign)
        AllSumNoatHF = int64_tmp(3+lenof_sign:2+2*lenof_sign)
        tot_parts_new_all = int64_tmp(3+2*lenof_sign:2+3*lenof_sign)

        ! HElement_t values (Calculates the energy by summing all on HF and 
        ! doubles)
        call MPISum ((/ENumCyc, SumENum/), helem_tmp)
        AllENumCyc = helem_tmp(1)
        AllSumENum = helem_tmp(2)

        ! real(dp) values
        call MPISum((/cyc_proje_denominator, sum_proje_denominator/),real_tmp)
        all_cyc_proje_denominator = real_tmp(1)!(1:lenof_sign)
        all_sum_proje_denominator = real_tmp(2)!(lenof_sign+1:2*lenof_sign)

        ! Max/Min values (check load balancing)
        call MPIReduce (TotWalkers, MPI_MAX, MaxWalkersProc)
        call MPIReduce (TotWalkers, MPI_MIN, MinWalkersProc)

        ! We need the total number on the HF and SumWalkersCyc to be valid on
        ! ALL processors (n.b. one of these is 32bit, the other 64)
        call MPISumAll (NoatHF, AllNoatHF)
        call MPISumAll (SumWalkersCyc, AllSumWalkersCyc)

!        WRITE(6,*) "***",iter_data%update_growth_tot,AllTotParts-AllTotPartsOld
        
#ifdef __DEBUG
        !Write this 'ASSERTROOT' out explicitly to avoid line lengths problems
        if ((iProcIndex == root) .and. .not. tSpinProject .and. &
         (.not.all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld))) then
            call stop_all (this_routine, &
                "Assertation failed: all(iter_data%update_growth_tot.eq.AllTotParts-AllTotPartsOld)")
        endif
#endif
    
    end subroutine collate_iter_data

    subroutine update_shift (iter_data)
     
        type(fcimc_iter_data), intent(in) :: iter_data
        integer(int64) :: tot_walkers
        logical :: tReZeroShift
        real(dp) :: AllGrowRateRe, AllGrowRateIm
        real(dp), dimension(lenof_sign) :: denominator, all_denominator

        integer :: error, i, proc, sgn(lenof_sign), pos

!        call flush(6)
        CALL MPIBarrier(error)
        ! collate_iter_data --> The values used are only valid on Root
        if (iProcIndex == Root) then
            ! Calculate the growth rate
!            WRITE(6,*) "iter_data%nborn: ",iter_data%nborn(:)
!            WRITE(6,*) "iter_data%ndied: ",iter_data%ndied(:)
!            WRITE(6,*) "iter_data%nannihil: ",iter_data%nannihil(:)
!            WRITE(6,*) "iter_data%naborted: ",iter_data%naborted(:)
!            WRITE(6,*) "iter_data%update_growth: ",iter_data%update_growth(:)
!            WRITE(6,*) "iter_data%update_growth_tot: ",iter_data%update_growth_tot(:)
!            WRITE(6,*) "iter_data%tot_parts_old: ",iter_data%tot_parts_old(:)
!            WRITE(6,*) "iter_data%update_iters: ",iter_data%update_iters
!            CALL FLUSH(6)

            AllGrowRate = (sum(iter_data%update_growth_tot &
                           + iter_data%tot_parts_old)) &
                          / real(sum(iter_data%tot_parts_old), dp)

            ! For complex case, obtain both Re and Im parts
            if (lenof_sign == 2) then
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
!            write(6,*) "All GrowRate: ",AllGrowRate,TargetGrowRate

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
            if (TSinglePartPhase) then
! AJWT dislikes doing this type of if based on a (seeminly unrelated) input option, but can't see another easy way.
!  TODO:  Something to make it better
                if(.not.tCCMC) then
                    tot_walkers = int(InitWalkers, int64) * int(nNodes,int64)
                else
                    tot_walkers = int(InitWalkers, int64)
                endif
                if ( (sum(AllTotParts) > tot_walkers) .or. &
                     (abs_int_sign(AllNoatHF) > MaxNoatHF)) then
!                     WRITE(6,*) "AllTotParts: ",AllTotParts(1),AllTotParts(2),tot_walkers
                    write (6, '(a,i13,a)') 'Exiting the single particle growth phase on iteration: ',iter, &
                                 ' - Shift can now change'
                    VaryShiftIter = Iter
                    tSinglePartPhase = .false.
                    if(TargetGrowRate.ne.0.D0) then
                        write(6,"(A)") "Setting target growth rate to 1."
                        TargetGrowRate=0.D0
                    endif
                    if(tSpawn_Only_Init.and.tSpawn_Only_Init_Grow) then
                        !Remove the restriction that only initiators can spawn.
                        write(6,*) "All determinants now with the ability to spawn new walkers."
                        tSpawn_Only_Init=.false.
                    endif
                endif
            elseif (abs_int_sign(AllNoatHF) < (MaxNoatHF - HFPopThresh)) then
                write (6, '(a,i13,a)') 'No at HF has fallen too low - reentering the &
                             &single particle growth phase on iteration',iter,' - particle number &
                             &may grow again.'
                tSinglePartPhase = .true.
                tReZeroShift = .true.
            endif

            ! How should the shift change for the entire ensemble of walkers 
            ! over all processors.
            if ((.not. tSinglePartPhase).or.(TargetGrowRate.ne.0.D0)) then

                !In case we want to continue growing, TargetGrowRate > 0.D0
                ! New shift value
                if(TargetGrowRate.ne.0.D0) then
                    if(sum(AllTotParts).gt.TargetGrowRateWalk) then
                        !Only allow targetgrowrate to kick in once we have > TargetGrowRateWalk walkers.
                        DiagSft = DiagSft - (log(AllGrowRate-TargetGrowRate) * SftDamp) / &
                                            (Tau * StepsSft)
                    endif
                else
                    DiagSft = DiagSft - (log(AllGrowRate) * SftDamp) / &
                                        (Tau * StepsSft)
                endif

                if (lenof_sign == 2) then
                    DiagSftRe = DiagSftRe - (log(AllGrowRateRe-TargetGrowRate) * SftDamp) / &
                                            (Tau * StepsSft)
                    DiagSftIm = DiagSftIm - (log(AllGrowRateIm-TargetGrowRate) * SftDamp) / &
                                            (Tau * StepsSft)
                endif

                ! Update the shift averages
                if ((iter - VaryShiftIter) >= nShiftEquilSteps) then
                    if ((iter-VaryShiftIter-nShiftEquilSteps) < StepsSft) &
                        write (6, '(a,i14)') 'Beginning to average shift value on iteration: ',iter
                    VaryShiftCycles = VaryShiftCycles + 1
                    SumDiagSft = SumDiagSft + DiagSft
                    AvDiagSft = SumDiagSft / real(VaryShiftCycles, dp)
                endif

                ! Update DiagSftAbort for initiator algorithm
                if (tTruncInitiator) then
                    DiagSftAbort = DiagSftAbort - &
                              (log(real(AllGrowRateAbort-TargetGrowRate, dp)) * SftDamp) / &
                              (Tau * StepsSft)

                    if (iter - VaryShiftIter >= nShiftEquilSteps) then
                        SumDiagSftAbort = SumDiagSftAbort + DiagSftAbort
                        AvDiagSftAbort = SumDiagSftAbort / &
                                         real(VaryShiftCycles, dp)
                    endif
                endif
            endif

            ! Calculate the instantaneous 'shift' from the HF population
            HFShift = -1.d0 / real(abs_int_sign(AllNoatHF), dp) * &
                              (real(abs_int_sign(AllNoatHF) - abs_int_sign(OldAllNoatHF), dp) / &
                              (Tau * real(StepsSft, dp)))
            InstShift = -1.d0 / sum(AllTotParts) * &
                        ((sum(AllTotParts) - sum(AllTotPartsOld)) / &
                         (Tau * real(StepsSft, dp)))

            ! When using a linear combination, the denominator is summed
            ! directly.
            if (.not. (proje_linear_comb .and. nproje_sum > 1)) then
                all_sum_proje_denominator = ARR_RE_OR_CPLX(AllSumNoatHF)
                all_cyc_proje_denominator = AllHFCyc
            endif

            ! Calculate the projected energy.
            if (any(AllSumNoatHF /= 0) .or. &
                (proje_linear_comb .and. nproje_sum > 1)) then
                ProjectionE = AllSumENum / all_sum_proje_denominator 
!                              ARR_RE_OR_CPLX(all_sum_proje_denominator)
                proje_iter = AllENumCyc / all_cyc_proje_denominator 
!                              ARR_RE_OR_CPLX(all_cyc_proje_denominator)
            endif

            ! If we are re-zeroing the shift
            if (tReZeroShift) then
                DiagSft = 0
                VaryShiftCycles = 0
                SumDiagSft = 0
                AvDiagSft = 0
            endif

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
        if(.not.tSinglePartPhase) TargetGrowRate=0.D0

    end subroutine


    subroutine rezero_iter_stats_each_iter (iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        real*8 :: TempTotParts

        NoInitDets = 0
        NoNonInitDets = 0
        NoInitWalk = 0
        NoNonInitWalk = 0
        InitRemoved = 0

        NoAborted = 0

        iter_data%nborn = 0
        iter_data%ndied = 0
        iter_data%nannihil = 0
        iter_data%naborted = 0

        TempTotParts=REAL(TotParts(1),dp)
        CALL MPIAllReduce(TempTotParts,MPI_SUM,AllTotPartstemp)

    end subroutine


    subroutine rezero_iter_stats_update_cycle (iter_data, tot_parts_new_all)
        
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(int64), dimension(lenof_sign), intent(in) :: tot_parts_new_all
        
        ! Zero all of the variables which accumulate for each iteration.

        IterTime = 0
        SumWalkersCyc = 0
        Annihilated = 0
        Acceptances = 0
        NoBorn = 0
        SpawnFromSing = 0
        NoDied = 0
        ENumCyc = 0
        HFCyc = 0

        ! Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld = TotWalkers
        TotPartsOld = TotParts

        ! Save the number at HF to use in the HFShift
        OldAllNoatHF = AllNoatHF

        ! Also the cumulative global variables
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        AllNoAbortedOld = AllNoAborted

        ! Reset the counters
        iter_data%update_growth = 0
        iter_data%update_iters = 0
        iter_data%tot_parts_old = tot_parts_new_all

        ! Reset the linear combination coefficients
        ! TODO: Need to rethink how/when this is done. This is just for tests
        if (proje_linear_comb) &
            call update_linear_comb_coeffs()


    end subroutine

    subroutine calculate_new_shift_wrapper (iter_data, tot_parts_new)

        type(fcimc_iter_data) :: iter_data
        integer(int64), dimension(lenof_sign), intent(in) :: tot_parts_new
        integer(int64), dimension(lenof_sign) :: tot_parts_new_all

        call collate_iter_data (iter_data, tot_parts_new, tot_parts_new_all)
        call iter_diagnostics ()
        if(tRestart) return
        call population_check ()
        call update_shift (iter_data)
        call WriteFCIMCStats ()
        call rezero_iter_stats_update_cycle (iter_data, tot_parts_new_all)

    end subroutine calculate_new_shift_wrapper

!This routine flips the sign of all particles on the node
    SUBROUTINE FlipSign()
        INTEGER :: i
        INTEGER, DIMENSION(lenof_sign) :: TempSign

        do i=1,TotWalkers
            call extract_sign(CurrentDets(:,i),TempSign)
            TempSign(1)=-TempSign(1)
            call encode_sign(CurrentDets(:,i),TempSign)
        enddo
        
!Reverse the flag for whether the sign of the particles has been flipped so the ACF can be correctly calculated
        TFlippedSign=.not.TFlippedSign
        RETURN
    
    END SUBROUTINE FlipSign

    SUBROUTINE WriteFciMCStatsHeader()

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(6,*) ""
            IF(tTruncInitiator.or.tDelayTruncInit) THEN
                WRITE(initiatorstats_unit,"(A2,A10,10A20)") "# ","1.Step","2.TotWalk","3.Annihil","4.Died", &
                & "5.Born","6.TotUniqDets",&
&               "7.InitDets","8.NonInitDets","9.InitWalks","10.NonInitWalks","11.AbortedWalks"
            ENDIF
            IF(tLogComplexPops) THEN
                WRITE(complexstats_unit,"(A)") '#   1.Step  2.Shift     3.RealShift     4.ImShift   5.TotParts      " &
                & //"6.RealTotParts      7.ImTotParts'
            ENDIF

#ifdef __CMPLX
            if(tMCOutput) then
                write(6, '(a)') "       Step     Shift      WalkerCng(Re)  &
                       &WalkerCng(Im)    TotWalkers(Re)   TotWalkers(Im)    &
                       &Proj.E(Re)   ProjE(Im)     Proj.E.ThisCyc(Re)  &
                       &Proj.E.ThisCyc(Im)   NoatHF(Re)   NoatHF(Im)   &
                       &NoatDoubs      AccRat     UniqueDets     IterTime"
            endif
            write(fcimcstats_unit, "(a,i4,a,l,a,l,a,l)") &
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
                   &28.Inst S^2"
#else
            if(tMCOutput) then
                write(6,"(A)") "       Step     Shift      WalkerCng    &
                      &GrowRate       TotWalkers    Annihil    NoDied    &
                      &NoBorn    Proj.E          Av.Shift     Proj.E.ThisCyc   &
                      &NoatHF NoatDoubs      AccRat     UniqueDets     IterTime"
            endif
            write(fcimcstats_unit, "(a,i4,a,l,a,l,a,l)") &
                  "# FCIMCStats VERSION 2 - REAL : NEl=", nel, &
                  " HPHF=", tHPHF, ' Lz=', tFixLz, &
                  ' Initiator=', tTruncInitiator
            write(fcimcstats_unit, "(A)") &
                  "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     &
                  &5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  &
                  &9.Proj.E       10.Av.Shift 11.Proj.E.ThisCyc  12.NoatHF &
                  &13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime &
                  &17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  &
                  &20.ProjE.ThisIter  21.HFInstShift  22.TotInstShift  &
                  &23.Tot-Proj.E.ThisCyc   24.HFContribtoE  25.NumContribtoE &
                  &26.HF weight    27.|Psi|     28.Inst S^2"
#endif
            
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

    subroutine WriteFCIMCStats()

        ! What is the current value of S2
        ! TODO: This should probably be placed somewhere cleaner.
        if (tCalcInstantS2) then
            if (mod(iter / StepsSft, instant_s2_multiplier) == 0) then
                curr_S2 = calc_s_squared_star (.false.)
            !    curr_S2_init = calc_s_squared_star(.true.)
            endif
        else
            curr_S2 = -1
            curr_S2_init = -1
        endif

        if (iProcIndex == root) then

#ifdef __CMPLX
            write(fcimcstats_unit,"(I12,G16.7,2I10,2I12,4G17.9,3I10,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,8G17.9)") &
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
                real(AllSpawnFromSing, dp) / real(AllNoBorn), &     !17.
                WalkersDiffProc, &                           !18.
                TotImagTime, &                               !19.
                HFShift, &                                   !20.
                InstShift, &                                 !21.
                real((AllHFCyc * conjg(AllHFCyc)),dp), &     !24     |n0|^2  This is the denominator for both calcs
                real((AllENumCyc * conjg(AllHFCyc)),dp), &   !22.    Re[\sum njH0j] x Re[n0] + Im[\sum njH0j] x Im[n0]   No div by StepsSft
                aimag(AllENumCyc * conjg(AllHFCyc)), &       !23.    Im[\sum njH0j] x Re[n0] - Re[\sum njH0j] x Im[n0]   since no physicality
                sqrt(float(sum(AllNoatHF**2))) / norm_psi, &
                norm_psi, &
                curr_S2

            if(tMCOutput) then
                write (6, "(I12,G16.7,2I10,2I12,4G17.9,3I10,G13.5,I12,G13.5)") &
                    Iter + PreviousCycles, &
                    DiagSft, &
                    AllTotParts(1) - AllTotPartsOld(1), &
                    AllTotParts(2) - AllTotPartsOld(2), &
                    AllTotParts(1), AllTotParts(2), &
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
#else

            write(fcimcstats_unit,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,10G17.9)") &
                Iter + PreviousCycles, &
                DiagSft, &
                sum(AllTotParts) - sum(AllTotPartsOld), &
                AllGrowRate, &
                sum(AllTotParts), &
                AllAnnihilated, &
                AllNoDied, &
                AllNoBorn, &
                ProjectionE, &
                AvDiagSft, &
                proje_iter, &
                AllNoatHF, &
                AllNoatDoubs, &
                AccRat, &
                AllTotWalkers, &
                IterTime, &
                real(AllSpawnFromSing) / real(AllNoBorn), &
                WalkersDiffProc, &
                TotImagTime, &
                0.D0, &
                HFShift, &
                InstShift, &
                proje_iter + Hii, &
                AllHFCyc / StepsSft, &
                AllENumCyc / StepsSft, &
                real(AllNoatHF, dp) / norm_psi, &
                norm_psi, &
                curr_S2, curr_S2_init

            if(tMCOutput) then
                write (6, "(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,&
                          &G13.5)") &
                    Iter + PreviousCycles, &
                    DiagSft, &
                    sum(AllTotParts) - sum(AllTotPartsOld), &
                    AllGrowRate, &
                    sum(AllTotParts), &
                    AllAnnihilated, &
                    AllNoDied, &
                    AllNoBorn, &
                    ProjectionE, &
                    AvDiagSft, &
                    proje_iter, &
                    AllNoatHF, &
                    AllNoatDoubs, &
                    AccRat, &
                    AllTotWalkers, &
                    IterTime
            endif
#endif

            if (tTruncInitiator .or. tDelayTruncInit) then
               write(initiatorstats_unit,"(I12,10I20)")&
                   Iter + PreviousCycles, sum(AllTotParts), &
                   AllAnnihilated, AllNoDied, AllNoBorn, AllTotWalkers,&
                   AllNoInitDets, AllNoNonInitDets, AllNoInitWalk, &
                   AllNoNonInitWalk,AllNoAborted
            endif

            if (tLogComplexPops) then
                write (complexstats_unit,"(I12,3G16.7,3I12)") &
                    Iter + PreviousCycles, DiagSft, DiagSftRe, DiagSftIm, &
                    sum(AllTotParts), AllTotParts(1), AllTotParts(lenof_sign)
            endif

            if(tMCOutput) then
                call flush(6)
            endif
            call flush(fcimcstats_unit)
            
        endif

    end subroutine WriteFCIMCStats

    !Ensure that the new FCIMCStats file which is about to be opened does not overwrite any other FCIMCStats
    !files. If there is already an FCIMCStats file present, then move it to FCIMCStats.x, where x is a largest
    !free filename.
    subroutine MoveFCIMCStatsFiles()
        integer :: extension,stat
        logical :: exists
        character(len=22) :: abstr
        character(len=*), parameter :: t_r='MoveFCIMCStatsFiles'

        inquire(file='FCIMCStats',exist=exists)
        if(exists) then
            !We already have an FCIMCStats file - move it to the end of the list of FCIMCStats files.

            extension=0
            do while(.true.)
                abstr=''
                write(abstr,'(I12)') extension
                abstr='FCIMCStats.'//adjustl(abstr)
                inquire(file=abstr,exist=exists)
                if(.not.exists) exit
                extension=extension+1
                if(extension.gt.1000) then
                    call stop_all(t_r,"Error finding free FCIMCStats name")
                endif
            enddo
            
            !We have got a unique filename
            !Do not use system call
!            command = 'mv' // ' FCIMCStats ' // abstr
!            call system(command)

            call rename('FCIMCStats',abstr)
            !Doesn't like the stat argument
!            if(stat.ne.0) then
!                call stop_all(t_r,"Error with renaming FCIMCStats file")
!            endif
        endif

    end subroutine MoveFCIMCStatsFiles

    SUBROUTINE SetupParameters()
        use SystemData, only : tUseBrillouin,iRanLuxLev,tSpn,tHPHFInts,tRotateOrbs,tROHF,tFindCINatOrbs,nOccBeta,nOccAlpha,tUHF
        use SystemData, only : tBrillouinsDefault,ECore,tNoSingExcits
        USE dSFMT_interface , only : dSFMT_init
        use CalcData, only: G_VMC_Seed, &
                            MemoryFacPart, MemoryFacAnnihil, TauFactor, &
                            StepsSftImag, tCheckHighestPop, tSpatialOnlyHash,tStartCAS
        use Determinants , only : GetH0Element3,GetH0Element4
        use SymData , only : SymLabelList,SymLabelCounts,TwoCycleSymGens
        use Logging , only : tTruncRODump
        use DetCalcData, only : NMRKS,tagNMRKS,FCIDets
        use SymExcit3, only : CountExcitations3 
        use constants, only: bits_n_int
        use util_mod, only: get_free_unit
        use HPHFRandExcitMod, only: ReturnAlphaOpenDet
        use sym_mod
        use HElem
        INTEGER :: ierr,i,j,HFDetTest(NEl),Seed,alpha,beta,symalpha,symbeta,endsymstate
        INTEGER :: HFConn,LargestOrb,nBits,HighEDet(NEl)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        HElement_t :: TempHii
        TYPE(BasisFn) HFSym
        real(dp) :: TotDets,SymFactor,r,Gap
        CHARACTER(len=*), PARAMETER :: this_routine='SetupParameters'
        CHARACTER(len=12) :: abstr
        LOGICAL :: tSuccess,tFoundOrbs(nBasis),FoundPair,tSwapped
        INTEGER :: HFLz,ChosenOrb,KPnt(3), step,SymHF

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors
        WRITE(6,*) ""
        WRITE(6,*) "Performing Parallel FCIMC...."
        WRITE(6,*) ""
        
!Set timed routine names
        Walker_Time%timer_name='WalkerTime'
        Annihil_Time%timer_name='AnnihilTime'
        Sort_Time%timer_name='SortTime'
        Comms_Time%timer_name='CommsTime'
        ACF_Time%timer_name='ACFTime'
        AnnSpawned_time%timer_name='AnnSpawnedTime'
        AnnMain_time%timer_name='AnnMainTime'
        BinSearch_time%timer_name='BinSearchTime'

        IF(TDebug) THEN
!This will open a file called LOCALPOPS-"iprocindex" on unit number 11 on every node.
            abstr=''
            write(abstr,'(I2)') iProcIndex
            abstr='LOCALPOPS-'//adjustl(abstr)
            OPEN(11,FILE=abstr,STATUS='UNKNOWN')
        ENDIF

        IF(iProcIndex.eq.Root) THEN
            fcimcstats_unit = get_free_unit()
            if (tReadPops) then
                ! Restart calculation.  Append to stats file (if it exists).
                OPEN(fcimcstats_unit,file='FCIMCStats',status='unknown',position='append')
            else
                call MoveFCIMCStatsFiles()          !This ensures that FCIMCStats files are not overwritten
                OPEN(fcimcstats_unit,file='FCIMCStats',status='unknown')
            end if
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
        CALL LogMemAlloc('HFDet',NEl,4,this_routine,HFDetTag)
        ALLOCATE(iLutHF(0:NIfTot),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for iLutHF")

        do i=1,NEl
            HFDet(i)=FDet(i)
        enddo
        CALL EncodeBitDet(HFDet,iLutHF)

        !iLutRef is the reference determinant for the projected energy.
        !Initially, it is chosen to be the same as the inputted reference determinant
        ALLOCATE(iLutRef(0:NIfTot),stat=ierr)
        ALLOCATE(ProjEDet(NEl),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for iLutRef")
        
        ! The reference / projected energy determinants are the same as the
        ! HF determinant.
        ! TODO: Make these pointers rather than copies?
        iLutRef = iLutHF
        ProjEDet = HFDet

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
                    write(6,*) "HPHF used, and open shell determinant spin-flipped for consistency."
                endif
                write(6,*) "Two *different* determinants contained in initial HPHF"
                write(6,*) "Projected energy will be calculated as projection onto both of these"
                tSpinCoupProjE=.true.
            else
                tSpinCoupProjE=.false.
            endif
            if(tMomInv) then
                call stop_all(this_routine,"Cannot currently have MomInv and HPHF functions. If this is important, bug ghb")
            endif
        elseif(tMomInv) then
            if(tAntisym_MI) then
                write(6,*) "Using hilbert space of antisymmetric momentum-coupled determinants..."
            else
                write(6,*) "Using hilbert space of symmetric momentum-coupled determinants..."
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
                    write(6,*) "Momentum-coupled function used, and reference determinant momentum-flipped for consistency."
                endif
                write(6,*) "Two *different* determinants contained in initial reference function"
                write(6,*) "Projected energy will be calculated as projection onto both of these"
                tSpinCoupProjE=.true.
            endif
            tSpinCoupProjE=.false.
        else
            tSpinCoupProjE=.false.
        endif

!Init hash shifting data
        hash_iter=0

        IF(tKPntSym) THEN
            CALL DecomposeAbelianSym(HFSym%Sym%S,KPnt)
            WRITE(6,"(A,3I5)") "Crystal momentum of reference determinant is: ",KPnt(1),KPnt(2),KPnt(3)
        ENDIF
        IF(tFixLz) THEN
            CALL GetLz(HFDet,NEl,HFLz)
            WRITE(6,"(A,I5)") "Ml value of reference determinant is: ",HFLz
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
                WRITE(6,*) "HFDet: ",HFDet(:)
                WRITE(6,*) "HFDetTest: ",HFDetTest(:)
                CALL Stop_All(this_routine,"HF Determinant incorrectly decoded.")
            ENDIF
        enddo
        CALL LargestBitSet(iLutHF,NIfD,LargestOrb)
        IF(LargestOrb.ne.HFDet(NEl)) THEN
            CALL Stop_All(this_routine,"LargestBitSet FAIL")
        ENDIF
        nBits = CountBits(iLutHF, NIfD, NEl)
        IF(nBits.ne.NEl) THEN
            CALL Stop_All(this_routine,"CountBits FAIL")
        ENDIF

        IF(tCheckHighestPop) THEN
            ALLOCATE(HighestPopDet(0:NIfDBO),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for HighestPopDet")
            HighestPopDet(:)=0
        ENDIF

!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess=.true.
        tFoundOrbs(:)=.false.

        IF((.not.tHub).and.(.not.tUEG).and.TwoCycleSymGens) THEN
            do i=1,nSymLabels
!                WRITE(6,*) "NSymLabels: ",NSymLabels,i-1
                EndSymState=SymLabelCounts(1,i)+SymLabelCounts(2,i)-1
!                WRITE(6,*) "Number of states: ",SymLabelCounts(2,i)
                do j=SymLabelCounts(1,i),EndSymState

                    Beta=(2*SymLabelList(j))-1
                    Alpha=(2*SymLabelList(j))
                    SymAlpha=INT((G1(Alpha)%Sym%S),4)
                    SymBeta=INT((G1(Beta)%Sym%S),4)
!                    WRITE(6,*) "***",Alpha,Beta

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
                    WRITE(6,*) "Orbital: ",i, " not found."
                    CALL Stop_All("SetupParameters","Orbital not found")
                ENDIF
            enddo
        ENDIF
        IF(.not.tSuccess) THEN
            WRITE(6,*) "************************************************"
            WRITE(6,*) "**                 WARNING!!!                 **"
            WRITE(6,*) "************************************************"
            WRITE(6,*) "Symmetry information of orbitals not the same in alpha and beta pairs."
            WRITE(6,*) "Symmetry now set up in terms of spin orbitals"
            WRITE(6,*) "I strongly suggest you check that the reference energy is correct."
        ENDIF
! From now on, the orbitals are contained in symlabellist2 and symlabelcounts2 rather than the original arrays.
! These are stored using spin orbitals.
!Assume that if we want to use the non-uniform random excitation generator, we also want to use the NoSpinSym full excitation generators if they are needed. 

        CALL GetSym(HFDet,NEl,G1,NBasisMax,HFSym)
        
        WRITE(6,"(A,I10)") "Symmetry of reference determinant is: ",INT(HFSym%Sym%S,4)
        SymHF=0
        do i=1,NEl
            SymHF=IEOR(SymHF,INT(G1(HFDet(i))%Sym%S,4))
        enddo
        WRITE(6,"(A,I10)") "Symmetry of reference determinant from spin orbital symmetry info is: ",SymHF
        if(SymHF.ne.INT(HFSym%Sym%S,4)) then
            !When is this allowed to happen?! Comment!!
            call warning(this_routine,"Inconsistency in the symmetry arrays. Beware.")
        endif

!If using a CAS space truncation, write out this CAS space
        IF(tTruncCAS) THEN
            IF(tTruncInitiator) THEN
                WRITE(6,'(A)') " *********** INITIATOR METHOD IN USE ***********"
                WRITE(6,'(A)') " Fixed initiator space defined using the CAS method."
            ELSE
                WRITE(6,*) "Truncated CAS space detected. Writing out CAS space..."
            ENDIF
            WRITE(6,'(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (", &
                (OccCASOrbs+VirtCASOrbs)/2,",",OccCASOrbs,")"
            DO I=NEl-OccCASorbs+1,NEl
                WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
                WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
            WRITE(6,'(A)') " -------------------------------------------------------------------------------------------------"
            DO I=NEl+1,NEl+VirtCASOrbs
                WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
                WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
        ELSEIF(tTruncInitiator) THEN
            WRITE(6,'(A)') " *********** INITIATOR METHOD IN USE ***********"
            WRITE(6,'(A /)') " Starting with only the HF determinant in the fixed initiator space."
        ENDIF

        if(tSpawn_Only_Init.and.tSpawn_Only_Init_Grow) then
            write(6,"(A)") "Only allowing initiator determinants to spawn during growth phase, but allowing all to "&
                & //"spawn after fixed shift."
        elseif(tSpawn_Only_Init) then
            write(6,"(A)") "Only allowing initiator determinants to spawn"
        endif
 
!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        IF(tUEG.or.tHub.or.tNoSingExcits) THEN
            exflag=2
        ELSE
            exflag=3
        ENDIF
        IF(.not.tKPntSym) THEN
!Count all possible excitations - put into HFConn
!TODO: Get CountExcitations3 working with tKPntSym
            CALL CountExcitations3(HFDet,exflag,nSingles,nDoubles)
        ELSE
            !use Alex's old excitation generators to enumerate all excitations.
            CALL CountExcitsOld(HFDet,exflag,nSingles,nDoubles)
        ENDIF
        HFConn=nSingles+nDoubles

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        if(.not.tRestart) then
            Seed=abs(G_VMC_Seed-iProcIndex)
            WRITE(6,*) "Value for seed is: ",Seed
            !Initialise...
            CALL dSFMT_init(Seed)
        else
            !Reset the DiagSft to its original value
            DiagSft = InputDiagSft
        endif
        
        ! Option tRandomiseHashOrbs has now been removed.
        ! Its behaviour is now considered default
        ! --> Create a random mapping for the orbitals 
        ALLOCATE(RandomHash(nBasis),stat=ierr)
        IF(ierr.ne.0) THEN
            CALL Stop_All(this_routine,"Error in allocating RandomHash")
        ENDIF
        RandomHash(:)=0
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

!                WRITE(6,*) "Random Orbital Indexing for hash:"
!                WRITE(6,*) RandomHash(:)
            if (tSpatialOnlyHash) then
                step = 2
            else
                step = 1
            endif
            do i=1,nBasis
                IF((RandomHash(i).eq.0).or.(RandomHash(i).gt.nBasis*1000)) THEN
                    CALL Stop_All(this_routine,"Random Hash incorrectly calculated")
                ENDIF
                do j = i+step, nBasis, step
                    IF(RandomHash(i).eq.RandomHash(j)) THEN
                        CALL Stop_All(this_routine,"Random Hash incorrectly calculated")
                    ENDIF
                enddo
            enddo
        ENDIF
        !Now broadcast to all processors
        CALL MPIBCast(RandomHash,nBasis)

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
        WRITE(6,*) "Reference Energy set to: ",Hii
        if(tUEG) then
            !We require calculation of the sum of fock eigenvalues,
            !without knowing them - calculate from the full 1e matrix elements
            !of full hamiltonian removing two electron terms.
            TempHii=GetH0Element4(HFDet,HFDet)
        else
            !Know fock eigenvalues, so just use these.
            TempHii=GetH0Element3(HFDet)
        endif
        Fii=REAL(TempHii,dp)

!Find the highest energy determinant...
        IF(.not.tSpn) THEN
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
            WRITE(6,"(A,G25.15)") "Highest energy determinant is (approximately): ",REAL(TempHii,dp)
            WRITE(6,"(A,F25.15)") "This means tau should be no more than about ",-2.D0/REAL(TempHii,dp)
!            WRITE(6,*) "Highest energy determinant is: ", HighEDet(:)
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
                WRITE(6,*) "High spin calculation with single excitations also used to calculate energy."
            ELSEIF(tUHF) THEN
                WRITE(6,*) "High spin calculation - but single excitations will *NOT* be used to calculate energy as "&
                & //"this is an unrestricted calculation."
            ELSE
                CALL Stop_All("SetupParameters","High-spin, restricted calculation detected, but single excitations are "&
                & //"not being used to calculate the energy.  &
                  & Either use the UHF keyword, or turn off brillouins theorem using NOBRILLOUINS, ROHF or ROTATEDORBS.")
            ENDIF
!            tRotatedOrbs=.true.
!        ELSEIF(LMS.ne.0) THEN
!            CALL Stop_All(this_routine,"Ms not equal to zero, but tSpn is false. Error here")
        ENDIF

!Initialise variables for calculation on each node
        iter=0          !This is set so that calls to CalcParentFlag in the initialisation are ok with the logging.
        iPopsTimers=1   !Number of timed popsfiles written out
        IterTime=0.0
        ProjectionE=0.D0
        AvSign=0.D0
        AvSignHFD=0.D0
        SumENum=0.D0
        SumNoatHF=0
        NoatHF=0
        NoatDoubs=0
        Annihilated=0
        Acceptances=0
        PreviousCycles=0
        NoBorn=0
        SpawnFromSing=0
        NoDied=0
        HFCyc=0
        ENumCyc=0.D0
        VaryShiftCycles=0
        AvDiagSft=0.D0
        SumDiagSft=0.D0
        SumDiagSftAbort=0.D0
        AvDiagSftAbort=0.D0
        NoAborted=0
        NoAddedInitiators=0
        NoInitDets=0
        NoNonInitDets=0
        NoInitWalk=0
        NoNonInitWalk=0
        NoExtraInitDoubs=0
        InitRemoved=0
        TotImagTime=0.D0
        DiagSftRe=0.D0
        DiagSftIm=0.D0
        sum_proje_denominator = 0
        cyc_proje_denominator = 0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllNoatHF=0
        AllNoatDoubs=0
        AllSumNoatHF = 0
        AllGrowRate=0.D0
        AllGrowRateAbort=0.D0
!        AllMeanExcitLevel=0.D0
        AllSumWalkersCyc=0
        AllAvSign=0.D0
        AllAvSignHFD=0.D0
        AllNoBorn=0
        AllSpawnFromSing=0
        AllNoDied=0
        AllAnnihilated=0
        AllENumCyc=0.D0
        AllHFCyc=0.D0
!        AllDetsNorm=0.D0
        AllNoAborted=0
        AllNoAddedInitiators=0
        AllNoInitDets=0
        AllNoNonInitDets=0
        AllNoInitWalk=0
        AllNoNonInitWalk=0
        AllNoExtraInitDoubs=0
        AllInitRemoved=0

        ! Initialise the fciqmc counters
        iter_data_fciqmc%update_growth = 0
        iter_data_fciqmc%update_iters = 0
 

        IF(tHistSpawn.or.(tCalcFCIMCPsi.and.tFCIMC).or.tHistHamil) THEN
            ALLOCATE(HistMinInd(NEl))
            ALLOCATE(HistMinInd2(NEl))
            maxdet=0
            do i=1,nel
                maxdet=maxdet+2**(nbasis-i)
            enddo

            IF(.not.allocated(FCIDets)) THEN
                CALL Stop_All(this_routine,"A Full Diagonalization is required before histogramming can occur.")
            ENDIF

            IF(tHistHamil) THEN
                WRITE(6,*) "Histogramming total Hamiltonian, with Dets=", Det
                ALLOCATE(HistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays") 
                HistHamil(:,:)=0.D0
                ALLOCATE(AvHistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                AvHistHamil(:,:)=0.D0
                IF(iProcIndex.eq.0) THEN
                    ALLOCATE(AllHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                    AllHistHamil(:,:)=0.D0
                    ALLOCATE(AllAvHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                    AllAvHistHamil(:,:)=0.D0
                ENDIF
            ELSE
                WRITE(6,*) "Histogramming spawning wavevector, with Dets=", Det
                ALLOCATE(Histogram(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays ")
                ENDIF
                Histogram(:,:)=0.D0
                ALLOCATE(AllHistogram(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
            ENDIF
            IF(tHistSpawn) THEN
                ALLOCATE(InstHist(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
                InstHist(:,:)=0.D0
                ALLOCATE(AvAnnihil(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
                AvAnnihil(:,:)=0.D0
                ALLOCATE(InstAnnihil(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
                InstAnnihil(:,:)=0.D0
            ENDIF

            IF(iProcIndex.eq.0) THEN
                IF(tHistSpawn) THEN
                    ALLOCATE(AllInstHist(1:lenof_sign,1:det),stat=ierr)
                    ALLOCATE(AllInstAnnihil(1:lenof_sign,1:det),stat=ierr)
                    ALLOCATE(AllAvAnnihil(1:lenof_sign,1:det),stat=ierr)
                ENDIF
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays")
                ENDIF
            ENDIF
        ELSEIF(tHistEnergies) THEN
            WRITE(6,*) "Histogramming the energies of the particles, with iNoBins=",iNoBins, " and BinRange=", BinRange
            WRITE(6,*) "Histogramming spawning events from ",-OffDiagMax, " with BinRange = ", OffDiagBinRange
            iOffDiagNoBins=INT((2.D0*OffDiagMax)/OffDiagBinRange)+1
            WRITE(6,*) "This gives ",iOffDiagNoBins," bins to histogram the off-diagonal matrix elements."
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
            HistogramEnergy(:)=0.D0
            AttemptHist(:)=0.D0
            SpawnHist(:)=0.D0
            SinglesHist(:)=0.D0
            SinglesAttemptHist(:)=0.D0
            SinglesHistOccOcc(:)=0.D0
            SinglesHistOccVirt(:)=0.D0
            SinglesHistVirtOcc(:)=0.D0
            SinglesHistVirtVirt(:)=0.D0
            DoublesHist(:)=0.D0
            DoublesAttemptHist(:)=0.D0
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


        ! Initialise the spin distribution histogramming
        if (tHistSpinDist) then
            call init_hist_spin_dist ()
        endif

!Need to declare a new MPI type to deal with the long integers we use in the hashing, and when reading in from POPSFILEs
!        CALL MPI_Type_create_f90_integer(18,mpilongintegertype,error)
!        CALL MPI_Type_commit(mpilongintegertype,error)
        IF(tUseBrillouin) THEN
            WRITE(6,*) "Brillouin theorem specified, but this will not be in use with the non-uniform excitation generators."
        ENDIF
        WRITE(6,*) "Non-uniform excitation generators in use."
        CALL CalcApproxpDoubles()
        IF(TauFactor.ne.0.D0) THEN
            WRITE(6,*) "TauFactor detected. Resetting Tau based on connectivity of: ",HFConn
            Tau=TauFactor/REAL(HFConn,dp)
            WRITE(6,*) "Tau set to: ",Tau
        ENDIF
        IF(StepsSftImag.ne.0.D0) THEN
            WRITE(6,*) "StepsShiftImag detected. Resetting StepsShift."
            StepsSft=NINT(StepsSftImag/Tau)
            IF(StepsSft.eq.0) StepsSft=1
            WRITE(6,*) "StepsShift set to: ",StepsSft
        ENDIF
!        IF(tConstructNOs) THEN
!! This is the option for constructing the natural orbitals actually during a NECI calculation.  This is different (and probably a lot more complicated and doesn't 
!! currently work) from the FINDCINATORBS option which finds the natural orbitals given a final wavefunction.
!            ALLOCATE(OneRDM(nBasis,nBasis),stat=ierr)
!            CALL LogMemAlloc('OneRDM',nBasis*nBasis,8,this_routine,OneRDMTag,ierr)
!            OneRDM(:,:)=0.D0
!        ENDIF

        IF(tCCMC) then      !Once alex's CCMC/FCIMC unification project is finished we can remove this
                            !The problem is that ValidSpawnedList is now setup in InitFCIMCCalcPar.
                            !This is for compatibility with POPSFILE v.3, where MaxSpawned is calculated
                            !from the number in the POPSFILE.
                            !Set it up here for CCMC.
            Call SetupValidSpawned(InitWalkers)
        endif

        IF(TPopsFile) THEN
            IF(mod(iWritePopsEvery,StepsSft).ne.0) then
                CALL Warning(this_routine,"POPSFILE writeout should be a multiple of the update cycle length.")
            endif
        ENDIF

        if (TReadPops) then
            if (tStartSinglePart .and. .not. tReadPopsRestart) then
                call warning(this_routine, &
                               "ReadPOPS cannot work with StartSinglePart: ignoring StartSinglePart")
                tStartSinglePart = .false.
            end if
        endif

        IF(.not.TReadPops) THEN
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
        ENDIF
        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
        WRITE(6,*) "Maximum connectivity of HF determinant is: ",HFConn
        IF(TStartSinglePart) THEN
            TSinglePartPhase=.true.
        ELSE
            TSinglePartPhase=.false.
        ENDIF
        
        IF(ICILevel.ne.0) THEN
!We are truncating the excitations at a certain value
            TTruncSpace=.true.
            WRITE(6,'(A,I4)') "Truncating the S.D. space at determinants will an excitation level w.r.t. HF of: ",ICILevel
        ENDIF
        IF(tTruncCAS.or.tStartCAS) THEN
!We are truncating the FCI space by only allowing excitations in a predetermined CAS space.
!The following line has already been written out if we are doing a CAS calculation.
!            WRITE(6,'(A,I4,A,I5)') "Truncating the S.D. space as determinants must be within a CAS of ",OccCASOrbs," , ",VirtCASOrbs
!The SpinInvBRR array is required for the tTruncCAS option. Its properties are explained more fully in the subroutine. 

            CALL CreateSpinInvBRR()

            CASmax=NEl+VirtCASorbs
! CASmax is the max spin orbital number (when ordered energetically) within the chosen active space.
! Spin orbitals with energies larger than this maximum value must be unoccupied for the determinant
! to be in the active space.
            CASmin=NEl-OccCASorbs
! CASmin is the max spin orbital number below the active space.  As well as the above criteria, spin 
! orbitals with energies equal to, or below that of the CASmin orbital must be completely occupied for 
! the determinant to be in the active space.

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
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the lowest ',NPartFrozen,' spin orbitals so that no more than ', &
                NHolesFrozen,' holes exist within this core.'
            CALL CreateSpinInvBRR()
        ENDIF
        IF(tPartFreezeVirt) THEN
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the highest ',NVirtPartFrozen, &
                ' virtual spin orbitals so that no more than ',NElVirtFrozen,' electrons occupy these orbitals.'
            CALL CreateSpinInvBRR()
        ENDIF

        if (tTruncNOpen) then
            write(6, '("Truncating determinant space at a maximum of ",i3," &
			           &unpaired electrons.")') trunc_nopen_max
        endif

        SymFactor=(Choose(NEl,2)*Choose(nBasis-NEl,2))/(HFConn+0.D0)
        TotDets=1.D0
        do i=1,NEl
            WRITE(6,"(A,I5,I20)") "Approximate excitation level population: ",i,NINT((Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor)
            TotDets=TotDets+(Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor
        enddo
        WRITE(6,"(A,I20)") "Approximate size of determinant space is: ",NINT(TotDets)

    END SUBROUTINE SetupParameters

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
            WRITE(6,*) "Turning Brillouins theorem off since we are using non-canonical complex orbitals"
            tNoBrillouin=.true.
        ENDIF

!Special case of defining a det with LMS=0, but which is open shell. No Brillouins if it's a restricted HF calc.
        tSpinPair = .false.
        IF(tDefineDet.and.(LMS.eq.0).and.(.not.tUHF)) THEN
!If we are defining our own reference determinant, we want to find out if it is open shell or closed to know whether or not brillouins theorem holds.            
!If LMS/=0, then it is easy and must be open shell, otherwise we need to consider the occupied orbitals.
            do i=1,(NEl-1),2
!Assuming things will probably go alpha beta alpha beta, run through each alpha and see if there's a corresponding beta.
                tSpinPair=.false.
                IF(MOD(BRR(FDet(i)),2).ne.0) THEN
!Odd energy, alpha orbital.                    
                    IF(BRR(FDet(i+1)).ne.(BRR(FDet(i))+1)) THEN
!Check the next orbital to see if it's the beta (will be alpha+1 when ordered by energy). 
!If not, check the other orbitals for the beta, as it's possible the orbitals are ordered weird (?).
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
                WRITE(6,'(A)') " Using an open shell reference determinant in a basis of restricted HF orbitals; " &
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

!        CASmax=NEl+VirtCASorbs
! CASmax is the max spin orbital number (when ordered energetically) within the chosen active space.
! Spin orbitals with energies larger than this maximum value must be unoccupied for the determinant
! to be in the active space.
!        CASmin=NEl-OccCASorbs   (These have been moved to the InitCalc subroutine so they're not calculated
! each time.
! CASmin is the max spin orbital number below the active space.  As well as the above criteria, spin 
! orbitals with energies equal to, or below that of the CASmin orbital must be completely occupied for 
! the determinant to be in the active space.

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
! if at any stage an electron has an energy greater than the CASmax value, the determinant can be ruled out
! of the active space.  Upon identifying this, it is not necessary to check the remaining electrons.
            else
                if (SpinInvBRR(orb).le.CASmin) THEN
                    z=z+1
                endif
! while running over all electrons, the number that occupy orbitals equal to or below the CASmin cutoff are
! counted.
            endif
        enddo

        if(tElecInVirt.or.(z.ne.CASmin)) THEN
! if an electron is in an orbital above the active space, or the inactive orbitals are not full, the determinant is automatically ruled out.            
            TestifDETinCAS=.false.
        else
! no orbital in virtual and the inactive orbitals are completely full - det in active space.        
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
            nopen = count_open_orbs (ilutnJ)
            if (nopen > trunc_nopen_max) &
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



!This is the same as BinSearchParts1, but this time, the list to search is passed in as an argument. The list goes from 1 to Length, but only between MinInd and MaxInd is actually searched.
    SUBROUTINE BinSearchParts3(iLut,List,Length,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER :: Length
        INTEGER(KIND=n_int) :: iLut(0:NIfTot), List(0:NIfTot,Length)
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
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
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(List(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
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

    END SUBROUTINE BinSearchParts3
    
    SUBROUTINE CalcApproxpDoubles()
        use SystemData , only : tAssumeSizeExcitgen,tUseBrillouin,tNoSingExcits
        use CalcData , only : SinglesBias
        use SymData , only : SymClassSize
        use SymExcit3 , only : CountExcitations3
        INTEGER :: iTotal
        integer :: nSing, nDoub, ncsf, excitcount, ierr, iExcit
        integer :: nStore(6), iMaxExcit, nExcitMemLen, nJ(nel)
        integer, allocatable :: EXCITGEN(:)
        character(*), parameter :: this_routine = 'CalcApproxpDoubles'
        logical :: TempUseBrill

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
                WRITE(6,*) "Since we are using a real-space hubbard model, only single excitations are connected &
                &   and will be generated."
                pDoubles=0.D0
                RETURN
            ELSE
                WRITE(6,*) "Since we are using a momentum-space hubbard model/UEG, only double excitaitons &
     &                          are connected and will be generated."
                pDoubles=1.D0
                RETURN
            ENDIF
        elseif(tNoSingExcits) then
            write(6,*) "Only double excitations will be generated"
            return
        ENDIF

!NSing=Number singles from HF, nDoub=No Doubles from HF

        WRITE(6,"(A)") " Calculating approximate pDoubles for use with excitation generator by looking a excitations from HF."
        exflag=3
        IF(tKPntSym) THEN
            !use Alex's old excitation generators.
            !However, we have to ensure that brillouins theorem isn't on!
            IF(tUseBrillouin) THEN
                TempUseBrill=.true.
                tUseBrillouin=.false.
            ELSE
                TempUseBrill=.false.
            ENDIF

            iMaxExcit=0
            nStore(1:6)=0
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
            ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating excitation generator")
            EXCITGEN(:)=0
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
        !    CALL GetSymExcitCount(EXCITGEN,DetConn)
            excitcount=0

        lp2: do while(.true.)
                CALL GenSymExcitIt2(HFDet,nEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
                IF(nJ(1).eq.0) exit lp2
                IF(iExcit.eq.1) THEN
                    nSing=nSing+1
                ELSEIF(iExcit.eq.2) THEN
                    nDoub=nDoub+1
                ELSE
                    CALL Stop_All(this_routine,"Trying to generate more than doubles!")
                ENDIF
            enddo lp2
            tUseBrillouin=TempUseBrill
        ELSE
            CALL CountExcitations3(HFDet,exflag,nSing,nDoub)
        ENDIF
        iTotal=nSing + nDoub + ncsf

        WRITE(6,"(I7,A,I7,A)") NDoub, " double excitations, and ",NSing, &
            " single excitations found from HF. This will be used to calculate pDoubles."

        IF(SinglesBias.ne.1.D0) THEN
            WRITE(6,*) "Singles Bias detected. Multiplying single excitation connectivity of HF determinant by ", &
                SinglesBias," to determine pDoubles."
        ENDIF

        IF((NSing+nDoub+ncsf).ne.iTotal) THEN
            CALL Stop_All("CalcApproxpDoubles","Sum of number of singles and number of doubles does " &
            & //"not equal total number of excitations")
        ENDIF
        IF((NSing.eq.0).or.(NDoub.eq.0)) THEN
            WRITE(6,*) "Number of singles or doubles found equals zero. pDoubles will be set to 0.95. Is this correct?"
            pDoubles = 0.95
            pSingles = 0.05
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

        IF(SinglesBias.ne.1.D0) THEN
            write (6, '("pDoubles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pDoubles, real(nDoub,dp) / real(iTotal,dp)
            write (6, '("pSingles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pSingles, real(nSing,dp) / real(iTotal,dp)

            WRITE(6,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ", &
                & real(nDoub,dp)/real(iTotal,dp)
        ELSE
            write (6,'(A,F14.6)') " pDoubles set to: ", pDoubles
            write (6,'(A,F14.6)') " pSingles set to: ", pSingles
        ENDIF

        WRITE(6,'(A,F15.10)') " Assuming an average K_ij magnitude of approx 0.01, an appropriate tau is " &
            & //"predicted to be around: ",(0.02*(1.D0/(REAL(NSing)+REAL(NDoub))))/0.01
!This is a rough guesstimate of what tau might like to be, assuming K_ij is approx 0.01 on average, and we want a probability of spawning to be about 0.02.        
!These are just stats taken from one system... will investigate further...

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
!            WRITE(6,*) "================================"
!            WRITE(6,*) "BRR is "
!            WRITE(6,*) BRR(:)
!        ENDIF
        
        SpinInvBRR(:)=0
        
        t=0
        DO I=1,NBASIS
            t=t+1
            SpinInvBRR(BRR(I))=t
        ENDDO

!        IF(iProcIndex.eq.root) THEN
!            WRITE(6,*) "================================"
!            WRITE(6,*) "SpinInvBRR is "
!            WRITE(6,*) SpinInvBRR(:)
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
            WRITE(6,*) VecSlot,NoDoubs
            CALL Stop_All("StoreDoubs","Problem enumerating all double excitations")
        ENDIF

    END SUBROUTINE StoreDoubs

!Initialize the Histogramming searching arrays if necessary
    SUBROUTINE InitHistMin()
        IF(tHistSpawn.or.tCalcFCIMCPsi.and.(Iter.ge.NHistEquilSteps)) THEN
            IF(Iter.eq.NHistEquilSteps) THEN
                IF(iProcIndex.eq.Root) WRITE(6,*) 'The iteration is equal to HISTEQUILSTEPS.  Beginning to histogram.'
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
        integer :: nopen, nup, pos, sgn(lenof_sign), nfound
        real(dp) :: norm
        character(*), parameter :: t_r = 'setup_linear_comb'

        write(6,*) 'Initialising projection onto linear combination of &
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
            call MPISumAll_inplace (proje_ref_coeffs)
            norm = sqrt(sum(proje_ref_coeffs**2))
            if (norm == 0) norm = 1
            proje_ref_coeffs = proje_ref_coeffs / norm

        endif

    end subroutine setup_linear_comb

    subroutine update_linear_comb_coeffs ()

        integer :: i, pos, sgn(lenof_sign)
        real(dp) :: norm

        if (nproje_sum > 1) then

            proje_ref_coeffs = 0
            do i = 1, nproje_sum

                pos = binary_search (CurrentDets(:,1:TotWalkers), &
                                     proje_ref_iluts(:,i), NIfD+1)
                if (pos > 0) then
                    call extract_sign (CurrentDets(:,pos), sgn)
                    proje_ref_coeffs(i) = sgn(1)
                endif
            enddo

            call MPISumAll_inplace (proje_ref_coeffs)
            norm = sqrt(sum(proje_ref_coeffs**2))
            if (norm == 0) norm = 1
            proje_ref_coeffs = proje_ref_coeffs / norm
            
        endif

    end subroutine


    subroutine clean_linear_comb ()

        if (allocated(proje_ref_dets)) &
            deallocate(proje_ref_dets)
        if (allocated(proje_ref_iluts)) &
            deallocate(proje_ref_iluts)
        if (allocated(proje_ref_coeffs)) &
            deallocate(proje_ref_coeffs)

    end subroutine clean_linear_comb

!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
!AJWT added optional argument dProbFin which is a probability that whatever gave this contribution was generated.
!  It defaults to 1, and weights the contribution of this det (only in the projected energy) by dividing its contribution by this number 
    subroutine SumEContrib (nI, ExcitLevel, WSign, ilut, HDiagCurr, dProbFin)

        integer, intent(in) :: nI(nel), ExcitLevel
        integer, intent(in) :: wSign(lenof_sign)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: HDiagCurr, dProbFin

        integer :: i, bin, pos, ExcitLevel_local, ExcitLevelSpinCoup
        integer :: PartInd, OpenOrbs, spatial_ic
        integer(n_int) :: iLutSym(0:NIfTot)
        logical tSuccess
        integer :: iUEG1, iUEG2
        HElement_t :: HOffDiag
        HElement_t :: HDoubDiag
        integer :: DoubEx(2,2),DoubEx2(2,2),kDoub(3) ! For histogramming UEG doubles
        logical :: tDoubParity,tDoubParity2 ! As above

        ! Are we performing a linear sum over various determinants?
        ! TODO: If we use this, function pointer it.
        HOffDiag = 0
        if (proje_linear_comb .and. nproje_sum > 1) then

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

                if (iter > NEquilSteps) SumNoatHF = SumNoatHF + wSign
                NoatHF = NoatHF + wSign
                ! Number at HF * sign over course of update cycle
                HFCyc = HFCyc + wSign

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
                
                if (ExcitLevel == 2) NoatDoubs = NoatDoubs + sum(abs(wSign))

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
        if (iter > NEquilSteps) &
            SumENum = SumENum + (HOffDiag * ARR_RE_OR_CPLX(wSign)) / dProbFin
        ENumCyc = ENumCyc + (HOffDiag * ARR_RE_OR_CPLX(wSign)) / dProbFin

        ! -----------------------------------
        ! HISTOGRAMMING
        ! -----------------------------------

        if ((tHistSpawn .or. (tCalcFCIMCPsi .and. tFCIMC)) .and. &
            (iter >= NHistEquilSteps)) then
            ! Histogram particles by determinant
            call add_hist_spawn (ilut, wSign, ExcitLevel_local, dProbFin)
        elseif (tHistEnergies) then
            ! Histogram particles by energy
            call add_hist_energies (ilut, wSign, HDiagCurr)
        endif

        ! Are we doing a spin-projection histogram?
        if (tHistSpinDist) &
            call test_add_hist_spin_dist_det (ilut, wSign)

        ! Maintain a list of the degree of occupation of each orbital
        if (tPrintOrbOcc .and. (iter >= StartPrintOrbOcc)) then
            if ((tPrintOrbOccInit .and. test_flag(ilut,flag_is_initiator(1)))&
                .or. .not. tPrintOrbOccInit) then
                forall (i = 1:nel) OrbOccs(nI(i)) = OrbOccs(nI(i)) &
                                          + (real(wSign(1)) * real(wSign(1)))
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
                    iUEG1=DoubsUEGLookup(DoubEx(1,1))
                    iUEG2=DoubsUEGLookup(DoubEx(1,2))
                    if (iUEG1.eq.0.or.iUEG2.eq.0) call stop_all("SumEContrib","Array bounds issue")
                    DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),1)=DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),1)+REAL(WSign(1))
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
                        DoubsUEG(iUEG1,iUEG2,DoubEx(2,1),4)=REAL(kDoub(1))**2.D0+REAL(kDoub(2))**2.D0+REAL(kDoub(3))**2.D0
                    endif
                endif
            endif
        endif

    end subroutine SumEContrib


!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        use FciMCLoggingMOD , only : InitHistInitPops
        use SystemData , only : tRotateOrbs
        use CalcData , only : InitialPart,tstartmp1,tStartCAS
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil
        use constants , only : size_n_int
        use DeterminantData , only : write_det
        use nElRDMMod, only: InitRDM
        INTEGER :: ierr,iunithead
        LOGICAL :: formpops,binpops
        INTEGER :: error,MemoryAlloc,PopsVersion,WalkerListSize,j,iLookup
        INTEGER, DIMENSION(lenof_sign) :: InitialSign
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'
        integer :: ReadBatch    !This parameter determines the length of the array to batch read in walkers from a popsfile
        real(8) :: Gap
        !Variables from popsfile header...
        logical :: tPop64Bit,tPopHPHF,tPopLz
        integer :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(8) :: iPopAllTotWalkers
        real(8) :: PopDiagSft
        integer(8) , dimension(lenof_sign) :: PopSumNoatHF
        HElement_t :: PopAllSumENum

        if(tReadPops.and..not.tPopsAlreadyRead) then
            call open_pops_head(iunithead,formpops,binpops)
            PopsVersion=FindPopsfileVersion(iunithead)
            if(iProcIndex.eq.root) close(iunithead)
            write(6,*) "POPSFILE VERSION ",PopsVersion," detected."
        endif

        if(tPopsMapping.and.(PopsVersion.lt.3)) then
            write(6,*) "Popsfile mapping cannot work with old POPSFILEs"
            call stop_all("InitFCIMCCalcPar","Popsfile mapping cannot work with old POPSFILEs")
        endif

        ! Initialise measurement of norm, to avoid divide by zero
        norm_psi = 1.0_dp

        if (tReadPops .and. (PopsVersion.lt.3) .and..not.tPopsAlreadyRead) then
!Read in particles from multiple POPSFILES for each processor
            !Ugh - need to set up ValidSpawnedList here too...
            call SetupValidSpawned(InitWalkers)
            WRITE(6,*) "Reading in initial particle configuration from *OLD* POPSFILES..."
            CALL ReadFromPopsFilePar()
        ELSE
!initialise the particle positions - start at HF with positive sign
!Set the maximum number of walkers allowed
            if(tReadPops.and..not.tPopsAlreadyRead) then
                !We must have a v.3 popsfile. Read header.
                call open_pops_head(iunithead,formpops,binpops)
                call ReadPopsHeadv3(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                        iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                        PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)

                call CheckPopsParams(tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                        iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                        PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,WalkerListSize)

                if(iProcIndex.eq.root) close(iunithead)
            else
                WalkerListSize=InitWalkers
            endif

            MaxWalkersPart=NINT(MemoryFacPart*WalkerListSize)
            WRITE(6,"(A,I14)") " Memory allocated for a maximum particle number per node of: ",MaxWalkersPart
            Call SetupValidSpawned(WalkerListSize)

!Put a barrier here so all processes synchronise
            CALL MPIBarrier(error)
!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
            MemoryAlloc=(NIfTot+1)*MaxWalkersPart*size_n_int    !Memory Allocated in bytes

            IF(.not.tRegenDiagHEls) THEN
                ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                WalkVecH(:)=0.d0
                MemoryAlloc=MemoryAlloc+8*MaxWalkersPart
            ELSE
                WRITE(6,"(A,F14.6,A)") " Diagonal H-Elements will not be stored. This will *save* ", &
                    & REAL(MaxWalkersPart*8,dp)/1048576.D0," Mb/Processor"
            ENDIF
            
            WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned, &
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

!Allocate pointers to the correct walker arrays
            CurrentDets=>WalkVecDets
            IF(.not.tRegenDiagHEls) THEN
                CurrentH=>WalkVecH
            ENDIF
        
            ! Get the (0-based) processor index for the HF det.
            iHFProc = DetermineDetNode(HFDet,0)
            WRITE(6,*) "Reference processor is: ",iHFProc
            write(6,*) "Initial reference is: "
            call write_det(6,HFDet,.true.)

            TotParts(:)=0
            TotPartsOld(:)=0
            NoatHF=0

            if (tSpawnSpatialInit) then
                max_inits = int(MemoryFacInit * INitWalkers)
                no_spatial_init_dets = 0
                allocate(CurrentInits(0:nIfTot, max_inits), stat=ierr)
                call LogMemAlloc('CurrentInits', max_inits * (NIfTot+1), &
                                 size_n_int, this_routine, CurrentInitTag, &
                                 ierr)
            endif

!If we have a popsfile, read the walkers in now.
            if(tReadPops.and..not.tPopsAlreadyRead) then

                ReadBatch=MaxSpawned    !ReadBatch is the number of walkers to read in from the popsfile at one time.
                                        !The larger it is, the fewer communications will be needed to scatter the particles.
                                        !By default, the new array (which is only created on the root processors) is the
                                        !same length as the spawning arrays.

                !TotWalkers and TotParts are returned as the dets and parts on each processor.
                call ReadFromPopsfilev3(iPopAllTotWalkers,ReadBatch,TotWalkers,TotParts,NoatHF,CurrentDets,MaxWalkersPart)

                !Setup global variables
                TotWalkersOld=TotWalkers
                TotPartsOld = TotParts
                call MPISumAll(TotWalkers,AllTotWalkers)
                AllTotWalkersOld = AllTotWalkers
                call MPISumAll(TotParts,AllTotParts)
                AllTotPartsOld=AllTotParts
                call MPISumAll(NoatHF,AllNoatHF)
                OldAllNoatHF=AllNoatHF
                AllNoAbortedOld=0.D0
                iter_data_fciqmc%tot_parts_old = AllTotParts

                ! Calculate the projected energy for this iteration.
                if (any(AllSumNoatHF /= 0)) &
                    ProjectionE = AllSumENum / ARR_RE_OR_CPLX(AllSumNoatHF)
                
                if(iProcIndex.eq.iHFProc) then
                    !Need to store SumENum and SumNoatHF, since the global variable All... gets wiped each iteration. 
                    !Rather than POPSFILE v2, where the average values were scattered, just store the previous
                    !energy contributions on the root node.
                    SumNoatHF=AllSumNoatHF
                    SumENum=AllSumENum

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
                        WRITE(6,"(A,I15,A,F9.3,A,I15)") " Initial number of particles set to ",InitialPart, &
                            " and shift will be held at ",DiagSft," until particle number gets to ",InitWalkers*nNodes
                    else
                        write(6,"(A,I16)") "Initial number of walkers per processor chosen to be: ", InitWalkers
                    endif


                    !Setup initial walker local variables for HF walkers start
                    IF(iProcIndex.eq.iHFProc) THEN

                        ! Encode the reference determinant identification.
                        call encode_det(CurrentDets(:,1), iLutHF)

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

                        ! HF energy is equal to 0 (by definition)
                        if (.not. tRegenDiagHEls) CurrentH(1) = 0

                        ! Obtain the initial sign
                        InitialSign = 0
                        if (tStartSinglePart) then
                            InitialSign(1) = InitialPart
                        else
                            InitialSign(1) = InitWalkers
                        endif
                        call encode_sign (CurrentDets(:,1), InitialSign)

                        ! set initial values for global control variables.
                        TotWalkers = 1
                        TotWalkersOld = 1
                        TotParts = InitialSign
                        TotPartsOld = InitialSign
                        NoatHF = InitialSign

                    ELSE
                        NoatHF = 0
                        TotWalkers = 0
                        TotWalkersOld = 0
                    ENDIF

                    OldAllNoatHF=0
                    AllNoatHF=0
                    IF(TStartSinglePart) THEN
        !Initialise global variables for calculation on the root node
                        IF(iProcIndex.eq.root) THEN
                            OldAllNoatHF(1)=InitialPart
                            AllNoatHF(1)=InitialPart
                            AllTotWalkers = 1
                            AllTotWalkersOld = 1
                            iter_data_fciqmc%tot_parts_old(1) = InitialPart
                            AllTotParts(1)=InitialPart
                            AllTotPartsOld(1)=InitialPart
                            AllNoAbortedOld=0.D0
                        ENDIF
                    ELSE
        !In this, only one processor has initial particles.
                        IF(iProcIndex.eq.Root) THEN
                            AllTotWalkers = 1
                            AllTotWalkersOld = 1
                            iter_data_fciqmc%tot_parts_old(1) = InitWalkers
                            AllTotParts(1)=InitWalkers
                            AllTotPartsOld(1)=InitWalkers
                            AllNoAbortedOld=0.D0
                        ENDIF
                    ENDIF

                endif   !tStartmp1
            endif  
        
            WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens + temp arrays) consists of : ", &
                & REAL(MemoryAlloc,dp)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Only one array of memory to store main particle list allocated..."
            WRITE(6,*) "Initial memory allocation sucessful..."
            CALL FLUSH(6)

        ENDIF   !End if initial walkers method

        ! If we are projecting onto a linear combination to calculate projE,
        ! then do the setup
        if (proje_linear_comb) &
            call setup_linear_comb ()

            
!Put a barrier here so all processes synchronise
        CALL MPIBarrier(error)

        IF(tTruncInitiator.or.tDelayTruncInit) THEN
            IF(tDelayTruncInit) tTruncInitiator=.false.
        ENDIF

        IF(tPrintOrbOcc) THEN
            ALLOCATE(OrbOccs(nBasis),stat=ierr)
            CALL LogMemAlloc('OrbOccs',nBasis,8,this_routine,OrbOccsTag,ierr)
            OrbOccs(:)=0.D0
        ENDIF

        IF(tPrintDoubsUEG) THEN
            ALLOCATE(DoubsUEG(NEl,NEl,nBasis,4),stat=ierr)
            DoubsUEG(:,:,:,:)=0.D0
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
        MaxInitPopPos=0
        MaxInitPopNeg=0

        IF(MaxNoatHF.eq.0) THEN
            MaxNoatHF=InitWalkers*nNodes
            HFPopThresh=MaxNoatHF
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

        if (tSpinProject) call init_yama_store ()

        !This keyword (tRDMonFly) is on from the beginning if we eventually plan to calculate the RDM's.
        !But the tFilling keywords don't become true until we actually start calculating them.
        IF(tRDMonFly) THEN
            !Initialises RDM stuff for both explicit and stochastic calculations of RDM.
            CALL InitRDM()
        ENDIF
        tFillingStochRDMonFly = .false.      
        tFillingExplicRDMonFly = .false.      
        !One of these becomes true when we have reached the relevant iteration to begin filling the RDM.

    end subroutine InitFCIMCCalcPar

!Routine to initialise the particle distribution according to a CAS diagonalisation. 
!This hopefully will help with close-lying excited states of the same sym.
    subroutine InitFCIMC_CAS()
        use SystemData, only : tSpn,tHPHFInts
        use CalcData, only: InitialPart
        use DeterminantData, only : write_det,write_det_len 
        use DetCalcData, only : NKRY,NBLK,B2L,nCycle
        use sym_mod , only : Getsym, writesym
        use MomInv, only: IsAllowedMI 
        type(BasisFN) :: CASSym
        integer :: i,j,ierr,nEval,NKRY1,NBLOCK,LSCR,LISCR,DetIndex,iNode,NoWalkers
        integer :: CASSpinBasisSize,elec,nCASDet,ICMax,GC,LenHamil,iInit,nHPHFCAS
        integer , allocatable :: CASBrr(:),CASDet(:),CASFullDets(:,:),nRow(:),Lab(:),ISCR(:),INDEX(:)
        integer , pointer :: CASDetList(:,:) => null()
        integer(n_int) :: iLutnJ(0:NIfTot)
        logical :: tMC
        HElement_t :: HDiagTemp
        real(dp) , allocatable :: CK(:,:),W(:),CKN(:,:),Hamil(:),A(:,:),V(:),BM(:),T(:),WT(:),SCR(:),WH(:),Work2(:),V2(:,:),AM(:)
        integer(TagIntType) :: ATag=0,VTag=0,BMTag=0,TTag=0,WTTag=0,SCRTag=0,WHTag=0,Work2Tag=0,V2Tag=0,ISCRTag=0,IndexTag=0,AMTag=0
        real(dp) :: CASRefEnergy,TotWeight,PartFac,amp,rat,r,GetHElement
        integer , dimension(lenof_sign) :: temp_sign
        character(len=*) , parameter :: this_routine='InitFCIMC_CAS'
        
        if(lenof_sign.ne.1) call stop_all(this_routine,"StartCAS currently does not work with complex walkers")
        if(tReadPops) call stop_all(this_routine,"StartCAS cannot work with with ReadPops")
        if(tStartSinglePart) call stop_all(this_routine,"StartCAS cannot work with StartSinglePart")
        if(tRestartHighPop) call stop_all(this_routine,"StartCAS cannot with with dynamically restarting calculations")

        write(6,*) "Initialising walkers proportional to a CAS diagonalisation..."
        write(6,'(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (" &
            ,(OccCASOrbs+VirtCASOrbs)/2,",",OccCASOrbs,")"
        DO I=NEl-OccCASorbs+1,NEl
            WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
            WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
            WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
        ENDDO
        WRITE(6,'(A)') " -------------------------------------------------------------------------------------------------"
        DO I=NEl+1,NEl+VirtCASOrbs
            WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
            WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
            WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
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

        write(6,*) "CAS Det is: "
        call write_det_len(6,CASDet,OccCASOrbs,.true.)
        call GetSym(CASDet,OccCASOrbs,G1,nBasisMax,CASSym)
        write(6,*) "Spatial symmetry of CAS determinants: ",CASSym%Sym%S
        write(6,*) "Ms of CAS determinants: ",CASSym%Ms
        if(tFixLz) then
            write(6,*) "Ml of CAS determinants: ",CASSym%Ml
        endif

        if(CASSym%Ml.ne.LzTot) call stop_all(this_routine,"Ml of CAS ref det does not match Ml of full reference det")
        if(CASSym%Ms.ne.0) call stop_all(this_routine,"CAS diagonalisation can only work with closed shell CAS spaces initially")
        if(CASSym%Sym%S.ne.HFSym%Sym%S) then
            call stop_all(this_routine,"Sym of CAS ref det does not match Sym of fulll reference det")
        endif

        !First, we need to generate all the excitations.
        call gndts(OccCASorbs,CASSpinBasisSize,CASBrr,nBasisMax,CASDetList,.true.,G1,tSpn,LMS,.true.,CASSym,nCASDet,CASDet)

        if(nCASDet.eq.0) call stop_all(this_routine,"No CAS determinants found.")
        write(6,*) "Number of symmetry allowed CAS determinants found to be: ",nCASDet
        Allocate(CASDetList(OccCASorbs,nCASDet),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating CASDetList")
        CASDetList(:,:)=0

        !Now fill up CASDetList...
        call gndts(OccCASorbs,CASSpinBasisSize,CASBrr,nBasisMax,CASDetList,.false.,G1,tSpn,LMS,.true.,CASSym,nCASDet,CASDet)

        !We have a complication here. If we calculate the hamiltonian from these CAS determinants, then we are not
        !including the mean-field generated from the other occupied orbitals. We need to either 'freeze' the occupied
        !orbitals and modify the 1 & two electron integrals, or add the other electrons back into the list. We do the latter.
        allocate(CASFullDets(NEl,nCASDet),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating CASFullDets")
        CASFullDets(:,:)=0

        do i=1,nCASDet
            do j=1,NEl-OccCASorbs
                CASFullDets(j,i)=ProjEDet(j)
            enddo
            do j=NEl-OccCASorbs+1,NEl
                CASFullDets(j,i)=CASDetList(j-(NEl-OccCASorbs),i)
            enddo
        enddo
        deallocate(CASDetList)

        write(6,*) "First CAS determinant in list is: "
        call write_det(6,CASFullDets(:,1),.true.)

        nEval=4
        write(6,"(A,I4,A)") "Calculating lowest ",nEval," eigenstates of CAS Hamiltonian..."
        Allocate(CkN(nCASDet,nEval), stat=ierr)
        CkN=0.D0
        Allocate(Ck(nCASDet,nEval),stat=ierr)
        Ck=0.D0
        Allocate(W(nEval),stat=ierr)    !Eigenvalues
        W=0.D0
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating")
        
        write(6,*) "Calculating hamiltonian..."
        allocate(nRow(nCASDet),stat=ierr)
        nRow=0
        ICMax=1
        tMC=.false.

        !HACK ALERT!! Need to fill up array in space of determinants, not HPHF functions.
        !Turn of tHPHFInts and turn back on when hamiltonian constructed.
        tHPHFInts=.false.

        CALL Detham(nCASDet,NEl,CASFullDets,Hamil,Lab,nRow,.true.,ICMax,GC,tMC)
        LenHamil=GC
        write(6,*) "Allocating memory for hamiltonian: ",LenHamil*2
        Allocate(Hamil(LenHamil),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating Hamil")
        Hamil=0.D0
        Allocate(Lab(LenHamil),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating Lab")
        Lab=0
        call Detham(nCASDet,NEl,CASFullDets,Hamil,Lab,nRow,.false.,ICMax,GC,tMC)

        CASRefEnergy=GETHELEMENT(1,1,HAMIL,LAB,NROW,NCASDET)
        write(6,*) "Energy of first CAS det is: ",CASRefEnergy

        !Turn back on HPHF integrals if needed.
        if(tHPHF) tHPHFInts=.true.

        if(abs(CASRefEnergy-Hii).gt.1.D-7) then
            call stop_all(this_routine,"CAS reference energy does not match reference energy of full space")
        endif

        !Lanczos
        NKRY1=NKRY+1
        NBLOCK=MIN(NEVAL,NBLK)
        LSCR=MAX(nCASDet*NEVAL,8*NBLOCK*NKRY)
        LISCR=6*NBLOCK*NKRY
        ALLOCATE(A(NEVAL,NEVAL),stat=ierr)
        CALL LogMemAlloc('A',NEVAL**2,8,this_routine,ATag,ierr)
        A=0.d0
        ALLOCATE(V(nCASDet*NBLOCK*NKRY1),stat=ierr)
        CALL LogMemAlloc('V',nCASDet*NBLOCK*NKRY1,8,this_routine,VTag,ierr)
        V=0.d0
        ALLOCATE(AM(NBLOCK*NBLOCK*NKRY1),stat=ierr)
        CALL LogMemAlloc('AM',NBLOCK*NBLOCK*NKRY1,8,this_routine,AMTag,ierr)
        AM=0.d0
        ALLOCATE(BM(NBLOCK*NBLOCK*NKRY),stat=ierr)
        CALL LogMemAlloc('BM',NBLOCK*NBLOCK*NKRY,8,this_routine,BMTag,ierr)
        BM=0.d0
        ALLOCATE(T(3*NBLOCK*NKRY*NBLOCK*NKRY),stat=ierr)
        CALL LogMemAlloc('T',3*NBLOCK*NKRY*NBLOCK*NKRY,8,this_routine,TTag,ierr)
        T=0.d0
        ALLOCATE(WT(NBLOCK*NKRY),stat=ierr)
        CALL LogMemAlloc('WT',NBLOCK*NKRY,8,this_routine,WTTag,ierr)
        WT=0.d0
        ALLOCATE(SCR(LScr),stat=ierr)
        CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag,ierr)
        SCR=0.d0
        ALLOCATE(ISCR(LIScr),stat=ierr)
        CALL LogMemAlloc('IScr',LIScr,4,this_routine,IScrTag,ierr)
        ISCR(1:LISCR)=0
        ALLOCATE(INDEX(NEVAL),stat=ierr)
        CALL LogMemAlloc('INDEX',NEVAL,4,this_routine,INDEXTag,ierr)
        INDEX(1:NEVAL)=0
        ALLOCATE(WH(nCASDet),stat=ierr)
        CALL LogMemAlloc('WH',nCASDet,8,this_routine,WHTag,ierr)
        WH=0.d0
        ALLOCATE(WORK2(3*nCASDet),stat=ierr)
        CALL LogMemAlloc('WORK2',3*nCASDet,8,this_routine,WORK2Tag,ierr)
        WORK2=0.d0
        ALLOCATE(V2(nCASDet,NEVAL),stat=ierr)
        CALL LogMemAlloc('V2',nCASDet*NEVAL,8,this_routine,V2Tag,ierr)
        V2=0.d0
!C..Lanczos iterative diagonalising routine
        CALL NECI_FRSBLKH(nCASDet,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT, &
     &  SCR,ISCR,INDEX,NCYCLE,B2L,.false.,.false.,.false.)
!Multiply all eigenvalues by -1.
        CALL DSCAL(NEVAL,-1.D0,W,1)
        if(CK(1,1).lt.0.D0) then
            do i=1,nCASDet
                CK(i,1)=-CK(i,1)
            enddo
        endif

        write(6,*) "Diagonalisation complete. Lowest energy CAS eigenvalues/corr E are: "
        do i=1,NEval
            write(6,*) i,W(i),W(i)-CASRefEnergy
        enddo

        TotWeight=0.D0
        nHPHFCAS=0
        do i=1,nCASDet
            if(tHPHF) then
                !Only allow valid HPHF functions
                call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                if(IsAllowedHPHF(iLutnJ)) then
                    nHPHFCAS=nHPHFCAS+1
                    if(.not.TestClosedShellDet(iLutnJ)) then
                        !Open shell. Weight is sqrt(2) of det weight.
                        TotWeight=TotWeight+(abs(CK(i,1))*sqrt(2.D0))
                        !Return this new weight to the CK array, so that we do not need to do this a second time.
                        CK(i,1)=CK(i,1)*sqrt(2.D0)
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
                        TotWeight=TotWeight+(abs(CK(i,1))*sqrt(2.D0))
                        !Return this new weight to the CK array, so that we do not need to do this a second time.
                        CK(i,1)=CK(i,1)*sqrt(2.D0)
                    else
                        !Closed Shell
                        TotWeight=TotWeight+abs(CK(i,1))
                    endif
                endif
            else
                TotWeight=TotWeight+abs(CK(i,1))
            endif
        enddo
        write(6,*) "Total weight of lowest eigenfunction: ",TotWeight
        if(tMomInv) write(6,*) "Converting into momentum-coupled space. Total MI functions: ",nHPHFCAS
        if(tHPHF) write(6,*) "Converting into HPHF space. Total HPHF CAS functions: ",nHPHFCAS

        if((InitialPart.eq.1).or.(InitialPart.ge.(InitWalkers*nNodes)-50)) then
            !Here, all the walkers will be assigned to the CAS wavefunction.
            !InitialPart = 1 by default
            write(6,"(A)") "All walkers specified in input will be distributed according to the CAS wavefunction."
            write(6,"(A)") "Shift will be allowed to vary from the beginning"
            write(6,"(A)") "Setting initial shift to equal CAS correlation energy",W(1)-CASRefEnergy
            DiagSft=W(1)-CASRefEnergy
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac=(real(InitWalkers,dp)* real(nNodes,dp))/TotWeight
        else
            !Here, not all walkers allowed will be initialised to the CAS wavefunction.
            write(6,"(A,I15,A)") "Initialising ",InitialPart, " walkers according to the CAS distribution."
            write(6,"(A,I15)") "Shift will remain fixed until the walker population reaches ",InitWalkers*nNodes
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac=real(InitialPart,dp)/TotWeight
            tSinglePartPhase=.true.
        endif

        !Now generate all excitations again, creating the required number of walkers on each one.
        DetIndex=1
        NoatHF(1)=0
        TotParts=0
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
                NoWalkers=int(amp)
                rat=amp-real(NoWalkers,dp)
                r=genrand_real2_dSFMT()
                if(abs(rat).gt.r) then
                    if(amp.lt.0.D0) then
                        NoWalkers=NoWalkers-1
                    else
                        NoWalkers=NoWalkers+1
                    endif
                endif

                if(NoWalkers.ne.0) then
                    call EncodeBitDet(CASFullDets(:,i),iLutnJ)
                    if(DetBitEQ(iLutnJ,iLutRef,NIfDBO)) then
                        !Check if this determinant is reference determinant, so we can count number on hf.
                        NoatHF(1) = NoWalkers
                    endif
                    call encode_det(CurrentDets(:,DetIndex),iLutnJ)
                    call clear_all_flags(CurrentDets(:,DetIndex))
                    temp_sign(1)=NoWalkers
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
                        CurrentH(DetIndex)=real(HDiagTemp,dp)-Hii
                    endif
                    DetIndex=DetIndex+1
                    TotParts(1)=TotParts(1)+abs(NoWalkers)
                endif
            endif   !End if desired node
        enddo

        TotWalkers=DetIndex-1   !This is the number of occupied determinants on each node
        TotWalkersOld=TotWalkers
        call sort(CurrentDets(:,1:TotWalkers),CurrentH(1:TotWalkers))

        !Set local&global variables
        TotPartsOld=TotParts
        call mpisumall(TotParts,AllTotParts)
        call mpisumall(NoatHF,AllNoatHF)
        call mpisumall(TotWalkers,AllTotWalkers)
        OldAllNoatHF=AllNoatHF
        AllTotWalkersOld=AllTotWalkers
        AllTotPartsOld=AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld=0.D0

        !Deallocate all the lanczos arrays now.
        deallocate(CK,W,CKN,Hamil,A,V,BM,T,WT,SCR,WH,Work2,V2,CASBrr,CASDet,CASFullDets,nRow,Lab,iscr,index,AM)
        call logmemdealloc(this_routine,ATag)
        call logmemdealloc(this_routine,VTag)
        call logmemdealloc(this_routine,BMTag)
        call logmemdealloc(this_routine,TTag)
        call logmemdealloc(this_routine,WTTag)
        call logmemdealloc(this_routine,SCRTag)
        call logmemdealloc(this_routine,WHTag)
        call logmemdealloc(this_routine,Work2Tag)
        call logmemdealloc(this_routine,V2Tag)
        call logmemdealloc(this_routine,iscrTag)
        call logmemdealloc(this_routine,indexTag)
        call logmemdealloc(this_routine,AMTag)

    end subroutine InitFCIMC_CAS 

!Routine to initialise the particle distribution according to the MP1 wavefunction.
!This hopefully will help with close-lying excited states of the same sym.
    subroutine InitFCIMC_MP1()
        use MomInv, only: IsAllowedMI
        use Determinants, only: GetH0Element3,GetH0Element4
        use SymExcit3 , only : GenExcitations3
        use CalcData , only : InitialPart
        real(dp) :: TotMP1Weight,amp,MP2Energy,PartFac,H0tmp,rat,r
        HElement_t :: hel,HDiagtemp
        integer :: iExcits,exflag,Ex(2,2),nJ(NEl),ic,DetIndex,iNode,NoWalkers,iInit
        integer(n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: temp_sign
        logical :: tAllExcitsFound,tParity
        character(len=*), parameter :: this_routine="InitFCIMC_MP1"

        if(lenof_sign.ne.1) call stop_all(this_routine,"StartMP1 currently does not work with complex walkers")
        if(tReadPops) call stop_all(this_routine,"StartMP1 cannot work with with ReadPops")
        if(tStartSinglePart) call stop_all(this_routine,"StartMP1 cannot work with StartSinglePart")
        if(tRestartHighPop) call stop_all(this_routine,"StartMP1 cannot with with dynamically restarting calculations")

        write(6,*) "Initialising walkers proportional to the MP1 amplitudes..."

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
        TotMP1Weight=1.D0
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
            if(tHPHF) then
                !Working in HPHF Space. Check whether determinant generated is an 'HPHF'
                call EncodeBitDet(nJ,iLutnJ)
                if(.not.IsAllowedHPHF(iLutnJ)) cycle
            elseif(tMomInv) then
                !Working in MI space, Check whether determinant generated is allowed
                call EncodeBitDet(nJ,iLutnJ)
                if(.not.IsAllowedMI(nJ,iLutnJ)) cycle
            endif
            iExcits=iExcits+1
            if(Ex(1,2).eq.0) then
                ic=1
            else
                ic=2
            endif
            if(tHPHF) then
                !Assume since we are using HPHF that the alpha and
                !beta orbitals of the same spatial orbital have the same
                !fock energies, so can consider either.
                hel=hphf_off_diag_helement(HFDet,nJ,iLutHF,iLutnJ)
            elseif(tMomInv) then
                hel=MI_off_diag_helement(HFDet,nJ,iLutHF,iLutnJ)
            else
                hel=get_helement(HFDet,nJ,ic,Ex,tParity)
            endif
            if(tUEG) then
                !This will calculate the MP2 energies without having to use the fock eigenvalues.
                !This is done via the diagonal determinant hamiltonian energies.
                H0tmp=getH0Element4(nJ,HFDet)
            else
                H0tmp=getH0Element3(nJ)
            endif
            H0tmp=Fii-H0tmp
            amp=hel/H0tmp
            TotMP1Weight=TotMP1Weight+abs(amp)
            MP2Energy=MP2Energy+(hel**2)/H0tmp
        enddo

        if((.not.tHPHF).and.(.not.tMomInv).and.(iExcits.ne.(nDoubles+nSingles))) then
            write(6,*) nDoubles,nSingles,iExcits
            call stop_all(this_routine,"Not all excitations accounted for in StartMP1")
        endif

        write(6,"(A,2G25.15)") "MP2 energy calculated: ",MP2Energy,MP2Energy+Hii

        if((InitialPart.eq.1).or.(InitialPart.ge.(InitWalkers*nNodes)-50)) then
            !Here, all the walkers will be assigned to the MP1 wavefunction.
            !InitialPart = 1 by default
            write(6,"(A)") "All walkers specified in input will be distributed according to the MP1 wavefunction."
            write(6,"(A)") "Shift will be allowed to vary from the beginning"
            write(6,"(A)") "Setting initial shift to equal MP2 correlation energy"
            DiagSft=MP2Energy
            !PartFac is the number of walkers that should reside on the HF determinant
            !in an intermediate normalised MP1 wavefunction. 
            PartFac=(real(InitWalkers,dp)* real(nNodes,dp))/TotMP1Weight
        else
            !Here, not all walkers allowed will be initialised to the MP1 wavefunction.
            write(6,"(A,I15,A)") "Initialising ",InitialPart, " walkers according to the MP1 distribution."
            write(6,"(A,I15)") "Shift will remain fixed until the walker population reaches ",InitWalkers*nNodes
            !PartFac is the number of walkers that should reside on the HF determinant
            !in an intermediate normalised MP1 wavefunction. 
            PartFac=real(InitialPart,dp)/TotMP1Weight
            tSinglePartPhase=.true.
        endif


        !Now generate all excitations again, creating the required number of walkers on each one.
        DetIndex=1
        TotParts=0
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
            if(tHPHF) then
                call EncodeBitDet(nJ,iLutnJ)
                if(.not.IsAllowedHPHF(iLutnJ)) cycle
            elseif(tMomInv) then
                call EncodeBitDet(nJ,iLutnJ)
                if(.not.IsAllowedMI(nJ,iLutnJ)) cycle
            endif

            iNode=DetermineDetNode(nJ,0)
            if(iProcIndex.eq.iNode) then
                if(Ex(1,2).eq.0) then
                    ic=1
                else
                    ic=2
                endif
                if(tHPHF) then
                    hel=hphf_off_diag_helement(HFDet,nJ,iLutHF,iLutnJ)
                elseif(tMomInv) then
                    hel=MI_off_diag_helement(HFDet,nJ,iLutHF,iLutnJ)
                else
                    hel=get_helement(HFDet,nJ,ic,Ex,tParity)
                endif
                if(tUEG) then
                    !This will calculate the MP2 energies without having to use the fock eigenvalues.
                    !This is done via the diagonal determinant hamiltonian energies.
                    H0tmp=getH0Element4(nJ,HFDet)
                else
                    H0tmp=getH0Element3(nJ)
                endif
                H0tmp=Fii-H0tmp
                !No parts on this det = PartFac*Amplitude
                amp=(hel/H0tmp)*PartFac
                NoWalkers=int(amp)
                rat=amp-real(NoWalkers,dp)

                r=genrand_real2_dSFMT()
                if(abs(rat).gt.r) then
                    if(amp.lt.0.D0) then
                        NoWalkers=NoWalkers-1
                    else
                        NoWalkers=NoWalkers+1
                    endif
                endif

                if(NoWalkers.ne.0) then
                    call EncodeBitDet(nJ,iLutnJ)
                    call encode_det(CurrentDets(:,DetIndex),iLutnJ)
                    call clear_all_flags(CurrentDets(:,DetIndex))
                    temp_sign(1)=NoWalkers
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
                        CurrentH(DetIndex)=real(HDiagTemp,dp)-Hii
                    endif
                    DetIndex=DetIndex+1
                    TotParts(1)=TotParts(1)+abs(NoWalkers)
                endif
            endif   !End if desired node

            
        enddo

        !Now for the walkers on the HF det
        if(iHFProc.eq.iProcIndex) then
            NoWalkers=int(PartFac)  !This will always be positive
            rat=PartFac-real(NoWalkers,dp)
            if(rat.lt.0.D0) call stop_all(this_routine,"Should not have negative weight on HF")
            r=genrand_real2_dSFMT()
            if(abs(rat).gt.r) NoWalkers=NoWalkers+1
            if(NoWalkers.ne.0) then
                call encode_det(CurrentDets(:,DetIndex),iLutHF)
                call clear_all_flags(CurrentDets(:,DetIndex))
                temp_sign(1)=NoWalkers
                call encode_sign(CurrentDets(:,DetIndex),temp_sign)
                if(tTruncInitiator) then
                    !Set initiator flag (always for HF)
                    call set_flag(CurrentDets(:,DetIndex),flag_is_initiator(1))
                    call set_flag(CurrentDets(:,DetIndex),flag_is_initiator(2))
                endif
                if(.not.tRegenDiagHEls) CurrentH(DetIndex)=0.D0
                DetIndex=DetIndex+1
                TotParts(1)=TotParts(1)+abs(NoWalkers)
                NoatHF(1) = NoWalkers
            else
                call stop_all(this_routine,"No walkers initialised on the HF det with StartMP1")
            endif
        else
            NoatHF(1)=0
        endif
            
        TotWalkers=DetIndex-1   !This is the number of occupied determinants on each node
        TotWalkersOld=TotWalkers
        call sort(CurrentDets(:,1:TotWalkers),CurrentH(1:TotWalkers))

        !Set local&global variables
        TotPartsOld=TotParts
        call mpisumall(TotParts,AllTotParts)
        call mpisumall(NoatHF,AllNoatHF)
        call mpisumall(TotWalkers,AllTotWalkers)
        OldAllNoatHF=AllNoatHF
        AllTotWalkersOld=AllTotWalkers
        AllTotPartsOld=AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld=0.D0

    end subroutine InitFCIMC_MP1

    subroutine CalcUEGMP2()
        use SymExcitDataMod, only: kPointToBasisFn
        use SystemData, only: ElecPairs,NMAXX,NMAXY,NMAXZ,OrbECutOff,tGCutoff,GCutoff, &
                                tMP2UEGRestrict,kiRestrict,kiMsRestrict,kjRestrict,kjMsRestrict, &
                                Madelung,tMadelung,tUEGFreeze,FreezeCutoff
        use GenRandSymExcitNUMod, only: FindNewDet
        use Determinants, only: GetH0Element4, get_helement_excit
        integer :: Ki(3),Kj(3),Ka(3),LowLoop,HighLoop,X,i,Elec1Ind,Elec2Ind,K,Orbi,Orbj
        integer :: iSpn,FirstA,nJ(NEl),a,Ex(2,2),kx,ky,kz,OrbB,FirstB
        integer :: ki2,kj2
        logical :: tParity,tMom
        real(dp) :: Ranger,mp2,mp2all,length,length_g,length_g_2
        HElement_t :: hel,H0tmp

        !Divvy up the ij pairs
        Ranger=real(ElecPairs)/real(nProcessors)
        LowLoop=int(iProcIndex*Ranger)+1
        Highloop=int((iProcIndex+1)*Ranger)

        if((iProcIndex+1).eq.nProcessors) Highloop=ElecPairs
        if(iProcIndex.eq.0) then
            if(lowLoop.ne.1) write(6,*) "Error here!"
        endif
        write(6,*) "Total ij pairs: ",ElecPairs
        write(6,*) "Considering ij pairs from: ",LowLoop," to ",HighLoop  
!        write(6,*) "HFDet: ",HFDet(:)

        do i=LowLoop,HighLoop   !Looping over electron pairs on this processor

            X=ElecPairs-i
            K=INT((SQRT(8.D0*REAL(X,dp)+1.D0)-1.D0)/2.D0)
            Elec1Ind=NEl-1-K
            Elec2Ind=NEl-X+((K*(K+1))/2)
            Orbi=HFDet(Elec1Ind)
            Orbj=HFDet(Elec2Ind)
            Ki=G1(Orbi)%k
            Kj=G1(Orbj)%k
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
                write(6,*) "Restricting calculation to i,j pair: ",Orbi,Orbj
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

!            write(6,*) "ijpair: ",Orbi,Orbj

            if((iSpn.eq.3).or.(iSpn.eq.1)) then
                if(iSpn.eq.3) then
                    FirstA=2    !Loop over alpha
                else
                    FirstA=1    !Loop over beta
                endif

                do a=FirstA,nBasis,2
                    !Loop over all a

                    !Reject if a is occupied
                    if(IsOcc(iLutHF,a)) cycle

                    Ka=G1(a)%k

                    !Find k labels of b
                    kx=Ki(1)+Kj(1)-Ka(1)
                    if(abs(kx).gt.NMAXX) cycle
                    ky=Ki(2)+Kj(2)-Ka(2)
                    if(abs(ky).gt.NMAXY) cycle
                    kz=Ki(3)+Kj(3)-Ka(3)
                    if(abs(kz).gt.NMAXZ) cycle
                    if(tGCutoff) then
                        length_g=real((kx-kj(1))**2+(ky-kj(2))**2+(kz-kj(3))**2)
                        length_g_2=real((kx-ki(1))**2+(ky-ki(2))**2+(kz-ki(3))**2)
                        if(length_g.gt.gCutoff.and.length_g_2.gt.gCutoff) cycle
                    endif
                    length=real((kx**2)+(ky**2)+(kz**2))
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
                    if(OrbB.ge.a) cycle

                    !Find det
!                    write(6,*) "OrbB: ",OrbB
                    call FindNewDet(HFDet,nJ,Elec1Ind,Elec2Ind,a,OrbB,Ex,tParity)
                    !Sum in mp2 contrib
                    hel=get_helement_excit(HFDet,nJ,2,Ex,tParity)

                    H0tmp=getH0Element4(nJ,HFDet)
                    H0tmp=Fii-H0tmp
                    if(tMadelung) then
                        H0tmp=H0tmp+Madelung
                    endif
                    mp2=mp2+(hel**2)/H0tmp
!                    write(6,*) (hel**2),H0tmp
                enddo

            elseif(iSpn.eq.2) then
                do a=1,nBasis
                    !Loop over all a
!                    write(6,*) "a: ",a

                    !Reject if a is occupied
                    if(IsOcc(iLutHF,a)) cycle

                    Ka=G1(a)%k

                    !Find k labels of b
                    kx=Ki(1)+Kj(1)-Ka(1)
                    if(abs(kx).gt.NMAXX) cycle
                    ky=Ki(2)+Kj(2)-Ka(2)
                    if(abs(ky).gt.NMAXY) cycle
                    kz=Ki(3)+Kj(3)-Ka(3)
                    if(abs(kz).gt.NMAXZ) cycle
                    if(tGCutoff) then
                        length_g=real((kx-kj(1))**2+(ky-kj(2))**2+(kz-kj(3))**2)
                        length_g_2=real((kx-ki(1))**2+(ky-ki(2))**2+(kz-ki(3))**2)
                        if(length_g.gt.gCutoff.and.length_g_2.gt.gCutoff) cycle
                    endif
                    length=real((kx**2)+(ky**2)+(kz**2))
                    if(length.gt.OrbECutoff) cycle

                    !Find the actual k orbital
                    if(is_beta(a)) then
                        !want alpha b orbital
                        OrbB=kPointToBasisFn(kx,ky,kz,2)
                    else
                        !want beta
                        OrbB=kPointToBasisFn(kx,ky,kz,1)
                    endif

                    !Reject k orbital if it is occupied or gt a
                    if(IsOcc(iLutHF,OrbB)) cycle
                    if(OrbB.ge.a) cycle

!                    write(6,*) "OrbB: ",OrbB
                    !Find det
                    call FindNewDet(HFDet,nJ,Elec1Ind,Elec2Ind,a,OrbB,Ex,tParity)
                    !Sum in mp2 contrib
                    hel=get_helement_excit(HFDet,nJ,2,Ex,tParity)
                    H0tmp=getH0Element4(nJ,HFDet)
                    H0tmp=Fii-H0tmp
                    if(tMadelung) then
                        H0tmp=H0tmp+Madelung
                    endif
                    mp2=mp2+(hel**2)/H0tmp
!                    write(6,*) (hel**2),H0tmp
                enddo
            endif

        enddo

!        write(6,*) "mp2: ",mp2
        mp2all=0.D0
        
        !Sum contributions across nodes.
        call MPISumAll(mp2,mp2all)
        write(6,"(A,2G25.15)") "MP2 energy calculated: ",MP2All,MP2All+Hii
        call flush(6)
        call stop_all("CalcUEGMP2","Dying after calculation of MP2 energy...")

    end subroutine CalcUEGMP2
            
    !Count excitations using ajwt3's old excitation generators which can handle non-abelian and
    !k-point symmetry
    SUBROUTINE CountExcitsOld(HFDet,exflag,nSing,nDoub)
        use SystemData , only : tUseBrillouin
        INTEGER , INTENT(IN) :: HFDet(NEl),exflag
        INTEGER , INTENT(OUT) :: nSing,nDoub
        LOGICAL :: TempUseBrill
        INTEGER :: nExcitMemLen,nJ(NEl),iMaxExcit,nStore(6),ierr,excitcount,iExcit
        INTEGER , ALLOCATABLE :: EXCITGEN(:)
        character(len=*) , parameter :: this_routine='CountExcitsOld'

        nSing=0
        nDoub=0

        !However, we have to ensure that brillouins theorem isn't on!
        IF(tUseBrillouin) THEN
            TempUseBrill=.true.
            tUseBrillouin=.false.
        ELSE
            TempUseBrill=.false.
        ENDIF

        iMaxExcit=0
        nStore(1:6)=0
        nJ(:)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
        ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating excitation generator")
        EXCITGEN(:)=0
        nJ(:)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
    !    CALL GetSymExcitCount(EXCITGEN,DetConn)
        excitcount=0

    lp2: do while(.true.)
            CALL GenSymExcitIt2(HFDet,nEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
            IF(nJ(1).eq.0) exit lp2
            IF(iExcit.eq.1) THEN
                nSing=nSing+1
            ELSEIF(iExcit.eq.2) THEN
                nDoub=nDoub+1
            ELSE
                CALL Stop_All(this_routine,"Trying to generate more than doubles!")
            ENDIF
        enddo lp2
        tUseBrillouin=TempUseBrill

        if(iMaxExcit.ne.(nSing+nDoub)) then
            call stop_all("CountExcitsOld","Error in counting old excitations")
        endif

    END SUBROUTINE CountExcitsOld

    SUBROUTINE DeallocFCIMCMemPar()
        use nElRDMMod, only: DeallocateRDM
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocFciMCMemPar'


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
        
        DEALLOCATE(HFDet)
        CALL LogMemDealloc(this_routine,HFDetTag)
        DEALLOCATE(iLutHF)
        DEALLOCATE(iLutRef)
        DEALLOCATE(ProjEDet)
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

        ! Cleanup excitation generation storage
        call clean_excit_gen_store (fcimc_excit_gen_store)

        ! Cleanup linear combination projected energy
        call clean_linear_comb ()

        ! Cleanup storage for spin projection
        call clean_yama_store ()


!There seems to be some problems freeing the derived mpi type.
!        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!Free the mpi derived type that we have created for the hashes.
!            CALL MPI_Type_free(mpilongintegertype,error)
!            IF(error.ne.MPI_SUCCESS) THEN
!                CALL MPI_Error_string(error,message,length,temp)
!                IF(temp.ne.MPI_SUCCESS) THEN
!                    WRITE(6,*) "REALLY SERIOUS PROBLEMS HERE!",temp
!                    CALL FLUSH(6)
!                ENDIF
!                WRITE(6,*) message(1:length)
!                WRITE(6,*) "ERROR FOUND"
!                CALL FLUSH(6)
!            ENDIF
!        ENDIF

    END SUBROUTINE DeallocFCIMCMemPar

   subroutine SetupValidSpawned(WalkerListSize)
      use CalcData, only: MemoryFacSpawn
      implicit none
      integer, intent(in) :: WalkerListSize
      integer ierr,i,j
      real(dp) Gap
      MaxSpawned=NINT(MemoryFacSpawn*WalkerListSize)
!            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned
            
      WRITE(6,*) "*Direct Annihilation* in use...Explicit load-balancing disabled."
      ALLOCATE(ValidSpawnedList(0:nNodes-1),stat=ierr)
      !InitialSpawnedSlots is now filled later, once the number of particles wanted is known
      !(it can change according to the POPSFILE).
      ALLOCATE(InitialSpawnedSlots(0:nNodes-1),stat=ierr)
!InitialSpawnedSlots now holds the first free position in the newly-spawned list for each processor, so it does not need to be reevaluated each iteration.
      MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
      Gap=REAL(MaxSpawned)/REAL(nNodes)
      do j=0,nNodes-1
          InitialSpawnedSlots(j)=NINT(Gap*j)+1
      enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
      ValidSpawnedList(:)=InitialSpawnedSlots(:)

   end subroutine

END MODULE FciMCParMod

!This routine will change the reference determinant to DetCurr. It will also re-zero all the energy estimators, since they now correspond to
!projection onto a different determinant.
    SUBROUTINE ChangeRefDet(DetCurr)
        use FciMCParMod
        use SystemData , only : NEl
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),i

        do i=1,NEl
            FDet(i)=DetCurr(i)
        enddo

        WRITE(6,"(A)") "*** Changing the reference determinant ***"
        WRITE(6,"(A)") "Switching reference and zeroing energy counters - restarting simulation"
!        
!Initialise variables for calculation on each node
        Iter=1
        CALL DeallocFCIMCMemPar()
        IF(iProcIndex.eq.Root) THEN
            CLOSE(fcimcstats_unit)
            IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(initiatorstats_unit)
            IF(tLogComplexPops) CLOSE(complexstats_unit)
        ENDIF
        IF(TDebug) CLOSE(11)
        CALL SetupParameters()
        CALL InitFCIMCCalcPar()

    END SUBROUTINE ChangeRefDet
            

!This is the same as BinSearchParts1, but this time, it searches though the full list of determinants created by the full diagonalizer when the histogramming option is on.
!This is outside the module so it is accessible to AnnihilateMod
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

!    WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!    CALL FLUSH(6)
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
!        WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
        Comp=DetBitLT(FCIDets(:,N),iLut(:),NIfDBO)

        IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
            tSuccess=.true.
            PartInd=N
            RETURN
        ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
            i=N
        ELSEIF(i.eq.N) THEN


            IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
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
    
