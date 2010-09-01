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
                          tNoBrillouin,tKPntSym
    use bit_reps, only: NIfD, NIfTot, NIfDBO, NIfY, decode_bit_det, &
                        encode_bit_rep, encode_det, extract_bit_rep, &
                        test_flag, set_flag, extract_flags, &
                        flag_is_initiator, clear_all_flags, set_flag_general
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
                        tChangeProjEDet, tCheckHighestPop
    use HPHFRandExcitMod, only: FindExcitBitDetSym, GenRandHPHFExcit, &
                                gen_hphf_excit
    use Determinants, only: FDet, get_helement, write_det, &
                            get_helement_det_only
    USE DetCalcData , only : ICILevel,nDet,Det,FCIDetIndex
    use GenRandSymExcitNUMod, only: gen_rand_excit, GenRandSymExcitNU, &
                                    ScratchSize
    use GenRandSymExcitCSF, only: gen_csf_excit
    use IntegralsData , only : fck,NMax,UMat,tPartFreezeCore,NPartFrozen,NHolesFrozen,tPartFreezeVirt,NVirtPartFrozen,NElVirtFrozen
    USE Logging , only : iWritePopsEvery,TPopsFile,iPopsPartEvery,tBinPops,tHistSpawn,iWriteHistEvery,tHistEnergies,IterShiftBlock,AllHistInitPops
    USE Logging , only : BinRange,iNoBins,OffDiagBinRange,OffDiagMax,AllHistInitPopsTag,tLogComplexPops
    USE Logging , only : tPrintFCIMCPsi,tCalcFCIMCPsi,NHistEquilSteps,tPrintOrbOcc,StartPrintOrbOcc,tPrintOrbOccInit
    USE Logging , only : tHFPopStartBlock,tIterStartBlock,IterStartBlocking,HFPopStartBlocking,tInitShiftBlocking,tHistHamil,iWriteHamilEvery,HistInitPopsTag
    USE Logging , only : OrbOccs,OrbOccsTag,tPrintPopsDefault,iWriteBlockingEvery,tBlockEveryIteration,tHistInitPops,HistInitPopsIter,HistInitPops
    USE Logging , only : FCIMCDebug
    USE SymData , only : nSymLabels
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE Parallel
    USE FciMCData
    USE AnnihilationMod
    use DetBitops, only: EncodeBitDet, DetBitEQ, DetBitLT, FindExcitBitDet, &
                         FindBitExcitLevel, countbits
    use csf, only: get_csf_bit_yama, iscsf, csf_orbital_mask, get_csf_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_spawn_sign, hphf_off_diag_helement_spawn
    use util_mod, only: choose,abs_int_sign,abs_int8_sign
    use constants, only: dp, int64, n_int, lenof_sign
    use soft_exit, only: ChangeVars 
    use FciMCLoggingMod, only: FinaliseBlocking, FinaliseShiftBlocking, &
                               PrintShiftBlocking, PrintBlocking, &
                               SumInErrorContrib, WriteInitPops, &
                               InitErrorBlocking, InitShiftErrorBlocking, &
                               SumInShiftErrorContrib
    use RotateOrbsMod, only: RotateOrbs
    use NatOrbsMod, only: PrintOrbOccs
    use spin_project, only: tSpinProject, spin_proj_interval, &
                            spin_proj_gamma, get_spawn_helement_spin_proj, &
                            generate_excit_spin_proj, attempt_die_spin_proj, &
                            iter_data_spin_proj, test_spin_proj, &
                            spin_proj_shift
#ifdef __DEBUG                            
    use DeterminantData, only: write_det
#endif

    implicit none

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)

        real(dp) :: Weight, Energyxw
        INTEGER :: error
        LOGICAL :: TIncrement,tWritePopsFound,tSoftExitFound,tSingBiasChange
        REAL(4) :: s,etime,tstart(2),tend(2)
        INTEGER(int64) :: MaxWalkers,MinWalkers
        real*8 :: AllTotWalkers,MeanWalkers,Inpair(2),Outpair(2)
        integer, dimension(lenof_sign) :: tmp_sgn
        integer :: tmp_int(lenof_sign)
        real(dp) :: grow_rate
        logical :: truncinit_tmp

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

        ! Initial output
        call WriteFciMCStatsHeader()
        ! Prepend a # to the initial status line so analysis doesn't pick up
        ! repetitions in the FCIMCSTATS file from restarts.
!        write (6,'("#")', advance='no')
        write (fcimcstats_unit,'("#")', advance='no')
        call WriteFCIMCStats()

        ! Put a barrier here so all processes synchronise before we begin.
        call MPIBarrier(error)


!Start MC simulation...
        TIncrement=.true.   !If TIncrement is true, it means that when it comes out of the loop, it wants to subtract 1 from the Iteration count to get the true number of iterations
        Iter=1
        do while(Iter.le.NMCyc)   !Iter=1,NMCyc
!Main iteration loop...
!            WRITE(6,*) 'Iter',Iter

            s=etime(tstart)
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
!                    print*, 'FCIQMC', iter, iter_data_fciqmc%nborn, &
!                                            iter_data_fciqmc%ndied
                endif
            endif

            ! Are we projecting the spin out between iterations?
            if (tSpinProject .and. (mod(Iter, spin_proj_interval) == 0 .or. &
                                    spin_proj_interval == -1)) then
                truncinit_tmp = tTruncInitiator
                tTruncInitiator = .false.
                call PerformFciMCycPar (generate_excit_spin_proj, &
                                       attempt_create_normal, &
                                       get_spawn_helement_spin_proj, &
                                       null_encode_child, &
                                       new_child_stats_normal, &
                                       attempt_die_spin_proj, &
                                       iter_data_spin_proj)
                tTruncInitiator = truncinit_tmp
!                print*, 'SPIN', iter, iter_data_spin_proj%nborn, &
!                                      iter_data_spin_proj%ndied
!                print*, 'TOTWALKRS', totwalkers
!                do i = 1, TotWalkers
!                    call extract_sign (CurrentDets(:,i), tmp_sgn)
!                    print*, 'det: ', i, tmp_sgn, &
!                            real(tmp_sgn(1)) / real(totparts(1))
!                enddo
            endif

            s=etime(tend)
            IterTime=IterTime+(tend(1)-tstart(1))

!            IF(tBlockEveryIteration) THEN
!                Inpair(1)=REAL(HFIter,dp)
!                Inpair(2)=ENumIter
!                CALL MPISumRoot(Inpair,2,Outpair,Root)
!                IterEnergy=Outpair(2)/Outpair(1)
!                IF(tErrorBlocking.and.(iProcIndex.eq.Root)) CALL SumInErrorContrib(Iter,Outpair(2),Outpair(1))
!                ENumIter=0.D0
!                HFIter=0
!            ENDIF

            if (mod(Iter, StepsSft) == 0) then

                ! Calculate the a new value for the shift (amongst other
                ! things). Generally, collate information from all processors,
                ! update statistics and output them to the user.
!>>>!                iter_data_fciqmc%update_growth = iter_data_fciqmc%update_growth + iter_data_spin_proj%update_growth
                if (tCCMC) then
                    call calculate_new_shift_wrapper (iter_data_ccmc, &
                                                      TotParts)
                else
                    call calculate_new_shift_wrapper (iter_data_fciqmc, &
                                                      TotParts)
                endif

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

                !    call MPIBCast (spin_proj_shift, root)

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

!This routine will check for a CHANGEVARS file and change the parameters of the calculation accordingly.
                CALL ChangeVars(tSingBiasChange,tSoftExitFound,tWritePopsFound)
                IF(tSoftExitFound) THEN
                    TIncrement=.false.
                    EXIT
                ENDIF
                IF(tWritePopsFound) THEN
!We have explicitly asked to write out the POPSFILE from the CHANGEVARS file.
                    CALL WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
                ENDIF
                IF(tSingBiasChange) THEN
                    CALL CalcApproxpDoubles()
                ENDIF
            
            ENDIF

!            IF(mod(Iter,iWriteBlockingEvery).eq.0) THEN
!                !Every 100 update cycles, write out a new blocking file.
!                IF(tErrorBlocking.and.(Iter.gt.IterStartBlocking)) CALL PrintBlocking(Iter) 
!                IF(tShiftBlocking.and.(Iter.gt.(VaryShiftIter+IterShiftBlock))) CALL PrintShiftBlocking(Iter)
!            ENDIF

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
        Energyxw=(ProjectionE)

        IF(tHistEnergies) CALL WriteHistogramEnergies()

        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
        ENDIF

! Print out some load balancing stats nicely to end.
        CALL MPIReduce(TotWalkers,MPI_MAX,MaxWalkers)
        CALL MPIReduce(TotWalkers,MPI_MIN,MinWalkers)
        CALL MPIAllReduce(Real(TotWalkers,dp),MPI_SUM,AllTotWalkers)
        if (iProcIndex.eq.Root) then
            MeanWalkers=AllTotWalkers/nProcessors
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
                            pGen, HEl, tFilled, scratch1, scratch2, scratch3)

                use SystemData, only: nel
                use bit_reps, only: niftot
                use GenRandSymExcitNUMod, only: scratchsize
                use constants, only: n_int,dp
                implicit none

                integer, intent(in) :: nI(nel) 
                integer(kind=n_int), intent(in) :: iLutI(0:niftot)
                integer, intent(in) :: exFlag
                integer, intent(inout) :: scratch1(scratchsize)
                integer, intent(inout) :: scratch2(scratchsize)
                integer, intent(inout) :: scratch3(scratchsize)
                integer, intent(out) :: nJ(nel) 
                integer(kind=n_int), intent(out) :: iLutJ(0:niftot)
                integer, intent(out) :: ic, ex(2,2)
                real*8, intent(out) :: pGen
                logical, intent(inout) :: tFilled
                logical, intent(out) :: tParity
                HElement_t, intent(out) :: HEl
            end subroutine
        end interface

        call assign_proc (ptr_excit_generator, gen)
    end subroutine

    subroutine set_attempt_create (attempt_create)
        use, intrinsic :: iso_c_binding
        implicit none
        interface
            function attempt_create (get_spawn_helement, nI, iLutI, wSign, &
                                     nJ, iLutJ, prob, HElGen, ic, ex, tPar, exLevel, part_type)&
                                     result(child)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int, dp, lenof_sign
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel), part_type 
                integer(kind=n_int), intent(in) :: iLutI(0:nIfTot)
                integer(kind=n_int), intent(inout) :: iLutJ(0:nIfTot)
                integer, intent(in) :: ic, ex(2,2), exLevel
                integer, dimension(lenof_sign), intent(in) :: wSign
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
            subroutine new_child_stats (iter_data, iLutI, iLutJ, ic, &
                                        walkExLevel, child)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int, lenof_sign
                use FciMCData, only: fcimc_iter_data
                implicit none
                integer(kind=n_int), intent(in) :: ilutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, walkExLevel
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
                             exFlag, IC, ex, tParity, pGen, HEl, tFilled, &
                             scratch1, scratch2, scratch3)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use GenRandSymExcitNUMod, only: scratchsize
                use constants, only: dp,n_int
                implicit none
                integer, intent(in) :: nI(nel)
                integer(kind=n_int), intent(in) :: iLutI(0:niftot)
                integer, intent(in) :: exFlag
                integer, intent(inout) :: scratch1(scratchsize)
                integer, intent(inout) :: scratch2(scratchsize)
                integer, intent(inout) :: scratch3(scratchsize)
                integer, intent(out) :: nJ(nel) 
                integer(kind=n_int), intent(out) :: iLutJ(0:niftot)
                integer, intent(out) :: ic, ex(2,2)
                real(dp), intent(out) :: pGen
                logical, intent(inout) :: tFilled
                logical, intent(out) :: tParity
                HElement_t, intent(out) :: HEl
            end subroutine
            function attempt_create (get_spawn_helement, nI, iLutI, wSign, &
                                     nJ, iLutJ, prob, HElGen, ic, ex, tPar, exLevel, part_type)&
                                     result(child)
                use systemdata, only: nel
                use bit_reps, only: niftot
                use constants, only: dp, n_int, lenof_sign
                implicit none
                integer, intent(in) :: nI(nel), nJ(nel), part_type
                integer(kind=n_int), intent(in) :: iLutI(0:nifTot)
                integer(kind=n_int), intent(inout) :: iLutJ(0:nIfTot)
                integer, intent(in) :: ic, ex(2,2), exLevel
                integer, dimension(lenof_sign), intent(in) :: wSign
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
            subroutine new_child_stats (iter_data, iLutI, iLutJ, ic, &
                                        walkExLevel, child)
                use SystemData, only: nel
                use bit_reps, only: niftot
                use constants, only: n_int, lenof_sign
                use FciMCData, only: fcimc_iter_data
                implicit none
                integer(kind=n_int), intent(in) :: ilutI(0:niftot), iLutJ(0:niftot)
                integer, intent(in) :: ic, walkExLevel
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
        integer :: VecSlot, j, p, error
        integer :: DetCurr(nel), nJ(nel), FlagsCurr, parent_flags
        integer, dimension(lenof_sign) :: SignCurr, child
        integer(kind=n_int) :: iLutnJ(0:niftot)
        integer :: IC, walkExcitLevel, ex(2,2), TotWalkersNew, part_type
        integer, dimension(ScratchSize) :: scratch1, scratch2, scratch3
        integer(int64) :: tot_parts_tmp(lenof_sign)
        logical :: tFilled, tParity
        real(dp) :: prob, HDiagCurr
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
                    WRITE(6,'(A,F10.5)') 'Beginning truncated initiator calculation and reducing tau by a factor of 10. New tau is : ',Tau
                ENDIF
                CALL MPIBCast(Tau,root)
            ENDIF
            tTruncInitiator=.true.
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
        iPartBloom = 0 ! Max number spawned from an excitation
        ! Next free position in newly spawned list.
        ValidSpawnedList = InitialSpawnedSlots
        
        ! Reset iteration specific counts
        ! n.b. this can be extracted to a separate function. This could even
        !      be function-pointerised!
        iter_data%nborn = 0
        iter_data%ndied = 0
        iter_data%nannihil = 0
        iter_data%naborted = 0
        tot_parts_tmp = TotParts

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

            !Debug output.
            IFDEBUG(FCIMCDebug,3) write(6,"(A,I10)",advance='no') "TW:", j
            IFDEBUG(FCIMCDebug,3) call WriteBitDet(6,CurrentDets(:,j),.true.)
            IFDEBUG(FCIMCDebug,3) call Flush(6) 

            ! Decode determinant from (stored) bit-representation.
            call extract_bit_rep (CurrentDets(:,j), DetCurr, SignCurr, &
                                  FlagsCurr)

            ! TODO: The next couple of bits could be done automatically

            ! We only need to find out if determinant is connected to the
            ! reference (so no ex. level above 2 required, except for HPHF
            ! or truncated etc.)
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

            ! TODO: Ensure that the HF determinant has its flags setup
            !       correctly at the start of a run.
            if (tTruncInitiator) call CalcParentFlag (j, VecSlot, Iter, &
                                                      parent_flags)

            ! Indicate that the scratch storage used for excitation generation
            ! from the same walker has not been filled (it is filled when we
            ! excite from the first particle on a determinant).
            tfilled = .false.

            !print*, '*******************'
            !call test_spin_proj (DetCurr, CurrentDets(:,j))

            ! Loop over the 'type' of particle. 
            ! lenof_sign == 1 --> Only real particles
            ! lenof_sign == 2 --> complex walkers
            !                 --> part_type == 1, 2; real and complex walkers
            do part_type=1,lenof_sign

                ! Loop over all the particles of a given type on the 
                ! determinant. CurrentSign gives number of walkers. Multiply 
                ! up by noMCExcits if attempting multiple excitations from 
                ! each walker (default 1)
                do p = 1, abs(SignCurr(part_type)) * noMCExcits
                    ! Generate a (random) excitation
                    call generate_excitation (DetCurr, CurrentDets(:,j), nJ, &
                                   ilutnJ, exFlag, IC, ex, tParity, prob, HElGen,&
                                   tFilled, scratch1, scratch2, scratch3)

                    ! If a valid excitation, see if we should spawn children.
                    if (.not. IsNullDet(nJ)) then
                        child = attempt_create (get_spawn_helement, DetCurr, &
                                            CurrentDets(:,j), SignCurr, &
                                            nJ,iLutnJ, Prob, HElGen, IC, ex, tParity, &
                                            walkExcitLevel,part_type)
                    else
                        child = 0
                    endif

                    !Children have been chosen to be spawned.
#ifdef __CMPLX
                    if ((child(1).ne.0).or.(child(2).ne.0)) then
                        IFDEBUG(FCIMCDebug,3) write(6,"(A,2I4,A)",advance='no') "Creating ",child(1),child(2)," particles: "
                        IFDEBUG(FCIMCDebug,3) write(6,"(A,2I4,A)",advance='no') "Parent flag: ",parent_flags,part_type
#else
                    if (child(1).ne.0) then
                        IFDEBUG(FCIMCDebug,3) write(6,"(A,I4,A)",advance='no') "Creating ",child(1)," particles: "
                        IFDEBUG(FCIMCDebug,3) write(6,"(A,2I4,A)",advance='no') "Parent flag: ",parent_flags,part_type
#endif                        
                        IFDEBUG(FCIMCDebug,3) CALL Write_Det(6,nJ,.true.)
                        IFDEBUG(FCIMCDebug,3) CALL FLUSH(6) 
                        ! We know we want to create a particle of this type.
                        ! Encode the bit representation if it isn't already.
                        call encode_child (CurrentDets(:,j), iLutnJ, ic, ex)

                        call new_child_stats (iter_data, CurrentDets(:,j), &
                                              iLutnJ, ic, walkExcitLevel, &
                                              child)

                        call create_particle (iLutnJ, child, parent_flags, part_type)
                    
                    endif ! (child /= 0). Child created

                enddo ! Cycling over mulitple particles on same determinant.

            enddo   ! Cycling over 'type' of particle on a given determinant.

            ! DEBUG
            ! if (VecSlot > j) call stop_all (this_routine, 'vecslot > j')
            call walker_death (attempt_die, iter_data, DetCurr, &
                               CurrentDets(:,j), HDiagCurr, SignCurr, VecSlot)

        enddo ! Loop over determinants.

        IFDEBUG(FCIMCDebug,2) WRITE(6,*) "Finished loop over determinants"

        ! SumWalkersCyc calculates the total number of walkers over an update
        ! cycle on each process.
        SumWalkersCyc = SumWalkersCyc + int(sum(TotParts), int64)

        ! Since VecSlot holds the next vacant slot in the array, TotWalkers
        ! should be one less than this.
        TotWalkersNew = VecSlot - 1

        IF((tHistInitPops.and.(MOD(Iter,HistInitPopsIter).eq.0))    &
            .or.tPrintHighPop) THEN
            CALL FindHighPopDet(TotWalkersNew)
            IF(iProcIndex.eq.0) WRITE(6,'(A)') 'Writing out the spread of the initiator determinant populations.'
            CALL WriteInitPops(Iter+PreviousCycles)
        ENDIF

        ! Print bloom/memory warnings
        call end_iteration_print_warn (totWalkersNew)
        call halt_timer (walker_time)
        
        ! For the direct annihilation algorithm. The newly spawned 
        ! walkers should be in a seperate array (SpawnedParts) and the other 
        ! list should be ordered.
        call set_timer (annihil_time, 30)
        call DirectAnnihilation (totWalkersNew, iter_data,.false.) !.false. for not single processor
        TotWalkers=TotWalkersNew
        CALL halt_timer(Annihil_Time)
        IFDEBUG(FCIMCDebug,2) WRITE(6,*) "Finished Annihilation step"
        
        ! Update iteration data
        iter_data%update_growth = iter_data%update_growth + iter_data%nborn &
                                - iter_data%ndied - iter_data%nannihil &
                                - iter_data%naborted
        iter_data%update_iters = iter_data%update_iters + 1

    end subroutine

    subroutine new_child_stats_normal (iter_data, iLutI, iLutJ, ic, &
                                       walkExLevel, child)

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel
        integer, dimension(lenof_sign), intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(n_int) :: iUnused
       
        ! Count the number of children born
        NoBorn = NoBorn + sum(abs(child))
        iter_data%nborn = iter_data%nborn + abs(child)

        if (ic == 1) SpawnFromSing = SpawnFromSing + sum(abs(child))

        if (sum(abs(child)) > abs(iPartBloom)) then
            iPartBloom = sign(sum(abs(child)), 2*ic - 3)
        endif

        ! Avoid compiler warnings
        iUnused = iLutI(0); iUnused = iLutJ(0); iUnused = walkExLevel

    end subroutine
                        
    subroutine create_particle (iLutJ, child, parent_flags, part_type)

        ! Create a child in the spawned particles arrays. We spawn particles
        ! into a separate array, but non-contiguously. The processor that the
        ! newly-spawned particle is going to be sent to has to be determined,
        ! and then it will be put into the appropriate element determined by
        ! ValidSpawnedList

        integer(kind=n_int), intent(in) :: iLutJ(0:niftot)
        integer, dimension(lenof_sign), intent(in) :: child
        integer, intent(in) :: parent_flags
        integer, intent(in) :: part_type        !This is the 'type' of the parent particle (i.e. real/imag)
        logical :: parent_init
        integer :: proc,j

        proc = DetermineDetProc(iLutJ)    ! 0 -> nProcessors-1

        ! This will also set the flag of the walker(parentInitiator) to be either 0 or 1
        ! according to if its parent is inside or outside the active space.
        call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), iLutJ, &
                            child, parent_flags)

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

        ! Has there been a particle bloom this iteration?
        if ( abs(iPartBloom) > InitiatorWalkNo) then
            write (6, bloom_warn_string, advance='no') iter, abs(iPartBloom)
            if (iPartBloom > 0) then
                write (6, '("double excit.")')
            else
                write (6, '("single excit.")')
            endif
        endif

        ! Too many particles?
        rat = real(TotWalkersNew,dp) / real(MaxWalkersPart,dp)
        if (rat > 0.95) then
            write (6, '(a)') '*WARNING* - Number of particles/determinants &
                             &has increased to over 95% of MaxWakersPart.'
            call flush(6)
        end if

        ! Are ony of the sublists near the end of their alloted space?
        if (nProcessors > 1) then
            do i = 0, nProcessors-1
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
    end subroutine


    subroutine CalcParentFlag(j, VecSlot, Iter, parent_flags)
!In the CurrentDets array, the flag at NIfTot refers to whether that determinant *itself* is an initiator or not.    
!We need to decide if this willchange due to the determinant acquiring a certain population, or its population dropping
!below the threshold.
!The CurrentDets(:,j) is the determinant we are currently spawning from, so this determines the ParentInitiator flag
!which is passed to the SpawnedDets array and refers to whether or not the walkers *parent* is an initiator or not.
!A flag of 0 means the determinant is an initiator, and 1 it is a non-initiator.
        integer, intent(in) :: j, VecSlot, Iter
        integer, intent(out) :: parent_flags
        integer, dimension(lenof_sign) :: CurrentSign
        integer :: part_type
        logical :: tDetinCAS, parent_init

        call extract_sign (CurrentDets(:,j), CurrentSign)

        do part_type=1,lenof_sign
            ! By default, the parent_flags are the flags of the parent...
            parent_init = test_flag (CurrentDets(:,j), flag_parent_initiator(part_type))

            ! The default path through this section makes no changes, leaving the
            ! initiator status of each parent unchanged.
            ! If tAddToInitiator is set, then the state of the parent may change.
            if (tAddToInitiator) then

                if (.not. parent_init) then
                    ! Determinant wasn't previously initiator 
                    ! - want to test if it has now got a large enough population 
                    ! to become an initiator.
                    if (abs(CurrentSign(part_type)) > InitiatorWalkNo) then
                        parent_init = .true.
                        NoAddedInitiators = NoAddedInitiators + 1
                    endif
                elseif (tRetestAddToInit) then
                    ! The source determinant is already an initiator.            
                    ! If tRetestAddToInit is on, the determinants become 
                    ! non-initiators again if their population falls below n_add
                    ! (this is on by default).
                    tDetInCas = .false.
                    if (tTruncCAS) &
                        tDetInCas = TestIfDetInCASBit (CurrentDets(0:NIfD,j))
                   
                    ! If det. in fixed initiator space, or is the HF det, it must
                    ! remain an initiator.
                    if (.not. tDetInCas .and. &
                        .not. (DetBitEQ(CurrentDets(:,j), iLutHF, NIfDBO)) &
                        .and. abs(CurrentSign(part_type)) <= InitiatorWalkNo) then
                        ! Population has fallen too low. Initiator status removed.
                        parent_init = .false.
                        NoAddedInitiators = NoAddedInitiators - 1
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
        INTEGER :: HighPopInNeg(2),HighPopInPos(2),HighPopoutNeg(2),HighPopoutPos(2),TotWalkersNew
        INTEGER, DIMENSION(lenof_sign) :: TempSign

!        WRITE(6,*) 'HighPopPos',HighPopPos
!        WRITE(6,*) 'CurrentSign(HighPopPos)',CurrentSign(HighPopPos)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopNeg),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF
        HighPopInNeg(1)=TempSign(1)
        HighPopInNeg(2)=iProcIndex
        CALL MPIAllReduceDatatype(HighPopinNeg,1,MPI_MINLOC,MPI_2INTEGER,HighPopoutNeg)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopPos),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF
        HighPopInPos(1)=TempSign(1)
        HighPopInPos(2)=iProcIndex
        CALL MPIAllReduceDatatype(HighPopinPos,1,MPI_MAXLOC,MPI_2INTEGER,HighPopoutPos)

!Now, on all processors, HighPopoutPos(1) is the highest positive population, and HighPopoutNeg(1) is the highest negative population.
!HighPopoutPos(2) is the processor the highest population came from.


        IF(iProcIndex.eq.HighPopOutNeg(2)) DetNeg(:)=CurrentDets(:,HighPopNeg)
        IF(iProcIndex.eq.HighPopOutPos(2)) DetPos(:)=CurrentDets(:,HighPopPos)
        CALL MPIBcast(DetNeg,HighPopoutNeg(2))
        CALL MPIBcast(DetPos,HighPopoutPos(2))

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
            tPartFreezeCore .or. tPartFreezeVirt .or. tFixLz .or. tUEG) then
            if (tHPHF .or. tCSF) then
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
        else
            call set_get_spawn_helement (get_helement_det_only)
        endif

        ! Once we have generated the children, do we need to encode them?
        if (.not. (tCSF .or. tHPHF)) then
            call set_encode_child (FindExcitBitDet)
        else
            call set_encode_child (null_encode_child)
        endif

        ! What message should we display for a particle bloom?
        if (tAddToInitiator) then
            bloom_warn_string = '("Particle blooms of more than n_add in &
                                &iteration ", i14, ": A max of ", i8, &
                                &"particles created in one attempt from ")'
        else
            ! Use this variable to store the bloom cutoff level.
            InitiatorWalkNo = 25
            bloom_warn_string = '("Particle blooms of more than 25 in &
                                &iteration ", i14, ": A max of ", i8, &
                                &"particles created in one attempt from ")'
        endif

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
        
    

!This routine is the same as WriteToPopsfilePar, but does not require two main arrays to hold the data.
!The root processors data will be stored in a temporary array while it recieves the data from the other processors.
!This routine will write out to a popsfile. It transfers all walkers to the head node sequentially, so does not want to be called too often
    SUBROUTINE WriteToPopsfileParOneArr(Dets,nDets)
        use util_mod, only: get_unique_filename, get_free_unit
        use CalcData, only: iPopsFileNoWrite
        use Logging, only: tIncrementPops
        use constants, only: size_n_int,MpiDetInt,n_int
        integer(int64),intent(in) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(in) :: Dets(0:nIfTot,1:nDets)
        INTEGER :: error
        integer(int64) :: WalkersonNodes(0:nProcessors-1)
        INTEGER :: Tag,Total,i,j,k
        INTEGER(KIND=n_int), ALLOCATABLE :: Parts(:,:)
        INTEGER :: PartsTag=0
        INTEGER :: nMaxDets
        integer :: iunit
        CHARACTER(len=*) , PARAMETER :: this_routine='WriteToPopsfileParOneArr'
        character(255) :: popsfile
        INTEGER, DIMENSION(lenof_sign) :: TempSign


        IF(lenof_sign.ne.1) THEN
            WRITE(6,*) "Cannot write complex walkers out to POPSFILE yet..."
        ENDIF

        CALL MPIBarrier(error)  !sync

!First, make sure we have up-to-date information - again collect AllTotWalkers,AllSumNoatHF and AllSumENum...
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_Sum,root,MPI_COMM_WORLD,error)    
!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        CALL MPISumRoot(SumNoatHF,1,AllSumNoatHF,Root)
        CALL MPISumRoot(SumENum,1,AllSumENum,Root)

!We also need to tell the root processor how many particles to expect from each node - these are gathered into WalkersonNodes
        CALL MPIAllGather(nDets,WalkersonNodes,error)

        Tag=125
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)

        IF(iProcIndex.eq.root) THEN
!First, check that we are going to receive the correct number of particles...
            Total=0
            do i=0,nProcessors-1
                Total=Total+INT(WalkersonNodes(i)/iPopsPartEvery)
            enddo
            AllTotWalkers=REAL(Total,dp)
!            IF(Total.ne.AllTotWalkers) THEN
!                CALL Stop_All("WriteToPopsfilePar","Not all walkers accounted for...")
!            ENDIF

!Write header information
#ifdef __INT64
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a 64-bit binary reduced POPSFILEBIN, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a 64-bit reduced POPSFILE, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to 64-bit binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to 64-bit POPSFILE..."
                ENDIF
            ENDIF
#else
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a binary reduced POPSFILEBIN, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a reduced POPSFILE, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to POPSFILE..."
                ENDIF
            ENDIF
#endif
            IF(tBinPops) THEN
                call get_unique_filename('POPSFILEHEAD',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ELSE
                call get_unique_filename('POPSFILE',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ENDIF
            iunit = get_free_unit()
            OPEN(iunit,FILE=popsfile,Status='replace')
            WRITE(iunit,"(A)") "# POPSFILE VERSION 2"
#ifdef __INT64
            WRITE(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') '64BitDets=',.TRUE.,'HPHF=',tHPHF,'Lz=',tFixLz,'Initiator=',tTruncInitiator
#else
            WRITE(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') '64BitDets=',.FALSE.,'HPHF=',tHPHF,'Lz=',tFixLz,'Initiator=',tTruncInitiator
#endif
            WRITE(iunit,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
            WRITE(iunit,*) DiagSft,"   DIAG SHIFT"
            WRITE(iunit,*) AllSumNoatHF,"   SUMNOATHF (all nodes)"
            WRITE(iunit,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
            WRITE(iunit,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
            IF(tBinPops) THEN
                CLOSE(iunit)
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
                OPEN(iunit,FILE=popsfile,Status='replace',form='unformatted')
            ENDIF

            IF(tBinPops) THEN
                do j=1,nDets
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        call extract_sign(Dets(:,j),TempSign)
                        WRITE(iunit) Dets(0:NIfDBO,j),TempSign(:)
                    ENDIF
                enddo
            ELSE
                do j=1,nDets
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        do k=0,NIfDBO
                            WRITE(iunit,"(I24)",advance='no') Dets(k,j)
                        enddo
                        call extract_sign(Dets(:,j),TempSign)
                        WRITE(iunit,*) TempSign(:)
                    ENDIF
                enddo
            ENDIF
!            WRITE(6,*) "Written out own walkers..."
!            CALL FLUSH(6)

!Now, we copy the head nodes data to a new array...
            nMaxDets=maxval(WalkersOnNodes)
            ALLOCATE(Parts(0:NIfTot,nMaxDets),stat=error)
            CALL LogMemAlloc('Parts',int(nMaxDets,int32)*(NIfTot+1),size_n_int,this_routine,PartsTag,error)

!Now we need to receive the data from each other processor sequentially
!We can overwrite the head nodes information, since we have now stored it elsewhere.
            do i=1,nProcessors-1
!Run through all other processors...receive the data...
                j=WalkersonNodes(i)*(NIfTot+1)
                CALL MPIRecv(Parts(0:NIfTot,1:WalkersonNodes(i)),j,i,Tag,error)
!                WRITE(6,*) "Recieved walkers for processor ",i
!                CALL FLUSH(6)
                
!Then write it out...
                IF(tBinPops) THEN
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            call extract_sign(Parts(:,j),TempSign)
                            WRITE(iunit) Parts(0:NIfDBO,j),TempSign(:)
                        ENDIF
                    enddo
                ELSE
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            do k=0,NIfDBO
                                WRITE(iunit,"(I24)",advance='no') Parts(k,j)
                            enddo
                            call extract_sign(Parts(:,j),TempSign)
                            WRITE(iunit,*) TempSign(:)
                        ENDIF
                    enddo
                ENDIF
!                WRITE(6,*) "Writted out walkers for processor ",i
!                CALL FLUSH(6)

            enddo

            CLOSE(iunit)

!Deallocate memory for temporary storage of information.
            DEALLOCATE(Parts)
            CALL LogMemDealloc(this_routine,PartsTag)

        ELSE
!All other processors need to send their data to root...
            j=nDets*(NIfTot+1)
            CALL MPISend(Dets(0:NIfTot,1:nDets),j,root,Tag,error)
!            WRITE(6,*) "Have sent info to head node..."
!            CALL FLUSH(6)
        ENDIF

!Reset the values of the global variables
        AllSumNoatHF = 0
        AllSumENum=0.D0
        AllTotWalkers=0.D0

        RETURN

    END SUBROUTINE WriteToPopsfileParOneArr


!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        use util_mod, only: get_unique_filename, get_free_unit
        use CalcData, only: iPopsFileNoRead
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,iWeightPopRead
        use Logging, only: tIncrementPops,tZeroProjE
        use constants, only: size_n_int,bits_n_int
        LOGICAL :: exists,tBinRead
        INTEGER :: AvWalkers,WalkerstoReceive(nProcessors)
        INTEGER*8 :: NodeSumNoatHF(nProcessors)
        integer(int64) :: TempTotParts(lenof_sign),TempCurrWalkers
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        INTEGER*8 :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        REAL*8 :: r,FracPart,Gap,DiagSftTemp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        integer :: sign_largest

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile does not work with complex walkers")
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) for particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"



!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters are the same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) AllTotWalkers
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers," configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.D0) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.D0
                AllSumNoatHF = 0
            ENDIF

!Need to calculate the number of walkers each node will receive...
            AvWalkers=NINT(AllTotWalkers/real(nProcessors,dp))

!Divide up the walkers to receive for each node
            do i=1,nProcessors-1
                WalkerstoReceive(i)=AvWalkers
            enddo
!The last processor takes the 'remainder'
            WalkerstoReceive(nProcessors)=AllTotWalkers-(AvWalkers*(nProcessors-1))

!Quick check to ensure we have all walkers accounted for
            total=0
            do i=1,nProcessors
                total=total+WalkerstoReceive(i)
            enddo
            if (total /= AllTotWalkers) then
                CALL Stop_All("ReadFromPopsfilePar","All Walkers not accounted for when reading in from POPSFILE")
            endif
            
!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = AllSumNoatHF(1)/nProcessors !This is the average Sumnoathf
            do i=1,nProcessors-1
                NodeSumNoatHF(i)=INT(AvSumNoatHF,int64)
            enddo
            NodeSumNoatHF(nProcessors)=AllSumNoatHF(1)-INT((AvSumNoatHF*(nProcessors-1)),int64)

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)
                
!Reset the global variables
            AllSumENum=0.D0
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft,root)
        CALL MPIBCast(SumENum,root)
        CALL MPIBCast(InitWalkers,root)
        CALL MPIBCast(NEquilSteps,root)
        CALL MPIBCast(NShiftEquilSteps,root)
        CALL MPIBCast(TSinglePartPhase,root)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and the SumNoatHF for each node which is distributed approximatly equally
        CALL MPIScatter(WalkerstoReceive,TempInitWalkers,root,error)
        CALL MPIScatter(NodeSumNoatHF,SumNoatHF(1),root,error)

        IF(MemoryFacPart.le.1.D0) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF
        
!Now we want to allocate memory on all nodes.
        MaxWalkersPart=NINT(MemoryFacPart*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
        MaxSpawned=NINT(MemoryFacSpawn*(NINT(InitWalkers*ScaleWalkers)))

        Gap=REAL(MaxSpawned)/REAL(nProcessors)
        do i=0,nProcessors-1
            InitialSpawnedSlots(i)=NINT(Gap*i)+1
        enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

        CALL MPIBarrier(error)  !Sync

!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecDets array.')
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(:,:)=0
        MemoryAlloc=(NIfTot+1)*MaxWalkersPart*size_n_int    !Memory Allocated in bytes

!Just allocating this here, so that the SpawnParts arrays can be used for sorting the determinants when using direct annihilation.
        WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
        ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVecTag,ierr)
        SpawnVec(:,:)=0
        ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
        SpawnVec2(:,:)=0
        MemoryAlloc=MemoryAlloc+MaxSpawned*(1+NIfTot)*2*size_n_int

!Point at correct spawning arrays
        SpawnedParts=>SpawnVec
        SpawnedParts2=>SpawnVec2
!Allocate pointer to the correct walker array...
        CurrentDets=>WalkVecDets

!Need to now allocate other arrays
        IF(.not.tRegenDiagHEls) THEN
            ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
            WalkVecH(:)=0.d0
            MemoryAlloc=MemoryAlloc+8*MaxWalkersPart
        ELSE
            WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*8,dp)/1048576.D0," Mb/Processor"
        ENDIF

        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
        ENDIF

! The hashing will be different in the new calculation from the one where the POPSFILE was produced, this means we must recalculate the processor each determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing them in the same way as the spawning steps are done - by finding the determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        do i=1,AllTotWalkers
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif
        
#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc=DetermineDetProc(iLutTemp)   !This wants to return a value between 0 -> nProcessors-1
            IF((Proc.eq.iProcIndex).and.(abs(TempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                call encode_bit_rep(CurrentDets(:,CurrWalkers),iLutTemp(0:NIfDBO),TempSign,0)   !Do not need to send a flag here...
                                                                                                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(TempSign(1)) > sign_largest) then
                sign_largest = abs(TempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=REAL(CurrWalkers,dp)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort (currentdets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(CurrentDets(:,i),CurrentDets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, currentdets(:, i-1)
                print*, currentdets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                call extract_sign(CurrentDets(:,l),TempSign)
                TempSign=TempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(TempSign(1).lt.0) THEN
                        TempSign(1)=TempSign(1)-1
                    ELSE
                        TempSign(1)=TempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(CurrentDets(:,l),TempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                aLlTotWalkersOld=AllTotWalkers
                iter_data_fciqmc%tot_parts_old = AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                iter_data_fciqmc%tot_parts_old = AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL FLUSH(6)

        ! If we are changing the reference determinant to the largest
        ! weighted one in the file, do it here
        if (tReadPopsChangeRef .or. tReadPopsRestart) then
            if (.not. DetBitEq(ilut_largest, iLutRef, NIfDBO)) then
                
                ! Set new reference
                iLutRef = ilut_largest
                call decode_bit_det (ProjEDet, iLutRef)
                tNoBrillouin = .true.

                ! Recalculate the reference E
                if (tHPHF) then
                    HElemTemp = hphf_diag_helement (ProjEDet, iLutRef)
                else
                    HElemTemp = get_helement (ProjEDet, ProjEDet, 0)
                endif
                Hii = real(HElemTemp, dp)

                ! Output info on root node.
                if (iProcIndex == root) then
                    write(6, '(a)', advance='no') &
                        "Changing projected energy reference determinant to: "
                    call write_det (6, ProjEDet, .true.)
                    write (6, '(a)') &
                        "Ensuring that Brillouin's theorem is no longer used."
                    write (6, '(a,g25.15)') &
                        "Reference energy now set to: ", Hii
                endif
            endif
        endif

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,TotWalkers
            call decode_bit_det (TempnI, currentDets(:,j))
            Excitlevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
            IF(Excitlevel.eq.0) THEN
                IF(.not.tRegenDiagHEls) CurrentH(j)=0.D0
            ELSE
                IF(.not.tRegenDiagHEls) THEN
                    if (tHPHF) then
                        HElemTemp = hphf_diag_helement (TempnI, &
                                                        CurrentDets(:,j))
                    else
                        HElemTemp = get_helement (TempnI, TempnI, 0)
                    endif
                    CurrentH(j)=REAL(HElemTemp,dp)-Hii
                ENDIF

            ENDIF
            call extract_sign(CurrentDets(:,j),TempSign)
            TotParts=TotParts+abs(TempSign(1))

        enddo

        TempTotParts=REAL(TotParts,dp)

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TempTotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfilePar

!This routine reads in particle configurations from a POPSFILE.
! It's a bastardisation of ReadFromPopsfilePar but only reads a popsfile into an array and leaves the rest until later.
! nDets comes in with the max number of amplitudes possible in Dets.
! Dets is then filled with all the amplitudes in the POPSFILE and nDets is returned with that number.

! This will very likely still be tweaking the inner heart-strings of FciMCPar.  Caveatis stulti!  AJWT
! GHB says he will incorporate this functionality into a rewrite of ReadFromPopsfilePar. 19/8/2010
    SUBROUTINE ReadFromPopsfileOnly(Dets,nDets)
        use util_mod, only: get_unique_filename, get_free_unit
        use CalcData, only: iPopsFileNoRead
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,iWeightPopRead
        use Logging, only: tIncrementPops,tZeroProjE
        use constants, only: size_n_int,bits_n_int
        integer(int64),intent(inout) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(out) :: Dets(0:nIfTot,1:nDets)
        LOGICAL :: exists,tBinRead
        integer(int64) :: TempTotParts(lenof_sign),TempCurrWalkers
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        INTEGER*8 :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        REAL*8 :: r,FracPart,Gap,DiagSftTemp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        integer :: sign_largest

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile does not work with complex walkers")
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) for particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"



!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters are the same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) AllTotWalkers
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers," configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.D0) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.D0
                AllSumNoatHF = 0
            ENDIF

!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = AllSumNoatHF(1)/nProcessors !This is the average Sumnoathf

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)
                
!Reset the global variables
            AllSumENum=0.D0
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft,root)
        CALL MPIBCast(SumENum,root)
        CALL MPIBCast(InitWalkers,root)
        CALL MPIBCast(NEquilSteps,root)
        CALL MPIBCast(NShiftEquilSteps,root)
        CALL MPIBCast(TSinglePartPhase,root)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and the SumNoatHF for each node which is distributed approximatly equally

        IF(MemoryFacPart.le.1.D0) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF
        
        CALL MPIBarrier(error)  !Sync

        if(AllTotWalkers>nDets) CALL Stop_All(this_routine,'Not enough memory to read in POPSFILE.')

! The hashing will be different in the new calculation from the one where the POPSFILE was produced, this means we must recalculate the processor each determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing them in the same way as the spawning steps are done - by finding the determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        write(6,*) "Reading in ", AllTotWalkers, " walkers"
        do i=1,AllTotWalkers
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif
        
#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc=0  !DetermineDetProc(iLutTemp)   !This wants to return a value between 0 -> nProcessors-1
            IF((Proc.eq.iProcIndex).and.(abs(TempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                call encode_bit_rep(Dets(:,CurrWalkers),iLutTemp(0:NIfDBO),TempSign,0)   !Do not need to send a flag here...
                                                                                                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(TempSign(1)) > sign_largest) then
                sign_largest = abs(TempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=REAL(CurrWalkers,dp)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort(Dets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(Dets(:,i),Dets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, dets(:, i-1)
                print*, dets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                call extract_sign(Dets(:,l),TempSign)
                TempSign=TempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(TempSign(1).lt.0) THEN
                        TempSign(1)=TempSign(1)-1
                    ELSE
                        TempSign(1)=TempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(Dets(:,l),TempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL FLUSH(6)

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,nDets
            call extract_sign(Dets(:,j),TempSign)
            TotParts=TotParts+abs(TempSign(1))
        enddo

        TempTotParts=REAL(TotParts,dp)

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TempTotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfileOnly

    function attempt_create_trunc_spawn (get_spawn_helement, DetCurr,&
                                         iLutCurr, wSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type) &
                                         result(child)
        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        integer, dimension(lenof_sign), intent(in) :: wSign
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
                               tParity, walkExcitLevel, part_type)
        else
            child = 0
        endif
    end function

    function attempt_create_trunc_spawn_encode (get_spawn_helement, DetCurr,&
                                         iLutCurr, wSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type) &
                                         result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        integer, dimension(lenof_sign), intent(in) :: wSign
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
                               tParity, walkExcitLevel, part_type)
        else
            child = 0
        endif
    end function

    function attempt_create_normal (get_spawn_helement, DetCurr, iLutCurr, &
                                    wSign, nJ, iLutnJ, prob, HElGen, ic, ex, tparity,&
                                    walkExcitLevel, part_type) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! 1 = Real parent particle, 2 = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        integer, dimension(lenof_sign), intent(in) :: wSign
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
        
        real(dp) :: rat, r
        integer :: extracreate, iUnused
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

#endif
        ! Avoid compiler warnings
        iUnused = walkExcitLevel

    end function

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,iLutCurr,WSign,nJ,iLutnJ,Prob,IC,Ex,tParity)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        use Logging, only : CCMCDebug
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,ExtraCreate,Ex(2,2),Bin
        INTEGER(KIND=n_int) :: iLutCurr(0:NIfTot),iLutnJ(0:NIfTot)
        LOGICAL :: tParity
        REAL*8 :: Prob,r,rat
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

                rat=Tau/abs(Prob)

                rh=Prob ! to get the signs right for later on.
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
                CALL Stop_All("AttemptCreatePar","Trying to histogram off-diagonal matrix elements, but outside histogram array bounds.")
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
                CALL Stop_All("AttemptCreatePar","Histogramming energies higher than the arrays can cope with. Increase iNoBins or BinRange")
            ENDIF
            IF(AttemptCreatePar.ne.0) THEN
                SpawnHist(Bin)=SpawnHist(Bin)+real(abs(AttemptCreatePar),dp)
!                WRITE(6,*) "Get Here!", real(abs(AttemptCreatePar),dp),Bin
            ENDIF
            AttemptHist(Bin)=AttemptHist(Bin)+(Tau/Prob)
        ENDIF

    END FUNCTION AttemptCreatePar

    subroutine walker_death (attempt_die, iter_data, DetCurr, iLutCurr, Kii, &
                             wSign, VecSlot)
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
        integer, dimension(lenof_sign), intent(in) :: wSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        integer, intent(inout) :: VecSlot
        real(dp), intent(in) :: Kii
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, dimension(lenof_sign) :: iDie
        integer, dimension(lenof_sign) :: CopySign

        ! Do particles on determinant die? iDie can be both +ve (deaths), or
        ! -ve (births, if shift > 0)
        iDie = attempt_die (DetCurr, Kii, wSign)

!        IF(iDie.ne.0) WRITE(6,*) "Death: ",iDie
        
        IFDEBUGTHEN(FCIMCDebug,3) 
            if(sum(abs(iDie)).ne.0) write(6,"(A,2I4)") "Death: ",iDie(:)
        ENDIFDEBUG

        ! Update death counter
        iter_data%ndied = iter_data%ndied + min(iDie, abs(wSign))
        NoDied = NoDied + sum(min(iDie, abs(wSign)))

        ! Count any antiparticles
        iter_data%nborn = iter_data%nborn + max(iDie - abs(wSign), 0)
        NoBorn = NoBorn + sum(max(iDie - abs(wSign), 0))

        ! Calculate new number of signed particles on the det.
        CopySign = wSign - (iDie * sign(1, wSign))

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

        real(dp) :: r, rat, fac
        integer :: i, iUnused

        fac = tau * (Kii-DiagSft)

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
        REAL*8 :: HList(Length)
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




    SUBROUTINE WriteHistogramEnergies()
        use util_mod, only: get_free_unit
        INTEGER :: i, io(8)
        REAL*8 :: Norm,EnergyBin

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
        REAL*8 :: norm,norm1

        CALL MPISum(Histogram,AllHistogram)
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
        INTEGER :: i,IterRead, io1, io2, io3
        REAL*8 :: norm,norm1,norm2,norm3,ShiftRead,AllERead,NumParts
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
            enddo

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

    subroutine new_child_stats_hist_hamil (iter_data, iLutI, iLutJ, ic, &
                                           walkExLevel, child)
        ! Based on old AddHistHamilEl. Histograms the hamiltonian matrix, and 
        ! then calls the normal statistics routine.

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel
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
        call new_child_stats_normal (iter_data, iLutI, iLutJ, ic, &
                                     walkExLevel, child)

    end subroutine


    ! TODO: COMMENTING


    subroutine iter_diagnostics ()

        character(*), parameter :: this_routine = 'iter_diagnostics'
        integer(int64) :: walkers_diff
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
                write(6,*) AllTotWalkers, TotWalkers
                call stop_all (this_routine, 'All particles have died. &
                              &Consider choosing new seed, or raising shift &
                              &value.')
            endif

! AJWT dislikes doing this type of if based on a (seeminly unrelated) input option, but can't see another easy way.
!  TODO:  Something to make it better

            if(.not.tCCMC) then
               ! Check how balanced the load on each processor is (even though
               ! we cannot load balance with direct annihilation).
               walkers_diff = MaxWalkersProc - MinWalkersProc
               mean_walkers = AllTotWalkers / real(nProcessors,dp)
               if (walkers_diff > nint(mean_walkers / 10.d0) .and. &
                   sum(AllTotParts) > real(nProcessors * 500, dp)) then
                   root_write (6, '(a, f20.10, 2i12)') &
                       'Number of determinants assigned to each processor &
                       &unbalanced: ', (walkers_diff * 10.d0) / &
                       real(mean_walkers), MinWalkersProc, MaxWalkersProc
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

    end subroutine

    subroutine population_check ()
        
        integer :: int_tmp(2), pop_highest, proc_highest, pop_change
        integer :: det(nel), i
        HElement_t :: h_tmp

        if (tCheckHighestPop) then

            ! Obtain the determinant (and its processor) with the highest
            ! population.
            call MPIAllReduceDatatype ((/iHighestPop, iProcIndex/), 1, &
                                       MPI_MAXLOC, MPI_2INTEGER, int_tmp)
            pop_highest = int_tmp(1)
            proc_highest = int_tmp(2)

            ! How many walkers do we need to switch dets?
            pop_change = int(FracLargerDet * real(abs_int_sign(AllNoAtHF), dp))
            if (pop_change < pop_highest .and. sum(AllTotParts) > 10000) then

                ! Write out info!
                if (tHPHF) then
                    root_print 'Highest weighted CLOSED-SHELL determinant not&
                               & reference det: ', pop_highest, abs_int_sign(AllNoAtHF)
                else
                    root_print 'Highest weighted determinant not reference &
                               &det: ', pop_highest, abs_int_sign(AllNoAtHF)
                endif

                ! Are we changing the reference determinant?
                if (tChangeProjEDet) then
                    ! Communicate the change to all dets and print out.
                    call MPIBcast (HighestPopDet(0:NIfTot), proc_highest)
                    iLutRef = HighestPopDet
                    call decode_bit_det (ProjEDet, iLutRef)
                    write (6, '(a)', advance='no') 'Changing projected &
                                  &energy reference & &determinant to: '
                    call write_det (6, ProjEDet, .true.)

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
                    write (6, '(a, g25.15)') 'Reference energy now set to: ',&
                                             Hii

                    ! Reset averages
                    SumENum = 0
                    SumNoatHF = 0
                    HFPopCyc = 0
                    ProjEIterSum = 0
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
                    call MPIBcast (HighestPopDet, proc_highest)
                    iLutRef = HighestPopDet
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

        integer :: int_tmp(7+lenof_sign)
        HElement_t :: real_tmp(2)
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
        IF(lenof_sign.eq.1) THEN
            AllHFCyc = real(int_tmp(5), dp)
            AllSpawnFromSing = int_tmp(6)
            iter_data%update_growth_tot = int_tmp(7:6+lenof_sign)
        ELSE
            AllHFCyc = cmplx(int_tmp(5),int_tmp(6), dp)
            AllSpawnFromSing = int_tmp(7)
            iter_data%update_growth_tot = int_tmp(8:7+lenof_sign)
        ENDIF

        ! Integer summations required for the initiator method
        if (tTruncInitiator) then
            call MPIReduce ((/NoAborted, NoAddedInitiators, NoInitDets, &
                              NoNonInitDets, NoInitWalk, NoNonInitWalk, &
                              NoExtraInitdoubs, InitRemoved/), &
                            MPI_SUM, int64_tmp)
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
        call MPIReduce ((/TotWalkers, TotParts, SumNoatHF, tot_parts_new/), &
                        MPI_SUM, int64_tmp(1:1+3*lenof_sign))
        AllTotWalkers = int64_tmp(1)
        AllTotParts = int64_tmp(2:1+lenof_sign)
        AllSumNoatHF = int64_tmp(2+lenof_sign:1+2*lenof_sign)
        tot_parts_new_all = int64_tmp(2+2*lenof_sign:1+3*lenof_sign)

        ! real(dp) values (Calculates the energy by summing all on HF and 
        ! doubles)
        call MPIReduce ((/ENumCyc, SumENum/), MPI_SUM, real_tmp)
        AllENumCyc = real_tmp(1)
        AllSumENum = real_tmp(2)

        ! Max/Min values (check load balancing)
        call MPIReduce (TotWalkers, MPI_MAX, MaxWalkersProc)
        call MPIReduce (TotWalkers, MPI_MIN, MinWalkersProc)

        ! We need the total number on the HF and SumWalkersCyc to be valid on
        ! ALL processors (n.b. one of these is 32bit, the other 64)
        call MPISum (NoatHF, AllNoatHF)
        call MPISum (SumWalkersCyc, AllSumWalkersCyc)

!        WRITE(6,*) "***",iter_data%update_growth_tot,AllTotParts-AllTotPartsOld
        ASSERT(iter_data%update_growth_tot(1).eq.AllTotParts(1)-AllTotPartsOld(1))
!        ASSERT(iter_data%update_growth_tot(2).eq.AllTotParts(2)-AllTotPartsOld(2))
        
    end subroutine

    subroutine update_shift (iter_data)

        type(fcimc_iter_data), intent(in) :: iter_data
        integer(int64) :: tot_walkers
        logical :: tReZeroShift
        real(dp) :: AllGrowRateRe, AllGrowRateIm

        integer :: error

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

            ! Initiator abort growth rate
            if (tTruncInitiator) then
                AllGrowRateAbort = (sum(iter_data%update_growth_tot + &
                                    iter_data%tot_parts_old) + AllNoAborted) &
                                    / (sum(iter_data%tot_parts_old) &
                                       + AllNoAbortedOld)
            endif

            ! Exit the single particle phase if the number of walkers exceeds
            ! the value in the input file. If particle no has fallen, re-enter
            ! it.
            tReZeroShift = .false.
            if (TSinglePartPhase) then
! AJWT dislikes doing this type of if based on a (seeminly unrelated) input option, but can't see another easy way.
!  TODO:  Something to make it better
                if(.not.tCCMC) then
                    tot_walkers = int(InitWalkers, int64) * int(nProcessors,int64)
                else
                    tot_walkers = int(InitWalkers, int64)
                endif
                if ( (sum(AllTotParts) > tot_walkers) .or. &
                     (abs_int_sign(AllNoatHF) > MaxNoatHF)) then
!                     WRITE(6,*) "AllTotParts: ",AllTotParts(1),AllTotParts(2),tot_walkers
                    write (6, *) 'Exiting the single particle growth phase - &
                                 &shift can now change'
                    VaryShiftIter = Iter
                    tSinglePartPhase = .false.
                endif
            elseif (abs_int_sign(AllNoatHF) < (MaxNoatHF - HFPopThresh)) then
                write (6, *) 'No at HF has fallen too low - reentering the &
                             &single particle growth phase - particle number &
                             &may grow again.'
                tSinglePartPhase = .true.
                tReZeroShift = .true.
            endif

            ! How should the shift change for the entire ensemble of walkers 
            ! over all processors.
            if (.not. tSinglePartPhase) then
                
                ! New shift value
                DiagSft = DiagSft - (log(AllGrowRate) * SftDamp) / &
                                    (Tau * StepsSft)

                if (lenof_sign == 2) then
                    DiagSftRe = DiagSftRe - (log(AllGrowRateRe) * SftDamp) / &
                                            (Tau * StepsSft)
                    DiagSftIm = DiagSftIm - (log(AllGrowRateIm) * SftDamp) / &
                                            (Tau * StepsSft)
                endif

                ! Update the shift averages
                if ((iter - VaryShiftIter) >= nShiftEquilSteps) then
                    if ((iter-VaryShiftIter-nShiftEquilSteps) < StepsSft) &
                        write (6, *) 'Beginning to average shift value.'
                    VaryShiftCycles = VaryShiftCycles + 1
                    SumDiagSft = SumDiagSft + DiagSft
                    AvDiagSft = SumDiagSft / real(VaryShiftCycles, dp)
                endif

                ! Update DiagSftAbort for initiator algorithm
                if (tTruncInitiator) then
                    DiagSftAbort = DiagSftAbort - &
                              (log(real(AllGrowRateAbort, dp)) * SftDamp) / &
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


            ! AllSumNoatHF can be 0 if equilsteps is on.
#ifdef __CMPLX
            if (AllSumNoatHF(1).ne.0.or.AllSumNoatHF(2).ne.0) then
                ProjectionE = AllSumENum / CMPLX(AllSumNoatHF(1), AllSumNoatHF(2), dp) 
            endif
#else
            if (AllSumNoatHF(1) /= 0) ProjectionE = AllSumENum / AllSumNoatHF(1)
#endif

            ! Calculate the projected energy where each update cycle 
            ! contributes the same weight to the average for its estimator 
            ! for the energy.
            if (abs(AllHFCyc) /= 0.D0) then
                ProjEItersum = ProjEIterSum + (AllENumCyc / AllHFCyc)
                ! Count the number of interactions where we have a non-zero
                ! contribution from HF particles
                HFPopCyc = HFPopCyc + 1
                ProjEIter = ProjEIterSum / real(HFPopCyc, dp)
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
        call MPIBcast (tSinglePartPhase, Root)
        call MPIBcast (VaryShiftIter, Root)
        call MPIBcast (DiagSft, Root)

    end subroutine


    subroutine rezero_iter_stats (iter_data, tot_parts_new_all)
        
        type(fcimc_iter_data) :: iter_data
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

        ! We don't want to rezero projEIter, as PrintFCIMCStats will print
        ! out zero on next iteration if NoAtHF = 0.
        ! ProjEIter = 0

        HFCyc = 0
        NoAborted = 0
        NoInitDets = 0
        NoNonInitDets = 0
        NoInitWalk = 0
        NoNonInitWalk = 0
        InitRemoved = 0

        ! Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld = TotWalkers
        TotPartsOld = TotParts

        ! Save the number at HF to use in the HFShift
        OldAllNoatHF = AllNoatHF

        ! Also reinitialise the global variables
        !  --> should not necessarily need to do this...
        ! TODO: Remove these zeroes?
        AllSumENum = 0
        AllSumNoatHF = 0
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        AllNoAbortedOld = AllNoAborted
        AllTotWalkers = 0
        AllTotParts = 0
        AllGrowRate = 0
        AllSumWalkersCyc = 0
        AllAnnihilated = 0
        AllNoatHF = 0
        AllNoatDoubs = 0
        AllNoBorn = 0
        AllSpawnFromSing = 0
        AllNoDied = 0
        AllNoAborted = 0
        AllNoAddedInitiators = 0
        AllNoInitDets = 0
        AllNoNonInitDets = 0
        AllNoInitWalk = 0
        AllNoNonInitWalk = 0
        AllNoExtraInitDoubs = 0
        AllInitRemoved = 0

        ! Reset the counters
        iter_data%update_growth = 0
        iter_data%update_iters = 0
        iter_data%tot_parts_old = tot_parts_new_all

    end subroutine

    subroutine calculate_new_shift_wrapper (iter_data, tot_parts_new)

        type(fcimc_iter_data) :: iter_data
        integer(int64), dimension(lenof_sign), intent(in) :: tot_parts_new
        integer(int64), dimension(lenof_sign) :: tot_parts_new_all

        call collate_iter_data (iter_data, tot_parts_new, tot_parts_new_all)
        call iter_diagnostics ()
        call population_check ()
        call update_shift (iter_data)
        call WriteFCIMCStats ()
        call rezero_iter_stats (iter_data, tot_parts_new_all)

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
                WRITE(initiatorstats_unit,"(A2,A10,2A16,2A16,1A20,6A18)") "# ","1.Step","2.No Aborted","3.NoAddedtoInit","4.FracDetsInit","5.FracWalksInit","6.InstAbortShift",&
&               "7.NoInit","8.NoNonInit"
            ENDIF
            IF(tLogComplexPops) THEN
                WRITE(complexstats_unit,"(A)") '#   1.Step  2.Shift     3.RealShift     4.ImShift   5.TotParts      6.RealTotParts      7.ImTotParts'
            ENDIF

#ifdef __CMPLX
            WRITE(6,"(A)") "       Step     Shift      WalkerCng(Re)  WalkerCng(Im)    TotWalkers(Re)   TotWalkers(Im)    Proj.E(Re)   ProjE(Im)     Proj.E.ThisCyc(Re)  Proj.E.ThisCyc(Im)"&
&             //"   NoatHF(Re)   NoatHF(Im)   NoatDoubs      AccRat     UniqueDets     IterTime"
            WRITE(fcimcstats_unit,"(A,I4,A,L,A,L,A,L)") "# FCIMCStats VERSION 2 - COMPLEX : NEl=",NEl," HPHF=",tHPHF,' Lz=',tFixLz,' Initiator=',tTruncInitiator
            WRITE(fcimcstats_unit,"(A)") "#     1.Step   2.Shift    3.WalkerCng(Re)  4.WalkerCng(Im)   5.TotWalkers(Re)  6.TotWalkers(Im)  7.Proj.E(Re)   8.Proj.E(Im)  "&
&           // " 9.Proj.E.ThisCyc(Re)  10.Proj.E.ThisCyc(Im)  11.NoatHF(Re)   12.NoatHF(Im)  13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime 17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  "&
&           // " 20.HFInstShift  21.TotInstShift  22.Tot-Proj.E.ThisCyc(Re)  23.HFContribtoE(Re)  24.HFContribtoE(Im)   25.NumContribtoE(Re)  26.NumContribtoE(Im)"
#else
            WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          Av.Shift     Proj.E.ThisCyc   NoatHF NoatDoubs      AccRat     UniqueDets     IterTime"
            WRITE(fcimcstats_unit,"(A,I4,A,L,A,L,A,L)") "# FCIMCStats VERSION 2 - REAL : NEl=",NEl," HPHF=",tHPHF,' Lz=',tFixLz,' Initiator=',tTruncInitiator
            WRITE(fcimcstats_unit,"(A)") "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  9.Proj.E       10.Av.Shift"&
&           // " 11.Proj.E.ThisCyc  12.NoatHF 13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime 17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  20.ProjE.ThisIter "&
&           // " 21.HFInstShift  22.TotInstShift  23.Tot-Proj.E.ThisCyc   24.HFContribtoE  25.NumContribtoE"
#endif
            
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

    subroutine WriteFCIMCStats()

        if (iProcIndex == root) then
#ifdef __CMPLX
            write(fcimcstats_unit,"(I12,G16.7,2I10,2I12,4G17.9,3I10,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,7G17.9)") &
                Iter + PreviousCycles, DiagSft, &
                AllTotParts(1)-AllTotPartsOld(1), AllTotParts(2)-AllTotPartsOld(2), &
                AllTotParts(1),AllTotParts(2), &
                REAL(ProjectionE,dp),AIMAG(ProjectionE), REAL(AllENumCyc / AllHFCyc,dp),AIMAG(AllENumCyc / AllHFCyc), &
                AllNoatHF(1),AllNoatHF(2), &
                AllNoatDoubs, AccRat, AllTotWalkers, IterTime, &
                real(AllSpawnFromSing) / real(AllNoBorn), WalkersDiffProc, &
                TotImagTime, HFShift, InstShift, &
                REAL(AllENumCyc / AllHFCyc,dp) + Hii, REAL(AllHFCyc / StepsSft,dp), AIMAG(AllHFCyc / StepsSft), &
                REAL(AllENumCyc / StepsSft,dp),AIMAG(AllENumCyc / StepsSft)
            write (6, "(I12,G16.7,2I10,2I12,4G17.9,3I10,G13.5,I12,&
                      &G13.5)") Iter + PreviousCycles, DiagSft, &
                AllTotParts(1)-AllTotPartsOld(1),AllTotParts(2)-AllTotPartsOld(2), &
                AllTotParts(1),AllTotParts(2), &
                real(ProjectionE,dp),aimag(ProjectionE), real(AllENumCyc / AllHFCyc,dp),aimag(AllENumCyc / AllHFCyc), &
                AllNoatHF(1),AllNoatHF(2), &
                AllNoatDoubs, AccRat, AllTotWalkers, IterTime
#else
            write(fcimcstats_unit,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,6G17.9)") &
                Iter + PreviousCycles, DiagSft, &
                sum(AllTotParts) - sum(AllTotPartsOld), AllGrowRate, &
                sum(AllTotParts), AllAnnihilated, AllNoDied, AllNoBorn, &
                ProjectionE, AvDiagSft, AllENumCyc / AllHFCyc, AllNoatHF, &
                AllNoatDoubs, AccRat, AllTotWalkers, IterTime, &
                real(AllSpawnFromSing) / real(AllNoBorn), WalkersDiffProc, &
                TotImagTime, 0.D0, HFShift, InstShift, &
                AllENumCyc / AllHFCyc + Hii, AllHFCyc / StepsSft, &
                AllENumCyc / StepsSft
            write (6, "(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,&
                      &G13.5)") Iter + PreviousCycles, DiagSft, &
                sum(AllTotParts) - sum(AllTotPartsOld), AllGrowRate, &
                sum(AllTotParts), AllAnnihilated, AllNoDied, AllNoBorn, &
                ProjectionE, AvDiagSft, AllENumCyc / AllHFCyc, AllNoatHF, &
                AllNoatDoubs, AccRat, AllTotWalkers, IterTime
#endif

            if (tTruncInitiator .or. tDelayTruncInit) then
               write(initiatorstats_unit,"(I12,2I15,2G16.7,3F18.7)")&
                   Iter + PreviousCycles, AllNoAborted, AllNoAddedInitiators,&
                   real(AllNoInitDets) / real(AllNoNonInitDets), &
                   real(AllNoInitWalk) / real(AllNoNonInitWalk), &
                   DiagSftAbort, &
                   AllNoInitDets / real(StepsSft), &
                   AllNoNonInitDets / real(StepsSft)
            endif

            if (tLogComplexPops) then
                write (complexstats_unit,"(I12,3G16.7,3I12)") &
                    Iter + PreviousCycles, DiagSft, DiagSftRe, DiagSftIm, &
                    sum(AllTotParts), AllTotParts(1), AllTotParts(lenof_sign)
            endif

            call flush(6)
            call flush(fcimcstats_unit)
            
        endif

    end subroutine WriteFCIMCStats


    SUBROUTINE SetupParameters()
        use SystemData, only : tUseBrillouin,iRanLuxLev,tSpn,tHPHFInts,tRotateOrbs,tNoBrillouin,tROHF,tFindCINatOrbs,nOccBeta,nOccAlpha,tUHF
        use SystemData, only : tFixLz,LzTot,BasisFN,tBrillouinsDefault
        USE dSFMT_interface , only : dSFMT_init
        use CalcData, only : tFCIMC
        use CalcData , only : tRandomiseHashOrbs
        use CalcData, only : VirtCASorbs,OccCASorbs,G_VMC_Seed
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,TauFactor,StepsSftImag,tCheckHighestPop
        use Determinants , only : GetH0Element3
        use SymData , only : nSymLabels,SymLabelList,SymLabelCounts,TwoCycleSymGens
        use Logging , only : tTruncRODump
        use DetCalcData, only : NMRKS,tagNMRKS,FCIDets
        use SymExcit3, only : CountExcitations3 
        use DetBitOps, only: CountBits
        use constants, only: bits_n_int
        use util_mod, only: get_free_unit
        use sym_mod
        use HElem
        INTEGER :: ierr,i,j,HFDetTest(NEl),Seed,alpha,beta,symalpha,symbeta,endsymstate
        INTEGER :: HFConn,LargestOrb,nBits,HighEDet(NEl)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfDBO)
        HElement_t :: TempHii
        TYPE(BasisFn) HFSym
        REAL*8 :: TotDets,SymFactor,r
        CHARACTER(len=*), PARAMETER :: this_routine='SetupParameters'
        CHARACTER(len=12) :: abstr
        LOGICAL :: tSuccess,tFoundOrbs(nBasis),FoundPair
        REAL :: Gap
        INTEGER :: nSingles,nDoubles,HFLz,ChosenOrb,KPnt(3)

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
                OPEN(fcimcstats_unit,file='FCIMCStats',status='unknown')
            end if
            IF(tTruncInitiator.or.tDelayTruncInit) THEN
                initiatorstats_unit = get_free_unit()
                OPEN(initiatorstats_unit,file='INITIATORStats',status='unknown',form='formatted')
            ENDIF
            IF(tLogComplexPops) THEN
                ComplexStats_unit = get_free_unit()
                OPEN(ComplexStats_unit,file='COMPLEXStats',status='unknown')
            ENDIF
        ENDIF

!Store information specifically for the HF determinant
        ALLOCATE(HFDet(NEl),stat=ierr)
        CALL LogMemAlloc('HFDet',NEl,4,this_routine,HFDetTag)
        do i=1,NEl
            HFDet(i)=FDet(i)
        enddo
        HFHash=CreateHash(HFDet)
        CALL GetSym(HFDet,NEl,G1,NBasisMax,HFSym)
        WRITE(6,"(A,I10)") "Symmetry of reference determinant is: ",INT(HFSym%Sym%S,4)
        IF(tKPntSym) THEN
            CALL DecomposeAbelianSym(HFSym%Sym%S,KPnt)
            WRITE(6,"(A,3I5)") "Crystal momentum of reference determinant is: ",KPnt(1),KPnt(2),KPnt(3)
        ENDIF
        IF(tFixLz) THEN
            CALL GetLz(HFDet,NEl,HFLz)
            WRITE(6,"(A,I5)") "Ml value of reference determinant is: ",HFLz
            IF(HFLz.ne.LzTot) THEN
                CALL Stop_All("SetupParameters","Chosen reference determinant does not have the same Lz value as indicated in the input.")
            ENDIF
        ENDIF

!Do a whole lot of tests to see if we can use Brillouins theorem or not.
        IF(tBrillouinsDefault) CALL CheckforBrillouins() 
        
!test the encoding of the HFdet to bit representation.
        ALLOCATE(iLutHF(0:NIfTot),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for iLutHF")

        ! Test that the bit operations are working correctly...
        ! TODO: Move this to using the extract_bit_det routines to test those
        !       too...
        CALL EncodeBitDet(HFDet,iLutHF)
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

        !iLutRef is the reference determinant for the projected energy.
        ALLOCATE(iLutRef(0:NIfTot),stat=ierr)
        ALLOCATE(ProjEDet(NEl),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for iLutRef")
        
        ! The reference / projected energy determinants are the same as the
        ! HF determinant.
        ! TODO: Make these pointers rather than copies?
        iLutRef = iLutHF
        ProjEDet = HFDet

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


!If using a CAS space truncation, write out this CAS space
        IF(tTruncCAS) THEN
            IF(tTruncInitiator) THEN
                WRITE(6,'(A)') " *********** INITIATOR METHOD IN USE ***********"
                WRITE(6,'(A)') " Fixed initiator space defined using the CAS method."
            ELSE
                WRITE(6,*) "Truncated CAS space detected. Writing out CAS space..."
            ENDIF
            WRITE(6,'(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (",(OccCASOrbs+VirtCASOrbs)/2,",",OccCASOrbs,")"
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
 
!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        IF(tUEG.or.tHub) THEN
            exflag=2
        ELSE
            exflag=3
        ENDIF
        IF(.not.tKPntSym) THEN
!Count all possible excitations - put into HFConn
!TODO: Get CountExcitations3 working with tKPntSym
            CALL CountExcitations3(HFDet,exflag,nSingles,nDoubles)
            HFConn=nSingles+nDoubles
        ENDIF

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=abs(G_VMC_Seed-iProcIndex)
        WRITE(6,*) "Value for seed is: ",Seed
!Initialise...
        CALL dSFMT_init(Seed)
        
        IF(tRandomiseHashOrbs) THEN
            ALLOCATE(RandomHash(nBasis),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All(this_routine,"Error in allocating RandomHash")
            ENDIF
            RandomHash(:)=0
            IF(iProcIndex.eq.root) THEN
                do i=1,nBasis
                    FoundPair=.false.
                    do while(.not.FoundPair)
                        r = genrand_real2_dSFMT()
                        ChosenOrb=INT(nBasis*r*1000)+1
                        do j=1,nBasis
                            IF(RandomHash(j).eq.ChosenOrb) EXIT
                        enddo
                        IF(j.eq.nBasis+1) THEN
                            RandomHash(i)=ChosenOrb
                            FoundPair=.true.
                        ELSE
                            FoundPair=.false.
                        ENDIF
                    enddo
                enddo

!                WRITE(6,*) "Random Orbital Indexing for hash:"
!                WRITE(6,*) RandomHash(:)
                do i=1,nBasis
                    IF((RandomHash(i).eq.0).or.(RandomHash(i).gt.nBasis*1000)) THEN
                        CALL Stop_All(this_routine,"Random Hash incorrectly calculated")
                    ENDIF
                    do j=i+1,nBasis
                        IF(RandomHash(i).eq.RandomHash(j)) THEN
                            CALL Stop_All(this_routine,"Random Hash incorrectly calculated")
                        ENDIF
                    enddo
                enddo
            ENDIF
            !Now broadcast to all processors
            CALL MPIBCast(RandomHash,nBasis,Root)
        ENDIF

        IF(tHPHF) THEN
            !IF(tLatticeGens) CALL Stop_All("SetupParameters","Cannot use HPHF with model systems currently.")
            IF(tROHF.or.(LMS.ne.0)) CALL Stop_All("SetupParameters","Cannot use HPHF with high-spin systems.")
            tHPHFInts=.true.
        ENDIF

!Calculate Hii
        IF(tHPHF) THEN
            TempHii = hphf_diag_helement (HFDet, iLutHF)
        ELSE
            TempHii = get_helement (HFDet, HFDet, 0)
        ENDIF
        Hii=REAL(TempHii,dp)
        WRITE(6,*) "Reference Energy set to: ",Hii
        TempHii=GetH0Element3(HFDet)
        Fii=REAL(TempHii,dp)

!Find the highest energy determinant...
        IF(.not.tSpn) THEN
            do i=1,NEl
                HighEDet(i)=Brr(nBasis-(i-1))
            enddo
            IF(tHPHF) THEN
                call EncodeBitDet (HighEDet, iLutTemp)
                TempHii = hphf_diag_helement (HighEDet, iLutTemp)
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
                WRITE(6,*) "High spin calculation - but single excitations will *NOT* be used to calculate energy as this is an unrestricted calculation."
            ELSE
                CALL Stop_All("SetupParameters","High-spin, restricted calculation detected, but single excitations are not being used to calculate the energy.  &
                              & Either use the UHF keyword, or turn off brillouins theorem using NOBRILLOUINS, ROHF or ROTATEDORBS.")
            ENDIF
!            tRotatedOrbs=.true.
!        ELSEIF(LMS.ne.0) THEN
!            CALL Stop_All(this_routine,"Ms not equal to zero, but tSpn is false. Error here")
        ENDIF

!Initialise variables for calculation on each node
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
        HFPopCyc=0
        ENumCyc=0.D0
        ProjEIter=0.D0
        ProjEIterSum=0.D0
        VaryShiftCycles=0
        AvDiagSft=0.D0
        SumDiagSft=0.D0
        SumDiagSftAbort=0.D0
        AvDiagSftAbort=0.D0
        NoAborted=0.D0
        NoAddedInitiators=0.D0
        NoInitDets=0.D0
        NoNonInitDets=0.D0
        NoInitWalk=0.D0
        NoNonInitWalk=0.D0
        NoExtraInitDoubs=0.D0
        InitRemoved=0.D0
        TotImagTime=0.D0
        DiagSftRe=0.D0
        DiagSftIm=0.D0

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
        CullInfo(1:10,1:3)=0
        NoCulls=0
        AllNoAborted=0.D0
        AllNoAddedInitiators=0.D0
        AllNoInitDets=0.D0
        AllNoNonInitDets=0.D0
        AllNoInitWalk=0.D0
        AllNoNonInitWalk=0.D0
        AllNoExtraInitDoubs=0.D0
        AllInitRemoved=0.D0

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
                CALL Stop_All(this_routine,"A Full Diagonalization is required in the same calculation before histogramming can occur.")
            ENDIF

            IF(tHistHamil) THEN
                WRITE(6,*) "Histogramming total Hamiltonian, with Dets=", Det
                ALLOCATE(HistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                HistHamil(:,:)=0.D0
                ALLOCATE(AvHistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                AvHistHamil(:,:)=0.D0
                IF(iProcIndex.eq.0) THEN
                    ALLOCATE(AllHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                    AllHistHamil(:,:)=0.D0
                    ALLOCATE(AllAvHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                    AllAvHistHamil(:,:)=0.D0
                ENDIF
            ELSE
                WRITE(6,*) "Histogramming spawning wavevector, with Dets=", Det
                ALLOCATE(Histogram(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                Histogram(:,:)=0.D0
                ALLOCATE(AllHistogram(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
            ENDIF
            IF(tHistSpawn) THEN
                ALLOCATE(InstHist(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                InstHist(:,:)=0.D0
                ALLOCATE(AvAnnihil(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                AvAnnihil(:,:)=0.D0
                ALLOCATE(InstAnnihil(1:lenof_sign,1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
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
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
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

!Need to declare a new MPI type to deal with the long integers we use in the hashing, and when reading in from POPSFILEs
!        CALL MPI_Type_create_f90_integer(18,mpilongintegertype,error)
!        CALL MPI_Type_commit(mpilongintegertype,error)
        IF(tUseBrillouin) THEN
            WRITE(6,*) "Brillouin theorem specified, but this will not be in use with the non-uniform excitation generators."
        ENDIF
        WRITE(6,*) "Non-uniform excitation generators in use."
        IF(tKPntSym) THEN
            !Since we didn't calculate HFConn before, we calculate it in this routine and return it.
            CALL CalcApproxpDoubles(HFConn)
        ELSE
            CALL CalcApproxpDoubles()
        ENDIF
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

        IF(TPopsFile) THEN
            IF(mod(iWritePopsEvery,StepsSft).ne.0) CALL Warning(this_routine,"POPSFILE writeout should be a multiple of the update cycle length.")
        ENDIF

        WRITE(6,*) "*Direct Annihilation* in use...Explicit load-balancing disabled."
        WRITE(6,*) ""
        ALLOCATE(ValidSpawnedList(0:nProcessors-1),stat=ierr)
        ALLOCATE(InitialSpawnedSlots(0:nProcessors-1),stat=ierr)
!InitialSpawnedSlots now holds the first free position in the newly-spawned list for each processor, so it does not need to be reevaluated each iteration.
        MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
        Gap=REAL(MaxSpawned)/REAL(nProcessors)
        do j=0,nProcessors-1
            InitialSpawnedSlots(j)=NINT(Gap*j)+1
        enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

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
        IF(tTruncCAS) THEN
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

            IF(OccCASOrbs.gt.NEl) CALL Stop_All("SetupParameters","Occupied orbitals in CAS space specified is greater than number of electrons")
            IF(VirtCASOrbs.gt.(nBasis-NEl)) CALL Stop_All("SetupParameters","Virtual orbitals in CAS space specified is greater than number of unoccupied orbitals")

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
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the lowest ',NPartFrozen,' spin orbitals so that no more than ',NHolesFrozen,' holes exist within this core.'
            CALL CreateSpinInvBRR()
        ENDIF
        IF(tPartFreezeVirt) THEN
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the highest ',NVirtPartFrozen,' virtual spin orbitals so that no more than ',NElVirtFrozen,' electrons occupy these orbitals.'
            CALL CreateSpinInvBRR()
        ENDIF

        SymFactor=(Choose(NEl,2)*Choose(nBasis-NEl,2))/(HFConn+0.D0)
        TotDets=1.D0
        do i=1,NEl
            WRITE(6,*) "Approximate excitation level population: ",i,NINT((Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor)
            TotDets=TotDets+(Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor
        enddo
        WRITE(6,*) "Approximate size of determinant space is: ",NINT(TotDets)

        IF(TStartSinglePart) THEN
            WRITE(6,"(A,F9.3,A,I9)") " Initial number of particles set to 1, and shift will be held at ",DiagSft," until particle number on root node gets to ",InitWalkers
        ELSE
            WRITE(6,*) "Initial number of walkers per processor chosen to be: ", InitWalkers
        ENDIF

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
                WRITE(6,'(A)') " Using an open shell reference determinant in a basis of restricted HF orbitals; Brillouins theorem is being turned off. "
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

!This function will tell us whether we should allow attempted spawning at an excitation when we are truncating the space.
!We pass in the excitation level of the original particle, the two representations of the excitation (we only need the bit-representation of the excitation
!for HPHF) and the magnitude of the excitation (for determinant representation).
    LOGICAL FUNCTION CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)
        INTEGER :: nJ(NEl),WalkExcitLevel,ExcitLevel,IC,iGetExcitLevel_2,i,NoInFrozenCore,MinVirt
        INTEGER(KIND=n_int) :: iLutnJ(0:NIfTot)
        INTEGER :: kx,ky,kz ! For UEG

        CheckAllowedTruncSpawn=.true.

        IF(tTruncSpace) THEN
!We are truncating the space by excitation level
            IF(tHPHF) THEN
!With HPHF, we can't rely on this, since one excitation could be a single, and one a double. Also, IC is not returned.
                ExcitLevel = FindBitExcitLevel(iLutHF, iLutnJ, ICILevel)
                IF(ExcitLevel.gt.ICILevel) THEN
                    CheckAllowedTruncSpawn=.false.
                ENDIF

            ELSE
!Determinant representation.

                IF(WalkExcitLevel.eq.(ICILevel-1)) THEN
!The current walker is one below the excitation cutoff - if IC is a double, then could go over - we need to check
                    
                    IF(IC.eq.2) THEN
                        ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                    ELSE
!Always allow this - a single cannot put us over the truncated excitation level
                        ExcitLevel=0
!                        CheckAllowedTruncSpawn=.true.
!                        RETURN
                    ENDIF
                    IF(ExcitLevel.gt.ICILevel) THEN
                        CheckAllowedTruncSpawn=.false.
!                    ELSE
!                        CheckAllowedTruncSpawn=.true.
                    ENDIF

                ELSEIF(WalkExcitLevel.ge.ICILevel) THEN
!Walker is at the excitation cutoff level - all possible excitations could be disallowed - check the actual excitation level
                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                    IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                        CheckAllowedTruncSpawn=.false.
!                    ELSE
!                        CheckAllowedTruncSpawn=.true.
                    ENDIF
!                ELSE
!!Excitation cannot be in a dissallowed excitation level - allow it as normal
!                    CheckAllowedTruncSpawn=.true.
                ENDIF 

            ENDIF   !endif tHPHF

        ENDIF   !endif tTruncSpace

        IF((tTruncCAS.and.(.not.tTruncInitiator)).and.CheckAllowedTruncSpawn) THEN
!This flag determines if the FCI space is restricted by whether the determinants are in the predescribed CAS.
!            IF(.not.TestIfDetinCAS(nJ)) THEN
            IF(.not.TestIfDetinCASBit(iLutnJ(0:NIfD))) THEN
!Excitation not in allowed CAS space.
                CheckAllowedTruncSpawn=.false.
            ENDIF

        ENDIF

        IF(tPartFreezeCore) THEN
!Want to check if the determinant we're about to spawn on has more than the restricted number of holes in the partially frozen core.            

!Run through the electrons in nJ, count the number in the partially frozen core - ie those occupying orbitals with energy (from BRR) 
!less than that of the partially frozen core limit.
!If this is less than NPartFrozen-NHolesFrozen then spawning is forbidden.
            NoInFrozenCore=0
!BRR(i)=j: orbital i is the j-th lowest in energy  
            do i=1,NEl
                IF(SpinInvBRR(nJ(i)).le.NPartFrozen) NoInFrozenCore=NoInFrozenCore+1
                IF(NoInFrozenCore.eq.(NPartFrozen-NHolesFrozen)) EXIT   ! Can exit out of the loop if this is satisfied, since excitation will definitely be accepted.
            enddo
            IF(NoInFrozenCore.lt.(NPartFrozen-NHolesFrozen)) THEN
!There are more holes in the partially frozen core than has been specified as allowed.
                CheckAllowedTruncSpawn=.false.
!            ELSE
!Either the 'partially frozen core' is completely full, or it has the allowed number of holes or less.                
!Allowed to spawn, CheckAllowedTruncSpawn=.true.
!                CheckAllowedTruncSpawn=.true.
            ENDIF

        ENDIF

        IF(tPartFreezeVirt) THEN
!Want to check if the determinant we're about to spawn on has more than the restricted number of electrons in the partially frozen virtual orbitals.

!Run through the electrons in nJ, count the number in the partially frozen virtual orbitals - ie those occupying orbitals with energy (from BRR) 
!greater than that of the minimum unfrozen virtual.
!If this is greater than NElVirtFrozen then spawning is forbidden.
            NoInFrozenCore=0
            MinVirt=nBasis-NVirtPartFrozen
!BRR(i)=j: orbital i is the j-th lowest in energy  
            do i=1,NEl
                IF(SpinInvBRR(nJ(i)).gt.MinVirt) NoInFrozenCore=NoInFrozenCore+1
                IF(NoInFrozenCore.gt.NElVirtFrozen) THEN
!There are more electrons in the partially frozen virtual orbitals than has been specified as allowed.
                    CheckAllowedTruncSpawn=.false.
                    EXIT   ! Can exit out of the loop if this is satisfied, since excitation will definitely be accepted.
                ENDIF
            enddo
!Either the 'partially frozen virtual orbitals' are completely empty, or have the allowed number of electrons or less.                
!Allowed to spawn, CheckAllowedTruncSpawn=.true.
!                CheckAllowedTruncSpawn=.true.
!Want to be able to make it unable to spawn, but not able to spawn again.
        ENDIF

        IF(tUEG.and.(.not.tLatticeGens)) THEN
!Check to see if this is an allowed excitation
!by summing kx, ky and kz to zero over all the electrons.
            kx=0
            ky=0
            kz=0
            do i=1,NEl
                kx=kx+G1(nJ(i))%k(1)
                ky=ky+G1(nJ(i))%k(2)
                kz=kz+G1(nJ(i))%k(3)
            enddo
            if( .not.((kx.eq.0) .and. (ky.eq.0) .and. (kz.eq.0)) ) then
                CheckAllowedTruncSpawn=.false.
            endif
        ENDIF


    END FUNCTION CheckAllowedTruncSpawn



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
    
    SUBROUTINE CalcApproxpDoubles(HFConn)
        use SystemData , only : tAssumeSizeExcitgen,tUseBrillouin
        use CalcData , only : SinglesBias
        use SymData , only : SymClassSize
        use SymExcit3 , only : CountExcitations3
        INTEGER :: iTotal
        integer :: nSing, nDoub, ncsf, excitcount, ierr, iExcit
        integer :: nStore(6), iMaxExcit, nExcitMemLen, nJ(nel)
        integer, allocatable :: EXCITGEN(:)
        integer, intent(out), optional :: HFConn
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
                WRITE(6,*) "Since we are using a real-space hubbard model, only single excitations are connected."
                WRITE(6,*) "Setting pDoub to 0.D0"
                pDoubles=0.D0
                RETURN
            ELSE
                WRITE(6,*) "Since we are using a momentum-space hubbard model/UEG, only double excitaitons are connected."
                WRITE(6,*) "Setting pDoub to 1.D0"
                pDoubles=1.D0
                RETURN
            ENDIF
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
        HFConn=iTotal

        WRITE(6,"(I7,A,I7,A)") NDoub, " double excitations, and ",NSing," single excitations found from HF. This will be used to calculate pDoubles."

        IF(SinglesBias.ne.1.D0) THEN
            WRITE(6,*) "Singles Bias detected. Multiplying single excitation connectivity of HF determinant by ",SinglesBias," to determine pDoubles."
        ENDIF

        IF((NSing+nDoub+ncsf).ne.iTotal) THEN
            CALL Stop_All("CalcApproxpDoubles","Sum of number of singles and number of doubles does not equal total number of excitations")
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

            WRITE(6,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ",real(nDoub,dp)/real(iTotal,dp)
        ELSE
            write (6,'(A,F14.6)') " pDoubles set to: ", pDoubles
            write (6,'(A,F14.6)') " pSingles set to: ", pSingles
        ENDIF

        WRITE(6,'(A,F15.10)') " Assuming an average K_ij magnitude of approx 0.01, an appropriate tau is predicted to be around: ",(0.02*(1.D0/(REAL(NSing)+REAL(NDoub))))/0.01
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
            CALL GenExcitations3(HFDet,iLutHF,nJ,exflag,ExcitMat3,tParity,tAllExcitFound)
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

!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
!AJWT added optional argument dProbFin which is a probability that whatever gave this contribution was generated.
!  It defaults to 1, and weights the contribution of this det (only in the projected energy) by dividing its contribution by this number 
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,WSign,iLutCurr,HDiagCurr,dProbFin)
        use SystemData, only : tNoBrillouin
        use CalcData, only: tFCIMC
        INTEGER , intent(in) :: DetCurr(NEl),ExcitLevel
        INTEGER, DIMENSION(lenof_sign) , INTENT(IN) :: WSign
        INTEGER(KIND=n_int), intent(in) :: iLutCurr(0:NIfTot)
        INTEGER :: i,Bin
        INTEGER :: PartInd,OpenOrbs
        INTEGER(KIND=n_int) :: iLutSym(0:NIfTot)
        LOGICAL :: tSuccess
        REAL*8 , intent(in) :: HDiagCurr,dProbFin
        HElement_t :: HOffDiag

        IF(ExcitLevel.eq.0) THEN
            IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+WSign
            NoatHF=NoatHF+WSign
            HFCyc=HFCyc+WSign      !This is simply the number at HF*sign over the course of the update cycle 
            
        ELSEIF(ExcitLevel.eq.2) THEN
            NoatDoubs=NoatDoubs+sum(abs(WSign(:)))
!At double excit - find and sum in energy
            IF(tHPHF) THEN
                HOffDiag = hphf_off_diag_helement (ProjEDet, DetCurr, iLutRef, &
                                                   iLutCurr)
            ELSE
                HOffDiag = get_helement (ProjEDet, DetCurr, ExcitLevel, iLutRef, &
                                         iLutCurr)
            ENDIF
            IF(lenof_sign.eq.1) THEN
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+(REAL(HOffDiag,dp)*WSign(1)/dProbFin)
                ENumCyc=ENumCyc+(REAL(HOffDiag,dp)*WSign(1)/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
            ELSE
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+(HOffDiag*CMPLX(WSign(1),WSign(2),dp))/dProbFin
                ENumCyc=ENumCyc+(HOffDiag*CMPLX(WSign(1),WSign(2),dp))/dProbFin     !This is simply the Hij*sign summed over the course of the update cycle
            ENDIF

        ELSEIF(ExcitLevel.eq.1) THEN
          if(tNoBrillouin) then
!For the real-space hubbard model, determinants are only connected to excitations one level away, and brillouins theorem can not hold.
!For Rotated orbitals, brillouins theorem also cannot hold, and energy contributions from walkers on singly excited determinants must
!be included in the energy values along with the doubles.
            IF(tHPHF) THEN
                HOffDiag = hphf_off_diag_helement (ProjEDet, DetCurr, iLutRef, &
                                                   iLutCurr)
            ELSE
                HOffDiag = get_helement (ProjEDet, DetCurr, ExcitLevel, ilutRef, &
                                         iLutCurr)
            ENDIF
            IF(lenof_sign.eq.1) THEN
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+(REAL(HOffDiag,dp)*WSign(1)/dProbFin)
                ENumCyc=ENumCyc+(REAL(HOffDiag,dp)*WSign(1)/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
            ELSE
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+(HOffDiag*CMPLX(WSign(1),WSign(2),dp))/dProbFin
                ENumCyc=ENumCyc+(HOffDiag*CMPLX(WSign(1),WSign(2),dp))/dProbFin     !This is simply the Hij*sign summed over the course of the update cycle
            ENDIF
          endif 

        ENDIF

!Histogramming diagnostic options...
        IF((tHistSpawn.or.(tCalcFCIMCPsi.and.tFCIMC)).and.(Iter.ge.NHistEquilSteps)) THEN
            IF(ExcitLevel.eq.NEl) THEN
                CALL BinSearchParts2(iLutCurr,HistMinInd(ExcitLevel),Det,PartInd,tSuccess)
                if(tFCIMC) HistMinInd(ExcitLevel)=PartInd  !CCMC doesn't sum particle contributions in order, so we must search the whole space again
            ELSEIF(ExcitLevel.eq.0) THEN
                PartInd=1
                tSuccess=.true.
            ELSE
                CALL BinSearchParts2(iLutCurr,HistMinInd(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartInd,tSuccess)
                if(tFCIMC) HistMinInd(ExcitLevel)=PartInd  !CCMC doesn't sum particle contributions in order, so we must search the whole space again
            ENDIF
            IF(tSuccess) THEN
                IF(tHPHF) THEN
                    CALL FindExcitBitDetSym(iLutCurr,iLutSym)
                    IF(.not.DetBitEQ(iLutCurr,iLutSym)) THEN
                        IF(tFlippedSign) THEN
                            Histogram(1,PartInd)=Histogram(1,PartInd)-(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                            IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)-(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                            IF(tHistSpawn) THEN 
                                InstHist(1,PartInd)=InstHist(1,PartInd)-(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)-(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                            ENDIF
                        ELSE
                            Histogram(1,PartInd)=Histogram(1,PartInd)+(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                            IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)+(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                            IF(tHistSpawn) THEN
                                InstHist(1,PartInd)=InstHist(1,PartInd)+(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)+(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                            ENDIF
                        ENDIF
                    ELSE
                        IF(tFlippedSign) THEN
                            Histogram(1,PartInd)=Histogram(1,PartInd)-REAL(WSign(1),dp)/dProbFin
                            IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)-REAL(WSign(lenof_sign),dp)/dProbFin
                            IF(tHistSpawn) THEN
                                InstHist(1,PartInd)=InstHist(1,PartInd)-REAL(WSign(1),dp)/dProbFin
                                IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)-REAL(WSign(lenof_sign),dp)/dProbFin
                            ENDIF
                        ELSE
                            Histogram(1,PartInd)=Histogram(1,PartInd)+REAL(WSign(1),dp)/dProbFin
                            IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)+REAL(WSign(lenof_sign),dp)/dProbFin
                            IF(tHistSpawn) THEN
                                InstHist(1,PartInd)=InstHist(1,PartInd)+REAL(WSign(1),dp)/dProbFin
                                IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)+REAL(WSign(lenof_sign),dp)/dProbFin
                            ENDIF
                        ENDIF
                    ENDIF
                ELSE    !not HPHF
                    IF(tFlippedSign) THEN
                        Histogram(1,PartInd)=Histogram(1,PartInd)-REAL(WSign(1),dp)/dProbFin
                        IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)-REAL(WSign(lenof_sign),dp)/dProbFin
                        IF(tHistSpawn) THEN
                            InstHist(1,PartInd)=InstHist(1,PartInd)-REAL(WSign(1),dp)/dProbFin
                            IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)-REAL(WSign(lenof_sign),dp)/dProbFin
                        ENDIF
                    ELSE
                        Histogram(1,PartInd)=Histogram(1,PartInd)+REAL(WSign(1),dp)/dProbFin
                        IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)+REAL(WSign(lenof_sign),dp)/dProbFin
                        IF(tHistSpawn) THEN
                            InstHist(1,PartInd)=InstHist(1,PartInd)+REAL(WSign(1),dp)/dProbFin
                            IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)+REAL(WSign(lenof_sign),dp)/dProbFin
                        ENDIF
                    ENDIF
                ENDIF
                IF(tHPHF) THEN
!With HPHF space, we need to also include the spin-coupled determinant, which will have the same amplitude as the original determinant, unless it is antisymmetric.
                    IF(.not.DetBitEQ(iLutCurr,iLutSym)) THEN
                        IF(ExcitLevel.eq.NEl) THEN
                            CALL BinSearchParts2(iLutSym,FCIDetIndex(ExcitLevel),Det,PartInd,tSuccess)
                        ELSEIF(ExcitLevel.eq.0) THEN
                            PartInd=1
                            tSuccess=.true.
                        ELSE
                            CALL BinSearchParts2(iLutSym,FCIDetIndex(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartInd,tSuccess)
                        ENDIF
                        IF(tSuccess) THEN
                            CALL CalcOpenOrbs(iLutSym,OpenOrbs)
                            IF(tFlippedSign) THEN
                                IF(mod(OpenOrbs,2).eq.1) THEN
                                    Histogram(1,PartInd)=Histogram(1,PartInd)+(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                    IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)+(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) THEN
                                        InstHist(1,PartInd)=InstHist(1,PartInd)+(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                        IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)+(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    ENDIF
                                ELSE
                                    Histogram(1,PartInd)=Histogram(1,PartInd)-(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                    IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)-(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) THEN
                                        InstHist(1,PartInd)=InstHist(1,PartInd)-(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                        IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)-(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    ENDIF
                                ENDIF
                            ELSE
                                IF(mod(OpenOrbs,2).eq.1) THEN
                                    Histogram(1,PartInd)=Histogram(1,PartInd)-(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                    IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)-(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) THEN
                                        InstHist(1,PartInd)=InstHist(1,PartInd)-(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                        IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)-(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    ENDIF
                                ELSE
                                    Histogram(1,PartInd)=Histogram(1,PartInd)+(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                    IF(lenof_sign.eq.2) Histogram(lenof_sign,PartInd)=Histogram(lenof_sign,PartInd)+(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) THEN
                                        InstHist(1,PartInd)=InstHist(1,PartInd)+(REAL(WSign(1),dp)/SQRT(2.0))/dProbFin
                                        IF(lenof_sign.eq.2) InstHist(lenof_sign,PartInd)=InstHist(lenof_sign,PartInd)+(REAL(WSign(lenof_sign),dp)/SQRT(2.0))/dProbFin
                                    ENDIF
                                ENDIF
                            ENDIF
                        ELSE
                            WRITE(6,*) DetCurr(:)
                            WRITE(6,*) "***",iLutSym(0:NIfTot)
                            WRITE(6,*) "***",ExcitLevel,Det
                            CALL Stop_All("SumEContrib","Cannot find corresponding spin-coupled FCI determinant when histogramming")
                        ENDIF
                    ENDIF
                ENDIF
            ELSE
                WRITE(6,*) DetCurr(:)
                WRITE(6,*) "***",iLutCurr(0:NIfTot)
                WRITE(6,*) "***",ExcitLevel,HistMinInd(ExcitLevel),Det
                Call WriteBitDet(6,iLutCurr(0:NIfTot),.true.)
                CALL Stop_All("SumEContrib","Cannot find corresponding FCI determinant when histogramming")
            ENDIF
        ELSEIF(tHistEnergies) THEN
!This wil histogramm the energies of the particles, rather than the determinants themselves.
            Bin=INT(HDiagCurr/BinRange)+1
            IF(Bin.gt.iNoBins) THEN
                CALL Stop_All("SumEContrib","Histogramming energies higher than the arrays can cope with. Increase iNoBins or BinRange")
            ENDIF
            HistogramEnergy(Bin)=HistogramEnergy(Bin)+real(abs(WSign(1)),dp)
        ENDIF

        IF(tPrintOrbOcc.and.(Iter.ge.StartPrintOrbOcc)) THEN
            IF((tPrintOrbOccInit.and.(test_flag(iLutCurr,flag_is_initiator(1)))).or.(.not.tPrintOrbOccInit)) then
                do i=1,NEl
                    OrbOccs(DetCurr(i))=OrbOccs(DetCurr(i))+(REAL(WSign(1))*REAL(WSign(1)))
                enddo
            ENDIF
        ENDIF

        RETURN

    END SUBROUTINE SumEContrib


!This routine will change the reference determinant to DetCurr. It will also re-zero all the energy estimators, since they now correspond to
!projection onto a different determinant.
    SUBROUTINE ChangeRefDet(DetCurr)
        use Determinants , only : GetH0Element3
!        use FciMCLoggingMod , only : RestartBlocking, RestartShiftBlocking
        INTEGER :: DetCurr(NEl),i

!        CALL Stop_All("ChangeRefDet","This option does not currently work. Bug ghb24 if its needed")
!Problem is that we need to rerun the simulation from scratch, and particles currently in the simulation will keep on
!changing the reference since their diagonal K element will remain negative.

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
!            IF(TAutoCorr) CLOSE(44)
        ENDIF
        IF(TDebug) CLOSE(11)
        CALL SetupParameters()
        CALL InitFCIMCCalcPar()
!        IF(iProcIndex.eq.0) THEN
!            CALL RestartBlocking(Iter)
!            CALL RestartShiftBlocking(Iter)
!        ENDIF



    END SUBROUTINE ChangeRefDet
    
!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        use FciMCLoggingMOD , only : InitHistInitPops
        use SystemData , only : tRotateOrbs
        use CalcData , only : InitialPart
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn
        use constants , only : size_n_int
        INTEGER :: ierr
        INTEGER :: error,MemoryAlloc
        INTEGER, DIMENSION(lenof_sign) :: InitialSign
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'
            
        if (tReadPops .and. .not. tPopsAlreadyRead) then
!Read in particles from multiple POPSFILES for each processor
            WRITE(6,*) "Reading in initial particle configuration from POPSFILES..."
            CALL ReadFromPopsFilePar()
        ELSE
!initialise the particle positions - start at HF with positive sign
!Set the maximum number of walkers allowed
            MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)
            WRITE(6,"(A,I14)") " Memory allocated for a maximum particle number per node of: ",MaxWalkersPart
            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
!            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned

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
                WRITE(6,"(A,F14.6,A)") " Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*8,dp)/1048576.D0," Mb/Processor"
            ENDIF
            

            WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
            ALLOCATE(SpawnVec(0:NIftot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVecTag,ierr)
            SpawnVec(:,:)=0
            ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
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
        
            iHFProc=DetermineDetProc(iLutHF)   !This wants to return a value between 0 -> nProcessors-1
            WRITE(6,*) "HF processor is: ",iHFProc

            TotParts(:)=0
            TotPartsOld(:)=0
            NoatHF=0
!            AllTotPartsOld(:)=1     !So that the first update gives a meaningful number

!Setup initial walker local variables
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
                    iter_data_fciqmc%tot_parts_old = 1
                    AllTotParts(1)=REAL(InitialPart,dp)
                    AllTotPartsOld(1)=REAL(InitialPart,dp)
                    AllNoAbortedOld=0.D0
                ENDIF
            ELSE
!In this, only one processor has initial particles.
                IF(iProcIndex.eq.Root) THEN
                    AllTotWalkers=1.D0
                    AllTotWalkersOld = 1
                    iter_data_fciqmc%tot_parts_old = AllTotWalkers
                    AllTotParts(1)=REAL(InitWalkers,dp)
                    AllTotPartsOld(1)=REAL(InitWalkers,dp)
                    AllNoAbortedOld=0.D0
                ENDIF
            ENDIF
        
            WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens + temp arrays) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Only one array of memory to store main particle list allocated..."
            WRITE(6,*) "Initial memory allocation sucessful..."
            CALL FLUSH(6)

        ENDIF   !End if initial walkers method
            
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

        IF(tHistInitPops) THEN
            CALL InitHistInitPops()
        ENDIF
        tPrintHighPop=.false.
        MaxInitPopPos=0
        MaxInitPopNeg=0

        IF(MaxNoatHF.eq.0) THEN
            MaxNoatHF=InitWalkers*nProcessors
            HFPopThresh=MaxNoatHF
        ENDIF

        IF((NMCyc.ne.0).and.(tRotateOrbs.and.(.not.tFindCINatOrbs))) CALL Stop_All(this_routine,"Currently not set up to rotate and then go straight into a spawning &
                                                                                    & calculation.  Ordering of orbitals is incorrect.  This may be fixed if needed.")

    end subroutine InitFCIMCCalcPar

    SUBROUTINE DeallocFCIMCMemPar()
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

END MODULE FciMCParMod

            

!This is the same as BinSearchParts1, but this time, it searches though the full list of determinants created by the full diagonalizer when the histogramming option is on.
!This is outside the module so it is accessible to AnnihilateMod
SUBROUTINE BinSearchParts2(iLut,MinInd,MaxInd,PartInd,tSuccess)
    use DetCalcData , only : FCIDets
    use DetBitOps, only: DetBitLT
    use constants, only: n_int
    use bit_reps, only: NIfTot
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
    
