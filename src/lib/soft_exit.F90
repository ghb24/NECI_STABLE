! JSS, based almost entirely on work by GHB.

! During the calculation, test for the existence of the file CHANGEVARS.  

! Various input options can be given in CHANGEVARS, including the ability
! to modify parameters or cause the code to perform a "soft exit" as soon
! as possible rather than finishing the calculation.

! It is left to the programmer to act on the variables set/modified by
! CHANGEVARS.

! Two procedures are provided:
! 
! * CHANGEVARS reads the CHANGEVARS input file (if it exists) and sets the
!   relevant variables accordingly.
! * test_SoftExit is a function which is a simple wrapper around CHANGEVARS
!   and is true if a soft exit has been requested.

! If you wish to exit immediately, rather than trapping the value returned by
! test_SOFTEXIT and then tracing back through the calling stack and hence avoid
! any further unwanted calculations, use test_SOFTEXIT with the termination
! routines provided in NECICore.F90 and in error_handling.F90::
!     if (test_SOFTEXIT) then
!         call NECICalcEnd(0)
!         call NECICodeEnd(tCPMD,tVASP)
!         call quiet_stop()
!     end if

module soft_exit

    use SystemData, only: nel, nBasis, NIfTot
    use util_mod, only: binary_search
    use FciMCData, only: iter, CASMin, CASMax, tTruncSpace, tSinglePartPhase,&
                         SumENum, SumNoatHF, HFPopCyc, ProjEIterSum, &
                         Histogram, AvAnnihil, VaryShiftCycles, SumDiagSft, &
                         VaryShiftIter, CurrentDets, CurrentSign, iLutHF, &
                         TotWalkers,tPrintHighPop
    use CalcData, only: Tau, DiagSft, SftDamp, StepsSft, SinglesBias, &
                        OccCASOrbs, VirtCASOrbs, NMCyc, tTruncCAS, &
                        NEquilSteps, tTruncInitiator, InitiatorWalkNo
    use DetCalcData, only: ICILevel
    use IntegralsData, only: tPartFreezeCore, NPartFrozen, NHolesFrozen, &
                             NVirtPartFrozen, NelVirtFrozen, tPartFreezeVirt
    use Parallel
    use Input
    use Logging, only: tHistSpawn, tCalcFCIMCPsi, tIterStartBlock, &
                       IterStartBlocking, tHFPopStartBlock, NHistEquilSteps
    use FCIMCLoggingMOD, only: PrintBlocking, RestartBlocking, &
                               PrintShiftBlocking, RestartShiftBlocking
    use AnnihilationMod, only: DetermineDetProc
    implicit none

contains

    subroutine ChangeVars(tSingBiasChange,tSoftExitFound,tWritePopsFound)

        ! If a file is created called CHANGEVARS in one of the working directories of
        ! the run, multiple values can be changed.

        ! Out:
        !     tSingBiasChange: true if the single bias is changed.
        !     tSoftExitFound: true if a SOFTEXIT is requested.
        !     tWritePopsFound: true if the output of a POPSFILE has been requested.
        
        ! Other variables that can be changed are modified directly at the
        ! module-level.
        ! Ways that the simulation can be affected are:

        !   EXCITE  XXX      Will change the excitation level of the simulation (< 0 or > NEl sets it to the full space)
        !   TRUNCATECAS  XXX XXX    Will change the CAS of the simulation (< 0 or > NEl sets it to the full space)
        !   SOFTEXIT         Will exit cleanly from the program
        !   WRITEPOPS        Will write a current popsfile
        !   VARYSHIFT        Will exit out of fixed shift phase
        !   NMCYC XXX        Will change the number of monte carlo cycles to perform
        !   TAU XXX          Will change the value of tau for the simulation
        !   DIAGSHIFT XXX    Will change the shift
        !   SHIFTDAMP XXX    Will change the damping parameter
        !   STEPSSHIFT XXX   Will change the length of the update cycle
        !   SINGLESBIAS XXX  Will change the singles bias for the non-uniform random excitation generator
        !   ZEROPROJE        Will rezero the averaged energy estimators
        !   ZEROHIST         Will rezero the averaged histogramming vectors
        !   PARTIALLYFREEZE XXX XXX Will change the number of holes/electrons in the core valence region
        !   PARTIALLYFREEZEVIRT XXX XXX Will change the number of electrons in the partially frozen virtual region
        !   PRINTERRORBLOCKING  Will print the blocking analysis
        !   STARTERRORBLOCKING  Will start the blocking analysis
        !   RESTARTERRORBLOCKING    Will restart the blocking analysis
        !   PRINTSHIFTBLOCKING      Will print the shift blocking analysis
        !   RESTARTSHIFTBLOCKING    Will restart the shift blocking analysis
        !   EQUILSTEPS XXX          Will change the number of steps to ignore in the averaging of the energy and the shift.
        !   STARTHIST        Will begin histogramming the determinant populations if the tCalcFCIMCPsi is on and the histogramming has been set up.
        !   HISTEQUILSTEPS XXX      Will change the iteration at which the histogramming begins to the value specified.
        !   TRUNCINITIATOR   Will expand the CAS calculation to a TRUNCINITIATOR calculation if DELAYTRUNCINITIATOR is present in the input.
        !   ADDTOINIT XXX    Will change the cutt-off population for which walkers are added to the initiator space.  Pop must be *above* specified value.
        !   SCALEHF XXX      Will scale the number of walkers at HF by the specified factor
        !   PRINTHIGHPOPDET  Will print the determinant with the highest population of different sign to the HF.

       integer :: error,i,ios,NewNMCyc, pos
       logical :: tSoftExitFound,tWritePopsFound,exists,AnyExist,deleted_file
       logical :: tEof,any_deleted_file,tChangeParams(27),tSingBiasChange
       real*8 :: hfScaleFactor
       Character(len=100) :: w

       tSoftExitFound=.false.
       tWritePopsFound=.false.
       tSingBiasChange=.false.
       ios=0
       inquire(file='CHANGEVARS',exist=exists)
#ifdef PARALLEL
       !This collective will check the exists logical on all nodes, and perform a logical or operation,
       !before broadcasting the result back to all nodes.
       CALL MPI_AllReduce(exists,AnyExist,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)
#endif
       if(AnyExist) then
           IF(iProcIndex.eq.0) THEN
               WRITE(6,*) "CHANGEVARS file detected on iteration ",Iter
           ENDIF
!Set the defaults
           tChangeParams = .false.
           deleted_file=.false.
           do i=0,nProcessors-1
               ! This causes each processor to attempt to delete
               ! CHANGEVARS in turn (as each cycle of the loop involves waiting
               ! for all processors to reach the AllReduce before the next cycle 
               ! can start, and hence avoid race conditions between processors 
               ! sharing the same disk.
               if (i==iProcIndex.and.exists) then
                   open(13,file='CHANGEVARS',Status='OLD',err=99,iostat=ios)
                   ir=13
                   Call input_options(echo_lines=.true.,skip_blank_lines=.true.)
                   Do
                       Call read_line(tEof)
                       IF(tEof) Exit
                       Call ReadU(w)
                       Select Case(w)
                       CASE("TAU")
                           tChangeParams(1)=.true.
                           Call Readf(Tau)
                       CASE("DIAGSHIFT")
                           tChangeParams(2)=.true.
                           CALL Readf(DiagSft)
                       CASE("SHIFTDAMP")
                           tChangeParams(3)=.true.
                           CALL Readf(SftDamp)
                       CASE("STEPSSHIFT")
                           tChangeParams(4)=.true.
                           CALL Readi(StepsSft)
                       CASE("EXCITE")
                           tChangeParams(5)=.true.
                           CALL Readi(ICILevel)
                       CASE("SOFTEXIT")
                           tChangeParams(6)=.true.
                       CASE("WRITEPOPS")
                           tChangeParams(7)=.true.
                       CASE("VARYSHIFT")
                           tChangeParams(8)=.true.
                       CASE("SINGLESBIAS")
                           tChangeParams(9)=.true.
                           CALL Readf(SinglesBias)
                       CASE("TRUNCATECAS")
                           tChangeParams(10)=.true.
                           CALL Geti(OccCASOrbs)
                           CALL Geti(VirtCASOrbs)
                       CASE("NMCYC")
                           tChangeParams(11)=.true.
                           CALL Geti(NewNMCyc)
                       CASE("ZEROPROJE")
                           tChangeParams(12)=.true.
                       CASE("ZEROHIST")
                           tChangeParams(13)=.true.
                       CASE("PARTIALLYFREEZE")
                           tChangeParams(14)=.true.
                           CALL Readi(NPartFrozen)
                           CALL Readi(NHolesFrozen)
                       CASE("PRINTERRORBLOCKING")
                           tChangeParams(15)=.true.
                       CASE("STARTERRORBLOCKING")
                           tChangeParams(16)=.true.
                       CASE("RESTARTERRORBLOCKING")
                           tChangeParams(17)=.true.
                       CASE("PRINTSHIFTBLOCKING")
                           tChangeParams(18)=.true.
                       CASE("RESTARTSHIFTBLOCKING")
                           tChangeParams(19)=.true.
                       CASE("EQUILSTEPS")
                           tChangeParams(20)=.true.
                           CALL Readi(NEquilSteps)
                       CASE("STARTHIST")
                           tChangeParams(21)=.true.
                       CASE("HISTEQUILSTEPS")
                           tChangeParams(22)=.true.
                           CALL Readi(NHistEquilSteps)
                       CASE("TRUNCINITIATOR")
                           tChangeParams(23)=.true.
                       CASE("PARTIALLYFREEZEVIRT")                           
                           tChangeParams(24)=.true.
                           CALL Readi(NVirtPartFrozen)
                           CALL Readi(NElVirtFrozen)
                       CASE("ADDTOINIT")
                           tChangeParams(25)=.true.
                           CALL Readi(InitiatorWalkNo)
                       case("SCALEHF")
                           tChangeParams(26) = .true.
                           call readf(hfScaleFactor)
                       case("PRINTHIGHPOPDET")
                           tChangeParams(27) = .true.
                       END SELECT
                   End Do
                   close(13,status='delete')
                   deleted_file=.true.
               end if
#ifdef PARALLEL
               call MPI_AllReduce(deleted_file,any_deleted_file,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)
#endif
               if (any_deleted_file) exit
           end do
#ifdef PARALLEL
           CALL MPI_BCast(tChangeParams,27,MPI_LOGICAL,i,MPI_COMM_WORLD,error)
#endif

           IF(tChangeParams(1)) THEN
!Change Tau
#ifdef PARALLEL
               CALL MPI_BCast(Tau,1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Tau changed to a value of : ",Tau
               ENDIF
           ENDIF
           IF(tChangeParams(2)) THEN
!Change DiagSft
#ifdef PARALLEL
               CALL MPI_BCast(DiagSft,1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "DIAGSHIFT changed to a value of : ",DiagSft
               ENDIF
           ENDIF
           IF(tChangeParams(3)) THEN
!Change SftDamp
#ifdef PARALLEL
               CALL MPI_BCast(SftDamp,1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "SHIFTDAMP changed to a value of : ",SftDamp
               ENDIF
           ENDIF
           IF(tChangeParams(4)) THEN
!Change StepsSft
#ifdef PARALLEL
               CALL MPI_BCast(StepsSft,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "STEPSSHIFT changed to a value of : ",StepsSft
               ENDIF
           ENDIF
           IF(tChangeParams(5)) THEN
!Change Excite Level
               IF(.not.tTruncSpace) THEN
                   IF(iProcIndex.eq.0) THEN
                       WRITE(6,*) "The space is not truncated, so EXCITE keyword in the CHANGEVARS file will not affect run."
                   ENDIF
               ELSE
                   IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
                       IF(iProcIndex.eq.0) WRITE(6,*) "Cannot increase truncation level, since histogramming wavefunction..."
                   ELSE
#ifdef PARALLEL
                       CALL MPI_BCast(ICILevel,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
                       IF((ICILevel.le.0).or.(ICILevel.ge.NEl)) THEN
                           tTruncSpace=.false.
                           IF(iProcIndex.eq.0) THEN
                               WRITE(6,*) "Expanding the space to the full space."
                           ENDIF
                       ELSE
                           IF(iProcIndex.eq.0) THEN
                               WRITE(6,*) "Increasing truncation level of space to a value of ",ICILevel
                           ENDIF
                       ENDIF
                   ENDIF
               ENDIF
           ENDIF
           IF(tChangeParams(6)) THEN
!SoftExit detected
!We don't need to broadcast this as it can't go the other way!
               tSoftExitFound=.true.
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "SOFTEXIT triggered. Exiting out of run."
               ENDIF
           ENDIF
           IF(tChangeParams(7)) THEN
!Write Popsfile
               tWritePopsFound=.true.
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Asked to write out a POPSFILE..."
               ENDIF
           ENDIF
           IF(tChangeParams(8)) THEN
!Vary Shift
!We don't need to broadcast this as it can't go the other way!
               IF(.not.tSinglePartPhase) THEN
                   IF(iProcIndex.eq.0) THEN
                       WRITE(6,*) "Request to vary shift denied, since simulation already in variable shift mode..."
                   ENDIF
               ELSE
                   tSinglePartPhase=.false.
                   VaryShiftIter=Iter
                   IF(iProcIndex.eq.0) THEN
                       WRITE(6,*) "Request to vary the shift detected on a node..."
                   ENDIF
               ENDIF
           ENDIF
           IF(tChangeParams(9)) THEN
!Vary SinglesBias
#ifdef PARALLEL
               CALL MPI_BCast(SinglesBias,1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "SINGLESBIAS changed to a value of : ",SinglesBias
               ENDIF
               tSingBiasChange=.true.
           ENDIF
           IF(tChangeParams(10)) THEN
!Change CAS space
               IF(.not.tTruncCAS) THEN
                   IF(iProcIndex.eq.0) THEN
                       WRITE(6,*) "The space is not truncated by CAS, so TRUNCATECAS keyword in the CHANGEVARS file will not affect run."
                   ENDIF
               ELSE
#ifdef PARALLEL
                   CALL MPI_BCast(OccCASOrbs,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
                   CALL MPI_BCast(VirtCASOrbs,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
                   IF(((OccCASOrbs.ge.NEl).and.(VirtCASOrbs.ge.(nBasis-NEl))).or.(OccCASOrbs.le.0).or.(VirtCASOrbs.le.0)) THEN
!CAS space is equal to or greater than the full space, or one of the arguments is less than zero.
                       tTruncCAS=.false.
                       IF(iProcIndex.eq.0) THEN
                           WRITE(6,*) "Expanding CAS to the full space"
                       ENDIF
                   ELSE
                       CASMax=NEl+VirtCASOrbs
                       CASMin=NEl-OccCASOrbs
                       IF(iProcIndex.eq.0) THEN
                           WRITE(6,"(A,I5,A,I6)") "Increasing CAS space accessible to ",OccCASOrbs," , ",VirtCASOrbs
                       ENDIF
                   ENDIF
               ENDIF
           ENDIF
           IF(tChangeParams(11)) THEN
!Change number of MC steps
#ifdef PARALLEL
               CALL MPI_BCast(NewNMCyc,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
               IF((iProcIndex.eq.0).and.(NewNMCyc.lt.Iter)) THEN
                   WRITE(6,*) "New value of NMCyc is LESS than the current iteration number."
                   WRITE(6,*) "Therefore, the number of iterations has been left at ",NMCyc
               ELSEIF(iProcIndex.eq.0) THEN
!Only update if new number of cycles is greater than old set.
                   NMCyc=NewNMCyc
                   WRITE(6,*) "Total number of MC Cycles has been changed to ",NMCyc
               ELSEIF(NewNMCyc.ge.Iter) THEN
                   NMCyc=NewNMCyc
               ENDIF
           ENDIF
           IF(tChangeParams(12)) THEN
               SumENum=0.D0
               SumNoatHF=0
               HFPopCyc=0
               ProjEIterSum=0.D0
               VaryShiftCycles=0
               SumDiagSft=0.D0
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Zeroing all average energy estimators..."
               ENDIF
           ENDIF
           IF(tChangeParams(13)) THEN
               Histogram(:)=0.D0
               IF(tHistSpawn) AvAnnihil(:)=0.D0
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Zeroing all average histograms..."
               ENDIF
           ENDIF
           IF(tChangeParams(14)) THEN
#ifdef PARALLEL
               CALL MPI_BCast(NPartFrozen,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
               CALL MPI_BCast(NHolesFrozen,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) WRITE(6,'(A,I4,A,I4,A)') 'Allowing ',NHolesFrozen,' holes in ',NPartFrozen,' partially frozen orbitals.'
               IF(NHolesFrozen.eq.NPartFrozen) THEN
                   ! Allowing as many holes as there are orbitals - equivalent to not freezing at all.
                   tPartFreezeCore=.false.
#ifdef PARALLEL
                   CALL MPI_BCast(tPartFreezeCore,1,MPI_LOGICAL,i,MPI_COMM_WORLD,error)
#endif
                   IF(iProcIndex.eq.0) THEN
                       WRITE(6,*) 'Unfreezing any partially frozen core.'
                   ENDIF
               ELSE
                   tPartFreezeCore=.true.
               ENDIF
           ENDIF
           IF(tChangeParams(15)) THEN
               WRITE(6,*) 'Printing blocking analysis at this point.'
               IF(iProcIndex.eq.0) CALL PrintBlocking(Iter)
           ENDIF
           IF(tChangeParams(16)) THEN
               IF((.not.tHFPopStartBlock).and.(.not.tIterStartBlock)) THEN
                   WRITE(6,*) 'Error blocking already started'
                   ! Don't want to re-call the init routine.
                   ! If both have these have been turned off, then this is true.
               ELSE
                   tIterStartBlock=.true.
                   IterStartBlocking=Iter
               ENDIF
           ENDIF
           IF(tChangeParams(17)) THEN
               WRITE(6,*) 'Restarting the error calculations.  All blocking arrays are re-set to zero.'
               IF(iProcIndex.eq.0) CALL RestartBlocking(Iter)
           ENDIF
           IF(tChangeParams(18)) THEN
               WRITE(6,*) 'Printing shift error blocking at this point.'
               IF(iProcIndex.eq.0) CALL PrintShiftBlocking(Iter)
           ENDIF
           IF(tChangeParams(19)) THEN
               WRITE(6,*) 'Restarting the shift error calculations.  All shift blocking arrays set to zero.'
               IF(iProcIndex.eq.0) CALL RestartShiftBlocking(Iter)
           ENDIF
           IF(tChangeParams(20)) THEN
#ifdef PARALLEL
               CALL MPI_BCast(NEquilSteps,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Changing the number of equilibration steps to: ",NEquilSteps
               ENDIF
           ENDIF
           IF(tChangeParams(21)) THEN
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Beginning to histogram at the next update"
                   NHistEquilSteps=Iter+StepsSft
                   IF(.not.tCalcFCIMCPsi) WRITE(6,*) "This has no effect, as the histograms have not been set up at the beginning of the calculation."
               ENDIF
           ENDIF
           IF(tChangeParams(22)) THEN
               IF(NHistEquilSteps.le.Iter) NHistEquilSteps=Iter+StepsSft
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Changing the starting iteration for histogramming to: ",NHistEquilSteps
                   IF(.not.tCalcFCIMCPsi) WRITE(6,*) "This has no effect, as the histograms have not been set up at the beginning of the calculation."
               ENDIF
           ENDIF
           IF(tChangeParams(23)) THEN
               tTruncInitiator=.true.
               Tau=Tau/10.D0
#ifdef PARALLEL
               CALL MPI_BCast(Tau,1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,*) "Beginning to allow spawning into inactive space for a truncated initiator calculation."
                   WRITE(6,*) "Reducing tau by an order of magnitude.  The new tau is: ",Tau
               ENDIF
           ENDIF
           IF(tChangeParams(24)) THEN
#ifdef PARALLEL
               CALL MPI_BCast(NVirtPartFrozen,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
               CALL MPI_BCast(NElVirtFrozen,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) WRITE(6,'(A,I4,A,I4,A)') 'Allowing ',NElVirtFrozen,' electrons in ',NVirtPartFrozen,' partially frozen virtual orbitals.'
               IF(NElVirtFrozen.eq.NEl) THEN
                   ! Allowing as many holes as there are orbitals - equivalent to not freezing at all.
                   tPartFreezeVirt=.false.
#ifdef PARALLEL
                   CALL MPI_BCast(tPartFreezeVirt,1,MPI_LOGICAL,i,MPI_COMM_WORLD,error)
#endif
                   IF(iProcIndex.eq.0) THEN
                       WRITE(6,*) 'Unfreezing any partially frozen virtuals.'
                   ENDIF
               ELSE
                   tPartFreezeVirt=.true.
               ENDIF
           ENDIF           
           IF(tChangeParams(25)) THEN
#ifdef PARALLEL
               CALL MPI_BCast(InitiatorWalkNo,1,MPI_INTEGER,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) THEN
                   WRITE(6,'(A,I5)') "Cutoff population for determinants to be added to the initiator space changed to ",InitiatorWalkNo
               ENDIF
           ENDIF
           if (tChangeParams(26)) then
#ifdef PARALLEL
               call MPI_BCast(HFScaleFactor, 1, MPI_DOUBLE_PRECISION, i, &
                              MPI_COMM_WORLD, error)
#endif
               if (iProcIndex == 0) then
                   write (6, '(a,f10.4)') "Number at Hartree-Fock scaled by &
                          &factor: ", hfScaleFactor
               endif

               SumNoatHF = SumNoatHF * hfScaleFactor
               if (iProcIndex == DetermineDetProc(iLutHF)) then
                   pos = binary_search (CurrentDets, iLutHF, NIfTot+1, &
                                        TotWalkers)
                   CurrentSign(pos) = CurrentSign(pos) * hfScaleFactor
               endif
           endif
           IF(tChangeParams(27)) THEN
               tPrintHighPop=.true.
#ifdef PARALLEL
               CALL MPI_BCast(tPrintHighPop,1,MPI_LOGICAL,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) WRITE(6,'(A)') 'Request to print the determinant with the highest population of different sign to the HF.'
           ENDIF   
 
       endif

99     IF (ios.gt.0) THEN
          WRITE (6,*) 'Problem reading CHANGEVARS file '
!          call stop_all('ChangeVars','CHANGEVARS read error.')
       END IF
    
    end subroutine ChangeVars

logical function test_SoftExit()

    logical :: tdummy1, tdummy2

    call ChangeVars(tdummy1, test_SoftExit, tdummy2)
    if (test_SoftExit) write (6,'(1X,a30)') 'Request for SOFTEXIT detected.'

end function test_SoftExit

end module soft_exit
