! JSS, based almost entirely on work by GHB.
! Somewhat tidied by SDS, due to expansion of options.

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

    use SystemData, only: nel, nBasis
    use bit_reps, only: NIfTot
    use util_mod, only: binary_search
    use FciMCData, only: iter, CASMin, CASMax, tTruncSpace, tSinglePartPhase,&
                         SumENum, SumNoatHF, HFPopCyc, ProjEIterSum, &
                         Histogram, AvAnnihil, VaryShiftCycles, SumDiagSft, &
                         VaryShiftIter, CurrentDets, iLutHF, &
                         TotWalkers,tPrintHighPop
    use CalcData, only: Tau, DiagSft, SftDamp, StepsSft, SinglesBias, &
                        OccCASOrbs, VirtCASOrbs, NMCyc, tTruncCAS, &
                        NEquilSteps, tTruncInitiator, InitiatorWalkNo, &
                        tCheckHighestPop, tRestartHighPop, tChangeProjEDet, &
                        tCheckHighestPopOnce, FracLargerDet
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
    use constants, only: lenof_sign
    use bit_reps, only: extract_sign,encode_sign
    use spin_project, only: tSpinProject, spin_proj_gamma, &
                            spin_proj_interval, spin_proj_shift
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
        !   CHANGEREFDET     Will change the reference determinant to the
        !                    det with the highest population
        !   RESTARTHIGHPOP   Restart the calculation with same parameters but
        !                    a new reference determinant
        !   SPIN-PROJECT     Change the interval between applications of
        !                    stochastic spin projection. If 0, disable it.
        !   SPIN-PROJECT-GAMMA
        !                    Change the delta-gamma value used for stochastic
        !                    spin projection
        !   SPIN-PROJECT-SHIFT
        !                    Change the spin projection shift value.
        
        integer, parameter :: excite = 1, truncatecas = 2, softexit = 3, &
                              writepops = 4, varyshift = 5, nmcyc = 6, &
                              tau = 7, diagshift = 8, shiftdamp = 9, &
                              stepsshift = 10, singlesbias = 11, &
                              zeroproje = 12, zerohist = 13, &
                              partiallyfreeze = 14, partiallyfreezevirt = 15,&
                              printerrorblocking = 16, &
                              starterrorblocking = 17, &
                              restarterrorblocking = 18, &
                              printshiftblocking = 19, &
                              restartshiftblocking = 20, &
                              equilsteps = 21, starthist = 22, &
                              histequilsteps = 23, truncinitiator = 24, &
                              addtoinit = 25, scalehf = 26, &
                              printhighpopdet = 27, changerefdet = 28, &
                              restarthighpop = 29, spin_project = 30, &
                              spin_project_gamma = 31, spin_project_shift = 32

        integer, parameter :: last_item = spin_project_shift
        integer, parameter :: max_item_len = 25
        character(*), parameter :: option_list(last_item) &
                                = (/character(max_item_len) ::&
                                   "excite", &
                                   "truncatecas", &
                                   "softexit", &
                                   "writepops", &
                                   "varyshift", &
                                   "nmcyc", &
                                   "tau", &
                                   "diagshift", &
                                   "shiftdamp", &
                                   "stepsshift", &
                                   "singlesbias", &
                                   "zeroproje", &
                                   "zerohist", &
                                   "partiallyfreeze", &
                                   "partiallyfreezevirt", &
                                   "printerrorblocking", &
                                   "starterrorblocking", &
                                   "restarterrorblocking", &
                                   "printshiftblocking", &
                                   "restartshiftblocking", &
                                   "equilsteps", &
                                   "starthist", &
                                   "histequilsteps", &
                                   "truncinitiator", &
                                   "addtoinit", &
                                   "scalehf", &
                                   "printhighpopdet", &
                                   "changerefdet", &
                                   "restarthighpop", &
                                   "spin-project", &
                                   "spin-project-gamma", &
                                   "spin-project-shift"/)


        logical :: exists, any_exist, eof, deleted, any_deleted
        logical :: opts_selected(last_item)
        character(len=100) :: w
        integer :: i, proc

        ! Test if the changevars file exists, and broadcast to all nodes.
        inquire (file='CHANGEVARS', exist=exists)
        call MPIAllReduce (exists, 1, MPI_LOR, any_exists)

        opts_selected = .false.
        deleted = .false.
        if (any_exist) then
            if (iProcIndex == 0) &
                write (6, *) "CHANGEVARS file detected on iteration ", iter

            ! Each processor attemtps to delete changevars in turn. Wait for
            ! all processors to reach AllReduce on each cycle, to avoid race
            ! condition between processors sharing the same disk.
            do proc = 0, nProcessors - 1
                if (proc == iProcIndex .and. exists) then
                    ! Set unit for read_line routine
                    ir = get_free_unit
                    open (ir, file='CHANGEVARS', status='old', err=99, &
                          iostat=ios)
                    call input_options (echo_lines=.true., &
                                        skip_blank_lines=.true.)

                    ! Loop over all options specified in the file.
                    do
                        call read_line (eof)
                        if (eof) exit
                        call readl (w)

                        ! Mark any selected options.
                        do i = 1, last_item
                            if (trim(w) == trim(options_list(i))) then
                                opts_selected(i) = .true.
                                exit
                            endif
                        enddo
                    enddo

                    ! Close and delete the file
                    close (ir, status='delete')
                    deleted = .true.
                endif

                ! Once one node has found and deleted the file, it is gone.
                call MPIAllReduce (deleted, 1, MPI_LOR, any_deleted)
                if (any_deleted) exit
            enddo ! Loop to read CHANGEVARS

            ! Broadcast the selected options list to all processors
            call MPIBCast (opts_selected, last_item, proc)

            ! ***********************
            ! Now we need to deal with the specific options.
            ! ***********************
            ! TODO: Need to implement the readi/readf etc.

            ! Change excit level
            if (opts_selected(excite)) then
                if (.not. tTruncSpace) then
                    root_write 'The space is not truncated, so EXCITE &
                               &keyowrd in CHANGEVARS has no effect.'
                else
                    if (tHistSpawn .or. tCalcFCIMCPsi) then
                        root_write 'Cannot increase truncation level, since &
                                   &histogramming wavefunction.'
                    else
                        call MPIBcast (ICILevel, 1, proc)

                        if ((ICILevel < 0) .or. (ICILevel > nel)) then
                            tTruncSpace = .false.
                            root_write 'Expanding to the full space.'
                        else
                            root_write 'Increasing truncation level of space &
                                       &to ', ICILevel
                        endif
                    endif
                endif
            endif

            ! Change the CAS space
            if (opts_selected(truncatecas)) then
                if (.not. tTruncCAS) then
                    root_write 'The space is not truncated by CAS, so &
                               &TRUNCATECAS keyword in CHANGEVARS has no &
                               &effect'
                else

                    call MPIBCast (OccCASORbs, 1, proc)
                    call MPIBCast (VirtCASOrbs, 1, proc)

                    if ( ((occCASOrbs>nel) .and. (VirtCASOrbs>nBasis - nel)) &
                        .or. (occCASORbs < 0) .or. (VirtCASORbs < 0) ) then
                        ! CAS space is equal to or greater than the full 
                        ! space, or one of the arguments is less than zero.
                        tTruncCAS = .false.
                        root_write 'Expanding CAS to the full space.'
                    else
                        CASMax = nel + VirtCASOrbs
                        CASMin = nel - OccCASOrbs
                        root_write 'Increasing CAS space accessible to ', &
                                   OccCASORbs, ", ", VirtCASORbs
                    endif
                endif
            endif

            ! softexit
            if (opts_selected(softexit)) then
                tSoftExitFound = .true.
                root_write 'SOFTEXIT triggered. Exiting run.'
            endif

            ! Write POPS file
            if (opts_selected(writepops)) then
                tWritePopsFound = .true.
                root_write 'Asked to write out a popsfile'
            endif

            ! Enter variable shift mode
            if (opts_selected(varyshift)) then
                if (.not. tSinglePartPhase) then
                    root_write 'Request to vary shift denied. Already in &
                               &variable shift mode.'
                else
                    tSinglePartPhase = .false.
                    VaryShiftIter = iter
                    root_write 'Request to vary the shift detected on a node.'
                endif
            endif

            ! Change number of MC steps
            if (opts_selected(nmcyc)) then
                call MPIBCast (nmcyc_new, 1, proc)

                if (nmcyc_new < iter) then
                    root_write 'New value of NMCyc is LESS than the current &
                               & iteration number.'
                    root_write 'Therefore, the number of iterations has been &
                               & left at ', nmcyc_value
                else
                    nmcyc_value = nmcyc_new
                    root_write 'Total number of MC cycles set to ', &
                               nmcyc_value
                endif
            endif

            ! Change Tau
            if (opts_selected(tau)) then
                call MPIBCast (tau_value, 1, proc)
                root_write 'Tau changed to: ', tau_value
            endif

            ! Change the shift value
            if (opts_selected(diagshift)) then
                call MPIBCast (DiagSft, 1, proc)
                root_write 'DIAGSHIFT change to: ', DiagSft
            endif

            ! Change the shift damping parameter
            if (opts_selected(shiftdamp)) then
                call MPIBCast (SftDamp, 1, proc)
                root_write 'SHIFTDAMP changed to: ', SftDamp
            endif

            ! Change the shift update (and output) interval
            if (opts_selected(stepsshift)) then
                call MPIBCast (StepsSft, 1, proc)
                root_write 'STEPSSHIFT changed to: ', StepsSft
            endif

            ! Change the singles bias
            if (opts_selected(singlesbias)) then
                call MPIBcast (SinglesBias_value, 1, proc)
                tSingBiasChange = .true.
                root_write 'SINGLESBIAS changed to: ', SinglesBias
            endif

            ! Zero the average energy estimators
            if (opts_selected(zeroproje)) then
                SumENum = 0
                SumNoatHF = 0
                HFPopCyc = 0
                ProjEIterSum = 0
                VaryShiftCycles = 0
                SumDiagSft = 0
                root_write 'Zeroing all average energy estimators.'
            endif

            ! Zero average histograms
            if (opts_selected(zerohist)) then
                histogram = 0
                if (tHistSpawn) avAnnihil = 0
                root_write 'Zeroing all average histograms'
            endif

            ! Change the number of holes/electrons in the core valence region
            if (opts_selected(partiallyfreeze)) then
                call MPIBCast (NPartFrozen, 1, proc)
                call MPIBcast (NHolesFrozen, 1, proc)

                root_write 'Allowing ', nHolesFrozen, ' holes in ', &
                           nPartFrozen, ' partially frozen orbitals.'

                if (nHolesFrozen == nPartFrozen) then
                    ! Allowing as many holes as there are orbitals
                    !  --> equivalent to not freezing at all.
                    tPartFreezeCore = .false.
                    root_write 'Unfreezing any partially frozen core'
                else
                    tPartFreeze = .true.
                endif
            endif

            if (opts_selected(partiallyfreezevirt)) then
                call MPIBcast (nVirtPartFrozen, 1, proc)
                call MPIBcast (nelVirtFrozen, 1, proc)

                root_write 'Allowing ', nelVirtFrozen, ' electrons in ', &
                           nVirtPartFrozen, ' partially frozen virtual &
                          &orbitals.'
                if (nelVirtFrozen == nel) then
                    ! Allowing as many holes as there are orbitals
                    ! --> Equivalent ton not freezing at all
                    tPartFreezeVirt = .false.
                    root_write 'Unfreezing any partially frozen virtual &
                               &orbitals'
                else
                    tPartFreezeVirt = .true.
                endif
            endif
            
            ! Print blocking analysis here.
            if (opts_selected(printerrorblocking)) then
                root_write 'Printing blocking analysis at this point.'
                if (iprocindex == 0) call PrintBlocking (iter)
            endif

















        endif




       integer :: error,i,ios,NewNMCyc, pos
       integer, dimension(lenof_sign) :: HFSign
       logical :: tSoftExitFound,tWritePopsFound,exists,AnyExist,deleted_file
       logical :: tEof,any_deleted_file,tChangeParams(32),tSingBiasChange
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
                       case("CHANGEREFDET")
                           tChangeParams(28) = .true.
                       case("RESTARTHIGHPOP")
                           tChangeParams(29) = .true.
                       case("SPIN-PROJECT")
                           tChangeParams(30) = .true.
                           call readi(spin_proj_interval)
                           if (spin_proj_interval == 0) &
                               tSpinProject = .false.
                       case("SPIN-PROJECT-GAMMA")
                           tChangeParams(31) = .true.
                           call readf (spin_proj_gamma)
                       case("SPIN-PROJECT-SHIFT")
                           tChangeParams(32) = .true.
                           call readf (spin_proj_shift)
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
           CALL MPI_BCast(tChangeParams,32,MPI_LOGICAL,i,MPI_COMM_WORLD,error)
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
                   call extract_sign(CurrentDets(:,pos),HFSign)
                   do i=1,lenof_sign
                       !Multiply real, and if applicable, imaginary parts of the determinant.
                       HFSign(i)=HFSign(i) * hfScaleFactor
                   enddo
                   call encode_sign(CurrentDets(:,pos),HFSign)
               endif
           endif
           IF(tChangeParams(27)) THEN
               tPrintHighPop=.true.
#ifdef PARALLEL
               CALL MPI_BCast(tPrintHighPop,1,MPI_LOGICAL,i,MPI_COMM_WORLD,error)
#endif
               IF(iProcIndex.eq.0) WRITE(6,'(A)') 'Request to print the determinant with the highest population of different sign to the HF.'
           ENDIF   

           if (tChangeParams(28)) then
               tCheckHighestPopOnce = .true.
               tCheckHighestPop = .true.
               tChangeProjEDet = .true.
               FracLargerDet = 1.0
#ifdef PARALLEL
                call MPI_BCast(tCheckHighestPopOnce, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
                call MPI_BCast(tCheckHighestPop, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
                call MPI_BCast(tChangeProjEDet, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
                call MPI_BCast(FracLargerDet, 1, MPI_DOUBLE_PRECISION, i, &
                               MPI_COMM_WORLD, error)
#endif
                if (iProcIndex == 0) then
                    write (6, '(a)') "Request to change reference det."
                endif
           endif

           if (tChangeParams(29)) then
               tCheckHighestPopOnce = .true.
               tCheckHighestPop = .true.
               tRestartHighPop = .true.
               FracLargerDet = 1.0
#ifdef PARALLEL
                call MPI_BCast(tCheckHighestPopOnce, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
                call MPI_BCast(tCheckHighestPop, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
                call MPI_BCast(tRestartHighPop, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
                call MPI_BCast(FracLargerDet, 1, MPI_DOUBLE_PRECISION, i, &
                               MPI_COMM_WORLD, error)
#endif
                if (iProcIndex == 0) then
                    write (6, '(a)') "Request to restart with new reference det."
                endif
           endif

           if (tChangeParams(30)) then
               if (spin_proj_interval == 0) then
                   tSpinProject = .false.
               else
                   tSpinProject = .true.
               endif
#ifdef PARALLEL
               call MPI_BCast (spin_proj_interval, 1, MPI_INTEGER, i, &
                               MPI_COMM_WORLD, error)
               call MPI_BCast (tSPinProject, 1, MPI_LOGICAL, i, &
                               MPI_COMM_WORLD, error)
#endif
           endif

           if (tChangeParams(31)) then
#ifdef PARALLEL
               call MPI_BCast (spin_proj_gamma, 1, MPI_DOUBLE_PRECISION, &
                               i, MPI_COMM_WORLD, error)
#endif
           endif

           if (tChangeParams(32)) then
#ifdef PARALLEL
               call MPI_BCast (spin_proj_shift, 1, MPI_DOUBLE_PRECISION, &
                               i, MPI_COMM_WORLD, error)
#endif
           endif
 
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
