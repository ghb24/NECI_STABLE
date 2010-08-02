#include "macros.h"
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

! *******************************************************
! Supported options: (n.b. mutiple values may be changed at once)
! 
!   EXCITE  XXX          Will change the excitation level of the simulation 
!                        (< 0 or > NEl sets it to the full space)
!   TRUNCATECAS  XXX XXX Will change the CAS of the simulation (< 0 or > NEl 
!                        sets it to the full space)
!   SOFTEXIT             Exit cleanly from the program
!   WRITEPOPS            Write a current popsfile
!   VARYSHIFT            Exit fixed shift phase
!   NMCYC XXX            Change the number of monte carlo cycles to perform
!   TAU XXX              Change the value of tau for the simulation
!   DIAGSHIFT XXX        Change the shift
!   SHIFTDAMP XXX        Change the shift damping parameter
!   STEPSSHIFT XXX       Change the length of the update cycle
!   SINGLESBIAS XXX      Change the singles bias for the non-uniform random 
!                        excitation generator
!   ZEROPROJE            Re-zero the averaged energy estimators
!   ZEROHIST             Re-ezero the averaged histogramming vectors
!   PARTIALLYFREEZE XXX XXX 
!                        Change the number of holes/electrons in the core 
!                        valence region
!   PARTIALLYFREEZEVIRT XXX XXX
!                        Change the number of electrons in the partially 
!                        frozen virtual region
!   PRINTERRORBLOCKING   Print the blocking analysis
!   STARTERRORBLOCKING   Sart the blocking analysis
!   RESTARTERRORBLOCKING Restart the blocking analysis
!   PRINTSHIFTBLOCKING   Print the shift blocking analysis
!   RESTARTSHIFTBLOCKING Restart the shift blocking analysis
!   EQUILSTEPS XXX       Change the number of steps to ignore in the 
!                        averaging of the energy and the shift.
!   STARTHIST            Begin histogramming the determinant populations if 
!                        the tCalcFCIMCPsi is on and the histogramming has 
!                        been set up.
!   HISTEQUILSTEPS XXX   Change the iteration at which the histogramming 
!                        begins to the value specified.
!   TRUNCINITIATOR       Expand the CAS calculation to a TRUNCINITIATOR 
!                        calculation if DELAYTRUNCINITIATOR is present in 
!                        the input.
!   ADDTOINIT XXX        Change the cutt-off population for which walkers are
!                        added to the initiator space.  Pop must be *above* 
!                        specified value.
!   SCALEHF XXX          Scale the number of walkers at HF by the specified 
!                        factor
!   PRINTHIGHPOPDET      Print the determinant with the highest population of
!                        different sign to the HF.
!   CHANGEREFDET         Change the reference determinant to the det with the
!                        highest population
!   RESTARTHIGHPOP       Restart the calculation with same parameters but a 
!                        new reference determinant
!   SPIN-PROJECT         Change the interval between applications of 
!                        stochastic spin projection. If 0, disable it.
!                        If -1, disable FCIQMC propagation.
!   SPIN-PROJECT-GAMMA   Change the delta-gamma value used for stochastic
!                        spin projection
!   SPIN-PROJECT-SHIFT   Change the spin projection shift value.
! **********************************************************

module soft_exit

    use SystemData, only: nel, nBasis
    use bit_reps, only: NIfTot
    use util_mod, only: binary_search, get_free_unit
    use FciMCData, only: iter, CASMin, CASMax, tTruncSpace, tSinglePartPhase,&
                         SumENum, SumNoatHF, HFPopCyc, ProjEIterSum, &
                         Histogram, AvAnnihil, VaryShiftCycles, SumDiagSft, &
                         VaryShiftIter, CurrentDets, iLutHF, &
                         TotWalkers,tPrintHighPop
    use CalcData, only: DiagSft, SftDamp, StepsSft, OccCASOrbs, VirtCASOrbs, &
                        tTruncCAS,  NEquilSteps, tTruncInitiator, &
                        InitiatorWalkNo, tCheckHighestPop, tRestartHighPop, &
                        tChangeProjEDet, tCheckHighestPopOnce, FracLargerDet,&
                        SinglesBias_value => SinglesBias, tau_value => tau, &
                        nmcyc_value => nmcyc
    use DetCalcData, only: ICILevel
    use IntegralsData, only: tPartFreezeCore, NPartFrozen, NHolesFrozen, &
                             NVirtPartFrozen, NelVirtFrozen, tPartFreezeVirt
    use Input
    use Logging, only: tHistSpawn, tCalcFCIMCPsi, tIterStartBlock, &
                       IterStartBlocking, tHFPopStartBlock, NHistEquilSteps
    use FCIMCLoggingMOD, only: PrintBlocking, RestartBlocking, &
                               PrintShiftBlocking_proc => PrintShiftBlocking,&
                               RestartShiftBlocking_proc=>RestartShiftBlocking
    use AnnihilationMod, only: DetermineDetProc
    use constants, only: lenof_sign, int32, dp
    use bit_reps, only: extract_sign,encode_sign
    use spin_project, only: tSpinProject, spin_proj_gamma, &
                            spin_proj_interval, spin_proj_shift, &
                            spin_proj_cutoff
    use Parallel
    implicit none

contains

    subroutine ChangeVars (tSingBiasChange, tSoftExitFound, tWritePopsFound)

        ! Read CHANGEVARS file as described in module header.
        !
        ! Out: tSingBiasChange - true if the single bias is changed.
        !      tSoftExitFound  - true if a SOFTEXIT is requested.
        !      tWritePopsFound - true if the output of a POPSFILE has been 
        !                        requested.
        ! Other changes are made directly to the modules concerned
        
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
                              spin_project_gamma = 31, &
                              spin_project_shift = 32, &
                              spin_project_cutoff = 33

        integer, parameter :: last_item = spin_project_cutoff
        integer, parameter :: max_item_len = 25
        character(max_item_len), parameter :: option_list(last_item) &
                               = (/"excite                   ", &
                                   "truncatecas              ", &
                                   "softexit                 ", &
                                   "writepops                ", &
                                   "varyshift                ", &
                                   "nmcyc                    ", &
                                   "tau                      ", &
                                   "diagshift                ", &
                                   "shiftdamp                ", &
                                   "stepsshift               ", &
                                   "singlesbias              ", &
                                   "zeroproje                ", &
                                   "zerohist                 ", &
                                   "partiallyfreeze          ", &
                                   "partiallyfreezevirt      ", &
                                   "printerrorblocking       ", &
                                   "starterrorblocking       ", &
                                   "restarterrorblocking     ", &
                                   "printshiftblocking       ", &
                                   "restartshiftblocking     ", &
                                   "equilsteps               ", &
                                   "starthist                ", &
                                   "histequilsteps           ", &
                                   "truncinitiator           ", &
                                   "addtoinit                ", &
                                   "scalehf                  ", &
                                   "printhighpopdet          ", &
                                   "changerefdet             ", &
                                   "restarthighpop           ", &
                                   "spin-project             ", &
                                   "spin-project-gamma       ", &
                                   "spin-project-shift       ", &
                                   "spin-project-cutoff      "/)


        logical :: exists, any_exist, eof, deleted, any_deleted
        logical :: opts_selected(last_item)
        logical, intent(out) :: tSingBiasChange, tSoftExitFound
        logical, intent(out) :: tWritePopsFound
        integer :: i, proc, nmcyc_new, ios, pos
        integer, dimension(lenof_sign) :: hfsign
        real(dp) :: hfScaleFactor
        character(len=100) :: w

        ! Test if the changevars file exists, and broadcast to all nodes.
        inquire (file='CHANGEVARS', exist=exists)
        call MPIAllReduce (exists, 1, MPI_LOR, any_exist)

        ! Default values
        opts_selected = .false.
        deleted = .false.
        tSoftExitFound = .false.
        tWritePopsfound = .false.
        tSingBiasChange = .false.
        ios = 0

        if (any_exist) then
            if (iProcIndex == 0) &
                write (6, *) "CHANGEVARS file detected on iteration ", iter

            ! Each processor attemtps to delete changevars in turn. Wait for
            ! all processors to reach AllReduce on each cycle, to avoid race
            ! condition between processors sharing the same disk.
            do proc = 0, nProcessors - 1
                if (proc == iProcIndex .and. exists) then
                    ! Set unit for read_line routine
                    ir = get_free_unit ()
                    open (ir, file='CHANGEVARS', status='old', iostat=ios)
                    if (ios /= 0) then
                        write (6, *) 'Problem reading CHANGEVARS file.'
                        cycle
                    endif
                    call input_options (echo_lines=.true., &
                                        skip_blank_lines=.true.)

                    ! Loop over all options specified in the file.
                    do
                        call read_line (eof)
                        if (eof) exit
                        call readl (w)

                        ! Mark any selected options.
                        do i = 1, last_item
                            if (trim(w) == trim(option_list(i))) then
                                opts_selected(i) = .true.
                                exit
                            endif
                        enddo

                        ! Do we have any other items to read in?
                        if (i == tau) then
                            call readf (tau_value)
                        elseif (i == diagshift) then
                            call readf (DiagSft)
                        elseif (i == shiftdamp) then
                            call readf (SftDamp)
                        elseif (i == stepsshift) then
                            call readi (StepsSft)
                        elseif (i == excite) then
                            call readi (ICILevel)
                        elseif (i == singlesbias) then
                            call readf (singlesbias_value)
                        elseif (i == truncatecas) then
                            call readi (OccCASOrbs)
                            call readi (VirtCASOrbs)
                        elseif (i == nmcyc) then
                            call readi (nmcyc_new)
                        elseif (i == partiallyfreeze) then
                            call readi (nPartFrozen)
                            call readi (nHolesFrozen)
                        elseif (i == equilsteps) then
                            call readi (nEquilSteps)
                        elseif (i == histequilsteps) then
                            call readi (nHistEquilSteps)
                        elseif (i == partiallyfreezevirt) then
                            call readi (nVirtPartFrozen)
                            call readi (nElVirtFrozen)
                        elseif (i == addtoinit) then
                            call readi (InitiatorWalkNo)
                        elseif (i == scalehf) then
                            call readf (hfScaleFactor)
                        elseif (i == spin_project) then
                            call readi (spin_proj_interval)
                        elseif (i == spin_project_gamma) then
                            call readf (spin_proj_gamma)
                        elseif (i == spin_project_shift) then
                            call readf (spin_proj_shift)
                        elseif (i == spin_project_cutoff) then
                            call readi (spin_proj_cutoff)
                        endif
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

            ! Change excit level
            if (opts_selected(excite)) then
                if (.not. tTruncSpace) then
                    root_print 'The space is not truncated, so EXCITE &
                               &keyword in CHANGEVARS has no effect.'
                else
                    if (tHistSpawn .or. tCalcFCIMCPsi) then
                        root_print 'Cannot increase truncation level, since &
                                   &histogramming wavefunction.'
                    else
                        call MPIBcast (ICILevel, 1, proc)

                        if ((ICILevel < 0) .or. (ICILevel > nel)) then
                            tTruncSpace = .false.
                            root_print 'Expanding to the full space.'
                        else
                            root_print 'Increasing truncation level of space &
                                       &to ', ICILevel
                        endif
                    endif
                endif
            endif

            ! Change the CAS space
            if (opts_selected(truncatecas)) then
                if (.not. tTruncCAS) then
                    root_print 'The space is not truncated by CAS, so &
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
                        root_print 'Expanding CAS to the full space.'
                    else
                        CASMax = nel + VirtCASOrbs
                        CASMin = nel - OccCASOrbs
                        root_print 'Increasing CAS space accessible to ', &
                                   OccCASORbs, ", ", VirtCASORbs
                    endif
                endif
            endif

            ! softexit
            if (opts_selected(softexit)) then
                tSoftExitFound = .true.
                root_print 'SOFTEXIT triggered. Exiting run.'
            endif

            ! Write POPS file
            if (opts_selected(writepops)) then
                tWritePopsFound = .true.
                root_print 'Asked to write out a popsfile'
            endif

            ! Enter variable shift mode
            if (opts_selected(varyshift)) then
                if (.not. tSinglePartPhase) then
                    root_print 'Request to vary shift denied. Already in &
                               &variable shift mode.'
                else
                    tSinglePartPhase = .false.
                    VaryShiftIter = iter
                    root_print 'Request to vary the shift detected on a node.'
                endif
            endif

            ! Change number of MC steps
            if (opts_selected(nmcyc)) then
                call MPIBCast (nmcyc_new, 1, proc)

                if (nmcyc_new < iter) then
                    root_print 'New value of NMCyc is LESS than the current &
                               & iteration number.'
                    root_print 'Therefore, the number of iterations has been &
                               & left at ', nmcyc_value
                else
                    nmcyc_value = nmcyc_new
                    root_print 'Total number of MC cycles set to ', &
                               nmcyc_value
                endif
            endif

            ! Change Tau
            if (opts_selected(tau)) then
                call MPIBCast (tau_value, 1, proc)
                root_print 'Tau changed to: ', tau_value
            endif

            ! Change the shift value
            if (opts_selected(diagshift)) then
                call MPIBCast (DiagSft, 1, proc)
                root_print 'DIAGSHIFT change to: ', DiagSft
            endif

            ! Change the shift damping parameter
            if (opts_selected(shiftdamp)) then
                call MPIBCast (SftDamp, 1, proc)
                root_print 'SHIFTDAMP changed to: ', SftDamp
            endif

            ! Change the shift update (and output) interval
            if (opts_selected(stepsshift)) then
                call MPIBCast (StepsSft, 1, proc)
                root_print 'STEPSSHIFT changed to: ', StepsSft
            endif

            ! Change the singles bias
            if (opts_selected(singlesbias)) then
                call MPIBcast (SinglesBias_value, 1, proc)
                tSingBiasChange = .true.
                root_print 'SINGLESBIAS changed to: ', SinglesBias
            endif

            ! Zero the average energy estimators
            if (opts_selected(zeroproje)) then
                SumENum = 0
                SumNoatHF = 0
                HFPopCyc = 0
                ProjEIterSum = 0
                VaryShiftCycles = 0
                SumDiagSft = 0
                root_print 'Zeroing all average energy estimators.'
            endif

            ! Zero average histograms
            if (opts_selected(zerohist)) then
                histogram = 0
                if (tHistSpawn) avAnnihil = 0
                root_print 'Zeroing all average histograms'
            endif

            ! Change the number of holes/electrons in the core valence region
            if (opts_selected(partiallyfreeze)) then
                call MPIBCast (NPartFrozen, 1, proc)
                call MPIBcast (NHolesFrozen, 1, proc)

                root_print 'Allowing ', nHolesFrozen, ' holes in ', &
                           nPartFrozen, ' partially frozen orbitals.'

                if (nHolesFrozen == nPartFrozen) then
                    ! Allowing as many holes as there are orbitals
                    !  --> equivalent to not freezing at all.
                    tPartFreezeCore = .false.
                    root_print 'Unfreezing any partially frozen core'
                else
                    tPartFreezeCore = .true.
                endif
            endif

            if (opts_selected(partiallyfreezevirt)) then
                call MPIBcast (nVirtPartFrozen, 1, proc)
                call MPIBcast (nelVirtFrozen, 1, proc)

                root_print 'Allowing ', nelVirtFrozen, ' electrons in ', &
                           nVirtPartFrozen, ' partially frozen virtual &
                          &orbitals.'
                if (nelVirtFrozen == nel) then
                    ! Allowing as many holes as there are orbitals
                    ! --> Equivalent ton not freezing at all
                    tPartFreezeVirt = .false.
                    root_print 'Unfreezing any partially frozen virtual &
                               &orbitals'
                else
                    tPartFreezeVirt = .true.
                endif
            endif
            
            ! Print blocking analysis here.
            if (opts_selected(printerrorblocking)) then
                root_print 'Printing blocking analysis at this point.'
                if (iprocindex == 0) call PrintBlocking (iter)
            endif

            ! Start blocking analysis
            if (opts_selected(starterrorblocking)) then
                if ((.not.tHFPopStartBlock) .and. (.not.tIterStartBlock)) then
                    root_print 'Error blocking already started'
                else
                    tIterStartBlock = .true.
                    IterStartBlocking = iter
                endif
            endif

            ! Restart error blocking
            if (opts_selected(restarterrorblocking)) then
                root_print 'Restarting the error calculations. All blocking &
                           &arrays are re-set to zero.'
                if (iProcIndex == 0) call RestartBlocking (iter)
            endif

            ! Print shift blocking analysis here
            if (opts_selected(printshiftblocking)) then
                root_print 'Printing shift error blocking.'
                if (iProcIndex == 0) call PrintShiftBlocking_proc (iter)
            endif

            ! Restart shift blocking analysis
            if (opts_selected(restartshiftblocking)) then
                root_print 'Restarting the shift error calculations. All &
                           &shift blocking arrays set to zero.'
                if (iProcIndex == 0) call RestartShiftBlocking_proc (iter)
            endif

            ! Change the number of equilibration steps
            if (opts_selected(equilsteps)) then
                call MPIBcast (nEquilSteps, 1, proc)
                root_print 'Changing the number of equilibration steps to ', &
                           nEquilSteps
            endif

            ! Start histogramming
            if (opts_selected(starthist)) then
                root_print 'Beginning to histogram at the next update'
                if (iProcIndex == 0) nHistEquilSteps = iter + StepsSft
                if (.not. tCalcFCIMCPsi) then
                    root_print 'This has no effect, as the histograms have &
                               &not been set up at the beginning of the &
                               &calculation.'
                endif
            endif

            ! Change the starting iteration for histogramming
            if (opts_selected(histequilsteps)) then
                if (nHistEquilSteps < iter) nHistEquilSteps = iter + StepsSft
                root_print 'Changing the starting iteration for &
                           &histogramming to ', nHistEquilSteps
                if (.not. tCalcFCIMCPsi) then
                    root_print 'This has no effec, as the histograms have &
                               &not been set up at the beginning of the &
                               &calculation.'
                endif
            endif

            ! Enable initiator truncation scheme
            if (opts_selected(truncinitiator)) then
                tTruncInitiator = .true.
                tau_value = tau_value / 10 ! Done by all. No need to BCast...
                root_print 'Beginning to allow spawning into inactive space &
                           &for a truncated initiator calculation.'
                root_print 'Reducing tau by an order of magnitude. The new &
                           &tau is ', tau_value
            endif

            ! Change the initiator cutoff parameter
            if (opts_selected(addtoinit)) then
                call MPIBCast (InitiatorWalkNo, 1, proc)
                root_print 'Cutoff propulation for determinants to be added &
                           &to the initiator space changed to ', &
                           InitiatorWalkNo
            endif
            
            ! Scale the number of walkers on the HF det
            if (opts_selected(scalehf)) then
                call MPIBcast (HFScaleFactor, 1, proc)
                root_print 'Number at Hartree-Fock scaled by factor: ', &
                           hfScaleFactor

                SumNoatHF = SumNoatHF * hfScaleFactor
                if (iProcIndex == DetermineDetProc(ilutHF)) then
                    pos = binary_search (CurrentDets, iLutHF, NIfTot+1, &
                                         int(TotWalkers,int32))
                    call extract_sign (CurrentDets(:,pos), hfsign)
                    do i = 1, lenof_sign
                        hfsign(i) = hfsign(i) * hfScaleFactor
                    enddo
                    call encode_sign (CurrentDets(:,pos), HFSign)
                endif
            endif

            ! Print the determinants with the largest +- populations
            if (opts_selected(printhighpopdet)) then
                tPrintHighPop = .true.
                root_print 'Request to print the determinants with the &
                           &largest populations detected.'
            endif

            ! Change the reference determinant on the fly
            if (opts_selected(changerefdet)) then
                tCheckHighestPopOnce = .true.
                tCheckHighestPop = .true.
                tChangeProjEDet = .true.
                FracLargerDet = 1.0
                root_print 'Changing the reference determinant to the most &
                           &highly weighted determinants.'
            endif

            ! Restart with new reference determinant
            if (opts_selected(restarthighpop)) then
                tCheckHighestPopOnce = .true.
                tCheckHighestPop = .true.
                tRestartHighPop = .true.
                FracLargerDet = 1.0
                root_print 'Restarting the calculation with the most highly &
                           &weighted determinant as the reference determiant.'
            endif

            ! Enable spin projection, and change application interval
            if (opts_selected(spin_project)) then
                call MPIBcast (spin_proj_interval, 1, proc)
                if (spin_proj_interval == 0) then
                    tSpinProject = .false.
                    root_print 'Stochastic spin projection disabled'
                else
                    tSpinProject = .true.
                    root_print 'Stochastic spin projection applied every ', &
                               spin_proj_interval, ' iterations.'
                endif
            endif

            ! Change delta-gamma for spin projection
            if (opts_selected(spin_project_gamma)) then
                call MPIBcast (spin_proj_gamma, 1, proc)
                root_print 'Changed gamma value for spin projection to ', &
                           spin_proj_gamma
            endif

            ! Change shift value for spin projection
            if (opts_selected(spin_project_shift)) then
                call MPIBcast (spin_proj_shift, 1, proc)
                root_print 'Changed shift value for spin projection to ', &
                           spin_proj_shift
            endif

            ! Change walker number cutoff value for spin projection
            if (opts_selected(spin_project_shift)) then
                call MPIBcast (spin_proj_cutoff, 1, proc)
                root_print 'Changed walker number cutoff value for spin &
                           &projection to ', spin_proj_shift
            endif
        endif

    end subroutine ChangeVars

logical function test_SoftExit()

    logical :: tdummy1, tdummy2

    call ChangeVars(tdummy1, test_SoftExit, tdummy2)
    if (test_SoftExit) write (6,'(1X,a30)') 'Request for SOFTEXIT detected.'

end function test_SoftExit

end module soft_exit
