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
!   TARGETGROWRATE XXX   Change the target growthrate for the simulation
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
!   REFSHIFT             Change the default use of the shift to now keep HF populations constant.
!   CALCRDMONFLY XXX XXX XXX
!                        Stochastically calculate the reduced density
!                        matrices.  The first integer specifies the
!                        XXX-electron RDM (3 for both 1 and 2).  The second
!                        is the number of iterations after the shift starts
!                        changing, to start filling the RDM, and the
!                        third is the frequency the energy is calculated
!                        and printed.
!   CALCEXPLICITRDM XXX XXX XXX
!                        Same as above, but the RDM is filled using the
!                        explicit algorithm.
!   FILLRDMITER XXX      Change the number of iterations after the shift has
!                        changed that the RDM are filled from.
!   DIAGFLYONERDM        Requests to diagonalise the 1-RDM at the end.
!   REFSHIFT             Change the default use of the shift to now keep HF
!                        populations constant.
!   PREPARE_REAL_TIME n m   start the print out of the current walker population
!                        n times with m cycles between them, which in turn
!                        will be used to start a subsequent real-time
!                        calculations with these popsfile as the groundstate
!   TIME                 Specify a total elapsed-time before the calculation
!                        performs an automatic soft-exit. If specified as -1,
!                        we don't stop automatically.
! **********************************************************

module soft_exit

    use SystemData, only: nel, nBasis, tHPHF
    use bit_reps, only: NIfTot
    use util_mod, only: binary_search, get_free_unit
    use FciMCData, only: iter, CASMin, CASMax, tTruncSpace, tSinglePartPhase,&
                         SumENum, SumNoatHF, tTimeExit, &
                         AvAnnihil, VaryShiftCycles, SumDiagSft, &
                         VaryShiftIter, CurrentDets, iLutHF, HFDet, &
                         TotWalkers,tPrintHighPop, MaxTimeExit, &
                         proje_iter
    use CalcData, only: DiagSft, SftDamp, StepsSft, OccCASOrbs, VirtCASOrbs, &
                        tTruncCAS,  NEquilSteps, tTruncInitiator, &
                        InitiatorWalkNo, tCheckHighestPop, tRestartHighPop, &
                        tChangeProjEDet, tCheckHighestPopOnce, FracLargerDet,&
                        SinglesBias_value => SinglesBias, &
                        nmcyc_value => nmcyc, tTruncNOpen, trunc_nopen_max, &
                        target_grow_rate => TargetGrowRate, tShiftonHFPop, &
                        tAllRealCoeff, tRealSpawnCutoff, tJumpShift
    use tau_search, only: tau_search_method, possible_tau_search_methods, &
        tau_value => tau, assign_value_to_tau
    use tau_search_hist, only: frq_ratio_cutoff, t_fill_frequency_hists
    use DetCalcData, only: ICILevel
    use IntegralsData, only: tPartFreezeCore, NPartFrozen, NHolesFrozen, &
                             NVirtPartFrozen, NelVirtFrozen, tPartFreezeVirt
    use input_parser_mod, only: ManagingFileReader_t, TokenIterator_t
    Use LoggingData, only: tCalcFCIMCPsi, tIterStartBlock, &
                       IterStartBlocking, tHFPopStartBlock, NHistEquilSteps, &
                       IterRDMonFly_value => IterRDMonFly, RDMExcitLevel, &
                       tExplicitAllRDM, tRDMonFly, tChangeVarsRDM, &
                       RDMEnergyIter, tDiagRDM, tPopsFile, tPrintPopsDefault, &
                       tIncrementPops, iWritePopsEvery
    use FCIMCLoggingMOD, only: PrintBlocking, RestartBlocking, &
                               PrintShiftBlocking_proc => PrintShiftBlocking,&
                               RestartShiftBlocking_proc=>RestartShiftBlocking
    use constants, only: lenof_sign, int32, dp
    use bit_rep_data, only: extract_sign
    use bit_reps, only: encode_sign
    use load_balance_calcnodes, only: DetermineDetNode
    use hist_data, only: Histogram, tHistSpawn
    use Parallel_neci
    use fortran_strings, only: to_lower, to_int, to_realdp

    implicit none

    logical, volatile :: tSoftExitFound = .false.

contains

    subroutine ChangeVars (tSingBiasChange, tWritePopsFound)

        ! Read CHANGEVARS file as described in module header.
        !
        ! Out: tSingBiasChange - true if the single bias is changed.
        !      tSoftExitFound  - true if a SOFTEXIT is requested.
        !      tWritePopsFound - true if the output of a POPSFILE has been
        !                        requested.
        ! Other changes are made directly to the modules concerned

        integer, parameter :: excite                =  1, &
                              truncatecas           =  2, &
                              softexit              =  3, &
                              writepops             =  4, &
                              varyshift             =  5, &
                              nmcyc                 =  6, &
                              tau                   =  7, &
                              diagshift             =  8, &
                              shiftdamp             =  9, &
                              stepsshift            = 10, &
                              singlesbias           = 11, &
                              zeroproje             = 12, &
                              zerohist              = 13, &
                              partiallyfreeze       = 14, &
                              partiallyfreezevirt   = 15, &
                              printerrorblocking    = 16, &
                              starterrorblocking    = 17, &
                              restarterrorblocking  = 18, &
                              printshiftblocking    = 19, &
                              restartshiftblocking  = 20, &
                              equilsteps            = 21, &
                              starthist             = 22, &
                              histequilsteps        = 23, &
                              truncinitiator        = 24, &
                              addtoinit             = 25, &
                              scalehf               = 26, &
                              printhighpopdet       = 27, &
                              changerefdet          = 28, &
                              restarthighpop        = 29, &
                              trunc_nopen           = 30, &
                              targetgrowrate        = 31, &
                              refshift              = 32, &
                              calc_rdm              = 33, &
                              calc_explic_rdm       = 34, &
                              fill_rdm_iter         = 35, &
                              diag_one_rdm          = 36, &
                              frequency_cutoff      = 37, & !for the histogram integration
                              time                  = 38

        integer, parameter :: last_item = time
        integer, parameter :: max_item_len = 30
        character(max_item_len), parameter :: option_list_molp(last_item) &
                               = (/"truncate                     ", &
                                   "not_option                   ", &
                                   "exit                         ", &
                                   "writepopsfile                ", &
                                   "varyshift                    ", &
                                   "iterations                   ", &
                                   "timestep                     ", &
                                   "shift                        ", &
                                   "shiftdamping                 ", &
                                   "interval                     ", &
                                   "singlesbias                  ", &
                                   "zeroproje                    ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "initiator_thresh             ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "changeref                    ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   ", &
                                   "not_option                   "/)

        character(max_item_len), parameter :: option_list(last_item) &
                               = (/"excite                       ", &
                                   "truncatecas                  ", &
                                   "softexit                     ", &
                                   "writepops                    ", &
                                   "varyshift                    ", &
                                   "nmcyc                        ", &
                                   "tau                          ", &
                                   "diagshift                    ", &
                                   "shiftdamp                    ", &
                                   "stepsshift                   ", &
                                   "singlesbias                  ", &
                                   "zeroproje                    ", &
                                   "zerohist                     ", &
                                   "partiallyfreeze              ", &
                                   "partiallyfreezevirt          ", &
                                   "printerrorblocking           ", &
                                   "starterrorblocking           ", &
                                   "restarterrorblocking         ", &
                                   "printshiftblocking           ", &
                                   "restartshiftblocking         ", &
                                   "equilsteps                   ", &
                                   "starthist                    ", &
                                   "histequilsteps               ", &
                                   "truncinitiator               ", &
                                   "addtoinit                    ", &
                                   "scalehf                      ", &
                                   "printhighpopdet              ", &
                                   "changerefdet                 ", &
                                   "restarthighpop               ", &
                                   "trunc-nopen                  ", &
                                   "targetgrowrate               ", &
                                   "refshift                     ", &
                                   "calcrdmonfly                 ", &
                                   "calcexplicitrdm              ", &
                                   "fillrdmiter                  ", &
                                   "diagflyonerdm                ", &
                                   "frequency-cutoff             ", &
                                   "time                         "/)

        ! Logical(4) datatypes for compilation with builds of openmpi that don't
        ! have support for logical(8). Gah.
        logical :: deleted, any_deleted, opts_selected(last_item)
        logical :: exists, any_exist

        logical :: tSource
        logical, intent(out) :: tSingBiasChange
        logical, intent(out) :: tWritePopsFound
        real(dp), dimension(lenof_sign) :: hfsign
        integer :: i, proc, nmcyc_new, ios, pos, trunc_nop_new, IterRDMonFly_new, run
        real(dp) :: hfScaleFactor
        character(len=100) :: w
        character(*), parameter :: file_name = 'CHANGEVARS'
        type(ManagingFileReader_t) :: file_reader
        type(TokenIterator_t) :: tokens

        ! Test if the changevars file exists, and broadcast to all nodes.
        any_exist=.false.
        inquire (file=file_name, exist=exists)
        call MPIAllLorLogical(exists, any_exist)

        ! Default values
        opts_selected = .false.
        deleted = .false.
        any_deleted=.false.
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
                    ! Instantiate the file-reader object and flag whether
                    ! there is a file-opening error by unsetting "exists".
                    file_reader = ManagingFileReader_t(file_name, err=ios)
                    if (ios /= 0) then
                        write(stdout, *) 'Problem reading CHANGEVARS file.'
                        write(stderr, *) 'Problem reading CHANGEVARS file.'
                        exists = .false.
                    endif
                end if

                ! Skip parsing if there was a file-opening error, but still
                ! run the later MPI call so we don't get stuck.
                if (proc == iProcIndex .and. exists) then

                    ! Loop over all options specified in the file.
                    do while (file_reader%nextline(tokens, skip_empty=.true.))
                        w = to_lower(tokens%next())

                        ! Mark any selected options.
                        do i = 1, last_item
                            if ((trim(w) == trim(option_list(i))).or.(trim(w).eq.trim(option_list_molp(i)))) then
                                opts_selected(i) = .true.
                                exit
                            endif
                        enddo

                        ! Do we have any other items to read in?
                        if (i == tau) then
                            call assign_value_to_tau(to_realdp(tokens%next()), 'Manual change via `CHANGEVARS` file.')
                        elseif (i == TargetGrowRate) then
                            target_grow_rate(1) = to_realdp(tokens%next())
                            if(inum_runs == 2) target_grow_rate(inum_runs)=target_grow_rate(1)
                        elseif (i == diagshift) then
                            DiagSft(1) = to_realdp(tokens%next())
                            if(inum_runs == 2) DiagSft(inum_runs)=DiagSft(1)
                        elseif (i == shiftdamp) then
                            SftDamp = to_realdp(tokens%next())
                        elseif (i == stepsshift) then
                            StepsSft = to_int(tokens%next())
                        elseif (i == excite) then
                            ICILevel = to_int(tokens%next())
                        elseif (i == singlesbias) then
                            singlesbias_value = to_realdp(tokens%next())
                        elseif (i == truncatecas) then
                            OccCASOrbs = to_int(tokens%next())
                            VirtCASOrbs = to_int(tokens%next())
                        elseif (i == nmcyc) then
                            nmcyc_new = to_int(tokens%next())
                        elseif (i == partiallyfreeze) then
                            nPartFrozen = to_int(tokens%next())
                            nHolesFrozen = to_int(tokens%next())
                        elseif (i == equilsteps) then
                            nEquilSteps = to_int(tokens%next())
                        elseif (i == histequilsteps) then
                            nHistEquilSteps = to_int(tokens%next())
                        elseif (i == partiallyfreezevirt) then
                            nVirtPartFrozen = to_int(tokens%next())
                            nElVirtFrozen = to_int(tokens%next())
                        elseif (i == addtoinit) then
                            InitiatorWalkNo = to_realdp(tokens%next())
                        elseif (i == scalehf) then
                            hfScaleFactor = to_realdp(tokens%next())
                        elseif (i == trunc_nopen) then
                            trunc_nop_new = to_int(tokens%next())
                        elseif (i == calc_rdm) then
                            RDMExcitLevel = to_int(tokens%next())
                            IterRDMonFly_new = to_int(tokens%next())
                            RDMEnergyIter = to_int(tokens%next())
                        elseif (i == calc_explic_rdm) then
                            RDMExcitLevel = to_int(tokens%next())
                            IterRDMonFly_new = to_int(tokens%next())
                            RDMEnergyIter = to_int(tokens%next())
                        elseif (i == fill_rdm_iter) then
                            IterRDMonFly_new = to_int(tokens%next())
                        elseif (i == frequency_cutoff) then
                            frq_ratio_cutoff = to_realdp(tokens%next())
                        elseif (i == time) then
                            MaxTimeExit = to_realdp(tokens%next())
                        endif
                    enddo
                    call file_reader%close(delete=.true.)
                    deleted = .true.
                endif

                ! Once one node has found and deleted the file, it is gone.
                any_deleted=.false.
                call MPIAllLORLogical(deleted, any_deleted)
                if (any_deleted) exit
            enddo ! Loop to read CHANGEVARS

            ! Do not proceed further if read errors have prevented loading
            ! the contents of the file on all processes.
            if (.not.any_deleted) return

            ! Relabel 'deleted' as 'tSource' for clarity
            ! --> If we have had the file, we should be the source node
            tSource = deleted

            ! Broadcast the selected options list to all processors
            call MPIBCast (opts_selected, tSource)

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
                        call MPIBcast (ICILevel, tSource)

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

                    call MPIBCast (OccCASORbs, tSource)
                    call MPIBCast (VirtCASOrbs, tSource)

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
                write(stdout,*) 'Asked to write out a popsfile on iteration: ',iter
            endif

            ! Enter variable shift mode
            if (opts_selected(varyshift)) then
                do run=1,inum_runs
                    if (.not. tSinglePartPhase(run)) then
                        root_print 'Request to vary shift denied. Already in &
                                   &variable shift mode.'
                    else
                        tSinglePartPhase(run) = .false.
                        VaryShiftIter(run) = iter
                        write(stdout,*) 'Request to vary the shift detected on a node on iteration: ',iter

                        ! If specified, jump the value of the shift to that
                        ! predicted by the projected energy
                        if (tJumpShift) &
                            DiagSft(run) = proje_iter(run)
                    endif
                enddo
            endif

            ! Change number of MC steps
            if (opts_selected(nmcyc)) then
                call MPIBCast (nmcyc_new, tSource)

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
                block
                    real(dp) :: local_tau
                    local_tau = tau_value
                    call MPIBCast (local_tau, tSource)
                    call assign_value_to_tau(local_tau, 'Manual change via `CHANGEVARS` file.')
                end block
                if (tau_search_method /= possible_tau_search_methods%off) then
                    if (tau_search_method == possible_tau_search_methods%HISTOGRAMMING) then
                        t_fill_frequency_hists = .false.
                    end if
                    write(stdout, *)
                    write(stdout, *) 'The current tau search method is: ', trim(tau_search_method%str)
                    write(stdout, *) 'It is switched off now.'
                    write(stdout, *)
                    tau_search_method = possible_tau_search_methods%OFF
                end if
            endif

            if(opts_selected(targetgrowrate)) then
                call MPIBCast(target_grow_rate, tSource)
                write(stdout,*) "TARGETGROWRATE changed to: ",target_grow_rate, "on iteration: ",iter
            endif

            ! Change the shift value
            if (opts_selected(diagshift)) then
                call MPIBCast (DiagSft, tSource)
                write(stdout,*) 'DIAGSHIFT changed to: ', DiagSft, 'on iteration: ',iter
            endif

            ! Change the shift damping parameter
            if (opts_selected(shiftdamp)) then
                call MPIBCast (SftDamp, tSource)
                write(stdout,*) 'SHIFTDAMP changed to: ', SftDamp, 'on iteration: ',iter
            endif

            ! Change the shift update (and output) interval
            if (opts_selected(stepsshift)) then
                call MPIBCast (StepsSft, tSource)
                write(stdout,*) 'STEPSSHIFT changed to: ', StepsSft, 'on iteration: ',iter
            endif

            ! Change the singles bias
            if (opts_selected(singlesbias)) then
                call MPIBcast (SinglesBias_value, tSource)
                tSingBiasChange = .true.
                write(stdout,*) 'SINGLESBIAS changed to: ', SinglesBias, 'on iteration: ',iter
            endif

            ! Zero the average energy estimators
            if (opts_selected(zeroproje)) then
                SumENum = 0
                SumNoatHF = 0
                VaryShiftCycles = 0
                SumDiagSft = 0
                write(stdout,*) 'Zeroing all average energy estimators on iteration: ',iter
            endif

            ! Zero average histograms
            if (opts_selected(zerohist)) then
                histogram = 0
                if (tHistSpawn) avAnnihil = 0
                root_print 'Zeroing all average histograms'
            endif

            ! Change the number of holes/electrons in the core valence region
            if (opts_selected(partiallyfreeze)) then
                call MPIBCast (NPartFrozen, tSource)
                call MPIBcast (NHolesFrozen, tSource)

                write(stdout,*) 'Allowing ', nHolesFrozen, ' holes in ', &
                           nPartFrozen, ' partially frozen orbitals on iteration: ',iter

                if (nHolesFrozen == nPartFrozen) then
                    ! Allowing as many holes as there are orbitals
                    !  --> equivalent to not freezing at all.
                    tPartFreezeCore = .false.
                    write(stdout,*) 'Unfreezing any partially frozen core on iteration: ',iter
                else
                    tPartFreezeCore = .true.
                endif
            endif

            if (opts_selected(partiallyfreezevirt)) then
                call MPIBcast (nVirtPartFrozen, tSource)
                call MPIBcast (nelVirtFrozen, tSource)

                write(stdout,*) 'Allowing ', nelVirtFrozen, ' electrons in ', &
                           nVirtPartFrozen, ' partially frozen virtual &
                          &orbitals on iteration: ',iter
                if (nelVirtFrozen == nel) then
                    ! Allowing as many holes as there are orbitals
                    ! --> Equivalent ton not freezing at all
                    tPartFreezeVirt = .false.
                    write(stdout,*) 'Unfreezing any partially frozen virtual &
                               &orbitals on iteration: ',iter
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
                write(stdout,*) 'Restarting the error calculations. All blocking &
                           &arrays are re-set to zero on iteration: ',iter
                if (iProcIndex == 0) call RestartBlocking (iter)
            endif

            ! Print shift blocking analysis here
            if (opts_selected(printshiftblocking)) then
                write(stdout,*) 'Printing shift error blocking on iteration: ',iter
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
                call MPIBcast (nEquilSteps, tSource)
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
                call assign_value_to_tau( &
                    tau_value / 10, &
                    'Beginning to allow spawning into inactive space &
                     &for a truncated initiator calculation. &
                    &Reducing tau by an order of magnitude.')
            endif

            ! Change the initiator cutoff parameter
            if (opts_selected(addtoinit)) then
                call MPIBCast (InitiatorWalkNo, tSource)
                root_print 'Cutoff propulation for determinants to be added &
                           &to the initiator space changed to ', &
                           InitiatorWalkNo
            endif

            ! Scale the number of walkers on the HF det
            if (opts_selected(scalehf)) then
                call MPIBcast (HFScaleFactor, tSource)
                root_print 'Number at Hartree-Fock scaled by factor: ', &
                           hfScaleFactor

                SumNoatHF = nint(real(SumNoatHF,dp) * hfScaleFactor,int64)
                if (iNodeIndex == DetermineDetNode(nel,HFDet,0).and. bNodeRoot) then
                    pos = binary_search (CurrentDets(:,1:TotWalkers), &
                                         iLutHF)
                    call extract_sign (CurrentDets(:,pos), hfsign)
                    do i = 1, lenof_sign
                        hfsign(i) = hfsign(i) * hfScaleFactor
                        if (.not. (tAllRealCoeff .or. tRealSpawnCutoff)) &
                            hfsign = nint(hfsign)
                    enddo
                    call encode_sign (CurrentDets(:,pos), HFSign)
                endif
            endif

            ! Print the determinants with the largest +- populations
            if (opts_selected(printhighpopdet)) then
                tPrintHighPop = .true.
                write(stdout,*) 'Request to print the determinants with the &
                           &largest populations detected on iteration: ',iter
            endif

            ! Change the reference determinant on the fly
            if (opts_selected(changerefdet)) then
                tCheckHighestPopOnce = .true.
                tCheckHighestPop = .true.
                tChangeProjEDet = .true.
                FracLargerDet = 1.0
                write(stdout,*) 'Changing the reference determinant to the most &
                           &highly weighted determinant on iteration: ',iter
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

            ! Change the maximum nopen truncation level
            if (opts_selected(trunc_nopen)) then
                if (tTruncNOpen) then
                    call MPIBcast (trunc_nop_new, tSource)
                    if (trunc_nop_new < 0 .or. trunc_nop_new > nel) then
                        tTruncNOpen = .false.
                        root_print 'Truncation by number of unpaired &
                                   &electrons disabled.'
                    elseif (trunc_nop_new >= trunc_nopen_max) then
                        trunc_nopen_max = trunc_nop_new
                        root_print 'Truncating space to a maximum of ', &
                                   trunc_nopen_max, ' unpaired electrons per &
                                   &determinant.'
                    else
                        root_print 'Cannot decrease truncation level for &
                                   &truncation by number of unpaired &
                                   &electrons during a run.'
                    endif
                else
                    root_print 'WARNING: Cannot enable truncation by number &
                               &of unpaired electrons during a run.'
                endif
            endif

            ! varyshift according to reference population
            if (opts_selected(refshift)) then
                tShiftonHFPop = .true.
                write(stdout,*) 'Request to change default shift action to REFSHIFT &
                &detected on a node on iteration: ',iter
            endif

            ! Initialise calculation of the stochastic RDM.
            if (opts_selected(calc_rdm)) then
                tChangeVarsRDM = .true.
                call MPIBCast (tChangeVarsRDM, tSource)
                call MPIBCast (RDMExcitLevel, tSource)
                call MPIBCast (IterRDMonFly_new, tSource)
                call MPIBCast (RDMEnergyIter, tSource)

                if (IterRDMonFly_new .le. (Iter - maxval(VaryShiftIter))) then
                    root_print 'Request to initialise the STOCHASTIC &
                               &calculation of the density matrices.'
                    root_print 'However the iteration specified to start &
                               &filling has already been.'
                    root_print 'Beginning to fill RDMs in the next iteration.'
                    IterRDMonFly_value = (Iter - maxval(VaryShiftIter)) + 1

                else
                    root_print 'Initialising the STOCHASTIC calculation of &
                               &the reduced density matrices'
                    IterRDMonFly_value = IterRDMonFly_new
                endif


            endif

            ! Initialise calculation of the explicit RDM.
            if (opts_selected(calc_explic_rdm)) then
                if(tHPHF) then
                    root_print 'Trying to set up calculation of the &
                               &EXPLICIT RDM.'
                    root_print 'But the EXPLICIT method does not work with &
                               &HPHF.'
                    root_print 'Ignoring request.'
                else
                    tChangeVarsRDM = .true.
                    tExplicitAllRDM = .true.
                    call MPIBCast (tChangeVarsRDM, tSource)
                    call MPIBCast (tExplicitAllRDM, tSource)
                    call MPIBCast (RDMExcitLevel, tSource)
                    call MPIBCast (IterRDMonFly_new, tSource)
                    call MPIBCast (RDMEnergyIter, tSource)

                    if (IterRDMonFly_new .le. (Iter - maxval(VaryShiftIter))) then
                        root_print 'Request to initialise the EXPLICIT &
                                   &calculation of the density matrices.'
                        root_print 'However the iteration specified to start &
                                   &filling has already been.'
                        root_print 'Beginning to fill RDMs in the next iteration.'
                        IterRDMonFly_value = (Iter - maxval(VaryShiftIter)) + 1

                    else
                        root_print 'Initialising the EXPLICIT calculation of &
                                   &the reduced density matrices'
                        IterRDMonFly_value = IterRDMonFly_new
                    endif
                endif
            endif

            ! Change the starting iteration for filling the rdm.
            if (opts_selected(fill_rdm_iter)) then
                call MPIBCast (IterRDMonFly_new, tSource)

                if (IterRDMonFly_new .le. (Iter - maxval(VaryShiftIter))) then
                    root_print 'New value of IterRDMonFly is LESS than or EQUAL TO &
                               &the current iteration number'
                    root_print 'The number of iterations after the shift change &
                               &to start filling the RDM has been left at ', IterRDMonFly_value
                elseif(tRDMonFly) then
                    IterRDMonFly_value = IterRDMonFly_new
                    if(tExplicitAllRDM) then
                        if(RDMExcitLevel.eq.3) then
                            root_print 'The 1 and 2 electron reduced density matrices &
                                      &will be EXPLICITLY filled '
                            root_print 'from the following number of iterations after the &
                                      &shift changes ', IterRDMonFly_value
                        else
                            root_print 'The ',RDMExcitLevel,' electron reduced density &
                                      &matrices will be EXPLICITLY filled '
                            root_print 'from the following number of iterations after the &
                                      &shift changes ', IterRDMonFly_value
                        endif
                    else
                        if(RDMExcitLevel.eq.3) then
                            root_print 'The 1 and 2 electron reduced density matrices &
                                      &will be STOCHASTICALLY filled '
                            root_print 'from the following number of iterations after the &
                                      &shift changes ', IterRDMonFly_value
                        else
                            root_print 'The ',RDMExcitLevel,' electron reduced density &
                                      &matrices will be STOCHASTICALLY filled '
                            root_print 'from the following number of iterations after the &
                                      &shift changes ', IterRDMonFly_value
                        endif
                    endif
                else
                    root_print 'Attempt to start filling the reduced density matrices.'
                    root_print 'This cannot be done, because the arrays have not &
                               &been set up.'
                    root_print 'Try CALC(EXPLICIT)RDM RDMExcitLevel RDMIter EnergyIter'
                endif
            endif

            if (opts_selected(diag_one_rdm)) then
                tDiagRDM = .true.
                root_print 'Requesting to diagonalise the 1-RDM at the end of the &
                           &calculation.'
            endif

            if (opts_selected(time)) then
                call MPIBcast(MaxTimeExit, tSource)
                if (MaxTimeExit <= 0) then
                    tTimeExit = .false.
                    root_print "Automatic time-based soft-exit disabled."
                else
                    tTimeExit = .true.
                    root_print "Automatic soft-exit set to ", MaxTimeExit, &
                               "mins"
                    MaxTimeExit = MaxTimeExit * 60.0_dp
                end if
            end if

        endif

    end subroutine ChangeVars

logical function test_SoftExit()

    logical :: tdummy1, tdummy2

    tSoftExitFound = .false.
    call ChangeVars(tdummy1, tdummy2)
    if (tSoftExitFound) write (6,'(1X,a30)') 'Request for SOFTEXIT detected.'
    test_SoftExit = tSoftExitFound
    tSoftExitFound = .false.

end function test_SoftExit

end module soft_exit
