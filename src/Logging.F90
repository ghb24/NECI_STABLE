#include "macros.h"

MODULE Logging

    use constants, only: dp, int64, nreplicas
    use MemoryManager, only: LogMemAlloc, LogMemDealloc, TagIntType
    use SystemData, only: nel, LMS, nbasis, tGUGA
    use CalcData, only: tCheckHighestPop, semistoch_shift_iter, trial_shift_iter, &
                        tPairedReplicas, tReplicaEstimates, iSampleRDMIters, tMoveGlobalDetData
    use constants, only: n_int, size_n_int, bits_n_int
    use bit_rep_data, only: NIfTot, NIfD
    use DetBitOps, only: EncodeBitDet
    use hist_data, only: iNoBins, tHistSpawn, BinRange
    use errors, only: Errordebug
    use LoggingData
    use spectral_data, only: tPrint_sl_eigenvecs

! RT_M_Merge: There seems to be no conflict here, so use both
    use real_time_data, only: n_real_time_copies, t_prepare_real_time, &
                              cnt_real_time_copies

    use rdm_data, only: nrdms_transition_input, states_for_transition_rdm, tApplyLC
    use rdm_data, only: rdm_main_size_fac, rdm_spawn_size_fac, rdm_recv_size_fac

    use analyse_wf_symmetry, only: t_symmetry_analysis, t_symmetry_mirror, &
                                   t_symmetry_rotation, symmetry_rotation_angle, &
                                   t_symmetry_inversion, symmertry_mirror_axis, &
                                   t_read_symmetry_states, n_symmetry_states, &
                                   t_pop_symmetry_states, symmetry_states, &
                                   symmetry_weights, symmetry_states_ilut

    use cc_amplitudes, only: t_plot_cc_amplitudes

    use fcimcdata, only: tFillingStochRDMonFly

    use input_parser_mod, only: FileReader_t, TokenIterator_t

    use fortran_strings, only: to_upper, to_lower, to_int, to_realdp

    IMPLICIT NONE

    logical, public :: RDMlinspace_in_inp, calcrdmonfly_in_inp

contains

    subroutine SetLogDefaults()
        != Set defaults for Logging data items.

        use default_sets
        implicit none

        ! real-time implementation changes:
        n_real_time_copies = 1
        cnt_real_time_copies = 1
        t_prepare_real_time = .false.

        ! By default, the output is given by the shift cycle
        StepsPrint = 10
        tCoupleCycleOutput = .true.

        tDipoles = .false.
        tPrintInitiators = .false.
        tDiagAllSpaceEver = .false.
        tCalcVariationalEnergy = .false.
        tJustBlocking = .false.
        iBlockEquilShift = 0
        iBlockEquilProjE = 0
        ErrorDebug = 0
        iHighPopWrite = 15    !How many highest weighted determinants to write out at the end of an FCIQMC calc.
        t_force_replica_output = .false.
        tDiagWalkerSubspace = .false.
        iDiagSubspaceIter = 1
        PopsfileTimer = 0.0_dp
        tMCOutput = .true.
        tLogComplexPops = .false.
        iWriteBlockingEvery = 1000
        tSaveBlocking = .false.
        OffDiagBinRange = 0.001_dp
        OffDiagMax = 1.0_dp
        BinRange = 0.001_dp
        iNoBins = 100000
        tHistEnergies = .false.
        tHistSpawn = .false.
        iWriteHistEvery = -1
        NoACDets(:) = 0
        TAutoCorr = .false.
        MaxHistE = 50.0_dp
        NoHistBins = 200
        iWritePopsEvery = 100000
        TCalcWavevector = .false.
        WavevectorPrint = 100
        TPopsFile = .true.
        tIncrementPops = .false.
        tPrintPopsDefault = .true.
        TDistrib = .false.
        ILOGGINGDef = 0
        iGlobalTimerLevel = 40
        nPrintTimer = 10
        HFLOGLEVEL = 0
        PreVarLogging = 0
        TDetPops = .false.
        TZeroProjE = .false.
        TWriteDetE = .false.
        iPopsPartEvery = 1
        tBinPops = .false.
        tROHistogramAll = .false.
        tROFciDump = .true.
        tTruncRODump = .false.
        tTruncDumpbyVal = .false.
        tROHistER = .false.
        tROHistDoubExc = .false.
        tROHistOffDiag = .false.
        tROHistSingExc = .false.
        tROHistOnePartOrbEn = .false.
        tROHistOneElInts = .false.
        tPrintInts = .false.
        tPrintSpinCoupHEl = .false.
        tPrintFCIMCPsi = .false.
        tCalcFCIMCPsi = .false.
        NHistEquilSteps = 0
        tPrintOrbOcc = .false.
        StartPrintOrbOcc = 0
        tPrintOrbOccInit = .false.
        FCIMCDebug = 0
        tHFPopStartBlock = .false.
        tIterStartBlock = .false.
        IterStartBlocking = 0
        HFPopStartBlocking = 100
        tInitShiftBlocking = .false.
        NoDumpTruncs = 0
        tWriteTransMat = .false.
        tHistInitPops = .false.
        HistInitPopsIter = 100000
        tLogDets = .false.
        tLogEXLEVELStats = .false.
        tCalcInstantS2 = .false.
        tCalcInstantS2Init = .false.
        tCalcInstSCpts = .false.
        tCalcPropEst = .false.
        iNumPropToEst = 0
        instant_s2_multiplier = 1
        tRDMonFly = .false.
        tFillingStochRDMonFly = .false.
        tChangeVarsRDM = .false.
        RDMEnergyIter = 100
        tDiagRDM = .false.
        tPrint1RDM = .false.
        tNoNOTransform = .false.
        tPrintRODump = .false.
        IterRDMonFly = 0
        RDMExcitLevel = 1
        tDo_Not_Calc_2RDM_est = .false.
        tExplicitAllRDM = .false.
        twrite_normalised_RDMs = .true.
        tWriteSpinFreeRDM = .false.
        twrite_RDMs_to_read = .false.
        tno_RDMs_to_read = .false.
        tReadRDMs = .false.
        tNoNewRDMContrib = .false.
        IterWriteRDMs = 10000
        tWriteMultRDMs = .false.
        tThreshOccRDMDiag = .false.
        ThreshOccRDM = 2.0_dp
        tDumpForcesInfo = .false.
        tPrintLagrangian = .false.
        instant_s2_multiplier_init = 1
        binarypops_min_weight = 0
        tSplitPops = .false.
        tWriteCore = .false.
        tWriteCoreEnd = .false.
        write_end_core_size = 0
        tWriteTrial = .false.
        tCompareTrialAmps = .false.
        compare_amps_period = 0
        tHistExcitToFrom = .false.
        tForceCauchySchwarz = .false.
        tBrokenSymNOs = .false.
        occ_numb_diff = 0.001_dp
        tBreakSymNOs = .false.
        local_cutoff = 0
        rottwo = 0
        rotthree = 0
        rotfour = 0
        tRDMInstEnergy = .true.
        tFullHFAv = .false.
        tPrintDataTables = .true.
        tOutputLoadDistribution = .false.
        tHDF5PopsRead = .false.
        tHDF5PopsWrite = .false.
        tReduceHDF5Pops = .false.
        HDF5PopsMin = 1.0_dp
        iHDF5PopsMinEx = 4
        tPopsProjE = .false.
        tHDF5TruncPopsWrite = .false.
        iHDF5TruncPopsEx = 0
        iHDF5TruncPopsIter = 0
        tAccumPops = .false.
        tAccumPopsActive = .false.
        iAccumPopsIter = 0
        iAccumPopsMaxEx = 2
        iAccumPopsExpireIters = 0
        AccumPopsExpirePercent = 0.9_dp
        iAccumPopsCounter = 0
        tWriteRefs = .false.
        maxInitExLvlWrite = 8
#ifdef PROG_NUMRUNS_
        tFCIMCStats2 = .true.
#else
        tFCIMCStats2 = .false.
#endif
        tFvalEnergyHist = .false.
        FvalEnergyHist_EnergyBins = 100
        FvalEnergyHist_FValBins = 10
        tFvalPopHist = .false.
        FvalPopHist_PopBins = 100
        FvalPopHist_FValBins = 10

! Feb08 defaults
        IF (Feb08) THEN
            !Mcpaths set
            ILOGGINGDef = 2
        end if

        ref_filename = "REFERENCES"

    end subroutine SetLogDefaults

    subroutine LogReadInput(file_reader)

        ! Read the logging section from the input file

        logical :: eof
        logical tUseOnlySingleReplicas
        integer :: i, line, ierr
        integer :: n_samples
        character(100) :: w
        character(100) :: PertFile(3)
        class(FileReader_t), intent(inout) :: file_reader
        character(*), parameter :: t_r = 'LogReadInput'
        character(*), parameter :: this_routine = 'LogReadInput'
        type(TokenIterator_t) :: tokens

        tUseOnlySingleReplicas = .false.

        ILogging = iLoggingDef

        PertFile(:) = ''
        logging: do while (file_reader%nextline(tokens, skip_empty=.true.))
            w = to_upper(tokens%next())
            select case (w)

            case('CI-COEFFICIENTS')
                ! collects ci coefficients of the wave function up to 3rd excitation level over a number of iter
                t_store_ci_coeff = .true.
                if (item < nitems) then
                   call readi(n_iter_ci_coeff)
                end if
                if (item < nitems) then
                   call readi(n_store_ci_level)
                end if

            case ("PRINT-MOLCAS-RDMS")
                ! output density matrices also in Molcas format in the GUGA RDM
                ! implementation
                t_print_molcas_rdms = .true.

            case ("PRINT-FREQUENCY-HISTOGRAMS")
                ! in this case print the frequency histograms to analyze the
                ! matrix element vs. pgen ratios
                t_print_frq_histograms = .true.

            case ("REBLOCKSHIFT")
                !Abort all other calculations, and just block data again with given equilibration time (in iterations)
                tJustBlocking = .true.
                iBlockEquilShift = to_int(tokens%next())
            case ("REBLOCKPROJE")
                !Abort all other calculations, and just block data again with given equilibration time (in iterations)
                tJustBlocking = .true.
                iBlockEquilProjE = to_int(tokens%next())
            case ("HIGHLYPOPWRITE")
                !At the end of an FCIMC calculation, how many highly populated determinants should we write out?
                iHighPopWrite = to_int(tokens%next())
            case("REPLICAS-POPWRITE")
                ! Print out the highest populated determinants from all replicas
                t_force_replica_output = .true.
                if( tokens%remaining_items() > 0) iHighPopWrite = to_int(tokens%next())
            case ("DIAGWALKERSUBSPACE")
                !Diagonalise walker subspaces every iDiagSubspaceIter iterations
                tDiagWalkerSubspace = .true.
                iDiagSubspaceIter = to_int(tokens%next())
            case ("DIAGALLSPACEEVER")
                !Diagonalise all space ever visited in the fciqmc dynamic. This will be written out each time HistSpawn is
                tDiagAllSpaceEver = .true.
            case ("CALCVARIATIONALENERGY")
                !Calculate the variational energy of the FCIQMC dynamic each time Histspawn is calculated
                tCalcVariationalEnergy = .true.

            case ("SPLITPROJE", "SPLITPROJE-G", "SPLITPROJE-K3")
                call stop_all(t_r, 'Option (SPLITPROJE*) deprecated')

            case ("NOMCOUTPUT")
                !No output to stdout from the fcimc iterations
                tMCOutput = .false.

            case ("STEPSOUTPUT")
                ! This is the number of steps taken between two lines in the output
                ! The default is equal to the update cycle length of the shift, since
                ! this saves some communication
                ! This clearly indicates that we do not want to have output and shift update
                ! going hand in hand
                tCoupleCycleOutput = .false.
                StepsPrint = to_int(tokens%next())

            case ("LOGCOMPLEXWALKERS")
                !This means that the complex walker populations are now logged.
                tLogComplexPops = .true.

            case ("PRINTNEWBLOCKING")
!This is the iteration interval period to write out the blocking files.
                iWriteBlockingEvery = to_int(tokens%next())
            case ("SAVEBLOCKING")
!In this case, blocking files are not overwritten each time they are printed out, but
                tSaveBlocking = .true.
            case ("ERRORBLOCKING")
!Performs blocking analysis on the errors in the instantaneous projected energy to get the error involved.
!This is default on, but can be turned off with this keyword followed by OFF.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tHFPopStartBlock = .false.
                    end select
                ELSE
                    tHFPopStartBlock = .true.
                end if

            case ("SHIFTERRORBLOCKING")
!Performs blocking analysis on the errors in the instantaneous projected energy to get the error involved.
!This is default on, but can be turned off with this keyword followed by OFF.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tInitShiftBlocking = .false.
                    end select
                ELSE
                    tInitShiftBlocking = .true.
                end if

            case ("BLOCKINGSTARTITER")
!This keyword can be used if we want to start the blocking error analysis at a particular iteration.
!If it is a negative integer, then this means that the blocking will start when we come out of fixed shift mode.
                tIterStartBlock = .true.
                tHFPopStartBlock = .false.
                IterStartBlocking = to_int(tokens%next())

            case ("SHIFTBLOCKINGSTARTITER")
!This keyword can be used if we want to start the blocking error analysis of the shift at a particular
!iteration after the shift begins to change.
                call stop_all(t_r, 'SHIFTBLOCKINGSTARTITER option deprecated')

            case ("BLOCKINGSTARTHFPOP")
!This keyword can be used if we want to start the blocking error analysis at a particular HF population.
!The current default is 100.
                tHFPopStartBlock = .true.
                HFPopStartBlocking = to_int(tokens%next())

            case ("ROFCIDUMP")
!Turning this option on prints out a new FCIDUMP file at the end of the orbital rotation.  At the moment, the rotation is very slow
!so this will prevent us having to do the transformation every time we run a calculation on a particular system
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROFCIDUMP = .false.
                    end select
                ELSE
                    tROFCIDUMP = .true.
                end if

            case ("TRUNCROFCIDUMP")
!This options truncates the rotated FCIDUMP file by removing the specified number of virtual orbitals, based on the occupation
!numbers given by diagonalisation of the MP2 variational density matrix.
                tTruncRODump = .true.
                NoDumpTruncs = 1
                allocate(NoTruncOrbs(NoDumpTruncs), stat=ierr)
                CALL LogMemAlloc('NoTruncOrbs', NoDumpTruncs, 4, 'Logging', NoTruncOrbsTag, ierr)
                NoTruncOrbs(:) = 0
                do i = 1, NoDumpTruncs
                    NoTruncOrbs(i) = to_int(tokens%next())
                end do

            case ("MULTTRUNCROFCIDUMP")
!This option allows us to specify multiple truncations, so that one calculation will print out multiple
!ROFCIDUMP files with different
!levels of truncation - prevents us from having to do multiple identical CISD calculations to get the different truncations.
                tTruncRODump = .true.
                NoDumpTruncs = to_int(tokens%next())
                allocate(NoTruncOrbs(NoDumpTruncs), stat=ierr)
                CALL LogMemAlloc('NoTruncOrbs', NoDumpTruncs, 4, 'Logging', NoTruncOrbsTag, ierr)
                NoTruncOrbs(:) = 0
                do i = 1, NoDumpTruncs
                    NoTruncOrbs(i) = to_int(tokens%next())
                end do

            case ("MULTTRUNCVALROFCIDUMP")
!This option allows us to specify particular cutoffs values for the eigenvalues - and print out multiply
!ROFCIDUMP files with orbitals
!with eigenvalues below these removed.
                tTruncRODump = .true.
                tTruncDumpbyVal = .true.
                NoDumpTruncs = to_int(tokens%next())
                allocate(TruncEvalues(NoDumpTruncs), stat=ierr)
                CALL LogMemAlloc('TruncEvalues', NoDumpTruncs, 8, 'Logging', TruncEvaluesTag, ierr)
                TruncEvalues(:) = 0.0_dp
                do i = 1, NoDumpTruncs
                    TruncEvalues(i) = to_realdp(tokens%next())
                end do

            case ("WRITETRANSFORMMAT")
!This option writes out the transformation matrix used to convert the HF orbitals into the natural orbitals.
!This can then be read into
!QChem to produce the natural orbital cube files and then visualise them using VMD.  Note : Currently,
!because of Fortran 90's weird dealings
!with writing and reading binary - this option is only compatible with QChem if the code is compiled using
!PGI - this will be fixed at
!some stage.  Also - QChem INTDUMP files must be used to be compatible.
                tWriteTransMat = .true.

            case ("HIST-INTEGRALS")
                tHistLMat = .true.

            case ("ROHISTOGRAMALL")
!This option goes with the orbital rotation routine.  If this keyword is included, all possible histograms are included.
!These maybe also turned off/on with individual keywords.
!As it stands, the bins run from -1 to 1 with increments of 0.05. These parameters may be made options in the future.
                tROHistogramAll = .true.
                tROHistOffDiag = .true.
                tROHistDoubExc = .true.
                tROHistSingExc = .true.
                tROHistOneElInts = .true.
                tROHistOnePartOrbEn = .true.
                tROHistER = .true.
                tROHistVirtCoulomb = .true.

            case ("ROHISTOFFDIAG")
!This option creates a histogram of the <ij|kl> terms where i<k and j<l.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    tROHistOffDiag = (w /= 'OFF')
                    select case (w)
                    case ("OFF")
                        tROHistOffDiag = .false.
                    end select
                ELSE
                    tROHistOffDiag = .true.
                end if

            case ("ROHISTDOUBEXC")
!This option creates a histogram of the 2<ij|kl>-<ij|lk> terms, the off diagonal hamiltonian elements for double excitations.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROHistDoubExc = .false.
                    end select
                ELSE
                    tROHistDoubExc = .true.
                end if

            case ("ROHISTSINGEXC")
!This option creates a histogram of the single excitation hamiltonian elements.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROHistSingExc = .false.
                    end select
                ELSE
                    tROHistSingExc = .true.
                end if

            case ("ROHISTER")
!This option creates a histogram of the <ii|ii> terms, the ones that are maximised in the edmiston-reudenberg localisation.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROHistER = .false.
                    end select
                ELSE
                    tROHistER = .true.
                end if

            case ("ROHISTONEElINTS")
!This option creates a histogram of the one electron integrals, the <i|h|i> terms.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROHistOneElInts = .false.
                    end select
                ELSE
                    tROHistOneElInts = .true.
                end if

            case ("ROHISTONEPARTORBEN")
!This option creates a histogram of the one particle orbital energies, epsilon_i = <i|h|i> + sum_j [<ij||ij>].
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROHistOnePartOrbEn = .false.
                    end select
                ELSE
                    tROHistOnePartOrbEn = .true.
                end if

            case ("ROHISTVIRTCOULOMB")
!This option creates a histogram of the coulomb integrals, <ij|ij>, where i and j are both virtual and i<j.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tROHistVirtCoulomb = .false.
                    end select
                ELSE
                    tROHistVirtCoulomb = .true.
                end if

            case ("PRINTINTEGRALS")
!This option prints 2 files containing the values of certain integrals at each rotation iteration.  This is so that we can see the
!effect the rotation is having on all values, other than just the one we are max/minimising.
                tPrintInts = .true.

            case ("PRINTTRICONNECTIONS")
!This option takes each generated pair of determinant and excitation and finds 3rd determinant to make up a triangular connection.
!The product of the three connecting elements are then histogrammed in two separate files.
!In one, the triangular connections that combine
!to be sign coherent are recorded, and in the other, those which are sign incoherent.
                CALL Stop_All(t_r, "PRINTTRICONNECTIONS option depreciated")
!            TriConMax = to_realdp(tokens%next())
!            NoTriConBins = to_int(tokens%next())
!            tPrintTriConnections=.true.

            case ("HISTTRICONNELEMENTS")
!This keyword takes the above triangles of connected determinants and histograms each connecting
!element that contributes to the triangle.
!It then prints these according to whether they are single or double connecting elements.
!It also prints a histogram and the average size of the Hjk elements (regardless of whether or not they are zero).
                CALL Stop_All(t_r, "HISTTRICONNELEMENTS option depreciated")
!            TriConHElSingMax = to_realdp(tokens%next())
!            TriConHElDoubMax = to_realdp(tokens%next())
!            NoTriConHElBins = to_int(tokens%next())
!            tHistTriConHEls=.true.

            case ("PRINTHELACCEPTSTATS")
!This keyword prints out an extra file that keeps track of the H elements involved in spawning attempts
!that are accepted or not accepted.
!It prints out the average H elements where spawning is accepted and the average where it is not accepted.
                CALL Stop_All(t_r, "PRINTHELACCEPTSTATS option depreciated")
!            tPrintHElAccept=.true.

            case ("PRINTSPINCOUPHELS")
!This option prints out the number of positive and negative (and their sums) H elements connecting two spin coupled determinants.
                tPrintSpinCoupHEl = .true.

            case ("HISTINITIATORPOPS")
!This option prints out a file (at every HistInitPopsIter iteration) containing the
!natural log of the populations of the initiator determinants
!and the number with this population. The range of populations histogrammed goes from ln(N_add) -> ln(1,000,000) with 50,000 bins.
                tHistInitPops = .true.
                HistInitPopsIter = to_int(tokens%next())

#if defined(PROG_NUMRUNS_)
            case ("PAIRED-REPLICAS")
                tPairedReplicas = .true.
                nreplicas = 2
#endif

            case ("UNPAIRED-REPLICAS")
                tUseOnlySingleReplicas = .true.
#if defined(PROG_NUMRUNS_)
                tPairedReplicas = .false.
                nreplicas = 1
#elif defined(DOUBLERUN_)
                call stop_all(t_r, "The unpaired-replicas option cannot be used with the dneci.x executable.")
#endif

#if defined(PROG_NUMRUNS_)
            case ("REPLICA-ESTIMATES")
                tReplicaEstimates = .true.
                tPairedReplicas = .true.
                nreplicas = 2
#endif

            case ("FULL-CORE-RDMS")
                ! samples the rdms within the core space between ALL
                ! states and not only the one connected by H
                t_full_core_rdms = .true.

            case ("CALCRDMONFLY")
!This keyword sets the calculation to calculate the reduced density matrix on the fly.
!This starts at IterRDMonFly iterations after the shift changes.
!If RDMExcitLevel = 1, only the 1 electron RDM is found, if RDMExcitLevel = 2,
! only the 2 electron RDM is found and if RDMExcitLevel = 3, both are found.
                calcrdmonfly_in_inp = .true.
                tRDMonFly = .true.
                tCheckHighestPop = .true.
                RDMExcitLevel = to_int(tokens%next())
                IterRDMonFly = to_int(tokens%next())
                RDMEnergyIter = to_int(tokens%next())

#if defined(PROG_NUMRUNS_)
                ! With this option, we want to use pairs of replicas.
                if (.not. tUseOnlySingleReplicas) then
                    tPairedReplicas = .true.
                    nreplicas = 2
                end if
#elif defined(DOUBLERUN_)
                tPairedReplicas = .true.
#endif
                if (IterRDMOnFly < semistoch_shift_iter) then
                    call stop_all(t_r, "Semi-stochastic needs to be &
                        &turned on before RDMs are turned on.")
                else if (IterRDMOnFly < trial_shift_iter) then
                    call stop_all(t_r, "Trial wavefunctions needs to be &
                        &turned on before RDMs are turned on.")
                end if

            case ("RDMLINSPACE")
!> This keyword is an improved way to specify the RDM sampling intervals.
!> The syntax is ``RDMlinspace  start n_samples  step``.
!> The RDMExcitLevel is set to three in this routine.
                RDMlinspace_in_inp = .true.
                tRDMonFly = .true.
                tCheckHighestPop = .true.

                RDMExcitLevel = 3
                IterRDMonFly = to_int(tokens%next())
                n_samples = to_int(tokens%next())
                RDMEnergyIter = to_int(tokens%next())

                iSampleRDMIters = n_samples * RDMEnergyIter
#if defined(PROG_NUMRUNS_)
                ! With this option, we want to use pairs of replicas.
                if (.not. tUseOnlySingleReplicas) then
                    tPairedReplicas = .true.
                    nreplicas = 2
                end if
#elif defined(DOUBLERUN_)
                tPairedReplicas = .true.
#endif
                if (IterRDMOnFly < semistoch_shift_iter) then
                    call stop_all(t_r, "Semi-stochastic needs to be &
                        &turned on before RDMs are turned on.")
                else if (IterRDMOnFly < trial_shift_iter) then
                    call stop_all(t_r, "Trial wavefunctions needs to be &
                        &turned on before RDMs are turned on.")
                end if

            case ("BIASED-RDMS")
                ! Only relevant for (k)-neci runs.
                ! By default the calculation stops with an error if RDMs are sampled with (k)-neci
                ! to prevent user error.
                ! With this keyword the user can explicitly say that they want to sample RDMs without replica.
                tUserKnowsBiasedRDMS = .true.

            case ("OLDRDMS")
                call stop_all(t_r, "OLDRDMS not supported anymore.")

            case ("RDM-MAIN-SIZE-FAC")
                rdm_main_size_fac = to_realdp(tokens%next())
            case ("RDM-SPAWN-SIZE-FAC")
                rdm_spawn_size_fac = to_realdp(tokens%next())
            case ("RDM-RECV-SIZE-FAC")
                rdm_recv_size_fac = to_realdp(tokens%next())

            case ("TRANSITION-RDMS")
                tTransitionRDMs = .true.
                nrdms_transition_input = to_int(tokens%next())
                allocate(states_for_transition_rdm(2, nrdms_transition_input), stat=ierr)

                do line = 1, nrdms_transition_input
                    if (file_reader%nextline(tokens, skip_empty=.false.)) then
                        do i = 1, 2
                            states_for_transition_rdm(i, line) = to_int(tokens%next())
                        end do
                    else
                        call stop_all(t_r, 'Unexpected end of file reached.')
                    end if
                end do

            case ("PRINT-1RDMS-FROM-2RDM-POPS")
                tPrint1RDMsFrom2RDMPops = .true.
                tReadRDMs = .true.

            case ("PRINT-1RDMS-FROM-SPINFREE")
                tPrint1RDMsFromSpinfree = .true.

            case ("NO-APPEND-STATS")
                t_no_append_stats = .true.

            case ("DIAGFLYONERDM")
!This sets the calculation to diagonalise the *1* electron reduced density matrix.
!The eigenvalues give the occupation numbers of the natural orbitals (eigenfunctions).
                tDiagRDM = .true.

            case ("FULLHFAV")
                !Continue to accumulate the average of N_HF even when it goes to zero
                !Necessary for good RDM accumulation in systems with small N_HF (i.e. multireference)
                tFullHFAv = .true.

            case ("NONOTRANSFORM")
! This tells the calc that we don't want to print the NO_TRANSFORM matrix.
! i.e. the diagonalisation is just done to get the correlation entropy.
                tNoNOTransform = .true.

            case ("PRINTONERDM")
! This prints the OneRDM, regardless of whether or not we are calculating just the 1-RDM, or the 2-RDM.
                tPrint1RDM = .true.

            case ("PRINTRODUMP")
                tPrintRODump = .true.
                tROFciDump = .true.
! This is to do with the calculation of the MP2 or CI natural orbitals.
!This should be used if we want the transformation matrix of the
! natural orbitals to be found, but no ROFCIDUMP file to be printed (i.e.
!the integrals don't need to be transformed).  This is so that at the end
! of a calculation, we may get the one body reduced density matrix from the
!wavefunction we've found, and then use the MOTRANSFORM file printed to
! visualise the natural orbitals with large occupation numbers.

            case ("FORCECAUCHYSCHWARZ")
                tForceCauchySchwarz = .true.
!This forces the inequality gamma_pq <= sqrt(gamma_pp * gamma_qq) is obeyed.
!we choose min(gamma_pq, sqrt(gamma_pp * gamma_qq) to ensure a positive-definite matrix
!May be of use for getting orbitals from an approximate initial FCIQMC calc.

            case ("BROKENSYMNOS")
                tBrokenSymNOs = .true.
                rottwo = to_int(tokens%next())
                rotthree = to_int(tokens%next())
                rotfour = to_int(tokens%next())
                local_cutoff = to_int(tokens%next())
                occ_numb_diff = to_realdp(tokens%next())
! This is to rotate the obtained natural orbitals (NOs) again in order to obtain
! symmetry broken NOs: pairs of NOs whose occupation numbers differ by less
! than the specified threshold occ_numb_diff (relative difference, i.e. difference
! divided by absolute value) will be rotated so as to
! maximally localise them using and Edminston Ruedenberg type localisation
! local_cutoff is the index of the spatial orbital which is the borderline
! for performing localisation or delocalisation of the NO: all chosen NOs pairs with
! orbital index less or equal to (and hence occupation numbers larger than)
! this orbital will be delocalised while the others will be localised
! If BREAKSYMNOS is specified rottwo gives number of pairs to rotate, rotthree the
! number of triples to rotate and rotfour the number of quadruples to rotate
! If BREAKSYMNOS is not present these will be ignored
! A new FCIDUMP file (BSFCIDUMP) with the rotated NOs is printed out

            case ("BREAKSYMNOS")
                tBreakSymNOs = .true.
! This is another option for BROKENSYMNOS p1 p2... t1 t2 t3... q1 q2 q3 q4...
! This contains just an ordered list of the spatial orbital indices of the NOs
! to rotate, firstly the doubles, then triples, then quadruples (the number of
! pairs, triples and quadruples is specified with the BROKENSYMNOS options), e.g.
! BROKENSYMNOS 2 1 1 0 0.1
! BREAKSYMNOS 3 4 4 5 1 5 6 7 8 9 10
! will rotate (2,3) (4,5) (1,5,6) (7,8,9,10)
                allocate(RotNOs((2 * rottwo) + (3 * rotthree) + (4 * rotfour)), stat=ierr)
                tagRotNOs = 0
                call LogMemAlloc('RotNOs', ((2 * rottwo) + (3 * rotthree) + (4 * rotfour))&
                    &, 4, t_r, tagRotNOs, ierr)
                RotNOs(:) = 0
                do i = 1, ((2 * rottwo) + (3 * rotthree) + (4 * rotfour))
                    RotNOs(i) = to_int(tokens%next())
                end do

            case ("DIPOLE_MOMENTS")
                !Calculate the dipole moments if we are in molpro
                tDipoles = .true.

            case ("CALCRDMENERGY")
                call stop_all(t_r, "The CALCRDMENERGY option has been replaced by CALC-2RDM-ESTIMATES. &
                                   &The 2-RDM energy is calculated by default when 2-RDMs are being &
                                   &sampled, so this option is only needed if one wants to turn this off.")

            case ("CALC-2RDM-ESTIMATES")
!This takes the 1 and 2 electron RDM and calculates the energy using the RDM expression.
!For this to be calculated, RDMExcitLevel must be = 3, so there is a check to make sure this
!is so if the CALCRDMENERGY keyword is present.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tDo_Not_Calc_2RDM_est = .true.
                    end select
                ELSE
                    tDo_Not_Calc_2RDM_est = .false.
                end if

            case ("CALC-PROP-ESTIMATES")
!Calculate the estimates of the one-electron properties using 1 electron RDM and 1 electron
!property integrals. It uses all the different RDMs that have been estimated and get the
!corresponding property estimations.
                tCalcPropEst = .true.
                if (tokens%remaining_items() == 0) then
                    call stop_all(t_r, "Please specify the name of the integral file corresponding the property")
                end if
                ! iNumPropToEst is the total number of properties to be estimated
                iNumPropToEst = iNumPropToEst + 1
                if (iNumPropToEst > 3) then
                    call stop_all(t_r, 'Only 3 different property integrals allowed')
                end if
                PertFile(iNumPropToEst) = to_upper(tokens%next())

            case ("NORDMINSTENERGY")
!Only calculate and print out the RDM energy (from the 2-RDM) at the end of the simulation
!This saves memory by only having to store one set of RDMs on the headnode rather than two
                tRDMInstEnergy = .false.

            case ("EXPLICITALLRDM")
!Explicitly calculates all the elements of the RDM.
                tExplicitAllRDM = .true.

            case ("WRITEINITIATORS")
                ! Requires a popsfile to be written out.  Writes out the initiator
                ! populations.
                tPrintInitiators = .true.

            case ("WRITERDMSTOREAD")
                ! Writes out the unnormalised RDMs (in binary), so they can be read
                ! back in, and the calculations restarted at a later point This is
                ! also tied to the POPSFILE/BINARYPOPS keyword - so if we're
                ! writing a normal POPSFILE, we'll write this too, unless
                ! **WRITERDMSTOREAD** OFF is used.
                IF (tokens%remaining_items() > 0) THEN
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("OFF")
                        tno_RDMs_to_read = .true.
                        twrite_RDMs_to_read = .false.
                    end select
                ELSE
                    twrite_RDMs_to_read = .true.
                    tno_RDMs_to_read = .false.
                end if

            case ("NONORMRDMS")
                ! Does not print out the normalised (final) RDMs - to be used if
                ! you know the calculation will not be converged, and don't want to
                ! take up disk space.
                twrite_normalised_RDMs = .false.

            case ("WRITE-SPIN-FREE-RDM")
                tWriteSpinFreeRDM = .true.

            case ("READRDMS")
! Read in the RDMs from a previous calculation, and continue accumulating the RDMs from the very beginning of this restart.
                tReadRDMs = .true.

            case ("NONEWRDMCONTRIB")
                ! To be used with READRDMs.  This option makes sure that we don't add in any
                ! new contributions to the RDM if filling stochastically
                ! This is useful if we want to read in an RDM from another calculation and then
                ! just print out the analysis, without adding in any more information.
                tNoNewRDMContrib = .true.

            case ("WRITERDMSEVERY")
! Write out the normalised, hermitian RDMs every IterWriteRDMs iterations.
                tWriteMultRDMs = .true.
                IterWriteRDMs = to_int(tokens%next())

            case ("THRESHOCCONLYRDMDIAG")
                ! Only add in a contribution to the diagonal elements of the RDM if the average sign
                ! of the determinant is greater than [ThreshOccRDM]
                tThreshOccRDMDiag = .true.
                ThreshOccRDM = to_realdp(tokens%next())

            case ("DUMPFORCESINFO")
! Using the finalised 2RDM, calculate the Lagrangian X used for the calculation of the forces,
!and dump all these in Molpro-friendly format
                tDumpForcesInfo = .true.

            case ("PRINTLAGRANGIAN")
                ! Print out the Lagrangian X to file (Only works in conjuction with DUMPFORCESINFO: otherwise, this option does nothing)
                tPrintLagrangian = .true.

            case ("AUTOCORR")
!This is a Parallel FCIMC option - it will calculate the largest weight MP1 determinants and histogramm them
!HF Determinant is always histogrammed. NoACDets(2) is number of doubles. NoACDets(3) is number of triples and NoACDets(4) is
!number of quads to histogram.
                TAutoCorr = .true.
                CALL Stop_All("LogReadInput", "The ACF code has been commented out in the FCIMCPar module")
                do i = 2, 4
                    IF (tokens%remaining_items() > 0) NoACDets(i) = to_int(tokens%next())
                end do
            case ("DETPOPS")
!This option no longer works...
                TDetPops = .true.
            case ("DISTRIBS")
                TDistrib = .true.
            case ("HISTPARTENERGIES")
!This will histogram the hamiltonian matrix elements of the particles in the parallel FCIMC algorithm.
                tHistEnergies = .true.
                BinRange = to_realdp(tokens%next())
                iNoBins = to_int(tokens%next())
                OffDiagBinRange = to_realdp(tokens%next())
                OffDiagMax = to_realdp(tokens%next())
                IF (OffDiagMax < 0.0_dp) THEN
                    OffDiagMax = -OffDiagMax
                end if
            case ("HISTSPAWN")
!This option will histogram the spawned wavevector, averaged over all previous iterations.
!It scales horrifically and can only be done for small systems
!which can be diagonalized. It requires a diagonalization initially to work.
!It can write out the average wavevector every iWriteHistEvery.
                tHistSpawn = .true.
                IF (tokens%remaining_items() > 0) iWriteHistEvery = to_int(tokens%next())
            case ("HISTHAMIL")
!This option will histogram the spawned hamiltonian, averaged over all previous iterations. It scales horrifically
!and can only be done for small systems
!which can be diagonalized. It will write out the hamiltonian every iWriteHamilEvery.
                call stop_all(t_r, 'HISTHAMIL option deprecated')
            case ("BLOCKEVERYITER")
!This will block the projected energy every iteration with the aim of achieving accurate error estimates.
!However, this does require a small amount of additional communication.
                tBlockEveryIteration = .true.
            case ("PRINTFCIMCPSI")
                tPrintFCIMCPsi = .true.
                tCalcFCIMCPsi = .true.
            case ("HISTEQUILSTEPS")
!This option sets the histogramming to only be done after the specified number of iterations.
                NHistEquilSteps = to_int(tokens%next())
            case ("PRINTORBOCCS")
!This option initiates the above histogramming of determinant populations and then at the end of the
!spawning uses these to find the normalised
!contribution of each orbital to the total wavefunction.
                tPrintOrbOcc = .true.
                IF (tokens%remaining_items() > 0) StartPrintOrbOcc = to_int(tokens%next())

            case ("PRINTDOUBSUEG")
                call stop_all(t_r, 'This option (PRINTDOUBSUEG) has been deprecated')

            case ("PRINTORBOCCSINIT")
!This option initiates the above histogramming of determinant populations and then
!at the end of the spawning uses these to find the normalised
!contribution of each orbital to the total wavefunction.
                tPrintOrbOcc = .true.
                tPrintOrbOccInit = .true.
                IF (tokens%remaining_items() > 0) StartPrintOrbOcc = to_int(tokens%next())
            case ("POPSFILE")
! This is so that the determinants at the end of the MC run are written
! out, to enable them to be read back in using READPOPS in the Calc section,
! if you want to restart the simulation at a later date.  !iWritePopsEvery
! will write the configuration of particles out each time the iteration
! passes that many.
                TPopsFile = .true.
                IF (tokens%remaining_items() > 0) THEN
                    iWritePopsEvery = to_int(tokens%next())
                    IF (iWritePopsEvery < 0) THEN
!If a negative argument is supplied to iWritePopsEvery, then the POPSFILE will
!never be written out, even at the end of a simulation.
!If it is exactly zero, this will be the same as without any argument, and a
!popsfile will only be written out in the instance of a clean exit
                        TPopsFile = .false.
                        tPrintPopsDefault = .false.
                    else if (iWritePopsEvery > 0) THEN
                        tPrintPopsDefault = .false.
                    end if
                end if
            case ("REDUCEDPOPSFILE")
!A reduced popsfile works in exactly the same way as a normal popsfile, but only every iPopsPartEvery particle is printed out.
                TPopsFile = .true.
                iWritePopsEvery = to_int(tokens%next())
                iPopsPartEvery = to_int(tokens%next())
            case ("POPSFILETIMER")
                PopsfileTimer = to_realdp(tokens%next())   !Write out a POPSFILE every "PopsfileTimer" hours.

            case ("BINARYPOPS")
                ! This means that the popsfile (full or reduced) will now be
                ! written out in binary format. This should now take up less
                ! space, and be written quicker.
                !
                ! By default, all particles are written into the popsfile. If
                ! a minimum weight is proveded, only those particles with at ]east
                ! that weight are included.
                tBinPops = .true.
                if (tokens%remaining_items() > 0) then
                    binarypops_min_weight = to_realdp(tokens%next())
                end if

            case ("HDF5-POPS")
                ! Use the new HDF5 popsfile format
                tHDF5PopsRead = .true.
                tHDF5PopsWrite = .true.

            case ("REDUCE-HDF5-POPS")

                ! Avoid writing a determinant to HDF5-popsfiles when its population
                ! is below or equal iHDF5PopsMin and its excitation is above iHDF5PopsMinEx
                ! Default values are 1.0 and 4, respectively.

                tReduceHDF5Pops = .true.

                if (tokens%remaining_items() > 0) then
                    HDF5PopsMin = to_realdp(tokens%next())
                    if (HDF5PopsMin < 0.0_dp) then
                        call stop_all(t_r, 'Minimum population should be greater than or equal zero')
                    end if
                end if

                if (tokens%remaining_items() > 0) then
                    iHDF5PopsMinEx = to_int(tokens%next())
                    if (iHDF5PopsMinEx < 2) then
                        call stop_all(t_r, 'Excitation of minimum population should be greater than one')
                    end if
                end if

            case ("HDF5-POPS-READ")
                ! Use the new HDF5 popsfile format just for reading
                tHDF5PopsRead = .true.

            case ("HDF5-POPS-WRITE")
                ! Use the new HDF5 popsfile format just for writing
                tHDF5PopsWrite = .true.

            case ("POPS-PROJE")
                ! Calculate and print the projected energy of
                ! popsfile wavefunction - instantaneous and accumulated (if available)
                tPopsProjE = .true.

            case ("HDF5-TRUNC-POPS-WRITE")
                ! Write another HDF5 popsfile with dets restricted to a maximum
                ! exitation level and/or minimum population
                tHDF5TruncPopsWrite = .true.
                iHDF5TruncPopsEx = to_int(tokens%next())
                if (iHDF5TruncPopsEx < 2) then
                    call stop_all(t_r, 'Maximum excitation level should be greater than 1')
                end if

                ! Number of iterations for the periodic writing of truncated popsfiles.
                ! The default value of zero indicates no periodic writing but
                ! only once at the end.
                if (tokens%remaining_items() > 0) then
                    iHDF5TruncPopsIter = to_int(tokens%next())

                    if (iHDF5TruncPopsEx < 0) then
                        call stop_all(t_r, 'Number of iterations should be greater than or equal zero')
                    end if
                end if

            case ("ACCUM-POPS")
                ! Accumulate the population of determinants and write them
                ! to the popsfile
                tAccumPops = .true.
                ! When to start accumulating the populations
                iAccumPopsIter = to_int(tokens%next())

                ! Normally, when dets become empty, they are removed from CurrentDets
                ! and any associated info (global_det_data) is lost. Therefore,
                ! when accumlating populations is active, (some) empty dets are kept alive.

                if (tokens%remaining_items() > 0) then
                    ! This parameter represents the maximum excitation level to consider
                    ! when keeping empty dets alive.
                    ! We keep up to double excitations indefinitely anyway.
                    ! Therefore, this should be greater than or equal two.
                    ! Default value: 2
                    iAccumPopsMaxEx = to_int(tokens%next())
                    if (iAccumPopsMaxEx < 2) then
                        call stop_all(t_r, 'iAccumPopsMaxEx should be greater than or equal two')
                    end if
                end if

                if (tokens%remaining_items() > 0) then
                    ! This parameter represents the maximum number of iterations,
                    ! empty dets are kept before being removed. The removal happens
                    ! when CurrentDets is almost full (see iAccumPopsExpirePercent).
                    ! A value of zero means keeping accumulated empty dets indefinitely.
                    ! Default value: 0
                    iAccumPopsExpireIters = to_int(tokens%next())
                    if (iAccumPopsExpireIters < 0) then
                        call stop_all(t_r, 'iAccumPopsExpireIters should be greater than or equal zero')
                    end if
                end if

                if (tokens%remaining_items() > 0) then
                    ! This parameter represents how full CurrentDets should be before
                    ! removing accumulated empty dets according to the above criteria.
                    ! Default value: 0.9
                    AccumPopsExpirePercent = to_realdp(tokens%next())
                    if (AccumPopsExpirePercent < 0.0_dp .or. AccumPopsExpirePercent > 1.0) then
                        call stop_all(t_r, 'iAccumPopsExpirePercent should be between zero and one.')
                    end if
                end if

                ! Accumlated populations are stored in global det data, so we need
                ! to preserve them when the dets change processors during load balancing
                tMoveGlobalDetData = .true.

                ! Print popsfile projected energy at the end
                tPopsProjE = .true.

            case ("INCREMENTPOPS")
! Don't overwrite existing POPSFILES.
                tIncrementPops = .true.
            case ("FCIMCDEBUG")
!FCIQMC debugging level. Takes an integer 0-6
                FCIMCDebug = to_int(tokens%next())
            case ("ERRORDEBUG")
!Error analysus debugging level. Takes an integer 0-6
                ErrorDebug = to_int(tokens%next())
            case ("WRITEDETE")
!This logging option will write out the energies of all determinants which have been spawned at in the simulation
! The two input options are the number of bins, and the maximum determinant energy to be histogrammed.
                TWriteDetE = .true.
                IF (tokens%remaining_items() > 0) NoHistBins = to_int(tokens%next())
                IF (tokens%remaining_items() > 0) MaxHistE = to_realdp(tokens%next())
            case ("ZEROPROJE")
! This is for FCIMC when reading in from a POPSFILE. If this is on, then the energy
! estimator will be restarted.
                TZeroProjE = .true.
            case ("WAVEVECTORPRINT")
! This is for FCIMC - if on, it will calculate the exact eigenvector and
! values initially, and then print out the running wavevector every
! WavevectorPrint MC steps. However, this is slower.
                TCalcWavevector = .true.
                WavevectorPrint = to_int(tokens%next())
            case ("BLOCKING")
                ILOGGING = IOR(ILOGGING, 2**13)
            case ("PREVAR")
                ILOGGING = IOR(ILOGGING, 2**14)
            case ("FMCPR")
!  We log the value
                ILOGGING = IOR(ILOGGING, 2**0)
                do while (tokens%remaining_items() > 0)
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("LABEL")
                        ILOGGING = IOR(ILOGGING, 2**2)
                    case ("RHO")
                        ILOGGING = IOR(ILOGGING, 2**3)
                    case ("1000")
                        ILOGGING = IOR(ILOGGING, 2**9)
                    case ("EXCITATION")
                        ILOGGING = IOR(ILOGGING, 2**12)
                    case ("XIJ")
                        ILOGGING = IOR(ILOGGING, 2**6)
                    case ("")
                        ILOGGING = IOR(ILOGGING, 2**2)
                    case default
                        call stop_all(this_routine, "Logging keyword FMCPR "//trim(w)       &
           &               //" not recognised", .true.)
                    end select
                end do
            case ("CALCPATH")
                do while (tokens%remaining_items() > 0)
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("LABEL")
                        ILOGGING = IOR(ILOGGING, 2**4)
                    case ("RHO")
                        ILOGGING = IOR(ILOGGING, 2**5)
                    case ("")
                        ILOGGING = IOR(ILOGGING, 2**4)
                    case default
                        call stop_all(this_routine, "Logging keyword CALCPATH "//trim(w)    &
           &               //" not recognised", .true.)
                    end select
                end do
            case ("XIJ")
                ILOGGING = IOR(ILOGGING, 2**6)
            case ("HAMILTONIAN")
                ILOGGING = IOR(ILOGGING, 2**7)
            case ("PSI")
                ILOGGING = IOR(ILOGGING, 2**8)
            case ("TIMING")
                do while (tokens%remaining_items() > 0)
                    w = to_upper(tokens%next())
                    select case (w)
                    case ("LEVEL")
                        iGlobalTimerLevel = to_int(tokens%next())
                    case ("PRINT")
                        nPrintTimer = to_int(tokens%next())
                    case default
                        call tokens%reset(-1)
                        iGlobalTimerLevel = to_int(tokens%next())
                    end select
                end do
            case ("VERTEX")
                do while (tokens%remaining_items() > 0)
                    w = to_upper(tokens%next())
                    select case (w)
                        ! case("1000")
                        !     ILOGGING = IOR(ILOGGING,2**9)
                    case ("EVERY")
                        ILOGGING = IOR(ILOGGING, 2**10)
                    case default
                        call tokens%reset(-1)
                        G_VMC_LOGCOUNT = to_int(tokens%next())
                        ILOGGING = IOR(ILOGGING, 2**9)
!                  call stop_all(this_routine, "Logging keyword VERTEX "//trim(w)    &
!     &               //" not recognised",.true.)
                    end select
                end do
            case ("HFBASIS")
                ILOGGING = IOR(ILOGGING, 2**11)
            case ("HFLOGLEVEL")
                HFLOGLEVEL = to_int(tokens%next())
            case ("SAVEPREVARLOGGING")
                PreVarLogging = iLogging
                iLogging = iLoggingDef
            case ("DETS")
                tLogDets = .true.
            case ("DETERMINANTS")
                tLogDets = .true.
            case ("EXLEVEL")
                tLogEXLEVELStats = .true.
            case ("INITS-EXLVL-WRITE")
                ! up to which excitation level the number of initiators is written out
                if (tokens%remaining_items() > 0) maxInitExLvlWrite = to_int(tokens%next())
            case ("INSTANT-S2-FULL")
                ! Calculate an instantaneous value for S^2, and output it to the
                ! relevant column in the FCIMCStats file.
                !
                ! The second parameter is a multiplier such that we only calculate
                ! S^2 once for every n update cycles (it must be on an update
                ! cycle such that norm_psi_squared is correct)
                tCalcInstantS2 = .true.
                if (tokens%remaining_items() > 0) instant_s2_multiplier = to_int(tokens%next())

            case ("PLOT-CC-AMPLITUDES")
                t_plot_cc_amplitudes = .true.

            case ("INSTANT-S2-INIT")
                ! Calculate an instantaneous value ofr S^2 considering only the
                ! initiators, and output it to the relevant column in the
                ! FCIMCStats file.
                !
                ! The second parameter is a multiplier such that we only calculate
                ! S^2 once for every n update cycles (it must be an update
                ! cycle such that norm_psi_squared is correct).
                tCalcInstantS2Init = .true.
                if (tokens%remaining_items() > 0) &
                    instant_s2_multiplier_init = to_int(tokens%next())

            case ("INSTANT-S-CPTS")
                ! Calculate components of the wavefunction with each value of S.
                ! n.b. This is NOT quantitatively correct.
                !      --> Only of QUALITATIVE utility.
                tCalcInstSCpts = .true.

            case ("SPLIT-POPS")
                ! Do we want to split a popfile up into multiple parts which are
                ! output on each of the nodes, and need to be combined/split-up and
                ! distributed to the nodes on our (sequential) time rather than
                ! on multi-processor time?
                tSplitPops = .true.
                tBinPops = .true.

            case ("WRITE-CORE")
                ! Output the semi-stochastic core space to a file.
                tWriteCore = .true.

            case ("PRINT-CORE-INFO")
                ! print core info, like energy, and maybe also the gs vector
                t_print_core_info = .true.

            case ("PRINT-CORE-HAMIL")
                t_print_core_info = .true.
                t_print_core_hamil = .true.

            case ("PRINT-CORE-VEC")
                t_print_core_vec = .true.

            case ("WRITE-MOST-POP-CORE-END")
                ! At the end of a calculation, find the write_end_core_size most
                ! populated determinants and write them to a CORESPACE file.
                tWriteCoreEnd = .true.
                write_end_core_size = to_int(tokens%next())

            case ("WRITE-TRIAL")
                ! Output the trial wavefunction space to a file.
                tWriteTrial = .true.

            case ("COMPARE-TRIAL-AND-FCIQMC-AMPS")
                tCompareTrialAmps = .true.
                compare_amps_period = to_int(tokens%next())

            case ("HIST-EXCIT-TOFROM")
                ! Histogram how many particles are spawned between sites with
                ! varying excitation levels from the Hartree--Fock.
                tHistExcitToFrom = .true.

            case ("FVAL-ENERGY-HIST")
                ! When using auto-adaptive shift, print a histogram of the shift factors over
                ! the energy
                tFValEnergyHist = .true.
                if (tokens%remaining_items() > 0) FValEnergyHist_EnergyBins = to_int(tokens%next())
                if (tokens%remaining_items() > 0) FValEnergyHist_FvalBins = to_int(tokens%next())

            case ("FVAL-POP-HIST")
                ! When using auto-adaptive shift, print a histogram of the shift factors over
                ! the population
                tFValPopHist = .true.
                if (tokens%remaining_items() > 0) FValPopHist_PopBins = to_int(tokens%next())
                if (tokens%remaining_items() > 0) FValPopHist_FvalBins = to_int(tokens%next())

            case ("ENDLOG")
                exit logging

            case ("TAU-SEARCH")
                ! Log the output of tau searching
                call stop_all(t_r, 'TAU-SEARCH option now deprecated')

            case ("POPS-MAX-TAU")
                ! If we were using the full enumeration excitation generator,
                ! what would the maximum acceptable value of tau be for the
                ! read-in walker distribution?
                call stop_all(t_r, 'POPS-MAX-TAU Deprecated')

            case ("FCIMCSTATS-2")
                ! Use the new-style FCIMCStats output.
                tFCIMCStats2 = .true.

            case ("PRINT-SL-EIGENVECS")
                tPrint_sl_eigenvecs = .true.

            case ("DONT-PRINT-DATA-TABLES")
                tPrintDataTables = .false.

            case ("LOAD-DISTRIBUTION")
                ! By default we don't output the load balancing distribution of
                ! particles between blocks, as for any reasonable sized system
                ! there are _many_ blocks.
                tOutputLoadDistribution = .true.

            case ("PRINT-UMAT")
                ! output umat also in the momentum space hubbard to be able to
                ! create a FCIDUMP file to compare GUGA matrix elements with
                ! DMRG results!
                t_umat_output = .true.

            case ("DOUBLE-OCCUPANCY")
                ! new functionality to measure the mean double occupancy
                ! as this is a only diagonal quantitity I decided to detach it
                ! from the RDM calculation, although it could be calculated
                ! from the RDMs and this should be used to test this functionality!
                ! Also, as it is a diagonal quantity, we need to unbias the
                ! quantitiy by using the replica trick, just like for the
                ! RDMs! Also this should be tested, to what extend the
                ! quantity differs in a biased and unbiased calculation

                t_calc_double_occ = .true.
                t_calc_double_occ_av = .true.

                if (tokens%remaining_items() > 0) then
                    t_calc_double_occ_av = .false.
                    equi_iter_double_occ = to_int(tokens%next())
                end if

            case ("LOCAL-SPIN")
                t_measure_local_spin = .true.

                if (.not. tGUGA) then
                    call stop_all(this_routine, &
                        "Guga required for local spin measurement!")
                end if

            case ("PRINT-SPIN-RESOLVED-RDMS")
                ! for giovanni enable the output of the spin-resolved rdms not
                ! only for ROHF calculations
                t_spin_resolved_rdms = .true.

            case ("LOG-IJA")
                t_log_ija = .true.

                if (tokens%remaining_items() > 0) then
                    ija_thresh = to_realdp(tokens%next())
                end if

            case ("WRITE-REFERENCES")
                ! Output the reference dets to a file
                tWriteRefs = .true.

            case ("SPIN-MEASUREMENTS")
                ! combine all the spatially resolved double occupancy and
                ! spin-difference measurements into one functionality to
                ! have a better overview
                ! this also includes the "standard" double occupancy measurement
                ! although leave the option to only do the old double occ meas.
                t_calc_double_occ = .true.
                t_calc_double_occ_av = .true.
                t_spin_measurements = .true.

                if (tokens%remaining_items() > 0) then
                    t_calc_double_occ_av = .false.
                    equi_iter_double_occ = to_int(tokens%next())
                end if

            case ('SYMMETRY-ANALYSIS')
                ! investigate the symmetry of the important part of the
                ! wavefuntion, by applying point-group symmetry operations on
                ! a certain number of determinants and check the sign change to
                ! the original wavefunction
                ! if we want multiple symmetries on has to specify this keyword
                ! multiple times with the according keywords
                t_symmetry_analysis = .true.

                if (tokens%remaining_items() > 0) then
                    w = to_lower(tokens%next())

                    select case (w)

                    case ('rot', 'rotation')
                        t_symmetry_rotation = .true.

                        if (tokens%remaining_items() > 0) then
                            symmetry_rotation_angle = to_realdp(tokens%next())
                        else
                            symmetry_rotation_angle = 90.0_dp
                        end if

                    case ('mirror')
                        t_symmetry_mirror = .true.

                        ! available mirror axes are : 'x','y','d' and 'o'
                        if (tokens%remaining_items() > 0) then
                            symmertry_mirror_axis = to_lower(tokens%next())
                        else
                            symmertry_mirror_axis = 'x'
                        end if

                    case ('inversion')
                        t_symmetry_inversion = .true.

                    case default
                        call stop_all(this_routine, "Logging keyword "//trim(w)//" not recognised", .true.)

                    end select
                else
                    ! default is 90° rotation:
                    t_symmetry_rotation = .true.
                    symmetry_rotation_angle = 90.0_dp

                end if

            case ('SYMMETRY-STATES')
                ! required input to determine which states to consider.
                ! Two options:
                if (tokens%remaining_items() > 0) then
                    w = to_lower(tokens%next())
                    select case (w)
                    case ('input', 'read')
                        t_read_symmetry_states = .true.

                        if (tokens%remaining_items() > 0) then
                            n_symmetry_states = to_int(tokens%next())
                        else
                            call stop_all(t_r, &
                                          "symmetry-states input need number of states!")
                        end if

                        allocate(symmetry_states(nel, n_symmetry_states))
                        symmetry_states = 0

                        do line = 1, n_symmetry_states
                            if (file_reader%nextline(tokens, skip_empty=.false.)) then
                                do i = 1, nel
                                    symmetry_states(i, line) = to_int(tokens%next())
                                end do
                            else
                                call stop_all(t_r, 'Unexpected end of file reached.')
                            end if
                        end do

                    case ('pop', 'most-populated', 'pops')
                        ! take the N most populated states
                        t_pop_symmetry_states = .true.

                        if (tokens%remaining_items() > 0) then
                            n_symmetry_states = to_int(tokens%next())
                        else
                            call stop_all(t_r, &
                                          "symmetry-states input need number of states!")
                        end if

                        allocate(symmetry_states(nel, n_symmetry_states))
                        symmetry_states = 0

                    end select
                else
                    ! default is: take the 6 most populated ones
                    t_pop_symmetry_states = .true.
                    n_symmetry_states = 6
                    allocate(symmetry_states(nel, n_symmetry_states))
                    symmetry_states = 0
                end if

                allocate(symmetry_weights(n_symmetry_states))
                symmetry_weights = 0.0_dp

                allocate(symmetry_states_ilut(0:niftot, n_symmetry_states))
                symmetry_states = 0_n_int

            case default
                call stop_all(this_routine, "Logging keyword "//trim(w)//" not recognised", .true.)
            end select
        end do logging

        if (tCalcPropEst) then
            !Save the name of the integral files in the proper place
            if (iNumPropToEst == 0) call stop_all(t_r, 'Error in the property estimations')
            if (allocated(EstPropFile)) deallocate(EstPropFile)
            allocate(EstPropFile(iNumPropToEst))
            EstPropFile(:) = PertFile(1:iNumPropToEst)
        end if
    END SUBROUTINE LogReadInput

END MODULE Logging
