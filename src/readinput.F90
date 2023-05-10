! ReadInput is called to read the input.
!   Filename is a Character(255) string which contains a filename to read.
!               If Filename=="" then we check to see if there's a filename on the command line.
!               Failing that, we use stdin
MODULE ReadInput_neci
    use constants, only: stdout, stdin
    use SystemData, only: tUHF, t_fci_pchb_excitgen, tStoreSpinOrbs, tMolpro, &
                         tROHF
    use pchb_excitgen, only: FCI_PCHB_options
    use gasci_pchb_main, only: GAS_PCHB_options
    use gasci_pchb_doubles_main, only: possible_PCHB_hole_selection
    use util_mod, only: operator(.implies.), stop_all
    Use Determinants, only: tDefineDet, DefDet
    use SystemData, only: lms, user_input_m_s, t_k_space_hubbard, t_trans_corr_2body
    use input_parser_mod, only: TokenIterator_t, FileReader_t, ManagingFileReader_t, AttachedFileReader_t
    use fortran_strings, only: to_upper, to_lower, to_int, to_realdp
    use tau_main, only: tau_start_val, possible_tau_start, &
        min_tau, max_tau, tau, readpops_but_tau_not_from_popsfile, MaxWalkerBloom, &
        tau_search_method, possible_tau_search_methods
    use CalcData, only: tTruncInitiator, InitiatorWalkNo, max_allowed_spawn, tScaleBlooms

    Implicit none
!   Used to specify which default set of inputs to use
!    An enum would be nice, but is sadly not supported
    integer, parameter :: idDefault = 0
    integer, parameter :: idFeb08 = 1
    integer, parameter :: idNov11 = 2

contains

    Subroutine ReadInputMain(cFilename, tOverride_input, kp)
        use SystemData, only: tMolpro
        use System, only: SysReadInput, SetSysDefaults
        use Calc, only: CalcReadInput, SetCalcDefaults
        use CalcData, only: tKP_FCIQMC, tUseProcsAsNodes
        use kp_fciqmc_data_mod, only: kp_fciqmc_data
        use kp_fciqmc_init, only: kp_fciqmc_read_inp
        use Integrals_neci, only: IntReadInput, SetIntDefaults
        Use Logging, only: LogReadInput, SetLogDefaults
        use Parallel_neci, only: iProcIndex
        use default_sets
        use util_mod, only: get_free_unit
        use real_time_read_input_module, only: real_time_read_input
!#ifdef NAGF95
!    !  USe doesn't get picked up by the make scripts
!        USe f90_unix_env, ONLY: getarg,iargc
!#endif
        Implicit none
!#ifndef NAGF95
!        Integer :: iargc
!#endif
        !  INPUT/OUTPUT params
        character(*), parameter :: this_routine = 'ReadInputMain'
        Character(len=*) cFilename    !Input  filename or "" if we check arg list or stdin

        Character(len=255) cInp         !temp storage for command line params
        Character(len=32) cTitle
!  Predeclared
        Character(len=100) w, x         !strings for input storage
        Logical tEof        !set when read_line runs out of lines
        logical tExists     !test for existence of input file.
        Integer idDef       !What default set do we use
        integer neci_iargc
        logical, intent(in) :: tOverride_input  !If running through molpro, is this an override input?
        integer, allocatable :: tmparr(:)
        type(kp_fciqmc_data), intent(inout) :: kp
        integer, parameter :: id_scratch_file = 7
        class(FileReader_t), allocatable :: file_reader
        type(TokenIterator_t) :: tokens

        cTitle = ""
        idDef = idDefault                 !use the Default defaults (pre feb08)
        If (trim(adjustl(cFilename)) /= '') Then
            file_reader = ManagingFileReader_t(trim(adjustl(cFilename)))
        else if (neci_iArgC() > 0) then
            ! We have some arguments we can process instead
#ifdef BLUEGENE_HACKS
            call neci_GetArg(neci_iArgC(), cInp)
#else
            Call neci_GetArg(1, cInp)      !Read argument 1 into inp
#endif
            write(stdout, *) "Processing arguments", cinp
            write(stdout, *) "Reading from file: ", Trim(cInp)
            file_reader = ManagingFileReader_t(trim(adjustl(cInp)))
        Else
            write(stdout, *) "Reading from STDIN"
            ! Save the input to a temporary file so we can scan for the
            ! defaults option and then re-read it for all other options.
            open(id_scratch_file, status='scratch')
            file_reader = AttachedFileReader_t(file_id=stdin, echo_lines=id_scratch_file)
        end if

        Do while (file_reader%nextline(tokens, skip_empty=.true.))
            w = to_upper(tokens%next())
            Select case (w)
            Case ("DEFAULTS")
                x = to_upper(tokens%next())
                select case (x)
!Add default options here
                case ("DEFAULT")
                    idDef = idDefault
                case ("FEB08")
                    idDef = idFeb08
                case ("NOV11")
                    idDef = idNov11
                case default
                    write(stdout, *) "No defaults selected - using 'default' defaults"
                    idDef = idDefault
                end select
            case ("END")
                exit
            end select
        End Do

        select case (idDef)
        case (0)
            write(stdout, *) 'Using the default set of defaults.'
        case (idFeb08)
            Feb08 = .true.
            write(stdout, *) 'Using the Feb08 set of defaults.'
        case (idNov11)
            Nov11 = .true.
            write(stdout, *) 'Using the November 2011 set of defaults'
        end select

        ! Set up defaults.
        call SetSysDefaults
        call SetCalcDefaults
        call SetIntDefaults
        call SetLogDefaults

        ! Now return to the beginning and process the whole input file
        select type(file_reader)
        type is (AttachedFileReader_t)
            file_reader = AttachedFileReader_t(file_id=id_scratch_file)
        end select

        call file_reader%rewind()

!Molpro writes out its own input file
        if (.not. tMolpro .or. tOverride_input) then
            if (iProcIndex == 0) then
                call file_reader%set_echo_lines(stdout)
            else
                call file_reader%set_echo_lines()
            end if
            write(stdout, '(/,64("*"),/)')
        end if

        Do while (file_reader%nextline(tokens, skip_empty=.true.))
            w = to_upper(tokens%next())
            select case (w)
            case ("TITLE")
                do while (tokens%remaining_items() > 0)
                    w = tokens%next()
                    cTitle = trim(cTitle)//" "//trim(w)
                end do
            case ("DEFAULTS")
                CONTINUE
            case ("SYSTEM")
                call SysReadInput(file_reader, tokens)
            case ("CALC")
                call CalcReadInput(file_reader)
            case ("INTEGRAL")
                call IntReadInput(file_reader)
            case ("LOGGING")
                call LogReadInput(file_reader)
            case ("KP-FCIQMC")
                tKP_FCIQMC = .true.
                tUseProcsAsNodes = .true.
                call kp_fciqmc_read_inp(file_reader, kp)
            case ("REALTIME")
                call real_time_read_input(file_reader)
            case ("END")
                exit
            case default
                call stop_all(this_routine, "Keyword "//trim(w)//" not recognized")
            end select
        end do
        write(stdout, '(/,64("*"),/)')
        call file_reader%close()
        call sanitize_input()
        RETURN
    END SUBROUTINE ReadInputMain

    subroutine sanitize_input()
        call evaluate_depending_keywords()
        call checkinput()
    end subroutine

    !> @brief
    !>   Certain keywords are optional and/or depend on others.
    !>   Evaluate this dependency here.
    subroutine evaluate_depending_keywords()
        use SystemData, only: tGAS
        use gasci, only: GAS_specification, GAS_exc_gen, &
            possible_GAS_exc_gen, user_input_GAS_exc_gen
        use CalcData, only: user_input_seed, G_VMC_SEED
        character(*), parameter :: this_routine = 'evaluate_depending_keywords'

        if (tGAS .and. allocated(user_input_GAS_exc_gen)) then
            GAS_exc_gen = user_input_GAS_exc_gen
            ! set fast weighting in case of indeterminate setting
            if (GAS_exc_gen == possible_GAS_exc_gen%PCHB) then
                if (GAS_PCHB_options%doubles%hole_selection &
                    == possible_PCHB_hole_selection%INDETERMINATE_FAST_FAST) then
                    if (tUHF) then
                        GAS_PCHB_options%doubles%hole_selection = possible_PCHB_hole_selection%SPINORB_FAST_FAST
                    else
                        GAS_PCHB_options%doubles%hole_selection = possible_PCHB_hole_selection%SPATORB_FAST_FAST
                    end if
                end if ! gasci pchb
            end if ! indeterminate fast-weighted
        end if

        if (allocated(user_input_seed)) then
            G_VMC_SEED = user_input_seed
        end if

        if (allocated(user_input_m_s)) then
            lms = user_input_m_s
        else if (tDefinedet) then
            lms = sum(merge(1, -1, mod(DefDet, 2) == 0))
        else
            lms = 0
        end if
        ! move from `src/readint.F90::INITFROMFCID` which overwrote results here
        ! note this comes before setting tstorespinorbs based on the hole selection
        ! algorithm below, since even if these requirements are not satisfied,
        ! we still want tstorespinorbs = .true.

        tStoreSpinOrbs = (tMolpro .and. tUHF) .or. (tUHF .and. (.not. tROHF))
        ! set fci pchb hole selection in case of indeterminate setting
        if (t_fci_pchb_excitgen) then
            if (FCI_PCHB_options%doubles%hole_selection &
                    == possible_PCHB_hole_selection%INDETERMINATE_FAST_FAST) then
                if (tUHF) then
                    FCI_PCHB_options%doubles%hole_selection = possible_PCHB_hole_selection%SPINORB_FAST_FAST
                else
                    FCI_PCHB_options%doubles%hole_selection = possible_PCHB_hole_selection%SPATORB_FAST_FAST
                end if
            end if
        end if

    end subroutine

    subroutine checkinput()

        ! Check that the specified runtime options are consistent and valid

        use SystemData, only: nel, tUseBrillouin, beta, tFixLz, &
                              tFindCINatOrbs, tNoRenormRandExcits, LMS, STOT, &
                              tSpn, tUHF, tGenHelWeighted, tHPHF, &
                              tGen_4ind_weighted, tGen_4ind_reverse, &
                              tMultiReplicas, tGen_4ind_part_exact, &
                              tGUGA, tgen_guga_weighted, &
                              tGen_4ind_lin_exact, tGen_4ind_2, tGAS, &
                              tComplexOrbs_RealInts, tLatticeGens
        use CalcData, only: I_VMAX, NPATHS, G_VMC_EXCITWEIGHT, &
                            G_VMC_EXCITWEIGHTS, EXCITFUNCS, TMCDIRECTSUM, &
                            TDIAGNODES, TSTARSTARS, TBiasing, TMoveDets, &
                            TNoSameExcit, TInitStar, tMP2Standalone, &
                            MemoryFacPart, &
                            tSemiStochastic, semistoch_shift_iter, ss_space_in, &
                            tSpatialOnlyHash, InitWalkers, tUniqueHFNode, &
                            tCheckHighestPop, &
                            tKP_FCIQMC, tReplicaEstimates, &
                            tRealCoeffByExcitLevel, &
                            tAllRealCoeff, tUseRealCoeffs, tChangeProjEDet, &
                            tOrthogonaliseReplicas, tReadPops, tStartMP1, &
                            tStartCAS, tUniqueHFNode, tContTimeFCIMC, &
                            tContTimeFull, tFCIMC, tPreCond, tOrthogonaliseReplicas, &
                            tMultipleInitialStates, pgen_unit_test_spec, &
                            user_input_seed, tTargetShiftdamp, tFixedN0, t_core_inits, &
                            tWalkContGrow
        use real_time_data, only: tInfInit
        use Calc, only : RDMsamplingiters_in_inp
        Use Determinants, only: SpecDet, tagSpecDet, tDefinedet, DefDet
        use IntegralsData, only: nFrozen, tDiscoNodes, tQuadValMax, &
                                 tQuadVecMax, tCalcExcitStar, tJustQuads, &
                                 tNoDoubs
        use IntegralsData, only: tDiagStarStars, tExcitStarsRootChange, &
                                 tRmRootExcitStarsRootChange, tLinRootChange
        use LoggingData, only: iLogging, tCalcFCIMCPsi, tRDMOnFly, &
                               tCalcInstantS2, tDiagAllSpaceEver, &
                               tCalcVariationalEnergy, tCalcInstantS2Init, &
                               tPopsFile, tRDMOnFly, tExplicitAllRDM, &
                               tHDF5PopsRead, tHDF5PopsWrite, tCalcFcimcPsi, &
                               tHistEnergies, tPrintOrbOcc, tUserKnowsBiasedRDMS
        use Logging, only: calcrdmonfly_in_inp, RDMlinspace_in_inp
        use real_time_data, only: t_real_time_fciqmc
        use DetCalc, only: tEnergy, tCalcHMat, tFindDets, tCompressDets
        use load_balance_calcnodes, only: tLoadBalanceBlocks
        use constants
        use global_utilities
        use FciMCData, only: nWalkerHashes, HashLengthFrac, InputDiagSft, t_global_core_space
        use hist_data, only: tHistSpawn
        use Parallel_neci, only: nNodes, nProcessors
        use UMatCache, only: tDeferred_Umat2d
        use gasci, only: GAS_specification, GAS_exc_gen, possible_GAS_exc_gen, user_input_GAS_exc_gen

        use guga_init, only: checkInputGUGA
        implicit none

        integer :: vv, kk, cc, ierr
        real(dp) :: InputDiagSftSingle
        logical :: check
        character(*), parameter :: t_r = 'checkinput', this_routine = t_r

        if (tDiagAllSpaceEver .and. .not. tHistSpawn) then
            call stop_all(t_r, "DIAGALLSPACEEVER requires HISTSPAWN option")
        end if
        if (tCalcVariationalEnergy .and. .not. tHistSpawn) then
            call stop_all(t_r, "CALCVARIATIONALENERGY requires HISTSPAWN option")
        end if
        if (tCalcVariationalEnergy .and. .not. tEnergy) then
            call stop_all(t_r, "CALCVARIATIONALENERGY requires initial FCI calculation")
        end if

        nWalkerHashes = nint(HashLengthFrac * InitWalkers)

        ! ================ GUGA implementation ===============================
        ! design convention to store as many guga related functionality in
        ! guga_*.F90 files and just call the routines in the main level modules
        ! checkInputGUGA() is found in guga_init.F90
        if (tGUGA) call checkInputGUGA()

        ! Turn on histogramming of fcimc wavefunction in order to find density
        ! matrix, or the orbital occupations
        if (tFindCINatOrbs) tCalcFCIMCPsi = .true.

        ! Used in the FCIMC. We find dets and compress them for later use
        if (tCalcFCIMCPsi .or. tHistSpawn) then
            tFindDets = .true.
            tCompressDets = .true.
        end if

        ! We need to have found the dets before calculating the H mat.
        if (tCalcHMat) tFindDets = .true.

        ! If we are using TNoSameExcit, then we have to start with the star -
        ! the other random graph algorithm cannot remove same excitation
        ! links yet.
        if (tNoSameExcit .and. .not. tInitStar) then
            call stop_all(this_routine, "If we are using TNoSameExcit, then we have to start&
                         & with the star. The other random graph algorithm &
                         &cannot remove same excitation links yet.")
        end if

        ! The MoveDets and Biasing algorithms cannot both be used in the
        ! GraphMorph Algorithm.
        if (tBiasing .and. tMoveDets) then
            call stop_all(this_routine, "Biasing algorithm and MoveDets algorithm cannot both&
                        & be used")
        end if

        ! ..RmRootExcitStarsRootChange must be used with DiagStarStars, and not
        ! with ExcitStarsRootChange
        if (tRmRootExcitStarsRootChange .and. .not. tDiagStarStars) then
            call stop_all(this_routine, "RmRootExcitStarsRootChange can only with used with &
                        &DiagStarStars currently")
        end if

        if (TRmRootExcitStarsRootChange .and. TExcitStarsRootChange) then
            call stop_all(this_routine, "RmRootExcitStarsRootChange and ExcitStarsRootChange &
                        &cannot both be used as they are both different &
                        &options with diagstarstars")
        end if

        !..ExcitStarsRootChange must be used with TDiagStarStars
        if (tExcitStarsRootChange .and. .not. tDiagStarStars) then
            call stop_all(this_routine, "ExcitStarsRootChange can only with used with &
                        &DiagStarStars currently")
        end if

        ! ..TDiagStarStars must be used with TStarStars, and cannot be used
        ! with TCalcExcitStar
        if (tDiagStarStars .and. .not. tStarStars) then
            call stop_all(this_routine, "DiagStarStars must be used with StarStars")
        end if
        if (tDiagStarStars .and. tCalcExcitStar) then
            call stop_all(this_routine, "DiagStarStars is incompatable with CalcExcitStar")
        end if
        if (tDiagStarStars .and. (tNoDoubs .or. tJustQuads)) then
            call stop_all(this_routine, "NoDoubs/JustQuads cannot be used with DiagStarStars &
                        &- try CalcExcitStar")
        end if

        ! ..TNoDoubs is only an option which applied to TCalcExcitStar, and
        ! cannot occurs with TJustQuads.
        if (tNoDoubs .and. .not. tCalcExcitStar) then
            call stop_all(this_routine, "STARNODOUBS is only an option which applied to &
                        &TCalcExcitStar")
        end if

        if (tNoDoubs .and. tJustQuads) then
            call stop_all(this_routine, "STARNODOUBS and STARQUADEXCITS cannot be applied &
                        &together!")
        end if

        ! .. TJustQuads is only an option which applies to TCalcExcitStar
        if (tJustQuads .and. .not. tCalcExcitStar) then
            call stop_all(this_routine, "STARQUADEXCITS is only an option which applies to &
                        &tCalcExcitStar")
        end if

        !.. tCalcExcitStar can only be used with tStarStars
        if (tCalcExcitStar .and. .not. tStarStars) then
            call stop_all(this_routine, "CalcExcitStar can only be used with StarStars set")
        end if

        !.. Brillouin Theorem must be applied when using TStarStars
        if (tStarStars .and. .not. tUseBrillouin) then
            call stop_all(this_routine, "Brillouin Theorem must be used when using &
                        &CalcExcitStar")
        end if

        !.. TQuadValMax and TQuadVecMax can only be used if TLINESTARSTARS set
        if ((tQuadValMax .or. tQuadVecMax) .and. .not. tStarStars) then
            call stop_all(this_routine, "TQuadValMax or TQuadVecMax can only be specified if &
                        &STARSTARS specified in method line")
        end if

        !.. TQuadValMax and TQuadVecMax cannot both be set
        if (tQuadValMax .and. tQuadVecMax) then
            call stop_all(this_routine, "TQuadValMax and TQuadVecMax cannot both be set")
        end if

        !.. TDISCONODES can only be set if NODAL is set in the star methods
        ! section
        if (tDiscoNodes .and. .not. tDiagNodes) then
            call stop_all(this_routine, "DISCONNECTED NODES ONLY POSSIBLE IF NODAL SET IN &
                        &METHOD")
        end if

        if (tMultipleInitialStates .or. tOrthogonaliseReplicas .or. &
            tPreCond) then
            if (tHistSpawn .or. &
                (tCalcFCIMCPsi .and. tFCIMC) .or. tHistEnergies .or. tPrintOrbOcc) then
                call stop_all(this_routine, "HistSpawn and PrintOrbOcc not yet supported for multi-replica with different references")
            end if
        end if

        if (.not. t_global_core_space) then
            if (t_real_time_fciqmc) call stop_all(this_routine, "Real-time FCIQMC requires a global core space")
            if (tKP_FCIQMC) call stop_all(this_routine, "KP-FCIQMC requires a global core space")
            if (tReplicaEstimates) call stop_all(this_routine, "Replica estimates require a global core space")
        end if

        !.. We still need a specdet space even if we don't have a specdet.
        if (.not. associated(SPECDET)) then
            allocate(SPECDET(nel - nFrozen), stat=ierr)
            call LogMemAlloc('SPECDET', nel - nFrozen, 4, t_r, tagSPECDET, ierr)
        end if

        !..   Testing ILOGGING
        !     ILOGGING = 0771
        if (I_VMAX == 0 .and. nPaths /= 0 .and. (.not. tKP_FCIQMC)) then
            call stop_all(this_routine, 'NPATHS!=0 and I_VMAX=0.  VERTEX SUM max level not &
                         &set')
        end if

        !Ensure beta is set.
        if (beta < 1.0e-6_dp .and. .not. tMP2Standalone) then
            call stop_all(this_routine, "No beta value provided.")
        end if

        do vv = 2, I_VMAX
            g_VMC_ExcitWeights(:, vv) = g_VMC_ExcitWeights(:, 1)
            G_VMC_EXCITWEIGHT(vv) = G_VMC_EXCITWEIGHT(1)
        end do

        !IF THERE IS NO WEIGHTING FUNCTION, ExcitFuncs(10)=.true.
        do vv = 1, 9
            IF (EXCITFUNCS(vv)) EXCITFUNCS(10) = .false.
        end do

        if (tNoRenormRandExcits .and. (.not. ExcitFuncs(10))) then
            write(stdout, *) "Random excitations WILL have to be renormalised, &
                       &since an excitation weighting has been detected."
        end if

        ! if the LMS value specified is not reachable with the number of electrons,
        ! fix this
        if (mod(abs(lms), 2) /= mod(nel, 2)) then
            call stop_all(t_r, "LMS Value is not reachable with the given number of electrons.")
        end if

        if (tCalcInstantS2 .or. tCalcInstantS2Init) then
            if (tUHF) then
                call stop_all(t_r, 'Cannot calculate instantaneous values of&
                              & S^2 with UF enabled.')
            end if
            write(stdout, *) 'Enabling calculation of instantaneous S^2 each &
                       &iteration.'
        end if

        if (tUniqueHFNode .and. nProcessors < 2) then
            write(stdout, *) "nNodes: ", nNodes
            write(stdout, *) 'nProcessors: ', nProcessors
            call stop_all(t_r, 'At least two nodes required to designate &
                          &a node uniquely to the HF determinant')
        end if

        if (tGenHelWeighted) then
            write(stdout, *)
            write(stdout, *) '*** WARNING ***'
            write(stdout, *) 'Slow HElement biased excitation generators in use.'
            write(stdout, *) 'NOT FOR PRODUCTION RUNS'
            write(stdout, *) '***************'
            write(stdout, *)
        end if

        if (tGen_4ind_weighted .or. tGen_4ind_reverse .or. tGen_4ind_2 &
            .or. tgen_guga_weighted) then

            ! We want to use UMAT2D...
            tDeferred_Umat2d = .true.

        end if

        if (tHPHF .and. tUHF) then
            call stop_all(t_r, 'HPHF functions cannot work with UHF')
        end if

#if PROG_NUMRUNS_
        if (tKP_FCIQMC .and. .not. tMultiReplicas) then

            write(stdout, *) 'Using KPFCIQMC without explicitly specifying the &
                       &number of replica simulations'
            write(stdout, *) 'Defaulting to using 2 replicas'
            tMultiReplicas = .true.
#ifdef CMPLX_
            lenof_sign = 4
#else
            lenof_sign = 2
#endif
            inum_runs = 2

            ! Correct the size of InputDiagSft:
            InputDiagSftSingle = InputDiagSft(1)
            deallocate(InputDiagSft)
            allocate(InputDiagSft(inum_runs))
            InputDiagSft = InputDiagSftSingle
        end if
#endif

#if PROG_NUMRUNS_
        if (tRDMonFly) then
            write(stdout, *) 'RDM on fly'

            if (.not. tMultiReplicas) then
                write(stdout, *) 'unspecified'
                write(stdout, *) 'Filling RDMs without explicitly specifying the &
                           &number of replica simplations'
                write(stdout, *) 'Defaulting to using 2 replicas'
                tMultiReplicas = .true.
                lenof_sign = 2
                inum_runs = 2

                ! Correct the size of InputDiagSft:
                InputDiagSftSingle = InputDiagSft(1)
                deallocate(InputDiagSft)
                allocate(InputDiagSft(inum_runs))
                InputDiagSft = InputDiagSftSingle
            end if
        end if
#endif


#if ! (defined(PROG_NUMRUNS_) || defined(DOUBLERUN_))
        if (tRDMonFly .and. .not. tUserKnowsBiasedRDMS) then
            write(stdout, *) 'RDM sampling is specified, but this version of neci'
            write(stdout, *) 'is not compiled with the replica trick.'
            write(stdout, *) 'You probably want to use dneci or mneci &
                    &(or their complex counterparts).'
            write(stdout, *)
            write(stdout, *) 'If you know what you do and really want to sample biased RDMS'
            write(stdout, *) 'you can also add `BIASED-RDMS` to the Logging block.'
            call stop_all(t_r, "Compiled version does not support RDM sampling.")
        end if
#endif

        if (tRDMOnFly .and. .not. tCheckHighestPop) then
            write(stdout, *) 'Highest population checking required for calculating &
                       &RDMs on the fly'
            write(stdout, *) 'If you are seeing this, it is an input parsing error'
            call stop_all(t_r, 'RDMs without CheckHighestPop')
        end if

        if (tSemiStochastic .and. .not. (tAllRealCoeff .and. tUseRealCoeffs)) then
            write(stdout, *) 'Semi-stochastic simulations only supported when using &
                       &ALLREALCOEFF option'
            call stop_all(t_r, 'Semistochastic without ALLREALCOEFF')
        end if

        if (ss_space_in%tPopsProportion &
                .and. .not. tSemiStochastic &
                .and. semistoch_shift_iter < 1) then
            call stop_all(t_r, 'POPS-CORE-PROPORTION requires SEMI-STOCHASTIC option')
        end if

        if (ss_space_in%tPopsCore .and. ss_space_in%tPopsProportion) then
            call stop_all(t_r, 'POPS-CORE and POPS-CORE-PROPORTION cannot be used&
                               & at the same time')
        end if

        if (tAllRealCoeff .and. tRealCoeffByExcitLevel) then
            call stop_all(t_r, 'Options ALLREALCOEFF and REALCOEFFBYEXCITLEVEL&
                               & are incompatibile')
        end if

        if (RDMlinspace_in_inp .and. (RDMsamplingiters_in_inp .or. calcrdmonfly_in_inp)) then
            call stop_all(t_r, 'RDMlinspace and (RDMsamplingiters + calcrdmonfly) &
                               &are mutually exclusive')
        end if

        if (tOrthogonaliseReplicas) then
            if (.not. tMultiReplicas) then
                call stop_all(t_r, 'Replica orthogonalisation requires &
                                   &SYSTEM-REPLICAS to determine the number &
                                   &of simulations')
            end if

            if (tStartMP1 .or. tStartCAS) then
                call stop_all(t_r, "MP1 or CAS starting not implemented for &
                                   &orthogonalised calculations")
            end if
        end if

        if (tLoadBalanceBlocks) then
            if (tUniqueHFNode) then
                call stop_all(t_r, "UNIQUE-HF-NODE requires disabling &
                                   &LOAD-BALANCE-BLOCKS")
            end if

            ! If there is only one node, then load balancing doesn't make
            ! a great deal of sense, and only slows things down...
            if (nNodes == 1) then
                write(stdout, *) 'Disabling load balancing for single node calculation'
                tLoadBalanceBlocks = .false.
            end if

            if (tContTimeFCIMC .and. tContTimeFull) then
                call stop_all(t_r, 'Load balancing not yet usable for &
                             &calculations requiring accumulated determinant &
                             &specific global data')
            end if
        end if

#ifndef USE_HDF_
        if (tHDF5PopsRead .or. tHDF5PopsWrite) then
            call stop_all(t_r, 'Support for HDF5 files disabled at compile time')
        end if
#endif

        if (tFixLz .and. tComplexOrbs_RealInts) then
            write(stdout, *) 'Options LZTOT and COMPLEXORBS_REALINTS incompatible'
            write(stdout, *)
            write(stdout, *) '1. Using multiple options that filter integrals at runtime is unsupported.'
            write(stdout, *) '   Only one integral filter may be used at once.'
            write(stdout, *)
            write(stdout, *) '2. This is almost certainly not what you intended to do.  LZTOT works using'
            write(stdout, *) '   abelian symmetries combined with momentum information. COMPLEXORBS_REALINTS'
            write(stdout, *) '   provides support for non-abelian symmetries in FCIDUMP files produced'
            write(stdout, *) '   using VASP'
            write(stdout, *)
            call stop_all(t_r, 'Options incompatible')
        end if

        if (tLatticeGens) then
            if (tGen_4ind_2 .or. tGen_4ind_weighted .or. tGen_4ind_reverse) then
                call stop_all(t_r, "Invalid excitation options")
            end if
        end if

        if (tGAS .neqv. allocated(user_input_GAS_exc_gen)) then
            call stop_all(this_routine, 'GAS-CI and GAS-SPEC required.')
        end if

        if (tGAS) then
            if (.not. tDefineDet) then
                call stop_all(t_r, "Running GAS requires a user-defined reference via definedet.")
            endif
            if (.not. GAS_specification%contains_conf(DefDet)) then
                call stop_all(t_r, "Reference determinant has to be contained in GAS space.")
            endif
            if (.not. GAS_specification%is_valid()) then
                call stop_all(t_r, "GAS specification not valid.")
            end if
            if (.not. GAS_specification%recoupling() .and. all(GAS_exc_gen /= [possible_GAS_exc_gen%DISCONNECTED, possible_GAS_exc_gen%PCHB])) then
                call stop_all(t_r, "Running GAS without spin-recoupling requires {DISCONNECTED, GENERAL_PCHB} implementations.")
            end if
            if (GAS_exc_gen == possible_GAS_exc_gen%DISCONNECTED .and.  GAS_specification%is_connected()) then
                call stop_all(t_r, "Running GAS-CI = DISCONNECTED requires disconnected spaces.")
            end if
        end if

        if (tDefineDet .and. allocated(user_input_m_s)) then
            if (sum(merge(1, -1, mod(DefDet, 2) == 0)) /= user_input_m_s) then
                call stop_all(t_r, "Spin of Definedet and Spin-restrict is not consistent.")
            end if
        end if

        if (allocated(pgen_unit_test_spec) .and. .not. tReadPops) then
            call stop_all(t_r, "UNIT-TEST-PGEN requires READPOPS.")
        end if

        if (tTargetShiftdamp .and. tFixedN0) then
            call stop_all(t_r, "TARGET-SHIFTDAMP and FIXED-N0 not compatible.")
        end if

        if (tTargetShiftdamp .and. tWalkContGrow) then
            call stop_all(t_r, "TARGET-SHIFTDAMP and WALKCONTGROW not compatible.")
        end if

        if (.not. (tInfInit .implies. t_core_inits)) then
            call stop_all(t_r, 'INFINITE-INIT requires CORE-INITS OFF.')
        end if

        block
            use load_balance_calcnodes, only: &
                tLoadBalanceBlocks, loadBalanceInterval
            logical :: deterministic
            if (allocated(user_input_seed)) then
                deterministic = loadBalanceInterval /= 0
                if (tLoadBalanceBlocks .and. .not. deterministic) then
        write(stdout, *) 'Seed was specified in input.'
        write(stdout, *) 'Please note that because of load-balancing the calculation is not fully deterministic (CPU load).'
        write(stdout, *) 'If a fully deterministic calculation is required use the `load-balance-interval` keyword.'
                end if
            end if
        end block


        time_step: block
            if (.not. allocated(tau_start_val)) then
                call stop_all(t_r, 'Start value for tau is required.')
            end if

            if (tau_start_val == possible_tau_start%from_popsfile .neqv. tReadPops) then
                if (tau_start_val == possible_tau_start%from_popsfile) then
                    call stop_all(t_r, 'Starting tau from popsfile requires readpops.')
                else if (.not. readpops_but_tau_not_from_popsfile) then
                    write(stdout, *) 'Using readpops while not reading tau from popsfile is &
                        &very likely an input error.'
                    write(stdout, *) 'If you think that your input is correct and you know what you do &
                        &add `readpops-but-tau-not-from-popsfile` to the tau-values keywords'
                    call stop_all(t_r, 'Readpops requires `tau-values start from-popsfile`.')
                end if
            end if

            if (tau_start_val == possible_tau_start%tau_factor &
                    .and. t_trans_corr_2body .and. t_k_space_hubbard) then
                call stop_all(this_routine, &
                              "finding the number of excits from HF breaks for too large lattice.")
            end if

            if (tau_start_val == possible_tau_start%refdet_connections .and. tGUGA) then
                call stop_all(this_routine, &
                              "tau-values start refdet-connections is not compatible with GUGA calculations!")
            end if

            if (tau_search_method /= possible_tau_search_methods%off) then
                if (tTruncInitiator .and. MaxWalkerBloom > InitiatorWalkNo) then
                    call stop_all(this_routine, &
                                  "MaxWalkerBloom has to be smaller equal than InitiatorWalkNo.")
                end if
                if (tScaleBlooms .and. MaxWalkerBloom > max_allowed_spawn) then
                    call stop_all(this_routine, &
                                  "MaxWalkerBloom has to be smaller equal than max_allowed_spawn.")
                end if
            end if
        end block time_step

    end subroutine checkinput

end Module ReadInput_neci
