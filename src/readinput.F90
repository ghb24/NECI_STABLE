! ReadInput is called to read the input.
!   Filename is a Character(255) string which contains a filename to read.
!               If Filename=="" then we check to see if there's a filename on the command line.
!               Failing that, we use stdin
!   ios is an Integer which is set to 0 on a successful return, or is non-zero if a file error has occurred, where it is the iostat.
MODULE ReadInput_neci
    Implicit none
!   Used to specify which default set of inputs to use
!    An enum would be nice, but is sadly not supported
    integer, parameter :: idDefault = 0
    integer, parameter :: idFeb08 = 1
    integer, parameter :: idNov11 = 2

    contains

    Subroutine ReadInputMain(cFilename,ios,tOverride_input, kp)
        USE input_neci
        use SystemData, only : tMolpro
        use System,     only : SysReadInput,SetSysDefaults
        use Calc,       only : CalcReadInput,SetCalcDefaults
        use CalcData, only: tKP_FCIQMC, tUseProcsAsNodes
        use kp_fciqmc_data_mod, only: kp_fciqmc_data
        use kp_fciqmc_init, only: kp_fciqmc_read_inp
        use Integrals_neci,  only : IntReadInput,SetIntDefaults
        Use Logging,    only : LogReadInput,SetLogDefaults
        use Parallel_neci,   only : iProcIndex
        use default_sets
        use util_mod, only: get_free_unit
!#ifdef NAGF95
!    !  USe doesn't get picked up by the make scripts
!        USe f90_unix_env, ONLY: getarg,iargc
!#endif
        Implicit none
!#ifndef NAGF95
!        Integer :: iargc
!#endif
    !  INPUT/OUTPUT params
        Character(len=*)  cFilename    !Input  filename or "" if we check arg list or stdin
        Integer             ios         !Output 0 if no error or nonzero iostat if error

        Character(len=255)  cInp         !temp storage for command line params
        Character(len=32)   cTitle
!  Predeclared
!        Integer             ir         !The file descriptor we are reading from
        Character(len=100)  w,x         !strings for input storage
        Logical             tEof        !set when read_line runs out of lines
        logical             tExists     !test for existence of input file.
        Integer             idDef       !What default set do we use
        integer neci_iargc
        logical, intent(in) :: tOverride_input  !If running through molpro, is this an override input?
        integer, allocatable :: tmparr(:)
        type(kp_fciqmc_data), intent(inout) :: kp

        cTitle=""
        idDef=idDefault                 !use the Default defaults (pre feb08)
        ir=get_free_unit()              !default to a free unit which we'll open below
        If(trim(adjustl(cFilename)) .ne. '') Then
            Write(6,*) "Reading from file: ", Trim(cFilename)
            inquire(file=cFilename,exist=tExists)
            if (.not.tExists) call stop_all('ReadInputMain','File '//Trim(cFilename)//' does not exist.')
            Open(ir,File=cFilename,Status='OLD',err=99,iostat=ios)
        ElseIf(neci_iArgC().gt.0) then
    ! We have some arguments we can process instead
#ifdef BLUEGENE_HACKS
            call neci_GetArg(neci_iArgC(), cInp)
#else
            Call neci_GetArg(1,cInp)      !Read argument 1 into inp
#endif
            write(6,*) "Processing arguments", cinp
            Write(6,*) "Reading from file: ", Trim(cInp)
            inquire(file=cInp,exist=tExists)
            if (.not.tExists) call stop_all('ReadInputMain','File '//Trim(cInp)//' does not exist.')
            Open(ir,File=cInp,Status='OLD',FORM="FORMATTED",err=99,iostat=ios)
        Else
            ir=5                    !file descriptor 5 is stdin
            Write(6,*) "Reading from STDIN"
            ! Save the input to a temporary file so we can scan for the
            ! defaults option and then re-read it for all other options.
            open(7,status='scratch',iostat=ios)
        Endif
        Call input_options(echo_lines=.false.,skip_blank_lines=.true.)

    !Look to find default options (line can be added anywhere in input)
        Do
            Call read_line(tEof)
            if (ir.eq.5) write (7,*) trim(char) ! Dump line from STDIN to temporary file.
            If(tEof) Exit
            Call readu(w)
            Select case(w)
            Case("DEFAULTS")
                call readu(x)
                select case(x)
!Add default options here
                case("DEFAULT")
                    idDef=idDefault
                case("FEB08")
                    idDef=idFeb08
                case("NOV11")
                    idDef=idNov11
                case default
                    write(6,*) "No defaults selected - using 'default' defaults"
                    idDef=idDefault
                end select
            case("END")
                exit
            end select
        End Do

        select case(idDef)
        case(0)
            write (6,*) 'Using the default set of defaults.'
        case(idFeb08)
            Feb08 = .true.
            write (6,*) 'Using the Feb08 set of defaults.'
        case(idNov11)
            Nov11 = .true.
            write(6,*) 'Using the November 2011 set of defaults'
        end select

        ! Set up defaults.
        call SetSysDefaults
        call SetCalcDefaults
        call SetIntDefaults
        call SetLogDefaults

!Now return to the beginning and process the whole input file
        if (ir.eq.5) ir=7 ! If read from STDIN, re-read from our temporary scratch file.
        Rewind(ir)
        if(tMolpro.and.(.not.tOverride_input)) then
!Molpro writes out its own input file
            Call input_options(echo_lines=.false.,skip_blank_lines=.true.)
        else
            Call input_options(echo_lines=iProcIndex.eq.0,skip_blank_lines=.true.)
            Write (6,'(/,64("*"),/)')
        endif


        Do
            Call read_line(tEof)
            If (tEof) exit
            call readu(w)
            select case(w)
            case("TITLE")
                do while (item.lt.nitems)
                    call reada(w)
                    cTitle = trim(cTitle)//" "//trim(w)
                enddo
            case("DEFAULTS")
                CONTINUE
            case("SYSTEM")
                call SysReadInput()
            case("CALC")
                call CalcReadInput()
            case("INTEGRAL")
                call IntReadInput()
            case("LOGGING")
                call LogReadInput()
            case("KP-FCIQMC")
                tKP_FCIQMC = .true.
                tUseProcsAsNodes = .true.
                call kp_fciqmc_read_inp(kp)

            case("END")
                exit
            case default
                call report ("Keyword "//trim(w)//" not recognized",.true.)
            end select
        end do
        write (6,'(/,64("*"),/)')
!        IF(IR.EQ.1.or.IR.EQ.7) CLOSE(ir)
        CLOSE(ir)
   99   IF (ios.gt.0) THEN
            WRITE (6,*) 'Problem reading input file ',TRIM(cFilename)
            call stop_all('ReadInputMain','Input error.')
        END IF
        call checkinput()
        RETURN
    END SUBROUTINE ReadInputMain



    subroutine checkinput()

        ! Check that the specified runtime options are consistent and valid

        use SystemData, only: nel, tUseBrillouin, beta, tFixLz, &
                              tFindCINatOrbs, tNoRenormRandExcits, LMS, STOT,&
                              tCSF, tSpn, tUHF, tGenHelWeighted, tHPHF, &
                              tGen_4ind_weighted, tGen_4ind_reverse, &
                              tMultiReplicas, tGen_4ind_part_exact, &
                              tGen_4ind_lin_exact, tGen_4ind_2, tGAS, tGASSpinRecoupling, &
                              tComplexOrbs_RealInts, tLatticeGens, tHistSpinDist
        use CalcData, only: I_VMAX, NPATHS, G_VMC_EXCITWEIGHT, &
                            G_VMC_EXCITWEIGHTS, EXCITFUNCS, TMCDIRECTSUM, &
                            TDIAGNODES, TSTARSTARS, TBiasing, TMoveDets, &
                            TNoSameExcit, TInitStar, tMP2Standalone, &
                            MemoryFacPart, tSemiStochastic, &
                            tSpatialOnlyHash, InitWalkers, tUniqueHFNode, &
                            tCheckHighestPop, &
                            tKP_FCIQMC, &
                            tRealCoeffByExcitLevel, &
                            tAllRealCoeff, tUseRealCoeffs, tChangeProjEDet, &
                            tOrthogonaliseReplicas, tReadPops, tStartMP1, &
                            tStartCAS, tUniqueHFNode, tContTimeFCIMC, &
                            tContTimeFull, tFCIMC, tPreCond, tOrthogonaliseReplicas, tMultipleInitialStates, tSpinProject, &
                            pgen_unit_test_spec
        use Calc, only : RDMsamplingiters_in_inp
        Use Determinants, only: SpecDet, tagSpecDet, tDefinedet
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
                           tHistEnergies, tPrintOrbOcc
        use Logging, only : calcrdmonfly_in_inp, RDMlinspace_in_inp
        use DetCalc, only: tEnergy, tCalcHMat, tFindDets, tCompressDets
        use load_balance_calcnodes, only: tLoadBalanceBlocks
        use input_neci
        use constants
        use global_utilities
        use spin_project, only: spin_proj_nopen_max
        use FciMCData, only: nWalkerHashes, HashLengthFrac, InputDiagSft
        use hist_data, only: tHistSpawn
        use Parallel_neci, only: nNodes,nProcessors
        use UMatCache, only: tDeferred_Umat2d
        use gasci, only: GAS_specification, is_connected, is_valid, GAS_exc_gen, possible_GAS_exc_gen, operator(==)

        implicit none

        integer :: vv, kk, cc, ierr
        real(dp) :: InputDiagSftSingle
        logical :: check
        character(*), parameter :: t_r='checkinput'

        if(tDiagAllSpaceEver.and..not.tHistSpawn) then
            call stop_all(t_r,"DIAGALLSPACEEVER requires HISTSPAWN option")
        endif
        if(tCalcVariationalEnergy.and..not.tHistSpawn) then
            call stop_all(t_r,"CALCVARIATIONALENERGY requires HISTSPAWN option")
        endif
        if(tCalcVariationalEnergy.and..not.tEnergy) then
           call stop_all(t_r,"CALCVARIATIONALENERGY requires initial FCI calculation")
        endif

        nWalkerHashes=nint(HashLengthFrac*InitWalkers)


        ! Turn on histogramming of fcimc wavefunction in order to find density
        ! matrix, or the orbital occupations
        if (tFindCINatOrbs) tCalcFCIMCPsi = .true.

        ! Used in the FCIMC. We find dets and compress them for later use
        if (tCalcFCIMCPsi .or. tHistSpawn) then
           tFindDets = .true.
           tCompressDets = .true.
        endif

        ! We need to have found the dets before calculating the H mat.
        if (tCalcHMat) tFindDets = .true.

        ! If we are using TNoSameExcit, then we have to start with the star -
        ! the other random graph algorithm cannot remove same excitation
        ! links yet.
        if (tNoSameExcit .and. .not. tInitStar) then
            call report ("If we are using TNoSameExcit, then we have to start&
                         & with the star. The other random graph algorithm &
                         &cannot remove same excitation links yet.", .true.)
        endif

        ! The MoveDets and Biasing algorithms cannot both be used in the
        ! GraphMorph Algorithm.
        if (tBiasing .and. tMoveDets) then
            call report("Biasing algorithm and MoveDets algorithm cannot both&
                        & be used",.true.)
        endif

        ! ..RmRootExcitStarsRootChange must be used with DiagStarStars, and not
        ! with ExcitStarsRootChange
        if (tRmRootExcitStarsRootChange .and. .not. tDiagStarStars) then
            call report("RmRootExcitStarsRootChange can only with used with &
                        &DiagStarStars currently",.true.)
        endif

        if (TRmRootExcitStarsRootChange .and. TExcitStarsRootChange) then
            call report("RmRootExcitStarsRootChange and ExcitStarsRootChange &
                        &cannot both be used as they are both different &
                        &options with diagstarstars", .true.)
        endif

        !..ExcitStarsRootChange must be used with TDiagStarStars
        if (tExcitStarsRootChange .and. .not. tDiagStarStars) then
            call report("ExcitStarsRootChange can only with used with &
                        &DiagStarStars currently", .true.)
        endif

        ! ..TDiagStarStars must be used with TStarStars, and cannot be used
        ! with TCalcExcitStar
        if (tDiagStarStars .and. .not. tStarStars) then
            call report("DiagStarStars must be used with StarStars", .true.)
        endif
        if (tDiagStarStars .and. tCalcExcitStar) then
            call report("DiagStarStars is incompatable with CalcExcitStar", &
                        .true.)
        endif
        if(tDiagStarStars .and. (tNoDoubs .or. tJustQuads)) then
            call report("NoDoubs/JustQuads cannot be used with DiagStarStars &
                        &- try CalcExcitStar")
        endif

        ! ..TNoDoubs is only an option which applied to TCalcExcitStar, and
        ! cannot occurs with TJustQuads.
        if (tNoDoubs .and. .not. tCalcExcitStar) then
            call report("STARNODOUBS is only an option which applied to &
                        &TCalcExcitStar", .true.)
        endif

        if (tNoDoubs .and. tJustQuads) then
            call report("STARNODOUBS and STARQUADEXCITS cannot be applied &
                        &together!", .true.)
        endif

        ! .. TJustQuads is only an option which applies to TCalcExcitStar
        if (tJustQuads.and..not.tCalcExcitStar) then
            call report("STARQUADEXCITS is only an option which applies to &
                        &tCalcExcitStar",.true.)
        endif

        !.. tCalcExcitStar can only be used with tStarStars
        if (tCalcExcitStar.and..not.tStarStars) then
            call report("CalcExcitStar can only be used with StarStars set", &
                        .true.)
        endif

        !.. Brillouin Theorem must be applied when using TStarStars
        if (tStarStars.and..not.tUseBrillouin) then
            call report("Brillouin Theorem must be used when using &
                        &CalcExcitStar", .true.)
        endif

        !.. TQuadValMax and TQuadVecMax can only be used if TLINESTARSTARS set
        if ((tQuadValMax .or. tQuadVecMax) .and. .not. tStarStars) then
            call report("TQuadValMax or TQuadVecMax can only be specified if &
                        &STARSTARS specified in method line", .true.)
        endif

        !.. TQuadValMax and TQuadVecMax cannot both be set
        if (tQuadValMax.and.tQuadVecMax) then
            call report("TQuadValMax and TQuadVecMax cannot both be set", &
                        .true.)
        endif

        !.. TDISCONODES can only be set if NODAL is set in the star methods
        ! section
        if (tDiscoNodes .and. .not. tDiagNodes) then
            call report("DISCONNECTED NODES ONLY POSSIBLE IF NODAL SET IN &
                        &METHOD",.true.)
        endif

        if(tMultipleInitialStates .or. tOrthogonaliseReplicas .or. &
             tPreCond) then
           if (tHistSpawn .or. &
                (tCalcFCIMCPsi .and. tFCIMC) .or. tHistEnergies .or. &
                tHistSpinDist .or. tPrintOrbOcc) &
                call report("HistSpawn and PrintOrbOcc not yet supported for multi-replica with different references"&
                ,.true.)
        endif

        !.. We still need a specdet space even if we don't have a specdet.
        if (.not. associated(SPECDET)) then
            allocate(SPECDET(nel - nFrozen), stat=ierr)
            call LogMemAlloc('SPECDET', nel-nFrozen, 4, t_r, tagSPECDET, ierr)
        endif

        !..   Testing ILOGGING
        !     ILOGGING = 0771
        if (I_VMAX == 0 .and. nPaths /= 0 .and. (.not. tKP_FCIQMC)) then
            call report ('NPATHS!=0 and I_VMAX=0.  VERTEX SUM max level not &
                         &set', .true.)
        endif

        !Ensure beta is set.
        if (beta < 1.0e-6_dp .and. .not. tMP2Standalone) then
            call report("No beta value provided.", .true.)
        endif

        do vv=2,I_VMAX
            g_VMC_ExcitWeights(:,vv)=g_VMC_ExcitWeights(:,1)
            G_VMC_EXCITWEIGHT(vv)=G_VMC_EXCITWEIGHT(1)
        enddo

        !IF THERE IS NO WEIGHTING FUNCTION, ExcitFuncs(10)=.true.
        do vv=1,9
            IF(EXCITFUNCS(vv)) EXCITFUNCS(10)=.false.
        enddo

        if (tNoRenormRandExcits .and. (.not.ExcitFuncs(10))) then
            write(6,*) "Random excitations WILL have to be renormalised, &
                       &since an excitation weighting has been detected."
        ENDIF

        ! if the LMS value specified is not reachable with the number of electrons,
        ! fix this
        if(mod(abs(lms),2).ne.mod(nel,2)) then
           write(6,*) "WARNING: LMS Value is not reachable with the given number of electrons."
           write(6,*) "Resetting LMS"
           LMS = -mod(nel,2)
        endif

        ! Check details for spin projection
        if (tSpinProject) then
            if (tCSF) &
                call stop_all (t_r, "Spin projection must not be used with &
                                    &CSFs")

            if (.not. tSpn) &
                call stop_all (t_r, "SPIN-RESTRICT must be used with SPIN-&
                                    &PROJECT to set the value of S, Ms")

            ! Unless specified, apply spin projection to ALL determinants.
            if (spin_proj_nopen_max == -1) &
                spin_proj_nopen_max = nel

            ! Set the value of STOT as required
            STOT = LMS
        endif

        if (tCalcInstantS2 .or. tCalcInstantS2Init) then
            if (tUHF) &
                call stop_all (t_r, 'Cannot calculate instantaneous values of&
                              & S^2 with UF enabled.')
            write(6,*) 'Enabling calculation of instantaneous S^2 each &
                       &iteration.'
        endif

        if (tUniqueHFNode .and. nProcessors < 2) then
            write(6,*) "nNodes: ",nNodes
            write(6,*) 'nProcessors: ', nProcessors
            call stop_all (t_r, 'At least two nodes required to designate &
                          &a node uniquely to the HF determinant')
        end if

        if (tGenHelWeighted) then
            write(6,*)
            write(6,*) '*** WARNING ***'
            write(6,*) 'Slow HElement biased excitation generators in use.'
            write(6,*) 'NOT FOR PRODUCTION RUNS'
            write(6,*) '***************'
            write(6,*)
        end if

        if (tGen_4ind_weighted .or. tGen_4ind_reverse .or. tGen_4ind_2) then

            ! We want to use UMAT2D...
            tDeferred_Umat2d = .true.

            if (tCSF) &
                call stop_all (t_r, 'Integral weighted excitation generators &
                              &not yet implemented with these keywords')
        end if

        if (tHPHF .and. tUHF) then
            call stop_all(t_r, 'HPHF functions cannot work with UHF')
        end if

#if PROG_NUMRUNS_
        if (tKP_FCIQMC .and. .not. tMultiReplicas) then

            write(6,*) 'Using KPFCIQMC without explicitly specifying the &
                       &number of replica simulations'
            write(6,*) 'Defaulting to using 2 replicas'
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
            write(6,*) 'RDM on fly'

            if (.not. tMultiReplicas) then
                write(6,*) 'unspecified'
                write(6,*) 'Filling RDMs without explicitly specifying the &
                           &number of replica simplations'
                write(6,*) 'Defaulting to using 2 replicas'
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

        if (tRDMOnFly .and. .not. tCheckHighestPop) then
            write(6,*) 'Highest population checking required for calculating &
                       &RDMs on the fly'
            write(6,*) 'If you are seeing this, it is an input parsing error'
            call stop_all(t_r, 'RDMs without CheckHighestPop')
        end if

        if (tSemiStochastic .and. .not.(tAllRealCoeff.and.tUseRealCoeffs)) then
            write(6,*) 'Semi-stochastic simulations only supported when using &
                       &ALLREALCOEFF option'
            call stop_all(t_r, 'Semistochastic without ALLREALCOEFF')
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
                write(6,*) 'Disabling load balancing for single node calculation'
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
            write(6,*) 'Options LZTOT and COMPLEXORBS_REALINTS incompatible'
            write(6,*)
            write(6,*) '1. Using multiple options that filter integrals at runtime is unsupported.'
            write(6,*) '   Only one integral filter may be used at once.'
            write(6,*)
            write(6,*) '2. This is almost certainly not what you intended to do.  LZTOT works using'
            write(6,*) '   abelian symmetries combined with momentum information. COMPLEXORBS_REALINTS'
            write(6,*) '   provides support for non-abelian symmetries in FCIDUMP files produced'
            write(6,*) '   using VASP'
            write(6,*)
            call stop_all(t_r, 'Options incompatible')
        end if

        if (tLatticeGens) then
           if (tGen_4ind_2 .or. tGen_4ind_weighted .or. tGen_4ind_reverse) then
              call stop_all(t_r, "Invalid excitation options")
           end if
        end if


        if (tGAS) then
            if(.not. tDefineDet) then
                call stop_all(t_r, "Running GAS requires a user-defined reference via definedet.")
            endif
            if (.not. is_valid(GAS_specification)) then
                call stop_all(t_r, "GAS specification not valid.")
            end if
            if (is_connected(GAS_specification) .and. .not. tGASSpinRecoupling) then
                call stop_all(t_r, "Running GAS without spin-recoupling requires disconnected spaces.")
            end if
            if (GAS_exc_gen == possible_GAS_exc_gen%DISCONNECTED .and. is_connected(GAS_specification)) then
                call stop_all(t_r, "Running GAS-CI = ONLY_DISCONNECTED requires disconnected spaces.")
            end if
        end if


        if (allocated(pgen_unit_test_spec)) then
            if (.not. tReadPops) then
                call stop_all(t_r, "UNIT-TEST-PGEN requires READPOPS.")
            end if
        end if
    end subroutine checkinput

end Module ReadInput_neci
