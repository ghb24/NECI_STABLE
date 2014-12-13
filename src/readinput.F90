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
        use kp_fciqmc_procs, only: kp_fciqmc_read_inp, kp_fciqmc_data
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
        Character(len=64)  cFilename    !Input  filename or "" if we check arg list or stdin
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
        type(kp_fciqmc_data), intent(inout) :: kp
        
        cTitle=""
        idDef=idDefault                 !use the Default defaults (pre feb08)
        ir=get_free_unit()              !default to a free unit which we'll open below
        If(cFilename.ne.'') Then
            Write(6,*) "Reading from file: ", Trim(cFilename)
            inquire(file=cFilename,exist=tExists)
            if (.not.tExists) call stop_all('ReadInputMain','File '//Trim(cFilename)//' does not exist.')
            Open(ir,File=cFilename,Status='OLD',err=99,iostat=ios)
        ElseIf(neci_iArgC().gt.0) then
    ! We have some arguments we can process instead
            Call neci_GetArg(1,cInp)      !Read argument 1 into inp
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

        use SystemData, only: nel, tStarStore, tUseBrillouin, beta, tFixLz, &
                              tFindCINatOrbs, tNoRenormRandExcits, LMS, STOT,&
                              tCSF, tSpn, tUHF, tGenHelWeighted, tHPHF, &
                              tGen_4ind_weighted, tGen_4ind_reverse
        use CalcData, only: I_VMAX, NPATHS, G_VMC_EXCITWEIGHT, &
                            G_VMC_EXCITWEIGHTS, EXCITFUNCS, TMCDIRECTSUM, &
                            TDIAGNODES, TSTARSTARS, TBiasing, TMoveDets, &
                            TNoSameExcit, TInitStar, tMP2Standalone, &
                            MemoryFacPart, tTruncInitiator, tSemiStochastic, &
                            tSpatialOnlyHash, InitWalkers, tUniqueHFNode, &
                            InitiatorCutoffEnergy, tCCMC, &
                            tSurvivalInitiatorThreshold, tKP_FCIQMC, &
                            tSurvivalInitMultThresh, tAddToInitiator, &
                            tMultiReplicaInitiators
        Use Determinants, only: SpecDet, tagSpecDet
        use IntegralsData, only: nFrozen, tDiscoNodes, tQuadValMax, &
                                 tQuadVecMax, tCalcExcitStar, tJustQuads, &
                                 tNoDoubs
        use IntegralsData, only: tDiagStarStars, tExcitStarsRootChange, &
                                 tRmRootExcitStarsRootChange, tLinRootChange
        use LoggingData, only: iLogging, tCalcFCIMCPsi, &
                           tCalcInstantS2, tDiagAllSpaceEver, &
                           tCalcVariationalEnergy, tCalcInstantS2Init, &
                           tPopsFile
        use DetCalc, only: tEnergy, tCalcHMat, tFindDets, tCompressDets
        use input_neci
        use constants
        use global_utilities
        use spin_project, only: tSpinProject, spin_proj_nopen_max
        use FciMCData, only: nWalkerHashes,HashLengthFrac,tHashWalkerList
        use hist_data, only: tHistSpawn
        use Parallel_neci, only: nNodes,nProcessors
        use UMatCache, only: tDeferred_Umat2d

        implicit none

        integer :: vv, kk, cc, ierr
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

        if(tHashWalkerList) then
            nWalkerHashes=nint(HashLengthFrac*InitWalkers)
        endif


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
  
        !IF FINDD or USED specified without using Excitweighting option
        if ((I_VMAX >= 3) .and. (tStarStore)) then 
            call report("Error - can only use STARSTOREREAD with double &
                        &excitations of HF",.true.)
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

        if (tGen_4ind_weighted .or. tGen_4ind_reverse) then

            ! We want to use UMAT2D...
            tDeferred_Umat2d = .true.

            if (tCSF) &
                call stop_all (t_r, 'Integral weighted excitation generators &
                              &not yet implemented with these keywords')
        end if

        if (tHPHF .and. tUHF) then
            call stop_all(t_r, 'HPHF functions cannot work with UHF')
        end if

        if (tCCMC .and. .not. (InitiatorCutoffEnergy > 1.0e99_dp)) then
            call stop_all(t_r, 'Initiator cutoff not implemented for CCMC')
        end if

        if (tPopsFile .and. (tSurvivalInitiatorThreshold .or. &
                             tSurvivalInitMultThresh)) then
            write(6,*) 'The initiator initial iteration details have not yet &
                       &been added to the POPSFILE reading/writing routines.'
            write(6,*) '--> Simulations will not display consistent behaviour &
                       &over restarts until this is implemented'
            call stop_all(t_r, 'POPSFILES cannot be used with initiator &
                               &survival criterion')
        end if

        if ((tSurvivalInitiatorThreshold .or. tSurvivalInitMultThresh) .and. &
            .not. tAddToInitiator) then

            write(6,*) 'Without the ADDTOINITIATOR option, the survival based &
                       &initiator thresholds do nothing'
            write(6,*) 'If ONLY survival based thresholds are desired, set &
                       &initiator walker number absurdly high'
            call stop_all(t_r, 'Inconsistent options')
        end if

        if (tMultiReplicaInitiators) then
#ifndef __PROG_NUMRUNS
            call stop_all(t_r, 'Aggregated initiator thresholds require &
                               &multiple simulations')
#endif
#ifdef __CMPLX
            call stop_all(t_r, 'Aggregated initator thresholds are not (yet) &
                               &implemented for complex particles')
#endif
            if (lenof_sign == 1) &
                call stop_all(t_r, 'Aggregated initator thresholds make no &
                                   &sense with only one system replica')
        end if

    end subroutine checkinput

end Module ReadInput_neci

        

