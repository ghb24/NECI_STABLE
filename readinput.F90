! ReadInput is called to read the input.
!   Filename is a Character(255) string which contains a filename to read.
!               If Filename=="" then we check to see if there's a filename on the command line.
!               Failing that, we use stdin
!   ios is an Integer which is set to 0 on a successful return, or is non-zero if a file error has occurred, where it is the iostat.
MODULE ReadInput 
    Implicit none
!   Used to specify which default set of inputs to use
!    An enum would be nice, but is sadly not supported
    integer, parameter :: idDefault=0
    integer, parameter :: idFeb08=1

    contains

    Subroutine ReadInputMain(cFilename,ios)
        USE input
        use System,     only : SysReadInput,SetSysDefaults
        USE Precalc,    only : PrecalcReadInput,SetPrecalcDefaults
        use Calc,       only : CalcReadInput,SetCalcDefaults
        use Integrals,  only : IntReadInput,SetIntDefaults
        Use Logging,    only : LogReadInput,SetLogDefaults
        use Parallel,   only : iProcIndex
        use default_sets
#ifdef NAGF95
    !  USe doesn't get picked up by the make scripts
        USe f90_unix_env, ONLY: getarg,iargc
#endif
        Implicit none
#ifndef NAGF95
        Integer :: iargc
#endif
    !  INPUT/OUTPUT params
        Character(len=255)  cFilename    !Input  filename or "" if we check arg list or stdin
        Integer             ios         !Output 0 if no error or nonzero iostat if error

        Character(len=255)  cInp         !temp storage for command line params
        Character(len=32)   cTitle
!  Predeclared
!        Integer             ir         !The file descriptor we are reading from
        Character(len=100)  w,x         !strings for input storage
        Logical             tEof        !set when read_line runs out of lines
        Integer             idDef       !What default set do we use
        
        cTitle=""
        idDef=idDefault                 !use the Default defaults (pre feb08)
        ir=1                            !default to file descriptor 1 which we'll open below
        If(cFilename.ne.'') Then
            Write(6,*) "Reading from file: ", Trim(cFilename)
            Open(1,File=cFilename,Status='OLD',err=99,iostat=ios)
        ElseIf(iArgC().gt.0) then
    ! We have some arguments we can process instead
            Call GetArg(1,cInp)      !Read argument 1 into inp
            Write(6,*) "Reading from file: ", Trim(cInp)
            Open(1,File=cInp,Status='OLD',FORM="FORMATTED",err=99,iostat=ios)
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
        case(1)
            Feb08=.true.
            write (6,*) 'Using the Feb08 set of defaults.'
        end select

        ! Set up defaults.
        call SetSysDefaults
        call SetPrecalcDefaults
        call SetCalcDefaults
        call SetIntDefaults
        call SetLogDefaults

!Now return to the beginning and process the whole input file
        if (ir.eq.5) ir=7 ! If read from STDIN, re-read from our temporary scratch file.
        Rewind(ir)
        Call input_options(echo_lines=iProcIndex.eq.0,skip_blank_lines=.true.)

        Write (6,'(/,64("*"),/)')

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
            case("PRECALC")
                call PrecalcReadInput()
            case("CALC")
                call CalcReadInput()
            case("INTEGRAL")
                call IntReadInput()
            case("LOGGING")
                call LogReadInput()
            case("END")
                exit
            case default
                call report ("Keyword "//trim(w)//" not recognized",.true.)
            end select
        end do
        write (6,'(/,64("*"),/)')
        IF(IR.EQ.1.or.IR.EQ.7) CLOSE(ir)
   99   IF (ios.gt.0) THEN
            WRITE (6,*) 'Problem reading input file ',TRIM(cFilename)
        END IF
        call checkinput()
        RETURN
    END SUBROUTINE ReadInputMain



      subroutine checkinput()
      use SystemData , only : NEL,TSTARSTORE,TUseBrillouin, Beta
      USE PrecalcData , only : PREIV_MAX,USEVAR,PRE_TAYLOG,             &
     &  TGRIDVAR,TLINEVAR,TOTALERROR,TRUECYCLES
      use CalcData , only : I_VMAX,NPATHS,                 &
     &  G_VMC_EXCITWEIGHT,G_VMC_EXCITWEIGHTS,EXCITFUNCS,TMCDIRECTSUM,   &
     & TDIAGNODES,TSTARSTARS,TBiasing,TMoveDets,TNoSameExcit,TInitStar,tMP2Standalone, &
     & GrowMaxFactor,MemoryFac
      Use Determinants, only : SpecDet,tagSpecDet
      use IntegralsData , only : NFROZEN,TDISCONODES,TQuadValMax,TQuadVecMax,TCalcExcitStar,TJustQuads,TNoDoubs,TDiagStarStars,TExcitStarsRootChange,TRmRootExcitStarsRootChange,TLinRootChange
      USE Logging , only : ILOGGING
      use SystemData, only : TNoRenormRandExcits
      USE input
      use global_utilities
      IMPLICIT NONE
      INTEGER :: vv,kk,cc,ierr
      LOGICAL :: CHECK
      character(*), parameter :: t_r='checkinput'

      IF(GrowMaxFactor.gt.MemoryFac) THEN
          CALL report("GrowMaxFactor is larger than MemoryFac - there will not be enough memory allocated if the walker number grows large. Think about increasing MemoryFac or reducing GrowMaxFactor.",.true.)
      ENDIF

!If we are using TNoSameExcit, then we have to start with the star - the other random graph algorithm cannot remove same excitation links yet.
      IF(TNoSameExcit.and..not.TInitStar) THEN
          CALL report("If we are using TNoSameExcit, then we have to start with the star - the other random graph algorithm cannot remove same excitation links yet.",.true.)
      ENDIF

!The MoveDets and Biasing algorithms cannot both be used in the GraphMorph Algorithm.
      IF(TBiasing.and.TMoveDets) THEN
          CALL report("Biasing algorithm and MoveDets algorithm cannot both be used",.true.)
      ENDIF

!..RmRootExcitStarsRootChange must be used with DiagStarStars, and not with ExcitStarsRootChange
      IF(TRmRootExcitStarsRootChange.and..not.TDiagStarStars) THEN
          CALL report("RmRootExcitStarsRootChange can only with used with DiagStarStars currently",.true.)
      ENDIF

      IF(TRmRootExcitStarsRootChange.and.TExcitStarsRootChange) THEN
          CALL report("RmRootExcitStarsRootChange and ExcitStarsRootChange cannot both be used as they are both different options with diagstarstars",.true.)
      ENDIF
      
!..ExcitStarsRootChange must be used with TDiagStarStars
      IF(TExcitStarsRootChange.and..not.TDiagStarStars) THEN
          CALL report("ExcitStarsRootChange can only with used with DiagStarStars currently",.true.)
      ENDIF
      
!..TDiagStarStars must be used with TStarStars, and cannot be used with TCalcExcitStar
      IF(TDiagStarStars.and..not.TStarStars) THEN
          CALL report("DiagStarStars must be used with StarStars",.true.)
      ENDIF
      IF(TDiagStarStars.and.TCalcExcitStar) THEN
          CALL report("DiagStarStars is incompatable with CalcExcitStar",.true.)
      ENDIF
      IF(TDiagStarStars.and.(TNoDoubs.or.TJustQuads)) THEN
          CALL report("NoDoubs/JustQuads cannot be used with DiagStarStars - try CalcExcitStar")
      ENDIF

!.. TNoDoubs is only an option which applied to TCalcExcitStar, and cannot occurs with TJustQuads.
      IF(TNoDoubs.and..not.TCalcExcitStar) THEN
          CALL report("STARNODOUBS is only an option which applied to TCalcExcitStar",.true.)
      ENDIF

      IF(TNoDoubs.and.TJustQuads) THEN
          CALL report("STARNODOUBS and STARQUADEXCITS cannot be applied together!",.true.)
      ENDIF
      
!.. TJustQuads is only an option which applies to TCalcExcitStar
      IF(TJustQuads.and..not.TCalcExcitStar) THEN
          CALL report("STARQUADEXCITS is only an option which applies to TCalcExcitStar",.true.)
      ENDIF
      
!.. TCalcExcitStar can only be used with TSTARSTARS
      IF(TCalcExcitStar.and..not.TSTARSTARS) THEN
          CALL report("CalcExcitStar can only be used with StarStars set",.true.)
      ENDIF

!.. Brillouin Theorem must be applied when using TStarStars
      IF(TStarStars.and..not.TUseBrillouin) THEN
          CALL report("Brillouin Theorem must be used when using CalcExcitStar",.true.)
      ENDIF

!.. TQuadValMax and TQuadVecMax can only be used if TLINESTARSTARS set
      IF((TQuadValMax.or.TQuadVecMax).and..not.TSTARSTARS) THEN
          CALL report("TQuadValMax or TQuadVecMax can only be specified if STARSTARS specified in method line",.true.)
      ENDIF

!.. TQuadValMax and TQuadVecMax cannot both be set
      IF(TQuadValMax.and.TQuadVecMax) THEN
          CALL report("TQuadValMax and TQuadVecMax cannot both be set",.true.)
      ENDIF
      
!.. TDISCONODES can only be set if NODAL is set in the star methods section
      IF(TDISCONODES.AND..NOT.TDIAGNODES) THEN
          CALL report("DISCONNECTED NODES ONLY POSSIBLE IF NODAL SET IN METHOD",.true.)
      ENDIF
      
!.. We still need a specdet space even if we don't have a specdet.
      IF(.NOT.ASSOCIATED(SPECDET)) THEN
          ALLOCATE(SPECDET(NEL-NFROZEN),STAT=ierr)
          CALL LogMemAlloc('SPECDET',NEL-NFROZEN,4,t_r,tagSPECDET,ierr)
      ENDIF
!      IF(IP_SPECDET.EQ.0) call MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')

!..   Testing ILOGGING
!     ILOGGING = 0771
      IF(I_VMAX.EQ.0.AND.NPATHS.NE.0)                                   &
     &   STOP 'NPATHS!=0 and I_VMAX=0.  VERTEX SUM max level not set'
      WRITE (6,"(A,Z4)") 'ILOGGING after input routine', ILOGGING

      !Ensure beta is set.
      if (beta.lt.1.d-6.and..not.tMP2Standalone) call report("No beta value provided.",.true.)
      
      !Make sure there aren't more precalc levels than true vertex levels - problems otherwise
      IF(preIV_MAX.gt.I_VMAX) THEN
          CALL report("There cannot be more precalc vertex levels "     &
     &     //"than vertex levels in the main program",.true.)
      ENDIF
      
      !We make sure that in precalc, a use statement is not specified more than once for any vertex level
      IF(preIV_MAX.ne.0) THEN
          do vv=2,I_VMAX
            check=.false.
            do kk=2,preIV_MAX
                do cc=1,8
                    IF((USEVAR(kk,cc).eq.vv).and.(check)) THEN
                     CALL report("Can only specify to use precalc "     &
     &               //"parameters on a given vertex level once",.true.)
                    ENDIF
                    IF(USEVAR(kk,cc).eq.vv) check=.true.
                enddo
            enddo
            !If use isn't specified for a vertex level, use the values given in the input file
!           Done later now
!            IF(.not.check) THEN
!                g_VMC_ExcitWeights(:,vv)=g_VMC_ExcitWeights(:,1)
!                G_VMC_EXCITWEIGHT(vv)=G_VMC_EXCITWEIGHT(1)
!            ENDIF
        enddo
      ENDIF
      
      !If not doing precalc, set all weighting parameters to the ones in the input file
      IF(preIV_MAX.eq.0) THEN
        do vv=2,I_VMAX
            g_VMC_ExcitWeights(:,vv)=g_VMC_ExcitWeights(:,1)
            G_VMC_EXCITWEIGHT(vv)=G_VMC_EXCITWEIGHT(1)
        enddo
      ENDIF

      !IF THERE IS NO WEIGHTING FUNCTION, ExcitFuncs(10)=.true.
      do vv=1,9
          IF(EXCITFUNCS(vv)) EXCITFUNCS(10)=.false.
      enddo

      IF(TNoRenormRandExcits.and.(.not.ExcitFuncs(10))) THEN
          WRITE(6,*) "Random excitations WILL have to be renormalised, "&
     &      //"since an excitation weighting has been detected."
      ENDIF

      !IF FINDD or USED specified without using Excitweighting option
      do vv=2,preIV_MAX
!    IF((pre_TAY(1,vv).eq.-20).and.((NWHTAY(1,vv).eq.-7).or.        &
!&    (NWHTAY(1,vv).eq.-19))) THEN
!    CALL report("Full precalc cannot be used on a vertex level"
!&   //" which is only sampled using MC in the main program",.true.)
!     ENDIF
          IF((pre_TAYLOG(5,vv).or.pre_TAYLOG(6,vv)).and.                &
     &         (.not.EXCITFUNCS(1))) THEN
               CALL report("Logging keyword FINDD and USED"             &
     &         //" can only be used with EXCITWEIGHTING",.true.)
          ENDIF
          IF(TGRIDVAR(vv).and.((.not.EXCITFUNCS(1)).and.(.not.          &
     &         EXCITFUNCS(4)))) THEN
               CALL report("GRIDVAR option only available with"         &
     &          //" excitation functions with two variables",.true.)
          ENDIF
          IF(TLINEVAR(vv).and.((.not.pre_TAYLOG(3,vv)).and.             &
     &              (.not.pre_TAYLOG(2,vv)))) THEN
              CALL report("LINEVAR option only available with"          &
     &          //" FINDIMPORT or FINDC",.true.)
          ENDIF
      ENDDO
      IF((preIV_MAX.ne.0).AND.(.NOT.TMCDIRECTSUM)) THEN
          CALL report("Precalculation can only work with the"           &
     &          //" MCDIRECTSUM option enabled",.true.)
      ENDIF
      IF((TOTALERROR.ne.0.D0).AND.(TRUECYCLES.ne.0)) THEN
          CALL report("Only TRUECYCLES or TOTALERROR can be"            &
     &      //" specified in precalc block",.true.)
      ENDIF
      IF((TOTALERROR.ne.0.D0).AND.(preIV_MAX.ne.I_VMAX)) THEN
          CALL report("TOTALERROR can only be used if the precalc"      &
     &    //" levels are equal to the main block vertex levels",.true.)
      ENDIF
      IF((I_VMAX.ge.3).and.(TSTARSTORE)) THEN 
          call report("Error - can only use STARSTOREREAD with "        &
     &    //"double excitations of HF",.true.)
      ENDIF

      end subroutine checkinput
End Module ReadInput

        

