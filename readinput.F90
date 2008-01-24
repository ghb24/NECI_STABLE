      SUBROUTINE READINPUT(FILENAME)
      USE input
      USE SYSREAD , only : readinputsys
      USE PRECALCREAD , only : readinputprecalc
      USE CALCREAD , only : readinputcalc
      USE INTREAD , only : readinputint
      USE LOGREAD , only : readinputlog
#ifdef NAGF95
USE f90_unix_env, ONLY: getarg,iargc
#endif
      IMPLICIT NONE
#ifndef NAGF95
      integer :: iargc
#endif
      CHARACTER(LEN=255) FILENAME,INP
      CHARACTER(LEN=32) TITLE
      CHARACTER(LEN=100) w
      LOGICAL eof 
      INTEGER ios

!     --------------------------------------------
      title = ""

      ir=1
      IF(FILENAME.NE.'') THEN
         WRITE(6,*) "Reading from file: ",TRIM(FILENAME)
         OPEN(1,FILE=FILENAME,STATUS="OLD",err=99,iostat=ios)
      ELSEIF(IARGC().GT.0) THEN
         CALL GETARG(1,INP)
         WRITE(6,*) "Reading from file: ",INP
         OPEN(1,FILE=INP,STATUS="OLD")
      ELSE
         ir=5
         WRITE(6,*) "Reading from STDIN"
      ENDIF

      write (6,'(/,64("*"),/)')
      call input_options(echo_lines=.true.,skip_blank_lines=.true.)
      main: do
          call read_line(eof)
          if (eof) exit
          call readu(w)
          select case(w)
            case("TITLE")
              do while (item.lt.nitems)
                  call reada(w)
                  title = trim(title)//" "//trim(w)
              enddo
            case("SYSTEM")
              call readinputsys()
            case("PRECALC")
              call readinputprecalc()
            case("CALC")
              call readinputcalc()
            case("INTEGRAL")
              call readinputint()
            case("LOGGING")
              call readinputlog()
            case("END")
              exit
            case default
             call report ("Keyword "//trim(w)//" not recognized",.true.)
          end select
      end do main
      write (6,'(/,64("*"),/)')
      IF(IR.EQ.1) CLOSE(1)
   99 IF (ios.gt.0) THEN
          WRITE (6,*) 'Problem reading input file ',TRIM(FILENAME)
      END IF
      call checkinput()
      RETURN
      END SUBROUTINE readinput
      
      subroutine inpgetprecalc(preIH)
         use input
         USE SYSREAD , only : TUSEBRILLOUIN
         implicit none
         integer preIH
         CHARACTER(LEN=16) w
                do while ( item .lt. nitems )
                  call readu(w)
                  select case(w)
                  case("HDIAG")
                      call readu(w)
                      select case(w)
                      case("FULL")
                          preIH=-20
                      case("MC")
                          preIH=-19
                          IF(.NOT.TUSEBRILLOUIN) THEN
                              write(6,*) "Warning  USEBRILLOUINTHEOREM" &
     &                       //" might need to be specified in system " &
     &                       //"block to use MC-PRECALC"
                          ENDIF
                      case default
                         call report("Error - must specify FULL"        &
     &                   //" or MC after HDIAG in PRECALC block",.true.)
                      end select
                   case("RHODIAG")
                         call readu(w)
                         select case(w)
                         case("FULL")
                             preIH=-8
                         case("MC")
                             preIH=-7
                         case default
                           call report("Error - must specify FULL or "  &
     &                     //"MC after RHODIAG in PRECALC block",.true.)
                         end select
                     case default
                           call report("Keyword error with "//trim(w),  &
     &                     .true.)
                     end select
                 end do
        end

      subroutine inpgetmethod(I_HMAX,NWHTAY,I_V)
         use input
         use UMatCache , only : TSTARSTORE
         USE CALCREAD , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,         &
     &          TMCDIRECTSUM,g_Multiweight,G_VMC_FAC,TMPTHEORY,         &
     &          STARPROD
         implicit none
         integer I_HMAX,NWHTAY,I_V
         CHARACTER(LEN=16) w
                  do while ( item .lt. nitems )
                    call readu(w)
                    select case(w)
                    case("VERTEX")
                        call readu(w)
                        select case(w)
                        case("SUM")
                           do while(item.lt.nitems)
                            call readu(w)
                            select case(w)
                            case("OLD")
                                I_HMAX = -1
                            case("NEW")
                                I_HMAX = -8
                            case("HDIAG")
                                I_HMAX = -20
                            case("READ")
                                I_HMAX=-14
                            case("SUB2VSTAR")
                                CALCP_SUB2VSTAR=.TRUE.
                            case("LOGWEIGHT")
                                CALCP_LOGWEIGHT=.TRUE.
                            case default
                                call report("Error - must specify OLD"  &
     &                         //" or NEW vertex sum method",.true.)
                            end select
                           enddo
                        case("MC","MCMETROPOLIS")
                           I_HMAX = -7
                            call readu(w)
                            select case(w)
                            case("HDIAG")
                                I_HMAX = -19
                            end select
                           tMCDirectSum=.FALSE.
                           IF(I_V.GT.0) g_MultiWeight(I_V)=1.D0
                        case("MCDIRECT")
                           I_HMAX = -7
                           tMCDirectSum=.TRUE.
                            call readu(w)
                            select case(w)
                            case("HDIAG")
                                I_HMAX = -19
                            end select
                           G_VMC_FAC=0.D0
                        case("MCMP")
                           tMCDirectSum=.TRUE.
                           G_VMC_FAC=0.D0
                           TMPTHEORY=.TRUE.
                        case("STAR")
                           I_HMAX=0
                           do while(item.lt.nitems)
                              call readu(w)
                              select case(w)
                              case("NEW")
                                 I_HMAX=-21
                              case("OLD")
                                 I_HMAX=-9
                              case("STARPROD")
                                 STARPROD=.TRUE.
                              case("COUNTEXCITS")
                                 NWHTAY=IBSET(NWHTAY,8)
                              case("ADDSINGLES")
                                 NWHTAY=IBSET(NWHTAY,7)
                                 IF(I_HMAX.NE.-21)  call report(        &
     &                              "Error - cannot use ADDSINGLES"     &
     &                              //" without STAR NEW",.true.)
                                 IF(TSTARSTORE) call report("Error - "  &
     &                            //"can only use STARSTOREREAD with "  &
     &                            //"double excitations of HF",.true.)
                              case("DIAG")
                                  NWHTAY=IBCLR(NWHTAY,0)
                              case("POLY")
                                  NWHTAY=IBSET(NWHTAY,0)
                              case("POLYMAX")
                                  NWHTAY=IBSET(NWHTAY,0)
                                  NWHTAY=IBSET(NWHTAY,1)
                              case("POLYCONVERGE")
                                  NWHTAY=IBSET(NWHTAY,0)
                                  NWHTAY=IBSET(NWHTAY,2)
                              case("POLYCONVERGE2")
                                  NWHTAY=IBSET(NWHTAY,0)
                                  NWHTAY=IBSET(NWHTAY,6)
                              case("H0")
                                  NWHTAY=IBSET(NWHTAY,5)
                                  if(I_HMAX.ne.-21) call report ("H0 "  &
     &                       //"can only be specified with POLY... NEW")
                              case default
                                call report("Error - must specify DIAG" &
     &                        //" or POLY vertex star method",.true.)
                               end select
                           enddo
                           IF(STARPROD.and.BTEST(NWHTAY,0)) THEN
                               call report("STARPROD can only be "      &
     &                        //"specified with HDIAG option",.true.)
                            ENDIF
                           if(i_hmax.eq.0)                              &
     &                   call report("OLD/NEW not specified for STAR",  &
     &                          .true.)
                        case default
                        call report("Keyword error with "//trim(w),     &
     &                          .true.)
                        end select
                    case default
                        call report("Error.  Method not specified."     &
     &                    //" Stopping.",.true.)
                    end select
               end do
      end

        
      subroutine inpgetexcitations(NWHTAY,w)
         use input
         IMPLICIT NONE
         INTEGER NWHTAY
         CHARACTER(LEN=16) w
!                  call readu(w)
                  select case(w)
                  case("FORCEROOT")
                     NWHTAY=IOR(NWHTAY,1)
                  case("FORCETREE")
                     NWHTAY=IOR(NWHTAY,2)
                  case("SINGLES")
                     NWHTAY=IOR(NWHTAY,8)
                  case("DOUBLES")
                     NWHTAY=IOR(NWHTAY,16)
                  case("ALL")
                     NWHTAY=0
                  case default
                        call report("Keyword error with EXCITATIONS "   &
     &                     //trim(w),                                   &
     &                          .true.)
                  end select
      end

      subroutine checkinput()
      USE SYSREAD , only : NEL,TSTARSTORE
      USE PRECALCREAD , only : PREIV_MAX,USEVAR,PRE_TAYLOG,             &
     &  TGRIDVAR,TLINEVAR,TOTALERROR,TRUECYCLES
      USE CALCREAD , only : BETA,I_VMAX,NPATHS,SPECDET,                 &
     &  G_VMC_EXCITWEIGHT,G_VMC_EXCITWEIGHTS,EXCITFUNCS,TMCDIRECTSUM
      USE INTREAD , only : NFROZEN
      USE LOGREAD , only : ILOGGING
      USE input
      IMPLICIT NONE
      INTEGER :: vv,kk,cc,ierr
      LOGICAL :: CHECK
      
!.. We still need a specdet space even if we don't have a specdet.
      IF(.NOT.ALLOCATED(SPECDET)) THEN
          ALLOCATE(SPECDET(NEL-NFROZEN),STAT=ierr)
          CALL MemAlloc(ierr,SPECDET,NEL-NFROZEN,'SPECDET')
      ENDIF
!      IF(IP_SPECDET.EQ.0) call MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')

!..   Testing ILOGGING
!     ILOGGING = 0771
      IF(I_VMAX.EQ.0.AND.NPATHS.NE.0)                                   &
     &   STOP 'NPATHS!=0 and I_VMAX=0.  VERTEX SUM max level not set'
      WRITE (6,"(A,Z4)") 'ILOGGING after input routine', ILOGGING

      !Ensure beta is set.
      if (beta.lt.1.d-6) call report("No beta value provided.",.true.)
      
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

      !IF THERE IS NO WEIGHTING FUNCTION
      do vv=1,9
          IF(EXCITFUNCS(vv)) EXCITFUNCS(10)=.false.
      enddo
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
