
      MODULE SYSREAD
        USE input
        IMPLICIT NONE

        LOGICAL TSTARBIN,TREADINT,THFORDER,TDFREAD,TPBC,TUEG,TCPMD,THUB
        LOGICAL TSPN,TCSF,TPARITY,TUSEBRILLOUIN,TEXCH,TREAL,TTILT
        LOGICAL TSTARSTORE
        LOGICAL TALPHA,TSTAR,TSTOREASEXCITATIONS,TBIN
        INTEGER LMS,STOT,IPARITY(5),NMAXX,NMAXY,NMAXZ,NMSH,COULDAMPORB
        INTEGER iPeriodicDampingType,ISTATE,NEL,ITILTX,ITILTY
        REAL*8 BOX,BOA,COA,FUEGRS,fRc,FCOUL,OrbECutoff,UHUB,BHUB
        REAL*8 ALPHA,FCOULDAMPBETA,FCOULDAMPMU
! Defaults stored in this module
        LOGICAL :: defaults,Feb08

        contains

        SUBROUTINE readinputsys()
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
        INTEGER I
        
!       SYSTEM defaults - leave these as the default defaults
!       Any further addition of defaults should change these after.
          TSTARBIN=.false.
          TREADINT=.false.
          THFORDER=.false.
          TDFREAD=.false.
          TPBC=.false.
          TUEG=.false.
          TCPMD=.false.
          THUB=.false.
          TUEG=.false.
          LMS=0
          TSPN=.false.
          TCSF=.false.
          STOT=0
          TPARITY = .false.
          IParity(:)=0
          NMAXX = 0
          NMAXY = 0
          NMAXZ = 0
          NMSH=32
          BOX=1.d0
          BOA=1.d0
          COA=1.d0
          TUSEBRILLOUIN=.false. 
          FUEGRS=0.D0
          iPeriodicDampingType=0
          fRc=0.D0
          TEXCH=.true.
          FCOUL=1.D0
          UHUB = 4
          BHUB = -1
          TREAL = .false.
          TTILT = .false.
          TALPHA = .false.
          ISTATE = 1
          OrbECutoff=-1e20
          tStar=.false.
          tStoreAsExcitations=.false.
          TBIN=.false.

!Feb08 defaults:
          IF(Feb08) THEN
              !...add defaults...
          ENDIF

! Coulomb damping function currently removed.
!      FCOULDAMPBETA=-1.D0
!      COULDAMPORB=0
          
        call readu(w)
        select case(w)
        case("STARSTOREREAD")
            TSTARSTORE = .true.
            TREADINT = .true.
            call readu(w)
            select case(w)
            case("ORDER")
                THFORDER = .true.
            end select
        case("STARBINREAD")
            TSTARSTORE=.true.
            TBIN=.true.
            TREADINT=.true.
            call readu(w)
            select case(w)
            case("ORDER")
                THFORDER=.true.
            end select
        case("DFREAD")
            TREADINT = .true.
            TDFREAD = .true.
            call readu(w)
            select case(w)
            case("ORDER")
                THFORDER = .true.
            end select
        case("BINREAD")
            TREADINT=.true.
            TBIN=.true.
            call readu(w)
            select case(w)
            case("ORDER")
                THFORDER=.true.
            end select
        case("READ","GENERIC")
            TREADINT = .true.
            call readu(w)
            select case(w)
            case("ORDER")
                THFORDER = .true.
            end select
        case("HUBBARD")
            THUB = .true.
            TPBC=.true.
        case("UEG")
            TUEG = .true.
        case("CPMD")
            TCPMD = .true.
            call readu(w)
            select case(w)
        case("ORDER")
                THFORDER = .true.
            end select
        case default
            call report ("System type "//trim(w)                        &
     &               //" not valid",.true.)
        end select
        

        system: do
          call read_line(eof)
          if (eof) then
              call report("Incomplete input file",.true.)
          end if
          call readu(w)
          select case(w)
          case("ELECTRONS","NEL")
              call geti(NEL)
          case("SPIN-RESTRICT")
              if(item.lt.nitems) then
                 call geti(LMS)
              else
                 LMS=0
              endif
              TSPN = .true.
          case("CSF")
              if(item.lt.nitems) then
                 call geti(STOT)
              else
                 STOT=0
              endif
              TCSF = .true.
          case("SYM")
              TPARITY = .true.
              do I = 1,4
                call geti(IPARITY(I))
              end do
! the last number is the symmetry specification - and is placed in position 5
              IPARITY(5)=IPARITY(4)
              IPARITY(4)=0
          case("CELL")
              call geti(NMAXX)
              call geti(NMAXY)
              call geti(NMAXZ)
          case("MESH")
              call geti(NMSH)
          case("BOXSIZE")
              call getf(BOX)
              if(item.lt.nitems) then
                 call getf(BOA)
                 call getf(COA)
              else
                 BOA=1.D0
                 COA=1.D0
              endif
          case("USEBRILLOUINTHEOREM")
            TUSEBRILLOUIN=.TRUE. 
          case("RS")
              call getf(FUEGRS)
          case("EXCHANGE-CUTOFF")
              iPeriodicDampingType=2
              if(item.lt.nitems) then
                 call getf(fRc)
              endif
          case("EXCHANGE-ATTENUATE")
              iPeriodicDampingType=1
              if(item.lt.nitems) then
                 call getf(fRc)
              endif
          case("EXCHANGE")
              call readu(w)
              select case(w)
                 case("ON")
                    TEXCH=.TRUE.
                 case("OFF")
                    TEXCH=.FALSE.
                 case default
                    call report ("EXCHANGE "//trim(w)                   &
     &                 //" not valid",.true.)
              end select
          case("COULOMB")
              call getf(FCOUL)
          case("COULOMB-DAMPING")
                  call report("Coulomb damping feature removed"         &
     &            ,.true.)
              
!                  call readu(w)
!                  select case(w)
!                  case("ENERGY")
!                     call getf(FCOULDAMPMU)
!                     call getf(FCOULDAMPBETA)
!                  case("ORBITAL")
!                     call geti(COULDAMPORB)
!                     call getf(FCOULDAMPBETA)
!                  end select
          case("U")
              call getf(UHUB)
          case("B")
              call getf(BHUB)
          case("REAL")
              TREAL = .true.
          case("APERIODIC")
              TPBC = .false.
          case("TILT")
              TTILT = .true.
              call geti(ITILTX)
              call geti(ITILTY)
          case("ALPHA")
              TALPHA = .true.
              call getf(ALPHA)
          case("STATE")
              call geti(ISTATE)
              if ( ISTATE /= 1 ) then
                  call report("Require ISTATE to be left set as"        &
     &            //" 1",.true.)
              end if
          case("ENERGY-CUTOFF")
              call getf(OrbECutoff)
          case("STORE-AS-EXCITATIONS")
             tStoreAsExcitations=.true.  

          case("ENDSYS") 
              exit system
          case default
              call report("Keyword "                                    &
     &          //trim(w)//" not recognized in SYSTEM block",.true.)
          end select
        end do system
        if(NEL.eq.0)                                                    &
     &     call report("Number of electrons cannot be zero.",.true.)
        if(THUB.OR.TUEG.OR..NOT.(TREADINT.OR.TCPMD)) then
           if(NMAXX.EQ.0)                                               &
     &        call report("Must specify CELL "                          &
     &        //"- the number of basis functions in each dim.",         &
     &        .true.)
           if(.NOT.THUB.AND.BOX.EQ.0.D0)                                &
     &        call report("Must specify BOX size.",.true.)
           if(TTILT.AND..NOT.THUB)                                      &
     &        call report("TILT can only be specified with "            &
     &     //"HUBBARD.",.true.)
        endif

        END SUBROUTINE readinputsys

      END MODULE SYSREAD
