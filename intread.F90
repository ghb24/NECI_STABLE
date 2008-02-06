      MODULE INTREAD
        USE input
        USE SYSREAD , only : NEL,defaults,Feb08,TUSEBRILLOUIN
        use UMatCache, only: tReadInCache,nSlotsInit,nMemInit,iDumpCacheFlag,iDFMethod
        IMPLICIT NONE

        LOGICAL TQUADRHO,TEXPRHO,THFBASIS,THFCALC,TCALCREALPROD
        LOGICAL TRHF,TReadTUMat,TReadHF,TQUASIEXCIT
        LOGICAL TSUMPROD,TCALCRHOPROD
        
        INTEGER NTAY(2),nHFit,NFROZEN,NTFROZEN,ORBORDER(8,2)
        INTEGER NRSTEPSMAX,IHFMETHOD
        
        REAL*8 NRCONV,RFCONV,ORBORDER2(8)
        REAL*8 HFMix,HFEDelta,HFCDelta

        contains

        SUBROUTINE readinputint()
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
        INTEGER :: i
           
! Integral defaults
      NRCONV=1.D-13
      RFCONV=1.D-8
      NRSTEPSMAX=50
      TQUADRHO=.false.
      TEXPRHO=.false.
      NTAY(1:2) = 1
      THFBASIS=.false.
      THFCALC = .false.
      nHFit=0
      HFMix=0.d0
      HFEDelta=0.d0
      HFCDelta=0.d0
      IHFMETHOD = 1
      TRHF=.true.
      TReadTUMat=.false.
      TReadHF=.false.
      NFROZEN=0
      NTFROZEN=0
      OrbOrder(:,:)=0
      OrbOrder2(:)=0.d0
      nSlotsInit=1024
      nMemInit=0
      iDumpCacheFlag=0
      tReadInCache=.false.
      iDFMethod=0

!Feb 08 defaults
      IF(Feb08) THEN
         NTAY(2)=3
      ENDIF
      
      
        integral: do
          call read_line(eof)
          if (eof) then
              exit
          end if
          call readu(w)
          select case(w)
          case("NRCONV")
              call readf(NRCONV)
          case("RFCONV")
              call readf(RFCONV)
          case("NRSTEPSMAX")
              call readi(NRSTEPSMAX)
          case("INCLUDEQUADRHO")
              TQUADRHO=.true.
          case("EXPRHO")
              TEXPRHO=.true.
          case("RHO-1STORDER")
              NTAY(2)=4
          case("FOCK-PARTITION")
              NTAY(2)=2
          case("FOCK-PARTITION-LOWDIAG")
              NTAY(2)=3
          case("FOCK-PARTITION-DCCORRECT-LOWDIAG")
              NTAY(2)=5
          case("DIAG-PARTITION")
              NTAY(2)=1
          case("CALCREALPROD")
              TCALCREALPROD=.TRUE.
              IF(.NOT.TUSEBRILLOUIN) THEN
                call report(trim(w)//" will not work unless "           &
     &          //"USEBRILLOUINTHEOREM set",.true.)
              ENDIF
          case("QUADEXCITSTAR")
              TQUASIEXCIT=.TRUE.
          case("CALCRHOPROD")
              TCALCRHOPROD=.TRUE.
          case("SUMPRODII")
              TSUMPROD=.TRUE.
          case("HF")
              THFBASIS = .true.
          case("CALCULATE")
              THFCALC = .true.
          case("MAXITERATIONS")
              call geti(NHFIT)
          case("MIX")
              call getf(HFMIX)
          case("THRESHOLD")
              do while ( item .lt. nitems )
                call readu(w)
                select case(w)
                case("ENERGY")
                    call readf(HFEDELTA)
                case("ORBITAL")
                    call readf(HFCDELTA)
                case default
                    call report(trim(w)//" not valid THRESHOLD"         &
     &           //"OPTION.  Specify ENERGY or ORBITAL convergence"     &
     &             //" threshold.",.true.)
                end select
              end do
          case("RHF")
              TRHF = .true.
          case("UHF")
              TRHF = .false.
          case("HFMETHOD")
              call readu(w)
              select case(w)
              case("DESCENT")
                  call readu(w)
                  select case(w)
                  case("OTHER")
                      IHFMETHOD = 2
                  case("SINGLES")
                      IHFMETHOD = 1
                  case default
                      call report(trim(w)//" not valid DESCENT"         &
     &                //" option",.true.)
                  end select
              case("STANDARD")
                  IHFMETHOD = 0
              case("MODIFIED")
                 IHFMETHOD=3
              case default
                  call report(trim(w)//" not valid HF method",          &
     &             .true.)
              end select
          case("READ")
              do while ( item .lt. nitems )
                call readu(w)
                select case(w)
                case("MATRIX")
                    TREADTUMAT = .true.
                case("BASIS")
                    TREADHF = .true.
                case default
                    call report(trim(w)//" is an invalid HF read"       &
     &              //" option.",.true.)
                end select
              end do
          case("FREEZE")
              call readi(NFROZEN)
              call readi(NTFROZEN)
              if ( mod(NFROZEN,2).ne.0 .or.                             &
     &         (NTFROZEN.GT.0 .and. mod(NTFROZEN,2).ne.0) ) then
                  call report("NFROZEN and (+ve) NTFROZEN must be"      &
     &            //"multiples of 2",.true.)
              end if
              if (                                                      &
     &         (NTFROZEN.LT.0 .and. mod(NEL-NTFROZEN,2).ne.0) ) then
                  call report("-ve NTFROZEN must be same parity  "      &
     &            //"as NEL",.true.)
              end if
          case("ORDER")
              I = 1
              do while ( item .lt. nitems )
                call readf(ORBORDER2(I))
                I = I + 1
              end do
              DO I=1,8
! two ways of specifying open orbitals
! if orborder2(I,1) is integral, then if it's odd, we have a single
! open orbital
                 IF(ORBORDER2(I).EQ.INT(ORBORDER2(I))) THEN
                    ORBORDER(I,1)=IAND(INT(ORBORDER2(I)),65534)
                    IF((INT(ORBORDER2(I))-ORBORDER(I,1)).GT.0) THEN
! we have an open orbital
                       ORBORDER(I,2)=2
                    ELSE
                       ORBORDER(I,2)=0
                    ENDIF
                 ELSE
! non-integral.  The integral part is the number of closed oribtals,
! and the fractional*1000 is the number of open orbitals.
! e.g. 6.002 would mean 6 closed and 2 open
! which would have orborder(I,1)=6, orborder(I,2)=4
! but say 5.002 would be meaningless as the integral part must be a
! multiple of 2
                    ORBORDER(I,1)=INT(ORBORDER2(I)+0.000001)
                    ORBORDER(I,2)=INT((ORBORDER2(I)-ORBORDER(I,1)+      &
     &                             0.000001)*1000)*2
                 ENDIF
              ENDDO
          case("UMATCACHE")
              call readu(w)
              select case(w)
                 case("SLOTS")
                    call geti(NSLOTSINIT)
                 case("MB")
                    call geti(NMEMINIT)
                    nSlotsInit=1
                 case("READ")
                     tReadInCache=.true.
                 case("DUMP")
                     if (iDumpCacheFlag.eq.0) iDumpCacheFlag=1
                 case("FORCE")
                     iDumpCacheFlag=2
                 case default
                    read(w,"(I)") NSLOTSINIT
              end select
          case("NOUMATCACHE")
              NSLOTSINIT=-1
          case ("DFMETHOD")
              call readu(w)
              select case(w)
              case ("DFOVERLAP")
                 iDFMethod=1
              case ("DFOVERLAP2NDORD")
                 iDFMethod=2
              case ("DFOVERLAP2")
                 iDFMethod=3
              case ("DFCOULOMB")
                 iDFMethod=4
             case default
                 call report("keyword "                                 &
     &      //trim(w)//" not recognized in DFMETHOD block",.true.)
              end select
          case("ENDINT")
               exit integral
          case default
              call report("keyword "                                    &
     &      //trim(w)//" not recognized in integral block",.true.)
          end select
        end do integral
      
        END SUBROUTINE 

      END MODULE INTREAD
