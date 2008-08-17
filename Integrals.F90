#include "macros.h"
MODULE Integrals
        Use MemoryManager, only: LogMemAlloc, LogMemDealloc
        USE HElem
        USE System , only : NEL,defaults,Feb08,TUSEBRILLOUIN,tStarStore,OrbOrder,NMSH,BasisFN
        use UMatCache, only: tReadInCache,nSlotsInit,nMemInit,iDumpCacheFlag,iDFMethod
        IMPLICIT NONE
        save

        LOGICAL TQUADRHO,TEXPRHO,THFBASIS,THFCALC,TCALCREALPROD
        LOGICAL TRHF,TReadTUMat,TReadHF,TQuadValMax,TQuadVecMax
        LOGICAL TSUMPROD,TCALCRHOPROD,TDISCONODES,TCalcExcitStar
        LOGICAL TJustQuads,TNoDoubs,TDiagStarStars,TExcitStarsRootChange
        LOGICAL TRmRootExcitStarsRootChange,TLinRootChange
        
        INTEGER NTAY(2),nHFit,NFROZEN,NTFROZEN
        INTEGER NRSTEPSMAX,IHFMETHOD
        
        REAL*8 NRCONV,RFCONV,OrbOrder2(8)
        REAL*8 HFMix,HFEDelta,HFCDelta
        REAL*8 HFRand
        REAL*8 DMatEpsilon !  The cutoff for density matrix elements
        Logical tPostFreezeHF ! Do we do HF after freezing

! // Transferred from system
!        LOGICAL TSTARSTORE

!  From NECI.F
        TYPE(HElement), pointer :: UMAT(:)      
        INTEGER tagUMat
        COMPLEX*16,pointer :: FCK(:)
        INTEGER tagFCK
        INTEGER NMAX
        REAL*8 CST

! from Calc      
        real*8 ChemPot
        Logical tNeedsVirts  ! Set if we need virtual orbitals  (usually set).  Will be unset (by Calc readinput) if I_VMAX=1 and TENERGY is false

        contains

        SUBROUTINE IntReadInput()
        USE input
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
        INTEGER :: i
           
! Integral defaults
      TLinRootChange=.false.
      TRmRootExcitStarsRootChange=.false.
      TExcitStarsRootChange=.false.
      TDiagStarStars=.false.    
      TJustQuads=.false.
      TNoDoubs=.false.
      TCalcExcitStar=.false.
      TQuadVecMax=.false.
      TQuadValMax=.false.
      TDISCONODES=.FALSE.
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
      HFRand=0.01
      DMatEpsilon=0
      tPostFreezeHF=.false.

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
          case("LINROOTCHANGE")
              TLinRootChange=.true.
          case("RMROOTEXCITSTARSROOTCHANGE")
              TRmRootExcitStarsRootChange=.true.
          case("EXCITSTARSROOTCHANGE")
              TExcitStarsRootChange=.true.
          case("DIAGSTARSTARS")
              TDiagStarStars=.true.
          case("STARQUADEXCITS")
              TJustQuads=.true.
          case("STARNODOUBS")
              TNoDoubs=.true.
          case("CALCEXCITSTAR")
              TCalcExcitStar=.true.
          case("QUADVECMAX")
              TQuadVecMax=.true.
          case("QUADVALMAX")
              TQuadValMax=.true.
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
          case("CALCRHOPROD")
              TCALCRHOPROD=.TRUE.
          case("SUMPRODII")
              TSUMPROD=.TRUE.
          case("DISCONNECTNODES")
              TDISCONODES=.TRUE.
          case("HF")
              THFBASIS = .true.
          case("CALCULATE")
              THFCALC = .true.
          case("MAXITERATIONS")
              call geti(NHFIT)
          case("MIX")
              call getf(HFMIX)
          case("RAND")
              call getf(HFRAND)
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
                    call reread(-1)
                    call geti(NSLOTSINIT)
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
          case("POSTFREEZEHF")
            tPostFreezeHF=.true.
          case("DMATEPSILON")
            call readf(DMatEpsilon)
          case("ENDINT")
               exit integral
          case default
              call report("keyword "                                    &
     &      //trim(w)//" not recognized in integral block",.true.)
          end select
        end do integral
      
        END SUBROUTINE IntReadInput


        Subroutine IntInit(iCacheFlag)
!who knows what for
            Use MemoryManager, only: LogMemAlloc, LogMemDealloc
            Use OneEInts, only: SetupTMat
            USE UMatCache, only : FreezeTransfer, CreateInvBRR, GetUMatSize, SetupUMat2D_df
            Use UMatCache, only: InitStarStoreUMat
            Use System, only : nBasisMax, Alpha,BHub, BRR,nmsh
            Use System, only : Ecore,G1,iSpinSkip,nBasis,nMax,nMaxZ
            Use System, only: Omega,tAlpha,TBIN,tCPMD,tDFread,THFORDER
            Use System, only: thub,tpbc,treadint,ttilt,TUEG,tVASP
            Use System, only: uhub, arr,alat,treal
            INTEGER iCacheFlag
            COMPLEX*16 ZIA(*)
            POINTER (IP_ZIA,ZIA)
            INTEGER tagZIA
            COMPLEX*16 COEFF(*)
            POINTER (IP_COEFF,COEFF)
            INTEGER i
            INTEGER TmatInt,UMatInt
            integer iErr
            INCLUDE 'cons.inc'
            character(25), parameter :: this_routine='IntInit'

          FREEZETRANSFER=.false.
                
          IF(THFBASIS) THEN
             WRITE(6,*) "Using Hartree-Fock Basis"
             IF(.NOT.THFCALC) WRITE(6,*) "Reading Hartree-Fock Basis"
          ENDIF
!I_VMAX.EQ.1.AND..NOT.TENERGY.AND.
          IF(.not.tNeedsVirts.and.NTFROZEN.EQ.0) THEN
             WRITE(6,*) "MaxVertices=1 and NTFROZEN=0."
             IF(ABS(ARR(NEL,1)-ARR(NEL+1,1)).GT.1.d-3) THEN
               WRITE(6,*) "NEL spinorbitals completely fills up degenerate set."
               WRITE(6,*) "Only calculating vertex sum, so freezing all virtuals."
               NTFROZEN=NBASIS-NEL
             ELSE
               WRITE(6,*) "NEL spinorbitals does not completely fill up degenerate set."
               WRITE(6,*) "NOT freezing all virtuals."
             ENDIF
          ENDIF

!                IF(.NOT.(TREAD.OR.TREADRHO)) THEN
             IF(NTFROZEN.LT.0) THEN
                I=NEL+(-NTFROZEN)
             ELSE
                I=nBasis-NTFROZEN
             ENDIF
             IF(TCPMD) THEN
    !.. We don't need to do init any 4-index integrals, but we do need to init the 2-index
                WRITE(6,*) " *** INITIALIZING CPMD 2-index integrals ***"
                Allocate(UMat(1), stat=ierr)
                LogAlloc(ierr, 'UMat', 1,HElementSizeB, tagUMat)
                !CALL MEMORY(IP_TMAT,HElementSize*nBasis*nBasis,'TMAT')
                !CALL AZZERO(TMAT,HElementSize*nBasis*nBasis)
                CALL GENSymStatePairs(nBasis/2,.false.)
                CALL SetupTMAT(nBasis,2,TMATINT)
                CALL CPMDINIT2INDINT(nBasis,I,NBASISMAX,ISPINSKIP,G1,NEL,ECORE,THFORDER,ARR,BRR,iCacheFlag)
             ELSEIF(tVASP) THEN
                Allocate(UMat(1), stat=ierr)
                LogAlloc(ierr, 'UMat', 1,HElementSizeB, tagUMat)
                CALL GENSymStatePairs(nBasis/2,.false.)
                CALL SetupTMAT(nBasis,2,TMATINT)
                CALL VASPInitIntegrals(I,ECore,tHFOrder)
             ELSEIF(TREADINT.AND.TDFREAD) THEN
                Allocate(UMat(1), stat=ierr)
                LogAlloc(ierr, 'UMat', 1,HElementSizeB, tagUMat)
                !CALL MEMORY(IP_TMAT,HElementSize*nBasis*nBasis,'TMAT')
                !CALL AZZERO(TMAT,HElementSize*nBasis*nBasis)
                CALL SetupTMAT(nBasis,2,TMATINT)
                Call ReadDalton1EIntegrals(G1,nBasis,Arr,Brr,ECore)
                Call ReadDF2EIntegrals(nBasis,I)
             ELSEIF(TREADINT.AND.TSTARSTORE) THEN
                WRITE(6,*) ' *** READING DOUBLES 2-VERTEX INTEGRALS FROM FCIDUMP *** '
                NBASISMAX(2,3)=2
                !NBASIS/nBasis is number of spin-orbitals
                CALL GENSymStatePairs(nBasis/2,.false.)
                !CALL WRITESYMORBS(nBasis,G1)
                !NBASISMAX(5,2)+1 is the number of irreps
                CALL SetupTMAT(nBasis,2,TMATINT)
                !CALL MEMORY(IP_TMAT,HElementSize*TMATINT,'TMAT')
                !CALL AZZERO(TMAT,HElementSize*TMATINT)
                CALL CREATEINVBRR(BRR,nBasis)
                Call InitStarStoreUMat(nEl/2, nBasis/2)
                CALL GetUMatSize(nBasis,nEl,2,UMATINT)
                Allocate(UMat(UMatInt), stat=ierr)
                LogAlloc(ierr, 'UMat', UMatInt,HElementSizeB, tagUMat)
                Call AZZERO(UMat,HElementSize*UMatInt)
                CALL SETUPUMAT2D_DF()
                IF(TBIN) THEN
                    CALL READFCIINTBIN(UMAT,NBASIS,ECORE,ARR,BRR,G1)
                ELSE
                    CALL READFCIINT(UMAT,NBASIS,ECORE,ARR,BRR,G1)
                ENDIF
                WRITE(6,*) ' ECORE=',ECORE
                ISPINSKIP=2
             ELSEIF(TREADINT) THEN
                WRITE(6,*) ' *** READING PRIMITIVE INTEGRALS FROM FCIDUMP *** '
    !.. Generate the 2e integrals (UMAT)
                ISPINSKIP=NBasisMax(2,3)
                IF(ISPINSKIP.eq.0) STOP 'NBASISMAX(2,3) ISpinSkip unset'
    !!C.. We can only do restricted systems, so alpha and beta spinorbital
    !!C... integrals will be the same 
    !nBasisMax(2,3) is iSpinSkip = 1 if UHF and 2 if RHF
                CALL GetUMatSize(nBasis,nEl,iSpinSkip,UMATINT)
                WRITE(6,*) "UMatSize: ",UMATINT
                Allocate(UMat(UMatInt), stat=ierr)
                LogAlloc(ierr, 'UMat', UMatInt,HElementSizeB, tagUMat)
                Call AZZERO(UMat,HElementSize*UMatInt)
    !nBasisMax(2,3) is iSpinSkip = 1 if UHF and 2 if RHF
                CALL SetupTMAT(nBasis,iSpinSkip,TMATINT)
    !            CALL MEMORY(IP_TMAT,HElementSize*nBasis*nBasis,'TMAT')
    !            CALL AZZERO(TMAT,HElementSize*nBasis*nBasis)
                IF(TBIN) THEN
                    CALL READFCIINTBIN(UMAT,NBASIS,ECORE,ARR,BRR,G1)
                ELSE
                    CALL READFCIINT(UMAT,NBASIS,ECORE,ARR,BRR,G1)
                ENDIF
                WRITE(6,*) ' ECORE=',ECORE
             ELSE
               ISPINSKIP=NBASISMAX(2,3)
               IF(NBASISMAX(1,3).GE.0) THEN
                   IF(TUEG.OR.THUB) THEN
                      IF(THUB.AND.TREAL) THEN
    !!C.. Real space hubbard
    !!C.. we pre-compute the 2-e integrals
                         WRITE(6,*) "Generating 2e integrals"
    !!C.. Generate the 2e integrals (UMAT)
                         CALL GetUMatSize(nBasis,nEl,iSpinSkip,UMATINT)
                         Allocate(UMat(UMatInt), stat=ierr)
                         LogAlloc(ierr, 'UMat', UMatInt,HElementSizeB, tagUMat)
                         Call AZZERO(UMat,HElementSize*UMatInt)
                         CALL CALCUMATHUBREAL(NEL,NBASIS,NBASISMAX,G1,UHUB,UMAT)
                      ELSEIF(THUB.AND..NOT.TPBC) THEN
    !!C.. we pre-compute the 2-e integrals
                         WRITE(6,*) "Generating 2e integrals"
    !!C.. Generate the 2e integrals (UMAT)
                         CALL GetUMatSize(nBasis,nEl,iSpinSkip,UMATINT)
                         Allocate(UMat(UMatInt), stat=ierr)
                         LogAlloc(ierr, 'UMat', UMatInt,HElementSizeB, tagUMat)
                         Call AZZERO(UMat,HElementSize*UMatInt)
    !!C.. Non-periodic hubbard (mom space)
                         CALL GEN_COUL_HUBNPBC(NEL,NBASISMAX,nBasis,G1,NMSH,NMAX,FCK,UMAT,ISPINSKIP,THUB,UHUB,OMEGA)
                      ELSE
    !!C.. Most normal Hubbards
                         IF(.NOT.TUEG) THEN
    !                        CALL GEN_COUL_UEG(NEL,NBASISMAX,nBasis,G1,
    !     &          NMSH,NMAX,FCK,UMAT,ISPINSKIP,THUB,UHUB,OMEGA,ALAT)
                            ISPINSKIP=-1
                            NBASISMAX(2,3)=-1
                         WRITE(6,*) "Not precomputing HUBBARD 2-e integrals"
                            Allocate(UMat(1), stat=ierr)
                            LogAlloc(ierr, 'UMat', 1,HElementSizeB, tagUMat)
                            UMAT(1)=UHUB/OMEGA
                         ENDIF
    !!C.. The UEG doesn't store coul integrals
                      ENDIF
                   ELSE
    !!C.. We need to init the arrays regardless of whether we're storing H
    !!C..Need to initialise the Fourier arrays
                     Allocate(Fck(nMsh**3),stat=ierr)
                     LogAlloc(ierr,'FCK',NMSH**3,16,tagFCK)
                     CALL MEMORY(IP_COEFF,2*(3*NMSH+48),'COEFF')
                     CALL MEMORY(IP_ZIA,2*(NMSH+1)*NMAX*NMAX,'ZIA')
                     WRITE(6,*) NMSH,NMAX
                     LogAlloc(ierr,'ZIA',2*(NMSH+1)*NMAX*NMAX,16,tagZIA)
    !!C..
                     CALL MEMORY_CHECK()
                     IF(NMAXZ.EQ.0) THEN
    !!C..   We're doing a 2D simulation
                        CALL INITFOU2D(NMSH,FCK,COEFF,NMAX,ALAT,TALPHA,ALPHA,OMEGA,ZIA)
                     ELSE
                        CALL INITFOU(NMSH,FCK,COEFF,NMAX,ALAT,TALPHA,ALPHA,OMEGA,ZIA)
                     ENDIF
                     CALL MEMORY_CHECK()
    !!C.. we pre-compute the 2-e integrals
                     WRITE(6,*) "Generating 2e integrals"
    !!C.. Generate the 2e integrals (UMAT)
                     CALL GetUMatSize(nBasis,nEl,iSpinSkip,UMATINT)
                     Allocate(UMat(UMatInt), stat=ierr)
                     LogAlloc(ierr, 'UMat', UMatInt,HElementSizeB, tagUMat)
                     Call AZZERO(UMat,HElementSize*UMatInt)
                     CALL GEN_COUL(NEL,NBASISMAX,nBasis,G1,NMSH,NMAX,FCK,UMAT,ISPINSKIP,ZIA)
                     LogDealloc(tagZIA)
                     Call FreeM(IP_ZIA)
                     LogDealloc(tagFCK)
                     Deallocate(FCK)
                   ENDIF
                ELSE
                   WRITE(6,*) "Not precomputing 2-e integrals"
                   Allocate(UMat(1), stat=ierr)
                   LogAlloc(ierr, 'UMat', 1,HElementSizeB, tagUMat)
                ENDIF
                CALL MEMORY_CHECK()
    !!C.. we need to generate TMAT - Now setup in individual routines
                !CALL MEMORY(IP_TMAT,HElementSize*nBasis*nBasis,'TMAT')
                !CALL AZZERO(TMAT,HElementSize*nBasis*nBasis)
                IF(THUB) THEN
                 CALL CALCTMATHUB(NBASIS,NBASISMAX,BHUB,TTILT,G1,TREAL,TPBC)
                ELSE
    !!C..Cube multiplier
                   CST=PI*PI/(2.D0*ALAT(1)*ALAT(1))
    !!C.. the UEG has k=2pi n/L rather than pi n/L, so we need to multiply the
    !!C.. KE up by 4
                   IF(NBASISMAX(1,1).LE.0) CST=CST*4
                   CALL CALCTMATUEG(NBASIS,ALAT,G1,CST,NBASISMAX(1,1).LE.0,OMEGA)
                ENDIF
             ENDIF
            !ENDIF
        End Subroutine IntInit
        
        
        Subroutine IntFreeze
            Use System, only: Alat,Brr,CoulDampOrb,ECore,fCoulDampMu
            Use System, only: G1,iSpinSkip
            use System, only: nBasis,nEl,arr,nbasismax
            use UMatCache, only: GetUMatSize
            character(25), parameter ::this_routine='IntFreeze'            
!//Locals
            TYPE(HElement), pointer :: UMAT2(:)
            INTEGER tagUMat2, ierr
            INTEGER nOcc
            integer UMATInt
            integer nHG
            nHG=nBasis
            
         CHEMPOT=(ARR(NEL,1)+ARR(NEL+1,1))/2.D0
         WRITE(6,*) "Chemical Potential: ",CHEMPOT
         IF(NTFROZEN.LT.0) THEN
            WRITE(6,*) "NTFROZEN<0.  Leaving ", -NTFROZEN," unfrozen virtuals."
            NTFROZEN=NTFROZEN+nBasis-NEL
         ENDIF
         IF(NFROZEN.GT.0.OR.NTFROZEN.GT.0) THEN
!!C.. At this point, we transform the UMAT and TMAT into a new UMAT and
!!C.. TMAT and Ecore with the frozen orbitals factored in
!!C..
!!C.. a,b are frozen spinorbitals
!!C.. E'core = Ecore+sum_a t_aa + sum_(a<b) (<ab|ab>-<ab|ba>)
!!C.. t'_ii = t_ii+ sum_a ( <ai|ai> - <ai|ia> ) 
!!C.. NHG contains the old number of orbitals
!!C.. NBASIS contains the new
            NBASIS=NBASIS-NFROZEN-NTFROZEN
!!C.. We need to transform some integrals
            !CALL MEMORY(IP_TMAT2,HElementSize*(NBASIS)**2,'TMAT2')
            !CALL AZZERO(TMAT2,HElementSize*(NBASIS)**2)
            IF(NBASISMAX(1,3).GE.0.AND.ISPINSKIP.NE.0) THEN
               CALL GetUMatSize(nBasis,(nEl-NFROZEN),iSpinSkip,UMATINT)
                   Allocate(UMat2(UMatInt), stat=ierr)
                   LogAlloc(ierr, 'UMat2', UMatInt,HElementSizeB, tagUMat)
               CALL AZZERO(UMAT2,HElementSize*UMATINT)
            ELSE
!!C.. we don't precompute 4-e integrals, so don't allocate a large UMAT
                   Allocate(UMat2(1), stat=ierr)
                   LogAlloc(ierr, 'UMat2', 1,HElementSizeB, tagUMat)
            ENDIF 
            CALL MEMORY_CHECK()
            CALL IntFREEZEBASIS(NHG,NBASIS,UMAT,UMAT2,ECORE, G1,NBASISMAX,ISPINSKIP,BRR,NFROZEN,NTFROZEN,NEL,ALAT)
            CALL MEMORY_CHECK()
            WRITE(6,*) "Freezing ",NFROZEN," core orbitals."
            WRITE(6,*) "Freezing ",NTFROZEN," virtual orbitals."
            WRITE(6,*) "ECORE now",ECORE
            WRITE(6,*) NBASIS," orbitals remain."
            NEL=NEL-NFROZEN
            NOCC=NEL/2
!!C.. NEL now only includes active electrons
            WRITE(6,*) "Number of active electrons:",NEL
            
            !CALL FREEM(IP_TMAT)
            !IP_TMAT=IP_TMAT2
            !IP_TMAT2=NULL
!!C.. Now we can remove the old UMATRIX, and set the pointer UMAT to point
!!C.. to UMAT2
            LogDealloc(tagUMat)
            Deallocate(UMat)
            UMat=>UMat2
            nullify(UMat2)
            tagUMat=tagUMat2
            tagUMat2=0
            CALL MEMORY_CHECK()
            WRITE(6,*) "Active basis functions:",NHG
            CALL WRITEBASIS(6,G1,NHG,ARR,BRR)
         ENDIF
!         CALL WRITETMAT(NBASIS)
!         CALL WRITESYMCLASSES(NBASIS)
         
         IF(COULDAMPORB.GT.0) THEN
            FCOULDAMPMU=(ARR(COULDAMPORB,1)+ARR(COULDAMPORB+1,1))/2
            WRITE(6,*) "Setting Coulomb damping mu between orbitals ",ARR(COULDAMPORB,1)," and ",ARR(COULDAMPORB+1,1)
            WRITE(6,*) "MU=",FCOULDAMPMU
         ENDIF
        End Subroutine
        Subroutine IntCleanup(iCacheFlag)
            Use System, only: G1, nBasis
            Use UMatCache, only: iDumpCacheFlag,tReadInCache,nStates, nStatesDump, DumpUMatCache, DestroyUMatCache, WriteUMatCacheStats
            INTEGER iCacheFlag
            if ((btest(iDumpCacheFlag,0).and.(nStatesDump.lt.nStates.or..not.tReadInCache)).or.btest(iDumpCacheFlag,1)) then
            call DumpUMatCache(nBasis,G1)
            end if
!If we're told explicitly not to destroy the cache, we don't
              if(.NOT.BTEST(iCacheFlag,0)) THEN
                 CALL DESTROYUMATCACHE()
              ELSE
                 Call WriteUMatCacheStats()
              ENDIF
        END Subroutine IntCleanup
      SUBROUTINE IntFREEZEBASIS(NHG,NBASIS,UMAT,UMAT2,ECORE,           &
     &         G1,NBASISMAX,ISS,BRR,NFROZEN,NTFROZEN,NEL,ALAT)
         USE HElem
         use System, only: Symmetry,BasisFN,BasisFNSize,arr,tagarr
         use OneEInts
         USE UMatCache, only: FreezeTransfer,UMatCacheData,UMatInd,TUMat2D
         Use UMatCache, only: FreezeUMatCache, CreateInvBrr2,FreezeUMat2D, SetupUMatTransTable
         use UMatCache, only: GTID
         IMPLICIT NONE
         INTEGER NHG,NBASIS,nBasisMax(5,*),ISS
         TYPE(BASISFN) G1(NHG)
         TYPE(HElement) UMAT(*)
!!C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         TYPE(HElement) UMAT2(*)
         REAL*8 ECORE
!!C.. was (NBASIS/ISS,NBASIS/ISS,NBASIS/ISS,NBASIS/ISS)
         REAL*8 ARR2(NBASIS,2)
         INTEGER NFROZEN,BRR(NHG),BRR2(NBASIS),GG(NHG)
         TYPE(BASISFN) G2(NHG)
         INTEGER NTFROZEN
         INTEGER I,J,K,L,A,B,IP,JP,KP,LP,IDI,IDJ,IDK,IDL
         INTEGER IB,JB,KB,LB,AB,BB,IPB,JPB,KPB,LPB
         INTEGER IDIP,IDJP,IDKP,IDLP
         INTEGER IDA,IDB,SGN
         INTEGER iSize,ierr
         INTEGER FDET(NEL),NEL
         TYPE(Symmetry) KSYM
         REAL*8 ALAT(3)
         character(*), parameter :: this_routine='IntFreezeBasis'

!!C.. Just check to see if we're not in the middle of a degenerate set with the same sym
         IF(NFROZEN.GT.0) THEN
            IF(ABS(ARR(NFROZEN,1)-ARR(NFROZEN+1,1)).LT.1.D-6.AND.       &
     &        G1(BRR(NFROZEN))%SYM%s.EQ.G1(BRR(NFROZEN+1))%SYM%s) THEN
               STOP "Cannot freeze in the middle of a degenerate set"
            ELSE IF (ABS(ARR(NFROZEN,1)-ARR(NFROZEN+1,1)).LT.1.D-6) THEN
               write (6,'(a)') 'WARNING: Freezing in the middle of a degenerate set.'
               write (6,'(a)') 'This should only be done for debugging purposes.'
            ENDIF
         ENDIF
         IF(NTFROZEN.GT.0) THEN
            IF(ABS(ARR(NHG-NTFROZEN,1)-ARR(NHG-NTFROZEN+1,1)).LT.1.D-6  &
     &         .AND.G1(BRR(NHG-NTFROZEN))%SYM%s                         &
     &               .EQ.G1(BRR(NHG-NTFROZEN+1))%SYM%s) THEN
          STOP "Cannot freeze in the middle of a degenerate virtual set"
            ELSE IF (ABS(ARR(NHG-NTFROZEN,1)-ARR(NHG-NTFROZEN+1,1)).LT.1.D-6) THEN
               write (6,'(a)') 'WARNING: Freezing in the middle of a degenerate set.'
               write (6,'(a)') 'This should only be done for debugging purposes.'
            ENDIF
         ENDIF

!!C.. At this point, we transform the UMAT and TMAT into a new UMAT and
!!C.. TMAT and Ecore with the frozen orbitals factored in
!!C..
!!C.. a,b are frozen spinorbitals
!!C.. E'core = Ecore+sum_a t_aa + sum_(a<b) (<ab|ab>-<ab|ba>)
!!C.. t'_ii = t_ii+ sum_a ( <ai|ai> - <ai|ia> ) 
!!C.. NHG contains the old number of orbitals
!!C.. NBASIS contains the new
!!C.. We first need to work out where each of the current orbitals will
!!C.. end up in the new set
         K=0
         DO I=1,NHG
            L=1
            DO J=1,NFROZEN
!C.. if orb I is to be frozen, L will become 0
               IF(BRR(J).EQ.I) L=0
            ENDDO
            DO J=NBASIS+NFROZEN+1,NHG
!C.. if orb I is to be frozen, L will become 0
               IF(BRR(J).EQ.I) L=0
            ENDDO
            IF(L.EQ.0) THEN
               GG(I)=0
            ELSE
!C.. we've got an orb which is not to be frozen 
               K=K+L
!C.. GG(I) is the new position in G of the (old) orb I
               GG(I)=K
!C.. copy the eigenvalue table to the new location
               CALL ICOPY(BasisFNSize,G1(I),1,G2(K),1)
            ENDIF
         ENDDO
!C.. Now construct the new BRR and ARR
         DO I=1,NBASIS
            BRR2(I)=GG(BRR(I+NFROZEN))
            ARR2(I,1)=ARR(I+NFROZEN,1)
         ENDDO
         DO I=1,NHG
            IF(GG(I).NE.0) ARR2(GG(I),2)=ARR(I,2)
         ENDDO

        IF(TSTARSTORE.or.tCPMDSymTMat) THEN
!C.. Now setup the default symmetry to include the frozen electrons
!C.. BRR(1:NFROZEN) is effectively the det of the frozens, so we get its sym
            CALL GETSYM(BRR,NFROZEN,G1,NBASISMAX,KSym)
            CALL SetupFreezeAllSym(KSYM)
!C.. Freeze the sym labels
            !SYMCLASSES2 gives the new symmetry of the frozen set of orbitals
            CALL FREEZESYMLABELS(NHG,NBASIS,GG,.true.)
!C.. Copy the new G1 over the old ones
            CALL ICOPY(BasisFNSize*NBASIS,G2,1,G1,1)
            !Redo SYMLABELCOUNTS
            CALL GENSymStatePairs(NBASIS/2,.true.)
         ENDIF
         CALL SetupTMAT2(NBASIS,2,iSize)

!C.. First deal with Ecore
         DO A=1,NFROZEN
            AB=BRR(A)
            ECORE=ECORE+DREAL(GetTMATEl(AB,AB))
            DO B=A+1,NFROZEN
               BB=BRR(B)
               CALL GTID(NBASISMAX,AB,IDA)
               CALL GTID(NBASISMAX,BB,IDB)
!C.. No sign problems from permuations here as all perms even
               ECORE=ECORE+DREAL(GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDB,IDA,IDB))
!C.. If we have spin-independent integrals, or 
!C.. if the spins are the same
               IF(G1(AB)%MS.EQ.G1(BB)%MS)                               &
     &            ECORE=ECORE-DREAL(GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDB,IDB,IDA))
            ENDDO
         ENDDO
!C.. now deal with the new TMAT
         FREEZETRANSFER=.true.
         DO I=1,NBASIS
            IP=I+NFROZEN
            IB=BRR(IP)
            IPB=GG(IB)
            CALL GTID(NBASISMAX,IB,IDI)
            DO J=1,NBASIS
               JP=J+NFROZEN
               JB=BRR(JP)
               JPB=GG(JB)
               CALL GTID(NBASISMAX,JB,IDJ)
               IF(TSTARSTORE.or.tCPMDSymTMat) THEN
                   TMATSYM2(NEWTMATInd(IPB,JPB))=GetTMATEl(IB,JB)
               ELSE
                   TMAT2D2(IPB,JPB)=GetTMATEl(IB,JB)
               ENDIF
               DO A=1,NFROZEN
                  AB=BRR(A)
                  CALL GTID(NBASISMAX,AB,IDA)
!C.. SGN takes into account permutationnness.
!C                  SGN=1
!C                  IF(IB.GT.AB) SGN=-SGN
!C                  IF(JB.GT.AB) SGN=-SGN
                  IF(G1(IB)%MS.EQ.G1(JB)%MS) THEN
                      IF(TSTARSTORE.or.tCPMDSymTMat) THEN
                        TMATSYM2(NEWTMATInd(IPB,JPB))=                  &
     &                  GetNEWTMATEl(IPB,JPB)+GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDI,IDA,IDJ)
                      ELSE
                         TMAT2D2(IPB,JPB)=TMAT2D2(IPB,JPB)+GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDI,IDA,IDJ)
                      ENDIF
                  ENDIF
!C.. If we have spin-independent integrals, ISS.EQ.2.OR
!C.. if the spins are the same
                  IF(G1(IB)%MS.EQ.G1(AB)%MS.AND.G1(AB)%MS.EQ.G1(JB)%MS) THEN
                        IF(TSTARSTORE.or.tCPMDSymTMat) THEN
                            TMATSYM2(NEWTMATInd(IPB,JPB))=GetNEWTMATEl(IPB,JPB) &
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDI,IDJ,IDA)        
                        ELSE
                            TMAT2D2(IPB,JPB)=GetNEWTMATEl(IPB,JPB)              &
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDI,IDJ,IDA)
                        ENDIF
                    ENDIF
               ENDDO
!               WRITE(6,*) "T",TMAT(IB,JB),I,J,TMAT2(IPB,JPB)
!            IF(TMAT(IPB,JPB).AGT.1.D-9) WRITE(16,*) I,J,TMAT2(IPB,JPB)
            ENDDO
         ENDDO
         IF(NBASISMAX(1,3).GE.0.AND.ISS.NE.0) THEN

             CALL CREATEINVBRR2(BRR2,NBASIS)
!CC Only do the below if we've a stored UMAT
!C.. Now copy the relevant matrix elements of UMAT across
!C.. the primed (...P) are the new versions
         DO I=1,NBASIS
            IB=BRR(I+NFROZEN)
            IPB=GG(IB)
            IF(ISS.NE.0.OR.G1(I)%MS.EQ.1) THEN
             CALL GTID(NBASISMAX,IB,IDI)
             CALL GTID(NBASISMAX,IPB,IDIP)
             DO J=1,NBASIS
               JB=BRR(J+NFROZEN)
               JPB=GG(JB)
               IF(ISS.NE.0.OR.G1(I)%MS.EQ.1) THEN
                CALL GTID(NBASISMAX,JB,IDJ)
                CALL GTID(NBASISMAX,JPB,IDJP)
                DO K=I,NBASIS
                  KB=BRR(K+NFROZEN)
                  KPB=GG(KB)
                  IF(ISS.NE.0.OR.G1(I)%MS.EQ.1) THEN
                   CALL GTID(NBASISMAX,KB,IDK)
                   CALL GTID(NBASISMAX,KPB,IDKP)
                   DO L=J,NBASIS
                     IF((K*(K-1))/2+I.GE.(L*(L-1))/2+J) THEN
                        LB=BRR(L+NFROZEN)
                        LPB=GG(LB)
                        IF(ISS.NE.0.OR.G1(I)%MS.EQ.1) THEN
                         CALL GTID(NBASISMAX,LB,IDL)
                         CALL GTID(NBASISMAX,LPB,IDLP)

                         IF(TSTARSTORE) THEN
                             IF(.NOT.TUMAT2D) STOP 'UMAT2D should be on'
                             IF((IDI.eq.IDJ.and.IDI.eq.IDK.and.IDI.eq.IDL).or.  &
     &                          (IDI.eq.IDK.and.IDJ.eq.IDL).or.                 &
     &                          (IDI.eq.IDL.and.IDJ.eq.IDK).or.                 &
     &                          (IDI.eq.IDJ.and.IDK.eq.IDL)) THEN
                                CONTINUE
                             ELSE
                                 IF(((I+NFROZEN).gt.NEL).or.((J+NFROZEN).gt.NEL)) THEN
                                     CONTINUE
                                 ELSE
                                    UMAT2(UMatInd(IDIP,IDJP,IDKP,IDLP,0,(NEL-NFROZEN)/2)) &
     &                          =UMAT(UMatInd(IDI,IDJ,IDK,IDL,NHG/2,0))
                                ENDIF
                             ENDIF
                         ELSE
                         
                            UMAT2(UMatInd(IDIP,IDJP,IDKP,IDLP,0,0))             &
     &                        =UMAT(UMatInd(IDI,IDJ,IDK,IDL,NHG/2,0))
                         ENDIF
                        ENDIF
                     ENDIF
                   ENDDO
                  ENDIF
                ENDDO
               ENDIF
             ENDDO
            ENDIF
         ENDDO
         CALL FLUSH(11)
         CALL FLUSH(12)
            IF(TSTARSTORE) CALL FreezeUMAT2D(NHG,NBASIS,GG,2)
         ELSEIF(Associated(UMatCacheData)) THEN
!.. We've a UMAT2D and a UMATCACHE.  Go and Freeze them
!C.. NHG contains the old number of orbitals
!C.. NBASIS contains the new
!C.. GG(I) is the new position in G of the (old) orb I
            CALL FreezeUMatCache(GG,NHG,NBASIS)
         ENDIF
         IF(ISS.EQ.0) CALL SetupUMatTransTable(GG,nHG,nBasis)
         IF(.NOT.TSTARSTORE) THEN
             CALL GETSYM(BRR,NFROZEN,G1,NBASISMAX,KSym)
             CALL SetupFREEZEALLSYM(KSym)
             CALL FREEZESYMLABELS(NHG,NBASIS,GG,.false.)
             CALL ICOPY(BasisFNSize*NBASIS,G2,1,G1,1)
         ENDIF 
         FREEZETRANSFER=.false.
!C.. Copy the new BRR and ARR over the old ones
         CALL SWAPTMAT(NBASIS,NHG,GG)
         deallocate(arr)
         LogDealloc(tagarr)
         allocate(arr(nBasis,2),stat=ierr)
         LogAlloc(ierr,'Arr',2*nBasis,8,tagArr)

         CALL ICOPY(NBASIS,BRR2,1,BRR,1)
         CALL DCOPY(NBASIS*2,ARR2,1,ARR,1)
!C.. Now reset the total number of orbitals
         NHG=NBASIS
         Call DetFreezeBasis(GG)
         
         RETURN
      end subroutine intfreezebasis



      function GetUMatEl2(I,J,A,B)
      ! A wrapper for GetUMatEl, now everything is available via modules.
      ! In:
      !    I,J,A,B: indices of integral
      ! Returns <ij|ab>
         use System, only: ALAT,G1,iSpinSkip,nBasis,nBasisMax
         implicit none
         TYPE(HElement) GetUMatEl2
         integer :: I,J,A,B
         
         GetUMatEl2=GetUMatEl(nBasisMax,UMat,ALAT,nBasis,iSpinSkip,G1,I,J,A,B)

      end function GetUMatEl2



      FUNCTION GetUMatEl(NBASISMAX,UMATstore,ALAT,NHG,ISS,G1,IDI,IDJ,IDK,IDL)
      ! Get a U matrix element <ij|u|kl> in multifarious ways, where orbitals
      ! are spatial orbitals.  Either from a passed-in UMAT, or ALAT parameters
      ! (Hubbard, UEG, particle in a box), or from UMatcache (CPMD, DF).
      ! In:
      !    nBasisMax: legacy.  Contains some information that Alex is not clear on (rant at will).
      !    UMatstore: Store of <ij|u|kl> integrals.
      !    ALAT: Size of cell/box for Hubbard/UEG/particle in a box.
      !    iSS: as above in GetUMatSize.
      !    NHG: # basis functions.`
      !    G1: symmetry and momentum information on the basis functions.
      !    IDI,IDJ,IDK,IDL: indices for integral.
         use System, only: Symmetry,BasisFN,tVASP
         use UMatCache
         use vasp_neci_interface, only: CONSTRUCT_IJAB_one
         IMPLICIT NONE
         TYPE(HElement) GetUMatEl
         INTEGER nBasisMax(5,*),I,J,K,L,NHG,ISS
         TYPE(BasisFN) G1(NHG)
         REAL*8 ALAT(3),GetNan
         TYPE(HElement) UMATstore(*)
         TYPE(HElement) UElems(0:nTypes-1)
         complex*16 vasp_int(1,0:1)
         INTEGER A,B,C,XXX
         INTEGER IDI,IDJ,IDK,IDL
         REAL*8 SUM
         real*8, PARAMETER :: PI=3.14159265358979323846264338327950288419716939937510D0
         INTEGER ICACHE,ICACHEI,ITYPE
         LOGICAL LSYMSYM
         TYPE(Symmetry) SYM,SYMPROD,SYMCONJ
         INTEGER ISUB,ISUB2
         Type(Symmetry) TotSymRep
         logical GetCachedUMatEl,calc2ints
!   IF NBASISMAX(1,3) is less than zero, we directly give the integral.
!   Otherwise we just look it up in umat
         IF(NBASISMAX(1,3).GE.0) THEN

!   See if we need to calculate on the fly
            IF(ISS.EQ.0) THEN

!  JSS - store <ij|ij> and <ij|ji> in UMAT2D.
!   Remember permutations.  
!   <ij|ji> = <ii|jj>

!.. Complex case is more difficult.

!  <ii|ii> is always real and allowed
!  <ij|ij> is always real and allowed (densities i*i and j*j)
!  <ij|ji> is always real and allowed (codensities i*j and j*i = (i*j)*)
!  <ii|jj> is not stored in UMAT2D, and may not be allowed by symmetry.  It can be complex.
           IF (IDI.eq.IDJ.and.IDI.eq.IDK.and.IDI.eq.IDL.AND.TUMAT2D) THEN
!    <ii|ii>
             GETUMATEL=UMAT2D(IDI,IDI)
           ELSE IF (IDI.eq.IDK.and.IDJ.eq.IDL.AND.TUMAT2D) THEN
!   <ij|ij>
             I=MIN(IDI,IDJ)
             J=MAX(IDI,IDJ)
             GETUMATEL=UMAT2D(I,J)
           ELSE IF (IDI.eq.IDL.and.IDJ.eq.IDK.AND.TUMAT2D) THEN
!   <ij|ji>
             I=MAX(IDI,IDJ)
             J=MIN(IDI,IDJ)
             GETUMATEL=UMAT2D(I,J)
!          ELSE IF (IDI.eq.IDJ.and.IDK.eq.IDL.AND.TUMAT2D.AND.HElementSize.EQ.1) THEN
!   <ii|jj> = <ij|ji> Only for real systems (and not for the local exchange
!   scheme.)
!            I=MAX(IDI,IDK)
!            J=MIN(IDI,IDK)
!            GETUMATEL=UMAT2D(I,J)
           ELSE
!   Check to see if the umat element is in the cache
               I=IDI
               J=IDJ
               K=IDK
               L=IDL
               SYM=TotSymRep()
               SYM=SYMPROD(SYM,SYMCONJ(G1(I*2-1)%Sym))
               SYM=SYMPROD(SYM,SYMCONJ(G1(J*2-1)%Sym))
               SYM=SYMPROD(SYM,G1(K*2-1)%Sym)
               SYM=SYMPROD(SYM,G1(L*2-1)%Sym)
!   Check the symmetry of the 4-index integrals
              IF(.NOT.LSYMSYM(SYM)) THEN
                  GETUMATEL=0.D0
                  RETURN
              ELSE
               
!First check whether we can reduce a set of k-points to a simpler symmetry related one
              If(HasKPoints()) THEN
                IF(TTRANSFINDX) THEN
                 I=TransTable(I)
                 J=TransTable(J)
                 K=TransTable(K)
                 L=TransTable(L)
                ENDIF
                ! As we're not looping over i,j,k,l, it's safe to return the
                ! k-pnt related labels in the same variables.
                call KPntSymInt(I,J,K,L,I,J,K,L)
                IF(TTRANSFINDX) THEN
                 I=InvTransTable(I)
                 J=InvTransTable(J)
                 K=InvTransTable(K)
                 L=InvTransTable(L)
                ENDIF
               ENDIF
!   This will rearrange I,J,K,L into the correct order
!   (i,k)<=(j,l) and i<=k, j<=l.
               IF(GETCACHEDUMATEL(I,J,K,L,GETUMATEL,ICACHE,ICACHEI,A,B,ITYPE)) THEN
!   We don't have a stored UMAT - we call to generate it.
                  IF(tDFInts) THEN
!   We're using density fitting
                     Call GetDF2EInt(I,J,K,L,UElems)
                     GetUMatEl=UElems(0)
                  ELSE IF (tVASP) then
                     IF(TTRANSFINDX) THEN
                        CALL CONSTRUCT_IJAB_one(TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L),vasp_int(1,0))
                        CALL CONSTRUCT_IJAB_one(TRANSTABLE(I),TRANSTABLE(L),TRANSTABLE(K),TRANSTABLE(J),vasp_int(1,1))
                     ELSE
                        CALL CONSTRUCT_IJAB_one(I,J,K,L,vasp_int(1,0))
                        CALL CONSTRUCT_IJAB_one(I,L,K,J,vasp_int(1,1))
                     END IF
                     UElems(0)=HElement(vasp_int(1,0))
                     UElems(1)=HElement(vasp_int(1,1))
                     GetUMatEl=UElems(0)
!  Bit 0 tells us which integral in the slot we need
                     GETUMATEL=UElems(IAND(ITYPE,1))
!  Bit 1 tells us whether we need to complex conj the integral
                     IF(BTEST(ITYPE,1)) GETUMATEL=DCONJG(GETUMATEL)
                  ELSE
!   Otherwise we call CPMD
!   Only need <IJ|KL> if we're doing a 2-vertex calculation unless the integral
!   is for a single excitation, in which case we need <IL|JK> as well.
                     calc2ints=gen2CPMDInts.or.(IDI.eq.IDJ.or.IDI.eq.IDK.or.IDI.eq.IDL.or.IDJ.eq.IDK.or.IDJ.eq.IDL)
                     IF(TTRANSFINDX) THEN
                        CALL INITFINDXI(TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L),UElems,calc2ints)
                     ELSE
                        CALL INITFINDXI(I,J,K,L,UElems,calc2ints)
!InitFindxI returns up to two integrals in UElems
!  <ij|u|kl> and <kj|u|il> (which are distinct when complex orbitals are used).
!  TYPE 0          TYPE 1

                     ENDIF
!  Bit 0 tells us which integral in the slot we need
                     GETUMATEL=UElems(IAND(ITYPE,1))
!  Bit 1 tells us whether we need to complex conj the integral
                     IF(BTEST(ITYPE,1)) GETUMATEL=DCONJG(GETUMATEL)
                  ENDIF
!  Because we've asked for the integral in the form to be stored, we store as iType=0
                  IF(ICACHE.NE.0) CALL CACHEUMATEL(A,B,UElems,ICACHE,ICACHEI,0)
                  NMISSES=NMISSES+1
               ELSE
                  NHITS=NHITS+1
              ENDIF
             ENDIF
           ENDIF

          ELSEIF(ISS.EQ.-1) THEN

!  A  non-stored hubbard integral.
          CALL GetHubUMatEl(IDI,IDJ,IDK,IDL,UMatstore,nBasisMax,G1,GetUMatEl)

          ELSE
             IF(TSTARSTORE) THEN
                 IF(.not.TUMAT2D) STOP 'UMAT2D should be on'
                 IF(IDI.eq.IDJ.and.IDI.eq.IDK.and.IDI.eq.IDL) THEN
!    <ii|ii>
                     GETUMATEL=UMAT2D(IDI,IDI)
                 ELSEIF (IDI.eq.IDK.and.IDJ.eq.IDL) THEN
!   <ij|ij> - coulomb
                     I=MIN(IDI,IDJ)
                     J=MAX(IDI,IDJ)
                     GETUMATEL=UMAT2D(I,J)
                 ELSEIF (IDI.eq.IDL.and.IDJ.eq.IDK) THEN
!   <ij|ji> - exchange
                     I=MAX(IDI,IDJ)
                     J=MIN(IDI,IDJ)
                     GETUMATEL=UMAT2D(I,J)
                 ELSEIF (IDI.eq.IDJ.and.IDK.eq.IDL) THEN
                     I=MAX(IDI,IDK)
                     J=MIN(IDI,IDK)
                     GETUMATEL=UMAT2D(I,J)
                 ELSE
                     XXX=UMatInd(IDI,IDJ,IDK,IDL,NHG/2,0)
                     IF(XXX.ne.-1) THEN
                         GETUMATEL=UMATstore(XXX)
                     ELSE
                         GETUMATEL=HElement(GETNAN())
                     ENDIF
                 ENDIF
             ELSE
                 GETUMATEL=UMATstore(UMatInd(IDI,IDJ,IDK,IDL,0,0))
             ENDIF
          ENDIF

         ELSEIF(NBASISMAX(1,3).EQ.-1) THEN
            CALL GetUEGUmatEl(IDI,IDJ,IDK,IDL,ISS,G1,ALAT,iPeriodicDampingType,GetUMatEl)
         ENDIF

         RETURN
      END FUNCTION GetUMatEl



      SUBROUTINE WRITESYMCLASSES(NBASIS)
        USE HElem
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        USE UMatCache
        use SymData, only: SymClasses,SymLabelCounts,nSymLabels
        IMPLICIT NONE
        INTEGER I,NBASIS
        
        DO I=1,NBASIS/2
            WRITE(13,*) I,SYMCLASSES(I)
            CALL FLUSH(13)
        ENDDO
        DO I=1,NSYMLABELS
            WRITE(13,*) I,SYMLABELCOUNTS(2,I)
        ENDDO
        WRITE(13,*) "******************"
        CALL FLUSH(13)
      END subroutine writesymclasses

END MODULE Integrals


! Calculate the diagonal matrix elements for the Uniform electron gas
!  CST is PI*PI/2L*L for the non-periodic case, and
!  CST is 4*PI*PI/2L*L for the periodic case, and

!  For the periodic case we must also add in a periodic images correction
!  which is 1/2 (<ii|ii> - <ii|ii>cell) for each orbital
!  We calculate <ii|ii>cell with a potenial v(r)=1/r (r<Rc) and 0 (r>=Rc)
!  Rc=ALAT(4).

      SUBROUTINE CALCTMATUEG(NBASIS,ALAT,G1,CST,TPERIODIC,OMEGA)
         USE HElem
         use System, only: BasisFN
         USE OneEInts, only : SetupTMAT,TMAT2D,TSTARSTORE
         IMPLICIT NONE
         INTEGER NBASIS
         TYPE(BASISFN) G1(NBASIS)
         REAL*8 ALAT(4),HFBASIS(NBASIS,NBASIS),CST
         INTEGER I,J
         INTEGER iSIZE
         REAL*8 SUM,S1,OMEGA
         LOGICAL TPERIODIC
         REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795029D0
         IF(TPERIODIC) WRITE(6,*) "Periodic UEG"
         OPEN(10,FILE='TMAT',STATUS='UNKNOWN')
         IF(TSTARSTORE) STOP 'Cannot use TSTARSTORE with UEG'
         CALL SetupTMAT(NBASIS,2,iSIZE)
          DO I=1,NBASIS
           TMAT2D(I,I)=((ALAT(1)**2)*((G1(I)%K(1)**2)/(ALAT(1)**2)+        &
     &         (G1(I)%K(2)**2)/(ALAT(2)**2)+(G1(I)%K(3)**2)/(ALAT(3)**2)))
           TMAT2D(I,I)=TMAT2D(I,I)*HElement(CST)
!..  The G=0 component is explicitly calculated for the cell interactions as 2 PI Rc**2 .
!     we *1/2 as we attribute only half the interaction to this cell.
           IF(TPERIODIC) TMAT2D(I,I)=TMAT2D(I,I)-HElement(PI*ALAT(4)**2/OMEGA)
            WRITE(10,*) I,I,TMAT2D(I,I)
         ENDDO
         CLOSE(10)
         RETURN
      END

