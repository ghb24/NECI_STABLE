#include "macros.h"
      MODULE System
        USE input
        
        IMPLICIT NONE
        save

        LOGICAL TSTARBIN,TREADINT,THFORDER,TDFREAD,TPBC,TUEG,TCPMD,THUB
        LOGICAL TSPN,TCSF,TPARITY,TUSEBRILLOUIN,TEXCH,TREAL,TTILT
        LOGICAL TALPHA,TSTAR,TSTOREASEXCITATIONS,TBIN,tStarStore
        INTEGER LMS,STOT,IPARITY(5),NMAXX,NMAXY,NMAXZ,NMSH,COULDAMPORB
        INTEGER iPeriodicDampingType,ISTATE,NEL,ITILTX,ITILTY
        REAL*8 BOX,BOA,COA,FUEGRS,fRc,FCOUL,OrbECutoff,UHUB,BHUB
        REAL*8 ALPHA,FCOULDAMPBETA,FCOULDAMPMU
! Defaults stored in this module
        LOGICAL :: defaults,Feb08

! Used to be stored in Integrals
        INTEGER ORBORDER(8,2)

! From NECICB
        integer lmsBasis

        TYPE Symmetry
           SEQUENCE
           INTEGER*8 S
        END TYPE
     
        integer, PARAMETER :: SymmetrySize=2
        integer, PARAMETER :: SymmetrySizeB=SymmetrySize*8
        TYPE BasisFN
           SEQUENCE
           INTEGER k(3)
           INTEGER Ms
           TYPE(Symmetry) sym
        END TYPE
  
        integer, PARAMETER :: BasisFNSize=SymmetrySize+4
        integer, PARAMETER :: BasisFNSizeB=BasisFNSize*8


        TYPE(BASISFN) SymRestrict
        INTEGER NBASISMAX(5,7)
        REAL*8 ALAT(5)
        REAL*8 ECore
        INTEGER nBasis
        integer nMax
        integer nnr
        integer nocc
        REAL*8 OMEGA
        logical tSpinPolar
        INTEGER iSpinSkip

!From Calc  
        REAL*8 Beta
        
!Renewed for compile

        REAL*8, pointer         :: Arr(:,:)        !List of orbital energies.  (:,1) is ordered, (:,2) is the energy of orbital :
        INTEGER tagArr
        INTEGER, pointer        :: BRR(:)          !Lists orbitals in energy order. i.e. Brr(1) is the lowest energy orbital
        INTEGER tagBrr
        Type(BasisFN), pointer  :: G1(:)           !Info about the basis functions
        INTEGER tagG1
        
        INTEGER LMS2
        
        contains

        SUBROUTINE SysReadInput()
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
        INTEGER I
        
!       SYSTEM defaults - leave these as the default defaults
!       Any further addition of defaults should change these after.
          TSTARSTORE=.false.
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

        END SUBROUTINE SysReadInput








        
        Subroutine SysInit
            Use MemoryManager, only: LogMemAlloc, LogMemDealloc
            character(25), parameter :: this_routine='SysInit'
            integer ierr

           CHARACTER CPAR(3)*1,CPARITY*3
! For init of mom
        TYPE(BasisFN) G
        INCLUDE 'csf.inc'
        INCLUDE 'cons.inc'
        
!  For the UEG
        REAL*8 FKF,Rs
        
! General variables
        INTEGER i,j,k,l,iG
        INTEGER len
        INTEGER iSub
        REAL*8 SUM
! Called functions
        type(Symmetry) TotSymRep
        TYPE(BasisFN) FrzSym
        logical kallowed

!      write (6,*)
!      call TimeTag()
!      if (.not.TCPMD) call Envir()
!      write (6,*)

      ECORE=0.D0
      
! //AJWT TBR
!      IFDET=0
!      TRHOIJND=.false.
      CALL TISET('SysInit   ',ISUB)


!C ==-------------------------------------------------------------------==
!C..Input parameters
      WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
      WRITE(6,*) ' TSPN : ' , TSPN
      CALL IAZZERO(NBASISMAX,5*7)
      TSPINPOLAR=.FALSE.
      DO I=1,3
         SymRestrict%k(I)=IPARITY(I)
      ENDDO
      SymRestrict%Ms=IPARITY(4)
      SymRestrict%Sym%s=IPARITY(5)

      if(tUseBrillouin) THEN
         WRITE(6,*) "Using Brillouin's Theorem to ignore single excitations"
      endif
      if(tStoreAsExcitations) THEN
         write(6,*) "Storing determinants as excitations from the HF determinant.  WARNING this may not work!"
         IF(nEL.lt.8) STOP 'tStoreAsExcitations requires nEl>=8.'
      endif
      IF(TSPN) THEN
!C.. If we're doing monte carlo without generating a list of
!C.. determinants, we cannot as yet force spin or parity, except by
!C.. restricting the basis set.  This will only work for Ms=NEL/2
!C         IF((.NOT.TBLOCK).AND.(.NOT.TCALCHMAT.OR.NTAY.LE.0)) THEN
!C            WRITE(6,*) 'TSPN set to TRUE.  Determinant list not being',
!C     &         ' used for MC.  Forcing MS=Nel/2'
!C            IF(MOD(NEL,2).EQ.0) THEN
!C               LMS=NEL/2
!C            ELSE
!C               LMS=NEL
!C            ENDIF
!C            TSPINPOLAR=.TRUE.
!C         ENDIF
         IF(MOD(LMS+NEL*2,2).NE.MOD(NEL,2)) THEN
           WRITE(6,*) 'LMS=',LMS,' not achievable with',NEL,' electrons'
           WRITE(6,*) 'Resetting LMS'
           LMS=MOD(NEL,2)
         ENDIF
         LMS2=LMS
      ENDIF
      WRITE(6,*) ' GLOBAL MS : ' , LMS
      IF(TCSF) THEN
         WRITE(6,*) "Using CSFs."
         IF(TSPN) THEN
            WRITE(6,*) "Restricting total spin*2 to ",STOT
            IF(LMS.GT.STOT) STOP "Cannot have LMS>STOT"
!C.. Encode the symmetry for the total spin in LMS
            LMS2=LMS+STOT*CSF_NBSTART
         ENDIF
         NBASISMAX(4,7)=1
      ENDIF
      

      IF(TCPMD) THEN
         WRITE(6,*) ' *** GENERIC SYSTEM USING KOHN-SHAM ORBITALS *** '
         CALL CPMDSYSTEMINIT(LEN)   
         IF(TPARITY) THEN
            WRITE(6,"(A,$)") ' SYMMETRIES : '
            CALL WRITEALLSYM(5,SymRestrict)
         ENDIF
         IF(THFORDER) WRITE(6,*)      "Ordering according to 1-electron energies."
      ELSEIF(TREADINT) THEN
!C.. we read in the integrals from FCIDUMP and ignore most config
!C..   
         WRITE(6,*) ' *** GENERIC SYSTEM *** '
         IF(THUB) THEN
            THUB=.FALSE.
            WRITE(6,*) "Setting THUB=.FALSE."
         ENDIF
         IF(TDFREAD) THEN
            WRITE(6,*) "Reading Density fitted integrals."
            LMSBASIS=LMS
            CALL InitDFBasis(nEl,nBasisMax,Len,LMsBasis)
         ELSE
            IF(TSTARSTORE) THEN
                WRITE(6,*) "Reading 2-vertex integrals of double excitations only"
            ENDIF
            LMSBASIS=LMS
            WRITE(6,*) "TBIN:",tBin
            CALL INITFROMFCID(NEL,NBASISMAX,LEN,LMSBASIS,TBIN)
!C.. say we're a UHF det so all singles are 0
            IF(LMS.EQ.0) THEN
               NBASISMAX(4,5)=1
            ELSE
               NBASISMAX(4,5)=2
            ENDIF
            IF(NBASISMAX(2,3).EQ.1) then
                WRITE(6,*) "Unrestricted calculation.  Cave Arthropodia"
            ENDIF
         ENDIF 
      ELSE   

      IF(TUEG) THEN
         WRITE(6,*) ' *** UNIFORM ELECTRON GAS CALCULATION ***' 
         IF(FUEGRS.NE.0.D0) THEN
            WRITE(6,*) 'Electron Gas Rs set to ',FUEGRS
            OMEGA=BOX*BOX*BOX*BOA*COA
!C.. required density is (3/(4 pi rs^3))
!C.. need omega to be (NEL* 4 pi rs^3 / 3)
!C.. need box to be (NEL*4 pi/(3 BOA COA))^(1/3) rs
            BOX=(NEL*4.D0*PI/(3.D0*BOA*COA))**(1.D0/3.D0)
            BOX=BOX*FUEGRS
            WRITE(6,*) "Resetting box size to ", BOX
         ENDIF
      ENDIF
      IF(THUB) WRITE(6,*) ' *** HUBBARD MODEL ***' 
!C..
      IF(.NOT.THUB.AND..NOT.TUEG) THEN
         WRITE(6,*) "Electron in cubic box."
         IF(TPARITY) THEN
            WRITE(6,*) ' ******************************* '
            WRITE(6,*) ' PARITY IS ON '
            DO I=1,3
               IF(IPARITY(I).EQ.1) THEN
                  CPAR(I)='G'
                ELSEIF(IPARITY(I).EQ.-1) THEN
                  CPAR(I)='U'
                ELSE 
                  STOP ' !!! PROBLEM WITH PARITY !!! '
                ENDIF
            ENDDO
            CPARITY=CPAR(1)//CPAR(2)//CPAR(3)
            WRITE(6,*) ' PARITY : ' , CPARITY
         ELSE
            WRITE(6,*) ' PARITY IS OFF '
         ENDIF
         WRITE(6,*) ' ******************************* '
         
!  //TBR
!         IF((.NOT.TBLOCK).AND.(.NOT.TCALCHMAT.OR.NTAY.LT.0)) STOP 'CANNOT USE PARITY WITHOUT LIST OF DETS' 
      ELSE
         IF(TPARITY) THEN
            WRITE(6,*) ' MOMENTUM : ',(IPARITY(I),I=1,3)
         ENDIF
      ENDIF
!C..
      NMAX=MAX(NMAXX,NMAXY,NMAXZ)
      NNR=NMSH*NMSH*NMSH
      WRITE(6,*) ' NMAXX : ' , NMAXX
      WRITE(6,*) ' NMAXY : ' , NMAXY
      WRITE(6,*) ' NMAXZ : ' , NMAXZ
      WRITE(6,*) ' NMSH : ' , NMSH 
!C.. 2D check
      IF(NMAXZ.EQ.0) THEN
         WRITE(6,*) 'NMAXZ=0.  2D calculation using C/A=1/A'
         COA=1/BOX
      ENDIF
      PI=ACOS(-1.D0)

!C..
      IF(THUB) THEN
         WRITE(6,'(1X,A,F19.5)') ' HUBBARD T : ' , BHUB
         WRITE(6,'(1X,A,F19.5)') ' HUBBARD U : ' , UHUB
         IF(TTILT) WRITE(6,*) 'TILTED LATTICE: ',ITILTX, ",",ITILTY
         IF(TTILT.AND.ITILTX.GT.ITILTY) STOP 'ERROR: ITILTX>ITILTY'
      ELSE
         WRITE(6,'(1X,A,F19.5)') ' BOX LENGTH : ' , BOX
         WRITE(6,'(1X,A,F19.5)') ' B/A : ' , BOA
         WRITE(6,'(1X,A,F19.5)') ' C/A : ' , COA
         TTILT=.FALSE.
      ENDIF
      ALAT(1)=BOX
      ALAT(2)=BOX*BOA
      ALAT(3)=BOX*COA
      IF(fRc.EQ.0.0.AND.iPeriodicDampingType.NE.0) THEN
         ALAT(4)=BOX*((BOA*COA)/(4*PI/3))**THIRD
      ELSE
         ALAT(4)=fRc
      ENDIF
!      ALAT(4)=2*BOX*(BOA*COA)**(1/3.D0)
      
      IF(THUB) THEN
         WRITE(6,*) ' X-LENGTH OF HUBBARD CHAIN:',(NMAXX)
         WRITE(6,*) ' Y-LENGTH OF HUBBARD CHAIN:',(NMAXY)
         WRITE(6,*) ' Z-LENGTH OF HUBBARD CHAIN:',(NMAXZ)
         WRITE(6,*) ' Periodic Boundary Conditions:',TPBC
         WRITE(6,*) ' Real space basis:',TREAL
         IF(TTILT.AND.THUB) THEN
            OMEGA=DFLOAT(NMAXX)*NMAXY*(ITILTX*ITILTX+ITILTY*ITILTY)
         ELSE
            OMEGA=DFLOAT(NMAXX)*(NMAXY)*(NMAXZ)
         ENDIF
         RS=1.D0
      ELSE
         OMEGA=ALAT(1)*ALAT(2)*ALAT(3)
         RS=(3.D0*OMEGA/(4.D0*PI*NEL))**THIRD
         ALAT(5)=RS
         IF(iPeriodicDampingType.NE.0) THEN
            IF(iPeriodicDampingType.EQ.1) THEN
               WRITE(6,*) "Using attenuated Coulomb potential for exchange interactions."
            ELSEIF(iPeriodicDampingType.EQ.2) THEN
               WRITE(6,*) "Using cut-off Coulomb potential for exchange interactions."
            ENDIF
      
            WRITE(6,*) "Rc cutoff: ",ALAT(4)
         ENDIF
         WRITE(6,*) "Wigner-Seitz radius Rs=",RS
         FKF=(9*PI/4)**THIRD/RS
         WRITE(6,*) "Fermi vector kF=",FKF
         WRITE(6,*) "Fermi Energy EF=",FKF*FKF/2
         WRITE(6,*) "Unscaled Fermi Energy nmax**2=",(FKF*FKF/2)/(0.5*(2*PI/ALAT(5))**2)
      ENDIF
      IF(OrbECutoff.gt.-1e-20) WRITE(6,*) "Orbital Energy Cutoff:",OrbECutoff
      WRITE(6,'(1X,A,F19.5)') ' VOLUME : ' , OMEGA
      WRITE(6,*) ' TALPHA : ' , TALPHA
      WRITE(6,'(1X,A,F19.5)') ' ALPHA : ' , ALPHA
      ALPHA=MIN(ALAT(1),ALAT(2),ALAT(3))*ALPHA
      WRITE(6,'(1X,A,F19.5)') ' SCALED ALPHA : ' , ALPHA

!C..
!C..Calculate number of basis functions
!C.. UEG allows from -NMAX->NMAX      
      IF(TSPINPOLAR) THEN
         NBASISMAX(4,1)=1 
!C.. spinskip
         NBASISMAX(2,3)=1
      ELSE
         NBASISMAX(4,1)=-1
!C.. spinskip
!  If spinskip is unset
         IF(NBASISMAX(2,3).EQ.0) NBASISMAX(2,3)=2
      ENDIF
      NBASISMAX(4,2)=1
      IF(THUB) THEN
         IF(TTILT) THEN
            CALL SETBASISLIM_HUBTILT(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,TPBC,ITILTX,ITILTY)
            IF(TREAL) STOP 'REAL TILTED HUBBARD NOT SUPPORTED'
          ELSE
            CALL SETBASISLIM_HUB(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,TPBC,TREAL)
         ENDIF
      ELSEIF(TUEG) THEN
         NBASISMAX(1,1)=-NMAXX
         NBASISMAX(1,2)=NMAXX
         NBASISMAX(2,1)=-NMAXY
         NBASISMAX(2,2)=NMAXY
         NBASISMAX(3,1)=-NMAXZ
         NBASISMAX(3,2)=NMAXZ
         NBASISMAX(1,3)=-1
         LEN=(2*NMAXX+1)*(2*NMAXY+1)*(2*NMAXZ+1)*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
!C.. UEG
         NBASISMAX(3,3)=-1
      ELSE
         NBASISMAX(1,1)=1
         NBASISMAX(1,2)=NMAXX
         NBASISMAX(2,1)=1
         NBASISMAX(2,2)=NMAXY
         NBASISMAX(3,1)=1
         NBASISMAX(3,2)=NMAXZ
         NBASISMAX(1,3)=0
         LEN=NMAXX*NMAXY*NMAXZ*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
         NBASISMAX(1,3)=0
!C.. particle in box
         NBASISMAX(3,3)=-2
      ENDIF
      ENDIF
!C..         (.NOT.TREADINT)


!C.. we actually store twice as much in arr as we need.
!C.. the ARR(1:LEN) are the energies of the orbitals ordered according to
!C.. BRR.  ARR(LEN+1:2*LEN) are the energies of the orbitals with default 
!C.. ordering.
      WRITE(6,*) "# basis", Len
      Allocate(Arr(LEN,2),STAT=ierr)
      LogAlloc(ierr,'Arr',2*LEN,8,tagArr)
! // TBR
!      IP_ARRSTORE=IP_ARR
      CALL AZZERO(ARR,2*LEN)
      Allocate(Brr(LEN),STAT=ierr)
      LogAlloc(ierr,'Brr',LEN,4,tagBrr)
      CALL IAZZERO(BRR,LEN)
      Allocate(G1(Len),STAT=ierr)
      LogAlloc(ierr,'G1',LEN,BasisFNSizeB,tagG1)
      CALL IAZZERO(G1,BasisFNSize*LEN)
      IF(TCPMD) THEN
         WRITE(6,*) ' *** INITIALIZING BASIS FNs FROM CPMD *** '
         CALL CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN) 
         NBASIS=LEN
         iSpinSkip=NBasisMax(2,3)
      ELSEIF(TREADINT.AND.TDFREAD) THEN
         WRITE(6,*) ' *** Creating Basis Fns from Dalton output ***'
         call InitDaltonBasis(nBasisMax,Arr,Brr,G1,Len)
         nBasis=Len
         call GenMolpSymTable(1,G1,nBasis,Arr,Brr)
      ELSEIF(TREADINT) THEN
         WRITE(6,*) ' *** CREATING BASIS FNs FROM FCIDUMP *** '
         CALL GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN,TBIN) 
         NBASIS=LEN
         CALL GENMOLPSYMTABLE(NBASISMAX(5,2)+1,G1,NBASIS,ARR,BRR)
      ELSE
!C.. Create plane wave basis functions
         IG=0
         DO I=NBASISMAX(1,1),NBASISMAX(1,2)
           DO J=NBASISMAX(2,1),NBASISMAX(2,2)
             DO K=NBASISMAX(3,1),NBASISMAX(3,2)
               DO L=NBASISMAX(4,1),NBASISMAX(4,2),2
                  G%k(1)=I
                  G%k(2)=J
                  G%k(3)=K
                  G%Ms=L
                  IF((THUB.AND.(TREAL.OR..NOT.TPBC)).OR.KALLOWED(G,NBASISMAX)) THEN
                    IF(THUB) THEN
!C..Note for the Hubbard model, the t is defined by ALAT(1)!
                       IF(TPBC) THEN
                       CALL HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ELSE
                      CALL HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ENDIF
                    ELSE
                       SUM=(BOX**2)*((I*I/ALAT(1)**2)+(J*J/ALAT(2)**2)+(K*K/ALAT(3)**2))
                    ENDIF
                    IF(SUM.GT.OrbECutoff) CYCLE
                    IG=IG+1
                    ARR(IG,1)=SUM
                    ARR(IG,2)=SUM
                    BRR(IG)=IG
!C..These are the quantum numbers: n,l,m and sigma
                    G1(IG)%K(1)=I
                    G1(IG)%K(2)=J
                    G1(IG)%K(3)=K
                    G1(IG)%MS=L
                    G1(IG)%Sym=TotSymRep()
                  ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
!C..Check to see if all's well
         WRITE(6,*) ' NUMBER OF BASIS FUNCTIONS : ' , IG 
         NBASIS=IG
         IF(LEN.NE.IG) THEN
            IF(OrbECutoff.gt.-1e20) then
               write(6,*) "Have removed ", LEN-IG, " high energy orbitals"
            else
               WRITE(6,*) "LEN=",LEN,"IG=",IG
               STOP ' LEN NE IG ' 
            endif
         ENDIF
         CALL GENMOLPSYMTABLE(1,G1,NBASIS,ARR,BRR)
      ENDIF
!C..        (.NOT.TREADINT)
!C.. Set the initial symmetry to be totally symmetric
      CALL IAZZERO(FrzSym,SymmetrySize)
      FrzSym%Sym=TotSymRep()
      CALL SetupFreezeSym(FrzSym)
!C..Now we sort them using SORT2 and then SORT

!C.. This sorts ARR and BRR into order of ARR [AJWT]
      CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
      CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)
      IF(NEL.GT.NBASIS) STOP 'MORE ELECTRONS THAN BASIS FUNCTIONS'
      CALL FLUSH(6)
      IF(TREAL.AND.THUB) THEN
!C.. we need to allow integrals between different spins
         NBASISMAX(2,3)=1
      ENDIF      
      
      !This is used in a test in UMatInd
      NOCC=NEl/2 
      IF(TREADINT) THEN
!C.. we're reading in integrals and have a molpro symmetry table
         CALL GENMOLPSYMREPS(NBASISMAX(5,2)+1,G1,NBASIS,ARR,BRR) 
      ELSEIF(TCPMD) THEN
!C.. If TCPMD, then we've generated the symmetry table earlier,
!C.. but we still need the sym reps table.
         CALL GENCPMDSYMREPS(G1,NBASIS,ARR,BRR,1.d-5)
      ELSEIF(THUB.AND..NOT.TREAL) THEN
         CALL GenHubMomIrrepsSymTable(G1,nBasis,nBasisMax)
         CALL GENHUBSYMREPS(nBasis/2,G1,NBASIS,ARR,BRR)
         CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)
      ELSE
!C.. no symmetry, so a simple sym table
         CALL GENMOLPSYMREPS(1,G1,NBASIS,ARR,BRR) 
      ENDIF

!C..
!// TBR
!      WRITE(6,*) ' TREAD : ' , TREAD


!// TBR
!      WRITE(6,*) ' ETRIAL : ',ETRIAL
      IF(FCOUL.NE.1.D0)  WRITE(6,*) "WARNING: FCOUL is not 1.D0. FCOUL=",FCOUL
      IF(FCOULDAMPBETA.GT.0) WRITE(6,*) "FCOUL Damping.  Beta ",FCOULDAMPBETA," Mu ",FCOULDAMPMU
         CALL TiHALT('SysInit   ',ISUB)
        End Subroutine SysInit

    Subroutine SysCleanup()
      CALL ENDSYM()
    End Subroutine SysCleanup
      END MODULE System
