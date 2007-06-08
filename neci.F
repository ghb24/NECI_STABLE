C ==------------------------------------------------------------------==
      SUBROUTINE NECICB(NEL,TSPN,LMS,TPARITY,IPARITY,NMAXX,NMAXY,
     &	NMAXZ,NMSH,BOX,BOA,COA,TALPHA,ALPHA,ISTATE,NEVAL,B2L,NCYCLE,
     &	NBLK,NKRY,TRHOOFR,TREAD,TCORR,TFODM,IOBS,JOBS,KOBS,TENERGY,TMC,
     &	NHISTBOXES,BETA,I_P,TRHOIJ,RHOEPSILON,TREADRHO,NPATHS,TMONTE,
     &  IMCSTEPS,IEQSTEPS,TTROT,I_HMAX,TCALCHMAT,NTAY,NWHTAY,I_VMAX,
     &  ILOGGING,TUEG,TBEGRAPH,THFBASIS,THUB,UHUB,BHUB,ETRIAL,TNPDERIV,
     & BETAP,BETAEQ,TMCDET,MDK,DETINV,THFCALC,NHFIT,TTILT,ITILTX,ITILTY,
     &   TBLOCK,HFMIX,HFEDELTA,TREADINT,TPBC,TREAL,NDETWORK,ORBORDER,
     &   NFROZEN,NTFROZEN,TRHF,HFCDELTA,IHFMETHOD,TREADHF,TSPECDET,
     &   SPECDET,DBETA,TREADTUMAT,THFORDER,TCSF,STOT,ICILEVEL,FUEGRS,
     &   TCPMD,TDFREAD,fRc)
      USE HElement
      USE UMatCache , only : GetUMatSize, DestroyUMatCache              &
     &   ,iPeriodicDampingType
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'cons.inc'
      INCLUDE 'irat.inc'
      INCLUDE 'vmc.inc'
      INCLUDE 'fcoul.inc'
C..Pointers
      POINTER (IP_NMRKS,NMRKS),(IP_ARR,ARR),(IP_BRR,BRR),
     &	(IP_HAMIL,HAMIL),(IP_CKN,CKN),(IP_G1,G1),
     &	(IP_LAB,LAB),(IP_NROW,NROW),(IP_CK,CK),(IP_A,A),
     &	(IP_W,W),(IP_V,V),(IP_AM,AM),(IP_BM,BM),(IP_T,T),
     &	(IP_WT,WT),(IP_SCR,SCR),(IP_ISCR,ISCR),(IP_INDEX,INDEX),
     &	(IP_WH,WH),(IP_WORK2,WORK2),(IP_V2,V2),(IP_IKILL,IKILL),
     &	(IP_TKE,TKE),(IP_UMAT,UMAT),(IP_WORK,WORK)
C..These are the CUBE specific pointers
      POINTER (IP_FCK,FCK),(IP_COEFF,COEFF),(IP_ZIA,ZIA),
     &	(IP_RHO,RHO),(IP_SITAB,SITAB),(IP_XCHOLE,XCHOLE),
     &	(IP_PSIR,PSIR),(IP_DLINE,DLINE),(IP_SUMA,SUMA),
     &	(IP_SCRTCH,SCRTCH),(IP_TMAT,TMAT)
C..
      POINTER (IP_RHOIJ,RHOIJ),(IP_RHOIJLAB,RHOIJLAB),
     &   (IP_RHOIJROWS,RHOIJROWS),(IP_RHOIJ1,RHOIJ1),
     &   (IP_RHOIJ2,RHOIJ2)
C..
      SAVE IP_NMRKS,IP_ARR,IP_BRR,IP_HAMIL,IP_CKN,
     &	IP_G1,IP_LAB,IP_NROW,IP_CK,IP_A,IP_W,IP_V,IP_AM,IP_BM,
     &	IP_T,IP_WT,IP_SCR,IP_ISCR,IP_INDEX,IP_WH,IP_WORK2,IP_V2,
     &	IP_IKILL,IP_TKE,IP_UMAT,IP_WORK
C..Save the Cube specific pointers
      SAVE IP_FCK,IP_COEFF,IP_ZIA,IP_RHO,IP_SITAB,IP_XCHOLE,
     &	IP_PSIR,IP_DLINE,IP_SUMA,IP_SCRTCH
C..
      SAVE IP_RHOIJ,IP_RHOIJLAB,IP_RHOIJROWS,IP_RHOIJ1,IP_RHOIJ2
#ifdef POINTER8
      INTEGER*8 GET_ADDR
#else
      INTEGER GET_ADDR
#endif
C.. NEL now includes all electrons, and we'll only want active later
      INTEGER NMRKS(NEL-NFROZEN,*),LAB(*),NROW(*),IKILL(*)
      INTEGER IPARITY(5),GC,IMCSTEPS,IEQSTEPS,I_HMAX,NTAY,NWHTAY
      REAL*8 ARR(*)
      TYPE(HElement) CKN(*),CK(*)
      INTEGER BRR(*)
      REAL*8 A(*),W(*),V(*),AM(*),BM(*),T(*),WT(*),SCR(*)
      TYPE(HElement) TMAT(*)
      REAL*8 WH(*),WORK2(*),V2(*),RHO(*),SITAB(*),TKE(*)
      TYPE(HElement) UMAT(*)
      TYPE(HElement) HAMIL(*)
      INCLUDE 'basis.inc'

      TYPE(Symmetry) FrzSym
      TYPE(BASISFN) G1(*),SymRestrict

C..Cube arrays
      REAL*8 ALAT(5),XCHOLE(*),SUMA(*)
      REAL*8 PSIR(*),DLINE(*),SCRTCH(*)
      COMPLEX*16 FCK(*)
      COMPLEX*16 COEFF(*)
      COMPLEX*16 ZIA(*)
      CHARACTER CPAR(3)*1,CPARITY*3
C..
      INTEGER INDEX(*),ISCR(*)
      LOGICAL TSPN,TALPHA,TPARITY,TCORR,TFODM,TREAD,TRHOOFR,TENERGY,TMC,
     &   TRHOIJ,TREADRHO,TMONTE,TTROT,TCALCHMAT,TRHOIJND,TUEG,TDFREAD
      REAL*8 RH,ETRIAL,ECORE
      REAL*8 RHOIJ(*),RHOIJ1(*),RHOIJ2(*)
      INTEGER RHOIJLAB(*),RHOIJROWS(*)
      INTEGER NBASISMAX(5,7),I_VMAX,NPATHS,NBASIS
      LOGICAL TSPINPOLAR,TBEGRAPH,THFBASIS,THUB,TNPDERIV,THFCALC
      POINTER (IP_HFBASIS,HFBASIS),(IP_UMAT2,UMAT2),(IP_HFE,HFE)
      POINTER (IP_TMAT2,TMAT2)
      REAL*8 HFBASIS(*),HFE(*)
      TYPE(HElement) TMAT2(*)

      TYPE(HElement) UMAT2(*)
      REAL*8 UHUB,BHUB,HFCDELTA
      LOGICAL TWARN


      POINTER (IP_NBLOCKSTARTS,NBLOCKSTARTS),(IP_BLOCKSYM,BLOCKSYM)
      INTEGER NBLOCKS,NBLOCKSTARTS(*)
      TYPE(BasisFN) BLOCKSYM(*)
      INTEGER FDET(NEL),IFDET
      LOGICAL TBLOCK

      LOGICAL TMCDET,TRHF
      INTEGER MCDET(NEL),DETINV
      TYPE(BasisFN) MDK,G,iSym
      LOGICAL TALLOW,TTILT,TREADINT,TPBC,TREAL,TREADHF
      LOGICAL KALLOWED
      INTEGER ORBORDER(8,2),ITILTX,ITILTY,NDETWORK
      INTEGER NFROZEN,IHFMETHOD
      INTEGER SPECDET(NEL-NFROZEN)
      LOGICAL TSPECDET,TREADTUMAT
      INTEGER IDEG
      LOGICAL THFORDER,TCSF
C.. The spin including the total spin
      INTEGER LMS2,STOT,LMSBASIS
      INTEGER ICILEVEL
      REAL*8 FUEGRS
      LOGICAL TCPMD
      REAL*8 fRc
      INCLUDE 'csf.inc'
      INCLUDE 'uhfdet.inc'
C ==-------------------------------------------------------------------==
      CALL TIMETAG()
      ECORE=0.D0

      CALL TISET('NECICUBE  ',ISUB)

C..      CALL DOZSUM()

C ==-------------------------------------------------------------------==
C..Input parameters
      WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
      WRITE(6,*) ' TSPN : ' , TSPN
      CALL IAZZERO(NBASISMAX,5*7)
      TSPINPOLAR=.FALSE.
      DO I=1,3
         SymRestrict%k(I)=IPARITY(I)
      ENDDO
      SymRestrict%Ms=IPARITY(4)
      SymRestrict%Sym%s=IPARITY(5)
      IF(TSPN) THEN
C.. If we're doing monte carlo without generating a list of
C.. determinants, we cannot as yet force spin or parity, except by
C.. restricting the basis set.  This will only work for Ms=NEL/2
C         IF((.NOT.TBLOCK).AND.(.NOT.TCALCHMAT.OR.NTAY.LE.0)) THEN
C            WRITE(6,*) 'TSPN set to TRUE.  Determinant list not being',
C     &         ' used for MC.  Forcing MS=Nel/2'
C            IF(MOD(NEL,2).EQ.0) THEN
C               LMS=NEL/2
C            ELSE
C               LMS=NEL
C            ENDIF
C            TSPINPOLAR=.TRUE.
C         ENDIF
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
C.. Encode the symmetry for the total spin in LMS
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
         IF(THFORDER) WRITE(6,*)
     &      "Ordering according to 1-electron energies."
      ELSEIF(TREADINT) THEN
C.. we read in the integrals from FCIDUMP and ignore most config
C..   
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
            LMSBASIS=LMS
            CALL INITFROMFCID(NEL,NBASISMAX,LEN,LMSBASIS)
C.. say we're a UHF det so all singles are 0
            IF(LMS.EQ.0) THEN
               NBASISMAX(4,5)=1
            ELSE
               NBASISMAX(4,5)=2
            ENDIF
         ENDIF 
      ELSE   

      IF(TUEG) THEN
         WRITE(6,*) ' *** UNIFORM ELECTRON GAS CALCULATION ***' 
         IF(FUEGRS.NE.0.D0) THEN
            WRITE(6,*) 'Electron Gas Rs set to ',FUEGRS
            OMEGA=BOX*BOX*BOX*BOA*COA
C.. required density is (3/(4 pi rs^3))
C.. need omega to be (NEL* 4 pi rs^3 / 3)
C.. need box to be (NEL*4 pi/(3 BOA COA))^(1/3) rs
            BOX=(NEL*4.D0*PI/(3.D0*BOA*COA))**(1.D0/3.D0)
            BOX=BOX*FUEGRS
            WRITE(6,*) "Resetting box size to ", BOX
         ENDIF
      ENDIF
      IF(THUB) WRITE(6,*) ' *** HUBBARD MODEL ***' 
C..
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
         IF((.NOT.TBLOCK).AND.(.NOT.TCALCHMAT.OR.NTAY.LT.0))
     &      STOP 'CANNOT USE PARITY WITHOUT LIST OF DETS' 
      ELSE
         IF(TPARITY) THEN
            WRITE(6,*) ' MOMENTUM : ',(IPARITY(I),I=1,3)
         ENDIF
      ENDIF
C..
      NMAX=MAX(NMAXX,NMAXY,NMAXZ)
      NNR=NMSH*NMSH*NMSH
      WRITE(6,*) ' NMAXX : ' , NMAXX
      WRITE(6,*) ' NMAXY : ' , NMAXY
      WRITE(6,*) ' NMAXZ : ' , NMAXZ
      WRITE(6,*) ' NMSH : ' , NMSH 
C.. 2D check
      IF(NMAXZ.EQ.0) THEN
         WRITE(6,*) 'NMAXZ=0.  2D calculation using C/A=1/A'
         COA=1/BOX
      ENDIF
      PI=ACOS(-1.D0)

C..
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
               WRITE(6,*) "Using attenuated Coulomb potential for " //
     &            "exchange interactions."
            ELSEIF(iPeriodicDampingType.EQ.2) THEN
               WRITE(6,*) "Using cut-off Coulomb potential for " //
     &            "exchange interactions."
            ENDIF
      
            WRITE(6,*) "Rc cutoff: ",ALAT(4)
         ENDIF
         WRITE(6,*) "Wigner-Seitz radius Rs=",RS
         FKF=(9*PI/4)**THIRD/RS
         WRITE(6,*) "Fermi vector kF=",FKF
         WRITE(6,*) "Fermi Energy EF=",FKF*FKF/2
         WRITE(6,*) "Unscaled Fermi Energy nmax**2=",
     &      (FKF*FKF/2)/(0.5*(2*PI/ALAT(5))**2)
      ENDIF
      WRITE(6,'(1X,A,F19.5)') ' VOLUME : ' , OMEGA
      WRITE(6,*) ' TALPHA : ' , TALPHA
      WRITE(6,'(1X,A,F19.5)') ' ALPHA : ' , ALPHA
      ALPHA=MIN(ALAT(1),ALAT(2),ALAT(3))*ALPHA
      WRITE(6,'(1X,A,F19.5)') ' SCALED ALPHA : ' , ALPHA
      IF(THFBASIS) THEN
         WRITE(6,*) "Using Hartree-Fock Basis"
         IF(.NOT.THFCALC) WRITE(6,*) "Reading Hartree-Fock Basis"
      ENDIF
      IF(THFCALC) THEN
         WRITE(6,*) "Calculating Hartree-Fock Basis"
         WRITE(6,*) "Max Iterations:",NHFIT
         WRITE(6,*) "FMIX,EDELTA",HFMIX,HFEDELTA
      ENDIF
      IF(G_VMC_FAC.LE.0) THEN
         WRITE(6,*) "G_VMC_FAC=",G_VMC_FAC
         STOP "G_VNC_FAC LE 0"
      ENDIF
C..
C..Calculate number of basis functions
C.. UEG allows from -NMAX->NMAX      
      IF(TSPINPOLAR) THEN
         NBASISMAX(4,1)=1 
C.. spinskip
         NBASISMAX(2,3)=1
      ELSE
         NBASISMAX(4,1)=-1
C.. spinskip
         NBASISMAX(2,3)=2
      ENDIF
      NBASISMAX(4,2)=1
      IF(THUB) THEN
         IF(TTILT) THEN
            CALL SETBASISLIM_HUBTILT(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,
     &         TPBC,ITILTX,ITILTY)
            IF(TREAL) STOP 'REAL TILTED HUBBARD NOT SUPPORTED'
          ELSE
            CALL SETBASISLIM_HUB(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,TPBC,
     &            TREAL)
         ENDIF
      ELSEIF(TUEG) THEN
         NBASISMAX(1,1)=-NMAXX
         NBASISMAX(1,2)=NMAXX
         NBASISMAX(2,1)=-NMAXY
         NBASISMAX(2,2)=NMAXY
         NBASISMAX(3,1)=-NMAXZ
         NBASISMAX(3,2)=NMAXZ
         NBASISMAX(1,3)=-1
         LEN=(2*NMAXX+1)*(2*NMAXY+1)*(2*NMAXZ+1)
     &      *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
C.. UEG
         NBASISMAX(3,3)=-1
      ELSE
         NBASISMAX(1,1)=1
         NBASISMAX(1,2)=NMAXX
         NBASISMAX(2,1)=1
         NBASISMAX(2,2)=NMAXY
         NBASISMAX(3,1)=1
         NBASISMAX(3,2)=NMAXZ
         NBASISMAX(1,3)=0
         LEN=NMAXX*NMAXY*NMAXZ
     &      *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
         NBASISMAX(1,3)=0
C.. particle in box
         NBASISMAX(3,3)=-2
      ENDIF
      ENDIF
C..         (.NOT.TREADINT)
C.. Set the initial symmetry to be totally symmetric
      CALL IAZZERO(FrzSym,SymmetrySize)
      FrzSym%s=1
      CALL SetupFreezeSym(FrzSym)


C.. we actually store twice as much in arr as we need.
C.. the ARR(1:LEN) are the energies of the orbitals ordered according to
C.. BRR.  ARR(LEN+1:2*LEN) are the energies of the orbitals with default 
C.. ordering.
      CALL MEMORY(IP_ARR,2*LEN,'ARR')
      IP_ARRSTORE=IP_ARR
      CALL AZZERO(ARR,2*LEN)
      CALL MEMORY(IP_BRR,LEN,'BRR')
      CALL AZZERO(BRR,LEN)
      CALL MEMORY(IP_G1,BASISFNSIZE*LEN,'G1')
      CALL IAZZERO(G1,5*LEN)
      IF(TCPMD) THEN
         WRITE(6,*) ' *** INITIALIZING BASIS FNs FROM CPMD *** '
         CALL CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN) 
         NBASIS=LEN
         NHG=LEN
      ELSEIF(TREADINT.AND.TDFREAD) THEN
         WRITE(6,*) ' *** Creating Basis Fns from Dalton output ***'
         call InitDaltonBasis(nBasisMax,Arr,Brr,G1,Len)
         nBasis=Len
         nHG=Len
         call GenMolpSymTable(1,G1,nBasis,Arr,Brr)
      ELSEIF(TREADINT) THEN
         WRITE(6,*) ' *** CREATING BASIS FNs FROM FCIDUMP *** '
         CALL GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN) 
         NBASIS=LEN
         NHG=LEN
         CALL GENMOLPSYMTABLE(NBASISMAX(5,2)+1,G1,NBASIS,ARR,BRR)
      ELSE
C.. Create plane wave basis functions
         IG=0
         DO I=NBASISMAX(1,1),NBASISMAX(1,2)
           DO J=NBASISMAX(2,1),NBASISMAX(2,2)
             DO K=NBASISMAX(3,1),NBASISMAX(3,2)
               DO L=NBASISMAX(4,1),NBASISMAX(4,2),2
                  G%k(1)=I
                  G%k(2)=J
                  G%k(3)=K
                  G$Ms=L
                  IF((THUB.AND.(TREAL.OR..NOT.TPBC))
     &               .OR.KALLOWED(G,NBASISMAX)) THEN
                    IG=IG+1
                    IF(THUB) THEN
C..Note for the Hubbard model, the t is defined by ALAT(1)!
                       IF(TPBC) THEN
                       CALL HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ELSE
                      CALL HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ENDIF
                    ELSE
                       SUM=(BOX**2)*((I*I/ALAT(1)**2)+(J*J/ALAT(2)**2)+
     &                   (K*K/ALAT(3)**2))
                    ENDIF
                    ARR(IG)=SUM
                    BRR(IG)=IG
C..These are the quantum numbers: n,l,m and sigma
                    G1(IG)%K(1)=I
                    G1(IG)%K(2)=J
                    G1(IG)%K(3)=K
                    G1(IG)%MS=L
                    G1(IG)%Sym%s=1
                  ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
C..Check to see if all's well
         WRITE(6,*) ' NUMBER OF BASIS FUNCTIONS : ' , IG 
         NBASIS=IG
         IF(LEN.NE.IG) THEN
            WRITE(6,*) "LEN=",LEN,"IG=",IG
            STOP ' LEN NE IG ' 
         ENDIF
         CALL GENMOLPSYMTABLE(1,G1,NBASIS,ARR,BRR)
      ENDIF
C..        (.NOT.TREADINT)
C..Now we sort them using SORT2 and then SORT

C.. This sorts ARR and BRR into order of ARR [AJWT]
      CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
      NHG=NBASIS
      CALL WRITEBASIS(6,G1,NHG,ARR,BRR)
      IF(NEL.GT.NBASIS) STOP 'MORE ELECTRONS THAN BASIS FUNCTIONS'
      CALL FLUSH(6)
      IF(TREAL.AND.THUB) THEN
C.. we need to allow integrals between different spins
         NBASISMAX(2,3)=1
      ENDIF      

C.. Calculate relevant symmetry tables
      IF(TREADINT) THEN
C.. we're reading in integrals and have a molpro symmetry table
         CALL GENMOLPSYMREPS(NBASISMAX(5,2)+1,G1,NBASIS,ARR,BRR) 
      ELSEIF(TCPMD) THEN
C.. If TCPMD, then we've generated the symmetry table earlier,
C.. but we still need the sym reps table
         CALL GENCPMDSYMREPS(G1,NBASIS,ARR,BRR)
      ELSEIF(THUB.AND..NOT.TREAL) THEN
         CALL GenHubMomIrrepsSymTable(G1,nBasis,nBasisMax)
         CALL GENHUBSYMREPS(nBasis/2,G1,NBASIS,ARR,BRR)
         CALL WRITEBASIS(6,G1,NHG,ARR,BRR)
      ELSE
C.. no symmetry, so a simple sym table
         CALL GENMOLPSYMREPS(1,G1,NBASIS,ARR,BRR) 
      ENDIF

C..
      WRITE(6,*) ' ISTATE : ' , ISTATE 
      WRITE(6,*) ' NEVAL : ' , NEVAL
      WRITE(6,'(1X,A,E19.3)') ' B2LIMIT : ' , B2L
      WRITE(6,*) ' NCYCLE : ' , NCYCLE
      WRITE(6,*) ' NBLK : ' , NBLK 
      WRITE(6,*) ' NKRY : ' , NKRY
      WRITE(6,*) ' TREAD : ' , TREAD
      WRITE(6,*) ' TCORR : ' , TCORR
      WRITE(6,*) ' TENERGY : ' , TENERGY
      IF(TCORR) THEN
        WRITE(6,*) 
     &	' *** EXCHANGE-CORRELATION HOLE WILL BE CALCULATED *** ' 
      ENDIF
      WRITE(6,*) ' IOBS : ' , IOBS 
      WRITE(6,*) ' JOBS : ' , JOBS 
      WRITE(6,*) ' KOBS : ' , KOBS 

      IF(BETAP.NE.0) THEN 
         I_P=NINT(BETA/BETAP)
         WRITE(6,*) ' BETAP=',BETAP
         WRITE(6,*) ' RESETTING P '
         IF(I_P.GT.100000) WRITE(6,*) ' *** WARNING I_P=',I_P
      ENDIF

      WRITE(6,*) ' BETA, P :',BETA,I_P

      IF(IOBS.GT.NMSH.OR.IOBS.LE.0.OR.JOBS.GT.NMSH.OR.JOBS.LE.0.
     &	OR.KOBS.GT.NMSH.OR.KOBS.LE.0) THEN
        STOP ' !!! REFERENCE PARTICLE NOT IN BOX !!! '
      ENDIF
      WRITE(6,*) ' ETRIAL : ',ETRIAL
      IF(TMONTE) THEN 
         WRITE(6,*) 'MC Determinant Symmetry:'
         WRITE(6,*) (MDK(I),I=1,4)
      ENDIF
      IF(FCOUL.NE.1.D0)
     &   WRITE(6,*) "WARNING: FCOUL is not 1.D0. FCOUL=",FCOUL
      IF(FCOULDAMPBETA.GT.0)
     &   WRITE(6,*) "FCOUL Damping.  Beta ",FCOULDAMPBETA,
     &      " Mu ",FCOULDAMPMU
      IF(.NOT.(TREAD.OR.TREADRHO)) THEN
         IF(NTFROZEN.LT.0) THEN
            I=NEL+(-NTFROZEN)
         ELSE
            I=NHG-NTFROZEN
         ENDIF
         IF(TCPMD) THEN
C.. We don't need to do init any 4-index integrals, but we do need to init the 2-index
            WRITE(6,*) " *** INITIALIZING CPMD 2-index integrals ***"
            CALL MEMORY(IP_UMAT,1,'UMAT')
            CALL MEMORY(IP_TMAT,HElementSize*NHG*NHG,'TMAT')
            CALL AZZERO(TMAT,HElementSize*NHG*NHG)
            CALL CPMDINIT2INDINT(TMAT,NHG,I,NBASISMAX,
     &         ISPINSKIP,G1,NEL,ECORE,THFORDER,ARR,BRR)
         ELSEIF(TREADINT.AND.TDFREAD) THEN
            CALL MEMORY(IP_UMAT,HElementSize,'UMAT')
            CALL MEMORY(IP_TMAT,HElementSize*NHG*NHG,'TMAT')
            CALL AZZERO(TMAT,HElementSize*NHG*NHG)
            Call ReadDalton1EIntegrals(G1,nBasis,TMat,Arr,Brr,ECore)
            Call ReadDF2EIntegrals(nBasis,I)
         ELSEIF(TREADINT) THEN
            WRITE(6,*)
     &         ' *** READING PRIMITIVE INTEGRALS FROM FCIDUMP *** '
C.. Generate the 2e integrals (UMAT)
            NBASISMAX(2,3)=2
C.. We can only do restricted systems, so alpha and beta spinorbital
C... integrals will be the same 
            CALL GetUMatSize(nHG,2,I)
            CALL MEMORY(IP_UMAT,HElementSize*I,'UMAT')
            CALL AZZERO(UMAT,HElementSize*I)
            CALL MEMORY(IP_TMAT,HElementSize*NHG*NHG,'TMAT')
            CALL AZZERO(TMAT,HElementSize*NHG*NHG)
            CALL READFCIINT(UMAT,TMAT,NBASIS,ECORE,ARR,BRR,G1)
            WRITE(6,*) ' ECORE=',ECORE
            ISPINSKIP=2
         ELSE
C.. We need to init the arrays regardless of whether we're storing H
C..Need to initialise the Fourier arrays
            CALL MEMORY(IP_FCK,2*NMSH**3,'FCK')
            CALL MEMORY(IP_COEFF,2*(3*NMSH+48),'COEFF')
            CALL MEMORY(IP_ZIA,2*(NMSH+1)*NMAX*NMAX,'ZIA')
C..
            CALL MEMORY_CHECK()
            IF(NMAXZ.EQ.0) THEN
C..   We're doing a 2D simulation
               CALL INITFOU2D(NMSH,FCK,COEFF,NMAX,ZIA,ALAT,TALPHA,
     &            ALPHA,OMEGA)
            ELSE
               CALL INITFOU(NMSH,FCK,COEFF,NMAX,ZIA,ALAT,TALPHA,
     &            ALPHA,OMEGA)
            ENDIF
            CALL MEMORY_CHECK()
            ISPINSKIP=NBASISMAX(2,3)
            IF(NBASISMAX(1,3).GE.0) THEN
C.. we pre-compute the 2-e integrals
               WRITE(6,*) "Generating 2e integrals"
C.. Generate the 2e integrals (UMAT)
               CALL GetUMatSize(nHG,iSpinSkip,I)
               CALL MEMORY(IP_UMAT,HElementSize*I,
     &            'UMAT')
               CALL AZZERO(UMAT,HElementSize*I)
               IF(TUEG.OR.THUB) THEN
                  IF(THUB.AND.TREAL) THEN
C.. Real space hubbard
                     CALL CALCUMATHUBREAL(NEL,NBASIS,NBASISMAX,G1,UHUB,
     &                  UMAT)
                  ELSEIF(THUB.AND..NOT.TPBC) THEN
C.. Non-periodic hubbard (mom space)
                     CALL GEN_COUL_HUBNPBC(NEL,NBASISMAX,NHG,G1,
     &                 NMSH,NMAX,FCK,ZIA,UMAT,ISPINSKIP,THUB,UHUB,OMEGA)
                  ELSE
C.. Most normal Hubbards
                     IF(.NOT.TUEG) THEN
!                        CALL GEN_COUL_UEG(NEL,NBASISMAX,NHG,G1,
!     &          NMSH,NMAX,FCK,ZIA,UMAT,ISPINSKIP,THUB,UHUB,OMEGA,ALAT)
                        ISPINSKIP=-1
                        NBASISMAX(2,3)=-1
                     WRITE(6,*) "Not precomputing HUBBARD 2-e integrals"
                        CALL MEMORY(IP_UMAT,1,'UMAT')
                        UMAT(1)=UHUB/OMEGA
                     ENDIF
C.. The UEG doesn't store coul integrals
                  ENDIF
               ELSE
                  CALL GEN_COUL(NEL,NBASISMAX,NHG,G1,
     &               NMSH,NMAX,FCK,ZIA,UMAT,ISPINSKIP)
               ENDIF
            ELSE
               WRITE(6,*) "Not precomputing 2-e integrals"
               CALL MEMORY(IP_UMAT,1,'UMAT')
            ENDIF
            CALL MEMORY_CHECK()
C.. we need to generate TMAT
            CALL MEMORY(IP_TMAT,HElementSize*NHG*NHG,'TMAT')
            CALL AZZERO(TMAT,HElementSize*NHG*NHG)
            IF(THUB) THEN
             CALL CALCTMATHUB(TMAT,NBASIS,NBASISMAX,BHUB,TTILT,G1,TREAL,
     &         TPBC)
            ELSE
C..Cube multiplier
               CST=PI*PI/(2.D0*ALAT(1)*ALAT(1))
C.. the UEG has k=2pi n/L rather than pi n/L, so we need to multiply the
C.. KE up by 4
               IF(NBASISMAX(1,1).LE.0) CST=CST*4
               CALL CALCTMATUEG(TMAT,NBASIS,ALAT,G1,CST,
     &            NBASISMAX(1,1).LE.0,OMEGA)
            ENDIF
         ENDIF
C..           (.NOT.TREADINT)
         FDET(1)=0
C.. If we are using an HF basis instead of our primitive basis, we need
C.. to load in the coeffs of the HF eigenfunctions in terms of the
C.. primitive basis.
C.. We load the coeffs from a file HFBASIS
         IF(THFBASIS.OR.THFCALC.OR.(THFORDER.AND..NOT.TCPMD)) THEN
            CALL MEMORY(IP_HFBASIS,NHG*NHG,'HFBASIS')
            CALL AZZERO(HFBASIS,NHG*NHG)
C.. Allocate an array to store the HF Energies
            CALL MEMORY(IP_HFE,NHG,'HFE')
            CALL AZZERO(HFE,NHG)
            IF(THFORDER.AND..NOT.THFBASIS) THEN
C.. If we're not using HF, but just calculating the HF order
C.. We generate the HF energies (this has no mixing or randomisation, so should jsut
C.. re-order the orbitals and give us some energy)
C.. HF basis is NOT using the LMS value set in the input
              CALL CALCHFBASIS(NBASIS,ISPINSKIP,NBASISMAX,G1,ARR,BRR,
     &         ECORE,TMAT,UMAT,HFE,HFBASIS,1,NEL,LMSBASIS,1.D0,HFEDELTA,
     &            HFCDELTA,.TRUE.,0,TREADHF,0.D0,FDET,ILOGGING)
               CALL ORDERBASISHF(ARR,BRR,HFE,HFBASIS,G1,NBASIS,FDET,NEL)
            ELSEIF(THFCALC) THEN
              CALL CALCHFBASIS(NBASIS,ISPINSKIP,NBASISMAX,G1,ARR,BRR,
     &         ECORE,TMAT,UMAT,HFE,HFBASIS,NHFIT,NEL,LMS,HFMIX,HFEDELTA,
     &            HFCDELTA,TRHF,IHFMETHOD,TREADHF,0.01D0,FDET,ILOGGING)
               CALL SETUPHFBASIS(NBASISMAX,G1,NBASIS,HFE,ARR,BRR)
            ELSEIF(THFBASIS) THEN
               CALL READHFBASIS(HFBASIS,HFE,G1,NBASIS)
               CALL SETUPHFBASIS(NBASISMAX,G1,NBASIS,HFE,ARR,BRR)
            ENDIF
            WRITE(6,*) "FINAL HF BASIS"
            CALL WRITEBASIS(6,G1,NHG,ARR,BRR)

            WRITE(6,"(A,$)") "Fermi det (D0):"
            CALL WRITEDET(6,FDET,NEL,.TRUE.)
            CALL FLUSH(6)
C.. If in Hubbard, we generate site-spin occupations
            IF(THUB) THEN
               CALL GENSITESPINOCC(NBASIS,NBASIS/ISPINSKIP,ISPINSKIP,
     &            NBASISMAX,G1,NEL,LMS,BRR,HFBASIS)
            ENDIF
C.. We now generate a new U matrix corresponding to the HF basis fns
C.. This requires a new matrix UMAT2 to be set up, as our HF basis is
C.. unrestricted
            IF(THFBASIS) THEN
               WRITE(6,*) "Allocating TMAT2"
               CALL MEMORY(IP_TMAT2,NHG**2,'TMAT2')
               CALL AZZERO(TMAT2,NHG**2)
               NORBUSED=NBASIS-NTFROZEN
               IF(TREADTUMAT) THEN
                  CALL READHFTMAT(TMAT,TMAT2,NBASIS,HFBASIS)
               ELSE
                  CALL CALCHFTMAT(TMAT,TMAT2,NBASIS,HFBASIS,NORBUSED)
               ENDIF
               CALL FREEM(IP_TMAT)
               IP_TMAT=IP_TMAT2
               IP_TMAT2=NULL
C.. Allocate the new matrix
               CALL GetUMatSize(nHG,1,I)
               CALL MEMORY(IP_UMAT2,HElementSize*I,'UMAT2')
               CALL AZZERO(UMAT2,HElementSize*I)
C.. We need to pass the TMAT to CALCHFUMAT as TMAT is no longer diagona
C.. This also modified G1, ARR, BRR
               IF(TREADTUMAT) THEN
                 CALL READHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,
     &            ISPINSKIP,HFE,ARR,BRR)
               ELSE
                 CALL CALCHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,
     &            ISPINSKIP,HFE,ARR,BRR,NORBUSED)
               ENDIF
C.. Now we can remove the old UMATRIX, and set the pointer UMAT to point
C.. to UMAT2
               CALL FREEM(IP_UMAT)
               IP_UMAT=IP_UMAT2
               IP_UMAT2=NULL
C.. spinskip values
               ISPINSKIP=1
               NBASISMAX(2,3)=1
C.. Indicate that we're in UHF and that <D0|H|D1>=0 for D1 being a
C.. single excitation
               IF(LMS.EQ.0) THEN
                  NBASISMAX(4,5)=1
                  DO I=1,NEL
                     NUHFDET(I)=BRR(I)
                  ENDDO
                  CALL SORTI(NEL,UHFDET)
               ELSE
                  NBASISMAX(4,5)=2
               ENDIF
            ENDIF
C.. Now deallocate the HF arrays
            CALL FREEM(IP_HFBASIS)
            CALL FREEM(IP_HFE)
            CALL FLUSH(6)
         ENDIF

         CHEMPOT=(ARR(NEL)+ARR(NEL+1))/2.D0
         WRITE(6,*) "Chemical Potential: ",CHEMPOT
         IF(NTFROZEN.LT.0) THEN
            WRITE(6,*) "NTFROZEN<0.  Leaving ", -NTFROZEN,
     &         " unfrozen virtuals."
            NTFROZEN=NTFROZEN+NHG-NEL
         ENDIF
         IF(NFROZEN.GT.0.OR.NTFROZEN.GT.0) THEN
C.. At this point, we transform the UMAT and TMAT into a new UMAT and
C.. TMAT and Ecore with the frozen orbitals factored in
C..
C.. a,b are frozen spinorbitals
C.. E'core = Ecore+sum_a t_aa + sum_(a<b) (<ab|ab>-<ab|ba>)
C.. t'_ii = t_ii+ sum_a ( <ai|ai> - <ai|ia> ) 
C.. NHG contains the old number of orbitals
C.. NBASIS contains the new
            NBASIS=NBASIS-NFROZEN-NTFROZEN
C.. We need to transform some integrals
            CALL MEMORY(IP_TMAT2,HElementSize*(NBASIS)**2,'TMAT2')
            CALL AZZERO(TMAT2,HElementSize*(NBASIS)**2)
            IF(NBASISMAX(1,3).GE.0.AND.ISPINSKIP.NE.0) THEN
C.. allocate a UMAT
               CALL GetUMatSize(nBasis,iSpinSkip,I) 
               CALL MEMORY(IP_UMAT2,HElementSize*I,
     &            'UMAT2')
               CALL AZZERO(UMAT2,HElementSize*I)
            ELSE
C.. we don't precompute 4-e integrals, so don't allocate a large UMAT
               CALL MEMORY(IP_UMAT2,1,'UMAT2')
            ENDIF 
            CALL FREEZEBASIS(NHG,NBASIS,UMAT,TMAT,UMAT2,TMAT2,ECORE,
     &         G1,NBASISMAX,ISPINSKIP,ARR,BRR,NFROZEN,NTFROZEN,FDET,NEL,
     &         ALAT)
            WRITE(6,*) "Freezing ",NFROZEN," core orbitals."
            WRITE(6,*) "Freezing ",NTFROZEN," virtual orbitals."
            WRITE(6,*) "ECORE now",ECORE
            WRITE(6,*) NBASIS," orbitals remain."
            NEL=NEL-NFROZEN
C.. NEL now only includes active electrons
            WRITE(6,*) "Number of active electrons:",NEL
            CALL FREEM(IP_TMAT)
            IP_TMAT=IP_TMAT2
            IP_TMAT2=NULL
C.. Now we can remove the old UMATRIX, and set the pointer UMAT to point
C.. to UMAT2
            CALL FREEM(IP_UMAT)
            IP_UMAT=IP_UMAT2
            IP_UMAT2=NULL
            CALL MEMORY_CHECK
            WRITE(6,*) "Active basis functions:",NHG
            CALL WRITEBASIS(6,G1,NHG,ARR,BRR)
         ENDIF
         IF(COULDAMPORB.GT.0) THEN
            FCOULDAMPMU=(ARR(COULDAMPORB)+ARR(COULDAMPORB+1))/2
            WRITE(6,*) "Setting Coulomb damping mu between orbitals ",
     &         ARR(COULDAMPORB)," and ",ARR(COULDAMPORB+1)
            WRITE(6,*) "MU=",FCOULDAMPMU
         ENDIF
         IF(THFBASIS.OR.THFCALC.OR.THFORDER) THEN
C.. We now generate a new U matrix corresponding to the HF basis fns
C.. This requires a new matrix UMAT2 to be set up, as our HF basis is
C.. unrestricted
            IF(THFBASIS) THEN
               CALL MEMORY(IP_TMAT2,HElementSize*NHG**2,'TMAT2')
               CALL AZZERO(TMAT2,HElementSize*NHG**2)
               IF(TREADTUMAT) THEN
                  CALL READHFTMAT(TMAT,TMAT2,NBASIS,HFBASIS)
               ELSE
                  CALL CALCHFTMAT(TMAT,TMAT2,NBASIS,HFBASIS)
               ENDIF
               CALL FREEM(IP_TMAT)
               IP_TMAT=IP_TMAT2
               IP_TMAT2=NULL
C.. Allocate the new matrix
               CALL GetUMatSize(nHG,1,I)
               CALL MEMORY(IP_UMAT2,HElementSize*I,'UMAT2')
               CALL AZZERO(UMAT2,HElementSize*I)
C.. We need to pass the TMAT to CALCHFUMAT as TMAT is no longer diagona
C.. This also modified G1, ARR, BRR
               IF(TREADTUMAT) THEN
                 CALL READHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,
     &            ISPINSKIP,HFE,ARR,BRR)
               ELSE
                 CALL CALCHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,
     &            ISPINSKIP,HFE,ARR,BRR)
               ENDIF
C.. Now we can remove the old UMATRIX, and set the pointer UMAT to point
C.. to UMAT2
               CALL FREEM(IP_UMAT)
               IP_UMAT=IP_UMAT2
               IP_UMAT2=NULL
C.. spinskip values
               ISPINSKIP=1
               NBASISMAX(2,3)=1
C.. Indicate that we're in UHF and that <D0|H|D1>=0 for D1 being a
C.. single excitation
               IF(LMS.EQ.0) THEN
                  NBASISMAX(4,5)=1
                  DO I=1,NEL
                     NUHFDET(I)=BRR(I)
                  ENDDO
                  CALL SORTI(NEL,UHFDET)
               ELSE
                  NBASISMAX(4,5)=2
               ENDIF
            ENDIF
         ENDIF
      ELSE
C.. if reading in FCI solutions, then 
C.. no symmetry, so a simple sym table
            CALL GENMOLPSYMREPS(1,G1,NBASIS,ARR,BRR) 
      ENDIF
C           (.NOT.(TREAD.OR.TREADRHO))

      WRITE(6,*) "SYMMETRY MULTIPLICATION TABLE"
      CALL WRITESYMTABLE(6)
   
      CALL GENSymStatePairs(NBASIS/2)


      IF(TCSF.AND.TSPECDET) THEN
         WRITE(6,*) "TSPECDET set.  SPECDET is"
         CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
         CALL ICOPY(NEL,SPECDET,1,FDET,1)
         CALL GETCSFFROMDET(FDET,SPECDET,NEL,STOT,LMS)
         WRITE(6,*) "CSF with 2S=",STOT," and 2Sz=",LMS,
     &         " now in SPECDET is"
         CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
      ENDIF

C.. Work out a preliminary Fermi det
      IF(IFDET.EQ.0.AND.FDET(1).EQ.0) THEN
         CALL GENFDET(BRR,G1,NBASIS,LMS,NEL,FDET)
      ENDIF
      WRITE(6,"(A,$)") "Fermi det (D0):"
      CALL WRITEDET(6,FDET,NEL,.TRUE.)
      CALL ICOPY(NEL,FDET,1,NUHFDET,1)
      IF(TSPECDET.AND.SPECDET(1).EQ.0) THEN
         WRITE(6,*) "TSPECDET set, but invalid.  using FDET"
         CALL ICOPY(NEL,FDET,1,SPECDET,1)
      ENDIF


C.. Check if we're blocking the hamiltonian
C      IF(THFBASIS.AND.TBLOCK) THEN
C         WRITE(6,*) "THFBASIS set and NBLK=0.  ",
C     &         "Cannot block diagonalize in HF Basis."
C         STOP
C      ENDIF
C      CALL SYMGENEXCITS(FDET,NEL,G1,NBASIS,NBASISMAX)
C      CALL PMAXMEM
C      CALL TIPRI
C      STOP

      IF(.NOT.TCALCHMAT) THEN
         WRITE(6,*) "Not storing the H matrix."
         IF(TENERGY.AND..NOT.TBLOCK) THEN
            WRITE(6,*) "Cannot calculate energies without blocking the"
     &            //" Hamiltonian."
            TENERGY=.FALSE.
         ENDIF
         IF(TENERGY.AND.NBLK.NE.0) THEN
C.. We're doing a Lanczos without calculating the H mat
            WRITE(6,*) "Cannot perform Lanczos without Hamiltonian"
            TENERGY=.FALSE.
         ENDIF
         IF(TRHOIJ) THEN
            IF(NPATHS.NE.0.OR.DETINV.GT.0) THEN
               TRHOIJND=.TRUE.
            ELSE
               WRITE(6,*) "Cannot calculate RHOs without NPATHS set."
            ENDIF
         ENDIF 
         TRHOIJ=.FALSE.
      ENDIF

C.. in order to calculate the H matrix, we need to work out all the determinants
C.. beware with NPATHS - it calcs the list of dets even if we don't calc H
C.. Could be big.
C..Now we see how many determinants we need
C      IF(NHG.GT.170) THEN
C..This fix is to stop floating overflow as FACTRL(NHG.GT.170) crashes
         NDET=1
         DNDET=1.D0
         DO I=0,NEL-1
            NDET=(NDET*(NHG-I))/(I+1)
            DNDET=(DNDET*DFLOAT(NHG-I))/DFLOAT(I+1)
         ENDDO
C      ELSE
C         NDET=FACTRL(NHG)/(FACTRL(NEL)*FACTRL(NHG-NEL)) 
C      ENDIF
      WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , NDET,DNDET
      
C      CALL TC(I_HMAX,I_P,NWHTAY)


C      IF(TCALCHMAT.OR.NPATHS.NE.0.OR.DETINV.GT.0.OR.TBLOCK) THEN
      IF(TENERGY.OR.TCALCHMAT) THEN
C..Need to determine the determinants
         IF(ICILEVEL.NE.0) THEN
            WRITE(6,*) "Performing truncated CI at level ",ICILEVEL
            IF(TSPECDET) THEN
               WRITE(6,*) "Using SPECDET:"
               CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
               CALL ICOPY(NEL,SPECDET,1,FDET,1)
            ELSE
               WRITE(6,*) "Using Fermi DET:"
               CALL WRITEDET(6,FDET,NEL,.TRUE.)
            ENDIF 
C.. if we're doing a truncated CI expansion
            CALL GENEXCIT(FDET,ICILEVEL,NBASIS,NEL,0,0,NDET,1,
     &            G1,.TRUE.,NBASISMAX,.TRUE.)
C.. We need to add in the FDET
            NDET=NDET+1
            II=NDET
            NBLOCKS=1
         ELSEIF(TBLOCK) THEN
            WRITE(6,*) "Determining determinants and blocks."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            CALL GNDTS_BLK(NEL,NHG,BRR,NBASISMAX,NMRKS, .TRUE.,
     &            NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS2,TPARITY,
     &           SymRestrict,IFDET,TMAT,.NOT.TREAD,NDETTOT,BLOCKSYM)
            WRITE(6,*) "NBLOCKS:",NBLOCKS
         ELSEIF(TCSF) THEN
            WRITE(6,*) "Determining CSFs."
C determinants.
            WRITE(6,*) "Determining determinants and blocks."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            CALL GNCSFS(NEL,NHG,BRR,NBASISMAX,
     &         NMRKS,.TRUE.,G1,TSPN,LMS2,TPARITY,
     &         SymRestrict,NDET,IFDET,.FALSE.,0,0,.FALSE.,0)
            NBLOCKS=1
            II=NDET
         ELSE
            WRITE(6,*) "Determining determinants."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            CALL GNDTS(NEL,NHG,BRR,NBASISMAX,NMRKS,.TRUE.,
     &            G1,TSPN,LMS,TPARITY,SymRestrict,II,IFDET)
            NBLOCKS=1
            NDET=II
         ENDIF
C..
         IF(II.EQ.0) THEN
            WRITE(6,*) "No determinants found.  Cannot continue"
            STOP "No determinants found.  Cannot continue"
         ENDIF
C.. NEL now only includes active electrons
         CALL MEMORY(IP_NMRKS,(NEL)*II,'NMRKS')
         CALL IAZZERO(NMRKS,(NEL)*II)
         CALL MEMORY(IP_NBLOCKSTARTS,(NBLOCKS+1)/IRAT+1,"NBLOCKSTARTS")
         CALL IAZZERO(NBLOCKSTARTS,NBLOCKS)
         CALL MEMORY(IP_BLOCKSYM,((NBLOCKS+1)*BasisFNSize)/IRAT+1,
     &         "BLOCKSYM")

         CALL IAZZERO(BLOCKSYM,NBLOCKS*BasisFNSize)
C..



         NDET=II   
         IF(ICILEVEL.NE.0) THEN
C.. Use HAMIL to temporarily hold a list of excitation levels
            CALL ICOPY(NEL,FDET,1,NMRKS,1)
            CALL MEMORY(IP_HAMIL,II*HElementSize,'HAMIL')
            CALL GENEXCIT(FDET,ICILEVEL,NBASIS,NEL,NMRKS(1,2),HAMIL,
     &            NDET,1,
     &            G1,.TRUE.,NBASISMAX,.FALSE.)
            CALL FREEM(IP_HAMIL)
            NDET=NDET+1
            NBLOCKSTARTS(1)=1
            NBLOCKSTARTS(2)=II+1
            IFDET=1
         ELSEIF(TBLOCK) THEN 
            CALL GNDTS_BLK(NEL,NHG,BRR,NBASISMAX,NMRKS, .FALSE.,
     &         NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS2,TPARITY,
     &           SymRestrict,IFDET,TMAT,.NOT.TREAD,NDETTOT,BLOCKSYM,
     &            TCSF)
         ELSEIF(TCSF) THEN
            CALL GNCSFS(NEL,NHG,BRR,NBASISMAX,
     &         NMRKS,.FALSE.,G1,TSPN,LMS2,TPARITY,
     &         SymRestrict,NDET,IFDET,.FALSE.,0,0,.FALSE.,0)
               NBLOCKSTARTS(1)=1
               NBLOCKSTARTS(2)=II+1
         ELSE
            CALL GNDTS(NEL,NHG,BRR,NBASISMAX,NMRKS,.FALSE.,
     &         G1,TSPN,LMS,TPARITY,SymRestrict,II,IFDET)
               NBLOCKSTARTS(1)=1
               NBLOCKSTARTS(2)=II+1
         ENDIF
         OPEN(8,FILE='DETS',STATUS='UNKNOWN')
         DO I=1,NDET
            CALL WRITEDET(8,NMRKS(1,I),NEL,.FALSE.)
            CALL GETSYM(NMRKS(1,I),NEL,G1,NBASISMAX,ISYM)
            CALL WRITESYM(8,ISym%Sym,.TRUE.)
         ENDDO
         CLOSE(8)
C..
C.. Now generate the fermi determiant
C.. Work out the fermi det
         DO I=1,NEL
            FDET(I)=NMRKS(I,IFDET)
         ENDDO
         WRITE(6,*) "Fermi Determinant:",IFDET

         WRITE(6,*) ' NUMBER OF SYMMETRY UNIQUE DETS ' , NDET

C         WRITE(6,*) ' TOTAL NUMBER OF DETS.' , NDETTOT
         IF(NEVAL.EQ.0) THEN
            WRITE(6,*) 'NEVAL=0.  Setting NEVAL=NDET'
            NEVAL=NDET
         ENDIF
         IF(NEVAL.GT.NDET) THEN
            WRITE(6,*) 'NEVAL>NDET.  Setting NEVAL=NDET'
            NEVAL=NDET
         ENDIF

         IF(NPATHS.EQ.-1) THEN
            WRITE(6,*) 'NPATHS=-1.  SETTING NPATHS to NDET'
            NPATHS=NDET
         ENDIF

         IF(ABS(DETINV).GT.NDET) THEN
            WRITE(6,*) 'DETINV=',DETINV,'>NDET=',NDET
            WRITE(6,*) 'Setting DETINV and NPATHS to 0'
            DETINV=0
            NPATHS=0
         ENDIF
         IF(ABS(DETINV)+NPATHS.GT.NDET) THEN
            WRITE(6,*) 'DETINV+NPATHS=',ABS(DETINV)+NPATHS,'>NDET=',NDET
            WRITE(6,*) 'Setting DETINV and NPATHS to 0'
            DETINV=0
            NPATHS=0
         ENDIF
         CALL FLUSH(6)
C ==----------------------------------------------------------------==
C..Set up memory for c's, nrow and the label
         IF(TCALCHMAT) THEN
            WRITE(6,*) "CK Size",NDET*NEVAL*HElementSize
            CALL MEMORY(IP_CKN,NDET*NEVAL*HElementSize,'CKN')
            CALL AZZERO(CKN,NDET*NEVAL*HElementSize)
C..
            CALL MEMORY(IP_CK,NDET*NEVAL*HElementSize,'CK')
            CALL AZZERO(CK,NDET*NEVAL*HElementSize)
C..
            CALL MEMORY(IP_W,NEVAL,'W')
            CALL AZZERO(W,NEVAL)
         ENDIF
         IF(TREADRHO.AND..NOT.TREAD) THEN
            WRITE(10,*)
     &         "TREADRHO specified, but not TREAD.  Setting TREAD=.T."
            TREAD=.TRUE.
         ENDIF
C..
         IF(TREAD) THEN
           CALL READ_PSI(BOX,BOA,COA,NDET,NEVAL,
     &             NBASISMAX,NEL,CK,W)
            GOTO 100
         ENDIF
         IF(TREADRHO) GOTO 100
      ENDIF
C.. we need to calculate a value for RHOEPS, so we approximate that
C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this
C.. we have put TMAT instead of ZIA
      IF(I_HMAX.NE.-20) THEN
C.. If we're using rhos,
         RHOEPS=GETRHOEPS(RHOEPSILON,BETA,NEL,NBASISMAX,G1,NHG,BRR,
     &               NMSH,FCK,TMAT,NMAX,ALAT,UMAT,I_P,ECORE)

         WRITE(6,*) "RHOEPS:",RHOEPS
      ELSE
C.. we're acutally diagonalizing H's, so we just leave RHOEPS as RHOEPSILON
         RHOEPS=RHOEPSILON
      ENDIF





      TMC=TCALCHMAT.AND.(.NOT.TENERGY)
C.. now back to the storing H
      IF(TCALCHMAT) THEN
         WRITE(6,*) "Calculating H matrix"
C..We need to measure HAMIL and LAB first 
         CALL MEMORY(IP_NROW,NDET,'NROW')
         CALL IAZZERO(NROW,NDET)
         ICMAX=1
         CALL DETHAM(NDET,NEL,NMRKS,NBASISMAX,
     &         NHG,HAMIL,G1,LAB,NROW,.TRUE.,
     &         NMSH,FCK,TMAT,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BRR)
         WRITE(6,*) ' FINISHED COUNTING '
         CALL FLUSH(6)
C..Now we know size, allocate memory to HAMIL and LAB
         LENHAMIL=GC
         CALL MEMORY(IP_HAMIL,LENHAMIL*HElementSize,'HAMIL')
         CALL AZZERO(HAMIL,LENHAMIL*HElementSize)
C..
         CALL MEMORY(IP_LAB,LENHAMIL/IRAT+1,'LAB')
         CALL IAZZERO(LAB,LENHAMIL)
C..Now we store HAMIL and LAB 
         CALL DETHAM(NDET,NEL,NMRKS,NBASISMAX,
     &      NHG,HAMIL,G1,LAB,NROW,.FALSE.,
     &      NMSH,FCK,TMAT,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BRR)

         IF(BTEST(ILOGGING,7)) THEN
C.. we write out H now
            OPEN(8,FILE='HAMIL',STATUS='UNKNOWN')
            J=0
            JR=0
C            HMAX=-dflmax()
C            HMIN=dflmax()
            DO I=1,LENHAMIL
               DO WHILE(I.GT.J)
                  JR=JR+1
                  J=J+NROW(JR)
               ENDDO
               WRITE(8,"(2I12,$)") JR,LAB(I)
               IF(HElementSize.EQ.1) THEN
                  WRITE(8,*) HAMIL(I)
               ELSE
                  WRITE(8,*) HAMIL(I),ABS(HAMIL(I))
               ENDIF
C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
C               WRITE(14,"(A,$)"),"|"
C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
C              WRITE(14,"(F27.20)") HAMIL(I)
C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
C               WRITE(14,"(A,$)"),"|"
C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
C               WRITE(14,"(F27.20)") HAMIL(I)

C               IF(HAMIL(I).GT.HMAX) HMAX=HAMIL(I)
C               IF(HAMIL(I).LT.HMIN) HMIN=HAMIL(I)
            ENDDO
            CLOSE(8)
         ENDIF
        WRITE(6,*) '<D0|H|D0>=',
     &            GETHELEMENT(IFDET,IFDET,HAMIL,LAB,NROW,NDET)
        WRITE(6,*) '<D0|T|D0>=',CALCT(NMRKS(1:NEL,IFDET),NEL,
     &     G1,TMAT,NBASIS) 
CC         CALL HAMHIST(HMIN,HMAX,LENHAMIL,NHISTBOXES)
      ENDIF
C.. We've now finished calculating H if we were going to.
C.. IF ENERGY CALC (for which we need to have calced H)

      IF(TENERGY) THEN
         IF(NBLK.NE.0) THEN
C..Things needed for Friesner-Pollard diagonalisation
            IF(TMC) STOP 'TMC and TENERGY set - Stopping'
            IF(HElementSize.NE.1)
     &         STOP 'Cannot do Lanczos on Complex orbitals.'
            NKRY1=NKRY+1
            NBLOCK=MIN(NEVAL,NBLK)
            LSCR=MAX(NDET*NEVAL,8*NBLOCK*NKRY)
            LISCR=6*NBLOCK*NKRY
C..
            WRITE(*,'(/,/,8X,64(1H*))')
            WRITE(*,'(7X," *",62X,"*")')
          WRITE(*,'(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
            WRITE(*,'(7X," *",62X,"*")')
            WRITE(*,'(7X,1X,64(1H*))')
C..Set up memory for FRSBLKH
            CALL MEMORY(IP_A,NEVAL*NEVAL,'A')
            CALL AZZERO(A,NEVAL*NEVAL)
C..
C,, W is now allocated with CK
C..
            CALL MEMORY(IP_V,NDET*NBLOCK*NKRY1,'V')
            CALL AZZERO(V,NDET*NBLOCK*NKRY1)
C..   
            CALL MEMORY(IP_AM,NBLOCK*NBLOCK*NKRY1,'AM')
            CALL AZZERO(AM,NBLOCK*NBLOCK*NKRY1)
C..
            CALL MEMORY(IP_BM,NBLOCK*NBLOCK*NKRY,'BM')
            CALL AZZERO(BM,NBLOCK*NBLOCK*NKRY)
C..
            CALL MEMORY(IP_T,3*NBLOCK*NKRY*NBLOCK*NKRY,'T')
            CALL AZZERO(T,3*NBLOCK*NKRY*NBLOCK*NKRY)
C..
            CALL MEMORY(IP_WT,NBLOCK*NKRY,'WT')
            CALL AZZERO(WT,NBLOCK*NKRY)
C..
            CALL MEMORY(IP_SCR,LSCR,'SCR')
            CALL AZZERO(SCR,LSCR)
            CALL MEMORY(IP_ISCR,LISCR,'ISCR')
            CALL IAZZERO(ISCR,LISCR)
            CALL MEMORY(IP_INDEX,NEVAL,'INDEX')
            CALL IAZZERO(INDEX,NEVAL)
C..
            CALL MEMORY(IP_WH,NDET,'WH')
            CALL AZZERO(WH,NDET)
            CALL MEMORY(IP_WORK2,3*NDET,'WORK2')
            CALL AZZERO(WORK2,3*NDET)
            CALL MEMORY(IP_V2,NDET*NEVAL,'V2')
            CALL AZZERO(V2,NDET*NEVAL)
C..Lanczos iterative diagonalising routine
            CALL FRSBLKH(NDET,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,
     &  NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT,
     &  SCR,ISCR,INDEX,WH,WORK2,V2,NCYCLE,B2L)
            CALL DSCAL(NEVAL,-1.D0,W,1)
         ELSE
C.. We splice in a non-Lanczos diagonalisin routine if NBLOCK=0
            IF(NEVAL.NE.NDET) THEN
               WRITE(6,*) "NEVAL.NE.NDET.",NEVAL,NDET,
     &           " Cannot exactly diagonalize."
               STOP
            ENDIF
            WRITE(6,*) "NBLK=0.  Doing exact diagonalization."
            IF(TCALCHMAT) THEN
               CALL MEMORY(IP_WORK,4*NDET*HElementSize,'WORK')
               CALL MEMORY(IP_WORK2,3*NDET,'WORK2')
               CALL HDIAG(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,LENHAMIL,
     &                  NBLOCKSTARTS,NBLOCKS,BLOCKSYM)
            ELSE
               CALL HDIAG_NH(NDET,NBLOCKSTARTS,NBLOCKS,
     &            NEL,NMRKS,NBASISMAX,NBASIS,G1,NMSH,BRR,TMAT,
     &            FCK,TMAT,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,
     &            BETA,I_P,ILOGGING,IFDET,ARR,BLOCKSYM)
C.. We're not storing the energies, so we pretend we weren't asked for
C.. them
               TENERGY=.FALSE.
            ENDIF
         ENDIF
C..
         CALL MEMORY(IP_TKE,NEVAL,'TKE')
C.. END ENERGY CALC
      ENDIF
      IF(TENERGY) THEN
         EXEN=CALCMCEN(NDET,NEVAL,CK,W,BETA,0.D0)
         WRITE(6,"(A,F19.9)") "EXACT E(BETA)=",EXEN
         GSEN=CALCDLWDB(IFDET,NDET,NEVAL,CK,W,BETA,0.D0)
         WRITE(6,"(A,F19.9)") "EXACT DLWDB(D0)=",GSEN
         WRITE(6,"(A,F19.9)") "GROUND E=",W(1)
      ENDIF

      call FLUSH(6)
C.. If we're calculating rhos (for which we have to have calced H
      IF(TRHOIJ) THEN
         IF((.NOT.TENERGY).AND.(.NOT.TREADRHO)) THEN
            WRITE(6,*) "Calculating approx RHOs"
            WRITE(6,*) "Using Trotter decomposition? ",TTROT
            WRITE(6,*) "Order of Taylor: ",ABS(NTAY)
            CALL CALCAPPROXRHOIJ(BETA,I_P,HAMIL,LAB,NROW,NDET,
     &            RHOMIN,RHOMAX,NRHOS,RHOEPS,TTROT,NTAY)
         ENDIF
      ENDIF
C..Free HAMIL AND LAB memory if we no longer need them
      IF(TCALCHMAT.AND..NOT.(TMONTE.AND.TMC)) THEN
         CALL FREEM(IP_HAMIL)
         CALL FREEM(IP_LAB)
      ENDIF

C.. IF ENERGY CALC
      IF (TENERGY) THEN
C..
         IF(.NOT.TCSF) THEN
            CALL CFF_CHCK(NDET,NEVAL,NMRKS,NBASISMAX,NEL,G1,CK,
     &         ALAT,TKE,TMAT,NHG,ILOGGING)
         ELSE
            DO I=1,NEVAL
               TKE(I)=0.D0
            ENDDO 
         ENDIF
         IF(BTEST(ILOGGING,7))
     &     CALL WRITE_PSI(BOX,BOA,COA,NDET,NEVAL,
     &      NBASISMAX,NEL,CK,W)
         IF(BTEST(ILOGGING,8))
     &     CALL WRITE_PSI_COMP(BOX,BOA,COA,NDET,NEVAL,
     &      NBASISMAX,NEL,CK,W)
         WRITE(6,*) 
     &  '       ==--------------------------------------------------== '
         WRITE(6,'(A5,5X,A15,1X,A18,1x,A20)') 'STATE',
     &  'KINETIC ENERGY', 'COULOMB ENERGY', 'TOTAL ENERGY'
         OPEN(15,FILE='ENERGIES',STATUS='UNKNOWN')
         DO IN=1,NEVAL
            WRITE(6,'(I5,2X,3(F19.11,2x))')
     &          IN,TKE(IN),W(IN)-TKE(IN),W(IN)
            WRITE(15,"(I7,$)") IN
            CALL WRITEDET(15,NMRKS(1,IN),NEL,.FALSE.)
            WRITE(15,"(F19.11)") W(IN)
         ENDDO
         CLOSE(15)
         WRITE(6,*) 
     &  '       ==--------------------------------------------------== '
C., END energy calc
      ENDIF

C.. Jump to here if just read Psi in
100   CONTINUE


C.. we need to calculate a value for RHOEPS, so we approximate that
C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this 
C.. If we haven't already calced RHOEPS, do it now
      IF(TENERGY) THEN
         RHOEPS=RHOEPSILON*EXP(-BETA*(W(1))/I_P)
         WRITE(6,*) "RHOEPS:",RHOEPS
         IF(TREAD) THEN
            EXEN=CALCMCEN(NDET,NEVAL,CK,W,BETA,0.D0)
            WRITE(6,"(A,F19.5)") "EXACT E(BETA)=",EXEN
            GSEN=CALCDLWDB(1,NDET,NEVAL,CK,W,BETA,0.D0)
            WRITE(6,"(A,F19.5)") "EXACT DLWDB(D0)=",GSEN
         ENDIF
         OPEN(14,FILE='RHOPIIex',STATUS='UNKNOWN')
         IF(NDETWORK.EQ.0.OR.NDETWORK.GT.NDET) NDETWORK=NDET
         DO III=1,NDETWORK
            CALL CALCRHOPII(III,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,
     &            0.D0,FLRI,FLSI,TWARN)
            IF(TWARN) THEN
               IF(III.EQ.1) THEN
                  WRITE(6,*) "Warning received from CALCRHOPII."
                  IF(TREAD) THEN
                     WRITE(6,*) "TREAD set. Cannot calculate RHOII."
                  ELSE
                  WRITE(6,*) "Calculating RHOII using 1st order Taylor."
                  ENDIF
               ENDIF
               IF(.NOT.TREAD) THEN
               CALL CALCRHO2(NMRKS(1,III),NMRKS(1,III),BETA,I_P,NEL,
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,
     &            UMAT,RH,1,0,ECORE)
C               WRITE(6,*) RH
               FLRI=LOG(RH)
               FLSI=FLSI-I_P*FLRI
               ENDIF
            ENDIF
            CALL WRITEDET(14,NMRKS(1,III),NEL,.FALSE.)
            GSEN=CALCDLWDB(III,NDET,NEVAL,CK,W,BETA,0.D0)
            CALL GETSYM(NMRKS(1,III),NEL,G1,NBASISMAX,ISYM)
            CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
            WRITE(14,"(4G25.16,I5)") 
     &         EXP(FLSI+I_P*FLRI),FLRI*I_P,FLSI,GSEN,IDEG
         ENDDO
C         CLOSE(17)
         CLOSE(14)
      ENDIF

C      DBRAT=0.001
C      DBETA=DBRAT*BETA
      WRITE(6,*) "DBETA=",DBETA

      IF(.NOT.TREAD) THEN
         WRITE(6,*) '<D0|H|D0>=',GETHELEMENT2(FDET,FDET,NEL,NBASISMAX,
     &             G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,IC,ECORE)
         WRITE(6,*) '<D0|T|D0>=',CALCT(FDET,NEL,G1,TMAT,NBASIS)
         IF(TUEG) THEN
!  The actual KE rather than the one-electron part of the Hamiltonian
            WRITE(6,*) 'Kinetic=',CALCT2(FDET,NEL,G1,ALAT,NBASIS,CST)
         ENDIF
      ENDIF
      IF(TMCDET) THEN
C.. Generate the determinant from which we start the MC
         NLIST=1
         CALL GENSYMDETSS(MDK,NEL,G1,BRR,NBASIS,MCDET,NLIST,NBASISMAX)
         IF(NLIST.EQ.0) THEN
C.. we couldn't find a det of that symmetry
            STOP 'Cannot find MC start determinant of correct symmetry'
         ENDIF
      ELSE
C         CALL GENRANDOMDET(NEL,NBASIS,MCDET)
         DO I=1,NEL
            MCDET(I)=FDET(I)
         ENDDO
      ENDIF
      IF(TMONTE) THEN
         WRITE(6,"(A,$)") 'MC Start Det: '
         CALL WRITEDET(6,MCDET,NEL,.TRUE.)
      ENDIF

      IF(TRHOIJ) THEN
C.. Approx RHOS have been generated above, and written to the file RHOS
C.. Exact RHOS are generated below

C..   First allocate memory to store our labels in our rho matrix
C..   This will store the index of each row within RHOIJ, the last element
C..   therefore being the total number of elements in RHOIJ
         WRITE(6,*) "CREATING RHO MATRIX"
         CALL MEMORY(IP_RHOIJROWS,(NDET+1)/IRAT+1,'RHOIJROWS')
         CALL IAZZERO(RHOIJROWS,NDET+1)

C.. Calculate RHO_IJ as Sum_K <I|K><K|J>exp(-(Beta/p)epsilon_K
         IF(TENERGY.AND.(.NOT.TREADRHO)) THEN
            WRITE(6,*) "Calculating exact RHOs"
            CALL CALCRHOIJ(NDET,NEVAL,CK,W,BETA,I_P,RHOMIN,RHOMAX,NRHOS,
     &               RHOEPS,0.D0)
C.. This writes the RHOS to the file RHOS
         ENDIF

         IF(TREADRHO) THEN
C.. If we're reading in RHOs instead of calculating, we need to know
C.. how many there are.
            OPEN(10,FILE="NRHOS",STATUS="OLD")
            READ(10,*) NRHOS
            CLOSE(10)
            WRITE(6,*) "NRHOS:",NRHOS
         ENDIF
         CALL MEMORY(IP_RHOIJLAB,NRHOS/IRAT+1,'RHOIJLAB')
         CALL IAZZERO(RHOIJLAB,NRHOS)

         CALL MEMORY(IP_RHOIJ,NRHOS,'RHOIJ')
         CALL AZZERO(RHOIJ,NRHOS)

C.. Read in the RHOS from the file RHOS
         CALL STORERHOIJ(NDET,NRHOS,RHOIJROWS,RHOIJ,RHOIJLAB,
     &            NRHOMAXWIDTH)
C               CALL MCPATHSR(482,I_P,NDET,RHOIJROWS,RHOIJ,RHOIJLAB,
C     &               NRHOMAXWIDTH,I_HMAX,NWHTAY,FLNRII,FLNSUM)

C.. Init the Chebyshev integration tables for the Zsums MMAX=511
C         IF((NPATHS.GT.0.OR.TMONTE).AND.(NWHTAY.EQ.-1)) CALL INITAL(511)

         
C.. NB - still in TRHOIJ
         IF(NPATHS.GT.0.OR.DETINV.GT.0) THEN
            WRITE(6,*) "Calculating ",NPATHS," W_Is..."
            WRITE(6,*) "I_HMAX:",I_HMAX
            IF(I_HMAX.EQ.-10) THEN
               OPEN(11,FILE="MCSUMMARY",STATUS="UNKNOWN")
               WRITE(11,*) "Calculating ",NPATHS," W_Is..."
               CLOSE(11)
            ELSE
               OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
               WRITE(11,*) "Calculating ",NPATHS," W_Is..."
               CLOSE(11)
            ENDIF
             IF(I_VMAX.EQ.0.AND.I_HMAX.NE.0) THEN
               WRITE(6,*) "I_VMAX=0.  Using paths limited by hop number"
               OPEN(15,FILE='RHOPII',STATUS='UNKNOWN')
               DO III=1,NPATHS
C.. we use the old MCPATHSR, which doesn't do vertices etc.
                  CALL MCPATHSR(III,I_P,NDET,RHOIJROWS,RHOIJ,RHOIJLAB,
     &               NRHOMAXWIDTH,I_HMAX,NWHTAY,ILOGGING,0.D0,
     &               FLNRII,FLNSUM)
                  WRITE(15,"(I12,$)") III
                  CALL WRITEDET(15,NMRKS(1,III),NEL,.FALSE.)
                  WRITE(15,"(A,3G25.16)") " ",
     &                EXP(FLNSUM+I_P*FLNRII),I_P*FLNRII,FLNSUM
               ENDDO
            ELSE
               WRITE(6,*) "I_VMAX<>0.  Using paths limited by ",I_VMAX,
     &                  " vertices."
C..     NTAY<0 => Exact new method
C..     FCK->RHOIJ
C..     ZIA->RHOIJLAB
C..     UMAT->RHOIJROWS
C..     NMSH->NMRKS
C..     BETA->NDET
               IF(TNPDERIV)
     &            WRITE(6,*) " *** CANNOT CALCULATE RHOPII derivs **"
               IF(NTAY.GE.0) THEN
                 WRITE(6,*) "NTAY must be <0 to use NPATHS and TRHOIJ ",
     &           "and TENERGY together.  Setting NTAY=-NEL."
                  NTAY=-NEL
               ENDIF
CCC..AAAA
CC.. Instead of NMAX, we put ARR
               CALL CALCRHOPII2(NMRKS,NDET,I_P,I_HMAX,I_VMAX,NEL,NDET,
     &           NBASISMAX,G1,NHG,BRR,NMRKS,RHOIJ,RHOIJLAB,ARR,ALAT,
     &           RHOIJROWS,NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,
     &            .FALSE.,0.D0,DETINV,TSPECDET,SPECDET)
C               CALL MCPATHSR3(NMRKS(1:NEL,III),NHG,I_P,I_HMAX,I_VMAX,
C     &               NEL,NBASISMAX,G1,NHG,BRR,NMRKS,RHOIJ,RHOIJLAB,
C     &               NMAX,ALAT,RHOIJROWS,NTAY,RHOEPS,
C     &               LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,FLNRII,FLNSUM)
            ENDIF
            
         ENDIF
C.. Still in IF(TRHOIJ)
         IF(TMONTE) THEN
            CALL MEMORY(IP_RHOIJ1,NRHOS,'RHOIJ1')
            CALL AZZERO(RHOIJ1,NRHOS)
            CALL MEMORY(IP_RHOIJ2,NRHOS,'RHOIJ2')
            CALL AZZERO(RHOIJ2,NRHOS)
            WRITE(6,*) "Calculating RHOS for BETA+-DBETA"
            IF(TMC) THEN
               WRITE(6,*) "Approx RHOS.  NTAY=",NTAY
               CALL CALCSTOREAPRHOIJ(NDET,BETA-DBETA,I_P,HAMIL,LAB,
     &            NROW,RHOIJROWS,RHOIJ1,RHOIJLAB,TTROT,NTAY,0.D0)
               CALL CALCSTOREAPRHOIJ(NDET,BETA+DBETA,I_P,HAMIL,LAB,
     &            NROW,RHOIJROWS,RHOIJ2,RHOIJLAB,TTROT,NTAY,0.D0)
            ELSEIF(TENERGY) THEN
               WRITE(6,*) "Exact RHOS"
               CALL CALCSTOREEXRHOIJ(NDET,BETA-DBETA,I_P,NEVAL,CK,
     &            W,0.D0,RHOIJROWS,RHOIJ1,RHOIJLAB)
               CALL CALCSTOREEXRHOIJ(NDET,BETA+DBETA,I_P,NEVAL,CK,
     &            W,0.D0,RHOIJROWS,RHOIJ2,RHOIJLAB)
            ELSE
               STOP "Neither approx nor exact. Cannot calc rhos"
            ENDIF
            WRITE(6,*) "I_HMAX:",I_HMAX
            WRITE(6,*) "Calculating MC Energy..."
            call FLUSH(6)
            IF(I_VMAX.EQ.0.AND.I_HMAX.NE.0) THEN
               WRITE(6,*) "I_VMAX=0.  Using paths limited by hops"
C.. NTAY>=0 uses old nonvertex based MCPATHSR
               EN=DMONTECARLO(I_P,NDET,RHOIJROWS,RHOIJ,RHOIJ1,RHOIJ2,
     &            BETA,DBETA,RHOIJLAB,NRHOMAXWIDTH,I_HMAX,IMCSTEPS,
     &            G1,NMRKS,NEL,NBASISMAX,NHG,BRR,IEQSTEPS,NWHTAY,
     &            ILOGGING,0.D0)
            ELSE
               WRITE(6,*) "I_VMAX<>0.  Using paths limited by ",I_VMAX,
     &            " vertices."
C..     NTAY<0 => new method, sending the RHOIJ matrices along
C..     FCK->RHOIJ
C..     ZIA->RHOIJLAB
C..     UMAT->RHOIJROWS
C..     NMSH->NMRKS
C..     NMAX->RHOIJ1
C..     ALAT->RHOIJ2
C..     BETA->NDET
               EN=DMONTECARLO2(MCDET,I_P,NDET,DBETA,I_HMAX,I_VMAX,
     &               IMCSTEPS,G1,NEL,NBASISMAX,NHG,BRR,IEQSTEPS,
     &               NMRKS,RHOIJ,RHOIJLAB,RHOIJ1,RHOIJ2,RHOIJROWS,NTAY,
     &               RHOEPS,NWHTAY,ILOGGING,ECORE,BETAEQ)
            ENDIF
            WRITE(6,*) "MC Energy=",EN
            IF(TMC) THEN
C.. Now we can get rid of the hamiltonian
               CALL FREEM(IP_HAMIL)
               CALL FREEM(IP_LAB)
            ENDIF
            CALL FREEM(IP_RHOIJ1)
            CALL FREEM(IP_RHOIJ2)
         ENDIF
C         WRITE(6,*) "Creating RHO Histogram"
C         CALL RHOHIST(RHOMIN,RHOMAX,NRHOS,NHISTBOXES)
         CALL FREEM(IP_RHOIJ)
         CALL FREEM(IP_RHOIJLAB)
         CALL FREEM(IP_RHOIJROWS)        
      ENDIF
C..   ENDIF(TRHOIJ)
      IF(TMONTE.AND.TENERGY.AND.NTAY.EQ.-1) THEN
         WRITE(6,*) "Calculating Exact MC Energy..."
         EN=DOEXMC(NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,0.D0,
     &      IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,NHG,BRR,IEQSTEPS)
      ENDIF
      IF(TBEGRAPH) THEN
         IF(TENERGY) THEN
            IF(NTAY.NE.0) THEN
               CALL DOBEGRAPH(NDET,NEVAL,CK,W,I_P,ILOGGING,
     &            G1,NMRKS,NEL,NBASISMAX,NHG,BRR)
            ELSE
C.. NTAY=0 signifying we're going to calculate the RHO values when we
C.. need them from the list of eigenvalues.   
C.. Hide NMSH=NEVAL
C..      FCK=W
C..      ZIA=CK
C..      UMAT=NDET
C..      ALAT=NMRKS        
               CALL DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,
     &            NBASISMAX,G1,NHG,BRR,NEVAL,W,CK,NMAX,NMRKS,NDET,
     &            NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
            ENDIF
         ELSE
            CALL DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,
     &            NBASISMAX,G1,NHG,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,
     &            NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
         ENDIF
      ENDIF
C.. We haven't calculated the RHOIJ matrix but can still do MC
      IF(TRHOIJND.OR..NOT.TRHOIJ) THEN 
         IF(NPATHS.NE.0.OR.DETINV.GT.0) THEN
            IF(TRHOIJND) THEN
C.. We're calculating the RHOs for interest's sake, and writing them,
C.. but not keeping them in memory
               WRITE(6,*) "Calculating RHOS..."
               WRITE(6,*) "Using approx NTAY=",NTAY
               CALL CALCRHOSD(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,
     &            NBASISMAX,G1,NHG,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,
     &            NTAY,RHOEPS,NWHTAY,ECORE)
            ENDIF
            WRITE(6,*) "Calculating ",NPATHS," W_Is..."
            IF(BTEST(ILOGGING,1)) THEN
               IF(I_HMAX.EQ.-10) THEN
                  OPEN(11,FILE="MCSUMMARY",STATUS="UNKNOWN")
                  WRITE(11,*) "Calculating ",NPATHS," W_Is..."
                  CLOSE(11)
               ELSE
                  OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
                  WRITE(11,*) "Calculating ",NPATHS," W_Is..."
                  CLOSE(11)
               ENDIF
            ENDIF
            IF(NTAY.EQ.0.AND.TENERGY) THEN
               WRITE(6,*) "Using exact RHOs generated on the fly"
C.. we've calculated energies, and we're passing them through to
C.. calculate the exact RHOS
C.. NTAY=0 signifying we're going to calculate the RHO values when we
C.. need them from the list of eigenvalues.   
C.. Hide NMSH=NEVAL
C..      FCK=W
C..      ZIA=CK
C..      UMAT=NDET
C..      ALAT=NMRKS
C..      NMAX=ARR
              CALL CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,
     &           NBASISMAX,G1,NHG,BRR,NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,
     &            DETINV,TSPECDET,SPECDET)

            ELSEIF(NTAY.GT.0) THEN
               WRITE(6,*) "Using list of determinants."
               WRITE(6,*) "Using approx RHOs generated on the fly,",
     &            "NTAY=",NTAY
C.. we haven't calculated the energy, so we're calculating the weights
C.. with approx RHOs
               IF(TENERGY) THEN
C.. If we've generated a list of dets
C.. Instead of NMAX, we put ARR
              CALL CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,
     &           NBASISMAX,G1,NHG,BRR,NMSH,FCK,TMAT,ARR,ALAT,UMAT,NTAY,
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,
     &            DETINV,TSPECDET,SPECDET)
               ELSE
                  IF(TCSF) THEN
                     IF(.NOT.TSPECDET) THEN
                        WRITE(6,*) "SPECDET not specified. Using Fermi",
     &                     " determinant ONLY"
                        TSPECDET=.TRUE.
                        CALL ICOPY(NEL,FDET,1,SPECDET,1)
                     ENDIF
                  ENDIF
C.. Instead of NMAX we have ARR
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
              CALL CALCRHOPII3(BETA,I_P,I_HMAX,I_VMAX,NEL,
     &           NBASISMAX,G1,NHG,BRR,NMSH,FCK,TMAT,ARR,ALAT,UMAT,NTAY,
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,
     &            DETINV,TSPN,LMS2,TPARITY,SymRestrict,TSPECDET,SPECDET)
               ENDIF
            ELSE
               WRITE(6,*)
     &  "Invalid combination of NTAY and TENERGY.  No NPATHS calculated"
            ENDIF
         ENDIF
         IF(TMONTE) THEN
C            DBRAT=0.01
C            DBETA=DBRAT*BETA
            WRITE(6,*) "I_HMAX:",I_HMAX
            WRITE(6,*) "Calculating MC Energy..."
            CALL FLUSH(6)
            IF(NTAY.GT.0) THEN
               WRITE(6,*) "Using approx RHOs generated on the fly,",
     &            "NTAY=",NTAY
C.. NMAX is now ARR
               EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,
     &               IMCSTEPS,
     &            G1,NEL,NBASISMAX,NHG,BRR,IEQSTEPS,
     &         NMSH,FCK,TMAT,ARR,ALAT,UMAT,NTAY,RHOEPS,NWHTAY,ILOGGING,
     &            ECORE,BETAEQ)
            ELSEIF(NTAY.EQ.0) THEN
               IF(TENERGY) THEN
                  WRITE(6,*) "Using exact RHOs generated on the fly"
C.. NTAY=0 signifying we're going to calculate the RHO values when we
C.. need them from the list of eigenvalues.   
C.. Hide NMSH=NEVAL
C..      FCK=W
C..      ZIA=CK
C..      UMAT=NDET
C..      ALAT=NMRKS
C..      NMAX=ARR
                  EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,
     &                     IMCSTEPS,
     &               G1,NEL,NBASISMAX,NHG,BRR,IEQSTEPS,
     &         NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,RHOEPS,NWHTAY,ILOGGING,
     &               ECORE,BETAEQ)
               ELSE
                  
                  STOP "TENERGY not set, but NTAY=0" 
               ENDIF
            ENDIF
            WRITE(6,*) "MC Energy:",EN
CC            WRITE(12,*) DBRAT,EN
         ENDIF
         
      ENDIF
C.. /AJWT
      IF(TRHOOFR) THEN
      
C..Generate memory for RHO and SITAB
        CALL MEMORY(IP_RHO,NMSH*NMSH*NMSH,'RHO')
        CALL MEMORY(IP_SITAB,NMSH*NMAX,'SITAB')
        CALL MEMORY(IP_SCRTCH,NMSH*NMSH*NMSH,'SCRTCH')
C..Calculate RHOOFR
        CALL RHOOFR(NHG,CK,G1,RHO,NMSH,SITAB,NMAX,NMRKS,NEL,
     &    NDET,NEVAL,RS,ALAT,OMEGA,SCRTCH)
C..
        CALL MEMORY(IP_DLINE,NMSH,'DLINE')
        CALL AZZERO(DLINE,NMSH)
C..Calculate RHOOFR in certain directions
C..001
        CALL PLANARAV(RHO,NMSH,DLINE,0,0,1,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV001',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
C..100
        CALL PLANARAV(RHO,NMSH,DLINE,1,0,0,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV100',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
C..010
        CALL PLANARAV(RHO,NMSH,DLINE,0,1,0,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV010',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
        IF(TCORR) THEN
C..Now generate memory for XCHOLE
          CALL MEMORY(IP_XCHOLE,NMSH*NMSH*NMSH,'XCHOLE')
          CALL MEMORY(IP_PSIR,2*NMSH+1,'PSIR')
          CALL AZZERO(PSIR,2*NMSH+1)
C..
          IXD=1
          IYD=0
          IZD=0
          SPAC=0.D0
          CALL GEN_XCHOLE(CK,PSIR,IOBS,JOBS,KOBS,G1,SITAB,
     &          NMAX,NMSH,NHG,IXD,IYD,IZD,RHO,.TRUE.,
     &          XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,NEL)
          CALL WRITE_RHO(10,'COMPXCHOLE',XCHOLE,NMSH,NMSH,NMSH,
     &          ALAT,.FALSE.,.TRUE.,RS)
C..
          CALL XCHOLES(CK,PSIR,IOBS,JOBS,KOBS,G1,SITAB,
     &          NMAX,NMSH,NHG,RHO,XCHOLE,SPAC,RS,ALAT,DLINE,
     &          OMEGA,NMRKS,NDET,NEL,NEVAL)
C          OPEN(11,FILE='XCEDENS',STATUS='UNKNOWN')
C          DO IZ=1,NMSH
C            CALL XCENERGY_DENSITY(CK,PSIR,DLINE,IOBS,JOBS,IZ,G1,
C     &       SITAB,NMAX,NMSH,NHG,RHO,XCHOLE,ALAT,OMEGA,
C     &       NMRKS,NDET,NEVAL,NEL)
C          ENDDO
C          CLOSE(11)
        ENDIF
      ENDIF
      IF(TFODM) THEN
        CALL MEMORY(IP_SUMA,NMAXX*NMAXY*NMAXZ,'SUMA')
        CALL AZZERO(SUMA,NMAXX*NMAXY*NMAXZ)
        CALL FODMAT(NEL,NHG,NDET,NEVAL,ISTATE,NMRKS,G1,CK,
     &          NMAXX,NMAXY,NMAXZ,SUMA)
      ENDIF

C..
      CALL ENDSYM()
      CALL DESTROYUMATCACHE()
      CALL MEMORY_CHECK
C ==-------------------------------------------------------------------==
      CALL TIHALT('NECICUBE  ',ISUB)
C ==-------------------------------------------------------------------==
      CALL PMAXMEM
      CALL TIPRI
      RETURN
      END
C ==-------------------------------------------------------------------==
      SUBROUTINE CFF_CHCK(NDET,NEVAL,NM,NBASISMAX,NEL,G1,CG,
&	ALAT,TKE,TMAT,NHG,ILOGGING)
      USE HElement
      IMPLICIT NONE
      TYPE(HElement) CG(NDET,NEVAL)
      INTEGER NM(NEL,*),NDET,NEL,NEVAL
      REAL*8 ALAT(3),TKE(NEVAL)
      INTEGER NBASISMAX(3,2),NHG,ILOGGING
      INCLUDE 'basis.inc'
      TYPE(BASISFN) G1(*)
      CHARACTER*255 STR
      TYPE(HElement) TMAT(NHG,NHG)
      REAL*8 PI,S,SUM
      TYPE(HDElement) AUX
      INTEGER I,J,IN,IEL,L
C..Calculate the expectation value of the kinetic energy
C..<Psi|T|Psi>
      PI=ACOS(-1.D0)
      DO IN=1,NEVAL
        TKE(IN)=0.D0
        DO I=1,NDET
          SUM=0.D0
          DO J=1,NEL
            AUX=TMAT(NM(J,I),NM(J,I))
C((ALAT(1)**2)*((G1(1,NM(J,I))**2)/(ALAT(1)**2)+
C     &		(G1(2,NM(J,I))**2)/(ALAT(2)**2)+
C     &		(G1(3,NM(J,I))**2)/(ALAT(3)**2)))
            SUM=SUM+DREAL(AUX)
          ENDDO
C..Cube multiplier
C          CST=PI*PI/(2.D0*ALAT(1)*ALAT(1))
C.. Deal with the UEG
C          IF(NBASISMAX(1,1).LE.0) CST=CST*4.D0
C          SUM=CST*SUM 
          TKE(IN)=TKE(IN)+SUM*SQ(CG(I,IN))
        ENDDO
      ENDDO
C ==--------------------------------------------------------------==
      IF(.F.) THEN
C      IF(BTEST(ILOGGING,7)) THEN
      OPEN(10,FILE='PSI',STATUS='UNKNOWN')
      DO J=1,NEVAL
        IF(J.EQ.1) THEN
          WRITE(10,*) ' GROUND STATE COEFFICIENTS  ' 
        ELSE
          WRITE(10,*) ' COEFFICIENTS FOR EXCITED STATE NUMBER : ' , J
        ENDIF
        S=0.D0
        DO I=1,NDET
         IF(CG(I,J).AGT.1.D-15) THEN
            DO IEL=1,NEL
               WRITE(10,"(I3,I3,2I3,2X,$)") (G1(NM(1,IEL))%K(L),L=1,5)
            ENDDO
            IF(HElementSize.EQ.1) THEN
               WRITE(10,"(F19.9,1X,I7)") CG(I,J),I
            ELSE
               WRITE(10,"(F19.9,1X,I7)") CG(I,J),I
            ENDIF
         ENDIF
         S=S+SQ(CG(I,J))
        ENDDO
        WRITE(10,'(A,F19.5)') ' SQUARE OF COEFFICIENTS : ' , S
        WRITE(10,*)
      ENDDO
      CLOSE(10)
      ENDIF
      RETURN
      END

      SUBROUTINE RHOHIST(RHOMIN,RHOMAX,NRHOS,NHISTBOXES)
         INTEGER IHIST(NHISTBOXES)
         WRITE(6,'(A,F19.5,F19.5)') "MIN,MAX RHO", RHOMIN,RHOMAX
         WRITE(6,'(A,I12)') "TOTAL # RHO ELEMENTS", NRHOS 
         OPEN(1,FILE="RHOS",STATUS="OLD")

C..Read in input parameters
         BOXSPACING=(RHOMAX-RHOMIN)/NHISTBOXES
         DO I=1,NRHOS
            READ(1,*) RH
            IBOX=INT((RH-RHOMIN)/BOXSPACING)+1
C.. This should only happen for RHOMAX
            IF(IBOX>NHISTBOXES) THEN IBOX=NHISTBOXES
            IHIST(IBOX)=IHIST(IBOX)+1
         ENDDO
         CLOSE(1)
         OPEN(1,FILE="Histogram.dat")
         DO I=1,NHISTBOXES
            WRITE(1,"(F19.7,I12)") RHOMIN+BOXSPACING*(I-1),IHIST(I)
         ENDDO
         CLOSE(1)
         RETURN
      END

C.. Note if TWARN becomes set, the RHII sum did not converge.
C.. FLSI will remain usable, but will be equal to log RHII(P), so the 
C.. sum I_P*FLRI+FLSI will still retain the correct value.
      SUBROUTINE CALCRHOPII(I,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ETRIAL,
     &               FLRI,FLSI,TWARN)
         USE HElement
         IMPLICIT NONE
         INTEGER NDET,NEVAL
         TYPE(HElement) CK(NDET,NEVAL)
         REAL*8 W(NEVAL)
         REAL*8 RHII,FLRI,FLSI,ETRIAL,BETA,RH,R
         INTEGER I_P,I,IK,ILOGGING
         LOGICAL LISNAN,TWARN
         RH=0.D0
         RHII=0.D0
         TWARN=.FALSE.
C.. We decompose ln(RHO^(P)_II) = p ln RHO_II + ln sI

C.. First we work out RHO_II
C         WRITE(6,*) BETA,I_P
         DO IK=1,NEVAL
            R=SQ(CK(I,IK))
            R=R*EXP(-(W(IK)-W(1))*BETA/I_P)
            RHII=RHII+R
         ENDDO
         IF(NEVAL.NE.NDET
     &      .AND.EXP(-(W(NEVAL)-W(1))*BETA/I_P).GT.1.D-2) THEN
C.. If we don't have all our eigenvectors and we think our sum has not
C.. converged, we print a warning the first time.
C.. we cannot calculate a proper RHII so we just guess at 1
            IF(I.EQ.1) THEN
             WRITE(6,*) 'WARNING: For Det 1 RHII sum has not converged.'
             WRITE(6,*) ' Setting RHII=1'
            ENDIF
            TWARN=.TRUE.
            RHII=1.D0
            FLRI=0.D0
         ELSE

C.. and Log it
           FLRI=LOG(RHII)-(W(1))*BETA/I_P
         ENDIF
C.. Now we work out RHO^(P)_II/RHO_II^P = sI         
         DO IK=1,NEVAL
            R=SQ(CK(I,IK))
            RH=RH+R*EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         FLSI=LOG(RH)-W(1)*BETA-I_P*FLRI
         IF(LISNAN((RH+1)-RH)) THEN
            RH=0
            FLSI=0
         ENDIF
C         WRITE(17,*) CK(1,I),W(I)
         RETURN
      END

      REAL*8 FUNCTION CALCMCEN(NDET,NEVAL,CK,W,BETA,ETRIAL)
         USE HElement
         IMPLICIT NONE
         INTEGER NDET,NEVAL,IK
         TYPE(HElement) CK(NDET,NEVAL)
         REAL*8  W(NEVAL),BETA,DNORM,EN,ETRIAL
         EN=0.D0
         DNORM=0.D0
         DO IK=1,NEVAL
            EN=EN+(W(IK))*EXP(-(W(IK)-W(1))*BETA)
            DNORM=DNORM+EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         CALCMCEN=EN/DNORM
         RETURN
      END

      REAL*8 FUNCTION CALCDLWDB(I,NDET,NEVAL,CK,W,BETA,ETRIAL)
         USE HElement
         IMPLICIT NONE
         INTEGER NDET,NEVAL,IK,I
         TYPE(HElement) CK(NDET,NEVAL)
         REAL*8  W(NEVAL),BETA,DNORM,EN,ETRIAL
         EN=0.D0
         DNORM=0.D0
         DO IK=1,NEVAL
            EN=EN+SQ(CK(I,IK))*(W(IK))
     &            *EXP(-(W(IK)-W(1))*BETA)
            DNORM=DNORM+SQ(CK(I,IK))*EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         CALCDLWDB=EN/DNORM
         RETURN
      END

      SUBROUTINE LISTEXC(NI,NEL,NBASIS,NBASISMAX)
         INTEGER LIST(NEL,NBASIS*NBASIS*NEL*NEL),NI(NEL)
         INTEGER I,J,ICLIST(NBASIS*NBASIS*NEL*NEL)
            DO J=1,NEL
               WRITE(6,"(I12,$)") NI(J)
            ENDDO
         WRITE(6,*)
         WRITE(6,*)
         CALL GENEXCIT(NI,2,NBASIS,NEL,LIST,ICLIST,NLIST,
     &      .FALSE.,NBASISMAX,.FALSE.)
         DO I=1,NLIST
            CALL WRITEDET(6,LIST(1,I),NEL,.FALSE.)
            WRITE(6,*) ICLIST(I)
         ENDDO
         RETURN
      END
      
      SUBROUTINE CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &   NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,
     &   DETINV,TSPECDET,SPECDET)
         include 'irat.inc'
         POINTER (IP_LSTE,LSTE),(IP_ICE,ICE),(IP_RIJLIST,RIJLIST)
         INTEGER LSTE(NEL,NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER ICE(NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         REAL*8 RIJLIST(*)
         INTEGER NMRKS(NEL,NDET),NPATHS
         INTEGER III,NWHTAY,I,IMAX,ILMAX
         REAL*8 WLRI,WLSI,ECORE,DBETA,WLRI1,WLRI2,WLSI1,WLSI2,WI
         REAL*8 TOT, NORM,WLRI0,WLSI0,WINORM
         LOGICAL TNPDERIV
         INTEGER DETINV
         INTEGER ISTART,IEND
         LOGICAL TSPECDET
         INTEGER SPECDET(NEL)
         TOT=0.D0
         NORM=0.D0
         IMAX=I_HMAX
         IF(I_VMAX.GT.IMAX) IMAX=I_VMAX
         CALL TISET('CLCRHOPII2',ISUB)
         ILMAX=NDET
C.. we don't need lists for I_HMAX=8
         IF((I_HMAX.GE.-10.AND.I_HMAX.LE.-7)
     &      .OR.I_HMAX.LE.-12) ILMAX=1
C         ILMAX=(NBASIS-NEL)**2*NEL*NEL/4
         CALL MEMORY(IP_LSTE,(1+ILMAX)*NEL*IMAX/IRAT,"LSTE")
         CALL MEMORY(IP_ICE,(1+ILMAX)*IMAX/IRAT,"ICE")
        CALL MEMORY(IP_RIJLIST,(1+ILMAX)*IMAX*2,"RIJLIST")
         IF(I_VMAX.NE.0) THEN
            WRITE(6,*) "Using Vertex approximation.  I_VMAX=",I_VMAX
            IF(I_HMAX.EQ.0)
     &  WRITE(6,*) "I_HMAX=0.  Summing all I_HMAX up to P using contour"
            IF(I_HMAX.GT.0) WRITE(6,*) "I_HMAX=",I_HMAX
         ELSEIF(I_HMAX.NE.0) THEN
            WRITE(6,*) "Using hop-restricted paths. I_HMAX:",I_HMAX
         ELSE
            WRITE(6,*) "I_HMAX=I_VMAX=0. Using rho diagonalisation."
         ENDIF
         IF(I_HMAX.EQ.-10) THEN
            OPEN(11,FILE="MCSUMMARY",STATUS="UNKNOWN")
            WRITE(11,*) "Calculating ",NPATHS," W_Is..."
            CLOSE(11)
         ELSE
            OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
            WRITE(11,*) "Calculating ",NPATHS," W_Is..."
            CLOSE(11)
         ENDIF
         OPEN(15,FILE='RHOPII',STATUS='UNKNOWN')
         IF(DETINV.NE.0) THEN
            ISTART=ABS(DETINV)
            IEND=ABS(DETINV)
         ELSEIF(TSPECDET) THEN
            ISTART=0
            IEND=0
            WRITE(6,*) "Using specified det:"
            CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
         ELSE
            ISTART=1
            IEND=NPATHS
         ENDIF
         DO III=ISTART,IEND
            IF(III.NE.0) THEN
              IF(NPATHS.EQ.1)  CALL WRITEDET(6,NMRKS(1,III),NEL,.TRUE.) 
               CALL MCPATHSR3(NMRKS(1,III),BETA,I_P,I_HMAX,I_VMAX,NEL,
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,
     &         WLRI,WLSI,DBETA,DLWDB2)
            ELSE
               CALL MCPATHSR3(SPECDET,BETA,I_P,I_HMAX,I_VMAX,NEL,
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,
     &         WLRI,WLSI,DBETA,DLWDB2)
            ENDIF
            WRITE(15,"(I12,$)") III
            IF(TSPECDET) THEN
               CALL WRITEDET(15,SPECDET,NEL,.FALSE.)
            ELSE
               CALL WRITEDET(15,NMRKS(1,III),NEL,.FALSE.)
            ENDIF
            WRITE(15,"(A,3G25.16,$)") " ",EXP(WLSI+I_P*WLRI),
     &         WLRI*I_P,WLSI
            IF(III.EQ.1) THEN
               WLRI0=WLRI
               WLSI0=WLSI
            ENDIF  
            IF(TNPDERIV) THEN
C.. if we're calculating the derivatives too
               IF(III.NE.0) THEN
               CALL MCPATHSR3(NMRKS(1,III),BETA+DBETA,I_P,I_HMAX,
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,
     &            ILOGGING,ECORE,ILMAX,WLRI1,WLSI1,DBETA,DLWDB3)
               CALL MCPATHSR3(NMRKS(1,III),BETA-DBETA,I_P,I_HMAX,
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,
     &            ILOGGING,ECORE,ILMAX,WLRI2,WLSI2,DBETA,DLWDB4)
               ELSE
               CALL MCPATHSR3(SPECDET,BETA+DBETA,I_P,I_HMAX,
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,
     &            ILOGGING,ECORE,ILMAX,WLRI1,WLSI1,DBETA,DLWDB3)
               CALL MCPATHSR3(SPECDET,BETA-DBETA,I_P,I_HMAX,
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,
     &            ILOGGING,ECORE,ILMAX,WLRI2,WLSI2,DBETA,DLWDB4)
               ENDIF
               DLWDB=-(I_P*(WLRI1-WLRI2)+(WLSI1-WLSI2))/(2*DBETA)
            ELSE
               DLWDB=DLWDB2
            ENDIF
C.. we calculate the energy with weightings normalized to the weight of
C.. the Fermi determinant, otherwise the numbers blow up
            WINORM=EXP(I_P*(WLRI-WLRI0)+(WLSI-WLSI0))
            NORM=NORM+WINORM
            TOT=TOT+WINORM*DLWDB
            WRITE(15,*) DLWDB
            IF(DETINV.EQ.III.AND.III.NE.0) THEN
               CALL FLUSH(15)
               WRITE(6,*) "Investigating det ",DETINV
               CALL FLUSH(6)
               CALL WIRD_SUBSET(NMRKS(1,DETINV),BETA,I_P,NEL,
     &       NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &       RHOEPS,ILOGGING,TSYM,ECORE)
            ENDIF
          ENDDO
         CLOSE(15)
         WRITE(6,*) "Summed approx E(Beta)=",TOT/NORM
         CALL FREEM(IP_RIJLIST)
         CALL FREEM(IP_LSTE)
         CALL FREEM(IP_ICE)
         CALL TIHALT('CLCRHOPII2',ISUB)
         RETURN
      END    

      SUBROUTINE CALCRHOSD(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &            NTAY,RHOEPS,NWHTAY,ECORE)
         POINTER (IP_LSTE,LSTE),(IP_ICE,ICE)
         INTEGER LSTE(NEL,NBASIS*NBASIS*NEL*NEL,0:I_HMAX-1)
         INTEGER ICE(NBASIS*NBASIS*NEL*NEL,0:I_HMAX-1)
         INTEGER NMRKS(NEL,NDET)
         INTEGER III,NWHTAY,III2,K
         REAL*8 WLRI,WLSI,RH
         CALL TISET('CALCRHOSD ',ISUB)
         OPEN(15,FILE='RHOS',STATUS='UNKNOWN')
         DO III=1,NDET
           DO III2=1,NDET
CNDET
            CALL CALCRHO2(NMRKS(1:NEL,III),NMRKS(1:NEL,III2),
     &                  BETA,I_P,
     &                  NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,
     &                  ALAT,UMAT,RH,NTAY,-1,ECORE)
            IF(ABS(RH).GT.RHOEPS) THEN
              WRITE(15,"(2I12,A,$)") III,III2,"("
              DO K=1,NEL
                 WRITE(15,"(I3,A,$)") NMRKS(K,III),","
              ENDDO
              WRITE(15,"(A,$)") "),("
              DO K=1,NEL
                 WRITE(15,"(I3,A,$)") NMRKS(K,III2),","
              ENDDO
              WRITE(15,"(A,$)") ")"
              WRITE(15,*) RH
            ENDIF
           ENDDO
         ENDDO
         CLOSE(15)
         CALL TIHALT('CALCRHOSD ',ISUB)
         RETURN
      END    


      SUBROUTINE TC(I_H,I_P,NWHTAY)
         REAL*8 RHOII(0:I_H-1)
         INTEGER IPATH(I_H)
         REAL*8 RAN2
         RHOII(0)=1.D0
         RHOII(1)=0.9
         ISEED=-7
         DO I=2,I_H-1
            RHOII(I)=1-0.1*RAN2(ISEED)
         ENDDO
         DO IH=2,I_P
         WRITE(6,*) "H:",IH," P:",I_P
         WRITE(6,*) "TC:",CALCPATH6_2(IPATH,RHOII,1.D0,I_H,I_P,
     &      IH,1,0)
C         WRITE(6,*) "TC2:",CALCPATH5(IPATH,RHOII,1.D0,IH,I_P,0)
         ENDDO
         RETURN
      END

      SUBROUTINE HAMHIST(HMIN,HMAX,NHS,NHISTBOXES)
         INTEGER IHIST(NHISTBOXES)
         WRITE(6,'(A,F19.5,F19.5)') "MIN,MAX H", HMIN,HMAX
         WRITE(6,'(A,I12)') "TOTAL # H ELEMENTS", NHS 
         OPEN(1,FILE="HAMIL",STATUS="OLD")

C..Read in input parameters
         BOXSPACING=(HMAX-HMIN)/NHISTBOXES
         DO I=1,NHS
            READ(1,*) II,IJ,H
            IBOX=INT((H-HMIN)/BOXSPACING)+1
C.. This should only happen for HMAX
            IF(IBOX.GT.NHISTBOXES) IBOX=NHISTBOXES
            IHIST(IBOX)=IHIST(IBOX)+1
         ENDDO
         CLOSE(1)
         OPEN(1,FILE="Histogram.dat")
         DO I=1,NHISTBOXES
            WRITE(1,"(F19.7,I12)") HMIN+BOXSPACING*(I-1),IHIST(I)
         ENDDO
         CLOSE(1)
         RETURN
      END

      REAL*8 FUNCTION GETRHOEPS(RHOEPSILON,BETA,NEL,NBASISMAX,G1,NHG,
     &      BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,I_P,ECORE)
         USE HElement
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         INTEGER NEL,NI(NEL),I,NBASISMAX(5,2),I_P
         INTEGER BRR(*),NMSH,NMAX,NHG
         COMPLEX*16 FCK(*)
         REAL*8 RHOEPSILON,BETA,ECORE,ALAT(*)
         TYPE(HElement) BP,UMat(*),TMat(*)
         TYPE(BasisFN) G1(*)
         DO I=1,NEL
            NI(I)=BRR(I)
         ENDDO
         CALL SORTI(NEL,NI)
         BP=HElement(-BETA/I_P)
        GETRHOEPS=DSQRT(SQ(HElement(RHOEPSILON)*EXP(BP*GETHELEMENT2(NI,
     &      NI,NEL,NBASISMAX,G1,NHG,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,
     &         0,ECORE))))
         RETURN
      END

      SUBROUTINE DOZSUM()
         REAL*8 XIS(4)
         REAL*8 SUM,SUMI
         INTEGER I
         XIS(1)=1.0
         XIS(2)=10.9
         XIS(3)=0.9
         XIS(4)=0.9
         DO I=2,2000,2
            CALL ZSUM(XIS,4,10,I,1.1D0*10.9,SUM,SUMI)
            WRITE(16,*) I,SUM,SUMI
         ENDDO
         STOP
      END

      REAL*8 FUNCTION CALCT(NI,NEL,G1,TMAT,NBASIS)
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),G1(*),NBASIS,I
         TYPE(HElement) TMAT(NBASIS,NBASIS)
         LOGICAL ISCSF
         CALCT=0.D0
         IF(ISCSF(NI,NEL)) RETURN
         DO I=1,NEL
            CALCT=CALCT+DREAL(TMAT(NI(I),NI(I)))
         ENDDO
         RETURN
      END

      REAL*8 FUNCTION CALCT2(NI,NEL,G1,ALAT,NBASIS,CST)
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,I,J
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         REAL*8 ALAT(4),CST,TMAT
         LOGICAL ISCSF
         CALCT2=0.D0
         IF(ISCSF(NI,NEL)) RETURN
         DO J=1,NEL
            I=NI(J)
           TMAT=((ALAT(1)**2)*((G1(I)%K(1)**2)/(ALAT(1)**2)+
     &         (G1(I)%K(2)**2)/(ALAT(2)**2)+
     &         (G1(I)%K(3)**2)/(ALAT(3)**2)))
           TMAT=TMAT*CST
           CALCT2=CALCT2+TMAT
         ENDDO
         RETURN
      END
       REAL*8 FUNCTION DOEXMC(NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ECORE,
     &     IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,NBASIS,BRR,IEQSTEPS)
         INTEGER NDET,NEVAL,I_P,ILOGGING
         REAL*8 CK(NEVAL),W(NEVAL),BETA,ECORE

         REAL*8 DLWDBS(NDET),WLRIS(NDET),WLSIS(NDET),EN
         INTEGER I
         LOGICAL TWARN
         
         DO I=1,NDET
            CALL CALCRHOPII(I,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ECORE,
     &               WLRIS(I),WLSIS(I),TWARN)
            DLWDBS(I)=CALCDLWDB(I,NDET,NEVAL,CK,W,BETA,ECORE)
         ENDDO
         EN=DMONTECARLOEXWI(NDET,WLRIS,WLSIS,DLWDBS,I_P,
     &         IMCSTEPS,
     &         G1,NMRKS,NEL,NBASISMAX,NBASIS,BRR,IEQSTEPS,
     &         ILOGGING)
         WRITE(6,*) "EXACT MC RESULT=",EN
         DOEXMC=EN
         RETURN
      END

      SUBROUTINE DOBEGRAPH(NDET,NEVAL,CK,W,I_P,ILOGGING,
     &   G1,NMRKS,NEL,NBASISMAX,NHG,BRR)
         INTEGER NDET,NEVAL,I_P,ILOGGING
         REAL*8 CK(NEVAL),W(NEVAL),BETA,ETRIAL

         REAL*8 WLRI,WLSI
         LOGICAL TWARN
         DO BETA=1,51,5
            WRITE(13,"(E25.16,$)") BETA
            DO ETRIAL=15.5,16.5,0.1
               
            CALL CALCRHOPII(1,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ETRIAL,
     &               WLRI,WLSI,TWARN)
           WRITE(13,"(E25.16,$)") WLRI*I_P+WLSI
            ENDDO
            WRITE(13,*)
         ENDDO
         RETURN
      END   
      
      SUBROUTINE DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &            NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
         POINTER (IP_LSTE,LSTE),(IP_ICE,ICE),(IP_RIJLIST,RIJLIST)
         INTEGER LSTE(NEL,NBASIS*NBASIS*NEL*NEL,0:I_HMAX-1)
         INTEGER ICE(NBASIS*NBASIS*NEL*NEL,0:I_HMAX-1)
         REAL*8 RIJLIST(*)
         INTEGER NMRKS(NEL,NDET),NPATHS
         INTEGER III,NWHTAY,I,IMAX,NI(NEL),BRR(NBASIS)
         REAL*8 WLRI,WLSI,ETRIAL
         IMAX=I_HMAX
         IF(I_VMAX.GT.IMAX) IMAX=I_VMAX
         CALL TISET('BEGRAPHAP ',ISUB)
         CALL MEMORY(IP_LSTE,NEL*NBASIS*NBASIS*NEL*NEL*IMAX,"LSTE")
         CALL MEMORY(IP_ICE,NBASIS*NBASIS*NEL*NEL*IMAX,"ICE")
         CALL MEMORY(IP_RIJLIST,NBASIS*NBASIS*NEL*NEL*IMAX*2,"RIJLIST")
         DO I=1,NEL
            NI(I)=BRR(I)
         ENDDO
         CALL SORTI(NEL,NI)
         DO BETA=1,51,2
            WRITE(13,"(E25.16,$)") BETA
            DO ETRIAL=15.8,16.3,0.01
                          
               CALL MCPATHSR3(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ETRIAL,WLRI,WLSI)
               WRITE(13,"(E25.16,$)") WLRI*I_P+WLSI
            ENDDO
            WRITE(13,*)
         ENDDO
         CALL FREEM(IP_RIJLIST)
         CALL FREEM(IP_LSTE)
         CALL FREEM(IP_ICE)
         CALL TIHALT('BEGRAPHAP ',ISUB)
         RETURN
      END    
! Calculate the diagonal matrix elements for the Uniform electron gas
!  CST is PI*PI/2L*L for the non-periodic case, and
!  CST is 4*PI*PI/2L*L for the periodic case, and

!  For the periodic case we must also add in a periodic images correction
!  which is 1/2 (<ii|ii> - <ii|ii>cell) for each orbital
!  We calculate <ii|ii>cell with a potenial v(r)=1/r (r<Rc) and 0 (r>=Rc)
!  Rc=ALAT(4).

      SUBROUTINE CALCTMATUEG(TMAT,NBASIS,ALAT,G1,CST,TPERIODIC,OMEGA)
         IMPLICIT NONE
         INTEGER NBASIS
         INCLUDE 'basis.inc'
         TYPE(BASISFN) G1(NBASIS)
         REAL*8 TMAT(NBASIS,NBASIS),ALAT(4),HFBASIS(NBASIS,NBASIS),CST
         INTEGER I,J
         REAL*8 SUM,S1,OMEGA
         LOGICAL TPERIODIC
         REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795029D0
         IF(TPERIODIC) WRITE(6,*) "Periodic UEG"
         OPEN(10,FILE='TMAT',STATUS='UNKNOWN')
          DO I=1,NBASIS
           TMAT(I,I)=((ALAT(1)**2)*((G1(I)%K(1)**2)/(ALAT(1)**2)+
     &         (G1(I)%K(2)**2)/(ALAT(2)**2)+
     &         (G1(I)%K(3)**2)/(ALAT(3)**2)))
           TMAT(I,I)=TMAT(I,I)*CST
!..  The G=0 component is explicitly calculated for the cell interactions as 2 PI Rc**2 .
!     we *1/2 as we attribute only half the interaction to this cell.
           IF(TPERIODIC) TMAT(I,I)=TMAT(I,I)-PI*ALAT(4)**2/OMEGA
            WRITE(10,*) I,I,TMAT(I,I)
         ENDDO
         CLOSE(10)
         RETURN
      END

      SUBROUTINE WRITEDET(NUNIT,NI,NEL,LTERM)
         IMPLICIT NONE
         INTEGER NUNIT,NEL,NI(NEL),I
         LOGICAL LTERM
         INTEGER IEL
         CHARACTER*2 SUFF
         INCLUDE 'csf.inc'
         WRITE(NUNIT,"(A,$)") "("
         DO I=1,NEL
            IEL=NI(I)
            IF(IEL.GE.CSF_NBSTART) THEN
               WRITE(NUNIT,"(I3,$)"),(IEL-CSF_NBSTART)/4+1
               IEL=IAND(IEL-CSF_NBSTART,3)
               IF(IEL.EQ.0) THEN
                  WRITE(NUNIT,"(A,$)") "-B,"
               ELSEIF(IEL.EQ.1) THEN
                  WRITE(NUNIT,"(A,$)") "-A,"
               ELSEIF(IEL.EQ.2) THEN
                  WRITE(NUNIT,"(A,$)") "+B,"
               ELSE
                  WRITE(NUNIT,"(A,$)") "+A,"
               ENDIF
            ELSE
               WRITE(NUNIT,"(I5,A,$)") IEL,","
            ENDIF
         ENDDO
         WRITE(NUNIT,"(A,$)") ")"
         IF(LTERM) WRITE(NUNIT,*)
         RETURN
      END
      LOGICAL FUNCTION KALLOWED(G,NBASISMAX)
         IMPLICIT NONE
         INTEGER G(5),NBASISMAX(5,5),NMAXX,I,J,AX,AY
         INTEGER KX,KY
         REAL*8 MX,MY,XX,YY
         LOGICAL TALLOW
         TALLOW=.TRUE.
         IF(NBASISMAX(3,3).EQ.1) THEN
C.. spatial symmetries
            IF(G(1).NE.0) TALLOW=.FALSE.
         ELSEIF(NBASISMAX(3,3).EQ.0) THEN
C.. Hubbard
            IF(NBASISMAX(1,3).EQ.1) THEN
C.. Tilted hubbard
               NMAXX=NBASISMAX(1,5)
C            NMAXY=
               AX=NBASISMAX(1,4)
               AY=NBASISMAX(2,4)
C.. (XX,YY) is the position of the bottom right corner of the unit cell
               XX=((AX+AY)/2.D0)*NMAXX
               YY=((AY-AX)/2.D0)*NMAXX
               MX=XX*AX+YY*AY
               MY=XX*AY-YY*AX
               I=G(1)
               J=G(2)
               KX=I*AX+J*AY
               KY=I*AY-J*AX
               IF(KX.GT.MX) TALLOW=.FALSE.
               IF(KY.GT.MY) TALLOW=.FALSE.
               IF(KX.LE.-MX) TALLOW=.FALSE.
               IF(KY.LE.-MY) TALLOW=.FALSE.
            ELSEIF(NBASISMAX(1,3).GE.4.OR.NBASISMAX(1,3).EQ.2) THEN
C.. Real space Hubbard
               IF(G(1).EQ.0.AND.G(2).EQ.0.AND.G(3).EQ.0) THEN
                  TALLOW=.TRUE.
               ELSE
                  TALLOW=.FALSE.
               ENDIF
C            ELSEIF(NBASISMAX(1,3).EQ.2) THEN
C.. mom space non-pbc non-tilt hub - parity sym
C               IF(  (G(1).EQ.0.OR.G(1).EQ.1)
C     &         .AND.(G(2).EQ.0.OR.G(2).EQ.1)
C     &         .AND.(G(3).EQ.0.OR.G(3).EQ.1)) THEN
C                  TALLOW=.TRUE.
C               ELSE
C                  TALLOW=.FALSE.
C               ENDIF
C            ELSEIF(NBASISMAX(1,3).EQ.2) THEN
C.. non-pbc hubbard
C               TALLOW=.TRUE.
C               IF(G(1).GT.NBASISMAX(1,2).OR.G(1).LT.NBASISMAX(1,1))
C     &            TALLOW=.FALSE.
C               IF(G(2).GT.NBASISMAX(2,2).OR.G(2).LT.NBASISMAX(2,1))
C     &            TALLOW=.FALSE.
C               IF(G(3).GT.NBASISMAX(3,2).OR.G(3).LT.NBASISMAX(3,1))
C     &            TALLOW=.FALSE.
            ELSE
C.. Normal Hubbard
               TALLOW=.TRUE.
               IF(G(1).GT.NBASISMAX(1,2).OR.G(1).LT.NBASISMAX(1,1))
     &            TALLOW=.FALSE.
               IF(G(2).GT.NBASISMAX(2,2).OR.G(2).LT.NBASISMAX(2,1))
     &            TALLOW=.FALSE.
               IF(G(3).GT.NBASISMAX(3,2).OR.G(3).LT.NBASISMAX(3,1))
     &            TALLOW=.FALSE.
            ENDIF
         ENDIF         
         KALLOWED=TALLOW
         RETURN
      END

      SUBROUTINE GENRANDOMDET(NEL,NBASIS,MCDET)
         IMPLICIT NONE
         INTEGER NEL,NBASIS,MCDET(NEL)
         INTEGER I,J,EL,SEED
         LOGICAL BR
         REAL*8 RAN2
         SEED=-7
         DO I=1,NEL
            BR=.TRUE.
            DO WHILE (BR)
               BR=.FALSE.
               EL=INT(RAN2(SEED)*NBASIS+1)
               DO J=1,I-1
                  IF(MCDET(J).EQ.EL) BR=.TRUE.
               ENDDO
            ENDDO
            MCDET(I)=EL
         ENDDO
         CALL SORTI(NEL,MCDET)
         RETURN
      END

      SUBROUTINE CALCRHOPII3(BETA,I_P,I_HMAX,I_VMAX,NEL,
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,NTAY,
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,
     &            DETINV,TSPN,LMS,TPARITY,SymRestrict,TSPECDET,SPECDET)
         USE HElement
         IMPLICIT NONE
         include 'irat.inc'
         POINTER (IP_LSTE,LSTE),(IP_ICE,ICE),(IP_RIJLIST,RIJLIST)
         INTEGER I_HMAX,NEL,NBASIS,I_VMAX
         INTEGER LSTE(NEL,NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER ICE(NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         TYPE(HElement)  RIJLIST(*),UMAT(*),TMat(*)
         REAL*8 BETA,FCK(*),ALAT(*),RHOEPS
         INTEGER NPATHS,NI(NEL),I_P,NBASISMAX(5,5)
         INCLUDE 'gndwork.inc'
         INTEGER Work(GNDWorkSize+2*NEL)
         TYPE(BASISFN) G1(NBASIS)
         INTEGER BRR(NBASIS),NMSH,NMAX(*),NTAY,ILOGGING
         INTEGER III,NWHTAY,I,IMAX,ILMAX,LMS,ISUB
         TYPE(BasisFN) ISYM,SymRestrict
         LOGICAL TSPN,TPARITY,TSYM
         REAL*8 DBETA,ECORE
         TYPE(HDElement) WLRI,WLSI,WLRI1,WLRI2,WLSI1,WLSI2,WI,DLWDB
         TYPE(HDElement) TOT,WLRI0,WLSI0,WINORM,HElP,NORM
         LOGICAL TNPDERIV,TDONE,TFIRST
         INTEGER DETINV
         INTEGER ISTART,IEND,IDEG
         LOGICAL TSPECDET
         INTEGER SPECDET(NEL)
         TYPE(HDElement) DLWDB2,DLWDB3,DLWDB4,TOT2
         HElP=HDElement(I_P)
         TSYM=.TRUE.
         TOT=0.D0
         TOT2=0.D0
         NORM=0.D0
         IMAX=I_HMAX
         IF(I_VMAX.GT.IMAX) IMAX=I_VMAX
         CALL TISET('CLCRHOPII3',ISUB)
         WRITE(6,*) "Entering CALCRHOPII3..."
C         ILMAX=NDET
C.. We don't need to store lists for I_HMAX=-8
         ILMAX=(NBASIS-NEL)**2*NEL*NEL/4
         IF((I_HMAX.GE.-10.AND.I_HMAX.LE.-7)
     &     .OR.I_HMAX.LE.-12) ILMAX=1
         CALL MEMORY(IP_LSTE,(1+ILMAX)*NEL*IMAX/IRAT,"LSTE")
         CALL MEMORY(IP_ICE,(1+ILMAX)*IMAX/IRAT,"ICE")
        CALL MEMORY(IP_RIJLIST,(1+ILMAX)*IMAX*2,"RIJLIST")
C:         CALL PRINT_MEMORY()
         IF(I_VMAX.NE.0) THEN
            WRITE(6,*) "Using Vertex approximation.  I_VMAX=",I_VMAX
            IF(I_HMAX.EQ.0)
     &  WRITE(6,*) "I_HMAX=0.  Summing all I_HMAX up to P using contour"
            IF(I_HMAX.GT.0) WRITE(6,*) "I_HMAX=",I_HMAX
         ELSEIF(I_HMAX.NE.0) THEN
            WRITE(6,*) "Using hop-restricted paths. I_HMAX:",I_HMAX
         ELSE
            WRITE(6,*) "I_HMAX=I_VMAX=0. Using rho diagonalisation."
         ENDIF
         IF(I_HMAX.EQ.-10) THEN
            OPEN(11,FILE="MCSUMMARY",STATUS="UNKNOWN")
            WRITE(11,*) "Calculating ",NPATHS," W_Is..."
            CLOSE(11)
         ELSE
            OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
            WRITE(11,*) "Calculating ",NPATHS," W_Is..."
            CLOSE(11)
         ENDIF
         OPEN(15,FILE='RHOPII',STATUS='UNKNOWN')
         IF(DETINV.NE.0) THEN
            ISTART=ABS(DETINV)
            IEND=ABS(DETINV)
         ELSEIF(TSPECDET) THEN
            WRITE(6,*) "Calculating vertex series for specific det:"
            CALL WRITEDET(6,SPECDET,NEL,.TRUE.) 
            ISTART=-1
            IEND=1
         ELSE
            ISTART=1
            IEND=NPATHS
         ENDIF

         III=0
         TDONE=.FALSE.
         TFIRST=.TRUE.
         IF(.NOT.TSPECDET)
     &       CALL GENNEXTDET(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &         TPARITY,SymRestrict,ISYM,NI,.TRUE.,TDONE,WORK)
         DO WHILE(III.NE.IEND.AND..NOT.TDONE)
          III=III+1
          IF(TSPECDET) THEN
             TDONE=.FALSE.
             CALL ICOPY(NEL,SPECDET,1,NI,1)
             IDEG=1
          ELSE
             CALL GENNEXTDET(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &           TPARITY,SymRestrict,ISYM,NI,.FALSE.,TDONE,WORK)
             CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
          ENDIF
          IF(III.GE.ISTART.AND..NOT.TDONE) THEN
            IF(NPATHS.EQ.1.AND..NOT.TSPECDET)
     &         CALL WRITEDET(6,NI,NEL,.TRUE.) 
            CALL MCPATHSR3(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,
     &         WLRI,WLSI,DBETA,DLWDB2)
            WRITE(15,"(I12,$)") III
            CALL WRITEDET(15,NI,NEL,.FALSE.)
            WRITE(15,"(3G25.16,$)") EXP(WLSI+HElP*WLRI),
     &         WLRI*HElP,WLSI
            IF(TFIRST) THEN
               TFIRST=.FALSE.
               WLRI0=WLRI
               WLSI0=WLSI
            ENDIF  
            IF(TNPDERIV) THEN
C.. if we're calculating the derivatives too
               CALL MCPATHSR3(NI,BETA+DBETA,I_P,I_HMAX,
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,
     &            ILOGGING,ECORE,ILMAX,WLRI1,WLSI1,DBETA,DLWDB3)
               CALL MCPATHSR3(NI,BETA-DBETA,I_P,I_HMAX,
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,
     &            ILOGGING,ECORE,ILMAX,WLRI2,WLSI2,DBETA,DLWDB4)
               DLWDB=-(HElP*(WLRI1-WLRI2)+(WLSI1-WLSI2))
     &                  /HDElement(2*DBETA)
            ELSE
               DLWDB=DLWDB2
            ENDIF
C.. we calculate the energy with weightings normalized to the weight of
C.. the Fermi determinant, otherwise the numbers blow up
            WINORM=EXP(HElP*(WLRI-WLRI0)+(WLSI-WLSI0))
            NORM=NORM+HDElement(IDEG)*WINORM
            TOT=TOT+HDElement(IDEG)*WINORM*DLWDB
            WRITE(15,"(G25.16,I5)") DLWDB,IDEG
            IF(DETINV.EQ.III) THEN
               CALL FLUSH(15)
               WRITE(6,*) "Investigating det ",DETINV
               CALL FLUSH(6)
               CALL WIRD_SUBSET(NI,BETA,I_P,NEL,
     &       NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,NTAY,
     &       RHOEPS,ILOGGING,TSYM,ECORE)
            ENDIF
           ENDIF
          ENDDO
         CLOSE(15)
         IF(TFIRST) THEN
            WRITE(6,*) "*** NO determinants found to calculate***"
         ENDIF
         WRITE(6,*) "Summed approx E(Beta)=",TOT/NORM
         CALL FREEM(IP_RIJLIST)
         CALL FREEM(IP_LSTE)
         CALL FREEM(IP_ICE)
         CALL TIHALT('CLCRHOPII3',ISUB)
         RETURN
      END    

      SUBROUTINE ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
         INTEGER NBASIS,BRR(NBASIS),ORBORDER(8,2),NBASISMAX(5,3)
         INTEGER BRR2(NBASIS)
         INCLUDE 'basis.inc'
         TYPE(BASISFN) G1(NBASIS)
         REAL*8 ARR(NBASIS,2),ARR2(NBASIS,2)
         INTEGER IDONE,I,J,IBFN,ITOT,ITYPE,K
         REAL*8 OEN
         IDONE=0
         ITOT=0
C.. copy the default ordered energies.
         CALL DCOPY(NBASIS,ARR(1,1),1,ARR(1,2),1)
         CALL DCOPY(NBASIS,ARR(1,1),1,ARR2(1,2),1)
         WRITE(6,"(A,8I4)") "Ordering Basis (Closed): ",
     &      (ORBORDER(I,1),I=1,8)
         WRITE(6,"(A,8I4)") "Ordering Basis (Open  ): ",
     &      (ORBORDER(I,2),I=1,8)
         IF(NBASISMAX(3,3).EQ.1) THEN
C.. we use the symmetries of the spatial orbitals
            DO ITYPE=1,2
             IBFN=1
             DO I=1,8
               DO J=1,ORBORDER(I,ITYPE)
                  DO WHILE(IBFN.LE.NBASIS
     &              .AND.(G1(IBFN)%SYM%s.LT.2**(I-1).OR.BRR(IBFN).EQ.0))
                     IBFN=IBFN+1
                  ENDDO
                  IF(IBFN.GT.NBASIS) THEN
                     STOP
     &                "Cannot find enough basis fns of correct symmetry"
                  ENDIF
                  IDONE=IDONE+1 
                  BRR2(IDONE)=IBFN
                  BRR(IBFN)=0
                  ARR2(IDONE,1)=ARR(IBFN,1)
                  IBFN=IBFN+1
               ENDDO
             ENDDO
             CALL SORT2(IDONE-ITOT,ARR2(ITOT+1,1),BRR2(ITOT+1))
             ITOT=IDONE
            ENDDO
            DO I=1,NBASIS
               IF(BRR(I).NE.0) THEN
                  ITOT=ITOT+1
                  BRR2(ITOT)=BRR(I)
                  ARR2(ITOT,1)=ARR(I,1)
               ENDIF
            ENDDO
            CALL ICOPY(NBASIS,BRR2,1,BRR,1)
            CALL DCOPY(NBASIS,ARR2,1,ARR,1) 
         ENDIF
         CALL SORT2(NBASIS-IDONE,ARR(IDONE+1,1),BRR(IDONE+1))
C.. We need to now go through each set of degenerate orbitals, and make
C.. the correct ones are paired together in BRR otherwise bad things
C.. happen in FREEZEBASIS
C.. We do this by ensuring that within a degenerate set, the BRR are in
C.. ascending order
         OEN=ARR(1,1)
         J=1
C         IF(NBASISMAX(3,3).EQ.1) G1(3,BRR(1))=J
         ITOT=1
         DO I=2,NBASIS
            IF(ABS(ARR(I,1)-OEN).GT.1.D-4) THEN
C.. We don't have degenerate orbitals
C.. First deal with the last set of degenerate orbitals
C.. We sort them into order of BRR
               CALL SORT2_(ITOT,BRR(I-ITOT),ARR(I-ITOT,1))
C.. now setup the new degenerate set.
               J=J+1
               ITOT=1
            ELSE
               ITOT=ITOT+1
            ENDIF
            OEN=ARR(I,1)
            IF(NBASISMAX(3,3).EQ.1) THEN
C.. If we've got a generic spatial sym or hf we mark degeneracies
C               G(3,BRR(I))=J
            ENDIF
         ENDDO
         CALL SORT2_(ITOT,BRR(I-ITOT),ARR(I-ITOT,1))
      END

      SUBROUTINE FREEZEBASIS(NHG,NBASIS,UMAT,TMAT,UMAT2,TMAT2,ECORE,
     &         G1,NBASISMAX,ISS,ARR,BRR,NFROZEN,NTFROZEN,FDET,NEL,ALAT)
         USE HElement
         USE UMatCache
         IMPLICIT NONE
         INTEGER NHG,NBASIS,NBASISMAX(5,6),ISS
         INCLUDE 'basis.inc'
         TYPE(BASISFN) G1(NHG)
         TYPE(HElement) UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         TYPE(HElement) TMAT(NHG,NHG),UMAT2(*)
         REAL*8 ECORE
C.. was (NBASIS/ISS,NBASIS/ISS,NBASIS/ISS,NBASIS/ISS)
         TYPE(HElement) TMAT2(NBASIS,NBASIS)
         REAL*8 ARR(NHG,2),ARR2(NBASIS,2)
         INTEGER NFROZEN,BRR(NHG),BRR2(NBASIS),GG(NHG)
         TYPE(BASISFN) G2(NHG)
         INTEGER NTFROZEN
         INTEGER I,J,K,L,A,B,IP,JP,KP,LP,IDI,IDJ,IDK,IDL
         INTEGER IB,JB,KB,LB,AB,BB,IPB,JPB,KPB,LPB
         INTEGER IDIP,IDJP,IDKP,IDLP
         INTEGER IDA,IDB,SGN
         INTEGER FDET(NEL),NEL
         TYPE(Symmetry) KSYM
         REAL*8 ALAT(3)
C.. Just check to see if we're not in the middle of a degenerate set with the same sym
         IF(NFROZEN.GT.0) THEN
            IF(ABS(ARR(NFROZEN,1)-ARR(NFROZEN+1,1)).LT.1.D-6.AND.
     &        G1(BRR(NFROZEN))%SYM%s.EQ.G1(BRR(NFROZEN+1))%SYM%s) THEN
               STOP "Cannot freeze in the middle of a degenerate set"
            ENDIF
         ENDIF
         IF(NTFROZEN.GT.0) THEN
            IF(ABS(ARR(NHG-NTFROZEN,1)-ARR(NHG-NTFROZEN+1,1)).LT.1.D-6
     &         .AND.G1(BRR(NHG-NTFROZEN))%SYM%s
     &               .EQ.G1(BRR(NHG-NTFROZEN+1))%SYM%s)
     &       THEN
          STOP "Cannot freeze in the middle of a degenerate virtual set"
            ENDIF
         ENDIF

C.. At this point, we transform the UMAT and TMAT into a new UMAT and
C.. TMAT and Ecore with the frozen orbitals factored in
C..
C.. a,b are frozen spinorbitals
C.. E'core = Ecore+sum_a t_aa + sum_(a<b) (<ab|ab>-<ab|ba>)
C.. t'_ii = t_ii+ sum_a ( <ai|ai> - <ai|ia> ) 
C.. NHG contains the old number of orbitals
C.. NBASIS contains the new
C.. We first need to work out where each of the current orbitals will
C.. end up in the new set
         K=0
         DO I=1,NHG
            L=1
            DO J=1,NFROZEN
C.. if orb I is to be frozen, L will become 0
               IF(BRR(J).EQ.I) L=0
            ENDDO
            DO J=NBASIS+NFROZEN+1,NHG
C.. if orb I is to be frozen, L will become 0
               IF(BRR(J).EQ.I) L=0
            ENDDO
            IF(L.EQ.0) THEN
               GG(I)=0
            ELSE
C.. we've got an orb which is not to be frozen 
               K=K+L
C.. GG(I) is the new position in G of the (old) orb I
               GG(I)=K
C.. copy the eigenvalue table to the new location
               CALL ICOPY(BasisFNSize,G1(I),1,G2(K),1)
            ENDIF
         ENDDO
C.. Now construct the new BRR and ARR
         DO I=1,NBASIS
            BRR2(I)=GG(BRR(I+NFROZEN))
            ARR2(I,1)=ARR(I+NFROZEN,1)
         ENDDO
         DO I=1,NHG
            IF(GG(I).NE.0) ARR2(GG(I),2)=ARR(I,2)
         ENDDO
C.. First deal with Ecore
         DO A=1,NFROZEN
            AB=BRR(A)
            ECORE=ECORE+DREAL(TMAT(AB,AB))
            DO B=A+1,NFROZEN
               BB=BRR(B)
               CALL GTID(NBASISMAX,AB,IDA)
               CALL GTID(NBASISMAX,BB,IDB)
C.. No sign problems from permuations here as all perms even
               ECORE=ECORE+DREAL(
     &        GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDB,IDA,IDB))
CCCC.. If we have spin-independent integrals, or 
C.. if the spins are the same
               IF(G1(AB)%MS.EQ.G1(BB)%MS)
     &            ECORE=ECORE-DREAL(
     &      GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDB,IDB,IDA))
            ENDDO
         ENDDO
C.. now deal with the new TMAT
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
               TMAT2(IPB,JPB)=TMAT(IB,JB)
               DO A=1,NFROZEN
                  AB=BRR(A)
                  CALL GTID(NBASISMAX,AB,IDA)
C.. SGN takes into account permutationnness.
C                  SGN=1
C                  IF(IB.GT.AB) SGN=-SGN
C                  IF(JB.GT.AB) SGN=-SGN
                  IF(G1(IB)%MS.EQ.G1(JB)%MS)
     &               TMAT2(IPB,JPB)=
     &                  TMAT2(IPB,JPB)
     &        +GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDI,IDA,IDJ)
CCCC.. If we have spin-independent integrals, ISS.EQ.2.OR
C.. if the spins are the same
                  IF(G1(IB)%MS.EQ.G1(AB)%MS.AND.G1(AB)%MS.EQ.G1(JB)%MS)
     &               TMAT2(IPB,JPB)=TMAT2(IPB,JPB)
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDA,IDI,IDJ,IDA)
               ENDDO
!               WRITE(6,*) "T",TMAT(IB,JB),I,J,TMAT2(IPB,JPB)
!            IF(TMAT(IPB,JPB).AGT.1.D-9) WRITE(16,*) I,J,TMAT2(IPB,JPB)
            ENDDO
         ENDDO
         IF(NBASISMAX(1,3).GE.0.AND.ISS.NE.0) THEN
CC Only do the below if we've a stored UMAT
C.. Now copy the relevant matrix elements of UMAT across
C.. the primed (...P) are the new versions
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
                         UMAT2(UMatInd(IDIP,IDJP,IDKP,IDLP))
     &                     =UMAT(UMatInd(IDI,IDJ,IDK,IDL))
                        ENDIF
                     ENDIF
                   ENDDO
                  ENDIF
                ENDDO
               ENDIF
             ENDDO
            ENDIF
         ENDDO
         ELSEIF(Allocated(UMatCacheData)) THEN
!.. We've a UMAT2D and a UMATCACHE.  Go and Freeze them
C.. NHG contains the old number of orbitals
C.. NBASIS contains the new
C.. GG(I) is the new position in G of the (old) orb I
            CALL FreezeUMatCache(GG,NHG,NBASIS)
         ENDIF
         IF(ISS.EQ.0) CALL SetupUMatTransTable(GG,nHG,nBasis)
C.. Now setup the default symmetry to include the frozen electrons
C.. BRR(1:NFROZEN) is effectively the det of the frozens, so we get its sym
         CALL GETSYM(BRR,NFROZEN,G1,NBASISMAX,KSym)
         CALL SetupFreezeAllSym(KSYM)
C.. Freeze the sym labels
         CALL FREEZESYMLABELS(NHG,NBASIS,GG)
C.. Copy the new G1 and BRR and ARR over the old ones
         CALL ICOPY(BasisFNSize*NBASIS,G2,1,G1,1)
         CALL ICOPY(NBASIS,BRR2,1,BRR,1)
         CALL DCOPY(NBASIS*2,ARR2,1,ARR,1)
C.. Now reset the total number of orbitals
         NHG=NBASIS
C.. Deal with FDET
C.. GG(I) is the new position in G of the (old) orb I
         IF(FDET(1).NE.0) THEN
            J=0
            DO I=1,NEL
               FDET(I)=GG(FDET(I))
C.. any orbitals which no longer exist, we move outside the basis
               IF(FDET(I).EQ.0) THEN
                  FDET(I)=NHG+1
               ELSE
                  J=J+1
               ENDIF
            ENDDO
            CALL SORTI(NEL,FDET)
            IF(J.NE.NEL-NFROZEN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               CALL WRITEDET(6,NEL,FDET,.TRUE.)
               STOP "After Freezing, FDET has wrong number of electrons"
            ENDIF
         ENDIF
         RETURN
      END   
      SUBROUTINE GENSITESPINOCC(NBASIS,NSBASIS,NSPINS,NBASISMAX,G1,NEL,
     &      LMS,BRR,HFBASIS)
         IMPLICIT NONE
         INTEGER NBASIS,NBASISMAX(5,5),NSPINS,NSBASIS
         INTEGER NEL,LMS,BRR(NBASIS)
         INCLUDE 'basis.inc'
         TYPE(BASISFN) G1(NBASIS)
         REAL*8 HFBASIS(NBASIS,NBASIS)
         
         REAL*8 DMAT(NBASIS,NBASIS)
         COMPLEX*16 STOT,II
         INTEGER I,J,K,KP
         INTEGER ISPN,NELS(NSPINS),NDET(NEL),NSPN(NSPINS)
         REAL*8 TOT,RJ(3),PI,STT
         INTEGER AX,AY,NMAXX,NMAXY,NMAXZ,BR(NSPINS)
         REAL*8 AA2
         REAL*8 P,Q
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         PI=ACOS(-1.D0)
         AX=NBASISMAX(1,4)
         AY=NBASISMAX(2,4)
         AA2=AX*AX+AY*AY
         NMAXX=NBASISMAX(1,5)
         NMAXY=NBASISMAX(2,5)
         NMAXZ=NBASISMAX(3,5)

         II=(0.D0,1.D0)
         IF(NSPINS.EQ.2) THEN
            NELS(2)=(LMS+NEL)/2
            NELS(1)=NEL-NELS(2)
         ELSE
            NELS(1)=NEL
         ENDIF
         DO I=1,NSPINS
            NSPN(I)=0
            BR(I)=0
         ENDDO
         I=1
         DO WHILE (I.LE.NEL)
            DO ISPN=1,NSPINS
               IF(NSPN(ISPN).LT.NELS(ISPN)) THEN
                  BR(ISPN)=BR(ISPN)+1
                  DO WHILE(G1(BRR(BR(ISPN)))%MS.NE.ISPN*2-3)
                     BR(ISPN)=BR(ISPN)+1
                  ENDDO
                  NDET(I)=BRR(BR(ISPN))
                  NSPN(ISPN)=NSPN(ISPN)+1
                  I=I+1
               ENDIF
            ENDDO
         ENDDO 

         DO K=1,NBASIS
            DO KP=1,NBASIS
               TOT=0.D0
               DO I=1,NEL
                  TOT=TOT+HFBASIS(NDET(I),K)*HFBASIS(NDET(I),KP)
               ENDDO
               DMAT(K,KP)=TOT
            ENDDO
         ENDDO
         STT=0.D0
         WRITE(6,*) "Calculate site spin density..."
         DO J=1,NBASIS
            STOT=(0.D0,0.D0)         
            P=(G1(J)%K(1)*AX+G1(J)%K(2)*AY)/AA2
            Q=(-G1(J)%K(1)*AY+G1(J)%K(2)*AX)/AA2
            DO K=1,NBASIS
               DO KP=1,NBASIS
C.. Check spin
                  IF(MOD(K+J,2).EQ.0) 
     &             STOT=STOT+DMAT(K,KP)*EXP(2*PI*II*
     &               ((G1(KP)%K(1)-G1(K)%K(1))*P/(NMAXX+0.D0)+
     &               (G1(KP)%K(2)-G1(K)%K(2))*Q/(NMAXY+0.D0)))
               ENDDO
            ENDDO
            STOT=STOT/(NMAXX*NMAXY*NMAXZ*AA2)
            WRITE(6,*) (G1(J)%K(I),I=1,5),REAL(STOT)
            STT=STT+REAL(STOT)
         ENDDO
         WRITE(6,*) NEL,STT
         
      END
      SUBROUTINE WRITEBASIS(NUNIT,G1,NHG,ARR,BRR)
         IMPLICIT NONE
         INTEGER NUNIT,NHG,BRR(NHG),I
         INCLUDE 'sym.inc'
         TYPE(BASISFN) G1(NHG)
         REAL*8 ARR(NHG,2)
         DO I=1,NHG
           WRITE(NUNIT,'(6I7,Z10,2F19.9)') I,BRR(I),G1(BRR(I))%K(1),
     &         G1(BRR(I))%K(2),G1(BRR(I))%K(3),
     &         G1(BRR(I))%MS,G1(BRR(I))%SYM,
C     &         SymLabels(SymClasses((BRR(I)+1)/2)),
     &         ARR(I,1),ARR(BRR(I),2)
         ENDDO
         RETURN
      END
