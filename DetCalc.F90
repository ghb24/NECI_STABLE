#include "macros.h"
MODULE DetCalc
        Use HElem
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
    IMPLICIT NONE
     save

!From input
      REAL*8 B2L  ! From Calc
      INTEGER NEVAL,NBLK,NKRY,NCYCLE
      INTEGER DETINV,ICILEVEL
      INTEGER IOBS,JOBS,KOBS
      LOGICAL TRHOOFR,TCORR,TFODM


      INTEGER NDET



      LOGICAL TCALCHMAT,TENERGY,TREAD,TBLOCK
      POINTER (IP_NBLOCKSTARTS,NBLOCKSTARTS)
      TYPE(BasisFN), pointer :: BLOCKSYM(:)
      INTEGER tagBlockSym
      INTEGER NBLOCKSTARTS(*)
      INTEGER NBLOCKS
      Type(HElement), pointer :: HAMIL(:)
      INTEGER tagHamil
      INTEGER, pointer :: NMRKS(:,:)      !NEL-NFROZEN
      INTEGER tagNMRKS
      INTEGER iFDEt
      TYPE(HElement), pointer :: CK(:), CKN(:)
      INTEGER tagCK, tagCKN
      REAL*8, pointer :: W(:)
      INTEGER tagW
      INTEGER LenHamil
    
CONTAINS
    Subroutine DetCalcInit

        Use global_utilities
        Use Determinants, only:  FDet, specdet, tSpecDet
        use SystemData, only : tCSF,lms, lms2, nBasis, nBasisMax, nEl, SymRestrict
        use SystemData, only : Alat, arr, brr, boa, box, coa, ecore, g1,Beta
        use SystemData, only : tParity, tSpn
        Type(BasisFn) ISym

        integer i,ii
        integer ierr
        include 'irat.inc'
        integer nDetTot
        
        character(25), parameter :: this_routine='DetCalcInit'
        

      IF(.NOT.TCALCHMAT) THEN
         WRITE(6,*) "Not storing the H matrix."
         IF(TENERGY.AND..NOT.TBLOCK) THEN
            WRITE(6,*) "Cannot calculate energies without blocking the Hamiltonian."
            TENERGY=.FALSE.
         ENDIF
         IF(TENERGY.AND.NBLK.NE.0) THEN
!C.. We're doing a Lanczos without calculating the H mat
            WRITE(6,*) "Cannot perform Lanczos without Hamiltonian"
            TENERGY=.FALSE.
         ENDIF
      ENDIF
!C      IF(TCALCHMAT.OR.NPATHS.NE.0.OR.DETINV.GT.0.OR.TBLOCK) THEN
      IF(TENERGY.OR.TCALCHMAT) THEN
!C..Need to determine the determinants
         IF(ICILEVEL.NE.0) THEN
            WRITE(6,*) "Performing truncated CI at level ",ICILEVEL
            IF(TSPECDET) THEN
               WRITE(6,*) "Using SPECDET:"
               CALL WRITEDET(6,SPECDET,NEL,.TRUE.)!
               CALL ICOPY(NEL,SPECDET,1,FDET,1)
            ELSE
               WRITE(6,*) "Using Fermi DET:"
               CALL WRITEDET(6,FDET,NEL,.TRUE.)
            ENDIF 
!C.. if we're doing a truncated CI expansion
            CALL GENEXCIT(FDET,ICILEVEL,NBASIS,NEL,0,0,NDET,1,G1,.TRUE.,NBASISMAX,.TRUE.)
!C.. We need to add in the FDET
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
            CALL GNDTS_BLK(NEL,nBasis,BRR,NBASISMAX,NMRKS, .TRUE.,             &
     &            NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS2,TPARITY,        &
     &           SymRestrict,IFDET,.NOT.TREAD,NDETTOT,BLOCKSYM)
            WRITE(6,*) "NBLOCKS:",NBLOCKS
         ELSEIF(TCSF) THEN
            WRITE(6,*) "Determining CSFs."
!C determinants.
            WRITE(6,*) "Determining determinants and blocks."
            IF(TPARITY) THEN
               WRITE(6,*) "Using symmetry restriction:"
               CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
            ENDIF
            IF(TSPN) THEN
               WRITE(6,*) "Using spin restriction:",LMS
            ENDIF
            CALL GNCSFS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.TRUE.,G1,TSPN,LMS2,TPARITY,        &
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
            CALL GNDTS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.TRUE.,G1,TSPN,LMS,TPARITY,SymRestrict,II,IFDET)
            NBLOCKS=1
            NDET=II
         ENDIF
!C..
         IF(II.EQ.0) THEN
            WRITE(6,*) "No determinants found.  Cannot continue"
            STOP "No determinants found.  Cannot continue"
         ENDIF
!C.. NEL now only includes active electrons
         Allocate(NMrks(nEl,II),stat=ierr)
         LogAlloc(ierr,'NMRKS',NEL*II,4,tagNMRKS)
         CALL IAZZERO(NMRKS,(NEL)*II)
         CALL MEMORY(IP_NBLOCKSTARTS,(NBLOCKS+1)/IRAT+1,"NBLOCKSTARTS")
         CALL IAZZERO(NBLOCKSTARTS,NBLOCKS)
         Allocate(BlockSym(NBLOCKS+1),stat=ierr)
         LogAlloc(ierr, 'BLOCKSYM', NBLOCKS+1,BasisFNSizeB, tagBlockSym)

         CALL IAZZERO(BLOCKSYM,NBLOCKS*BasisFNSize)
!C..



         NDET=II   
         IF(ICILEVEL.NE.0) THEN
!C.. Use HAMIL to temporarily hold a list of excitation levels
            CALL ICOPY(NEL,FDET,1,NMRKS,1)
            Allocate(Hamil(II), stat=ierr)
            LogAlloc(ierr, 'HAMIL', II, HElementSizeB, tagHamil)
            CALL GENEXCIT(FDET,ICILEVEL,NBASIS,NEL,NMRKS(1,2),HAMIL,NDET,1,G1,.TRUE.,NBASISMAX,.FALSE.)
            Deallocate(Hamil)
            LogDealloc(tagHamil)
            NDET=NDET+1
            NBLOCKSTARTS(1)=1
            NBLOCKSTARTS(2)=II+1
            IFDET=1
         ELSEIF(TBLOCK) THEN 
            CALL GNDTS_BLK(NEL,nBasis,BRR,NBASISMAX,NMRKS, .FALSE.,NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS2,TPARITY, &
     &           SymRestrict,IFDET,.NOT.TREAD,NDETTOT,BLOCKSYM,TCSF)
         ELSEIF(TCSF) THEN
            CALL GNCSFS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.FALSE.,G1,TSPN,LMS2,TPARITY, &
     &         SymRestrict,NDET,IFDET,.FALSE.,0,0,.FALSE.,0)
               NBLOCKSTARTS(1)=1
               NBLOCKSTARTS(2)=II+1
         ELSE
            CALL GNDTS(NEL,nBasis,BRR,NBASISMAX,NMRKS,.FALSE.,G1,TSPN,LMS,TPARITY,SymRestrict,II,IFDET)
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
!C..
!C.. Now generate the fermi determiant
!C.. Work out the fermi det
         DO I=1,NEL
            FDET(I)=NMRKS(I,IFDET)
         ENDDO
         WRITE(6,*) "Fermi Determinant:",IFDET

         WRITE(6,*) ' NUMBER OF SYMMETRY UNIQUE DETS ' , NDET

!C         WRITE(6,*) ' TOTAL NUMBER OF DETS.' , NDETTOT
         IF(NEVAL.EQ.0) THEN
            WRITE(6,*) 'NEVAL=0.  Setting NEVAL=NDET'
            NEVAL=NDET
         ENDIF
         IF(NEVAL.GT.NDET) THEN
            WRITE(6,*) 'NEVAL>NDET.  Setting NEVAL=NDET'
            NEVAL=NDET
         ENDIF

         IF(ABS(DETINV).GT.NDET) THEN
            WRITE(6,*) 'DETINV=',DETINV,'>NDET=',NDET
            WRITE(6,*) 'Setting DETINV to 0'
            DETINV=0
         ENDIF
         CALL FLUSH(6)
    
!C ==----------------------------------------------------------------==
!C..Set up memory for c's, nrow and the label
         IF(TCALCHMAT) THEN
            WRITE(6,*) "CK Size",NDET*NEVAL*HElementSize
            Allocate(CkN(nDet*nEval), stat=ierr)
            LogAlloc(ierr,'CKN',nDet*nEval, HElementSizeB, tagCKN)
            CALL AZZERO(CKN,NDET*NEVAL*HElementSize)
!C..
            Allocate(Ck(nDet*nEval), stat=ierr)
            LogAlloc(ierr,'CK',nDet*nEval, HElementSizeB, tagCK)
            CALL AZZERO(CK,NDET*NEVAL*HElementSize)
!C..
            allocate(W(nEval), stat=ierr)
            LogAlloc(ierr, 'W', nEval,8,tagW)
            CALL AZZERO(W,NEVAL)
         ENDIF
!         IF(TREADRHO.AND..NOT.TREAD) THEN
!            WRITE(10,*) "TREADRHO specified, but not TREAD.  Setting TREAD=.T."
!            TREAD=.TRUE.
!         ENDIF
!C..
         IF(TREAD) THEN
           CALL READ_PSI(BOX,BOA,COA,NDET,NEVAL,NBASISMAX,NEL,CK,W)
         ENDIF
      ENDIF

      



!      TMC=TCALCHMAT.AND.(.NOT.TENERGY)


    
    End Subroutine DetCalcInit
    
    Subroutine DoDetCalc
      Use global_utilities
      Use HElem
      use SystemData, only : Alat, arr, brr, boa, box, coa, ecore, g1,Beta
      use SystemData, only : nBasis, nBasisMax,nEl,nMsh
      use IntegralsData, only: FCK,NMAX, UMat
      Use Logging, only: iLogging
      use SystemData, only  : tCSF
!*      POINTER (IP_TKE,TKE)
!*      REAL*8 TKE(*)
!*      REAL*8 A(*),V(*),AM(*),BM(*),T(*),WT(*),SCR(*),V2(*),WORK(*),WORK2(*),WH(*)
!*      POINTER (IP_A,A),(IP_V,V),(IP_AM,AM),(IP_BM,BM),(IP_T,T),(IP_WT,WT),(IP_SCR,SCR),(IP_ISCR,ISCR)
!*      POINTER (IP_INDEX,INDEX),(IP_LAB,LAB), (IP_NROW,NROW), (IP_V2,V2), (IP_WORK,WORK), (IP_WH,WH), (IP_WORK2,WORK2)
!*      INTEGER LAB(*),NROW(*),INDEX(*),ISCR(*)

      REAL*8 , ALLOCATABLE :: TKE(:),A(:,:),V(:),AM(:),BM(:),T(:),WT(:),SCR(:),WH(:),WORK2(:),V2(:,:),WORK(:)
      INTEGER , ALLOCATABLE :: LAB(:),NROW(:),INDEX(:),ISCR(:)

      integer :: LabTag=0,NRowTag=0,TKETag=0,ATag=0,VTag=0,AMTag=0,BMTag=0,TTag=0
      INTEGER :: WTTag=0,SCRTag=0,ISCRTag=0,INDEXTag=0,WHTag=0,Work2Tag=0,V2Tag=0,WorkTag=0
      integer :: ierr
      character(25), parameter :: this_routine = 'DoDetCalc'
      REAL*8 EXEN,GSEN,FLRI,FLSI

        INTEGER GC,I,ICMAX
        INTEGER iDeg,III,IN
        INCLUDE 'irat.inc'
        INTEGER NBLOCK
        INTEGER nKry1
        
        INTEGER J,JR
        INTEGER LSCR,LISCR
        LOGICAL tMC
        real*8 GetHElement, calct, calcmcen, calcdlwdb
! Doesn't seem to have been inited

      WRITE(6,'(1X,A,E19.3)') ' B2LIMIT : ' , B2L
      WRITE(6,*) ' NBLK : ' , NBLK 
      WRITE(6,*) ' NKRY : ' , NKRY
      WRITE(6,*) ' NEVAL : ' , NEVAL

              WRITE(6,*) ' NCYCLE : ' , NCYCLE
            WRITE(6,*) ' TCORR : ' , TCORR
      WRITE(6,*) ' TENERGY : ' , TENERGY
      IF(TCORR) THEN
        WRITE(6,*) ' *** EXCHANGE-CORRELATION HOLE WILL BE CALCULATED *** ' 
      ENDIF
      WRITE(6,*) ' IOBS : ' , IOBS 
      WRITE(6,*) ' JOBS : ' , JOBS 
      WRITE(6,*) ' KOBS : ' , KOBS 
      WRITE(6,*) ' NMSH : ' , NMSH 
      IF(IOBS.GT.NMSH.OR.IOBS.LE.0.OR.JOBS.GT.NMSH.OR.JOBS.LE.0.OR.KOBS.GT.NMSH.OR.KOBS.LE.0) THEN
        STOP ' !!! REFERENCE PARTICLE NOT IN BOX !!! '
      ENDIF

!C.. now back to the storing H
      IF(TCALCHMAT) THEN
         WRITE(6,*) "Calculating H matrix"
!C..We need to measure HAMIL and LAB first 
         ALLOCATE(NROW(NDET),stat=ierr)
         CALL LogMemAlloc('NROW',NDET,4,this_routine,NROWTag,ierr)
!*         CALL MEMORY(IP_NROW,NDET,'NROW')
         CALL IAZZERO(NROW,NDET)
         ICMAX=1
!Falsify tMC
         TMC=.FALSE.
         CALL DETHAM(NDET,NEL,NMRKS,NBASISMAX,nBasis,HAMIL,G1,LAB,NROW,.TRUE.,NMSH,FCK,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BRR)
         WRITE(6,*) ' FINISHED COUNTING '
         CALL FLUSH(6)
!C..Now we know size, allocate memory to HAMIL and LAB
         LENHAMIL=GC
         Allocate(Hamil(LenHamil), stat=ierr)
         LogAlloc(ierr, 'HAMIL', LenHamil, HElementSizeB, tagHamil)
         CALL AZZERO(HAMIL,LENHAMIL*HElementSize)
!C..
         ALLOCATE(LAB(LENHAMIL),stat=ierr)
         CALL LogMemAlloc('LAB',LenHamil,4,this_routine,LabTag,ierr)

!*         CALL MEMORY(IP_LAB,LENHAMIL/IRAT+1,'LAB')
         CALL IAZZERO(LAB,LENHAMIL)
!C..Now we store HAMIL and LAB 
         CALL DETHAM(NDET,NEL,NMRKS,NBASISMAX,nBasis,HAMIL,G1,LAB,NROW,.FALSE.,NMSH,FCK,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BRR)

         IF(BTEST(ILOGGING,7)) THEN
!C.. we write out H now
            OPEN(8,FILE='HAMIL',STATUS='UNKNOWN')
            J=0
            JR=0
!C            HMAX=-dflmax()
!C            HMIN=dflmax()
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
!C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
!C               WRITE(14,"(A,$)"),"|"
!C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
!C              WRITE(14,"(F27.20)") HAMIL(I)
!C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
!C               WRITE(14,"(A,$)"),"|"
!C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
!C               WRITE(14,"(F27.20)") HAMIL(I)

!C               IF(HAMIL(I).GT.HMAX) HMAX=HAMIL(I)
!C               IF(HAMIL(I).LT.HMIN) HMIN=HAMIL(I)
            ENDDO
            CLOSE(8)
         ENDIF
        WRITE(6,*) '<D0|H|D0>=',GETHELEMENT(IFDET,IFDET,HAMIL,LAB,NROW,NDET)
        WRITE(6,*) '<D0|T|D0>=',CALCT(NMRKS(1,IFDET),NEL,G1,NBASIS) 
!CC         CALL HAMHIST(HMIN,HMAX,LENHAMIL,NHISTBOXES)
      ENDIF
!C.. We've now finished calculating H if we were going to.
!C.. IF ENERGY CALC (for which we need to have calced H)
! 
      IF(TENERGY) THEN
         IF(NBLK.NE.0) THEN
!C..Things needed for Friesner-Pollard diagonalisation
            IF(TMC) STOP 'TMC and TENERGY set - Stopping'
            IF(HElementSize.NE.1)  STOP 'Cannot do Lanczos on Complex orbitals.'
            NKRY1=NKRY+1
            NBLOCK=MIN(NEVAL,NBLK)
            LSCR=MAX(NDET*NEVAL,8*NBLOCK*NKRY)
            LISCR=6*NBLOCK*NKRY
!C..
            WRITE(*,'(/,/,8X,64(1H*))')
            WRITE(*,'(7X," *",62X,"*")')
          WRITE(*,'(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
            WRITE(*,'(7X," *",62X,"*")')
            WRITE(*,'(7X,1X,64(1H*))')
!C..Set up memory for FRSBLKH

!*            CALL MEMORY(IP_A,NEVAL*NEVAL,'A')
            ALLOCATE(A(NEVAL,NEVAL),stat=ierr)
            CALL LogMemAlloc('A',NEVAL**2,8,this_routine,ATag,ierr)
            CALL AZZERO(A,NEVAL*NEVAL)
!C..
!C,, W is now allocated with CK
!C..
!*            CALL MEMORY(IP_V,NDET*NBLOCK*NKRY1,'V')
            ALLOCATE(V(NDET*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('V',NDET*NBLOCK*NKRY1,8,this_routine,VTag,ierr)
            CALL AZZERO(V,NDET*NBLOCK*NKRY1)
!C..   
!*            CALL MEMORY(IP_AM,NBLOCK*NBLOCK*NKRY1,'AM')
            ALLOCATE(AM(NBLOCK*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('AM',NBLOCK*NBLOCK*NKRY1,8,this_routine,AMTag,ierr)
            CALL AZZERO(AM,NBLOCK*NBLOCK*NKRY1)
!C..
!*            CALL MEMORY(IP_BM,NBLOCK*NBLOCK*NKRY,'BM')
            ALLOCATE(BM(NBLOCK*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('BM',NBLOCK*NBLOCK*NKRY,8,this_routine,BMTag,ierr)
            CALL AZZERO(BM,NBLOCK*NBLOCK*NKRY)
!C..
!*            CALL MEMORY(IP_T,3*NBLOCK*NKRY*NBLOCK*NKRY,'T')
            ALLOCATE(T(3*NBLOCK*NKRY*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('T',3*NBLOCK*NKRY*NBLOCK*NKRY,8,this_routine,TTag,ierr)
            CALL AZZERO(T,3*NBLOCK*NKRY*NBLOCK*NKRY)
!C..
!*            CALL MEMORY(IP_WT,NBLOCK*NKRY,'WT')
            ALLOCATE(WT(NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('WT',NBLOCK*NKRY,8,this_routine,WTTag,ierr)
            CALL AZZERO(WT,NBLOCK*NKRY)
!C..
!*            CALL MEMORY(IP_SCR,LSCR,'SCR')
            ALLOCATE(SCR(LScr),stat=ierr)
            CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag,ierr)
            CALL AZZERO(SCR,LSCR)
!*            CALL MEMORY(IP_ISCR,LISCR,'ISCR')
            ALLOCATE(ISCR(LIScr),stat=ierr)
            CALL LogMemAlloc('IScr',LIScr,4,this_routine,IScrTag,ierr)
            CALL IAZZERO(ISCR,LISCR)
!*            CALL MEMORY(IP_INDEX,NEVAL,'INDEX')
            ALLOCATE(INDEX(NEVAL),stat=ierr)
            CALL LogMemAlloc('INDEX',NEVAL,4,this_routine,INDEXTag,ierr)
            CALL IAZZERO(INDEX,NEVAL)
!C..
!*            CALL MEMORY(IP_WH,NDET,'WH')
            ALLOCATE(WH(NDET),stat=ierr)
            CALL LogMemAlloc('WH',NDET,8,this_routine,WHTag,ierr)
            CALL AZZERO(WH,NDET)
!*            CALL MEMORY(IP_WORK2,3*NDET,'WORK2')
            ALLOCATE(WORK2(3*NDET),stat=ierr)
            CALL LogMemAlloc('WORK2',3*NDET,8,this_routine,WORK2Tag,ierr)
            CALL AZZERO(WORK2,3*NDET)
!*            CALL MEMORY(IP_V2,NDET*NEVAL,'V2')
            ALLOCATE(V2(NDET,NEVAL),stat=ierr)
            CALL LogMemAlloc('V2',NDET*NEVAL,8,this_routine,V2Tag,ierr)
            CALL AZZERO(V2,NDET*NEVAL)
!C..Lanczos iterative diagonalising routine
            CALL FRSBLKH(NDET,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT, &
     &  SCR,ISCR,INDEX,WH,WORK2,V2,NCYCLE,B2L,.true.,.false.,.false.)

!Multiply all eigenvalues by -1.
            CALL DSCAL(NEVAL,-1.D0,W,1)
         ELSE
!C.. We splice in a non-Lanczos diagonalisin routine if NBLOCK=0
            IF(NEVAL.NE.NDET) THEN
               WRITE(6,*) "NEVAL.NE.NDET.",NEVAL,NDET," Cannot exactly diagonalize."
               STOP
            ENDIF
            WRITE(6,*) "NBLK=0.  Doing exact diagonalization."
            IF(TCALCHMAT) THEN
               ALLOCATE(WORK(4*NDET),stat=ierr)
               CALL LogMemAlloc('WORK',4*NDET,8,this_routine,WorkTag,ierr)
               ALLOCATE(WORK2(3*NDET),stat=ierr)
               CALL LogMemAlloc('WORK2',3*NDET,8,this_routine,WORK2Tag,ierr)
!*               CALL MEMORY(IP_WORK,4*NDET*HElementSize,'WORK')
!*               CALL MEMORY(IP_WORK2,3*NDET,'WORK2')
               CALL HDIAG(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,LENHAMIL,NBLOCKSTARTS,NBLOCKS,BLOCKSYM)
            ELSE
!I_P we've replaced by 0
               CALL HDIAG_NH(NDET,NBLOCKSTARTS,NBLOCKS,NEL,NMRKS,NBASISMAX,NBASIS,G1,NMSH,BRR, &
     &            FCK,NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BETA,0,ILOGGING,IFDET,ARR,BLOCKSYM)
!C.. We're not storing the energies, so we pretend we weren't asked for
!C.. them
               TENERGY=.FALSE.
            ENDIF
         ENDIF
!C..
!  Since we no longer use HAMIL or LAB, we deallocate
      LogDealloc(tagHamil)
      Deallocate(Hamil)
      DEALLOCATE(LAB)
      CALL LogMemDealloc(this_routine,LabTag)
!*      CALL FREEM(IP_LAB)
      ALLOCATE(TKE(NEVAL),stat=ierr)
      CALL LogMemAlloc('TKE',NEVAL,8,this_routine,TKETag,ierr)
!*      CALL MEMORY(IP_TKE,NEVAL,'TKE')
!C.. END ENERGY CALC
      ENDIF
      IF(TENERGY) THEN
         EXEN=CALCMCEN(NDET,NEVAL,CK,W,BETA,0.D0)
         WRITE(6,"(A,F19.9)") "EXACT E(BETA)=",EXEN
         GSEN=CALCDLWDB(IFDET,NDET,NEVAL,CK,W,BETA,0.D0)
         WRITE(6,"(A,F19.9)") "EXACT DLWDB(D0)=",GSEN
         WRITE(6,"(A,F19.9)") "GROUND E=",W(1)
      ENDIF

      call FLUSH(6)
!C.. If we're calculating rhos (for which we have to have calced H
!No longer used
!      IF(TRHOIJ) THEN
!         IF((.NOT.TENERGY).AND.(.NOT.TREADRHO)) THEN
!            WRITE(6,*) "Calculating approx RHOs"
!            WRITE(6,*) "Using Trotter decomposition? ",TTROT
!            WRITE(6,*) "Order of Taylor: ",ABS(NTAY)
!            CALL CALCAPPROXRHOIJ(BETA,I_P,HAMIL,LAB,NROW,NDET,RHOMIN,RHOMAX,NRHOS,RHOEPS,TTROT,NTAY)
!         ENDIF
!      ENDIF
!C..Free HAMIL AND LAB memory if we no longer need them
!      IF(TCALCHMAT.AND..NOT.(TMONTE.AND.TMC)) THEN
!         CALL FREEM(IP_HAMIL)
!         CALL FREEM(IP_LAB)
!      ENDIF

!C.. IF ENERGY CALC
      IF (TENERGY) THEN
!C..
         IF(.NOT.TCSF) THEN
            CALL CFF_CHCK(NDET,NEVAL,NMRKS,NBASISMAX,NEL,G1,CK,ALAT,TKE,nBasis,ILOGGING)
         ELSE
            DO I=1,NEVAL
               TKE(I)=0.D0
            ENDDO 
         ENDIF
         IF(BTEST(ILOGGING,7)) CALL WRITE_PSI(BOX,BOA,COA,NDET,NEVAL,NBASISMAX,NEL,CK,W)
         IF(BTEST(ILOGGING,8)) CALL WRITE_PSI_COMP(BOX,BOA,COA,NDET,NEVAL,NBASISMAX,NEL,CK,W)
         WRITE(6,*) '       ==--------------------------------------------------== '
         WRITE(6,'(A5,5X,A15,1X,A18,1x,A20)') 'STATE','KINETIC ENERGY', 'COULOMB ENERGY', 'TOTAL ENERGY'
         OPEN(15,FILE='ENERGIES',STATUS='UNKNOWN')
         DO IN=1,NEVAL
            WRITE(6,'(I5,2X,3(F19.11,2x))') IN,TKE(IN),W(IN)-TKE(IN),W(IN)
!            WRITE(15,"(I7,$)") IN
!            CALL WRITEDET(15,NMRKS(1,IN),NEL,.FALSE.)
            WRITE(15,"(F19.11)") W(IN)
         ENDDO
         CLOSE(15)
         WRITE(6,*)   '       ==--------------------------------------------------== '
!C., END energy calc
      ENDIF

!C.. Jump to here if just read Psi in
100   CONTINUE

      IF(TRHOOFR) THEN
        Call CalcRhoOfR()
      ENDIF
      IF(TFODM) THEN
        Call CalcFoDM()
      ENDIF


        End Subroutine DoDetCalc
    Subroutine CalcRhoOfR()
        use SystemData, only: Alat, G1, nBasis, Omega, nEl,nMsh
        use IntegralsData, only: nMax
        REAL*8 DLINE(*)
        POINTER (IP_DLINE, DLINE)
        REAL*8 PSIR(*)
        POINTER (IP_PSIR, PSIR)
        REAL*8 RHO(*)
        POINTER (IP_RHO, RHO)
        REAL*8 SITAB(*)
        POINTER (IP_SITAB,SITAB)
        REAL*8 SCRTCH(*)
        POINTER (IP_SCRTCH, SCRTCH)
        REAL*8 XCHOLE(*)
        POINTER (IP_XCHOLE,XCHOLE)
        INTEGER iXD, iYD, iZD
        REAL*8 SPAC, Rs
!C..Generate memory for RHO and SITAB
        CALL MEMORY(IP_RHO,NMSH*NMSH*NMSH,'RHO')
        CALL MEMORY(IP_SITAB,NMSH*NMAX,'SITAB')
        CALL MEMORY(IP_SCRTCH,NMSH*NMSH*NMSH,'SCRTCH')
!C..Calculate RHOOFR
        CALL RHOOFR(nBasis,CK,G1,RHO,NMSH,SITAB,NMAX,NMRKS,NEL,NDET,NEVAL,RS,ALAT,OMEGA,SCRTCH)
!C..
        CALL MEMORY(IP_DLINE,NMSH,'DLINE')
        CALL AZZERO(DLINE,NMSH)
!C..Calculate RHOOFR in certain directions
!C..001
        CALL PLANARAV(RHO,NMSH,DLINE,0,0,1,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV001',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
!C..100
        CALL PLANARAV(RHO,NMSH,DLINE,1,0,0,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV100',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
!C..010
        CALL PLANARAV(RHO,NMSH,DLINE,0,1,0,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV010',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
        IF(TCORR) THEN
!C..Now generate memory for XCHOLE
          CALL MEMORY(IP_XCHOLE,NMSH*NMSH*NMSH,'XCHOLE')
          CALL MEMORY(IP_PSIR,2*NMSH+1,'PSIR')
          CALL AZZERO(PSIR,2*NMSH+1)
!C..
          IXD=1
          IYD=0
          IZD=0
          SPAC=0.D0
          CALL GEN_XCHOLE(CK,PSIR,IOBS,JOBS,KOBS,G1,SITAB,NMAX,NMSH,nBasis,IXD,IYD,IZD,RHO,.TRUE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,NEL)
          CALL WRITE_RHO(10,'COMPXCHOLE',XCHOLE,NMSH,NMSH,NMSH,ALAT,.FALSE.,.TRUE.,RS)
!C..
          CALL XCHOLES(CK,PSIR,IOBS,JOBS,KOBS,G1,SITAB,NMAX,NMSH,nBasis,RHO,XCHOLE,SPAC,RS,ALAT,DLINE,OMEGA,NMRKS,NDET,NEL,NEVAL)
        ENDIF
    End Subroutine CalcRhoOfR
    Subroutine CalcFoDM()
        use SystemData, only: G1, nBasis, nMaxX, nMaxY, nMaxZ, nEl
        REAL*8 SUMA
        POINTER (IP_SUMA, SUMA)
        INTEGER ISTATE
        ISTATE=0
        WRITE(6,*) ' ISTATE : ' , ISTATE 
        CALL MEMORY(IP_SUMA,NMAXX*NMAXY*NMAXZ,'SUMA')
        CALL AZZERO(SUMA,NMAXX*NMAXY*NMAXZ)
        CALL FODMAT(NEL,NBasis,NDET,NEVAL,ISTATE,NMRKS,G1,CK,NMAXX,NMAXY,NMAXZ,SUMA)
    End Subroutine CalcFoDM
END MODULE DetCalc

! If we have a list of determinants NMRKS calculate 'PATHS' for NPATHS of them.
      SUBROUTINE CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET, &
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,  &
     &   NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,   &
     &   DETINV,TSPECDET,SPECDET)
         use HElem
         use SystemData, only: BasisFN
         use global_utilities
         implicit none
         include 'irat.inc'
         INTEGER NEL,I_P,I_HMAX,I_VMAX,NDET,nBasisMax(5,*),nBasis
         INTEGER BRR(*),NMSH,NMAX(*),NTAY,ILOGGING
         type(timer), save :: proc_timer
         TYPE(HElement) UMat(*)
         TYPE(HDElement) DLWDB, DLWDB2, DLWDB3, DLWDB4
         TYPE(BasisFN) g1(*),ALAT(*)
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),RHOEPS
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
         DLWDB2=0.D0
         IMAX=I_HMAX
         IF(I_VMAX.GT.IMAX) IMAX=I_VMAX
         proc_timer%timer_name='CLCRHOPII2'
         call set_timer(proc_timer)
         ILMAX=NDET
!.. we don't need lists for I_HMAX=8
         IF((I_HMAX.GE.-10.AND.I_HMAX.LE.-7)      .OR.I_HMAX.LE.-12) ILMAX=1
!         ILMAX=(NBASIS-NEL)**2*NEL*NEL/4
         CALL MEMORY(IP_LSTE,(1+ILMAX)*NEL*IMAX/IRAT,"LSTE")
         CALL MEMORY(IP_ICE,(1+ILMAX)*IMAX/IRAT,"ICE")
        CALL MEMORY(IP_RIJLIST,(1+ILMAX)*IMAX*2,"RIJLIST")
         IF(I_VMAX.NE.0) THEN
            WRITE(6,*) "Using Vertex approximation.  I_VMAX=",I_VMAX
            IF(I_HMAX.EQ.0) WRITE(6,*) "I_HMAX=0.  Summing all I_HMAX up to P using contour"
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
               CALL MCPATHSR3(NMRKS(1,III),BETA,I_P,I_HMAX,I_VMAX,NEL, NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT, &
     &         UMAT,NTAY, RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX, WLRI,WLSI,DBETA,DLWDB2)
            ELSE
               CALL MCPATHSR3(SPECDET,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
     &         RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX, WLRI,WLSI,DBETA,DLWDB2)
            ENDIF
            WRITE(15,"(I12,$)") III
            IF(TSPECDET) THEN
               CALL WRITEDET(15,SPECDET,NEL,.FALSE.)
            ELSE
               CALL WRITEDET(15,NMRKS(1,III),NEL,.FALSE.)
            ENDIF
            WRITE(15,"(A,3G25.16,$)") " ",EXP(WLSI+I_P*WLRI),WLRI*I_P,WLSI
            IF(III.EQ.1) THEN
               WLRI0=WLRI
               WLSI0=WLSI
           else
               WLRI0=0.d0
               WLSI0=0.d0
            ENDIF  
            IF(TNPDERIV) THEN
!.. if we're calculating the derivatives too
               IF(III.NE.0) THEN
               CALL MCPATHSR3(NMRKS(1,III),BETA+DBETA,I_P,I_HMAX,    &
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,       &
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,&
     &            ILOGGING,ECORE,ILMAX,WLRI1,WLSI1,DBETA,DLWDB3)
               CALL MCPATHSR3(NMRKS(1,III),BETA-DBETA,I_P,I_HMAX,    &
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,       &
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,&
     &            ILOGGING,ECORE,ILMAX,WLRI2,WLSI2,DBETA,DLWDB4)
               ELSE
               CALL MCPATHSR3(SPECDET,BETA+DBETA,I_P,I_HMAX,         &
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,       &
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,&
     &            ILOGGING,ECORE,ILMAX,WLRI1,WLSI1,DBETA,DLWDB3)
               CALL MCPATHSR3(SPECDET,BETA-DBETA,I_P,I_HMAX,         &
     &            I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,       &
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,&
     &            ILOGGING,ECORE,ILMAX,WLRI2,WLSI2,DBETA,DLWDB4)
               ENDIF
               DLWDB=-(I_P*(WLRI1-WLRI2)+(WLSI1-WLSI2))/(2*DBETA)
            ELSE
               DLWDB=DLWDB2
            ENDIF
!.. we calculate the energy with weightings normalized to the weight of
!.. the Fermi determinant, otherwise the numbers blow up
            WINORM=EXP(I_P*(WLRI-WLRI0)+(WLSI-WLSI0))
            NORM=NORM+WINORM
            TOT=TOT+WINORM*DREAL(DLWDB)
            WRITE(15,*) DLWDB
            IF(DETINV.EQ.III.AND.III.NE.0) THEN
               CALL FLUSH(15)
               WRITE(6,*) "Investigating det ",DETINV
               CALL FLUSH(6)
               CALL WIRD_SUBSET(NMRKS(1,DETINV),BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
     &       RHOEPS,ILOGGING,TSYM,ECORE)
            ENDIF
          ENDDO
         CLOSE(15)
         WRITE(6,*) "Summed approx E(Beta)=",TOT/NORM
         CALL FREEM(IP_RIJLIST)
         CALL FREEM(IP_LSTE)
         CALL FREEM(IP_ICE)
         call halt_timer(proc_timer)
         RETURN
      END    

!   Using the exact eigenvectors and -values calculate the exact value of rho^P_ii

!.. Note if TWARN becomes set, the RHII sum did not converge.
!.. FLSI will remain usable, but will be equal to log RHII(P), so the 
!.. sum I_P*FLRI+FLSI will still retain the correct value.
      SUBROUTINE CALCRHOPII(I,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ETRIAL,FLRI,FLSI,TWARN)
         USE HElem
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
!.. We decompose ln(RHO^(P)_II) = p ln RHO_II + ln sI

!.. First we work out RHO_II
!         WRITE(6,*) BETA,I_P
         DO IK=1,NEVAL
            R=SQ(CK(I,IK))
            R=R*EXP(-(W(IK)-W(1))*BETA/I_P)
            RHII=RHII+R
         ENDDO
         IF(NEVAL.NE.NDET.AND.EXP(-(W(NEVAL)-W(1))*BETA/I_P).GT.1.D-2) THEN
!.. If we don't have all our eigenvectors and we think our sum has not
!.. converged, we print a warning the first time.
!.. we cannot calculate a proper RHII so we just guess at 1
            IF(I.EQ.1) THEN
             WRITE(6,*) 'WARNING: For Det 1 RHII sum has not converged.'
             WRITE(6,*) ' Setting RHII=1'
            ENDIF
            TWARN=.TRUE.
            RHII=1.D0
            FLRI=0.D0
         ELSE

!.. and Log it
           FLRI=LOG(RHII)-(W(1))*BETA/I_P
         ENDIF
!.. Now we work out RHO^(P)_II/RHO_II^P = sI         
         DO IK=1,NEVAL
            R=SQ(CK(I,IK))
            RH=RH+R*EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         FLSI=LOG(RH)-W(1)*BETA-I_P*FLRI
         IF(LISNAN((RH+1)-RH)) THEN
            RH=0
            FLSI=0
         ENDIF
!         WRITE(17,*) CK(1,I),W(I)
         RETURN
      END

!  Given an exact calculation of eigen-vectors and -values, calculate the expectation value of E(Beta)
      REAL*8 FUNCTION CALCMCEN(NDET,NEVAL,CK,W,BETA,ETRIAL)
         USE HElem
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

!  Given an exact calculation of eigen-vectors and -values, calculate the expectation value of E~(Beta)_I for det I
      REAL*8 FUNCTION CALCDLWDB(I,NDET,NEVAL,CK,W,BETA,ETRIAL)
         USE HElem
         IMPLICIT NONE
         INTEGER NDET,NEVAL,IK,I
         TYPE(HElement) CK(NDET,NEVAL)
         REAL*8  W(NEVAL),BETA,DNORM,EN,ETRIAL
         EN=0.D0
         DNORM=0.D0
         DO IK=1,NEVAL
            EN=EN+SQ(CK(I,IK))*(W(IK))*EXP(-(W(IK)-W(1))*BETA)
            DNORM=DNORM+SQ(CK(I,IK))*EXP(-(W(IK)-W(1))*BETA)
         ENDDO
         CALCDLWDB=EN/DNORM
         RETURN
      END

      SUBROUTINE CFF_CHCK(NDET,NEVAL,NM,NBASISMAX,NEL,G1,CG,ALAT,TKE,NHG,ILOGGING)
      USE HElem
      USE OneEInts, only : GetTMATEl
      use SystemData, only: BasisFN
      IMPLICIT NONE
      TYPE(HElement) CG(NDET,NEVAL)
      INTEGER NM(NEL,*),NDET,NEL,NEVAL
      REAL*8 ALAT(3),TKE(NEVAL)
      INTEGER NBASISMAX(3,2),NHG,ILOGGING
      TYPE(BASISFN) G1(*)
      CHARACTER*255 STR
      REAL*8 PI,S,SUM
      TYPE(HDElement) AUX
      INTEGER I,J,IN,IEL,L
!..Calculate the expectation value of the kinetic energy
!..<Psi|T|Psi>
      PI=ACOS(-1.D0)
      DO IN=1,NEVAL
        TKE(IN)=0.D0
        DO I=1,NDET
          SUM=0.D0
          DO J=1,NEL
            AUX=GetTMATEl(NM(J,I),NM(J,I))
!((ALAT(1)**2)*((G1(1,NM(J,I))**2)/(ALAT(1)**2)+
!     &		(G1(2,NM(J,I))**2)/(ALAT(2)**2)+
!     &		(G1(3,NM(J,I))**2)/(ALAT(3)**2)))
            SUM=SUM+DREAL(AUX)
          ENDDO
!..Cube multiplier
!          CST=PI*PI/(2.D0*ALAT(1)*ALAT(1))
!.. Deal with the UEG
!          IF(NBASISMAX(1,1).LE.0) CST=CST*4.D0
!          SUM=CST*SUM 
          TKE(IN)=TKE(IN)+SUM*SQ(CG(I,IN))
        ENDDO
      ENDDO
! ==--------------------------------------------------------------==
      IF(.FALSE.) THEN
!      IF(BTEST(ILOGGING,7)) THEN
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


!Given exact eigenvalues and vectors, do monte carlo in det space with exact weights and E~
       REAL*8 FUNCTION DOEXMC(NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ECORE,IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,NBASIS,BRR,IEQSTEPS)
         use HElem, only: HElement
         INTEGER NDET,NEVAL,I_P,ILOGGING
         type(HElement) CK(NEVAL)
         REAL*8 W(NEVAL),BETA,ECORE

         REAL*8 DLWDBS(NDET),WLRIS(NDET),WLSIS(NDET),EN
         INTEGER I
         LOGICAL TWARN
         
         DO I=1,NDET
            CALL CALCRHOPII(I,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,ECORE,WLRIS(I),WLSIS(I),TWARN)
            DLWDBS(I)=CALCDLWDB(I,NDET,NEVAL,CK,W,BETA,ECORE)
         ENDDO
         EN=DMONTECARLOEXWI(NDET,WLRIS,WLSIS,DLWDBS,I_P,IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,NBASIS,BRR,IEQSTEPS,ILOGGING)
         WRITE(6,*) "EXACT MC RESULT=",EN
         DOEXMC=EN
         RETURN
      END

