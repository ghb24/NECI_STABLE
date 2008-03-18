#include "macros.h"
MODULE Determinants
    Use HElem
    implicit none
    save
! Set by Calc on input
      INTEGER nActiveSpace(2)
        INTEGER, DIMENSION(:), POINTER :: SPECDET
        Logical TSPECDET

      INTEGER, pointer :: FDet(:)
      integer tagFDet
!nActiveBasis(1) is the lowest non-active orbital
!nActiveBasis(2) is the highest active orbital.  There can be virtuals above this.
!  Active orbitals are used for generating the determinants whose energy/weight is to be found
      Integer nActiveBasis(2)
!  Set by input to indicate which type of active basis we need
      Integer iActiveBasis
!Not really Locals - needed for the DetCalc
      INCLUDE 'basis.inc'
      TYPE(BasisFN) ISym
!Used to be from uhfdet.inc
      INTEGER nUHFDet(5000)
      REAL*8  E0HFDet

     save FDet
      
    contains
    Subroutine DetPreFreezeInit()
        Use MemoryManager, only: LogMemAlloc, LogMemDealloc
        Use System, only : nEl, ECore, Arr, Brr, G1, nBasis, LMS
        integer ierr
        integer i
   
        character(25), parameter :: this_routine='DetPreFreezeInit'
        Allocate(FDet(nEl), stat=ierr)
        LogAlloc(ierr, 'FDet', nEl, 4, tagFDet)
         CALL GENFDET(BRR,G1,NBASIS,LMS,NEL,FDET)
!      ENDIF
      WRITE(6,"(A,$)") "Fermi det (D0):"
      CALL WRITEDET(6,FDET,NEL,.TRUE.)
      CALL ICOPY(NEL,FDET,1,NUHFDET,1)
      E0HFDET=ECORE
      DO I=1,NEL
         E0HFDET=E0HFDET+ARR(NBASIS+NUHFDET(i),1)
      ENDDO     
      WRITE(6,*) "Fock operator energy:",E0HFDET
    End Subroutine DetPreFreezeInit
    Subroutine DetInit()
        Use MemoryManager, only: LogMemAlloc, LogMemDealloc
        Use HElem
        Use System, only: nel, Alat, Boa, Coa, BOX, BRR, ECore
        Use System, only: G1, LMS, nBasis, STot, tCSF, Arr
        use integrals, only: nfrozen
      
      real*8 DNDET
      integer i,ii
      integer*8 nDet
      integer ierr
      character(25), parameter :: this_routine='DetInit'

      WRITE(6,*) "SYMMETRY MULTIPLICATION TABLE"
      CALL WRITESYMTABLE(6)
   
      CALL GENSymStatePairs(NBASIS/2,.false.)


!iActiveBasis is a copy of nPaths
      IF(iActiveBasis.eq.-2) then
!  PATHS ACTIVE SETS
         Call GenActiveBasis(ARR,BRR,G1,nBasis,LMS,nEl,nActiveBasis,nActiveSpace(1),nActiveSpace(2))
      elseif(iActiveBasis.eq.-3) then
!  PATHS ACTIVE ORBITALS
         nActiveBasis(1)=nEl+1-nActiveSpace(1)
         nActiveBasis(2)=nEl+nActiveSpace(2)
         WRITE(6,*) "Active space:", nActiveBasis(1)," TO ", nActiveBasis(2)," (ordered labels)."
         WRITE(6,*) "Active space electrons:",nEl-nActiveBasis(1)+1
      else
         nActiveBasis(1)=1
         nActiveBasis(2)=nBasis
      endif
!C.. Work out a preliminary Fermi det
!      IF(FDET(1).EQ.0) THEN

 

!C.. Check if we're blocking the hamiltonian
!C      IF(THFBASIS.AND.TBLOCK) THEN
!C         WRITE(6,*) "THFBASIS set and NBLK=0.  ",
!C     &         "Cannot block diagonalize in HF Basis."
!C         STOP
!C      ENDIF
!C      CALL SYMGENEXCITS(FDET,NEL,G1,NBASIS,NBASISMAX)
!C      CALL LeaveMemoryManager
!C      CALL TIPRI
!C      STOP


!C.. in order to calculate the H matrix, we need to work out all the determinants
!C.. beware with NPATHS - it calcs the list of dets even if we don't calc H
!C.. Could be big.
!C..Now we see how many determinants we need
!C      IF(nBasis.GT.170) THEN
!C..This fix is to stop floating overflow as FACTRL(nBasis.GT.170) crashes
         NDET=1
         DNDET=1.D0
         DO I=0,NEL-1
            NDET=(NDET*(nBasis-I))/(I+1)
            DNDET=(DNDET*DFLOAT(nBasis-I))/DFLOAT(I+1)
         ENDDO
!C      ELSE
!C         NDET=FACTRL(nBasis)/(FACTRL(NEL)*FACTRL(nBasis-NEL)) 
!C      ENDIF
        IF(NDET.ne.DNDET) THEN
         WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , DNDET
         NDET=-1
        ELSE
         WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , NDET
        ENDIF
      
!C      CALL TC(I_HMAX,I_P,NWHTAY)


    End Subroutine DetInit
    

!.. GETHELEMENT2
!.. Get matrix element of the hamiltonian
!.. IC is the number of basis fns that differ in NI and NJ (or -1 if not known)
!.. ECORE is the uniform background energy

      TYPE(HElement) FUNCTION GetHElement2(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,iC2,ECore)
         Use HElem
         USE System , only : TSTOREASEXCITATIONS
         IMPLICIT NONE
         INTEGER NMSH,NMAX
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(*)
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS,BRR(*)
         TYPE(HElement) UMat(*)
         INTEGER I,nEl,NI(nEl),NJ(nEl),iC,nBasisMax(5,5),iC2
         REAL*8 ECore
         TYPE(HElement) Sum,Sum2
         INTEGER IGETEXCITLEVEL
         LOGICAL ISCSF
         INTEGER ISUB
         IF(ISCSF(NI,NEL).OR.ISCSF(NJ,NEL)) THEN
            CALL CSFGETHELEMENT(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,ECore,Sum2)
            GETHELEMENT2=SUM2
            RETURN
         ENDIF
         IF(tStoreAsExcitations.AND.nI(1).eq.-1.and.nJ(1).eq.-1) then
            if(ic2.ne.2) stop 'tStoreAsExcitations in GetHElement2 requires ic=2 (doubles).'
            Call SCR2Excit(nBasisMax,nJ,G1,nBasis,UMat,Alat,nBasisMax(2,3),Sum)
            GetHElement2=Sum
            RETURN
         endif
         IC=IC2
         GetHElement2%v=0.D0
         IF(IC.LT.0) IC=IGETEXCITLEVEL(NI,NJ,NEL)
!.. if we differ by more than 2 spin orbital, then the hamiltonian element is 0         
         IF(IC.GT.2) RETURN
!.. SLTCND has IC is # electrons the same in 2 dets
         CALL TISETL('GETHELEM2 ',ISUB,60)
         CALL SltCnd(nEl,nBasisMax,nBasis,NI,NJ,G1,nEl-iC,NMSH,FCK,NMAX,ALAT,UMat,Sum)
         GetHElement2=Sum
         IF(iC.EQ.0) GetHElement2%v=GetHElement2%v+ECore
!         CALL WRITEDET(6,NI,NEL,.FALSE.)
!         CALL WRITEDET(6,NJ,NEL,.FALSE.)
!         WRITE(6,*) GetHElement2
         CALL TIHALTL('GETHELEM2 ',ISUB,60)
         RETURN
      END FUNCTION
      Subroutine DetCleanup()
      End Subroutine DetCleanup
END MODULE Determinants
      subroutine GetH0Element(nI,nEl,Arr,nBasis,ECore,hEl)
         USE System , only : TSTOREASEXCITATIONS
         USE HElem
         implicit none
         integer nI(nEl),nEl,nBasis
         type(HElement) hEl
         real*8 Arr(nBasis,2),ECore
         integer i
         if(tStoreAsExcitations.and.nI(1).eq.-1) then
!The excitation storage starts with -1.  The next number is the excitation level,L .  
!Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
            hEl=0.d0
            do i=4,nI(2)+4-1
               hEl=hEl-HElement(Arr(nI(i),2))
            enddo
            do i=i,i+nI(2)-1
               hEl=hEl+HElement(Arr(nI(i),2))
            enddo
         else
            hEl=ECore
            do i=1,nEl
               hEl=hEl+HElement(Arr(nI(i),2))
            enddo
         endif
!         call writedet(77,nI,nel,.false.)
!         write(77,*) "H0",hEl
!         call flush(77)
      end subroutine
!  Get a matrix element of the unperturbed Hamiltonian.  This is just the sum of the Hartree-Fock eigenvalues
      subroutine DetFreezeBasis(GG)
        Use Determinants, only: FDet
        Use System, only : nEl, nBasis
        Use Integrals, only : nFrozen
        implicit none
        integer i,j
        INTEGER GG(*)
!C.. Deal with FDET
!C.. GG(I) is the new position in G of the (old) orb I
         IF(FDET(1).NE.0) THEN
            J=0
            DO I=1,NEL
               FDET(I)=GG(FDET(I))
!C.. any orbitals which no longer exist, we move outside the basis
               IF(FDET(I).EQ.0) THEN
                  FDET(I)=nBasis+1
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
      end subroutine
