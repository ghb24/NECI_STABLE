#include "macros.h"
MODULE Determinants
    Use HElem
    use SystemData, only: BasisFN
    implicit none
    save
! Set by Calc on input
      INTEGER nActiveSpace(2)
        INTEGER, DIMENSION(:), POINTER :: SPECDET
        INTEGER :: tagSPECDET=0
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
      TYPE(BasisFN) ISym
!Used to be from uhfdet.inc
      INTEGER nUHFDet(5000)
      REAL*8  E0HFDet

    contains
    Subroutine DetPreFreezeInit()
        Use global_utilities
        use SystemData, only : nEl, ECore, Arr, Brr, G1, nBasis, LMS
        integer ierr
        integer i
   
        character(25), parameter :: this_routine='DetPreFreezeInit'
        Allocate(FDet(nEl), stat=ierr)
        LogAlloc(ierr, 'FDet', nEl, 4, tagFDet)
         CALL GENFDET(BRR,G1,NBASIS,LMS,NEL,FDET)
!      ENDIF
      WRITE(6,"(A)",advance='no') "Fermi det (D0):"
      CALL WRITEDET(6,FDET,NEL,.TRUE.)
      CALL NECI_ICOPY(NEL,FDET,1,NUHFDET,1)
      E0HFDET=ECORE
      DO I=1,NEL
         E0HFDET=E0HFDET+ARR(NBASIS+NUHFDET(i),1)
      ENDDO     
      WRITE(6,*) "Fock operator energy:",E0HFDET
    End Subroutine DetPreFreezeInit
    Subroutine DetInit()
        Use global_utilities
        Use HElem
        use SystemData, only: nel, Alat, Boa, Coa, BOX, BRR, ECore
        use SystemData, only: G1, LMS, nBasis, STot, tCSF, Arr
        use IntegralsData, only: nfrozen
      
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
!C      STOP


!C.. in order to calculate the H matrix, we need to work out all the determinants
!C.. beware with NPATHS - it calcs the list of dets even if we don't calc H
!C.. Could be big.
!C..Now we see how many determinants we need
!C      IF(nBasis.GT.170) THEN
!C..This fix is to stop floating overflow as taking the factorial of (nBasis.GT.170) crashes
!C  using the old FACTRL routine.
         NDET=1
         DNDET=1.D0
         DO I=0,NEL-1
            NDET=(NDET*(nBasis-I))/(I+1)
            DNDET=(DNDET*DFLOAT(nBasis-I))/DFLOAT(I+1)
         ENDDO
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
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: BasisFN
         use global_utilities
         IMPLICIT NONE
         INTEGER NMSH,NMAX
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(*)
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS,BRR(*)
         TYPE(HElement) UMat(*)
         INTEGER I,nEl,NI(nEl),NJ(nEl),iC,nBasisMax(5,*),iC2
         REAL*8 ECore
         TYPE(HElement) Sum,Sum2
         INTEGER IGETEXCITLEVEL_2
         LOGICAL ISCSF
         type(timer), save :: proc_timer
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
         IF(IC.LT.0) THEN
             IC=IGETEXCITLEVEL_2(NI,NJ,NEL,2)   !Calculate whether connected or not
         ENDIF
!.. if we differ by more than 2 spin orbital, then the hamiltonian element is 0         
         IF(IC.GT.2) RETURN
!.. SLTCND has IC is # electrons the same in 2 dets
         proc_timer%timer_name='GETHELEM2 '
         call set_timer(proc_timer,60)
         CALL SltCnd(nEl,nBasisMax,nBasis,NI,NJ,G1,nEl-iC,NMSH,FCK,NMAX,ALAT,UMat,Sum)
         GetHElement2=Sum
         IF(iC.EQ.0) GetHElement2%v=GetHElement2%v+ECore
!         CALL WRITEDET(6,NI,NEL,.FALSE.)
!         CALL WRITEDET(6,NJ,NEL,.FALSE.)
!         WRITE(6,*) GetHElement2
         call halt_timer(proc_timer)
         RETURN
      END FUNCTION
!Call GetHElement2 without needing so many arguments
      TYPE(HElement) FUNCTION GetHElement3(NI,NJ,iC)
         Use HElem
         use SystemData, only : nEl,nBasisMax,G1,nBasis,Brr
         use SystemData, only : ECore,ALat,NMSH
         use IntegralsData, only : UMat,FCK,NMAX
         INTEGER NI(nEl),NJ(nEl),iC
         GetHElement3=GetHElement2(NI,NJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,iC,ECore)
      END Function GetHElement3

      type(HElement) function GetH0Element3(nI)
         ! Wrapper for GetH0Element.
         ! Returns the matrix element of the unperturbed Hamiltonian, which is
         ! just the sum of the eigenvalues of the occupied orbitals and the core
         ! energy.
         !  Note that GetH0Element{1,2} don't exist. The name is to be
         !  consistent with GetHElement3, i.e. offer the most abstraction possible.
         ! In: 
         !    nI(nEl)  list of occupied spin orbitals in the determinant.
         use HElem
         use SystemData, only: nEl,nBasis,Arr,ECore
         integer nI(nEl)
         type(HElement) hEl
         call GetH0Element(nI,nEl,Arr(1:nBasis,1:2),nBasis,ECore,hEl)
         GetH0Element3=hEl
      end function

      Subroutine DetCleanup()
      End Subroutine DetCleanup
END MODULE Determinants

      subroutine GetH0Element(nI,nEl,Arr,nBasis,ECore,hEl)
         !  Get a matrix element of the unperturbed Hamiltonian.  This is just
         !  the sum of the Hartree-Fock eigenvalues and the core energy.
         !  In:
         !     nI(nEl)  list of occupied spin orbitals in the determinant.
         !     nEl      # of electrons.
         !     Arr      array containing the eigenvalues of the spin-orbitals.
         !              (See System for how it's defined/stored).
         !     nBasis   # spin orbitals.
         !     ECore    Core energy.
         !  Out:
         !     hEl      <D_i|H_0|D_i>, the unperturbed Hamiltonian matrix element.
         use SystemData , only : TSTOREASEXCITATIONS
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

      subroutine DetFreezeBasis(GG)
        Use Determinants, only: FDet, nUHFDet
        use SystemData, only : nEl, nBasis
        use IntegralsData, only : nFrozen
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
            CALL NECI_SORTI(NEL,FDET)
            IF(J.NE.NEL-NFROZEN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               CALL WRITEDET(6,NEL,FDET,.TRUE.)
               STOP "After Freezing, FDET has wrong number of electrons"
            ENDIF
         ENDIF
         IF(nUHFDet(1).NE.0) THEN
            J=0
            DO I=1,NEL
               nUHFDET(I)=GG(nUHFDET(I))
!C.. any orbitals which no longer exist, we move outside the basis
               IF(nUHFDET(I).EQ.0) THEN
                  nUHFDET(I)=nBasis+1
               ELSE
                  J=J+1
               ENDIF
            ENDDO
            CALL NECI_SORTI(NEL,nUHFDET)
            IF(J.NE.NEL-NFROZEN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               CALL WRITEDET(6,NEL,nUHFDET,.TRUE.)
               STOP "After Freezing, UHFDET has wrong number of electrons"
            ENDIF
         ENDIF
         WRITE(6,"(A)",advance='no') "Post-Freeze Fermi det (D0):"
         CALL WRITEDET(6,FDET,NEL-NFROZEN,.TRUE.)
      end subroutine


      LOGICAL FUNCTION ISUHFDET(NI,NEL)
         use SystemData , only : TUSEBRILLOUIN
         Use Determinants, only : NUHFDET
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         INTEGER I
         ISUHFDET=.FALSE.
         IF(.NOT.TUSEBRILLOUIN) RETURN
            ISUHFDET=.TRUE.
            DO I=1,NEL
               IF(NI(I).NE.NUHFDET(I)) ISUHFDET=.FALSE.
            ENDDO
!         ISUHFDET=.FALSE.
         RETURN
      END Function


! Generate the active space from a basis.
! The Active basis can be used to in PATHS calculations and (?) as a CASCI

! nActiveBasis(1:2) contains (First Active Basis Fn, Last Active Basis Fn)
! nDown is the number of orbital sets  below the Fermi level
! nUp is the number of orbital sets  above the Fermi level

      SUBROUTINE GenActiveBasis(ARR,BRR,G1,nBasis,LMS,nEl,nActiveBasis, nDown,nUp)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         REAL*8 ARR(nBasis)
         INTEGER BRR(nBasis)
         TYPE(BasisFN) G1(nBasis)
         INTEGER LMS,nEl,nActiveBasis(2),nBasis
         INTEGER I,nDown,nUp,nLeft
         I=nEl+1
         nLeft=1+nUp
         IF(nDown.NE.0.AND.nUp.NE.0) WRITE(6,*) "Including ",-nDown,",",nUp," extra degenerate sets in active space."
         DO WHILE (nLeft.GT.0.AND.I.LT.nBasis)
            DO WHILE (I.LT.nBasis.AND.ABS(ARR(I)-ARR(I-1)).LT.1.d-5)
               I=I+1
            ENDDO
            nLeft=nLeft-1
            IF(nLeft.EQ.nUp.AND.I.NE.nEl+1) WRITE(6,*) "Fermi determinant degenerate.  "
            IF(nLeft.ne.0) I=I+2
         ENDDO
         IF(I.EQ.nEl+1.and.nDown.eq.0) THEN
!We have no degeneracies at the Fermi Energy
            WRITE(6,*) "Fermi determinant non-degenerate.  "
            IF(nDown.eq.0) THEN
               WRITE(6,*) "Active space empty."
               nActiveBasis(1)=nEl+1
               nActiveBasis(2)=nEl
               RETURN
            ENDIF
         ENDIF
         nActiveBasis(2)=I-1
         I=nEl-1
         nLeft=nDown
         Do WHILE(nLeft.GT.0.AND.I.Gt.0)
      
            DO WHILE (I.GT.0.AND.ABS(ARR(I)-ARR(I+1)).LT.1.d-5)
               I=I-1
            ENDDO
            nLeft=nLeft-1
         ENDDO
         nActiveBasis(1)=I+1
         WRITE(6,*) "Active space:", nActiveBasis(1)," TO ",nActiveBasis(2)," (ordered labels)."
         WRITE(6,*) "Active space electrons:",nEl-nActiveBasis(1)+1
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
         CALL NECI_SORTI(NEL,MCDET)
         RETURN
      END

! Write determinant NI(NEL) to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
      SUBROUTINE WRITEDET(NUNIT,NI,NEL,LTERM)
         IMPLICIT NONE
         INTEGER NUNIT,NEL,NI(NEL),I
         LOGICAL LTERM
         INTEGER IEL
         CHARACTER*2 SUFF
         INCLUDE 'csf.inc'
         WRITE(NUNIT,"(A)",advance='no') "("
         DO I=1,NEL
            IEL=NI(I)
            IF(IEL.GE.CSF_NBSTART) THEN
               WRITE(NUNIT,"(I3)",advance='no'),(IEL-CSF_NBSTART)/4+1
               IEL=IAND(IEL-CSF_NBSTART,3)
               IF(IEL.EQ.0) THEN
                  WRITE(NUNIT,"(A)",advance='no') "-B,"
               ELSEIF(IEL.EQ.1) THEN
                  WRITE(NUNIT,"(A)",advance='no') "-A,"
               ELSEIF(IEL.EQ.2) THEN
                  WRITE(NUNIT,"(A)",advance='no') "+B,"
               ELSE
                  WRITE(NUNIT,"(A)",advance='no') "+A,"
               ENDIF
            ELSE
               WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
            ENDIF
         ENDDO
         WRITE(NUNIT,"(A)",advance='no') ")"
         IF(LTERM) WRITE(NUNIT,*)
         RETURN
      END

! Calculate the one-electron part of the energy of a det
      REAL*8 FUNCTION CALCT(NI,NEL,G1,NBASIS)
         USE HElem
         USE OneEInts, only : GetTMatEl
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),G1(*),NBASIS,I
         LOGICAL ISCSF
         CALCT=0.D0
         IF(ISCSF(NI,NEL)) RETURN
         DO I=1,NEL
            CALCT=CALCT+DREAL(GetTMATEl(NI(I),NI(I)))
         ENDDO
         RETURN
      END


