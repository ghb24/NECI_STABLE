#include "macros.h"
MODULE Determinants
    use constants, only: dp
    use SystemData, only: BasisFN, tCSF, nel, G1, Brr, ECore, ALat, NMSH, &
                          nBasis, nBasisMax, tStoreAsExcitations, tHPHFInts, &
                          NIfToT, tCSF
    use IntegralsData, only: UMat, FCK, NMAX
    use csf, only: det_to_random_csf, iscsf, csf_orbital_mask, &
                   csf_yama_bit, CSFGetHelement
    use sltcnd_mod, only: sltcnd, sltcnd_excit, sltcnd_2, sltcnd_compat, &
                          sltcnd_knowIC
    use global_utilities
    use DetBitOps, only: EncodeBitDet
    use DeterminantData
    implicit none

    interface get_helement
        module procedure get_helement_compat
        module procedure get_helement_excit
        module procedure get_helement_normal
    end interface

    save
! Set by Calc on input
      INTEGER nActiveSpace(2)
        INTEGER, DIMENSION(:), POINTER :: SPECDET
        INTEGER :: tagSPECDET=0
        Logical TSPECDET

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

      INTEGER, allocatable :: DefDet(:)
      Logical :: tDefineDet
      integer :: tagDefDet=0

contains


    Subroutine DetPreFreezeInit()
        Use global_utilities
        use SystemData, only : nEl, ECore, Arr, Brr, G1, nBasis, LMS, nBasisMax,tFixLz, tUEGSpecifyMomentum
        integer ierr
        integer i,Lz
        type(BasisFn) s
   
        character(25), parameter :: this_routine='DetPreFreezeInit'
        Allocate(FDet(nEl), stat=ierr)
        LogAlloc(ierr, 'FDet', nEl, 4, tagFDet)
        IF(tDefineDet) THEN
            WRITE(6,*) 'Defining FDet according to input'
            do i=1,NEl
                FDet(i)=DefDet(i)
            enddo
        ELSE
             CALL GENFDET(BRR,G1,NBASIS,LMS,NEL,FDET)
             IF(tUEGSpecifyMomentum) THEN
                WRITE(6,*) 'Defining FDet according to a momentum input'
                CALL ModifyMomentum(FDET)
            ENDIF
        ENDIF
!      ENDIF
      WRITE(6,"(A)",advance='no') " Fermi det (D0):"
      call write_det (6, FDET, .true.)
      Call GetSym(FDet,nEl,G1,nBasisMax,s)
      WRITE(6,"(A)",advance='no') " Symmetry: "
      Call WriteSym(6,s%Sym,.true.)
      IF(tFixLz) THEN
         Call GetLz(FDet,nEl,Lz)
         WRITE(6,"(A,I5)") "Lz of Fermi det:",Lz
      ENDIF
      CALL NECI_ICOPY(NEL,FDET,1,NUHFDET,1)
      E0HFDET=ECORE
      DO I=1,NEL
         E0HFDET=E0HFDET+ARR(NUHFDET(i),2)
      ENDDO     
      WRITE(6,*) "Fock operator energy:",E0HFDET
    End Subroutine DetPreFreezeInit

    
    Subroutine DetInit()
        Use global_utilities
        use constants, only: dp
        use SystemData, only: nel, Alat, Boa, Coa, BOX, BRR, ECore
        use SystemData, only: G1, LMS, nBasis, STot, tCSFOLD, Arr,tHub,tUEG
        use SymData , only : nSymLabels,SymLabelList,SymLabelCounts,TwoCycleSymGens
        use IntegralsData, only: nfrozen
      
      real*8 DNDET
      integer i,ii,j
      integer*8 nDet
      integer ierr
      integer :: alpha,beta,symalpha,symbeta,endsymstate
      character(25), parameter :: this_routine='DetInit'
      LOGICAL :: tSuccess,tFoundOrbs(nBasis)
      integer :: ncsf

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
!         WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , DNDET
         NDET=-1
        ELSE
!         WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , NDET
        ENDIF
      
!C      CALL TC(I_HMAX,I_P,NWHTAY)

        
!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess=.true.
        tFoundOrbs(:)=.false.

        IF((.not.tHub).and.(.not.tUEG).and.TwoCycleSymGens) THEN
            do i=1,nSymLabels
!                WRITE(6,*) "NSymLabels: ",NSymLabels,i-1
                EndSymState=SymLabelCounts(1,i)+SymLabelCounts(2,i)-1
!                WRITE(6,*) "Number of states: ",SymLabelCounts(2,i)
                do j=SymLabelCounts(1,i),EndSymState

                    Beta=(2*SymLabelList(j))-1
                    Alpha=(2*SymLabelList(j))
                    SymAlpha=INT((G1(Alpha)%Sym%S),4)
                    SymBeta=INT((G1(Beta)%Sym%S),4)
!                    WRITE(6,*) "***",Alpha,Beta

                    IF(.not.tFoundOrbs(Beta)) THEN
                        tFoundOrbs(Beta)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF
                    IF(.not.tFoundOrbs(Alpha)) THEN
                        tFoundOrbs(Alpha)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF

                    IF(G1(Beta)%Ms.ne.-1) THEN
                        tSuccess=.false.
                    ELSEIF(G1(Alpha)%Ms.ne.1) THEN
                        tSuccess=.false.
                    ELSEIF((SymAlpha.ne.(i-1)).or.(SymBeta.ne.(i-1))) THEN
                        tSuccess=.false.
                    ENDIF
                enddo
            enddo
            do i=1,nBasis
                IF(.not.tFoundOrbs(i)) THEN
                    WRITE(6,*) "Orbital: ",i, " not found."
                    CALL Stop_All("SetupParameters","Orbital not found")
                ENDIF
            enddo
        ENDIF
        IF(.not.tSuccess) THEN
            WRITE(6,*) "************************************************"
            WRITE(6,*) "**                 WARNING!!!                 **"
            WRITE(6,*) "************************************************"
            WRITE(6,*) "Symmetry information not set up correctly in NECI initialisation"
            WRITE(6,*) "Will attempt to set up the symmetry again, but now in terms of spin orbitals"
            WRITE(6,*) "Old excitation generators will not work"
            WRITE(6,*) "I strongly suggest you check that the reference energy is correct."
            CALL SpinOrbSymSetup(.true.) 
        ELSE
            WRITE(6,*) "Symmetry and spin of orbitals correctly set up for excitation generators."
            WRITE(6,*) "Simply transferring this into a spin orbital representation."
            CALL SpinOrbSymSetup(.false.) 
        ENDIF
! From now on, the orbitals are also contained in symlabellist2 and symlabelcounts2.
! These are stored using spin orbitals.

        ! If we are using CSFs, then we need to convert this into a csf
        if (tCSF) then
            ncsf = det_to_random_csf (FDET)
            write (6, '("Generated starting CSF: ")', advance='no')
            call write_det (6, FDET, .true.)
        endif

    End Subroutine DetInit

    function get_helement_compat (nI, nJ, IC, iLutI, iLutJ) result (hel)
        use constants, only: n_int
       
        ! Get the matrix element of the hamiltonian. This assumes that we
        ! already know IC. We do not need to know iLutI, iLutJ (although
        ! they are helpful). This better fits the requirements of existing
        ! code than get_helement_normal.
        !
        ! In:  nI, nJ       - The determinants to consider
        !      iLutI, iLutJ - Bit representations of I,J (optional, helpful)
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The desired matrix element.

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in), optional :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: IC
        HElement_t :: hel

        character(*), parameter :: this_routine = 'get_helement_compat'

        if (tHPHFInts) &
            call stop_all (this_routine, "Should not be calling HPHF &
                          &integrals from here.")

        if (tCSF) then
            if (iscsf(nI) .or. iscsf(nJ)) then
                hel = CSFGetHelement (nI, nJ)
                return
            endif
        endif

        if (tStoreAsExcitations) &
            call stop_all(this_routine, "tStoreExcitations not supported")

        if (present(iLutJ)) then
            hel = sltcnd_knowIC (nI, nJ, iLutI, iLutJ, IC)
        else
            hel = sltcnd_compat (nI, nJ, IC)
        endif

        ! Add in ECore if for a diagonal element
        if (IC == 0) hel = hel + (ECore)
    end function
    
    function get_helement_normal (nI, nJ, iLutI, iLutJ, ICret) result(hel)
        use constants, only: n_int

        ! Get the matrix element of the hamiltonian.
        !
        ! In:  nI, nJ       - The determinants to consider
        !      iLutI, iLutJ - Bit representations of I,J (optional, helpful)
        ! Out: ICret        - The number of orbitals I,J differ by
        ! Ret: hel          - The desired matrix element.
        
        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in), optional :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(out), optional :: ICret
        HElement_t :: hel

        character(*), parameter :: this_routine = 'get_helement_normal'
        integer :: ex(2,2), IC
        integer(kind=n_int) :: ilut(0:NIfTot,2)

        if (tHPHFInts) &
            call stop_all (this_routine, "Should not be calling HPHF &
                          &integrals from here.")

        if (tCSF) then
            if (iscsf(nI) .or. iscsf(nJ)) then
                hel = CSFGetHelement (nI, nJ)
                return
            endif
        endif
         
        if (tStoreAsExcitations .and. nI(1) == -1 .and. nJ(1) == -1) then
            ! TODO: how to express requirement for double?
            !if (IC /= 2) &
            !    call stop_all (this_routine, "tStoreAsExcitations in &
            !                  &get_helement requires IC=2 (doubles)")

            ex(1,:) = nJ(4:5)
            ex(2,:) = nJ(6:7)
            hel = sltcnd_2 (ex, .false.)
        endif

        if (present(iLutJ)) then
            hel = sltcnd (nI, nJ, iLutI, iLutJ, IC)
        else
            call EncodeBitDet (nI, iLut(:,1))
            call EncodeBitdet (nJ, iLut(:,2))
            ! TODO: This is not an ideal place to end up...
            hel = sltcnd (nI, nJ, iLut(:,1), ilut(:,2), IC)
        endif

        ! Add in ECore for a diagonal element
        if (IC == 0) hel = hel + (ECore)

        ! If requested, return IC
        if (present(ICret)) then
            ICret = IC
        endif

    end function get_helement_normal

    function get_helement_excit (NI, NJ, IC, ExcitMat, TParity) result(hel)

        ! Calculate the Hamiltonian Matrix Element for a given determinant (or
        ! csf), when we have the excitation matrix and parity of the
        ! excitation.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        !      ex           - The excitation matrix
        !      tParity      - The parity of the excitation
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), IC
        integer, intent(in) :: ExcitMat(2,2)
        logical, intent(in) :: tParity
        HElement_t :: hel

        character(*), parameter :: this_routine = 'get_helement_excit'

        ! If we are using CSFs, then call the csf routine.
        ! TODO: Passing through of ExcitMat to CSFGetHelement
        if (tCSF) then
            if (iscsf(NI) .or. iscsf(NJ)) then
                hel = CSFGetHelement (nI, nJ)
                return
            endif
        endif
         
        if (IC < 0) &
            call stop_all(this_routine, "get_helement_excit should only be &
                         &used if we know the number of excitations and the &
                         &excitation matrix")

        hel = sltcnd_excit (nI, nJ, IC, ExcitMat, tParity)

        if (IC == 0)  hel = hel + (ECore)
    end function get_helement_excit

    function get_helement_det_only (nI, nJ, iLutI, iLutJ, ic, ex, tParity, &
                                    prob) result (hel)
        
        ! Calculate the Hamiltonian Matrix Element for a determinant as above.
        ! This function assumes that we have got it correct for determinants
        ! (i.e. no error checking), and no conditionals. It also has extra
        ! arguments for compatibility with the function pointer methods.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      iLutI, iLutJ - Bit representations (unused)
        !      ic           - The number of orbitals i,j differ by
        !      ex           - Excitation matrix
        !      tParity      - Parity of the excitation
        ! Ret: hel          - The H matrix element

        use constants, only: n_int
        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2,2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: prob
        HElement_t :: hel

        hel = sltcnd_excit (nI, nJ, IC, ex, tParity)

        if (IC == 0) hel = hel + ECore
    end function


      HElement_t function GetH0Element3(nI)
         ! Wrapper for GetH0Element.
         ! Returns the matrix element of the unperturbed Hamiltonian, which is
         ! just the sum of the eigenvalues of the occupied orbitals and the core
         ! energy.
         !  Note that GetH0Element{1,2} don't exist. The name is to be
         !  consistent with GetHElement3, i.e. offer the most abstraction possible.
         ! In: 
         !    nI(nEl)  list of occupied spin orbitals in the determinant.
         use constants, only: dp
         use SystemData, only: nEl,nBasis,Arr,ECore
         integer nI(nEl)
         HElement_t hEl
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
         use constants, only: dp
         implicit none
         integer nI(nEl),nEl,nBasis
         HElement_t hEl
         real*8 Arr(nBasis,2),ECore
         integer i
         if(tStoreAsExcitations.and.nI(1).eq.-1) then
!The excitation storage starts with -1.  The next number is the excitation level,L .  
!Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
            hEl=0.d0
            do i=4,nI(2)+4-1
               hEl=hEl-(Arr(nI(i),2))
            enddo
            do i=i,i+nI(2)-1
               hEl=hEl+(Arr(nI(i),2))
            enddo
         else
            hEl=ECore
            do i=1,nEl
               hEl=hEl+(Arr(nI(i),2))
            enddo
         endif
!         call writedet(77,nI,nel,.false.)
!         write(77,*) "H0",hEl
!         call flush(77)
      end subroutine

      subroutine DetFreezeBasis(GG)
        Use Determinants, only: FDet, nUHFDet, write_det, write_det_len
        use SystemData, only : nEl, nBasis, nBasisMax,BasisFN,G1,tFixLz
        use IntegralsData, only : nFrozen,nFrozenIn
        implicit none
        integer i,j
        INTEGER GG(*),Lz
        Type(BasisFn) s
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
            IF(J.NE.NEL-NFROZEN-NFROZENIN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               call write_det (6, FDET, .true.)
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
            IF(J.NE.NEL-NFROZEN-NFROZENIN) THEN
               WRITE(6,*) "Failed Freezing Det:"
               call write_det (6, nUHFDET, .true.)
               STOP "After Freezing, UHFDET has wrong number of electrons"
            ENDIF
         ENDIF
         WRITE(6,"(A)",advance='no') " Post-Freeze Fermi det (D0):"
         call write_det_len (6, fDet, nel-nfrozen-nfrozenin, .true.)
         WRITE(6,"(A)",advance='no') " Symmetry: "
         Call GetSym(FDet,nEl-nFrozen-nFrozenIn,G1,nBasisMax,s)
         Call WriteSym(6,s%Sym,.true.)
         IF(tFixLz) THEN
             Call GetLz(FDet,nEl-nFrozen-nFrozenIn,Lz)
             WRITE(6,"(A,I5)") " Lz of Fermi det:",Lz
         ENDIF

      end subroutine


      LOGICAL FUNCTION ISUHFDET(NI,NEL)
         use SystemData , only : TUSEBRILLOUIN
         Use Determinants, only : NUHFDET
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         INTEGER I
         IF(.NOT.TUSEBRILLOUIN) THEN
             ISUHFDET=.FALSE.
             RETURN
         ENDIF
         ISUHFDET=.TRUE.
         DO I=NEL,1,-1
            IF(NI(I).NE.NUHFDET(I)) THEN
                ISUHFDET=.FALSE.
                EXIT
            ENDIF
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
      !SUBROUTINE WRITEDET(NUNIT,NI,NEL,LTERM)
      !   use legacy_data, only: CSF_NBSTART
      !   IMPLICIT NONE
      !   INTEGER NUNIT,NEL,NI(NEL),I
      !   LOGICAL LTERM
      !   INTEGER IEL
      !   CHARACTER*2 SUFF
      !   WRITE(NUNIT,"(A)",advance='no') "("
      !   DO I=1,NEL
      !      IEL=NI(I)
      !      IF(IEL.GE.CSF_NBSTART) THEN
      !         WRITE(NUNIT,"(I3)",advance='no'),(IEL-CSF_NBSTART)/4+1
      !         IEL=IAND(IEL-CSF_NBSTART,3)
      !         IF(IEL.EQ.0) THEN
      !            WRITE(NUNIT,"(A)",advance='no') "-B,"
      !         ELSEIF(IEL.EQ.1) THEN
      !            WRITE(NUNIT,"(A)",advance='no') "-A,"
      !         ELSEIF(IEL.EQ.2) THEN
      !            WRITE(NUNIT,"(A)",advance='no') "+B,"
      !         ELSE
      !            WRITE(NUNIT,"(A)",advance='no') "+A,"
      !         ENDIF
      !      ELSE
      !         WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
      !      ENDIF
      !   ENDDO
      !   WRITE(NUNIT,"(A)",advance='no') ")"
      !   IF(LTERM) WRITE(NUNIT,*)
      !   RETURN
      !END
    subroutine writedet_oldcsf (nunit, nI, nel, lTerm)
        use systemdata, only: tCSF, tCSFOLD
        
        ! Write a human readable determinant to specified file unit. For use
        ! with old csf routines.
        ! Not easy to test, as both iscsf routines will return true sometimes.
        ! This is here for use if it becomes necessary (eg debugging)
        !
        ! In: nunit    - File unit 
        !     nI (nel) - Determinant to print
        !     nel      - Number of electrons
        !     lTerm    - Do we write an end-of-line character
        
        use legacy_data, only: CSF_NBSTART
        implicit none
        integer, intent(in) :: nunit, nel, nI(nel)
        logical, intent(in) :: lTerm
        integer :: i, orb
        logical iscsf, bCSF

        ! Is this a csf? Note use of old (non-modularised) iscsf
        bCSF = (tCSF .or. tCSFOLD) .and. iscsf(nI, nel)

        ! Begin with an open bracket
        write(nunit,'("(")',advance='no')
        do i=1,nel
            orb = nI(i)
            if (bCSF .and. orb > CSF_NBSTART) then
                ! Output the component orbital
                orb = (orb - CSF_NBSTART) / 4 + 1
                write(nunit,'(i4)',advance='no') orb

                ! Output components of Yamanouchi symbol
                orb = nI(i) - csf_nbstart
                if (btest(orb, 1)) then
                    write(nunit,'("+")',advance='no')
                else
                    write(nunit,'("-")',advance='no')
                endif

                ! Output components of Ms
                if (btest(orb, 0)) then
                    write(nunit,'("A")',advance='no')
                else
                    write(nunit,'("B")',advance='no')
                endif
            else
                write(nunit,'(i6)',advance='no') orb
            endif
            if (i /= nel) write(nunit,'(",")',advance='no')
        enddo

        ! Close the written determinant off
        write(nunit,'(")")',advance='no')
        if (lTerm) write(nunit,*)

    end subroutine


! Calculate the one-electron part of the energy of a det
      REAL*8 FUNCTION CALCT(NI,NEL,G1,NBASIS)
         use constants, only: dp
         USE SystemData, only : BasisFN
         USE OneEInts, only : GetTMatEl
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,I
         TYPE(BasisFN) :: G1(*)
         LOGICAL ISCSF
         CALCT=0.D0
         IF(ISCSF(NI,NEL)) RETURN
         DO I=1,NEL
            CALCT=CALCT+GetTMATEl(NI(I),NI(I))
         ENDDO
         RETURN
      END

! Write bit-determinant NI to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
      SUBROUTINE WriteBitDet(nUnit,iLutnI,lTerm)
         use SystemData, only : nEl, nIfTot
         use DetBitops, only: DecodeBitDet
         use Determinants, only: write_det
         use constants, only: n_int
         implicit none
         integer nUnit,nI(nEl)
         integer(kind=n_int) :: iLutnI(0:nIfTot)
         logical lTerm
         CALL DecodeBitDet(nI,iLutnI)
         call write_det (nUnit, nI, lTerm)
      END

! Write bit-determinant NI to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
      SUBROUTINE WriteBitEx(nUnit,iLutRef,iLutnI,lTerm)
         use SystemData, only : nEl, NIfTot
         use constants, only: n_int
         implicit none
         integer nUnit,nExpI(nEl)
         integer(kind=n_int) :: iLutRef(0:nIfTot),iLutnI(0:nIfTot)
         integer Ex(2,nEl)
         logical lTerm
         logical tSign
         INTEGER iEl,I
         EX(1,1)=nEl  !Indicate the length of EX
         CALL GetBitExcitation(iLutRef,iLutnI,Ex,tSign)
         WRITE(NUNIT,"(A)",advance='no') "("
! First the excit from
         DO I=1,NEL
            IEL=EX(1,I)
            if(iEl.eq.0) EXIT
            WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
         ENDDO
         IF(tSign) THEN
            WRITE(NUNIT,"(A)",advance='no') ")->-("
         ELSE
            WRITE(NUNIT,"(A)",advance='no') ")->+("
         ENDIF
!Now the excit to
         DO I=1,NEL
            IEL=EX(2,I)
            if(iEl.eq.0) EXIT
            WRITE(NUNIT,"(I5,A)",advance='no') IEL,","
         ENDDO
         WRITE(NUNIT,"(A)",advance='no') ")"
         IF(LTERM) WRITE(NUNIT,*)
         RETURN
      END

! This takes a the ground state FDet generated for the UEG and changes its total momentum
! According to input options
      subroutine ModifyMomentum(FDet)
        use SystemData, only : nEl,G1,k_momentum,nBasis,tUEG
        implicit none
        integer :: i,j ! Loop variables
        integer, intent(inout) :: FDet(NEl)
        integer :: k_total(3) ! Stores the total momentum of FDet
        integer :: delta_k(3) ! Stores the difference between the current FDet and the FDet we're aiming for
        integer, allocatable :: kPointToBasisFn(:,:,:,:) ! Look up table for kPoints to basis functions
        integer :: kmaxX,kminX,kmaxY,kminY,kmaxZ,kminZ,iSpinIndex ! Stores the limits of kPointToBasisFn
        integer :: det_sorted(NEl), e_store ! Storage for the sorting routine
        logical :: sorted ! As above
        integer :: wrapped_index
        integer :: k_old, k_new

        IF(.not.tUEG) call stop_all("ModifyMomentum", "Only works for UEG")

        ! Finds current momentum, and finds the difference between this and the target momentum
        ! most commonly will be zero to start with
        k_total(1)=0
        k_total(2)=0
        k_total(3)=0
        do j=1,NEl
            k_total(1)=k_total(1)+G1(FDet(j))%k(1)
            k_total(2)=k_total(2)+G1(FDet(j))%k(2)
            k_total(3)=k_total(3)+G1(FDet(j))%k(3)
        enddo
        delta_k=k_momentum-k_total

        if (delta_k(1).eq.0.and.delta_k(2).eq.0.and.delta_k(3).eq.0) write(6,*) "WARNING: specified momentum is ground state"

        ! Creates a look-up table for k-points (this was taken from symrandexcit2.F90)
        kmaxX=0
        kminX=0
        kmaxY=0
        kminY=0
        kminZ=0
        kmaxZ=0
        do i=1,nBasis 
            IF(G1(i)%k(1).gt.kmaxX) kmaxX=G1(i)%k(1)
            IF(G1(i)%k(1).lt.kminX) kminX=G1(i)%k(1)
            IF(G1(i)%k(2).gt.kmaxY) kmaxY=G1(i)%k(2)
            IF(G1(i)%k(2).lt.kminY) kminY=G1(i)%k(2)
            IF(G1(i)%k(3).gt.kmaxZ) kmaxZ=G1(i)%k(3)
            IF(G1(i)%k(3).lt.kminZ) kminZ=G1(i)%k(3)
        enddo
        allocate(kPointToBasisFn(kminX:kmaxX,kminY:kmaxY,kminZ:kmaxZ,2))
        do i=1,nBasis
            iSpinIndex=(G1(i)%Ms+1)/2+1 ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
            kPointToBasisFn(G1(i)%k(1),G1(i)%k(2),G1(i)%k(3),iSpinIndex)=i
        enddo

        ! For each of the three dimensions, nudge electrons one at a time by one momentum unit until delta_k is reached
        det_sorted=FDet

        ! Bubble sort to order det_sorted in order of kx of the corresponding electron
        do 
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (G1(det_sorted(j))%k(1).gt.G1(det_sorted(i))%k(1)) then
                    sorted=.false.
                    e_store=det_sorted(i)
                    det_sorted(i)=det_sorted(j)
                    det_sorted(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        ! Nudge momenta one at a time
        if (delta_k(1).gt.0) then
            do i=1,delta_k(1)
                wrapped_index=mod(i,NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index.eq.0) then    ! Deal with the i=NEl case
                    wrapped_index=NEl
                endif
                j=wrapped_index                 ! For convenience asign this to j
                k_new=G1(det_sorted(j))%k(1)+1  ! Find the new momentum of this electron
                if (k_new.gt.kmaxX) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1     ! Spin of the new orbital is the same as the old
                det_sorted(j)=kPointToBasisFn(k_new,G1(det_sorted(j))%k(2),G1(det_sorted(j))%k(3),iSpinIndex) ! Finds basis number for the new momentum
            enddo
        else if (delta_k(1).lt.0) then ! For the negative case, i must run through negative numbers
            do i=-1,delta_k(1),-1
                wrapped_index=mod(i,NEl)
                if (wrapped_index.eq.0) then
                    wrapped_index=-NEl
                endif
                j=NEl+wrapped_index+1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new=G1(det_sorted(j))%k(1)-1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new.lt.kminX) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1 ! Spin of the new orbital is the same as the old
                det_sorted(j)=kPointToBasisFn(k_new,G1(det_sorted(j))%k(2),G1(det_sorted(j))%k(3),iSpinIndex) ! Finds basis number for the new momentum
            enddo
        endif

        FDet=det_sorted
        
        !====ky treated as kx above
        do 
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (G1(det_sorted(j))%k(2).gt.G1(det_sorted(i))%k(2)) then
                    sorted=.false.
                    e_store=det_sorted(i)
                    det_sorted(i)=det_sorted(j)
                    det_sorted(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        ! Nudge momenta one at a time
        if (delta_k(2).gt.0) then
            do i=1,delta_k(2)
                wrapped_index=mod(i,NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index.eq.0) then    ! Deal with the i=NEl case
                    wrapped_index=NEl
                endif
                j=wrapped_index                 ! For convenience asign this to j
                k_new=G1(det_sorted(j))%k(2)+1  ! Find the new momentum of this electron
                if (k_new.gt.kmaxY) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1     ! Spin of the new orbital is the same as the old
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),k_new,G1(det_sorted(j))%k(3),iSpinIndex) ! Finds basis number for the new momentum
            enddo
        else if (delta_k(2).lt.0) then ! For the negative case, i must run through negative numbers
            do i=-1,delta_k(2),-1
                wrapped_index=mod(i,NEl)
                if (wrapped_index.eq.0) then
                    wrapped_index=-NEl
                endif
                j=NEl+wrapped_index+1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new=G1(det_sorted(j))%k(2)-1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new.lt.kminY) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1 ! Spin of the new orbital is the same as the old
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),k_new,G1(det_sorted(j))%k(3),iSpinIndex) ! Finds basis number for the new momentum
            enddo
        endif

        FDet=det_sorted
        
        !====kz treated as kx and ky above
        do 
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (G1(det_sorted(j))%k(3).gt.G1(det_sorted(i))%k(3)) then
                    sorted=.false.
                    e_store=det_sorted(i)
                    det_sorted(i)=det_sorted(j)
                    det_sorted(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        ! Nudge momenta one at a time
        if (delta_k(3).gt.0) then
            do i=1,delta_k(3)
                wrapped_index=mod(i,NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index.eq.0) then    ! Deal with the i=NEl case
                    wrapped_index=NEl
                endif
                j=wrapped_index                 ! For convenience asign this to j
                k_new=G1(det_sorted(j))%k(3)+1  ! Find the new momentum of this electron
                if (k_new.gt.kmaxZ) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1     ! Spin of the new orbital is the same as the old
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),G1(det_sorted(j))%k(2),k_new,iSpinIndex) ! Finds basis number for the new momentum
            enddo
        else if (delta_k(3).lt.0) then ! For the negative case, i must run through negative numbers
            do i=-1,delta_k(3),-1
                wrapped_index=mod(i,NEl)
                if (wrapped_index.eq.0) then
                    wrapped_index=-NEl
                endif
                j=NEl+wrapped_index+1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new=G1(det_sorted(j))%k(3)-1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new.lt.kminZ) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                endif
                iSpinIndex=(G1(j)%Ms+1)/2+1 ! Spin of the new orbital is the same as the old
                det_sorted(j)=kPointToBasisFn(G1(det_sorted(j))%k(1),G1(det_sorted(j))%k(2),k_new,iSpinIndex) ! Finds basis number for the new momentum
            enddo
        endif

        FDet=det_sorted
        
        ! Bubble sort to order FDet back into increasing order by number
        do 
            sorted=.true.
            do i=1,NEl-1
                j=i+1
                if (FDet(j).lt.FDet(i)) then
                    sorted=.false.
                    e_store=FDet(i)
                    FDet(i)=FDet(j)
                    FDet(j)=e_store
                endif
            enddo
            if (sorted) exit
        enddo

        write(6,*) "Total momentum set to", k_momentum

      end subroutine ModifyMomentum
