! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module OneEInts

! For calculating, storing and retrieving the one electron integrals, <i|h|j>, where 
! h is the one-particle Hamiltonian operator and i and j are spin-orbitals.
! Variables associated with this usually have TMAT in the name.

! The 2 or NEW versions of routines and variables are used for the set of
! orbitals after freezing has been done. The "original" versions are used in the
! pre-freezing stage.

use constants, only: dp
use SystemData, only: TSTARSTORE
use MemoryManager, only: TagIntType
use util_mod, only: get_free_unit

implicit none

save
public

! For storing the <i|h|j> elements which are non-zero by symmetry.
! Ranges from -1 to N, where N is the # of non-zero elements. The -1
! element is used to return the zero value.  Used for Abelian symmetries,
! where i and j must span the same representation in order for <i|h|j> 
! to be non-zero.
! The symmetries used in Dalton and Molpro are all Abelian.
HElement_t, dimension(:), POINTER :: TMATSYM
! For non-Abelian symmetries, we store the entire <i|h|j> matrix, as 
! the direct products of the representations can contain the totally
! symmetric representation (and are not necessarily the same).  We could
! compress this in a similar fashion at some point.
HElement_t, dimension(:,:), POINTER :: TMAT2D

HElement_t, dimension(:), POINTER :: TMATSYM2
HElement_t, dimension(:,:), POINTER :: TMAT2D2

! True if using TMatSym in CPMD (currently only if using k-points, which form
! an Abelian group).
logical tCPMDSymTMat
logical tOneElecDiag    !Indicates that the one-electron integral matrix is diagonal - 
                        !basis functions are eigenstates of KE operator.

! Memory book-keeping tags
integer(TagIntType) :: tagTMat2D=0
integer(TagIntType) :: tagTMat2D2=0
integer(TagIntType) :: tagTMATSYM=0,tagTMATSYM2=0

contains

      pure INTEGER FUNCTION TMatInd(I,J)
        ! In: 
        !    i,j: spin orbitals.
        ! Return the index of the <i|h|j> element in TMatSym2.
        ! This is only used with TSTARSTORE, where the TMAT is compressed to only store states, not spin-orbitals,
        ! Added compression supplied by only storing symmetry allowed integrals - therefore needs SymData info.
        ! We assume a restricted calculation.  We note that TMat is a Hermitian matrix.
        ! For the TMat(i,j) to be non-zero, i and j have to belong to the same symmetry, as we're working in Abelian groups.
        ! We store the non-zero elements in TMatSym(:).
        !
        ! The indexing scheme is:
        !   Arrange basis into blocks of the same symmetry, ie so TMat is of a block
        !   (but not necessarily block diagonal) nature.  Order the blocks by their
        !   symmetry index.
        !   Within each block, convert i and j to a spatial index and label by the 
        !   (energy) order in which they come within that symmetry: i->k,j->l.
        !   We only need to store the upper diagonal of each block:
        !        blockind=k*(k-1)/2+l, l<k.
        !   The overall index depends on the number of non-zero integrals in the blocks
        !   preceeding the block i and j are in.  This is given by SYMLABELINTSCUM(symI-1).
        !        TMatInd=k*(k-1)/2+l + SymLabelIntsCum(symI-1).
        !   If the element is zero by symmetry, return -1 (TMatSym(-1) is set to 0). 
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymClasses,StateSymMap,SymLabelIntsCum
        IMPLICIT NONE
        integer, intent(in) :: i, j
        integer A,B,symI,symJ,Block,ind,K,L
        A=mod(I,2)
        B=mod(J,2)
        !If TMatInd = -1, then the spin-orbitals have different spins, or are symmetry disallowed 
!therefore have a zero integral (apart from in UHF - might cause problem if we want this)
        IF(A.ne.B) THEN
            TMatInd=-1
            RETURN
        ENDIF
        !To convert to spatial orbitals
        K=(I+1)/2
        L=(J+1)/2
        
        symI=SYMCLASSES(K)
        symJ=SYMCLASSES(L)

        IF(symI.ne.symJ) THEN
            ! <i|t|j> is symmetry forbidden.
            TMatInd=-1
            RETURN
        ELSE
            ! Convert K and L into their index within their symmetry class.
            K=StateSymMap(K)
            L=StateSymMap(L)
            ! Block index: how many symmetry allowed integrals are there in the
            ! preceeding blocks?
            IF(symI.eq.1) THEN
                Block=0
            ELSE
                Block=SYMLABELINTSCUM(symI-1)
            ENDIF
            ! Index within symmetry block.
            IF(K.ge.L) THEN
                ind=(K*(K-1))/2+L
            ELSE
                ind=(L*(L-1))/2+K
            ENDIF
            TMatInd=Block+ind
            RETURN
        ENDIF
      END FUNCTION TMatInd
     


      INTEGER FUNCTION NEWTMatInd(I,J)
      ! In: 
      !    i,j: spin orbitals.
      ! Return the index of the <i|h|j> element in TMatSym2.
      ! See notes for TMatInd. Used post-freezing.
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymClasses2,StateSymMap,SymLabelIntsCum2
        IMPLICIT NONE
        INTEGER I,J,A,B,symI,symJ,Block,ind,K,L
        A=mod(I,2)
        B=mod(J,2)
        ! If TMatInd = -1, then the spin-orbitals have different spins, or are symmetry 
!disallowed therefore have a zero integral (apart from in UHF - might cause problem if we want this)
        IF(A.ne.B) THEN
            NEWTMatInd=-1
            RETURN
        ENDIF
        !To convert to spatial orbitals
        K=(I+1)/2
        L=(J+1)/2
        
        symI=SYMCLASSES2(K)
        symJ=SYMCLASSES2(L)

        IF(symI.ne.symJ) THEN
            NEWTMatInd=-1
            RETURN
        ELSE
            IF(symI.eq.1) THEN
                Block=0
            ELSE
                Block=SYMLABELINTSCUM2(symI-1)
            ENDIF
            ! Convert K and L into their index within their symmetry class.
            K=StateSymMap(K)
            L=StateSymMap(L)
            IF(K.ge.L) THEN
                ind=(K*(K-1))/2+L
            ELSE
                ind=(L*(L-1))/2+K
            ENDIF
            NEWTMatInd=Block+ind
            RETURN
        ENDIF
      END FUNCTION NewTMatInd

    
    elemental function GetTMatEl (i, j) result(ret)

        ! Return the one electron integral <i|h|j>
        !
        ! In: i,j - Spin orbitals

        integer, intent(in) :: i, j
        HElement_t :: ret
#ifdef __CMPLX
        HElement_t :: t
#endif

        if (tStarStore) then
            ret = TMatSym(TMatInd(i, j))
        else if (tCPMDSymTMat) then
            ! TMat is Hermitian, rather than symmetric.
            ! Only the upper diagonal of each symmetry block is stored
            if (j .ge. i) then
                ret = TMatSym(TMatInd(i,j))
            else
                ! Work around a bug in gfortran's parser: it doesn't like
                ! doing conjg(TMatSym).
#ifdef __CMPLX
                t = TMatSym(TmatInd(i,j))
                ret = conjg(t)
#else
                ret = TMatSym(TmatInd(i,j))
#endif
            endif
        else
            if(tOneElecDiag) then
                if(i.ne.j) then
                    ret = 0.0_dp
                else
                    ret = TMat2D(i,1)
                endif
            else
                ret = TMat2D(i, j)
            endif
        endif
    end function GetTMatEl

      FUNCTION GetNewTMatEl(I,J)
      ! In: 
      !    i,j: spin orbitals.
      ! Return <i|h|j>, the "TMat" element.
      ! Used post-freezing. See also GetTMatEl.
        IMPLICIT NONE
        INTEGER I,J
        HElement_t GetNEWTMATEl

        IF(TSTARSTORE) THEN
            GetNEWTMATEl=TMATSYM2(NEWTMATInd(I,J))
        else if (tCPMDSymTMat) then
#ifdef __CMPLX
            if (j.ge.i) then
                GetNewTMatEl=TMATSYM2(TMatInd(I,J))
            else
                GetNewTMatEl=Conjg(TMATSYM2(TMatInd(I,J)))
            end if
#else
            GetNewTMatEl=TMATSYM2(TMatInd(I,J))
#endif
        ELSE
            if(tOneElecDiag) then
                if(I.ne.J) then
                    GetNEWTMATEl=0.0_dp
                else
                    GetNewTMATEl=TMAT2D2(I,1)
                endif
            else
                GetNEWTMATEl=TMAT2D2(I,J)
            endif
        ENDIF
      END FUNCTION GetNewTMatEl


      SUBROUTINE WriteTMat(NBASIS)
        ! In:
        !    nBasis: size of basis (# orbitals).
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymLabelCounts,SymLabelCountsCum,nSymLabels
        use SymData, only: SymLabelIntsCum,SymLabelIntsCum2,SymLabelCountsCum2
        IMPLICIT NONE
        INTEGER II,I,J,NBASIS,iunit
        
        iunit = get_free_unit()
        open(iunit, file="TMATSYMLABEL", status="unknown")
        IF(associated(SYMLABELINTSCUM)) THEN
            write(iunit,*) "SYMLABELCOUNTS,SYMLABELCOUNTSCUM,SYMLABELINTSCUM:"
            DO I=1,NSYMLABELS
                WRITE(iunit,"(I5)",advance='no') SYMLABELCOUNTS(2,I)
                CALL neci_flush(iunit)
            ENDDO
            WRITE(iunit,*) ""
            DO I=1,NSYMLABELS
                WRITE(iunit,"(I5)",advance='no') SYMLABELCOUNTSCUM(I)
                CALL neci_flush(iunit)
            ENDDO
            WRITE(iunit,*) ""
            DO I=1,NSYMLABELS
                WRITE(iunit,"(I5)",advance='no') SYMLABELINTSCUM(I)
                CALL neci_flush(iunit)
            ENDDO
            WRITE(iunit,*) ""
            WRITE(iunit,*) "**********************************"
        ENDIF
        IF(associated(SYMLABELINTSCUM2)) THEN
            write(iunit,*) "SYMLABELCOUNTS,SYMLABELCOUNTSCUM2,SYMLABELINTSCUM2:"
            DO I=1,NSYMLABELS
                WRITE(iunit,"(I5)",advance='no') SYMLABELCOUNTS(2,I)
                CALL neci_flush(iunit)
            ENDDO
            WRITE(iunit,*) ""
            DO I=1,NSYMLABELS
                WRITE(iunit,"(I5)",advance='no') SYMLABELCOUNTSCUM2(I)
                CALL neci_flush(iunit)
            ENDDO
            WRITE(iunit,*) ""
            DO I=1,NSYMLABELS
                WRITE(iunit,"(I5)",advance='no') SYMLABELINTSCUM2(I)
                CALL neci_flush(iunit)
            ENDDO
            WRITE(iunit,*) ""
            WRITE(iunit,*) "**********************************"
            CALL neci_flush(iunit)
        ENDIF
        WRITE(iunit,*) "TMAT:"
        IF(TSTARSTORE) THEN
            DO II=1,NSYMLABELS
                DO I=SYMLABELCOUNTSCUM(II-1)+1,SYMLABELCOUNTSCUM(II)
                    DO J=SYMLABELCOUNTSCUM(II-1)+1,I
                        WRITE(iunit,*) I,J,GetTMATEl((2*I),(2*J))
                        CALL neci_flush(iunit)
                    ENDDO
                ENDDO
            ENDDO
        ELSE
            DO I=1,NBASIS,2
                DO J=1,NBASIS,2
                    WRITE(iunit,*) (I+1)/2,(J+1)/2, GetTMATEl(I,J)
                ENDDO
            ENDDO
        ENDIF
        WRITE(iunit,*) "**********************************"
        CALL neci_flush(iunit)
        IF(ASSOCIated(TMATSYM2).or.ASSOCIated(TMAT2D2)) THEN
            WRITE(iunit,*) "TMAT2:"
            DO II=1,NSYMLABELS
                DO I=SYMLABELCOUNTSCUM(II-1)+1,SYMLABELCOUNTSCUM(II)
                    DO J=SYMLABELCOUNTSCUM(II-1)+1,I
                        WRITE(iunit,*) I,J,GetNEWTMATEl((2*I),(2*J))
                        CALL neci_flush(iunit)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF
        WRITE(iunit,*) "*********************************"
        WRITE(iunit,*) "*********************************"
        CALL neci_flush(iunit)
      END SUBROUTINE WriteTMat
        
!Routine to calculate number of elements allocated for TMAT matrix
      SUBROUTINE CalcTMATSize(nBasis,iSize)
      use SymData, only: SymLabelCounts,nSymLabels
      INTEGER :: iSize,nBasis,basirrep,i,Nirrep

          IF(TSTARSTORE.or.tCPMDSymTMat) THEN 
              iSize=0
              Nirrep=NSYMLABELS
              do i=1,Nirrep
                  basirrep=SYMLABELCOUNTS(2,i)
                  ! Block diagonal.
                  iSize=iSize+(basirrep*(basirrep+1))/2
              enddo
              iSize=iSize+2     !lower index is -1
          ELSE
              if(tOneElecDiag) then
                  iSize=nBasis
              else
                  iSize=nBasis*nBasis
              endif
          ENDIF

      END SUBROUTINE CalcTMATSize

      SUBROUTINE SetupTMAT(nBASIS,iSS,iSize)   
        ! In:
        !    nBasis: number of basis functions (orbitals).
        !    iSS: ratio of nBasis to the number of spatial orbitals.
        !    iSize: number of elements in TMat/TMatSym.
        ! Initial allocation of TMat2D or TMatSym (if using symmetry-compressed
        ! storage of the <i|h|j> integrals).
        use CPMDData, only: tKP
        use SystemData, only: tCPMD,tVASP
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymLabelCounts,SymLabelCountsCum,SymClasses
        use SymData, only: SymLabelIntsCum,nSymLabels,StateSymMap
        use SymData, only: tagSymLabelIntsCum,tagStateSymMap,tagSymLabelCountsCum
        use HElem, only: HElement_t_size
        use global_utilities
        IMPLICIT NONE
        integer Nirrep,nBasis,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer iSize, iunit
        character(len=*),parameter :: thisroutine='SetupTMAT'
        
        ! If this is a CPMD k-point calculation, then we're operating
        ! under Abelian symmetry: can use George's memory efficient
        ! TMAT.  
        if (tCPMD) tCPMDSymTMat=tKP
        if (tVASP) tCPMDSymTMat=.true.
        IF(TSTARSTORE.or.tCPMDSymTMat) THEN 
            ! Set up info for indexing scheme (see TMatInd for full description).
            Nirrep=NSYMLABELS
            nBi=nBasis/iSS
            iSize=0

            if (associated(SymLabelIntsCum)) then
                deallocate(SymLabelIntsCum)
                call LogMemDealloc(thisroutine,tagSymLabelIntsCum)
            end if
            if (allocated(StateSymMap)) then
                deallocate(StateSymMap)
                call LogMemDealloc(thisroutine,tagStateSymMap)
            end if
            if (associated(SymLabelCountsCum)) then
                deallocate(SymLabelCountsCum)
                call LogMemDealloc(thisroutine,tagSymLabelCountsCum)
            end if

            allocate(SymLabelIntsCum(nIrrep))
            call LogMemAlloc('SymLabelIntsCum',nIrrep,4,thisroutine,tagSymLabelIntsCum)
            allocate(SymLabelCountsCum(nIrrep))
            call LogMemAlloc('SymLabelCountsCum',nIrrep,4,thisroutine,tagSymLabelCountsCum)
            allocate(StateSymMap(nBi))
            call LogMemAlloc('StateSymMap',nBi,4,thisroutine,tagStateSymMap)
            SYMLABELCOUNTSCUM(1:Nirrep)=0
            SYMLABELINTSCUM(1:Nirrep)=0
            do i=1,Nirrep
                basirrep=SYMLABELCOUNTS(2,i)
                ! Block diagonal.
                iSize=iSize+(basirrep*(basirrep+1))/2
                ! # of integrals in this symmetry class and all preceeding
                ! symmetry classes.
                SYMLABELINTSCUM(i)=iSize
                ! JSS: we no longer need SymLabelCountsCum, but it's a useful test.
                ! SymLabelCountsCum is the cumulative total # of basis functions
                ! in the preceeding symmetry classes.
                IF(i.eq.1) THEN
                    SYMLABELCOUNTSCUM(i)=0
                ELSE
                    DO t=1,(i-1)
                        SYMLABELCOUNTSCUM(i)=SYMLABELCOUNTSCUM(i)+SYMLABELCOUNTS(2,t)
                    ENDDO
                ENDIF
                ! JSS: Label states of symmetry i by the order in which they come.
                nStateIrrep=0
                do iState=1,nBi
                    if (SymClasses(iState).eq.i) then
                        nStateIrrep=nStateIrrep+1
                        StateSymMap(iState)=nStateIrrep
                    end if
                end do
            enddo
            ! The number of basis functions before the last irrep
            ! plus the number of basis functions in the last irrep
            ! (which is in basirrep on exiting the above do loop)
            ! must equal the total # of basis functions.
            IF((SYMLABELCOUNTSCUM(Nirrep)+basirrep).ne.nBI) THEN
                iunit = get_free_unit()
                open(iunit, file="TMATSYMLABEL", status="unknown")
                DO i=1,Nirrep
                    WRITE(iunit,*) SYMLABELCOUNTSCUM(i)
                ENDDO
                write(iunit,*) "***************"
                write(iunit,*) NBI
                CALL neci_flush(iunit)
                close(iunit)
                STOP 'Not all basis functions found while setting up TMAT'
            ENDIF
            !iSize=iSize+2
            !This is to allow the index of '-1' in the array to give a zero value
            !Refer to TMatSym(-1) for when <i|h|j> is zero by symmetry.
            
            Allocate(TMATSYM(-1:iSize),STAT=ierr)
            Call LogMemAlloc('TMATSym',iSize+2,HElement_t_size*8,thisroutine,tagTMATSYM,ierr)
            TMATSYM=(0.0_dp)

        ELSE

            if(tOneElecDiag) then
                ! In the UEG, the orbitals are eigenfunctions of KE operator, so TMAT is diagonal. 
                iSize=nBasis
                Allocate(TMAT2D(nBasis,1),STAT=ierr)
                call LogMemAlloc('TMAT2D',nBasis,HElement_t_size*8,thisroutine,tagTMat2D)
                TMAT2D=(0.0_dp)
            else
                ! Using a square array to hold <i|h|j> (incl. elements which are
                ! zero by symmetry).
                iSize=nBasis*nBasis
                Allocate(TMAT2D(nBasis,nBasis),STAT=ierr)
                call LogMemAlloc('TMAT2D',nBasis*nBasis,HElement_t_size*8,thisroutine,tagTMat2D)
                TMAT2D=(0.0_dp)
            endif
        
        ENDIF
    
      END SUBROUTINE SetupTMAT


     
      SUBROUTINE SetupTMAT2(nBASISFRZ,iSS,iSize)
        ! In:
        !    nBasisFRZ: number of active basis functions (orbitals).
        !    iSS: ratio of nBasisFRZ to the number of active spatial orbitals.
        !    iSize: number of elements in TMat2/TMatSym2.
        ! Initial allocation of TMat2D2 or TMatSym2 (if using symmetry-compressed
        ! storage of the <i|h|j> integrals) for post-freezing.
        ! See also notes in SetupTMat.
        use CPMDData, only: tKP
        use SystemData, only: tCPMD,tVASP
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymLabelCounts,SymClasses2,SymLabelCountsCum2
        use SymData, only: SymLabelIntsCum2,nSymLabels,StateSymMap2
        use SymData, only: tagSymLabelIntsCum2,tagStateSymMap2,tagSymLabelCountsCum2
        use global_utilities
        use HElem, only: HElement_t_size
        IMPLICIT NONE
        integer Nirrep,nBasisfrz,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer iSize, iunit
        character(len=*),parameter :: thisroutine='SetupTMAT2'
        
        ! If this is a CPMD k-point calculation, then we're operating
        ! under Abelian symmetry: can use George's memory efficient
        ! TMAT.
        if (tCPMD) tCPMDSymTMat=tKP
        if (tVASP) tCPMDSymTMat=.true.
        IF(TSTARSTORE.or.tCPMDSymTMat) THEN 
            ! Set up info for indexing scheme (see TMatInd for full description).
            Nirrep=NSYMLABELS
            nBi=nBasisFRZ/iSS
            iSize=0
            if (associated(SymLabelIntsCum2)) then
                deallocate(SymLabelIntsCum2)
                call LogMemDealloc(thisroutine,tagSymLabelIntsCum2)
            end if
            if (allocated(StateSymMap2)) then
                deallocate(StateSymMap2)
                call LogMemDealloc(thisroutine,tagStateSymMap2)
            end if
            allocate(SymLabelIntsCum2(nIrrep))
            call LogMemAlloc('SymLabelIntsCum2',nIrrep,4,thisroutine,tagSymLabelIntsCum2)
            allocate(SymLabelCountsCum2(nIrrep))
            call LogMemAlloc('SymLabelCountsCum2',nIrrep,4,thisroutine,tagSymLabelCountsCum2)
            allocate(StateSymMap2(nBi))
            call LogMemAlloc('StateSymMap2',nBi,4,thisroutine,tagStateSymMap2)
            SYMLABELINTSCUM2(1:Nirrep)=0
            SYMLABELCOUNTSCUM2(1:Nirrep)=0
            do i=1,Nirrep
            !SYMLABELCOUNTS is now mbas only for the frozen orbitals
                basirrep=SYMLABELCOUNTS(2,i)
                iSize=iSize+(basirrep*(basirrep+1))/2
                ! # of integrals in this symmetry class and all preceeding
                ! symmetry classes.
                SYMLABELINTSCUM2(i)=iSize
                ! JSS: we no longer need SymLabelCountsCum, but it's a useful test.
                ! SymLabelCountsCum is the cumulative total # of basis functions
                ! in the preceeding symmetry classes.
                IF(i.eq.1) THEN
                    SYMLABELCOUNTSCUM2(i)=0
                ELSE
                    DO t=1,(i-1)
                        SYMLABELCOUNTSCUM2(i)=SYMLABELCOUNTSCUM2(i)+SYMLABELCOUNTS(2,t)
                    ENDDO
                ENDIF
!                write(6,*) basirrep,SYMLABELINTSCUM(i),SYMLABELCOUNTSCUM(i)
!                call neci_flush(6)
                ! JSS: Label states of symmetry i by the order in which they come.
                nStateIrrep=0
                do iState=1,nBi
                    if (SymClasses2(iState).eq.i) then
                        nStateIrrep=nStateIrrep+1
                        StateSymMap2(iState)=nStateIrrep
                    end if
                end do
            enddo
            IF((SYMLABELCOUNTSCUM2(Nirrep)+basirrep).ne.nBI) THEN
                iunit = get_free_unit()
                open(iunit, file="SYMLABELCOUNTS", status="unknown")
                DO i=1,Nirrep
                    WRITE(iunit,*) SYMLABELCOUNTS(2,i)
                    CALL neci_flush(iunit)
                ENDDO
                STOP 'Not all basis functions found while setting up TMAT2'
            ENDIF
            !iSize=iSize+2
            !This is to allow the index of '-1' in the array to give a zero value
            !Refer to TMatSym(-1) for when <i|h|j> is zero by symmetry.
            
            Allocate(TMATSYM2(-1:iSize),STAT=ierr)
            CALL LogMemAlloc('TMatSym2',iSize+2,HElement_t_size*8,thisroutine,tagTMATSYM2,ierr)
            TMATSYM2=(0.0_dp)

        ELSE

            if(tOneElecDiag) then
                ! In the UEG, the orbitals are eigenfunctions of KE operator, so TMAT is diagonal. 
                iSize=nBasisFRZ
                Allocate(TMAT2D2(nBasisFRZ,1),STAT=ierr)
                call LogMemAlloc('TMAT2D2',nBasisFRZ,HElement_t_size*8,thisroutine,tagTMat2D2)
                TMAT2D2=(0.0_dp)
            else
                ! Using a square array to hold <i|h|j> (incl. elements which are
                ! zero by symmetry).
                iSize=nBasisFRZ*nBasisFRZ
                Allocate(TMAT2D2(nBasisFRZ,nBasisFRZ),STAT=ierr)
                call LogMemAlloc('TMAT2D2',nBasisFRZ*nBasisFRZ,HElement_t_size*8,thisroutine,tagTMat2D2)
                TMAT2D2=(0.0_dp)
            endif
        
        ENDIF
      END SUBROUTINE SetupTMat2
    


      SUBROUTINE DestroyTMat(NEWTMAT)
        ! In:
        !    NewTMat : if true, destroy arrays used in storing "new" (post-freezing) TMat,
        !              else destroy those used for the pre-freezing TMat.
        use global_utilities
        IMPLICIT NONE
        LOGICAL :: NEWTMAT
        character(len=*), parameter :: thisroutine='DestroyTMat'

        IF(TSTARSTORE) THEN
            IF(NEWTMAT) THEN
                IF(ASSOCIATED(TMATSYM2)) THEN
                    CALL LogMemDealloc(thisroutine,tagTMatSym2)
                    Deallocate(TMATSYM2)
                    NULLIFY(TMATSYM2)
                ENDIF
            ELSE
                IF(ASSOCIATED(TMATSYM)) THEN
                    CALL LogMemDealloc(thisroutine,tagTMatSym)
                    Deallocate(TMATSYM)
                    NULLIFY(TMATSYM)
                ENDIF
            ENDIF
        ELSE
            IF(NEWTMAT) THEN
                IF(ASSOCIATED(TMAT2D2)) THEN
                    call LogMemDealloc(thisroutine,tagTMat2D2)
                    Deallocate(TMAT2D2)
                    NULLIFY(TMAT2D2)
                ENDIF
            ELSE
                IF(ASSOCIATED(TMAT2D)) THEN
                    Deallocate(TMAT2D)
                    call LogMemDealloc(thisroutine,tagTMat2D)
                    NULLIFY(TMAT2D)
                ENDIF
            ENDIF
        ENDIF
      END SUBROUTINE DestroyTMat

   


      SUBROUTINE SwapTMat(NBASIS,NHG,GG)
        ! In:
        !    nBasis: the number of active orbitals post-freezing.
        !    NHG: the number of orbitals used initially (pre-freezing).
        !    GG: GG(I) is the new (post-freezing) index of the old
        !        (pre-freezing) orbital I.
        ! During freezing, we need to know the TMat arrays both pre- and
        ! post-freezing.  Once freezing is done, clear all the pre-freezing
        ! arrays and point them to the post-freezing arrays, so the code
        ! referencing pre-freezing arrays can be used post-freezing.
        use constants, only: dp
        USE UMatCache
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymLabelCountsCum,SymLabelIntsCum
        use SymData, only: SymLabelCountsCum2,SymLabelIntsCum2
        use SymData, only: tagSymLabelCountsCum,tagSymLabelIntsCum
        use SymData, only: SymClasses2,tagSymClasses2
        use sym_mod
        use global_utilities
        IMPLICIT NONE
        integer NBASIS,NHG,GG(NHG)
        character(*),parameter :: this_routine='SwapTMat'

         IF(TSTARSTORE) THEN
            !Transfer across the vectors for indexing TMAT
            deallocate(SymClasses2)
            call LogMemDealloc(this_routine,tagSymClasses2)
            !copy across the new frozen symclasses and symlabelstuff
            CALL FREEZESYMLABELS(NHG,NBASIS,GG,.false.)
            !Redo SYMLABELCOUNTS
            CALL GENSymStatePairs(NBASIS/2,.false.)
            deallocate(SYMLABELCOUNTSCUM)
            call LogMemDealloc(this_routine,tagSYMLABELCOUNTSCUM)
            deallocate(SYMLABELINTSCUM)
            call LogMemDealloc(this_routine,tagSYMLABELINTSCUM)
            SYMLABELCOUNTSCUM=>SYMLABELCOUNTSCUM2
            nullify(SymLabelCountsCum2)
            SYMLABELINTSCUM=>SYMLABELINTSCUM2
            nullify(SymLabelIntsCum2)
             !Deallocate TMAT & reallocate with right size
             CALL DestroyTMAT(.false.)
             !CALL SetupTMAT(NBASIS,2,TMATINT)
             !DO LG=1,TMATINT
             !    TMATSYM(LG)=TMATSYM2(LG)
             !ENDDO
             !CALL DestroyTMAT(.true.)
             TMATSYM => TMATSYM2
             NULLIFY(TMATSYM2)
             !Swap the INVBRR and deallocate old one
             CALL LogMemDealloc(this_routine,tagINVBRR)
             Deallocate(INVBRR)
             INVBRR => INVBRR2
             NULLIFY(INVBRR2)
         ELSE
             !Deallocate TMAT & reallocate with right size
             CALL DestroyTMAT(.false.)
             !CALL SetupTMAT(NBASIS,2,TMATINT)
             !IF(TMATINT.ne.(NBASIS*NBASIS)) STOP 'Errorfrz'
             TMAT2D => TMAT2D2
             NULLIFY(TMAT2D2)
!             DO I=1,NBASIS
!                 DO J=1,NBASIS
!                     TMAT2D(NBASIS,NBASIS)=
!     &               TMAT2D2(NBASIS,NBASIS)
!                 ENDDO
!             ENDDO
             !CALL DestroyTMAT(.true.)
         ENDIF
      END SUBROUTINE SwapTMat
      
      

end module OneEInts
