module OneEInts

! For calculating, storing and retrieving the one electron integrals, <i|h|j>, where 
! h is the one-particle Hamiltonian operator and i and j are spin-orbitals.
! Variables associated with this usually have TMAT in the name.

USE HElem
USE System, only: TSTARSTORE

implicit none

save
public

TYPE(HElement), dimension(:,:), POINTER :: TMAT2D
TYPE(HElement), dimension(:), POINTER :: TMATSYM
TYPE(HElement), dimension(:), POINTER :: TMATSYM2
TYPE(HElement), dimension(:,:), POINTER :: TMAT2D2
logical tCPMDSymTMat

! Memory book-keeping tags
integer :: tagTMat2D=0
integer :: tagTMat2D2=0

contains

      !Get the index of TMAT element h_IJ (I&J are spin-orbs). This is only used with TSTARSTORE, where the TMAT is compressed to only store states, not spin-orbitals,
      !Added compression supplied by only storing symmetry allowed integrals - therefore needs sym.inc info.
      ! We assume a restricted calculation.  We note that TMat is a Hermitian matrix.
      ! For the TMat(i,j) to be non-zero, i and j have to belong to the same symmetry.
      ! We store the non-zero elements in TMatSym(:).
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
      !   If the element is zero, return -1 (TMatSym(-1) is set to 0). 
      INTEGER FUNCTION TMatInd(I,J)
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        IMPLICIT NONE
        INTEGER I,J,A,B,symI,symJ,Block,ind,cumulative,K,L
        include 'sym.inc'
        A=mod(I,2)
        B=mod(J,2)
        !If TMatInd = -1, then the spin-orbitals have different spins, or are symmetry disallowed therefore have a zero integral (apart from in UHF - might cause problem if we want this)
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
     


      ! See notes for TMatInd.
      INTEGER FUNCTION NEWTMatInd(I,J)
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        IMPLICIT NONE
        INTEGER I,J,A,B,symI,symJ,Block,ind,cumulative,K,L
        include 'sym.inc'
        A=mod(I,2)
        B=mod(J,2)
        !If TMatInd = -1, then the spin-orbitals have different spins, or are symmetry disallowed therefore have a zero integral (apart from in UHF - might cause problem if we want this)
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



      !Input at spin orbitals
      FUNCTION GetTMatEl(I,J)
        IMPLICIT NONE
        INTEGER I,J
        TYPE(HElement) GetTMatEl

        IF(TSTARSTORE) THEN
            GetTMatEl=TMATSYM(TMatInd(I,J))
        else if (tCPMDSymTMat) then
            ! TMat is Hermitian, rather than symmetric.
            ! Only the upper diagonal of each symmetry block stored.
            if (j.ge.i) then
                GetTMatEl=TMATSYM(TMatInd(I,J))
            else
                GetTMatEl=dConjg(TMATSYM(TMatInd(I,J)))
            end if
        ELSE
            GetTMatEl=TMAT2D(I,J)
        ENDIF
      END FUNCTION GetTMatEl



      ! See GetTMatEl.
      FUNCTION GetNewTMatEl(I,J)
        IMPLICIT NONE
        INTEGER I,J
        TYPE(HElement) GetNEWTMATEl

        IF(TSTARSTORE) THEN
            GetNEWTMATEl=TMATSYM2(NEWTMATInd(I,J))
        else if (tCPMDSymTMat) then
            if (j.ge.i) then
                GetNewTMatEl=TMATSYM2(TMatInd(I,J))
            else
                GetNewTMatEl=dConjg(TMATSYM2(TMatInd(I,J)))
            end if
        ELSE
            GetNEWTMATEl=TMAT2D2(I,J)
        ENDIF
      END FUNCTION GetNewTMatEl


      SUBROUTINE WriteTMat(NBASIS)
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        IMPLICIT NONE
        include 'sym.inc'
        INTEGER II,I,J,NBASIS
        
        IF(IP_SYMLABELINTSCUM.ne.0) THEN
            write(12,*) "SYMLABELCOUNTS,SYMLABELCOUNTSCUM,SYMLABELINTSCUM:"
            DO I=1,NSYMLABELS
                WRITE(12,"(I5,$)") SYMLABELCOUNTS(2,I)
                CALL FLUSH(12)
            ENDDO
            WRITE(12,*) ""
            DO I=1,NSYMLABELS
                WRITE(12,"(I5,$)") SYMLABELCOUNTSCUM(I)
                CALL FLUSH(12)
            ENDDO
            WRITE(12,*) ""
            DO I=1,NSYMLABELS
                WRITE(12,"(I5,$)") SYMLABELINTSCUM(I)
                CALL FLUSH(12)
            ENDDO
            WRITE(12,*) ""
            WRITE(12,*) "**********************************"
        ENDIF
        IF(IP_SYMLABELINTSCUM2.ne.0) THEN
            write(12,*) "SYMLABELCOUNTS,SYMLABELCOUNTSCUM2,SYMLABELINTSCUM2:"
            DO I=1,NSYMLABELS
                WRITE(12,"(I5,$)") SYMLABELCOUNTS(2,I)
                CALL FLUSH(12)
            ENDDO
            WRITE(12,*) ""
            DO I=1,NSYMLABELS
                WRITE(12,"(I5,$)") SYMLABELCOUNTSCUM2(I)
                CALL FLUSH(12)
            ENDDO
            WRITE(12,*) ""
            DO I=1,NSYMLABELS
                WRITE(12,"(I5,$)") SYMLABELINTSCUM2(I)
                CALL FLUSH(12)
            ENDDO
            WRITE(12,*) ""
            WRITE(12,*) "**********************************"
            CALL FLUSH(12)
        ENDIF
        WRITE(12,*) "TMAT:"
        IF(TSTARSTORE) THEN
            DO II=1,NSYMLABELS
                DO I=SYMLABELCOUNTSCUM(II-1)+1,SYMLABELCOUNTSCUM(II)
                    DO J=SYMLABELCOUNTSCUM(II-1)+1,I
                        WRITE(12,*) I,J,DREAL(GetTMATEl((2*I),(2*J)))
                        CALL FLUSH(12)
                    ENDDO
                ENDDO
            ENDDO
        ELSE
            DO I=1,NBASIS,2
                DO J=1,NBASIS,2
                    WRITE(12,*) (I+1)/2,(J+1)/2,DREAL(GetTMATEl(I,J))
                ENDDO
            ENDDO
        ENDIF
        WRITE(12,*) "**********************************"
        CALL FLUSH(12)
        IF(ASSOCIated(TMATSYM2).or.ASSOCIated(TMAT2D2)) THEN
            WRITE(12,*) "TMAT2:"
            DO II=1,NSYMLABELS
                DO I=SYMLABELCOUNTSCUM(II-1)+1,SYMLABELCOUNTSCUM(II)
                    DO J=SYMLABELCOUNTSCUM(II-1)+1,I
                        WRITE(12,*) I,J,DREAL(GetNEWTMATEl((2*I),(2*J)))
                        CALL FLUSH(12)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF
        WRITE(12,*) "*********************************"
        WRITE(12,*) "*********************************"
        CALL FLUSH(12)
      END SUBROUTINE WriteTMat
        


      SUBROUTINE SetupTMAT(nBASIS,iSS,iSize)   
        use System, only: tCPMD
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        use MemoryManager, only: LogMemAlloc
        IMPLICIT NONE
        include 'cpmddata.inc'
        include 'sym.inc'
        integer Nirrep,nBasis,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer iSize
        character(len=*),parameter :: thisroutine='SetupTMAT'
        
        ! If this is a CPMD k-point calculation, then we're operating
        ! under Abelian symmetry: can use George's memory efficient
        ! TMAT.  Bit of a hack for now: should place this line somewhere more
        ! relevant.
        if (tCPMD) tCPMDSymTMat=tKP
        IF(TSTARSTORE.or.tCPMDSymTMat) THEN 
            ! Set up info for indexing scheme (see TMatInd for full description).
            Nirrep=NSYMLABELS
            nBi=nBasis/iSS
            iSize=0
            IF(IP_SYMLABELINTSCUM.ne.0) CALL FREEM(IP_SYMLABELINTSCUM)
            IF(IP_SYMLABELCOUNTSCUM.ne.0) CALL FREEM(IP_SYMLABELCOUNTSCUM)
            IF(IP_StateSymMap.ne.0) call FreeM(IP_StateSymMap) ! I feel dirty doing this.
            CALL MEMORY(IP_SYMLABELINTSCUM,Nirrep,'SYMLABELINTSCUM')
            CALL IAZZERO(SYMLABELINTSCUM,Nirrep)
            CALL MEMORY(IP_SYMLABELCOUNTSCUM,Nirrep,'SYMLABELCOUNTSCUM')
            CALL IAZZERO(SYMLABELCOUNTSCUM,Nirrep)
            call Memory(IP_StateSymMap,nBi,'StateSymMap')
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
                DO i=1,Nirrep
                    WRITE(12,*) SYMLABELCOUNTSCUM(i)
                ENDDO
                write(12,*) "***************"
                write(12,*) NBI
                CALL FLUSH(12)
                STOP 'Not all basis functions found while setting up TMAT'
            ENDIF
            !iSize=iSize+2
            !This is to allow the index of '-1' in the array to give a zero value
            !Refer to TMatSym(-1) for when <i|h|j> is zero by symmetry.
            
            Allocate(TMATSYM(-1:iSize),STAT=ierr)
            Call MemAlloc(ierr,TMATSYM,HElementSize*(iSize+2),'TMATSYM')
            Call AZZERO(TMATSYM,HElementSize*(iSize+2))

        ELSE

            ! Using a square array to hold <i|h|j> (incl. elements which are
            ! zero by symmetry).
            iSize=nBasis*nBasis
            Allocate(TMAT2D(nBasis,nBasis),STAT=ierr)
            call LogMemAlloc('TMAT2D',nBasis*nBasis,HElementSize*8,thisroutine,tagTMat2D)
            Call AZZERO(TMAT2D,HElementSize*iSize)
        
        ENDIF
    
      END SUBROUTINE SetupTMAT


     
      ! See notes in SetupTMat as well.
      SUBROUTINE SetupTMAT2(nBASISFRZ,iSS,iSize)
        use System, only: tCPMD
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        use MemoryManager, only: LogMemAlloc
        IMPLICIT NONE
        include 'cpmddata.inc'
        include 'sym.inc'
        integer Nirrep,nBasisfrz,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer iSize
        character(len=*),parameter :: thisroutine='SetupTMAT2'
        
        ! If this is a CPMD k-point calculation, then we're operating
        ! under Abelian symmetry: can use George's memory efficient
        ! TMAT.
        if (tCPMD) tCPMDSymTMat=tKP
        IF(TSTARSTORE.or.tCPMDSymTMat) THEN 
            ! Set up info for indexing scheme (see TMatInd for full description).
            Nirrep=NSYMLABELS
            nBi=nBasisFRZ/iSS
            iSize=0
            IF(IP_SYMLABELINTSCUM2.ne.0) CALL FREEM(IP_SYMLABELINTSCUM2)
            IF(IP_StateSymMap2.ne.0) call FreeM(IP_StateSymMap2) ! I feel dirty doing this.
            CALL MEMORY(IP_SYMLABELINTSCUM2,Nirrep,'SYMLABELINTSCUM2')
            CALL IAZZERO(SYMLABELINTSCUM2,Nirrep)
            CALL MEMORY(IP_SYMLABELCOUNTSCUM2,Nirrep,'SYMLABELCOUNTSCUM2')
            CALL IAZZERO(SYMLABELCOUNTSCUM2,Nirrep)
            call Memory(IP_StateSymMap2,nBi,'StateSymMap2')
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
!                call flush(6)
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
                DO i=1,Nirrep
                    WRITE(14,*) SYMLABELCOUNTS(2,i)
                    CALL FLUSH(14)
                ENDDO
                STOP 'Not all basis functions found while setting up TMAT2'
            ENDIF
            !iSize=iSize+2
            !This is to allow the index of '-1' in the array to give a zero value
            !Refer to TMatSym(-1) for when <i|h|j> is zero by symmetry.
            
            Allocate(TMATSYM2(-1:iSize),STAT=ierr)
            CALL MemAlloc(ierr,TMATSYM2,HElementSize*(iSize+2),'TMATSYM2')
            CALL AZZERO(TMATSYM2,HElementSize*(iSize+2))

        ELSE

            ! Using a square array to hold <i|h|j> (incl. elements which are
            ! zero by symmetry).
            iSize=nBasisFRZ*nBasisFRZ
            Allocate(TMAT2D2(nBasisFRZ,nBasisFRZ),STAT=ierr)
            call LogMemAlloc('TMAT2D2',nBasisFRZ*nBasisFRZ,HElementSize*8,thisroutine,tagTMat2D2)
            Call AZZERO(TMAT2D2,HElementSize*iSize)
        
        ENDIF
      END SUBROUTINE SetupTMat2
    


      SUBROUTINE DestroyTMat(NEWTMAT)
        use MemoryManager, only: LogMemDealloc
        IMPLICIT NONE
        LOGICAL :: NEWTMAT
        character(len=*), parameter :: thisroutine='DestroyTMat'

        IF(TSTARSTORE) THEN
            IF(NEWTMAT) THEN
                IF(ASSOCIATED(TMATSYM2)) THEN
                    CALL MemDealloc(TMATSYM2)
                    Deallocate(TMATSYM2)
                    NULLIFY(TMATSYM2)
                ENDIF
            ELSE
                IF(ASSOCIATED(TMATSYM)) THEN
                    CALL MemDealloc(TMATSYM)
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
                    call LogMemDealloc(thisroutine,tagTMat2D)
                    Deallocate(TMAT2D)
                    NULLIFY(TMAT2D)
                ENDIF
            ENDIF
        ENDIF
      END SUBROUTINE DestroyTMat

   


      SUBROUTINE SwapTMat(NBASIS,NHG,GG)
        USE HElem
        USE UMatCache
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB
        use System, only: BasisFN,BasisFNSize,BasisFNSizeB
        IMPLICIT NONE
        include 'sym.inc'
        integer I,J,NBASIS,NHG,GG(NHG)
        integer*8 TMATINT,LG
        
         IF(TSTARSTORE) THEN
            !Transfer across the vectors for indexing TMAT
            CALL FREEM(IP_SYMCLASSES2)
            !copy across the new frozen symclasses and symlabelstuff
            CALL FREEZESYMLABELS(NHG,NBASIS,GG,.false.)
            !Redo SYMLABELCOUNTS
            CALL GENSymStatePairs(NBASIS/2,.false.)
            CALL FREEM(IP_SYMLABELCOUNTSCUM)
            IP_SYMLABELCOUNTSCUM=IP_SYMLABELCOUNTSCUM2
            IP_SYMLABELCOUNTSCUM2=0
            CALL FREEM(IP_SYMLABELINTSCUM)
            IP_SYMLABELINTSCUM=IP_SYMLABELINTSCUM2
            IP_SYMLABELINTSCUM2=0
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
             CALL MemDealloc(INVBRR)
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
