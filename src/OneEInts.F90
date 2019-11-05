module OneEInts

! For calculating, storing and retrieving the one electron integrals, <i|h|j>, where
! h is the one-particle Hamiltonian operator and i and j are spin-orbitals.
! Variables associated with this usually have TMAT in the name.

! The 2 or NEW versions of routines and variables are used for the set of
! orbitals after freezing has been done. The "original" versions are used in the
! pre-freezing stage.

use constants, only: dp
use MemoryManager, only: TagIntType, LogMemalloc, LogMemDealloc
use util_mod, only: get_free_unit
use SystemData, only: Symmetry, BasisFN, tCPMD, tVASP
use LoggingData, only: iNumPropToEst
use CPMDData, only: tKP
use HElem, only: HElement_t_size
use global_utilities
use UMatCache, only: nStates
use SymData

implicit none

save
public

! For storing the <i|h|j> elements which are non-zero by symmetry.
! Ranges from -1 to N, where N is the # of non-zero elements. The -1
! element is used to return the zero value.  Used for Abelian symmetries,
! where i and j must span the same representation in order for <i|h|j>
! to be non-zero.
! The symmetries used in Dalton and Molpro are all Abelian.
HElement_t(dp), dimension(:), POINTER :: TMATSYM
! For non-Abelian symmetries, we store the entire <i|h|j> matrix, as
! the direct products of the representations can contain the totally
! symmetric representation (and are not necessarily the same).  We could
! compress this in a similar fashion at some point.
HElement_t(dp), dimension(:,:), POINTER :: TMAT2D

HElement_t(dp), dimension(:), POINTER :: TMATSYM2
HElement_t(dp), dimension(:,:), POINTER :: TMAT2D2

! One electron integrals corresponding to properties that would be estimated
! using 1 body RDM
HElement_t(dp), dimension(:,:,:), pointer :: OneEPropInts
real(dp), dimension(:), pointer :: PropCore

! A second pointer to get the integrals after freezing orbitals
HElement_t(dp), dimension(:,:,:), pointer :: OneEPropInts2

! True if using TMatSym in CPMD (currently only if using k-points, which form
! an Abelian group).
logical tCPMDSymTMat
logical tOneElecDiag    !Indicates that the one-electron integral matrix is diagonal -
                        !basis functions are eigenstates of KE operator.

! Memory book-keeping tags
integer(TagIntType) :: tagTMat2D=0
integer(TagIntType) :: tagTMat2D2=0
integer(TagIntType) :: tagTMATSYM=0,tagTMATSYM2=0
integer(TagIntType) :: tagOneEPropInts=0
integer(TagIntType) :: tagOneEPropInts2=0
integer(TagIntType) :: tagPropCore=0


contains

      pure INTEGER FUNCTION TMatInd(I,J)
        ! In:
        !    i,j: spin orbitals.
        ! Return the index of the <i|h|j> element in TMatSym2.
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
        HElement_t(dp) :: ret
#ifdef __CMPLX
        HElement_t(dp) :: t
#endif

        if (tCPMDSymTMat) then
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

    function GetPropIntEl(i,j,iprop) result(integral)

!       use OneEInts, only: OneEPropInts
        integer, intent(in) :: i, j, iprop
        real(dp) :: integral

        integral = OneEPropInts(i,j,iprop)

    end function

      FUNCTION GetNewTMatEl(I,J)
      ! In:
      !    i,j: spin orbitals.
      ! Return <i|h|j>, the "TMat" element.
      ! Used post-freezing. See also GetTMatEl.
        IMPLICIT NONE
        INTEGER I,J
        HElement_t(dp) GetNEWTMATEl

        if (tCPMDSymTMat) then
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
        DO I=1,NBASIS,2
            DO J=1,NBASIS,2
                WRITE(iunit,*) (I+1)/2,(J+1)/2, GetTMATEl(I,J)
            ENDDO
        ENDDO
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
      INTEGER :: iSize,nBasis,basirrep,i,Nirrep

          IF(tCPMDSymTMat) THEN
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
        IMPLICIT NONE
        integer Nirrep,nBasis,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer iSize, iunit
        character(*), parameter :: thisroutine = 'SetupTMAT'

        ! If this is a CPMD k-point calculation, then we're operating
        ! under Abelian symmetry: can use George's memory efficient
        ! TMAT.
        if (tCPMD) tCPMDSymTMat=tKP
        if (tVASP) tCPMDSymTMat=.true.
        IF(tCPMDSymTMat) THEN
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

            ! Ensure no compiler warnings
            basirrep=SYMLABELCOUNTS(2, 1)

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
                call stop_all(thisroutine, 'Not all basis functions found while setting up TMAT')
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

      subroutine SetupPropInts(nBasis)
        implicit none
        integer, intent(in) :: nBasis
        integer :: ierr,iSize
        character(*),parameter :: t_r = 'SetupPropertyInts'

        ! Using a square array to hold <i|h|j> (incl. elements which are
        ! zero by symmetry).
        Allocate(OneEPropInts(nBasis,nBasis,iNumPropToEst),STAT=ierr)
        iSize = NBasis*NBasis*iNumPropToEst
        call LogMemAlloc('OneEPropInts',nBasis*nBasis*iNumPropToEst,HElement_t_size*8,t_r,tagOneEPropInts)
        OneEPropInts = (0.0_dp)
        Allocate(PropCore(iNumPropToEst),STAT=ierr)
        call LogMemAlloc('PropCore',iNumPropToEst,dp,t_r,tagPropCore)
        PropCore = 0.0d0

      end subroutine SetupPropInts

      subroutine SetupPropInts2(nBasisFrz)
        implicit none
        integer, intent(in) :: nBasisFrz
        integer :: ierr,iSize
        character(*),parameter :: t_r = 'SetupPropertyInts2'

        ! Using a square array to hold <i|h|j> (incl. elements which are
        ! zero by symmetry).
        Allocate(OneEPropInts2(nBasisFrz,nBasisFrz,iNumPropToEst),STAT=ierr)
        iSize = NBasisFrz*NBasisFrz*iNumPropToEst
        call LogMemAlloc('OneEPropInts2',iSize,HElement_t_size*8,t_r,tagOneEPropInts2)
        OneEPropInts2 = (0.0_dp)

      end subroutine SetupPropInts2

      SUBROUTINE SetupTMAT2(nBASISFRZ,iSS,iSize)
        ! In:
        !    nBasisFRZ: number of active basis functions (orbitals).
        !    iSS: ratio of nBasisFRZ to the number of active spatial orbitals.
        !    iSize: number of elements in TMat2/TMatSym2.
        ! Initial allocation of TMat2D2 or TMatSym2 (if using symmetry-compressed
        ! storage of the <i|h|j> integrals) for post-freezing.
        ! See also notes in SetupTMat.
        IMPLICIT NONE
        integer Nirrep,nBasisfrz,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer iSize, iunit
        character(*),parameter :: thisroutine='SetupTMAT2'

        ! If this is a CPMD k-point calculation, then we're operating
        ! under Abelian symmetry: can use George's memory efficient
        ! TMAT.
        if (tCPMD) tCPMDSymTMat=tKP
        if (tVASP) tCPMDSymTMat=.true.
        IF(tCPMDSymTMat) THEN
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

            ! Ensure no compiler warnings
            basirrep=SYMLABELCOUNTS(2, 1)

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
                call stop_all(thisroutine, 'Not all basis functions found while setting up TMAT2')
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
        IMPLICIT NONE
        LOGICAL :: NEWTMAT
        character(len=*), parameter :: thisroutine='DestroyTMat'

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
      END SUBROUTINE DestroyTMat

      SUBROUTINE DestroyPropInts()
        implicit none
        character(*),parameter :: t_r = 'DestroyPropInts'

        Deallocate(OneEPropInts)
        call LogMemDealloc(t_r,tagOneEPropInts)
        if (associated(OneEPropInts2)) then
            call LogMemDealloc(t_r,tagOneEPropInts2)
            Deallocate(OneEPropInts2)
        end if
        Deallocate(PropCore)

      END SUBROUTINE DestroyPropInts

      SUBROUTINE SwapTMat()
        ! In:
        !    nBasis: the number of active orbitals post-freezing.
        !    NHG: the number of orbitals used initially (pre-freezing).
        !    GG: GG(I) is the new (post-freezing) index of the old
        !        (pre-freezing) orbital I.
        ! During freezing, we need to know the TMat arrays both pre- and
        ! post-freezing.  Once freezing is done, clear all the pre-freezing
        ! arrays and point them to the post-freezing arrays, so the code
        ! referencing pre-freezing arrays can be used post-freezing.
!         use sym_mod
        IMPLICIT NONE
        character(*),parameter :: this_routine='SwapTMat'

        ! Deallocate TMAT & reallocate with right size
        CALL DestroyTMAT(.false.)
        TMAT2D => TMAT2D2
        NULLIFY(TMAT2D2)
      END SUBROUTINE SwapTMat

      SUBROUTINE SwapOneEPropInts(nBasisFrz,iNum)

        ! IN: iNum is the number of perturbation operator used in the calculation
        ! During freezing, we need to know the OneEPropInts arrays both pre- and
        ! post-freezing.  Once freezing is done, clear all the pre-freezing
        ! arrays and point them to the post-freezing arrays, so the code
        ! referencing pre-freezing arrays can be used post-freezing.
        implicit none
        integer, intent(in) :: nBasisFrz,iNum
        integer :: iSize, ierr
        character(*),parameter :: t_r = 'SwapOneEPropInts'

        Deallocate(OneEPropInts)
        call LogMemDealloc(t_r,tagOneEPropInts)
        NULLIFY(OneEPropInts)

        Allocate(OneEPropInts(nBasisFrz,nBasisFrz,iNum),STAT=ierr)
        iSize = nBasisFrz*nBasisFrz*iNum
        call LogMemAlloc('OneEPropInts',iSize,HElement_t_size*8,t_r,tagOneEPropInts)

        OneEPropInts => OneEPropInts2

        NULLIFY(OneEPropInts2)

      END SUBROUTINE SwapOneEPropInts

end module OneEInts
