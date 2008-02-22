!Based on umatcache.F, this modulises the umat cache, and retrieval functions.
MODULE UMatCache
      USE HElem
      USE System , only : TSTARSTORE
      IMPLICIT NONE
! The cache is stored in UMatCacheData
!  For real systems, nTypes=1, and we have a single value for each unique ordered pair of ordered pairs (i,k), (j,l)
!    out of <ij|u|kl> (where i,j,k,l are state labels).
!  For complex systems, nTypes=2, giving two possible values.
!  For each nPairs (which corresponds to the greater of the two pairs), there  are nSlots, with the other pair label in !  UMatLabels, and the values in UMatCache.
!  nStates is the maximum number of states stored in the cache (which may not be all the states if there are frozen virtuals).
!
      TYPE(HElement), Pointer :: UMatCacheData(:,:,:) !(0:nTypes-1,nSlots,nPairs)
      INTEGER, Pointer :: UMatLabels(:,:) !(nSlots,nPairs)
      INTEGER nSlots,nPairs,nTypes
      INTEGER nStates
!tSmallUMat is set if we have nStates slots per pair for storing the <ik|jk> integrals.
!  This should only be used prior to freezing to store precalculated integrals.
!  For each pair (i,j), we store the <ik|jk> integral in slot k.
      LOGICAL tSmallUMat
         SAVE nSlots,nPairs,nTypes,nStates,tSmallUMat,UMatCacheData,UMatLabels
!  iDumpCacheFlag: Dumps the cache if we're told to.
!    =0: No cache dumping.
!    =1: Dump cache unless we've read in a cache dump that's from a larger
!    calculation than the current one.
!    =2: Force dumping of cache (over-writing previous cache level).
!  nStatesDump: number of states used in calculation which produced dump file.
!  tReadInCache: does what it says on the tin.
      integer iDumpCacheFlag,nStatesDump
      logical tReadInCache
         save iDumpCacheFlag,nStatesDump,tReadInCache


! For the more requently used <ij|u|ij> and <ij|u|ji> integrals, we store them in a separate cache (if TUMat2D is true)
      TYPE(HElement), Pointer :: UMat2D(:,:) !(nStates,nStates)
      LOGICAL tUMat2D
         SAVE tUMat2D ,UMat2D

! This is for the storage of the one-electron integrals
      TYPE(HElement), dimension(:,:), POINTER :: TMAT2D
      TYPE(HElement), dimension(:), POINTER :: TMATSYM
      TYPE(HElement), dimension(:), POINTER :: TMATSYM2
      TYPE(HElement), dimension(:,:), POINTER :: TMAT2D2
      logical tCPMDSymTMat
      SAVE TMAT2D,TMATSYM,TMATSYM2,TMAT2D2

! This vector stores the energy ordering for each spatial orbital, which is the inverse of the BRR vector
! This is needed for the memory saving star indexing system.
! E.g. Element 4 will give the the order in the energy of element 4
     INTEGER, DIMENSION(:), POINTER :: INVBRR
     INTEGER, DIMENSION(:), POINTER :: INVBRR2
     SAVE INVBRR,INVBRR2

!NOCC is number of occupied spatial orbitals - needed for test in UMATInd, thought would be quicker than passing it in each time.
!Freezetransfer is a temporary measure to tell UMATIND when the freezing of orbitals is occuring.
      INTEGER,SAVE :: NOCC
      LOGICAL,SAVE :: FREEZETRANSFER
! For the UEG, we damp the exchange interactions.
!    0 means none
!    1 means attenuated (using an erfc)
!    2 means cut-off    (at a distance Rc=ALAT(4))
      INTEGER iPeriodicDampingType

!  Book-keeping information
!  nSlotsInit is the number of slots requested on input.  If the number required is less, then the lower value is allocated
!  If nSlotsInit is set to 0, then general <ij|u|kl> element caching is not performed, but UMat2D <ij|u|ij> and <ij|u|ji> is.  For nSlotsInit=-1 neither is performed.
      INTEGER nSlotsInit,nMemInit
! nHits and nMisses are the number of cache hits and misses.
      INTEGER nHits,nMisses
!.. UMatCacheFlag=0 is normal operation
!.. UMatCacheFlag=1 means cache from bottom up
!..      This is useful when lots of sequential pieces of data are being stored.
!..  When UMatCacheFlag is reset to 0, the data which are present are spread evenly around the slots for a given Pair.
      INTEGER UMatCacheFlag
!.. The number of cache overwrites
      INTEGER iCacheOvCount
         SAVE nSlotsInit,nMemInit,nHits,nMisses,UMatCacheFlag,iCacheOvCount

!.. Some various translation tables to convert between different orderings of states.
      LOGICAL tTransGTID,tTransFindx
      INTEGER, Pointer :: TransTable(:) !(NSTATES)
      INTEGER, Pointer :: InvTransTable(:) !(NSTATES)
         SAVE tTransGTID,tTransFindx,TransTable,InvTransTable


!.. Density fitting cache information.  If we're generating integrals on the fly from density fitting.
      integer nAuxBasis,nBasisPairs
      logical tDFInts      
      Real*8,Pointer :: DFCoeffs(:,:) !(nAuxBasis,nBasisPairs)
      Real*8,Pointer :: DFInts(:,:) !(nAuxBasis,nBasisPairs)
      Real*8,Pointer :: DFFitInts(:,:) !(nAuxBasis,nAuxBasis)
      Real*8,Pointer :: DFInvFitInts(:,:) !(nAuxBasis,nAuxBasis)
      INTEGER iDFMethod
!Some possible DFMethods sums over P, Q implied.  All precontracted to run in order(X) except DFOVERLAP2NDORD
! 0 - no DF
! DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)
! DFOVERLAP2NDORD  2 - (ij|u|ab)= (ij|u|P)(P|ab)+(ij|P)(P|u|ab)-(ij|P)(P|u|Q)(Q|ab)
! DFOVERLAP2       3 - (ij|u|ab)= (ij|P)(P|u|Q)(Q|ab)
! DFCOULOMB        4 - (ij|u|ab)= (ij|u|P)[(P|u|Q)^-1](Q|u|ab)

         SAVE nAuxBasis,nBasisPairs,tDFInts,DFCoeffs,DFInts,DFFitInts,DFInvFitInts,iDFMethod

      Contains

      !Create new INVBRR for the freezing process
      SUBROUTINE CREATEINVBRR2(BRR2,NBASIS)
        IMPLICIT NONE
        INTEGER BRR2(NBASIS),NBASIS,ierr,I,t

        ALLOCATE(INVBRR2(NBASIS/2),STAT=ierr)
        CALL MemAlloc(ierr,INVBRR2,NBASIS/2,'INVBRR2')
        CALL IAZZERO(INVBRR2,NBASIS/2)
        t=0
        DO I=2,NBASIS,2
            t=t+1
            INVBRR2(BRR2(I)/2)=t
        ENDDO
        RETURN
      END subroutine

      
      SUBROUTINE CREATEINVBRR(BRR,NBASIS)
        IMPLICIT NONE
        INTEGER BRR(NBASIS),NBASIS,ierr,I,t

        IF(ASSOCIATED(INVBRR)) THEN
            CALL MemDealloc(INVBRR)
            DEALLOCATE(INVBRR)
        ENDIF
        ALLOCATE(INVBRR(NBASIS/2),STAT=ierr)
        CALL MemAlloc(ierr,INVBRR,NBASIS/2,'INVBRR')
        CALL IAZZERO(INVBRR,NBASIS/2)
        t=0
        DO I=2,NBASIS,2
            t=t+1
            INVBRR(BRR(I)/2)=t
        ENDDO
!        DO I=1,NBASIS/2
!            WRITE(11,*) INVBRR(I)
!        ENDDO
        return
      END subroutine

      
! Get the index of physical order UMAT element <IJ|KL>.  Indices are internally reordered such that I>K, J>L,(I,K)>(J,L) 
!NB This is a different order from UMatCache. Orbitals are passed as arguments for spatial orbitals
!If NBASIS or NOCCUPIED is passed in as zero, the values defined in the module are used, which are set elsewhere.
      INTEGER FUNCTION UMatInd(I,J,K,L,NBASIS,NOCCUPIED)
         IMPLICIT NONE
         INTEGER I,J,K,L,AA,BB,NBASIS
         INTEGER R,S,T,U,A,B,C,D,NOCCUPIED
         IF(TSTARSTORE) THEN
            !Rearrange, so that orbitals ordered over energy, and first two indices are occupied
            !Could be a problem in the future r.e. partially filled degenerate fermi levels - is BRR then the best way to determine if an orbital is occupied or not??
            IF(NOCCUPIED.EQ.0) THEN
            
                R=INVBRR(I)
                S=INVBRR(J)
                T=INVBRR(K)
                U=INVBRR(L)
            ELSE
                R=INVBRR2(I)
                S=INVBRR2(J)
                T=INVBRR2(K)
                U=INVBRR2(L)
            ENDIF

            !Need to create unique pairings from all permutations of indices
            IF(R.le.T) THEN
                A=R
                C=T
            ELSE
                A=T
                C=R
            ENDIF
            IF(S.le.U) THEN
                B=S
                D=U
            ELSE
                B=U
                D=S
            ENDIF
            IF((A.lt.B).or.((A.eq.B).and.(C.le.D))) THEN
                R=A
                S=B
                T=C
                U=D
            ELSE
                R=B
                S=A
                T=D
                U=C
            ENDIF
            !During the freezing routine, it tries to lookup <ia|ib>, for the h_ab integrals where a & b are distinct and virtual, but these are unneeded if just considering double excitations.
            IF(FREEZETRANSFER.and.(S.gt.NOCC)) THEN
                UMatInd=-1
                RETURN
            ENDIF
            IF(NOCCUPIED.EQ.0) THEN
                IF((R.gt.NOCC).or.(S.gt.NOCC)) THEN
                    WRITE(6,*) "NO OCCUPIED ORBITAL PAIR REQUESTED - STARBINREAD CANNOT BE USED"
                    WRITE(6,*) "USING ORIGINAL UMAT AND STORED NOCC"
                    WRITE(6,*) "NOCC is: ",NOCC
                    WRITE(6,*) "SPIN-ORBITALS ",I*2,", ", J*2,", ",K*2,", ",L*2," requested."
                    CALL FLUSH(6)
                    STOP 'NO OCCUPIED ORBITAL PAIR REQUESTED'
                ENDIF
            ELSE
                IF((R.gt.NOCCUPIED).or.(S.gt.NOCCUPIED)) THEN
                    WRITE(6,*) "NO OCCUPIED ORBITAL PAIR REQUESTED - STARBINREAD CANNOT BE USED"
                    WRITE(6,*) "USING UMAT2 AND NOCC FROM ARGUMENT"
                    WRITE(6,*) "NOCC is: ",NOCC
                    WRITE(6,*) "SPIN-ORBITALS ",I*2,", ", J*2,", ",K*2,", ",L*2," requested."
                    CALL FLUSH(6)
                    STOP 'NO OCCUPIED ORBITAL PAIR REQUESTED'
                ENDIF
            ENDIF
                
                
            IF(NBASIS.ne.0) THEN
                BB=((S-1)*NBASIS)+U
                AA=((R-1)*NBASIS)+T
            ELSE
                BB=((S-1)*nStates)+U
                AA=((R-1)*nStates)+T
            ENDIF
                
            UMatInd=((BB*(BB-1))/2)+AA
            RETURN
         ELSE
            IF(I.GT.K) THEN
                A=(I*(I-1))/2+K
            ELSE
                A=(K*(K-1))/2+I
            ENDIF
         
            IF(J.GT.L) THEN
                B=(J*(J-1))/2+L 
            ELSE
                B=(L*(L-1))/2+J
            ENDIF
            IF(A.GT.B) THEN
                UMatInd=(A*(A-1))/2+B
                RETURN
            ELSE
                UMatInd=(B*(B-1))/2+A
                RETURN
            ENDIF
         ENDIF
      END function

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
      END function
      
     ! See notes for TMatInd.
     INTEGER FUNCTION NEWTMatInd(I,J)
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
      END function
     
      
!Get the prospective size of a UMat (not a UMatCache) for completely storing FCIDUMP 2-e integrals
!  The UMat is currently passed as a parameter, but in future be absorbed into UMatCache.
      SUBROUTINE GetUMatSize(nBasis,nEl,iSS,iSize)
         IMPLICIT NONE
         INTEGER nBasis,iSS
         INTEGER iPairs,nBi,nEl,noccup
         INTEGER*8 :: iSize
         nBi=nBasis/iSS
         IF(TSTARSTORE) THEN
            IF(MOD(nel,2).ne.0) THEN
                noccup=(nel+1)/iSS
            ELSE
                noccup=nel/iSS
            ENDIF
            iPairs=noccup*nBi
            iSize=(iPairs*(iPairs+1))/2 
         ELSE
            iPairs=(nBi*(nBi+1))/2
            iSize=(iPairs*(iPairs+1))/2
         ENDIF
      END subroutine
      
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
      END function

      ! See GetTMatEl.
      FUNCTION GetNEWTMATEl(I,J)
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
      END function
     
      ! See notes in SetupTMat as well.
      SUBROUTINE SetupTMAT2(nBASISFRZ,iSS,iSize)
        use System, only: tCPMD
        IMPLICIT NONE
        include 'cpmddata.inc'
        include 'sym.inc'
         integer Nirrep,nBasisfrz,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer*8 iSize
        
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
            Call MemAlloc(ierr,TMAT2D2,HElementSize*iSize,'TMAT2D2')
            Call AZZERO(TMAT2D2,HElementSize*iSize)
        
        ENDIF
      END subroutine
    
      SUBROUTINE DestroyTMAT(NEWTMAT)
        IMPLICIT NONE
        LOGICAL :: NEWTMAT

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
                    CALL MemDealloc(TMAT2D2)
                    Deallocate(TMAT2D2)
                    NULLIFY(TMAT2D2)
                ENDIF
            ELSE
                IF(ASSOCIATED(TMAT2D)) THEN
                    CALL MemDealloc(TMAT2D)
                    Deallocate(TMAT2D)
                    NULLIFY(TMAT2D)
                ENDIF
            ENDIF
        ENDIF
      END subroutine
      SUBROUTINE WRITETMAT(NBASIS)
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
      END subroutine
        
      SUBROUTINE SetupTMAT(nBASIS,iSS,iSize)   
        use System, only: tCPMD
        IMPLICIT NONE
        include 'cpmddata.inc'
        include 'sym.inc'
        integer Nirrep,nBasis,iSS,nBi,i,basirrep,t,ierr,iState,nStateIrrep
        integer*8 iSize
        
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
            Call MemAlloc(ierr,TMAT2D,HElementSize*iSize,'TMAT2D')
            Call AZZERO(TMAT2D,HElementSize*iSize)
        
        ENDIF
    
      END SUBROUTINE SetupTMAT

      !Get a U matrix element <ij|u|kl> in multifarious ways, where orbitals are spatial orbitals.  Either from a passed-in UMAT, or ALAT parameters, 
! or from UMatcache.
      FUNCTION GetUMatEl(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDI,IDJ,IDK,IDL)
         IMPLICIT NONE
         TYPE(HElement) GetUMatEl
         INTEGER NBASISMAX(5,3),I,J,K,L,NHG,ISS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(NHG)
         REAL*8 ALAT(3),GetNan
         TYPE(HElement) UMAT(*)
         TYPE(HElement) UElems(0:nTypes-1)
         INTEGER A,B,C,XXX
         INTEGER IDI,IDJ,IDK,IDL
         REAL*8 SUM
         real*8, PARAMETER :: PI=3.14159265358979323846264338327950288419716939937510D0
         INTEGER ICACHE,ICACHEI,ITYPE
         LOGICAL LSYMSYM
         TYPE(Symmetry) SYM,SYMPROD,SYMCONJ
         INTEGER ISUB,ISUB2
         LOGICAL GetCachedUMatEl,HasKPoints
         integer*8 TotSymRep
!         CALL TISET(' GETUMATEL',ISUB)
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
               SYM%s=TotSymRep()
               SYM=SYMPROD(SYM,SYMCONJ(G1(I*2-1)%Sym))
               SYM=SYMPROD(SYM,SYMCONJ(G1(J*2-1)%Sym))
               SYM=SYMPROD(SYM,G1(K*2-1)%Sym)
               SYM=SYMPROD(SYM,G1(L*2-1)%Sym)
!         WRITE(6,"(A,5I5)") "NN",IDI,IDJ,IDK,IDL,SYM%s
!   Check the symmetry of the 4-index integrals
              IF(.NOT.LSYMSYM(SYM)) THEN
                  GETUMATEL=0.D0
! JSS --- comment out the following TIHALT.  (TISET commented out above.)
!                 CALL TIHALT(' GETUMATEL',ISUB)
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
!            if (idi.eq.59.and.idj.eq.1.and.idk.eq.67.and.idl.eq.1)then
!                write (6,*) 
!            end if
!            if (idi.eq.59.and.idj.eq.2.and.idk.eq.67.and.idl.eq.2)then
!                write (6,*) 
!            end if
!            if (idi.eq.60.and.idj.eq.1.and.idk.eq.68.and.idl.eq.1)then
!                write (6,*) 
!            end if
                call KPntSymInt(I,J,K,L,I,J,K,L)
!                if (i.eq.63.and.j.eq.75.and.k.eq.47.and.l.eq.63) then
!                   write (6,*) 
!                end if
                IF(TTRANSFINDX) THEN
                 I=InvTransTable(I)
                 J=InvTransTable(J)
                 K=InvTransTable(K)
                 L=InvTransTable(L)
                ENDIF
               ENDIF
!    This will rearrange I,J,K,L into the correct order
!   (i,k)<=(j,l) and i<=k, j<=l.
!               WRITE(6,*) "J",I,J,K,L
!                        CALL INITFINDXI(I,J,K,L,UElems)
!               GetUMatEl=UElems(0)
!         WRITE(6,"(4I5,$)") IDI,IDJ,IDK,IDL
!         WRITE(6,*) GETUMATEL,ABS(GETUMATEL)
!               return
               IF(GETCACHEDUMATEL(I,J,K,L,GETUMATEL,ICACHE,ICACHEI,A,B,ITYPE)) THEN
!   We don't have a stored UMAT - we call to generate it.
                  IF(tDFInts) THEN
!   We're using density fitting
                     Call GetDF2EInt(I,J,K,L,UElems)
                     GetUMatEl=UElems(0)
                  ELSE
!   Otherwise we call CPMD
!                    write (6,*) TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L)
                     IF(TTRANSFINDX) THEN
!         WRITE(6,"(A,4I5,$)") "MM",TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L)
                        CALL INITFINDXI(TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L),UElems)
                     ELSE
                        CALL INITFINDXI(I,J,K,L,UElems)
!InitFindxI returns up to two integrals in UElems
!  <ij|u|kl> and <kj|u|il> (which are distinct when complex orbitals are used).
!  TYPE 0          TYPE 1

                     ENDIF
!  Bit 0 tells us which integral in the slot we need
                     GETUMATEL=UElems(IAND(ITYPE,1))
!  Bit 1 tells us whether we need to complex conj the integral
                     IF(BTEST(ITYPE,1)) GETUMATEL=DCONJG(GETUMATEL)
!                     WRITE(6,"(A3,6I6,$)") "I",I,J,K,L,A,B
!                     WRITE(6,*) UElems
!                     IF(TTRANSFINDX) THEN
!                        CALL INITFINDXI(TRANSTABLE(K),TRANSTABLE(J),TRANSTABLE(I),TRANSTABLE(L),UElems)
!                     ELSE
!                        CALL INITFINDXI(K,J,I,L,UElems)
!                     ENDIF
!                     WRITE(6,*) "I",K,J,I,L,UElems
!                     CALL INITFINDXI(I,L,K,J,UElems)
!                     WRITE(6,*) "I",I,L,K,J,UElems
!                     CALL INITFINDXI(IDI,IDJ,IDK,IDL,UElems)
!                     WRITE(6,*) "I",IDI,IDJ,IDK,IDL,UElems
!                     CALL INITFINDXI(8,4,5,8,UElems)
!                     WRITE(6,*) "I",8,4,5,8,UElems
                  ENDIF
!  Because we've asked for the integral in the form to be stored, we store as iType=0
                  IF(ICACHE.NE.0) CALL CACHEUMATEL(A,B,UElems,ICACHE,ICACHEI,0)
                  NMISSES=NMISSES+1
!                  WRITE(6,*) "MISS",I,J,K,L,A,B,GETUMATEL
               ELSE
!                  WRITE(6,*) "HIT ",IDI,IDJ,IDK,IDL,GETUMATEL
!                  WRITE(6,*) A,B
                  NHITS=NHITS+1
              ENDIF
             ENDIF
           ENDIF
          ELSEIF(ISS.EQ.-1) THEN
!  A  non-stored hubbard integral.
          CALL GetHubUMatEl(IDI,IDJ,IDK,IDL,UMat,nBasisMax,G1,GetUMatEl)
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
                         GETUMATEL=UMAT(XXX)
                     ELSE
                         !GETUMATEL=HElement(0.D0)
                         GETUMATEL=HElement(GETNAN())
                     ENDIF
                 ENDIF
             ELSE
                 GETUMATEL=UMAT(UMatInd(IDI,IDJ,IDK,IDL,0,0))
             ENDIF
          ENDIF
         ELSEIF(NBASISMAX(1,3).EQ.-1) THEN
            CALL GetUEGUmatEl(IDI,IDJ,IDK,IDL,ISS,G1,ALAT,iPeriodicDampingType,GetUMatEl)
         ENDIF
!         write (6,*) idi,idj,idk,idl,GetUMatEl
!         WRITE(6,"(4I5,$)") IDI,IDJ,IDK,IDL
!         WRITE(6,*) GETUMATEL,ABS(GETUMATEL)
         RETURN
!         CALL TIHALT(' GETUMATEL',ISUB)
      END function


!  TSMALL is used if we create a pre-freezing cache to hold just the <ij|kj> integrals
      SUBROUTINE SETUPUMATCACHE(NSTATE,TSMALL)
         IMPLICIT NONE
         INTEGER NSTATE
         LOGICAL TSMALL
         INTEGER ierr
         INCLUDE 'irat.inc'
         NTYPES=HElementSize
         NHITS=0
         NMISSES=0
         iCacheOvCount=0
         NSTATES=NSTATE
         IF(NSLOTSINIT.LE.0) THEN
            NSLOTS=0
            WRITE(6,*) "Not using UMATCACHE."
         ELSE
            NPAIRS=NSTATES*(NSTATES+1)/2
            IF(TSMALL) THEN
               NSLOTS=NSTATES
               tSmallUMat=.TRUE.
               WRITE(6,*) "Using small pre-freezing UMat Cache."
            ELSE
               IF(nMemInit.NE.0) THEN
                  WRITE(6,*) "Allocating ",nMemInit,"Mb for UMatCache+Labels."
                  nSlotsInit=(nMemInit*1048576/8)/(nPairs*(nTypes*HElementSize+1.D0/irat))
               ENDIF
               NSLOTS=MIN(NPAIRS, NSLOTSINIT)
               tSmallUMat=.FALSE.
            ENDIF
            UMATCACHEFLAG=0
            WRITE(6,"(A,I3,2I7,I10)") "UMAT NTYPES,NSLOTS,NPAIRS,TOT",NTYPES,NSLOTS,NPAIRS,NSLOTS*NPAIRS*NTYPES
            TUMAT2D=.FALSE.
            Allocate(UMatCacheData(0:nTypes-1,nSlots,nPairs), STAT=ierr)
            CALL MemAlloc(ierr,UMatCacheData,nTypes*HElementSize*nSlots*nPairs,'UMATCACHE')
            Allocate(UMatLabels(nSlots,nPairs), STAT=ierr)
            CALL MemAlloc(ierr,UMatLabels,nSlots*nPairs/irat+1,'UMATLABELS')
            CALL AZZERO(UMatCacheData,nTypes*HElementSize*nPairs*nSlots)
            CALL IAZZERO(UMATLABELS,nPairs*nSlots)
            if (.not.tSmallUMat.and.tReadInCache) then
                write (6,*) 'reading in cache'
                call ReadInUMatCache
            end if
         ENDIF
      END subroutine

      SUBROUTINE SETUPUMAT2D(G1,HarInt)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER ierr
         complex*16 HarInt(nStates,nStates)
         IF((NSLOTSINIT.LT.0).AND.(.not.TSTARSTORE)) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            Allocate(UMat2D(nStates,nStates),STAT=ierr)
            Call MemAlloc(ierr,UMat2D,HElementSize*NSTATES*NSTATES,'UMAT2D')
            CALL CPMDANTISYMINTEL(G1,UMAT2D,HarInt,NSTATES)
         ENDIF
      END subroutine


      SUBROUTINE SETUPUMAT2D_DF()
         IMPLICIT NONE
         INTEGER ierr
         IF(NSLOTSINIT.LT.0) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            Allocate(UMat2D(nStates,nStates),STAT=ierr)
            Call MemAlloc(ierr,UMat2D,HElementSize*NSTATES*NSTATES,'UMAT2D')
            IF(TSTARSTORE) THEN
                RETURN
            ELSE
                CALL ReadDalton2EIntegrals(nStates,UMat2D)
            ENDIF
         ENDIF
      END subroutine

     
      SUBROUTINE SETUMATTRANS(TRANS)
         IMPLICIT NONE
         INTEGER TRANS(NSTATES),ierr
         Allocate(TransTable(nStates),STAT=ierr)
         CALL MemAlloc(ierr,TransTable,NSTATES,'TransTable')
         CALL ICOPY(NSTATES,TRANS,1,TransTable,1)
         TTRANSGTID=.TRUE.
      END subroutine
 
      SUBROUTINE SetupUMatTransTable(OldNew,nOld,nNew)
         IMPLICIT NONE
         INTEGER nNew,nOld,I
         INTEGER OldNew(*),ierr
         LOGICAL tDiff
         Allocate(TransTable(nNew/2),STAT=ierr)
         CALL MemAlloc(ierr,TransTable,nNew/2,'TransTable')
         Allocate(InvTransTable(nOld/2),STAT=ierr)
         CALL MemAlloc(ierr,InvTransTable,nOld/2,'InvTransTable')
         CALL IAZZERO(InvTransTable, nOld/2)
         tDiff=.FALSE.
         DO I=2,nOld,2
            IF(OldNew(I).NE.0) THEN
               TransTable(OldNew(I)/2)=I/2
               InvTransTable(I/2)=OldNew(I)/2
               IF(OldNew(I)/2.NE.I/2) tDiff=.TRUE.
            ENDIF
         ENDDO
         IF(tDiff) THEN
            Write(6,*) "New->Old State Translation Table"
            DO I=1,nNew/2
               WRITE(6,*) I,TransTable(I)
            ENDDO
         ENDIF
         TTRANSFINDX=.TRUE.
      END subroutine
      SUBROUTINE DESTROYUMATCACHE
         IMPLICIT NONE
         CALL WriteUMatCacheStats()
         IF(ASSOCIated(UMatCacheData)) THEN
            WRITE(6,*) "Destroying UMatCache"
            CALL MemDealloc(UMatCacheData)
            Deallocate(UMatCacheData)
            CALL MemDealloc(UMATLABELS)
            Deallocate(UMatLabels)
            IF(ASSOCIated(UMat2D)) THEN
               CALL MemDealloc(UMat2D)
               Deallocate(UMat2D) 
            ENDIF
            IF(ASSOCIated(TransTable)) THEN
               CALL MemDealloc(TransTable)
               Deallocate(TransTable)
            ENDIF
            IF(ASSOCIated(InvTRANSTABLE)) THEN
               CALL MemDealloc(InvTransTable)
               Deallocate(InvTRANSTABLE)
            ENDIF
         ENDIF
      END subroutine

      SUBROUTINE WriteUMatCacheStats
         IMPLICIT NONE
         IF(ASSOCIated(UMatCacheData)) THEN
            WRITE(6,*) "UMAT Cache Statistics"
            WRITE(6,*) NHITS, " hits"
            WRITE(6,*) NMISSES, " misses"
            WRITE(6,*) iCacheOvCount, " overwrites"
            WRITE(6,"(F6.2,A)") (NHITS/(NHITS+NMISSES+0.D0))*100,"% success"
         ENDIF
      END SUBROUTINE
      SUBROUTINE SETUMATCACHEFLAG(NEWFLAG)
         IMPLICIT NONE
         INTEGER NEWFLAG,NF
         SELECT CASE(UMATCACHEFLAG)
         CASE(1)
!  We were in direct cache mode where values were distributed correctly throughout the cache.
            IF(NEWFLAG.EQ.0.AND..NOT.tSmallUMat) THEN
!  We need to fill the cache properly with values
               CALL FILLUPCACHE()
            ENDIF
         ENDSELECT
         UMATCACHEFLAG=NEWFLAG
!         WRITE(69,*) "FLAG",NEWFLAG
         SELECT CASE(NEWFLAG)
         CASE(1)
            IF(NSLOTS.EQ.NPAIRS) THEN ! we're storing every element, so we don't need to deal with different cacheing
               UMATCACHEFLAG=0
            ELSE
               CALL IAZZERO(UMATLABELS,NSLOTS*NPAIRS)
!Turn on the direct caching, and clear the cache.
            ENDIF
         ENDSELECT
         RETURN
      END subroutine

!The cache consists of an unordered set of labels and elements.
!We must order this, and then distribute the elements throughout each set of SLOTS.
      SUBROUTINE FillUpCache()
         IMPLICIT NONE
         INTEGER I,J,K,N,nK
         DO I=1,nPairs
!Find the last value in the cache
            CALL BinarySearch(nPairs+1,UMatLabels(1:nSlots,I),1,nSlots,N,J,K)
            N=J-1
!  N is now the last element and thus number of elements.
!  Sort according to label
            CALL SortIRN(N,UMatLabels(1,I),UMatCacheData(0,1,I),nTypes*HElementSize)
            K=nSlots
!  Now disperse among the whole array, from the end
            DO J=N,1,-1
               nK=(nSlots*(J-1))/N+1
               UMatLabels(nK:K,I)=UMatLabels(J,I)
               DO K=K,nK,-1
                  UMatCacheData(:,K,I)=UMatCacheData(:,J,I)
               ENDDO
               K=nK
            ENDDO
         ENDDO
      END subroutine

!   A binary search to find VAL in TAB.  TAB is sorted, but can have
!   multiple entries being the same.  If the search terminated unsuccessfully, 
!   the entry indicated is one after half-way through the set of entries which 
!   would be immediately prior to it.  From here until the label changes
!   should be filled with VAL if it is to be entered into the table.
!   A and B are the limits of the table.
      SUBROUTINE BINARYSEARCH(VAL,TAB,A,B,LOC,LOC1,LOC2)
         IMPLICIT NONE
         INTEGER VAL,A,B,LOC,LOC1,LOC2
         INTEGER TAB(A:B)
         INTEGER I,J,IFIRST,N,ILAST
!         DO I=A,B
!            WRITE(6,*) I,TAB(I)
!         ENDDO
         I=A
         J=B
         IFIRST=I
         ILAST=J
         DO WHILE(J-I.GE.1)
            N=(I+J)/2
!            WRITE(6,"(A,5I3)") "TN",I,J,N,TAB(N),VAL
            IF(TAB(N).LT.VAL.AND.TAB(N).NE.0.AND.I.NE.N) THEN
               IF(TAB(N).NE.TAB(IFIRST)) IFIRST=N
!   reset the lower limit
               I=N
            ELSEIF(TAB(N).GT.VAL.OR.TAB(N).EQ.0) THEN
               IF(TAB(N).NE.TAB(ILAST)) ILAST=N
!   reset the upper limit
               J=N
            ELSEIF(TAB(N).EQ.VAL) THEN
!   bingo, we've got it!
               LOC=N
!         DO I=A,B
!            WRITE(6,*) I,TAB(I),I.EQ.LOC
!         ENDDO
               LOC1=N
               LOC2=N
               RETURN
            ELSE
!   we've reached a situation where I and J's entries have the same value, and it's
!   not the one we want.  Leave the loop.
               I=J
            ENDIF
         ENDDO
!   We've failed.  However, the new value should sit between I and J.
!   Split whichever of the prior or after slots which has the most duplicates
!         WRITE(6,*) "FAIL:",IFIRST,I,J,ILAST
         LOC1=IFIRST+1
         LOC2=ILAST-1
         IF(TAB(IFIRST).EQ.TAB(ILAST)) THEN
            LOC=(IFIRST+ILAST)/2
            LOC1=IFIRST
            LOC2=ILAST
         ELSEIF(I-IFIRST.GE.ILAST-J) THEN
            LOC=(IFIRST+I)/2
         ELSE
            LOC=(ILAST+J)/2
         ENDIF
!         DO I=A,B
!            WRITE(6,*) I,TAB(I),I.EQ.LOC
!         ENDDO
      END subroutine

!   Get a unique index corresponding to pair (I,J), and return in RET.
!   I<=J<=N.  12/5/06
!  Old Scheme
!   e.g. 11 12 13 14 15 corresponds to    1  2  3  4  5
!           22 23 24 25                      6  7  8  9
!              33 34 35                        10 11 12
!                 44 45                           13 14
!                    55                              15
!  New Scheme 1/2/07
!   e.g. 11 12 13 14 15 corresponds to    1  2  4  7 11
!           22 23 24 25                      3  5  8 12
!              33 34 35                         6  9 13
!                 44 45                           10 14
!                    55                              15
      SUBROUTINE GETCACHEINDEX(I,J,N,RET)
         IMPLICIT NONE
         INTEGER I,J,RET,N
         RET=J*(J-1)/2+I
!         RET=N*(I-1)-I*(I-1)/2+J
      END subroutine


!      2n              (n)   sqrt(2n)
!     2  4  8 14 22... 92  -> 1  2  2  3  4 ... 9
!        6 10 16 24    94        2  3  4  4     9
!          12 18 26    96           3  4  5     9
!             20 28    98              4  5     9
!                30   100                 5    10
!                     102
!                     104
!                     106
!                     108
!                     110                      10
!   return I<=J
      SUBROUTINE GETCACHEINDEXSTATES(IND,N,I,J)
         IMPLICIT NONE
         INTEGER I,J,IND,N
         J=SQRT(2.0*IND)
         IF(J*(J+1)/2.LT.IND) J=J+1
         I=IND-J*(J-1)/2
      END subroutine
      SUBROUTINE SWAP(A,B)
         IMPLICIT NONE
         INTEGER A,B,C
         C=A
         A=B
         B=C
         RETURN
      END subroutine



!  We're in the middle of freezing some orbitals.
!  OrbTrans(i) will give us the new position of the old orbital i.
      SUBROUTINE FreezeUMatCache(OrbTrans,nOld,nNew)
         IMPLICIT NONE
         INTEGER nOld,nNew,OrbTrans(nOld)
         INTEGER onSlots,onPairs
         INTEGER I,J
         if(nNew/2.NE.nStates.OR.tSmallUMat) THEN
            WRITE(6,*) "Reordering UMatCache for freezing"
            onSlots=nSlots
            onPairs=nPairs
            CALL FreezeUMatCacheInt(OrbTrans,nOld,nNew,onSlots,onPairs)
         else
            WRITE(6,*) "UMatCache size not changing.  Not reordering."
         endif
      END subroutine

      SUBROUTINE FreezeUMAT2D(OldBasis,NewBasis,OrbTrans,iSS)
         IMPLICIT NONE
         INTEGER NewBasis,OldBasis,iSS,ierr,OrbTrans(OldBasis),i,j
         TYPE(HElement),POINTER :: NUMat2D(:,:)

         Allocate(NUMat2D(NewBasis/iSS,NewBasis/iSS),STAT=ierr)
         CALL MemAlloc(ierr,NUMat2D,HElementSize*(NewBasis/iSS)**2,'NUMat2D')
         NUMat2D(:,:)=HElement(0.D0)
         DO i=1,OldBasis/2
            IF(OrbTrans(i*2).NE.0) THEN
                DO j=1,OldBasis/2
                    IF(OrbTrans(j*2).NE.0) THEN
                        NUMat2D(OrbTrans(i*2)/2,OrbTrans(j*2)/2)=UMat2D(i,j)
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        CALL MemDealloc(UMat2D)
        Deallocate(UMat2D)
        UMat2D=>NUMat2D
        NULLIFY(NUMat2D)
        RETURN
      END subroutine
                
      SUBROUTINE FreezeUMatCacheInt(OrbTrans,nOld,nNew,onSlots,onPairs)
         IMPLICIT NONE
         INTEGER nOld,nNew,OrbTrans(nOld)
         TYPE(HElement),Pointer :: NUMat2D(:,:) !(nNew/2,nNew/2)
         TYPE(HElement) El(0:nTypes-1)
         INTEGER i,j,k,l,m,n
         INTEGER ni,nj,nk,nl,nm,nn,A,B,iType
         TYPE(HElement),Pointer :: OUMatCacheData(:,:,:) !(0:nTypes-1,onSlots,onPairs)
         INTEGER,Pointer :: OUMatLabels(:,:) !(onSlots,onPairs)
         
         INTEGER onSlots,onPairs,ierr
         LOGICAL toSmallUMat,tlog,toUMat2D
         LOGICAL GetCachedUMatEl
                  
         toUMat2D=tUMat2D
         IF(tUMat2D) then
            Allocate(NUMat2D(nNew/2,nNew/2),STAT=ierr)
            CALL MemAlloc(ierr,NUMat2D,HElementSize*(nNew/2)**2,'NUMat2D')
! /2 because UMat2D works in states, not in orbitals
            DO i=1,nOld/2
               IF(OrbTrans(i*2).NE.0) THEN
                  DO j=1,nOld/2
                     IF(OrbTrans(j*2).NE.0) THEN
                       NUMat2D(OrbTrans(i*2)/2,OrbTrans(j*2)/2)=UMat2D(i,j)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO    
            CALL MemDealloc(UMat2D)
            Deallocate(UMat2D)
            UMat2D=>NUMat2D
         endif
! Now go through the other cache.
! First save the memory used for it.
!         onSlots=nSlots
!         onPairs=nPairs
         OUMatCacheData=>UMatCacheData
         OUMatLabels=>UMatLabels
         toSmallUMat=tSmallUMat
         Nullify(UMatCacheData)
         Nullify(UMatLabels)
!Now reinitialize the cache.
         CALL SetupUMatCache(nNew/2,.FALSE.)
         TUMAT2D=toUMat2D
         CALL SetUMatcacheFlag(1)
         DO i=1,nOld/2
          IF(OrbTrans(i*2).NE.0) THEN
           DO k=i,nOld/2
            IF(OrbTrans(k*2).NE.0) THEN
             CALL GetCacheIndex(i,k,nOld/2,m)
             DO n=1,onSlots
              IF((onSlots.EQ.onPairs.OR.toSmallUMat)                    &
     &            .OR.((onSlots.NE.onPairs).AND.((n.EQ.1).OR.           &
     &            OUMatLabels(n,m).NE.OUMatLabels(n-1,m)))) THEN
               IF(OUMatLabels(n,m).NE.0) THEN
                ni=OrbTrans(i*2)/2
                nk=OrbTrans(k*2)/2
!Now get the label of the slot and convert to orbitals
                IF(onSlots.EQ.onPairs) THEN
                 CALL GetCacheIndexStates(n,nOld/2,j,l)
                ELSEIF(toSmallUMat) THEN
                  j=n
                  l=n
                ELSE
                 CALL GetCacheIndexStates(OUMatLabels(n,m),nOld/2,j,l)
                ENDIF
                nj=OrbTrans(j*2)/2
                nl=OrbTrans(l*2)/2
                IF(nj.NE.0.AND.nl.NE.0) THEN
!                  WRITE(6,"(A,4I5,A,2I5,$)") "NC",I,J,K,L,"L",m,n
                   Tlog=GetCachedUmatEl(ni,nj,nk,nl,El,nm,nn,A,B,iType)
!                   WRITE(6,"(A,4I5,2I5,$)") "->",NI,NJ,NK,NL,A,B
!OUMatLabels(n,m)
!                   WRITE(6,*) OUMatCacheData(0,n,m)
                   CALL CacheUMatEl(A,B,OUMatCacheData(0,n,m),nm,nn,iType)
                ENDIF
               ENDIF
              ENDIF
             ENDDO
            ENDIF
           ENDDO
          ENDIF
         ENDDO
         CALL MemDealloc(OUMatLabels) 
         Deallocate(OUMatLabels)
         CALL MemDealloc(OUMatCacheData) 
         Deallocate(OUMatCacheData)
         CALL SetUMatCacheFlag(0)               
      END subroutine


! JSS: Read in cache file.
      subroutine ReadInUMatCache()
      implicit none
      integer  i,j,k,l,iCache1,iCache2,A,B,readerr,iType
      integer  iSlot,iPair
      type(HElement) UMatEl(0:nTypes-1),DummyUMatEl(0:nTypes-1)
      logical  tDummy,GetCachedUMatEl,testfile
      inquire(file="CacheDump",exist=testfile)
      if (.not.testfile) then
          write (6,*) 'CacheDump does not exist.'
          return
      end if
      open (21,file="CacheDump",status="old",iostat=readerr)
      if (readerr.ne.0) then 
          write (6,*) 'Error reading CacheDump.'
          return
      end if
      read (21,*) nStatesDump
      readerr=0
      do while (readerr.eq.0)
        read (21,*,iostat=readerr) i,j,k,l,UMatEl
        DummyUMatEl=UMatEl
        if (TTRANSFINDX) then
            i=TransTable(i)
            j=TransTable(j)
            k=TransTable(k)
            l=TransTable(l)
        end if
        if (min(i,j,k,l).gt.0.and.max(i,j,k,l).le.nStates) then
            tDummy=GetCachedUMatEl(i,j,k,l,DummyUmatEl,iCache1,iCache2,A,B,iType)
            call CacheUMatEl(A,B,UMatEl,iCache1,iCache2,iType)
        end if
      end do
      close(21)
      return
      end subroutine

!  JSS: Print out the cache contents so they can be read back in for a future
!  calculation.  Need to print out the full set of indices, as the number of
!  states may change with the next calculation.
      subroutine DumpUMatCache(NHG,G1)
      implicit none
      include  'basis.inc'
      integer  NHG
      type(BasisFN) G1(NHG)
      ! Variables
      integer iPair,iSlot,i,j,k,l,iCache1,iCache2,A,B,iType
      logical GetCachedUMatEl,LSymSym
      integer*8 TotSymRep
      type(HElement) UMatEl(0:nTypes-1)
      type(Symmetry) Sym,Symprod,SymConj
      ! 1. test read in.
      ! 2. binary file.
      open (21,file="CacheDump",status="unknown")
      write (21,*) nStates
      do iPair=1,nPairs
        do iSlot=iPair,nSlots
          call GetCacheIndexStates(iPair,nStates,i,k)
          call GetCacheIndexStates(iSlot,nStates,j,l)
          Sym%s=TotSymRep()
!          Sym=SymProd(Sym,SymConj(G1(I*2-1)%Sym))
!          Sym=SymProd(Sym,SymConj(G1(J*2-1)%Sym))
!          Sym=SymProd(Sym,G1(K*2-1)%Sym)
!          Sym=SymProd(Sym,G1(L*2-1)%Sym)
          if (LSymSym(Sym)) then
              if (.not.GetCachedUMatEl(i,j,k,l,UMatEl,iCache1,iCache2,A,B,iType)) then
                  if (TTRANSFINDX) then
                      i=InvTransTable(i)
                      j=InvTransTable(j)
                      k=InvTransTable(k)
                      l=InvTransTable(l)
                  end if
                  write (21,*) i,j,k,l,UMatCacheData(:,ICACHE2,ICACHE1)
              end if
          end if
        end do
      end do
      close(21,status="keep")
      return
      end subroutine DumpUMatCache
      
END MODULE UMatCache

      logical function HasKPoints()
         IMPLICIT NONE
         include 'cpmddata.inc'
         IF(NKPS.GT.1) THEN
            HasKPoints=.TRUE.
         ELSE
            HasKPoints=.FALSE.
         ENDIF
      end function
 

      SUBROUTINE GTID(NBASISMAX,GIND,ID)
         USE UMatCache
         IMPLICIT NONE
         INTEGER GIND,NBASISMAX(5,3),ID
            IF(NBASISMAX(2,3).GT.0) THEN
               ID=(GIND-1)/NBASISMAX(2,3)+1
            ELSE
               ID=(GIND-1)/2+1
               IF(TTRANSGTID) ID=TRANSTABLE(ID)
            ENDIF
      RETURN
      END subroutine

!  CacheUMatEl and GetCachedUMatEl are outside the module so that they can interact with CPMD
!      REAL*8  FUNCTION WR
!         WRITE(70,*) UMatCacheData(0,2,7)
!         CALL FLUSH(70)
!         WR=sq(UMatCacheData(0,2,7))
!      END
!     Set an element in the cache.  All the work has been done for us before
!   as the element we have to set is in (ICACHEI,ICACHE)
!    iType tells us whether we need to swap/conjugate the nTypes integrals within the slot
!   We still need to fill out the space before or after  us if we've been put in the 
!   middle of a block of duplicates
      SUBROUTINE CACHEUMATEL(A,B,UMATEL,ICACHE,ICACHEI,iType)
         USE HElem
         USE UMatCache
         IMPLICIT NONE
         INTEGER A,B,ICACHE,ICACHEI
         TYPE(HElement) UMATEL(0:NTYPES-1),TMP(0:NTYPES-1)
         INTEGER OLAB,IC1,I,J,ITOTAL
         INTEGER iType
         INTEGER iIntPos
         SAVE ITOTAL
         DATA ITOTAL /0/
         if (nSlots.eq.0) return
!         WRITE(6,*) "CU",A,B,UMATEL,iType
         if(nTypes.gt.1) then
! A number of different cases to deal with depending on the order the integral came in (see GetCachedUMatEl for details)
!  First get which pos in the slot will be the new first pos
            iIntPos=iand(iType,1) 
!  If bit 1 is set we must conjg the (to-be-)first integral
            if(btest(iType,1)) then
               Tmp(0)=dconjg(UMatEl(iIntPos))
            else
               Tmp(0)=UMatEl(iIntPos)
            endif
!  If bit 2 is set we must conjg the (to-be-)second integral
            if(btest(iType,2)) then
               Tmp(1)=dconjg(UMatEl(1-iIntPos))
            else
               Tmp(1)=UMatEl(1-iIntPos)
            endif
            UMatEl=Tmp
         endif
!         WRITE(6,*) "CU",A,B,UMATEL,iType
!         WRITE(69,*) NSLOTS,A,B,UMATEL,ICACHE,ICACHEI
         IF(NSLOTS.EQ.NPAIRS.OR.UMATCACHEFLAG.EQ.1.OR.tSmallUMat) THEN
!   small system.  only store a single element
            UMATLABELS(ICACHEI,ICACHE)=B
            UMatCacheData(:,ICACHEI,ICACHE)=UMATEL
            ITOTAL=ITOTAL+1
            RETURN
         ENDIF
         IC1=ICACHEI
!         WRITE(6,*) "ICI",ICACHEI,ICACHE
         OLAB=UMATLABELS(ICACHEI,ICACHE)
!   If we're in a block of prior, fill after
         DO WHILE(OLAB.LT.B.AND.ICACHEI.LE.NSLOTS)
            UMatCacheData(:,ICACHEI,ICACHE)=UMATEL
            UMATLABELS(ICACHEI,ICACHE)=B
!            IF(ICACHEI.LT.1.OR.ICACHE.LT.1.OR.ICACHEI.GT.NSLOTS.OR.ICACHE.GT.NPAIRS) THEN
!               WRITE(6,*) ICACHEI,ICACHE
!               STOP "a"
!            ENDIF
            ICACHEI=ICACHEI+1
            IF(ICACHEI.LE.NSLOTS) THEN
               OLAB=UMATLABELS(ICACHEI,ICACHE)
            ELSE
               OLAB=0
            ENDIF
         ENDDO
         IF(OLAB.EQ.0) ICACHEI=IC1
!        WRITE(6,*) "ICI2",ICACHEI,ICACHE
         DO WHILE((OLAB.GT.B.OR.OLAB.EQ.0).AND.ICACHEI.GT.0)
            UMatCacheData(:,ICACHEI,ICACHE)=UMATEL
            UMATLABELS(ICACHEI,ICACHE)=B
!            IF(ICACHEI.LT.1.OR.ICACHE.LT.1.OR.ICACHEI.GT.NSLOTS.OR.ICACHE.GT.NPAIRS) THEN
!               WRITE(6,*) ICACHEI,ICACHE
!               STOP "b"
!            ENDIF
            ICACHEI=ICACHEI-1
            OLAB=UMATLABELS(ICACHEI,ICACHE)
         ENDDO
      END subroutine
!   Lookup in the cache to see if there's a stored element.  If not, return TRUE.
!   If there is, return the stored element in UMatEl.
!    This will rearrange IDI,IDJ,IDK,IDL into the correctkorder
!   (i,k)<=(j,l) and i<=k, j<=l.  ICACHE corresponds to the pair (i,j), and
!   ICACHEI is the index in that cache where the cache should be located.
      LOGICAL FUNCTION GETCACHEDUMATEL(IDI,IDJ,IDK,IDL,UMATEL,ICACHE,ICACHEI,A,B,ITYPE)
         USE HElem
         USE UMatCache
         IMPLICIT NONE
         INTEGER IDI,IDJ,IDK,IDL,ICACHE,ICACHEI
         INTEGER ICACHEI1,ICACHEI2
         TYPE(HElement) UMATEL
         INTEGER I,A,B,ITYPE,ISTAR,ISWAP
!         WRITE(6,"(A,4I5)") "GCUI",IDI,IDJ,IDK,IDL
         IF(NSLOTS.EQ.0) THEN
!We don't have a cache so signal failure.
            GETCACHEDUMATEL=.TRUE.
            ICACHE=0
            ITYPE=0
            RETURN
         ENDIF
!   First ensure the indices are in the correct order
         IF(tSmallUMat) THEN
!tSmallUMat is set if we have nStates slots per pair for storing the <ik|jk> integrals.
            ITYPE=0
            IF(IDI.EQ.IDK) THEN
               B=IDI
               IF(IDL.LT.IDJ) THEN
                  CALL SWAP(IDJ,IDL)
                  ITYPE=2
               ENDIF
               CALL GETCACHEINDEX(IDJ,IDL,NSTATES,A)
            ELSEIF(IDJ.EQ.IDL) THEN
               B=IDJ
               IF(IDK.LT.IDI) THEN
                  CALL SWAP(IDI,IDK)
                  ITYPE=2
               ENDIF
               CALL GETCACHEINDEX(IDI,IDK,NSTATES,A)
            ELSE !Can't find a pair the same
!We don't have a cache for this type of element so signal failure.
               GETCACHEDUMATEL=.TRUE.
               ICACHE=0
               ITYPE=0
               RETURN
            ENDIF
         ELSE
!Otherwise Normal caching
            ITYPE=0
            ISTAR=0
            ISWAP=0
            IF(IDK.LT.IDI) THEN
               CALL SWAP(IDI,IDK)
               ISTAR=IOR(ISTAR,1)
            ENDIF
            IF(IDL.LT.IDJ) THEN
               CALL SWAP(IDJ,IDL)
               ISTAR=IOR(ISTAR,2)
            ENDIF
            CALL GETCACHEINDEX(IDI,IDK,NSTATES,A)
            CALL GETCACHEINDEX(IDJ,IDL,NSTATES,B)
            IF(A.GT.B) THEN
               CALL SWAP(A,B)
               CALL SWAP(IDI,IDJ)
               CALL SWAP(IDK,IDL)
               ISWAP=1
            ENDIF
            IF(HElementSize.EQ.1) THEN
!  Eight integrals from ijkl are the same.
               ITYPE=0
            ELSE
!  Complex orbitals and integrals, so we need to consider different types
!  Using notation abcd rather than ijkl.  6/2/07 and 19/2/06
!
!  <> mean swap pairs <ab|cd> -> <ba|dc>
!  *.  means complex conjugate of first codensity i.e. <ab|cd> -> <ca|bd>
!  .* for second and ** for both.
!
!  abcd   -> badc <>
!  |  |
!  | \|/       |-> cdab ** -> dcba **<>
!  | cbad *.  -|
!  |           |-> bcda *.<>
! \|/
!  adcb .* -> dabc .*<>

!Now consider what must occur to the other integral in the slot to recover pair (abcd,cbad).  0 indicates 1st in slot, 1 means 2nd.  * indicated conjg.
!
!  ..    abcd  cbad  0  1
!  *.    cbad  abcd  1  0
!  .*    adcb  cbad  1* 0*
!  **    cdab  adcb  0* 1*
!    <>  badc  dabc  0  1*
!  *.<>  bcda  dcba  1* 0
!  .*<>  dabc  badc  1  0*
!  **<>  dcba  bcda  0* 1


! Of the type, bit zero indicates which of the two integrals in a slot to use.  Bit 1 is set if the integral should be complex conjugated.
!   Bit 2 is set if the other integral in the slot should be complex conjugated if we are to have the structure (<ij|kl>,<kj|il>) in the slot.
!      This is only used by CacheUMatEl when adding a slot to the cache.
!
!  <ij|u|kl> and <kj|u|il> (which are distinct when complex orbitals are used).
!  TYPE 0          TYPE 1
!
!
               IF(ISTAR.EQ.0) THEN
                  IF(ISWAP.EQ.0) THEN
                     ITYPE=0  !0  1
                  ELSE
                     ITYPE=4  !0  1*
                  ENDIF
               ELSEIF(ISTAR.EQ.1) THEN
!  If we star the first pair, that corresponds to the plain TYPE 1.  If we swap too, then we complex conj.
                  IF(ISWAP.EQ.0) THEN
                     ITYPE=1  !1  0
                  ELSE
                     ITYPE=3  !1* 0
                  ENDIF
               ELSEIF(ISTAR.EQ.2) THEN
!  If we star the second pair, that corresponds to TYPE 1.
!  If there's no swap, it's complex conjugated, otherwise it's not.
                  IF(ISWAP.EQ.0) THEN
                     ITYPE=7  !1* 0*
                  ELSE
                     ITYPE=1  !1  0*
                  ENDIF
               ELSEIF(ISTAR.EQ.3) THEN
! We've starred both pairs
!  We complex conjg setting bit 1 but using type 0
                  IF(ISWAP.EQ.0) THEN
                     ITYPE=6  !0* 1*
                  ELSE
                     ITYPE=2  !0* 1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!         WRITE(6,"(A,7I5)") "GCUE",IDI,IDJ,IDK,IDL,A,B,iType
         ICACHE=A
!               IF(ISTAR.EQ.1) THEN
!!  If we star the first pair, that corresponds to the plain TYPE 1
!                  ITYPE=1
!               ELSEIF(ISTAR.EQ.2) THEN
!!  If we star the second pair, that corresponds to TYPE 1.
!!  If there's no swap, it's complex conjugated, otherwise it's not.
!                  IF(ISWAP.EQ.0) THEN
!                     ITYPE=3
!                  ELSE
!                     ITYPE=1
!                  ENDIF
!               ELSEIF(ISTAR.EQ.3) THEN
!! We've starred both pairs
!!  We complex conjg setting bit 1 but using type 0
!                  ITYPE=2
!               ENDIF
!            ENDIF
!         ENDIF
!         ICACHE=A
         IF(NSLOTS.EQ.NPAIRS.OR.tSmallUMat) THEN
!   we've a small enough system to store everything.
            ICACHEI=B
         ELSE
            IF(UMATCACHEFLAG.EQ.1) THEN
!  UMATCACHEFLAG=1 means we are storing a sequence of cache elements,
!  in a blank cache
!  We store them linearly in the cache, and distribute them around later

!Find the last value in the cache
            CALL BINARYSEARCH(NPAIRS+1,UMATLABELS(1:NSLOTS,A),1,NSLOTS,ICACHEI,ICACHEI1,ICACHEI2)
               ICACHEI=ICACHEI1
               ICACHEI2=ICACHEI1
               IF(UMatLabels(iCacheI,A).NE.0) iCacheOvCount=iCacheOvCount+1
!WRITE(6,*) "Cache Overwrite", A,B
   
!                  WRITE(6,*) IDI,IDJ,IDK,IDL
!                  WRITE(6,*) A,B,NSLOTS,NPAIRS
!                  WRITE(6,*) ICACHEI1,ICACHEI2
!                  WRITE(6,*) ICACHEI
            ELSE
            CALL BINARYSEARCH(B,UMATLABELS(1:NSLOTS,A),1,NSLOTS,ICACHEI,ICACHEI1,ICACHEI2)
            ENDIF
         ENDIF
         IF(UMATLABELS(ICACHEI,ICACHE).EQ.B) THEN
            !WRITE(6,*) "C",IDI,IDJ,IDK,IDL,ITYPE,UMatCacheData(0:nTypes-1,ICACHEI,ICACHE)
            UMATEL=UMatCacheData(IAND(ITYPE,1),ICACHEI,ICACHE)
            IF(BTEST(ITYPE,1)) UMATEL=DCONJG(UMATEL)  ! Bit 1 tells us whether we need to complex conjg the integral
!   signal success
            GETCACHEDUMATEL=.FALSE.
         ELSE
!   signal failure
            GETCACHEDUMATEL=.TRUE.
!            WRITE(68,*) A,B,ICACHEI1,ICACHEI2,ICACHEI
         ENDIF
         RETURN
      END function

      
