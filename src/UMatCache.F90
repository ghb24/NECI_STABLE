
MODULE UMatCache

    use constants, only: dp,sizeof_int,int64

    use SystemData, only: tROHF,tStoreSpinOrbs, tComplexWalkers_RealInts, &
                          Symmetry, BasisFN, UMatEps, tROHF

    use SystemData, only: tRIIntegrals,tCacheFCIDUMPInts, t_non_hermitian

    use util_mod, only: swap, get_free_unit, NECI_ICOPY, near_zero

    use sort_mod

    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc

    use HElem, only: HElement_t_size

    use CPMDData, only: NKPS

    use sym_mod, only: TotSymRep, LSymSym

    use HElem, only: HElement_t_size

    use legacy_data, only: irat

    use procedure_pointers, only: get_umat_el

      IMPLICIT NONE

      SAVE

! Integrals are cached for CPMD and density fitting calculations, where we might
! not be able to store all possible integrals in the available memory.

! The integral cache is stored in UMatCacheData.
! For real systems, nTypes=1, and we have a single value for each unique ordered pair of ordered pairs (i,k), (j,l)
! out of <ij|u|kl> (where i,j,k,l are state labels).
! For complex systems, nTypes=2, giving two possible values.
! For each nPairs (which corresponds to the greater of the two pairs), there  are nSlots, with the other pair label in
! UMatLabels, and the values in UMatCache.
! nStates is the maximum number of states stored in the cache (which may not be all the states if there are frozen virtuals).

      HElement_t(dp), Pointer :: UMatCacheData(:,:,:) => null() !(0:nTypes-1,nSlots,nPairs)
      INTEGER, Pointer :: UMatLabels(:,:) => null() !(nSlots,nPairs)
      INTEGER :: nSlots, nPairs, nTypes
      INTEGER :: nStates

! tSmallUMat is set if we have nStates slots per pair for storing the <ik|jk> integrals.
! This should only be used prior to freezing to store precalculated integrals.
! For each pair (i,j), we store the <ik|jk> integral in slot k.
      LOGICAL tSmallUMat

!  iDumpCacheFlag: Dump the cache to disk if we're told to.
!    =0: No cache dumping.
!    =1: Dump cache unless we've read in a cache dump that's from a larger
!    calculation than the current one.
!    =2: Force dumping of cache (over-writing previous cache level).
!  nStatesDump: number of states used in calculation which produced dump file.
!  tReadInCache: does what it says on the tin.
      integer iDumpCacheFlag,nStatesDump
      logical tReadInCache

! For the more frequently used <ij|u|ij> and <ij|u|ji> integrals, we store them
! in a separate cache (if TUMat2D is true).
!     UMAT2D : Stores the integrals of the form <ij|ij> and <ij|ji>
!     <ij|ij> is stored in the upper diagaonal, <ij|ji> in the
!     off-diagonal elements of the lower triangle.
      HElement_t(dp), Pointer :: UMat2D(:,:) => null() !(nStates,nStates)
      LOGICAL :: tUMat2D, tDeferred_Umat2d

! This vector stores the energy ordering for each spatial orbital, which is the inverse of the BRR vector
! This is needed for the memory saving star indexing system.
! E.g. Element 4 will give the the order in the energy of element 4
     INTEGER, DIMENSION(:), POINTER :: INVBRR => null()
     INTEGER, DIMENSION(:), POINTER :: INVBRR2 => null()

!NOCC is number of occupied spatial orbitals - needed for test in UMATInd, thought would be quicker
!than passing it in each time.
!Freezetransfer is a temporary measure to tell UMATIND when the freezing of orbitals is occuring.
      INTEGER :: NOCC
      LOGICAL :: FREEZETRANSFER


! Book-keeping information
! nSlotsInit is the number of slots requested on input.  If the number required is less,
!then the lower value is allocated
! If nSlotsInit is set to 0, then general <ij|u|kl> element caching is not performed, but
!UMat2D <ij|u|ij> and <ij|u|ji> is.  For nSlotsInit=-1 neither is performed.
      INTEGER nSlotsInit,nMemInit

! UMatCacheFlag=0 is normal operation
! UMatCacheFlag=1 means cache from bottom up
! This is useful when lots of sequential pieces of data are being stored.
! When UMatCacheFlag is reset to 0, the data which are present are spread evenly around the slots for a given Pair.
      INTEGER UMatCacheFlag

! If the max vertex level is 2 or less, then we just need to calculate
! <ij|ab> and never need <ib|aj>, unless the integral is for a single
! excitation.
      LOGICAL gen2CPMDInts

! nHits and nMisses are the number of cache hits and misses.
      INTEGER :: nHits, nMisses
! The number of cache overwrites
      INTEGER :: iCacheOvCount

! Some various translation tables to convert between different orderings of states.
      LOGICAL :: tTransGTID, tTransFindx
      INTEGER, Pointer :: TransTable(:) => null() !(NSTATES)
      INTEGER, Pointer :: InvTransTable(:) => null() !(NSTATES)

! Density fitting cache information: for generating integrals on the fly from density fitting.
      integer nAuxBasis,nBasisPairs
      logical tDFInts
      real(dp),Pointer :: DFCoeffs(:,:) => null() !(nAuxBasis,nBasisPairs)
      real(dp),Pointer :: DFInts(:,:) => null() !(nAuxBasis,nBasisPairs)
      real(dp),Pointer :: DFFitInts(:,:) => null() !(nAuxBasis,nAuxBasis)
      real(dp),Pointer :: DFInvFitInts(:,:) => null() !(nAuxBasis,nAuxBasis)
      INTEGER iDFMethod
!Some possible DFMethods sums over P, Q implied.  All precontracted to run in order(X) except DFOVERLAP2NDORD
! 0 - no DF
! DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)
! DFOVERLAP2NDORD  2 - (ij|u|ab)= (ij|u|P)(P|ab)+(ij|P)(P|u|ab)-(ij|P)(P|u|Q)(Q|ab)
! DFOVERLAP2       3 - (ij|u|ab)= (ij|P)(P|u|Q)(Q|ab)
! DFCOULOMB        4 - (ij|u|ab)= (ij|u|P)[(P|u|Q)^-1](Q|u|ab)

      ! Memory book-keeping tags
      integer(TagIntType) :: tagUMatCacheData=0
      integer(TagIntType) :: tagUMatLabels=0
      integer(TagIntType) :: tagOUMatCacheData=0
      integer(TagIntType) :: tagOUMatLabels=0
      integer(TagIntType) :: tagUMat2D=0
      ! [W.D]
      ! those two tags are also defined in OneEInts..
!       integer(TagIntType) :: tagTMat2D=0
!       integer(TagIntType) :: tagTMat2D2=0
      integer(TagIntType) :: tagTransTable=0
      integer(TagIntType) :: tagInvTransTable=0
      integer(TagIntType) :: tagDFCoeffs=0
      integer(TagIntType) :: tagDFInts=0
      integer(TagIntType) :: tagDFFitInts=0
      integer(TagIntType) :: tagDFInvFitInts=0
      integer(TagIntType) :: tagInvBRR=0
      integer(TagIntType) :: tagInvBRR2=0

      Contains

      SUBROUTINE CreateInvBRR2(BRR2,NBASIS)
      ! Create new INVBRR for the freezing process
      ! In:
      !    BRR(i)=j: orbital i is the j-th lowest in energy.
      !    nBasis: size of bais
      ! InvBRR is the inverse of BRR.  InvBRR(j)=i: the j-th lowest energy
      ! orbital corresponds to the i-th orbital in the original basis.
        IMPLICIT NONE
        INTEGER NBASIS
        INTEGER BRR2(NBASIS),ierr,I,t
        character(*), parameter :: t_r='CreateInvBRR2'

!        WRITE(6,*) "================================"
!        WRITE(6,*) "BRR2 is "
!        WRITE(6,*) BRR2(:)

        ALLOCATE(INVBRR2(NBASIS/2),STAT=ierr)
        CALL LogMemAlloc('INVBRR2',NBASIS/2,4,t_r,tagINVBRR2,ierr)
        INVBRR2(1:NBASIS/2)=0
        t=0
        DO I=2,NBASIS,2
            t=t+1
            INVBRR2(BRR2(I)/2)=t
        ENDDO

!        WRITE(6,*) "================================"
!        WRITE(6,*) "InvBRR2 is "
!        WRITE(6,*) INVBRR2(:)

        RETURN
      END SUBROUTINE CreateInvBRR2


      SUBROUTINE CreateInvBRR(BRR,NBASIS)
      ! Create new INVBRR for the freezing process
      ! In:
      !    BRR(i)=j: orbital i is the j-th lowest in energy.
      !    nBasis: size of bais
      ! InvBRR is the inverse of BRR.  InvBRR(j)=i: the j-th lowest energy
      ! orbital corresponds to the i-th orbital in the original basis.
        IMPLICIT NONE
        INTEGER NBASIS
        INTEGER BRR(NBASIS),ierr,I,t
        character(*), parameter :: t_r='CreateInvBRR'

        IF(ASSOCIATED(INVBRR)) THEN
            CALL LogMemDealloc(t_r,tagINVBRR)
            DEALLOCATE(INVBRR)
        ENDIF
        ALLOCATE(INVBRR(NBASIS/2),STAT=ierr)
        CALL LogMemAlloc('INVBRR',NBASIS/2,4,t_r,tagINVBRR,ierr)
        INVBRR(1:NBASIS/2)=0
        t=0
        DO I=2,NBASIS,2
            t=t+1
            INVBRR(BRR(I)/2)=t
        ENDDO
        RETURN
      END SUBROUTINE CreateInvBRR



      INTEGER(int64) FUNCTION UMatInd(I,J,K,L)
         ! Get the index of physical order UMAT element <IJ|KL>.
         ! Indices are internally reordered such that I>K, J>L,(I,K)>(J,L)
         ! Note: (i,k)>(j,l) := (k>l) || ((k==l)&&(i>j))
         ! In:
         !    I,J,K,L: orbital indices. These refer to spin orbitals in
         !      unrestricted calculations and spatial orbitals in restricted
         !      calculations.
         !    nBasis: size of basis. If =0, use nStates instead.
         !    nOccupied: # of occupied orbitals.  If =0, then nOcc is used.
         !    Should only be passed as non-zero during the freezing process.
         use SystemData, only: nbasis
         IMPLICIT NONE
         INTEGER, intent(in) :: I,J,K,L
         INTEGER A,B, nbi, iss

         if (t_non_hermitian) then
             IF(tStoreSpinOrbs) THEN
                 iSS=1
             ELSE
                 iSS=2
             ENDIF

             nBi=nBasis/iSS

             A=(I-1)*nBi+K
             B=(J-1)*nBi+L
         else

             !Combine indices I and K, ensuring I>K
             IF(I.GT.K) THEN
                 A=(I*(I-1))/2+K
             ELSE
                 A=(K*(K-1))/2+I
             ENDIF

             !Combine indices J and L, ensuring J>L
             IF(J.GT.L) THEN
                 B=(J*(J-1))/2+L
             ELSE
                 B=(L*(L-1))/2+J
             ENDIF
         end if

         !Combine (IK) and (JL) in a unique way  (k > l or if k = l then i > j)
         IF(A.GT.B) THEN
             UMatInd=(int(A,int64)*int(A-1,int64))/2+int(B,int64)
         ELSE
             UMatInd=(int(B,int64)*int(B-1,int64))/2+int(A,int64)
         ENDIF
#ifdef __CMPLX
         if(.not. tComplexWalkers_RealInts) then
             UMatInd = (UmatInd-1)*2 + 1
             !We need to test whether we have swapped i and k or j and l independantly of each other
             !If we have done this, it is one of the 'other' integrals - add one.
             if (((I.gt.K).and.(J.lt.L)) .or. ((I.lt.K).and.(J.gt.L))) then
                UMatInd = UMatInd + 1
             endif
         endif
#endif
      END FUNCTION UMatInd

      HElement_t(dp) function UMatConj(I,J,K,L,val)
         integer, intent(in) :: I,J,K,L
         HElement_t(dp), intent(in) :: val
#ifdef __CMPLX
         INTEGER :: IDI,IDJ,IDK,IDL,NewA,A

         !Changing index ordering for the real ordering.
         IDI=I
         IDJ=J
         IDK=K
         IDL=L

         !First find rearranged indices.
         IF(idi.lt.idk) then
             !swap idi and idk
             call swap(idi,idk)
         ENDIF

         IF(idj.lt.idl) then
             !swap idj and idl
             call swap(idj,idl)
         ENDIF

         IF((idl.lt.idk).or.((idl.eq.idk).and.(idi.lt.idj))) THEN
             !We would want to swap the (ik) and (jl) pair.
             call swap(idi,idj)
             call swap(idk,idl)
         ENDIF

         !Indices now permuted to the real case ordering. Is this now the same integral?
         if (((I.gt.K).and.(J.lt.L)) .or. ((I.lt.K).and.(J.gt.L))) then
             !Type II integral - reducing to lowest ordering will give 'other'
             !integral, where one of (ik) and (jl) pairs have been swapped independantly.
             !If i = k, or j = l, we do not go through here.
             call swap(idi,idk)  !Revert back to the correct integral by swapping just i and k.
         endif

         !Want to see if the pairs of indices have swapped sides.
         !Make unique index from the ij and kl pairs
         IF(IDI.gt.IDJ) THEN
             A=IDI*(IDI-1)/2+IDJ
         ELSE
             A=IDJ*(IDJ-1)/2+IDI
         ENDIF

         !Create uniques indices from the original pairs of indices.
         !We only need to consider whether the (ij) pair has swapped sides, since
         !the <ij|ij> and <ij|ji> integrals are real by construction, and so we do not
         !need to consider what happens if the (ij) pair = (kl) pair.
         IF(I.gt.J) THEN
             NewA=I*(I-1)/2+J
         ELSE
             NewA=J*(J-1)/2+I
         ENDIF

         !Check whether pairs of indices have swapped sides.
         IF(NewA.ne.A) THEN
             UMatConj=CONJG(val) !Index pair i and j have swapped sides - take CC.
         ELSE
             UMatConj=val
         ENDIF

#else
        integer :: tmp
         UMatConj = val

         ! Eliminate warnings
         tmp=i; tmp=j; tmp=k; tmp=l
#endif
      end function UMatConj


      SUBROUTINE GetUMatSize(nBasis,iSize)
        use SystemData, only: tStoreSpinOrbs
      ! Get the prospective size of a UMat (not a UMatCache) for completely
      ! storing FCIDUMP 2-e integrals
      ! In:
      !    nBasis: as above.
      !    nEl: # electrons.
      !    iSS: ratio of spatial orbitals to spin orbitals.
      !         iSS=0 integrals not stored in UMAT.
      !         iSS=1 unrestricted calculation
      !         iSS=2 restricted calculation
      !         iSS=-1 flag for the Hubbard model.
      ! Out:
      !    iSize: size of UMAT.
         IMPLICIT NONE
         INTEGER nBasis,iSS
         INTEGER iPairs,nBi
         INTEGER(int64), intent(out) :: iSize
         IF(tStoreSpinOrbs) THEN
             iSS=1
         ELSE
             iSS=2
         ENDIF

         nBi=nBasis/iSS
         if (t_non_hermitian) then
             iPairs = nbi**2
         else
             iPairs=(nBi*(nBi+1))/2
         end if
         iSize=(int(iPairs,int64)*int(iPairs+1,int64))/2
#ifdef __CMPLX
         !Since we now only have 4-fold symmetry, rather than 8-fold.
         iSize = iSize * 2
#endif
      END SUBROUTINE GetUMatSize



      SUBROUTINE SETUPUMATCACHE(NSTATE,TSMALL)
         ! nState: # states.
         ! TSMALL is used if we create a pre-freezing cache to hold just the <ij|kj> integrals.
         IMPLICIT NONE
         INTEGER NSTATE
         real(dp) Memory
         LOGICAL TSMALL
         INTEGER ierr
         character(len=*),parameter :: thisroutine='SETUPUMATCACHE'
         NTYPES=HElement_t_size
         NHITS=0
         NMISSES=0
         iCacheOvCount=0
         NSTATES=NSTATE
         IF(NSLOTSINIT.LE.0) THEN
            NSLOTS=0
            WRITE(6,*) "Not using UMATCACHE."
         ELSE
            NPAIRS=NSTATES*(NSTATES+1)/2
            WRITE(6,*) "NPairs: ",NSTATES,NPAIRS
            IF(TSMALL) THEN
               NSLOTS=NSTATES
               tSmallUMat=.TRUE.
               WRITE(6,*) "Using small pre-freezing UMat Cache."
            ELSE
               IF(nMemInit.NE.0) THEN
                  WRITE(6,*) "Allocating ",nMemInit,"Mb for UMatCache+Labels."
                  nSlotsInit=nint((nMemInit*1048576/8)/(nPairs*(nTypes*HElement_t_size+1.0_dp/irat)),sizeof_int)
               ENDIF
               NSLOTS=MIN(NPAIRS, NSLOTSINIT)
               tSmallUMat=.FALSE.
            ENDIF
            UMATCACHEFLAG=0
            WRITE(6,"(A,I3,2I7,I10)") "UMAT NTYPES,NSLOTS,NPAIRS,TOT",NTYPES,NSLOTS,NPAIRS,NSLOTS*NPAIRS*NTYPES
            TUMAT2D=.FALSE.
            ! Each cache element stores <ij|ab> and <ib|aj>.  If real orbitals
            ! then these are identical and we can use this to halve the storage
            ! space (setting nTypes=1).  If not, we must store both explicitly
            ! (nTypes=2).
            Allocate(UMatCacheData(0:nTypes-1,nSlots,nPairs), STAT=ierr)
            call LogMemAlloc('UMatCache',nTypes*nSlots*nPairs,8*HElement_t_size,thisroutine,tagUMatCacheData)
            Allocate(UMatLabels(nSlots,nPairs), STAT=ierr)
            CALL LogMemAlloc('UMATLABELS',nSlots*nPairs,4,thisroutine,tagUMatLabels)
            Memory=(REAL(nTypes*nSlots,dp)*nPairs*8.0_dp*HElement_t_size+nSlots*nPairs*4.0_dp)*9.536743316e-7_dp
            WRITE(6,"(A,G20.10,A)") "Total memory allocated for storage of integrals in cache is: ",Memory,"Mb/Processor"

            UMatCacheData=(0.0_dp)
            UMATLABELS(1:nSlots,1:nPairs)=0
!If tSmallUMat is set here, and have set tCacheFCIDUMPInts, then we need to read in
!the <ik|u|jk> integrals from the FCIDUMP file, then disperse them using the
!FillUMatCache routine. Otherwise, we need to read in all the integrals.
            if (.not.tSmallUMat.and.tReadInCache) then
                write (6,*) 'reading in cache'
                call ReadInUMatCache
            end if
         ENDIF
      END SUBROUTINE SetupUMatCache



      SUBROUTINE SETUPUMAT2D(G1,HarInt)
         ! Set up UMat2D for storing the <ij|u|ij> and <ij|u|ji> integrals,
         ! and pre-calculate the common integrals (<ij|u|ij>, <ij|u|ji>,
         ! <i|v_har|j>) for CPMD calculations.
         ! In:
         !    G1: symmetry and momentum information on the basis functions.
         ! Out:
         !    HarInt(i,j)=<i|v_har|j>, where v_har is the Hartree potential.
         IMPLICIT NONE
         TYPE(BasisFN) G1(*)
         INTEGER ierr
         complex(dp) HarInt(nStates,nStates)
         character(len=*),parameter :: thisroutine='SETUPUMAT2D'
         IF((NSLOTSINIT.LT.0)) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            Allocate(UMat2D(nStates,nStates),STAT=ierr)
            call LogMemAlloc('UMat2D',nStates**2,8*HElement_t_size,thisroutine,tagUMat2D,ierr)
            CALL CPMDANTISYMINTEL(G1,UMAT2D,HarInt,NSTATES)
         ENDIF
      END SUBROUTINE SetupUMat2D



      SUBROUTINE SETUPUMAT2D_DF()
         ! Set up UMat2D for storing the <ij|u|ij> and <ij|u|ji> integrals for
         ! density fitting calculations.
         IMPLICIT NONE
         INTEGER ierr
         character(len=*),parameter :: thisroutine='SETUPUMAT2D_DF'
         IF(NSLOTSINIT.LT.0) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            Allocate(UMat2D(nStates,nStates),STAT=ierr)
!            WRITE(6,*) "nStates for UMat2D: ",nStates
            call LogMemAlloc('UMat2D',nStates**2,8*HElement_t_size,thisroutine,tagUMat2D,ierr)
            IF(tRIIntegrals.or.tCacheFCIDUMPInts) THEN
        !        CALL ReadRI2EIntegrals(nStates,UMat2D,tUMat2D)
         !  This happens later
            ELSE
                CALL ReadDalton2EIntegrals(nStates,UMat2D,tUMat2D)
            ENDIF
         ENDIF
      END SUBROUTINE SetupUMat2D_DF



      SUBROUTINE SETUMATTRANS(TRANS)
         ! In:
         !    Trans: Translation list of orbitrals from one ordering to a new one.
         ! Currently only called in cpmdinit to re-order states by the
         ! one-particle energies (option is rarely used).
         ! Copy to UMatCache's translation table.
         IMPLICIT NONE
         INTEGER TRANS(NSTATES),ierr
         character(*), parameter :: thisroutine='SetupUMatTrans'
         Allocate(TransTable(nStates),STAT=ierr)
         call LogMemAlloc('TransTable',nStates,4,thisroutine,tagTransTable,ierr)
         CALL NECI_ICOPY(NSTATES,TRANS,1,TransTable,1)
         TTRANSGTID=.TRUE.
      END SUBROUTINE SetUMatTrans



      SUBROUTINE SetupUMatTransTable(OldNew,nOld,nNew)
         ! Set up translational table for freezing.
         ! In:
         !    nOld: # of old states.
         !    nNew: # of new states.
         !    OldNew: convert index in the old (pre-freezing) indexing scheme to
         !            the new (post-freezing) indexing scheme.
         IMPLICIT NONE
         INTEGER nNew,nOld,I
         INTEGER OldNew(*),ierr
         LOGICAL tDiff
         character(*), parameter :: thisroutine='SetupUMatTransTable'
         Allocate(TransTable(nNew/2),STAT=ierr)
         call LogMemAlloc('TransTable',nNew/2,4,thisroutine,tagTransTable,ierr)
         Allocate(InvTransTable(nOld/2),STAT=ierr)
         call LogMemAlloc('InvTransTable',nOld/2,4,thisroutine,tagInvTransTable,ierr)
         InvTransTable(1: nOld/2)=0
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
      END SUBROUTINE SetupUMatTransTable



      SUBROUTINE DESTROYUMATCACHE
         IMPLICIT NONE
         character(len=*), parameter :: thisroutine='DESTROYUMATCACHE'
         CALL WriteUMatCacheStats()
         IF(ASSOCIated(UMatCacheData)) THEN
            WRITE(6,*) "Destroying UMatCache"
            CALL LogMemDealloc(thisroutine,tagUMatCacheData)
            Deallocate(UMatCacheData)
            CALL LogMemDealloc(thisroutine,tagUMATLABELS)
            Deallocate(UMatLabels)
         end if
         IF(ASSOCIated(UMat2D)) THEN
            CALL LogMemDealloc(thisroutine,tagUMat2D)
            Deallocate(UMat2D)
         ENDIF
         IF(ASSOCIated(TransTable)) THEN
            CALL LogMemDealloc(thisroutine,tagTransTable)
            Deallocate(TransTable)
         ENDIF
         IF(ASSOCIated(InvTRANSTABLE)) THEN
            CALL LogMemDealloc(thisroutine,tagInvTransTable)
            Deallocate(InvTRANSTABLE)
         ENDIF
      END SUBROUTINE DESTROYUMATCACHE



      SUBROUTINE WriteUMatCacheStats
         IMPLICIT NONE
         IF(ASSOCIated(UMatCacheData)) THEN
            WRITE(6,*) "UMAT Cache Statistics"
            WRITE(6,*) NHITS, " hits"
            WRITE(6,*) NMISSES, " misses"
            WRITE(6,*) iCacheOvCount, " overwrites"
            if(NHITS+NMISSES.gt.0) then
                WRITE(6,"(F6.2,A)") (NHITS/(NHITS+NMISSES+0.0_dp))*100,"% success"
            endif
         ENDIF
      END SUBROUTINE WriteUMatCacheStats



      SUBROUTINE SETUMATCACHEFLAG(NEWFLAG)
         ! Change caching mode of UMatCache,
         ! In:
         !    NewFlag [0,1]: new value for UMatCacheFlag.
         !  flag=1: Storing just the <ik|u|jk> integrals in order they arrive
         !  in (and only have room to do so).
         !  flag=0: Distribute integrals throughout the cache in the scheme
         !  described at the top.
         IMPLICIT NONE
         INTEGER NEWFLAG
         SELECT CASE(UMATCACHEFLAG)
         CASE(1)
!  We were in direct cache mode where values were distributed correctly throughout the cache.
            IF(NEWFLAG.EQ.0.AND..NOT.tSmallUMat) THEN
!  We need to fill the cache properly with values from the small cache.
               CALL FILLUPCACHE()
            ENDIF
         ENDSELECT
         UMATCACHEFLAG=NEWFLAG
         SELECT CASE(NEWFLAG)
         CASE(1)
            IF(NSLOTS.EQ.NPAIRS) THEN ! we're storing every element, so we don't need to deal with different cacheing
               UMATCACHEFLAG=0
            ELSE
               UMATLABELS(1:NSLOTS,1:NPAIRS)=0
!Turn on the direct caching, and clear the cache.
            ENDIF
         ENDSELECT
         RETURN
      END SUBROUTINE SETUMATCACHEFLAG



      SUBROUTINE FillUpCache()
         ! Disperse the (pre-calculated) <ik|u|jk> integrals throughout the cache.
         ! The cache consists of an unordered set (in the standard UMatCache
         ! sense) of labels and elements.
         ! We must order this, and then distribute the elements throughout each set of SLOTS.
         IMPLICIT NONE
         INTEGER I,J,K,N,nK
         DO I=1,nPairs
! Find the last value in the cache
            CALL BinarySearch(nPairs+1,UMatLabels(1:nSlots,I),1,nSlots,N,J,K)
            N=J-1
! N is now the last element and thus number of elements.
! Sort according to label
            call sort (UMatLabels(1:N, i), UMatCacheData(:, 1:N, i))
            K=nSlots
! Now copy element among the whole array for this, from the end.
! (Multiple copies of the same integral in an unfilled cache make adding
! elements ! substantially faster.)
            DO J=N,1,-1
               nK=(nSlots*(J-1))/N+1
               UMatLabels(nK:K,I)=UMatLabels(J,I)
               DO K=K,nK,-1
                  UMatCacheData(:,K,I)=UMatCacheData(:,J,I)
               ENDDO
               K=nK
            ENDDO
         ENDDO
      END SUBROUTINE FillUpCache



      SUBROUTINE BINARYSEARCH(VAL,TAB,A,B,LOC,LOC1,LOC2)
!   A binary search to find VAL in TAB.  TAB is sorted, but can have
!   multiple entries being the same.  If the search terminated unsuccessfully,
!   the entry indicated is one after half-way through the set of entries which
!   would be immediately prior to it.  From here until the label changes
!   should be filled with VAL if it is to be entered into the table.
!   A and B are the limits of the table.
!   If the search is successful, the location of VAL in TAB is returned in LOC
!   (and LOC1,LOC2).
!   If the search fails, then VAL should fit between LOC1 and LOC2 in TAB.
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
!Finally, check the last element of the array, as it can still be there.
         IF(TAB(B).eq.VAL) THEN
             LOC=B
             LOC1=B
             LOC2=B
             RETURN
         ENDIF
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
      END SUBROUTINE BinarySearch



!   Get a unique index corresponding to pair (I,J), and return in RET.
!  New Scheme 1/2/07
!   e.g. 11 12 13 14 15 corresponds to    1  2  4  7 11
!           22 23 24 25                      3  5  8 12
!              33 34 35                         6  9 13
!                 44 45                           10 14
!                    55                              15
      SUBROUTINE GETCACHEINDEX(I,J,RET)
         ! In:
         !    I,J (I<=J): state indices
         ! Out:
         !    Cache indexing scheme.
         IMPLICIT NONE
         INTEGER I,J,RET
         RET=J*(J-1)/2+I
      END SUBROUTINE GetCacheIndex



! Example: calculating int(sqrt(2*ind)) from 2*ind.
!    2ind                    sqrt(2ind)
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
! But indexing scheme:
!   2ind                      J
!     2  4  8 14 22... 92  -> 1  2  3  4  5 ...10
!        6 10 16 24    94        2  3  4  5    10
!          12 18 26    96           3  4  5    10
!             20 28    98              4  5    10
!                30   100                 5    10
! Hence need to +1 to J in some cases.
      SUBROUTINE GETCACHEINDEXSTATES(IND,I,J)
         ! In:
         !   Ind: Cache index.
         ! Out:
         !   I,J (I<=J): states corresponding to cache index.
         ! Reverse of GetCacheIndex.
         IMPLICIT NONE
         INTEGER I,J,IND
         J=int(SQRT(2.0d0*IND))
         IF(J*(J+1)/2.LT.IND) J=J+1
         I=IND-J*(J-1)/2
      END SUBROUTINE GetCacheIndexStates

      SUBROUTINE FreezeUMatCache(OrbTrans,nOld,nNew)
         ! We're in the middle of freezing some orbitals.
         ! OrbTrans(i) will give us the new position of the old orbital i.
         IMPLICIT NONE
         INTEGER nOld,nNew,OrbTrans(nOld)
         INTEGER onSlots,onPairs
         if(nNew/2.NE.nStates.OR.tSmallUMat) THEN
            WRITE(6,*) "Reordering UMatCache for freezing"
            onSlots=nSlots
            onPairs=nPairs
            CALL FreezeUMatCacheInt(OrbTrans,nOld,nNew,onSlots,onPairs)
         else
            WRITE(6,*) "UMatCache size not changing.  Not reordering."
         endif
      END SUBROUTINE FreezeUMatCache



      SUBROUTINE FreezeUMAT2D(OldBasis,NewBasis,OrbTrans,iSS)
         IMPLICIT NONE
         INTEGER NewBasis,OldBasis,iSS,ierr,OrbTrans(OldBasis),i,j
         HElement_t(dp),POINTER :: NUMat2D(:,:)
         integer(TagIntType) :: tagNUMat2D=0
         character(len=*),parameter :: thisroutine='FreezeUMat2D'

         Allocate(NUMat2D(NewBasis/iSS,NewBasis/iSS),STAT=ierr)
         call LogMemAlloc('UMat2D',(NewBasis/iSS)**2,8*HElement_t_size,thisroutine,tagNUMat2D,ierr)
         NUMat2D(:,:)=(0.0_dp)
         DO i=1,OldBasis/2
            IF(OrbTrans(i*2).NE.0) THEN
                DO j=1,OldBasis/2
                    IF(OrbTrans(j*2).NE.0) THEN
                        NUMat2D(OrbTrans(i*2)/2,OrbTrans(j*2)/2)=UMat2D(i,j)
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        call LogMemDealloc(thisroutine,tagUMat2D)
        Deallocate(UMat2D)
        UMat2D=>NUMat2D
        NULLIFY(NUMat2D)
        tagUMat2D=tagNUMat2D
        RETURN
      END SUBROUTINE FreezeUMAT2D



      SUBROUTINE FreezeUMatCacheInt(OrbTrans,nOld,nNew,onSlots,onPairs)
         IMPLICIT NONE
         INTEGER nOld,nNew,OrbTrans(nOld)
         HElement_t(dp),Pointer :: NUMat2D(:,:) !(nNew/2,nNew/2)
         integer(TagIntType) :: tagNUMat2D=0
         HElement_t(dp) El(0:nTypes-1)
         INTEGER i,j,k,l,m,n
         INTEGER ni,nj,nk,nl,nm,nn,A,B,iType
         HElement_t(dp),Pointer :: OUMatCacheData(:,:,:) !(0:nTypes-1,onSlots,onPairs)
         INTEGER,Pointer :: OUMatLabels(:,:) !(onSlots,onPairs)

         INTEGER onSlots,onPairs,ierr
         LOGICAL toSmallUMat,tlog,toUMat2D,tmpl
         character(len=*),parameter :: thisroutine='FreezeUMatCacheInt'

         toUMat2D=tUMat2D
         IF(tUMat2D) then
            Allocate(NUMat2D(nNew/2,nNew/2),STAT=ierr)
            call LogMemAlloc('UMat2D',(nNew/2)**2,8*HElement_t_size,thisroutine,tagNUMat2D,ierr)
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
            CALL LogMemDealloc(thisroutine,tagUMat2D)
            Deallocate(UMat2D)
            UMat2D=>NUMat2D
            nullify(NUMat2D)
            tagUMat2D=tagNUMat2D
         endif
! Now go through the other cache.
! First save the memory used for it.
!         onSlots=nSlots
!         onPairs=nPairs
         OUMatCacheData=>UMatCacheData
         tagOUMatCacheData=tagUMatCacheData
         OUMatLabels=>UMatLabels
         tagOUMatLabels=tagUMatLabels
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
             CALL GetCacheIndex(i,k,m)
             DO n=1,onSlots


              tmpl = .true.
              if (n /= 1) then
                  if (OUMatLabels(n,m) /= OUMatLabels(n-1,m)) tmpl = .false.
              endif
              if ( (onSlots == onPairs .or. toSmallUMat) .or. &
                   (onSlots /= onPairs .and. tmpl) ) then
               IF(OUMatLabels(n,m).NE.0) THEN
                ni=OrbTrans(i*2)/2
                nk=OrbTrans(k*2)/2
!Now get the label of the slot and convert to orbitals
                IF(onSlots.EQ.onPairs) THEN
                 CALL GetCacheIndexStates(n,j,l)
                ELSEIF(toSmallUMat) THEN
                  j=n
                  l=n
                ELSE
                 CALL GetCacheIndexStates(OUMatLabels(n,m),j,l)
                ENDIF
                nj=OrbTrans(j*2)/2
                nl=OrbTrans(l*2)/2
                IF(nj.NE.0.AND.nl.NE.0) THEN
                   ! JSS: Alex, please check!
                   ! GetCachedUMatEl called for cache indices: El(O) is
                   ! a dummy argument.
                   Tlog=GetCachedUmatEl(ni,nj,nk,nl,El(0),nm,nn,A,B,iType)
                   CALL CacheUMatEl(B,OUMatCacheData(0,n:n+nTypes,m:m+nTypes),nm,nn,iType)
                ENDIF
               ENDIF
              ENDIF
             ENDDO
            ENDIF
           ENDDO
          ENDIF
         ENDDO
         CALL LogMemDealloc(thisroutine,tagOUMatLabels)
         Deallocate(OUMatLabels)
         CALL LogMemDealloc(thisroutine,tagOUMatCacheData)
         Deallocate(OUMatCacheData)
         CALL SetUMatCacheFlag(0)
      END SUBROUTINE FreezeUMatCacheInt


      SUBROUTINE CacheFCIDUMP(I,J,K,L,Z,CacheInd,ZeroedInt,NonZeroInt)
          IMPLICIT NONE
          INTEGER :: I,J,K,L,CacheInd(nPairs)
          INTEGER(int64) :: ZeroedInt,NonZeroInt
          INTEGER :: A,B
          HElement_t(dp) :: Z

          IF(abs(Z).lt.UMatEps) THEN
!We have an epsilon cutoff for the size of the two-electron integrals - UMatEps
              ZeroedInt=ZeroedInt+1
              RETURN
          ELSE
              NonZeroInt=NonZeroInt+1
          ENDIF

!Find unique indices within permutational symmetry.
          IF(K.lt.I) THEN
              CALL SWAP(I,K)
          ENDIF
          IF(L.lt.J) THEN
              CALL SWAP(J,L)
          ENDIF
          CALL GETCACHEINDEX(I,K,A)
          CALL GETCACHEINDEX(J,L,B)
          IF(A.gt.B) THEN
              CALL SWAP(A,B)
              CALL SWAP(I,J)
              CALL SWAP(K,L)
          ENDIF

!          WRITE(6,*) "Final Phys ordering: ",I,J,K,L
!          WRITE(6,*) "Pair indices: ",A,B

          IF(A.gt.nPairs) THEN
              WRITE(6,*) "Final Phys ordering: ",I,J,K,L
              WRITE(6,*) "Pair indices: ",A,B
              WRITE(6,*) "nPairs,nSlots: ", nPairs,nSlots
              CALL Stop_All("CacheFCIDUMP","Error in caching")
          ENDIF

!Store the integral in a contiguous fashion. A is the index for the i,k pair
          IF(UMATLABELS(CacheInd(A),A).ne.0) THEN
              IF((abs(REAL(UMatCacheData(nTypes-1,CacheInd(A),A),dp)-Z)).gt.1.0e-7_dp) THEN
                  WRITE(6,*) i,j,k,l,z,UMatCacheData(nTypes-1,CacheInd(A),A)
                  CALL Stop_All("CacheFCIDUMP","Same integral cached in same place with different value")
              ENDIF

              CALL Stop_All("CacheFCIDUMP","Overwriting UMATLABELS")
          ENDIF
          UMATLABELS(CacheInd(A),A)=B
          IF(.not. near_zero(REAL(UMatCacheData(nTypes-1,CacheInd(A),A),dp))) THEN
              CALL Stop_All("CacheFCIDUMP","Overwriting when trying to fill cache.")
          ENDIF
          UMatCacheData(nTypes-1,CacheInd(A),A)=Z

          CacheInd(A)=CacheInd(A)+1
          IF(CacheInd(A).gt.nSlots) THEN
              CALL Stop_All("CacheFCIDUMP","Error in filling cache")
          ENDIF

!After we have filled all of the cache, this will want to be sorted.
      END SUBROUTINE CacheFCIDUMP


!This is a routine to calculate the maximum size needed for nSlotsInit to hold all the needed two-electron
!integrals in the cache.
!We are assuming that there is no more than one integral per i,j pair and so permutational
!symmetry is taken into account when determining islotsmax.
      SUBROUTINE CalcNSlotsInit(I,J,K,L,Z,nPairs2,MaxSlots)
          IMPLICIT NONE
          INTEGER :: I,J,K,L,nPairs2,MaxSlots(1:nPairs2),A,B,C,D,X,Y
          HElement_t(dp) :: Z

!The (ii|jj) and (ij|ij) integrals are not stored in the cache (they are stored in UMAT2D, so
!we do not want to include them in the consideration of the size of the cache.
          IF((I.eq.J).and.(K.eq.L)) THEN
              RETURN
          ELSEIF((I.eq.K).and.(J.eq.L)) THEN
              RETURN
          ELSEIF(min(I,J,K,L).eq.0) THEN
              RETURN
          ENDIF
          A=I
          B=J
          C=K
          D=L
          IF(B.lt.A) THEN
              CALL SWAP(B,A)
          ENDIF
          IF(D.lt.C) THEN
              CALL SWAP(C,D)
          ENDIF
          CALL GETCACHEINDEX(A,B,X)
          CALL GETCACHEINDEX(C,D,Y)
          IF(X.gt.Y) THEN
              CALL SWAP(X,Y)
              CALL SWAP(A,B)
              CALL SWAP(C,D)
          ENDIF

          IF(abs(Z).gt.UMatEps) THEN
              MaxSlots(X)=MaxSlots(X)+1
              IF(X.gt.nPairs2) THEN
                  CALL Stop_All("CalcNSlotsInit","Problem since X > nPairs")
              ENDIF
              IF((MaxSlots(X).gt.nPairs2).and.(.not.tROHF)) THEN
                  WRITE(6,*) "Final Phys ordering: ",I,J,K,L
                  WRITE(6,*) "Pair indices: ",X,Y
                  WRITE(6,*) "nPairs,nSlots: ", nPairs2,MaxSlots(X)
                  CALL Stop_All("CalcNSlotsInit","Problem since more integrals for a given ik pair found than possible.")
              ENDIF
          ENDIF

      END SUBROUTINE CalcNSlotsInit


      subroutine ReadInUMatCache()
      ! Read in cache file from CacheDump.
      implicit none
      integer  i,j,k,l,iCache1,iCache2,A,B,readerr,iType, iunit
      HElement_t(dp) UMatEl(0:nTypes-1),DummyUMatEl(0:nTypes-1)
      logical  tDummy,testfile
      inquire(file="CacheDump",exist=testfile)
      if (.not.testfile) then
          write (6,*) 'CacheDump does not exist.'
          return
      end if
      iunit = get_free_unit()
      open (iunit,file="CacheDump",status="old",iostat=readerr)
      if (readerr.ne.0) then
          write (6,*) 'Error reading CacheDump.'
          return
      end if
      read (iunit,*) nStatesDump
      readerr=0
      do while (readerr.eq.0)
        read (iunit,*,iostat=readerr) i,j,k,l,UMatEl
        DummyUMatEl=UMatEl
        if (TTRANSFINDX) then
            i=TransTable(i)
            j=TransTable(j)
            k=TransTable(k)
            l=TransTable(l)
        end if
        if (min(i,j,k,l).gt.0.and.max(i,j,k,l).le.nStates) then
            ! Need to get cache indices before we cache the integral:
            ! a dummy call to GetCachedUMatEl returns the needed indices and
            ! integral type information.
            tDummy=GetCachedUMatEl(i,j,k,l,DummyUmatEl(0),iCache1,iCache2,A,B,iType)
            call CacheUMatEl(B,UMatEl,iCache1,iCache2,iType)
        end if
      end do
      close(iunit)
      return
      end subroutine ReadInUMatCache



      subroutine DumpUMatCache()
      ! Print out the cache contents so they can be read back in for a future
      ! calculation.  Need to print out the full set of indices, as the number of
      ! states may change with the next calculation.
      implicit none
      ! Variables
      integer iPair,iSlot,i,j,k,l,iCache1,iCache2,A,B,iType
      HElement_t(dp) UMatEl
      type(Symmetry) Sym
      integer iunit
      iunit = get_free_unit()
      open (iunit,file="CacheDump",status="unknown")
!      do i=1,nPairs !Run through ik pairs
!          do j=1,nSlots !Run through all pairs (unordered in the list)
!              WRITE(iunit,*) i,j,UMatLabels(j,i),UMatCacheData(:,j,i)  !ik label, slot value, jl label, integral
!          enddo
!      enddo
!      WRITE(iunit,*) "*****"
      write (iunit,*) nStates
      do iPair=1,nPairs
        do iSlot=iPair,nSlots
          call GetCacheIndexStates(iPair,i,k)
          call GetCacheIndexStates(iSlot,j,l)
          Sym=TotSymRep()
          ! All integrals stored in the cache are non-zero by symmetry.
          if (LSymSym(Sym)) then
              if (.not.GetCachedUMatEl(i,j,k,l,UMatEl,iCache1,iCache2,A,B,iType)) then
                  if (TTRANSFINDX) then
                      i=InvTransTable(i)
                      j=InvTransTable(j)
                      k=InvTransTable(k)
                      l=InvTransTable(l)
                  end if
                  ! Print out UmatCacheData as UMatEl holds just a single
                  ! integral.
                  write (iunit,*) i,j,k,l,UMatCacheData(:,ICACHE2,ICACHE1)!,A,B
              end if
          end if
        end do
      end do
      close(iunit,status="keep")
      return
      end subroutine DumpUMatCache



      logical function HasKPoints()
         IMPLICIT NONE
         IF(NKPS.GT.1) THEN
            HasKPoints=.TRUE.
         ELSE
            HasKPoints=.FALSE.
         ENDIF
      end function HasKPoints

    elemental function GTID (gInd) result(id)

        ! Convert spin orbital index to spacial orbital index if required
        ! (e.g. for restricted calculation)
        !
        ! In:  gInd - Spin orbital index
        ! Ret: id   - Overall index

        integer, intent(in) :: gInd
        integer :: id

        if (tStoreSpinOrbs) then
            ! Storing as spin-orbitals (UHF/default ROHF)
            id = gInd
        else
            ! Storing as spatial orbitals (RHF or explicit input option ROHF)
            id = (gInd-1)/2 + 1
        endif
        if (tTransGTID) id = TransTable(id)
    end function

    elemental function spatial(spin_orb) result(spat_orb)

        ! Convert a spin orbital label into a spatial one. This is more
        ! appropriate than the above function, gtid, when one always wants the
        ! spatial label, regardless of the way integrals are being stored.

        integer, intent(in) :: spin_orb
        integer :: spat_orb

        spat_orb = (spin_orb-1)/2 + 1

    end function spatial

      LOGICAL FUNCTION GETCACHEDUMATEL(IDI,IDJ,IDK,IDL,UMATEL,ICACHE,ICACHEI,A,B,ITYPE)
         ! In:
         !    IDI,IDJ,IDK,IDL: orbitals indices of the integral(s). These are
         !    rearranged, as disucssed below.
         ! Out:
         !    UMatEl: <ij|u|kl> integral (if found in the cache).
         !    A,B: Cache indices of (I,K) and (J,L) pairs (see GetCacheIndex).
         !    ICache,ICacheI: location of integral within the cache (or where it
         !                    should be stored if not found).
         !    iType: gives information (bit-wise) on where the requested
         !           integral lies within the cache slot.  Useful if the
         !           integral has to be computed (using the reordered indices).

! Lookup in the cache to see if there's a stored element.  If not, return TRUE,
! together with the information on how to compute the integral (re-ordered
! indices, cache indices and a "type" index).

! If the integral exists in the cache, return false and return the stored element in UMatEl.

! NOTE: This will rearrange IDI,IDJ,IDK,IDL into the correct order
! (i,k)<=(j,l) and i<=k, j<=l.  ICACHE corresponds to the pair (i,j), and
! ICACHEI is the index in that cache where the cache should be located.
! Note: (i,k)>(j,l) := (k>l) || ((k==l)&&(i>j))

!         use SystemData, only : nBasis,G1
         IMPLICIT NONE
         INTEGER IDI,IDJ,IDK,IDL,ICACHE,ICACHEI
         INTEGER ICACHEI1,ICACHEI2
         HElement_t(dp) UMATEL
         INTEGER A,B,ITYPE,ISTAR,ISWAP
!         LOGICAL tDebug
!         IF(IDI.eq.14.and.IDJ.eq.17.and.IDK.eq.23.and.IDL.eq.6) THEN
!             WRITE(6,*) "Setting tDebug!"
!             tDebug=.true.
!         ELSE
!             tDebug=.false.
!         ENDIF
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
               CALL GETCACHEINDEX(IDJ,IDL,A)
            ELSEIF(IDJ.EQ.IDL) THEN
               B=IDJ
               IF(IDK.LT.IDI) THEN
                  CALL SWAP(IDI,IDK)
                  ITYPE=2
               ENDIF
               CALL GETCACHEINDEX(IDI,IDK,A)
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
            CALL GETCACHEINDEX(IDI,IDK,A)
            CALL GETCACHEINDEX(IDJ,IDL,B)
            IF(A.GT.B) THEN
               CALL SWAP(A,B)
               CALL SWAP(IDI,IDJ)
               CALL SWAP(IDK,IDL)
               ISWAP=1
            ENDIF
            IF(HElement_t_size.EQ.1) THEN
!  Eight integrals from ijkl are the same.
               ITYPE=0
            ELSE
!  Complex orbitals and integrals, so we need to consider different types
!  Using notation abcd rather than ijkl.  6/2/07 and 19/2/06
!
!  <> mean swap pairs <ab|cd> -> <ba|dc>
!  *.  means complex conjugate of first codensity i.e. <ab|cd> -> <cb|ad>
!  .* for second and ** for both.
!
!  abcd   -> badc <>
!  |  |
!  | \|/       |-> cdab ** -> dcba **<>
!  | cbad *.  -|
!  |           |-> bcda *.<>
! \|/
!  adcb .* -> dabc .*<>

!Now consider what must occur to the other integral in the slot to recover pair
!(abcd,cbad).  0 indicates 1st in slot, 1 means 2nd.  * indicated conjg.
!
!  ..    abcd  cbad  0  1
!  *.    cbad  abcd  1  0
!  .*    adcb  cbad  1* 0*
!  **    cdab  adcb  0* 1*
!    <>  badc  dabc  0  1*
!  *.<>  bcda  dcba  1* 0
!  .*<>  dabc  badc  1  0*
!  **<>  dcba  bcda  0* 1


! Of the type, bit zero indicates which of the two integrals in a slot to use.
!Bit 1 is set if the integral should be complex conjugated.
!   Bit 2 is set if the other integral in the slot should be complex conjugated
!if we are to have the structure (<ij|kl>,<kj|il>) in the slot.
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
!         IF(tDebug) WRITE(6,"(A,8I5)") "GCUE",IDI,IDJ,IDK,IDL,A,B,iType,UMatCacheFlag
         ICACHE=A
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
!                IF(tDebug) THEN
!                     CALL DumpUMatCache(nBasis)
!                     WRITE(8,*) B,NSLOTS
!                     WRITE(8,*) UMATLABELS(1:NSLOTS,A)
!                 ENDIF
                CALL BINARYSEARCH(B,UMATLABELS(1:NSLOTS,A),1,NSLOTS,ICACHEI,ICACHEI1,ICACHEI2)
!                IF(tDebug) WRITE(8,*) "***",UMATLABELS(ICACHEI,A),ICACHEI,ICACHEI1,ICACHEI2
            ENDIF
         ENDIF
         IF(UMATLABELS(ICACHEI,ICACHE).EQ.B) THEN
            !WRITE(6,*) "C",IDI,IDJ,IDK,IDL,ITYPE,UMatCacheData(0:nTypes-1,ICACHEI,ICACHE)
            UMATEL=UMatCacheData(IAND(ITYPE,1),ICACHEI,ICACHE)
#ifdef __CMPLX
            IF(BTEST(ITYPE,1)) UMATEL=CONJG(UMATEL)  ! Bit 1 tells us whether we need to complex conjg the integral
#endif
!   signal success
            GETCACHEDUMATEL=.FALSE.
         ELSE
!   signal failure
            GETCACHEDUMATEL=.TRUE.
!            WRITE(68,*) A,B,ICACHEI1,ICACHEI2,ICACHEI
         ENDIF
         RETURN
      END FUNCTION GETCACHEDUMATEL

      !------------------------------------------------------------------------------------------!

      function numBasisIndices(nBasis) result(nBI)
        implicit none
        integer, intent(in) :: nBasis
        integer :: iSS
        integer :: nBI

        if(tStoreSpinOrbs) then
           iSS = 1
        else
           iSS = 2
        endif

        ! number of distinct indices of the integrals
        nBI = nBasis / iSS
      end function numBasisIndices

      !------------------------------------------------------------------------------------------!

      subroutine SetupUMat2d_dense(nBasis)
        implicit none
        integer, intent(in) :: nBasis
        integer :: nBI, i, j, idX, idN

        nBI = numBasisIndices(nBasis)
        allocate(UMat2D(nBI,nBI))

        do i = 1, nBI
           do j = 1, nBI
              idX = max(i,j)
              idN = min(i,j)
              ! here, we introduce a cheap redundancy in memory to allow
              ! for faster access (no need to get max/min of indices
              ! and have contiguous access)
              ! store the integrals <ij|ij> in UMat2D
              UMat2D(j,i) = get_umat_el(idN,idX,idN,idX)
           end do
        end do
      end subroutine SetupUMat2d_dense

      !------------------------------------------------------------------------------------------!

      subroutine freeUmat2d_dense()
        implicit none

        ! deallocate auxiliary arrays storing the integrals <ij|ij> and <ij|ji>
        if(associated(UMat2d)) deallocate(UMat2d)
      end subroutine freeUmat2d_dense


END MODULE UMatCache
! Still useful to keep CacheUMatEl and GetCachedUMatEl outside of the module for
! CPMD interaction (though this should be fixed: the problem lies with the type
! mismatch (GKElement and ) in the argument list).



! Set an element in the cache.  All the work has been done for us before as the
! element we have to set is in (ICACHEI,ICACHE) iType tells us whether we need
! to swap/conjugate the nTypes integrals within the slot We still need to fill
! out the space before or after  us if we've been put in the middle of a block
! of duplicates.
      SUBROUTINE CACHEUMATEL(B,UMATEL,ICACHE,ICACHEI,iType)
         ! In:
         !    A,B: cache indices of the element.
         !    UMatEl: element being stored.  For calculations involving real
         !            orbitals, this is a array of size 1 containing the
         !            <ij|u|kl> integral (nTypes=1).  For calculations involving
         !            complex orbtials, this is an array of size 2 containing
         !            the <ij|u|kl> and <il|u|jk> integrals (nTypes=2).
         !    ICache: Segment index of the cache for storing integrals involving
         !            index A (often equal to A).
         !    ICacheI: Slot within ICache segment for storing UMatEl involving
         !             B.
         !    iType: See notes below.
         use constants, only: dp
         use UMatCache
         IMPLICIT NONE
         INTEGER B,ICACHE,ICACHEI
         HElement_t(dp) UMATEL(0:NTYPES-1),TMP(0:NTYPES-1)
         INTEGER OLAB,IC1,ITOTAL
         INTEGER iType
         INTEGER iIntPos
         SAVE ITOTAL
         DATA ITOTAL /0/
         if (nSlots.eq.0) return
!         WRITE(6,*) "CU",A,B,UMATEL,iType
!         WRITE(6,*) A,ICache,B,ICacheI
         if(nTypes.gt.1) then
! A number of different cases to deal with depending on the order the integral came in (see GetCachedUMatEl for details)
!  First get which pos in the slot will be the new first pos
            iIntPos=iand(iType,1)
!  If bit 1 is set we must conjg the (to-be-)first integral
#ifdef __CMPLX
            if(btest(iType,1)) then
               Tmp(0)=conjg(UMatEl(iIntPos))
            else
               Tmp(0)=UMatEl(iIntPos)
            endif
#else
            Tmp(0)=UMatEl(iIntPos)
#endif
!  If bit 2 is set we must conjg the (to-be-)second integral
#ifdef __CMPLX
            if(btest(iType,2)) then
               Tmp(1)=conjg(UMatEl(1-iIntPos))
            else
               Tmp(1)=UMatEl(1-iIntPos)
            endif
#else
            Tmp(1)=UMatEl(1-iIntPos)
#endif
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
            if (icachei > 0) OLAB=UMATLABELS(ICACHEI,ICACHE)
         ENDDO
      END SUBROUTINE CacheUMatEl


