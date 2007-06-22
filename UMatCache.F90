!Based on umatcache.F, this modulises the umat cache, and retrieval functions.
MODULE UMatCache
      USE HElement
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


! For the more requently used <ij|u|ij> and <ij|u|ji> integrals, we store them in a separate cache (if TUMat2D is true)
      TYPE(HElement), Pointer :: UMat2D(:,:) !(nStates,nStates)
      LOGICAL tUMat2D
         SAVE tUMat2D ,UMat2D
      

! For the UEG, we damp the exchange interactions.
!    0 means none
!    1 means attenuated (using an erfc)
!    2 means cut-off    (at a distance Rc=ALAT(4))
      INTEGER iPeriodicDampingType

!  Book-keeping information
!  nSlotsInit is the number of slots requested on input.  If the number required is less, then the lower value is allocated
!  If nSlotsInit is set to 0, then general <ij|u|kl> element caching is not performed, but UMat2D <ij|u|ij> and <ij|u|ji> is.  For nSlotsInit=-1 neither is performed.
      INTEGER nSlotsInit
! nHits and nMisses are the number of cache hits and misses.
      INTEGER nHits,nMisses
!.. UMatCacheFlag=0 is normal operation
!.. UMatCacheFlag=1 means cache from bottom up
!..      This is useful when lots of sequential pieces of data are being stored.
!..  When UMatCacheFlag is reset to 0, the data which are present are spread evenly around the slots for a given Pair.
      INTEGER UMatCacheFlag
         SAVE nSlotsInit,nHits,nMisses,UMatCacheFlag

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
         SAVE nAuxBasis,nBasisPairs,tDFInts,DFCoeffs,DFInts

      Contains

! Get the index of physical order UMAT element <IJ|KL>.  Indices are internally reordered such that I>K, J>L,(I,K)>(K,L) 
!NB This is a different order from UMatCache
      INTEGER FUNCTION UMatInd(I,J,K,L)
         IMPLICIT NONE
         INTEGER I,J,K,L,A,B
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
         ELSE
            UMatInd=(B*(B-1))/2+A
         ENDIF
      END

!Get the prospective size of a UMat (not a UMatCache) for completely storing FCIDUMP 2-e integrals
!  The UMat is currently passed as a parameter, but in future be absorbed into UMatCache.
      SUBROUTINE GetUMatSize(nBasis,iSS,iSize)
         IMPLICIT NONE
         INTEGER nBasis,iSS
         INTEGER iPairs,nBi,iSize
         nBi=nBasis/iSS
         iPairs=(nBi*(nBi+1))/2
         iSize=(iPairs*(iPairs+1))/2
      END

!Get a U matrix element <ij|u|kl> in multifarious ways.  Either from a passed-in UMAT, or ALAT parameters, 
! or from UMatcache.
      FUNCTION GetUMatEl(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,IDI,IDJ,IDK,IDL)
!         USE HElement
         IMPLICIT NONE
         TYPE(HElement) GetUMatEl
         INTEGER NBASISMAX(5,3),I,J,K,L,NHG,ISS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(NHG)
         REAL*8 ALAT(3)
         TYPE(HElement) UMAT(*)
         TYPE(HElement) UElems(0:nTypes-1)
         INTEGER A,B,C
         INTEGER IDI,IDJ,IDK,IDL
         REAL*8 SUM
         PARAMETER PI=3.14159265358979323846264338327950288419716939937510D0
         INTEGER ICACHE,ICACHEI,ITYPE
         LOGICAL LSYMSYM
         TYPE(Symmetry) SYM,SYMPROD,SYMCONJ
         INTEGER ISUB,ISUB2
         LOGICAL GetCachedUMatEl
!         CALL TISET(' GETUMATEL',ISUB)
!   IF NBASISMAX(1,3) is less than zero, we directly give the integral.
!   Otherwise we just look it up in umat
         IF(NBASISMAX(1,3).GE.0) THEN
!   See if we need to calculate on the fly
          IF(ISS.EQ.0) THEN

!  JSS - store <ij|ij> and <ij|ji> in UMAT2D.
!   Remember permutations.  Complex case is the same as real.
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
           ELSE IF (IDI.eq.IDJ.and.IDK.eq.IDL.AND.TUMAT2D.AND.HElementSize.EQ.1) THEN
!   <ii|jj> = <ij|ji> Only for real systems
             I=MAX(IDI,IDK)
             J=MIN(IDI,IDK)
             GETUMATEL=UMAT2D(I,J)
           ELSE
!   Check to see if the umat element is in the cache
               I=IDI
               J=IDJ
               K=IDK
               L=IDL
               SYM%s=1
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
!                CALL GetKPInd(I,J,K,L)
                IF(TTRANSFINDX) THEN
                 I=InvTransTable(I)
                 J=InvTransTable(J)
                 K=InvTransTable(K)
                 L=InvTransTable(L)
                ENDIF
               ENDIF
!    This will rearrange I,J,K,L into the correct order
!   (i,k)<=(j,l) and i<=k, j<=l.
               IF(GETCACHEDUMATEL(I,J,K,L,GETUMATEL,ICACHE,ICACHEI,A,B,ITYPE)) THEN
!   We don't have a stored UMAT - we call to generate it.
                  IF(tDFInts) THEN
!   We're using density fitting
                     Call GetDF2EInt(I,J,K,L,UElems)
                     GetUMatEl=UElems(0)
                  ELSE
!   Otherwise we call CPMD
                     IF(TTRANSFINDX) THEN
!         WRITE(6,"(A,4I5,$)") "MM",TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L)
                        CALL INITFINDXI(TRANSTABLE(I),TRANSTABLE(J),TRANSTABLE(K),TRANSTABLE(L),UElems)
                     ELSE
                        CALL INITFINDXI(I,J,K,L,UElems)
                     ENDIF
                     GETUMATEL=UElems(IAND(ITYPE,1))
                     IF(ITYPE.GT.1) GETUMATEL=DCONJG(GETUMATEL)
!                     WRITE(6,*) "I",I,J,K,L,UElems
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
                  IF(ICACHE.NE.0) CALL CACHEUMATEL(A,B,UElems,ICACHE,ICACHEI)
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
             GETUMATEL=UMAT(UMatInd(IDI,IDJ,IDK,IDL))
          ENDIF
         ELSEIF(NBASISMAX(1,3).EQ.-1) THEN
            CALL GetUEGUmatEl(IDI,IDJ,IDK,IDL,ISS,G1,ALAT,iPeriodicDampingType,GetUMatEl)
         ENDIF
!         WRITE(6,"(4I5,$)") IDI,IDJ,IDK,IDL
!         WRITE(6,*) GETUMATEL,ABS(GETUMATEL)
         RETURN
!         CALL TIHALT(' GETUMATEL',ISUB)
      END

      logical function HasKPoints()
         IMPLICIT NONE
         include 'cpmddata.inc'
         INTEGER I,J,K,L
         IF(NKPS.GT.1) THEN
            HasKPoints=.TRUE.
         ELSE
            HasKPoints=.FALSE.
         ENDIF
      end


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
         IF(NSLOTSINIT.LE.0) THEN
            NSLOTS=0
            WRITE(6,*) "Not using UMATCACHE."
         ELSE
            NSTATES=NSTATE
            NPAIRS=NSTATES*(NSTATES+1)/2
            IF(TSMALL) THEN
               NSLOTS=NSTATES
               tSmallUMat=.TRUE.
               WRITE(6,*) "Using small pre-freezing UMat Cache."
            ELSE
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
         ENDIF
      END

      SUBROUTINE SETUPUMAT2D(G1)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER ierr
         IF(NSLOTSINIT.LT.0) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            Allocate(UMat2D(nStates,nStates),STAT=ierr)
            Call MemAlloc(ierr,UMat2D,HElementSize*NSTATES*NSTATES,'UMAT2D')
            CALL CPMDANTISYMINTEL(G1,UMAT2D,NSTATES)
         ENDIF
      END


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
            CALL ReadDalton2EIntegrals(nStates,UMat2D)
         ENDIF
      END

     
      SUBROUTINE SETUMATTRANS(TRANS)
         IMPLICIT NONE
         INTEGER TRANS(NSTATES),ierr
         Allocate(TransTable(nStates),STAT=ierr)
         CALL MemAlloc(ierr,TransTable,NSTATES,'TransTable')
         CALL ICOPY(NSTATES,TRANS,1,TransTable,1)
         TTRANSGTID=.TRUE.
      END
 
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
      END
      SUBROUTINE DESTROYUMATCACHE
         IMPLICIT NONE
         IF(Allocated(UMatCacheData)) THEN
            CALL MemDealloc(UMatCacheData)
            Deallocate(UMatCacheData)
            CALL MemDealloc(UMATLABELS)
            Deallocate(UMatLabels)
            IF(Allocated(UMat2D)) THEN
               CALL MemDealloc(UMat2D)
               Deallocate(UMat2D) 
            ENDIF
            IF(Allocated(TransTable)) THEN
               CALL MemDealloc(TransTable)
               Deallocate(TransTable)
            ENDIF
            IF(Allocated(InvTRANSTABLE)) THEN
               CALL MemDealloc(InvTransTable)
               Deallocate(InvTRANSTABLE)
            ENDIF
            WRITE(6,*) "UMAT Cache Statistics"
            WRITE(6,*) NHITS, " hits"
            WRITE(6,*) NMISSES, " misses"
            WRITE(6,"(F6.2,A)") (NHITS/(NHITS+NMISSES+0.D0))*100,"% success"
         ENDIF
      END

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
      END

!The cache consists of an unordered set of labels and elements.
!We must order this, and then distribute the elements throughout each set of SLOTS.
      SUBROUTINE FillUpCache()
         IMPLICIT NONE
         INTEGER I,J,K,N,nK
         DO I=1,nPairs
!Find the last value in the cache
            CALL BinarySearch(nPairs+1,UMatLabels(1,I),1,nSlots,N,J,K)
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
      END

!   A binary search to find VAL in TAB.  TAB is sorted, but can have
!   multiple entries being the same.  If the search terminated unsuccessfully, 
!   the entry indicated is one after half-way through the set of entries which 
!   would be immediately prior to it.  From here until the label changes
!   should be filled with VAL if it is to be entered into the table.
!   A and B are the limits of the table.
      SUBROUTINE BINARYSEARCH(VAL,TAB,A,B,LOC,LOC1,LOC2)
         IMPLICIT NONE
         INTEGER VAL,A,B,LOC,LOC1,LOC2
         INTEGER TAB(B)
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
      END

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
      END


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
      END
      SUBROUTINE SWAP(A,B)
         IMPLICIT NONE
         INTEGER A,B,C
         C=A
         A=B
         B=C
         RETURN
      END



!  We're in the middle of freezing some orbitals.
!  OrbTrans(i) will give us the new position of the old orbital i.
      SUBROUTINE FreezeUMatCache(OrbTrans,nOld,nNew)
         IMPLICIT NONE
         INTEGER nOld,nNew,OrbTrans(nOld)
         INTEGER onSlots,onPairs
         INTEGER I,J
         onSlots=nSlots
         onPairs=nPairs
         CALL FreezeUMatCacheInt(OrbTrans,nOld,nNew,onSlots,onPairs)
      END
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
                   CALL CacheUMatEl(A,B,OUMatCacheData(0,n,m),nm,nn)
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
      END
END MODULE UMatCache

      SUBROUTINE GTID(NBASISMAX,GIND,ID)
         USE UMatCache
         IMPLICIT NONE
         INTEGER GIND,NBASISMAX(5,2),ID
            IF(NBASISMAX(2,3).GT.0) THEN
               ID=(GIND-1)/NBASISMAX(2,3)+1
            ELSE
               ID=(GIND-1)/2+1
               IF(TTRANSGTID) ID=TRANSTABLE(ID)
            ENDIF
      RETURN
      END

!  CacheUMatEl and GetCachedUMatEl are outside the module so that they can interact with CPMD
!      REAL*8  FUNCTION WR
!         WRITE(70,*) UMatCacheData(0,2,7)
!         CALL FLUSH(70)
!         WR=sq(UMatCacheData(0,2,7))
!      END
!     Set an element in the cache.  All the work has been done for us before
!   as the element we have to set is in (ICACHEI,ICACHE)
!   We still need to fill out the space before or after  us if we've been put in the 
!   middle of a block of duplicates
      SUBROUTINE CACHEUMATEL(A,B,UMATEL,ICACHE,ICACHEI)
         USE HElement
         USE UMatCache
         IMPLICIT NONE
         INTEGER A,B,ICACHE,ICACHEI
         TYPE(HElement) UMATEL(0:NTYPES-1)
         INTEGER OLAB,IC1,I,J,ITOTAL
         SAVE ITOTAL
         DATA ITOTAL /0/
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
      END
!   Lookup in the cache to see if there's a stored element.  If not, return TRUE.
!    This will rearrange IDI,IDJ,IDK,IDL into the correctkorder
!   (i,k)<=(j,l) and i<=k, j<=l.  ICACHE corresponds to the pair (i,j), and
!   ICACHEI is the index in that cache where the cache should be located.
      LOGICAL FUNCTION GETCACHEDUMATEL(IDI,IDJ,IDK,IDL,UMATEL,ICACHE,ICACHEI,A,B,ITYPE)
         USE HElement
         USE UMatCache
         IMPLICIT NONE
         INTEGER IDI,IDJ,IDK,IDL,ICACHE,ICACHEI
         INTEGER ICACHEI1,ICACHEI2
         TYPE(HElement) UMATEL
         INTEGER I,A,B,ITYPE,ISTAR,ISWAP
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
!         WRITE(6,"(6I3)") IDI,IDJ,IDK,IDL,A,B
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
!  abcd   -> badc <>
!->cbad *.-> cdab ** -> dcba **<>
!         -> bcda *.<>
!->adcb .*-> dabc .*<>

               IF(ISTAR.EQ.1) THEN
!  If we star the first pair, that corresponds to the plain TYPE 1
                  ITYPE=1
               ELSEIF(ISTAR.EQ.2) THEN
!  If we star the second pair, that corresponds to TYPE 1.
!  If there's no swap, it's starred, otherwise it's not starred.
                  IF(ISWAP.EQ.0) THEN
                     ITYPE=3
                  ELSE
                     ITYPE=1
                  ENDIF
               ELSEIF(ISTAR.EQ.3) THEN
!  We complex conjg if we've swapped within each of the pairs, setting bit 1 if ISTAR=2 (ISTAR!=1)
                  ITYPE=2
               ENDIF
            ENDIF
         ENDIF
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
            CALL BINARYSEARCH(NPAIRS+1,UMATLABELS(1,A),1,NSLOTS,ICACHEI,ICACHEI1,ICACHEI2)
               ICACHEI=ICACHEI1
               ICACHEI2=ICACHEI1
               IF(UMatLabels(iCacheI,A).NE.0) WRITE(6,*) "Cache Overwrite", A,B
   
!                  WRITE(6,*) IDI,IDJ,IDK,IDL
!                  WRITE(6,*) A,B,NSLOTS,NPAIRS
!                  WRITE(6,*) ICACHEI1,ICACHEI2
!                  WRITE(6,*) ICACHEI
            ELSE
            CALL BINARYSEARCH(B,UMATLABELS(1,A),1,NSLOTS,ICACHEI,ICACHEI1,ICACHEI2)
            ENDIF
         ENDIF
         IF(UMATLABELS(ICACHEI,ICACHE).EQ.B) THEN
!            WRITE(6,*) "C",IDI,IDJ,IDK,IDL,ITYPE,UMATCACHE(0:nTypes-1,ICACHEI,ICACHE)
            UMATEL=UMatCacheData(IAND(ITYPE,1),ICACHEI,ICACHE)
            IF(ITYPE.GT.1) UMATEL=DCONJG(UMATEL)
!   signal success
            GETCACHEDUMATEL=.FALSE.
         ELSE
!   signal failure
            GETCACHEDUMATEL=.TRUE.
!            WRITE(68,*) A,B,ICACHEI1,ICACHEI2,ICACHEI
         ENDIF
         RETURN
      END

      
