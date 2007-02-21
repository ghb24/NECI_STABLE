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
      SUBROUTINE GetUMatSize(nBasis,iSS,iSize)
         IMPLICIT NONE
         INTEGER nBasis,iSS
         INTEGER iPairs,nBi,iSize
         nBi=nBasis/iSS
         iPairs=(nBi*(nBi+1))/2
         iSize=(iPairs*(iPairs+1))/2
      END
      FUNCTION GetUMatEl(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,
     &      IDI,IDJ,IDK,IDL)
         USE HElement
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
         LOGICAL GETCACHEDUMATEL,LSYMSYM
         TYPE(Symmetry) SYM,SYMPROD,SYMCONJ
         INCLUDE 'umatcache.inc'
         INTEGER ISUB,ISUB2
         INTEGER UMatInd
!         CALL TISET(' GETUMATEL',ISUB)
C.. IF NBASISMAX(1,3) is less than zero, we directly give the integral.
C.. Otherwise we just look it up in umat
         IF(NBASISMAX(1,3).GE.0) THEN
C.. See if we need to calculate on the fly
          IF(ISS.EQ.0) THEN

C  JSS - store <ij|ij> and <ij|ji> in UMAT2D.
C.. Remember permutations.  Complex case is the same as real.
C.. <ij|ji> = <ii|jj>

!.. Complex case is more difficult.

!  <ii|ii> is always real and allowed
!  <ij|ij> is always real and allowed (densities i*i and j*j)
!  <ij|ji> is always real and allowed (codensities i*j and j*i = (i*j)*)
!  <ii|jj> is not stored in UMAT2D, and may not be allowed by symmetry.  It can be complex.
           IF (IDI.eq.IDJ.and.IDI.eq.IDK.and.IDI.eq.IDL
     &      .AND.TUMAT2D) THEN
C..  <ii|ii>
             GETUMATEL=UMAT2D(IDI,IDI)
           ELSE IF (IDI.eq.IDK.and.IDJ.eq.IDL
     &      .AND.TUMAT2D) THEN
C.. <ij|ij>
             I=MIN(IDI,IDJ)
             J=MAX(IDI,IDJ)
             GETUMATEL=UMAT2D(I,J)
           ELSE IF (IDI.eq.IDL.and.IDJ.eq.IDK
     &      .AND.TUMAT2D) THEN
C.. <ij|ji>
             I=MAX(IDI,IDJ)
             J=MIN(IDI,IDJ)
             GETUMATEL=UMAT2D(I,J)
           ELSE IF (IDI.eq.IDJ.and.IDK.eq.IDL
     &      .AND.TUMAT2D.AND.HElementSize.EQ.1) THEN
C.. <ii|jj> = <ij|ji> Only for real systems
             I=MAX(IDI,IDK)
             J=MIN(IDI,IDK)
             GETUMATEL=UMAT2D(I,J)
           ELSE
C.. Check to see if the umat element is in the cache
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
C.. Check the symmetry of the 4-index integrals
              IF(.NOT.LSYMSYM(SYM)) THEN
                  GETUMATEL=0.D0
! JSS --- comment out the following TIHALT.  (TISET commented out above.)
!                 CALL TIHALT(' GETUMATEL',ISUB)
                  RETURN
              ELSE
               
C..  This will rearrange I,J,K,L into the correct order
C.. (i,k)<=(j,l) and i<=k, j<=l.
               IF(GETCACHEDUMATEL(I,J,K,L,GETUMATEL,ICACHE,ICACHEI,A,B,
     &            ITYPE)) THEN
C.. We don't have a stored UMAT - we call to generate it.
                  IF(IP_DFInts.NE.0) THEN
C.. We're using density fitting
                     Call GetDF2EInt(I,J,K,L,GETUMATEL)
                  ELSE
C.. Otherwise we call CPMD
                     IF(TTRANSFINDX) THEN
!         WRITE(6,"(A,4I5,$)") "MM",TRANSTABLE(I),TRANSTABLE(J),
!     &                     TRANSTABLE(K),TRANSTABLE(L)
                        CALL INITFINDXI(TRANSTABLE(I),TRANSTABLE(J),
     &                     TRANSTABLE(K),TRANSTABLE(L),UElems)
                     ELSE
                        CALL INITFINDXI(I,J,K,L,UElems)
                     ENDIF
                     GETUMATEL=UElems(IAND(ITYPE,1))
                     IF(ITYPE.GT.1) GETUMATEL=DCONJG(GETUMATEL)
!                     WRITE(6,*) "I",I,J,K,L,UElems
!                     IF(TTRANSFINDX) THEN
!                        CALL INITFINDXI(TRANSTABLE(K),TRANSTABLE(J),
!     &                     TRANSTABLE(I),TRANSTABLE(L),UElems)
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
                  IF(ICACHE.NE.0) CALL CACHEUMATEL(A,B,UElems,ICACHE,
     &                  ICACHEI)
                  NMISSES=NMISSES+1
!                  WRITE(6,*) "MISS",I,J,K,L,A,B,GETUMATEL
               ELSE
C                  WRITE(6,*) "HIT ",IDI,IDJ,IDK,IDL,GETUMATEL
C                  WRITE(6,*) A,B
                  NHITS=NHITS+1
              ENDIF
             ENDIF
           ENDIF
          ELSE
!         CALL TIHALT(' GETUMATEL',ISUB)
!         CALL TISET('GETUMATEL2',ISUB2)
             GETUMATEL=UMAT(UMatInd(IDI,IDJ,IDK,IDL))
!         CALL TIHALT('GETUMATEL2',ISUB2)
          ENDIF
!          CALL TIHALT(' GETUMATEL',ISUB)
         ELSEIF(NBASISMAX(1,3).EQ.-1) THEN
            I=(IDI-1)*ISS+1
            J=(IDJ-1)*ISS+1
            K=(IDK-1)*ISS+1
            L=(IDL-1)*ISS+1
C.. The Uniform electron gas
                   A=G1(I)%k(1)-G1(K)%k(1)
                   B=G1(I)%k(2)-G1(K)%k(2)
                   C=G1(I)%k(3)-G1(K)%k(3)
                   IF(    (G1(L)%k(1)-G1(J)%k(1)).EQ.A
     &               .AND.(G1(L)%k(2)-G1(J)%k(2)).EQ.B
     &               .AND.(G1(L)%k(3)-G1(J)%k(3)).EQ.C
     &               .AND.((A.NE.0).OR.(B.NE.0).OR.(C.NE.0))) THEN
C                     WRITE(6,*) "(",I,J,"|",K,L,")",A,B,C

CC  AJWT  <IJ|r_12^-1|KL> = v_(G_I-G_K) delta_((G_I-G_K)-(G_L-G_J)
CC  v_G = 4 Pi/ G**2.  G=2 Pi/L(nx,ny,nx) etc.
                     SUM=((A/ALAT(1))**2+(B/ALAT(2))**2)
                     IF(ALAT(3).NE.0.D0) SUM=SUM+(C/ALAT(3))**2
                     SUM=1/(PI*SUM*ALAT(1)*ALAT(2)*ALAT(3))
                   ELSE
                     SUM=0.D0
                   ENDIF
            GETUMATEL=SUM
         ENDIF
!         WRITE(6,"(4I5,$)") IDI,IDJ,IDK,IDL
!         WRITE(6,*) GETUMATEL,ABS(GETUMATEL)
         RETURN
!         CALL TIHALT(' GETUMATEL',ISUB)
      END


      SUBROUTINE GTID(NBASISMAX,GIND,ID)
         USE HElement
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER GIND,NBASISMAX(5,2),ID
            IF(NBASISMAX(2,3).GT.0) THEN
               ID=(GIND-1)/NBASISMAX(2,3)+1
            ELSE
               ID=(GIND-1)/2+1
               IF(TTRANSGTID) ID=TRANSTABLE(ID)
            ENDIF
      RETURN
      END

      SUBROUTINE SETUPUMATCACHE(NSTATE)
         USE HElement
         IMPLICIT NONE
         INTEGER NSTATE
         INCLUDE 'umatcache.inc'
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
c        NSLOTS=MIN(NPAIRS, 1024)
            NSLOTS=MIN(NPAIRS, NSLOTSINIT)
            UMATCACHEFLAG=0
            WRITE(6,"(A,I3,2I7,I10)") "UMAT NTYPES,NSLOTS,NPAIRS,TOT",
     &         NTYPES,NSLOTS,NPAIRS,NSLOTS*NPAIRS*NTYPES
            TUMAT2D=.FALSE.
            CALL MEMORY(IP_UMATCACHE,
     &         NTYPES*HElementSize*NSLOTS*NPAIRS,'UMATCACHE')
            CALL MEMORY(IP_UMATLABELS,NSLOTS*NPAIRS/IRAT+1,'UMATLABELS')
            CALL AZZERO(UMATCACHE,NTYPES*HElementSize*NPAIRS*NSLOTS)
            CALL IAZZERO(UMATLABELS,NPAIRS*NSLOTS)
         ENDIF
      END

      SUBROUTINE SETUPUMAT2D(G1)
         USE HElement
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER G1(*)
         IF(NSLOTSINIT.LT.0) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            CALL MEMORY(IP_UMAT2D,HElementSize*NSTATES*NSTATES,'UMAT2D')
            CALL CPMDANTISYMINTEL(G1,UMAT2D,NSTATES)
         ENDIF
      END


      SUBROUTINE SETUPUMAT2D_DF()
         USE HElement
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         IF(NSLOTSINIT.LT.0) THEN
            TUMAT2D=.FALSE.
            WRITE(6,*) "Not using UMAT2D."
         ELSE
            TUMAT2D=.TRUE.
            CALL MEMORY(IP_UMAT2D,HElementSize*NSTATES*NSTATES,'UMAT2D')
            CALL ReadDalton2EIntegrals(nStates,UMat2D)
         ENDIF
      END

     
      SUBROUTINE SETUMATTRANS(TRANS)
         USE HElement
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER TRANS(NSTATES)
         CALL MEMORY(IP_TRANSTABLE,NSTATES,'TRANSTABLE')
         CALL ICOPY(NSTATES,TRANS,1,TRANSTABLE,1)
         TTRANSGTID=.TRUE.
      END
 
      SUBROUTINE SetupUMatTransTable(OldNew,nOld,nNew)
         USE HElement
         IMPLICIT NONE
         INTEGER nNew,nOld,I
         INTEGER OldNew(*)
         INCLUDE 'umatcache.inc'
         CALL MEMORY(IP_TRANSTABLE,nNew/2,'TRANSTABLE')
         DO I=2,nOld,2
            IF(OldNew(I).NE.0) TransTable(OldNew(I)/2)=I/2
         ENDDO
         Write(6,*) "New->Old State Translation Table"
         DO I=1,nNew/2
            WRITE(6,*) I,TransTable(I)
         ENDDO
         TTRANSFINDX=.TRUE.
      END
      SUBROUTINE DESTROYUMATCACHE
         USE HElement
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         IF(IP_UMATCACHE.NE.0) THEN
            CALL FREEM(IP_UMATCACHE)
            IP_UMATCACHE=0
            CALL FREEM(IP_UMATLABELS)
            IP_UMATLABELS=0
            IF(IP_UMAT2D.NE.0) CALL FREEM(IP_UMAT2D)
            IP_UMAT2D=0
            IF(IP_TRANSTABLE.NE.0) THEN
               CALL FREEM(IP_TRANSTABLE)
               IP_TRANSTABLE=0
            ENDIF
            WRITE(6,*) "UMAT Cache Statistics"
            WRITE(6,*) NHITS, " hits"
            WRITE(6,*) NMISSES, " misses"
            WRITE(6,"(F6.2,A)") (NHITS/(NHITS+NMISSES+0.D0))*100,
     &      "% success"
         ENDIF
      END

      SUBROUTINE SETUMATCACHEFLAG(NEWFLAG)
         USE HElement
         IMPLICIT NONE
         INTEGER NEWFLAG,NF
         INCLUDE 'umatcache.inc'
         SELECT CASE(UMATCACHEFLAG)
         CASE(1)
!  We were in direct cache mode where values were distributed correctly throughout the cache.
            IF(NEWFLAG.EQ.0) THEN
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
         USE HElement
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER I,J,K,N,nK
         DO I=1,nPairs
!Find the last value in the cache
            CALL BinarySearch(nPairs+1,UMatLabels(1,I),1,nSlots,N,
     &         J,K)
            N=J-1
!  N is now the last element and thus number of elements.
!  Sort according to label
            CALL SortIRN(N,UMatLabels(1,I),UMatCache(0,1,I),
     &         nTypes*HElementSize)
            CALL SortIRN(N,UMatLabels(1,I),UMatCache(0,1,I),
     &         nTypes*HElementSize)
            K=nSlots
!  Now disperse among the whole array, from the end
            DO J=N,1,-1
               nK=(nSlots*(J-1))/N+1
               UMatLabels(nK:K,I)=UMatLabels(J,I)
               DO K=K,nK,-1
                  UMatCache(:,K,I)=UMatCache(:,J,I)
               ENDDO
               K=nK
            ENDDO
         ENDDO
      END

C.. Lookup in the cache to see if there's a stored element.  If not, return TRUE.
C..  This will rearrange IDI,IDJ,IDK,IDL into the correctkorder
C.. (i,k)<=(j,l) and i<=k, j<=l.  ICACHE corresponds to the pair (i,j), and
C.. ICACHEI is the index in that cache where the cache should be located.
      LOGICAL FUNCTION GETCACHEDUMATEL(IDI,IDJ,IDK,IDL,UMATEL,ICACHE,
     &      ICACHEI,A,B,ITYPE)
         USE HElement
         IMPLICIT NONE
         INTEGER IDI,IDJ,IDK,IDL,ICACHE,ICACHEI
         INTEGER ICACHEI1,ICACHEI2
         TYPE(HElement) UMATEL
         INCLUDE 'umatcache.inc'      
         INTEGER I,A,B,ITYPE,ISTAR,ISWAP
         IF(NSLOTS.EQ.0) THEN
!We don't have a cache so signal failure.
            GETCACHEDUMATEL=.TRUE.
            ICACHE=0
            ITYPE=0
            RETURN
         ENDIF
C.. First ensure the indices are in the correct order
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
         ICACHE=A
         IF(NSLOTS.EQ.NPAIRS) THEN
C.. we've a small enough system to store everything.
            ICACHEI=B
         ELSE
            IF(UMATCACHEFLAG.EQ.1) THEN
!  UMATCACHEFLAG=1 means we are storing a sequence of cache elements,
!  in a blank cache
!  We store them linearly in the cache, and distribute them around later

!Find the last value in the cache
            CALL BINARYSEARCH(NPAIRS+1,UMATLABELS(1,A),1,NSLOTS,ICACHEI,
     &         ICACHEI1,ICACHEI2)
               ICACHEI=ICACHEI1
               ICACHEI2=ICACHEI1
               IF(UMatLabels(iCacheI,A).NE.0)
     &            WRITE(6,*) "Cache Overwrite", A,B
   
!                  WRITE(6,*) IDI,IDJ,IDK,IDL
!                  WRITE(6,*) A,B,NSLOTS,NPAIRS
!                  WRITE(6,*) ICACHEI1,ICACHEI2
!                  WRITE(6,*) ICACHEI
            ELSE
            CALL BINARYSEARCH(B,UMATLABELS(1,A),1,NSLOTS,ICACHEI,
     &         ICACHEI1,ICACHEI2)
            ENDIF
         ENDIF
         IF(UMATLABELS(ICACHEI,ICACHE).EQ.B) THEN
!            WRITE(6,*) "C",IDI,IDJ,IDK,IDL,ITYPE,
!     &         UMATCACHE(0:nTypes-1,ICACHEI,ICACHE)
            UMATEL=UMATCACHE(IAND(ITYPE,1),ICACHEI,ICACHE)
            IF(ITYPE.GT.1) UMATEL=DCONJG(UMATEL)
C.. signal success
            GETCACHEDUMATEL=.FALSE.
         ELSE
C.. signal failure
            GETCACHEDUMATEL=.TRUE.
!            WRITE(68,*) A,B,ICACHEI1,ICACHEI2,ICACHEI
         ENDIF
         RETURN
      END

C.. A binary search to find VAL in TAB.  TAB is sorted, but can have
C.. multiple entries being the same.  If the search terminated unsuccessfully, 
C.. the entry indicated is one after half-way through the set of entries which 
C.. would be immediately prior to it.  From here until the label changes
C.. should be filled with VAL if it is to be entered into the table.
C.. A and B are the limits of the table.
      SUBROUTINE BINARYSEARCH(VAL,TAB,A,B,LOC,LOC1,LOC2)
         IMPLICIT NONE
         INTEGER VAL,A,B,LOC,LOC1,LOC2
         INTEGER TAB(B)
         INTEGER I,J,IFIRST,N,ILAST
C         DO I=A,B
C            WRITE(6,*) I,TAB(I)
C         ENDDO
         I=A
         J=B
         IFIRST=I
         ILAST=J
         DO WHILE(J-I.GE.1)
            N=(I+J)/2
!            WRITE(6,"(A,5I3)") "TN",I,J,N,TAB(N),VAL
            IF(TAB(N).LT.VAL.AND.TAB(N).NE.0.AND.I.NE.N) THEN
               IF(TAB(N).NE.TAB(IFIRST)) IFIRST=N
C.. reset the lower limit
               I=N
            ELSEIF(TAB(N).GT.VAL.OR.TAB(N).EQ.0) THEN
               IF(TAB(N).NE.TAB(ILAST)) ILAST=N
C.. reset the upper limit
               J=N
            ELSEIF(TAB(N).EQ.VAL) THEN
C.. bingo, we've got it!
               LOC=N
C         DO I=A,B
C            WRITE(6,*) I,TAB(I),I.EQ.LOC
C         ENDDO
               LOC1=N
               LOC2=N
               RETURN
            ELSE
C.. we've reached a situation where I and J's entries have the same value, and it's
C.. not the one we want.  Leave the loop.
               I=J
            ENDIF
         ENDDO
C.. We've failed.  However, the new value should sit between I and J.
C.. Split whichever of the prior or after slots which has the most duplicates
C         WRITE(6,*) "FAIL:",IFIRST,I,J,ILAST
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
C         DO I=A,B
C            WRITE(6,*) I,TAB(I),I.EQ.LOC
C         ENDDO
      END

C.. Get a unique index corresponding to pair (I,J), and return in RET.
C.. I<=J<=N.  12/5/06
!  Old Scheme
C.. e.g. 11 12 13 14 15 corresponds to    1  2  3  4  5
C..         22 23 24 25                      6  7  8  9
C..            33 34 35                        10 11 12
C..               44 45                           13 14
C..                  55                              15
!  New Scheme 1/2/07
C.. e.g. 11 12 13 14 15 corresponds to    1  2  4  7 11
C..         22 23 24 25                      3  5  8 12
C..            33 34 35                         6  9 13
C..               44 45                           10 14
C..                  55                              15
      SUBROUTINE GETCACHEINDEX(I,J,N,RET)
         IMPLICIT NONE
         INTEGER I,J,RET,N
         RET=J*(J-1)/2+I
!         RET=N*(I-1)-I*(I-1)/2+J
      END


!      2n              (n)   sqrt(2n)
C..   2  4  8 14 22... 92  -> 1  2  2  3  4 ... 9
C..      6 10 16 24    94        2  3  4  4     9
C..        12 18 26    96           3  4  5     9
C..           20 28    98              4  5     9
C..              30   100                 5    10
C..                   102
c..                   104
c..                   106
c..                   108
c..                   110                      10
C.. return I<=J
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

C..   Set an element in the cache.  All the work has been done for us before
C.. as the element we have to set is in (ICACHEI,ICACHE)
C.. We still need to fill out the space before or after  us if we've been put in the 
C.. middle of a block of duplicates
      SUBROUTINE CACHEUMATEL(A,B,UMATEL,ICACHE,ICACHEI)
         USE HElement
         IMPLICIT NONE
         INTEGER A,B,ICACHE,ICACHEI
         TYPE(HElement) UMATEL(0:NTYPES-1)
         INCLUDE 'umatcache.inc'      
         INTEGER OLAB,IC1
!         WRITE(69,*) "CACHE",A,B,UMATEL,ICACHE,ICACHEI
         IF(NSLOTS.EQ.NPAIRS.OR.UMATCACHEFLAG.EQ.1) THEN
C.. small system.  only store a single element
            UMATLABELS(ICACHEI,ICACHE)=B
            UMATCACHE(:,ICACHEI,ICACHE)=UMATEL
!            WRITE(6,*),
            RETURN
         ENDIF
         IC1=ICACHEI
C         WRITE(6,*) "ICI",ICACHEI,ICACHE
         OLAB=UMATLABELS(ICACHEI,ICACHE)
C.. If we're in a block of prior, fill after
         DO WHILE(OLAB.LT.B.AND.ICACHEI.LE.NSLOTS)
            UMATCACHE(:,ICACHEI,ICACHE)=UMATEL
            UMATLABELS(ICACHEI,ICACHE)=B
C            IF(ICACHEI.LT.1.OR.ICACHE.LT.1
C     &         .OR.ICACHEI.GT.NSLOTS.OR.ICACHE.GT.NPAIRS) THEN
C               WRITE(6,*) ICACHEI,ICACHE
C               STOP "a"
C            ENDIF
            ICACHEI=ICACHEI+1
            IF(ICACHEI.LE.NSLOTS) THEN
               OLAB=UMATLABELS(ICACHEI,ICACHE)
            ELSE
               OLAB=0
            ENDIF
         ENDDO
         IF(OLAB.EQ.0) ICACHEI=IC1
C        WRITE(6,*) "ICI2",ICACHEI,ICACHE
         DO WHILE((OLAB.GT.B.OR.OLAB.EQ.0).AND.ICACHEI.GT.0)
            UMATCACHE(:,ICACHEI,ICACHE)=UMATEL
            UMATLABELS(ICACHEI,ICACHE)=B
C            IF(ICACHEI.LT.1.OR.ICACHE.LT.1
C     &         .OR.ICACHEI.GT.NSLOTS.OR.ICACHE.GT.NPAIRS) THEN
C               WRITE(6,*) ICACHEI,ICACHE
C               STOP "b"
C            ENDIF
            ICACHEI=ICACHEI-1
            OLAB=UMATLABELS(ICACHEI,ICACHE)
         ENDDO
      END


!  We're in the middle of freezing some orbitals.
!  OrbTrans(i) will give us the new position of the old orbital i.
      SUBROUTINE FreezeUMatCache(OrbTrans,nOld,nNew)
         USE HElement
         IMPLICIT NONE
         INTEGER nOld,nNew,OrbTrans(nOld)
         INCLUDE 'umatcache.inc'
         TYPE(HElement) NUMat2D(nNew/2,nNew/2),El(0:nTypes-1)
         POINTER (IP_NUMat2D,NUMat2D)
         INTEGER i,j,k,l,m,n
         INTEGER ni,nj,nk,nl,nm,nn,A,B,iType
         TYPE(HElement) OUMatCache(0:nTypes-1,nSlots,nPairs)
         INTEGER OUMatLabels(nSlots,nPairs)
         
         POINTER (IP_OUmatCache,OUmatCache)
         POINTER (IP_OUmatLabels,OUmatLabels)
         INTEGER onSlots,onPairs
                  
         CALL MEMORY(IP_NUMat2D,HElementSize*(nNew/2)**2,'NUMat2D')
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
         CALL FREEM(IP_UMat2D)      
         IP_UMat2D=IP_NUMat2D
! Now go through the other cache.
! First save the memory used for it.
         onSlots=nSlots
         onPairs=nPairs
         IP_OUMatCache=IP_UMatCache
         IP_OUMatLabels=IP_UMatLabels
         IP_UMatCache=0
         IP_UMatLabels=0
!Now reinitialize the cache.
         CALL SetupUMatCache(nNew/2)
         TUMAT2D=.TRUE.
         CALL SetUMatcacheFlag(1)
         DO i=1,nOld/2
          IF(OrbTrans(i*2).NE.0) THEN
           DO k=i,nOld/2
            IF(OrbTrans(k*2).NE.0) THEN
             CALL GetCacheIndex(i,k,nOld/2,m)
             DO n=1,onSlots
              IF(n.EQ.1.OR.
     &            (onSlots.EQ.onPairs.AND.OUMatLabels(n,m).NE.0)
     &            .OR.(onSlots.NE.onPairs.AND.
     &            OUMatLabels(n,m).NE.OUMatLabels(n-1,m))) THEN
                ni=OrbTrans(i*2)/2
                nk=OrbTrans(k*2)/2
                IF(onSlots.EQ.onPairs) THEN
                 CALL GetCacheIndexStates(n,nNew/2,j,l)
                ELSE
                 CALL GetCacheIndexStates(OUMatLabels(n,m),nNew/2,j,l)
                ENDIF
                nj=OrbTrans(j*2)/2
                nl=OrbTrans(l*2)/2
                IF(nj.NE.0.AND.nl.NE.0) THEN
!                  WRITE(6,"(A,4I5,A,2I5,$)") "NC",I,J,K,L,"L",m,n
                   CALL GetCachedUmatEl(ni,nj,nk,nl,El,nm,nn,A,B,iType)
!                   WRITE(6,"(A,4I5,I5,$)") "->",NI,NJ,NK,NL,
!     &   OUMatLabels(n,m)
!                   WRITE(6,*) OUMatCache(0,n,m)
                   CALL CacheUMatEl(A,B,OUMatCache(0,n,m),nm,nn)
                ENDIF
              ENDIF
             ENDDO
            ENDIF
           ENDDO
          ENDIF
         ENDDO
         CALL FREEM(IP_OUMatLabels) 
         CALL FREEM(IP_OUMatCache) 
         CALL SetUMatCacheFlag(0)               
      END
      
