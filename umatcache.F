      REAL*8 FUNCTION GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,
     &      IDI,IDJ,IDK,IDL)
         IMPLICIT NONE
         INTEGER NBASISMAX(5,3),I,J,K,L,NHG,ISS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(NHG)
         REAL*8 ALAT(3)
         REAL*8 UMAT(*)
         INTEGER A,B,C
         INTEGER IDI,IDJ,IDK,IDL
         REAL*8 SUM
         PARAMETER PI=3.14159265358979323846264338327950288419716939937510D0
         INTEGER ICACHE,ICACHEI
         LOGICAL GETCACHEDUMATEL,LSYMSYM
         TYPE(Symmetry) SYM,SYMPROD
         INCLUDE 'umatcache.inc'
         INTEGER ISUB,ISUB2
!         CALL TISET(' GETUMATEL',ISUB)
C.. IF NBASISMAX(1,3) is less than zero, we directly give the integral.
C.. Otherwise we just look it up in umat
         IF(NBASISMAX(1,3).GE.0) THEN
C.. See if we need to calculate on the fly
          IF(ISS.EQ.0) THEN

C  JSS - store <ij|ij> and <ij|ji> in UMAT2D.
C.. Remember permutations
C.. <ij|ji> = <ii|jj>
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
     &      .AND.TUMAT2D) THEN
C.. <ii|jj> = <ij|ji>
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
               SYM=SYMPROD(SYM,G1(I*2-1)%Sym)
               SYM=SYMPROD(SYM,G1(J*2-1)%Sym)
               SYM=SYMPROD(SYM,G1(K*2-1)%Sym)
               SYM=SYMPROD(SYM,G1(L*2-1)%Sym)
C.. Check the symmetry of the 4-index integrals
               IF(.NOT.LSYMSYM(SYM)) THEN
                  GETUMATEL=0.D0
                  CALL TIHALT(' GETUMATEL',ISUB)
                  RETURN
               ENDIF
               
C..  This will rearrange I,J,K,L into the correct order
C.. (i,k)<=(j,l) and i<=k, j<=l.
               IF(GETCACHEDUMATEL(I,J,K,L,GETUMATEL,ICACHE,ICACHEI,A,B))
     &             THEN
C.. We don't have a stored UMAT - we call to generate it.
                  IF(IP_DFInts.NE.0) THEN
C.. We're using density fitting
                     Call GetDF2EInt(I,J,K,L,GETUMATEL)
                  ELSE
C.. Otherwise we call CPMD
                     CALL INITFINDXI(IDI,IDJ,IDK,IDL,GETUMATEL)
                  ENDIF
                  CALL CACHEUMATEL(A,B,GETUMATEL,ICACHE,
     &                  ICACHEI)
                  NMISSES=NMISSES+1
                  WRITE(63,*) I,J,K,L,GETUMATEL
               ELSE
C                  WRITE(6,*) "HIT ",IDI,IDJ,IDK,IDL,GETUMATEL
C                  WRITE(6,*) A,B
                  NHITS=NHITS+1
              ENDIF
           ENDIF
          ELSE
!         CALL TIHALT(' GETUMATEL',ISUB)
!         CALL TISET('GETUMATEL2',ISUB2)
             GETUMATEL=UMAT(IDI+(NHG/ISS)*(IDJ-1+
     &                     (NHG/ISS)*(IDK-1+(NHG/ISS)*(IDL-1))))
!         CALL TIHALT('GETUMATEL2',ISUB2)
         RETURN
          ENDIF
!          CALL TIHALT(' GETUMATEL',ISUB)
          RETURN
         END IF   
         IF(NBASISMAX(1,3).EQ.-1) THEN
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
!         CALL TIHALT(' GETUMATEL',ISUB)
      END


      SUBROUTINE GTID(NBASISMAX,GIND,ID)
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER GIND,NBASISMAX(5,2),ID
            IF(NBASISMAX(2,3).GT.0) THEN
               ID=(GIND-1)/NBASISMAX(2,3)+1
            ELSE
               ID=(GIND-1)/2+1
               IF(IP_TRANSTABLE.NE.0) ID=TRANSTABLE(ID)
            ENDIF
      RETURN
      END

      SUBROUTINE SETUPUMATCACHE(NSTATE)
         IMPLICIT NONE
         INTEGER NSTATE
         INCLUDE 'umatcache.inc'
         INCLUDE 'irat.inc'
         NSTATES=NSTATE
         NPAIRS=NSTATES*(NSTATES+1)/2
c        NSLOTS=MIN(NPAIRS, 1024)
         NSLOTS=MIN(NPAIRS, NSLOTSINIT)
         UMATCACHEFLAG=0
         NHITS=0
         NMISSES=0
         WRITE(6,*) "UMAT NSLOTS,NPAIRS,TOT",NSLOTS,NPAIRS,NSLOTS*NPAIRS
         TUMAT2D=.FALSE.
         CALL MEMORY(IP_UMATCACHE,NSLOTS*NPAIRS,'UMATCACHE')
         CALL MEMORY(IP_UMATLABELS,NSLOTS*NPAIRS/IRAT+1,'UMATLABELS')
         CALL AZZERO(UMATCACHE,NPAIRS*NSLOTS)
         CALL IAZZERO(UMATLABELS,NPAIRS*NSLOTS)
      END

      SUBROUTINE SETUPUMAT2D(G1)
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER G1(*)
         TUMAT2D=.TRUE.
         CALL MEMORY(IP_UMAT2D,NSTATES*NSTATES,'UMAT2D')
         CALL CPMDANTISYMINTEL(G1,UMAT2D,NSTATES)
      END
     
      SUBROUTINE SETUMATTRANS(TRANS)
         IMPLICIT NONE
         INCLUDE 'umatcache.inc'
         INTEGER TRANS(NSTATES)
         CALL MEMORY(IP_TRANSTABLE,NSTATES,'TRANSTABLE')
         CALL ICOPY(NSTATES,TRANS,1,TRANSTABLE,1)
      END
 
      SUBROUTINE DESTROYUMATCACHE
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
         IMPLICIT NONE
         INTEGER NEWFLAG
         INCLUDE 'umatcache.inc'
         UMATCACHEFLAG=NEWFLAG
         RETURN
      END

C.. Lookup in the cache to see if there's a stored element.  If not, return TRUE.
C..  This will rearrange IDI,IDJ,IDK,IDL into the correctkorder
C.. (i,k)<=(j,l) and i<=k, j<=l.  ICACHE corresponds to the pair (i,j), and
C.. ICACHEI is the index in that cache where the cache should be located.
      LOGICAL FUNCTION GETCACHEDUMATEL(IDI,IDJ,IDK,IDL,UMATEL,ICACHE,
     &      ICACHEI,A,B)
         IMPLICIT NONE
         INTEGER IDI,IDJ,IDK,IDL,ICACHE,ICACHEI
         INTEGER ICACHEI1,ICACHEI2
         REAL*8 UMATEL
         INCLUDE 'umatcache.inc'      
         INTEGER I,A,B
C.. First ensure the indices are in the correct order
         IF(IDK.LT.IDI) CALL SWAP(IDI,IDK)
         IF(IDL.LT.IDJ) CALL SWAP(IDJ,IDL)
         CALL GETCACHEINDEX(IDI,IDK,NSTATES,A)
         CALL GETCACHEINDEX(IDJ,IDL,NSTATES,B)
C         WRITE(6,"(6I3)") IDI,IDJ,IDK,IDL,A,B
         IF(A.GT.B) THEN
            CALL SWAP(A,B)
            CALL SWAP(IDI,IDJ)
            CALL SWAP(IDK,IDL)
         ENDIF
         ICACHE=A
         IF(NSLOTS.EQ.NPAIRS) THEN
C.. we've a small enough system to store everything.
            ICACHEI=B
         ELSE
            CALL BINARYSEARCH(B,UMATLABELS(1,A),1,NSLOTS,ICACHEI,
     &         ICACHEI1,ICACHEI2)
C.. UMATCACHEFLAG=1 means we store at the lowest possible position
            IF(UMATCACHEFLAG.EQ.1) ICACHEI=ICACHEI1
         ENDIF
         IF(UMATLABELS(ICACHEI,ICACHE).EQ.B) THEN
            UMATEL=UMATCACHE(ICACHEI,ICACHE)
C.. signal success
            GETCACHEDUMATEL=.FALSE.
         ELSE
C.. signal failure
            GETCACHEDUMATEL=.TRUE.
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
C            WRITE(6,"(A,5I3)") "TN",I,J,N,TAB(N),VAL
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
C.. e.g. 11 12 13 14 15 corresponds to    1  2  3  4  5
C..         22 23 24 25                      6  7  8  9
C..            33 34 35                        10 11 12
C..               44 45                           13 14
C..                  55                              15
      SUBROUTINE GETCACHEINDEX(I,J,N,RET)
         IMPLICIT NONE
         INTEGER I,J,RET,N
         RET=N*(I-1)-I*(I-1)/2+J
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
         IMPLICIT NONE
         INTEGER A,B,ICACHE,ICACHEI
         REAL*8 UMATEL
         INCLUDE 'umatcache.inc'      
         INTEGER OLAB,IC1
C         WRITE(6,*) "CACHE",A,B,UMATEL,ICACHE,ICACHEI
         IF(NSLOTS.EQ.NPAIRS) THEN
C.. small system.  only store a single element
            UMATLABELS(ICACHEI,ICACHE)=B
            UMATCACHE(ICACHEI,ICACHE)=UMATEL
            RETURN
         ENDIF
         IC1=ICACHEI
C         WRITE(6,*) "ICI",ICACHEI,ICACHE
         OLAB=UMATLABELS(ICACHEI,ICACHE)
C.. If we're in a block of prior, fill after
         DO WHILE(OLAB.LT.B.AND.ICACHEI.LE.NSLOTS)
            UMATCACHE(ICACHEI,ICACHE)=UMATEL
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
            UMATCACHE(ICACHEI,ICACHE)=UMATEL
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
