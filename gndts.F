C.. A recursive version of GNDTS which doesn't require #defines
C.. This sets up variables, and calls the recursion 
      SUBROUTINE GNDTS(NEL,NHG,BRR,NBASISMAX,
     &  NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &  IPARITY,II,IFDET)
         IMPLICIT NONE
         INTEGER NEL,NHG,LMS,II
         INTEGER NMRKS(NEL,*)
         INTEGER BRR(NHG)
         INTEGER G1(5,NHG)
         INTEGER IPARITY(5)
         LOGICAL TCOUNT,TSPN,TPARITY,TPRTY
         INTEGER NBASISMAX(5,2),IFDET,IDETS,ITDETS
C..  This will hold which determinant each electron is on at the moment
         INTEGER IELECS(NEL)
C..  This will hold the spin of each electron
         INTEGER ISPINS(NEL)
         INTEGER KJ(5),I
            
C.. Total number of determined determinants
         II=0
         DO I=1,5
            IF(TPARITY) THEN
               KJ(I)=IPARITY(I)
            ELSE
               KJ(I)=0
            ENDIF
         ENDDO
         IF(TSPN) THEN
            KJ(4)=LMS
         ELSE
            KJ(4)=0 
         ENDIF
C.. Start with the first electron
         CALL GNDTS_R(NEL,NHG,BRR,NBASISMAX,
     &      NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &      IPARITY,II,IELECS,ISPINS,1,KJ)
         IFDET=1
         RETURN
      END

C.. The recursive routine in GNDTS.  IELEC is the current electron
      RECURSIVE SUBROUTINE GNDTS_R(NEL,NHG,BRR,NBASISMAX,
     &  NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &  IPARITY,II,IELECS,ISPINS,IELEC,KJ)
      IMPLICIT NONE
      INTEGER NEL,NHG,LMS,II
      INTEGER NMRKS(NEL,*)
      INTEGER BRR(NHG)
      INTEGER G1(5,NHG)
      INTEGER IPARITY(3)
      LOGICAL TCOUNT,TSPN,TPARITY,TPRTY
      INTEGER NBASISMAX(5,2),IFDET,IDETS,ITDETS

      INTEGER IELECS(NEL)
      INTEGER ISPINS(NEL)

      INTEGER NI(NEL)
      INTEGER IELEC,IST,IEL,MSS,I,IP
C.. This will hold the parity of the basis fn
      INTEGER IPAR(3)
      INTEGER KI(5),KJ(5)
      LOGICAL LCHKSYM
C.. Iterate of the the spin basis functions available for this electron
C.. I_(N+1) = I_N +1 .. NHG
      IST=1
      IF(IELEC.GT.1) IST=IELECS(IELEC-1)+1
      DO IEL=IST,NHG
         IELECS(IELEC)=IEL
         ISPINS(IELEC)=G1(4,BRR(IELECS(IELEC)))
C.. If we're on the last electron
         IF(NEL.EQ.IELEC) THEN
            DO I=1,NEL
               NI(I)=BRR(IELECS(I))
            ENDDO
            CALL SORTI(NEL,NI)
            CALL GETSYM(NI,NEL,G1,NBASISMAX,KI)
            IF(.NOT.TSPN) KI(4)=0
            IF(.NOT.TPARITY) THEN
               DO I=1,3
                  KI(I)=0
               ENDDO
               KI(5)=0
            ENDIF
            IF(LCHKSYM(KI,KJ)) THEN
               II=II+1
C.. If we're generating rather than counting
               IF(.NOT.TCOUNT) THEN
                  CALL ICOPY(NEL,NI,1,NMRKS(1:NEL,II),1)
               ENDIF
            ENDIF
         ELSE
C..If we're not on the last electron
C..We recurse over the next electrons
            CALL GNDTS_R(NEL,NHG,BRR,NBASISMAX,
     &            NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &            IPARITY,II,IELECS,ISPINS,IELEC+1,KJ)
         ENDIF
      ENDDO      

      RETURN
      END

C ==--------------------------------------------------==
      SUBROUTINE PARITY(IPAR,K,TPARITY)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL TPARITY
      IF(IABS(K).NE.1) THEN
         WRITE(6,*) K
         CALL STOPGM('PARITY',' IABS(K) NE 1')
      ENDIF
      IF(MOD(IPAR,2).EQ.0) THEN
         IF(K.EQ.1) THEN
            TPARITY=.TRUE.
         ELSEIF(K.EQ.-1) THEN
            TPARITY=.FALSE.
         ENDIF
      ELSE
         IF(K.EQ.1) THEN
            TPARITY=.FALSE.
         ELSEIF(K.EQ.-1) THEN
            TPARITY=.TRUE.
         ENDIF
      ENDIF
      RETURN
      END

