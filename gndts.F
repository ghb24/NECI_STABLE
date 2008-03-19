C.. A recursive version of GNDTS which doesn't require #defines
C.. This sets up variables, and calls the recursion 
      SUBROUTINE GNDTS(NEL,NHG,BRR,NBASISMAX,
     &  NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &  SymRestrict,II,IFDET)
         use System, only: BasisFN,BasisFNSize
         IMPLICIT NONE
         INTEGER NEL,NHG,LMS,II
         INTEGER NMRKS(NEL,*)
         INTEGER BRR(NHG)
         TYPE(BasisFN) G1(NHG),SymRestrict,kJ
         LOGICAL TCOUNT,TSPN,TPARITY,TPRTY
         INTEGER NBASISMAX(5,2),IFDET,IDETS,ITDETS
C..  This will hold which determinant each electron is on at the moment
         INTEGER IELECS(NEL)
C..  This will hold the spin of each electron
         INTEGER ISPINS(NEL)
         INTEGER I
            
C.. Total number of determined determinants
         II=0
         IF(TPARITY) THEN
            kJ=SymRestrict
         ELSE
            CALL IAZZERO(kJ,BasisFNSize)
         ENDIF
         IF(TSPN) THEN
            kJ%Ms=LMS
         ELSE
            kJ%Ms=0 
         ENDIF
C.. Start with the first electron
         CALL GNDTS_R(NEL,NHG,BRR,NBASISMAX,
     &      NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &      SymRestrict,II,IELECS,ISPINS,1,KJ)
         IFDET=1
         RETURN
      END

C.. The recursive routine in GNDTS.  IELEC is the current electron
      RECURSIVE SUBROUTINE GNDTS_R(NEL,NHG,BRR,NBASISMAX,
     &  NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &  SymRestrict,II,IELECS,ISPINS,IELEC,KJ)
      use System, only: BasisFN
      IMPLICIT NONE
      INTEGER NEL,NHG,LMS,II
      INTEGER NMRKS(NEL,*)
      INTEGER BRR(NHG)
      TYPE(BasisFN) G1(nHG),kI,kJ,SymRestrict
      LOGICAL TCOUNT,TSPN,TPARITY,TPRTY
      INTEGER NBASISMAX(5,2),IFDET,IDETS,ITDETS

      INTEGER IELECS(NEL)
      INTEGER ISPINS(NEL)

      INTEGER NI(NEL)
      INTEGER IELEC,IST,IEL,MSS,I,IP
C.. This will hold the parity of the basis fn
      LOGICAL LCHKSYM
C.. Iterate of the the spin basis functions available for this electron
C.. I_(N+1) = I_N +1 .. NHG
      IST=1
      IF(IELEC.GT.1) IST=IELECS(IELEC-1)+1
      DO IEL=IST,NHG
         IELECS(IELEC)=IEL
         ISPINS(IELEC)=G1(BRR(IELECS(IELEC)))%Ms
C.. If we're on the last electron
         IF(NEL.EQ.IELEC) THEN
            DO I=1,NEL
               NI(I)=BRR(IELECS(I))
            ENDDO
            CALL SORTI(NEL,NI)
            CALL GETSYM(NI,NEL,G1,NBASISMAX,KI)
            IF(.NOT.TPARITY) THEN
               DO I=1,3
                  KI%k(I)=0
               ENDDO
               KI%Sym%s=0 ! Ignore symmetry in generation of determinants.
            ENDIF
            IF(.NOT.TSPN) KI%Ms=0
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
     &            SymRestrict,II,IELECS,ISPINS,IELEC+1,KJ)
         ENDIF
      ENDDO      

      RETURN
      END
