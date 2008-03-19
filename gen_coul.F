      SUBROUTINE GEN_COUL(NEL,NBASISMAX,
     &   NHG,G1,NMSH,NMAX,FCK,UMAT,ISS)
      use System, only: BasisFN
      IMPLICIT REAL*8 (A-H,O-Z)
      TYPE(BASISFN) G1(NHG)
      INTEGER NBASISMAX(5,3),ISS
      DIMENSION UMAT(NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
C..Cube arrays
      COMPLEX*16 FCK(NMSH,NMSH,NMSH)
!      COMPLEX*16 ZIA(-NMSH/2:NMSH/2,NMAX,NMAX)
C..This routine generates *ALL* possible combinations of Coulomb integrals
C..stored in the form (u1 u2 | U | u1' u2') = UMAT(n1 n2 n3 n4)
C..
C..This first call calculates the inner integral 
C..The call to SCOUL calculates the outer integral 
C ==--------------------------------------------------------------------==
      CALL TISET('GEN_COUL  ',ISUB)
C ==--------------------------------------------------------------------==
      OPEN(10,FILE='UMAT',STATUS='UNKNOWN')
      II=0
      ISPINSKIP=NBASISMAX(2,3)
      DO I=1,NHG,ISPINSKIP
        DO J=1,NHG,ISPINSKIP
          DO K=1,NHG,ISPINSKIP
            DO L=1,NHG,ISPINSKIP
              SUM=0.D0
C..Original call
C              CALL SLATCOULFOU(G1(1,I),G1(1,J),
C     &		G1(1,K),G1(1,L),NMSH,FCK,NMAX,SUM)
              CALL GTID(NBASISMAX,I,ID1)
              CALL GTID(NBASISMAX,J,ID2)
              CALL GTID(NBASISMAX,K,ID3)
              CALL GTID(NBASISMAX,L,ID4)
              IF(UMAT(ID1,ID2,ID3,ID4).NE.0.D0) GOTO 100
              CALL SLATCOULFOU(G1(I),G1(J),
     &		G1(K),G1(L),NMSH,FCK,NMAX,SUM)
c              SUM=SUM*COUPLE 
C..Original call
              UMAT(ID1,ID2,ID3,ID4)=SUM
C..Symmetries
              UMAT(ID3,ID4,ID1,ID2)=SUM
              UMAT(ID2,ID1,ID4,ID3)=SUM
              UMAT(ID4,ID3,ID2,ID1)=SUM
              IF(SUM.LE.1.D-10) GOTO 100 
              WRITE(10,'(4I7,F19.9)') ID1,ID2,ID3,ID4,SUM
100           CONTINUE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(10)
      WRITE(6,*) ' !!! FINISHED CALCULATING ALL 2E INTEGRALS !!! '
C ==--------------------------------------------------------------------==
      CALL TIHALT('GEN_COUL  ',ISUB)
C ==--------------------------------------------------------------------==
      RETURN
      END
C ==------------------------------------------------------------------==
