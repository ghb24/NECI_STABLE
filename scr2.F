      SUBROUTINE SCR2(NEL,NBASISMAX,NDET1,NDET2,
     &     G1,NHG,UMAT,ALAT,ISS,COULCOUPLE,SUM)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
      IMPLICIT NONE
      INTEGER NDET1(NEL)
      INTEGER NDET2(NEL)
      INTEGER G1(5,*),ISS
      INTEGER ND1(2),ND2(2),IDS(4)
      REAL*8 UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3)
      REAL*8 GETUMATEL
      INTEGER NEL,NHG
      INTEGER NBASISMAX(5,3)
      REAL*8 COULCOUPLE,SUM

      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

C..Here we zero ND1 and ND2
      DO I=1,2
         ND1(I)=0
         ND2(I)=0
      ENDDO
C..IQ ensures that only two differences are recorded
C..This first routine compares NDET1 to NDET2
      IQ=0
      DO I=1,NEL
         K=0
         DO J=1,NEL
            IF(NDET1(I).EQ.NDET2(J)) CONTINUE
            IF(NDET1(I).NE.NDET2(J)) K=K+1
         ENDDO
         IF(K.EQ.NEL) THEN
            IQ=IQ+1
            ND1(IQ)=NDET1(I)
            IF(IQ.GE.2) GOTO 100
         ENDIF
      ENDDO
 100  CONTINUE
C..This routine compares NDET2 to NDET1
      IQ=0
      DO I=1,NEL
         K=0
         DO J=1,NEL
            IF(NDET2(I).EQ.NDET1(J)) CONTINUE
            IF(NDET2(I).NE.NDET1(J)) K=K+1
         ENDDO
         IF(K.EQ.NEL) THEN
            IQ=IQ+1
            ND2(IQ)=NDET2(I)
            IF(IQ.GE.2) GOTO 200
         ENDIF
      ENDDO
 200  CONTINUE
      CALL GTID(NBASISMAX,ND1(1),ID2)
      CALL GTID(NBASISMAX,ND1(2),ID1)
C..
      CALL GTID(NBASISMAX,ND2(1),ID4)
      CALL GTID(NBASISMAX,ND2(2),ID3)
C..Check the spin combinations of N', N'-1, N and N-1 
      CALL IAZZERO(IDS,4)
      CALL GENKRON(G1,NBASISMAX,ND1,ND2,IDS)
C..N'>N'-1 and N>N-1 *See Levine for rule 3*
      SUM=DFLOAT(IDS(1)*IDS(2))*
     &   GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
     & -DFLOAT(IDS(3)*IDS(4))*
     &   GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      IF(SUM.NE.0.D0) CALL CHK_MS(NDET1,NDET2,NBASISMAX,
     &   NEL,3,G1)
      SUM=SUM*COULCOUPLE
      RETURN
      END 
C ==----------------------------------------------------------------==
      SUBROUTINE GENKRON(G1,NBASISMAX,ND1,ND2,IDS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER G1(5,*)
      INTEGER ND1(2),ND2(2),IDS(4)
C..This routine generates all the requird Kronecker deltas
C..See notes for full explanation
      CALL CHK(G1(4,ND1(2)),G1(4,ND2(2)),IDS(1))
      CALL CHK(G1(4,ND1(1)),G1(4,ND2(1)),IDS(2))
      CALL CHK(G1(4,ND1(2)),G1(4,ND2(1)),IDS(3))
      CALL CHK(G1(4,ND1(1)),G1(4,ND2(2)),IDS(4))
C..
      RETURN
      END
