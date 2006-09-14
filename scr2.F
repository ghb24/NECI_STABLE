      SUBROUTINE SCR2(NEL,NBASISMAX,NDET1,NDET2,
     &     G1,NHG,UMAT,ALAT,ISS,COULCOUPLE,SUM)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
      IMPLICIT NONE
      INTEGER NDET1(NEL)
      INTEGER NDET2(NEL)
      INCLUDE 'basis.inc'
      TYPE(BasisFN) G1(*)
      INTEGER ISS
      INTEGER EX(2,2)
      REAL*8 UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3)
      REAL*8 GETUMATEL
      INTEGER NEL,NHG
      INTEGER NBASISMAX(5,3)
      REAL*8 COULCOUPLE,SUM,SUM2
      REAL*8 COULDAMP
      LOGICAL TSIGN
      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

      INTEGER ISUB
!      CALL TISET('      SCR2',ISUB)
      EX(1,1)=2
      CALL GETEXCITATION(NDET1,NDET2,NEL,EX,TSIGN)
      CALL GTID(NBASISMAX,EX(1,1),ID2)
      CALL GTID(NBASISMAX,EX(1,2),ID1)
C..
      CALL GTID(NBASISMAX,EX(2,1),ID4)
      CALL GTID(NBASISMAX,EX(2,2),ID3)
      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,2))%Ms) THEN
         SUM=GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
      ELSE
         SUM=0.D0
      ENDIF

      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,2))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,1))%Ms) THEN
         SUM=SUM
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      ENDIF
         
      SUM=SUM*COULCOUPLE
      IF(TSIGN) SUM=-SUM
C     &      COULDAMP(ND1(1),NHG)*COULDAMP(ND1(2),NHG)*
C     &      COULDAMP(ND2(1),NHG)*COULDAMP(ND2(2),NHG)

!      CALL TIHALT('      SCR2',ISUB)
      RETURN
      END 
C ==----------------------------------------------------------------==
