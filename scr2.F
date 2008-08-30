      SUBROUTINE SCR2(NEL,NBASISMAX,NDET1,NDET2,
     &     G1,NHG,UMAT,ALAT,ISS,COULCOUPLE,TEXCH,SUM)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
      USE HElem
      use Integrals, Only : GetUMatEl
      USE UMatCache , Only : GTID
      use SystemData, only: BasisFN
      IMPLICIT NONE
      INTEGER NDET1(NEL)
      INTEGER NDET2(NEL)
      TYPE(BasisFN) G1(*)
      INTEGER ISS
      INTEGER EX(2,2)
      TYPE(HElement) UMAT(*),SUM,SUM2
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3)
      INTEGER NEL,NHG
      INTEGER nBasisMax(5,*)
      REAL*8 COULCOUPLE
      REAL*8 COULDAMP
      LOGICAL TSIGN,TEXCH
      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

      INTEGER,SAVE :: ISUB=0
!      call set_timer('      SCR2',ISUB)
      EX(1,1)=2
      CALL GETEXCITATION(NDET1,NDET2,NEL,EX,TSIGN)
      CALL GTID(NBASISMAX,EX(1,1),ID1)
      CALL GTID(NBASISMAX,EX(1,2),ID2)
C..
      CALL GTID(NBASISMAX,EX(2,1),ID3)
      CALL GTID(NBASISMAX,EX(2,2),ID4)
      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,2))%Ms) THEN
         SUM=GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
      ELSE
         SUM=0.D0
      ENDIF

      IF(TEXCH.AND.G1(EX(1,1))%Ms.EQ.G1(EX(2,2))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,1))%Ms) THEN
         SUM=SUM
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      ENDIF
         
      SUM=SUM*HElement(COULCOUPLE)
      IF(TSIGN) SUM=-SUM
C     &      COULDAMP(ND1(1),NHG)*COULDAMP(ND1(2),NHG)*
C     &      COULDAMP(ND2(1),NHG)*COULDAMP(ND2(2),NHG)

!      call halt_timer(ISUB)
      RETURN
      END 
C ==----------------------------------------------------------------==

!  Work out a double-excitation matrix element in excitation format
      SUBROUTINE SCR2Excit(NBASISMAX,nJ,
     &     G1,NHG,UMAT,ALAT,ISS,SUM)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
      USE HElem
      use Integrals, Only : GetUMatEl
      USE UMatCache , Only : GTID
      use SystemData, only: BasisFN
      IMPLICIT NONE
      INTEGER nJ(7)
      TYPE(BasisFN) G1(*)
      INTEGER ISS
      TYPE(HElement) UMAT(*),SUM,SUM2
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3),COULCOUPLE
      INTEGER NHG,EX(2,2)
      INTEGER nBasisMax(5,*)
      LOGICAL TSIGN,TEXCH
      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

      INTEGER,SAVE :: ISUB=0
!      call set_timer('      SCR2',ISUB)
      EX(1,1)=nJ(4)
      EX(1,2)=nJ(5)
      EX(2,1)=nJ(6)
      EX(2,2)=nJ(7)
      TSIGN=.FALSE.
      TEXCH=.TRUE.
      COULCOUPLE=1
      CALL GTID(NBASISMAX,EX(1,1),ID1)
      CALL GTID(NBASISMAX,EX(1,2),ID2)
C..
      CALL GTID(NBASISMAX,EX(2,1),ID3)
      CALL GTID(NBASISMAX,EX(2,2),ID4)
      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,2))%Ms) THEN
         SUM=GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
      ELSE
         SUM=0.D0
      ENDIF

      IF(TEXCH.AND.G1(EX(1,1))%Ms.EQ.G1(EX(2,2))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,1))%Ms) THEN
         SUM=SUM
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      ENDIF
         
      SUM=SUM*HElement(COULCOUPLE)
      IF(TSIGN) SUM=-SUM
C     &      COULDAMP(ND1(1),NHG)*COULDAMP(ND1(2),NHG)*
C     &      COULDAMP(ND2(1),NHG)*COULDAMP(ND2(2),NHG)

!      call halt_timer(ISUB)
      RETURN
      END 
