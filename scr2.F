      SUBROUTINE SCR2(NEL,NBASISMAX,NDET1,NDET2,
     &     G1,NHG,UMAT,ALAT,ISS,COULCOUPLE,TEXCH,hel)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
      USE HElem
      use Integrals, Only : GetUMatEl
      USE UMatCache , Only : GTID
      use SystemData, only: BasisFN
!      use global_utilities
      IMPLICIT NONE
      INTEGER NDET1(NEL)
      INTEGER NDET2(NEL)
      TYPE(BasisFN) G1(*)
      INTEGER ISS
      INTEGER EX(2,2)
      TYPE(HElement) UMAT(*),hel,SUM2
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3)
      INTEGER NEL,NHG
      INTEGER nBasisMax(5,*)
      REAL*8 COULCOUPLE
      REAL*8 COULDAMP
      LOGICAL TSIGN,TEXCH
      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

!      type(timer), save :: proc_timer
!      proc_timer%timer_name='      SCR2'
!      call set_timer(proc_timer)
      EX(1,1)=2
      CALL GETEXCITATION(NDET1,NDET2,NEL,EX,TSIGN)
!      WRITE(6,*) "TSIGN: ",TSIGN
      ID1 = GTID(EX(1,1))
      ID2 = GTID(EX(1,2))
C..
      ID3 = GTID(EX(2,1))
      ID4 = GTID(EX(2,2))
      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,2))%Ms) THEN
         hel=GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
      ELSE
         hel=0.D0
      ENDIF

      IF(TEXCH.AND.G1(EX(1,1))%Ms.EQ.G1(EX(2,2))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,1))%Ms) THEN
         hel=hel
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      ENDIF
 
!      WRITE(6,*) EX(1:2,1:2),hel
        
      hel=hel*HElement(COULCOUPLE)
      IF(TSIGN) hel=-hel
C     &      COULDAMP(ND1(1),NHG)*COULDAMP(ND1(2),NHG)*
C     &      COULDAMP(ND2(1),NHG)*COULDAMP(ND2(2),NHG)

!      call halt_timer(proc_timer)
      RETURN
      END 
C ==----------------------------------------------------------------==

!  Work out a double-excitation matrix element in excitation format
      SUBROUTINE SCR2Excit(NBASISMAX,nJ,
     &     G1,NHG,UMAT,ALAT,ISS,hel)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
C..NB. Only used for tStoreAsExcitations. Parity is not required.
      USE HElem
      use Integrals, Only : GetUMatEl
      USE UMatCache , Only : GTID
      use SystemData, only: BasisFN
!      use global_utilities
      IMPLICIT NONE
      INTEGER nJ(7)
      TYPE(BasisFN) G1(*)
      INTEGER ISS
      TYPE(HElement) UMAT(*),hel,SUM2
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3),COULCOUPLE
      INTEGER NHG,EX(2,2)
      INTEGER nBasisMax(5,*)
      LOGICAL TSIGN,TEXCH
      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

!      type(timer), save :: proc_timer
!      proc_timer%timer_name='      SCR2'
!      call set_timer(proc_timer)
      EX(1,1)=nJ(4)
      EX(1,2)=nJ(5)
      EX(2,1)=nJ(6)
      EX(2,2)=nJ(7)
      TSIGN=.FALSE.
      TEXCH=.TRUE.
      COULCOUPLE=1
      ID1 = GTID(EX(1,1))
      ID2 = GTID(EX(1,2))
C..
      ID3 = GTID(EX(2,1))
      ID4 = GTID(EX(2,2))
      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,2))%Ms) THEN
         hel=GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
      ELSE
         hel=0.D0
      ENDIF

      IF(TEXCH.AND.G1(EX(1,1))%Ms.EQ.G1(EX(2,2))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,1))%Ms) THEN
         hel=hel
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      ENDIF
         
      hel=hel*HElement(COULCOUPLE)
      IF(TSIGN) hel=-hel
C     &      COULDAMP(ND1(1),NHG)*COULDAMP(ND1(2),NHG)*
C     &      COULDAMP(ND2(1),NHG)*COULDAMP(ND2(2),NHG)

!      call halt_timer(proc_timer)
      RETURN
      END 

C ==----------------------------------------------------------------==

!  Work out a double-excitation matrix element in excitation format
! We actually pass in the excitation matrix, rather than the STORE-AS-EXCITATIONS
! determinant string.
      SUBROUTINE SCR2Excit2(NBASISMAX,EX,
     &     G1,NHG,UMAT,ALAT,ISS,hel,TSIGN)
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.
      USE HElem
      use Integrals, Only : GetUMatEl
      USE UMatCache , Only : GTID 
      use SystemData, only: BasisFN
!      use global_utilities
      IMPLICIT NONE
      INTEGER nJ(7)
      TYPE(BasisFN) G1(*)
      INTEGER ISS
      TYPE(HElement) UMAT(*),hel,SUM2
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 ALAT(3)
      INTEGER NHG,EX(2,2)
      INTEGER nBasisMax(5,*)
      LOGICAL TSIGN
      INTEGER I,IQ,J,K,ID1,ID2,ID3,ID4

!      type(timer), save :: proc_timer
!      proc_timer%timer_name='      SCR2'
!      call set_timer(proc_timer)

      ID1 = GTID(EX(1,1))
      ID2 = GTID(EX(1,2))
      ID3 = GTID(EX(2,1))
      ID4 = GTID(EX(2,2))

      IF(    G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,2))%Ms) THEN
         hel=GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID3,ID4)
      ELSE
         hel=0.D0
      ENDIF

      IF(G1(EX(1,1))%Ms.EQ.G1(EX(2,2))%Ms
     &    .AND.G1(EX(1,2))%Ms.EQ.G1(EX(2,1))%Ms) THEN
         hel=hel
     &        -GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID2,ID4,ID3)
      ENDIF
         
!      hel=hel*HElement(COULCOUPLE)
      IF(TSIGN) hel=-hel
C     &      COULDAMP(ND1(1),NHG)*COULDAMP(ND1(2),NHG)*
C     &      COULDAMP(ND2(1),NHG)*COULDAMP(ND2(2),NHG)

!      call halt_timer(proc_timer)
      RETURN
      END 
