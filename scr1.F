      SUBROUTINE SCR1(NEL,NBASISMAX,NDET1,NDET2,G1,
     &	NHG,UMAT,ALAT,ISS,COULCOUPLE,TEXCH,SUM)
      USE HElem
      use Integrals, only : GetUMatEl
      USE UMatCache, only : GTID
      use OneEInts, only: GetTMatEl
      use SystemData, only: BasisFN
!      use global_utilities
      IMPLICIT NONE
      INTEGER NDET1(NEL),NDET2(NEL)
      TYPE(BasisFN) G1(*)
      INTEGER nBasisMax(5,*)
      INTEGER ID,EX(2)
      TYPE(HElement) UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      TYPE(HElement) SUM
      LOGICAL TSIGN
      INTEGER NEL,NHG,ISS
      REAL*8 COULCOUPLE,CD,ALAT(3)
      REAL*8 COULDAMP

      INTEGER IQ,I,J,K,ID1,ID2
      TYPE(HElement) SUM2,SUM3
      LOGICAL LCOUL,TEXCH
!      type(timer), save :: proc_timer
!      proc_timer%timer_name='      SCR1'
!      call set_timer(proc_timer) 
C.. NB in HF we need to include TMAT element


C..In this routine the determinants differ by only 1 spin orbital
      EX(1)=1
      CALL GETEXCITATION(NDET1,NDET2,NEL,EX,TSIGN)
C.. ND2 is the deifferent basis fn in NDET2
C..The two electron integrals are to be done 
      CALL GTID(NBASISMAX,EX(1),ID1)
      CALL GTID(NBASISMAX,EX(2),ID2)
C      CD=COULDAMP(ND1,NHG)*COULDAMP(ND2,NHG)
C..
      LCOUL=G1(EX(1))%Ms.EQ.G1(EX(2))%Ms
      SUM2=0.D0
      DO I=1,NEL
         IF(EX(1).NE.NDET1(I)) THEN
            CALL GTID(NBASISMAX,NDET1(I),ID)
            IF(LCOUL) THEN 
               SUM3=
     &     GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID,ID2,ID)
!               WRITE(6,*) ID1,ID,ID2,ID,SUM3
            ELSE
               SUM3=0.D0
            ENDIF
            IF(TEXCH.AND.G1(EX(1))%Ms.EQ.G1(NDET1(I))%Ms
     &         .AND.     G1(EX(2))%Ms.EQ.G1(NDET1(I))%Ms) THEN
               SUM3=SUM3-
     &     GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID,ID,ID2)
!               WRITE(6,*) ID1,ID,ID,ID2,SUM3
            ENDIF
            SUM2=SUM2+SUM3
         ENDIF
C*(COULDAMP(NDET1(I),NHG)**2)*CD
      ENDDO
C..  We then add in the non-diagonal part of the kinetic energy -
C..  <ph_a|T|ph_a'> where a and a' are the only basis fns that differ in
C..  NDET1 and NDET2
      SUM=SUM2*HElement(COULCOUPLE)+GetTMATEl(EX(1),EX(2))
!      WRITE(6,*) EX(1),EX(2)
!      WRITE(6,*) SUM2,GetTMATEl(EX(1),EX(2))
      IF(TSIGN) SUM=-SUM
!      call halt_timer(proc_timer) 
      RETURN
      END
      
!Same as above, but calculates the one excit level hamiltonian matrix element
!from the determinants and their excitation matrix.
      
      SUBROUTINE SCR1Excit2(NEL,NBASISMAX,NDET1,NDET2,G1,
     &	NHG,UMAT,ALAT,ISS,COULCOUPLE,TEXCH,SUM,EX,TParity)
      USE HElem
      use Integrals, only : GetUMatEl
!      USE UMatCache, only : GTID
      use OneEInts, only: GetTMatEl
      use SystemData, only: BasisFN
!      use global_utilities
      IMPLICIT NONE
      INTEGER NDET1(NEL),NDET2(NEL)
      TYPE(BasisFN) G1(*)
      INTEGER nBasisMax(5,*)
      INTEGER ID,EX(2,2)
      TYPE(HElement) UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      TYPE(HElement) SUM
      LOGICAL TParity
      INTEGER NEL,NHG,ISS
      REAL*8 COULCOUPLE,CD,ALAT(3)
      REAL*8 COULDAMP

      INTEGER IQ,I,J,K,ID1,ID2
      TYPE(HElement) SUM2,SUM3
      LOGICAL LCOUL,TEXCH
!      type(timer), save :: proc_timer
!      proc_timer%timer_name='      SCR1'
!      call set_timer(proc_timer) 
C.. NB in HF we need to include TMAT element


C..In this routine the determinants differ by only 1 spin orbital
!      EX(1)=1
!      CALL GETEXCITATION(NDET1,NDET2,NEL,EX,TParity)
C.. ND2 is the deifferent basis fn in NDET2
C..The two electron integrals are to be done 
!      CALL GTID(NBASISMAX,EX(1,1),ID1)
!      CALL GTID(NBASISMAX,EX(2,1),ID2)
!Explicitly inline GTID here.

       IF(NBASISMAX(2,3).gt.0) THEN
           ID1=(EX(1,1)-1)/NBASISMAX(2,3)+1
           ID2=(EX(2,1)-1)/NBASISMAX(2,3)+1
       ELSE
           ID1=(EX(1,1)-1)/2+1
           ID2=(EX(2,1)-1)/2+1
!BEWARE that the tTransGTID option is not here.
       ENDIF


C      CD=COULDAMP(ND1,NHG)*COULDAMP(ND2,NHG)
C..
      LCOUL=G1(EX(1,1))%Ms.EQ.G1(EX(2,1))%Ms
      SUM2=0.D0
      DO I=1,NEL
         IF(EX(1,1).NE.NDET1(I)) THEN
!Again, explicitly inline gtid.
            IF(NBASISMAX(2,3).gt.0) THEN
                ID=(NDET1(I)-1)/NBASISMAX(2,3)+1
            ELSE
                ID=(NDET1(I)-1)/2+1
            ENDIF
!            CALL GTID(NBASISMAX,NDET1(I),ID)
            IF(LCOUL) THEN 
               SUM3=
     &     GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID,ID2,ID)
            ELSE
               SUM3=0.D0
            ENDIF
            IF(TEXCH.AND.G1(EX(1,1))%Ms.EQ.G1(NDET1(I))%Ms
     &         .AND.     G1(EX(2,1))%Ms.EQ.G1(NDET1(I))%Ms) THEN
               SUM3=SUM3-
     &     GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID,ID,ID2)
            ENDIF
            SUM2=SUM2+SUM3
         ENDIF
C*(COULDAMP(NDET1(I),NHG)**2)*CD
      ENDDO
C..  We then add in the non-diagonal part of the kinetic energy -
C..  <ph_a|T|ph_a'> where a and a' are the only basis fns that differ in
C..  NDET1 and NDET2
!      WRITE(6,*) EX(1,1),EX(2,1)
!      WRITE(6,*) SUM2,GetTMATEl(EX(1,1),EX(2,1))
      SUM=SUM2*HElement(COULCOUPLE)+GetTMATEl(EX(1,1),EX(2,1))
      IF(TParity) SUM=-SUM
!      call halt_timer(proc_timer) 
      RETURN
      END
