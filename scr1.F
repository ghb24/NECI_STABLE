      SUBROUTINE SCR1(NEL,NBASISMAX,NDET1,NDET2,G1,
     &	NHG,UMAT,ALAT,ISS,TMAT,COULCOUPLE,SUM)
      IMPLICIT NONE
      INTEGER NDET1(NEL),NDET2(NEL)
      INCLUDE 'basis.inc'
      TYPE(BasisFN) G1(*)
      INTEGER NBASISMAX(5,2)
      INTEGER ID(NEL-1),IDS1(NEL-1),IDS2(NEL-1)
      REAL*8 UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 TMAT(NHG,NHG),ALAT(3)
      REAL*8 GETUMATEL
      INTEGER NEL,NHG,ISS
      REAL*8 COULCOUPLE,SUM,CD
      REAL*8 COULDAMP

      INTEGER ND1,ND2,IQ,I,J,K,ID1,ID2,IDS
      REAL*8 SUM2,SUM3
      INTEGER ISUB
!      CALL TISET('      SCR1',ISUB) 
C.. NB in HF we need to include TMAT element


C..In this routine the determinants differ by only 1 spin orbital
C..The integers ND1 and ND2 store the spin orbitals that are different 
C..in each determinant
      ND1=0
      ND2=0
C..Compare NDET1 to NDET2
      IQ=0
      DO I=1,NEL
        K=0
        DO J=1,NEL
          IF(NDET1(I).EQ.NDET2(J)) CONTINUE
          IF(NDET1(I).NE.NDET2(J)) K=K+1
        ENDDO
        IF(K.EQ.NEL) THEN
          IQ=IQ+1
          ND1=NDET1(I)
          IF(IQ.GE.1) GOTO 100
        ENDIF
      ENDDO
100   CONTINUE
C.. ND1 is the basis fn which differs in NDET1 to NDET2
C..Now we compare NDET2 to NDET1
      IQ=0
      DO I=1,NEL
        K=0
        DO J=1,NEL
          IF(NDET2(I).EQ.NDET1(J)) CONTINUE
          IF(NDET2(I).NE.NDET1(J)) K=K+1
        ENDDO
        IF(K.EQ.NEL) THEN
          IQ=IQ+1
          ND2=NDET2(I)
          IF(IQ.GE.1) GOTO 200
        ENDIF
      ENDDO
200   CONTINUE 
C.. ND2 is the deifferent basis fn in NDET2
C..The two electron integrals are to be done 
      IQ=0
      DO I=1,NEL
        IF(NDET1(I).NE.ND1) THEN
          IQ=IQ+1
          CALL GTID(NBASISMAX,NDET1(I),ID(IQ))
          CALL CHK(G1(NDET1(I))%Ms,G1(ND1)%Ms,IDS1(IQ))
          CALL CHK(G1(NDET1(I))%Ms,G1(ND2)%Ms,IDS2(IQ))
          IF(IQ.GE.(NEL-1)) GOTO 500
        ENDIF 
      ENDDO 
500   CONTINUE
      IF(ND1.EQ.ND2) STOP ' !!! PROBLEM IN SCR1 !!! '
      CALL GTID(NBASISMAX,ND1,ID1)
      CALL GTID(NBASISMAX,ND2,ID2)
C      CD=COULDAMP(ND1,NHG)*COULDAMP(ND2,NHG)
C..
      CALL CHK(G1(ND1)%Ms,G1(ND2)%Ms,IDS)
C..
      SUM2=0.D0
      DO I=1,NEL-1
        SUM3=DFLOAT(IDS)*
     &  GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID(I),ID2,ID(I))
     &           -DFLOAT(IDS1(I)*IDS2(I))*
     &  GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID1,ID(I),ID(I),ID2)
        SUM2=SUM2+SUM3
C*(COULDAMP(NDET1(I),NHG)**2)*CD
      ENDDO
C      IF(SUM2.NE.0.D0) CALL CHK_MS(NDET1,NDET2,NBASISMAX,NEL,2,G1)
C..  We then add in the non-diagonal part of the kinetic energy -
C..  <ph_a|T|ph_a'> where a and a' are the only basis fns that differ in
C..  NDET1 and NDET2
      SUM=SUM2*COULCOUPLE+TMAT(ND1,ND2)
!      CALL TIHALT('      SCR1',ISUB) 
      RETURN
      END
