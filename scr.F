C.. TMAT has sometimes been hidden in ZIA
      SUBROUTINE SCR0(NEL,NBASISMAX,NDET1,NDET2,
     &	G1,NHG,UMAT,ALAT,ISS,TMAT,COULCOUPLE,SUM) 
      IMPLICIT NONE
      INTEGER G1(5,*),ISS
      REAL*8 UMAT(*)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      INTEGER NDET1(NEL),NDET2(NEL)
      INTEGER ID(NEL)
      INTEGER NBASISMAX(5,2)
C..Cube array
      REAL*8 ALAT(3)
      REAL*8 KP
      LOGICAL TVAB
      REAL*8 TMAT(NHG,NHG)
      REAL*8 GETUMATEL
      INTEGER NEL,NHG
      REAL*8 SUM,COULCOUPLE
      REAL*8 SUM1,SUM2,PI,SUM3,SUM4
      INTEGER IDS,I,J
      INTEGER ISUB
      REAL*8 COULDAMP
      CALL TISET('      SCR0',ISUB) 
      PI=ACOS(-1.D0)
C..This is the Slater-Condon rule when the two determinants are the *SAME*
      SUM1=0.D0
      SUM2=0.D0
C..We must first do the one electron integrals
      DO I=1,NEL
C.. we extract the KEs from TMAT
         SUM1=SUM1+TMAT(NDET1(I),NDET1(I))
C        SUM1=SUM1+((ALAT(1)**2)*((G(1,NDET1(I))**2)/(ALAT(1)**2)+
C     &      (G(2,NDET1(I))**2)/(ALAT(2)**2)+
C     &      (G(3,NDET1(I))**2)/(ALAT(3)**2)))
        CALL GTID(NBASISMAX,NDET1(I),ID(I))
      ENDDO
C..Cube multiplier
C      CST=PI*PI/(2.D0*ALAT(1)*ALAT(1))
C.. the UEG has k=2pi n/L rather than pi n/L, so we need to multiply the
C.. KE up by 4
C      IF(NBASISMAX(1,1).LE.0) CST=CST*4
C      SUM1=CST*SUM1
C.. If we're in the UEG, we need to add in the part of the background
C.. contribution which doesn't cancel
C         IF(NBASISMAX(1,1).LT.1) THEN
C            SUM=SUM-NEL*UMAT(1,1,1,1)/2
C         ENDIF
C..Now for the two electron part
      DO I=1,NEL-1
C        SUM3=COULDAMP(NDET1(I),NHG)**2
        DO J=I+1,NEL
          CALL CHK(G1(4,NDET1(I)),G1(4,NDET2(J)),IDS)
          SUM4=
     & GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID(I),ID(J),ID(I),ID(J))
     &   -DFLOAT(IDS)*
     & GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID(I),ID(J),ID(J),ID(I))
          SUM2=SUM2+SUM4
C  *SUM3*(COULDAMP(NDET2(J),NHG)**2)
        ENDDO
      ENDDO
C..
      IF(SUM2.NE.0.D0) CALL CHK_MS(NDET1,NDET2,NBASISMAX,
     &   NEL,1,G1)
      SUM=SUM1+SUM2*COULCOUPLE
      CALL TIHALT('      SCR0',ISUB) 
      RETURN
      END
C=============================================================
      SUBROUTINE CHK(N,L,IDS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IDS=N*L
      IF(IDS.LT.0) IDS=0 
      RETURN
      END
