C.. TMAT has sometimes been hidden in ZIA
      SUBROUTINE SCR0(NEL,NBASISMAX,NDET1,NDET2,
     &	G1,NHG,UMAT,ALAT,ISS,TMAT,COULCOUPLE,SUM)
      USE HElement
      IMPLICIT NONE
      INCLUDE 'umatcache.inc'
      INCLUDE 'basis.inc'
      TYPE(BasisFn) G1(*)
      INTEGER ISS
      TYPE(HElement) UMAT(*),TMAT(NHG,NHG)

C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      INTEGER NDET1(NEL),NDET2(NEL)
      INTEGER ID(NEL)
      INTEGER NBASISMAX(5,2)
C..Cube array
      REAL*8 ALAT(3)
      REAL*8 KP
      LOGICAL TVAB
      TYPE(HElement) GETUMATEL
      INTEGER NEL,NHG
      TYPE(HElement) SUM,SUM1,SUM2,SUM3,SUM4
      REAL*8 COULCOUPLE
      INTEGER IDS,I,J,IG1I,INDET
      INTEGER ISUB
      REAL*8 COULDAMP
c      PARAMETER PI=3.14159265358979323846264338327950288419716939937510D0
!      CALL TISET('      SCR0',ISUB) 
C      PI=ACOS(-1.D0)
C..This is the Slater-Condon rule when the two determinants are the *SAME*
      SUM1=0.D0
      SUM2=0.D0
C..We must first do the one electron integrals
      DO I=1,NEL
C.. we extract the KEs from TMAT
         INDET=NDET1(I)
         SUM1=SUM1+TMAT(INDET,INDET)
        CALL GTID(NBASISMAX,NDET1(I),ID(I))
      ENDDO
C..Now for the two electron part
C  JSS --access the <ij|ij> and <ij|ji> integrals directly.  Note J>I and how the integrals are
C  stored in UMAT2D (see note in cpmddata.inc).   
      DO I=1,NEL-1
C        SUM3=COULDAMP(NDET1(I),NHG)**2
         IG1I=G1(NDET1(I))%Ms
        DO J=I+1,NEL
          CALL CHK(IG1I,G1(NDET2(J))%Ms,IDS)
          IF(TUMAT2D) THEN
             SUM4=UMAT2D(ID(I),ID(J))
          ELSE
            SUM4=
     & GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID(I),ID(J),ID(I),ID(J))
          ENDIF
          IF(IDS.NE.0) THEN 
            IF(TUMAT2D) THEN
               SUM4=SUM4-HElement(IDS)*UMAT2D(ID(J),ID(I))
            ELSE
               SUM4=SUM4-HElement(IDS)*
     &GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,ID(I),ID(J),ID(J),ID(I))
            ENDIF
          ENDIF
          SUM2=SUM2+SUM4
C  *SUM3*(COULDAMP(NDET2(J),NHG)**2)
        ENDDO
      ENDDO
C..
C      IF(SUM2.NE.0.D0) CALL CHK_MS(NDET1,NDET2,NBASISMAX,
C     &   NEL,1,G1)
      SUM=SUM1+SUM2*HElement(COULCOUPLE)
!      CALL TIHALT('      SCR0',ISUB) 
      RETURN
      END
C=============================================================
      SUBROUTINE CHK(N,L,IDS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IDS=N*L
      IF(IDS.LT.0) IDS=0 
      RETURN
      END
