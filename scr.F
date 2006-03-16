C.. TMAT has sometimes been hidden in ZIA
      SUBROUTINE SCR0(NEL,NBASISMAX,NDET1,NDET2,
     &	G,NHG,UMAT,ALAT,ISS,TMAT,COULCOUPLE,SUM) 
      IMPLICIT NONE
      INTEGER G(4,*),ISS
      REAL*8 UMAT(NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      INTEGER NDET1(NEL),NDET2(NEL)
      INTEGER ID(NEL)
      INTEGER NBASISMAX(4,2)
C..Cube array
      REAL*8 ALAT(3)
      REAL*8 KP
      LOGICAL TVAB
      REAL*8 TMAT(NHG,NHG)
      REAL*8 GETUMATEL
      INTEGER NEL,NHG
      REAL*8 SUM,COULCOUPLE
      REAL*8 SUM1,SUM2,PI
      INTEGER IDS,I,J
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
        DO J=I+1,NEL
          CALL CHK(G(4,NDET1(I)),G(4,NDET2(J)),IDS)
          SUM2=SUM2+
     & GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G,ID(I),ID(J),ID(I),ID(J))
     &   -DFLOAT(IDS)*
     & GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G,ID(I),ID(J),ID(J),ID(I))
        ENDDO
      ENDDO
C..
      IF(SUM2.NE.0.D0) CALL CHK_MS(NDET1,NDET2,NBASISMAX,
     &   NEL,1,G)
      SUM=SUM1+SUM2*COULCOUPLE
      RETURN
      END
C=============================================================
      SUBROUTINE CHK(N,L,IDS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IDS=N*L
      IF(IDS.LT.0) IDS=0 
      RETURN
      END

      REAL*8 FUNCTION GETUMATEL(NBASISMAX,UMAT,ALAT,NHG,ISS,G1,
     &      IDI,IDJ,IDK,IDL)
         IMPLICIT NONE
         INTEGER NBASISMAX(4,3),I,J,K,L,NHG,ISS
         INTEGER G1(4,NHG)
         REAL*8 ALAT(3)
         REAL*8 UMAT(NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         INTEGER A,B,C
         INTEGER IDI,IDJ,IDK,IDL
         REAL*8 SUM
         PARAMETER PI=3.14159265358979323846264338327950288419716939937510D0

C.. IF NBASISMAX(1,3) is less than zero, we directly give the integral.
C.. Otherwise we just look it up in umat
         IF(NBASISMAX(1,3).GE.0) THEN
            GETUMATEL=UMAT(IDI,IDJ,IDK,IDL)
            RETURN
         ENDIF
         IF(NBASISMAX(1,3).EQ.-1) THEN
            I=(IDI-1)*ISS+1
            J=(IDJ-1)*ISS+1
            K=(IDK-1)*ISS+1
            L=(IDL-1)*ISS+1
C.. The Uniform electron gas
                   A=G1(1,I)-G1(1,K)
                   B=G1(2,I)-G1(2,K)
                   C=G1(3,I)-G1(3,K)
                   IF(    (G1(1,L)-G1(1,J)).EQ.A
     &               .AND.(G1(2,L)-G1(2,J)).EQ.B
     &               .AND.(G1(3,L)-G1(3,J)).EQ.C
     &               .AND.((A.NE.0).OR.(B.NE.0).OR.(C.NE.0))) THEN
C                     WRITE(6,*) "(",I,J,"|",K,L,")",A,B,C

CC  AJWT  <IJ|r_12^-1|KL> = v_(G_I-G_K) delta_((G_I-G_K)-(G_L-G_J)
CC  v_G = 4 Pi/ G**2.  G=2 Pi/L(nx,ny,nx) etc.
                     SUM=((A/ALAT(1))**2+(B/ALAT(2))**2)
                     IF(ALAT(3).NE.0.D0) SUM=SUM+(C/ALAT(3))**2
                     SUM=1/(PI*SUM*ALAT(1)*ALAT(2)*ALAT(3))
                   ELSE
                     SUM=0.D0
                   ENDIF
            GETUMATEL=SUM
         ENDIF
            
      END
