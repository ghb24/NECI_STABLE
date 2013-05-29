! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
      SUBROUTINE RLSXC2_SC(NEL,NMSH,NMAX,SITAB,IOBS,JOBS,KOBS,
     &  NDET1,NDET2,G1,NHG,AUXC)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER IDS(4),NMSH,NMAX,NEL
      integer :: NDET1(NEL),NDET2(NEL)
      real(dp) :: AUXC(NMSH,NMSH,NMSH)
C..Cube array
      real(dp) :: SITAB(NMSH,NMAX)
      INTEGER G1(*),ND1(2),ND2(2)
      EXTERNAL U1
      INTEGER :: IQ,I,J,K,IOBS,JOBS,KOBS,NHG,IX,IY,IZ
      real(dp) :: U2ND22R2,U2ND21R2,U2ND12R2,U2ND11R2
      real(dp) :: U1ND22R1,U1ND21R1,U1ND12R1,U1ND11R1,U1

C..Zero ND1 and ND2
      AUXC=0.0_dp
      ND1(1:2)=0
      ND2(1:2)=0
C..IQ ensures that only two differences are recorded
C..First routine compares NDET1 to NDET2
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
C.. 
      IDS(1:4)=0
      CALL GENKRON2(G1,NHG,ND1,ND2,IDS)
      AUXC=0.0_dp
C..
C..electron in orbital ND1(1) at position R1Z, 0 
      U1ND11R1=U1(ND1(1),NMAX,NMSH,NHG,G1,SITAB,IOBS,JOBS,KOBS)
C..electron in orbital ND1(2) at position R1Z, 0
      U1ND12R1=U1(ND1(2),NMAX,NMSH,NHG,G1,SITAB,IOBS,JOBS,KOBS)
C..electron in orbital ND2(1) at position R1Z, 0
      U1ND21R1=U1(ND2(1),NMAX,NMSH,NHG,G1,SITAB,IOBS,JOBS,KOBS)
C..electron in orbital ND2(2) at position R1Z, 0
      U1ND22R1=U1(ND2(2),NMAX,NMSH,NHG,G1,SITAB,IOBS,JOBS,KOBS)
C..
      DO IZ=1,NMSH
        DO IY=1,NMSH
          DO IX=1,NMSH
C..electron in orbital ND1(1) at position R2, THT
            U2ND11R2=U1(ND1(1),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ) 
C..electron in orbital ND1(2) at position R2, THT
            U2ND12R2=U1(ND1(2),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
C..electron in orbital ND2(1) at position R2, THT
            U2ND21R2=U1(ND2(1),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
C..electron in orbital ND2(2) at position R2, THT
            U2ND22R2=U1(ND2(2),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
C..
            AUXC(IX,IY,IZ)=AUXC(IX,IY,IZ)+(
     &    U1ND12R1*U1ND22R1*U2ND11R2*U2ND21R2*real(
     &    IDS(1)*IDS(2),dp)+
     &          U1ND11R1*U1ND21R1*U2ND12R2*U2ND22R2*real(
     &          IDS(1)*IDS(2),dp)-
     &    U1ND12R1*U1ND21R1*U2ND11R2*U2ND22R2*real(
     &    IDS(3)*IDS(4),dp)-
     &    U1ND11R1*U1ND22R1*U2ND12R2*U2ND21R2*real(
     &    IDS(3)*IDS(4),dp))    
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C ==--------------------------------------------------------==
      SUBROUTINE GENKRON2(G,NHG,ND1,ND2,IDS)
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER NHG,G(4,NHG)
      INTEGER ND1(2),ND2(2),IDS(4)
C..
      CALL CHK(G(4,ND1(2)),G(4,ND2(2)),IDS(1))
      CALL CHK(G(4,ND1(1)),G(4,ND2(1)),IDS(2))
      CALL CHK(G(4,ND1(2)),G(4,ND2(1)),IDS(3))
      CALL CHK(G(4,ND1(1)),G(4,ND2(2)),IDS(4))
C..
      RETURN
      END
