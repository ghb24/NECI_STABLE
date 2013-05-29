! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
      SUBROUTINE RLSXC_SC(NEL,NMSH,NMAX,SITAB,IOBS,JOBS,KOBS,
     &    NDET1,G1,NHG,AUXC)
      use SystemData, only: BasisFN
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: NHG,NMAX,NMSH,NEL,KI,KJ,IX,IY,IZ,IDS
      integer :: NDET1(NEL)
      real(dp) :: AUXC(NMSH,NMSH,NMSH)
      real(dp) :: SITAB(NMSH,NMAX)
      TYPE(BasisFN) G1(NHG)
      real(dp) :: U1KIR1,U2KJR2,U1KJR1,U2KIR2,U1
      integer :: IOBS,JOBS,KOBS
C..
      AUXC=0.0_dp
      DO KI=1,NEL-1
        DO KJ=KI+1,NEL
C..
          CALL CHK(G1(NDET1(KI))%Ms,G1(NDET1(KJ))%Ms,IDS)
          DO IZ=1,NMSH
          DO IY=1,NMSH
          DO IX=1,NMSH
C..
C..U1KIR1=U1(NDET1(KI),NMAX,LMAX,G1,RAD,R1Z)
        U1KIR1=U1(NDET1(KI),NMAX,NMSH,NHG,G1,SITAB,
     &      IOBS,JOBS,KOBS)
C..U2KJR2=U2(NDET1(KJ),NMAX,LMAX,G1,RAD,R2,THT)
        U2KJR2=U1(NDET1(KJ),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
C..U1KJR1=U1(NDET1(KJ),NMAX,LMAX,G1,RAD,R1Z)
        U1KJR1=U1(NDET1(KJ),NMAX,NMSH,NHG,G1,SITAB,
     &      IOBS,JOBS,KOBS)
C..U2KIR2=U2(NDET1(KI),NMAX,LMAX,G1,RAD,R2,THT)
        U2KIR2=U1(NDET1(KI),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
C..
        AUXC(IX,IY,IZ)=AUXC(IX,IY,IZ)+
     &      U1KIR1*U1KIR1*U2KJR2*U2KJR2+
     &      U1KJR1*U1KJR1*U2KIR2*U2KIR2-
     &       U1KIR1*U1KJR1*U2KJR2*U2KIR2*real(IDS,dp)-
     &      U1KJR1*U1KIR1*U2KIR2*U2KJR2*real(IDS,dp)
          ENDDO
            ENDDO
          ENDDO
C..
        ENDDO
      ENDDO
      RETURN
      END
C ==-----------------------------------------------------------------==
      FUNCTION U1(ND,NX,NMSH,NHG,G1,SITAB,I,J,K)
      use SystemData, only: BasisFN
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: NMSH,NX,NHG,I,J,K,ND
      real(dp) :: SITAB(NMSH,NX),U1
      TYPE(BasisFN) G1(NHG)
C..
      U1=
     &   SITAB(I,G1(ND)%k(1))*SITAB(J,G1(ND)%k(2))*SITAB(K,G1(ND)%k(3))
      RETURN
      END
C =====================================
C We want to get rid of this really ... it is a bit opaque
      SUBROUTINE CHK(N,L,IDS)
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: N,L,IDS,DS
      DS=N*L
      IF(IDS.LT.0) IDS=0
      RETURN
      END
