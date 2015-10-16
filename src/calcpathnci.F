! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 

C.. Generate a real that is (P H) (the binomial coefficient)
C.. Best for small H
      FUNCTION FNCR(I_P,I_H)
         use constants, only: dp
         implicit none
         integer :: i,i_h,i_p
         real(dp) :: FNCR
         FNCR=1.0_dp
         DO I=0,I_H-1
            FNCR=FNCR*(I_P-I)/(I+1)
         ENDDO
         RETURN
      END

      LOGICAL FUNCTION LSAME2(A,B)
        use constants, only: dp
        IMPLICIT NONE
        REAL(dp) :: A,B
!         IMPLICIT real(dp)(A-H,O-Z)
!         IMPLICIT LOGICAL(L)
         LSAME2=.FALSE.
         IF((ABS((A/B)-1.0_dp).LT.1.0e-10_dp).OR.
     &      ABS(A-B).LT.1.0e-10_dp) LSAME2=.TRUE.
C         IF(LSAME2) WRITE(10,*) "LSAME2",A,B,LSAME2
         RETURN
      END

C.. GETRHOELEMENT
C.. Get matrix element of rho
      FUNCTION GETRHOELEMENT(II,IJ,RHOIJ,RHOIJLAB,
     &              RHOIJROWS,NDET)
      use constants, only: dp
      implicit none
      real(dp) RHOIJ(*),GETRHOELEMENT
      INTEGER RHOIJLAB(*),NDET
      INTEGER RHOIJROWS(NDET),IJ,II,I,J,INDXROW,IMAX,K
C.. We only have half of RHO, so if J<I, return the symmetrical (J,I) element
C.. Or if we have the whole RHO, it's quicker to look closer to its beginning
      IF(IJ<II) THEN
         I=IJ
         J=II
      ELSE
         I=II
         J=IJ
      ENDIF
C      IF(I.GT.0) THEN
      GETRHOELEMENT=0.0_dp
      INDXROW=RHOIJROWS(I)
      IMAX=RHOIJROWS(I+1)-1
      DO K=INDXROW,IMAX
         IF(RHOIJLAB(K).GT.J) RETURN
         IF(RHOIJLAB(K).EQ.J) THEN
            GETRHOELEMENT=RHOIJ(K)
            RETURN
         ENDIF
      ENDDO
C      ELSE
C         JK=1
C      ENDIF
      GETRHOELEMENT=0.0_dp
      RETURN
      END
