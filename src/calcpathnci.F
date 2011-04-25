
C.. Generate a real that is (P H) (the binomial coefficient)
C.. Best for small H
      real(dp) FUNCTION FNCR(I_P,I_H)
         use constants, only: dp
         implicit none
         integer :: i,i_h,i_p
         FNCR=1.D0
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
         IF((ABS((A/B)-1.D0).LT.1.D-10).OR.
     &      ABS(A-B).LT.1.D-10) LSAME2=.TRUE.
C         IF(LSAME2) WRITE(10,*) "LSAME2",A,B,LSAME2
         RETURN
      END

C.. GETRHOELEMENT
C.. Get matrix element of rho
      real(dp) FUNCTION GETRHOELEMENT(II,IJ,RHOIJ,RHOIJLAB,             &
     &              RHOIJROWS,NDET)
      use constants, only: dp
      implicit none
      real(dp) RHOIJ(*)
      INTEGER RHOIJLAB(*)
      INTEGER RHOIJROWS(NDET),IJ,II,I,J,INDXROW,IMAX,K,NDET
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
      GETRHOELEMENT=0.D0
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
      GETRHOELEMENT=0.D0
      RETURN
      END
