C ========================================================
      SUBROUTINE SLATCOULFOU(G1,G2,G1P,G2P,N,
     &     CK,NMAX,ZIA,OUT)
C=========================================================
C.. Returns the Couloumb integral between the Slater 
C.. determinants of Sin function basis, using the Fourier 
C.. method. 
C=========================================================
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER :: N
      INTEGER G1(3),G2(3),G1P(3),G2P(3),NMAX
      complex(dp) ZIA(-N/2:N/2,NMAX,NMAX) 
      complex(dp) CK(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1)
      real(dp) :: SUM1,OUT

      SUM1=0.0_dp
      CALL VCOULFOU(N,G1,G2,G1P,G2P,SUM1,CK,NMAX,ZIA)
      OUT=SUM1*(2**6)
      RETURN
      END
C ===================================================================
      SUBROUTINE VCOULFOU(N,G1,G2,G1P,G2P,SUM,CK,NMAX,ZIA)
C =============================================================
C Returns the Coulomb integral between states (g1g2) & (g1pg2p)
C using the Fourier method 
C =============================================================
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER N,NMAX
      complex(dp) ZIA(-N/2:N/2,NMAX,NMAX) 
      INTEGER G1(3),G2(3),G1P(3),G2P(3),I,K,J
      complex(dp) CK(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1)
      complex(dp) CAUX,CAUX1,CAUX2,CAUX3,ZZERO
      real(dp) :: SUM
C..Sum over k
      SUM=0.0_dp
      ZZERO=(0.0_dp,0.0_dp)
      DO K=-N/2,N/2-1
        CAUX1=ZIA(K,G1(3),G1P(3))*CONJG(ZIA(K,G2(3),G2P(3)))
        IF(CAUX1.EQ.ZZERO) GOTO 100
        DO J=-N/2,N/2-1
         CAUX2=CAUX1*ZIA(J,G1(2),G1P(2))*CONJG(ZIA(J,G2(2),G2P(2)))
          IF(CAUX2.EQ.ZZERO) GOTO 90
          DO I=-N/2,N/2-1
            CAUX3=ZIA(I,G1(1),G1P(1))*CONJG(ZIA(I,G2(1),G2P(1)))
            CAUX=CK(I,J,K)*CAUX2*CAUX3
            SUM=SUM+REAL(CAUX,dp)
          ENDDO
 90       CONTINUE
        ENDDO
 100    CONTINUE
      ENDDO
C..
      RETURN 
      END
