#include "macros.h"
module fcoul_mod
    use constants, only: dp, EPS
    better_implicit_none
    private
    public :: SlatCoulFou

contains

! ========================================================
    SUBROUTINE SLATCOULFOU(G1, G2, G1P, G2P, N, CK, NMAX, ZIA, OUT)
!=========================================================
!.. Returns the Couloumb integral between the Slater
!.. determinants of Sin function basis, using the Fourier
!.. method.
!=========================================================
        INTEGER :: N
        INTEGER G1(3), G2(3), G1P(3), G2P(3), NMAX
        complex(dp) ZIA(-N / 2:N / 2, NMAX, NMAX)
        complex(dp) CK(-N / 2:N / 2 - 1, -N / 2:N / 2 - 1, -N / 2:N / 2 - 1)
        real(dp) :: SUM1, OUT

        SUM1 = 0.0_dp
        CALL VCOULFOU(N, G1, G2, G1P, G2P, SUM1, CK, NMAX, ZIA)
        OUT = SUM1 * (2**6)
    END
! ===================================================================
    SUBROUTINE VCOULFOU(N, G1, G2, G1P, G2P, SUM, CK, NMAX, ZIA)
! =============================================================
! Returns the Coulomb integral between states (g1g2) & (g1pg2p)
! using the Fourier method
! =============================================================
        INTEGER N, NMAX
        complex(dp) ZIA(-N / 2:N / 2, NMAX, NMAX)
        INTEGER G1(3), G2(3), G1P(3), G2P(3), I, K, J
        complex(dp) CK(-N / 2:N / 2 - 1, -N / 2:N / 2 - 1, -N / 2:N / 2 - 1)
        complex(dp) CAUX, CAUX1, CAUX2, CAUX3, ZZERO
        real(dp) :: SUM
!..Sum over k
        SUM = 0.0_dp
        ZZERO = (0.0_dp, 0.0_dp)
        DO K = -N / 2, N / 2 - 1
            CAUX1 = ZIA(K, G1(3), G1P(3)) * CONJG(ZIA(K, G2(3), G2P(3)))
            if (abs(CAUX1) <= EPS) goto 100
            DO J = -N / 2, N / 2 - 1
                CAUX2 = CAUX1 * ZIA(J, G1(2), G1P(2)) * CONJG(ZIA(J, G2(2), G2P(2)))
                if (abs(CAUX2) <= EPS) goto 90
                DO I = -N / 2, N / 2 - 1
                    CAUX3 = ZIA(I, G1(1), G1P(1)) * CONJG(ZIA(I, G2(1), G2P(1)))
                    CAUX = CK(I, J, K) * CAUX2 * CAUX3
                    SUM = SUM + REAL(CAUX, dp)
                END DO
90              CONTINUE
            END DO
100         CONTINUE
        END DO
    END
end module fcoul_mod
