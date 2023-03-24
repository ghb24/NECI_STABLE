#include "macros.h"
module init_coul2D_mod
    use constants, only: dp, sp, pi, int64
    use util_mod, only: stop_all, near_zero, neci_etime
    use init_coul_mod, only: gen_zia, gen_ck_fft
    better_implicit_none
    private
    public :: INITFOU2D

contains

!==========================================================
!.. Routines needed for the initialisation of the Fourier
!.. Method of evaluating the Coulomb integrals
!  ========================================================
    SUBROUTINE INITFOU2D(NMESHX, CK, NMAX, A, TALPHA, ALPHA, OMEGA, ZIA)
        integer :: NMESHX, NMAX
        complex(dp) CK(NMESHX, NMESHX, NMESHX)
        complex(dp) ZIA(-NMESHX / 2:NMESHX / 2, NMAX, NMAX)
        real(dp) A(3), OMEGA, ALPHA
        real(dp) t1, t2, t3, t(2)
        LOGICAL TALPHA
        INTEGER, SAVE :: IFIRST = 0
        IF (IFIRST == 1) RETURN
        IFIRST = 1
        T1 = neci_etime(t)
        CALL GEN_CK_FFT2D(NMESHX, CK, A, TALPHA, ALPHA, OMEGA)
        CALL GEN_ZIA(NMESHX, NMAX, ZIA)
        T2 = neci_etime(t)
        T3 = (T2 - T1)
        WRITE(stdout, *) 'V0=', CK(NMESHX / 2 + 1, NMESHX / 2 + 1, NMESHX / 2 + 1)
        WRITE(stdout, *) ' TIME FOR INITIALISATION:', T3 / 1000.0_dp
    end subroutine
! =========================================================
    SUBROUTINE GEN_CK_FFT2D(N, DIST, A, TALPHA, ALPHA, OMEGA)
        integer :: i, j, k, n, n1, n2, n3
        complex(dp) DIST(-N / 2:N / 2 - 1, -N / 2:N / 2 - 1, -N / 2:N / 2 - 1)
        real(dp) A(3), x, y, z, aux, r, ALPHA, ALPHA2, OMEGA, OMEGAP
        LOGICAL TALPHA
        complex(dp) D(-N / 2:N / 2 - 1, -N / 2:N / 2 - 1)
        real(dp) :: HSTEPX, HSTEPY, HSTEPZ, gx, gy, gz, SUM, G
        character(len=*), parameter :: t_r = "GEN_CK_FFT2D"
        ALPHA2 = ALPHA * ALPHA
        HSTEPX = 2.0_dp * A(1) / N
        HSTEPY = 2.0_dp * A(2) / N
        HSTEPZ = 0

#if defined(NAGF95) || defined(GFORTRAN_) || defined(BLUEGENE_HACKS)
        call stop_all(t_r, "No ERF in NAG?")
#endif

        DO I = -N / 2, N / 2 - 1
        DO J = -N / 2, N / 2 - 1
            K = 0
            X = real(I, dp) * HSTEPX
            Y = real(J, dp) * HSTEPY
            Z = real(K, dp) * HSTEPZ
            AUX = X * X + Y * Y + Z * Z
            IF (.not. near_zero(AUX)) THEN
            IF (TALPHA) THEN
                R = SQRT(AUX)
#if !defined(NAGF95) && !defined(GFORTRAN_) && !defined(BLUEGENE_HACKS)
                AUX = 1.0_dp / R * (-1)**(I + J) * DERF(R / ALPHA)
#endif
            ELSE
                AUX = 1.0_dp / SQRT(AUX) * (-1)**(I + J)
            END IF
            ELSE
            IF (TALPHA) THEN
                AUX = (2.0_dp / SQRT(PI) / ALPHA) * (-1)**(I + J)
            ELSE
                AUX = real(N / 2, dp) * (-1)**(I + J)
            END IF
            END IF
            D(I, J) = CMPLX(AUX, 0.0_dp, dp)

        END DO
        END DO
#ifdef __SGI
!..FFT parameters
        N1 = N
        N2 = N
        N3 = N
        LA1 = N1
        LA2 = N2
!..Initialise FFT
        CALL ZFFT3DI(N1, N2, N3, COEFF)
!..FORWARD TRANSFORM
        CALL ZFFT3D(-1, N1, N2, N3, DIST, LA1, LA2, COEFF)
#elif __alpha
        N1 = N
        N2 = N
        N3 = N
        LA1 = N1
        LA2 = N2
        CALL ZFFT_3D('C', 'C', 'F', DIST, DIST, N1, N2, N3, LA1, LA2, 1, 1, 1)
#else
        N1 = N
        N2 = N
        N3 = N
#ifndef DISABLE_FFTW
        CALL DFFTW_PLAN_DFT_2D(PLAN, N1, N2, D, DIST(-N / 2, -N / 2, 0), FFTW_FORWARD, FFTW_ESTIMATE)
        CALL DFFTW_EXECUTE(PLAN)
        CALL DFFTW_DESTROY_PLAN(PLAN)
#else
        call stop_all("gen_ck_fft2d", "FFTW disabled")
#endif
#endif
!..Shift origin and normalise
        DO I = -N / 2, N / 2 - 1
        DO J = -N / 2, N / 2 - 1
            K = 0
            DIST(I, J, 0) = DIST(I, J, 0) * (-1)**(i + j + k) / real(n1 * n2, dp)
        end do
        end do
!..the short-range correction for alpha
        IF (TALPHA) THEN

            OMEGAP = 4.0_dp * OMEGA
            DO I = -N / 2, N / 2 - 1
            DO J = -N / 2, N / 2 - 1
                K = 0
                GX = PI * I / A(1)
                GY = PI * J / A(2)
                GZ = 0
                G = SQRT(GX * GX + GY * GY + GZ * GZ)
                IF (.not. near_zero(G)) THEN
#if !defined(NAGF95) && !defined(GFORTRAN_) && !defined(BLUEGENE_HACKS)
                    SUM = 2.0_dp * PI * DERF(ALPHA * G / 2) / OMEGAP / G
#endif
                ELSE
                    SUM = ALPHA / (SQRT(PI) * OMEGAP)
                END IF
                DIST(I, J, K) = DIST(I, J, K) + CMPLX(SUM, 0.0_dp, dp)
            END DO
            END DO
        END IF
    end subroutine

end module
