!==========================================================
!.. Routines needed for the initialisation of the Fourier
!.. Method of evaluating the Coulomb integrals
!  ========================================================
SUBROUTINE INITFOU(NMESHX, CK, NMAX, A, TALPHA, ALPHA, OMEGA, ZIA)
    use constants, only: dp, int64, sp
    use util_mod, only: neci_etime
    IMPLICIT NONE
    integer :: NMESHX, NMAX
    complex(dp) CK(NMESHX, NMESHX, NMESHX)
    complex(dp) ZIA(-NMESHX / 2:NMESHX / 2, NMAX, NMAX)
    real(dp) :: ALPHA, OMEGA
    real(dp) A(3), t3
    real(dp) t(2), t1, t2
    LOGICAL TALPHA
    INTEGER, SAVE :: IFIRST = 0
    IF (IFIRST == 1) RETURN
    IFIRST = 1
    T1 = neci_etime(t)
    CALL GEN_CK_FFT(NMESHX, CK, A, TALPHA, ALPHA, OMEGA)
    CALL GEN_ZIA(NMESHX, NMAX, ZIA)
    T2 = neci_etime(t)
    T3 = (T2 - T1)
    WRITE(6, *) 'V0=', CK(NMESHX / 2, NMESHX / 2, NMESHX / 2)
    WRITE(6, *) ' TIME FOR INITIALISATION:', T3 / 1000.
END
! =========================================================
SUBROUTINE GEN_CK_FFT(N, DIST, A, TALPHA, ALPHA, OMEGA)
    use constants, only: dp, int64, pi
    use util_mod, only: near_zero, stop_all

    IMPLICIT NONE
#if !defined(__alpha) || !defined(__SGI)
    INTEGER, parameter :: FFTW_FORWARD = -1, FFTW_ESTIMATE = 64
#endif
    integer :: N, I, J, K, N1, N2, N3
    complex(dp) DIST(-N / 2:N / 2 - 1, -N / 2:N / 2 - 1, -N / 2:N / 2 - 1)
    integer(int64) PLAN
    real(dp) A(3), ALPHA2, ALPHA, HSTEPX, HSTEPY, HSTEPZ, DERF
    real(dp) :: X, Y, Z, AUX, SUM, G2, GX, GY, GZ, OMEGA, OMEGAP, R
    LOGICAL TALPHA
    character(len=*), parameter :: t_r = "GEN_CK_FFT"
    ALPHA2 = ALPHA * ALPHA
    HSTEPX = 2.0_dp * A(1) / N
    HSTEPY = 2.0_dp * A(2) / N
    HSTEPZ = 2.0_dp * A(3) / N
!..
#if defined(NAGF95) || defined(GFORTRAN_) || defined(BLUEGENE_HACKS)
    call stop_all(t_r, "No ERF in NAG/GFortran?")
#endif
    DO I = -N / 2, N / 2 - 1
    DO J = -N / 2, N / 2 - 1
    DO K = -N / 2, N / 2 - 1
        X = real(I, dp) * HSTEPX
        Y = real(J, dp) * HSTEPY
        Z = real(K, dp) * HSTEPZ
        AUX = X * X + Y * Y + Z * Z
        IF (.not. near_zero(AUX)) THEN
        IF (TALPHA) THEN
            R = SQRT(AUX)
#if !defined(NAGF95) && !defined(GFORTRAN_) && !defined(BLUEGENE_HACKS)
            AUX = 1.0_dp / R * (-1)**(I + J + K) * DERF(R / ALPHA)
#endif
        ELSE
            AUX = 1.0_dp / SQRT(AUX) * (-1)**(I + J + K)
        END IF
        ELSE
        IF (TALPHA) THEN
            AUX = 2.0_dp / SQRT(PI) / ALPHA
        ELSE
            AUX = real(N / 2, dp) * (-1)**(I + J + K)
        END IF
        END IF
        DIST(I, J, K) = CMPLX(AUX, 0.0_dp, dp)
    END DO
    END DO
    END DO
#ifdef __SGI
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
    CALL DFFTW_PLAN_DFT_3D(PLAN, N1, N2, N3, DIST, DIST, FFTW_FORWARD, FFTW_ESTIMATE)
    CALL DFFTW_EXECUTE(PLAN)
    CALL DFFTW_DESTROY_PLAN(PLAN)
#else
    call stop_all("gen_ck_fft", "FFTW disabled")
#endif
#endif
!..Shift origin and normalise
    DO I = -N / 2, N / 2 - 1
    DO J = -N / 2, N / 2 - 1
    DO K = -N / 2, N / 2 - 1
        DIST(I, J, K) = DIST(I, J, K) * (-1)**(i + j + k) / real(n1 * n2 * n3, dp)
    end do
    end do
    end do
!..the short-range correction for alpha
    IF (TALPHA) THEN
        OMEGAP = 8.0_dp * OMEGA
        DO I = -N / 2, N / 2 - 1
        DO J = -N / 2, N / 2 - 1
        DO K = -N / 2, N / 2 - 1
            GX = PI * I / A(1)
            GY = PI * J / A(2)
            GZ = PI * K / A(3)
            G2 = GX * GX + GY * GY + GZ * GZ
            IF (.not. near_zero(G2)) THEN
                SUM = 4.0_dp * PI * (1.0_dp - AUX) / OMEGAP / G2
            ELSE
                SUM = PI * ALPHA2 / OMEGAP
            END IF
            DIST(I, J, K) = DIST(I, J, K) + CMPLX(SUM, 0.0_dp, dp)
        END DO
        END DO
        END DO
    END IF
END
! ============================================================
SUBROUTINE GEN_ZIA(KMAX, NMAX, ZIA)
    use constants, only: dp
    IMPLICIT NONE
    integer :: K, N, M, KMAX, NMAX
    complex(dp) ZIA(-KMAX / 2:KMAX / 2, NMAX, NMAX), ZIO
    DO K = -KMAX / 2, KMAX / 2
    DO N = 1, NMAX
    DO M = 1, NMAX
        ZIA(K, N, M) = 0.25_dp * (ZIO(N - M + K) + ZIO(-N + M + K) - ZIO(N + M + K) - ZIO(-N - M + K))
    END DO
    END DO
    END DO
    RETURN
END
FUNCTION ZIO(K)
    use constants, only: dp, pi
    IMPLICIT NONE
    complex(dp) ZIO
    INTEGER :: K
    IF (K == 0) THEN
        ZIO = (1.0_dp, 0.0_dp)
    ELSE
        ZIO = CMPLX(0.0_dp, -((-1.0_dp)**K - 1.0_dp) / (REAL(K, dp) * PI), dp)
    END IF
    RETURN
END
