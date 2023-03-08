#include "macros.h"

module lineup_mod
    use util_mod, only: NECI_ICOPY
    better_implicit_none
    private
    public :: lineup

contains

    SUBROUTINE LINEUP(N, LST1, LST2, NBASIS, L1E, L2E, IFLAG, IPTOT)
        INTEGER :: N, LST1(N), LST2(N), nBasis, I, iC, iFlag, iP1, iP2, iPTot, j, iP
        INTEGER :: LIST1(N), LIST2(N), L1E(NBASIS), L2E(NBASIS), LD(NBASIS)

        CALL NECI_ICOPY(N, LST1, 1, LIST1, 1)
        CALL NECI_ICOPY(N, LST2, 1, LIST2, 1)
        CALL SCATTER_neci(N, NBASIS, LIST1, L1E)
        CALL SCATTER_neci(N, NBASIS, LIST2, L2E)
        IC = 0
        DO I = 1, NBASIS
            LD(I) = L1E(I) - L2E(I)
            IF (LD(I) > 0) IC = IC + 1
        END DO
        IF (IC > 2) THEN
            IFLAG = 3
            RETURN
        END IF
        IF (IC == 0) THEN
            CALL ASCENDING_ORDER(N, LIST1, IP1)
            CALL ASCENDING_ORDER(N, LIST2, IP2)
            IPTOT = IP1 * IP2
            IFLAG = 0
            RETURN
        END IF
        IF (IC == 1) THEN
            IPTOT = 1
!..find out which one in list one differ
            DO I = 1, NBASIS
            IF (LD(I) == 1) THEN
            DO J = 1, N
            IF (LIST1(J) == I) THEN
                CALL CYCLIC_REORDER(N - J + 1, LIST1(J), IP)
                IPTOT = IPTOT * IP
            END IF
            END DO
            END IF
            CONTINUE
            IF (LD(I) == -1) THEN
            DO J = 1, N
            IF (LIST2(J) == I) THEN
                CALL CYCLIC_REORDER(N - J + 1, LIST2(J), IP)
                IPTOT = IPTOT * IP
            END IF
            END DO
            END IF
            END DO
            CONTINUE
            CALL ASCENDING_ORDER(N - 1, LIST1, IP1)
            CALL ASCENDING_ORDER(N - 1, LIST2, IP2)
            IPTOT = IPTOT * IP1 * IP2
            IFLAG = 1
            RETURN
        END IF
        IF (IC == 2) THEN
            IPTOT = 1
            IFLAG = 2
!..find out which two in list one differ
            DO I = 1, NBASIS
            IF (LD(I) == 1) THEN
            DO J = 1, N
            IF (LIST1(J) == I) THEN
                CALL CYCLIC_REORDER(N - J + 1, LIST1(J), IP)
                IPTOT = IPTOT * IP
            END IF
            END DO
            END IF
            IF (LD(I) == -1) THEN
            DO J = 1, N
            IF (LIST2(J) == I) THEN
                CALL CYCLIC_REORDER(N - J + 1, LIST2(J), IP)
                IPTOT = IPTOT * IP
            END IF
            END DO
            END IF
            END DO
            CALL ASCENDING_ORDER(N - 2, LIST1, IP1)
            CALL ASCENDING_ORDER(N - 2, LIST2, IP2)
            IPTOT = IPTOT * IP1 * IP2
            CALL ASCENDING_ORDER(2, LIST1(N - 1), IP1)
            CALL ASCENDING_ORDER(2, LIST2(N - 1), IP2)
            IPTOT = IPTOT * IP1 * IP2
        END IF
    END

    SUBROUTINE SCATTER_neci(N, NBASIS, LIST, LE)
        integer :: n, nBasis, LIST(N), LE(NBASIS)

        integer :: i
        LE(1:NBASIS) = 0
        DO I = 1, N
            LE(LIST(I)) = 1
        END DO
        RETURN
    END

    SUBROUTINE CYCLIC_REORDER(N, L, IP)
        integer :: N, L(N), iP
        integer :: iS, I
        IS = L(1)
        DO I = 2, N
            L(I - 1) = L(I)
        END DO
        L(N) = IS
        IP = (-1)**(N - 1)
        RETURN
    END

    SUBROUTINE ASCENDING_ORDER(N, L, IP)
        integer :: N, l(n), iP
        integer :: i, iT
        IP = 1
100     CONTINUE
        DO I = 1, N - 1
        IF (L(I) > L(I + 1)) THEN
            IT = L(I)
            L(I) = L(I + 1)
            L(I + 1) = IT
            IP = IP * (-1)
            GOTO 100
        END IF
        END DO
        RETURN
    END

end module
