#include "macros.h"

module read_psi_mod
    use constants, only: dp, EPS
    use util_mod, only: stop_all
    better_implicit_none
    private
    public :: read_psi, write_psi_comp

contains

    SUBROUTINE READ_PSI(BOX, BOA, COA, NDET, NEVAL, NBASISMAX, NEL, CK, W)
        INTEGER nBasisMax(5, *)
        integer :: NEVAL, NDET, I, J, NEL, NDETCK, NEVALCK, NELCK
        real(dp) :: CK(NDET, NEVAL)
        real(dp) :: W(NEVAL)
        integer :: NMAXCK, NMAYCK, NMAZCK
        real(dp) :: BOX, BOA, COA, C, BOXCK, BOACK, COACK
        character(*), parameter :: this_routine = 'read_psi'
        OPEN(10, FILE='PSI_LONG', STATUS='OLD', ERR=10)
        READ(10, *) NMAXCK, NMAYCK, NMAZCK, NELCK, NEVALCK, NDETCK
        READ(10, *) BOXCK, BOACK, COACK
        IF (NMAXCK /= NBASISMAX(1, 2) .OR. NMAYCK /= NBASISMAX(2, 2) .OR. NMAZCK /= NBASISMAX(3, 2)) then
            call stop_all(this_routine, ' !!! DIFFERENT SIZE BASIS !!! ')
        end if
        IF (NELCK /= NEL) call stop_all(this_routine, ' !!! DIFFERENT NUMBER OF ELECTRONS !!! ')
        IF (NEVALCK /= NEVAL) call stop_all(this_routine, ' !!! NEVAL DIFFERENT !!! ')
        IF (NDETCK /= NDET) call stop_all(this_routine, ' !!! DIFFERENT NUMBER OF DETS !!! ')
        IF (abs(BOXCK - BOX) > EPS .or. abs(BOACK - BOA) > EPS .or. abs(COACK - COA) > EPS) then
            call stop_all(this_routine, ' !!! DIFFERENT SIZE CUBES !!! ')
        end if
        DO I = 1, NEVAL
            DO J = 1, NDET
                READ(10, *) CK(J, I)
            END DO
        END DO
        DO I = 1, NEVAL
            READ(10, *) W(I)
        END DO
        CLOSE(10)
        WRITE(6, '(A)') ' HAVE READ IN PSI FROM FILE ./PSI_LONG '
        RETURN
10      CONTINUE
!.. PSI_LONG doesn't exist.  what about PSI_COMP
        OPEN(10, FILE='PSI_COMP', STATUS='OLD')
        READ(10, *) NMAXCK, NMAYCK, NMAZCK, NELCK, NEVALCK, NDETCK
        READ(10, *) BOXCK, BOACK, COACK
        IF (NMAXCK /= NBASISMAX(1, 2) .OR. NMAYCK /= NBASISMAX(2, 2) .OR. NMAZCK /= NBASISMAX(3, 2)) then
            call stop_all(this_routine, ' !!! DIFFERENT SIZE BASIS !!! ')
        end if
        IF (NELCK /= NEL) call stop_all(this_routine, ' !!! DIFFERENT NUMBER OF ELECTRONS !!! ')
        IF (NEVALCK /= NEVAL) call stop_all(this_routine, ' !!! NEVAL DIFFERENT !!! ')
        IF (NDETCK /= NDET) call stop_all(this_routine, ' !!! DIFFERENT NUMBER OF DETS !!! ')
        IF (abs(BOXCK - BOX) > EPS .or. abs(BOACK - BOA) > EPS .or. abs(COACK - COA) > EPS) then
            call stop_all(this_routine, ' !!! DIFFERENT SIZE CUBES !!! ')
        end if
        I = 1
        DO WHILE (I /= 0)
            READ(10, *) I, J, C
            IF (I /= 0) CK(J, I) = C
        END DO
        DO I = 1, NEVAL
            READ(10, *) W(I)
        END DO
        CLOSE(10)
    end subroutine


    SUBROUTINE WRITE_PSI(BOX, BOA, COA, NDET, NEVAL, NBASISMAX, NEL, CK, W)
        INTEGER nBasisMax(5, *)
        INTEGER I, NEVAL, NDET, NEL, J
        HElement_t(dp) CK(NDET, NEVAL)
        real(dp) W(NEVAL), BOA, BOX, COA
        OPEN(10, FILE='PSI_LONG', STATUS='UNKNOWN')
        WRITE(10, *) NBASISMAX(1, 2), NBASISMAX(2, 2), NBASISMAX(3, 2), NEL, NEVAL, NDET
        WRITE(10, *) BOX, BOA, COA
        DO I = 1, NEVAL
            DO J = 1, NDET
                WRITE(10, *) CK(J, I), ABS(CK(J, I))
            END DO
        END DO
        DO I = 1, NEVAL
            WRITE(10, *) W(I)
        END DO
        CLOSE(10)
    end subroutine

    SUBROUTINE WRITE_PSI_COMP(BOX, BOA, COA, NDET, NEVAL, NBASISMAX, NEL, CK, W)
        INTEGER I, NEVAL, NDET, NEL, J
        INTEGER nBasisMax(5, *)
        HElement_t(dp) CK(NDET, NEVAL)
        real(dp) W(NEVAL), BOA, BOX, COA

        OPEN(10, FILE='PSI_COMP', STATUS='UNKNOWN')
        WRITE(10, *) NBASISMAX(1, 2), NBASISMAX(2, 2), NBASISMAX(3, 2), NEL, NEVAL, NDET
        WRITE(10, *) BOX, BOA, COA
        DO I = 1, NEVAL
            DO J = 1, NDET
                IF (abs(CK(J, I)) > 1.0e-2_dp) THEN
                    WRITE(10, "(2I5)", advance='no') I, J
                    WRITE(10, *) CK(J, I), ABS(CK(J, I))
                END IF
            END DO
        END DO
        WRITE(10, *) 0, 0, 0
        DO I = 1, NEVAL
            WRITE(10, *) W(I)
        END DO
        CLOSE(10)
    end subroutine

end module
