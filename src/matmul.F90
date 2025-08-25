#include "macros.h"
module matmul_mod
    use constants, only: dp
    use global_utilities, only: timer, set_timer, halt_timer
    use blas_interface_mod, only: dgemm, dscal
    better_implicit_none
    private
    public :: my_hpsi

contains
    SUBROUTINE MY_HPSI(NDET, NEVAL, NROW, LAB, HAMIL, CK, CKN, TLargest)
        integer :: NEVAL, NDET, I, J, K, L, IBEG
        real(dp) :: HAMIL(*), CKN(NDET, NEVAL)
        INTEGER LAB(*), NROW(NDET)
        real(dp) :: CK(NDET, NEVAL), AUX
        LOGICAL :: TLargest
        type(timer), save :: proc_timer
! ==-----------------------------------------------------------------==
        proc_timer%timer_name = ' MY_HPSI  '
        call set_timer(proc_timer)
! ==-----------------------------------------------------------------==
!..Run over rows
        CKN = 0.0_dp
        DO I = 1, NDET
!..Run over columns
            IF (I == 1) THEN
                IBEG = 0
            ELSE
                IBEG = IBEG + NROW(I - 1)
            END IF
            DO K = 1, NEVAL
!..
                J = LAB(IBEG + 1)
                CKN(I, K) = CKN(I, K) + HAMIL(IBEG + 1) * CK(J, K)
                DO L = 2, NROW(I)
                    J = LAB(IBEG + L)
                    AUX = HAMIL(IBEG + L)
                    CKN(I, K) = CKN(I, K) + AUX * CK(J, K)
                    CKN(J, K) = CKN(J, K) + AUX * CK(I, K)
                END DO
!..
            END DO
        END DO
!..Need to do the following mult. by -1 so that eigenvalues are calc. in
!..ascending order
        IF (.NOT. TLargest) THEN
            CALL DSCAL(NEVAL * NDET, -1.0_dp, CKN, 1)
        END IF
! ==-----------------------------------------------------------------==
        call halt_timer(proc_timer)
! ==-----------------------------------------------------------------==
    END
! ==-----------------------------------------------------------------==
end module
