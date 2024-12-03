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

        do k = 1, neval
!..
            j = lab(1)
            ckn(1, k) = ckn(1, k) + hamil(1) * ck(j, k)
            do l = 2, nrow(1)
                j = lab(l)
                aux = hamil(l)
                ckn(1, k) = ckn(1, k) + aux * ck(j, k)
                ckn(j, k) = ckn(j, k) + aux * ck(1, k)
            end do
!..
        end do

        ibeg = 0
        do i = 2, ndet
!..run over columns
            ibeg = ibeg + nrow(i - 1)
            do k = 1, neval
!..
                j = lab(ibeg + 1)
                ckn(i, k) = ckn(i, k) + hamil(ibeg + 1) * ck(j, k)
                do l = 2, nrow(i)
                    j = lab(ibeg + l)
                    aux = hamil(ibeg + l)
                    ckn(i, k) = ckn(i, k) + aux * ck(j, k)
                    ckn(j, k) = ckn(j, k) + aux * ck(i, k)
                end do
!..
            end do
        end do
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
