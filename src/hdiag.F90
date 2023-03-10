#include "macros.h"

module hdiag_mod
    use constants, only: dp, int32, stdout
    use SystemData, only: BasisFN
    use global_utilities, only: timer, set_timer, halt_timer
    use HElem, only: HElement_t_size
    use util_mod, only: stop_all
    use blas_interface_mod, only: dsyev, zheev
    better_implicit_none
    private
    public :: HDIAG_neci

contains

!Diagonalize a compressed hamiltonian (in HAMIL), returning eigenvalues in W and eigenvectors in CK.
    SUBROUTINE HDIAG_neci(NDET, HAMIL, LAB, NROW, CK, W, WORK2, WORK, NBLOCKSTARTS, NBLOCKS)
        INTEGER, intent(in) :: NDET, NROW(NDET), LAB(*)
        HElement_t(dp), intent(in) :: HAMIL(*)
        HElement_t(dp), intent(out) :: CK(NDET, NDET), WORK(*), WORK2(3 * NDET)
        real(dp), intent(inout) :: W(NDET)
        INTEGER, intent(in) :: NBLOCKS, NBLOCKSTARTS(NBLOCKS + 1)
        integer :: i, j, ind, indz, nbs
        INTEGER(int32) :: INFO
        type(timer), save :: proc_timer
        real(dp) :: GSEN
        character(len=*), parameter :: t_r = "HDIAG_neci"
#ifndef CMPLX_
        associate(tmp => WORK(1:1)); end associate
#endif

        proc_timer%timer_name = 'HDIAG     '
        call set_timer(proc_timer)
        GSEN = 1.D100
        CK(:, :) = 0._dp
!.. Now we fill the RIJ array
        IND = 1
        INDZ = 1
        DO I = 1, NDET
            INDZ = INDZ + NROW(I)
            DO WHILE (IND < INDZ)
                CK(I, LAB(IND)) = HAMIL(IND)
                IND = IND + 1
            END DO
        END DO
        IF (HElement_t_size /= 1) THEN
!ensure hermitian
            do i = 1, ndet
                do j = 1, i - 1
#ifdef CMPLX_
                    !To avoid compiler warnings
                    CK(i, j) = CONJG(CMPLX(CK(j, i), kind=dp))
#endif
                end do
            end do
        END IF

!****************************
!.. Diagonalize
        DO I = 1, NBLOCKS
            NBS = NBLOCKSTARTS(I)
            IF (HElement_t_size == 1) THEN
#ifndef CMPLX_
                CALL DSYEV('V', 'U', NBLOCKSTARTS(I + 1) - NBS, CK(NBS, NBS), NDET, W(NBS), WORK2, 3 * NDET, INFO)
#else
                call stop_all(t_r, "This block for real calcs only")
#endif
            ELSE
#ifdef CMPLX_
                CALL ZHEEV('V', 'U', NBLOCKSTARTS(I + 1) - NBS, CK(NBS, NBS), NDET, W(NBS), WORK, 4 * NDET, WORK2, INFO)
#else
                call stop_all(t_r, "This block for complex calcs only")
#endif
            END IF
            IF (INFO /= 0) THEN
                WRITE(6, *) 'DYSEV error: ', INFO
                call stop_all(t_r, "DSYEV error")
            END IF
            IF (W(NBS) < GSEN) GSEN = W(NBS)
        END DO

!.. CK now contains the eigenvectors, and W the eigenvalues
        WRITE(stdout, "(A,F19.11,I4)") "GROUND E=", GSEN
        call halt_timer(proc_timer)
    end subroutine hdiag_neci

end module hdiag_mod
