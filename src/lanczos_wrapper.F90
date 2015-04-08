#include  "macros.h"

module lanczos_wrapper

    use constants

    implicit none

contains

    subroutine frsblk_wrapper(det_list, ndets, nexcit, evals, evecs)

        ! A wrapper to the frsblk Lanczos routine. This wrapper has a much
        ! smaller number of input parameters than frsblkh. The input parameters
        ! are:

        ! det_list: A list of all determinants in the space to be diagonalised,
        !     in orbital form, i.e. det_list(:,i) contains the spin orbitals
        !     occupied in the i'th determinant.
        ! ndets: The number of determinants in the space to be diagonalised,
        !     which should be equal to the size(det_list,2).
        ! nexcit: The number of eigenvectors to calculate. The nexcit lowest
        !     states will be calculated.

        ! On input, evals and evecs should be allocated already.

        ! On output, evals stores the nexict lowest eigenvalues and evecs
        ! stores the corresponding eigenvectors.

        use DetCalcData, only: nkry, nblk, b2l, ncycle
        use SystemData, only: nel, tHPHF

        integer, intent(in) :: det_list(:,:)
        integer, intent(in) :: ndets
        integer, intent(in) :: nexcit 
        real(dp), intent(out) :: evals(:)
        real(dp), intent(out) :: evecs(:,:)

        integer :: ierr, nkry1, nblock, len_scr, len_iscr, ICMax, LenHamil
        integer, allocatable :: nRow(:), Lab(:), iscr(:), ind(:)
        real(dp), allocatable :: evecs_space(:,:), Hamil(:), A_Arr(:,:)
        real(dp), allocatable :: V(:), BM(:), T(:), WT(:), scr(:), AM(:)
        logical :: tMC
        character(len=*), parameter :: t_r = 'frsblk_wrapper'

        allocate(nRow(ndets), stat=ierr)
        nRow = 0
        ICMax = 1
        tMC = .false.

        call Detham(ndets, nel, det_list, Hamil, Lab, nRow, .true., ICMax, LenHamil, tMC)

        allocate(Hamil(LenHamil), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error allocating Hamil.")
        allocate(Lab(LenHamil), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error allocating Lab.")

        Hamil = 0.0_dp
        Lab = 0

        call Detham(ndets, NEl, det_list, Hamil, Lab, nRow, .false., ICMax, LenHamil, tMC)

        nkry1 = nkry+1
        nblock = min(nexcit, nblk)
        len_scr = max(ndets*nexcit, 8*nblock*nkry)
        len_iscr = 6*nblock*nkry

        allocate(A_Arr(nexcit, nexcit), stat=ierr)
        allocate(V(ndets*nblock*nkry1), stat=ierr)
        allocate(AM(nblock*nblock*nkry1), stat=ierr)
        allocate(BM(nblock*nblock*nkry), stat=ierr)
        allocate(T(3*nblock*nkry*nblock*nkry), stat=ierr)
        allocate(WT(nblock*nkry), stat=ierr)
        allocate(scr(len_scr), stat=ierr)
        allocate(iscr(len_iscr), stat=ierr)
        allocate(ind(nexcit), stat=ierr)
        allocate(evecs_space(ndets, nexcit), stat=ierr)

        A_Arr = 0.0_dp
        V = 0.0_dp
        AM = 0.0_dp
        BM = 0.0_dp
        T = 0.0_dp
        WT = 0.0_dp
        scr = 0.0_dp
        iscr(1:len_iscr) = 0
        ind(1:nexcit) = 0
        evecs_space = 0.0_dp

        ! Perform Lanczos procedure.
        call neci_frsblkh(ndets, ICMax, nexcit, Hamil, Lab, evecs, evecs_space, nkry, nkry1, nblock, nrow, &
                           len_scr, len_iscr, A_Arr, evals, V, AM, BM, T, WT, scr, iscr, ind, ncycle, b2l, &
                           .false., .false., .false., .true.)

        ! The above routine returns *minus* the eigenvalues. Remove this factor:
        evals = -evals

        deallocate(evecs_space, &
                   A_Arr, &
                   V, &
                   BM, &
                   T, &
                   WT, &
                   scr, &
                   iscr, &
                   ind, &
                   AM, &
                   nrow, &
                   Lab, &
                   Hamil)

    end subroutine frsblk_wrapper

end module lanczos_wrapper
