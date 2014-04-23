#include "macros.h"

module exact_diag

    use constants

    implicit none

contains

    subroutine calculate_full_hamiltonian(ilut_list, local_hamil)

        use bit_reps, only: decode_bit_det
        use Determinants, only: get_helement
        use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
        use SystemData, only: tHPHF, nel

        integer(n_int), intent(in) :: ilut_list(0:,:)
        real(dp), intent(inout), allocatable :: local_hamil(:,:)

        integer :: ndets, i, j, ierr
        integer :: nI(nel), nJ(nel)
        character(*), parameter :: t_r = "calculate_full_hamiltonian"

        ! Initial checks that arrays passed in are consistent.
        ndets = size(ilut_list, 2)
        if (allocated(local_hamil)) then
            if(size(local_hamil,1) /= ndets) call stop_all(t_r, "Inconsistent sizes Hamiltonian and ilut arrays.")
        else
            allocate(local_hamil(ndets, ndets), stat=ierr)
            if (ierr /= 0) then
                write(6,'(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating Hamiltonian array.")
            end if
        end if

        ! Loop over every pair of determinants and calculate all elements.
        do i = 1, ndets
            call decode_bit_det(nI, ilut_list(:,i))
            do j = i, ndets
                call decode_bit_det(nJ, ilut_list(:,j))
                if (i == j) then
                    if (tHPHF) then
                        local_hamil(i,i) = hphf_diag_helement(nI, ilut_list(:,i))
                    else
                        local_hamil(i,i) = get_helement(nI, nI, 0)
                    end if
                else
                    if (tHPHF) then
                        local_hamil(i,j) = hphf_off_diag_helement(nI, nJ, ilut_list(:,i), ilut_list(:,j))
                    else
                        local_hamil(i,j) = get_helement(nI, nJ, ilut_list(:,i), ilut_list(:,j))
                    end if
                end if
                local_hamil(j,i) = local_hamil(i,j)
            end do
        end do

    end subroutine calculate_full_hamiltonian

end module exact_diag
