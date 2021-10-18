#include "macros.h"

module exact_diag

    use constants
    use FciMCData, only: hamiltonian
    use SystemData, only: t_non_hermitian, tGUGA, nSpatOrbs
    use guga_bitRepOps, only: CSF_Info_t

    implicit none

    integer :: ndets_ed
    real(dp), allocatable :: eigv_ed(:)

contains

    subroutine perform_exact_diag_all_symmetry()

        integer :: lwork, info
        real(dp), allocatable :: work(:)

        call init_exact_diag()

        ! Create the workspace for dsyev.
        lwork = max(1, 3 * ndets_ed - 1)
        allocate(work(lwork))

        write(stdout, '(1x,a28)', advance='no') "Diagonalising Hamiltonian..."
        call neci_flush(6)

        ! Perform the diagonalisation.
        call dsyev('V', 'U', ndets_ed, hamiltonian, ndets_ed, eigv_ed, work, lwork, info)

        write(stdout, '(1x,a9,/)') "Complete."
        call neci_flush(6)

        call output_exact_spectrum()

        deallocate(work)

        call end_exact_spectrum()

    end subroutine perform_exact_diag_all_symmetry

    subroutine init_exact_diag()

        use bit_rep_data, only: NIfTot
        use gndts_mod, only: gndts_all_sym_this_proc
        use SystemData, only: nbasis, nel
        use util_mod, only: choose

        integer :: expected_ndets_tot
        integer(n_int), allocatable :: ilut_list(:, :)
        character(len=*), parameter :: t_r = 'init_exact_diag'
        integer :: ierr

        write(stdout, '(/,1x,a57,/)') "Beginning exact diagonalisation in all symmetry sectors."
        call neci_flush(6)

        write(stdout, '(1x,a56)', advance='no') "Enumerating and storing all determinants in the space..."
        call neci_flush(6)

        ! Generate and count all the determinants on this processor, but don't store them.
        call gndts_all_sym_this_proc(ilut_list, .true., ndets_ed)
        allocate(ilut_list(0:NIfTot, ndets_ed))
        ! Now generate them again and store them this time.
        call gndts_all_sym_this_proc(ilut_list, .false., ndets_ed)

        write(stdout, '(1x,a9)') "Complete."
        call neci_flush(6)

        expected_ndets_tot = int(choose(nbasis, nel))
        if (ndets_ed /= expected_ndets_tot) then
            write(stdout, *) "ndets counted:", ndets_ed, "ndets expected:", expected_ndets_tot
            call stop_all('t_r', 'The number of determinants generated is not &
                                    &consistent with the expected number.')
        end if

        allocate(eigv_ed(ndets_ed), stat=ierr)
        if (ierr /= 0) then
            write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating eigenvalue array.")
        end if
        eigv_ed = 0.0_dp

        write(stdout, '(1x,a48)') "Allocating and calculating Hamiltonian matrix..."
        call neci_flush(6)
        allocate(hamiltonian(ndets_ed, ndets_ed), stat=ierr)
        if (ierr /= 0) then
            write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating Hamiltonian array.")
        end if
        write(stdout, '(1x,a46)') "Hamiltonian allocation completed successfully."
        call neci_flush(6)
        call calculate_full_hamiltonian(ilut_list, hamiltonian)
        write(stdout, '(1x,a33)') "Hamiltonian calculation complete."
        call neci_flush(6)

    end subroutine init_exact_diag

    subroutine calculate_full_hamiltonian(ilut_list, local_hamil)

        use bit_reps, only: decode_bit_det
        use Determinants, only: get_helement
        use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
        use SystemData, only: tHPHF, nel
        use guga_excitations, only: calc_guga_matrix_element
        use guga_data, only: ExcitationInformation_t
        use guga_bitrepops, only: new_CSF_Info_t, fill_csf_i
        use bit_rep_data, only: nifd
        integer(n_int), intent(in) :: ilut_list(0:, :)
        HElement_t(dp), intent(inout), allocatable :: local_hamil(:, :)

        integer :: ndets, i, j, ierr
        integer :: nI(nel), nJ(nel)
        character(*), parameter :: t_r = "calculate_full_hamiltonian"
        type(ExcitationInformation_t) :: excitInfo
        type(CSF_Info_t) :: csf_i

        ! Initial checks that arrays passed in are consistent.
        ndets = size(ilut_list, 2)
        if (allocated(local_hamil)) then
            if (size(local_hamil, 1) /= ndets) call stop_all(t_r, "Inconsistent sizes Hamiltonian and ilut arrays.")
        else
            allocate(local_hamil(ndets, ndets), stat=ierr)
            if (ierr /= 0) then
                write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating Hamiltonian array.")
            end if
        end if

        ! Loop over every pair of determinants and calculate all elements.
        if (t_non_hermitian) then
            call stop_all(t_r, "check this for non-hermitian hamil!")
        end if

        call new_CSF_Info_t(nSpatOrbs, csf_i)
        do i = 1, ndets
            call decode_bit_det(nI, ilut_list(:, i))
            if (tGUGA) call fill_csf_i(ilut_list(0:nifd, i), csf_i)

            do j = i, ndets
                call decode_bit_det(nJ, ilut_list(:, j))
                if (i == j) then
                    if (tHPHF) then
                        local_hamil(i, i) = hphf_diag_helement(nI, ilut_list(:, i))
                    else
                        local_hamil(i, i) = get_helement(nI, nI, 0)
                    end if
                else
                    if (tHPHF) then
                        local_hamil(i, j) = hphf_off_diag_helement(nI, nJ, ilut_list(:, i), ilut_list(:, j))
                    else if (tGUGA) then
                        call calc_guga_matrix_element(ilut_list(:, i), csf_i, ilut_list(:, j), &
                                                      excitInfo, local_hamil(i, j), .true., 1)
! #ifdef CMPLX_
!                         local_hamil(i,j) = conjg(local_hamil(i,j))
! #endif
                    else
                        local_hamil(i, j) = get_helement(nI, nJ, ilut_list(:, i), ilut_list(:, j))
                    end if
                end if
                local_hamil(j, i) = local_hamil(i, j)
            end do
        end do

    end subroutine calculate_full_hamiltonian

    subroutine output_exact_spectrum()

        use util_mod, only: get_free_unit

        integer :: i, temp_unit

        temp_unit = get_free_unit()
        open(temp_unit, file='EIGENVALUES', status='replace')
        write(temp_unit, '(1x,a12)') "Total energy"
        do i = 1, ndets_ed
            write(temp_unit, '(1x,i7,5x,f15.10)') i, eigv_ed(i)
        end do
        close(temp_unit)

    end subroutine output_exact_spectrum

    subroutine end_exact_spectrum()

        deallocate(eigv_ed)
        deallocate(hamiltonian)

    end subroutine end_exact_spectrum

end module exact_diag
