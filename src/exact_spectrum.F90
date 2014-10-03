#include "macros.h"

module exact_spectrum

    use constants
    use FciMCData, only: hamiltonian, perturbation
    use PopsfileMod, only: read_popsfile_wrapper
    use spectral_data
    use spectral_lanczos, only: output_spectrum, return_perturbed_ground_spec

    implicit none

    real(dp), allocatable :: eigv_es(:)

contains

    subroutine get_exact_spectrum()

        integer :: lwork, info, ndets
        real(dp), allocatable :: work(:)
        ! Data for the testsuite to use.
        real(dp) :: spec_low, spec_high

        call init_exact_spectrum(ndets)

        ! Create the workspace for dsyev.
        lwork = max(1,3*ndets-1)
        allocate(work(lwork))

        write(6,'(1x,a28)',advance='no') "Diagonalising Hamiltonian..."
        call neci_flush(6)

        ! Perform the diagonalisation.
        call dsyev('V', 'U', ndets, hamiltonian, ndets, eigv_es, work, lwork, info)

        write(6,'(1x,a9,/)') "Complete."
        call neci_flush(6)

        trans_amps_right = matmul(pert_ground_right, hamiltonian)
        trans_amps_left = matmul(pert_ground_left, hamiltonian)

        call output_spectrum(ndets, eigv_es, spec_low, spec_high)

        call write_exact_spec_testsuite_data(eigv_es(1), eigv_es(ndets), spec_low, spec_high)

        call end_exact_spectrum()

        deallocate(work)

    end subroutine get_exact_spectrum

    subroutine init_exact_spectrum(ndets)

        use bit_rep_data, only: NIfTot, NIfDBO
        use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt
        use exact_diag, only: calculate_full_hamiltonian
        use gndts_mod, only: gndts
        use sort_mod, only: sort
        use SystemData, only: nbasis, nel, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict

        integer, intent(out) :: ndets
        integer, allocatable :: nI_list(:,:)
        integer(n_int), allocatable :: ilut_list(:,:)
        integer(n_int) :: ilut(0:NIfTot)
        integer :: i, hf_ind, temp(1,1), ierr
        character(len=*), parameter :: t_r = 'init_exact_spectrum'

        write(6,'(/,1x,a48,/)') "Beginning calculation of exact spectral density."
        call neci_flush(6)

        write(6,'(1x,a56)',advance='yes') "Enumerating and storing all determinants in the space..."
        call neci_flush(6)

        ! Calculate the number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)
        allocate(nI_list(nel, ndets))
        ! Now generate them again and store them this time.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        write(6,'(1x,a33,i7)') "Number of determinants in space: ", ndets

        allocate(ilut_list(0:NIfTot, ndets))
        ilut_list = 0_n_int

        do i = 1, ndets
            call EncodeBitDet(nI_list(:,i), ilut)
            ilut_list(0:NIfDBO, i) = ilut(0:NIfDBO)
        end do

        call sort(ilut_list, ilut_lt, ilut_gt)

        allocate(eigv_es(ndets), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating eigenvalue array.")
        end if
        eigv_es = 0.0_dp

        allocate(trans_amps_left(ndets), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold left transition amplitudes.")
        end if
        trans_amps_left = 0.0_dp

        allocate(trans_amps_right(ndets), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold right transition amplitudes.")
        end if
        trans_amps_right = 0.0_dp

        allocate(pert_ground_left(ndets), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold left perturbed ground state components.")
        end if

        allocate(pert_ground_right(ndets), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold right perturbed ground state components.")
        end if

        write(6,'(1x,a48)') "Allocating and calculating Hamiltonian matrix..."
        call neci_flush(6)
        allocate(hamiltonian(ndets, ndets), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating Hamiltonian array.")
        end if
        write(6,'(1x,a46)') "Hamiltonian allocation completed successfully."
        call neci_flush(6)
        call calculate_full_hamiltonian(ilut_list, hamiltonian)
        write(6,'(1x,a33)') "Hamiltonian calculation complete."
        call neci_flush(6)

        call return_perturbed_ground_spec(left_perturb_spectral, ilut_list, pert_ground_left, left_pert_norm)
        call return_perturbed_ground_spec(right_perturb_spectral, ilut_list, pert_ground_right, right_pert_norm)

        deallocate(ilut_list)
        deallocate(nI_list)

    end subroutine init_exact_spectrum

    subroutine write_exact_spec_testsuite_data(eigv_low, eigv_high, spec_low, spec_high)

        real(dp), intent(in) :: eigv_low, eigv_high, spec_low, spec_high

        write(6,'(/,1X,64("="))')
        write(6,'(1X,"Exact spectrum testsuite data:")')
        write(6,'(1X,"Lowest eigenvalue of H:",21X,es20.13)') eigv_low
        write(6,'(1X,"Highest eigenvalue of H:",20X,es20.13)') eigv_high
        write(6,'(1X,"Spectral weight at the lowest omega value:",2X,es20.13)') spec_low
        write(6,'(1X,"Spectral weight at the highest omega value:",1X,es20.13)') spec_high
        write(6,'(1X,64("="))')

    end subroutine write_exact_spec_testsuite_data

    subroutine end_exact_spectrum()

        deallocate(pert_ground_left)
        deallocate(pert_ground_right)
        deallocate(trans_amps_left)
        deallocate(trans_amps_right)
        deallocate(eigv_es)
        deallocate(hamiltonian)

    end subroutine end_exact_spectrum

end module exact_spectrum
