#include "macros.h"

module exact_spectrum

    use constants
    use FciMCData, only: hamiltonian
    use spectral_lanczos, only: nomega_spectral, delta_omega_spectral, spectral_broadening
    use spectral_lanczos, only: spectral_ground_energy, tIncludeGroundSpectral, min_omega_sl

    implicit none

    integer :: ndets_es
    real(dp), allocatable :: perturbed_ground_es(:)
    real(dp), allocatable :: transition_amps_es(:)
    real(dp), allocatable :: eigv_es(:)

contains

    subroutine get_exact_spectrum()

        integer :: lwork, info
        real(dp), allocatable :: work(:)

        call init_exact_spectrum()

        ! Create the workspace for dsyev.
        lwork = max(1,3*ndets_es-1)
        allocate(work(lwork))

        write(6,'(1x,a28)',advance='no') "Diagonalising Hamiltonian..."
        call neci_flush(6)

        ! Perform the diagonalisation.
        call dsyev('V', 'U', ndets_es, hamiltonian, ndets_es, eigv_es, work, lwork, info)

        write(6,'(1x,a9,/)') "Complete."
        call neci_flush(6)

        transition_amps_es = matmul(perturbed_ground_es, hamiltonian)

        call output_exact_spectrum()

        deallocate(work)

        call end_exact_spectrum()

    end subroutine get_exact_spectrum

    subroutine init_exact_spectrum()

        use bit_rep_data, only: NIfTot, NIfDBO, extract_sign
        use bit_reps, only: decode_bit_det
        use CalcData, only: pops_norm
        use DetBitOps, only: DetBitEq, EncodeBitDet, ilut_lt, ilut_gt
        use exact_diag, only: calculate_full_hamiltonian
        use FciMCData, only: TotWalkers, CurrentDets
        use gndts_mod, only: gndts
        use sort_mod, only: sort
        use SystemData, only: nbasis, nel, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict

        integer, allocatable :: nI_list(:,:)
        integer(n_int), allocatable :: ilut_list(:,:)
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: i, j, hf_ind, temp(1,1), ierr
        real(dp) :: real_sign(lenof_sign)
        character(len=*), parameter :: t_r = 'init_exact_spectrum'

        write(6,'(/,1x,a48,/)') "Beginning calculation of exact spectral density."
        call neci_flush(6)

        write(6,'(1x,a56)',advance='yes') "Enumerating and storing all determinants in the space..."
        call neci_flush(6)

        ! Calculate the number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets_es, hf_ind)
        allocate(nI_list(nel, ndets_es))
        ! Now generate them again and store them this time.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets_es, hf_ind)

        write(6,'(1x,a33,i7)') "Number of determinants in space: ", ndets_es

        allocate(ilut_list(0:NIfTot, ndets_es))
        ilut_list = 0_n_int

        do i = 1, ndets_es
            call EncodeBitDet(nI_list(:,i), ilut)
            ilut_list(0:NIfDBO, i) = ilut(0:NIfDBO)
        end do

        call sort(ilut_list, ilut_lt, ilut_gt)

        allocate(eigv_es(ndets_es), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating eigenvalue array.")
        end if
        eigv_es = 0.0_dp

        allocate(transition_amps_es(ndets_es), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold transition amplitudes.")
        end if
        transition_amps_es = 0.0_dp

        allocate(perturbed_ground_es(ndets_es), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold perturbed ground state components.")
        end if
        perturbed_ground_es = 0.0_dp

        write(6,'(1x,a48)') "Allocating and calculating Hamiltonian matrix..."
        call neci_flush(6)
        allocate(hamiltonian(ndets_es, ndets_es), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating Hamiltonian array.")
        end if
        write(6,'(1x,a46)') "Hamiltonian allocation completed successfully."
        call neci_flush(6)
        call calculate_full_hamiltonian(ilut_list, hamiltonian)
        write(6,'(1x,a33)') "Hamiltonian calculation complete."
        call neci_flush(6)

        j = 0
        call sort(CurrentDets(:,1:TotWalkers), ilut_lt, ilut_gt)
        do i = 1, int(TotWalkers, sizeof_int)
            call extract_sign(CurrentDets(:,i), real_sign)
            ! If the POPSFILE was generated from a calculation with the linear-scaling
            ! algorithm then it is possible to have the same determinant twice in a
            ! row, once with zero weight. If this is the case, skip the second case.
            if (i > 1) then
                if (DetBitEq(CurrentDets(:,i-1), CurrentDets(:,i), NIfDBO)) cycle
            end if
            do
                j = j + 1
                if (DetBitEq(CurrentDets(:,i), ilut_list(:,j), NIfDBO)) then
                    perturbed_ground_es(j) = real_sign(1)
                    exit
                else if (j > ndets_es) then
                    write(6,*) "Determinant POPSFILE not found:", CurrentDets(:,i)
                    call decode_bit_det(nI, CurrentDets(:,i))
                    write(6,*) "Occupied orbitals in the above POPSFILE determinant:", nI
                    call stop_all(t_r, "Determinant in POPSFILE not found in the enumerated FCI space.")
                end if
            end do
        end do
        perturbed_ground_es = perturbed_ground_es/sqrt(pops_norm)

        deallocate(ilut_list)
        deallocate(nI_list)

    end subroutine init_exact_spectrum

    subroutine output_exact_spectrum()

        use util_mod, only: get_free_unit

        integer :: i, j, min_vec, temp_unit
        real(dp) :: omega, spectral_weight

        temp_unit = get_free_unit()
        open(temp_unit, file='SPECTRAL_DATA',status='replace')
        write(temp_unit,'(1x,a38)') "Eigenvalues and transition amplitudes:"
        do i = 1, ndets_es
            write(temp_unit,'(1x,i7,5x,f15.10,5x,f15.10)') i, eigv_es(i), transition_amps_es(i)
        end do
        close(temp_unit)

        write(6,'(1x,a5,18X,a15)') "Omega", "Spectral_weight"

        ! Do we include the ground state in the spectrum or not?
        if (tIncludeGroundSpectral) then
            min_vec = 1
        else
            min_vec = 2
        end if

        omega = min_omega_sl
        do i = 1, nomega_spectral + 1
            spectral_weight = 0.0_dp
            do j = min_vec, ndets_es
                spectral_weight = spectral_weight + &
                    (transition_amps_es(j)**2*spectral_broadening)/&
                    (pi*(spectral_broadening**2 + (spectral_ground_energy-eigv_es(j)+omega)**2))
            end do
            write(6,'(f18.12, 4x, f18.12)') omega, spectral_weight
            omega = omega + delta_omega_spectral
        end do

    end subroutine output_exact_spectrum

    subroutine end_exact_spectrum()

        deallocate(perturbed_ground_es)
        deallocate(transition_amps_es)
        deallocate(eigv_es)
        deallocate(hamiltonian)

    end subroutine end_exact_spectrum

end module exact_spectrum
