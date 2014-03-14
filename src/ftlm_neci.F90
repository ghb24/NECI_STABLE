#include "macros.h"

! This module contains the routines to perform the finite-temperature Lanczos
! method deterministically.

module ftlm_neci

    use bit_rep_data, only: NIfTot
    use bit_reps, only: test_flag
    use CalcData, only: NMCyc, tSemiStochastic
    use constants
    use DetBitOps, only: DetBitLT
    use FciMCData, only: HFDet, ilutHF, iHFProc, CurrentDets, core_hamiltonian, &
                         determ_proc_sizes, determ_space_size, partial_determ_vector, &
                         TotWalkers
    use Parallel_neci, only: iProcIndex, nProcessors, MPIAllGatherV, MPIAllGather, &
                             MPISumAll
    use ParallelHelper, only: root
    use sparse_arrays, only: calculate_sparse_hamiltonian_parallel, sparse_ham

    implicit none

    real(dp), allocatable :: ftlm_hamil(:,:), ftlm_vecs(:,:)
    real(dp), allocatable :: full_vec_ftlm(:)

    integer(MPIArg), allocatable :: ndets_ftlm(:), disps_ftlm(:)
    integer :: n_init_vecs_ftlm, n_lanc_vecs_ftlm, nbeta_ftlm

    real(dp), allocatable :: ftlm_h_eigv(:), ftlm_trace(:), ftlm_e_num(:)
    real(dp) :: delta_beta_ftlm

    integer :: eigv_unit

contains

    subroutine perform_ftlm()

        integer :: i, j

        call init_ftlm()

        do i = 1, n_init_vecs_ftlm

            write(6,'(1x,a28,1x,i3)') "Starting from initial vector", i
            call neci_flush(6)

            ftlm_vecs = 0.0_dp
            ftlm_hamil = 0.0_dp

            call gen_init_vec_ftlm()

            do j = 1, n_lanc_vecs_ftlm-1
                call subspace_expansion_ftlm(j)
                write(6,'(1x,a19,1x,i3)') "Iteration complete:", j
                call neci_flush(6)
            end do

            call calc_final_hamil_elem_ftlm()

            call subspace_extraction_ftlm()

            call add_in_contribs_to_energy()

            write(6,'(1x,a45,/)') "Calculation complete for this initial vector."
            call neci_flush(6)

        end do

        write(6,'(1x,a48,/)') "FTLM calculation complete. Outputting results..."
        call neci_flush(6)

        call output_ftlm()

        call end_ftlm()

    end subroutine perform_ftlm

    subroutine init_ftlm()

        use gndts_mod, only: gndts_all_sym_this_proc
        use SystemData, only: nbasis, nel
        use util_mod, only: choose, get_free_unit

        integer :: ndets_this_proc, ndets_tot, expected_ndets_tot
        integer(MPIArg) :: mpi_temp
        integer(n_int), allocatable :: ilut_list(:,:)
        character(len=*), parameter :: t_r = 'init_ftlm'
        integer :: i, ierr

        write(6,'(/,1x,a49,/)') "Beginning finite-temperature Lanczos calculation."
        call neci_flush(6)

        eigv_unit = get_free_unit()
        open(eigv_unit, file='FTLM_EIGV',status='replace')

        allocate(ndets_ftlm(0:nProcessors-1))
        allocate(disps_ftlm(0:nProcessors-1))

        write(6,'(1x,a56)',advance='no') "Enumerating and storing all determinants in the space..."
        call neci_flush(6)

        ! Generate and count all the determinants on this processor, but don't store them.
        call gndts_all_sym_this_proc(ilut_list, .true., ndets_this_proc)
        allocate(ilut_list(0:NIfTot, ndets_this_proc))
        ! Now generate them again and store them this time.
        call gndts_all_sym_this_proc(ilut_list, .false., ndets_this_proc)

        write(6,'(1x,a9)') "Complete."
        call neci_flush(6)

        mpi_temp = int(ndets_this_proc, MPIArg)
        call MPIAllGather(mpi_temp, ndets_ftlm, ierr)

        disps_ftlm(0) = 0
        do i = 1, nProcessors-1
            disps_ftlm(i) = sum(ndets_ftlm(:i-1))
        end do

        ndets_tot = int(sum(ndets_ftlm), sizeof_int)
        expected_ndets_tot = choose(nbasis, nel)
        if (ndets_tot /= expected_ndets_tot) then
            write(6,*) "ndets counted:", ndets_tot, "ndets expected:", expected_ndets_tot
            call stop_all('t_r', 'The number of determinants generated is not &
                                    &consistent with the expected number.')
        end if

        write(6,'(1x,a44)',advance='no') "Allocating arrays to hold Lanczos vectors..."
        call neci_flush(6)
        allocate(ftlm_vecs(ndets_this_proc, n_lanc_vecs_ftlm))
        allocate(full_vec_ftlm(ndets_tot))
        write(6,'(1x,a9)') "Complete."
        call neci_flush(6)

        allocate(ftlm_hamil(n_lanc_vecs_ftlm, n_lanc_vecs_ftlm))
        allocate(ftlm_trace(nbeta_ftlm+1))
        allocate(ftlm_e_num(nbeta_ftlm+1))
        allocate(ftlm_h_eigv(n_lanc_vecs_ftlm))
        ftlm_trace = 0.0_dp
        ftlm_e_num = 0.0_dp
        ftlm_h_eigv = 0.0_dp

        write(6,'(1x,a48)') "Allocating and calculating Hamiltonian matrix..."
        call neci_flush(6)
        call calculate_sparse_hamiltonian_parallel(ndets_ftlm, ilut_list, .true.)
        write(6,'(1x,a48,/)') "Hamiltonian allocation and calculation complete."
        call neci_flush(6)

    end subroutine init_ftlm

    subroutine gen_init_vec_ftlm()

        use dSFMT_interface, only: genrand_real2_dSFMT

        integer :: i
        real(dp) :: r, amp

        ! The magnitude of the amplitude required to make the initial vector normalised to 1.0.
        amp = 1/sqrt(real(sum(ndets_ftlm), dp))

        do i = 1, ndets_ftlm(iProcIndex)
            r = genrand_real2_dSFMT()
            if (r > 0.5_dp) then
                ftlm_vecs(i,1) = amp
            else
                ftlm_vecs(i,1) = -amp
            end if
        end do

    end subroutine gen_init_vec_ftlm

    subroutine subspace_expansion_ftlm(ivec)

        integer, intent(in) :: ivec
        integer :: i, j
        real(dp) :: overlap, tot_overlap, last_norm, norm, tot_norm

        call MPIAllGatherV(ftlm_vecs(:,ivec), full_vec_ftlm, ndets_ftlm, disps_ftlm)

        ! Multiply the last Lanczos vector by the Hamiltonian.
        ftlm_vecs(:,ivec+1) = 0.0_dp
        do i = 1, ndets_ftlm(iProcIndex)
            do j = 1, sparse_ham(i)%num_elements
                ftlm_vecs(i, ivec+1) = ftlm_vecs(i, ivec+1) + &
                    sparse_ham(i)%elements(j)*full_vec_ftlm(sparse_ham(i)%positions(j))
            end do
        end do

        overlap = dot_product(ftlm_vecs(:,ivec+1), ftlm_vecs(:,ivec))
        call MPISumAll(overlap, tot_overlap)

        if (ivec == 1) then
            ftlm_vecs(:,ivec+1) = ftlm_vecs(:,ivec+1) - tot_overlap*ftlm_vecs(:,ivec)
        else
            last_norm = ftlm_hamil(ivec, ivec-1)
            ftlm_vecs(:,ivec+1) = ftlm_vecs(:,ivec+1) - tot_overlap*ftlm_vecs(:,ivec) - last_norm*ftlm_vecs(:,ivec-1)
        end if

        norm = 0.0_dp
        do i = 1, ndets_ftlm(iProcIndex)
            norm = norm + ftlm_vecs(i,ivec+1)*ftlm_vecs(i,ivec+1)
        end do
        call MPISumAll(norm, tot_norm)
        tot_norm = sqrt(tot_norm)

        ! The final new Lanczos vector.
        ftlm_vecs(:,ivec+1) = ftlm_vecs(:,ivec+1)/tot_norm

        ! Update the subspace Hamiltonian. It can be shown that this has a tridiagonal form
        ! where the non-zero elements are equal to the following.
        ftlm_hamil(ivec, ivec) = tot_overlap
        ftlm_hamil(ivec, ivec+1) = tot_norm
        ftlm_hamil(ivec+1, ivec) = tot_norm

    end subroutine subspace_expansion_ftlm

    subroutine calc_final_hamil_elem_ftlm()

        integer :: i, j
        real(dp) :: temp, overlap, tot_overlap

        call MPIAllGatherV(ftlm_vecs(:,n_lanc_vecs_ftlm), full_vec_ftlm, ndets_ftlm, disps_ftlm)

        overlap = 0.0_dp
        do i = 1, ndets_ftlm(iProcIndex)
            ! If we denote the final Lanczos vector as V, then at the end of the following
            ! do loop, temp will hold (H*V)_i. That is, the element of H*V corresponding to
            ! the i'th basis vector, which we're currently on (where H is the Hamiltonian matrix).
            temp = 0.0_dp
            do j = 1, sparse_ham(i)%num_elements
                temp = temp + sparse_ham(i)%elements(j)*full_vec_ftlm(sparse_ham(i)%positions(j))
            end do
            overlap = overlap + temp*ftlm_vecs(i,n_lanc_vecs_ftlm)
        end do
        call MPISumAll(overlap, tot_overlap)

        ftlm_hamil(n_lanc_vecs_ftlm, n_lanc_vecs_ftlm) = tot_overlap

    end subroutine calc_final_hamil_elem_ftlm

    subroutine subspace_extraction_ftlm()

        integer :: lwork, info, i
        real(dp), allocatable :: work(:)

        if (iProcIndex /= root) return

        ! Scrap space for the diagonaliser.
        lwork = max(1,3*n_lanc_vecs_ftlm-1)
        allocate(work(lwork))

        ! This routine diagonalises a symmetric matrix, A.
        ! V tells the routine to calculate eigenvalues *and* eigenvectors.
        ! U tells the routine to get the upper half of A (it is symmetric).
        ! n_lanc_vecs_ftlm is the number of rows and columns in A.
        ! A = ftlm_hamil. This matrix stores the eigenvectors in its columns on output.
        ! n_lanc_vecs_ftlm is the leading dimension of A.
        ! ftlm_h_eigv stores the eigenvalues on output.
        ! work is scrap space.
        ! lwork is the length of the work array.
        ! info = 0 on output if diagonalisation is successful.
        call dsyev('V', 'U', n_lanc_vecs_ftlm, ftlm_hamil, n_lanc_vecs_ftlm, ftlm_h_eigv, work, lwork, info)

        deallocate(work)

        ! Output all of the eigenvalues.
        do i = 1, n_lanc_vecs_ftlm
            write(eigv_unit,'(1x,f15.10)',advance='no') ftlm_h_eigv(i)
        end do
        write(eigv_unit,'()')

    end subroutine subspace_extraction_ftlm
    
    subroutine add_in_contribs_to_energy()

        integer :: i, j
        real(dp) :: beta

        if (iProcIndex /= root) return

        do i = 1, n_lanc_vecs_ftlm
            beta = 0.0_dp
            do j = 1, nbeta_ftlm + 1
                ftlm_trace(j) = ftlm_trace(j) + (ftlm_hamil(1,i)**2)*exp(-beta*ftlm_h_eigv(i))
                ftlm_e_num(j) = ftlm_e_num(j) + (ftlm_hamil(1,i)**2)*ftlm_h_eigv(i)*exp(-beta*ftlm_h_eigv(i))
                beta = beta + delta_beta_ftlm
            end do
        end do

    end subroutine add_in_contribs_to_energy

    subroutine output_ftlm()

        integer :: i
        real(dp) :: beta

        if (iProcIndex /= root) return

        write(6,'(1x,a4,18X,a11,11x,a11)') "Beta", "E_numerator", "Denominator"

        beta = 0.0_dp
        do i = 1, nbeta_ftlm + 1
            write(6,'(es17.10,5x,es17.10,5x,es17.10)') beta, ftlm_e_num(i), ftlm_trace(i)
            beta = beta + delta_beta_ftlm
        end do

    end subroutine output_ftlm

    subroutine end_ftlm()

        close(eigv_unit)

        deallocate(ftlm_vecs)
        deallocate(full_vec_ftlm)
        deallocate(ftlm_hamil)
        deallocate(ftlm_trace)
        deallocate(ftlm_e_num)
        deallocate(ftlm_h_eigv)
        deallocate(ndets_ftlm)
        deallocate(disps_ftlm)

        write(6,'(/,1x,a18,/)') "Exiting FTLM code."

    end subroutine end_ftlm

end module ftlm_neci
