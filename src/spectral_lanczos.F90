#include "macros.h"

module spectral_lanczos

    use bit_rep_data, only: NIfTot
    use bit_reps, only: test_flag
    use CalcData, only: NMCyc, tSemiStochastic
    use constants
    use DetBitOps, only: DetBitLT
    use FciMCData, only: HFDet, ilutHF, iHFProc, CurrentDets, core_hamiltonian, &
                         determ_proc_sizes, determ_space_size, partial_determ_vector, &
                         TotWalkers
    use ftlm_neci, only: subspace_expansion_lanczos, calc_final_hamil_elem
    use Parallel_neci, only: iProcIndex, nProcessors, MPIAllGatherV, MPIAllGather, &
                             MPISumAll
    use ParallelHelper, only: root
    use sparse_arrays, only: calculate_sparse_hamiltonian_parallel, sparse_ham

    implicit none

    real(dp), allocatable :: sl_hamil(:,:), sl_vecs(:,:)
    real(dp), allocatable :: full_vec_sl(:)

    integer(MPIArg), allocatable :: ndets_sl(:), disps_sl(:)
    integer :: n_lanc_vecs_sl, nomega_sl
    real(dp) :: broadening_sl

    real(dp), allocatable :: sl_h_eigv(:), sl_overlaps(:), sl_trans_amps(:)
    real(dp) :: allnorm_sl, allnorm_pert_sl
    real(dp) :: delta_omega_sl

    integer :: sl_unit
    integer :: sl_elem, sl_bit

contains

    subroutine perform_spectral_lanczos()

        integer :: i

        call init_spectral_lanczos()

        do i = 1, n_lanc_vecs_sl-1
            call subspace_expansion_lanczos(i, sl_vecs, full_vec_sl, sl_hamil, ndets_sl, disps_sl)
            write(6,'(1x,a19,1x,i3)') "Iteration complete:", i
            call neci_flush(6)
        end do

        call calc_final_hamil_elem(sl_vecs, full_vec_sl, sl_hamil, ndets_sl, disps_sl)

        call subspace_extraction_sl()

        write(6,'(1x,a60,/)') "Spectral Lanczos calculation complete. Outputting results..."
        call neci_flush(6)

        call output_spectral_lanczos()

        call end_spectral_lanczos()

    end subroutine perform_spectral_lanczos

    subroutine init_spectral_lanczos()

        use bit_rep_data, only: NIfDBO, extract_sign
        use DetBitOps, only: DetBitEq, EncodeBitDet, ilut_lt, ilut_gt
        use gndts_mod, only: gndts
        use hash, only: DetermineDetNode
        use sort_mod, only: sort
        use SystemData, only: nbasis, nel, BRR, nBasisMax, G1, tSpn, LMS, tParity, SymRestrict
        use util_mod, only: choose, get_free_unit

        integer :: ndets_this_proc, ndets_tot
        integer(MPIArg) :: mpi_temp
        integer, allocatable :: nI_list(:,:)
        integer(n_int), allocatable :: ilut_list(:,:)
        integer(n_int) :: ilut(0:NIfTot)
        character(len=*), parameter :: t_r = 'init_spectral_lanczos'
        integer :: i, j, ndets, proc, ierr, temp(1,1), hf_ind
        real(dp) :: real_sign(lenof_sign)
        real(dp) :: norm, norm_pert

        write(6,'(/,1x,a39,/)') "Beginning spectral Lanczos calculation."
        call neci_flush(6)

        sl_unit = get_free_unit()
        open(sl_unit, file='SPECTRAL_DATA',status='replace')

        allocate(ndets_sl(0:nProcessors-1))
        allocate(disps_sl(0:nProcessors-1))

        write(6,'(1x,a56)',advance='no') "Enumerating and storing all determinants in the space..."
        call neci_flush(6)

        ! Determine the total number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)
        allocate(nI_list(nel, ndets))
        ! Now generate them again and store them this time.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        ! Count the number of determinants belonging to this processor.
        ndets_this_proc = 0
        do i = 1, ndets
            write(6,*) nI_list(:,i)
            proc = DetermineDetNode(nI_list(:,i), 0)
            if (proc == iProcIndex) ndets_this_proc = ndets_this_proc + 1
        end do

        allocate(ilut_list(0:NIfTot, ndets_this_proc))
        ilut_list = 0_n_int

        j = 0
        do i = 1, ndets
            proc = DetermineDetNode(nI_list(:,i), 0)
            if (proc == iProcIndex) then
                call EncodeBitDet(nI_list(:,i), ilut)
                j = j + 1
                ilut_list(0:NIfDBO, j) = ilut(0:NIfDBO)
            end if
        end do

        call sort(ilut_list, ilut_lt, ilut_gt)

        write(6,'(1x,a9)') "Complete."
        call neci_flush(6)

        mpi_temp = int(ndets_this_proc, MPIArg)
        call MPIAllGather(mpi_temp, ndets_sl, ierr)

        disps_sl(0) = 0
        do i = 1, nProcessors-1
            disps_sl(i) = sum(ndets_sl(:i-1))
        end do

        ndets_tot = int(sum(ndets_sl), sizeof_int)

        write(6,'(1x,a44)',advance='no') "Allocating arrays to hold Lanczos vectors..."
        call neci_flush(6)
        allocate(sl_vecs(ndets_this_proc, n_lanc_vecs_sl))
        sl_vecs = 0.0_dp
        allocate(full_vec_sl(ndets_tot))
        write(6,'(1x,a9)') "Complete."
        call neci_flush(6)

        allocate(sl_hamil(n_lanc_vecs_sl, n_lanc_vecs_sl))
        allocate(sl_h_eigv(n_lanc_vecs_sl))
        allocate(sl_overlaps(n_lanc_vecs_sl))
        allocate(sl_trans_amps(n_lanc_vecs_sl))
        sl_hamil = 0.0_dp
        sl_h_eigv = 0.0_dp
        sl_overlaps = 0.0_dp
        sl_trans_amps = 0.0_dp

        write(6,'(1x,a48)') "Allocating and calculating Hamiltonian matrix..."
        call neci_flush(6)
        call calculate_sparse_hamiltonian_parallel(ndets_sl, ilut_list, .true.)
        write(6,'(1x,a48,/)') "Hamiltonian allocation and calculation complete."
        call neci_flush(6)

        sl_elem = (BRR(1)-1)/bits_n_int
        sl_bit = mod(BRR(1)-1, bits_n_int)

        write(6,*) "CurrentDets:"
        do i = 1, int(TotWalkers, sizeof_int)
            call extract_sign(CurrentDets(:, i), real_sign)
            write(6,'(i7, i12, 4x, f18.7, 4x, f18.7)') i, CurrentDets(0, i), real_sign
        end do

        write(6,*) "ilut_list:"
        do i = 1, ndets_this_proc
            call extract_sign(ilut_list(:, i), real_sign)
            write(6,'(i7, i12, 4x, f18.7, 4x, f18.7)') i, ilut_list(0, i), real_sign
        end do

        ! Copy the the initial vector to sl_vecs(:,1), apply the perturbation
        ! operator and normalise it.
        j = 0
        norm = 0.0_dp
        norm_pert = 0.0_dp
        do i = 1, int(TotWalkers, sizeof_int)
            call extract_sign(CurrentDets(:,i), real_sign)
            ! If the POPSFILE was generated from a calculation with the linear-scaling
            ! algorithm then it is possible to have the same determinant twice in a
            ! row, once with zero weight. If this is the case, skip the second case.
            if (i > 1) then
                if (DetBitEq(CurrentDets(:,i-1), CurrentDets(:,i), NIfDBO)) cycle
            end if
            norm = norm + real_sign(1)*real_sign(1)
            do
                j = j + 1
                !write(6,*) "j:", j, "CurrentDets:", CurrentDets(:,i), "ilut_list:", ilut_list(:,j)
                call neci_flush(6)
                if (DetBitEq(CurrentDets(:,i), ilut_list(:,j), NIfDBO)) then
                    if (btest(ilut_list(sl_elem, j), sl_bit)) then
                        sl_vecs(j, 1) = real_sign(1)
                        norm_pert = norm_pert + real_sign(1)*real_sign(1)
                    end if
                    exit
                end if
            end do
        end do
        call MPISumAll(norm, allnorm_sl)
        allnorm_sl = sqrt(allnorm_sl)
        call MPISumAll(norm_pert, allnorm_pert_sl)
        allnorm_pert_sl = sqrt(allnorm_pert_sl)
        sl_vecs(:,1) = sl_vecs(:,1)/allnorm_pert_sl

        deallocate(ilut_list)

    end subroutine init_spectral_lanczos

    subroutine subspace_extraction_sl()

        integer :: lwork, info, i
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: temp_overlaps(:)

        if (iProcIndex /= root) return

        ! Scrap space for the diagonaliser.
        lwork = max(1,3*n_lanc_vecs_sl-1)
        allocate(work(lwork))

        ! This routine diagonalises a symmetric matrix, A.
        ! V tells the routine to calculate eigenvalues *and* eigenvectors.
        ! U tells the routine to get the upper half of A (it is symmetric).
        ! n_lanc_vecs_sl is the number of rows and columns in A.
        ! A = sl_hamil. This matrix stores the eigenvectors in its columns on output.
        ! n_lanc_vecs_sl is the leading dimension of A.
        ! sl_h_eigv stores the eigenvalues on output.
        ! work is scrap space.
        ! lwork is the length of the work array.
        ! info = 0 on output if diagonalisation is successful.
        call dsyev('V', 'U', n_lanc_vecs_sl, sl_hamil, n_lanc_vecs_sl, sl_h_eigv, work, lwork, info)

        deallocate(work)

        ! Calculate the overlap of the perturbed ground state with each Lanczos vector.
        allocate(temp_overlaps(ndets_sl(iProcIndex)))
        ! temp_overlaps will hold the overlaps between the first Lanczos vector and all the
        ! other Lanczos vectors from this processor only.
        temp_overlaps = matmul(sl_vecs(:,1), sl_vecs)
        ! ...And now from all processors.
        call MPISumAll(temp_overlaps, sl_overlaps)
        ! The first Lanczos vector differs from the perturbed ground state only by a
        ! normalisation factor, which is accounted for here.
        sl_overlaps = sl_overlaps*allnorm_pert_sl/allnorm_sl
        deallocate(temp_overlaps)

        ! Now calculate the overlaps between the final Lanczos eigenvectors and the original
        ! perturbed ground state, which are the transition amplitudes which decide the
        ! height of the peaks in the spectrum.
        sl_trans_amps = matmul(sl_overlaps, sl_hamil)

        ! Output all of the eigenvalues and transition amplitudes.
        do i = 1, n_lanc_vecs_sl
            write(sl_unit,'(1x,f15.10,5x,f15.10)') sl_h_eigv(i), sl_trans_amps(i)
        end do
        write(sl_unit,'()')

    end subroutine subspace_extraction_sl
    
    subroutine output_spectral_lanczos()

        integer :: i, j
        real(dp) :: omega, spectral_weight

        if (iProcIndex /= root) return

        write(6,'(1x,a5,18X,a15)') "Omega", "Spectral_weight"

        omega = 0.0_dp
        do i = 1, nomega_sl
            spectral_weight = 0.0_dp
            do j = 2, n_lanc_vecs_sl
                spectral_weight = spectral_weight + &
                    (sl_trans_amps(j)**2*broadening_sl)/(pi*(broadening_sl**2 + (sl_h_eigv(1)-sl_h_eigv(j)+omega)**2))
            end do
            write(6,'(f18.12, 4x, f18.12)') omega, spectral_weight
            omega = omega + delta_omega_sl
        end do

    end subroutine output_spectral_lanczos

    subroutine end_spectral_lanczos()

        close(sl_unit)

        deallocate(sl_vecs)
        deallocate(full_vec_sl)
        deallocate(sl_overlaps)
        deallocate(sl_trans_amps)
        deallocate(sl_hamil)
        deallocate(sl_h_eigv)
        deallocate(ndets_sl)
        deallocate(disps_sl)

        write(6,'(/,1x,a30,/)') "Exiting Spectral Lanczos code."

    end subroutine end_spectral_lanczos

end module spectral_lanczos
