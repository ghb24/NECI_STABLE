#include "macros.h"

module spectral_lanczos

    use bit_rep_data, only: NIfTot
    use bit_reps, only: test_flag
    use CalcData, only: NMCyc, tSemiStochastic, pops_norm
    use constants
    use DetBitOps, only: DetBitLT
    use FciMCData, only: HFDet, ilutHF, iHFProc, CurrentDets, determ_proc_sizes, &
                         determ_space_size, partial_determ_vector, TotWalkers
    use ftlm_neci, only: subspace_expansion_lanczos, calc_final_hamil_elem
    use Parallel_neci, only: iProcIndex, nProcessors, MPIAllGather, MPISumAll, MPIBCast
    use ParallelHelper, only: root
    use PopsfileMod, only: read_popsfile_wrapper
    use sparse_arrays, only: calculate_sparse_hamiltonian_parallel, sparse_ham
    use spectral_data

    implicit none

    real(dp), allocatable :: sl_hamil(:,:), sl_vecs(:,:)
    real(dp), allocatable :: full_vec_sl(:)

    integer(MPIArg), allocatable :: ndets_sl(:), disps_sl(:)
    integer :: n_lanc_vecs_sl

    real(dp), allocatable :: sl_h_eigv(:), sl_overlaps(:)

contains

    subroutine perform_spectral_lanczos()

        integer :: i
        ! Data for the testsuite to use.
        real(dp) :: spec_low, spec_high

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

        call output_spectrum(n_lanc_vecs_sl, sl_h_eigv, spec_low, spec_high)

        call write_spec_lanc_testsuite_data(spec_low, spec_high)

        call end_spectral_lanczos()

    end subroutine perform_spectral_lanczos

    subroutine init_spectral_lanczos()

        use bit_rep_data, only: NIfDBO, extract_sign
        use DetBitOps, only: DetBitEq, EncodeBitDet, ilut_lt, ilut_gt
        use gndts_mod, only: gndts
        use hash, only: DetermineDetNode
        use sort_mod, only: sort
        use SystemData, only: nbasis, nel, BRR, nBasisMax, G1, tSpn, LMS, tParity, SymRestrict
        use util_mod, only: choose

        integer :: ndets_this_proc, ndets_tot
        integer(MPIArg) :: mpi_temp
        integer, allocatable :: nI_list(:,:)
        integer(n_int), allocatable :: ilut_list(:,:)
        integer(n_int) :: ilut(0:NIfTot)
        character(len=*), parameter :: t_r = 'init_spectral_lanczos'
        integer :: i, j, ndets, proc, ierr, temp(1,1), hf_ind
        real(dp) :: real_sign(lenof_sign)
        real(dp) :: norm_pert

        write(6,'(/,1x,a39,/)') "Beginning spectral Lanczos calculation."
        call neci_flush(6)

        allocate(ndets_sl(0:nProcessors-1))
        allocate(disps_sl(0:nProcessors-1))

        write(6,'(1x,a56)',advance='yes') "Enumerating and storing all determinants in the space..."
        call neci_flush(6)

        ! Determine the total number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)
        allocate(nI_list(nel, ndets))
        ! Now generate them again and store them this time.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        ! Count the number of determinants belonging to this processor.
        ndets_this_proc = 0
        do i = 1, ndets
            proc = DetermineDetNode(nel,nI_list(:,i), 0)
            if (proc == iProcIndex) ndets_this_proc = ndets_this_proc + 1
        end do

        allocate(ilut_list(0:NIfTot, ndets_this_proc))
        ilut_list = 0_n_int

        j = 0
        do i = 1, ndets
            proc = DetermineDetNode(nel,nI_list(:,i), 0)
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

        allocate(pert_ground_left(ndets_this_proc), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array to hold left perturbed ground state components.")
        end if

        allocate(sl_hamil(n_lanc_vecs_sl, n_lanc_vecs_sl))
        allocate(sl_h_eigv(n_lanc_vecs_sl))
        allocate(sl_overlaps(n_lanc_vecs_sl))
        allocate(trans_amps_left(n_lanc_vecs_sl))
        allocate(trans_amps_right(n_lanc_vecs_sl))
        sl_hamil = 0.0_dp
        sl_h_eigv = 0.0_dp
        sl_overlaps = 0.0_dp
        trans_amps_left = 0.0_dp
        trans_amps_right = 0.0_dp

        write(6,'(1x,a48)') "Allocating and calculating Hamiltonian matrix..."
        call neci_flush(6)
        call calculate_sparse_hamiltonian_parallel(ndets_sl, ilut_list, .true.)
        write(6,'(1x,a48,/)') "Hamiltonian allocation and calculation complete."
        call neci_flush(6)

        call return_perturbed_ground_spectral(left_perturb_spectral, ilut_list, pert_ground_left, left_pert_norm)
        call return_perturbed_ground_spectral(right_perturb_spectral, ilut_list, sl_vecs(:,1), right_pert_norm)
        
        ! Normalise the perturbed wave function, currently stored in sl_vecs(:,1),
        ! so that it can be used as the first lanczos vector.
        sl_vecs(:,1) = sl_vecs(:,1)*sqrt(pops_norm)/right_pert_norm

        deallocate(ilut_list)

    end subroutine init_spectral_lanczos

    subroutine return_perturbed_ground_spectral(perturbs, ilut_list, pert_ground_local, all_norm_pert)

        use bit_rep_data, only: NIfDBO, extract_sign
        use bit_reps, only: decode_bit_det
        use CalcData, only: pops_norm
        use DetBitOps, only: DetBitEq, EncodeBitDet, ilut_lt, ilut_gt
        use FciMCData, only: TotWalkers, CurrentDets
        use sort_mod, only: sort
        use SystemData, only: nel

        type(perturbation), intent(in), allocatable :: perturbs(:)
        integer(n_int), intent(in) :: ilut_list(:,:)
        real(dp), intent(inout) :: pert_ground_local(:)
        real(dp), intent(out) :: all_norm_pert

        integer :: i, j
        integer :: nI(nel)
        real(dp) :: real_sign(lenof_sign)
        real(dp) :: norm_pert
        character(len=*), parameter :: t_r = 'return_perturbed_ground_spectral'

        pert_ground_local = 0.0_dp

        ! Read in the POPSFILE and apply the perturbation operator in perturbs,
        ! as we do so. Afterwards, the perturbed wave function will be stored in
        ! CurrentDets.
        call read_popsfile_wrapper(perturbs)

        j = 0
        norm_pert = 0.0_dp
        call sort(CurrentDets(:,1:TotWalkers), ilut_lt, ilut_gt)
        do i = 1, int(TotWalkers, sizeof_int)
            call extract_sign(CurrentDets(:,i), real_sign)
            ! If the POPSFILE was generated from a calculation with the linear-scaling
            ! algorithm then it is possible to have the same determinant twice in a
            ! row, once with zero weight. If this is the case, skip the second case.
            if (i > 1) then
                if (DetBitEq(CurrentDets(:,i-1), CurrentDets(:,i), NIfDBO)) cycle
            end if
            norm_pert = norm_pert + real_sign(1)*real_sign(1)
            do
                j = j + 1
                if (DetBitEq(CurrentDets(:,i), ilut_list(:,j), NIfDBO)) then
                    pert_ground_local(j) = real_sign(1)
                    exit
                end if
            end do
        end do

        ! The contributions to pops_norm are added in as the determinants are
        ! read in from the popsfile, which is only done on the root process, so
        ! brooadcast the value to other processors.
        call MPIBCast(pops_norm)
        pert_ground_local = pert_ground_local/sqrt(pops_norm)

        call MPISumAll(norm_pert, all_norm_pert)
        all_norm_pert = sqrt(all_norm_pert)

    end subroutine return_perturbed_ground_spectral

    subroutine subspace_extraction_sl()

        integer :: lwork, info, i
        real(dp), allocatable :: work(:), left_pert_overlaps(:), left_pert_overlaps_all(:)

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

        ! The first Lanczos vector differs from the perturbed ground state only
        ! by a normalisation factor, which is accounted for here.
        trans_amps_right = sl_hamil(1,:)*right_pert_norm/sqrt(pops_norm)

        ! Calculate the overlaps between the left-perturbed ground state and all
        ! of the Lanczos vectors.
        allocate(left_pert_overlaps(n_lanc_vecs_sl))
        allocate(left_pert_overlaps_all(n_lanc_vecs_sl))
        do i = 1, n_lanc_vecs_sl
            left_pert_overlaps = matmul(pert_ground_left, sl_vecs)
        end do
        ! Sum the contributions to the overlaps from all processes.
        call MPISumAll(left_pert_overlaps, left_pert_overlaps_all)
        ! And then use these overlaps to calculate the overlaps of the
        ! left-perturbed ground state with the final eigenvectors.
        trans_amps_left = matmul(left_pert_overlaps_all, sl_hamil)
        deallocate(left_pert_overlaps)
        deallocate(left_pert_overlaps_all)

    end subroutine subspace_extraction_sl
    
    subroutine output_spectrum(neigv, eigv, spec_low, spec_high)

        use util_mod, only: get_free_unit

        integer, intent(in) :: neigv
        real(dp), intent(in) :: eigv(neigv)
        ! Data for the testsuite to use.
        real(dp), intent(out) :: spec_low, spec_high

        integer :: i, j, min_vec, temp_unit
        real(dp) :: omega, spectral_weight

        if (iProcIndex /= root) return

        temp_unit = get_free_unit()
        open(temp_unit, file='SPECTRAL_DATA',status='replace')
        write(temp_unit,'(1x,a53)') "Eigenvalues and left and right transition amplitudes:"
        do i = 1, neigv
            write(temp_unit,'(1x,i7,5x,f15.10,5x,f15.10,5x,f15.10)') &
                i, eigv(i), trans_amps_left(i), trans_amps_right(i)
        end do
        close(temp_unit)

        write(6,'(1x,a5,18X,a15)') "Omega", "Spectral_weight"

        ! Do we include the ground state in the spectrum or not?
        if (tIncludeGroundSpectral) then
            min_vec = 1
        else
            min_vec = 2
        end if

        omega = min_omega_spectral
        do i = 1, nomega_spectral + 1
            spectral_weight = 0.0_dp
            do j = min_vec, neigv
                spectral_weight = spectral_weight + &
                    (trans_amps_left(j)*trans_amps_right(j)*spectral_broadening)/&
                    (pi*(spectral_broadening**2 + (spectral_ground_energy-eigv(j)+omega)**2))
            end do
            write(6,'(f18.12, 4x, f18.12)') omega, spectral_weight
            omega = omega + delta_omega_spectral

            ! Store the values of the spectrum for the highest and lowest
            ! values of omega for the testsuite to use.
            if (i == 1) spec_low = spectral_weight
            if (i == nomega_spectral + 1) spec_high = spectral_weight
        end do

    end subroutine output_spectrum

    subroutine write_spec_lanc_testsuite_data(spec_low, spec_high)

        real(dp), intent(in) :: spec_low, spec_high

        write(6,'(/,1X,64("="))')
        write(6,'(1X,"Spectral Lanczos testsuite data:")')
        write(6,'(1X,"Lowest eigenvalue of H from the last Lanczos space:",2X,es20.13)') sl_h_eigv(1)
        write(6,'(1X,"Highest eigenvalue of H from the last Lanczos space:",1X,es20.13)') sl_h_eigv(n_lanc_vecs_sl)
        write(6,'(1X,"Spectral weight at the lowest omega value:",11X,es20.13)') spec_low
        write(6,'(1X,"Spectral weight at the highest omega value:",10X,es20.13)') spec_high
        write(6,'(1X,64("="))')

    end subroutine write_spec_lanc_testsuite_data

    subroutine end_spectral_lanczos()

        deallocate(sl_vecs)
        deallocate(full_vec_sl)
        deallocate(sl_overlaps)
        deallocate(trans_amps_left)
        deallocate(trans_amps_right)
        deallocate(pert_ground_left)
        deallocate(sl_hamil)
        deallocate(sl_h_eigv)
        deallocate(ndets_sl)
        deallocate(disps_sl)

        write(6,'(/,1x,a30,/)') "Exiting Spectral Lanczos code."
        call neci_flush(6)

    end subroutine end_spectral_lanczos

end module spectral_lanczos
