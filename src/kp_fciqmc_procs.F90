#include "macros.h"

module kp_fciqmc_procs

    use bit_rep_data
    use bit_reps, only: decode_bit_det, encode_sign
    use constants
    use util_mod, only: near_zero
    use Parallel_neci, only: iProcIndex, MPISum, nProcessors
    use kp_fciqmc_data_mod
    use SystemData, only: tGUGA

    implicit none

contains

    subroutine store_krylov_vec(ivec)

        use FciMCData, only: TotWalkers, CurrentDets
        use global_det_data, only: det_diagH
        use hash, only: hash_table_lookup, add_hash_table_entry
        use Loggingdata, only: tPrintDataTables
        use semi_stoch_procs, only: check_determ_flag
        use SystemData, only: nel
        use util_mod, only: int_fmt

        integer, intent(in) :: ivec

        integer :: idet, iamp, sign_ind, flag_ind, hash_val, det_ind
        integer :: nI(nel)
        integer(n_int) :: temp, int_sign(lenof_sign_kp)
        logical :: tDetFound, tCoreDet
        real(dp) :: amp_fraction, real_sign(lenof_sign_kp)
        character(len=*), parameter :: t_r = "store_krylov_vec"

        if (tPrintDataTables) then
            write(stdout, '(a71)', advance='no') "# Adding the current walker configuration to the Krylov vector array..."
            call neci_flush(6)
        end if

        ! The index of the first element referring to the sign, for this ivec.
        sign_ind = nifd + lenof_sign_kp * (ivec - 1) + 1
        flag_ind = nifd + lenof_all_signs + 1

        ! Loop over all occupied determinants for this new Krylov vector.
        do idet = 1, int(TotWalkers)
            int_sign = CurrentDets(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign_kp - 1, idet)
            call extract_sign(CurrentDets(:, idet), real_sign)
            tCoreDet = check_determ_flag(CurrentDets(:, idet))
            ! Don't add unoccpied determinants, unless they are core determinants.
            if (IsUnoccDet(real_sign) .and. (.not. tCoreDet)) cycle

            tDetFound = .false.
            call decode_bit_det(nI, CurrentDets(:, idet))
            ! Search the hash table for this determinant.
            call hash_table_lookup(nI, CurrentDets(:, idet), nifd, krylov_vecs_ht, krylov_vecs, det_ind, hash_val, tDetFound)
            if (tDetFound) then
                ! In this case the determinant is already in the Krylov vector
                ! array.
                krylov_vecs(sign_ind:sign_ind + lenof_sign_kp - 1, det_ind) = int_sign
            else
                ! A new determiant needs to be added.
                TotWalkersKP = TotWalkersKP + 1
                det_ind = TotWalkersKP
                if (TotWalkersKP > krylov_vecs_length) then
                    call stop_all(t_r, "There are no slots left in the krylov_vecs array for the next determinant. &
                                       &You can increase the size of this array using the memory-factor option in &
                                       &the kp-fciqmc block of the input file.")
                end if

                ! Add the determinant's index to the hash table.
                call add_hash_table_entry(krylov_vecs_ht, det_ind, hash_val)

                ! Copy determinant data across.
                krylov_vecs(0:nifd, det_ind) = CurrentDets(0:nifd, idet)
                krylov_vecs(sign_ind:sign_ind + lenof_sign_kp - 1, det_ind) = int_sign
                krylov_helems(det_ind) = det_diagH(idet)
                krylov_vecs(flag_ind, det_ind) = CurrentDets(IlutBits%ind_flag, idet)
            end if

            ! Update information about how much of the hash table is filled.
            do iamp = 1, lenof_sign_kp
                if (.not. near_zero(real_sign(iamp))) then
                    nkrylov_amp_elems_used = nkrylov_amp_elems_used + 1
                end if
            end do

        end do

        if (tPrintDataTables) then
            write(stdout, '(1x,a5)', advance='yes') "Done."
            write(stdout, '(a56,'//int_fmt(TotWalkersKP, 1)//',1x,a17,'//int_fmt(krylov_vecs_length, 1)//')') &
                "# Number unique determinants in the Krylov vector array:", TotWalkersKP, "out of a possible", krylov_vecs_length
            amp_fraction = real(nkrylov_amp_elems_used, dp) / real(nkrylov_amp_elems_tot, dp)
            write(stdout, '(a69,1x,es10.4)') "# Fraction of the amplitude elements used in the Krylov vector array:", amp_fraction
            call neci_flush(6)
        end if

    end subroutine store_krylov_vec

    subroutine calc_overlap_matrix(nvecs, krylov_array, array_len, s_matrix)

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:, :)
        integer, intent(in) :: array_len
        real(dp), intent(out) :: s_matrix(:, :)

        integer :: idet, ivec, jvec, ind(nvecs)
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign_1(lenof_sign_kp), real_sign_2(lenof_sign_kp)

        ! Just in case!
        s_matrix(:, :) = 0.0_dp

        do jvec = 1, nvecs
            ! The first index of the sign in krylov_array, for each vector.
            ind(jvec) = nifd + lenof_sign_kp * (jvec - 1) + 1
        end do

        ! Loop over all determinants in the Krylov array.
        do idet = 1, array_len
            ! Loop over all Krylov vectors.
            do ivec = 1, nvecs
                int_sign = krylov_array(ind(ivec):ind(ivec) + lenof_sign_kp - 1, idet)
                real_sign_1 = transfer(int_sign, real_sign_1)
                if (IsUnoccDet(real_sign_1)) cycle

                do jvec = 1, ivec
                    int_sign = krylov_array(ind(jvec):ind(jvec) + lenof_sign_kp - 1, idet)
                    real_sign_2 = transfer(int_sign, real_sign_1)
                    if (IsUnoccDet(real_sign_2)) cycle

                    s_matrix(jvec, ivec) = s_matrix(jvec, ivec) + &
                                           (real_sign_1(kp_ind_1(1)) * real_sign_2(kp_ind_2(1)) + &
                                            real_sign_1(kp_ind_2(1)) * real_sign_2(kp_ind_1(1))) / 2.0_dp

                    ! Fill in the lower-half of the overlap matrix.
                    s_matrix(ivec, jvec) = s_matrix(jvec, ivec)
                end do
            end do
        end do

    end subroutine calc_overlap_matrix

    subroutine calc_perturbation_overlap(ivec)

        use hash, only: hash_table_lookup
        use SystemData, only: nel

        integer, intent(in) :: ivec

        integer :: idet, sign_ind, hash_val, det_ind
        integer :: nI(nel)
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign_1(lenof_sign_kp), real_sign_2(lenof_sign_kp)
        real(dp) :: overlap
        logical :: tDetFound

        overlap = 0.0_dp

        sign_ind = nifd + lenof_sign_kp * (ivec - 1) + 1

        ! Loop over all determinants in perturbed_ground
        do idet = 1, size(perturbed_ground, dim=2)
            call extract_sign(perturbed_ground(:, idet), real_sign_1)
            if (IsUnoccDet(real_sign_1)) cycle

            call decode_bit_det(nI, perturbed_ground(:, idet))
            ! Search to see if this determinant is in any Krylov vector.
            call hash_table_lookup(nI, perturbed_ground(:, idet), nifd, krylov_vecs_ht, krylov_vecs, &
                                   det_ind, hash_val, tDetFound)
            if (tDetFound) then
                ! If here then this determinant was found in the hash table.
                ! Add in the contributions to the overlap.
                int_sign = krylov_vecs(sign_ind:sign_ind + lenof_sign_kp - 1, det_ind)
                real_sign_2 = transfer(int_sign, real_sign_1)
                overlap = overlap + (real_sign_1(kp_ind_1(1)) * real_sign_2(kp_ind_2(1)) + &
                                     real_sign_1(kp_ind_2(1)) * real_sign_2(kp_ind_1(1))) / 2.0_dp
            end if
        end do

        pert_overlaps(ivec) = pert_overlaps(ivec) + overlap

    end subroutine calc_perturbation_overlap

    subroutine calc_hamil_exact(nvecs, krylov_array, array_len, h_matrix, h_diag)

        use DetBitOps, only: FindBitExcitLevel
        use Determinants, only: get_helement
        use FciMCData, only: Hii, exact_subspace_h_time
        use global_det_data, only: det_diagH
        use hphf_integrals, only: hphf_off_diag_helement
        use SystemData, only: tHPHF, nel
        use timing_neci, only: set_timer, halt_timer

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:, :)
        integer, intent(in) :: array_len
        real(dp), intent(out) :: h_matrix(:, :)
        real(dp), intent(in), optional :: h_diag(:)

        integer :: i, j, idet, jdet, ic
        integer(n_int) :: ilut_1(0:NIfTot), ilut_2(0:NIfTot)
        integer(n_int) :: int_sign(lenof_all_signs)
        integer :: nI(nel), nJ(nel)
        real(dp) :: real_sign_1(lenof_all_signs), real_sign_2(lenof_all_signs)
        real(dp) :: h_elem
        logical :: any_occ, occ_1, occ_2
        integer(4), allocatable :: occ_flags(:)

        call set_timer(exact_subspace_h_time)

        h_matrix = 0.0_dp

        allocate(occ_flags(array_len))
        occ_flags = 0

        ilut_1 = 0_n_int
        ilut_2 = 0_n_int

        ! Check to see if there are any replica 1 or 2 walkers on this determinant.
        do idet = 1, array_len
            int_sign = krylov_array(IlutBits%ind_pop:IlutBits%ind_pop + lenof_all_signs - 1, idet)

            any_occ = .false.
            do i = 1, nvecs
                any_occ = any_occ .or. (int_sign(kp_ind_1(i)) /= 0)
            end do
            if (any_occ) occ_flags = ibset(occ_flags(idet), 0)

            any_occ = .false.
            do i = 1, nvecs
                any_occ = any_occ .or. (int_sign(kp_ind_2(i)) /= 0)
            end do
            if (any_occ) occ_flags = ibset(occ_flags(idet), 1)
        end do

        ! Loop over all determinants in krylov_array.
        do idet = 1, array_len
            ilut_1(0:nifd) = krylov_array(0:nifd, idet)
            call decode_bit_det(nI, ilut_1)
            int_sign = krylov_array(IlutBits%ind_pop:IlutBits%ind_pop + lenof_all_signs - 1, idet)
            real_sign_1 = transfer(int_sign, real_sign_1)
            occ_1 = btest(occ_flags(idet), 0)
            occ_2 = btest(occ_flags(idet), 1)

            do jdet = idet, array_len
                ! If one of the two determinants are unoccupied in all vectors,
                ! then cycle.
                if (.not. ((occ_1 .and. btest(occ_flags(jdet), 1)) .or. &
                           (occ_2 .and. btest(occ_flags(jdet), 0)))) cycle

                ilut_2(0:nifd) = krylov_array(0:nifd, jdet)
                ic = FindBitExcitLevel(ilut_1, ilut_2)
                if (ic > 2) cycle

                call decode_bit_det(nJ, ilut_2)
                int_sign = krylov_array(IlutBits%ind_pop:IlutBits%ind_pop + lenof_all_signs - 1, jdet)
                real_sign_2 = transfer(int_sign, real_sign_1)

                if (idet == jdet) then
                    if (present(h_diag)) then
                        h_elem = h_diag(idet) + Hii
                    else
                        h_elem = det_diagH(idet) + Hii
                    end if
                else
                    if (tHPHF) then
                        h_elem = hphf_off_diag_helement(nI, nJ, ilut_1, ilut_2)
                    else
                        if (tGUGA) then
                            call stop_all("calc_hamil_exact", "modify for GUGA")
                        end if
                        h_elem = get_helement(nI, nJ, ic, ilut_1, ilut_2)
                    end if
                end if

                ! Finally, add in the contribution to all of the Hamiltonian elements.
                do i = 1, nvecs
                    do j = i, nvecs
                        if (idet == jdet) then
                            h_matrix(i, j) = h_matrix(i, j) + &
                                             h_elem * (real_sign_1(kp_ind_1(i)) * real_sign_2(kp_ind_2(j)) + &
                                                       real_sign_1(kp_ind_2(i)) * real_sign_2(kp_ind_1(j))) / 2
                        else
                            h_matrix(i, j) = h_matrix(i, j) + &
                                             h_elem * (real_sign_1(kp_ind_1(i)) * real_sign_2(kp_ind_2(j)) + &
                                                       real_sign_1(kp_ind_2(i)) * real_sign_2(kp_ind_1(j)) + &
                                                       real_sign_1(kp_ind_1(j)) * real_sign_2(kp_ind_2(i)) + &
                                                       real_sign_1(kp_ind_2(j)) * real_sign_2(kp_ind_1(i))) / 2
                        end if
                    end do
                end do

            end do
        end do

        do i = 1, nvecs
            do j = 1, i - 1
                h_matrix(i, j) = h_matrix(j, i)
            end do
        end do

        deallocate(occ_flags)

        call halt_timer(exact_subspace_h_time)

    end subroutine calc_hamil_exact

    subroutine communicate_kp_matrices(overlap_matrix, hamil_matrix, spin_matrix)

        ! Add all the overlap and projected Hamiltonian matrices together, with
        ! the result being held only on the root node.

        use MPI_wrapper, only: root

        real(dp), intent(inout) :: overlap_matrix(:, :)
        real(dp), intent(inout) :: hamil_matrix(:, :)
        real(dp), optional, intent(inout) :: spin_matrix(:, :)

        real(dp), allocatable :: mpi_mat_in(:, :), mpi_mat_out(:, :)
        integer :: nrow, ncol

        nrow = size(hamil_matrix, 1)

        if (present(spin_matrix)) then
            ncol = 3 * nrow
        else
            ncol = 2 * nrow
        end if

        allocate(mpi_mat_in(nrow, ncol))
        allocate(mpi_mat_out(nrow, ncol))

        mpi_mat_in(1:nrow, 1:nrow) = overlap_matrix
        mpi_mat_in(1:nrow, nrow + 1:2 * nrow) = hamil_matrix
        if (present(spin_matrix)) mpi_mat_in(1:nrow, 2 * nrow + 1:3 * nrow) = spin_matrix

        call MPISum(mpi_mat_in, mpi_mat_out)

        if (iProcIndex == root) then
            overlap_matrix = mpi_mat_out(1:nrow, 1:nrow)
            hamil_matrix = mpi_mat_out(1:nrow, nrow + 1:2 * nrow)
            if (present(spin_matrix)) spin_matrix = mpi_mat_out(1:nrow, 2 * nrow + 1:3 * nrow)
        end if

        deallocate(mpi_mat_in)
        deallocate(mpi_mat_out)

    end subroutine communicate_kp_matrices

    subroutine output_kp_matrices_wrapper(config_label, overlap_mats, hamil_mats)

        integer, intent(in) :: config_label
        real(dp), intent(in) :: overlap_mats(:, :, :)
        real(dp), intent(in) :: hamil_mats(:, :, :)

        call output_kp_matrices(config_label, 'overlap', overlap_mats)
        call output_kp_matrices(config_label, 'hamil  ', hamil_mats)

    end subroutine output_kp_matrices_wrapper

    subroutine output_kp_matrices(config_label, stem, matrices)

        use util_mod, only: int_fmt, get_free_unit

        integer, intent(in) :: config_label
        character(7), intent(in) :: stem
        real(dp), intent(in) :: matrices(:, :, :)
        character(25) :: ind1, filename
        integer :: i, j, k, temp_unit, nrepeats

        write(ind1, '(i15)') config_label

        filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
        nrepeats = size(matrices, 3)
        temp_unit = get_free_unit()

        open(temp_unit, file=trim(filename), status='replace')

        ! Write all the components of the various estimates of the matrix, above and including the
        ! diagonal, one after another on separate lines.
        do i = 1, size(matrices, 1)
            do j = i, size(matrices, 2)
                ! Write the index of the matrix element.
                write(temp_unit, '('//int_fmt(i, 0)//', "," ,'//int_fmt(j, 0)//')', advance='no') i, j
                do k = 1, nrepeats
                    write(temp_unit, '(1x,es19.12)', advance='no') matrices(i, j, k)
                end do
                write(temp_unit, '()', advance='yes')
            end do
        end do

        close(temp_unit)

    end subroutine output_kp_matrices

    subroutine average_kp_matrices_wrapper(config_label, nrepeats, overlap_mats, hamil_mats, &
                                           overlap_mean, hamil_mean, overlap_se, hamil_se)

        integer, intent(in) :: config_label, nrepeats
        real(dp), intent(in) :: overlap_mats(:, :, :), hamil_mats(:, :, :)
        real(dp), intent(out) :: overlap_mean(:, :), hamil_mean(:, :)
        real(dp), intent(out) :: overlap_se(:, :), hamil_se(:, :)

        call average_kp_matrices(nrepeats, hamil_mats, hamil_mean, hamil_se)
        call average_kp_matrices(nrepeats, overlap_mats, overlap_mean, overlap_se)
        if (tOutputAverageKPMatrices) then
            call output_average_kp_matrix(config_label, nrepeats, 'av_hamil  ', hamil_mean, hamil_se)
            call output_average_kp_matrix(config_label, nrepeats, 'av_overlap', overlap_mean, overlap_se)
        end if

    end subroutine average_kp_matrices_wrapper

    subroutine average_kp_matrices(nrepeats, matrices, mean, se)

        integer, intent(in) :: nrepeats
        real(dp), intent(in) :: matrices(:, :, :)
        real(dp), intent(out) :: mean(:, :), se(:, :)

        integer :: irepeat

        mean = 0.0_dp
        se = 0.0_dp

        do irepeat = 1, nrepeats
            mean = mean + matrices(:, :, irepeat)
        end do
        mean = mean / nrepeats

        if (nrepeats > 1) then
            do irepeat = 1, nrepeats
                se = se + (matrices(:, :, irepeat) - mean)**2
            end do
            se = se / ((nrepeats - 1) * nrepeats)
            se = sqrt(se)
        end if

    end subroutine average_kp_matrices

    subroutine output_average_kp_matrix(config_label, nrepeats, stem, mean, se)

        use util_mod, only: int_fmt, get_free_unit

        integer, intent(in) :: config_label, nrepeats
        character(10), intent(in) :: stem
        real(dp), intent(in) :: mean(:, :), se(:, :)

        integer :: irepeat
        character(25) :: ind1, filename
        integer :: i, j, k, temp_unit, repeat_ind

        write(ind1, '(i15)') config_label
        filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        ! Write all the components of the various estimates of the matrix, above and including the
        ! diagonal, one after another on separate lines.
        do i = 1, size(mean, 1)
            do j = i, size(mean, 2)
                ! Write the index of the matrix element.
                write(temp_unit, '('//int_fmt(i, 0)//', "1x" ,'//int_fmt(j, 0)//')', advance='no') i, j
                if (nrepeats > 1) then
                    ! Write the mean and standard error.
                    write(temp_unit, '(1x,es19.12,1x,a3,es19.12)') mean(i, j), "+/-", se(i, j)
                else if (nrepeats == 1) then
                    ! If we only have one sample then a standard error was not calculated, so
                    ! only output the mean.
                    write(temp_unit, '(1x,es19.12)') mean(i, j)
                end if
            end do
        end do

        close(temp_unit)

    end subroutine output_average_kp_matrix

    subroutine average_and_comm_pert_overlaps(nrepeats)

        ! Average and perform MPI communication.

        integer, intent(in) :: nrepeats

        pert_overlaps = pert_overlaps / nrepeats
        call MPISum(pert_overlaps, kp_all_pert_overlaps)

    end subroutine average_and_comm_pert_overlaps

    subroutine find_and_output_lowdin_eigv(config_label, nvecs, overlap_matrix, hamil_matrix, npositive, &
                                           all_evals, tOutput, spin_matrix, all_spin)

        use CalcData, only: tWritePopsNorm, pops_norm
        use util_mod, only: int_fmt, get_free_unit

        integer, intent(in) :: config_label, nvecs
        real(dp), intent(in) :: overlap_matrix(:, :), hamil_matrix(:, :)
        integer, intent(out) :: npositive
        real(dp), intent(out) :: all_evals(:, :)
        logical, intent(in) :: tOutput
        real(dp), optional, intent(in) :: spin_matrix(:, :)
        real(dp), optional, intent(out) :: all_spin(:, :)

        integer :: lwork, counter, i, nkeep, nkeep_len
        integer :: info, ierr, temp_unit
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: kp_final_hamil(:, :), kp_hamil_eigv(:), rotated_spin(:, :)
        real(dp) :: kp_pert_energy_overlaps(nvecs)
        character(7) :: string_fmt
        character(25) :: ind1, filename
        character(len=*), parameter :: stem = "lowdin"

        if (tOutput) then
            write(ind1, '(i15)') config_label
            filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
            temp_unit = get_free_unit()
            open(temp_unit, file=trim(filename), status='replace')

            if (tWritePopsNorm) then
                write(temp_unit, '(4("-"),a41,25("-"))') "Norm of unperturbed initial wave function"
                write(temp_unit, '(1x,es19.12,/)') sqrt(pops_norm)
            end if
        end if

        ! Create the workspace for the diagonaliser.
        lwork = max(1, 3 * nvecs - 1)
        allocate(work(lwork), stat=ierr)

        kp_overlap_eigenvecs = overlap_matrix

        ! Now perform the diagonalisation.
        call dsyev('V', 'U', nvecs, kp_overlap_eigenvecs, nvecs, kp_overlap_eigv, work, lwork, info)

        npositive = 0
        all_evals = 0.0_dp
        if (present(spin_matrix)) all_spin = 0.0_dp

        if (tOutput) write(temp_unit, '(4("-"),a26,40("-"))') "Overlap matrix eigenvalues"
        do i = 1, nvecs
            if (tOutput) write(temp_unit, '(1x,es19.12)') kp_overlap_eigv(i)
            if (kp_overlap_eigv(i) > 0.0_dp) npositive = npositive + 1
        end do

        if (tOutput .and. tOverlapPert) then
            write (temp_unit, '(/,"# The overlap of each final Hamiltonian eigenstate with each &
                              &requested perturbed ground state will be printed. The first &
                              &printed overlap is that with the first Krylov vector. The second &
                              &printed overlap is that with the vector specified with the &
                              &OVERLAP-PERTURB-ANNIHILATE and OVERLAP-PERTURB-CREATION &
                              &options.")')
        end if

        do nkeep = 1, npositive

            allocate(kp_final_hamil(nkeep, nkeep))
            allocate(kp_hamil_eigv(nkeep))
            if (present(spin_matrix)) allocate(rotated_spin(nkeep, nkeep))

            associate(transform_matrix => kp_transform_matrix(1:nvecs, 1:nkeep), &
                       inter_matrix => kp_inter_matrix(1:nvecs, 1:nkeep), &
                       eigenvecs_krylov => kp_eigenvecs_krylov(1:nvecs, 1:nkeep), &
                       init_overlaps => kp_init_overlaps(1:nkeep))

                counter = 0
                do i = nvecs - nkeep + 1, nvecs
                    counter = counter + 1
                    transform_matrix(:, counter) = kp_overlap_eigenvecs(:, i) / sqrt(kp_overlap_eigv(i))
                end do

                inter_matrix = matmul(hamil_matrix, transform_matrix)
                kp_final_hamil = matmul(transpose(transform_matrix), inter_matrix)

                call dsyev('V', 'U', nkeep, kp_final_hamil, nkeep, kp_hamil_eigv, work, lwork, info)

                eigenvecs_krylov = matmul(transform_matrix, kp_final_hamil)
                init_overlaps = matmul(overlap_matrix(1, :), eigenvecs_krylov) / scaling_factor

                if (tOverlapPert) kp_pert_energy_overlaps(1:nkeep) = matmul(kp_all_pert_overlaps, eigenvecs_krylov)

                ! The spin matrix in the final eigenvalue basis.
                if (present(spin_matrix)) then
                    ! Use inter_matrix as temporary space.
                    inter_matrix = matmul(spin_matrix, eigenvecs_krylov)
                    rotated_spin = matmul(transpose(eigenvecs_krylov), inter_matrix)
                end if

                if (tOutput) then
                    nkeep_len = ceiling(log10(real(abs(nkeep) + 1, dp)))
                    write(string_fmt, '(i2,a5)') 15 - nkeep_len, '("-")'
                    write(temp_unit, '(/,4("-"),a37,'//int_fmt(nkeep, 1)//',1x,a12,'//string_fmt//')') &
                        "Eigenvalues and overlaps when keeping", nkeep, "eigenvectors"
                    do i = 1, nkeep
                        write(temp_unit, '(1x,es19.12,1x,es19.12)', advance='no') kp_hamil_eigv(i), init_overlaps(i)
                        if (tOverlapPert) write(temp_unit, '(1x,es19.12)', advance='no') kp_pert_energy_overlaps(i)
                        write(temp_unit, '()', advance='yes')
                    end do
                end if

                all_evals(1:nkeep, nkeep) = kp_hamil_eigv(1:nkeep)
                if (present(spin_matrix)) then
                    do i = 1, nkeep
                        all_spin(i, nkeep) = rotated_spin(i, i)
                    end do
                end if

            end associate

            deallocate(kp_final_hamil)
            deallocate(kp_hamil_eigv)
            if (allocated(rotated_spin)) deallocate(rotated_spin)

        end do

        deallocate(work)

        if (tOutput) close(temp_unit)

    end subroutine find_and_output_lowdin_eigv

    subroutine find_and_output_gs_eigv(config_label, nvecs)

        ! Use the Gram-Schmidt approach rather than the Lowdin approach above
        ! (see Phys. Rev. B. 85, 205119).

        use util_mod, only: int_fmt, get_free_unit

        integer, intent(in) :: config_label, nvecs

        integer :: lwork, counter, nkeep, nkeep_len, temp_unit
        integer :: npositive, info, ierr
        integer :: i, j, n, m
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: kp_final_hamil(:, :), kp_hamil_eigv(:)
        character(7) :: string_fmt
        character(25) :: ind1, filename
        character(len=*), parameter :: stem = "gram_schmidt"

        write(ind1, '(i15)') config_label
        filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        ! Create the workspace for the diagonaliser.
        lwork = max(1, 3 * nvecs - 1)
        allocate(work(lwork), stat=ierr)

        ! Use the following allocated arrays as work space for the following routine. Not ideal, I know...
        allocate(kp_final_hamil(nvecs, nvecs))
        associate(S_tilde => kp_inter_matrix, N => kp_init_overlaps)
            call construct_gs_transform_matrix(kp_overlap_mean, kp_transform_matrix, S_tilde, kp_final_hamil, N, nvecs)
            npositive = 0
            do i = 1, nvecs
                if (N(i) > 0.0_dp) then
                    npositive = npositive + 1
                else
                    exit
                end if
            end do

            write(temp_unit, '(4("-"),a21,45("-"))') "Normalisation factors"
            do i = 1, npositive
                write(temp_unit, '(1x,es19.12)') N(i)
            end do
        end associate
        deallocate(kp_final_hamil)

        do nkeep = 1, npositive

            allocate(kp_final_hamil(nkeep, nkeep))
            allocate(kp_hamil_eigv(nkeep))

            associate(S => kp_transform_matrix(1:nvecs, 1:nkeep), &
                       init_overlaps => kp_init_overlaps(1:nkeep))

                kp_final_hamil = 0.0_dp
                do n = 1, nkeep
                    do m = 1, nkeep
                        do i = 1, n
                            do j = 1, m
                                kp_final_hamil(n, m) = kp_final_hamil(n, m) + S(i, n) * kp_hamil_mean(i, j) * S(j, m)
                            end do
                        end do
                    end do
                end do

                call dsyev('V', 'U', nkeep, kp_final_hamil, nkeep, kp_hamil_eigv, work, lwork, info)

                init_overlaps = kp_final_hamil(:, 1) * sqrt(kp_overlap_mean(1, 1))

                nkeep_len = ceiling(log10(real(abs(nkeep) + 1, dp)))
                write(string_fmt, '(i2,a5)') 15 - nkeep_len, '("-")'
                write(temp_unit, '(/,4("-"),a37,'//int_fmt(nkeep, 1)//',1x,a12,'//string_fmt//')') &
                    "Eigenvalues and overlaps when keeping", nkeep, "eigenvectors"
                do i = 1, nkeep
                    write(temp_unit, '(1x,es19.12,1x,es19.12)') kp_hamil_eigv(i), init_overlaps(i)
                end do

            end associate

            deallocate(kp_final_hamil)
            deallocate(kp_hamil_eigv)

        end do

        deallocate(work)

        close(temp_unit)

    end subroutine find_and_output_gs_eigv

    subroutine construct_gs_transform_matrix(overlap, S, S_tilde, k, N, matrix_size)

        ! Construct the matrix, S, which transforms the Krylov vectors to a set of orthogonal vectors.

        ! We use the notation from the Appendix of Phys. Rev. B. 85, 205119 for the matrices and the
        ! indices, except we use 'm' instead of 'n' and 'l' instead of 'k' for the indices.

        ! Usage: overlap on input should contain the overlap matrix of the Krylov vectors that
        ! you want to produce a transformation matrix for. All other matrices input should be
        ! allocated on input, and should be the same size as overlap.

        real(dp), intent(inout) :: overlap(:, :), S(:, :), S_tilde(:, :), k(:, :), N(:)
        integer, intent(in) :: matrix_size
        integer :: m, i, l, p, q

        S = 0.0_dp
        S_tilde = 0.0_dp
        k = 0.0_dp
        N = 0.0_dp

        do m = 1, matrix_size
            S(m, m) = 1.0_dp
            S_tilde(m, m) = 1.0_dp
        end do

        N(1) = overlap(1, 1)
        S(1, 1) = 1 / sqrt(N(1))

        do m = 2, matrix_size
            do i = 1, m - 1
                do l = 1, i
                    k(i, m) = k(i, m) + S(l, i) * overlap(m, l)
                end do
            end do
            do p = 1, m - 1
                do q = p, m - 1
                    S_tilde(p, m) = S_tilde(p, m) - k(q, m) * S(p, q)
                end do
            end do
            do q = 1, m
                do p = 1, m
                    N(m) = N(m) + S_tilde(p, m) * S_tilde(q, m) * overlap(p, q)
                end do
            end do
            if (N(m) < 0.0_dp) return
            S(:, m) = S_tilde(:, m) / sqrt(N(m))
        end do

    end subroutine construct_gs_transform_matrix

    subroutine print_populations_kp()

        ! A useful test routine which will output the total walker population on both
        ! replicas, for each Krylov vector.

        integer :: ihash
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: real_sign(lenof_all_signs), total_pop(lenof_all_signs)
        type(ll_node), pointer :: temp_node

        int_sign = 0_n_int
        total_pop = 0.0_dp
        real_sign = 0.0_dp

        do ihash = 1, nhashes_kp
            temp_node => krylov_vecs_ht(ihash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    int_sign = krylov_vecs(IlutBits%ind_pop:IlutBits%ind_pop + lenof_all_signs - 1, temp_node%ind)
                    real_sign = transfer(int_sign, real_sign)
                    total_pop = total_pop + abs(real_sign)
                    temp_node => temp_node%next
                end do
            end if
        end do

        nullify (temp_node)

        write(stdout, *) "krylov_vecs populations:", total_pop

    end subroutine print_populations_kp

    subroutine print_amplitudes_kp(irepeat)

        ! A (*very* slow and memory intensive) test routine to print the current amplitudes (as stored
        ! in CurrentDets) of *all* determinants to a file. The amplitude of each replica will be printed
        ! one after the other. Since this is intended to be used with kp-fciqmc, irepeat is the number of
        ! the current repeat, but it will simply be used in naming the output file.

        use DetBitOps, only: EncodeBitDet, IsAllowedHPHF
        use FciMCData, only: nWalkerHashes, HashIndex, CurrentDets, HFSym
        use hash, only: FindWalkerHash
        use gndts_mod, only: gndts
        use sym_mod, only: getsym
        use SystemData, only: nel, nbasis, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict, BasisFn
        use util_mod, only: int_fmt, get_free_unit

        integer, intent(in) :: irepeat
        integer, allocatable :: nI_list(:, :)
        integer :: temp(1, 1), hf_ind, ndets
        integer :: i, j, counter, temp_unit, DetHash
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign(lenof_sign)
        type(ll_node), pointer :: temp_node
        type(BasisFn) :: iSym
        character(15) :: ind, filename

        ! Determine the total number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)
        allocate(nI_list(nel, ndets))
        ! Generate the determinants and move them to nI_list.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        write(ind, '(i15)') irepeat
        filename = trim('amps.'//adjustl(ind))

        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        counter = 0

        do i = 1, ndets
            call getsym(nI_list(:, i), nel, G1, nBasisMax, iSym)
            ! Only carry on if the symmetry of this determinant is correct.
            if (iSym%Sym%S /= HFSym%Sym%S .or. iSym%Ms /= HFSym%Ms .or. iSym%Ml /= HFSym%Ml) cycle
            call EncodeBitDet(nI_list(:, i), ilut)
            if (.not. IsAllowedHPHF(ilut(0:nifd))) cycle
            counter = counter + 1
            real_sign = 0.0_dp
            DetHash = FindWalkerHash(nI_list(:, i), nWalkerHashes)
            temp_node => HashIndex(DetHash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if (all(ilut(0:nifd) == CurrentDets(0:nifd, temp_node%ind))) then
                        int_sign = CurrentDets(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign - 1, temp_node%ind)
                        real_sign = transfer(int_sign, real_sign)
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            do j = 1, lenof_sign
                write(temp_unit, '(a1,'//int_fmt(counter, 0)//',a1,'//int_fmt(j, 0)//',a1,1x,es19.12)') &
                    "(", counter, ",", j, ")", real_sign(j)
            end do
        end do

        close(temp_unit)

        deallocate(nI_list)

    end subroutine print_amplitudes_kp

    subroutine write_ex_state_header(nvecs, irepeat)

        use util_mod, only: int_fmt, get_free_unit

        integer, intent(in) :: nvecs, irepeat
        integer :: ivec, icolumn
        integer :: temp_unit
        character(22) :: column_label
        character(len=*), parameter :: filename = "EIGV_DATA"

        temp_unit = get_free_unit()

        if (irepeat == 1) then
            open(temp_unit, file=trim(filename), status='replace')
        else
            open(temp_unit, file=trim(filename), status='old', position='append')
        end if

        ! Write header.
        if (irepeat == 1) then
            ! The number of the column.
            icolumn = 1

            write(temp_unit, '("#",1X,"1. Iteration")', advance='no')

            ! Energy estimates.
            do ivec = 1, nvecs
                icolumn = icolumn + 1
                write(column_label, '('//int_fmt(icolumn, 0)//',".",1X,"Energy",1X,'//int_fmt(ivec, 0)//')') icolumn, ivec
                column_label = adjustr(column_label)
                write(temp_unit, '(a22)', advance='no') column_label
            end do
            do ivec = 1, nvecs
                icolumn = icolumn + 1
                write(column_label, '('//int_fmt(icolumn, 0)//',".",1X,"Diag. energy",1X,'//int_fmt(ivec, 0)//')') icolumn, ivec
                column_label = adjustr(column_label)
                write(temp_unit, '(a22)', advance='no') column_label
            end do

            ! Spin estimates.
            if (tCalcSpin) then
                do ivec = 1, nvecs
                    icolumn = icolumn + 1
                    write(column_label, '('//int_fmt(icolumn, 0)//',".",1X,"Spin^2",1X,'//int_fmt(ivec, 0)//')') icolumn, ivec
                    column_label = adjustr(column_label)
                    write(temp_unit, '(a22)', advance='no') column_label
                end do
                do ivec = 1, nvecs
                    icolumn = icolumn + 1
                    write(column_label, '('//int_fmt(icolumn, 0)//',".",1X,"Diag spin^2",1X,'//int_fmt(ivec, 0)//')') icolumn, ivec
                    column_label = adjustr(column_label)
                    write(temp_unit, '(a22)', advance='no') column_label
                end do
            end if

            write(temp_unit, '()')
        end if

        write(temp_unit, '("#",1X,"Repeat",'//int_fmt(irepeat, 1)//')') irepeat

        close(temp_unit)

    end subroutine write_ex_state_header

    subroutine write_ex_state_data(niters, nlowdin, lowdin_evals, hamil_matrix, overlap_matrix, spin_matrix, lowdin_spin)

        use SystemData, only: nel
        use util_mod, only: get_free_unit

        integer, intent(in) :: niters
        integer, intent(in) :: nlowdin
        real(dp), intent(in) :: lowdin_evals(:, :)
        real(dp), intent(in) :: hamil_matrix(:, :), overlap_matrix(:, :)
        real(dp), optional, intent(in) :: spin_matrix(:, :), lowdin_spin(:, :)

        integer :: ivec, nvecs
        integer :: temp_unit
        character(len=*), parameter :: filename = "EIGV_DATA"

        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='old', position='append')

        nvecs = size(lowdin_evals, 1)

        write(temp_unit, '(5X,i9)', advance='no') niters

        do ivec = 1, nlowdin
            write(temp_unit, '(3X,es19.12)', advance='no') lowdin_evals(ivec, nlowdin)
        end do
        do ivec = nlowdin + 1, nvecs
            write(temp_unit, '(12X,"NaN",7X)', advance='no')
        end do

        ! Diagonal energies.
        do ivec = 1, nvecs
            write(temp_unit, '(3X,es19.12)', advance='no') hamil_matrix(ivec, ivec) / overlap_matrix(ivec, ivec)
        end do

        if (present(lowdin_spin)) then
            do ivec = 1, nlowdin
                write(temp_unit, '(3X,es19.12)', advance='no') lowdin_spin(ivec, nlowdin) + (0.75_dp * nel)
            end do
            do ivec = nlowdin + 1, nvecs
                write(temp_unit, '(12X,"NaN",7X)', advance='no')
            end do
        end if

        ! Diagonal spin entries.
        if (present(spin_matrix)) then
            do ivec = 1, nvecs
                write(temp_unit, '(3X,es19.12)', advance='no') spin_matrix(ivec, ivec) / overlap_matrix(ivec, ivec) + (0.75_dp * nel)
            end do
        end if

        close(temp_unit)

    end subroutine write_ex_state_data

    subroutine write_kpfciqmc_testsuite_data(s_sum, h_sum)

        ! Write out information about the completed KP-FCIQMC calculation, for
        ! use in the testsuite.

        real(dp), intent(in) :: s_sum, h_sum

        write(stdout, '(/,1X,64("="))')
        write(stdout, '(1X,"KP-FCIQMC testsuite data:")')
        write(stdout, '(1X,"Sum of overlap matrix elements:",12X,es20.13)') s_sum
        write(stdout, '(1X,"Sum of H elements:",25X,es20.13)') h_sum
        write(stdout, '(1X,64("="))')

    end subroutine write_kpfciqmc_testsuite_data

end module kp_fciqmc_procs
