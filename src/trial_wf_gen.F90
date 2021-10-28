#include "macros.h"

module trial_wf_gen

    use bit_rep_data, only: NIfTot, nifd, flag_trial, flag_connected, IlutBits
    use CalcData
    use Parallel_neci
    use semi_stoch_gen
    use semi_stoch_procs
    use sparse_arrays
    use SystemData, only: nel, tHPHF, t_non_hermitian

    use guga_data, only: ExcitationInformation_t
    use guga_matrixElements, only: calc_guga_matrix_element
    use guga_bitrepops, only: write_det_guga, fill_csf_i

    use util_mod, only: get_free_unit, binary_search_custom, operator(.div.)
    use FciMCData, only: con_send_buf, NConEntry

    implicit none

contains

    subroutine init_trial_wf(trial_in, nexcit_calc, nexcit_keep, replica_pairs)

        use DetBitOps, only: ilut_lt, ilut_gt
        use enumerate_excitations, only: generate_connected_space
        use FciMCData, only: trial_space, trial_space_size, con_space, con_space_size, trial_wfs, tot_trial_space_size
        use FciMCData, only: trial_energies, ConTag, ConVecTag, TempTag, TrialTag, TrialWFTag
        use FciMCData, only: TrialTempTag, ConTempTag, OccTrialTag, Trial_Init_Time
        use FciMCData, only: OccConTag, CurrentTrialTag, current_trial_amps
        use FciMCData, only: MaxWalkersPart, tTrialHash, tIncCancelledInitEnergy
        use FciMCData, only: con_space_vecs, ntrial_excits, trial_numerator, trial_denom
        use FciMCData, only: tot_trial_numerator, tot_trial_denom, HashIndex
        use FciMCData, only: tot_init_trial_numerator, tot_init_trial_denom
        use FciMCData, only: init_trial_numerator, init_trial_denom
        use initial_trial_states, only: calc_trial_states_lanczos, calc_trial_states_qmc, calc_trial_states_direct
        use LoggingData, only: tWriteTrial, tCompareTrialAmps
        use MemoryManager, only: LogMemAlloc, LogMemDealloc
        use MPI_wrapper, only: root
        use ras_data, only: trial_ras
        use searching, only: remove_repeated_states
        use sort_mod, only: sort
        use SystemData, only: tAllSymSectors

        type(subspace_in) :: trial_in
        integer, intent(in) :: nexcit_calc, nexcit_keep
        logical, intent(in) :: replica_pairs

        integer :: i, ierr, num_states_on_proc, con_space_size_old
        integer :: excit, tot_con_space_size
        integer :: con_counts(0:nProcessors - 1)
        integer :: min_elem, max_elem, num_elem
        integer(MPIArg) :: trial_counts(0:nProcessors - 1), trial_displs(0:nProcessors - 1)
        integer(MPIArg) :: con_sendcounts(0:nProcessors - 1), con_recvcounts(0:nProcessors - 1)
        integer(MPIArg) :: con_senddispls(0:nProcessors - 1), con_recvdispls(0:nProcessors - 1)
        integer(n_int), allocatable, dimension(:, :) :: temp_space
        HElement_t(dp), allocatable :: trial_wfs_all_procs(:, :), temp_wfs(:, :)
        real(dp) :: temp_energies(nexcit_calc)
        character(len=*), parameter :: t_r = "init_trial_wf"

        ! Perform checks.
        if (tIncCancelledInitEnergy .and. (.not. tTrialHash)) &
            call stop_all(t_r, "The inc-cancelled-init-energy option cannot be used with the &
                               &trial-bin-search option.")
        if (.not. tUseRealCoeffs) call stop_all(t_r, "To use a trial wavefunction you must also &
            &use real coefficients.")
        if (nexcit_keep > nexcit_calc) call stop_all(t_r, "The number of required trial wave functions &
                               &is more than the number that has been requested to be calculated.")

        call set_timer(Trial_Init_Time)

        write(stdout, '(/,11("="),1X,"Trial wavefunction initialisation",1X,10("="))')

        ntrial_excits = nexcit_keep

        ! Simply allocate the trial vector to have up to 1 million elements for now...
        allocate(trial_space(0:NIfTot, 1000000), stat=ierr)
        call LogMemAlloc('trial_space', 1000000 * (NIfTot + 1), size_n_int, t_r, TrialTag, ierr)

        allocate(trial_energies(nexcit_keep))

        trial_energies = 0.0_dp
        trial_space = 0_n_int

        write(stdout, '("Generating the trial space...")'); call neci_flush(6)

        if (qmc_trial_wf) then
#ifdef CMPLX_
            call stop_all(t_r, "QMC trial state initiation not supported for complex wavefunctions.")
#else
            call calc_trial_states_qmc(trial_in, nexcit_keep, CurrentDets, HashIndex, replica_pairs, &
                                       trial_space_size, trial_space, trial_wfs, trial_counts, trial_displs)
#endif
        else
            if (allocated(trial_est_reorder)) then
                call calc_trial_states_lanczos(trial_in, nexcit_calc, trial_space_size, trial_space, temp_wfs, &
                                               temp_energies, trial_counts, trial_displs, trial_est_reorder)
            else
                call calc_trial_states_lanczos(trial_in, nexcit_calc, trial_space_size, trial_space, temp_wfs, &
                                               temp_energies, trial_counts, trial_displs)
            end if
        end if

        write(stdout, '("Size of trial space on this processor:",1X,i8)') trial_space_size; call neci_flush(6)

        if (.not. qmc_trial_wf) then
            ! Allocate the array to hold the final trial wave functions which we
            ! decide to keep, in the correct order.
            allocate(trial_wfs(nexcit_keep, trial_space_size), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error allocating trial_wfs.")
            ! Go through each replica and find which trial state matches it best.
            if (nexcit_calc > 1) then
                call assign_trial_states(replica_pairs, CurrentDets, HashIndex, trial_space, temp_wfs, &
                                         trial_wfs, temp_energies, trial_energies)
            else
                trial_wfs = temp_wfs
                trial_energies = temp_energies
                root_print "energy eigenvalue(s): ", trial_energies(1:nexcit_keep)
            end if
            deallocate(temp_wfs, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error deallocating temp_wfs.")
        end if

        ! At this point, each processor has only those states which reside on them, and
        ! have only counted those states. Send all states to all processors for the next bit.
        tot_trial_space_size = int(sum(trial_counts), sizeof_int)
        write(stdout, '("Total size of the trial space:",1X,i8)') tot_trial_space_size; call neci_flush(6)

        ! Use SpawnedParts as temporary space:
        call MPIAllGatherV(trial_space(0:NIfTot, 1:trial_space_size), &
                           SpawnedParts(0:NIfTot, 1:tot_trial_space_size), trial_counts, trial_displs)

        call sort(SpawnedParts(0:NIfTot, 1:tot_trial_space_size), ilut_lt, ilut_gt)

        call assign_elements_on_procs(tot_trial_space_size, min_elem, max_elem, num_elem)

        ! set the size of the entries in con_ht
#ifdef CMPLX_
        NConEntry = nifd + 2 * nexcit_keep
#else
        NConEntry = nifd + nexcit_keep
#endif

        if (num_elem > 0) then

            ! Find the states connected to the trial space. This typically takes a long time, so
            ! it is done in parallel by letting each processor find the states connected to a
            ! portion of the trial space.
            write(stdout, '("Calculating the number of states in the connected space...")'); call neci_flush(6)

            call generate_connected_space(num_elem, SpawnedParts(0:NIfTot, min_elem:max_elem), con_space_size)

            write(stdout, '("Attempting to allocate con_space. Size =",1X,F12.3,1X,"Mb")') &
                real(con_space_size, dp) * (NIfTot + 1.0_dp) * 7.629392e-06_dp; call neci_flush(6)
            allocate(con_space(0:NIfTot, con_space_size), stat=ierr)
            call LogMemAlloc('con_space', con_space_size * (NIfTot + 1), size_n_int, t_r, ConTag, ierr)
            con_space = 0_n_int

            write(stdout, '("States found on this processor, including repeats:",1X,i8)') con_space_size

            write(stdout, '("Generating and storing the connected space...")'); call neci_flush(6)

            call generate_connected_space(num_elem, SpawnedParts(0:NIfTot, min_elem:max_elem), &
                                          con_space_size, con_space)

            write(stdout, '("Removing repeated states and sorting by processor...")'); call neci_flush(6)

            call remove_repeated_states(con_space, con_space_size)

            call sort_space_by_proc(con_space(:, 1:con_space_size), con_space_size, con_sendcounts)

        else
            con_space_size = 0
            con_sendcounts = 0
            allocate(con_space(0, 0), stat=ierr)
            write(stdout, '("This processor will not search for connected states.")'); call neci_flush(6)
            !Although the size is zero, we should allocate it, because the rest of the code use it.
            !Otherwise, we get segmentation fault later.
            allocate(con_space(0:NIfTot, con_space_size), stat=ierr)
        end if

        write(stdout, '("Performing MPI communication of connected states...")'); call neci_flush(6)

        ! Send the connected states to their processors.
        ! con_sendcounts holds the number of states to send to other processors from this one.
        ! con_recvcounts will hold the number of states to be sent to this processor from the others.
        call MPIAlltoAll(con_sendcounts, 1, con_recvcounts, 1, ierr)
        con_space_size_old = con_space_size
        con_space_size = sum(con_recvcounts)
        ! The displacements necessary for mpi_alltoall.
        con_sendcounts = con_sendcounts * int(NIfTot + 1, MPIArg)
        con_recvcounts = con_recvcounts * int(NIfTot + 1, MPIArg)
        con_senddispls(0) = 0
        con_recvdispls(0) = 0
        do i = 1, nProcessors - 1
            con_senddispls(i) = con_senddispls(i - 1) + con_sendcounts(i - 1)
            con_recvdispls(i) = con_recvdispls(i - 1) + con_recvcounts(i - 1)
        end do

        write(stdout, '("Attempting to allocate temp_space. Size =",1X,F12.3,1X,"Mb")') &
            real(con_space_size, dp) * (NIfTot + 1.0_dp) * 7.629392e-06_dp; call neci_flush(6)
        allocate(temp_space(0:NIfTot, con_space_size), stat=ierr)
        call LogMemAlloc('temp_space', con_space_size * (NIfTot + 1), size_n_int, t_r, TempTag, ierr)

        call MPIAlltoAllV(con_space(:, 1:con_space_size_old), con_sendcounts, con_senddispls, &
                          temp_space(:, 1:con_space_size), con_recvcounts, con_recvdispls, ierr)

        if (allocated(con_space)) then
            deallocate(con_space, stat=ierr)
            call LogMemDealloc(t_r, ConTag, ierr)
        end if
        write(stdout, '("Attempting to allocate con_space. Size =",1X,F12.3,1X,"Mb")') &
            real(con_space_size, dp) * (NIfTot + 1.0_dp) * 7.629392e-06_dp; call neci_flush(6)
        allocate(con_space(0:NIfTot, 1:con_space_size), stat=ierr)
        call LogMemAlloc('con_space', con_space_size * (NIfTot + 1), size_n_int, t_r, ConTag, ierr)
        con_space = temp_space
        deallocate(temp_space, stat=ierr)
        call LogMemDealloc(t_r, TempTag, ierr)
        ! Finished sending to states to their processors.

        ! This will also sort the connected space.
        call remove_repeated_states(con_space, con_space_size)

        ! Remove states in the connected space which are also in the trial space.
        ! We don't use this anymore, but may want to try using it again at
        ! some point, so just leave it commented out.
!        call remove_list1_states_from_list2(SpawnedParts, con_space, tot_trial_space_size, con_space_size)

        call MPISumAll(con_space_size, tot_con_space_size)
        ! get the sizes of the connected/trial space on the other procs, for
        ! estimating the max size of the buffer
        call MPIAllGather(con_space_size, con_counts, ierr)
        ! allocate buffer for communication of con_ht
        ! it is normally also allocated upon initialization, so deallocate
        ! the dummy version
        if (allocated(con_send_buf)) deallocate(con_send_buf)
        ! the most we can communicate at once is the full size of the largest
        ! trial/connected space on a single proc
        allocate(con_send_buf(0:NConEntry, max(maxval(con_counts), maxval(trial_counts))))

        write(stdout, '("Total size of connected space:",1X,i10)') tot_con_space_size
        write(stdout, '("Size of connected space on this processor:",1X,i10)') con_space_size
        call neci_flush(6)

        ! Create the trial wavefunction from all processors, on all processors.
        allocate(trial_wfs_all_procs(nexcit_keep, tot_trial_space_size), stat=ierr)
        call MPIAllGatherV(trial_wfs, trial_wfs_all_procs, trial_counts, trial_displs)

        call sort_space_by_proc(SpawnedParts(0:NIfTot, 1:tot_trial_space_size), tot_trial_space_size, trial_counts)

        write(stdout, '("Generating the vector \sum_j H_{ij} \psi^T_j...")'); call neci_flush(6)
        allocate(con_space_vecs(nexcit_keep, con_space_size), stat=ierr)
        call LogMemAlloc('con_space_vecs', con_space_size, 8, t_r, ConVecTag, ierr)
        call generate_connected_space_vector(SpawnedParts, trial_wfs_all_procs, con_space, con_space_vecs)

        call MPIBarrier(ierr)

        if (tWriteTrial) call write_trial_space()
        if (tCompareTrialAmps) call update_compare_trial_file(.true.)

        allocate(current_trial_amps(nexcit_keep, MaxWalkersPart), stat=ierr)
        call LogMemAlloc('current_trial_amps', nexcit_keep * MaxWalkersPart, 8, t_r, CurrentTrialTag, ierr)
        call init_current_trial_amps()

        if (tTrialHash) call create_trial_hashtables(nexcit_keep)

        ! Set these to zero, to prevent junk being printed in the initial report.
        trial_numerator = 0.0_dp
        tot_trial_numerator = 0.0_dp
        trial_denom = 0.0_dp
        tot_trial_denom = 1.0_dp

        init_trial_numerator = 0.0_dp
        tot_init_trial_numerator = 0.0_dp
        init_trial_denom = 0.0_dp
        tot_init_trial_denom = 0.0_dp

        call halt_timer(Trial_Init_Time)

        if (.not. qmc_trial_wf) then
            root_print "Energy eigenvalue(s) of the trial space:", trial_energies(1:nexcit_keep)
        end if
        root_print "Trial wavefunction initialisation complete."
        root_print "Total time (seconds) taken for trial wavefunction initialisation:", &
            get_total_time(Trial_Init_Time)
        call neci_flush(6)

        if (tAS_TrialOffset) then
            call set_AS_TrialOffset(nexcit_keep, replica_pairs)
        end if

    end subroutine init_trial_wf

    subroutine assign_trial_states(replica_pairs, ilut_list, ilut_ht, trial_dets, trial_amps, &
                                   trials_kept, energies, energies_kept)

        ! Calculate the overlaps between each trial state and FCIQMC replica
        ! pair. For each replica, keep the trial state which has the largest
        ! overlap (by magnitude).

        use bit_reps, only: extract_sign
        use FciMCData, only: ll_node
        use hash, only: hash_table_lookup

        logical, intent(in) :: replica_pairs
        integer(n_int), intent(in) :: ilut_list(0:, :)
        type(ll_node), pointer, intent(inout) :: ilut_ht(:)
        integer(n_int), intent(in) :: trial_dets(0:, :)
        HElement_t(dp), intent(in) :: trial_amps(:, :)
        HElement_t(dp), intent(out) :: trials_kept(:, :)
        real(dp), intent(in) :: energies(:)
        real(dp), intent(out) :: energies_kept(:)

        integer :: i, idet, itrial, ireplica, det_ind, hash_val
        integer :: nI(nel), best_trial(1)
        real(dp) :: fciqmc_amps_real(size(energies_kept)), all_fciqmc_amps(lenof_sign)
        real(dp) :: overlaps_real(size(energies_kept), size(trial_amps, 1))
        real(dp) :: all_overlaps_real(size(energies_kept), size(trial_amps, 1))
#ifdef CMPLX_
        real(dp) :: overlaps_imag(size(energies_kept), size(trial_amps, 1))
        real(dp) :: all_overlaps_imag(size(energies_kept), size(trial_amps, 1))
        real(dp) :: fciqmc_amps_imag(size(energies_kept))
#endif
        logical :: tDetFound

        overlaps_real = 0.0_dp
        all_overlaps_real = 0.0_dp
#ifdef CMPLX_
        overlaps_imag = 0.0_dp
        all_overlaps_imag = 0.0_dp
        unused_var(replica_pairs)
#endif

        ! Loop over all basis states (determinants) in the trial space.
        ! For each, add the overlap for each trial amplitude to a running
        ! total for replica-trial state combinations.

        do idet = 1, size(trial_amps, 2)
            ! Find if this determinant is occupied in any of the FCIQMC wave
            ! functions.
            call decode_bit_det(nI, trial_dets(0:NIfTot, idet))
            ! Search the hash table for this determinant.
            call hash_table_lookup(nI, trial_dets(:, idet), nifd, ilut_ht, ilut_list, det_ind, hash_val, tDetFound)
            if (tDetFound) then
                call extract_sign(ilut_list(:, det_ind), all_fciqmc_amps)
#ifdef CMPLX_
                do i = 1, inum_runs
                    fciqmc_amps_real(i) = all_fciqmc_amps(min_part_type(i))
                    fciqmc_amps_imag(i) = all_fciqmc_amps(max_part_type(i))
                end do
#else
                if (replica_pairs) then
#if defined(PROG_NUMRUNS_) || defined(DOUBLERUN_)
                    do i = 1, lenof_sign.div.2
                        ! When using pairs of replicas, average their amplitudes.
                        fciqmc_amps_real(i) = sum(all_fciqmc_amps(2 * i - 1:2 * i)) / 2.0_dp
                    end do
#endif
                else
                    fciqmc_amps_real = all_fciqmc_amps
                end if
#endif
                ! Add in the outer product between fciqmc_amps and the trial
                ! state amplitudes.
                do itrial = 1, size(trial_amps, 1)
#ifdef CMPLX_
                    ! (a+ib)(c+id) = ac-bd +i(ad+bc)
                    overlaps_real(:, itrial) = overlaps_real(:, itrial) &
                                    + real(trial_amps(itrial, idet)) * fciqmc_amps_real - aimag(trial_amps(itrial, idet)) * fciqmc_amps_imag
                    overlaps_imag(:, itrial) = overlaps_imag(:, itrial) &
                                    + real(trial_amps(itrial, idet)) * fciqmc_amps_imag + aimag(trial_amps(itrial, idet)) * fciqmc_amps_real
#else
                    overlaps_real(:, itrial) = overlaps_real(:, itrial) + trial_amps(itrial, idet) * fciqmc_amps_real
#endif
                end do
            end if
        end do

        call MPISumAll(overlaps_real, all_overlaps_real)
#ifdef CMPLX_
        call MPISumAll(overlaps_imag, all_overlaps_imag)
#endif

        ! Now, find the best trial state for each FCIQMC replica:
        if (t_choose_trial_state) then

#ifdef CMPLX_
            do ireplica = 1, inum_runs
                trials_kept(ireplica, :) = trial_amps(trial_excit_choice(ireplica), :)
                energies_kept(ireplica) = energies(trial_excit_choice(ireplica))
            end do
#else
            if (replica_pairs) then
                do ireplica = 1, lenof_sign.div.2
                    trials_kept(ireplica, :) = trial_amps(trial_excit_choice(ireplica), :)
                    energies_kept(ireplica) = energies(trial_excit_choice(ireplica))

                    root_print "trial state: ", trial_excit_choice(ireplica), &
                        " chosen for replica: ", ireplica, &
                        " chosen by input, with energy: ", energies(trial_excit_choice(ireplica))
                end do
            else
                do ireplica = 1, lenof_sign
                    trials_kept(ireplica, :) = trial_amps(trial_excit_choice(ireplica), :)
                    energies_kept(ireplica) = energies(trial_excit_choice(ireplica))

                    root_print "trial state: ", trial_excit_choice(ireplica), &
                        " chosen for replica: ", ireplica, &
                        " chosen by input, with energy: ", energies(trial_excit_choice(ireplica))

                end do
            end if
#endif
        else
#ifdef CMPLX_
            do ireplica = 1, inum_runs
                best_trial = maxloc(abs(all_overlaps_real(ireplica, :)**2 + all_overlaps_imag(ireplica, :)**2))
                trials_kept(ireplica, :) = trial_amps(best_trial(1), :)
                energies_kept(ireplica) = energies(best_trial(1))
            end do
#else
            if (replica_pairs) then
                do ireplica = 1, lenof_sign.div.2
                    best_trial = maxloc(abs(all_overlaps_real(ireplica, :)))
                    trials_kept(ireplica, :) = trial_amps(best_trial(1), :)
                    energies_kept(ireplica) = energies(best_trial(1))
#ifdef DEBUG_
                    root_print "trial state: ", best_trial, " kept for replica ", ireplica, &
                        " based on overlap, with energy: ", energies(best_trial(1))
#endif
                end do
            else
                do ireplica = 1, lenof_sign
                    best_trial = maxloc(abs(all_overlaps_real(ireplica, :)))
                    trials_kept(ireplica, :) = trial_amps(best_trial(1), :)
                    energies_kept(ireplica) = energies(best_trial(1))

                    root_print "trial state: ", best_trial, " kept for replica ", ireplica, &
                        " based on overlap, with energy: ", energies(best_trial(1))

                end do
            end if
#endif
        end if

    end subroutine assign_trial_states

    subroutine remove_states_not_on_proc(ilut_list, ilut_list_size, update_trial_vector)

        use FciMCData, only: trial_wfs
        use load_balance_calcnodes, only: DetermineDetNode

        integer, intent(inout) :: ilut_list_size
        integer(n_int), intent(inout) :: ilut_list(0:, :)
        logical, intent(in) :: update_trial_vector
        integer :: i, counter, proc
        integer :: nI(nel)

        counter = 0
        do i = 1, ilut_list_size
            call decode_bit_det(nI, ilut_list(0:NIfTot, i))
            proc = DetermineDetNode(nel, nI, 0)
            ! If this state and the previous one were identical, don't add this state to the
            ! list so that repeats aren't included.
            if (proc == iProcIndex) then
                counter = counter + 1
                ilut_list(:, counter) = ilut_list(:, i)
                if (update_trial_vector) trial_wfs(:, counter) = trial_wfs(:, i)
            end if
        end do

        ilut_list_size = counter

    end subroutine remove_states_not_on_proc

    subroutine remove_list1_states_from_list2(list_1, list_2, list_1_size, list_2_size)

        use util_mod, only: binary_search

        integer, intent(in) :: list_1_size
        integer, intent(inout) :: list_2_size
        integer(n_int), intent(in) :: list_1(0:NIfTot, list_1_size)
        integer(n_int), intent(inout) :: list_2(0:NIfTot, list_2_size)
        integer :: i, counter, pos, min_ind

        min_ind = 1

        do i = 1, list_2_size
            ! Binary search list_1 to see if list_2(:,i) is in it.
            pos = binary_search(list_1(:, min_ind:list_1_size), list_2(:, i), NifD + 1)
            ! If it is in list 1, remove the state by setting it to 0.
            ! If it isn't in list 1 (pos < 0) then we can still search a smaller list next time.
            if (pos > 0) then
                list_2(:, i) = 0
                min_ind = min_ind + pos
            else
                min_ind = min_ind - pos - 1
            end if
            if (min_ind > list_1_size) exit
        end do

        ! Now compress the new list by overwriting the removed states:
        counter = 0
        do i = 1, list_2_size
            ! If the state wasn't set to 0:
            if (.not. all(list_2(:, i) == 0)) then
                counter = counter + 1
                list_2(:, counter) = list_2(:, i)
            end if
        end do

        list_2_size = counter

    end subroutine remove_list1_states_from_list2

    subroutine generate_connected_space_vector(trial_space, trial_vecs, con_space, con_vecs)

        ! Calculate the vector
        ! \sum_j H_{ij} \psi_j,
        ! where \psi is the trial vector, j runs over all trial space
        ! states and i runs over all connected space states. This is output
        ! in con_vecs.

        use hphf_integrals, only: hphf_off_diag_helement
        use MemoryManager, only: LogMemAlloc

        integer(n_int), intent(in) :: trial_space(0:, :)
        HElement_t(dp), intent(in) :: trial_vecs(:, :)
        integer(n_int), intent(in) :: con_space(0:, :)
        HElement_t(dp), intent(out) :: con_vecs(:, :)

        integer :: i, j, ierr
        integer :: nI(nel), nJ(nel)
        HElement_t(dp) :: H_ij
        character(len=*), parameter :: this_routine = "generate_connected_space_vector"
        type(ExcitationInformation_t) :: excitInfo
        type(CSF_Info_t) :: csf_i, csf_j

        con_vecs = 0.0_dp

        ! do i need to change this here for the non-hermitian transcorrelated
        ! hamiltonians?
        do i = 1, size(con_vecs, 2)
            call decode_bit_det(nI, con_space(0:NIfTot, i))

            ! i am only here in the guga case if i use the new way to calc
            ! the off-diagonal elements..
            if (tGUGA) csf_i = CSF_Info_t(con_space(0 : nifd, i))

            do j = 1, size(trial_vecs, 2)

                call decode_bit_det(nJ, trial_space(0:NIfTot, j))
                if (tGUGA) csf_j = CSF_Info_t(trial_space(0:NIfTot, j))

                if (all(con_space(0:nifd, i) == trial_space(0:nifd, j))) then
                    if (tHPHF) then
                        H_ij = hphf_diag_helement(nI, trial_space(:, j))
                    else if (tGUGA) then
                        H_ij = calcDiagMatEleGuga_nI(nI)
                    else
                        H_ij = get_helement(nI, nJ, 0)
                    end if
                else
                    ! need guga changes here!
                    ! and need
                    if (tHPHF) then
                        ! maybe i need a non-hermitian keyword here..
                        ! since I am not sure if this breaks the kneci
                        H_ij = hphf_off_diag_helement(nJ, nI, trial_space(:, j), con_space(:, i))
                        ! H_ij = hphf_off_diag_helement(nI, nJ, con_space(:,i), trial_space(:,j))
                    else if (tGUGA) then
                        ASSERT(.not. t_non_hermitian)
                        call calc_guga_matrix_element(&
                            con_space(:, i), csf_i, trial_space(:, j), csf_j, &
                            excitInfo, H_ij, .true.)
#ifdef CMPLX_
                        H_ij = conjg(H_ij)
#endif
!                         call calc_guga_matrix_element(con_space(:,i), trial_space(:,j), &
!                             excitInfo, H_ij, .true., 1)
                    else
                        ! maybe i need a non-hermitian keyword here..
                        ! since I am not sure if this breaks the kneci
                        H_ij = get_helement(nJ, nI, trial_space(:, j), con_space(:, i))
                        ! H_ij = get_helement(nI, nJ, con_space(:,i), trial_space(:,j))
                    end if
                end if
                ! workaround for complex matrix elements here:
#ifdef CMPLX_
                H_ij = conjg(H_ij)
#endif
                con_vecs(:, i) = con_vecs(:, i) + H_ij * trial_vecs(:, j)
            end do
        end do

    end subroutine generate_connected_space_vector

    subroutine assign_elements_on_procs(list_length, min_elem, max_elem, num_elem)

        ! Split list_length into nProcessor parts. Note that this is not done based on any hash.

        integer, intent(in) :: list_length
        integer, intent(out) :: min_elem, max_elem, num_elem
        integer :: floor_div_list_length, mod_list_length
        integer :: num_elem_all_procs(0:nProcessors - 1)
        integer :: i

        mod_list_length = mod(list_length, nProcessors)
        floor_div_list_length = (list_length - mod_list_length) / nProcessors

        do i = 0, nProcessors - 1
            num_elem_all_procs(i) = floor_div_list_length
            if (i < mod_list_length) num_elem_all_procs(i) = num_elem_all_procs(i) + 1
        end do

        num_elem = num_elem_all_procs(iProcIndex)

        if (num_elem == 0) then
            if (iProcIndex == 0) call stop_all("assign_elements_on_procs", "There are no states &
                                               &in the trial space.")
            min_elem = 0
            max_elem = 0
            return
        end if

        if (iProcIndex == 0) then
            min_elem = 1
            max_elem = num_elem
        else
            min_elem = 0
            do i = 0, iProcIndex - 1
                min_elem = min_elem + num_elem_all_procs(i)
            end do
            max_elem = min_elem + num_elem_all_procs(iProcIndex)
            min_elem = min_elem + 1
        end if

    end subroutine assign_elements_on_procs

    subroutine write_trial_space()

        use FciMCData, only: trial_space, trial_space_size

        integer :: i, j, k, iunit, ierr
        logical :: texist
        character(len=*), parameter :: t_r = 'write_trial_space'

        write(stdout, '("Writing the trial space to a file...")');
        iunit = get_free_unit()

        ! Let each processor write its trial states to the file. Each processor waits for
        ! the processor before it to finish before starting.
        do i = 0, nProcessors - 1

            if (iProcIndex == i) then

                if (i == 0) then
                    open(iunit, file='TRIALSPACE', status='replace')
                else
                    inquire (file='TRIALSPACE', exist=texist)
                    if (.not. texist) call stop_all(t_r, '"TRIALSPACE" file not found')
                    open(iunit, file='TRIALSPACE', status='old', position='append')
                end if

                do j = 1, trial_space_size
                    do k = 0, nifd
                        write(iunit, '(i24)', advance='no') trial_space(k, j)
                    end do
                    write(iunit, *)
                end do

                close(iunit)

            end if

            call MPIBarrier(ierr)

            if (i == 0) call MPIBCast(iunit, 0)

        end do

    end subroutine write_trial_space

    subroutine update_compare_trial_file(tFirstCall)

        ! Routine to output the trial wavefunction amplitudes and FCIQMC amplitudes in the trial
        ! space. This is a test routine and is very unoptimised.

        use bit_reps, only: extract_sign
        use FciMCData, only: trial_space, trial_space_size, trial_wfs
        use MPI_wrapper, only: root
        use searching, only: BinSearchParts

        logical, intent(in) :: tFirstCall
        logical :: tSuccess
        integer :: i, j, k, comp, MinInd, PartInd, n_walkers, iunit, ierr
        real(dp) :: temp_sign(lenof_sign)
        real(dp) :: all_norm_squared, norm_squared, norm

        MinInd = 1
        n_walkers = int(TotWalkers, sizeof_int)

        if (tFirstCall) then
            if (iProcIndex == root) then
                iunit = get_free_unit()
                open(iunit, file='TRIALCOMPARE', status='replace')
                close(iunit)
            end if
            call MPIBCast(iunit, root)
        else
            ! Calculate the norm of the wavefunction.
            norm_squared = 0.0_dp
            do i = 1, n_walkers
                call extract_sign(CurrentDets(:, i), temp_sign)
                norm_squared = norm_squared + sum(temp_sign**2)
            end do
            call MPIAllReduce(norm_squared, MPI_SUM, all_norm_squared)
            norm = sqrt(all_norm_squared)

            ! Let each processor write its amplitudes to the file. Each processor waits for
            ! the processor before it to finish before starting.
            do i = 0, nProcessors - 1

                if (iProcIndex == i) then

                    if (iProcIndex == 0) then
                        open(iunit, file='TRIALCOMPARE', status='old', position='rewind')
                        write(iunit, '(a14,5x,a6)') "Trial function", "FCIQMC"
                    else
                        open(iunit, file='TRIALCOMPARE', status='old', position='append')
                    end if

                    do j = 1, trial_space_size
                        write(iunit, '(es15.8,3x)', advance='no') trial_wfs(1, j)

                        call BinSearchParts(trial_space(:, j), MinInd, n_walkers, PartInd, tSuccess)

                        if (tSuccess) then
                            call extract_sign(CurrentDets(:, PartInd), temp_sign)
                            write(iunit, '(es15.8)', advance='yes') temp_sign / norm
                        else
                            write(iunit, '(es15.8)', advance='yes') 0.0_dp
                        end if

                        MinInd = PartInd
                    end do

                    close(iunit)

                end if

                call MPIBarrier(ierr)

            end do
        end if

    end subroutine update_compare_trial_file

    subroutine init_current_trial_amps()

        use bit_reps, only: set_flag, decode_bit_det
        use FciMCData, only: ll_node, trial_space, trial_space_size, con_space, con_space_size
        use FciMCData, only: con_space_vecs, current_trial_amps, HashIndex, trial_wfs, nWalkerHashes
        use FciMCData, only: CurrentDets
        use hash, only: FindWalkerHash
        use SystemData, only: nel

        integer :: i, hash_val
        integer :: nI(nel)
        type(ll_node), pointer :: temp_node

        ! Don't do anything is this is called before the trial wave function
        ! initialisation.
        if (.not. allocated(current_trial_amps)) return

        current_trial_amps = 0.0_dp

        do i = 1, trial_space_size
            call decode_bit_det(nI, trial_space(:, i))
            hash_val = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(hash_val)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if (all(trial_space(0:nifd, i) == CurrentDets(0:nifd, temp_node%ind))) then
                        call set_flag(CurrentDets(:, temp_node%ind), flag_trial)
                        current_trial_amps(:, temp_node%ind) = trial_wfs(:, i)
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            nullify (temp_node)
        end do

        do i = 1, con_space_size
            call decode_bit_det(nI, con_space(:, i))
            hash_val = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(hash_val)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if (all(con_space(0:nifd, i) == CurrentDets(0:nifd, temp_node%ind))) then
                        ! If not also in the trial space. If it is, then we
                        ! don't want the connected flag to be set, or the
                        ! connected vector amplitude to be used.
                        if (.not. test_flag(CurrentDets(:, temp_node%ind), flag_trial)) then
                            call set_flag(CurrentDets(:, temp_node%ind), flag_connected)
                            current_trial_amps(:, temp_node%ind) = con_space_vecs(:, i)
                        end if
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            nullify (temp_node)
        end do

    end subroutine init_current_trial_amps

    subroutine create_trial_hashtables(nexcit)

        use FciMCData, only: trial_space, trial_space_size, con_space, con_space_size
        use FciMCData, only: con_space_vecs, TrialTag, ConTag, TrialWFTag, ConVecTag
        use FciMCData, only: trial_wfs
        use hash, only: FindWalkerHash
        use MemoryManager, only: LogMemDealloc

        integer, intent(in) :: nexcit

        integer :: i, nclash, hash_val, mode, ierr
        integer :: nI(nel)
#ifdef CMPLX_
        integer(n_int) :: temp(2 * nexcit)
#else
        integer(n_int) :: temp(nexcit)
#endif
        character(len=*), parameter :: t_r = "create_trial_hashtables"

        ! Create the trial space hash table.

        allocate(trial_ht(trial_space_size), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error allocating trial_ht.")

        do i = 1, trial_space_size
            trial_ht(i)%nclash = 0
        end do

        ! When mode = 1, count the number of clashes.
        ! Alllocate arrays at the end of the mode = 1 loop.
        ! When mode = 2, fill in arrays.
        do mode = 1, 2
            do i = 1, trial_space_size
                call decode_bit_det(nI, trial_space(:, i))
                hash_val = FindWalkerHash(nI, trial_space_size)

                if (mode == 1) then
                    trial_ht(hash_val)%nclash = trial_ht(hash_val)%nclash + 1
                else
                    nclash = trial_ht(hash_val)%nclash + 1
                    trial_ht(hash_val)%nclash = nclash
                    trial_ht(hash_val)%states(0:nifd, nclash) = trial_space(0:nifd, i)
                    trial_ht(hash_val)%states(IlutBits%ind_pop:, nclash) = &
                        transfer(trial_wfs(:, i), temp)
                end if
            end do

            if (mode == 1) then
                do i = 1, size(trial_ht)
                    nclash = trial_ht(i)%nclash
                    allocate(trial_ht(i)%states(0:NConEntry, nclash))
                    ! Set this back to zero to use it as a counter next time
                    ! around (when mode == 2).
                    trial_ht(i)%nclash = 0
                end do
            end if
        end do

        ! No longer need these arrays in this form.
        if (allocated(trial_space)) then
            deallocate(trial_space, stat=ierr)
            call LogMemDealloc(t_r, TrialTag, ierr)
        end if
        if (allocated(trial_wfs)) then
            deallocate(trial_wfs, stat=ierr)
        end if

        ! Create the connected space hash table.

        allocate(con_ht(con_space_size), stat=ierr)
        if (ierr /= 0) then
            write(stdout, '("ierr:")') ierr
            call neci_flush(6)
            call stop_all("t_r", "Error in allocating con_ht array.")
        end if

        do i = 1, con_space_size
            con_ht(i)%nclash = 0
        end do

        do mode = 1, 2
            do i = 1, con_space_size
                call decode_bit_det(nI, con_space(:, i))
                hash_val = FindWalkerHash(nI, con_space_size)

                if (mode == 1) then
                    con_ht(hash_val)%nclash = con_ht(hash_val)%nclash + 1
                else
                    nclash = con_ht(hash_val)%nclash + 1
                    con_ht(hash_val)%nclash = nclash
                    con_ht(hash_val)%states(0:nifd, nclash) = con_space(0:nifd, i)
                    con_ht(hash_val)%states(IlutBits%ind_pop:, nclash) &
                        = transfer(con_space_vecs(:, i), temp)
                end if
            end do

            if (mode == 1) then
                do i = 1, size(con_ht)
                    nclash = con_ht(i)%nclash
                    allocate(con_ht(i)%states(0:NConEntry, nclash))
                    ! Set this back to zero to use it as a counter next time
                    ! around (when mode == 2).
                    con_ht(i)%nclash = 0
                end do
            end if
        end do

        ! No longer need these arrays in this form.
        if (allocated(con_space)) then
            deallocate(con_space, stat=ierr)
            call LogMemDealloc(t_r, ConTag, ierr)
        end if
        if (allocated(con_space_vecs)) then
            deallocate(con_space_vecs, stat=ierr)
            call LogMemDealloc(t_r, ConVecTag, ierr)
        end if

    end subroutine create_trial_hashtables

    subroutine end_trial_wf()

        use FciMCData, only: trial_space, con_space, con_space_vecs, TrialTag, ConTag
        use FciMCData, only: TrialWFTag, ConVecTag, CurrentTrialTag, ConTempTag, OccTrialTag
        use FciMCData, only: OccConTag, TrialTempTag, current_trial_amps
        use FciMCData, only: trial_wfs, trial_energies
        use MemoryManager, only: LogMemDealloc
        use sparse_arrays, only: deallocate_trial_hashtable

        character(len=*), parameter :: t_r = "end_trial_wf"
        integer :: ierr

        call deallocate_trial_hashtable(trial_ht)
        call deallocate_trial_hashtable(con_ht)

        if (allocated(trial_space)) then
            deallocate(trial_space, stat=ierr)
            call LogMemDealloc(t_r, TrialTag, ierr)
        end if
        if (allocated(trial_wfs)) then
            deallocate(trial_wfs, stat=ierr)
            call LogMemDealloc(t_r, TrialWFTag, ierr)
        end if
        if (allocated(trial_energies)) then
            deallocate(trial_energies, stat=ierr)
            if (ierr /= 0) write(stdout, '("Error when deallocating trial_energies:",1X,i8)') ierr
        end if
        if (allocated(con_space)) then
            deallocate(con_space, stat=ierr)
            call LogMemDealloc(t_r, ConTag, ierr)
        end if
        if (allocated(con_space_vecs)) then
            deallocate(con_space_vecs, stat=ierr)
            call LogMemDealloc(t_r, ConVecTag, ierr)
        end if
        if (allocated(current_trial_amps)) then
            deallocate(current_trial_amps, stat=ierr)
            call LogMemDealloc(t_r, CurrentTrialTag, ierr)
        end if
        if (allocated(con_send_buf)) deallocate(con_send_buf)

    end subroutine end_trial_wf

!------------------------------------------------------------------------------------------!

    subroutine refresh_trial_wf(trial_in, nexcit_calc, nexcit_keep, replica_pairs)
        implicit none
        type(subspace_in) :: trial_in
        integer, intent(in) :: nexcit_calc, nexcit_keep
        logical, intent(in) :: replica_pairs

        ! first, clear the current trial wavefunction
        call end_trial_wf()
        ! then, remove the flag from all determinants
        call reset_trial_space()
        ! now, generate the new trial wavefunction
        call init_trial_wf(trial_in, nexcit_calc, nexcit_keep, replica_pairs)
    end subroutine refresh_trial_wf

!------------------------------------------------------------------------------------------!

    subroutine reset_trial_space()
        use bit_reps, only: clr_flag
        implicit none
        integer(int64) :: i

        do i = 1_int64, TotWalkers
            ! remove the trial flag from all determinants
            call clr_flag(CurrentDets(:, i), flag_trial)
        end do
    end subroutine reset_trial_space

    !> Set the offset of the adaptive shift equal to the eigen energy(s)
    !> of the trial space.
    !> @param[in] nexcit_keep  number of wave functions/energies kept during trial-wf initialization
    !> @param[in] replica_pairs  whether replicas are assumed to be paired during trial-wf initialization
    subroutine Set_AS_TrialOffset(nexcit_keep, replica_pairs)
        use FciMCData, only: trial_energies
        integer, intent(in) :: nexcit_keep
        logical, intent(in) :: replica_pairs
        integer :: i
        character(*), parameter :: this_routine = 'Set_AS_TrialOffset'

        if (replica_pairs) then
            ASSERT(nexcit_keep == int(inum_runs / 2.0))
            do i = 1, nexcit_keep
                ShiftOffset(2 * i - 1) = trial_energies(i)
                ShiftOffset(2 * i) = trial_energies(i)
            end do
        else
            ASSERT(nexcit_keep == inum_runs)
            ShiftOffset(1:inum_runs) = trial_energies(1:inum_runs)
        end if

        tAS_Offset = .true.
        write(stdout, *) "The adaptive shift is offset by the eigen energy(s) of this trial-space."

    end subroutine

end module trial_wf_gen
