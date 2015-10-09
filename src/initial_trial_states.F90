#include "macros.h"

module initial_trial_states

    use bit_rep_data
    use constants
    use kp_fciqmc_data_mod

    implicit none

contains

    subroutine calc_trial_states_lanczos(space_in, nexcit, ndets_this_proc, trial_iluts, evecs_this_proc, evals, &
                                         space_sizes, space_displs, reorder)

        use bit_reps, only: decode_bit_det
        use CalcData, only: subspace_in
        use DetBitOps, only: ilut_lt, ilut_gt
        use FciMCData, only: ilutHF
        use lanczos_wrapper, only: frsblk_wrapper
        use Parallel_neci, only: MPIScatterV, MPIGatherV, MPIBCast, MPIArg, iProcIndex
        use Parallel_neci, only: nProcessors
        use ParallelHelper, only: root
        use semi_stoch_gen
        use sort_mod, only: sort
        use SystemData, only: nel, tAllSymSectors

        type(subspace_in) :: space_in
        integer, intent(in) :: nexcit
        integer, intent(out) :: ndets_this_proc
        integer(n_int), intent(out) :: trial_iluts(0:,:)
        real(dp), allocatable, intent(out) :: evecs_this_proc(:,:)
        real(dp), intent(out) :: evals(:)
        integer(MPIArg), intent(out) :: space_sizes(0:nProcessors-1), space_displs(0:nProcessors-1)
        integer, optional, intent(in) :: reorder(nexcit)

        integer(n_int), allocatable :: ilut_list(:,:)
        integer, allocatable :: det_list(:,:)
        integer :: i, j, max_elem_ind(1), ierr
        integer :: temp_reorder(nexcit)
        integer(MPIArg) :: ndets_all_procs, ndets_this_proc_mpi
        integer(MPIArg) :: sndcnts(0:nProcessors-1), displs(0:nProcessors-1)
        integer(MPIArg) :: rcvcnts
        integer, allocatable :: evec_abs(:)
        real(dp), allocatable :: evecs(:,:), evecs_transpose(:,:)
        character(len=*), parameter :: t_r = "calc_trial_states_lanczos"

        ndets_this_proc = 0
        trial_iluts = 0_n_int

        ! Choose the correct generating routine.
        if (space_in%tHF) call add_state_to_space(ilutHF, trial_iluts, ndets_this_proc)
        if (space_in%tPops) call generate_space_most_populated(space_in%npops, trial_iluts, ndets_this_proc)
        if (space_in%tRead) call generate_space_from_file(space_in%read_filename, trial_iluts, ndets_this_proc)
        if (space_in%tDoubles) call generate_sing_doub_determinants(trial_iluts, ndets_this_proc, .false.)
        if (space_in%tCAS) call generate_cas(space_in%occ_cas, space_in%virt_cas, trial_iluts, ndets_this_proc)
        if (space_in%tRAS) call generate_ras(space_in%ras, trial_iluts, ndets_this_proc)
        if (space_in%tOptimised) call generate_optimised_space(space_in%opt_data, space_in%tLimitSpace, &
                                                         trial_iluts, ndets_this_proc, space_in%max_size)
        if (space_in%tMP1) call generate_using_mp1_criterion(space_in%mp1_ndets, trial_iluts, ndets_this_proc)
        if (space_in%tFCI) then
            if (tAllSymSectors) then
                call gndts_all_sym_this_proc(trial_iluts, .true., ndets_this_proc)
            else
                call generate_fci_core(trial_iluts, ndets_this_proc)
            end if
        end if

        if (.not. (space_in%tPops .or. space_in%tRead .or. space_in%tDoubles .or. space_in%tCAS .or. &
                   space_in%tRAS .or. space_in%tOptimised .or. space_in%tMP1 .or. space_in%tFCI)) then
            call stop_all(t_r, "A space for the trial functions was not chosen.")
        end if

        ndets_this_proc_mpi = int(ndets_this_proc, MPIArg)
        call MPIAllGather(ndets_this_proc_mpi, space_sizes, ierr)
        ndets_all_procs = sum(space_sizes)

        if (ndets_all_procs < nexcit) call stop_all(t_r, "The number of excited states that you have asked &
            &for is larger than the size of the trial space used to create the excited states. Since this &
            &routine generates trial states that are orthogonal, this is not possible.")

        space_displs(0) = 0_MPIArg
        do i = 1, nProcessors-1
            space_displs(i) = sum(space_sizes(:i-1))
        end do

        call sort(trial_iluts(:,1:ndets_this_proc), ilut_lt, ilut_gt)

        if (iProcIndex == root) then
            allocate(ilut_list(0:NIfTot, ndets_all_procs))
        else
            ! On these other processes ilut_list and evecs_transpose are not
            ! needed, but we need them to be allocated for the MPI wrapper
            ! function to work, so just allocate them to be small.
            allocate(ilut_list(1,1))
            allocate(evecs_transpose(1,1))
        end if

        call MPIGatherV(trial_iluts(:,1:space_sizes(iProcIndex)), ilut_list, &
                        space_sizes, space_displs, ierr)

        ! Only perform the diagonalisation on the root process.
        if (iProcIndex == root) then
            allocate(det_list(nel, ndets_all_procs))

            do i = 1, ndets_all_procs
                call decode_bit_det(det_list(:,i), ilut_list(:,i))
            end do

            deallocate(ilut_list)

            allocate(evecs(ndets_all_procs, nexcit), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error allocating eigenvectors array.")
            evecs = 0.0_dp

            allocate(evec_abs(ndets_all_procs), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error allocating evec_abs array.")
            evec_abs = 0.0_dp
            
            ! Perform the Lanczos procedure.
            call frsblk_wrapper(det_list, int(ndets_all_procs, sizeof_int), nexcit, evals, evecs)

            ! For consistency between compilers, enforce a rule for the sign of
            ! the eigenvector. To do this, make sure that the largest component
            ! of each vector is positive. The largest component is found using
            ! maxloc. If there are multiple determinants with the same weight
            ! then the FORTRAN standard says that the first such component will
            ! be used, so this will hopefully be consistent across compilers.
            ! To avoid numerical errors bwteen compilers for such elements, we
            ! also round the array components to integers. Because evecs will
            ! be normalised, also multiply by 100000 before rounding.
            do i = 1, nexcit
                ! First find the maximum element.
                evec_abs = nint(100000.0_dp*abs(evecs(:,i)))
                max_elem_ind = maxloc(evec_abs)
                if (evecs(max_elem_ind(1),i) < 0.0_dp) then
                    evecs(:,i) = -evecs(:,i)
                end if
            end do

            deallocate(det_list)
            deallocate(evec_abs)

            ! Reorder the trial states, if the user has asked for this to be done.
            if (present(reorder)) then
                temp_reorder = reorder
                call sort(temp_reorder, evals)
                temp_reorder = reorder
                call sort(temp_reorder, evecs)
            end if

            ! Unfortunately to perform the MPIScatterV call we need the transpose
            ! of the eigenvector array.
            allocate(evecs_transpose(nexcit, ndets_all_procs), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error allocating transposed eigenvectors array.")
            evecs_transpose = transpose(evecs)
        else
            deallocate(ilut_list)
        end if

        call MPIBCast(evals, size(evals), root)

        ndets_this_proc_mpi = space_sizes(iProcIndex)
        ! The number of elements to send and receive in the MPI call, and the
        ! displacements.
        sndcnts = space_sizes*int(nexcit, MPIArg)
        rcvcnts = ndets_this_proc_mpi*int(nexcit, MPIArg)
        displs = space_displs*int(nexcit, MPIArg)

        ! Send the components to the correct processors using the following
        ! array as temporary space.
        allocate(evecs_this_proc(nexcit, ndets_this_proc), stat=ierr)
        call MPIScatterV(evecs_transpose, sndcnts, displs, evecs_this_proc, rcvcnts, ierr)
        if (ierr /= 0) call stop_all(t_r, "Error in MPIScatterV call.")

        ! Clean up.
        if (iProcIndex == root) deallocate(evecs)
        deallocate(evecs_transpose)

    end subroutine calc_trial_states_lanczos

    subroutine calc_trial_states_qmc(space_in, nexcit, qmc_iluts, qmc_ht, paired_replicas, ndets_this_proc, &
                                     trial_iluts, evecs_this_proc, space_sizes, space_displs)

        use CalcData, only: subspace_in
        use DetBitOps, only: ilut_lt, ilut_gt
        use FciMCData, only: ilutHF
        use Parallel_neci, only: nProcessors
        use semi_stoch_gen
        use sort_mod, only: sort
        use SystemData, only: tAllSymSectors

        type(subspace_in) :: space_in
        integer, intent(in) :: nexcit
        integer(n_int), intent(in) :: qmc_iluts(0:,:)
        type(ll_node), pointer, intent(inout) :: qmc_ht(:)
        logical, intent(in) :: paired_replicas
        integer, intent(out) :: ndets_this_proc
        integer(n_int), intent(out) :: trial_iluts(0:,:)
        real(dp), allocatable, intent(out) :: evecs_this_proc(:,:)
        integer(MPIArg), intent(out) :: space_sizes(0:nProcessors-1), space_displs(0:nProcessors-1)

        integer(n_int), allocatable :: ilut_list(:,:)
        integer, allocatable :: det_list(:,:)
        integer :: i, j, ierr
        integer(MPIArg) :: ndets_all_procs, ndets_this_proc_mpi
        character(len=*), parameter :: this_routine = "calc_trial_states_qmc"

        if (paired_replicas) then
            ASSERT(nexcit == lenof_sign/2)
        else
            ASSERT(nexcit == lenof_sign)
        end if

        ndets_this_proc = 0
        trial_iluts = 0_n_int

        ! Choose the correct generating routine.
        if (space_in%tHF) call add_state_to_space(ilutHF, trial_iluts, ndets_this_proc)
        if (space_in%tPops) call generate_space_most_populated(space_in%npops, trial_iluts, ndets_this_proc)
        if (space_in%tRead) call generate_space_from_file(space_in%read_filename, trial_iluts, ndets_this_proc)
        if (space_in%tDoubles) call generate_sing_doub_determinants(trial_iluts, ndets_this_proc, .false.)
        if (space_in%tCAS) call generate_cas(space_in%occ_cas, space_in%virt_cas, trial_iluts, ndets_this_proc)
        if (space_in%tRAS) call generate_ras(space_in%ras, trial_iluts, ndets_this_proc)
        if (space_in%tOptimised) call generate_optimised_space(space_in%opt_data, space_in%tLimitSpace, &
                                                         trial_iluts, ndets_this_proc, space_in%max_size)
        if (space_in%tMP1) call generate_using_mp1_criterion(space_in%mp1_ndets, trial_iluts, ndets_this_proc)
        if (space_in%tFCI) then
            if (tAllSymSectors) then
                call gndts_all_sym_this_proc(trial_iluts, .true., ndets_this_proc)
            else
                call generate_fci_core(trial_iluts, ndets_this_proc)
            end if
        end if

        if (.not. (space_in%tPops .or. space_in%tRead .or. space_in%tDoubles .or. space_in%tCAS .or. &
                   space_in%tRAS .or. space_in%tOptimised .or. space_in%tMP1 .or. space_in%tFCI)) then
            call stop_all(this_routine, "A space for the trial functions was not chosen.")
        end if

        ndets_this_proc_mpi = int(ndets_this_proc, MPIArg)
        call MPIAllGather(ndets_this_proc_mpi, space_sizes, ierr)
        ndets_all_procs = sum(space_sizes)

        space_displs(0) = 0_MPIArg
        do i = 1, nProcessors-1
            space_displs(i) = sum(space_sizes(:i-1))
        end do

        call sort(trial_iluts(:,1:ndets_this_proc), ilut_lt, ilut_gt)

        allocate(evecs_this_proc(nexcit, ndets_this_proc), stat=ierr)

        call get_qmc_trial_weights(trial_iluts, ndets_this_proc, qmc_iluts, qmc_ht, nexcit, paired_replicas, evecs_this_proc)

    end subroutine calc_trial_states_qmc

    subroutine get_qmc_trial_weights(trial_iluts, ntrial, qmc_iluts, qmc_ht, nexcit, paired_replicas, evecs_this_proc)

        use bit_reps, only: decode_bit_det
        use hash, only: hash_table_lookup
        use Parallel_neci, only: MPISumAll
        use SystemData, only: nel

        integer(n_int), intent(in) :: trial_iluts(0:,:)
        integer, intent(in) :: ntrial
        integer(n_int), intent(in) :: qmc_iluts(0:,:)
        type(ll_node), pointer, intent(inout) :: qmc_ht(:)
        integer, intent(in) :: nexcit
        logical, intent(in) :: paired_replicas
        real(dp), intent(out) :: evecs_this_proc(:,:)

        integer :: i, j, ind, hash_val
        integer :: nI(nel)
        real(dp) :: qmc_sign(lenof_sign), trial_sign(nexcit)
        real(dp) :: norm(nexcit), tot_norm(nexcit)
        logical :: found
        character(*), parameter :: this_routine = 'get_qmc_trial_weights'

        if (paired_replicas) then
            ASSERT(nexcit == lenof_sign/2)
        else
            ASSERT(nexcit == lenof_sign)
        end if

        norm = 0.0_dp

        ! Loop over all trial states.
        do i = 1, ntrial
            ! Check the QMC hash table to see if this state exists in the
            ! QMC list. 
            call decode_bit_det(nI, trial_iluts(:,i)) 
            call hash_table_lookup(nI, trial_iluts(:,i), NIfDBO, qmc_ht, &
                                   qmc_iluts, ind, hash_val, found)
            if (found) then
                call extract_sign(qmc_iluts(:,ind), qmc_sign)

                ! If using paired replicas then average the sign on each pair.
                if (paired_replicas) then
                    do j = 2, lenof_sign, 2
                        trial_sign(j/2) = sum(qmc_sign(j-1:j)/2.0_dp)
                    end do
                else
                    trial_sign = qmc_sign
                end if

                norm = norm + trial_sign**2
                evecs_this_proc(:,i) = trial_sign
            else
                evecs_this_proc(:,i) = 0.0_dp
            end if
        end do

        call MPISumAll(norm, tot_norm)

        do i = 1, ntrial
            evecs_this_proc(:,i) = evecs_this_proc(:,i)/sqrt(tot_norm)
        end do

    end subroutine get_qmc_trial_weights

    subroutine set_trial_populations(nexcit, ndets_this_proc, trial_vecs)

        use CalcData, only: InitialPart, InitWalkers, tStartSinglePart
        use Parallel_neci, only: MPISumAll

        integer, intent(in) :: nexcit, ndets_this_proc
        real(dp), intent(inout) :: trial_vecs(:,:)

        real(dp) :: eigenvec_pop, tot_eigenvec_pop
        integer :: i, j

        ! We need to normalise all of the vectors to have the correct number of
        ! walkers.
        do j = 1, nexcit
            eigenvec_pop = 0.0_dp
            do i = 1, ndets_this_proc
                eigenvec_pop = eigenvec_pop + abs(trial_vecs(j,i))
            end do

            call MPISumAll(eigenvec_pop, tot_eigenvec_pop)

            if (tStartSinglePart) then
                trial_vecs(j,:) = trial_vecs(j,:)*InitialPart/tot_eigenvec_pop
            else
                trial_vecs(j,:) = trial_vecs(j,:)*InitWalkers/tot_eigenvec_pop
            end if
        end do

    end subroutine set_trial_populations

    subroutine set_trial_states(ndets_this_proc, init_vecs, trial_iluts, &
                                semistoch_started, paired_replicas)

        use bit_reps, only: encode_sign
        use CalcData, only: tSemiStochastic, tTrialWavefunction
        use FciMCData, only: CurrentDets, TotWalkers, HashIndex, nWalkerHashes
        use replica_data, only: set_initial_global_data
        use hash, only: clear_hash_table, fill_in_hash_table
        use semi_stoch_procs, only: fill_in_diag_helements, copy_core_dets_to_spawnedparts
        use semi_stoch_procs, only: add_core_states_currentdet_hash, reinit_current_trial_amps

        integer, intent(in) :: ndets_this_proc
        real(dp), intent(in) :: init_vecs(:,:)
        integer(n_int), intent(in) :: trial_iluts(0:,:) 
        logical, intent(in) :: semistoch_started
        logical, intent(in), optional :: paired_replicas

        real(dp) :: real_sign(lenof_sign)
        integer :: i, j
        logical :: paired_reps_local

        if (present(paired_replicas)) then
            paired_reps_local = paired_replicas
        else
            paired_reps_local = .true.
        end if

        ! Now copy the amplitudes across to the CurrentDets array:
        ! First, get the correct states in CurrentDets.
        CurrentDets(0:NIfTot, 1:ndets_this_proc) = 0_n_int
        CurrentDets(0:NIfDBO, 1:ndets_this_proc) = trial_iluts(0:NIfDBO, 1:ndets_this_proc)

        ! Set signs.
        do i = 1, ndets_this_proc
            ! Construct the sign array to be encoded.
            if (paired_reps_local) then
                do j = 2, lenof_sign, 2
                    real_sign(j-1:j) = init_vecs(j/2,i)
                end do
            else
                do j = 1, lenof_sign
                    real_sign(j) = init_vecs(j, i)
                end do
            end if
            call encode_sign(CurrentDets(:,i), real_sign)
        end do

        TotWalkers = int(ndets_this_proc, int64)

        ! Reset and fill in the hash table.
        call clear_hash_table(HashIndex)
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets_this_proc, .true.)

        if (tSemiStochastic .and. semistoch_started) then
            ! core_space stores all core determinants from all processors. Move those on this
            ! processor to trial_iluts, which add_core_states_currentdet_hash uses.
            call copy_core_dets_to_spawnedparts()
            ! Any core space determinants which are not already in CurrentDets will be added
            ! by this routine.
            call add_core_states_currentdet_hash()
        end if

        if (tTrialWavefunction) call reinit_current_trial_amps()

        ! Calculate and store the diagonal elements of the Hamiltonian for
        ! determinants in CurrentDets.
        call fill_in_diag_helements()

        call set_initial_global_data(TotWalkers, CurrentDets)

    end subroutine set_trial_states

end module initial_trial_states
