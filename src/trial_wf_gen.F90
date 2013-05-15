module trial_wf_gen

    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: encode_det
    use CalcData, only: tSortDetermToTop
    use davidson, only: perform_davidson, davidson_eigenvalue, davidson_eigenvector, &
                        sparse_hamil_type
    use DetBitOps, only: FindBitExcitLevel, ilut_lt, ilut_gt
    use DeterminantData, only: write_det
    use FciMCData, only: trial_space, trial_space_size, con_space, &
                         con_space_size, trial_wf, trial_energy, &
                         con_space_vector, ilutHF, Hii, nSingles, nDoubles, &
                         ConTag, DavidsonTag, TempTag, TrialTag, TrialWFTag, iHFProc
    use hphf_integrals, only: hphf_off_diag_helement
    use LoggingData, only: tWriteTrial
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier, MPIAllToAll, MPISumAll, &
                             MPIAllToAllV, MPIArg, MPIBCast
    use ParallelHelper, only: root
    use semi_stoch_gen
    use semi_stoch_procs
    use sparse_hamil
    use SystemData, only: nel, tDoublesTrial, tOptimisedTrial, tCASTrial, tPopsTrial, &
                          tLowETrial, tHPHF
    use util_mod, only: get_free_unit

    implicit none

contains

    subroutine init_trial_wf()

        integer :: i, ierr, num_states_on_proc
        integer :: excit, tot_trial_space_size, tot_con_space_size
        integer :: min_element, max_element, num_elements
        integer(MPIArg) :: sendcounts(0:nProcessors-1), recvcounts(0:nProcessors-1)
        integer(MPIArg) :: senddisps(0:nProcessors-1), recvdisps(0:nProcessors-1)
        integer(n_int), allocatable, dimension(:,:) :: temp_space
        integer :: nI(nel)
        character (len=*), parameter :: t_r = "init_trial_wf"

        write(6,'()')
        write(6,'(a56)') "=========== Trial wavefunction initialisation =========="

        ! Simply allocate the trial vector to have up to 1 million elements for now...
        allocate(trial_space(0:NIfTot, 1000000), stat=ierr)
        call LogMemAlloc('trial_space', 1000000*(NIfTot+1), size_n_int, t_r, TrialTag, ierr)

        trial_space = 0
        trial_space_size = 0

        write(6,'(a29)') "Generating the trial space..."
        call neci_flush(6)

        ! Generate the trial space and place the corresponding states in trial_space.
        if (tDoublesTrial) then
            call generate_sing_doub_determinants(called_from_trial)
        elseif (tCASTrial) then
            call generate_cas(called_from_trial)
        elseif (tOptimisedTrial) then
            call generate_optimised_core(called_from_trial)
        elseif (tPopsTrial) then
            call generate_space_from_pops(called_from_trial)
        elseif (tReadTrial) then
            call generate_space_from_file(called_from_trial)
        elseif (tLowETrial) then
            call generate_low_energy_core(called_from_trial)
        end if

        if (tLimitTrialSpace) call remove_high_energy_orbs&
                                       (trial_space(:, 1:trial_space_size), &
                                        trial_space_size, max_trial_size, .true.)

        write(6,'(a38,1X,i8)') "Size of trial space on this processor:", trial_space_size
        call neci_flush(6)

        ! At this point, each processor has only those states which reside on them, and
        ! have only counted those states. Send all states to all processors for the next bit.
        call MPIAllGather(int(trial_space_size,MPIArg), recvcounts, ierr)
        tot_trial_space_size = int(sum(recvcounts),4)
        recvdisps(0) = 0
        do i = 1, nProcessors-1
            recvdisps(i) = sum(recvcounts(:i-1))
        end do

        write(6,'(a30,1X,i8)') "Total size of the trial space:", tot_trial_space_size
        call neci_flush(6)

        ! Use SpawnedParts as temporary space:
        SpawnedParts(:, 1:trial_space_size) = trial_space(:, 1:trial_space_size)
        call MPIAllGatherV(SpawnedParts(:, 1:trial_space_size), &
                           trial_space(:, 1:tot_trial_space_size), recvcounts, recvdisps)
        SpawnedParts = 0

        ! Temporarily set this to be the case.
        trial_space_size = tot_trial_space_size

        call sort(trial_space(0:NIfTot, 1:trial_space_size), ilut_lt, ilut_gt)

        allocate(trial_wf(trial_space_size), stat=ierr)
        call LogMemAlloc('trial_wf', trial_space_size, 8, t_r, TrialWFTag, ierr)
        
        call assign_elements_on_procs(trial_space_size, min_element, max_element, num_elements)

        if (num_elements > 0) then

            ! Find the states connected to the trial space. This typically takes a long time, so
            ! it is done in parallel by letting each processor find the states connected to a
            ! portion of the trial space.
            write(6,'(a58)') "Calculating the number of states in the connected space..."
            call neci_flush(6)

            call generate_connected_space(num_elements, trial_space(:, min_element:max_element), &
                                          con_space_size)

            allocate(con_space(0:NIfTot, con_space_size), stat=ierr)
            call LogMemAlloc('con_space', con_space_size*(NIfTot+1), size_n_int, t_r, ConTag, ierr)
            con_space = 0

            write(6,'(a50,1X,i8)') "States found on this processor, including repeats:", con_space_size

            write(6,'(a45)') "Generating and storing the connected space..."
            call neci_flush(6)

            call generate_connected_space(num_elements, trial_space(:, min_element:max_element), &
                                          con_space_size, con_space)

            write(6,'(a52)') "Removing repeated states and sorting by processor..."
            call neci_flush(6)

            call remove_repeated_states(con_space, con_space_size)

            call sort_space_by_proc(con_space(:, 1:con_space_size), con_space_size, sendcounts)

        else
            con_space_size = 0
            sendcounts = 0
        end if

        write(6,'(a51)') "Performing MPI communication of connected states..."
        call neci_flush(6)

        ! Send the connected states to their processors.
        ! sendcounts holds the number of states to send to other processors from this one.
        ! recvcounts will hold the number of states to be sent to this processor from the others.
        call MPIAlltoAll(sendcounts, 1, recvcounts, 1, ierr)
        con_space_size = sum(recvcounts)
        ! The displacements necessary for mpi_alltoall.
        sendcounts = sendcounts*int(NIfTot+1,MPIArg)
        recvcounts = recvcounts*int(NIfTot+1,MPIArg)
        senddisps(0) = 0
        recvdisps(0) = 0
        do i = 1, nProcessors-1
            senddisps(i) = sum(sendcounts(:i-1))
            recvdisps(i) = sum(recvcounts(:i-1))
        end do

        allocate(temp_space(0:NIfTot, con_space_size), stat=ierr)
        call LogMemAlloc('temp_space', con_space_size*(NIfTot+1), size_n_int, t_r, TempTag, ierr)

        call MPIAlltoAllV(con_space, sendcounts, senddisps, temp_space, recvcounts, recvdisps, ierr)

        deallocate(con_space, stat=ierr)
        call LogMemDealloc(t_r, ConTag, ierr)
        allocate(con_space(0:NIfTot, 1:con_space_size), stat=ierr)
        call LogMemAlloc('con_space', con_space_size*(NIfTot+1), size_n_int, t_r, ConTag, ierr)
        con_space = temp_space
        deallocate(temp_space, stat=ierr)
        call LogMemDealloc(t_r, TempTag, ierr)
        ! Finished sending to states to their processors.

        call remove_repeated_states(con_space, con_space_size)

        ! Remove states in the connected space which are also in the trial space.
        call remove_list1_states_from_list2(trial_space, con_space, &
                                            trial_space_size, con_space_size)

        call MPISumAll(con_space_size, tot_con_space_size)

        write(6,'(a30,1X,i8)') "Total size of connected space:", tot_con_space_size
        write(6,'(a42,1X,i8)') "Size of connected space on this processor:", con_space_size
        call neci_flush(6)

        ! Now the correct states in both the trial space and connected space have been fully
        ! generated. We want both lists to be sorted in the same order as CurrentDets (to use
        ! binary searching later). If tSortDetermToTop is true then deterministic states are
        ! sorted to the top in CurrentDets. Hence, we must find any deterministic states and
        ! sort them as such in the trial and connected spaces.
        if (tSortDetermToTop) then
            call find_determ_states_and_sort(trial_space, trial_space_size)
            call find_determ_states_and_sort(con_space, con_space_size)
        end if

        ! Only perform this on the root processor, due to memory demands.
        if (iProcindex == root ) then
            write(6,'(a50)') "Calculating the Hamiltonian for the trial space..."
            call neci_flush(6)
            call calculate_sparse_hamiltonian(trial_space_size, &
                                              trial_space(0:NIfTot, 1:trial_space_size))

            write(6,'(a50)') "Calculating the ground state in the trial space..."
            call neci_flush(6)
            call perform_davidson(sparse_hamil_type, .false.)

            trial_wf = davidson_eigenvector
            trial_energy = davidson_eigenvalue

            call deallocate_sparse_ham()
            deallocate(hamil_diag, stat=ierr)
            call LogMemDealloc(t_r, HDiagTag, ierr)
        else
            if (allocated(davidson_eigenvector)) then
                deallocate(davidson_eigenvector, stat=ierr)
                call LogMemDealloc(t_r, DavidsonTag, ierr)
            end if
            allocate(davidson_eigenvector(trial_space_size), stat=ierr)
            call LogMemAlloc("davidson_eigenvector", trial_space_size, 8, t_r, DavidsonTag, ierr)
        end if

        call MPIBCast(trial_energy, 1, root)
        call MPIBCast(trial_wf, size(trial_wf), root)

        write(6,'(a47)') "Generating the vector \sum_j H_{ij} \psi^T_j..."
        call neci_flush(6)
        call generate_connected_space_vector()

        ! Remove all states in trial_space that do not belong to this processor (and move the
        ! amplitudes in trial_wf at the same time).
        allocate(temp_space(0:NIfTot, trial_space_size), stat=ierr)
        call LogMemAlloc('temp_space', trial_space_size*(NIfTot+1), size_n_int, t_r, TempTag, ierr)
        temp_space = trial_space(0:NIfTot, 1:trial_space_size)
        call remove_states_not_on_proc(temp_space, trial_space_size, .true.)
        deallocate(trial_space, stat=ierr)
        call LogMemDealloc(t_r, TrialTag, ierr)
        allocate(trial_space(0:NIfTot, trial_space_size), stat=ierr)
        call LogMemAlloc('trial_space', trial_space_size*(NIfTot+1), size_n_int, t_r, TrialTag, ierr)
        trial_space = temp_space(0:NIfTot, 1:trial_space_size)

        ! Finally, correct the size of the trial_wf array. Use Davidson eigenvector
        ! as temporary space.
        davidson_eigenvector = trial_wf
        deallocate(trial_wf, stat=ierr)
        call LogMemDealloc(t_r, TrialWFTag, ierr)
        allocate(trial_wf(trial_space_size), stat=ierr)
        call LogMemAlloc('trial_wf', trial_space_size, 8, t_r, TrialWFTag, ierr)
        trial_wf = davidson_eigenvector(1:trial_space_size)

        call MPIBarrier(ierr)

        if (tWriteTrial) call write_trial_space()

        ! Deallocate remaining arrays.
        deallocate(temp_space, stat=ierr)
        call LogMemDealloc(t_r, TempTag, ierr)
        deallocate(davidson_eigenvector, stat=ierr)
        call LogMemDealloc(t_r, DavidsonTag, ierr)

        write(6,'(a37,1X,f13.7)') "Energy eigenvalue of the trial space:", trial_energy
        write(6,'(a43)') "Trial wavefunction initialisation complete."
        call neci_flush(6)

    end subroutine init_trial_wf

    subroutine remove_states_not_on_proc(ilut_list, ilut_list_size, update_trial_vector)

        integer, intent(inout) :: ilut_list_size
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, ilut_list_size)
        logical, intent(in) :: update_trial_vector
        integer :: i, counter, proc
        integer :: nI(nel)

        counter = 0
        do i = 1, ilut_list_size
            call decode_bit_det(nI, ilut_list(0:NIfTot, i))
            proc = DetermineDetNode(nI, 0)
            ! If this state and the previous one were identical, don't add this state to the
            ! list so that repeats aren't included.
            if (proc == iProcIndex) then
                counter = counter + 1
                ilut_list(:, counter) = ilut_list(:, i)
                if (update_trial_vector) trial_wf(counter) = trial_wf(i)
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
            pos = binary_search(list_1(:, min_ind:list_1_size), list_2(:,i), NifD+1)
            ! If it is in list 1, remove the state by setting it to 0.
            ! If it isn't in list 1 (pos < 0) then we can still search a smaller list next time.
            if (pos > 0) then
                list_2(:,i) = 0
                min_ind = min_ind + pos
            else
                min_ind = min_ind - pos - 1
            end if
        end do

        ! Now compress the new list by overwriting the removed states:
        counter = 0
        do i = 1, list_2_size
            ! If the state wasn't set to 0:
            if (.not. all(list_2(:,i) == 0)) then
                counter = counter + 1
                list_2(:, counter) = list_2(:, i)
            end if
        end do

     list_2_size = counter

    end subroutine remove_list1_states_from_list2

    subroutine generate_connected_space_vector()

        ! Calculate the vector
        ! \sum_j H_{ij} \psi^T_j,
        ! where \psi^T is the trial vector, j runs over all trial space vectors and i runs over
        ! all connected space vectors. This quantity is stored in the con_space_vector array.

        integer :: i, j
        integer :: nI(nel), nJ(nel)
        real(dp) :: H_ij

        allocate(con_space_vector(con_space_size))
        con_space_vector = 0.0_dp

        do i = 1, con_space_size
            call decode_bit_det(nI, con_space(0:NIfTot, i))
            do j = 1, trial_space_size
                call decode_bit_det(nJ, trial_space(0:NIfTot, j))
                ! Note that, because the connected and trial spaces do not contain any common
                ! states, we never have diagonal Hamiltonian elements.
                if (.not. tHPHF) then
                    H_ij = get_helement(nI, nJ, con_space(:,i), trial_space(:,j))
                else
                    H_ij = hphf_off_diag_helement(nI, nJ, con_space(:,i), trial_space(:,j))
                end if
                con_space_vector(i) = con_space_vector(i) + H_ij*trial_wf(j)
            end do
        end do

    end subroutine generate_connected_space_vector

    subroutine assign_elements_on_procs(list_length, min_element, max_element, num_elements)

        ! Split list_length into nProcessor parts. Note that this is not done based on any hash.

        integer, intent(in) :: list_length
        integer, intent(out) :: min_element, max_element, num_elements
        integer :: floor_div_list_length, mod_list_length
        integer :: num_elem_all_procs(0:nProcessors-1)
        integer :: i

        mod_list_length = mod(list_length, nProcessors)
        floor_div_list_length = (list_length-mod_list_length)/nProcessors

        do i = 0, nProcessors-1
            num_elem_all_procs(i) = floor_div_list_length
            if (i < mod_list_length) num_elem_all_procs(i) = num_elem_all_procs(i) + 1
        end do

        num_elements = num_elem_all_procs(iProcIndex)

        if (num_elements == 0) then
            if (iProcIndex == 0) call stop_all("assign_elements_on_procs","There are no states &
                                               &in the trial space.")
            min_element = 0
            max_element = 0
            return
        end if

        if (iProcIndex == 0) then
            min_element = 1
            max_element = num_elements
        else
            min_element = sum(num_elem_all_procs(0:iProcIndex-1)) + 1
            max_element = sum(num_elem_all_procs(0:iProcIndex))
        end if

    end subroutine assign_elements_on_procs

    subroutine write_trial_space()

        integer :: i, j, k, iunit, ierr

        write(6,'(a36)') "Writing the trial space to a file..."

        ! Let each processor write its trial states to the file. Each processor waits for
        ! the processor before it to finish before starting.
        do i = 0, nProcessors-1

            if (iProcIndex == i) then

                if (i == 0) then
                    iunit = get_free_unit()
                    open(iunit, file='TRIALSPACE', status='replace')
                else
                    open(iunit, file='TRIALSPACE', status='old', position='append')
                end if
                
                do j = 1, trial_space_size 
                    do k = 0, NIfDBO
                        write(iunit, '(i24)', advance = 'no') trial_space(k,j)
                    end do
                    write(iunit, *)
                end do

                close(iunit)

            end if

            call MPIBarrier(ierr)

            if (i == 0) call MPIBCast(iunit, 0)

        end do

    end subroutine write_trial_space

end module trial_wf_gen
