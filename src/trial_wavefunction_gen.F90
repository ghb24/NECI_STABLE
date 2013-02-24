module trial_wavefunction_gen

    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: encode_det
    use CalcData, only: tSortDetermToTop
    use davidson, only: perform_davidson, davidson_eigenvalue, davidson_eigenvector, &
                        sparse_hamil_type
    use DetBitOps, only: FindBitExcitLevel
    use DeterminantData, only: write_det
    use FciMCData, only: trial_space, trial_space_size, connected_space, &
                         connected_space_size, trial_wavefunction, trial_energy, &
                         connected_space_vector, ilutHF, Hii, nSingles, nDoubles
    use hphf_integrals, only: hphf_off_diag_helement
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier, MPIAllToAll, MPISumAll, &
                             MPIAllToAllV, MPIArg
    use ParallelHelper, only: root
    use semi_stochastic
    use SystemData, only: nel, tDoublesTrial, tOptimisedTrial, tCASTrial, tHPHF

    implicit none

contains

    subroutine init_trial_wavefunction()

        integer :: i, ierr, num_states_on_proc
        integer :: excit, tot_trial_space_size, tot_connected_space_size
        integer :: connected_storage_space_size, min_element, max_element, num_elements
        integer(MPIArg) :: sendcounts(0:nProcessors-1), recvcounts(0:nProcessors-1)
        integer(MPIArg) :: senddisps(0:nProcessors-1), recvdisps(0:nProcessors-1)
        integer(n_int), allocatable, dimension(:,:) :: temp_space
        integer :: nI(nel)
        integer :: old, new

        write(6,'()')
        write(6,'(a56)') "=========== Trial wavefunction initialisation =========="

        ! Simply allocate the trial vector to have up to 1 million elements for now...
        allocate(trial_space(0:NIfTot, 1000000))
        trial_space = 0

        ! Encode the Hartree-Fock state first.
        call encode_det(trial_space(:,1), ilutHF)
        trial_space_size = 1

        ! Generate the trial space and place the corresponding states in trial_space.
        if (tDoublesTrial) then
            call generate_sing_doub_determinants(called_from_trial)
        elseif (tCASTrial) then
            call generate_cas(called_from_trial)
        elseif (tOptimisedTrial) then
            call generate_optimised_core(called_from_trial)
        end if

        call sort(trial_space(0:NIfTot, 1:trial_space_size), ilut_lt, ilut_gt)
        
        write(6,'(a50)') "Calculating the Hamiltonian for the trial space..."
        call neci_flush(6)
        call calculate_sparse_hamiltonian(trial_space_size, &
            trial_space(0:NIfTot, 1:trial_space_size))

        ! To allocate storage space for the connected space states, assume that each state
        ! in the trial space has roughly the same number as connected states as the HF state
        ! (nSingles+nDoubles). Also divide by the number of processors.
        connected_storage_space_size = &
            ceiling(real(trial_space_size)/real(nProcessors))*(nSingles+nDoubles)
        allocate(connected_space(0:NIfTot, connected_storage_space_size))
        connected_space = 0
        
        call find_min_max_elements_on_proc(trial_space_size, min_element, max_element)
        num_elements = max_element-min_element+1

        ! Find the states connected to the trial space. This typically takes a long time, so
        ! it is done in parallel by letting each processor find the states connected to a
        ! portion of the trial space.
        write(6,'(a33)') "Generating the connected space..."
        call neci_flush(6)
        call generate_connected_space(num_elements, trial_space(:, min_element:max_element), &
                                      connected_space_size, connected_space, &
                                      connected_storage_space_size)

        call remove_repeated_states(connected_space, connected_space_size)

        call sort_space_by_proc(connected_space(:, 1:connected_space_size), connected_space_size, &
                                sendcounts)

        ! Send the connected states to their processors.
        ! sendcounts holds the number of states to send to other processors from this one.
        ! recvcounts will hold the number of states to be sent to this processor from the others.
        call MPIAlltoAll(sendcounts, 1, recvcounts, 1, ierr)
        connected_space_size = sum(recvcounts)
        ! The displacements necessary for mpi_alltoall.
        sendcounts = sendcounts*int(NIfTot+1,MPIArg)
        recvcounts = recvcounts*int(NIfTot+1,MPIArg)
        senddisps(0) = 0
        recvdisps(0) = 0
        do i = 1, nProcessors-1
            senddisps(i) = sum(sendcounts(:i-1))
            recvdisps(i) = sum(recvcounts(:i-1))
        end do

        allocate(temp_space(0:NIfTot, connected_space_size))
        call MPIAlltoAllV(connected_space, sendcounts, senddisps, temp_space, recvcounts, recvdisps, ierr)
        deallocate(connected_space)
        allocate(connected_space(0:NIfTot, 1:connected_space_size))
        connected_space = temp_space
        deallocate(temp_space)

        call remove_repeated_states(connected_space, connected_space_size)

        ! Remove states in the connected space which are also in the trial space:
        call remove_list1_states_from_list2(trial_space, connected_space, &
            trial_space_size, connected_space_size)

        tot_trial_space_size = trial_space_size
        call MPISumAll(connected_space_size, tot_connected_space_size)

        ! Now the correct states in both the trial space and connected space have been fully
        ! generated. We want both lists to be sorted in the same order as CurrentDets (to use
        ! binary searching later). If tSortDetermToTop is true then deterministic states are
        ! sorted to the top in CurrentDets. Hence, we must find any deterministic states and
        ! sort them as such in the trial and connected spaces.
        if (tSortDetermToTop) then
            call find_determ_states_and_sort(trial_space, trial_space_size)
            call find_determ_states_and_sort(connected_space, connected_space_size)
        end if

        write(6,'(a50)') "Calculating the ground state in the trial space..."
        call neci_flush(6)
        call perform_davidson(sparse_hamil_type, .false.)
        allocate(trial_wavefunction(trial_space_size))
        trial_wavefunction = davidson_eigenvector
        trial_energy = davidson_eigenvalue

        write(6,'(a47)') "Generating the vector \sum_j H_{ij} \psi^T_j..."
        call neci_flush(6)
        call generate_connected_space_vector()

        ! Remove all states in trial_space that do not belong to this processor (and move the
        ! amplitudes in trial_wavefunction at the same time).
        allocate(temp_space(0:NIfTot, trial_space_size))
        temp_space = trial_space(0:NIfTot, 1:trial_space_size)
        call remove_states_not_on_proc(temp_space, trial_space_size, .true.)
        deallocate(trial_space)
        allocate(trial_space(0:NIfTot, trial_space_size))
        trial_space = temp_space(0:NIfTot, 1:trial_space_size)
        deallocate(temp_space)

        ! Finally, correct the size of the trial_wavefunction array. Use Davidson eigenvector
        ! as temporary space.
        davidson_eigenvector = trial_wavefunction
        deallocate(trial_wavefunction)
        allocate(trial_wavefunction(trial_space_size))
        trial_wavefunction = davidson_eigenvector(1:trial_space_size)
        deallocate(davidson_eigenvector)

        call MPIBarrier(ierr)

        write(6,'(a30,1X,i8)') "Total size of the trial space:", tot_trial_space_size
        write(6,'(a38,1X,i8)') "Size of trial space on this processor:", trial_space_size
        write(6,'(a30,1X,i8)') "Total size of connected space:", tot_connected_space_size
        write(6,'(a42,1X,i8)') "Size of connected space on this processor:", connected_space_size
        write(6,'(a37,1X,f13.7)') "Energy eigenvalue of the trial space:", trial_energy
        write(6,'(a42)') "Trial wavefunction initialisation complete"
        call neci_flush(6)

    end subroutine init_trial_wavefunction

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
                if (update_trial_vector) trial_wavefunction(counter) = trial_wavefunction(i)
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
        ! all connected space vectors. This quantity is stored in the connected_space_vector array.

        integer :: i, j
        integer :: nI(nel), nJ(nel)
        real(dp) :: H_ij

        allocate(connected_space_vector(connected_space_size))
        connected_space_vector = 0.0_dp

        do i = 1, connected_space_size
            call decode_bit_det(nI, connected_space(0:NIfTot, i))
            do j = 1, trial_space_size
                call decode_bit_det(nJ, trial_space(0:NIfTot, j))
                ! Note that, because the connected and trial spaces do not contain any common
                ! states, we never have diagonal Hamiltonian elements.
                if (.not. tHPHF) then
                    H_ij = get_helement(nI, nJ, connected_space(:,i), trial_space(:,j))
                else
                    H_ij = hphf_off_diag_helement(nI, nJ, connected_space(:,i), trial_space(:,j))
                end if
                connected_space_vector(i) = connected_space_vector(i) + H_ij*trial_wavefunction(j)
            end do
        end do

    end subroutine generate_connected_space_vector

    subroutine find_min_max_elements_on_proc(list_length, min_element, max_element)

        ! Split list_length into nProcessor parts. Not that this is not done based on any hash.

        integer, intent(in) :: list_length
        integer, intent(out) :: min_element, max_element
        integer :: floor_divided_list_length, mod_list_length

        mod_list_length = mod(list_length, nProcessors)
        floor_divided_list_length = (list_length-mod_list_length)/nProcessors

        min_element = floor_divided_list_length*iProcIndex + 1
        max_element = floor_divided_list_length*(iProcIndex+1)
        if ((iProcIndex-1) < mod_list_length .and. iProcIndex > 0) min_element = min_element + 1
        if (iProcIndex < mod_list_length) max_element = max_element + 1

    end subroutine find_min_max_elements_on_proc

    subroutine sort_space_by_proc(ilut_list, ilut_list_size, num_states_procs)

        ! And also output the number of states on each processor in the space...

        integer, intent(in) :: ilut_list_size
        integer(n_int) :: ilut_list(0:NIfTot, 1:ilut_list_size)
        integer(MPIArg), intent(out) :: num_states_procs(0:nProcessors-1)
        integer :: proc_list(ilut_list_size)
        integer :: nI(nel)
        integer :: i

        num_states_procs = 0

        ! Create a list, proc_list, with the processor numbers of the corresponding iluts.
        do i = 1, ilut_list_size
            call decode_bit_det(nI, ilut_list(:,i))
            proc_list(i) = DetermineDetNode(nI,0)
            num_states_procs(proc_list(i)) = num_states_procs(proc_list(i)) + 1
        end do

        ! Now sort the list according to this processor number.
        call sort(proc_list, ilut_list(0:NIfTot, 1:ilut_list_size))

    end subroutine sort_space_by_proc

end module trial_wavefunction_gen
