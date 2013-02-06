module trial_wavefunction_gen

    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: encode_det
    use CalcData, only: tSortDetermToTop
    use davidson, only: perform_davidson, davidson_eigenvalue, davidson_eigenvector, &
                        sparse_hamil_type
    use DetBitOps, only: FindBitExcitLevel
    use DeterminantData, only: write_det
    use FciMCData, only: trial_vector_space, trial_vector_space_size, connected_space, &
                         connected_space_size, trial_wavefunction, trial_energy, &
                         connected_space_vector, ilutHF, Hii
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier
    use ParallelHelper, only: root
    use semi_stochastic
    use SystemData, only: nel, tDoublesTrial, tOptimisedTrial, tCASTrial

    implicit none

contains

    subroutine init_trial_wavefunction()

        integer :: i, ierr, counter, comp, num_states_on_proc
        integer(n_int), allocatable, dimension(:,:) :: temp_space
        integer :: nI(nel)
        integer :: excit, total_trial_space_size, total_connected_space_size

        write(6,'()')
        write(6,'(a56)') "=========== Trial wavefunction initialisation =========="

        ! Simply allocate the trial vector to have up to 1 million elements for now...
        allocate(trial_vector_space(0:NIfTot, 1000000))
        allocate(connected_space(0:NIfTot, 1000000))
        trial_vector_space = 0
        connected_space = 0

        ! Encode the Hartree-Fock state first.
        call encode_det(trial_vector_space(:,1), ilutHF)
        trial_vector_space_size = 1

        ! Generate the trial space and place the corresponding states in trial_vector_space.
        if (tDoublesTrial) then
            call generate_sing_doub_determinants(called_from_trial)
        elseif (tCASTrial) then
            call generate_cas(called_from_trial)
        elseif (tOptimisedTrial) then
            call generate_optimised_core(called_from_trial)
        end if

        call sort(trial_vector_space(0:NIfTot, 1:trial_vector_space_size), ilut_lt, ilut_gt)
        
        write(6,'(a50)') "Calculating the Hamiltonian for the trial space..."
        call neci_flush(6)
        call calculate_sparse_hamiltonian(trial_vector_space_size, &
            trial_vector_space(0:NIfTot, 1:trial_vector_space_size))

        ! Find the states connected to the trial space.
        write(6,'(a33)') "Generating the connected space..."
        call neci_flush(6)
        call generate_connected_space(trial_vector_space_size, trial_vector_space, &
            connected_space_size, connected_space)

        ! Annihilation-like steps to remove repeated states.
        call sort(connected_space(0:NIfTot, 1:connected_space_size), ilut_lt, ilut_gt)
        counter = 1
        do i = 2, connected_space_size
            comp = DetBitLT(connected_space(:, i-1), connected_space(:, i), NIfD, .false.)
            ! If this state and the previous one were identical, don't add this state to the
            ! list so that repeats aren't included.
            if (comp /= 0) then
                counter = counter + 1
                connected_space(:, counter) = connected_space(:, i)
            end if
        end do

        connected_space_size = counter

        ! Now we want to remove states in the connected space which are also in the trial space:
        call remove_list_1_states_from_list_2(trial_vector_space, connected_space, &
            trial_vector_space_size, connected_space_size)

        total_trial_space_size = trial_vector_space_size
        total_connected_space_size = connected_space_size

        ! Remove all states in connected_space that do not belong to this processor.
        allocate(temp_space(0:NIfTot, connected_space_size))
        temp_space = connected_space(0:NIfTot, 1:connected_space_size)
        call remove_states_not_on_proc(temp_space, connected_space_size, .false.)
        deallocate(connected_space)
        ! Reallocate connected_space with the new number of states, now that those states on
        ! different processors have been removed.
        allocate(connected_space(0:NIfTot, connected_space_size))
        connected_space = temp_space(0:NIfTot, 1:connected_space_size)

        ! Now the correct states in both the trial space and connected space have been fully
        ! generated. We want both lists to be sorted in the same order as CurrentDets, to use
        ! binary searching later. If tSortDetermToTop is true then deterministic states are
        ! sorted to the top in CurrentDets. Hence, we must find any deterministic states and
        ! sort them as such in the trial and connected spaces.
        if (tSortDetermToTop) then
            call find_determ_states_and_sort(trial_vector_space, trial_vector_space_size)
            call find_determ_states_and_sort(connected_space, connected_space_size)
        end if

        ! Use the Davidson method to find the ground state of the trial space.
        write(6,'(a50)') "Calculating the ground state in the trial space..."
        call neci_flush(6)
        call perform_davidson(sparse_hamil_type, .false.)
        allocate(trial_wavefunction(trial_vector_space_size))
        trial_wavefunction = davidson_eigenvector
        trial_energy = davidson_eigenvalue

        write(6,'(a47)') "Generating the vector \sum_j H_{ij} \psi^T_j..."
        call neci_flush(6)
        call generate_connected_space_vector()

        ! Remove all states in trial_vector_space that do not belong to this processor (and
        ! move the amplitudes in trial_wavefunction at the same time).
        deallocate(temp_space)
        allocate(temp_space(0:NIfTot, trial_vector_space_size))
        temp_space = trial_vector_space(0:NIfTot, 1:trial_vector_space_size)
        call remove_states_not_on_proc(temp_space, trial_vector_space_size, .true.)
        deallocate(trial_vector_space)
        allocate(trial_vector_space(0:NIfTot, trial_vector_space_size))
        trial_vector_space = temp_space(0:NIfTot, 1:trial_vector_space_size)

        ! Finally, correct the size of trial_wavefunction, too. Use Davidson eigenvector as
        ! temporary space.
        davidson_eigenvector = trial_wavefunction
        deallocate(trial_wavefunction)
        allocate(trial_wavefunction(trial_vector_space_size))
        trial_wavefunction = davidson_eigenvector(1:trial_vector_space_size)
        deallocate(davidson_eigenvector)

        call MPIBarrier(ierr)

        write(6,'(a30,1X,i8)') "Total size of the trial space:", total_trial_space_size
        write(6,'(a38,1X,i8)') "Size of trial space on this processor:", trial_vector_space_size
        write(6,'(a30,1X,i8)') "Total size of connected space:", total_connected_space_size
        write(6,'(a42,1X,i8)') "Size of connected space on this processor:", connected_space_size
        write(6,'(a37,1X,f12.8)') "Energy eigenvalue of the trial space:", trial_energy
        write(6,'(a42)') "Trial wavefunction initialisation complete"
        call neci_flush(6)

    end subroutine init_trial_wavefunction

    subroutine count_states_on_this_proc(num_states_on_proc, ilut_list_size, ilut_list)

        integer, intent(out) :: num_states_on_proc
        integer, intent(in) :: ilut_list_size
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, ilut_list_size)
        integer :: nI(nel)
        integer :: i, proc

        num_states_on_proc = 0

        do i = 1, ilut_list_size

            call decode_bit_det(nI, trial_vector_space(0:NIfTot, i))
            ! Find the processor which this current state belongs to.
            proc = DetermineDetNode(nI, 0)
            ! If on this processor.
            if (proc == iProcIndex) num_states_on_proc = num_states_on_proc + 1

        end do

    end subroutine count_states_on_this_proc

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

    subroutine remove_list_1_states_from_list_2(ilut_list_1, ilut_list_2, list_1_size, list_2_size)

        use util_mod, only: binary_search

        integer, intent(in) :: list_1_size
        integer, intent(inout) :: list_2_size
        integer(n_int), intent(in) :: ilut_list_1(0:NIfTot, list_1_size)
        integer(n_int), intent(inout) :: ilut_list_2(0:NIfTot, list_2_size)
        integer :: i, counter, pos, min_ind

        min_ind = 1

        do i = 1, list_2_size
            ! Binary search ilut_list_1 to see if ilut_list_2(:,i) is in it.
            pos = binary_search(ilut_list_1(:, min_ind:list_1_size), ilut_list_2(:,i), NifD+1)
            ! If it is in list 1, remove the state by setting it to 0.
            ! If it isn't in list 1 (pos < 0) then we can still search a smaller list next time.
            if (pos > 0) then
                ilut_list_2(:,i) = 0
                min_ind = min_ind + pos
            else
                min_ind = min_ind - pos - 1
            end if
        end do

        ! Now compress the new list by overwriting the removed states:
        counter = 0
        do i = 1, list_2_size
            ! If the state wasn't set to 0:
            if (.not. all(ilut_list_2(:,i) == 0)) then
                counter = counter + 1
                ilut_list_2(:, counter) = ilut_list_2(:, i)
            end if
        end do

     list_2_size = counter

    end subroutine remove_list_1_states_from_list_2

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
            do j = 1, trial_vector_space_size
                call decode_bit_det(nJ, trial_vector_space(0:NIfTot, j))
                ! Note that, because the connected and trial spaces do not contain any common
                ! states, we never have diagonal Hamiltonian elements.
                H_ij = get_helement(nI, nJ, connected_space(:,i), trial_vector_space(:,j))
                connected_space_vector(i) = connected_space_vector(i) + H_ij*trial_wavefunction(j)
            end do
        end do

    end subroutine generate_connected_space_vector

end module trial_wavefunction_gen
