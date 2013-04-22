#include "macros.h"

! Some general procedures, created for the semi-stochastic (and trial wavefunction) code.
! Some are just used for initialisation and others in the main FCIQMC loop.

module semi_stoch_procs

    use bit_rep_data, only: flag_deterministic, nIfDBO, NIfD, NIfTot
    use bit_reps, only: decode_bit_det, set_flag, extract_part_sign
    use CalcData, only: tRegenDiagHEls, tau, DiagSft
    use constants
    use DetBitOps, only: ilut_lt, ilut_gt, FindBitExcitLevel, DetBitLT, &
                         count_set_bits, DetBitEq
    use Determinants, only: get_helement
    use FciMCData, only: ilutHF, Hii, CurrentH, determ_proc_sizes, determ_proc_indices, &
                         full_determ_vector, partial_determ_vector, core_hamiltonian, &
                         determ_space_size, SpawnedParts, SemiStoch_Comms_Time, &
                         SemiStoch_Multiply_Time
    use hash, only: DetermineDetNode
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIBarrier, MPIArg, &
                             MPIAllGatherV, MPISum, MPISumAll
    use ras, only: core_ras
    use sort_mod, only: sort
    use SystemData, only: tSemiStochastic, tCSFCore, tDeterminantCore, tDoublesCore, &
                          tCASCore, tRASCore, cas_determ_not_bitmask, core_ras1_bitmask, &
                          core_ras3_bitmask, nel, OccDetermCASOrbs, tHPHF, nBasis, BRR, &
                          ARR
    use timing_neci

    implicit none

contains

    subroutine deterministic_projection()

        ! This subroutine gathers together the partial_determ_vectors from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        integer :: i, info

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_determ_vector, full_determ_vector, determ_proc_sizes, &
                            determ_proc_indices)

        call halt_timer(SemiStoch_Comms_Time)

        call set_timer(SemiStoch_Multiply_Time)

        if (determ_proc_sizes(iProcIndex) >= 1) then

            ! This function performs y := alpha*A*x + beta*y
            ! N specifies not to use the transpose of A.
            ! determ_proc_sizes(iProcIndex) is the number of rows in A.
            ! determ_space_size is the number of columns of A.
            ! alpha = -1.0_dp.
            ! A = core_hamiltonian.
            ! determ_proc_sizes(iProcIndex) is the first dimension of A.
            ! input x = full_determ_vector.
            ! 1 is the increment of the elements of x.
            ! beta = 0.0_dp.
            ! output y = partial_determ_vector.
            ! 1 is the incremenet of the elements of y.
            call dgemv('N', &
                       determ_proc_sizes(iProcIndex), &
                       determ_space_size, &
                       -1.0_dp, &
                       core_hamiltonian, &
                       determ_proc_sizes(iProcIndex), &
                       full_determ_vector, &
                       1, &
                       0.0_dp, &
                       partial_determ_vector, &
                       1)

            ! Now add shift*full_determ_vector, to account for the shift, not stored in
            ! core_hamiltonian.
            partial_determ_vector = partial_determ_vector + &
               DiagSft * full_determ_vector(determ_proc_indices(iProcIndex)+1:&
                 determ_proc_indices(iProcIndex)+determ_proc_sizes(iProcIndex))

            ! Now multiply the vector by tau to get the final projected vector.
            partial_determ_vector = partial_determ_vector * tau

            call halt_timer(SemiStoch_Multiply_Time)

        end if

    end subroutine deterministic_projection

    subroutine check_if_in_determ_space(ilut, tInDetermSpace)

        ! In: ilut - The state that we will check, in a bitstring form.
        ! Out: tInDetermSpace - This will be set to true if the state (ilut) is in the
        !                       deterministic space, and false if not.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer(n_int) :: ilut_temp1(0:NIfD), ilut_temp2(0:NIfD)
        logical, intent(out) :: tInDetermSpace
        integer :: IC, bits_set

        tInDetermSpace = .false.

        if (tDoublesCore) then
            if (tDeterminantCore) then
                IC = FindBitExcitLevel(ilut, ilutHF)
                if (IC <= 2) tInDetermSpace = .true.
            else if (tCSFCore) then
                IC = FindBitExcitLevel(ilut, ilutHF)
                if (IC <= 2) tInDetermSpace = .true.
            end if
        else if (tCASCore) then
            ! The state ilut, but with all active space orbitals unoccupied.
            ilut_temp1 = iand(cas_determ_not_bitmask, ilut(0:NIfD))

            ! The state ilut, but with all active space orbitals and all orbitals not occupied
            ! in the HF determinant unoccupied.
            ilut_temp2 = iand(ilut_temp1, ilutHF(0:NIfD))
            ! All these orbitals should be occupied if the state is in the deterministic space.
            bits_set = sum(count_set_bits(ilut_temp2))
            if (bits_set /= nel - OccDetermCASOrbs) return

            ! The state ilut, but with all active space orbitals and all orbitals in the HF
            ! determinant unoccupied.
            ilut_temp2 = iand(ilut_temp1, not(ilutHF(0:NIfD)))
            ! All these orbitals should be unoccupied if the state is in the deterministic space.
            bits_set = sum(count_set_bits(ilut_temp2))
            if (bits_set /= 0) return

            ! If both of the above occupation criteria are met, then the state is in the
            ! deterministic space.
            tInDetermSpace = .true.
        else if (tRASCore) then
            ! Check that the minimum number of electrons in RAS1 is obeyed.
            ilut_temp1 = iand(core_ras1_bitmask, ilut(0:NIfD))
            bits_set = sum(count_set_bits(ilut_temp1))
            if (bits_set < core_ras%min_1) return

            ! Check that the maximum number of electrons in RAS3 is not exceeded.
            ilut_temp2 = iand(core_ras3_bitmask, ilut(0:NIfD))
            bits_set = sum(count_set_bits(ilut_temp2))
            if (bits_set > core_ras%max_3) return

            ! If both of the above conditions were met.
            tInDetermSpace = .true.
        end if

    end subroutine check_if_in_determ_space

    subroutine calculate_det_hamiltonian_normal()

        integer :: i, j, iproc, col_index, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        real(dp), allocatable, dimension(:,:) :: test_hamiltonian
        integer(TagIntType) :: TempStoreTag
        character (len=*), parameter :: this_routine = "calculate_det_hamiltonian_normal"

        ! temp_store is storage space for bitstrings so that the Hamiltonian matrix
        ! elements can be calculated.
        allocate(temp_store(0:NIfTot, maxval(determ_proc_sizes)), stat=ierr)
        call LogMemAlloc('temp_store', maxval(determ_proc_sizes)*(NIfTot+1), 8, &
                         this_routine, TempStoreTag, ierr)

        ! Loop over each processor in order and broadcast the sorted vector of bitstrings
        ! from the processor to all other processors so that all matrix elements can be found.
        do iproc = 0, nProcessors-1

            ! If we are broadcasting from this processor next, transfer the bitstrings
            ! to the array temp_store.
            if (iproc == iProcIndex) temp_store(:,1:determ_proc_sizes(iproc)) = &
                                      SpawnedParts(:,1:determ_proc_sizes(iproc))

            ! Perform the broadcasting to other all other processors.
            call MPIBCast(temp_store, size(temp_store), iproc)

            ! The index of the column before the first column of the block of the
            ! Hamiltonian currently being calculated.
            col_index = determ_proc_indices(iproc)

            ! Loop over all the elements in the block of the Hamiltonian corresponding
            ! to these two prcoessors.
            do i = 1, determ_proc_sizes(iProcIndex)
                do j = 1, determ_proc_sizes(iproc)

                    ! Get nI form of the basis functions.
                    call decode_bit_det(nI, SpawnedParts(:, i))
                    call decode_bit_det(nJ, temp_store(:, j))

                    ! If on the diagonal of the Hamiltonian.
                    if ((iProcIndex == iproc) .and. (i == j)) then
                        if (tHPHF) then
                            core_hamiltonian(i, col_index + j) = &
                                hphf_diag_helement(nI, SpawnedParts(:,i)) - Hii
                        else
                            core_hamiltonian(i, col_index + j) = &
                                get_helement(nI, nJ, 0) - Hii
                        end if
                        ! We calculate and store CurrentH at this point for ease.
                        if (.not.tRegenDiagHEls) CurrentH(1,i) = &
                                             core_hamiltonian(i, col_index + j)
                    else
                        if (tHPHF) then
                            core_hamiltonian(i, col_index + j) = &
                                hphf_off_diag_helement(nI, nJ, SpawnedParts(:,i), &
                                temp_store(:,j))
                        else
                            core_hamiltonian(i, col_index + j) = &
                                get_helement(nI, nJ, SpawnedParts(:, i), temp_store(:, j))
                        end if
                    end if

                end do
            end do

        end do

        call MPIBarrier(ierr)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(this_routine, TempStoreTag, ierr)

    end subroutine calculate_det_hamiltonian_normal

    subroutine check_if_connected_to_old_space(ilut_store_old, nI, ilut, num_states, &
                                               counter, connected)

        ! This subroutine loops over all states in the old space and sees if any of them is
        ! connected to the state ilut. If so, it also adds one to counter.

        integer(n_int), intent(in) :: ilut_store_old(0:NIfTot, 1000000)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: num_states
        integer, intent(inout) :: counter
        logical, intent(inout) :: connected
        integer :: nJ(nel)
        real(dp) :: matrix_element
        integer :: i

        connected = .false.

        ! Loop over all states in the old space and see if any of them is connected
        ! to the state just generated.
        do i = 1, num_states
            ! If this is true then the state is *in* the old space. These states
            ! aren't required in the defined space generated, so we can skip them
            ! (although they can be included too, as most will be, and removed afterwards
            ! with an annihilation-like step, for simplicity).
            if (DetBitEQ(ilut, ilut_store_old(:, i), NIfDBO)) return

            call decode_bit_det(nJ, ilut_store_old(0:NIfD, i))

            if (.not. tHPHF) then
                matrix_element = get_helement(nI, nJ, ilut, ilut_store_old(:, i))
            else
                matrix_element = hphf_off_diag_helement(nI, nJ, ilut, ilut_store_old(:, i))
            end if

            ! If true, then this state should be added to the new ilut store, so
            ! return connected as true.
            if (abs(matrix_element) > 0.0_dp) then
                counter = counter + 1 
                connected = .true.
                return
            end if
        end do

    end subroutine check_if_connected_to_old_space

    subroutine remove_repeated_states(list, list_size)

        integer, intent(inout) :: list_size
        integer(n_int), intent(inout) :: list(0:NIfTot, list_size)
        integer :: i, counter, comp

        ! Annihilation-like steps to remove repeated states.
        call sort(list(:, 1:list_size), ilut_lt, ilut_gt)
        counter = 1
        do i = 2, list_size
            comp = DetBitLT(list(:, i-1), list(:, i), NIfD, .false.)
            ! If this state and the previous one were identical, don't add this state to the
            ! list so that repeats aren't included.
            if (comp /= 0) then
                counter = counter + 1
                list(:, counter) = list(:, i)
            end if
        end do

        list_size = counter

    end subroutine remove_repeated_states

    subroutine remove_high_energy_orbs(ilut_list, num_states, target_num_states, tParallel)

        integer, intent(inout) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:num_states)
        integer, intent(in) :: target_num_states
        logical, intent(in) :: tParallel
        integer, allocatable, dimension(:) :: orbitals_rmvd
        integer :: i, j, orb, num_orbs_rmvd
        integer :: bit, elem
        integer :: states_rmvd_this_proc, counter
        integer(MPIArg) :: tot_num_states, states_rmvd_all_procs
        logical :: occupied

        states_rmvd_this_proc = 0

        if (tParallel) then
            call MPISumAll(int(num_states, MPIArg), tot_num_states)
        else
            tot_num_states = num_states
        end if

        if (tot_num_states <= target_num_states) return

        write(6,'(a32)') "Removing high energy orbitals..."

        ! Loop through all orbitals, from highest energy to lowest.
        do i = nBasis, 1, -1

            num_orbs_rmvd = nBasis-i+1

            orb = BRR(i)

            bit = mod(orb-1, bits_n_int)
            elem = (orb-1-bit)/bits_n_int

            ! Loop through all states and remove those states with orbital orb
            ! occupied.
            do j = 1, num_states
                occupied = btest(ilut_list(elem, j), bit)
                if (occupied) then
                    ilut_list(:, j) = 0
                    states_rmvd_this_proc = states_rmvd_this_proc + 1
                end if
            end do

            if (tParallel) then
                call MPISumAll(int(states_rmvd_this_proc, MPIArg), states_rmvd_all_procs)
            else
                states_rmvd_all_procs = states_rmvd_this_proc
            end if

            ! If there are degenerate orbitals, then cycle to remove the
            ! degenerate orbitals too, before giving the program chance to quit.
            if (i > 1) then
                ! If the orbitals energies are the same:
                if ( Arr(i, 1) == Arr((i-1), 1) ) cycle
            end if

            if (tot_num_states-states_rmvd_all_procs <= target_num_states) exit

        end do

        ! Loop through the list and shuffle states down to fill in the gaps
        ! created above.
        counter = 0
        do i = 1, num_states
            ! If the state wasn't set to 0:
            if (.not. all(ilut_list(:,i) == 0)) then
                counter = counter + 1
                ilut_list(:, counter) = ilut_list(:, i)
            end if
        end do
        if (counter /= num_states-states_rmvd_this_proc) &
            call stop_all("remove_high_energy_orbs", &
                          "Wrong number of states found.")
        num_states = counter

        ! Finally, output information:
        write(6,'(a36)',advance='no') "The following orbitals were removed:"
        allocate(orbitals_rmvd(num_orbs_rmvd))
        orbitals_rmvd = BRR(nBasis-num_orbs_rmvd+1:nBasis)
        call sort(orbitals_rmvd)
        do i = 1, num_orbs_rmvd
            write(6,'(i5)',advance='no') orbitals_rmvd(i)
            if (i == num_orbs_rmvd) write(6,'()',advance='yes')
        end do
        deallocate(orbitals_rmvd)
        write(6,'(i6,1x,a29)') states_rmvd_all_procs, "states were removed in total."
        write(6,'(i6,1x,a17)') tot_num_states-states_rmvd_all_procs, "states were kept."
        call neci_flush(6)

    end subroutine remove_high_energy_orbs

    subroutine sort_states_by_energy(ilut_list, num_states, tSortDoubles)

        ! Note: If requested, keep all doubles at the top, then sort by energy.

        use bit_reps, only: decode_bit_det
        use Determinants, only: get_helement
        use hphf_integrals, only: hphf_diag_helement
        use SystemData, only: tHPHF

        integer, intent(in) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:num_states)
        logical, intent(in) :: tSortDoubles
        integer(n_int) :: temp_ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: i, excit_level, num_sing_doub, block_size
        real(dp) :: energies(num_states)

        do i = 1, num_states
            call decode_bit_det(nI, ilut_list(:,i))
            if (tHPHF) then
                energies(i) = hphf_diag_helement(nI, ilut_list(:,i))
            else
                energies(i) = get_helement(nI, nI, 0)
            end if
        end do

        ! Sort the states in order of energy, from smallest to largest.
        call sort(energies, ilut_list(0:NIfTot, 1:num_states))

        ! If requested, sort singles and doubles to the top, keeping the rest of
        ! the ordering the same.
        if (tSortDoubles) then
            num_sing_doub = 0
            block_size = 0
            do i = 1, num_states
                excit_level = FindBitExcitLevel(ilut_list(:,i), ilutHF)
                if (excit_level <= 2) then
                    num_sing_doub = num_sing_doub + 1
                    temp_ilut = ilut_list(:,i)

                    ! Move block of non-singles or doubles down one in ilut_list.
                    ilut_list(:, num_sing_doub+1:num_sing_doub+block_size) = &
                        ilut_list(:, num_sing_doub:num_sing_doub+block_size-1)

                    ! Then move the single or double just found into the space
                    ! created above this block.
                    ilut_list(:, num_sing_doub) = temp_ilut
                else
                    block_size = block_size + 1
                end if
            end do
        end if

    end subroutine sort_states_by_energy

    subroutine find_determ_states_and_sort(ilut_list, ilut_list_size)

        integer, intent(in) :: ilut_list_size
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, ilut_list_size)
        integer :: i
        logical :: in_determ_space

        do i = 1, ilut_list_size
            call check_if_in_determ_space(ilut_list(0:NIfTot, i), in_determ_space)
            if (in_determ_space) call set_flag(ilut_list(0:NIfTot, i), flag_deterministic)
        end do

        call sort(ilut_list, ilut_lt, ilut_gt)

    end subroutine find_determ_states_and_sort

    subroutine sort_space_by_proc(ilut_list, ilut_list_size, num_states_procs)

        ! And also output the number of states on each processor in the space.

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

end module semi_stoch_procs
