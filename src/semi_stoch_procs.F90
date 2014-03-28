! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

! Some general procedures, created for the semi-stochastic (and trial wavefunction) code.
! Some are just used for initialisation and others in the main FCIQMC loop.

module semi_stoch_procs

    use AnnihilationMod, only: RemoveDetHashIndex
    use bit_rep_data, only: flag_deterministic, nIfDBO, NIfD, NIfTot, test_flag, &
                            flag_is_initiator, NOffSgn, NIfSgn
    use bit_reps, only: decode_bit_det, set_flag, extract_part_sign, extract_sign, &
                        encode_sign
    use CalcData
    use constants
    use davidson_neci, only: davidson_eigenvector, parallel_sparse_hamil_type, perform_davidson
    use DetBitOps, only: ilut_lt, ilut_gt, FindBitExcitLevel, DetBitLT, &
                         count_set_bits, DetBitEq, sign_lt, sign_gt, IsAllowedHPHF, &
                         EncodeBitDet
    use Determinants, only: get_helement, GetH0Element3, GetH0Element4
    use FciMCData, only: ilutHF, Hii, CurrentH, determ_proc_sizes, determ_proc_indices, &
                         full_determ_vector, partial_determ_vector, core_hamiltonian, &
                         determ_space_size, SpawnedParts, SemiStoch_Comms_Time, &
                         SemiStoch_Multiply_Time, TotWalkers, CurrentDets, CoreTag, &
                         PDetermTag, FDetermTag, IDetermTag, indices_of_determ_states, &
                         HashIndex, core_space, CoreSpaceTag, ll_node, nWalkerHashes, &
                         tFill_RDM, IterLastRDMFill, full_determ_vector_av, &
                         tFillingStochRDMonFly, Iter, IterRDMStart, CoreHashIndex, &
                         core_ham_diag, DavidsonTag, Fii, HFDet
    use hash, only: DetermineDetNode, FindWalkerHash
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use MI_integrals, only: MI_off_diag_helement
    use MomInv, only: IsAllowedMI
    use nElRDMMod, only: fill_RDM_offdiag_deterministic
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIBarrier, MPIArg, &
                             MPIAllGatherV, MPISum, MPISumAll, MPIScatterV
    use ParallelHelper, only: root
    use ras, only: core_ras
    use searching, only: BinSearchParts
    use sort_mod, only: sort
    use sparse_arrays, only: sparse_core_ham, SparseCoreHamilTags, deallocate_sparse_ham, &
                            core_connections, sparse_ham, hamil_diag, HDiagTag, &
                            SparseHamilTags, allocate_sparse_ham_row
    use SystemData, only: nel, tHPHF, nBasis, BRR, ARR, tUEG, tMomInv
    use timing_neci
    use util_mod, only: get_free_unit

    implicit none

contains

    subroutine deterministic_projection()

        ! This subroutine gathers together the partial_determ_vectors from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        integer :: i, j, info, ierr

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_determ_vector, full_determ_vector, determ_proc_sizes, &
                            determ_proc_indices)

        call halt_timer(SemiStoch_Comms_Time)

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Multiply_Time)

        if(tFillingStochRDMonFly) then !Update the average signs in full_determ_vector_av
            full_determ_vector_av(:)=(((real(Iter,dp)-IterRDMStart)*full_determ_vector_av(:)) &
                                      + full_determ_vector(:))/(real(Iter,dp) - IterRDMStart + 1.0_dp)
        endif
            
        if (determ_proc_sizes(iProcIndex) >= 1) then

            ! Perform the multiplication. This can be done in two ways depending on
            ! whether the the core Hamiltonian uses a sparse representation or not.
            if (tSparseCoreHamil) then
                
                if(tFill_RDM) call fill_RDM_offdiag_deterministic()

                    !For the moment, we're only adding in these contributions when we need the energy
                    !This will need refinement if we want to continue with the option of inst vs true full RDMs
                    ! (as in another CMO branch).

                partial_determ_vector = 0.0_dp

                do i = 1, determ_proc_sizes(iProcIndex)
                    do j = 1, sparse_core_ham(i)%num_elements
                        partial_determ_vector(i) = partial_determ_vector(i) - &
                            sparse_core_ham(i)%elements(j)*full_determ_vector(sparse_core_ham(i)%positions(j))
                    end do
                end do

            else

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

            end if

            ! Now add shift*full_determ_vector, to account for the shift, not stored in
            ! core_hamiltonian.
            partial_determ_vector = partial_determ_vector + &
               DiagSft * full_determ_vector(determ_proc_indices(iProcIndex)+1:&
                 determ_proc_indices(iProcIndex)+determ_proc_sizes(iProcIndex))

            ! Now multiply the vector by tau to get the final projected vector.
            partial_determ_vector = partial_determ_vector * tau

        end if

        call halt_timer(SemiStoch_Multiply_Time)

    end subroutine deterministic_projection

    function is_core_state(ilut) result (core_state)

        use FciMCData, only: ll_node

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: DetHash, PartInd
        type(ll_node), pointer :: temp_node
        logical :: core_state

        core_state = .false.

        call decode_bit_det(nI, ilut)
        DetHash = FindWalkerHash(nI, int(determ_space_size,sizeof_int))
        temp_node => CoreHashIndex(DetHash)

        if (temp_node%ind == 0) then
            ! If there are no core states with this hash value.
            nullify(temp_node)
            return
        else
            do while (associated(temp_node))
                if (DetBitEQ(ilut, core_space(:,temp_node%ind),NIfDBO)) then
                    core_state = .true.
                    nullify(temp_node)
                    return
                end if
                temp_node => temp_node%next
            end do
        end if

    end function is_core_state

    subroutine calc_determ_hamil_normal()

        integer :: i, j, iproc, col_index, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer(TagIntType) :: TempStoreTag
        character(len=*), parameter :: t_r = "calc_determ_hamil_normal"

        ! Allocate the core hamiltonian.
        allocate(core_hamiltonian(determ_proc_sizes(iProcIndex), determ_space_size), stat=ierr)
        call LogMemAlloc('core_hamiltonian', int(determ_space_size*&
                         &determ_proc_sizes(iProcIndex),sizeof_int), 8, t_r, CoreTag, ierr)
        allocate(core_ham_diag(determ_proc_sizes(iProcIndex)), stat=ierr)

        ! temp_store is storage space for bitstrings so that the Hamiltonian matrix
        ! elements can be calculated.
        allocate(temp_store(0:NIfTot, maxval(determ_proc_sizes)), stat=ierr)
        call LogMemAlloc('temp_store', maxval(determ_proc_sizes)*(NIfTot+1), 8, t_r, &
                         TempStoreTag, ierr)

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
            ! to these two processors.
            do i = 1, determ_proc_sizes(iProcIndex)

                call decode_bit_det(nI, SpawnedParts(:, i))

                do j = 1, determ_proc_sizes(iproc)

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
                        core_ham_diag(i) = core_hamiltonian(i, col_index + j)
                        ! We calculate and store CurrentH at this point for ease.
                        if ((.not. tRegenDiagHEls) .and. (.not. tReadPops)) &
                            CurrentH(1,i) = core_hamiltonian(i, col_index + j)
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
        call LogMemDealloc(t_r, TempStoreTag, ierr)

    end subroutine calc_determ_hamil_normal

    subroutine recalc_core_hamil_diag(old_Hii, new_Hii)

        real(dp) :: old_Hii, new_Hii
        real(dp) :: Hii_shift
        integer :: i, j

        write(6,'(a56)') "Recalculating diagonal elements of the core Hamiltonian."

        Hii_shift = old_Hii - new_Hii

        if (tSparseCoreHamil) then
            do i = 1, determ_proc_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    if (sparse_core_ham(i)%positions(j) == i + determ_proc_indices(iProcIndex)) then
                        sparse_core_ham(i)%elements(j) = sparse_core_ham(i)%elements(j) + Hii_shift
                    end if
                end do
            end do
        else
        end if

        core_ham_diag = core_ham_diag + Hii_shift

    end subroutine recalc_core_hamil_diag

    subroutine generate_core_connections()

        if (tSparseCoreHamil) then
            call gen_core_connections_sparse()
        else
            ! To be implemented.
            ! call generate_core_connections_full()
        end if

    end subroutine generate_core_connections

    subroutine gen_core_connections_sparse()

        integer :: i, j, ic, counter, ierr
        integer :: Ex(2,nel)
        logical :: tSign
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer(TagIntType) :: TempStoreTag
        character(len=*), parameter :: t_r = "calculate_det_hamiltonian_sparse"

        integer :: nI(nel), nJ(nel)

        allocate(core_connections(determ_proc_sizes(iProcIndex)))

        allocate(temp_store(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', maxval(determ_proc_sizes)*(NIfTot+1), 8, t_r, &
                         TempStoreTag, ierr)

        ! Stick together the deterministic states from all processors, on all processors.
        call MPIAllGatherV(SpawnedParts(:,1:determ_proc_sizes(iProcIndex)), temp_store, &
                       determ_proc_sizes, determ_proc_indices)

        ! Over all core states on this processor.
        do i = 1, determ_proc_sizes(iProcIndex)

            ! The number of non-zero elements in this array will be almost the same as in
            ! the core Hamiltonian array, except the diagonal element is not considered,
            ! so there will actually be one less.
            allocate(core_connections(i)%elements(sparse_core_ham(i)%num_elements-1)) 
            allocate(core_connections(i)%positions(sparse_core_ham(i)%num_elements-1))

            ! The total number of non-zero elements in row i.
            core_connections(i)%num_elements = sparse_core_ham(i)%num_elements-1

            counter = 0
            do j = 1, sparse_core_ham(i)%num_elements
                ! If not the diagonal element.
                if (sparse_core_ham(i)%positions(j) /= i + determ_proc_indices(iProcIndex)) then
                    Ex = 0
                    Ex(1,1) = nel
                    counter = counter + 1
                    ! The positions of the non-zero and non-diagonal elements in this row i.
                    core_connections(i)%positions(counter) = sparse_core_ham(i)%positions(j)

                    ic = FindBitExcitLevel(SpawnedParts(:,i), temp_store(:, sparse_core_ham(i)%positions(j)))
                    call GetBitExcitation(SpawnedParts(0:NIfD,i), temp_store(0:NIfD, &
                                          sparse_core_ham(i)%positions(j)),Ex,tSign)
                    if (tSign) then
                        ! Odd number of permutations. Minus the excitation level.
                        core_connections(i)%elements(counter) = -ic
                    else
                        ! Even number of permutations. The excitation level.
                        core_connections(i)%elements(counter) = ic
                    end if
                end if
            end do

        end do

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)

    end subroutine gen_core_connections_sparse

    subroutine store_whole_core_space()

        integer :: ierr
        character(len=*), parameter :: t_r = "store_whole_core_space"

        allocate(core_space(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('core_space', maxval(determ_proc_sizes)*(NIfTot+1), 8, t_r, &
                         CoreSpaceTag, ierr)

        call MPIAllGatherV(SpawnedParts(:,1:determ_proc_sizes(iProcIndex)), core_space, &
                       determ_proc_sizes, determ_proc_indices)

    end subroutine store_whole_core_space

    subroutine initialise_core_hash_table()

        use bit_reps, only: decode_bit_det
        use FciMCData, only: core_space
        use SystemData, only: nel

        integer :: nI(nel)
        integer :: i, ierr, DetHash, counter, total
        type(ll_node), pointer :: temp_node
        character(len=*), parameter :: t_r = "initialise_core_hash_table"

        allocate(CoreHashIndex(determ_space_size), stat=ierr)

        do i = 1, determ_space_size
            CoreHashIndex(i)%ind = 0
            nullify(CoreHashIndex(i)%next)
        end do

        do i = 1, determ_space_size
            call decode_bit_det(nI, core_space(:,i))
            DetHash = FindWalkerHash(nI, int(determ_space_size,sizeof_int))
            temp_node => CoreHashIndex(DetHash)
            ! If the first element in the list has not been used.
            if (temp_node%ind == 0) then
                temp_node%ind = i
            else
                do while (associated(temp_node%next))
                    temp_node => temp_node%next
                end do
                allocate(temp_node%next)
                nullify(temp_node%next%next)
                temp_node%next%ind = i
            end if
        end do

        nullify(temp_node)

    end subroutine initialise_core_hash_table

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
            tot_num_states = int(num_states,MPIArg)
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
                states_rmvd_all_procs = int(states_rmvd_this_proc,MPIArg)
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

    subroutine sort_space_by_proc(ilut_list, ilut_list_size, num_states_procs)

        ! And also output the number of states on each processor in the space.

        integer, intent(in) :: ilut_list_size
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:ilut_list_size)
        integer(MPIArg), intent(out) :: num_states_procs(0:nProcessors-1)
        integer(n_int), allocatable, dimension(:,:) :: temp_list
        integer, allocatable, dimension(:) :: proc_list
        integer :: nI(nel)
        integer :: i, ierr
        integer :: counter(0:nProcessors-1)
        integer(TagIntType) :: TempConTag, ProcListTag
        character(len=*), parameter :: t_r = "sort_space_by_proc"

        allocate(proc_list(ilut_list_size), stat=ierr)
        call LogMemAlloc('proc_list', ilut_list_size, sizeof_int, t_r, ProcListTag, ierr)

        allocate(temp_list(0:NIfTot, ilut_list_size), stat=ierr)
        call LogMemAlloc('con_space_temp', ilut_list_size*(NIfTot+1), size_n_int, t_r, &
                         TempConTag, ierr)

        num_states_procs = 0

        ! Create a list, proc_list, with the processor numbers of the corresponding iluts.
        do i = 1, ilut_list_size
            call decode_bit_det(nI, ilut_list(:,i))
            proc_list(i) = DetermineDetNode(nI,0)
            num_states_procs(proc_list(i)) = int(num_states_procs(proc_list(i)) + 1,MPIArg)
        end do

        counter(0) = 0
        do i = 1, nProcessors-1
            counter(i) = sum(num_states_procs(:i-1))
        end do

        do i = 1, ilut_list_size
            counter(proc_list(i)) = counter(proc_list(i)) + 1 
            temp_list(:, counter(proc_list(i))) = ilut_list(:,i)
        end do

        ilut_list = temp_list

        deallocate(temp_list, stat=ierr)
        deallocate(proc_list, stat=ierr)
        call LogMemDealloc(t_r, TempConTag, ierr)
        call LogMemDealloc(t_r, ProcListTag, ierr)

    end subroutine sort_space_by_proc

    subroutine fill_in_CurrentH()

        integer(int64) :: i
        integer :: nI(nel)

        CurrentH = 0

        do i = 1, TotWalkers

            call decode_bit_det(nI, CurrentDets(:,i))

            if (tHPHF) then
                CurrentH(1,i) = hphf_diag_helement(nI, CurrentDets(:,i)) - Hii
            else
                CurrentH(1,i) = get_helement(nI, nI, 0) - Hii
            end if

        end do

    end subroutine fill_in_CurrentH

    subroutine write_core_space()

        integer :: i, k, iunit, ierr
        integer(int64) :: j
        logical :: texist
        character(len=*), parameter :: t_r='write_core_space'

        write(6,'(a35)') "Writing the core space to a file..."

        iunit = get_free_unit()

        ! Let each processor write its core states to the file. Each processor waits for
        ! the processor before it to finish before starting.
        do i = 0, nProcessors-1

            if (iProcIndex == i) then

                if (i == 0) then
                    open(iunit, file='CORESPACE', status='replace')
                else
                    inquire(file='CORESPACE',exist=texist)
                    if(.not.texist) call stop_all(t_r,'"CORESPACE" file cannot be found')
                    open(iunit, file='CORESPACE', status='old', position='append')
                end if
                
                do j = 1, TotWalkers 
                    if (test_flag(CurrentDets(:,j), flag_deterministic)) then
                        do k = 0, NIfDBO
                            write(iunit, '(i24)', advance='no') CurrentDets(k,j)
                        end do
                        write(iunit, *)
                    end if
                end do

                close(iunit)

            end if

            call MPIBarrier(ierr)

        end do

    end subroutine write_core_space

    subroutine add_core_states_currentdets()

        ! And if the state is already present, simply set its flag.
        ! Also sort the states afterwards.

        integer :: i, comp, MinInd, PartInd, nwalkers
        logical :: tSuccess

        integer :: j

        MinInd = 1
        nwalkers = int(TotWalkers,sizeof_int)

        do i = 1, determ_proc_sizes(iProcIndex)

            if (nwalkers > 0) then
                ! If there is only one state in CurrentDets to check then BinSearchParts doesn't
                ! return the desired value for PartInd, so do this separately...
                if (MinInd == nwalkers) then
                    comp = DetBitLT(CurrentDets(:,MinInd), SpawnedParts(:,i), NIfDBO, .false.)
                    if (comp == 0) then
                        tSuccess = .true.
                        PartInd = MinInd
                    else if (comp == 1) then
                        tSuccess = .false.
                        PartInd = MinInd
                    else if (comp == -1) then
                        tSuccess = .false.
                        PartInd = MinInd - 1
                    end if
                else
                    call BinSearchParts(SpawnedParts(:,i), MinInd, nwalkers, PartInd, tSuccess)
                end if
            else
                tSuccess = .false.
                PartInd = 0
            end if

            if (tSuccess) then
                call set_flag(CurrentDets(:,PartInd), flag_deterministic)
                if (tTruncInitiator) then
                    call set_flag(CurrentDets(:,PartInd), flag_is_initiator(1))
                    call set_flag(CurrentDets(:,PartInd), flag_is_initiator(2))
                end if
                MinInd = PartInd
            else
                ! Move all states below PartInd down one and insert the new state in the slot.
                CurrentDets(:, PartInd+2:nwalkers+1) = CurrentDets(:, PartInd+1:nwalkers)
                CurrentDets(:, PartInd+1) = SpawnedParts(:,i)
                nwalkers = nwalkers + 1
                MinInd = PartInd + 1
            end if

        end do

        call sort(CurrentDets(:,1:nwalkers), ilut_lt, ilut_gt)

        TotWalkers = int(nwalkers, int64)

    end subroutine add_core_states_currentdets

    subroutine add_core_states_currentdet_hash

        ! This routine adds the core states in SpawnedParts into CurrentDets. For all
        ! such states already in CurrentDets, we want to keep the amplitude (which
        ! may have come from a popsfile).

        ! This routine is for when the tHashWalkerList option is used. In this case,
        ! as all core states are always kept in the list, it is beneficial to keep
        ! them at the top always. So, in this routine, we move the non-core states
        ! in CurrentDets to the end and add the new core states in the gaps.

        integer :: i, DetHash, PartInd, nwalkers, i_non_core
        integer :: nI(nel)
        real(dp) :: walker_sign(lenof_sign)
        type(ll_node), pointer :: temp_node, curr, prev
        logical :: tSuccess

        nwalkers = int(TotWalkers,sizeof_int)

        ! First find which CurrentDet states are in the core space.
        do i = 1, determ_proc_sizes(iProcIndex)

            tSuccess = .false.
            call decode_bit_det (nI, SpawnedParts(:,i))
            DetHash = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(DetHash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if (DetBitEQ(SpawnedParts(:,i), CurrentDets(:,temp_node%ind),NIfDBO)) then
                        tSuccess = .true.
                        PartInd = temp_node%ind
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            nullify(temp_node)

            ! Core state i is in CurrentDets.
            if (tSuccess) then
                call set_flag(CurrentDets(:,PartInd), flag_deterministic)
                ! Copy the amplitude of the state across to SpawnedParts.
                call extract_sign(CurrentDets(:,PartInd), walker_sign)
                call encode_sign(SpawnedParts(:,i), walker_sign)
            else
                ! This will be a new state added to CurrentDets.
                nwalkers = nwalkers + 1
            end if

        end do

        ! Next loop through CurrentDets and move all non-core states to after the last
        ! core state slot in SpawnedParts.
        i_non_core = determ_proc_sizes(iProcIndex)
        do i = 1, int(TotWalkers,sizeof_int)
            if (.not. test_flag(CurrentDets(:,i), flag_deterministic)) then
                i_non_core = i_non_core + 1
                SpawnedParts(:,i_non_core) = CurrentDets(:,i)
            end if
        end do

        ! Now copy all the core states in SpawnedParts into CurrentDets.
        ! Note that the amplitude in CurrentDets was copied across, so this is fine.
        do i = 1, nwalkers
            CurrentDets(:,i) = SpawnedParts(:,i)
        end do

        ! Reset the hash index array.
        do i = 1, nWalkerHashes
            curr => HashIndex(i)%next
            prev => HashIndex(i)
            prev%ind = 0
            nullify(prev%next)
            do while (associated(curr))
                prev => curr
                curr => curr%next
                deallocate(prev)
            end do
        end do
        nullify(curr)
        nullify(prev)

        ! Finally, add the indices back into the hash index array.
        do i = 1, nwalkers
            call decode_bit_det(nI, CurrentDets(:,i))
            DetHash = FindWalkerHash(nI,nWalkerHashes)
            temp_node => HashIndex(DetHash)
            ! If the first element in the list has not been used.
            if (temp_node%ind == 0) then
                temp_node%ind = i
            else
                do while (associated(temp_node%next))
                    temp_node => temp_node%next
                end do
                allocate(temp_node%next)
                nullify(temp_node%next%next)
                temp_node%next%ind = i
            end if
            nullify(temp_node)

            ! These core states will always stay in the same position.
            if (i <= determ_proc_sizes(iProcIndex)) indices_of_determ_states(i) = i
        end do

        TotWalkers = int(nwalkers, int64)

    end subroutine add_core_states_currentdet_hash

    subroutine return_most_populated_states(n_keep, largest_walkers, norm)

        ! Return the most populated states in CurrentDets on *this* processor only. 
        ! Also return the norm of these states, if requested.

        integer, intent(in) :: n_keep
        integer(n_int), intent(out) :: largest_walkers(0:NIfTot, n_keep)
        real(dp), intent(out), optional :: norm
        integer :: i, j, smallest_pos
        real(dp) :: smallest_sign, sign_curr_real
        real(dp), dimension(lenof_sign) :: sign_curr, low_sign

        largest_walkers = 0
        smallest_sign = 0.0_dp
        smallest_pos = 1
        if (present(norm)) norm = 0.0_dp

        ! Run through all walkers on process.
        do i = 1, int(TotWalkers,sizeof_int)
            call extract_sign(CurrentDets(:,i), sign_curr)

            if (lenof_sign == 1) then
                sign_curr_real = real(abs(sign_curr(1)),dp)
            else
                sign_curr_real = sqrt(real(sign_curr(1),dp)**2 + real(sign_curr(lenof_sign),dp)**2)
            endif

            if (present(norm)) norm = norm + (sign_curr_real**2.0)

            ! Is this determinant more populated than the smallest. First in the list is always
            ! the smallest.
            if (sign_curr_real > smallest_sign) then
                largest_walkers(:,smallest_pos) = CurrentDets(:,i)

                ! Instead of resorting, just find new smallest sign and position.
                call extract_sign(largest_walkers(:,1),low_sign)

                if (lenof_sign == 1) then
                    smallest_sign = real(abs(low_sign(1)),dp)
                else
                    smallest_sign = sqrt(real(low_sign(1),dp)**2+real(low_sign(lenof_sign),dp)**2)
                endif

                smallest_pos = 1
                do j = 2, n_keep
                    if (smallest_sign < 1.0e-7_dp) exit
                    call extract_sign(largest_walkers(:,j), low_sign)
                    if (lenof_sign == 1) then
                        sign_curr_real = real(abs(low_sign(1)), dp)
                    else
                        sign_curr_real = sqrt(real(low_sign(1),dp)**2 + real(low_sign(lenof_sign),dp)**2)
                    end if

                    if (sign_curr_real < smallest_sign) then
                        smallest_pos = j
                        smallest_sign = sign_curr_real
                    end if
                end do

            endif

        end do

        call sort(largest_walkers(:,1:n_keep), sign_lt, sign_gt)

    end subroutine return_most_populated_states

    subroutine return_largest_indices(n_keep, list_size, list, largest_indices)

        ! Return the indices of the largest elements in list.

        integer, intent(in) :: n_keep, list_size
        real(dp), intent(in) :: list(list_size)
        integer, intent(out) :: largest_indices(n_keep)
        integer :: i, j, ind, smallest_pos
        real(dp) :: smallest_sign, sign_curr, sign_curr_abs, low_sign

        largest_indices = 0
        smallest_sign = 0.0_dp
        smallest_pos = 1

        do i = 1, list_size
            sign_curr = list(i)
            sign_curr_abs = abs(sign_curr)

            if (sign_curr_abs > smallest_sign) then
                largest_indices(smallest_pos) = i

                low_sign = list(largest_indices(1))

                smallest_sign = abs(low_sign)

                smallest_pos = 1
                do j = 2, n_keep
                    if (smallest_sign < 1.0e-7_dp) exit
                    ind = largest_indices(j)
                    if (ind == 0) then
                        low_sign = 0.0_dp
                    else
                        low_sign = list(ind)
                    end if

                    sign_curr_abs = abs(low_sign)

                    if (sign_curr_abs < smallest_sign) then
                        smallest_pos = j
                        smallest_sign = sign_curr_abs
                    end if
                end do
            endif
        end do

    end subroutine return_largest_indices

    subroutine start_walkers_from_core_ground()

        integer :: i, counter, ierr
        real(dp) :: eigenvec_pop
        character(len=*), parameter :: t_r = "start_walkers_from_core_ground"

        ! Create the arrays used by the Davidson routine.
        ! First, the whole Hamiltonian in sparse form.
        allocate(sparse_ham(determ_proc_sizes(iProcIndex)))
        allocate(SparseHamilTags(2, determ_proc_sizes(iProcIndex)))
        do i = 1, determ_proc_sizes(iProcIndex)
            call allocate_sparse_ham_row(sparse_ham, i, sparse_core_ham(i)%num_elements, "sparse_ham", SparseHamilTags(:,i)) 
            sparse_ham(i)%elements = sparse_core_ham(i)%elements
            sparse_ham(i)%positions = sparse_core_ham(i)%positions
            sparse_ham(i)%num_elements = sparse_core_ham(i)%num_elements
        end do

        ! Next create the diagonal used by Davidson by copying the core one.
        allocate(hamil_diag(determ_proc_sizes(iProcIndex)),stat=ierr)
        call LogMemAlloc('hamil_diag', int(determ_proc_sizes(iProcIndex),sizeof_int), 8, t_r, HDiagTag, ierr)
        hamil_diag = core_ham_diag

        write(6,'(a69)') "Using the deterministic ground state as initial walker configuration."
        write(6,'(a34)') "Performing Davidson calculation..."
        call neci_flush(6)

        ! Call the Davidson routine to find the ground state of the core space. 
        call perform_davidson(parallel_sparse_hamil_type, .false.)

        write(6,'(a30)') "Davidson calculation complete."
        call neci_flush(6)

        ! The ground state compnents are now stored in davidson_eigenvector on the root.
        ! First, we need to normalise this vector to have the correct 'number of walkers'.
        if (iProcIndex == root) then
            eigenvec_pop = 0.0_dp
            do i = 1, determ_space_size
                eigenvec_pop = eigenvec_pop + abs(davidson_eigenvector(i))
            end do
            if (tStartSinglePart) then
                davidson_eigenvector = davidson_eigenvector*InitialPart/eigenvec_pop
            else
                davidson_eigenvector = davidson_eigenvector*InitWalkers/eigenvec_pop
            end if
        end if

        ! Send the components to the correct processors and use partial_determ_vector as
        ! temporary space.
        call MPIScatterV(davidson_eigenvector, determ_proc_sizes, determ_proc_indices, &
                         partial_determ_vector, determ_proc_sizes(iProcIndex), ierr)

        ! Finally, copy these amplitudes across to the corresponding states in CurrentDets.
        counter = 0
        do i = 1, int(TotWalkers)
            if (test_flag(CurrentDets(:,i), flag_deterministic)) then
                counter = counter + 1
                call encode_sign(CurrentDets(:,i), partial_determ_vector(counter))
            end if
        end do

        partial_determ_vector = 0.0_dp
        deallocate(davidson_eigenvector)
        call LogMemDealloc(t_r, DavidsonTag, ierr)
        deallocate(hamil_diag, stat=ierr)
        call LogMemDealloc(t_r, HDiagTag, ierr)
        call deallocate_sparse_ham(sparse_ham, 'sparse_ham', SparseHamilTags)

    end subroutine start_walkers_from_core_ground

    subroutine return_mp1_amp_and_mp2_energy(nI, ilut, ex, tParity, amp, energy_contrib)

        ! For a given determinant (input as nI), find the amplitude of it in the MP1 wavefunction.
        ! Also return the contribution from this determinant in the MP2 energy.

        ! To use this routine, generate an excitation from the Hartree-Fock determinant using the
        ! GenExcitations3 routine. This will return nI, ex and tParity which can be input into this
        ! routine.

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(out) :: amp, energy_contrib
        integer :: ic
        real(dp) :: hel, H0tmp, denom

        amp = 0.0_dp
        energy_contrib = 0.0_dp

        if (ex(1,2) == 0) then
            ic = 1
        else
            ic = 2
        end if

        if (tHPHF) then
            ! Assume since we are using HPHF that the alpha and
            ! beta orbitals of the same spatial orbital have the same
            ! fock energies, so can consider either.
            hel = hphf_off_diag_helement(HFDet, nI, iLutHF, ilut)
        else if (tMomInv) then
            hel = MI_off_diag_helement(HFDet, nI, iLutHF, ilut)
        else
            hel = get_helement(HFDet, nI, ic, ex, tParity)
        end if

        if (tUEG) then
            ! This will calculate the MP2 energies without having to use the fock eigenvalues.
            ! This is done via the diagonal determinant hamiltonian energies.
            H0tmp = getH0Element4(nI, HFDet)
        else
            H0tmp = getH0Element3(nI)
        end if

        ! If the relevant excitation from the Hartree-Fock takes electrons from orbitals
        ! (i,j) to (a,b), then denom will be equal to 
        ! \epsilon_a + \epsilon_b - \epsilon_i - \epsilon_j
        ! as required in the denominator of the MP1 amplitude and MP2 energy.
        denom = Fii - H0tmp

        if (.not. (abs(denom) > 0.0_dp)) then
            call warning_neci("return_mp1_amp_and_mp2_energy", &
            "One of the determinants under consideration for the MP1 wave function is degenerate &
            &with the Hartree-Fock determinant. Degenerate perturbation theory has not been &
            &considered, but the amplitude of this determinant will be returned as huge(0.0_dp) &
            &so that it should be included in the space.")
            amp = huge(0.0_dp)
        else
            amp = hel/denom
            energy_contrib = (hel**2)/denom
        end if

    end subroutine return_mp1_amp_and_mp2_energy

    subroutine end_semistoch()

        character(len=*), parameter :: t_r = "end_semistoch"
        integer :: ierr

        if (tSparseCoreHamil) call deallocate_sparse_ham(sparse_core_ham, 'sparse_core_ham', SparseCoreHamilTags)

        if (allocated(core_hamiltonian)) then
            deallocate(core_hamiltonian, stat=ierr)
            call LogMemDealloc(t_r, CoreTag, ierr)
        end if
        if (allocated(full_determ_vector)) then
            deallocate(full_determ_vector, stat=ierr)
            call LogMemDealloc(t_r, FDetermTag, ierr)
        end if
        if (allocated(partial_determ_vector)) then
            deallocate(partial_determ_vector, stat=ierr)
            call LogMemDealloc(t_r, PDetermTag, ierr)
        end if
        if (allocated(indices_of_determ_states)) then
            deallocate(indices_of_determ_states, stat=ierr)
            call LogMemDealloc(t_r, IDetermTag, ierr)
        end if

    end subroutine end_semistoch
    
end module semi_stoch_procs
