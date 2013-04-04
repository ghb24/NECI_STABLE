#include "macros.h"

module semi_stochastic

    use bit_rep_data, only: flag_deterministic, nIfDBO, nOffY, nIfY, nOffFlag, &
                            flag_determ_parent, deterministic_mask, determ_parent_mask, &
                            flag_is_initiator, flag_bit_offset, NOffSgn
    use bit_reps, only: NIfD, NIfTot, decode_bit_det, encode_bit_rep, set_flag
    use CalcData, only: tRegenDiagHEls, tau, tSortDetermToTop, tTruncInitiator, DiagSft, &
                        tStartCAS, tReadPops
    use csf, only: csf_get_yamas, get_num_csfs, get_csf_bit_yama, csf_apply_yama
    use csf_data, only: iscsf, csf_orbital_mask
    use constants
    use DetBitOps, only: ilut_lt, ilut_gt, count_open_orbs, FindBitExcitLevel, DetBitLT, &
                         count_set_bits, IsAllowedHPHF, DetBitEq, TestClosedShellDet
    use DeterminantData, only: write_det
    use Determinants, only: get_helement
    use enumerate_excitations
    use FciMCData, only: HFDet, ilutHF, iHFProc, Hii, CurrentDets, CurrentH, &
                         determ_proc_sizes, determ_proc_indices, &
                         full_determ_vector, partial_determ_vector, core_hamiltonian, &
                         determ_space_size, TotWalkers, TotWalkersOld, &
                         indices_of_determ_states, ProjEDet, SpawnedParts, SparseHamilTags, &
                         HDiagTag, CoreTag, FDetermTag, PDetermTag, &
                         trial_space, trial_space_size, SemiStoch_Comms_Time, &
                         SemiStoch_Multiply_Time
    use gndts_mod, only: gndts
    use hash , only: DetermineDetNode
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIBarrier, MPIArg, &
                             MPIAllGatherV, MPIAllGather, MPIScatter, MPIScatterV
    use ParallelHelper, only: root
    use sort_mod, only: sort
    use sym_mod, only: getsym
    use SystemData
    use timing_neci

    implicit none

    ! Some of the following subroutines are used for generating both the trial space and the
    ! deterministic space, but must be used slightly differently for each case. Hence, the
    ! following integers may be input to these subroutines to tell them what to do.
    integer :: called_from_semistoch = 1
    integer :: called_from_trial = 2

contains

    subroutine init_semi_stochastic()

        ! Initialise the semi-stochastic information. This includes enumerating a list
        ! of all determinants or CSFs in the deterministic space and calculating and
        ! storing the resulting Hamiltonian matrix elements. The lists which will store
        ! the psip vectors in the deterministic space are also allocated.

        integer :: i, j, IC, ierr, comp
        integer :: nI(nel)
        integer :: space_size
        character (len=*), parameter :: this_routine = "init_semi_stochastic"

        write(6,'()')
        write(6,'(a56)') "============ Semi-stochastic initialisation ============"
        call neci_flush(6)

        allocate(determ_proc_sizes(0:nProcessors-1))
        allocate(determ_proc_indices(0:nProcessors-1))
        determ_proc_sizes = 0
        determ_proc_indices = 0

        ! Count the HF state and add it to SpawnedParts.
        determ_proc_sizes(iHFProc) = 1
        if (.not. tReadPops) then
            SpawnedParts(:,1) = CurrentDets(:,1)
        else
            if (iProcIndex == iHFProc) then
                SpawnedParts(:,1) = 0
                SpawnedParts(0:NIfTot,1) = ilutHF(0:NIfTot)
                call set_flag(SpawnedParts(:,1), flag_deterministic)
            end if
        end if

        if (.not. (tDeterminantCore .or. tCSFCore)) then
            call warning_neci("init_semi_stochastic", "You have not selected to use either &
                              &determinants or CSFs for the deterministic space. Determinants &
                              &will be used.")
            tDeterminantCore = .true.
        end if

        ! The following subroutines call the enumerating subroutines to create all excitations
        ! and add these states to the main list, CurrentDets, on the correct processor. As
        ! they do this, they count the size of the deterministic space on each processor.
        write(6,'(a37)') "Generating the deterministic space..."
        call neci_flush(6)
        if (tStartCAS) then
            do i = 1, TotWalkers
                call set_flag(CurrentDets(:, i), flag_deterministic)
                call set_flag(CurrentDets(:, i), flag_is_initiator(1))
                call set_flag(CurrentDets(:, i), flag_is_initiator(2))
            end do
            call MPIAllGather(int(TotWalkers, MPIArg), determ_proc_sizes, ierr)
        else if (tDeterminantCore) then
            if (tDoublesCore) then
                call generate_sing_doub_determinants(called_from_semistoch)
            else if (tCASCore) then
                call generate_cas(called_from_semistoch)
            else if (tOptimisedCore) then
                call generate_optimised_core(called_from_semistoch)
            else if (tLowECore) then
                call generate_low_energy_core(called_from_semistoch)
            end if
        else if (tCSFCore) then
            if (tDoublesCore) then
                call generate_sing_doub_csfs(called_from_semistoch)
            else if (tCASCore) then
                call stop_all("init_semi_stochastic", "CAS core space with CSFs is not &
                              &currently implemented.")
            else if (tOptimisedCore) then
                call stop_all("init_semi_stochastic", "Optimised core space with CSFs is not &
                              &currently implemented.")
            else if (tLowECore) then
                call stop_all("init_semi_stochastic", "Low energy core space with CSFs is not &
                              &currently implemented.")
            end if
        end if

        if (tLimitDetermSpace) then
            tSortDetermToTop = .false.
            space_size = determ_proc_sizes(iProcIndex)
            call remove_high_energy_orbs(SpawnedParts(:, 1:space_size), space_size, &
                                           max_determ_size, .true.)
        else
            space_size = determ_proc_sizes(iProcIndex)
        end if
        call MPIAllGather(int(space_size, MPIArg), determ_proc_sizes, ierr)

        ! We now know the size of the deterministic space on each processor, so now find
        ! the total size of the space and also allocate vectors to store psip amplitudes
        ! and the deterministic Hamiltonian.
        determ_space_size = sum(determ_proc_sizes)
        allocate(full_determ_vector(determ_space_size), stat=ierr)
        call LogMemAlloc('full_determ_vector', int(determ_space_size,sizeof_int), 8, this_routine, &
                         FDetermTag, ierr)
        allocate(partial_determ_vector(determ_proc_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('partial_determ_vector', int(determ_proc_sizes(iProcIndex), &
                         sizeof_int), 8, this_routine, PDetermTag, ierr)
        allocate(core_hamiltonian(determ_proc_sizes(iProcIndex), determ_space_size), stat=ierr)
        call LogMemAlloc('core_hamiltonian', int(determ_space_size*&
                         &determ_proc_sizes(iProcIndex),sizeof_int), 8, this_routine, CoreTag, ierr)
        full_determ_vector = 0.0_dp
        partial_determ_vector = 0.0_dp

        ! Write space sizes to screen.
        write(6,'(a34,1X,i8)') "Total size of deterministic space:", determ_space_size
        write(6,'(a46,1X,i8)') "Size of deterministic space on this processor:", &
                    determ_proc_sizes(iProcIndex)
        call neci_flush(6)

        ! Also allocate the vector to store the positions of the deterministic states in
        ! CurrentDets.
        allocate(indices_of_determ_states(determ_proc_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('indices_of_determ_states', int(determ_proc_sizes(iProcIndex), &
                         sizeof_int), bytes_int, this_routine, FDetermTag, ierr)

        ! Calculate the indices in the full vector at which the various processors take
        ! over, relative to the first index position in the vector (ie, the first value
        ! in this vector will be 0, *not* 1 - this is required for mpiallgatherv later).
        do i = 0, nProcessors-1
            if (i == 0) then
                determ_proc_indices(i) = 0
            else
                determ_proc_indices(i) = determ_proc_indices(i-1) + &
                                                  determ_proc_sizes(i-1)
            end if
        end do

        ! Sort the list of basis states so that it is correctly ordered in the predefined
        ! order which is always kept throughout the simulation.
        call sort(SpawnedParts(:,1:determ_proc_sizes(iProcIndex)), ilut_lt, ilut_gt)

        ! Do a check that no states are in the deterministic space twice. The list is sorted
        ! already so simply check states next to each other in the list:
        do i = 2, determ_proc_sizes(iProcIndex)
            comp = DetBitLT(SpawnedParts(:, i-1), SpawnedParts(:, i), NIfD, .false.)
            if (comp == 0) then
                call decode_bit_det(nI, SpawnedParts(:,i))
                write(6,'(a18)') "State found twice:"
                write(6,*) SpawnedParts(:,i)
                call write_det(6, nI, .true.)
                call stop_all("init_semi_stochastic", &
                    "The same state has been found twice in the deterministic space.")
            end if
        end do

        ! Calculate and store the deterministic Hamiltonian.
        write(6,'(a56)') "Generating the Hamiltonian in the deterministic space..."
        call neci_flush(6)
        call calculate_det_hamiltonian_normal()

        ! If not reading from a popsfile, move the deterministic states to CurrentDets ready for
        ! FCIQMC spawning to start. If reading from a popsfile then CurrenDets will have been
        ! initialised already.
        if (.not. tReadPops) then
            CurrentDets(:,1:determ_proc_sizes(iProcIndex)) = &
                SpawnedParts(:,1:determ_proc_sizes(iProcIndex))
            SpawnedParts = 0
            TotWalkers = int(determ_proc_sizes(iProcIndex), int64)
            TotWalkersOld = int(determ_proc_sizes(iProcIndex), int64)
        else
            SpawnedParts = 0
        end if

        write(6,'(a40)') "Semi-stochastic initialisation complete."
        call neci_flush(6)

    end subroutine init_semi_stochastic

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

        ! Calcuation of the Hamiltonian matrix elements to be stored. Loop over each
        ! processor in order and broadcast the sorted vector of bitstrings from the
        ! processor to all other processors so that all matrix elements can be found.
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

    subroutine calculate_sparse_hamiltonian(num_states, ilut_list)

        use FciMCData, only: sparse_matrix_info, sparse_hamil, hamil_diag

        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        real(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character (len=*), parameter :: this_routine = "calculate_sparse_hamiltonian"

        ! Allocate arrays needed for the Hamiltonian construction and for finding the ground
        ! state
        allocate(sparse_hamil(num_states))
        allocate(hamiltonian_row(num_states), stat=ierr)
        call LogMemAlloc('hamiltonian_row', num_states, 8, this_routine, HRTag, ierr)
        allocate(hamil_diag(num_states), stat=ierr)
        call LogMemAlloc('hamil_diag', num_states, 8, this_routine, HDiagTag, ierr)
        allocate(sparse_row_sizes(num_states), stat=ierr)
        call LogMemAlloc('sparse_row_sizes', num_states, bytes_int, this_routine, SRTag, ierr)
        allocate(sparse_diag_positions(num_states), stat=ierr)
        call LogMemAlloc('sparse_diag_positions', num_states, bytes_int, this_routine, SDTag, ierr)
        allocate(indices(num_states), stat=ierr)
        call LogMemAlloc('indices', num_states, bytes_int, this_routine, ITag, ierr)

        ! Allocate the tags for the various rows of data in sparse_hamil.
        allocate(SparseHamilTags(2, num_states))

        ! Set each element to one to count the digonal elements straight away.
        sparse_row_sizes = 1

        ! In the following, the Hamiltonian for the new deterministic space is generated in a
        ! form which takes advantage of its sparse nature (see the type sparse_matrix_info in
        ! the FciMCData module).
        do i = 1, num_states

            hamiltonian_row = 0.0_dp

            ! Get nI form of the basis function.
            call decode_bit_det(nI, ilut_list(:, i))

            ! sparse_diag_positions(j) stores the number of non-zero elements in row j of the
            ! Hamiltonian, up to and including the diagonal element.
            ! sparse_row_sizes stores this number currently as all non-zero elements before
            ! the diagonal have been counted, as the Hamiltonian is symmetric.
            sparse_diag_positions(i) = sparse_row_sizes(i)

            do j = i, num_states

                call decode_bit_det(nJ, ilut_list(:, j))

                ! If on the diagonal of the Hamiltonian.
                if (i == j) then
                    if (.not. tHPHF) then
                        hamiltonian_row(j) = get_helement(nI, nJ, 0)
                    else
                        hamiltonian_row(j) = hphf_diag_helement(nI, ilut_list(:, i))
                    end if
                    hamil_diag(j) = hamiltonian_row(j)
                else
                    if (.not. tHPHF) then
                        hamiltonian_row(j) = get_helement(nI, nJ, ilut_list(:, i), &
                                                                 ilut_list(:, j))
                    else
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, ilut_list(:, i), &
                                                                           ilut_list(:, j))
                    end if
                    if (abs(hamiltonian_row(j)) > 0.0_dp) then
                        ! If element is nonzero, update the following sizes.
                        sparse_row_sizes(i) = sparse_row_sizes(i) + 1
                        sparse_row_sizes(j) = sparse_row_sizes(j) + 1
                    end if
                end if

            end do

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so
            ! allocate it.
            call allocate_sparse_hamil_row(i, sparse_row_sizes(i))

            sparse_hamil(i)%elements = 0.0_dp
            sparse_hamil(i)%positions = 0
            sparse_hamil(i)%num_elements = sparse_row_sizes(i)

            ! Now fill in all matrix elements beyond and including the diagonal, as these are
            ! stored in hamiltonian_row.
            counter = sparse_diag_positions(i)
            do j = i, num_states
                if (abs(hamiltonian_row(j)) > 0.0_dp) then
                    sparse_hamil(i)%positions(counter) = j
                    sparse_hamil(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
            end do

        end do

        ! At this point, sparse_hamil has been allocated with the correct sizes, but only the
        ! matrix elements in above and including the diagonal have been filled in. Now we must
        ! fill in the other elements. To do this, cycle through every element above the
        ! diagonal and fill in every corresponding element below it:
        indices = 1
        do i = 1, num_states
            do j = sparse_diag_positions(i) + 1, sparse_row_sizes(i)

                sparse_hamil(sparse_hamil(i)%positions(j))%&
                    &positions(indices(sparse_hamil(i)%positions(j))) = i
                sparse_hamil(sparse_hamil(i)%positions(j))%&
                    &elements(indices(sparse_hamil(i)%positions(j))) = sparse_hamil(i)%elements(j)

                indices(sparse_hamil(i)%positions(j)) = indices(sparse_hamil(i)%positions(j))+ 1

            end do
        end do

        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(this_routine, HRTag, ierr)
        deallocate(sparse_row_sizes, stat=ierr)
        call LogMemDealloc(this_routine, SRTag, ierr)
        deallocate(sparse_diag_positions, stat=ierr)
        call LogMemDealloc(this_routine, SDTag, ierr)
        deallocate(indices, stat=ierr)
        call LogMemDealloc(this_routine, ITag, ierr)

    end subroutine calculate_sparse_hamiltonian

    subroutine allocate_sparse_hamil_row(row, sparse_row_size)

        use FciMCData, only: sparse_hamil

        integer, intent(in) :: row, sparse_row_size
        integer :: ierr
        character (len=1024) :: string_row
        character (len=1024) :: var_name
        character (len=*), parameter :: this_routine = "allocate_sparse_hamil_row"

        write (string_row, '(I10)') row

        var_name = "sparse_hamil_"//trim(string_row)//"_elements"
        allocate(sparse_hamil(row)%elements(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, 8, this_routine, &
                         SparseHamilTags(1,row), ierr)

        var_name = "sparse_hamil_"//trim(string_row)//"_positions"
        allocate(sparse_hamil(row)%positions(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, bytes_int, this_routine, &
                         SparseHamilTags(2,row), ierr)

    end subroutine allocate_sparse_hamil_row

    subroutine deallocate_sparse_hamil()

        use FciMCData, only: sparse_matrix_info, sparse_hamil, hamil_diag

        integer :: sparse_hamil_size, i, ierr
        character (len=*), parameter :: this_routine = "allocate_sparse_hamil"

        sparse_hamil_size = size(sparse_hamil)

        do i = sparse_hamil_size, 1, -1

            deallocate(sparse_hamil(i)%elements, stat=ierr)
            call LogMemDealloc(this_routine, SparseHamilTags(1,i), ierr)

            deallocate(sparse_hamil(i)%positions, stat=ierr)
            call LogMemDealloc(this_routine, SparseHamilTags(2,i), ierr)
        end do

        deallocate(SparseHamilTags)
        deallocate(sparse_hamil)

    end subroutine deallocate_sparse_hamil

    subroutine generate_optimised_core(called_from)

        ! This routine generates a deterministic space by diagonalising a small fraction
        ! of the whole space, and choosing the basis states with the largest weights in
        ! the ground state for the deterministic states. This therefore aims to find
        ! some of the basis states with the largest weights in the true ground state.

        ! More specifically, we start with the Hartree-Fock state, and generate a first
        ! space by finding all states connected to it. We then find the ground state of
        ! the Hamiltonian in this space. Using this ground state, we pick out the basis
        ! states with the largest amplitudes (up to some user-specified criteria), and
        ! define these basis states as our new space. We then find all states connected
        ! to the states in this space, and diagonalise the Hamiltonian in this new space
        ! and again pick out the basis states with the largest weights. This process is
        ! iterated for as many time as the user requests.

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        use davidson, only: perform_davidson, davidson_eigenvalue, davidson_eigenvector, &
                            sparse_hamil_type
        use FciMCData, only: sparse_matrix_info, sparse_hamil, hamil_diag

        integer, intent(in) :: called_from
        integer(n_int), allocatable, dimension(:,:) :: ilut_store, temp_space
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: counter, i, j, k, ierr
        integer :: old_num_states, new_num_states, max_space_size, comp
        integer :: num_generation_loops, array_size
        integer(MPIArg) :: proc_space_sizes(0:nProcessors-1), disps(0:nProcessors-1), &
                           sendcounts(0:nProcessors-1), recvcount, this_proc_size
        integer, allocatable, dimension(:) :: space_cutoff_num
        real(dp), allocatable, dimension(:) :: space_cutoff_amp
        logical :: tAmplitudeCutoff, tLimitSpace
        integer(TagIntType) :: IlutTag, TempTag, FinalTag
        character (len=*), parameter :: this_routine = "generate_optimised_core"

        if (iProcIndex == root) then

            ! Use the correct set of parameters, depending this function was called
            ! from for generating the deterministic space or the trial space:
            if (called_from == called_from_semistoch) then
                num_generation_loops = num_det_generation_loops
                tAmplitudeCutoff = tDetermAmplitudeCutoff
                max_space_size = max_determ_size
                tLimitSpace = tLimitDetermSpace
                if (tAmplitudeCutoff) then
                    array_size = size(determ_space_cutoff_amp,1)
                    allocate(space_cutoff_amp(array_size))
                    space_cutoff_amp = determ_space_cutoff_amp
                else
                    array_size = size(determ_space_cutoff_num,1)
                    allocate(space_cutoff_num(array_size))
                    space_cutoff_num = determ_space_cutoff_num
                end if

            elseif (called_from == called_from_trial) then
                num_generation_loops = num_trial_generation_loops
                tAmplitudeCutoff = tTrialAmplitudeCutoff
                max_space_size = max_trial_size
                tLimitSpace = tLimitTrialSpace
                if (tAmplitudeCutoff) then
                    array_size = size(trial_space_cutoff_amp,1)
                    allocate(space_cutoff_amp(array_size))
                    space_cutoff_amp = trial_space_cutoff_amp
                else
                    array_size = size(trial_space_cutoff_num,1)
                    allocate(space_cutoff_num(array_size))
                    space_cutoff_num = trial_space_cutoff_num
                end if
            end if

            ! Allocate the stores of ilut's that will hold these deterministic states.
            ! For now, assume that we won't have a deterministic space of more than one
            ! million states. Could make this user-specified later.
            allocate(ilut_store(0:NIfTot, 1000000), stat=ierr)
            call LogMemAlloc("ilut_store", 1000000*(NIfTot+1), size_n_int, this_routine, &
                             IlutTag, ierr)
            allocate(temp_space(0:NIfTot, 1000000), stat=ierr)
            call LogMemAlloc("temp_store", 1000000*(NIfTot+1), size_n_int, this_routine, &
                             TempTag, ierr)
            ilut_store = 0
            temp_space = 0

            ! Put the Hartree-Fock state in the list first.
            ilut_store(0:NIfTot, 1) = ilutHF(0:NIfTot)

            ! old_num_states will hold the number of deterministic states in the current
            ! space. This is just 1 for now, with only the Hartree-Fock.
            old_num_states = 1

            ! Now we start the iterating loop. Find all states which are either a single or
            ! double excitation from each state in the old ilut store, and then see if they
            ! have a non-zero Hamiltonian matrix element with any state in the old ilut store:

            ! Over the total number of iterations.
            do i = 1, num_generation_loops

                write(6,'(a37,1X,i2)') "Optimised space generation: Iteration", i
                call neci_flush(6)

                ! Find all states connected to the states currently in ilut_store.
                write(6,'(a29)') "Generating connected space..."
                call neci_flush(6)
                call generate_connected_space(old_num_states, ilut_store(:, 1:old_num_states), &
                                              new_num_states, temp_space(:, 1:1000000), 1000000)
                write(6,'(a26)') "Connected space generated."
                call neci_flush(6)

                ! Add these states to the ones already in the ilut stores.
                ilut_store(:, old_num_states+1:old_num_states+new_num_states) = &
                    temp_space(:, 1:new_num_states)

                new_num_states = new_num_states + old_num_states

                call remove_repeated_states(ilut_store(:, 1:new_num_states), new_num_states)

                write(6,'(i8,1X,a13)') new_num_states, "states found."
                call neci_flush(6)

                if (tLimitSpace) call remove_high_energy_orbs(ilut_store(:, 1:new_num_states), &
                                                              new_num_states, max_space_size, .false.)

                write(6,'(a27)') "Constructing Hamiltonian..."
                call neci_flush(6)

                call calculate_sparse_hamiltonian(new_num_states, ilut_store(:,1:new_num_states))

                write (6,'(a29)') "Performing diagonalisation..."
                call neci_flush(6)

                ! Now that the Hamiltonian is generated, we can finally find the ground state of it:
                call perform_davidson(sparse_hamil_type, .false.)

                ! davidson_eigenvector now stores the ground state eigenvector. We want to use the
                ! vector whose components are the absolute values of this state:
                davidson_eigenvector = abs(davidson_eigenvector)
                ! Multiply by -1.0_dp so that the sort operation later is the right way around.
                davidson_eigenvector = -1.0_dp*davidson_eigenvector

                ! Now decide which states to keep for the next iteration. There are two ways of
                ! doing this, as decided by the user. Either all basis states with an amplitude
                ! in the ground state greater than a given value are kept (tAmplitudeCutoff =
                ! .true.), or a number of states to keep is specified and we pick the states with
                ! the biggest amplitudes (tAmplitudeCutoff = .false.).
                if (tAmplitudeCutoff) then
                    counter = 0
                    do j = 1, new_num_states
                        if (davidson_eigenvector(j) > space_cutoff_amp(i)) then
                            counter = counter + 1
                            ilut_store(:, counter) = ilut_store(:, j)
                        end if
                    end do
                    old_num_states = counter
                else
                    ! Sort the list in order of the amplitude of the states in the ground state,
                    ! from largest to smallest.
                    call sort(davidson_eigenvector(:), ilut_store(:, 1:new_num_states))

                    ! Set old_num_states to specify the number of states which should be kept.
                    old_num_states = space_cutoff_num(i)
                    if (old_num_states > new_num_states) old_num_states = new_num_states
                end if

                write(6,'(i8,1X,a12)') old_num_states, "states kept."
                call neci_flush(6)

                call deallocate_sparse_hamil()
                deallocate(hamil_diag, stat=ierr)
                call LogMemDealloc(this_routine, HDiagTag, ierr)

            end do

        end if ! If on root.

        ! When used for generating the deterministic space, send each processor the
        ! states which belong to that processor only.
        ! When used for generating the trial space, send each processor every state.
        if (called_from == called_from_semistoch) then

            if (iProcIndex == root) then
                ! Find which core each state belongs to and sort accordingly.
                call sort_space_by_proc(ilut_store, old_num_states, proc_space_sizes)

                ! Create displacement and sendcount arrays for MPIScatterV later:
                sendcounts = proc_space_sizes*(NIfTot+1)
                disps(0) = 0
                do i = 1, nProcessors
                    disps(i) = sum(proc_space_sizes(0:i-1))*(NIfTot+1)
                end do
            end if

            ! Send the number of states on each processor to the corresponding processor.
            call MPIScatter(proc_space_sizes, this_proc_size, ierr)
            recvcount = this_proc_size*(NIfTot+1)
            ! Finally send the actual states to the SpawnedParts array.
            call MPIScatterV(ilut_store, sendcounts, disps, &
                             SpawnedParts(:, 1:this_proc_size), recvcount, ierr)

            determ_proc_sizes(iProcIndex) = this_proc_size

            ! Set the flags and the amplitude on the HF.
            do i = 1, this_proc_size
                call set_flag(SpawnedParts(:,i), flag_deterministic)
                if (tTruncInitiator) then
                    call set_flag(SpawnedParts(:,i), flag_is_initiator(1))
                    call set_flag(SpawnedParts(:,i), flag_is_initiator(2))
                end if
                comp = DetBitLT(SpawnedParts(:,i), ilutHF, NIfD, .false.)
                if (comp == 0) SpawnedParts(nOffSgn,i) = CurrentDets(nOffSgn,1)
            end do

        else if (called_from == called_from_trial) then

           if (iProcIndex == root) then
               trial_space(:, 1:old_num_states) = ilut_store(:, 1:old_num_states)
               trial_space_size = old_num_states
           end if

           call MPIBCast(trial_space_size, 1, root)
           call MPIBCast(trial_space(:, 1:trial_space_size), &
                         trial_space_size*(NIfTot+1), root)

        end if

        ! Finally, deallocate arrays.
        if (allocated(space_cutoff_num)) deallocate(space_cutoff_num)
        if (allocated(space_cutoff_amp)) deallocate(space_cutoff_amp)
        if (iProcIndex == root) then
            deallocate(temp_space, stat=ierr)
            call LogMemDealloc(this_routine, TempTag, ierr)
            deallocate(ilut_store, stat=ierr)
            call LogMemDealloc(this_routine, IlutTag, ierr)
        end if

    end subroutine generate_optimised_core

    subroutine generate_connected_space(original_space_size, original_space, connected_space_size, &
        connected_space, storage_space_size, tSinglesOnly)

        integer, intent(in) :: original_space_size
        integer(n_int), intent(in) :: original_space(0:NIfTot, original_space_size)
        integer, intent(in) :: storage_space_size
        integer, intent(out) :: connected_space_size
        integer(n_int), intent(out) :: connected_space(0:NIfTot, storage_space_size)
        logical, intent(in), optional :: tSinglesOnly

        integer(n_int) :: ilut(0:NIfTot), ilut_tmp(0:NIfTot)
        integer :: nI(nel)
        integer :: i, counter
        logical :: first_loop, connected, tSkipDoubles
        type(excit_store), target :: gen_store
        character (len=1024) :: storage_space_string

        ! By default, generate both singles and doubles (the whole connected space).
        if (present(tSinglesOnly)) then
            tSkipDoubles = tSinglesOnly
        else
            tSkipDoubles = .false.
        end if

        connected_space_size = 0

        ! Over all the states in the original list:
        do i = 1, original_space_size

            call decode_bit_det(nI, original_space(:,i))

            ! The singles:

            ! Ensure ilut(0) \= -1 so that the loop can be entered.
            ilut(0) = 0
            first_loop = .true.

            ! When no more basis functions are found, this value is returned and the
            ! loop is exited.
            do while(ilut(0) /= -1)

                ! The first time the enumerating subroutine is called, setting ilut(0)
                ! to -1 tells it to  initialise everything first.
                if (first_loop) then
                    ilut(0) = -1
                    first_loop = .false.
                end if

                ! Find the next state.
                call enumerate_all_single_excitations (original_space(:,i), nI, ilut, &
                                                                             gen_store)

                ! If a new state was not found.
                if (ilut(0) == -1) exit

                ! If using HPHFs, and if this isn't the allowed state for this HPHF, find
                ! the the allowed state and continue with this as ilut.
                if (tHPHF .and. (.not. TestClosedShellDet(ilut))) then
                    if (.not. IsAllowedHPHF(ilut(0:NIfD))) then
                        ilut_tmp = ilut
                        call FindExcitBitDetSym(ilut_tmp, ilut)
                    end if
                end if

                ! Check if any of the states in the old space is connected to the state
                ! just generated. If so, this function returns the value true for the
                ! logical connected. In this case, this state is then added to the new ilut
                ! store.
                call check_if_connected_to_old_space(original_space, nI, ilut, &
                    original_space_size, connected_space_size, connected)
                if (connected_space_size > storage_space_size) then
                    write (storage_space_string, '(I10)') storage_space_size
                    call stop_all("generate_connected_space","No space left in storage array &
                    &for the next connected space state. "//trim(storage_space_string)//" &
                    &elements were allocated and this number has been exceeded.")
                end if
                if (connected) connected_space(0:NIfD, connected_space_size) = ilut(0:NIfD)

            end do

            ! If only generating the singles space.
            if (tSkipDoubles) cycle

            ! The doubles:

            ilut(0) = 0
            first_loop = .true.

            do while(ilut(0) /= -1)

                if (first_loop) then
                    ilut(0) = -1
                    first_loop = .false.
                end if

                call enumerate_all_double_excitations (original_space(:,i), nI, ilut, &
                                                                             gen_store)

                if (ilut(0) == -1) exit

                ! If using HPHFs, and if this isn't the allowed state for this HPHF, find
                ! the the allowed state and continue with this as ilut.
                if (tHPHF .and. (.not. TestClosedShellDet(ilut))) then
                    if (.not. IsAllowedHPHF(ilut(0:NIfD))) then
                        ilut_tmp = ilut
                        call FindExcitBitDetSym(ilut_tmp, ilut)
                    end if
                end if

                call check_if_connected_to_old_space(original_space, nI, ilut, &
                    original_space_size, connected_space_size, connected)
                if (connected_space_size > storage_space_size) then
                    write (storage_space_string, '(I10)') storage_space_size
                    call stop_all("generate_connected_space","No space left in storage array &
                    &for the next connected space state. "//trim(storage_space_string)//" &
                    &elements were allocated and this number has been exceeded.")
                end if
                if (connected) connected_space(0:NIfD, connected_space_size) = ilut(0:NIfD)

            end do

        end do

    end subroutine generate_connected_space

    subroutine check_if_connected_to_old_space(ilut_store_old, nI, ilut, num_states, &
                                                                         counter, connected)

        ! This subroutine loops over all states in the old space and see if any of them is
        ! connected to the state just generated. If so, it also adds one to counter.

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

    subroutine generate_sing_doub_determinants(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        type(excit_store), target :: gen_store
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: space_size
        integer :: i, ierr

        ! This condition tells the enumerating subroutines to initialise the loop.
        ilut(0) = -1

        ! Find the first single excitation.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        ! Subroutine to add this state to the spawned list if on this processor.
        call add_basis_state_to_list(ilut, called_from)

        ! When no more basis functions are found, this value is returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)
            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)

            ! If a determinant is returned (if we did not find the final one last time.)
            if (ilut(0) /= -1) call add_basis_state_to_list(ilut, called_from)
        end do

        ! Now generate the double excitations...

        ilut(0) = -1
        ! The first double excitation.
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_basis_state_to_list(ilut, called_from)

        do while(ilut(0) /= -1)
            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)

            if (ilut(0) /= -1) call add_basis_state_to_list(ilut, called_from)
        end do

        if (called_from == called_from_trial) then
            if (tLimitTrialSpace) call remove_high_energy_orbs(trial_space(:, 1:trial_space_size), &
                                                             trial_space_size, max_trial_size, .false.)
        end if

    end subroutine generate_sing_doub_determinants

    subroutine generate_sing_doub_csfs(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        type(excit_store), target :: gen_store
        integer(n_int) :: ilutHF_loc(0:NIfTot), ilut(0:NIfTot)
        integer :: HFDet_loc(nel), nI(nel)
        integer :: exflag
        logical :: first_loop

        if (tFixLz) call stop_all("generate_sing_doub_csfs", "The CSF generating routine &
            &does not work when Lz symmetry is applied.")

        ! Setting the first two bits of this flag tells the generating subroutine to
        ! generate both single and double spatial excitations.
        exflag = 0
        exflag = ibset(exflag, 0)
        exflag = ibset(exflag, 1)

        ! For Stot /= 0, the HF state will be a CSF. For the purpose of
        ! generating all spatial orbitals, we just want a determinant, so use a
        ! state without the CSF information.
        HFdet_loc = iand(HFDet, csf_orbital_mask)
        ilutHF_loc = ilutHF
        ilutHF_loc(NOffY) = 0

        ! Ensure ilut(0) \= -1 so that the loop can be entered.
        ilut(0) = 0
        first_loop = .true.
        do while (ilut(0) /= -1)

            if (first_loop) then
                ilut(0) = -1
                first_loop = .false.
            end if

            ! Generate the next spatial excitation (orbital configuration).
            call enumerate_spatial_excitations(ilutHF_loc, HFDet_loc, ilut, exflag, gen_store)

            if (ilut(0) == -1) exit

            ! Find the nI representation.
            call decode_bit_det(nI, ilut)

            call generate_all_csfs_from_orb_config(ilut, nI, called_from)

        end do

    end subroutine generate_sing_doub_csfs

    subroutine generate_cas(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        type(BasisFN) :: CASSym
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int), allocatable, dimension(:,:) :: ilut_store
        integer :: HFdet_loc(nel)
        integer :: num_active_orbs, elec, nCASDet, i, j, counter, comp, ierr
        integer, allocatable :: CASBrr(:), CASRef(:)
        integer(n_int) :: cas_bitmask(0:NIfD), cas_not_bitmask(0:NIfD)
        integer, pointer :: CASDets(:,:) => null()
        integer :: OccOrbs, VirtOrbs, iCASDet
        integer(TagIntType) :: CASDetsTag, IlutTag
        character (len=*), parameter :: this_routine = "generate_cas"

        if (called_from == called_from_semistoch) then
            OccOrbs = OccDetermCASOrbs
            VirtOrbs = VirtDetermCASOrbs
        elseif (called_from == called_from_trial) then
            OccOrbs = OccTrialCASOrbs
            VirtOrbs = VirtTrialCASOrbs
        end if

        ! This option should be true. It tells the subroutine gndts to only consider states
        ! with an Ms value in the correct spin subspace.
        if (.not. tSpn) call stop_all("generate_cas", "tSpn is not set to true.")

        ! The total number of orbitals in the active space:
        num_active_orbs = OccOrbs + VirtOrbs
        allocate(CASBrr(1:num_active_orbs))
        allocate(CASRef(1:OccOrbs))
        do i = 1, num_active_orbs
            ! Run through the cas space, and create an array which will map these orbtials to the
            ! orbitals they actually represent.
            ! i.e. CASBRR(1) will store the lowest energy orbital in the CAS space and
            ! CASBRR(num_active_orbs) will store the highest energy orbital in the CAS space.
            CASBrr(i) = BRR(i + (nel - OccOrbs))
        end do

        ! Create a bit mask which has 1's in the bits which represent active orbitals and 0's in
        ! all other orbitals.
        cas_bitmask = 0
        do i = 1, num_active_orbs
            set_orb(cas_bitmask, CASBrr(i))
        end do
        ! Create a bit mask which has 0's in the bits which represent active orbitals and 1's in
        ! all other orbitals.
        cas_not_bitmask = not(cas_bitmask)

        ! For Stot /= 0, the HF state will be a CSF. For the purpose of generating all spatial
        ! orbitals, we just want a determinant, so use a state without the CSF information.
        HFdet_loc = iand(HFDet, csf_orbital_mask)

        elec = 1
        do i = nel-OccOrbs+1, nel
            ! CASRef(elec) will store the orbital number of the electron elec in the reference
            ! state, HFDet. elec runs from 1 to the number of electrons in the active space.
            CASRef(elec) = HFDet_loc(i)
            elec = elec + 1
        end do

        call GetSym(CASRef, OccOrbs, G1, nBasisMax, CASSym)

        ! First, generate all excitations so we know how many there are, to allocate the arrays.
        call gndts(OccOrbs, num_active_orbs, CASBrr, nBasisMax, CASDets, &
                              .true., G1, tSpn, LMS, .true., CASSym, nCASDet, iCASDet)

        if (nCASDet == 0) call stop_all("generate_cas","No CAS determinants found.")

        ! Now allocate the array CASDets. CASDets(:,i) will store the orbitals in the active space
        ! which are occupied in the i'th determinant generated.
        allocate(CASDets(OccOrbs, nCASDet), stat=ierr)
        call LogMemAlloc("CASDets", OccOrbs*nCASDet, 8, this_routine, CASDetsTag, ierr)
        CASDets(:,:) = 0

        if (tCASCore) then
            allocate(ilut_store(nCASDet-1, 0:NIfTot), stat=ierr)
            call LogMemAlloc("ilut_store", (nCASDet-1)*(NIfTot+1), size_n_int, this_routine, &
                             IlutTag, ierr)
            ilut_store = 0
            counter = 1
        end if

        ! Now fill up CASDets...
        call gndts(OccOrbs, num_active_orbs, CASBrr, nBasisMax, CASDets, &
                              .false., G1, tSpn, LMS, .true., CASSym, nCASDet, iCASDet)

        do i = 1, nCASDet
            ! First, create the bitstring representing this state:
            ! Start from the HF determinant and apply cas_not_bitmask to clear all active space
            ! orbitals.
            ilut(0:NIfD) = iand(ilutHF(0:NIfD), cas_not_bitmask)
            ! Then loop through the occupied orbitals in the active space, stored in CASDets(:,i),
            ! and set the corresponding bits.
            do j = 1, OccOrbs
                set_orb(ilut, CASDets(j,i))
            end do

            ! The function gndts always outputs the reference state. This is already included, so
            ! we want to cycle when we get to this state to ignore it.
            ! comp will be 0 if ilut and ilutHF are the same.
            comp = DetBitLT(ilut, ilutHF, NIfD, .false.)
            if (comp == 0) cycle

            if (tDeterminantCore) then
                ! Now that we have fully generated the determinant, add it to the main list.
                call add_basis_state_to_list(ilut, called_from)
            else if (tCSFCore) then
                ilut_store(counter, 0:NifD) = ilut(0:NifD)
                counter = counter + 1
            end if

        end do

        if (tCSFCore) call create_CAS_csfs(ilut_store, nCASDet-1, called_from)

        if (called_from == called_from_semistoch) then
            allocate(cas_determ_bitmask(0:NIfD))
            allocate(cas_determ_not_bitmask(0:NIfD))
            cas_determ_bitmask = cas_bitmask
            cas_determ_not_bitmask = cas_not_bitmask
        end if

        deallocate(CASBrr)
        deallocate(CASRef)
        deallocate(CASDets, stat=ierr)
        call LogMemDealloc(this_routine, CASDetsTag, ierr)
        deallocate(ilut_store, stat=ierr)
        call LogMemDealloc(this_routine, IlutTag, ierr)

    end subroutine generate_cas

    subroutine generate_low_energy_core(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        integer :: i, num_loops, num_states_keep, comp, ierr
        integer :: old_num_states, new_num_states, max_space_size, low_e_excit
        logical :: tAllDoubles, tSinglesOnly
        integer(n_int), allocatable, dimension(:,:) :: ilut_store, temp_space
        integer(TagIntType) :: IlutTag, TempTag
        character (len=*), parameter :: this_routine = "generate_low_energy_core"

        ! Use the correct set of parameters, depending this function was called
        ! from for generating the deterministic space or the trial space:
        if (called_from == called_from_semistoch) then
            low_e_excit = low_e_core_excit
            tAllDoubles = tLowECoreAllDoubles
            num_states_keep = low_e_core_num_keep
            max_space_size = max_determ_size
        elseif (called_from == called_from_trial) then
            low_e_excit = low_e_trial_excit
            tAllDoubles = tLowETrialAllDoubles
            num_states_keep = low_e_trial_num_keep
            max_space_size = max_trial_size
        end if

        ! low_e_excit holds the maximum excitation level to go to. In the first loop,
        ! generate both single and double excitations.
        num_loops = max(low_e_excit-1, 1)

        if (low_e_excit == 1) call warning_neci("generate_low_energy_core", "You asked for &
                                                &singles only, but both singles and doubles &
                                                &will be generated in the first iteration.")

        allocate(ilut_store(0:NIfTot, 1000000), stat=ierr)
        call LogMemAlloc("ilut_store", 1000000*(NIfTot+1), size_n_int, this_routine, &
                         IlutTag, ierr)
        allocate(temp_space(0:NIfTot, 1000000), stat=ierr)
        call LogMemAlloc("temp_store", 1000000*(NIfTot+1), size_n_int, this_routine, &
                         TempTag, ierr)
        ilut_store = 0
        temp_space = 0

        ! Put the Hartree-Fock state in the list first.
        ilut_store(0:NIfTot, 1) = ilutHF(0:NIfTot)

        ! old_num_states will hold the number of deterministic states in the current
        ! space. This is just 1 for now, with only the Hartree-Fock.
        old_num_states = 1
        new_num_states = 1

        do i = 1, num_loops

            write(6,'(a38,1X,i2)') "Low energy space generation: Iteration", i
            call neci_flush(6)

            ! The number of states in the list to work with, after the last iteration.
            old_num_states = min(new_num_states, num_states_keep)

            ! In the first iteration, generate all singles and doubles.
            if (i == 1) then
                tSinglesOnly = .false.
            else
                tSinglesOnly = .true.
            end if

            ! Find all *single* excitations (not doubles) to the states in ilut_store.
            call generate_connected_space(old_num_states, ilut_store(:, 1:old_num_states), &
                                          new_num_states, temp_space(:, 1:1000000), &
                                          1000000, tSinglesOnly)

            ! Add these states to the ones already in the ilut stores.
            ilut_store(:, old_num_states+1:old_num_states+new_num_states) = &
                temp_space(:, 1:new_num_states)

            new_num_states = new_num_states + old_num_states

            call remove_repeated_states(ilut_store(:, 1:new_num_states), new_num_states)

            write(6,'(i8,1X,a13)') new_num_states, "states found."
            call neci_flush(6)

            call sort_states_by_energy(ilut_store, new_num_states, tAllDoubles)

            ! If user has asked to keep all singles and doubles but has not asked for enough
            ! states to be kept.
            if (i == 1 .and. new_num_states > num_states_keep .and. tAllDoubles) &
                call warning_neci("generate_low_energy_core", "You have asked to keep all &
                                   &singles and doubles, but the maximum number of states &
                                   &you have asked to keep is to small for this. Some singles &
                                   &or doubles will not be kept.")

        end do

        if (max_space_size /= 0) new_num_states = min(new_num_states, max_space_size)

        ! Add all states (except for the HF state, included already) to the appropriate
        ! space.
        do i = 1, new_num_states
            comp = DetBitLT(ilut_store(:, i), ilutHF, NIfD, .false.)
            if (comp == 0) cycle
            call add_basis_state_to_list(ilut_store(:, i), called_from)
        end do

        deallocate(temp_space, stat=ierr)
        call LogMemDealloc(this_routine, TempTag, ierr)
        deallocate(ilut_store, stat=ierr)
        call LogMemDealloc(this_routine, IlutTag, ierr)

    end subroutine generate_low_energy_core

    subroutine sort_states_by_energy(ilut_list, num_states, tSortDoubles)

        ! Note: If requested, sort doubles to the top, then sort by energy.

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

    subroutine create_CAS_csfs(ilut_store, num_dets, called_from)

        integer, intent(in) :: num_dets
        integer(n_int), intent(inout) :: ilut_store(num_dets, 0:NIfTot)
        integer, intent(in) :: called_from
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: i, j, comp
        integer :: orb1, orb2

        do i = 1, num_dets

            nI = 0
            ilut = ilut_store(i,:)
            call decode_bit_det(nI, ilut)
            nI = iand(nI, csf_orbital_mask)
            print *, "nI:", nI
            print *, "ilut:", ilut

            do j = 1, nel
                orb1 = nI(j)
                orb2 = ab_pair(orb1)

                if (is_alpha(orb1) .and. IsNotOcc(ilut, orb2)) then
                    clr_orb(ilut, orb1)
                    set_orb(ilut, orb2)
                end if
            end do

        end do

        call sort(ilut_store(:,:), ilut_lt, ilut_gt)

        do i = 1, num_dets
            if (i > 1) then
                comp = DetBitLT(ilut_store(i,:), ilut_store(i-1,:), NIfD, .false.)
                if (comp == 0) cycle
            end if

            call decode_bit_det(nI, ilut_store(i,:))

            call generate_all_csfs_from_orb_config(ilut_store(i,:), nI, called_from)
        end do

    end subroutine create_CAS_csfs

    subroutine generate_all_csfs_from_orb_config(ilut, nI, called_from)

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(inout) :: nI(nel)
        integer, intent(in) :: called_from
        integer :: i, n_open, num_csfs, max_num_csfs
        integer, allocatable, dimension(:,:) :: yama_symbols

        ! Can't possibly have more open orbitals than nel.
        ! TODO: Check that num_csfs *always* increases with n_open.
        max_num_csfs = get_num_csfs(nel, Stot)

        ! Use max_num_csfs to allocate the array of Yamanouchi symbols to be large enough.
        allocate(yama_symbols(max_num_csfs, nel))
        yama_symbols = 0

        ! Find the number of open orbitals.
        n_open = count_open_orbs(ilut)
        ! Find the number of CSFs (and hence Yamanouchi symbols) with this value of Stot for
        ! this orbital configuration.
        num_csfs = get_num_csfs(n_open, Stot)

        ! Enumerate the list of all possible Yamanouchi symbols for this orbital
        ! configuration.
        call csf_get_yamas(n_open, Stot, yama_symbols(1:num_csfs,1:n_open), num_csfs)

        ! Loop over all Yamanouchi symbols.
        do i = 1, num_csfs
            ! If n_open = 0 then we just have a determinant.
            if (n_open > 0) then
                ! Encode the Yamanouchi symbol in nI representation.
                call csf_apply_yama(nI, yama_symbols(i,1:n_open))
                ! Encode the Yamanouchi symbol in the ilut representation.
                call get_csf_bit_yama(nI, ilut(nOffY:nOffY+nIfY-1))
            end if
            ! Finally add the CSF to the SpawnedParts list.
            call add_basis_state_to_list(ilut, called_from, nI)
        end do

    end subroutine generate_all_csfs_from_orb_config

    subroutine add_basis_state_to_list(ilut, called_from, nI_in)

        ! This subroutine, when called from init_semi_stochastic, takes a bitstring,
        ! finds the nI representation, decides if the basis state lives on this
        ! processor and, if so, adds it to the spawned list. It also increases the
        ! vector sizes being calculated.

        ! When called from the trial wavefunction generation code, only the root
        ! processor accesses this code, and each state is added to the list of iluts,
        ! trial_space on this one processor.

        ! In: ilut - The determinant in a bitstring form.
        ! In: called_from - Integer to specify which part of the code this routine was
        !     called from, and hence which space this state should be added to.
        ! In (optional) : nI_in - A list of the occupied orbitals in the determinant.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: called_from
        integer, optional :: nI_in(nel)
        integer :: nI(nel)
        integer :: flags
        integer :: proc
        real(dp) :: sgn(lenof_sign)

        if (called_from == called_from_semistoch) then

            ! If using HPHFs then only allow the correct HPHFs to be added to the list.
            if (tHPHF) then
                if (.not. IsAllowedHPHF(ilut(0:NIfD))) return
            end if

            sgn = 0.0_dp
            ! Flag to specify that these basis states are in the deterministic space.
            flags = 0
            flags = ibset(flags, flag_deterministic)
            if (tTruncInitiator) then
                flags = ibset(flags, flag_is_initiator(1))
                flags = ibset(flags, flag_is_initiator(2))
            end if

            ! Find the nI representation of determinant.
            if (present(nI_in)) then
                nI = nI_in
            else
                call decode_bit_det(nI, ilut)
            end if

            ! Find the processor which this state belongs to.
            proc = DetermineDetNode(nI,0)

            ! Keep track of the size of the deterministic space on this processor.
            if (proc == iProcIndex) determ_proc_sizes(proc) = determ_proc_sizes(proc) + 1

            ! If this determinant belongs to this processor, add it to the main list.
            if (proc == iProcIndex) call encode_bit_rep(SpawnedParts(0:NIfTot, &
                determ_proc_sizes(iProcIndex)), ilut(0:nIfDBO), sgn, flags)

        else if (called_from == called_from_trial) then

            if (tHPHF) then
                if (.not. IsAllowedHPHF(ilut)) return
            end if

            ! Keep track of the size of the trial space on this processor.
            trial_space_size = trial_space_size + 1

            trial_space(0:NIfTot, trial_space_size) = ilut(0:NIfTot)

        end if

    end subroutine add_basis_state_to_list

    subroutine check_if_in_determ_space(ilut, tInDetermSpace)

        ! This subroutine takes in a state in its bitstring representation and finds out
        ! whether it is in the deterministic space or not.

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
        end if

    end subroutine check_if_in_determ_space

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

end module semi_stochastic
