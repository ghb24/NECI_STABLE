#include "macros.h"

module semi_stochastic

    use bit_rep_data, only: flag_deterministic, nIfDBO, nOffY, nIfY, nOffFlag, &
                            flag_determ_parent, deterministic_mask, determ_parent_mask, &
                            flag_is_initiator
    use bit_reps, only: NIfD, NIfTot, decode_bit_det, encode_bit_rep, set_flag
    use CalcData, only: tRegenDiagHEls, tau, tSortDetermToTop, tTruncInitiator, DiagSft, &
                        occCASorbs, virtCASorbs
    use csf, only: csf_get_yamas, get_num_csfs, get_csf_bit_yama, csf_apply_yama
    use csf_data, only: iscsf, csf_orbital_mask
    use constants
    use DetBitOps, only: ilut_lt, ilut_gt, count_open_orbs, FindBitExcitLevel, DetBitLT, &
                         count_set_bits
    use Determinants, only: get_helement
    use enumerate_excitations
    use FciMCData, only: HFDet, ilutHF, iHFProc, Hii, CurrentDets, CurrentH, &
                         deterministic_proc_sizes, deterministic_proc_indices, &
                         full_determ_vector, partial_determ_vector, core_hamiltonian, &
                         determ_space_size, TotWalkers, TotWalkersOld, &
                         indices_of_determ_states, ProjEDet
    use hash , only: DetermineDetNode
    use hphf_integrals, only: hphf_diag_helement
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIBarrier, &
                             MPIAllGatherV
    use sort_mod, only: sort
    use sym_mod, only: getsym
    use SystemData, only: nel, tHPHF, tCSFCore, tDeterminantCore, tDoublesCore, tCASCore, &
                          tOptimisedCore, Stot, lms, BasisFN, nBasisMax, BRR, tSpn, &
                          cas_bitmask, cas_not_bitmask, num_det_generation_loops, &
                          determ_space_cutoff_amp, determ_space_cutoff_num, tAmplitudeCutoff, &
                          determ_ground_state

    implicit none

contains

    subroutine init_semi_stochastic()

        ! Initialise the semi-stochastic information. This includes enumerating a list
        ! of all determinants or CSFs in the deterministic space and calculating and
        ! storing the resulting Hamiltonian matrix elements. The lists which will store
        ! the psip vectors in the deterministic space are also allocated.

        integer :: i, j, IC

        ! Initialise the deterministic masks.
        deterministic_mask = 0
        determ_parent_mask = 0
        deterministic_mask = ibset(deterministic_mask, flag_deterministic)
        determ_parent_mask = ibset(determ_parent_mask, flag_determ_parent)

        allocate(deterministic_proc_sizes(0:nProcessors-1))
        allocate(deterministic_proc_indices(0:nProcessors-1))
        deterministic_proc_sizes = 0
        deterministic_proc_indices = 0

        ! Count the Hartree-Fock determinant (The Hartree-Fock determinant has already been
        ! added to the main list as usual in the subroutine InitFCIMCCalcPar).
        deterministic_proc_sizes(iHFProc) = 1

        ! The following subroutines call the enumerating subroutines to create all excitations
        ! and add these states to the main list, CurrentDets, on the correct processor. As
        ! they do this, they count the size of the deterministic space on each processor.
        if (tDeterminantCore) then
            tSortDetermToTop = .true.
            if (tDoublesCore) then
                call generate_sing_doub_determinants()
            else if (tCASCore) then
                call generate_cas()
            end if
        else if (tCSFCore) then
            tSortDetermToTop = .true.
            if (tDoublesCore) then
                call generate_sing_doub_csfs()
            else if (tCASCore) then
                call generate_cas()
            end if
        else if (tOptimisedCore) then
            tSortDetermToTop = .false.
            call generate_optimised_core()
        end if

        ! We now know the size of the deterministic space on each processor, so now find
        ! the total size of the space and also allocate vectors to store psip amplitudes
        ! and the deterministic Hamiltonian.
        determ_space_size = sum(deterministic_proc_sizes)
        allocate(full_determ_vector(determ_space_size))
        allocate(partial_determ_vector(deterministic_proc_sizes(iProcIndex)))
        allocate(core_hamiltonian(deterministic_proc_sizes(iProcIndex), determ_space_size))
        full_determ_vector = 0.0_dp
        partial_determ_vector = 0.0_dp

        ! Write space sizes to screen.
        write(6,*) "Total size of deterministic space:", determ_space_size
        write(6,*) "Size of deterministic space on this processor:", &
                    deterministic_proc_sizes(iProcIndex)

        ! Also allocate the vector to store the positions of the deterministic states in
        ! CurrentDets.
        allocate(indices_of_determ_states(deterministic_proc_sizes(iProcIndex)))

        TotWalkers = int(deterministic_proc_sizes(iProcIndex), int64)
        TotWalkersOld = int(deterministic_proc_sizes(iProcIndex), int64)

        ! Calculate the indices in the full vector at which the various processors take
        ! over, relative to the first index position in the vector (ie, the first value
        ! in this vector will be 0, *not* 1 - this is required for mpiallgatherv later).
        do i = 0, nProcessors-1
            if (i == 0) then
                deterministic_proc_indices(i) = 0
            else
                deterministic_proc_indices(i) = deterministic_proc_indices(i-1) + &
                                                  deterministic_proc_sizes(i-1)
            end if
        end do

        ! Sort the list of basis states so that it is correctly ordered in the predefined
        ! order which is always kept throughout the simulation.
        call sort(CurrentDets(:,1:deterministic_proc_sizes(iProcIndex)), ilut_lt, ilut_gt)

        ! Calculate and store the deterministic Hamiltonian.
        call calculate_det_hamiltonian_normal()

    end subroutine init_semi_stochastic

    subroutine calculate_det_hamiltonian_normal()

        integer :: i, j, iproc, col_index, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        real(dp), allocatable, dimension(:,:) :: test_hamiltonian

        ! temp_store is storage space for bitstrings so that the Hamiltonian matrix
        ! elements can be calculated.
        allocate(temp_store(0:NIfTot, maxval(deterministic_proc_sizes)))

        ! Calcuation of the Hamiltonian matrix elements to be stored. Loop over each
        ! processor in order and broadcast the sorted vector of bitstrings from the
        ! processor to all other processors so that all matrix elements can be found.
        do iproc = 0, nProcessors-1

            ! If we are broadcasting from this processor next, transfer the bitstrings
            ! to the array temp_store.
            if (iproc == iProcIndex) temp_store(:,1:deterministic_proc_sizes(iproc)) = &
                                      CurrentDets(:,1:deterministic_proc_sizes(iproc))

            ! Perform the broadcasting to other all other processors.
            call MPIBCast(temp_store, size(temp_store), iproc)

            ! The index of the column before the first column of the block of the
            ! Hamiltonian currently being calculated.
            col_index = deterministic_proc_indices(iproc)

            ! Loop over all the elements in the block of the Hamiltonian corresponding
            ! to these two prcoessors.
            do i = 1, deterministic_proc_sizes(iProcIndex)
                do j = 1, deterministic_proc_sizes(iproc)

                    ! Get nI form of the basis functions.
                    call decode_bit_det(nI, CurrentDets(:, i))
                    call decode_bit_det(nJ, temp_store(:, j))

                    ! If on the diagonal of the Hamiltonian.
                    if ((iProcIndex == iproc) .and. (i == j)) then
                        core_hamiltonian(i, col_index + j) = &
                                             get_helement(nI, nJ, 0) - Hii
                        ! We calculate and store CurrentH at this point for ease.
                        if (.not.tRegenDiagHEls) CurrentH(1,i) = &
                                             core_hamiltonian(i, col_index + j) - Hii
                    else
                        core_hamiltonian(i, col_index + j) = &
                                         get_helement(nI, nJ, CurrentDets(:, i), &
                                                            temp_store(:, j))
                    end if

                end do
            end do

        end do

        call MPIBarrier(ierr)

        deallocate(temp_store)

    end subroutine calculate_det_hamiltonian_normal

    subroutine generate_optimised_core()

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

        type(excit_store), target :: gen_store
        integer(n_int), allocatable, dimension(:,:) :: ilut_store_old, ilut_store_new
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel), nJ(nel)
        real(dp), allocatable, dimension(:,:) :: hamiltonian
        real(dp), allocatable, dimension(:) :: work
        integer :: lwork, info
        integer :: counter, i, j, k
        integer :: old_num_det_states, new_num_det_states, comp
        logical :: connected, first_loop

        ! First, allocate the stores of ilut's that will hold these deterministic states.
        ! For now, assume that we won't have a deterministic space of more than one
        ! million states. Could make this user-specified later.
        ! ilut_store_old will store the states in the previously generated deterministic
        ! space from which we want to find all connected states.
        allocate(ilut_store_old(0:NIfTot, 1000000))
        allocate(ilut_store_new(0:NIfTot, 1000000))
        ilut_store_old = 0
        ilut_store_new = 0

        ! Put the Hartree-Fock state in the old list first.
        ilut_store_old(0:NIfTot, 1) = ilutHF(0:NIfTot)
        ! Also put it in the new list, as it is connected to itself.
        ilut_store_new(0:NIfTot, 1) = ilutHF(0:NIfTot)

        counter = 1

        ! old_num_det_states will hold the number of deterministic states in the current
        ! space. This is just 1 for now, with only the Hartree-Fock.
        old_num_det_states = 1

        ! Now we start the iterating loop. Find all states which are either a single or
        ! double connection to each state in the old ilut store, and then see if they
        ! have a non-zero Hamiltonian matrix element with any state in the old
        ! ilut store:

        ! Over the total number of iterations.
        do i = 1, num_det_generation_loops

            write(6,*) "Iteration:", i
            call neci_flush(6)

            ! Over all the states in the old ilut store.
            do j = 1, old_num_det_states

                call decode_bit_det(nI, ilut_store_old(:,j))

                ! First, the singles:

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
                    call enumerate_all_single_excitations (ilut_store_old(:,j), nI, ilut, &
                                                                                 gen_store)

                    ! If a new state was not found.
                    if (ilut(0) == -1) exit

                    ! Check if any of the states in the old space is connected to the state
                    ! just generated. If so, this function adds one to counter and adds the
                    ! generated state to the new ilut store.
                    call check_if_connected_to_old_space(ilut_store_old, nI, ilut, &
                                                     old_num_det_states, counter, connected)
                    if (connected) ilut_store_new(0:NIfD, counter) = ilut(0:NIfD)

                end do

                ! Now for the doubles:

                ilut(0) = 0
                first_loop = .true.

                do while(ilut(0) /= -1)

                    if (first_loop) then
                        ilut(0) = -1
                        first_loop = .false.
                    end if

                    call enumerate_all_double_excitations (ilut_store_old(:,j), nI, ilut, &
                                                                                 gen_store)

                    call check_if_connected_to_old_space(ilut_store_old, nI, ilut, &
                                                    old_num_det_states, counter, connected)
                    if (connected) ilut_store_new(0:NIfD, counter) = ilut(0:NIfD)

                end do

            end do


            ! By this point we should have all states connected to the states in the old ilut
            ! store now stored in ilut_store_new, and the number of such states is counter.
            new_num_det_states = counter
            counter = 1

            ! Perform annihilation-like steps to remove repeated states.
            call sort(ilut_store_new(:, 1:new_num_det_states), ilut_lt, ilut_gt)

            do j = 2, new_num_det_states
                comp = DetBitLT(ilut_store_new(:, j-1), ilut_store_new(:, j), NIfD, .false.)
                if (comp /= 0) then
                    counter = counter + 1
                    ilut_store_new(:, counter) = ilut_store_new(:, j)
                end if
            end do

            ! After deleting any states which are in the new ilut store twice, the new total
            ! number of states in ilut_store_new.
            new_num_det_states = counter

            write(6,*) new_num_det_states, "states found. Constructing Hamiltonian..."
            call neci_flush(6)

            allocate(hamiltonian(new_num_det_states, new_num_det_states))
            allocate(determ_ground_state(new_num_det_states))
            hamiltonian = 0.0_dp
            determ_ground_state = 0.0_dp

            do j = 1, new_num_det_states

                ! Get nI form of the basis functions.
                call decode_bit_det(nI, ilut_store_new(:, j))

                do k = j, new_num_det_states

                    call decode_bit_det(nJ, ilut_store_new(:, k))

                    ! If on the diagonal of the Hamiltonian.
                    if (j == k) then
                        hamiltonian(j, k) = get_helement(nI, nJ, 0)
                    else
                        hamiltonian(j, k) = get_helement(nI, nJ, ilut_store_new(:, j), &
                                                                 ilut_store_new(:, k))
                        hamiltonian(k, j) = hamiltonian(j, k)
                    end if

                end do

            end do

            ! Allocate work space for the diagonaliser, dsyev.
            lwork = max(1,3*new_num_det_states-1)
            allocate(work(lwork))

            write (6,*) "Performing diagonalisation..."
            call neci_flush(6)

            ! Now that the Hamiltonian is generated, we can finally diagonalise it:
            call dsyev('V', 'U', new_num_det_states, hamiltonian, new_num_det_states, &
                       determ_ground_state, work, lwork, info)

            determ_ground_state = hamiltonian(:,1)

            write(6,*) "Diagonalisation complete."
            call neci_flush(6)

            ! Hamiltonian now stores the eigenvectors. determ_ground_state(:) stores the
            ! ground-state eigenvector.

            if (tAmplitudeCutoff) then
                counter = 0
                do j = 1, new_num_det_states
                    if (abs(determ_ground_state(j)) > determ_space_cutoff_amp(i)) then
                        counter = counter + 1
                        ilut_store_old(:, counter) = ilut_store_new(:, j)
                        ilut_store_new(:, counter) = ilut_store_old(:, counter)
                    end if
                end do
            else
                do j = 1, new_num_det_states
                    ilut_store_new(NOffFlag, j) = j
                end do

                call sort(ilut_store_new(:, 1:new_num_det_states), ilut_lt_eigen, &
                                                                        ilut_gt_eigen)
                do j = 1, determ_space_cutoff_num(i)
                    ilut_store_old(:,j) = ilut_store_new(:,j)
                end do

                counter = determ_space_cutoff_num(i)
            end if

            old_num_det_states = counter

            write(6,*) old_num_det_states, "states kept."
            call neci_flush(6)

            deallocate(work)
            deallocate(hamiltonian)
            deallocate(determ_ground_state)

        end do

        do i = 1, old_num_det_states
            comp = DetBitLT(ilut_store_old(:, i), ilutHF, NIfD, .false.)
            if (comp == 0) cycle
            call add_basis_state_to_list(ilut_store_old(:, i))
        end do

        deallocate(ilut_store_old)
        deallocate(ilut_store_new)

    end subroutine generate_optimised_core

    subroutine check_if_connected_to_old_space(ilut_store_old, nI, ilut, num_det_states, &
                                                                         counter, connected)

        ! This subroutine loops over all states in the old space and see if any of them is
        ! connected to the state just generated. If so, it also adds one to counter.

        integer(n_int), intent(in) :: ilut_store_old(0:NIfTot, 1000000)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(inout) :: num_det_states, counter
        logical, intent(inout) :: connected
        integer :: nJ(nel)
        real(dp) :: matrix_element
        integer :: i

        connected = .false.

        ! Loop over all states in the old space and see if any of them is connected
        ! to the state just generated.
        do i = 1, num_det_states
            call decode_bit_det(nJ, ilut_store_old(0:NIfD, i))

            matrix_element = get_helement(nI, nJ, ilut, ilut_store_old(0:NIfD, i))

            ! If true, then this state should be added to the new ilut store.
            if (abs(matrix_element) > 0.0_dp) then
                counter = counter + 1 
                connected = .true.
                return
            end if
        end do

    end subroutine check_if_connected_to_old_space

    subroutine generate_sing_doub_determinants()

        type(excit_store), target :: gen_store
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)

        ! This condition tells the enumerating subroutines to initialise the loop.
        ilut(0) = -1

        ! Find the first single excitation.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        ! Subroutine to add this state to the spawned list if on this processor.
        call add_basis_state_to_list(ilut)

        ! When no more basis functions are found, this value is returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)
            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)

            ! If a determinant is returned (if we did not find the final one last time.)
            if (ilut(0) /= -1) call add_basis_state_to_list(ilut)
        end do

        ! Now generate the double excitations...

        ilut(0) = -1
        ! The first double excitation.
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_basis_state_to_list(ilut)

        do while(ilut(0) /= -1)
            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)

            if (ilut(0) /= -1) call add_basis_state_to_list(ilut)
        end do

    end subroutine generate_sing_doub_determinants

    subroutine generate_sing_doub_csfs()

        type(excit_store), target :: gen_store
        integer(n_int) :: ilutHF_loc(0:NIfTot), ilut(0:NIfTot)
        integer :: HFDet_loc(nel), nI(nel)
        integer :: exflag
        logical :: first_loop

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

            call generate_all_csfs_from_orb_config(ilut, nI)

        end do

    end subroutine generate_sing_doub_csfs

    subroutine generate_cas()

        type(BasisFN) :: CASSym
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int), allocatable, dimension(:,:) :: ilut_store
        integer :: HFdet_loc(nel)
        integer :: num_active_orbs, elec, nCASDet, i, j, counter, comp
        integer, allocatable :: CASBrr(:), CASRef(:)
        integer, pointer :: CASDets(:,:) => null()

        ! This option should be true. It tells the subroutine gndts to only consider states
        ! with an Ms value in the correct spin subspace.
        if (.not. tSpn) call stop_all("generate_cas", "tSpn is not set to true.")

        ! The total number of orbitals in the active space:
        num_active_orbs = OccCASorbs + VirtCASorbs
        allocate(CASBrr(1:num_active_orbs))
        allocate(CASRef(1:OccCasOrbs))
        do i = 1, num_active_orbs
            ! Run through the cas space, and create an array which will map these orbtials to the
            ! orbitals they actually represent.
            ! i.e. CASBRR(1) will store the lowest energy orbital in the CAS space and
            ! CASBRR(num_active_orbs) will store the highest energy orbital in the CAS space.
            CASBrr(i) = BRR(i + (nel - OccCasorbs))
        end do

        ! Create a bit mask which has 1's in the bits which represent active orbitals and 0's in
        ! all other orbitals.
        allocate(cas_bitmask(0:NIfD))
        allocate(cas_not_bitmask(0:NIfD))
        cas_bitmask = 0
        do i = 1, num_active_orbs
            set_orb(cas_bitmask, CASBRR(i))
        end do
        ! Create a bit mask which has 0's in the bits which represent active orbitals and 1's in
        ! all other orbitals.
        cas_not_bitmask = not(cas_bitmask)

        ! For Stot /= 0, the HF state will be a CSF. For the purpose of generating all spatial
        ! orbitals, we just want a determinant, so use a state without the CSF information.
        HFdet_loc = iand(HFDet, csf_orbital_mask)

        elec = 1
        do i = nel-OccCasOrbs+1, nel
            ! CASRef(elec) will store the orbital number of the electron elec in the reference
            ! state, HFDet. elec runs from 1 to the number of electrons in the active space.
            CASRef(elec) = HFDet_loc(i)
            elec = elec + 1
        end do

        call GetSym(CASRef, OccCASOrbs, G1, nBasisMax, CASSym)

        ! First, generate all excitations so we know how many there are, to allocate the arrays.
        call gndts(OccCASorbs, num_active_orbs, CASBrr, nBasisMax, CASDets, &
                              .true., G1, tSpn, LMS, .true., CASSym, nCASDet, CASRef)

        if (nCASDet == 0) call stop_all("generate_cas","No CAS determinants found.")

        ! Now allocate the array CASDets. CASDets(:,i) will store the orbitals in the active space
        ! which are occupied in the i'th determinant generated.
        allocate(CASDets(OccCASorbs, nCASDet))
        CASDets(:,:) = 0

        if (tCASCore) then
            allocate(ilut_store(nCASDet-1, 0:NIfTot))
            ilut_store = 0
            counter = 1
        end if

        ! Now fill up CASDets...
        call gndts(OccCASorbs, num_active_orbs, CASBrr, nBasisMax, CASDets, &
                              .false., G1, tSpn, LMS, .true., CASSym, nCASDet, CASRef)

        do i = 1, nCASDet
            ! First, create the bitstring representing this state:
            ! Start from the HF determinant and apply cas_not_bitmask to clear all active space
            ! orbitals.
            ilut(0:NIfD) = iand(ilutHF(0:NIfD), cas_not_bitmask)
            ! Then loop through the occupied orbitals in the active space, stored in CASDets(:,i),
            ! and set the corresponding bits.
            do j = 1, OccCASorbs
                set_orb(ilut, CASDets(j,i))
            end do

            ! The function gndts always outputs the reference state. This is already included, so
            ! we want to cycle when we get to this state to ignore it.
            ! comp will be 0 if ilut and ilutHF are the same.
            comp = DetBitLT(ilut, ilutHF, NIfD, .false.)
            if (comp == 0) cycle

            if (tDeterminantCore) then
                ! Now that we have fully generated the determinant, add it to the main list.
                call add_basis_state_to_list(ilut)
            else if (tCSFCore) then
                ilut_store(counter, 0:NifD) = ilut(0:NifD)
                counter = counter + 1
            end if

        end do

        if (tCSFCore) call create_CAS_csfs(ilut_store, nCASDet-1)

    end subroutine generate_cas

    subroutine create_CAS_csfs(ilut_store, num_dets)

        integer, intent(in) :: num_dets
        integer(n_int), intent(inout) :: ilut_store(num_dets, 0:NIfTot)
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

                if (is_alpha(orb1) .and. IsNotOcc(ilut, orb2) ) then
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

            call generate_all_csfs_from_orb_config(ilut_store(i,:), nI)
        end do

    end subroutine create_CAS_csfs

    subroutine generate_all_csfs_from_orb_config(ilut, nI)

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(inout) :: nI(nel)
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
            ! Finally add the CSF to the CurrentDets list.
            call add_basis_state_to_list(ilut, nI)
        end do

    end subroutine generate_all_csfs_from_orb_config

    pure function ilut_lt_eigen(ilutI, ilutJ) result (bLt)

        integer(n_int), intent(in) :: ilutI(0:), ilutJ(0:)
        logical :: bLt

        bLt = .false.
        
        if (ilutI(NOffFlag) == 0) then
            bLt = .true.
            return
        end if

        bLt = abs(determ_ground_state(ilutI(NOffFlag))) > &
                 abs(determ_ground_state(ilutJ(NOffFlag)))

    end function

    pure function ilut_gt_eigen(ilutI, ilutJ) result (bGt)

        integer(n_int), intent(in) :: ilutI(0:), ilutJ(0:)
        logical :: bGt

        bGt = .false.

        if (ilutI(NOffFlag) == 0) then
            bGt = .true.
            return
        end if

        bGt = abs(determ_ground_state(ilutI(NOffFlag))) < &
                 abs(determ_ground_state(ilutJ(NOffFlag)))

    end function

    subroutine add_basis_state_to_list(ilut, nI_in)

        ! This subroutine, called from init_semi_stochastic, takes a bitstring, finds
        ! the nI representation, decides if the basis state lives on this processor
        ! and, if so, adds it to the spawned list. It also increases the vector sizes
        ! being calculated.

        ! In: ilut - The determinant in a bitstring form.
        ! In (optional) : nI_in - A list of the occupied orbitals in the determinant.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, optional :: nI_in(nel)
        integer :: nI(nel)
        integer(n_int) :: flags
        integer :: proc
        real(dp) :: sgn(lenof_sign)

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

        ! Increase the size of the deterministic space on the correct processor.
        deterministic_proc_sizes(proc) = deterministic_proc_sizes(proc) + 1

        ! If this determinant belongs to this processor, add it to the main list.
        if (proc == iProcIndex) call encode_bit_rep(CurrentDets(:, &
                       deterministic_proc_sizes(iProcIndex)), ilut(0:nIfDBO), sgn, flags)

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
            ilut_temp1 = iand(cas_not_bitmask, ilut(0:NIfD))

            ! The state ilut, but with all active space orbitals and all orbitals not occupied
            ! in the HF determinant unoccupied.
            ilut_temp2 = iand(ilut_temp1, ilutHF(0:NIfD))
            ! All these orbitals should be occupied if the state is in the deterministic space.
            bits_set = sum(count_set_bits(ilut_temp2))
            if (bits_set /= nel - occCASorbs) return

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

    subroutine deterministic_projection()

        ! This subroutine gathers together the partial_determ_vectors from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        integer :: i, info

        call MPIAllGatherV(partial_determ_vector, full_determ_vector, deterministic_proc_sizes, &
                            deterministic_proc_indices)

        ! This function performs y := alpha*A*x + beta*y
        ! N specifies not to use the transpose of A.
        ! deterministic_proc_sizes(iProcIndex) is the number of rows in A.
        ! determ_space_size is the number of columns of A.
        ! alpha = -1.0_dp.
        ! A = core_hamiltonian.
        ! deterministic_proc_sizes(iProcIndex) is the first dimension of A.
        ! input x = full_determ_vector.
        ! 1 is the increment of the elements of x.
        ! beta = 0.0_dp.
        ! output y = partial_determ_vector.
        ! 1 is the incremenet of the elements of y.
        call dgemv('N', &
                   deterministic_proc_sizes(iProcIndex), &
                   determ_space_size, &
                   -1.0_dp, &
                   core_hamiltonian, &
                   deterministic_proc_sizes(iProcIndex), &
                   full_determ_vector, &
                   1, &
                   0.0_dp, &
                   partial_determ_vector, &
                   1)

        ! Now add shift*full_determ_vector, to account for the shift, not stored in core_hamiltonian.
        partial_determ_vector = partial_determ_vector + &
           DiagSft * full_determ_vector(deterministic_proc_indices(iProcIndex)+1:&
             deterministic_proc_indices(iProcIndex)+deterministic_proc_sizes(iProcIndex))

        ! Now multiply the vector by tau to get the final projected vector.
        partial_determ_vector = partial_determ_vector * tau

    end subroutine deterministic_projection

end module semi_stochastic
