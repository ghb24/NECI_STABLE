#include  "macros.h"

module semi_stoch_gen

    use bit_rep_data, only: flag_deterministic, nIfDBO, nOffY, nIfY, NIfD, &
                            flag_is_initiator, flag_determ_parent, NOffSgn, NIfTot
    use bit_reps, only: decode_bit_det, encode_bit_rep, set_flag, extract_sign, &
                        clr_flag
    use CalcData, only: tSortDetermToTop, tTruncInitiator, &
                        tStartCAS, tReadPops, tRegenDiagHels
    use csf, only: csf_get_yamas, get_num_csfs, get_csf_bit_yama, csf_apply_yama
    use csf_data, only: csf_orbital_mask
    use constants
    use DetBitOps, only: ilut_lt, ilut_gt, count_open_orbs, DetBitLT, IsAllowedHPHF
    use DeterminantData, only: write_det
    use enumerate_excitations
    use FciMCData, only: HFDet, ilutHF, iHFProc, CurrentDets, determ_proc_sizes, &
                         determ_proc_indices, full_determ_vector, partial_determ_vector, &
                         core_hamiltonian, determ_space_size, TotWalkers, TotWalkersOld, &
                         indices_of_determ_states, SpawnedParts, CoreTag, FDetermTag, &
                         PDetermTag, trial_space, trial_space_size
    use gndts_mod, only: gndts
    use hash, only: DetermineDetNode
    use LoggingData, only: tWriteCore
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIArg, MPIAllGatherV, &
                             MPIAllGather, MPIScatter, MPIScatterV
    use ParallelHelper, only: root
    use ras
    use semi_stoch_procs
    use sort_mod, only: sort
    use sparse_hamil
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

        ! Initialise the semi-stochastic information. This includes enumerating a list of all
        ! determinants or CSFs in the deterministic space and calculating and storing the resulting
        ! Hamiltonian matrix elements. The lists which will store the walker amplitude vectors in
        ! the deterministic space are also allocated.

        integer :: i, j, IC, ierr, comp
        integer :: nI(nel)
        integer(MPIArg) :: mpi_temp
        character (len=*), parameter :: t_r = "init_semi_stochastic"

        write(6,'()')
        write(6,'(a56)') "============ Semi-stochastic initialisation ============"
        call neci_flush(6)

        allocate(determ_proc_sizes(0:nProcessors-1))
        allocate(determ_proc_indices(0:nProcessors-1))
        determ_proc_sizes = 0
        determ_proc_indices = 0

        if (.not. (tStartCAS .or. tPopsCore .or. tDoublesCore .or. tCASCore .or. tRASCore .or. &
                   tOptimisedCore .or. tLowECore .or. tReadCore)) then
            call stop_all("init_semi_stochastic", "You have not selected a semi-stochastic core &
                          &space to use.")
        end if
        if (.not. (tDeterminantCore .or. tCSFCore .or. tStartCAS .or. tPopsCore .or. tReadCore)) then
            call warning_neci("init_semi_stochastic", "You have not selected to use either &
                              &determinants or CSFs for the deterministic space. Determinants &
                              &will be used.")
            tDeterminantCore = .true.
        end if

        ! Call the enumerating subroutines to create all excitations and add these states to
        ! SpawnedParts on the correct processor. As they do this, they count the size of the
        ! deterministic space (on their own processor only).
        write(6,'(a37)') "Generating the deterministic space..."
        call neci_flush(6)
        call generate_space()

        ! So that all procs store the size of the deterministic spaces on all procs.
        mpi_temp = determ_proc_sizes(iProcIndex)
        call MPIAllGather(mpi_temp, determ_proc_sizes, ierr)

        determ_space_size = sum(determ_proc_sizes)

        write(6,'(a34,1X,i8)') "Total size of deterministic space:", determ_space_size
        write(6,'(a46,1X,i8)') "Size of deterministic space on this processor:", &
                    determ_proc_sizes(iProcIndex)
        call neci_flush(6)

        ! Allocate the vectors to store the psip amplitudes and the deterministic Hamiltonian.
        allocate(full_determ_vector(determ_space_size), stat=ierr)
        call LogMemAlloc('full_determ_vector', int(determ_space_size,sizeof_int), 8, t_r, &
                         FDetermTag, ierr)
        allocate(partial_determ_vector(determ_proc_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('partial_determ_vector', int(determ_proc_sizes(iProcIndex), &
                         sizeof_int), 8, t_r, PDetermTag, ierr)
        allocate(core_hamiltonian(determ_proc_sizes(iProcIndex), determ_space_size), stat=ierr)
        call LogMemAlloc('core_hamiltonian', int(determ_space_size*&
                         &determ_proc_sizes(iProcIndex),sizeof_int), 8, t_r, CoreTag, ierr)

        full_determ_vector = 0.0_dp
        partial_determ_vector = 0.0_dp

        ! This array will hold the positions of the deterministic states in CurrentDets.
        allocate(indices_of_determ_states(determ_proc_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('indices_of_determ_states', int(determ_proc_sizes(iProcIndex), &
                         sizeof_int), bytes_int, t_r, FDetermTag, ierr)

        ! Calculate the indices in the full vector at which the various processors take over, relative
        ! to the first index position in the vector (i.e. the array disps in MPI routines).
        determ_proc_indices(0) = 0
        do i = 1, nProcessors-1
            determ_proc_indices(i) = sum(determ_proc_sizes(:i-1))
        end do

        ! Sort the states to the order that is kept throughout the simulation.
        call sort(SpawnedParts(:,1:determ_proc_sizes(iProcIndex)), ilut_lt, ilut_gt)

        ! Do a check that no states are in the deterministic space twice. The list is sorted
        ! already so simply check states next to each other in the list.
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

        write(6,'(a56)') "Generating the Hamiltonian in the deterministic space..."
        call neci_flush(6)
        call calculate_det_hamiltonian_normal()

        ! Finally, move the states to CurrentDets.
        call add_semistoch_states_to_currentdets()
        ! If starting from a popsfile then CurrentH won't have been filled in yet.
        if ((.not. tRegenDiagHels) .and. tReadPops) call fill_in_CurrentH()
        SpawnedParts = 0
        TotWalkersOld = TotWalkers

        if (tWriteCore) call write_core_space()

        write(6,'(a40)') "Semi-stochastic initialisation complete."
        call neci_flush(6)

    end subroutine init_semi_stochastic

    subroutine generate_space()

        ! A wrapper to call the correct generating routine.

        integer :: i, space_size, ierr

        ! Choose the correct generating routine.
        if (tStartCAS) then
            do i = 1, TotWalkers
                call set_flag(CurrentDets(:, i), flag_deterministic)
                call set_flag(CurrentDets(:, i), flag_is_initiator(1))
                call set_flag(CurrentDets(:, i), flag_is_initiator(2))
            end do
            call MPIAllGather(int(TotWalkers, MPIArg), determ_proc_sizes, ierr)
        else if (tPopsCore) then
            call generate_space_from_pops(called_from_semistoch)
        else if (tReadCore) then
            call generate_space_from_file(called_from_semistoch)
        else if (tDeterminantCore) then
            if (tDoublesCore) then
                call generate_sing_doub_determinants(called_from_semistoch)
            else if (tCASCore) then
                call generate_cas(called_from_semistoch)
            else if (tRASCore) then
                call generate_ras(called_from_semistoch)
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
            else if (tCASCore) then
                call stop_all("init_semi_stochastic", "Cannot use a RAS core space with &
                              &CSFs.")
            else if (tOptimisedCore) then
                call stop_all("init_semi_stochastic", "Optimised core space with CSFs is not &
                              &currently implemented.")
            else if (tLowECore) then
                call stop_all("init_semi_stochastic", "Low energy core space with CSFs is not &
                              &currently implemented.")
            end if
        end if

        ! If requested, remove high energy orbitals so that the space size is below some max.
        if (tLimitDetermSpace) then
            space_size = determ_proc_sizes(iProcIndex)
            call remove_high_energy_orbs(SpawnedParts(:, 1:space_size), space_size, &
                                           max_determ_size, .true.)
            determ_proc_sizes(iProcIndex) = space_size
            tSortDetermToTop = .false.
        end if

    end subroutine generate_space

    subroutine add_state_to_space(ilut, called_from, nI_in)

        ! This subroutine, when called from the semi-stochastic generation code,
        ! takes a state, decides if it lives on this processor and, if so, adds it to
        ! SpawnedParts. It also adds one to the deterministic space size. on this proc.

        ! When called from the trial wavefunction generation code, only the root
        ! processor accesses this code, and each state is added to the list of iluts,
        ! trial_space.

        ! In: ilut - The determinant in a bitstring form.
        ! In: called_from - Integer to specify which part of the code this routine was
        !     called from, and hence which space this state should be added to.
        ! In (optional) : nI_in - A list of the occupied orbitals in the determinant.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: called_from
        integer, optional :: nI_in(nel)
        integer :: nI(nel)
        integer :: flags, proc, comp
        real(dp) :: sgn(lenof_sign)

        ! If using HPHFs then only allow the correct HPHFs to be added to the list.
        if (tHPHF) then
            if (.not. IsAllowedHPHF(ilut(0:NIfD))) return
        end if

        ! Find the nI representation of determinant.
        if (present(nI_in)) then
            nI = nI_in
        else
            call decode_bit_det(nI, ilut)
        end if

        proc = DetermineDetNode(nI,0)

        if (.not. (proc == iProcIndex)) return

        if (called_from == called_from_semistoch) then

            sgn = 0.0_dp
            ! Flag to specify that these basis states are in the deterministic space.
            flags = 0
            flags = ibset(flags, flag_deterministic)
            if (tTruncInitiator) then
                flags = ibset(flags, flag_is_initiator(1))
                flags = ibset(flags, flag_is_initiator(2))
            end if

            ! Keep track of the size of the deterministic space on this processor.
            determ_proc_sizes(proc) = determ_proc_sizes(proc) + 1

            ! Add the state to the space.
            call encode_bit_rep(SpawnedParts(0:NIfTot, determ_proc_sizes(iProcIndex)), &
                                ilut(0:nIfDBO), sgn, flags)

        else if (called_from == called_from_trial) then

            ! Keep track of the size of the trial space on this processor.
            trial_space_size = trial_space_size + 1

            trial_space(0:NIfTot, trial_space_size) = ilut(0:NIfTot)

        end if

    end subroutine add_state_to_space

    subroutine generate_sing_doub_determinants(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        type(excit_store), target :: gen_store
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: space_size
        integer :: i, ierr

        write(6,*) "ilutHF:", ilutHF
 
        ! Start by adding the HF state.
        call add_state_to_space(ilutHF, called_from)

        ! This condition tells the enumerating subroutines to initialise the loop.
        ilut(0) = -1
        ! Find the first single excitation.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_state_to_space(ilut, called_from)

        ! When no more basis functions are found, this value is returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)
            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)

            ! If a determinant is returned (if we did not find the final one last time.)
            if (ilut(0) /= -1) call add_state_to_space(ilut, called_from)
        end do

        ! Now generate the double excitations...

        ilut(0) = -1
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_state_to_space(ilut, called_from)

        do while(ilut(0) /= -1)
            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)

            if (ilut(0) /= -1) call add_state_to_space(ilut, called_from)
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
        ! generating all spatial orbitals we just want a determinant, so use a
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

            call decode_bit_det(nI, ilut)

            call generate_all_csfs_from_orb_config(ilut, nI, called_from)

        end do

    end subroutine generate_sing_doub_csfs

    subroutine generate_ras(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        type(ras_class_data), allocatable, dimension(:) :: core_classes
        integer(n_int), allocatable, dimension(:,:) :: ilut_list
        integer :: nI(nel)
        integer :: space_size, i

        tot_nelec = nel/2
        tot_norbs = nbasis/2

        ! Do a check that the RAS parameters are possible.
        if (core_ras%size_1+core_ras%size_2+core_ras%size_3 /= tot_norbs .or. &
            core_ras%min_1 > core_ras%size_1*2 .or. &
            core_ras%max_3 > core_ras%size_3*2) &
            call stop_all("generate_ras", "RAS parameters are not possible.")

        if (mod(nel, 2) /= 0) call stop_all("generate_ras", "RAS-core only implmented for &
                                            & closed shell molecules.")

        ! Create bitmasks, used in check_if_in_determ_space.
        allocate(core_ras1_bitmask(0:NIfD))
        core_ras1_bitmask = 0
        do i = 1, core_ras%size_1*2
            set_orb(core_ras1_bitmask, Brr(i))
        end do

        allocate(core_ras3_bitmask(0:NIfD))
        core_ras3_bitmask = 0
        do i = (core_ras%size_1+core_ras%size_2)*2+1, nbasis
            set_orb(core_ras3_bitmask, Brr(i))
        end do

        call initialise_ras_space(core_ras, core_classes)

        call find_ras_size(core_ras, core_classes, space_size)

        allocate(ilut_list(0:NIfTot, space_size))

        call generate_entire_ras_space(core_ras, core_classes, space_size, ilut_list)

        do i = 1, space_size
            call add_state_to_space(ilut_list(:,i), called_from)
        end do

        deallocate(core_classes)
        deallocate(ilut_list)

    end subroutine generate_ras

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
        character (len=*), parameter :: t_r = "generate_cas"

        if (called_from == called_from_semistoch) then
            OccOrbs = OccDetermCASOrbs
            VirtOrbs = VirtDetermCASOrbs
        elseif (called_from == called_from_trial) then
            OccOrbs = OccTrialCASOrbs
            VirtOrbs = VirtTrialCASOrbs
        end if

        ! Start by adding the HF state.
        call add_state_to_space(ilutHF, called_from)

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
        call LogMemAlloc("CASDets", OccOrbs*nCASDet, 8, t_r, CASDetsTag, ierr)
        CASDets(:,:) = 0

        if (tCASCore) then
            allocate(ilut_store(nCASDet-1, 0:NIfTot), stat=ierr)
            call LogMemAlloc("ilut_store", (nCASDet-1)*(NIfTot+1), size_n_int, t_r, IlutTag, ierr)
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
                call add_state_to_space(ilut, called_from)
            else if (tCSFCore) then
                ilut_store(counter, 0:NifD) = ilut(0:NifD)
                counter = counter + 1
            end if

        end do

        if (called_from == called_from_semistoch) then
            allocate(cas_determ_bitmask(0:NIfD))
            allocate(cas_determ_not_bitmask(0:NIfD))
            cas_determ_bitmask = cas_bitmask
            cas_determ_not_bitmask = cas_not_bitmask
        end if

        deallocate(CASBrr)
        deallocate(CASRef)
        deallocate(CASDets, stat=ierr)
        call LogMemDealloc(t_r, CASDetsTag, ierr)
        deallocate(ilut_store, stat=ierr)
        call LogMemDealloc(t_r, IlutTag, ierr)

    end subroutine generate_cas

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
        use sparse_hamil, only: sparse_matrix_info, sparse_ham, hamil_diag

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
        character (len=*), parameter :: t_r = "generate_optimised_core"

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
            call LogMemAlloc("ilut_store", 1000000*(NIfTot+1), size_n_int, t_r, IlutTag, ierr)
            allocate(temp_space(0:NIfTot, 1000000), stat=ierr)
            call LogMemAlloc("temp_store", 1000000*(NIfTot+1), size_n_int, t_r, TempTag, ierr)
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
                ! Allow for up to 1 million connected states to be created.
                new_num_states = 1000000
                call generate_connected_space(old_num_states, ilut_store(:, 1:old_num_states), &
                                              new_num_states, temp_space(:, 1:1000000))
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

                call deallocate_sparse_ham()
                deallocate(hamil_diag, stat=ierr)
                call LogMemDealloc(t_r, HDiagTag, ierr)

            end do

        end if ! If on root.

        ! At this point, the space has been generated on the root processor. The rest is
        ! just sending the right info to the right processors...

        if (iProcIndex == root) then
            ! Find which core each state belongs to and sort accordingly.
            call sort_space_by_proc(ilut_store, old_num_states, proc_space_sizes)

            ! Create displacement and sendcount arrays for MPIScatterV later:
            sendcounts = proc_space_sizes*(NIfTot+1)
            disps(0) = 0
            do i = 1, nProcessors-1
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

        ! If being called for the trial space, we want the states to be stored in the
        ! trial_space array, not SpawnedParts.
        if (called_from == called_from_trial) then
            trial_space(:, 1:this_proc_size) = SpawnedParts(:, 1:this_proc_size)
            SpawnedParts = 0
        end if

        ! Finally, deallocate arrays.
        if (allocated(space_cutoff_num)) deallocate(space_cutoff_num)
        if (allocated(space_cutoff_amp)) deallocate(space_cutoff_amp)
        if (iProcIndex == root) then
            deallocate(temp_space, stat=ierr)
            call LogMemDealloc(t_r, TempTag, ierr)
            deallocate(ilut_store, stat=ierr)
            call LogMemDealloc(t_r, IlutTag, ierr)
        end if

    end subroutine generate_optimised_core

    subroutine generate_space_from_pops(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        real(dp), allocatable, dimension(:) :: amps_this_proc, amps_all_procs
        integer(MPIArg) :: length_this_proc, total_length
        integer(MPIArg) :: lengths(0:nProcessors-1), disps(0:nProcessors-1)
        integer(n_int) :: temp_ilut(0:NIfTot)
        integer(n_int), dimension(:,:), allocatable :: largest_states
        integer, allocatable, dimension(:) :: procs_list
        integer :: i, ierr, comp, n_pops_keep, min_ind, max_ind, n_states_this_proc

        if (called_from == called_from_semistoch) then
            n_pops_keep = n_core_pops
        else if (called_from == called_from_trial) then
            n_pops_keep = n_trial_pops
        end if

        if (.not. tReadPops) then
            if (called_from == called_from_semistoch) then
                call stop_all("generate_space_from_pops", "NECI must be started &
                               &from a popsfile to use the pops-core option.")
            else if (called_from == called_from_trial) then
                call stop_all("generate_space_from_pops", "NECI must be started &
                               &from a popsfile to use the pops-trial option.")
            end if
        end if

        length_this_proc = min(n_pops_keep, TotWalkers)

        allocate(amps_this_proc(length_this_proc))
        allocate(amps_all_procs(n_pops_keep*nProcessors))
        allocate(procs_list(n_pops_keep*nProcessors))
        allocate(largest_states(0:NIfTot, length_this_proc))

        call MPIAllGather(length_this_proc, lengths, ierr)
        total_length = sum(lengths)
        if (total_length < n_pops_keep) then
            call warning_neci("generate_space_from_pops", "The number of states in &
                              &POPSFILE is less than the number you requested. All states &
                              &will be used.")
            n_pops_keep = total_length
        end if

        disps(0) = 0
        do i = 0, nProcessors-1
            disps(i) = sum(lengths(:i-1))
        end do

        ! Return the most populated states in CurrentDets on *this* processor.
        call return_most_populated_states(int(length_this_proc,sizeof_int), largest_states)

        ! Store the amplitudes in their real form.
        do i = 1, length_this_proc
            call extract_sign(largest_states(:,i), amps_this_proc(i))
            ! Take the negative absolute value for the sort operation later.
            amps_this_proc(i) = -abs(amps_this_proc(i))
        end do

        ! Now we want to combine all the most populated states from each processor to find
        ! how many states to keep from each processor.

        ! Take the top length_this_proc states from each processor.
        call MPIAllGatherV(amps_this_proc(1:length_this_proc), &
                           amps_all_procs(1:total_length), lengths, disps)

        ! Get procs_list to store the processor labels of the amplitudes in amps_all_procs.
        max_ind = 0
        do i = 0, nProcessors-1
            min_ind = max_ind + 1
            max_ind = min_ind + lengths(i) - 1
            procs_list(min_ind:max_ind) = i
        end do

        ! Sort the list from all processors along with the processor labels, and then
        ! count how many states from each processor are in the top n_pops_keep states.
        call sort(amps_all_procs, procs_list)

        n_states_this_proc = 0
        do i = 1, n_pops_keep
            if (procs_list(i) == iProcIndex) n_states_this_proc = n_states_this_proc + 1
        end do

        ! Add the states to the SpawnedParts array so that they can be processed for the
        ! semi-stochastic space in the standard, consistent way.
        temp_ilut = 0
        do i = 1, n_states_this_proc
            ! The states in largest_states are sorted from smallest to largest.
            temp_ilut(0:NIfD) = largest_states(0:NIfD, length_this_proc-i+1)
            call add_state_to_space(temp_ilut, called_from)
        end do

        deallocate(amps_this_proc)
        deallocate(amps_all_procs)
        deallocate(procs_list)
        deallocate(largest_states)

    end subroutine generate_space_from_pops

    subroutine generate_space_from_file(called_from)

        integer, intent(in) :: called_from
        integer :: iunit, stat
        integer(n_int) :: ilut(0:NIfTot)
        logical :: does_exist
        character(10) :: filename

        if (called_from == called_from_semistoch) then
            filename = 'CORESPACE'
        else if (called_from == called_from_trial) then
            filename = 'TRIALSPACE'
        end if

        inquire(file=filename, exist=does_exist)
        if (.not. does_exist) call stop_all("generate_space_from_file", &
                                            "No "//trim(filename)//" file detected.")

        iunit = get_free_unit()
        open(iunit, file=filename, status='old')

        ilut = 0

        do
            read(iunit, *, iostat=stat) ilut(0:NIfDBO)

            ! If the end of the file.
            if (stat < 0) exit

            call add_state_to_space(ilut, called_from)
        end do

    end subroutine generate_space_from_file

    subroutine generate_low_energy_core(called_from)

        ! In: called_from - Integer to specify whether this routine was called from the
        !     the semi-stochastic generation code or the trial vector generation code.

        integer, intent(in) :: called_from
        integer :: i, num_loops, num_states_keep, ierr
        integer :: old_num_states, new_num_states, max_space_size, low_e_excit
        logical :: tAllDoubles, tSinglesOnly
        integer(n_int), allocatable, dimension(:,:) :: ilut_store, temp_space
        integer(TagIntType) :: IlutTag, TempTag
        character (len=*), parameter :: t_r = "generate_low_energy_core"

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
        call LogMemAlloc("ilut_store", 1000000*(NIfTot+1), size_n_int, t_r, IlutTag, ierr)
        allocate(temp_space(0:NIfTot, 1000000), stat=ierr)
        call LogMemAlloc("temp_store", 1000000*(NIfTot+1), size_n_int, t_r, TempTag, ierr)
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
            ! Allow for up to 1 million connected states.
            new_num_states = 1000000
            call generate_connected_space(old_num_states, ilut_store(:, 1:old_num_states), &
                                          new_num_states, temp_space(:, 1:1000000), tSinglesOnly)

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

        do i = 1, new_num_states
            call add_state_to_space(ilut_store(:, i), called_from)
        end do

        deallocate(temp_space, stat=ierr)
        call LogMemDealloc(t_r, TempTag, ierr)
        deallocate(ilut_store, stat=ierr)
        call LogMemDealloc(t_r, IlutTag, ierr)

    end subroutine generate_low_energy_core

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
            call add_state_to_space(ilut, called_from, nI)
        end do

    end subroutine generate_all_csfs_from_orb_config

end module semi_stoch_gen
