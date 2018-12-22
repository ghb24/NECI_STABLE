#include "macros.h"

module semi_stoch_gen

    use bit_rep_data, only: nIfDBO, NIfD, NIfTot
    use bit_reps, only: decode_bit_det
    use CalcData
    use constants
    use DetBitOps, only: EncodeBitDet
    use FciMCData, only: HFDet, ilutHF
    use gndts_mod, only: gndts, gndts_all_sym_this_proc
    use LoggingData, only: tWriteCore, tRDMonFly
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg, MPIAllGatherV, &
                             MPIAllGather, MPIScatter, MPIScatterV, MPIBarrier
    use ParallelHelper, only: root
    use semi_stoch_procs
    use sparse_arrays
    use timing_neci

    implicit none

contains

    subroutine init_semi_stochastic(core_in, tStartedFromCoreGround)

        ! Initialise the semi-stochastic information. This includes enumerating a list of all
        ! determinants or CSFs in the deterministic space and calculating and storing the resulting
        ! Hamiltonian matrix elements. The lists which will store the walker amplitude vectors in
        ! the deterministic space are also allocated.

        use DetBitOps, only: ilut_lt, ilut_gt
        use DeterminantData, only: write_det
        use FciMCData, only: determ_sizes, determ_displs, full_determ_vecs, full_determ_vecs_av
        use FciMCData, only: partial_determ_vecs, determ_space_size, determ_space_size_int
        use FciMCData, only: TotWalkers, TotWalkersOld, indices_of_determ_states, SpawnedParts
        use FciMCData, only: FDetermTag, FDetermAvTag, PDetermTag, IDetermTag
        use FciMCData, only: tStartCoreGroundState, iter_data_fciqmc, SemiStoch_Init_Time
        use FciMCData, only: tFillingStochRdmOnFly, core_space
        use load_balance, only: adjust_load_balance
        use load_balance_calcnodes, only: tLoadBalanceBlocks
        use sort_mod, only: sort
        use SystemData, only: nel

        type(subspace_in) :: core_in
        logical, intent(out) :: tStartedFromCoreGround

        integer :: i, j, ierr
        integer :: nI(nel)
        integer(MPIArg) :: mpi_temp
        character (len=*), parameter :: t_r = "init_semi_stochastic"

        ! If we are load balancing, this gets disabled once semi stochastic
        ! has been initialised. Therefore we should do a last-gasp load
        ! adjustment at this point.
        if (tLoadBalanceBlocks .and. .not. tFillingStochRDMOnFly) then
            call adjust_load_balance(iter_data_fciqmc)
        end if

!#ifdef __CMPLX
!        call stop_all(t_r, "Semi-stochastic has not been implemented with complex coefficients.")
!#endif

        call MPIBarrier(ierr, tTimeIn=.false.)

        call set_timer(SemiStoch_Init_Time)
        
        write(6,'(/,12("="),1x,a30,1x,12("="))') "Semi-stochastic initialisation"; call neci_flush(6)

        allocate(determ_sizes(0:nProcessors-1))
        allocate(determ_displs(0:nProcessors-1))
        determ_sizes = 0_MPIArg
        determ_displs = 0_MPIArg

        if (.not. (tStartCAS .or. core_in%tPops .or. core_in%tDoubles .or. core_in%tCAS .or. core_in%tRAS .or. &
                   core_in%tOptimised .or. core_in%tLowE .or. core_in%tRead .or. core_in%tMP1 .or. &
                   core_in%tFCI .or. core_in%tHeisenbergFCI .or. core_in%tHF)) then
            call stop_all("init_semi_stochastic", "You have not selected a semi-stochastic core space to use.")
        end if
        if (.not. tUseRealCoeffs) call stop_all(t_r, "To use semi-stochastic you must also use real coefficients.")

        ! Call the enumerating subroutines to create all excitations and add these states to
        ! SpawnedParts on the correct processor. As they do this, they count the size of the
        ! deterministic space (on their own processor only).
        write(6,'("Generating the deterministic space...")'); call neci_flush(6)
        if (core_in%tApproxSpace) write(6,'(" ... approximately using the factor of",1X,i5)') core_in%nApproxSpace
        call generate_space(core_in)

        ! So that all procs store the size of the deterministic spaces on all procs.
        mpi_temp = determ_sizes(iProcIndex)
        call MPIAllGather(mpi_temp, determ_sizes, ierr)

        determ_space_size = sum(determ_sizes)
        determ_space_size_int = int(determ_space_size,sizeof_int)

        write(6,'("Total size of deterministic space:",1X,i8)') determ_space_size
        write(6,'("Size of deterministic space on this processor:",1X,i8)') determ_sizes(iProcIndex)
        call neci_flush(6)

        ! Allocate the vectors to store the walker amplitudes and the deterministic Hamiltonian.
        allocate(full_determ_vecs(lenof_sign,determ_space_size), stat=ierr)
        call LogMemAlloc('full_determ_vecs', determ_space_size_int*lenof_sign, &
                         8, t_r, FDetermTag, ierr)
        allocate(full_determ_vecs_av(lenof_sign,determ_space_size), stat=ierr)
        call LogMemAlloc('full_determ_vecs_av', determ_space_size_int*lenof_sign, &
                         8, t_r, FDetermAvTag, ierr)
        allocate(partial_determ_vecs(lenof_sign,determ_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('partial_determ_vecs', int(determ_sizes(iProcIndex), &
                         sizeof_int)*lenof_sign, 8, t_r, PDetermTag, ierr)

        full_determ_vecs = 0.0_dp
        full_determ_vecs_av = 0.0_dp
        partial_determ_vecs = 0.0_dp

        ! This array will hold the positions of the deterministic states in CurrentDets.
        allocate(indices_of_determ_states(determ_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('indices_of_determ_states', int(determ_sizes(iProcIndex), &
                         sizeof_int), bytes_int, t_r, IDetermTag, ierr)

        ! Calculate the indices in the full vector at which the various processors take over, relative
        ! to the first index position in the vector (i.e. the array disps in MPI routines).
        determ_displs(0) = 0
        do i = 1, nProcessors-1
            determ_displs(i) = determ_displs(i-1) + determ_sizes(i-1)
        end do

        call sort(spawnedparts(0:NIfTot,1:determ_sizes(iprocindex)), ilut_lt, ilut_gt)

        ! Do a check that no states are in the deterministic space twice. The list is sorted
        ! already so simply check states next to each other in the list.
        do i = 2, determ_sizes(iProcIndex)
            if (all(SpawnedParts(0:NIfDBO, i-1) == SpawnedParts(0:NIfDBO, i))) then
                call decode_bit_det(nI, SpawnedParts(:,i))
                write(6,'("State found twice:")')
                write(6,*) SpawnedParts(:,i)
                call write_det(6, nI, .true.)
                call stop_all("init_semi_stochastic", &
                    "The same state has been found twice in the deterministic space.")
            end if
        end do

        ! Store every core determinant from all processors on all processors, in core_space.
        call store_whole_core_space()
        ! Create the hash table to address the core determinants.
        call initialise_core_hash_table(core_space, determ_space_size, core_ht)

        if (tWriteCore) call write_core_space()

        write(6,'("Generating the Hamiltonian in the deterministic space...")'); call neci_flush(6)
        call calc_determ_hamil_sparse()

        if (tRDMonFly) call generate_core_connections()

        ! Move the states to CurrentDets.
        call add_core_states_currentdet_hash()

        ! If using a trial wavefunction, and that initialisation has already
        ! been performed, then the current_trial_amps array needs correcting
        ! after the core states were added and sorted into CurrentDets.
        call reinit_current_trial_amps()

        ! If starting from a popsfile then global_determinant_data will not
        ! have been initialised, or if in the middle of a calculation then new
        ! determinants may have been added.
        if (tReadPops .or. (semistoch_shift_iter /= 0)) call fill_in_diag_helements()

        SpawnedParts = 0_n_int
        TotWalkersOld = TotWalkers

        tStartedFromCoreGround = .false.
        if (tStartCoreGroundState .and. (.not. tReadPops) .and. tStaticCore .and. (.not. tTrialInit)) then
            call start_walkers_from_core_ground(tPrintInfo = .true.)
            tStartedFromCoreGround = .true.
        end if

        ! Call MPIBarrier here so that Semistoch_Init_Time will give the
        ! initialisation time for all processors to finish.
        call MPIBarrier(ierr, tTimeIn=.false.)

        call halt_timer(SemiStoch_Init_Time)

        write(6,'("Semi-stochastic initialisation complete.")')
        write(6,'("Total time (seconds) taken for semi-stochastic initialisation:", f9.3, /)') &
           get_total_time(SemiStoch_Init_Time)
        call neci_flush(6)

    end subroutine init_semi_stochastic

    subroutine generate_space(core_in)

        ! A wrapper to call the correct generating routine.

        use bit_rep_data, only: flag_deterministic, flag_initiator, flag_static_init
        use bit_reps, only: set_flag, encode_sign
        use FciMCData, only: determ_sizes, SpawnedParts
        use searching, only: remove_repeated_states
        use SystemData, only: tAllSymSectors

        type(subspace_in) :: core_in

        integer :: space_size, i, ierr, run
        real(dp) :: zero_sign(lenof_sign)
        character (len=*), parameter :: t_r = "generate_space"

        space_size = 0

        ! Call the requested generating routines.
        if (core_in%tHF) call add_state_to_space(ilutHF, SpawnedParts, space_size)
        if (core_in%tPops) call generate_space_most_populated(core_in%npops, & 
                                    core_in%tApproxSpace, core_in%nApproxSpace, SpawnedParts, space_size)
        if (core_in%tRead) call generate_space_from_file(core_in%read_filename, SpawnedParts, space_size)
        if (.not. tCSFCore) then
            if (core_in%tDoubles) call generate_sing_doub_determinants(SpawnedParts, space_size, core_in%tHFConn)
            if (core_in%tCAS) call generate_cas(core_in%occ_cas, core_in%virt_cas, SpawnedParts, space_size)
            if (core_in%tRAS) call generate_ras(core_in%ras, SpawnedParts, space_size)
            if (core_in%tOptimised) call generate_optimised_space(core_in%opt_data, core_in%tLimitSpace, &
                                                             SpawnedParts, space_size, core_in%max_size)
            if (core_in%tMP1) call generate_using_mp1_criterion(core_in%mp1_ndets, SpawnedParts, space_size)
            if (core_in%tFCI) then
                if (tAllSymSectors) then
                    call gndts_all_sym_this_proc(SpawnedParts, .false., space_size)
                else
                    call generate_fci_core(SpawnedParts, space_size)
                end if
            else if (core_in%tHeisenbergFCI) then
                call generate_heisenberg_fci(SpawnedParts, space_size)
            end if
        else if (tCSFCore) then
            if (core_in%tDoubles) then
                call generate_sing_doub_csfs(SpawnedParts, space_size)
            else if (core_in%tCAS) then
                call stop_all("init_semi_stochastic", "CAS core space with CSFs is not &
                              &currently implemented.")
            else if (core_in%tCAS) then
                call stop_all("init_semi_stochastic", "Cannot use a RAS core space with &
                              &CSFs.")
            else if (core_in%tOptimised) then
                call stop_all("init_semi_stochastic", "Optimised core space with CSFs is not &
                              &currently implemented.")
            else if (core_in%tLowE) then
                call stop_all("init_semi_stochastic", "Low energy core space with CSFs is not &
                              &currently implemented.")
            else if (core_in%tMP1) then
                call stop_all("init_semi_stochastic", "The use of the MP1 wave function criterion &
                              &with CSFs is not implemented.")
            end if
        end if

        ! If two different deterministic spaces have been called then there may
        ! be some repeated states. We don't want repeats, so remove them and
        ! update space_size accordingly.
        call remove_repeated_states(SpawnedParts, space_size)

        zero_sign = 0.0_dp
        do i = 1, space_size
            call encode_sign(SpawnedParts(:,i), zero_sign)

            call set_flag(SpawnedParts(:,i), flag_deterministic)
            if (tTruncInitiator) then
                do run = 1, inum_runs
                    call set_flag(SpawnedParts(:,i), get_initiator_flag_by_run(run))
                    call set_flag(CurrentDets(:,i), flag_static_init(run))
                end do
            end if
        end do

        ! Set the deterministic space size for this process.
        determ_sizes(iProcIndex) = int(space_size, MPIArg)

        ! If requested, remove high energy orbitals so that the space size is below some max.
        if (core_in%tLimitSpace) then
            call remove_high_energy_orbs(SpawnedParts(0:NIfTot, 1:space_size), space_size, &
                                           core_in%max_size, .true.)
            determ_sizes(iProcIndex) = int(space_size, MPIArg)
        end if

    end subroutine generate_space

    subroutine add_state_to_space(ilut, ilut_list, space_size, nI_in)

        ! This subroutine, takes a state, decides if it lives on this processor and,
        ! if so, adds it to ilut_list. It also adds one to the space size on this proc.

        ! In: ilut - The determinant in a bitstring form.
        ! In/Out: ilut_list - List of determinants to which ilut will be added.
        ! In/Out: space_size - The number of determinants belonging to this process.
        ! In (optional): nI_in - A list of the occupied orbitals in the determinant.

        use DetBitOps, only: IsAllowedHPHF
        use load_balance_calcnodes, only: DetermineDetNode
        use SystemData, only: nel

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size
        integer, optional, intent(in) :: nI_in(nel)

        integer :: nI(nel)
        integer :: proc

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

        proc = DetermineDetNode(nel, nI, 0)

        if (.not. (proc == iProcIndex)) return

        space_size = space_size + 1
        ilut_list(:, space_size) = 0_n_int
        ilut_list(0:NIfTot, space_size) = ilut(0:NIfTot)

    end subroutine add_state_to_space

    subroutine generate_sing_doub_determinants(ilut_list, space_size, only_keep_conn)

        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use determinants, only: get_helement
        use SymExcit3, only: GenExcitations3
        use SymExcit4, only: GenExcitations4, ExcitGenSessionType
        use SystemData, only: nel, tKPntSym, tReltvy

        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size
        logical, intent(in) :: only_keep_conn

        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: excit(2,2)
        integer :: nsing, ndoub, ex_flag
        logical :: tAllExcitFound, tParity
        HElement_t(dp) :: HEl

        type(ExcitGenSessionType) :: session

        ! Always generate both the single and double excitations.
        ex_flag = 3

        ! Start by adding the HF state.
        call add_state_to_space(ilutHF, ilut_list, space_size)

        if (tKPntSym) then
            call enumerate_sing_doub_kpnt(ex_flag, only_keep_conn, nsing, ndoub, .true., ilut_list, space_size)
        else

            tAllExcitFound = .false.
            excit = 0

            do while(.true.)
                ! Generate the next determinant.
                if (tReltvy) then
                    call GenExcitations4(session, HFDet, nI, ex_flag, excit, tParity, tAllExcitFound, .false.)
                else
                    call GenExcitations3(HFDet, ilutHF, nI, ex_flag, excit, tParity, tAllExcitFound, .false.)
                endif
                if (tAllExcitFound) exit

                call EncodeBitDet(nI, ilut)
                ! If using a deterministic space connected to the Hartree-Fock
                ! then check that this determinant is actually connected to it!
                if (only_keep_conn) then
                    HEl = get_helement(HFDet, nI, ilutHF, ilut)
                    ! [W.D. 15.5.2017:]
                    ! is this still enough, even for Hamiltonians containing
                    ! complex entries?? 
                    ! and why is the cast to real(dp) done here?? 
!                     if (abs(real(HEl,dp)) < 1.e-12_dp) cycle
                    if (abs(HEl) < 1.e-12_dp) cycle
                end if
                call add_state_to_space(ilut, ilut_list, space_size, nI)
            end do
        end if

    end subroutine generate_sing_doub_determinants

    subroutine generate_sing_doub_csfs(ilut_list, space_size)

        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use bit_rep_data, only: nOffY
        use csf_data, only: csf_orbital_mask
        use enumerate_excitations, only: excit_store, enumerate_spatial_excitations
        use SystemData, only: nel, tFixLz

        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        type(excit_store), target :: gen_store
        integer(n_int) :: ilutHF_loc(0:NIfTot), ilut(0:NIfTot)
        integer :: HFDet_loc(nel), nI(nel)
        integer :: ex_flag
        logical :: first_loop

        if (tFixLz) call stop_all("generate_sing_doub_csfs", "The CSF generating routine &
            &does not work when Lz symmetry is applied.")

        ! Setting the first two bits of this flag tells the generating subroutine to
        ! generate both single and double spatial excitations.
        ex_flag = 0
        ex_flag = ibset(ex_flag, 0)
        ex_flag = ibset(ex_flag, 1)

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
            call enumerate_spatial_excitations(ilutHF_loc, HFDet_loc, ilut, ex_flag, gen_store)

            if (ilut(0) == -1) exit

            call decode_bit_det(nI, ilut)

            call gen_all_csfs_from_orb_config(ilut, nI, ilut_list, space_size)

        end do

    end subroutine generate_sing_doub_csfs

    subroutine generate_ras(ras_info, ilut_list, space_size)

        ! In: ras - Parameters for the RAS space.
        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use ras
        use SystemData, only: nel

        type(ras_parameters), intent(inout) :: ras_info
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        type(ras_class_data), allocatable, dimension(:) :: ras_classes
        integer(n_int), allocatable, dimension(:,:) :: temp_list
        integer :: nI(nel)
        integer :: temp_size, i

        tot_nelec = nel/2
        tot_norbs = nbasis/2

        ! Do a check that the RAS parameters are possible.
        if (ras_info%size_1+ras_info%size_2+ras_info%size_3 /= tot_norbs .or. &
            ras_info%min_1 > ras_info%size_1*2 .or. &
            ras_info%max_3 > ras_info%size_3*2) &
            call stop_all("generate_ras", "RAS parameters are not possible.")

        if (mod(nel, 2) /= 0) call stop_all("generate_ras", "RAS space only implmented for &
                                            & closed shell molecules.")

        call initialise_ras_space(ras_info, ras_classes)

        call find_ras_size(ras_info, ras_classes, temp_size)

        allocate(temp_list(0:NIfTot, temp_size))

        call generate_entire_ras_space(ras_info, ras_classes, temp_size, temp_list)

        do i = 1, temp_size
            call add_state_to_space(temp_list(:,i), ilut_list, space_size)
        end do

        deallocate(ras_classes)
        deallocate(temp_list)

    end subroutine generate_ras

    subroutine generate_cas(occ_orbs, virt_orbs, ilut_list, space_size)

        ! In: occ_orbs - The number of electrons in the CAS space.
        ! In: virt_orbs - The number of virtual spin orbitals in the CAS space.
        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use csf_data, only: csf_orbital_mask
        use DetBitOps, only: DetBitLT
        use sort_mod, only: sort
        use sym_mod, only: getsym
        use SystemData, only: nel, tSpn, G1, nBasisMax, LMS, BasisFn, BRR

        integer, intent(in) :: occ_orbs, virt_orbs
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        type(BasisFN) :: CASSym
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int), allocatable, dimension(:,:) :: ilut_store
        integer :: HFdet_loc(nel), iCASDet
        integer :: num_active_orbs, nCASDet, i, j, counter, comp, ierr
        integer, allocatable :: CASBrr(:), CASRef(:)
        integer(n_int) :: cas_bitmask(0:NIfD), cas_not_bitmask(0:NIfD)
        integer, pointer :: CASDets(:,:) => null()
        integer(TagIntType) :: CASDetsTag, IlutTag
        character (len=*), parameter :: t_r = "generate_cas"

        ! Start by adding the HF state.
        call add_state_to_space(ilutHF, ilut_list, space_size)

        ! This option should be true. It tells the subroutine gndts to only consider states
        ! with an Ms value in the correct spin subspace.
        if (.not. tSpn) call stop_all("generate_cas", "tSpn is not set to true.")

        ! The total number of orbitals in the active space:
        num_active_orbs = occ_orbs + virt_orbs
        allocate(CASBrr(1:num_active_orbs))
        allocate(CASRef(1:occ_orbs))
        do i = 1, num_active_orbs
            ! Run through the cas space, and create an array which will map these orbtials to the
            ! orbitals they actually represent.
            ! i.e. CASBRR(1) will store the lowest energy orbital in the CAS space and
            ! CASBRR(num_active_orbs) will store the highest energy orbital in the CAS space.
            CASBrr(i) = BRR(i + (nel - occ_orbs))
        end do

        ! Create a bit mask which has 1's in the bits which represent active orbitals and 0's in
        ! all other orbitals.
        cas_bitmask = 0_n_int
        do i = 1, num_active_orbs
            set_orb(cas_bitmask, CASBrr(i))
        end do
        ! Create a bit mask which has 0's in the bits which represent active orbitals and 1's in
        ! all other orbitals.
        cas_not_bitmask = not(cas_bitmask)

        ! For Stot /= 0, the HF state will be a CSF. For the purpose of generating all spatial
        ! orbitals, we just want a determinant, so use a state without the CSF information.
        HFdet_loc = iand(HFDet, csf_orbital_mask)

        ! CASRef holds the part of the HF determinant in the active space.
        CASRef = CasBRR(1:occ_orbs)
        call sort(CasRef)
        call GetSym(CASRef, occ_orbs, G1, nBasisMax, CASSym)

        ! First, generate all excitations so we know how many there are, to allocate the arrays.
        call gndts(occ_orbs, num_active_orbs, CASBrr, nBasisMax, CASDets, &
                              .true., G1, tSpn, LMS, .true., CASSym, nCASDet, iCASDet)

        if (nCASDet == 0) call stop_all("generate_cas","No CAS determinants found.")

        ! Now allocate the array CASDets. CASDets(:,i) will store the orbitals in the active space
        ! which are occupied in the i'th determinant generated.
        allocate(CASDets(occ_orbs, nCASDet), stat=ierr)
        call LogMemAlloc("CASDets", occ_orbs*nCASDet, 8, t_r, CASDetsTag, ierr)
        CASDets(:,:) = 0

        if (tCSFCore) then
            allocate(ilut_store(nCASDet-1, 0:NIfTot), stat=ierr)
            call LogMemAlloc("ilut_store", (nCASDet-1)*(NIfTot+1), size_n_int, t_r, IlutTag, ierr)
            ilut_store = 0_n_int
            counter = 1
        end if

        ! Now fill up CASDets...
        call gndts(occ_orbs, num_active_orbs, CASBrr, nBasisMax, CASDets, &
                              .false., G1, tSpn, LMS, .true., CASSym, nCASDet, iCASDet)

        do i = 1, nCASDet
            ! First, create the bitstring representing this state:
            ! Start from the HF determinant and apply cas_not_bitmask to clear all active space
            ! orbitals.
            ilut(0:NIfD) = iand(ilutHF(0:NIfD), cas_not_bitmask)
            ! Then loop through the occupied orbitals in the active space, stored in CASDets(:,i),
            ! and set the corresponding bits.
            do j = 1, occ_orbs
                set_orb(ilut, CASDets(j,i))
            end do

            ! The function gndts always outputs the reference state. This is already included, so
            ! we want to cycle when we get to this state to ignore it.
            ! comp will be 0 if ilut and ilutHF are the same.
            comp = DetBitLT(ilut, ilutHF, NIfD)
            if (comp == 0) cycle

            if (tCSFCore) then
                ilut_store(counter, 0:NifD) = ilut(0:NifD)
                counter = counter + 1
            else
                ! Now that we have fully generated the determinant, add it to the main list.
                call add_state_to_space(ilut, ilut_list, space_size)
            end if

        end do

        deallocate(CASBrr)
        deallocate(CASRef)
        deallocate(CASDets, stat=ierr)
        call LogMemDealloc(t_r, CASDetsTag, ierr)
        if (tCSFCore) then
            deallocate(ilut_store, stat=ierr)
            call LogMemDealloc(t_r, IlutTag, ierr)
        end if

    end subroutine generate_cas

    subroutine generate_optimised_space(opt_data, tLimitSpace, ilut_list, space_size, max_space_size)

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

        ! In: opt_data: Derived type containing the parameters for the algorithm to
        !     generate the space.
        ! In: tLimitSpace: If true then, every iteration of generating algorithm, after
        !     all connections have generated, the space is cut off to have a space with
        !     a maximum size of max_space_size. This is done by removing the highest
        !     energy determinants.
        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.
        ! In (optional): max_space_size - Only used if tLimitSpace is true. See
        !     tLimitSpace for an explanation of use.

        use davidson_neci, only: perform_davidson, DestroyDavidsonCalc, DavidsonCalcType
        use hamiltonian_linalg, only: sparse_hamil_type
        use enumerate_excitations, only: generate_connected_space
        use searching, only: remove_repeated_states
        use sort_mod, only: sort
        use sparse_arrays, only: sparse_ham, hamil_diag
        use SystemData, only: nel

        type(opt_space_data), intent(in) :: opt_data
        logical :: tLimitSpace
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size
        integer, optional, intent(in) :: max_space_size

        integer(n_int), allocatable, dimension(:,:) :: ilut_store, temp_space
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: counter, i, j, ierr
        integer :: old_num_states, new_num_states
        integer(MPIArg) :: proc_space_sizes(0:nProcessors-1), disps(0:nProcessors-1), &
                           sendcounts(0:nProcessors-1), recvcount, this_proc_size
        integer(TagIntType) :: IlutTag, TempTag, FinalTag
        character (len=*), parameter :: t_r = "generate_optimised_space"

        type(DavidsonCalcType) :: davidsonCalc

        if (iProcIndex /= root) then
            ! Allocate some space so that the MPIScatterV call does not crash.
            allocate(ilut_store(0:1, 1), stat=ierr)
            call LogMemAlloc("ilut_store", 1000000*(NIfTot+1), size_n_int, t_r, IlutTag, ierr)
        else if (iProcIndex == root) then
            ! Allocate the stores of ilut's that will hold these deterministic states.
            ! For now, assume that we won't have a deterministic space of more than one
            ! million states. Could make this user-specified later.
            allocate(ilut_store(0:NIfTot, 1000000), stat=ierr)
            call LogMemAlloc("ilut_store", 1000000*(NIfTot+1), size_n_int, t_r, IlutTag, ierr)
            allocate(temp_space(0:NIfTot, 1000000), stat=ierr)
            call LogMemAlloc("temp_store", 1000000*(NIfTot+1), size_n_int, t_r, TempTag, ierr)
            ilut_store = 0_n_int
            temp_space = 0_n_int

            ! Put the Hartree-Fock state in the list first.
            ilut_store(0:NIfTot, 1) = ilutHF(0:NIfTot)

            ! old_num_states will hold the number of deterministic states in the current
            ! space. This is just 1 for now, with only the Hartree-Fock.
            old_num_states = 1

            ! Now we start the iterating loop. Find all states which are either a single or
            ! double excitation from each state in the old ilut store, and then see if they
            ! have a non-zero Hamiltonian matrix element with any state in the old ilut store:

            ! Over the total number of iterations.
            do i = 1, opt_data%ngen_loops

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

                call remove_repeated_states(ilut_store, new_num_states)

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
                call perform_davidson(davidsonCalc, sparse_hamil_type, .false.)

                associate( &
                    davidson_eigenvector => davidsonCalc%davidson_eigenvector &
                )

                ! davidson_eigenvector now stores the ground state eigenvector. We want to use the
                ! vector whose components are the absolute values of this state:
                davidson_eigenvector = abs(davidson_eigenvector)
                ! Multiply by -1.0_dp so that the sort operation later is the right way around.
                davidson_eigenvector = -1.0_dp*davidson_eigenvector

                ! Now decide which states to keep for the next iteration. There are two ways of
                ! doing this, as decided by the user. Either all basis states with an amplitude
                ! in the ground state greater than a given value are kept (tAmpCutoff = .true.),
                ! or a number of states to keep is specified and we pick the states with the
                ! biggest amplitudes (tAmpCutoff = .false.).
                if (opt_data%tAmpCutoff) then
                    counter = 0
                    do j = 1, new_num_states
                        if (abs(davidson_eigenvector(j)) > opt_data%cutoff_amps(i)) then
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
                    old_num_states = opt_data%cutoff_nums(i)
                    if (old_num_states > new_num_states) old_num_states = new_num_states
                end if

                write(6,'(i8,1X,a12)') old_num_states, "states kept."
                call neci_flush(6)

                call deallocate_sparse_ham(sparse_ham, 'sparse_ham', SparseHamilTags)
                deallocate(hamil_diag, stat=ierr)
                call LogMemDealloc(t_r, HDiagTag, ierr)

                end associate
            end do

        end if ! If on root.

        ! At this point, the space has been generated on the root processor. The rest is
        ! just sending the right info to the right processors...

        if (iProcIndex == root) then
            ! Find which core each state belongs to and sort accordingly.
            call sort_space_by_proc(ilut_store, old_num_states, proc_space_sizes)

            ! Create displacement and sendcount arrays for MPIScatterV later:
            sendcounts = int(proc_space_sizes*(NIfTot+1),MPIArg)
            disps(0) = 0
            do i = 1, nProcessors-1
                disps(i) = int(disps(i-1) + proc_space_sizes(i-1)*(NIfTot+1),MPIArg)
!               disps(i) = int(sum(proc_space_sizes(0:i-1))*(NIfTot+1),MPIArg)
            end do
        end if

        ! Send the number of states on each processor to the corresponding processor.
        call MPIScatter(proc_space_sizes, this_proc_size, ierr)
        recvcount = int(this_proc_size*(NIfTot+1),MPIArg)

        ! Finally send the actual determinants to the ilut_list array.
        call MPIScatterV(ilut_store, sendcounts, disps, ilut_list(:, space_size+1:space_size+this_proc_size), recvcount, ierr)

        space_size = space_size + int(this_proc_size, sizeof_int)

        ! Finally, deallocate arrays.
        if (allocated(ilut_store)) then
            deallocate(ilut_store, stat=ierr)
            call LogMemDealloc(t_r, IlutTag, ierr)
        end if
        if (iProcIndex == root) then
            deallocate(temp_space, stat=ierr)
            call LogMemDealloc(t_r, TempTag, ierr)
        end if

    end subroutine generate_optimised_space

    subroutine generate_space_most_populated(target_space_size, tApproxSpace, nApproxSpace, ilut_list, space_size)

        ! In: target_space_size - The number of determinants to attempt to keep
        !         from if less determinants are present then use all of them.
        ! In: tApproxSpace - If true then only find *approximately* the best
        !         space, to save memory (although in many cases it will end up
        !         finding the best space).
        ! In: nApproxSpace - factor that defines how many states are kept on each 
        !         process if tApproxSpace is true, 1 =< nApproxSpace =< nProcessors. 
        !         The larger nApproxSpace, the more memory is consumed and the slower  
        !         (but more accurate) the semi-stochastic initialisation is.
        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use bit_reps, only: extract_sign
        use FciMCData, only: TotWalkers

        integer, intent(in) :: target_space_size, nApproxSpace
        logical, intent(in) :: tApproxSpace
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        real(dp), allocatable, dimension(:) :: amps_this_proc, amps_all_procs
        real(dp) :: real_sign(lenof_sign)
        integer(MPIArg) :: length_this_proc, total_length
        integer(MPIArg) :: lengths(0:nProcessors-1), disps(0:nProcessors-1)
        integer(n_int) :: temp_ilut(0:NIfTot)
        integer(n_int), dimension(:,:), allocatable :: largest_states
        integer, allocatable, dimension(:) :: indices_to_keep
        integer :: i, j, ierr, ind, n_pops_keep, min_ind, max_ind, n_states_this_proc
        integer(TagIntType) :: TagA, TagB, TagC, TagD
        character (len=*), parameter :: t_r = "generate_space_most_populated"
        integer :: nzero_dets

        n_pops_keep = target_space_size


        ! Quickly loop through and find the number of determinants with
        ! zero sign.
        nzero_dets = 0
        do i = 1, TotWalkers
            call extract_sign(CurrentDets(:,i), real_sign)
            if (sum(abs(real_sign)) < 1.e-8_dp) nzero_dets = nzero_dets + 1
        end do

        if (tApproxSpace .and. nProcessors > nApproxSpace) then
            ! Look at nApproxSpace  times the number of states that we expect to keep on
            ! this process. This is done instead of sending the best
            ! target_space_size states to all processes, which is often
            ! overkill and uses up too much memory.
            length_this_proc = min( ceiling(real(nApproxSpace*n_pops_keep)/real(nProcessors), MPIArg), &
                                   int(TotWalkers-nzero_dets,MPIArg) )
        else
            length_this_proc = min( int(n_pops_keep, MPIArg), int(TotWalkers-nzero_dets,MPIArg) )
        end if


        call MPIAllGather(length_this_proc, lengths, ierr)
        total_length = sum(lengths)
        if (total_length < n_pops_keep) then
            call warning_neci(t_r, "The number of states in the walker list is less &
                                   &than the number you requested. All states &
                                   &will be used.")
            n_pops_keep = total_length
        end if


        ! Allocate necessary arrays and log the memory used.
        allocate(amps_this_proc(length_this_proc), stat = ierr)
        call LogMemAlloc("amps_this_proc", int(length_this_proc, sizeof_int), 8, t_r, TagA, ierr)
        allocate(amps_all_procs(total_length), stat = ierr)
        call LogMemAlloc("amps_all_procs", int(total_length, sizeof_int), 8, t_r, TagB, ierr)
        allocate(indices_to_keep(n_pops_keep), stat = ierr)
        call LogMemAlloc("indices_to_keep", n_pops_keep, sizeof_int, t_r, TagC, ierr)
        allocate(largest_states(0:NIfTot, length_this_proc), stat = ierr)
        call LogMemAlloc("largest_states", int(length_this_proc,sizeof_int)*(NIfTot+1), &
                         size_n_int, t_r, TagD, ierr)

        disps(0) = 0_MPIArg
        do i = 1, nProcessors-1
            disps(i) = disps(i-1) + lengths(i-1)
        end do
        

        ! Return the most populated states in CurrentDets on *this* processor.
        call return_most_populated_states(int(length_this_proc,sizeof_int), largest_states)

        ! Store the amplitudes in their real form.
        do i = 1, length_this_proc
            call extract_sign(largest_states(:,i), real_sign)
            ! We are interested in the absolute values of the ampltiudes.
            amps_this_proc(i) = sum(abs(real_sign))
        end do

        ! Now we want to combine all the most populated states from each processor to find
        ! how many states to keep from each processor.
        ! Take the top length_this_proc states from each processor.
        call MPIAllGatherV(amps_this_proc(1:length_this_proc), amps_all_procs(1:total_length), lengths, disps)
        ! This routine returns indices_to_keep, which will store the indices in amps_all_procs
        ! of those amplitudes which are among the n_pops_keep largest (but not sorted).
        call return_largest_indices(n_pops_keep, int(total_length, sizeof_int), amps_all_procs, indices_to_keep)

        n_states_this_proc = 0
        do i = 1, n_pops_keep

            ind = indices_to_keep(i)

            ! Find the processor label for this state. The states of the ampltidues in
            ! amps_all_procs are together in blocks, in order of the corresponding processor
            ! label. Hence, if a state has index <= lengths(0) then it is on processor 0.
            do j = 0, nProcessors-1
                ind = ind - lengths(j)
                if (ind <= 0) then
                    ! j now gives the processor label. If the state is on *this* processor,
                    ! update n_states_this_proc.
                    if (j == iProcIndex) n_states_this_proc = n_states_this_proc + 1
                    exit
                end if
            end do

        end do
        ! Add the states to the ilut_list array.
        temp_ilut = 0_n_int
        do i = 1, n_states_this_proc
            ! The states in largest_states are sorted from smallest to largest.
            temp_ilut(0:NIfTot) = largest_states(0:NIfTot, length_this_proc-i+1)
            call add_state_to_space(temp_ilut, ilut_list, space_size)
        end do
        deallocate(amps_this_proc)
        call LogMemDealloc(t_r, TagA, ierr)
        deallocate(amps_all_procs)
        call LogMemDealloc(t_r, TagB, ierr)
        deallocate(indices_to_keep)
        call LogMemDealloc(t_r, TagC, ierr)
        deallocate(largest_states)
        call LogMemDealloc(t_r, TagD, ierr)

    end subroutine generate_space_most_populated

    subroutine generate_space_from_file(filename, ilut_list, space_size)

        ! In: filename - Name of file to read for determinants.
        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use DetBitOps, only: IsAllowedHPHF, spin_sym_ilut
        use util_mod, only: get_free_unit

        character(255), intent(in) :: filename
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        integer :: iunit, stat
        integer(n_int) :: ilut(0:NIfTot), ilut_tmp(0:NIfTot)
        logical :: does_exist

        inquire(file=trim(filename), exist=does_exist)
        if (.not. does_exist) call stop_all("generate_space_from_file", &
                                            "No "//trim(filename)//" file detected.")

        iunit = get_free_unit()
        open(iunit, file=trim(filename), status='old')

        ilut = 0_n_int
        ilut_tmp = 0_n_int

        do
            read(iunit, *, iostat=stat) ilut(0:NIfDBO)

            ! If the end of the file.
            if (stat < 0) exit

            ! If this determinant isn't the correct determinant for the time
            ! reversal symmetry, then try the spin-flipped version.
            if (tHPHF) then
                if (.not. IsAllowedHPHF(ilut(0:NIfD))) then
                    call spin_sym_ilut (ilut(0:NIfD), ilut_tmp(0:NIfD))
                    if (.not. IsAllowedHPHF(ilut_tmp(0:NIfD))) then
                        cycle
                    else
                        ilut(0:NIfD) = ilut_tmp(0:NIfD)
                    end if
                end if
            end if

            call add_state_to_space(ilut, ilut_list, space_size)
        end do

        close(iunit)

    end subroutine generate_space_from_file

    subroutine gen_all_csfs_from_orb_config(ilut, nI, ilut_list, space_size)

        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use bit_rep_data, only: NIfY, NOffY
        use csf, only: csf_get_yamas, get_num_csfs, get_csf_bit_yama, csf_apply_yama
        use DetBitOps, only: count_open_orbs
        use SystemData, only: nel, Stot

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(inout) :: nI(nel)
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

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
            ! Finally add the CSF to ilut_list.
            call add_state_to_space(ilut, ilut_list, space_size, nI)
        end do

    end subroutine gen_all_csfs_from_orb_config

    subroutine generate_using_mp1_criterion(target_ndets, ilut_list, space_size)

        ! In: target_ndets - The number of determinants to keep, unless there are less
        !     singles and doubles than this value, in which case all singles and doubles
        !     will be kept.
        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use DetBitOps, only: IsAllowedHPHF
        use SymExcit3, only: GenExcitations3
        use SystemData, only: nel, tHub, tUEG, tReal, tNoSingExcits
        use util_mod, only: binary_search_real

        integer, intent(in) :: target_ndets
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        integer(n_int), allocatable :: temp_list(:,:)
        real(dp), allocatable :: amp_list(:)
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: ex(2,2), ex_flag, ndets
        integer :: pos, i
        real(dp) :: amp, energy_contrib
        logical :: tAllExcitFound, tParity

        allocate(amp_list(target_ndets))
        allocate(temp_list(0:NIfD, target_ndets))
        amp_list = 0.0_dp
        temp_list = 0_n_int

        ! Should we generate just singles (1), just doubles (2), or both (3)?
        if (tUEG .or. tNoSingExcits) then
            ex_flag = 2
        else if (tHub) then
            if (tReal) then
                ex_flag = 1
            else
                ex_flag = 2
            end if
        else
            ! Generate both the single and double excitations.
            ex_flag = 3
        end if

        tAllExcitFound = .false.
        ! Count the HF determinant.
        ndets = 1
        ex = 0
        ilut = 0_n_int

        ! Start by adding the HF state.
        temp_list(0:NIfD, 1) = ilutHF(0:NIfD)

        amp_list = huge(amp)
        ! Set this amplitude to be the lowest possible number so that it doesn't get removed from
        ! the list in the following selection - we always want to keep the HF determinant.
        amp_list(1) = -huge(amp)

        ! Loop through all connections to the HF determinant and keep the required number which
        ! have the largest MP1 weights.

        do while (.true.)
            call GenExcitations3(HFDet, ilutHF, nI, ex_flag, ex, tParity, tAllExcitFound, .false.)
            ! When no more basis functions are found, this value is returned and the loop is exited.
            if (tAllExcitFound) exit

            call EncodeBitDet(nI, ilut)
            if (tHPHF) then
                if (.not. IsAllowedHPHF(ilut(0:NIfD))) cycle
            end if
            ndets = ndets + 1

            ! If a determinant is returned (if we did not find the final one last time.)
            if (.not. tAllExcitFound) then
                call return_mp1_amp_and_mp2_energy(nI, ilut, ex, tParity, amp, energy_contrib)
                
                pos = binary_search_real(amp_list, -abs(amp), 1.0e-8_dp)

                ! If pos is less then there isn't another determinant with the same amplitude
                ! (which will be common), but -pos specifies where in the list it should be
                ! inserted to keep amp_list in order.
                if (pos < 0) pos = -pos

                if (pos > 0 .and. pos <= target_ndets) then
                    ! Shuffle all less significant determinants down one slot, and throw away the
                    ! previous least significant determinant.
                    temp_list(0:NIfD, pos+1:target_ndets) = temp_list(0:NIfD, pos:target_ndets-1)
                    amp_list(pos+1:target_ndets) = amp_list(pos:target_ndets-1)

                    ! Add in the new ilut and amplitude in the correct position.
                    temp_list(0:NIfD, pos) = ilut(0:NIfD)
                    ! Store the negative absolute value, because the binary search sorts from
                    ! lowest (most negative) to highest.
                    amp_list(pos) = -abs(amp)
                end if

            end if
        end do

        call warning_neci("generate_using_mp1_criterion", &
            "Note that there are less connections to the Hartree-Fock than the requested &
            &size of the space. The space will therefore be smaller than requested, &
            &containing all connections.")

        ! Now that the correct determinants have been selected, add them to the desired space.
        do i = 1, min(ndets, target_ndets)
            ilut(0:NIfD) = temp_list(0:NIfD, i)
            call add_state_to_space(ilut, ilut_list, space_size)
        end do

    end subroutine generate_using_mp1_criterion

    subroutine enumerate_sing_doub_kpnt(ex_flag, only_keep_conn, nSing, nDoub, tStore, ilut_list, space_size)

        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use csf_data, only: csf_orbital_mask
        use neci_intfce
        use SystemData, only: nel, G1, tUseBrillouin, nBasis

        integer, intent(in) :: ex_flag
        logical, intent(in) :: only_keep_conn
        integer, intent(out) :: nSing, nDoub
        logical, intent(in) :: tStore
        integer(n_int), optional, intent(inout) :: ilut_list(0:,:)
        integer, optional, intent(inout) :: space_size

        integer, allocatable :: excit_gen(:)
        integer(n_int) :: ilut(0:NIfTot)
        integer :: iExcit, iMaxExcit, ierr
        integer :: nJ(nel), hfdet_loc(nel), nStore(6), nExcitMemLen(1)
        logical :: tTempUseBrill
        character(*), parameter :: t_r = 'enumerate_doubles_kpnt'
        HElement_t(dp) :: HEl

        ! A quick hack. Count excitations as though we were a determinant.
        ! We could fix this later...
        hfdet_loc = iand(hfdet, csf_orbital_mask)

        nSing = 0
        nDoub = 0
        iMaxExcit = 0
        nStore(1:6) = 0

        ! Use Alex's old excitation generators. However, we have to ensure
        ! that brillouins theorem isn't on!
        if (tUseBrillouin) then
            tTempUseBrill = .true.
            tUseBrillouin = .false.
        else
            tTempUseBrill = .false.
        end if

        call GenSymExcitIt2(HFDet_loc, nel, G1, nBasis, .true., nExcitMemLen, &
                nJ, iMaxExcit, nStore, ex_flag)

        allocate(excit_gen(nExcitMemLen(1)), stat=ierr)
        if (ierr .ne. 0) call Stop_All(t_r, "Problem allocating excitation generator")
        excit_gen = 0

        call GenSymExcitIt2(HFDet_loc, nel, G1, nBasis, .true., excit_gen, nJ, &
                iMaxExcit, nStore, ex_flag)

        do while(.true.)
            call GenSymExcitIt2(HFDet_loc, nel, G1, nBasis, .false., excit_gen, &
                    nJ, iExcit, nStore, ex_flag)

            if (nJ(1).eq.0) exit

            if (tStore) then
                call EncodeBitDet(nJ,ilut)
                ! If using a deterministic space connected to the Hartree-Fock
                ! then check that this determinant is actually connected to it!
                if (only_keep_conn) then
                    HEl = get_helement(HFDet_loc, nJ, ilutHF, ilut)
                    if (abs(real(HEl,dp)) < 1.e-12_dp) cycle
                end if
                call add_state_to_space(ilut, ilut_list, space_size, nJ)
            end if

            if (iExcit.eq.1) then
                nSing = nSing + 1
            else if (iExcit.eq.2) then
                nDoub = nDoub + 1
            else
                call stop_all(t_r, "Trying to generate more than doubles!")
            end if
        end do

        tUseBrillouin = tTempUseBrill

    end subroutine enumerate_sing_doub_kpnt

    subroutine generate_heisenberg_fci(ilut_list, space_size)

        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use SystemData, only: nel, nbasis

        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        integer :: nsites, nup
        integer :: up_spins(nel/2+1)

        nsites = nbasis/2
        nup = nel/2
        call generate_heisenberg_fci_r(1, up_spins, nsites, nup, ilut_list, space_size)

    end subroutine generate_heisenberg_fci

    recursive subroutine generate_heisenberg_fci_r(ispin, up_spins, nsites, nup, ilut_list, space_size)

        use SystemData, only: nel

        integer, value :: ispin
        integer, intent(inout) :: up_spins(nel/2+1)
        integer, intent(in) :: nsites, nup
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        integer :: i, isite, counter, starting_site
        integer :: alpha_ind, beta_ind, pos
        integer(n_int) :: ilut(0:NIfTot)

        starting_site = 1
        if (ispin > 1) starting_site = up_spins(ispin-1) + 1

        do isite = starting_site, nsites
            counter = 1
            ilut = 0_n_int
            up_spins(ispin) = isite
            ! If we're on the last spin.
            if (ispin == nup) then
                do i = 1, nsites
                    ! If this site has been chosen to have an up spin on it.
                    if (i == up_spins(counter)) then
                        alpha_ind = 2*i
                        pos = (alpha_ind - 1)/bits_n_int
                        ilut(pos) = ibset(ilut(pos), mod(alpha_ind-1, bits_n_int))
                        ! Consider the next up spin on the next loop.
                        counter = counter + 1
                    else
                        beta_ind = 2*i-1
                        pos = (beta_ind - 1)/bits_n_int
                        ilut(pos) = ibset(ilut(pos), mod(beta_ind-1, bits_n_int))
                    end if
                end do
                call add_state_to_space(ilut, ilut_list, space_size)
            else
                call generate_heisenberg_fci_r(ispin+1, up_spins, nsites, nup, ilut_list, space_size)
            end if
        end do 
        
    end subroutine generate_heisenberg_fci_r

    subroutine generate_fci_core(ilut_list, space_size)

        ! In/Out: ilut_list - List of determinants generated.
        ! In/Out: space_size - Number of determinants in the generated space.
        !             If ilut_list is not empty on input and you want to keep
        !             the states already in it, then on input space_size should
        !             be equal to the number of states to be kept in ilut_list,
        !             and new states will be added in from space_size+1.
        !             Otherwise, space_size must equal 0 on input.
        !             On output space_size will equal the total number of
        !             generated plus what space_size was on input.

        use SystemData, only: nel, nbasis, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict, tSymSet

        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(inout) :: space_size

        integer, allocatable :: nI_list(:,:)
        integer(n_int) :: ilut(0:NIfTot)
        integer :: proc, temp(1,1), hf_ind, ndets, i
        character(*), parameter :: t_r = "generate_fci_core"

        if (.not. tSymSet) call stop_all(t_r, "To use the 'FCI-CORE' option you must also choose the symmetry sector of &
                                              &space to be generated by using the 'SYM' option.")

        ! Count the total number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)
        allocate(nI_list(nel, ndets))
        ! Generate and store all the determinants in nI_list.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        do i = 1, ndets
            call EncodeBitDet(nI_list(:,i), ilut)
            call add_state_to_space(ilut, ilut_list, space_size, nI_list(:,i))
        end do

    end subroutine generate_fci_core

    subroutine write_most_pop_core_at_end(target_space_size)

        use FciMCData, only: SpawnedParts
        use util_mod, only: get_free_unit

        ! Write the most populated states in CurrentDets to a DETFILE file,
        ! using the routine generate_space_most_populated, which is the same
        ! routine used by the pops-core semi-stochastic input option. So this
        ! routine basically generates a pops-core space, but can be used at the
        ! end of a calculation, rather than at the start.

        integer, intent(in) :: target_space_size
        integer :: i, j, k, ierr, iunit
        integer :: space_size
        logical :: texist
        character(*), parameter :: t_r = "write_most_pop_core_at_end"

        write(6,'(/,"Finding most populated states...")'); call neci_flush(6)

        space_size = 0

        ! Calling this routine will find the most populated states in
        ! CurrentDets and copy them across to SpawnedParts, to the first
        ! space_size slots in it (overwriting anything which was there before,
        ! which presumably won't be needed now).
        call generate_space_most_populated(target_space_size, .false., 0, SpawnedParts, space_size)

        write(6,'("Writing the most populated states to DETFILE...")'); call neci_flush(6)

        iunit = get_free_unit()

        ! Let each process write its states to the file. Each process waits for
        ! the process before it to finish before starting.
        do i = 0, nProcessors-1

            if (iProcIndex == i) then

                if (i == 0) then
                    open(iunit, file='DETFILE', status='replace')
                else
                    inquire(file='DETFILE',exist=texist)
                    if(.not.texist) call stop_all(t_r,'"DETFILE" file cannot be found')
                    open(iunit, file='DETFILE', status='old', position='append')
                end if
                
                do j = 1, space_size
                    do k = 0, NIfDBO
                        write(iunit, '(i24)', advance='no') SpawnedParts(k,j)
                    end do
                    write(iunit, *)
                end do

                close(iunit)

            end if

            call MPIBarrier(ierr)

        end do

    end subroutine write_most_pop_core_at_end

!------------------------------------------------------------------------------------------!

    subroutine refresh_semistochastic_space()

      use FciMCData, only: iter_data_fciqmc
      use semi_stoch_procs, only: end_semistoch

      implicit none

      logical :: tStartedFromCoreGround

      ! The reinitialization of the semistochastic space can affect population
      ! because of stochastic rounds. To log this correctly, set the iter_data to 0 here
      iter_data_fciqmc%nborn = 0.0_dp
      iter_data_fciqmc%nremoved = 0.0_dp
      tStaticCore = .false.

      ! as the important determinants might change over time, this
      ! resets the semistochastic space taking the current population to get a new one
      call end_semistoch()
      ! the flag_deterministic flag has to be cleared from all determinants as it is
      ! assumed that no states have that flag when init_semi_stochastic starts
      call reset_core_space()
      ! Now, generate the new deterministic space
      call init_semi_stochastic(ss_space_in, tStartedFromCoreGround)

      ! Changing the semi-stochastic space can involve some roundings
      ! if determinants with population < realSpawnCutoff stop being 
      ! in the corespace. Then, we need to log these events.
      iter_data_fciqmc%update_growth = iter_data_fciqmc%update_growth + iter_data_fciqmc%nborn &
           - iter_data_fciqmc%nremoved

    end subroutine refresh_semistochastic_space

!------------------------------------------------------------------------------------------!

    subroutine reset_core_space()

      use bit_reps, only: clr_flag
      use FciMCData, only: MaxWalkersPart

      implicit none

      integer :: i
      
      do i = 1, MaxWalkersPart
         call clr_flag(CurrentDets(:,i),flag_deterministic)
      end do
      
    end subroutine reset_core_space

!------------------------------------------------------------------------------------------!

end module semi_stoch_gen
