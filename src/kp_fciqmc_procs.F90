#include "macros.h"
 
module kp_fciqmc_procs
 
    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use bit_rep_data
    use bit_reps, only: decode_bit_det, encode_sign
    use CalcData, only: tTruncInitiator, tStartSinglePart, InitialPart, InitWalkers
    use CalcData, only: tSemiStochastic, tReadPops, tUseRealCoeffs, tau, DiagSft
    use CalcData, only: AvMCExcits
    use constants
    use DetBitOps, only: DetBitEq, EncodeBitDet, IsAllowedHPHF, FindBitExcitLevel
    use Determinants, only: get_helement
    use dSFMT_interface , only: dSFMT_init, genrand_real2_dSFMT
    use FciMCData, only: ilutHF, HFDet, CurrentDets, SpawnedParts, SpawnedParts2, TotWalkers
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, HashIndex, nWalkerHashes
    use FciMCData, only: fcimc_iter_data, ll_node, MaxWalkersPart, tStartCoreGroundState
    use FciMCData, only: tPopsAlreadyRead, tHashWalkerList, CurrentH, determ_proc_sizes
    use FciMCData, only: core_ham_diag, InputDiagSft, Hii, max_spawned_ind, SpawnedPartsKP
    use FciMCData, only: partial_determ_vector, full_determ_vector, determ_proc_indices, HFSym
    use FciMCData, only: TotParts, TotPartsOld, AllTotParts, AllTotPartsOld, iter, MaxSpawned
    use FciMCData, only: SpawnedPartsKP2, SpawnVecKP, SpawnVecKP2, partial_determ_vecs_kp
    use FciMCData, only: full_determ_vecs_kp, determ_space_size, OldAllAvWalkersCyc
    use FciMCParMod, only: create_particle, InitFCIMC_HF, SetupParameters, InitFCIMCCalcPar
    use FciMCParMod, only: init_fcimc_fn_pointers, WriteFciMCStats, WriteFciMCStatsHeader
    use FciMCParMod, only: rezero_iter_stats_each_iter, tSinglePartPhase
    use gndts_mod, only: gndts
    use hash, only: FindWalkerHash, init_hash_table, reset_hash_table, fill_in_hash_table
    use hash, only: DetermineDetNode, remove_node
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias, create_rand_heisenberg_det
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use Parallel_neci, only: MPIBarrier, iProcIndex, MPISum, MPIReduce, nProcessors
    use ParallelHelper, only: root
    use procedure_pointers
    use semi_stoch_procs, only: copy_core_dets_this_proc_to_spawnedparts, fill_in_CurrentH
    use semi_stoch_procs, only: add_core_states_currentdet_hash, start_walkers_from_core_ground
    use sym_mod, only: getsym
    use SystemData, only: nel, nbasis, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict
    use SystemData, only: BasisFn, tHeisenberg, tHPHF
    use util_mod, only: get_free_unit

    implicit none

    type kp_fciqmc_data
        ! If true then a ground-state calculation is being performed.
        logical :: tGround
        ! If true then a finite-temperature calculation is being performed.
        logical :: tFiniteTemp
        ! The number of different initial walker configurations to start
        ! calculations from.
        integer :: nconfigs
        ! The number of simulations to perform for each initial walker
        ! configuration.
        integer :: nrepeats
        ! The number of different Krylov vectors to sample (the number of
        ! vectors which form the Krylov subspace at the end of a calculation).
        integer :: nvecs
        ! The number of iterations to perform *between each Krylov vector being
        ! sampled*. niters(i) holds the number to be performed between the
        ! i-th and (i+1)th vectors. The final element is not set by the user,
        ! it is always set to 1 (because we don't want to go any further
        ! beyond the final Krylov vector).
        integer, allocatable :: niters(:)
        ! For finite-temperature calculations, when creating the inital vector,
        ! this variable specifies how many walkers should be added to each
        ! chosen site.
        real(dp) :: nwalkers_per_site_init
        ! On iterations where the spawned walkers are used to estiamate the
        ! Hamiltonian, this is used instead of AvMCExcits.
        real(dp) :: av_mc_excits_kp

        real(dp), allocatable :: overlap_matrix(:,:)
        real(dp), allocatable :: hamil_matrix(:,:)
    end type

    type(kp_fciqmc_data) :: kp

    integer :: nhashes_kp
    integer :: TotWalkersKP
    integer(n_int), allocatable :: krylov_vecs(:,:)
    type(ll_node), pointer :: krylov_vecs_ht(:) 

    integer :: TotWalkersInit
    integer :: TotPartsInit(lenof_sign)
    integer :: AllTotPartsInit(lenof_sign)
    integer(n_int), allocatable :: init_kp_config(:,:)
    logical :: vary_niters

    ! The total sign length for all Krylov vectors together.
    integer :: lenof_sign_kp

    ! If true then calculate the projected Hamiltonian exactly (useful
    ! for testing only, in practice).
    logical :: tExactHamil
    ! If true, use the spawning from the main FCIQMC iterations to
    ! calculate the projected Hamiltonian.
    logical :: tHamilOnFly

    ! If true then use generate_init_config_this_proc to generate the initial
    ! walker distribution for finite-temperature calculations. This will always
    ! get the request walker population.
    logical :: tInitCorrectNWalkers

    integer :: MaxSpawnedEachProc

    logical :: tUseInitConfigSeeds
    integer, allocatable :: init_config_seeds(:)

contains

    subroutine kp_fciqmc_read_inp()

        use input_neci

        logical :: eof
        character(len=100) :: w
        integer :: i, niters_temp
        character (len=*), parameter :: t_r = "kp_fciqmc_read_inp"

        ! Default values.
        kp%nconfigs = 1
        kp%nrepeats = 1
        kp%nvecs = 0
        niters_temp = 0
        kp%nwalkers_per_site_init = 1.0_dp
        kp%av_mc_excits_kp = 0.0_dp
        kp%tGround = .false.
        kp%tFiniteTemp = .false.
        tExactHamil = .false.
        tHamilOnFly = .false.
        tInitCorrectNWalkers = .false.
        vary_niters = .false.
        tUseInitConfigSeeds = .false.

        read_inp: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)
            case("END-KP-FCIQMC")
                exit read_inp
            case("GROUND-STATE")
                kp%tGround = .true.
            case("FINITE-TEMPERATURE")
                kp%tFiniteTemp = .true.
            case("NUM-INIT-CONFIGS")
                call geti(kp%nconfigs)
            case("NUM-REPEATS-PER-INIT-CONFIG")
                call geti(kp%nrepeats)
            case("NUM-KRYLOV-VECS")
                call geti(kp%nvecs)
                allocate(kp%niters(kp%nvecs))
                kp%niters = 0
            case("NUM-ITERS-BETWEEN-VECS")
                call geti(niters_temp)
            case("NUM-ITERS-BETWEEN-VECS-VARY")
                if (kp%nvecs == 0) call stop_all(t_r, 'Please enter the number of krylov vectors before &
                                                          &specifying the number of iterations between vectors.')
                vary_niters = .true.
                do i = 1, kp%nvecs-1
                    call geti(kp%niters(i))
                end do
                kp%niters(kp%nvecs) = 1
            case("NUM-WALKERS-PER-SITE-INIT")
                call getf(kp%nwalkers_per_site_init)
            case("AVERAGEMCEXCITS-HAMIL")
                call getf(kp%av_mc_excits_kp)
            case("EXACT-HAMIL")
                tExactHamil = .true.
            case("HAMIL-ON-FLY")
                tHamilOnFly = .true.
            case("INIT-CORRECT-WALKER-POP")
                tInitCorrectNWalkers = .true.
            case("INIT-CONFIG-SEEDS")
                if (kp%nconfigs == 0) call stop_all(t_r, 'Please use the num-init-configs option to enter the number of &
                                                          &initial configurations before using the init-vec-seeds option.')
                tUseInitConfigSeeds = .true.
                allocate(init_config_seeds(kp%nconfigs))
                do i = 1, kp%nconfigs
                    call geti(init_config_seeds(i))
                end do
            case default
                call report("Keyword "//trim(w)//" not recognized in kp-fciqmc block", .true.)
            end select
        end do read_inp

        if (kp%nvecs == 0 .and. niters_temp /= 0) then
            ! If nvecs was not set but niters was.
            kp%nvecs = 1
            allocate(kp%niters(kp%nvecs))
            kp%niters = niters_temp
            kp%niters(kp%nvecs) = 1
        else if (kp%nvecs /= 0 .and. niters_temp /= 0 .and. (.not. vary_niters)) then
            ! If both were set (but not through the variable niters option).
            kp%niters = niters_temp
            kp%niters(kp%nvecs) = 1
        else if (kp%nvecs /= 0 .and. niters_temp == 0 .and. (.not. vary_niters)) then
            ! If nvecs was set but niters was not.
            kp%niters = 1
        else if (kp%nvecs == 0) then
            ! If nvecs was not set and neither was niters.
            kp%nvecs = 1
            allocate(kp%niters(kp%nvecs))
            kp%niters = 1
        end if

    end subroutine kp_fciqmc_read_inp

    subroutine init_kp_fciqmc(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: ierr
        character (len=*), parameter :: t_r = "init_kp_fciqmc"

        if (.not. tHashWalkerList) call stop_all('t_r','kp-fciqmc can only be run using &
            &the linscalefcimcalgo option (the linear scaling algorithm).')

        if (.not. tUseRealCoeffs) call stop_all('t_r','kp-fciqmc can only be run using &
            &real coefficients).')

        if (tExactHamil .and. nProcessors /= 1) call stop_all('t_r','The exact-hamil &
            &option can only be used when running with one processor.')

        tPopsAlreadyRead = .false.
        call SetupParameters()
        call InitFCIMCCalcPar()
        call init_fcimc_fn_pointers() 

        call WriteFciMCStatsHeader()

        if (n_int == 4) call stop_all('t_r', 'Use of RealCoefficients does not work with 32 bit &
             &integers due to the use of the transfer operation from dp reals to 64 bit integers.')

        lenof_sign_kp = lenof_sign*kp%nvecs
        NIfLan = NIfDBO + lenof_sign_kp + 1 + NIfFlag

        nhashes_kp = nWalkerHashes
        TotWalkersKP = 0
        allocate(krylov_vecs(0:NIfLan, MaxWalkersPart), stat=ierr)
        krylov_vecs = 0_n_int
        allocate(krylov_vecs_ht(nhashes_kp), stat=ierr)
        call init_hash_table(krylov_vecs_ht)

        allocate(kp%overlap_matrix(kp%nvecs, kp%nvecs), stat=ierr)
        allocate(kp%hamil_matrix(kp%nvecs, kp%nvecs), stat=ierr)

        ! If performing a finite-temperature calculation with more than one run for each initial
        ! configuration, we store this walker configuration so that we can restart from it later.
        if (kp%tFiniteTemp .and. kp%nrepeats > 1) then
            allocate(init_kp_config(0:NIfTot, MaxWalkersPart), stat=ierr)
        end if

        ! If av_mc_excits_kp hasn't been set by the user, just use AvMCExcits.
        if (kp%av_mc_excits_kp == 0.0_dp) kp%av_mc_excits_kp = AvMCExcits

        if (tHamilOnFly) then
            allocate(SpawnVecKP(0:NIfTot,MaxSpawned),stat=ierr)
            SpawnVecKP(:,:) = 0
            SpawnedPartsKP => SpawnVecKP
        else
            allocate(SpawnVecKP(0:NOffSgn+lenof_sign_kp-1,MaxSpawned),stat=ierr)
            allocate(SpawnVecKP2(0:NOffSgn+lenof_sign_kp-1,MaxSpawned),stat=ierr)
            SpawnVecKP(:,:) = 0
            SpawnVecKP2(:,:) = 0
            SpawnedPartsKP => SpawnVecKP
            SpawnedPartsKP2 => SpawnVecKP2
            if (tSemiStochastic) then
                allocate(partial_determ_vecs_kp(lenof_sign_kp,determ_proc_sizes(iProcIndex)), stat=ierr)
                allocate(full_determ_vecs_kp(lenof_sign_kp,determ_space_size), stat=ierr)
            end if
        end if

        MaxSpawnedEachProc = int(0.85*real(MaxSpawned,dp)/nProcessors)

        call MPIBarrier(ierr)

    end subroutine init_kp_fciqmc

    subroutine init_kp_fciqmc_repeat(kp, iconfig, irepeat)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer, intent(in) :: iconfig, irepeat

        call create_initial_config(kp, iconfig, irepeat)

        call reset_hash_table(krylov_vecs_ht)
        krylov_vecs = 0_n_int

        kp%overlap_matrix = 0.0_dp
        kp%hamil_matrix = 0.0_dp
        TotWalkersKP = 0
        iter = 0

        DiagSft = InputDiagSft
        if (tStartSinglePart) then
            OldAllAvWalkersCyc = InitialPart
            AllTotPartsOld=InitialPart
        else
            OldAllAvWalkersCyc = InitWalkers
            AllTotPartsOld=InitWalkers
        end if
        ! Setting this variable to true stops the shift from varying instantly.
        tSinglePartPhase = .true.

    end subroutine init_kp_fciqmc_repeat

    subroutine init_kp_fciqmc_iter(iter_data, determ_index)

        use FciMCData, only: fcimc_excit_gen_store, FreeSlot, iStartFreeSlot, iEndFreeSlot

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(out) :: determ_index

        ! Reset positions to spawn into in the spawning array.
        ValidSpawnedList = InitialSpawnedSlots

        ! Reset the array which holds empty slots in CurrentDets.
        FreeSlot(1:iEndFreeSlot) = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0

        ! Index for counting deterministic states.
        determ_index = 1

        call rezero_iter_stats_each_iter(iter_data)

    end subroutine init_kp_fciqmc_iter

    subroutine create_initial_config(kp, iconfig, irepeat)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer, intent(in) :: iconfig, irepeat
        integer :: DetHash, i, nwalkers, nwalkers_target

        call reset_hash_table(HashIndex)

        if (kp%tGround) then
            ! Put a walker on the Hartree-Fock again, with the requested amplitude.
            call InitFCIMC_HF()
            if (tSemiStochastic) then
                ! core_space stores all core determinants from all processors. Move those on this
                ! processor to SpawnedParts, which add_core_states_currentdet_hash uses.
                call copy_core_dets_this_proc_to_spawnedparts()
                call add_core_states_currentdet_hash()
                SpawnedParts = 0_n_int
                if (tStartCoreGroundState .and. (.not. tReadPops)) &
                    call start_walkers_from_core_ground(tPrintInfo = .false.)
                ! Reset the diagonal Hamiltonian elements.
                CurrentH(1, 1:determ_proc_sizes(iProcIndex)) = core_ham_diag
            end if
        else if (kp%tFiniteTemp) then
            if (irepeat == 1) then
                ! Convert the initial number of walkers to an integer.
                if (tStartSinglePart) then
                    nwalkers_target = ceiling(real(InitialPart,dp)/real(nProcessors,dp))
                else
                    nwalkers_target = ceiling(real(InitWalkers,dp)/real(nProcessors,dp))
                end if
                ! Finally, call the routine to create the walker distribution.
                if (tUseInitConfigSeeds) call dSFMT_init((iProcIndex+1)*init_config_seeds(iconfig))
                if (tInitCorrectNWalkers) then
                    call generate_init_config_this_proc(nwalkers_target, kp%nwalkers_per_site_init, nwalkers)
                else
                    call generate_init_config_basic(nwalkers_target, kp%nwalkers_per_site_init, nwalkers)
                end if
                TotWalkersInit = TotWalkers
                TotPartsInit = TotParts
                AllTotPartsInit = AllTotParts
                ! If starting from this configuration more than once, store it.
                if (kp%nrepeats > 1) init_kp_config(:,1:TotWalkers) = CurrentDets(:,1:TotWalkers)
            else if (irepeat > 1) then
                TotWalkers = TotWalkersInit
                TotParts = TotPartsInit
                TotPartsOld = TotPartsInit
                AllTotParts = AllTotPartsInit
                AllTotPartsOld = AllTotPartsInit

                CurrentDets(:,1:TotWalkers) = init_kp_config(:,1:TotWalkers)
                call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, TotWalkers)
            end if
            call fill_in_CurrentH()
        end if

    end subroutine create_initial_config

    subroutine generate_init_config_basic(nwalkers, nwalkers_per_site_init, ndets)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        integer, intent(in) :: nwalkers
        real(dp) :: nwalkers_per_site_init
        integer, intent(out) :: ndets
        integer :: i, ireplica, excit, nattempts, DetHash
        integer :: nspawns
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_amp, walker_sign(lenof_sign)
        logical :: tInitiatorTemp
        type(fcimc_iter_data) :: unused_data
        integer(n_int), pointer :: PointTemp(:,:)

        ! Turn off the initiator method for the annihilation steps to be used here.
        tInitiatorTemp = tTruncInitiator
        tTruncInitiator = .false.

        ! Set the spawning slots to their starting positions.
        ValidSpawnedList = InitialSpawnedSlots

        ilut = 0_n_int
        nspawns = ceiling(real(nwalkers,dp)/nwalkers_per_site_init)

        do i = 1, nspawns
            ! Generate the determinant (output to ilut).
            if (tHeisenberg) then
                call create_rand_heisenberg_det(ilut)
            else
                call CreateRandomExcitLevDetUnbias(nel, HFDet, ilutHF, ilut, excit, nattempts)
            end if
            call decode_bit_det(nI, ilut)

            ! Choose whether the walker should have a positive or negative amplitude, with
            ! 50% chance of each.
            walker_amp = nwalkers_per_site_init
            r = genrand_real2_dSFMT()
            if (r < 0.5) walker_amp = -1.0_dp*walker_amp

            do ireplica = 1, inum_runs
                walker_sign = 0.0_dp
                walker_sign(ireplica) = walker_amp
                call create_particle(nI, ilut, walker_sign, 0, ireplica)
            end do
        end do

        ! Perform annihilation steps:
        ! Send the walkers to their correct processors. The resulting walkers will be stored in
        ! SpawnedParts2.
        call SendProcNewParts(ndets, tSingleProc = .false.)
        ! CompressSpawnedList works on SpawnedParts, not SpawnedParts2, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp
        call CompressSpawnedList(ndets, unused_data) 

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        TotParts = 0.0_dp
        do i = 1, ndets
            CurrentDets(:,i) = SpawnedParts(:,i)
            walker_sign = transfer(CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, i), walker_sign)
            TotParts = TotParts + abs(walker_sign)
        end do
        TotPartsOld = TotParts

        ! Add the entries into the hash table.
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets)

        call MPIReduce(TotParts, MPI_SUM, AllTotParts)
        AllTotPartsOld = AllTotParts
        TotWalkers = int(ndets, int64)

        if (tSemiStochastic) then
            ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
            ! These routines will do this.
            call copy_core_dets_this_proc_to_spawnedparts()
            call add_core_states_currentdet_hash()
        end if

        ValidSpawnedList = InitialSpawnedSlots
        SpawnedParts = 0_n_int

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

    end subroutine generate_init_config_basic

    subroutine generate_init_config_this_proc(nwalkers, nwalkers_per_site_init, ndets)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        integer, intent(in) :: nwalkers
        real(dp) :: nwalkers_per_site_init
        integer, intent(out) :: ndets
        integer :: proc, excit, nattempts, idet, DetHash, det_ind, nI(nel)
        integer(n_int) :: ilut(0:NIfTot), int_sign(lenof_sign)
        type(ll_node), pointer :: prev, temp_node
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        real(dp) :: new_sign(lenof_sign), r
        logical :: tDetFound

        ilut = 0_n_int
        ndets = 0
        TotParts = 0.0_dp

        do
            ! Generate the determinant (output to ilut).
            if (tHeisenberg) then
                call create_rand_heisenberg_det(ilut)
            else
                call CreateRandomExcitLevDetUnbias(nel, HFDet, ilutHF, ilut, excit, nattempts)
            end if
            call decode_bit_det(nI, ilut)
            proc = DetermineDetNode(nI, 0)
            if (proc /= iProcIndex) cycle

            ! Choose whether the walker should have a positive or negative amplitude, with
            ! 50% chance of each.
            real_sign_1 = nwalkers_per_site_init
            r = genrand_real2_dSFMT()
            if (r < 0.5) real_sign_1 = -1.0_dp*real_sign_1
            int_sign = transfer(real_sign_1, int_sign)

            tDetFound = .false.
            DetHash = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(DetHash)
            prev => null()
            ! If the first element in the list for this hash value has been used.
            if (.not. temp_node%ind == 0) then
                ! Loop over all determinants with this hash value which are already in the list.
                do while (associated(temp_node))
                    if (DetBitEQ(CurrentDets(:,temp_node%ind), ilut, NIfDBO)) then
                        ! This determinant is already in the list.
                        tDetFound = .true.
                        int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, temp_node%ind)
                        real_sign_2 = transfer(int_sign, real_sign_2)
                        new_sign = real_sign_1 + real_sign_2
                        call encode_sign(CurrentDets(:, temp_node%ind), new_sign)
                        ! If the walkers have annihilated completley, remove the determinant.
                        !if (IsUnoccDet(new_sign)) call remove_node(prev, temp_node)
                        TotParts = TotParts - abs(real_sign_2) + abs(new_sign)
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    prev => temp_node
                    temp_node => temp_node%next
                end do

                if (.not. tDetFound) then
                    ! We need to add a new determinant in the next position in the list.
                    ! So create that next position!
                    allocate(prev%next)
                    temp_node => prev%next
                    nullify(temp_node%next)
                end if
            end if

            if (.not. tDetFound) then
                ! A new determiant needs to be added.
                ndets = ndets + 1
                TotParts = TotParts + abs(real_sign_1)
                det_ind = ndets
                temp_node%ind = det_ind
                ! Copy determinant data across.
                CurrentDets(0:NIfDBO, det_ind) = ilut(0:NIfDBO)
                CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, det_ind) = int_sign
                if (tUseFlags) CurrentDets(NOffFlag, det_ind) = 0_n_int

                nullify(temp_node)
                nullify(prev)
            end if

            if (TotParts(1) >= nwalkers) exit

        end do

        ! Remove the nodes of all unoccupied determinants from the hash table.
        !do idet = 1, ndets
        !    int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, idet)
        !    if (.not. IsUnoccDet(int_sign)) cycle
        !    tDetFound = .false.
        !    call decode_bit_det(nI, CurrentDets(:, idet))
        !    DetHash = FindWalkerHash(nI, nWalkerHashes)
        !    temp_node => HashIndex(DetHash)
        !    prev => null()
        !    if (.not. temp_node%ind == 0) then
        !        ! Loop over all determinants with this hash value.
        !        do while (associated(temp_node))
        !            if (temp_node%ind == idet) then
        !                tDetFound = .true.
        !                call remove_node(prev, temp_node)
        !                exit
        !            end if
        !            ! Move on to the next determinant with this hash value.
        !            prev => temp_node
        !            temp_node => temp_node%next
        !        end do
        !    end if
        !    ASSERT(tDetFound)
        !end do

        !write(6,*) "Dets present before:"
        !do idet = 1, ndets
        !    int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, idet)
        !    real_sign_2 = transfer(int_sign, real_sign_2)
        !    if (tUseFlags) then
        !        write(6,'(i7, i12, 4x, f18.7, 4x, f18.7, 4x, l1)') idet, CurrentDets(0, idet), real_sign_2, &
        !            test_flag(CurrentDets(:, idet), flag_deterministic)
        !    else
        !        write(6,'(i7, i12, 4x, f18.7, 4x, f18.7)') idet, CurrentDets(0, idet), real_sign_2
        !    end if
        !end do

        TotPartsOld = TotParts

        call MPIReduce(TotParts, MPI_SUM, AllTotParts)
        AllTotPartsOld = AllTotParts
        TotWalkers = int(ndets, int64)

        if (tSemiStochastic) then
            ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
            ! These routines will do this.
            call copy_core_dets_this_proc_to_spawnedparts()
            call add_core_states_currentdet_hash()
        end if

        SpawnedParts = 0_n_int

    end subroutine generate_init_config_this_proc

    subroutine store_krylov_vec(ivec, nvecs)

        integer, intent(in) ::  ivec, nvecs
        integer :: idet, sign_ind, hdiag_ind, flag_ind, DetHash, det_ind
        integer :: nI(nel)
        integer(n_int) :: temp
        logical :: tDetFound
        type(ll_node), pointer :: temp_node, prev

        ! The index of the first element referring to the sign, for this ivec.
        sign_ind = NIfDBO + lenof_sign*(ivec-1) + 1
        hdiag_ind = NIfDBO + lenof_sign_kp + 1
        if (tUseFlags) flag_ind = NIfDBO + lenof_sign_kp + 2

        ! Loop over all occupied determinants for this new Krylov vector.
        do idet = 1, TotWalkers
            tDetFound = .false.

            call decode_bit_det(nI, CurrentDets(:,idet))
            DetHash = FindWalkerHash(nI, nhashes_kp)
            temp_node => krylov_vecs_ht(DetHash)

            ! If the first element in the list for this hash value has been used.
            if (.not. temp_node%ind == 0) then
                ! Loop over all determinants with this hash value which are already in the list.
                do while (associated(temp_node))
                    if (DetBitEQ(CurrentDets(:,idet), krylov_vecs(:,temp_node%ind), NIfDBO)) then
                        ! This determinant is already in the list.
                        det_ind = temp_node%ind
                        ! Add the sign for the new Krylov vector. The determinant and flag are
                        ! there already.
                        krylov_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = &
                              CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,idet)
                        tDetFound = .true.
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    prev => temp_node
                    temp_node => temp_node%next
                end do

                if (.not. tDetFound) then
                    ! We need to add a new determinant in the next position in the list.
                    ! So create that next position!
                    allocate(prev%next)
                    temp_node => prev%next
                    nullify(temp_node%next)
                end if
            end if

            if (.not. tDetFound) then
                ! A new determiant needs to be added.
                TotWalkersKP = TotWalkersKP + 1
                det_ind = TotWalkersKP
                temp_node%ind = det_ind

                ! Copy determinant data across.
                krylov_vecs(0:NIfDBO,det_ind) = CurrentDets(0:NIfDBO,idet)
                krylov_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,idet)
                krylov_vecs(hdiag_ind,det_ind) = transfer(CurrentH(1,idet), temp)
                if (tUseFlags) krylov_vecs(flag_ind,det_ind) = CurrentDets(NOffFlag,idet)
            end if

            nullify(temp_node)
            nullify(prev)
        end do

    end subroutine store_krylov_vec

    subroutine calc_overlap_matrix_elems(kp, ivec)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec)
        integer(n_int) :: sgn(lenof_sign)
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign)

        associate(s_matrix => kp%overlap_matrix)

            ! Just in case!
            s_matrix(1:ivec, ivec) = 0.0_dp
            s_matrix(ivec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in krylov_vecs, for each Krylov vector.
                ind(jvec) = NIfDBO + lenof_sign*(jvec-1) + 1
            end do

            ! Loop over all determinants in krylov_vecs.
            do idet = 1, TotWalkersKP
                sgn = krylov_vecs(ind(ivec):ind(ivec)+1, idet)
                if (IsUnoccDet(sgn)) cycle
                sign1 = transfer(sgn, sign1)
                ! Loop over all Krylov vectors currently stored.
                do jvec = 1, ivec
                    sgn = krylov_vecs(ind(jvec):ind(jvec)+1, idet)
                    if (IsUnoccDet(sgn)) cycle
                    sign2 = transfer(sgn, sign1)
                    s_matrix(jvec,ivec) = s_matrix(jvec,ivec) + &
                        (sign1(1)*sign2(2) + sign1(2)*sign2(1))/2.0_dp
                end do
            end do

            ! Fill in the lower-half of the overlap matrix.
            do jvec = 1, ivec
                s_matrix(ivec,jvec) = s_matrix(jvec,ivec)
            end do

        end associate

    end subroutine calc_overlap_matrix_elems

    subroutine calc_hamil_on_fly(kp, ivec)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec), nI(nel)
        integer :: det_ind, hdiag_ind, flag_ind, ideterm, DetHash
        integer(n_int) :: sgn(lenof_sign)
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign)
        real(dp) :: temp
        type(ll_node), pointer :: temp_node
        logical :: tDetFound, tDeterm
        character(len=*), parameter :: t_r = "calc_hamil_elems_on_fly"

        associate(h_matrix => kp%hamil_matrix)

            h_matrix(1:ivec, ivec) = 0.0_dp
            h_matrix(ivec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in krylov_vecs, for each Krylov vector.
                ind(jvec) = NIfDBO + lenof_sign*(jvec-1) + 1
            end do
            hdiag_ind = NIfDBO + lenof_sign_kp + 1
            if (tUseFlags) flag_ind = NIfDBO + lenof_sign_kp + 2

            ideterm = 0

            ! Loop over all determinants in SpawnedPartsKP.
            do idet = 1, max_spawned_ind
                sgn = SpawnedPartsKP(NOffSgn:NOffSgn+1, idet)
                sign1 = transfer(sgn, sign1)
                call decode_bit_det(nI, SpawnedPartsKP(:,idet))
                DetHash = FindWalkerHash(nI, nhashes_kp)
                ! Point to the first node with this hash value in krylov_vecs.
                temp_node => krylov_vecs_ht(DetHash)
                if (temp_node%ind == 0) then
                    ! If there are no determinants at all with this hash value in krylov_vecs.
                    cycle
                else
                    tDetFound = .false.
                    do while (associated(temp_node))
                        if (DetBitEQ(SpawnedPartsKP(:,idet), krylov_vecs(:,temp_node%ind), NIfDBO)) then
                            ! If this determinant has been found in krylov_vecs.
                            det_ind = temp_node%ind
                            tDetFound = .true.
                            exit
                        end if
                        ! Move on to the next determinant with this hash value.
                        temp_node => temp_node%next
                    end do
                    if (tDetFound) then
                        ! Add in the contribution to the projected Hamiltonian, for each Krylov vector.
                        do jvec = 1, ivec
                            sgn = krylov_vecs(ind(jvec):ind(jvec)+1, det_ind)
                            if (IsUnoccDet(sgn)) cycle
                            sign2 = transfer(sgn, sign1)
                            h_matrix(jvec,ivec) = h_matrix(jvec,ivec) - &
                                (sign1(1)*sign2(2) + sign1(2)*sign2(1))/2.0_dp
                        end do
                    end if
                end if
            end do

            ! Loop over all determinants in krylov_vecs.
            do idet = 1, TotWalkersKP

                sign1 = 0.0_dp
                tDeterm = .false.
                if (tUseFlags) then
                    tDeterm = btest(krylov_vecs(flag_ind, idet), flag_deterministic + flag_bit_offset)
                end if

                if (tDeterm) then
                    ideterm  = ideterm + 1
                    sign1 = - partial_determ_vector(:,ideterm) + &
                             (DiagSft+Hii) * tau * full_determ_vector(:, ideterm + determ_proc_indices(iProcIndex))
                else
                    sgn = krylov_vecs(ind(ivec):ind(ivec)+1, idet)
                    sign1 = transfer(sgn, sign1)
                    sign1 = tau * sign1 * (transfer(krylov_vecs(hdiag_ind, idet), temp) + Hii)
                end if
                if (IsUnoccDet(sign1)) cycle

                ! Loop over all Krylov vectors currently stored.
                do jvec = 1, ivec
                    sgn = krylov_vecs(ind(jvec):ind(jvec)+1, idet)
                    if (IsUnoccDet(sgn)) cycle
                    sign2 = transfer(sgn, sign1)
                    h_matrix(jvec,ivec) = h_matrix(jvec,ivec) + &
                        (sign1(1)*sign2(2) + sign1(2)*sign2(1))/2.0_dp
                end do
            end do

            do jvec = 1, ivec
                h_matrix(jvec,ivec) = h_matrix(jvec,ivec)/tau
                h_matrix(ivec,jvec) = h_matrix(jvec,ivec)
            end do

            if (tSemiStochastic) then
                if (ideterm /= determ_proc_sizes(iProcIndex)) then
                    write(6,*) "determ_proc_sizes(iProcIndex):", determ_proc_sizes(iProcIndex)
                    write(6,*) "ideterm:", ideterm
                    call neci_flush(6)
                    call stop_all(t_r, "An incorrect number of core determinants have been counted.")
                end if
            end if

        end associate

    end subroutine calc_hamil_on_fly

    subroutine calc_hamil_exact(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: i, j, idet, jdet, ic, hdiag_ind
        integer(n_int) :: ilut_1(0:NIfTot), ilut_2(0:NIfTot)
        integer(n_int) :: int_sign(lenof_sign_kp)
        integer :: nI(nel), nJ(nel)
        real(dp) :: real_sign_1(lenof_sign_kp), real_sign_2(lenof_sign_kp)
        real(dp) :: h_elem
        logical :: any_occ, occ_1, occ_2
        integer(4), allocatable :: occ_flags(:)

        hdiag_ind = NIfDBO + lenof_sign_kp + 1

        associate(h_matrix => kp%hamil_matrix)

            h_matrix = 0.0_dp

            allocate(occ_flags(TotWalkersKP))
            occ_flags = 0

            ilut_1 = 0_n_int
            ilut_2 = 0_n_int

            ! Check to see if there are any replica 1 or 2 walkers on this determinant.
            do idet = 1, TotWalkersKP
                int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, idet)

                any_occ = .false.
                do i = 1, kp%nvecs
                    any_occ = any_occ .or. (int_sign(2*i-1) /= 0)
                end do
                if (any_occ) occ_flags = ibset(occ_flags(idet), 0)

                any_occ = .false.
                do i = 1, kp%nvecs
                    any_occ = any_occ .or. (int_sign(2*i) /= 0)
                end do
                if (any_occ) occ_flags = ibset(occ_flags(idet), 1)
            end do

            ! Loop over all determinants in krylov_vecs.
            do idet = 1, TotWalkersKP
                ilut_1 = krylov_vecs(0:NIfDBO, idet)
                call decode_bit_det(nI, ilut_1)
                int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, idet)
                real_sign_1 = transfer(int_sign, real_sign_1)
                occ_1 = btest(occ_flags(idet),0)
                occ_2 = btest(occ_flags(idet),1)

                do jdet = idet, TotWalkersKP
                    if (.not. ((occ_1 .and. btest(occ_flags(jdet),1)) .or. &
                        (occ_2 .and. btest(occ_flags(jdet),0)))) cycle

                    ilut_2 = krylov_vecs(0:NIfDBO, jdet)
                    ic = FindBitExcitLevel(ilut_1, ilut_2)
                    if (ic > 2) cycle

                    call decode_bit_det(nJ, ilut_2)
                    int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, jdet)
                    real_sign_2 = transfer(int_sign, real_sign_1)
                    
                    if (idet == jdet) then
                        h_elem = transfer(krylov_vecs(hdiag_ind, idet), h_elem) + Hii
                    else
                        if (tHPHF) then
                            h_elem = hphf_off_diag_helement(nI, nJ, ilut_1, ilut_2)
                        else
                            h_elem = get_helement(nI, nJ, ic, ilut_1, ilut_2)
                        end if
                    end if

                    ! Finally, add in the contribution to all of the Hamiltonian elements.
                    do i = 1, kp%nvecs
                        do j = i, kp%nvecs
                            if (idet == jdet) then
                                h_matrix(i,j) = h_matrix(i,j) + &
                                    h_elem*(real_sign_1(2*i-1)*real_sign_2(2*j) + &
                                    real_sign_1(2*i)*real_sign_2(2*j-1))/2
                            else
                                h_matrix(i,j) = h_matrix(i,j) + &
                                    h_elem*(real_sign_1(2*i-1)*real_sign_2(2*j) + &
                                    real_sign_1(2*i)*real_sign_2(2*j-1) + &
                                    real_sign_1(2*j-1)*real_sign_2(2*i) + &
                                    real_sign_1(2*j)*real_sign_2(2*i-1))/2
                            end if
                        end do
                    end do

                end do
            end do

            do i = 1, kp%nvecs
                do j = 1, i-1
                    h_matrix(i,j) = h_matrix(j,i)
                end do
            end do

            deallocate(occ_flags)

        end associate

    end subroutine calc_hamil_exact

    subroutine communicate_kp_matrices(kp)

        ! Add all the overlap and projected Hamiltonian matrices together, with the result being
        ! held only on the root node.

        type(kp_fciqmc_data), intent(inout) :: kp
        real(dp) :: inp_matrices(2*kp%nvecs, kp%nvecs)
        real(dp) :: out_matrices(2*kp%nvecs, kp%nvecs)

        inp_matrices(1:kp%nvecs, 1:kp%nvecs) = kp%overlap_matrix
        inp_matrices(kp%nvecs+1:2*kp%nvecs, 1:kp%nvecs) = kp%hamil_matrix

        call MPISum(inp_matrices, out_matrices)

        if (iProcIndex == root) then
            kp%overlap_matrix = out_matrices(1:kp%nvecs, 1:kp%nvecs)
            kp%hamil_matrix = out_matrices(kp%nvecs+1:2*kp%nvecs, 1:kp%nvecs)
        end if

    end subroutine communicate_kp_matrices

    subroutine output_kp_matrices(kp, iconfig, irepeat)

        type(kp_fciqmc_data), intent(in) :: kp
        integer, intent(in) :: iconfig, irepeat
        real(dp) :: average_h_matrix(kp%nvecs, kp%nvecs)
        real(dp) :: average_s_matrix(kp%nvecs, kp%nvecs)

        if (iProcIndex == root) then
            call output_matrix(kp, iconfig, irepeat, 'hamil  ', kp%hamil_matrix)
            call output_matrix(kp, iconfig, irepeat, 'overlap', kp%overlap_matrix)
        end if

    end subroutine output_kp_matrices

    subroutine output_matrix(kp, iconfig, irepeat, stem, matrix)

        type(kp_fciqmc_data), intent(in) :: kp
        integer, intent(in) :: iconfig, irepeat
        character(7), intent(in) :: stem
        character(2) :: ifmt, jfmt
        real(dp), intent(in) :: matrix(kp%nvecs, kp%nvecs)
        character(25) :: ind1, ind2, filename
        integer :: i, j, ilen, jlen, new_unit

        ! Create the filename.
        write(ind2,'(i15)') irepeat
        if (kp%tGround) then
            filename = trim(trim(stem)//'.'//trim(adjustl(ind2)))
        else if (kp%tFiniteTemp) then
            write(ind1,'(i15)') iconfig
            filename = trim(trim(stem)//'.'//trim(adjustl(ind1))//'.'//trim(adjustl(ind2)))
        end if

        new_unit = get_free_unit()
        open(new_unit, file=trim(filename), status='replace')

        ! Write all the components of the matrix, above and including the diagonal, one
        ! after another on separate lines. Each element is preceeded by the indices involved.
        do i = 1, kp%nvecs
            ilen = ceiling(log10(real(abs(i)+1)))
            ! ifmt will hold the correct integer length so that there will be no spaces printed out.
            ! Note that this assumes that ilen < 10, which is very reasonable!
            write(ifmt,'(a1,i1)') "i", ilen
            do j = i, kp%nvecs
                jlen = ceiling(log10(real(abs(j)+1)))
                write(jfmt,'(a1,i1)') "i", jlen

                ! Finally write the line.
                write(new_unit,'(a1,'//ifmt//',a1,'//jfmt//',a1,1x,es19.12)') &
                    "(",i,",",j,")", matrix(i,j)
            end do
        end do
        close(new_unit)

    end subroutine output_matrix

    subroutine print_populations_kp(kp)
    
        ! A useful test routine which will output the total walker population on both
        ! replicas, for each Krylov vector.

        type(kp_fciqmc_data), intent(in) :: kp
        integer :: ihash
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign(lenof_sign_kp), total_pop(lenof_sign_kp)
        type(ll_node), pointer :: temp_node

        int_sign = 0_n_int
        total_pop = 0.0_dp
        real_sign = 0.0_dp
        
        do ihash = 1, nhashes_kp
            temp_node => krylov_vecs_ht(ihash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, temp_node%ind)
                    real_sign = transfer(int_sign, real_sign)
                    total_pop = total_pop + abs(real_sign)
                    temp_node => temp_node%next
                end do
            end if
        end do

        nullify(temp_node)

        write(6,*) "krylov_vecs populations:", total_pop

    end subroutine print_populations_kp

    subroutine print_amplitudes_kp(irepeat)

        ! A (*very* slow and memory intensive) test routine to print the current amplitudes (as stored
        ! in CurrentDets) of *all* determinants to a file. The amplitude of each replica will be printed
        ! one after the other. Since this is intended to be used with kp-fciqmc, irepeat is the number of
        ! the current repeat, but it will simply be used in naming the output file.

        ! Note that this routine will only work when using the tHashWalkerList option.

        integer, intent(in) :: irepeat
        integer, allocatable :: nI_list(:,:)
        integer :: temp(1,1), hf_ind, ndets
        integer :: i, j, ilen, counter, new_unit, DetHash
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int) :: sgn(lenof_sign)
        real(dp) :: real_sign(lenof_sign)
        type(ll_node), pointer :: temp_node
        type(BasisFn) :: iSym
        character(2) :: ifmt
        character(15) :: ind, filename

        ! Determine the total number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        allocate(nI_list(nel, ndets))

        ! Generate the determinants and move them to nI_list.
        ! Important: the above routine does not take symmetry into account. It returns all possible combinations.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, LMS, tParity, SymRestrict, ndets, hf_ind)

        write(ind,'(i15)') irepeat
        filename = trim('amps.'//adjustl(ind))

        new_unit = get_free_unit()
        open(new_unit, file=trim(filename), status='replace')

        counter = 0

        do i = 1, ndets
            call getsym(nI_list(:,i), nel, G1, nBasisMax, iSym)
            ! Only carry on if the symmetry of this determinant is correct.
            if (iSym%Sym%S /= HFSym%Sym%S .or. iSym%Ms /= HFSym%Ms .or. iSym%Ml /= HFSym%Ml) cycle
            call EncodeBitDet(nI_list(:,i), ilut)
            if (.not. IsAllowedHPHF(ilut(0:NIfDBO))) cycle
            counter = counter + 1
            real_sign = 0.0_dp
            DetHash = FindWalkerHash(nI_list(:,i), nWalkerHashes)
            temp_node => HashIndex(DetHash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if (DetBitEQ(ilut, CurrentDets(:,temp_node%ind), NIfDBO)) then
                        sgn = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,temp_node%ind)
                        real_sign = transfer(sgn, real_sign)
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            ilen = ceiling(log10(real(abs(counter)+1)))
            ! ifmt will hold the correct integer length so that there will be no spaces printed out.
            ! Note that this assumes that ilen < 10, which is very reasonable!
            write(ifmt,'(a1,i1)') "i", ilen
            do j = 1, lenof_sign
                ! This assumes that lenof_sign < 10. Probably will always be 2.
                write(new_unit,'(a1,'//ifmt//',a1,i1,a1,1x,es19.12)') "(", counter,",", j, ")", real_sign(j)
            end do
        end do

        close(new_unit)

        deallocate(nI_list)

    end subroutine print_amplitudes_kp

end module kp_fciqmc_procs
