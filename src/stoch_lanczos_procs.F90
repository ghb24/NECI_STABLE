#include "macros.h"

module stoch_lanczos_procs

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use bit_rep_data
    use bit_reps, only: decode_bit_det
    use CalcData, only: tTruncInitiator, tStartSinglePart, InitialPart, InitWalkers
    use CalcData, only: tSemiStochastic, tReadPops, tUseRealCoeffs
    use constants
    use DetBitOps, only: DetBitEq
    use dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData, only: ilutHF, HFDet, CurrentDets, SpawnedParts, SpawnedParts2, TotWalkers
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, HashIndex, nWalkerHashes
    use FciMCData, only: fcimc_iter_data, ll_node, MaxWalkersPart, tStartCoreGroundState
    use FciMCData, only: tPopsAlreadyRead, tHashWalkerList, CurrentH, determ_proc_sizes
    use FciMCData, only: core_ham_diag
    use FciMCParMod, only: create_particle, InitFCIMC_HF, SetupParameters, InitFCIMCCalcPar
    use FciMCParMod, only: init_fcimc_fn_pointers, WriteFciMCStats, WriteFciMCStatsHeader
    use FciMCParMod, only: rezero_iter_stats_each_iter
    use hash, only: FindWalkerHash, init_hash_table, reset_hash_table, fill_in_hash_table
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias
    use Parallel_neci, only: MPIBarrier, iProcIndex
    use procedure_pointers
    use semi_stoch_procs, only: copy_core_dets_this_proc_to_spawnedparts
    use semi_stoch_procs, only: add_core_states_currentdet_hash, start_walkers_from_core_ground
    use SystemData, only: nel

    implicit none

    type stoch_lanczos_data
        ! If true then a ground-state calculation is being performed.
        logical :: tGround
        ! If true then a finite-temperature calculation is being performed.
        logical :: tFiniteTemp
        ! The number of different initial walker configurations to start
        ! Lanczos calculations from.
        integer :: nconfigs
        ! The number of Lanczos calculations to perform for each initial walker
        ! configuration.
        integer :: nrepeats
        ! The number of different Lanczos vectors to sample (the number of
        ! vectors which form the Krylov subspace at the end of a calculation).
        integer :: nvecs
        ! The number of iterations to perform *between each Lanczos vector being
        ! sampled*.
        integer :: niters
        real(dp), allocatable :: overlap_matrix(:,:)
        real(dp), allocatable :: hamil_matrix(:,:)
    end type

    type(stoch_lanczos_data) :: lanczos

    integer :: nhashes_lanczos
    integer :: TotWalkersLanczos
    integer(n_int), allocatable :: lanczos_vecs(:,:)
    type(ll_node), pointer :: lanczos_hash_table(:) 

    integer :: TotWalkers_Lanc
    integer(n_int), allocatable :: init_lanczos_config(:,:)

contains

    subroutine stoch_lanczos_read_inp()

        use input_neci

        logical :: eof
        character(len=100) :: w

        lanczos%nconfigs = 1
        lanczos%nrepeats = 1
        lanczos%tGround = .false.
        lanczos%tFiniteTemp = .false.

        read_inp: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)
            case("END-LANCZOS")
                exit read_inp
            case("GROUND-STATE")
                lanczos%tGround = .true.
            case("FINITE-TEMPERATURE")
                lanczos%tFiniteTemp = .true.
            case("NUM-INIT-CONFIGS")
                call geti(lanczos%nconfigs)
            case("NUM-REPEATS-PER-CONFIG")
                call geti(lanczos%nrepeats)
            case("NUM-LANCZOS-VECS")
                call geti(lanczos%nvecs)
            case("NUM-ITERS-PER-VEC")
                call geti(lanczos%niters)
            case default
                call report("Keyword "//trim(w)//" not recognized in stoch-lanczos block", .true.)
            end select
        end do read_inp

    end subroutine stoch_lanczos_read_inp

    subroutine init_stoch_lanczos(lanczos)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer :: ierr
        character (len=*), parameter :: t_r = "init_stoch_lanczos"

        associate(nvecs => lanczos%nvecs, s_matrix => lanczos%overlap_matrix, hamil => lanczos%hamil_matrix)

        if (.not. tHashWalkerList) call stop_all('t_r','Stochastic Lanczos can only be run using &
            &the linscalefcimcalgo option (the linear scaling algorithm).')

        if (.not. tUseRealCoeffs) call stop_all('t_r','Stochastic Lanczos can only be run using &
            &real coefficients).')

        tPopsAlreadyRead = .false.
        call SetupParameters()
        call InitFCIMCCalcPar()
        call init_fcimc_fn_pointers() 

        call WriteFciMCStatsHeader()
        call WriteFCIMCStats()

        if (n_int == 4) call stop_all('t_r', 'Use of RealCoefficients does not work with 32 bit &
             &integers due to the use of the transfer operation from dp reals to 64 bit integers.')

        ! Assuming Yamanouchi symbols (and CSFs) not used.
        NIfLan = NIfD + lenof_sign*lanczos%nvecs + NIfFlag

        nhashes_lanczos = nWalkerHashes
        TotWalkersLanczos = 0
        allocate(lanczos_vecs(0:NIfLan, MaxWalkersPart), stat=ierr)
        lanczos_vecs = 0
        allocate(lanczos_hash_table(nhashes_lanczos), stat=ierr)
        call init_hash_table(lanczos_hash_table)

        allocate(lanczos%overlap_matrix(lanczos%nvecs, lanczos%nvecs), stat=ierr)
        allocate(lanczos%hamil_matrix(lanczos%nvecs, lanczos%nvecs), stat=ierr)

        ! If performing a finite-temperature calculation with more than one run for each initial
        ! configuration, we store this walker configuration so that we can restart from it later.
        if (lanczos%tFiniteTemp .and. lanczos%nrepeats > 1) then
            allocate(init_lanczos_config(0:NIfTot, MaxWalkersPart), stat=ierr)
        end if

        call MPIBarrier(ierr)

        end associate

    end subroutine init_stoch_lanczos

    subroutine init_stoch_lanczos_iter(iter_data, determ_index)

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

    end subroutine init_stoch_lanczos_iter

    subroutine create_initial_config(lanczos, irepeat)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer, intent(in) :: irepeat
        integer :: DetHash, i, nwalkers, nwalkers_target

        call reset_hash_table(HashIndex)

        if (lanczos%tGround) then
            ! Put a walker on the Hartree-Fock again, with the requested amplitude.
            call InitFCIMC_HF()
            if (tSemiStochastic) then
                ! core_space stores all core determinants from all processors. Move those on this
                ! processor to SpawnedParts, which add_core_states_currentdet_hash uses.
                call copy_core_dets_this_proc_to_spawnedparts()
                call add_core_states_currentdet_hash()
                if (tStartCoreGroundState .and. (.not. tReadPops)) &
                    call start_walkers_from_core_ground(tPrintInfo = .false.)
                ! Reset the diagonal Hamiltonian elements.
                CurrentH(1, 1:determ_proc_sizes(iProcIndex)) = core_ham_diag
            end if
        else if (lanczos%tFiniteTemp) then
            if (irepeat == 1) then
                ! Convert the initial number of walkers to an integer.
                if (tStartSinglePart) then
                    nwalkers_target = int(InitialPart)
                else
                    nwalkers_target = int(InitWalkers)
                end if
                ! Finally, call the routine to create the walker distribution.
                call generate_init_config_basic(nwalkers_target, nwalkers)
                TotWalkers = int(nwalkers, int64)
                TotWalkers_Lanc = TotWalkers
                ! If starting from this configuration more than once, store it.
                if (lanczos%nrepeats > 1) init_lanczos_config(:, 1:TotWalkers) = CurrentDets(:, 1:TotWalkers)
            else if (irepeat > 1) then
                TotWalkers = TotWalkers_Lanc
                CurrentDets(:, 1:TotWalkers) = init_lanczos_config(:, 1:TotWalkers)
                call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, TotWalkers)
            end if
        end if

    end subroutine create_initial_config

    subroutine generate_init_config_basic(nwalkers, ndets)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        integer, intent(in) :: nwalkers
        integer, intent(out) :: ndets
        integer :: i, ireplica, excit, nattempts, DetHash
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_sign(lenof_sign)
        logical :: tInitiatorTemp
        type(fcimc_iter_data) :: temp_data
        integer(n_int), pointer :: PointTemp(:,:)

        ! Turn off the initiator method for the annihilation steps to be used here.
        tInitiatorTemp = tTruncInitiator
        tTruncInitiator = .false.

        ! Set the spawning slots to their starting positions.
        ValidSpawnedList = InitialSpawnedSlots

        ilut = 0

        do ireplica = 1, inum_runs
            walker_sign = 0.0_dp
            do i = 1, nwalkers
                ! Generate the determinant (output to ilut).
                call CreateRandomExcitLevDetUnbias(nel, HFDet, ilutHF, ilut, excit, nattempts)
                call decode_bit_det(nI, ilut)

                ! Choose whether the walker to be added has an amplitude of plus or minus one, with
                ! 0.5 chance of each.
                walker_sign(ireplica) = 1.0_dp
                r = genrand_real2_dSFMT()
                if (r < 0.5) walker_sign(ireplica) = -1.0_dp*walker_sign(ireplica)

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
        call CompressSpawnedList(ndets, temp_data) 

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        CurrentDets = SpawnedParts

        ! Add the entries into the hash table.
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets)

        ValidSpawnedList = InitialSpawnedSlots
        SpawnedParts = 0

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

    end subroutine generate_init_config_basic

    subroutine store_lanczos_vec(ivec)

        integer, intent(in) ::  ivec
        integer :: idet, sign_ind, DetHash, det_ind
        integer :: nI(nel)
        logical :: tDetFound
        type(ll_node), pointer :: temp_node

        ! The index of the first element referring to the sign, for this ivec.
        sign_ind = NIfD + lenof_sign*(ivec-1) + 1

        ! Loop over all occupied determinants for this new Lanczos vector.
        do idet = 1, TotWalkers
            tDetFound = .false.

            call decode_bit_det(nI, CurrentDets(:,idet))
            DetHash = FindWalkerHash(nI, nhashes_lanczos)
            temp_node => lanczos_hash_table(DetHash)

            ! If the first element in the list for this hash value has been used.
            if (.not. temp_node%ind == 0) then
                ! Loop over all determinants with this hash value which are already in the list.
                do while (associated(temp_node%next))
                    if (DetBitEQ(CurrentDets(:,idet), lanczos_vecs(:,TotWalkersLanczos), NIfDBO)) then
                        ! This determinant is already in the list.
                        det_ind = temp_node%ind
                        ! Just add the sign for the new Lanczos vector. The determinant and flag are
                        ! there already.
                        lanczos_vecs(sign_ind:sign_ind+NIfSgn,det_ind) = CurrentDets(NOffSgn:NOffSgn+NIfSgn,idet)
                        tDetFound = .true.
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    temp_node => temp_node%next
                end do

                if (.not. tDetFound) then
                    ! We need to add a new determinant in the next position in the list.
                    ! So create that next position!
                    allocate(temp_node%next)
                    temp_node => temp_node%next
                    nullify(temp_node%next)
                end if
            end if

            if (.not. tDetFound) then
                ! A new determiant needs to be added.
                TotWalkersLanczos = TotWalkersLanczos + 1
                det_ind = TotWalkersLanczos
                temp_node%ind = det_ind

                ! Copy determinant across.
                lanczos_vecs(0:NIfD,det_ind) = CurrentDets(0:NIfD,idet)
                ! Copy sign across.
                lanczos_vecs(sign_ind:sign_ind+NIfSgn,det_ind) = CurrentDets(NOffSgn:NOffSgn+NIfSgn,idet)
                ! Copy Flags across.
                lanczos_vecs(NOffFlag,det_ind) = CurrentDets(NOffFlag,idet)
            end if

            nullify(temp_node)
        end do

    end subroutine store_lanczos_vec

    subroutine calc_overlap_matrix_elems(lanczos, ivec)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec), sgn(lenof_sign)
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign)

        associate(s_matrix => lanczos%overlap_matrix)

            s_matrix(1:ivec, jvec) = 0.0_dp
            s_matrix(jvec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in lanczos_vecs, for each Lanczos vector.
                ind(jvec) = NIfD + lenof_sign*(jvec-1) + 1
            end do

            do idet = 1, TotWalkersLanczos
                sgn = lanczos_vecs(ind(ivec):ind(ivec)+1, idet)
                sign1 = transfer(sgn, sign1)
                do jvec = 1, ivec
                    sgn = lanczos_vecs(ind(jvec):ind(jvec)+1, idet)
                    sign2 = transfer(sgn, sign1)
                    s_matrix(jvec,ivec) = s_matrix(jvec,ivec) + (sign1(1)*sign2(2) + sign1(2)*sign2(1))/2.0_dp
                end do
            end do

        end associate

    end subroutine calc_overlap_matrix_elems

end module stoch_lanczos_procs
