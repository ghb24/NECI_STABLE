#include "macros.h"

module stoch_lanczos_procs

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use bit_rep_data
    use bit_reps, only: decode_bit_det
    use CalcData, only: tTruncInitiator, tStartSinglePart, InitialPart, InitWalkers
    use CalcData, only: tSemiStochastic, tReadPops, tUseRealCoeffs, tau, DiagSft
    use constants
    use DetBitOps, only: DetBitEq, EncodeBitDet, IsAllowedHPHF
    use dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData, only: ilutHF, HFDet, CurrentDets, SpawnedParts, SpawnedParts2, TotWalkers
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, HashIndex, nWalkerHashes
    use FciMCData, only: fcimc_iter_data, ll_node, MaxWalkersPart, tStartCoreGroundState
    use FciMCData, only: tPopsAlreadyRead, tHashWalkerList, CurrentH, determ_proc_sizes
    use FciMCData, only: core_ham_diag, InputDiagSft, Hii, max_spawned_ind, SpawnedPartsLanc
    use FciMCData, only: partial_determ_vector, full_determ_vector, determ_proc_indices, HFSym
    use FciMCData, only: TotParts, TotPartsOld, AllTotParts, AllTotPartsOld
    use FciMCParMod, only: create_particle, InitFCIMC_HF, SetupParameters, InitFCIMCCalcPar
    use FciMCParMod, only: init_fcimc_fn_pointers, WriteFciMCStats, WriteFciMCStatsHeader
    use FciMCParMod, only: rezero_iter_stats_each_iter, tSinglePartPhase
    use gndts_mod, only: gndts
    use hash, only: FindWalkerHash, init_hash_table, reset_hash_table, fill_in_hash_table
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias
    use Parallel_neci, only: MPIBarrier, iProcIndex, MPISum, MPIReduce
    use ParallelHelper, only: root
    use procedure_pointers
    use semi_stoch_procs, only: copy_core_dets_this_proc_to_spawnedparts, fill_in_CurrentH
    use semi_stoch_procs, only: add_core_states_currentdet_hash, start_walkers_from_core_ground
    use sym_mod, only: getsym
    use SystemData, only: nel, nbasis, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict
    use SystemData, only: BasisFn
    use util_mod, only: get_free_unit

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
        real(dp), allocatable :: overlap_matrix_1(:,:), overlap_matrix_2(:,:)
        real(dp), allocatable :: hamil_matrix_1(:,:), hamil_matrix_2(:,:)
    end type

    type(stoch_lanczos_data) :: lanczos

    integer :: nhashes_lanczos
    integer :: TotWalkersLanczos
    integer(n_int), allocatable :: lanczos_vecs(:,:)
    type(ll_node), pointer :: lanczos_hash_table(:) 

    integer :: TotWalkersInit
    integer :: TotPartsInit(lenof_sign)
    integer :: AllTotPartsInit(lenof_sign)
    integer(n_int), allocatable :: init_lanczos_config(:,:)

contains

    subroutine stoch_lanczos_read_inp()

        use input_neci

        logical :: eof
        character(len=100) :: w

       ! Default values.
        lanczos%nconfigs = 1
        lanczos%nrepeats = 1
        lanczos%nvecs = 1
        lanczos%niters = 1
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
            case("NUM-REPEATS-PER-INIT-CONFIG")
                call geti(lanczos%nrepeats)
            case("NUM-LANCZOS-VECS")
                call geti(lanczos%nvecs)
            case("NUM-ITERS-BETWEEN-VECS")
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
        NIfLan = NIfD + lenof_sign*lanczos%nvecs + 1 + NIfFlag

        nhashes_lanczos = nWalkerHashes
        TotWalkersLanczos = 0
        allocate(lanczos_vecs(0:NIfLan, MaxWalkersPart), stat=ierr)
        lanczos_vecs = 0_n_int
        allocate(lanczos_hash_table(nhashes_lanczos), stat=ierr)
        call init_hash_table(lanczos_hash_table)

        allocate(lanczos%overlap_matrix_1(lanczos%nvecs, lanczos%nvecs), stat=ierr)
        allocate(lanczos%overlap_matrix_2(lanczos%nvecs, lanczos%nvecs), stat=ierr)
        allocate(lanczos%hamil_matrix_1(lanczos%nvecs, lanczos%nvecs), stat=ierr)
        allocate(lanczos%hamil_matrix_2(lanczos%nvecs, lanczos%nvecs), stat=ierr)

        ! If performing a finite-temperature calculation with more than one run for each initial
        ! configuration, we store this walker configuration so that we can restart from it later.
        if (lanczos%tFiniteTemp .and. lanczos%nrepeats > 1) then
            allocate(init_lanczos_config(0:NIfTot, MaxWalkersPart), stat=ierr)
        end if

        call MPIBarrier(ierr)

    end subroutine init_stoch_lanczos

    subroutine init_stoch_lanczos_repeat(lanczos, irepeat)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer :: irepeat

        call create_initial_config(lanczos, irepeat)

        call reset_hash_table(lanczos_hash_table)
        lanczos_vecs = 0_n_int

        lanczos%overlap_matrix_1 = 0.0_dp
        lanczos%overlap_matrix_2 = 0.0_dp
        lanczos%hamil_matrix_1 = 0.0_dp
        lanczos%hamil_matrix_2 = 0.0_dp
        TotWalkersLanczos = 0

        DiagSft = InputDiagSft
        ! Setting this variable to true stops the shift from varying instantly.
        tSinglePartPhase = .true.

    end subroutine init_stoch_lanczos_repeat

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
                SpawnedParts = 0_n_int
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
                TotWalkersInit = TotWalkers
                TotPartsInit = TotParts
                AllTotPartsInit = AllTotParts
                ! If starting from this configuration more than once, store it.
                if (lanczos%nrepeats > 1) init_lanczos_config(:,1:TotWalkers) = CurrentDets(:,1:TotWalkers)
            else if (irepeat > 1) then
                TotWalkers = TotWalkersInit
                TotParts = TotPartsInit
                TotPartsOld = TotPartsInit
                AllTotParts = AllTotPartsInit
                AllTotPartsOld = AllTotPartsInit

                CurrentDets(:,1:TotWalkers) = init_lanczos_config(:,1:TotWalkers)
                call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, TotWalkers)
            end if
            call fill_in_CurrentH()
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
        type(fcimc_iter_data) :: unused_data
        integer(n_int), pointer :: PointTemp(:,:)

        ! Turn off the initiator method for the annihilation steps to be used here.
        tInitiatorTemp = tTruncInitiator
        tTruncInitiator = .false.

        ! Set the spawning slots to their starting positions.
        ValidSpawnedList = InitialSpawnedSlots

        ilut = 0_n_int

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
        call CompressSpawnedList(ndets, unused_data) 

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        TotParts = 0.0_dp
        do i = 1, ndets
            CurrentDets(:,i) = SpawnedParts(:,i)
            walker_sign = transfer(CurrentDets(NOffSgn:NOffSgn+lenof_sign, i), walker_sign)
            TotParts = TotParts + abs(walker_sign)
        end do
        TotPartsOld = TotParts

        ! Add the entries into the hash table.
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets)

        call MPIReduce(TotParts, MPI_SUM, AllTotParts)
        AllTotPartsOld = AllTotParts

        ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
        ! These routines will do this.
        TotWalkers = int(ndets, int64)
        call copy_core_dets_this_proc_to_spawnedparts()
        call add_core_states_currentdet_hash()

        ValidSpawnedList = InitialSpawnedSlots
        SpawnedParts = 0_n_int

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

    end subroutine generate_init_config_basic

    subroutine store_lanczos_vec(ivec, nvecs)

        integer, intent(in) ::  ivec, nvecs
        integer :: idet, sign_ind, hdiag_ind, flag_ind, DetHash, det_ind
        integer :: nI(nel)
        integer(n_int) :: temp
        logical :: tDetFound
        type(ll_node), pointer :: temp_node, prev

        ! The index of the first element referring to the sign, for this ivec.
        sign_ind = NIfD + lenof_sign*(ivec-1) + 1
        hdiag_ind = NIfD + lenof_sign*nvecs + 1
        if (tUseFlags) flag_ind = NIfD + lenof_sign*nvecs + 2

        ! Loop over all occupied determinants for this new Lanczos vector.
        do idet = 1, TotWalkers
            tDetFound = .false.

            call decode_bit_det(nI, CurrentDets(:,idet))
            DetHash = FindWalkerHash(nI, nhashes_lanczos)
            temp_node => lanczos_hash_table(DetHash)

            ! If the first element in the list for this hash value has been used.
            if (.not. temp_node%ind == 0) then
                ! Loop over all determinants with this hash value which are already in the list.
                do while (associated(temp_node))
                    if (DetBitEQ(CurrentDets(:,idet), lanczos_vecs(:,temp_node%ind), NIfDBO)) then
                        ! This determinant is already in the list.
                        det_ind = temp_node%ind
                        ! Add the sign for the new Lanczos vector. The determinant and flag are
                        ! there already.
                        lanczos_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = &
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
                TotWalkersLanczos = TotWalkersLanczos + 1
                det_ind = TotWalkersLanczos
                temp_node%ind = det_ind

                ! Copy determinant data across.
                lanczos_vecs(0:NIfD,det_ind) = CurrentDets(0:NIfD,idet)
                lanczos_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,idet)
                lanczos_vecs(hdiag_ind,det_ind) = transfer(CurrentH(1,idet), temp)
                if (tUseFlags) lanczos_vecs(flag_ind,det_ind) = CurrentDets(NOffFlag,idet)
            end if

            nullify(temp_node)
            nullify(prev)
        end do

    end subroutine store_lanczos_vec

    subroutine calc_overlap_matrix_elems(lanczos, ivec)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec)
        integer(n_int) :: sgn(lenof_sign)
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign)

        associate(s_matrix1 => lanczos%overlap_matrix_1, s_matrix2 => lanczos%overlap_matrix_2)

            ! Just in case!
            s_matrix1(1:ivec, ivec) = 0.0_dp
            s_matrix1(ivec, 1:ivec) = 0.0_dp
            s_matrix2(1:ivec, ivec) = 0.0_dp
            s_matrix2(ivec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in lanczos_vecs, for each Lanczos vector.
                ind(jvec) = NIfD + lenof_sign*(jvec-1) + 1
            end do

            ! Loop over all determinants in lanczos_vecs.
            do idet = 1, TotWalkersLanczos
                sgn = lanczos_vecs(ind(ivec):ind(ivec)+1, idet)
                if (IsUnoccDet(sgn)) cycle
                sign1 = transfer(sgn, sign1)
                ! Loop over all lanczos vectors currently stored.
                do jvec = 1, ivec
                    sgn = lanczos_vecs(ind(jvec):ind(jvec)+1, idet)
                    if (IsUnoccDet(sgn)) cycle
                    sign2 = transfer(sgn, sign1)
                    s_matrix1(jvec,ivec) = s_matrix1(jvec,ivec) + sign1(2)*sign2(1)
                    s_matrix2(jvec,ivec) = s_matrix2(jvec,ivec) + sign1(1)*sign2(2)
                end do
            end do

            ! Fill in the lower-half of the overlap matrices.
            do jvec = 1, ivec
                s_matrix1(ivec,jvec) = s_matrix1(jvec,ivec)
                s_matrix2(ivec,jvec) = s_matrix2(jvec,ivec)
            end do

        end associate

    end subroutine calc_overlap_matrix_elems

    subroutine calc_hamil_elems(lanczos, ivec)

        ! Note: this doesn't work!

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec)
        integer :: det_ind, DetHash, nI(nel)
        integer(n_int) :: sgn(lenof_sign)
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign), full_shift(lenof_sign)
        type(ll_node), pointer :: temp_node
        logical :: tDetFound

        associate(h_matrix => lanczos%hamil_matrix_1, &
                  s_matrix1 => lanczos%overlap_matrix_1, s_matrix2 => lanczos%overlap_matrix_2)

            h_matrix(1:ivec, ivec) = 0.0_dp
            h_matrix(ivec, 1:ivec) = 0.0_dp
            ! The full shift, including the Hartree-Fock energy, *not* relative to it.
            full_shift = DiagSft + Hii

            do jvec = 1, ivec
                ! The first index of the sign in lanczos_vecs, for each Lanczos vector.
                ind(jvec) = NIfD + lenof_sign*(jvec-1) + 1
            end do

            ! Loop over all determinants in CurrentDets.
            do idet = 1, TotWalkers
                sgn = CurrentDets(NOffSgn:NOffSgn+1, idet)
                sign1 = transfer(sgn, sign1)
                call decode_bit_det(nI, CurrentDets(:,idet))
                DetHash = FindWalkerHash(nI, nhashes_lanczos)
                ! Point to the first node with this hash value in lanczos_vecs.
                temp_node => lanczos_hash_table(DetHash)
                if (temp_node%ind == 0) then
                    ! If there are no determinants at all with this hash value in lanczos_vecs.
                    cycle
                else
                    tDetFound = .false.
                    do while (associated(temp_node))
                        if (DetBitEQ(CurrentDets(:,idet), lanczos_vecs(:,temp_node%ind), NIfDBO)) then
                            ! If this CurrentDets determinant has been found in lanczos_vecs.
                            det_ind = temp_node%ind
                            tDetFound = .true.
                            exit
                        end if
                        ! Move on to the next determinant with this hash value.
                        temp_node => temp_node%next
                    end do
                    if (tDetFound) then
                        ! Add in the contribution to the projected Hamiltonian, for each lanczos vector.
                        do jvec = 1, ivec
                            sgn = lanczos_vecs(ind(jvec):ind(jvec)+1, det_ind)
                            if (IsUnoccDet(sgn)) cycle
                            sign2 = transfer(sgn, sign1)
                            h_matrix(jvec,ivec) = h_matrix(jvec,ivec) + (sign1(1)*sign2(2) + sign1(2)*sign2(1))/2.0_dp
                        end do
                    end if
                end if
            end do

            ! The above actually holds the subspace version of (1-\tau(H-S)) acting on the most recent
            ! Lanczos vector. But we want just H acting on it. We therefore have to add some extra contributions
            ! to get this. Note, we have to be careful about the two different overlap matrices to be included.
            do jvec = 1, ivec
                h_matrix(jvec,ivec) = -h_matrix(jvec,ivec) + 0.5_dp*(1+tau*full_shift(1))*s_matrix1(jvec,ivec) + &
                                                            0.5_dp*(1+tau*full_shift(2))*s_matrix2(jvec,ivec)
                h_matrix(jvec,ivec) = h_matrix(jvec,ivec)/tau
                h_matrix(ivec,jvec) = h_matrix(jvec,ivec)
            end do

            lanczos%hamil_matrix_2 = h_matrix

        end associate

    end subroutine calc_hamil_elems

    subroutine calc_hamil_elems_direct(lanczos, ivec)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec), nI(nel)
        integer :: det_ind, hdiag_ind, flag_ind, ideterm, DetHash
        integer(n_int) :: sgn(lenof_sign)
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign)
        real(dp) :: temp
        type(ll_node), pointer :: temp_node
        logical :: tDetFound, tDeterm
        character(len=*), parameter :: t_r = "calc_hamil_elems_direct"

        associate(h_matrix_1 => lanczos%hamil_matrix_1, h_matrix_2 => lanczos%hamil_matrix_2)

            h_matrix_1(1:ivec, ivec) = 0.0_dp
            h_matrix_1(ivec, 1:ivec) = 0.0_dp
            h_matrix_2(1:ivec, ivec) = 0.0_dp
            h_matrix_2(ivec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in lanczos_vecs, for each Lanczos vector.
                ind(jvec) = NIfD + lenof_sign*(jvec-1) + 1
            end do
            hdiag_ind = NIfD + lenof_sign*lanczos%nvecs + 1
            if (tUseFlags) flag_ind = NIfD + lenof_sign*lanczos%nvecs + 2

            ideterm = 0

            ! Loop over all determinants in SpawnedPartsLanc.
            do idet = 1, max_spawned_ind
                sgn = SpawnedPartsLanc(NOffSgn:NOffSgn+1, idet)
                sign1 = transfer(sgn, sign1)
                call decode_bit_det(nI, SpawnedPartsLanc(:,idet))
                DetHash = FindWalkerHash(nI, nhashes_lanczos)
                ! Point to the first node with this hash value in lanczos_vecs.
                temp_node => lanczos_hash_table(DetHash)
                if (temp_node%ind == 0) then
                    ! If there are no determinants at all with this hash value in lanczos_vecs.
                    cycle
                else
                    tDetFound = .false.
                    do while (associated(temp_node))
                        if (DetBitEQ(SpawnedPartsLanc(:,idet), lanczos_vecs(:,temp_node%ind), NIfDBO)) then
                            ! If this determinant has been found in lanczos_vecs.
                            det_ind = temp_node%ind
                            tDetFound = .true.
                            exit
                        end if
                        ! Move on to the next determinant with this hash value.
                        temp_node => temp_node%next
                    end do
                    if (tDetFound) then
                        ! Add in the contribution to the projected Hamiltonian, for each lanczos vector.
                        do jvec = 1, ivec
                            sgn = lanczos_vecs(ind(jvec):ind(jvec)+1, det_ind)
                            if (IsUnoccDet(sgn)) cycle
                            sign2 = transfer(sgn, sign1)
                            h_matrix_1(jvec,ivec) = h_matrix_1(jvec,ivec) - sign1(1)*sign2(2)
                            h_matrix_2(jvec,ivec) = h_matrix_2(jvec,ivec) - sign1(2)*sign2(1)
                        end do
                    end if
                end if
            end do

            ! Loop over all determinants in lanczos_vecs.
            do idet = 1, TotWalkersLanczos

                sign1 = 0.0_dp
                tDeterm = .false.
                if (tUseFlags) then
                    tDeterm = btest(lanczos_vecs(flag_ind, idet), flag_deterministic + flag_bit_offset)
                end if

                if (tDeterm) then
                    ideterm  = ideterm + 1
                    sign1 = - partial_determ_vector(:,ideterm) + &
                             (DiagSft+Hii) * tau * full_determ_vector(:, ideterm + determ_proc_indices(iProcIndex))
                else
                    sgn = lanczos_vecs(ind(ivec):ind(ivec)+1, idet)
                    sign1 = transfer(sgn, sign1)
                    sign1 = tau * sign1 * (transfer(lanczos_vecs(hdiag_ind, idet), temp) + Hii)
                end if
                if (IsUnoccDet(sign1)) cycle

                ! Loop over all lanczos vectors currently stored.
                do jvec = 1, ivec
                    sgn = lanczos_vecs(ind(jvec):ind(jvec)+1, idet)
                    if (IsUnoccDet(sgn)) cycle
                    sign2 = transfer(sgn, sign1)
                    h_matrix_1(jvec,ivec) = h_matrix_1(jvec,ivec) + sign1(1)*sign2(2)
                    h_matrix_2(jvec,ivec) = h_matrix_2(jvec,ivec) + sign1(2)*sign2(1)
                end do
            end do

            do jvec = 1, ivec
                h_matrix_1(jvec,ivec) = h_matrix_1(jvec,ivec)/tau
                h_matrix_2(jvec,ivec) = h_matrix_2(jvec,ivec)/tau
                h_matrix_1(ivec,jvec) = h_matrix_1(jvec,ivec)
                h_matrix_2(ivec,jvec) = h_matrix_2(jvec,ivec)
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

    end subroutine calc_hamil_elems_direct

    subroutine communicate_lanczos_matrices(lanczos)

        ! Add all the overlap and projected Hamiltonian matrices together, with the result being
        ! held only on the root node.

        type(stoch_lanczos_data), intent(inout) :: lanczos
        real(dp) :: inp_matrices(4*lanczos%nvecs, lanczos%nvecs)
        real(dp) :: out_matrices(4*lanczos%nvecs, lanczos%nvecs)

        inp_matrices(1:lanczos%nvecs, 1:lanczos%nvecs) = lanczos%overlap_matrix_1
        inp_matrices(lanczos%nvecs+1:2*lanczos%nvecs, 1:lanczos%nvecs) = lanczos%overlap_matrix_2
        inp_matrices(2*lanczos%nvecs+1:3*lanczos%nvecs, 1:lanczos%nvecs) = lanczos%hamil_matrix_1
        inp_matrices(3*lanczos%nvecs+1:4*lanczos%nvecs, 1:lanczos%nvecs) = lanczos%hamil_matrix_2

        call MPISum(inp_matrices, out_matrices)

        if (iProcIndex == root) then
            lanczos%overlap_matrix_1 = out_matrices(1:lanczos%nvecs, 1:lanczos%nvecs)
            lanczos%overlap_matrix_2 = out_matrices(lanczos%nvecs+1:2*lanczos%nvecs, 1:lanczos%nvecs)
            lanczos%hamil_matrix_1 = out_matrices(2*lanczos%nvecs+1:3*lanczos%nvecs, 1:lanczos%nvecs)
            lanczos%hamil_matrix_2 = out_matrices(3*lanczos%nvecs+1:4*lanczos%nvecs, 1:lanczos%nvecs)
        end if

    end subroutine communicate_lanczos_matrices

    subroutine output_lanczos_matrices(lanczos, iconfig, irepeat)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer, intent(in) :: iconfig, irepeat
        real(dp) :: average_h_matrix(lanczos%nvecs, lanczos%nvecs)
        real(dp) :: average_s_matrix(lanczos%nvecs, lanczos%nvecs)

        if (iProcIndex == root) then
            average_h_matrix = (lanczos%hamil_matrix_1 + lanczos%hamil_matrix_2)/2.0_dp
            call output_matrix(lanczos, iconfig, irepeat, 'hamil  ', average_h_matrix)
            ! We have two overlap matrices as we have two replicas. So average them for
            ! better statistics.
            average_s_matrix = (lanczos%overlap_matrix_1 + lanczos%overlap_matrix_2)/2.0_dp
            call output_matrix(lanczos, iconfig, irepeat, 'overlap', average_s_matrix)
        end if

    end subroutine output_lanczos_matrices

    subroutine output_matrix(lanczos, iconfig, irepeat, stem, matrix)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer, intent(in) :: iconfig, irepeat
        character(7), intent(in) :: stem
        character(2) :: ifmt, jfmt
        real(dp), intent(in) :: matrix(lanczos%nvecs, lanczos%nvecs)
        character(25) :: ind1, ind2, filename
        integer :: i, j, ilen, jlen, new_unit

        ! Create the filename.
        write(ind2,'(i15)') irepeat
        if (lanczos%tGround) then
            filename = trim(trim(stem)//'.'//trim(adjustl(ind2)))
        else if (lanczos%tFiniteTemp) then
            write(ind1,'(i15)') iconfig
            filename = trim(trim(stem)//'.'//trim(adjustl(ind1))//'.'//trim(adjustl(ind2)))
        end if

        new_unit = get_free_unit()
        open(new_unit, file=trim(filename), status='replace')

        ! Write all the components of the matrix, above and including the diagonal, one
        ! after another on separate lines. Each element is preceeded by the indices involved.
        do i = 1, lanczos%nvecs
            ilen = ceiling(log10(real(abs(i)+1)))
            ! ifmt will hold the correct integer length so that there will be no spaces printed out.
            ! Note that this assumes that ilen < 10, which is very reasonable!
            write(ifmt,'(a1,i1)') "i", ilen
            do j = i, lanczos%nvecs
                jlen = ceiling(log10(real(abs(j)+1)))
                write(jfmt,'(a1,i1)') "i", jlen

                ! Finally write the line.
                write(new_unit,'(a1,'//ifmt//',a1,'//jfmt//',a1,1x,es19.12)') &
                    "(",i,",",j,")", matrix(i,j)
            end do
        end do
        close(new_unit)

    end subroutine output_matrix

    subroutine print_populations_sl(lanczos)
    
        ! A useful test routine which will output the total walker population on both
        ! replicas, for each Lanczos vector.

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer :: ihash
        integer(n_int) :: int_sign(lenof_sign*lanczos%nvecs)
        real(dp) :: real_sign(lenof_sign*lanczos%nvecs), total_pop(lenof_sign*lanczos%nvecs)
        type(ll_node), pointer :: temp_node

        int_sign = 0_n_int
        total_pop = 0.0_dp
        real_sign = 0.0_dp
        
        do ihash = 1, nhashes_lanczos
            temp_node => lanczos_hash_table(ihash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    int_sign = lanczos_vecs(NIfD+1:NIfD+lenof_sign*lanczos%nvecs, temp_node%ind)
                    real_sign = transfer(int_sign, real_sign)
                    total_pop = total_pop + abs(real_sign)
                    temp_node => temp_node%next
                end do
            end if
        end do

        nullify(temp_node)

        write(6,*) "lanczos_vec populations:", total_pop

    end subroutine print_populations_sl

    subroutine print_amplitudes_sl(irepeat)

        ! A (*very* slow and memory intensive) test routine to print the current amplitudes (as stored
        ! in CurrentDets) of *all* determinants to a file. The amplitude of each replica will be printed
        ! one after the other. Since this is intended to be used with stochastic Lanczos, irepeat
        ! is the number of the current repeat, but it will simply be used in naming the output file.

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
            if (.not. IsAllowedHPHF(ilut(0:NIfD))) cycle
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

    end subroutine print_amplitudes_sl

end module stoch_lanczos_procs
