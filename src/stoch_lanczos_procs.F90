#include "macros.h"

module stoch_lanczos_procs

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use bit_rep_data
    use bit_reps, only: decode_bit_det
    use CalcData, only: tTruncInitiator, tStartSinglePart, InitialPart, InitWalkers
    use CalcData, only: tSemiStochastic, tReadPops, tUseRealCoeffs, tau, DiagSft
    use constants
    use DetBitOps, only: DetBitEq
    use dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData, only: ilutHF, HFDet, CurrentDets, SpawnedParts, SpawnedParts2, TotWalkers
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, HashIndex, nWalkerHashes
    use FciMCData, only: fcimc_iter_data, ll_node, MaxWalkersPart, tStartCoreGroundState
    use FciMCData, only: tPopsAlreadyRead, tHashWalkerList, CurrentH, determ_proc_sizes
    use FciMCData, only: core_ham_diag, InputDiagSft, Hii
    use FciMCParMod, only: create_particle, InitFCIMC_HF, SetupParameters, InitFCIMCCalcPar
    use FciMCParMod, only: init_fcimc_fn_pointers, WriteFciMCStats, WriteFciMCStatsHeader
    use FciMCParMod, only: rezero_iter_stats_each_iter, tSinglePartPhase
    use hash, only: FindWalkerHash, init_hash_table, reset_hash_table, fill_in_hash_table
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias
    use Parallel_neci, only: MPIBarrier, iProcIndex, MPISum
    use ParallelHelper, only: root
    use procedure_pointers
    use semi_stoch_procs, only: copy_core_dets_this_proc_to_spawnedparts
    use semi_stoch_procs, only: add_core_states_currentdet_hash, start_walkers_from_core_ground
    use SystemData, only: nel
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
        NIfLan = NIfD + lenof_sign*lanczos%nvecs + NIfFlag

        nhashes_lanczos = nWalkerHashes
        TotWalkersLanczos = 0
        allocate(lanczos_vecs(0:NIfLan, MaxWalkersPart), stat=ierr)
        lanczos_vecs = 0_n_int
        allocate(lanczos_hash_table(nhashes_lanczos), stat=ierr)
        call init_hash_table(lanczos_hash_table)

        allocate(lanczos%overlap_matrix_1(lanczos%nvecs, lanczos%nvecs), stat=ierr)
        allocate(lanczos%overlap_matrix_2(lanczos%nvecs, lanczos%nvecs), stat=ierr)
        allocate(lanczos%hamil_matrix(lanczos%nvecs, lanczos%nvecs), stat=ierr)

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
        lanczos%overlap_matrix_1 = 0.0_dp
        lanczos%overlap_matrix_2 = 0.0_dp
        lanczos%hamil_matrix = 0.0_dp
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

    subroutine store_lanczos_vec(ivec, nvecs)

        integer, intent(in) ::  ivec, nvecs
        integer :: idet, sign_ind, flag_ind, DetHash, det_ind
        integer :: nI(nel)
        logical :: tDetFound
        type(ll_node), pointer :: temp_node, prev

        ! The index of the first element referring to the sign, for this ivec.
        sign_ind = NIfD + lenof_sign*(ivec-1) + 1
        ! The index of the flag, which is *not* NOffFlag!
        flag_ind = NIfD + lenof_sign*nvecs + 1

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
                        ! Just add the sign for the new Lanczos vector. The determinant and flag are
                        ! there already.
                        lanczos_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,idet)
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

                ! Copy determinant across.
                lanczos_vecs(0:NIfD,det_ind) = CurrentDets(0:NIfD,idet)
                ! Copy sign across.
                lanczos_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,idet)
                ! Copy Flags across.
                if (tUseFlags) lanczos_vecs(flag_ind,det_ind) = CurrentDets(NOffFlag,idet)
            end if

            nullify(temp_node)
            nullify(prev)
        end do

    end subroutine store_lanczos_vec

    subroutine calc_overlap_matrix_elems(lanczos, ivec)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec), sgn(lenof_sign)
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

            ! Loop over all walkers in lanczos_vecs.
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

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: ivec
        integer :: idet, jvec, ind(ivec), sgn(lenof_sign), nI(nel)
        integer :: det_ind, DetHash
        real(dp) :: sign1(lenof_sign), sign2(lenof_sign), full_shift(lenof_sign)
        type(ll_node), pointer :: temp_node
        logical :: tDetFound

        associate(h_matrix => lanczos%hamil_matrix, &
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

        end associate

    end subroutine calc_hamil_elems

    subroutine communicate_lanczos_matrices(lanczos)

        ! Add all the overlap and projected Hamiltonian matrices together, with the result being
        ! held only on the root node.

        type(stoch_lanczos_data), intent(inout) :: lanczos
        real(dp) :: inp_matrices(3*lanczos%nvecs, lanczos%nvecs)
        real(dp) :: out_matrices(3*lanczos%nvecs, lanczos%nvecs)

        inp_matrices(1:lanczos%nvecs, 1:lanczos%nvecs) = lanczos%overlap_matrix_1
        inp_matrices(lanczos%nvecs+1:2*lanczos%nvecs, 1:lanczos%nvecs) = lanczos%overlap_matrix_2
        inp_matrices(2*lanczos%nvecs+1:3*lanczos%nvecs, 1:lanczos%nvecs) = lanczos%hamil_matrix

        call MPISum(inp_matrices, out_matrices)

        if (iProcIndex == root) then
            lanczos%overlap_matrix_1 = out_matrices(1:lanczos%nvecs, 1:lanczos%nvecs)
            lanczos%overlap_matrix_2 = out_matrices(lanczos%nvecs+1:2*lanczos%nvecs, 1:lanczos%nvecs)
            lanczos%hamil_matrix = out_matrices(2*lanczos%nvecs+1:3*lanczos%nvecs, 1:lanczos%nvecs)
        end if

    end subroutine communicate_lanczos_matrices

    subroutine output_lanczos_matrices(lanczos, irepeat)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer, intent(in) :: irepeat
        real(dp) :: average_s_matrix(lanczos%nvecs, lanczos%nvecs)

        if (iProcIndex == root) then
            call output_matrix(lanczos, irepeat, 'hamil  ', lanczos%hamil_matrix)
            ! We have two overlap matrices as we have two replicas. So average them for
            ! better statistics.
            average_s_matrix = (lanczos%overlap_matrix_1 + lanczos%overlap_matrix_2)/2.0_dp
            call output_matrix(lanczos, irepeat, 'overlap', average_s_matrix)
        end if

    end subroutine output_lanczos_matrices

    subroutine output_matrix(lanczos, irepeat, stem, matrix)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer, intent(in) :: irepeat
        character(7), intent(in) :: stem
        character(2) :: ifmt, jfmt
        real(dp), intent(in) :: matrix(lanczos%nvecs, lanczos%nvecs)
        character(15) :: ind, filename
        integer :: i, j, ilen, jlen, new_unit

        write(ind,'(i15)') irepeat
        filename = trim(trim(stem)//'.'//adjustl(ind))

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
                write(new_unit,'(a1,'//ifmt//',a1,'//jfmt//',a1,1x,es15.8)') &
                    "(",i,",",j,")", matrix(i,j)
            end do
        end do
        close(new_unit)

    end subroutine output_matrix

    subroutine print_populations(lanczos)
    
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

        write(6,*) "lanczos_vec populations", total_pop

    end subroutine print_populations

end module stoch_lanczos_procs
