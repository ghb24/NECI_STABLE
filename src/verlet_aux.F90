#include "macros.h"

! module for auxiliary routines used in the third order verlet algorithm

module verlet_aux

    use constants, only: n_int, lenof_sign, dp, EPS, inum_runs, null_part, maxExcit, stdout

    use AnnihilationMod, only: DirectAnnihilation, SendProcNewParts, CompressSpawnedList

    use hash, only: clear_hash_table, hash_table_lookup, add_hash_table_entry

    use bit_rep_data, only: niftot, nifd, extract_sign, IlutBits

    use bit_reps, only: decode_bit_det, set_flag, get_initiator_flag_by_run, encode_sign, &
                        add_ilut_lists, extract_bit_rep, test_flag, encode_bit_rep

    use real_time_data, only: spawnBuf, spawnBufSize, dpsi_cache, dpsi_size, max_cache_size, &
                              backup_size, temp_det_list, tau_real, tau_imag, iterInit, tDynamicAlpha, tVerletSweep, &
                              runge_kutta_step

    use real_time_procs, only: attempt_die_realtime, real_time_determ_projection

    use SystemData, only: nel

    use FciMCData, only: CurrentDets, HashIndex, maxSpawned, iStartFreeSlot, iEndFreeSlot, &
                         inum_runs, SpawnedParts, TotWalkers, spawn_ht, iter_data_fciqmc, &
                         InitialSpawnedSlots, ValidSpawnedList, fcimc_excit_gen_store, &
                         FreeSlot, ll_node, &
                         popsfile_dets, WalkVecDets, exFlag, max_calc_ex_level, ilutRef, Hii, &
                         fcimc_iter_data, core_run
    use core_space_util, only: cs_replicas
    use CalcData, only: tTruncInitiator, AvMCExcits, tPairedReplicas, &
                        tSemiStochastic, tInitCoherentRule

    use procedure_pointers, only: attempt_create, attempt_die, generate_excitation, &
                                  encode_child

    use DetBitOps, only: FindBitExcitLevel

    use matel_getter, only: get_diagonal_matel, get_off_diagonal_matel

    use fcimc_pointed_fns, only: attempt_create_normal

    use fcimc_helper, only: CalcParentFlag, decide_num_to_spawn, create_particle_with_hash_table, &
                            SumEContrib, check_semistoch_flags, checkValidSpawnedList, rezero_iter_stats_each_iter

    use semi_stoch_procs, only: check_determ_flag

    use load_balance_calcnodes, only: DetermineDetNode

    use MPI_wrapper, only: iProcIndex

    use tau_search, only: tau, assign_value_to_tau

    implicit none

contains

    subroutine init_verlet_sweep()
        ! here, we initiate the first dpsi_cache from a runge-kutta calculation
        character(*), parameter :: this_routine = "init_verlet_sweep"

        write(stdout, *) "Prepared initial delta_psi, starting verlet calculation"
        call build_initial_delta_psi()
        ! rescale the timestep
        call assign_value_to_tau(iterInit * tau, this_routine)
        tau_imag = iterInit * tau_imag
        tau_real = iterInit * tau_real
        ! There is only one step now (we might log the second spawns as quasi-second
        ! step later on
        runge_kutta_step = 1

    end subroutine init_verlet_sweep

!-----------------------------------------------------------------------------------------------!

    subroutine check_verlet_sweep(iterRK)
        integer, intent(inout) :: iterRK
        ! if iterInit iterations were done in the runge-Kutta, switch to verlet scheme
        if (iterRK == iterInit) then
            tVerletSweep = .true.
            iterRK = 0
            call init_verlet_sweep()
        end if
    end subroutine check_verlet_sweep

!-----------------------------------------------------------------------------------------------!

    subroutine end_verlet_sweep()
        ! switch back to runge-kutta after adjusting alpha to get a new delta_psi
        write(stdout, *) "Switching to runge-kutta for update of alpha"
        tVerletSweep = .false.
        call assign_value_to_tau(tau / iterInit, 'end_verlet_sweep')
    end subroutine end_verlet_sweep

!-----------------------------------------------------------------------------------------------!

    subroutine build_initial_delta_psi()
        character(*), parameter :: this_routine = "build_initial_delta_psi"
        ! build delta_psi as  psi(delta_t) - psi(0), where psi(delta_t) is the current population
        ! and psi(0) the backup stored in popsfile_dets
        if (allocated(popsfile_dets)) then
            call add_ilut_lists(int(TotWalkers), backup_size, .true., CurrentDets, popsfile_dets, &
                                dpsi_cache, dpsi_size, -1.0_dp)

            ! we do not need popsfile dets anymore
            deallocate(popsfile_dets)
        else
            call stop_all(this_routine, "Backup buffer not allocated")
        end if
    end subroutine build_initial_delta_psi

!-----------------------------------------------------------------------------------------------!

    subroutine setup_delta_psi()
        ! temp_det_list is sufficient as a cache
        ! note that therefore, dpsi_cache also has a hashtable (temp_det_hash)
        dpsi_cache => temp_det_list
        max_cache_size = size(WalkVecDets, dim=2)
    end subroutine setup_delta_psi

!-----------------------------------------------------------------------------------------------!

    subroutine backup_initial_state()
        character(*), parameter :: this_routine = "backup_initial_state"

        ! popsfile_dets is certainly large enough to house CurrentDets, so use it if available
        if (.not. allocated(popsfile_dets)) allocate(popsfile_dets(0:niftot, TotWalkers))
        popsfile_dets(:, 1:TotWalkers) = CurrentDets(:, 1:TotWalkers)
        ! number of initial walkers
        backup_size = int(TotWalkers)
    end subroutine backup_initial_state

    subroutine init_verlet_iteration()
        use rdm_data, only: rdm_definitions_t
        type(rdm_definitions_t) :: dummy

        FreeSlot(1:iEndFreeSlot) = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0

        ! verlet always uses a spawn hashtable
        call clear_hash_table(spawn_ht)
        call rezero_iter_stats_each_iter(iter_data_fciqmc, dummy)

    end subroutine init_verlet_iteration

!-----------------------------------------------------------------------------------------------!

    ! applies the hamiltonian twice to the current population and stores the result
    ! in spawnedParts
    subroutine obtain_h2_psi()
        real(dp) :: dummy_sign(lenof_sign)

        ! apply H once, we now have the spawnedParts from a single iteration
        call apply_hamiltonian(CurrentDets, int(TotWalkers), .true., tTruncInitiator, .true.)

        ! communicate the spawns between processors and store the compressed spawns into a buffer
        call generate_spawn_buf()

        ! apply H to the buffer to get H^2 on the original population. The result is stored
        ! in spawnedParts (uncommunicated)
        call apply_hamiltonian(spawnBuf, spawnBufSize, .false., .false., .false.)

        ! communicate the result and compress the population (such that each determinant
        ! only occurs once)
        call SendProcNewParts(spawnBufSize, .false.)
        call CompressSpawnedList(spawnBufSize, iter_data_fciqmc)

    end subroutine obtain_h2_psi

!-----------------------------------------------------------------------------------------------!

    subroutine apply_hamiltonian(population, popsize, tGetFreeSlots, tGetInitFlags, tSumE)
        ! this subroutine performs spawning from population to spawnVec by applying
        ! delta t H once. No annihilation is performed and no other steps performed
        ! keep it minimalistic and stick to the SRP principle

        ! by default writes into SpawnedParts
        integer(n_int), intent(in) :: population(0:, 1:)
        integer, intent(in) :: popsize
        ! store the free slots in population
        logical, intent(in) :: tGetInitflags, tGetFreeSlots, tSumE
        integer :: idet, nI(nel), determ_index, unused_flags, ex_level
        real(dp) :: parent_sign(lenof_sign), hdiag
        logical :: tEmptyDet, tCoreDet
        HElement_t(dp) :: hoffdiag

        ! where to put this?
        ! attempt_create => attempt_create_normal

        ValidSpawnedList = InitialSpawnedSlots
        ! we do not know what is in spawn_ht, but we clearly want it to be empty now
        call clear_hash_table(spawn_ht)

        determ_index = 1

        do idet = 1, popsize
            ! apply spawn and death for each walker
            fcimc_excit_gen_store%tFilled = .false.
            unused_flags = 0

            call extract_bit_rep(population(:, idet), nI, parent_sign, unused_flags, idet, &
                                 fcimc_excit_gen_store)

            tEmptyDet = IsUnoccDet(parent_sign)
            ! we collect free slots only in the first application
            if (tEmptyDet) then
                if (tGetFreeSlots) then
                    iEndFreeSlot = iEndFreeSlot + 1
                    FreeSlot(iEndFreeSlot) = idet
                end if
                cycle
            end if

            ! also, initiator flags are only reset in the first iteration
            ! when using the initiator criterium
            ! else, we just let the flags be
            hdiag = 0.0_dp
            if (tGetInitFlags) call CalcParentFlag(idet, unused_flags)
            hdiag = get_diagonal_matel(nI, population(:, idet)) - Hii
            ! if desired, sum in the energy contribution
            if (tSumE) then
                ex_level = FindBitExcitLevel(ilutRef, population(:, idet), max_calc_ex_level)
                hoffdiag = get_off_diagonal_matel(nI, population(:, idet))
                call SumEContrib(nI, ex_level, parent_sign, population(:, idet), &
                                 hdiag, hoffdiag, 1.0_dp, tPairedReplicas, idet)
            end if

            tCoreDet = check_determ_flag(population(:, idet))

            if (tCoreDet) then
                associate(rep => cs_replicas(core_run))
                    rep%indices_of_determ_states(determ_index) = idet
                    rep%partial_determ_vecs(:, determ_index) = parent_sign
                end associate
                determ_index = determ_index + 1
                if (IsUnoccDet(parent_sign)) cycle
            end if

            ! the initiator flags are set upon the first iteration

            call perform_spawn(population(:, idet), nI, parent_sign, hdiag, tCoreDet)
        end do

        if (tSemiStochastic) call real_time_determ_projection()
    end subroutine apply_hamiltonian

!-----------------------------------------------------------------------------------------------!

    subroutine perform_spawn(iLut_parent, nI, parent_sign, hdiag, tCoreDet)
        real(dp), intent(in) :: parent_sign(lenof_sign), hdiag
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut_parent(0:niftot)
        logical, intent(in) :: tCoreDet
        integer :: part, nspawn, ispawn, nI_child(nel), ic, ex(2, maxExcit), unused_ex_level
        integer(n_int) :: ilut_child(0:niftot)!, ilut_parent_init(0:niftot)
        real(dp) :: prob, child_sign(lenof_sign), unused_rdm_real, unused_sign(nel)
        real(dp) :: unused_precond_fac, unused_avEx
        logical :: tParity, break
        integer :: err
        HElement_t(dp) :: HElGen

        unused_ex_level = 0
        do part = 1, lenof_sign
            call decide_num_to_spawn(parent_sign(part), AvMCExcits, nspawn)
            do ispawn = 1, nspawn
                ilut_child = 0_n_int

                call generate_excitation(nI, iLut_parent, nI_child, ilut_child, exFlag, ic, &
                                         ex, tParity, prob, HElGen, fcimc_excit_gen_store)

                if (.not. IsNullDet(nI_child)) then
                    call encode_child(ilut_parent, ilut_child, ic, ex)
                    ilut_child(IlutBits%ind_flag) = 0_n_int

                    ! treating semi-stochastic space
                    ! note that diagonal event either are not in the core space at all
                    ! or from the core space to the core space, i.e. they are covered fully
                    ! by this branching
                    if (tSemiStochastic) then
                        break = check_semistoch_flags(ilut_child, nI_child, &
                                                      part_type_to_run(part), tCoredet)
                        if (break) cycle
                    end if

                    child_sign = attempt_create(nI, iLut_parent, parent_sign, nI_child, iLut_child, &
                                                prob, HElGen, ic, ex, tParity, unused_ex_level, part, &
                                                unused_sign, unused_avEx, unused_rdm_real, unused_precond_fac)
                else
                    child_sign = 0.0_dp
                end if

                if ((any(abs(child_sign) > EPS)) .and. (ic /= 0) .and. (ic <= 2)) then
                    call create_particle_with_hash_table(nI_child, ilut_child, child_sign, &
                                                         part, ilut_parent, iter_data_fciqmc, err)
                    if (err /= 0) return
                end if ! If a child was spawned.

            end do ! Over mulitple particles on same determinant.
        end do
        if (.not. tCoreDet) then
            child_sign = -attempt_die_realtime(hdiag, parent_sign, unused_ex_level)
            if (any(abs(child_sign) > EPS)) then
                ! diagonal events are treated the same way as offdiagonal ones,
                ! except for an extra flag indicating the diagonal spawn
                !ilut_parent_init = ilut_parent
                !call set_flag(ilut_parent_init,get_initiator_flag(part_type))
                ! we have to think about whether diagonal spawns get special treatment
                ! I think they should be subject to the same rules as offdiagonal spawns
                call create_diagonal_with_hashtable(nI, iLut_parent, child_sign, err)
                if (err /= 0) return
            end if
        end if

    end subroutine perform_spawn

!-----------------------------------------------------------------------------------------------!

    subroutine create_diagonal_with_hashtable(nI, iLut, sign, err)
        ! this subroutine is somewhat a variant of create_particle_with_hashtable
        ! that takes a global sign with many entries as it appears in diagonal spawning events
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: iLut(0:niftot)
        real(dp), intent(in) :: sign(lenof_sign)
        integer, intent(out) :: err
        integer :: ind, hash_val, run, proc
        integer, parameter :: flags = 0
        logical :: tSuccess
        real(dp) :: old_sign(lenof_sign)
        character(*), parameter :: this_routine = "create_diagonal_with_hashtable"

        err = 0

        call hash_table_lookup(nI, iLut, nifd, spawn_ht, SpawnedParts, ind, hash_val, tSuccess)
        if (tSuccess) then
            ! if it already exists, add in the
            call extract_sign(SpawnedParts(:, ind), old_sign)
            call encode_sign(SpawnedParts(:, ind), old_sign + sign)

            ! check for initiator criterium
            if (tTruncInitiator) then
                do run = 1, inum_runs
                    if (tInitCoherentRule) then
                        if (.not. is_run_unnocc(old_sign, run) .or. test_flag( &
                            ilut, get_initiator_flag_by_run(run))) then
                            call set_flag(SpawnedParts(:, ind), get_initiator_flag_by_run(run))
                        end if
                    else
                        if (test_flag(ilut, get_initiator_flag_by_run(run))) then
                            call set_flag(SpawnedParts(:, ind), get_initiator_flag_by_run(run))
                        end if
                    end if
                end do
            end if
        else
            proc = DetermineDetNode(nel, nI, 0)

            if (checkValidSpawnedList(proc)) then
                err = 1
                return
            end if

            call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), ilut(0:nifd), &
                                sign, flags)

            if (tTruncInitiator) then
                do run = 1, inum_runs
                    if (test_flag(ilut, get_initiator_flag_by_run(run))) call set_flag( &
                        SpawnedParts(:, ValidSpawnedList(proc)), get_initiator_flag_by_run(run))
                end do
            end if

            call add_hash_table_entry(spawn_ht, ValidSpawnedList(proc), hash_val)
            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        end if
    end subroutine create_diagonal_with_hashtable

!-----------------------------------------------------------------------------------------------!

    subroutine generate_spawn_buf()
        character(*), parameter :: this_routine = "generate_spawn_buf"

        call SendProcNewParts(spawnBufSize, .false.)
        call CompressSpawnedList(spawnBufSize, iter_data_fciqmc)

        if (tSemiStochastic) call add_semistoch_spawns(spawnedParts, spawnBufSize, spawn_ht, &
                                                       maxSpawned, CurrentDets)

        spawnBuf(:, 1:spawnBufSize) = SpawnedParts(:, 1:spawnBufSize)
    end subroutine generate_spawn_buf

!-----------------------------------------------------------------------------------------------!

    subroutine add_semistoch_spawns(population, populationSize, hashTable, maxSize, &
                                    sourcePopulation)
        integer(n_int), intent(inout) :: population(:, :)
        integer, intent(inout) :: populationSize
        type(ll_node), pointer, intent(inout) :: hashTable(:)
        integer, intent(in) :: maxSize
        integer(n_int), intent(in) :: sourcePopulation(:, :)

        integer :: i, hashValue, ilutindex, nI(nel)
        real(dp) :: sign(lenof_sign)
        logical :: tSuccess
        character(*), parameter :: this_routine = "add_semistoch_spawns"

        associate(rep => cs_replicas(core_run))
        do i = 1, rep%determ_sizes(iProcIndex)
            ! check if the core-space determinant was already spawned upon
            call decode_bit_det(nI, sourcePopulation(:, rep%indices_of_determ_states(i)))
            call hash_table_lookup(nI, sourcePopulation(:, rep%indices_of_determ_states(i)), nifd, &
                                   hashTable, population, ilutindex, hashValue, tSuccess)

            if (tSuccess) then
                ! if it is found, add the signs
                call extract_sign(population(:, ilutIndex), sign)
                call encode_sign(population(:, ilutIndex), sign + rep%partial_determ_vecs(:, i))
            else
                ! the spawn is new, add it to population
                populationSize = populationSize + 1
                if (populationSize > maxSize) call stop_all(this_routine, &
                                                            "Out of memory for adding semistochastic spawns")
                population(:, populationSize) = sourcePopulation(:, rep%indices_of_determ_states(i))
                call encode_sign(population(:, populationSize), rep%partial_determ_vecs(:, i))
                ! add the hash table entry for the new determinant
                call add_hash_table_entry(hashTable, populationSize, hashValue)
            end if
        end do
        end associate

    end subroutine add_semistoch_spawns

!-----------------------------------------------------------------------------------------------!

    subroutine merge_ilut_lists(listA, listB, hashTable, sizeA, sizeB, maxSizeA)
        integer(n_int), intent(in) :: listB(0:, 1:)
        integer(n_int), intent(inout) :: listA(0:, 1:)
        integer, intent(in) :: sizeB, maxSizeA
        integer, intent(inout) :: sizeA
        type(ll_node), pointer, intent(inout) :: hashTable(:)
        integer :: nJ(nel), ilutIndex, hashValue, i, insertPos
        real(dp) :: signA(lenof_sign), signB(lenof_sign)
        logical :: tSuccess
        character(*), parameter :: this_routine = "merge_ilut_lists"

        ! this merges listB into listA. In the current form, empty slots in listA are
        ! not exploited, because it should not make a difference currently (optimization follows)
        do i = 1, sizeB
            call decode_bit_det(nJ, listB(:, i))

            call hash_table_lookup(nJ, listB(:, i), nifd, hashTable, listA, ilutIndex, &
                                   hashValue, tSuccess)

            if (tSuccess) then
                ! the i-th determinant in listB is already present in listA
                ! -> add up the signs
                call extract_sign(listA(:, ilutIndex), signA)
                call extract_sign(listB(:, i), signB)

                ! we do not fill up empty slots, so we do not care if signA == 0
                ! this differs from AnnihilateSpawnedParts
                call encode_sign(listA(:, ilutIndex), signA + signB)

                ! the initiator criterium is checked upon annihilation, no need to do so here
            else
                ! if the entry in listB is not in listA, add it as the last entry
                sizeA = sizeA + 1
                insertPos = sizeA
                if (sizeA > maxSizeA) &
                    call stop_all(this_routine, "Out of memory for merging ilut lists")
                listA(:, insertPos) = listB(:, i)
                call add_hash_table_entry(hashTable, insertPos, hashValue)
            end if
        end do
    end subroutine merge_ilut_lists

!-----------------------------------------------------------------------------------------------!

    subroutine update_delta_psi()
        ! here, we add H^2 psi to delta_psi to generate the new delta_psi
        ! this might cause trouble as stochastic error is not averaged out,
        ! but we just carry over the error from the last iteration
        character(*), parameter :: this_routine = "update_delta_psi"

        ! add up delta_psi from the last iteration and spawnedParts (i.e. H^2 psi)
        call merge_ilut_lists(spawnedParts, dpsi_cache, spawn_ht, spawnBufSize, &
                              dpsi_size, maxSpawned)
        ! if semistochastic mode is run, we add the semistochastic spawns now
        ! because they have to be included in delta psi

        if (tSemiStochastic) call add_semistoch_spawns(spawnedParts, spawnBufSize, &
                                                       spawn_ht, maxSpawned, spawnBuf)

        ! cache delta_psi for the next iteration
        if (spawnBufSize > max_cache_size) call stop_all(this_routine, &
                                                         "Insufficient memory for creating delta_psi")

        dpsi_cache(:, 1:spawnBufSize) = spawnedParts(:, 1:spawnBufSize)
        dpsi_size = spawnBufSize

        ! for now, consider all entries in delta_psi as safe spawns for the next iteration
        if (tTruncInitiator) call set_initiator_flags_array(dpsi_cache, dpsi_size)
    end subroutine update_delta_psi

!-----------------------------------------------------------------------------------------------!

    subroutine set_initiator_flags_array(list, listSize)
        integer(n_int), intent(inout) :: list(0:, 1:)
        integer, intent(in) :: listSize
        integer :: i, run

        do i = 1, listSize
            do run = 1, inum_runs
                call set_flag(list(:, i), get_initiator_flag_by_run(run))
            end do
        end do
    end subroutine set_initiator_flags_array

end module verlet_aux
