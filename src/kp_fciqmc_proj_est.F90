! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

module kp_fciqmc_proj_est

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use kp_fciqmc_data_mod, only: lenof_all_signs, kp_ind_1, kp_ind_2

    implicit none

contains

    subroutine calc_projected_hamil(nvecs, krylov_array, krylov_ht, ndets, h_matrix, partial_vecs, full_vecs, h_diag)

        use bit_rep_data, only: NOffFlag, tUseFlags
        use bit_reps, only: decode_bit_det
        use CalcData, only: ss_space_in, tDetermHFSpawning
        use constants
        use DetBitOps, only: return_ms, FindBitExcitLevel
        use enumerate_excitations, only: generate_connection_normal, generate_connection_kpnt
        use enumerate_excitations, only: init_generate_connected_space
        use FciMCData, only: fcimc_excit_gen_store, exFlag, SpawnVecKP, SpawnVecKP2
        use FciMCData, only: SpawnVec, SpawnVec2, determ_sizes, determ_displs, SpawnedParts
        use FciMCData, only: SpawnedParts2, InitialSpawnedSlots, ValidSpawnedList, ll_node
        use FciMCData, only: spawn_ht, subspace_hamil_time, iLutHF_True, max_calc_ex_level
        use FCIMCData, only: HFConn, Hii
        use global_det_data, only: det_diagH
        use hash, only: clear_hash_table
        use kp_fciqmc_data_mod, only: tSemiStochasticKPHamil, tExcitedStateKP, av_mc_excits_kp
        use kp_fciqmc_data_mod, only: tExactHamilSpawning, kp_hamil_exact_frac
        use procedure_pointers, only: generate_excitation, encode_child, get_spawn_helement
        use semi_stoch_procs, only: is_core_state, check_determ_flag
        use semi_stoch_procs, only: determ_projection_kp_hamil
        use SystemData, only: nbasis, tAllSymSectors, nOccAlpha, nOccBeta, nel
        use SystemData, only: tKPntSym, tUseBrillouin
        use timing_neci, only: set_timer, halt_timer
        use util_mod, only: stochastic_round
        
        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        type(ll_node), pointer, intent(in) :: krylov_ht(:)
        integer, intent(in) :: ndets
        real(dp), intent(out) :: h_matrix(:,:)
        real(dp), allocatable, intent(inout), optional :: partial_vecs(:,:), full_vecs(:,:)
        real(dp), intent(in), optional :: h_diag(:)

        integer :: idet, ispawn, nspawn, i, j, ex_level_to_hf
        integer :: determ_ind, flag_ind, ic, ex(2,2), ms_parent
        integer :: ex_flag, nStore(6)
        integer, allocatable :: excit_gen(:)
        integer :: nI_parent(nel), nI_child(nel)
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        real(dp) :: prob, tot_pop, h_diag_elem
        real(dp), dimension(lenof_all_signs) :: child_sign, parent_sign
        integer(n_int) :: int_sign(lenof_all_signs)
        logical :: tChildIsDeterm, tParentIsDeterm, tParentUnoccupied, tParity
        logical :: tNearlyFull, tFinished, tAllFinished, all_excits_found, tTempUseBrill
        HElement_t :: HElGen, HEl

        call set_timer(subspace_hamil_time)

        h_matrix(:,:) = 0.0_dp
        ilut_parent = 0_n_int
        if (tUseFlags) flag_ind = NIfDBO + lenof_all_signs + 1
        ValidSpawnedList = InitialSpawnedSlots
        call clear_hash_table(spawn_ht)
        tNearlyFull = .false.
        tFinished = .false.
        determ_ind = 1

        if (.not. tExcitedStateKP) then
            ! We want SpawnedParts and SpawnedParts2 to point to the 'wider' spawning arrays which
            ! hold all the signs of *all* the Krylov vectors. This is because SendProcNewParts uses
            ! SpawnedParts and SpawnedParts2.
            SpawnedParts => SpawnVecKP
            SpawnedParts2 => SpawnVecKP2
        end if

        do idet = 1, ndets

            ! The 'parent' determinant from which spawning is to be attempted.
            ilut_parent(0:NIfDBO) = krylov_array(0:NIfDBO,idet)
            if (tUseFlags) ilut_parent(NOffFlag) = krylov_array(flag_ind,idet)

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call decode_bit_det(nI_parent, ilut_parent)
            int_sign = krylov_array(NIfDBO+1:NIfDBO+lenof_all_signs, idet)
            parent_sign = transfer(int_sign, parent_sign)
            tot_pop = sum(abs(parent_sign))

            tParentIsDeterm = check_determ_flag(ilut_parent)
            tParentUnoccupied = all(abs(parent_sign) < 1.e-16_dp)

            ! If this determinant is in the deterministic space then store the relevant
            ! data in arrays for later use.
            if (tParentIsDeterm) then
                ! Add the amplitude to the deterministic vector.
                partial_vecs(:,determ_ind) = parent_sign
                determ_ind = determ_ind + 1
            end if

            if (tParentUnoccupied) cycle

            if (ss_space_in%tDoubles .and. tDetermHFSpawning) then
                ex_level_to_hf = FindBitExcitLevel (iLutHF_true, ilut_parent, max_calc_ex_level)
                if (ex_level_to_hf == 0) cycle
            end if

            if (tAllSymSectors) then
                ms_parent = return_ms(ilut_parent)
                ! If this condition is met (if all electrons have spin up or all have spin
                ! down) then there will be no determinants to spawn to, so don't attempt spawning.
                if (abs(ms_parent) == nbasis/2) cycle
                nOccAlpha = (nel+ms_parent)/2
                nOccBeta = (nel-ms_parent)/2
            end if

            if (tExactHamilSpawning .and. av_mc_excits_kp*tot_pop >= real(HFConn,dp)*kp_hamil_exact_frac) then

                ! We no longer want to treat these deterministic states by the
                ! usual deterministic projection. So set the corresponding
                ! amplitude in the deterministic vector to 0 and add in the
                ! diagonal element separately.
                if (tParentIsDeterm) then
                    partial_vecs(:,determ_ind-1) = 0.0_dp

                    if (present(h_diag)) then
                        h_diag_elem = h_diag(idet) + Hii
                    else
                        h_diag_elem = det_diagH(idet) + Hii
                    end if

                    child_sign = calc_amp_kp_hamil(parent_sign, 1.0_dp, 1.0_dp, h_diag_elem)
                    call create_particle_kp_estimates(nI_parent, ilut_parent, child_sign, tNearlyFull)
                end if

                call init_generate_connected_space(nI_parent, ex_flag, all_excits_found, ex, excit_gen, &
                                                    nstore, tTempUseBrill)

                do while(.true.)
                    if (tKPntSym) then
                        call generate_connection_kpnt(nI_parent, ilut_parent, nI_child, ilut_child, ex_flag, &
                                                         all_excits_found, nstore, excit_gen, hel=HEl)
                    else
                        call generate_connection_normal(nI_parent, ilut_parent, nI_child, ilut_child, ex_flag, &
                                                         ex, all_excits_found, hel=HEl)
                    end if

                    if (all_excits_found) exit

                    child_sign = calc_amp_kp_hamil(parent_sign, 1.0_dp, 1.0_dp, real(HEl,dp))
                    call create_particle_kp_estimates(nI_child, ilut_child, child_sign, tNearlyFull)

                    if (tNearlyFull) then
                        call add_in_contribs(nvecs, krylov_array, krylov_ht, tFinished, tAllFinished, 1.0_dp, h_matrix)
                        call clear_hash_table(spawn_ht)
                        tNearlyFull = .false.
                    end if
                end do

                if (tKPntSym) then
                    tUseBrillouin = tTempUseBrill
                    deallocate(excit_gen)
                end if

            else

                nspawn = stochastic_round(av_mc_excits_kp*tot_pop)
                
                do ispawn = 1, nspawn

                    ! Zero the bit representation, to ensure no extraneous data gets through.
                    ilut_child = 0_n_int

                    call generate_excitation (nI_parent, ilut_parent, nI_child, &
                                        ilut_child, exFlag, ic, ex, tParity, prob, &
                                        HElGen, fcimc_excit_gen_store)

                    ! If a valid excitation.
                    if (.not. IsNullDet(nI_child)) then

                        call encode_child (ilut_parent, ilut_child, ic, ex)

                        ! If spawning is both too and from the core space, abort it.
                        if (tSemiStochasticKPHamil .and. tParentisDeterm) then
                            if(is_core_state(ilut_child, nI_child)) cycle
                        end if

                        HEl = get_spawn_helement(nI_parent, nI_child, ilut_parent, ilut_child, ic, ex, &
                                                 tParity, HElGen)

                        child_sign = calc_amp_kp_hamil(parent_sign, prob, av_mc_excits_kp*tot_pop, real(HEl,dp))
                    else
                        child_sign = 0.0_dp
                    end if

                    ! If any (valid) children have been spawned.
                    if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                        call create_particle_kp_estimates(nI_child, ilut_child, child_sign, tNearlyFull)

                        if (tNearlyFull) then
                            call add_in_contribs(nvecs, krylov_array, krylov_ht, tFinished, tAllFinished, 1.0_dp, h_matrix)
                            call clear_hash_table(spawn_ht)
                            tNearlyFull = .false.
                        end if

                    end if ! If a child was spawned.

                end do ! Over mulitple spawns from the same determinant.
            end if

        end do ! Over all determinants.

        tFinished = .true.
        do
            call add_in_contribs(nvecs, krylov_array, krylov_ht, tFinished, tAllFinished, 1.0_dp, h_matrix)
            if (tAllFinished) exit
        end do

        call calc_hamil_contribs_diag(nvecs, krylov_array, ndets, h_matrix, h_diag)

        if (tSemiStochasticKPHamil) then
            call determ_projection_kp_hamil(partial_vecs, full_vecs, determ_sizes, determ_displs)
            call calc_hamil_contribs_semistoch(nvecs, krylov_array, h_matrix, partial_vecs)
        end if

        ! Symmetrise the projected Hamiltonian.
        do i = 1, nvecs
            do j = 1, i-1
                h_matrix(i,j) = h_matrix(j,i)
            end do
        end do

        if (.not. tExcitedStateKP) then
            ! Now let SpawnedParts and SpawnedParts2 point back to their
            ! original arrays.
            SpawnedParts => SpawnVec
            SpawnedParts2 => SpawnVec2
        end if

        call halt_timer(subspace_hamil_time)

    end subroutine calc_projected_hamil

    subroutine calc_projected_spin(nvecs, krylov_array, krylov_ht, ndets, spin_matrix)

        use bit_rep_data, only: NOffSgn
        use bit_reps, only: decode_bit_det
        use DetBitOps, only: EncodeBitDet
        use FciMCData, only: spawn_ht, SpawnVecKP, SpawnVecKP2, SpawnVec, SpawnVec2
        use FciMCData, only: ll_node, SpawnedParts, SpawnedParts2, InitialSpawnedSlots
        use FciMCData, only: ValidSpawnedList, subspace_spin_time
        use get_excit, only: make_double
        use hash, only: clear_hash_table
        use kp_fciqmc_data_mod, only: tExcitedStateKP
        use SystemData, only: nel
        use timing_neci, only: set_timer, halt_timer

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        type(ll_node), pointer, intent(in) :: krylov_ht(:)
        integer, intent(in) :: ndets
        real(dp), intent(out) :: spin_matrix(:,:)

        integer :: idet, iel, jel, iorb, jorb, i, j
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        integer :: nI_parent(nel), nI_child(nel)
        logical :: ialpha, jalpha, ibeta, jbeta
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: parent_sign(lenof_all_signs), child_sign(lenof_all_signs)
        logical :: tNearlyFull, tFinished, tAllFinished, tParity
        integer :: ex(2,2)

        call set_timer(subspace_spin_time)

        spin_matrix(:,:) = 0.0_dp
        ilut_parent = 0_n_int
        ValidSpawnedList = InitialSpawnedSlots
        call clear_hash_table(spawn_ht)
        tNearlyFull = .false.
        tFinished = .false.

        if (.not. tExcitedStateKP) then
            ! We want SpawnedParts and SpawnedParts2 to point to the 'wider' spawning arrays which
            ! hold all the signs of *all* the Krylov vectors. This is because SendProcNewParts uses
            ! SpawnedParts and SpawnedParts2.
            SpawnedParts => SpawnVecKP
            SpawnedParts2 => SpawnVecKP2
        end if

        ! Loop over all states and add all contributions.
        do idet = 1, ndets

            ! The 'parent' determinant to consider 'spawning' from.
            ilut_parent(0:NIfDBO) = krylov_array(0:NIfDBO,idet)

            ! Get the real walker amplitudes.
            int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, idet)
            parent_sign = transfer(int_sign, parent_sign)

            call decode_bit_det(nI_parent, ilut_parent)

            do iel = 1, nel
                iorb = nI_parent(iel)
                ibeta = is_beta(iorb)
                ialpha = .not. ibeta

                do jel = iel+1, nel

                    jorb = nI_parent(jel)
                    jbeta = is_beta(jorb)
                    jalpha = .not. jbeta

                    ! If both are beta orbitals, or both are alpha orbitals,
                    ! then add the contributions to spin matrix.
                    if ((ibeta .and. jbeta) .or. (ialpha .and. jalpha)) then
                        do i = 1, nvecs
                            do j = i, nvecs
                                spin_matrix(i,j) = spin_matrix(i,j) + &
                                    0.5_dp*(parent_sign(kp_ind_1(i))*parent_sign(kp_ind_2(j)) + &
                                            parent_sign(kp_ind_2(i))*parent_sign(kp_ind_1(j)))/2.0_dp
                            end do
                        end do
                    else
                        ! One orbital has alpha spin, the other has beta spin.

                        ! First add the elements corresponding to 'spawning'
                        ! from and to the same determinant.
                        do i = 1, nvecs
                            do j = i, nvecs
                                spin_matrix(i,j) = spin_matrix(i,j) - &
                                    0.5_dp*(parent_sign(kp_ind_1(i))*parent_sign(kp_ind_2(j)) + &
                                            parent_sign(kp_ind_2(i))*parent_sign(kp_ind_1(j)))/2.0_dp
                            end do
                        end do

                        ! Now for the more tricky bit - 'spawning' to
                        ! different determinants.
                        
                        ! First check if the two electrons are in the same
                        ! spatial orbital. If they are then this will actually
                        ! 'spawn' to the same determinant, so just add the
                        ! contribution in.
                        if (is_in_pair(iorb, jorb)) then
                            ! Only take each contribution once.
                            do i = 1, nvecs
                                do j = i, nvecs
                                    spin_matrix(i,j) = spin_matrix(i,j) - &
                                        (parent_sign(kp_ind_1(i))*parent_sign(kp_ind_2(j)) + &
                                         parent_sign(kp_ind_2(i))*parent_sign(kp_ind_1(j)))/2.0_dp
                                end do
                            end do
                        else
                            ! If one of the new orbitals is already occupied
                            ! then there will be no contribution.
                            if (IsOcc(ilut_parent, ab_pair(iorb)) .or. &
                                IsOcc(ilut_parent, ab_pair(jorb))) cycle
                            ! Otherwise, we need to spawn to a different
                            ! determinant. Create it and find the parity.
                            call make_double(nI_parent, nI_child, iel, jel, ab_pair(iorb), ab_pair(jorb), ex, tParity)
                            call EncodeBitDet(nI_child, ilut_child)

                            if (tParity) then
                                child_sign = parent_sign
                            else
                                child_sign = -parent_sign
                            end if

                            call create_particle_kp_estimates(nI_child, ilut_child, child_sign, tNearlyFull)

                            if (tNearlyFull) then
                                call add_in_contribs(nvecs, krylov_array, krylov_ht, tFinished, tAllFinished, -1.0_dp, spin_matrix)
                                call clear_hash_table(spawn_ht)
                                tNearlyFull = .false.
                            end if

                        end if
                    end if

                end do ! Over all occupied electrons > iel (jel).
            end do ! Over all occupied electrons in this determiant (iel).

        end do ! Over all occupied determinants.

        tFinished = .true.
        do
            call add_in_contribs(nvecs, krylov_array, krylov_ht, tFinished, tAllFinished, -1.0_dp, spin_matrix)
            if (tAllFinished) exit
        end do

        ! Symmetrise the spin matrix.
        do i = 1, nvecs
            do j = 1, i-1
                spin_matrix(i,j) = spin_matrix(j,i)
            end do
        end do

        if (.not. tExcitedStateKP) then
            ! Now let SpawnedParts and SpawnedParts2 point back to their
            ! original arrays.
            SpawnedParts => SpawnVec
            SpawnedParts2 => SpawnVec2
        end if

        call halt_timer(subspace_spin_time)

    end subroutine calc_projected_spin

    function calc_amp_kp_hamil(parent_sign, prob, av_nspawn, HEl) result(child_sign)

        real(dp), intent(in) :: parent_sign(lenof_all_signs)
        real(dp), intent(in) :: prob
        real(dp), intent(in) :: av_nspawn
        real(dp), intent(in) :: HEl

        real(dp) :: child_sign(lenof_all_signs)
        real(dp) :: corrected_prob

        corrected_prob = prob*av_nspawn
        child_sign = HEl*parent_sign/corrected_prob

    end function calc_amp_kp_hamil

    subroutine create_particle_kp_estimates (nI_child, ilut_child, child_sign, tNearlyFull)

        use bit_rep_data, only: NOffSgn
        use load_balance_calcnodes, only: DetermineDetNode
        use hash, only: hash_table_lookup, add_hash_table_entry
        use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, SpawnedParts, spawn_ht
        use kp_fciqmc_data_mod, only: MaxSpawnedEachProc
        use SystemData, only: nel

        integer, intent(in) :: nI_child(nel)
        integer(n_int), intent(in) :: ilut_child(0:NIfTot)
        real(dp), intent(in) :: child_sign(lenof_all_signs)
        logical, intent(inout) :: tNearlyFull

        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: real_sign(lenof_all_signs)
        integer :: proc, ind, hash_val
        logical :: tSuccess

        call hash_table_lookup(nI_child, ilut_child, NIfDBO, spawn_ht, SpawnedParts, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! If this determinant is already in the spawned array.
            ! Extract the old sign.
            real_sign = transfer(SpawnedParts(NOffSgn:NOffSgn+lenof_all_signs-1, ind), real_sign)
            ! Find the total new sign.
            real_sign = real_sign + child_sign
            ! Encode the new sign.
            int_sign = transfer(real_sign, int_sign)
            SpawnedParts(NOffSgn:NOffSgn+lenof_all_signs-1, ind) = int_sign
        else
            ! If this determinant is a new entry to the spawned array.
            proc = DetermineDetNode(nel, nI_child, 0)

            SpawnedParts(0:NIfDBO, ValidSpawnedList(proc)) = ilut_child(0:NIfDBO)
            int_sign = transfer(child_sign, int_sign)
            SpawnedParts(NOffSgn:NOffSgn+lenof_all_signs-1, ValidSpawnedList(proc)) = int_sign
            call add_hash_table_entry(spawn_ht, ValidSpawnedList(proc), hash_val)

            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1

            if (ValidSpawnedList(proc)-InitialSpawnedSlots(proc) > MaxSpawnedEachProc) tNearlyFull = .true.
        end if

    end subroutine create_particle_kp_estimates

    subroutine distribute_spawns_kp_estimates(nspawns_this_proc)

        use AnnihilationMod, only: SendProcNewParts
        use FciMCData, only: SpawnedParts, SpawnedParts2

        integer, intent(out) :: nspawns_this_proc
        integer(n_int), pointer :: PointTemp(:,:)

        call SendProcNewParts(nspawns_this_proc, .false.)

        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

        nullify(PointTemp)

    end subroutine distribute_spawns_kp_estimates

    subroutine add_in_contribs(nvecs, krylov_array, krylov_ht, tFinished, tAllFinished, fac, est_matrix)

        use FciMCData, only: InitialSpawnedSlots, ValidSpawnedList, ll_node
        use Parallel_neci, only: MPIAllGather, nProcessors

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        type(ll_node), pointer, intent(in) :: krylov_ht(:)
        real(dp), intent(in) :: fac
        real(dp), intent(inout) :: est_matrix(:,:)

        logical, intent(in) :: tFinished
        logical, intent(out) :: tAllFinished
        logical :: tFinished_AllProcs(nProcessors)
        integer :: nspawns_this_proc, ierr

        call distribute_spawns_kp_estimates(nspawns_this_proc)

        call calc_contribs_spawn(nvecs, krylov_array, krylov_ht, nspawns_this_proc, fac, est_matrix)

        ValidSpawnedList = InitialSpawnedSlots

        call MPIAllGather(tFinished, 1, tFinished_AllProcs, 1, ierr)

        tAllFinished = all(tFinished_AllProcs)

    end subroutine add_in_contribs

    subroutine calc_contribs_spawn(nvecs, krylov_array, krylov_ht, nspawns_this_proc, fac, est_matrix)

        use bit_rep_data, only: NOffSgn
        use bit_reps, only: decode_bit_det
        use FciMCData, only: SpawnedParts, ll_node
        use hash, only: FindWalkerHash
        use SystemData, only: nel

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        type(ll_node), pointer, intent(in) :: krylov_ht(:)
        integer, intent(in) :: nspawns_this_proc
        real(dp), intent(in) :: fac
        real(dp), intent(inout) :: est_matrix(:,:)

        integer :: idet, det_hash, det_ind, i, j
        integer(n_int) :: ilut_spawn(0:NIfTot)
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: real_sign_1(lenof_all_signs), real_sign_2(lenof_all_signs)
        type(ll_node), pointer :: temp_node
        logical :: tDetFound

        ilut_spawn = 0_n_int

        do idet = 1, nspawns_this_proc

            ilut_spawn(0:NIfDBO) = SpawnedParts(0:NIfDBO, idet)
            int_sign = SpawnedParts(NOffSgn:NOffSgn+lenof_all_signs-1, idet)
            real_sign_1 = transfer(int_sign, real_sign_1)
            call decode_bit_det(nI_spawn, ilut_spawn)
            det_hash = FindWalkerHash(nI_spawn, size(krylov_ht))
            ! Point to the first node with this hash value in krylov_array.
            temp_node => krylov_ht(det_hash)

            if (temp_node%ind == 0) then
                ! If there are no determinants at all with this hash value in krylov_array.
                cycle
            else
                tDetFound = .false.
                do while (associated(temp_node))
                    if ( all(ilut_spawn(0:NIfDBO) == krylov_array(0:NIfDBO,temp_node%ind)) ) then
                        ! If this CurrentDets determinant has been found in krylov_array.
                        det_ind = temp_node%ind
                        tDetFound = .true.
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    temp_node => temp_node%next
                end do
                if (tDetFound) then
                    int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, det_ind)
                    real_sign_2 = transfer(int_sign, real_sign_1)
                    if (IsUnoccDet(real_sign_2)) cycle

                    ! Finally, add in the contribution to the projected Hamiltonian for each pair of Krylov vectors.
                    do i = 1, nvecs
                        do j = i, nvecs
                            est_matrix(i,j) = est_matrix(i,j) + &
                                fac*(real_sign_1(kp_ind_1(i))*real_sign_2(kp_ind_2(j)) + &
                                     real_sign_1(kp_ind_2(i))*real_sign_2(kp_ind_1(j)))/2.0_dp
                        end do
                    end do
                end if
            end if
        end do

    end subroutine calc_contribs_spawn

    subroutine calc_hamil_contribs_diag(nvecs, krylov_array, ndets, h_matrix, h_diag)
    
        use bit_rep_data, only: NOffSgn
        use FciMCData, only: determ_sizes, Hii
        use global_det_data, only: det_diagH
        use kp_fciqmc_data_mod, only: tSemiStochasticKPHamil
        use Parallel_neci, only: iProcIndex
        use SystemData, only: nel

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        integer, intent(in) :: ndets
        real(dp), intent(inout) :: h_matrix(:,:)
        real(dp), intent(in), optional :: h_diag(:)

        integer :: idet, i, j, min_idet
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: real_sign(lenof_all_signs)
        real(dp) :: h_diag_elem

        ! In semi-stochastic calculations the diagonal elements of the core space
        ! are taken care of in the core Hamiltonian calculation, so skip them here.
        ! Core determinants are always kept at the top of the list, so they're simple
        ! to skip.
        if (tSemiStochasticKPHamil) then
            min_idet = determ_sizes(iProcIndex) + 1
        else
            min_idet = 1
        end if

        do idet = min_idet, ndets
            int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, idet)
            real_sign = transfer(int_sign, real_sign)
            if (present(h_diag)) then
                h_diag_elem = h_diag(idet) + Hii
            else
                h_diag_elem = det_diagH(idet) + Hii
            end if

            ! Finally, add in the contribution to the projected Hamiltonian for each pair of Krylov vectors.
            do i = 1, nvecs
                do j = i, nvecs
                    h_matrix(i,j) = h_matrix(i,j) + &
                        h_diag_elem*(real_sign(kp_ind_1(i))*real_sign(kp_ind_2(j)) + &
                                     real_sign(kp_ind_2(i))*real_sign(kp_ind_1(j)))/2.0_dp
                end do
            end do

        end do

    end subroutine calc_hamil_contribs_diag

    subroutine calc_hamil_contribs_semistoch(nvecs, krylov_array, h_matrix, partial_vecs)

        use bit_rep_data, only: NOffSgn
        use FciMCData, only: determ_sizes
        use Parallel_neci, only: iProcIndex
        use SystemData, only: nel
    
        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        real(dp), intent(inout) :: h_matrix(:,:)
        real(dp), intent(in) :: partial_vecs(:,:)

        integer :: idet, i, j
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: real_sign(lenof_all_signs)

        do idet = 1, determ_sizes(iProcIndex)
            int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, idet)
            real_sign = transfer(int_sign, real_sign)

            do i = 1, nvecs
                do j = i, nvecs
                    h_matrix(i,j) = h_matrix(i,j) + &
                        (real_sign(kp_ind_2(j))*partial_vecs(kp_ind_1(i), idet) + &
                         real_sign(kp_ind_1(j))*partial_vecs(kp_ind_2(i), idet))/2.0_dp
                end do
            end do

        end do

    end subroutine calc_hamil_contribs_semistoch

end module kp_fciqmc_proj_est
