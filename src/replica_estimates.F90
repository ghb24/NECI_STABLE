#include "macros.h"

module replica_estimates

    use AnnihilationMod, only: test_abort_spawn
    use bit_rep_data
    use bit_reps, only: encode_sign, set_flag, &
                        get_initiator_flag, extract_spawn_hdiag
    use constants, only: lenof_sign, sizeof_int
    use FciMCData
    use hash
    use Parallel_neci
    use searching
    use sort_mod
    use SystemData, only: NEl

    implicit none

    contains

    subroutine get_ests_from_spawns(ValidSpawned, proj_energy)

        use CalcData, only: tPreCond, tReplicaEstimates, tTruncInitiator
        use CalcData, only: tSetInitFlagsBeforeDeath
        use fcimc_helper, only: set_init_flag_spawns_to_occ

        integer, intent(in) :: ValidSpawned
        real(dp), intent(out) :: proj_energy(lenof_sign)

        if (tPreCond .or. tReplicaEstimates) then
            ! The preconditioned energy is used in perturbative estimates
            ! (and also when performing preconditioned FCIQMC).
            call set_timer(proj_e_time, 30)
            call get_proj_e_for_preconditioner(ValidSpawned, proj_energy)
            call halt_timer(proj_e_time)
        end if

        if (tReplicaEstimates) then
            call set_timer(precond_e_time, 30)
            call calc_ests_and_set_init_flags(ValidSpawned, proj_energy)
            call halt_timer(precond_e_time)
        end if

        ! With preconditiong and a time step of 1, death will kill all
        ! walkers entirely, so the initiator criterion will not be applied
        ! unless we set flags now. Do this now, unless done already in
        ! calc_ests_and_set_init_flags (for efficiency improvement).
        if (tSetInitFlagsBeforeDeath .and. (.not. tReplicaEstimates)) then
            if (tTruncInitiator) then
                call set_timer(init_flag_time, 30)
                call set_init_flag_spawns_to_occ(ValidSpawned)
                call halt_timer(init_flag_time)
            end if
        end if

    end subroutine get_ests_from_spawns

    subroutine get_ests_from_spawns_simple(ValidSpawned, proj_energy)

        use CalcData, only: tPreCond, tReplicaEstimates, tTruncInitiator
        use CalcData, only: tSetInitFlagsBeforeDeath
        use fcimc_helper, only: set_init_flag_spawns_to_occ

        integer, intent(in) :: ValidSpawned
        real(dp), intent(out) :: proj_energy(lenof_sign)

        if (tPreCond .or. tReplicaEstimates) then
            ! The preconditioned energy is used in perturbative estimates
            ! (and also when performing preconditioned FCIQMC).
            call set_timer(proj_e_time, 30)
            call get_proj_e_for_preconditioner(ValidSpawned, proj_energy)
            call halt_timer(proj_e_time)
        end if

        if (tReplicaEstimates) then
            call set_timer(precond_e_time, 30)
            call calc_ests_simple_initiator(ValidSpawned, proj_energy)
            call halt_timer(precond_e_time)
        end if

    end subroutine get_ests_from_spawns_simple

    subroutine get_proj_e_for_preconditioner(ValidSpawned, proj_energy)

        use CalcData, only: tau

        integer, intent(in) :: ValidSpawned
        real(dp), intent(out) :: proj_energy(lenof_sign)

        integer :: i, run, ierr
        real(dp) :: SignTemp(lenof_sign), ref_pop(lenof_sign)
        logical :: ref_found(lenof_sign), tSuccess
        integer :: PartInd, DetHash

        proj_energy = 0.0_dp
        ref_found = .false.

        ! Find the weight spawned on the Hartree--Fock determinant.
        if (tSemiStochastic) then
            do run = 1, lenof_sign
                do i = 1, determ_sizes(iProcIndex)
                    if (DetBitEQ(core_space(0:NIfDBO, determ_displs(iProcIndex)+i), iLutRef(:,run), NIfDBO)) then
                        proj_energy(run) = -partial_determ_vecs(run,i)
                        ref_found(run) = .true.
                    end if
                end do
            end do
        end if

        do i = 1, ValidSpawned
            do run = 1, lenof_sign
                if (DetBitEQ(SpawnedParts(:,i), iLutRef(:,run), NIfDBO)) then
                    call extract_sign(SpawnedParts(:,i), SignTemp)
                    proj_energy(run) = proj_energy(run) - SignTemp(run)
                    ref_found(run) = .true.
                end if
            end do
        end do

        do run = 1, lenof_sign
            if (iProcIndex == iRefProc(run)) then
                if ( (.not. ref_found(run)) ) then
                    proj_energy(run) = -0.01_dp
                    write(6,'("WARNING: The reference determinant was not spawned to in the last iteration.")')
                else if (abs(proj_energy(run)) < 1.e-12_dp) then
                    proj_energy(run) = -0.01_dp
                    write(6,'("WARNING: The projected energy from the last iteration was zero. Setting to -0.1.")')
                end if

                call hash_table_lookup(ProjEDet(:,run), ilutRef(:,run), NIfDBO, HashIndex, &
                                       CurrentDets, PartInd, DetHash, tSuccess)

                if (tSuccess) then
                    call extract_sign(CurrentDets(:,PartInd), ref_pop)
                    proj_energy(run) = proj_energy(run)/ref_pop(run)
                else
                    write(6,'("WARNING: Reference determinant not found in main walker list.")')
                end if
            end if

            call MPIBarrier(ierr)
            call MPIBCast(proj_energy(run), 1, iRefProc(run))
        end do

        ! Remove time step
        proj_energy = proj_energy / tau

    end subroutine  get_proj_e_for_preconditioner

    subroutine calc_ests_and_set_init_flags(ValidSpawned, proj_energy)

        ! This routine calculates the various energy estimates to be printed
        ! to the file for the preconditioned approach.

        ! *IMPORTANT* It also sets the initiator flag for all spawnings to
        ! already-occupied determinants. This is because if tau=1, then
        ! all walkers will die in the death step to follow, which causes
        ! the initiator criteria to not be properly applied (as the
        ! AnnihilateSpawnedParts routine thinks these are spawnings to
        ! unoccupied determinants). So setting the initiator flag here
        ! allows everything to proceed as normal.

        ! This second point means that this routine *does* affect the
        ! walker dynamics - it doesn't just calculate some energies
        ! without affecting walkers. This step is done here for efficiency.

        use CalcData, only: tau, tEN2Init, tEN2Rigorous, tTruncInitiator
        use global_det_data, only: det_diagH, replica_est_len
        use semi_stoch_procs, only: core_space_pos, check_determ_flag

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, j, PartInd, DetHash, determ_pos, ierr
        integer :: nI_spawn(nel)
        real(dp) :: spwnsign(lenof_sign), cursign(lenof_sign)
        real(dp) :: hdiag, mean_energy(replica_est_len)
        logical :: tCoreDet, tSuccess, abort(lenof_sign)

        ! Allow room to send up to 1000 elements.
        real(dp) :: send_arr(6*replica_est_len)
        ! Allow room to receive up to 1000 elements.
        real(dp) :: recv_arr(6*replica_est_len)

        var_e_num = 0.0_dp
        rep_est_overlap = 0.0_dp
        var_e_num_all = 0.0_dp
        rep_est_overlap_all(1) = 0.0_dp

        e_squared_num = 0.0_dp
        e_squared_num_all = 0.0_dp

        en2_pert = 0.0_dp
        en2_pert_all = 0.0_dp

        precond_e_num = 0.0_dp
        precond_denom = 0.0_dp
        precond_e_num_all = 0.0_dp
        precond_denom_all = 0.0_dp

        do i = 1, replica_est_len
            mean_energy(i) = Hii + ( proj_energy(2*i-1) + proje_ref_energy_offsets(2*i-1) + &
                                     proj_energy(2*i) + proje_ref_energy_offsets(2*i) ) / 2.0_dp
        end do

        tCoreDet = .false.
        tSuccess = .false.
        abort = .false.

        tDetermSpawnedTo = .false.

        ! Contributions from diagonal where necessary
        do i = 1, int(TotWalkers, sizeof_int)
            hdiag = det_diagH(i) + Hii
            call extract_sign(CurrentDets(:, i), cursign)
            do j = 1, replica_est_len
                rep_est_overlap(j) = rep_est_overlap(j) + cursign(2*j-1) * cursign(2*j)
                var_e_num(j) = var_e_num(j) + hdiag * cursign(2*j-1) * cursign(2*j)
                e_squared_num(j) = e_squared_num(j) + (hdiag**2) * cursign(2*j-1) * cursign(2*j)
            end do
        end do

        do i = 1, ValidSpawned
            call decode_bit_det(nI_spawn, SpawnedParts(:,i))
            call extract_sign(SpawnedParts(:,i), spwnsign)

            ! Now add in the diagonal elements
            call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            if (tSuccess) then
                call extract_sign(CurrentDets(:,PartInd), cursign)

                ! Set initiator flags for the spawning, before the currently
                ! occupied determinant is potentially killed in the death step.
                if (tTruncInitiator) then
                    do j = 1, lenof_sign
                        if (abs(cursign(j)) > 1.e-12_dp) then
                            call set_flag(SpawnedParts(:,i), get_initiator_flag(j))
                        end if
                    end do
                end if

                hdiag = det_diagH(PartInd) + Hii

                if (tSemiStochastic) then
                    tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                    if (tCoreDet) then
                        determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                        tDetermSpawnedTo(determ_pos) = .true.
                        spwnsign = spwnsign + partial_determ_vecs(:, determ_pos)
                    end if
                end if

                do j = 1, replica_est_len
                    var_e_num(j) = var_e_num(j) - &
                     (spwnsign(2*j-1)*cursign(2*j) + spwnsign(2*j)*cursign(2*j-1)) / (2.0_dp*tau)

                    precond_e_num(j) = precond_e_num(j) - &
                     (spwnsign(2*j-1)*hdiag*cursign(2*j) + spwnsign(2*j)*hdiag*cursign(2*j-1)) / &
                      (2.0_dp*(mean_energy(j) - hdiag))

                    precond_denom(j) = precond_denom(j) - &
                     (spwnsign(2*j-1)*cursign(2*j) + spwnsign(2*j) * cursign(2*j-1)) / (2.0_dp*(mean_energy(j) - hdiag))

                    e_squared_num(j) = e_squared_num(j) - &
                     (spwnsign(2*j-1)*hdiag*cursign(2*j) + spwnsign(2*j)*hdiag*cursign(2*j-1)) / tau
                end do
            end if

            ! Add in the contributions corresponding to off-diagonal
            ! elements of the Hamiltonian
            do j = 1, replica_est_len
                if (abs(spwnsign(2*j-1)) > 1.e-12_dp .and. abs(spwnsign(2*j)) > 1.e-12_dp) then
                    hdiag = extract_spawn_hdiag(SpawnedParts(:,i))

                    precond_e_num(j) = precond_e_num(j) + &
                      spwnsign(2*j-1) * spwnsign(2*j) / (tau*(mean_energy(j) - hdiag))

                    e_squared_num(j) = e_squared_num(j) + spwnsign(2*j-1) * spwnsign(2*j) / (tau**2)

                    ! Only get EN2 contribution if we're due to cancel this
                    ! spawning on both replicas
                    if (tEN2Init .and. (.not. tSuccess)) then
                        abort(2*j-1) = test_abort_spawn(SpawnedParts(:, i), 2*j-1)
                        abort(2*j)   = test_abort_spawn(SpawnedParts(:, i), 2*j)

                        if (abort(2*j-1) .and. abort(2*j)) then
                            en2_pert(j) = en2_pert(j) + &
                              spwnsign(2*j-1)*spwnsign(2*j) / ( (tau**2) * (mean_energy(j) - hdiag ) )
                        end if
                    end if

                end if
            end do

        end do

        ! Contribution from deterministic states that were not spawned to
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                if (.not. tDetermSpawnedTo(i)) then
                    call extract_sign(CurrentDets(:, indices_of_determ_states(i)), cursign)

                    do j = 1, replica_est_len
                        ! Variational energy terms
                        var_e_num(j) = var_e_num(j) - &
                            (partial_determ_vecs(2*j-1, i) * cursign(2*j) + &
                             partial_determ_vecs(2*j, i) * cursign(2*j-1)) / (2.0_dp * tau)

                        ! Preconditioned estimator
                        precond_e_num(j) = precond_e_num(j) + &
                          partial_determ_vecs(2*j-1, i) * partial_determ_vecs(2*j, i) / &
                            (tau * (mean_energy(j) - core_ham_diag(i) - Hii))

                        precond_e_num(j) = precond_e_num(j) - &
                          (partial_determ_vecs(2*j-1, i) * (core_ham_diag(i) + Hii) * cursign(2*j) + &
                           partial_determ_vecs(2*j, i) * (core_ham_diag(i) + Hii) * cursign(2*j-1)) / &
                            (2.0_dp*(mean_energy(j) - core_ham_diag(i) - Hii))

                        precond_denom(j) = precond_denom(j) - &
                          (partial_determ_vecs(2*j-1, i) * cursign(2*j) + partial_determ_vecs(2*j, i) * cursign(2*j-1)) / &
                            (2.0_dp*(mean_energy(j) - core_ham_diag(i) - Hii))

                        ! E squared terms
                        e_squared_num(j) = e_squared_num(j) + &
                          partial_determ_vecs(2*j-1, i) * partial_determ_vecs(2*j, i) / (tau**2)

                        e_squared_num(j) = e_squared_num(j) - &
                         (partial_determ_vecs(2*j-1, i) * (core_ham_diag(i) + Hii) * cursign(2*j) + &
                          partial_determ_vecs(2*j, i) * (core_ham_diag(i) + Hii) * cursign(2*j-1)) / tau
                    end do

                end if
            end do
        end if

        ! ---- MPI communication --------------------------------
        send_arr(0*replica_est_len+1:1*replica_est_len) = var_e_num
        send_arr(1*replica_est_len+1:2*replica_est_len) = rep_est_overlap
        send_arr(2*replica_est_len+1:3*replica_est_len) = en2_pert
        send_arr(3*replica_est_len+1:4*replica_est_len) = precond_e_num
        send_arr(4*replica_est_len+1:5*replica_est_len) = precond_denom
        send_arr(5*replica_est_len+1:6*replica_est_len) = e_squared_num

        call MPIBarrier(ierr)
        call MPISumAll(send_arr, recv_arr)

        var_e_num_all       = recv_arr(0*replica_est_len+1:1*replica_est_len)
        rep_est_overlap_all = recv_arr(1*replica_est_len+1:2*replica_est_len)
        en2_pert_all        = recv_arr(2*replica_est_len+1:3*replica_est_len)
        precond_e_num_all   = recv_arr(3*replica_est_len+1:4*replica_est_len)
        precond_denom_all   = recv_arr(4*replica_est_len+1:5*replica_est_len)
        e_squared_num_all   = recv_arr(5*replica_est_len+1:6*replica_est_len)
        ! -------------------------------------------------------

        ! Use a slightly more rigorous expression for the variational plus
        ! EN2 energy.
        if (tEN2Rigorous) then
            en2_new = 0.0_dp
            en2_new_all = 0.0_dp
            mean_energy = var_e_num_all / rep_est_overlap_all

            do i = 1, int(TotWalkers, sizeof_int)
                hdiag = det_diagH(i) + Hii
                call extract_sign(CurrentDets(:, i), cursign)
                do j = 1, replica_est_len
                    en2_new(j) = en2_new(j) + hdiag * cursign(2*j-1) * cursign(2*j)
                end do
            end do

            do i = 1, ValidSpawned
                call decode_bit_det(nI_spawn, SpawnedParts(:,i))
                call extract_sign(SpawnedParts(:,i), spwnsign)

                ! Now add in the diagonal elements
                call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                       CurrentDets, PartInd, DetHash, tSuccess)

                if (tSuccess .and. tSemiStochastic) then
                    tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                    if (tCoreDet) then
                        determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                        tDetermSpawnedTo(determ_pos) = .true.
                        spwnsign = spwnsign + partial_determ_vecs(:, determ_pos)
                    end if
                end if

                do j = 1, replica_est_len
                    if (abs(spwnsign(2*j-1)) > 1.e-12_dp .and. abs(spwnsign(2*j)) > 1.e-12_dp) then
                        hdiag = extract_spawn_hdiag(SpawnedParts(:,i))

                        en2_new(j) = en2_new(j) + &
                          spwnsign(2*j-1) * spwnsign(2*j) / ((tau**2)*(mean_energy(j) - hdiag))
                    end if
                end do
            end do

            if (tSemiStochastic) then
                do i = 1, determ_sizes(iProcIndex)
                    if (.not. tDetermSpawnedTo(i)) then

                        do j = 1, replica_est_len
                            en2_new(j) = en2_new(j) + &
                              partial_determ_vecs(2*j-1, i) * partial_determ_vecs(2*j, i) / &
                                ((tau)**2 * (mean_energy(j) - core_ham_diag(i) - Hii))
                        end do

                    end if
                end do
            end if

            call MPISumAll(en2_new, en2_new_all)
        end if

        if (iProcIndex == 0) then
            write(replica_est_unit, '(1x,i13)', advance='no') Iter + PreviousCycles

            do j = 1, replica_est_len
                write(replica_est_unit, '(2(3x,es20.13))', advance='no') var_e_num_all(j), e_squared_num_all(j)
                if (tEN2Init) then
                    write(replica_est_unit, '(2(3x,es20.13))', advance='no') en2_pert_all(j), var_e_num_all(j)+en2_pert_all(j)
                end if
                if (tEN2Rigorous) then
                    write(replica_est_unit, '(1(3x,es20.13))', advance='no') en2_new_all(j)
                end if
                write(replica_est_unit, '(1(3x,es20.13))', advance='no') rep_est_overlap_all(j)
                write(replica_est_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(j), precond_denom_all(j)
            end do
            write(replica_est_unit,'()')
        end if

    end subroutine calc_ests_and_set_init_flags

    subroutine calc_ests_simple_initiator(ValidSpawned, proj_energy)

        ! This routine calculates the various energy estimates to be printed
        ! to the file for the preconditioned approach.

        use CalcData, only: tau, tEN2Init, tEN2Rigorous, tTruncInitiator
        use global_det_data, only: det_diagH, replica_est_len
        use semi_stoch_procs, only: core_space_pos, check_determ_flag

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, j, PartInd, DetHash, determ_pos, nrepeats, ierr
        integer :: nI_spawn(nel)
        real(dp) :: spwnsign_init(lenof_sign), spwnsign_non(lenof_sign)
        real(dp) :: spwnsign(lenof_sign), cursign(lenof_sign)
        real(dp) :: hdiag, mean_energy(replica_est_len)
        logical :: tCoreDet, tSuccess, abort(lenof_sign)
        character(len=*), parameter :: t_r = 'calc_ests_simple_initiator'

        real(dp) :: b_term(replica_est_len), b_term_all(replica_est_len)
        real(dp) :: c_term(replica_est_len), c_term_all(replica_est_len)

        ! Allow room to send up to 1000 elements.
        real(dp) :: send_arr(6*replica_est_len)
        ! Allow room to receive up to 1000 elements.
        real(dp) :: recv_arr(6*replica_est_len)

        var_e_num = 0.0_dp
        rep_est_overlap = 0.0_dp
        var_e_num_all = 0.0_dp
        rep_est_overlap_all(1) = 0.0_dp

        en2_pert = 0.0_dp
        en2_pert_all = 0.0_dp

        e_squared_num = 0.0_dp
        e_squared_num_all = 0.0_dp

        precond_e_num = 0.0_dp
        precond_denom = 0.0_dp
        precond_e_num_all = 0.0_dp
        precond_denom_all = 0.0_dp

        do i = 1, replica_est_len
            mean_energy(i) = Hii + ( proj_energy(2*i-1) + proje_ref_energy_offsets(2*i-1) + &
                                     proj_energy(2*i) + proje_ref_energy_offsets(2*i) ) / 2.0_dp
        end do

        tCoreDet = .false.
        tSuccess = .false.
        abort = .false.

        tDetermSpawnedTo = .false.

        ! Contributions from diagonal where necessary
        do i = 1, int(TotWalkers, sizeof_int)
            hdiag = det_diagH(i) + Hii
            call extract_sign(CurrentDets(:, i), cursign)
            do j = 1, replica_est_len
                rep_est_overlap(j) = rep_est_overlap(j) + cursign(2*j-1) * cursign(2*j)
                var_e_num(j) = var_e_num(j) + hdiag * cursign(2*j-1) * cursign(2*j)
                e_squared_num(j) = e_squared_num(j) + (hdiag**2) * cursign(2*j-1) * cursign(2*j)
            end do
        end do

        i = 1
        do
            call decode_bit_det(nI_spawn, SpawnedParts(:,i))

            ! How many times is this particular determinant repeated in the
            ! spawning array (either 0 or 1)?
            nrepeats = 0
            spwnsign_init = 0.0_dp
            spwnsign_non = 0.0_dp

            ! Only need to check initiator status for one replica - if
            ! one is an initiator then they all are, when using the simple
            ! initiator approximation.
            if (test_flag(SpawnedParts(:,i), get_initiator_flag(1))) then
                ! This is an initiator
                call extract_sign(SpawnedParts(:,i), spwnsign_init)
                if (i+1 <= ValidSpawned) then
                    if (DetBitEq(SpawnedParts(:,i), SpawnedParts(:,i+1), NIfDBO)) then
                        ! This next state is a different state, and so will
                        ! be a non-initiator
                        call extract_sign(SpawnedParts(:,i+1), spwnsign_non)
                        nrepeats = 1
                    end if
                end if
            else
                call extract_sign(SpawnedParts(:,i), spwnsign_non)
            end if

            spwnsign = spwnsign_init + spwnsign_non

            ! Now add in the diagonal elements
            call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            if (tSuccess) then
                call extract_sign(CurrentDets(:,PartInd), cursign)

                hdiag = det_diagH(PartInd) + Hii

                if (tSemiStochastic) then
                    tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                    if (tCoreDet) then
                        determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                        tDetermSpawnedTo(determ_pos) = .true.
                        spwnsign = spwnsign + partial_determ_vecs(:, determ_pos)
                    end if
                end if

                do j = 1, replica_est_len
                    var_e_num(j) = var_e_num(j) - &
                     (spwnsign(2*j-1)*cursign(2*j) + spwnsign(2*j)*cursign(2*j-1)) / (2.0_dp*tau)

                    precond_e_num(j) = precond_e_num(j) - &
                     (spwnsign(2*j-1)*hdiag*cursign(2*j) + spwnsign(2*j)*hdiag*cursign(2*j-1)) / &
                      (2.0_dp*(mean_energy(j) - hdiag))

                    precond_denom(j) = precond_denom(j) - &
                     (spwnsign(2*j-1)*cursign(2*j) + spwnsign(2*j) * cursign(2*j-1)) / (2.0_dp*(mean_energy(j) - hdiag))

                    e_squared_num(j) = e_squared_num(j) - &
                     (spwnsign(2*j-1)*hdiag*cursign(2*j) + spwnsign(2*j)*hdiag*cursign(2*j-1)) / tau

                    !if (tEN2Init) then
                    !    en2_pert(j) = en2_pert(j) - &
                    !      (spwnsign_non(2*j-1)*cursign(2*j) + spwnsign_non(2*j)*cursign(2*j-1)) / (2.0_dp*tau)
                    !end if
                end do
            end if

            ! Add in the contributions corresponding to off-diagonal
            ! elements of the Hamiltonian
            do j = 1, replica_est_len
                hdiag = extract_spawn_hdiag(SpawnedParts(:,i))

                precond_e_num(j) = precond_e_num(j) + &
                  spwnsign(2*j-1) * spwnsign(2*j) / (tau*(mean_energy(j) - hdiag))

                e_squared_num(j) = e_squared_num(j) + spwnsign(2*j-1) * spwnsign(2*j) / (tau**2)

                ! Only get EN2 contribution if we're due to cancel this
                ! spawning on both replicas
                if (tEN2Init) then
                    en2_pert(j) = en2_pert(j) + &
                      spwnsign_non(2*j-1)*spwnsign_non(2*j) / ( (tau**2) * (mean_energy(j) - hdiag ) )
                end if
            end do

            i = i + 1 + nrepeats
            if (i > ValidSpawned + 1) call stop_all(t_r, "The spawning array index is larger than it should be.")
            if (i == ValidSpawned + 1) exit
        end do

        ! Contribution from deterministic states that were not spawned to
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                if (.not. tDetermSpawnedTo(i)) then
                    call extract_sign(CurrentDets(:, indices_of_determ_states(i)), cursign)

                    do j = 1, replica_est_len
                        ! Variational energy terms
                        var_e_num(j) = var_e_num(j) - &
                            (partial_determ_vecs(2*j-1, i) * cursign(2*j) + &
                             partial_determ_vecs(2*j, i) * cursign(2*j-1)) / (2.0_dp * tau)

                        ! Preconditioned estimator
                        precond_e_num(j) = precond_e_num(j) + &
                          partial_determ_vecs(2*j-1, i) * partial_determ_vecs(2*j, i) / &
                            (tau * (mean_energy(j) - core_ham_diag(i) - Hii))

                        precond_e_num(j) = precond_e_num(j) - &
                          (partial_determ_vecs(2*j-1, i) * (core_ham_diag(i) + Hii) * cursign(2*j) + &
                           partial_determ_vecs(2*j, i) * (core_ham_diag(i) + Hii) * cursign(2*j-1)) / &
                            (2.0_dp*(mean_energy(j) - core_ham_diag(i) - Hii))

                        precond_denom(j) = precond_denom(j) - &
                          (partial_determ_vecs(2*j-1, i) * cursign(2*j) + partial_determ_vecs(2*j, i) * cursign(2*j-1)) / &
                            (2.0_dp*(mean_energy(j) - core_ham_diag(i) - Hii))

                        ! E squared terms
                        e_squared_num(j) = e_squared_num(j) + &
                          partial_determ_vecs(2*j-1, i) * partial_determ_vecs(2*j, i) / (tau**2)

                        e_squared_num(j) = e_squared_num(j) - &
                         (partial_determ_vecs(2*j-1, i) * (core_ham_diag(i) + Hii) * cursign(2*j) + &
                          partial_determ_vecs(2*j, i) * (core_ham_diag(i) + Hii) * cursign(2*j-1)) / tau
                    end do

                end if
            end do
        end if

        ! ---- MPI communication --------------------------------
        send_arr(0*replica_est_len+1:1*replica_est_len) = var_e_num
        send_arr(1*replica_est_len+1:2*replica_est_len) = rep_est_overlap
        send_arr(2*replica_est_len+1:3*replica_est_len) = en2_pert
        send_arr(3*replica_est_len+1:4*replica_est_len) = precond_e_num
        send_arr(4*replica_est_len+1:5*replica_est_len) = precond_denom
        send_arr(5*replica_est_len+1:6*replica_est_len) = e_squared_num

        call MPIBarrier(ierr)
        call MPISumAll(send_arr, recv_arr)

        var_e_num_all       = recv_arr(0*replica_est_len+1:1*replica_est_len)
        rep_est_overlap_all = recv_arr(1*replica_est_len+1:2*replica_est_len)
        en2_pert_all        = recv_arr(2*replica_est_len+1:3*replica_est_len)
        precond_e_num_all   = recv_arr(3*replica_est_len+1:4*replica_est_len)
        precond_denom_all   = recv_arr(4*replica_est_len+1:5*replica_est_len)
        e_squared_num_all   = recv_arr(5*replica_est_len+1:6*replica_est_len)
        ! -------------------------------------------------------

        ! Use a slightly more rigorous expression for the variational plus
        ! EN2 energy.
        if (tEN2Rigorous) then
            en2_new = 0.0_dp
            en2_new_all = 0.0_dp
            b_term = 0.0_dp
            b_term_all = 0.0_dp
            c_term = 0.0_dp
            c_term_all = 0.0_dp

            mean_energy = var_e_num_all / rep_est_overlap_all

            do i = 1, int(TotWalkers, sizeof_int)
                hdiag = det_diagH(i) + Hii
                call extract_sign(CurrentDets(:, i), cursign)
                do j = 1, replica_est_len
                    en2_new(j) = en2_new(j) + hdiag * cursign(2*j-1) * cursign(2*j)
                    b_term(j)  = b_term(j)  + cursign(2*j-1) * cursign(2*j)
                    c_term(j)  = c_term(j)  - cursign(2*j-1) * cursign(2*j) / (mean_energy(j) - hdiag)
                end do
            end do

            i = 1
            do
                call decode_bit_det(nI_spawn, SpawnedParts(:,i))

                ! How many times is this particular determinant repeated in the
                ! spawning array (either 0 or 1)?
                nrepeats = 0
                spwnsign_init = 0.0_dp
                spwnsign_non = 0.0_dp

                ! Only need to check initiator status for one replica - if
                ! one is an initiator then they all are, when using the simple
                ! initiator approximation.
                if (test_flag(SpawnedParts(:,i), get_initiator_flag(1))) then
                    ! This is an initiator
                    call extract_sign(SpawnedParts(:,i), spwnsign_init)
                    if (i+1 <= ValidSpawned) then
                        if (DetBitEq(SpawnedParts(:,i), SpawnedParts(:,i+1), NIfDBO)) then
                            ! This next state is a different state, and so will
                            ! be a non-initiator
                            call extract_sign(SpawnedParts(:,i+1), spwnsign_non)
                            nrepeats = 1
                        end if
                    end if
                else
                    call extract_sign(SpawnedParts(:,i), spwnsign_non)
                end if

                spwnsign = spwnsign_init + spwnsign_non

                ! Now add in the diagonal elements
                call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                       CurrentDets, PartInd, DetHash, tSuccess)

                if (tSuccess) then
                    call extract_sign(CurrentDets(:,PartInd), cursign)

                    hdiag = det_diagH(PartInd) + Hii

                    if (tSemiStochastic) then
                        tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                        if (tCoreDet) then
                            determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                            tDetermSpawnedTo(determ_pos) = .true.
                            spwnsign = spwnsign + partial_determ_vecs(:, determ_pos)
                        end if
                    end if

                    do j = 1, replica_est_len
                        b_term(j) = b_term(j) + &
                         (spwnsign(2*j-1)*cursign(2*j) + spwnsign(2*j)*cursign(2*j-1)) / (2.0_dp*tau*(mean_energy(j) - hdiag))
                    end do
                end if

                do j = 1, replica_est_len
                    hdiag = extract_spawn_hdiag(SpawnedParts(:,i))

                    en2_new(j) = en2_new(j) + &
                      spwnsign(2*j-1) * spwnsign(2*j) / ((tau**2)*(mean_energy(j) - hdiag))
                end do

                i = i + 1 + nrepeats
                if (i > ValidSpawned + 1) call stop_all(t_r, "The spawning array index is larger than it should be.")
                if (i == ValidSpawned + 1) exit
            end do

            if (tSemiStochastic) then
                do i = 1, determ_sizes(iProcIndex)
                    if (.not. tDetermSpawnedTo(i)) then
                        call extract_sign(CurrentDets(:, indices_of_determ_states(i)), cursign)

                        do j = 1, replica_est_len
                            en2_new(j) = en2_new(j) + &
                              partial_determ_vecs(2*j-1, i) * partial_determ_vecs(2*j, i) / &
                                ((tau)**2 * (mean_energy(j) - core_ham_diag(i) - Hii))

                            b_term(j) = b_term(j) + &
                              (partial_determ_vecs(2*j-1, i) * cursign(2*j) + partial_determ_vecs(2*j, i) * cursign(2*j-1)) / &
                                (2.0_dp*tau*(mean_energy(j) - core_ham_diag(i) - Hii))
                        end do

                    end if
                end do
            end if

            call MPISumAll(en2_new, en2_new_all)
            call MPISumAll(b_term,  b_term_all)
            call MPISumAll(c_term,  c_term_all)
        end if

        if (iProcIndex == 0) then
            write(replica_est_unit, '(1x,i13)', advance='no') Iter + PreviousCycles

            do j = 1, replica_est_len
                write(replica_est_unit, '(2(3x,es20.13))', advance='no') var_e_num_all(j), e_squared_num_all(j)
                if (tEN2Init) then
                    write(replica_est_unit, '(2(3x,es20.13))', advance='no') en2_pert_all(j), var_e_num_all(j)+en2_pert_all(j)
                end if
                if (tEN2Rigorous) then
                    write(replica_est_unit, '(1(3x,es20.13))', advance='no') en2_new_all(j)
                end if
                write(replica_est_unit, '(1(3x,es20.13))', advance='no') rep_est_overlap_all(j)
                if (tEN2Rigorous) then
                    write(replica_est_unit, '(2(3x,es20.13))', advance='no') b_term_all(j), c_term_all(j)
                end if
                write(replica_est_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(j), precond_denom_all(j)
            end do
            write(replica_est_unit,'()')
        end if

    end subroutine calc_ests_simple_initiator

    subroutine time_hash(ValidSpawned)

        integer, intent(in) :: ValidSpawned

        integer :: i, PartInd, DetHash, nI(nel)
        logical :: tSuccess

        do i = 1, ValidSpawned
            call decode_bit_det(nI, SpawnedParts(:,i))
            call hash_table_lookup(nI, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)
        end do

    end subroutine time_hash

    subroutine open_replica_est_file()

        use CalcData, only: tEN2Init, tEN2Rigorous
        use fcimc_output, only: open_create_stats
        use global_det_data, only: replica_est_len
        use util_mod, only: get_free_unit

        integer :: j

        if (iProcIndex == 0) then
            replica_est_unit = get_free_unit()
            call open_create_stats('replica_estimates', replica_est_unit)
            write(replica_est_unit, '("#", 4X, "Iteration")', advance='no')
            do j = 1, replica_est_len
                write(replica_est_unit, '(7x,"Energy numerator")', advance='no')
                write(replica_est_unit, '(5x,"Energy^2 numerator")', advance='no')
                if (tEN2Init) then
                    write(replica_est_unit, '(10x,"EN2 numerator")', advance='no')
                    write(replica_est_unit, '(9x,"Var + EN2 num.")', advance='no')
                end if
                if (tEN2Rigorous) then
                    write(replica_est_unit, '(10x,"Var + EN2 new")', advance='no')
                end if
                write(replica_est_unit, '(10x,"Normalisation")', advance='no')
                if (tEN2Rigorous) then
                    write(replica_est_unit, '(11x,"B correction")', advance='no')
                    write(replica_est_unit, '(10x,"C denominator")', advance='no')
                end if
                write(replica_est_unit, '(8x,"Precond. energy")', advance='no')
                write(replica_est_unit, '(9x,"Precond. norm.")', advance='no')
            end do
            write(replica_est_unit,'()')
        end if

    end subroutine open_replica_est_file

end module replica_estimates
