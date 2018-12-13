#include "macros.h"

module precond_annihilation_mod

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use AnnihilationMod, only: test_abort_spawn
    use AnnihilationMod, only: AnnihilateSpawnedParts, deterministic_annihilation
    use bit_rep_data
    use bit_reps, only: encode_sign, encode_part_sign, set_flag, &
                        get_initiator_flag, extract_spawn_hdiag
    use CalcData, only: OccupiedThresh, tTruncInitiator, PrecondSpawnCutoff
    use constants, only: n_int, lenof_sign, null_part, sizeof_int
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData
    use hash
    use Parallel_neci
    use load_balance, only: CalcHashTableStats
    use searching
    use sort_mod
    use SystemData, only: NEl, tHPHF

    implicit none

    logical, allocatable :: tSpawnedTo(:)

    contains

    subroutine precond_annihilation(TotWalkersNew, proj_energy, iter_data, tSingleProc)

        integer, intent(inout) :: TotWalkersNew
        real(dp), intent(inout) :: proj_energy(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        logical, intent(in) :: tSingleProc

        integer :: MaxIndex
        integer(n_int), pointer :: PointTemp(:,:)
        type(timer), save :: Compress_time
        integer :: ref_positions(lenof_sign)

        ! This routine will send all the newly-spawned particles to their
        ! correct processor. 
        call SendProcNewParts(MaxIndex, tSingleProc)

        ! CompressSpawnedList works on SpawnedParts arrays, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

        Compress_time%timer_name='Compression interface'
        call set_timer(Compress_time, 20)
        ! Now we want to order and compress the spawned list of particles. 
        ! This will also annihilate the newly spawned particles amongst themselves.
        ! MaxIndex will change to reflect the final number of unique determinants in the newly-spawned list, 
        ! and the particles will end up in the spawnedSign/SpawnedParts lists.
        call CompressSpawnedList(MaxIndex, iter_data)
        call halt_timer(Compress_time)

        !call set_timer(hash_test_time, 30)
        !call time_hash(MaxIndex)
        !call halt_timer(hash_test_time)

        call set_timer(proj_e_time, 30)
        call get_proj_energy(MaxIndex, proj_energy, ref_positions)
        call halt_timer(proj_e_time)

        call set_timer(precond_e_time, 30)
        call calc_e_and_set_init_flags(MaxIndex, proj_energy)
        call halt_timer(precond_e_time)

        call set_timer(rescale_time, 30)
        call rescale_spawns(MaxIndex, proj_energy)
        call halt_timer(rescale_time)

        call set_timer(precond_death_time, 30)
        call perform_death_with_precond(iter_data)
        call halt_timer(precond_death_time)

        ! If the semi-stochastic approach is being used then the following routine performs the
        ! annihilation of the deterministic states. These states are subsequently skipped in the
        ! AnnihilateSpawnedParts routine.
        if (tSemiStochastic) call deterministic_annihilation(iter_data)

        ! Binary search the main list and copy accross/annihilate determinants which are found.
        ! This will also remove the found determinants from the spawnedparts lists.
        call AnnihilateSpawnedParts(MaxIndex, TotWalkersNew, iter_data)

        call set_timer(Sort_Time, 30)
        call CalcHashTableStats(TotWalkersNew, iter_data)
        call halt_timer(Sort_Time)

    end subroutine precond_annihilation

    subroutine get_proj_energy(ValidSpawned, proj_energy, ref_positions)

        use CalcData, only: tau

        integer, intent(in) :: ValidSpawned
        real(dp), intent(out) :: proj_energy(lenof_sign)
        integer, intent(out) :: ref_positions(lenof_sign)

        integer :: i, run, ierr
        real(dp) :: SignTemp(lenof_sign), ref_pop(lenof_sign)
        logical :: ref_found(lenof_sign), tSuccess
        integer :: PartInd, DetHash

        proj_energy = 0.0_dp
        ref_positions = 0
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
                    ref_positions(run) = i
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

    end subroutine get_proj_energy

    subroutine calc_e_and_set_init_flags(ValidSpawned, proj_energy)

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

        use CalcData, only: tau, tEN2Init, tEN2Rigorous
        use global_det_data, only: det_diagH, replica_est_len
        use semi_stoch_procs, only: core_space_pos, check_determ_flag

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, j, PartInd, DetHash, determ_pos, ierr
        integer :: nI_spawn(nel)
        real(dp) :: spwnsign(lenof_sign), cursign(lenof_sign)
        real(dp) :: hdiag, mean_energy(replica_est_len)
        logical :: tCoreDet, tSuccess, abort(lenof_sign)

        real(dp) :: var_e_num(replica_est_len),         overlap(replica_est_len)
        real(dp) :: var_e_num_all(replica_est_len),     overlap_all(replica_est_len)
        real(dp) :: e_squared_num(replica_est_len),     e_squared_num_all(replica_est_len)
        real(dp) :: en2_pert(replica_est_len),          en2_pert_all(replica_est_len)
        real(dp) :: en2_new(replica_est_len),           en2_new_all(replica_est_len)
        real(dp) :: precond_e_num(replica_est_len),     precond_denom(replica_est_len)
        real(dp) :: precond_e_num_all(replica_est_len), precond_denom_all(replica_est_len)

        ! Allow room to send up to 1000 elements.
        real(dp) :: send_arr(6*replica_est_len)
        ! Allow room to receive up to 1000 elements.
        real(dp) :: recv_arr(6*replica_est_len)

        var_e_num = 0.0_dp
        overlap = 0.0_dp
        var_e_num_all = 0.0_dp
        overlap_all(1) = 0.0_dp

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

        tSpawnedTo = .false.

        ! Contributions from diagonal where necessary
        do i = 1, int(TotWalkers, sizeof_int)
            hdiag = det_diagH(i) + Hii
            call extract_sign(CurrentDets(:, i), cursign)
            do j = 1, replica_est_len
                overlap(j) = overlap(j) + cursign(2*j-1) * cursign(2*j)
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
                do j = 1, lenof_sign
                    if (abs(cursign(j)) > 1.e-12_dp) then
                        call set_flag(SpawnedParts(:,i), get_initiator_flag(j))
                    end if
                end do

                hdiag = det_diagH(PartInd) + Hii

                tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                if (tCoreDet) then
                    determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                    tSpawnedTo(determ_pos) = .true.
                    spwnsign = spwnsign + partial_determ_vecs(:, determ_pos)
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

                    ! Only get EN2 contribution is we're due to cancel this
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
                if (.not. tSpawnedTo(i)) then
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
        send_arr(1*replica_est_len+1:2*replica_est_len) = overlap
        send_arr(2*replica_est_len+1:3*replica_est_len) = en2_pert
        send_arr(3*replica_est_len+1:4*replica_est_len) = precond_e_num
        send_arr(4*replica_est_len+1:5*replica_est_len) = precond_denom
        send_arr(5*replica_est_len+1:6*replica_est_len) = e_squared_num

        call MPIBarrier(ierr)
        call MPISumAll(send_arr, recv_arr)

        var_e_num_all     = recv_arr(0*replica_est_len+1:1*replica_est_len)
        overlap_all       = recv_arr(1*replica_est_len+1:2*replica_est_len)
        en2_pert_all      = recv_arr(2*replica_est_len+1:3*replica_est_len)
        precond_e_num_all = recv_arr(3*replica_est_len+1:4*replica_est_len)
        precond_denom_all = recv_arr(4*replica_est_len+1:5*replica_est_len)
        e_squared_num_all = recv_arr(5*replica_est_len+1:6*replica_est_len)
        ! -------------------------------------------------------

        ! Use a slightly more rigorous expression for the variational plus
        ! EN2 energy.
        if (tEN2Rigorous) then
            en2_new = 0.0_dp
            en2_new_all = 0.0_dp
            mean_energy = var_e_num_all / overlap_all

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

                if (tSuccess) then
                    tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                    if (tCoreDet) then
                        determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                        tSpawnedTo(determ_pos) = .true.
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
                    if (.not. tSpawnedTo(i)) then

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
                write(replica_est_unit, '(1(3x,es20.13))', advance='no') overlap_all(j)
                write(replica_est_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(j), precond_denom_all(j)
            end do
            write(replica_est_unit,'()')
        end if

    end subroutine calc_e_and_set_init_flags

    subroutine rescale_spawns(ValidSpawned, proj_energy)

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i
        real(dp) :: spwnsign(lenof_sign), hdiag

        ! Find the weight spawned on the Hartree--Fock determinant.
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                partial_determ_vecs(:,i) = partial_determ_vecs(:,i) / &
                  (core_ham_diag(i) - proj_energy - proje_ref_energy_offsets)
            end do
        end if

        do i = 1, ValidSpawned
            hdiag = extract_spawn_hdiag(SpawnedParts(:,i))

            call extract_sign(SpawnedParts(:,i), spwnsign)
            spwnsign = spwnsign / (hdiag - proj_energy - proje_ref_energy_offsets - Hii)
            call encode_sign(SpawnedParts(:,i), spwnsign)
        end do

    end subroutine rescale_spawns

    subroutine perform_death_with_precond(iter_data)

        use DetBitOps, only: FindBitExcitLevel
        use fcimc_helper, only: walker_death
        use global_det_data, only: det_diagH

        type(fcimc_iter_data), intent(inout) :: iter_data

        integer :: ex_level, nI(nel), j
        real(dp) :: sgn(lenof_sign), hdiag

        do j = 1, int(TotWalkers, sizeof_int)

            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

            ex_level = FindBitExcitLevel(iLutRef(:,1), CurrentDets(:,j))
            hdiag = det_diagH(j)

            call decode_bit_det(nI, CurrentDets(:,j))

            call walker_death(iter_data, nI, CurrentDets(:,j), hdiag, &
                              sgn, j, ex_level)
        end do

    end subroutine perform_death_with_precond

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
            call open_create_stats('var_estimates', replica_est_unit)
            write(replica_est_unit, '("#", 4X, "Iteration")', advance='no')
            do j = 1, replica_est_len
                write(replica_est_unit, '(7x,"Energy numerator")', advance='no')
                write(replica_est_unit, '(5x,"Energy^2 numerator")', advance='no')
                if (tEN2Init) then
                    write(replica_est_unit, '(10x,"EN2 numerator")', advance='no')
                    write(replica_est_unit, '(9x,"Var + EN2 num.")', advance='no')
                end if
                if (tEN2Rigorous) then
                    write(replica_est_unit, '(8x,"Var + EN2 new")', advance='no')
                end if
                write(replica_est_unit, '(10x,"Normalisation")', advance='no')
                write(replica_est_unit, '(8x,"Precond. energy")', advance='no')
                write(replica_est_unit, '(9x,"Precond. norm.")', advance='no')
            end do
            write(replica_est_unit,'()')
        end if

    end subroutine open_replica_est_file

end module precond_annihilation_mod
