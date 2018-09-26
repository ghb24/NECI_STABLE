#include "macros.h"

module precond_annihilation_mod

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use AnnihilationMod, only: test_abort_spawn
    use AnnihilationMod, only: AnnihilateSpawnedParts, deterministic_annihilation
    use bit_rep_data
    use bit_reps, only: encode_sign, encode_part_sign, set_flag, get_initiator_flag
    use CalcData, only: OccupiedThresh, tTruncInitiator, PrecondSpawnCutoff
    use constants, only: n_int, lenof_sign, null_part, sizeof_int
    use determinants, only: get_helement
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData
    use hash
    use hphf_integrals, only: hphf_diag_helement
    use Parallel_neci
    use load_balance, only: CalcHashTableStats
    use searching
    use sort_mod
    use SystemData, only: NEl, tHPHF

    implicit none

    logical, allocatable :: tSpawnedTo(:)

    contains

    subroutine precond_annihilation(TotWalkersNew, iter_data, tSingleProc)

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        logical, intent(in) :: tSingleProc

        integer :: MaxIndex
        integer(n_int), pointer :: PointTemp(:,:)
        type(timer), save :: Compress_time
        real(dp) :: proj_energy(lenof_sign), var_e_num_all, overlap_all

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

        !call set_timer(hii_test_time, 30)
        !call time_hii(MaxIndex)
        !call halt_timer(hii_test_time)

        call set_timer(var_e_time, 30)
        call calc_var_energy(MaxIndex, var_e_num_all, overlap_all)
        call halt_timer(var_e_time)

        call set_timer(proj_e_time, 30)
        call get_proj_energy(MaxIndex, proj_energy)
        call halt_timer(proj_e_time)

        call set_timer(precond_e_time, 30)
        call calc_precond_energy(MaxIndex, proj_energy)
        call halt_timer(precond_e_time)

        call set_timer(precond_round_time, 30)
        call round_spawns(MaxIndex, precondSpawnCutoff)
        call halt_timer(precond_round_time)

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

    subroutine calc_var_energy(ValidSpawned, var_e_num_all, overlap_all)

        use CalcData, only: tau, tEN2Init
        use global_det_data, only: det_diagH

        integer, intent(in) :: ValidSpawned
        real(dp), intent(out) :: var_e_num_all, overlap_all

        integer :: i, j, PartInd, DetHash, ierr
        integer :: nI_spawn(nel)
        real(dp) :: SpawnedSign(lenof_sign), CurrentSign(lenof_sign)
        real(dp) :: h_diag, var_e_num, overlap, en2_pert, en2_pert_all
        logical :: tSuccess, abort(lenof_sign), pert_contrib

        var_e_num = 0.0_dp
        overlap = 0.0_dp
        var_e_num_all = 0.0_dp
        overlap_all = 0.0_dp
        en2_pert = 0.0_dp
        en2_pert_all = 0.0_dp

        tSuccess = .false.
        abort = .false.

        ! Contribution from diagonal
        do i = 1, int(TotWalkers, sizeof_int)
            h_diag = det_diagH(i) + Hii
            call extract_sign(CurrentDets(:, i), CurrentSign)
            overlap = overlap + CurrentSign(1) * CurrentSign(2)
            var_e_num = var_e_num + h_diag * CurrentSign(1) * CurrentSign(2)
        end do

        ! Contribution from deterministic space
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                call extract_sign(CurrentDets(:, indices_of_determ_states(i)), CurrentSign)
                var_e_num = var_e_num - &
                    (partial_determ_vecs(1, i) * CurrentSign(2) + &
                     partial_determ_vecs(2, i) * CurrentSign(1)) / (2.0_dp * tau)
            end do
        end if

        ! Contribution from stochastic spawnings
        do i = 1, ValidSpawned
            ! Check if this spawned determinant is already in the walker list
            call decode_bit_det(nI_spawn, SpawnedParts(:,i))
            call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            if (tSuccess) then
                call extract_sign(SpawnedParts(:,i), SpawnedSign)
                call extract_sign(CurrentDets(:,PartInd), CurrentSign)
                var_e_num = var_e_num - &
                    (SpawnedSign(1) * CurrentSign(2) + SpawnedSign(2) * CurrentSign(1))/(2.0_dp * tau)
            end if
        end do

        call MPIBarrier(ierr)
        call MPISumAll(var_e_num, var_e_num_all)
        call MPISumAll(overlap, overlap_all)

        ! Contribution from stochastic spawnings
        if (tEN2Init) then
            do i = 1, ValidSpawned
                ! Check if this spawned determinant is already in the walker list
                call decode_bit_det(nI_spawn, SpawnedParts(:,i))
                call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                       CurrentDets, PartInd, DetHash, tSuccess)

                if (.not. tSuccess) then
                    ! Not already in the main list - are we going to abort it?
                    ! If so, add in a contribution to the EN2 correction
                    do j = 1, lenof_sign
                        abort(j) = test_abort_spawn(SpawnedParts(:, i), j)
                    end do
                    call extract_sign(SpawnedParts(:,i), SpawnedSign)

                    pert_contrib = abort(1) .and. abort(2) .and. &
                      abs(SpawnedSign(1)) > 1.e-12_dp .and. abs(SpawnedSign(2)) > 1.e-12_dp

                    if (pert_contrib) then
                        if (tHPHF) then
                            h_diag = hphf_diag_helement(nI_spawn, SpawnedParts(:, i))
                        else
                            h_diag = get_helement(nI_spawn, nI_spawn, 0)
                        end if

                        en2_pert = en2_pert + &
                          SpawnedSign(1)*SpawnedSign(2) / ( (tau**2) * ((var_e_num_all / overlap_all) - h_diag ) )
                    end if

                end if
            end do

            call MPIBarrier(ierr)
            call MPISumAll(en2_pert, en2_pert_all)
        end if

        if (iProcIndex == 0) then
            write(var_unit, '(1x,i13)', advance='no') Iter + PreviousCycles

            if (tEN2Init) then
                write(var_unit, '(4(3x,es20.13))', advance='no') var_e_num_all, en2_pert_all, &
                                                                  var_e_num_all+en2_pert_all, overlap_all
            else
                write(var_unit, '(2(3x,es20.13))', advance='no') var_e_num_all, overlap_all
            end if
        end if

    end subroutine calc_var_energy

    subroutine get_proj_energy(ValidSpawned, proj_energy)

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

    end subroutine get_proj_energy

    subroutine calc_precond_energy(ValidSpawned, proj_energy)

        use CalcData, only: tau, tEN2Init
        use global_det_data, only: det_diagH
        use semi_stoch_procs, only: core_space_pos, check_determ_flag

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, PartInd, DetHash, determ_pos, ierr
        integer :: nI_spawn(nel)
        real(dp) :: SpawnedSign(lenof_sign), CurrentSign(lenof_sign)
        real(dp) :: h_diag
        real(dp) :: precond_e_num(lenof_sign), precond_denom(lenof_sign)
        real(dp) :: precond_e_num_all(lenof_sign), precond_denom_all(lenof_sign)
        real(dp) :: precond_e_num_av, precond_denom_av
        logical :: tCoreDet, tSuccess, abort(lenof_sign)

        precond_e_num = 0.0_dp
        precond_denom = 0.0_dp
        precond_e_num_all = 0.0_dp
        precond_denom_all = 0.0_dp

        tCoreDet = .false.
        tSuccess = .false.
        abort = .false.

        tSpawnedTo = .false.

        do i = 1, ValidSpawned
            call decode_bit_det(nI_spawn, SpawnedParts(:,i))
            call extract_sign(SpawnedParts(:,i), SpawnedSign)

            ! Now add in the diagonal elements
            call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            if (tSuccess) then
                call extract_sign(CurrentDets(:,PartInd), CurrentSign)
                h_diag = det_diagH(PartInd) + Hii

                tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                if (tCoreDet) then
                    determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                    tSpawnedTo(determ_pos) = .true.
                    SpawnedSign = SpawnedSign + partial_determ_vecs(:, determ_pos)
                end if

                precond_e_num(1) = precond_e_num(1) - &
                    SpawnedSign(1) * h_diag * CurrentSign(2) / (tau * (proj_energy(1) + Hii - h_diag ))

                precond_e_num(2) = precond_e_num(2) - &
                    SpawnedSign(2) * h_diag * CurrentSign(1) / (tau * (proj_energy(2) + Hii - h_diag ))

                precond_denom(1) = precond_denom(1) - &
                    SpawnedSign(1) * CurrentSign(2) / (tau * (proj_energy(1) + Hii - h_diag ))

                precond_denom(2) = precond_denom(2) - &
                    SpawnedSign(2) * CurrentSign(1) / (tau * (proj_energy(2) + Hii - h_diag ))
            end if

            ! Add in the contributions corresponding to off-diagonal
            ! elements of the Hamiltonian
            if (abs(SpawnedSign(1)) > 1.e-12_dp .and. abs(SpawnedSign(2)) > 1.e-12_dp) then
                ! If we don't have h_diag already...
                if (.not. tSuccess) then
                    if (tHPHF) then
                        h_diag = hphf_diag_helement(nI_spawn, SpawnedParts(:, i))
                    else
                        h_diag = get_helement(nI_spawn, nI_spawn, 0)
                    end if
                end if

                precond_e_num(1) = precond_e_num(1) + &
                  SpawnedSign(1) * SpawnedSign(2) / ((tau**2) * ( proj_energy(1) + Hii - h_diag ))

                precond_e_num(2) = precond_e_num(2) + &
                  SpawnedSign(1) * SpawnedSign(2) / ((tau**2) * ( proj_energy(2) + Hii - h_diag ))
            end if
        end do

        ! Contribution from deterministic states that were not spawned to
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                if (.not. tSpawnedTo(i)) then
                    call extract_sign(CurrentDets(:, indices_of_determ_states(i)), CurrentSign)

                    precond_e_num(1) = precond_e_num(1) + &
                      partial_determ_vecs(1, i) * partial_determ_vecs(2, i) / &
                        ( (tau**2) * (proj_energy(1) - core_ham_diag(i)) )

                    precond_e_num(2) = precond_e_num(2) + &
                      partial_determ_vecs(1, i) * partial_determ_vecs(2, i) / &
                       ( (tau**2) * (proj_energy(2) - core_ham_diag(i)) )

                    precond_e_num(1) = precond_e_num(1) - &
                      partial_determ_vecs(1, i) * (core_ham_diag(i) + Hii) * CurrentSign(2) / &
                        ( tau * (proj_energy(1) - core_ham_diag(i)) )

                    precond_e_num(2) = precond_e_num(2) - &
                      partial_determ_vecs(2, i) * (core_ham_diag(i) + Hii) * CurrentSign(1) / &
                       ( tau * (proj_energy(2) - core_ham_diag(i)) )

                    precond_denom(1) = precond_denom(1) - &
                      partial_determ_vecs(1, i) * CurrentSign(2) / &
                        ( tau * (proj_energy(1) - core_ham_diag(i)) )

                    precond_denom(2) = precond_denom(2)- &
                      partial_determ_vecs(2, i) * CurrentSign(1) / &
                       ( tau * (proj_energy(2) - core_ham_diag(i)) )
                end if
            end do
        end if

        call MPIBarrier(ierr)
        call MPISumAll(precond_e_num, precond_e_num_all)
        call MPISumAll(precond_denom, precond_denom_all)

        precond_e_num_av = ( precond_e_num_all(1) + precond_e_num_all(2) ) / 2.0_dp
        precond_denom_av = ( precond_denom_all(1) + precond_denom_all(2) ) / 2.0_dp

        if (iProcIndex == 0) then
            write(var_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(1), precond_denom_all(1)
            write(var_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(2), precond_denom_all(2)
            write(var_unit, '(2(3x,es20.13))', advance='no') precond_e_num_av, precond_denom_av
            write(var_unit,'()')
        end if

    end subroutine calc_precond_energy

    subroutine round_spawns(ValidSpawned, cutoff)

        integer, intent(inout) :: ValidSpawned
        real(dp), intent(in) :: cutoff

        integer :: i, j, new_length
        real(dp) :: SpawnedSign(lenof_sign), prob_remove, r

        write(6,*) "Number of spawns before round:", ValidSpawned

        new_length = 0

        do i = 1, ValidSpawned
            call extract_sign(SpawnedParts(:,i), SpawnedSign)

            do j = 1, lenof_sign
                if ( abs(SpawnedSign(j)) > 1.e-12_dp .and. abs(SpawnedSign(j)) < cutoff ) then
                    prob_remove = (cutoff - abs(SpawnedSign(j)))/cutoff
                    r = genrand_real2_dSFMT()

                    if (prob_remove > r) then
                        SpawnedSign(j) = 0.0_dp
                    else
                        SpawnedSign(j) = sign(cutoff, SpawnedSign(j))
                    end if

                    call encode_part_sign(SpawnedParts(:,i), SpawnedSign(j), j)
                end if
            end do

            if (.not. IsUnoccDet(SpawnedSign)) then
                new_length = new_length + 1
                SpawnedParts(:,new_length) = SpawnedParts(:,i)
            end if

        end do

        ValidSpawned = new_length

        write(6,*) "Number of spawns after round:", ValidSpawned

    end subroutine round_spawns

    subroutine rescale_spawns(ValidSpawned, proj_energy)

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, nJ(nel)
        real(dp) :: SpawnedSign(lenof_sign), hdiag

        ! Find the weight spawned on the Hartree--Fock determinant.
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                partial_determ_vecs(:,i) = partial_determ_vecs(:,i) / (core_ham_diag(i) - proj_energy)
            end do
        end if

        do i = 1, ValidSpawned
            call decode_bit_det(nJ, SpawnedParts(:,i))

            if (tHPHF) then
                HDiag = hphf_diag_helement (nJ, SpawnedParts(:,i))
            else
                HDiag = get_helement (nJ, nJ, 0)
            end if

            call extract_sign(SpawnedParts(:,i), SpawnedSign)
            SpawnedSign = SpawnedSign / (hdiag - proj_energy - Hii)
            call encode_sign(SpawnedParts(:,i), SpawnedSign)
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

    subroutine time_hii(ValidSpawned)

        integer, intent(in) :: ValidSpawned

        integer :: i, nI(nel)
        real(dp) :: hdiag

        do i = 1, ValidSpawned
            call decode_bit_det(nI, SpawnedParts(:,i))
            if (tHPHF) then
                HDiag = hphf_diag_helement (nI, SpawnedParts(:,i))
            else
                HDiag = get_helement (nI, nI, 0)
            end if
        end do

    end subroutine time_hii

end module precond_annihilation_mod
