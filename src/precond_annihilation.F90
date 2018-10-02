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
        real(dp) :: proj_energy(lenof_sign)
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

        !call set_timer(hii_test_time, 30)
        !call time_hii(MaxIndex)
        !call halt_timer(hii_test_time)

        call set_timer(proj_e_time, 30)
        call get_proj_energy(MaxIndex, proj_energy, ref_positions)
        call halt_timer(proj_e_time)

        spawn_hii(1:MaxIndex) = 0.0_dp

        call set_timer(precond_e_time, 30)
        call calc_e_and_set_init_flags(MaxIndex, proj_energy)
        call halt_timer(precond_e_time)

        call set_timer(precond_round_time, 30)
        call round_spawns(MaxIndex, precondSpawnCutoff, ref_positions)
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

        use CalcData, only: tau, tEN2Init
        use global_det_data, only: det_diagH
        use semi_stoch_procs, only: core_space_pos, check_determ_flag

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, j, PartInd, DetHash, determ_pos, ierr
        integer :: nI_spawn(nel)
        real(dp) :: SpawnedSign(lenof_sign), CurrentSign(lenof_sign)
        real(dp) :: h_diag, mean_energy(lenof_sign)
        logical :: tCoreDet, tSuccess, tCalcHii, abort(lenof_sign)

        real(dp) :: var_e_num(lenof_sign), overlap(lenof_sign)
        real(dp) :: var_e_num_all(lenof_sign), overlap_all(lenof_sign)
        real(dp) :: en2_pert(lenof_sign), en2_pert_all(lenof_sign)
        real(dp) :: precond_e_num(lenof_sign), precond_denom(lenof_sign)
        real(dp) :: precond_e_num_all(lenof_sign), precond_denom_all(lenof_sign)
        real(dp) :: precond_e_num_av, precond_denom_av

        ! Allow room to send up to 1000 elements.
        real(dp) :: send_arr(5*lenof_sign)
        ! Allow room to receive up to 1000 elements.
        real(dp) :: recv_arr(5*lenof_sign)

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

        mean_energy(1) = Hii + ( proj_energy(1) + proje_ref_energy_offsets(1) + &
                                 proj_energy(2) + proje_ref_energy_offsets(2) ) / 2.0_dp

        tCoreDet = .false.
        tSuccess = .false.
        abort = .false.

        tSpawnedTo = .false.

        ! Contribution from diagonal for variational energy
        do i = 1, int(TotWalkers, sizeof_int)
            h_diag = det_diagH(i) + Hii
            call extract_sign(CurrentDets(:, i), CurrentSign)
            overlap(1) = overlap(1) + CurrentSign(1) * CurrentSign(2)
            var_e_num(1) = var_e_num(1) + h_diag * CurrentSign(1) * CurrentSign(2)
        end do

        do i = 1, ValidSpawned
            tCalcHii = .false.
            call decode_bit_det(nI_spawn, SpawnedParts(:,i))
            call extract_sign(SpawnedParts(:,i), SpawnedSign)

            ! Now add in the diagonal elements
            call hash_table_lookup(nI_spawn, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            if (tSuccess) then
                call extract_sign(CurrentDets(:,PartInd), CurrentSign)

                ! Set initiator flags for the spawning, before the currently
                ! occupied determinant is potentially killed in the death step.
                do j = 1, lenof_sign
                    if (abs(CurrentSign(j)) > 1.e-12_dp) then
                        call set_flag(SpawnedParts(:,i), get_initiator_flag(j))
                    end if
                end do

                h_diag = det_diagH(PartInd) + Hii
                tCalcHii = .true.

                tCoreDet = check_determ_flag(CurrentDets(:,PartInd))
                if (tCoreDet) then
                    determ_pos = core_space_pos(SpawnedParts(:,i), nI_spawn) - determ_displs(iProcIndex)
                    tSpawnedTo(determ_pos) = .true.
                    SpawnedSign = SpawnedSign + partial_determ_vecs(:, determ_pos)
                end if

                var_e_num(1) = var_e_num(1) - &
                    (SpawnedSign(1) * CurrentSign(2) + SpawnedSign(2) * CurrentSign(1))/(2.0_dp * tau)

                precond_e_num(1) = precond_e_num(1) - &
                    SpawnedSign(1) * h_diag * CurrentSign(2) / ( (proj_energy(1) + proje_ref_energy_offsets(1) + Hii - h_diag ))

                precond_e_num(2) = precond_e_num(2) - &
                    SpawnedSign(2) * h_diag * CurrentSign(1) / ( (proj_energy(2) + proje_ref_energy_offsets(2) + Hii - h_diag ))

                precond_denom(1) = precond_denom(1) - &
                    SpawnedSign(1) * CurrentSign(2) / ( (proj_energy(1) + proje_ref_energy_offsets(1) + Hii - h_diag ))

                precond_denom(2) = precond_denom(2) - &
                    SpawnedSign(2) * CurrentSign(1) / ( (proj_energy(2) + proje_ref_energy_offsets(2) + Hii - h_diag ))
            end if

            ! Add in the contributions corresponding to off-diagonal
            ! elements of the Hamiltonian
            if (abs(SpawnedSign(1)) > 1.e-12_dp .and. abs(SpawnedSign(2)) > 1.e-12_dp) then
                tCalcHii = .true.
                ! If we don't have h_diag already...
                if (.not. tSuccess) then
                    if (tHPHF) then
                        h_diag = hphf_diag_helement(nI_spawn, SpawnedParts(:, i))
                    else
                        h_diag = get_helement(nI_spawn, nI_spawn, 0)
                    end if
                end if

                precond_e_num(1) = precond_e_num(1) + &
                  SpawnedSign(1) * SpawnedSign(2) / (tau * ( proj_energy(1) + proje_ref_energy_offsets(1) + Hii - h_diag ))

                precond_e_num(2) = precond_e_num(2) + &
                  SpawnedSign(1) * SpawnedSign(2) / (tau * ( proj_energy(2) + proje_ref_energy_offsets(2) + Hii - h_diag ))

                ! Only get EN2 contribution is we're due to cancel this
                ! spawning on both replicas
                if (tEN2Init .and. (.not. tSuccess)) then
                    do j = 1, lenof_sign
                        abort(j) = test_abort_spawn(SpawnedParts(:, i), j)
                    end do

                    if (abort(1) .and. abort(2)) then
                        en2_pert(1) = en2_pert(1) + &
                          SpawnedSign(1)*SpawnedSign(2) / ( (tau**2) * (mean_energy(1) - h_diag ) )
                    end if
                end if

            end if

            ! Store, to avoid recalculating later
            if (tCalcHii) spawn_hii(i) = h_diag
        end do

        ! Contribution from deterministic states that were not spawned to
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                if (.not. tSpawnedTo(i)) then
                    call extract_sign(CurrentDets(:, indices_of_determ_states(i)), CurrentSign)

                    var_e_num(1) = var_e_num(1) - &
                        (partial_determ_vecs(1, i) * CurrentSign(2) + &
                         partial_determ_vecs(2, i) * CurrentSign(1)) / (2.0_dp * tau)

                    precond_e_num(1) = precond_e_num(1) + &
                      partial_determ_vecs(1, i) * partial_determ_vecs(2, i) / &
                        ( tau * (proj_energy(1) + proje_ref_energy_offsets(1) - core_ham_diag(i)) )

                    precond_e_num(2) = precond_e_num(2) + &
                      partial_determ_vecs(1, i) * partial_determ_vecs(2, i) / &
                       ( tau * (proj_energy(2) + proje_ref_energy_offsets(2) - core_ham_diag(i)) )

                    precond_e_num(1) = precond_e_num(1) - &
                      partial_determ_vecs(1, i) * (core_ham_diag(i) + Hii) * CurrentSign(2) / &
                        ( (proj_energy(1) + proje_ref_energy_offsets(1) - core_ham_diag(i)) )

                    precond_e_num(2) = precond_e_num(2) - &
                      partial_determ_vecs(2, i) * (core_ham_diag(i) + Hii) * CurrentSign(1) / &
                       ( (proj_energy(2) + proje_ref_energy_offsets(2) - core_ham_diag(i)) )

                    precond_denom(1) = precond_denom(1) - &
                      partial_determ_vecs(1, i) * CurrentSign(2) / &
                        ( (proj_energy(1) + proje_ref_energy_offsets(1) - core_ham_diag(i)) )

                    precond_denom(2) = precond_denom(2) - &
                      partial_determ_vecs(2, i) * CurrentSign(1) / &
                       ( (proj_energy(2) + proje_ref_energy_offsets(2) - core_ham_diag(i)) )
                end if
            end do
        end if

        ! ---- MPI communication --------------------------------
        send_arr(0*lenof_sign+1:1*lenof_sign) = var_e_num
        send_arr(1*lenof_sign+1:2*lenof_sign) = overlap
        send_arr(2*lenof_sign+1:3*lenof_sign) = en2_pert
        send_arr(3*lenof_sign+1:4*lenof_sign) = precond_e_num
        send_arr(4*lenof_sign+1:5*lenof_sign) = precond_denom

        call MPIBarrier(ierr)
        call MPISumAll(send_arr, recv_arr)

        var_e_num_all     = recv_arr(0*lenof_sign+1:1*lenof_sign)
        overlap_all       = recv_arr(1*lenof_sign+1:2*lenof_sign)
        en2_pert_all      = recv_arr(2*lenof_sign+1:3*lenof_sign)
        precond_e_num_all = recv_arr(3*lenof_sign+1:4*lenof_sign)
        precond_denom_all = recv_arr(4*lenof_sign+1:5*lenof_sign)
        ! -------------------------------------------------------

        precond_e_num_av = ( precond_e_num_all(1) + precond_e_num_all(2) ) / 2.0_dp
        precond_denom_av = ( precond_denom_all(1) + precond_denom_all(2) ) / 2.0_dp

        if (iProcIndex == 0) then
            write(var_unit, '(1x,i13)', advance='no') Iter + PreviousCycles
            write(var_unit, '(1(3x,es20.13))', advance='no') var_e_num_all(1)
            if (tEN2Init) then
                write(var_unit, '(2(3x,es20.13))', advance='no') en2_pert_all(1), var_e_num_all(1)+en2_pert_all(1)
            end if
            write(var_unit, '(1(3x,es20.13))', advance='no') overlap_all(1)
            write(var_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(1), precond_denom_all(1)
            write(var_unit, '(2(3x,es20.13))', advance='no') precond_e_num_all(2), precond_denom_all(2)
            write(var_unit, '(2(3x,es20.13))', advance='no') precond_e_num_av, precond_denom_av
            write(var_unit,'()')
        end if

    end subroutine calc_e_and_set_init_flags

    subroutine round_spawns(ValidSpawned, cutoff, ref_positions)

        integer, intent(inout) :: ValidSpawned
        real(dp), intent(in) :: cutoff
        integer, intent(in) :: ref_positions(lenof_sign)

        integer :: i, j, new_length
        real(dp) :: SpawnedSign(lenof_sign), prob_remove, r

        write(6,*) "Number of spawns before round:", ValidSpawned

        new_length = 0

        do i = 1, ValidSpawned
            call extract_sign(SpawnedParts(:,i), SpawnedSign)

            ! Don't perform rounding on reference determinants
            if (.not. any(ref_positions == i)) then

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

            end if

            if (.not. IsUnoccDet(SpawnedSign)) then
                new_length = new_length + 1
                SpawnedParts(:,new_length) = SpawnedParts(:,i)
                spawn_hii(new_length) = spawn_hii(i)
            end if
        end do

        ValidSpawned = new_length

        write(6,*) "Number of spawns after round:", ValidSpawned

    end subroutine round_spawns

    subroutine rescale_spawns(ValidSpawned, proj_energy)

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)

        integer :: i, nJ(nel)
        real(dp) :: SpawnedSign(lenof_sign), h_diag

        ! Find the weight spawned on the Hartree--Fock determinant.
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                partial_determ_vecs(:,i) = partial_determ_vecs(:,i) / &
                  (core_ham_diag(i) - proj_energy -  proje_ref_energy_offsets)
            end do
        end if

        do i = 1, ValidSpawned
            if (abs(spawn_hii(i)) < 1e-12_dp) then
                call decode_bit_det(nJ, SpawnedParts(:,i))

                if (tHPHF) then
                    h_diag = hphf_diag_helement (nJ, SpawnedParts(:,i))
                else
                    h_diag = get_helement (nJ, nJ, 0)
                end if
            else
                h_diag = spawn_hii(i)
            end if

            call extract_sign(SpawnedParts(:,i), SpawnedSign)
            SpawnedSign = SpawnedSign / (h_diag - proj_energy - proje_ref_energy_offsets - Hii)
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
