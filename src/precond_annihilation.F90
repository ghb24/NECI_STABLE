#include "macros.h"

module precond_annihilation_mod

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use AnnihilationMod, only: test_abort_spawn
    use AnnihilationMod, only: AnnihilateSpawnedParts, deterministic_annihilation
    use bit_rep_data
    use bit_reps, only: encode_sign, set_flag, get_initiator_flag
    use CalcData, only: OccupiedThresh, tTruncInitiator
    use constants, only: n_int, lenof_sign, null_part, sizeof_int
    use determinants, only: get_helement
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData
    use hash
    use hphf_integrals, only: hphf_diag_helement
    use Parallel_neci
    use load_balance, only: CalcHashTableStats
    use searching
    use semi_stoch_procs, only: fill_in_diag_helements, is_core_state
    use sort_mod
    use SystemData, only: NEl, tHPHF

    implicit none

    contains

    subroutine precond_annihilation(TotWalkersNew, iter_data, tSingleProc)

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        logical, intent(in) :: tSingleProc

        integer :: MaxIndex
        integer(n_int), pointer :: PointTemp(:,:)
        type(timer), save :: Compress_time
        real(dp) :: proj_energy(lenof_sign)

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

        call get_proj_energy(MaxIndex, proj_energy)

        call rescale_spawns(MaxIndex, proj_energy)

        call perform_death_with_precond(iter_data)

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

    subroutine get_proj_energy(ValidSpawned, proj_energy)

        use CalcData, only: tau

        integer, intent(in) :: ValidSpawned
        real(dp), intent(out) :: proj_energy(lenof_sign)

        integer :: i
        real(dp) :: SignTemp(lenof_sign)
        logical :: ref_found, tSuccess
        integer :: PartInd, DetHash

        proj_energy = 0.0_dp
        ref_found = .false.

        ! Find the weight spawned on the Hartree--Fock determinant.
        if (tSemiStochastic) then
            do i = 1, determ_sizes(iProcIndex)
                if (DetBitEQ(core_space(0:NIfDBO,determ_displs(iProcIndex)+i), iLutHF, NIfDBO)) then
                    proj_energy = -partial_determ_vecs(:,i)
                    ref_found = .true.
                    exit
                end if
            end do
        end if

        do i = 1, ValidSpawned
            if (DetBitEQ(SpawnedParts(:,i), iLutHF, NIfDBO)) then
                call extract_sign(SpawnedParts(:,i), SignTemp)
                proj_energy = proj_energy - SignTemp
                ref_found = .true.
                exit
            end if
        end do

        if (iProcIndex == iRefProc(1)) then
            if ( (.not. ref_found) ) then
                proj_energy = -0.01_dp
                write(6,'("WARNING: The reference determinant was not spawned to in the last iteration.")')
            else if (abs(proj_energy(1)) < 1.e-12_dp) then
                proj_energy = -0.01_dp
                write(6,'("WARNING: The projected energy from the last iteration was zero. Setting to -0.1.")')
            end if

            call hash_table_lookup(HFDet, ilutHF, NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            if (tSuccess) then
                call extract_sign(CurrentDets(:,PartInd), InstNoAtHF)
                proj_energy = proj_energy/InstNoAtHF
            else
                write(6,'("WARNING: HF determinant not found in main walker list.")')
            end if
        end if

        ! Remove time step
        proj_energy = proj_energy / tau

        call MPIBCast(proj_energy, lenof_sign, iRefProc(1))

    end subroutine get_proj_energy

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

        ! Calculate an instantaneous value of the projected energy for the
        ! given walkers distributions

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

end module precond_annihilation_mod
