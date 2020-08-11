#include "macros.h"

! module containing useful functions and subroutines needed in the real-time
! implementation of the FCIQMC algotrithm

module real_time_procs
    use hash, only: hash_table_lookup, init_hash_table, clear_hash_table, &
                    add_hash_table_entry, fill_in_hash_table
    use SystemData, only: nel, nBasis, tHPHF
    use real_time_data, only: gf_overlap, TotWalkers_orig, overlap_states, tInfInit, &
                              t_complex_ints, real_time_info, temp_freeslot, dyn_norm_red, &
                              temp_det_list, temp_det_pointer, temp_iendfreeslot, &
                              temp_det_hash, temp_totWalkers, pert_norm, allGfs, &
                              valid_diag_spawns, DiagParts, n_diag_spawned, tOverpopulate, &
                              gf_count, tVerletSweep, t_real_time_fciqmc, &
                              t_rotated_time, tau_imag, tau_real, gs_energy, TotPartsLastAlpha, &
                              shift_damping, normsize, tStabilizerShift, dyn_norm_psi, &
                              TotPartsPeak, numCycShiftExcess, shiftLimit, t_kspace_operators, &
                              tDynamicAlpha, tDynamicDamping, stepsAlpha, phase_factors, &
                              elapsedImagTime, elapsedRealTime, tStaticShift, asymptoticShift, &
                              iunitCycLog, trajFile, tauCache, alphaCache, tNewOverlap, &
                              alphaLog, alphaLogSize, alphaLogPos, tOnlyPositiveShift, &
                              tHFOverlap
    use real_time_aux, only: write_overlap_state, write_overlap_state_serial

    use kp_fciqmc_data_mod, only: perturbed_ground, overlap_pert
    use constants, only: dp, lenof_sign, int64, n_int, EPS, iout, null_part, &
                         sizeof_int, MPIArg
    use bit_reps, only: decode_bit_det, test_flag, encode_sign, &
                        set_flag, encode_bit_rep, extract_bit_rep, &
                        flag_deterministic, encode_part_sign, &
                        get_initiator_flag, get_initiator_flag_by_run, &
                        clr_flag, test_flag_multi
    use util_mod, only: get_free_unit, get_unique_filename, near_zero, &
                        operator(.isclose.)
    use bit_rep_data, only: extract_sign, nifd, niftot
    use FciMCData, only: CurrentDets, HashIndex, popsfile_dets, MaxWalkersPart, &
                         WalkVecDets, freeslot, spawn_ht, nhashes_spawn, MaxSpawned, &
                         iStartFreeSlot, iEndFreeSlot, ValidSpawnedList, &
                         InitialSpawnedSlots, iLutRef, inum_runs, max_cyc_spawn, core_run, &
                         tSearchTau, tFillingStochRDMonFly, fcimc_iter_data, &
                         NoAddedInitiators, SpawnedParts, acceptances, TotWalkers, &
                         nWalkerHashes, iter, fcimc_excit_gen_store, NoDied, &
                         NoBorn, NoAborted, NoRemoved, HolesInList, TotParts, Hii, &
                         tSinglePartPhase, perturbation, alloc_popsfile_dets
    use core_space_util, only: cs_replicas
    use perturbations, only: apply_perturbation, init_perturbation_creation, &
                             init_perturbation_annihilation, apply_perturbation_array
    use util_mod, only: int_fmt
    use CalcData, only: AvMCExcits, tAllRealCoeff, tRealCoeffByExcitLevel, &
                        tRealSpawnCutoff, RealSpawnCutoff, tau, RealCoeffExcitThresh, &
                        DiagSft, tTruncInitiator, OccupiedThresh, tReadPops, InitiatorWalkNo, &
                        tSpinProject
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use procedure_pointers, only: get_spawn_helement
    use util_mod, only: stochastic_round
    use tau_search, only: log_spawn_magnitude
    use rdm_general, only: calc_rdmbiasfac
    use global_det_data, only: global_determinant_data
! RT_M_Merge: Disabled rdms
!    use rdm_data, only: nrdms, rdms
    use hash, only: remove_hash_table_entry
    use dSFMT_interface, only: genrand_real2_dSFMT
    use load_balance_calcnodes, only: DetermineDetNode
    use Determinants, only: tDefineDet, DefDet
    use ParallelHelper, only: nNodes, bNodeRoot, ProcNode, NodeRoots, MPIBarrier, &
                              iProcIndex, MPI_SUM, root
    use Parallel_neci
    use LoggingData, only: tNoNewRDMContrib
    use AnnihilationMod, only: test_abort_spawn
    use load_balance, only: AddNewHashDet, CalcHashTableStats, get_diagonal_matel
    use semi_stoch_gen, only: generate_space_most_populated, reset_core_space
    use semi_stoch_procs, only: GLOBAL_RUN
    implicit none

    type(timer) :: calc_gf_time

contains

!------------------------------------------------------------------------------------------!

    subroutine DirectAnnihilation_diag(TotWalkersNew, iter_data)
        ! new direct annihilation routine to mimick the diagonal death step
        ! in the y(n) + k2 combination between reloaded CurrentDets and the
        ! DiagParts list
        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "DirectAnnihilation_diag"
        integer :: numSpawns

        ! As only diagonal events are considered here, no communication
        ! is required
        ! Also, this eliminates the need for compression as all dets are
        ! already stored contigously and in the right orde
        ! (no annihilation inside DiagParts can occur)

        ! valid_diag_spawns gets increased after the spawn
        ! to DiagParts(:,valid_diag_spawns) -> highest index is
        ! actually valid_diag_spawns-1

        numSpawns = valid_diag_spawns - 1
        call AnnihilateDiagParts(numSpawns, TotWalkersNew, iter_data)

        ! also should update the hashtable stats, specific for this diagonal
        ! spawning event, but the original one should work also for this
        ! since it only takes CurrentDets into account!

        call CalcHashTableStats(TotWalkersNew, iter_data)

        ! this should be it, deterministic annihilation is carried out in the next
        ! step, within the 'regular' annihilation

    end subroutine DirectAnnihilation_diag

!------------------------------------------------------------------------------------------!

    subroutine AnnihilateDiagParts(ValidSpawned, TotWalkersNew, iter_data)
        ! this is the new "annihilation" routine which mimics the actual
        ! diagonal death-step, between the reloaded CurrentDets y(n) and the
        ! diagonal parts in k2, DiagParts
        integer, intent(inout) :: ValidSpawned, TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "AnnihilateDiagParts"

        integer :: PartInd, i, j
        real(dp), dimension(lenof_sign) :: CurrentSign, SpawnedSign, SignTemp
        real(dp), dimension(lenof_sign) :: SignProd
        real(dp) :: HDiag
        integer :: DetHash, nJ(nel)
        logical :: tSuccess, tDetermState
        integer :: run, err, pos

        ! rewrite the original Annihilation routine to fit the new
        ! requirements here

        ! this routine updated the NoDied variables etc.. for the 2nd RK step
        ! so i dont think i need to change much here
        ! Only node roots to do this.
        if (.not. bNodeRoot) return

        do i = 1, ValidSpawned

            call decode_bit_det(nJ, DiagParts(:, i))

            ! Search the hash table HashIndex for the determinant defined by
            ! nJ and DiagParts(:,i). If it is found, tSuccess will be
            ! returned .true. and PartInd will hold the position of the
            ! determinant in CurrentDets. Else, tSuccess will be returned
            ! .false. (and PartInd shouldn't be accessed).
            ! Also, the hash value, DetHash, is returned by this routine.
            ! tSuccess will determine whether the particle has been found or not.
            call hash_table_lookup(nJ, DiagParts(:, i), nifd, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            tDetermState = .false.
            CurrentSign = 0.0_dp

!            write(6,*) 'i,DiagParts(:,i)',i,DiagParts(:,i)

            if (tSuccess) then

                ! Our DiagParts determinant is found in CurrentDets.

                call extract_sign(CurrentDets(:, PartInd), CurrentSign)
                call extract_sign(DiagParts(:, i), SpawnedSign)

                SignProd = CurrentSign * SpawnedSign

                tDetermState = test_flag_multi(CurrentDets(:, PartInd), flag_deterministic)

                if (sum(abs(CurrentSign)) >= 1.e-12_dp .or. tDetermState) then
                    ! Transfer new sign across.
                    call encode_sign(CurrentDets(:, PartInd), SpawnedSign + CurrentSign)
                    ! I dont see why DiagParts has to be nullified. In the verlet scheme
                    ! warmup, we might want to use DiagParts afterwards
                    call encode_sign(DiagParts(:, i), null_part)
                    do j = 1, lenof_sign
                        run = part_type_to_run(j)
                        if (is_run_unnocc(CurrentSign, run)) then
                            ! This determinant is actually *unoccupied* for the
                            ! walker type/set we're considering. We need to
                            ! decide whether to abort it or not.
                            ! rt_iter_adapt : also count those spawned onto an initiator
                            ! this is not necessary in the normal version as walkers are
                            ! counted there on spawn
                            ! also, they are born for any sign
                            iter_data%nborn(j) = iter_data%nborn(j) + &
                                                 abs(SpawnedSign(j))
                            NoBorn(run) = NoBorn(run) + abs(SpawnedSign(j))

                        else if (SignProd(j) < 0) then
                            ! in the real-time for the final combination
                            ! y(n) + k2 i have to check if the "spawned"
                            ! particle is actually a diagonal death/born
                            ! walker
                            ! remember the - sign when filling up the DiagParts
                            ! list -> so opposite sign here means a death!
                            iter_data%ndied(j) = iter_data%ndied(j) + &
                                                 min(abs(CurrentSign(j)), abs(SpawnedSign(j)))

                            NoDied(run) = NoDied(run) + min(abs(CurrentSign(j)), &
                                                            abs(SpawnedSign(j)))
                            ! and if the Spawned sign magnitude is even higher
                            ! then the currentsign -> born anti-particles

                            iter_data%nborn(j) = iter_data%nborn(j) + &
                                                 max(abs(SpawnedSign(j)) - abs(CurrentSign(j)), 0.0_dp)

                            NoBorn(run) = NoBorn(run) + max(abs(SpawnedSign(j)) - &
                                                            abs(CurrentSign(j)), 0.0_dp)
                        else
                            ! if it has the same sign i have to keep track of
                            ! the born particles, or as in the original
                            ! walker_death, reduce the number of died parts
                            iter_data%ndied(j) = iter_data%ndied(j) - &
                                                 abs(SpawnedSign(j))

                        end if

                    end do ! Over all components of the sign.

                    if (.not. tDetermState) then
                        call extract_sign(CurrentDets(:, PartInd), SignTemp)
                        if (IsUnoccDet(SignTemp)) then
                            ! All walkers in this main list have been annihilated
                            ! away. Remove it from the hash index array so that
                            ! no others find it (it is impossible to have another
                            ! spawned walker yet to find this determinant).
                            call remove_hash_table_entry(HashIndex, nJ, PartInd)
                            ! Add to "freeslot" list so it can be filled in.
                            iEndFreeSlot = iEndFreeSlot + 1
                            FreeSlot(iEndFreeSlot) = PartInd
                        end if
                    end if

                    ! RT_M_Merge: Disabled rdm functionality
                    ! We must use the instantaneous value for the off-diagonal
                    ! contribution. However, we can't just use CurrentSign from
                    ! the previous iteration, as this has been subject to death
                    ! but not the new walkers. We must add on SpawnedSign, so
                    ! we're effectively taking the instantaneous value from the
                    ! next iter. This is fine as it's from the other population,
                    ! and the Di and Dj signs are already strictly uncorrelated.

                end if

            end if

            if (((.not. tSuccess) .or. (tSuccess .and. sum(abs(CurrentSign)) < 1.e-12_dp .and. (.not. tDetermState)))) then

                ! Running the full, non-initiator scheme.
                ! Determinant in newly spawned list is not found in
                ! CurrentDets. If coeff <1, apply removal criterion.
                call extract_sign(DiagParts(:, i), SignTemp)

                if (.not. IsUnoccDet(SignTemp)) then
                    ! Walkers have not been aborted and so we should copy the
                    ! determinant straight over to the main list. We do not
                    ! need to recompute the hash, since this should be the
                    ! same one as was generated at the beginning of the loop.
                    ! also here treat those new walkers as born particles

                    iter_data%nborn = iter_data%nborn + abs(SignTemp)
                    do run = 1, inum_runs
                        NoBorn(run) = NoBorn(run) + sum(abs(SignTemp( &
                                                            min_part_type(run):max_part_type(run))))
                    end do
                    HDiag = get_diagonal_matel(nJ, DiagParts(:, i))
                    call AddNewHashDet(TotWalkersNew, DiagParts(:, i), DetHash, nJ, HDiag, pos, err)
                    if (err /= 0) exit

                end if
            end if
        end do
        ! Update remaining number of holes in list for walkers stats.
        if (iStartFreeSlot > iEndFreeSlot) then
            ! All slots filled
            HolesInList = 0
        else
            HolesInList = iEndFreeSlot - (iStartFreeSlot - 1)
        end if

    end subroutine AnnihilateDiagParts

!------------------------------------------------------------------------------------------!

    function count_holes_in_currentDets() result(holes)
        integer :: holes
        integer(n_int), pointer :: ilut_parent(:)
        integer :: nI_parent(nel), unused_flags, idet
        real(dp) :: parent_sign(lenof_sign)

        holes = 0

        do idet = 1, int(TotWalkers, sizeof_int)

            ilut_parent => CurrentDets(:, idet)

            call extract_bit_rep(ilut_parent, nI_parent, parent_sign, unused_flags, idet, &
                                 fcimc_excit_gen_store)

            if (IsUnoccDet(parent_sign)) then
                holes = holes + 1
            end if

        end do

    end function count_holes_in_currentDets

!------------------------------------------------------------------------------------------!

    subroutine create_diagonal_as_spawn(ilut, diag_sign, iter_data)
        ! new routine to create diagonal particles into new DiagParts
        ! array to distinguish between spawns and diagonal events in the
        ! combination y(n) + k2
        use Parallel_neci, only: iProcIndex
        integer(n_int), intent(in) :: ilut(0:niftot)
        real(dp), intent(in) :: diag_sign(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "create_diagonal_as_spawn"

        logical :: list_full
        integer, parameter :: flags = 0

        unused_var(iter_data)

        ! Note that this is a diagonal event, no communication is needed

        ! Check that the position described by valid_diag_spawn_list is acceptable.
        ! If we have filled up the memory that would be acceptable, then
        ! kill the calculation hard (i.e. stop_all) with a descriptive
        ! error message.
        list_full = .false.
        if (valid_diag_spawns > MaxWalkersPart) list_full = .true.
        if (list_full) then
            print *, "Attempting to spawn particle onto processor: ", iProcIndex
            print *, "No memory slots available for this spawn."
            print *, "Please increase MEMORYFACSPAWN"
            print *, "VALID DIAG SPAWN LIST", valid_diag_spawns
            print *, "NUMBER OF OCC DETERMINANTS", TotWalkers
            call stop_all(this_routine, "Out of memory for spawned particles")
        end if

        call encode_bit_rep(DiagParts(:, valid_diag_spawns), ilut, &
                            diag_sign, flags)

        if (tFillingStochRDMonFly) then
            call stop_all(this_routine, "RDM not yet implemented in the rt-fciqmc!")
            ! We are spawning from ilutI to
            ! SpawnedParts(:,valid_diag_spawn_list(proc)). We want to store the
            ! parent (D_i) with the spawned child (D_j) so that we can add in
            ! Ci.Cj to the RDM later.
            ! The parent is nifd integers long, and stored in the second
            ! part of the SpawnedParts array from NIfTot+1 --> NIfTot+1+nifd

        end if

        valid_diag_spawns = valid_diag_spawns + 1

    end subroutine create_diagonal_as_spawn

!------------------------------------------------------------------------------------------!

    function attempt_die_realtime(Kii, RealwSign, walkExcitLevel) &
        result(ndie)
        ! also write a function, which calculates the new "signs"(weights) of
        ! the real and complex walkers for the diagonal death/cloning step
        ! since i need that for both 1st and 2nd loop of RK, but at different
        ! points
        implicit none
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        real(dp), intent(in) :: Kii
        real(dp), dimension(lenof_sign) :: ndie
        integer, intent(in) :: walkExcitLevel
        integer :: run
        real(dp) :: fac(lenof_sign), rat, r

        character(*), parameter :: this_routine = "attempt_die_realtime"

        ! do it dependent on if there is damping, since otherwise i can save
        ! some effort..

        ! do i need a minus sign here?? just from the convention
        ! in the rest of the neci code? yes!
        ! rmneci_setup: there is no reason to use an imaginary shift
        ! except if when dealing with rotated times)
        ! TODO: Energies should be taken with respect to the N-particle ground state energy

        ! important: the matrix element Kii does not contain the reference energy,
        ! therefore it has to be added manually
        do run = 1, inum_runs
            fac(min_part_type(run)) = tau_real * (Kii + Hii - gs_energy(run))
            fac(max_part_type(run)) = 0.0_dp
        end do

        if (abs(real_time_info%damping) < EPS .and. .not. t_rotated_time) then
            if (any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTau) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do run = 1, lenof_sign
                            fac(run) = min(2.0_dp, fac(run))
                        end do
                    else
                        call stop_all(this_routine, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    do run = 1, inum_runs
                        write(iout, '(1X,f13.7)', advance='no') fac(run)
                    end do
                    write(iout, '()')
                end if
            end if
            do run = 1, inum_runs
                if (((tRealCoeffByExcitLevel .and. (WalkExcitLevel <= RealCoeffExcitThresh)) &
                     .or. tAllRealCoeff)) then
                    ndie(min_part_type(run)) = -fac(min_part_type(run)) &
                                               * realwSign(max_part_type(run))
                    ! and - from Re -> Im
                    ! does this give the correct sign compared to the parent sign?
                    ! additional -1 is added in postprocessing to convert ndie -> nborn
                    ndie(max_part_type(run)) = fac(min_part_type(run)) &
                                               * realwSign(min_part_type(run))

                else
                    ! if not exact i have to round stochastically
                    ! Im -> Re
                    rat = -fac(min_part_type(run)) * RealwSign(max_part_type(run))

                    ndie(min_part_type(run)) = real(int(rat), dp)
                    rat = rat - ndie(min_part_type(run))

                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) ndie(min_part_type(run)) = &
                        ndie(min_part_type(run)) + real(nint(sign(1.0_dp, rat)), dp)

                    ! Re -> Im
                    rat = fac(min_part_type(run)) * RealwSign(min_part_type(run))
                    ndie(max_part_type(run)) = real(int(rat), dp)
                    rat = rat - ndie(max_part_type(run))
                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) ndie(max_part_type(run)) = &
                        ndie(max_part_type(run)) + real(nint(sign(1.0_dp, rat)), dp)
                end if
            end do

            ! here i only have influence from the other type of walkers, which
            ! is already stored in the ndie() array
        else
            ! there is an imaginary energy damping factor.. -> rotated time
            ! so i have mixed influences for the diagonal step
            ! n_i' = -e n_i + H_ii c_i
            ! c_i' = -e c_i - H_ii n_i

            ! temporarily use the 2 entries of fac to store both influences
            ! not sure about the sign of this fac factor but the 2nd entry
            ! as it is implemented right now has to have the opposite sign
            ! of the first on above
            ! the diagnonal matrix elements are always real, so there is
            ! no contribution from tau_real via Kii (and no influence of
            ! tau_imag on fac(min_part_type(run))
            do run = 1, inum_runs
                ! tau_real and tau_imag have opposite signs -> shift changed sign
                ! when moved from real to imaginary part
                fac(max_part_type(run)) = tau_real * (real_time_info%damping) &
                                          + tau_imag * (Kii + Hii - gs_energy(run) - DiagSft(run))
            end do

            ! and also about the fac restrictions.. for now but it here anyway..
            if (any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTau) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do run = 1, lenof_sign
                            fac(run) = min(2.0_dp, fac(run))
                        end do
                    else
                        call stop_all(this_routine, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    write(iout, *) "Warning, diagonal spawn probability > 1. Reduce timestep."
                    do run = 1, inum_runs
                        write(iout, '(1X,f13.7)', advance='no') fac(run)
                    end do
                    write(iout, '()')
                end if
            end if
            do run = 1, inum_runs
                if (((tRealCoeffByExcitLevel .and. (WalkExcitLevel <= RealCoeffExcitThresh)) &
                     .or. tAllRealCoeff) .and. .not. tVerletSweep) then
                    ! this gives a huge overhead in the verlet scheme, as all diagonal
                    ! events are classified as valid -> #spawns ~ #dets
                    ! i exact.. just get the weights, this should also still work.
                    ! but have to exchange the weights to come from the other
                    ! type of particles..
                    ! the number of deaths has +sign from Im -> Re
                    ! can i just add the other contribution here?
                    ! and also include the sign of the parent occupation here
                    ! already.
                    ndie(min_part_type(run)) = -fac(min_part_type(run)) &
                                               * realwSign(max_part_type(run)) - &
                                               fac(max_part_type(run)) * realwSign(min_part_type(run))
                    ! and - from Re -> Im
                    ! does this give the correct sign compared to the parent sign?
                    ndie(max_part_type(run)) = fac(min_part_type(run)) * &
                                               (realwSign(min_part_type(run))) - &
                                               fac(max_part_type(run)) * (RealwSign(max_part_type(run)))

                else
                    ! if not exact i have to round stochastically
                    ! is this ok here to just add the second contribution? todo
                    !  -> Re
                    rat = -fac(min_part_type(run)) * RealwSign(max_part_type(run)) &
                          - fac(max_part_type(run)) * RealwSign(min_part_type(run))

                    ndie(min_part_type(run)) = real(int(rat), dp)
                    rat = rat - ndie(min_part_type(run))

                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) ndie(min_part_type(run)) = ndie(min_part_type(run)) &
                                                                 + real(nint(sign(1.0_dp, rat)), dp)

                    !  -> Im
                    rat = fac(min_part_type(run)) * RealwSign(min_part_type(run)) &
                          - fac(max_part_type(run)) * RealwSign(max_part_type(run))

                    ndie(max_part_type(run)) = real(int(rat), dp)
                    rat = rat - ndie(max_part_type(run))

                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) ndie(max_part_type(run)) = &
                        ndie(max_part_type(run)) + real(nint(sign(1.0_dp, rat)), dp)

                end if
            end do
        end if

    end function attempt_die_realtime

!------------------------------------------------------------------------------------------!

    subroutine walker_death_realtime(iter_data, DetCurr, iLutCurr, Kii, RealwSign, &
                                     DetPosition, walkExcitLevel)
        ! need new walker_death routine for the real-time fciqmc, since i
        ! have complex diagonal influence: due to real-time formulation and
        ! the use of a imaginary energy damping factot -ie
        ! (n_i + ic_i)' = -i(H_ii -E0 - ie)(n_i + ic_i)
        ! n_i' = -e n_i + (H_ii - E0) c_i
        ! c_i' = -e c_i - (H_ii - E0) n_i
        ! so the complex and imaginary walkers mix in the diagonal death
        ! step also! but with opposite sign between Re <-> Im spawn
        ! so the 'death' as annihilation only happens after 2 steps..

        ! i may also need a second version of this, where the diagonal step
        ! is no directly executed on the current list, but builds it into
        ! the spawned list array, since i need the original y(n) list to
        ! combine with k2: y(n+1) = y(n) + k2
        ! i guess that could be done..
        ! this routine is now only used in the first loop of the 2nd order
        ! RK! in the 2nd loop i store the diagonal events in a specific
        ! list, which later gets merged with the original y(n) list with
        ! the Annihilation routine

        integer, intent(in) :: DetCurr(nel)
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        real(dp), intent(in) :: Kii
        integer, intent(in) :: DetPosition
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign) :: ndie
        real(dp), dimension(lenof_sign) :: CopySign
        integer, intent(in) :: walkExcitLevel
        integer :: i, j
        real(dp) :: r

        character(len=*), parameter :: t_r = "walker_death_realtime"

        integer(n_int), pointer :: ilut(:)

        ilut => CurrentDets(:, DetPosition)

        ndie = attempt_die_realtime(Kii, realwSign, walkExcitLevel)

        ! this routine only gets called in the first runge-kutta step ->
        ! so only update the stats for the first here!
        do i = 1, lenof_sign
            ! check if the parent and ndie have the same sign
            if (sign(1.0_dp, RealwSign(i)) .isclose.sign(1.0_dp, ndie(i))) then
                ! then the entries in ndie kill the parent, but only maximally
                ! the already occupying walkers can get killed
                iter_data%ndied(i) = iter_data%ndied(i) + &
                                     abs(min(abs(RealwSign(i)), abs(ndie(i))))

                ! if ndie is bigger than the original occupation i am actually
                ! spawning 'anti-particles' which i have to count as born
                ! and reduce the ndied number.. or not?
                ! hm the old code is actually not counting births, due to
                ! the shift.. interesting.. but just subtracts that from
                ! the ndied quantity...
                iter_data%nborn(i) = iter_data%nborn(i) + &
                                     max(abs(ndie(i)) - abs(RealwSign(i)), 0.0_dp)

            else
                ! if they have opposite sign, as in the original algorithm
                ! reduce the number of ndied by that amount
                iter_data%ndied(i) = iter_data%ndied(i) - abs(ndie(i))

            end if
        end do

        CopySign = RealwSign - nDie

        if (any(.not. near_zero(CopySign))) then
            ! For the hashed walker main list, the particles don't move.
            ! Therefore just adjust the weight.
            call encode_sign(CurrentDets(:, DetPosition), CopySign)
            ! rotated_time_setup: there is only one initiator flag for
            ! both complex and real walkers -> initiator flag is kept
        else

            if (tTruncInitiator) then
                ! All particles on this determinant have gone. If the determinant was an initiator, update the stats
                do j = 1, inum_runs
                    if (test_flag(iLutCurr, get_initiator_flag_by_run(j))) then
                        NoAddedInitiators(j) = NoAddedInitiators(j) - 1
                    end if
                end do
            end if

            ! Remove the determinant from the indexing list
            call remove_hash_table_entry(HashIndex, DetCurr, DetPosition)
            ! Add to the "freeslot" list
            iEndFreeSlot = iEndFreeSlot + 1
            FreeSlot(iEndFreeSlot) = DetPosition
            ! Encode a null det to be picked up
            call encode_sign(CurrentDets(:, DetPosition), null_part)
        end if

    end subroutine walker_death_realtime

!------------------------------------------------------------------------------------------!

    subroutine walker_death_spawn()
        ! this routine is for the 2nd RK step, in which the list k2, which
        ! has to be combined with the original y(n) walker list, is created
        ! since the diagonal death/cloning step cannot be done on the currently
        ! iterated on list y(n) + k1/2, treat the diagonal step, like a spawning
        ! step and store it also in the spawned array
        ! possible considerations: maybe the original spawned list is then
        ! too small, if we have a really big y(n) + k1/2 list, from which
        ! essentially all determinants get stored into the spawned list
        ! during the death/cloning step..
        ! talk about that with ali!

        character(*), parameter :: this_routine = "walker_death_spawn"

    end subroutine walker_death_spawn

    function attempt_create_realtime(DetCurr, iLutCurr, RealwSign, nJ, iLutnJ, &
                                     prob, HElGen, ic, ex, tParity, walkExcitLevel, part_type, AvSignCurr, &
                                     AvExPerWalker, RDMBiasFacCurr, precond_fac) result(child)

        ! create a specific attempt_create function for the real-time fciqmc
        ! to avoid preprocessor flag jungle..

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! 1 = Real parent particle, 2 = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2, ic), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp), dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp), intent(in) :: AvExPerWalker
        real(dp), intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: precond_fac
        HElement_t(dp), intent(inout) :: HElGen
        character(*), parameter :: this_routine = 'attempt_create_realtime'

        real(dp) :: walkerweight, pSpawn, nSpawn, MatEl, p_spawn_rdmfac, &
                    sepSign, fac_unused
        integer :: extracreate, tgt_cpt, component, iUnused
        integer :: TargetExcitLevel, tmp_ex(2, ic)
        logical :: tRealSpawning
        real(dp) :: rh_imag
        HElement_t(dp) :: rh_used

        unused_var(precond_fac)
        unused_var(AvSignCurr)

        ! This is crucial
        child = 0.0_dp

        ! If each walker does not have exactly one spawning attempt
        ! (if AvMCExcits /= 1.0_dp) then the probability of an excitation
        ! having been chosen, prob, must be altered accordingly.
        prob = prob * AvExPerWalker

        ! In the case of using HPHF, and when tGenMatHEl is on, the matrix
        ! element is calculated at the time of the excitation generation,
        ! and returned in HElGen. In this case, get_spawn_helement simply
        ! returns HElGen, rather than recomputing the matrix element.
        tmp_ex(1, :) = ex(2, :)
        tmp_ex(2, :) = ex(1, :)
        rh_used = get_spawn_helement(nJ, DetCurr, iLutnJ, ilutCurr, ic, tmp_ex, &
                                     tParity, HElGen)

        tRealSpawning = .false.
        if (tAllRealCoeff) then
            tRealSpawning = .true.
        else if (tRealCoeffByExcitLevel) then
            TargetExcitLevel = FindBitExcitLevel(iLutRef, iLutnJ)
            if (TargetExcitLevel <= RealCoeffExcitThresh) &
                tRealSpawning = .true.
        end if

        ! do the whole real-time shabang here, since it also depends on the
        ! fact if it is a pure real-hamiltonian(or atleast we can save some effort)
        ! change in code here for the real-time fciqmc
        ! for now itsonly implemented for pure real hamiltonians
        ! but the -i * H in the diff. eq. makes it kind of all imaginary
        ! for the walker dynamics(except no influence on the conjugation
        ! of H for H_ij -> H_ji*
        ! todo: the case of a paritally imaginary Hamiltonian as input
        ! has also to be considered later on!
        ! here i have then contributions from both real and imaginary
        ! walker populations: when writing the Hamiltonian as: H = H + iJ
        ! and a determinant weight as: n + ic:
        ! n_j + ic_j = -i(H_ij + iJ_ij)^+ (n_i + ic_i)
        ! n_j + ic_j = -i(H_ji - iJ_ji) (n_i + ic_i) =>
        ! leads to =>
        ! n_j = H_ji c_i - J_ji n_i
        ! c_j = -H_ji n_i - J_ji c_i
        ! ... but this should also be with the "normal" complex Hamiltonians
        ! or? talk to ali..
        ! yes this is exactly what is done already.. now implement that
        ! also for the additional -i in the real-time walker dynamics
        if (.not. t_complex_ints .and. .not. t_rotated_time) then
            ! if it is a pure real-hamiltonian there is only spawing from
            ! real to complex walkers and v.v.
            tgt_cpt = rotate_part(part_type)
            walkerweight = sign(1.0_dp, RealwSign(part_type))
            MatEl = real(rh_used, dp)

            ! spawn from real-parent to imaginary child: no sign change
            ! from imaginary to real -> sign change
            if (mod(part_type, 2) == 0) walkerweight = -walkerweight

            nSpawn = -tau * MatEl * walkerweight / prob

            ! n.b. if we ever end up with |walkerweight| /= 1, then this
            !      will need to ffed further through.
            if (tSearchTau .and. (.not. tFillingStochRDMonFly)) &
                call log_spawn_magnitude(ic, ex, matel, prob)

            ! Keep track of the biggest spawn this cycle
            max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)

            if (tRealSpawning) then
                ! Continuous spawning. Add in acceptance probabilities.

                if (tRealSpawnCutoff .and. &
                    abs(nSpawn) < RealSpawnCutoff) then
                    p_spawn_rdmfac = abs(nSpawn) / RealSpawnCutoff
                    nSpawn = RealSpawnCutoff &
                             * stochastic_round(nSpawn / RealSpawnCutoff)
                else
                    p_spawn_rdmfac = 1.0_dp !The acceptance probability of some kind of child was equal to 1
                end if
            else
                if (abs(nSpawn) >= 1) then
                    p_spawn_rdmfac = 1.0_dp !We were certain to create a child here.
                    ! This is the special case whereby if P_spawn(j | i) > 1,
                    ! then we will definitely spawn from i->j.
                    ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                    ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                    ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                else
                    p_spawn_rdmfac = abs(nSpawn)
                end if

                ! How many children should we spawn?

                ! And round this to an integer in the usual way
                ! HACK: To use the same number of random numbers for the tests.
                nSpawn = real(stochastic_round(nSpawn), dp)

            end if
            ! And create the parcticles
            child(tgt_cpt) = nSpawn

        else

            ! have to loop over the tgt_cpt similar to the complex impl
            ! if the Hamiltonian has real and imaginary components do it
            ! similarily to complex implementation with H <-> J switched
            ! rmneci_setup: adjusted for multirun, fixed complex -> real spawns
            do component = 1, (lenof_sign / inum_runs)
                tgt_cpt = min_part_type(part_type_to_run(part_type)) - 1 + component
                ! keep track of the sign due to the kind of spawn event
                sepSign = 1.0_dp
                ! if (part_type == 2 .and. inum_runs == 1) component = 3 - tgt_cpt !?

                walkerweight = sign(1.0_dp, RealwSign(part_type))
                if (mod(part_type, 2) == 0 .and. component == 1) &
                    sepSign = (-1.0_dp)
                if (t_real_time_fciqmc) then
#ifdef CMPLX_
                    rh_imag = real(aimag(rh_used), dp)
#else
                    rh_imag = 0.0_dp
#endif
                    ! part_type is given as input, for that part_type, the real part of
                    ! the HElement is used if rotation occurs and the imaginary part if not
                    if (mod(component, 2) == mod(part_type, 2)) then
                        ! spawn part_type -> part_type
                        MatEl = -rh_imag * tau_real - real(rh_used, dp) * tau_imag
                    else
                        ! spawn part_type -> rotate_part(part_type)
                        MatEl = real(rh_used, dp) * tau_real - rh_imag * tau_imag
                    end if
                end if
                nSpawn = -sepSign * MatEl * walkerweight / prob

                ! n.b. if we ever end up with |walkerweight| /= 1, then this
                !      will need to ffed further through.
                if (tSearchTau .and. (.not. tFillingStochRDMonFly)) &
                    call log_spawn_magnitude(ic, ex, matel, prob)

                ! Keep track of the biggest spawn this cycle
                max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)

                if (tRealSpawning) then
                    ! Continuous spawning. Add in acceptance probabilities.

                    if (tRealSpawnCutoff .and. &
                        abs(nSpawn) < RealSpawnCutoff) then
                        p_spawn_rdmfac = abs(nSpawn) / RealSpawnCutoff
                        nSpawn = RealSpawnCutoff &
                                 * stochastic_round(nSpawn / RealSpawnCutoff)
                    else
                        p_spawn_rdmfac = 1.0_dp !The acceptance probability of some kind of child was equal to 1
                    end if
                else
                    if (abs(nSpawn) >= 1) then
                        p_spawn_rdmfac = 1.0_dp !We were certain to create a child here.
                        ! This is the special case whereby if P_spawn(j | i) > 1,
                        ! then we will definitely spawn from i->j.
                        ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                        ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                        ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                    else
                        p_spawn_rdmfac = abs(nSpawn)
                    end if

                    ! How many children should we spawn?

                    ! And round this to an integer in the usual way
                    ! HACK: To use the same number of random numbers for the tests.
                    nSpawn = real(stochastic_round(nSpawn), dp)

                end if
                ! And create the parcticles (in the correct run)
                child(tgt_cpt) = nSpawn
            end do
        end if

        if (tFillingStochRDMonFly) then
            if (.not. near_zero(child(part_type))) then
                !Only add in contributions for spawning events within population 1
                !(Otherwise it becomes tricky in annihilation as spawnedparents doesn't tell you which population
                !the event came from at present)
                call calc_rdmbiasfac(p_spawn_rdmfac, prob, realwSign(part_type), RDMBiasFacCurr)
            else
                RDMBiasFacCurr = 0.0_dp
            end if
        else
            ! Not filling the RDM stochastically, bias is zero.
            RDMBiasFacCurr = 0.0_dp
        end if

        ! Avoid compiler warnings
        iUnused = walkExcitLevel

    end function attempt_create_realtime

!------------------------------------------------------------------------------------------!

    subroutine update_elapsed_time()
        implicit none
        integer :: run
        ! bookkeeping of timestats
        ! each iteration step constist of two tau-steps -> factor of 2
        elapsedRealTime = elapsedRealTime + tau_real
        elapsedImagTime = elapsedImagTime + tau_imag
        if (tStaticShift) then
            do run = 1, inum_runs
                if (.not. tSinglePartPhase(run)) DiagSft(run) = asymptoticShift
            end do
        end if

        ! normally, we strictly forbid negative shifts
        if (tOnlyPositiveShift) then
            do run = 1, inum_runs
                if (DiagSft(run) < 0.0_dp) DiagSft(run) = 0.0_dp
            end do
        end if

        ! cheap way of removing all initiators
        if (tInfInit) InitiatorWalkNo = TotWalkers + 1

    end subroutine update_elapsed_time

!------------------------------------------------------------------------------------------!

    subroutine save_current_dets()
        use real_time_data, only: TotPartsStorage
        ! routine to copy the currentDets array and all the associated
        ! pointers an hashtable related quantities to the 2nd temporary
        ! list, from which the first spawn and y(n) + k1/2 addition is done
        ! and the k2 spawing list is created to then use CurrentDets to go to
        ! the next time-step y(n+1) = y(n) + k2
        character(*), parameter :: this_routine = "save_current_dets"

        ! save the WalkVecDets variable, i think thats the only necessary
        ! variable, the pointers don't count
        temp_det_list(:, 1:TotWalkers) = WalkVecDets(:, 1:TotWalkers)

        ! for now also store the pointer, but thats not needed i guess
        temp_det_pointer => temp_det_list

        ! and the freeslot.. although this one gets reinitialized to 0
        ! every iteration or not? yeah it is.. so i only have to reset it
        ! twice in the rt-fciqmc before the y(n) + k2 combination
        ! do that in the reload_current_dets routine!
        ! same with n_determ_states var.

        ! also have to save current number of determinants! (maybe totparts too?)
        temp_totWalkers = TotWalkers

        ! And save the old TotParts value, as this might have changed and iter_data is reset
        ! (some weird scenario in which CalcHashTableStats is called at the end of the
        ! time-step and and then modifies both TotParts and iter_data correctly, but iter_data
        ! is reset at the beginning of the iteration, so TotParts also has to)
        TotPartsStorage = TotParts

    end subroutine save_current_dets

!------------------------------------------------------------------------------------------!

    subroutine reload_current_dets()
        ! routine to reload the saved y(n) CurrentDets array for the final
        ! y(n) + k/2 combination to move to the next time step y(n+1)
        ! have also to think about the death-step and annihilation step
        ! influences.. maybe have to write new death/born routines to split
        ! that from the spawned list creation..
        character(*), parameter :: this_routine = "reload_current_dets"

        ! copy the list
        WalkVecDets(:, 1:temp_totWalkers) = temp_det_list(:, 1:temp_totWalkers)

        ! and point to it
        CurrentDets => WalkVecDets

        ! also have to reset the number of determinants to the original
        ! value?!
        TotWalkers = temp_totWalkers

        ! and the hash
        ! cant just copy the hash-table like that have to associate all
        ! entries correctly
        ! here i 'just' have to reassign the original HashIndex with the
        ! stored CurrentDets list!
        call clear_hash_table(HashIndex)
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, &
                                int(TotWalkers, sizeof_int), .true.)

        ! for correct load Balancin i also have to reset the the freeslot var.
        ! here i have to reset the freeslot values to the values after the
        ! first spawn loop, as the empty entries get determined there!
        ! and i only want to do an additional annihilation step to combine
        ! y(n) + k2
        ! also reload the positions of empty slots in the ensemble
        iStartFreeSlot = 1
        iEndFreeSlot = temp_iendfreeslot
        FreeSlot = temp_freeslot

        call reset_tot_parts()

    end subroutine reload_current_dets

!------------------------------------------------------------------------------------------!

    subroutine reset_spawned_list()
        ! also need a routine to reset the spawned lists before the second
        ! spawning step for the 2nd order RK method in the rt-fciqmc
        character(*), parameter :: this_routine = "reset_spawned_list"

        ! Reset positions to spawn into in the spawning array.
        ValidSpawnedList = InitialSpawnedSlots

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        ! also reset the diagonal specific valid spawn list.. i think i
        ! can just reuse the InitialSpawnedSlots also
        valid_diag_spawns = 1

    end subroutine reset_spawned_list

!------------------------------------------------------------------------------------------!

    subroutine setup_temp_det_list()
        ! setup the second list to temporaly store the list of determinants
        ! necessary in the real-time fciqmc list
        ! determine the necessary size from the already setup CurrentDets
        character(*), parameter :: this_routine = "setup_temp_det_list"
        integer :: ierr, tmp_siz1, tmp_siz2, i, spawn_ht_mem

        tmp_siz1 = size(WalkVecDets, dim=1)
        tmp_siz2 = size(WalkVecDets, dim=2)

        ! allocate the array
        allocate(temp_det_list(0:tmp_siz1 - 1, tmp_siz2), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error in allocation")

        ! and init it
        temp_det_list(0:tmp_siz1 - 1, 1:tmp_siz2) = 0

        ! and point to it
        temp_det_pointer => temp_det_list

        ! and also allocate the hash-table
        tmp_siz1 = size(HashIndex)
        allocate(temp_det_hash(tmp_siz1), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error in allocation")

        ! and initialize it to 0
        do i = 1, tmp_siz1
            temp_det_hash(i)%ind = 0
        end do

        ! also use the spawn_ht hash table, so also allocate it here!

        ! Allocate the hash table to the spawning array.
        ! The number of MB of memory required to allocate spawn_ht.
        ! Each node requires 16 bytes.
        nhashes_spawn = int(0.8_dp * real(MaxSpawned, dp))
        spawn_ht_mem = nhashes_spawn * 16 / 1000000
        write(6, '(a78,'//int_fmt(spawn_ht_mem, 1)//')') "About to allocate hash table to the spawning array. &
                                       &Memory required (MB):", spawn_ht_mem
        write(6, '(a13)', advance='no') "Allocating..."; call neci_flush(6)
        allocate(spawn_ht(nhashes_spawn), stat=ierr)
        if (ierr /= 0) then
            write(6, '(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(this_routine, "Error allocating spawn_ht array.")
        else
            write(6, '(1x,a5)') "Done."
            write(6, '(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                              &increase as further nodes are added."
        end if

        call init_hash_table(spawn_ht)

    end subroutine setup_temp_det_list

!------------------------------------------------------------------------------------------!

    subroutine update_gf_overlap()
        ! subroutine to calculate the overlap of the current y(t) = a_j(a^+_j)(t)y(0)>
        ! time evolved wavefunction to the saved <y(0)|a^+_i(a_i)
        use timing_neci, only: timer, get_total_time
        implicit none
        integer :: idet, nI(nel), det_ind, hash_val, runA, runB, iGf
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        complex(dp) :: overlap(normsize)
        logical :: tDetFound
        real(dp) :: gf_time

        call set_timer(calc_gf_time)

        do iGf = 1, gf_count
            overlap = cmplx(0.0_dp, 0.0_dp, dp)
            do idet = 1, overlap_states(iGf)%nDets

                call extract_sign(overlap_states(iGf)%dets(:, idet), real_sign_1)

                if (IsUnoccDet(real_sign_1)) cycle

                call decode_bit_det(nI, overlap_states(iGf)%dets(:, idet))

                ! search for the hash table associated with the time evolved
                ! wavefunction -> is this already initialized correctly?
                call hash_table_lookup(nI, overlap_states(iGf)%dets(:, idet), nifd, &
                                       HashIndex, CurrentDets, det_ind, hash_val, tDetFound)

                if (tDetFound) then
                    ! both real and imaginary part of the time-evolved wf are required
                    call extract_sign(CurrentDets(:, det_ind), real_sign_2)

                    do runA = 1, inum_runs
                        do runB = 1, inum_runs
                            ! overlap is now treated as complex type
                            ! this only works for the complex code
                            overlap(overlap_index(runA, runB)) = overlap(overlap_index(runA, runB)) &
                                                                 + conjg(cmplx(real_sign_1(min_part_type(runA)), &
                                                                               real_sign_1(max_part_type(runA)), dp)) &
                                                             * cmplx(real_sign_2(min_part_type(runB)), real_sign_2(max_part_type(runB)), dp)
                        end do
                    end do
                end if
            end do

            ! rmneci_setup: the overlap has to be reduced as each proc
            ! only computes its own part
            call MPIReduce(overlap, MPI_SUM, gf_overlap(:, iGf))
        end do

        call halt_timer(calc_gf_time)
        gf_time = get_total_time(calc_gf_time)

    end subroutine update_gf_overlap

!------------------------------------------------------------------------------------------!

    subroutine normalize_gf_overlap(overlapList, avReal, avImag)
        implicit none
        complex(dp), intent(out) :: overlapList(normsize, gf_count)
        real(dp), intent(out) :: avReal(gf_count), avImag(gf_count)
        integer :: i, j
        complex(dp), allocatable :: overlap_buf(:)

        allocate(overlap_buf(gf_count))
        call update_gf_overlap()

        do j = 1, gf_count
            do i = 1, normsize
                overlapList(i, j) = gf_overlap(i, j) / dyn_norm_red(i, j) * &
                                    exp(shift_damping(((i - 1) / inum_runs + 1)))
            end do

            !normalize the greens function
            overlap_buf(j) = sum(gf_overlap(:, j)) / sum(dyn_norm_red(:, j)) * &
                             sum(exp(shift_damping)) / inum_runs

            avReal(j) = real(overlap_buf(j))
            avImag(j) = aimag(overlap_buf(j))
        end do

        deallocate(overlap_buf)
    end subroutine normalize_gf_overlap

!------------------------------------------------------------------------------------------!

    function calc_norm(dets, num_dets) result(cd_norm)
        ! the first dimension of dets has to be niftot
        ! function to calculate the norm of a state and
        ! the overlap between replicas(general function)
        complex(dp) :: cd_norm(normsize)
        integer(dp) :: dets(0:, 1:)
        integer, intent(in) :: num_dets
        character(*), parameter :: this_routine = "calc_perturbed_norm"

        integer :: idet, run, targetRun
        real(dp) :: tmp_sign(lenof_sign)

        cd_norm = 0.0_dp

        do idet = 1, num_dets

            call extract_sign(dets(:, idet), tmp_sign)
            do run = 1, inum_runs
                ! we calculate the overlap between any two replicas, including the norm
                ! of each individually
                do targetRun = 1, run
! this only works for complex builds, it would not even compile else
                    cd_norm(overlap_index(run, targetRun)) = cd_norm(overlap_index(run, targetRun)) &
                                                            + conjg(cmplx(tmp_sign(min_part_type(run)), tmp_sign(max_part_type(run)), dp)) &
                                                             * cmplx(tmp_sign(min_part_type(targetRun)), tmp_sign( &
                                                                     max_part_type(targetRun)), dp)
                end do
            end do
            do run = 1, inum_runs
                do targetRun = run + 1, inum_runs
                    cd_norm(overlap_index(run, targetRun)) = conjg(cd_norm(overlap_index(targetRun, run)))
                end do
            end do
        end do

    end function calc_norm

!------------------------------------------------------------------------------------------!

    subroutine makePopSnapshot(i)
        use real_time_data, only: popSnapshot, snapshotOrbs, numSnapshotOrbs

        implicit none
        integer, intent(in) :: i
        integer :: iOrb, nI(nel), iEl, part
        real(dp) :: avPop, tmpSign(lenof_sign)

        call decode_bit_det(nI, CurrentDets(:, i))
        do iOrb = 1, numSnapshotOrbs
            do iEl = 1, nel
                if (nI(iEl) == snapshotOrbs(iOrb)) then
                    avPop = 0
                    call extract_sign(CurrentDets(:, i), tmpSign)
                    do part = 1, lenof_sign
                        avPop = avPop + abs(tmpSign(part))
                    end do
                    avPop = avPop / inum_runs
                    popSnapshot(iOrb) = popSnapshot(iOrb) + avPop
                end if
            end do
        end do

    end subroutine makePopSnapshot

!------------------------------------------------------------------------------------------!

    subroutine real_time_determ_projection()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use FciMCData, only: SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV

        integer :: i, j, ierr, run

        associate(rep => cs_replicas(core_run))
            call MPIBarrier(ierr)

            call set_timer(SemiStoch_Comms_Time)

            call MPIAllGatherV(rep%partial_determ_vecs, rep%full_determ_vecs, &
                               rep%determ_sizes, rep%determ_displs)

            call halt_timer(SemiStoch_Comms_Time)

            call set_timer(SemiStoch_Multiply_Time)

#ifdef CMPLX_

            if (rep%determ_sizes(iProcIndex) >= 1) then

                ! Perform the multiplication.

                rep%partial_determ_vecs = 0.0_dp

                do i = 1, rep%determ_sizes(iProcIndex)
                    do j = 1, rep%sparse_core_ham(i)%num_elements
                        do run = 1, inum_runs
                            ! real part of the 'spawn'
                            rep%partial_determ_vecs(min_part_type(run), i) = &
                                rep%partial_determ_vecs(min_part_type(run), i) + &
                                tau_real * Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                min_part_type(run), rep%sparse_core_ham(i)%positions(j)) + &
                                tau_real * Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                max_part_type(run), rep%sparse_core_ham(i)%positions(j)) + &
                                tau_imag * Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                min_part_type(run), rep%sparse_core_ham(i)%positions(j)) - &
                                tau_imag * Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                max_part_type(run), rep%sparse_core_ham(i)%positions(j))

                            ! imaginary part
                            rep%partial_determ_vecs(max_part_type(run), i) = &
                                rep%partial_determ_vecs(max_part_type(run), i) - &
                                tau_real * Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                min_part_type(run), rep%sparse_core_ham(i)%positions(j)) + &
                                tau_real * Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                max_part_type(run), rep%sparse_core_ham(i)%positions(j)) + &
                                tau_imag * Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                max_part_type(run), rep%sparse_core_ham(i)%positions(j)) + &
                                tau_imag * Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs( &
                                min_part_type(run), rep%sparse_core_ham(i)%positions(j))

                        end do
                    end do
                end do

                ! Now add shift*rep%full_determ_vecs to account for the shift, not stored in
                ! rep%sparse_core_ham.
                do i = 1, rep%determ_sizes(iProcIndex)
                    do run = 1, inum_runs
                        ! real part
                        rep%partial_determ_vecs(min_part_type(run), i) = &
                            rep%partial_determ_vecs(min_part_type(run), i) + &
                            (Hii - gs_energy(run)) * rep%full_determ_vecs(max_part_type(run), i + &
                                                                          rep%determ_displs(iProcIndex)) * &
                            tau_real + (tau_imag * (Hii - gs_energy(run) - DiagSft(run)) + &
                                        tau_real * real_time_info%damping) * rep%full_determ_vecs( &
                            min_part_type(run), i + rep%determ_displs(iProcIndex))

                        ! imaginary part
                        rep%partial_determ_vecs(max_part_type(run), i) = &
                            rep%partial_determ_vecs(max_part_type(run), i) + (tau_imag * &
                                                                              (Hii - gs_energy(run) - DiagSft(run)) + tau_real &
                                                                              * real_time_info%damping) * rep%full_determ_vecs( &
                            max_part_type(run), i + rep%determ_displs(iProcIndex)) - tau_real &
                            * (Hii - gs_energy(run)) * rep%full_determ_vecs(min_part_type(run), i &
                                                                            + rep%determ_displs(iProcIndex))
                    end do
                end do
            end if

#endif
            call halt_timer(SemiStoch_Multiply_Time)
        end associate

    end subroutine real_time_determ_projection

!------------------------------------------------------------------------------------------!

    subroutine refresh_semistochastic_space()
        use CalcData, only: ss_space_in
        use semi_stoch_gen, only: init_semi_stochastic
        use semi_stoch_procs, only: end_semistoch
        implicit none
        logical :: tStartedFromCoreGround
        ! as the important determinants might change during time evolution, this
        ! resets the semistochastic space taking the current population to get a new one
        call end_semistoch()
        ! the flag_deterministic flag has to be cleared from all determinants as it is
        ! assumed that no states have that flag when init_semi_stochastic starts
        call reset_core_space()
        call init_semi_stochastic(ss_space_in, tStartedFromCoreGround)

    end subroutine refresh_semistochastic_space

!------------------------------------------------------------------------------------------!

    subroutine create_perturbed_ground()
        ! routine to create from the already read in popsfile info in
        ! popsfile_dets the left hand <y(0)| by applying the corresponding
        ! creation or annihilation operator
        implicit none
        character(*), parameter :: this_routine = "create_perturbed_ground"
        integer :: tmp_totwalkers, totwalkers_backup, TotWalkers_orig_max
        integer :: ierr, i, totNOccDets, iProc, nPertRefs
        integer(n_int), allocatable :: perturbed_buf(:, :)
        logical :: t_use_perturbed_buf

        if (tReadPops .and. .not. tNewOverlap) then
            tmp_totwalkers = int(TotWalkers_orig)
        else
            tmp_totwalkers = int(TotWalkers)
        end if

        write(iout, *) "Creating the wavefunction to projected on!"
        write(iout, *) "Initial number of walkers: ", tmp_totwalkers

        if (allocated(overlap_pert)) then
            call MPISumAll(tmp_totwalkers, TotWalkers_orig_max)
        else
            TotWalkers_orig_max = MaxWalkersPart
        end if
        t_use_perturbed_buf = allocated(overlap_pert) .and. tNewOverlap

        if (.not. allGfs == 0) call setup_pert_array(allGfs)

        allocate(overlap_states(gf_count), stat=ierr)
        if (t_use_perturbed_buf) &
            allocate(perturbed_buf(0:niftot, TotWalkers_orig_max), stat=ierr)
        write(iout, *) "Read-in dets", TotWalkers_orig
        do i = 1, gf_count
            totwalkers_backup = tmp_totwalkers
            if (t_use_perturbed_buf) then
                if (.not. t_kspace_operators) then
                    if (tReadPops) then
                        ! if the perturbation is to be created from the read-in population
                        ! explicitly
                        perturbed_buf = 0.0_dp
                        call apply_perturbation(overlap_pert(i), tmp_totwalkers, popsfile_dets, &
                                                perturbed_buf)

                        ! The HF-Overlap option makes us use a reference for projection
                        ! instead of the full wavefunction
                        if (tHFOverlap) call create_perturbed_ref(perturbed_buf, TotWalkers_orig_max)
                    else
                        ! else, perturb the current population (this is for starting from
                        ! a defined determinant)
                        call apply_perturbation(overlap_pert(i), tmp_totwalkers, CurrentDets, &
                                                perturbed_buf)
                    end if
                else
                    ! a less useful feature for fourier-transforming the perturbation
                    if (gf_count > 1) call stop_all("create_perturbed_ground", &
                                                    "Unable to use momentum operators for multiple correlation functions")
                    if (tReadPops) then
                        perturbed_buf = 0.0_dp
                        call apply_perturbation_array(overlap_pert, tmp_totwalkers, popsfile_dets, &
                                                      perturbed_buf, phase_factors)
                    else
                        call apply_perturbation_array(overlap_pert, tmp_totwalkers, CurrentDets, &
                                                      perturbed_buf, phase_factors)
                    end if
                end if

                call write_overlap_state_serial(perturbed_buf, TotWalkers_orig_max, i)
            else
                write(6, *) "Generated overlap state"
                call write_overlap_state_serial(CurrentDets, int(TotWalkers), i)
                write(6, *) "Written overlap state to array"
            end if
            call MPISumAll(overlap_states(i)%nDets, totNOccDets)
            if (totNOccDets == 0) then
                if (gf_count == 1) then
                    call stop_all('create_perturbed_ground', 'No walkers survived perturbation')
                else
                    write(6, *) "WARNING, EMPTY PERTURBED STATE WITH INDEX", i
                end if
            end if
            tmp_totwalkers = totwalkers_backup
        end do

        if (allocated(overlap_states)) write(iout, *) &
            "Determinants remaining in perturbed ground state:", overlap_states(1)%nDets
        if (allocated(perturbed_buf)) deallocate(perturbed_buf, stat=ierr)

    end subroutine create_perturbed_ground

!------------------------------------------------------------------------------------------!

    subroutine create_perturbed_ref(perturbed_buf, tmp_totwalkers)
        implicit none
        integer(n_int), intent(out) :: perturbed_buf(:, :)
        integer, intent(out) :: tmp_totwalkers

        integer :: nPertRefs
        integer(n_int) :: tmpRef(0:NIfTot, 1)
        character(*), parameter :: t_r = "create_perturbed_reference"

        ! create a perturbed reference determinant as overlap state
        nPertRefs = 0
        ! i.e. take the most populated determinant (over all procs)
        call generate_space_most_populated(1, .false., 1, tmpRef, nPertRefs, &
                                           GLOBAL_RUN, perturbed_buf, int(tmp_totwalkers, n_int))

        ! from now on, we treat the perturbed_buf as 1-sized
        perturbed_buf = 0
        perturbed_buf(:, 1) = tmpRef(:, 1)
        tmp_totwalkers = nPertRefs
    end subroutine create_perturbed_ref

!------------------------------------------------------------------------------------------!

    subroutine check_update_growth(iter_data, message)
        use real_time_data, only: TotPartsStorage
        implicit none
        character(len=*), intent(in) :: message
        type(fcimc_iter_data), intent(in) :: iter_data
        real(dp) :: growth(lenof_sign), growth_tot(lenof_sign)
        real(dp) :: allWalkers(lenof_sign), allWalkersOld(lenof_sign)
        real(dp) :: allDied(lenof_sign), allBorn(lenof_sign), allAnnihil(lenof_sign), &
                    allAbrt(lenof_sign), allRmv(lenof_sign)
        growth = iter_data%nborn &
                 - iter_data%ndied - iter_data%nannihil &
                 - iter_data%naborted - iter_data%nremoved

        call MPISumAll(growth, growth_tot)
        call MPISumAll(TotParts, allWalkers)
        call MPIsumAll(TotPartsStorage, allWalkersOld)

        call MPISumAll(iter_data%nborn, allBorn)
        call MPISumAll(iter_data%ndied, allDied)
        call MPISumAll(iter_data%nannihil, allAnnihil)
        call MPISumAll(iter_data%naborted, allAbrt)
        call MPISumAll(iter_data%nremoved, allRmv)
        TotPartsStorage = TotParts
        if ((iProcIndex == root) .and. .not. tSpinProject .and. &
            any(abs(growth_tot - (allWalkers - allWalkersOld)) > 1.0e-4_dp)) then
            write(iout, *) message
            write(iout, *) "update_growth: ", growth_tot
            write(iout, *) "AllTotParts: ", allWalkers
            write(iout, *) "AllTotPartsOld: ", allWalkersOld
            write(6, *) "nborn", allBorn
            write(6, *) "ndied", allDied
            write(6, *) "nannihil", allAnnihil
            write(6, *) "naborted", allAbrt
            write(6, *) "nremoved", allRmv

            call stop_all("check_update_growth", &
                          "Assertation failed: all(iter_data_fciqmc%update_growth_tot.eq.AllTotParts_1-AllTotPartsOld_1)")
        end if

    end subroutine check_update_growth

!------------------------------------------------------------------------------------------!

    subroutine update_shift_damping()
        implicit none
        ! sign convention for imaginary and real time differ
        shift_damping = shift_damping + tau_imag * DiagSft
        if (tDynamicDamping) shift_damping = shift_damping - tau_real &
                                             * real_time_info%damping

    end subroutine update_shift_damping

!------------------------------------------------------------------------------------------!

    subroutine setup_pert_array(ctrn_index)
        implicit none
        integer, intent(in) :: ctrn_index
        integer :: i

        gf_count = nBasis
        allocate(overlap_pert(nBasis))
        do i = 1, nBasis
            if (ctrn_index == 2) then
                overlap_pert(i)%ncreate = 1
                allocate(overlap_pert(i)%crtn_orbs(1))
                overlap_pert(i)%crtn_orbs(1) = i
                call init_perturbation_creation(overlap_pert(i))
            else
                overlap_pert(i)%nannihilate = 1
                allocate(overlap_pert(i)%ann_orbs(1))
                overlap_pert(i)%ann_orbs(1) = i
                call init_perturbation_annihilation(overlap_pert(i))
            end if
        end do

    end subroutine setup_pert_array

!------------------------------------------------------------------------------------------!

    subroutine merge_spawn(nspawn, prefactor)
        use real_time_data, only: nspawnMax
        implicit none
        integer :: nspawn
        real(dp) :: prefactor
        ! truncate the number of spawns from a single determinant
        ! for now, use as a threshold a multiple of the average population
        ! for a full SpawnVec
        if (nspawn > nspawnMax) then
            prefactor = nspawn / real(nspawnMax, dp)
            nspawn = nspawnMax
        else
            prefactor = 1.0_dp
        end if
        ! the prefactor is used to unbias therefor
    end subroutine merge_spawn

!------------------------------------------------------------------------------------------!

    subroutine trunc_shift()
        implicit none
        integer :: run

        do run = 1, inum_runs
            ! remember that shiftLimit is the absolute value, but we are only
            ! interested in shifts that are too small
            if (DiagSft(run) < -1.0_dp * shiftLimit) then
                numCycShiftExcess(run) = numCycShiftExcess(run) + 1
                if (numCycShiftExcess(run) > 20) DiagSft(run) = -0.95_dp * shiftLimit
            else
                numCycShiftExcess(run) = 0
            end if
        end do
    end subroutine trunc_shift

!------------------------------------------------------------------------------------------!

    subroutine adjust_decay_channels()
        use FciMCData, only: AllTotParts
        use CalcData, only: InitWalkers
        use Parallel_neci, only: nProcessors
        use real_time_data, only: alphaDamping, etaDamping, tStartVariation, rotThresh, &
                                  numSnapShotOrbs, tLowerThreshold
        implicit none
        real(dp) :: allWalkersOld(lenof_sign), walkersOld(lenof_sign)
        real(dp) :: deltaAlpha, deltaEta

        ! once the walker number exceeds the total walkers set in the input, start
        ! adjusting the damping and the real/imag timestep ratio
        if (.not. tStartVariation) then
            if (sum(AllTotParts) / inum_runs > rotThresh) tStartVariation = .true.
            if (tLowerThreshold) then
                if (tStartVariation) then
                    tStartVariation = .false.
                else
                    tStartVariation = .true.
                end if
            end if
        end if
        ! once started, we have to do so forever, else we might kill all walkers

        if (tStartVariation) then
            call MPIReduce(TotPartsLastAlpha, MPI_Sum, allWalkersOld)
            ! as AllTotParts is only reduced for shift computation, we need to
            ! do it here manually (TotParts is recomputed every iteration as
            ! a part of the RK-Scheme)
            call MPIReduce(TotParts, MPI_Sum, AllTotParts)
            if (iProcIndex == root) then
                ! compare the walker number the last time the angle was adjusted to
                ! the walker number now
                if (tDynamicAlpha) then
                    deltaAlpha = alphaDamping * atan(sum(AllTotParts) / real(sum(allWalkersOld), dp) - 1)
                    real_time_info%time_angle = real_time_info%time_angle + deltaAlpha
                end if
                ! if the damping is also to be adjusted on the fly, do so here
                if (tDynamicDamping) then
                    deltaEta = etaDamping * log(sum(AllTotParts) / real(sum(allWalkersOld), dp)) / &
                               (tau_real * stepsAlpha)
                    real_time_info%damping = real_time_info%damping - deltaEta
                end if
            end if
            ! communicate the updated quantities
            if (tDynamicAlpha) then
                call MPIBCast(real_time_info%time_angle)
                ! Store the value of alpha in a log, overwriting an old value
                alphaLog(alphaLogPos) = real_time_info%time_angle
                alphaLogPos = alphaLogPos + 1
                ! set back the position if exceeding the end
                if (alphaLogPos > alphaLogSize) alphaLogPos = 1
                !possibly change the stepsize
                call adjust_stepsAlpha()
            end if
            if (tDynamicDamping) call MPIBCast(real_time_info%damping)
        end if
        TotPartsLastAlpha = TotParts

    end subroutine adjust_decay_channels

!------------------------------------------------------------------------------------------!

    subroutine adjust_stepsAlpha()
        implicit none
        integer :: newStep, maxAlpha, minAlpha
        real(dp), parameter :: ratioThreshold = 0.01

        maxAlpha = int(maxval(alphaLog))
        minAlpha = int(minval(alphaLog))
    end subroutine adjust_stepsAlpha

!------------------------------------------------------------------------------------------!

    subroutine reset_tot_parts()
        ! if the second RK step is to be compared, the reference has to be reset
        ! -> recount the TotParts from the restored data
        use real_time_data, only: TotPartsStorage
        implicit none
        TotParts = get_tot_parts()
        TotPartsStorage = TotParts
    end subroutine reset_tot_parts

!------------------------------------------------------------------------------------------!

    function get_tot_parts() result(allWalkersSummed)
        ! if the second RK step is to be compared, the reference has to be reset
        ! -> recount the TotParts from the restored data
        implicit none
        integer(int64) :: i
        real(dp) :: CurrentSign(lenof_sign)
        real(dp) :: allWalkersSummed(lenof_sign)
        allWalkersSummed = 0.0_dp
        do i = 1, TotWalkers
            call extract_sign(CurrentDets(:, i), CurrentSign)
            allWalkersSummed = allWalkersSummed + abs(CurrentSign)
        end do
        if (tStabilizerShift) then
            do i = 1, inum_runs
                TotPartsPeak(i) = max(allWalkersSummed(min_part_type(i)) + &
                                      allWalkersSummed(max_part_type(i)), TotPartsPeak(i))
            end do
        end if
    end function get_tot_parts

!------------------------------------------------------------------------------------------!

    subroutine update_peak_walker_number()
        use FciMCData, only: AllTotParts
        implicit none
        integer :: run

        do run = 1, inum_runs
            TotPartsPeak(run) = max(AllTotParts(min_part_type(run)) + AllTotParts(max_part_type(run)) &
                                    , TotPartsPeak(run))
        end do
    end subroutine update_peak_walker_number

!------------------------------------------------------------------------------------------!

    subroutine clean_overlap_states()
        implicit none
        integer :: i, ierr
        if (allocated(overlap_states)) then
            do i = 1, gf_count
                if (allocated(overlap_states(i)%dets)) deallocate(overlap_states(i)%dets, stat=ierr)
            end do
            deallocate(overlap_states, stat=ierr)
        end if
    end subroutine clean_overlap_states

!------------------------------------------------------------------------------------------!

    subroutine logTimeCurve
        implicit none
        character(*), parameter :: this_routine = "read_in_trajectory"

        if (iProcIndex == root) write(iunitCycLog, *) tau, real_time_info%time_angle

    end subroutine logTimeCurve

!------------------------------------------------------------------------------------------!

    subroutine openTauContourFile
        implicit none

        ! Only root writes out the trajectory
        if (iProcIndex == root) then
            iunitCycLog = get_free_unit()
            call get_unique_filename('tauContour', .true., .true., 0, trajFile)

            open(iunitCycLog, file=trajFile, status='new')
        end if

    end subroutine openTauContourFile

!------------------------------------------------------------------------------------------!

    subroutine closeTauContourFile
        implicit none

        if (iProcIndex == root) then
            close(iunitCycLog)
        end if
    end subroutine closeTauContourFile

!------------------------------------------------------------------------------------------!

    subroutine get_current_alpha_from_cache
        implicit none

        ! the logging and reading are done before iter is updated
        real_time_info%time_angle = alphaCache(iter + 1)
        tau = tauCache(iter + 1)
    end subroutine get_current_alpha_from_cache

!------------------------------------------------------------------------------------------!

    subroutine expand_corespace_buf(buffer, buffer_size)
        use real_time_data, only: ssht, wn_threshold
        use real_time_aux, only: add_semistochastic_state
        use FciMCData, only: AllTotParts
        implicit none
        integer(n_int), intent(inout) :: buffer(0:, 1:)
        integer, intent(inout) :: buffer_size
        integer(int64) :: i
        real(dp) :: sgn(lenof_sign)

        do i = 1, TotWalkers
            call extract_sign(CurrentDets(:, i), sgn)
            if (sum(abs(sgn)) / inum_runs > wn_threshold * sum(abs(AllTotParts)) / inum_runs) then
                call add_semistochastic_state(buffer, buffer_size, ssht, CurrentDets(:, i))
            end if
        end do

    end subroutine expand_corespace_buf

!------------------------------------------------------------------------------------------!

    subroutine get_corespace_from_buf(buffer, buffer_size)
        use core_space_util, only: cs_replicas
        use CalcData, only: ss_space_in
        use semi_stoch_gen, only: generate_space_most_populated
        use semi_stoch_procs, only: store_whole_core_space, write_core_space
        implicit none
        integer, intent(in) :: buffer_size
        integer(n_int), intent(in), pointer :: buffer(:, :)
        integer :: space_size, ierr, i
        integer(MPIArg) :: mpi_buf

        associate(rep => cs_replicas(core_run))
            ! Get the most populated from those determinants that exceeded the threshold once
            space_size = 0
            allocate(rep%determ_sizes(0:nProcessors - 1))
            call generate_space_most_populated(ss_space_in%npops, ss_space_in%tApproxSpace, &
                                               ss_space_in%nApproxSpace, SpawnedParts, space_size, GLOBAL_RUN, buffer, &
                                               int(buffer_size, n_int))

            ! Then, communicate the number of states per core
            mpi_buf = int(space_size, MPIArg)
            call MPIAllGather(mpi_buf, rep%determ_sizes, ierr)
            rep%determ_space_size = sum(rep%determ_sizes)

            ! Prepare the store_whole_core_space communication routine
            allocate(rep%determ_displs(0:nProcessors - 1))

            ! Communciate the newly built corespace and output it
            call store_whole_core_space(rep)

            call write_core_space(rep)
        end associate

    end subroutine get_corespace_from_buf

end module real_time_procs
