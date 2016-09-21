#include "macros.h"

! module containing useful functions and subroutines needed in the real-time
! implementation of the FCIQMC algotrithm

module real_time_procs
    use hash, only: hash_table_lookup, init_hash_table, clear_hash_table, &
                    add_hash_table_entry, fill_in_hash_table
    use SystemData, only: nel
    use real_time_data, only: gf_overlap, TotWalkers_orig, TotWalkers_pert, &
                              t_complex_ints, real_time_info, temp_freeslot, & 
                              temp_iendfreeslot, temp_det_list, temp_det_pointer, &
                              temp_det_hash, temp_totWalkers, pert_norm, &
                              valid_diag_spawn_list, DiagParts, n_diag_spawned, &
                              DiagParts2, NoDied_1, NoBorn_1, SumWalkersCyc_1, &
                              current_overlap
    use kp_fciqmc_data_mod, only: perturbed_ground, overlap_pert
    use constants, only: dp, lenof_sign, int64, n_int, EPS, iout, null_part, &
                         sizeof_int, MPIArg
    use bit_reps, only: decode_bit_det, test_flag, flag_initiator, encode_sign, &
                        set_flag, encode_bit_rep, extract_bit_rep, &
                        flag_has_been_initiator, flag_deterministic, encode_part_sign, &
                        nullify_ilut_part

                        
    use bit_rep_data, only: extract_sign, nifdbo, niftot
    use FciMCData, only: CurrentDets, HashIndex, popsfile_dets, MaxWalkersPart, &
                         WalkVecDets, freeslot, spawn_ht, nhashes_spawn, MaxSpawned, &
                         iStartFreeSlot, iEndFreeSlot, ValidSpawnedList, &
                         InitialSpawnedSlots, iLutRef, inum_runs, max_cyc_spawn, &
                         tSearchTau, tFillingStochRDMonFly, fcimc_iter_data, &
                         NoAddedInitiators, SpawnedParts, acceptances, TotWalkers, &
                         nWalkerHashes, iter, fcimc_excit_gen_store, NoDied, &
                         NoBorn, NoAborted, NoRemoved, HolesInList, TotParts
    use perturbations, only: apply_perturbation
    use util_mod, only: int_fmt
    use CalcData, only: AvMCExcits, tAllRealCoeff, tRealCoeffByExcitLevel, &
                        tRealSpawnCutoff, RealSpawnCutoff, tau, RealCoeffExcitThresh, &
                        DiagSft, tTruncInitiator, OccupiedThresh
    use DetBitOps, only: FindBitExcitLevel
    use procedure_pointers, only: get_spawn_helement
    use util_mod, only: stochastic_round
    use tau_search, only: log_spawn_magnitude
    use rdm_general, only: calc_rdmbiasfac
    use global_det_data, only: global_determinant_data
    use rdm_filling, only: det_removed_fill_diag_rdm, check_fillRDM_DiDj
! RT_M_Merge: Disabled rdms
!    use rdm_data, only: nrdms, rdms
    use hash, only: remove_hash_table_entry
    use dSFMT_interface, only: genrand_real2_dSFMT
    use load_balance_calcnodes, only: DetermineDetNode
    use ParallelHelper, only: nNodes, bNodeRoot, ProcNode, NodeRoots, MPIBarrier, iProcIndex
    use Parallel_neci, only: nProcessors, MPIAlltoAll, MPIAlltoAllv
    use LoggingData, only: tNoNewRDMContrib
    use AnnihilationMod, only: test_abort_spawn
    use load_balance, only: AddNewHashDet, CalcHashTableStats, adjust_load_balance, test_hash_table

    implicit none

contains

    subroutine DirectAnnihilation_diag(TotWalkersNew, iter_data, tSingleProc)
        ! new direct annihilation routine to mimick the diagonal death step
        ! in the y(n) + k2 combination between reloaded CurrentDets and the 
        ! DiagParts list
        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        logical, intent(in) :: tSingleProc
        character(*), parameter :: this_routine = "DirectAnnihilation_diag"
        integer :: numSpawns

        ! As only diagonal events are considered here, no communication
        ! is required
        ! Also, this eliminates the need for compression as all dets are
        ! already stored contigously and in the right orde
        ! (no annihilation inside DiagParts can occur)

        numSpawns = valid_diag_spawn_list(iProcIndex)
        call AnnihilateDiagParts(numSpawns, TotWalkersNew, iter_data)

        ! also should update the hashtable stats, specific for this diagonal 
        ! spawning event, but the original one should work also for this 
        ! since it only takes CurrentDets into account! 
        
        call CalcHashTableStats(TotWalkersNew,iter_data)
        
        ! this should be it..
        
    end subroutine DirectAnnihilation_diag

    subroutine AnnihilateDiagParts(ValidSpawned, TotWalkersNew, iter_data)
        ! this is the new "annihilation" routine which mimics the actual 
        ! diagonal death-step, between the reloaded CurrentDets y(n) and the 
        ! diagonal parts in k2, DiagParts
        integer, intent(inout) :: ValidSpawned, TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data 
        character(*), parameter :: this_routine = "AnnihilateDiagParts"

        integer :: PartInd, i, j, PartIndex
        real(dp), dimension(lenof_sign) :: CurrentSign, SpawnedSign, SignTemp
        real(dp), dimension(lenof_sign) :: TempCurrentSign, SignProd
        real(dp) :: pRemove, r
        integer :: ExcitLevel, DetHash, nJ(nel)
        logical :: tSuccess, tSuc, tPrevOcc, tDetermState
        integer :: run       
        integer :: nI(nel)
        
        ! rewrite the original Annihilation routine to fit the new 
        ! requirements here

        ! this routine updated the NoDied variables etc.. for the 2nd RK step
        ! so i dont think i need to change much here
        ! Only node roots to do this.
        if (.not. bNodeRoot) return

!         call set_timer(AnnMain_time, 30)

!         if (tHistSpawn) HistMinInd2(1:nEl) = FCIDetIndex(1:nEl)

!         call set_timer(BinSearch_time,45)

        do i = InitialSpawnedSlots(iProcIndex), ValidSpawned

           call decode_bit_det(nJ, DiagParts(:,i)) 

            ! Search the hash table HashIndex for the determinant defined by
            ! nJ and DiagParts(:,i). If it is found, tSuccess will be
            ! returned .true. and PartInd will hold the position of the
            ! determinant in CurrentDets. Else, tSuccess will be returned
            ! .false. (and PartInd shouldn't be accessed).
            ! Also, the hash value, DetHash, is returned by this routine.
            ! tSuccess will determine whether the particle has been found or not.
            call hash_table_lookup(nJ, DiagParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            tDetermState = .false.

!            WRITE(6,*) 'i,DiagParts(:,i)',i,DiagParts(:,i)
            
            if (.true. .and. tSuccess) then

                ! Our DiagParts determinant is found in CurrentDets.

                call extract_sign(CurrentDets(:,PartInd),CurrentSign)
                call extract_sign(DiagParts(:,i),SpawnedSign)

                SignProd = CurrentSign*SpawnedSign

                tDetermState = test_flag(CurrentDets(:,PartInd), flag_deterministic)

                if (sum(abs(CurrentSign)) >= 1.e-12_dp .or. tDetermState) then
                    ! Transfer new sign across.
                   call encode_sign(CurrentDets(:,PartInd), SpawnedSign+CurrentSign)
                   call encode_sign(DiagParts(:,i), null_part)

                    ! If we are spawning onto a site and growing it, then
                    ! count that spawn for initiator purposes.
                    ! not sure if this function is even used ever
                    ! if (any(signprod > 0)) call inc_spawn_count(PartInd)

                    ! stick with the old convention on how to count deaths / 
                    ! births in the walker_death routine
                    ! i changed the sign of the child_sign when filling up 
                    ! the diag_parts list.. so i have to change that here 
                    ! again.. 
                    ! rt_iter_adapt : this will be counted below
                    ! iter_data%ndied = iter_data%ndied + &
                    !    min(-SpawnedSign,abs(CurrentSign))
                    
                    ! how does that combine with the nAborted below..? todo
                    do j = 1, lenof_sign
                        run = part_type_to_run(j)
                        ! here it seems, it treats both the real and complex 
                        ! walker occupation for the same initiator criteria
                        ! meaning if, any of the real or imaginary occupations
                        ! is an initiator, the whole determinant is an 
                        ! initiator.. hm.. do i want this? if yes i have 
                        ! to change a lot of the previous implemented code 
                        ! which distinguished between real and imaginary 
                        ! initiators..
                        ! for now, stick to the distinction! but ask ali! 
                        if (abs(CurrentSign(j)) < 1.e-12_dp) then
                           ! This determinant is actually *unoccupied* for the
                           ! walker type/set we're considering. We need to
                           ! decide whether to abort it or not.
                           if (tTruncInitiator .and. .not. test_flag (DiagParts(:,i), &
                                flag_initiator(j)) .and. .not. tDetermState) then
                              ! Walkers came from outside initiator space.
                              NoAborted(run) = NoAborted(run) + abs(SpawnedSign(j))
                              ! rt_iter_adapt:
                              ! these must not be counted as they never were considered born
                              ! iter_data%naborted(j) = iter_data%naborted(j) + abs(SpawnedSign(j))
                              call encode_part_sign (CurrentDets(:,PartInd), 0.0_dp, j)
                           else
                              ! rt_iter_adapt : also count those spawned onto an initiator
                              ! this is not necessary in the normal version as walkers are
                              ! counted there on spawn
                              ! also, they are born for any sign
                              iter_data%nborn(j) = iter_data%nborn(j) + &
                                   abs(SpawnedSign(j))
                              NoBorn(run) = NoBorn(run) + abs(SpawnedSign(j))
                           end if

                        else if (SignProd(j) < 0) then
                            ! in the real-time for the final combination
                            ! y(n) + k2 i have to check if the "spawned" 
                            ! particle is actually a diagonal death/born
                            ! walker
                            ! UPDATE! we changed that, so that we deal with 
                            ! the diagonal particles specifically, so i know
                            ! ever death/born here is a death or spawn..
                            ! This indicates that the particle has found the
                            ! same particle of opposite sign to annihilate with.
                            ! In this case we just need to update some statistics:
                            ! stick here with the convention in the original 
                            ! walker_death routine on how to count the deaths
                            ! and borns
                            ! remember the - sign when filling up the DiagParts
                            ! list -> so opposite sign here means a death!
                            iter_data%ndied(j) = iter_data%ndied(j) + &
                                min(abs(CurrentSign(j)),abs(SpawnedSign(j)))

                            NoDied(run) = NoDied(run) + min(abs(CurrentSign(j)), &
                                abs(SpawnedSign(j)))

                            ! and if the Spawned sign magnitude is even higher
                            ! then the currentsign -> born anti-particles

                            iter_data%nborn(j) = iter_data%nborn(j) + &
                                max(abs(SpawnedSign(j)) - abs(CurrentSign(j)), 0.0_dp) 

                            NoBorn(run) = NoBorn(run) + max(abs(SpawnedSign(j)) - &
                                abs(CurrentSign(j)), 0.0_dp)

!                             Annihilated(run) = Annihilated(run) + 2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j))))
!                             iter_data%nannihil(j) = iter_data%nannihil(j) + &
!                                 2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j))))

!                             if (tHistSpawn) then
!                                 ! We want to histogram where the particle
!                                 ! annihilations are taking place.
!                                 ExcitLevel = FindBitExcitLevel(DiagParts(:,i), iLutHF, nel)
!                                 if (ExcitLevel == NEl) then
!                                     call BinSearchParts2(SpawnedParts(:,i), HistMinInd2(ExcitLevel), Det, PartIndex, tSuc)
!                                 else if (ExcitLevel == 0) then
!                                     PartIndex = 1
!                                     tSuc = .true.
!                                 else
!                                     call BinSearchParts2(SpawnedParts(:,i), HistMinInd2(ExcitLevel), &
!                                             FCIDetIndex(ExcitLevel+1)-1, PartIndex, tSuc)
!                                 end if
!                                 HistMinInd2(ExcitLevel) = PartIndex
!                                 if (tSuc) then
!                                     AvAnnihil(j,PartIndex) = AvAnnihil(j,PartIndex)+ &
!                                     real(2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j)))), dp)
!                                     InstAnnihil(j,PartIndex) = InstAnnihil(j,PartIndex)+ &
!                                     real(2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j)))), dp)
!                                 else
!                                     write(6,*) "***",SpawnedParts(0:NIftot,i)
!                                     Call WriteBitDet(6,SpawnedParts(0:NIfTot,i), .true.)
!                                     call stop_all("AnnihilateSpawnedParts","Cannot find corresponding FCI "&
!                                         & //"determinant when histogramming")
!                                 end if
!                             end if
                        else
                            ! if it has the same sign i have to keep track of 
                            ! the born particles, or as in the original 
                            ! walker_death, reduce the number of died parts
                           iter_data%ndied(j) = iter_data%ndied(j) - &
                                abs(SpawnedSign(j))
                            
                        end if

                    end do ! Over all components of the sign.
                    
                    if (.not. tDetermState) then
                        call extract_sign (CurrentDets(:,PartInd), SignTemp)
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
                    !if (tFillingStochRDMonFly .and. (.not.tNoNewRDMContrib)) then
                    !call extract_sign(CurrentDets(:,PartInd), TempCurrentSign)
                        ! We must use the instantaneous value for the off-diagonal
                        ! contribution. However, we can't just use CurrentSign from
                        ! the previous iteration, as this has been subject to death
                        ! but not the new walkers. We must add on SpawnedSign, so
                        ! we're effectively taking the instantaneous value from the
                        ! next iter. This is fine as it's from the other population,
                        ! and the Di and Dj signs are already strictly uncorrelated.
                    !    call check_fillRDM_DiDj(rdms, i, CurrentDets(:,PartInd), TempCurrentSign)
                    !end if 

                end if

            end if
                
            if (((.not.tSuccess) .or. (tSuccess .and. sum(abs(CurrentSign)) < 1.e-12_dp .and. (.not. tDetermState)))) then

                ! Determinant in newly spawned list is not found in CurrentDets.
                ! Usually this would mean the walkers just stay in this list and
                ! get merged later - but in this case we want to check where the
                ! walkers came from - because if the newly spawned walkers are
                ! from a parent outside the active space they should be killed,
                ! as they have been spawned on an unoccupied determinant.

                ! this can also happen in the rt-fciqmc, if this is a diagonal
                ! spawn from a excitation in the first loop.. 
                ! can think about specific rules here.. maybe i dont want 
                ! to do that here or use some different kind of initiator
                ! criteria, to do this spawn... todo and talk with ali! 
                if (tTruncInitiator) then

                    call extract_sign (DiagParts(:,i), SignTemp)

                    tPrevOcc = .false.
                    if (.not. IsUnoccDet(SignTemp)) tPrevOcc=.true.   
                        
                    do j = 1, lenof_sign
                        run = part_type_to_run(j)

                        if (test_abort_spawn(DiagParts(:, i), j)) then

                            ! If this option is on, include the walker to be
                            ! cancelled in the trial energy estimate.
!                             if (tIncCancelledInitEnergy) call add_trial_energy_contrib(DiagParts(:,i), SignTemp(j), j)

                            ! Walkers came from outside initiator space.
                            NoAborted(run) = NoAborted(run) + abs(SignTemp(j))
                            ! rt_iter_adapt : As above, aborted walkers must not be counted
                            ! iter_data%naborted(j) = iter_data%naborted(j) + abs(SignTemp(j))
                            ! We've already counted the walkers where SpawnedSign
                            ! become zero in the compress, and in the merge, all
                            ! that's left is those which get aborted which are
                            ! counted here only if the sign was not already zero
                            ! (when it already would have been counted).
                            SignTemp(j) = 0.0_dp
                            call encode_part_sign (DiagParts(:,i), SignTemp(j), j)

                        end if

                        
                        ! RT_M_Merge: Removed initiator case
                        ! Either the determinant has never been an initiator,
                        ! or we want to treat them all the same, as before.
                        ! this is for the real-coefficient.. if its 
                        ! below 1.0_dp eg.. also not of concern yet 
                        if ((abs(SignTemp(j)) > 1.e-12_dp) .and. (abs(SignTemp(j)) < OccupiedThresh)) then
                           ! We remove this walker with probability 1-RealSignTemp
                           pRemove=(OccupiedThresh-abs(SignTemp(j)))/OccupiedThresh
                           r = genrand_real2_dSFMT ()
                           if (pRemove > r) then
                              ! Remove this walker.
                              NoRemoved(run) = NoRemoved(run) + abs(SignTemp(j))
                              !Annihilated = Annihilated + abs(SignTemp(j))
                              !iter_data%nannihil = iter_data%nannihil + abs(SignTemp(j))
                              
                              !rt_iter_adapt : see above
                              !iter_data%nremoved(j) = iter_data%nremoved(j) &
                              !     + abs(SignTemp(j))
                              SignTemp(j) = 0.0_dp
                              call nullify_ilut_part (DiagParts(:,i), j)
                           else! if (tEnhanceRemainder) then
                              NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(SignTemp(j))
                              iter_data%nborn(j) = iter_data%nborn(j) &
                                   + OccupiedThresh - abs(SignTemp(j))
                              SignTemp(j) = sign(OccupiedThresh, SignTemp(j))
                              call encode_part_sign (DiagParts(:,i), SignTemp(j), j)
                           end if
                        end if
                    end do
                    if (.not. IsUnoccDet(SignTemp)) then
                        ! Walkers have not been aborted and so we should copy the
                        ! determinant straight over to the main list. We do not
                        ! need to recompute the hash, since this should be the
                        ! same one as was generated at the beginning of the loop.
                        ! but here i have to keep track of the number of born
                        ! and died particles... or? 
                        ! hm.. in the first spawn.. i count them as spawns or..
                        ! so can i here count them as born.. deaths..
                        ! and if it is not found here .. 
                        ! atleast it can never be a death, as it is not 
                        ! found in the main list..
                        ! for now, just count them as birth

                       iter_data%nborn = iter_data%nborn + abs(SignTemp)
                       
                       NoBorn(run) = NoBorn(run) + sum(abs(SignTemp))

                       call AddNewHashDet(TotWalkersNew, DiagParts(:,i), DetHash, nJ)
                    end if

                 else
                    ! Running the full, non-initiator scheme.
                    ! Determinant in newly spawned list is not found in
                    ! CurrentDets. If coeff <1, apply removal criterion.
                    call extract_sign (DiagParts(:,i), SignTemp)
                    
                    tPrevOcc = .false.
                    if (.not. IsUnoccDet(SignTemp)) tPrevOcc = .true. 
                    
                    do j = 1, lenof_sign
                        run = part_type_to_run(j)
                        if ((abs(SignTemp(j)) > 1.e-12_dp) .and. (abs(SignTemp(j)) < OccupiedThresh)) then
                            ! We remove this walker with probability 1-RealSignTemp.
                            pRemove = (OccupiedThresh-abs(SignTemp(j)))/OccupiedThresh
                            r = genrand_real2_dSFMT ()
                            if (pRemove  >  r) then
                                ! Remove this walker.
                                NoRemoved(run) = NoRemoved(run) + abs(SignTemp(j))
                                ! rt_iter_adapt : see above
                                !Annihilated = Annihilated + abs(SignTemp(j))
                                !iter_data%nannihil = iter_data%nannihil + abs(SignTemp(j))
                                !iter_data%nremoved(j) = iter_data%nremoved(j) &
                                !+ abs(SignTemp(j))
                                SignTemp(j) = 0
                                call nullify_ilut_part (DiagParts(:,i), j)
                             else !if (tEnhanceRemainder) then
                                NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(SignTemp(j))
                                iter_data%nborn(j) = iter_data%nborn(j) &
                                            + OccupiedThresh - abs(SignTemp(j))
                                SignTemp(j) = sign(OccupiedThresh, SignTemp(j))
                                call encode_part_sign (DiagParts(:,i), SignTemp(j), j)
                            end if
                        end if
                    end do
                    
                    if (.not. IsUnoccDet(SignTemp)) then
                        ! Walkers have not been aborted and so we should copy the
                        ! determinant straight over to the main list. We do not
                        ! need to recompute the hash, since this should be the
                        ! same one as was generated at the beginning of the loop.
                        ! also here treat those new walkers as born particles

                       iter_data%nborn = iter_data%nborn + abs(SignTemp)
                       
                       NoBorn(run) = NoBorn(run) + sum(abs(SignTemp))
                       call AddNewHashDet(TotWalkersNew, DiagParts(:,i), DetHash, nJ)

                    end if
                end if

!                if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) then
                    ! We must use the instantaneous value for the off-diagonal contribution.
!                call check_fillRDM_DiDj(rdms, i, SpawnedParts(0:NifTot,i), SignTemp)
!                end if 
            end if

        end do

!         call halt_timer(BinSearch_time)

        ! Update remaining number of holes in list for walkers stats.
        if (iStartFreeSlot > iEndFreeSlot) then
            ! All slots filled
            HolesInList = 0
        else
            HolesInList = iEndFreeSlot - (iStartFreeSlot-1)
        endif

!         call halt_timer(AnnMain_time)

    end subroutine AnnihilateDiagParts

    function count_holes_in_currentDets() result(holes)
        integer :: holes 
        integer(n_int), pointer :: ilut_parent(:)
        integer :: nI_parent(nel), unused_flags, idet
        real(dp) :: parent_sign(lenof_sign)

        holes = 0

        do idet = 1, int(TotWalkers, sizeof_int)

            ilut_parent => CurrentDets(:,idet) 

            call extract_bit_rep(ilut_parent, nI_parent, parent_sign, unused_flags, &
                fcimc_excit_gen_store)

            if (IsUnoccDet(parent_sign)) then
                holes = holes + 1
            end if

        end do

    end function count_holes_in_currentDets

    subroutine create_diagonal_as_spawn(nI, ilut, diag_sign, iter_data)
        ! new routine to create diagonal particles into new DiagParts 
        ! array to distinguish between spawns and diagonal events in the 
        ! combination y(n) + k2
      use Parallel_neci, only :iProcIndex
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:niftot)
        real(dp), intent(in) :: diag_sign(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "create_diagonal_as_spawn"

        integer :: proc, i
        logical :: list_full
        integer, parameter :: flags = 0

        ! Determine which processor the particle should end up on in the
        ! DirectAnnihilation algorithm.
        ! Note that this is a diagonal event, no communication is needed
        proc = iProcIndex    ! (0 -> nNodes-1)

        ! Check that the position described by valid_diag_spawn_list is acceptable.
        ! If we have filled up the memory that would be acceptable, then
        ! kill the calculation hard (i.e. stop_all) with a descriptive
        ! error message.
        list_full = .false.
        if (proc == nNodes - 1) then
            if (valid_diag_spawn_list(proc) > MaxSpawned) list_full = .true.
        else
            if (valid_diag_spawn_list(proc) >= InitialSpawnedSlots(proc+1)) &
                list_full=.true.
        end if
        if (list_full) then
            write(6,*) "Attempting to spawn particle onto processor: ", proc
            write(6,*) "No memory slots available for this spawn."
            write(6,*) "Please increase MEMORYFACSPAWN"
            call stop_all(this_routine, "Out of memory for spawned particles")
        end if

        call encode_bit_rep(DiagParts(:, valid_diag_spawn_list(proc)), ilut, &
                            diag_sign, flags)

        ! If the parent was an initiator then set the initiator flag for the
        ! child, to allow it to survive.
        ! in the real-time, depending if damping is present or not both 
        ! particle types can or can't spawn to both types of child-particles
        if (tTruncInitiator) then
           ! rmneci_setup: changed 1,2 -> max/min_part_type(run)
           ! for rotated time: real_time_info%damping -> tRotateTime
            if (real_time_info%damping > EPS) then
                ! both particle spawns possible -> check for both initiator flags! 
                do i = 1, lenof_sign
                    if (abs(diag_sign(i)) > EPS) then
                       ! rmneci_setup: these are the flags for the corresponding run
                        if (test_flag(ilut, flag_initiator((1+i)/2)) .or. &
                            test_flag(ilut, flag_initiator(1+(1+i)/2))) then
                            call set_flag(DiagParts(:,valid_diag_spawn_list(proc)), &
                                flag_initiator(i))
                        end if
                    end if
                end do
            else
                ! in this case only spawns to the other type of particle are 
                ! possible! check for the init of the other parent type 
                do i = 1, lenof_sign
                    if (abs(diag_sign(i)) > EPS .and. test_flag(ilut, &
                        flag_initiator(rotate_part(i)))) then
                        call set_flag(DiagParts(:, valid_diag_spawn_list(proc)), &
                            flag_initiator(i))
                    end if
                end do
            end if
        end if

        if (tFillingStochRDMonFly) then
            call stop_all(this_routine, "RDM not yet implemented in the rt-fciqmc!")
            ! We are spawning from ilutI to 
            ! SpawnedParts(:,valid_diag_spawn_list(proc)). We want to store the
            ! parent (D_i) with the spawned child (D_j) so that we can add in
            ! Ci.Cj to the RDM later.
            ! The parent is NIfDBO integers long, and stored in the second
            ! part of the SpawnedParts array from NIfTot+1 --> NIfTot+1+NIfDBO
!             call store_parent_with_spawned (RDMBiasFacCurr, WalkerNo, &
!                                             ilutI, WalkersToSpawn, ilutJ, &
!                                             proc, part_type)

        end if

        ! If we are storing the parent coefficient with the particle, then
        ! do that it this point
!         if (tBroadcastParentCoeff) then
!             ! n.b. SignCurr(part_type) --> this breaks with CPLX
!             call stop_all(this_routine, 'Not implemented (yet)')
!             call set_parent_coeff(DiagParts(:, valid_diag_spawn_list(proc)), &
!                                   SignCurr(part_type))
!         end if

! #ifdef __CMPLX
!         if (tTruncInitiator) then
!             ! With complex walkers, things are a little more tricky.
!             ! We want to transfer the flag for all particles created (both
!             ! real and imag) from the specific type of parent particle. This
!             ! can mean real walker flags being transfered to imaginary
!             ! children and vice versa.
!             ! This is unneccesary for real walkers.
!             ! Test the specific flag corresponding to the parent, of type
!             ! 'part_type'
!             ! ok.. here it is obvious that the childs get all the flags from
!             ! the parent and not only one type..
!             ! but i essentially do that above, by checking both! 
!             ! so this is redundant
!             parent_init = test_flag(SpawnedParts(:,valid_diag_spawn_list(proc)), &
!                                     flag_initiator(part_type))
!             ! Assign this flag to all spawned children.
!             do j=1,lenof_sign
!                 if (child(j) /= 0) then
!                     call set_flag (SpawnedParts(:,valid_diag_spawn_list(proc)), &
!                                    flag_initiator(j), parent_init)
!                 endif
!             enddo
!         end if
! #endif

        valid_diag_spawn_list(proc) = valid_diag_spawn_list(proc) + 1
        
        ! Sum the number of created children to use in acceptance ratio.
        ! i dont need to do acceptances here, since diagonal events are not 
        ! counted in the acceptance ratio..
!         acceptances(1) = acceptances(1) + sum(abs(diag_sign))

    end subroutine create_diagonal_as_spawn

    subroutine create_diagonal_as_spawn_old(nI, ilut, diag_sign, iter_data)
        ! routine to treat the diagonal death/cloning step in the 2nd 
        ! RK loop as spawned particles, with using hash table too
        ! have to write this routine, since i also need the old one
        ! in the "normal" spawnings of the rt-fciqmc
        ! only need the currently looped over determinant! and probably 
        ! also just have to copy all the sign(or just store the original 
        ! parent_ilut in the hash table
        ! it does a simultanious spawning of real and complex walkers! 
        ! UPDATE: change the way this "diagonal spawning" is done to store
        ! it in a seperate SpawnedPartsDiag array! so its easier to keep 
        ! track of the 2 seperate processes death/cloning and spawning
        ! the good thing here is that i do not need to check if a diagonal 
        ! spawn is already in the hash-list of "spawned" parts. since its 
        ! only diagonal! and when i write a new "annihilation" routine i 
        ! can treat those particles just as "normal" death or birth processes
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:niftot)
        real(dp), intent(in) :: diag_sign(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "create_diagonal_as_spawn_old"

        integer :: proc, ind, hash_val, i
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_old(lenof_sign), real_sign_new(lenof_sign)
        real(dp) :: sgn_prod(lenof_sign)
        logical :: list_full, tSuccess
        integer, parameter :: flags = 0       

        ! i have to work with the DiagParts array here, not the Spawned Parts
        ! i do not need to use a hash, as only diagonal events happen
        
        call hash_table_lookup(nI, ilut, nifdbo, spawn_ht, SpawnedParts, ind, &
            hash_val, tSuccess)

        if (tSuccess) then
            ! If the spawned child is already in the spawning array.
            ! Extract the old sign.
            call extract_sign(SpawnedParts(:,ind), real_sign_old)

            ! If the new child has an opposite sign to that of walkers already
            ! on the site, then annihilation occurs. The stats for this need
            ! accumulating.
            sgn_prod = real_sign_old * diag_sign
            ! this annihilation is kinda right.. it annihilates previously 
            ! spawned particles in the spawned array...
            do i = 1, lenof_sign
                if (sgn_prod(i) < 0.0_dp) then
                    iter_data%nannihil(i) = iter_data%nannihil(i) + &
                        2*min( abs(real_sign_old(i)), abs(diag_sign(i)) )
                end if
            end do

            ! Find the total new sign.
            real_sign_new = real_sign_old + diag_sign
            ! Encode the new sign.
            call encode_sign(SpawnedParts(:,ind), real_sign_new)

            ! i also have to keep track of the nDied and nBorn in here if 
            ! i treat these diagonal terms as "spawns" ...
            ! but right now, as calced above they are treated as annihilations!
            ! to deal with it later in the annihilation step as nborns and 
            ! ndied mark them as diagonal particles here.. 
            ! Set the initiator flags appropriately.
            ! If this determinant (on this replica) has already been spawned to
            ! then set the initiator flag. Also if this child was spawned from
            ! an initiator, set the initiator flag.
            ! in this routine i do a simultanious spawning of real and complex 
            ! walkers -> and here i just considere the diagonal "spawn" of 
            ! already occupied dets..
            ! shouldnt i also check if this didt annihilate all the walkers? 
            ! or is this handled later in the annihilation/calchashtable 
            ! routines?
            if (tTruncInitiator) then
                ! i can still check if only a real hamiltonian is used, so 
                ! i know all the diagonal "spawns" come from the other particle
                ! type

               ! rotated_time: this will need adaption
                if (real_time_info%damping > EPS) then
                    ! update: i am stupid! there are no complex diagonal 
                    ! elements! duh.. hermiticity!!
                    ! BUT there is still the damping factor, which leads to 
                    ! 'spawns' to both particle type from each particle
                    ! both spawns possible
                    ! can i do that with the simultanious spawn? since i have 
                    ! no information left here what the parent was ..
                    ! fuck.. -> so for now dont consider those pesky complex
                    ! Hamiltonians! 
                    ! see update above! 
                    do i = 1, lenof_sign
                        if (abs(real_sign_new(i)) > EPS) then
                            if (abs(real_sign_old(i)) > EPS .or. &
                                test_flag(ilut,flag_initiator(min_part_type(i))) .or. &
                                test_flag(ilut,flag_initiator(max_part_type(i)))) then

                                call set_flag(SpawnedParts(:,ind),flag_initiator(i))
                            end if
                        end if
                    end do

!                     call stop_all(this_routine, &
!                         "i loose track of what the parent det is, if i do the death step simultaniously!")
                else
                    ! this routine gets called only if some part of the 
                    ! child is non-zero, so i am sure there was some entry
                    ! non-zero in the child.. but since the spawns of 
                    ! both particle types are done simultaniously i need to 
                    ! check still.. 
                    do i = 1, lenof_sign 
                        if (abs(real_sign_new(i)) > EPS) then
                            ! then i know there is something left on the det
                            if (abs(real_sign_old(i)) > EPS .or. &
                                test_flag(ilut, flag_initiator(rotate_part(i)))) then
                                ! if 2 spawns to the same, or parent (of 
                                ! opposite part. type) was initiator
                                call set_flag(SpawnedParts(:,ind), flag_initiator(i))
                            end if
                        end if
                    end do
                end if
            end if
        else
            ! Determine which processor the particle should end up on in the
            ! DirectAnnihilation algorithm.
            proc = DetermineDetNode(nel, nI, 0)

            ! i also have to count the born/died walkers here! 
            ! i will later on annihilate this spawned list with the original
            ! y(n) list.. where it would be treated as an annihilation again
            ! but here i know that this is an diagonal step actually... 

            ! fuck! by the way: also with damping i have the spawn to both 
            ! particle types from each type of walker...
            ! todo: i have to rework this diagonal step.. because already 
            ! with damping only i loose track of the parent type.. 
            ! maybe also the 'regular' death step in the first loop of the 
            ! real-time fciqmc...

            ! Check that the position described by ValidSpawnedList is acceptable.
            ! If we have filled up the memory that would be acceptable, then
            ! kill the calculation hard (i.e. stop_all) with a descriptive
            ! error message.
            list_full = .false.
            if (proc == nNodes - 1) then
                if (ValidSpawnedList(proc) > MaxSpawned) list_full = .true.
            else
                if (ValidSpawnedList(proc) > InitialSpawnedSlots(proc+1)) &
                    list_full=.true.
            end if
            if (list_full) then
                write(6,*) "Attempting to spawn particle onto processor: ", proc
                write(6,*) "No memory slots available for this spawn."
                write(6,*) "Please increase MEMORYFACSPAWN"
                call stop_all(this_routine, "Out of memory for spawned particles")
            end if

            call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), &
                ilut(0:NIfDBO), diag_sign, flags)

            ! If the parent was an initiator then set the initiator flag for the
            ! child, to allow it to survive.
            if (tTruncInitiator) then
                ! i have to consider the 'other particle type as parent.. 
                ! so check that initiator flag! 
                if (real_time_info%damping > EPS) then
                    ! for now assume that both of the parents contribute some
                    ! progeny and check initiator flags of both parent parts.
                    do i = 1, lenof_sign
                        if (abs(diag_sign(i)) > EPS) then
                            if (test_flag(ilut,flag_initiator(max_part_type(i))) .or. &
                                test_flag(ilut,flag_initiator(min_part_type(i)))) then
                                call set_flag(SpawnedParts(:,ValidSpawnedList(proc)), &
                                    flag_initiator(i))
                            end if
                        end if
                    end do
!                     call stop_all(this_routine, &
!                         "i loose info on which the parent is if i do the death_as_spawn simultaniously for both types!")
                else
                    do i = 1, lenof_sign
                        ! i have to check if the specific element is > 0
                        if (abs(diag_sign(i)) > EPS .and. &
                            test_flag(ilut, flag_initiator(rotate_part(i)))) then
                            
                            call set_flag(SpawnedParts(:,ValidSpawnedList(proc)),&
                                flag_initiator(i))
                        end if
                    end do
                end if
            end if

            call add_hash_table_entry(spawn_ht, ValidSpawnedList(proc), hash_val)

            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        end if
        
        ! Sum the number of created children to use in acceptance ratio.
        ! how is this done in the rt-fciqmc??
        acceptances(1) = acceptances(1) + sum(abs(diag_sign))

    end subroutine create_diagonal_as_spawn_old

    function attempt_die_realtime(DetCurr, Kii, RealwSign, walkExcitLevel) &
            result(ndie)
        ! also write a function, which calculates the new "signs"(weights) of 
        ! the real and complex walkers for the diagonal death/cloning step 
        ! since i need that for both 1st and 2nd loop of RK, but at different 
        ! points 
      implicit none
        integer, intent(in) :: DetCurr(nel) 
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        real(dp), intent(in) :: Kii
        real(dp), dimension(lenof_sign) :: ndie
        real(dp), dimension(lenof_sign) :: CopySign
        integer, intent(in) :: walkExcitLevel
        integer :: run, irdm
        ! each run needs its own fac, since the shift varies
        real(dp) :: fac(lenof_sign), rat, r

        character(*), parameter :: this_routine = "attempt_die_realtime"

        ! do it dependent on if there is damping, since otherwise i can save
        ! some effort..

        ! do i need a minus sign here?? just from the convention
        ! in the rest of the neci code? yes!
        do run = 1, inum_runs
           fac(min_part_type(run)) = -tau * (Kii - DiagSft(run))
        enddo

        if (real_time_info%damping < EPS) then
            ! in this case there is only Re <-> Im 
            ! n_i' =  H_ii c_i 
            ! c_i' = -H_ii n_i
            ! how should i log death_magnitudes here.. for tau-search.. 
            ! that probably has to change too..
            ! an actual death is only done when the damping factor is present
            ! as otherwise we 'only' spawn to the other particle type
!             call log_death_magnitude(Kii - DiagSft(1))

            ! and also not sure about this criteria.. if that applies in the 
            ! rt-fciqmc..
            if(any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTau) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do run = 1, inum_runs
                            fac(run) = min(2.0_dp, fac(run))
                        end do
                    else
                        call stop_all(this_routine, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    write(iout,'("** WARNING ** Death probability > 1: Creating Antiparticles. "&
                        & //"Timestep errors possible: ")',advance='no')
                    do run = 1, inum_runs
                        write(iout,'(1X,f13.7)',advance='no') fac(run)
                    end do
                    write(iout,'()')
                endif
            endif
            do run = 1, inum_runs
               if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel .le. RealCoeffExcitThresh)) &
                    .or. tAllRealCoeff ) then

                  ! i exact.. just get the weights, this should also still work.
                  ! but have to exchange the weights to come from the other 
                  ! type of particles..
                  ! the number of deaths has +sign from Im -> Re
                  ! dont use abs() here compared to the old code, but already 
                  ! here include the sign of the walkers on the determinant
                  ndie(min_part_type(run)) = fac(min_part_type(run)) * realwSign(max_part_type(run))
                  ! and - from Re -> Im
                  ! does this give the correct sign compared to the parent sign?
                  ndie(max_part_type(run)) = -fac(min_part_type(run)) * realwSign(min_part_type(run))

               else 
                  ! if not exact i have to round stochastically
                  ! Re -> Im
                  rat = fac(min_part_type(run)) * RealwSign(max_part_type(run))

                  ndie(min_part_type(run)) = real(int(rat), dp) 
                  rat = rat - ndie(min_part_type(run))

                  r = genrand_real2_dSFMT() 
                  if (abs(rat) > r) ndie(min_part_type(run)) = &
                       ndie(min_part_type(run)) + real(nint(sign(1.0_dp,rat)),dp)

                  ! Im -> Re
                  rat = -fac(min_part_type(run)) * RealwSign(max_part_type(run))
                  ndie(max_part_type(run)) = real(int(rat), dp)
                  rat = rat - ndie(max_part_type(run))
                  r = genrand_real2_dSFMT()
                  if (abs(rat) > r) ndie(max_part_type(run)) = &
                       ndie(max_part_type(run)) + real(nint(sign(1.0_dp,rat)),dp)
               end if
            enddo

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
           do run = 1, inum_runs
              fac(max_part_type(run)) = tau * real_time_info%damping 
           end do
            ! here i am definetly not sure about the logging of the death
            ! magnitude.. 
            ! but also with damping.. the death is only related to the 
            ! damping and not on shift or diagonal matrix element..
            ! so i dont think i need it here! 

            ! and also about the fac restrictions.. for now but it here anyway..
            if(any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTau) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do run = 1, inum_runs
                            fac(run) = min(2.0_dp, fac(run))
                        end do
                    else
                        call stop_all(this_routine, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    write(iout,'("** WARNING ** Death probability > 1: Creating Antiparticles. "&
                        & //"Timestep errors possible: ")',advance='no')
                    do run = 1, inum_runs
                        write(iout,'(1X,f13.7)',advance='no') fac(run)
                    end do
                    write(iout,'()')
                endif
            endif
            do run = 1, inum_runs
               if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel .le. RealCoeffExcitThresh)) &
                    .or. tAllRealCoeff ) then

                  ! i exact.. just get the weights, this should also still work.
                  ! but have to exchange the weights to come from the other 
                  ! type of particles..
                  ! the number of deaths has +sign from Im -> Re
                  ! can i just add the other contribution here? 
                  ! and also include the sign of the parent occupation here
                  ! already.
                  ndie(min_part_type(run)) = fac(min_part_type(run)) &
                       * realwSign(max_part_type(run)) + &
                       fac(max_part_type(run)) * realwSign(min_part_type(run))
                  ! and - from Re -> Im
                  ! does this give the correct sign compared to the parent sign?
                  ndie(max_part_type(run)) = -fac(min_part_type(run)) * &
                       abs(realwSign(min_part_type(run))) + &
                       fac(max_part_type(run)) * abs(RealwSign(max_part_type(run)))

               else 
                  ! if not exact i have to round stochastically
                  ! is this ok here to just add the second contribution? todo
                  ! Re -> Im
                  rat = fac(min_part_type(run)) * RealwSign(max_part_type(run)) &
                       + fac(max_part_type(run)) * RealwSign(1)

                  ndie(min_part_type(run)) = real(int(rat), dp) 
                  rat = rat - ndie(min_part_type(run))

                  r = genrand_real2_dSFMT() 
                  if (abs(rat) > r) ndie(min_part_type(run)) = ndie(min_part_type(run)) &
                       + real(nint(sign(1.0_dp,rat)),dp)

                  ! Im -> Re
                  rat = -fac(min_part_type(run)) * RealwSign(min_part_type(run)) &
                  + fac(max_part_type(run)) * RealwSign(max_part_type(run))

                  ndie(max_part_type(run)) = real(int(rat), dp)
                  rat = rat - ndie(max_part_type(run))

                  r = genrand_real2_dSFMT()
                  if (abs(rat) > r) ndie(max_part_type(run)) = &
                       ndie(max_part_type(run)) + real(nint(sign(1.0_dp,rat)),dp)
               end if
            enddo
        end if

    end function attempt_die_realtime

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
        integer :: i, irdm, j
        real(dp) :: fac(lenof_sign), rat, r

        character(len=*), parameter :: t_r = "walker_death_realtime"

        integer(n_int), pointer :: ilut(:)
        logical :: set_init(lenof_sign)

        ilut => CurrentDets(:,DetPosition)

        ndie = attempt_die_realtime(DetCurr, Kii, RealwSign, walkExcitLevel)
        ! combine in here old walker_death() and attempt_die() routines..
        ! update: no dont! because i need the attempt_die in the 2nd RK loop 
        ! alone
        ! update deatch counter..
        ! how to do that in rt-fciqmc? use 2 different infos or only do 
        ! that at the end of the full loop? so i wouldnt have to do 
        ! anything here..

        ! the only particles, which definetly die are the diagonal
        ! damped particles.. the other could have the same sign as the 
        ! walker which already occupy the determinant..
        ! but i have to integrate that above, since have unmix the both 
        ! influences.. 
!         if (real_time_info%damping > EPS) then
            ! and only if there is damping -> otherwise no 'death' 
            ! or do i have to check if the spawn from Re <-> Im caused 
            ! annihilated particles?!
            ! the previous counter for died particles assumed an opposite sing
            ! between the two.. but in the rt-fciqmc i have a mix because 
            ! of the Re <-> Im 'diagonal' spawning.. should i count that 
            ! as spawing or as death/cloning.. how to set intitator? 

            ! different to the old algorithm i include the sign of the 
            ! parent walker already in the ndie variable, so i have to 
            ! consider that here.. 
            ! and if 'anti-particles' get born the number of died particles
            ! naturally are reduced

            ! also do it independent of the damping and also consider the 
            ! diagonal switch between Re <-> Im as born or died particles

        ! rt_iter_adapt: this will be counted below
        !iter_data%ndied = iter_data%ndied + min(ndie, abs(RealwSign))

        ! this routine only gets called in the first runge-kutta step -> 
        ! so only update the stats for the first here!
        do i = 1, lenof_sign
            ! check if the parent and ndie have the same sign 
            if (sign(1.0_dp,RealwSign(i)) == sign(1.0_dp, ndie(i))) then
                ! then the entries in ndie kill the parent, but only maximally
                ! the already occupying walkers can get killed 
                iter_data%ndied(i) = iter_data%ndied(i) + &
                    abs(min(abs(RealwSign(i)),abs(ndie(i))))

                NoDied_1(part_type_to_run(i)) = NoDied_1(part_type_to_run(i)) &
                     + abs(min(abs(ndie(i)),abs(RealwSign(i))))
                ! if ndie is bigger than the original occupation i am actually
                ! spawning 'anti-particles' which i have to count as born
                ! and reduce the ndied number.. or not? 
                ! hm the old code is actually not counting births, due to 
                ! the shift.. interesting.. but just subtracts that from 
                ! the ndied quantity... 
                iter_data%nborn(i) = iter_data%nborn(i) + &
                    max(abs(ndie(i)) - abs(RealwSign(i)), 0.0_dp)

                ! not sure if i even need this quantity..
                NoBorn_1(part_type_to_run(i)) = NoBorn_1(part_type_to_run(i)) + &
                     max(abs(ndie(i)) - abs(RealwSign(i)), 0.0_dp)
            else
                ! if they have opposite sign, as in the original algorithm
                ! reduce the number of ndied by that amount 
                iter_data%ndied(i) = iter_data%ndied(i) - abs(ndie(i))

                NoDied_1(part_type_to_run(i)) = NoDied_1(part_type_to_run(i)) &
                     - abs(ndie(i))

            end if
        end do

!         NoDied_1(1) = NoDied_1(1) + sum(min(iDie, abs(RealwSign)))

        ! Count any antiparticles
!         iter_data%nborn = iter_data%nborn + max(iDie - abs(RealwSign), 0.0_dp)

!         end if
! 
        ! try to 'just' apply the rest of the code here.. 

        CopySign = RealwSign - nDie
        ! Calculate new number of signed particles on the det.
!         CopySign = RealwSign + (nDie * sign(1.0_dp, RealwSign))
!         CopySign = RealwSign - (nDie * sign(1.0_dp, RealwSign))
!         CopySign = RealwSign + ([nDie(1)* sign(1.0_dp, RealwSign(2)), &
!             ndie(2)*sign(1.0_dp,RealwSign(1))])

        ! In the initiator approximation, abort any anti-particles.
        ! do i want to do that also in the rt-fciqmc? 
        ! thats all a different .. and no i dont think this should be 
        ! aborted in this case..
        if (tTruncInitiator .and. any(CopySign /= 0)) then
            do i = 1, lenof_sign
                if (sign(1.0_dp, CopySign(i)) /= &
                        sign(1.0_dp, RealwSign(i))) then
!                     print *, "anti-particle even! but do not cancel this in the rt-fciqmc!"
!                     NoAborted = NoAborted + abs(CopySign(i))
!                     iter_data%naborted(i) = iter_data%naborted(i) &
!                                           + abs(CopySign(i))
!                     if (test_flag(ilutCurr, flag_initiator(i))) &
!                         NoAddedInitiators = NoAddedInitiators - 1
!                     CopySign(i) = 0
                end if
            end do
        end if

        ! i also have to update initiator information, since there is a lot 
        ! of population/depopulation between the Re <-> Im parts of an 
        ! occupied det, and i don't think this is done in any different part 
        ! of algorithm.. todo
        if (any(CopySign /= 0)) then
            ! For the hashed walker main list, the particles don't move.
            ! Therefore just adjust the weight.
            call encode_sign (CurrentDets(:, DetPosition), CopySign)
            ! i also have to update the initiator information, due to the 
            ! big (de)population beteen real and complex walker populations..
            ! and i also have to deal with allowed anti-particles in the 
            ! rt-fciqmc even with the initiator method in use! 
            if (tTruncInitiator) then
                ! if the other particle type was an initiator, the newly 
                ! spawned child should also become an initiator
                ! or only if the occupation rises above the initiator 
                ! threshhold it should also become one
                do i = 1, lenof_sign
                    ! other particle type
                    j = rotate_part(i)
                    if (abs(CopySign(j)) > EPS .and. test_flag(ilut, flag_initiator(i))) then
                        set_init(j) = .true.
                    end if
                end do
                do i = 1, lenof_sign
                    if (set_init(i)) then
                        call set_flag(ilut,flag_initiator(i))
                    end if
                end do
            end if
        else
            ! All walkers died.
!            if(tFillingStochRDMonFly) then
!                do irdm = 1, nrdms
!                    call det_removed_fill_diag_rdm(rdms(irdm), irdm, CurrentDets(:,DetPosition), DetPosition)
!                end do
                ! Set the average sign and occupation iteration to zero, so
                ! that the same contribution will not be added in in
                ! CalcHashTableStats, if this determinant is not overwritten
                ! before then
!                global_determinant_data(:, DetPosition) = 0.0_dp
!            endif

            if (tTruncInitiator) then
                ! All particles on this determinant have gone. If the determinant was an initiator, update the stats
                if (test_flag(iLutCurr,flag_initiator(1))) then
                    NoAddedInitiators = NoAddedInitiators - 1
                else if (test_flag(iLutCurr,flag_initiator(lenof_sign))) then
                    NoAddedInitiators(inum_runs) = NoAddedInitiators(inum_runs) - 1
                end if
            end if

            ! Remove the determinant from the indexing list
            call remove_hash_table_entry(HashIndex, DetCurr, DetPosition)
            ! Add to the "freeslot" list
            iEndFreeSlot = iEndFreeSlot + 1
            FreeSlot(iEndFreeSlot) = DetPosition
            ! Encode a null det to be picked up
            call encode_sign(CurrentDets(:,DetPosition), null_part)
        end if


    end subroutine walker_death_realtime

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

    function attempt_create_realtime(DetCurr, iLutCurr, RealwSign, nJ, iLutnJ,&
            prob, HElGen, ic, ex, tParity, walkExcitLevel, part_type, AvSignCurr,&
            RDMBiasFacCurr) result(child)

        ! create a specific attempt_create function for the real-time fciqmc
        ! to avoid preprocessor flag jungle..

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! 1 = Real parent particle, 2 = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t(dp) , intent(in) :: HElGen
        character(*), parameter :: this_routine = 'attempt_create_realtime'

        real(dp) :: rat, r, walkerweight, pSpawn, nSpawn, MatEl, p_spawn_rdmfac
        integer :: extracreate, tgt_cpt, component, i, iUnused
        integer :: TargetExcitLevel
        logical :: tRealSpawning
        HElement_t(dp) :: rh, rh_used

        ! Just in case
        child = 0.0_dp

        ! If each walker does not have exactly one spawning attempt
        ! (if AvMCExcits /= 1.0_dp) then the probability of an excitation
        ! having been chosen, prob, must be altered accordingly.
        prob = prob * AvMCExcits

        ! In the case of using HPHF, and when tGenMatHEl is on, the matrix
        ! element is calculated at the time of the excitation generation, 
        ! and returned in HElGen. In this case, get_spawn_helement simply
        ! returns HElGen, rather than recomputing the matrix element.
        rh = get_spawn_helement (DetCurr, nJ, iLutCurr, iLutnJ, ic, ex, &
                                 tParity, HElGen)

        !write(6,*) 'p,rh', prob, rh

        ! The following is useful for debugging the contributions of single
        ! excitations, and double excitations of spin-paired/opposite
        ! electron pairs to the value of tau.
!        if (ic == 2) then
!            if (G1(ex(1,1))%Ms /= G1(ex(1,2))%Ms) then
!                write(6,*) 'OPP', rh, prob
!            else
!                write(6,*) 'SAM', rh, prob
!            end if
!        else
!            write(6,*) 'IC1', rh, prob
!        end if

        ! Are we doing real spawning?
        
        tRealSpawning = .false.
        if (tAllRealCoeff) then
            tRealSpawning = .true.
        elseif (tRealCoeffByExcitLevel) then
            TargetExcitLevel = FindBitExcitLevel (iLutRef, iLutnJ)
            if (TargetExcitLevel <= RealCoeffExcitThresh) &
                tRealSpawning = .true.
        endif

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
        if (.not. t_complex_ints) then
            ! if it is a pure real-hamiltonian there is only spawing from 
            ! real to complex walkers and v.v.
            tgt_cpt = rotate_part(part_type)
            walkerweight = sign(1.0_dp,RealwSign(part_type))
            MatEl = real(rh,dp)

            ! spawn from real-parent to imaginary child: no sign change
            ! from imaginary to real -> sign change
            if (part_type == 2) walkerweight = -walkerweight

            nSpawn = - tau * MatEl * walkerweight / prob

            ! n.b. if we ever end up with |walkerweight| /= 1, then this
            !      will need to ffed further through.
            if (tSearchTau .and. (.not. tFillingStochRDMonFly)) &
                call log_spawn_magnitude (ic, ex, matel, prob)

            ! Keep track of the biggest spawn this cycle
            max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)
            
            if (tRealSpawning) then
                ! Continuous spawning. Add in acceptance probabilities.
                
                if (tRealSpawnCutoff .and. &
                    abs(nSpawn) < RealSpawnCutoff) then
                    p_spawn_rdmfac=abs(nSpawn)/RealSpawnCutoff
                    nSpawn = RealSpawnCutoff &
                           * stochastic_round (nSpawn / RealSpawnCutoff)
               else
                    p_spawn_rdmfac=1.0_dp !The acceptance probability of some kind of child was equal to 1
               endif
            else
                if(abs(nSpawn).ge.1) then
                    p_spawn_rdmfac=1.0_dp !We were certain to create a child here.
                    ! This is the special case whereby if P_spawn(j | i) > 1, 
                    ! then we will definitely spawn from i->j.
                    ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                    ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                    ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                else
                    p_spawn_rdmfac=abs(nSpawn)
                endif
                
                ! How many children should we spawn?

                ! And round this to an integer in the usual way
                ! HACK: To use the same number of random numbers for the tests.
                nSpawn = real(stochastic_round (nSpawn), dp)
                
            endif
            ! And create the parcticles
            child(tgt_cpt) = nSpawn

        else

            ! We actually want to calculate Hji - take the complex conjugate, 
            ! rather than swap around DetCurr and nJ.
#ifdef __REALTIME
            rh_used = conjg(rh)
#endif

            ! have to loop over the tgt_cpt similar to the complex implo
            ! if the Hamiltonian has real and imaginary components do it 
            ! similarily to complex implementation with H <-> J switched
            ! change this (lenof/inum) below .. todo
            do tgt_cpt = 1, (lenof_sign/inum_runs)
                component = tgt_cpt
                if (part_type == 2 .and. inum_runs == 1) component = 3 - tgt_cpt

                walkerweight = sign(1.0_dp,RealwSign(part_type)) 
#ifdef __REALTIME
                if (component == 1) then
                    MatEl = real(aimag(rh_used),dp)
                else 
                    MatEl = real(rh_used,dp)
                    if (part_type == 2) walkerweight = -walkerweight
                end if
#endif

                nSpawn = - tau * MatEl * walkerweight / prob

                ! n.b. if we ever end up with |walkerweight| /= 1, then this
                !      will need to ffed further through.
                if (tSearchTau .and. (.not. tFillingStochRDMonFly)) &
                    call log_spawn_magnitude (ic, ex, matel, prob)

                ! Keep track of the biggest spawn this cycle
                max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)
                
                if (tRealSpawning) then
                    ! Continuous spawning. Add in acceptance probabilities.
                    
                    if (tRealSpawnCutoff .and. &
                        abs(nSpawn) < RealSpawnCutoff) then
                        p_spawn_rdmfac=abs(nSpawn)/RealSpawnCutoff
                        nSpawn = RealSpawnCutoff &
                               * stochastic_round (nSpawn / RealSpawnCutoff)
                   else
                        p_spawn_rdmfac=1.0_dp !The acceptance probability of some kind of child was equal to 1
                   endif
                else
                    if(abs(nSpawn).ge.1) then
                        p_spawn_rdmfac=1.0_dp !We were certain to create a child here.
                        ! This is the special case whereby if P_spawn(j | i) > 1, 
                        ! then we will definitely spawn from i->j.
                        ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                        ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                        ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                    else
                        p_spawn_rdmfac=abs(nSpawn)
                    endif
                    
                    ! How many children should we spawn?

                    ! And round this to an integer in the usual way
                    ! HACK: To use the same number of random numbers for the tests.
                    nSpawn = real(stochastic_round (nSpawn), dp)
                    
                endif
                ! And create the parcticles
                child(tgt_cpt) = nSpawn
            end do
        end if
       
        if(tFillingStochRDMonFly) then
            if (child(part_type).ne.0) then
                !Only add in contributions for spawning events within population 1
                !(Otherwise it becomes tricky in annihilation as spawnedparents doesn't tell you which population
                !the event came from at present)
                call calc_rdmbiasfac(p_spawn_rdmfac, prob, realwSign(part_type), RDMBiasFacCurr) 
            else
                RDMBiasFacCurr = 0.0_dp
            endif
        else
            ! Not filling the RDM stochastically, bias is zero.
            RDMBiasFacCurr = 0.0_dp
        endif

        ! Avoid compiler warnings
        iUnused = walkExcitLevel

    end function attempt_create_realtime

    subroutine save_current_dets() 
        ! routine to copy the currentDets array and all the associated 
        ! pointers an hashtable related quantities to the 2nd temporary 
        ! list, from which the first spawn and y(n) + k1/2 addition is done 
        ! and the k2 spawing list is created to then use CurrentDets to go to
        ! the next time-step y(n+1) = y(n) + k2
        ! new idea to reduce the amount of rewriting of old routines:
        ! save the CurrentDets and all associated quantities in the temporary
        ! variables, then normally work on CurrentDets to create the k1 
        ! spawning and also the CurrentDets = y(n) + k1 / 2 combination
        ! then spawn again from this list to create k2 spawning array
        ! but then before combing reload the saved CurrentDets from the 
        ! temporary variable
        ! probably have to think about parallelism issues here..
        character(*), parameter :: this_routine = "save_current_dets"

        ! save the WalkVecDets variable, i think thats the only necessary 
        ! variable, the pointers don't count 
        temp_det_list = WalkVecDets
        
        ! for now also store the pointer, but thats not needed i guess
        temp_det_pointer => temp_det_list

        ! also need the hash table and the freeslot i guess
        ! cannot just copy the hash like that .. have to loop over it 
        ! and initialize it correctly.. that takes some effort
        ! do not actually have to do that here... just in the reload! 
        ! there i have to reassign HashIndex to the stored det-list
!         call clear_hash_table(temp_det_hash) 
!         call fill_in_hash_table(temp_det_hash, nWalkerHashes, CurrentDets, &
!             int(TotWalkers,sizeof_int), .true.)

!         call copy_hash_table(HashIndex, temp_det_hash)
!         temp_det_hash = HashIndex

        ! and the freeslot.. although this one gets reinitialized to 0 
        ! every iteration or not? yeah it is.. so i only have to reset it 
        ! twice in the rt-fciqmc before the y(n) + k2 combination
        ! do that in the reload_current_dets routine!
        ! same with n_determ_states var.
        
        ! also have to save current number of determinants! (maybe totparts too?)
        temp_totWalkers = TotWalkers

    end subroutine save_current_dets

!     subroutine copy_hash_table()
!         ! routine to copy the hash-table to a temporary 
! 
!     end subroutine copy_hash_table

    subroutine reload_current_dets()
        ! routine to reload the saved y(n) CurrentDets array for the final
        ! y(n) + k/2 combination to move to the next time step y(n+1)
        ! have also to think about the death-step and annihilation step 
        ! influences.. maybe have to write new death/born routines to split 
        ! that from the spawned list creation..
        character(*), parameter :: this_routine = "reload_current_dets" 

        ! copy the list
        WalkVecDets = temp_det_list

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

!         call copy_hash_table(temp_det_hash, HashIndex)
!         HashIndex = temp_det_hash

        ! for correct load Balancin i also have to reset the the freeslot var.
        ! here i have to reset the freeslot values to the values after the 
        ! first spawn loop, as the empty entries get determined there! 
        ! and i only want to do an additional annihilation step to combine 
        ! y(n) + k2
        FreeSlot = temp_freeslot
        ! setting this to one is enough
        iStartFreeSlot = 1
        iEndFreeSlot = temp_iendfreeslot
        ! not sure if i have to reset this here, or in the routine to reset 
        ! the spawned list
        ! i think that has to be done in the reset_spawned_list
!         n_determ_states = 1

#ifdef __DEBUG
        call reset_tot_parts()
#endif        

    end subroutine reload_current_dets

    subroutine reset_spawned_list()
        ! also need a routine to reset the spawned lists before the second 
        ! spawning step for the 2nd order RK method in the rt-fciqmc
!         integer, intent(out) :: n_determ_states
        character(*), parameter :: this_routine = "reset_spawned_list"

        ! Reset positions to spawn into in the spawning array.
        ValidSpawnedList = InitialSpawnedSlots

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        ! think i have to reset the deterministic counter here
!         n_determ_states = 1

        ! also reset the diagonal specific valid spawn list.. i think i 
        ! can just reuse the InitialSpawnedSlots also
        valid_diag_spawn_list = InitialSpawnedSlots

        print *, "valid_diag_spawn_list", valid_diag_spawn_list

        ! also save the number of particles from this spawning to calc. 
        ! first step specific acceptance rate
        ! changed that! this has to be done BEFORE the annihilation step! 
!         SumWalkersCyc_1 = SumWalkersCyc_1 + sum(TotParts)

    end subroutine reset_spawned_list

    subroutine setup_temp_det_list()
        ! setup the second list to temporaly store the list of determinants
        ! necessary in the real-time fciqmc list
        ! determine the necessary size from the already setup CurrentDets
        character(*), parameter :: this_routine = "setup_temp_det_list"
        integer :: ierr, tmp_siz1, tmp_siz2, i, spawn_ht_mem

        tmp_siz1 = size(WalkVecDets,dim=1)
        tmp_siz2 = size(WalkVecDets,dim=2)

        ! allocate the array
        allocate(temp_det_list(0:tmp_siz1-1,tmp_siz2), stat = ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")

        ! and init it
        temp_det_list(0:tmp_siz1-1,1:tmp_siz2) = 0
        
        ! and point to it 
        temp_det_pointer => temp_det_list

        ! and also allocate the hash-table
        tmp_siz1 = size(HashIndex)
        allocate(temp_det_hash(tmp_siz1), stat = ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")

        ! and initialize it to 0
        do i = 1, tmp_siz1
            temp_det_hash(i)%ind = 0
        end do

        ! also need a corresponding freeslot array
        tmp_siz1 = size(freeslot)
        allocate(temp_freeslot(tmp_siz1), stat = ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")

        ! and init it
        temp_freeslot = 0

        ! also use the spawn_ht hash table, so also allocate it here! 

        ! Allocate the hash table to the spawning array.
        ! The number of MB of memory required to allocate spawn_ht.
        ! Each node requires 16 bytes.
        nhashes_spawn = int(0.8_dp * real(MaxSpawned,dp))
        spawn_ht_mem = nhashes_spawn*16/1000000
        write(6,'(a78,'//int_fmt(spawn_ht_mem,1)//')') "About to allocate hash table to the spawning array. &
                                       &Memory required (MB):", spawn_ht_mem
        write(6,'(a13)',advance='no') "Allocating..."; call neci_flush(6)
        allocate(spawn_ht(nhashes_spawn), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(this_routine, "Error allocating spawn_ht array.")
        else
            write(6,'(1x,a5)') "Done."
            write(6,'(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                              &increase as further nodes are added."
        end if

        call init_hash_table(spawn_ht)

    end subroutine setup_temp_det_list

    subroutine new_child_stats_realtime(iter_data, ilutI, nJ, ilutJ, ic, &
                walkExLevel, child, parent_flags, part_type) 
        ! new realtime fciqmc specific spawn statistics gathering routine 
        ! due to the different algorithm in the real-time fciqmc, with 2 
        ! application of the Hamiltonian in one cycle, the bookkeeping 
        ! function also have to be changed majorly! 
        ! the change is not so much, mainly that the particles spawn to the 
        ! other type of particles (Re <-> Im) all the time if a pure real
        ! Hamiltonian is used. and both if the Hamiltonian has complex entries
        ! but i guess for the complex run this is also wrongly implemented 
        ! in the old code, as this is not considered..
        ! NO! its correctly considered in the old routine.. doh! 
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel) 
        integer, intent(in) :: part_type
        real(dp), intent(in) :: child(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "new_child_stats_realtime"

!         NoBorn(1) = NoBorn(1) + sum(abs(child))
!         if (ic == 1) SpawnFromSing(1) = SpawnFromSing(1) + sum(abs(child))

        !... nah its correct in the original routine... duh..
    end subroutine new_child_stats_realtime

    ! subroutine to calculate the overlap of the current y(t) = a_j(a^+_j)(t)y(0)>
    ! time evolved wavefunction to the saved <y(0)|a^+_i(a_i) 
    subroutine update_gf_overlap() 
        ! this routine only deals with globally defined variables
        integer :: idet, nI(nel), det_ind, hash_val, i
        real(dp) :: overlap(lenof_sign), real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        logical :: tDetFound

        overlap = 0.0_dp

        do idet = 1, size(perturbed_ground, dim = 2)

            call extract_sign(perturbed_ground(:,idet), real_sign_1) 

            ! why should i store zero valued dets in the perturbed_gs? 
            ! and i have to expand that to complex code later on..
            ! although.. the perturbed ground state wavefunction is always
            ! real-valued -> so i only need the actual real-part of the 
            ! time evolved wavefunction! 
            if (IsUnoccDet(real_sign_1)) cycle

            call decode_bit_det(nI, perturbed_ground(:,idet))

            ! search for the hash table associated with the time evolved 
            ! wavefunction -> is this already initialized correctly? 
            call hash_table_lookup(nI, perturbed_ground(:,idet), nifdbo, &
                HashIndex, CurrentDets, det_ind, hash_val, tDetFound)

            if (tDetFound) then
                ! since the original gound state is pure real... (atleast thats
                ! what we assume for now -> we only need the real part of 
                ! the time-evolved wf
                call extract_sign(CurrentDets(:,det_ind), real_sign_2)

                do i = 1, lenof_sign
                    overlap(i) = overlap(i) + real_sign_1(1) * real_sign_2(i)
                end do
            end if
        end do

        ! need the timestep here... or the cycle of the current real-time loop
        gf_overlap(:,iter) = overlap 


    end subroutine update_gf_overlap

    function calc_perturbed_norm() result(pert_norm) 
        ! function to calculate the norm of the left-hand <y(0)|
        real(dp) :: pert_norm
        character(*), parameter :: this_routine = "calc_perturbed_norm"

        integer :: idet 
        real(dp) :: tmp_sign(lenof_sign)

        pert_norm = 0.0_dp

        do idet = 1, TotWalkers_pert

            call extract_sign(perturbed_ground(:,idet), tmp_sign)

            ! for now assume real-only groundstates from which we start 
            pert_norm = pert_norm + tmp_sign(1)**2

        end do

    end function calc_perturbed_norm

    subroutine create_perturbed_ground()
        ! routine to create from the already read in popsfile info in 
        ! popsfile_dets the left hand <y(0)| by applying the corresponding 
        ! creation or annihilation operator
        character(*), parameter :: this_routine = "create_perturbed_ground"
        integer :: tmp_totwalkers
        integer :: ierr

        tmp_totwalkers = TotWalkers_orig

        print *, "Creating the wavefunction to projected on!"
        print *, "Initial number of walkers: ", tmp_totwalkers


        allocate(perturbed_ground(0:niftot,TotWalkers_orig), stat = ierr)
        call apply_perturbation(overlap_pert(1), tmp_totwalkers, popsfile_dets,&
            perturbed_ground)
        TotWalkers_pert = int(tmp_totwalkers, int64)

        print *, "Walkers remaining in perturbed ground state:" , TotWalkers_pert

        ! also need to create and associated hash table to this list 
!         call clear_hash_table(perturbed_ground_ht)
        ! or maybe not... lets see later on..


    end subroutine create_perturbed_ground
    
    subroutine check_update_growth(iter_data, message)
      use Parallel_neci, only : iProcIndex, MPISumAll, root
      use spin_project, only : tSpinProject
      use real_time_data, only : TotPartsStorage
      implicit none
      character(len=*), intent(in) :: message
      type(fcimc_iter_data), intent(in) :: iter_data
      real(dp) :: growth(lenof_sign), growth_tot(lenof_sign)
      real(dp) :: allWalkers(lenof_sign), allWalkersOld(lenof_sign)
      growth = iter_data%nborn &
         - iter_data%ndied - iter_data%nannihil &
         - iter_data%naborted - iter_data%nremoved

      call MPISumAll(growth,growth_tot)
      call MPISumAll(TotParts,allWalkers)
      call MPIsumAll(TotPartsStorage,allWalkersOld)
      TotPartsStorage = TotParts
         write(iout,*) "update_growth: ", growth_tot
         write(iout,*) "AllTotParts: ", allWalkers
         write(iout,*) "AllTotPartsOld: ", allWalkersOld
      if((iProcIndex == root) .and. .not. tSpinProject .and. &
           any(abs(growth_tot - (allWalkers - allWalkersOld)) > 1.0e-5_dp)) then
         write(iout,*) message

         call stop_all("check_update_growth", &
              "Assertation failed: all(iter_data_fciqmc%update_growth_tot.eq.AllTotParts_1-AllTotPartsOld_1)")
      end if

    end subroutine check_update_growth

    subroutine reset_tot_parts()
      ! if the second RK step is to be compared, the reference has to be reset
      ! -> recount the TotParts from the restored data
      use real_time_data, only : TotPartsStorage
      implicit none
      integer :: i
      real(dp) :: CurrentSign(lenof_sign)
      TotParts = 0.0_dp
      do i=1, TotWalkers
         call extract_sign(CurrentDets(:,i),CurrentSign)
         TotParts = TotParts + abs(CurrentSign)
      end do
      TotPartsStorage = TotParts
    end subroutine reset_tot_parts
end module real_time_procs
