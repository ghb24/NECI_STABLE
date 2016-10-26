#include "macros.h"

! module containing useful functions and subroutines needed in the real-time
! implementation of the FCIQMC algotrithm

module real_time_procs
    use hash, only: hash_table_lookup, init_hash_table, clear_hash_table, &
                    add_hash_table_entry, fill_in_hash_table
    use SystemData, only: nel
    use real_time_data, only: gf_overlap, TotWalkers_orig, TotWalkers_pert, &
                              t_complex_ints, real_time_info, temp_freeslot, & 
                              temp_det_list, temp_det_pointer,  temp_iendfreeslot, &
                              temp_det_hash, temp_totWalkers, pert_norm, &
                              valid_diag_spawn_list, DiagParts, n_diag_spawned, &
                              NoDied_1, NoBorn_1, SumWalkersCyc_1, &
                              t_rotated_time, wf_norm, tau_imag, tau_real, gs_energy, &
                              dyn_norm_psi, shift_damping, MaxSpawnedDiag, & 
                              InitialSpawnedSlotsDiag
    use kp_fciqmc_data_mod, only: perturbed_ground, overlap_pert
    use constants, only: dp, lenof_sign, int64, n_int, EPS, iout, null_part, &
                         sizeof_int, MPIArg
    use bit_reps, only: decode_bit_det, test_flag, encode_sign, &
                        set_flag, encode_bit_rep, extract_bit_rep, &
                        flag_has_been_initiator, flag_deterministic, encode_part_sign, &
                        nullify_ilut_part, get_initiator_flag, get_initiator_flag_by_run

                        
    use bit_rep_data, only: extract_sign, nifdbo, niftot
    use FciMCData, only: CurrentDets, HashIndex, popsfile_dets, MaxWalkersPart, &
                         WalkVecDets, freeslot, spawn_ht, nhashes_spawn, MaxSpawned, &
                         iStartFreeSlot, iEndFreeSlot, ValidSpawnedList, &
                         InitialSpawnedSlots, iLutRef, inum_runs, max_cyc_spawn, &
                         tSearchTau, tFillingStochRDMonFly, fcimc_iter_data, &
                         NoAddedInitiators, SpawnedParts, acceptances, TotWalkers, &
                         nWalkerHashes, iter, fcimc_excit_gen_store, NoDied, &
                         NoBorn, NoAborted, NoRemoved, HolesInList, TotParts, Hii, &
                         determ_sizes, determ_displs, determ_space_size, core_space
    use sparse_arrays, only: sparse_core_ham
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
    use ParallelHelper, only: nNodes, bNodeRoot, ProcNode, NodeRoots, MPIBarrier, &
         iProcIndex, MPI_SUM, root
    use Parallel_neci
    use LoggingData, only: tNoNewRDMContrib
    use AnnihilationMod, only: test_abort_spawn
    use load_balance, only: AddNewHashDet, CalcHashTableStats, test_hash_table

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

        ! valid_diag_spawn_list gets increased after the spawn
        ! to DiagParts(:,valid_diag_spawn_list) -> highest index is
        ! actually valid_diag_spawn_list-1 

        numSpawns = valid_diag_spawn_list(iProcIndex)-1
        call AnnihilateDiagParts(numSpawns, TotWalkersNew, iter_data)

        ! also should update the hashtable stats, specific for this diagonal 
        ! spawning event, but the original one should work also for this 
        ! since it only takes CurrentDets into account! 
        
        call CalcHashTableStats(TotWalkersNew,iter_data)
        
        ! this should be it, deterministic annihilation is carried out in the next
        ! step, within the 'regular' annihilation
        
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

        do i = InitialSpawnedSlotsDiag(iProcIndex), ValidSpawned

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
            
            if (tSuccess) then

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
                        if (is_run_unnocc(CurrentSign,run)) then
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

               ! Running the full, non-initiator scheme.
               ! Determinant in newly spawned list is not found in
               ! CurrentDets. If coeff <1, apply removal criterion.
               call extract_sign (DiagParts(:,i), SignTemp)

               if (.not. IsUnoccDet(SignTemp)) then
                  ! Walkers have not been aborted and so we should copy the
                  ! determinant straight over to the main list. We do not
                  ! need to recompute the hash, since this should be the
                  ! same one as was generated at the beginning of the loop.
                  ! also here treat those new walkers as born particles

                  iter_data%nborn = iter_data%nborn + abs(SignTemp)
                  do run=1, inum_runs
                     NoBorn(run) = NoBorn(run) + sum(abs(SignTemp(&
                          min_part_type(run):max_part_type(run))))
                  enddo
                  call AddNewHashDet(TotWalkersNew, DiagParts(:,i), DetHash, nJ)

               end if
            end if

            !                if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) then
            ! We must use the instantaneous value for the off-diagonal contribution.
            !                call check_fillRDM_DiDj(rdms, i, SpawnedParts(0:NifTot,i), SignTemp)
            !                end if 

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
            if (valid_diag_spawn_list(proc) > MaxSpawnedDiag) list_full = .true.
        else
            if (valid_diag_spawn_list(proc) >= InitialSpawnedSlotsDiag(proc+1)) &
                list_full=.true.
        end if
        if (list_full) then
            write(6,*) "Attempting to spawn particle onto processor: ", proc
            write(6,*) "No memory slots available for this spawn."
            write(6,*) "Please increase MEMORYFACSPAWN"
            write(6,*) "VALID DIAG SPAWN LIST", valid_diag_spawn_list
            write(6,*) "NUMBER OF OCCUPIED DETS", TotWalkers
            call stop_all(this_routine, "Out of memory for spawned particles")
        end if

        call encode_bit_rep(DiagParts(:, valid_diag_spawn_list(proc)), ilut, &
                            diag_sign, flags)

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

        valid_diag_spawn_list(proc) = valid_diag_spawn_list(proc) + 1
        
        ! Sum the number of created children to use in acceptance ratio.
        ! i dont need to do acceptances here, since diagonal events are not 
        ! counted in the acceptance ratio..
!         acceptances(1) = acceptances(1) + sum(abs(diag_sign))

    end subroutine create_diagonal_as_spawn

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
           fac(min_part_type(run)) = tau_real * (Kii + Hii - gs_energy(run) )
           fac(max_part_type(run)) = 0.0_dp
        enddo

        if (real_time_info%damping < EPS .and. .not. t_rotated_time) then

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
                        do run = 1, lenof_sign
                           fac(run) = min(2.0_dp, fac(run))
                        end  do
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
                  ndie(min_part_type(run)) =  - fac(min_part_type(run)) &
                       * realwSign(max_part_type(run))
                  ! and - from Re -> Im
                  ! does this give the correct sign compared to the parent sign?
                  ! additional -1 is added in postprocessing to convert ndie -> nborn
                  ndie(max_part_type(run)) = fac(min_part_type(run)) &
                       * realwSign(min_part_type(run))

               else 
                  ! if not exact i have to round stochastically
                  ! Im -> Re
                  rat =  - fac(min_part_type(run)) * RealwSign(max_part_type(run))

                  ndie(min_part_type(run)) = real(int(rat), dp) 
                  rat = rat - ndie(min_part_type(run))

                  r = genrand_real2_dSFMT() 
                  if (abs(rat) > r) ndie(min_part_type(run)) = &
                       ndie(min_part_type(run)) + real(nint(sign(1.0_dp,rat)),dp)

                  ! Re -> Im
                  rat =   fac(min_part_type(run)) * RealwSign(min_part_type(run))
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
            ! the diagnonal matrix elements are always real, so there is
            ! no contribution from tau_real via Kii (and no influence of
            ! tau_imag on fac(min_part_type(run))
           do run = 1, inum_runs
              ! tau_real and tau_imag have opposite signs -> shift changed sign
              ! when moved from imaginary to real part
              fac(max_part_type(run)) =  tau_real * (real_time_info%damping) &
                  + tau_imag*(Kii + Hii - gs_energy(run) - DiagSft(run) )
           enddo

            ! and also about the fac restrictions.. for now but it here anyway..
            if(any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTau) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do run = 1, 2
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
                  ndie(min_part_type(run)) =  - fac(min_part_type(run)) &
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
                       + real(nint(sign(1.0_dp,rat)),dp)

                  !  -> Im
                  rat = fac(min_part_type(run)) * RealwSign(min_part_type(run)) &
                  - fac(max_part_type(run)) * RealwSign(max_part_type(run))

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

        ndie = attempt_die_realtime(DetCurr, Kii, realwSign, walkExcitLevel)
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
            ! rotated_time_setup: there is only one initiator flag for
            ! both complex and real walkers -> initiator flag is kept
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
               do j = 1, inum_runs
                  if (test_flag(iLutCurr,get_initiator_flag_by_run(j))) then
                     NoAddedInitiators(j) = NoAddedInitiators(j) - 1
                  end if
               enddo
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

        real(dp) :: rat, r, walkerweight, pSpawn, nSpawn, MatEl, p_spawn_rdmfac, &
             sepSign
        integer :: extracreate, tgt_cpt, component, i, iUnused
        integer :: TargetExcitLevel
        logical :: tRealSpawning
        HElement_t(dp) :: rh, rh_used

        ! This is crucial
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
!        print *, "Spawn event with mat el", rh, "sign", RealwSign, &
!             "timestep", tau, "probability", prob
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
        if (.not. t_complex_ints .and. .not. t_rotated_time) then
            ! if it is a pure real-hamiltonian there is only spawing from 
            ! real to complex walkers and v.v.
            tgt_cpt = rotate_part(part_type)
            walkerweight = sign(1.0_dp,RealwSign(part_type))
            MatEl = real(rh,dp)

            ! spawn from real-parent to imaginary child: no sign change
            ! from imaginary to real -> sign change
            if (mod(part_type,2) == 0) walkerweight = -walkerweight

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

            ! have to loop over the tgt_cpt similar to the complex impl
            ! if the Hamiltonian has real and imaginary components do it 
            ! similarily to complex implementation with H <-> J switched
            ! rmneci_setup: adjusted for multirun, fixed complex -> real spawns
            do component = 1, (lenof_sign/inum_runs)
               tgt_cpt = min_part_type(part_type_to_run(part_type)) - 1 + component
                ! keep track of the sign due to the kind of spawn event
                sepSign = 1.0_dp
                ! if (part_type == 2 .and. inum_runs == 1) component = 3 - tgt_cpt !?

                walkerweight = sign(1.0_dp,RealwSign(part_type)) 
                if (mod(part_type,2) == 0 .and. component == 1) &
                     sepSign = (-1.0_dp)
#ifdef __REALTIME
                ! part_type is given as input, for that part_type, the real part of
                ! the HElement is used if rotation occurs and the imaginary part if not
                if (mod(component,2) == mod(part_type,2)) then
                   ! spawn part_type -> part_type
                    MatEl = - real(aimag(rh_used),dp)*tau_real - real(rh_used,dp)*tau_imag
                else 
                   ! spawn part_type -> rotate_part(part_type)
                    MatEl = real(rh_used,dp)*tau_real - real(aimag(rh_used),dp)*tau_imag
                end if
#endif
                nSpawn = - sepSign * MatEl * walkerweight / prob

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
                ! And create the parcticles (in the correct run)
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
      use real_time_data, only: TotPartsStorage
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
        
        ! And save the old TotParts value, as this might have changed and iter_data is reset
        ! (some weird scenario in which CalcHashTableStats is called at the end of the 
        ! time-step and and then modifies both TotParts and iter_data correctly, but iter_data
        ! is reset at the beginning of the iteration, so TotParts also has to)
        TotPartsStorage = TotParts

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
        ! also reload the positions of empty slots in the ensemble
        iStartFreeSlot = 1
        iEndFreeSlot = temp_iendfreeslot
        FreeSlot = temp_freeslot

        ! i think that has to be done in the reset_spawned_list
!         n_determ_states = 1

        call reset_tot_parts()

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
        valid_diag_spawn_list = InitialSpawnedSlotsDiag

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

    ! subroutine to calculate the overlap of the current y(t) = a_j(a^+_j)(t)y(0)>
    ! time evolved wavefunction to the saved <y(0)|a^+_i(a_i) 
    subroutine update_gf_overlap() 
        ! this routine only deals with globally defined variables
      implicit none
        integer :: idet, nI(nel), det_ind, hash_val, i, run
        real(dp) :: overlap(lenof_sign), real_sign_1(lenof_sign), real_sign_2(lenof_sign), &
             tmp_norm(inum_runs), norm_buf(inum_runs)
        logical :: tDetFound

        overlap = 0.0_dp
        tmp_norm = 0.0_dp

        do idet = 1, size(perturbed_ground, dim = 2)

            call extract_sign(perturbed_ground(:,idet), real_sign_1) 

            if (IsUnoccDet(real_sign_1)) cycle

            call decode_bit_det(nI, perturbed_ground(:,idet))

            ! search for the hash table associated with the time evolved 
            ! wavefunction -> is this already initialized correctly? 
            call hash_table_lookup(nI, perturbed_ground(:,idet), nifdbo, &
                HashIndex, CurrentDets, det_ind, hash_val, tDetFound)

            if (tDetFound) then
               ! both real and imaginary part of the time-evolved wf are required
                call extract_sign(CurrentDets(:,det_ind), real_sign_2)
                
                do i = 1, lenof_sign
                   run = part_type_to_run(i)
                   if(mod(i,2)==0) then
                      ! imaginary part of the overlap
                    overlap(i) = overlap(i) + &
                         real_sign_1(min_part_type(run)) * real_sign_2(max_part_type(run)) - &
                         real_sign_2(min_part_type(run)) * real_sign_1(max_part_type(run))
                    else 
                       ! real part of the overlap
                       overlap(i) = overlap(i) + &
                         real_sign_1(min_part_type(run)) * real_sign_2(min_part_type(run)) + &
                         real_sign_2(max_part_type(run)) * real_sign_1(max_part_type(run)) 
                    endif
                    tmp_norm(run) = tmp_norm(run) + real_sign_2(i)**2
                end do
            end if
        end do

        ! rmneci_setup: the overlap has to be reduced as each proc
        ! only computes its own part
        call MPIReduce(overlap,MPI_SUM,gf_overlap(:,iter))
        ! communicate the norm as it is the sum over all walkers
        call MPIReduce(tmp_norm,MPI_SUM,norm_buf)
        ! also store the new norm information 
        do run = 1, inum_runs
           wf_norm(run,iter) = sqrt(norm_buf(run) * pert_norm(run))
        end do

    end subroutine update_gf_overlap

    function calc_norm(dets, num_dets) result(cd_norm) 
      ! the first dimension of dets has to be lenof_sign
        ! function to calculate the norm of the left-hand <y(0)| (general function)
        real(dp) :: cd_norm(inum_runs)
        integer(dp) :: dets(:,:)
        integer, intent(in) :: num_dets
        character(*), parameter :: this_routine = "calc_perturbed_norm"

        integer :: idet, run
        real(dp) :: tmp_sign(lenof_sign)

        cd_norm = 0.0_dp


        do idet = 1, num_dets

           call extract_sign(dets(:,idet), tmp_sign)
           do run = 1,inum_runs

              cd_norm(run) = cd_norm(run) + tmp_sign(min_part_type(run))**2 &
                   + tmp_sign(max_part_type(run))**2

           end do
        end do 

    end function calc_norm

    subroutine real_time_determ_projection()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use FciMCData, only: partial_determ_vecs, full_determ_vecs, SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV

        integer :: i, j, ierr, run, part_type

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_determ_vecs, full_determ_vecs, &
                            determ_sizes, determ_displs)

        call halt_timer(SemiStoch_Comms_Time)

        call set_timer(SemiStoch_Multiply_Time)

        if (determ_sizes(iProcIndex) >= 1) then

            ! For the moment, we're only adding in these contributions when we need the energy
            ! This will need refinement if we want to continue with the option of inst vs true full RDMs
            ! (as in another CMO branch).

            ! Perform the multiplication.

            partial_determ_vecs = 0.0_dp

            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    do run = 1, inum_runs
                       ! real part of the 'spawn'
                       partial_determ_vecs(min_part_type(run),i) = &
                            partial_determ_vecs(min_part_type(run),i) + &

                            tau_real*Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            min_part_type(run),sparse_core_ham(i)%positions(j)) +&

                            tau_real*Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            max_part_type(run),sparse_core_ham(i)%positions(j)) +&
                            
                            tau_imag*Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            min_part_type(run),sparse_core_ham(i)%positions(j)) -&

                            tau_imag*Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            max_part_type(run),sparse_core_ham(i)%positions(j))

                       ! imaginary part
                       partial_determ_vecs(max_part_type(run),i) = &
                            partial_determ_vecs(max_part_type(run),i) - &

                            tau_real*Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            min_part_type(run),sparse_core_ham(i)%positions(j)) +&

                            tau_real*Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            max_part_type(run),sparse_core_ham(i)%positions(j)) +&

                            tau_imag*Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            max_part_type(run),sparse_core_ham(i)%positions(j)) +&

                            tau_imag*Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(&
                            min_part_type(run),sparse_core_ham(i)%positions(j))
                        
                    end do
                end do
            end do

            ! Now add shift*full_determ_vecs to account for the shift, not stored in
            ! sparse_core_ham.
            do i = 1, determ_sizes(iProcIndex)
                do run  = 1, inum_runs 
                   ! real part
                   partial_determ_vecs(min_part_type(run),i) = &
                        partial_determ_vecs(min_part_type(run),i) + &
                        (Hii - gs_energy(run)) * full_determ_vecs(max_part_type(run),i + &
                        determ_displs(iProcIndex)) * &
                        tau_real + tau_imag * full_determ_vecs(min_part_type(run),i + &
                        determ_displs(iProcIndex)) * (Hii - gs_energy(run) - DiagSft(run))

                   ! imaginary part
                   partial_determ_vecs(max_part_type(run),i) = &
                        partial_determ_vecs(max_part_type(run),i) + tau_imag * &
                        (Hii - gs_energy(run) - DiagSft(run)) * full_determ_vecs(& 
                        max_part_type(run),i + determ_displs(iProcIndex)) - tau_real &
                        * (Hii - gs_energy(run)) * full_determ_vecs(min_part_type(run),i &
                        + determ_displs(iProcIndex)) 
                enddo
            end do
        end if

        call halt_timer(SemiStoch_Multiply_Time)

    end subroutine real_time_determ_projection

    subroutine refresh_semistochastic_space()
      use CalcData, only: ss_space_in
      use semi_stoch_gen, only: init_semi_stochastic
      use semi_stoch_procs, only: end_semistoch
      implicit none
      ! as the important determinants might change during time evolution, this
      ! resets the semistochastic space taking the current population to get it
      call end_semistoch()
      call init_semi_stochastic(ss_space_in)
    end subroutine refresh_semistochastic_space

    subroutine create_perturbed_ground()
        ! routine to create from the already read in popsfile info in 
        ! popsfile_dets the left hand <y(0)| by applying the corresponding 
        ! creation or annihilation operator
      implicit none
        character(*), parameter :: this_routine = "create_perturbed_ground"
        integer :: tmp_totwalkers
        integer(int64) :: TotWalkers_orig_max
        integer :: ierr

        tmp_totwalkers = TotWalkers_orig

        print *, "Creating the wavefunction to projected on!"
        print *, "Initial number of walkers: ", tmp_totwalkers

        call MPISumAll(TotWalkers_orig,TotWalkers_orig_max)
        allocate(perturbed_ground(0:niftot,TotWalkers_orig_max), stat = ierr)
        print *, "Read-in dets", TotWalkers_orig_max
        if(allocated(overlap_pert)) then
           call apply_perturbation(overlap_pert(1), tmp_totwalkers, popsfile_dets,&
                perturbed_ground)
        else
           perturbed_ground = CurrentDets
        endif
        TotWalkers_pert = int(tmp_totwalkers, int64)

        print *, "Walkers remaining in perturbed ground state:" , TotWalkers_pert

        ! also need to create and associated hash table to this list 
!         call clear_hash_table(perturbed_ground_ht)
        ! or maybe not... lets see later on..


    end subroutine create_perturbed_ground

    
    subroutine check_update_growth(iter_data, message)
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
      if((iProcIndex == root) .and. .not. tSpinProject .and. &
           any(abs(growth_tot - (allWalkers - allWalkersOld)) > 1.0e-5_dp)) then
         write(iout,*) message
         write(iout,*) "update_growth: ", growth_tot
         write(iout,*) "AllTotParts: ", allWalkers
         write(iout,*) "AllTotPartsOld: ", allWalkersOld
         call stop_all("check_update_growth", &
              "Assertation failed: all(iter_data_fciqmc%update_growth_tot.eq.AllTotParts_1-AllTotPartsOld_1)")
      end if

    end subroutine check_update_growth

    subroutine update_shift_damping()
      implicit none
! sign convention for imaginary and real time differ
      shift_damping = shift_damping + tau_imag * DiagSft

    end subroutine update_shift_damping

    subroutine reset_tot_parts()
      ! if the second RK step is to be compared, the reference has to be reset
      ! -> recount the TotParts from the restored data
      use real_time_data, only : TotPartsStorage
      implicit none
      TotParts = get_tot_parts()
      TotPartsStorage = TotParts
    end subroutine reset_tot_parts

    function get_tot_parts() result(allWalkersSummed)
      ! if the second RK step is to be compared, the reference has to be reset
      ! -> recount the TotParts from the restored data
      implicit none
      integer :: i
      real(dp) :: CurrentSign(lenof_sign)
      real(dp) :: allWalkersSummed(lenof_sign)
      allWalkersSummed = 0.0_dp
      do i=1, TotWalkers
         call extract_sign(CurrentDets(:,i),CurrentSign)
         allWalkersSummed = allWalkersSummed + abs(CurrentSign)
      end do
    end function get_tot_parts
end module real_time_procs
