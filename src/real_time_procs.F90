#include "macros.h"

! module containing useful functions and subroutines needed in the real-time
! implementation of the FCIQMC algotrithm

module real_time_procs
    use hash, only: hash_table_lookup, init_hash_table, clear_hash_table, &
                    add_hash_table_entry
    use SystemData, only: nel
    use real_time_data, only: gf_overlap, TotWalkers_orig, TotWalkers_pert, &
                              t_complex_ints, real_time_info, temp_freeslot, & 
                              temp_iendfreeslot, temp_det_list, temp_det_pointer, &
                              temp_det_hash

    use kp_fciqmc_data_mod, only: perturbed_ground, overlap_pert
    use constants, only: dp, lenof_sign, int64, n_int, EPS, iout, null_part
    use bit_reps, only: decode_bit_det, test_flag, flag_initiator, encode_sign, &
                        set_flag, encode_bit_rep
    use bit_rep_data, only: extract_sign, nifdbo, niftot
    use FciMCData, only: CurrentDets, HashIndex, popsfile_dets, MaxWalkersPart, &
                         WalkVecDets, freeslot, spawn_ht, nhashes_spawn, MaxSpawned, &
                         iStartFreeSlot, iEndFreeSlot, ValidSpawnedList, &
                         InitialSpawnedSlots, iLutRef, inum_runs, max_cyc_spawn, &
                         tSearchTau, tFillingStochRDMonFly, fcimc_iter_data, &
                         NoAddedInitiators, SpawnedParts, acceptances
    use perturbations, only: apply_perturbation
    use util_mod, only: int_fmt
    use CalcData, only: AvMCExcits, tAllRealCoeff, tRealCoeffByExcitLevel, &
                        tRealSpawnCutoff, RealSpawnCutoff, tau, RealCoeffExcitThresh, &
                        DiagSft, tTruncInitiator
    use DetBitOps, only: FindBitExcitLevel
    use procedure_pointers, only: get_spawn_helement
    use util_mod, only: stochastic_round
    use tau_search, only: log_spawn_magnitude
    use rdm_general, only: calc_rdmbiasfac
    use global_det_data, only: global_determinant_data
    use rdm_filling, only: det_removed_fill_diag_rdm
    use rdm_data, only: nrdms, rdms
    use hash, only: remove_hash_table_entry
    use dSFMT_interface, only: genrand_real2_dSFMT
    use load_balance_calcnodes, only: DetermineDetNode
    use ParallelHelper, only: nNodes

    implicit none

contains

    subroutine create_diagonal_as_spawn(nI, ilut, diag_sign, iter_data)
        ! routine to treat the diagonal death/cloning step in the 2nd 
        ! RK loop as spawned particles, with using hash table too
        ! have to write this routine, since i also need the old one
        ! in the "normal" spawnings of the rt-fciqmc
        ! only need the currently looped over determinant! and probably 
        ! also just have to copy all the sign(or just store the original 
        ! parent_ilut in the hash table
        ! it does a simultanious spawning of real and complex walkers! 
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:niftot)
        real(dp), intent(in) :: diag_sign(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = "create_diagonal_as_spawn"

        integer :: proc, ind, hash_val, i
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_old(lenof_sign), real_sign_new(lenof_sign)
        real(dp) :: sgn_prod(lenof_sign)
        logical :: list_full, tSuccess
        integer, parameter :: flags = 0       
        
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

            ! Set the initiator flags appropriately.
            ! If this determinant (on this replica) has already been spawned to
            ! then set the initiator flag. Also if this child was spawned from
            ! an initiator, set the initiator flag.
            ! in this routine i do a simultanious spawning of real and complex 
            ! walkers -> and here i just considere the diagonal "spawn" of 
            ! already occupied dets..
            if (tTruncInitiator) then
                do i = 1, lenof_sign
                    if (abs(real_sign_old(i)) > EPS .or. &
                        test_flag(ilut, flag_initiator(i))) then
                        call set_flag(SpawnedParts(:,ind), flag_initiator(i))
                    end if
                end do
            end if
        else
            ! Determine which processor the particle should end up on in the
            ! DirectAnnihilation algorithm.
            proc = DetermineDetNode(nel, nI, 0)

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
                do i = 1, lenof_sign
                    if (test_flag(ilut, flag_initiator(i))) then
                        call set_flag(SpawnedParts(:, ValidSpawnedList(proc)),&
                            flag_initiator(i))
                    end if
                end do
            end if

            call add_hash_table_entry(spawn_ht, ValidSpawnedList(proc), hash_val)

            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        end if
        
        ! Sum the number of created children to use in acceptance ratio.
        ! how is this done in the rt-fciqmc??
        acceptances(1) = acceptances(1) + sum(abs(diag_sign))

    end subroutine create_diagonal_as_spawn

    function attempt_die_realtime(DetCurr, RealwSign, ilutCurr, Kii, DetPosition, &
            iter_data, walkExcitLevel) result(ndie)
        ! also write a function, which calculates the new "signs"(weights) of 
        ! the real and complex walkers for the diagonal death/cloning step 
        ! since i need that for both 1st and 2nd loop of RK, but at different 
        ! points 
        integer, intent(in) :: DetCurr(nel) 
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        real(dp), intent(in) :: Kii
        integer, intent(in) :: DetPosition
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign) :: ndie
        real(dp), dimension(lenof_sign) :: CopySign
        integer, intent(in) :: walkExcitLevel
        integer :: i, irdm
        real(dp) :: fac(lenof_sign), rat, r

        character(*), parameter :: this_routine = "attempt_die_realtime"

        ! do it dependent on if there is damping, since otherwise i can save
        ! some effort..

        ! do i need a minus sign here?? just from the convention
        ! in the rest of the neci code? yes!
        fac(1) = -tau * (Kii - DiagSft(1))

        if (real_time_info%damping < EPS) then
            ! in this case there is only Re <-> Im 
            ! n_i' =  H_ii c_i 
            ! c_i' = -H_ii n_i
            ! how should i log death_magnitudes here.. for tau-search.. 
            ! that probably has to change too..
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
                        do i = 1, inum_runs
                            fac(i) = min(2.0_dp, fac(i))
                        end do
                    else
                        call stop_all(this_routine, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    write(iout,'("** WARNING ** Death probability > 1: Creating Antiparticles. "&
                        & //"Timestep errors possible: ")',advance='no')
                    do i = 1, inum_runs
                        write(iout,'(1X,f13.7)',advance='no') fac(i)
                    end do
                    write(iout,'()')
                endif
            endif

            if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel .le. RealCoeffExcitThresh)) &
                .or. tAllRealCoeff ) then

                ! i exact.. just get the weights, this should also still work.
                ! but have to exchange the weights to come from the other 
                ! type of particles..
                ! the number of deaths has +sign from Im -> Re
                ! dont use abs() here compared to the old code, but already 
                ! here include the sign of the walkers on the determinant
                ndie(1) = fac(1) * realwSign(2)
                ! and - from Re -> Im
                ! does this give the correct sign compared to the parent sign?
                ndie(2) = -fac(1) * realwSign(1)

            else 
                ! if not exact i have to round stochastically
                ! Re -> Im
                rat = fac(1) * RealwSign(2)

                ndie(1) = real(int(rat), dp) 
                rat = rat - ndie(1)

                r = genrand_real2_dSFMT() 
                if (abs(rat) > r) ndie(1) = ndie(1) + real(nint(sign(1.0_dp,rat)),dp)

                ! Im -> Re
                rat = -fac(1) * RealwSign(1)
                ndie(2) = real(int(rat), dp)
                rat = rat - ndie(2)
                r = genrand_real2_dSFMT()
                if (abs(rat) > r) ndie(2) = ndie(2) + real(nint(sign(1.0_dp,rat)),dp)
            end if

            ! here i only have influence from the other type of walkers, which 
            ! is already stored in the ndie() array
        else
            ! there is an imaginary energy damping factor.. 
            ! so i have mixed influences for the diagonal step
            ! n_i' = -e n_i + H_ii c_i
            ! c_i' = -e c_i - H_ii n_i

            ! temporarily use the 2 entries of fac to store both influences 
            ! not sure about the sign of this fac factor but the 2nd entry
            ! as it is implemented right now has to have the opposite sign 
            ! of the first on above
            fac(2) = tau * real_time_info%damping 

            ! here i am definetly not sure about the logging of the death
            ! magnitude.. 
            ! and also about the fac restrictions.. for now but it here anyway..
            if(any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTau) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do i = 1, inum_runs
                            fac(i) = min(2.0_dp, fac(i))
                        end do
                    else
                        call stop_all(this_routine, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    write(iout,'("** WARNING ** Death probability > 1: Creating Antiparticles. "&
                        & //"Timestep errors possible: ")',advance='no')
                    do i = 1, inum_runs
                        write(iout,'(1X,f13.7)',advance='no') fac(i)
                    end do
                    write(iout,'()')
                endif
            endif

            if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel .le. RealCoeffExcitThresh)) &
                .or. tAllRealCoeff ) then

                ! i exact.. just get the weights, this should also still work.
                ! but have to exchange the weights to come from the other 
                ! type of particles..
                ! the number of deaths has +sign from Im -> Re
                ! can i just add the other contribution here? 
                ! and also include the sign of the parent occupation here
                ! already.
                ndie(1) = fac(1) * realwSign(2) + fac(2) * realwSign(1)
                ! and - from Re -> Im
                ! does this give the correct sign compared to the parent sign?
                ndie(2) = -fac(1) * abs(realwSign(1)) + fac(2) * abs(RealwSign(2))

            else 
                ! if not exact i have to round stochastically
                ! is this ok here to just add the second contribution? todo
                ! Re -> Im
                rat = fac(1) * RealwSign(2) + fac(2) * RealwSign(1)

                ndie(1) = real(int(rat), dp) 
                rat = rat - ndie(1)

                r = genrand_real2_dSFMT() 
                if (abs(rat) > r) ndie(1) = ndie(1) + real(nint(sign(1.0_dp,rat)),dp)

                ! Im -> Re
                rat = -fac(1) * RealwSign(1) + fac(2) * RealwSign(2)

                ndie(2) = real(int(rat), dp)
                rat = rat - ndie(2)

                r = genrand_real2_dSFMT()
                if (abs(rat) > r) ndie(2) = ndie(2) + real(nint(sign(1.0_dp,rat)),dp)
            end if
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

        integer, intent(in) :: DetCurr(nel) 
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        real(dp), intent(in) :: Kii
        integer, intent(in) :: DetPosition
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign) :: ndie
        real(dp), dimension(lenof_sign) :: CopySign
        integer, intent(in) :: walkExcitLevel
        integer :: i, irdm
        real(dp) :: fac(lenof_sign), rat, r

        character(len=*), parameter :: t_r = "walker_death_realtime"

        ndie = attempt_die_realtime(DetCurr, RealwSign, ilutCurr, Kii, DetPosition, &
            iter_data, walkExcitLevel)
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
        if (real_time_info%damping > EPS) then
            ! and only if there is damping -> otherwise no 'death' 
            ! or do i have to check if the spawn from Re <-> Im caused 
            ! annihilated particles?!
!             
!         iter_data%ndied = iter_data%ndied + min(iDie, abs(RealwSign))
! #ifdef __CMPLX
!         NoDied(1) = NoDied(1) + sum(min(iDie, abs(RealwSign)))
! #else
!         NoDied = NoDied + min(iDie, abs(RealwSign))
! #endif
! 
!         ! Count any antiparticles
!         iter_data%nborn = iter_data%nborn + max(iDie - abs(RealwSign), 0.0_dp)
! #ifdef __CMPLX
!         NoBorn(1) = NoBorn(1) + sum(max(iDie - abs(RealwSign), 0.0_dp))
! #else
!         NoBorn = NoBorn + max(iDie - abs(RealwSign), 0.0_dp)
! #endif
        end if
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
                    print *, "anti-particle even! but do not cancel this in the rt-fciqmc!"
!                     NoAborted = NoAborted + abs(CopySign(i))
!                     iter_data%naborted(i) = iter_data%naborted(i) &
!                                           + abs(CopySign(i))
!                     if (test_flag(ilutCurr, flag_initiator(i))) &
!                         NoAddedInitiators = NoAddedInitiators - 1
!                     CopySign(i) = 0
                end if
            end do
        end if

        if (any(CopySign /= 0)) then
            ! For the hashed walker main list, the particles don't move.
            ! Therefore just adjust the weight.
            call encode_sign (CurrentDets(:, DetPosition), CopySign)
        else
            ! All walkers died.
            if(tFillingStochRDMonFly) then
                do irdm = 1, nrdms
                    call det_removed_fill_diag_rdm(rdms(irdm), irdm, CurrentDets(:,DetPosition), DetPosition)
                end do
                ! Set the average sign and occupation iteration to zero, so
                ! that the same contribution will not be added in in
                ! CalcHashTableStats, if this determinant is not overwritten
                ! before then
                global_determinant_data(:, DetPosition) = 0.0_dp
            endif

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
            tgt_cpt = 3 - part_type
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
            rh_used = conjg(rh)

            ! have to loop over the tgt_cpt similar to the complex implo
            ! if the Hamiltonian has real and imaginary components do it 
            ! similarily to complex implementation with H <-> J switched
            ! change this (lenof/inum) below .. todo
            do tgt_cpt = 1, (lenof_sign/inum_runs)
                component = tgt_cpt
                if (part_type == 2 .and. inum_runs == 1) component = 3 - tgt_cpt

                walkerweight = sign(1.0_dp,RealwSign(part_type)) 

                if (component == 1) then
                    MatEl = real(aimag(rh_used),dp)
                else 
                    MatEl = real(rh_used,dp)
                    if (part_type == 2) walkerweight = -walkerweight
                end if

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
        temp_det_hash = HashIndex

        ! and the freeslot.. although this one gets reinitialized to 0 
        ! every iteration or not? yeah it is.. so i only have to reset it 
        ! twice in the rt-fciqmc before the y(n) + k2 combination
        ! do that in the reload_current_dets routine!
        ! same with n_determ_states var.
        
    end subroutine save_current_dets

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

        ! and the hash
        HashIndex = temp_det_hash

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
        temp_det_list(0:tmp_siz1,1:tmp_siz2) = 0
        
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

    ! subroutine to calculate the overlap of the current y(t) = a_j(a^+_j)(t)y(0)>
    ! time evolved wavefunction to the saved <y(0)|a^+_i(a_i) 
    subroutine update_gf_overlap() 
        ! this routine only deals with globally defined variables
        integer :: idet, nI(nel), det_ind, hash_val
        real(dp) :: overlap, real_sign_1(lenof_sign), real_sign_2(lenof_sign)
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
                call extract_sign(CurrentDets(:,det_ind), real_sign_2)

                overlap = overlap + real_sign_1(1) * real_sign_2(1)
                
            end if
        end do

        ! need the timestep here... or the cycle of the current real-time loop
        gf_overlap(0) = overlap 

    end subroutine update_gf_overlap

    subroutine create_perturbed_ground()
        ! routine to create from the already read in popsfile info in 
        ! popsfile_dets the left hand <y(0)| by applying the corresponding 
        ! creation or annihilation operator
        character(*), parameter :: this_routine = "create_perturbed_ground"
        integer(int64) :: tmp_totwalkers
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

end module real_time_procs
