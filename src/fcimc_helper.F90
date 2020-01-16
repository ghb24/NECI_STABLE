#include "macros.h"

module fcimc_helper

    use constants
    use util_mod
    use systemData, only: nel, tHPHF, tNoBrillouin, G1, tUEG, &
                          tLatticeGens, nBasis, tHistSpinDist, tRef_Not_HF, &
                          t_3_body_excits, t_non_hermitian, t_ueg_3_body, t_mol_3_body
    use HPHFRandExcitMod, only: ReturnAlphaOpenDet
    use semi_stoch_procs, only: recalc_core_hamil_diag, is_core_state
    use bit_reps, only: NIfTot, test_flag, extract_flags, &
                        encode_bit_rep, NIfD, set_flag_general, NIfDBO, &
                        extract_sign, set_flag, encode_sign, &
                        flag_trial, flag_connected, flag_deterministic, &
                        extract_part_sign, encode_part_sign, decode_bit_det, &
                        get_initiator_flag, get_initiator_flag_by_run, &

                        log_spawn, increase_spawn_counter, all_runs_are_initiator, &
                        encode_spawn_hdiag, extract_spawn_hdiag, flag_static_init
    use bit_rep_data, only: flag_determ_parent

    use DetBitOps, only: FindBitExcitLevel, FindSpatialBitExcitLevel, &
                         DetBitEQ, count_open_orbs, EncodeBitDet, &
                         TestClosedShellDet
    use Determinants, only: get_helement, write_det
    use FciMCData
    use hist, only: test_add_hist_spin_dist_det, add_hist_spawn, &
                    add_hist_energies, HistMinInd
    use hist_data, only: tHistSpawn
    use hphf_integrals, only: hphf_off_diag_helement
    use LoggingData, only: OrbOccs, tPrintOrbOcc, tPrintOrbOccInit, &
                           tHistEnergies, tExplicitAllRDM, RDMExcitLevel, &
                           RDMEnergyIter, tFullHFAv, tLogComplexPops, &
                           nHistEquilSteps, tCalcFCIMCPsi, StartPrintOrbOcc, &
                           HistInitPopsIter, tHistInitPops, iterRDMOnFly, &
                           FciMCDebug, tLogEXLEVELStats, maxInitExLvlWrite, &
                           initsPerExLvl
    use CalcData, only: NEquilSteps, tFCIMC, tTruncCAS, &
                        InitiatorWalkNo, &
                        tTruncInitiator, tTruncNopen, trunc_nopen_max, &
                        tRealCoeffByExcitLevel, tGlobalInitFlag, tInitsRDM, &
                        tSemiStochastic, tTrialWavefunction, DiagSft, &
                        MaxWalkerBloom, tEN2, tEN2Started, &
                        NMCyc, iSampleRDMIters, ErrThresh, &
                        tOrthogonaliseReplicas, tPairedReplicas, t_back_spawn, &
                        t_back_spawn_flex, tau, DiagSft, &
                        tSeniorInitiators, SeniorityAge, tInitCoherentRule, &
                        tLogAverageSpawns, &
                        spawnSgnThresh, minInitSpawns, &
                        tAutoAdaptiveShift, tAAS_MatEle, tAAS_MatEle2,&
                        tAAS_MatEle3, tAAS_MatEle4, AAS_DenCut, &
                        tPrecond, &
                        tReplicaEstimates, tInitiatorSpace, tPureInitiatorSpace, tSimpleInit, &
                        allowedSpawnSign, tAS_Offset, FullShiftOffset
    use adi_data, only: tSignedRepAv
    use IntegralsData, only: tPartFreezeVirt, tPartFreezeCore, NElVirtFrozen, &
                             nPartFrozen, nVirtPartFrozen, nHolesFrozen
    use procedure_pointers, only: attempt_die, extract_bit_rep_avsign, &
         scaleFunction
    use DetCalcData, only: FCIDetIndex, ICILevel, det
    use hash, only: remove_hash_table_entry, add_hash_table_entry, hash_table_lookup
    use load_balance_calcnodes, only: DetermineDetNode, tLoadBalanceBlocks
    use load_balance, only: adjust_load_balance, RemoveHashDet, get_diagonal_matel
    use rdm_filling, only: det_removed_fill_diag_rdm
    use rdm_general, only: store_parent_with_spawned, extract_bit_rep_avsign_norm
    use Parallel_neci
    use FciMCLoggingMod, only: HistInitPopulations, WriteInitPops
    use csf_data, only: csf_orbital_mask
    use csf, only: iscsf
    use hphf_integrals, only: hphf_diag_helement
    use global_det_data, only: get_av_sgn_tot, set_av_sgn_tot, set_det_diagH, &
                               global_determinant_data, det_diagH, &
                               get_spawn_pop, get_tau_int, get_shift_int, &
                               get_neg_spawns, get_pos_spawns
    use searching, only: BinSearchParts2
    use back_spawn, only: setup_virtual_mask
    use initiator_space_procs, only: is_in_initiator_space
    use fortran_strings, only: str

    implicit none

    save

    interface CalcParentFlag
       module procedure CalcParentFlag_normal
       module procedure CalcParentFlag_det
    end interface CalcParentFlag

    interface TestInitiator
       module procedure TestInitiator_ilut
       module procedure TestInitiator_explicit
    end interface TestInitiator

contains

    function TestMCExit(Iter,RDMSamplingIter) result(ExitCriterion)
        implicit none
        integer, intent(in) :: Iter,RDMSamplingIter
        logical :: ExitCriterion

        ExitCriterion = .false.
        if((Iter.gt.NMCyc).and.(NMCyc.ne.-1)) then
            write(6,"(A)") "Total iteration number limit reached. Finishing FCIQMC loop..."
            ExitCriterion = .true.
        elseif((RDMSamplingIter.gt.iSampleRDMIters).and.(iSampleRDMIters.ne.-1)) then
            write(6,"(A)") "RDM Sampling iteration number limit reached. Finishing FCIQMC loop..."
            ExitCriterion = .true.
        endif

    end function TestMCExit

    subroutine create_particle (nJ, iLutJ, child, part_type, hdiag_spawn, err, ilutI, SignCurr, &
                                WalkerNo, RDMBiasFacCurr, WalkersToSpawn, matel, ParentPos)
        ! Create a child in the spawned particles arrays. We spawn particles
        ! into a separate array, but non-contiguously. The processor that the
        ! newly-spawned particle is going to be sent to has to be determined,
        ! and then it will be put into the appropriate element determined by
        ! ValidSpawnedList

        ! 'type' of the particle - i.e. real/imag
        integer, intent(in) :: nJ(nel), part_type
        integer(n_int), intent(in) :: iLutJ(0:niftot)
        real(dp), intent(in) :: child(lenof_sign)
        integer, intent(out) :: err
        HElement_t(dp), intent(in) :: hdiag_spawn
        integer(n_int), intent(in), optional :: ilutI(0:niftot)
        real(dp), intent(in), optional :: SignCurr(lenof_sign)
        integer, intent(in), optional :: WalkerNo
        real(dp), intent(in), optional :: RDMBiasFacCurr
        integer, intent(in), optional :: WalkersToSpawn
        real(dp), intent(in), optional :: matel
        integer, intent(in), optional :: ParentPos
        integer :: proc, j, run
        real(dp) :: r
        integer, parameter :: flags = 0
        logical :: list_full, allowed_child
        character(*), parameter :: this_routine = 'create_particle'

        logical :: parent_init
        real(dp)  :: weight_acc, weight_rej, weight_rev, weight_den, weight_den2

        ! error flag to indicate if the attempt was successful
        err = 0
        !Ensure no cross spawning between runs - run of child same as run of
        !parent
        run = part_type_to_run(part_type)
#ifdef DEBUG_
        ASSERT(sum(abs(child))-sum(abs(child(min_part_type(run):max_part_type(run)))) < 1.0e-12_dp)
#endif

        ! Determine which processor the particle should end up on in the
        ! DirectAnnihilation algorithm.
        proc = DetermineDetNode(nel,nJ,0)    ! (0 -> nNodes-1)

        if(checkValidSpawnedList(proc)) then
           err = 1
           return
        endif
        !We initially encode no flags
        call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), iLutJ, &
                            child, flags)

        ! If the parent was an initiator then set the initiator flag for the
        ! child, to allow it to survive.
        if (tTruncInitiator) then
           allowed_child = .false.
           ! optionally: allow all spawns with a given sign
           if(allowedSpawnSign.ne.0) then
              if(allowedSpawnSign * child(part_type) * SignCurr(part_type)> 0) &
                   allowed_child = .true.
           endif
            if (allowed_child .or. test_flag(ilutI, get_initiator_flag(part_type))) then
                call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), get_initiator_flag(part_type))
            endif
        end if

        if(tAutoAdaptiveShift)then
            SpawnInfo(SpawnParentIdx, ValidSpawnedList(proc)) = ParentPos
            SpawnInfo(SpawnRun, ValidSpawnedList(proc)) = run
            if(tAAS_MatEle)then
                weight_acc = abs(matel)
                weight_rej = abs(matel)
            else if(tAAS_MatEle2) then
                weight_den = abs((get_diagonal_matel(nJ, ilutJ)-Hii) - DiagSft(run))
                if(weight_den<AAS_DenCut)then
                    weight_den = AAS_DenCut
                end if
                weight_acc = abs(matel)/weight_den
                weight_rej = abs(matel)/weight_den
            else if(tAAS_MatEle3) then
                weight_den = abs((get_diagonal_matel(nJ, ilutJ)-Hii) - DiagSft(run))
                if(weight_den<AAS_DenCut)then
                    weight_den = AAS_DenCut
                end if
                weight_acc = 1.0_dp
                weight_rej = abs(matel)/weight_den
            else if(tAAS_MatEle4) then
                weight_den = abs((get_diagonal_matel(nJ, ilutJ)-Hii) - DiagSft(run))
                if(weight_den<AAS_DenCut)then
                    weight_den = AAS_DenCut
                end if
                weight_den2 = abs((get_diagonal_matel(nJ, ilutJ)-Hii))
                if(weight_den2<AAS_DenCut)then
                    weight_den2 = AAS_DenCut
                end if
                weight_rej = abs(matel)/weight_den
                weight_acc = abs(matel)/weight_den2
            else
                weight_acc = 1.0_dp
                weight_rej = 1.0_dp
            end if

            !Enocde weight, which is real, as an integer
            SpawnInfo(SpawnWeightAcc, ValidSpawnedList(proc)) = transfer(weight_acc, SpawnInfo(SpawnWeightAcc, ValidSpawnedList(proc)))
            SpawnInfo(SpawnWeightRej, ValidSpawnedList(proc)) = transfer(weight_rej, SpawnInfo(SpawnWeightRej, ValidSpawnedList(proc)))

        end if

        ! store global data - number of spawns
        ! Is using the pure initiator space option, then if this spawning
        ! occurs to within the defined initiator space (regardless of
        ! whether or not it is occupied already), then it should never be
        ! rejected by the initiator criterion, so set the initiator
        ! flag now if not done already.
        if (tPureInitiatorSpace) then
            if ( .not. test_flag(SpawnedParts(:, ValidSpawnedList(proc)), get_initiator_flag(part_type)) ) then
                if (is_in_initiator_space(SpawnedParts(:, ValidSpawnedList(proc)), nJ)) then
                    call set_flag(SpawnedParts(:,ValidSpawnedList(proc)), get_initiator_flag(part_type))
                end if
            end if
        end if
        if(tLogNumSpawns) call log_spawn(SpawnedParts(:,ValidSpawnedList(proc) ) )

        if (tFillingStochRDMonFly) then
            ! We are spawning from ilutI to
            ! SpawnedParts(:,ValidSpawnedList(proc)). We want to store the
            ! parent (D_i) with the spawned child (D_j) so that we can add in
            ! Ci.Cj to the RDM later.
            ! The parent is NIfDBO integers long, and stored in the second
            ! part of the SpawnedParts array from NIfTot+1 --> NIfTot+1+NIfDBO
            call store_parent_with_spawned (RDMBiasFacCurr, WalkerNo, &
                                            ilutI, WalkersToSpawn, ilutJ, &
                                            proc)
        end if

        if (tPreCond .or. tReplicaEstimates) then
            call encode_spawn_hdiag(SpawnedParts(:, ValidSpawnedList(proc)), hdiag_spawn)
        end if

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1

        ! Sum the number of created children to use in acceptance ratio.
        ! Note that if child is an array, it should only have one non-zero
        ! element which has changed.
        ! rmneci_setup: clarified dependence of run on part_type
        run = part_type_to_run(part_type)
        acceptances(run) = acceptances(run) + sum(abs(child(min_part_type(run):max_part_type(run))))

    end subroutine create_particle


    subroutine create_particle_with_hash_table (nI_child, ilut_child, child_sign, part_type, ilut_parent, iter_data, err)
        use hash, only: hash_table_lookup, add_hash_table_entry
        integer, intent(in) :: nI_child(nel), part_type
        integer(n_int), intent(in) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        real(dp), intent(in) :: child_sign(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(out) :: err

        integer :: proc, ind, hash_val_cd, hash_val, i, run
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_old(lenof_sign), real_sign_new(lenof_sign)
        real(dp) :: sgn_prod(lenof_sign)
        logical :: list_full, tSuccess, allowed_child
        integer :: global_position
        integer, parameter :: flags = 0
        character(*), parameter :: this_routine = 'create_particle_with_hash_table'

        err = 0
        ! Only one element of child should be non-zero
        ASSERT((sum(abs(child_sign))-maxval(abs(child_sign)))<1.0e-12_dp)

        if (tSimpleInit) then
            call stop_all(this_routine, "Cannot use a hash table to the spawned list when using the &
                                        &simple-initiator option.")
        end if

        call hash_table_lookup(nI_child, ilut_child, NIfDBO, spawn_ht, SpawnedParts, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! If the spawned child is already in the spawning array.
            ! Extract the old sign.
            call extract_sign(SpawnedParts(:,ind), real_sign_old)

            ! If the new child has an opposite sign to that of walkers already
            ! on the site, then annihilation occurs. The stats for this need
            ! accumulating.
            if (.not. tPrecond) then
                sgn_prod = real_sign_old * child_sign
                do i = 1, lenof_sign
                    if (sgn_prod(i) < 0.0_dp) then
                        iter_data%nannihil(i) = iter_data%nannihil(i) + 2*min( abs(real_sign_old(i)), abs(child_sign(i)) )
                    end if
                end do
            end if

            ! Find the total new sign.
            real_sign_new = real_sign_old + child_sign
            ! Encode the new sign.
            call encode_sign(SpawnedParts(:,ind), real_sign_new)


            ! this is not correctly considered for the real-time or complex
            ! code .. probably nobody thought about using this in the __cmplx
            ! implementation..

            ! (There is now an option (tInitCoherentRule = .false.) to turn this
            ! coherent spawning rule off, mainly for testing purposes).
            if (tTruncInitiator) then
                if (tInitCoherentRule) then
                    if (abs(real_sign_old(part_type)) > 1.e-12_dp .or. test_flag(ilut_parent, get_initiator_flag(part_type))) &
                        call set_flag(SpawnedParts(:,ind), get_initiator_flag(part_type))
                else
                    if (test_flag(ilut_parent, get_initiator_flag(part_type))) &
                        call set_flag(SpawnedParts(:,ind), get_initiator_flag(part_type))
                end if
            end if


            ! log the spawn
            global_position = ind
        else
            ! Determine which processor the particle should end up on in the
            ! DirectAnnihilation algorithm.
            proc = DetermineDetNode(nel, nI_child, 0)

            ! Check that the position described by ValidSpawnedList is acceptable.
            ! If we have filled up the memory that would be acceptable, then
            ! kill the calculation hard (i.e. stop_all) with a descriptive
            ! error message.
            if(checkValidSpawnedList(proc)) then
               err = 1
               return
            endif

            call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), ilut_child(0:NIfDBO), child_sign, flags)

            ! If the parent was an initiator then set the initiator flag for the
            ! child, to allow it to survive.

            if (tTruncInitiator) then
               if (test_flag(ilut_parent, get_initiator_flag(part_type))) &
                    call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), &
                    get_initiator_flag(part_type))
            end if

             ! where to store the global data
             global_position = ValidSpawnedList(proc)

            call add_hash_table_entry(spawn_ht, ValidSpawnedList(proc), hash_val)
            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        end if

        ! store global data
        if(tLogNumSpawns) call increase_spawn_counter(SpawnedParts(:,global_position))

        ! Sum the number of created children to use in acceptance ratio.
        ! in the rt-fciqmc i have to track the stats of the 2 RK steps
        ! seperately
        ! RT_M_Merge: Merge
        ! rmneci_setup: introduced multirun support, fixed issue in non
        ! real-time scheme
        run = part_type_to_run(part_type)
        acceptances(run) = &
             acceptances(run) + maxval(abs(child_sign))

      end subroutine create_particle_with_hash_table

    function checkValidSpawnedList(proc) result(list_full)
        ! Check that the position described by ValidSpawnedList is acceptable.
        ! If we have filled up the memory that would be acceptable, then
        ! end the calculation, i.e. throw an error
        implicit none
        logical :: list_full
        integer, intent(in) :: proc
        list_full = .false.
        if (proc == nNodes - 1) then
            if (ValidSpawnedList(proc) > MaxSpawned) list_full = .true.
        else
            if (ValidSpawnedList(proc) >= InitialSpawnedSlots(proc+1)) &
                list_full=.true.
        end if
        if (list_full) then
#ifdef DEBUG_
            write(6,*) "Attempting to spawn particle onto processor: ", proc
            write(6,*) "No memory slots available for this spawn."
            write(6,*) "Please increase MEMORYFACSPAWN"
            write(6,*) ValidSpawnedList
            write(6,*) InitialSpawnedSlots
            if(MaxSpawned / nProcessors < 0.1_dp * TotWalkers) write(6,*) &
                "Memory available for spawns is too low, number of processes might be too high for the given walker number"            
#else
            print *, "Attempting to spawn particle onto processor: ", proc
            print *, "No memory slots available for this spawn, terminating calculation"
            print *, "Please increase MEMORYFACSPAWN"
            print *, ValidSpawnedList
            print *, InitialSpawnedSlots

            ! give a note on the counter-intuitive scaling behaviour
            if(MaxSpawned / nProcessors < 0.1_dp * TotWalkers) print *, &
                "Memory available for spawns is too low, number of processes might be too high for the given walker number"    
#endif
        end if
        
    end function checkValidSpawnedList

    ! This routine sums in the energy contribution from a given walker and
    ! updates stats such as mean excit level AJWT added optional argument
    ! dProbFin which is a probability that whatever gave this contribution
    ! was generated. It defaults to 1, and weights the contribution of this
    ! det (only in the projected energy) by dividing its contribution by
    ! this number
    subroutine SumEContrib (nI, ExcitLevel, RealWSign, ilut, HDiagCurr, &
                            dProbFin, tPairedReplicas, ind)

        use CalcData, only: qmc_trial_wf
        use searching, only: get_con_amp_trial_space

        integer, intent(in) :: nI(nel), ExcitLevel
        real(dp), intent(in) :: RealwSign(lenof_sign)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: HDiagCurr, dProbFin
        logical, intent(in) :: tPairedReplicas
        integer, intent(in), optional :: ind

        integer :: i, ExcitLevel_local, ExcitLevelSpinCoup
        integer :: run
        HElement_t(dp) :: HOffDiag(inum_runs)
        character(*), parameter :: this_routine = 'SumEContrib'

#ifdef CMPLX_
        complex(dp) :: CmplxwSign
#endif

        real(dp) :: amps(size(current_trial_amps,1))
        real(dp) :: w(0:2)

        if (tReplicaReferencesDiffer) then
            call SumEContrib_different_refs(nI, realWSign, ilut, dProbFin, tPairedReplicas, ind)
            return
        endif

        HOffDiag = 0.0_dp

        ! Add in the contributions to the numerator and denominator of the trial
        ! estimator, if it is being used.
#ifdef CMPLX_
        CmplxwSign = ARR_RE_OR_CPLX(realwsign, 1)

        if (tTrialWavefunction .and. present(ind)) then
            if (test_flag(ilut,flag_trial)) then
                if(ntrial_excits == 1) then
                   trial_denom = trial_denom + conjg(current_trial_amps(1,ind))*CmplxwSign
                   ! this does somehow not support kmneci
                   if(test_flag(ilut, get_initiator_flag_by_run(1))) &
                        init_trial_denom = init_trial_denom + conjg(&
                        current_trial_amps(1,ind))*CmplxwSign
                else if(ntrial_excits == lenof_sign) then
                   call stop_all(this_routine, 'ntrial_excits has to be 1 currently for complex')
                end if

                if(qmc_trial_wf) then
                   call stop_all(this_routine, 'qmc_trial_wf currently not implemented for complex')
                end if
            else if (test_flag(ilut,flag_connected)) then
                if(ntrial_excits == 1) then
                   trial_numerator = trial_numerator + conjg(current_trial_amps(1,ind))*cmplxwsign
                   if(test_flag(ilut, get_initiator_flag_by_run(1))) &
                        init_trial_numerator = init_trial_numerator + conjg(&
                        current_trial_amps(1,ind))*CmplxwSign
                else if(ntrial_excits == lenof_sign) then
                   call stop_all(this_routine, 'ntrial_excits has to be 1 currently for complex')
                end if
            end if
        end if
#else
        if (tTrialWavefunction .and. present(ind)) then
            if (test_flag(ilut, flag_trial)) then
                if (ntrial_excits == 1) then
                    trial_denom = trial_denom + current_trial_amps(1,ind)*RealwSign
                    trial_denom_inst = trial_denom_inst + current_trial_amps(1,ind)*RealwSign
                    do run = 1, inum_runs
                       if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                            init_trial_denom(run) = init_trial_denom(run) + &
                            current_trial_amps(1,ind) * RealwSign(run)
                    end do
                else if (ntrial_excits == lenof_sign) then
                    trial_denom = trial_denom + current_trial_amps(:,ind)*RealwSign
                    trial_denom_inst = trial_denom_inst + current_trial_amps(:,ind)*RealwSign
                    do run = 1, inum_runs
                       if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                            init_trial_denom(run) = init_trial_denom(run) + &
                            current_trial_amps(run,ind) * RealwSign(run)
                    end do
                end if

                if (qmc_trial_wf) then
                    call get_con_amp_trial_space(ilut, amps)

                    if (ntrial_excits == 1) then
                        trial_numerator = trial_numerator + amps(1)*RealwSign
                        do run = 1, inum_runs
                           if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                                init_trial_numerator(run) = init_trial_numerator(run) + &
                                amps(1) * RealwSign(run)
                        end do
                    else if (ntrial_excits == lenof_sign) then
                        trial_numerator = trial_numerator + amps*RealwSign
                        do run = 1, inum_runs
                           if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                                init_trial_numerator(run) = init_trial_numerator(run) + &
                                amps(run) * RealwSign(run)
                        end do
                    end if
                end if

            else if (test_flag(ilut, flag_connected)) then
                ! Note, only attempt to add in a contribution from the
                ! connected space if we're not also in the trial space.
                if (ntrial_excits == 1) then
                    trial_numerator = trial_numerator + current_trial_amps(1,ind)*RealwSign
                    trial_num_inst = trial_num_inst + current_trial_amps(1,ind)*RealwSign
                    do run = 1, inum_runs
                       if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                            ! this is the real case, so inum_runs == lenof_sign
                            init_trial_numerator(run) = init_trial_numerator(run) + &
                            current_trial_amps(1,ind) * RealwSign(run)
                    end do
                else if (ntrial_excits == lenof_sign) then
                    trial_numerator = trial_numerator + current_trial_amps(:,ind)*RealwSign
                    trial_num_inst = trial_num_inst + current_trial_amps(:,ind)*RealwSign
                    do run = 1, inum_runs
                       if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                            init_trial_numerator(run) = init_trial_numerator(run) + &
                            current_trial_amps(run,ind)*RealwSign(run)
                    end do
                end if
            end if
        end if
#endif

        ! ExcitLevel indicates the excitation level between the det and
        ! *one* of the determinants in an HPHF/MomInv function. If needed,
        ! calculate the connection between it and the other one. If either
        ! is connected, then it has to be counted. Since the excitation
        ! level is the same to either det, we don't need to consider the
        ! spin-coupled det of both reference and current HPHFs.
        !
        ! For determinants, set ExcitLevel_local == ExcitLevel.
        ExcitLevel_local = ExcitLevel
        if (tSpinCoupProjE(1) .and. (ExcitLevel /= 0)) then
            ExcitLevelSpinCoup = FindBitExcitLevel (iLutRefFlip(:, 1), &
                                                    ilut, 2)
            if (ExcitLevelSpinCoup <= 2 .or. ExcitLevel <= 2) &
                ExcitLevel_local = 2
        endif

        ! this is the first important change to make the triples run!!
        ! consider the matrix elements of triples!

        ! Perform normal projection onto reference determinant
        if (ExcitLevel_local == 0) then


            ! for the real-time i have to distinguish between the first and
            ! second RK step, if i want to keep track of the statistics
            ! seperately: in the first loop i analyze the the wavefunction
            ! from on step behind.. so store it in the "normal" noathf var
            if (iter > NEquilSteps) &
                SumNoatHF(1:lenof_sign) = SumNoatHF(1:lenof_sign) + RealwSign
            NoatHF(1:lenof_sign) = NoatHF(1:lenof_sign) + RealwSign
            ! Number at HF * sign over course of update cycle
            HFCyc(1:lenof_sign) = HFCyc(1:lenof_sign) + RealwSign
            HFOut(1:lenof_sign) = HFOut(1:lenof_sign) + RealwSign

        elseif (ExcitLevel_local == 2 .or. &
                (ExcitLevel_local == 1 .and. tNoBrillouin)) then

            ! For the real-space Hubbard model, determinants are only
            ! connected to excitations one level away, and Brillouins
            ! theorem cannot hold.
            !
            ! For Rotated orbitals, Brillouins theorem also cannot hold,
            ! and energy contributions from walkers on singly excited
            ! determinants must also be included in the energy values
            ! along with the doubles
           ! RT_M_Merge: Adjusted to kmneci
           ! rmneci_setup: Added multirun functionality for real-time
            if (ExcitLevel_local == 2) then
            do run = 1, inum_runs
#ifdef CMPLX_
                NoatDoubs(run) = NoatDoubs(run) + sum(abs(RealwSign(min_part_type(run):max_part_type(run))))
#else
                NoatDoubs(run) = NoatDoubs(run) + abs(RealwSign(run))
#endif
            enddo
            end if
            ! Obtain off-diagonal element
            if (tHPHF) then
                HOffDiag(1:inum_runs) = hphf_off_diag_helement (ProjEDet(:,1), nI, &
                                                                iLutRef(:,1), ilut)
            else
                HOffDiag(1:inum_runs) = get_helement (ProjEDet(:,1), nI, &
                                                      ExcitLevel, ilutRef(:,1), ilut)
            endif

        else if (ExcitLevel_local == 3 .and. (t_3_body_excits .or. t_ueg_3_body .or. t_mol_3_body)) then
            ! the new 3-body terms in the transcorrelated momentum space hubbard
            ! hphf not yet implemented!
            ASSERT(.not. tHPHF)
            HOffDiag(1:inum_runs) = get_helement( ProjEDet(:,1), nI, ilutRef(:,1), ilut)

        endif ! ExcitLevel_local == 1, 2, 3

        ! L_{0,1,2} norms of walker weights by excitation level.
        if (tLogEXLEVELStats) then
            do run = 1, inum_runs
#ifdef CMPLX_
                w(0) = real(1 + max_part_type(run) - min_part_type(run), dp)
                w(1) = sum(abs(RealwSign(min_part_type(run):&
                                         max_part_type(run))))
                w(2) = sum(RealwSign(min_part_type(run):&
                                     max_part_type(run))**2)
#else
                w(0) = 1_dp
                w(1) = abs(RealwSign(run))
                w(2) = RealwSign(run)**2
#endif
                EXLEVEL_WNorm(0:2,ExcitLevel,run) = &
                      EXLEVEL_WNorm(0:2,ExcitLevel,run) + w(0:2)
            enddo ! run
        endif ! tLogEXLEVELStats

        ! if in the real-time fciqmc: when we are in the 2nd loop
        ! return here since, the energy got already calculated in the
        ! first RK step, and doing it on the intermediate step would
        ! be meaningless

        ! Sum in energy contribution
        do run=1, inum_runs
           
           if (iter > NEquilSteps) &
                SumENum(run) = SumENum(run) + enum_contrib()

           ENumCyc(run) = ENumCyc(run) + enum_contrib()
           ENumOut(run) = ENumOut(run) + enum_contrib()
           ENumCycAbs(run) = ENumCycAbs(run) + abs(enum_contrib() )
           if(test_flag(ilut, get_initiator_flag_by_run(run))) &
                InitsENumCyc(run) = InitsENumCyc(run) + enum_contrib()
        end do

        ! -----------------------------------
        ! HISTOGRAMMING
        ! -----------------------------------

        if ((tHistSpawn .or. (tCalcFCIMCPsi .and. tFCIMC)) .and. &
            (iter >= NHistEquilSteps)) then
            ! Histogram particles by determinant
            call add_hist_spawn (ilut, RealwSign, ExcitLevel_local, dProbFin)
        elseif (tHistEnergies) then
            ! Histogram particles by energy
            call add_hist_energies (ilut, RealwSign, HDiagCurr)
        endif

        ! Are we doing a spin-projection histogram?
        if (tHistSpinDist) then
            if (tRealCoeffByExcitLevel) &
                call stop_all(this_routine, 'Not set up to use real coeffs &
                                            &with tHistSpindist')
            call test_add_hist_spin_dist_det (ilut, RealwSign)
        end if

        ! Maintain a list of the degree of occupation of each orbital
        if (tPrintOrbOcc .and. (iter >= StartPrintOrbOcc)) then
            if (iter == StartPrintOrbOcc .and. &
                 DetBitEq(ilut, ilutHF, NIfDBO)) then
                write(6,*) 'Beginning to fill the HF orbital occupation list &
                           &during iteration', iter
                if (tPrintOrbOccInit) &
                    write(6,*) 'Only doing so for initiator determinants'
            end if
            if ((tPrintOrbOccInit .and. test_flag(ilut,get_initiator_flag(1)))&
                .or. .not. tPrintOrbOccInit) then
                forall (i = 1:nel) OrbOccs(iand(nI(i), csf_orbital_mask)) &
                        = OrbOccs(iand(nI(i), csf_orbital_mask)) &
                                   + (RealwSign(1) * RealwSign(1))
            endif
         endif

       contains

         function enum_contrib() result(dE)
           implicit none
           HElement_t(dp) :: dE

           dE = (HOffDiag(run) * ARR_RE_OR_CPLX(RealwSign,run)) / dProbFin
         end function enum_contrib
        
    end subroutine SumEContrib


    subroutine SumEContrib_different_refs(nI, sgn, ilut, dProbFin, tPairedReplicas, ind)

        use CalcData, only: qmc_trial_wf
        use searching, only: get_con_amp_trial_space

        ! This is a modified version of SumEContrib for use where the
        ! projected energies need to be calculated relative to differing
        ! references.
        !
        ! Some of the arguments to SumeEContrib have been dropped, as they
        ! only refer to the first of the replicas, and thus cannot be relied on

        ! n.b. We don't want to just modify SumEContrib for this, as performing
        !      the calculations _inside_ a sum over runs would radically slow
        !      down any simulations that were not using differing references

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NifTot)
        real(dp), intent(in) :: sgn(lenof_sign), dProbFin
        logical, intent(in) :: tPairedReplicas
        integer, intent(in), optional :: ind

        integer :: run, exlevel
        HElement_t(dp) :: sgn_run
        HElement_t(dp) :: hoffdiag
        character(*), parameter :: this_routine = 'SumEContrib_different_refs'

        HElement_t(dp) :: amps(size(current_trial_amps,1))

        if (tHistSpawn .or. &
            (tCalcFCIMCPsi .and. tFCIMC) .or. tHistEnergies .or. &
            tHistSpinDist .or. tPrintOrbOcc) &
            call stop_all(this_routine, "Not yet supported: Turn off HISTSPAWN,&
                   & PRINTFCIMCPSI, PRINTORBOCCS, HISTPARTENERGIES, HIST-SPIN-DIST")

        ! Add in the contributions to the numerator and denominator of the trial
        ! estimator, if it is being used.
        if (tTrialWavefunction .and. present(ind)) then
            if (test_flag(ilut, flag_trial)) then
                if (ntrial_excits == 1) then
                    trial_denom = trial_denom + current_trial_amps(1,ind)*sgn
                else
                    if (tPairedReplicas) then
#if defined(PROG_NUMRUNS_) || defined(DOUBLERUN_)
                        do run = 2, inum_runs, 2
#ifdef CMPLX_
                            trial_denom(run-1) = trial_denom(run-1) + current_trial_amps(run/2,ind)* &
                                cmplx(sgn(min_part_type(run-1)),sgn(max_part_type(run-1)),dp)
                            trial_denom(run) = trial_denom(run) + current_trial_amps(run/2,ind)* &
                                cmplx(sgn(min_part_type(run)),sgn(max_part_type(run)),dp)
#else
                            trial_denom(run-1:run) = trial_denom(run-1:run) + current_trial_amps(run/2,ind)*sgn(run-1:run)
#endif
                        end do
#else
                        call stop_all(this_routine, "INVALID")
#endif
                    else
#ifdef CMPLX_
                        do run=1,inum_runs
                            trial_denom(run) = trial_denom(run) + current_trial_amps(run,ind)* &
                                cmplx(sgn(min_part_type(run)),sgn(max_part_type(run)),dp)
                        enddo
#else
                        trial_denom = trial_denom + current_trial_amps(:,ind)*sgn(:)
#endif
                    end if
                end if

                if (qmc_trial_wf) then
                    call get_con_amp_trial_space(ilut, amps)

                    if (ntrial_excits == 1) then
                        trial_numerator = trial_numerator + amps(1)*sgn
                    else
                        if (tPairedReplicas) then
#if defined(PROG_NUMRUNS_) || defined(DOUBLERUN_)
                            do run = 2, inum_runs, 2
                                trial_numerator(run-1:run) = trial_numerator(run-1:run) + amps(run/2)*sgn(run-1:run)
                            end do
#else
                            call stop_all(this_routine, "INVALID")
#endif
                        else
                            trial_numerator = trial_numerator + amps*sgn
                        end if
                    end if
                end if

            else if (test_flag(ilut, flag_connected)) then
                ! Note, only attempt to add in a contribution from the
                ! connected space if we're not also in the trial space.
                if (ntrial_excits == 1) then
                    trial_numerator = trial_numerator + current_trial_amps(1,ind)*sgn
                else
                    if (tPairedReplicas) then
#if defined(PROG_NUMRUNS_) || defined(DOUBLERUN_)
                        do run = 2, inum_runs, 2
#ifdef CMPLX_
                            trial_numerator(run-1) = trial_numerator(run-1) + current_trial_amps(run/2,ind)* &
                                cmplx(sgn(min_part_type(run-1)),sgn(max_part_type(run-1)),dp)
                            trial_numerator(run) = trial_numerator(run) + current_trial_amps(run/2,ind)* &
                                cmplx(sgn(min_part_type(run)),sgn(max_part_type(run)),dp)
#else
                            trial_numerator(run-1:run) = trial_numerator(run-1:run) + current_trial_amps(run/2,ind)*sgn(run-1:run)
#endif
                        end do
#else
                        call stop_all(this_routine, "INVALID")
#endif
                    else
#ifdef CMPLX_
                        do run=1,inum_runs
                            trial_numerator(run) = trial_numerator(run) + current_trial_amps(run,ind)* &
                                cmplx(sgn(min_part_type(run)),sgn(max_part_type(run)),dp)
                        enddo
#else
                        trial_numerator = trial_numerator + current_trial_amps(:,ind)*sgn
#endif
                    end if
                end if
            end if
        end if


        ! This is the normal projected energy calculation, but split over
        ! multiple runs, rather than done in one go.
        do run = 1, inum_runs

            ! We need to use the excitation level relevant for this run
            exlevel = FindBitExcitLevel(ilut, ilutRef(:, run), t_hphf_ic = .true.)
            if (tSpinCoupProjE(run) .and. exlevel /= 0) then
                if (exlevel <= 2) then
                    exlevel = 2
                else if (FindBitExcitLevel(ilut, ilutRefFlip(:,run)) <= 2) then
                    exlevel = 2
                end if
            end if
#ifdef CMPLX_
            sgn_run = cmplx(sgn(min_part_type(run)),sgn(max_part_type(run)),dp)
#else
            sgn_run = sgn(run)
#endif

            hoffdiag = 0.0_dp
            if (exlevel == 0) then

               if (iter > nEquilSteps) then
                  call add_sign_on_run(SumNoatHF)
                  call add_sign_on_run(NoatHF)
                  call add_sign_on_run(HFCyc)
                  call add_sign_on_run(HFOut)
                endif

            else if (exlevel == 2 .or. (exlevel == 1 .and. tNoBrillouin)) then

                ! n.b. Brillouins theorem cannot hold for real-space Hubbard
                ! model or for rotated orbitals.

                if (exlevel == 2) then
                    NoatDoubs(run) = NoatDoubs(run) + mag_of_run(sgn, run)
                endif

                ! Obtain the off-diagonal elements
                if (tHPHF) then
                    hoffdiag = hphf_off_diag_helement(ProjEDet(:,run), nI, &
                                                      iLutRef(:,run), ilut)
                else
                    hoffdiag = get_helement (ProjEDet(:,run), nI, exlevel, &
                                             ilutRef(:,run), ilut)
                endif

            else if (exlevel == 3 .and. (t_3_body_excits .or. t_ueg_3_body .or. t_mol_3_body) ) then
                ASSERT(.not. tHPHF)
                hoffdiag = get_helement(ProjEDet(:,run), nI, exlevel, &
                    iLutRef(:,run), ilut)

            end if


            ! Sum in energy contributions
            if (iter > nEquilSteps) &
                SumENum(run) = SumENum(run) + enum_contrib()
            ENumCyc(run) = ENumCyc(run) + enum_contrib()
            ENumOut(run) = ENumOut(run) + enum_contrib()
            ENumCycAbs(run) = ENumCycAbs(run) + abs(enum_contrib() )

        end do

      contains

        subroutine add_sign_on_run(var)
          implicit none
          real(dp), intent(inout) :: var(lenof_sign)
          
#ifdef CMPLX_
          var(min_part_type(run)) = var(min_part_type(run)) + real(sgn_run)
          var(max_part_type(run)) = var(max_part_type(run)) + aimag(sgn_run)
#else
          var(run) = var(run) + sgn_run
#endif
        end subroutine add_sign_on_run

        function enum_contrib() result(dE)
          implicit none
          HElement_t(dp) :: dE

          dE = (hoffdiag * sgn_run) / dProbFin
        end function enum_contrib
        
    end subroutine SumEContrib_different_refs

    subroutine CalcParentFlag_normal(j, parent_flags)
      implicit none
      integer, intent(in) :: j
      integer, intent(out) :: parent_flags
      integer :: nI(nel), exLvl

      ! If we do not supply nI externally, get it now.
      ! This routine mainly exists for compatibility
      call decode_bit_det(nI, CurrentDets(:,j))
      exLvl = FindBitExcitLevel(ilutRef(:,1), CurrentDets(:,j), t_hphf_ic = .true.)
      call CalcParentFlag_det(j, nI, exLvl, parent_flags)
    end subroutine CalcParentFlag_normal

    subroutine CalcParentFlag_det(j, nI,  exLvl, parent_flags)

        ! In the CurrentDets array, the flag at NIfTot refers to whether that
        ! determinant *itself* is an initiator or not. We need to decide if
        ! this willchange due to the determinant acquiring a certain
        ! population, or its population dropping below the threshold.
        ! The CurrentDets(:,j) is the determinant we are currently spawning
        ! from, so this determines the ParentInitiator flag which is passed to
        ! the SpawnedDets array and refers to whether or not the walkers
        ! *parent* is an initiator or not.

        integer, intent(in) :: j, nI(nel), exLvl
        integer, intent(out) :: parent_flags
        real(dp) :: CurrentSign(lenof_sign)
        integer :: run, nopen
        logical :: tDetinCAS, parent_init
        real(dp) :: init_tm, expected_lifetime, hdiag
        character(*), parameter :: this_routine = 'CalcParentFlag'

        call extract_sign (CurrentDets(:,j), CurrentSign)

        if (tTruncInitiator) then

            ! Now loop over the particle types, and update the flags
            do run = 1, inum_runs

                ! By default, the parent_flags are the flags of the parent.
                parent_init = test_flag (CurrentDets(:,j), get_initiator_flag_by_run(run))

                ! Should this particle be considered to be an initiator
                ! for spawning purposes.
                if (tPureInitiatorSpace) then
                    parent_init = TestInitiator_pure_space(CurrentDets(:,j), nI, j, parent_init, run)
                else
                    parent_init = TestInitiator_explicit(CurrentDets(:,j), nI, j, parent_init, &
                                                CurrentSign, exLvl, run)
                end if

                ! log the initiator
                if(parent_init) then
                   if(exLvl <= maxInitExLvlWrite .and. exLvl >0) &
                        initsPerExLvl(exLvl) = initsPerExLvl(exLvl) + 1
                endif

                ! Update counters as required.
                if (parent_init) then
                    NoInitDets(run) = NoInitDets(run) + 1_int64
                    NoInitWalk(run) = NoInitWalk(run) + mag_of_run(CurrentSign, run)
                else
                    NoNonInitDets(run) = NoNonInitDets(run) + 1_int64
                    NoNonInitWalk(run) = NoNonInitWalk(run) + mag_of_run(CurrentSign, run)
                endif

                ! Update the parent flag as required.
                call set_flag (CurrentDets(:,j), get_initiator_flag_by_run(run), parent_init)

            enddo

        endif

        ! Store this flag for use in the spawning routines...
        parent_flags = extract_flags (CurrentDets(:,j))

        ! We don't want the deterministic flag to be set in parent_flags, as
        ! that would set the same flag in the child in create_particle, which
        ! we don't want in general.
        parent_flags = ibclr(parent_flags, flag_deterministic)

        ! We don't want the child to have trial or connected flags.
        parent_flags = ibclr(parent_flags, flag_trial)
        parent_flags = ibclr(parent_flags, flag_connected)

        if ((tHistInitPops .and. mod(iter, histInitPopsIter) == 0) &
            .or. tPrintHighPop) then
             call HistInitPopulations (CurrentSign(1), j)
        endif

      end subroutine CalcParentFlag_det

      function TestInitiator_ilut(ilut, site_idx, is_init, run) result(initiator)

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(in) :: run, site_idx
        logical, intent(in) :: is_init
        integer :: nI(nel), exLvl
        real(dp) :: sgn(lenof_sign)
        logical :: initiator

        exLvl = FindBitExcitLevel(ilut, ilutRef(:,run),t_hphf_ic = .true.)
        call decode_bit_det(nI,ilut)
        call extract_sign(ilut, sgn)

        if (tPureInitiatorSpace) then
            initiator = TestInitiator_pure_space(ilut, nI, site_idx, is_init, run)
        else
            initiator = TestInitiator_explicit(ilut, nI, site_idx, is_init, sgn, exLvl, run)
        end if

      end function TestInitiator_ilut


      function TestInitiator_explicit(ilut, nI, det_idx,is_init, sgn, exLvl, run) result(initiator)
        use adi_initiators, only: check_static_init
        implicit none
        ! For a given particle (with its given particle type), should it
        ! be considered as an initiator for the purposes of spawning.
        !
        ! Inputs: The ilut, the determinant, the particles sign, and if the particle
        !         is currently considered to be an initiator.

        ! N.B. This intentionally DOES NOT directly reference part_type.
        !      This means we can call it for individual, or aggregate,
        !      particles.

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        logical, intent(in) :: is_init
        integer, intent(in) :: run, nI(nel), exLvl, det_idx
        real(dp), intent(in) :: sgn(lenof_sign)

        logical :: initiator, staticInit, popInit, spawnInit
        integer :: i

        logical :: Senior
        real(dp) :: DetAge, HalfLife, AvgShift, diagH

        ! initiator flag according to population/spawn coherence
        popInit = initiator_criterium(sgn, det_diagH(det_idx), run) .or. &
             spawn_criterium(det_idx)

        ! initiator flag according to SI or a static initiator space
        staticInit = check_static_init(ilut, nI, sgn, exLvl, run)

        if (tInitiatorSpace) then
            staticInit = test_flag(ilut, flag_static_init(run)) .or. test_flag(ilut, flag_deterministic)
            if (.not. staticInit) then
                if (is_in_initiator_space(ilut, nI)) then
                    staticInit = .true.
                    call set_flag(CurrentDets(:, det_idx), flag_static_init(run))
                end if
            end if
        end if

        ! By default the particles status will stay the same
        initiator = is_init

        ! SI-caused initiators also have the initiator flag
        if(staticInit) then
           initiator = .true.
           ! thats it, we never remove static initiators
           return
        endif

        Senior = .false.
        if (tSeniorInitiators .and. .not. is_run_unnocc(sgn, run) ) then
            DetAge = get_tau_int(det_idx, run)
            diagH = det_diagH(det_idx)
            AvgShift = get_shift_int(det_idx, run)/DetAge
            HalfLife = log(2.0_dp) / (diagH  - AvgShift)
            !Usually the shift is negative, so the HalfLife is always positive.
            !In some cases, however, the shift is set to positive (to increase the birth at HF).
            !This will to a negative HalfLife for some determinants.
            if (HalfLife>0.0) then
                Senior = DetAge>HalfLife*SeniorityAge
            end if
        end if
        if(Senior) initiator = .true.

        if (.not. initiator) then

           ! Determinant wasn't previously initiator
           ! - want to test if it has now got a large enough
           !   population to become an initiator.
           if (popInit) then
              initiator = .true.
              NoAddedInitiators = NoAddedInitiators + 1_int64
           endif

        else

           ! The determinants become
           ! non-initiators again if their population falls below
           ! n_add (this is on by default).

           ! All of the references stay initiators
           if(DetBitEQ(ilut, ilutRef(:,run),NIfDBO)) staticInit = .true.
           ! If det. is the HF det, or it
           ! is in the deterministic space, then it must remain an initiator.
           if ( .not. (staticInit) &
                .and. .not. test_flag(ilut, flag_deterministic) &
                .and. .not. Senior &
                .and. (.not. popInit )) then
              ! Population has fallen too low. Initiator status
              ! removed.
              initiator = .false.
              NoAddedInitiators = NoAddedInitiators - 1_int64
           endif

        end if

      end function TestInitiator_explicit

      function TestInitiator_pure_space(ilut, nI, site_idx, initiator_before, run) result(initiator)

          integer(n_int), intent(inout) :: ilut(0:NIfTot)
          integer, intent(in) :: nI(nel), site_idx, run
          logical, intent(in) :: initiator_before

          logical :: initiator

          ! Has this already been marked as a determinant in the static space?
          initiator = test_flag(ilut, flag_static_init(run)) .or. test_flag(ilut, flag_deterministic)

          ! If not, then it may be new, so check.
          ! Deterministic states are always in CurrentDets, so don't need to
          ! check if it's a new state in the deterministic space.
          if (.not. initiator) then
              if (is_in_initiator_space(ilut, nI)) then
                  initiator = .true.
                  call set_flag(CurrentDets(:, site_idx), flag_static_init(run))
              end if
          end if

          if (initiator .and. (.not. initiator_before)) then
              NoAddedInitiators = NoAddedInitiators + 1_int64
          else if ((.not. initiator) .and. initiator_before) then
              NoAddedInitiators = NoAddedInitiators - 1_int64
          end if

      end function TestInitiator_pure_space

      function initiator_criterium(sign,hdiag,run) result(init_flag)
        implicit none
        real(dp), intent(in) :: sign(lenof_sign), hdiag
        integer, intent(in) :: run
        ! variance of sign and either a single value or an aggregate
        real(dp) :: sigma, tot_sgn
        integer :: crun, nOcc
        real(dp) :: scaledInitiatorWalkNo
        logical :: init_flag

        if(tEScaleWalkers) then
           scaledInitiatorWalkNo = InitiatorWalkNo * scaleFunction(hdiag)
        else
           scaledInitiatorWalkNo = InitiatorWalkNo
        endif
        ! option to use the average population instead of the local one
        ! for purpose of initiator threshold
        if(tGlobalInitFlag) then
           ! we can use a signed or unsigned sum
           if(tSignedRepAv) then
              tot_sgn = real(abs(sum(sign)),dp)/inum_runs
           else
              tot_sgn = av_pop(sign)
           endif
        else
           tot_sgn = mag_of_run(sign,run)
        endif
        ! make it an initiator
        init_flag = (tot_sgn > scaledInitiatorWalkNo)

        ! option to use the average population instead of the local one
        ! for purpose of initiator threshold
        if(tGlobalInitFlag) then
           ! we can use a signed or unsigned sum
           if(tSignedRepAv) then
              tot_sgn = real(abs(sum(sign)),dp)/inum_runs
           else
              tot_sgn = av_pop(sign)
           endif
        else
           tot_sgn = mag_of_run(sign,run)
        endif
        ! make it an initiator
        init_flag = (tot_sgn > scaledInitiatorWalkNo)

      end function initiator_criterium

      function spawn_criterium(idx) result(spawnInit)
        implicit none
        ! makes something an initiator if the sign of spawns is sufficiently unique
        integer, intent(in) :: idx
        logical :: spawnInit

        real(dp) :: negSpawn(lenof_sign), posSpawn(lenof_sign)

        if(tLogAverageSpawns) then
           negSpawn = get_neg_spawns(idx)
           posSpawn = get_pos_spawns(idx)
           if(any((negSpawn + posSpawn) .ge. minInitSpawns)) then
              if(all(min(negSpawn,posSpawn) > eps)) then
                 spawnInit = all(max(negSpawn,posSpawn)/min(negSpawn,posSpawn) > spawnSgnThresh)
              else
                 spawnInit = .true.
              endif
           else
              spawnInit = .false.
           endif
        else
           spawnInit = .false.
        endif

      end function spawn_criterium

      subroutine rezero_iter_stats_each_iter(iter_data, rdm_defs)

        use global_det_data, only: len_av_sgn_tot
        use rdm_data, only: rdm_definitions_t

        type(fcimc_iter_data), intent(inout) :: iter_data
        type(rdm_definitions_t), intent(in) :: rdm_defs

        real(dp) :: prev_AvNoatHF(len_av_sgn_tot), AllInstNoatHF(lenof_sign)
        integer :: irdm, av_ind_1, av_ind_2, part_type

        NoInitDets = 0_int64
        NoNonInitDets = 0_int64
        NoInitWalk = 0.0_dp
        NoNonInitWalk = 0.0_dp
        InitRemoved = 0_int64

        ! replica-initiator info

        NoAborted = 0.0_dp
        NoRemoved = 0.0_dp
        NoatHF = 0.0_dp
        NoatDoubs = 0.0_dp
        if (tLogEXLEVELStats) EXLEVEL_WNorm = 0.0_dp

        ! for the real-time fciqmc also rezero the info on the intermediate
        ! RK step

        iter_data%nborn = 0.0_dp
        iter_data%ndied = 0.0_dp
        iter_data%nannihil = 0.0_dp
        iter_data%naborted = 0.0_dp
        iter_data%nremoved = 0.0_dp


        call InitHistMin()

        if (tFillingStochRDMonFly) then
            associate(ind => rdm_defs%sim_labels)

            call MPISumAll(InstNoatHF, AllInstNoAtHF)
            InstNoAtHF = AllInstNoAtHF

            if (tFullHFAv) then
                Prev_AvNoatHF = AvNoatHF

                do irdm = 1, rdm_defs%nrdms
                    if (abs(IterRDM_HF(irdm)) > 1.0e-12_dp) then
                        AvNoatHF(irdm) = ( (real((Iter+PreviousCycles - IterRDM_HF(irdm)),dp) * Prev_AvNoatHF(irdm)) &
                                                + InstNoatHF(irdm) ) / real((Iter+PreviousCycles - IterRDM_HF(irdm)) + 1, dp)
                    end if
                end do

            else
                if (((Iter+PreviousCycles-IterRDMStart) > 0) .and. &
                    & (mod(((Iter-1)+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0)) then
                    ! The previous iteration was one where we added in diagonal
                    ! elements To keep things unbiased, we need to set up a new
                    ! averaging block now.
                    do irdm = 1, rdm_defs%nrdms
                        av_ind_1 = irdm*2-1
                        av_ind_2 = irdm*2

                        AvNoAtHF(av_ind_1) = InstNoAtHF(ind(1,irdm))
                        IterRDM_HF(av_ind_1) = Iter + PreviousCycles
                        AvNoAtHF(av_ind_2) = InstNoAtHF(ind(2,irdm))
                        IterRDM_HF(av_ind_2) = Iter + PreviousCycles
                    end do
                else
                    do irdm = 1, rdm_defs%nrdms
                        ! The indicies of the first and second replicas in this
                        ! particular pair, in the *average* sign arrays.
                        av_ind_1 = irdm*2-1
                        av_ind_2 = irdm*2

                        if ((abs(InstNoAtHF(ind(1,irdm))) < 1.0e-12_dp .and. abs(IterRDM_HF(av_ind_1)) > 1.0e-12_dp) .or.&
                            (abs(InstNoAtHF(ind(2,irdm))) < 1.0e-12_dp .and. abs(IterRDM_HF(av_ind_2)) > 1.0e-12_dp)) then
                            ! At least one of the populations has just become
                            ! zero. Start a new averaging block.
                            IterRDM_HF(av_ind_1) = Iter + PreviousCycles
                            IterRDM_HF(av_ind_2) = Iter + PreviousCycles
                            AvNoatHF(av_ind_1) = InstNoAtHF(ind(1,irdm))
                            AvNoatHF(av_ind_2) = InstNoAtHF(ind(2,irdm))
                            if (abs(InstNoAtHF(ind(1,irdm))) < 1.0e-12_dp) IterRDM_HF(av_ind_1) = 0
                            if (abs(InstNoAtHF(ind(2,irdm))) < 1.0e-12_dp) IterRDM_HF(av_ind_2) = 0

                        else if ((abs(InstNoAtHF(ind(1,irdm))) > 1.0e-12_dp .and. abs(IterRDM_HF(av_ind_1)) < 1.0e-12_dp) .or.&
                                 (abs(InstNoAtHF(ind(2,irdm))) > 1.0e-12_dp .and. abs(IterRDM_HF(av_ind_2)) < 1.0e-12_dp)) then
                            ! At least one of the populations has just
                            ! become occupied. Start a new block here.
                            IterRDM_HF(av_ind_1) = Iter + PreviousCycles
                            IterRDM_HF(av_ind_2) = Iter + PreviousCycles
                            AvNoAtHF(av_ind_1) = InstNoAtHF(ind(1,irdm))
                            AvNoAtHF(av_ind_2) = InstNoAtHF(ind(2,irdm))
                            if (abs(InstNoAtHF(ind(1,irdm))) < 1.0e-12_dp) IterRDM_HF(av_ind_1) = 0
                            if (abs(InstNoAtHF(ind(2,irdm))) < 1.0e-12_dp) IterRDM_HF(av_ind_2) = 0
                        else
                            Prev_AvNoatHF = AvNoatHF

                            if (abs(IterRDM_HF(av_ind_1)) > 1.0e-12_dp) then
                                 AvNoatHF(av_ind_1) = ((real((Iter+PreviousCycles - IterRDM_HF(av_ind_1)), dp) &
                                                            * Prev_AvNoatHF(av_ind_1)) + InstNoatHF(ind(1,irdm)) ) &
                                                            / real((Iter+PreviousCycles - IterRDM_HF(av_ind_1)) + 1, dp)
                            end if
                            if (abs(IterRDM_HF(av_ind_2)) > 1.0e-12_dp) then
                                 AvNoatHF(av_ind_2) = ((real((Iter+PreviousCycles - IterRDM_HF(av_ind_2)), dp) &
                                                            * Prev_AvNoatHF(av_ind_2)) + InstNoatHF(ind(2,irdm)) ) &
                                                            / real((Iter+PreviousCycles - IterRDM_HF(av_ind_2)) + 1, dp)
                            end if

                        end if
                    end do

                end if
            end if

            end associate
        end if

        HFInd = 0

        min_trial_ind = 1
        min_conn_ind = 1

        trial_num_inst = 0.0_dp
        trial_denom_inst = 0.0_dp

    end subroutine rezero_iter_stats_each_iter

    subroutine InitHistMin ()

        ! Initialize the Histogramming searching arrays if necessary
        if (tHistSpawn .or. tCalcFCIMCPsi) then
            if (iter == nHistEquilSteps) then
                 root_print 'The iteration is equal to HISTEQUILSTEPS. &
                            &Beginning to histogram.'
            end if

            ! Initialise start of binary searches.
            if (iter >= nHistEquilSteps) &
                HistMinInd(1:Nel) = FCIDetIndex(1:nel)
        end if

    end subroutine

    subroutine end_iter_stats (TotWalkersNew)

        implicit none
        integer, intent(in) :: TotWalkersNew
        integer :: proc, pos, i, k
        real(dp) :: sgn(lenof_sign)
        integer :: run

        ! SumWalkersCyc calculates the total number of walkers over an update
        ! cycle on each process. (used for shift update)
        call add_part_number(SumWalkersCyc)
        ! SumWalkersOut is the total walker number over an output cycle (used for acc. rate)
        call add_part_number(SumWalkersOut)

        ! Write initiator histograms if on the correct iteration.
        ! Why is this done here - before annihilation!
        if ((tHistInitPops .and. mod(iter, HistInitPopsIter) == 0) &
            .or. tPrintHighPop) then
            call FindHighPopDet (TotWalkersNew)
            if (tHistInitPops) then
                root_write(iout,'(a)') 'Writing out the spread of the &
                                       &initiator determinant populations.'
                call WriteInitPops (iter + PreviousCycles)
            end if
         endif

       contains

         subroutine add_part_number(var)
           real(dp), intent(inout) :: var(inum_runs)

#ifdef CMPLX_
           do run = 1, inum_runs
              var(run) = var(run) + sum(TotParts(min_part_type(run):max_part_type(run)))
           enddo
#else
           var = var + TotParts
#endif           
         end subroutine add_part_number

    end subroutine end_iter_stats


    logical function TestIfDETinCASBit (ilutnI)

        ! In:
        !    iLutNI: bit string representation of a determinant.
        ! Returns:
        !    true if the determinant is in the complete active space.

        integer(n_int), intent(in) :: iLutnI(0:NIfD)

        ! A determinant is in the CAS iff
        !  a) all orbitals in the core space are occupied;
        !  b) no orbitals in the external space are occupied;
        ! Thus ANDing the determinant with CASMask (containing set bits for
        ! the core and external orbitals) will give precisely the core
        ! orbitals if the determinant is in the CAS.

        TestifDETinCASBit = all(iand(iLutNI,CASMask) == CoreMask)

    end function TestIfDETinCASBit

    LOGICAL FUNCTION TestifDETinCAS(CASDet)
        INTEGER :: k,z,CASDet(NEl), orb
        LOGICAL :: tElecInVirt, bIsCsf

        ! CASmax is the max spin orbital number (when ordered energetically)
        ! within the chosen active space. Spin orbitals with energies larger
        ! than this maximum value must be unoccupied for the determinant to
        ! be in the active space.
!        CASmax=NEl+VirtCASorbs

        ! CASmin is the max spin orbital number below the active space.  As
        ! well as the above criteria, spin orbitals with energies equal to,
        ! or below that of the CASmin orbital must be completely occupied for
        ! the determinant to be in the active space.
        !
        ! (These have been moved to the InitCalc subroutine so they're not
        !  calculated each time).
!        CASmin=NEl-OccCASorbs

        bIsCsf = iscsf(CASDet)

        z=0
        tElecInVirt=.false.
        do k=1,NEl      ! running over all electrons
            ! TODO: is it reasonable to just apply the orbital mask anyway?
            !       it is probably faster than running iscsf...
            if (bIsCsf) then
                orb = iand(CASDet(k), csf_orbital_mask)
            else
                orb = CASDet(k)
            endif

            if (SpinInvBRR(orb).gt.CASmax) THEN
                tElecInVirt=.true.
                EXIT
                ! if at any stage an electron has an energy greater than the
                ! CASmax value, the determinant can be ruled out of the active
                ! space.  Upon identifying this, it is not necessary to check
                ! the remaining electrons.
            else
                if (SpinInvBRR(orb).le.CASmin) THEN
                    z=z+1
                endif
                ! while running over all electrons, the number that occupy
                ! orbitals equal to or below the CASmin cutoff are counted.
            endif
        enddo

        if(tElecInVirt.or.(z.ne.CASmin)) THEN
            ! if an electron is in an orbital above the active space, or the
            ! inactive orbitals are not full, the determinant is automatically
            ! ruled out.
            TestifDETinCAS=.false.
        else
            ! no orbital in virtual and the inactive orbitals are completely
            ! full - det in active space.
            TestifDETinCAS=.true.
        endif

        RETURN

    END FUNCTION TestifDETinCAS



    SUBROUTINE FindHighPopDet(TotWalkersNew)

        ! Found the highest population on each processor, need to find out
        ! which of these has the highest of all.

        INTEGER(n_int) :: DetPos(0:NIfTot),DetNeg(0:NIfTot)
        INTEGER :: TotWalkersNew
        real(dp) :: HighPopInNeg(2),HighPopInPos(2),HighPopoutNeg(2),HighPopoutPos(2)
        real(dp) :: TempSign(lenof_sign)

!        WRITE(6,*) 'HighPopPos',HighPopPos
!        WRITE(6,*) 'CurrentSign(HighPopPos)',CurrentSign(HighPopPos)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopNeg),TempSign)
        ELSE
            TempSign(:)=0.0_dp
        ENDIF

        HighPopInNeg(1) = TempSign(1)
        ! [W.D. 15.5.2017:]
        ! why cast to int32?? and not real(dp)
!         HighPopInNeg(2)= int(iProcIndex,int32)
        HighPopInNeg(2)= real(iProcIndex,dp)

        CALL MPIAllReduceDatatype(HighPopinNeg,1,MPI_MINLOC,MPI_2DOUBLE_PRECISION,HighPopoutNeg)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopPos),TempSign)
        ELSE
            TempSign(:)=0.0_dp
        ENDIF

        HighPopInPos(1) = TempSign(1)
        ! [W.D. 15.5.2017:]
        ! why cast to int32?? and not real(dp)
!         HighPopInPos(2)=int(iProcIndex,int32)
        HighPopInPos(2)=real(iProcIndex,dp)

        CALL MPIAllReduceDatatype(HighPopinPos,1,MPI_MAXLOC,MPI_2DOUBLE_PRECISION,HighPopoutPos)

        ! Now, on all processors, HighPopoutPos(1) is the highest positive
        ! population, and HighPopoutNeg(1) is the highest negative population.
        ! HighPopoutPos(2) is the processor the highest population came from.

        if (abs(iProcIndex - HighPopOutNeg(2)) < 1.0e-12_dp) DetNeg(:)=CurrentDets(:,HighPopNeg)
        if (abs(iProcIndex - HighPopOutPos(2)) < 1.0e-12_dp) DetPos(:)=CurrentDets(:,HighPopPos)

        ! This is a horrible hack, because the process argument should be of
        ! type 'integer' - whatever that is, but the highpopoutneg is
        ! explicitly an int(4), so that it works with MPI_2INTEGER. Because
        ! of the explicit interfaces, we need to do this.
        CALL MPIBcast(DetNeg ,NIfTot+1, int(HighPopOutNeg(2)))
        CALL MPIBcast(DetPos, NIfTot+1, int(HighPopOutPos(2)))

        if (iProcIndex == 0) then
            write (iout, '(a,f12.5,a)') 'The most highly populated determinant &
                                  & with the opposite sign to the HF has ', &
                                  HighPopoutNeg(1), ' walkers.'
            call WriteBitDet (iout, DetNeg, .true.)

            write (iout,'(a,f12.5,a9)') 'The most highly populated determinant &
                                  & with the same sign as the HF has ', &
                                  HighPopoutPos(1), ' walkers.'
            call WriteBitDet (iout, DetPos, .true.)
        endif

        tPrintHighPop=.false.


    END SUBROUTINE FindHighPopDet


    function CheckAllowedTruncSpawn (WalkExcitLevel, nJ, ilutnJ, IC) &
                                    result(bAllowed)
        ! Under any currently applied truncation schemes, is an excitation to
        ! this determinant allowed?
        !
        ! In:  WalkExcitLevel - Current excitation level relative to HF
        !      nJ             - Natural integer representation of det
        !                       (not Needed for HPHF/tTruncNOpen/MomInv)
        !      ilutnJ         - Bit representation of det
        !      IC             - Excitation level relative to parent
        ! Ret: bAllowed       - .true. if excitation is allowed

        integer, intent(in) :: nJ(nel), WalkExcitLevel, IC
        integer(n_int), intent(in) :: ilutnJ(0:NIfTot)
        logical :: bAllowed

        integer :: NoInFrozenCore, MinVirt, ExcitLevel, i
        integer :: k(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "CheckAllowedTruncSpawn"
#endif

        bAllowed = .true.

        ! Truncate space by excitation level
        if (tTruncSpace) then
            ! If parent walker is one below excitation cutoff, could be
            ! disallowed if double. If higher, then all excits could
            ! be disallowed. If HPHF, excit could be single or double,
            ! and IC not returned --> Always test.
            if (tHPHF .or. WalkExcitLevel >= ICILevel .or. &
                (WalkExcitLevel == (ICILevel-1) .and. IC == 2)) then
                ExcitLevel = FindBitExcitLevel (iLutHF, ilutnJ, ICILevel, .true.)
                if (ExcitLevel > ICILevel) &
                    bAllowed = .false.
            endif
        endif

        ! Is the number of unpaired electrons too high?
        if (tTruncNOpen .and. bAllowed) then
            if (count_open_orbs(ilutnJ) > trunc_nopen_max) &
                bAllowed = .false.
        endif


        ! If the FCI space is restricted by a predetermined CAS space
        if (tTruncCAS .and. .not. tTruncInitiator .and. bAllowed) then
            if (.not. TestIfdetinCASBit(ilutnJ(0:NIfD))) &
                bAllowed = .false.
        endif


        ! Does the spawned determinant have more than the restricted number
        ! of holes in the partially frozen core?
        !
        ! --> Run through the e- in nJ, count the number in the partially
        !     frozen core (i.e. with energy, from BRR, less than the frozen
        !     core limit). If too few, then forbidden.
        if (tPartFreezeCore .and. bAllowed) then
            NoInFrozenCore = 0
            bAllowed = .false.
            do i = 1, nel
                if (SpinInvBRR(nJ(i)) <= NPartFrozen) &
                    NoInFrozenCore = NoInFrozenCore + 1
                if (NoInFrozenCore == (NPartFrozen - NHolesFrozen)) then
                    bAllowed = .true.
                    exit
                endif
            enddo
        endif


        ! Does the spawned determinant have more than the restricted number
        ! of e- in the partially frozen virtual orbitals?
        !
        ! --> Run through the e- in nJ, count the number in the partially
        !     frozen orbitals (i.e. with energy, from BRR, greater than
        !     minumum unfrozen virtual). If too many, then forbidden
        if (tPartFreezeVirt .and. bAllowed) then
            NoInFrozenCore = 0
            MinVirt = nBasis - NVirtPartFrozen
            ! BRR(i) = j: orbital i is the j-th lowest in energy
            do i = 1, nel
                if (SpinInvBRR(nJ(i)) > MinVirt) &
                    NoInFrozenCore = NoInFrozenCore + 1
                if (NoInFrozenCore > NElVirtFrozen) then
                    ! Too many e- in part-frozen orbs
                    bAllowed = .false.
                    exit
                endif
            enddo
        endif


        ! Check to see if UEG excitation is allowed, by summing kx, ky, kz
        ! over all the electrons
        if (tUEG .and. .not. tLatticeGens .and. bAllowed) then
            k = 0
            do i = 1, nel
                k = k + G1(nJ(i))%k
            enddo
            if (.not. all(k == 0)) &
                bAllowed = .false.
        endif


      end function CheckAllowedTruncSpawn
      
!Routine which takes a set of determinants and returns the ground state energy
    subroutine LanczosFindGroundE(Dets,DetLen,GroundE,ProjGroundE,tInit)
        use DetCalcData, only : NKRY,NBLK,B2L,nCycle
        implicit none
        real(dp) , intent(out) :: GroundE,ProjGroundE
        logical , intent(in) :: tInit
        integer, intent(in) :: DetLen
        integer, intent(in) :: Dets(NEl,DetLen)
        integer :: gc,ierr,icmax,lenhamil,nKry1,nBlock,LScr,LIScr
        integer(n_int) :: ilut(0:NIfTot)
        logical :: tmc,tSuccess
        real(dp) , allocatable :: A_Arr(:,:),V(:),AM(:),BM(:),T(:),WT(:),SCR(:),WH(:),WORK2(:),V2(:,:),W(:)
        HElement_t(dp), allocatable :: Work(:)
        HElement_t(dp), allocatable :: CkN(:,:),Hamil(:),TruncWavefunc(:)
        HElement_t(dp), allocatable :: Ck(:,:)  !This holds eigenvectors in the end
        HElement_t(dp) :: Num,Denom,HDiagTemp
        integer , allocatable :: LAB(:),NROW(:),INDEX(:),ISCR(:)
        integer(TagIntType) :: LabTag=0,NRowTag=0,ATag=0,VTag=0,AMTag=0,BMTag=0,TTag=0,tagCKN=0,tagCK=0
        integer(TagIntType) :: WTTag=0,SCRTag=0,ISCRTag=0,INDEXTag=0,WHTag=0,Work2Tag=0,V2Tag=0,tagW=0
        integer(TagIntType) :: tagHamil=0,WorkTag=0
        integer :: nEval,nBlocks,nBlockStarts(2),ExcitLev,iGetExcitLevel,i,nDoubles
        character(*), parameter :: t_r='LanczosFindGroundE'
        real(dp) :: norm
        integer :: PartInd,ioTrunc
        character(24) :: abstr

        if(DetLen.gt.1300) then
            nEval = 3
        else
            nEval = DetLen
        endif

        write(iout,'(A,I10)') 'DetLen = ',DetLen
!C..
        Allocate(Ck(DetLen,nEval), stat=ierr)
        call LogMemAlloc('CK',DetLen*nEval, 8,t_r, tagCK,ierr)
        CK=(0.0_dp)
!C..
        allocate(W(nEval), stat=ierr)
        call LogMemAlloc('W', nEval,8,t_r,tagW,ierr)
        W=0.0_dp
!        write(iout,*) "Calculating H matrix"
!C..We need to measure HAMIL and LAB first
        allocate(NROW(DetLen),stat=ierr)
        call LogMemAlloc('NROW',DetLen,4,t_r,NROWTag,ierr)
        NROW(:)=0
        ICMAX=1
        TMC=.FALSE.
        call DETHAM(DetLen,NEL,Dets,HAMIL,LAB,NROW,.TRUE.,ICMAX,GC,TMC)
!        WRITE(iout,*) ' FINISHED COUNTING '
        WRITE(iout,*) "Allocating memory for hamiltonian: ",GC*2
        CALL neci_flush(iout)
!C..Now we know size, allocate memory to HAMIL and LAB
        LENHAMIL=GC
        Allocate(Hamil(LenHamil), stat=ierr)
        call LogMemAlloc('HAMIL', LenHamil, 8, t_r,tagHamil,ierr)
        HAMIL=(0.0_dp)
!C..
        ALLOCATE(LAB(LENHAMIL),stat=ierr)
        CALL LogMemAlloc('LAB',LenHamil,4,t_r,LabTag,ierr)

        LAB(1:LENHAMIL)=0
!C..Now we store HAMIL and LAB
        CALL DETHAM(DetLen,NEL,Dets,HAMIL,LAB,NROW,.FALSE.,ICMAX,GC,TMC)

        if(DetLen.gt.1300) then
            !Lanczos diag

            Allocate(CkN(DetLen,nEval), stat=ierr)
            call LogMemAlloc('CKN',DetLen*nEval, 8,t_r, tagCKN,ierr)
            CKN=(0.0_dp)

            NKRY1=NKRY+1
            NBLOCK=MIN(NEVAL,NBLK)
            LSCR=MAX(DetLen*NEVAL,8*NBLOCK*NKRY)
            LISCR=6*NBLOCK*NKRY
    !C..
    !       write (iout,'(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
    !C..Set up memory for FRSBLKH

            ALLOCATE(A_Arr(NEVAL,NEVAL),stat=ierr)
            CALL LogMemAlloc('A_Arr',NEVAL**2,8,t_r,ATag,ierr)
            A_Arr=0.0_dp
    !C..
    !C,, W is now allocated with CK
    !C..
            ALLOCATE(V(DetLen*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('V',DetLen*NBLOCK*NKRY1,8,t_r,VTag,ierr)
            V=0.0_dp
    !C..
            ALLOCATE(AM(NBLOCK*NBLOCK*NKRY1),stat=ierr)
            CALL LogMemAlloc('AM',NBLOCK*NBLOCK*NKRY1,8,t_r,AMTag,ierr)
            AM=0.0_dp
    !C..
            ALLOCATE(BM(NBLOCK*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('BM',NBLOCK*NBLOCK*NKRY,8,t_r,BMTag,ierr)
            BM=0.0_dp
    !C..
            ALLOCATE(T(3*NBLOCK*NKRY*NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('T',3*NBLOCK*NKRY*NBLOCK*NKRY,8,t_r,TTag,ierr)
            T=0.0_dp
    !C..
            ALLOCATE(WT(NBLOCK*NKRY),stat=ierr)
            CALL LogMemAlloc('WT',NBLOCK*NKRY,8,t_r,WTTag,ierr)
            WT=0.0_dp
    !C..
            ALLOCATE(SCR(LScr),stat=ierr)
            CALL LogMemAlloc('SCR',LScr,8,t_r,SCRTag,ierr)
            SCR=0.0_dp
            ALLOCATE(ISCR(LIScr),stat=ierr)
            CALL LogMemAlloc('IScr',LIScr,4,t_r,IScrTag,ierr)
            ISCR(1:LISCR)=0
            ALLOCATE(INDEX(NEVAL),stat=ierr)
            CALL LogMemAlloc('INDEX',NEVAL,4,t_r,INDEXTag,ierr)
            INDEX(1:NEVAL)=0
    !C..
            ALLOCATE(WH(DetLen),stat=ierr)
            CALL LogMemAlloc('WH',DetLen,8,t_r,WHTag,ierr)
            WH=0.0_dp
            ALLOCATE(WORK2(3*DetLen),stat=ierr)
            CALL LogMemAlloc('WORK2',3*DetLen,8,t_r,WORK2Tag,ierr)
            WORK2=0.0_dp
            ALLOCATE(V2(DetLen,NEVAL),stat=ierr)
            CALL LogMemAlloc('V2',DetLen*NEVAL,8,t_r,V2Tag,ierr)
            V2=0.0_dp
    !C..Lanczos iterative diagonalising routine
            if (t_non_hermitian) then
                call stop_all(t_r, &
                    "NECI_FRSBLKH not adapted for non-hermitian Hamiltonians!")
            end if
            CALL NECI_FRSBLKH(DetLen,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A_Arr,W,V,AM,BM,T,WT, &
             &  SCR,ISCR,INDEX,NCYCLE,B2L,.true.,.false.,.false.,.false.)

            !Eigenvalues may come out wrong sign - multiply by -1
            if(W(1).gt.0.0_dp) then
                GroundE = -W(1)
            else
                GroundE = W(1)
            endif

            deallocate(A_Arr,V,AM,BM,T,WT,SCR,WH,V2,CkN,Index,IScr)
            call LogMemDealloc(t_r,ATag)
            call LogMemDealloc(t_r,VTag)
            call LogMemDealloc(t_r,AMTag)
            call LogMemDealloc(t_r,BMTag)
            call LogMemDealloc(t_r,TTag)
            call LogMemDealloc(t_r,tagCKN)
            call LogMemDealloc(t_r,WTTag)
            call LogMemDealloc(t_r,SCRTag)
            call LogMemDealloc(t_r,ISCRTag)
            call LogMemDealloc(t_r,IndexTag)
            call LogMemDealloc(t_r,WHTag)
            call LogMemDealloc(t_r,V2Tag)
        else
            !Full diag
            allocate(Work(4*DetLen),stat=ierr)
            call LogMemAlloc('Work',4*DetLen,8,t_r,WorkTag,ierr)
            allocate(Work2(3*DetLen),stat=ierr)
            call logMemAlloc('Work2',3*DetLen,8,t_r,Work2Tag,ierr)
            nBlockStarts(1) = 1
            nBlockStarts(2) = DetLen+1
            nBlocks = 1
            if (t_non_hermitian) then
                call stop_all(t_r, &
                    "HDIAG_neci is not set up for non-hermitian Hamiltonians!")
            end if
            call HDIAG_neci(DetLen,Hamil,Lab,nRow,CK,W,Work2,Work,nBlockStarts,nBlocks)
            GroundE = W(1)
            deallocate(Work)
            call LogMemDealloc(t_r,WorkTag)
        endif

        !Calculate proje.
        nDoubles = 0
        Num = 0.0_dp
        do i=1,DetLen
            ExcitLev = iGetExcitLevel(HFDet,Dets(:,i),NEl)
            if((ExcitLev.eq.1).or.(ExcitLev.eq.2)) then
                HDiagTemp = get_helement(HFDet,Dets(:,i),ExcitLev)
                Num = Num + (HDiagTemp * CK(i,1))
                nDoubles = nDoubles + 1
            elseif(ExcitLev.eq.0) then
                Denom = CK(i,1)
            endif
        enddo
        ProjGroundE = (Num / Denom) + Hii
!        write(iout,*) "***", nDoubles

        if(tHistSpawn.and..not.tInit) then
            !Write out the ground state wavefunction for this truncated calculation.
            !This needs to be sorted, into the order given in the original FCI in DetCalc.
            allocate(TruncWavefunc(Det),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc error")
            TruncWavefunc(:) = 0.0_dp

            Norm = 0.0_dp
            do i=1,DetLen
            !Loop over iluts in instantaneous list.

                call EncodeBitDet(Dets(:,i),iLut)
                !Calculate excitation level (even if hf isn't in list)
                ExcitLev = FindBitExcitLevel(iLutHF,iLut,NEl)

                if(ExcitLev.eq.nel) then
                    call BinSearchParts2(iLut,FCIDetIndex(ExcitLev),Det,PartInd,tSuccess)
                elseif(ExcitLev.eq.0) then
                    PartInd = 1
                    tSuccess = .true.
                else
                    call BinSearchParts2(iLut,FCIDetIndex(ExcitLev),FCIDetIndex(ExcitLev+1)-1,PartInd,tSuccess)
                endif
                if(.not.tSuccess) then
                    call stop_all(t_r,"Error here as cannot find corresponding determinant in FCI expansion")
                endif
                TruncWavefunc(PartInd) = CK(i,1)

                !Check normalisation
                Norm = Norm + CK(i,1)**2.0_dp
            enddo
            Norm = sqrt(Norm)
            if(abs(Norm-1.0_dp).gt.1.0e-7_dp) then
                write(iout,*) "***",norm
                call warning_neci(t_r,"Normalisation not correct for diagonalised wavefunction!")
            endif

            !Now write out...
            abstr = 'TruncWavefunc-'//str(Iter)
            ioTrunc = get_free_unit()
            open(ioTrunc,file=abstr,status='unknown')
            do i=1,Det
                write(ioTrunc,"(I13,F25.16)") i,TruncWavefunc(i)/norm
            enddo
            close(ioTrunc)

            deallocate(TruncWavefunc)

        endif


        deallocate(Work2,W,Hamil,Lab,nRow,CK)
        call LogMemDealloc(t_r,Work2Tag)
        call LogMemDealloc(t_r,tagW)
        call LogMemDealloc(t_r,tagHamil)
        call LogMemDealloc(t_r,LabTag)
        call LogMemDealloc(t_r,NRowTag)
        call LogMemDealloc(t_r,tagCK)

    end subroutine LanczosFindGroundE

    subroutine FlipSign(part_type)

        ! This routine flips the sign of all particles on the node with a
        ! given particle type. This allows us to keep multiple parallel
        ! simulations sign coherent.

        integer, intent(in) :: part_type
        character(*), parameter :: t_r = 'FlipSign'
        integer :: i
        real(dp) :: sgn

        do i = 1, int(TotWalkers, sizeof_int)

            sgn = extract_part_sign(CurrentDets(:,i), part_type)
            sgn = -sgn
            call encode_part_sign(CurrentDets(:,i), sgn, part_type)

            ! Flip average signs too.
            if (tFillingStochRDMOnFly) then
                if (lenof_sign /= 1) &
                    call stop_all(t_r, 'Not yet implemented')
                call set_av_sgn_tot(i, -get_av_sgn_tot(i))
            end if

        enddo

        ! Reverse the flag for whether the sign of the particles has been
        ! flipped so the ACF can be correctly calculated
        tFlippedSign = .not. tFlippedSign

    end subroutine FlipSign

!This routine takes the walkers from all subspaces, constructs the hamiltonian, and diagonalises it.
!Currently, this only works in serial.
    subroutine DiagWalkerSubspace()
        implicit none
        integer :: i,iSubspaceSize,ierr,iSubspaceSizeFull
        real(dp) :: CurrentSign(lenof_sign)
        integer, allocatable :: ExpandedWalkerDets(:,:)
        integer(TagIntType) :: ExpandedWalkTag=0
        real(dp) :: GroundEFull,GroundEInit,CreateNan,ProjGroundEFull,ProjGroundEInit
        character(*), parameter :: t_r='DiagWalkerSubspace'

        if(nProcessors.gt.1) call stop_all(t_r,"Walker subspace diagonalisation only works in serial")
        if(lenof_sign.ne.1) call stop_all(t_r,'Cannot do Lanczos on complex orbitals.')

        if(tTruncInitiator) then
            !First, diagonalise initiator subspace
            write(iout,'(A)') 'Diagonalising initiator subspace...'

            iSubspaceSize = 0
            do i=1,int(TotWalkers,sizeof_int)
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if((abs(CurrentSign(1)) > InitiatorWalkNo) .or. &
                        (DetBitEQ(CurrentDets(:,i),iLutHF,NIfDBO))) then
                    !Is allowed initiator. Add to subspace.
                    iSubspaceSize = iSubspaceSize + 1
                endif
            enddo

            write(iout,'(A,I12)') "Number of initiators found to diagonalise: ",iSubspaceSize
            allocate(ExpandedWalkerDets(NEl,iSubspaceSize),stat=ierr)
            call LogMemAlloc('ExpandedWalkerDets',NEl*iSubspaceSize,4,t_r,ExpandedWalkTag,ierr)

            iSubspaceSize = 0
            do i=1,int(TotWalkers,sizeof_int)
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if((abs(CurrentSign(1)) > InitiatorWalkNo) .or. &
                        (DetBitEQ(CurrentDets(:,i),iLutHF,NIfDBO))) then
                    !Is allowed initiator. Add to subspace.
                    iSubspaceSize = iSubspaceSize + 1
                    call decode_bit_det(ExpandedWalkerDets(:,iSubspaceSize),CurrentDets(:,i))
                endif
            enddo

            if(iSubspaceSize.gt.0) then
!Routine to diagonalise a set of determinants, and return the ground state energy
                call LanczosFindGroundE(ExpandedWalkerDets,iSubspaceSize,GroundEInit,ProjGroundEInit,.true.)
                write(iout,'(A,G25.10)') 'Ground state energy of initiator walker subspace = ',GroundEInit
            else
                CreateNan=-1.0_dp
                write(iout,'(A,G25.10)') 'Ground state energy of initiator walker subspace = ',sqrt(CreateNan)
            endif


            deallocate(ExpandedWalkerDets)
            call LogMemDealloc(t_r,ExpandedWalkTag)

        endif

        iSubspaceSizeFull = int(TotWalkers,sizeof_int)

        !Allocate memory for walker list.
        write(iout,'(A)') "Allocating memory for diagonalisation of full walker subspace"
        write(iout,'(A,I12,A)') "Size = ",iSubspaceSizeFull," walkers."

        allocate(ExpandedWalkerDets(NEl,iSubspaceSizeFull),stat=ierr)
        call LogMemAlloc('ExpandedWalkerDets',NEl*iSubspaceSizeFull,4,t_r,ExpandedWalkTag,ierr)
        do i=1,iSubspaceSizeFull
            call decode_bit_det(ExpandedWalkerDets(:,i),CurrentDets(:,i))
        enddo

!Routine to diagonalise a set of determinants, and return the ground state energy
        call LanczosFindGroundE(ExpandedWalkerDets,iSubspaceSizeFull,GroundEFull,ProjGroundEFull,.false.)

        write(iout,'(A,G25.10)') 'Ground state energy of full walker subspace = ',GroundEFull

        if(tTruncInitiator) then
            write(unitWalkerDiag,'(3I14,4G25.15)') Iter,iSubspaceSize,iSubspaceSizeFull,GroundEInit-Hii,    &
                    GroundEFull-Hii,ProjGroundEInit-Hii,ProjGroundEFull-Hii
        else
            write(unitWalkerDiag,'(2I14,2G25.15)') Iter,iSubspaceSizeFull,GroundEFull-Hii,ProjGroundEFull-Hii
        endif

        deallocate(ExpandedWalkerDets)
        call LogMemDealloc(t_r,ExpandedWalkTag)

    end subroutine DiagWalkerSubspace


    subroutine decide_num_to_spawn(parent_pop, av_spawns_per_walker, nspawn)

        real(dp), intent(in) :: parent_pop
        real(dp), intent(in) :: av_spawns_per_walker
        integer, intent(out) :: nspawn
        real(dp) :: prob_extra_walker, r

        nspawn = abs(int(parent_pop*av_spawns_per_walker))
        if (abs(abs(parent_pop*av_spawns_per_walker) - real(nspawn,dp)) > 1.e-12_dp) then
            prob_extra_walker = abs(parent_pop*av_spawns_per_walker) - real(nspawn,dp)
            r = genrand_real2_dSFMT()
            if (prob_extra_walker > r) nspawn = nspawn + 1
        end if

    end subroutine decide_num_to_spawn

    subroutine rescale_spawns(ValidSpawned, proj_energy, iter_data)

        integer, intent(in) :: ValidSpawned
        real(dp), intent(in) :: proj_energy(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data

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

            iter_data%nborn = iter_data%nborn + abs(spwnsign)
        end do

    end subroutine rescale_spawns

    subroutine set_init_flag_spawns_to_occ(ValidSpawned)

        ! Loop through the SpawnedParts array and set the initiator flag for
        ! any spawnings to determinants already occupied in CurrenDets.

        ! Usually this is done in AnnihilateSpawnedParts, but with
        ! preconditioning and a time step of exactly 1, all walkers are
        ! killed and removed from CurrenDets before then.

        ! IMPORTANT: This should only be used after spawnings have been
        ! sent to their parent process. And preferably should not be
        ! called until repeated spawnings ahve been compressed, for the
        ! sake of efficiency.

        integer, intent(in) :: ValidSpawned

        integer :: i, j, PartInd, DetHash
        integer :: nI_spawn(nel)
        real(dp) :: cursign(lenof_sign)
        logical :: tSuccess

        do i = 1, ValidSpawned
            call decode_bit_det(nI_spawn, SpawnedParts(:,i))

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
            end if

        end do

    end subroutine set_init_flag_spawns_to_occ

    subroutine perform_death_all_walkers(iter_data)

        use DetBitOps, only: FindBitExcitLevel
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

    end subroutine perform_death_all_walkers

    subroutine walker_death (iter_data, DetCurr, iLutCurr, Kii, RealwSign, &
                             DetPosition, walkExcitLevel)

        use global_det_data, only: get_iter_occ_tot, get_av_sgn_tot
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: one_rdms, two_rdm_spawn, rdm_definitions, &
                    inits_one_rdms, two_rdm_inits_spawn
        use semi_stoch_procs, only: check_determ_flag

        integer, intent(in) :: DetCurr(nel)
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        real(dp), intent(in) :: Kii
        integer, intent(in) :: DetPosition
        type(fcimc_iter_data), intent(inout) :: iter_data

        real(dp) :: iDie(lenof_sign), CopySign(lenof_sign)
        real(dp) :: av_sign(len_av_sgn_tot), iter_occ(len_iter_occ_tot)
        integer, intent(in) :: walkExcitLevel
        integer :: i, irdm, run
        logical :: tCoreDet
        character(len=*), parameter :: t_r = "walker_death"

        ! Do particles on determinant die? iDie can be both +ve (deaths), or
        ! -ve (births, if shift > 0)
        iDie = attempt_die (DetCurr, Kii, realwSign, WalkExcitLevel, DetPosition)

        IFDEBUG(FCIMCDebug,3) then
            if (sum(abs(iDie)) > 1.0e-10_dp) then
                write(iout,"(A)",advance='no') "Death: "
                do i = 1,lenof_sign-1
                    write(iout,"(f10.5)",advance='no') iDie(i)
                enddo
                write(iout,"(f10.5)") iDie(i)
            endif
        endif

        ! Update death counter
        iter_data%ndied = iter_data%ndied + min(iDie, abs(RealwSign))
#ifdef CMPLX_
        do run = 1, inum_runs
            NoDied(run) = NoDied(run) &
                + sum(min(iDie(min_part_type(run):max_part_type(run)), abs(RealwSign(min_part_type(run):max_part_type(run)) )))
        enddo
#else
        NoDied = NoDied + min(iDie, abs(RealwSign))
#endif

        ! Count any antiparticles
        iter_data%nborn = iter_data%nborn + max(iDie - abs(RealwSign), 0.0_dp)
#ifdef CMPLX_
        do run = 1, inum_runs
            NoBorn(run) = NoBorn(run) &
                + sum(max(iDie(min_part_type(run):max_part_type(run)) &
                - abs(RealwSign(min_part_type(run):max_part_type(run))), 0.0_dp))
        enddo
#else
        NoBorn = NoBorn + max(iDie - abs(RealwSign), 0.0_dp)
#endif

        ! Calculate new number of signed particles on the det.
        CopySign = RealwSign - (iDie * sign(1.0_dp, RealwSign))

        ! In the initiator approximation, abort any anti-particles.
        if (tTruncInitiator .and. any(abs(CopySign) > 1.0e-12_dp)) then
            do i = 1, lenof_sign
                if (CopySign(i) > 0.0_dp .neqv. RealwSign(i) > 0.0_dp) then
                    NoAborted(i) = NoAborted(i) + abs(CopySign(i))
                    iter_data%naborted(i) = iter_data%naborted(i) &
                                          + abs(CopySign(i))
                    if (test_flag(ilutCurr, get_initiator_flag(i))) &
                        NoAddedInitiators = NoAddedInitiators - 1_int64
                    CopySign(i) = 0
                end if
            end do
        end if

        tCoreDet = check_determ_flag(iLutCurr)

        if (any(abs(CopySign) > 1.0e-12_dp) .or. tCoreDet) then
            ! For the hashed walker main list, the particles don't move.
            ! Therefore just adjust the weight.
            call encode_sign (CurrentDets(:, DetPosition), CopySign)
        else
            ! All walkers died.
            if (tFillingStochRDMonFly) then
                av_sign = get_av_sgn_tot(DetPosition)
                iter_occ = get_iter_occ_tot(DetPosition)
                call det_removed_fill_diag_rdm(two_rdm_spawn, one_rdms, CurrentDets(:,DetPosition), av_sign, iter_occ)
                if(tInitsRDM .and. all_runs_are_initiator(CurrentDets(:,DetPosition))) &
                     call det_removed_fill_diag_rdm(two_rdm_inits_spawn, inits_one_rdms, &
                     CurrentDets(:,DetPosition), av_sign, iter_occ, .false.)
                ! Set the average sign and occupation iteration to zero, so
                ! that the same contribution will not be added in in
                ! CalcHashTableStats, if this determinant is not overwritten
                ! before then
                global_determinant_data(:, DetPosition) = 0.0_dp
            end if

            if (tTruncInitiator) then
                ! All particles on this determinant have gone. If the determinant was an initiator, update the stats
                do i = 1, lenof_sign
                    if (test_flag(iLutCurr,get_initiator_flag(i))) then
                        NoAddedInitiators(i) = NoAddedInitiators(i) - 1_int64
                    endif
                enddo
            end if

            ! Remove the determinant from the indexing list
            call RemoveHashDet(HashIndex, DetCurr, DetPosition)
            ! Encode a null det to be picked up
            call encode_sign(CurrentDets(:,DetPosition), null_part)
        end if

        ! Test - testsuite, RDM still work, both still work with Linscalealgo (all in debug)
        ! Null particle not kept if antiparticles aborted.
        ! When are the null particles removed?

    end subroutine walker_death

    subroutine check_start_rdm()

        ! This routine checks if we should start filling the RDMs -
        ! and does so if we should.

        use rdm_general, only: realloc_SpawnedParts
        use LoggingData, only: tReadRDMs, tTransitionRDMs
        implicit none
        logical :: tFullVaryshift

        tFullVaryShift=.false.

        if (all(.not. tSinglePartPhase)) tFullVaryShift = .true.

        ! If we're reading in the RDMs we've already started accumulating them in a previous calculation
        ! We don't want to put in an arbitrary break now!
        if (tReadRDMs) IterRDMonFly = 0

        if (tFullVaryShift .and. ((Iter - maxval(VaryShiftIter)).eq.(IterRDMonFly+1))) then
        ! IterRDMonFly is the number of iterations after the shift has changed that we want
        ! to fill the RDMs.  If this many iterations have passed, start accumulating the RDMs!

            IterRDMStart = Iter + PreviousCycles
            IterRDM_HF = Iter + PreviousCycles

            if (tEN2) tEN2Started = .true.

            ! We have reached the iteration where we want to start filling the RDM.
            if (tExplicitAllRDM) then
                ! Explicitly calculating all connections - expensive...
                if (tPairedReplicas) call stop_all('check_start_rdm',"Cannot yet do replica RDM sampling with explicit RDMs. &
                    & e.g Hacky bit in Gen_Hist_ExcDjs to make it compile.")

                tFillingExplicRDMonFly = .true.
                if(tHistSpawn) NHistEquilSteps = Iter
            else

                ! If we are load balancing, this will disable the load balancer
                ! so we should do a last-gasp balance at this point.
                if (tLoadBalanceBlocks .and. .not. tSemiStochastic) &
                    call adjust_load_balance(iter_data_fciqmc)

                extract_bit_rep_avsign => extract_bit_rep_avsign_norm
                ! By default - we will do a stochastic calculation of the RDM.
                tFillingStochRDMonFly = .true.
                if (tTransitionRDMs) tTransitionRDMsStarted = .true.

                call realloc_SpawnedParts()
                ! The SpawnedParts array now needs to carry both the spawned parts Dj, and also it's
                ! parent Di (and it's sign, Ci). - We deallocate it and reallocate it with the larger size.
                ! Don't need any of this if we're just doing HF_Ref_Explicit calculation.
                ! This is all done in the add_rdm_hfconnections routine.
            end if

            if (RDMExcitLevel .eq. 1) then
                write(6,'(A)') 'Calculating the 1 electron density matrix on the fly.'
            else
                write(6,'(A)') 'Calculating the 2 electron density matrix on the fly.'
            end if
            write(6,'(A,I10)') 'Beginning to fill the RDMs during iteration', Iter
        ENDIF

    end subroutine check_start_rdm

    subroutine update_run_reference(ilut, run)
        use adi_references, only: update_first_reference
        ! Update the reference used for a particular run to the one specified.
        ! Update the HPHF flipped arrays, and adjust the stored diagonal
        ! energies to account for the change if necessary.
        use SystemData, only: BasisFn, nBasisMax
        use sym_mod, only: writesym, getsym
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: run
        character(*), parameter :: this_routine = 'update_run_reference'

        HElement_t(dp) :: h_tmp
        real(dp) :: old_hii
        integer :: i, det(nel)
        logical :: tSwapped
        Type(BasisFn) :: isym

        iLutRef(:, run) = 0_n_int
        iLutRef(0:NIfDBO, run) = ilut(0:NIfDBO)
        call decode_bit_det (ProjEDet(:, run), iLutRef(:, run))
        write (iout, '(a,i3,a)', advance='no') 'Changing projected &
              &energy reference determinant for run', run, &
              ' on the next update cycle to: '
        call write_det (iout, ProjEDet(:, run), .true.)
        call GetSym(ProjEDet(:,run),nEl,G1,nBasisMax,isym)
        write(6,"(A)",advance='no') " Symmetry: "
        call writeSym(6,isym%sym,.true.)

        if(tHPHF) then
            if(.not.Allocated(RefDetFlip)) then
                allocate(RefDetFlip(NEl, inum_runs), &
                         ilutRefFlip(0:NifTot, inum_runs))
                RefDetFlip = 0
                iLutRefFlip = 0_n_int
            endif
            if(.not. TestClosedShellDet(iLutRef(:, run))) then
                ! Complications. We are now effectively projecting
                ! onto a LC of two dets. Ensure this is done correctly.
                call ReturnAlphaOpenDet(ProjEDet(:,run), &
                                        RefDetFlip(:, run), &
                                        iLutRef(:,run), &
                                        iLutRefFlip(:, run), &
                                        .true., .true., tSwapped)
                if(tSwapped) then
                    ! The iLutRef should already be the correct
                    ! one, since it was obtained by the normal
                    ! calculation!
                    call stop_all(this_routine, &
                        "Error in changing reference determinant &
                        &to open shell HPHF")
                endif
                write(iout,"(A,i3)") "Now projecting onto open-shell &
                    &HPHF as a linear combo of two determinants...&
                    & for run", run
                tSpinCoupProjE(run) = .true.
            endif
        else
            ! In case it was already on, and is now projecting
            ! onto a CS HPHF.
            tSpinCoupProjE(run) = .false.
        endif

        ! We can't use Brillouin's theorem if not a converged,
        ! closed shell, ground state HF det.
        tNoBrillouin = .true.
        tRef_Not_HF = .true.
        root_print "Ensuring that Brillouin's theorem is no &
                   &longer used."

        ! If this is the first replica, update the global reference
        ! energy.
        if (run == 1) then

            old_Hii = Hii
            if(tZeroRef) then
               h_tmp = 0.0_dp
            else if (tHPHF) then
                h_tmp = hphf_diag_helement (ProjEDet(:,1), &
                                            iLutRef(:,1))
            else
                h_tmp = get_helement (ProjEDet(:,1), &
                                      ProjEDet(:,1), 0)
            endif
            Hii = real(h_tmp, dp)
            write (iout, '(a, g25.15)') &
                'Reference energy now set to: ', Hii

            ! Regenerate all the diagonal elements relative to the
            ! new reference det.
            write (iout,*) 'Regenerating the stored diagonal &
                           &HElements for all walkers.'
            do i = 1, int(Totwalkers,sizeof_int)
                call decode_bit_det (det, CurrentDets(:,i))
                if (tHPHF) then
                    h_tmp = hphf_diag_helement (det, &
                                                CurrentDets(:,i))
                else
                    h_tmp = get_helement (det, det, 0)
                endif
                call set_det_diagH(i, real(h_tmp, dp) - Hii)
            enddo
            if (tSemiStochastic) &
                call recalc_core_hamil_diag(old_Hii, Hii)

            if (tReplicaReferencesDiffer) then
                ! Ensure that the energy references for all of the runs are
                ! relative to the new Hii
                do i = 1, inum_runs
                    proje_ref_energy_offsets(i) = proje_ref_energy_offsets(i) &
                                                + old_hii - hii
                end do
            end if

            ! All of the shift energies are relative to Hii, so they need to
            ! be offset
            DiagSft = DiagSft + old_hii - hii

        end if ! run == 1

        ! Ensure that our energy offsets for outputting the correct
        ! data have been updated correctly.
        if (tHPHF) then
            h_tmp = hphf_diag_helement (ProjEDet(:,run), &
                                        ilutRef(:,run))
        else
            h_tmp = get_helement (ProjEDet(:,run), &
                                  ProjEDet(:,run), 0)
        endif
        proje_ref_energy_offsets(run) = real(h_tmp, dp) - Hii

        ! Update the processor on which the reference is held
        iRefProc(run) = DetermineDetNode(nel, ProjEDet(:, run), 0)

        ! [W.D] need to also change the virtual mask
        if (t_back_spawn .or. t_back_spawn_flex) then
            call setup_virtual_mask()
        end if

        ! Also update ilutRefAdi - this has to be done completely
        call update_first_reference()

      end subroutine update_run_reference

    subroutine calc_inst_proje()

        ! Calculate an instantaneous value of the projected energy for the
        ! given walkers distributions

        integer :: ex_level, det(nel), j, run
        real(dp), dimension(max(lenof_sign,inum_runs)) :: RealAllHFCyc
        real(dp) :: sgn(lenof_sign)

        ! Reset the accumulators
        HFCyc = 0.0_dp
        ENumCyc = 0.0_dp
        NoatDoubs = 0.0_dp

        ! Main loop
        do j = 1, int(TotWalkers, sizeof_int)

            ! n.b. non-contiguous list
            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

            ex_level = FindBitExcitLevel (iLutRef(:,1), CurrentDets(:,j), t_hphf_ic = .true.)

            call decode_bit_det(det, CurrentDets(:,j))
            call SumEContrib(det, ex_level, sgn, CurrentDets(:,j), 0.0_dp, &
                             1.0_dp, tPairedReplicas, j)
        end do

        ! Accumulate values over all processors
        call MPISumAll(HFCyc, RealAllHFCyc)
        call MPISumAll(ENumCyc, AllENumCyc)

        do run = 1, inum_runs
            AllHFCyc(run) = ARR_RE_OR_CPLX(RealAllHFCyc, run)
        end do

        proje_iter = AllENumCyc / AllHFCyc + proje_ref_energy_offsets


    end subroutine

    !> Set the offset of the adaptive shift equal to the ground state energy
    !> of the trial space.
    subroutine Set_AS_TrialOffset()


        integer :: run

        FullShiftOffset = trial_energies(1)
        !Find the lowest trial energy
        do run=2, inum_runs
            if(trial_energies(run)<FullShiftOffset) &
                FullShiftOffset = trial_energies(run)
        enddo
        tAS_Offset = .true.
        write(6,*) "The adaptive shift is offset by the correlation energy of trial-wavefunction: ", FullShiftOffset-Hii
    end subroutine

end module
