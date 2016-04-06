#include "macros.h"

module fcimc_helper

    use constants
    use util_mod
    use systemData, only: nel, tHPHF, tNoBrillouin, G1, tUEG, &
                          tLatticeGens, nBasis, tHistSpinDist, tRef_Not_HF
    use HPHFRandExcitMod, only: ReturnAlphaOpenDet
    use semi_stoch_procs, only: recalc_core_hamil_diag
    use bit_reps, only: NIfTot, flag_initiator, test_flag, extract_flags, &
                        encode_bit_rep, NIfD, set_flag_general, NIfDBO, &
                        extract_sign, set_flag, encode_sign, &
                        flag_trial, flag_connected, flag_deterministic, &
                        extract_part_sign, encode_part_sign, decode_bit_det, &
                        set_has_been_initiator, flag_has_been_initiator, &
                        set_parent_coeff, flag_weak_initiator
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
                           FciMCDebug
    use CalcData, only: NEquilSteps, tFCIMC, tTruncCAS, &
                        tAddToInitiator, InitiatorWalkNo, &
                        tTruncInitiator, tTruncNopen, trunc_nopen_max, &
                        tRealCoeffByExcitLevel, tSurvivalInitiatorThreshold, &
                        tSemiStochastic, tTrialWavefunction, DiagSft, &
                        InitiatorCutoffEnergy, InitiatorCutoffWalkNo, &
                        im_time_init_thresh, tSurvivalInitMultThresh, &
                        init_survival_mult, MaxWalkerBloom, &
                        tMultiReplicaInitiators, NMCyc, iSampleRDMIters, &
                        tSpawnCountInitiatorThreshold, init_spawn_thresh, &
                        tOrthogonaliseReplicas, tPairedReplicas, &
                        tBroadcastParentCoeff, tWeakInitiators, weakthresh
    use IntegralsData, only: tPartFreezeVirt, tPartFreezeCore, NElVirtFrozen, &
                             nPartFrozen, nVirtPartFrozen, nHolesFrozen
    use procedure_pointers, only: attempt_die, extract_bit_rep_avsign
    use DetCalcData, only: FCIDetIndex, ICILevel, det
    use hash, only: remove_hash_table_entry
    use load_balance_calcnodes, only: DetermineDetNode, tLoadBalanceBlocks
    use load_balance, only: adjust_load_balance
    use rdm_filling, only: det_removed_fill_diag_rdm
    use rdm_general, only: store_parent_with_spawned, extract_bit_rep_avsign_norm
    use Parallel_neci
    use FciMCLoggingMod, only: HistInitPopulations, WriteInitPops
    use csf_data, only: csf_orbital_mask
    use csf, only: iscsf
    use hphf_integrals, only: hphf_diag_helement
    use global_det_data, only: get_av_sgn, set_av_sgn, set_det_diagH, &
                               global_determinant_data, set_iter_occ, &
                               get_part_init_time, det_diagH, get_spawn_count
    use searching, only: BinSearchParts2
    use rdm_data, only: nrdms
    implicit none
    save

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


    subroutine create_particle (nJ, iLutJ, child, part_type, ilutI, SignCurr, &
                                WalkerNo, RDMBiasFacCurr, WalkersToSpawn)

        ! Create a child in the spawned particles arrays. We spawn particles
        ! into a separate array, but non-contiguously. The processor that the
        ! newly-spawned particle is going to be sent to has to be determined,
        ! and then it will be put into the appropriate element determined by
        ! ValidSpawnedList

        ! 'type' of the particle - i.e. real/imag

        integer, intent(in) :: nJ(nel), part_type
        integer(n_int), intent(in) :: iLutJ(0:niftot)
        real(dp), intent(in) :: child(lenof_sign)
        integer(n_int), intent(in), optional :: ilutI(0:niftot)
        real(dp), intent(in), optional :: SignCurr(lenof_sign)
        integer, intent(in), optional :: WalkerNo
        real(dp), intent(in), optional :: RDMBiasFacCurr
        integer, intent(in), optional :: WalkersToSpawn
        integer :: proc
        real(dp) :: r
        integer, parameter :: flags = 0
        logical :: list_full
        character(*), parameter :: this_routine = 'create_particle'

#ifdef __CMPLX
        integer :: j
        logical :: parent_init
#endif

        ! Determine which processor the particle should end up on in the
        ! DirectAnnihilation algorithm.
        proc = DetermineDetNode(nel,nJ,0)    ! (0 -> nNodes-1)

        ! Check that the position described by ValidSpawnedList is acceptable.
        ! If we have filled up the memory that would be acceptable, then
        ! kill the calculation hard (i.e. stop_all) with a descriptive
        ! error message.
        list_full = .false.
        if (proc == nNodes - 1) then
            if (ValidSpawnedList(proc) > MaxSpawned) list_full = .true.
        else
            if (ValidSpawnedList(proc) >= InitialSpawnedSlots(proc+1)) &
                list_full=.true.
        end if
        if (list_full) then
            write(6,*) "Attempting to spawn particle onto processor: ", proc
            write(6,*) "No memory slots available for this spawn."
            write(6,*) "Please increase MEMORYFACSPAWN"
            call stop_all(this_routine, "Out of memory for spawned particles")
        end if

        call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), iLutJ, &
                            child, flags)

        ! If the parent was an initiator then set the initiator flag for the
        ! child, to allow it to survive.
        if (tTruncInitiator) then
            if (test_flag(ilutI, flag_initiator(part_type)).or.test_flag(ilutI, flag_weak_initiator(part_type))) &
                call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), flag_initiator(part_type))
            if (tWeakInitiators) then
              if(test_flag(ilutI, flag_initiator(part_type))) then
                r = genrand_real2_dSFMT()
                        if(weakthresh > r) then
                             call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), flag_weak_initiator(part_type))
                        else
                             call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), flag_weak_initiator(part_type),.false.)
                        endif
              else
                call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), flag_weak_initiator(part_type),.false.)
              endif
            endif

        end if

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

        ! If we are storing the parent coefficient with the particle, then
        ! do that it this point
        if (tBroadcastParentCoeff) then
#ifdef __CMPLX
            ! n.b. SignCurr(part_type) --> this breaks with CPLX
            call stop_all(this_routine, 'Not implemented (yet)')
#endif
            call set_parent_coeff(SpawnedParts(:, ValidSpawnedList(proc)), &
                                  SignCurr(part_type))
        end if

#ifdef __CMPLX
        if (tTruncInitiator) then
            ! With complex walkers, things are a little more tricky.
            ! We want to transfer the flag for all particles created (both
            ! real and imag) from the specific type of parent particle. This
            ! can mean real walker flags being transfered to imaginary
            ! children and vice versa.
            ! This is unneccesary for real walkers.
            ! Test the specific flag corresponding to the parent, of type
            ! 'part_type'
            parent_init = test_flag(SpawnedParts(:,ValidSpawnedList(proc)), &
                                    flag_initiator(part_type))
            ! Assign this flag to all spawned children.
            do j=1,lenof_sign
                if (child(j) /= 0) then
                    call set_flag (SpawnedParts(:,ValidSpawnedList(proc)), &
                                   flag_initiator(j), parent_init)
                endif
            enddo
        end if
#endif

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        
        ! Sum the number of created children to use in acceptance ratio.
#ifdef __CMPLX
        acceptances(1) = acceptances(1) + sum(abs(child))
#else
        acceptances = acceptances + abs(child)
#endif

    end subroutine create_particle

    subroutine create_particle_with_hash_table (nI_child, ilut_child, child_sign, part_type, ilut_parent, iter_data)

        use hash, only: hash_table_lookup, add_hash_table_entry

        integer, intent(in) :: nI_child(nel), part_type
        integer(n_int), intent(in) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        real(dp), intent(in) :: child_sign(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data

        integer :: proc, ind, hash_val, i
        real(dp) :: real_sign_old(lenof_sign), real_sign_new(lenof_sign)
        real(dp) :: sgn_prod(lenof_sign)
        logical :: list_full, tSuccess
        integer, parameter :: flags = 0
        character(*), parameter :: this_routine = 'create_particle_with_hash_table'

        call hash_table_lookup(nI_child, ilut_child, NIfDBO, spawn_ht, SpawnedParts, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! If the spawned child is already in the spawning array.
            ! Extract the old sign.
            call extract_sign(SpawnedParts(:,ind), real_sign_old)

            ! If the new child has an opposite sign to that of walkers already
            ! on the site, then annihilation occurs. The stats for this need
            ! accumulating.
            sgn_prod = real_sign_old * child_sign
            do i = 1, lenof_sign
                if (sgn_prod(i) < 0.0_dp) then
                    iter_data%nannihil(i) = iter_data%nannihil(i) + 2*min( abs(real_sign_old(i)), abs(child_sign(i)) )
                end if
            end do

            ! Find the total new sign.
            real_sign_new = real_sign_old + child_sign
            ! Encode the new sign.
            call encode_sign(SpawnedParts(:,ind), real_sign_new)

            ! Set the initiator flags appropriately.
            ! If this determinant (on this replica) has already been spawned to
            ! then set the initiator flag. Also if this child was spawned from
            ! an initiator, set the initiator flag.
            if (tTruncInitiator) then
                if (abs(real_sign_old(part_type)) > 1.e-12_dp .or. test_flag(ilut_parent, flag_initiator(part_type))) &
                    call set_flag(SpawnedParts(:,ind), flag_initiator(part_type))
            end if
        else
            ! Determine which processor the particle should end up on in the
            ! DirectAnnihilation algorithm.
            proc = DetermineDetNode(nel, nI_child, 0)

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

            call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), ilut_child(0:NIfDBO), child_sign, flags)
            ! If the parent was an initiator then set the initiator flag for the
            ! child, to allow it to survive.
            if (tTruncInitiator) then
                if (test_flag(ilut_parent, flag_initiator(part_type))) &
                    call set_flag(SpawnedParts(:, ValidSpawnedList(proc)), flag_initiator(part_type))
            end if

            call add_hash_table_entry(spawn_ht, ValidSpawnedList(proc), hash_val)

            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        end if
        
        ! Sum the number of created children to use in acceptance ratio.
#ifdef __CMPLX
        acceptances(1) = acceptances(1) + sum(abs(child_sign))
#else
        acceptances = acceptances + abs(child_sign)
#endif

    end subroutine create_particle_with_hash_table

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

#ifdef __CMPLX
        complex(dp) :: CmplxwSign
#endif

        real(dp) :: amps(size(current_trial_amps,1))

        if (tReplicaReferencesDiffer) then
            call SumEContrib_different_refs(nI, realWSign, ilut, dProbFin, tPairedReplicas, ind)
            return
        endif

        HOffDiag = 0

        ! Add in the contributions to the numerator and denominator of the trial
        ! estimator, if it is being used.
#ifdef __CMPLX
        CmplxwSign = ARR_RE_OR_CPLX(realwsign, 1)

        if (tTrialWavefunction .and. present(ind)) then
            if (test_flag(ilut,flag_trial)) then
                if(ntrial_excits == 1) then
                   trial_denom = trial_denom + conjg(current_trial_amps(1,ind))*CmplxwSign
                else if(ntrial_excits == lenof_sign) then
                   call stop_all(this_routine, 'ntrial_excits has to be 1 currently for complex')
                end if

                if(qmc_trial_wf) then
                   call stop_all(this_routine, 'qmc_trial_wf currently not implemented for complex')
                end if
            else if (test_flag(ilut,flag_connected)) then
                if(ntrial_excits == 1) then
                   trial_numerator = trial_numerator + conjg(current_trial_amps(1,ind))*cmplxwsign
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
                else if (ntrial_excits == lenof_sign) then
                    trial_denom = trial_denom + current_trial_amps(:,ind)*RealwSign
                end if

                if (qmc_trial_wf) then
                    call get_con_amp_trial_space(ilut, amps)

                    if (ntrial_excits == 1) then
                        trial_numerator = trial_numerator + amps(1)*RealwSign
                    else if (ntrial_excits == lenof_sign) then
                        trial_numerator = trial_numerator + amps*RealwSign
                    end if
                end if

            else if (test_flag(ilut, flag_connected)) then
                ! Note, only attempt to add in a contribution from the
                ! connected space if we're not also in the trial space.
                 if (ntrial_excits == 1) then
                    trial_numerator = trial_numerator + current_trial_amps(1,ind)*RealwSign
                else if (ntrial_excits == lenof_sign) then
                    trial_numerator = trial_numerator + current_trial_amps(:,ind)*RealwSign
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

        ! Perform normal projection onto reference determinant
        if (ExcitLevel_local == 0) then

            if (iter > NEquilSteps) &
                SumNoatHF(1:lenof_sign) = SumNoatHF(1:lenof_sign) + RealwSign
            NoatHF(1:lenof_sign) = NoatHF(1:lenof_sign) + RealwSign
            ! Number at HF * sign over course of update cycle
            HFCyc(1:lenof_sign) = HFCyc(1:lenof_sign) + RealwSign

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
            
            if (ExcitLevel_local == 2) then
#ifdef __CMPLX
                NoatDoubs(1) = NoatDoubs(1) + sum(abs(RealwSign))
#else
                do run = 1, inum_runs
                    NoatDoubs(run) = NoatDoubs(run) + abs(RealwSign(run))
                end do
#endif
            end if

            ! Obtain off-diagonal element
            if (tHPHF) then
                HOffDiag(1:inum_runs) = hphf_off_diag_helement (ProjEDet(:,1), nI, &
                                                                iLutRef(:,1), ilut)
            else
                HOffDiag(1:inum_runs) = get_helement (ProjEDet(:,1), nI, &
                                                      ExcitLevel, ilutRef(:,1), ilut)
            endif

        endif ! ExcitLevel_local == 1, 2, 3


        ! Sum in energy contribution
        do run=1, inum_runs
            if (iter > NEquilSteps) &
                SumENum(run) = SumENum(run) + (HOffDiag(run) * ARR_RE_OR_CPLX(RealwSign,run)) &
                                  / dProbFin

            ENumCyc(run) = ENumCyc(run) + (HOffDiag(run) * ARR_RE_OR_CPLX(RealwSign,run)) / dProbFin
            ENumCycAbs(run) = ENumCycAbs(run) + abs(HoffDiag(run) * ARR_RE_OR_CPLX(RealwSign,run)) &
                                      / dProbFin
            
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
            if ((tPrintOrbOccInit .and. test_flag(ilut,flag_initiator(1)))&
                .or. .not. tPrintOrbOccInit) then
                forall (i = 1:nel) OrbOccs(iand(nI(i), csf_orbital_mask)) &
                        = OrbOccs(iand(nI(i), csf_orbital_mask)) &
                                   + (RealwSign(1) * RealwSign(1))
            endif
        endif
        
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
        real(dp) :: sgn_run
        HElement_t(dp) :: hoffdiag
        character(*), parameter :: this_routine = 'SumEContrib_different_refs'

        HElement_t(dp) :: amps(size(current_trial_amps,1))

        ASSERT(inum_runs == lenof_sign)
        ASSERT(tReplicaReferencesDiffer)
#ifdef __CMPLX
        call stop_all(this_routine, "Complex not supported")
#endif
        if (tHistSpawn .or. &
            (tCalcFCIMCPsi .and. tFCIMC) .or. tHistEnergies .or. &
            tHistSpinDist .or. tPrintOrbOcc) &
            call stop_all(this_routine, "Not yet supported")


        ! Add in the contributions to the numerator and denominator of the trial
        ! estimator, if it is being used.
        if (tTrialWavefunction .and. present(ind)) then
            if (test_flag(ilut, flag_trial)) then
                if (ntrial_excits == 1) then
                    trial_denom = trial_denom + current_trial_amps(1,ind)*sgn
                else
                    if (tPairedReplicas) then
#if defined(__PROG_NUMRUNS) || defined(__DOUBLERUN)
                        do run = 2, inum_runs, 2
                            trial_denom(run-1:run) = trial_denom(run-1:run) + current_trial_amps(run/2,ind)*sgn(run-1:run)
                        end do
#else
                        call stop_all(this_routine, "INVALID")
#endif
                    else
                        trial_denom = trial_denom + current_trial_amps(:,ind)*sgn
                    end if
                end if

                if (qmc_trial_wf) then
                    call get_con_amp_trial_space(ilut, amps)

                    if (ntrial_excits == 1) then
                        trial_numerator = trial_numerator + amps(1)*sgn
                    else
                        if (tPairedReplicas) then
#if defined(__PROG_NUMRUNS) || defined(__DOUBLERUN)
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
#if defined(__PROG_NUMRUNS) || defined(__DOUBLERUN)
                        do run = 2, inum_runs, 2
                            trial_numerator(run-1:run) = trial_numerator(run-1:run) + current_trial_amps(run/2,ind)*sgn(run-1:run)
                        end do
#else
                        call stop_all(this_routine, "INVALID")
#endif
                    else
                        trial_numerator = trial_numerator + current_trial_amps(:,ind)*sgn
                    end if
                end if
            end if
        end if


        ! This is the normal projected energy calculation, but split over
        ! multiple runs, rather than done in one go.
        do run = 1, inum_runs

            ! We need to use the excitation level relevant for this run
            exlevel = FindBitExcitLevel(ilut, ilutRef(:, run))
            if (tSpinCoupProjE(run) .and. exlevel /= 0) then
                if (exlevel <= 2) then
                    exlevel = 2
                else if (FindBitExcitLevel(ilut, ilutRefFlip(:,run)) <= 2) then
                    exlevel = 2
                end if
            end if
            sgn_run = sgn(run)

            hoffdiag = 0
            if (exlevel == 0) then

                if (iter > nEquilSteps) &
                    SumNoatHF(run) = SumNoatHF(run) + sgn_run
                NoatHF(run) = NoatHF(run) + sgn_run
                HFCyc(run) = HFCyc(run) + sgn_run

            else if (exlevel == 2 .or. (exlevel == 1 .and. tNoBrillouin)) then

                ! n.b. Brillouins theorem cannot hold for real-space Hubbard
                ! model or for rotated orbitals.

                if (exlevel == 2) &
                    NoatDoubs(run) = NoatDoubs(run) + sgn_run

                ! Obtain the off-diagonal elements
                if (tHPHF) then
                    hoffdiag = hphf_off_diag_helement(ProjEDet(:,run), nI, &
                                                      iLutRef(:,run), ilut)
                else
                    hoffdiag = get_helement (ProjEDet(:,run), nI, exlevel, &
                                             ilutRef(:,run), ilut)
                endif

            end if

            ! Sum in energy contributions
            if (iter > nEquilSteps) &
                SumENum(run) = SumENum(run) + (hoffdiag * sgn_run) / dProbFin
            ENumCyc(run) = ENumCyc(run) + (hoffdiag * sgn_run) / dProbFin
            ENumCycAbs(run) = ENumCycAbs(run) + abs(hoffdiag * sgn_run) / dProbFin

        end do

    end subroutine SumEContrib_different_refs


    subroutine CalcParentFlag(j, parent_flags, diagH)

        ! In the CurrentDets array, the flag at NIfTot refers to whether that
        ! determinant *itself* is an initiator or not. We need to decide if 
        ! this willchange due to the determinant acquiring a certain 
        ! population, or its population dropping below the threshold.
        ! The CurrentDets(:,j) is the determinant we are currently spawning 
        ! from, so this determines the ParentInitiator flag which is passed to
        ! the SpawnedDets array and refers to whether or not the walkers 
        ! *parent* is an initiator or not.

        integer, intent(in) :: j
        integer, intent(out) :: parent_flags
        real(dp) :: CurrentSign(lenof_sign)
        real(dp), intent(in) :: diagH
        integer :: part_type
        logical :: parent_init

        call extract_sign (CurrentDets(:,j), CurrentSign)

        ! The default path through this section makes no changes, leaving
        ! the initiator status of each parent unchanged.  If 
        ! tAddToInitiator is set, then the state of the parent may change.
        if (tAddToInitiator) then

            ! If we are considering initiators on a whole-site basis, then
            ! do that here.
            !
            ! n.b. these tests are not replicated in SortMerge.F90, as for
            !      there to be initiators on multiple particle types, there
            !      must be some merging happening.
            if (tMultiReplicaInitiators) then
                parent_init = test_flag (CurrentDets(:,j), flag_initiator(1))
                ! We sum the sign. Given that the reference site always has
                ! positive walkers on it, this is a sensible test.
                parent_init = TestInitiator(CurrentDets(:,j), parent_init, &
                                            sum(CurrentSign), diagH, &
                                            j, part_type)
            end if

            ! Now loop over the particle types, and update the flags
            do part_type = 1, lenof_sign

                if (.not. tMultiReplicaInitiators) then
                    ! By default, the parent_flags are the flags of the parent.
                    parent_init = test_flag (CurrentDets(:,j), flag_initiator(part_type))

                    ! Should this particle be considered to be an initiator
                    ! for spawning purposes.
                    parent_init = TestInitiator(CurrentDets(:,j), parent_init, &
                                                CurrentSign(part_type), diagH, &
                                                j, part_type)
                end if

                ! Update counters as required.
                if (parent_init) then
                    NoInitDets = NoInitDets + 1
                    NoInitWalk = NoInitWalk + abs(CurrentSign(part_type))
                else
                    NoNonInitDets = NoNonInitDets + 1
                    NoNonInitWalk = NoNonInitWalk + abs(CurrentSign(part_type))
                endif

                ! Update the parent flag as required.
                call set_flag (CurrentDets(:,j), flag_initiator(part_type), parent_init)

                if (parent_init) &
                    call set_has_been_initiator(CurrentDets(:,j), &
                                                flag_has_been_initiator(1))
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

    end subroutine CalcParentFlag


    function TestInitiator(ilut, is_init, sgn, diagH, site_idx, part_type) result(initiator)

        ! For a given particle (with its given particle type), should it
        ! be considered as an initiator for the purposes of spawning.
        !
        ! Inputs: The determinant, the particles sign, and if the particle
        !         is currently considered to be an initiator.

        ! N.B. This intentionally DOES NOT directly reference part_type.
        !      This means we can call it for individual, or aggregate,
        !      particles.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        logical, intent(in) :: is_init
        real(dp), intent(in) :: sgn, diagH
        integer, intent(in) :: site_idx, part_type
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'TestInitiator'
#endif

        logical :: initiator, tDetInCAS
        real(dp) :: init_thresh, low_init_thresh, init_tm, expected_lifetime
        real(dp) :: hdiag
        integer :: spwn_cnt, run

        ! By default the particles status will stay the same
        initiator = is_init

        ! Nice numbers
        init_thresh = InitiatorWalkNo
        low_init_thresh = InitiatorCutoffWalkNo
        run = part_type_to_run(part_type)

        if (.not. is_init) then

            ! Determinant wasn't previously initiator 
            ! - want to test if it has now got a large enough 
            !   population to become an initiator.
            if ((diagH > InitiatorCutoffEnergy &
                 .and. abs(sgn) > low_init_thresh) &
                .or. abs(sgn) > init_thresh) then
                initiator = .true.
                NoAddedInitiators = NoAddedInitiators + 1
            endif

        else

            ! The source determinant is already an initiator.            
            ! If tRetestAddToInit is on, the determinants become 
            ! non-initiators again if their population falls below 
            ! n_add (this is on by default).

            tDetInCas = .false.
            if (tTruncCAS) &
                tDetInCas = TestIfDetInCASBit (ilut)
           
            ! If det. in fixed initiator space, or is the HF det, or it
            ! is in the deterministic space, then it must remain an initiator.
            if (.not. tDetInCas .and. &
                .not. (DetBitEQ(ilut, iLutRef(:,run), NIfDBO)) &
                .and. .not. test_flag(ilut, flag_deterministic) &
                .and. ((diagH <= InitiatorCutoffEnergy .and. abs(sgn) <= init_thresh) .or. &
                       (diagH > InitiatorCutoffEnergy .and. abs(sgn) <= low_init_thresh))) then
                ! Population has fallen too low. Initiator status 
                ! removed.
                initiator = .false.
                NoAddedInitiators = NoAddedInitiators - 1
            endif

        end if

        ! If this site has survived for a long time, but otherwise
        ! would not be an initiator, then it is possible we ought
        ! to be considering it as well.
        if (.not. initiator .and. tSurvivalInitiatorThreshold) then
            init_tm = get_part_init_time(site_idx)
            if ((TotImagTime - init_tm) > im_time_init_thresh) &
                initiator = .true.
        end if

#ifdef __DEBUG
        if (tSurvivalInitiatorThreshold) then
            init_tm = get_part_init_time(site_idx)
            ASSERT(init_tm >= 0.0_dp)
        end if
#endif

        if (.not. initiator .and. tSurvivalInitMultThresh) then
            init_tm = get_part_init_time(site_idx)
            hdiag = det_diagH(site_idx) - DiagSft(part_type)
            if (hdiag > 0) then
                expected_lifetime = &
                    log(2.0_dp * max(MaxWalkerBloom, 1.0_dp)) / hdiag
                if ((TotImagTime - init_tm) > & !0.5_dp) then !&
                        init_survival_mult * expected_lifetime) then
                    initiator = .true.
                !    write(6,*) 'allowing', j, totimagtime-init_tm, &
                !        expected_lifetime, init_survival_mult*expected_lifetime, &
                !        init_survival_mult
                end if
            end if
        end if

        if (.not. initiator .and. tSpawnCountInitiatorThreshold) then
            spwn_cnt = get_spawn_count(site_idx)
            if (spwn_cnt >= init_spawn_thresh) &
                initiator = .true.
        end if

    end function TestInitiator

    subroutine rezero_iter_stats_each_iter(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data

        real(dp) :: prev_AvNoatHF(lenof_sign), AllInstNoatHF(lenof_sign)
        integer :: irdm, ind1, ind2, part_type

        NoInitDets = 0
        NoNonInitDets = 0
        NoInitWalk = 0.0_dp
        NoNonInitWalk = 0.0_dp
        InitRemoved = 0

        NoAborted = 0.0_dp
        NoRemoved = 0.0_dp
        NoatHF = 0.0_dp
        NoatDoubs = 0.0_dp

        iter_data%nborn = 0
        iter_data%ndied = 0
        iter_data%nannihil = 0
        iter_data%naborted = 0
        iter_data%nremoved = 0

        call InitHistMin()

        if (tFillingStochRDMonFly) then
            call MPISumAll(InstNoatHF, AllInstNoAtHF)
            InstNoAtHF = AllInstNoAtHF

            if (tFullHFAv) then
                Prev_AvNoatHF = AvNoatHF

                do ind1 = 1, lenof_sign
                    if (abs(IterRDM_HF(ind1)) > 1.0e-12_dp) then
                        AvNoatHF(ind1) = ( (real((Iter+PreviousCycles - IterRDM_HF(ind1)),dp) * Prev_AvNoatHF(ind1)) &
                                                + InstNoatHF(ind1) ) / real((Iter+PreviousCycles - IterRDM_HF(ind1)) + 1, dp)
                    end if
                end do
            else
                if (((Iter+PreviousCycles-IterRDMStart) .gt. 0) .and. &
                    & (mod(((Iter-1)+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) .eq. 0)) then 
                    ! The previous iteration was one where we added in diagonal
                    ! elements To keep things unbiased, we need to set up a new
                    ! averaging block now.
                    AvNoAtHF = InstNoAtHF
                    IterRDM_HF = Iter + PreviousCycles
                else
                    if (tPairedReplicas) then
                        do irdm = 1, nrdms

                            ! The indicies of the first and second replicas in this
                            ! particular pair, in the sign arrays.
                            ind1 = irdm*2-1
                            ind2 = irdm*2

                            if ((abs(InstNoAtHF(ind1)) < 1.0e-12_dp .and. abs(IterRDM_HF(ind1)) > 1.0e-12_dp) .or.&
                                (abs(InstNoAtHF(ind2)) < 1.0e-12_dp .and. abs(IterRDM_HF(ind2)) > 1.0e-12_dp)) then
                                ! At least one of the populations has just become
                                ! zero. Start a new averaging block.
                                IterRDM_HF(ind1) = Iter + PreviousCycles
                                IterRDM_HF(ind2) = Iter + PreviousCycles
                                AvNoatHF(ind1) = InstNoAtHF(ind1)
                                AvNoatHF(ind2) = InstNoAtHF(ind2)
                                if (abs(InstNoAtHF(ind1)) < 1.0e-12_dp) IterRDM_HF(ind1) = 0
                                if (abs(InstNoAtHF(ind2)) < 1.0e-12_dp) IterRDM_HF(ind2) = 0

                            else if ((abs(InstNoAtHF(ind1)) > 1.0e-12_dp .and. abs(IterRDM_HF(ind1)) < 1.0e-12_dp) .or. &
                                     (abs(InstNoAtHF(ind2)) > 1.0e-12_dp .and. abs(IterRDM_HF(ind2)) < 1.0e-12_dp)) then
                                ! At least one of the populations has just
                                ! become occupied. Start a new block here.
                                IterRDM_HF(ind1) = Iter + PreviousCycles
                                IterRDM_HF(ind2) = Iter + PreviousCycles
                                AvNoAtHF(ind1) = InstNoAtHF(ind1)
                                AvNoAtHF(ind2) = InstNoAtHF(ind2)
                                if (abs(InstNoAtHF(ind1)) < 1.0e-12_dp) IterRDM_HF(ind1) = 0
                                if (abs(InstNoAtHF(ind2)) < 1.0e-12_dp) IterRDM_HF(ind2) = 0
                            else
                                Prev_AvNoatHF = AvNoatHF

                                do part_type = 2*irdm-1, 2*irdm
                                    if (abs(IterRDM_HF(part_type)) > 1.0e-12_dp) then
                                         AvNoatHF(part_type) = ((real((Iter+PreviousCycles - IterRDM_HF(part_type)), dp) &
                                                                    * Prev_AvNoatHF(part_type)) + InstNoatHF(part_type) ) &
                                                                    / real((Iter+PreviousCycles - IterRDM_HF(part_type)) + 1, dp)
                                    end if
                                end do
                            end if
                        end do

                    else
                        do irdm = 1, nrdms

                            if (abs(InstNoAtHF(irdm)) < 1.0e-12_dp .and. abs(IterRDM_HF(irdm)) > 1.0e-12_dp) then
                                ! At least one of the populations has just become
                                ! zero. Start a new averaging block.
                                IterRDM_HF(irdm) = 0
                                AvNoatHF(irdm) = 0.0_dp

                            else if (abs(InstNoAtHF(irdm)) > 1.0e-12_dp .and. abs(IterRDM_HF(irdm)) < 1.0e-12_dp) then
                                ! At least one of the populations has just
                                ! become occupied. Start a new block here.
                                IterRDM_HF(irdm) = Iter + PreviousCycles
                                AvNoAtHF(irdm) = InstNoAtHF(irdm)
                            else
                                Prev_AvNoatHF = AvNoatHF

                                if (IterRDM_HF(irdm) /= 0) then
                                     AvNoatHF(irdm) = ((real((Iter+PreviousCycles - IterRDM_HF(irdm)), dp) &
                                                                * Prev_AvNoatHF(irdm)) + InstNoatHF(irdm) ) &
                                                                / real((Iter+PreviousCycles - IterRDM_HF(irdm)) + 1, dp)
                                end if
                            end if
                        end do

                    end if
                end if
            end if
        end if
        HFInd = 0

        min_trial_ind = 1
        min_conn_ind = 1

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

        integer, intent(in) :: TotWalkersNew

        ! SumWalkersCyc calculates the total number of walkers over an update
        ! cycle on each process.
#ifdef __CMPLX
        SumWalkersCyc = SumWalkersCyc + sum(TotParts)
#else
        SumWalkersCyc = SumWalkersCyc + TotParts
#endif

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
            TempSign(:)=0
        ENDIF

        HighPopInNeg(1) = TempSign(1)
        HighPopInNeg(2)=int(iProcIndex,int32)

        CALL MPIAllReduceDatatype(HighPopinNeg,1,MPI_MINLOC,MPI_2DOUBLE_PRECISION,HighPopoutNeg)

        IF(TotWalkersNew.gt.0) THEN
            call extract_sign(CurrentDets(:,HighPopPos),TempSign)
        ELSE
            TempSign(:)=0
        ENDIF

        HighPopInPos(1) = TempSign(1)
        HighPopInPos(2)=int(iProcIndex,int32)

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

        bAllowed = .true.

        ! Truncate space by excitation level
        if (tTruncSpace) then
            ! If parent walker is one below excitation cutoff, could be
            ! disallowed if double. If higher, then all excits could
            ! be disallowed. If HPHF, excit could be single or double,
            ! and IC not returned --> Always test.
            if (tHPHF .or. WalkExcitLevel >= ICILevel .or. &
                (WalkExcitLevel == (ICILevel-1) .and. IC == 2)) then
                ExcitLevel = FindBitExcitLevel (iLutHF, ilutnJ, ICILevel)
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
            abstr=''
            write(abstr,'(I12)') Iter
            abstr='TruncWavefunc-'//adjustl(abstr)
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
                call set_av_sgn(i, -get_av_sgn(i))
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
            write(iout,'(A)') 'Diagonalising initator subspace...'

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

    subroutine walker_death (iter_data, DetCurr, iLutCurr, Kii, RealwSign, &
                             DetPosition, walkExcitLevel)

        use rdm_data, only: one_rdms, rdms, two_rdm_spawn

        integer, intent(in) :: DetCurr(nel) 
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        integer(kind=n_int), intent(in) :: iLutCurr(0:niftot)
        real(dp), intent(in) :: Kii
        integer, intent(in) :: DetPosition
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp), dimension(lenof_sign) :: iDie
        real(dp), dimension(lenof_sign) :: CopySign
        integer, intent(in) :: walkExcitLevel
        integer :: i, irdm

        ! Do particles on determinant die? iDie can be both +ve (deaths), or
        ! -ve (births, if shift > 0)
        iDie = attempt_die (DetCurr, Kii, realwSign, WalkExcitLevel)

        IFDEBUG(FCIMCDebug,3) then 
            if (sum(abs(iDie)) > 1.0e-10_dp) write(iout,"(A,2f10.5)") "Death: ",iDie(:)
        endif

        ! Update death counter
        iter_data%ndied = iter_data%ndied + min(iDie, abs(RealwSign))
#ifdef __CMPLX
        NoDied(1) = NoDied(1) + sum(min(iDie, abs(RealwSign)))
#else
        NoDied = NoDied + min(iDie, abs(RealwSign))
#endif

        ! Count any antiparticles
        iter_data%nborn = iter_data%nborn + max(iDie - abs(RealwSign), 0.0_dp)
#ifdef __CMPLX
        NoBorn(1) = NoBorn(1) + sum(max(iDie - abs(RealwSign), 0.0_dp))
#else
        NoBorn = NoBorn + max(iDie - abs(RealwSign), 0.0_dp)
#endif

        ! Calculate new number of signed particles on the det.
        CopySign = RealwSign - (iDie * sign(1.0_dp, RealwSign))

        ! In the initiator approximation, abort any anti-particles.
        if (tTruncInitiator .and. any(abs(CopySign) > 1.0e-12_dp)) then
            do i = 1, lenof_sign
                if (CopySign(i) > 0.0_dp .neqv. RealwSign(i) > 0.0_dp) then
                    NoAborted = NoAborted + abs(CopySign(i))
                    iter_data%naborted(i) = iter_data%naborted(i) &
                                          + abs(CopySign(i))
                    if (test_flag(ilutCurr, flag_initiator(i))) &
                        NoAddedInitiators = NoAddedInitiators - 1
                    CopySign(i) = 0
                end if
            end do
        end if

        if (any(abs(CopySign) > 1.0e-12_dp)) then
            ! For the hashed walker main list, the particles don't move.
            ! Therefore just adjust the weight.
            call encode_sign (CurrentDets(:, DetPosition), CopySign)
        else
            ! All walkers died.
            if(tFillingStochRDMonFly) then
                do irdm = 1, nrdms
                    call det_removed_fill_diag_rdm(two_rdm_spawn, one_rdms(irdm), rdms(irdm), irdm, &
                                                    CurrentDets(:,DetPosition), DetPosition)
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

        ! Test - testsuite, RDM still work, both still work with Linscalealgo (all in debug)
        ! Null particle not kept if antiparticles aborted.
        ! When are the null particles removed?

    end subroutine walker_death

    subroutine check_start_rdm()

        ! This routine checks if we should start filling the RDMs - 
        ! and does so if we should. 

        use rdm_general, only: DeAlloc_Alloc_SpawnedParts
        use LoggingData, only: tReadRDMs
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

                call DeAlloc_Alloc_SpawnedParts()
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

        ! Update the reference used for a particular run to the one specified.
        ! Update the HPHF flipped arrays, and adjust the stored diagonal
        ! energies to account for the change if necessary.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: run
        character(*), parameter :: this_routine = 'update_run_reference'

        HElement_t(dp) :: h_tmp
        real(dp) :: old_hii
        integer :: i, det(nel)
        logical :: tSwapped

        iLutRef(:, run) = 0
        iLutRef(0:NIfDBO, run) = ilut(0:NIfDBO)
        call decode_bit_det (ProjEDet(:, run), iLutRef(:, run))
        write (iout, '(a,i3,a)', advance='no') 'Changing projected &
              &energy reference determinant for run', run, &
              ' on the next update cycle to: '
        call write_det (iout, ProjEDet(:, run), .true.)

        if(tHPHF) then
            if(.not.Allocated(RefDetFlip)) then
                allocate(RefDetFlip(NEl, inum_runs), &
                         ilutRef(0:NifTot, inum_runs))
                RefDetFlip = 0
                iLutRefFlip = 0
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
                write(iout,"(A)") "Now projecting onto open-shell &
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
            if (tHPHF) then
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

    end subroutine update_run_reference

    subroutine calc_inst_proje()

        ! Calculate an instantaneous value of the projected energy for the
        ! given walkers distributions

        integer :: ex_level, det(nel), j, run
        real(dp), dimension(max(lenof_sign,inum_runs)) :: RealAllHFCyc
        real(dp) :: sgn(lenof_sign)

        ! Reset the accumulators
        HFCyc = 0
        ENumCyc = 0

        ! Main loop
        do j = 1, int(TotWalkers, sizeof_int)

            ! n.b. non-contiguous list
            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

            ex_level = FindBitExcitLevel (iLutRef, CurrentDets(:,j))

            call decode_bit_det(det, CurrentDets(:,j))
            call SumEContrib(det, ex_level, sgn, CurrentDets(:,j), 0.0_dp, &
                             1.0_dp, tPairedReplicas, j)
        end do

        ! Accumulate values over all processors
        call MPISum(HFCyc, RealAllHFCyc)
        call MPISum(ENumCyc, AllENumCyc)

        do run = 1, inum_runs
            AllHFCyc(run) = ARR_RE_OR_CPLX(RealAllHFCyc, run)
        end do

        proje_iter = AllENumCyc / AllHFCyc + proje_ref_energy_offsets

        write(6,*) 'Calculated instantaneous projected energy', proje_iter

    end subroutine

end module
