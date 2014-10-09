#include "macros.h"

module fcimc_helper

    use constants
    use util_mod
    use systemData, only: nel, tHPHF, tNoBrillouin, G1, tUEG, &
                          tLatticeGens, nBasis
    use bit_reps, only: NIfTot, flag_is_initiator, test_flag, extract_flags, &
                        flag_parent_initiator, encode_bit_rep, NIfD, &
                        set_flag_general, flag_make_initiator, NIfDBO, &
                        extract_sign, set_flag, extract_first_iter, &
                        flag_trial, flag_connected, flag_deterministic
    use spatial_initiator, only: add_initiator_list, rm_initiator_list
    use DetBitOps, only: FindBitExcitLevel, FindSpatialBitExcitLevel, &
                         DetBitEQ, count_open_orbs
    use Determinants, only: get_helement
    use FciMCData
    use hist, only: test_add_hist_spin_dist_det, add_hist_spawn, &
                    add_hist_energies, HistMinInd, tHistSpawn
    use hphf_integrals, only: hphf_off_diag_helement
    use Logging, only: OrbOccs, tPrintOrbOcc, tPrintOrbOccInit, &
                       tHistSpinDist, tHistSpawn, tHistEnergies, &
                       RDMEnergyIter, tFullHFAv, &
                       nHistEquilSteps, tCalcFCIMCPsi, StartPrintOrbOcc, &
                       HistInitPopsIter, tHistInitPops
    use CalcData, only: NEquilSteps, tFCIMC, tSpawnSpatialInit, tTruncCAS, &
                        tRetestAddToInit, tAddToInitiator, InitiatorWalkNo, &
                        tTruncInitiator, tTruncNopen, trunc_nopen_max, &
                        tRealCoeffByExcitLevel, tSurvivalInitiatorThreshold, &
                        tSemiStochastic, tTrialWavefunction, nItersInitiator, &
                        InitiatorCutoffEnergy, InitiatorCutoffWalkNo, &
                        tLogComplexPops
    use DeterminantData, only: FDet
    use IntegralsData, only: tPartFreezeVirt, tPartFreezeCore, NElVirtFrozen, &
                             nPartFrozen, nVirtPartFrozen, nHolesFrozen
    use DetCalcData, only: FCIDetIndex, ICILevel
    use hash, only: DetermineDetNode
    use nElRDMMod, only: store_parent_with_spawned
    use Parallel_neci
    use FciMCLoggingMod, only: HistInitPopulations, WriteInitPops
    use csf_data, only: csf_orbital_mask
    use fcimc_initialisation, only: DeallocFCIMCMemPar, SetupParameters, &
                                    InitFCIMCCalcPar
    implicit none
    save

contains

    subroutine create_particle (nJ, iLutJ, child, parent_flags, part_type, &
                                ilutI, SignCurr, WalkerNo, RDMBiasFacCurr, &
                                WalkersToSpawn)

        ! Create a child in the spawned particles arrays. We spawn particles
        ! into a separate array, but non-contiguously. The processor that the
        ! newly-spawned particle is going to be sent to has to be determined,
        ! and then it will be put into the appropriate element determined by
        ! ValidSpawnedList

        ! 'type' of the particle - i.e. real/imag

        integer, intent(in) :: nJ(nel), parent_flags, part_type
        integer(n_int), intent(in) :: iLutJ(0:niftot)
        real(dp), intent(in) :: child(lenof_sign)
        integer(n_int), intent(in), optional :: ilutI(0:niftot)
        real(dp), intent(in), optional :: SignCurr(lenof_sign)
        integer, intent(in), optional :: WalkerNo
        real(dp), intent(in), optional :: RDMBiasFacCurr
        integer, intent(in), optional :: WalkersToSpawn
        integer :: proc, flags, j
        logical :: parent_init

        proc = DetermineDetNode(nel,nJ,0)    ! 0 -> nNodes-1)
        ! We need to include any flags set both from the parent and from the
        ! spawning steps. No we don't! - ghb
        ! This is highly yucky and needs cleaning up.
        ! Potentially, ilutJ can be given the flag of its parent in
        ! FindExcitBitDet routine. I don't think it should be. To make things
        ! confusing, this only happens for non-HPHF/CSF runs.
        ! Things are even more confusing given the fact that CCMC is using
        ! tihs routine. Who knows whether theyrequire the flag there or not...
        ! TODO: CLEAN THIS UP. Make it clear, and transparent, with one way
        !       to change the flag. Otherwise, this will trip many people up
        !       in the future.
        flags = ior(parent_flags, extract_flags(ilutJ))

        call encode_bit_rep(SpawnedParts(:, ValidSpawnedList(proc)), iLutJ, &
                            child, flags)

        if (tFillingStochRDMonFly) then
            ! We are spawning from ilutI to 
            ! SpawnedParts(:,ValidSpawnedList(proc)). We want to store the
            ! parent (D_i) with the spawned child (D_j) so that we can add in
            ! Ci.Cj to the RDM later.
            ! The parent is NIfDBO integers long, and stored in the second
            ! part of the SpawnedParts array from NIfTot+1 --> NIfTot+1+NIfDBO
            call store_parent_with_spawned (RDMBiasFacCurr, WalkerNo, &
                                            ilutI, WalkersToSpawn, ilutJ, &
                                            proc, part_type)
        end if

        if (tTruncInitiator) then
            IF(lenof_sign.eq.2) THEN
                ! With complex walkers, things are a little more tricky.
                ! We want to transfer the flag for all particles created (both
                ! real and imag) from the specific type of parent particle. This
                ! can mean real walker flags being transfered to imaginary
                ! children and vice versa.
                ! This is unneccesary for real walkers.
                ! Test the specific flag corresponding to the parent, of type
                ! 'part_type'
                parent_init = test_flag(SpawnedParts(:,ValidSpawnedList(proc)), &
                                        flag_parent_initiator(part_type))
                !Assign this flag to all spawned children
                do j=1,lenof_sign
                    if (child(j) /= 0) then
                        call set_flag (SpawnedParts(:,ValidSpawnedList(proc)), &
                                       flag_parent_initiator(j), parent_init)
                    endif
                enddo
            ENDIF
        end if

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
        
        ! Sum the number of created children to use in acceptance ratio.
        acceptances = acceptances + int(sum(abs(child)), kind(acceptances))
    end subroutine


    ! This routine sums in the energy contribution from a given walker and 
    ! updates stats such as mean excit level AJWT added optional argument 
    ! dProbFin which is a probability that whatever gave this contribution 
    ! was generated. It defaults to 1, and weights the contribution of this
    ! det (only in the projected energy) by dividing its contribution by 
    ! this number 
    subroutine SumEContrib (nI, ExcitLevel, RealWSign, ilut, HDiagCurr, &
                            dProbFin, ind)

        integer, intent(in) :: nI(nel), ExcitLevel
        real(dp), intent(in) :: RealwSign(lenof_sign)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: HDiagCurr, dProbFin
        integer, intent(in), optional :: ind

        integer :: i, k, bin, pos, ExcitLevel_local, ExcitLevelSpinCoup
        integer :: PartInd, OpenOrbs, spatial_ic, run
        integer(n_int) :: iLutSym(0:NIfTot)
        real(dp) :: ovp
        logical tSuccess
        integer :: iUEG1, iUEG2, ProjEBin
        HElement_t :: HOffDiag
        HElement_t :: HDoubDiag
        integer :: DoubEx(2,2),DoubEx2(2,2),kDoub(3) ! For histogramming UEG doubles
        integer :: ExMat(2,2), nopen
        integer :: doub_parity, doub_parity2, parity
        character(*), parameter :: this_routine = 'SumEContrib'

        HOffDiag = 0

        ! Add in the contributions to the numerator and denominator of the trial
        ! estimator, if it is being used.
        if (tTrialWavefunction .and. present(ind)) then
            if (tHashWalkerlist) then
                if (test_flag(ilut, flag_trial)) then
                    do run = 1, inum_runs
                        trial_denom(run) = trial_denom(run) &
                                    + current_trial_amps(ind)*RealwSign(run)
                    end do
                else if (test_flag(ilut, flag_connected)) then
                    do run = 1, inum_runs
                        trial_numerator(run) = trial_numerator(run) &
                                    + current_trial_amps(ind) * RealwSign(run)
                    end do
                end if
            else
                if (test_flag(ilut, flag_trial)) then
                    ! Take the next element in the occupied trial vector.
                    trial_ind = trial_ind + 1
                    do run = 1, inum_runs
                        trial_denom(run) = trial_denom(run) &
                                + occ_trial_amps(trial_ind) * RealwSign(run)
                    end do
                else if (test_flag(ilut, flag_connected)) then
                    ! Take the next element in the occupied connected vector.
                    con_ind = con_ind + 1
                    do run = 1, inum_runs
                        trial_numerator(run) = trial_numerator(run) &
                                    + occ_con_amps(con_ind) * RealwSign(run)
                    end do
                end if
            end if
        end if

        ! ExcitLevel indicates the excitation level between the det and
        ! *one* of the determinants in an HPHF/MomInv function. If needed,
        ! calculate the connection between it and the other one. If either
        ! is connected, then it has to be counted. Since the excitation
        ! level is the same to either det, we don't need to consider the
        ! spin-coupled det of both reference and current HPHFs.
        !
        ! For determinants, set ExcitLevel_local == ExcitLevel.
        ExcitLevel_local = ExcitLevel
        if (tSpinCoupProjE .and. (ExcitLevel /= 0)) then
            ExcitLevelSpinCoup = FindBitExcitLevel (iLutRefFlip, &
                                                    ilut, 2)
            if (ExcitLevelSpinCoup <= 2 .or. ExcitLevel <= 2) &
                ExcitLevel_local = 2
        endif

        ! Perform normal projection onto reference determinant
        if (ExcitLevel_local == 0) then

            if (iter > NEquilSteps) &
                SumNoatHF = SumNoatHF + RealwSign
            NoatHF = NoatHF + RealwSign
            ! Number at HF * sign over course of update cycle
            HFCyc = HFCyc + RealwSign

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
                HOffDiag = hphf_off_diag_helement (ProjEDet, nI, iLutRef,&
                                                   ilut)
            else
                HOffDiag = get_helement (ProjEDet, nI, ExcitLevel, &
                                         ilutRef, ilut)
            endif

        endif ! ExcitLevel_local == 1, 2, 3


        ! Sum in energy contribution
        do run=1, inum_runs
            if (iter > NEquilSteps) &
                SumENum(run) = SumENum(run) + (HOffDiag * ARR_RE_OR_CPLX(RealwSign,run)) &
                                  / dProbFin

        ENumCyc(run) = ENumCyc(run) + (HOffDiag * ARR_RE_OR_CPLX(RealwSign,run)) / dProbFin
        ENumCycAbs(run) = ENumCycAbs(run) + abs(HoffDiag * ARR_RE_OR_CPLX(RealwSign,run)) &
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
            if ((tPrintOrbOccInit .and. test_flag(ilut,flag_is_initiator(1)))&
                .or. .not. tPrintOrbOccInit) then
                forall (i = 1:nel) OrbOccs(iand(nI(i), csf_orbital_mask)) &
                        = OrbOccs(iand(nI(i), csf_orbital_mask)) &
                                   + (RealwSign(1) * RealwSign(1))
            endif
        endif
        
    end subroutine SumEContrib


    subroutine CalcParentFlag(j, VecSlot, parent_flags, diagH)

        ! In the CurrentDets array, the flag at NIfTot refers to whether that
        ! determinant *itself* is an initiator or not. We need to decide if 
        ! this willchange due to the determinant acquiring a certain 
        ! population, or its population dropping below the threshold.
        ! The CurrentDets(:,j) is the determinant we are currently spawning 
        ! from, so this determines the ParentInitiator flag which is passed to
        ! the SpawnedDets array and refers to whether or not the walkers 
        ! *parent* is an initiator or not.

        integer, intent(in) :: j, VecSlot
        integer, intent(out) :: parent_flags
        real(dp) :: CurrentSign(lenof_sign), init_thresh, low_init_thresh
        real(dp), intent(in) :: diagH
        integer :: part_type, nopen, first_iter
        logical :: tDetinCAS, parent_init
        character(*), parameter :: this_routine = 'CalcParentFlag'

        call extract_sign (CurrentDets(:,j), CurrentSign)

        init_thresh = InitiatorWalkNo
        low_init_thresh = InitiatorCutoffWalkNo

        tcurr_initiator = .false.
        do part_type=1,lenof_sign

            ! By default, the parent_flags are the flags of the parent...
            parent_init = test_flag (CurrentDets(:,j), &
                                     flag_parent_initiator(part_type))

            ! The default path through this section makes no changes, leaving
            ! the initiator status of each parent unchanged.  If 
            ! tAddToInitiator is set, then the state of the parent may change.
            if (tAddToInitiator) then

                if (.not. parent_init) then
                    ! Determinant wasn't previously initiator 
                    ! - want to test if it has now got a large enough 
                    !   population to become an initiator.
                    if ((diagH > InitiatorCutoffEnergy &
                         .and. abs(CurrentSign(part_type)) > low_init_thresh) &
                        .or. abs(CurrentSign(part_type)) > init_thresh) then
                        parent_init = .true.
                        NoAddedInitiators = NoAddedInitiators + 1
                        if (tSpawnSpatialInit) &
                            call add_initiator_list (CurrentDets(:,j))
                    endif
                elseif (tRetestAddToInit) then
                    ! The source determinant is already an initiator.            
                    ! If tRetestAddToInit is on, the determinants become 
                    ! non-initiators again if their population falls below 
                    ! n_add (this is on by default).
                    tDetInCas = .false.
                    if (tTruncCAS) &
                        tDetInCas = TestIfDetInCASBit (CurrentDets(0:NIfD,j))
                   
                    ! If det. in fixed initiator space, or is the HF det, or it
                    ! is in the deterministic space, then it must remain an initiator.
                    if (.not. tDetInCas .and. &
                        .not. (DetBitEQ(CurrentDets(:,j), iLutHF, NIfDBO)) &
                        .and. .not. test_flag(CurrentDets(:,j), flag_deterministic) &
                        .and. abs(CurrentSign(part_type)) <= init_thresh &
                        .and. diagH <= InitiatorCutoffEnergy &
                        .and. .not. test_flag(CurrentDets(:,j), &
                        flag_make_initiator(part_type))) then
                        ! Population has fallen too low. Initiator status 
                        ! removed.
                        parent_init = .false.
                        NoAddedInitiators = NoAddedInitiators - 1
                        if (tSpawnSpatialInit) &
                            call rm_initiator_list (CurrentDets(:,j))
                    endif
                endif

                ! If this site has survived for a long time, but otherwise
                ! would not be an initiator, then it is possible we ought
                ! to be considering it as well.
                if (.not. parent_init .and. tSurvivalInitiatorThreshold) then
                    first_iter = extract_first_iter(CurrentDets(:,j))
                    if ((iter - first_iter) > nItersInitiator) &
                        parent_init = .true.
                end if
#ifdef __DEBUG
                if (tSurvivalInitiatorThreshold) then
                    first_iter = extract_first_iter(CurrentDets(:,j))
                    ASSERT(first_iter <= iter .and. first_iter > 0)
                end if
#endif
            endif

            ! Update counters as required.
            if (parent_init) then
                NoInitDets = NoInitDets + 1
                NoInitWalk = NoInitWalk + abs(CurrentSign(part_type))
            else
                NoNonInitDets = NoNonInitDets + 1
                NoNonInitWalk = NoNonInitWalk + abs(CurrentSign(part_type))
            endif

            ! Update the parent flag as required.
            call set_flag (CurrentDets(:,j), flag_parent_initiator(part_type),&
                           parent_init)
            
            if(parent_init) then                           
                tcurr_initiator = .true.
            endif

        enddo

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
             call HistInitPopulations (CurrentSign(1), VecSlot)
        endif

    end subroutine CalcParentFlag


    subroutine rezero_iter_stats_each_iter (iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp) :: prev_AvNoatHF(lenof_sign), AllInstNoatHF(lenof_sign)
        integer :: j

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

        if(tFillingStochRDMonFly) then
            call MPISumAll(InstNoatHF, AllInstNoAtHF)
            InstNoAtHF=AllInstNoAtHF
            if (tFullHFAv) then
                Prev_AvNoatHF(1) = AvNoatHF(1)
                if (IterRDM_HF(1).ne.0) AvNoatHF(1) = ( (real((Iter+PreviousCycles - IterRDM_HF(1)),dp) * Prev_AvNoatHF(1)) &
                                            + InstNoatHF(1) ) / real((Iter+PreviousCycles - IterRDM_HF(1)) + 1,dp)
                if(inum_runs.eq.2) then
                    Prev_AvNoatHF(inum_runs) = AvNoatHF(inum_runs)
                   if(IterRDM_HF(inum_runs).ne.0) AvNoatHF(inum_runs) = &
                               & ( (real((Iter+PreviousCycles - IterRDM_HF(inum_runs)),dp) * &
                                 &       Prev_AvNoatHF(inum_runs)) + InstNoatHF(inum_runs) ) &
                                 &   / real((Iter+PreviousCycles - IterRDM_HF(inum_runs)) + 1,dp)
                endif
            else
                if(((Iter+PreviousCycles-IterRDMStart).gt.0) .and. &
                    & (mod(((Iter-1)+PreviousCycles - IterRDMStart + 1),RDMEnergyIter).eq.0)) then 
                ! The previous iteration was one where we added in diagonal elements
                ! To keep things unbiased, we need to set up a new averaging block now.
                    AvNoAtHF=InstNoAtHF
                    IterRDM_HF(1)=real(Iter+PreviousCycles,dp)
                    IterRDM_HF(inum_runs)=real(Iter+PreviousCycles,dp)
                else
                    if((InstNoatHF(1).eq.0.0).and.(InstNoAtHF(lenof_sign).eq.0.0) &
                        .and. (.not. tSemiStochastic)) then
                        !The HF determinant won't be in currentdets, so the CurrentH averages will have been wiped.
                        !NB - there will be a small issue here if the HF determinant isn't in the core space
                        IterRDM_HF(1) = 0.0_dp
                        AvNoatHF(1) = 0.0_dp
                        if(inum_runs.eq.2) then
                            IterRDM_HF(inum_runs) = 0.0_dp 
                            AvNoatHF(inum_runs) = 0.0_dp
                        endif
                    elseif(((InstNoAtHF(1).eq.0.0).and.(IterRDM_HF(1).ne.0)) .or. &
                       &  ((InstNoAtHF(inum_runs).eq.0.0).and.(IterRDM_HF(inum_runs).ne.0))) then
                        !At least one of the populations has just become zero
                        !Start a new averaging block
                        IterRDM_HF(1) = Iter+PreviousCycles  
                        AvNoatHF(1) = InstNoAtHF(1)
                        IterRDM_HF(inum_runs) = Iter+PreviousCycles  
                        AvNoatHF(inum_runs) = InstNoAtHF(inum_runs)
                        do j=1,inum_runs
                            if(InstNoAtHF(j).eq.0) then
                                IterRDM_HF(j)=0
                            endif
                        enddo
                    elseif(((InstNoAtHF(1).ne.0).and.(IterRDM_HF(1).eq.0)) .or. &
                           ((InstNoAtHF(inum_runs).ne.0).and.(IterRDM_HF(inum_runs).eq.0))) then
                            !At least one of the populations has just become occupied
                            !Start a new block here
                            IterRDM_HF(1)=real(Iter+PreviousCycles,dp)
                            IterRDM_HF(inum_runs)=real(Iter+PreviousCycles,dp)
                            AvNoAtHF(1)=InstNoAtHF(1)
                            AvNoAtHF(inum_runs)=InstNoAtHF(inum_runs)
                            do j=1,inum_runs
                                if(InstNoAtHF(j).eq.0) then
                                    IterRDM_HF(j)=0
                                endif
                            enddo
                    else
                        Prev_AvNoatHF(1) = AvNoatHF(1)
                        if (IterRDM_HF(1).ne.0) AvNoatHF(1) =((real((Iter+PreviousCycles - IterRDM_HF(1)),dp) &
                                                        * Prev_AvNoatHF(1)) + InstNoatHF(1) ) &
                                                        / real((Iter+PreviousCycles - IterRDM_HF(1)) + 1,dp)
                        if(inum_runs.eq.2) then
                            Prev_AvNoatHF(inum_runs) = AvNoatHF(inum_runs)
                           if(IterRDM_HF(inum_runs).ne.0) AvNoatHF(inum_runs)=&
                                                ((real((Iter+PreviousCycles-IterRDM_HF(inum_runs)),dp) * &
                                                Prev_AvNoatHF(inum_runs)) + InstNoatHF(inum_runs) ) &
                                                / real((Iter+PreviousCycles - IterRDM_HF(inum_runs)) + 1,dp)
                        endif
                    endif
                endif
            endif
        endif
        HFInd = 0

        trial_ind = 0
        con_ind = 0
        min_trial_ind = 1
        min_conn_ind = 1

    end subroutine

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
        integer :: proc, pos, i, k
        real(dp) :: sgn(lenof_sign)

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

        IF(iProcIndex.eq.HighPopOutNeg(2)) DetNeg(:)=CurrentDets(:,HighPopNeg)
        IF(iProcIndex.eq.HighPopOutPos(2)) DetPos(:)=CurrentDets(:,HighPopPos)

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

        integer :: NoInFrozenCore, MinVirt, ExcitLevel, nopen, i
        ! For UEG
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


    ! 
    ! This is a null routine for encoding spawned sites
    ! --> DOES NOTHING!!!
    subroutine null_encode_child (ilutI, ilutJ, ic, ex)
        implicit none
        integer(kind=n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: ic, ex(2,2)
        integer(kind=n_int), intent(inout) :: ilutj(0:niftot)

        ! Avoid compiler warnings
        integer :: iUnused
        integer(n_int) :: iUnused2
        iLutJ(0) = iLutJ(0); iUnused = IC; iUnused = ex(2,2)
        iUnused2 = iLutI(0)
    end subroutine
end module
