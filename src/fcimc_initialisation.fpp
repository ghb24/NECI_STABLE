#include "macros.h"
#:include "macros.fpph"

module fcimc_initialisation

    use SystemData, only: tUseBrillouin, iRanLuxLev, tSpn, tHPHFInts, tHPHF, &
                          tRotateOrbs, tROHF, tFindCINatOrbs, nOccBeta, tHUB, &
                          nOccAlpha, tUHF, tBrillouinsDefault, ECore, tUEG, &
                          tNoSingExcits, tOddS_HPHF, tSpn, tNoBrillouin, G1, &
                          tAssumeSizeExcitgen, tMolproMimic, tMolpro, tFixLz, &
                          tRef_Not_HF, LzTot, LMS, tKPntSym, tReal, nBasisMax, &
                          tRotatedOrbs, MolproID, nBasis, arr, brr, nel, &
                          tPickVirtUniform, tGen_4ind_reverse, &
                          tGenHelWeighted, tGen_4ind_weighted, tLatticeGens, &
                          tGUGA, tGen_nosym_guga, ref_stepvector, ref_b_vector_int, &
                          ref_occ_vector, ref_b_vector_real, tUEGNewGenerator, &
                          tGen_4ind_2, tReltvy, t_new_real_space_hubbard, &
                          t_lattice_model, t_tJ_model, t_heisenberg_model, &
                          t_k_space_hubbard, t_3_body_excits, breathingCont, &
                          momIndexTable, t_trans_corr_2body, t_non_hermitian, &
                          tgen_guga_crude, &
                          t_uniform_excits, t_mol_3_body, t_ueg_transcorr, t_ueg_3_body, tLatticeGens, &
                          irrepOrbOffset, nIrreps, &
                          tTrcorrExgen, nClosedOrbs, irrepOrbOffset, nIrreps, &
                          nOccOrbs, tNoSinglesPossible, t_pcpp_excitgen, &
                          t_pchb_excitgen, tGAS
    use tc_three_body_data, only: ptriples
    use SymExcitDataMod, only: tBuildOccVirtList, tBuildSpinSepLists
    use core_space_util, only: cs_replicas
    use dSFMT_interface, only: dSFMT_init

    use CalcData, only: G_VMC_Seed, MemoryFacPart, TauFactor, StepsSftImag, &
                        tCheckHighestPop, tSpatialOnlyHash, tStartCAS, tau, &
                        MaxWalkerBloom, InitialPart, tStartMP1, tReadPops, &
                        InitialPartVec, iReadWalkersRoot, SinglesBias, NMCYC, &
                        tTruncCAS, tTruncInitiator, DiagSft, tFCIMC, &
                        tTrialWavefunction, tSemiStochastic, OccCASOrbs, &
                        VirtCASOrbs, StepsSft, tStartSinglePart, InitWalkers, &
                        tShiftOnHFPop, tReadPopsRestart, tTruncNOpen, &
                        trunc_nopen_max, MemoryFacInit, MaxNoatHF, HFPopThresh, &
                        tAddToInitiator, InitiatorWalkNo, tRestartHighPop, &
                        tAllRealCoeff, tRealCoeffByExcitLevel, tTruncInitiator, &
                        RealCoeffExcitThresh, aliasStem, tPopsAlias, &
                        tDynamicCoreSpace, TargetGrowRate, &
                        TargetGrowRateWalk, InputTargetGrowRate, semistoch_shift_iter, &
                        InputTargetGrowRateWalk, tOrthogonaliseReplicas, &
                        use_spawn_hash_table, tReplicaSingleDetStart, RealSpawnCutoff, &
                        ss_space_in, trial_space_in, init_trial_in, trial_shift_iter, &
                        tContTimeFCIMC, tContTimeFull, tMultipleInitialRefs, &
                        initial_refs, trial_init_reorder, tStartTrialLater, tTrialInit, &
                        ntrial_ex_calc, tPairedReplicas, tMultiRefShift, tPreCond, &
                        tMultipleInitialStates, initial_states, t_hist_tau_search, &
                        t_direct_guga_ref, &
                        t_previous_hist_tau, t_fill_frequency_hists, t_back_spawn, &
                        t_trunc_nopen_diff, t_guga_back_spawn, tExpAdaptiveShift, &
                        t_back_spawn_option, t_back_spawn_flex_option, &
                        t_back_spawn_flex, back_spawn_delay, ScaleWalkers, tfixedN0, &
                        tReplicaEstimates, tDeathBeforeComms, pSinglesIn, pParallelIn, &
                        tSetInitFlagsBeforeDeath, tSetInitialRunRef, tEN2Init, &
                        tAutoAdaptiveShift, &
                        tInitializeCSF, S2Init, tWalkContgrow, tSkipRef, &
                        AAS_Cut, &
                        tInitiatorSpace, i_space_in, tLinearAdaptiveShift, &
                        tExpAdaptiveShift, tAS_TrialOffset, ShiftOffset, &
                        tSpinProject

    use adi_data, only: tReferenceChanged, tAdiActive, nExChecks, nExCheckFails, &
                        nRefUpdateInterval, SIUpdateInterval

    use Determinants, only: GetH0Element3, GetH0Element4, tDefineDet, &
                            get_helement, get_helement_det_only

    use hphf_integrals, only: hphf_diag_helement, hphf_spawn_sign, &
                              hphf_off_diag_helement_spawn

    use SymData, only: SymLabelList, SymLabelCounts, TwoCycleSymGens, &
                       SymClassSize, nSymLabels, sym_psi

    use DeterminantData, only: write_det, write_det_len, FDet

    use LoggingData, only: tTruncRODump, tCalcVariationalEnergy, tReadRDMs, &
                           tDiagAllSpaceEver, tFCIMCStats2, tCalcFCIMCPsi, &
                           tLogComplexPops, tHistExcitToFrom, tPopsFile, &
                           iWritePopsEvery, tRDMOnFly, tCoupleCycleOutput, StepsPrint, &
                           tDiagWalkerSubspace, tPrintOrbOcc, OrbOccs, &
                           tHistInitPops, OrbOccsTag, tHistEnergies, tMCOutput, &
                           HistInitPops, AllHistInitPops, OffDiagMax, &
                           OffDiagBinRange, iDiagSubspaceIter, &
                           AllHistInitPopsTag, HistInitPopsTag, tHDF5PopsRead, &
                           tTransitionRDMs, tLogEXLEVELStats, t_no_append_stats, &
                           t_spin_measurements, &
                           maxInitExLvlWrite, initsPerExLvl, AllInitsPerExLvl

    use DetCalcData, only: NMRKS, tagNMRKS, FCIDets, NKRY, NBLK, B2L, nCycle, &
                           ICILevel, det
    use IntegralsData, only: tPartFreezeCore, nHolesFrozen, tPartFreezeVirt, &
                             nVirtPartFrozen, nPartFrozen, nelVirtFrozen

    use bit_rep_data, only: NIfTot, NIfD, IlutBits, flag_initiator, &
                            flag_deterministic, extract_sign

    use bit_reps, only: encode_det, clear_all_flags, set_flag, encode_sign, &
                        decode_bit_det, nullify_ilut, encode_part_sign, &
                        extract_run_sign, tBuildSpinSepLists, nifguga, &
                        get_initiator_flag, &
                        get_initiator_flag_by_run
    use hist_data, only: tHistSpawn, HistMinInd, HistMinInd2, Histogram, &
                         BeforeNormHist, InstHist, iNoBins, AllInstHist, &
                         HistogramEnergy, AllHistogramEnergy, AllHistogram, &
                         BinRange
    use hist, only: init_hist_excit_tofrom, clean_hist_excit_tofrom
    use PopsfileMod, only: FindPopsfileVersion, initfcimc_pops, &
                           ReadFromPopsfilePar, ReadPopsHeadv3, &
                           ReadPopsHeadv4, open_pops_head, checkpopsparams
    use HPHFRandExcitMod, only: gen_hphf_excit, FindDetSpinSym
    use GenRandSymExcitNUMod, only: gen_rand_excit, init_excit_gen_store, &
                                    clean_excit_gen_store
    use replica_estimates, only: open_replica_est_file
    use procedure_pointers, only: generate_excitation, attempt_create, &
                                  get_spawn_helement, encode_child, &
                                  attempt_die, extract_bit_rep_avsign, &
                                  fill_rdm_diag_currdet, &
                                  new_child_stats, get_conn_helement, scaleFunction, &
                                  generate_two_body_excitation, shiftFactorFunction, gen_all_excits
    use symrandexcit3, only: gen_rand_excit3
    use symrandexcit_Ex_Mag, only: gen_rand_excit_Ex_Mag
    use excit_gens_int_weighted, only: gen_excit_hel_weighted, &
                                       gen_excit_4ind_weighted, &
                                       gen_excit_4ind_reverse
    use hash, only: FindWalkerHash, add_hash_table_entry, init_hash_table, &
                    hash_table_lookup
    use load_balance_calcnodes, only: DetermineDetNode, RandomOrbIndex
    use load_balance, only: tLoadBalanceBlocks, addNormContribution, get_diagonal_matel, &
                            AddNewHashDet
    use SymExcit3, only: CountExcitations3, GenExcitations3
    use SymExcit4, only: CountExcitations4, GenExcitations4
    use HPHFRandExcitMod, only: ReturnAlphaOpenDet
    use FciMCLoggingMOD, only: InitHistInitPops
    use SymExcitDataMod, only: SymLabelList2, OrbClassCount, SymLabelCounts2
    use rdm_general, only: init_rdms, dealloc_global_rdm_data, &
                           extract_bit_rep_avsign_no_rdm
    use rdm_filling, only: fill_rdm_diag_currdet_norm
    use DetBitOps, only: FindBitExcitLevel, CountBits, TestClosedShellDet, &
                         FindExcitBitDet, IsAllowedHPHF, DetBitEq, &
                         EncodeBitDet, DetBitLT
    use fcimc_pointed_fns, only: att_create_trunc_spawn_enc, &
                                 attempt_create_normal, &
                                 attempt_create_trunc_spawn, &
                                 new_child_stats_hist_hamil, &
                                 new_child_stats_normal, &
                                 null_encode_child, attempt_die_normal, attempt_die_precond, &
                                 powerScaleFunction, expScaleFunction, negScaleFunction, &
                                 expCOScaleFunction, &
                                 expShiftFactorFunction, constShiftFactorFunction, &
                                 linearShiftFactorFunction, autoShiftFactorFunction

    use initial_trial_states, only: calc_trial_states_lanczos, &
                                    set_trial_populations, set_trial_states, calc_trial_states_direct
    use global_det_data, only: global_determinant_data, set_det_diagH, &
                               clean_global_det_data, init_global_det_data, &
                               set_spawn_rate, store_decoding
    use semi_stoch_gen, only: init_semi_stochastic, end_semistoch, &
                              enumerate_sing_doub_kpnt
    use semi_stoch_procs, only: return_mp1_amp_and_mp2_energy
    use initiator_space_procs, only: init_initiator_space
    use kp_fciqmc_data_mod, only: tExcitedStateKP
    use sym_general_mod, only: ClassCountInd
    use trial_wf_gen, only: init_trial_wf, end_trial_wf
    use load_balance, only: clean_load_balance, init_load_balance
    use ueg_excit_gens, only: gen_ueg_excit
    use gndts_mod, only: gndts
    use excit_gen_5, only: gen_excit_4ind_weighted2
    use tc_three_body_excitgen, only: gen_excit_mol_tc, setup_mol_tc_excitgen
    use pcpp_excitgen, only: gen_rand_excit_pcpp, init_pcpp_excitgen, finalize_pcpp_excitgen

    use tau_search, only: init_tau_search, max_death_cpt

    use fcimc_helper, only: CalcParentFlag, update_run_reference

    use cont_time_rates, only: spawn_rate_full, oversample_factors, &
                               secondary_gen_store, ostag

    use soft_exit, only: tSoftExitFound

    use get_excit, only: make_double

    use excitation_types, only: NoExc_t
    use sltcnd_mod, only: sltcnd_excit
    use rdm_data, only: nrdms_transition_input, rdmCorrectionFactor, InstRDMCorrectionFactor, &
                        ThisRDMIter
    use rdm_data, only: nrdms_transition_input
    use Parallel_neci

    use FciMCData

    use util_mod

    use fortran_strings, only: str

    use sort_mod

    use sym_mod

    use HElem

    use constants

    use guga_tausearch, only: init_tau_search_guga_nosym, log_spawn_magnitude_guga_nosym, &
                              update_tau_guga_nosym, init_hist_tau_search_guga_nosym, &
                              update_hist_tau_guga_nosym
    use guga_data, only: bVectorRef_ilut, bVectorRef_nI, projE_replica
    use guga_bitRepOps, only: calcB_vector_nI, calcB_vector_ilut, convert_ilut_toNECI, &
                              convert_ilut_toGUGA, getDeltaB, write_det_guga, write_guga_list, &
                              calc_csf_info

    use guga_excitations, only: generate_excitation_guga, create_projE_list, &
                                actHamiltonian
    use guga_matrixElements, only: calcDiagMatEleGUGA_ilut, calcDiagMatEleGUGA_nI

    use real_time_data, only: t_real_time_fciqmc

    use real_time_procs, only: attempt_create_realtime

    use adi_references, only: setup_reference_space, clean_adi

    use tau_search_hist, only: init_hist_tau_search

    use double_occ_mod, only: init_spin_measurements

    use back_spawn, only: init_back_spawn, setup_virtual_mask

    use real_space_hubbard, only: init_real_space_hubbard, init_get_helement_hubbard

    use back_spawn_excit_gen, only: gen_excit_back_spawn, gen_excit_back_spawn_ueg, &
                                    gen_excit_back_spawn_hubbard, gen_excit_back_spawn_ueg_new
    use gasci, only: GAS_exc_gen, possible_GAS_exc_gen, operator(==), GAS_specification
    use gasci_disconnected, only: gen_GASCI_disconnected, init_disconnected_GAS, clearGAS
    use gasci_general, only: gen_GASCI_general, gen_all_excits_GAS => gen_all_excits
    use gasci_discarding, only: gen_GASCI_discarding, init_GASCI_discarding, finalize_GASCI_discarding

    use cepa_shifts, only: t_cepa_shift, init_cepa_shifts

    use tj_model, only: init_get_helement_tj, init_get_helement_heisenberg, &
                        init_get_helement_heisenberg_guga, init_get_helement_tj_guga

    use k_space_hubbard, only: init_get_helement_k_space_hub, init_k_space_hubbard, &
                               gen_excit_k_space_hub_transcorr, gen_excit_uniform_k_space_hub_transcorr

    use OneEInts, only: tmat2d

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard, gen_all_excits_r_space_hubbard

    use pchb_excitgen, only: gen_rand_excit_pchb, init_pchb_excitgen, finalize_pchb_excitgen

    use symexcit3, only: gen_all_excits_default => gen_all_excits
    implicit none

contains

    SUBROUTINE SetupParameters()

        INTEGER :: ierr, i, j, HFDetTest(NEl), Seed, alpha, beta, symalpha, symbeta, endsymstate
        INTEGER :: LargestOrb, nBits, HighEDet(NEl), orb
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        HElement_t(dp) :: TempHii
        real(dp) :: UpperTau, r
        CHARACTER(len=*), PARAMETER :: t_r = 'SetupParameters'
        CHARACTER(*), parameter :: this_routine = t_r
        CHARACTER(len=12) :: abstr
        character(len=40) :: filename, filename2
        LOGICAL :: tSuccess, tFoundOrbs(nBasis), FoundPair, tSwapped, tAlreadyOcc
        INTEGER :: HFLz, ChosenOrb, step, SymFinal, run
        integer(int64) :: SymHF
        integer(n_int), allocatable :: dummy_list(:, :)

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors
        write(iout, *)
        if (nProcessors > 1) then
            write(iout, *) "Performing Parallel FCIQMC...."
        else
            write(iout, *) "Performing single-core FCIQMC...."
        end if
        write(iout, *)

!Set timed routine names
        Walker_Time%timer_name = 'WalkerTime'
        Annihil_Time%timer_name = 'AnnihilTime'
        Sort_Time%timer_name = 'SortTime'
        Comms_Time%timer_name = 'CommsTime'
        ACF_Time%timer_name = 'ACFTime'
        AnnSpawned_time%timer_name = 'AnnSpawnedTime'
        AnnMain_time%timer_name = 'AnnMainTime'
        BinSearch_time%timer_name = 'BinSearchTime'
        Sync_Time%timer_name = 'SyncTime'
        SemiStoch_Comms_Time%timer_name = 'SemiStochCommsTime'
        SemiStoch_Multiply_Time%timer_name = 'SemiStochMultiplyTime'
        Trial_Search_Time%timer_name = 'TrialSearchTime'
        SemiStoch_Init_Time%timer_name = 'SemiStochInitTime'
        SemiStoch_Hamil_Time%timer_name = 'SemiStochHamilTime'
        SemiStoch_Davidson_Time%timer_name = 'SemiStochDavidsonTime'
        Trial_Init_Time%timer_name = 'TrialInitTime'
        InitSpace_Init_Time%timer_name = 'InitSpaceTime'
        kp_generate_time%timer_name = 'KPGenerateTime'
        Stats_Comms_Time%timer_name = 'StatsCommsTime'
        subspace_hamil_time%timer_name = 'SubspaceHamilTime'
        exact_subspace_h_time%timer_name = 'ExactSubspace_H_Time'
        subspace_spin_time%timer_name = 'SubspaceSpinTime'
        var_e_time%timer_name = 'VarEnergyTime'
        precond_e_time%timer_name = 'PreCondEnergyTime'
        proj_e_time%timer_name = 'ProjEnergyTime'
        rescale_time%timer_name = 'RescaleTime'
        death_time%timer_name = 'DeathTime'
        hash_test_time%timer_name = 'HashTestTime'
        hii_test_time%timer_name = 'HiiTestTime'
        init_flag_time%timer_name = 'InitFlagTime'

        ! Initialise allocated arrays with input data
        TargetGrowRate(:) = InputTargetGrowRate
        TargetGrowRateWalk(:) = InputTargetGrowRateWalk

        ! Initialize
        AllTotParts = 0.0_dp
        AllTotPartsOld = 0.0_dp
        AllTotPartsLastOutput = 0.0_dp

        IF (TDebug) THEN
!This will open a file called LOCALPOPS-"iprocindex" on unit number 11 on every node.
            abstr = 'LOCALPOPS-'//str(iProcIndex)
            open(11, FILE=abstr, STATUS='UNKNOWN')
        end if

        IF (iProcIndex == Root) THEN
            if (.not. tFCIMCStats2) then
                fcimcstats_unit = get_free_unit()
                if (tReadPops .and. .not. t_no_append_stats) then
                    ! Restart calculation.  Append to stats file (if it exists).
                    if (tMolpro .and. .not. tMolproMimic) then
                        filename = 'FCIQMCStats_'//adjustl(MolproID)
                        open(fcimcstats_unit, file=filename, status='unknown', position='append')
                    else
                        open(fcimcstats_unit, file='FCIMCStats', status='unknown', position='append')
                    end if
                else
                    call MoveFCIMCStatsFiles()          !This ensures that FCIMCStats files are not overwritten
                    if (tMolpro .and. .not. tMolproMimic) then
                        filename = 'FCIQMCStats_'//adjustl(MolproID)
                        open(fcimcstats_unit, file=filename, status='unknown')
                    else
                        open(fcimcstats_unit, file='FCIMCStats', status='unknown')
                    end if
                end if
            end if
#ifndef PROG_NUMRUNS_
            if (inum_runs == 2) then
                fcimcstats_unit2 = get_free_unit()
                if (tReadPops .and. .not. t_no_append_stats) then
                    ! Restart calculation.  Append to stats file (if it exists).
                    if (tMolpro .and. .not. tMolproMimic) then
                        filename2 = 'FCIQMCStats2_'//adjustl(MolproID)
                        open(fcimcstats_unit2, file=filename2, status='unknown', position='append')
                    else
                        open(fcimcstats_unit2, file='FCIMCStats2', status='unknown', position='append')
                    end if
                else
                    if (tMolpro .and. .not. tMolproMimic) then
                        filename2 = 'FCIQMCStats2_'//adjustl(MolproID)
                        open(fcimcstats_unit2, file=filename2, status='unknown')
                    else
                        open(fcimcstats_unit2, file='FCIMCStats2', status='unknown')
                    end if
                end if
            end if
#endif

            IF (tTruncInitiator .and. (.not. tFCIMCStats2)) THEN
                initiatorstats_unit = get_free_unit()
                if (tReadPops .and. .not. t_no_append_stats) then
! Restart calculation.  Append to stats file (if it exists)
                    open(initiatorstats_unit, file='INITIATORStats', status='unknown', form='formatted', position='append')
                else
                    open(initiatorstats_unit, file='INITIATORStats', status='unknown', form='formatted')
                end if
            end if
            IF (tLogComplexPops) THEN
                ComplexStats_unit = get_free_unit()
                open(ComplexStats_unit, file='COMPLEXStats', status='unknown')
            end if
            if (tLogEXLEVELStats) then
                EXLEVELStats_unit = get_free_unit()
                open (EXLEVELStats_unit, file='EXLEVELStats', &
                      status='unknown', position='append')
            end if
        end if

!Store information specifically for the HF determinant
        allocate(HFDet(NEl), stat=ierr)
        CALL LogMemAlloc('HFDet', NEl, 4, t_r, HFDetTag)
        allocate(iLutHF(0:NIfTot), stat=ierr)
        IF (ierr /= 0) CALL Stop_All(t_r, "Cannot allocate memory for iLutHF")

        do i = 1, NEl
            HFDet(i) = FDet(i)
        end do
        CALL EncodeBitDet(HFDet, iLutHF)

        if (tHPHF .and. TestClosedShellDet(iLutHF) .and. tOddS_HPHF .and. TwoCycleSymGens) then
            !This is not a compatible reference function.
            !Create single excitation of the correct symmetry
            !Use this as the reference.
            write(6, "(A)") "Converging to ODD S eigenstate"

            SymFinal = int((G1(HFDet(nel))%Sym%S), sizeof_int) + 1

            tAlreadyOcc = .false.
            do i = SymLabelCounts(1, SymFinal), SymLabelCounts(1, SymFinal) + SymLabelCounts(2, SymFinal) - 1
                if (G1(HFDet(nel))%Ms == -1) then
                    !Choose beta ones
                    orb = (2 * SymLabelList(i)) - 1
                else
                    orb = (2 * SymLabelList(i))
                end if
                tAlreadyOcc = .false.
                do j = 1, nel
                    if (orb == HFDet(j)) then
                        tAlreadyOcc = .true.
                        exit
                    end if
                end do
                if (.not. tAlreadyOcc) then
                    HFDet(nel) = orb
                    call sort(HFDet)
                    exit
                end if
            end do
            if (tAlreadyOcc) &
                call stop_all(t_r, "Cannot automatically detect open-shell determinant for reference to use with odd S")
            call EncodeBitDet(HFDet, iLutHF)
            if (TestClosedShellDet(iLutHF)) &
                call stop_all(t_r, "Fail to create open-shell determinant for reference to use with odd S")
            write(6, *) "Reference determinant changed to the open-shell:"
            call write_det(iout, HFDet, .true.)
        end if

        !iLutRef is the reference determinant for the projected energy.
        !Initially, it is chosen to be the same as the inputted reference determinant
        call setup_adi()
        allocate(iLutRef(0:NIfTot, inum_runs), stat=ierr)
        ilutRef = 0
        allocate(ProjEDet(NEl, inum_runs), stat=ierr)

        IF (ierr /= 0) CALL Stop_All(t_r, "Cannot allocate memory for iLutRef")

        ! The reference / projected energy determinants are the same as the
        ! HF determinant.
        call assign_reference_dets()

        allocate(iLutHF_True(0:NIfTot), stat=ierr)
        IF (ierr /= 0) CALL Stop_All(t_r, "Cannot allocate memory for iLutHF_True")
        allocate(HFDet_True(NEl), stat=ierr)
        IF (ierr /= 0) CALL Stop_All(t_r, "Cannot allocate memory for HFDet_True")

        !RDM uses HFDet_True in some parts but ilutRef in others
        !Sorry here we make them the same to avoid errors there.
        !Let's hope that unkonwn parts of the code do not depend on HFDet_True

        if (tRef_Not_HF) then
            do i = 1, NEl
                HFDet_True(i) = BRR(i)
            end do
            call sort(HFDet_True(1:NEl))
            CALL EncodeBitDet(HFDet_True, iLutHF_True)
        else
            iLutHF_True = iLutHF
            HFDet_True = HFDet
        end if

        ! for GUGA calculations also save the b vector of the reference det
        ! also initialize the persistently stored list of H|ref> to
        ! calculate the projected energy
        if (tGUGA) then
            allocate(bVectorRef_nI(nEl), stat=ierr)
            allocate(bVectorRef_ilut(nBasis / 2), stat=ierr)

            ! store more information of the reference determinant, like
            ! stepvector and stuff

            ASSERT(allocated(ref_stepvector))

            call calc_csf_info(ilutRef, ref_stepvector, ref_b_vector_int, &
                               ref_occ_vector)

            ref_b_vector_real = real(ref_b_vector_int, dp)

            ! for multiple runs i have to initialize all the necessary
            ! projected energy lists
            ! have to first allocate the type proje_replica
            allocate(projE_replica(inum_runs), stat=ierr)

            ! only initialize that if we use the old way to calc the
            ! reference energy!
            ! for testing purposes initiaize both
            if (.not. t_direct_guga_ref) then
                do run = 1, inum_runs
                    call create_projE_list(run)
                end do
            end if
        end if

        if (tHPHF) then
            allocate(RefDetFlip(NEl, inum_runs), &
                      ilutRefFlip(0:NifTot, inum_runs))
            do run = 1, inum_runs
                if (.not. TestClosedShellDet(ilutRef(:, run))) then

                    ! If the reference determinant corresponds to an open shell
                    ! HPHF, then we need to specify the paired determinant and
                    ! mark that this needs to be considered in calculating
                    ! the projected energy.

                    tSpinCoupProjE(run) = .true.
                    call ReturnAlphaOpenDet(ProjEDet(:, run), &
                                            RefDetFlip(:, run), &
                                            ilutRef(:, run), &
                                            ilutRefFlip(:, run), &
                                            .true., .true., tSwapped)
                    if (tSwapped) &
                        write(iout, *) 'HPHF used, and open shell determinant &
                                      &for run ', run, ' spin-flippd for &
                                      &consistency.'
                    write(iout, *) "Two *different* determinants contained in &
                                  &initial HPHF for run ", run
                    write(iout, *) "Projected energy will be calculated as &
                                  &projection onto both of these."

                else
                    tSpinCoupProjE(run) = .false.
                end if
            end do

        else
            tSpinCoupProjE(:) = .false.
        end if

!Init hash shifting data
        hash_iter = 0

        IF (tFixLz) THEN
            CALL GetLz(HFDet, NEl, HFLz)
            write(iout, "(A,I5)") "Ml value of reference determinant is: ", HFLz
            IF (HFLz /= LzTot) THEN
                CALL Stop_All("SetupParameters", "Chosen reference determinant does not have the " &
                & //"same Lz value as indicated in the input.")
            end if
        end if

!Do a whole lot of tests to see if we can use Brillouins theorem or not.
        IF (tBrillouinsDefault) CALL CheckforBrillouins()

!test the encoding of the HFdet to bit representation.
        ! Test that the bit operations are working correctly...
        ! TODO: Move this to using the extract_bit_det routines to test those
        !       too...
        call decode_bit_det(HFDetTest, iLutHF)
        do i = 1, NEl
            IF (HFDetTest(i) /= HFDet(i)) THEN
                write(iout, *) "HFDet: ", HFDet(:)
                write(iout, *) "HFDetTest: ", HFDetTest(:)
                CALL Stop_All(t_r, "HF Determinant incorrectly decoded.")
            end if
        end do
        CALL LargestBitSet(iLutHF, NIfD, LargestOrb)
        IF (LargestOrb /= hfdet(nel)) then
            write(6, *) 'ilut HF', ilutHF
            write(6, *) 'largest orb', largestorb
            write(6, *) 'HFDet', hfdet
            CALL Stop_All(t_r, "LargestBitSet FAIL")
        end if
        nBits = CountBits(iLutHF, NIfD, NEl)
        IF (nBits /= NEl) THEN
            CALL Stop_All(t_r, "CountBits FAIL")
        end if

        allocate(HighestPopDet(0:NIfTot, inum_runs), stat=ierr)
        IF (ierr /= 0) CALL Stop_All(t_r, "Cannot allocate memory for HighestPopDet")
        HighestPopDet(:, :) = 0
        iHighestPop = 0

!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess = .true.
        tFoundOrbs(:) = .false.

        IF ((.not. tHub) .and. (.not. tUEG) .and. TwoCycleSymGens) THEN
            do i = 1, nSymLabels
!                write(iout,*) "NSymLabels: ",NSymLabels,i-1
                EndSymState = SymLabelCounts(1, i) + SymLabelCounts(2, i) - 1
!                write(iout,*) "Number of states: ",SymLabelCounts(2,i)
                do j = SymLabelCounts(1, i), EndSymState

                    Beta = (2 * SymLabelList(j)) - 1
                    Alpha = (2 * SymLabelList(j))
                    SymAlpha = INT((G1(Alpha)%Sym%S), sizeof_int)
                    SymBeta = INT((G1(Beta)%Sym%S), sizeof_int)
!                    write(iout,*) "***",Alpha,Beta

                    IF (.not. tFoundOrbs(Beta)) THEN
                        tFoundOrbs(Beta) = .true.
                    ELSE
                        CALL Stop_All("SetupParameters", "Orbital specified twice")
                    end if
                    IF (.not. tFoundOrbs(Alpha)) THEN
                        tFoundOrbs(Alpha) = .true.
                    ELSE
                        CALL Stop_All("SetupParameters", "Orbital specified twice")
                    end if

                    IF (G1(Beta)%Ms /= -1) THEN
                        tSuccess = .false.
                    else if (G1(Alpha)%Ms /= 1) THEN
                        tSuccess = .false.
                    else if ((SymAlpha /= (i - 1)) .or. (SymBeta /= (i - 1))) THEN
                        tSuccess = .false.
                    end if
                end do
            end do
            do i = 1, nBasis
                IF (.not. tFoundOrbs(i)) THEN
                    write(iout, *) "Orbital: ", i, " not found."
                    CALL Stop_All("SetupParameters", "Orbital not found")
                end if
            end do
        end if
        IF (.not. tSuccess) THEN
            write(iout, *) "************************************************"
            write(iout, *) "**                 WARNING!!!                 **"
            write(iout, *) "************************************************"
            write(iout, *) "Symmetry information of orbitals not the same in alpha and beta pairs."
            write(iout, *) "Symmetry now set up in terms of spin orbitals"
            write(iout, *) "I strongly suggest you check that the reference energy is correct."
        end if
        ! From now on, the orbitals are contained in symlabellist2 and
        ! symlabelcounts2 rather than the original arrays. These are stored
        ! using spin orbitals. Assume that if we want to use the non-uniform
        ! random excitation generator, we also want to use the NoSpinSym full
        ! excitation generators if they are needed.

        CALL GetSym(HFDet, NEl, G1, NBasisMax, HFSym)

        Sym_Psi = INT(HFSym%Sym%S, sizeof_int)  !Store the symmetry of the wavefunction for later
        write(iout, "(A,I10)") "Symmetry of reference determinant is: ", INT(HFSym%Sym%S, sizeof_int)

        if (TwoCycleSymGens) then
            SymHF = 0
            do i = 1, NEl
                SymHF = IEOR(SymHF, G1(HFDet(i))%Sym%S)
            end do
            write(iout, "(A,I10)") "Symmetry of reference determinant from spin orbital symmetry info is: ", SymHF
            if (SymHF /= HFSym%Sym%S) then
                call stop_all(t_r, "Inconsistency in the symmetry arrays.")
            end if
        end if

!If using a CAS space truncation, write out this CAS space
        IF (tTruncCAS) THEN
            write(iout, *) "Truncated CAS space detected. Writing out CAS space..."
            write(iout, '(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (", &
                (OccCASOrbs + VirtCASOrbs) / 2, ",", OccCASOrbs, ")"
            do I = NEl - OccCASorbs + 1, NEl
                write(iout, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(iout, G1(BRR(I))%SYM, .FALSE.)
                write(iout, '(I4)', advance='no') G1(BRR(I))%Ml
                write(iout, '(2F19.9)') ARR(I, 1), ARR(BRR(I), 2)
            end do
            write(iout, '(A)') " ================================================================================================="
            do I = NEl + 1, NEl + VirtCASOrbs
                write(iout, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(iout, G1(BRR(I))%SYM, .FALSE.)
                write(iout, '(I4)', advance='no') G1(BRR(I))%Ml
                write(iout, '(2F19.9)') ARR(I, 1), ARR(BRR(I), 2)
            end do
        else if (tTruncInitiator) THEN
            write(iout, '(A)') "*********** INITIATOR METHOD IN USE ***********"
            write(iout, '(A)') "Starting with only the reference determinant in the fixed initiator space."
        end if

        ! Setup excitation generator for the HF determinant. If we are using
        ! assumed sized excitgens, this will also be assumed size.
        IF (tUEG .or. tHub .or. tNoSingExcits) THEN
            exflag = 2
        ELSE
            exflag = 3
        end if
        IF (.not. tKPntSym) THEN
!Count all possible excitations - put into HFConn
!TODO: Get CountExcitations3 working with tKPntSym
            CALL CountExcitations3(HFDet, exflag, nSingles, nDoubles)
        ELSE
            if (t_k_space_hubbard) then
                ! use my gen_all_excits_k_space_hubbard routine from the
                ! unit tests.. but i might have to set up some other stuff
                ! for this to work and also make sure this works with my
                ! new symmetry implementation
                if (.not. t_trans_corr_2body) then
                    call gen_all_excits_k_space_hubbard(HFDet, nDoubles, dummy_list)
                end if
                nSingles = 0
            else
                ! Use Alex's old excitation generators to enumerate all excitations.
                call enumerate_sing_doub_kpnt(exflag, .false., nSingles, nDoubles, .false.)
            end if
        end if
        HFConn = nSingles + nDoubles

        if (tAutoAdaptiveShift .and. AAS_Cut < 0.0) then
            !The user did not specify the value, use this as a default
            if (HFConn > 0) then
                AAS_Cut = 1.0_dp / HFConn
            else
                ! if the HF is disconnected (can happen in rare corner cases), set it to 0
                AAS_Cut = 0.0_dp
            end if
        end if

        ! Initialise random number seed - since the seeds need to be different
        ! on different processors, subract processor rank from random number
        if (.not. tRestart) then
            Seed = abs(G_VMC_Seed - iProcIndex)
            write(iout, "(A,I12)") "Value for seed is: ", Seed
            !Initialise...
            CALL dSFMT_init(Seed)
            if (tMolpro) then
                if ((NMCyc == -1) .and. (.not. tTimeExit)) then
                    !No iteration number, or TIME option has been specified.
                    call warning_neci(t_r, &
                                 "No iteration number specified. Only running for 100 iterations initially. Change with ITERATIONS option.")
                    NMCyc = 100   !Only run for 100 iterations.
                else if (tTimeExit .and. (NMCyc == -1)) then
                    write(iout, "(A,F10.3,A)") "Running FCIQMC for ", MaxTimeExit / 60.0_dp, " minutes."
                else if (tTimeExit .and. (NMCyc /= -1)) then
                    write(iout, "(A,F10.3,A,I15,A)") "Running FCIQMC for ", MaxTimeExit / 60.0_dp, " minutes OR ", NMCyc, " iterations."
                else if ((.not. tTimeExit) .and. (NMCyc > 0)) then
                    write(iout, "(A,I15,A)") "Running FCIQMC for ", NMCyc, " iterations."
                else
                    call stop_all(t_r, "Iteration number/Time unknown for simulation - contact ghb")
                end if
            end if
        end if

        ! Option tRandomiseHashOrbs has now been removed.
        ! Its behaviour is now considered default
        ! --> Create a random mapping for the orbitals
        allocate(RandomOrbIndex(nBasis), stat=ierr)
        IF (ierr /= 0) THEN
            CALL Stop_All(t_r, "Error in allocating RandomOrbIndex")
        end if
        RandomOrbIndex(:) = 0

        ! We want another independent randomizing array for the hash table, so
        ! we do not introduce correlations between the two
        allocate(RandomHash2(nBasis), stat=ierr)
        if (ierr /= 0) &
            call stop_all(t_r, "Error in allocating RandomHash2")
        RandomHash2(:) = 0

        IF (iProcIndex == root) THEN
            do i = 1, nBasis
                ! If we want to hash only by spatial orbitals, then the
                ! spin paired orbitals must be set equal
                if (tSpatialOnlyHash) then
                    if (.not. btest(i, 0)) then
                        RandomOrbIndex(i) = RandomOrbIndex(i - 1)
                        cycle
                    end if
                end if

                ! Ensure that we don't set two values to be equal accidentally
                FoundPair = .false.
                do while (.not. FoundPair)
                    r = genrand_real2_dSFMT()
                    ChosenOrb = INT(nBasis * r * 1000) + 1

                    ! Check all values which have already been set.
                    do j = 1, nBasis
                        IF (RandomOrbIndex(j) == ChosenOrb) EXIT
                    end do

                    ! If not already used, then we can move on
                    if (j == nBasis + 1) FoundPair = .true.
                    RandomOrbIndex(i) = ChosenOrb
                end do
            end do

            !Do again for RandomHash2
            do i = 1, nBasis
                ! If we want to hash only by spatial orbitals, then the
                ! spin paired orbitals must be set equal
                if (tSpatialOnlyHash) then
                    if (.not. btest(i, 0)) then
                        RandomHash2(i) = RandomHash2(i - 1)
                        cycle
                    end if
                end if
                ! Ensure that we don't set two values to be equal accidentally
                FoundPair = .false.
                do while (.not. FoundPair)
                    r = genrand_real2_dSFMT()
                    ChosenOrb = INT(nBasis * r * 1000) + 1

                    ! Check all values which have already been set.
                    do j = 1, nBasis
                        IF (RandomHash2(j) == ChosenOrb) EXIT
                    end do

                    ! If not already used, then we can move on
                    if (j == nBasis + 1) FoundPair = .true.
                    RandomHash2(i) = ChosenOrb
                end do
            end do

!            write(iout,*) "Random Orbital Indexing for hash:"
!            write(iout,*) RandomHash2(:)
            if (tSpatialOnlyHash) then
                step = 2
            else
                step = 1
            end if
            do i = 1, nBasis
                IF ((RandomOrbIndex(i) == 0) .or. (RandomOrbIndex(i) > nBasis * 1000)) THEN
                    CALL Stop_All(t_r, "Random Hash incorrectly calculated")
                end if
                IF ((RandomHash2(i) == 0) .or. (RandomHash2(i) > nBasis * 1000)) THEN
                    CALL Stop_All(t_r, "Random Hash 2 incorrectly calculated")
                end if
                do j = i + step, nBasis, step
                    IF (RandomOrbIndex(i) == RandomOrbIndex(j)) THEN
                        CALL Stop_All(t_r, "Random Hash incorrectly calculated")
                    end if
                    IF (RandomHash2(i) == RandomHash2(j)) THEN
                        CALL Stop_All(t_r, "Random Hash 2 incorrectly calculated")
                    end if
                end do
            end do
        end if

        !Now broadcast to all processors
        CALL MPIBCast(RandomOrbIndex, nBasis)
        call MPIBCast(RandomHash2, nBasis)

        t_initialized_roi = .true.

        call init_load_balance()

        IF (tHPHF) THEN
            !IF(tLatticeGens) CALL Stop_All("SetupParameters","Cannot use HPHF with model systems currently.")
            IF (tROHF .or. (LMS /= 0)) CALL Stop_All("SetupParameters", "Cannot use HPHF with high-spin systems.")
            tHPHFInts = .true.
        end if

        if (t_lattice_model) then
            if (t_tJ_model) then
                if (tGUGA) then
                    call init_get_helement_tj_guga()
                else
                    call init_get_helement_tj()
                end if
            else if (t_heisenberg_model) then
                if (tGUGA) then
                    call init_get_helement_heisenberg_guga
                else
                    call init_get_helement_heisenberg()
                end if
            else if (t_new_real_space_hubbard) then
                call init_get_helement_hubbard()
            else if (t_k_space_hubbard) then
                call init_get_helement_k_space_hub()
            end if
        end if

!Calculate Hii (unless suppressed)
        if (tZeroRef) then
            TempHii = 0.0_dp
        else IF (tHPHF) THEN
            TempHii = hphf_diag_helement(HFDet, iLutHF)
        ELSE
            TempHii = get_helement(HFDet, HFDet, 0)
        end if
        Hii = REAL(TempHii, dp)
        write(iout, "(A,F20.10)") "Reference Energy set to: ", Hii
        if (tUEG) then
            !We require calculation of the sum of fock eigenvalues,
            !without knowing them - calculate from the full 1e matrix elements
            !of full hamiltonian removing two electron terms.
            TempHii = GetH0Element4(HFDet, HFDet)
        else
            !Know fock eigenvalues, so just use these.
            TempHii = GetH0Element3(hfdet)
        end if
        Fii = REAL(TempHii, dp)

!Find the highest energy determinant...
        IF (LMS == 0) THEN
            do i = 1, NEl
                HighEDet(i) = Brr(nBasis - (i - 1))
            end do
            call sort(HighEDet)
            IF (tHPHF) THEN
                call EncodeBitDet(HighEDet, iLutTemp)
                TempHii = hphf_diag_helement(HighEDet, iLutTemp)
            ELSE
                TempHii = get_helement(HighEDet, HighEDet, 0)
            end if
            if (abs(TempHii - Hii) > EPS) then
                UpperTau = 1.0_dp / REAL(TempHii - Hii, dp)
            else
                UpperTau = 0.0_dp
            end if
            write(iout, "(A,G25.15)") "Highest energy determinant is (approximately): ", REAL(TempHii, dp)
            write(iout, "(a,g25.15)") "Corresponding to a correlation energy of: ", real(temphii - hii, dp)
            write(iout, "(A,F25.15)") "This means tau should be no more than about ", UpperTau
            write(iout, *) "Highest energy determinant is: ", HighEDet(:)
        else
            UpperTau = 0.0_dp
        end if

        ! Initialise DiagSft according to the input parameters. If we have
        ! multiple projected-energy references, then the shift on each of the
        ! runs should be adjusted so that it is still relative to the first
        ! replica, but is offset by the replica's reference's diagonal energy.
        DiagSft = InputDiagSft
        proje_ref_energy_offsets = 0.0_dp
        if (tOrthogonaliseReplicas) then
            do run = 1, inum_runs
                if (tHPHF) then
                    TempHii = hphf_diag_helement(ProjEDet(:, run), ilutRef(:, run))
                else
                    TempHii = get_helement(ProjEDet(:, run), ProjEDet(:, run), 0)
                end if
                proje_ref_energy_offsets(run) = real(TempHii, dp) - Hii

                if (tMultiRefShift) DiagSft(run) = proje_ref_energy_offsets(run)
            end do
        end if

        IF (tHub) THEN
            IF (tReal) THEN
!We also know that in real-space hubbard calculations, there are only single excitations.
                exFlag = 1
            ELSE
!We are doing a momentum space hubbard calculation - set exFlag to 2 since only doubles are connected for momentum conservation.
                exFlag = 2
            end if
        end if

        IF (LMS /= 0) THEN
            IF (tNoBrillouin .or. (tHub .and. tReal) .or. tRotatedOrbs) THEN
                write(iout, *) "No brillouin theorem assumed. Single excitations also used to calculate projected energy."
            else if (tUHF) THEN
                write(iout, *) "High spin calculation - but single excitations will *NOT* be used to calculate energy as "&
                & //"this is an unrestricted calculation."
            ELSE
                CALL Stop_All("SetupParameters", "High-spin, restricted calculation detected, but single excitations are "&
                & //"not being used to calculate the energy.  &
                  & Either use the UHF keyword, or turn off brillouins theorem using NOBRILLOUINS, ROHF or ROTATEDORBS.")
            end if
!            tRotatedOrbs=.true.
!        else if(LMS.ne.0) THEN
!            CALL Stop_All(t_r,"Ms not equal to zero, but tSpn is false. Error here")
        end if

!Initialise variables for calculation on each node
        iter = 0          !This is set so that calls to CalcParentFlag in the initialisation are ok with the logging.
        iPopsTimers = 1   !Number of timed popsfiles written out
        iBlockingIter = 0
        IterTime = 0.0
        ProjectionE(:) = 0.0_dp
        AvSign = 0.0_dp
        AvSignHFD = 0.0_dp
        SumENum(:) = 0.0_dp
        InitsENumCyc(:) = 0.0_dp
        SumNoatHF(:) = 0.0_dp
        NoatHF(:) = 0.0_dp
        InstNoatHF(:) = 0.0_dp
        Annihilated(:) = 0.0_dp
        Acceptances(:) = 0.0_dp
        PreviousCycles = 0
        NoBorn = 0.0_dp
        SpawnFromSing = 0
        max_cyc_spawn = 0.0_dp
        ! in case tau-search is off
        max_death_cpt = 0.0_dp
        NoDied = 0
        HFCyc = 0.0_dp
        HFOut = 0.0_dp
        ENumCyc = 0.0_dp
        ENumOut = 0.0_dp
        ENUmCycAbs = 0.0_dp
        VaryShiftCycles = 0
        AvDiagSft(:) = 0.0_dp
        SumDiagSft(:) = 0.0_dp
        SumWalkersCyc(:) = 0.0_dp
        SumWalkersOut(:) = 0.0_dp
!        SumDiagSftAbort=0.0_dp
!        AvDiagSftAbort=0.0_dp
        NoAborted(:) = 0.0_dp
        NoRemoved(:) = 0.0_dp
        NoAddedInitiators = 0
        NoInitDets = 0
        NoNonInitDets = 0
        NoInitWalk(:) = 0.0_dp
        NoNonInitWalk(:) = 0.0_dp
        NoExtraInitDoubs = 0
        InitRemoved = 0
        TotImagTime = 0.0_dp
        DiagSftRe = 0.0_dp
        DiagSftIm = 0.0_dp
        sum_proje_denominator = 0
        cyc_proje_denominator = 0
!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum(:) = 0.0_dp
        AllInitsENumCyc(:) = 0.0_dp
        AllNoatHF(:) = 0.0_dp
        AllNoatDoubs(:) = 0.0_dp
        NoAtDoubs = 0.0_dp
        if (tLogEXLEVELStats) AllEXLEVEL_WNorm(:, :, :) = 0.0_dp
        AllSumNoatHF(:) = 0.0_dp
        AllGrowRate(:) = 0.0_dp
        AllGrowRateAbort(:) = 0
!        AllMeanExcitLevel=0.0_dp
        AllSumWalkersCyc(:) = 0
        AllSumWalkersOut(:) = 0.0_dp
        AllAvSign = 0.0_dp
        AllAvSignHFD = 0.0_dp
        AllNoBorn(:) = 0
        AllSpawnFromSing(:) = 0
        AllNoDied(:) = 0
        AllAnnihilated(:) = 0
        AllENumCyc(:) = 0.0_dp
        AllENumOut(:) = 0.0_dp
        AllENumCycAbs = 0.0_dp
        AllHFCyc(:) = 0.0_dp
        all_cyc_proje_denominator = 1.0_dp
        AllHFOut(:) = 0.0_dp
!        AllDetsNorm=0.0_dp
        AllNoAborted = 0
        AllNoRemoved = 0
        AllNoAddedInitiators = 0
        AllNoInitDets = 0
        AllNoNonInitDets = 0
        AllNoInitWalk = 0.0_dp
        AllNoNonInitWalk = 0.0_dp
        AllNoExtraInitDoubs = 0
        AllInitRemoved = 0
        all_n_core_non_init = 0
        n_core_non_init = 0
        proje_iter = 0
        inits_proje_iter = 0.0_dp
        AccRat = 0
        HFShift = 0
        bloom_count = 0
        InstShift = 0
        AbsProjE = 0
        norm_semistoch = 0
        norm_psi = 0
        bloom_sizes = 0
        proje_iter_tot = 0.0_dp
        ! initialize as one (kind of makes sense for a norm)
        norm_psi_squared = 1.0_dp
        all_norm_psi_squared = 1.0_dp
        norm_semistoch_squared = 1.0_dp
        tSoftExitFound = .false.
        tReferenceChanged = .false.

        ! Set the flag to indicate that no shift adjustment has been made
        tfirst_cycle = .true.

        ! Initialise the fciqmc counters
        iter_data_fciqmc%update_growth = 0.0_dp
        iter_data_fciqmc%update_iters = 0

        nExChecks = 0
        nExCheckFails = 0

        ! 0-initialize truncated weight
        truncatedWeight = 0.0_dp
        AllTruncatedWeight = 0.0_dp

        ! RDMs are taken as they are until we have some data on the f-function
        ! of the adaptive shift
        rdmCorrectionFactor = 0.0_dp
        InstRDMCorrectionFactor = 0.0_dp
        ThisRDMIter = 0.0_dp
        ! initialize excitation number trackers
        nInvalidExcits = 0
        nValidExcits = 0
        allNInvalidExcits = 0
        allNValidExcits = 0

        if (tEScaleWalkers) then
            if (abs(RealSpawnCutoff - sFBeta) > eps) then
                write(iout, *) &
                    "Warning: Overriding RealSpawnCutoff with scale function parameter"
                RealSpawnCutoff = sFBeta
            end if
        end if
        tNoSinglesPossible = t_k_space_hubbard .or. tUEG .or. tNoSinglesPossible

        if (.not. allocated(allInitsPerExLvl)) allocate(allInitsPerExLvl(maxInitExLvlWrite))
        if (.not. allocated(initsPerExLvl)) allocate(initsPerExLvl(maxInitExLvlWrite))
        initsPerExlvl = 0
        allInitsPerExLvl = 0

        IF (tHistSpawn .or. (tCalcFCIMCPsi .and. tFCIMC)) THEN
            allocate(HistMinInd(NEl))
            allocate(HistMinInd2(NEl))
            maxdet = 0
            do i = 1, nel
                maxdet = maxdet + 2**(nbasis - i)
            end do

            IF (.not. allocated(FCIDets)) THEN
                CALL Stop_All(t_r, "A Full Diagonalization is required before histogramming can occur.")
            end if

            write(iout, *) "Histogramming spawning wavevector, with Dets=", Det
            allocate(Histogram(1:lenof_sign, 1:det), stat=ierr)
            IF (ierr /= 0) THEN
                CALL Stop_All("SetupParameters", "Error assigning memory for histogramming arrays ")
            end if
            Histogram(:, :) = 0.0_dp
            allocate(AllHistogram(1:lenof_sign, 1:det), stat=ierr)
            allocate(BeforeNormHist(1:det), stat=ierr)
            IF (ierr /= 0) CALL Stop_All("SetupParameters", "Error assigning memory for histogramming arrays")
            IF (tHistSpawn) THEN
                allocate(InstHist(1:lenof_sign, 1:det), stat=ierr)
                IF (ierr /= 0) THEN
                    CALL Stop_All("SetupParameters", "Error assigning memory for histogramming arrays")
                end if
                InstHist(:, :) = 0.0_dp
                allocate(AvAnnihil(1:lenof_sign, 1:det), stat=ierr)
                IF (ierr /= 0) THEN
                    CALL Stop_All("SetupParameters", "Error assigning memory for histogramming arrays")
                end if
                AvAnnihil(:, :) = 0.0_dp
                allocate(InstAnnihil(1:lenof_sign, 1:det), stat=ierr)
                IF (ierr /= 0) THEN
                    CALL Stop_All("SetupParameters", "Error assigning memory for histogramming arrays")
                end if
                InstAnnihil(:, :) = 0.0_dp
            end if
            IF (iProcIndex == 0) THEN
                IF (tHistSpawn) THEN
                    Tot_Unique_Dets_Unit = get_free_unit()
                    open(Tot_Unique_Dets_Unit, FILE='TOTUNIQUEDETS', STATUS='UNKNOWN')
                    if (tCalcVariationalEnergy .and. tDiagAllSpaceEver) then
                        write(Tot_Unique_Dets_Unit, "(A)") "# Iter  UniqueDetsEver  AvVarEnergy   InstVarEnergy   GroundE_Ever"
                    else if (tCalcVariationalEnergy) then
                        write(Tot_Unique_Dets_Unit, "(A)") "# Iter  UniqueDetsEver  AvVarEnergy   InstVarEnergy"
                    else if (tDiagAllSpaceEver) then
                        write(Tot_Unique_Dets_Unit, "(A)") "# Iter  UniqueDetsEver  GroundE_Ever"
                    else
                        write(Tot_Unique_Dets_Unit, "(A)") "# Iter  UniqueDetsEver"
                    end if
                    allocate(AllInstHist(1:lenof_sign, 1:det), stat=ierr)
                    allocate(AllInstAnnihil(1:lenof_sign, 1:det), stat=ierr)
                    allocate(AllAvAnnihil(1:lenof_sign, 1:det), stat=ierr)
                end if
                IF (ierr /= 0) THEN
                    CALL Stop_All("SetupParameters", "Error assigning memory for histogramming arrays")
                end if
            end if
        else if (tHistEnergies) THEN
            write(iout, *) "Histogramming the energies of the particles, with iNoBins=", iNoBins, " and BinRange=", BinRange
            write(iout, *) "Histogramming spawning events from ", -OffDiagMax, " with BinRange = ", OffDiagBinRange
            iOffDiagNoBins = INT((2.0_dp * OffDiagMax) / OffDiagBinRange) + 1
            write(iout, *) "This gives ", iOffDiagNoBins, " bins to histogram the off-diagonal matrix elements."
            allocate(HistogramEnergy(1:iNoBins))
            allocate(AttemptHist(1:iNoBins))
            allocate(SpawnHist(1:iNoBins))
            allocate(SinglesHist(1:iOffDiagNoBins))
            allocate(SinglesAttemptHist(1:iOffDiagNoBins))
            allocate(SinglesHistOccOcc(1:iOffDiagNoBins))
            allocate(SinglesHistOccVirt(1:iOffDiagNoBins))
            allocate(SinglesHistVirtOcc(1:iOffDiagNoBins))
            allocate(SinglesHistVirtVirt(1:iOffDiagNoBins))
            allocate(DoublesHist(1:iOffDiagNoBins))
            allocate(DoublesAttemptHist(1:iOffDiagNoBins))
            HistogramEnergy(:) = 0.0_dp
            AttemptHist(:) = 0.0_dp
            SpawnHist(:) = 0.0_dp
            SinglesHist(:) = 0.0_dp
            SinglesAttemptHist(:) = 0.0_dp
            SinglesHistOccOcc(:) = 0.0_dp
            SinglesHistOccVirt(:) = 0.0_dp
            SinglesHistVirtOcc(:) = 0.0_dp
            SinglesHistVirtVirt(:) = 0.0_dp
            DoublesHist(:) = 0.0_dp
            DoublesAttemptHist(:) = 0.0_dp
            IF (iProcIndex == Root) THEN
                allocate(AllHistogramEnergy(1:iNoBins))
                allocate(AllAttemptHist(1:iNoBins))
                allocate(AllSpawnHist(1:iNoBins))
                allocate(AllSinglesHist(1:iOffDiagNoBins))
                allocate(AllDoublesHist(1:iOffDiagNoBins))
                allocate(AllSinglesAttemptHist(1:iOffDiagNoBins))
                allocate(AllDoublesAttemptHist(1:iOffDiagNoBins))
                allocate(AllSinglesHistOccOcc(1:iOffDiagNoBins))
                allocate(AllSinglesHistOccVirt(1:iOffDiagNoBins))
                allocate(AllSinglesHistVirtOcc(1:iOffDiagNoBins))
                allocate(AllSinglesHistVirtVirt(1:iOffDiagNoBins))
            end if
        end if

        if ((iProcIndex == Root) .and. tDiagWalkerSubspace) then
            write(iout, '(A,I9,A)') "Diagonalising walker subspace every ", iDiagSubspaceIter, " iterations"
            unitWalkerDiag = get_free_unit()
            open(unitWalkerDiag, file='WalkerSubspaceDiag', status='unknown')
            if (tTruncInitiator) then
                write(unitWalkerDiag, '(A)') "# Iter    NoInitDets   NoOccDets  InitiatorSubspaceEnergy   &
                        & FullSubspaceEnergy   ProjInitEnergy   ProjFullEnergy"
            else
                write(unitWalkerDiag, '(A)') "# Iter    NoOccDets    InitiatorSubspaceEnergy         &
                        & FullSubspaceEnergy    ProjFullEnergy"
            end if
        end if

        if (tHistExcitToFrom) &
            call init_hist_excit_tofrom()

        IF (tUseBrillouin) THEN
            write(iout, "(A)") "Brillouin theorem in use for calculation of projected energy."
        end if
        if (.not. (t_k_space_hubbard .and. t_trans_corr_2body)) then
            ! for too big lattices my implementation breaks, due to
            ! memory limitations.. but i think we do not actually need it.
            CALL CalcApproxpDoubles()
        end if
        IF (abs(TauFactor) > 1.0e-12_dp) THEN
            if (t_trans_corr_2body .and. t_k_space_hubbard) then
                call Stop_All(this_routine, &
                              "finding the number of excits from HF breaks for too large lattice")
            end if
            write(iout, *) "TauFactor detected. Resetting Tau based on connectivity of: ", HFConn
            Tau = TauFactor / REAL(HFConn, dp)
            write(iout, *) "Timestep set to: ", Tau
        end if

        ! [W.D.] I guess I want to initialize that before the tau-search,
        ! or otherwise some pgens get calculated incorrectly
        if (t_back_spawn .or. t_back_spawn_flex) then
            call init_back_spawn()
        end if

        if (t_guga_back_spawn) then
            call setup_virtual_mask()
        end if

        ! also i should warn the user if this is a restarted run with a
        ! set delay in the back-spawning method:
        ! is there actually a use-case where someone really wants to delay
        ! a back-spawn in a restarted run?
        if (tReadPops .and. back_spawn_delay /= 0) then
            call Warning_neci(t_r, &
                              "Do you really want a delayed back-spawn in a restarted run?")
        end if

        ! can i initialize the k-space hubbard here already?
        ! because we need information for the tau-search already..
        if (t_k_space_hubbard) then
            call init_k_space_hubbard()
        end if

        if (t_new_real_space_hubbard) then
            call init_real_space_hubbard
        end if

        if (t_spin_measurements) then
            call init_spin_measurements()
        end if

        IF (abs(StepsSftImag) > 1.0e-12_dp) THEN
            write(iout, *) "StepsShiftImag detected. Resetting StepsShift."
            StepsSft = NINT(StepsSftImag / Tau)
            IF (StepsSft == 0) StepsSft = 1
            write(iout, *) "StepsShift set to: ", StepsSft
        end if

        ! StepsPrint < 1 while not coupling update and output cycle means no std output
        if (.not. tCoupleCycleOutput .and. StepsPrint < 1) then
            tMCOutput = .false.
            ! But there shall be output in the FCIMCStats file
            ! There is no specified output cycle, so we default to the shift cycle -> couple them
            tCoupleCycleOutput = .true.
        end if
        ! Coupling output and shift update means these two are the same
        if (tCoupleCycleOutput) StepsPrint = StepsSft
        if (StepsPrint == StepsSft) tCoupleCycleOutput = .true.

        IF (TPopsFile) THEN
            IF (mod(iWritePopsEvery, StepsSft) /= 0) then
                CALL Warning_neci(t_r, "POPSFILE writeout should be a multiple of the update cycle length.")
            end if
        end if

        if (TReadPops) then
            if (tStartSinglePart .and. .not. tReadPopsRestart) then
                if (iProcIndex == root) call warning_neci(t_r, &
                                                          "ReadPOPS cannot work with StartSinglePart: ignoring StartSinglePart")
                tStartSinglePart = .false.
            end if
        end if

        IF (.not. TReadPops) THEN
            write(iout, "(A,F20.10)") "Initial Diagonal Shift: ", DiagSft(1)
        end if

        if (tShiftonHFPop) then
            write(iout, *) "Shift will be varied in order to keep the population on the reference determinant fixed"
        end if
        write(iout, *) "Connectivity of HF determinant is: ", HFConn
        IF (TStartSinglePart) THEN
            TSinglePartPhase(:) = .true.
        ELSE
            TSinglePartPhase(:) = .false.
        end if

        IF (ICILevel /= 0) THEN
!We are truncating the excitations at a certain value
            TTruncSpace = .true.
            write(iout, '(A,I4)') "Truncating the S.D. space at determinants will an excitation level w.r.t. HF of: ", ICILevel
        end if
        IF (tTruncCAS .or. tStartCAS) THEN
            ! We are truncating the FCI space by only allowing excitations
            ! in a predetermined CAS space.
            ! The following line has already been written out if we are doing
            ! a CAS calculation.

!            write(iout,'(A,I4,A,I5)') "Truncating the S.D. space as &
!                                   &determinants must be within a CAS of ", &
!                                   OccCASOrbs, " , ", VirtCASOrbs
            ! The SpinInvBRR array is required for the tTruncCAS option. Its
            ! properties are explained more fully in the subroutine.

            CALL CreateSpinInvBRR()

            ! CASmax is the max spin orbital number (when ordered
            ! energetically) within the chosen active space.
            ! Spin orbitals with energies larger than this maximum value must
            ! be unoccupied for the determinant to be in the active space.
            CASmax = NEl + VirtCASorbs

            ! CASmin is the max spin orbital number below the active space.
            ! As well as the above criteria, spin orbitals with energies
            ! equal to, or below that of the CASmin orbital must be completely
            ! occupied for the determinant to be in the active space.
            CASmin = NEl - OccCASorbs

            IF (OccCASOrbs > NEl) then
                CALL Stop_All("SetupParameters", "Occupied orbitals in CAS space specified is greater than number of electrons")
            end if
            IF (VirtCASOrbs > (nBasis - NEl)) then
                CALL Stop_All("SetupParameters", "Virtuals in CAS space specified greater than number of unoccupied orbitals")
            end if

!Create the bit masks for the bit calculation of these properties.
            allocate(CASMask(0:NIfD))
            allocate(CoreMask(0:NIfD))
            CASMask(:) = 0
            CoreMask(:) = 0
            do i = 1, nBasis
                IF (SpinInvBRR(i) > CASmax) THEN
                    !Orbital is in external space
                    CASMask((SpinInvBRR(i) - 1) / bits_n_int) = ibset(CASMask((i - 1) / bits_n_int), MOD((i - 1), bits_n_int))

                else if (SpinInvBRR(i) <= CASmin) THEN
                    !Orbital is in core space
                    CoreMask((SpinInvBRR(i) - 1) / bits_n_int) = ibset(CoreMask((i - 1) / bits_n_int), MOD((i - 1), bits_n_int))
                    CASMask((SpinInvBRR(i) - 1) / bits_n_int) = ibset(CASMask((i - 1) / bits_n_int), MOD((i - 1), bits_n_int))
                end if
            end do

        end if
        IF (tPartFreezeCore) THEN
            write(iout, '(A,I4,A,I5)') 'Partially freezing the lowest ', NPartFrozen, ' spin orbitals so that no more than ', &
                NHolesFrozen, ' holes exist within this core.'
            CALL CreateSpinInvBRR()
        end if
        IF (tPartFreezeVirt) THEN
            write(iout, '(A,I4,A,I5)') 'Partially freezing the highest ', NVirtPartFrozen, &
                ' virtual spin orbitals so that no more than ', NElVirtFrozen, ' electrons occupy these orbitals.'
            CALL CreateSpinInvBRR()
        end if

        if (tTruncNOpen) then
            write (iout, '("Truncating determinant space at a maximum of ",i3," &
                    &unpaired electrons.")') trunc_nopen_max
        end if

        ! for the (uniform) 3-body excitgen, the generation probabilities are uniquely given
        ! by the number of alpha and beta electrons and the number of orbitals
        ! and can hence be precomputed
        if (t_mol_3_body .or. t_ueg_3_body) call setup_mol_tc_excitgen()

        if (tGAS) then
            if (GAS_exc_gen == possible_GAS_exc_gen%DISCONNECTED) then
                call init_disconnected_GAS(GAS_specification)
            else if (GAS_exc_gen == possible_GAS_exc_gen%DISCARDING) then
                call init_GASCI_discarding(projEDet(:, 1))
            end if
        end if
    END SUBROUTINE SetupParameters

    ! This initialises the calculation, by allocating memory, setting up the
    ! initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        INTEGER :: ierr, iunithead
        logical :: formpops, binpops, tStartedFromCoreGround
        INTEGER :: error, MemoryAlloc, PopsVersion
        character(*), parameter :: t_r = 'InitFCIMCPar', this_routine = t_r
        integer :: PopBlockingIter
        real(dp) :: ExpectedMemWalk, read_tau, read_psingles, read_pparallel, read_ptriples
        integer(int64) :: read_walkers_on_nodes(0:nProcessors - 1)
        integer :: read_nnodes
        !Variables from popsfile header...
        logical :: tPop64Bit, tPopHPHF, tPopLz
        integer :: iPopLenof_sign, iPopNel, iPopIter, PopNIfD, PopNIfSgn, PopNIfFlag, PopNIfTot, Popinum_runs
        integer :: PopRandomHash(2056), PopBalanceBlocks
        integer(int64) :: iPopAllTotWalkers
        integer :: i, run
        real(dp) :: PopDiagSft(1:inum_runs)
        real(dp), dimension(lenof_sign) :: PopSumNoatHF
        HElement_t(dp) :: PopAllSumENum(1:inum_runs)
        integer :: perturb_ncreate, perturb_nannihilate
        integer :: nrdms_standard, nrdms_transition
        character(255) :: identifier
        ! set default pops version, this should not have any functional impact,
        ! but prevents using it uninitialized
        !default
        Popinum_runs = 1
        ! default version for popsfiles, this does not have any functional effect,
        ! but prevents it from using uninitialized
        PopsVersion = 4

        if (tPopsAlias) then
            identifier = aliasStem
        else
            identifier = 'POPSFILE'
        end if

        if (tReadPops .and. .not. (tPopsAlreadyRead .or. tHDF5PopsRead)) then
            call open_pops_head(iunithead, formpops, binpops, identifier)
            PopsVersion = FindPopsfileVersion(iunithead)
            if (iProcIndex == root) close(iunithead)
            write(iout, *) "POPSFILE VERSION ", PopsVersion, " detected."
        end if

        if (tReadPops .and. (PopsVersion < 3) .and. &
            .not. (tPopsAlreadyRead .or. tHDF5PopsRead)) then
!Read in particles from multiple POPSFILES for each processor
            !Ugh - need to set up ValidSpawnedList here too...
            call SetupValidSpawned(int(InitWalkers, int64))
            write(iout, *) "Reading in initial particle configuration from *OLD* POPSFILES..."
            CALL ReadFromPopsFilePar()
        ELSE
            !Scale walker number
            !This is needed to be done here rather than later,
            !because the arrays should be allocated with appropariate sizes
            if (tReadPops .and. .not. tPopsAlreadyRead) then
                InitWalkers = InitWalkers * ScaleWalkers
            end if

!initialise the particle positions - start at HF with positive sign
!Set the maximum number of walkers allowed
            if (tReadPops .and. .not. (tPopsAlreadyRead .or. tHDF5PopsRead)) then
                !Read header.
                call open_pops_head(iunithead, formpops, binpops, identifier)
                if (PopsVersion == 3) then
                    call ReadPopsHeadv3(iunithead, tPop64Bit, tPopHPHF, tPopLz, iPopLenof_Sign, iPopNel, &
                                        iPopAllTotWalkers, PopDiagSft, PopSumNoatHF, PopAllSumENum, iPopIter, &
                                        PopNIfD, PopNIfSgn, PopNIfFlag, PopNIfTot)

                    ! The following values were not read in...
                    read_tau = 0.0_dp
                    read_nnodes = 0
                    PopBalanceBlocks = -1
                else if (PopsVersion == 4) then
                    ! The only difference between 3 & 4 is just that 4 reads
                    ! in via a namelist, so that we can add more details
                    ! whenever we want.
                    call ReadPopsHeadv4(iunithead, tPop64Bit, tPopHPHF, tPopLz, iPopLenof_Sign, iPopNel, &
                                        iPopAllTotWalkers, PopDiagSft, PopSumNoatHF, PopAllSumENum, iPopIter, &
                                        PopNIfD, PopNIfSgn, Popinum_runs, PopNIfFlag, PopNIfTot, &
                                        read_tau, PopBlockingIter, PopRandomHash, read_psingles, &
                                        read_pparallel, read_ptriples, read_nnodes, read_walkers_on_nodes, &
                                        PopBalanceBlocks)
                    ! Use the random hash from the Popsfile. This is so that,
                    ! if we are using the same number of processors as the job
                    ! which produced the Popsfile, we can send the determinants
                    ! to the correct processor directly (see read_pops_nnodes).
                    ! This won't work if a different random number seed has
                    ! been used, so we copy the random hash across for safety.
                    ! This will overwrite current value of RandomHash.
                    !RandomHash = PopRandomHash(1:nBasis)
                    !call MPIBcast(RandomHash)
                else
                    call stop_all(this_routine, "Popsfile version invalid")
                end if

                ! Check the number of electrons created and annihilated by the
                ! perturbation operators, if any are being used.
                if (allocated(pops_pert)) then
                    perturb_ncreate = pops_pert(1)%ncreate
                    perturb_nannihilate = pops_pert(1)%nannihilate
                else
                    perturb_ncreate = 0
                    perturb_nannihilate = 0
                end if

                call CheckPopsParams(tPop64Bit, tPopHPHF, tPopLz, iPopLenof_Sign, iPopNel, &
                                     iPopAllTotWalkers, PopDiagSft, PopSumNoatHF, PopAllSumENum, iPopIter, &
                                     PopNIfD, PopNIfSgn, PopNIfTot, &
                                     MaxWalkersUncorrected, read_tau, PopBlockingIter, read_psingles, read_pparallel, &
                                     read_ptriples, perturb_ncreate, perturb_nannihilate)

                if (iProcIndex == root) close(iunithead)
            else
                MaxWalkersUncorrected = int(InitWalkers, sizeof_int)
            end if

            MaxWalkersPart = NINT(MemoryFacPart * MaxWalkersUncorrected)
            ExpectedMemWalk = real((NIfTot + 1) * MaxWalkersPart * size_n_int + 8 * MaxWalkersPart, dp) / 1048576.0_dp
            if (ExpectedMemWalk < 20.0) then
                !Increase memory allowance for small runs to a min of 20mb
                MaxWalkersPart = int(20.0 * 1048576.0 / real((NIfTot + 1) * size_n_int + 8, dp), sizeof_int)
                write(iout, "(A)") "Low memory requested for walkers, so increasing memory to 20Mb to avoid memory errors"
            end if
            write(iout, "(A,I14)") "Memory allocated for a maximum particle number per node of: ", MaxWalkersPart
            !Here is where MaxSpawned is set up - do we want to set up a minimum allocation here too?
            Call SetupValidSpawned(int(InitWalkers, int64))

!Put a barrier here so all processes synchronise
            CALL MPIBarrier(error)
!Allocate memory to hold walkers
            allocate(WalkVecDets(0:NIfTot, MaxWalkersPart), stat=ierr)
            CALL LogMemAlloc('WalkVecDets', MaxWalkersPart * (NIfTot + 1), size_n_int, this_routine, WalkVecDetsTag, ierr)
            WalkVecDets(0:NIfTot, 1:MaxWalkersPart) = 0
            MemoryAlloc = (NIfTot + 1) * MaxWalkersPart * size_n_int    !Memory Allocated in bytes

            if (alloc_popsfile_dets) allocate(popsfile_dets(0:NIfTot, MaxWalkersPart), stat=ierr)

            nrdms_standard = 0
            nrdms_transition = 0

            if (tRDMOnFly) then
                if (tPairedReplicas) then
                    nrdms_standard = lenof_sign.div.2
                else
                    nrdms_standard = lenof_sign
                end if
                if (tTransitionRDMs) then
                    ! nrdms_transition_input will have been read in from the user
                    ! input. But if we have two different replicas for each state
                    ! sampled, then there are two ways to form each transition RDMs.
                    nrdms_transition = nrdms_transition_input * nreplicas
                end if
            end if

            ! Allocate storage for persistent data to be stored alongside
            ! the current determinant list (particularly diagonal matrix
            ! elements, i.e. CurrentH; now global_determinant_data).
            call init_global_det_data(nrdms_standard, nrdms_transition)

            ! If we are doing cont time, then initialise it here
            call init_cont_time()

            ! set the dummies for trial wavefunction connected space
            ! load balancing before trial wf initialization
            if (tTrialWavefunction) then
                allocate(con_send_buf(0, 0))
                con_space_size = 0
            end if

            write(iout, "(A,I12,A)") "Spawning vectors allowing for a total of ", MaxSpawned, &
                " particles to be spawned in any one iteration per core."
            write(iout, *) "Memory requirement ", IlutBits%len_bcast * 8.0_dp * ( &
                MaxSpawned / 1048576.0_dp), "MB"

            allocate(SpawnVec(0:IlutBits%len_bcast, MaxSpawned), &
                      stat=ierr, source=0_n_int)
            @:log_alloc(SpawnVec, SpawnVecTag, ierr)
            allocate(SpawnVec2(0:IlutBits%len_bcast, MaxSpawned), &
                      stat=ierr, source=0_n_int)
            @:log_alloc(SpawnVec2, SpawnVec2Tag, ierr)

            if (use_spawn_hash_table) then
                nhashes_spawn = ceiling(0.8 * MaxSpawned)
                allocate(spawn_ht(nhashes_spawn), stat=ierr)
                call init_hash_table(spawn_ht)
            end if

!Point at correct spawning arrays
            SpawnedParts => SpawnVec
            SpawnedParts2 => SpawnVec2

            MemoryAlloc = MemoryAlloc + (NIfTot + 1) * MaxSpawned * 2 * size_n_int

            if (tAutoAdaptiveShift) then
                allocate(SpawnInfoVec(1:SpawnInfoWidth, MaxSpawned), stat=ierr)
                @:log_alloc(SpawnInfoVec, SpawnInfoVecTag, ierr)
                allocate(SpawnInfoVec2(1:SpawnInfoWidth, MaxSpawned), stat=ierr)
                @:log_alloc(SpawnInfoVec2, SpawnInfoVec2Tag, ierr)
                SpawnInfoVec(:, :) = 0
                SpawnInfoVec2(:, :) = 0
                SpawnInfo => SpawnInfoVec
                SpawnInfo2 => SpawnInfoVec2
                MemoryAlloc = MemoryAlloc + (SpawnInfoWidth) * MaxSpawned * 2 * size_n_int
            end if

            write(iout, "(A)") "Storing walkers in hash-table. Algorithm is now formally linear scaling with walker number"
            write(iout, "(A,I15)") "Length of hash-table: ", nWalkerHashes
            write(iout, "(A,F20.5)") "Length of hash-table as a fraction of targetwalkers: ", HashLengthFrac
            ! TODO: Correct the memory usage.
            !MemTemp=2*(8*(nClashMax+1)*nWalkerHashes)+8*MaxWalkersPart
            !write(iout,"(A,F10.3,A)") "This will use ",real(MemTemp,dp)/1048576.0_dp,&
            !    " Mb of memory per process, although this is likely to increase as it expands"
            allocate(HashIndex(nWalkerHashes), stat=ierr)
            if (ierr /= 0) call stop_all(this_routine, "Error in allocation")
            do i = 1, nWalkerHashes
                HashIndex(i)%ind = 0
            end do
            !Also need to allocate memory for the freeslot array
            allocate(FreeSlot(MaxWalkersPart), stat=ierr)
            if (ierr /= 0) call stop_all(this_routine, "Error in allocation")
            freeslot(:) = 0
            !MemoryAlloc=MemoryAlloc+MemTemp

            ! Allocate pointers to the correct walker arrays
            CurrentDets => WalkVecDets

            ! Get the (0-based) processor index for the HF det.
            do run = 1, inum_runs
                iRefProc(run) = DetermineDetNode(nel, ProjEDet(:, run), 0)
            end do
            write(iout, "(A,I8)") "Reference processor is: ", iRefProc(1)
            write(iout, "(A)", advance='no') "Initial reference is: "
            call write_det(iout, ProjEDet(:, 1), .true.)

            TotParts(:) = 0.0
            TotPartsOld(:) = 0.0
            NoatHF(:) = 0.0_dp
            InstNoatHF(:) = 0.0_dp

            ! Has been moved to guarantee initialization before first load balancing
            ! Initialises RDM stuff for both explicit and stochastic calculations of RDM.

            tFillingStochRDMonFly = .false.
            tFillingExplicRDMonFly = .false.
            !One of these becomes true when we have reached the relevant iteration to begin filling the RDM.

            ! If we have a popsfile, read the walkers in now.
            if (tReadPops .and. .not. tPopsAlreadyRead) then
                call InitFCIMC_pops(iPopAllTotWalkers, PopNIfSgn, iPopNel, read_nnodes, &
                                    read_walkers_on_nodes, pops_pert, &
                                    PopBalanceBLocks, PopDiagSft)
            else
                if (tStartMP1) then
                    !Initialise walkers according to mp1 amplitude.
                    call InitFCIMC_MP1()

                else if (tStartCAS) then
                    !Initialise walkers according to a CAS diagonalisation.
                    call InitFCIMC_CAS()

                else if (tTrialInit .or. (tOrthogonaliseReplicas .and. &
                                          .not. tReplicaSingleDetStart)) then

                    call InitFCIMC_trial()

                else if (tInitializeCSF) then

                    call InitFCIMC_CSF()

                else !Set up walkers on HF det

                    if (tStartSinglePart) then
                        write(iout, "(A,I10,A,F9.3,A,I15)") "Initial number of particles set to ", int(InitialPart), &
                            " and shift will be held at ", DiagSft(1), " until particle number gets to ", int(InitWalkers * nNodes)
                    else
                        write(iout, "(A,I16)") "Initial number of walkers per processor chosen to be: ", nint(InitWalkers)
                    end if

                    call InitFCIMC_HF()

                end if   !tStartmp1
            end if
            if (t_ueg_3_body .and. tTrcorrExgen) tLatticeGens = .false.
            write(iout, "(A,F14.6,A)") " Initial memory (without excitgens + temp arrays) consists of : ", &
                & REAL(MemoryAlloc, dp) / 1048576.0_dp, " Mb/Processor"
            write(iout, *) "Only one array of memory to store main particle list allocated..."
            write(iout, *) "Initial memory allocation sucessful..."
            write(iout, *) "============================================="
            CALL neci_flush(iout)

        end if   !End if initial walkers method

        ! There was no last output, use the same value as for the shift update
        AllTotPartsLastOutput = AllTotPartsOld
!Put a barrier here so all processes synchronise
        CALL MPIBarrier(error)

        call init_norm()

        IF (tPrintOrbOcc) THEN
            allocate(OrbOccs(nBasis), stat=ierr)
            CALL LogMemAlloc('OrbOccs', nBasis, 8, this_routine, OrbOccsTag, ierr)
            OrbOccs(:) = 0.0_dp
        end if

        IF (tHistInitPops) THEN
            CALL InitHistInitPops()
        end if
        tPrintHighPop = .false.
        MaxInitPopPos = 0.0
        MaxInitPopNeg = 0.0

        IF (abs(MaxNoatHF) < 1.0e-12_dp) THEN
            MaxNoatHF = InitWalkers * nNodes
            HFPopThresh = int(MaxNoatHF, int64)
        end if

        ! Initialise excitation generation storage
        call init_excit_gen_store(fcimc_excit_gen_store)

        ! initialize excitation generator
        if (t_pcpp_excitgen) call init_pcpp_excitgen()
        if (t_pchb_excitgen) call init_pchb_excitgen(projEDet(:, 1))
        ! [W.D.] I guess I want to initialize that before the tau-search,
        ! or otherwise some pgens get calculated incorrectly
        if (t_back_spawn .or. t_back_spawn_flex) then
            call init_back_spawn()
        end if
        ! also i should warn the user if this is a restarted run with a
        ! set delay in the back-spawning method:
        ! is there actually a use-case where someone really wants to delay
        ! a back-spawn in a restarted run?
        if (tReadPops .and. back_spawn_delay /= 0) then
            call Warning_neci(t_r, &
                              "Do you really want a delayed back-spawn in a restarted run?")
        end if

        ! [Werner Dobrautz 5.5.2017:]
        ! if this is a continued run from a histogramming tau-search
        ! and a restart of the tau-search is not forced by input, turn
        ! both the new and the old tau-search off!
        ! i cannot do it here, since this is called before the popsfile read-in
        if (t_previous_hist_tau) then
            ! i have to check for tau-search option and stuff also, so that
            ! the death tau adaption is still used atleast! todo!
            tSearchTau = .false.
            t_hist_tau_search = .false.
            t_fill_frequency_hists = .false.
            write(iout, *) "Turning OFF the tau-search, since continued run!"
        end if

        if (tSearchTau) then
            if (tGen_nosym_guga) then
                call init_tau_search_guga_nosym()
            else
                call init_tau_search()
            end if
            ! set!
            ! [Werner Dobrautz 4.4.2017:]
            if (t_hist_tau_search) then
                ! some setup went wrong!
                call Stop_All(t_r, &
                              "Input error! both standard AND Histogram tau-search chosen!")
            end if

        else if (t_hist_tau_search) then
            if (tGen_nosym_guga) then
                call init_hist_tau_search_guga_nosym()
            else
                call init_hist_tau_search()
            end if

        else if (t_hist_tau_search) then
            call init_hist_tau_search()

        else
            ! Add a couple of checks for sanity
            if (nOccAlpha == 0 .or. nOccBeta == 0) then
                pParallel = 1.0_dp
            end if
            if (nOccAlpha == 1 .and. nOccBeta == 1) then
                pParallel = 0.0_dp
            end if
        end if

        IF ((NMCyc /= 0) .and. (tRotateOrbs .and. (.not. tFindCINatOrbs))) then
            CALL Stop_All(this_routine, "Currently not set up to rotate and then go straight into a spawning &
            & calculation.  Ordering of orbitals is incorrect.  This may be fixed if needed.")
        end if

        if (tRDMonFly) then
            call init_rdms(nrdms_standard, nrdms_transition)
            !If the iteration specified to start filling the RDM has already been, want to
            !start filling as soon as possible.
            do run = 1, inum_runs
                if (.not. tSinglePartPhase(run)) VaryShiftIter(run) = 0
            end do
        end if

        ! Perform all semi-stochastic initialisation. This includes generating all the states in the
        ! deterministic space, finding their processors, ordering them, inserting them into
        ! CurrentDets, calculating and storing all Hamiltonian matrix elements and initalising all
        ! arrays required to store and distribute the vectors in the deterministic space later.

        ! in the real-time application, this is done after the initial state is set up
        if (tSemiStochastic .and. .not. t_real_time_fciqmc) then
            if (tDynamicCoreSpace .and. tRDMonFly) then
                tSemiStochastic = .false.
                semistoch_shift_iter = 1
            else
                call init_semi_stochastic(ss_space_in, tStartedFromCoreGround)
                if (tStartedFromCoreGround .and. tSetInitialRunRef) call set_initial_run_references()
            end if
        end if

        ! If the number of trial states to calculate hasn't been set by the
        ! user, then simply use the minimum number
        if ((tTrialWavefunction .or. tStartTrialLater) .and. (ntrial_ex_calc == 0)) then
            ntrial_ex_calc = inum_runs
        end if

        ! Initialise the trial wavefunction information which can be used for the energy estimator.
        ! This includes generating the trial space, generating the space connected to the trial space,
        ! diagonalising the trial space to find the trial wavefunction and calculating the vector
        ! in the connected space, required for the energy estimator.
        if (tRDMonFly .and. tDynamicCoreSpace .and. tTrialWavefunction) then
            tTrialWavefunction = .false.
            tStartTrialLater = .true.
            trial_shift_iter = 1
        end if
        if (tTrialWavefunction) then
            if (tPairedReplicas) then
                call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs.div.2, .true.)
            else
                call init_trial_wf(trial_space_in, ntrial_ex_calc, inum_runs, .false.)
            end if
        else if (tStartTrialLater) then
            ! If we are going to turn on the use of a trial wave function
            ! later in the calculation, then zero the trial estimate arrays
            ! for now, to prevent junk being printed before then.
            trial_numerator = 0.0_dp
            tot_trial_numerator = 0.0_dp
            trial_denom = 0.0_dp
            tot_trial_denom = 0.0_dp

            init_trial_numerator = 0.0_dp
            tot_init_trial_numerator = 0.0_dp
            init_trial_denom = 0.0_dp
            tot_init_trial_denom = 0.0_dp

            trial_num_inst = 0.0_dp
            tot_trial_num_inst = 0.0_dp
            trial_denom_inst = 0.0_dp
            tot_trial_denom_inst = 0.0_dp
        end if

        replica_overlaps_real(:, :) = 0.0_dp
#ifdef CMPLX_
        replica_overlaps_imag(:, :) = 0.0_dp
#endif

        if (t_cepa_shift) call init_cepa_shifts()

        ! Set up the reference space for the adi-approach
        ! in real-time, we do this in the real-time init
        if (.not. t_real_time_fciqmc) call setup_reference_space(tReadPops)

        ! in fixed-n0, the variable shift mode and everything connected is
        ! controlled over the reference population
        if (tFixedN0) then
            if (tReadPops .and. .not. tWalkContGrow) then
                tSkipRef = .true.
                tSinglePartPhase = .false.
            else
                tSkipRef = .false.
                tSinglePartPhase = .true.
            end if
        end if

        if (tRDMonFly .and. tDynamicCoreSpace) call sync_rdm_sampling_iter()

        ! for the (uniform) 3-body excitgen, the generation probabilities are uniquely given
        ! by the number of alpha and beta electrons and the number of orbitals
        ! and can hence be precomputed
        if (t_mol_3_body .or. t_ueg_3_body) call setup_mol_tc_excitgen()

        if (tInitiatorSpace) call init_initiator_space(i_space_in)

        if (tReplicaEstimates) then
            if (.not. tPairedReplicas) then
                call stop_all(this_routine, "The paired-replicas option must be used the logging &
                                            &block, in order to calculate replica estimates.)")
            end if

            if (tSemiStochastic) allocate (tDetermSpawnedTo( &
                                           cs_replicas(core_run)%determ_sizes(iProcIndex)))

            call open_replica_est_file()
        end if

        ! When should we perform death before communication?
        ! For integer walkers, do death before comms just so the tests don't fail (or need updating).
        if (.not. tAllRealCoeff) then
            tDeathBeforeComms = .true.
        end if
        if (t_back_spawn .or. t_back_spawn_flex) then
            tDeathBeforeComms = .true.
        end if

        ! For FCIQMC with preconditioning and a time step of 1, death will
        ! kill all walkers and remove them from the hash table. In this
        ! case, we must set the initiator flags for spawning to occupied
        ! determinants before this occurs.
        if (tPreCond .and. (tau.isclose.1.0_dp)) tSetInitFlagsBeforeDeath = .true.

        ! Make sure we are performing death *after* communication, in cases
        ! where this is essential.
        if (tPreCond .and. tDeathBeforeComms) then
            call stop_all(this_routine, "With preconditioning, death must &
                                        &be performed after communication.")
        end if
        if (tReplicaEstimates .and. tDeathBeforeComms) then
            call stop_all(this_routine, "In order to calculate replica estimates, &
                                        &death must be performed after communication.")
        end if
        if (tEN2Init .and. (.not. tTruncInitiator)) then
            call stop_all(this_routine, "Cannot calculate the EN2 correction to initiator &
                                        &error as the initiator method is not in use.")
        end if

    end subroutine InitFCIMCCalcPar

    subroutine init_fcimc_fn_pointers()
        character(*), parameter :: this_routine = "init_fcimc_fn_pointers"
        character(*), parameter :: t_r = 'init_fcimc_fn_pointers'

        ! Select the excitation generator.
        if (tHPHF) then
            generate_excitation => gen_hphf_excit
        else if (tGAS) then
            if (GAS_exc_gen == possible_GAS_exc_gen%GENERAL) then
                generate_excitation => gen_GASCI_general
            else if (GAS_exc_gen == possible_GAS_exc_gen%DISCONNECTED) then
                generate_excitation => gen_GASCI_disconnected
            else if (GAS_exc_gen == possible_GAS_exc_gen%DISCARDING) then
                generate_excitation => gen_GASCI_discarding
            else
                call stop_all(t_r, 'Invalid GAS excitation generator')
            end if
        else if (t_3_body_excits .and. .not. (t_mol_3_body .or. t_ueg_3_body)) then
            if (t_uniform_excits) then
                generate_excitation => gen_excit_uniform_k_space_hub_transcorr
            else
                generate_excitation => gen_excit_k_space_hub_transcorr
            end if
        else if (tHPHF) then
            generate_excitation => gen_hphf_excit
        else if (t_ueg_3_body) then
            if (tTrcorrExgen) then
                generate_two_body_excitation => gen_ueg_excit
            else if (TLatticeGens) then
                generate_two_body_excitation => gen_rand_excit
            end if
            generate_excitation => gen_excit_mol_tc
        else if ((t_back_spawn_option .or. t_back_spawn_flex_option)) then
            if (tHUB .and. tLatticeGens) then
                ! for now the hubbard + back-spawn still uses the old
                ! genrand excit gen
                generate_excitation => gen_excit_back_spawn_hubbard
            else if (tUEGNewGenerator .and. tLatticeGens) then
                generate_excitation => gen_excit_back_spawn_ueg_new
            else if (tUEG .and. tLatticeGens) then
                generate_excitation => gen_excit_back_spawn_ueg
            else
                generate_excitation => gen_excit_back_spawn
            end if
        else if (tUEGNewGenerator) then
            generate_excitation => gen_ueg_excit
        else if (tPickVirtUniform) then
            ! pick-uniform-random-mag is on
            if (tReltvy) then
                generate_excitation => gen_rand_excit_Ex_Mag
            else
                call stop_all(t_r, "Excitation generator has not been set!")
            end if
        else if (tGenHelWeighted) then
            generate_excitation => gen_excit_hel_weighted
        else if (tGen_4ind_2) then
            generate_excitation => gen_excit_4ind_weighted2
        else if (tGen_4ind_weighted) then
            generate_excitation => gen_excit_4ind_weighted
        else if (tGen_4ind_reverse) then
            generate_excitation => gen_excit_4ind_reverse

        else if (tGUGA) then
            if (tgen_guga_crude) then
                generate_excitation => gen_excit_4ind_weighted2
            else
                generate_excitation => generate_excitation_guga
            end if

        else if (t_pcpp_excitgen) then
            generate_excitation => gen_rand_excit_pcpp
        else if (t_pchb_excitgen) then
            generate_excitation => gen_rand_excit_pchb
        else
            generate_excitation => gen_rand_excit
        end if

        ! if we are using the 3-body excitation generator, embed the chosen excitgen
        ! in the three-body one
        if (t_mol_3_body) then
            ! yes, fortran pointers work this way
            generate_two_body_excitation => generate_excitation
            generate_excitation => gen_excit_mol_tc
        end if

        ! In the main loop, we only need to find out if a determinant is
        ! connected to the reference det or not (so no ex. level above 2 is
        ! required). Except in some cases where we need to know the maximum
        ! excitation level
        if (tTruncSpace .or. tHistSpawn .or. tCalcFCIMCPsi) then
            max_calc_ex_level = nel
        else
            if (t_3_body_excits) then
                max_calc_ex_level = 3
            else
                max_calc_ex_level = 2
            end if
        end if

        ! How many children should we spawn given an excitation?
        if (t_real_time_fciqmc) then
            attempt_create => attempt_create_realtime
        else if (tTruncCas .or. tTruncSpace .or. &
                 tPartFreezeCore .or. tPartFreezeVirt .or. tFixLz .or. &
                 (tUEG .and. .not. tLatticeGens) .or. tTruncNOpen .or. t_trunc_nopen_diff) then
            if (tHPHF .or. tSemiStochastic) then
                attempt_create => attempt_create_trunc_spawn
            else
                attempt_create => att_create_trunc_spawn_enc
            end if
        else
            attempt_create => attempt_create_normal
        end if

        ! In attempt create, how should we evaluate the off diagonal matrix
        ! elements between a parent and its (potentially) spawned offspring?
        if (tHPHF) then
            if (tGenMatHEL) then
                get_spawn_helement => hphf_spawn_sign
            else
                get_spawn_helement => hphf_off_diag_helement_spawn
            end if

            ! new guga addition: do not need to recalculate Helement
        else if (tGUGA) then
            ! use hphf_routine also, since it does exactly what needed
            get_spawn_helement => hphf_spawn_sign
        else
            get_spawn_helement => get_helement_det_only
        end if

        ! When calling routines to generate all possible connections, this
        ! routine is called to generate the corresponding Hamiltonian matrix
        ! elements.
        if (tHPHF) then
            get_conn_helement => hphf_off_diag_helement_spawn
        else
            get_conn_helement => get_helement_det_only
        end if

        ! Once we have generated the children, do we need to encode them?
        if (.not. (tHPHF .or. tGen_4ind_weighted .or. tGUGA)) then
            encode_child => FindExcitBitDet
        else
            encode_child => null_encode_child
        end if

        ! What message should we display for a particle bloom?
        if (tAddToInitiator) then
            bloom_warn_string = '("Bloom of more than n_add on ", a, " excit: &
                                &A max of ", f10.2, " particles created. ", &
                                &i8, " blooms occurred.")'
        else
            ! Use this variable to store the bloom cutoff level.
            InitiatorWalkNo = 3.0_dp
            bloom_warn_string = '("Bloom of more than 3 on ", a, " excit: &
                                &A max of ", f10.2, " particles created. ", &
                                &i8, " blooms occurred.")'
        end if
        bloom_max = 0

        ! Perform the correct statistics on new child particles
        new_child_stats => new_child_stats_normal

        if (tPreCond) then
            attempt_die => attempt_die_precond
        else
            attempt_die => attempt_die_normal
        end if

        extract_bit_rep_avsign => extract_bit_rep_avsign_no_rdm

        fill_rdm_diag_currdet => fill_rdm_diag_currdet_norm

        select case (sfTag)
        case (0)
            scaleFunction => powerScaleFunction
        case (1)
            scaleFunction => expScaleFunction
        case (2)
            scaleFunction => negScaleFunction
        case (3)
            scaleFunction => expCOScaleFunction
        case default
            call stop_all(t_r, "Invalid scale function specified")
        end select

        if (tExpAdaptiveShift) then
            shiftFactorFunction => expShiftFactorFunction
        else if (tLinearAdaptiveShift) then
            shiftFactorFunction => linearShiftFactorFunction
        else if (tAutoAdaptiveShift) then
            shiftFactorFunction => autoShiftFactorFunction
        else
            shiftFactorFunction => constShiftFactorFunction
        end if

        ! select the procedure that returns all connected determinants.
        if (tGAS) then
            gen_all_excits => gen_all_excits_GAS
        else if (t_k_space_hubbard) then
            gen_all_excits => gen_all_excits_k_space_hubbard
        else if (t_new_real_space_hubbard) then
            gen_all_excits => gen_all_excits_r_space_hubbard
        else
            gen_all_excits => gen_all_excits_default
        end if

    end subroutine init_fcimc_fn_pointers

    subroutine DeallocFCIMCMemPar()

        CHARACTER(len=*), PARAMETER :: this_routine = 'DeallocFciMCMemPar'
        type(ll_node), pointer :: Curr, Prev
        integer :: i, ierr

        deallocate(RandomHash2, stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Err deallocating")

        ! Deallocate the linked list
        do i = 1, nWalkerHashes
            Curr => HashIndex(i)%Next
            Prev => HashIndex(i)
            nullify (Prev%Next)
            do while (associated(Curr))
                Prev => Curr
                Curr => Curr%Next
                deallocate(Prev)
                if (ierr /= 0) call stop_all(this_routine, "Err deallocating")
            end do
        end do
        deallocate(HashIndex, stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Err deallocating")
        nullify (Curr)
        nullify (Prev)

        deallocate(FreeSlot, stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Err deallocating")

        IF (tHistSpawn .or. tCalcFCIMCPsi) THEN
            DEallocate(Histogram)
            DEallocate(AllHistogram)
            IF (tHistSpawn) THEN
                DEallocate(InstHist)
                DEallocate(InstAnnihil)
                DEallocate(AvAnnihil)
            end if
            IF (iProcIndex == 0) THEN
                IF (tHistSpawn) THEN
                    DEallocate(AllInstHist)
                    DEallocate(AllAvAnnihil)
                    DEallocate(AllInstAnnihil)
                end if
            end if
        else if (tHistEnergies) THEN
            DEallocate(HistogramEnergy)
            DEallocate(AttemptHist)
            DEallocate(SpawnHist)
            DEallocate(SinglesHist)
            DEallocate(DoublesHist)
            DEallocate(DoublesAttemptHist)
            DEallocate(SinglesAttemptHist)
            DEallocate(SinglesHistOccOcc)
            DEallocate(SinglesHistVirtOcc)
            DEallocate(SinglesHistOccVirt)
            DEallocate(SinglesHistVirtVirt)
            IF (iProcIndex == Root) THEN
                DEallocate(AllHistogramEnergy)
                DEallocate(AllAttemptHist)
                DEallocate(AllSpawnHist)
                DEallocate(AllSinglesAttemptHist)
                DEallocate(AllSinglesHist)
                DEallocate(AllDoublesAttemptHist)
                DEallocate(AllDoublesHist)
                DEallocate(AllSinglesHistOccOcc)
                DEallocate(AllSinglesHistVirtOcc)
                DEallocate(AllSinglesHistOccVirt)
                DEallocate(AllSinglesHistVirtVirt)
            end if
        end if
        if (tHistExcitToFrom) call clean_hist_excit_tofrom()
        DEallocate(WalkVecDets)
        CALL LogMemDealloc(this_routine, WalkVecDetsTag)
        DEallocate(SpawnVec)
        CALL LogMemDealloc(this_routine, SpawnVecTag)
        DEallocate(SpawnVec2)
        CALL LogMemDealloc(this_routine, SpawnVec2Tag)
        if (tAutoAdaptiveShift) then
            DEallocate(SpawnInfoVec)
            CALL LogMemDealloc(this_routine, SpawnInfoVecTag)
            DEallocate(SpawnInfoVec2)
            CALL LogMemDealloc(this_routine, SpawnInfoVec2Tag)
        end if

        if (allocated(TempSpawnedParts)) then
            deallocate(TempSpawnedParts)
            log_dealloc(TempSpawnedPartsTag)
        end if
        DEallocate(HFDet)
        CALL LogMemDealloc(this_routine, HFDetTag)
        DEallocate(iLutHF)
        DEallocate(iLutRef)
        DEallocate(ProjEDet)
        DEallocate(iLutHF_True)
        DEallocate(HFDet_True)
        IF (ALLOCATED(HighestPopDet)) DEallocate(HighestPopDet)
        IF (ALLOCATED(RandomOrbIndex)) DEallocate(RandomOrbIndex)

        IF (ALLOCATED(SpinInvBrr)) THEN
            CALL LogMemDealloc(this_routine, SpinInvBRRTag)
            DEallocate(SpinInvBRR)
        end if
        IF (ALLOCATED(CoreMask)) THEN
            DEallocate(CoreMask)
            DEallocate(CASMask)
        end if
        IF (tPrintOrbOcc) THEN
            DEallocate(OrbOccs)
            CALL LogMemDeAlloc(this_routine, OrbOccsTag)
        end if

        IF (tHistInitPops) THEN
            if (allocated(HistInitPops)) then
                deallocate(HistInitPops)
                call LogMemDeAlloc(this_routine, HistInitPopsTag)
            end if
            IF (iProcIndex == 0) THEN
                if (allocated(AllHistInitPops)) then
                    deallocate(AllHistInitPops)
                    call LogMemDeAlloc(this_routine, AllHistInitPopsTag)
                end if
            end if
        end if

        if (tHub) then
            if (allocated(momIndexTable)) deallocate(momIndexTable)
            deallocate(breathingCont)
        end if

        if (tRDMonFly) call dealloc_global_rdm_data()

        if (allocated(refdetflip)) deallocate(refdetflip)
        if (allocated(ilutrefflip)) deallocate(ilutrefflip)
        if (allocated(ValidSpawnedList)) deallocate(ValidSpawnedList)
        if (allocated(InitialSpawnedSlots)) deallocate(InitialSpawnedSlots)

        ! Cleanup global storage
        call clean_global_det_data()

        ! Cleanup excitation generation storage
        call clean_excit_gen_store(fcimc_excit_gen_store)

        ! Cleanup cont time
        call clean_cont_time()

        ! Cleanup the load balancing
        call clean_load_balance()

        ! Cleanup adi caches
        call clean_adi()


        if (tGAS) then
            if (GAS_exc_gen == possible_GAS_exc_gen%DISCONNECTED) then
                call clearGAS()
            else if (GAS_exc_gen == possible_GAS_exc_gen%DISCARDING) then
                call finalize_GASCI_discarding()
            end if
        end if

        ! Cleanup excitation generator
        if (t_pcpp_excitgen) call finalize_pcpp_excitgen()
        if (t_pchb_excitgen) call finalize_pchb_excitgen()

        if (tSemiStochastic) call end_semistoch()

        if (tTrialWavefunction) call end_trial_wf()

!There seems to be some problems freeing the derived mpi type.
!        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!Free the mpi derived type that we have created for the hashes.
!            CALL MPI_Type_free(mpilongintegertype,error)
!            IF(error.ne.MPI_SUCCESS) THEN
!                CALL MPI_Error_string(error,message,length,temp)
!                IF(temp.ne.MPI_SUCCESS) THEN
!                    write(iout,*) "REALLY SERIOUS PROBLEMS HERE!",temp
!                    CALL neci_flush(iout)
!                end if
!                write(iout,*) message(1:length)
!                write(iout,*) "ERROR FOUND"
!                CALL neci_flush(iout)
!            end if
!        end if

    end subroutine DeallocFCIMCMemPar

    subroutine InitFCIMC_CSF()
        implicit none
        integer(n_int), allocatable ::  initSpace(:, :)
        integer :: count, nUp, nOpen
        integer :: i, j, lwork, proc
        integer :: DetHash, pos, TotWalkersTmp, nI(nel), nJ(nel)
        integer(n_int) :: ilutJ(0:NIfTot)
        integer(n_int), allocatable :: openSubspace(:)
        real(dp), allocatable :: S2(:, :), eigsImag(:), eigs(:), evs(:, :), void(:, :), work(:)
        real(dp) :: normalization, rawWeight, HDiag, tmpSgn(lenof_sign)
        integer :: err
        real(dp) :: HFWeight(inum_runs)
        logical :: tSuccess
        character(*), parameter :: t_r = "InitFCIMC_CSF"

        ! get the number of open orbitals
        nOpen = sum(nOccOrbs) - sum(nClosedOrbs)
        ! in a closed shell system, nothing to do
        if (nOpen == 0) then
            call InitFCIMC_HF()
            return
        end if
        ! first, set up the space considered for the CSF
        call generateInitSpace()
        if (allocated(openSubspace)) deallocate(openSubspace)
        ! we now have initSpace(:,:) with iluts belonging to all possible initial
        ! dets (i.e. all dets contributing to the target CSF) -> construct S2
        allocate(S2(count, count))
        do i = 1, count
            do j = 1, count
                S2(i, j) = S2Matel(initSpace(:, i), initSpace(:, j))
            end do
        end do

        ! prepare the diagonalization
        allocate(eigs(count))
        allocate(evs(count, count))
        allocate(eigsImag(count))
        allocate(work(1))
        allocate(void(0, 0))
        ! workspace query, get how much tmp memory we need
        call dgeev('N', 'V', count, S2, count, eigs, eigsImag, void, count, evs, count, work, -1, err)
        ! allocate work array
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        ! diagonalize S2
        call dgeev('N', 'V', count, S2, count, eigs, eigsImag, void, count, evs, count, work, lwork, err)
        deallocate(void)
        deallocate(work)
        deallocate(eigsImag)

        ! transfer the eigenvector
        do i = 1, count
            if (abs(S2Init * (S2Init + 1) - eigs(i)) < eps) exit
        end do
        if (i > count) then
            call stop_all(t_r, "Requested S2 eigenvalue does not exist")
        end if
        eigs = evs(:, i)
!      normalization = minval(eigs)
        rawWeight = sum(abs(eigs))
        normalization = InitialPart / rawWeight

!      refLoc = maxloc(eigs)
        eigs = eigs * normalization

        TotWalkers = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0
        do i = 1, count
            call decode_bit_det(nI, initSpace(:, i))
            proc = DetermineDetNode(nel, nI, 0)
            if (iProcIndex == proc) then
                HDiag = get_diagonal_matel(nI, initSpace(:, i))
                DetHash = FindWalkerHash(nI, size(HashIndex))
                TotWalkersTmp = int(TotWalkers)
                tmpSgn = eigs(i)
                call encode_sign(initSpace(:, i), tmpSgn)
                if (tHPHF) then
                    call FindDetSpinSym(nI, nJ, nel)
                    call encodebitdet(nJ, ilutJ)
                    ! if initSpace(:,i) is not the det of the HPHF pair we are storing,
                    ! skip this - the correct contribution will be stored once
                    ! the spin-flipped version is stored
                    if (DetBitLT(initSpace(:, i), ilutJ, NIfD) == 1) cycle
                end if
                call AddNewHashDet(TotWalkersTmp, initSpace(:, i), DetHash, nI, HDiag, pos, err)
                TotWalkers = TotWalkersTmp
            end if
            ! reset the reference?
        end do

        call hash_table_lookup(HFDet, ilutHF, nifd, HashIndex, CurrentDets, i, DetHash, tSuccess)
        if (tSuccess) then
            call extract_sign(CurrentDets(:, i), tmpSgn)
            do i = 1, inum_runs
                HFWeight(i) = mag_of_run(tmpSgn, i)
            end do
        else
            HFWeight = 0.0_dp
        end if

        AllTotParts = InitialPart
        AllTotPartsOld = InitialPart
        OldAllAvWalkersCyc = InitialPart
        OldAllHFCyc = HFWeight
        OldAllNoatHF = HFWeight

        ! cleanup
        deallocate(evs)
        deallocate(eigs)
        deallocate(S2)
        deallocate(initSpace)
    contains

        !------------------------------------------------------------------------------------------!

        subroutine generateOpenOrbIluts()
            use IntegralsData, only: nfrozen
            implicit none

            count = 0
            nUp = (nel + lms) / 2 - sum(nClosedOrbs) + nfrozen / 2
            do i = 1, 2**nOpen - 1
                if (popcnt(i) == nUp) then
                    count = count + 1
                    if (allocated(openSubspace)) openSubspace(count) = i
                end if
            end do
        end subroutine generateOpenOrbIluts

        !------------------------------------------------------------------------------------------!

        subroutine generateInitSpace()
            implicit none

            integer :: nI(nel), nIBase(nel)
            integer :: iEl, iElBase, iOpen
            integer :: openOrbList(nel)
            logical :: previousCont, nextCont, tClosed
            character(*), parameter :: t_r = "generateInitSpace"

            ! create a list of all open-shell determinants with the correct spin+orbs
            nUp = (nel + lms) / 2
            call generateOpenOrbIluts()
            allocate(openSubspace(count))
            call generateOpenOrbIluts()

            ! convert open-shell-only iluts into full iluts
            ! use the reference to determine which orbitals shall participate
            iElBase = 1
            nIBase = 0
            ! generate a list of open orbitals
            iOpen = 0
            openOrbList = 0
            do i = 1, nel
                if (i == 1) then
                    previousCont = .false.
                else
                    previousCont = FDet(i - 1) == FDet(i) - 1
                end if
                if (i == nel) then
                    nextCont = .false.
                else
                    nextCont = FDet(i + 1) == FDet(i) + 1
                end if
                tClosed = .true.
                ! identify open orbitals, using the ordering: does the next one belong
                ! to the same spatial orb?
                if (is_beta(FDet(i)) .and. .not. nextCont) then
                    iOpen = iOpen + 1
                    openOrbList(iOpen) = FDet(i)
                    tClosed = .false.
                end if
                ! or the previous one?
                if (is_alpha(FDet(i)) .and. .not. previousCont) then
                    iOpen = iOpen + 1
                    openOrbList(iOpen) = FDet(i)
                    tClosed = .false.
                end if
                ! if the orbital is not open, it is closed
                if (tClosed) then
                    nIBase(iElBase) = FDet(i)
                    iElBase = iElBase + 1
                end if
            end do
            if (iOpen /= nOpen) then
                write(iout, *) "nOpen/iOpen conflict", FDet, openOrbList, iOpen, nOpen
                call stop_all(t_r, "Error in determining open shell orbitals")
            end if

            allocate(initSpace(0:NIfTot, count))
            ! now, add the open-shell contribution
            do i = 1, count
                ! start from the closed-shell base
                nI = nIBase
                iEl = iElBase
                ! for each open orb, add the electron
                do j = 1, nOpen
                    ! based on openSubspace(i), we are considering alpha/beta for this
                    ! orb
                    if (btest(openSubspace(i), j - 1)) then
                        nI(iEl) = get_alpha(openOrbList(j))
                    else
                        nI(iEl) = get_beta(openOrbList(j))
                    end if
                    iEl = iEl + 1
                end do
                ! encode the determinant to the initial space
                call sort(nI)
                call EncodeBitDet(nI, initSpace(:, i))
            end do

        end subroutine generateInitSpace

        function S2Matel(ilutA, ilutB) result(matel)
            implicit none
            integer(n_int), intent(in) :: ilutA(0:NIfTot), ilutB(0:NIfTot)
            integer(n_int) :: splus(0:NIfTot), sminus(0:NIfTot)
            real(dp) :: matel
            integer :: k, m, nI(nel), upOrb, downOrb

            matel = 0.0_dp
            if (DetBitEq(ilutA, ilutB, NIfD)) then
                matel = matel + real(lms * (lms + 2), dp) / 4.0_dp
            end if
            ! get the offdiag part of S2: S-S+
            call decode_bit_det(nI, ilutA)
            do k = 1, nel
                ! first, apply S- to all electrons (additively)
                if (is_beta(nI(k))) then
                    sminus = ilutA
                    downOrb = get_alpha(nI(k))
                    clr_orb(sminus, nI(k))
                    set_orb(sminus, downOrb)
                    ! now, apply S+ to all electrons (again, sum)
                    do m = 1, nel
                        ! check if it yields 0
                        if (is_alpha(nI(m)) .or. m == k) then
                            splus = sminus
                            ! if not, go on
                            if (m == k) then
                                clr_orb(splus, downOrb)
                            else
                                clr_orb(splus, nI(m))
                            end if
                            upOrb = get_beta(nI(m))
                            set_orb(splus, upOrb)
                            if (DetBitEq(splus, ilutB, NIFD)) matel = matel + 1.0_dp
                        end if
                    end do
                end if
            end do
        end function S2Matel
    end subroutine InitFCIMC_CSF

    !------------------------------------------------------------------------------------------!

    subroutine InitFCIMC_HF()

        integer :: run, DetHash
        real(dp), dimension(lenof_sign) :: InitialSign
        real(dp) :: h_temp

        if (tOrthogonaliseReplicas) then
            call InitFCIMC_HF_orthog()
            return
        end if

        InitialPartVec = 0.0_dp
        do run = 1, inum_runs
            InitialPartVec(min_part_type(run)) = InitialPart
#ifdef CMPLX_
            InitialPartVec(max_part_type(run)) = 0.0_dp
#endif
        end do

        !Setup initial walker local variables for HF walkers start
        IF (iProcIndex == iRefProc(1)) THEN

            ! Encode the reference determinant identification.
            call encode_det(CurrentDets(:, 1), iLutHF)

            !Point at the correct position for the first walker
            DetHash = FindWalkerHash(HFDet, nWalkerHashes)    !Find det hash position
            HashIndex(DetHash)%Ind = 1

            ! Clear the flags
            call clear_all_flags(CurrentDets(:, 1))

            ! Set reference determinant as an initiator if
            ! tTruncInitiator is set, for both imaginary and real flags
            ! in real-time calculations, the reference does not have any special role
            if (tTruncInitiator) then
                do run = 1, inum_runs
                    call set_flag(CurrentDets(:, 1), get_initiator_flag_by_run(run))
                end do
            end if

            ! If running a semi-stochastic simulation, set flag to specify the Hartree-Fock is in the
            ! deterministic space.
            if (tSemiStochastic) then
                do run = 1, inum_runs
                    call set_flag(CurrentDets(:, 1), flag_deterministic(run))
                end do
            end if

            ! if no reference energy is used, explicitly get the HF energy
            if (tZeroRef) then
                if (tHPHF) then
                    h_temp = hphf_diag_helement(HFDet, ilutHF)
                else
                    h_temp = get_helement(HFDet, HFDet, 0)
                end if
            else
                ! HF energy is equal to 0 (when used as reference energy)
                h_temp = 0.0_dp
            end if
            call set_det_diagH(1, h_temp)
            HFInd = 1

            call store_decoding(1, HFDet)

            if (tContTimeFCIMC .and. tContTimeFull) &
                call set_spawn_rate(1, spawn_rate_full(HFDet, ilutHF))

            ! Obtain the initial sign
            InitialSign = 0.0_dp
            if (tStartSinglePart) then
                InitialSign(:) = InitialPartVec(:)
                TotParts(:) = InitialPartVec(:)
                TotPartsOld(:) = InitialPartVec(:)
            else
                do run = 1, inum_runs
                    InitialSign(min_part_type(run)) = InitWalkers
                    TotParts(min_part_type(run)) = real(InitWalkers, dp)
                    TotPartsOld(min_part_type(run)) = real(InitWalkers, dp)
#ifdef CMPLX_
                    TotParts(max_part_type(run)) = 0.0_dp
                    TotPartsOld(max_part_type(run)) = 0.0_dp
#endif
                end do
            end if

            ! set initial values for global control variables.

            TotWalkers = 1
            TotWalkersOld = 1
            NoatHF(:) = InitialSign(:)
            call encode_sign(CurrentDets(:, 1), InitialSign)
        ELSE
            NoatHF(:) = 0.0_dp
            TotWalkers = 0
            TotWalkersOld = 0
        end if

        OldAllNoatHF(:) = 0.0_dp
        AllNoatHF(:) = 0.0_dp
        IF (TStartSinglePart) THEN
            !Initialise global variables for calculation on the root node
            IF (iProcIndex == root) THEN
                OldAllNoatHF = InitialPartVec
                do run = 1, inum_runs
                    OldAllAvWalkersCyc(run) = sum(InitialPartVec( &
                                                  min_part_type(run):max_part_type(run)))
                end do
                AllNoatHF = InitialPartVec
                InstNoatHF = InitialPartVec
                AllTotParts = InitialPartVec
                AllTotPartsOld = InitialPartVec
                AllNoAbortedOld(:) = 0.0_dp
                iter_data_fciqmc%tot_parts_old = InitialPartVec
                AllTotWalkers = 1
                AllTotWalkersOld = 1
                do run = 1, inum_runs
                    OldAllHFCyc(run) = ARR_RE_OR_CPLX(InitialPartVec, run)
                end do
            end if
        ELSE
            !In this, only one processor has initial particles.
            IF (iProcIndex == Root) THEN
                AllTotWalkers = 1
                AllTotWalkersOld = 1
                do run = 1, inum_runs
                    iter_data_fciqmc%tot_parts_old(run) = real(InitWalkers, dp)
                    AllTotParts(run) = InitWalkers
                    AllTotPartsOld(run) = InitWalkers
                    AllNoAbortedOld(run) = 0.0_dp
                end do
            end if
        end if

    end subroutine InitFCIMC_HF

    subroutine InitFCIMC_HF_orthog()

        ! This is a reimplementation of InitFCIMC_HF to work with multiple
        ! different reference states.
        !
        ! In the end, we expect to be able to just substitute it back in
        ! for the original. It should give the same results
        ! TODO: Substitute back

        integer :: run, site, hash_val, i
        logical :: repeated
        HElement_t(dp) :: hdiag
        character(*), parameter :: this_routine = 'InitFCIMC_HF_orthog'

        ! Add some implementation guards

        ! Default values, unless overridder for individual procs
        NoatHF = 0.0_dp
        TotWalkers = 0
        TotWalkersOld = 0
        tRef_Not_HF = .true.
        tNoBrillouin = .true.

        !
        ! Initialise each of the runs separately.
        site = 0
        do run = 1, inum_runs

            ! If this run should have the reference on this site, then
            ! initialise it as appropriate.
            if (iProcIndex == iRefProc(run)) then

                ! Check if this reference is the same as any of the previous
                ! ones. If it is not, then at the end of the loop (i == site+1)
                repeated = .false.
                do i = 1, site
                    if (DetBitEQ(CurrentDets(:, i), ilutRef(:, run))) then
                        repeated = .true.
                        exit
                    end if
                end do
                site = i

                if (.not. repeated) then
                    ! Add the site to the main list (unless it is already there)
                    call encode_det(CurrentDets(:, site), ilutRef(:, run))
                    hash_val = FindWalkerHash(ProjEDet(:, run), nWalkerHashes)
                    call add_hash_table_entry(HashIndex, site, hash_val)

                    ! Clear all the flags and sign
                    call clear_all_flags(CurrentDets(:, site))
                    call nullify_ilut(CurrentDets(:, site))
                end if

                ! Set reference determinant as an initiator if tTruncInitiator
                if (tTruncInitiator) then
                    call set_flag(CurrentDets(:, site), get_initiator_flag_by_run(run))
                end if

                ! The global reference is the HF and is primary for printed
                ! energies.
                if (run == 1) HFInd = site
                if (tHPHF) then
                    hdiag = hphf_diag_helement(ProjEDet(:, run), ilutRef(:, run))
                else
                    hdiag = get_helement(ProjEDet(:, run), ProjEDet(:, run), 0)
                end if
                call set_det_diagH(site, real(hdiag, dp) - Hii)

                ! store the determinant
                call store_decoding(site, ProjEDet(:, run))

                ! Obtain the initial sign
                if (.not. tStartSinglePart) &
                    call stop_all(this_routine, "Only startsinglepart supported")
                call encode_part_sign(CurrentDets(:, site), InitialPart, min_part_type(run))

                ! Initial control values
                TotWalkers = site
                TotWalkersOld = site
                NoatHF(min_part_type(run)) = InitialPart
                TotParts(min_part_type(run)) = real(InitialPart, dp)
                TotPartsOld(min_part_type(run)) = real(InitialPart, dp)
            end if
        end do

        ! Check to ensure that the following code is valid
        if (.not. tStartSinglePart) &
            call stop_all(this_routine, "Only startsinglepart supported")

        ! Initialise global variabes for calculation on the root node
        OldAllNoatHF = 0.0_dp
        AllNoatHF = 0.0_dp
        call MPISum(TotWalkers, AllTotWalkers)
        if (iProcIndex == root) then
            OldAllNoatHF(:) = InitialPart
            OldAllAvWalkersCyc(:) = InitialPart
            AllNoatHF(:) = InitialPart
            InstNoatHF(:) = InitialPart
            AllTotParts(:) = InitialPart
            AllTotPartsOld(:) = InitialPart
            AllNoAbortedOld(:) = InitialPart
            OldAllHFCyc(:) = InitialPart

            TotWalkersOld = TotWalkers
        end if

    end subroutine InitFCIMC_HF_orthog

    subroutine InitFCIMC_trial()

        ! Use the code generated for the KPFCIQMC excited state calculations
        ! to initialise the FCIQMC simulation.

!         integer :: nexcit, ndets_this_proc, i, det(nel)
        integer :: nexcit, ndets_this_proc, det(nel)
        integer(int64) :: i

        type(basisfn) :: sym
        real(dp) :: evals(inum_runs / nreplicas)
        HElement_t(dp), allocatable :: evecs_this_proc(:, :)
        integer(MPIArg) :: space_sizes(0:nProcessors - 1), space_displs(0:nProcessors - 1)
        character(*), parameter :: this_routine = 'InitFCIMC_trial'

        nexcit = inum_runs / nreplicas

        ! Create the trial excited states
        if (allocated(trial_init_reorder)) then
            call calc_trial_states_lanczos(init_trial_in, nexcit, ndets_this_proc, &
                                           SpawnedParts, evecs_this_proc, evals, &
                                           space_sizes, space_displs, trial_init_reorder)
        else
            call calc_trial_states_lanczos(init_trial_in, nexcit, ndets_this_proc, &
                                           SpawnedParts, evecs_this_proc, evals, &
                                           space_sizes, space_displs)
        end if
        ! Determine the walker populations associated with these states
        call set_trial_populations(nexcit, ndets_this_proc, evecs_this_proc)
        ! Set the trial excited states as the FCIQMC wave functions
        call set_trial_states(ndets_this_proc, evecs_this_proc, SpawnedParts, &
                              .false., tPairedReplicas)

        deallocate(evecs_this_proc)

        if (tSetInitialRunRef) call set_initial_run_references()

        ! Add an initialisation check on symmetries.
        if ((.not. tHub) .and. (.not. tUEG)) then
            do i = 1, TotWalkers
                call decode_bit_det(det, CurrentDets(:, i))
                call getsym_wrapper(det, sym)
                if (sym%sym%S /= HFSym%sym%S .or. sym%ml /= HFSym%Ml) &
                    call stop_all(this_routine, "Invalid det found")
            end do
        end if

    end subroutine InitFCIMC_trial

    subroutine set_initial_run_references()

        ! Analyse each of the runs, and set the reference determinant to the
        ! det with the largest coefficient, rather than the currently guessed
        ! one...

        HElement_t(dp) :: largest_coeff, sgn
        integer(n_int) :: largest_det(0:NIfTot)
!         integer :: run, j
        integer(int64) :: j
        integer :: run
        integer(int32) :: proc_highest
        integer(n_int) :: ilut(0:NIfTot)
        integer(int32) :: int_tmp(2)
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'set_initial_run_references'
#endif

        do run = 1, inum_runs

            if (tMultipleInitialRefs) then
                ! Use user specified reference states.
                call EncodeBitDet(initial_refs(:, run), ilut)
                call update_run_reference(ilut, run)
            else
                ! Find the largest det on this processor
                largest_coeff = h_cast(0.0_dp)
                do j = 1, TotWalkers
                    sgn = extract_run_sign(CurrentDets(:, j), run)
                    if (abs(sgn) > abs(largest_coeff)) then
                        largest_coeff = sgn
                        largest_det = CurrentDets(:, j)
                    end if
                end do

                ! Find the largest det on any processor (n.b. discard the
                ! non-integer part. This isn't all that important).
                ! [W.D. 15.5.2017:]
                ! for the test suite problems, maybe it is important..
                ! because there seems to be some compiler dependent
                ! differences..
                call MPIAllReduceDatatype( &
                    (/int(abs(largest_coeff), int32), int(iProcIndex, int32)/), 1, &
                    MPI_MAXLOC, MPI_2INTEGER, int_tmp)
                proc_highest = int_tmp(2)
                call MPIBCast(largest_det, NIfTot + 1, int(proc_highest, sizeof_int))
!                 call MPIBCast(largest_det, NIfTot+1, proc_highest)

                write(6, *) 'Setting ref', run
                call writebitdet(6, largest_det, .true.)

                ! Set this det as the reference
                call update_run_reference(largest_det, run)

            end if
        end do

        if (tMultiRefShift) then
            do run = 1, inum_runs
                DiagSft(run) = proje_ref_energy_offsets(run)
            end do
        end if

    end subroutine set_initial_run_references

    subroutine InitFCIMC_CAS()

        ! Routine to initialise the particle distribution according to a CAS diagonalisation.
        ! This hopefully will help with close-lying excited states of the same sym.

        type(BasisFN) :: CASSym
        integer :: i, ierr, nEval, NKRY1, NBLOCK, LSCR, LISCR, DetIndex
        integer :: iNode, nBlocks, nBlockStarts(2), DetHash
        integer :: CASSpinBasisSize, nCASDet, ICMax, GC, LenHamil, iInit
        integer :: nHPHFCAS, iCasDet, ExcitLevel
        real(dp) :: NoWalkers
        integer, allocatable :: CASBrr(:), CASDet(:), CASFullDets(:, :), nRow(:), Lab(:), ISCR(:), INDEX(:)
        integer, pointer :: CASDetList(:, :) => null()
        integer(n_int) :: iLutnJ(0:NIfTot)
        logical :: tMC, tHPHF_temp, tHPHFInts_temp
        HElement_t(dp) :: HDiagTemp
        real(dp), allocatable :: CK(:, :), W(:), CKN(:, :), Hamil(:), A_Arr(:, :), V(:), BM(:), T(:), WT(:)
        real(dp), allocatable :: SCR(:), WH(:), Work2(:), V2(:, :), AM(:)
        real(dp), allocatable :: Work(:)
        integer(TagIntType) :: ATag = 0, VTag = 0, BMTag = 0, TTag = 0, WTTag = 0, SCRTag = 0, WHTag = 0, Work2Tag = 0, V2Tag = 0
        integer(TagIntType) :: ISCRTag = 0, IndexTag = 0, AMTag = 0
        integer(TagIntType) :: WorkTag = 0
        real(dp) :: CASRefEnergy, TotWeight, PartFac, amp, rat, r, GetHElement
        real(dp), dimension(lenof_sign) :: temp_sign
        real(dp) :: energytmp(nel), max_wt
        integer  :: tmp_det(nel), det_max, run
        type(ll_node), pointer :: TempNode
        character(len=*), parameter :: this_routine = 'InitFCIMC_CAS'
#ifdef CMPLX_
        call stop_all(this_routine, "StartCAS currently does not work with complex walkers")
#endif
        if (tReadPops) call stop_all(this_routine, "StartCAS cannot work with with ReadPops")
        if (tStartSinglePart) call stop_all(this_routine, "StartCAS cannot work with StartSinglePart")
        if (tRestartHighPop) call stop_all(this_routine, "StartCAS cannot with with dynamically restarting calculations")

        write(iout, *) "Initialising walkers proportional to a CAS diagonalisation..."
        write(iout, '(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (" &
            , (OccCASOrbs + VirtCASOrbs) / 2, ",", OccCASOrbs, ")"
        do I = NEl - OccCASorbs + 1, NEl
            write(iout, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(iout, G1(BRR(I))%SYM, .FALSE.)
            write(iout, '(I4)', advance='no') G1(BRR(I))%Ml
            write(iout, '(2F19.9)') ARR(I, 1), ARR(BRR(I), 2)
        end do
        write(iout, '(A)') " ================================================================================================="
        do I = NEl + 1, NEl + VirtCASOrbs
            write(iout, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(iout, G1(BRR(I))%SYM, .FALSE.)
            write(iout, '(I4)', advance='no') G1(BRR(I))%Ml
            write(iout, '(2F19.9)') ARR(I, 1), ARR(BRR(I), 2)
        end do

        CASSpinBasisSize = OccCASorbs + VirtCASorbs
        allocate(CASBrr(1:CASSpinBasisSize))
        allocate(CASDet(1:OccCasOrbs))
        do i = 1, CASSpinBasisSize
            !Run through the cas space, and create an array which will map these orbtials to the
            !orbitals they actually represent.
            CASBrr(i) = BRR(i + (NEl - OccCasorbs))
        end do

        !Calculate symmetry of CAS determinants, and check that this will be the same as the reference determinant
        !for the rest of the FCIMC calculations.
        CASDet = CasBRR(1:OccCasOrbs)
        call sort(CasDet)

        write(iout, *) "CAS Det is: "
        call write_det_len(iout, CASDet, OccCASOrbs, .true.)
        call GetSym(CASDet, OccCASOrbs, G1, nBasisMax, CASSym)
        write(iout, *) "Spatial symmetry of CAS determinants: ", CASSym%Sym%S
        write(iout, *) "Ms of CAS determinants: ", CASSym%Ms
        if (tFixLz) then
            write(iout, *) "Ml of CAS determinants: ", CASSym%Ml
        end if
        call neci_flush(iout)

        if (CASSym%Ml /= LzTot) call stop_all(this_routine, "Ml of CAS ref det does not match Ml of full reference det")
        if (CASSym%Ms /= 0) call stop_all(this_routine, "CAS diagonalisation can only work with closed shell CAS spaces initially")
        if (CASSym%Sym%S /= HFSym%Sym%S) then
            call stop_all(this_routine, "Sym of CAS ref det does not match Sym of full reference det")
        end if

        !First, we need to generate all the excitations.
        call gndts(OccCASorbs, CASSpinBasisSize, CASBrr, nBasisMax, CASDetList, .true., G1, tSpn, LMS, .true., CASSym, nCASDet, iCASDet)

        if (nCASDet == 0) call stop_all(this_routine, "No CAS determinants found.")
        write(iout, *) "Number of symmetry allowed CAS determinants found to be: ", nCASDet
        allocate(CASDetList(OccCASorbs, nCASDet), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating CASDetList")
        CASDetList(:, :) = 0

        !Now fill up CASDetList...
        call gndts(OccCASorbs, CASSpinBasisSize, CASBrr, nBasisMax, &
                   CASDetList, .false., G1, tSpn, LMS, .true., CASSym, &
                   nCASDet, iCASDet)

        !We have a complication here. If we calculate the hamiltonian from these CAS determinants, then we are not
        !including the mean-field generated from the other occupied orbitals. We need to either 'freeze' the occupied
        !orbitals and modify the 1 & two electron integrals, or add the other electrons back into the list. We do the latter.
        allocate(CASFullDets(NEl, nCASDet), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating CASFullDets")
        CASFullDets(:, :) = 0

        ! Get the first part of a determinant with the lowest energy, rather
        ! than lowest index number orbitals
        energytmp = ARR(ProjEDet(:, 1), 2)
        tmp_det = ProjEDet(:, 1)
        call sort(energytmp, tmp_det)

        ! Construct the determinants resulting from the CAS expansion.
        do i = 1, nCASDet
            CASFullDets(1:nel - OccCASorbs, i) = tmp_det(1:nel - OccCASOrbs)
            CASFullDets(nel - OccCASorbs + 1:nel, i) = CASDetList(1:OccCASorbs, i)
            call sort(CASFullDets(:, i))
        end do
        deallocate(CASDetList)

        write(iout, *) "First CAS determinant in list is: "
        call write_det(iout, CASFullDets(:, 1), .true.)

        if (nCASDet > 1300) then
            !Do lanczos
            nEval = 4
        else
            nEval = nCASDet
        end if
        write(iout, "(A,I4,A)") "Calculating lowest ", nEval, " eigenstates of CAS Hamiltonian..."
        allocate(Ck(nCASDet, nEval), stat=ierr)
        Ck = 0.0_dp
        allocate(W(nEval), stat=ierr)    !Eigenvalues
        W = 0.0_dp
        if (ierr /= 0) call stop_all(this_routine, "Error allocating")

        write(iout, *) "Calculating hamiltonian..."
        allocate(nRow(nCASDet), stat=ierr)
        nRow = 0
        ICMax = 1
        tMC = .false.

        !HACK ALERT!! Need to fill up array in space of determinants, not HPHF functions.
        !Turn off tHPHFInts and tHPHF and turn back on after the hamiltonian constructed.
        tHPHF_temp = tHPHF
        tHPHFInts_temp = tHPHFInts
        tHPHF = .false.
        tHPHFInts = .false.

        ! do not pass an unallocated array, so dummy-allocate
        allocate(Hamil(0))
        allocate(Lab(0))
        CALL Detham(nCASDet, NEl, CASFullDets, Hamil, Lab, nRow, .true., ICMax, GC, tMC)
        deallocate(Lab)
        deallocate(Hamil)
        LenHamil = GC
        write(iout, *) "Allocating memory for hamiltonian: ", LenHamil * 2
        allocate(Hamil(LenHamil), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating Hamil")
        Hamil = 0.0_dp
        allocate(Lab(LenHamil), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating Lab")
        Lab = 0
        call Detham(nCASDet, NEl, CASFullDets, Hamil, Lab, nRow, .false., ICMax, GC, tMC)

        CASRefEnergy = GETHELEMENT(1, 1, HAMIL, LAB, NROW, NCASDET)
        write(iout, *) "Energy of first CAS det is: ", CASRefEnergy

        ! Turn back on HPHFs if needed.
        tHPHF = tHPHF_temp
        tHPHFInts = tHPHFInts_temp

!        if(abs(CASRefEnergy-Hii).gt.1.0e-7_dp) then
!            call stop_all(this_routine,"CAS reference energy does not match reference energy of full space")
!        end if

        if (nCASDet > 1300) then
            !Lanczos
            NKRY1 = NKRY + 1
            NBLOCK = MIN(NEVAL, NBLK)
            LSCR = MAX(nCASDet * NEVAL, 8 * NBLOCK * NKRY)
            LISCR = 6 * NBLOCK * NKRY
            allocate(A_Arr(NEVAL, NEVAL), stat=ierr)
            CALL LogMemAlloc('A_Arr', NEVAL**2, 8, this_routine, ATag, ierr)
            A_Arr = 0.0_dp
            allocate(V(nCASDet * NBLOCK * NKRY1), stat=ierr)
            CALL LogMemAlloc('V', nCASDet * NBLOCK * NKRY1, 8, this_routine, VTag, ierr)
            V = 0.0_dp
            allocate(AM(NBLOCK * NBLOCK * NKRY1), stat=ierr)
            CALL LogMemAlloc('AM', NBLOCK * NBLOCK * NKRY1, 8, this_routine, AMTag, ierr)
            AM = 0.0_dp
            allocate(BM(NBLOCK * NBLOCK * NKRY), stat=ierr)
            CALL LogMemAlloc('BM', NBLOCK * NBLOCK * NKRY, 8, this_routine, BMTag, ierr)
            BM = 0.0_dp
            allocate(T(3 * NBLOCK * NKRY * NBLOCK * NKRY), stat=ierr)
            CALL LogMemAlloc('T', 3 * NBLOCK * NKRY * NBLOCK * NKRY, 8, this_routine, TTag, ierr)
            T = 0.0_dp
            allocate(WT(NBLOCK * NKRY), stat=ierr)
            CALL LogMemAlloc('WT', NBLOCK * NKRY, 8, this_routine, WTTag, ierr)
            WT = 0.0_dp
            allocate(SCR(LScr), stat=ierr)
            CALL LogMemAlloc('SCR', LScr, 8, this_routine, SCRTag, ierr)
            SCR = 0.0_dp
            allocate(ISCR(LIScr), stat=ierr)
            CALL LogMemAlloc('IScr', LIScr, 4, this_routine, IScrTag, ierr)
            ISCR(1:LISCR) = 0
            allocate(INDEX(NEVAL), stat=ierr)
            CALL LogMemAlloc('INDEX', NEVAL, 4, this_routine, INDEXTag, ierr)
            INDEX(1:NEVAL) = 0
            allocate(WH(nCASDet), stat=ierr)
            CALL LogMemAlloc('WH', nCASDet, 8, this_routine, WHTag, ierr)
            WH = 0.0_dp
            allocate(WORK2(3 * nCASDet), stat=ierr)
            CALL LogMemAlloc('WORK2', 3 * nCASDet, 8, this_routine, WORK2Tag, ierr)
            WORK2 = 0.0_dp
            allocate(V2(nCASDet, NEVAL), stat=ierr)
            CALL LogMemAlloc('V2', nCASDet * NEVAL, 8, this_routine, V2Tag, ierr)
            V2 = 0.0_dp
            allocate(CkN(nCASDet, nEval), stat=ierr)
            CkN = 0.0_dp
            !C..Lanczos iterative diagonalising routine
            if (t_non_hermitian) then
                call stop_all(this_routine, &
                              "NECI_FRSBLKH not adapted for non-hermitian Hamiltonians!")
            end if
         CALL NECI_FRSBLKH(nCASDet, ICMAX, NEVAL, HAMIL, LAB, CK, CKN, NKRY, NKRY1, NBLOCK, NROW, LSCR, LISCR, A_Arr, W, V, AM, BM, T, WT, &
         &  SCR, ISCR, INDEX, NCYCLE, B2L, .false., .false., .false., .true.)
            !Multiply all eigenvalues by -1.
            CALL DSCAL(NEVAL, -1.0_dp, W, 1)
            if (CK(1, 1) < 0.0_dp) then
                do i = 1, nCASDet
                    CK(i, 1) = -CK(i, 1)
                end do
            end if

            deallocate(CKN, A_Arr, V, BM, T, WT, SCR, WH, V2, iscr, index, AM)
            call logmemdealloc(this_routine, ATag)
            call logmemdealloc(this_routine, VTag)
            call logmemdealloc(this_routine, BMTag)
            call logmemdealloc(this_routine, TTag)
            call logmemdealloc(this_routine, WTTag)
            call logmemdealloc(this_routine, SCRTag)
            call logmemdealloc(this_routine, WHTag)
            call logmemdealloc(this_routine, V2Tag)
            call logmemdealloc(this_routine, iscrTag)
            call logmemdealloc(this_routine, indexTag)
            call logmemdealloc(this_routine, AMTag)
        else
            !complete diagonalisation
            allocate(Work(4 * nCASDet), stat=ierr)
            call LogMemAlloc('Work', 4 * nCASDet, 8, this_routine, WorkTag, ierr)
            allocate(Work2(3 * nCASDet), stat=ierr)
            call logMemAlloc('Work2', 3 * nCASDet, 8, this_routine, Work2Tag, ierr)
            nBlockStarts(1) = 1
            nBlockStarts(2) = nCASDet + 1
            nBlocks = 1
            if (t_non_hermitian) then
                call stop_all(this_routine, &
                              "HDIAG_neci is not set up for non-hermitian Hamiltonians!")
            end if
            call HDIAG_neci(nCASDet, Hamil, Lab, nRow, CK, W, Work2, Work, nBlockStarts, nBlocks)
            deallocate(Work)
            call LogMemDealloc(this_routine, WorkTag)
        end if
        !Deallocate all the lanczos arrays now.
        deallocate(nrow, lab, work2)
        call logmemdealloc(this_routine, Work2Tag)

        write(iout, *) "Diagonalisation complete. Lowest energy CAS eigenvalues/corr E are: "
        do i = 1, min(NEval, 10)
            write(iout, *) i, W(i), W(i) - CASRefEnergy
        end do

        TotWeight = 0.0_dp
        nHPHFCAS = 0
        max_wt = 0
        det_max = 0
        do i = 1, nCASDet
            if (tHPHF) then
                !Only allow valid HPHF functions
                call EncodeBitDet(CASFullDets(:, i), iLutnJ)
                if (IsAllowedHPHF(iLutnJ)) then
                    nHPHFCAS = nHPHFCAS + 1
                    if (.not. TestClosedShellDet(iLutnJ)) then
                        !Open shell. Weight is sqrt(2) of det weight.
                        TotWeight = TotWeight + (abs(CK(i, 1)) * sqrt(2.0_dp))
                        !Return this new weight to the CK array, so that we do not need to do this a second time.
                        CK(i, 1) = CK(i, 1) * sqrt(2.0_dp)
                    else
                        !Closed Shell
                        TotWeight = TotWeight + abs(CK(i, 1))
                    end if
                end if
            else
                TotWeight = TotWeight + abs(CK(i, 1))
            end if

            ! Find the maximum weighted determinant.
            if (abs(ck(i, 1)) > max_wt) then
                max_wt = abs(ck(i, 1))
                det_max = i
            end if
        end do

        ! Output details
        write(iout, *) "Total weight of lowest eigenfunction: ", TotWeight
        write(iout, *) "Maximum weighted det: ", det_max, max_wt

        !If the reference det is not the maximum weighted det, suggest that
        !we change it!
        if (.not. all(CASFullDets(:, det_max) == ProjEDet(:, 1))) then
            write(iout, *) 'The specified reference determinant is not the &
                       &maximum weighted determinant in the CAS expansion'
            write(iout, *) 'Use following det as reference:'
            call write_det(6, CASFullDets(:, det_max), .true.)
            call warning_neci(this_routine, "Poor reference chosen")
        end if

        if (tHPHF) write(iout, *) "Converting into HPHF space. Total HPHF CAS functions: ", nHPHFCAS

        if ((InitialPart.isclose.1._dp) .or. (InitialPart >= (InitWalkers * nNodes) - 50)) then
            !Here, all the walkers will be assigned to the CAS wavefunction.
            !InitialPart = 1 by default
            write(iout, "(A)") "All walkers specified in input will be distributed according to the CAS wavefunction."
            write(iout, "(A)") "Shift will be allowed to vary from the beginning"
            write(iout, "(A,F20.9)") "Setting initial shift to equal CAS correlation energy", W(1) - CASRefEnergy
            DiagSft = W(1) - CASRefEnergy
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac = (real(InitWalkers, dp) * real(nNodes, dp)) / TotWeight
        else
            !Here, not all walkers allowed will be initialised to the CAS wavefunction.
            write(iout, "(A,I15,A)") "Initialising ", int(InitialPart), " walkers according to the CAS distribution."
            write(iout, "(A,I15)") "Shift will remain fixed until the walker population reaches ", int(InitWalkers * nNodes)
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac = real(InitialPart, dp) / TotWeight
            tSinglePartPhase(:) = .true.
        end if

        !Now generate all excitations again, creating the required number of walkers on each one.
        DetIndex = 1
        NoatHF(:) = 0.0
        TotParts(:) = 0.0
        do i = 1, nCASDet
            if (tHPHF) then
                call EncodeBitDet(CASFullDets(:, i), iLutnJ)
                if (.not. IsAllowedHPHF(iLutnJ)) cycle
            end if
            iNode = DetermineDetNode(nel, CASFullDets(:, i), 0)
            if (iProcIndex == iNode) then
                !Number parts on this det = PartFac*Amplitude
                amp = CK(i, 1) * PartFac

                if (tRealCoeffByExcitLevel) ExcitLevel = FindBitExcitLevel(iLutnJ, iLutRef, nEl)
                if (tAllRealCoeff .or. &
                    & (tRealCoeffByExcitLevel .and. (ExcitLevel <= RealCoeffExcitThresh))) then
                    NoWalkers = amp
                else
                    NoWalkers = int(amp)
                    rat = amp - real(NoWalkers, dp)
                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) then
                        if (amp < 0.0_dp) then
                            NoWalkers = NoWalkers - 1
                        else
                            NoWalkers = NoWalkers + 1
                        end if
                    end if
                end if

                if (abs(NoWalkers) > 1.0e-12_dp) then
                    call EncodeBitDet(CASFullDets(:, i), iLutnJ)
                    if (DetBitEQ(iLutnJ, iLutRef(:, 1), nifd)) then
                        !Check if this determinant is reference determinant, so we can count number on hf.
                        do run = 1, inum_runs
                            NoatHF(run) = NoWalkers
                        end do
                    end if
                    call encode_det(CurrentDets(:, DetIndex), iLutnJ)
                    call clear_all_flags(CurrentDets(:, DetIndex))
                    do run = 1, inum_runs
                        temp_sign(run) = NoWalkers
                    end do
                    call encode_sign(CurrentDets(:, DetIndex), temp_sign)

                    ! Store the diagonal matrix elements
                    if (tHPHF) then
                        HDiagTemp = hphf_diag_helement(CASFullDets(:, i), iLutnJ)
                    else
                        HDiagTemp = get_helement(CASFullDets(:, i), CASFullDets(:, i), 0)
                    end if
                    call set_det_diagH(DetIndex, real(HDiagTemp, dp) - Hii)
                    call store_decoding(DetIndex, CASFullDets(:, i))

                    if (tTruncInitiator) then
                        !Set initiator flag if needed (always for HF)
                        call CalcParentFlag(DetIndex, iInit)
                    end if

                    DetHash = FindWalkerHash(CASFullDets(:, i), nWalkerHashes)
                    TempNode => HashIndex(DetHash)
                    ! If the first element in the list has not been used.
                    if (TempNode%Ind == 0) then
                        TempNode%Ind = DetIndex
                    else
                        do while (associated(TempNode%Next))
                            TempNode => TempNode%Next
                        end do
                        allocate(TempNode%Next)
                        nullify (TempNode%Next%Next)
                        TempNode%Next%Ind = DetIndex
                    end if
                    nullify (TempNode)

                    DetIndex = DetIndex + 1
                    do run = 1, inum_runs
                        TotParts(run) = TotParts(run) + abs(NoWalkers)
                    end do
                end if
            end if   !End if desired node
        end do

        TotWalkers = DetIndex - 1   !This is the number of occupied determinants on each node
        TotWalkersOld = TotWalkers

        !Set local&global variables
        TotPartsOld = TotParts
        call mpisumall(TotParts, AllTotParts)
        call mpisumall(NoatHF, AllNoatHF)
        call mpisumall(TotWalkers, AllTotWalkers)
        OldAllNoatHF = AllNoatHF
        do run = 1, inum_runs
            OldAllHFCyc(run) = AllNoatHF(run)
            OldAllAvWalkersCyc(run) = AllTotParts(run)
        end do
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld = 0.0_dp

        deallocate(CK, W, Hamil, CASBrr, CASDet, CASFullDets)

    end subroutine InitFCIMC_CAS

    !Routine to initialise the particle distribution according to the MP1 wavefunction.
    !This hopefully will help with close-lying excited states of the same sym.
    subroutine InitFCIMC_MP1()
        real(dp) :: TotMP1Weight, amp, MP2Energy, PartFac, rat, r, energy_contrib
        HElement_t(dp) :: HDiagtemp
        integer :: iExcits, exflag, Ex(2, maxExcit), nJ(NEl), DetIndex, iNode
        integer :: iInit, DetHash, ExcitLevel, run, part_type
        integer(n_int) :: iLutnJ(0:NIfTot)
        real(dp) :: NoWalkers, temp_sign(lenof_sign)
        logical :: tAllExcitsFound, tParity
        type(ll_node), pointer :: TempNode
        character(len=*), parameter :: this_routine = "InitFCIMC_MP1"
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:, :)
        integer :: i

#ifdef CMPLX_
        call stop_all(this_routine, "StartMP1 currently does not work with complex walkers")
#endif
        if (tReadPops) call stop_all(this_routine, "StartMP1 cannot work with with ReadPops")
        if (tStartSinglePart) call stop_all(this_routine, "StartMP1 cannot work with StartSinglePart")
        if (tRestartHighPop) call stop_all(this_routine, "StartMP1 cannot with with dynamically restarting calculations")

        write(iout, *) "Initialising walkers proportional to the MP1 amplitudes..."

        if (tHPHF) then
            if (.not. TestClosedShellDet(iLutHF)) then
                call stop_all(this_routine, "Cannot use HPHF with StartMP1 if your reference is open-shell")
            end if
        end if

        if (tUEG) then
            !Parallel N^2M implementation of MP2 for UEG
            call CalcUEGMP2()
        end if

        !First, calculate the total weight - TotMP1Weight
        mp2energy = 0.0_dp
        TotMP1Weight = 1.0_dp
        iExcits = 0
        tAllExcitsFound = .false.
        if (tUEG) then
            exflag = 2
        else
            exflag = 3
        end if
        Ex(:, :) = 0
        ! for GUGA, this whole excitation generation is different!
        if (tGUGA) then
            ! should finally do some general routine, which does all this
            ! below...
            call convert_ilut_toGUGA(ilutHF, ilutG)
            call actHamiltonian(ilutG, excitations, iExcits)
            do i = 1, iExcits
                call convert_ilut_toNECI(excitations(:, i), ilutnJ)
                call decode_bit_det(nJ, iLutnJ)
                call return_mp1_amp_and_mp2_energy(nJ, iLutnJ, Ex, tParity, amp, energy_contrib)
                TotMP1Weight = TotMP1Weight + abs(amp)
                MP2Energy = MP2Energy + energy_contrib
            end do

        else
            do while (.true.)
                call GenExcitations3(HFDet, iLutHF, nJ, exflag, Ex, tParity, tAllExcitsFound, .false.)
                if (tAllExcitsFound) exit !All excits found

                call EncodeBitDet(nJ, iLutnJ)
                if (tHPHF) then
                    !Working in HPHF Space. Check whether determinant generated is an 'HPHF'
                    if (.not. IsAllowedHPHF(iLutnJ)) cycle
                end if
                iExcits = iExcits + 1

                call return_mp1_amp_and_mp2_energy(nJ, iLutnJ, Ex, tParity, amp, energy_contrib)
                TotMP1Weight = TotMP1Weight + abs(amp)
                MP2Energy = MP2Energy + energy_contrib
            end do
        end if

        if ((.not. tHPHF .and. .not. tGUGA) .and. (iExcits /= (nDoubles + nSingles))) then
            write(iout, *) nDoubles, nSingles, iExcits
            call stop_all(this_routine, "Not all excitations accounted for in StartMP1")
        end if

        write(iout, "(A,2G25.15)") "MP2 energy calculated: ", MP2Energy, MP2Energy + Hii

        if ((InitialPart.isclose.1._dp) .or. (InitialPart >= (InitWalkers * nNodes) - 50)) then
            !Here, all the walkers will be assigned to the MP1 wavefunction.
            !InitialPart = 1 by default
            write(iout, "(A)") "All walkers specified in input will be distributed according to the MP1 wavefunction."
            write(iout, "(A)") "Shift will be allowed to vary from the beginning"
            write(iout, "(A)") "Setting initial shift to equal MP2 correlation energy"
            DiagSft = MP2Energy
            !PartFac is the number of walkers that should reside on the HF determinant
            !in an intermediate normalised MP1 wavefunction.
            PartFac = (real(InitWalkers, dp) * real(nNodes, dp)) / TotMP1Weight
        else
            !Here, not all walkers allowed will be initialised to the MP1 wavefunction.
            write(iout, "(A,G15.5,A)") "Initialising ", InitialPart, " walkers according to the MP1 distribution."
            write(iout, "(A,G15.5)") "Shift will remain fixed until the walker population reaches ", InitWalkers * nNodes
            !PartFac is the number of walkers that should reside on the HF determinant
            !in an intermediate normalised MP1 wavefunction.
            PartFac = real(InitialPart, dp) / TotMP1Weight
            tSinglePartPhase(:) = .true.
        end if

        !Now generate all excitations again, creating the required number of walkers on each one.
        ! puh... for GUGA this gets messy to change below.. damn

        DetIndex = 1
        TotParts = 0.0
        tAllExcitsFound = .false.
        if (tUEG) then
            exflag = 2
        else
            exflag = 3
        end if
        Ex(:, :) = 0
        ! figure out if the HF det gets store in the excitation list too?
        ! if yes I have to modify that all a bit, and maybe also in
        ! other parts of the NECI code ... todo
        if (tGUGA) call stop_all(this_routine, "deprecated option with GUGA!")

        do while (.true.)
            call GenExcitations3(HFDet, iLutHF, nJ, exflag, Ex, tParity, tAllExcitsFound, .false.)
            if (tAllExcitsFound) exit !All excits found
            call EncodeBitDet(nJ, iLutnJ)
            if (tHPHF) then
                !Working in HPHF Space. Check whether determinant generated is an 'HPHF'
                if (.not. IsAllowedHPHF(iLutnJ)) cycle
            end if

            iNode = DetermineDetNode(nel, nJ, 0)
            if (iProcIndex == iNode) then
                call return_mp1_amp_and_mp2_energy(nJ, iLutnJ, Ex, tParity, amp, energy_contrib)
                amp = amp * PartFac

                if (tRealCoeffByExcitLevel) ExcitLevel = FindBitExcitLevel(iLutnJ, iLutRef(:, 1), nEl)
                if (tAllRealCoeff .or. &
                    & (tRealCoeffByExcitLevel .and. (ExcitLevel <= RealCoeffExcitThresh))) then
                    NoWalkers = amp
                else
                    NoWalkers = int(amp)
                    rat = amp - real(NoWalkers, dp)
                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) then
                        if (amp < 0.0_dp) then
                            NoWalkers = NoWalkers - 1
                        else
                            NoWalkers = NoWalkers + 1
                        end if
                    end if
                end if

                if (abs(NoWalkers) > 1.0e-12_dp) then
                    call encode_det(CurrentDets(:, DetIndex), iLutnJ)
                    call clear_all_flags(CurrentDets(:, DetIndex))
                    do run = 1, inum_runs
                        temp_sign(run) = NoWalkers
                    end do
                    call encode_sign(CurrentDets(:, DetIndex), temp_sign)

                    ! Store the diagonal matrix elements
                    if (tHPHF) then
                        HDiagTemp = hphf_diag_helement(nJ, iLutnJ)
                    else
                        HDiagTemp = get_helement(nJ, nJ, 0)
                    end if
                    call set_det_diagH(DetIndex, real(HDiagtemp, dp) - Hii)
                    ! store the determinant
                    call store_decoding(DetIndex, nJ)

                    if (tTruncInitiator) then
                        !Set initiator flag if needed (always for HF)
                        call CalcParentFlag(DetIndex, iInit)
                    end if

                    DetHash = FindWalkerHash(nJ, nWalkerHashes)
                    TempNode => HashIndex(DetHash)
                    ! If the first element in the list has not been used.
                    if (TempNode%Ind == 0) then
                        TempNode%Ind = DetIndex
                    else
                        do while (associated(TempNode%Next))
                            TempNode => TempNode%Next
                        end do
                        allocate(TempNode%Next)
                        nullify (TempNode%Next%Next)
                        TempNode%Next%Ind = DetIndex
                    end if
                    nullify (TempNode)

                    DetIndex = DetIndex + 1
                    do part_type = 1, lenof_sign
                        TotParts(part_type) = TotParts(part_type) + abs(NoWalkers)
                    end do
                end if
            end if   !End if desired node

        end do

        !Now for the walkers on the HF det
        if (iRefProc(1) == iProcIndex) then
            ! dont have to change this below, should also work for the GUGAc
            if (tAllRealCoeff .or. tRealCoeffByExcitLevel) then
                NoWalkers = PartFac
            else
                NoWalkers = int(PartFac)
                rat = PartFac - real(NoWalkers, dp)
                if (rat < 0.0_dp) &
                    call stop_all(this_routine, "Should not have negative weight on HF")
                r = genrand_real2_dSFMT()
                if (abs(rat) > r) NoWalkers = NoWalkers + 1
            end if
            if (abs(NoWalkers) > 1.0e-12_dp) then
                call encode_det(CurrentDets(:, DetIndex), iLutHF)
                call clear_all_flags(CurrentDets(:, DetIndex))
                do run = 1, inum_runs
                    temp_sign(run) = NoWalkers
                end do
                call encode_sign(CurrentDets(:, DetIndex), temp_sign)
                if (tTruncInitiator) then
                    !Set initiator flag (always for HF)
                    do run = 1, inum_runs
                        call set_flag(CurrentDets(:, DetIndex), get_initiator_flag(run))
                    end do
                end if
                call set_det_diagH(DetIndex, 0.0_dp)

                ! store the determinant
                call store_decoding(DetIndex, HFDet)

                ! Now add the Hartree-Fock determinant (not with index 1).
                DetHash = FindWalkerHash(HFDet, nWalkerHashes)
                TempNode => HashIndex(DetHash)
                ! If the first element in the list has not been used.
                if (TempNode%Ind == 0) then
                    TempNode%Ind = DetIndex
                else
                    do while (associated(TempNode%Next))
                        TempNode => TempNode%Next
                    end do
                    allocate(TempNode%Next)
                    nullify (TempNode%Next%Next)
                    TempNode%Next%Ind = DetIndex
                end if
                nullify (TempNode)

                DetIndex = DetIndex + 1
                do run = 1, inum_runs
                    TotParts(run) = TotParts(run) + abs(NoWalkers)
                    NoatHF(run) = NoWalkers
                end do
            else
                call stop_all(this_routine, "No walkers initialised on the HF det with StartMP1")
            end if
        else
            NoatHF(:) = 0.0_dp
        end if

        TotWalkers = DetIndex - 1   !This is the number of occupied determinants on each node
        TotWalkersOld = TotWalkers

        !Set local&global variables
        TotPartsOld = TotParts
        call mpisumall(TotParts, AllTotParts)
        call mpisumall(NoatHF, AllNoatHF)
        call mpisumall(TotWalkers, AllTotWalkers)
        OldAllNoatHF = AllNoatHF
        do run = 1, inum_runs
            OldAllHFCyc(run) = AllNoatHF(run)
            OldAllAvWalkersCyc(run) = AllTotParts(run)
        end do
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld = 0.0_dp

    end subroutine InitFCIMC_MP1

    SUBROUTINE CheckforBrillouins()
        INTEGER :: i, j
        LOGICAL :: tSpinPair

!Standard cases.
        IF ((tHub .and. tReal) .or. (tRotatedOrbs) .or. ((LMS /= 0) .and. (.not. tUHF)) .or. tReltvy) THEN
!Open shell, restricted.
            tNoBrillouin = .true.
        ELSE
!Closed shell restricted, or open shell unrestricted are o.k.
            tNoBrillouin = .false.
            tUseBrillouin = .true.
        end if

!Special case of complex orbitals.
        IF (tFixLz .and. (.not. tNoBrillouin)) THEN
            write(iout, *) "Turning Brillouins theorem off since we are using non-canonical complex orbitals"
            tNoBrillouin = .true.
        end if

        ! Special case of defining a det with LMS=0, but which is open shell.
        ! No Brillouins if it's a restricted HF calc.
        tSpinPair = .false.
        IF (tDefineDet .and. (LMS == 0) .and. (.not. tUHF)) THEN
            ! If we are defining our own reference determinant, we want to
            ! find out if it is open shell or closed to know whether or not
            ! brillouins theorem holds.
            !
            ! If LMS/=0, then it is easy and must be open shell, otherwise
            ! we need to consider the occupied orbitals.
            do i = 1, (NEl - 1), 2
                ! Assuming things will probably go alpha beta alpha beta,
                ! run through each alpha and see if there's a corresponding
                ! beta.
                tSpinPair = .false.
                IF (MOD(BRR(FDet(i)), 2) /= 0) THEN
!Odd energy, alpha orbital.
                    IF (BRR(FDet(i + 1)) /= (BRR(FDet(i)) + 1)) THEN
                        ! Check the next orbital to see if it's the beta (will
                        ! be alpha+1 when ordered by energy). If not, check
                        ! the other orbitals for the beta, as it's possible
                        ! the orbitals are ordered weird (?).
                        do j = 1, NEl
                            IF (BRR(FDet(j)) == (BRR(FDet(i)) + 1)) tSpinPair = .true.
                        end do
                    ELSE
                        tSpinPair = .true.
                    end if
                ELSE
!Even energy, beta orbital. The corresponding alpha will be beta-1.
                    IF (BRR(FDet(i + 1)) /= (BRR(FDet(i)) - 1)) THEN
                        do j = 1, NEl
                            IF (BRR(FDet(j)) == (BRR(FDet(i)) - 1)) tSpinPair = .true.
                        end do
                    ELSE
                        tSpinPair = .true.
                    end if
                end if
                IF (.not. tSpinPair) EXIT
            end do
            IF (.not. tSpinPair) THEN
!Open shell LMS=0 determinant.
!If restricted HF orbitals are being used, brillouins theorem does not hold.
                tNoBrillouin = .true.
                tUseBrillouin = .false.
                write(iout, '(A)') " Using an open shell reference determinant in a basis of restricted HF orbitals; " &
                & //"Brillouins theorem is being turned off. "
            end if
        end if

    ENDSUBROUTINE CheckforBrillouins

    SUBROUTINE CalcApproxpDoubles()
        implicit none
        real(dp) :: denom
        INTEGER :: iTotal
        integer :: nSingles, nDoubles, nSing_spindiff1, nDoub_spindiff1, nDoub_spindiff2
        integer :: nTot
        character(*), parameter :: this_routine = "CalcApproxpDoubles"
        integer(n_int), allocatable :: dummy_list(:, :)

        nSingles = 0
        nDoubles = 0
        if (tReltvy) then
            nSing_spindiff1 = 0
            nDoub_spindiff1 = 0
            nDoub_spindiff2 = 0
            pDoub_spindiff1 = 0.0_dp
            pDoub_spindiff2 = 0.0_dp
        end if

!NSing=Number singles from HF, nDoub=No Doubles from HF

        write(iout, "(A)") " Calculating approximate pDoubles for use with &
                       &excitation generator by looking a excitations from &
                       &reference."
        exflag = 3
        if (tReltvy) then
            write(iout, *) "Counting magnetic excitations"
            ! subroutine CountExcitations4(nI, minRank, maxRank, minSpinDiff, maxSpinDiff, tot)
            call CountExcitations4(hfdet, 1, 1, 0, 0, nSingles)
            call CountExcitations4(hfdet, 1, 1, 1, 1, nSing_spindiff1)
            call CountExcitations4(hfdet, 2, 2, 0, 0, nDoubles)
            call CountExcitations4(hfdet, 2, 2, 1, 1, nDoub_spindiff1)
            call CountExcitations4(hfdet, 2, 2, 2, 2, nDoub_spindiff2)
            call CountExcitations4(hfdet, 1, 2, 0, 2, nTot)
            ASSERT(nTot == (nSingles + nSing_spindiff1 + nDoubles + nDoub_spindiff1 + nDoub_spindiff2))

            iTotal = nSingles + nDoubles + nSing_spindiff1 + nDoub_spindiff1 + nDoub_spindiff2

        else
            if (tKPntSym) THEN
                if (t_k_space_hubbard) then
                    ! change this to the new implementation
                    call gen_all_excits_k_space_hubbard(HFDet, nDoubles, dummy_list)
                else
                    call enumerate_sing_doub_kpnt(exFlag, .false., nSingles, nDoubles, .false.)
                end if
            else
                call CountExcitations3(hfdet, exflag, nSingles, nDoubles)
            end if
            iTotal = nSingles + nDoubles
        end if

        IF (tHub .or. tUEG) THEN
            IF (tReal) THEN
                write(iout, *) "Since we are using a real-space hubbard model, only single excitations are connected &
                &   and will be generated."
                pDoubles = 0.0_dp
                if (tReltvy) then
                    pDoub_spindiff1 = 0.0_dp
                    pDoub_spindiff2 = 0.0_dp
                    pSingles = real(nSingles, dp) / real(nSingles + nSing_spindiff1, dp)
                    pSing_spindiff1 = 1.0_dp - pSingles
                else
                    pSingles = 1.0_dp
                end if
                return
            ELSE
                write(iout, *) "Since we are using a momentum-space hubbard model/UEG, only double excitaitons &
     &                          are connected and will be generated."
                pSingles = 0.0_dp
                if (tReltvy) then
                    pSing_spindiff1 = 0.0_dp
                    pDoubles = real(nDoubles, dp) / real(nDoubles + nDoub_spindiff1 + nDoub_spindiff2, dp)
                    pDoub_spindiff1 = real(nDoub_spindiff1, dp) / real(nDoubles + nDoub_spindiff1 + nDoub_spindiff2, dp)
                    pDoub_spindiff2 = 1.0_dp - pDoubles - pDoub_spindiff1
                else
                    pDoubles = 1.0_dp
                end if
                return
            end if

        else if (tNoSingExcits) then
            pSingles = 0.0_dp
            if (tReltvy) then
                pSing_spindiff1 = 0.0_dp
                pDoubles = real(nDoubles, dp) / real(nDoubles + nDoub_spindiff1 + nDoub_spindiff2, dp)
                pDoub_spindiff1 = real(nDoub_spindiff1, dp) / real(nDoubles + nDoub_spindiff1 + nDoub_spindiff2, dp)
                pDoub_spindiff2 = 1.0_dp - pDoubles - pDoub_spindiff1
            else
                pDoubles = 1.0_dp
            end if
            write(iout, *) "Only double excitations will be generated"
            return
        end if

        write(iout, "(I7,A,I7,A)") nDoubles, " double excitations, and ", nSingles, &
            " single excitations found from reference. This will be used to calculate pDoubles."

        IF (abs(SinglesBias - 1.0_dp) > 1.0e-12_dp) THEN
            write(iout, *) "Singles Bias detected. Multiplying single excitation connectivity of HF determinant by ", &
                SinglesBias, " to determine pDoubles."
        end if

        IF ((nSingles == 0) .or. (nDoubles == 0)) THEN
            write(iout, *) "Number of singles or doubles found equals zero. pDoubles will be set to 0.95. Is this correct?"
            pDoubles = 0.95_dp
            pSingles = 0.05_dp
            return
        else if ((nSingles < 0) .or. (nDoubles < 0)) then
            call stop_all("CalcApproxpDoubles", &
                          "Number of singles, doubles or Yamanouchi symbols &
                          &found to be a negative number. Error here.")
        end if
        if (tReltvy) then
            denom = real(nSingles + nSing_spindiff1, dp) * SinglesBias &
                    + real(nDoubles + nDoub_spindiff1 + nDoub_spindiff2, dp)
            pSingles = real(nSingles, dp) * SinglesBias / denom
            pSing_spindiff1 = real(nSing_spindiff1, dp) * SinglesBias / denom
            pDoubles = real(nDoubles, dp) / denom
            pDoub_spindiff1 = real(nDoub_spindiff1, dp) / denom
            pDoub_spindiff2 = 1.0_dp - pSingles - pSing_spindiff1 &
                              - pDoubles - pDoub_spindiff1 - pDoub_spindiff2
        else
            denom = real(nSingles, dp) * SinglesBias + real(nDoubles, dp)
            pSingles = real(nSingles, dp) * SinglesBias / denom
            pDoubles = 1.0_dp - pSingles
        end if

        IF (abs(SinglesBias - 1.0_dp) > 1.0e-12_dp) THEN
            write (iout, '("pDoubles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pDoubles, real(nDoubles, dp) / real(iTotal, dp)
            write (iout, '("pSingles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pSingles, real(nSingles, dp) / real(iTotal, dp)

!            write(iout,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ", &
!                & real(nDoub,dp)/real(iTotal,dp)
        ELSE
            if (tReltvy) then
                write(iout, '(A)') " Where s and t are alpha or beta spin function labels: "
                write(iout, '(A30,F14.6)') " pSingles(s->s) set to: ", pSingles
                write(iout, '(A30,F14.6)') " pSingles(s->s') set to: ", pSing_spindiff1
                write(iout, '(A30,F14.6)') " pDoubles(st->st) set to: ", pDoubles
                write(iout, '(A30,F14.6)') " pDoubles(st->s't) set to: ", pDoub_spindiff1
                write(iout, '(A30,F14.6)') " pDoubles(st->s't') set to: ", pDoub_spindiff2
            else
                write(iout, '(A,F14.6)') " pDoubles set to: ", pDoubles
                write(iout, '(A,F14.6)') " pSingles set to: ", pSingles
            end if
        end if

        if (pSinglesIn > 1.e-12_dp) then
            pSingles = pSinglesIn
            pDoubles = 1.0_dp - pSinglesIn
            write(iout, '(" Using the input value of pSingles:",1x, f14.6)') pSingles
            write(iout, '(" Using the input value of pSingles:",1x, f14.6)') pDoubles
        end if
        if (pParallelIn > 1.e-12_dp) then
            write(iout, '(" Using the input value of pSingles:",1x, f14.6)') pParallelIn
            pParallel = pParallelIn
        end if

    END SUBROUTINE CalcApproxpDoubles

    SUBROUTINE CreateSpinInvBRR()

        ! Create an SpinInvBRR containing spin orbitals,
        ! unlike 'createInvBRR' which only has spatial orbitals.
        ! This is used for the FixCASshift option in establishing whether or not
        ! a determinant is in the complete active space.
        ! In:
        !    BRR(i)=j: orbital i is the j-th lowest in energy.
        !    nBasis: size of basis
        ! SpinInvBRR is the inverse of BRR.  SpinInvBRR(j)=i: the j-th lowest energy
        ! orbital corresponds to the i-th orbital in the original basis.
        ! i.e the position in SpinInvBRR now corresponds to the orbital number and
        ! the value to the relative energy of this orbital.

        IMPLICIT NONE
        INTEGER :: I, t, ierr
        CHARACTER(len=*), PARAMETER :: this_routine = 'CreateSpinInvBrr'

        IF (ALLOCATED(SpinInvBRR)) return

        allocate(SpinInvBRR(NBASIS), STAT=ierr)
        CALL LogMemAlloc('SpinInvBRR', NBASIS, 4, this_routine, SpinInvBRRTag, ierr)

        SpinInvBRR(:) = 0

        t = 0
        do I = 1, NBASIS
            t = t + 1
            SpinInvBRR(BRR(I)) = t
        end do

        return

    END SUBROUTINE CreateSpinInvBRR

    subroutine SetupValidSpawned(WalkerListSize)

        use CalcData, only: MemoryFacSpawn

        implicit none
        integer(int64), intent(in) :: WalkerListSize
        integer ierr, j
        real(dp) Gap

        !When running normally, WalkerListSize will be equal to totalwalkers
        !However, when reading in (and not continuing to grow) it should be equal to the number of dets in the popsfile
        MaxSpawned = NINT(MemoryFacSpawn * WalkerListSize * inum_runs)
!            write(iout,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned

!      write(iout,"(A)") "*Direct Annihilation* in use...Explicit load-balancing disabled."
        allocate(ValidSpawnedList(0:nNodes - 1), stat=ierr)
        ! InitialSpawnedSlots is now filled later, once the number of particles
        ! wanted is known
        !(it can change according to the POPSFILE).
        allocate(InitialSpawnedSlots(0:nNodes - 1), stat=ierr)
        ! InitialSpawnedSlots now holds the first free position in the
        ! newly-spawned list for each processor, so it does not need to be
        ! reevaluated each iteration.
!      MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
        Gap = REAL(MaxSpawned, dp) / REAL(nNodes, dp)
        do j = 0, nNodes - 1
            InitialSpawnedSlots(j) = NINT(Gap * j) + 1
        end do
        ! ValidSpawndList now holds the next free position in the newly-spawned
        ! list, but for each processor.
        ValidSpawnedList(:) = InitialSpawnedSlots(:)

    end subroutine SetupValidSpawned

    subroutine sync_rdm_sampling_iter()
        use LoggingData, only: RDMEnergyIter, IterRDMOnfly
        use CalcData, only: coreSpaceUpdateCycle, semistoch_shift_iter
        implicit none
        integer :: frac
        ! first, adjust the offset to make the rdm sampling start right when a semistochastic
        ! update cycle ends
        IterRDMOnFly = IterRDMOnFly - &
                       mod(IterRDMOnFly - semistoch_shift_iter, coreSpaceUpdateCycle) - 1
        ! The -1 is just because the sampling starts one iteration after IterRDMOnFly
        ! If we subtracted too much, jump one cycle backwards
        if (IterRDMOnFly < semistoch_shift_iter) IterRDMOnFly = IterRDMOnFly + coreSpaceUpdateCycle
        write(6, *) "Adjusted starting iteration of RDM sampling to ", IterRDMOnFly

        ! Now sync the update cycles
        if (RDMEnergyIter > coreSpaceUpdateCycle) then
            RDMEnergyIter = coreSpaceUpdateCycle
            write(6, *) "The RDM sampling interval cannot be larger than the update " &
                //"interval of the semi-stochastic space. Reducing it to ", RDMEnergyIter
        end if
        if (mod(coreSpaceUpdateCycle, RDMEnergyIter) /= 0) then
            ! first, try to ramp up the RDMEnergyIter to meet the coreSpaceUpdateCycle
            frac = coreSpaceUpdateCycle / RDMEnergyIter
            RDMEnergyIter = coreSpaceUpdateCycle / frac
            write(6, *) "Update cycle of semi-stochastic space and RDM sampling interval" &
                //" out of sync. "
            write(6, *) "Readjusting RDM sampling interval to ", RDMEnergyIter

            ! now, if this did not succeed, adjust the coreSpaceUpdateCycle
            if (mod(coreSpaceUpdateCycle, RDMEnergyIter) /= 0) then
                coreSpaceUpdateCycle = coreSpaceUpdateCycle - &
                                       mod(coreSpaceUpdateCycle, RDMEnergyIter)
                write(6, *) "Adjusted update cycle of semi-stochastic space to ", &
                    coreSpaceUpdateCycle
            end if
        end if
    end subroutine sync_rdm_sampling_iter

    subroutine CalcUEGMP2()
        use SymExcitDataMod, only: kPointToBasisFn
        use SystemData, only: ElecPairs, NMAXX, NMAXY, NMAXZ, OrbECutOff, &
                              tMP2UEGRestrict, kiRestrict, kiMsRestrict, kjRestrict, kjMsRestrict, &
                              Madelung, tMadelung, tUEGFreeze, FreezeCutoff, kvec, tUEG2
        use Determinants, only: GetH0Element4, get_helement_excit
        integer :: Ki(3), Kj(3), Ka(3), LowLoop, HighLoop, X, i, Elec1Ind, Elec2Ind, K, Orbi, Orbj
        integer :: iSpn, FirstA, nJ(NEl), a_loc, Ex(2, maxExcit), kx, ky, kz, OrbB
        integer :: ki2, kj2
        logical :: tParity
        real(dp) :: Ranger, mp2, mp2all, length
        HElement_t(dp) :: hel, H0tmp

        !Divvy up the ij pairs
        Ranger = real(ElecPairs, dp) / real(nProcessors, dp)
        LowLoop = int(iProcIndex * Ranger) + 1
        Highloop = int((iProcIndex + 1) * Ranger)

        if ((iProcIndex + 1) == nProcessors) Highloop = ElecPairs
        if (iProcIndex == 0) then
            if (lowLoop /= 1) write(iout, *) "Error here!"
        end if
        write(iout, *) "Total ij pairs: ", ElecPairs
        write(iout, *) "Considering ij pairs from: ", LowLoop, " to ", HighLoop
!        write(iout,*) "HFDet: ",HFDet(:)

        do i = LowLoop, HighLoop   !Looping over electron pairs on this processor

            X = ElecPairs - i
            K = INT((SQRT(8.0_dp * REAL(X, dp) + 1.0_dp) - 1.0_dp) / 2.0_dp)
            Elec1Ind = NEl - 1 - K
            Elec2Ind = NEl - X + ((K * (K + 1)) / 2)
            Orbi = HFDet(Elec1Ind)
            Orbj = HFDet(Elec2Ind)
            Ki = G1(Orbi)%k
            Kj = G1(Orbj)%k
            !=======================================
            if (tUEG2) then
                Ki = kvec(Orbi, 1:3)
                Kj = kvec(Orbj, 1:3)
            end if
            !=======================================
            if (tUEGFreeze) then
                ki2 = ki(1)**2 + ki(2)**2 + ki(3)**2
                kj2 = kj(1)**2 + kj(2)**2 + kj(3)**2
                if (.not. (ki2 > FreezeCutoff .and. kj2 > FreezeCutoff)) cycle
            end if
            if (tMP2UEGRestrict) then
                if (.not. ( &
                    (kiRestrict(1) == ki(1) .and. kiRestrict(2) == ki(2) .and. kiRestrict(3) == ki(3) .and. &
                     kjRestrict(1) == kj(1) .and. kjRestrict(2) == kj(2) .and. kjRestrict(3) == kj(3) .and. &
                     kjMsRestrict == G1(Orbi)%Ms .and. kiMsRestrict == G1(Orbj)%Ms) .or. &
                    ! the other way round
                    (kiRestrict(1) == kj(1) .and. kiRestrict(2) == kj(2) .and. kiRestrict(3) == kj(3) .and. &
                     kjRestrict(1) == ki(1) .and. kjRestrict(2) == ki(2) .and. kjRestrict(3) == ki(3) .and. &
                     kiMsRestrict == G1(Orbi)%Ms .and. kjMsRestrict == G1(Orbj)%Ms)) &
                    ) cycle
                write(iout, *) "Restricting calculation to i,j pair: ", Orbi, Orbj
            end if

            IF ((G1(Orbi)%Ms) * (G1(Orbj)%Ms) == -1) THEN
!We have an alpha beta pair of electrons.
                iSpn = 2
            ELSE
                IF (G1(Orbi)%Ms == 1) THEN
!We have an alpha alpha pair of electrons.
                    iSpn = 3
                ELSE
!We have a beta beta pair of electrons.
                    iSpn = 1
                end if
            end if

!            write(iout,*) "ijpair: ",Orbi,Orbj

            if ((iSpn == 3) .or. (iSpn == 1)) then
                if (iSpn == 3) then
                    FirstA = 2    !Loop over alpha
                else
                    FirstA = 1    !Loop over beta
                end if

                do a_loc = FirstA, nBasis, 2
                    !Loop over all a

                    !Reject if a is occupied
                    if (IsOcc(iLutHF, a_loc)) cycle

                    Ka = G1(a_loc)%k
                    !=======================================
                    if (tUEG2) then
                        Ka = kvec(a_loc, 1:3)
                    end if
                    !=======================================

                    !Find k labels of b
                    kx = Ki(1) + Kj(1) - Ka(1)
                    if (abs(kx) > NMAXX) cycle
                    ky = Ki(2) + Kj(2) - Ka(2)
                    if (abs(ky) > NMAXY) cycle
                    kz = Ki(3) + Kj(3) - Ka(3)
                    if (abs(kz) > NMAXZ) cycle
                    length = real((kx**2) + (ky**2) + (kz**2), dp)
                    if (length > OrbECutoff) cycle

                    !Find the actual k orbital
                    if (iSpn == 3) then
                        !want alpha
                        OrbB = kPointToBasisFn(kx, ky, kz, 2)
                    else
                        !want beta
                        OrbB = kPointToBasisFn(kx, ky, kz, 1)
                    end if

                    !Reject k orbital if it is occupied or gt a
                    if (IsOcc(iLutHF, OrbB)) cycle
                    if (OrbB >= a_loc) cycle

                    !Find det
!                    write(iout,*) "OrbB: ",OrbB
                    call make_double(HFDet, nJ, elec1ind, elec2ind, a_loc, &
                                     orbB, ex, tParity)
                    !Sum in mp2 contrib
                    hel = get_helement_excit(HFDet, nJ, 2, Ex, tParity)

                    H0tmp = getH0Element4(nJ, HFDet)
                    H0tmp = Fii - H0tmp
                    if (tMadelung) then
                        H0tmp = H0tmp + 2.0_dp * Madelung
                    end if
                    mp2 = mp2 + (hel**2) / H0tmp
!                    write(iout,*) (hel**2),H0tmp
                end do

            else if (iSpn == 2) then
                do a_loc = 1, nBasis
                    !Loop over all a_loc
!                    write(iout,*) "a_loc: ",a_loc

                    !Reject if a is occupied
                    if (IsOcc(iLutHF, a_loc)) cycle

                    Ka = G1(a_loc)%k
                    !=======================================
                    if (tUEG2) then
                        Ka = kvec(a_loc, 1:3)
                    end if
                    !=======================================

                    !Find k labels of b
                    kx = Ki(1) + Kj(1) - Ka(1)
                    if (abs(kx) > NMAXX) cycle
                    ky = Ki(2) + Kj(2) - Ka(2)
                    if (abs(ky) > NMAXY) cycle
                    kz = Ki(3) + Kj(3) - Ka(3)
                    if (abs(kz) > NMAXZ) cycle
                    length = real((kx**2) + (ky**2) + (kz**2), dp)
                    if (length > OrbECutoff) cycle

                    !Find the actual k orbital
                    if (is_beta(a_loc)) then
                        !want alpha b orbital
                        OrbB = kPointToBasisFn(kx, ky, kz, 2)
                    else
                        !want beta
                        OrbB = kPointToBasisFn(kx, ky, kz, 1)
                    end if

                    !Reject k orbital if it is occupied or gt a
                    if (IsOcc(iLutHF, OrbB)) cycle
                    if (OrbB >= a_loc) cycle

!                    write(iout,*) "OrbB: ",OrbB
                    !Find det
                    call make_double(HFDet, nJ, elec1ind, elec2ind, a_loc, &
                                     orbB, ex, tParity)
                    !Sum in mp2 contrib
                    hel = get_helement_excit(HFDet, nJ, 2, Ex, tParity)
                    H0tmp = getH0Element4(nJ, HFDet)
                    H0tmp = Fii - H0tmp
                    if (tMadelung) then
                        H0tmp = H0tmp + 2.0_dp * Madelung
                    end if
                    mp2 = mp2 + (hel**2) / H0tmp
!                    write(iout,*) (hel**2),H0tmp
                end do
            end if

        end do

!        write(iout,*) "mp2: ",mp2
        mp2all = 0.0_dp

        !Sum contributions across nodes.
        call MPISumAll(mp2, mp2all)
        write(iout, "(A,2G25.15)") "MP2 energy calculated: ", MP2All, MP2All + Hii
        call neci_flush(iout)

    end subroutine CalcUEGMP2

    !Ensure that the new FCIMCStats file which is about to be opened does not overwrite any other FCIMCStats
    !files. If there is already an FCIMCStats file present, then move it to FCIMCStats.x, where x is a largest
    !free filename.
    subroutine MoveFCIMCStatsFiles()
#ifdef NAGF95
        USe f90_unix_dir, only: rename
#endif
        integer :: extension
        logical :: exists
        character(len=22) :: abstr
!        character(len=36) :: command
        character(len=*), parameter :: t_r = 'MoveFCIMCStatsFiles'

        if (tMolpro) then
            inquire (file='FCIQMCStats', exist=exists)
        else
            inquire (file='FCIMCStats', exist=exists)
        end if
        if (exists) then
            !We already have an FCIMCStats file - move it to the end of the list of FCIMCStats files.
            extension = 1
            do while (.true.)
                if (tMolpro) then
                    abstr = 'FCIQMCStats.'//str(extension)
                else
                    abstr = 'FCIMCStats.'//str(extension)
                end if
                inquire (file=trim(adjustl(abstr)), exist=exists)
                if (.not. exists) exit
                extension = extension + 1
                if (extension > 10000) then
                    call stop_all(t_r, "Error finding free FCIMCStats name")
                end if
            end do

            !We have got a unique filename
            !Do not use system call
!            command = 'mv' // ' FCIMCStats ' // abstr
!            stat = neci_system(trim(command))

            if (tMolpro) then
                call rename('FCIQMCStats', trim(adjustl(abstr)))
            else
                call rename('FCIMCStats', trim(adjustl(abstr)))
            end if
            !Doesn't like the stat argument
!            if(stat.ne.0) then
!                call stop_all(t_r,"Error with renaming FCIMCStats file")
!            end if
        end if

    end subroutine MoveFCIMCStatsFiles

    subroutine assign_reference_dets()

        ! Depending on the configuration we may have one, or multiple,
        ! reference determinants.

        integer :: det(nel), orbs(nel), orb, orb2, norb
        integer :: run, cc_idx, label_idx, i, j, found_orbs(inum_runs)
        real(dp) :: energies(nel), hdiag

        ! for now guga only works with non-complex code
        integer(n_int), pointer :: excitations(:, :)
        integer :: n_excits, ierr
        real(dp), allocatable :: diag_energies(:)
        logical, allocatable :: found_mask(:)
        character(*), parameter :: this_routine = "assign_reference_dets"

        ! If the user has specified all of the (multiple) reference states,
        ! then just copy these across to the ilutRef array:
        if (tMultipleInitialStates) then

            tReplicaReferencesDiffer = .true.

            do run = 1, inum_runs
                ProjEDet(:, run) = initial_states(:, run)
                call EncodeBitDet(ProjEDet(:, run), ilutRef(:, run))
            end do

        else if (tOrthogonaliseReplicas) then

            tReplicaReferencesDiffer = .true.

            ! The first replica is just a normal FCIQMC simulation.
            ilutRef(:, 1) = ilutHF
            ProjEDet(:, 1) = HFDet

            found_orbs = 0

            ! in the GUGA approach we have to change that simple
            ! single excitation of course or otherwise we get non-
            ! allowed CSF or the wrong STOT symmetry sector.
            if (tGUGA) then
                ! run the exact single excitations on the HF det
                ! and find inum_runs lowest energetically ones..

                ! create all excitations from the HF
                ! do i remain in the same symmetry sector with this
                ! info??
                ! i think my exact hamil application routine does not
                ! deal with symmetry at all...
                ! i guess i have to change that for the real-space
                ! hubbard model implementation!
                if (tHUB .or. tUEG .or. .not. (tNoBrillouin)) then
                    call actHamiltonian(ilutHF, excitations, n_excits)
                else
                    call actHamiltonian(ilutHF, excitations, n_excits, .true.)
                end if

                ! if no excitations possible... there is something wrong
                if (.not. (n_excits > 0) .or. n_excits < inum_runs - 1) then
                    print *, "n_excits: ", n_excits
                    print *, "requested excited states:", inum_runs - 1
                    call write_guga_list(6, excitations(:, 1:n_excits))
                    call stop_all(this_routine, "not enough excitations from HF!")
                end if

                ! and the choose the inum_runs - 1 lowest energetically
                ! excitations

                ! also have to calculate the correct diagonal element to
                ! compare energies
                allocate(diag_energies(n_excits), stat=ierr)
                allocate(found_mask(n_excits), stat=ierr)
                found_mask = .true.

                do i = 1, n_excits
                    diag_energies(i) = &
                        calcDiagMatEleGUGA_ilut(excitations(:, i))
                end do

                ! can i sort the excitation list, according to energies?
                do run = 2, inum_runs
                    ! find the minimum energy
                    i = minloc(diag_energies, 1, found_mask)
                    ! dont find that specific one again
                    found_mask(i) = .false.
                    ! assign the reference determinant(transform between
                    ! guga and actual iluts!)
                    call convert_ilut_toNECI(excitations(:, i), ilutRef(:, run))
                    ! and calc. the nI representation
                    call decode_bit_det(ProjEDet(:, run), ilutRef(:, run))
                    ! that should be it.. should the diagonal element also
                    ! be already stored within the ilut? is that done here?
                end do

            else
                do run = 2, inum_runs

                    ! Now we want to find the lowest energy single excitation with
                    ! the same symmetry as the reference site.

                    do i = 1, nel
                        ! Find the excitations, and their energy
                        orb = HFDet(i)
                        cc_idx = ClassCountInd(orb)
                        label_idx = SymLabelCounts2(1, cc_idx)
                        norb = OrbClassCount(cc_idx)

                        ! nb. sltcnd_0 does not depend on the ordering of the det,
                        !     so we don't need to do any sorting here.
                        energies(i) = 9999999.9_dp
                        do j = 1, norb
                            orb2 = SymLabelList2(label_idx + j - 1)
                            if ((.not. any(orb2 == HFDet)) .and. &
                                (.not. any(orb2 == found_orbs))) then
                                det = HFDet
                                det(i) = orb2
                                hdiag = real(sltcnd_excit(det, NoExc_t()), dp)
                                if (hdiag < energies(i)) then
                                    energies(i) = hdiag
                                    orbs(i) = orb2
                                end if
                            end if
                        end do
                    end do

                    ! Which of the electrons that is excited gives the lowest energy?
                    i = minloc(energies, 1)
                    found_orbs(run) = orbs(i)

                    ! Construct that determinant, and set it as the reference.
                    ProjEDet(:, run) = HFDet
                    ProjEDet(i, run) = orbs(i)
                    call sort(ProjEDet(:, run))
                    call EncodeBitDet(ProjEDet(:, run), ilutRef(:, run))

                end do
            end if

        else if (tPreCond) then

            do run = 1, inum_runs
                ilutRef(:, run) = ilutHF
                ProjEDet(:, run) = HFDet
            end do

            ! And make sure that the rest of the code knows this
            tReplicaReferencesDiffer = .true.

        else

            ! This is the normal case. All simultions are essentially doing
            ! the same thing...

            do run = 1, inum_runs
                ilutRef(:, run) = ilutHF
                ProjEDet(:, run) = HFDet
            end do

            ! And make sure that the rest of the code knows this
            tReplicaReferencesDiffer = .false.

        end if

        write(6, *) 'Generated reference determinants:'
        do run = 1, inum_runs
            call write_det(6, ProjEDet(:, run), .false.)
            hdiag = real(get_helement(ProjEDet(:, run), ProjEDet(:, run), 0), dp)
            write(6, '(" E = ", f16.9)') hdiag
        end do

    end subroutine assign_reference_dets

    subroutine init_cont_time()

        integer :: ierr
        character(*), parameter :: this_routine = 'init_cont_time'
        character(*), parameter :: t_r = this_routine

        call clean_cont_time()

        allocate(oversample_factors(1:2, LMS:nel), stat=ierr)
        @:log_alloc(oversample_factors, ostag, ierr)
        oversample_factors = 1.0_dp

        ! We need somewhere for our nested excitation generators to call home
        call init_excit_gen_store(secondary_gen_store)

    end subroutine

    subroutine clean_cont_time()

        character(*), parameter :: this_routine = 'clean_cont_time'

        if (allocated(oversample_factors)) then
            deallocate(oversample_factors)
            log_dealloc(ostag)
        end if

    end subroutine

!------------------------------------------------------------------------------------------!

    subroutine setup_adi()
        ! We initialize the flags for the adi feature
        use adi_data, only: tSetDelayAllDoubsInits, tDelayAllDoubsInits, &
                            tAllDoubsInitiators, tDelayGetRefs, &
                            NoTypeN, tReadRefs, maxNRefs, SIUpdateOffset
        use CalcData, only: InitiatorWalkNo
        use adi_references, only: enable_adi, reallocate_ilutRefAdi, &
                                  reset_coherence_counter
        implicit none
        call reallocate_ilutRefAdi(maxNRefs)

        ! If using adi with dynamic SIs, also use a dynamic corespace by default
        call setup_dynamic_core()

        ! Check if one of the keywords is specified as delayed
        if (tSetDelayAllDoubsInits .and. tAllDoubsInitiators) then
            tAllDoubsInitiators = .false.
            tDelayAllDoubsInits = .true.
        end if

        ! Check if we want to get the references right away
        if (.not. (tReadRefs .or. tReadPops)) tDelayGetRefs = .true.
        if (tDelayAllDoubsInits) tDelayGetRefs = .true.
        ! Give a status message
        if (tAllDoubsInitiators) then
            call enable_adi()
            tAdiActive = .true.
        end if

        ! there is a minimum cycle lenght for updating the number of SIs, as the reference population
        ! needs some time to equilibrate
        nRefUpdateInterval = max(SIUpdateInterval, 500)
        SIUpdateOffset = 0

        ! Initialize the logging variables
        call reset_coherence_counter()
    end subroutine setup_adi

!------------------------------------------------------------------------------------------!

    subroutine setup_dynamic_core()
        use CalcData, only: tDynamicCoreSpace, coreSpaceUpdateCycle, tIntervalSet
        use adi_data, only: tAllDoubsInitiators
        implicit none

        ! Enable dynamic corespace if both
        ! a) using adi with dynamic SIs (default)
        ! b) no other keywords regarding the dynamic corespace are given
        if (SIUpdateInterval > 0 .and. .not. tIntervalSet .and. tAllDoubsInitiators) then
            tDynamicCoreSpace = .true.
            coreSpaceUpdateCycle = SIUpdateInterval
        end if
    end subroutine setup_dynamic_core

!------------------------------------------------------------------------------------------!

    subroutine init_norm()
        use bit_rep_data, only: test_flag_multi
        ! initialize the norm_psi, norm_psi_squared
        implicit none
        integer(int64) :: j
        real(dp) :: sgn(lenof_sign)
        logical :: tIsStateDeterm

        norm_psi_squared = 0.0_dp
        norm_semistoch_squared = 0.0_dp
        ! has to be set only once, if it changes in one iteration, it is reset in every iteration
        tIsStateDeterm = .false.
        do j = 1, TotWalkers
            ! get the sign
            call extract_sign(CurrentDets(:, j), sgn)
            if (tSemiStochastic) tIsStateDeterm = test_flag_multi(CurrentDets(:, j), flag_deterministic)
            call addNormContribution(sgn, tIsStateDeterm)
        end do

        ! sum up the norm over the procs
        call MPISumAll(norm_psi_squared, all_norm_psi_squared)

        ! assign the sqrt norm
#ifdef CMPLX_
        norm_psi = sqrt(sum(all_norm_psi_squared))
        norm_semistoch = sqrt(sum(norm_semistoch_squared))
#else
        norm_psi = sqrt(all_norm_psi_squared)
        norm_semistoch = sqrt(norm_semistoch_squared)
#endif

        old_norm_psi = norm_psi

        all_norms = all_norm_psi_squared
    end subroutine init_norm

!------------------------------------------------------------------------------------------!

end module fcimc_initialisation

! This routine will change the reference determinant to DetCurr. It will
! also re-zero all the energy estimators, since they now correspond to
! projection onto a different determinant.
!
! n.b. NOT MODULARISED. This is a little evil, but there is an unbreakable
!      circular dependency otherwise.
!
! **** See interface in Popsfile.F90 ****
subroutine ChangeRefDet(DetCurr)
    use fcimc_initialisation
    INTEGER :: DetCurr(NEl), i

    do i = 1, NEl
        FDet(i) = DetCurr(i)
    end do

    write(iout, "(A)") "*** Changing the reference determinant ***"
    write(iout, "(A)") "Switching reference and zeroing energy counters - restarting simulation"
!
!Initialise variables for calculation on each node
    Iter = 1
    CALL DeallocFCIMCMemPar()
    IF (iProcIndex == Root) THEN
        close(fcimcstats_unit)
        if (inum_runs == 2) close(fcimcstats_unit2)
        IF (tTruncInitiator) close(initiatorstats_unit)
        IF (tLogComplexPops) close(complexstats_unit)
        if (tLogEXLEVELStats) close(EXLEVELStats_unit)
    end if
    IF (TDebug) close(11)
    CALL SetupParameters()
    CALL InitFCIMCCalcPar()

end subroutine ChangeRefDet
