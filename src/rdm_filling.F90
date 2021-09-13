#include "macros.h"

module rdm_filling

    ! This module contains routines used to perform filling of the RDM arrays,
    ! as done on-the-fly during an FCIQMC simulation.

    use bit_rep_data, only: NIfTot, nifd, test_flag, IlutBits, IlutBitsParent
    use bit_reps, only: get_initiator_flag_by_run
    use constants
    use SystemData, only: tGUGA, nbasis
    use guga_bitRepOps, only: extract_stochastic_rdm_info, extract_2_rdm_ind, &
                              init_csf_information, identify_excitation
    use rdm_data, only: rdm_spawn_t, rdmCorrectionFactor
    use CalcData, only: tAdaptiveShift, tNonInitsForRDMs, tInitsRDMRef, &
                        tNonVariationalRDMs
    use FciMCData, only: projEDet, ilutRef
    use DetBitOps, only: DetBitEq, GetBitExcitation
    use guga_rdm, only: Add_RDM_From_IJ_Pair_GUGA, fill_diag_1rdm_guga, &
                        fill_spawn_rdm_diag_guga, Add_RDM_HFConnections_GUGA, &
                        fill_sings_1rdm_guga, fill_sings_2rdm_guga, &
                        add_rdm_from_ij_pair_guga_exact
    use util_mod, only: near_zero
    use guga_data, only: excit_type, ExcitationInformation_t
    use LoggingData, only: t_full_core_rdms

    implicit none

contains

    subroutine fill_rdm_diag_wrapper(rdm_defs, spawn, one_rdms, ilut_list, ndets, tNonInit, &
                                     tLagrCorr)

        ! Loop over all states in ilut_list and see if any signs have just
        ! become unoccupied or become reoccupied. In which case, we have
        ! started a new averaging block, so we need to add in the
        ! contributions from the last block to the corresponding RDMs.

        use bit_rep_data, only: extract_sign
        use bit_reps, only: all_runs_are_initiator
        use CalcData, only: tPairedReplicas
        use global_det_data, only: get_iter_occ_tot, get_av_sgn_tot
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: one_rdm_t, rdm_definitions_t

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: ilut_list(:, :)
        integer, intent(in) :: ndets
        logical, intent(in), optional :: tNonInit, tLagrCorr

        integer :: idet, irdm, av_ind_1, av_ind_2
        real(dp) :: curr_sign(lenof_sign), adapted_sign(len_av_sgn_tot)
        real(dp) :: av_sign(len_av_sgn_tot), iter_occ(len_iter_occ_tot)
        logical :: tAllContribs, tLC

        def_default(tAllContribs, tNonInit, .true.)
        def_default(tLC, tLagrCorr, .true.)

        associate(ind => rdm_defs%sim_labels)

            do idet = 1, ndets

                call extract_sign(ilut_list(:, idet), curr_sign)
                ! All average sign from all RDMs.
                av_sign = get_av_sgn_tot(idet)
                ! The iteration on which each replica became occupied.
                iter_occ = get_iter_occ_tot(idet)

                adapted_sign = 0.0_dp

                do irdm = 1, rdm_defs%nrdms
                    ! The indicies of the first and second replicas in this
                    ! particular pair, in the *average* sign arrays (and
                    ! therefore also for the iter_occ array).
                    av_ind_1 = irdm * 2 - 1
                    av_ind_2 = irdm * 2

                    if ((abs(curr_sign(ind(1, irdm))) < 1.0e-10_dp .and. abs(iter_occ(av_ind_1)) > 1.0e-10_dp) .or. &
                        (abs(curr_sign(ind(2, irdm))) < 1.0e-10_dp .and. abs(iter_occ(av_ind_2)) > 1.0e-10_dp) .or. &
                        (abs(curr_sign(ind(1, irdm))) > 1.0e-10_dp .and. abs(iter_occ(av_ind_1)) < 1.0e-10_dp) .or. &
                        (abs(curr_sign(ind(2, irdm))) > 1.0e-10_dp .and. abs(iter_occ(av_ind_2)) < 1.0e-10_dp)) then

                        ! In this case we want to include this diagonal element,
                        ! so transfer the sign.
                        adapted_sign(av_ind_1:av_ind_2) = av_sign(av_ind_1:av_ind_2)
                    end if
                end do

                ! At least one of the signs has just gone to zero or just become
                ! reoccupied, so we need to add in diagonal elements and connections to HF
                if (any(abs(adapted_sign) > 1.e-12_dp)) then
                    if (tAllContribs .or. all_runs_are_initiator(ilut_list(:, idet))) &
                        call det_removed_fill_diag_rdm(spawn, one_rdms, ilut_list(:, idet), adapted_sign, iter_occ, tLC)
                end if

            end do

        end associate

    end subroutine fill_rdm_diag_wrapper

    subroutine fill_rdm_diag_currdet_norm(spawn, one_rdms, iLutnI, nI, ExcitLevelI, av_sign, iter_occ, tCoreSpaceDet, tLagrCorr)

        ! This routine calculates the diagonal RDM contribution, and explicit
        ! connections to the HF, from the current determinant.

        ! j --> Which element of the main list CurrentDets are we considering?
        ! IterLastRDMFill is the number of iterations since the last time the
        ! RDM contributions were added in (often the frequency of the RDM
        ! energy calculation).

        ! For the instantaneous RDMs we need to multiply the RDM contributions
        ! by either this, or the number of iterations the determinant has been
        ! occupied, which ever is fewer. For the full RDMs we need to multiply
        ! the RDM contributions by the number of iterations the determinant has
        ! been occupied.

        use bit_reps, only: decode_bit_det, all_runs_are_initiator
        use DetBitOps, only: TestClosedShellDet, FindBitExcitLevel
        use FciMCData, only: Iter, IterRDMStart, PreviousCycles, AvNoAtHF
        use global_det_data, only: len_iter_occ_tot
        use hphf_integrals, only: hphf_sign
        use HPHFRandExcitMod, only: FindExcitBitDetSym
        use LoggingData, only: RDMEnergyIter, RDMExcitLevel
        use rdm_data, only: one_rdm_t
        use SystemData, only: nel, tHPHF
        use CalcData, only: tNonInitsForRDMs

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: nI(nel), ExcitLevelI
        real(dp), intent(in) :: av_sign(:), iter_occ(:)
        logical, intent(in), optional :: tCoreSpaceDet, tLagrCorr

        real(dp) :: full_sign(spawn%rdm_send%sign_length), IterDetOcc_all(len_iter_occ_tot)
        integer(n_int) :: SpinCoupDet(0:nIfTot)
        integer :: nSpinCoup(nel), SignFac, HPHFExcitLevel
        integer :: IterLastRDMFill, AvSignIters, IterRDM
        integer :: AvSignIters_new(spawn%rdm_send%sign_length), IterRDM_new(spawn%rdm_send%sign_length)
        logical :: tUseDet, tInit, tLC
        integer :: run

        tUseDet = tNonInitsForRDMs
        if (.not. tUseDet) then
            tUseDet = all_runs_are_initiator(ilutnI)
        end if

        if (present(tLagrCorr)) then
            tLC = tLagrCorr
        else
            tLC = .true.
        end if

        if (tUseDet) then
            ! This is the number of iterations this determinant has been occupied,
            ! over *all* replicas.
            IterDetOcc_all = real(Iter + PreviousCycles, dp) - iter_occ + 1.0_dp

            ! IterLastRDMFill is the number of iterations from the last time the
            ! energy was calculated.
            IterLastRDMFill = mod((Iter + PreviousCycles - IterRDMStart + 1), RDMEnergyIter)

            AvSignIters_new = int(min(IterDetOcc_all(1::2), IterDetOcc_all(2::2)))

            ! The number of iterations we want to weight this RDM contribution by is:
            if (IterLastRDMFill > 0) then
                IterRDM_new = min(AvSignIters_new, IterLastRDMFill)
            else
                IterRDM_new = AvSignIters_new
            end if

            full_sign = 0.0_dp

            if (tHPHF) then
                if (.not. TestClosedShellDet(iLutnI)) then
                    if (RDMExcitLevel == 1) then
                        call fill_diag_1rdm(one_rdms, nI, av_sign / sqrt(2.0_dp), tCoreSpaceDet, IterRDM_new, tLC)
                    else
                        full_sign = IterRDM_new*av_sign(1::2)*av_sign(2::2) / 2.0_dp
                        call fill_spawn_rdm_diag(spawn, nI, full_sign)
                    end if

                    ! C_X D_X = C_X / sqrt(2) [ D_I +/- D_I'] - for open shell dets,
                    ! divide stored C_X by sqrt(2).
                    ! Add in I.
                    call FindExcitBitDetSym(iLutnI, SpinCoupDet)
                    call decode_bit_det(nSpinCoup, SpinCoupDet)
                    ! Find out if it's + or - in the above expression
                    SignFac = hphf_sign(iLutnI)

                    if (RDMExcitLevel == 1) then
                        call fill_diag_1rdm(one_rdms, nSpinCoup, real(SignFac, dp) * av_sign / sqrt(2.0_dp), &
                                            tCoreSpaceDet, IterRDM_new, tLC)
                    else
                        full_sign = IterRDM_new*av_sign(1::2)*av_sign(2::2) / 2.0_dp
                        call applyRDMCorrection()
                        call fill_spawn_rdm_diag(spawn, nI, full_sign)
                    end if

                    ! For HPHF we're considering < D_I + D_I' | a_a+ a_b+ a_j a_i | D_I + D_I' >
                    ! Not only do we have diagonal < D_I | a_a+ a_b+ a_j a_i | D_I > terms, but also cross terms
                    ! < D_I | a_a+ a_b+ a_j a_i | D_I' > if D_I and D_I' can be connected by a single or double
                    ! excitation. Find excitation level between D_I and D_I' and add in the contribution if connected.
                    HPHFExcitLevel = FindBitExcitLevel(iLutnI, SpinCoupDet, 2)
                    if (HPHFExcitLevel <= 2) then
                        call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nSpinCoup, &
                                                  IterRDM_new*av_sign(1::2)/sqrt(2.0_dp), &
                                                  (real(SignFac, dp)*av_sign(2::2)) / sqrt(2.0_dp))
                        call Add_RDM_From_IJ_Pair(spawn, one_rdms, nSpinCoup, nI, &
                                                  IterRDM_new*av_sign(2::2)/sqrt(2.0_dp), &
                                                  (real(SignFac, dp)*av_sign(1::2)) / sqrt(2.0_dp))
                    end if
                else

                    ! HPHFs on, but determinant closed shell.
                    if (RDMExcitLevel == 1) then
                        call fill_diag_1rdm(one_rdms, nI, av_sign, tCoreSpaceDet, &
                                            IterRDM_new, tLC)
                    else
                        full_sign = IterRDM_new*av_sign(1::2)*av_sign(2::2)
                        call applyRDMCorrection()
                        call fill_spawn_rdm_diag(spawn, nI, full_sign)
                    end if

                end if
                call Add_RDM_HFConnections_HPHF(spawn, one_rdms, iLutnI, nI, &
                                                av_sign, AvNoAtHF, ExcitLevelI, IterRDM_new)

            else if (tGUGA) then
                if (any(abs(av_sign(1::2)*av_sign(2::2)) > 1.0e-10_dp)) then
                    if (RDMExcitLevel == 1) then
                        call fill_diag_1rdm_guga(one_rdms, nI, av_sign)
                    else
                        full_sign = IterRDM_new*av_sign(1::2)*av_sign(2::2)
                        call applyRDMCorrection()
                        call fill_spawn_rdm_diag_guga(spawn, nI, full_sign)
                    end if
                end if

                if ((.not. all(near_zero(AvNoatHF))) .and. &
                    (.not. all(near_zero(av_sign)))) then
                    call Add_RDM_HFConnections_GUGA(spawn, one_rdms, ilutnI, av_sign, &
                                                    AvNoAtHF, ExcitLevelI, IterRDM_new)
                end if
            else
                ! Not using HPHFs.
                if (any(abs(av_sign(1::2)*av_sign(2::2)) > 1.0e-10_dp)) then
                    if (RDMExcitLevel == 1) then
                        call fill_diag_1rdm(one_rdms, nI, av_sign, tCoreSpaceDet, IterRDM_new)
                    else
                        full_sign = IterRDM_new*av_sign(1::2)*av_sign(2::2)
                        ! in adaptive shift mode, the reference contribution is rescaled
                        ! projEDet has to be the same on all runs
                        call applyRDMCorrection()
                        call fill_spawn_rdm_diag(spawn, nI, full_sign)
                    end if
                end if
                call Add_RDM_HFConnections_Norm(spawn, one_rdms, iLutnI, nI, av_sign, AvNoAtHF, ExcitLevelI, IterRDM_new)
            end if
        end if

    contains

        subroutine applyRDMCorrection()
            implicit none
            if (tAdaptiveShift .and. DetBitEq(ilutRef(:, 1), ilutnI) .and. &
                tNonInitsForRDMs .and. .not. tInitsRDMRef .and. tLC) &
                full_sign = full_sign + IterRDM_new * rdmCorrectionFactor
        end subroutine applyRDMCorrection

    end subroutine fill_rdm_diag_currdet_norm

    subroutine det_removed_fill_diag_rdm(spawn, one_rdms, iLutnI, av_sign, iter_occ, tLagrCorr)

        ! This routine is called if a determinant is removed from the list of
        ! currently occupied. At this point we need to add in its diagonal
        ! contribution for the number of iterations it has been occupied (or
        ! since the contribution was last included).

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: FindBitExcitLevel
        use CalcData, only: NMCyc
        use FciMCData, only: iLutHF_True, Iter, IterRDMStart, PreviousCycles
        use LoggingData, only: RDMEnergyIter
        use rdm_data, only: one_rdm_t
        use SystemData, only: nel, max_ex_level

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp), intent(in) :: av_sign(:), iter_occ(:)
        logical, intent(in), optional :: tLagrCorr

        integer :: nI(nel), ExcitLevel
        logical :: tLC

        if (present(tLagrCorr)) then
            tLC = tLagrCorr
        else
            tLC = .true.
        end if

        ! If the determinant is removed on an iteration that the diagonal
        ! RDM elements are  already being calculated, it will already have
        ! been counted. So check this isn't the case first.
        if (.not. ((Iter == NMCyc) .or. (mod((Iter + PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0))) then
            call decode_bit_det(nI, iLutnI)
            ExcitLevel = FindBitExcitLevel(iLutHF_True, iLutnI, max_ex_level)

            call fill_rdm_diag_currdet_norm(spawn, one_rdms, iLutnI, nI, ExcitLevel, av_sign, iter_occ, .false., tLC)
        end if

    end subroutine det_removed_fill_diag_rdm

    subroutine Add_RDM_HFConnections_Norm(spawn, one_rdms, iLutJ, nJ, AvSignJ, AvSignHF, walkExcitLevel, IterRDM)

        ! This is called when we run over all TotWalkers in CurrentDets.
        ! It is called for each CurrentDet which is a single or double of the HF.
        ! It explicitly adds in the HF - S/D connection, as if the HF were D_i and
        ! the single or double D_j. This is the standard full space RDM calc (No HPHF).
        ! In this case the diagonal elements wll already be taken care of.

        use FciMCData, only: HFDet_True
        use rdm_data, only: one_rdm_t
        use SystemData, only: nel, t_3_body_excits

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(nel)
        real(dp), intent(in) :: AvSignJ(:), AvSignHF(:)
        integer, intent(in) :: walkExcitLevel
        integer, intent(in) :: IterRDM(:)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "Add_RDM_HFConnections_Norm"
#endif

        integer(n_int) :: iUnused

        ! If we have a single or double, add in the connection to the HF,
        ! symmetrically.
        ASSERT(.not. t_3_body_excits)
        if ((walkExcitLevel == 1) .or. (walkExcitLevel == 2)) then
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, HFDet_True, nJ, AvSignHF(2::2), IterRDM*AvSignJ(1::2))
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, nJ, HFDet_True, AvSignJ(2::2), IterRDM*AvSignHF(1::2))
        end if

        ! Eliminate compiler warnings.
        iUnused = ilutJ(0)

    end subroutine Add_RDM_HFConnections_Norm

    subroutine Add_RDM_HFConnections_HPHF(spawn, one_rdms, iLutJ, nJ, AvSignJ, AvSignHF, walkExcitLevel, IterRDM)

        ! This is called when we run over all TotWalkers in CurrentDets.
        ! It is called for each CurrentDet which is a single or double of the HF.
        ! It adds in the HF - S/D connection. The diagonal elements will already
        ! have been taken care of by the extract routine.

        use FciMCData, only: HFDet_True, iLutHF_True
        use rdm_data, only: one_rdm_t
        use SystemData, only: nel, t_3_body_excits

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(nel)
        real(dp), intent(in) :: AvSignJ(:), AvSignHF(:)
        integer, intent(in) :: walkExcitLevel
        integer, intent(in) :: IterRDM(:)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "Add_RDM_HFConnections_HPHF"
#endif

        ! Now if the determinant is connected to the HF (i.e. single or double),
        ! add in the diagonal elements of this connection as well -
        ! symmetrically because no probabilities are involved.
        ASSERT(.not. t_3_body_excits)
        if ((walkExcitLevel == 1) .or. (walkExcitLevel == 2)) then
            call Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutHF_True, iLutJ, HFDet_True, nJ, AvSignHF(2::2), IterRDM*AvSignJ(1::2))

            call Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutJ, iLutHF_True, nJ, HFDet_True, AvSignJ(2::2), IterRDM*AvSignHF(1::2))
        end if

    end subroutine Add_RDM_HFConnections_HPHF

    subroutine check_fillRDM_DiDj(rdm_defs, spawn, one_rdms, Spawned_No, iLutJ, realSignJ, &
                                  tNonInits)

        ! The spawned parts contain the Dj's spawned by the Di's in CurrentDets.
        ! If the SpawnedPart is found in the CurrentDets list, it means that
        ! the Dj has a non-zero cj - and therefore the Di.Dj pair will have a
        ! non-zero ci.cj to contribute to the RDM. The index i tells us where
        ! to look in the parent array, for the Di's to go with this Dj.

        use FciMCData, only: iLutHF_True
        use rdm_data, only: one_rdm_t, rdm_definitions_t

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: Spawned_No
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), intent(in) :: realSignJ(lenof_sign)
        logical, intent(in), optional :: tNonInits
        logical :: tAllContribs

        ! optionally only sum in initiator contributions
        if (present(tNonInits)) then
            tAllContribs = tNonInits
        else
            tAllContribs = .true.
        end if

        if (.not. DetBitEQ(iLutHF_True, iLutJ, nifd)) then
            call DiDj_Found_FillRDM(rdm_defs, spawn, one_rdms, Spawned_No, &
                                    iLutJ, realSignJ, tAllContribs)
        end if

    end subroutine check_fillRDM_DiDj

    subroutine DiDj_Found_FillRDM(rdm_defs, spawn, one_rdms, Spawned_No, iLutJ, real_sign_j_all, &
                                  tNonInits)

        ! This routine is called when we have found a Di (or multiple Di's)
        ! spawning onto a Dj with sign /= 0 (i.e. occupied). We then want to
        ! run through all the Di, Dj pairs and add their coefficients
        ! (with appropriate de-biasing factors) into the 1 and 2 electron RDM.

        use bit_reps, only: decode_bit_det
        use FciMCData, only: Spawned_Parents, Spawned_Parents_Index, iLutHF_True
        use rdm_data, only: one_rdm_t, rdm_definitions_t
        use SystemData, only: nel, tHPHF
        use bit_reps, only: all_runs_are_initiator

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: Spawned_No
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), intent(in) :: real_sign_j_all(lenof_sign)
        logical, intent(in) :: tNonInits

        integer :: i, irdm, rdm_ind, nI(nel), nJ(nel)
        real(dp) :: realSignI
        real(dp) :: input_sign_i(rdm_defs%nrdms), input_sign_j(rdm_defs%nrdms)
        integer :: dest_part_type, source_part_type, run
        integer(n_int) :: source_flags
        logical :: spawning_from_ket_to_bra, tNonInitParent
        integer(int_rdm) :: guga_rdm_ind
        real(dp) :: x0, x1

        ! Spawning from multiple parents, to iLutJ, which has SignJ.

        ! We are at position Spawned_No in the SpawnedParts array.
        ! Spawned_Parents_Index(1,Spawned_No) is therefore the start position
        ! of the list of parents (Di's) which spawned on the Dj in
        ! SpawnedParts(Spawned_No). There are Spawned_Parents_Index(2,Spawned_No)
        ! of these parent Di's. Spawned_Parents(0:IlutBitsParent%len_orb,x)
        ! is the determinant Di, Spawned_Parents(IlutBitsParent%ind_pop,x)
        ! is the un-biased ci. Spawned_Parents(IlutBitsParent%ind_flag,x)
        ! is the flag (only initiator i believe. WD) and
        ! Spawned_Parents(IlutBitsParent%ind_source,x) is the simulation index,
        ! where the spawn was coming from

        ! Run through all Di's.

        ! in the GUGA RDMs it somehow can happen that a all 0 real_sign_j_all
        ! gets in here.. to abort if that happens
        if (all(near_zero(real_sign_j_all))) return

        do i = Spawned_Parents_Index(1, Spawned_No), &
            Spawned_Parents_Index(1, Spawned_No) + Spawned_Parents_Index(2, Spawned_No) - 1

            if (DetBitEQ(iLutHF_True, &
                         Spawned_Parents(0:IlutBitsParent%len_orb, i), IlutBits%len_orb)) then
                ! We've already added HF - S, and HF - D symmetrically.
                ! Any connection with the HF has therefore already been added.
                cycle
            end if

            call decode_bit_det(nI, Spawned_Parents(0:IlutBits%len_orb, i))
            call decode_bit_det(nJ, iLutJ)

            realSignI = transfer(Spawned_Parents(IlutBits%ind_pop, i), realSignI)

            ! if the sign is 0, there is nothing to do, i.e. all contributions will be 0
            if (abs(realSignI) < eps) cycle

            ! The original spawning event (and the RealSignI) came from this
            ! population.
            source_part_type = int(Spawned_Parents(IlutBitsParent%ind_source, i))

            ! WD: due to some weird convention in Annihilation (TODO check that!)
            ! it can happen that a source_part_type = 0 can end up here..
            ! so cycle if that happens
            if (source_part_type == 0) cycle

            ! if we only sum in initiator contriubtions, check the flags here
            ! This block is entered if
            ! a) tNonInits == .false. -> we are looking at the inits-rdms
            ! b) tNonInitsForRDMs == .false. -> we are calculating inits-only rdms
            ! In both cases, an initiator check is required
            if (.not. (tNonInits .and. tNonInitsForRDMs)) then
                tNonInitParent = .false.
                ! Only initiators are used
                if (.not. (all_runs_are_initiator(ilutJ) .or. &
                           ! or the non-variational inits-only version (only require an init in the ket)
                           (tNonVariationalRDMs .and. .not. tNonInits))) return
                do run = 1, inum_runs
                    if (.not. btest(Spawned_Parents(IlutBitsParent%ind_flag, i), &
                                    get_initiator_flag_by_run(run))) tNonInitParent = .true.
                    ! if a non-initiator is participating, do not sum in
                    ! that contribution
                end do
                if (tNonInitParent) cycle
            end if

            ! Loop over all RDMs to which the simulation with label
            ! source_part_type contributes to.
            do irdm = 1, rdm_defs%nrdms_per_sim(source_part_type)
                ! Get the label of the simulation that is paired with this,
                ! replica, for this particular RDM.
                dest_part_type = rdm_defs%sim_pairs(irdm, source_part_type)

                ! The label of the RDM that this is contributing to.
                rdm_ind = rdm_defs%rdm_labels(irdm, source_part_type)

                input_sign_i = 0.0_dp
                input_sign_j = 0.0_dp
                input_sign_i(rdm_ind) = realSignI
                input_sign_j(rdm_ind) = real_sign_j_all(dest_part_type)

                ! If the two FCIQMC simulations involved are different then we
                ! need a factor of a half, since spawning can occur from
                ! from either of the two simulations.
                if (dest_part_type /= source_part_type) input_sign_i = 0.5_dp * input_sign_i

                ! sim_labels(1,irdm) holds the state in the 'bra' of the RDM.
                ! sim_labels(2,irdm) holds the state in the 'ket' of the RDM.
                ! If source_part_type is equal to the part type of the ket,
                ! then we can say that we are spawning from the ket to the bra.
                spawning_from_ket_to_bra = source_part_type == rdm_defs%sim_labels(2, rdm_ind)

                ! If we are spawning from the ket of the RDM, then the RDM
                ! element to be added in is:
                !
                ! \psi_j^b* \psi_i^a < nJ | ... | nI > (1)
                !
                ! (where a and b are state labels), which is what we want.
                !
                ! If we are spawning from the bra to the ket then we cannot
                ! just pass input_sign_i, input_sign_j, nI and nJ into the
                ! following routines in the same order, because then we would
                ! be adding in the same element as in Eq. (1), but we actually
                ! want the Hermitian conjugate of that term - <nJ| and |nI>
                ! will be in the wrong order, and so the RDM indices will
                ! be calculated the wrong way around, and the element will
                ! be added in on the 'wrong side' of the RDM's diagonal.
                ! As such, we instead swap the order of things so that we are
                ! adding in the term
                !
                ! \psi_i^a* \psi_j^b < nI | ... | nJ > (2)
                !
                ! which is clearly as it should be if the determinant we are
                ! spawning from (nI) is the bra of the RDM.

                if (tGUGA) then
                    call extract_stochastic_rdm_info(IlutBitsParent, &
                             Spawned_Parents(:, i), guga_rdm_ind, x0, x1)
                end if

                if (spawning_from_ket_to_bra) then
                    if (tHPHF) then
                        call Fill_Spin_Coupled_RDM(spawn, one_rdms, &
                                                   Spawned_Parents(0:IlutBitsParent%len_orb, i), &
                                                   iLutJ, nI, nJ, input_sign_i, input_sign_j)

                    else if (tGUGA) then
                        call Add_RDM_From_IJ_Pair_GUGA(spawn, one_rdms, &
                           nI, nJ, input_sign_i, input_sign_j, spawning_from_ket_to_bra, &
                           rdm_ind_in=guga_rdm_ind, x0=x0, x1=x1)

                    else
                        call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, &
                                                  input_sign_i, input_sign_j)
                    end if
                else
                    ! Spawning from the bra to the ket - swap the order of the
                    ! signs and states in the routine interface.
                    if (tHPHF) then
                        call Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutJ, &
                                                   Spawned_Parents(0:IlutBitsParent%len_orb, i), &
                                                   nJ, nI, input_sign_j, input_sign_i)

                    else if (tGUGA) then
                        call Add_RDM_From_IJ_Pair_GUGA(spawn, one_rdms, &
                           nJ, nI, input_sign_j, input_sign_i, spawning_from_ket_to_bra, &
                           rdm_ind_in=guga_rdm_ind, x0=x0, x1=x1)

                    else
                        call Add_RDM_From_IJ_Pair(spawn, one_rdms, nJ, nI, &
                                                  input_sign_j, input_sign_i)
                    end if
                end if
            end do

        end do

    end subroutine DiDj_Found_FillRDM

    subroutine Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutnI, iLutnJ, nI, nJ, realSignI, realSignJ)

        ! It takes to HPHF functions, and calculates what needs to be summed
        ! into the RDMs.

        ! If the two HPHF determinants we're considering consist of I + I' and
        ! J + J', where X' is the spin coupled (all spins flipped) version of X,
        ! then we have already considered the I -> J excitation. And if I and
        ! J are connected by a double excitation, tDoubleConnection is true
        ! and we have also considered I' -> J'. But we need to also account
        ! for I -> J' and I' -> J.

        use DetBitOps, only: TestClosedShellDet, FindBitExcitLevel
        use hphf_integrals, only: hphf_sign
        use HPHFRandExcitMod, only: FindExcitBitDetSym, FindDetSpinSym
        use rdm_data, only: one_rdm_t
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutnI(0:NIfTot), iLutnJ(0:NIfTot)
        real(dp), intent(in) :: realSignI(:), realSignJ(:)
        integer, intent(in) :: nI(nel), nJ(nel)

        integer(n_int) :: iLutnI2(0:NIfTot)
        integer :: nI2(nel), nJ2(nel)
        real(dp) :: NewSignJ(size(realSignJ)), NewSignI(size(realSignI))
        real(dp) :: PermSignJ(size(realSignJ)), PermSignI(size(realSignI))
        integer :: I_J_ExcLevel, ICoup_J_ExcLevel

        if (TestClosedShellDet(iLutnI)) then

            if (TestClosedShellDet(iLutnJ)) then
                ! Closed shell -> Closed shell - just as in determinant case
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, realSignI, realSignJ)
            else
                ! Closed shell -> open shell.
                call FindDetSpinSym(nJ, nJ2, nel)
                NewSignJ = realSignJ / Root2
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, realSignI, NewSignJ)
                ! What is the permutation between Di and Dj'
                NewSignJ = NewSignJ * hphf_sign(iLutnJ)
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ2, realSignI, NewSignJ)
            end if

        else if (TestClosedShellDet(iLutnJ)) then
            ! Open shell -> closed shell
            call FindDetSpinSym(nI, nI2, nel)
            NewSignI = realSignI / Root2
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, NewSignI, realSignJ)
            ! What is the permutation between Di' and Dj?
            NewSignI = NewSignI * hphf_sign(iLutnI)
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI2, nJ, NewSignI, realSignJ)

        else
            ! Open shell -> open shell
            NewSignI = realSignI / Root2
            NewSignJ = realSignJ / Root2
            PermSignJ = NewSignJ * real(hphf_sign(iLutnJ), dp)
            PermSignI = NewSignI * real(hphf_sign(iLutnI), dp)
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            call FindDetSpinSym(nI, nI2, nel)
            call FindDetSpinSym(nJ, nJ2, nel)
            I_J_ExcLevel = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            ICoup_J_ExcLevel = FindBitExcitLevel(iLutnI2, iLutnJ, 2)

            if (I_J_ExcLevel <= 2) then
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, NewSignI, NewSignJ)
                ! Di -> Dj
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI2, nJ2, PermSignI, PermSignJ)
                ! Di' -> Dj'  (both permuted sign)
            end if

            if (ICoup_J_ExcLevel <= 2) then
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI2, nJ, PermSignI, NewSignJ)
                ! Di' -> Dj  (i permuted sign)
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ2, NewSignI, PermSignJ)
                ! Di  -> Dj'  (j permuted sign)
            end if
        end if

    end subroutine Fill_Spin_Coupled_RDM

    subroutine Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, realSignI, realSignJ)

        ! This routine takes a pair of different determinants Di and Dj, and
        ! figures out which type of elements need to be added in to the RDM.

        ! ------------------------ *IMPORTANT* --------------------------------
        ! nI refers to the ket determinant and nJ refers to the bra,
        ! determinant. Similarly, realSignI must be the sign corresponding to
        ! the ket determinant, and realSignJ to the bra, i.e. the element
        ! added in is equal to:
        !   realSignJ^* realSignI < nJ | ... | nI >
        ! Getting nI and nJ the wrong way around will result in RDM elements
        ! being added to the wrong side of the diagonal, which is a problem for
        ! transition RDMs which are *NOT* hermitian, so this is important.
        ! Similarly, getting nI and nJ right but realSignI and realSignJ the
        ! wrong way around will cause problems for complex RDMs, since the
        ! complex conjugate will be taken of the wrong sign.
        ! ---------------------------------------------------------------------

        use LoggingData, only: RDMExcitLevel
        use rdm_data, only: one_rdm_t
        use rdm_data_utils, only: add_to_rdm_spawn_t
        use SystemData, only: nel, t_3_body_excits

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(in) :: realSignI(:), realSignJ(:)

        integer :: Ex(2, maxExcit), Ex_symm(2, maxExcit)
        logical :: tParity
        real(dp) :: full_sign(spawn%rdm_send%sign_length)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "Add_RDM_From_IJ_Pair"
#endif
        Ex(:, :) = 0
        ! Maximum excitation level - we know they are connected by a double
        ! or a single excitation.
        Ex(1, 1) = 2
        tParity = .false.
        ASSERT(.not. t_3_body_excits)

        ! Ex(1,:) comes out as the orbital(s) excited from, i.e. i,j.
        ! Ex(2,:) comes out as the orbital(s) excited to, i.e. a,b.
        call GetExcitation(nI, nJ, nel, Ex, tParity)

        full_sign = 0.0_dp
        if (tParity) then
            full_sign = -realSignI * realSignJ
        else
            full_sign = realSignI * realSignJ
        end if

        ASSERT(all(ex(:, 1) > 0) .and. all(ex(:, 1) <= nbasis))

        if ((Ex(1, 2) == 0) .and. (Ex(2, 2) == 0)) then

            ! Di and Dj are separated by a single excitation.
            ! Add in the contribution from this pair into the 1-RDM.

            if (RDMExcitLevel == 1) then
                call fill_sings_1rdm(one_rdms, Ex, tParity, realSignI, realSignJ, .false.)
            else
                call fill_spawn_rdm_singles(spawn, nI, Ex, full_sign)
            end if

        else if (RDMExcitLevel /= 1) then

            ! Otherwise Di and Dj are connected by a double excitation.
            ! Add in this contribution to the 2-RDM, if being calculated.
            call add_to_rdm_spawn_t(spawn, Ex(2, 1), Ex(2, 2), Ex(1, 1), Ex(1, 2), full_sign, .false.)
        end if

    end subroutine Add_RDM_From_IJ_Pair

    subroutine fill_spawn_rdm_diag(spawn, nI, full_sign)

        use rdm_data_utils, only: add_to_rdm_spawn_t
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: full_sign(spawn%rdm_send%sign_length)

        integer :: iel, jel

        ! Looking at elements of the type Gamma(i,j,i,j).
        do iel = 1, nel - 1
            do jel = iel + 1, nel
                associate(i => nI(iel), j => nI(jel))
                    call add_to_rdm_spawn_t(spawn, i, j, i, j, full_sign, .false.)
                end associate
            end do
        end do

    end subroutine fill_spawn_rdm_diag

    subroutine fill_spawn_rdm_singles(spawn, nI, Ex, full_sign)

        use rdm_data_utils, only: add_to_rdm_spawn_t
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nI(nel), Ex(2, maxExcit)
        real(dp), intent(in) :: full_sign(spawn%rdm_send%sign_length)

        integer :: iel

        ! Looking at elements of the type Gamma(a,k,i,k).
        do iel = 1, nel
            associate(k => nI(iel))
                if (k < Ex(1, 1) .and. k < Ex(2, 1)) then
                    call add_to_rdm_spawn_t(spawn, k, Ex(2, 1), k, Ex(1, 1), full_sign, .false.)
                else if (k < Ex(1, 1) .and. k > Ex(2, 1)) then
                    call add_to_rdm_spawn_t(spawn, Ex(2, 1), k, k, Ex(1, 1), -full_sign, .false.)
                else if (k > Ex(1, 1) .and. k < Ex(2, 1)) then
                    call add_to_rdm_spawn_t(spawn, k, Ex(2, 1), Ex(1, 1), k, -full_sign, .false.)
                else if (k > Ex(1, 1) .and. k > Ex(2, 1)) then
                    call add_to_rdm_spawn_t(spawn, Ex(2, 1), k, Ex(1, 1), k, full_sign, .false.)
                end if
            end associate
        end do

    end subroutine fill_spawn_rdm_singles

! =======================================================================================
! THESE NEXT ROUTINES ARE GENERAL TO BOTH STOCHASTIC AND EXPLICIT
! =======================================================================================

    subroutine fill_diag_1rdm(one_rdms, nI, contrib_sign, tCoreSpaceDetIn, RDMItersIn, tLagrCorr)

        ! Add the contribution to the diagonal elements of the 1-RDM from
        ! determinant nI, with the corresponding sign(s) from an FCIQMC
        ! simulation gievn by contrib_sign.

        ! The (spinned) 1-RDM is defined by
        !
        ! \gamma_{i,j} = < \Psi | a^{\dagger}_i a_j | \Psi >.
        !
        ! If we have a closed shell system then we instead add to the spin-free
        ! 1-RDM, obtained by tracing over the spin component.

        use rdm_data, only: one_rdm_t, tOpenShell
        use LoggingData, only: ThreshOccRDM, tThreshOccRDMDiag, tExplicitAllRDM
        use RotateOrbsData, only: SymLabelListInv_rot
        use UMatCache, only: gtID

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(:)
        real(dp), intent(in) :: contrib_sign(:)
        logical, optional, intent(in) :: tCoreSpaceDetIn
        integer, optional, intent(in) :: RDMItersIn(:)
        logical, optional, intent(in) :: tLagrCorr

        integer :: i, ind, irdm
        integer :: RDMIters(size(one_rdms))
        real(dp) :: ScaleContribFac, final_contrib(size(one_rdms))
        logical :: tCoreSpaceDet, tLC

        if (present(tLagrCorr)) then
            tLC = tLagrCorr
        else
            tLC = .true.
        end if

        ScaleContribFac = 1.0_dp

        RDMIters = 1
        if (present(RDMItersIn)) RDMIters = RDMItersIn

        tCoreSpaceDet = .false.
        if (present(tCoreSpaceDetIn)) tCoreSpaceDet = tCoreSpaceDetIn

        ! This is the single-run cutoff being applied (do not use in DR mode):
        if (.not. tCoreSpaceDet) then
            ! Dets in the core space are never removed from main list, so
            ! strictly do not require corrections.
            if (tThreshOccRDMDiag .and. (abs(contrib_sign(1)) <= ThreshOccRDM)) ScaleContribFac = 0.0_dp
        end if

        do i = 1, size(nI)
            ! The SymLabelListInv_rot array is used to index the 1-RDM so that
            ! it will be filled in block diagonal order, where each block holds
            ! one symmetry sector of the 1-RDM (i.e., orbitals are ordered
            ! first by symmetry, then by the standard order).
            if (tOpenShell) then
                ind = SymLabelListInv_rot(nI(i))
            else
                ind = SymLabelListInv_rot(gtID(nI(i)))
            end if
            if (size(contrib_sign) == 1) then
                final_contrib = contrib_sign**2 * RDMIters * ScaleContribFac
            else
                final_contrib = contrib_sign(1::2)*contrib_sign(2::2) * RDMIters * ScaleContribFac
            end if
            ! in adaptive shift mode, the reference contribution is rescaled
            ! we assume that projEDet is the same on all runs, else there is no point
            if (tAdaptiveShift .and. all(nI == projEDet(:, 1)) &
                .and. tNonInitsForRDMs .and. .not. tInitsRDMRef .and. tLC) &
                final_contrib = final_contrib + RDMIters * ScaleContribFac * rdmCorrectionFactor
            do irdm = 1, size(one_rdms)
                one_rdms(irdm)%matrix(ind, ind) = one_rdms(irdm)%matrix(ind, ind) + final_contrib(irdm)
            end do
        end do

    end subroutine fill_diag_1rdm

    subroutine fill_sings_1rdm(one_rdms, Ex, tParity, contrib_sign_i, contrib_sign_j, fill_symmetric)

        use rdm_data, only: one_rdm_t, tOpenShell
        use RotateOrbsData, only: SymLabelListInv_rot
        use UMatCache, only: gtID

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: Ex(2, maxExcit)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: contrib_sign_i(:), contrib_sign_j(:)
        logical, intent(in) :: fill_symmetric

        integer :: i, a, ind_i, ind_a, irdm
        real(dp) :: ParityFactor

        ParityFactor = 1.0_dp
        if (tParity) ParityFactor = -1.0_dp

        ! Get the orbital labels involved in the excitation.
        if (tOpenShell) then
            i = Ex(1, 1)
            a = Ex(2, 1)
        else
            i = gtID(Ex(1, 1))
            a = gtID(Ex(2, 1))
        end if

        ! The SymLabelListInv_rot array is used to index the 1-RDM so that it
        ! will be filled in block diagonal order, where each block holds one
        ! symmetry sector of the 1-RDM (i.e., orbitals are ordered first by
        ! symmetry, then by the standard order).
        ind_i = SymLabelListInv_rot(i)
        ind_a = SymLabelListInv_rot(a)

        do irdm = 1, size(one_rdms)
            one_rdms(irdm)%matrix(ind_a, ind_i) = one_rdms(irdm)%matrix(ind_a, ind_i) + &
                                                  (ParityFactor * contrib_sign_i(irdm) * contrib_sign_j(irdm))

            if (fill_symmetric) then
                one_rdms(irdm)%matrix(ind_i, ind_a) = one_rdms(irdm)%matrix(ind_i, ind_a) + &
                                                      (ParityFactor * contrib_sign_i(irdm) * contrib_sign_j(irdm))
            end if
        end do

    end subroutine fill_sings_1rdm

    subroutine fill_RDM_offdiag_deterministic(rdm_defs, spawn, one_rdms)

        use bit_rep_data, only: NIfD
        use bit_reps, only: decode_bit_det
        use DetBitOps, only: get_bit_excitmat, FindBitExcitLevel
        use FciMCData, only: Iter, IterRDMStart, PreviousCycles, iLutHF_True, core_run
        use core_space_util, only: cs_replicas
        use LoggingData, only: RDMExcitLevel, RDMEnergyIter
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: one_rdm_t, rdm_definitions_t
        use rdm_data_utils, only: add_to_rdm_spawn_t
        use SystemData, only: nel, tHPHF

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)

        integer :: i, j, irdm
        integer :: SingEx(2, 1), Ex(2, maxExcit)
        real(dp) :: AvSignI(spawn%rdm_send%sign_length), AvSignJ(spawn%rdm_send%sign_length)
        real(dp) :: full_sign(spawn%rdm_send%sign_length)
        logical :: tParity
        integer(n_int) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer :: nI(nel), nJ(nel), IC, n
        integer :: IterRDM, connect_elem, num_j
        type(ExcitationInformation_t) :: excitInfo
        character(*), parameter :: this_routine = "fill_rdm_offdiag_deterministic"

        ! IterRDM will be the number of iterations that the contributions are
        ! ech weighted by.
        if (mod((iter + PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0) then
            ! This numer of iterations is how regularly the energy is printed
            ! out.
            IterRDM = RDMEnergyIter
        else
            ! This must be the final iteration, as we've got tFill_RDM=.true.
            ! for an iteration where we wouldn't normally need the energy
            IterRDM = mod((Iter + PreviousCycles - IterRDMStart + 1), RDMEnergyIter)
        end if

        Ex = 0
        associate(rep => cs_replicas(core_run))
            if (t_full_core_rdms) then
                ! num_j = rep%determ_sizes(iProcIndex)
                num_j = rep%determ_space_size
            end if

            ! Cycle over all core dets on this process.
            do i = 1, rep%determ_sizes(iProcIndex)
                if (.not. t_full_core_rdms) then
                    num_j = rep%sparse_core_ham(i)%num_elements - 1
                end if
                iLutI = rep%core_space(:, rep%determ_displs(iProcIndex) + i)

                ! Connections to the HF are added in elsewhere, so skip them here.
                if (DetBitEq(iLutI, iLutHF_True, nifd)) cycle

                do irdm = 1, rdm_defs%nrdms
                    AvSignI(irdm) = rep%full_determ_vecs_av(rdm_defs%sim_labels(2, irdm), &
                        rep%determ_displs(iProcIndex) + i)
                end do

                call decode_bit_det(nI, iLutI)

                if (tGUGA) call init_csf_information(ilutI(0:nifd))

                do j = 1, num_j
                    ! Running over all non-zero off-diag matrix elements
                    ! Connections to whole space (1 row), excluding diagonal elements

                    ! Note: determ_displs holds sum(determ_sizes(0:proc-1))
                    ! Core space holds all the core determinants on every processor,
                    ! so we need to shuffle up to the range of indices corresponding
                    ! to this proc (using determ_displs) and then select the
                    ! correct one, i.

                    if (t_full_core_rdms) then
                        iLutJ = rep%core_space(:, j)
                    else
                        iLutJ = rep%core_space(:, rep%core_connections(i)%positions(j))
                    end if

                    ! Connections to the HF are added in elsewhere, so skip them here.
                    if (DetBitEq(iLutJ, iLutHF_True, nifd)) cycle
                    if (DetBitEq(iLutJ, iLutI, nifd)) cycle


                    if (t_full_core_rdms) then
                        if (.not. tGUGA) then
                            ic = FindBitExcitLevel(iLutI, ilutJ)
                            call GetBitExcitation(ilutI, ilutJ, ex, tParity)
                        end if

                        do irdm = 1, rdm_defs%nrdms
                            AvSignJ(irdm) = rep%full_determ_vecs_av(rdm_defs%sim_labels(1, irdm), j)
                        end do
                    else
                        do irdm = 1, rdm_defs%nrdms
                            AvSignJ(irdm) = rep%full_determ_vecs_av(rdm_defs%sim_labels(1, irdm), &
                                                rep%core_connections(i)%positions(j))
                        end do
                        connect_elem = rep%core_connections(i)%elements(j)
                        IC = abs(connect_elem)
                        if (sign(1, connect_elem) > 0) then
                            tParity = .false.
                        else
                            tParity = .true.
                        end if
                    end if

                    ! if (ic > 2) cycle

                    if (tParity) then
                        full_sign = -AvSignI * AvSignJ * IterRDM
                    else
                        full_sign = AvSignI * AvSignJ * IterRDM
                    end if

                    if (tHPHF) then
                        call decode_bit_det(nJ, iLutJ)
                        call Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutI, iLutJ, &
                            nI, nJ, AvSignI * IterRDM, AvSignJ)

                    else if (tGUGA) then

                        call add_rdm_from_ij_pair_guga_exact(spawn, one_rdms, &
                                  ilutI, ilutJ, AvSignI * IterRDM, AvSignJ, calc_type=1)
                    else
                        if (IC == 1) then
                            ! Single excitation - contributes to 1- and 2-RDM
                            ! (if calculated).

                            ! Note: get_bit_excitmat may be buggy (DetBitOps),
                            ! but will do for now as we need the Ex...
                            call get_bit_excitmat(iLutI(0:NIfD), iLutJ(0:NIfD), SingEx, IC)
                            Ex(:, 1) = SingEx(:, 1)

                            ! No need to explicitly fill symmetrically as we'll
                            ! generate pairs of determinants both ways around using
                            ! the connectivity matrix.
                            if (RDMExcitLevel /= 1) call fill_spawn_rdm_singles(spawn, nI, Ex, full_sign)

                            if (RDMExcitLevel == 1) then
                                call fill_sings_1rdm(one_rdms, Ex, tParity, AvSignI * IterRDM, AvSignJ, .false.)
                            end if

                        else if ((IC == 2) .and. (RDMExcitLevel /= 1)) then
                            ! Note: get_bit_excitmat may be buggy (DetBitOps),
                            ! but will do for now as we need the Ex...
                            call get_bit_excitmat(iLutI(0:NIfD), iLutJ(0:NIfD), Ex, IC)

                            call add_to_rdm_spawn_t(spawn, Ex(2, 1), Ex(2, 2), Ex(1, 1), Ex(1, 2), full_sign, .false.)
                        end if
                    end if
                end do
            end do
        end associate

    end subroutine fill_RDM_offdiag_deterministic

end module rdm_filling
