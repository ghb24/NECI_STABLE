#include "macros.h"

module rdm_filling

    ! This module contains routines used to perform filling of the RDM arrays,
    ! as done on-the-fly during an FCIQMC simulation.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use rdm_data, only: rdm_spawn_t

    implicit none

contains

    subroutine fill_rdm_diag_wrapper(spawn, one_rdms, ilut_list, ndets)

        ! Loop over all states in ilut_list and see if any signs have just
        ! become unoccupied or become reoccupied. In which case, we have
        ! started a new averaging block, so we need to add in the
        ! contributions from the last block to the corresponding RDMs.

        use bit_rep_data, only: extract_sign
        use CalcData, only: tPairedReplicas
        use global_det_data, only: get_iter_occ_tot, get_av_sgn_tot
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: one_rdm_t, nrdms, signs_for_rdm

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: ilut_list(:,:)
        integer, intent(in) :: ndets

        integer :: idet, irdm, av_ind_1, av_ind_2
        real(dp) :: curr_sign(lenof_sign), adapted_sign(len_av_sgn_tot)
        real(dp) :: av_sign(len_av_sgn_tot), iter_occ(len_iter_occ_tot)

        associate(ind => signs_for_rdm)

            do idet = 1, ndets

                call extract_sign(ilut_list(:,idet), curr_sign)
                ! All average sign from all RDMs.
                av_sign = get_av_sgn_tot(idet)
                ! The iteration on which each replica became occupied.
                iter_occ = get_iter_occ_tot(idet)

                adapted_sign = 0.0_dp

                if (tPairedReplicas) then
                    do irdm = 1, nrdms

                        ! The indicies of the first and second replicas in this
                        ! particular pair, in the *average* sign arrays (and
                        ! therefore also for the iter_occ array).
                        av_ind_1 = irdm*2-1
                        av_ind_2 = irdm*2

                        if ((abs(curr_sign(ind(1,irdm))) < 1.0e-10_dp .and. abs(iter_occ(av_ind_1)) > 1.0e-10_dp) .or. &
                            (abs(curr_sign(ind(2,irdm))) < 1.0e-10_dp .and. abs(iter_occ(av_ind_2)) > 1.0e-10_dp) .or. &
                            (abs(curr_sign(ind(1,irdm))) > 1.0e-10_dp .and. abs(iter_occ(av_ind_1)) < 1.0e-10_dp) .or. &
                            (abs(curr_sign(ind(2,irdm))) > 1.0e-10_dp .and. abs(iter_occ(av_ind_2)) < 1.0e-10_dp)) then 

                            ! In this case we want to include this diagonal element,
                            ! so transfer the sign.
                            adapted_sign(av_ind_1:av_ind_2) = av_sign(av_ind_1:av_ind_2)
                        end if
                    end do

                else

                    do irdm = 1, nrdms
                        if (abs(curr_sign(ind(1,irdm))) < 1.0e-10_dp) then
                            ! If this RDM sign has gone to zero, then we want to add
                            ! the contribution for this RDM.
                            adapted_sign(irdm) = av_sign(irdm)
                        end if
                    end do

                end if

                ! At least one of the signs has just gone to zero or just become
                ! reoccupied, so we need to add in diagonal elements and connections to HF
                if (any(abs(adapted_sign) > 1.e-12_dp)) then
                    call det_removed_fill_diag_rdm(spawn, one_rdms, ilut_list(:,idet), adapted_sign, iter_occ)
                end if

            end do

        end associate

    end subroutine fill_rdm_diag_wrapper

    subroutine fill_rdm_diag_currdet_norm(spawn, one_rdms, iLutnI, nI, ExcitLevelI, av_sign, iter_occ, tCoreSpaceDet)

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

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: TestClosedShellDet, FindBitExcitLevel
        use FciMCData, only: Iter, IterRDMStart, PreviousCycles, AvNoAtHF
        use global_det_data, only: len_iter_occ_tot
        use hphf_integrals, only: hphf_sign
        use HPHFRandExcitMod, only: FindExcitBitDetSym
        use LoggingData, only: RDMEnergyIter, RDMExcitLevel
        use rdm_data, only: one_rdm_t
        use SystemData, only: nel, tHPHF

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: nI(nel), ExcitLevelI
        real(dp), intent(in) :: av_sign(:), iter_occ(:)
        logical, intent(in), optional :: tCoreSpaceDet

        real(dp) :: full_sign(spawn%rdm_send%sign_length), IterDetOcc_all(len_iter_occ_tot)
        integer(n_int) :: SpinCoupDet(0:nIfTot)
        integer :: nSpinCoup(nel), SignFac, HPHFExcitLevel
        integer :: IterLastRDMFill, AvSignIters, IterRDM
        integer :: AvSignIters_new(spawn%rdm_send%sign_length), IterRDM_new(spawn%rdm_send%sign_length)

        ! This is the number of iterations this determinant has been occupied,
        ! over *all* replicas.
        IterDetOcc_all = real(Iter+PreviousCycles,dp) - iter_occ + 1.0_dp

        ! IterLastRDMFill is the number of iterations from the last time the
        ! energy was calculated.
        IterLastRDMFill = mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter)

        AvSignIters_new = int(min(IterDetOcc_all(1::nreplicas), IterDetOcc_all(nreplicas::nreplicas)))

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
                    call fill_diag_1rdm(one_rdms, nI, av_sign/sqrt(2.0_dp), tCoreSpaceDet, IterRDM_new)
                else
                    full_sign = IterRDM_new*av_sign(1::nreplicas)*av_sign(nreplicas::nreplicas)/2.0_dp
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
                    call fill_diag_1rdm(one_rdms, nSpinCoup, real(SignFac,dp)*av_sign/sqrt(2.0_dp), &
                                       tCoreSpaceDet, IterRDM_new)
                else
                    full_sign = IterRDM_new*av_sign(1::nreplicas)*av_sign(nreplicas::nreplicas)/2.0_dp
                    call fill_spawn_rdm_diag(spawn, nI, full_sign)
                end if

                ! For HPHF we're considering < D_I + D_I' | a_a+ a_b+ a_j a_i | D_I + D_I' >
                ! Not only do we have diagonal < D_I | a_a+ a_b+ a_j a_i | D_I > terms, but also cross terms
                ! < D_I | a_a+ a_b+ a_j a_i | D_I' > if D_I and D_I' can be connected by a single or double 
                ! excitation. Find excitation level between D_I and D_I' and add in the contribution if connected.
                HPHFExcitLevel = FindBitExcitLevel(iLutnI, SpinCoupDet, 2)
                if (HPHFExcitLevel <= 2) then 
                    call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nSpinCoup, &
                                              IterRDM_new*av_sign(1::nreplicas)/sqrt(2.0_dp), &
                                              (real(SignFac,dp)*av_sign(nreplicas::nreplicas))/sqrt(2.0_dp), .true.)
                end if
            else

                ! HPHFs on, but determinant closed shell.
                if (RDMExcitLevel == 1) then
                    call fill_diag_1rdm(one_rdms, nI, av_sign, tCoreSpaceDet, IterRDM_new)
                else
                    full_sign = IterRDM_new*av_sign(1::nreplicas)*av_sign(nreplicas::nreplicas)
                    call fill_spawn_rdm_diag(spawn, nI, full_sign)
                end if

            end if
            call Add_RDM_HFConnections_HPHF(spawn, one_rdms, iLutnI, nI, av_sign, AvNoAtHF, ExcitLevelI, IterRDM_new)

        else
            ! Not using HPHFs.
            if (any(abs(av_sign(1::nreplicas) * av_sign(nreplicas::nreplicas)) > 1.0e-10_dp)) then
                if (RDMExcitLevel == 1) then
                    call fill_diag_1rdm(one_rdms, nI, av_sign, tCoreSpaceDet, IterRDM_new)
                else
                    full_sign = IterRDM_new*av_sign(1::nreplicas)*av_sign(nreplicas::nreplicas)
                    call fill_spawn_rdm_diag(spawn, nI, full_sign)
                end if
            end if

            call Add_RDM_HFConnections_Norm(spawn, one_rdms, iLutnI, nI, av_sign, AvNoAtHF, ExcitLevelI, IterRDM_new)
        end if

    end subroutine fill_rdm_diag_currdet_norm

    subroutine det_removed_fill_diag_rdm(spawn, one_rdms, iLutnI, av_sign, iter_occ)

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
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp), intent(in) :: av_sign(:), iter_occ(:)

        integer :: nI(nel), ExcitLevel

        ! If the determinant is removed on an iteration that the diagonal
        ! RDM elements are  already being calculated, it will already have
        ! been counted. So check this isn't the case first.
        if (.not. ((Iter == NMCyc) .or. (mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0))) then
            call decode_bit_det (nI, iLutnI)
            ExcitLevel = FindBitExcitLevel(iLutHF_True, iLutnI, 2)

            call fill_rdm_diag_currdet_norm(spawn, one_rdms, iLutnI, nI, ExcitLevel, av_sign, iter_occ, .false.)
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
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(nel)
        real(dp), intent(in) :: AvSignJ(:), AvSignHF(:)
        integer, intent(in) :: walkExcitLevel
        integer, intent(in) :: IterRDM(:)

        integer(n_int) :: iUnused

        ! If we have a single or double, add in the connection to the HF,
        ! symmetrically.
        if ((walkExcitLevel == 1) .or. (walkExcitLevel == 2)) then
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, HFDet_True, nJ, AvSignHF(1::nreplicas), &
                                      (1.0_dp/real(nreplicas,dp))*IterRDM*AvSignJ(nreplicas::nreplicas), .true.)

            call Add_RDM_From_IJ_Pair(spawn, one_rdms, HFDet_True, nJ, AvSignHF(nreplicas::nreplicas), &
                                      (1.0_dp/real(nreplicas,dp))*IterRDM*AvSignJ(1::nreplicas), .true.)
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
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(nel)
        real(dp), intent(in) :: AvSignJ(:), AvSignHF(:)
        integer, intent(in) :: walkExcitLevel
        integer, intent(in) :: IterRDM(:)

        ! Now if the determinant is connected to the HF (i.e. single or double),
        ! add in the diagonal elements of this connection as well -
        ! symmetrically because no probabilities are involved.
        if ((walkExcitLevel == 1) .or. (walkExcitLevel == 2)) then
            call Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutHF_True, iLutJ, HFDet_True, nJ, AvSignHF(1::nreplicas), &
                                          IterRDM*AvSignJ(nreplicas::nreplicas), .true.)
        end if

    end subroutine Add_RDM_HFConnections_HPHF

    subroutine check_fillRDM_DiDj(spawn, one_rdms, Spawned_No, iLutJ, realSignJ)

        ! The spawned parts contain the Dj's spawned by the Di's in CurrentDets.
        ! If the SpawnedPart is found in the CurrentDets list, it means that
        ! the Dj has a non-zero cj - and therefore the Di.Dj pair will have a
        ! non-zero ci.cj to contribute to the RDM. The index i tells us where
        ! to look in the parent array, for the Di's to go with this Dj.

        use DetBitOps, only: DetBitEq
        use FciMCData, only: iLutHF_True
        use rdm_data, only: one_rdm_t

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: Spawned_No
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), intent(in) :: realSignJ(lenof_sign)

        if (.not. DetBitEQ(iLutHF_True, iLutJ, NIfDBO)) then
                call DiDj_Found_FillRDM(spawn, one_rdms, Spawned_No, iLutJ, realSignJ)
        end if

    end subroutine check_fillRDM_DiDj
 
    subroutine DiDj_Found_FillRDM(spawn, one_rdms, Spawned_No, iLutJ, real_sign_j_all)

        ! This routine is called when we have found a Di (or multiple Di's)
        ! spawning onto a Dj with sign /= 0 (i.e. occupied). We then want to
        ! run through all the Di, Dj pairs and add their coefficients 
        ! (with appropriate de-biasing factors) into the 1 and 2 electron RDM.

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: DetBitEq
        use FciMCData, only: Spawned_Parents, Spawned_Parents_Index, iLutHF_True
        use rdm_data, only: one_rdm_t, nrdms
        use rdm_data, only: nrdms_each_simulation, rdm_replica_pairs, rdm_labels_for_sims
        use SystemData, only: nel, tHPHF

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: Spawned_No
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), intent(in) :: real_sign_j_all(lenof_sign)

        integer :: i, irdm, rdm_ind, nI(nel), nJ(nel)
        real(dp) :: realSignI
        real(dp) :: input_sign_i(nrdms), input_sign_j(nrdms)
        integer :: dest_part_type, source_part_type

        ! Spawning from multiple parents, to iLutJ, which has SignJ.        

        ! We are at position Spawned_No in the SpawnedParts array.
        ! Spawned_Parents_Index(1,Spawned_No) is therefore the start position
        ! of the list of parents (Di's) which spawned on the Dj in
        ! SpawnedParts(Spawned_No). There are Spawned_Parents_Index(2,Spawned_No)
        ! of these parent Di's. Spawned_Parents(0:NIfDBO,x) is the determinant Di,
        ! Spawned_Parents(NIfDBO+1,x) is the un-biased ci.

        ! Run through all Di's.

        do i = Spawned_Parents_Index(1,Spawned_No), &
                Spawned_Parents_Index(1,Spawned_No) + Spawned_Parents_Index(2,Spawned_No) - 1 

            if (DetBitEQ(iLutHF_True, Spawned_Parents(0:NIfDBO,i), NIfDBO)) then
                ! We've already added HF - S, and HF - D symmetrically.
                ! Any connection with the HF has therefore already been added.
                cycle
            end if
            
            call decode_bit_det (nI, Spawned_Parents(0:NIfDBO,i))
            call decode_bit_det (nJ, iLutJ)

            realSignI = transfer( Spawned_Parents(NIfDBO+1,i), realSignI )
            ! The original spawning event (and the RealSignI) came from this
            ! population.
            source_part_type = Spawned_Parents(NIfDBO+2,i)

            ! Loop over all RDMs to which the simulation with label
            ! source_part_type contributes to.
            do irdm = 1, nrdms_each_simulation(source_part_type)
                ! Get the label of the simulation that is paired with this, 
                ! replica, for this particular RDM.
                dest_part_type = rdm_replica_pairs(irdm,source_part_type)
                ! The label of the RDM that this be contributing to.
                rdm_ind = rdm_labels_for_sims(irdm, source_part_type)

                input_sign_i = 0.0_dp
                input_sign_j = 0.0_dp
                input_sign_i(rdm_ind) = realSignI
                input_sign_j(rdm_ind) = real_sign_j_all(dest_part_type)

                ! Given the Di,Dj and Ci,Cj - find the orbitals involved in the
                ! excitation, and therefore the RDM elements we want to add the
                ! Ci.Cj to. We have to halve the contributions for DR as we're
                ! summing in pairs that originated from spawning events in both
                ! pop 1 and pop 2 -- i.e., double counted wrt diagonal elements.
                if (tHPHF) then
                    call Fill_Spin_Coupled_RDM(spawn, one_rdms, Spawned_Parents(0:NIfDBO,i), iLutJ, &
                                               nI, nJ, (1.0_dp/real(nreplicas,dp))*input_sign_i, input_sign_j, .false.)
                else
                    call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, &
                                              (1.0_dp/real(nreplicas,dp))*input_sign_i, input_sign_j, .false.)
                end if
            end do

        end do

    end subroutine DiDj_Found_FillRDM

    subroutine Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutnI, iLutnJ, nI, nJ, realSignI, realSignJ, tFill_CiCj_Symm)

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
        logical, intent(in) :: tFill_CiCj_Symm

        integer(n_int) :: iLutnI2(0:NIfTot)
        integer :: nI2(nel), nJ2(nel)
        real(dp) :: NewSignJ(size(realSignJ)), NewSignI(size(realSignI))
        real(dp) :: PermSignJ(size(realSignJ)), PermSignI(size(realSignI))
        integer :: I_J_ExcLevel, ICoup_J_ExcLevel

        if (TestClosedShellDet(iLutnI)) then

            if (TestClosedShellDet(iLutnJ)) then
                ! Closed shell -> Closed shell - just as in determinant case
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, realSignI, realSignJ, tFill_CiCj_Symm)
            else
                ! Closed shell -> open shell.
                call FindDetSpinSym(nJ, nJ2, nel)
                NewSignJ = realSignJ/Root2
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, realSignI, NewSignJ, tFill_CiCj_Symm)
                ! What is the permutation between Di and Dj'
                NewSignJ = NewSignJ * hphf_sign(iLutnJ)
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ2, realSignI, NewSignJ, tFill_CiCj_Symm)
            end if

        else if (TestClosedShellDet(iLutnJ)) then
            ! Open shell -> closed shell
            call FindDetSpinSym(nI,nI2,nel)
            NewSignI = realSignI/Root2
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, NewSignI, realSignJ, tFill_CiCj_Symm)
            ! What is the permutation between Di' and Dj?
            NewSignI = NewSignI * hphf_sign(iLutnI)
            call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI2, nJ, NewSignI, realSignJ, tFill_CiCj_Symm)

        else
            ! Open shell -> open shell
            NewSignI = realSignI/Root2
            NewSignJ = realSignJ/Root2
            PermSignJ = NewSignJ * real(hphf_sign(iLutnJ),dp)
            PermSignI = NewSignI * real(hphf_sign(iLutnI),dp)
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            call FindDetSpinSym(nI, nI2, nel)
            call FindDetSpinSym(nJ, nJ2, nel)
            I_J_ExcLevel = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            ICoup_J_ExcLevel = FindBitExcitLevel(iLutnI2, iLutnJ, 2)

            if (I_J_ExcLevel .le. 2) then
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, NewSignI, NewSignJ, tFill_CiCj_Symm)
                ! Di -> Dj
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI2, nJ2, PermSignI, PermSignJ, tFill_CiCj_Symm)
                ! Di' -> Dj'  (both permuted sign)
            end if

            if (ICoup_J_ExcLevel .le. 2) then
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI2, nJ, PermSignI, NewSignJ, tFill_CiCj_Symm)
                ! Di' -> Dj  (i permuted sign)
                call Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ2, NewSignI, PermSignJ, tFill_CiCj_Symm)
                ! Di  -> Dj'  (j permuted sign)
            end if
        end if

    end subroutine Fill_Spin_Coupled_RDM

    subroutine Add_RDM_From_IJ_Pair(spawn, one_rdms, nI, nJ, realSignI, realSignJ, tFill_CiCj_Symm)

        ! This routine takes a pair of different determinants Di and Dj, and
        ! figures out which type of elements need to be added in to the RDM.

        use LoggingData, only: RDMExcitLevel
        use rdm_data, only: one_rdm_t
        use rdm_data_utils, only: add_to_rdm_spawn_t
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(in) :: realSignI(:), realSignJ(:)
        logical, intent(in) :: tFill_CiCj_Symm

        integer :: Ex(2,2), Ex_symm(2,2)
        logical :: tParity
        real(dp) :: full_sign(spawn%rdm_send%sign_length)

        Ex(:,:) = 0
        ! Maximum excitation level - we know they are connected by a double
        ! or a single excitation.
        Ex(1,1) = 2
        tParity = .false.

        ! Ex(1,:) comes out as the orbital(s) excited from, i.e. i,j.
        ! Ex(2,:) comes out as the orbital(s) excited to, i.e. a,b.
        call GetExcitation(nI, nJ, nel, Ex, tParity)

        full_sign = 0.0_dp
        if (tParity) then
            full_sign = -realSignI*realSignJ
        else
            full_sign = realSignI*realSignJ
        end if

        if ((Ex(1,2) .eq. 0) .and. (Ex(2,2) .eq. 0)) then
            
            ! Di and Dj are separated by a single excitation.
            ! Add in the contribution from this pair into the 1-RDM.
            
            if (RDMExcitLevel == 1) then
                call fill_sings_1rdm(one_rdms, Ex, tParity, realSignI, realSignJ, tFill_CiCj_Symm)
            else
                call fill_spawn_rdm_singles(spawn, nI, Ex, full_sign)
                if (tFill_CiCj_Symm) then
                    Ex_symm(1,:) = Ex(2,:)
                    Ex_symm(2,:) = Ex(1,:)
                    call fill_spawn_rdm_singles(spawn, nI, Ex_symm, full_sign)
                end if
            end if
    
        else if (RDMExcitLevel /= 1) then

            ! Otherwise Di and Dj are connected by a double excitation.
            ! Add in this contribution to the 2-RDM, if being calculated.
            call add_to_rdm_spawn_t(spawn, Ex(2,1), Ex(2,2), Ex(1,1), Ex(1,2), full_sign, .false.)
            if (tFill_CiCj_Symm) call add_to_rdm_spawn_t(spawn, Ex(1,1), Ex(1,2), Ex(2,1), Ex(2,2), full_sign, .false.)
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
        do iel = 1, nel-1
            do jel = iel+1, nel
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
        integer, intent(in) :: nI(nel), Ex(2,2)
        real(dp), intent(in) :: full_sign(spawn%rdm_send%sign_length)

        integer :: iel

        ! Looking at elements of the type Gamma(a,k,i,k).
        do iel = 1, nel
            associate(k => nI(iel))
                if (k < Ex(1,1) .and. k < Ex(2,1)) then
                    call add_to_rdm_spawn_t(spawn, k, Ex(2,1), k, Ex(1,1), full_sign, .false.)
                else if (k < Ex(1,1) .and. k > Ex(2,1)) then
                    call add_to_rdm_spawn_t(spawn, Ex(2,1), k, k, Ex(1,1), -full_sign, .false.)
                else if (k > Ex(1,1) .and. k < Ex(2,1)) then
                    call add_to_rdm_spawn_t(spawn, k, Ex(2,1), Ex(1,1), k, -full_sign, .false.)
                else if (k > Ex(1,1) .and. k > Ex(2,1)) then
                    call add_to_rdm_spawn_t(spawn, Ex(2,1), k, Ex(1,1), k, full_sign, .false.)
                end if
            end associate
        end do
        
    end subroutine fill_spawn_rdm_singles

! =======================================================================================    
! THESE NEXT ROUTINES ARE GENERAL TO BOTH STOCHASTIC AND EXPLICIT    
! =======================================================================================    

    subroutine fill_diag_1rdm(one_rdms, nI, contrib_sign, tCoreSpaceDetIn, RDMItersIn)

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
        use LoggingData, only: ThreshOccRDM, tThreshOccRDMDiag
        use RotateOrbsData, only: SymLabelListInv_rot
        use UMatCache, only: gtID

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(:)
        real(dp), intent(in) :: contrib_sign(:)
        logical, optional, intent(in) :: tCoreSpaceDetIn
        integer, optional, intent(in) :: RDMItersIn(:)

        integer :: i, ind, irdm
        integer :: RDMIters(size(one_rdms))
        real(dp) :: ScaleContribFac, final_contrib(size(one_rdms))
        logical :: tCoreSpaceDet

        ScaleContribFac = 1.0_dp

        RDMIters = 1
        if (present(RDMItersIn)) RDMIters = RDMItersIn

        tCoreSpaceDet = .false.
        if (present(tCoreSpaceDetIn)) tCoreSpaceDet = tCoreSpaceDetIn

        ! This is the single-run cutoff being applied (do not use in DR mode):
        if (.not. tCoreSpaceDet) then
            ! Dets in the core space are never removed from main list, so
            ! strictly do not require corrections.
            if (tThreshOccRDMDiag .and. (abs(contrib_sign(1)) .le. ThreshOccRDM)) ScaleContribFac = 0.0_dp
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

            final_contrib = contrib_sign(1::nreplicas) * contrib_sign(nreplicas::nreplicas) * RDMIters * ScaleContribFac

            do irdm = 1, size(one_rdms)
                one_rdms(irdm)%matrix(ind,ind) = one_rdms(irdm)%matrix(ind,ind) + final_contrib(irdm)
            end do
        end do

    end subroutine fill_diag_1rdm

    subroutine fill_sings_1rdm(one_rdms, Ex, tParity, contrib_sign_i, contrib_sign_j, fill_symmetric)

        use rdm_data, only: one_rdm_t, tOpenShell
        use RotateOrbsData, only: SymLabelListInv_rot
        use UMatCache, only: gtID

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: Ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: contrib_sign_i(:), contrib_sign_j(:)
        logical, intent(in) :: fill_symmetric

        integer :: i, a, ind_i, ind_a, irdm
        real(dp) :: ParityFactor

        ParityFactor = 1.0_dp
        if (tParity) ParityFactor = -1.0_dp

        ! Get the orbital labels involved in the excitation.
        if (tOpenShell) then
            i = Ex(1,1)
            a = Ex(2,1)
        else
            i = gtID(Ex(1,1))
            a = gtID(Ex(2,1))
        end if

        ! The SymLabelListInv_rot array is used to index the 1-RDM so that it
        ! will be filled in block diagonal order, where each block holds one
        ! symmetry sector of the 1-RDM (i.e., orbitals are ordered first by
        ! symmetry, then by the standard order).
        ind_i = SymLabelListInv_rot(i)
        ind_a = SymLabelListInv_rot(a)

        do irdm = 1, size(one_rdms)
            one_rdms(irdm)%matrix( ind_i, ind_a ) = one_rdms(irdm)%matrix( ind_i, ind_a ) + &
                                                    (ParityFactor * contrib_sign_i(irdm) * contrib_sign_j(irdm))
            if (fill_symmetric) then
                one_rdms(irdm)%matrix( ind_a, ind_i ) = one_rdms(irdm)%matrix( ind_a, ind_i ) + &
                                                        (ParityFactor * contrib_sign_i(irdm) * contrib_sign_j(irdm))
            end if
        end do

    end subroutine fill_sings_1rdm

    subroutine fill_RDM_offdiag_deterministic(spawn, one_rdms)

        use bit_rep_data, only: NIfD
        use bit_reps, only: decode_bit_det
        use DetBitOps, only: get_bit_excitmat, DetBitEq
        use FciMCData, only: Iter, IterRDMStart, PreviousCycles, iLutHF_True
        use FciMCData, only: core_space, determ_sizes, determ_displs, full_determ_vecs_av
        use LoggingData, only: RDMExcitLevel, RDMEnergyIter
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: one_rdm_t, nrdms, signs_for_rdm
        use rdm_data_utils, only: add_to_rdm_spawn_t
        use sparse_arrays, only: sparse_core_ham, core_connections
        use SystemData, only: nel, tHPHF

        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)

        integer :: i, j, irdm
        integer :: SingEx(2,1), Ex(2,2)
        real(dp) :: AvSignI(spawn%rdm_send%sign_length), AvSignJ(spawn%rdm_send%sign_length)
        real(dp) :: full_sign(spawn%rdm_send%sign_length)
        logical :: tParity
        integer(n_int) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer :: nI(nel), nJ(nel), IC
        integer :: IterRDM, connect_elem

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

        ! Cycle over all core dets on this process.
        do i = 1, determ_sizes(iProcIndex)
            iLutI = core_space(:,determ_displs(iProcIndex)+i)
                         
            ! Connections to the HF are added in elsewhere, so skip them here.
            if (DetBitEq(iLutI, iLutHF_True, NifDBO)) cycle
           
            do irdm = 1, nrdms
                AvSignI(irdm) = full_determ_vecs_av(signs_for_rdm(1,irdm), determ_displs(iProcIndex)+i)
            end do

            call decode_bit_det(nI,iLutI)
 
            do j = 1, sparse_core_ham(i)%num_elements-1
                ! Running over all non-zero off-diag matrix elements
                ! Connections to whole space (1 row), excluding diagonal elements

                ! Note: determ_displs holds sum(determ_sizes(0:proc-1))
                ! Core space holds all the core determinants on every processor,
                ! so we need to shuffle up to the range of indices corresponding
                ! to this proc (using determ_displs) and then select the
                ! correct one, i.
                
                iLutJ = core_space(:,core_connections(i)%positions(j))

                ! Connections to the HF are added in elsewhere, so skip them here.
                if (DetBitEq(iLutJ, iLutHF_True, NifDBO)) cycle
                
                do irdm = 1, nrdms
                    AvSignJ(irdm) = full_determ_vecs_av(signs_for_rdm(nreplicas, irdm), core_connections(i)%positions(j))
                end do

                connect_elem = core_connections(i)%elements(j)

                IC = abs(connect_elem)

                if (sign(1, connect_elem) > 0) then
                    tParity = .false.
                    ! For nreplicas=2, this will multiply the odd indices of AvSignI_Par
                    ! (representing replica 1) by the even indices of AvSignJ (representing
                    ! replica 2).
                    full_sign = AvSignI*AvSignJ * IterRDM
                else
                    tParity = .true.
                    full_sign = -AvSignI*AvSignJ * IterRDM
                end if

                if (tHPHF) then
                    call decode_bit_det(nJ, iLutJ)
                    call Fill_Spin_Coupled_RDM(spawn, one_rdms, iLutI, iLutJ, nI, nJ, AvSignI*IterRDM, AvSignJ, .false.)
                else
                    if (IC == 1) then
                        ! Single excitation - contributes to 1- and 2-RDM
                        ! (if calculated).
                         
                        ! Note: get_bit_excitmat may be buggy (DetBitOps),
                        ! but will do for now as we need the Ex...
                        call get_bit_excitmat(iLutI(0:NIfD),iLutJ(0:NIfD), SingEx, IC)
                        Ex(:,1) = SingEx(:,1)
                       
                        ! No need to explicitly fill symmetrically as we'll
                        ! generate pairs of determinants both ways around using
                        ! the connectivity matrix.
                        if (RDMExcitLevel /= 1) call fill_spawn_rdm_singles(spawn, nI, Ex, full_sign)

                        if (RDMExcitLevel == 1) then
                            call fill_sings_1rdm(one_rdms, Ex, tParity, AvSignI*IterRDM, AvSignJ, .false.)
                        end if

                    else if ((IC == 2) .and. (RDMExcitLevel /= 1)) then
                        ! Note: get_bit_excitmat may be buggy (DetBitOps),
                        ! but will do for now as we need the Ex...
                        call get_bit_excitmat(iLutI(0:NIfD), iLutJ(0:NIfD), Ex, IC)

                        call add_to_rdm_spawn_t(spawn, Ex(2,1), Ex(2,2), Ex(1,1), Ex(1,2), full_sign, .false.)
                    end if
                end if
            end do
        end do

    end subroutine fill_RDM_offdiag_deterministic 

end module rdm_filling
