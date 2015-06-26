#include "macros.h"

module rdm_filling

    ! This module contains routines used to perform filling of the RDM arrays,
    ! as done on-the-fly during an FCIQMC simulation.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants

    implicit none

contains

    subroutine fill_rdm_diag_currdet_norm(rdm, irdm, iLutnI, nI, j, ExcitLevelI, tCoreSpaceDet)

        ! This routine calculates the diagonal RDM contribution, and explicit
        ! connections to the HF, from the current determinant, for *every* RDM
        ! passed in through the rdms array.

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
        use FciMCData, only: Iter, IterRDMStart, PreviousCycles
        use global_det_data, only: get_iter_occ, get_av_sgn
        use hphf_integrals, only: hphf_sign
        use HPHFRandExcitMod, only: FindExcitBitDetSym
        use LoggingData, only: RDMEnergyIter
        use rdm_data, only: rdm_t
        use SystemData, only: nel, tHPHF

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: irdm
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: nI(nel), ExcitLevelI, j
        logical, intent(in), optional :: tCoreSpaceDet

        real(dp) :: IterDetOcc_all(lenof_sign), IterDetOcc_sing(nreplicas)
        real(dp) :: AvSignCurr_all(lenof_sign), AvSignCurr_sing(nreplicas)
        integer(n_int) :: SpinCoupDet(0:nIfTot)
        integer :: nSpinCoup(nel), SignFac, HPHFExcitLevel, part_type
        integer :: IterLastRDMFill, AvSignIters, IterRDM, sign_ind_1, sign_ind_2

        ! This is the number of iterations this determinant has been occupied,
        ! over *all* replicas.
        IterDetOcc_all(1:lenof_sign) = real(Iter+PreviousCycles,dp) - get_iter_occ(j) + 1.0_dp

        ! All signs from all RDMs.
        AvSignCurr_all = get_av_sgn(j)

        ! IterLastRDMFill is the number of iterations from the last time the
        ! energy was calculated.
        IterLastRDMFill = mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter)

            ! The indices of the signs for the RDM that we are considering.
            sign_ind_1 = nreplicas*irdm-nreplicas+1
            sign_ind_2 = nreplicas*irdm

            ! This is the number of iterations this determinant has been occupied,
            ! over the replicas relevant for the requested RDM.
            IterDetOcc_sing(1:nreplicas) = IterDetOcc_all(sign_ind_1:sign_ind_2)

            AvSignIters = min(IterDetOcc_sing(1), IterDetOcc_sing(nreplicas))

            ! The number of iterations we want to weight this RDM contribution by is:
            if (IterLastRDMFill .gt. 0) then
                IterRDM = min(AvSignIters, IterLastRDMFill)
            else
                IterRDM = AvSignIters
            end if

            ! The signs corresponding to this RDM.
            AvSignCurr_sing = AvSignCurr_all(sign_ind_1:sign_ind_2)

            if (tHPHF) then
                if (.not. TestClosedShellDet(iLutnI)) then
                    call Fill_Diag_RDM(rdm, nI, AvSignCurr_sing/sqrt(2.0_dp), tCoreSpaceDet, IterRDM)

                    ! C_X D_X = C_X / sqrt(2) [ D_I +/- D_I'] - for open shell dets,
                    ! divide stored C_X by sqrt(2). 
                    ! Add in I.
                    call FindExcitBitDetSym(iLutnI, SpinCoupDet)
                    call decode_bit_det(nSpinCoup, SpinCoupDet)
                    ! Find out if it's + or - in the above expression
                    SignFac = hphf_sign(iLutnI)

                    call Fill_Diag_RDM(rdm, nSpinCoup, real(SignFac,dp)*AvSignCurr_sing/sqrt(2.0_dp), &
                                       tCoreSpaceDet, IterRDM)

                    ! For HPHF we're considering < D_I + D_I' | a_a+ a_b+ a_j a_i | D_I + D_I' >
                    ! Not only do we have diagonal < D_I | a_a+ a_b+ a_j a_i | D_I > terms, but also cross terms
                    ! < D_I | a_a+ a_b+ a_j a_i | D_I' > if D_I and D_I' can be connected by a single or double 
                    ! excitation. Find excitation level between D_I and D_I' and add in the contribution if connected.
                    HPHFExcitLevel = FindBitExcitLevel(iLutnI, SpinCoupDet, 2)
                    if (HPHFExcitLevel .le. 2) then 
                        call Add_RDM_From_IJ_Pair(rdm, nI, nSpinCoup, IterRDM*AvSignCurr_sing(1)/sqrt(2.0_dp), &
                                                  (real(SignFac,dp)*AvSignCurr_sing(nreplicas))/sqrt(2.0_dp), .true.)
                    end if
                else

                    ! HPHFs on, but determinant closed shell.
                    call Fill_Diag_RDM(rdm, nI, AvSignCurr_sing, tCoreSpaceDet, IterRDM)

                end if
                call Add_RDM_HFConnections_HPHF(rdm, iLutnI, nI, AvSignCurr_sing, ExcitLevelI, IterRDM)

            else
                ! Not using HPHFs.
                if (AvSignCurr_sing(1)*AvSignCurr_sing(nreplicas) .ne. 0) &
                    call Fill_Diag_RDM(rdm, nI, AvSignCurr_sing, tCoreSpaceDet, IterRDM)

                call Add_RDM_HFConnections_Norm(rdm, iLutnI, nI, AvSignCurr_sing, ExcitLevelI, IterRDM)

            end if

    end subroutine fill_rdm_diag_currdet_norm

    subroutine det_removed_fill_diag_rdm(rdm, irdm, iLutnI, j)

        ! This routine is called if a determinant is removed from the list of
        ! currently occupied. At this point we need to add in its diagonal
        ! contribution for the number of iterations it has been occupied (or
        ! since the contribution was last included).

        ! This is done for *all* RDMs passed in through the rdms array.

        ! j --> which element of the main list CurrentDets are we considering.

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: FindBitExcitLevel
        use CalcData, only: NMCyc
        use FciMCData, only: iLutRef, iLutHF_True, Iter, IterRDMStart, PreviousCycles
        use LoggingData, only: RDMEnergyIter
        use rdm_data, only: rdm_t
        use SystemData, only: nel, tRef_Not_HF

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: irdm
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: j

        integer :: nI(nel), ExcitLevel, IterLastRDMFill

        ! If the determinant is removed on an iteration that the diagonal
        ! RDM elements are  already being calculated, it will already have
        ! been counted. So check this isn't the case first.
        if (.not. ((Iter .eq. NMCyc) .or. (mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) .eq. 0))) then
            call decode_bit_det (nI, iLutnI)
            if (tRef_Not_HF) then
                ExcitLevel = FindBitExcitLevel(iLutHF_True, iLutnI, 2)
            else
                ExcitLevel = FindBitExcitLevel(iLutRef, iLutnI, 2)
            end if

            call fill_rdm_diag_currdet_norm(rdm, irdm, iLutnI, nI, j, ExcitLevel, .false.)

        end if

    end subroutine det_removed_fill_diag_rdm

    subroutine Add_RDM_HFConnections_Norm(rdm, iLutJ, nJ, AvSignJ, walkExcitLevel, IterRDM)

        ! This is called when we run over all TotWalkers in CurrentDets.    
        ! It is called for each CurrentDet which is a single or double of the HF.
        ! It explicitly adds in the HF - S/D connection, as if the HF were D_i and 
        ! the single or double D_j. This is the standard full space RDM calc (No HPHF).
        ! In this case the diagonal elements wll already be taken care of.

        use FciMCData, only: InstNoatHF, HFDet_True, AvNoatHF
        use LoggingData, only: tFullHFAv
        use rdm_data, only: rdm_t
        use SystemData, only: nel

        type(rdm_t), intent(inout) :: rdm
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(nel)
        real(dp), intent(in) :: AvSignJ(nreplicas)
        integer, intent(in) :: IterRDM
        integer, intent(in) :: walkExcitLevel

        integer(n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(nel), HPHFExcitLevel, part_type

        !! Quick check that the HF population is being calculated correctly.
        !if (.not. tFullHFAv) then
        !    ! If tFullHFAv, we continue the accumulation of AvNoAtHF even when
        !    ! InstNoAtHF is zero. Therefore, AvNoAtHF is allowed to be different
        !    ! to the AvSignJ stored in CurrentH for this det.
        !    if (walkExcitLevel .eq. 0) then
        !        do part_type = 1, lenof_sign
        !            if (abs(AvSignJ(part_type)-AvNoatHF(part_type)) .gt. 1.0e-10_dp) then
        !                write(6,*) 'HFDet_True', HFDet_True
        !                write(6,*) 'nJ', nJ
        !                write(6,*) 'iLutJ', iLutJ
        !                write(6,*) 'AvSignJ', AvSignJ
        !                write(6,*) 'AvNoatHF', AvNoatHF
        !                write(6,*) "instnoathf", instnoathf
        !                call Stop_All('Add_RDM_HFConnections_Norm','Incorrect average HF population.')
        !            end if
        !        end do
        !    end if
        !end if

        ! If we have a single or double, add in the connection to the HF,
        ! symmetrically.
        if ((walkExcitLevel .eq. 1) .or. (walkExcitLevel .eq. 2)) then
            call Add_RDM_From_IJ_Pair(rdm, HFDet_True, nJ, AvNoatHF(1), &
                                      (1.0_dp/real(nreplicas,dp))*IterRDM*AvSignJ(nreplicas), .true.)

            call Add_RDM_From_IJ_Pair(rdm, HFDet_True, nJ, AvNoatHF(nreplicas), &
                                      (1.0_dp/real(nreplicas,dp))*IterRDM*AvSignJ(1), .true.)
        end if

    end subroutine Add_RDM_HFConnections_Norm

    subroutine Add_RDM_HFConnections_HPHF(rdm, iLutJ, nJ, AvSignJ, walkExcitLevel, IterRDM)

        ! This is called when we run over all TotWalkers in CurrentDets.
        ! It is called for each CurrentDet which is a single or double of the HF.
        ! It adds in the HF - S/D connection. The diagonal elements will already
        ! have been taken care of by the extract routine.

        use FciMCData, only: HFDet_True, iLutHF_True, AvNoatHF
        use LoggingData, only: tFullHFAv
        use rdm_data, only: rdm_t
        use SystemData, only: nel

        type(rdm_t), intent(inout) :: rdm
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(nel)
        integer, intent(in) :: IterRDM
        real(dp), intent(in) :: AvSignJ(nreplicas)
        integer, intent(in) :: walkExcitLevel

        integer(n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(nel), HPHFExcitLevel, part_type

        !if (.not. tFullHFAv) then
        !    ! If tFullHFAv, we continue the accumulation of AvNoAtHF even
        !    ! when InstNoAtHF is zero. Therefore, AvNoAtHF is allowed to be
        !    ! different to the AvSignJ stored in CurrentH for this det.
        !    if (walkExcitLevel .eq. 0) then
        !        do part_type = 1, lenof_sign
        !            if (AvSignJ(part_type) .ne. AvNoatHF(part_type)) then
        !                write(6,*) 'AvSignJ', AvSignJ
        !                write(6,*) 'AvNoatHF', AvNoatHF
        !                call Stop_All('Add_RDM_HFConnections_HPHF','Incorrect average HF population.')
        !            end if
        !        end do
        !    end if
        !end if

        ! Now if the determinant is connected to the HF (i.e. single or double),
        ! add in the diagonal elements of this connection as well -
        ! symmetrically because no probabilities are involved.
        if ((walkExcitLevel .eq. 1) .or. (walkExcitLevel .eq. 2)) &
            call Fill_Spin_Coupled_RDM(rdm, iLutHF_True, iLutJ, HFDet_True, nJ, AvNoatHF(1), &
                                          IterRDM*AvSignJ(nreplicas), .true.)

    end subroutine Add_RDM_HFConnections_HPHF

    subroutine check_fillRDM_DiDj(rdms, Spawned_No, iLutJ, realSignJ)

        ! The spawned parts contain the Dj's spawned by the Di's in CurrentDets.
        ! If the SpawnedPart is found in the CurrentDets list, it means that
        ! the Dj has a non-zero cj - and therefore the Di.Dj pair will have a
        ! non-zero ci.cj to contribute to the RDM. The index i tells us where
        ! to look in the parent array, for the Di's to go with this Dj.

        use DetBitOps, only: DetBitEq
        use FciMCData, only: iLutHF_True
        use rdm_data, only: rdm_t

        type(rdm_t), intent(inout) :: rdms(:)
        integer, intent(in) :: Spawned_No
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), intent(in) :: realSignJ(lenof_sign)

        integer :: ExcitLevel
    
        if (.not. DetBitEQ(iLutHF_True, iLutJ, NIfDBO)) then
                call DiDj_Found_FillRDM(rdms, Spawned_No, iLutJ, realSignJ)
        end if

    end subroutine check_fillRDM_DiDj
 
    subroutine DiDj_Found_FillRDM(rdms, Spawned_No, iLutJ, realSignJ)

        ! This routine is called when we have found a Di (or multiple Di's)
        ! spawning onto a Dj with sign /= 0 (i.e. occupied). We then want to
        ! run through all the Di, Dj pairs and add their coefficients 
        ! (with appropriate de-biasing factors) into the 1 and 2 electron RDM.

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: DetBitEq
        use FciMCData, only: Spawned_Parents, Spawned_Parents_Index, iLutHF_True
        use rdm_data, only: rdm_t
        use SystemData, only: nel, tHPHF

        type(rdm_t), intent(inout) :: rdms(:)
        integer, intent(in) :: Spawned_No
        integer(n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), intent(in) :: realSignJ(lenof_sign)

        integer :: i, j, rdm_ind, nI(nel), nJ(nel), walkExcitLevel
        real(dp) :: part_realSignI
        integer :: dest_part_type, source_part_type
        logical :: tParity, tDetAdded

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

            part_realSignI = transfer( Spawned_Parents(NIfDBO+1,i), part_realSignI )

            ! The original spawning event (and the RealSignI) came from this
            ! population.
            source_part_type = Spawned_Parents(NIfDBO+2,i)

            ! Get the index of the sign of the replica that is paired with
            ! this replica. Also get the label of the RDM to which this
            ! spawning event is contributing.
            if (nreplicas == 1) then
                ! This is the biased case - just use the same sign.
                dest_part_type = source_part_type
                rdm_ind = source_part_type
            else if (nreplicas == 2) then
                ! The unbiased case - get the sign from the replica simulation.
                dest_part_type = paired_replica(source_part_type)
                rdm_ind = (source_part_type+1)/2
            end if

            ! Given the Di,Dj and Ci,Cj - find the orbitals involved in the
            ! excitation, and therefore the RDM elements we want to add the
            ! Ci.Cj to. We have to halve the contributions for DR as we're
            ! summing in pairs that originated from spawning events in both
            ! pop 1 and pop 2 -- i.e., double counted wrt diagonal elements.
            if (tHPHF) then
                call Fill_Spin_Coupled_RDM(rdms(rdm_ind), Spawned_Parents(0:NIfDBO,i), iLutJ, nI, nJ, &
                           (1.0_dp/real(nreplicas,dp))*part_realSignI, realSignJ(dest_part_type), .false.)
            else
                call Add_RDM_From_IJ_Pair(rdms(rdm_ind), nI, nJ, (1.0_dp/real(nreplicas,dp))*part_realSignI, &
                                          realSignJ(dest_part_type), .false.)
            end if

        end do

    end subroutine DiDj_Found_FillRDM

    subroutine Fill_Spin_Coupled_RDM(rdm, iLutnI, iLutnJ, nI, nJ, realSignI, realSignJ, tFill_CiCj_Symm)

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
        use rdm_data, only: rdm_t
        use SystemData, only: nel, tOddS_hphf

        type(rdm_t), intent(inout) :: rdm
        integer(n_int), intent(in) :: iLutnI(0:NIfTot), iLutnJ(0:NIfTot)
        real(dp), intent(in) :: realSignI, realSignJ
        integer, intent(in) :: nI(nel),nJ(nel)
        logical, intent(in) :: tFill_CiCj_Symm

        integer(n_int) :: iLutnI2(0:NIfTot)
        integer :: nI2(nel), nJ2(nel)
        real(dp) :: NewSignJ, NewSignI, PermSignJ, PermSignI
        integer :: I_J_ExcLevel, ICoup_J_ExcLevel
        character(*), parameter :: t_r = 'Fill_Spin_Coupled_RDM'

        if (TestClosedShellDet(iLutnI)) then
            if (tOddS_HPHF) then
                call stop_all(t_r, "Should not be any closed shell determinants in high S states")
            end if

            if (TestClosedShellDet(iLutnJ)) then
                ! Closed shell -> Closed shell - just as in determinant case
                call Add_RDM_From_IJ_Pair(rdm, nI, nJ, realSignI, realSignJ, tFill_CiCj_Symm)
            else
                ! Closed shell -> open shell.
                call FindDetSpinSym(nJ, nJ2, nel)
                NewSignJ = realSignJ/Root2
                call Add_RDM_From_IJ_Pair(rdm, nI, nJ, realSignI, NewSignJ, tFill_CiCj_Symm)
                ! What is the permutation between Di and Dj'
                NewSignJ = NewSignJ * hphf_sign(iLutnJ)
                call Add_RDM_From_IJ_Pair(rdm, nI, nJ2, realSignI, NewSignJ, tFill_CiCj_Symm)
            end if

        else if (TestClosedShellDet(iLutnJ)) then
            ! Open shell -> closed shell
            call FindDetSpinSym(nI,nI2,nel)
            NewSignI = realSignI/Root2
            call Add_RDM_From_IJ_Pair(rdm, nI, nJ, NewSignI, realSignJ, tFill_CiCj_Symm)
            ! What is the permutation between Di' and Dj?
            NewSignI = NewSignI * hphf_sign(iLutnI)
            call Add_RDM_From_IJ_Pair(rdm, nI2, nJ, NewSignI, realSignJ, tFill_CiCj_Symm)

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
                call Add_RDM_From_IJ_Pair(rdm, nI, nJ, NewSignI, NewSignJ, tFill_CiCj_Symm)
                ! Di -> Dj
                call Add_RDM_From_IJ_Pair(rdm, nI2, nJ2, PermSignI, PermSignJ, tFill_CiCj_Symm)
                ! Di' -> Dj'  (both permuted sign)
            end if

            if (ICoup_J_ExcLevel .le. 2) then
                call Add_RDM_From_IJ_Pair(rdm, nI2, nJ, PermSignI, NewSignJ, tFill_CiCj_Symm)
                ! Di' -> Dj  (i permuted sign)
                call Add_RDM_From_IJ_Pair(rdm, nI, nJ2, NewSignI, PermSignJ, tFill_CiCj_Symm)
                ! Di  -> Dj'  (j permuted sign)
            end if
        end if

    end subroutine Fill_Spin_Coupled_RDM

    subroutine Add_RDM_From_IJ_Pair(rdm, nI, nJ, realSignI, realSignJ, tFill_CiCj_Symm)

        ! This routine takes a pair of different determinants Di and Dj, and
        ! figures out which type of elements need to be added in to the RDM.

        use LoggingData, only: RDMExcitLevel
        use rdm_data, only: rdm_t
        use SystemData, only: nel

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(in) :: realSignI, realSignJ
        logical, intent(in) :: tFill_CiCj_Symm

        integer :: Ex(2,2), j
        logical :: tParity

        Ex(:,:) = 0
        Ex(1,1) = 2         ! Maximum excitation level - we know they are connected by
                            ! a double or single.
        tParity = .false.

        ! Ex(1,:) comes out as the orbital(s) excited from, i.e. i,j.
        ! Ex(2,:) comes out as the orbital(s) excited to, i.e. a,b.
        call GetExcitation(nI, nJ, nel, Ex, tParity)

        if (Ex(1,1) .le. 0) then
            ! Error.
            write(6,*) '*'
            write(6,*) 'nI', nI
            write(6,*) 'nJ', nJ
            write(6,*) 'Ex(:,:)', Ex(1,1), Ex(1,2), Ex(2,1), Ex(2,2)
            write(6,*) 'tParity', tParity
            write(6,*) 'realSignI', realSignI
            write(6,*) 'realSignJ', realSignJ
            write(6,*) '*'
            call neci_flush(6)
            call stop_all('Add_RDM_From_IJ_Pair', 'Excitation level between pair not 1 or 2 as it should be.')
        end if

        if ((Ex(1,2) .eq. 0) .and. (Ex(2,2) .eq. 0)) then
            
            ! Di and Dj are separated by a single excitation.
            ! Add in the contribution from this pair into the 1-RDM.
            
            call Fill_Sings_RDM(rdm, nI, Ex, tParity, realSignI, realSignJ, tFill_CiCj_Symm)
    
        else if (RDMExcitLevel .ne. 1) then

            ! Otherwise Di and Dj are connected by a double excitation.
            ! Add in this contribution to the 2-RDM (as long as we're
            ! calculating this obv).
            call Fill_Doubs_RDM(rdm, Ex, tParity, realSignI, realSignJ, tFill_CiCj_Symm)

        end if

    end subroutine Add_RDM_From_IJ_Pair

! =======================================================================================    
! THESE NEXT ROUTINES ARE GENERAL TO BOTH STOCHASTIC AND EXPLICIT    
! =======================================================================================    

    subroutine Fill_Diag_RDM(rdm, nI, realSignDi, tCoreSpaceDetIn, RDMItersIn)

        ! Fill diagonal elements of 1- and 2-RDM.
        ! These are < Di | a_i+ a_i | Di > and < Di | a_i+ a_j+ a_j a_i | Di >.

        use rdm_data, only: rdm_t, tOpenShell
        use LoggingData, only: RDMExcitLevel, ThreshOccRDM, tThreshOccRDMDiag
        use NatOrbsMod, only: NatOrbMat
        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: nel
        use UMatCache, only: gtID

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: realSignDi(nreplicas)
        logical, intent(in), optional :: tCoreSpaceDetIn
        integer, intent(in), optional :: RDMItersIn

        integer :: i, j, iSpat, jSpat, Ind, iInd
        real(dp) :: ScaleContribFac
        integer :: RDMIters
        logical :: tCoreSpaceDet

        ScaleContribFac = 1.0
        
        RDMIters = 1.0_dp
        if (present(RDMItersIn)) RDMIters = RDMItersIn

        tCoreSpaceDet = .false.
        if (present(tCoreSpaceDetIn)) tCoreSpaceDet = tCoreSpaceDetIn

        ! This is the single-run cutoff being applied (do not use in DR mode):
        if (.not. tCoreSpaceDetIn) then
            ! Dets in the core space are never removed from main list, so
            ! strictly do not require corrections.
            if (tThreshOccRDMDiag .and. (abs(RealSignDi(1)) .le. ThreshOccRDM)) ScaleContribFac = 0.0_dp
        end if
        
        if (RDMExcitLevel .eq. 1) then
            do i = 1, nel
                if (tOpenShell) then
                    iInd = SymLabelListInv_rot(nI(i))
                else 
                    ! SymLabelListInv_rot will be in spat orbitals too.
                    iInd = SymLabelListInv_rot(gtID(nI(i)))
                end if
                NatOrbMat(iInd,iInd) = NatOrbMat(iInd,iInd) &
                                          + ( realSignDi(1) * realSignDi(nreplicas) * RDMIters )*ScaleContribFac 
            end do
        else
            ! Only calculating 2-RDM.
            ! nI(i) - spin orbital label. Odd=beta, even=alpha.
            do i = 1, nel - 1
                iSpat = gtID(nI(i))
                if (tOpenShell) iSpat = (nI(i)-1)/2 + 1

                ! Orbitals in nI ordered lowest to highest so nI(j) > nI(i),
                ! and jSpat >= iSpat (can only be equal if different spin).
                do j = i+1, nel
                    jSpat = gtID(nI(j))
                     if (tOpenShell) jSpat = (nI(j)-1)/2 + 1 
               
                    ! Either alpha alpha or beta beta -> aaaa/bbbb arrays.
                    if ( ((mod(nI(i),2) .eq. 1) .and. (mod(nI(j),2) .eq. 1)) .or. &
                        ((mod(nI(i),2) .eq. 0) .and. (mod(nI(j),2) .eq. 0)) ) then

                        ! Ind doesn't include diagonal terms (when iSpat == jSpat).
                        Ind = ( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                        if (( mod(nI(i),2) .eq. 0) .or. (.not. tOpenShell))then
                            ! nI(i) is even --> aaaa.
                            rdm%aaaa( Ind, Ind ) = rdm%aaaa( Ind, Ind ) &
                                          + ( realSignDi(1) * realSignDi(nreplicas) * RDMIters)*ScaleContribFac

                        else if ( mod(nI(i),2) .eq. 1)then
                            ! nI(i) is odd --> bbbb.
                            rdm%bbbb( Ind, Ind ) = rdm%bbbb( Ind, Ind ) &
                                          + ( realSignDi(1) * realSignDi(nreplicas) * RDMIters)*ScaleContribFac

                        end if

                    ! Either alpha beta or beta alpha -> abab/baba arrays.
                    else

                        ! Ind does include diagonal terms (when iSpat == jSpat).
                        Ind = ( ( (jSpat-1) * jSpat ) / 2 ) + iSpat

                        if (jSpat .eq. iSpat) then
                                ! aSpat == bSpat == iSpat == jSpat terms are
                                ! saved in abab only.
                                rdm%abab( Ind, Ind ) = rdm%abab( Ind, Ind ) &
                                                + ( realSignDi(1) * realSignDi(nreplicas) * RDMIters)*ScaleContribFac

                        else

                            if ((mod(nI(i),2) .eq. 0) .or. (.not. tOpenShell)) then
                                ! nI(i) is even ---> abab.
                                rdm%abab( Ind, Ind ) = rdm%abab( Ind, Ind ) &
                                          + ( realSignDi(1) * realSignDi(nreplicas) * RDMIters)*ScaleContribFac
                            else if (mod(nI(i),2) .eq. 1) then
                                ! nI(i) is odd ---> baba.
                                rdm%baba( Ind, Ind ) = rdm%baba( Ind, Ind ) &
                                               + ( realSignDi(1) * realSignDi(nreplicas) * RDMIters)*ScaleContribFac
                            end if

                       end if 

                    end if

                end do
            end do
        end if

    end subroutine Fill_Diag_RDM

    subroutine Fill_Sings_RDM(rdm, nI, Ex, tParity, realSignDi, realSignDj, tFill_CiCj_Symm)

        ! This routine adds in the contribution to the 1- and 2-RDM from
        ! determinants connected by a single excitation.

        use rdm_data, only: rdm_t, tOpenShell
        use LoggingData, only: RDMExcitLevel
        use NatOrbsMod, only: NatOrbMat
        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: nel
        use UMatCache, only: gtID

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: nI(nel), Ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: realSignDi, realSignDj
        logical, intent(in) :: tFill_CiCj_Symm

        integer :: k, Indik, Indak, iSpat, aSpat, kSpat, iInd, aInd
        real(dp) :: ParityFactor, ParityFactor2

        ParityFactor = 1.0_dp
        if (tParity) ParityFactor = -1.0_dp

        if (RDMExcitLevel .eq. 1) then

            ! SymLabelList2_rot(i) gives the orbital in position i
            ! SymLabelListInv_rot(i) gives the position orbital i should go in.
            if (tOpenShell) then
                iInd = Ex(1,1)
                aInd = Ex(2,1)
            else
                iInd = gtID(Ex(1,1))
                aInd = gtID(Ex(2,1))   ! These two must have the same spin.
            end if
            Indik = SymLabelListInv_rot(iInd)    ! Position of i 
            Indak = SymLabelListInv_rot(aInd)    ! Position of a.
            
            ! Adding to 1-RDM(i,a), ci.cj effectively.
            NatOrbMat( Indik, Indak ) = NatOrbMat( Indik, Indak ) + (ParityFactor * realSignDi * realSignDj)

            if (tFill_CiCj_Symm) then                                
                NatOrbMat( Indak, Indik ) = NatOrbMat( Indak, Indik ) + (ParityFactor * realSignDi * realSignDj)
            end if
        else
            ! Looking at elements of the type Gamma(i,k,a,k).

            ! The two determinants Di and Dj will have the same occupations
            ! except for the i and a. Any of the N-1 other electrons can be
            ! annihilated and created in the same orbital. So we run over
            ! all k = all N-1 other occupied orbitals.
            
            iSpat = gtID(Ex(1,1))
            aSpat = gtID(Ex(2,1))  ! These two must have the same spin.
            if (tOpenShell) then
                iSpat = (Ex(1,1)-1)/2 + 1
                aSpat = (Ex(2,1)-1)/2 + 1
            end if

            do k = 1, nel

                kSpat = gtID(nI(k))
                if (tOpenShell) kSpat = (nI(k)-1)/2 + 1

                if (nI(k) .ne. Ex(1,1)) then

                    if ((iSpat .eq. kSpat) .or. (aSpat .eq. kSpat)) then
                        ! It is possible for i = k or a = k if they 
                        ! have different spins. the only arrays with
                        ! i = j or a = b are abab/baba.
                        ! -> abba/baab terms must be reordered to
                        ! become abab/baba

                        Indik = ( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                        Indak = ( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                        if ((iSpat .eq. aSpat) .or. (.not. tOpenShell) ) then

                            ! The iSpat == jSpat == aSpat == bSpat term is
                            ! saved in the abab array only (not in baba).
                            rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                            if (tFill_CiCj_Symm) then
                                rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                            end if

                        else if (iSpat .eq. kSpat) then

                            ! We get the term k i -> k a, which is abab or baba.
                            ! If iSpat == kSpat and aSpat > kSpat, the indeces are
                            ! already ordered correctly. If they are not ordered
                            ! correctly (aSpat<kSpat), then we need to swap a and k
                            ! This results into an abba/baab term, but for equal
                            ! spatial orbitals, there is no abba/baab array, and
                            ! so we save the equivalent abab/baba term. Therefore
                            ! i and k have to be swapped as well, giving an abab/baba
                            ! term. Because there are none or two swaps, the parity does
                            ! not change!

                            if (aSpat .gt. kSpat) then

                                if ( (mod(Ex(2,1),2) .eq. 1) .or. (.not. tOpenShell) )then ! ki, ka -> last index beta
                                    rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if

                                else if (mod(Ex(2,1),2) .eq. 0)then ! ki, ka -> last index alpha
                                    rdm%baba( Indik, Indak ) = rdm%baba( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%baba( Indak, Indik ) = rdm%baba( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            else if  (aSpat .lt. kSpat) then

                                if ( (mod(Ex(2,1),2) .eq. 0) .or. (.not. tOpenShell) )then ! ik, ak -> third index alpha
                                    rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(2,1),2) .eq. 1)then ! ik, ak -> third index beta
                                    rdm%baba( Indik, Indak ) = rdm%baba( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%baba( Indak, Indik ) = rdm%baba( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            end if

                        else if (aSpat .eq. kSpat) then

                            if (iSpat .gt. kSpat) then 
                                if ( (mod(Ex(1,1),2) .eq. 1) .or. (.not. tOpenShell) )then ! ki, ka -> second index beta
                                    rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(1,1),2) .eq. 0)then ! ki, ka -> second index alpha
                                    rdm%baba( Indik, Indak ) = rdm%baba( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%baba( Indak, Indik ) = rdm%baba( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            else if  (iSpat .lt. kSpat) then

                                if ( (mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell) )then ! ik, ak -> first index alpha
                                    rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(1,1),2) .eq. 1)then ! ik, ak -> first index beta
                                    rdm%baba( Indik, Indak ) = rdm%baba( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%baba( Indak, Indik ) = rdm%baba( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            end if
                        end if  ! a=k

                    else ! not (iSpat.eq.kSpat) .or. (aSpat.eq.kSpat))
                        ! Checking spins of i and k.
                        ! If same, i.e alpha alpha or beta beta -> aaaa array.
                        if ( ((mod(Ex(1,1),2) .eq. 1) .and. (mod(nI(k),2) .eq. 1)) .or. &
                             ((mod(Ex(1,1),2) .eq. 0) .and. (mod(nI(k),2) .eq. 0)) ) then

                            ! 2-RDM(i,j,a,b) is defined to have i < j and a < b, as that is how the unique 
                            ! indices are defined for i,j and a,b.
                            ! But the parity is defined so that the i -> a excitation is aligned.

                            ! I.e. we're adding these as nI(k),Ex(1,1) -> nI(k), Ex(2,1)
                            ! So if Ex(1,1) < nI(k), or Ex(2,1) < nI(k) then we need 
                            ! to switch the parity.
                            ParityFactor2 = ParityFactor
                            if ((Ex(1,1) .lt. nI(k)) .and. (Ex(2,1) .gt. nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.0_dp)
                            if ((Ex(1,1) .gt. nI(k)) .and.( Ex(2,1) .lt. nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.0_dp)

                            ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                            Indik = ( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                            Indak = ( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                            ! nI(k) even or odd. odd=bbbb, even=aaaa.                          
                            if ((mod(nI(k),2) .eq. 0) .or. (.not. tOpenShell)) then
                                rdm%aaaa( Indik, Indak ) = rdm%aaaa( Indik, Indak ) + ( ParityFactor2 * &
                                                                                    realSignDi * realSignDj )
                                if (tFill_CiCj_Symm) then
                                    rdm%aaaa( Indak, Indik ) = rdm%aaaa( Indak, Indik ) + ( ParityFactor2 * &
                                                                                    realSignDi * realSignDj )
                                end if

                            else if (mod(nI(k),2) .eq. 1) then
                                rdm%bbbb( Indik, Indak ) = rdm%bbbb( Indik, Indak ) + ( ParityFactor2 * &
                                                                                 realSignDi * realSignDj )
                                if (tFill_CiCj_Symm) then
                                        rdm%bbbb( Indak, Indik ) = rdm%bbbb( Indak, Indik ) + ( ParityFactor2 * &
                                                                                 realSignDi * realSignDj )
                                end if
                            end if

                        ! either abab/baba or abba/baab array. 
                        ! we distinguish between these because i<j and a<b.
                        else   ! abab/baba or abba/baab

                            if ( (Ex(1,1) .lt. nI(k)) .and. (Ex(2,1) .lt. nI(k)) ) then
                            !i k a k -> abab/baba 

                                ! It is possible for i = k or j = k if they are spat orbitals 
                                ! and have different spins.
                                Indik = ( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                                Indak = ( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                                ! Ex(1,1) (i spin orb): first index even or odd 
                                if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell)) then
                                    rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(1,1),2).eq.1)then
                                    rdm%baba( Indik, Indak ) = rdm%baba( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                            rdm%baba( Indak, Indik ) = rdm%baba( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj ) 
                                    end if
                                end if

                            else if ( (Ex(1,1) .gt. nI(k)) .and. (Ex(2,1) .gt. nI(k)) ) then
                            ! k i k a -> abab/baba 

                                Indik = ( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                                Indak = ( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                                !Ex(1,1) (i spin orb): second index even or odd 
                                if ((mod(Ex(1,1),2).eq.1) .or. (.not. tOpenShell)) then
                                    rdm%abab( Indik, Indak ) = rdm%abab( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                    if (tFill_CiCj_Symm) then
                                        rdm%abab( Indak, Indik ) = rdm%abab( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj ) 
                                    end if
                                else if (mod(Ex(1,1),2).eq.0) then
                                    rdm%baba( Indik, Indak ) = rdm%baba( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                    if (tFill_CiCj_Symm) then
                                            rdm%baba( Indak, Indik ) = rdm%baba( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj )  
                                    end if
                                end if


                            else if ( (Ex(1,1) .gt. nI(k)) .and. (Ex(2,1) .lt. nI(k)) ) then 
                            ! k i a k -> abba/baab

                                ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                                Indik = ( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                                Indak = ( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                                !i spin orb: second odd or even. 
                                if ((mod(Ex(1,1),2) .eq. 1) .or. (.not. tOpenShell))then
                                    rdm%abba( Indik, Indak ) = rdm%abba( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                    if (tFill_CiCj_Symm ) then
                                        if (.not. tOpenShell) then
                                            rdm%abba( Indak, Indik ) = rdm%abba( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        else
                                            rdm%baab( Indak, Indik ) = rdm%baab( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        end if
                                    end if

                                else if (mod(Ex(1,1),2) .eq. 0)then
                                    rdm%baab( Indik, Indak ) = rdm%baab( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abba( Indak, Indik ) = rdm%abba( Indak, Indik ) - ( ParityFactor * &
                                                                                     realSignDi * realSignDj )
                                    end if
                                end if

                            else if ( (Ex(1,1) .lt. nI(k)) .and. (Ex(2,1) .gt. nI(k)) ) then 
                            ! i k k a -> abba/baab

                                ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                                Indik = ( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                                Indak = ( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                                !i spin orb: first index odd or even.
                                if ((mod(Ex(1,1),2) .eq. 0) .or. (.not. tOpenShell)) then
                                    rdm%abba( Indik, Indak ) = rdm%abba( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        if (.not. tOpenShell) then
                                            rdm%abba( Indak, Indik ) = rdm%abba( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        else
                                            rdm%baab( Indak, Indik ) = rdm%baab( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        end if
                                    end if

                                else if (mod(Ex(1,1),2) .eq. 1) then
                                    rdm%baab( Indik, Indak ) = rdm%baab( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        rdm%abba( Indak, Indik ) = rdm%abba( Indak, Indik ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    end if
                                end if

                            end if ! order of k i k j

                        end if !abab/baba or abba/baab
                    end if  ! not (iSpat.eq.kSpat).or.(aSpat.eq.kSpat))

                end if 
            end do 

        end if  

    end subroutine Fill_Sings_RDM

    subroutine Fill_Doubs_RDM(rdm, Ex, tParity, realSignDi, realSignDj, tFill_CiCj_Symm)

        ! This routine adds in the contribution to the 2-RDM from determinants
        ! connected by a double excitation.

        use rdm_data, only: rdm_t, tOpenShell
        use UMatCache, only: gtID

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: Ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: realSignDi, realSignDj
        logical, intent(in) :: tFill_CiCj_Symm

        integer :: Indij, Indab, iSpat, jSpat, aSpat, bSpat
        real(dp) :: ParityFactor

        ! Adding to elements Gamma(i,j,a,b)

        ParityFactor = 1.0_dp
        if (tParity) ParityFactor = -1.0_dp

        iSpat = gtID(Ex(1,1))
        jSpat = gtID(Ex(1,2))       ! Ex(1,1) < Ex(1,2)
        aSpat = gtID(Ex(2,1)) 
        bSpat = gtID(Ex(2,2))       ! Ex(2,1) < Ex(2,2)

        if (tOpenShell)then 
            iSpat = (Ex(1,1)-1)/2 + 1
            jSpat = (Ex(1,2)-1)/2 + 1
            aSpat = (Ex(2,1)-1)/2 + 1
            bSpat = (Ex(2,2)-1)/2 + 1
        end if         

        if ((iSpat .eq. jSpat) .or. (aSpat .eq. bSpat)) then

            ! if i and a are different spin -> abba (but adding as abab - mult by -1).
            if ( ((mod(Ex(1,1),2) .eq. 0) .and. (mod(Ex(2,1),2) .eq. 1)) .or. &
                ((mod(Ex(1,1),2) .eq. 1) .and. (mod(Ex(2,1),2) .eq. 0)) ) &
                    ParityFactor = ParityFactor * (-1.0_dp)

            Indij=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
            Indab=( ( (bSpat-1) * bSpat ) / 2 ) + aSpat

            if ((iSpat.eq.jSpat).and.(aSpat.eq.bSpat)) then
                ! abab and baba terms are equal and are saved in abab only.

                rdm%abab( Indij, Indab ) = rdm%abab( Indij, Indab ) + ( ParityFactor * &
                                                         realSignDi * realSignDj )
                if (tFill_CiCj_Symm) then
                    rdm%abab( Indab, Indij ) = rdm%abab( Indab, Indij ) + ( ParityFactor * &
                                                            realSignDi * realSignDj )
                end if

            else if (iSpat .eq. jSpat) then
               ! i and j may have to be swapped to get abab/baba term -> get
               ! spin from a (third index).
               if ((mod(Ex(2,1),2) .eq. 0) .or. (.not. tOpenShell)) then
                    rdm%abab( Indij, Indab ) = rdm%abab( Indij, Indab ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        rdm%abab( Indab, Indij ) = rdm%abab( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               else if (mod(Ex(2,1),2) .eq. 1) then
                    rdm%baba( Indij, Indab ) = rdm%baba( Indij, Indab ) + ( ParityFactor * &
                                                                  realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        rdm%baba( Indab, Indij ) = rdm%baba( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               end if

            else if (aSpat .eq. bSpat) then

               ! a and b may have to be swapped to get abab/baba term -> get
               ! spin from i (first index)
               if ((mod(Ex(1,1),2) .eq. 0) .or. (.not. tOpenShell)) then
                    rdm%abab( Indij, Indab ) = rdm%abab( Indij, Indab ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        rdm%abab( Indab, Indij ) = rdm%abab( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               else if (mod(Ex(1,1),2) .eq. 1) then
                    rdm%baba( Indij, Indab ) = rdm%baba( Indij, Indab ) + ( ParityFactor * &
                                                                  realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        rdm%baba( Indab, Indij ) = rdm%baba( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               end if

            end if

          else ! not ((iSpat.eq.jSpat) .or. (aSpat.eq.bSpat))

            ! Checking spins of i and j (these must be same combination as a and b).
            ! If alpha alpha or beta beta -> aaaa array.
            if ( ((mod(Ex(1,1),2) .eq. 1) .and. (mod(Ex(1,2),2) .eq. 1)) .or. &
                 ((mod(Ex(1,1),2) .eq. 0) .and. (mod(Ex(1,2),2) .eq. 0)) ) then

                ! Don't need to worry about diagonal terms, i can't equal j.
                ! jSpat > iSpat and bSpat > aSpat
                Indij = ( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                Indab = ( ( (bSpat-2) * (bSpat-1) ) / 2 ) + aSpat

                if ((mod(Ex(1,1),2) .eq. 0) .or. (.not. tOpenShell)) then
                    rdm%aaaa( Indij, Indab ) = rdm%aaaa( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                    if (tFill_CiCj_Symm) then
                            rdm%aaaa( Indab, Indij ) = rdm%aaaa( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                    end if

                else if (mod(Ex(1,1),2) .eq. 1)then
                    rdm%bbbb( Indij, Indab ) = rdm%bbbb( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                    if (tFill_CiCj_Symm) then
                          rdm%bbbb( Indab, Indij ) = rdm%bbbb( Indab, Indij ) + (ParityFactor * &
                                                                    realSignDi * realSignDj )
                    end if
                end if
                
            ! Either alpha beta or beta alpha -> abab array.
            else
                
                ! If when ordering i < j and a < b, is it abab or abba.

                ! i and a are the same spin -> abab/baba.
                if ( ((mod(Ex(1,1),2) .eq. 0) .and. (mod(Ex(2,1),2) .eq. 0)) .or. &
                     ((mod(Ex(1,1),2) .eq. 1) .and. (mod(Ex(2,1),2) .eq. 1)) ) then

                    Indij = ( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
                    Indab = ( ( (bSpat-1) * bSpat ) / 2 ) + aSpat

                    if ((mod(Ex(1,1), 2) .eq. 0) .or. (.not. tOpenShell)) then
                        rdm%abab( Indij, Indab ) = rdm%abab( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                            rdm%abab( Indab, Indij ) = rdm%abab( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                        end if
                    else if (mod(Ex(1,1), 2) .eq. 1) then
                        rdm%baba( Indij, Indab ) = rdm%baba( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                                rdm%baba( Indab, Indij ) = rdm%baba( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                        end if
                    end if

                ! i and a are different spin -> abba
                ! the only double excitation case with Indij = Indab will go
                ! in here.
                else
                    ! Don't need to worry about diagonal terms, i can't equal j.
                    ! jSpat > iSpat and bSpat > aSpat
                    Indij = ( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                    Indab = ( ( (bSpat-2) * (bSpat-1) ) / 2 ) + aSpat

                    if ((mod(Ex(1,1),2) .eq. 0) .or. (.not. tOpenShell))then
                        rdm%abba( Indij, Indab ) = rdm%abba( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                            if (.not. tOpenShell)then
                                rdm%abba( Indab, Indij ) = rdm%abba( Indab, Indij ) + ( ParityFactor * &
                                                                        realSignDi * realSignDj )
                            else
                                rdm%baab( Indab, Indij ) = rdm%baab( Indab, Indij ) + ( ParityFactor * &
                                                                        realSignDi * realSignDj )
                            end if
                        end if

                    else if (mod(Ex(1,1),2) .eq. 1) then

                        rdm%baab( Indij, Indab ) = rdm%baab( Indij, Indab ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                                rdm%abba( Indab, Indij ) = rdm%abba( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                        end if

                    end if

                end if
            end if
        end if

    end subroutine Fill_Doubs_RDM

    subroutine fill_RDM_offdiag_deterministic()

        use bit_rep_data, only: NIfD
        use bit_reps, only: decode_bit_det
        use DetBitOps, only: get_bit_excitmat, DetBitEq
        use FciMCData, only: Iter, IterRDMStart, PreviousCycles, iLutHF_True
        use FciMCData, only: core_space, determ_sizes, determ_displs, full_determ_vecs_av
        use LoggingData, only: RDMExcitLevel, RDMEnergyIter
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdms
        use sparse_arrays, only: sparse_core_ham, core_connections
        use SystemData, only: nel, tHPHF

        integer :: i, j
        integer :: SingEx(2,1), Ex(2,2)
        real(dp) :: InstSignI, InstSignJ
        real(dp) :: AvSignI, AvSignJ
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
        
        Ex(:,:) = 0

        ! Cycle over all core dets on this process.
        do i = 1, determ_sizes(iProcIndex)
            iLutI = core_space(:,determ_displs(iProcIndex)+i)
                         
            ! Connections to the HF are added in elsewhere, so skip them here.
            if (DetBitEq(iLutI, iLutHF_True, NifDBO)) cycle
           
            AvSignI = full_determ_vecs_av(1,determ_displs(iProcIndex)+i)

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
                 
                 AvSignJ = full_determ_vecs_av(inum_runs,core_connections(i)%positions(j))

                 connect_elem = core_connections(i)%elements(j)

                 IC = abs(connect_elem)

                 if (sign(1, connect_elem) .gt. 0) then
                     tParity = .false.
                 else
                     tParity = .true.
                 end if

                 if (tHPHF) then
                     call decode_bit_det(nJ, iLutJ)

                     call Fill_Spin_Coupled_RDM(rdms(1), iLutI, iLutJ, nI, nJ, AvSignI*IterRDM, AvSignJ, .false.)
                 else
                     if (IC .eq. 1) then
                         ! Single excitation - contributes to 1- and 2-RDM
                         ! (if calculated).
                          
                         ! Note: get_bit_excitmat may be buggy (DetBitOps),
                         ! but will do for now as we need the Ex...
                         call get_bit_excitmat(iLutI(0:NIfD),iLutJ(0:NIfD), SingEx, IC)
                         Ex(:,1) = SingEx(:,1)
                        
                         ! No need to explicitly fill symmetrically as we'll
                         ! generate pairs of determinants both ways around using
                         ! the connectivity matrix.
                         call Fill_Sings_RDM(rdms(1), nI, Ex, tParity, AvSignI*IterRDM, AvSignJ, .false.)

                     else if ((IC .eq. 2) .and. (RDMExcitLevel .ne. 1)) then
                         
                         ! Note: get_bit_excitmat may be buggy (DetBitOps),
                         ! but will do for now as we need the Ex...
                         call get_bit_excitmat(iLutI(0:NIfD), iLutJ(0:NIfD), Ex, IC)
                         call Fill_Doubs_RDM(rdms(1), Ex, tParity, AvSignI*IterRDM, AvSignJ, .false.)
                     end if
                 end if
             end do
        end do

    end subroutine fill_RDM_offdiag_deterministic 

end module rdm_filling
