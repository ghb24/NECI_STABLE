#include "macros.h"
#ifndef __CMPLX

module guga_rdm
    ! RDM module specifically for the GUGA spin-adapted implementation 
    
    use constants, only: n_int, dp, lenof_sign, EPS, sizeof_int, int_rdm
    use SystemData, only: nel, nSpatOrbs, current_stepvector, currentB_ilut
    use bit_reps, only: extract_bit_rep, decode_bit_det, niftot, nifdbo
    use LoggingData, only: RDMExcitLevel
    use rdm_data, only: one_rdms, two_rdm_spawn, rdmCorrectionFactor
    use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs, rdm_spawn_t, one_rdm_t
    use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, rdm_list_t
    use rdm_data, only: Sing_ExcList, Doub_ExcList, OneEl_Gap, TwoEl_Gap
    use DetBitOps, only: EncodeBitDet, count_open_orbs, DetBitEq
    use load_balance_calcnodes, only: DetermineDetNode
    use guga_excitations, only: init_csf_information, excitationIdentifier
    use guga_excitations, only: init_singleWeight, calcRemainingSwitches
    use guga_excitations, only: createSingleStart, singleUpdate, singleEnd
    use guga_excitations, only: checkCompatibility
    use guga_excitations, only: calcDoubleExcitation_withWeight, &
                                calcNonOverlapDouble, calcSingleOverlapLowering, &
                                calcSingleOverlapRaising, calcSingleOverlapMixed, &
                                calcDoubleLowering, calcDoubleRaising, &
                                calcDoubleL2R, calcDoubleR2L, calcFullstopLowering, &
                                calcFullstopRaising, calcFullStopL2R, &
                                calcFullStopR2L, calcFullStartLowering, &
                                calcFulLStartRaising, calcFullStartL2R, &
                                calcFullStartR2L, calcFullStartFullStopAlike, &
                                calcFullStartFullStopMixed, &
                                calcRemainingSwitches_excitInfo_double
    use guga_data, only: excitationInformation, tag_tmp_excits, tag_excitations
    use guga_types, only: weight_obj
    use guga_bitRepOps, only: update_matrix_element, setDeltaB, extract_matrix_element
    use guga_bitRepOps, only: isProperCSF_ilut, isDouble
    use guga_bitRepOps, only: write_guga_list, write_det_guga
    use guga_bitRepOps, only: convert_ilut_toGUGA, convert_ilut_toNECI, getDeltaB
    use MemoryManager, only: LogMemAlloc, LogMemDealloc
    use bit_reps, only: nifguga
    use FciMCData, only: projEDet, CurrentDets, TotWalkers
    use LoggingData, only: ThreshOccRDM, tThreshOccRDMDiag
    use UMatCache, only: gtID
    use RotateOrbsData, only: SymLabelListInv_rot
    use CalcData, only: tAdaptiveShift
    use Parallel_neci, only: nProcessors, MPIArg, MPIAlltoAll, MPIAlltoAllv
    use searching, only: BinSearchParts_rdm
    use rdm_data_utils, only: add_to_rdm_spawn_t, calc_combined_rdm_label, &
                              calc_separate_rdm_labels, extract_sign_rdm
    use OneEInts, only: GetTMatEl
    use procedure_pointers, only: get_umat_el
    use guga_matrixElements, only: calcDiagExchangeGUGA_nI

    implicit none

    ! test the symmetric filling of the GUGA-RDM, if the assumptions about 
    ! the hermiticity are correct..
    logical :: t_test_sym_fill = .true.
    logical :: t_test_diagonal = .false.

contains 

    subroutine gen_exc_djs_guga(ilutI)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        character(*), parameter :: this_routine = "gen_exc_djs_guga"

        integer :: nI(nel), flags_I, n_singles, n_doubles
        real(dp) :: sign_i(lenof_sign), full_sign(1)

        integer(n_int), pointer :: excits(:,:)
        integer(n_int) :: ilutG(0:nifguga)

        call extract_bit_rep(ilutI, nI, sign_I, flags_I)

        if (RDMExcitLevel == 1) then 
            call fill_diag_1rdm_guga(one_rdms, nI, sign_I)
        else
            full_sign = sign_i(1) * sign_I(lenof_sign)
            call fill_spawn_rdm_diag_guga(two_rdm_spawn, nI, full_sign)
        end if

        call convert_ilut_toGUGA(ilutI, ilutG)

        ! one-rdm is always calculated 
        ! calculate the excitations here
        ! but with my general two-body excitation routine i do not need 
        ! to calculate this in case of more than singles:
        if (RDMExcitLevel == 1) then
            call calc_explicit_1_rdm_guga(ilutG, n_singles, excits)

            ! and then sort them correctly in the communicated array
            call assign_excits_to_proc_guga(n_singles, excits, 1)

            deallocate(excits) 
            call LogMemDealloc(this_routine, tag_excitations)
        end if

        ! now to double excitations if requsted: 
        if (RDMExcitLevel /= 1) then 

            call calc_explicit_2_rdm_guga(ilutG, n_doubles, excits)

            call assign_excits_to_proc_guga(n_doubles, excits, 2)

            deallocate(excits) 
            call LogMemDealloc(this_routine, tag_excitations)
        end if

    end subroutine gen_exc_djs_guga

    subroutine fill_spawn_rdm_diag_guga(spawn, nI, full_sign) 
        ! i have to write a routine, which correctly takes into 
        ! account all the diagonal contributions from double excitations
        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: full_sign(spawn%rdm_send%sign_length)
        character(*), parameter :: this_routine = "fill_spawn_rdm_diag_guga"

        integer :: i, s_orbs(nel), s, inc_i, j, inc_j, p
        integer :: ssss, spsp, sp, sp2, a,b,c,d

        real(dp) :: occ_i, occ_j, x0, x1

        ! i have to figure out what exactly contributes to here..
        ! according to the paper about the two-RDM labels the 
        ! diagonal elements of the 2-RDM are given by 
        ! D_ij,kl = <psi| e_{ki,lj} | psi>
        ! which for the diagonal contribution kl = ij yields
        ! D_{ij,ij} = <psi| e_{ii,jj} |psi>
        ! which is definetly purely diagonal! 

        ! i could try to use the add_to_rdm_spawn_t routine in 
        ! rdm_data_utils with the according modified matrix elements. 
        i = 1
        s_orbs = gtID(nI)

        ! TODO: I could also try to add the diagonal exchange 
        ! contributions directly here! so this routine would be 
        ! applicable in the stochastic spawning already 
        do while (i <= nel) 

            s = s_orbs(i)

            if (isDouble(nI,i)) then 
                occ_i = 2.0_dp 
                inc_i = 2

                ! i think here i have to put the full diagonal 
                ! D_{ii,ii} entry, which counts double occupancies
                ! and i will try for now to reuse the following routine:
                ! put with the spatial orbital s
                call add_to_rdm_spawn_t(spawn, s, s, s, s, 2.0_dp * full_sign, .true.)

            else
                occ_i = 1.0_dp 
                inc_i = 1
            end if

            j = i + inc_i

            do while (j <= nel) 
                
                p = s_orbs(j)

                if (isDouble(nI,j)) then 
                    occ_j = 2.0_dp
                    inc_j = 2
                else
                    occ_j = 1.0_dp
                    inc_j = 1
                end if

                call add_to_rdm_spawn_t(spawn, s, s, p, p, & 
                    occ_i * occ_j * full_sign, .true.)

                ! i could also just multiply by 2 here, since this will 
                ! get strored in the same D_{ij,ij} RDM element! 
                ! but for now I decided to fill in all disting 2-RDM 
                ! elements
                call add_to_rdm_spawn_t(spawn, p, p, s, s, & 
                    occ_i * occ_j * full_sign, .true.)

                ! i can also add the fully diagonal exchange contributions 
                ! here. this is also necessary to do, if I want to use 
                ! this routine in the stochastic sampling

                ! but for open-shell to open-shell exchange excitations 
                ! I have to calculate the correct x1 matrix element..
                
                ! x0 matrix element is easy
                x0 = -occ_i * occ_j / 2.0_dp 

                if (inc_i == 1 .and. inc_j == 1) then 
                    ! if we have open-shell to open chell 
                    x1 = calcDiagExchangeGUGA_nI(i, j, nI) / 2.0_dp

                    call add_to_rdm_spawn_t(spawn, s, p, p, s, &
                       (x0 - x1) * full_sign, .true.)

                    call add_to_rdm_spawn_t(spawn, p, s, s, p, &
                       (x0 - x1) * full_sign, .true.)

                else
                    call add_to_rdm_spawn_t(spawn, s, p, p, s, &
                        x0 * full_sign, .true.)

                    ! and the symmetric version:
                    call add_to_rdm_spawn_t(spawn, p, s, s, p, &
                        x0 * full_sign, .true.)
                end if

                j = j + inc_j
            end do
            i = i + inc_i
        end do

    end subroutine fill_spawn_rdm_diag_guga

    subroutine fill_diag_1rdm_guga(one_rdms, nI, contrib_sign, &
                                    tCoreSpaceDetIn, RDMItersIn)
        ! I think I have to change these routines to also take into 
        ! account the possible coupling coefficient of doubly occupied 
        ! orbitals 
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: contrib_sign(:)
        logical, optional, intent(in) :: tCoreSpaceDetIn
        integer, optional, intent(in) :: RDMItersIn(:)
        character(*), parameter :: this_routine = "fill_diag_1rdm_guga"

        integer :: i, ind, irdm, s, s_orbs(nel), inc
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

        s_orbs = gtID(nI)

        i = 1

        do while (i <= nel) 
            
            s = s_orbs(i)
            
            ind = SymLabelListInv_rot(s)

            if (size(contrib_sign) == 1) then
                final_contrib = contrib_sign**2 * RDMIters * ScaleContribFac
            else
                final_contrib = contrib_sign(1::2) * contrib_sign(2::2) * RDMIters * ScaleContribFac
            end if

            if (isDouble(nI,i)) then
                inc = 2

                final_contrib = 2.0_dp * final_contrib

            else
                inc = 1

            end if

            if(tAdaptiveShift .and. all(nI == projEDet(:,1))) &
                 final_contrib = final_contrib + RDMIters * ScaleContribFac * rdmCorrectionFactor
            do irdm = 1, size(one_rdms)
                one_rdms(irdm)%matrix(ind,ind) = one_rdms(irdm)%matrix(ind,ind) + final_contrib(irdm)
            end do

            i = i + inc

        end do

    end subroutine fill_diag_1rdm_guga

    subroutine send_proc_ex_djs()
        ! i also need a specific routint to send the excitations for 
        ! the GUGA RDMs, otherwise this clutters up the det-based routines

        ! this routine is very similar to SendProcExcDjs in the rdm_explicit 
        ! module. see there for comments
        character(*), parameter :: this_routine = "send_proc_ex_djs"

        integer :: i, j
        integer :: error, MaxSendIndex,MaxIndex
        integer(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        integer(MPIArg) :: sing_recvcounts(nProcessors)
        integer(MPIArg) :: sing_recvdisps(nProcessors)
        integer(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        if (RDMExcitLevel == 1) then
            do i = 0, nProcessors - 1
                sendcounts(i+1) = int(Sing_ExcList(i)-(nint(OneEl_Gap*i)+1), MPIArg)
                disps(i+1) = nint(OneEl_Gap*i, MPIArg)
            end do

            MaxSendIndex = Sing_ExcList(nProcessors-1) - 1

            sing_recvcounts(1:nProcessors) = 0
            call MPIAlltoAll(sendcounts, 1, sing_recvcounts, 1, error)

            sing_recvdisps(1) = 0
            do i = 2, nProcessors
                sing_recvdisps(i) = sing_recvdisps(i-1) + sing_recvcounts(i-1)
            end do

            MaxIndex = sing_recvdisps(nProcessors) + sing_recvcounts(nProcessors)

            do i = 1, nProcessors
                sendcounts(i) = sendcounts(i)*(int(nifguga+1,MPIArg))
                disps(i) = disps(i)*(int(nifguga+1,MPIArg))
                sing_recvcounts(i) = sing_recvcounts(i)*(int(nifguga+1,MPIArg))
                sing_recvdisps(i) = sing_recvdisps(i)*(int(nifguga+1,MPIArg))
            end do


#ifdef PARALLEL
            call MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex), sendcounts, disps,&
                                Sing_ExcDjs2, sing_recvcounts, sing_recvdisps, error)
#else
            Sing_ExcDjs2(0:nifguga,1:MaxIndex) = Sing_ExcDjs(0:nifguga,1:MaxSendIndex)
#endif

            ! and also write a new routine for the search of occ. dets
            call singles_search_guga(sing_recvcounts, sing_recvdisps)
        end if

        if (RDMExcitLevel /= 1) then

            do i = 0, nProcessors-1
                ! Sendcounts is the number of singly excited determinants we
                ! want to send for each processor (but goes from 1, not 0).
                sendcounts(i+1) = int(Doub_ExcList(i)-(nint(TwoEl_Gap*i)+1),MPIArg)

                ! disps is the first slot for each processor - 1.            
                disps(i+1) = nint(TwoEl_Gap*i,MPIArg)
            end do

            MaxSendIndex = Doub_ExcList(nProcessors-1) - 1

            ! We now need to calculate the recvcounts and recvdisps -
            ! this is a job for AlltoAll
            doub_recvcounts(1:nProcessors) = 0
            call MPIAlltoAll(sendcounts, 1, doub_recvcounts, 1, error)

            ! We can now get recvdisps from recvcounts, since we want the data
            ! to be contiguous after the move.
            doub_recvdisps(1) = 0
            do i = 2, nProcessors
                doub_recvdisps(i) = doub_recvdisps(i-1) + doub_recvcounts(i-1)
            end do

            MaxIndex = doub_recvdisps(nProcessors)+doub_recvcounts(nProcessors)
            ! But the actual number of integers we need to send is the
            ! calculated values * NIfTot+1.
            do i = 1, nProcessors
                sendcounts(i) = sendcounts(i)*(int(nifguga+1,MPIArg))
                disps(i) = disps(i)*(int(nifguga+1,MPIArg))
                doub_recvcounts(i) = doub_recvcounts(i)*(int(nifguga+1,MPIArg))
                doub_recvdisps(i) = doub_recvdisps(i)*(int(nifguga+1,MPIArg))
            end do

            ! This is the main send of all the single excitations to the
            ! corresponding processors.
#ifdef PARALLEL
            call MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex), sendcounts, disps,&
                                    Doub_ExcDjs2, doub_recvcounts, doub_recvdisps, error)
#else
            Doub_ExcDjs2(0:nifguga,1:MaxIndex) = Doub_ExcDjs(0:nifguga,1:MaxSendIndex)
#endif

            call doubles_search_guga(doub_recvcounts, doub_recvdisps)

        end if

    end subroutine send_proc_ex_djs

    subroutine doubles_search_guga(recvcounts, recvdisps)
        ! this is the adapted explicit doubles search for the 
        ! GUGA-RDMs. 
        integer(MPIArg), intent(in) :: recvcounts(nProcessors), recvdisps(nProcessors)
        character(*), parameter :: this_routine = "doubles_search_guga"

        integer :: i, NoDets, StartDets, nI(nel), nJ(nel), PartInd, FlagsDj
        integer :: rdm_ind, j, ab, cd, a, b, c, d
        integer(n_int) :: ilutJ(0:niftot)
        real(dp) :: mat_ele, sign_i(lenof_sign), sign_j(lenof_sign)
        logical :: tDetFound

        do i = 1, nProcessors

            NoDets = recvcounts(i) / (nifguga + 1)
            StartDets = (recvdisps(i) / (nifguga + 1)) + 1

            if (NoDets > 1) then 

                sign_i = extract_matrix_element(Doub_ExcDjs2(:,StartDets), 2)

                do j = StartDets + 1, (NoDets + StartDets - 1)

                    call convert_ilut_toNECI(Doub_ExcDjs2(:,j), ilutJ)

                    call BinSearchParts_rdm(ilutJ, 1, int(TotWalkers, sizeof_int), &
                        PartInd, tDetFound)

                    if (tDetFound) then 

                        call extract_bit_rep(CurrentDets(:,PartInd), nJ, sign_j, FlagsDj)

                        mat_ele = extract_matrix_element(Doub_ExcDjs2(:,j), 1)
                        rdm_ind = getDeltaB(Doub_ExcDjs2(:,j))

                        call calc_separate_rdm_labels(rdm_ind, ab, cd, a, b, c, d)

                        call add_to_rdm_spawn_t(two_rdm_spawn, a, b, c, d, & 
                            sign_i * sign_j * mat_ele, .true.)

                        if (t_test_sym_fill) then 
                            ! i only calculate excitations with ab < cd in this case
                            ! so fill in the missing ones
                            call add_to_rdm_spawn_t(two_rdm_spawn, b, a, d, c, & 
                                sign_i * sign_j * mat_ele, .true.)


                        end if
                    end if
                end do
            end if
        end do

    end subroutine doubles_search_guga

    subroutine singles_search_guga(recvcounts, recvdisps)
        ! make a tailored GUGA search for explicit single excitations
        
        ! this routine is very similar to Sing_SearchOccDets in the rdm_explicit
        ! module. see there for comments
        integer(MPIArg), intent(in) :: recvcounts(nProcessors), recvdisps(nProcessors)
        character(*), parameter :: this_routine = "singles_search_guga"

        integer :: i, NoDets, StartDets, nI(nel), nJ(nel), PartInd, FlagsDj
        integer :: rdm_ind, j
        integer(n_int) :: ilutJ(0:niftot), ilutG
        real(dp) :: mat_ele, sign_i(lenof_sign), sign_j(lenof_sign)
        logical :: tDetFound


        do i = 1, nProcessors

            NoDets = recvcounts(i)/(nifguga + 1)
            StartDets = (recvdisps(i) / (nifguga + 1)) + 1

            if (NoDets > 1) then 
                
                ! extract the determinant, the matrix element and the 
                ! RDM index from the single excitation array

                ! the convention is: i encode the coupling coefficient 
                ! in the x0 matrix element and the sign of Di in the 
                ! x1 element. and the combined RDM index in the 
                ! Delta B value position! 
                sign_i = extract_matrix_element(Sing_ExcDjs2(:,StartDets),2)

                do j = StartDets + 1, (NoDets + StartDets - 1)

                    ! apparently D_i is in the first spot and all 
                    ! excitations come here afterwards.. 

                    call convert_ilut_toNECI(Sing_ExcDjs2(:,j), ilutJ)

                    call BinSearchParts_rdm(ilutJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    if (tDetFound) then 

                        call extract_bit_rep(CurrentDets(:,PartInd), nJ, sign_j, FlagsDj)

                        mat_ele = extract_matrix_element(Sing_ExcDjs2(:,j), 1)
                        rdm_ind = getDeltaB(Sing_ExcDjs2(:,j))

                        if (RDMExcitLevel == 1) then 
                            call fill_sings_1rdm_guga(one_rdms, sign_i, sign_j, &
                                mat_ele, rdm_ind, t_test_sym_fill)
                        end if
                    end if
                end do
            end if
        end do

    end subroutine singles_search_guga

    subroutine fill_sings_1rdm_guga(one_rdms, sign_i, sign_j, mat_ele, rdm_ind, &
            t_fill_symmetric)
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        integer, intent(in) :: rdm_ind
        logical, intent(in) :: t_fill_symmetric
        character(*), parameter :: this_routine = "fill_sings_1rdm_guga"

        integer :: i, a, ind_i, ind_a, irdm

        call extract_1_rdm_ind(rdm_ind, i, a)

        ind_i = SymLabelListInv_rot(i)
        ind_a = SymLabelListInv_rot(a)

        do irdm = 1, size(one_rdms)

            one_rdms(irdm)%matrix(ind_a, ind_i) = one_rdms(irdm)%matrix(ind_a, ind_i) & 
                + sign_i(irdm) * sign_j(irdm) * mat_ele

            if (t_fill_symmetric) then 

                one_rdms(irdm)%matrix(ind_i, ind_a) = one_rdms(irdm)%matrix(ind_i, ind_a) & 
                    + sign_i(irdm) * sign_j(irdm) * mat_ele
            end if
        end do

    end subroutine fill_sings_1rdm_guga

    subroutine extract_1_rdm_ind(rdm_ind, i, a)
        ! the converstion routine between the combined and explicit rdm 
        ! indices for the 1-RDM
        integer, intent(in) :: rdm_ind
        integer, intent(out) :: i, a
        character(*), parameter :: this_routine = "extract_matrix_element"

        a = mod(rdm_ind - 1, nSpatOrbs)  + 1
        i = (rdm_ind - 1)/nSpatOrbs + 1

    end subroutine extract_1_rdm_ind

    function contract_1_rdm_ind(i,a) result(rdm_ind) 
        ! the inverse function of the routine above, to give the combined 
        ! rdm index of two explicit ones
        integer, intent(in) :: i, a
        integer :: rdm_ind
        character(*), parameter :: this_routine = "contract_1_rdm_ind"

        rdm_ind = nSpatOrbs * (i - 1) + a

    end function contract_1_rdm_ind

    function contract_2_rdm_ind(i,j,k,l) result(ijkl)
        ! since I only ever have spatial orbitals in the GUGA-RDM make 
        ! the definition of the RDM-index combination differently
        integer, intent(in) :: i,j,k,l
        integer(int_rdm) :: ijkl
        character(*), parameter :: this_routine = "contract_2_rdm_ind"

        integer :: ij, kl

        ij = (i - 1) * nSpatOrbs + j 
        kl = (k - 1) * nSpatOrbs + l 

        ijkl = (ij - 1) * (int(nSpatOrbs, int_rdm)**2) + kl

    end function contract_2_rdm_ind

    subroutine extract_2_rdm_ind(ijkl, i, j, k, l, ij_out, kl_out)
        ! the inverse routine of the function above. 
        ! it is actually practical to have ij and kl also available at 
        ! times, since it can identify diagonal entries of the two-RDM
        integer(int_rdm), intent(in) :: ijkl
        integer, intent(out) :: i,j,k,l
        integer, intent(out), optional :: ij_out, kl_out
        character(*), parameter :: this_routine = "extract_2_rdm_ind"

        integer :: ij, kl

        kl = mod(ijkl - 1, int(nSpatOrbs, int_rdm)**2) + 1
        ij = (ijkl - kl)/(nSpatOrbs ** 2) + 1

        j = mod(ij - 1, nSpatOrbs) + 1
        i = (ij - j) / nSpatOrbs + 1 

        l = mod(kl - 1, nSpatOrbs) + 1 
        k = (kl - l) / nSpatOrbs + 1

        if (present(ij_out)) ij_out = ij
        if (present(kl_out)) kl_out = kl

    end subroutine extract_2_rdm_ind


    subroutine assign_excits_to_proc_guga(n_tot, excits, excit_lvl)
        integer, intent(in) :: n_tot, excit_lvl
        integer(n_int), intent(in), pointer :: excits(:,:)
        character(*), parameter :: this_routine = "assign_excits_to_proc_guga"

        integer :: i, proc, nJ(nel)
        integer(n_int) :: ilut(0:niftot)

        ASSERT(excit_lvl == 1 .or. excit_lvl == 2)

        if (excit_lvl == 1) then 
            do i = 1, n_tot

                call convert_ilut_toNECI(excits(:,i), ilut)

                call decode_bit_det(nJ, ilut)

                proc = DetermineDetNode(nel, nJ, 0)

                Sing_ExcDjs(:, Sing_ExcList(proc)) = excits(:,i)
                Sing_ExcList(proc) = Sing_ExcList(proc) + 1

#ifdef __DEBUG
                if (Sing_ExcList(Proc) .gt. nint(OneEl_Gap*(Proc+1))) then
                    write(6,*) 'Proc', Proc
                    write(6,*) 'Sing_ExcList', Sing_ExcList
                    write(6,*) 'No. spaces for each proc', nint(OneEl_Gap)
                    call Stop_All('GenExcDjs', 'Too many excitations for space available.')
                end if
#endif
            end do

        else if (excit_lvl == 2) then 
            do i = 1, n_tot

                call convert_ilut_toNECI(excits(:,i), ilut)

                call decode_bit_det(nJ, ilut)
                proc = DetermineDetNode(nel, nJ, 0)

                Doub_ExcDjs(:, Doub_ExcList(proc)) = excits(:,i)
                Doub_ExcList(proc) = Doub_ExcList(proc) + 1

#ifdef __DEBUG
                if (Doub_ExcList(Proc) .gt. nint(TwoEl_Gap*(Proc+1))) then
                    write(6,*) 'Proc', Proc
                    write(6,*) 'Doub_ExcList', Doub_ExcList
                    write(6,*) 'No. spaces for each proc', nint(TwoEl_Gap)
                    call Stop_All('GenExcDjs','Too many excitations for space available.')
                end if
#endif
            end do
        end if

    end subroutine assign_excits_to_proc_guga

    subroutine calc_explicit_2_rdm_guga(ilut, n_tot, excitations)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), pointer :: excitations(:,:)
        character(*), parameter :: this_routine = "calc_explicit_2_rdm_guga"

        integer :: i, j, k, l, nMax, ierr, n, n_excits, jl,ik
        integer(n_int), pointer :: temp_excits(:,:), tmp_all_excits(:,:)
        integer(int_rdm) :: ijkl

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**3 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:nifguga,nMax), stat = ierr)
        call LogMemAlloc('tmp_all_excits',(nifguga+1)*nMax,8,this_routine,tag_tmp_excits)

        n_tot = 0

        if (t_test_diagonal) then 
            do i = 1, nSpatOrbs 
                do j = 1, nSpatOrbs

                    if (i == j) cycle

                    call calc_all_excits_guga_rdm_doubles(ilut, j, i, i, j, &
                        temp_excits, n_excits)

#ifdef __DEBUG
                    do n = 1, n_excits
                        ASSERT(isProperCSF_ilut(temp_excits(:,n), .true.))
                    end do
#endif
                    if (n_excits > 1) then 
                        call add_guga_lists_rdm(n_tot, n_excits-1, & 
                            tmp_all_excits, temp_excits(:,2:))
                    end if

                    deallocate(temp_excits)
                end do
            end do

        else if (.not. t_test_sym_fill) then
            do i = 1, nSpatOrbs
                do j = 1, nSpatOrbs
                    do k = 1, nSpatOrbs
                        do l = 1, nSpatOrbs

                            ! pure diagonal contribution:
                            if (i == j .and. k == l) cycle

                            call calc_all_excits_guga_rdm_doubles(ilut, i, j, k, l, &
                                temp_excits, n_excits)

#ifdef __DEBUG
                            do n = 1, n_excits
                                ASSERT(isProperCSF_ilut(temp_excits(:,n), .true.))
                            end do
#endif
                            if (i == l .and. j == k) then 
                                ! exclude the diagonal exchange here, 
                                ! as it is already accounted for in the 
                                ! diagonal contribution routine
#ifdef __DEBUG
                                if (n_excits > 0) then
                                    if (.not. DetBitEQ(ilut, temp_excits(:,1), nifdbo)) then
                                        print *, "not equal!"
                                    end if
                                end if
#endif

                                if (n_excits > 1) then
                                    call add_guga_lists_rdm(n_tot, n_excits - 1, &
                                        tmp_all_excits, temp_excits(:,2:))
                                end if
                            else
                                if (n_excits > 0) then 
                                    call add_guga_lists_rdm(n_tot, n_excits, tmp_all_excits, temp_excits)
                                end if
                            end if

                            deallocate(temp_excits)

                        end do
                    end do
                end do
            end do
        else
            do i = 1, nSpatOrbs
                do j = 1, nSpatOrbs
                    do k = 1, nSpatOrbs
                        do l = 1, nSpatOrbs

                            if (i == j .and. k == l) cycle

                            call calc_combined_rdm_label(j,l,i,k,ijkl,jl,ik)

                            if (jl > ik) cycle

                            call calc_all_excits_guga_rdm_doubles(ilut, i, j, k, l, &
                                temp_excits, n_excits)

#ifdef __DEBUG
                            do n = 1, n_excits
                                ASSERT(isProperCSF_ilut(temp_excits(:,n), .true.))
                            end do
#endif
                            if (i == l .and. j == k) then 
#ifdef __DEBUG
                                if (n_excits > 0) then
                                    if (.not. DetBitEQ(ilut, temp_excits(:,1), nifdbo)) then
                                        print *, "not equal!"
                                    end if
                                end if
#endif
                                if (n_excits > 1) then

                                    call add_guga_lists_rdm(n_tot, n_excits - 1, &
                                        tmp_all_excits, temp_excits(:,2:))
                                end if
                            else
                                if (n_excits > 0) then 
                                    call add_guga_lists_rdm(n_tot, n_excits, tmp_all_excits, temp_excits)
                                end if
                            end if

                            deallocate(temp_excits)

                        end do
                    end do
                end do
            end do
        end if

        j = 1
        do i = 1, n_tot
            if (abs(extract_matrix_element(tmp_all_excits(:,i),1)) < EPS) cycle

            tmp_all_excits(:,j) = tmp_all_excits(:,i)

            j = j + 1

        end do

        n_tot = j - 1

        allocate(excitations(0:nifguga,n_tot), stat = ierr) 
        ! hm to log that does not make so much sense.. since it gets called 
        ! more than once and is only a temporary array..
        call LogMemAlloc('excitations',n_tot,8,this_routine,tag_excitations)

        excitations = tmp_all_excits(:,1:n_tot)

        deallocate(tmp_all_excits) 
        call LogMemDealloc(this_routine, tag_tmp_excits)

    end subroutine calc_explicit_2_rdm_guga

    subroutine calc_explicit_1_rdm_guga(ilut, n_tot, excitations)
        ! routine to calculate the one-RDM explicitly in the GUGA formalism
        ! the one RDM is given by rho_ij = <Psi|E_ij|Psi> 
        ! so, similar to the actHamiltonian routine I need to loop over all 
        ! the orbital indices (i,j) and calculate the action of a given E_ij 
        ! on the sampled wavefunction |Psi> and the resulting overlap with <Psi|
        ! make this routine similar to GenExcDjs() in rdm_explicit.F90
        ! to insert it there to calculate the GUGA RDMs in this case
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), pointer :: excitations(:,:)
        character(*), parameter :: this_routine = "calc_explicit_1_rdm_guga"

        integer :: i, j, nMax, ierr, n, n_excits
        integer(n_int), pointer :: temp_excits(:,:), tmp_all_excits(:,:)

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**2 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:nifguga,nMax), stat = ierr)
        call LogMemAlloc('tmp_all_excits',(nifguga+1)*nMax,8,this_routine,tag_tmp_excits)

        n_tot = 0

        if (.not. t_test_sym_fill) then 
            do i = 1, nSpatOrbs
                do j = 1, nSpatOrbs 
                    if (i == j) cycle

                    call calc_all_excits_guga_rdm_singles(ilut, i, j, temp_excits, &
                        n_excits)

#ifdef __DEBUG
                    do n = 1, n_excits
                        ASSERT(isProperCSF_ilut(temp_excits(:,n), .true.))
                    end do
#endif

                    if (n_excits > 0) then 
                        call add_guga_lists_rdm(n_tot, n_excits, tmp_all_excits, temp_excits)
                    end if

                    deallocate(temp_excits)

                end do
            end do

        else
            do i = 1, nSpatOrbs - 1
                do j = i + 1, nSpatOrbs

                    call calc_all_excits_guga_rdm_singles(ilut, i, j, temp_excits, &
                        n_excits)

#ifdef __DEBUG
                    do n = 1, n_excits
                        ASSERT(isProperCSF_ilut(temp_excits(:,n), .true.))
                    end do
#endif

                    if (n_excits > 0) then 
                        call add_guga_lists_rdm(n_tot, n_excits, tmp_all_excits, temp_excits)
                    end if

                    deallocate(temp_excits)

                end do
            end do
        end if

        j = 1
        do i = 1, n_tot
            if (abs(extract_matrix_element(tmp_all_excits(:,i),1)) < EPS) cycle

            tmp_all_excits(:,j) = tmp_all_excits(:,i)

            j = j + 1

        end do

        n_tot = j - 1

        allocate(excitations(0:nifguga,n_tot), stat = ierr) 
        ! hm to log that does not make so much sense.. since it gets called 
        ! more than once and is only a temporary array..
        call LogMemAlloc('excitations',n_tot,8,this_routine,tag_excitations)

        excitations = tmp_all_excits(:,1:n_tot)

        deallocate(tmp_all_excits) 
        call LogMemDealloc(this_routine, tag_tmp_excits)


    end subroutine calc_explicit_1_rdm_guga


    subroutine calc_all_excits_guga_rdm_doubles(ilut, i, j, k, l, excits, n_excits)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: i, j, k, l
        integer(n_int), intent(out), pointer :: excits(:,:)
        integer, intent(out) :: n_excits
        character(*), parameter :: this_routine = "calc_all_excits_guga_rdm_doubles"

        real(dp) :: umat, posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer :: ierr, n, exlevel
        type(excitationInformation) :: excitInfo
        logical :: compFlag
        integer(n_int) :: tmp_ilut(0:niftot)
        integer(int_rdm) :: ijkl

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        ! if called with k = l = 0 -> call single version of function
        if (k == 0 .and. l == 0) then
            call calc_all_excits_guga_rdm_singles(ilut, i, j, excits, n_excits)
            return
        end if

        ASSERT(k > 0 .and. k <= nSpatOrbs)
        ASSERT(l > 0 .and. l <= nSpatOrbs)

        ! default:
        n_excits = 0

        ! otherwise get excitation information
        excitInfo = excitationIdentifier(i, j, k, l)
       
        ! screw it. for now write a function which checks if indices and ilut
        ! are compatible, and not initiate the excitation right away
        ! but check cases again 
        ! in checkCompatibility the number of switches is already 
        ! calulated to check if the probabilistic weights fit... maybe but 
        ! that out and reuse.. to not waste any effort.
        call checkCompatibility(ilut, excitInfo, compFlag, posSwitches, negSwitches)


        ! for mixed type full starts and/or full stops i have to consider 
        ! the possible diagonal/single excitations here!
        ! although I think the full-stops are already correctly taken into 
        ! account..
        ! and also the the full-starts are maybe correct already.. 
        ! so it was just the full-start into full-stop mixed!
        if (.not.compFlag) then !.and. .not. excitInfo%typ == 23) then
            allocate(excits(0,0), stat = ierr)
            return
        end if

        select case(excitInfo%typ)
        case(0)
            ! shouldnt be here.. onyl single excits and full weight gens
            allocate(excits(0,0), stat = ierr)
            return

        case(1) ! weight + lowering gen.
            ! can be treated almost like a single excitation 
            ! essentially the same, except if d(w) == 3 in the excitaton regime
            call calcDoubleExcitation_withWeight(ilut, excitInfo, excits,&
                n_excits, posSwitches, negSwitches)

            exlevel = 1

        case(2) ! weight + raising gen
            call calcDoubleExcitation_withWeight(ilut, excitInfo, excits,&
                n_excits, posSwitches, negSwitches)

            exlevel = 1

        case(3) ! non overlap 
            call calcNonOverlapDouble(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case(4) ! single overlap two lowering
            ! how can i efficiently adress that?
            ! can i write that efficiently in one function or do i need more?
            ! probably need more... i already determined
            call calcSingleOverlapLowering(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 1

        case(5) ! single overlap raising 
            call calcSingleOverlapRaising(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 1

        case (6) ! single overlap lowering into raising
            call calcSingleOverlapMixed(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (7) ! single overlap raising into lowering
            call calcSingleOverlapMixed(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (8) ! normal double overlap two lowering
            call calcDoubleLowering(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (9) ! normal double overlap two raising
            call calcDoubleRaising(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (10) ! lowering into raising into lowering
            call calcDoubleRaising(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (11) ! raising into lowering into raising
            call calcDoubleLowering(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (12) ! lowering into raising double
            call calcDoubleL2R(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (13) ! raising into lowering double
            call calcDoubleR2L(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (14) ! full stop 2 lowering
            ! can i write a function for both alike generator combinations
            ! i think i can
            call calcFullstopLowering(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (15) ! full stop 2 raising
            call calcFullstopRaising(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (16) ! full stop lowering into raising 
            call calcFullStopL2R(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            ! in this case there is also the possibility for one single-like 
            ! excitation if there is no change in the double overlap region! 
            ! todo! how to fix that? is that so important? its only max. 1 
            exlevel = 2

        case (17) ! full stop raising into lowering 
            call calcFullStopR2L(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        case (18) ! full start 2 lowering
            call calcFullStartLowering(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (19) ! full start 2 raising
            call calcFulLStartRaising(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            exlevel = 2

        case (20) ! full start lowering into raising
            call calcFullStartL2R(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        case (21) ! full start raising into lowering
            call calcFullStartR2L(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        case (22) ! full start into full stop alike
            call calcFullStartFullStopAlike(ilut, excitInfo, excits)
            n_excits = 1

            exlevel = 2

        case (23) ! full start into full stop mixed
            call calcFullStartFullStopMixed(ilut, excitInfo, excits, n_excits, &
                posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        end select


        ! indicate the level of excitation IC for the remaining NECI code
        if (n_excits > 0) then 
            call calc_combined_rdm_label(i,j,k,l, ijkl)
            excits(nifguga,:) = ijkl
        end if

    end subroutine calc_all_excits_guga_rdm_doubles

    subroutine calc_all_excits_guga_rdm_singles(ilut, i, j, excits, n_excits)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: i, j
        integer(n_int), intent(out), pointer :: excits(:,:)
        integer, intent(out) :: n_excits
        character(*), parameter :: this_routine = "calc_all_excits_guga_rdm_singles"

        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer :: st, iEx, iOrb, ierr
        integer(n_int), pointer :: tempExcits(:,:)
        real(dp) :: minusWeight, plusWeight

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        n_excits = 0

        excitInfo = excitationIdentifier(i,j)

        if (excitInfo%gen1 /= 0 ) then
            if (current_stepvector(i) == 3 .or. current_stepvector(j) == 0) then
                allocate(excits(0,0), stat = ierr)
                return
            end if
        end if

        call calcRemainingSwitches(ilut, excitInfo, posSwitches, negSwitches)

        weights = init_singleWeight(ilut, excitInfo%fullEnd)
        plusWeight = weights%proc%plus(posSwitches(excitInfo%fullStart), &
            currentB_ilut(excitInfo%fullStart), weights%dat)
        minusWeight = weights%proc%minus(negSwitches(excitInfo%fullStart),&
            currentB_ilut(excitInfo%fullStart), weights%dat)

        st = excitInfo%fullStart
        ! check compatibility of chosen indices

        if ((current_stepvector(st) == 1 .and. plusWeight < EPS) .or.&
            (current_stepvector(st) == 2 .and. minusWeight < EPS).or.&
            (minusWeight + plusWeight < EPS)) then
            allocate(excits(0,0), stat = ierr)
            return
        end if


        ! have to give probabilistic weight object as input, to deal 
        call createSingleStart(ilut, excitInfo, posSwitches, &
            negSwitches, weights, tempExcits, n_excits)

        do iOrb = excitInfo%fullStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, iOrb, excitInfo, posSwitches, &
                negSwitches, weights, tempExcits, n_excits)
        end do

        call singleEnd(ilut, excitInfo, tempExcits, &
            n_excits, excits)

        ! encode the combined RDM-ind in the deltaB position for 
        ! communication purposes
        excits(nifguga,:) = contract_1_rdm_ind(i,j)

    end subroutine calc_all_excits_guga_rdm_singles

    subroutine add_guga_lists_rdm(n_dets_tot, n_dets, list_tot, list)
        ! i need a new guga list addition for RDMs as I must not add 
        ! identical excitations coming from different RDM index combinations
        ! otherwise I loose the information of those excitations!
        integer, intent(inout) ::        n_dets_tot
        integer, intent(in) ::           n_dets
        integer(n_int), intent(inout) :: list_tot(0:,1:), list(0:,1:)
        character(*), parameter :: this_routine = "add_guga_lists_rdm"

        ! here I essentially only need to append the list to the total 
        ! list and add the number of new elements to n_dets_tot
 
        list_tot(:,n_dets_tot+1:n_dets_tot+n_dets) = list
        n_dets_tot = n_dets_tot + n_dets

    end subroutine add_guga_lists_rdm

    subroutine calc_rdm_energy_guga(rdm, rdm_energy_1, rdm_energy_2)
        ! the GUGA RDM version of the energy calculator

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(out) :: rdm_energy_1(rdm%sign_length)
        real(dp), intent(out) :: rdm_energy_2(rdm%sign_length)
        character(*), parameter :: this_routine = "calc_rdm_energy_guga"

        integer(int_rdm) :: ijkl
        integer :: ielem, ij, kl, i, j, k, l
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_energy_1 = 0.0_dp
        rdm_energy_2 = 0.0_dp

        ! Loop over all elements in the 2-RDM.
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

            ! Decode pqrs label into p, q, r and s labels.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

            ! D_{ij,kl} corresponds to V_{ki,lj} * e_{ki,lj} i believe..
            ! i have to figure out if the problem is here or somewhere else..
            rdm_energy_2 = rdm_energy_2 + rdm_sign * get_umat_el(k,l,i,j)

            if (.not. t_test_sym_fill) then
                if (i == k) then
                    rdm_energy_1 = rdm_energy_1 + rdm_sign * GetTMatEl(2*j,2*l) / real(nel-1,dp)
                end if
                if (j == l) then 
                    rdm_energy_1 = rdm_energy_1 + rdm_sign * GetTMatEl(2*i,2*k) / real(nel-1,dp)
                end if
            else
                if (j == l) then 
                    rdm_energy_1 = rdm_energy_1 + 2.0_dp * rdm_sign * GetTMatEl(2*i,2*k) / real(nel-1,dp)
                end if
            end if
        end do

    end subroutine calc_rdm_energy_guga

end module guga_rdm

#endif
