#include "macros.h"
#ifndef __CMPLX

module guga_rdm
    ! RDM module specifically for the GUGA spin-adapted implementation 
    
    use constants, only: n_int, dp, lenof_sign, EPS
    use bit_rep_data, only: niftot
    use SystemData, only: nel, nSpatOrbs, current_stepvector, currentB_ilut
    use bit_reps, only: extract_bit_rep
    use LoggingData, only: RDMExcitLevel
    use rdm_filling, only: fill_diag_1rdm, fill_spawn_rdm_diag
    use rdm_data, only: one_rdms, two_rdm_spawn
    use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs
    use rdm_data, only: Sing_ExcList, Doub_ExcList, OneEl_Gap, TwoEl_Gap
    use DetBitOps, only: EncodeBitDet, count_open_orbs
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
                                calcFullStartFullStopMixed
    use guga_data, only: excitationInformation, tag_tmp_excits, tag_excitations
    use guga_types, only: weight_obj
    use guga_bitRepOps, only: update_matrix_element, setDeltaB, extract_matrix_element
    use guga_bitRepOps, only: add_guga_lists, isProperCSF_ilut
    use MemoryManager, only: LogMemAlloc, LogMemDealloc
    use bit_reps, only: nifguga

    implicit none

contains 

    subroutine gen_exc_djs_guga(ilutI)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        character(*), parameter :: this_routine = "gen_exc_djs_guga"

        integer :: nI(nel), flags_I, n_singles, n_doubles
        real(dp) :: sign_i(lenof_sign), full_sign(1)

        integer(n_int), pointer :: excits(:,:)

        call extract_bit_rep(ilutI, nI, sign_I, flags_I)

        if (RDMExcitLevel == 1) then 
            call fill_diag_1rdm(one_rdms, nI, sign_I)
        else
            full_sign = sign_i(1) * sign_I(lenof_sign)
            call fill_spawn_rdm_diag(two_rdm_spawn, nI, full_sign)
        end if

        ! one-rdm is always calculated 
        ! calculate the excitations here
        call calc_explicit_1_rdm_guga(ilutI, n_singles, excits)

        ! and then sort them correctly in the communicated array
        call assign_excits_to_proc_guga(n_singles, excits, 1)

        ! now to double excitations if requsted: 
        if (RDMExcitLevel /= 1) then 

            deallocate(excits)

            call calc_explicit_2_rdm_guga(ilutI, n_doubles, excits)

            call assign_excits_to_proc_guga(n_doubles, excits, 2)

        end if

    end subroutine gen_exc_djs_guga

    subroutine assign_excits_to_proc_guga(n_tot, excits, excit_lvl)
        integer, intent(in) :: n_tot, excit_lvl
        integer(n_int), intent(in), pointer :: excits(:,:)
        character(*), parameter :: this_routine = "assign_excits_to_proc_guga"

        integer :: i, proc, nJ(nel)

        ASSERT(excit_lvl == 1 .or. excit_lvl == 2)

        if (excit_lvl == 1) then 
            do i = 1, n_tot

                call EncodeBitDet(nJ, excits(:,i))

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

                call EncodeBitDet(nJ, excits(:,i))
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
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), pointer :: excitations(:,:)
        character(*), parameter :: this_routine = "calc_explicit_2_rdm_guga"

        integer :: i, j, k, l, nMax, ierr, n, n_excits
        integer(n_int), pointer :: temp_excits(:,:), tmp_all_excits(:,:)

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**4 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:nifguga,nMax), stat = ierr)
        call LogMemAlloc('tmp_all_excits',(nifguga+1)*nMax,8,this_routine,tag_tmp_excits)

        n_tot = 0

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

                        ! exclude possible diagonal terms, but do not 
                        ! forget to account for those in the diagonal 
                        ! contributions! 
                        if (i == l .and. k == j) then 
                            if (n_excits > 1) then 
                                call add_guga_lists(n_tot, n_excits - 1, &
                                    tmp_all_excits, temp_excits(:,2:))
                            end if
                        else
                            if (n_excits > 0) then 
                                call add_guga_lists(n_tot, n_excits, tmp_all_excits, temp_excits)
                            end if
                        end if

                    end do
                end do
            end do
        end do

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
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), pointer :: excitations(:,:)
        character(*), parameter :: this_routine = "calc_explicit_1_rdm_guga"

        integer :: i, j, nMax, ierr, n, n_excits
        integer(n_int), pointer :: temp_excits(:,:), tmp_all_excits(:,:)

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**4 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:nifguga,nMax), stat = ierr)
        call LogMemAlloc('tmp_all_excits',(nifguga+1)*nMax,8,this_routine,tag_tmp_excits)

        n_tot = 0

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
                    call add_guga_lists(n_tot, n_excits, tmp_all_excits, temp_excits)
                end if
            end do
        end do

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
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: i, j, k, l
        integer(n_int), intent(out), pointer :: excits(:,:)
        integer, intent(out) :: n_excits
        character(*), parameter :: this_routine = "calc_all_excits_guga_rdm_doubles"

        real(dp) :: umat, posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer :: ierr, n, exlevel
        type(excitationInformation) :: excitInfo
        logical :: compFlag

        integer(n_int) :: tmp_ilut(0:niftot)

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

        if (.not.compFlag) then
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
        excits(nifguga,:) = exlevel 

    end subroutine calc_all_excits_guga_rdm_doubles

    subroutine calc_all_excits_guga_rdm_singles(ilut, i, j, excits, n_excits)
        integer(n_int), intent(in) :: ilut(0:niftot)
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

        ! encode IC = 1 in the deltB information of the GUGA 
        ! excitation to handle it in the remaining NECI code 
        ! correctly 
        do iEx = 1, n_excits
            call setDeltaB(1, excits(:,iEx))
        end do

    end subroutine calc_all_excits_guga_rdm_singles


end module guga_rdm

#endif
