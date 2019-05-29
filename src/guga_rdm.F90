#include "macros.h"
#ifndef __CMPLX

module guga_rdm
    ! RDM module specifically for the GUGA spin-adapted implementation 
    
    use constants, only: n_int, dp
    use bit_rep_data, only: niftot
    use SystemData, only:

    implicit none

contains 


    subroutine gen_exc_djs_guga(ilutI)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        character(*), parameter :: this_routine = "gen_exc_djs_guga"

        integer(n_int) :: ilutJ(0:niftot)
        integer :: nI(nel), nJ(nel), proc, flags_I
        real(dp) :: sign_i(lenof_sign), full_sign(1)

        call extract_bit_rep(ilutI, nI, sign_I, flags_I)

        if (RDMExcitLevel == 1) then 
            call fill_diag_1rdm(one_rdms, nI, sign_I)
        else
            full_sign = sign_i(1) * sign_I(lenof_sign)
            call fill_spawn_rdm_diag(two_rdm_spawn, nI, full_sign)
        end if

        ! one-rdm is always calculated 
        ! calculate the excitations here
        call calc_explicit_1_rdm(ilutI, n_tot, excits)

        ! and then sort them correctly in the communicated array
        call assign_excits_to_proc_guga(n_tot, excits)

    end subroutine gen_exc_djs_guga

    subroutine assign_excits_to_proc_guga(n_tot, excits, excit_lvl)
        integer, intent(in) :: n_tot, excit_lvl
        integer(n_int), intent(in), pointer :: excits(:,:)
        character(*), parameter :: this_routine = "assign_excits_to_proc_guga"

        integer :: i, proc

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

    end subroutine assign_excits_to_proc_guga


    subroutine calc_explicit_1_rdm(ilut, n_tot, excitations)
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

        character(*), parameter :: this_routine = "calc_explicit_1_rdm"

        integer :: i, j

        call init_csf_information(ilut)

        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs 
                if (i == j) cycle

                call calc_all_excits_guga_rdm_singles(ilut, i, j, temp_excits, &
                    n_excits)

#ifdef __DEBUG
                do n = 1, n_excits
                    ASSERT(isProperCSF_ilut(temp_excits(:,n), .true.))
#endif

                if (n_excits > 0) then 
                    call add_guga_lists(n_tot, n_excits, tmp_all_excits, temp_excits)
                end if
            end do
        end do

    end subroutine calc_explicit_1_rdm

    subroutine calc_all_excits_guga_rdm_singles(ilut, i, j, excits, n_excits)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: i, j
        integer(n_int), intent(out), pointer :: excits(:,:)
        integer, intent(out) :: n_excits
        character(*), parameter :: this_routine = "calc_all_excits_guga_rdm_singles"

        type(excitationInformation) :: excitInfo

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        n_excits = 0

        excitInfo = excitationIdentifier(i,j)

        if (excitInfo%gen1 /= 0 ) then
            if (current_stepvector(i) == 3 .or. current_stepvector(j) == 0) then
                allocate(excitations(0,0), stat = ierr)
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
            allocate(excitations(0,0), stat = ierr)
            return
        end if


        ! have to give probabilistic weight object as input, to deal 
        call createSingleStart(ilut, excitInfo, posSwitches, &
            negSwitches, weights, tempExcits, nExcits)

        ! to not call getTmatEl again in createSingleStart loop over 
        ! the atmost two excitations here and multiply with tmat

        do iEx = 1, nExcits
            call update_matrix_element(tempExcits(:,iEx), tmat, 1)
        end do


        do iOrb = excitInfo%fullStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, iOrb, excitInfo, posSwitches, &
                negSwitches, weights, tempExcits, nExcits)
        end do

        call singleEnd(ilut, excitInfo, tempExcits, &
            nExcits, excitations)

        ! encode IC = 1 in the deltB information of the GUGA 
        ! excitation to handle it in the remaining NECI code 
        ! correctly 
        do iEx = 1, nExcits
            call setDeltaB(1, excitations(:,iEx))
        end do

    end subroutine calc_all_excits_guga_rdm_singles


end module guga_rdm

#endif
