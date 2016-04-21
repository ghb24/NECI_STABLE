#include "macros.h"
! GUGA bit-representation operations:
! contains functions, concerning the representation of CSFs and 
! associated information in NECI and operations thereon. 
! some of these functions should be replaced by bit-ops defined in the 
! macros.h header file, but for clarity and test purposes write them down
! explicetly too
#ifndef __CMPLX
module guga_bitRepOps

    use SystemData, only: nEl, Stot, nSpatOrbs
    use guga_data ! get explicit here too!
    use bit_reps, only: niftot, set_flag, nIfGUGA, nIfD, nifdbo
    use constants, only: dp, n_int, bits_n_int, bni_, bn2_
    use DetBitOps, only: return_ms, count_set_bits, MaskAlpha, &
                    count_open_orbs, ilut_lt, ilut_gt, MaskAlpha, MaskBeta, &
                    CountBits
    use bit_rep_data, only: test_flag, flag_deltaB_single, &
        flag_deltaB_double, flag_deltaB_sign
    use util_mod, only: binary_search
    use sort_mod, only: sort

    implicit none

    ! interfaces
    
!     interface convert_ilut
!         module procedure convert_ilut_toGUGA
!         module procedure convert_ilut_toNECI
!     end interface convert_ilut
! 
    interface isProperCSF_ilut
        module procedure isProperCSF_b
        module procedure isProperCSF_sys
    end interface isProperCSF_ilut

    interface find_switches
        module procedure find_switches_ilut
        module procedure find_switches_stepvector
    end interface find_switches
!     interface calcB_vector
!         module procedure calcB_vector_nI
!         module procedure calcB_vector_ilut
!     end interface calcB_vector

!     interface isProperCSF
!         ! used to check if it is a proper csf (b>=0) etc. 
!         module procedure isProperCSF_nI
!         !module procedure isProperCSF_ilut
!     end interface isProperCSF

contains
    subroutine getExcitation_guga(nI, nJ, ex) 
        ! routine to determine excitation in guga basis
        ! for now to it very naively and unellegant by converting to ilut 
        ! format and calculating occupation vectors
        integer, intent(in) :: nI(nEl), nJ(nEl) 
        integer, intent(out) :: ex(2,2)
        character(*), parameter :: this_routine = "getExcitation_guga"

        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot) 
        integer :: first, last, cnt_e, cnt_h, occ_diff(nSpatOrbs), i

        call EncodeBitDet_guga(nI, ilutI)
        call EncodeBitDet_guga(nJ, ilutJ) 

        occ_diff = calcOcc_vector_int(ilutI(0:nifd)) - calcOcc_vector_int(ilutJ(0:nifd))

        select case(sum(abs(occ_diff)))

            case (0)
                ! there is a different definition of RDMs in the guga case or? 
                ! because for _RL_ -> ^RL^ excitations or nI = nJ i end up here
                ! but a lot of orbital index combinations can lead to 
                ! this type of excitation.. where should i assign it to..? 
                ! read the R.Shephard paper and think of a new calculation of 
                ! the RDM calculation in the GUGA case.. 

                ! ... for now, but in the indices of the first and last switch..

                first = findFirstSwitch(ilutI, ilutJ, 1, nSpatOrbs)

                last = findLastSwitch(ilutI, ilutJ, first, nSpatOrbs)

                ex(1,1) = 2*first
                ex(1,2) = 2*last-1

                ex(2,1) = 2*first-1
                ex(2,2) = 2*last

            case (2) 

                ! this is a "normal" double excitation 
                ! find the electron in nI which gets excited

                ex(1,2) = 0
                ex(2,2) = 0

                do i = 1, nSpatOrbs
                    if (occ_diff(i) == 1) ex(1,1) = 2 * i
                    if (occ_diff(i) == -1) ex(2,1) = 2 * i
                end do

            case (4)

                ! keep count of the already found electrons and holes 
                cnt_e = 1
                cnt_h = 1

                do i = 1, nSpatOrbs

                    select case(occ_diff(i))
                        
                        ! this choice of default spin-orbitals below 
                        ! makes certain two_rdm samplings as default alos..
                        ! not sure if this choice alone is valid..
                        ! also have to ensure i get the "spins" right so the 
                        ! rest of the NECI RDM routines can handle that..
                        case (2)
                            ! two eletrons get excited from orb i: 
                            ex(1,1) = 2 * i - 1
                            ex(1,2) = 2 * i

                        case (1) 
                            ! one electron gets excited from i 
                            ! at the first encountered electron cnt_e = 1
                            ! -> so this below gives me an alpha electron! 
                            ! -> at the second it will be an beta to ensure 
                            ! i get "correct" spins..
                            ex(1,cnt_e) = 2 * i - cnt_e + 1

                            cnt_e = cnt_e + 1 

                        case (-1) 
                            ! one hole found 
                            ex(2,cnt_h) = 2 * i - cnt_h + 1

                            cnt_h = cnt_h + 1

                        case (-2)
                            ! two electron get excited to orb i
                            ex(2,1) = 2 * i - 1
                            ex(2,2) = 2 * i 

                    end select 

                end do
                
        end select

    end subroutine getExcitation_guga
        
    subroutine find_switches_ilut(ilut, ind, lower, upper)
        ! for single excitations this checks for available switches around an 
        ! already chosen index
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: ind
        integer, intent(out) :: lower, upper
        character(*), parameter :: this_routine = "find_switches_ilut"

        integer :: i
        ! set defaults if no such switches are available
        lower = 1
        upper = nSpatOrbs

        if (isOne(ilut,ind)) then
            do i = ind - 1, 2, -1
                if (isTwo(ilut,i)) then
                    lower = i
                    exit
                end if
            end do
            do i = ind + 1, nSpatOrbs - 1
                if (isTwo(ilut,i)) then
                    upper = i
                    exit
                end if
            end do
        else if (isTwo(ilut,ind)) then
            do i = ind - 1, 2, -1
                if (isOne(ilut,i)) then
                    lower = i
                    exit
                end if
            end do
            do i = ind + 1, nSpatOrbs - 1
                if (isOne(ilut,i)) then
                    upper = i
                    exit
                end if
            end do
        end if

    end subroutine find_switches_ilut

    subroutine find_switches_stepvector(ind, lower, upper) 
        ! same as above but using the already calculated stepvector 
        integer, intent(in) :: ind
        integer, intent(out) :: lower, upper
        character(*), parameter :: this_routine = "find_switches_stepvector"

        integer :: switch, i
        ! set defaults
        lower = 1
        upper = nSpatOrbs

        if (current_stepvector(ind) == 1) then
            switch = 2

        else if (current_stepvector(ind) == 2) then
            switch = 1

        else
            ! wrong input! 
            call stop_all(this_routine, "wrong input! stepvalue /= {1,2}!")
        end if

        do i = ind - 1, 2, -1
            if (current_stepvector(i) == switch) then
                lower = i
                exit
            end if
        end do
        do i = ind + 1, nSpatOrbs - 1
            if (current_stepvector(i) == switch) then
                upper = i
                exit
            end if
        end do

    end subroutine find_switches_stepvector

    function findFirstSwitch(iI, iJ, start, semi) result(orb)
        ! write a scratch implementation to find the first change in 
        ! stepvector for two given CSFs. do it inefficiently for now
        ! improve later on
        integer(n_int), intent(in) :: iI(0:nifguga), iJ(0:nifguga)
        integer, intent(in) :: start, semi
        integer :: orb, a, b
        character(*), parameter :: this_routine = "findFirstSwitch"

        integer :: i
        ! with the fortran 2008 intrinsic funcitons it would be easy...
        ! for now just do a loop over double overlap region and compare
        ! stepvalues

        ! i could also use the current_stepvector quantity here.. or? 
        ! if it is always called for the current looked at ilut..
        ! i guess it does..

        orb = 0
        do i = start, semi - 1
!             if (getStepvalue(iI,i) /= getStepvalue(iJ,i)) then
            a = getStepvalue(iI,i)
!             a = current_stepvector(i)
            b = getStepvalue(iJ,i)
            if (a /= b) then
                orb = i
                return
            end if
        end do

    end function findFirstSwitch

    function findLastSwitch(ilutI, ilutJ, semi, ende) result(orb)
        ! function to find last switch in a mixed fullstop excitation
        integer(n_int), intent(in) :: ilutI(0:nifguga), ilutJ(0:nifguga)
        integer, intent(in) :: ende, semi
        integer :: orb, a, b
        character(*), parameter :: this_routine = "findLastSwitch"

        integer :: iOrb

        ! set it to impossible value, so contribution does not get 
        ! calculated if no switch happened, (which shouldnt be reached anyway)
        orb = nSpatOrbs + 1

        do iOrb = ende, semi + 1, -1
            a = getStepvalue(ilutI,iOrb)
!             a = current_stepvector(iOrb)
            b = getStepvalue(ilutJ,iOrb)
            if (a /= b) then
!             if (getStepvalue(ilutI,iOrb) /= getStepvalue(ilutJ,iOrb)) then
                orb = iOrb
                return
            end if
        end do
    
    end function findLastSwitch

    ! write custom add_ilut_lists
    subroutine add_guga_lists(nDets1, nDets2, list1, list2)
        use DetBitOps, only: ilut_lt, ilut_gt
        use sort_mod, only: sort
        use util_mod, only: binary_search_custom
        integer, intent(inout) :: nDets1
        integer, intent(in) :: nDets2 
        integer(n_int), intent(inout) :: list1(0:,1:), list2(0:,1:)
        character(*), parameter :: this_routine = "add_guga_lists"

        integer :: i, min_ind, pos, abs_pos

        ! first sort lists to use binary search 
        call sort(list1(:,1:nDets1), ilut_lt, ilut_gt)
        call sort(list2(:,1:nDets2), ilut_lt, ilut_gt)

        abs_pos = 0
        min_ind = 1

        do i = 1, nDets2
           
            pos = binary_search(list1(:,min_ind:ndets1), list2(:,i), nifd + 1)

            if (pos > 0) then
                ! try new implementation of that without the need of an extra
                ! output list
                
                ! need the absolute position after binary searches in only 
                ! sublists
                abs_pos = abs_pos + pos

                ! when element found just update the matrix element and update
                ! the indices 
                call encode_matrix_element(list1(:,abs_pos), &
                    extract_matrix_element(list1(:,abs_pos),1) + &
                    extract_matrix_element(list2(:,i),1), 1)

                ! min_ind to search next element is then 
                min_ind = min_ind + pos


!                 listOut(:,nDetsOut+1:nDetsOut + pos -1) = list1(:,min_ind:min_ind+pos-2)
! 
!                 nDetsOut = nDetsOut + pos
! 
!                 weight = extract_matrix_element(list1(:,min_ind+pos-1), 1)
! 
!                 weight = weight + extract_matrix_element(list2(:,i), 1)
! 
!                 listOut(:, nDetsOut) = list2(:,i)
! 
!                 call encode_matrix_element(listOut(:,nDetsOut), weight, 1)
! 
!                 min_ind = min_ind + pos

            else
                ! if the entry is not in list1 i have to move down all 
                ! subsequent elements after the found position and insert the 
                ! new entry at the indicated absolute position
                abs_pos = abs_pos - pos

                list1(:,abs_pos + 1:nDets1+1) = list1(:,abs_pos:nDets1)
                ! and add the new entry
                list1(:,abs_pos) = list2(:,i)

                ! the minimum index to search from now on should not include 
                ! the inserted element, since there should be no repetitions 
                ! in both lists
                min_ind = min_ind - pos

                ! and i have to update the number of determinants in list1
                nDets1 = nDets1 + 1

!                 listOut(:,nDetsOut+1:nDetsOut-pos-1) = list1(:,min_ind:min_ind-pos-2)
! 
!                 nDetsOut = nDetsOut - pos
! 
!                 listOut(:,nDetsOut) = list2(:,i)
! 
!                 min_ind = min_ind - pos - 1

            end if
        end do
                
    end subroutine add_guga_lists

    ! have to write custom x0 and x1 matrix encoding functions for the 
    ! GUGA excitation integer lists, as i need an additional entry for the 
    ! x1 matrix element


    subroutine write_guga_list(nunit, ilut)
        integer(n_int), intent(in) :: ilut(:,:)
        integer, intent(in) :: nunit

        integer :: i

        print *, " ilut list: "
        print *, " ==========="
        do i = 1, size(ilut,2)
            call write_det_guga(nunit, ilut(:,i))
        end do
        print *, " ==========="

    end subroutine write_guga_list

    subroutine write_det_guga(nunit, ilut, flag)
        ! subroutine which prints the stepvector representation of an ilut
        integer(n_int), intent(in) :: ilut(0:nIfGUGA)
        integer, intent(in) :: nunit
        logical, intent(in), optional :: flag

        integer :: step(nSpatOrbs), i

        step = calcStepvector(ilut(0:nifd))

        write(nunit,'("(")', advance='no')

        do i = 1, nSpatOrbs
            write(nunit, '(i3)', advance = 'no') step(i)
            if (i /= nSpatOrbs) write(nunit, '(",")', advance = 'no')
        end do
        write(nunit,'(")")', advance='no')
        

        write(nunit,'("(")', advance='no')
        do i = 1, 2
            write(nunit, "(f16.7)", advance = 'no') extract_matrix_element(ilut,i)
            if (i /= 2) write(nunit, "(A)", advance = 'no') ","
        end do

        if (present(flag)) then
            if (flag) then
                write(nunit, "(A,i3)", advance = 'yes') ") ", ilut(nifguga)
            else
                write(nunit, "(A,i3)", advance = 'no') ") ", ilut(nifguga)
            end if
        else
            write(nunit, "(A,i3)", advance = 'yes') ") ", ilut(nifguga)
        end if

    end subroutine write_det_guga

    subroutine encode_matrix_element(ilut, mat_ele, mat_type)
        ! encodes the x0 or x1 matrix element needed during the excitation
        ! creation. 
        ! mat_ele   ... x0 or x1 matrix element
        ! mat_type  ... 1...x0, 2...x1 
        integer(n_int), intent(inout) :: ilut(0:nIfGUGA)
        real(dp), intent(in) :: mat_ele
        integer, intent(in) :: mat_type
        character(*), parameter :: this_routine = "encode_matrix_element"

        integer(n_int) :: mat_int ! integer version of real

        ASSERT(mat_type == 1 .or. mat_type == 2)

        mat_int = transfer(mat_ele, mat_int)

        ilut(nIfD + mat_type) = mat_int

    end subroutine encode_matrix_element

    function extract_matrix_element(ilut, mat_type) result(mat_ele)
        ! function to extract matrix element of a GUGA ilut
        integer(n_int), intent(in) :: ilut(0:nIfGUGA)
        integer, intent(in) :: mat_type
        real(dp) :: mat_ele
        character(*), parameter :: this_routine = "extract_matrix_element"

        ASSERT(mat_type == 1 .or. mat_type == 2)

        mat_ele = transfer(ilut(nIfD + mat_type), mat_ele)

    end function extract_matrix_element

    subroutine update_matrix_element(ilut, mat_ele, mat_type)
        ! function to update already encoded matrix element multiplicative
        integer(n_int), intent(inout) :: ilut(0:nIfGUGA)
        real(dp), intent(in) :: mat_ele
        integer, intent(in) :: mat_type
        character(*), parameter :: this_routine = "update_matrix_element"

        integer(n_int) :: mat_int
        real(dp) :: temp_ele

        ASSERT(mat_type == 1 .or. mat_type == 2)

        temp_ele = transfer(ilut(nIfD + mat_type), temp_ele)

        mat_int = transfer(temp_ele * mat_ele, mat_int)

        ilut(nIfD + mat_type) = mat_int

    end subroutine update_matrix_element

    function count_beta_orbs_ij(ilut, i, j) result(nOpen)
        ! function to count the number of 1s in a CSF det between spatial 
        ! orbitals i and j
        integer(n_int), intent(in) :: ilut(0:nifd)
        integer, intent(in) :: i, j
        integer :: nOpen
        character(*), parameter :: this_routine = "count_beta_orbs_ij"

        integer(n_int) :: mask(0:nifd), beta(0:nifd), alpha(0:nifd)
        integer :: k

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        nOpen = 0

        ! quick and dirty fix to deal with the excitation range mask probs:
        ! do i always call that for the current det the excitation is 
        ! calculated for?  i think so..
        do k = i, j
!             if (isOne(ilut,k)) then 
            if (current_stepvector(k) == 1) then
                nOpen = nOpen + 1
            end if
        end do
! 
!         if (i < j) then
!             mask = getExcitationRangeMask(i, j)
!             mask = iand(ilut, mask)
!             beta = iand(mask, MaskBeta)
!             alpha = iand(mask, MaskAlpha)
! 
!             alpha = ishft(alpha,-1_n_int)
! 
!             beta = iand(ieor(alpha,beta),beta)
! 
!             nOpen = CountBits(beta, nifd)
!         else if (i == j) then
!             if (isOne(ilut,i)) then
!                 nOpen = 1
!             end if
!         end if

    end function count_beta_orbs_ij

    function count_alpha_orbs_ij(ilut, i, j) result(nOpen)
        ! function to count the number of 2s in a CSF det between spatial 
        ! orbitals i and j
        integer(n_int), intent(in) :: ilut(0:nifd)
        integer, intent(in) :: i, j
        integer :: nOpen
        character(*), parameter :: this_routine = "count_alpha_orbs_ij"

        integer(n_int) :: mask(0:nifd), alpha(0:nifd), beta(0:nifd)
        integer :: k

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        nOpen = 0

        ! quick fix for now to see if thats the problem: loop and check! 
        do k = i, j
!             if (isTwo(ilut,k)) then
            if (current_stepvector(k) == 2) then
                nOpen = nOpen + 1
            end if 
        end do
! 
!         if (i < j) then
!             mask = getExcitationRangeMask(i, j)
!             if ( i == 2 .and. j == 16) then
!                 print *, "mask bef: ", mask
!             end if
!             mask = iand(ilut, mask)
!             alpha = iand(mask, MaskAlpha)
!             beta = iand(mask, MaskBeta)
! 
!             beta = ishft(beta,+1_n_int)
! 
!             alpha = iand(ieor(beta,alpha),alpha)
! 
!             nOpen = CountBits(alpha, nifd)
!         else if (i == j) then
!             if (isTwo(ilut,i)) then
!                 nOpen = 1
!             end if
!         end if
! ! 
!         if (i == 2 .and. j == 16) then 
!             print *, "nOpen: ", nOpen
!             print *, "ilut: ", ilut
!             print *, "mask:" , mask
!             print *, "alpha: ", alpha
!             print *, "beta: ", beta
!         end if

    end function count_alpha_orbs_ij

    function count_open_orbs_ij(i, j, L) result(nOpen)
        ! function to calculate the number of open orbitals between spatial
        ! orbitals i and j in ilut. i and j have to be given ordered i<j
        integer, intent(in) :: i, j
        integer(n_int), intent(in), optional :: L(0:nifd)
        integer :: nOpen
        character(*), parameter :: this_routine = "count_open_orbs_ij"
        
        logical :: flag
        integer(n_int) :: mask
        integer :: k

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)
!         ASSERT(i < j)
        ! scrap this assert and change in that way to output 0 if the indices
        ! dont fit or are reversed. to deal with to short overlap ranges

        nOpen = 0

        ! also here a quick fix do deal with excitrangemask probs:

        ! if the ilut input is present use it otherwise just look at the 
        ! current_stepvector 
        if (present(L)) then
            do k = i, j
                flag = isOne(L,k)
                if (flag .or. isTwo(L,k)) then 
                    nOpen = nOpen + 1
                end if
            end do
        else
            do k = i, j 
                if (current_stepvector(k) == 1 .or. current_stepvector(k) == 2) then
                    nOpen = nOpen + 1
                end if
            end do
        end if

!         if (i < j) then
! 
!             ! first have to create a integer mask were every bit between 
!             ! correspongin spin orbitals for i and j is set. and then i can give
!             ! it to the already provided open orbital counting function
!             mask = getExcitationRangeMask(i, j)
! 
!             ! use it to indicate only these orbitals
!             nOpen = count_open_orbs(iand(L, mask))
!         else if (i == j) then
!             ! have to do this to avoid too long lines...
!             flag = isOne(L,i)
!             if (flag .or. isTwo(L,i)) then
!                 nOpen = 1
!             end if
!         end if

    end function count_open_orbs_ij

    function getExcitationRangeMask(i, j) result(mask)
        ! function to create an integer corresponding to a mask were every 
        ! bits between the spin orbitals corresponding to spatial orbitals 
        ! i and j are set to 1 and 0 else. used to access the excitation 
        ! range for GUGA excitation calculations
        integer, intent(in) :: i, j
        integer(n_int) :: mask
        character(*), parameter :: this_routine = "getExcitationRangeMask"

        integer(n_int) :: k
        integer(n_int) :: tmp_i, tmp_j

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)
        ASSERT(j > i)

        tmp_i = int(i,n_int)
        tmp_j = int(j,n_int)

        ! not quite sure about LMSB or RMSB... todo
        mask = 0_n_int

        do k = 2_n_int*tmp_i - 1_n_int, 2_n_int*tmp_j ! convert to spin orbitals
            mask = mask + 2_n_int**(k-1_n_int)
        end do

    end function getExcitationRangeMask

    subroutine setDeltaB(deltaB, ilut)
        ! routine to encode the deltaB value in a given CSF in ilut bit
        ! representation, by using the newly defined flags:
        ! flag_deltaB_sign   ... 7
        ! flag_deltaB_single ... 5 
        ! and if necessary
        ! flag_deltaB_double ... 6
        integer, intent(in) :: deltaB
        integer(n_int), intent(inout) :: ilut(0:nIfGUGA)
        character(*), parameter :: this_routine = "set_delta_b"
    

        ASSERT(deltaB <= 2 .and. deltaB >= -2)

        ! should no just be:
        ilut(nIfGUGA) = deltaB

!         select case(deltaB)
!             case(-2)
!                 call set_flag(ilut, flag_deltaB_double, .true.)
!                 call set_flag(ilut, flag_deltaB_sign, .true.)
!                 call set_flag(ilut, flag_deltaB_single, .false.)
! 
!             case (-1)
!                 call set_flag(ilut, flag_deltaB_double, .false.)
!                 call set_flag(ilut, flag_deltaB_sign, .true.)
!                 call set_flag(ilut, flag_deltaB_single, .true.)
! 
!             case (0)
!                 call set_flag(ilut, flag_deltaB_double, .false.)
!                 call set_flag(ilut, flag_deltaB_sign, .false.)
!                 call set_flag(ilut, flag_deltaB_single, .false.)
! 
!             case (1)
!                 call set_flag(ilut, flag_deltaB_double, .false.)
!                 call set_flag(ilut, flag_deltaB_sign, .false.)
!                 call set_flag(ilut, flag_deltaB_single, .true.)
! 
!             case (2)
!                 call set_flag(ilut, flag_deltaB_double, .true.)
!                 call set_flag(ilut, flag_deltaB_sign, .false.)
!                 call set_flag(ilut, flag_deltaB_single, .false.)
! 
!             case default
!                 call stop_all(this_routine, "invalid deltaB value!")
! 
!         end select
    end subroutine setDeltaB

    function getDeltaB(ilut) result(deltaB)
        ! function to get the deltaB value encoded in the flag-byte in ilut
        integer(n_int), intent(in) :: ilut(0:nIfGUGA)
        integer :: deltaB
        character(*), parameter :: this_routine = "getDeltaB"

        ! check if flags are correctly set
!         ASSERT(.not.(test_flag(ilut, flag_deltaB_double) .and. test_flag(ilut, flag_deltaB_single)))

        
        ! and this should now jsut be:
        deltaB = ilut(nIfGUGA)

!         if (test_flag(ilut, flag_deltaB_double)) then
!             if (test_flag(ilut, flag_deltaB_sign)) then
!                 deltaB = -2
!             else
!                 deltaB = 2
!             end if
! 
!         else 
!             if (test_flag(ilut, flag_deltaB_single)) then
!                 if (test_flag(ilut, flag_deltaB_sign)) then
!                     deltaB = -1
!                 else
!                     deltaB = 1
!                 end if
! 
!             else 
!                 deltaB = 0
!             end if
!         end if

    end function getDeltaB

    subroutine convert_ilut_toNECI(ilutG, ilutN, HElement)
        integer(n_int), intent(in) :: ilutG(0:nifguga)
        integer(n_int), intent(inout) :: ilutN(0:niftot)
        HElement_t(dp), intent(out), optional :: HElement
        character(*), parameter :: this_routine = "convert_ilut_toNECI"

        ASSERT(isProperCSF_ilut(ilutG))

        ! i think i just need to copy over the det part again
        ilutN(0:nifd) = ilutG(0:nifd)

        ! and then extract the matrix element 
        ! here i need to check what type of matrix element is necessary 
        ! dependent on which type of compilation, 
        ! extract_matrix_element always gives a real(dp)!
        if (present(HElement)) then
#ifdef __CMPLX
            HElement = complex(extract_matrix_element(ilutG,1), 0.0_dp)
#else
            HElement = extract_matrix_element(ilutG, 1)
#endif
        end if

    end subroutine convert_ilut_toNECI

    subroutine convert_ilut_toGUGA(ilutN, ilutG)
        integer(n_int), intent(in) :: ilutN(0:niftot)
        integer(n_int), intent(out) :: ilutG(0:nifguga)
        character(*), parameter :: this_routine = "convert_ilut_toGUGA"

        ! need only the det part essentially..
        ilutG(0:nifd) = ilutN(0:nifd)

        ! and set matrix elements to 1 and delta b to 0
        call encode_matrix_element(ilutG,1.0_dp,1)
        call encode_matrix_element(ilutG,1.0_dp,2)

        call setDeltaB(0,ilutG)

    end subroutine convert_ilut_toGUGA


    function isProperCSF_nI(nI) result(flag)
        ! function to check if provided CSF in nI(nEl) format is a proper CSF
        integer, intent(in) :: nI(nEl)
        logical :: flag

        real(dp) :: bVector(nEl)
        integer :: calcS

        flag = .true.

        bVector = calcB_vector_nI(nI)

        ! check if the b value drops below zero
        if (any(bVector < 0.0_dp)) flag = .false.

        ! check if total spin is same as input
        calcS = abs(sum(get_spin_pn(nI(:))))

        if (calcS /= STOT) flag = .false.

    end function isProperCSF_nI

    function isProperCSF_b(ilut) result(flag)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        logical :: flag

        flag = .true.

        ! check if b value drops below zero
        if (any(calcB_vector_int(ilut(0:nifd)) < 0)) flag = .false.

    end function isProperCSF_b

    function isProperCSF_sys(ilut, sysFlag) result(flag)
        ! function to check if provided CSF in ilut format is a proper CSF
        ! checks b vector positivity and is total S is correct
        integer(n_int), intent(in) :: ilut(0:nifguga)
        logical, intent(in):: sysFlag
        logical :: flag

        integer(n_int) :: tmp_ilut(0:niftot)

        flag = .true.

        ! check if b value drops below zero
        if (any(calcB_vector_int(ilut(0:nifd)) < 0)) flag = .false.

        ! check if total spin is same as input, cant avoid loop here i think..
        !calcS = 0
        !do iOrb = 1, nBasis/2
        !    if (isOne(ilut,iOrb)) calcS = calcS + 1
        !    if (isTwo(ilut,iOrb)) calcS = calcS - 1
        !end do
        
        tmp_ilut = 0_n_int
        tmp_ilut(0:nifd) = ilut(0:nifd)

        ! if system flag is also given as input also check if the CSF fits 
        ! concerning total S and the number of electrons
        if (sysFlag) then
            if (abs(return_ms(tmp_ilut)) /= STOT) then
                print *, "CSF does not have correct total spin!:"
                call write_det_guga(6,ilut)
                print *, "System S: ", STOT
                print *, "CSF S: ", abs(return_ms(tmp_ilut))
                flag = .false.
            end if

            if (int(sum(calcOcc_vector_ilut(ilut(0:nifd)))) /= nEl) then
                print *, "CSF does not have right number of electrons!:"
                call write_det_guga(6,ilut)
                print *, "System electrons: ", nEl
                print *, "CSF electrons: ", int(sum(calcOcc_vector_ilut(ilut(0:nifd))))
                flag = .false.
            end if
        end if

    end function isProperCSF_sys

    function calcB_vector_nI(nI) result(bVector)
        ! function to calculate the bVector from a CSF given in nI 
        ! representation. Gives b vector also of length nEl. 
        ! not yet quite sure if i should output b as integers or real
        integer, intent(in) :: nI(nEl)
        real(dp) :: bVector(nEl), bValue

        integer :: iOrb, i, inc

        ! init
        iOrb = 1
        bVector = 0.0_dp
        bValue = 0.0_dp
        inc = 0

        do while (iOrb <= nEl)
            if (isDouble(nI, iOrb)) then
                bVector(iOrb) = bValue
                bVector(iOrb + 1) = bValue
                inc = 2
    
            else
                inc = 1
                if is_alpha(nI(iOrb)) then
                    bValue = bValue - 1.0_dp

                else 
                    bValue = bValue + 1.0_dp
                
                end if
                bVector(iOrb) = bValue

            end if
            
            iOrb = iOrb + inc

        end do

    end function calcB_vector_nI


    function isDouble(nI, iOrb) result(flag)
        ! returns a logical .true. if spinorbital iOrb is doubly occupied 
        ! and .false. elsewise.
        integer, intent(in) :: nI(nEl), iOrb
        logical :: flag

        integer :: pair
        
        flag = .false.
        pair = ab_pair(nI(iOrb))

        ! ask simon about the ordering in nI. If its not always ordered I 
        ! have to to a quicksearch for pair in nI. If its ordered I can 
        ! just check adjacent entries in nI

        ! assume ordered for now! 
        if (iOrb == 1) then
            ! special conditions if iOrb is 1
            flag = (nI(2) == pair)

        else if (iOrb == nEl) then
            flag = (nI(nEl-1) == pair)

        else
            flag = ( (nI(iOrb-1) == pair) .or. (nI(iOrb+1) == pair))

        end if

    end function isDouble


    function calcStepvector(ilut) result(stepVector)
        ! function to calculate stepvector of length nReps, corresponding
        ! to the ilut bit-representation, if each stepvalue is needed 
        ! often within a function.
        ! there is probably a very efficient way of programming that! 
        !TODO ask simon for improvements.
        integer(n_int), intent(in) :: ilut(0:nIfD)
        integer :: stepVector(nSpatOrbs)
        integer :: iOrb

        do iOrb = 1, nSpatOrbs
            stepVector(iOrb) = getStepvalue(ilut, iOrb)
        end do

    end function calcStepvector

    function calcOcc_vector_ilut(ilut) result(occVector)
        ! probably more efficiently implemented by simon already... 
        ! but for now do it in this stupid way todo -> ask simon
        integer(n_int), intent(in) :: ilut(0:nIfD)
        real(dp) :: occVector(nSpatOrbs)

        integer :: iOrb

        do iOrb = 1, nSpatOrbs
            occVector(iOrb) = getSpatialOccupation(ilut,iOrb)
        end do

    end function calcOcc_vector_ilut

    function calcOcc_vector_int(ilut) result(occVector)
        ! function which gives the occupation vector in integer form
        integer(n_int), intent(in) :: ilut(0:nifd)
        integer :: occVector(nSpatOrbs)

        integer :: i

        do i = 1, nSpatOrbs

            occVector(i) = getSpatialOccupation(ilut,i)
! 
!             if (isZero(ilut,i)) then
!                 occVector(i) = 0
!             else if (isThree(ilut,i)) then
!                 occVector(i) = 2
!             else
!                 occVector(i) = 1
!             end if
        end do

    end function calcOcc_vector_int

    function calcB_vector_ilut(ilut) result(bVector)
        ! function to calculate bVector of length (nBasis) for a given
        ! CSF bitRep ilut of length (2*nBasis), save this b-vector 
        ! in the nI array, normally used to store occupied orbitals
        ! update: changed convention to not use nI for bVector but also for 
        ! occupied orbitals as in normal NECI implementation
        ! but nethertheless need a b-vector of lenght of spatial orbitals 
        ! for excitaiton and matrix element calculation
        ! UPDATE: change b vector calculation to update b vector directly on
        ! the spot of the corresponding stepvector, eg:
        ! d = 0 1 2 3
        ! b = 0 1 0 0
        ! because otherwise i actually only needed the bvector of one orbital
        ! ahead. and never the first entry otherwise. and that would cause 
        ! problems when accessing b vector at the end of an excitaiton if 
        ! that is the last orbital..
        integer(n_int), intent(in) :: ilut(0:nIfD)
        real(dp) :: bVector(nSpatOrbs)        ! b-vector stored in bVector
        integer :: i
        real(dp) :: bValue 

        ! loop over CSF entries and increment, decrement bValue 
        ! accordingly, have to correctly access ilut entries and 
        ! check if they are 0,1,2 or 3... -> how to for multiple 
        ! integers?
        bVector = 0.0_dp
        bValue = 0.0_dp

        do i = 1, nSpatOrbs

            if (isOne(ilut,i)) then
                bValue = bValue + 1.0_dp
            else if (isTwo(ilut,i)) then
                bValue = bValue - 1.0_dp
            end if
            ! define bvalue to always only get updated for the next
            ! UPDATE: changed definition to update on the spot.
            bVector(i) = bValue
        end do

    end function calcB_vector_ilut

    function calcB_vector_int(ilut) result(bVector)
        ! function to calculate the bvector in integer form
        integer(n_int), intent(in) :: ilut(0:nifd) 
        integer :: bVector(nSpatOrbs)

        integer :: i, bValue

        bVector = 0
        bValue = 0
        
        do i = 1, nSpatOrbs
            
            if (isOne(ilut, i)) then
                bValue = bValue + 1
                
            else if (isTwo(ilut, i)) then
                bValue = bValue - 1

            end if

            bVector(i) = bValue

        end do

    end function calcB_vector_int


    function getStepvalueExp(ilut, sOrb) result(stepValue)
        ! function to get stepvector value of a given spatial orbital
        ! sOrb -> has to later be included in "macros.h" for efficiency
        
        integer(n_int), intent(in) :: ilut(0:nIfD)
        integer, intent(in) :: sOrb     ! spatial orbital
        integer :: stepValue            ! resulting stepvector value
        ! if the number of orbitals is too large, multiple integers 
        ! in bit-representation are used to encode a single determinant
        ! have to figure out this access to the iluts
        ! determinants are stores as a list of integer in form of iluts
        ! first have to figure out which integer to take: 
        integer :: indInt, offset, mask

        ! integer division to get the necessary ilut entry, remember
        ! the occupation of the spatial orbitals are asked for
        indInt = int((sOrb-1)/(bits_n_int/2))


        ! now i have to pick out the spatial alpha-beta combination, 
        ! corresponding to the asked for sOrb, with a mask
        ! the offset, where this mask has to be shifted to:
        offset = int(2*mod((sOrb-1), bits_n_int/2))   
        
        ! shift (11) = 3 mask left to corresponding spinorbitals 
        mask = ishft(3,offset)

        ! now pick out corresponding spin orbitals with iand() and 
        ! shift back to the first postion to give integer
        stepValue = ishft(iand(ilut(indInt),mask),-offset)

        
        ! already have some function in macros.h

    end function getStepvalueExp

    pure subroutine EncodeBitDet_guga(nI, ilut) 
        ! special function to encode bit dets for the use in the guga 
        ! excitation generation
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(out) :: ilut(0:nifguga)

        integer :: i, pos

        ilut = 0_n_int

        do i = 1, nEl
            pos = (nI(i) - 1) / bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(nI(i)-1, bits_n_int))
        end do

    end subroutine EncodeBitDet_guga
            
    function getSpatialOccupation(iLut, s) result(nOcc)

        integer(n_int), intent(in) :: ilut(0:nIfD)
        integer, intent(in) :: s
        real(dp) :: nOcc

        if (isZero(ilut, s)) then
            nOcc = 0.0_dp

        else if (isThree(ilut, s)) then
            nOcc = 2.0_dp

        else
            nOcc = 1.0_dp
        end if
        
    end function getSpatialOccupation
! 
! 
!     function calcMeanB(nI) result(meanB)
!         ! function to calculate mean value of the b-vector, used in 
!         ! the branching probabilities for future switches
!         integer, intent(in) :: nI(nReps)
!         real(dp) :: meanB
! 
!         meanB = sum(nI)/(max(1,nReps))
! 
!     end function calcMeanB
! 
end module guga_bitRepOps
#endif
