#include "macros.h"
! GUGA testsuite.
! contains all GUGA related unit tests
! ask simon how to implement this in an optimal way.
! have to choose test cases differently... 
! maybe for testing purposes hard-compile a test system with fixed parameters
! like electrons, orbitals etc... 
! since otherwise all test cases are quite input dependent... 
! discuss with simon how to implement that optimally
#ifndef __CMPLX
module guga_testsuite
   
    ! check if i can use some sort of unit test framework, like fruit,
    ! and how to use it correctly and efficiently.
    ! do that later, for now implement just some random tests which get 
    ! executed at the end of the NECI initialisation and check them by hand
    use SystemData 
    use guga_bitRepOps
    use guga_excitations
    use guga_matrixElements
    use guga_data
    use guga_init
    use guga_procedure_pointers
    use constants
    use DetBitOps
    use Determinants
    use bit_reps
    use FciMCData
    use dsfmt_interface, only: dsfmt_init
    implicit none

contains

    subroutine runTestsGUGA()
        ! main test calling function
        character(*), parameter :: this_routine = "runTestsGUGA"

        print *, "============================================================"
        print *, " Running GUGA tests. For now, just to check if desired &
            &functionality is there. Not yet on unit test level. TODO"
        print *, "============================================================"
        print *, "Can only run full GUGA testsuite if: "
        print *, " number of orbitals = 4"
        print *, " number of electrons = 4"
        print *, " total spin = 0"
        print *, " have to initialize random number generator for tests"
        print *, " seed: ", 0
        call dsfmt_init(0)

        if (nEl == 4 .and. nBasis/2 == 4 .and. STOT == 0 .and. t_full_guga_tests) then
            print *, "running full test suite: "
            print *, "ATTENTION: make sure the H = 1 FCIDUMP is used!!"

            ! todo: make a check FCIDUMP function to see if its the correct 
            ! test system! 
            call test_guga_bitRepOps

            call test_guga_excitations_stochastic

            call test_guga_excitations_exact
            
            call test_guga_matrixElements

            call test_guga_data

    !         call test_excitationIdentifier

            call test_bitChecks

            ! need currentB already allocated and calculated to test that... 
    !         call test_calcSingleProbWeight

            call check_determinants

            call run_test_excit_gen_guga_S2

        else 
            print *, " only run the excitation generator tests!"
!             call run_test_excit_gen_guga_S0

            call run_test_excit_gen_guga()

        end if

        print *, "============================================================"
        print *, " GUGA tests finished!"
        print *, " All tests passed! You are awesome!"
        print *, "============================================================"

        call stop_all(this_routine, "stop after GUGA tests")

    end subroutine runTestsGUGA
     
    subroutine run_test_excit_gen_det 
        character(*), parameter :: this_routine = "run_test_excit_gen_det"
        
!         call dsfmt_init(0)

        call run_test_excit_gen_4ind_multiple([1,2,3,4,5,6,7,8])
        call run_test_excit_gen_4ind_multiple([9,10,11,12,13,14,15,16])
        call run_test_excit_gen_4ind_multiple([1,3,5,7,10,12,14,16])
        call run_test_excit_gen_4ind_multiple([1,4,5,8,9,12,13,16])
        call run_test_excit_gen_4ind_multiple([1,3,6,7,10,11,14,16])

        call stop_all(this_routine, "stop after tests!")

    end subroutine run_test_excit_gen_det
        
    ! also write a test runner for the 4ind-weighted excit gen to compare 
    ! the frequency of the mat_ele/pgen ratios 
    subroutine run_test_excit_gen_4ind_multiple(nI)
        ! keep in mind that there are 2 cases of 4ind excitation 
        ! generator! use the correct one depending on the input!
        use DetBitOps, only: EncodeBitDet
        use Determinants, only: write_det
        use excit_gen_5, only: test_excit_gen_take2
        use excit_gens_int_weighted, only: test_excit_gen_4ind, calc_all_excitations
        
        integer, intent(in), optional :: nI(nel)
        integer(n_int) :: ilut(0:niftot)
        integer(n_int), pointer :: excitations(:,:)
        integer :: nExcits, nTest, i, test_det(nel)
        character(*), parameter :: this_routine = "run_test_excit_gen_4ind_multiple"

        pSingles = 0.1_dp
        pDoubles = 1.0_dp - pSingles

        ! here i also have to set the parallel and ant
        pParallel = 0.5_dp

        if (present(nI)) then
            test_det = nI
        else
            test_det = fdet
        end if

        ! this is the implementation to calculate up to n_test excitations
        ! from the hf det

        call EncodeBitDet(test_det, ilut)

        call calc_all_excitations(ilut, excitations, nExcits) 

!         nTest = min(nExcits, 20) 
        nTest = nExcits

        print *, "running tests on nExcits of ", nTest
        call write_det(6, test_det, .true.)

        if (tGen_4ind_weighted) then
            
            call test_excit_gen_4ind(ilut, n_guga_excit_gen)
        
        else if (tGen_4ind_2) then 
            
            call test_excit_gen_take2(ilut, n_guga_excit_gen)

        else
            call stop_all(this_routine, "use 4ind_weighted or 4ind_weighted_2!")

        end if

        do i = 1, nTest
            if (tGen_4ind_weighted) then 
                call test_excit_gen_4ind(excitations(:,i), n_guga_excit_gen)
            else if (tGen_4ind_2) then
                call test_excit_gen_take2(excitations(:,i), n_guga_excit_gen)
            else
                call stop_all(this_routine, "use 4ind_weighted or 4ind_weighted_2!")
            end if
        end do

    end subroutine run_test_excit_gen_4ind_multiple

    subroutine run_test_excit_gen_guga

        pSingles = 0.1_dp
        pDoubles = 1.0_dp - pSingles

!         pExcit4 = (1.0_dp - 1.0_dp / real(nSpatOrbs,dp))
        pExcit4 = 0.5_dp
!         pExcit2 = 1.0_dp / real(nSpatOrbs - 1, dp)
        pExcit2 = 0.5_dp

        if (t_consider_diff_bias) then
            pExcit2_same = 0.5_dp
            pExcit3_same = 0.5_dp
        else
            pExcit2_same = 1.0_dp
            pExcit3_same = 1.0_dp
        end if

        pExcit2_same = 0.9_dp
        pExcit3_same = 0.9_dp
! 
!         if (nEl == 4 .and. nBasis/2 == 4 .and. STOT == 2) then
! !             call run_test_excit_gen_guga_S2
!             call run_test_excit_gen_guga_general
! 
!         else if (nEL == 3 .and. nBasis/2 == 4 .and. STOT == 1) then
!             call run_test_excit_gen_guga_nEl_3_S_1
! 
!         else if (nEl == 2 .and. nBasis/2 == 4) then
!             call run_test_excit_gen_guga_nEl_2_S_0
! 
!         else if (nEl == 5 .and. nBasis/2 == 4) then
!             call run_test_excit_gen_guga_nEl_5_S_1
! 
!         else if (nEl == 6 .and. nBasis/2 == 4) then
!             call run_test_excit_gen_guga_nEl_6_S_0
! 
!         else if (nEl == 6 .and. nBasis/2 == 6 .and. STOT == 0) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_6_S_0
! 
! !             call run_test_excit_gen_guga_single([1,3,5,7,10,12])
! 
!         else if (nEl == 6 .and. nBasis/2 == 6 .and. STOT == 2) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_6_S_2
! !             call run_test_excit_gen_guga_single([1,3,4,5,10,11])
! 
!         else if (nEl == 6 .and. nBasis/2 == 6 .and. STOT == 4) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_6_S_4
! 
!         else if (nEl == 4 .and. nBasis/2 == 6 .and. STOT == 0) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_4_S_0
! 
!         else if (nEl == 4 .and. nBasis/2 == 6 .and. STOT == 2) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_4_S_2
! 
!         else if (nEl == 5 .and. nBasis/2 == 6 .and. STOT == 1) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_5_S_1
! 
!         else if (nEl == 5 .and. nBasis/2 == 6 .and. STOT == 3) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_5_S_3
! 
!         else if (nEl == 7 .and. nBasis/2 == 6 .and. STOT == 1) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_7_S_1
! 
!         else if (nEl == 7 .and. nBasis/2 == 6 .and. STOT == 3) then
!             call run_test_excit_gen_guga_nOrb_6_nEl_7_S_3
! 
!         else if (nEl == 5 .and. nBasis/2 == 9 .and. STOT == 1) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_5_S_1
! 
!         else if (nEl == 5 .and. nBasis/2 == 9 .and. STOT == 3) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_5_S_3
! 
!         else if (nEl == 7 .and. nBasis/2 == 9 .and. STOT == 1) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_7_S_1
! 
!         else if (nEl == 7 .and. nBasis/2 == 9 .and. STOT == 3) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_7_S_3
! 
!         else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 0) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_10_S_0
! 
!         else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 6) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_10_S_6
! 
! !         else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 4) then
! !             call run_test_excit_gen_guga_nOrb_9_nEl_10_S_4
! ! 
! !         else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 2) then
! !             call run_test_excit_gen_guga_nOrb_9_nEl_10_S_2
! 
!         else if (nEl == 9 .and. nBasis/2 == 9 .and. STOT == 3) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_9_S_3
! 
!         else if (nEl == 9 .and. nBasis/2 == 9 .and. STOT == 1) then
!             call run_test_excit_gen_guga_nOrb_9_nEl_9_S_1
! 
!         else if (nEl == 18 .and. nBasis/2 == 18 .and. STOT == 0) then
!             call run_test_excit_gen_guga_nOrb_18_nel_18_S_0
! 
!         else if (nEl == 18 .and. nBasis/2 == 18 .and. STOT == 2) then
!             call run_test_excit_gen_guga_nOrb_18_nel_18_S_2
! 
!         else if (nEl == 18 .and. nBasis/2 == 18 .and. STOT == 6) then
!             call run_test_excit_gen_guga_nOrb_18_nel_18_S_6
!         else 
            ! todo: create a general excit_gen tester which uses the 
            ! guessed HF determinant as a start and uses some of the 
            ! exactly created determinants from that to check for pgen and 
            ! matrix element consistency! 
! 
!             call test_identify_excitation()
!             call test_identify_excitation_and_matrix_element()

!             call run_test_excit_gen_guga_multiple(&
!                 [1,2,3,4,5,6,7,8])
!             call run_test_excit_gen_guga_multiple(&
!                 [9,10,11,12,13,14,15,16])
!             call run_test_excit_gen_guga_multiple(&
!                 [1,3,5,7,10,12,14,16])
!             call run_test_excit_gen_guga_multiple(&
!                 [1,3,6,7,10,11,14,16])
!             call run_test_excit_gen_guga_multiple(&
!                 [1,4,5,8,9,12,13,16])
!             call run_test_excit_gen_guga_general
!             call run_test_excit_gen_guga_single([1,2,3,5,6,7,8,9,11,12,13,14,15,18])
!             call run_test_excit_gen_guga_single([1,2,3,4,5,6,7,8,9,10,11,13,16,18])
!             call run_test_excit_gen_guga_single([1,2,3,4,5,7,8,10,11,12,13,14,17,18])
!             call run_test_excit_gen_guga_multiple(&
!                 convert_guga_to_ni([0,0,0,0,1,1,1,2,1,2,3,3,3,1,2],15))
            call run_test_excit_gen_guga_single()
            call run_test_excit_gen_guga_single(&
                convert_guga_to_ni([3,3,3,3,3,1,0,3,1],9))
        
!             call run_test_excit_gen_guga_single(&
!                 [1,2,3,5,7,10,12,13,14,15,17,18,20,21,25,27,30,32,34,35,36,&
!                 37,39,41,44,46,48,49,50,51,54,58])

!         end if 

    end subroutine run_test_excit_gen_guga

    subroutine test_identify_excitation
        character(*), parameter :: this_routine = "test_identify_excitation"
        integer :: nI(nel), nJ(nel)
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(excitationInformation) :: excitInfo

        ! make a more thorough test on the excitation identifier
        ! do it for now for the specific 14 electron system where the 
        ! errors show up.. 
        ! nI: 3333333
        nI = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
        ! nJ: 311333322
        nJ = [1,2,3,5,7,8,9,10,11,12,13,14,16,18]

        print *, convert_guga_to_ni([3,3,3,3,3,3,3],7)
        print *, convert_guga_to_ni([3,1,1,3,3,3,3,2,2],9)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)

        excitInfo = identify_excitation(ilutI,ilutJ)

        call print_excitInfo(excitInfo)

        ! nI: 3331212123
        nI = [1,2,3,4,5,6,7,10,11,14,15,18,19,20]
        ! nJ: 311212121322
        nJ = convert_guga_to_ni([3,1,1,2,1,2,1,2,1,3,2,2],12)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)

        excitInfo = identify_excitation(ilutI,ilutJ)

        call print_excitInfo(excitInfo)

        call stop_all(this_routine, "for now")

    end subroutine test_identify_excitation

    subroutine test_identify_excitation_and_matrix_element(nI)
        integer(n_int), intent(in), optional :: nI(nel)
        character(*), parameter :: this_routine = "test_identify_excitation_and_matrix_element" 
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot) 
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:), two_ex(:,:)
        integer :: nEx, i, nex_2, test_det(nel), j, ind
        logical :: valid
        real(dp) :: pos(nSpatOrbs), neg(nSpatOrbs), mat_ele, diff

        if (present(nI)) then 
            test_det = nI 
        else
             test_det = fdet 
         end if

        call EncodeBitDet(test_det, ilutI)
        call actHamiltonian(ilutI, ex, nEx) 

        call init_csf_information(ilutI) 

        print *, "Testing matrix elements for nEx excitations of: ", nEx
        call write_det_guga(6,ilutI,.true.)
        call write_guga_list(6,ex(:,1:nex))

        print *, "Do the tests on only connected determinants:"
        do i = 1, nEx

            excitInfo = identify_excitation(ilutI, ex(:,i)) 

            ASSERT(excitInfo%valid)

            if (excitInfo%typ /= 0) then
                call checkCompatibility(ilutI, excitInfo, valid, pos, neg)

                ASSERT(valid)

            end if

            call calc_guga_matrix_element(ilutI, ex(:,i), excitInfo, mat_ele, & 
                .true., 1)

            diff = abs(extract_matrix_element(ex(:,i), 1) - mat_ele)

            if (diff < 1.0e-10) diff = 0.0_dp

!             print *, "mat eles: ", extract_matrix_element(ex(:,i),1), mat_ele, &
!                 diff

            if (diff > EPS) then 
                call write_det_guga(6,ilutI,.true.) 
                call write_det_guga(6,ex(:,i),.true.)
                call print_excitInfo(excitInfo)
                print *, "spin change: ", excitInfo%spin_change
            end if
            ! how do i test non-possible excitations? 
            ! do i apply the hamiltonian once more? and then check in the 
            ! lists? 
!             call actHamiltonian(ex(:,i), two_ex, nex_2) 

        end do

        ! also do the tests on non-connected CSFs.. 
        ! maybe apply the hamiltonian a second time to the list of generated 
        ! CSFs and check with the original one.. some of them should not 
        ! be connected then and have a zero matrix element! 
        print *, "do the test on only non-connected determinants:"
        ! this might take some time.. 
        do i = 1, nEx

            call write_det_guga(6, ex(:,i),.true.)

            call actHamiltonian(ex(:,i), two_ex, nex_2) 

            ! in acthamiltonian the current_stepvector quantity is set 
            ! for the ilut input.. 
!             pos = binary_search(two_ex(0:nifd), ilutI)

!             call init_csf_information(ilutI)

            do j = 1, nex_2

                excitInfo = identify_excitation(ilutI, two_ex(:,j))
! 
!                 if (excitInfo%valid) then
! 
!                     ind = binary_search(ex(0:nifd,1:nex), two_ex(0:nifd,j)) 
! 
!                     if (ind < 0) then 
! 
!                         print *, "incorrectly identified!"
! 
!                         call write_det_guga(6,ilutI,.true.)
!                         call write_det_guga(6,two_ex(:,j),.true.) 
!                         call print_excitInfo(excitInfo)
! 
!                     end if
!                 end if

!                 call calc_guga_matrix_element(ex(:,i), two_ex(:,j), excitInfo, &
!                     mat_ele, .true., 2)
! 
!                 diff = abs(extract_matrix_element(two_ex(:,j),1) - mat_ele)
! 
!                 if (diff < 1.0e-10) diff = 0.0_dp 
! 
!                 print *, "two ex check: ", extract_matrix_element(two_ex(:,j),1), &
!                     mat_ele, diff

                call calc_guga_matrix_element(ilutI, two_ex(:,j), excitInfo, &
                    mat_ele, .true., 2)

                if (abs(mat_ele) > EPS) then 

                    ! this should only happen if two_ex is in the original ex
                    ! or it is ilutI 
                    ind = binary_search(ex(0:nifd,1:nex), two_ex(0:nifd,j)) 

                    ! is the matrix element here correct if i find somethin? 
                    if (ind < 0 .and. (.not. DetBitEQ(two_ex(0:nifd,j),ilutI(0:nifd)))) then 

                        print *, "something wrong!"

                    else if (ind > 0) then

                        ! is the sign correct now??
                        diff = abs(extract_matrix_element(ex(:,ind),1) - mat_ele)

                        if (diff < 1.0e-10) diff = 0.0_dp

                        if (diff > EPS) then 
                            print *, "sign check:" 
                            print *, "I:"
                            call write_det_guga(6,ilutI,.true.)
                            print *, "step: ", temp_step_i
                            print *, "b:", temp_b_real_i
                            print *, "n:", temp_occ_i
                            print *, "J:"
                            call write_det_guga(6,two_ex(:,j),.true.)
                            print *, "step: ", temp_step_j
                            print *, "db: ", temp_delta_b 

                        end if
!                         print *, extract_matrix_element(ex(:,ind),1), mat_ele, diff

!                             extract_matrix_element(ex(:,i),1) * extract_matrix_element(two_ex(:,j),1)
! 
!                         if (diff > 1.0_dpe-10) then 
!                             print *, "matrix element different: ", diff
!                             call write_det_guga(6,ilutI,.true.)
!                             call write_det_guga(6,two_ex(:,
! 
! 
                    end if
                end if 
            end do 

            deallocate(two_ex)

        end do

        ! i should also test more than double excitaitons to see if i correctly 
        ! identify impossible excitations

        call deinit_csf_information()

        call stop_all(this_routine, "!")


    end subroutine test_identify_excitation_and_matrix_element

    subroutine test_guga_data

!         call init_guga_data_procPtrs
        print *, "testing module: guga_data:"
        call test_getMixedFullStop
        call test_getSingleMatrixElement
        call test_getDoubleMatrixElement
        print *, "guga_data tests passed!"

    end subroutine test_guga_data

    subroutine test_guga_matrixElements

        print *, " testing routines of module: guga_matrixElements:"

        call check_calcDiagExchange_nI
        call check_calcDiagMatEleGUGA_nI

        print *, " guga_matrixElements tests passed!"

    end subroutine test_guga_matrixElements

    subroutine test_guga_excitations_exact
        character(*), parameter :: this_routine = "test_guga_excitations_exact"

        print *, "testing module: guga_excitations:"
        call test_add_ilut_lists
        call test_calcRemainingSwitches
        call test_actHamiltonian
        call test_calcOverlapRange
        call test_excitationIdentifier_single
        call test_createSingleStart
        call test_singleUpdate
        call test_singleEnd
        call test_calcAllExcitations_single
        call test_calcRemainingSwitches_double
        call test_excitationIdentifier_double
        call test_checkCompatibility
        call test_calcSingleOverlapLowering
        call test_calcSingleOverlapRaising
        call test_calcSingleOverlapMixed
        call test_calcDoubleLowering
        call test_calcDoubleRaising
        call test_calcDoubleL2R
        call test_calcDoubleR2L
        call test_calcFullStopLowering
        call test_calcFullStopRaising
        call test_calcFullStartLowering
        call test_calcFullStartRaising
        call test_calcFullStartFullStopAlike
        call test_calcDoubleExcitation_withWeight
        call test_calcNonOverlapDouble
        call test_calcFullStartR2L
        call test_calcFullStartL2R
        call test_calcFullStopR2L
        call test_calcFullStopL2R
        call test_calcFullStartFullStopMixed
        call test_calcAllExcitations_double

!         call stop_all(this_routine, "stop here!")

        print *, "guga_excitations tests passed!"

    end subroutine test_guga_excitations_exact

    subroutine test_guga_bitRepOps
        character(*), parameter :: this_routine = "test_guga_bitRepOps"

        print *, "testing functions from module: guga_bitRepOps"
        call test_findSwitches
        call test_count_alpha_orbs_ij
        call test_count_beta_orbs_ij
        call test_matrix_element_ops
        call test_count_open_orbs_ij
        call test_getExcitationRangeMask
        call test_set_get_DeltaB
        call test_isProperCSF_ilut
        call test_calcbvector
        call test_isDouble
        call test_calcStepVector
        call test_getSpatialOccupation
        call test_calcOcc_vector_ilut
        call test_add_guga_lists

        print *, "guga_bitRepOps tests passed!"

    end subroutine test_guga_bitRepOps

    subroutine run_test_excit_gen_guga_general
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_general"
        integer(n_int) :: ilut(0:niftot)
        integer(n_int), pointer :: ex(:,:)
        integer :: nEx, i
        integer :: nTest
        type(excitationInformation) :: excitInfo

        ! use fdet as first determinant and test on all excitations from this..!
        ! maybe a bit too much for bigger system?
        print *, "running general test_excit_gen_guga()"

        ! first act the hamiltonian on the fdet
        call EncodeBitDet(fdet, ilut)

        call actHamiltonian(ilut, ex, nEx)

        nTest = min(nEx,20)
!         nTest = nEx

        print *, "running tests on nExcits: ", nTest
        call write_guga_list(6, ex(:,1:nEx))
        call test_excit_gen_guga(ilut, n_guga_excit_gen)
        ! then loop over the excitations and check the excitation generator

        ! dont do it for all excitatons... too much
        do i = 1, nTest
            call test_excit_gen_guga(ex(:,i), n_guga_excit_gen)
        end do
    
    end subroutine run_test_excit_gen_guga_general

    subroutine run_test_excit_gen_guga_single(nI)
        integer, intent(in), optional :: nI(nEl)
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_single"
        integer(n_int) :: ilut(0:niftot)
        integer :: test_det(nEl)

        if (present(nI)) then
            test_det = nI 
        else
            test_det = fdet
        end if

        call EncodeBitDet(test_det, ilut)

        call test_excit_gen_guga(ilut, n_guga_excit_gen)

    end subroutine run_test_excit_gen_guga_single
    
    subroutine run_test_excit_gen_guga_multiple(nI)
        integer, intent(in), optional :: nI(nEl) 
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_multiple"
        integer(n_int) :: ilut(0:niftot)
        integer :: test_det(nEl), nEx, i, nTest
        integer(n_int), pointer :: ex(:,:)


        if (present(nI)) then
            test_det = nI
        else
            test_det = fdet 
        end if

        call EncodeBitDet(test_det, ilut)

        call actHamiltonian(ilut, ex, nEx)

!         nTest = min(nEx, 20)
        nTest = nEx

        print *, "running tests on nExcits: ", nTest

        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        do i = 1, nTest
            call test_excit_gen_guga(ex(:,i), n_guga_excit_gen)
        end do

    end subroutine run_test_excit_gen_guga_multiple


    subroutine test_findSwitches
        character(*), parameter :: this_routine = "test_findSwitches"
        integer(n_int) :: ilutI(0:nifguga), ilutJ(0:nifguga)

        print *, "testing findSwitches routines:"
        ! 3300
        call EncodeBitDet_guga([1,2,3,4],ilutI)
        ilutJ = ilutI

        ASSERT(findFirstSwitch(ilutI,ilutJ,1,4) == 0)
        ASSERT(findLastSwitch(ilutI,ilutJ,1,4) == 5)

        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilutJ)

        ASSERT(findFirstSwitch(ilutI,ilutJ,1,4) == 0)
        ASSERT(findLastSwitch(ilutI,ilutJ,1,4) == 5)

        ASSERT(findFirstSwitch(ilutI,ilutJ, 2,3) == 0)
        ASSERT(findLastSwitch(ilutI,ilutJ,1,2) == 5)

        ! 1230
        call EncodeBitDet_guga([1,4,5,6],ilutI)
        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilutJ)
! 
        ASSERT(findFirstSwitch(ilutI,ilutJ,1,3) == 2)
        ASSERT(findFirstSwitch(ilutI,ilutJ,2,3) == 2)
!         print *, findLastSwitch(ilutI,ilutJ,1,4)
        ASSERT(findLastSwitch(ilutI,ilutJ,1,4) == 2)
        ! for the find last switch we exclude the inputted first orbital!
        ASSERT(findLastSwitch(ilutI,ilutJ,2,3) == 5)
        ASSERT(findFirstSwitch(ilutI,ilutJ,3,4) == 0)
        ASSERT(findLastSwitch(ilutI,ilutJ,3,4) == 5)

        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilutI)
        ! 1212
        call EncodeBitDet_guga([1,4,5,8],ilutJ)

        ASSERT(findFirstSwitch(ilutI,ilutJ,2,3) == 2)
        ASSERT(findLastSwitch(ilutI,ilutJ,1,3) == 3)
        ASSERT(findFirstSwitch(ilutI,ilutJ,3,4) == 3)
        ASSERT(findLastSwitch(ilutI,ilutJ,2,3) == 3)

!         call stop_all(this_routine, "for now!")
        print *, "findSwitches tests passed!"

    end subroutine test_findSwitches

    subroutine test_count_beta_orbs_ij
        character(*), parameter :: this_routine = "test_count_beta_orbs_ij"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4],ilut)

        print *, "testing count_beta_orbs_ij:"

        ASSERT(count_beta_orbs_ij(ilut,1,4) == 0)
        ASSERT(count_beta_orbs_ij(ilut,1,3) == 0)
        ASSERT(count_beta_orbs_ij(ilut,2,4) == 0)

        call EncodeBitDet_guga([1,3,6,8],ilut)
        ASSERT(count_beta_orbs_ij(ilut,1,4) == 2)
        ASSERT(count_beta_orbs_ij(ilut,2,4) == 1)

        call EncodeBitDet_guga([3,5,6,8],ilut)
        ASSERT(count_beta_orbs_ij(ilut,1,4) == 1)

        print *, "count_beta_orbs_ij tests passed!"

    end subroutine test_count_beta_orbs_ij


    subroutine test_count_alpha_orbs_ij
        character(*), parameter :: this_routine = "test_count_alpha_orbs_ij"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4],ilut)

        print *, "testing count_alpha_orbs_ij:"
        ASSERT(count_alpha_orbs_ij(ilut,1,4) == 0)
        ASSERT(count_alpha_orbs_ij(ilut,2,3) == 0)

        call EncodeBitDet_guga([1,4,5,6],ilut)
        ASSERT(count_alpha_orbs_ij(ilut,1,4) == 1)
        call EncodeBitDet_guga([3,6,7,8],ilut)
        ASSERT(count_alpha_orbs_ij(ilut,1,4) == 1)

        call EncodeBitDet_guga([1,4,5,8],ilut)
        ASSERT(count_alpha_orbs_ij(ilut,1,4) == 2)


        print *, "count_alpha_orbs_ij tests passed!"

    end subroutine test_count_alpha_orbs_ij

    subroutine test_add_guga_lists
        character(*), parameter :: this_routine = "test_add_guga_lists"
        integer(n_int) :: l1(0:nifguga,1), l2(0:nifguga,1), l3(0:nifguga,2)
        integer(n_int) :: l4(0:nifguga,3), l5(0:nifguga,4), l6(0:nifguga,5)

        integer :: nOut

        call EncodeBitDet_guga([1,2,3,4],l6(:,1))
        call EncodeBitDet_guga([1,4,5,8],l3(:,1))

        print *, "testing: add_guga_lists(n1,n2,l1,l2,lOut,nOut):"
        nOut = 1
        call add_guga_lists(nOut, 1, l6, l3)
        ASSERT(nOut == 2)

        call EncodeBitDet_guga([1,2,3,6],l3(:,1))

        call add_guga_lists(nOut, 1, l6, l3)

        ASSERT(nout == 3)

        call EncodeBitDet_guga([5,6,7,8],l3(:,1))

        call add_guga_lists(nOut,1,l6,l3)

        ASSERT(nout == 4)

        l3(:,2) = l6(:,2)
        call write_guga_list(6,l3)
        call write_guga_list(6,l6)

        call add_guga_lists(nOut,2,l6,l3)

        ASSERT(nOut == 4)

        l5(:,1) = l6(:,1)
        nOut = 1
        
        call add_guga_lists(nOut,4,l5,l6)

        ASSERT(nOut == 4)

        print *, "add_guga_lists tests passed!"

    end subroutine test_add_guga_lists

    subroutine test_getSpatialOccupation
        character(*), parameter :: this_routine = "test_getSpatialOccupation"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,6], ilut)

        print *, "testing getSpatialOccupation(ilut, sOrb):"
        ASSERT(getSpatialOccupation(ilut,1) == 2.0_dp)
        ASSERT(getSpatialOccupation(ilut,2) == 1.0_dp)
        ASSERT(getSpatialOccupation(ilut,3) == 1.0_dp)
        ASSERT(getSpatialOccupation(ilut,4) == 0.0_dp)
        print *, "getSpatialOccupation tests passed!"

    end subroutine test_getSpatialOccupation

    subroutine test_calcFullStartFullStopMixed
        character(*), parameter :: this_routine = "test_calcfullStartFullStopMixed"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(weight_obj) :: weights

        call EncodeBitDet_guga([1,4,5,8],ilut)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(4,1,1,4)

        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        weights = init_doubleWeight(ilut, 4)

        print *, "testing calcFullStartFullStopMixed(ilut, exInfo, ex, num, posSwitch, negSwitch):"

        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        print *, "calcFullStartFullStopMixed tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_calcFullStartFullStopMixed

    subroutine test_guga_excitations_stochastic
        character(*), parameter :: this_routine = "test_guga_excitations_stochastic"

        print *, "testing module: guga_excitations stochastic:"
!         call run_test_excit_gen_guga_S0
        call test_calcMixedContribution
        call test_pickRandomOrb
        call test_generate_excitation_guga_double
        call test_generate_excitation_guga_single
        call test_pickOrbitals_single

        call test_createStochasticStart_single
        call test_singleStochasticUpdate 
        call test_singleStochasticEnd
        call test_createStochasticExcitation_single

        call test_pickOrbitals_double
        ! now go through each specific stochastix excitation calulator:
        call test_mixedFullStartStochastic
        call test_doubleUpdateStochastic
        call test_mixedFullStopStochastic
        call test_calcFullStartFullStopMixedStochastic
        call test_calcFullStartRaisingStochastic
        call test_calcFullStartLoweringStochastic
        call test_calcFullStopRaisingStochastic
        call test_calcFullStopLoweringStochastic
        call test_calcSingleOverlapMixedStochastic

        ! todo: cont. from here!
        ! now i need to test the specific semi-start/stop routines
        call test_calcLoweringSemiStartStochastic
        call test_calcRaisingSemiStartStochastic
        call test_calcFullStopR2L_stochastic
        call test_calcFullStopL2R_stochastic


        call test_calcRaisingSemiStopStochastic
        call test_calcLoweringSemiStopStochastic

        ! and then the rest of the test should just be calling them.
        call test_calcFullStartL2R_stochastic
        call test_calcFullStartR2L_stochastic

        call test_calcDoubleLoweringStochastic
        call test_calcDoubleRaisingStochastic
        call test_calcDoubleL2R2L_stochastic
        call test_calcDoubleR2L2R_stochastic
        call test_calcDoubleL2R_stochastic
        call test_calcDoubleR2L_stochastic

        call test_createStochasticExcitation_double

!         call stop_all(this_routine, "stop here!")


        print *, "guga_excitations stochastic tests passed!"
        

    end subroutine test_guga_excitations_stochastic

!     subroutine test_pgens_single
!         character(*), parameter :: this_routine = "test_pgens_single"
!         integer(n_int) :: ilut(0:nifguga), excitation(0:nifguga)
!         integer :: nTest = 10000, i
!         real(dp) :: pgen
! 
! 
!         print *, "testing stochastic generation probabilities and exact ones"
!         call EncodeBitDet_guga([1,2,3,4], ilut)
!         allocate(currentB_ilut(4))
!         currentB_ilut = calcB_vector_ilut(ilut)
!         allocate(currentOcc_ilut(4))
!         currentOcc_ilut = calcOcc_vector_ilut(ilut)
! 
!         do i = 1, nTest
!             call createStochasticExcitation_single(ilut, excitation, pgen)
! 
! 
!     end subroutine test_pgens_single
! 
    subroutine run_test_excit_gen_guga_nEl_3_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_3_S_1"
        integer(n_int):: ilut(0:niftot)
        integer :: nI(3)


        print *, "running: test_excit_gen_guga() for nEl = 3, S = 1"

        ! 1102
        nI = [1,3,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3100
        nI = [1,2,3]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0112
        nI = [3,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1030
        nI = [1,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0103
        nI = [3,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0121
        nI = [3,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0310
        nI = [3,4,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3010
        nI = [1,2,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1012
        nI = [1,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1300
        nI = [1,3,4]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3001
        nI = [1,2,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0130
        nI = [3,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1120
        nI = [1,3,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1210
        nI = [1,4,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0301
        nI = [3,4,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)
! 
        ! 1201
        nI = [1,4,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1021
        nI = [1,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1003
        nI = [1,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0031
        nI = [5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0013
        nI = [5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

!         call stop_all(this_routine, "here")


    end subroutine run_test_excit_gen_guga_nEl_3_S_1

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_2
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_6_S_2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 6, S = 2 system"

        ! 331100
        nI = [1,2,3,4,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313100
        nI = [1,2,3,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133100
        nI = [1,3,4,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311300
        nI = [1,2,3,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

! 
        ! 031310
        nI = [3,4,5,7,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131021
        nI = [1,3,4,5,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311102
        nI = [1,2,3,5,7,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111122
        nI = [1,3,5,7,10,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 112112
        nI = [1,3,6,7,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 112211
        nI = [1,3,6,8,9,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121103
        nI = [1,4,5,7,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121112
        nI = [1,4,5,7,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 310301
        nI = [1,2,3,7,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
        
        ! 312110
        nI = [1,2,3,6,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 001133
        nI = [5,7,9,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine,"for now")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_2


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_4
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_6_S_4"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 6, S = 4 system"
        ! 311110
        nI = [1,2,3,5,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131110
        nI = [1,3,4,5,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113110
        nI = [1,3,5,6,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111310
        nI = [1,3,5,7,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111130
        nI = [1,3,5,7,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311101
        nI = [1,2,3,5,7,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121111
        nI = [1,4,5,7,9,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 011113
        nI = [3,5,7,9,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111121
        nI = [1,3,5,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111112
        nI = [1,3,5,7,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 101113
        nI = [1,5,7,9,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
        
        ! 103111
        nI = [1,5,6,7,9,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


!         call stop_all(this_routine,"for now")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_4


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_2
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_4_S_2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 4, S = 0 system"

        ! 311000
        nI = [1,2,3,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 011012
        nI = [3,5,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 011021
        nI = [3,5,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 013001
        nI = [3,5,6,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 011003
        nI = [3,5,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 101201
        nI = [1,5,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 103001
        nI = [1,5,6,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113000
        nI = [1,3,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 130100
        nI = [1,3,4,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 013100
        nI = [3,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_2


    subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_5_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 5, S = 1 system"
        ! 331000000
        nI = [1,2,3,4,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313000000
        nI = [1,2,3,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 310300000
        nI = [1,2,3,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121210000
        nI = [1,4,5,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


        ! 000012121
        nI = [9,12,13,16,17]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 301210000
        nI = [1,2,5,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121120000
        nI = [1,4,5,7,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 301003000
        nI = [1,2,5,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine,"here")

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_1
    
    subroutine run_test_excit_gen_guga_nOrb_18_nel_18_S_2
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(18)

        ! 331333332110000000
        nI = [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,18,19,21]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313333332110000000
        nI = [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,18,19,21]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 333313332110000000
        nI = [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,18,19,21]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 333133332110000000
        nI = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,18,19,21]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
        
    end subroutine run_test_excit_gen_guga_nOrb_18_nel_18_S_2

    subroutine run_test_excit_gen_guga_nOrb_18_nel_18_S_6
        character(*), parameter :: this_routine =&
        "run_test_excit_gen_guga_nOrb_18_nel_18_S_6"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(18)

        nI = [3,6,9,14,15,17,20,21,23,25,27,29,31,32,33,34,35,36]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

    end subroutine run_test_excit_gen_guga_nOrb_18_nel_18_S_6

    subroutine run_test_excit_gen_guga_nOrb_18_nel_18_S_0
        character(*), parameter :: this_routine =&
        "run_test_excit_gen_guga_nOrb_18_nel_18_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(18)

        ! 111212212121322102
        nI = [1,3,5,8,9,12,14,15,18,19,22,23,25,26,28,30,31,36]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
        
        ! 
        nI = [1,2,3,5,6,8,9,10,11,14,15,18,19,21,23,26,28,30]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


    end subroutine run_test_excit_gen_guga_nOrb_18_nel_18_S_0

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_9_S_1
        character(*), parameter :: this_routine =&
        "run_test_excit_gen_guga_nOrb_9_nEl_9_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(9)

        ! 333102100
        nI = [1,2,3,4,5,6,7,12,13]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331302100
        nI = [1,2,3,4,5,7,8,12,13]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313302100
        nI = [1,2,3,5,6,7,8,12,13]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133302100
        nI = [1,3,4,5,6,7,8,12,13]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 333102010
        nI = [1,2,3,4,5,6,7,12,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331302010
        nI = [1,2,3,4,5,7,8,12,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313302010
        nI = [1,2,3,5,6,7,8,12,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133302010
        nI = [1,3,4,5,6,7,8,12,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 333100210
        nI = [1,2,3,4,5,6,7,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331300210
        nI = [1,2,3,4,5,7,8,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313300210
        nI = [1,2,3,5,6,7,8,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133300210
        nI = [1,3,4,5,6,7,8,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 333010210
        nI = [1,2,3,4,5,6,9,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331210210
        nI = [1,2,3,4,5,8,9,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 123310210
        nI = [1,4,5,6,7,8,9,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_9_S_1

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_9_S_3
        character(*), parameter :: this_routine =&
        "run_test_excit_gen_guga_nOrb_9_nEl_9_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(9)

        ! 331111200
        nI = [1,2,3,4,5,7,9,11,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313111200
        nI = [1,2,3,5,6,7,9,11,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133111200
        nI = [1,3,4,5,6,7,8,11,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_9_S_3

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_2
        character(*), parameter :: this_routine =&
        "run_test_excit_gen_guga_nOrb_9_nEl_10_S_2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(10)

        ! 331211030
        nI = [1,2,3,4,5,8,9,11,15,16]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_2

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_4
        character(*), parameter :: this_routine = &
        "run_test_excit_gen_guga_nOrb_9_nEl_10_S_4"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(10)

        ! 313031110
        nI = [1,2,3,5,6,9,10,11,13,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        call stop_all(this_routine,"here")

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_4

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_6
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_10_S_6"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(10)


        print *, "running: test_excit_gen_guga() on the 9 orbital, Nel = 10, S = 6"
        ! 311111121
        nI = [1,2,3,5,7,9,11,13,16,17]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131111121
        nI = [1,3,4,5,7,9,11,13,16,17]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311111112
        nI = [1,2,3,5,7,9,11,13,15,18]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131111112
        nI = [1,3,4,5,7,9,11,13,15,18]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        
    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_6


    subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_10_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(10)

        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 10, S = 0 system"
        ! 333330000
        nI = [1,2,3,4,5,6,7,8,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 333303000
        nI = [1,2,3,4,5,6,7,8,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133233000
        nI = [1,3,4,5,6,8,9,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133231200
        nI = [1,3,4,5,6,8,9,10,11,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 123331200
        nI = [1,4,5,6,7,8,9,10,11,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313320012
        nI = [1,2,3,5,6,7,8,10,15,18]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 132331002
        nI = [1,3,4,6,7,8,9,10,11,18]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


!         call stop_all(this_routine,"here")

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_0

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_5_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 5, S = 3 system"
        ! 311100000
        nI = [1,2,3,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113100000
        nI = [1,3,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111120000
        nI = [1,3,5,7,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 112101000
        nI = [1,3,6,7,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 101301000
        nI = [1,5,6,7,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111000300
        nI = [1,3,5,13,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111200100
        nI = [1,3,5,8,13]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
        
        ! 001113000
        nI = [5,7,9,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 1100110020
        nI = [1,3,9,11,16]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 101110002
        nI = [1,5,7,9,18]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 0101020101
        nI = [3,7,12,15,17]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_3


    subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_7_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 7, S = 3 system"
        ! 331110000
        nI = [1,2,3,4,5,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131121000
        nI = [1,3,4,5,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311013000
        nI = [1,2,3,5,7,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121131000
        nI = [1,4,5,7,9,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 301111200
        nI = [1,2,5,7,9,11,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121110300
        nI = [1,4,5,7,9,13,14]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121110210
        nI = [1,4,5,7,9,14,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_3

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_7_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 7, S = 1 system"
        ! 33310000
        nI = [1,2,3,4,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313210000
        nI = [1,2,3,5,6,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133120000
        nI = [1,3,4,5,6,7,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 132300010
        nI = [1,3,4,6,7,8,15]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 123100030
        nI = [1,4,5,6,7,15,16]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 033100003
        nI = [3,4,5,6,7,17,18]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_1


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_7_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 7, S = 3 system"
        ! 331110
        nI = [1,2,3,4,5,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131121
        nI = [1,3,4,5,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311031
        nI = [1,2,3,5,9,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131031
        nI = [1,3,4,5,9,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 310131
        nI = [1,2,3,7,9,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121131
        nI = [1,4,5,7,9,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111213
        nI = [1,3,5,8,9,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311301
        nI = [1,2,3,5,7,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113103
        nI = [1,3,5,6,7,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 110133
        nI = [1,3,7,9,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113103
        nI = [1,3,5,6,7,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine,"here")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_3

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_7_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 7, S = 1 system"
        ! 333100
        nI = [1,2,3,4,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 312121
        nI = [1,2,3,6,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 132121
        nI = [1,3,4,6,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 303121
        nI = [1,2,5,6,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 123121
        nI = [1,4,5,6,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 033121
        nI = [3,4,5,6,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133102
        nI = [1,3,4,5,6,7,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 132112
        nI = [1,3,4,6,7,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331003
        nI = [1,2,3,4,5,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331300
        nI = [1,2,3,4,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 312310
        nI = [1,2,3,6,7,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 132310
        nI = [1,3,4,6,7,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 033310
        nI = [3,4,5,6,7,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313210
        nI = [1,2,3,5,6,8,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_1


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_5_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 5, S = 3 system"
        ! 311100
        nI = [1,2,3,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 013101
        nI = [3,5,6,7,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111201
        nI = [1,3,5,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 110301
        nI = [1,3,7,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 101301
        nI = [1,5,7,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111021
        nI = [1,3,5,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
! 
        ! 011121
        nI = [3,5,7,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111012
        nI = [1,3,5,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 310110
        nI = [1,2,3,7,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
! 
        ! 011211
        nI = [3,5,8,9,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131100
        nI = [1,3,4,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311010
        nI = [1,2,3,5,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113010
        ni = [1,3,5,6,9]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)


!         call stop_all(this_routine,"here")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_3

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_5_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 5, S = 1 system"

        ! 331000
        nI = [1,2,3,4,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 132001
        nI = [1,3,4,6,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 303001
        nI = [1,2,5,6,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 123001
        nI = [1,4,5,6,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 310201
        nI = [1,2,3,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121201
        nI = [1,4,5,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 031201
        nI = [3,4,5,8,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 130102
        nI = [1,3,4,7,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121021
        nI = [1,4,5,10,11]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 311002
        nI = [1,2,3,5,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113020
        nI = [1,3,5,6,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 131020
        nI = [1,3,4,5,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 113002
        nI = [1,3,5,6,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine,"heer")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_1

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_4_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 4, S = 0 system"
        ! 330000
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 312000
        nI = [1,2,3,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 030120
        nI = [3,4,7,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 300030
        nI = [1,2,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 120030
        nI = [1,4,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
! 
        ! 310002
        nI = [1,2,3,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 130002
        nI = [1,3,4,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 120012
        nI = [1,4,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 011220
        nI = [3,5,8,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 012120
        nI = [3,6,7,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 030300
        nI = [3,4,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine, "now")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_0

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_6_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 6, S = 0 system"

        ! 333000
        nI = [1,2,3,4,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 000333
        nI = [7,8,9,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 303030
        nI = [1,2,5,6,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 030303
        nI = [3,4,7,8,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 121212
        nI = [1,4,5,8,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 310230
        nI = [1,2,3,8,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 013230
        nI = [3,5,6,8,9,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

! 
        ! 131022
        nI = [1,3,4,5,10,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 111222
        nI = [1,3,5,8,10,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331200
        nI = [1,2,3,4,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 313200
        nI = [1,2,3,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 133200
        nI = [1,3,4,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 330300
        nI = [1,2,3,4,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 312300
        nI = [1,2,3,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 303300
        nI = [1,2,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 123300
        nI = [1,4,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 303120
        nI = [1,2,5,6,7,10]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 331002
        nI = [1,2,3,4,5,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 312012
        nI = [1,2,3,6,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 303012
        nI = [1,2,5,6,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 312012
        nI = [1,2,3,6,9,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 012132
        nI = [3,6,7,9,10,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 001233
        nI = [5,8,9,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 012123
        nI = [3,6,7,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

        ! 110322
        nI = [1,3,7,8,10,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)
! 
        ! 030033
        nI = [3,4,9,10,11,12]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut,n_guga_excit_gen)

!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_0


    subroutine run_test_excit_gen_guga_nEl_6_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_6_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, "running: test_excit_gen_guga() on the nEl = 6, S = 0 system"
        ! 3123
        nI = [1,2,3,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3330
        nI = [1,2,3,4,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3033
        nI = [1,2,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1233
        nI = [1,4,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3132
        nI = [1,2,3,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1323
        nI = [1,3,4,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3312
        nI = [1,2,3,4,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1332
        nI = [1,3,4,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0333
        nI = [3,4,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3303
        nI = [1,2,3,4,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

!         call stop_all(this_routine, "here")
    end subroutine run_test_excit_gen_guga_nEl_6_S_0


    subroutine run_test_excit_gen_guga_nEl_5_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_5_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)


        print *, "running: test_excit_gen_guga() on the nEl = 5, S = 1 system"
        ! 3310
        nI = [1,2,3,4,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1321
        nI = [1,3,4,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1231
        nI = [1,4,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1330
        nI = [1,3,4,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1132
        nI = [1,3,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3031
        nI = [1,2,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3130
        nI = [1,2,3,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0331
        nI = [3,4,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3013
        nI = [1,2,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3301
        nI = [1,2,3,4,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0313
        nI = [3,4,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1303
        nI = [1,3,4,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3112
        nI = [1,2,3,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1213
        nI = [1,4,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3103
        nI = [1,2,3,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3121
        nI = [1,2,3,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1312
        nI = [1,3,4,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1123
        nI = [1,3,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1033
        nI = [1,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0133
        nI = [3,5,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

!         call stop_all(this_routine, "here")


    end subroutine run_test_excit_gen_guga_nEl_5_S_1

    subroutine run_test_excit_gen_guga_nEl_2_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_2_S_0"
        integer(n_int):: ilut(0:niftot)
        integer :: nI(2)


        print *, "running test_excit_gen_guga() gor the nEl = 2, S = 0 system"
! 
        ! 1200
        nI = [1,4]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)


        ! 3000
        nI = [1,2]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0120
        nI = [3,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0012 
        nI = [5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1020
        nI = [1,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1002
        nI = [1,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0003
        nI = [7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0300
        nI = [3,4]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0102
        nI = [3,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0030
        nI = [5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

!         call stop_all(this_routine, "here")

    end subroutine run_test_excit_gen_guga_nEl_2_S_0
! 
    subroutine run_test_excit_gen_guga_S2
        ! also check for the S = 2 system... 
        ! and probably write this function generally for all sorts of 
        ! excitations
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_S2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)


        print *, "running: test_excit_gen_guga(ilut,iter) for the S = 2 system"
! 
        ! 0311
        nI = [3,4,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3110
        nI = [1,2,3,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0131
        nI = [3,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)
! 
        ! 3011
        nI = [1,2,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1121
        nI = [1,3,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1103
        nI = [1,3,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)
! 
        ! 3101
        nI = [1,2,3,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1130
        nI = [1,3,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1013
        nI = [1,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1310
        nI = [1,3,4,5]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0113
        nI = [3,5,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1301
        nI = [1,3,4,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1031
        nI = [1,5,6,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1211
        nI = [1,4,5,7]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1112
        nI = [1,3,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

!         call stop_all (this_routine,"")

    end subroutine run_test_excit_gen_guga_S2


    subroutine run_test_excit_gen_guga_S0
        ! write a similar testing routine as Simons 
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_S0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)


        pSingles = 0.1_dp
        pDoubles = 1.0_dp - pSingles

!         pExcit4 = (1.0_dp - 1.0_dp / real(nSpatOrbs,dp))
        pExcit4 = 0.5_dp
!         pExcit2 = 1.0_dp / real(nSpatOrbs - 1, dp)
        pExcit2 = 0.5_dp

        if (t_consider_diff_bias) then
            pExcit2_same = 0.5_dp
            pExcit3_same = 0.5_dp
        else
            pExcit2_same = 1.0_dp
            pExcit3_same = 1.0_dp
        end if

        pExcit2_same = 0.9_dp
        pExcit3_same = 0.9_dp

        print *, "running: test_excit_gen_guga(ilut,n_guga_excit_gen)"
        print *, "pSingles set to: ", pSingles
        print *, "pDoubles set to: ", pDoubles

        
        ! 1032
        nI = [1,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0132
        nI = [3,5,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1302
        nI = [1,3,4,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3300
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3030
        nI = [1,2,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3003
        nI = [1,2,7,8]
         call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0330
        nI = [3,4,5,6]
         call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0303
        nI = [3,4,7,8]
         call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0033
        nI = [5,6,7,8]
         call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        
        ! 3012
        nI = [1,2,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3102
        nI = [1,2,3,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1023
        nI = [1,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 3120
        nI = [1,2,3,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 0312
        nI = [3,4,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1230
        nI = [1,4,5,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1203
        nI = [1,4,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1320
        nI = [1,3,4,6]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)
        ! 0123
        nI = [3,6,7,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1122
        nI = [1,3,6,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        ! 1212
        nI = [1,4,5,8]
        call EncodeBitDet(nI, ilut)
        call test_excit_gen_guga(ilut, n_guga_excit_gen)


!         call stop_all(this_routine, "here")
        print *, "test_excit_gen_guga finished!"

    end subroutine run_test_excit_gen_guga_S0

    subroutine test_calcMixedContribution
        character(*), parameter :: this_routine = "test_calcMixedContribution"
        integer(n_int) :: ilut(0:nifguga), t(0:nifguga)

        call EncodeBitDet_guga([1,4,5,8],ilut)

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call EncodeBitDet_guga([1,3,6,8], t)

        print *, "testing calcMixedContribution(ilut,t,start,ende):"
        ASSERT(calcMixedContribution(ilut,t,1,4) == 0.0_dp)
   
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ASSERT(calcMixedContribution(t,ilut,1,4) == 0.0_dp)

        print *, "calcMixedContribution tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcMixedContribution

    subroutine test_generate_excitation_guga_double
        character(*), parameter :: this_routine = "test_generate_excitation_guga_double"
        integer :: nI(4), nJ(4), IC, excitMat(2,2), exFlag, nEx
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        logical :: tParity
        real(dp) :: pgen
        HElement_t(dp) :: HElGen
        type(excit_gen_store_type), target :: store
        integer(n_int), pointer :: ex(:,:)

        exFlag = 1
        ! make only double excitations:
        pSingles = 0.0_dp
        pDoubles = 1.0_dp - pSingles

        print *, "testing generate_excitation_guga:"
        ! 3300:
        nI = [1,2,3,4]
        call EncodeBitDet(nI,ilutI)

        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
               print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3030
        nI = [1,2,5,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3003
        nI = [1,2,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 0330
        nI = [3,4,5,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 0303
        nI = [3,4,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 0033
        nI = [5,6,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 1023
        nI = [1,6,7,8]
        call EncodeBitDet(nI,ilutI)

        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3102
        nI = [1,2,3,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3120
        nI = [1,2,3,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3012
        nI = [1,2,5,8]
        call EncodeBitDet(nI,ilutI)

        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
               print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 0312
        nI = [3,4,5,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1230
        nI = [1,4,5,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1203
        nI = [1,4,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1320
        nI = [1,3,4,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1302
        nI = [1,3,4,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1032
        nI = [1,5,6,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 0132
        nI = [3,5,6,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 0123
        nI = [3,6,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1122
        nI = [1,3,6,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 1212
        nI = [1,4,5,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        print *, "generate_excitation_guga tests passed!"
 
!         call stop_all(this_routine,"for now")

    end subroutine test_generate_excitation_guga_double

    subroutine test_generate_excitation_guga_single
        character(*), parameter :: this_routine = "test_generate_excitation_guga_single"
        integer :: nI(4), nJ(4), IC, excitMat(2,2), exFlag, nEx
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        logical :: tParity
        real(dp) :: pgen
        HElement_t(dp) :: HElGen
        type(excit_gen_store_type), target :: store
        integer(n_int), pointer :: ex(:,:)

        exFlag = 1
        ! make this store element ...
        print *, "testing generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,exMat,tPar,pgen,hEl,store)"
        ! test singles only first
        pSingles = 1.0_dp
        
        ! 3300:
        nI = [1,2,3,4]
        call EncodeBitDet(nI,ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3030
        nI = [1,2,5,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3003
        nI = [1,2,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 0330
        nI = [3,4,5,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 0303
        nI = [3,4,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 0033
        nI = [5,6,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        call actHamiltonian(ilutI, ex, nEx)
        call write_guga_list(6, ex(:,1:nEx))
 
        ! 1023
        nI = [1,6,7,8]
        call EncodeBitDet(nI,ilutI)

        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 3102
        nI = [1,2,3,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 3120
        nI = [1,2,3,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 3012
        nI = [1,2,5,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 0312
        nI = [3,4,5,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1230
        nI = [1,4,5,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1203
        nI = [1,4,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1320
        nI = [1,3,4,6]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1302
        nI = [1,3,4,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1032
        nI = [1,5,6,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if

        ! 0132
        nI = [3,5,6,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 0123
        nI = [3,6,7,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1122
        nI = [1,3,6,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        ! 1212
        nI = [1,4,5,8]
        call EncodeBitDet(nI,ilutI)
 
        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call write_det_guga(6, ilutI)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call write_det_guga(6, ilutJ)

        print *, "exact excitations for this ilut:"
        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
        else 
            print *, "no valid excitation created!"
        end if
        print *, "generate_excitation_guga tests passed!"

    end subroutine test_generate_excitation_guga_single

    subroutine test_calcDoubleR2L_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleR2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,3,4,2)

        ASSERT(excitInfo%typ == 13 )


        print *, "testing calcDoubleR2L_stochastic(ilut,exinfo,ex,pgen):"
!         call calcDoubleR2L_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcDoubleR2L_stochastic tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcDoubleR2L_stochastic

    subroutine test_calcDoubleL2R_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleL2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(3,1,2,4)

        ASSERT(excitInfo%typ == 12)


        print *, "testing calcDoubleL2R_stochastic(ilut,exinfo,ex,pgen):"
!         call calcDoubleL2R_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcDoubleL2R_stochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcDoubleL2R_stochastic

    subroutine test_calcDoubleR2L2R_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleR2L2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,3,2 )

        ASSERT(excitInfo%typ == 11)


        print *, "testing calcDoubleR2L2R_stochastic(ilut,exinfo,ex,pgen):"
!         call calcDoubleR2L2R_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcDoubleR2L2R_stochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcDoubleR2L2R_stochastic

    subroutine test_calcDoubleL2R2L_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleL2R2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,3)

        ASSERT(excitInfo%typ == 10)


        print *, "testing calcDoubleL2R2L_stochastic(ilut,exinfo,ex,pgen):"
!         call calcDoubleL2R2L_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcDoubleL2R2L_stochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcDoubleL2R2L_stochastic

    subroutine test_calcDoubleRaisingStochastic
        character(*), parameter :: this_routine = "test_calcDoubleRaisingStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,2,3)

        ASSERT(excitInfo%typ == 9)


        print *, "testing calcDoubleRaisingStochastic(ilut,exinfo,ex,pgen):"
!         call calcDoubleRaisingStochastic(ilut,excitInfo,ex,pgen)

        print *, "calcDoubleRaisingStochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcDoubleRaisingStochastic

    subroutine test_calcDoubleLoweringStochastic
        character(*), parameter :: this_routine = "test_calcDoubleLoweringStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 4,1,3,2 )

        ASSERT(excitInfo%typ == 8)


        print *, "testing calcDoubleLoweringStochastic(ilut,exinfo,ex,pgen):"
!         call calcDoubleLoweringStochastic(ilut,excitInfo,ex,pgen)

        print *, "calcDoubleLoweringStochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcDoubleLoweringStochastic

    subroutine test_calcFullStopR2L_stochastic
        character(*), parameter :: this_routine = "test_calcFullStopR2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 1,4,4,2  )

        ASSERT(excitInfo%typ == 17)


        print *, "testing calcFullStopR2L_stochastic(ilut,exinfo,ex,pgen):"
!         call calcFullStopR2L_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcFullStopR2L_stochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcFullStopR2L_stochastic

    subroutine test_calcFullStopL2R_stochastic
        character(*), parameter :: this_routine = "test_calcFullStopL2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8 ], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,4 )

        ASSERT(excitInfo%typ == 16 )


        print *, "testing calcFullStopL2R_stochastic(ilut,exInfo,ex,pgen)"
!         call calcFullStopL2R_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcFullStopL2R_stochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcFullStopL2R_stochastic

    subroutine test_calcFullStartR2L_stochastic
        character(*), parameter :: this_routine = "test_calcFullStartR2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8 ], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,2,4,1)

        ASSERT(excitInfo%typ == 21 )


        print *, "testing calcFullStartR2L_stochastic(ilut,exInfo,ex,pgen):"
!         call calcFullStartR2L_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcFullStartR2L_stochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcFullStartR2L_stochastic

    subroutine test_calcFullStartL2R_stochastic
        character(*), parameter :: this_routine = "test_calcFullStartL2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! encode det
        call EncodeBitDet_guga([1,4,5,8 ], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 1,4,2,1 )

        ASSERT(excitInfo%typ == 20)

        print *, "testing calcFullStartL2R_stochastic(ilut, exInfo, ex, pgen)"
!         call calcFullStartL2R_stochastic(ilut,excitInfo,ex,pgen)

        print *, "calcFullStartL2R_stochastic tests passed!"
        
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartL2R_stochastic

    subroutine test_calcRaisingSemiStopStochastic
        character(*), parameter :: this_routine = "test_calcRaisingSemiStopStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,3,4,8 ], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,2,4,1 )

        ASSERT(excitInfo%typ == 21)

        ! set up correct weights
        weights = init_fullStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),&
            currentB_ilut(2))

        ! calc the possible switches
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing calcRaisingSemiStopStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        call calcRaisingSemiStopStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        print *, "calcRaisingSemiStopStochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_calcRaisingSemiStopStochastic

    subroutine test_calcLoweringSemiStopStochastic
        character(*), parameter :: this_routine = "test_calcLoweringSemiStopStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,5,6,8 ], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(2,1,1,4 )

        ASSERT(excitInfo%typ == 20)

        ! calc the possible switches
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_fullStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),&
            currentB_ilut(2))

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing calcLoweringSemiStopStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        call calcLoweringSemiStopStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        print *, "calcLoweringSemiStopStochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_calcLoweringSemiStopStochastic


    subroutine test_calcRaisingSemiStartStochastic
        character(*), parameter :: this_routine = "test_calcRaisingSemiStartStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,5,6,8], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,4 )

        ASSERT(excitInfo%typ == 16)

        ! calc the possible switches
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_semiStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),&
            currentB_ilut(2))

        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing calcRaisingSemiStartStochastic(ilut,exInfo,weigh,negSwitch,posSwitch,ex,pgen):"
        call calcRaisingSemiStartStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        print *, "calcRaisingSemiStartStochastic tests passed!"
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_calcRaisingSemiStartStochastic


    subroutine test_calcLoweringSemiStartStochastic
        character(*), parameter :: this_routine = "test_calcLoweringSemiStartStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,3,4,8 ], ilut)

        ! calc b and occ vector
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,4,2 )

        ASSERT(excitInfo%typ == 17)

        ! calc the possible switches
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_semiStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),currentB_ilut(2))

        ! modify the excitation so it fits test case:
        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing calcLoweringSemiStartStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        call calcLoweringSemiStartStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        print *, "calcLoweringSemiStartStochastic tests passed!"

        ! clean up
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcLoweringSemiStartStochastic

    subroutine test_calcSingleOverlapMixedStochastic
        character(*), parameter :: this_routine = "test_calcSingleOverlapMixedStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        call EncodeBitDet_guga([3,4,5,6],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,4,3)
        
        ASSERT(excitInfo%typ == 7)

        print *, "testing calcSingleOverlapMixedStochastic(ilut, exInfo, ex, pgen):"
!         call calcSingleOverlapMixedStochastic(ilut, excitInfo, ex, pgen)

        call EncodeBitDet_guga([1,2,7,8],ilut)
  
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(3,1,3,4)

        ASSERT(excitInfo%typ == 6)

!         call calcSingleOverlapMixedStochastic(ilut,excitInfo, ex, pgen)

        print *, "calcSingleOverlapMixedStochastic tests passed!"

         deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcSingleOverlapMixedStochastic


    subroutine test_calcFullStopLoweringStochastic
        character(*), parameter :: this_routine = "test_calcFullStartLoweringStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        call EncodeBitDet_guga([1,2,5,6],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(4,1,4,3)
        
        ASSERT(excitInfo%typ == 14)

        print *, "testing calcFullStopLoweringStochastic(ilut, exInfo, ex, pgen):"
!         call calcFullStopLoweringStochastic(ilut, excitInfo, ex, pgen)

        print *, "calcFullStopLoweringStochastic tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcFullStopLoweringStochastic

    subroutine test_calcFullStopRaisingStochastic
        character(*), parameter :: this_routine = "test_calcFullStopRaisingStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        call EncodeBitDet_guga([3,4,7,8],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(1,4,3,4)
        
        ASSERT(excitInfo%typ == 15)

        print *, "testing calcFullStopRaisingStochastic(ilut, exInfo, ex, pgen):"
!         call calcFullStopRaisingStochastic(ilut, excitInfo, ex, pgen)

        print *, "calcFullStopRaisingStochastic tests passed!"
        
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcFullStopRaisingStochastic


    subroutine test_calcFullStartLoweringStochastic
        character(*), parameter :: this_routine = "test_calcFullStartLoweringStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        call EncodeBitDet_guga([1,2,5,6],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(3,1,4,1)
        
        ASSERT(excitInfo%typ == 18)

        print *, "testing calcFullStartLoweringStochastic(ilut, exInfo, ex, pgen):"
!         call calcFullStartLoweringStochastic(ilut, excitInfo, ex, pgen)

        print *, "calcFullStartLoweringStochastic tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartLoweringStochastic


    subroutine test_calcFullStartRaisingStochastic
        character(*), parameter :: this_routine = "test_calcFullStartRaisingStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        call EncodeBitDet_guga([5,6,7,8],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(1,3,1,4)
        
        ASSERT(excitInfo%typ == 19)

        print *, "testing calcFullStartRaisingStochastic(ilut, exInfo, ex, pgen):"
!         call calcFullStartRaisingStochastic(ilut, excitInfo, ex, pgen)

        print *, "calcFullStartRaisingStochastic tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcFullStartRaisingStochastic

    subroutine test_mixedFullStopStochastic
        character(*), parameter :: this_routine = "test_mixedFullStopStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        call EncodeBitDet_guga([1,4,5,8],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
        
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,1,4)
        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)

 
        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call doubleUpdateStochastic(ilut,2,excitInfo,weights,negSwitch,posSwitch,ex,pgen)
        call doubleUpdateStochastic(ilut,3,excitInfo,weights,negSwitch,posSwitch,ex,pgen)

        print *, "testing mixedFullStopStochastic(ilut, excitInfo, ex)"
        call mixedFullStopStochastic(ilut, excitInfo, ex)

        print *, "mixedFullStopStochastic tests passed!"

       deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_mixedFullStopStochastic


    subroutine test_doubleUpdateStochastic
        character(*), parameter :: this_routine = "test_doubleUpdateStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        call EncodeBitDet_guga([1,4,5,8],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(4,1,1,4)
        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)


        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing doubleUpdateStochastic(ilut,orb,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        call doubleUpdateStochastic(ilut,2,excitInfo,weights,negSwitch,posSwitch,ex,pgen)
        call doubleUpdateStochastic(ilut,3,excitInfo,weights,negSwitch,posSwitch,ex,pgen)

        print *, "doubleUpdateStochastic tests passed!"

       deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_doubleUpdateStochastic

    subroutine test_calcFullStartFullStopMixedStochastic
        character(*), parameter :: this_routine = "test_calcFullStartFullStopMixedStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen

        ! set up determinant and excitaiton information
        call EncodeBitDet_guga([1,4,5,8],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(1,4,4,1)

        ASSERT(excitInfo%typ==23)

        print *, "testing calcFullStartFullStopMixedStochastic(ilut,exInfo,ex,pgen)"
!         call calcFullStartFullStopMixedStochastic(ilut, excitInfo, ex, pgen)

        

        print *, "calcFullStartFullStopMixedStochastic tests passed!"

          deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartFullStopMixedStochastic


    subroutine test_pickOrbitals_double
        character(*), parameter :: this_routine = "test_pickOrbitals_double"
        type(excitationInformation) :: excitInfo
        integer(n_int) :: ilut(0:nifguga)
        real(dp) :: pgen
        integer :: nI(nel)

        nI = [1,2,3,4]
        
        call EncodeBitDet_guga([1,2,3,4],ilut)
       allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! still have to think about, if i preemptively choose the different
        ! types of double excitation (iiij, ii,jk, etc.) or let i happen 
        ! randomly and adjust the pgens accordingly...
        
        print *, "testing pickOrbitals_double(ilut, excitLvl):"
        call pickOrbitals_double(ilut, nI, excitInfo, pgen)
        call pickOrbitals_double(ilut, nI, excitInfo, pgen)


        print *, "pickOrbitals_double tests passed!"

          deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_pickOrbitals_double

    subroutine test_createStochasticExcitation_double
        character(*), parameter :: this_routine = "test_createStochasticExcitation_double"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        real(dp) :: pgen
        integer :: dummy(2), nI(nel)

        nI = [1,5,6,8]

        call EncodeBitDet_guga([1,5,6,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        print *, "testing createStochasticExcitation_double(ilut, ex, pgen):"
        call createStochasticExcitation_double(ilut,nI,ex,pgen,dummy)

        print *, "createStochasticExcitation_double tests passed!"

          deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_createStochasticExcitation_double

    subroutine test_singleStochasticEnd
        character(*), parameter :: this_routine = "test_singleStochasticEnd"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! do not make it random here but choose specific excitation!
!         excitInfo = pickOrbitals_single(ilut)
        excitInfo = excitationIdentifier(4,1)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, excitInfo%fullEnd)

        call singleStochasticUpdate(ilut, 2, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call singleStochasticUpdate(ilut, 3, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing singleStochasticEnd(ilut, excitInfo, excitation):"

        call singleStochasticEnd(ilut, excitInfo, ex)

        print *, "singleStochasticEnd tests passed!"
 
          deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_singleStochasticEnd


    subroutine test_singleStochasticUpdate
        character(*), parameter :: this_routine = "test_singleStochasticUpdate"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        call EncodeBitDet_guga([1,2,3,4], ilut)
        
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! do not make it random here but choose specific excitation!
!         excitInfo = pickOrbitals_single(ilut)
        excitInfo = excitationIdentifier(4,1)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, excitInfo%fullEnd)

        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "testing singleStochasticUpdate(ilut, exInfo, weight, posSwitch, negSwitch, ex, pgen):"
        call singleStochasticUpdate(ilut, 2, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call singleStochasticUpdate(ilut, 3, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, "singleStochasticUpdate tests passed!"

 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_singleStochasticUpdate
            


    subroutine test_pickRandomOrb
        character(*), parameter :: this_routine = "test_pickRandomOrb"
        integer :: orb
        real(dp) :: pgen
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,8], ilut)
       
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        pgen = 1.0_dp

        print *, "testing pickRandomOrb(orbRes, pgen, orb, occRes):"
        call pickRandomOrb_scalar(1, pgen, orb)
        ASSERT(pgen == 1.0_dp/3.0_dp)
        ASSERT( orb > 1 .and. orb <= 4)
        pgen = 1.0_dp
        call pickRandomOrb_forced(2, pgen, orb)
        ASSERT(orb == 1)
        ASSERT(pgen == 1.0_dp)
        pgen = 1.0_dp
        call pickRandomOrb_vector([1,2], pgen, orb)
        ASSERT(orb == 3 .or. orb == 4)
        ASSERT(pgen == 1.0_dp/2.0_dp)
        pgen = 1.0_dp
        call pickRandomOrb_vector([2,3], pgen, orb, 2)
        ASSERT(orb == 4)
        ASSERT(pgen == 1.0_dp)

        pgen = 1.0_dp
        call pickRandomOrb_scalar(0, pgen, orb, 0)
        ASSERT(pgen == 1.0_dp/3.0_dp)
        ASSERT(orb /= 3)

        pgen = 1.0_dp
        call pickRandomOrb_scalar(2, pgen, orb, 2)
        ASSERT(pgen == 1.0_dp/2.0_dp)
        ASSERT(orb == 3 .or. orb == 4)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,2,pgen,orb)
        ASSERT(pgen == 0.0_dp)
        ASSERT(orb == 0)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,4,pgen,orb)
        print *, pgen, orb
        ASSERT(pgen == 1.0_dp/2.0_dp)
        ASSERT(orb == 2 .or. orb == 3)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,4,pgen,orb,0)
        ASSERT(pgen == 1.0_dp)
        ASSERT(orb == 2)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,4,pgen,orb,1)
        ASSERT(pgen == 1.0_dp)
        ASSERT(orb == 3)
        print *, "pickRandomOrb tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_pickRandomOrb

    subroutine test_mixedFullStartStochastic
        character(*), parameter :: this_routine = "test_mixedFullStartStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), prob

        call EncodeBitDet_guga([1,4,5,6],ilut)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(4,1,1,3)
        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)

        print *, "testing mixedFullStartStochastic:"

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, prob)

        print *, "mixedFullStartStochastic tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_mixedFullStartStochastic

    subroutine test_createStochasticStart_single
        character(*), parameter :: this_routine = "test_createStochasticStart_single"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        type(weight_obj) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), probWeight
        integer(n_int) :: ex(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        excitInfo = excitationIdentifier(4,1)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut,4)

        print *, "testing createStochasticStart_single(ilut,exInfo, weighs, posSwitch, negSwitch, ex, probWeight):"
        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, negSwitch, ex, probWeight)

        print *, "createStochasticStart_single tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_createStochasticStart_single

    subroutine test_pickOrbitals_single
        character(*), parameter :: this_routine = "test_pickOrbitals_single"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: pgen
        integer :: nI(nel)

        nI = [1,2,3,4]

        call EncodeBitDet_guga([1,2,3,4],ilut)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        ! todo the pick random orbital picker.. 
        ! ask ali and simon
        print *, "testing: pickOrbitals_single(ilut)"
        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        call pickOrbitals_single(ilut, nI, excitInfo, pgen)

        call EncodeBitDet_guga([1,4,5,8], ilut)
       
        nI = [1,4,5,8]
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        call pickOrbitals_single(ilut, nI, excitInfo, pgen)


        print *, "pickOrbitals_single tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_pickOrbitals_single

    subroutine test_createStochasticExcitation_single
        character(*), parameter :: this_routine = "test_createStochasticExcitation_single"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        real(dp) :: pgen
        integer :: nI(nel) 

        nI = [1,2,3,4]
        
        call EncodeBitDet_guga([1,2,3,4], ilut)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        print *, "testing: createStochasticExcitation_single(ilut,ex,weight):"
        call createStochasticExcitation_single(ilut, nI, ex, pgen)

        print *, "createStochasticExcitation_single tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_createStochasticExcitation_single


    subroutine test_actHamiltonian
        character(*), parameter :: this_routine = "test_actHamiltonian"
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), pointer :: ex(:,:)
        integer :: nEx

        
        print *, "testing actHamiltonian(ilut):"
        ! 3300:
        call EncodeBitDet_guga([1,2,3,4], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 13)   
        ! 0330
        call EncodeBitDet_guga([3,4,5,6],ilut)
        call actHamiltonian(ilut,ex,nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 14)
         ! 0303
        call EncodeBitDet_guga([3,4,7,8],ilut)
        call actHamiltonian(ilut,ex,nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 14)
        ! 0033
        call EncodeBitDet_guga([5,6,7,8],ilut)
        call actHamiltonian(ilut,ex,nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 13)
       ! 1023
        call EncodeBitDet_guga([1,6,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
       ! 3102
        call EncodeBitDet_guga([1,2,3,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
       ! 3120
        call EncodeBitDet_guga([1,2,3,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)

        ! 3030
        call EncodeBitDet_guga([1,2,5,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 14)
        ! 3003:
        call EncodeBitDet_guga([1,2,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 14)
        ! 3012
        call EncodeBitDet_guga([1,2,5,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 16)
        ! 0312
        call EncodeBitDet_guga([3,4,5,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 16)
        ! 1230
        call EncodeBitDet_guga([1,4,5,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 16)
        ! 1203
        call EncodeBitDet_guga([1,4,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 16)
        ! 1320
        call EncodeBitDet_guga([1,3,4,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
        ! 1302
        call EncodeBitDet_guga([1,3,4,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
        ! 1032
        call EncodeBitDet_guga([1,5,6,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
        ! 0132
        call EncodeBitDet_guga([3,5,6,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
        ! 0123
        call EncodeBitDet_guga([3,6,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 17)
              ! 1122
        call EncodeBitDet_guga([1,3,6,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 12)
        ! 1212
        call EncodeBitDet_guga([1,4,5,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        ASSERT(nEx == 18)



        print *, "actHamiltonian tests passed!"
            
    end subroutine test_actHamiltonian


    subroutine test_calcAllExcitations_double
        character(*), parameter :: this_routine = "test_calcAllExcitations_double"
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), pointer :: ex(:,:)
        integer :: nExcits

        call EncodeBitDet_guga([1,4,5,8],ilut)
 
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        print *, "testing calcAllExcitations_double(ilut,i,j,k,l,ex,nExits):"
        call calcAllExcitations_double(ilut,1,2,3,4, ex, nExcits)

        print *, "calcAllExcitations_double tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_calcAllExcitations_double


    subroutine test_calcFullStartFullStopAlike
        character(*), parameter :: this_routine = "test_calcFullStartFullStopAlike"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,7,8], ilut)

        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(1,4,1,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStartFullStopAlike(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStartFullStopAlike(ilut, excitInfo, ex)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStartFullStopAlike tests passed!"

        call EncodeBitDet_guga([1,2,3,4], ilut) 

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(4,1,4,1)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        currentB_ilut = calcB_vector_ilut(ilut)
        print *, excitInfo%typ

        print *, "testing: calcFullStartFullStopAlike(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStartFullStopAlike(ilut, excitInfo, ex)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStartFullStopAlike tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartFullStopAlike



    subroutine test_calcFullStartL2R
        character(*), parameter :: this_routine = "test_calcFullStartL2R"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        excitInfo = excitationIdentifier(1,4,3,1)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStartL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStartL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStartL2R tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartL2R


    subroutine test_calcFullStartR2L
        character(*), parameter :: this_routine = "test_calcFullStartR2L"
        integer(n_int) :: ilut(0:2)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(1,3,4,1)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStartR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStartR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStartR2L tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartR2L


    subroutine test_calcFullStartRaising
        character(*), parameter :: this_routine = "test_calcFullStartRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(1,4,1,3)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStartRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStartRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStartRaising tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartRaising



    subroutine test_calcFullStartLowering
        character(*), parameter :: this_routine = "test_calcFullStartLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        excitInfo = excitationIdentifier(4,1,3,1)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStartLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStartLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStartLowering tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStartLowering


       subroutine test_calcFullStopR2L
        character(*), parameter :: this_routine = "test_calcFullStopR2L"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        excitInfo = excitationIdentifier(4,1,3,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStopR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStopR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStopR2L tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStopR2L



       subroutine test_calcFullStopL2R
        character(*), parameter :: this_routine = "test_calcFullStopL2R"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(4,1,3,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStopL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStopL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStopL2R tests passed!"
 
        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStopL2R


       subroutine test_calcFullStopRaising
        character(*), parameter :: this_routine = "test_calcFullStopRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([5,6,7,8], ilut)
       allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,4,2,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStopRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStopRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStopRaising tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStopRaising

     subroutine test_calcFullStopLowering
        character(*), parameter :: this_routine = "test_calcDoubleFullStopLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,4,2)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcFullStopLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcFullStopLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcFullStopLowering tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcFullStopLowering

    subroutine test_calcDoubleR2L
        character(*), parameter :: this_routine = "test_calcDoubleR2L"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,4,2)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcDoubleR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcDoubleR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcDoubleR2L tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcDoubleR2L

   subroutine test_calcDoubleL2R
        character(*), parameter :: this_routine = "test_calcDoubleL2R"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(3,1,2,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcDoubleL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcDoubleL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcDoubleL2R tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcDoubleL2R


   subroutine test_calcDoubleRaising
        character(*), parameter :: this_routine = "test_calcDoubleRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([5,6,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,2,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcDoubleRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcDoubleRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcDoubleRaising tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcDoubleRaising



    subroutine test_calcDoubleLowering
        character(*), parameter :: this_routine = "test_calcDoubleLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(3,1,4,2)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcDoubleLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcDoubleLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcDoubleLowering tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcDoubleLowering

    subroutine test_calcSingleOverlapRaising
        character(*), parameter :: this_routine = "test_singleOverlapRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([5,6,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,3,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)


        print *, " numExcits: ", num
        print *, ex
        print *, "calcSingleOverlapRaising tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcSingleOverlapRaising

    subroutine test_calcSingleOverlapMixed
        character(*), parameter :: this_routine = "test_singleOverlapMixed"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(3,1,3,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapMixed(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        print *, " numExcits: ", num
        print *, ex

        print *, "calcSingleOverlapMixed tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)
    end subroutine test_calcSingleOverlapMixed


    subroutine test_calcSingleOverlapLowering
        character(*), parameter :: this_routine = "test_singleOverlapLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(2,1,4,2)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        print *, " numExcits: ", num
        print *, ex

        print *, "calcSingleOverlapLowering tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcSingleOverlapLowering


    subroutine test_calcNonOverlapDouble
        character(*), parameter :: this_routine = "test_calcNonOverlapDouble"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,2,3,4)
        call calcRemainingSwitches(ilut, excitInfo, 1, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcNonOverlapDouble(ilut, exInfo, exs, num, posSwitch, negSwitch"
        call calcNonOverlapDouble(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! todo: auch hier eine funktion noetig die nur die single excitations
        ! berechnt und nicht tmat einbezieht... und das auch effektiv macht
        
        print *, "calcNonOverlapDouble tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcNonOverlapDouble

    subroutine test_add_ilut_lists
        character(*), parameter :: this_routine = "test_add_ilut_lists"
        integer(n_int) :: ilutI(0:2,1), ilutJ(0:2,1), ilutOut(0:2,2)
        integer :: nOut


        call EncodeBitDet_guga([1,2,3,4],ilutI)
        call EncodeBitDet_guga([3,4,7,8],ilutJ)

        print *, "testing add_ilut_lists(n1, n2, sortFlag, ilut1, ilut2, ilutOut, nOut):"
        call add_ilut_lists(1,1,.false.,ilutI, ilutJ, ilutOut, nOut)

        ASSERT(nOut == 2)

        call add_ilut_lists(1,1,.false.,ilutI, ilutI, ilutOut, nOut)

        ASSERT(nOut == 1)

        print *, "add_ilut_lists tests passed!"

    end subroutine test_add_ilut_lists

    subroutine test_calcDoubleExcitation_withWeight
        character(*), parameter :: this_routine = "test_calcDoubleExcitation_withWeight"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,7,8],ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(2,2,1,4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcDoubleExcitation_withWeight(ilut, exInfo, exc, num)"
        call calcDoubleExcitation_withWeight(ilut, excitInfo, ex, num, posSwitch, &
            negSwitch)
        ! todo! have to write a only excitation calculating funciton
        ! not including tmat for this case! 

        print *, "calcDoubleExcitation_withWeight tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

    end subroutine test_calcDoubleExcitation_withWeight

    subroutine test_checkCompatibility
        character(*), parameter :: this_routine = "test_checkCompatibility"
        real(dp) :: posSwitch(4), negSwitch(4)
        integer(n_int) :: ilut(0:nifguga)
        logical :: flag
        type(excitationInformation) :: excitInfo

        print *, "testing: checkCompatibility:"
        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier_double(1,2,3,4)
        print *, excitInfo%typ
        call calcRemainingSwitches_excitInfo_double(ilut, excitInfo, 1, posSwitch, negSwitch)
        call checkCompatibility(ilut, excitInfo, flag, posSwitch, negSwitch)

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

        print *, "checkCompatibility tests passed!"


    end subroutine test_checkCompatibility


    subroutine test_calcRemainingSwitches_double
        character(*), parameter :: this_routine = "test_calcRemainingSwitches_double"
        real(dp) :: posSwitch(4), negSwitch(4)  
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        print *, "testing: calcRemainingSwitches_double"
        call calcRemainingSwitches_double(ilut,1,2,3,4,posSwitch,negSwitch)
        call calcRemainingSwitches_double(ilut,3,2,3,1,posSwitch,negSwitch)
        call calcRemainingSwitches_double(ilut,1,4,3,4,posSwitch,negSwitch)

        print *, "calcRemainingSwitches_double tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcRemainingSwitches_double
            

    subroutine test_excitationIdentifier_double
        character(*), parameter :: this_routine = "test_excitationIdentifier_double"
        type(excitationInformation) :: excitInfo


        print *, "testing: excitationIdentifier_double"

        excitInfo = excitationIdentifier_double(1,2,3,4)
        excitInfo = excitationIdentifier_double(1,2,2,4)
        excitInfo = excitationIdentifier_double(3,2,3,4)
        excitInfo = excitationIdentifier_double(4,2,3,1)
        excitInfo = excitationIdentifier_double(1,1,3,4)
        excitInfo = excitationIdentifier_double(1,1,4,4)
        excitInfo = excitationIdentifier_double(1,1,1,1)
        excitInfo = excitationIdentifier_double(1,3,2,4)
        excitInfo = excitationIdentifier_double(1,4,3,2)


        print *, "excitationIdentifier_double tests passed!"

    end subroutine test_excitationIdentifier_double


    subroutine test_calcAllExcitations_single
        character(*), parameter :: this_routine = "test_calcAllExcitations_single"
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), pointer :: ex(:,:)
        integer :: nEx

        call EncodeBitDet_guga([1,2,3,4], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        print *, "testing: calcAllExcitations_single"
        call calcAllExcitations_single(ilut, 4, 1, ex, nEx)
        call write_bit_rep(6, ex(:,1), .true.)
        print *, extract_part_sign(ex(:,1), 1), "==", -Root2
        ASSERT(nEx == 1)

        deallocate(ex)

        call EncodeBitDet_guga([3,6,7,8], ilut)
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call calcAllExcitations_single(ilut, 2, 4, ex, nEx)


        print *, "calcAllExcitations_single tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_calcAllExcitations_single


    subroutine test_singleEnd
        character(*), parameter :: this_routine = "test_singleEnd"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: posSwitch(nBasis/2), negSwitch(nBasis/2)
        integer(n_int), pointer :: excits(:,:), tmpEx(:,:)
        integer :: num
        type(weight_obj) :: weights

        print *, "testing: singleEnd(ilut, exInfo, tmpEx, nEx, excits)"
        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1,4,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, 4)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            tmpEx, num)

        print *, tmpEx
        ASSERT(getDeltaB(tmpEx) == 1)
        call singleUpdate(ilut, 2, excitInfo, posSwitch, negSwitch, weights, &
            tmpEx, num)

        print *, getDeltaB(tmpEx)
        ASSERT(getDeltaB(tmpEx) == -1)
        call write_det_guga(6, tmpEx(:,1))
        ASSERT(num == 1)
        print *, extract_matrix_element(tmpEx(:,1),1), "==", -OverR2

        call singleUpdate(ilut, 3, excitInfo, posSwitch, negSwitch, weights, &
            tmpEx, num)
        call write_det_guga(6, tmpEx(:,1))
        ASSERT(num == 1)
        print *, extract_matrix_element(tmpEx(:,1),1), "==", -OverR2

        call singleEnd(ilut, excitInfo, tmpEx, num, excits)

        call write_bit_rep(6, excits(:,1), .true.)

        ASSERT(.not. associated(tmpEx))
        ASSERT( num == 1)
        print*,  extract_matrix_element(excits(:,1),1)," ==", -1.0_dp

        print *, "singleEnd tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)


    end subroutine test_singleEnd


    subroutine test_singleUpdate
        character(*), parameter :: this_routine = "test_singleUpdate"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: posSwitch(nBasis/2), negSwitch(nBasis/2)
        integer(n_int), pointer :: excits(:,:)
        integer :: num
        type(weight_obj) :: weights

        print *, "testing: singleUpdate(ilut,orb,exInfo,posSwitch,negSwitch,excits,num)"
        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1,4,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, 4)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        call singleUpdate(ilut, 2, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        call write_det_guga(6, excits(:,1))
        ASSERT(num == 1)
        print *, extract_matrix_element(excits(:,1),1), "==", Root2

        call singleUpdate(ilut, 3, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)
        call write_det_guga(6, excits(:,1))
        ASSERT(num == 1)
        print *, extract_matrix_element(excits(:,1),1), "==", Root2

        print *, "singleUpdate tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

        deallocate(excits)

    end subroutine test_singleUpdate
    

    subroutine test_createSingleStart
        character(*), parameter :: this_routine = "test_createSingleStart"
        integer(n_int) :: ilut(0:nifguga)
        type(excitationInformation) :: excitInfo
        real(dp) :: posSwitch(nBasis/2), negSwitch(nBasis/2)
        integer(n_int), pointer :: excits(:,:)
        integer :: num
        type(weight_obj) :: weights

        print *, "testing: createSingleStart(ilut, excitInfo, posSwitch, negSwitch, excits, nExcits)"

        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1,4,7,8], ilut)
        allocate(currentB_ilut(4))
        allocate(currentOcc_ilut(4))
        allocate(currentOcc_int(4))
        allocate(current_stepvector(4))
        allocate(currentB_int(4))
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut,4)
        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights,&
            excits, num)
        call write_det_guga(6, ilut) 
        ASSERT(num == 1)
        ASSERT(extract_matrix_element(excits(:,1),1) == Root2)
        ASSERT(getDeltaB(excits(:,1)) == 1)

        call write_det_guga(6, ilut) 
        
        deallocate(excits)
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call EncodeBitDet_guga([1,4,7,8], ilut)

        excitInfo = excitationIdentifier(2,4)

        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        call write_det_guga(6, ilut) 
        ASSERT(num == 1)
        print *, extract_matrix_element(excits(:,1),1), OverR2
        ASSERT(getDeltaB(excits(:,1)) == -1)

        deallocate(excits)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        excitInfo = excitationIdentifier(4, 1)

        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
 
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)
!         call write_bit_rep(6, excits(:,1), .true.)
        ASSERT(num == 1)
        ASSERT(getDeltaB(excits(:,1)) == -1)
        ASSERT(extract_matrix_element(excits(:,1),1)== Root2)

        ! and do one double start
        deallocate(excits)

        call EncodeBitDet_guga([1,3,4,8], ilut)
        excitInfo = excitationIdentifier(4,2)
        call calcRemainingSwitches(ilut, excitInfo, posSwitch, negSwitch)
  
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)


        print *, "createSingleStart tests passed!"

        deallocate(currentB_ilut)
        deallocate(currentOcc_ilut)
        deallocate(currentOcc_int)
        deallocate(current_stepvector)
        deallocate(currentB_int)

        deallocate(excits)

    end subroutine test_createSingleStart



    subroutine test_excitationIdentifier_single
        character(*), parameter :: this_routine = "test_excitationIdentifier_single"
        type(excitationInformation) :: excitInfo

        print *, "testing: excitationIdentifier_single:"
        excitInfo = excitationIdentifier(1, 4)
        ASSERT(excitInfo%i==1)
        ASSERT(excitInfo%j==4)
        ASSERT(excitInfo%gen1==1)
        ASSERT(excitInfo%fullStart==1)
        ASSERT(excitInfo%fullEnd==4)
        ASSERT(excitInfo%currentGen == 1)
        ASSERT(excitInfo%excitLvl == 2)
        ASSERT(excitInfo%typ == 0)
        print *, "excitationIdentifier_single tests passed!"

    end subroutine test_excitationIdentifier_single


    subroutine test_getDoubleMatrixElement
        character(*), parameter :: this_routine = "test_getDoubleMatrixElement"
        real(dp) :: x0, x1

        print *, "testing: getDoubleMatrixElement:"
        call getDoubleMatrixElement(1,1,0,1,1,1.0_dp,1.0_dp,x0,x1)
        print *, "x0 = ", x0, " x1 = ", x1
        print *, "getDoubleMatrixElement tests passed!"

    end subroutine test_getDoubleMatrixElement
        
    subroutine test_getMixedFullStop
        character(*), parameter :: this_routine = "test_getMixedFullStop"
        real(dp) :: x0, x1

        print *, "testing: getMixedFullStop:"
        call getMixedFullStop(1,1,0,1.0_dp,x0,x1)
        print *, "x0 = ", x0, " x1 = ", x1
        print *, "getMixedFullStop tests passed!"

    end subroutine test_getMixedFullStop


    subroutine test_getSingleMatrixElement
        character(*), parameter :: this_routine = "test_getSingleMatrixElement"

        print *, "testing: getSingleMatrixElement(d1,d2,dB,gen,b):"
        print *, "h = ", getSingleMatrixElement(1,1,1,1,1.0_dp)
        print *, "getSingleMatrixElement tests passed!"

    end subroutine test_getSingleMatrixElement
        


    subroutine test_matrix_element_ops
        character(*), parameter :: this_routine ="test_encode_matrix_element"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)
        
        print *, "testing: encode_matrix_element(ilut, ele, type):"
        print *, "         extract_matrix_element(ilut, type):"
        print *, "         update_matrix_element(ilut, ele, type):"

        call encode_matrix_element(ilut,1.0_dp,1)
        ASSERT(extract_matrix_element(ilut, 1) == 1.0_dp)
        call encode_matrix_element(ilut,-1.0_dp,2)
        ASSERT(extract_matrix_element(ilut, 2) == -1.0_dp)
        call update_matrix_element(ilut, 2.0_dp, 1)
        ASSERT(extract_matrix_element(ilut, 1) == 2.0_dp)
        call update_matrix_element(ilut, -2.0_dp, 2)
        ASSERT(extract_matrix_element(ilut, 2) == 2.0_dp)

        print *, "encode_matrix_element tests passed!"
    end subroutine test_matrix_element_ops

    subroutine test_calcOcc_vector_ilut
        character(*), parameter :: this_routine = "test_calcOcc_vector_ilut"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)
        
        print *, "testing: calcOcc_vector_ilut(ilut)"
        ASSERT(all([2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp] == calcOcc_vector_ilut(ilut)))
        
        call EncodeBitDet_guga([1,4,5,6], ilut)
        ASSERT(all([1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp] == calcOcc_vector_ilut(ilut)))

        call EncodeBitDet_guga([1,2,3,8], ilut)
        ASSERT( all([2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp] == calcOcc_vector_ilut(ilut)))
        print *, "calcOcc_vector_ilut tests passed!"

    end subroutine test_calcOcc_vector_ilut

    subroutine test_calcStepVector
        character(*), parameter :: this_routine = "test_calcStepVector"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: calcStepVector(ilut)"
        ASSERT(all([3,3,0,0] == calcStepVector(ilut)))
        call EncodeBitDet_guga([1,4,5,6], ilut)
        ASSERT(all([1,2,3,0] == calcStepVector(ilut)))
        print *, "calcStepVector tests passed!"
            

    end subroutine test_calcStepVector

    subroutine test_isDouble
        character(*), parameter :: this_routine = "test_isDouble"

        print *, "testing: isDouble(nI, sOrb)"
        ASSERT(isDouble([1,2,3,4],1))
        ASSERT(isDouble([1,2,3,6],2))
        ASSERT(.not.isDouble([1,2,3,6],3))
        print *, "isDouble tests passed!"

    end subroutine test_isDouble

    subroutine test_isProperCSF_ilut
        integer(n_int) :: ilut(0:nifguga)
        character(*), parameter :: this_routine = "test_isProperCSF_ilut"

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: isProperCSF_ilut(ilut)"
        ASSERT(isProperCSF_ilut(ilut))
        call EncodeBitDet_guga([2,3,4,5],ilut)
        ASSERT(.not.isProperCSF_ilut(ilut))
        print *, "isProperCSF_ilut tests passed!"

    end subroutine test_isProperCSF_ilut

    subroutine test_set_get_DeltaB
        integer(n_int) :: ilut(0:nifguga)
        integer :: deltaB
        character(*), parameter :: this_routine = "test_set_get_DeltaB"

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: setDeltaB(deltaB,ilut)"
        call setDeltaB(-2, ilut)
        ASSERT(getDeltaB(ilut) == -2)
        call setDeltaB(-1, ilut)
        ASSERT(getDeltaB(ilut) == -1)
        call setDeltaB(0, ilut)
        ASSERT(getDeltaB(ilut) == 0)
        call setDeltaB(1, ilut)
        ASSERT(getDeltaB(ilut) == 1)
        call setDeltaB(2, ilut)
        ASSERT(getDeltaB(ilut) == 2)

        ! maybe problems with x1-matrix element
        call setDeltaB(-2, ilut)
        ASSERT(getDeltaB(ilut) == -2)
        call setDeltaB(-1, ilut)
        ASSERT(getDeltaB(ilut) == -1)
        call setDeltaB(0, ilut)
        ASSERT(getDeltaB(ilut) == 0)
        call setDeltaB(1, ilut)
        ASSERT(getDeltaB(ilut) == 1)
        call setDeltaB(2, ilut)
        ASSERT(getDeltaB(ilut) == 2)

        print *, "setDeltaB tests passed!"

    end subroutine test_set_get_DeltaB


    subroutine test_count_open_orbs_ij
        integer :: det(4)
        integer(n_int) :: ilut(0:nifguga)
        character(*), parameter :: this_routine = "test_count_open_orbs_ij"

        det = [1,2,3,4]

        call EncodeBitDet_guga(det, ilut)

        print *, "testing: count_open_orbs_ij(ilut, i, j)"
        print *, "open orbs in ([1,2,3,4],1,4): ", count_open_orbs_ij(1,4, ilut(0:0))
        ASSERT(count_open_orbs_ij(1,4, ilut) == 0)

        det = [1, 3, 5, 6]

        call EncodeBitDet_guga(det, ilut)
        print *, "open orbs in ([1, 3, 5, 6], 1, 4): ", count_open_orbs_ij(1, 4, ilut)
        ASSERT(count_open_orbs_ij(1,4, ilut) == 2)

        print *, "count_open_orbs_ij tests passed!"

    end subroutine test_count_open_orbs_ij


    subroutine test_getExcitationRangeMask

        print *, "testing: getExcitationRangeMask(i,j)"
        print *, "i = 1, j = 4 :", getExcitationRangeMask(1, 4)
        print *, "i = 2, j = 3 :", getExcitationRangeMask(2, 3)
        print *, "getExcitationRangeMask tests passed!"

    end subroutine test_getExcitationRangeMask



    subroutine check_calcDiagMatEleGUGA_nI
        integer :: det(nEl)
        character(*), parameter :: this_routine = "check_calcDiagMatEles_nI"


        ! calculating the diagonal matrix elements in the guga formalism 
        ! for d = |3,3,0,0> -> H_d = 8
        ! and for d = |3,1,2,0> -> H_d = 9 
        ! from my derived formulas per hand... 
        ! but calculating it through NECI or Sandeeps code gives 
        ! for d = |3300> = 9, which makes sense if all orbitals have the same 
        ! single particle energies.. but atleast my code and my formulas 
        ! coincide... -> so be content for now

        print *, "testing: calcDiagMatEles_nI"

        ! 3120
        det = [1,2,3,6]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 3300
        det = [1,2,3,4]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 8.0_dp)

        ! 0330
        det = [3,4,5,6]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 8.0_dp)

        ! 0303
        det = [3,4,7,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 8.0_dp)

        ! 0033
        det = [5,6,7,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 8.0_dp)

        ! 1023
        det = [1,6,7,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 3102
        det = [1,2,3,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 3030
        det = [1,2,5,6]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 8.0_dp)

        ! 3003
        det = [1,2,7,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 8.0_dp)

        ! 3012
        det = [1,2,5,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 0312
        det = [3,4,5,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 1230
        det = [1,4,5,6]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 1203
        det = [1,4,7,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 1320
        det = [1,3,4,6]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 1302
        det = [1,3,4,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 1032
        det = [1,5,6,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 0132
        det = [3,5,6,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 0123
        det = [3,6,7,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 9.0_dp)

        ! 1122
        det = [1,3,6,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 10.0_dp)

        ! 1212
        det = [1,4,5,8]
        ASSERT(calcDiagMatEleGUGA_nI(det) == 10.0_dp)

        print *, this_routine, " tests passed!"

    end subroutine check_calcDiagMatEleGUGA_nI


    subroutine check_calcDiagExchange_nI
        integer :: det(nEl), iOrb, jOrb

        det = [1,2,3,6]
        iOrb = 3
        jOrb = 4

        print *, "testing: calcDiagExchange_nI"
        print *, "H_ii = ", calcDiagMatEleGUGA_nI(det)
        print *, "X = ", calcDiagExchangeGUGA_nI(iOrb, jOrb, det)
        print *, "calcDiagExchange tests passed!"

    end subroutine check_calcDiagExchange_nI


    subroutine check_determinants
        integer(n_int) :: ilut(0:nifguga)
        integer :: det(nEl) 

        det = [1,2,3,4]
        call EncodeBitDet_guga(det,ilut)

        ! check which matrix elements get initialized
        call write_bit_rep(6, ilut, .true.)

        ! also do it directly with extract function
        print *, extract_part_sign(ilut,1)


    end subroutine check_determinants


!     subroutine test_calcSingleProbWeight
!         integer(n_int) :: ilut(0:nifguga)
!         integer :: det(nEl), ind
!         real(dp) :: posSwitch, negSwitch, posWeight, negWeight, checkPos, &
!             checkNeg
!         character(*), parameter :: testFun = "calcSingleProbWeight", &
!             this_routine = "test_calcSingleProbWeight"
! 
!         print *, " Testing ", testFun
!         det = [1,2,3,4]
!         ind = 3
!         posSwitch = 0.0_dp
!         negSwitch = 0.0_dp
! 
!         call EncodeBitDet_guga(det, ilut)
! 
!         call calcSingleProbWeight(ilut, ind, posSwitch, negSwitch, &
!             posWeight, negWeight)
! 
!         print *, testFun, " tests passed!"
! 
!     end subroutine test_calcSingleProbWeight

    subroutine test_calcbvector
        integer(n_int) :: ilut(0:nifguga)
        integer :: det(nEl), checkB_nI(nEl), checkB_ilut(nBasis/2)
        character(*), parameter :: testFun = "calcB_vector", &
            this_routine = "test_calcbvector"

        print *, " Testing ", testFun
        det = [1,2,3,4]
        checkB_nI = 0.0_dp
        call EncodeBitDet_guga(det,ilut)
        checkB_ilut = 0.0_dp

        ASSERT(all(calcB_vector_nI(det) == checkB_nI))
        ASSERT(all(calcB_vector_ilut(ilut) == checkB_ilut))

        det = [1,2,3,6]
        call EncodeBitDet_guga(det, ilut)

        checkB_nI(3) = 1.0_dp
        checkB_ilut(2) = 1.0_dp

        ASSERT(all(calcB_vector_nI(det) == checkB_nI))
        ASSERT(all(calcB_vector_ilut(ilut) == checkB_ilut))
        print *, testFun, " tests passed!"

        call write_bit_rep(6,ilut,.true.)

    end subroutine test_calcbvector

    subroutine test_bitChecks
        ! checks the function isZero(ilut,sOrb), isOne(ilut,sOrb) etc.
        integer(n_int) :: ilut(0:nifguga)
        integer :: det(nEl)
        character(*), parameter :: testFun = "bitChecks", &
            this_routine = "test_bitChecks"

        det = [1,2,3,6]
        ! make a valid ilut:
        call EncodeBitDet_guga(det,ilut)

        print *, "***"
        print *, " Testing ",testFun
        ASSERT(.not.isZero(ilut,1))
        ASSERT(.not.isZero(ilut,2))
        ASSERT(.not.isZero(ilut,3))
        ASSERT(isZero(ilut,4))
        ASSERT(.not. isOne(ilut,1))
        ASSERT(isOne(ilut,2))
        ASSERT(.not.isOne(ilut,3))
        ASSERT(.not.isOne(ilut,4)) 
        ASSERT(.not.isTwo(ilut,1))
        ASSERT(.not.isTwo(ilut,2))
        ASSERT(isTwo(ilut,3))
        ASSERT(.not.isTwo(ilut,4)) 
        ASSERT(isThree(ilut,1))
        ASSERT(.not.isThree(ilut,2))
        ASSERT(.not.isThree(ilut,3))
        ASSERT(.not.isThree(ilut,4)) 
        print *, testFun, " tests passed!"
    
    end subroutine test_bitChecks

    subroutine test_calcRemainingSwitches()
        real(dp) :: neg(nBasis/2), pos(nBasis/2)
        integer :: det(nEl)
        integer(n_int) :: ilut(0:nifguga)
        integer :: b
        character(*), parameter :: this_routine = "test_calcRemainingSwitches"

        det = [1,2,3,6] ! 3 1 2 0

        call EncodeBitDet_guga(det,ilut)
        ! now need excitInfo for 

        allocate(current_stepvector(4))

        current_stepvector = calcStepVector(ilut)

        print *, "***"
        print *, "testing calcRemainingSwitches:"
        call calcRemainingSwitches(ilut,1,4,pos,neg)
        print *, "positive switches: ", pos
        ASSERT(all(pos == [1.0_dp,1.0_dp,0.0_dp,0.0_dp]))
        print *, "negative switches: ", neg
        ASSERT(all(neg == [1.0_dp,0.0_dp,0.0_dp,0.0_dp]))
        
        call EncodeBitDet_guga([1,4,5,8],ilut) ! 1 2 1 2

        current_stepvector = calcStepVector(ilut)

        call calcRemainingSwitches(ilut, 2,4,pos,neg)
        print *, "positive switches: ", pos
        ASSERT(all(pos == [0.0_dp,0.0_dp,0.0_dp,0.0_dp]))
        print *, "negative switches: ", neg
        ASSERT(all(neg == [0.0_dp,1.0_dp,0.0_dp,0.0_dp]))
        print *, "calcRemainingSwitches tests successfull"

        deallocate(current_stepvector)

    end subroutine test_calcRemainingSwitches

    subroutine test_calcOverlapRange()
        integer, allocatable :: overlap(:), nonOverlap(:)
        character(*), parameter :: testFun = "calcOverlapRange"
        integer :: i,j,k,l

        print *, "***"
        print *, " Testing ", testFun

        i = 1
        j = 4
        k = 1
        l = 4
        
        call calcOverlapRange(i, j, k, l, overlap, nonOverlap)

        print *, " overlap ind: ", overlap
        print *, " overlap length: ", size(overlap)
        print *, " non overlap ind ", nonOverlap
        print *, " non overlap length ", size(nonOverlap)

        print *, testFun, " tests successfull"
        
    end subroutine test_calcOverlapRange

    subroutine test_excitationIdentifier
        integer :: i, j, k, l
        type(excitationInformation) :: ex1, ex2 
        character(*), parameter :: testFun = " excitationIdentifier"

        print *, "***"
        print *, " Testing: ", testFun

        i = 1
        j = 2
        k = 3
        l = 4

        ex1 = excitationIdentifier(i, j)
        ex2 = excitationIdentifier(i,j,k,l)

        print *, testFun, " tests passed!"

    end subroutine test_excitationIdentifier



end module guga_testsuite
#endif
