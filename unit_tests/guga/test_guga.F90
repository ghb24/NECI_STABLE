#include "macros.h"
! GUGA testsuite.
! contains all GUGA related unit tests
! maybe for testing purposes hard-compile a test system with fixed parameters
! like electrons, orbitals etc...
! since otherwise all test cases are quite input dependent...
! discuss with simon how to implement that optimally

program test_guga
    use fruit

    use SystemData
    use guga_bitRepOps
    use guga_bitRepOps, only: global_csf_i => current_csf_i
    use guga_excitations
    use guga_main, only: generate_excitation_guga, &
        createStochasticExcitation_single, createStochasticExcitation_double
    use guga_matrixElements
    use guga_data
    use guga_types
    use guga_init
    use guga_procedure_pointers
    use guga_plugin, only: init_guga_plugin
    use guga_rdm, only: calc_all_excits_guga_rdm_singles, calc_explicit_1_rdm_guga, &
                        calc_explicit_2_rdm_guga, &
                        combine_x0_x1, &
                        pure_rdm_ind, generator_sign, create_all_rdm_contribs, &
                        extract_molcas_1_rdm_index, contract_molcas_1_rdm_index, &
                        extract_molcas_2_rdm_index, contract_molcas_2_rdm_index, &
                        calc_all_excits_guga_rdm_doubles, &
                        conjugate_rdm_ind
    use constants
    use DetBitOps
    use Determinants
    use bit_reps
    use FciMCData
    use dsfmt_interface, only: dsfmt_init
    use util_mod, only: operator(.isclose.), near_zero, operator(.div.), &
                        binary_search, get_free_unit, stop_all, get_unique_filename
    use sort_mod, only: sort
    use rdm_data_utils, only: calc_combined_rdm_label, calc_separate_rdm_labels
    use fruit_extensions, only: my_run_test_case

    better_implicit_none

    real(dp), parameter :: tol = 1.0e-10_dp
    integer :: failed_count

    call init_fruit()
    call dsfmt_init(0)

    call guga_test_driver()

    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop - 1

contains

    subroutine guga_test_driver

        call init_guga_testsuite()
        call run_test_case(compare_rdm_all_excits_and_mat_eles, &
                           "compare_rdm_all_excits_and_mat_eles")

        call test_guga_bitRepOps
        call test_guga_excitations_stochastic
        call test_guga_excitations_exact
        call test_guga_matrixElements
        call test_guga_data
        call test_guga_explicit_rdms()

        call run_test_case(test_excitationIdentifier, "test_excitationIdentifier")
        call run_test_case(test_bitChecks, "test_bitChecks")
        call run_test_case(test_identify_excitation, "test_identify_excitation")
        call run_test_case(test_identify_excitation_and_matrix_element, &
                           "test_identify_excitation_and_matrix_element")

        !TODO maybe run the excit-gen test also!
        !call run_test_excit_gen_guga_S0

    end subroutine guga_test_driver

    subroutine test_contract_extract_1_rdm_molcas

        integer :: i, j, ij, ij_
        print *, ""
        print *, "testing: contract/extract 1-RDM molcas style rdm index"

        ij = contract_molcas_1_rdm_index(1, 1)
        call extract_molcas_1_rdm_index(ij, i, j)
        call assert_equals(1, i)
        call assert_equals(1, j)

        ij = int(contract_1_rdm_ind(1, 2))
        call extract_molcas_1_rdm_index(ij, i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)
        ij_ = contract_molcas_1_rdm_index(2, 1)

        call extract_molcas_1_rdm_index(ij_, i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(ij, ij_)

        print *, ""
        print *, "testing: contract/extract 1-RDM molcas style rdm index. DONE!"
    end subroutine test_contract_extract_1_rdm_molcas

    subroutine test_contract_extract_2_rdm_molcas

        integer :: i, j, k, l, ij, kl, ijkl
        print *, ""
        print *, "testing: contract/extract 2-RDM molcas style rdm index"
        ijkl = contract_molcas_2_rdm_index(1, 1, 1, 1)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(1, l)
        call assert_equals(1, ij)
        call assert_equals(1, kl)

        ijkl = contract_molcas_2_rdm_index(2, 1, 1, 1)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(1, l)
        call assert_equals(2, ij)
        call assert_equals(1, kl)

        ijkl = contract_molcas_2_rdm_index(1, 1, 1, 2)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(1, l)
        call assert_equals(2, ij)
        call assert_equals(1, kl)

        ijkl = contract_molcas_2_rdm_index(1, 2, 2, 1)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)
        call assert_equals(2, ij)
        call assert_equals(2, kl)

        ijkl = contract_molcas_2_rdm_index(2, 1, 2, 1)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)
        call assert_equals(2, ij)
        call assert_equals(2, kl)

        ijkl = contract_molcas_2_rdm_index(2, 1, 1, 2)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)
        call assert_equals(2, ij)
        call assert_equals(2, kl)

        ijkl = contract_molcas_2_rdm_index(2, 1, 2, 1)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)
        call assert_equals(2, ij)
        call assert_equals(2, kl)

        ijkl = contract_molcas_2_rdm_index(5, 1, 4, 5)
        call extract_molcas_2_rdm_index(ijkl, i, j, k, l, ij, kl)

        call assert_equals(5, i)
        call assert_equals(4, j)
        call assert_equals(5, k)
        call assert_equals(1, l)

        print *, ""
        print *, "testing: contract/extract 2-RDM molcas style rdm index. DONE!"

    end subroutine test_contract_extract_2_rdm_molcas

    subroutine test_create_all_rdm_contribs

        integer(int_rdm), allocatable :: rdm_inds(:), rdm_ind_ex(:)
        real(dp), allocatable :: rdm_mats(:), rdm_mat_ex(:)
        integer(int_rdm) :: rdm_ind, rdm_ex
        real(dp) :: x0, x1
        integer(n_int) :: ilut(0:GugaBits%len_tot), t(0:GugaBits%len_tot), &
                          ilutJ(0:GugaBits%len_tot)
        real(dp) :: pgen, mat_ex
        integer :: nI(4), nex, pos, i, j, cnt
        HElement_t(dp) :: mat_ele
        integer(n_int), allocatable :: ex(:, :)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        logical :: compFlag
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, "testing: create_all_rdm_contribs"

        nI = [1, 2, 3, 4]
        call EncodeBitDet_guga(nI, ilut)
        csf_i = CSF_Info_t(ilut)

        ! to test it fully, create a stochastic excitation and then use
        ! the obtained rdm_ind, x0 and x1 and then compare this to the
        ! calc_guga_matrix_element and calc_explicit_1/2_rdm_guga routines!!

        call createStochasticExcitation_single(ilut, nI, csf_i, t, pgen)

        if (pgen > EPS) then
            call extract_stochastic_rdm_info(GugaBits, t, rdm_ind, x0, x1)
            call create_all_rdm_contribs(rdm_ind, x0, x1, rdm_inds, rdm_mats)

            call assert_equals(1, size(rdm_inds))
            call assert_equals(rdm_ind, rdm_inds(1))
            call assert_equals(x0, rdm_mats(1))

            call calc_explicit_1_rdm_guga(ilut, csf_i, nEx, ex)
            do i = 1, nex
                if (DetBitEq(t(0:nifd), ex(0:nifd, i), nifd)) then
                    pos = i
                end if
            end do
            call assert_true(pos > 0)

            ilutJ = ex(:, pos)
            call assert_equals(extract_rdm_ind(ilutJ), pure_rdm_ind(rdm_inds(1)))
            call assert_equals(extract_matrix_element(ilutJ, 1), rdm_mats(1))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilut, csf_i, t, CSF_Info_t(t), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_ex, &
                                          rdm_mat=rdm_mat_ex)
            x0 = extract_stochastic_rdm_x0(GugaBits, t)
            call assert_equals(1, size(rdm_ind_ex))
            call assert_equals(rdm_ind_ex(1), pure_rdm_ind(rdm_inds(1)))
            call assert_equals(rdm_mat_ex(1), rdm_mats(1))

        end if

        nI = [1, 2, 3, 6]
        call EncodeBitDet_guga(nI, ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 2, 3)
        call assert_true(excitInfo%typ == excit_type%double_L_to_R_to_L)
        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call calcDoubleL2R2L_stochastic(ilut, csf_i, excitInfo, t, pgen, posSwitches, negSwitches)
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(t) == [1, 3, 0, 2]))

        call extract_stochastic_rdm_info(GugaBits, t, rdm_ind, x0, x1)
        call create_all_rdm_contribs(rdm_ind, x0, x1, rdm_inds, rdm_mats)

        call assert_true(size(rdm_inds) > 0)
        call assert_equals(rdm_ind, rdm_inds(1))

        call calc_explicit_2_rdm_guga(ilut, csf_i, nEx, ex)

        cnt = 0
        do i = 1, nex
            if (DetBitEq(t(0:nifd), ex(0:nifd, i))) then
                rdm_ex = extract_rdm_ind(ex(:, i))
                mat_ex = extract_matrix_element(ex(:, i), 1)
                do j = 1, size(rdm_inds)
                    if (pure_rdm_ind(rdm_inds(j)) == rdm_ex) then
                        cnt = cnt + 1
                        call assert_equals(rdm_mats(j), mat_ex, 1e-10_dp)
                    end if
                end do
            end if
        end do
        call assert_true(cnt > 0)

        ! also test with matrix element calculator!
        call calc_guga_matrix_element(ilut, csf_i, t, CSF_Info_t(t), excitInfo, &
                                      mat_ele, t_hamil=.true., rdm_ind=rdm_ind_ex, &
                                      rdm_mat=rdm_mat_ex)

        cnt = 0
        do i = 1, size(rdm_inds)
            do j = 1, size(rdm_ind_ex)
                if (pure_rdm_ind(rdm_inds(i)) == rdm_ind_ex(j) .or. &
                    pure_rdm_ind(rdm_inds(i)) == conjugate_rdm_ind(rdm_ind_ex(j), 2)) then
                    call assert_equals(rdm_mat_ex(j), rdm_mats(i))
                    cnt = cnt + 1
                end if
            end do
        end do
        call assert_true(cnt > 0)

        print *, ""
        print *, "testing: create_all_rdm_contribs. DONE"

    end subroutine test_create_all_rdm_contribs

    subroutine test_transfer_stochastic_rdm_info
        integer(n_int) :: ilutG(0:GugaBits%len_tot), ilutG2(0:GugaBits%len_tot)
        integer(n_int) :: ilutN(0:IlutBits%len_tot), ilutN2(0:IlutBits%len_tot)
        integer(n_int) :: ilutP(0:IlutBitsParent%len_tot)
        integer(int_rdm) :: rdm_ind
        real(dp) :: x0, x1

        print *, ""
        print *, "testing: transfer_stochastic_rdm_info"

        call encode_stochastic_rdm_info(GugaBits, ilutG, &
                                        rdm_ind=1_n_int, x0=-1.0_dp, x1=1.0_dp)

        call transfer_stochastic_rdm_info(ilutG, ilutg2, BitIndex_to=GugaBits)
        call extract_stochastic_rdm_info(GugaBits, ilutg2, rdm_ind, x0, x1)
        call assert_equals(1_n_int, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(1.0_dp, x1)

        call transfer_stochastic_rdm_info(ilutG, ilutN)
        call extract_stochastic_rdm_info(IlutBits, ilutN, rdm_ind, x0, x1)
        call assert_equals(1_n_int, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(1.0_dp, x1)

        call transfer_stochastic_rdm_info(ilutN, ilutN2, BitIndex_from=IlutBits)
        call extract_stochastic_rdm_info(IlutBits, ilutN2, rdm_ind, x0, x1)
        call assert_equals(1_n_int, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(1.0_dp, x1)

        call transfer_stochastic_rdm_info(ilutN, ilutP, &
                                          BitIndex_from=IlutBits, BitIndex_to=IlutBitsParent)
        call extract_stochastic_rdm_info(IlutBitsParent, ilutP, rdm_ind, x0, x1)
        call assert_equals(1_n_int, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(1.0_dp, x1)

        print *, ""
        print *, "testing: transfer_stochastic_rdm_info. DONE!"

    end subroutine test_transfer_stochastic_rdm_info

    subroutine test_encode_extract_stochastic_rdm_ind

        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer(n_int) :: ilutN(0:IlutBits%len_tot)
        integer(n_int) :: ilutP(0:IlutBitsParent%len_tot)
        integer(int_rdm) :: rdm_ind

        print *, ""
        print *, "testing: encode and extract stochastic rdm ind"

        call encode_stochastic_rdm_ind(GugaBits, ilut, 2_int_rdm)
        rdm_ind = extract_stochastic_rdm_ind(GugaBits, ilut)
        call assert_equals(2_int_rdm, rdm_ind)

        call encode_stochastic_rdm_ind(GugaBits, ilut, -1_int_rdm)
        rdm_ind = extract_stochastic_rdm_ind(GugaBits, ilut)
        call assert_equals(-1_int_rdm, rdm_ind)

        call encode_stochastic_rdm_ind(IlutBits, IlutN, 2_int_rdm)
        rdm_ind = extract_stochastic_rdm_ind(IlutBits, IlutN)
        call assert_equals(2_int_rdm, rdm_ind)

        call encode_stochastic_rdm_ind(IlutBits, IlutN, -1_int_rdm)
        rdm_ind = extract_stochastic_rdm_ind(IlutBits, IlutN)
        call assert_equals(-1_int_rdm, rdm_ind)

        call encode_stochastic_rdm_ind(IlutBitsParent, IlutP, 2_int_rdm)
        rdm_ind = extract_stochastic_rdm_ind(IlutBitsParent, IlutP)
        call assert_equals(2_int_rdm, rdm_ind)

        call encode_stochastic_rdm_ind(IlutBitsParent, IlutP, -1_int_rdm)
        rdm_ind = extract_stochastic_rdm_ind(IlutBitsParent, IlutP)
        call assert_equals(-1_int_rdm, rdm_ind)

        print *, ""
        print *, "testing: encode and extract stochastic rdm ind. DONE"

    end subroutine test_encode_extract_stochastic_rdm_ind

    subroutine test_encode_extract_stochastic_rdm_x0

        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer(n_int) :: ilutN(0:IlutBits%len_tot)
        integer(n_int) :: ilutP(0:IlutBitsParent%len_tot)
        real(dp) :: x0

        print *, ""
        print *, "testing: encode and extract stochastic rdm x0"

        call encode_stochastic_rdm_x0(GugaBits, ilut, 0.0_dp)
        x0 = extract_stochastic_rdm_x0(GugaBits, ilut)
        call assert_equals(0.0_dp, x0)

        call encode_stochastic_rdm_x0(GugaBits, ilut, 1.0_dp)
        x0 = extract_stochastic_rdm_x0(GugaBits, ilut)
        call assert_equals(1.0_dp, x0)

        call encode_stochastic_rdm_x0(GugaBits, ilut, -1.0_dp)
        x0 = extract_stochastic_rdm_x0(GugaBits, ilut)
        call assert_equals(-1.0_dp, x0)

        call encode_stochastic_rdm_x0(IlutBits, ilutN, 0.0_dp)
        x0 = extract_stochastic_rdm_x0(IlutBits, ilutN)
        call assert_equals(0.0_dp, x0)

        call encode_stochastic_rdm_x0(IlutBits, ilutN, 1.0_dp)
        x0 = extract_stochastic_rdm_x0(IlutBits, ilutN)
        call assert_equals(1.0_dp, x0)

        call encode_stochastic_rdm_x0(IlutBits, ilutN, -1.0_dp)
        x0 = extract_stochastic_rdm_x0(IlutBits, ilutN)
        call assert_equals(-1.0_dp, x0)

        call encode_stochastic_rdm_x0(IlutBitsParent, ilutP, 0.0_dp)
        x0 = extract_stochastic_rdm_x0(IlutBitsParent, ilutP)
        call assert_equals(0.0_dp, x0)

        call encode_stochastic_rdm_x0(IlutBitsParent, ilutP, 1.0_dp)
        x0 = extract_stochastic_rdm_x0(IlutBitsParent, ilutP)
        call assert_equals(1.0_dp, x0)

        call encode_stochastic_rdm_x0(IlutBitsParent, ilutP, -1.0_dp)
        x0 = extract_stochastic_rdm_x0(IlutBitsParent, ilutP)
        call assert_equals(-1.0_dp, x0)

        print *, ""
        print *, "testing: encode and extract stochastic rdm x0. DONE!"

    end subroutine test_encode_extract_stochastic_rdm_x0

    subroutine test_encode_extract_stochastic_rdm_x1

        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer(n_int) :: ilutN(0:IlutBits%len_tot)
        integer(n_int) :: ilutP(0:IlutBitsParent%len_tot)
        real(dp) :: x1
        print *, ""
        print *, "testing: encode and extract stochastic rmd x1"

        call encode_stochastic_rdm_x1(GugaBits, ilut, 0.0_dp)
        x1 = extract_stochastic_rdm_x1(GugaBits, ilut)
        call assert_equals(0.0_dp, x1)

        call encode_stochastic_rdm_x1(GugaBits, ilut, 1.0_dp)
        x1 = extract_stochastic_rdm_x1(GugaBits, ilut)
        call assert_equals(1.0_dp, x1)

        call encode_stochastic_rdm_x1(GugaBits, ilut, -1.0_dp)
        x1 = extract_stochastic_rdm_x1(GugaBits, ilut)
        call assert_equals(-1.0_dp, x1)

        call encode_stochastic_rdm_x1(IlutBits, ilutN, 0.0_dp)
        x1 = extract_stochastic_rdm_x1(IlutBits, ilutN)
        call assert_equals(0.0_dp, x1)

        call encode_stochastic_rdm_x1(IlutBits, ilutN, 1.0_dp)
        x1 = extract_stochastic_rdm_x1(IlutBits, ilutN)
        call assert_equals(1.0_dp, x1)

        call encode_stochastic_rdm_x1(IlutBits, ilutN, -1.0_dp)
        x1 = extract_stochastic_rdm_x1(IlutBits, ilutN)
        call assert_equals(-1.0_dp, x1)

        call encode_stochastic_rdm_x1(IlutBitsParent, ilutP, 0.0_dp)
        x1 = extract_stochastic_rdm_x1(IlutBitsParent, ilutP)
        call assert_equals(0.0_dp, x1)

        call encode_stochastic_rdm_x1(IlutBitsParent, ilutP, 1.0_dp)
        x1 = extract_stochastic_rdm_x1(IlutBitsParent, ilutP)
        call assert_equals(1.0_dp, x1)

        call encode_stochastic_rdm_x1(IlutBitsParent, ilutP, -1.0_dp)
        x1 = extract_stochastic_rdm_x1(IlutBitsParent, ilutP)
        call assert_equals(-1.0_dp, x1)

        print *, ""
        print *, "testing: encode and extract stochastic rmd x1. DONE!"

    end subroutine test_encode_extract_stochastic_rdm_x1

    subroutine test_encode_extract_stochastic_rdm_info

        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer(n_int) :: ilutN(0:IlutBits%len_tot)
        integer(n_int) :: ilutP(0:IlutBitsParent%len_tot)
        integer(int_rdm) :: rdm_ind
        real(dp) :: x0, x1

        print *, ""
        print *, "testing: encode and extract stochastic rdm info"

        call encode_stochastic_rdm_info(GugaBits, ilut, 0_int_rdm, 0.0_dp, 0.0_dp)
        call extract_stochastic_rdm_info(GugaBits, ilut, rdm_ind, x0, x1)
        call assert_equals(0_int_rdm, rdm_ind)
        call assert_equals(0.0_dp, x0)
        call assert_equals(0.0_dp, x1)

        call encode_stochastic_rdm_info(GugaBits, ilut, 1_int_rdm, -1.0_dp, 10.0_dp)
        call extract_stochastic_rdm_info(GugaBits, ilut, rdm_ind, x0, x1)
        call assert_equals(1_int_rdm, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(10.0_dp, x1)

        call encode_stochastic_rdm_info(IlutBits, ilutN, 0_int_rdm, 0.0_dp, 0.0_dp)
        call extract_stochastic_rdm_info(IlutBits, ilutN, rdm_ind, x0, x1)
        call assert_equals(0_int_rdm, rdm_ind)
        call assert_equals(0.0_dp, x0)
        call assert_equals(0.0_dp, x1)

        call encode_stochastic_rdm_info(IlutBits, ilutN, 1_int_rdm, -1.0_dp, 10.0_dp)
        call extract_stochastic_rdm_info(IlutBits, ilutN, rdm_ind, x0, x1)
        call assert_equals(1_int_rdm, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(10.0_dp, x1)

        call encode_stochastic_rdm_info(IlutBitsParent, ilutP, 0_int_rdm, 0.0_dp, 0.0_dp)
        call extract_stochastic_rdm_info(IlutBitsParent, ilutP, rdm_ind, x0, x1)
        call assert_equals(0_int_rdm, rdm_ind)
        call assert_equals(0.0_dp, x0)
        call assert_equals(0.0_dp, x1)

        call encode_stochastic_rdm_info(IlutBitsParent, ilutP, 1_int_rdm, -1.0_dp, 10.0_dp)
        call extract_stochastic_rdm_info(IlutBitsParent, ilutP, rdm_ind, x0, x1)
        call assert_equals(1_int_rdm, rdm_ind)
        call assert_equals(-1.0_dp, x0)
        call assert_equals(10.0_dp, x1)

        print *, ""
        print *, "testing: encode and extract stochastic rdm info. DONE"

    end subroutine test_encode_extract_stochastic_rdm_info

    subroutine compare_rdm_all_excits_and_mat_eles

        integer, allocatable :: nI(:)
        integer(n_int) :: ilutI(0:nifguga), ilutJ(0:nifguga)
        integer :: n_tot, n, m
        integer(n_int), allocatable :: excits(:, :)
        integer(int_rdm) :: rdm_ind_1
        real(dp) :: rdm_mat_1
        type(ExcitationInformation_t) :: excitInfo
        integer(int_rdm), allocatable :: rdm_ind(:)
        real(dp), allocatable :: rdm_mat(:)
        HElement_t(dp) :: mat_ele
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, "comparing coupling coeffs from exact and from calc_guga_matrix_element"

        nel = 2
        allocate (nI(nel))
        nI = [1, 2]
        call EncodeBitDet_guga(nI, ilutI)
        csf_i = CSF_Info_t(ilutI)

        call calc_explicit_2_rdm_guga(ilutI, csf_i, n_tot, excits)

        do n = 1, n_tot
            ilutJ = excits(:, n)
            rdm_ind_1 = extract_rdm_ind(ilutJ)
            rdm_mat_1 = real(extract_h_element(ilutJ), dp)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                          t_hamil=.false., rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            call assert_true(any(pure_rdm_ind(rdm_ind) == rdm_ind_1) .or. &
                             any(pure_rdm_ind(conjugate_rdm_ind(rdm_ind, 2)) == rdm_ind_1))

            do m = 1, size(rdm_ind)
                if (rdm_ind(m) == rdm_ind_1 .or. &
                    conjugate_rdm_ind(rdm_ind(m), 2) == rdm_ind_1) then
                    call assert_equals(rdm_mat(m), rdm_mat_1)
                end if
            end do
        end do
        deallocate (nI)

        allocate (nI(nel))
        nI = [1, 4]

        call EncodeBitDet_guga(nI, ilutI)
        csf_i = CSF_Info_t(ilutI)
        call calc_explicit_2_rdm_guga(ilutI, csf_i, n_tot, excits)

        do n = 1, n_tot
            ilutJ = excits(:, n)
            rdm_ind_1 = extract_rdm_ind(ilutJ)
            rdm_mat_1 = real(extract_h_element(ilutJ), dp)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                          t_hamil=.false., rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            call assert_true(any(pure_rdm_ind(rdm_ind) == rdm_ind_1) .or. &
                             any(pure_rdm_ind(conjugate_rdm_ind(rdm_ind, 2)) == rdm_ind_1))

            do m = 1, size(rdm_ind)
                if (rdm_ind(m) == rdm_ind_1 .or. &
                    conjugate_rdm_ind(rdm_ind(m), 2) == rdm_ind_1) then
                    call assert_equals(rdm_mat(m), rdm_mat_1, 1e-12_dp)
                end if
            end do
        end do
        deallocate (nI)

        nel = 4
        allocate (nI(nel))

        nI = [3, 4, 5, 6]

        call EncodeBitDet_guga(nI, ilutI)
        csf_i = CSF_Info_t(ilutI)
        call calc_explicit_2_rdm_guga(ilutI, csf_i, n_tot, excits)

        do n = 1, n_tot
            ilutJ = excits(:, n)
            rdm_ind_1 = extract_rdm_ind(ilutJ)
            rdm_mat_1 = real(extract_h_element(ilutJ), dp)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                          t_hamil=.false., rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            call assert_true(any(pure_rdm_ind(rdm_ind) == rdm_ind_1) .or. &
                             any(pure_rdm_ind(conjugate_rdm_ind(rdm_ind, 2)) == rdm_ind_1))

            do m = 1, size(rdm_ind)
                if (rdm_ind(m) == rdm_ind_1 .or. &
                    conjugate_rdm_ind(rdm_ind(m), 2) == rdm_ind_1) then
                    call assert_equals(rdm_mat(m), rdm_mat_1, 1e-12_dp)
                end if
            end do
        end do

        nI = [1, 2, 7, 8]

        call EncodeBitDet_guga(nI, ilutI)
        csf_i = CSF_Info_t(ilutI)
        call calc_explicit_2_rdm_guga(ilutI, csf_i, n_tot, excits)

        do n = 1, n_tot
            ilutJ = excits(:, n)
            rdm_ind_1 = extract_rdm_ind(ilutJ)
            rdm_mat_1 = real(extract_h_element(ilutJ), dp)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                          t_hamil=.false., rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            call assert_true(any(pure_rdm_ind(rdm_ind) == rdm_ind_1) .or. &
                             any(pure_rdm_ind(conjugate_rdm_ind(rdm_ind, 2)) == rdm_ind_1))

            do m = 1, size(rdm_ind)
                if (rdm_ind(m) == rdm_ind_1 .or. &
                    conjugate_rdm_ind(rdm_ind(m), 2) == rdm_ind_1) then
                    call assert_equals(rdm_mat(m), rdm_mat_1, 1e-12_dp)
                end if
            end do
        end do

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilutI)
        csf_i = CSF_Info_t(ilutI)
        call calc_explicit_2_rdm_guga(ilutI, csf_i, n_tot, excits)

        do n = 1, n_tot
            ilutJ = excits(:, n)
            rdm_ind_1 = extract_rdm_ind(ilutJ)
            rdm_mat_1 = real(extract_h_element(ilutJ), dp)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                          t_hamil=.false., rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            call assert_true(any(pure_rdm_ind(rdm_ind) == rdm_ind_1) .or. &
                             any(pure_rdm_ind(conjugate_rdm_ind(rdm_ind, 2)) == rdm_ind_1))

            do m = 1, size(rdm_ind)
                if (rdm_ind(m) == rdm_ind_1 .or. &
                    conjugate_rdm_ind(rdm_ind(m), 2) == rdm_ind_1) then
                    call assert_equals(rdm_mat(m), rdm_mat_1, 1e-12_dp)
                end if
            end do
        end do

        nI = [1, 3, 6, 8]
        call EncodeBitDet_guga(nI, ilutI)
        csf_i = CSF_Info_t(ilutI)
        call calc_explicit_2_rdm_guga(ilutI, csf_i, n_tot, excits)

        do n = 1, n_tot
            ilutJ = excits(:, n)
            rdm_ind_1 = extract_rdm_ind(ilutJ)
            rdm_mat_1 = real(extract_h_element(ilutJ), dp)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                          t_hamil=.false., rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            call assert_true(any(pure_rdm_ind(rdm_ind) == rdm_ind_1) .or. &
                             any(pure_rdm_ind(conjugate_rdm_ind(rdm_ind, 2)) == rdm_ind_1))

            do m = 1, size(rdm_ind)
                if (rdm_ind(m) == rdm_ind_1 .or. &
                    conjugate_rdm_ind(rdm_ind(m), 2) == rdm_ind_1) then
                    call assert_equals(rdm_mat(m), rdm_mat_1, 1e-12_dp)
                end if
            end do
        end do

        print *, ""
        print *, "comparing coupling coeffs from exact and from calc_guga_matrix_element. DONE"

    end subroutine compare_rdm_all_excits_and_mat_eles

    subroutine init_guga_testsuite
        call init_guga_plugin(t_testmode_=.true., nel_=4, nbasis_=8, &
                              nSpatOrbs_=4)
    end subroutine init_guga_testsuite

    subroutine test_guga_explicit_rdms
        print *, ""
        print *, "testing explicit GUGA RDM routines"
        print *, ""

        call run_test_case(test_calc_all_excits_guga_rdm_singles, &
                           "test_calc_all_excits_guga_rdm_singles")
        call run_test_case(test_calc_explicit_1_rdm_guga, "test_calc_explicit_1_rdm_guga")

        call run_test_case(test_calc_all_excits_guga_rdm_doubles, &
                           "test_calc_all_excits_guga_rdm_doubles")
        call run_test_case(test_calc_explicit_2_rdm_guga, "test_calc_explicit_2_rdm_guga")

        call run_test_case(test_compare_RDM_indexing, "test_compare_RDM_indexing")

        call run_test_case(test_create_all_rdm_contribs, &
                           "test_create_all_rdm_contribs")
        call run_test_case(test_generator_sign, "test_generator_sign")

        print *, ""
        print *, "explicit RDM routines passed!"
        print *, ""
    end subroutine test_guga_explicit_rdms

    subroutine test_generator_sign

        print *, ""
        print *, "testing: generator_sign"

        call assert_equals(1.0_dp, generator_sign(0, 0, 0, 0))
        call assert_equals(1.0_dp, generator_sign(1, 1, 1, 1))
        call assert_equals(1.0_dp, generator_sign(1, 1, 2, 2))
        call assert_equals(1.0_dp, generator_sign(1, 2, 3, 4))
        call assert_equals(1.0_dp, generator_sign(1, 3, 2, 4))
        call assert_equals(1.0_dp, generator_sign(4, 2, 1, 3))
        call assert_equals(-1.0_dp, generator_sign(2, 3, 1, 4))
        call assert_equals(-1.0_dp, generator_sign(1, 4, 2, 3))
        call assert_equals(1.0_dp, generator_sign(4, 1, 2, 3))
        call assert_equals(-1.0_dp, generator_sign(4, 1, 3, 2))
        print *, ""
        print *, "testing: generator_sign. DONE"

    end subroutine test_generator_sign

    subroutine test_compare_RDM_indexing

        integer(int_rdm) :: ijkl, abcd, ab, cd
        integer :: i, j, k, l, a, b, c, d, ij, kl

        print *, ""
        print *, " compare 'old' SD-based RDM indexing and GUGA convention"
        print *, ""

        call calc_combined_rdm_label(1, 1, 1, 1, ijkl)
        abcd = contract_2_rdm_ind(1, 1, 1, 1)

        call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
        call extract_2_rdm_ind(abcd, a, b, c, d, ab, cd)

        call assert_equals(i, a)
        call assert_equals(j, b)
        call assert_equals(k, c)
        call assert_equals(l, d)

        call calc_combined_rdm_label(1, 2, 3, 4, ijkl)
        abcd = contract_2_rdm_ind(1, 2, 3, 4)

        call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
        call extract_2_rdm_ind(abcd, a, b, c, d, ab, cd)

        call assert_equals(i, a)
        call assert_equals(j, b)
        call assert_equals(k, c)
        call assert_equals(l, d)

        call calc_combined_rdm_label(3, 2, 1, 4, ijkl)
        abcd = contract_2_rdm_ind(3, 2, 1, 4)

        call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
        call extract_2_rdm_ind(abcd, a, b, c, d, ab, cd)

        call assert_equals(i, a)
        call assert_equals(j, b)
        call assert_equals(k, c)
        call assert_equals(l, d)

        print *, ""
        print *, " compare 'old' SD-based RDM indexing and GUGA convention DONE"
        print *, ""

    end subroutine test_compare_RDM_indexing

    subroutine test_calc_explicit_2_rdm_guga

        integer :: n_tot
        integer(n_int), allocatable :: excits(:, :)
        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, cnt, n
        integer(n_int) :: ilut(0:nifguga)
        type(CSF_Info_t) :: csf_i

        nel = 2
        call EncodeBitDet_guga([5, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        print *, ""
        print *, "testing: calc_explicit_2_rdm_guga"
        print *, ""

        call calc_explicit_2_rdm_guga(ilut, csf_i, n_tot, excits)

        call assert_equals(9, n_tot)

        !  1 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(2.0_dp, real(extract_h_element(excits(:, 1)), dp), 1e-12_dp)

        ! 1 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:, 2))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:, 2)), dp), 1e-12_dp)

        ! 1 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:, 3))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:, 3)), dp), 1e-12_dp)

        ! 2 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:, 4))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:, 4)), dp), 1e-12_dp)

        ! 2 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:, 5))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(2.0_dp, real(extract_h_element(excits(:, 5)), dp), 1e-12_dp)

        ! 2 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:, 6))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:, 6)), dp), 1e-12_dp)

        ! 4 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:, 7))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:, 7)), dp), 1e-12_dp)

        ! 4 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:, 8))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:, 8)), dp), 1e-12_dp)

        ! 4 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:, 9))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(2.0_dp, real(extract_h_element(excits(:, 9)), dp), 1e-12_dp)

        nel = 4

        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_explicit_2_rdm_guga(ilut, csf_i, n_tot, excits)

        cnt = 0
        do n = 1, n_tot
            if (DetBitEQ(ilut, excits(:, n))) cnt = cnt + 1
        end do
        call assert_equals(0, cnt)

        print *, ""
        print *, "testing: calc_explicit_2_rdm_guga DONE"
        print *, ""

        nel = 4

    end subroutine test_calc_explicit_2_rdm_guga

    subroutine test_calc_all_excits_guga_rdm_doubles

        integer(n_int) :: ilut(0:nifguga)
        integer :: n_excits, i, j, k, l
        integer(n_int), allocatable :: excits(:, :)
        integer, allocatable :: nJ(:)
        real(dp) :: rdm_mat
        integer(int_rdm) :: rdm_ind
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, "testing: calc_all_excits_guga_rdm_doubles"
        print *, ""

        nel = 2
        call EncodeBitDet_guga([1, 2], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 2, 1, 0, 0, excits, n_excits)

        allocate (nJ(2))

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([1, 4], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 2, 1, 2, 1, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([3, 4], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(2.0_dp, rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i=i, j=j, k=k, l=l)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 2, 1, 3, 1, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([3, 6], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i=i, j=j, k=k, l=l)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(3, k)
        call assert_equals(1, l)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 3, 1, 2, 1, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([3, 6], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i=i, j=j, k=k, l=l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)

        call EncodeBitDet_guga([3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 1, 2, 1, 2, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([1, 2], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(2.0_dp, rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i=i, j=j, k=k, l=l)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, k)
        call assert_equals(2, l)

        call EncodeBitDet_guga([5, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 1, 3, 2, 3, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([1, 4], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i=i, j=j, k=k, l=l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 1, 3, 4, 3, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([1, 8], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:, 1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i=i, j=j, k=k, l=l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)

        nel = 3

        call EncodeBitDet_guga([1, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 3, 1, 2, 3, excits, n_excits)
        call assert_equals(0, n_excits)

        nel = 4
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        csf_i = CSF_Info_t(ilut)

        deallocate (nJ)
        allocate (nJ(nel))

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 1, 4, 4, 1, excits, n_excits)

        call assert_equals(2, n_excits)

        call decode_bit_det(nJ, excits(:, 1))
        call assert_equals([1, 3, 6, 8], nJ, 4)
        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(1, l)

        call decode_bit_det(nJ, excits(:, 2))
        call assert_equals([1, 4, 5, 8], nJ, 4)
        rdm_ind = extract_rdm_ind(excits(:, 2))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(1, l)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 4, 1, 1, 4, excits, n_excits)

        call assert_equals(2, n_excits)

        call decode_bit_det(nJ, excits(:, 1))
        call assert_equals([1, 3, 6, 8], nJ, 4)
        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(4, l)

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, csf_i, 1, 4, 4, 1, excits, n_excits)

        call decode_bit_det(nJ, excits(:, 1))
        call assert_equals([1, 4, 5, 8], nJ, 4)
        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(1, l)

        call decode_bit_det(nJ, excits(:, 2))
        call assert_equals([1, 3, 6, 8], nJ, 4)
        rdm_ind = extract_rdm_ind(excits(:, 2))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(1, l)

        print *, ""
        print *, "testing: calc_all_excits_guga_rdm_doubles DONE"
        print *, ""

        nel = 4

    end subroutine test_calc_all_excits_guga_rdm_doubles

    subroutine test_calc_explicit_1_rdm_guga

        integer :: n_tot, i, iex, j
        real(dp) :: rdm_mat
        integer(int_rdm) :: rdm_ind
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), allocatable :: excits(:, :)
        integer, allocatable :: nJ(:)
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, "testing: calc_explicit_1_rdm_guga"
        print *, ""

        nel = 2
        call EncodeBitDet_guga([1, 2], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_explicit_1_rdm_guga(ilut, csf_i, n_tot, excits)
        call assert_equals(3, n_tot)

        allocate (nJ(2))

        do iEx = 1, 3

            rdm_mat = real(extract_h_element(excits(:, iEx)), dp)
            call assert_equals(sqrt(2.0_dp), rdm_mat)

            rdm_ind = extract_rdm_ind(excits(:, iEx))
            call extract_1_rdm_ind(rdm_ind, i, j)

            call assert_equals(j, 1)
            call assert_equals(i, iex + 1)

            call decode_bit_det(nJ, excits(:, iEx))

            call assert_equals([1, 2 * (iEx + 1)], nJ, 2)
        end do

        call EncodeBitDet_guga([7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_explicit_1_rdm_guga(ilut, csf_i, n_tot, excits)
        call assert_equals(3, n_tot)

        do iEx = 1, 3

            rdm_mat = real(extract_h_element(excits(:, iEx)), dp)
            call assert_equals(sqrt(2.0_dp), rdm_mat)

            rdm_ind = extract_rdm_ind(excits(:, iEx))
            call extract_1_rdm_ind(rdm_ind, i, j)

            call assert_equals(j, 4)
            call assert_equals(i, iex)

            call decode_bit_det(nJ, excits(:, iEx))

            call assert_equals([2 * iEx - 1, 8], nJ, 2)
        end do

        print *, ""
        print *, "testing: calc_explicit_1_rdm_guga DONE"
        print *, ""

        nel = 4
    end subroutine test_calc_explicit_1_rdm_guga

    subroutine test_calc_all_excits_guga_rdm_singles

        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), allocatable :: excits(:, :)
        integer :: n_excits, i, j
        integer, allocatable :: nJ(:)
        real(dp) :: rdm_mat
        integer(int_rdm) :: rdm_ind
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, "testing: calc_all_excits_guga_rdm_singles"
        print *, ""

        nel = 2
        call EncodeBitDet_guga([1, 2], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_singles(ilut, csf_i, 1, 2, excits, n_excits)

        call assert_equals(0, n_excits)

        call calc_all_excits_guga_rdm_singles(ilut, csf_i, 2, 1, excits, n_excits)
        call assert_equals(1, n_excits)

        allocate (nJ(nel), source=0)

        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([1, 4], nJ, 2)

        rdm_mat = real(extract_h_element(excits(:, 1)), dp)

        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, j)
        call assert_equals(2, i)

        call EncodeBitDet_guga([3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_singles(ilut, csf_i, 2, 1, excits, n_excits)
        call assert_equals(0, n_excits)

        call calc_all_excits_guga_rdm_singles(ilut, csf_i, 1, 2, excits, n_excits)
        call assert_equals(1, n_excits)

        call decode_bit_det(nJ, excits(:, 1))

        call assert_equals([1, 4], nJ, 2)

        rdm_mat = real(extract_h_element(excits(:, 1)), dp)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(2, j)
        call assert_equals(1, i)

        call assert_equals(sqrt(2.0_dp), rdm_mat)

        nel = 3
        deallocate (nJ);
        call EncodeBitDet_guga([1, 2, 3], ilut)
        csf_i = CSF_Info_t(ilut)

        call calc_all_excits_guga_rdm_singles(ilut, csf_i, 3, 1, excits, n_excits)

        call assert_equals(2, n_excits)
        allocate (nJ(nel), source=0)
        call decode_bit_det(nJ, excits(:, 2))
        call assert_equals([1, 4, 5], nJ, 3)

        rdm_ind = extract_rdm_ind(excits(:, 2))
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, j)
        call assert_equals(3, i)

        call decode_bit_det(nJ, excits(:, 1))
        call assert_equals([1, 3, 6], nJ, 3)

        rdm_ind = extract_rdm_ind(excits(:, 1))
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, j)
        call assert_equals(3, i)

        print *, ""
        print *, " calc_all_excits_guga_rdm_singles passed!"
        print *, ""

        nel = 4

    end subroutine test_calc_all_excits_guga_rdm_singles

    subroutine test_guga_bitRepOps

        print *, ""
        print *, "testing functions from module: guga_bitRepOps"
        print *, ""
        call run_test_case(test_findSwitches, "test_findSwitches")
        call run_test_case(test_count_alpha_orbs_ij, "test_count_alpha_orbs_ij")
        call run_test_case(test_count_beta_orbs_ij, "test_count_beta_orbs_ij")
        call run_test_case(test_matrix_element_ops, "test_matrix_element_ops")
        call run_test_case(test_count_open_orbs_ij, "test_count_open_orbs_ij")
        call run_test_case(test_set_get_DeltaB, "test_set_get_DeltaB")
        call run_test_case(test_isProperCSF_ilut, "test_isProperCSF_ilut")
        call run_test_case(test_calcbvector, "test_calcbvector")
        call run_test_case(test_isDouble, "test_isDouble")
        call run_test_case(test_calcStepVector, "test_calcStepVector")
        call run_test_case(test_getSpatialOccupation, "test_getSpatialOccupation")
        call run_test_case(test_calcOcc_vector_ilut, "test_calcOcc_vector_ilut")
        call run_test_case(test_contract_extract_1_rdm, "test_contract_extract_1_rdm")
        call run_test_case(test_contract_extract_2_rdm, "test_contract_extract_2_rdm")
        call run_test_case(test_contract_extract_1_rdm_with_excitInfo, &
                           "test_contract_extract_1_rdm_with_excitInfo")
        call run_test_case(test_contract_extract_2_rdm_with_excitInfo, &
                           "test_contract_extract_2_rdm_with_excitInfo")
        call run_test_case(test_encode_extract_stochastic_rdm_ind, &
                           "test_encode_extract_stochastic_rdm_ind")
        call run_test_case(test_encode_extract_stochastic_rdm_x0, &
                           "test_encode_extract_stochastic_rdm_x0")
        call run_test_case(test_encode_extract_stochastic_rdm_x1, &
                           "test_encode_extract_stochastic_rdm_x1")
        call run_test_case(test_encode_extract_stochastic_rdm_info, &
                           "test_encode_extract_stochastic_rdm_info")
        call run_test_case(test_transfer_stochastic_rdm_info, &
                           "test_transfer_stochastic_rdm_info")

        call run_test_case(test_pure_rdm_ind, "test_pure_rdm_ind")
        call run_test_case(test_contract_extract_1_rdm_molcas, &
                           "test_contract_extract_1_rdm_molcas")
        call run_test_case(test_contract_extract_2_rdm_molcas, &
                           "test_contract_extract_2_rdm_molcas")

        call my_run_test_case(encode_excit_info_type_test, &
                              "encode_excit_info_type_test", "encode_excit_info_type")

        call my_run_test_case(extract_excit_info_type_test, &
                              "extract_excit_info_type_test", "extract_excit_info_type")

        call my_run_test_case(encode_extract_excit_info_indices, &
                              "encode_extract_excit_info_indices", &
                              "encode_excit_info_indices() and extract_excit_info_indices")

        call my_run_test_case(encode_and_extract_excit_info_test, &
            "encode_and_extract_excit_info_test", &
            "encode_excit_info_scalar(), encode_excit_info_vec(), &
            & encode_excit_info_type(), extract_excit_info_obj(), &
            & extract_excit_info_scalar() and extract_excit_info_vec")

        print *, ""
        print *, "guga_bitRepOps tests passed!"
        print *, ""

    end subroutine test_guga_bitRepOps

    subroutine encode_extract_excit_info_indices
        integer(int64) :: excit_info_int
        integer :: a, i, b, j, vec(4)

        call encode_excit_info_indices(excit_info_int, 1, 2, 3, 4)
        call extract_excit_info_indices(excit_info_int, a, i, b, j)
        call extract_excit_info_indices(excit_info_int, vec)
        call assert_equals(1, a)
        call assert_equals(2, i)
        call assert_equals(3, b)
        call assert_equals(4, j)
        call assert_equals([1, 2, 3, 4], vec, 4)

        call encode_excit_info_indices(excit_info_int, [4, 3, 2, 1])
        call extract_excit_info_indices(excit_info_int, a, i, b, j)
        call extract_excit_info_indices(excit_info_int, vec)
        call assert_equals(4, a)
        call assert_equals(3, i)
        call assert_equals(2, b)
        call assert_equals(1, j)
        call assert_equals([4, 3, 2, 1], vec, 4)

    end subroutine encode_extract_excit_info_indices

    subroutine encode_excit_info_type_test
        integer(int64) :: excit_info_int

        excit_info_int = 0_int64

        call encode_excit_info_type(excit_info_int, excit_type%single_overlap_R_to_L)
        call assert_equals(excit_type%single_overlap_R_to_L, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%double_lowering)
        call assert_equals(excit_type%double_lowering, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%double_R_to_L)
        call assert_equals(excit_type%double_R_to_L, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstop_R_to_L)
        call assert_equals(excit_type%fullstop_R_to_L, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstart_R_to_L)
        call assert_equals(excit_type%fullstart_R_to_L, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstop_raising)
        call assert_equals(excit_type%fullstop_raising, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstart_stop_mixed)
        call assert_equals(excit_type%fullstart_stop_mixed, int(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%single_overlap_L_to_R)
        call assert_equals(excit_type%single_overlap_L_to_R, int(excit_info_int))

    end subroutine encode_excit_info_type_test

    subroutine extract_excit_info_type_test
        integer(int64) :: excit_info_int

        excit_info_int = 0_int64
        call encode_excit_info_type(excit_info_int, excit_type%single_overlap_R_to_L)
        call assert_equals(excit_type%single_overlap_R_to_L, &
                           extract_excit_info_type(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstop_R_to_L)
        call assert_equals(excit_type%fullstop_R_to_L, &
                           extract_excit_info_type(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstart_R_to_L)
        call assert_equals(excit_type%fullstart_R_to_L, &
                           extract_excit_info_type(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%double_lowering)
        call assert_equals(excit_type%double_lowering, &
                           extract_excit_info_type(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%double_R_to_L)
        call assert_equals(excit_type%double_R_to_L, &
                           extract_excit_info_type(excit_info_int))

        call encode_excit_info_type(excit_info_int, excit_type%fullstart_stop_alike)
        call assert_equals(excit_type%fullstart_stop_alike, &
                           extract_excit_info_type(excit_info_int))

    end subroutine extract_excit_info_type_test

    subroutine encode_and_extract_excit_info_test
        integer(int64) :: excit_info_int
        integer :: typ, a, i, b, j, inds(4)
        type(ExcitationInformation_t) :: excitInfo

        excit_info_int = encode_excit_info(excit_type%double_R_to_L_to_R, &
                                           1, 4, 3, 2)

        call extract_excit_info(excit_info_int, typ, a, i, b, j)
        call assert_equals(excit_type%double_R_to_L_to_R, typ)
        call extract_excit_info(excit_info_int, typ, inds)
        call assert_equals(excit_type%double_R_to_L_to_R, typ)
        call extract_excit_info(excit_info_int, excitInfo)
        call assert_equals(1, a)
        call assert_equals(4, i)
        call assert_equals(3, b)
        call assert_equals(2, j)
        call assert_equals([1, 4, 3, 2], inds, 4)
        call assert_equals(excit_type%double_R_to_L_to_R, excitInfo%typ)
        call assert_equals(gen_type%R, excitInfo%gen1)
        call assert_equals(gen_type%L, excitInfo%gen2)

        excit_info_int = encode_excit_info(excit_type%double_L_to_R_to_L, &
                                           [2, 3, 4, 1])

        call extract_excit_info(excit_info_int, typ, a, i, b, j)
        call assert_equals(excit_type%double_L_to_R_to_L, typ)
        call extract_excit_info(excit_info_int, typ, inds)
        call assert_equals(excit_type%double_L_to_R_to_L, typ)
        call extract_excit_info(excit_info_int, excitInfo)
        call assert_equals(2, a)
        call assert_equals(3, i)
        call assert_equals(4, b)
        call assert_equals(1, j)
        call assert_equals([2, 3, 4, 1], inds, 4)
        call assert_equals(excit_type%double_L_to_R_to_L, excitInfo%typ)
        call assert_equals(gen_type%R, excitInfo%gen1)
        call assert_equals(gen_type%L, excitInfo%gen2)

        excitInfo%typ = excit_type%fullstart_stop_mixed
        excitInfo%i = 1
        excitInfo%j = 4
        excitInfo%k = 4
        excitInfo%l = 1

        excit_info_int = encode_excit_info(excitInfo)

        call extract_excit_info(excit_info_int, typ, a, i, b, j)
        call assert_equals(excit_type%fullstart_stop_mixed, typ)
        call extract_excit_info(excit_info_int, typ, inds)
        call assert_equals(excit_type%fullstart_stop_mixed, typ)
        call extract_excit_info(excit_info_int, excitInfo)
        call assert_equals(1, a)
        call assert_equals(4, i)
        call assert_equals(4, b)
        call assert_equals(1, j)
        call assert_equals([1, 4, 4, 1], inds, 4)
        call assert_equals(excit_type%fullstart_stop_mixed, excitInfo%typ)
        call assert_equals(gen_type%R, excitInfo%gen1)
        call assert_equals(gen_type%L, excitInfo%gen2)

    end subroutine encode_and_extract_excit_info_test

    subroutine test_guga_excitations_stochastic

        print *, ""
        print *, "testing module: guga_excitations stochastic:"
        print *, ""
        call run_test_case(test_calcMixedContribution, "test_calcMixedContribution")
        call run_test_case(test_pickRandomOrb, "test_pickRandomOrb")
        call run_test_case(test_generate_excitation_guga_double, "test_generate_excitation_guga_double")
        call run_test_case(test_generate_excitation_guga_single, "test_generate_excitation_guga_single")
        call run_test_case(test_pickOrbitals_single, "test_pickOrbitals_single")

        call run_test_case(test_createStochasticStart_single, "test_createStochasticStart_single")
        call run_test_case(test_singleStochasticUpdate, "test_singleStochasticUpdate")
        call run_test_case(test_singleStochasticEnd, "test_singleStochasticEnd")
        call run_test_case(test_createStochasticExcitation_single, "test_createStochasticExcitation_single")

        call run_test_case(test_pickOrbitals_double, "test_pickOrbitals_double")
        ! now go through each specific stochastix excitation calulator:
        call run_test_case(test_mixedFullStartStochastic, "test_mixedFullStartStochastic")
        call run_test_case(test_doubleUpdateStochastic, "test_doubleUpdateStochastic")
        call run_test_case(test_mixedFullStopStochastic, "test_mixedFullStopStochastic")
        call run_test_case(test_calcFullStartFullStopMixedStochastic, "test_calcFullStartFullStopMixedStochastic")
        call run_test_case(test_calcFullStartRaisingStochastic, "test_calcFullStopRaisingStochastic")
        call run_test_case(test_calcFullStartLoweringStochastic, "test_calcFullStartLoweringStochastic")
        call run_test_case(test_calcFullStopRaisingStochastic, "test_calcFullStopRaisingStochastic")
        call run_test_case(test_calcFullStopLoweringStochastic, "test_calcFullStopLoweringStochastic")
        call run_test_case(test_calcSingleOverlapMixedStochastic, "test_calcSingleOverlapMixedStochastic")

        ! now i need to test the specific semi-start/stop routines
        call run_test_case(test_calcLoweringSemiStartStochastic, "test_calcLoweringSemiStartStochastic")
        call run_test_case(test_calcRaisingSemiStartStochastic, "test_calcRaisingSemiStartStochastic")
        call run_test_case(test_calcFullStopR2L_stochastic, "test_calcFullStopR2L_stochastic")
        call run_test_case(test_calcFullStopL2R_stochastic, "test_calcFullStopL2R_stochastic")

        call run_test_case(test_calcRaisingSemiStopStochastic, "test_calcRaisingSemiStopStochastic")
        call run_test_case(test_calcLoweringSemiStopStochastic, "test_calcLoweringSemiStopStochastic")

        ! and then the rest of the test should just be calling them.
        call run_test_case(test_calcFullStartL2R_stochastic, "test_calcFullStartL2R_stochastic")
        call run_test_case(test_calcFullStartR2L_stochastic, "test_calcFullStartR2L_stochastic")

        call run_test_case(test_calcDoubleLoweringStochastic, "test_calcDoubleLoweringStochastic")
        call run_test_case(test_calcDoubleRaisingStochastic, "test_calcDoubleRaisingStochastic")
        call run_test_case(test_calcDoubleL2R2L_stochastic, "test_calcDoubleL2R2L_stochastic")
        call run_test_case(test_calcDoubleR2L2R_stochastic, "test_calcDoubleR2L2R_stochastic")
        call run_test_case(test_calcDoubleL2R_stochastic, "test_calcDoubleL2R_stochastic")
        call run_test_case(test_calcDoubleR2L_stochastic, "test_calcDoubleR2L_stochastic")

        call run_test_case(test_createStochasticExcitation_double, "test_createStochasticExcitation_double")

        print *, ""
        print *, "guga_excitations stochastic tests passed!"
        print *, ""

    end subroutine test_guga_excitations_stochastic

    subroutine test_guga_excitations_exact

        print *, ""
        print *, "testing module: guga_excitations:"
        print *, ""
        call run_test_case(test_calcRemainingSwitches, "test_calcRemainingSwitches")
        call run_test_case(test_actHamiltonian, "test_actHamiltonian")
        call run_test_case(test_excitationIdentifier_single, "test_excitationIdentifier_single")
        call run_test_case(test_createSingleStart, "test_createSingleStart")
        call run_test_case(test_singleUpdate, "test_singleUpdate")
        call run_test_case(test_singleEnd, "test_singleEnd")
        call run_test_case(test_calcAllExcitations_single, "test_calcAllExcitations_single")
        call run_test_case(test_excitationIdentifier_double, "test_excitationIdentifier_double")
        call run_test_case(test_checkCompatibility, "test_checkCompatibility")
        call run_test_case(test_calcSingleOverlapLowering, "test_calcSingleOverlapLowering")
        call run_test_case(test_calcSingleOverlapRaising, "test_calcSingleOverlapRaising")
        call run_test_case(test_calcSingleOverlapMixed, "test_calcSingleOverlapMixed")
        call run_test_case(test_calcDoubleLowering, "test_calcDoubleLowering")
        call run_test_case(test_calcDoubleRaising, "test_calcDoubleRaising")
        call run_test_case(test_calcDoubleL2R, "test_calcDoubleL2R")
        call run_test_case(test_calcDoubleR2L, "test_calcDoubleR2L")
        call run_test_case(test_calcFullStopLowering, "test_calcFullStopLowering")
        call run_test_case(test_calcFullStopRaising, "test_calcFullStopRaising")
        call run_test_case(test_calcFullStartLowering, "test_calcFullStartLowering")
        call run_test_case(test_calcFullStartRaising, "test_calcFullStartRaising")
        call run_test_case(test_calcFullStartFullStopAlike, "test_calcFullStartFullStopAlike")
        call run_test_case(test_calcDoubleExcitation_withWeight, "test_calcDoubleExcitation_withWeight")
        call run_test_case(test_calcNonOverlapDouble, "test_calcNonOverlapDouble")
        call run_test_case(test_calcFullStartR2L, "test_calcFullStartR2L")
        call run_test_case(test_calcFullStartL2R, "test_calcFullStartL2R")
        call run_test_case(test_calcFullStopR2L, "test_calcFullStopR2L")
        call run_test_case(test_calcFullStopL2R, "test_calcFullStopL2R")
        call run_test_case(test_calcFullStartFullStopMixed, "test_calcFullStartFullStopMixed")
        call run_test_case(test_calcAllExcitations_double, "test_calcAllExcitations_double")

        print *, ""
        print *, "guga_excitations tests passed!"
        print *, ""

    end subroutine test_guga_excitations_exact

    subroutine test_guga_matrixElements

        print *, ""
        print *, " =============================================================="
        print *, "  ===== testing routines of module: guga_matrixElements: ===== "
        print *, " =============================================================="
        print *, ""

        call run_test_case(check_calcDiagExchange_nI, "check_calcDiagExchange_nI")
        call run_test_case(check_calcDiagMatEleGUGA_nI, "check_calcDiagMatEleGUGA_nI")
        call run_test_case(test_coupling_coeffs, "test_coupling_coeffs")

        print *, ""
        print *, " guga_matrixElements tests passed!"
        print *, ""

    end subroutine test_guga_matrixElements

    subroutine test_guga_data

        print *, ""
        print *, "testing module: guga_data:"
        print *, ""
        call run_test_case(test_getMixedFullStop, "test_getMixedFullStop")
        call run_test_case(test_getSingleMatrixElement, "test_getSingleMatrixElement")
        call run_test_case(test_getDoubleMatrixElement, "test_getDoubleMatrixElement")
        print *, ""
        print *, "guga_data tests passed!"
        print *, ""

    end subroutine test_guga_data

    subroutine test_pure_rdm_ind
        integer(int_rdm) :: rdm_ind, rdm_ind_orig

        print *, ""
        print *, "testing: pure_rdm_ind"

        rdm_ind_orig = contract_1_rdm_ind(1, 2)
        call assert_equals(rdm_ind_orig, pure_rdm_ind(rdm_ind_orig))

        rdm_ind = contract_1_rdm_ind(1, 2, 1)
        call assert_equals(rdm_ind_orig, pure_rdm_ind(rdm_ind))

        rdm_ind = contract_1_rdm_ind(1, 2, 1, 2)
        call assert_equals(rdm_ind_orig, pure_rdm_ind(rdm_ind))

        rdm_ind_orig = contract_2_rdm_ind(1, 2, 3, 4)
        call assert_equals(rdm_ind_orig, pure_rdm_ind(rdm_ind_orig))

        rdm_ind = contract_2_rdm_ind(1, 2, 3, 4, 2)
        call assert_equals(rdm_ind_orig, pure_rdm_ind(rdm_ind))

        rdm_ind = contract_2_rdm_ind(1, 2, 3, 4, 2, 10)
        call assert_equals(rdm_ind_orig, pure_rdm_ind(rdm_ind))

        print *, ""
        print *, "testing: pure_rdm_ind. DONE."
    end subroutine test_pure_rdm_ind

    subroutine test_contract_extract_2_rdm
        integer(int_rdm) :: ijkl
        integer :: i, j, k, l
        integer(int_rdm) :: ij, kl

        print *, ""
        print *, "testing: contract and extract 2 rdm index: "
        print *, ""

        ijkl = contract_2_rdm_ind(1, 1, 1, 1)
        call assert_equals(1_int_rdm, ijkl)
        call extract_2_rdm_ind(ijkl, i=i, j=j, k=k, l=l, ij_out=ij, kl_out=kl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(1, l)
        call assert_equals(1_int_rdm, ij)
        call assert_equals(1_int_rdm, kl)

        ijkl = contract_2_rdm_ind(1, 2, 3, 4)
!         call assert_equals(1_int_rdm, ijkl)
        call extract_2_rdm_ind(ijkl, i=i, j=j, k=k, l=l, ij_out=ij, kl_out=kl)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2_int_rdm, ij)
        call assert_equals(12_int_rdm, kl)

        ijkl = contract_2_rdm_ind(i=1, j=1, k=2, l=2)
!         call assert_equals(1_int_rdm, ijkl)
        call extract_2_rdm_ind(ijkl, i=i, j=j, k=k, l=l, ij_out=ij, kl_out=kl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(2, l)
        call assert_equals(1_int_rdm, ij)
        call assert_equals(6_int_rdm, kl)

        print *, ""

    end subroutine test_contract_extract_2_rdm

    subroutine test_contract_extract_1_rdm_with_excitInfo
        integer(int_rdm) :: rdm_ind
        integer :: i, j, excit_lvl, excit_typ

        print *, ""
        print *, "testing: contract and extract 1 rdm index with the additional excitInfo"

        rdm_ind = contract_1_rdm_ind(1, 1, 0, excit_type%invalid)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_lvl)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_lvl)
        call assert_equals(excit_type%invalid, excit_typ)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(excit_type%invalid, excit_typ)

        rdm_ind = contract_1_rdm_ind(1, 1, 0)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_lvl)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_lvl)
        call assert_equals(0, excit_typ)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_typ)

        rdm_ind = contract_1_rdm_ind(1, 1, excit_typ=excit_type%invalid)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_lvl)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(0, excit_lvl)
        call assert_equals(excit_type%invalid, excit_typ)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(excit_type%invalid, excit_typ)

        rdm_ind = contract_1_rdm_ind(1, 2, 1, excit_type%single)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, excit_lvl)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, excit_lvl)
        call assert_equals(excit_type%single, excit_typ)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(excit_type%single, excit_typ)

        rdm_ind = contract_1_rdm_ind(1, 2, 1)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, excit_lvl)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, excit_lvl)
        call assert_equals(excit_type%invalid, excit_typ)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(excit_type%invalid, excit_typ)

        rdm_ind = contract_1_rdm_ind(1, 2, excit_typ=excit_type%single)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(0, excit_lvl)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(0, excit_lvl)
        call assert_equals(excit_type%single, excit_typ)
        call extract_1_rdm_ind(rdm_ind, i, j, excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(excit_type%single, excit_typ)

        print *, ""
        print *, "testing: contrat and extract 1 rdm index with the additional excitInfo. DONE"

    end subroutine test_contract_extract_1_rdm_with_excitInfo

    subroutine test_contract_extract_2_rdm_with_excitInfo

        integer(int_rdm) :: rdm_ind, ij, kl
        integer :: i, j, k, l, excit_lvl, excit_typ
        print *, ""
        print *, "testing: contract and exctract 2 rdm index with excit info"

        rdm_ind = contract_2_rdm_ind(1, 1, 1, 1, 0, excit_type%weight)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, ij, kl, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(1, l)
        call assert_equals(0, excit_lvl)
        call assert_equals(excit_typ, excit_type%weight)

        rdm_ind = contract_2_rdm_ind(1, 2, 3, 4, 2, excit_type%non_overlap)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, ij, kl, excit_lvl, excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2, excit_lvl)
        call assert_equals(excit_typ, excit_type%non_overlap)

        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)

        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=excit_lvl, &
                               excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2, excit_lvl)
        call assert_equals(excit_typ, excit_type%non_overlap)

        rdm_ind = contract_2_rdm_ind(1, 2, 3, 4, excit_lvl=2)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=excit_lvl, &
                               excit_typ=excit_typ)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2, excit_lvl)
        call assert_equals(excit_typ, excit_type%invalid)

        rdm_ind = contract_2_rdm_ind(1, 2, 3, 4, excit_typ=excit_type%double_lowering)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=excit_lvl, &
                               excit_typ=excit_typ)

        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(0, excit_lvl)
        call assert_equals(excit_typ, excit_type%double_lowering)

        print *, ""
        print *, "testing: contract and exctract 2 rdm index with excit info. DONE!"

    end subroutine test_contract_extract_2_rdm_with_excitInfo

    subroutine test_contract_extract_1_rdm
        integer(int_rdm) :: rdm_ind
        integer :: i, j

        print *, ""
        print *, " testing: contract and extract 1 rdm index"
        print *, ""

        rdm_ind = contract_1_rdm_ind(1, 1)
        call assert_equals(1_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(1, j)

        rdm_ind = contract_1_rdm_ind(1, 2)
        call assert_equals(2_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(2, j)

        rdm_ind = contract_1_rdm_ind(2, 1)
        call assert_equals(5_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)

        rdm_ind = contract_1_rdm_ind(2, 2)
        call assert_equals(6_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(2, i)
        call assert_equals(2, j)

        print *, ""

    end subroutine test_contract_extract_1_rdm

    subroutine test_excitationIdentifier
        integer :: i, j, k, l
        type(ExcitationInformation_t) :: ex1, ex2
        character(*), parameter :: testFun = " excitationIdentifier"

        print *, ""
        print *, " Testing: ", testFun
        print *, ""

        i = 1
        j = 2
        k = 3
        l = 4

        ex1 = excitationIdentifier(i, j)
        ex2 = excitationIdentifier(i, j, k, l)

        print *, ""
        print *, testFun, " tests passed!"
        print *, ""

    end subroutine test_excitationIdentifier

    subroutine test_bitChecks
        ! checks the function isZero(ilut,sOrb), isOne(ilut,sOrb) etc.
        integer(n_int) :: ilut(0:nifguga)
        integer :: det(4)
        integer :: i
        character(*), parameter :: testFun = "bitChecks"
        type(CSF_Info_t) :: csf_i

        nel = 4
        det = [1, 2, 3, 6]
        ! make a valid ilut:
        call EncodeBitDet_guga(det, ilut)
        csf_i = CSF_Info_t(ilut)

        print *, ""
        print *, " Testing ", testFun
        print *, ""
        ! use variable i to avoid compiler warning
        i = 1; call assert_true(.not. isZero(ilut, i))
        i = 2; call assert_true(.not. isZero(ilut, i))
        i = 3; call assert_true(.not. isZero(ilut, i))
        i = 4; call assert_true(isZero(ilut, i))
        i = 1; call assert_true(.not. isOne(ilut, i))
        i = 2; call assert_true(isOne(ilut, i))
        i = 3; call assert_true(.not. isOne(ilut, i))
        i = 4; call assert_true(.not. isOne(ilut, i))
        i = 1; call assert_true(.not. isTwo(ilut, i))
        i = 2; call assert_true(.not. isTwo(ilut, i))
        i = 3; call assert_true(isTwo(ilut, i))
        i = 4; call assert_true(.not. isTwo(ilut, i))
        i = 1; call assert_true(isThree(ilut, i))
        i = 2; call assert_true(.not. isThree(ilut, i))
        i = 3; call assert_true(.not. isThree(ilut, i))
        i = 4; call assert_true(.not. isThree(ilut, i))
        print *, ""
        print *, testFun, " tests passed!"
        print *, ""

    end subroutine test_bitChecks

    subroutine test_identify_excitation
        integer, allocatable :: nI(:), nJ(:)
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(ExcitationInformation_t) :: excitInfo

        ! make a more thorough test on the excitation identifier
        ! do it for now for the specific 14 electron system where the
        ! errors show up..
        ! nI: 3333333

        print *, ""
        print *, "testing: identify_excitation"
        print *, ""
        nel = 14
        allocate (nI(nel))
        nI = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
        ! nJ: 311333322
        allocate (nJ(nel))
        nJ = [1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18]

        print *, convert_guga_to_ni([3, 3, 3, 3, 3, 3, 3], 7)
        print *, convert_guga_to_ni([3, 1, 1, 3, 3, 3, 3, 2, 2], 9)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)

        excitInfo = identify_excitation(ilutI, ilutJ)

        call print_excitInfo(excitInfo)

        ! nI: 3331212123
        nI = [1, 2, 3, 4, 5, 6, 7, 10, 11, 14, 15, 18, 19, 20]
        ! nJ: 311212121322
        nJ = convert_guga_to_ni([3, 1, 1, 2, 1, 2, 1, 2, 1, 3, 2, 2], 12)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)

        excitInfo = identify_excitation(ilutI, ilutJ)

        call print_excitInfo(excitInfo)

        print *, ""
        print *, "identify_excitation tests passed!"
        print *, ""
        deallocate (nI)
        deallocate (nJ)
        nel = 4

    end subroutine test_identify_excitation

    subroutine test_identify_excitation_and_matrix_element
        character(*), parameter :: this_routine = "test_identify_excitation_and_matrix_element"
        integer(n_int) :: ilutI(0:niftot), ilutG(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :), two_ex(:, :)
        integer :: nEx, i, nex_2, test_det(4), j, ind
        logical :: valid
        real(dp) :: pos(nSpatOrbs), neg(nSpatOrbs), diff
        HElement_t(dp) :: mat_ele
        type(CSF_Info_t) :: csf_G

        test_det = [1, 2, 3, 4]

        call EncodeBitDet(test_det, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutG)
        csf_G = CSF_Info_t(ilutG)
        call actHamiltonian(ilutG, csf_G, ex, nEx)

        print *, ""
        print *, "Testing matrix elements for nEx excitations of: ", nEx
        print *, ""
        call write_det_guga(6, ilutG, .true.)

        print *, ""
        print *, "Do the tests on only connected determinants:"
        print *, ""
        do i = 1, nEx

            excitInfo = identify_excitation(ilutG, ex(:, i))

            call assert_true(excitInfo%valid)

            if (excitInfo%typ /= excit_type%single) then
                call checkCompatibility(csf_G, excitInfo, valid, pos, neg)

                call assert_true(valid)

            end if

            call calc_guga_matrix_element(ilutG, csf_G, ex(:, i), CSF_Info_t(ex(:, i)), excitInfo, mat_ele, &
                                          .true.)

            diff = abs(extract_matrix_element(ex(:, i), 1) - mat_ele)

            if (diff < 1.0e-10) diff = 0.0_dp

            if (diff > EPS) then
                call write_det_guga(6, ilutG, .true.)
                call write_det_guga(6, ex(:, i), .true.)
                call print_excitInfo(excitInfo)
                print *, "spin change: ", excitInfo%spin_change
                call stop_all(this_routine, "wrong matrix element!")
            end if

        end do

        print *, ""
        print *, "connected determinants correct!"
        print *, ""

        ! also do the tests on non-connected CSFs..
        ! maybe apply the hamiltonian a second time to the list of generated
        ! CSFs and check with the original one.. some of them should not
        ! be connected then and have a zero matrix element!
        print *, ""
        print *, "do the test on only non-connected determinants:"
        print *, ""
        ! this might take some time..
        do i = 1, nEx

            call write_det_guga(6, ex(:, i), .true.)

            call actHamiltonian(ex(:, i), CSF_Info_t(ex(:, i)), two_ex, nex_2)

            ! in acthamiltonian the current_stepvector quantity is set
            ! for the ilut input..

            do j = 1, nex_2

                excitInfo = identify_excitation(ilutG, two_ex(:, j))

                call calc_guga_matrix_element(ilutG, csf_G, two_ex(:, j), CSF_Info_t(two_ex(:, j)), &
                                              excitInfo, mat_ele, .true.)

                if (abs(mat_ele) > EPS) then

                    ! this should only happen if two_ex is in the original ex
                    ! or it is ilutI
                    ind = binary_search(ex(0:nifd, 1:nex), two_ex(0:nifd, j))

                    ! is the matrix element here correct if i find somethin?
                    if (ind < 0 .and. (.not. DetBitEQ(two_ex(0:nifd, j), ilutG(0:nifd)))) then

                        print *, "something wrong!"
                        call stop_all(this_routine, "matrix element should be 0!")

                    else if (ind > 0) then

                        ! is the sign correct now??
                        diff = abs(extract_matrix_element(ex(:, ind), 1) - mat_ele)

                        if (diff < 1.0e-10) diff = 0.0_dp

                        if (diff > EPS) then
                            call stop_all(this_routine, "incorrect sign!")
                        end if
                    end if
                end if
            end do

            deallocate (two_ex)

        end do

        print *, ""
        print *, "non-connected determinants correctly 0!"
        print *, ""

        ! i should also test more than double excitaitons to see if i correctly
        ! identify impossible excitations

        print *, ""
        print *, "test_identify_excitation_and_matrix_element passed!"
        print *, ""

    end subroutine test_identify_excitation_and_matrix_element

!
    subroutine test_coupling_coeffs
!
        integer(n_int) :: ilutI(0:nifguga), ilutJ(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        HElement_t(dp) :: mat_ele
        integer(int_rdm), allocatable :: rdm_ind(:)
        real(dp), allocatable :: rdm_mat(:)
        integer :: i, j, k, l
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, " =============================================================="
        print *, " ====== testing the coupling coefficient calculation =========="
        print *, " =============================================================="
        print *, ""
!
        nel = 1
        call EncodeBitDet_guga([1], ilutI)
        csf_i = CSF_Info_t(ilutI)
        call EncodeBitDet_guga([3], ilutJ)
!
        call calc_guga_matrix_element(ilutI, CSF_Info_t(ilutI), ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
!
        ! Single excitations:
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(5_int_rdm, rdm_ind(1))
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!

        call calc_guga_matrix_element(ilutJ, CSF_Info_t(ilutJ), ilutI, CSF_Info_t(ilutI), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
!
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(2_int_rdm, rdm_ind(1))
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        nel = 2
        call EncodeBitDet_guga([1, 2], ilutI)
        csf_i = CSF_Info_t(ilutI)
        call EncodeBitDet_guga([1, 4], ilutJ)
!
        call calc_guga_matrix_element(ilutI, CSF_Info_t(ilutI), ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call assert_equals(5_int_rdm, rdm_ind(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, CSF_Info_t(ilutJ), ilutI, CSF_Info_t(ilutI), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call assert_equals(2_int_rdm, rdm_ind(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call EncodeBitDet_guga([1, 6], ilutJ)
!
        call calc_guga_matrix_element(ilutI, CSF_Info_t(ilutI), ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, CSF_Info_t(ilutJ), ilutI, CSF_Info_t(ilutI), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        nel = 3
        call EncodeBitDet_guga([1, 3, 4], ilutI)
        csf_i = CSF_Info_t(ilutI)
        call EncodeBitDet_guga([3, 4, 5], ilutJ)
!
        call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(-1.0_dp), mat_ele)
        call assert_equals(-1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, CSF_Info_t(ilutJ), ilutI, CSF_Info_t(ilutI), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(-1.0_dp), mat_ele)
        call assert_equals(-1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call EncodeBitDet_guga([1, 3, 5], ilutI)
        csf_i = CSF_Info_t(ilutI)
        call EncodeBitDet_guga([3, 5, 7], ilutJ)
!
        call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutj, CSF_Info_t(ilutJ), iluti, CSF_Info_t(iluti), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(4, j)
        call assert_equals(1, i)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        nel = 2
        call EncodeBitDet_guga([1, 2], ilutI)
        csf_i = CSF_Info_t(ilutI)
        call EncodeBitDet_guga([3, 4], ilutJ)
!
        call calc_guga_matrix_element(iluti, CSF_Info_t(ilutI), ilutj, CSF_Info_t(ilutj), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(2.0_dp), mat_ele)
        call assert_equals(2.0_dp, rdm_mat(1))
        call extract_2_rdm_ind(rdm_ind(1), i=i, j=j, k=k, l=l)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, CSF_Info_t(ilutJ), ilutI, CSF_Info_t(ilutI), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(2.0_dp), mat_ele)
        call assert_equals(2.0_dp, rdm_mat(1))
        call extract_2_rdm_ind(rdm_ind(1), i=i, j=j, k=k, l=l)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, k)
        call assert_equals(2, l)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call EncodeBitDet_guga([3, 6], ilutJ)
!
        call calc_guga_matrix_element(ilutI, CSF_Info_t(ilutI), ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, &
                                      t_hamil=.false., rdm_ind=rdm_ind, rdm_mat=rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call extract_2_rdm_ind(rdm_ind(1), i=i, j=j, k=k, l=l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))

        nel = 4
    end subroutine test_coupling_coeffs

    subroutine test_findSwitches
        integer(n_int) :: ilutI(0:nifguga), ilutJ(0:nifguga)
        type(CSF_Info_t) :: csf_i

        print *, ""
        print *, "testing findSwitches routines:"
        print *, ""
        ! 3300
        call EncodeBitDet_guga([1, 2, 3, 4], ilutI)
        csf_i = CSF_Info_t(ilutI)
        ilutJ = ilutI

        call assert_true(findFirstSwitch(ilutI, ilutJ, 1, 4) == -1)
        call assert_equals(6, findLastSwitch(ilutI, ilutJ, 1, 4))

        ! 1122
        call EncodeBitDet_guga([1, 3, 6, 8], ilutJ)

        ! i changed the definition of the find-switches routines againe..
        ! so now it is the old way again, where every change in the
        ! stepvector is considered, not only spin-changes..
        ! this relies on the correct index input then! to only provide
        ! the correct range, where a spin-change is suspected!
        ! and dont forget the convention that for the findfirstswitch the
        ! end index (but end-1)is not considered and for the findlastswitch the
        ! first index is not considered (but start+1)
        ! 3300
        ! 1122
        call assert_true(findFirstSwitch(ilutI, ilutJ, 1, 4) == 1)
        call assert_true(findLastSwitch(ilutI, ilutJ, 1, 4) == 4)

        call assert_true(findFirstSwitch(ilutI, ilutJ, 2, 3) == 2)
        call assert_true(findLastSwitch(ilutI, ilutJ, 1, 2) == 2)

        ! 1230
        call EncodeBitDet_guga([1, 4, 5, 6], ilutI)
        csf_i = CSF_Info_t(ilutI)
        ! 1122
        call EncodeBitDet_guga([1, 3, 6, 8], ilutJ)
!
        ! 1230
        ! 1122
        call assert_true(findFirstSwitch(ilutI, ilutJ, 1, 3) == 2)
        call assert_true(findFirstSwitch(ilutI, ilutJ, 2, 3) == 2)
        call assert_true(findLastSwitch(ilutI, ilutJ, 1, 4) == 4)
        ! for the find last switch we exclude the inputted first orbital!
        call assert_true(findLastSwitch(ilutI, ilutJ, 2, 3) == 3)
        call assert_true(findFirstSwitch(ilutI, ilutJ, 3, 4) == 3)
        call assert_true(findLastSwitch(ilutI, ilutJ, 3, 4) == 4)

        ! 1122
        call EncodeBitDet_guga([1, 3, 6, 8], ilutI)
        csf_i = CSF_Info_t(ilutI)
        ! 1212
        call EncodeBitDet_guga([1, 4, 5, 8], ilutJ)

        ! 1122
        ! 1212
        call assert_true(findFirstSwitch(ilutI, ilutJ, 2, 3) == 2)
        call assert_true(findLastSwitch(ilutI, ilutJ, 1, 3) == 3)
        call assert_true(findFirstSwitch(ilutI, ilutJ, 3, 4) == 3)
        call assert_true(findLastSwitch(ilutI, ilutJ, 2, 3) == 3)

        call EncodeBitDet_guga([1, 2, 5, 6], ilutI)
        csf_i = CSF_Info_t(ilutI)
        call EncodeBitDet_guga([1, 2, 7, 8], ilutJ)

        ! 3030
        ! 3003
        call assert_true(findFirstSwitch(ilutI, ilutJ, 1, 4) == 3)
        call assert_true(findFirstSwitch(ilutI, ilutJ, 1, 3) == -1)
        call assert_true(findLastSwitch(ilutI, ilutJ, 1, 4) == 4)
        call assert_equals(6, (findLastSwitch(ilutI, IlutJ, 1, 2)))

        call EncodeBitDet_guga([3, 4, 7, 8], ilutI)
        csf_i = CSF_Info_t(ilutI)

        ! 0303
        ! 3003
        call assert_equals(6, findLastSwitch(ilutI, ilutJ, 2, 4))

        print *, ""
        print *, "findSwitches tests passed!"
        print *, ""

    end subroutine test_findSwitches

    subroutine test_count_beta_orbs_ij
        integer(n_int) :: ilut(0:nifguga)
        type(CSF_Info_t) :: csf_i

        ! these routine now need the current_stepvector quantitiy!

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)
        print *, ""
        print *, "testing count_beta_orbs_ij:"
        print *, ""

        call assert_true(count_beta_orbs_ij(csf_i, 1, 4) == 0)
        call assert_true(count_beta_orbs_ij(csf_i, 1, 3) == 0)
        call assert_true(count_beta_orbs_ij(csf_i, 2, 4) == 0)

        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        call assert_true(count_beta_orbs_ij(csf_i, 1, 4) == 2)
        call assert_true(count_beta_orbs_ij(csf_i, 2, 4) == 1)

        call EncodeBitDet_guga([3, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        call assert_true(count_beta_orbs_ij(csf_i, 1, 4) == 1)

        print *, ""
        print *, "count_beta_orbs_ij tests passed!"
        print *, ""

    end subroutine test_count_beta_orbs_ij

    subroutine test_count_alpha_orbs_ij
        integer(n_int) :: ilut(0:nifguga)
        type(CSF_Info_t) :: csf_i

        ! this routines now need the current_stepvector quantity!
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        print *, ""
        print *, "testing count_alpha_orbs_ij:"
        print *, ""
        ! 3300
        call assert_true(count_alpha_orbs_ij(csf_i, 1, 4) == 0)
        call assert_true(count_alpha_orbs_ij(csf_i, 2, 3) == 0)

        call EncodeBitDet_guga([1, 4, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        ! 1230
        call assert_true(count_alpha_orbs_ij(csf_i, 1, 4) == 1)

        call EncodeBitDet_guga([3, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call assert_true(count_alpha_orbs_ij(csf_i, 1, 4) == 1)

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call assert_true(count_alpha_orbs_ij(csf_i, 1, 4) == 2)

        print *, ""
        print *, "count_alpha_orbs_ij tests passed!"
        print *, ""

    end subroutine test_count_alpha_orbs_ij

    subroutine test_getSpatialOccupation
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1, 2, 3, 6], ilut)

        print *, ""
        print *, "testing getSpatialOccupation(ilut, sOrb):"
        print *, ""
        call assert_true(getSpatialOccupation(ilut, 1) .isclose.2.0_dp)
        call assert_true(getSpatialOccupation(ilut, 2) .isclose.1.0_dp)
        call assert_true(getSpatialOccupation(ilut, 3) .isclose.1.0_dp)
        call assert_true(getSpatialOccupation(ilut, 4) .isclose.0.0_dp)
        print *, ""
        print *, "getSpatialOccupation tests passed!"
        print *, ""

    end subroutine test_getSpatialOccupation

    subroutine test_calcFullStartFullStopMixed
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 1, 4)

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)

        print *, ""
        print *, "testing calcFullStartFullStopMixed(ilut, exInfo, ex, num, posSwitch, negSwitch):"
        print *, ""

        call calcFullStartFullStopMixed(ilut, CSF_Info_t(ilut), excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 1212
        ! 1122

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 1, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 1, 2, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) + 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3, 2, 2, 3)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 1, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 1, 2, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) + 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3, 1, 1, 3)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 1, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 1, 2, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) + 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) + sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(4, 2, 2, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 1, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 1, 2, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) + 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) + sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        ! 1122
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 1, 4)

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)

        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1122
        ! 1122
        ! 1212

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 1, 2, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 2, 1, 2]))
        ! -1/2 + 1 -> +1/2
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3, 2, 2, 3)

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)

        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 1, 2, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 2, 1, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3, 1, 1, 3)

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 1, 2, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 2, 1, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) + sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(4, 2, 2, 4)

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 1, 2, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 2, 1, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:, 2), 1) + sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartFullStopMixed tests passed!"
        print *, ""

    end subroutine test_calcFullStartFullStopMixed

    subroutine test_calcMixedContribution
        integer(n_int) :: ilut(0:nifguga), t(0:nifguga)
        HElement_t(dp) :: temp

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)

        call EncodeBitDet_guga([1, 3, 6, 8], t)

        print *, ""
        print *, "testing calcMixedContribution(ilut,t,start,ende):"
        print *, ""
        ! 1212
        ! 1122
        call calc_mixed_contr_integral(ilut, CSF_Info_t(ilut), t, 1, 4, temp)
        call assert_true(temp .isclose. h_cast(0.0_dp))


        call calc_mixed_contr_integral(t, CSF_Info_t(t), ilut, 1, 4, temp)
        call assert_true(temp .isclose. h_cast(0.0_dp))

        print *, ""
        print *, "calcMixedContribution tests passed!"
        print *, ""

    end subroutine test_calcMixedContribution

    subroutine test_generate_excitation_guga_double
        integer :: nI(4), nJ(4), IC, excitMat(2, maxExcit), exFlag, nEx, pos
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer(n_int) :: ilutGi(0:nifguga), ilutGj(0:nifguga)
        logical :: tParity
        real(dp) :: pgen
        HElement_t(dp) :: HElGen, mat_ele
        type(excit_gen_store_type), target :: store
        integer(n_int), allocatable :: ex(:, :)
        integer(int_rdm) :: rdm_ind, rdm_ind_, rdm_ind_1
        real(dp) :: x0, x1, rdm_mat_ex, rdm_comb
        integer(int_rdm), allocatable :: rdm_ind_v(:)
        real(dp), allocatable :: rdm_mat(:)
        type(ExcitationInformation_t) :: excitInfo
        integer :: i

        exFlag = 1
        ! make only double excitations:
        pSingles = 0.0_dp
        pDoubles = 1.0_dp - pSingles

        print *, ""
        print *, "testing generate_excitation_guga for doubles"
        print *, ""
        ! 3300:
        nI = [1, 2, 3, 4]; ilutI = 0_n_int
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)

        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)

            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_equals(helgen, extract_h_element(ex(:, pos)))

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do

            call assert_true(pos > 0)
            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 3030
        nI = [1, 2, 5, 6]
        call EncodeBitDet(nI, ilutI)

        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)

        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)


        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v) .or. &
                             any(rdm_ind_ == conjugate_rdm_ind(rdm_ind_v, ic)))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i) .or. &
                    rdm_ind_ == conjugate_rdm_ind(rdm_ind_v(i), ic)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do
        end if

        ! 3003
        nI = [1, 2, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)


        if (pgen > EPS) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 0330
        nI = [3, 4, 5, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)


        if (pgen > EPS) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 0303
        nI = [3, 4, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)


        if (pgen > EPS) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 0033
        nI = [5, 6, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > EPS) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 1023
        nI = [1, 6, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)


        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 3102
        nI = [1, 2, 3, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 3120
        nI = [1, 2, 3, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 3012
        nI = [1, 2, 5, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj),  excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v) .or. &
                             any(rdm_ind_ == conjugate_rdm_ind(rdm_ind_v, ic)))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i) .or. &
                    rdm_ind_ == conjugate_rdm_ind(rdm_ind_v(i), ic)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 0312
        nI = [3, 4, 5, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, &
                                          rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 1230
        nI = [1, 4, 5, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 1203
        nI = [1, 4, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 1320
        nI = [1, 3, 4, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 1302
        nI = [1, 3, 4, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 1032
        nI = [1, 5, 6, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 0132
        nI = [3, 5, 6, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        ! 0123
        nI = [3, 6, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                    call assert_equals(rdm_mat(i), rdm_mat_ex)
                end if
            end do

        end if

        ! 1122
        nI = [1, 3, 6, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, global_csf_i, nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, global_csf_i, ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                    call assert_equals(rdm_mat(i), rdm_mat_ex)
                end if
            end do

        end if

        ! 1212
        nI = [1, 4, 5, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, CSF_Info_t(ilutI), ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            call extract_stochastic_rdm_info(IlutBits, ilutJ, rdm_ind, x0, x1)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            call calc_explicit_2_rdm_guga(ilutGi, CSF_Info_t(ilutGi), nex, ex)

            do i = 1, nex
                if (DetBitEQ(ex(0:nifd, i), ilutJ(0:nifd))) then
                    rdm_ind_1 = extract_rdm_ind(ex(:, i))
                    if (rdm_ind_1 == rdm_ind_) then
                        pos = i
                    end if
                end if
            end do
            call assert_true(pos > 0)

            ilutGj = ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutGj)

            rdm_mat_ex = extract_matrix_element(ilutGj, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilutGi, CSF_Info_t(ilutGi), ilutGj, CSF_Info_t(ilutGj), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if
        print *, ""
        print *, "generate_excitation_guga tests passed!"

    end subroutine test_generate_excitation_guga_double

    subroutine test_generate_excitation_guga_single
        integer :: nI(4), nJ(4), IC, excitMat(2, maxExcit), exFlag, nEx, pos
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot), ilutGi(0:nifguga)
        integer(n_int) :: ilutGj(0:nifguga)
        logical :: tParity
        real(dp) :: pgen
        HElement_t(dp) :: HElGen, mat_ele
        type(excit_gen_store_type), target :: store
        integer(n_int), allocatable :: ex(:, :)
        integer(int_rdm) :: rdm_ind, rdm_ind_
        integer(int_rdm), allocatable :: rdm_ind_v(:)
        real(dp), allocatable :: rdm_mat(:)
        real(dp) :: x0
        type(ExcitationInformation_t) :: excitInfo

        exFlag = 1
        ! make this store element ...
        print *, ""
        print *, "testing generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,exMat,tPar,pgen,hEl,store)"
        print *, ""
        ! test singles only first
        pSingles = 1.0_dp
        pDoubles = 0.0_dp

        ! 3300:
        nI = [1, 2, 3, 4]; ilutI = 0_n_int
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)

            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), rdm_ind_)
            call assert_equals(extract_matrix_element(ex(:, pos), 1), x0)

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 3030
        nI = [1, 2, 5, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 3003
        nI = [1, 2, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 0330
        nI = [3, 4, 5, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 0303
        nI = [3, 4, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 0033
        nI = [5, 6, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1023
        nI = [1, 6, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 3102
        nI = [1, 2, 3, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 3120
        nI = [1, 2, 3, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 3012
        nI = [1, 2, 5, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 0312
        nI = [3, 4, 5, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1230
        nI = [1, 4, 5, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1203
        nI = [1, 4, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1320
        nI = [1, 3, 4, 6]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1302
        nI = [1, 3, 4, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1032
        nI = [1, 5, 6, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 0132
        nI = [3, 5, 6, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutI)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 0123
        nI = [3, 6, 7, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1122
        nI = [1, 3, 6, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutI)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        ! 1212
        nI = [1, 4, 5, 8]
        call EncodeBitDet(nI, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutGi)
        global_csf_i = CSF_Info_t(ilutGi)
        call generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, excitMat, &
                                      tParity, pgen, HElGen, store)
        call convert_ilut_toGUGA(ilutJ, ilutGj)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilutI, global_csf_i, ex, nEx)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)

            rdm_ind = extract_stochastic_rdm_ind(IlutBits, ilutJ)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilutGi, global_csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), ilutJ(0:nifd))
            call assert_true(pos > 0)

            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(IlutBits, ilutJ))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilutI, global_csf_i, ex(:, pos), CSF_Info_t(ex(:, pos)), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(IlutBits, ilutJ)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, helgen)
            call assert_equals(rdm_mat(1), x0)

        end if

        print *, ""
        print *, "generate_excitation_guga tests passed!"
        print *, ""

    end subroutine test_generate_excitation_guga_single

    subroutine test_calcDoubleR2L_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 3, 4, 2)

        print *, ""
        print *, "testing calcDoubleR2L_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call assert_true(excitInfo%typ == excit_type%double_R_to_L)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call calcDoubleR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 3003

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3, 0, 0, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_R_to_L, ex_typ)

        ! mixed: -1
        ! nonover: +2 -> +2

        ! 0132
        ! 1023
        call EncodeBitDet_guga([3, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 2, 1, 3)

        call assert_true(excitInfo%typ == excit_type%double_R_to_L)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 0, 2, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(2, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_R_to_L, ex_typ)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 2, 1, 3)

        call assert_true(excitInfo%typ == excit_type%double_R_to_L)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 0, 2, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(2, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_R_to_L, ex_typ)

        ! nonoverlap : -2
        ! mixed: +1 -> -1
        print *, ""
        print *, "calcDoubleR2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleR2L_stochastic

    subroutine test_calcDoubleL2R_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(3, 1, 2, 4)

        print *, ""
        print *, "testing calcDoubleL2R_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call assert_true(excitInfo%typ == excit_type%double_L_to_R)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 0330

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 3, 3, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_L_to_R, ex_typ)

        ! mixed matele: -1
        ! nonover: 2 -> +1

        ! 3102
        ! 1320
        call EncodeBitDet_guga([1, 2, 3, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(2, 4, 3, 1)

        call assert_true(excitInfo%typ == excit_type%double_L_to_R)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 2, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(2, i)
        call assert_equals(4, j)
        call assert_equals(3, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_L_to_R, ex_typ)

        ! mixed: -2
        ! nonover: +1 -> -1

        print *, ""
        print *, "calcDoubleL2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleL2R_stochastic

    subroutine test_calcDoubleR2L2R_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 3, 2)

        print *, ""
        print *, "testing calcDoubleR2L2R_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call assert_true(excitInfo%typ == excit_type%double_R_to_L_to_R)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleR2L2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 3030
        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3, 0, 3, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(3, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_R_to_L_to_R, ex_typ)

        ! 0123
        ! 1032
        call EncodeBitDet_guga([3, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(3, 2, 1, 4)

        call assert_true(excitInfo%typ == excit_type%double_R_to_L_to_R)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleR2L2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 0, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < EPS)
        ! mixed ele: -1/2 - 3/2 = -2
        ! nonover: 1 -> -1

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(2, j)
        call assert_equals(1, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_R_to_L_to_R, ex_typ)

        ! mixes matele: -1
        ! nonoverlp: 2 -> +1
        print *, ""
        print *, "calcDoubleR2L2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleR2L2R_stochastic

    subroutine test_calcDoubleL2R2L_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 1, 2, 3)

        print *, ""
        print *, "testing calcDoubleL2R2L_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call assert_true(excitInfo%typ == excit_type%double_L_to_R_to_L)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleL2R2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 0303

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 3, 0, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        ! argh i have to consider the non-overlap one..
        ! ok: the first is: -1
        ! the non-overlap: 2 -> so in total: +1
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_L_to_R_to_L, ex_typ)

        ! 1032
        ! 0123
        call EncodeBitDet_guga([1, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(2, 3, 4, 1)

        call assert_true(excitInfo%typ == excit_type%double_L_to_R_to_L)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleL2R2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 1, 2, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        ! the mixed is -2
        ! the non-overlap: + 1 -> so -1 in total!
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_L_to_R_to_L, ex_typ)

        print *, ""
        print *, "calcDoubleL2R2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleL2R2L_stochastic

    subroutine test_calcDoubleRaisingStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 2, 3)

        print *, ""
        print *, "testing calcDoubleRaisingStochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call assert_true(excitInfo%typ == excit_type%double_raising)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 3300
        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 2.0_dp) < EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_raising, ex_typ)

        excitInfo = excitationIdentifier(2, 4, 1, 3)

        call assert_true(excitInfo%typ == excit_type%double_raising)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call calcDoubleRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 3300
        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 2.0_dp) < EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(2, i)
        call assert_equals(4, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_raising, ex_typ)

        ! 0132
        ! 1320

        ! encode det
        call EncodeBitDet_guga([3, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 2, 3)

        call assert_true(excitInfo%typ == excit_type%double_raising)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 2, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_raising, ex_typ)

        excitInfo = excitationIdentifier(1, 3, 2, 4)

        call assert_true(excitInfo%typ == excit_type%double_raising)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 2, 0]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_raising, ex_typ)

        print *, ""
        print *, "calcDoubleRaisingStochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleRaisingStochastic

    subroutine test_calcDoubleLoweringStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        ! 1212
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 1, 3, 2)

        print *, ""
        print *, "testing calcDoubleLoweringStochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call assert_true(excitInfo%typ == excit_type%double_lowering)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 0033
        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 0, 3, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        ! have to think about the other index comb too!
        call assert_true(abs(extract_matrix_element(ex, 1) + 2.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(3, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_lowering, ex_typ)

        excitInfo = excitationIdentifier(3, 2, 4, 1)

        call assert_true(excitInfo%typ == excit_type%double_lowering)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call calcDoubleLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)
        ! 1212
        ! 0033
        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 0, 3, 3]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        ! have to think about the other index comb too!
        call assert_true(abs(extract_matrix_element(ex, 1) + 2.0_dp) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(2, j)
        call assert_equals(4, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_lowering, ex_typ)

        ! 3120
        ! 1032
        call EncodeBitDet_guga([1, 2, 3, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 1, 3, 2)

        call assert_true(excitInfo%typ == excit_type%double_lowering)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 0, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(3, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_lowering, ex_typ)

        excitInfo = excitationIdentifier(3, 1, 4, 2)

        call assert_true(excitInfo%typ == excit_type%double_lowering)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcDoubleLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 0, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp) < EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(4, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%double_lowering, ex_typ)

        print *, ""
        print *, "calcDoubleLoweringStochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleLoweringStochastic

    subroutine test_calcFullStopR2L_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 4, 2)

        print *, ""
        print *, "testing calcFullStopR2L_stochastic(ilut,csf_i, exinfo,ex,pgen):"
        print *, ""
        call assert_equals(excit_type%fullstop_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStopR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 3012

        ! since no switch in the overlap region happended this should be 0

        call assert_true(all(ex == 0_n_int))
        call assert_equals(0.0_dp, pgen)

        ! is there a valid possible?
        ! 1122
        ! 3012
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 4, 2)

        call assert_equals(excit_type%fullstop_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStopR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [3, 0, 1, 2], 4)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstop_R_to_L, ex_typ)
        call assert_equals(0.0_dp, x0)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 3, 3, 2)

        call assert_equals(excit_type%fullstop_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStopR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [3, 0, 1, 2], 4)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(2, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstop_R_to_L, ex_typ)
        call assert_equals(0.0_dp, x0)

        print *, ""
        print *, "calcFullStopR2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopR2L_stochastic

    subroutine test_calcFullStopL2R_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i
        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 1, 2, 4)

        print *, ""
        print *, "testing calcFullStopL2R_stochastic(ilut,exInfo,ex,pgen)"
        print *, ""
        call assert_equals(excit_type%fullstop_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStopL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 0312
        ! no possible excitation
        call assert_true(all(ex == 0_n_int))
        call assert_equals(0.0_dp, pgen)

        ! also do a possible one..

        ! encode det
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 1, 2, 4)

        call assert_equals(excit_type%fullstop_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStopL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1122
        ! 0312

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [0, 3, 1, 2], 4)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstop_L_to_R, ex_typ)
        call assert_equals(0.0_dp, x0)

        excitInfo = excitationIdentifier(2, 4, 4, 1)

        call assert_equals(excit_type%fullstop_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStopL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [0, 3, 1, 2], 4)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(2, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstop_L_to_R, ex_typ)
        call assert_equals(0.0_dp, x0)

        print *, ""
        print *, "calcFullStopL2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopL2R_stochastic

    subroutine test_calcFullStartR2L_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i
        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 2, 4, 1)

        print *, ""
        print *, "testing calcFullStartR2L_stochastic(ilut,exInfo,ex,pgen):"
        print *, ""
        call assert_equals(excit_type%fullstart_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! also should not yield a valid excitation
        ! 1212
        ! 1023
        call assert_equals(0.0_dp, pgen)
        call assert_true(all(ex == 0_n_int))

        ! 1122
        ! 1203
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 3, 4, 1)

        call assert_equals(excit_type%fullstart_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [1, 2, 0, 3], 4)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_R_to_L, ex_typ)
        call assert_equals(0.0_dp, x0)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 2, 2, 3)

        call assert_equals(excit_type%fullstart_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartR2L_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [1, 2, 0, 3], 4)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(2, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_R_to_L, ex_typ)
        call assert_equals(0.0_dp, x0)

        print *, ""
        print *, "calcFullStartR2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartR2L_stochastic

    subroutine test_calcFullStartL2R_stochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i
        ! encode det
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 2, 1)

        print *, ""
        print *, "testing calcFullStartL2R_stochastic(ilut, exInfo, ex, pgen)"
        print *, ""
        call assert_equals(excit_type%fullstart_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 1212
        ! 1320 -> normal single!
        call assert_equals(0.0_dp, pgen)
        call assert_true(all(ex == 0_n_int))

        ! 1122
        ! 1230 should work!
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 3, 1)

        call assert_equals(excit_type%fullstart_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_equals(calcStepVector(ex), [1, 2, 3, 0], 4)
        call assert_true(pgen > EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(3, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_L_to_R, ex_typ)
        call assert_equals(0.0_dp, x0)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(3, 2, 2, 4)

        call assert_equals(excit_type%fullstart_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartL2R_stochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_equals(calcStepVector(ex), [1, 2, 3, 0], 4)
        call assert_true(pgen > EPS)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(2, j)
        call assert_equals(2, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_L_to_R, ex_typ)
        call assert_equals(0.0_dp, x0)

        print *, ""
        print *, "calcFullStartL2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartL2R_stochastic

    subroutine test_calcRaisingSemiStopStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4), posSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 3, 4, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 2, 4, 1)

        call assert_true(excitInfo%typ == excit_type%fullstart_R_to_L)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        ! but switches are not yet set up... wtf
        weights = init_fullStartWeight(csf_i, 2, 4, negSwitch(2), posSwitch(2), &
                                       csf_i%B_real(2))

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "testing calcRaisingSemiStopStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitch, &
                                           posSwitch, ex, pgen)

        ! 1302: there should be 2 possible! why not?
        ! 1203.. ah yes.. since no overlap changes..
        call assert_true(pgen < EPS)
        call assert_true(all(ex == 0))

        ! 1122
        ! 1203
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(2, 3, 4, 2)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        ! but switches are not yet set up... wtf
        weights = init_fullStartWeight(csf_i, 3, 4, negSwitch(3), posSwitch(3), &
                                       csf_i%B_real(3))

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 2, 2, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        call calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitch, &
                                           posSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 2, 0, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcRaisingSemiStopStochastic tests passed!"
        print *, ""

    end subroutine test_calcRaisingSemiStopStochastic

    subroutine test_calcLoweringSemiStopStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4), posSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(2, 1, 1, 4)

        call assert_true(excitInfo%typ == excit_type%fullstart_L_to_R)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_fullStartWeight(csf_i, 2, 4, negSwitch(2), posSwitch(2), &
                                       csf_i%B_real(2))

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, pgen)

        print *, ""
        print *, "testing calcLoweringSemiStopStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitch, &
                                            posSwitch, ex, pgen)

        ! 1032
        ! 1230
        ! but this is again only a single..
        call assert_true(pgen < EPS)
        call assert_true(all(ex == 0_n_int))

        ! do:
        ! 1122
        ! 1230
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(2, 4, 3, 2)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        ! but switches are not yet set up... wtf
        weights = init_fullStartWeight(csf_i, 3, 4, negSwitch(3), posSwitch(3), &
                                       csf_i%B_real(3))

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 2, 2, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        call calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitch, &
                                            posSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 2, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 2) + sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "calcLoweringSemiStopStochastic tests passed!"
        print *, ""

    end subroutine test_calcLoweringSemiStopStochastic

    subroutine test_calcRaisingSemiStartStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4), posSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4, 1, 2, 4)
        ! 1032
        ! 0132
        ! is this even compatible??

        call assert_true(excitInfo%typ == excit_type%fullstop_L_to_R)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_semiStartWeight(csf_i, 2, 4, negSwitch(2), posSwitch(2), &
                                       csf_i%B_real(2))

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitch, &
                                          negSwitch, ex, pgen)

        ! 1032
        ! 0x
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 0, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1) - 1.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "testing calcRaisingSemiStartStochastic(ilut,exInfo,weigh,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcRaisingSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitch, &
                                            posSwitch, ex, pgen)

        ! 1032
        ! 01x
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0, 1, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "calcRaisingSemiStartStochastic tests passed!"
        print *, ""

    end subroutine test_calcRaisingSemiStartStochastic

    subroutine test_calcLoweringSemiStartStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4), posSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! encode det
        call EncodeBitDet_guga([1, 3, 4, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1, 4, 4, 2)

        call assert_true(excitInfo%typ == excit_type%fullstop_R_to_L)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_semiStartWeight(csf_i, 2, 4, negSwitch(2), &
                                       posSwitch(2), csf_i%B_real(2))

        ! modify the excitation so it fits test case:
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitch, &
                                          negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3, 3, 0, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1) - Root2) < 1.0e-10_dp)
        print *, ""
        print *, "testing calcLoweringSemiStartStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcLoweringSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitch, &
                                             posSwitch, ex, pgen)

        ! 1302
        ! 31x
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3, 1, 0, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex, 2) + sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)
        print *, ""
        print *, "calcLoweringSemiStartStochastic tests passed!"
        print *, ""

    end subroutine test_calcLoweringSemiStartStochastic

    subroutine test_calcSingleOverlapMixedStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! 0330
        call EncodeBitDet_guga([3, 4, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        excitInfo = excitationIdentifier(1, 3, 4, 3)

        print *, ""
        print *, "testing calcSingleOverlapMixedStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call assert_equals(excit_type%single_overlap_R_to_L, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcSingleOverlapMixedStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 0330
        ! 1302
        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(all(calcStepVector(ex) == [1, 3, 0, 2]))
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, abs(extract_matrix_element(ex, 2)))
        call assert_equals(-Root2, extract_matrix_element(ex, 1))

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%single_overlap_R_to_L, ex_typ)
        call assert_equals(0.0_dp, x1)
        call assert_equals(-Root2, x0)

        ! 3003
        call EncodeBitDet_guga([1, 2, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(3, 1, 3, 4)

        call assert_equals(excit_type%single_overlap_L_to_R, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcSingleOverlapMixedStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 3003
        ! 1032
        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(all(calcStepVector(ex) == [1, 0, 3, 2]))
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, abs(extract_matrix_element(ex, 2)))
        call assert_equals(Root2, extract_matrix_element(ex, 1), 1e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%single_overlap_L_to_R, ex_typ)
        call assert_equals(0.0_dp, x1)
        call assert_equals(Root2, x0, 1e-12_dp)

        print *, ""
        print *, "calcSingleOverlapMixedStochastic tests passed!"
        print *, ""

    end subroutine test_calcSingleOverlapMixedStochastic

    subroutine test_calcFullStopLoweringStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! 3030
        call EncodeBitDet_guga([1, 2, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 4, 3)

        print *, ""
        print *, "testing calcFullStopLoweringStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call assert_equals(excit_type%fullstop_lowering, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call calcFullStopLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)
        ! 3030
        ! 1023
        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_equals(calcStepVector(ex), [1, 0, 2, 3], 4)
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, abs(extract_matrix_element(ex, 2)))
        call assert_equals(-Root2, extract_matrix_element(ex, 1))

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstop_lowering, ex_typ)
        call assert_equals(0.0_dp, x1)
        call assert_equals(-Root2, x0)

        print *, ""
        print *, "calcFullStopLoweringStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopLoweringStochastic

    subroutine test_calcFullStopRaisingStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! 0303
        call EncodeBitDet_guga([3, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        excitInfo = excitationIdentifier(1, 4, 3, 4)

        print *, ""
        print *, "testing calcFullStopRaisingStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call assert_equals(excit_type%fullstop_raising, excitInfo%typ)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call calcFullStopRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! 0303
        ! 1320
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [1, 3, 2, 0]))
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, (extract_matrix_element(ex, 2)))
        call assert_equals(-Root2, extract_matrix_element(ex, 1))

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstop_raising, ex_typ)
        call assert_equals(0.0_dp, x1)
        call assert_equals(-Root2, x0)

        print *, ""
        print *, "calcFullStopRaisingStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopRaisingStochastic

    subroutine test_calcFullStartLoweringStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! 3030
        call EncodeBitDet_guga([1, 2, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        excitInfo = excitationIdentifier(3, 1, 4, 1)

        call assert_true(excitInfo%typ == excit_type%fullstart_lowering)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        call assert_true(.not. compFlag)
        print *, ""
        print *, "testing calcFullStartLoweringStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcFullStartLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(4, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_lowering, ex_typ)
        call assert_equals(0.0_dp, x1)

        excitInfo = excitationIdentifier(2, 1, 4, 1)

        call assert_true(excitInfo%typ == excit_type%fullstart_lowering)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call assert_true(compFlag)
        call calcFullStartLoweringStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)
        ! 3030
        ! 0132
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [0, 1, 3, 2]))
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + Root2) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(4, k)
        call assert_equals(1, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_lowering, ex_typ)
        call assert_equals(0.0_dp, x1)

        print *, ""
        print *, "calcFullStartLoweringStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartLoweringStochastic

    subroutine test_calcFullStartRaisingStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l, ex_lvl, ex_typ
        real(dp) :: x0, x1
        type(CSF_Info_t) :: csf_i

        ! 0033
        call EncodeBitDet_guga([5, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 3, 1, 4)

        call assert_true(excitInfo%typ == excit_type%fullstart_raising)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        print *, ""
        print *, "testing calcFullStartRaisingStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcFullStartRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! only result is: pgen should be 1..
        ! 0033
        ! 3012
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [3, 0, 1, 2]))
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + Root2) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_raising, ex_typ)
        call assert_equals(0.0_dp, x1)
        call assert_equals(-Root2, x0)

        excitInfo = excitationIdentifier(2, 3, 2, 4)

        call assert_true(excitInfo%typ == excit_type%fullstart_raising)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        call calcFullStartRaisingStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)
        ! 0033
        ! 0312
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [0, 3, 1, 2]))
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex, 1) + Root2) < 1.0e-10_dp)

        call extract_stochastic_rdm_info(GugaBits, ex, rdm_ind, x0, x1)
        call extract_2_rdm_ind(rdm_ind, i, j, k, l, excit_lvl=ex_lvl, excit_typ=ex_typ)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(4, l)
        call assert_equals(2, ex_lvl)
        call assert_equals(excit_type%fullstart_raising, ex_typ)
        call assert_equals(0.0_dp, x1)
        call assert_equals(-Root2, x0)

        print *, ""
        print *, "calcFullStartRaisingStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartRaisingStochastic

    subroutine test_mixedFullStopStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 1, 4)
        weights = init_doubleWeight(csf_i, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, pgen)

        call doubleUpdateStochastic(ilut, csf_i, 2, excitInfo, weights, negSwitch, posSwitch, ex, pgen)
        call doubleUpdateStochastic(ilut, csf_i, 3, excitInfo, weights, negSwitch, posSwitch, ex, pgen)

        ! i should never get the other matrix element.. due to the 0
        ! matrix element or?? hopefully!
        ! no! it is not 0!
        print *, ""
        print *, "testing mixedFullStopStochastic(ilut, excitInfo, ex)"
        print *, ""
        call mixedFullStopStochastic(ilut, csf_i, excitInfo, ex)

        if (isOne(ex, 3)) then
            call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
            call assert_true(abs(extract_matrix_element(ex, 1) + 1.0_dp / 2.0_dp) < 1.0e-10_dp)
        else if (isTwo(ex, 3)) then
            call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)
        end if

        print *, ""
        print *, "mixedFullStopStochastic tests passed!"
        print *, ""

    end subroutine test_mixedFullStopStochastic

    subroutine test_doubleUpdateStochastic
        character(*), parameter :: this_routine = "test_doubleUpdateStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! 1212
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 1, 4)
        weights = init_doubleWeight(csf_i, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, pgen)

        ! 1212
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 2, 1, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "testing doubleUpdateStochastic(ilut,orb,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call doubleUpdateStochastic(ilut, csf_i, 2, excitInfo, weights, negSwitch, posSwitch, ex, pgen)

        ! now there are 2 possibs.
        ! although.. do i exclude the "diagonal" excitation??
        ! because if yes, then there is only one possib here..
        call assert_true(pgen < 1.0_dp)
        if (isTwo(ex, 2)) then
            call assert_true(all(calcStepVector(ex) == [1, 2, 1, 2]))
            call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
            call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)
        else if (isOne(ex, 2)) then
            call assert_true(all(calcStepVector(ex) == [1, 1, 1, 2]))
            call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex, 2) + sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        else
            call stop_all(this_routine, "wrong stepvalue")
        end if

        call doubleUpdateStochastic(ilut, csf_i, 3, excitInfo, weights, negSwitch, posSwitch, ex, pgen)

        call assert_true(pgen < 1.0_dp)

        ! 121
        if (isOne(ex, 3)) then
            call assert_true(all(calcStepVector(ex) == [1, 2, 1, 2]))
            call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
            call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)

        else if (isTwo(ex, 3)) then
            call assert_true(all(calcStepVector(ex) == [1, 1, 2, 2]))
            call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex, 2) - OverR2) < 1.0e-10_dp)

        else
            call stop_all(this_routine, "wrong stepvalue!")

        end if

        print *, ""
        print *, "doubleUpdateStochastic tests passed!"
        print *, ""

    end subroutine test_doubleUpdateStochastic

    subroutine test_calcFullStartFullStopMixedStochastic
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(CSF_Info_t) :: csf_i

        ! set up determinant and excitaiton information
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 4, 1)

        call assert_true(excitInfo%typ == excit_type%fullstart_stop_mixed)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        print *, ""
        print *, "testing calcFullStartFullStopMixedStochastic(ilut,exInfo,ex,pgen)"
        print *, ""
        call calcFullStartFullStopMixedStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! in this constellation no excitaiton should be possible, due to 0
        ! matrix elements..
        call assert_true(all(ex == 0) .or. all(calcStepVector(ex) == [1, 1, 2, 2]))

        ! 1122
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 4, 1)

        call assert_true(excitInfo%typ == excit_type%fullstart_stop_mixed)

        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)
        print *, ""
        print *, "testing calcFullStartFullStopMixedStochastic(ilut,exInfo,ex,pgen)"
        print *, ""
        call calcFullStartFullStopMixedStochastic(ilut, csf_i, excitInfo, ex, pgen, posSwitches, negSwitches)

        ! in this constellation no excitaiton should be possible, due to 0
        ! matrix elements..
        call assert_true(all(ex == 0) .or. all(calcStepVector(ex) == [1, 2, 1, 2]))

        print *, ""
        print *, "calcFullStartFullStopMixedStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartFullStopMixedStochastic

    subroutine test_pickOrbitals_double
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int) :: ilut(0:nifguga)
        real(dp) :: pgen
        integer :: nI(4)
        type(CSF_Info_t) :: csf_i

        nI = [1, 2, 3, 4]
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        ! still have to think about, if i preemptively choose the different
        ! types of double excitation (iiij, ii,jk, etc.) or let i happen
        ! randomly and adjust the pgens accordingly...
        ! and what should i test here??

        print *, ""
        print *, "testing pickOrbitals_double(ilut, excitLvl):"
        print *, ""
        ! 3300
        call pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            ! what can i test here?
            ! only lowerings possible..
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(-1, excitInfo%gen1)
            call assert_equals(-1, excitInfo%gen2)
        end if

        call pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            ! what can i test here?
            ! only lowerings possible..
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(-1, excitInfo%gen1)
            call assert_equals(-1, excitInfo%gen2)
        end if

        nI = [5, 6, 7, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_i = CSF_Info_t(ilut)

        ! 0033
        call pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            ! what can i test here?
            ! only lowerings possible..
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(1, excitInfo%gen1)
            call assert_equals(1, excitInfo%gen2)
        end if

        call pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            ! what can i test here?
            ! only lowerings possible..
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(1, excitInfo%gen1)
            call assert_equals(1, excitInfo%gen2)
        end if

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_i = CSF_Info_t(ilut)

        call pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > EPS)
        end if
        print *, ""
        print *, "pickOrbitals_double tests passed!"
        print *, ""

    end subroutine test_pickOrbitals_double

    subroutine test_createStochasticExcitation_double
        integer(n_int) :: ilut(0:GugaBits%len_tot), ex(0:GugaBits%len_tot), &
                          ilutJ(0:GugaBits%len_tot)
        real(dp) :: pgen
        integer :: dummy(2), nI(4), pos, nex, i
        integer(n_int), allocatable :: all_ex(:, :)
        HElement_t(dp) :: helgen, mat_ele, mat_exact
        integer(int_rdm) :: rdm_ind, rdm_ind_, rdm_ind_1
        integer(int_rdm), allocatable :: rdm_ind_v(:)
        real(dp), allocatable :: rdm_mat(:)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: x0, x1, rdm_mat_ex, rdm_comb
        type(CSF_Info_t) :: csf_i

        nI = [1, 5, 6, 8]

        call EncodeBitDet_guga([1, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        print *, ""
        print *, "testing createStochasticExcitation_double(ilut, ex, pgen):"
        print *, ""
        call createStochasticExcitation_double(ilut, nI, csf_i, ex, pgen, dummy)

        ! what should i test here?
        if (pgen > EPS) then
            HElGen = extract_matrix_element(ex, 1)
            rdm_ind = extract_rdm_ind(ex)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            x0 = extract_stochastic_rdm_x0(GugaBits, ex)
            x1 = extract_stochastic_rdm_x1(GugaBits, ex)

            call actHamiltonian(ilut, csf_i, all_ex, nex)
            pos = binary_search(all_ex(0:nifd, 1:nex), ex(0:nifd))
            call assert_true(pos > 0)
            ilutJ = all_ex(:, pos)
            mat_exact = extract_matrix_element(ilutJ, 1)

            call assert_true(pos > 0)
            call assert_equals(HElGen, mat_exact)

            call calc_explicit_2_rdm_guga(ilut, csf_i, nex, all_ex)

            pos = binary_search(all_ex(0:nifd, 1:nex), ex(0:nifd))
            call assert_true(pos > 0)

            ilutJ = all_ex(:, pos)
            rdm_ind_1 = extract_rdm_ind(ilutJ)

            rdm_mat_ex = extract_matrix_element(ilutJ, 1)

            call assert_equals(rdm_ind_1, rdm_ind_)
            rdm_comb = combine_x0_x1(rdm_ind, x0, x1)
            call assert_equals(rdm_mat_ex, rdm_comb)

            call calc_guga_matrix_element(ilut, csf_i, ex, CSF_Info_t(ex), excitInfo, mat_ele, &
                                          t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)

            call assert_equals(mat_ele, HElGen)
            call assert_true(any(rdm_ind_ == rdm_ind_v))

            do i = 1, size(rdm_ind_v)
                if (rdm_ind_ == rdm_ind_v(i)) then
                    call assert_equals(rdm_mat(i), rdm_comb)
                end if
            end do

        end if

        print *, ""
        print *, "createStochasticExcitation_double tests passed!"
        print *, ""

    end subroutine test_createStochasticExcitation_double

    subroutine test_singleStochasticEnd
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! 3300
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)
        ! do not make it random here but choose specific excitation!
        excitInfo = excitationIdentifier(4, 1)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(csf_i, excitInfo%fullEnd)

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitch, &
                                          negSwitch, ex, pgen)

        call singleStochasticUpdate(ilut, csf_i, 2, excitInfo, weights, posSwitch, &
                                    negSwitch, ex, pgen)

        call singleStochasticUpdate(ilut, csf_i, 3, excitInfo, weights, posSwitch, &
                                    negSwitch, ex, pgen)

        print *, ""
        print *, "testing singleStochasticEnd(excitInfo, excitation):"
        print *, ""

        call singleStochasticEnd(csf_i, excitInfo, ex)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 0, 2]))
        call assert_true(abs(extract_matrix_element(ex, 1) + Root2) < 1.0e-10_dp)

        print *, ""
        print *, "singleStochasticEnd tests passed!"
        print *, ""

    end subroutine test_singleStochasticEnd

    subroutine test_singleStochasticUpdate
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen
        type(CSF_Info_t) :: csf_i

        ! 3300
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        ! do not make it random here but choose specific excitation!
!         excitInfo = pickOrbitals_single(ilut)
        excitInfo = excitationIdentifier(4, 1)

        call assert_true(excitInfo%typ == excit_type%single)
        call assert_true(excitInfo%fullStart == 1 .and. excitInfo%fullEnd == 4)
        call assert_true(excitInfo%gen1 == -1)

        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(csf_i, excitInfo%fullEnd)

        call assert_true(all(posSwitch < EPS))
        call assert_true(all(negSwitch < EPS))

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitch, &
                                          negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex, 1) - Root2) < 1.0e-10_dp)

        print *, ""
        print *, "testing singleStochasticUpdate(ilut, exInfo, weight, posSwitch, negSwitch, ex, pgen):"
        print *, ""
        call singleStochasticUpdate(ilut, csf_i, 2, excitInfo, weights, posSwitch, &
                                    negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex, 1) + Root2) < 1.0e-10_dp)

        call singleStochasticUpdate(ilut, csf_i, 3, excitInfo, weights, posSwitch, &
                                    negSwitch, ex, pgen)

        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex, 1) + Root2) < 1.0e-10_dp)

        print *, ""
        print *, "singleStochasticUpdate tests passed!"
        print *, ""

    end subroutine test_singleStochasticUpdate

    subroutine test_pickRandomOrb
        integer :: orb
        real(dp) :: pgen
        integer(n_int) :: ilut(0:nifguga)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 2, 3, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        pgen = 1.0_dp

        call pickRandomOrb_scalar(csf_i, 1, pgen, orb)
        call assert_true(pgen.isclose.1.0_dp / 3.0_dp)
        call assert_true(orb > 1 .and. orb <= 4)
        pgen = 1.0_dp
        call pickRandomOrb_forced(csf_i, 2, pgen, orb)
        call assert_true(orb == 1)
        call assert_true(pgen.isclose.1.0_dp)
        pgen = 1.0_dp
        call pickRandomOrb_vector(csf_i, [1, 2], pgen, orb)
        call assert_true(orb == 3 .or. orb == 4)
        call assert_true(pgen.isclose.1.0_dp / 2.0_dp)
        pgen = 1.0_dp
        call pickRandomOrb_vector(csf_i, [2, 3], pgen, orb, 2)
        call assert_true(orb == 4)
        call assert_true(pgen.isclose.1.0_dp)

        pgen = 1.0_dp
        call pickRandomOrb_scalar(csf_i, 0, pgen, orb, 0)
        call assert_true(pgen.isclose.1.0_dp / 3.0_dp)
        call assert_true(orb /= 3)

        pgen = 1.0_dp
        call pickRandomOrb_scalar(csf_i, 2, pgen, orb, 2)
        call assert_true(pgen.isclose.1.0_dp / 2.0_dp)
        call assert_true(orb == 3 .or. orb == 4)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(csf_i, 1, 2, pgen, orb)
        call assert_true(pgen.isclose.0.0_dp)
        call assert_true(orb == 0)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(csf_i, 1, 4, pgen, orb)
        call assert_true(pgen.isclose.1.0_dp / 2.0_dp)
        call assert_true(orb == 2 .or. orb == 3)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(csf_i, 1, 4, pgen, orb, 0)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(orb == 2)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(csf_i, 1, 4, pgen, orb, 1)
        call assert_true(pgen.isclose.1.0_dp)
        call assert_true(orb == 3)
        print *, ""
        print *, "pickRandomOrb tests passed!"
        print *, ""

    end subroutine test_pickRandomOrb

    subroutine test_mixedFullStartStochastic
        character(*), parameter :: this_routine = "test_mixedFullStartStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), prob
        type(CSF_Info_t) :: csf_i

        ! 1230
        call EncodeBitDet_guga([1, 4, 5, 6], ilut)

        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 1, 3)

        weights = init_doubleWeight(csf_i, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        print *, ""
        print *, "testing mixedFullStartStochastic:"
        print *, ""

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, prob)

        call assert_true(prob.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 2, 3, 0]))
        call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(4, 2, 2, 3)

        weights = init_doubleWeight(csf_i, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitch, &
                                      negSwitch, ex, prob)

        ! i think i have two possibs here..
        ! 1230
        ! 1212
        ! 1122
        call assert_true(prob < 1.0_dp)

        if (isTwo(ex, 2)) then
            call assert_true(all(calcStepVector(ex) == [1, 2, 3, 0]))
            call assert_true(abs(extract_matrix_element(ex, 1) + OverR2) < 1.0e-10_dp)
            call assert_true(abs(extract_matrix_element(ex, 2)) < EPS)

        else if (isOne(ex, 2)) then
            call assert_true(all(calcStepVector(ex) == [1, 1, 3, 0]))
            call assert_true(abs(extract_matrix_element(ex, 1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex, 2) - sqrt(3.0_dp) / 2.0_dp) < 1.0e-10_dp)

        else
            call stop_all(this_routine, "wrong stepvalue at fullstart!")
        endif

        print *, ""
        print *, "mixedFullStartStochastic tests passed!"
        print *, ""

    end subroutine test_mixedFullStartStochastic

    subroutine test_createStochasticStart_single
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), probWeight
        integer(n_int) :: ex(0:nifguga)
        type(CSF_Info_t) :: csf_i

        ! 3300
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1)

        call assert_true(excitInfo%typ == excit_type%single)
        call assert_true(excitInfo%fullStart == 1 .and. excitInfo%fullEnd == 4)
        call assert_true(excitInfo%gen1 == -1)

        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(csf_i, 4)

        call assert_true(all(posSwitch < EPS))
        call assert_true(all(negSwitch < EPS))

        print *, ""
        print *, "testing createStochasticStart_single(ilut,exInfo, weighs, posSwitch, negSwitch, ex, probWeight):"
        print *, ""
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitch, negSwitch, ex, probWeight)

        ! i should check the matrix element and the excitation to be sure
        ! about the effect!
        call assert_true(probWeight.isclose.1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex, 1) - Root2) < 1.0e-10_dp)

        print *, ""
        print *, "createStochasticStart_single tests passed!"
        print *, ""

    end subroutine test_createStochasticStart_single

    subroutine test_pickOrbitals_single
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        integer :: nI(4)
        type(CSF_Info_t) :: csf_i

        nI = [1, 2, 3, 4]

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        csf_i = CSF_Info_t(ilut)

        ! what should i test here?..
        print *, ""
        print *, "testing: pickOrbitals_single(ilut)"
        print *, ""

        ! 3300
        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == -1)
        end if

        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == -1)
        end if

        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == -1)
        end if

        ! 0033
        nI = [5, 6, 7, 8]
        call EncodeBitDet_guga([5, 6, 7, 8], ilut)

        csf_i = CSF_Info_t(ilut)

        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == 1)
        end if

        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == 1)
        end if

        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == 1)
        end if


        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_i = CSF_Info_t(ilut)

        ! 1212
        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
        end if

        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == excit_type%single)
        end if

        print *, ""
        print *, "pickOrbitals_single tests passed!"
        print *, ""

    end subroutine test_pickOrbitals_single

    subroutine test_createStochasticExcitation_single
        integer(n_int) :: ilut(0:nifguga), t(0:nifguga)
        real(dp) :: pgen
        integer :: nI(4), pos, nex
        integer(n_int), allocatable :: ex(:, :)
        integer(int_rdm) :: rdm_ind
        HElement_t(dp) :: mat_ele
        integer(int_rdm), allocatable :: rdm_ind_v(:)
        real(dp), allocatable :: rdm_mat(:)
        real(dp) :: x0
        integer(int_rdm) :: rdm_ind_
        type(ExcitationInformation_t) :: excitInfo
        type(CSF_Info_t) :: csf_i

        nI = [1, 2, 3, 4]
        call EncodeBitDet_guga(nI, ilut)
        csf_i = CSF_Info_t(ilut)

        print *, ""
        print *, "testing: createStochasticExcitation_single(ilut,t,weight):"
        print *, ""
        call createStochasticExcitation_single(ilut, nI, csf_i, t, pgen)

        if (pgen > 0.0_dp) then
            call actHamiltonian(ilut, csf_i, ex, nEx)

            rdm_ind = extract_rdm_ind(t)

            pos = binary_search(ex(0:nifd, 1:nex), t(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(extract_matrix_element(t, 1) - extract_matrix_element(ex(:, pos), 1)) < 1.0e-10_dp)
            call assert_equals(1, extract_excit_lvl_rdm(rdm_ind))
            call assert_equals(excit_type%single, extract_excit_type_rdm(rdm_ind))

            call calc_explicit_1_rdm_guga(ilut, csf_i, nEx, ex)
            pos = binary_search(ex(0:nifd, 1:nex), t(0:nifd), nifd)
            call assert_true(pos > 0)
            call assert_equals(extract_rdm_ind(ex(:, pos)), pure_rdm_ind(rdm_ind))
            call assert_equals(extract_matrix_element(ex(:, pos), 1), &
                               extract_stochastic_rdm_x0(GugaBits, t))

            ! also test with matrix element calculator!
            call calc_guga_matrix_element(ilut, csf_i, t, CSF_Info_t(t), excitInfo, &
                                          mat_ele, t_hamil=.true., rdm_ind=rdm_ind_v, rdm_mat=rdm_mat)
            x0 = extract_stochastic_rdm_x0(GugaBits, t)
            rdm_ind_ = pure_rdm_ind(rdm_ind)
            call assert_equals(1, size(rdm_ind_v))
            call assert_equals(rdm_ind_v(1), rdm_ind_)
            call assert_equals(mat_ele, extract_h_element(t))
            call assert_equals(rdm_mat(1), x0)

        end if

        print *, "createStochasticExcitation_single tests passed!"

    end subroutine test_createStochasticExcitation_single

    subroutine test_actHamiltonian
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), allocatable :: ex(:, :)
        integer :: nEx
        type(CSF_Info_t) :: csf_i

        nel = 4

        print *, ""
        print *, "testing actHamiltonian(ilut):"
        print *, ""
        ! 3300:
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(13, nEx)
        ! 0330
        call EncodeBitDet_guga([3, 4, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(14, nEx)
        ! 0303
        call EncodeBitDet_guga([3, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(14, nEx)
        ! 0033
        call EncodeBitDet_guga([5, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(13, nEx)
        ! 1023
        call EncodeBitDet_guga([1, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 3102
        call EncodeBitDet_guga([1, 2, 3, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 3120
        call EncodeBitDet_guga([1, 2, 3, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)

        ! 3030
        call EncodeBitDet_guga([1, 2, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(14, nEx)
        ! 3003:
        call EncodeBitDet_guga([1, 2, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(14, nEx)
        ! 3012
        call EncodeBitDet_guga([1, 2, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(16, nEx)
        ! 0312
        call EncodeBitDet_guga([3, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(16, nEx)
        ! 1230
        call EncodeBitDet_guga([1, 4, 5, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(16, nEx)
        ! 1203
        call EncodeBitDet_guga([1, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(16, nEx)
        ! 1320
        call EncodeBitDet_guga([1, 3, 4, 6], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 1302
        call EncodeBitDet_guga([1, 3, 4, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 1032
        call EncodeBitDet_guga([1, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 0132
        call EncodeBitDet_guga([3, 5, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 0123
        call EncodeBitDet_guga([3, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(17, nEx)
        ! 1122
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(12, nEx)
        ! 1212
        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        call actHamiltonian(ilut, csf_i, ex, nEx)
        print *, "number of excitations for: ", nEx
        call assert_equals(18, nEx)

        print *, ""
        print *, "actHamiltonian tests passed!"
        print *, ""

    end subroutine test_actHamiltonian

    subroutine test_calcAllExcitations_double
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), allocatable :: ex(:, :)
        integer :: nExcits
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)

        csf_i = CSF_Info_t(ilut)

        print *, ""
        print *, "testing calcAllExcitations_double(ilut,i,j,k,l,ex,nExits):"
        print *, ""
        call calcAllExcitations_double(ilut, csf_i, 1, 2, 3, 4, ex, nExcits)

        ! meh... was soll ich hier testen?
        ! 1212
        ! 3030
        call assert_true(nExcits == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [3, 0, 3, 0]))
        print *, ""
        print *, "calcAllExcitations_double tests passed!"
        print *, ""

    end subroutine test_calcAllExcitations_double

    subroutine test_calcFullStartFullStopAlike
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([3, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 1, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstart_stop_alike)

        print *, ""
        print *, "testing: calcFullStartFullStopAlike(ilut, exInfo, ex)"
        print *, ""
        call calcFullStartFullStopAlike(ilut, csf_i, excitInfo, ex)

        ! 0303
        ! 3300
        call assert_true(all(calcStepVector(ex(:, 1)) == [3, 3, 0, 0]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartFullStopAlike tests passed!"
        print *, ""

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 4, 1)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstart_stop_alike)

        print *, ""
        print *, "testing: calcFullStartFullStopAlike(ilut, exInfo, ex)"
        print *, ""
        call calcFullStartFullStopAlike(ilut, csf_i, excitInfo, ex)

        ! 3300
        ! 0303
        call assert_true(all(calcStepVector(ex(:, 1)) == [0, 3, 0, 3]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartFullStopAlike tests passed!"
        print *, ""

    end subroutine test_calcFullStartFullStopAlike

    subroutine test_calcFullStartL2R
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 3, 1)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstart_L_to_R)

        print *, ""
        print *, "testing: calcFullStartL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartL2R(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 1230

        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 3, 0]))

        excitInfo = excitationIdentifier(2, 4, 3, 2)
        call calcFullStartL2R(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 3, 0]))
        print *, ""
        print *, "calcFullStartL2R tests passed!"
        print *, ""

    end subroutine test_calcFullStartL2R

    subroutine test_calcFullStartR2L
        integer(n_int) :: ilut(0:2)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)

        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 3, 4, 1)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstart_R_to_L)

        print *, ""
        print *, "testing: calcFullStartR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartR2L(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 1203
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 0, 3]))

        excitInfo = excitationIdentifier(2, 3, 4, 2)
        call calcFullStartR2L(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 0, 3]))

        print *, ""
        print *, "calcFullStartR2L tests passed!"
        print *, ""

    end subroutine test_calcFullStartR2L

    subroutine test_calcFullStartRaising
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([3, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 1, 3)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstart_raising)

        print *, ""
        print *, "testing: calcFullStartRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartRaising(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0312
        ! 3300

        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [3, 3, 0, 0]))

        print *, ""
        print *, "calcFullStartRaising tests passed!"
        print *, ""

    end subroutine test_calcFullStartRaising

    subroutine test_calcFullStartLowering
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 2, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 3, 1)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstart_lowering)

        print *, ""
        print *, "testing: calcFullStartLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartLowering(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3012
        ! 0033
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [0, 0, 3, 3]))

        print *, ""
        print *, "calcFullStartLowering tests passed!"
        print *, ""

    end subroutine test_calcFullStartLowering

    subroutine test_calcFullStopR2L
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 4, 3)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstop_R_to_L)

        print *, ""
        print *, "testing: calcFullStopR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopR2L(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 3102
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [3, 1, 0, 2]))

        excitInfo = excitationIdentifier(1, 3, 3, 2)
        call calcFullStopR2L(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)
        ! 1212
        ! 3012
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [3, 0, 1, 2]))

        print *, ""
        print *, "calcFullStopR2L tests passed!"
        print *, ""

    end subroutine test_calcFullStopR2L

    subroutine test_calcFullStopL2R
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 3, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstop_L_to_R)

        print *, ""
        print *, "testing: calcFullStopL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopL2R(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 0132
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [0, 1, 3, 2]))

        excitInfo = excitationIdentifier(3, 1, 2, 3)
        call calcFullStopL2R(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 0312
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [0, 3, 1, 2]))

        print *, ""
        print *, "calcFullStopL2R tests passed!"
        print *, ""

    end subroutine test_calcFullStopL2R

    subroutine test_calcFullStopRaising
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([5, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4, 2, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstop_raising)

        print *, ""
        print *, "testing: calcFullStopRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopRaising(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0033
        ! 1230
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 3, 0]))

        print *, ""
        print *, "calcFullStopRaising tests passed!"
        print *, ""

    end subroutine test_calcFullStopRaising

    subroutine test_calcFullStopLowering
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1, 4, 2)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%fullstop_lowering)

        print *, ""
        print *, "testing: calcFullStopLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopLowering(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3300
        ! 1203
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 0, 3]))

        print *, ""
        print *, "calcFullStopLowering tests passed!"
        print *, ""

    end subroutine test_calcFullStopLowering

    subroutine test_calcDoubleR2L
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 3, 4, 2)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%double_R_to_L)

        print *, ""
        print *, "testing: calcDoubleR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleR2L(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 3003
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [3, 0, 0, 3]))

        print *, ""
        print *, " numExcits: ", num
        print *, ""
        print *, ex
        print *, ""
        print *, "calcDoubleR2L tests passed!"
        print *, ""

    end subroutine test_calcDoubleR2L

    subroutine test_calcDoubleL2R
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 4, 5, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(3, 1, 2, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%double_L_to_R)

        print *, ""
        print *, "testing: calcDoubleL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleL2R(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 0330
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [0, 3, 3, 0]))

        print *, ""
        print *, "calcDoubleL2R tests passed!"
        print *, ""

    end subroutine test_calcDoubleL2R

    subroutine test_calcDoubleRaising
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([5, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 3, 2, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%double_raising)

        print *, ""
        print *, "testing: calcDoubleRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleRaising(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0033
        ! 1122
        ! 1212

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 1, 2, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 2, 1, 2]))

        print *, "calcDoubleRaising tests passed!"

    end subroutine test_calcDoubleRaising

    subroutine test_calcDoubleLowering
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(3, 1, 4, 2)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)

        call assert_true(excitInfo%typ == excit_type%double_lowering)
        print *, ""
        print *, "testing: calcDoubleLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleLowering(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3300
        ! 1122
        ! 1212
        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 1, 2, 2]))
        call assert_true(all(calcStepVector(ex(:, 2)) == [1, 2, 1, 2]))

        print *, ""
        print *, "calcDoubleLowering tests passed!"
        print *, ""

    end subroutine test_calcDoubleLowering

    subroutine test_calcSingleOverlapRaising
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([5, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 3, 3, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapRaising(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0033
        ! 1032
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 0, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - Root2) < 1.0e-10_dp)

        print *, "calcSingleOverlapRaising tests passed!"

    end subroutine test_calcSingleOverlapRaising

    subroutine test_calcSingleOverlapMixed
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 2, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(3, 1, 3, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)

        print *, "testing: calcSingleOverlapMixed(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapMixed(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3003
        ! 1032
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 0, 3, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - Root2) < 1.0e-10_dp)

        print *, "calcSingleOverlapMixed tests passed!"

    end subroutine test_calcSingleOverlapMixed

    subroutine test_calcSingleOverlapLowering
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(2, 1, 4, 2)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapLowering(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3300
        ! 1302
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 3, 0, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - Root2) < 1.0e-10_dp)

        print *, "calcSingleOverlapLowering tests passed!"

    end subroutine test_calcSingleOverlapLowering

    subroutine test_calcNonOverlapDouble
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([3, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 2, 3, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%non_overlap)

        print *, "testing: calcNonOverlapDouble(ilut, exInfo, exs, num, posSwitch, negSwitch"
        call calcNonOverlapDouble(ilut, csf_i, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0303
        ! 1212

        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 2, 1, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) - 2.0_dp) < 1.0e-10_dp)

        ! todo: auch hier eine funktion noetig die nur die single excitations
        ! berechnt und nicht tmat einbezieht... und das auch effektiv macht

        print *, "calcNonOverlapDouble tests passed!"

    end subroutine test_calcNonOverlapDouble

    subroutine test_calcDoubleExcitation_withWeight
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: ex(:, :)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        logical :: compFlag
        type(CSF_Info_t) :: csf_i

        call EncodeBitDet_guga([3, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(2, 2, 1, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == excit_type%raising)

        print *, "testing: calcDoubleExcitation_withWeight(ilut, exInfo, exc, num)"
        call calcDoubleExcitation_withWeight(ilut, csf_i, excitInfo, ex, num, posSwitch, &
                                             negSwitch)

        ! 0303
        ! 1302
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:, 1)) == [1, 3, 0, 2]))
        call assert_true(abs(extract_matrix_element(ex(:, 1), 1) + 2.0_dp * Root2) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3, 3, 1, 4)
        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitch, negSwitch)
        call assert_true(.not. compFlag)

        ! todo! have to write a only excitation calculating funciton
        ! not including tmat for this case!

        print *, "calcDoubleExcitation_withWeight tests passed!"

    end subroutine test_calcDoubleExcitation_withWeight

    subroutine test_checkCompatibility
        real(dp) :: posSwitch(4), negSwitch(4)
        integer(n_int) :: ilut(0:nifguga)
        logical :: flag
        type(ExcitationInformation_t) :: excitInfo
        type(CSF_Info_t) :: csf_i

        print *, "testing: checkCompatibility:"
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier_double(1, 2, 3, 4)
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitch, negSwitch)
        call checkCompatibility(csf_i, excitInfo, flag, posSwitch, negSwitch)

        call assert_true(.not. flag)

        ! 3300

        excitInfo = excitationIdentifier_double(1, 2, 1, 4)
        call checkCompatibility(csf_i, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(.not. flag)

        excitInfo = excitationIdentifier_double(3, 2, 4, 1)
        call checkCompatibility(csf_i, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(flag)

        ! 0033
        call EncodeBitDet_guga([5, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier_double(1, 2, 1, 4)
        call checkCompatibility(csf_i, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(.not. flag)

        excitInfo = excitationIdentifier_double(1, 3, 1, 4)
        call checkCompatibility(csf_i, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(flag)

        print *, "checkCompatibility tests passed!"

    end subroutine test_checkCompatibility

    subroutine test_excitationIdentifier_double
        type(ExcitationInformation_t) :: excitInfo

        print *, "testing: excitationIdentifier_double"

        excitInfo = excitationIdentifier_double(1, 2, 3, 4)
        call assert_true(excitInfo%i == 1)
        call assert_true(excitInfo%j == 2)
        call assert_true(excitInfo%gen1 == 1)
        call assert_true(excitInfo%gen2 == 1)
        call assert_true(excitInfo%typ == excit_type%non_overlap)
        excitInfo = excitationIdentifier_double(1, 2, 2, 4)
        call assert_true(excitInfo%fullStart == 1)
        call assert_true(excitInfo%secondStart == 2)
        call assert_true(excitInfo%firstEnd == 2)
        call assert_true(excitInfo%fullEnd == 4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%typ == excit_type%single_overlap_raising)

        excitInfo = excitationIdentifier_double(3, 2, 3, 4)
        call assert_true(excitInfo%typ == excit_type%single_overlap_L_to_R)
        excitInfo = excitationIdentifier_double(4, 2, 3, 1)
        call assert_true(excitInfo%typ == excit_type%double_lowering)
        excitInfo = excitationIdentifier_double(1, 1, 3, 4)
        call assert_true(excitInfo%typ == excit_type%raising)
        excitInfo = excitationIdentifier_double(1, 1, 4, 4)
        call assert_true(excitInfo%typ == excit_type%weight)
        excitInfo = excitationIdentifier_double(1, 1, 1, 1)
        call assert_true(excitInfo%typ == excit_type%weight)
        excitInfo = excitationIdentifier_double(1, 3, 2, 4)
        call assert_true(excitInfo%typ == excit_type%double_raising)
        excitInfo = excitationIdentifier_double(1, 4, 3, 2)
        call assert_true(excitInfo%typ == excit_type%double_R_to_L_to_R)

        print *, "excitationIdentifier_double tests passed!"

    end subroutine test_excitationIdentifier_double

    subroutine test_calcAllExcitations_single
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), allocatable :: ex(:, :)
        integer :: nEx
        type(CSF_Info_t) :: csf_i

        ! 3300
        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        print *, "testing: calcAllExcitations_single"
        call calcAllExcitations_single(ilut, csf_i, 4, 1, ex, nEx)
        call assert_equals(-Root2, extract_matrix_element(ex(:, 1), 1))
        call assert_equals(1, nEx)

        deallocate (ex)

        ! 0123
        call EncodeBitDet_guga([3, 6, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        call calcAllExcitations_single(ilut, csf_i, 2, 4, ex, nEx)

        call assert_equals(1, nex)
        call assert_equals(-1.0_dp, extract_matrix_element(ex(:, 1), 1), 1e-10_dp)

        ! 0123
        ! 0312

        print *, "calcAllExcitations_single tests passed!"

    end subroutine test_calcAllExcitations_single

    subroutine test_singleEnd
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitch(nBasis / 2), negSwitch(nBasis / 2)
        integer(n_int), allocatable :: excits(:, :), tmpEx(:, :)
        integer :: num
        type(WeightObj_t) :: weights
        type(CSF_Info_t) :: csf_i

        print *, "testing: singleEnd(ilut, exInfo, tmpEx, nEx, excits)"
        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)
        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(csf_i, 4)

        call createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, weights, &
                               tmpEx, num)

        call assert_true(getDeltaB(tmpEx) == 1)
        call singleUpdate(ilut, csf_i, 2, excitInfo, posSwitch, negSwitch, weights, &
                          tmpEx, num)

        call assert_true(getDeltaB(tmpEx) == -1)
        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(tmpEx(:, 1), 1) + OverR2) < 1.0e-10_dp)

        call singleUpdate(ilut, csf_i, 3, excitInfo, posSwitch, negSwitch, weights, &
                          tmpEx, num)
        call assert_true(num == 1)
        call assert_true(getDeltaB(tmpEx) == -1)
        call assert_true(abs(extract_matrix_element(tmpEx(:, 1), 1) + OverR2) < 1.0e-10_dp)

        call singleEnd(ilut, csf_i, excitInfo, tmpEx, num, excits)

        call assert_true(.not. allocated(tmpEx))
        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(excits(:, 1), 1) + 1.0_dp) < 1.0e-10_dp)

        print *, "singleEnd tests passed!"

    end subroutine test_singleEnd

    subroutine test_singleUpdate
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitch(nBasis / 2), negSwitch(nBasis / 2)
        integer(n_int), allocatable :: excits(:, :)
        integer :: num
        type(WeightObj_t) :: weights
        type(CSF_Info_t) :: csf_i

        print *, "testing: singleUpdate(ilut,csf_i,orb,exInfo,posSwitch,negSwitch,excits,num)"
        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1, 4, 7, 8], ilut)

        csf_i = CSF_Info_t(ilut)
        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(csf_i, 4)

        call createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, weights, &
                               excits, num)

        call singleUpdate(ilut, csf_i, 2, excitInfo, posSwitch, negSwitch, weights, &
                          excits, num)

        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(excits(:, 1), 1) + OverR2) < 1.0e-10_dp)

        call singleUpdate(ilut, csf_i, 3, excitInfo, posSwitch, negSwitch, weights, &
                          excits, num)

        ! 1203
        ! 310x
        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(excits(:, 1), 1) + OverR2) < 1.0e-10_dp)

        print *, "singleUpdate tests passed!"

        deallocate (excits)

    end subroutine test_singleUpdate

    subroutine test_createSingleStart
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitch(nBasis / 2), negSwitch(nBasis / 2)
        integer(n_int), allocatable :: excits(:, :)
        integer :: num
        type(WeightObj_t) :: weights
        type(CSF_Info_t) :: csf_i

        print *, "testing: createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, excits, nExcits)"

        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(csf_i, 4)
        call createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, weights, &
                               excits, num)
        call assert_true(num == 1)
        call assert_true(extract_matrix_element(excits(:, 1), 1) .isclose.Root2)
        call assert_true(getDeltaB(excits(:, 1)) == 1)

        deallocate (excits)

        call EncodeBitDet_guga([1, 4, 7, 8], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(2, 4)

        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        weights = init_singleWeight(csf_i, 4)

        call createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, weights, &
                               excits, num)

        call assert_true(num == 1)
        call assert_true(getDeltaB(excits(:, 1)) == -1)

        deallocate (excits)

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 1)

        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        weights = init_singleWeight(csf_i, 4)

        call createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, weights, &
                               excits, num)
        call assert_true(num == 1)
        call assert_true(getDeltaB(excits(:, 1)) == -1)
        call assert_true(extract_matrix_element(excits(:, 1), 1) .isclose.Root2)

        ! and do one double start
        deallocate (excits)

        ! 1320
        ! 1122
        ! 1212
        call EncodeBitDet_guga([1, 3, 4, 6], ilut)
        csf_i = CSF_Info_t(ilut)

        excitInfo = excitationIdentifier(4, 2)

        weights = init_singleWeight(csf_i, 4)

        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitch, negSwitch)

        call createSingleStart(ilut, csf_i, excitInfo, posSwitch, negSwitch, weights, &
                               excits, num)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(excits(:, 1)) == [1, 1, 2, 0]))
        call assert_true(all(calcStepVector(excits(:, 2)) == [1, 2, 2, 0]))
        call assert_true(abs(extract_matrix_element(excits(:, 1), 1) - sqrt(3.0_dp / 2.0_dp)) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(excits(:, 2), 1) - OverR2) < 1.0e-10_dp)

        print *, "createSingleStart tests passed!"

        deallocate (excits)

    end subroutine test_createSingleStart

    subroutine test_excitationIdentifier_single
        type(ExcitationInformation_t) :: excitInfo

        print *, "testing: excitationIdentifier_single:"
        excitInfo = excitationIdentifier(1, 4)
        call assert_true(excitInfo%i == 1)
        call assert_true(excitInfo%j == 4)
        call assert_true(excitInfo%gen1 == 1)
        call assert_true(excitInfo%fullStart == 1)
        call assert_true(excitInfo%fullEnd == 4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == excit_type%single)

        excitInfo = excitationIdentifier(1, 4)
        call assert_true(excitInfo%i == 1)
        call assert_true(excitInfo%j == 4)
        call assert_true(excitInfo%gen1 == 1)
        call assert_true(excitInfo%fullStart == 1)
        call assert_true(excitInfo%fullEnd == 4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == excit_type%single)

        excitInfo = excitationIdentifier(2, 4)
        call assert_true(excitInfo%i == 2)
        call assert_true(excitInfo%j == 4)
        call assert_true(excitInfo%gen1 == 1)
        call assert_true(excitInfo%fullStart == 2)
        call assert_true(excitInfo%fullEnd == 4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == excit_type%single)

        excitInfo = excitationIdentifier(3, 2)
        call assert_true(excitInfo%i == 3)
        call assert_true(excitInfo%j == 2)
        call assert_true(excitInfo%gen1 == -1)
        call assert_true(excitInfo%fullStart == 2)
        call assert_true(excitInfo%fullEnd == 3)
        call assert_true(excitInfo%currentGen == -1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == excit_type%single)

        print *, "excitationIdentifier_single tests passed!"

    end subroutine test_excitationIdentifier_single

    subroutine test_getDoubleMatrixElement
        real(dp) :: x0, x1

        print *, "testing: getDoubleMatrixElement:"
        call getDoubleMatrixElement(1, 1, 0, 1, 1, 1.0_dp, 1.0_dp, x0, x1)
        ! todo more extensive tests here.. but not for now..
        print *, "x0 = ", x0, " x1 = ", x1
        print *, "getDoubleMatrixElement tests passed!"

    end subroutine test_getDoubleMatrixElement

    subroutine test_getMixedFullStop
        real(dp) :: x0, x1

        print *, "testing: getMixedFullStop:"
        call getMixedFullStop(1, 1, 0, 1.0_dp, x0, x1)
        call assert_true(abs(x0 - OverR2) < 1.0e-10_dp)
        call assert_true(abs(x1) < EPS)
        call getMixedFullStop(0, 0, 0, 1.0_dp, x0, x1)
        call assert_true(abs(x0) < EPS)
        call assert_true(abs(x1) < EPS)
        call getMixedFullStop(3, 3, 0, 2.0_dp, x0, x1)
        call assert_true(abs(x0 - Root2) < 1.0e-10_dp)
        call assert_true(abs(x1) < EPS)
        call getMixedFullStop(2, 2, 0, 3.0_dp, x0, x1)
        call assert_true(abs(x0 - OverR2) < 1.0e-10_dp)
        call assert_true(abs(x1 - sqrt(6.0_dp / 8.0_dp)) < 1.0e-10_dp)
        call getMixedFullStop(1, 2, 2, 2.0_dp, x0, x1)
        call assert_true(abs(x0) < EPS)
        call assert_true(abs(x1 - 1.0_dp) < 1.0e-10_dp)
        call getMixedFullStop(2, 1, -2, 1.0_dp, x0, x1)
        call assert_true(abs(x0) < EPS)
        call assert_true(abs(x1 - 1.0_dp) < 1.0e-10_dp)
        call getMixedFullStop(1, 1, 0, 2.0_dp, x0, x1)
        call assert_true(abs(x1 + 1.0_dp / sqrt(6.0_dp)) < 1.0e-10_dp)

        print *, "x0 = ", x0, " x1 = ", x1
        print *, "getMixedFullStop tests passed!"

    end subroutine test_getMixedFullStop

    subroutine test_getSingleMatrixElement

        print *, "testing: getSingleMatrixElement(d1,d2,dB,gen,b):"
        call assert_true(abs(getSingleMatrixElement(0, 0, -1, 1, 1.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0, 0, 1, 1, 1.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0, 0, 1, -1, 1.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0, 0, -1, -1, 1.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(3, 3, -1, 1, 1.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(3, 3, 1, 1, 1.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(3, 3, 1, -1, 1.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(3, 3, -1, -1, 1.0_dp) + 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(0, 1, -1, 1, 2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0, 1, 1, -1, 2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(0, 2, 1, 1, 2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0, 2, -1, -1, 2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1, 0, -1, 1, 2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1, 0, -1, -1, 2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2, 0, 1, 1, 2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2, 0, +1, -1, 2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1, 1, -1, 1, 3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1, 1, 1, 1, 3.0_dp) - sqrt(8.0_dp) / 3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1, 1, 1, -1, 3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1, 1, -1, -1, 2.0_dp) - sqrt(8.0_dp) / 3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2, 2, 1, 1, 3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2, 2, -1, 1, 1.0_dp) - sqrt(8.0_dp) / 3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2, 2, -1, -1, 3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2, 2, 1, -1, 2.0_dp) - sqrt(8.0_dp) / 3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1, 2, 1, 1, 2.0_dp) + 1.0_dp / 4.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1, 2, 1, -1, 2.0_dp) - 1.0_dp / 3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2, 1, -1, 1, 2.0_dp) - 1.0_dp / 2.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2, 1, -1, -1, 2.0_dp) + 1.0_dp / 3.0_dp) < tol)

        print *, "getSingleMatrixElement tests passed!"

    end subroutine test_getSingleMatrixElement

    subroutine test_matrix_element_ops
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        print *, "testing: encode_matrix_element(ilut, ele, type):"
        print *, "         extract_matrix_element(ilut, type):"
        print *, "         update_matrix_element(ilut, ele, type):"

        print *, "nifguga: ", nifguga
        print *, "niftot: ", niftot

        call encode_matrix_element(ilut, 1.0_dp, 1)
        call assert_true(extract_matrix_element(ilut, 1) .isclose.1.0_dp)
        call encode_matrix_element(ilut, -1.0_dp, 2)
        call assert_true(extract_matrix_element(ilut, 2) .isclose.-1.0_dp)
        call update_matrix_element(ilut, 2.0_dp, 1)
        call assert_true(extract_matrix_element(ilut, 1) .isclose.2.0_dp)
        call update_matrix_element(ilut, -2.0_dp, 2)
        call assert_true(extract_matrix_element(ilut, 2) .isclose.2.0_dp)

        print *, "encode_matrix_element tests passed!"
    end subroutine test_matrix_element_ops

    subroutine test_calcOcc_vector_ilut
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        print *, "testing: calcOcc_vector_ilut(ilut)"
        call assert_true(all([2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp] .isclose.calcOcc_vector_ilut(ilut)))

        call EncodeBitDet_guga([1, 4, 5, 6], ilut)
        call assert_true(all([1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp] .isclose.calcOcc_vector_ilut(ilut)))

        call EncodeBitDet_guga([1, 2, 3, 8], ilut)
        call assert_true(all([2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp] .isclose.calcOcc_vector_ilut(ilut)))
        print *, "calcOcc_vector_ilut tests passed!"

    end subroutine test_calcOcc_vector_ilut

    subroutine test_calcStepVector
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        print *, "testing: calcStepVector(ilut)"
        call assert_true(all([3, 3, 0, 0] == calcStepVector(ilut)))
        call EncodeBitDet_guga([1, 4, 5, 6], ilut)
        call assert_true(all([1, 2, 3, 0] == calcStepVector(ilut)))
        print *, "calcStepVector tests passed!"

    end subroutine test_calcStepVector

    subroutine test_isDouble

        print *, "testing: isDouble(nI, sOrb)"
        call assert_true(isDouble([1, 2, 3, 4], 1))
        call assert_true(isDouble([1, 2, 3, 6], 2))
        call assert_true(.not. isDouble([1, 2, 3, 6], 3))
        print *, "isDouble tests passed!"

    end subroutine test_isDouble

    subroutine test_isProperCSF_ilut
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        print *, "testing: isProperCSF_ilut(ilut)"
        call assert_true(isProperCSF_ilut(ilut))
        call EncodeBitDet_guga([2, 3, 4, 5], ilut)
        call assert_true(.not. isProperCSF_ilut(ilut))
        print *, "isProperCSF_ilut tests passed!"

    end subroutine test_isProperCSF_ilut

    subroutine test_set_get_DeltaB
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1, 2, 3, 4], ilut)

        print *, "testing: setDeltaB(deltaB,ilut)"
        call setDeltaB(-2, ilut)
        call assert_true(getDeltaB(ilut) == -2)
        call setDeltaB(-1, ilut)
        call assert_true(getDeltaB(ilut) == -1)
        call setDeltaB(0, ilut)
        call assert_true(getDeltaB(ilut) == 0)
        call setDeltaB(1, ilut)
        call assert_true(getDeltaB(ilut) == 1)
        call setDeltaB(2, ilut)
        call assert_true(getDeltaB(ilut) == 2)

        ! maybe problems with x1-matrix element
        call setDeltaB(-2, ilut)
        call assert_true(getDeltaB(ilut) == -2)
        call setDeltaB(-1, ilut)
        call assert_true(getDeltaB(ilut) == -1)
        call setDeltaB(0, ilut)
        call assert_true(getDeltaB(ilut) == 0)
        call setDeltaB(1, ilut)
        call assert_true(getDeltaB(ilut) == 1)
        call setDeltaB(2, ilut)
        call assert_true(getDeltaB(ilut) == 2)

        print *, "setDeltaB tests passed!"

    end subroutine test_set_get_DeltaB

    subroutine test_count_open_orbs_ij
        integer :: det(4)
        integer(n_int) :: ilut(0:nifguga)
        type(CSF_Info_t) :: csf_i

        det = [1, 2, 3, 4]
        call EncodeBitDet_guga(det, ilut)
        csf_i = CSF_Info_t(ilut)

        print *, "testing: count_open_orbs_ij(ilut, i, j)"
        print *, "open orbs in ([1,2,3,4],1,4): ", count_open_orbs_ij(csf_i, 1, 4)
        call assert_true(count_open_orbs_ij(csf_i, 1, 4) == 0)

        det = [1, 3, 5, 6]
        call EncodeBitDet_guga(det, ilut)
        csf_i = CSF_Info_t(ilut)

        print *, "open orbs in ([1, 3, 5, 6], 1, 4): ", count_open_orbs_ij(csf_i, 1, 4)
        call assert_true(count_open_orbs_ij(csf_i, 1, 4) == 2)

        print *, "count_open_orbs_ij tests passed!"

    end subroutine test_count_open_orbs_ij

    subroutine check_calcDiagMatEleGUGA_nI
        integer :: det(4)
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
        det = [1, 2, 3, 6]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3300
        det = [1, 2, 3, 4]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0330
        det = [3, 4, 5, 6]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0303
        det = [3, 4, 7, 8]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0033
        det = [5, 6, 7, 8]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1023
        det = [1, 6, 7, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3102
        det = [1, 2, 3, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3030
        det = [1, 2, 5, 6]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3003
        det = [1, 2, 7, 8]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3012
        det = [1, 2, 5, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0312
        det = [3, 4, 5, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1230
        det = [1, 4, 5, 6]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1203
        det = [1, 4, 7, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1320
        det = [1, 3, 4, 6]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1302
        det = [1, 3, 4, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1032
        det = [1, 5, 6, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0132
        det = [3, 5, 6, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0123
        det = [3, 6, 7, 8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1122
        det = [1, 3, 6, 8]
        call assert_equals(h_cast(10.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1212
        det = [1, 4, 5, 8]
        call assert_equals(h_cast(10.0_dp), calcDiagMatEleGUGA_nI(det))

        print *, this_routine, " tests passed!"

    end subroutine check_calcDiagMatEleGUGA_nI

    subroutine check_calcDiagExchange_nI
        integer :: det(4), iOrb, jOrb

        det = [1, 2, 3, 6]
        iOrb = 3
        jOrb = 4

        print *, "testing: calcDiagExchange_nI"
        print *, "H_ii = ", calcDiagMatEleGUGA_nI(det)
        print *, "X = ", calcDiagExchangeGUGA_nI(iOrb, jOrb, det)
        ! asserts would be nice..

        ! 3120

        print *, "calcDiagExchange tests passed!"

    end subroutine check_calcDiagExchange_nI

    subroutine test_calcbvector
        integer(n_int) :: ilut(0:nifguga)
        integer :: det(4)
        real(dp) :: checkB_nI(4), checkB_ilut(nBasis / 2)
        character(*), parameter :: testFun = "calcB_vector"

        print *, " Testing ", testFun
        det = [1, 2, 3, 4]
        checkB_nI = 0.0_dp
        call EncodeBitDet_guga(det, ilut)
        checkB_ilut = 0.0_dp

        call assert_true(all(calcB_vector_nI(det) .isclose.checkB_nI))
        call assert_true(all(calcB_vector_ilut(ilut) .isclose.0.0_dp))

        det = [1, 2, 3, 6]
        call EncodeBitDet_guga(det, ilut)

        checkB_nI(3) = 1.0_dp
        checkB_ilut(2) = 1.0_dp

        call assert_true(all(calcB_vector_nI(det) .isclose. [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]))
        call assert_true(all(calcB_vector_ilut(ilut) .isclose. [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]))
        print *, testFun, " tests passed!"

    end subroutine test_calcbvector

    subroutine test_calcRemainingSwitches()
        real(dp) :: neg(nBasis / 2), pos(nBasis / 2)
        integer :: det(4)
        integer(n_int) :: ilut(0:nifguga)
        type(CSF_Info_t) :: csf_i

        det = [1, 2, 3, 6] ! 3 1 2 0
        call EncodeBitDet_guga(det, ilut)
        csf_i = CSF_Info_t(ilut)

        print *, "***"
        print *, "testing calcRemainingSwitches:"
        call calcRemainingSwitches_single(csf_i, 1, 4, pos, neg)
        print *, "positive switches: ", pos
        call assert_true(all(pos.isclose. [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]))
        print *, "negative switches: ", neg
        call assert_true(all(neg.isclose. [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]))

        call EncodeBitDet_guga([1, 4, 5, 8], ilut) ! 1 2 1 2
        csf_i = CSF_Info_t(ilut)

        call calcRemainingSwitches_single(csf_i, 2, 4, pos, neg)
        print *, "positive switches: ", pos
        call assert_true(all(pos.isclose. [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]))
        print *, "negative switches: ", neg
        call assert_true(all(neg.isclose. [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]))
        print *, "calcRemainingSwitches tests successfull"

    end subroutine test_calcRemainingSwitches


end program test_guga
