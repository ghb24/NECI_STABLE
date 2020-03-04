#include "macros.h"
! GUGA testsuite.
! contains all GUGA related unit tests
! maybe for testing purposes hard-compile a test system with fixed parameters
! like electrons, orbitals etc...
! since otherwise all test cases are quite input dependent...
! discuss with simon how to implement that optimally
program test_guga

    ! check if i can use some sort of unit test framework, like fruit,
    ! and how to use it correctly and efficiently.

    use fruit

    use SystemData
    use guga_bitRepOps
    use guga_excitations
    use guga_matrixElements
    use guga_data
    use guga_types
    use guga_init
    use guga_procedure_pointers
    use guga_rdm, only: calc_all_excits_guga_rdm_singles, calc_explicit_1_rdm_guga, &
                        calc_all_excits_guga_rdm_doubles, t_mimic_stochastic, &
                        calc_explicit_diag_2_rdm_guga, calc_explicit_2_rdm_guga
    use constants
    use DetBitOps
    use Determinants
    use bit_reps
    use FciMCData
    use dsfmt_interface, only: dsfmt_init
    use genrandsymexcitnumod, only: testgenrandsymexcitnu
    use symrandexcit3, only: test_sym_excit3
    use util_mod, only: operator(.isclose.), near_zero, operator(.div.), &
                        binary_search
    use Integrals_neci, only: get_umat_el_normal
    use procedure_pointers, only: get_umat_el
    use read_fci, only: initfromfcid, fcidump_name, readfciint
    use UMatCache, only: tTransGTID, GetUMatSize, tumat2d, umat2d, tdeferred_umat2d
    use Parallel_neci, only: MPIInit
    use Calc, only: SetCalcDefaults, CalcInit
    use System, only: SetSysDefaults, SysInit
    use OneEInts, only: TMat2d
    use shared_memory_mpi, only: shared_allocate_mpi
    use IntegralsData, only: umat_win, umat
    use DetCalc, only: DetCalcInit
    use unit_test_helper_excitgen, only: generate_uniform_integrals
    use CalcData, only: t_guga_mat_eles

    use rdm_data_utils, only: calc_combined_rdm_label, calc_separate_rdm_labels

    implicit none

    real(dp), parameter :: tol = 1.0e-10_dp
    integer :: failed_count
    logical :: t_test_guga_excit_gen

    call init_fruit()
    call dsfmt_init(0)

    call guga_test_driver()

    ! change logical to test excit-gen
    t_test_guga_excit_gen = .false.
    if (t_test_guga_excit_gen) then
        nel = 0
        nbasis = 0
        call run_test_excit_gen_guga(nel, nbasis, stot)
    end if

    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

    subroutine guga_test_driver

        call init_guga_testsuite()

        call test_guga_bitRepOps
        call test_guga_excitations_stochastic
        call test_guga_excitations_exact
        call test_guga_matrixElements
        call test_guga_data
        call test_guga_explicit_rdms()

        call run_test_case(test_excitationIdentifier, "test_excitationIdentifier")
        call run_test_case(test_bitChecks, "test_bitChecks")
        call run_test_case(test_identify_excitation, "test_identify_excitation")
        call run_test_case(test_identify_excitation_and_matrix_element, "test_identify_excitation_and_matrix_element")

        !TODO maybe run the excit-gen test also!
        !call run_test_excit_gen_guga_S0

    end subroutine guga_test_driver

    subroutine init_guga_testsuite

        integer(int64) :: umatsize
        integer :: nBasisMax(5,3), lms, stot
        real(dp) :: ecore

        umatsize = 0
        nel = 4
        nbasis = 8
        nSpatOrbs = 4
        stot = 0
        lms = 0
        tGUGA = .true.

        call init_bit_rep()
        t_full_guga_tests = .true.

        tGen_sym_guga_mol = .true.
        tgen_guga_weighted = .true.
        tdeferred_umat2d = .true.
        tumat2d = .false.

        t_guga_mat_eles = .true.

        call init_guga()

        fcidump_name = "FCIDUMP"
        UMatEps = 1.0e-8
        tStoreSpinOrbs = .false.
        tTransGTID = .false.
        tReadFreeFormat = .true.

        call MPIInit(.false.)

        call dSFMT_init(8)

        call SetCalcDefaults()
        call SetSysDefaults()
        tReadInt = .true.

        call generate_uniform_integrals()

        get_umat_el => get_umat_el_normal

        call initfromfcid(nel,nbasismax,nBasis,lms,.false.)

        call GetUMatSize(nBasis, umatsize)

        allocate(TMat2d(nBasis,nBasis))

        call shared_allocate_mpi(umat_win, umat, (/umatsize/))

        call readfciint(UMat,umat_win,nBasis,ecore,.false.)
        call SysInit()
        ! required: set up the spin info

        call DetInit()
        call DetCalcInit()
        ! call SpinOrbSymSetup()

        call DetPreFreezeInit()

        call CalcInit()

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
        call run_test_case(test_calc_explicit_diag_2_rdm_guga, &
            "test_calc_explicit_diag_2_rdm_guga")
        call run_test_case(test_calc_explicit_2_rdm_guga, "test_calc_explicit_2_rdm_guga")

        call test_compare_RDM_indexing
        print *, ""
        print *, "explicit RDM routines passed!"
        print *, ""
    end subroutine test_guga_explicit_rdms

    subroutine test_compare_RDM_indexing

        integer(int_rdm) :: ijkl, abcd, ab, cd
        integer :: i, j, k, l, a, b, c, d, ij, kl

        print *, ""
        print *, " compare 'old' SD-based RDM indexing and GUGA convention"
        print *, ""

        call calc_combined_rdm_label(1,1,1,1, ijkl)
        abcd = contract_2_rdm_ind(1,1,1,1)

        call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
        call extract_2_rdm_ind(abcd, a, b, c, d, ab, cd)

        call assert_equals(i,a)
        call assert_equals(j,b)
        call assert_equals(k,c)
        call assert_equals(l,d)

        call calc_combined_rdm_label(1,2,3,4, ijkl)
        abcd = contract_2_rdm_ind(1,2,3,4)

        call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
        call extract_2_rdm_ind(abcd, a, b, c, d, ab, cd)

        call assert_equals(i,a)
        call assert_equals(j,b)
        call assert_equals(k,c)
        call assert_equals(l,d)


        call calc_combined_rdm_label(3,2,1,4, ijkl)
        abcd = contract_2_rdm_ind(3,2,1,4)

        call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
        call extract_2_rdm_ind(abcd, a, b, c, d, ab, cd)

        call assert_equals(i,a)
        call assert_equals(j,b)
        call assert_equals(k,c)
        call assert_equals(l,d)

        print *, ""
        print *, " compare 'old' SD-based RDM indexing and GUGA convention DONE"
        print *, ""

    end subroutine test_compare_RDM_indexing

    subroutine test_calc_explicit_diag_2_rdm_guga

        integer(n_int) :: ilut(0:nifguga)
        integer :: n_tot, i, j, k, l, iEx
        integer(n_int), pointer :: excits(:,:)
        real(dp) :: rdm_mat
        integer(int_rdm) :: rdm_ind
        integer, allocatable :: nJ(:)

        print *, ""
        print *, "testing: calc_explicit_diag_2_rdm_guga"
        print *, ""

        nel = 3
        call EncodeBitDet_guga([1,7,8], ilut)

        call calc_explicit_diag_2_rdm_guga(ilut, n_tot, excits)

        call assert_equals(0, n_tot)

        call EncodeBitDet_guga([1, 4, 5], ilut)
        call calc_explicit_diag_2_rdm_guga(ilut, n_tot, excits)

        call assert_equals(4, n_tot)

        ! 1 3 - 3 1
        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(1, l)

        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,1)), dp))


        ! 2 3 - 3 2
        rdm_ind = extract_rdm_ind(excits(:,2))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(2, l)

        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,2)), dp))

        ! 3 1 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,3))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,3)), dp))

        ! 3 2 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,4))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(2, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,4)), dp))



        call EncodeBitDet_guga([1, 3, 6], ilut)
        call calc_explicit_diag_2_rdm_guga(ilut, n_tot, excits)

        call assert_equals(4, n_tot)

        ! 1 3 - 3 1
        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(1, l)

        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,1)), dp),1e-12_dp)

        ! 2 3 - 3 2
        rdm_ind = extract_rdm_ind(excits(:,2))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(2, l)

        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,2)), dp),1e-12_dp)

        ! 3 1 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,3))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,3)), dp),1e-12_dp)

        ! 3 2 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,4))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(2, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,4)), dp),1e-12_dp)


        nel = 4
        call EncodeBitDet_guga([1, 3, 6, 8], ilut)

        call calc_explicit_diag_2_rdm_guga(ilut, n_tot, excits)

        ! 1 3 - 3 1
        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(1, l)

        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,1)), dp),1e-12_dp)

        ! 1 4 - 4 1
        rdm_ind = extract_rdm_ind(excits(:,2))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(1, l)

        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,2)), dp),1e-12_dp)

        ! 2 3 - 3 2
        rdm_ind = extract_rdm_ind(excits(:,3))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(2, l)
        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,3)), dp),1e-12_dp)

        ! 2 4 - 4 2
        rdm_ind = extract_rdm_ind(excits(:,4))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(2, i)
        call assert_equals(4, j)
        call assert_equals(4, k)
        call assert_equals(2, l)
        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,4)), dp),1e-12_dp)


        ! 3 1 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,5))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(3, l)

        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,5)), dp),1e-12_dp)

        ! 3 2 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,6))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(2, j)
        call assert_equals(2, k)
        call assert_equals(3, l)

        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,6)), dp),1e-12_dp)

        ! 4 1 - 1 4
        rdm_ind = extract_rdm_ind(excits(:,7))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(4, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(4, l)
        call assert_equals( sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,7)), dp),1e-12_dp)

        ! 4 2 - 2 4
        rdm_ind = extract_rdm_ind(excits(:,8))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(4, i)
        call assert_equals(2, j)
        call assert_equals(2, k)
        call assert_equals(4, l)
        call assert_equals( -sqrt(3.0_dp)/2.0_dp, real(extract_h_element(excits(:,8)), dp), 1e-12_dp)


        print *, ""
        print *, "testing: calc_explicit_diag_2_rdm_guga DONE"
        print *, ""

        nel = 4

    end subroutine test_calc_explicit_diag_2_rdm_guga

    subroutine test_calc_explicit_2_rdm_guga

        integer :: n_tot
        integer(n_int), pointer :: excits(:,:)
        integer, allocatable :: nJ(:)
        real(dp) :: rdm_mat
        integer(int_rdm) :: rdm_ind
        integer :: i, j, k, l
        integer(n_int) :: ilut(0:nifguga)


        nel = 2
        call EncodeBitDet_guga([5,6], ilut)

        print *, ""
        print *, "testing: calc_explicit_2_rdm_guga"
        print *, ""

        t_mimic_stochastic = .false.
        call calc_explicit_2_rdm_guga(ilut, n_tot, excits)

        call assert_equals(15, n_tot)
        !  1 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(2.0_dp, real(extract_h_element(excits(:,1)), dp), 1e-12_dp)

        ! 1 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,2))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,2)), dp), 1e-12_dp)

        ! 1 3 - 3 3
        rdm_ind = extract_rdm_ind(excits(:,3))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,3)), dp), 1e-12_dp)

        ! 1 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:,4))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,4)), dp), 1e-12_dp)

        ! 2 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,5))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,5)), dp), 1e-12_dp)

        ! 2 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,6))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(2.0_dp, real(extract_h_element(excits(:,6)), dp), 1e-12_dp)

        ! 2 3 - 3 3
        rdm_ind = extract_rdm_ind(excits(:,7))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,7)), dp), 1e-12_dp)

        ! 2 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:,8))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(2, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,8)), dp), 1e-12_dp)

        ! 3 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,9))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(3, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,9)), dp), 1e-12_dp)

        ! 3 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,10))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(3, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,10)), dp), 1e-12_dp)

        ! 3 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:,11))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(3, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,11)), dp), 1e-12_dp)

        ! 4 3 - 1 3
        rdm_ind = extract_rdm_ind(excits(:,12))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(1, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,12)), dp), 1e-12_dp)

        ! 4 3 - 2 3
        rdm_ind = extract_rdm_ind(excits(:,13))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,13)), dp), 1e-12_dp)

        ! 4 3 - 3 3
        rdm_ind = extract_rdm_ind(excits(:,14))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(3, k)
        call assert_equals(3, l)
        call assert_equals(sqrt(2.0_dp), real(extract_h_element(excits(:,14)), dp), 1e-12_dp)

        ! 4 3 - 4 3
        rdm_ind = extract_rdm_ind(excits(:,15))
        call extract_2_rdm_ind(rdm_ind, i, j, k, l)
        call assert_equals(4, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)
        call assert_equals(2.0_dp, real(extract_h_element(excits(:,15)), dp), 1e-12_dp)


        t_mimic_stochastic = .true.
        call calc_explicit_2_rdm_guga(ilut, n_tot, excits)
        call assert_equals(9, n_tot)

        print *, ""
        print *, "testing: calc_explicit_2_rdm_guga DONE"
        print *, ""

        t_mimic_stochastic = .false.
        nel = 4

    end subroutine test_calc_explicit_2_rdm_guga

    subroutine test_calc_all_excits_guga_rdm_doubles

        integer(n_int) :: ilut(0:nifguga)
        integer :: n_excits, i, j, k, l, iEx
        integer(n_int), pointer :: excits(:,:)
        integer, allocatable :: nJ(:)
        real(dp) :: rdm_mat
        integer(int_rdm) :: rdm_ind


        t_mimic_stochastic = .false.
        print *, ""
        print *, "testing: calc_all_excits_guga_rdm_doubles"
        print *, ""

        nel = 2
        call EncodeBitDet_guga([1,2], ilut)
        call init_csf_information(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, 2, 1, 0, 0, excits, n_excits)

        allocate(nJ(2))

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([1,4], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)

        call calc_all_excits_guga_rdm_doubles(ilut, 2, 1, 2, 1, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([3,4], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(2.0_dp, rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)

        call calc_all_excits_guga_rdm_doubles(ilut, 2, 1, 3, 1, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([3,6], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(2, i)
        call assert_equals(1, j)
        call assert_equals(3, k)
        call assert_equals(1, l)

        call calc_all_excits_guga_rdm_doubles(ilut, 3, 1, 2, 1, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([3,6], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(1, l)

        call EncodeBitDet_guga([3,4], ilut)
        call init_csf_information(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, 1, 2, 1, 2, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([1,2], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(2.0_dp, rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(1, k)
        call assert_equals(2, l)


        call EncodeBitDet_guga([5,6], ilut)
        call init_csf_information(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, 1, 3, 2, 3, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([1,4], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(2, k)
        call assert_equals(3, l)


        call calc_all_excits_guga_rdm_doubles(ilut, 1, 3, 4, 3, excits, n_excits)

        call assert_equals(1, n_excits)
        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([1,8], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(sqrt(2.0_dp), rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(1, i)
        call assert_equals(3, j)
        call assert_equals(4, k)
        call assert_equals(3, l)


        nel = 3

        call EncodeBitDet_guga([1,5,6], ilut)
        call init_csf_information(ilut)

        call calc_all_excits_guga_rdm_doubles(ilut, 3, 1, 2, 3, excits, n_excits)

        call assert_equals(1, n_excits)
        deallocate(nJ)
        allocate(nJ(nel))

        call decode_bit_det(nJ, excits(:,1))

        call assert_equals([3,5,6], nJ, 2)
        rdm_mat = real(extract_h_element(excits(:,1)), dp)
        call assert_equals(-1.0_dp, rdm_mat)

        rdm_ind = extract_rdm_ind(excits(:,1))
        call extract_2_rdm_ind(rdm_ind, i = i, j = j, k = k, l = l)
        call assert_equals(3, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(3, l)

        t_mimic_stochastic = .true.

        call calc_all_excits_guga_rdm_doubles(ilut, 3, 1, 2, 3, excits, n_excits)
        call assert_equals(0, n_excits)

        t_mimic_stochastic = .false.

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
        integer(n_int), pointer :: excits(:,:)
        integer, allocatable :: nJ(:)

        print *, ""
        print *, "testing: calc_explicit_1_rdm_guga"
        print *, ""

        nel = 2
        call EncodeBitDet_guga([1,2], ilut)

        call calc_explicit_1_rdm_guga(ilut, n_tot, excits)
        call assert_equals(3, n_tot)

        allocate(nJ(2))

        do iEx = 1, 3

            rdm_mat = real(extract_h_element(excits(:,iEx)), dp)
            call assert_equals(sqrt(2.0_dp), rdm_mat)

            rdm_ind = extract_rdm_ind(excits(:,iEx))
            call extract_1_rdm_ind(rdm_ind, i, j)

            call assert_equals(j, 1)
            call assert_equals(i, iex+1)

            call decode_bit_det(nJ, excits(:,iEx))

            call assert_equals([1,2*(iEx+1)], nJ, 2)
        end do

        call EncodeBitDet_guga([7,8], ilut)

        call calc_explicit_1_rdm_guga(ilut, n_tot, excits)
        call assert_equals(3, n_tot)

        do iEx = 1, 3

            rdm_mat = real(extract_h_element(excits(:,iEx)), dp)
            call assert_equals(sqrt(2.0_dp), rdm_mat)

            rdm_ind = extract_rdm_ind(excits(:,iEx))
            call extract_1_rdm_ind(rdm_ind, i, j)

            call assert_equals(j, 4)
            call assert_equals(i, iex)

            call decode_bit_det(nJ, excits(:,iEx))

            call assert_equals([2*iEx-1,8] , nJ, 2)
        end do


        print *, ""
        print *, "testing: calc_explicit_1_rdm_guga DONE"
        print *, ""

        nel = 4
    end subroutine test_calc_explicit_1_rdm_guga

    subroutine test_calc_all_excits_guga_rdm_singles

            integer(n_int) :: ilut(0:nifguga)
            integer(n_int), pointer :: excits(:,:)
            integer :: n_excits, i, j
            integer, allocatable :: nJ(:)
            real(dp) :: rdm_mat
            integer(int_rdm) :: rdm_ind


            print *, ""
            print *, "testing: calc_all_excits_guga_rdm_singles"
            print *, ""

            nel = 2
            call EncodeBitDet_guga([1,2], ilut)

            call init_csf_information(ilut)

            call calc_all_excits_guga_rdm_singles(ilut, 1, 2, excits, n_excits)

            call assert_equals(0, n_excits)

            call calc_all_excits_guga_rdm_singles(ilut, 2, 1, excits, n_excits)
            call assert_equals(1, n_excits)

            allocate(nJ(nel), source = 0)

            call decode_bit_det(nJ, excits(:,1))

            call assert_equals([1,4], nJ, 2)

            rdm_mat = real(extract_h_element(excits(:,1)), dp)

            call assert_equals(sqrt(2.0_dp), rdm_mat)

            rdm_ind = extract_rdm_ind(excits(:,1))
            call extract_1_rdm_ind(rdm_ind, i, j)
            call assert_equals(1, j)
            call assert_equals(2, i)

            call EncodeBitDet_guga([3,4],ilut)
            call init_csf_information(ilut)

            call calc_all_excits_guga_rdm_singles(ilut, 2, 1, excits, n_excits)
            call assert_equals(0, n_excits)

            call calc_all_excits_guga_rdm_singles(ilut, 1, 2, excits, n_excits)
            call assert_equals(1, n_excits)

            call decode_bit_det(nJ, excits(:,1))

            call assert_equals([1,4], nJ, 2)

            rdm_mat = real(extract_h_element(excits(:,1)), dp)

            rdm_ind = extract_rdm_ind(excits(:,1))
            call extract_1_rdm_ind(rdm_ind, i, j)
            call assert_equals(2, j)
            call assert_equals(1, i)

            call assert_equals(sqrt(2.0_dp), rdm_mat)

            nel = 3
            deallocate(nJ);

            call EncodeBitDet_guga([1,2,3], ilut)
            call init_csf_information(ilut)

            call calc_all_excits_guga_rdm_singles(ilut, 3, 1, excits, n_excits)

            call assert_equals(2, n_excits)
            allocate(nJ(nel), source = 0)
            call decode_bit_det(nJ, excits(:,2))
            call assert_equals([1,4,5], nJ, 3)

            rdm_ind = extract_rdm_ind(excits(:,2))
            call extract_1_rdm_ind(rdm_ind, i, j)
            call assert_equals(1, j)
            call assert_equals(3, i)

            call decode_bit_det(nJ, excits(:,1))
            call assert_equals([1,3,6], nJ, 3)

            rdm_ind = extract_rdm_ind(excits(:,1))
            call extract_1_rdm_ind(rdm_ind, i, j)
            call assert_equals(1, j)
            call assert_equals(3, i)

            print *, ""
            print *, " calc_all_excits_guga_rdm_singles passed!"
            print *, ""

            nel = 4

    end subroutine test_calc_all_excits_guga_rdm_singles

    subroutine test_guga_bitRepOps
        character(*), parameter :: this_routine = "test_guga_bitRepOps"

        print *, ""
        print *, "testing functions from module: guga_bitRepOps"
        print *, ""
        call run_test_case(test_findSwitches, "test_findSwitches")
        call run_test_case(test_count_alpha_orbs_ij, "test_count_alpha_orbs_ij")
        call run_test_case(test_count_beta_orbs_ij, "test_count_beta_orbs_ij")
        call run_test_case(test_matrix_element_ops, "test_matrix_element_ops")
        call run_test_case(test_count_open_orbs_ij, "test_count_open_orbs_ij")
        call run_test_case(test_getExcitationRangeMask, "test_getExcitationRangeMask")
        call run_test_case(test_set_get_DeltaB, "test_set_get_DeltaB")
        call run_test_case(test_isProperCSF_ilut, "test_isProperCSF_ilut")
        call run_test_case(test_calcbvector, "test_calcbvector")
        call run_test_case(test_isDouble, "test_isDouble")
        call run_test_case(test_calcStepVector, "test_calcStepVector")
        call run_test_case(test_getSpatialOccupation, "test_getSpatialOccupation")
        call run_test_case(test_calcOcc_vector_ilut, "test_calcOcc_vector_ilut")
        call run_test_case(test_contract_extract_1_rdm, "test_contract_extract_1_rdm")
        call run_test_case(test_contract_extract_2_rdm, "test_contract_extract_2_rdm")

        print *, ""
        print *, "guga_bitRepOps tests passed!"
        print *, ""

    end subroutine test_guga_bitRepOps

    subroutine test_guga_excitations_stochastic
        character(*), parameter :: this_routine = "test_guga_excitations_stochastic"

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
        character(*), parameter :: this_routine = "test_guga_excitations_exact"


        print *, ""
        print *, "testing module: guga_excitations:"
        print *, ""
        call  run_test_case(test_calcRemainingSwitches, "test_calcRemainingSwitches")
        call  run_test_case(test_actHamiltonian, "test_actHamiltonian")
        call  run_test_case(test_calcOverlapRange, "test_calcOverlapRange")
        call  run_test_case(test_excitationIdentifier_single, "test_excitationIdentifier_single")
        call  run_test_case(test_createSingleStart, "test_createSingleStart")
        call  run_test_case(test_singleUpdate, "test_singleUpdate")
        call  run_test_case(test_singleEnd, "test_singleEnd")
        call  run_test_case(test_calcAllExcitations_single, "test_calcAllExcitations_single")
        call  run_test_case(test_calcRemainingSwitches_double, "test_calcRemainingSwitches_double")
        call  run_test_case(test_excitationIdentifier_double, "test_excitationIdentifier_double")
        call  run_test_case(test_checkCompatibility, "test_checkCompatibility")
        call  run_test_case(test_calcSingleOverlapLowering, "test_calcSingleOverlapLowering")
        call  run_test_case(test_calcSingleOverlapRaising, "test_calcSingleOverlapRaising")
        call  run_test_case(test_calcSingleOverlapMixed, "test_calcSingleOverlapMixed")
        call  run_test_case(test_calcDoubleLowering, "test_calcDoubleLowering")
        call  run_test_case(test_calcDoubleRaising, "test_calcDoubleRaising")
        call  run_test_case(test_calcDoubleL2R, "test_calcDoubleL2R")
        call  run_test_case(test_calcDoubleR2L, "test_calcDoubleR2L")
        call  run_test_case(test_calcFullStopLowering, "test_calcFullStopLowering")
        call  run_test_case(test_calcFullStopRaising, "test_calcFullStopRaising")
        call  run_test_case(test_calcFullStartLowering, "test_calcFullStartLowering")
        call  run_test_case(test_calcFullStartRaising, "test_calcFullStartRaising")
        call  run_test_case(test_calcFullStartFullStopAlike, "test_calcFullStartFullStopAlike")
        call  run_test_case(test_calcDoubleExcitation_withWeight, "test_calcDoubleExcitation_withWeight")
        call  run_test_case(test_calcNonOverlapDouble, "test_calcNonOverlapDouble")
        call  run_test_case(test_calcFullStartR2L, "test_calcFullStartR2L")
        call  run_test_case(test_calcFullStartL2R, "test_calcFullStartL2R")
        call  run_test_case(test_calcFullStopR2L, "test_calcFullStopR2L")
        call  run_test_case(test_calcFullStopL2R, "test_calcFullStopL2R")
        call  run_test_case(test_calcFullStartFullStopMixed, "test_calcFullStartFullStopMixed")
        call  run_test_case(test_calcAllExcitations_double, "test_calcAllExcitations_double")

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

    subroutine test_contract_extract_2_rdm
        integer(int_rdm) :: ijkl
        integer :: i,j,k,l
        integer(int_rdm) :: ij, kl
        character(*), parameter :: this_routine = "test_contract_extract_2_rdm"

        print *, ""
        print *, "testing: contract and extract 2 rdm index: "
        print *, ""

        ijkl = contract_2_rdm_ind(1,1,1,1)
        call assert_equals(1_int_rdm, ijkl)
        call extract_2_rdm_ind(ijkl, i = i, j = j, k = k, l = l, ij_out = ij, kl_out = kl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(1, k)
        call assert_equals(1, l)
        call assert_equals(1_int_rdm, ij)
        call assert_equals(1_int_rdm, kl)

        ijkl = contract_2_rdm_ind(1,2,3,4)
!         call assert_equals(1_int_rdm, ijkl)
        call extract_2_rdm_ind(ijkl, i = i, j = j, k = k, l = l, ij_out = ij, kl_out = kl)
        call assert_equals(1, i)
        call assert_equals(2, j)
        call assert_equals(3, k)
        call assert_equals(4, l)
        call assert_equals(2_int_rdm, ij)
        call assert_equals(12_int_rdm, kl)

        ijkl = contract_2_rdm_ind(i = 1, j = 1, k = 2, l = 2)
!         call assert_equals(1_int_rdm, ijkl)
        call extract_2_rdm_ind(ijkl, i = i, j = j, k = k, l = l, ij_out = ij, kl_out = kl)
        call assert_equals(1, i)
        call assert_equals(1, j)
        call assert_equals(2, k)
        call assert_equals(2, l)
        call assert_equals(1_int_rdm, ij)
        call assert_equals(6_int_rdm, kl)

        print *, ""

    end subroutine test_contract_extract_2_rdm

    subroutine test_contract_extract_1_rdm
        integer(int_rdm) :: rdm_ind
        integer :: i, j
        character(*), parameter :: this_routine = "test_contract_extract_1_rdm"

        print *, ""
        print *, " testing: contract and extract 1 rdm index"
        print *, ""


        rdm_ind = contract_1_rdm_ind(1,1)
        call assert_equals(1_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(1, j)

        rdm_ind = contract_1_rdm_ind(1,2)
        call assert_equals(2_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(1, i)
        call assert_equals(2, j)

        rdm_ind = contract_1_rdm_ind(2,1)
        call assert_equals(5_int_rdm, rdm_ind)
        call extract_1_rdm_ind(rdm_ind, i, j)
        call assert_equals(2, i)
        call assert_equals(1, j)

        rdm_ind = contract_1_rdm_ind(2,2)
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
        ex2 = excitationIdentifier(i,j,k,l)

        print *, ""
        print *, testFun, " tests passed!"
        print *, ""

    end subroutine test_excitationIdentifier

    subroutine test_bitChecks
        ! checks the function isZero(ilut,sOrb), isOne(ilut,sOrb) etc.
        integer(n_int) :: ilut(0:nifguga)
        integer :: det(4)
        integer :: i
        character(*), parameter :: testFun = "bitChecks", &
            this_routine = "test_bitChecks"

        nel = 4
        det = [1,2,3,6]
        ! make a valid ilut:
        call EncodeBitDet_guga(det,ilut)

        print *, ""
        print *, " Testing ",testFun
        print *, ""
        ! use variable i to avoid compiler warning
        i = 1; call assert_true(.not.isZero(ilut,i))
        i = 2; call assert_true(.not.isZero(ilut,i))
        i = 3; call assert_true(.not.isZero(ilut,i))
        i = 4; call assert_true(isZero(ilut,i))
        i = 1; call assert_true(.not. isOne(ilut,i))
        i = 2; call assert_true(isOne(ilut,i))
        i = 3; call assert_true(.not.isOne(ilut,i))
        i = 4; call assert_true(.not.isOne(ilut,i))
        i = 1; call assert_true(.not.isTwo(ilut,i))
        i = 2; call assert_true(.not.isTwo(ilut,i))
        i = 3; call assert_true(isTwo(ilut,i))
        i = 4; call assert_true(.not.isTwo(ilut,i))
        i = 1; call assert_true(isThree(ilut,i))
        i = 2; call assert_true(.not.isThree(ilut,i))
        i = 3; call assert_true(.not.isThree(ilut,i))
        i = 4; call assert_true(.not.isThree(ilut,i))
        print *, ""
        print *, testFun, " tests passed!"
        print *, ""

    end subroutine test_bitChecks

    subroutine test_identify_excitation
        character(*), parameter :: this_routine = "test_identify_excitation"
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
        allocate(nI(nel))
        nI = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
        ! nJ: 311333322
        allocate(nJ(nel))
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

        print *, ""
        print *, "identify_excitation tests passed!"
        print *, ""
        deallocate(nI)
        deallocate(nJ)
        nel = 4

    end subroutine test_identify_excitation

    subroutine test_identify_excitation_and_matrix_element
        character(*), parameter :: this_routine = "test_identify_excitation_and_matrix_element"
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot), ilutG(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:), two_ex(:,:)
        integer :: nEx, i, nex_2, test_det(4), j, ind
        logical :: valid
        real(dp) :: pos(nSpatOrbs), neg(nSpatOrbs), diff
        HElement_t(dp) :: mat_ele

         test_det = [1,2,3,4]

        call EncodeBitDet(test_det, ilutI)
        call convert_ilut_toGUGA(ilutI, ilutG)
        call init_csf_information(ilutG)
        call actHamiltonian(ilutG, ex, nEx)


        print *, ""
        print *, "Testing matrix elements for nEx excitations of: ", nEx
        print *, ""
        call write_det_guga(6,ilutG,.true.)
        call write_guga_list(6,ex(:,1:nex))

        print *, ""
        print *, "Do the tests on only connected determinants:"
        print *, ""
        do i = 1, nEx

            excitInfo = identify_excitation(ilutG, ex(:,i))

            call assert_true(excitInfo%valid)

            if (excitInfo%typ /= 0) then
                call checkCompatibility(ilutG, excitInfo, valid, pos, neg)

                call assert_true(valid)

            end if

            call calc_guga_matrix_element(ilutG, ex(:,i), excitInfo, mat_ele, &
                .true., 1)

            diff = abs(extract_matrix_element(ex(:,i), 1) - mat_ele)

            if (diff < 1.0e-10) diff = 0.0_dp

            if (diff > EPS) then
                call write_det_guga(6,ilutG,.true.)
                call write_det_guga(6,ex(:,i),.true.)
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

            call write_det_guga(6, ex(:,i),.true.)

            call actHamiltonian(ex(:,i), two_ex, nex_2)

            ! in acthamiltonian the current_stepvector quantity is set
            ! for the ilut input..

            do j = 1, nex_2

                excitInfo = identify_excitation(ilutG, two_ex(:,j))

                call calc_guga_matrix_element(ilutG, two_ex(:,j), excitInfo, &
                    mat_ele, .true., 2)

                if (abs(mat_ele) > EPS) then

                    ! this should only happen if two_ex is in the original ex
                    ! or it is ilutI
                    ind = binary_search(ex(0:nifd,1:nex), two_ex(0:nifd,j))

                    ! is the matrix element here correct if i find somethin?
                    if (ind < 0 .and. (.not. DetBitEQ(two_ex(0:nifd,j),ilutG(0:nifd)))) then

                        print *, "something wrong!"
                        call stop_all(this_routine, "matrix element should be 0!")

                    else if (ind > 0) then

                        ! is the sign correct now??
                        diff = abs(extract_matrix_element(ex(:,ind),1) - mat_ele)

                        if (diff < 1.0e-10) diff = 0.0_dp

                        if (diff > EPS) then
                            print *, "sign check:"
                            print *, "I:"
                            call write_det_guga(6,ilutG,.true.)
                            print *, "step: ", temp_step_i
                            print *, "b:", temp_b_real_i
                            print *, "n:", temp_occ_i
                            print *, "J:"
                            call write_det_guga(6,two_ex(:,j),.true.)
                            print *, "step: ", temp_step_j
                            print *, "db: ", temp_delta_b

                            call stop_all(this_routine, "incorrect sign!")

                        end if
                    end if
                end if
            end do

            deallocate(two_ex)

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
        character(*), parameter :: this_routine = "test_coupling_coeffs"
!
        print *, ""
        print *, " =============================================================="
        print *, " ====== testing the coupling coefficient calculation =========="
        print *, " =============================================================="
        print *, ""
!
        nel = 1
        call EncodeBitDet_guga([1], ilutI)
        call EncodeBitDet_guga([3], ilutJ)
!
        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
!
        ! Single excitations:
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(5_int_rdm, rdm_ind(1))
        call assert_equals(2,i)
        call assert_equals(1,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!

        call calc_guga_matrix_element(ilutJ, ilutI, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
!
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(2_int_rdm, rdm_ind(1))
        call assert_equals(1,i)
        call assert_equals(2,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        nel = 2
        call EncodeBitDet_guga([1,2], ilutI)
        call EncodeBitDet_guga([1,4], ilutJ)
!
        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call assert_equals(5_int_rdm, rdm_ind(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(2,i)
        call assert_equals(1,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, ilutI, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call assert_equals(2_int_rdm, rdm_ind(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(1,i)
        call assert_equals(2,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call EncodeBitDet_guga([1,6], ilutJ)
!
        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(3,i)
        call assert_equals(1,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, ilutI, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(1,i)
        call assert_equals(3,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        nel = 3
        call EncodeBitDet_guga([1,3,4], ilutI)
        call EncodeBitDet_guga([3,4,5], ilutJ)
!
        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(-1.0_dp), mat_ele)
        call assert_equals(-1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(3,i)
        call assert_equals(1,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, ilutI, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(-1.0_dp), mat_ele)
        call assert_equals(-1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(1,i)
        call assert_equals(3,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call EncodeBitDet_guga([1,3,5], ilutI)
        call EncodeBitDet_guga([3,5,7], ilutJ)
!
        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(4,i)
        call assert_equals(1,j)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutj, iluti, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(1.0_dp), mat_ele)
        call assert_equals(1.0_dp, rdm_mat(1))
        call extract_1_rdm_ind(rdm_ind(1), i, j)
        call assert_equals(4,j)
        call assert_equals(1,i)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        nel = 2
        call EncodeBitDet_guga([1,2], ilutI)
        call EncodeBitDet_guga([3,4], ilutJ)
!
        call calc_guga_matrix_element(iluti, ilutj, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(2.0_dp), mat_ele)
        call assert_equals(2.0_dp, rdm_mat(1))
        call extract_2_rdm_ind(rdm_ind(1), i = i, j = j, k = k, l = l)
        call assert_equals(2,i)
        call assert_equals(1,j)
        call assert_equals(2,k)
        call assert_equals(1,l)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call calc_guga_matrix_element(ilutJ, ilutI, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(2.0_dp), mat_ele)
        call assert_equals(2.0_dp, rdm_mat(1))
        call extract_2_rdm_ind(rdm_ind(1), i = i, j = j, k = k, l = l)
        call assert_equals(1,i)
        call assert_equals(2,j)
        call assert_equals(1,k)
        call assert_equals(2,l)
        call assert_equals(1, size(rdm_mat))
        call assert_equals(1, size(rdm_ind))
!
        call EncodeBitDet_guga([3,6], ilutJ)
!
        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
            t_hamil = .false., calc_type = 2, rdm_ind = rdm_ind, rdm_mat = rdm_mat)
        call assert_equals(h_cast(sqrt(2.0_dp)), mat_ele)
        call assert_equals(sqrt(2.0_dp), rdm_mat(1))
        call assert_equals(sqrt(2.0_dp), rdm_mat(2))
        call extract_2_rdm_ind(rdm_ind(1), i = i, j = j, k = k, l = l)
        call assert_equals(3,i)
        call assert_equals(1,j)
        call assert_equals(2,k)
        call assert_equals(1,l)
        call extract_2_rdm_ind(rdm_ind(2), i = i, j = j, k = k, l = l)
        call assert_equals(2,i)
        call assert_equals(1,j)
        call assert_equals(3,k)
        call assert_equals(1,l)
        call assert_equals(2, size(rdm_mat))
        call assert_equals(2, size(rdm_ind))

        nel = 4
!
!
!
    end subroutine test_coupling_coeffs
!
    subroutine run_test_excit_gen_guga(nel_in, nbasis_in, stot_in)
        integer, intent(in) :: nel_in, nbasis_in, stot_in

        nel = nel_in
        nBasis = nbasis_in
        stot = stot_in

        if (tUEG .or. tHUB) then
            if (.not. treal) then
                pSingles = 0.0_dp
                pDoubles = 1.0_dp
            else
                pSingles = 1.0_dp
                pDoubles = 0.0_dp
            end if
        else if (t_heisenberg_model) then
            pSingles = 0.0_dp
            pDoubles = 1.0_dp

        else if (t_tJ_model) then
            pSingles = 0.1_dp
            pDoubles = 1.0 - pSingles

        else if (t_k_space_hubbard) then
            pSingles = 0.0_dp
            pDoubles = 1.0_dp

        else if (t_new_real_space_hubbard) then
            pSingles = 1.0_dp
            pDoubles = 0.0_dp

        else
            pSingles = 0.1_dp
            pDoubles = 1.0_dp - pSingles
        end if

        pExcit4 = 0.5_dp
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

        if (nEl == 4 .and. nBasis/2 == 4 .and. STOT == 2) then
            call run_test_excit_gen_guga_S2

        else if (nEL == 3 .and. nBasis/2 == 4 .and. STOT == 1) then
            call run_test_excit_gen_guga_nEl_3_S_1

        else if (nEl == 2 .and. nBasis/2 == 4) then
            call run_test_excit_gen_guga_nEl_2_S_0

        else if (nEl == 5 .and. nBasis/2 == 4) then
            call run_test_excit_gen_guga_nEl_5_S_1

        else if (nEl == 6 .and. nBasis/2 == 4) then
            call run_test_excit_gen_guga_nEl_6_S_0

        else if (nEl == 6 .and. nBasis/2 == 6 .and. STOT == 0) then
            call run_test_excit_gen_guga_nOrb_6_nEl_6_S_0

        else if (nEl == 6 .and. nBasis/2 == 6 .and. STOT == 2) then
            call run_test_excit_gen_guga_nOrb_6_nEl_6_S_2

        else if (nEl == 6 .and. nBasis/2 == 6 .and. STOT == 4) then
            call run_test_excit_gen_guga_nOrb_6_nEl_6_S_4

        else if (nEl == 4 .and. nBasis/2 == 6 .and. STOT == 0) then
            call run_test_excit_gen_guga_nOrb_6_nEl_4_S_0

        else if (nEl == 4 .and. nBasis/2 == 6 .and. STOT == 2) then
            call run_test_excit_gen_guga_nOrb_6_nEl_4_S_2

        else if (nEl == 5 .and. nBasis/2 == 6 .and. STOT == 1) then
            call run_test_excit_gen_guga_nOrb_6_nEl_5_S_1

        else if (nEl == 5 .and. nBasis/2 == 6 .and. STOT == 3) then
            call run_test_excit_gen_guga_nOrb_6_nEl_5_S_3

        else if (nEl == 7 .and. nBasis/2 == 6 .and. STOT == 1) then
            call run_test_excit_gen_guga_nOrb_6_nEl_7_S_1

        else if (nEl == 7 .and. nBasis/2 == 6 .and. STOT == 3) then
            call run_test_excit_gen_guga_nOrb_6_nEl_7_S_3

        else if (nEl == 5 .and. nBasis/2 == 9 .and. STOT == 1) then
            call run_test_excit_gen_guga_nOrb_9_nEl_5_S_1

        else if (nEl == 5 .and. nBasis/2 == 9 .and. STOT == 3) then
            call run_test_excit_gen_guga_nOrb_9_nEl_5_S_3

        else if (nEl == 7 .and. nBasis/2 == 9 .and. STOT == 1) then
            call run_test_excit_gen_guga_nOrb_9_nEl_7_S_1

        else if (nEl == 7 .and. nBasis/2 == 9 .and. STOT == 3) then
            call run_test_excit_gen_guga_nOrb_9_nEl_7_S_3

        else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 0) then
            call run_test_excit_gen_guga_nOrb_9_nEl_10_S_0

        else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 6) then
            call run_test_excit_gen_guga_nOrb_9_nEl_10_S_6

        else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 4) then
            call run_test_excit_gen_guga_nOrb_9_nEl_10_S_4

        else if (nEl == 10 .and. nBasis/2 == 9 .and. STOT == 2) then
            call run_test_excit_gen_guga_nOrb_9_nEl_10_S_2

        else if (nEl == 9 .and. nBasis/2 == 9 .and. STOT == 3) then
            call run_test_excit_gen_guga_nOrb_9_nEl_9_S_3

        else if (nEl == 9 .and. nBasis/2 == 9 .and. STOT == 1) then
            call run_test_excit_gen_guga_nOrb_9_nEl_9_S_1

        else if (nEl == 18 .and. nBasis/2 == 18 .and. STOT == 0) then
            call run_test_excit_gen_guga_nOrb_18_nel_18_S_0

        else if (nEl == 18 .and. nBasis/2 == 18 .and. STOT == 2) then
            call run_test_excit_gen_guga_nOrb_18_nel_18_S_2

        else if (nEl == 18 .and. nBasis/2 == 18 .and. STOT == 6) then
            call run_test_excit_gen_guga_nOrb_18_nel_18_S_6
        else
            ! todo: create a general excit_gen tester which uses the
            ! guessed HF determinant as a start and uses some of the
            ! exactly created determinants from that to check for pgen and
            ! matrix element consistency!

            call run_test_excit_gen_guga_multiple(&
                [1,2,3,4,5,6,7,8])
            call run_test_excit_gen_guga_multiple(&
                [9,10,11,12,13,14,15,16])
            call run_test_excit_gen_guga_multiple(&
                [1,3,5,7,10,12,14,16])
            call run_test_excit_gen_guga_multiple(&
                [1,3,6,7,10,11,14,16])
            call run_test_excit_gen_guga_multiple(&
                [1,4,5,8,9,12,13,16])
            call run_test_excit_gen_guga_general
            call run_test_excit_gen_guga_multiple([3,7,8,9,10,12,13,14])
            call run_test_excit_gen_guga_multiple([1,4,5,7])
            call run_test_excit_gen_guga_multiple([1,3,6,7])
            call run_test_excit_gen_guga_multiple([1,2,5,8])
            call run_test_excit_gen_guga_multiple([1,3,4,8])

            call run_test_excit_gen_guga_single([1,2,3,5,6,7,8,9,11,12,13,14,15,18])
            call run_test_excit_gen_guga_single([1,2,3,4,5,6,7,8,9,10,11,13,16,18])
            call run_test_excit_gen_guga_single([1,2,3,4,5,7,8,10,11,12,13,14,17,18])
            call run_test_excit_gen_guga_multiple(&
                convert_guga_to_ni([0,0,0,0,1,1,1,2,1,2,3,3,3,1,2],15))
            call run_test_excit_gen_guga_single()
            call run_test_excit_gen_guga_single(&
                convert_guga_to_ni([3,3,3,3,3,1,0,3,1],9))

            call run_test_excit_gen_guga_single(&
                [1,2,3,5,7,10,12,13,14,15,17,18,20,21,25,27,30,32,34,35,36,&
                37,39,41,44,46,48,49,50,51,54,58])

        end if

    end subroutine run_test_excit_gen_guga






    subroutine run_test_excit_gen_guga_general
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_general"
        integer(n_int) :: ilut(0:niftot)
        integer(n_int), pointer :: ex(:,:)
        integer :: nEx, i
        integer :: nTest
        type(ExcitationInformation_t) :: excitInfo

        ! use fdet as first determinant and test on all excitations from this..!
        ! maybe a bit too much for bigger system?
        print *, ""
        print *, "running general test_excit_gen_guga()"
        print *, ""

        ! first act the hamiltonian on the fdet
        call EncodeBitDet(fdet, ilut)

        call actHamiltonian(ilut, ex, nEx)

        if (t_full_guga_tests .or. t_guga_testsuite) then
            ! in this case also check if not too many n_guga_excits are
            ! asked for.. otherwise it takes forever!
            ! if 1B are asked for only do 1 CSF
            if (n_guga_excit_gen < 1000000) then
                ! under a million do all the excits
                nTest = nex
            else if (n_guga_excit_gen < 10000000) then
                ! for a million do 100 excits
                nTest = min(nex,1000)
            else if (n_guga_excit_gen < 100000000) then
                nTest = min(nex,100)
            else if (n_guga_excit_gen < 1000000000) then
                ! for 100 million do only 10
                nTest = min(nex,10)
            else
                ! for more do only the reference!
                nTest = 1
            end if
        else
            nTest = min(nEx,20)
        end if

        print *, ""
        print *, "running tests on nExcits: ", nTest
        print *, ""
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

        print *, ""
        print *, "running tests on nExcits: ", nTest
        print *, ""

        call test_excit_gen_guga(ilut, n_guga_excit_gen)

        do i = 1, nTest
            call test_excit_gen_guga(ex(:,i), n_guga_excit_gen)
        end do

    end subroutine run_test_excit_gen_guga_multiple


    subroutine test_findSwitches
        character(*), parameter :: this_routine = "test_findSwitches"
        integer(n_int) :: ilutI(0:nifguga), ilutJ(0:nifguga)

        print *, ""
        print *, "testing findSwitches routines:"
        print *, ""
        ! 3300
        call EncodeBitDet_guga([1,2,3,4],ilutI)
        ilutJ = ilutI

        call assert_true(findFirstSwitch(ilutI,ilutJ,1,4) == 0)
        call assert_equals(5, findLastSwitch(ilutI,ilutJ,1,4))

        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilutJ)

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
        call assert_true(findFirstSwitch(ilutI,ilutJ,1,4) == 1)
        call assert_true(findLastSwitch(ilutI,ilutJ,1,4) == 4)

        call assert_true(findFirstSwitch(ilutI,ilutJ,2,3) == 2)
        call assert_true(findLastSwitch(ilutI,ilutJ,1,2) == 2)

        ! 1230
        call EncodeBitDet_guga([1,4,5,6],ilutI)
        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilutJ)
!
        ! 1230
        ! 1122
        call assert_true(findFirstSwitch(ilutI,ilutJ,1,3) == 2)
        call assert_true(findFirstSwitch(ilutI,ilutJ,2,3) == 2)
        call assert_true(findLastSwitch(ilutI,ilutJ,1,4) == 4)
        ! for the find last switch we exclude the inputted first orbital!
        call assert_true(findLastSwitch(ilutI,ilutJ,2,3) == 3)
        call assert_true(findFirstSwitch(ilutI,ilutJ,3,4) == 3)
        call assert_true(findLastSwitch(ilutI,ilutJ,3,4) == 4)

        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilutI)
        ! 1212
        call EncodeBitDet_guga([1,4,5,8],ilutJ)

        ! 1122
        ! 1212
        call assert_true(findFirstSwitch(ilutI,ilutJ,2,3) == 2)
        call assert_true(findLastSwitch(ilutI,ilutJ,1,3) == 3)
        call assert_true(findFirstSwitch(ilutI,ilutJ,3,4) == 3)
        call assert_true(findLastSwitch(ilutI,ilutJ,2,3) == 3)

        call EncodeBitDet_guga([1,2,5,6],ilutI)
        call EncodeBitDet_guga([1,2,7,8],ilutJ)

        ! 3030
        ! 3003
        call assert_true(findFirstSwitch(ilutI,ilutJ,1,4) == 3)
        call assert_true(findFirstSwitch(ilutI,ilutJ,1,3) == 0)
        call assert_true(findLastSwitch(ilutI,ilutJ,1,4) == 4)
        call assert_equals(5, (findLastSwitch(ilutI,IlutJ,1,2)))

        call EncodeBitDet_guga([3,4,7,8],ilutI)

        ! 0303
        ! 3003
        call assert_equals(5, findLastSwitch(ilutI,ilutJ,2,4))


        print *, ""
        print *, "findSwitches tests passed!"
        print *, ""

    end subroutine test_findSwitches

    subroutine test_count_beta_orbs_ij
        character(*), parameter :: this_routine = "test_count_beta_orbs_ij"
        integer(n_int) :: ilut(0:nifguga)

        ! these routine now need the current_stepvector quantitiy!

        call EncodeBitDet_guga([1,2,3,4],ilut)

        current_stepvector = calcStepVector(ilut)

        print *, ""
        print *, "testing count_beta_orbs_ij:"
        print *, ""

        call assert_true(count_beta_orbs_ij(ilut,1,4) == 0)
        call assert_true(count_beta_orbs_ij(ilut,1,3) == 0)
        call assert_true(count_beta_orbs_ij(ilut,2,4) == 0)

        call EncodeBitDet_guga([1,3,6,8],ilut)

        current_stepvector = calcStepVector(ilut)

        call assert_true(count_beta_orbs_ij(ilut,1,4) == 2)
        call assert_true(count_beta_orbs_ij(ilut,2,4) == 1)

        call EncodeBitDet_guga([3,5,6,8],ilut)

        current_stepvector = calcStepVector(ilut)

        call assert_true(count_beta_orbs_ij(ilut,1,4) == 1)

        print *, ""
        print *, "count_beta_orbs_ij tests passed!"
        print *, ""

    end subroutine test_count_beta_orbs_ij


    subroutine test_count_alpha_orbs_ij
        character(*), parameter :: this_routine = "test_count_alpha_orbs_ij"
        integer(n_int) :: ilut(0:nifguga)

        ! this routines now need the current_stepvector quantity!
        call EncodeBitDet_guga([1,2,3,4],ilut)

        current_stepvector = calcStepVector(ilut)

        print *, ""
        print *, "testing count_alpha_orbs_ij:"
        print *, ""
        ! 3300
        call assert_true(count_alpha_orbs_ij(ilut,1,4) == 0)
        call assert_true(count_alpha_orbs_ij(ilut,2,3) == 0)

        call EncodeBitDet_guga([1,4,5,6],ilut)

        current_stepvector = calcStepVector(ilut)

        ! 1230
        call assert_true(count_alpha_orbs_ij(ilut,1,4) == 1)
        call EncodeBitDet_guga([3,6,7,8],ilut)

        current_stepvector = calcStepVector(ilut)

        call assert_true(count_alpha_orbs_ij(ilut,1,4) == 1)

        call EncodeBitDet_guga([1,4,5,8],ilut)

        current_stepvector = calcStepVector(ilut)

        call assert_true(count_alpha_orbs_ij(ilut,1,4) == 2)

        print *, ""
        print *, "count_alpha_orbs_ij tests passed!"
        print *, ""

    end subroutine test_count_alpha_orbs_ij

    subroutine test_add_guga_lists
        character(*), parameter :: this_routine = "test_add_guga_lists"
        integer(n_int) :: l1(0:nifguga,1), l2(0:nifguga,1), l3(0:nifguga,2)
        integer(n_int) :: l4(0:nifguga,3), l5(0:nifguga,4), l6(0:nifguga,5)

        integer :: nOut

        l1 = 0_n_int
        l2 = 0_n_int
        l3 = 0_n_int
        l4 = 0_n_int
        l5 = 0_n_int
        l6 = 0_n_int

        call EncodeBitDet_guga([1,2,3,4],l6(:,1))
        call EncodeBitDet_guga([1,4,5,8],l3(:,1))

        print *, "testing: add_guga_lists(n1,n2,l1,l2)"
        print *, ""
        nOut = 1
        call add_guga_lists(nOut, 1, l6, l3)
        call assert_true(nOut == 2)

        call EncodeBitDet_guga([1,2,3,6],l3(:,1))

        call add_guga_lists(nOut, 1, l6, l3)

        call assert_true(nout == 3)

        call EncodeBitDet_guga([5,6,7,8],l3(:,1))

        call add_guga_lists(nOut,1,l6,l3)

        call assert_true(nout == 4)

        l3(:,2) = l6(:,2)
        call write_guga_list(6,l3(:,1:2))
        call write_guga_list(6,l6(:,1:4))

        call add_guga_lists(nOut,2,l6,l3)

        call assert_true(nOut == 4)

        l5(:,1) = l6(:,1)
        nOut = 1

        call add_guga_lists(nOut,4,l5,l6)

        call assert_true(nOut == 4)

        print *, ""
        print *, "add_guga_lists tests passed!"
        print *, ""

    end subroutine test_add_guga_lists

    subroutine test_getSpatialOccupation
        character(*), parameter :: this_routine = "test_getSpatialOccupation"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,6], ilut)

        print *, ""
        print *, "testing getSpatialOccupation(ilut, sOrb):"
        print *, ""
        call assert_true(getSpatialOccupation(ilut,1) .isclose. 2.0_dp)
        call assert_true(getSpatialOccupation(ilut,2) .isclose. 1.0_dp)
        call assert_true(getSpatialOccupation(ilut,3) .isclose. 1.0_dp)
        call assert_true(getSpatialOccupation(ilut,4) .isclose. 0.0_dp)
        print *, ""
        print *, "getSpatialOccupation tests passed!"
        print *, ""

    end subroutine test_getSpatialOccupation

    subroutine test_calcFullStartFullStopMixed
        character(*), parameter :: this_routine = "test_calcfullStartFullStopMixed"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        type(WeightObj_t) :: weights

        call EncodeBitDet_guga([1,4,5,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,1,4)

        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        print *, ""
        print *, "testing calcFullStartFullStopMixed(ilut, exInfo, ex, num, posSwitch, negSwitch):"
        print *, ""

        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 1212
        ! 1122

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,1,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,1,2,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) + 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)


        excitInfo = excitationIdentifier(3,2,2,3)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,1,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,1,2,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) + 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3,1,1,3)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,1,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,1,2,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) + 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) + sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(4,2,2,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,1,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,1,2,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) + 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) + sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)


        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,1,4)

        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1122
        ! 1122
        ! 1212

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,1,2,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,2,1,2]))
        ! -1/2 + 1 -> +1/2
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3,2,2,3)

        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,1,2,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,2,1,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3,1,1,3)

        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,1,2,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,2,1,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) + sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(4,2,2,4)

        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        call calcFullStartFullStopMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,1,2,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,2,1,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 1.0_dp/2.0_dp) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex(:,2),1) + sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartFullStopMixed tests passed!"
        print *, ""

    end subroutine test_calcFullStartFullStopMixed


!
    subroutine run_test_excit_gen_guga_nEl_3_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_3_S_1"
        integer(n_int):: ilut(0:niftot)
        integer :: nI(3)


        print *, ""
        print *, "running: test_excit_gen_guga() for nEl = 3, S = 1"
        print *, ""

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



    end subroutine run_test_excit_gen_guga_nEl_3_S_1

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_2
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_6_S_2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 6, S = 2 system"
        print *, ""

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


    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_2


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_4
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_6_S_4"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 6, S = 4 system"
        print *, ""
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



    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_4


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_2
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_4_S_2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 4, S = 0 system"
        print *, ""

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


    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_2


    subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_5_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 5, S = 1 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_4

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_6
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_10_S_6"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(10)


        print *, ""
        print *, "running: test_excit_gen_guga() on the 9 orbital, Nel = 10, S = 6"
        print *, ""
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

        print *, ""
        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 10, S = 0 system"
        print *, ""
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



    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_10_S_0

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_5_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 5, S = 3 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_5_S_3


    subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_7_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 7, S = 3 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_3

    subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_9_nEl_7_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 9 orbital, nEl = 7, S = 1 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nOrb_9_nEl_7_S_1


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_7_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 7, S = 3 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_3

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_7_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(7)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 7, S = 1 system"
        print *, ""
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



    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_7_S_1


    subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_3
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_5_S_3"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 5, S = 3 system"
        print *, ""
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



    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_3

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_5_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 5, S = 1 system"
        print *, ""

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


    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_5_S_1

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_4_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 4, S = 0 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_4_S_0

    subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nOrb_6_nEl_6_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, ""
        print *, "running: test_excit_gen_guga() on the 6 orbital, nEl = 6, S = 0 system"
        print *, ""

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


    end subroutine run_test_excit_gen_guga_nOrb_6_nEl_6_S_0


    subroutine run_test_excit_gen_guga_nEl_6_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_6_S_0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(6)

        print *, ""
        print *, "running: test_excit_gen_guga() on the nEl = 6, S = 0 system"
        print *, ""
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

    end subroutine run_test_excit_gen_guga_nEl_6_S_0


    subroutine run_test_excit_gen_guga_nEl_5_S_1
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_5_S_1"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(5)


        print *, ""
        print *, "running: test_excit_gen_guga() on the nEl = 5, S = 1 system"
        print *, ""
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



    end subroutine run_test_excit_gen_guga_nEl_5_S_1

    subroutine run_test_excit_gen_guga_nEl_2_S_0
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_nEl_2_S_0"
        integer(n_int):: ilut(0:niftot)
        integer :: nI(2)


        print *, ""
        print *, "running test_excit_gen_guga() gor the nEl = 2, S = 0 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_nEl_2_S_0
!
    subroutine run_test_excit_gen_guga_S2
        ! also check for the S = 2 system...
        ! and probably write this function generally for all sorts of
        ! excitations
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_S2"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)


        print *, ""
        print *, "running: test_excit_gen_guga(ilut,iter) for the S = 2 system"
        print *, ""
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


    end subroutine run_test_excit_gen_guga_S2


    subroutine run_test_excit_gen_guga_S0
        ! write a similar testing routine as Simons
        character(*), parameter :: this_routine = "run_test_excit_gen_guga_S0"
        integer(n_int) :: ilut(0:niftot)
        integer :: nI(4)


        pSingles = 0.1_dp
        pDoubles = 1.0_dp - pSingles

        pExcit4 = 0.5_dp
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

        print *, ""
        print *, "running: test_excit_gen_guga_S0(ilut,n_guga_excit_gen)"
        print *, "pSingles set to: ", pSingles
        print *, "pDoubles set to: ", pDoubles
        print *, ""


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

        print *, ""
        print *, "test_excit_gen_guga finished!"
        print *, ""

    end subroutine run_test_excit_gen_guga_S0

    subroutine test_calcMixedContribution
        character(*), parameter :: this_routine = "test_calcMixedContribution"
        integer(n_int) :: ilut(0:nifguga), t(0:nifguga)

        call EncodeBitDet_guga([1,4,5,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call EncodeBitDet_guga([1,3,6,8], t)

        print *, ""
        print *, "testing calcMixedContribution(ilut,t,start,ende):"
        print *, ""
        ! 1212
        ! 1122
        call assert_true(calcMixedContribution(ilut,t,1,4) .isclose. h_cast(0.0_dp))

        currentB_ilut = calcB_vector_ilut(t)
        currentOcc_ilut = calcOcc_vector_ilut(t)
        currentOcc_int = calcOcc_vector_int(t)
        current_stepvector = calcStepVector(t)
        currentB_int = calcB_vector_int(t)

        call assert_true(calcMixedContribution(t,ilut,1,4) .isclose. h_cast(0.0_dp))

        print *, ""
        print *, "calcMixedContribution tests passed!"
        print *, ""

    end subroutine test_calcMixedContribution

    subroutine test_generate_excitation_guga_double
        character(*), parameter :: this_routine = "test_generate_excitation_guga_double"
        integer :: nI(4), nJ(4), IC, excitMat(2,maxExcit), exFlag, nEx, pos
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer(n_int) :: ilutGi(0:nifguga), ilutGj(0:nifguga)
        logical :: tParity
        real(dp) :: pgen
        HElement_t(dp) :: HElGen
        type(excit_gen_store_type), target :: store
        integer(n_int), pointer :: ex(:,:)

        exFlag = 1
        ! make only double excitations:
        pSingles = 0.0_dp
        pDoubles = 1.0_dp - pSingles

        print *, ""
        print *, "testing generate_excitation_guga:"
        print *, ""
        ! 3300:
        nI = [1,2,3,4]; ilutI = 0_n_int
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        print *, ""
        print *, "random double excitation for :"
        print *, ""
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, ""
        print *, "pgen: ", pgen, "matEle: ", HElGen
        print *, ""

        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, ""
            print *, "exact excitations for this ilut:"
            print *, ""
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))

            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)

        else
            print *, ""
            print *, "no valid excitation created!"
            print *, ""
        end if

        ! 3030
        nI = [1,2,5,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 3003
        nI = [1,2,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > EPS) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        end if

        ! 0330
        nI = [3,4,5,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > EPS) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        end if

        ! 0303
        nI = [3,4,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > EPS) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        end if

        ! 0033
        nI = [5,6,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > EPS) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        end if

        ! 1023
        nI = [1,6,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 3102
        nI = [1,2,3,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 3120
        nI = [1,2,3,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 3012
        nI = [1,2,5,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
               print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 0312
        nI = [3,4,5,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1230
        nI = [1,4,5,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1203
        nI = [1,4,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1320
        nI = [1,3,4,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1302
        nI = [1,3,4,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1032
        nI = [1,5,6,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 0132
        nI = [3,5,6,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 0123
        nI = [3,6,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1122
        nI = [1,3,6,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if

        ! 1212
        nI = [1,4,5,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random double excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)


        else
            print *, "no valid excitation created!"
        end if
        print *, "generate_excitation_guga tests passed!"


    end subroutine test_generate_excitation_guga_double

    subroutine test_generate_excitation_guga_single
        character(*), parameter :: this_routine = "test_generate_excitation_guga_single"
        integer :: nI(4), nJ(4), IC, excitMat(2,maxExcit), exFlag, nEx, pos
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot), ilutGi(0:nifguga)
        integer(n_int) :: ilutGj(0:nifguga)
        logical :: tParity
        real(dp) :: pgen
        HElement_t(dp) :: HElGen
        type(excit_gen_store_type), target :: store
        integer(n_int), pointer :: ex(:,:)

        exFlag = 1
        ! make this store element ...
        print *, ""
        print *, "testing generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,exMat,tPar,pgen,hEl,store)"
        print *, ""
        ! test singles only first
        pSingles = 1.0_dp
        pDoubles = 0.0_dp

        ! 3300:
        nI = [1,2,3,4]; ilutI = 0_n_int
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 3030
        nI = [1,2,5,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 3003
        nI = [1,2,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 0330
        nI = [3,4,5,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 0303
        nI = [3,4,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 0033
        nI = [5,6,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1023
        nI = [1,6,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 3102
        nI = [1,2,3,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 3120
        nI = [1,2,3,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 3012
        nI = [1,2,5,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 0312
        nI = [3,4,5,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1230
        nI = [1,4,5,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1203
        nI = [1,4,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1320
        nI = [1,3,4,6]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1302
        nI = [1,3,4,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1032
        nI = [1,5,6,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 0132
        nI = [3,5,6,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 0123
        nI = [3,6,7,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1122
        nI = [1,3,6,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        ! 1212
        nI = [1,4,5,8]
        call EncodeBitDet(nI,ilutI)

        call init_csf_information(ilutI)

        call generate_excitation_guga(nI,ilutI,nJ,ilutJ,exFlag,IC,excitMat,&
            tParity,pgen,HElGen,store)
        print *, "random single excitation for :"
        call convert_ilut_toGUGA(ilutI, ilutGi)
        call write_det_guga(6, ilutGi)

        print *, "pgen: ", pgen, "matEle: ", HElGen
        call convert_ilut_toGUGA(ilutJ, ilutGj)
        call write_det_guga(6, ilutGj)

        if (pgen > 0.0_dp) then
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilutI, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))
            pos = binary_search(ex(0:nifd,1:nex),ilutJ(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(helgen - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)
        else
            print *, "no valid excitation created!"
        end if

        print *, ""
        print *, "generate_excitation_guga tests passed!"
        print *, ""

    end subroutine test_generate_excitation_guga_single

    subroutine test_calcDoubleR2L_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleR2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,3,4,2)

        call assert_true(excitInfo%typ == 13 )

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        print *, ""
        print *, "testing calcDoubleR2L_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcDoubleR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 3003

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3,0,0,3]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < EPS)

        ! mixed: -1
        ! nonover: +2 -> +2

        ! 0132
        ! 1023
        call EncodeBitDet_guga([3,5,6,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,3,4,2)

        call assert_true(excitInfo%typ == 13 )

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,0,2,3]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp) < 1.0e-10_dp)

        ! nonoverlap : -2
        ! mixed: +1 -> -1
        print *, ""
        print *, "calcDoubleR2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleR2L_stochastic

    subroutine test_calcDoubleL2R_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleL2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(3,1,2,4)

        call assert_true(excitInfo%typ == 12)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcDoubleL2R_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcDoubleL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 0330

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,3,3,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < 1.0e-10_dp)

        ! mixed matele: -1
        ! nonover: 2 -> +1

        ! 3102
        ! 1320
        call EncodeBitDet_guga([1,2,3,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(3,1,2,4)

        call assert_true(excitInfo%typ == 12)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,2,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp) < 1.0e-10_dp)

        ! mixed: -2
        ! nonover: +1 -> -1


        print *, ""
        print *, "calcDoubleL2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleL2R_stochastic

    subroutine test_calcDoubleR2L2R_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleR2L2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,3,2 )

        call assert_true(excitInfo%typ == 11)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcDoubleR2L2R_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcDoubleR2L2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 3030
        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3,0,3,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < EPS)

        ! 0123
        ! 1032
        call EncodeBitDet_guga([3,6,7,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,3,2 )

        call assert_true(excitInfo%typ == 11)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleR2L2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,0,3,2]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp) < EPS)
        ! mixed ele: -1/2 - 3/2 = -2
        ! nonover: 1 -> -1



        ! mixes matele: -1
        ! nonoverlp: 2 -> +1
        print *, ""
        print *, "calcDoubleR2L2R_stochastic tests passed!"
        print *, ""


    end subroutine test_calcDoubleR2L2R_stochastic

    subroutine test_calcDoubleL2R2L_stochastic
        character(*), parameter :: this_routine = "test_calcDoubleL2R2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,3)

        call assert_true(excitInfo%typ == 10)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcDoubleL2R2L_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcDoubleL2R2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 0303

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,3,0,3]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        ! argh i have to consider the non-overlap one..
        ! ok: the first is: -1
        ! the non-overlap: 2 -> so in total: +1
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < 1.0e-10_dp)


        ! 1032
        ! 0123
        call EncodeBitDet_guga([1,5,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,3)

        call assert_true(excitInfo%typ == 10)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleL2R2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,1,2,3]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        ! the mixed is -2
        ! the non-overlap: + 1 -> so -1 in total!
        call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcDoubleL2R2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleL2R2L_stochastic

    subroutine test_calcDoubleRaisingStochastic
        character(*), parameter :: this_routine = "test_calcDoubleRaisingStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,2,3)

        call assert_true(excitInfo%typ == 9)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcDoubleRaisingStochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcDoubleRaisingStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 3300
        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 2.0_dp) < EPS)

        excitInfo = excitationIdentifier(1,3,2,4)

        call assert_true(excitInfo%typ == 9)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        call calcDoubleRaisingStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 3300
        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 2.0_dp) < EPS)

        ! 0132
        ! 1320

        ! encode det
        call EncodeBitDet_guga([3,5,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,2,3)

        call assert_true(excitInfo%typ == 9)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleRaisingStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,2,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(1,3,2,4)

        call assert_true(excitInfo%typ == 9)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleRaisingStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,2,0]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcDoubleRaisingStochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleRaisingStochastic

    subroutine test_calcDoubleLoweringStochastic
        character(*), parameter :: this_routine = "test_calcDoubleLoweringStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        ! 1212
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,3,2)

        call assert_true(excitInfo%typ == 8)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcDoubleLoweringStochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcDoubleLoweringStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 0033
        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,0,3,3]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        ! have to think about the other index comb too!
        call assert_true(abs(extract_matrix_element(ex,1) + 2.0_dp) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(3,2,4,1)

        call assert_true(excitInfo%typ == 8)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        call calcDoubleLoweringStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)
        ! 1212
        ! 0033
        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,0,3,3]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        ! have to think about the other index comb too!
        call assert_true(abs(extract_matrix_element(ex,1) + 2.0_dp) < 1.0e-10_dp)

        ! 3120
        ! 1032
        call EncodeBitDet_guga([1,2,3,6], ilut)

        ! calc b and occ vector

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,3,2)

        call assert_true(excitInfo%typ == 8)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleLoweringStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,0,3,2]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp) < EPS)

        excitInfo = excitationIdentifier(4,2,3,1)

        call assert_true(excitInfo%typ == 8)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcDoubleLoweringStochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,0,3,2]))
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp) < EPS)


        print *, ""
        print *, "calcDoubleLoweringStochastic tests passed!"
        print *, ""

    end subroutine test_calcDoubleLoweringStochastic

    subroutine test_calcFullStopR2L_stochastic
        character(*), parameter :: this_routine = "test_calcFullStopR2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 1,4,4,2  )

        call assert_equals(17, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcFullStopR2L_stochastic(ilut,exinfo,ex,pgen):"
        print *, ""
        call calcFullStopR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 3012

        ! since no switch in the overlap region happended this should be 0

        call assert_true(all(ex == 0_n_int))
        call assert_equals(0.0_dp, pgen)

        ! is there a valid possible?
        ! 1122
        ! 3012
        call EncodeBitDet_guga([1,3,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,4,2)

        call assert_equals(17, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStopR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [3,0,1,2], 4)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,3,3,2)

        call assert_equals(17, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStopR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [3,0,1,2], 4)

        print *, ""
        print *, "calcFullStopR2L_stochastic tests passed!"
        print *, ""


    end subroutine test_calcFullStopR2L_stochastic

    subroutine test_calcFullStopL2R_stochastic
        character(*), parameter :: this_routine = "test_calcFullStopL2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,4 )

        call assert_equals(16, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcFullStopL2R_stochastic(ilut,exInfo,ex,pgen)"
        print *, ""
        call calcFullStopL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 0312
        ! no possible excitation
        call assert_true(all(ex == 0_n_int))
        call assert_equals(0.0_dp, pgen)

        ! also do a possible one..

        ! encode det
        call EncodeBitDet_guga([1,3,6,8], ilut)

        ! calc b and occ vector

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,4 )

        call assert_equals(16, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStopL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1122
        ! 0312

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [0,3,1,2], 4)

        excitInfo = excitationIdentifier(4,1,2,4)

        call assert_equals(16, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStopL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [0,3,1,2], 4)

        print *, ""
        print *, "calcFullStopL2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopL2R_stochastic

    subroutine test_calcFullStartR2L_stochastic
        character(*), parameter :: this_routine = "test_calcFullStartR2L_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,2,4,1)

        call assert_equals(21, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcFullStartR2L_stochastic(ilut,exInfo,ex,pgen):"
        print *, ""
        call calcFullStartR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! also should not yield a valid excitation
        ! 1212
        ! 1023
        call assert_equals(0.0_dp, pgen)
        call assert_true(all(ex == 0_n_int))

        ! 1122
        ! 1203
        call EncodeBitDet_guga([1,3,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,3,4,1)

        call assert_equals(21, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStartR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [1,2,0,3], 4)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(2,3,4,2)

        call assert_equals(21, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStartR2L_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(pgen > EPS)
        call assert_equals(calcStepVector(ex), [1,2,0,3], 4)

        print *, ""
        print *, "calcFullStartR2L_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartR2L_stochastic

    subroutine test_calcFullStartL2R_stochastic
        character(*), parameter :: this_routine = "test_calcFullStartL2R_stochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! encode det
        call EncodeBitDet_guga([1,4,5,8 ], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 1,4,2,1 )

        call assert_equals(20, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcFullStartL2R_stochastic(ilut, exInfo, ex, pgen)"
        print *, ""
        call calcFullStartL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        ! 1212
        ! 1320 -> normal single!
        call assert_equals(0.0_dp, pgen)
        call assert_true(all(ex == 0_n_int))

        ! 1122
        ! 1230 should work!
        call EncodeBitDet_guga([1,3,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 1,4,3,1 )

        call assert_equals(20, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStartL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_equals(calcStepVector(ex), [1,2,3,0],4)
        call assert_true(pgen > EPS)

        ! set up correct excitation information
        excitInfo = excitationIdentifier( 2,4,3,2 )

        call assert_equals(20, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStartL2R_stochastic(ilut,excitInfo,ex,pgen,posSwitches,negSwitches)

        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_equals(calcStepVector(ex), [1,2,3,0],4)
        call assert_true(pgen > EPS)

        print *, ""
        print *, "calcFullStartL2R_stochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartL2R_stochastic

    subroutine test_calcRaisingSemiStopStochastic
        character(*), parameter :: this_routine = "test_calcRaisingSemiStopStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,3,4,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,2,4,1)

        call assert_true(excitInfo%typ == 21)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        ! but switches are not yet set up... wtf
        weights = init_fullStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),&
            currentB_ilut(2))

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "testing calcRaisingSemiStopStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcRaisingSemiStopStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        ! 1302: there should be 2 possible! why not?
        ! 1203.. ah yes.. since no overlap changes..
        call assert_true(pgen < EPS)
        call assert_true(all(ex == 0))

        ! 1122
        ! 1203
        call EncodeBitDet_guga([1,3,6,8],ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(2,3,4,2)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        ! but switches are not yet set up... wtf
        weights = init_fullStartWeight(ilut,3,4,negSwitch(3),posSwitch(3),&
            currentB_ilut(3))

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,2,2,2]))
        call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        call calcRaisingSemiStopStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,2,0,2]))
        call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcRaisingSemiStopStochastic tests passed!"
        print *, ""

    end subroutine test_calcRaisingSemiStopStochastic

    subroutine test_calcLoweringSemiStopStochastic
        character(*), parameter :: this_routine = "test_calcLoweringSemiStopStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,5,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(2,1,1,4)

        call assert_true(excitInfo%typ == 20)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_fullStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),&
            currentB_ilut(2))

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, ""
        print *, "testing calcLoweringSemiStopStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcLoweringSemiStopStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        ! 1032
        ! 1230
        ! but this is again only a single..
        call assert_true(pgen < EPS)
        call assert_true(all(ex == 0_n_int))

        ! do:
        ! 1122
        ! 1230
        call EncodeBitDet_guga([1,3,6,8],ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! set up correct excitation information
        excitInfo = excitationIdentifier(2,4,3,2)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        ! but switches are not yet set up... wtf
        weights = init_fullStartWeight(ilut,3,4,negSwitch(3),posSwitch(3),&
            currentB_ilut(3))

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,2,2,2]))
        call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        call calcLoweringSemiStopStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,2,3,2]))
        call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,2) + sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "calcLoweringSemiStopStochastic tests passed!"
        print *, ""

    end subroutine test_calcLoweringSemiStopStochastic


    subroutine test_calcRaisingSemiStartStochastic
        character(*), parameter :: this_routine = "test_calcRaisingSemiStartStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,5,6,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(4,1,2,4 )
        ! 1032
        ! 0132
        ! is this even compatible??

        call assert_true(excitInfo%typ == 16)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_semiStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),&
            currentB_ilut(2))

        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        ! 1032
        ! 0x
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,0,3,2]))
        call assert_true(abs(extract_matrix_element(ex,1) - 1.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "testing calcRaisingSemiStartStochastic(ilut,exInfo,weigh,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcRaisingSemiStartStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        ! 1032
        ! 01x
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [0,1,3,2]))
        call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "calcRaisingSemiStartStochastic tests passed!"
        print *, ""

    end subroutine test_calcRaisingSemiStartStochastic


    subroutine test_calcLoweringSemiStartStochastic
        character(*), parameter :: this_routine = "test_calcLoweringSemiStartStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: negSwitch(4),posSwitch(4),pgen

        ! encode det
        call EncodeBitDet_guga([1,3,4,8], ilut)

        ! calc b and occ vector
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! set up correct excitation information
        excitInfo = excitationIdentifier(1,4,4,2 )

        call assert_true(excitInfo%typ == 17)

        ! calc the possible switches
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        ! set up correct weights
        weights = init_semiStartWeight(ilut,2,4,negSwitch(2),posSwitch(2),currentB_ilut(2))

        ! modify the excitation so it fits test case:
        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3,3,0,2]))
        call assert_true(abs(extract_matrix_element(ex,1) - Root2) < 1.0e-10_dp)
        print *, ""
        print *, "testing calcLoweringSemiStartStochastic(ilut,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call calcLoweringSemiStartStochastic(ilut,excitInfo,weights,negSwitch,&
            posSwitch,ex,pgen)

        ! 1302
        ! 31x
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [3,1,0,2]))
        call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex,2) + sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)
        print *, ""
        print *, "calcLoweringSemiStartStochastic tests passed!"
        print *, ""

    end subroutine test_calcLoweringSemiStartStochastic

    subroutine test_calcSingleOverlapMixedStochastic
        character(*), parameter :: this_routine = "test_calcSingleOverlapMixedStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! 0330
        call EncodeBitDet_guga([3,4,5,6],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,4,3)

        call assert_equals(7, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        print *, ""
        print *, "testing calcSingleOverlapMixedStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcSingleOverlapMixedStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)

        ! 0330
        ! 1302
        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(all(calcStepVector(ex) == [1,3,0,2]))
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, abs(extract_matrix_element(ex,2)))
        call assert_equals(-Root2, extract_matrix_element(ex,1))

        ! 3003
        call EncodeBitDet_guga([1,2,7,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(3,1,3,4)

        call assert_equals(6, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcSingleOverlapMixedStochastic(ilut,excitInfo, ex, pgen,posSwitches,negSwitches)

        ! 3003
        ! 1032
        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_true(all(calcStepVector(ex) == [1,0,3,2]))
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, abs(extract_matrix_element(ex,2)))
        call assert_equals(Root2, extract_matrix_element(ex,1),1e-10_dp)

        print *, ""
        print *, "calcSingleOverlapMixedStochastic tests passed!"
        print *, ""

    end subroutine test_calcSingleOverlapMixedStochastic


    subroutine test_calcFullStopLoweringStochastic
        character(*), parameter :: this_routine = "test_calcFullStopLoweringStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! 3030
        call EncodeBitDet_guga([1,2,5,6],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,4,3)

        call assert_equals(14, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        print *, ""
        print *, "testing calcFullStopLoweringStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcFullStopLoweringStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)
        ! 3030
        ! 1023
        call print_excitInfo(excitInfo)
        call assert_true(compFlag)
        call assert_true(excitInfo%valid)
        call assert_equals(calcStepVector(ex), [1,0,2,3],4)
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, abs(extract_matrix_element(ex,2)))
        call assert_equals(-Root2, extract_matrix_element(ex,1))

        print *, ""
        print *, "calcFullStopLoweringStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopLoweringStochastic

    subroutine test_calcFullStopRaisingStochastic
        character(*), parameter :: this_routine = "test_calcFullStopRaisingStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! 0303
        call EncodeBitDet_guga([3,4,7,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(1,4,3,4)

        call assert_equals(15, excitInfo%typ)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        print *, ""
        print *, "testing calcFullStopRaisingStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcFullStopRaisingStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)

        ! 0303
        ! 1320
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [1,3,2,0]))
        call assert_equals(1.0_dp, pgen)
        call assert_equals(0.0_dp, (extract_matrix_element(ex,2)))
        call assert_equals(-Root2, extract_matrix_element(ex,1))

        print *, ""
        print *, "calcFullStopRaisingStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStopRaisingStochastic


    subroutine test_calcFullStartLoweringStochastic
        character(*), parameter :: this_routine = "test_calcFullStartLoweringStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! 3030
        call EncodeBitDet_guga([1,2,5,6],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(3,1,4,1)

        call assert_true(excitInfo%typ == 18)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        call assert_true(.not.compFlag)
        print *, ""
        print *, "testing calcFullStartLoweringStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcFullStartLoweringStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)

        excitInfo = excitationIdentifier(2,1,4,1)

        call assert_true(excitInfo%typ == 18)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call assert_true(compFlag)
        call calcFullStartLoweringStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)
        ! 3030
        ! 0132
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [0,1,3,2]))
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + Root2) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartLoweringStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartLoweringStochastic


    subroutine test_calcFullStartRaisingStochastic
        character(*), parameter :: this_routine = "test_calcFullStartRaisingStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! 0033
        call EncodeBitDet_guga([5,6,7,8],ilut)
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,1,4)

        call assert_true(excitInfo%typ == 19)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        print *, ""
        print *, "testing calcFullStartRaisingStochastic(ilut, exInfo, ex, pgen):"
        print *, ""
        call calcFullStartRaisingStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)

        ! only result is: pgen should be 1..
        ! 0033
        ! 3012
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [3,0,1,2]))
        call assert_true(pgen .isclose. 1.0_dp)
        ! umat is also stored in there.. so i hope i get it right
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + Root2) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(2,3,2,4)

        call assert_true(excitInfo%typ == 19)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)

        call calcFullStartRaisingStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)
        ! 0033
        ! 0312
        call assert_true(compFlag)
        call assert_true(all(calcStepVector(ex) == [0,3,1,2]))
        call assert_true(pgen .isclose. 1.0_dp)
        ! umat is also stored in there.. so i hope i get it right
        call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        call assert_true(abs(extract_matrix_element(ex,1) + Root2) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartRaisingStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartRaisingStochastic

    subroutine test_mixedFullStopStochastic
        character(*), parameter :: this_routine = "test_mixedFullStopStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        call EncodeBitDet_guga([1,4,5,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,1,4)
        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call doubleUpdateStochastic(ilut,2,excitInfo,weights,negSwitch,posSwitch,ex,pgen)
        call doubleUpdateStochastic(ilut,3,excitInfo,weights,negSwitch,posSwitch,ex,pgen)

        ! i should never get the other matrix element.. due to the 0
        ! matrix element or?? hopefully!
        ! no! it is not 0!
        print *, ""
        print *, "testing mixedFullStopStochastic(ilut, excitInfo, ex)"
        print *, ""
        call mixedFullStopStochastic(ilut, excitInfo, ex)

        if (isOne(ex,3)) then
            call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
            call assert_true(abs(extract_matrix_element(ex,1) + 1.0_dp/2.0_dp) < 1.0e-10_dp)
        else if (isTwo(ex,3)) then
            call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)
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

        ! 1212
        call EncodeBitDet_guga([1,4,5,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,1,4)
        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        ! 1212
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,2,1,2]))
        call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)

        print *, ""
        print *, "testing doubleUpdateStochastic(ilut,orb,exInfo,weight,negSwitch,posSwitch,ex,pgen):"
        print *, ""
        call doubleUpdateStochastic(ilut,2,excitInfo,weights,negSwitch,posSwitch,ex,pgen)

        ! now there are 2 possibs.
        ! although.. do i exclude the "diagonal" excitation??
        ! because if yes, then there is only one possib here..
        call assert_true(pgen < 1.0_dp)
        if (isTwo(ex,2)) then
            call assert_true(all(calcStepVector(ex) == [1,2,1,2]))
            call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
            call assert_true(abs(extract_matrix_element(ex,2)) < EPS)
        else if (isOne(ex,2)) then
            call assert_true(all(calcStepVector(ex) == [1,1,1,2]))
            call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex,2) + sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        else
            call stop_all(this_routine, "wrong stepvalue")
        end if

        call doubleUpdateStochastic(ilut,3,excitInfo,weights,negSwitch,posSwitch,ex,pgen)

        call assert_true(pgen < 1.0_dp)

        ! 121
        if (isOne(ex,3)) then
            call assert_true(all(calcStepVector(ex) == [1,2,1,2]))
            call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
            call assert_true(abs(extract_matrix_element(ex,2)) < EPS)

        else if (isTwo(ex,3)) then
            call assert_true(all(calcStepVector(ex) == [1,1,2,2]))
            call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex,2) - OverR2) < 1.0e-10_dp)

        else
            call stop_all(this_routine, "wrong stepvalue!")

        end if


        print *, ""
        print *, "doubleUpdateStochastic tests passed!"
        print *, ""

    end subroutine test_doubleUpdateStochastic

    subroutine test_calcFullStartFullStopMixedStochastic
        character(*), parameter :: this_routine = "test_calcFullStartFullStopMixedStochastic"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ! set up determinant and excitaiton information
        call EncodeBitDet_guga([1,4,5,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(1,4,4,1)

        call assert_true(excitInfo%typ==23)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        print *, ""
        print *, "testing calcFullStartFullStopMixedStochastic(ilut,exInfo,ex,pgen)"
        print *, ""
        call calcFullStartFullStopMixedStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)

        ! in this constellation no excitaiton should be possible, due to 0
        ! matrix elements..
        call assert_true(all(ex == 0) .or. all(calcStepVector(ex) == [1,1,2,2]))

        ! 1122
        call EncodeBitDet_guga([1,3,6,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier(1,4,4,1)

        call assert_true(excitInfo%typ==23)

        call checkCompatibility(ilut,excitInfo,compFlag,posSwitches,negSwitches)
        print *, ""
        print *, "testing calcFullStartFullStopMixedStochastic(ilut,exInfo,ex,pgen)"
        print *, ""
        call calcFullStartFullStopMixedStochastic(ilut, excitInfo, ex, pgen,posSwitches,negSwitches)

        ! in this constellation no excitaiton should be possible, due to 0
        ! matrix elements..
        call assert_true(all(ex == 0) .or. all(calcStepVector(ex) == [1,2,1,2]))

        print *, ""
        print *, "calcFullStartFullStopMixedStochastic tests passed!"
        print *, ""

    end subroutine test_calcFullStartFullStopMixedStochastic


    subroutine test_pickOrbitals_double
        character(*), parameter :: this_routine = "test_pickOrbitals_double"
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int) :: ilut(0:nifguga)
        real(dp) :: pgen
        integer :: nI(4)

        nI = [1,2,3,4]

        call EncodeBitDet_guga([1,2,3,4],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! still have to think about, if i preemptively choose the different
        ! types of double excitation (iiij, ii,jk, etc.) or let i happen
        ! randomly and adjust the pgens accordingly...
        ! and what should i test here??

        print *, ""
        print *, "testing pickOrbitals_double(ilut, excitLvl):"
        print *, ""
        ! 3300
        call pickOrbitals_double(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            ! what can i test here?
            ! only lowerings possible..
            call print_excitInfo(excitInfo)
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(-1, excitInfo%gen1)
            call assert_equals(-1, excitInfo%gen2)
        end if

        call pickOrbitals_double(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            ! what can i test here?
            ! only lowerings possible..
            call print_excitInfo(excitInfo)
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(-1, excitInfo%gen1)
            call assert_equals(-1, excitInfo%gen2)
        end if

        nI = [5,6,7,8]
        call EncodeBitDet_guga(nI,ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! 0033
        call pickOrbitals_double(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            call print_excitInfo(excitInfo)
            ! what can i test here?
            ! only lowerings possible..
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(1, excitInfo%gen1)
            call assert_equals(1, excitInfo%gen2)
        end if

        call pickOrbitals_double(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            call print_excitInfo(excitInfo)
            ! what can i test here?
            ! only lowerings possible..
            call assert_true(pgen > EPS)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_equals(1, excitInfo%gen1)
            call assert_equals(1, excitInfo%gen2)
        end if

        nI = [1,4,5,8]

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call pickOrbitals_double(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > EPS)
        end if
        print *, ""
        print *, "pickOrbitals_double tests passed!"
        print *, ""

    end subroutine test_pickOrbitals_double

    subroutine test_createStochasticExcitation_double
        character(*), parameter :: this_routine = "test_createStochasticExcitation_double"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        real(dp) :: pgen
        integer :: dummy(2), nI(4), pos, nex
        integer(n_int), pointer :: all_ex(:,:)

        nI = [1,5,6,8]

        call EncodeBitDet_guga([1,5,6,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        print *, ""
        print *, "testing createStochasticExcitation_double(ilut, ex, pgen):"
        print *, ""
        call createStochasticExcitation_double(ilut,nI,ex,pgen,dummy)

        ! what should i test here?
        if (pgen > EPS) then
            call actHamiltonian(ilut,all_ex,nex)

            pos = binary_search(all_ex(0:nifd,1:nex),ex(0:nifd))

            call assert_true(pos > 0)
            call assert_true(abs(extract_matrix_element(all_ex(:,pos),1) - extract_matrix_element(ex,1)) < 1.0e-10_dp)

        end if

        print *, ""
        print *, "createStochasticExcitation_double tests passed!"
        print *, ""

    end subroutine test_createStochasticExcitation_double

    subroutine test_singleStochasticEnd
        character(*), parameter :: this_routine = "test_singleStochasticEnd"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        ! 3300
        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        ! do not make it random here but choose specific excitation!
!         excitInfo = pickOrbitals_single(ilut)
        excitInfo = excitationIdentifier(4,1)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, excitInfo%fullEnd)

        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call singleStochasticUpdate(ilut, 2, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call singleStochasticUpdate(ilut, 3, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        print *, ""
        print *, "testing singleStochasticEnd(excitInfo, excitation):"
        print *, ""

        call singleStochasticEnd(excitInfo, ex)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,0,2]))
        call assert_true(abs(extract_matrix_element(ex,1) + Root2) < 1.0e-10_dp)


        print *, ""
        print *, "singleStochasticEnd tests passed!"
        print *, ""

    end subroutine test_singleStochasticEnd


    subroutine test_singleStochasticUpdate
        character(*), parameter :: this_routine = "test_singleStochasticUpdate"
        integer(n_int) :: ilut(0:nifguga), ex(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), pgen

        ! 3300
        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! do not make it random here but choose specific excitation!
!         excitInfo = pickOrbitals_single(ilut)
        excitInfo = excitationIdentifier(4,1)

        call assert_true(excitInfo%typ == 0)
        call assert_true(excitInfo%fullStart == 1 .and. excitInfo%fullEnd == 4)
        call assert_true(excitInfo%gen1 == -1)

        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, excitInfo%fullEnd)

        call assert_true( all(posSwitch < EPS))
        call assert_true( all(negSwitch < EPS))

        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex,1) - Root2) < 1.0e-10_dp)

        print *, ""
        print *, "testing singleStochasticUpdate(ilut, exInfo, weight, posSwitch, negSwitch, ex, pgen):"
        print *, ""
        call singleStochasticUpdate(ilut, 2, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex,1) + Root2) < 1.0e-10_dp)

        call singleStochasticUpdate(ilut, 3, excitInfo, weights, posSwitch, &
            negSwitch, ex, pgen)

        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex,1) + Root2) < 1.0e-10_dp)


        print *, ""
        print *, "singleStochasticUpdate tests passed!"
        print *, ""

    end subroutine test_singleStochasticUpdate

    subroutine test_pickRandomOrb
        character(*), parameter :: this_routine = "test_pickRandomOrb"
        integer :: orb
        real(dp) :: pgen
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        pgen = 1.0_dp

        call pickRandomOrb_scalar(1, pgen, orb)
        call assert_true(pgen .isclose. 1.0_dp/3.0_dp)
        call assert_true( orb > 1 .and. orb <= 4)
        pgen = 1.0_dp
        call pickRandomOrb_forced(2, pgen, orb)
        call assert_true(orb == 1)
        call assert_true(pgen .isclose. 1.0_dp)
        pgen = 1.0_dp
        call pickRandomOrb_vector([1,2], pgen, orb)
        call assert_true(orb == 3 .or. orb == 4)
        call assert_true(pgen .isclose. 1.0_dp/2.0_dp)
        pgen = 1.0_dp
        call pickRandomOrb_vector([2,3], pgen, orb, 2)
        call assert_true(orb == 4)
        call assert_true(pgen .isclose. 1.0_dp)

        pgen = 1.0_dp
        call pickRandomOrb_scalar(0, pgen, orb, 0)
        call assert_true(pgen .isclose. 1.0_dp/3.0_dp)
        call assert_true(orb /= 3)

        pgen = 1.0_dp
        call pickRandomOrb_scalar(2, pgen, orb, 2)
        call assert_true(pgen .isclose. 1.0_dp/2.0_dp)
        call assert_true(orb == 3 .or. orb == 4)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,2,pgen,orb)
        call assert_true(pgen .isclose. 0.0_dp)
        call assert_true(orb == 0)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,4,pgen,orb)
        call assert_true(pgen .isclose. 1.0_dp/2.0_dp)
        call assert_true(orb == 2 .or. orb == 3)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,4,pgen,orb,0)
        call assert_true(pgen .isclose. 1.0_dp)
        call assert_true(orb == 2)

        pgen = 1.0_dp
        call pickRandomOrb_restricted(1,4,pgen,orb,1)
        call assert_true(pgen .isclose. 1.0_dp)
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

        ! 1230
        call EncodeBitDet_guga([1,4,5,6],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,1,3)

        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        print *, ""
        print *, "testing mixedFullStartStochastic:"
        print *, ""

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, prob)

        call assert_true(prob .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,2,3,0]))
        call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)

        excitInfo = excitationIdentifier(4,2,2,3)

        weights = init_doubleWeight(ilut, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        call mixedFullStartStochastic(ilut, excitInfo, weights, posSwitch, &
            negSwitch, ex, prob)

        ! i think i have two possibs here..
        ! 1230
        ! 1212
        ! 1122
        call assert_true(prob < 1.0_dp)

        if (isTwo(ex,2)) then
            call assert_true(all(calcStepVector(ex) == [1,2,3,0]))
            call assert_true(abs(extract_matrix_element(ex,1) + OverR2) < 1.0e-10_dp)
            call assert_true(abs(extract_matrix_element(ex,2)) < EPS)

        else if (isOne(ex,2)) then
            call assert_true(all(calcStepVector(ex) == [1,1,3,0]))
            call assert_true(abs(extract_matrix_element(ex,1)) < EPS)
            call assert_true(abs(extract_matrix_element(ex,2) - sqrt(3.0_dp)/2.0_dp) < 1.0e-10_dp)

        else
            call stop_all(this_routine, "wrong stepvalue at fullstart!")
        endif

        print *, ""
        print *, "mixedFullStartStochastic tests passed!"
        print *, ""

    end subroutine test_mixedFullStartStochastic

    subroutine test_createStochasticStart_single
        character(*), parameter :: this_routine = "test_createStochasticStart_single"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitch(4), negSwitch(4), probWeight
        integer(n_int) :: ex(0:nifguga)

        ! 3300
        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(4,1)

        call assert_true(excitInfo%typ == 0)
        call assert_true(excitInfo%fullStart == 1 .and. excitInfo%fullEnd == 4)
        call assert_true(excitInfo%gen1 == -1)

        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut,4)

        call assert_true( all(posSwitch < EPS))
        call assert_true( all(negSwitch < EPS))

        print *, ""
        print *, "testing createStochasticStart_single(ilut,exInfo, weighs, posSwitch, negSwitch, ex, probWeight):"
        print *, ""
        call createStochasticStart_single(ilut, excitInfo, weights, posSwitch, negSwitch, ex, probWeight)

        ! i should check the matrix element and the excitation to be sure
        ! about the effect!
        call assert_true(probWeight .isclose. 1.0_dp)
        call assert_true(all(calcStepVector(ex) == [1,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex,1) - Root2) < 1.0e-10_dp)

        print *, ""
        print *, "createStochasticStart_single tests passed!"
        print *, ""

    end subroutine test_createStochasticStart_single

    subroutine test_pickOrbitals_single
        character(*), parameter :: this_routine = "test_pickOrbitals_single"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        integer :: nI(4)

        nI = [1,2,3,4]

        call EncodeBitDet_guga([1,2,3,4],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        ! what should i test here?..
        print *, ""
        print *, "testing: pickOrbitals_single(ilut)"
        print *, ""

        ! 3300
        call pickOrbitals_single(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == -1)
        end if

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == -1)
        end if

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == -1)
        end if

        ! 0033
        nI = [5,6,7,8]
        call EncodeBitDet_guga([5,6,7,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == 1)
        end if

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == 1)
        end if

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)

        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
            call assert_true(excitInfo%fullstart == 1 .or. excitInfo%fullstart == 2)
            call assert_true(excitInfo%fullEnd == 3 .or. excitInfo%fullEnd == 4)
            call assert_true(excitInfo%gen1 == 1)
        end if

        call EncodeBitDet_guga([1,4,5,8], ilut)

        nI = [1,4,5,8]
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        ! 1212
        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
        end if

        call pickOrbitals_single(ilut, nI, excitInfo, pgen)
        if (excitInfo%valid) then
            call assert_true(pgen > 0.0_dp)
            call assert_true(excitInfo%typ == 0)
        end if

        print *, ""
        print *, "pickOrbitals_single tests passed!"
        print *, ""

    end subroutine test_pickOrbitals_single

    subroutine test_createStochasticExcitation_single
        character(*), parameter :: this_routine = "test_createStochasticExcitation_single"
        integer(n_int) :: ilut(0:nifguga), t(0:nifguga)
        real(dp) :: pgen
        integer :: nI(4), pos, nex
        HElement_t(dp) :: HElGen
        integer(n_int), pointer :: ex(:,:)

        nI = [1,2,3,4]

        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        print *, ""
        print *, "testing: createStochasticExcitation_single(ilut,t,weight):"
        print *, ""
        call createStochasticExcitation_single(ilut, nI, t, pgen)

        if (pgen > 0.0_dp) then
            print *, "stochastic excitation: "
            call write_det_guga(6,t,.true.)
            print *, "exact excitations for this ilut:"
            call actHamiltonian(ilut, ex, nEx)
            call write_guga_list(6, ex(:,1:nEx))

            pos = binary_search(ex(0:nifd,1:nex),t(0:nifd))
            call assert_true(pos > 0)
            call assert_true(abs(extract_matrix_element(t,1) - extract_matrix_element(ex(:,pos),1)) < 1.0e-10_dp)

        else
            print *, "no valid excitation created!"
        end if

        print *, "createStochasticExcitation_single tests passed!"

    end subroutine test_createStochasticExcitation_single


    subroutine test_actHamiltonian
        character(*), parameter :: this_routine = "test_actHamiltonian"
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), pointer :: ex(:,:)
        integer :: nEx


        nel = 4

        print *, ""
        print *, "testing actHamiltonian(ilut):"
        print *, ""
        ! 3300:
        call EncodeBitDet_guga([1,2,3,4], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(13, nEx)
        ! 0330
        call EncodeBitDet_guga([3,4,5,6],ilut)
        call actHamiltonian(ilut,ex,nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(14, nEx)
         ! 0303
        call EncodeBitDet_guga([3,4,7,8],ilut)
        call actHamiltonian(ilut,ex,nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(14, nEx)
        ! 0033
        call EncodeBitDet_guga([5,6,7,8],ilut)
        call actHamiltonian(ilut,ex,nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(13, nEx)
       ! 1023
        call EncodeBitDet_guga([1,6,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
       ! 3102
        call EncodeBitDet_guga([1,2,3,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
       ! 3120
        call EncodeBitDet_guga([1,2,3,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)

        ! 3030
        call EncodeBitDet_guga([1,2,5,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(14, nEx)
        ! 3003:
        call EncodeBitDet_guga([1,2,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(14, nEx)
        ! 3012
        call EncodeBitDet_guga([1,2,5,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(16, nEx)
        ! 0312
        call EncodeBitDet_guga([3,4,5,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(16, nEx)
        ! 1230
        call EncodeBitDet_guga([1,4,5,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(16, nEx)
        ! 1203
        call EncodeBitDet_guga([1,4,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(16, nEx)
        ! 1320
        call EncodeBitDet_guga([1,3,4,6], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
        ! 1302
        call EncodeBitDet_guga([1,3,4,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
        ! 1032
        call EncodeBitDet_guga([1,5,6,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
        ! 0132
        call EncodeBitDet_guga([3,5,6,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
        ! 0123
        call EncodeBitDet_guga([3,6,7,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(17, nEx)
              ! 1122
        call EncodeBitDet_guga([1,3,6,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(12,nEx)
        ! 1212
        call EncodeBitDet_guga([1,4,5,8], ilut)
        call actHamiltonian(ilut, ex, nEx)
        print *, "number of excitations for: ", nEx
        call write_det_guga(6,ilut)
        call write_guga_list(6,ex(:,1:nEx))
        call assert_equals(18, nEx)

        print *, ""
        print *, "actHamiltonian tests passed!"
        print *, ""

    end subroutine test_actHamiltonian


    subroutine test_calcAllExcitations_double
        character(*), parameter :: this_routine = "test_calcAllExcitations_double"
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), pointer :: ex(:,:)
        integer :: nExcits

        call EncodeBitDet_guga([1,4,5,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        print *, ""
        print *, "testing calcAllExcitations_double(ilut,i,j,k,l,ex,nExits):"
        print *, ""
        call calcAllExcitations_double(ilut,1,2,3,4, ex, nExcits)

        ! meh... was soll ich hier testen?
        ! 1212
        ! 3030
        call assert_true(nExcits == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [3,0,3,0]))
        print *, ""
        print *, "calcAllExcitations_double tests passed!"
        print *, ""

    end subroutine test_calcAllExcitations_double


    subroutine test_calcFullStartFullStopAlike
        character(*), parameter :: this_routine = "test_calcFullStartFullStopAlike"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(1,4,1,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 22)

        print *, ""
        print *, "testing: calcFullStartFullStopAlike(ilut, exInfo, ex)"
        print *, ""
        call calcFullStartFullStopAlike(ilut, excitInfo, ex)

        ! 0303
        ! 3300
        call assert_true(all(calcStepVector(ex(:,1)) == [3,3,0,0]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartFullStopAlike tests passed!"
        print *, ""

        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(4,1,4,1)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        currentB_ilut = calcB_vector_ilut(ilut)
        call assert_true(excitInfo%typ == 22)

        print *, ""
        print *, "testing: calcFullStartFullStopAlike(ilut, exInfo, ex)"
        print *, ""
        call calcFullStartFullStopAlike(ilut, excitInfo, ex)

        ! 3300
        ! 0303
        call assert_true(all(calcStepVector(ex(:,1)) == [0,3,0,3]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 2.0_dp) < 1.0e-10_dp)

        print *, ""
        print *, "calcFullStartFullStopAlike tests passed!"
        print *, ""


    end subroutine test_calcFullStartFullStopAlike

    subroutine test_calcFullStartL2R
        character(*), parameter :: this_routine = "test_calcFullStartL2R"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        excitInfo = excitationIdentifier(1,4,3,1)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 20)

        print *, ""
        print *, "testing: calcFullStartL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 1230

        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,3,0]))

        excitInfo = excitationIdentifier(2,4,3,2)
        call calcFullStartL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,3,0]))
        print *, ""
        print *, "calcFullStartL2R tests passed!"
        print *, ""

    end subroutine test_calcFullStartL2R


    subroutine test_calcFullStartR2L
        character(*), parameter :: this_routine = "test_calcFullStartR2L"
        integer(n_int) :: ilut(0:2)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(1,3,4,1)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 21)

        print *, ""
        print *, "testing: calcFullStartR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 1203
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,0,3]))

        excitInfo = excitationIdentifier(2,3,4,2)
        call calcFullStartR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,0,3]))

        print *, ""
        print *, "calcFullStartR2L tests passed!"
        print *, ""

    end subroutine test_calcFullStartR2L


    subroutine test_calcFullStartRaising
        character(*), parameter :: this_routine = "test_calcFullStartRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(1,4,1,3)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 19)

        print *, ""
        print *, "testing: calcFullStartRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0312
        ! 3300

        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [3,3,0,0]))

        print *, ""
        print *, "calcFullStartRaising tests passed!"
        print *, ""

    end subroutine test_calcFullStartRaising

    subroutine test_calcFullStartLowering
        character(*), parameter :: this_routine = "test_calcFullStartLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)

        excitInfo = excitationIdentifier(4,1,3,1)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 18)

        print *, ""
        print *, "testing: calcFullStartLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStartLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3012
        ! 0033
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [0,0,3,3]))

        print *, ""
        print *, "calcFullStartLowering tests passed!"
        print *, ""

    end subroutine test_calcFullStartLowering


       subroutine test_calcFullStopR2L
        character(*), parameter :: this_routine = "test_calcFullStopR2L"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        currentB_ilut = calcB_vector_ilut(ilut)


        excitInfo = excitationIdentifier(1,4,4,3)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 17)

        print *, ""
        print *, "testing: calcFullStopR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 3102
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [3,1,0,2]))

        excitInfo = excitationIdentifier(1,3,3,2)
        call calcFullStopR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)
        ! 1212
        ! 3012
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [3,0,1,2]))

        print *, ""
        print *, "calcFullStopR2L tests passed!"
        print *, ""

    end subroutine test_calcFullStopR2L

    subroutine test_calcFullStopL2R
        character(*), parameter :: this_routine = "test_calcFullStopL2R"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(4,1,3,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 16)

        print *, ""
        print *, "testing: calcFullStopL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 0132
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [0,1,3,2]))

        excitInfo = excitationIdentifier(3,1,2,3)
        call calcFullStopL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 0312
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [0,3,1,2]))

        print *, ""
        print *, "calcFullStopL2R tests passed!"
        print *, ""

    end subroutine test_calcFullStopL2R


    subroutine test_calcFullStopRaising
        character(*), parameter :: this_routine = "test_calcFullStopRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([5,6,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,4,2,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 15)

        print *, ""
        print *, "testing: calcFullStopRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0033
        ! 1230
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,3,0]))

        print *, ""
        print *, "calcFullStopRaising tests passed!"
        print *, ""

    end subroutine test_calcFullStopRaising

    subroutine test_calcFullStopLowering
        character(*), parameter :: this_routine = "test_calcDoubleFullStopLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(4,1,4,2)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 14)

        print *, ""
        print *, "testing: calcFullStopLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcFullStopLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3300
        ! 1203
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,0,3]))

        print *, ""
        print *, "calcFullStopLowering tests passed!"
        print *, ""

    end subroutine test_calcFullStopLowering

    subroutine test_calcDoubleR2L
        character(*), parameter :: this_routine = "test_calcDoubleR2L"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,4,2)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 13)

        print *, ""
        print *, "testing: calcDoubleR2L(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleR2L(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 3003
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [3,0,0,3]))

        print *, ""
        print *, " numExcits: ", num
        print *, ""
        print *, ex
        print *, ""
        print *, "calcDoubleR2L tests passed!"
        print *, ""

    end subroutine test_calcDoubleR2L

    subroutine test_calcDoubleL2R
        character(*), parameter :: this_routine = "test_calcDoubleL2R"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,4,5,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(3,1,2,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 12)

        print *, ""
        print *, "testing: calcDoubleL2R(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleL2R(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 1212
        ! 0330
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [0,3,3,0]))

        print *, ""
        print *, "calcDoubleL2R tests passed!"
        print *, ""

    end subroutine test_calcDoubleL2R


    subroutine test_calcDoubleRaising
        character(*), parameter :: this_routine = "test_calcDoubleRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([5,6,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,2,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 9)

        print *, ""
        print *, "testing: calcDoubleRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0033
        ! 1122
        ! 1212

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,1,2,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,2,1,2]))

        print *, "calcDoubleRaising tests passed!"

    end subroutine test_calcDoubleRaising

    subroutine test_calcDoubleLowering
        character(*), parameter :: this_routine = "test_calcDoubleLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(3,1,4,2)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        call assert_true(excitInfo%typ==8)
        print *, ""
        print *, "testing: calcDoubleLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        print *, ""
        call calcDoubleLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3300
        ! 1122
        ! 1212
        call assert_true(num == 2)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,1,2,2]))
        call assert_true(all(calcStepVector(ex(:,2)) == [1,2,1,2]))

        print *, ""
        print *, "calcDoubleLowering tests passed!"
        print *, ""

    end subroutine test_calcDoubleLowering

    subroutine test_calcSingleOverlapRaising
        character(*), parameter :: this_routine = "test_singleOverlapRaising"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([5,6,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,3,3,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapRaising(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapRaising(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0033
        ! 1032
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,0,3,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - Root2) < 1.0e-10_dp)

        print *, "calcSingleOverlapRaising tests passed!"

    end subroutine test_calcSingleOverlapRaising

    subroutine test_calcSingleOverlapMixed
        character(*), parameter :: this_routine = "test_singleOverlapMixed"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(3,1,3,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)

        print *, "testing: calcSingleOverlapMixed(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapMixed(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3003
        ! 1032
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,0,3,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - Root2) < 1.0e-10_dp)

        print *, "calcSingleOverlapMixed tests passed!"

    end subroutine test_calcSingleOverlapMixed

    subroutine test_calcSingleOverlapLowering
        character(*), parameter :: this_routine = "test_singleOverlapLowering"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)


        excitInfo = excitationIdentifier(2,1,4,2)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        print *, excitInfo%typ

        print *, "testing: calcSingleOverlapLowering(ilut, exInfo, ex, num, posSwitch, negSwitch)"
        call calcSingleOverlapLowering(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 3300
        ! 1302
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,3,0,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - Root2) < 1.0e-10_dp)

        print *, "calcSingleOverlapLowering tests passed!"

    end subroutine test_calcSingleOverlapLowering

    subroutine test_calcNonOverlapDouble
        character(*), parameter :: this_routine = "test_calcNonOverlapDouble"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)

        call EncodeBitDet_guga([3,4,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1,2,3,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 3)

        print *, "testing: calcNonOverlapDouble(ilut, exInfo, exs, num, posSwitch, negSwitch"
        call calcNonOverlapDouble(ilut, excitInfo, ex, num, posSwitch, negSwitch)

        ! 0303
        ! 1212

        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,2,1,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) - 2.0_dp) < 1.0e-10_dp)

        ! todo: auch hier eine funktion noetig die nur die single excitations
        ! berechnt und nicht tmat einbezieht... und das auch effektiv macht

        print *, "calcNonOverlapDouble tests passed!"

    end subroutine test_calcNonOverlapDouble


    subroutine test_calcDoubleExcitation_withWeight
        character(*), parameter :: this_routine = "test_calcDoubleExcitation_withWeight"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), pointer :: ex(:,:)
        integer :: num
        real(dp) :: posSwitch(4), negSwitch(4)
        logical :: compFlag

        call EncodeBitDet_guga([3,4,7,8],ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(2,2,1,4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        call assert_true(excitInfo%typ == 1)

        print *, "testing: calcDoubleExcitation_withWeight(ilut, exInfo, exc, num)"
        call calcDoubleExcitation_withWeight(ilut, excitInfo, ex, num, posSwitch, &
            negSwitch)

        ! 0303
        ! 1302
        call assert_true(num == 1)
        call assert_true(all(calcStepVector(ex(:,1)) == [1,3,0,2]))
        call assert_true(abs(extract_matrix_element(ex(:,1),1) + 2.0_dp*Root2) < 1.0e-10_dp)


        excitInfo = excitationIdentifier(3,3,1,4)
        call checkCompatibility(ilut,excitInfo,compFlag,posSwitch,negSwitch)
        call assert_true(.not.compFlag)


        ! todo! have to write a only excitation calculating funciton
        ! not including tmat for this case!

        print *, "calcDoubleExcitation_withWeight tests passed!"

    end subroutine test_calcDoubleExcitation_withWeight

    subroutine test_checkCompatibility
        character(*), parameter :: this_routine = "test_checkCompatibility"
        real(dp) :: posSwitch(4), negSwitch(4)
        integer(n_int) :: ilut(0:nifguga)
        logical :: flag
        type(ExcitationInformation_t) :: excitInfo

        print *, "testing: checkCompatibility:"
        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)
        excitInfo = excitationIdentifier_double(1,2,3,4)
        call calcRemainingSwitches_excitInfo_double(excitInfo, posSwitch, negSwitch)
        call checkCompatibility(ilut, excitInfo, flag, posSwitch, negSwitch)

        call assert_true(.not.flag)

        ! 3300

        excitInfo = excitationIdentifier_double(1,2,1,4)
        call checkCompatibility(ilut, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(.not.flag)

        excitInfo = excitationIdentifier_double(3,2,4,1)
        call checkCompatibility(ilut, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(flag)

        ! 0033
        call EncodeBitDet_guga([5,6,7,8],ilut)
        call init_csf_information(ilut)

        excitInfo = excitationIdentifier_double(1,2,1,4)
        call checkCompatibility(ilut, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(.not.flag)

        excitInfo = excitationIdentifier_double(1,3,1,4)
        call checkCompatibility(ilut, excitInfo, flag, posSwitch, negSwitch)
        call assert_true(flag)

        print *, "checkCompatibility tests passed!"


    end subroutine test_checkCompatibility


    subroutine test_calcRemainingSwitches_double
        character(*), parameter :: this_routine = "test_calcRemainingSwitches_double"
        real(dp) :: posSwitch(4), negSwitch(4)
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        print *, "testing: calcRemainingSwitches_double"
        ! solche scheiß tests...
        call calcRemainingSwitches_double(1,2,3,4,posSwitch,negSwitch)
        call assert_true(all(posSwitch .isclose. 0.0_dp))
        call assert_true(all(negSwitch .isclose. 0.0_dp))
        call calcRemainingSwitches_double(3,2,3,1,posSwitch,negSwitch)
        call assert_true(all(posSwitch .isclose. 0.0_dp))
        call assert_true(all(negSwitch .isclose. 0.0_dp))
        call calcRemainingSwitches_double(1,4,3,4,posSwitch,negSwitch)
        call assert_true(all(posSwitch .isclose. 0.0_dp))
        call assert_true(all(negSwitch .isclose. 0.0_dp))

        call EncodeBitDet_guga([1,4,5,6], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        ! 1212 todo

        print *, "calcRemainingSwitches_double tests passed!"

    end subroutine test_calcRemainingSwitches_double


    subroutine test_excitationIdentifier_double
        character(*), parameter :: this_routine = "test_excitationIdentifier_double"
        type(ExcitationInformation_t) :: excitInfo


        print *, "testing: excitationIdentifier_double"

        excitInfo = excitationIdentifier_double(1,2,3,4)
        call assert_true(excitInfo%i==1)
        call assert_true(excitInfo%j==2)
        call assert_true(excitInfo%gen1==1)
        call assert_true(excitInfo%gen2==1)
        call assert_true(excitInfo%typ==3)
        excitInfo = excitationIdentifier_double(1,2,2,4)
        call assert_true(excitInfo%fullStart==1)
        call assert_true(excitInfo%secondStart==2)
        call assert_true(excitInfo%firstEnd==2)
        call assert_true(excitInfo%fullEnd==4)
        call assert_true(excitInfo%currentGen==1)
        call assert_true(excitInfo%typ==5)

        excitInfo = excitationIdentifier_double(3,2,3,4)
        call assert_true(excitInfo%typ==6)
        excitInfo = excitationIdentifier_double(4,2,3,1)
        call assert_true(excitInfo%typ==8)
        excitInfo = excitationIdentifier_double(1,1,3,4)
        call assert_true(excitInfo%typ==1)
        excitInfo = excitationIdentifier_double(1,1,4,4)
        call assert_true(excitInfo%typ==-2)
        excitInfo = excitationIdentifier_double(1,1,1,1)
        call assert_true(excitInfo%typ==-2)
        excitInfo = excitationIdentifier_double(1,3,2,4)
        call assert_true(excitInfo%typ==9)
        excitInfo = excitationIdentifier_double(1,4,3,2)
        call assert_true(excitInfo%typ==11)

        print *, "excitationIdentifier_double tests passed!"

    end subroutine test_excitationIdentifier_double

    subroutine test_calcAllExcitations_single
        character(*), parameter :: this_routine = "test_calcAllExcitations_single"
        integer(n_int) :: ilut(0:nifguga)
        integer(n_int), pointer :: ex(:,:)
        integer :: nEx

        ! 3300
        call EncodeBitDet_guga([1,2,3,4], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        print *, "testing: calcAllExcitations_single"
        call calcAllExcitations_single(ilut, 4, 1, ex, nEx)
        call assert_equals(-Root2, extract_matrix_element(ex(:,1),1))
        call assert_equals(1, nEx)

        deallocate(ex)

        ! 0123
        call EncodeBitDet_guga([3,6,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call calcAllExcitations_single(ilut, 2, 4, ex, nEx)

        call assert_equals(1, nex)
        call assert_equals(-1.0_dp, extract_matrix_element(ex(:,1),1),1e-10_dp)

        ! 0123
        ! 0312

        print *, "calcAllExcitations_single tests passed!"


    end subroutine test_calcAllExcitations_single


    subroutine test_singleEnd
        character(*), parameter :: this_routine = "test_singleEnd"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitch(nBasis/2), negSwitch(nBasis/2)
        integer(n_int), pointer :: excits(:,:), tmpEx(:,:)
        integer :: num
        type(WeightObj_t) :: weights

        print *, "testing: singleEnd(ilut, exInfo, tmpEx, nEx, excits)"
        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1,4,7,8], ilut)
        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, 4)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            tmpEx, num)

        call assert_true(getDeltaB(tmpEx) == 1)
        call singleUpdate(ilut, 2, excitInfo, posSwitch, negSwitch, weights, &
            tmpEx, num)

        call assert_true(getDeltaB(tmpEx) == -1)
        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(tmpEx(:,1),1) + OverR2) < 1.0e-10_dp)

        call singleUpdate(ilut, 3, excitInfo, posSwitch, negSwitch, weights, &
            tmpEx, num)
        call assert_true(num == 1)
        call assert_true(getDeltaB(tmpEx) == -1)
        call assert_true(abs(extract_matrix_element(tmpEx(:,1),1) + OverR2) < 1.0e-10_dp)

        call singleEnd(ilut, excitInfo, tmpEx, num, excits)

        call assert_true(.not. associated(tmpEx))
        call assert_true( num == 1)
        call assert_true(abs(extract_matrix_element(excits(:,1),1) + 1.0_dp) < 1.0e-10_dp)

        print *, "singleEnd tests passed!"

    end subroutine test_singleEnd

    subroutine test_singleUpdate
        character(*), parameter :: this_routine = "test_singleUpdate"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitch(nBasis/2), negSwitch(nBasis/2)
        integer(n_int), pointer :: excits(:,:)
        integer :: num
        type(WeightObj_t) :: weights

        print *, "testing: singleUpdate(ilut,orb,exInfo,posSwitch,negSwitch,excits,num)"
        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1,4,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut, 4)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)


        call singleUpdate(ilut, 2, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(excits(:,1),1) + OverR2) < 1.0e-10_dp)

        call singleUpdate(ilut, 3, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        ! 1203
        ! 310x
        call assert_true(num == 1)
        call assert_true(abs(extract_matrix_element(excits(:,1),1) + OverR2) < 1.0e-10_dp)

        print *, "singleUpdate tests passed!"

        deallocate(excits)

    end subroutine test_singleUpdate


    subroutine test_createSingleStart
        character(*), parameter :: this_routine = "test_createSingleStart"
        integer(n_int) :: ilut(0:nifguga)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitch(nBasis/2), negSwitch(nBasis/2)
        integer(n_int), pointer :: excits(:,:)
        integer :: num
        type(WeightObj_t) :: weights

        print *, "testing: createSingleStart(ilut, excitInfo, posSwitch, negSwitch, excits, nExcits)"

        ! test a [1,2,0,3] E_1,4 raising start:
        call EncodeBitDet_guga([1,4,7,8], ilut)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        excitInfo = excitationIdentifier(1, 4)
        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)
        weights = init_singleWeight(ilut,4)
        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights,&
            excits, num)
        call assert_true(num == 1)
        call assert_true(extract_matrix_element(excits(:,1),1) .isclose. Root2)
        call assert_true(getDeltaB(excits(:,1)) == 1)

        deallocate(excits)

        call EncodeBitDet_guga([1,4,7,8], ilut)

        excitInfo = excitationIdentifier(2,4)

        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        weights = init_singleWeight(ilut,4)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        call assert_true(num == 1)
        call assert_true(getDeltaB(excits(:,1)) == -1)

        deallocate(excits)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        excitInfo = excitationIdentifier(4, 1)


        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        weights = init_singleWeight(ilut,4)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)
        call assert_true(num == 1)
        call assert_true(getDeltaB(excits(:,1)) == -1)
        call assert_true(extract_matrix_element(excits(:,1),1) .isclose. Root2)

        ! and do one double start
        deallocate(excits)

        ! 1320
        ! 1122
        ! 1212
        call EncodeBitDet_guga([1,3,4,6], ilut)
        excitInfo = excitationIdentifier(4,2)

        currentB_ilut = calcB_vector_ilut(ilut)
        currentOcc_ilut = calcOcc_vector_ilut(ilut)
        currentOcc_int = calcOcc_vector_int(ilut)
        current_stepvector = calcStepVector(ilut)
        currentB_int = calcB_vector_int(ilut)

        weights = init_singleWeight(ilut,4)

        call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitch, negSwitch)

        call createSingleStart(ilut, excitInfo, posSwitch, negSwitch, weights, &
            excits, num)

        call assert_true(num == 2)
        call assert_true(all(calcStepVector(excits(:,1)) == [1,1,2,0]))
        call assert_true(all(calcStepVector(excits(:,2)) == [1,2,2,0]))
        call assert_true(abs(extract_matrix_element(excits(:,1),1) - sqrt(3.0_dp/2.0_dp)) < 1.0e-10_dp)
        call assert_true(abs(extract_matrix_element(excits(:,2),1) - OverR2) < 1.0e-10_dp)


        print *, "createSingleStart tests passed!"

        deallocate(excits)

    end subroutine test_createSingleStart

    subroutine test_excitationIdentifier_single
        character(*), parameter :: this_routine = "test_excitationIdentifier_single"
        type(ExcitationInformation_t) :: excitInfo

        print *, "testing: excitationIdentifier_single:"
        excitInfo = excitationIdentifier(1, 4)
        call assert_true(excitInfo%i==1)
        call assert_true(excitInfo%j==4)
        call assert_true(excitInfo%gen1==1)
        call assert_true(excitInfo%fullStart==1)
        call assert_true(excitInfo%fullEnd==4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == 0)

        excitInfo = excitationIdentifier(1, 4)
        call assert_true(excitInfo%i==1)
        call assert_true(excitInfo%j==4)
        call assert_true(excitInfo%gen1==1)
        call assert_true(excitInfo%fullStart==1)
        call assert_true(excitInfo%fullEnd==4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == 0)

        excitInfo = excitationIdentifier(2, 4)
        call assert_true(excitInfo%i==2)
        call assert_true(excitInfo%j==4)
        call assert_true(excitInfo%gen1==1)
        call assert_true(excitInfo%fullStart==2)
        call assert_true(excitInfo%fullEnd==4)
        call assert_true(excitInfo%currentGen == 1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == 0)

        excitInfo = excitationIdentifier(3, 2)
        call assert_true(excitInfo%i==3)
        call assert_true(excitInfo%j==2)
        call assert_true(excitInfo%gen1==-1)
        call assert_true(excitInfo%fullStart==2)
        call assert_true(excitInfo%fullEnd==3)
        call assert_true(excitInfo%currentGen == -1)
        call assert_true(excitInfo%excitLvl == 2)
        call assert_true(excitInfo%typ == 0)

        print *, "excitationIdentifier_single tests passed!"

    end subroutine test_excitationIdentifier_single


    subroutine test_getDoubleMatrixElement
        character(*), parameter :: this_routine = "test_getDoubleMatrixElement"
        real(dp) :: x0, x1

        print *, "testing: getDoubleMatrixElement:"
        call getDoubleMatrixElement(1,1,0,1,1,1.0_dp,1.0_dp,x0,x1)
        ! todo more extensive tests here.. but not for now..
        print *, "x0 = ", x0, " x1 = ", x1
        print *, "getDoubleMatrixElement tests passed!"

    end subroutine test_getDoubleMatrixElement

    subroutine test_getMixedFullStop
        character(*), parameter :: this_routine = "test_getMixedFullStop"
        real(dp) :: x0, x1

        print *, "testing: getMixedFullStop:"
        call getMixedFullStop(1,1,0,1.0_dp,x0,x1)
        call assert_true(abs(x0 - OverR2) < 1.0e-10_dp)
        call assert_true(abs(x1) < EPS)
        call getMixedFullStop(0,0,0,1.0_dp,x0,x1)
        call assert_true(abs(x0) < EPS)
        call assert_true(abs(x1) < EPS)
        call getMixedFullStop(3,3,0,2.0_dp,x0,x1)
        call assert_true(abs(x0 - Root2) < 1.0e-10_dp)
        call assert_true(abs(x1) < EPS)
        call getMixedFullStop(2,2,0,3.0_dp,x0,x1)
        call assert_true(abs(x0 - OverR2) < 1.0e-10_dp)
        call assert_true(abs(x1 - sqrt(6.0_dp/8.0_dp)) < 1.0e-10_dp)
        call getMixedFullStop(1,2,2,2.0_dp,x0,x1)
        call assert_true(abs(x0) < EPS)
        call assert_true(abs(x1 - 1.0_dp) < 1.0e-10_dp)
        call getMixedFullStop(2,1,-2,1.0_dp,x0,x1)
        call assert_true(abs(x0) < EPS)
        call assert_true(abs(x1 - 1.0_dp) < 1.0e-10_dp)
        call getMixedFullStop(1,1,0,2.0_dp,x0,x1)
        call assert_true(abs(x1 + 1.0_dp / sqrt(6.0_dp)) < 1.0e-10_dp)

        print *, "x0 = ", x0, " x1 = ", x1
        print *, "getMixedFullStop tests passed!"

    end subroutine test_getMixedFullStop


    subroutine test_getSingleMatrixElement
        character(*), parameter :: this_routine = "test_getSingleMatrixElement"

        print *, "testing: getSingleMatrixElement(d1,d2,dB,gen,b):"
        call assert_true(abs(getSingleMatrixElement(0,0,-1,1,1.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0,0,1,1,1.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0,0,1,-1,1.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0,0,-1,-1,1.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(3,3,-1,1,1.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(3,3,1,1,1.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(3,3,1,-1,1.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(3,3,-1,-1,1.0_dp) + 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(0,1,-1,1,2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0,1,1,-1,2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(0,2,1,1,2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(0,2,-1,-1,2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1,0,-1,1,2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1,0,-1,-1,2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2,0,1,1,2.0_dp) - 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2,0,+1,-1,2.0_dp) - 1.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1,1,-1,1,3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1,1,1,1,3.0_dp) - sqrt(8.0_dp)/3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1,1,1,-1,3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1,1,-1,-1,2.0_dp) - sqrt(8.0_dp)/3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2,2,1,1,3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2,2,-1,1,1.0_dp) - sqrt(8.0_dp)/3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2,2,-1,-1,3.0_dp) + 1.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2,2,1,-1,2.0_dp) - sqrt(8.0_dp)/3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(1,2,1,1,2.0_dp) + 1.0_dp/4.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(1,2,1,-1,2.0_dp) - 1.0_dp/3.0_dp) < tol)

        call assert_true(abs(getSingleMatrixElement(2,1,-1,1,2.0_dp) - 1.0_dp/2.0_dp) < tol)
        call assert_true(abs(getSingleMatrixElement(2,1,-1,-1,2.0_dp) + 1.0_dp/3.0_dp) < tol)

        print *, "getSingleMatrixElement tests passed!"

    end subroutine test_getSingleMatrixElement



    subroutine test_matrix_element_ops
        character(*), parameter :: this_routine ="test_encode_matrix_element"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: encode_matrix_element(ilut, ele, type):"
        print *, "         extract_matrix_element(ilut, type):"
        print *, "         update_matrix_element(ilut, ele, type):"

        print *, "nifguga: ", nifguga
        print *, "niftot: ", niftot

        call encode_matrix_element(ilut,1.0_dp,1)
        call assert_true(extract_matrix_element(ilut, 1)  .isclose. 1.0_dp)
        call encode_matrix_element(ilut,-1.0_dp,2)
        call assert_true(extract_matrix_element(ilut, 2)  .isclose. -1.0_dp)
        call update_matrix_element(ilut, 2.0_dp, 1)
        call assert_true(extract_matrix_element(ilut, 1)  .isclose. 2.0_dp)
        call update_matrix_element(ilut, -2.0_dp, 2)
        call assert_true(extract_matrix_element(ilut, 2)  .isclose. 2.0_dp)

        print *, "encode_matrix_element tests passed!"
    end subroutine test_matrix_element_ops

    subroutine test_calcOcc_vector_ilut
        character(*), parameter :: this_routine = "test_calcOcc_vector_ilut"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: calcOcc_vector_ilut(ilut)"
        call assert_true(all([2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp]  .isclose. calcOcc_vector_ilut(ilut)))

        call EncodeBitDet_guga([1,4,5,6], ilut)
        call assert_true(all([1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp]  .isclose. calcOcc_vector_ilut(ilut)))

        call EncodeBitDet_guga([1,2,3,8], ilut)
        call assert_true( all([2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp] .isclose. calcOcc_vector_ilut(ilut)))
        print *, "calcOcc_vector_ilut tests passed!"

    end subroutine test_calcOcc_vector_ilut

    subroutine test_calcStepVector
        character(*), parameter :: this_routine = "test_calcStepVector"
        integer(n_int) :: ilut(0:nifguga)

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: calcStepVector(ilut)"
        call assert_true(all([3,3,0,0] == calcStepVector(ilut)))
        call EncodeBitDet_guga([1,4,5,6], ilut)
        call assert_true(all([1,2,3,0] == calcStepVector(ilut)))
        print *, "calcStepVector tests passed!"


    end subroutine test_calcStepVector

    subroutine test_isDouble
        character(*), parameter :: this_routine = "test_isDouble"

        print *, "testing: isDouble(nI, sOrb)"
        call assert_true(isDouble([1,2,3,4],1))
        call assert_true(isDouble([1,2,3,6],2))
        call assert_true(.not.isDouble([1,2,3,6],3))
        print *, "isDouble tests passed!"

    end subroutine test_isDouble

    subroutine test_isProperCSF_ilut
        integer(n_int) :: ilut(0:nifguga)
        character(*), parameter :: this_routine = "test_isProperCSF_ilut"

        call EncodeBitDet_guga([1,2,3,4], ilut)

        print *, "testing: isProperCSF_ilut(ilut)"
        call assert_true(isProperCSF_ilut(ilut))
        call EncodeBitDet_guga([2,3,4,5],ilut)
        call assert_true(.not.isProperCSF_ilut(ilut))
        print *, "isProperCSF_ilut tests passed!"

    end subroutine test_isProperCSF_ilut

    subroutine test_set_get_DeltaB
        integer(n_int) :: ilut(0:nifguga)
        integer :: deltaB
        character(*), parameter :: this_routine = "test_set_get_DeltaB"

        call EncodeBitDet_guga([1,2,3,4], ilut)

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
        character(*), parameter :: this_routine = "test_count_open_orbs_ij"

        det = [1,2,3,4]

        call EncodeBitDet_guga(det, ilut)

        print *, "testing: count_open_orbs_ij(ilut, i, j)"
        print *, "open orbs in ([1,2,3,4],1,4): ", count_open_orbs_ij(1,4, ilut(0:0))
        call assert_true(count_open_orbs_ij(1,4, ilut) == 0)

        det = [1, 3, 5, 6]

        call EncodeBitDet_guga(det, ilut)
        print *, "open orbs in ([1, 3, 5, 6], 1, 4): ", count_open_orbs_ij(1, 4, ilut)
        call assert_true(count_open_orbs_ij(1,4, ilut) == 2)

        print *, "count_open_orbs_ij tests passed!"

    end subroutine test_count_open_orbs_ij


    subroutine test_getExcitationRangeMask

        print *, "testing: getExcitationRangeMask(i,j)"
        print *, "i = 1, j = 4 :", getExcitationRangeMask(1, 4)
        print *, "i = 2, j = 3 :", getExcitationRangeMask(2, 3)
        print *, "getExcitationRangeMask tests passed!"

    end subroutine test_getExcitationRangeMask



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
        det = [1,2,3,6]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3300
        det = [1,2,3,4]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0330
        det = [3,4,5,6]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0303
        det = [3,4,7,8]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0033
        det = [5,6,7,8]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1023
        det = [1,6,7,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3102
        det = [1,2,3,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3030
        det = [1,2,5,6]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3003
        det = [1,2,7,8]
        call assert_equals(h_cast(8.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 3012
        det = [1,2,5,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0312
        det = [3,4,5,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1230
        det = [1,4,5,6]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1203
        det = [1,4,7,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1320
        det = [1,3,4,6]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1302
        det = [1,3,4,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1032
        det = [1,5,6,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0132
        det = [3,5,6,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 0123
        det = [3,6,7,8]
        call assert_equals(h_cast(9.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1122
        det = [1,3,6,8]
        call assert_equals(h_cast(10.0_dp), calcDiagMatEleGUGA_nI(det))

        ! 1212
        det = [1,4,5,8]
        call assert_equals(h_cast(10.0_dp), calcDiagMatEleGUGA_nI(det))

        print *, this_routine, " tests passed!"

    end subroutine check_calcDiagMatEleGUGA_nI


    subroutine check_calcDiagExchange_nI
        integer :: det(4), iOrb, jOrb

        det = [1,2,3,6]
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
        real(dp) :: checkB_nI(4), checkB_ilut(nBasis/2)
        character(*), parameter :: testFun = "calcB_vector", &
            this_routine = "test_calcbvector"

        print *, " Testing ", testFun
        det = [1,2,3,4]
        checkB_nI = 0.0_dp
        call EncodeBitDet_guga(det,ilut)
        checkB_ilut = 0.0_dp

        call assert_true(all(calcB_vector_nI(det) .isclose. checkB_nI))
        call assert_true(all(calcB_vector_ilut(ilut) .isclose. 0.0_dp))

        det = [1,2,3,6]
        call EncodeBitDet_guga(det, ilut)

        checkB_nI(3) = 1.0_dp
        checkB_ilut(2) = 1.0_dp

        call assert_true(all(calcB_vector_nI(det) .isclose. [0.0_dp,0.0_dp,1.0_dp,0.0_dp]))
        call assert_true(all(calcB_vector_ilut(ilut) .isclose. [0.0_dp,1.0_dp,0.0_dp,0.0_dp]))
        print *, testFun, " tests passed!"


    end subroutine test_calcbvector


    subroutine test_calcRemainingSwitches()
        real(dp) :: neg(nBasis/2), pos(nBasis/2)
        integer :: det(4)
        integer(n_int) :: ilut(0:nifguga)
        integer :: b
        character(*), parameter :: this_routine = "test_calcRemainingSwitches"

        det = [1,2,3,6] ! 3 1 2 0

        call EncodeBitDet_guga(det,ilut)
        ! now need excitInfo for

        current_stepvector = calcStepVector(ilut)

        print *, "***"
        print *, "testing calcRemainingSwitches:"
        call calcRemainingSwitches_single(1,4,pos,neg)
        print *, "positive switches: ", pos
        call assert_true(all(pos .isclose. [1.0_dp,1.0_dp,0.0_dp,0.0_dp]))
        print *, "negative switches: ", neg
        call assert_true(all(neg .isclose. [1.0_dp,0.0_dp,0.0_dp,0.0_dp]))

        call EncodeBitDet_guga([1,4,5,8],ilut) ! 1 2 1 2

        current_stepvector = calcStepVector(ilut)

        call calcRemainingSwitches_single(2,4,pos,neg)
        print *, "positive switches: ", pos
        call assert_true(all(pos .isclose. [0.0_dp,0.0_dp,0.0_dp,0.0_dp]))
        print *, "negative switches: ", neg
        call assert_true(all(neg .isclose. [0.0_dp,1.0_dp,0.0_dp,0.0_dp]))
        print *, "calcRemainingSwitches tests successfull"

    end subroutine test_calcRemainingSwitches

    subroutine test_calcOverlapRange()
        integer, allocatable :: overlap(:), nonOverlap(:)
        character(*), parameter :: testFun = "calcOverlapRange"
        integer :: i,j,k,l

        print *, ""
        print *, " Testing ", testFun

        ! anyway not used anymore..
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




end program test_guga
