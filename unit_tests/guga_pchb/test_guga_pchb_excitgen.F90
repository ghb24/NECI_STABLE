#include "macros.h"
program test_guga_pchb_excitgen

    use constants
    use SystemData
    use CalcData
    use FciMCData
    use dSFMT_interface
    use procedure_pointers
    use read_fci
    use UMatCache
    use Parallel_neci
    use Calc
    use System
    use shared_memory_mpi
    use DetCalc
    use unit_test_helper_excitgen
    use LoggingData
    use bit_reps
    use bit_rep_data
    use guga_init
    use guga_bitRepOps
    use guga_data
    use guga_types
    use guga_pchb_excitgen
    use util_mod

    use fruit
    use fruit_extensions

    ! use template

    implicit none

    integer :: failed_count

    call init_fruit()
    call guga_pchb_test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

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

        tRDMonfly = .true.
        tFillingStochRDMOnFly = .true.
        call init_bit_rep()
        t_full_guga_tests = .true.

        t_pchb_weighted_singles = .false.
        t_pchb_excitgen = .true.
!         tGen_sym_guga_mol = .true.
!         tgen_guga_weighted = .true.
        tdeferred_umat2d = .true.
        tumat2d = .false.

        ! set this to false before the init to setup all the ilut variables
        tExplicitAllRDM = .false.

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

    subroutine guga_pchb_test_driver()

        call init_guga_testsuite()

        call my_run_test_case(pick_uniform_spatial_hole_test, &
            "pick_uniform_spatial_hole_test", "pick_uniform_spatial_hole")
        call my_run_test_case(pick_orbitals_pure_uniform_singles_test, &
            "pick_orbitals_pure_uniform_singles_test", &
            "pick_orbitals_pure_uniform_singles")
        call my_run_test_case(calc_orb_pgen_uniform_singles_test, &
            "calc_orb_pgen_uniform_singles_test", &
            "calc_orb_pgen_uniform_singles_ex() and &
            & calc_orb_pgen_uniform_singles_excitInfo")

        call my_run_test_case(guga_pchb_sampler_test, &
            "guga_pchb_sampler_test", "guga_pchb_sampler type")

        call my_run_test_case(setup_pchb_sampler_conditional_test, &
            "setup_pchb_sampler_conditional_test", &
            "setup_pchb_sampler_conditional")



    end subroutine guga_pchb_test_driver

    subroutine guga_pchb_sampler_test

    end subroutine guga_pchb_sampler_test

    subroutine setup_pchb_sampler_conditional_test

        integer :: i, ijMax, j

        call setup_pchb_sampler_conditional()

        ! the first index: (1,1) - > 1 and (1,1) -> 1 should be 0 since
        ! only weight
        print *, "fuse"
        do i = 1, 4
            do j = 1, 4
                print *, i, j, "fuse: ", fuseIndex(i,j)
            end do
        end do
        print *, "--"
        associate(sampler => guga_pchb_sampler, &
                  a_sampler => guga_pchb_sampler%alias_sampler)

            call assert_equals(0_int64, sampler%get_info(1,1))
            call assert_equals(0.0_dp, a_sampler%get_prob(1,1))

        end associate

        ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)

        do i = 1, ijMax
            do j = 1, ijMax
                print *, guga_pchb_sampler%get_info(i,j), &
                    guga_pchb_sampler%alias_sampler%get_prob(i,j)
            end do
        end do

    end subroutine setup_pchb_sampler_conditional_test

    subroutine calc_orb_pgen_uniform_singles_test

        integer :: ex(2,2), nI(4)
        integer(n_int) :: ilut(0:GugaBits%len_tot)
        real(dp) :: pgen
        type(ExcitationInformation_t) :: excitInfo

        nI = [1,2,3,4]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        ex(:,2) = 0
        ex(1,1) = 1
        ex(2,1) = 2

        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(0.0_dp, pgen)

        excitInfo%j = 1
        excitInfo%i = 1
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(0.0_dp, pgen)

        ex(2,1) = 3
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(0.0_dp, pgen)

        ex(2,1) = 4
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(0.0_dp, pgen)

        excitInfo%i = 2
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(0.0_dp, pgen)

        ex(2,1) = 5
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(1.0_dp/4.0_dp, pgen)

        excitInfo%i = 3
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(1.0_dp/4.0_dp, pgen)

        ex(2,1) = 8
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(1.0_dp/4.0_dp, pgen)

        excitInfo%i = 4
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(1.0_dp/4.0_dp, pgen)

        nI = [1,2,3,8]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        ex(2,1) = 8
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(1.0_dp/9.0_dp, pgen)

        excitInfo%i = 4
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(1.0_dp/9.0_dp, pgen)

        ex(1,1) = 3
        ex(2,1) = 8
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(1.0_dp/6.0_dp, pgen)

        excitInfo%j = 2
        excitInfo%i = 4
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(1.0_dp/6.0_dp, pgen)

        nI = [1,4,5,8]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        ex(1,1) = 8
        ex(2,1) = 3
        pgen = calc_orb_pgen_uniform_singles(ex)
        call assert_equals(1.0_dp/12.0_dp, pgen)

        excitInfo%j = 4
        excitInfo%i = 2
        pgen = calc_orb_pgen_uniform_singles(excitInfo)
        call assert_equals(1.0_dp/12.0_dp, pgen)

    end subroutine calc_orb_pgen_uniform_singles_test

    subroutine pick_orbitals_pure_uniform_singles_test

        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer :: nI(4)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen

        nI = [1,2,3,4]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        call pick_orbitals_pure_uniform_singles(ilut, nI, excitInfo, pgen)

        call assert_equals(excitInfo%typ, excit_type%single)
        call assert_equals(gen_type%L, excitInfo%gen1)
        call assert_true(excitInfo%j == 1 .or. excitInfo%j == 2)
        call assert_true(excitInfo%i == 3 .or. excitInfo%i == 4)
        call assert_equals(0.25_dp, pgen)


        nI = [1,4,5,8]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        call pick_orbitals_pure_uniform_singles(ilut, nI, excitInfo, pgen)

        call assert_equals(excitInfo%typ, excit_type%single)
        call assert_true(any(excitInfo%j == [1,2,3,4]))
        call assert_true(any(excitInfo%i == [1,2,3,4]))
        call assert_true(excitInfo%i /= excitInfo%j)
        call assert_equals(1.0_dp/12.0_dp, pgen)

    end subroutine pick_orbitals_pure_uniform_singles_test

    subroutine pick_uniform_spatial_hole_test

        integer :: nI(4), elec, s_orb
        real(dp) :: pgen
        integer(n_int) :: ilut(0:GugaBits%len_tot)


        nI = [1,2,3,4]

        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        elec = 1
        call pick_uniform_spatial_hole(elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4)
        call assert_equals(0.5_dp, pgen)

        elec = 2
        call pick_uniform_spatial_hole(elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4)
        call assert_equals(0.5_dp, pgen)

        nI = [1,2,3,8]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        elec = 1
        call pick_uniform_spatial_hole(elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4 .or. s_orb == 2)
        call assert_equals(1.0_dp/3.0_dp, pgen)

        elec = 2
        call pick_uniform_spatial_hole(elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4)
        call assert_equals(1.0_dp/2.0_dp, pgen)

        nI = [1,4,5,8]
        call EncodeBitDet_guga(nI, ilut)
        call init_csf_information(ilut)

        elec = 1
        call pick_uniform_spatial_hole(elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4 .or. s_orb == 2)
        call assert_equals(1.0_dp/3.0_dp, pgen)







    end subroutine pick_uniform_spatial_hole_test




end program test_guga_pchb_excitgen

