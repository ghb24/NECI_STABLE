#include "macros.h"
program test_guga_pchb_excitgen

    use constants, only: dp, int64, n_int
    use SystemData, only: tstorespinorbs, tReadInt, UMatEps, tReadFreeFormat, &
        tGUGA, t_guga_pchb_weighted_singles, t_guga_pchb, nSpatOrbs, &
        nel, nbasis
    use FciMCData, only: tFillingStochRDMOnFly
    use dSFMT_interface, only: dSFMT_init
    use procedure_pointers, only: get_umat_el
    use read_fci, only: initfromfcid, GetUMatSize, readfciint, fcidump_name, &
        tumat2d, TMat2d
    use UMatCache, only: tTransGTID, tdeferred_umat2d
    use Parallel_neci, only: MPIInit, MPIEnd
    use Calc, only: CalcInit, SetCalcDefaults
    use System, only: SetSysDefaults, SysInit
    use shared_memory_mpi, only: shared_allocate_mpi
    use DetCalc, only: DetCalcInit
    use unit_test_helper_excitgen, only: generate_uniform_integrals
    use IntegralsData, only: Umat, umat_win
    use Integrals_Neci, only: get_umat_el_normal
    use LoggingData, only: tExplicitAllRDM, tRDMonfly
    use bit_reps, only:  init_bit_rep
    use bit_rep_data, only: GugaBits
    use guga_init, only: init_guga
    use guga_bitRepOps, only: CSF_Info_t, encodebitdet_guga
    use guga_data, only: ExcitationInformation_t, excit_type, gen_type
    use guga_pchb_class, only: GugaAliasSampler_t, &
        calc_orb_pgen_uniform_singles, pick_uniform_spatial_hole, &
        pick_orbitals_pure_uniform_singles
    use Determinants, only: DetInit, DetPreFreezeInit

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case, assert_equals, assert_true
    use fruit_extensions, only: my_run_test_case

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
        integer :: nBasisMax(5, 3), lms, stot
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

        t_guga_pchb_weighted_singles = .false.
        t_guga_pchb = .true.
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

        call initfromfcid(nel, nbasismax, nBasis, lms, .false.)

        call GetUMatSize(nBasis, umatsize)

        allocate (TMat2d(nBasis, nBasis))

        call shared_allocate_mpi(umat_win, umat, (/umatsize/))

        call readfciint(UMat, umat_win, nBasis, ecore)
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

    end subroutine guga_pchb_test_driver

    subroutine guga_pchb_sampler_test

    end subroutine guga_pchb_sampler_test

    subroutine calc_orb_pgen_uniform_singles_test

        integer :: ex(2, 2), nI(4)
        integer(n_int) :: ilut(0:GugaBits%len_tot)
        real(dp) :: pgen
        type(ExcitationInformation_t) :: excitInfo
        type(CSF_Info_t) :: csf_info

        nI = [1, 2, 3, 4]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        ex(:, 2) = 0
        ex(1, 1) = 1
        ex(2, 1) = 2

        excitInfo%j = 1
        excitInfo%i = 1
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(0.0_dp, pgen)

        excitInfo%i = 2
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(0.0_dp, pgen)

        excitInfo%i = 3
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(1.0_dp / 4.0_dp, pgen)

        excitInfo%i = 4
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(1.0_dp / 4.0_dp, pgen)

        nI = [1, 2, 3, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        excitInfo%i = 4
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(1.0_dp / 9.0_dp, pgen)

        excitInfo%j = 2
        excitInfo%i = 4
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(1.0_dp / 6.0_dp, pgen)

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        excitInfo%j = 4
        excitInfo%i = 2
        pgen = calc_orb_pgen_uniform_singles(csf_info, excitInfo)
        call assert_equals(1.0_dp / 12.0_dp, pgen)

    end subroutine calc_orb_pgen_uniform_singles_test

    subroutine pick_orbitals_pure_uniform_singles_test
        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer :: nI(4)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        type(CSF_Info_t) :: csf_info

        nI = [1, 2, 3, 4]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        call pick_orbitals_pure_uniform_singles(ilut, nI, csf_info, excitInfo, pgen)

        call assert_equals(excitInfo%typ, excit_type%single)
        call assert_equals(gen_type%L, excitInfo%gen1)
        call assert_true(excitInfo%j == 1 .or. excitInfo%j == 2)
        call assert_true(excitInfo%i == 3 .or. excitInfo%i == 4)
        call assert_equals(0.25_dp, pgen)

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        call pick_orbitals_pure_uniform_singles(ilut, nI, csf_info, excitInfo, pgen)

        call assert_equals(excitInfo%typ, excit_type%single)
        call assert_true(any(excitInfo%j == [1, 2, 3, 4]))
        call assert_true(any(excitInfo%i == [1, 2, 3, 4]))
        call assert_true(excitInfo%i /= excitInfo%j)
        call assert_equals(1.0_dp / 12.0_dp, pgen)

    end subroutine pick_orbitals_pure_uniform_singles_test

    subroutine pick_uniform_spatial_hole_test

        integer :: nI(4), elec, s_orb
        real(dp) :: pgen
        integer(n_int) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t) :: csf_info

        nI = [1, 2, 3, 4]

        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        elec = 1
        call pick_uniform_spatial_hole(csf_info, elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4)
        call assert_equals(0.5_dp, pgen)

        elec = 2
        call pick_uniform_spatial_hole(csf_info, elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4)
        call assert_equals(0.5_dp, pgen)

        nI = [1, 2, 3, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        elec = 1
        call pick_uniform_spatial_hole(csf_info, elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4 .or. s_orb == 2)
        call assert_equals(1.0_dp / 3.0_dp, pgen)

        elec = 2
        call pick_uniform_spatial_hole(csf_info, elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4)
        call assert_equals(1.0_dp / 2.0_dp, pgen)

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        elec = 1
        call pick_uniform_spatial_hole(csf_info, elec, s_orb, pgen)
        call assert_true(pgen > 0.0_dp)
        call assert_true(s_orb == 3 .or. s_orb == 4 .or. s_orb == 2)
        call assert_equals(1.0_dp / 3.0_dp, pgen)
    end subroutine pick_uniform_spatial_hole_test

end program test_guga_pchb_excitgen
