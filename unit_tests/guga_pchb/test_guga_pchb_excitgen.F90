#include "macros.h"

module test_guga_pchb_mod
    use constants, only: dp, n_int
    use SystemData, only: nEl
    use FciMCData, only: pSingles
    use procedure_pointers, only: get_umat_el
    use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
    use IntegralsData, only: umat_win, umat
    use Integrals_neci, only: get_umat_el_normal
    USE OneEInts, only: TMAT2D
    use DetCalc, only: DetCalcInit
    use Determinants, only: DetInit, DetPreFreezeInit
    use unit_test_helper_excitgen, only: init_guga_testsuite, finalize_guga_testuite, &
            FciDumpWriter_t
    use orb_idx_mod, only: sum, calc_spin_raw
    use LoggingData, only: tExplicitAllRDM, tRDMonFly
    use bit_reps, only: init_bit_rep
    use bit_rep_data, only: GugaBits
    use guga_init, only: init_guga
    use guga_bitRepOps, only: CSF_Info_t, EncodeBitDet_guga
    use guga_data, only: ExcitationInformation_t, excit_type, gen_type
    use guga_pchb_class, only: PCHB_t, options, FCI_PCHB_options_t, &
        FCI_PCHB_options_vals_t, options_vals


    use guga_gas_pchb_singles_main, only: GAS_PCHB_SinglesOptions_t
    use guga_gas_pchb_doubles_main, only: PCHB_DoublesOptions_t

    use fruit, only: assert_true, assert_equals
    use fruit_extensions, only: my_run_test_case

    better_implicit_none
    private
    public :: guga_pchb_test_driver

    type(FCI_PCHB_options_t), parameter :: chosen_options = &
            FCI_PCHB_options_t(&
                GAS_PCHB_SinglesOptions_t(&
                    options_vals%singles%algorithm%UNIF_UNIF, &
                    options_vals%singles%weighting%UNIFORM&
                ), &
                PCHB_DoublesOptions_t(&
                    options_vals%doubles%particle_selection%UNIF_UNIF, &
                    options_vals%doubles%hole_selection%FAST_FAST, &
                    options_vals%doubles%weighting%INITIAL_DOBRAUTZ &
                ), &
                use_lookup=.false. &
            )

    type(PCHB_t) :: exc_generator
contains

    subroutine guga_pchb_test_driver()

        call init_guga_testsuite()
        call finalize_guga_testuite()

        call init_guga_testsuite()
        call my_run_test_case(pick_orbitals_pure_uniform_singles_test, &
                              "pick_orbitals_pure_uniform_singles_test", &
                              "pick_orbitals_pure_uniform_singles")
        call my_run_test_case(calc_orb_pgen_uniform_singles_test, &
            "calc_orb_pgen_uniform_singles_test", &
            "calc_orb_pgen_uniform_singles_ex() and &
            & calc_orb_pgen_uniform_singles_excitInfo")

        call my_run_test_case(guga_pchb_sampler_test, &
                              "guga_pchb_sampler_test", "guga_pchb_sampler type")

        call finalize_guga_testuite()
    end subroutine guga_pchb_test_driver

    subroutine guga_pchb_sampler_test
    end subroutine guga_pchb_sampler_test

    subroutine calc_orb_pgen_uniform_singles_test

        integer :: nI(4)
        integer(n_int) :: ilutI(0:GugaBits%len_tot)
        real(dp) :: pgen
        type(ExcitationInformation_t) :: excitInfo
        type(CSF_Info_t) :: csf_info

        pSingles = 1._dp

        nI = [1, 2, 3, 4]
        call exc_generator%init(chosen_options)
        call EncodeBitDet_guga(nI, ilutI)
        csf_info = CSF_Info_t(ilutI)

        excitInfo%j = 1
        excitInfo%i = 1
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(0.0_dp, pgen)

        excitInfo%i = 2
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(0.0_dp, pgen)

        excitInfo%i = 3
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(1.0_dp / 4.0_dp, pgen)

        excitInfo%i = 4
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(1.0_dp / 4.0_dp, pgen)

        nI = [1, 2, 3, 8]
        call EncodeBitDet_guga(nI, ilutI)
        csf_info = CSF_Info_t(ilutI)

        excitInfo%i = 4
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(1.0_dp / 9.0_dp, pgen)

        excitInfo%j = 2
        excitInfo%i = 4
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(1.0_dp / 6.0_dp, pgen)

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilutI)
        csf_info = CSF_Info_t(ilutI)

        excitInfo%j = 4
        excitInfo%i = 2
        pgen = exc_generator%calc_pgen(nI, ilutI, csf_info, excitInfo)
        call assert_equals(1.0_dp / 12.0_dp, pgen)

        call exc_generator%finalize()
    end subroutine calc_orb_pgen_uniform_singles_test

    subroutine pick_orbitals_pure_uniform_singles_test
        integer(n_int) :: ilut(0:GugaBits%len_tot)
        integer :: nI(4)
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        type(CSF_Info_t) :: csf_info

        call exc_generator%init(chosen_options)

        nI = [1, 2, 3, 4]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        call pick_orbitals_pure_uniform_singles(ilut, nI, csf_info, excitInfo, pgen)

        call assert_true(excitInfo%typ == excit_type%single)
        call assert_equals(gen_type%L, excitInfo%gen1)
        call assert_true(excitInfo%j == 1 .or. excitInfo%j == 2)
        call assert_true(excitInfo%i == 3 .or. excitInfo%i == 4)
        call assert_equals(0.25_dp, pgen)

        nI = [1, 4, 5, 8]
        call EncodeBitDet_guga(nI, ilut)
        csf_info = CSF_Info_t(ilut)

        call pick_orbitals_pure_uniform_singles(ilut, nI, csf_info, excitInfo, pgen)

        call assert_true(excitInfo%typ == excit_type%single)
        call assert_true(any(excitInfo%j == [1, 2, 3, 4]))
        call assert_true(any(excitInfo%i == [1, 2, 3, 4]))
        call assert_true(excitInfo%i /= excitInfo%j)
        call assert_equals(1.0_dp / 12.0_dp, pgen)

        call exc_generator%finalize()

    contains

        subroutine pick_orbitals_pure_uniform_singles(ilut, nI, csf_i, excitInfo, pgen)
            integer(n_int), dimension(0:gugabits%len_tot), intent(in) :: ilut
            integer, dimension(nel), intent(in) :: ni
            type(csf_info_t), intent(in) :: csf_i
            type(excitationinformation_t), intent(out) :: excitinfo
            real(dp), intent(out) :: pgen
            call exc_generator%pickOrbitals_single(nI, ilut, csf_i, excitInfo, pgen)
        end subroutine

    end subroutine pick_orbitals_pure_uniform_singles_test


end module


program test_guga_pchb_excitgen

    use util_mod, only: stop_all
    use test_guga_pchb_mod, only: guga_pchb_test_driver
    use Parallel_neci, only: MPIInit, MPIEnd

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, get_failed_count

    implicit none

    integer :: failed_count

    call MPIInit(.false.)
    call init_fruit()
    call guga_pchb_test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) call stop_all('test_guga_pchb_excitgen', 'Tests failed')

    call MPIEnd(.false.)

end program test_guga_pchb_excitgen
