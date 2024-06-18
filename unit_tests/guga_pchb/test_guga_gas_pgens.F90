#include "macros.h"

module test_guga_gas_pgens_mod
    use constants, only: dp, n_int, int64, stdout
    use fortran_strings, only: Token_t, join, replace, to_lower
    use SystemData, only: nEl, nSpatOrbs
    use FciMCData, only: pSingles
    use procedure_pointers, only: get_umat_el
    use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
    use IntegralsData, only: umat_win, umat
    use Integrals_neci, only: get_umat_el_normal
    USE OneEInts, only: TMAT2D
    use DetCalc, only: DetCalcInit
    use Determinants, only: DetInit, DetPreFreezeInit
    use unit_test_helper_excitgen, only: init_guga_testsuite, finalize_guga_testuite, init_excitgen_test, &
        finalize_excitgen_test

    use unit_test_helper_fcidumps, only: FciDumpWriter_t, RandomFciDumpWriter_t, BenzeneDimerFciDumpWriter_t
    use gasci, only: CAS_spec
    use gasci_supergroup_index, only: get_CAS_indexer, GAS_SuperGroupIndexer_t
    use orb_idx_mod, only: sum, calc_spin_raw
    use LoggingData, only: tExplicitAllRDM, tRDMonFly
    use bit_reps, only: init_bit_rep
    use bit_rep_data, only: GugaBits
    use guga_init, only: init_guga
    use guga_bitRepOps, only: CSF_Info_t, EncodeBitDet_guga, StepVector_t, stepvector_vals
    use guga_data, only: ExcitationInformation_t, excit_type, gen_type, DistinctDouble_t, distinct_doubles
    use guga_pchb_class, only: PCHB_t, options, FCI_PCHB_options_t, &
        FCI_PCHB_options_vals_t, options_vals
    use guga_base_class, only: test_pgen, Entry_t, analyze_entries, AnalyzedEntries_t, output
    use util_mod, only: near_zero, operator(.implies.)


    use guga_gas_pchb_doubles_main, only: PCHB_DoublesOptions_t
    use guga_gas_pchb_doubles_select_holes, only: GAS_PCHB_FastFast_t

    use fruit, only: assert_true, assert_equals, run_test_case

    better_implicit_none
    private
    public :: guga_pchb_test_driver

    type(StepVector_t), parameter :: s0 = stepvector_vals%s0, u = stepvector_vals%u, &
        d = stepvector_vals%d, s2 = stepvector_vals%s2

contains

    subroutine guga_pchb_test_driver()
        call run_test_case(test_guga_pchb_sampler, "test_guga_pchb_sampler")
    end subroutine guga_pchb_test_driver

    subroutine test_guga_pchb_sampler
        type(BenzeneDimerFciDumpWriter_t) :: fcidump_writer
        type(CSF_Info_t) :: csfs(2)
        type(PCHB_DoublesOptions_t) :: options(2)
        integer, parameter :: n_el = 12, ref_csf(n_el) = [1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18]
        integer(int64) :: n_iter = 5_int64 * 10_int64**4_int64
        integer :: i, j

        fcidump_writer = BenzeneDimerFciDumpWriter_t(filepath='FCIDUMP')

        call init_excitgen_test(ref_csf, fcidump_writer, .true.)

        csfs = [CSF_Info_t([u, u, d, d, u, u, d, d, u, u, d, d]), &
                CSF_Info_t([s2, s2, s2, s0, s0, s0, s2, s2, s2, s0, s0, s0])]

        options =  [ &
                    PCHB_DoublesOptions_t(&
                        options_vals%doubles%particle_selection%UNIF_UNIF, &
                        options_vals%doubles%hole_selection%FAST_FAST, &
                        options_vals%doubles%weighting%INITIAL_DOBRAUTZ &
                    ), &
                    PCHB_DoublesOptions_t(&
                        options_vals%doubles%particle_selection%FULL_FULL, &
                        options_vals%doubles%hole_selection%FAST_FAST, &
                        options_vals%doubles%weighting%INITIAL_DOBRAUTZ &
                    ) &
        ]

        do i = 1, size(csfs)
            do j = 1, size(options)
                call my_test_pgen(n_iter, csfs(i), options(j))
            end do
        end do


        call finalize_excitgen_test()

        contains
            logical elemental function is_problematic(analyzed_entry)
                type(AnalyzedEntries_t), intent(in) :: analyzed_entry
                is_problematic = analyzed_entry%p_gen > 1e-6_dp .and. abs(analyzed_entry%quality - 1._dp) > 0.03
            end function

#ifndef IFORT_
            logical pure function quality_comp(p1, p2)
                type(AnalyzedEntries_t), intent(in) :: p1, p2
                quality_comp = p1%quality >= p2%quality
            end function
#endif

            subroutine my_test_pgen(n_iter, csf_i, chosen_options)
                integer(int64), intent(in) :: n_iter
                type(CSF_Info_t), intent(in) :: csf_i
                type(PCHB_DoublesOptions_t), intent(in) :: chosen_options

                integer :: file_id
                type(GAS_PCHB_FastFast_t) :: orbital_picker
                type(Entry_t), allocatable :: entries(:)

                block
                    type(GAS_SuperGroupIndexer_t) :: indexer
                    indexer = get_CAS_indexer(n_el, nSpatOrbs)
                    call orbital_picker%init(indexer, &
                            chosen_options%particle_selection, &
                            chosen_options%weighting, .false., .false.)
                end block

                open(newunit=file_id, file=get_filename(csf_i, chosen_options))


                    write(file_id, *) 'Starting CSF: ', csf_i%to_str()
                    write(file_id, *) 'Chosen options: ', chosen_options%to_str()
                    write(file_id, *) 'n_iter: ', n_iter
                    write(file_id, *) '`calc_orbital_pgen_contr` with division of 2'

                    call test_pgen(orbital_picker, csf_i, n_iter, entries)
#ifdef IFORT_
                    associate(analyzed_entries => analyze_entries(entries, n_iter))
#else
                    associate(analyzed_entries => analyze_entries(entries, n_iter, quality_comp))
#endif
                        call output(file_id, pack(analyzed_entries, is_problematic(analyzed_entries)))
                    end associate

                close(file_id)

                call orbital_picker%finalize()
            end subroutine

            pure function get_filename(csf_i, options) result(res)
                type(CSF_Info_t), intent(in) :: csf_i
                type(PCHB_DoublesOptions_t), intent(in) :: options
                character(:), allocatable :: res

                res = to_lower(replace(options%to_str(), ' ', '_')) &
                        // '_' // &
                        join(csf_i%to_tokens(), '-') &
                        // '.dat'
            end function

    end subroutine

end module


program test_guga_gas_pgens

    use util_mod, only: stop_all
    use test_guga_gas_pgens_mod, only: guga_pchb_test_driver
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

end program
