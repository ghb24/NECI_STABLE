module test_guga_matrixelements_mod
    use fruit, only: assert_true, assert_false, assert_equals, run_test_case
    use constants, only: dp, n_int, maxExcit, stdout
    use guga_bitRepOps, only: convert_stepvector_to_nI, StepVectorVals_t, StepVector_t, stepvector_vals, &
        EncodeBitDet_guga, CSF_Info_t
    use SystemData, only: nEl, nSpatOrbs, nBasis
    use guga_matrixElements, only: calc_guga_matrix_element
    use guga_data, only: excit_type, gen_type, ExcitationInformation_t
    use bit_rep_data, only: GugaBits
    use guga_plugin, only: init_guga_plugin
    implicit none
    private
    public :: test_driver

    type(StepVector_t), parameter :: s0 = stepvector_vals%s0, u = stepvector_vals%u, &
        d = stepvector_vals%d, s2 = stepvector_vals%s2
contains

    subroutine test_driver()
        call run_test_case(test_input_stepvector, "test_input_stepvector")
    end subroutine


    subroutine test_input_stepvector()

        call init_guga_plugin('FCIDUMP', t_testmode_=.true., nel_=2, nbasis_=8, &
                              nSpatOrbs_=4)

        block
            type(StepVector_t), parameter :: step_v(9) = [s2, u, u, d, u, u, s0, d, s0]
            call assert_true(all(step_v%val == [3, 1, 1, 2, 1, 1, 0, 2, 0]))
        end block

        write(stdout, *) get_coupling_coefficients([s0, u, s0, u], [u, s0, u, s0])
    end subroutine

    pure function get_coupling_coefficients(stepv_i, stepv_j) result(res)
        type(StepVector_t), intent(in) :: stepv_i(:), stepv_j(:)
        real(dp) :: res

        type(ExcitationInformation_t) :: excitInfo
        integer(n_int) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)

        ilutI = get_ilut(stepv_i)
        ilutJ = get_ilut(stepv_j)
        call calc_guga_matrix_element(ilutI, CSF_Info_t(ilutI), ilutJ, CSF_Info_t(ilutJ), excitInfo, res, t_hamil=.false.)
    end function

    pure function get_ilut(stepv) result(ilut)
        type(StepVector_t), intent(in) :: stepv(:)
        integer(n_int) :: ilut(0:GugaBits%len_tot)
        call EncodeBitDet_guga(convert_stepvector_to_nI(stepv), ilut)
    end function
end module

program test_guga_matrixelements_prog

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_guga_matrixelements_mod, only: test_driver
    use util_mod, only: stop_all

    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

end program
