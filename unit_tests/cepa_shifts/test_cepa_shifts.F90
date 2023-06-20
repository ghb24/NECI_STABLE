! replace _template by whatever one needs
program test_cepa_shifts

    use constants, only: dp
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case, assert_equals
    use cepa_shifts, only: cepa_0, cepa_aqcc, cepa_shift_double, &
        cepa_acpf, init_cepa_shifts, cepa_shift_single, &
        cepa_shift, cepa_method, aqcc_factor
    use CalcData, only: diagsft

    implicit none

    integer :: failed_count

    call init_fruit()
    call cepa_shifts_test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

    subroutine cepa_shifts_test_driver()

        call run_test_case(cepa_0_test, "cepa_0_test")
        call run_test_case(cepa_aqcc_test, "cepa_aqcc_test")
        call run_test_case(cepa_acpf_test, "cepa_acpf_test")
        call run_test_case(init_cepa_shifts_test, "init_cepa_shifts_test")
        call run_test_case(cepa_shift_test, "cepa_shift_test")

    end subroutine cepa_shifts_test_driver

    subroutine cepa_shift_test
        use replica_data, only: diagsft

        allocate(diagsft(2))
        diagsft(1) = 1.0
        diagsft(2) = 2.0

        print *, ""
        print *, "testing: cepa_shift: "
        cepa_shift_single => cepa_0
        cepa_shift_double => cepa_0

        call assert_equals(1.0_dp, cepa_shift(1,1))
        call assert_equals(2.0_dp, cepa_shift(2,2))

        call assert_equals(0.0_dp, cepa_shift(1,0))
        call assert_equals(0.0_dp, cepa_shift(1,3))
        call assert_equals(0.0_dp, cepa_shift(2,-1))
        call assert_equals(0.0_dp, cepa_shift(2,100))

        ! the rest is tested in init_cepa_shifts_test..

    end subroutine cepa_shift_test

    subroutine cepa_aqcc_test
        use SystemData, only: nel
        use replica_data, only: diagsft

        nel = 4
        allocate(diagsft(2))

        diagsft(1) = 1.0
        diagsft(2) = 2.0

        aqcc_factor = (1.0_dp - real((nel - 3)*(nel - 2),dp)/real(nel*(nel - 1), dp))

        print *, ""
        print *, "testing: cepa_aqcc_test: "

        call assert_equals(1.0_dp * (1.0_dp - aqcc_factor), cepa_aqcc(1))
        call assert_equals(2.0_dp * (1.0_dp - aqcc_factor), cepa_aqcc(2))

        nel = -1
        deallocate(diagsft)

    end subroutine cepa_aqcc_test

    subroutine cepa_acpf_test
        use Systemdata, only: nel
        use replica_data, only: diagsft

        nel = 10
        allocate(diagsft(2))
        diagsft(1) = 1.0
        diagsft(2) = 2.0

        print *, ""
        print *, "testing: cepa_acpf"

        call assert_equals(1.0_dp * (1.0_dp - 2.0_dp / 10.0_dp), cepa_acpf(1))
        call assert_equals(2.0_dp * (1.0_dp - 2.0_dp / 10.0_dp), cepa_acpf(2))

        nel = -1
        deallocate(diagsft)

    end subroutine cepa_acpf_test

    subroutine cepa_0_test
        print *, ""
        print *, "testing: cepa_0"

        allocate(diagsft(2))
        diagsft = 0.0_dp

        call assert_equals(0.0_dp, cepa_0(1))
        call assert_equals(0.0_dp, cepa_0(2))

        deallocate(diagsft)

    end subroutine cepa_0_test

    subroutine init_cepa_shifts_test
        use SystemData, only: nel
        use replica_data, only: diagsft
        use FciMCData, only: ilutref
        use DetBitOps, only: EncodeBitDet
        use bit_rep_data, only: niftot

        integer :: i

        nel = 10
        niftot = 0
        allocate(diagsft(2))
        diagsft(1) = 1.0
        diagsft(2) = 2.0

        allocate(ilutref(0:niftot,2))
        call EncodeBitDet([(i, i = 1,nel)], ilutref(:,1))
        call EncodeBitDet([(i, i = 1,nel)], ilutref(:,2))

        print *, ""
        print *, "testing init_cepa_shifts"

        cepa_method = '0'

        call init_cepa_shifts()

        call assert_equals(cepa_0(1), cepa_shift_single(1))
        call assert_equals(cepa_0(2), cepa_shift_single(2))

        call assert_equals(cepa_0(1), cepa_shift_double(1))
        call assert_equals(cepa_0(2), cepa_shift_double(2))

        call assert_equals(cepa_0(1), cepa_shift(1,1))
        call assert_equals(cepa_0(1), cepa_shift(1,2))
        call assert_equals(0.0_dp, cepa_shift(1,0))
        call assert_equals(0.0_dp, cepa_shift(2,0))

        call assert_equals(0.0_dp, cepa_shift(1,3))
        call assert_equals(0.0_dp, cepa_shift(2,100))

        cepa_method = 'acpf'
        call init_cepa_shifts()

        call assert_equals(cepa_acpf(1), cepa_shift_single(1))
        call assert_equals(cepa_acpf(2), cepa_shift_single(2))

        call assert_equals(cepa_acpf(1), cepa_shift(1,1))
        call assert_equals(cepa_acpf(2), cepa_shift(2,1))

        call assert_equals(cepa_acpf(1), cepa_shift_double(1))
        call assert_equals(cepa_acpf(2), cepa_shift_double(2))

        call assert_equals(cepa_acpf(1), cepa_shift(1,2))
        call assert_equals(cepa_acpf(2), cepa_shift(2,2))

        call assert_equals(0.0_dp, cepa_shift(1,0))
        call assert_equals(0.0_dp, cepa_shift(1,3))
        call assert_equals(0.0_dp, cepa_shift(2,0))
        call assert_equals(0.0_dp, cepa_shift(2,100))

        nel = 4
        cepa_method = 'aqcc'

        call init_cepa_shifts()

        call assert_equals(cepa_aqcc(1), cepa_shift_single(1))
        call assert_equals(cepa_aqcc(2), cepa_shift_single(2))

        call assert_equals(cepa_aqcc(1), cepa_shift(1,1))
        call assert_equals(cepa_aqcc(2), cepa_shift(2,1))

        call assert_equals(0.0_dp, cepa_shift(1,0))
        call assert_equals(0.0_dp, cepa_shift(1,3))

        call assert_equals(cepa_aqcc(1), cepa_shift_double(1))
        call assert_equals(cepa_aqcc(2), cepa_shift_double(2))

        call assert_equals(cepa_aqcc(1), cepa_shift(1,2))
        call assert_equals(cepa_aqcc(2), cepa_shift(2,2))

        call assert_equals(0.0_dp, cepa_shift(2,0))
        call assert_equals(0.0_dp, cepa_shift(2,100))

        nel = -1
        deallocate(diagsft)
        niftot = -1
        deallocate(ilutref)

    end subroutine init_cepa_shifts_test

end program test_cepa_shifts

