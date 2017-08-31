! replace _template by whatever one needs 
program test_cepa_shifts

    use fruit
    use cepa_shifts

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

    end subroutine cepa_shifts_test_driver

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

        call assert_equals(10.0/12.0, cepa_aqcc(1))
        call assert_equals(20.0/12.0, cepa_aqcc(2))

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

        call assert_equals(2.0/real(nel,dp), cepa_acpf(1))
        call assert_equals(4.0/real(nel,dp), cepa_acpf(2))

        nel = -1 
        deallocate(diagsft)

    end subroutine cepa_acpf_test

    subroutine cepa_0_test
        print *, "" 
        print *, "testing: cepa_0"

        call assert_equals(0.0, cepa_0(1))
        call assert_equals(0.0, cepa_0(-1))
        call assert_equals(0.0, cepa_0(0))
        call assert_equals(0.0, cepa_0(100))

    end subroutine cepa_0_test

    subroutine init_cepa_shifts_test
        use SystemData, only: nel
        use replica_data, only: diagsft

        nel = 10
        allocate(diagsft(1))
        diagsft(1) = 1.0
        diagsft(2) = 2.0

        print *, "" 
        print *, "testing init_cepa_shifts" 

        cepa_method = '0'

        call init_cepa_shifts()

        call assert_equals(0.0, cepa_shift_single(0))
        call assert_equals(0.0, cepa_shift_single(-1))
        call assert_equals(0.0, cepa_shift_single(1))
        call assert_equals(0.0, cepa_shift_single(1000))

        call assert_equals(0.0, cepa_shift_double(0))
        call assert_equals(0.0, cepa_shift_double(-1))
        call assert_equals(0.0, cepa_shift_double(1))
        call assert_equals(0.0, cepa_shift_double(1000))

        cepa_method = 'acpf'
        call init_cepa_shifts()

        call assert_equals(2.0/real(nel,dp), cepa_shift_single(1))
        call assert_equals(4.0/real(nel,dp), cepa_shift_single(2))

        call assert_equals(2.0/real(nel,dp), cepa_shift_double(1))
        call assert_equals(4.0/real(nel,dp), cepa_shift_double(2))

        nel = 4 
        cepa_method = 'aqcc'

        call init_cepa_shifts()

        call assert_equals(10.0/12.0, cepa_shift_single(1))
        call assert_equals(20.0/12.0, cepa_shift_single(2))

        call assert_equals(10.0/12.0, cepa_shift_double(1))
        call assert_equals(20.0/12.0, cepa_shift_double(2))

    end subroutine init_cepa_shifts_test

end program test_cepa_shifts

