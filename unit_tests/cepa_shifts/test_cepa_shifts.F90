! replace _template by whatever one needs 
program test_cepa_shifts

    use fruit
    use cepa_shifts
    use cc_amplitudes
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
        call run_test_case(calc_n_single_excits_test, "calc_n_single_excits_test")
        call run_test_case(binomial_test, "binomial_test")
        call run_test_case(calc_n_parallel_excitations_test, "calc_n_parallel_excitations_test")
        call run_test_case(calc_number_of_excitations_test, "calc_number_of_excitations_test")

    end subroutine cepa_shifts_test_driver

    subroutine calc_number_of_excitations_test

        print *, ""
        print *, "calc_number_of_excitations" 

        call assert_equals([8,18], calc_number_of_excitations(2,2,2,4),2)
        call assert_equals([8,18,8], calc_number_of_excitations(2,2,3,4),2)
        call assert_equals([8,18,8,1], calc_number_of_excitations(2,2,4,4),2)

        call assert_equals([7,13,3,0], calc_number_of_excitations(1,2,4,4), 2)
        call assert_equals([7,13,3,0], calc_number_of_excitations(2,1,4,4), 2)

        call assert_equals([9,9,1,0], calc_number_of_excitations(0,3,4,6), 2)
        
        call assert_equals([10,27,12,0], calc_number_of_excitations(1,3,4,5), 2)
        call assert_equals([10,27,12,0], calc_number_of_excitations(3,1,4,5), 2)

        call assert_equals([12,42,36,9,0], calc_number_of_excitations(2,3,5,5), 2)

    end subroutine calc_number_of_excitations_test

    subroutine calc_n_parallel_excitations_test

        print *, ""
        print *, "testing: calc_n_parallel_excitations "

        call assert_equals(4, calc_n_parallel_excitations(2,4,1))
        call assert_equals(6, calc_n_parallel_excitations(2,5,1))
        call assert_equals(2, calc_n_parallel_excitations(2,3,1))
        call assert_equals(1, calc_n_parallel_excitations(1,2,1))
        call assert_equals(3, calc_n_parallel_excitations(3,4,1))

        call assert_equals(1, calc_n_parallel_excitations(2,4,2))
        call assert_equals(3, calc_n_parallel_excitations(2,5,2))

    end subroutine calc_n_parallel_excitations_test

    subroutine binomial_test 
        print *, ""
        print *, "testing: binomial" 
        
        call assert_equals(1, binomial(1,1))
        call assert_equals(1, binomial(2,0))
        call assert_equals(1, binomial(3,0))
        call assert_equals(1, binomial(4,0))
        call assert_equals(1, binomial(2,0))
        call assert_equals(1, binomial(3,0))
        call assert_equals(1, binomial(4,0))

        call assert_equals(2, binomial(2,1))
        call assert_equals(3, binomial(3,1))
        call assert_equals(4, binomial(4,1))

        call assert_equals(10, binomial(5,2))
        call assert_equals(10, binomial(5,3))

        call assert_equals(20, binomial(6,3))
        call assert_equals(35, binomial(7,4))

        call assert_equals(1, binomial(0,0))
        call assert_equals(0, binomial(0,1))
        call assert_equals(0, binomial(0,100))
        call assert_equals(0, binomial(0,2))

    end subroutine binomial_test 

    subroutine calc_n_single_excits_test

        print *, ""
        print *, "testing: calc_n_single_excits "

        call assert_equals(4, calc_n_single_excits(2,4))
        call assert_equals(6, calc_n_single_excits(2,5))
        call assert_equals(2, calc_n_single_excits(2,3))
        call assert_equals(1, calc_n_single_excits(1,2))
        call assert_equals(3, calc_n_single_excits(3,4))

    end subroutine calc_n_single_excits_test

    subroutine cepa_shift_test
        use replica_data, only: diagsft

        allocate(diagsft(2))
        diagsft(1) = 1.0
        diagsft(2) = 2.0

        print *, "" 
        print *, "testing: cepa_shift: "
        cepa_shift_single => cepa_0
        cepa_shift_double => cepa_0

        call assert_equals(0.0_dp, cepa_shift(1,1))
        call assert_equals(0.0_dp, cepa_shift(2,2))
        
        call assert_equals(1.0_dp, cepa_shift(1,0))
        call assert_equals(1.0_dp, cepa_shift(1,3))
        call assert_equals(2.0_dp, cepa_shift(2,-1))
        call assert_equals(2.0_dp, cepa_shift(2,100))

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

        call assert_equals(10.0_dp/12.0_dp, cepa_aqcc(1))
        call assert_equals(20.0_dp/12.0_dp, cepa_aqcc(2))

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

        call assert_equals(2.0_dp/real(nel,dp), cepa_acpf(1))
        call assert_equals(4.0_dp/real(nel,dp), cepa_acpf(2))

        nel = -1 
        deallocate(diagsft)

    end subroutine cepa_acpf_test

    subroutine cepa_0_test
        print *, "" 
        print *, "testing: cepa_0"

        call assert_equals(0.0_dp, cepa_0(1))
        call assert_equals(0.0_dp, cepa_0(-1))
        call assert_equals(0.0_dp, cepa_0(0))
        call assert_equals(0.0_dp, cepa_0(100))

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

        allocate(ilutref(0:niftot,1:2))
        call EncodeBitDet([(i, i = 1,nel)], ilutref(:,1))
        call EncodeBitDet([(i, i = 1,nel)], ilutref(:,2))

        print *, "" 
        print *, "testing init_cepa_shifts" 

        cepa_method = '0'

        call init_cepa_shifts()

        call assert_equals(0.0_dp, cepa_shift_single(0))
        call assert_equals(0.0_dp, cepa_shift_single(-1))
        call assert_equals(0.0_dp, cepa_shift_single(1))
        call assert_equals(0.0_dp, cepa_shift_single(1000))

        call assert_equals(0.0_dp, cepa_shift_double(0))
        call assert_equals(0.0_dp, cepa_shift_double(-1))
        call assert_equals(0.0_dp, cepa_shift_double(1))
        call assert_equals(0.0_dp, cepa_shift_double(1000))

        call assert_equals(0.0_dp, cepa_shift(1,1)) 
        call assert_equals(0.0_dp, cepa_shift(1,2)) 
        call assert_equals(1.0_dp, cepa_shift(1,0)) 
        call assert_equals(2.0_dp, cepa_shift(2,0)) 

        call assert_equals(1.0_dp, cepa_shift(1,3)) 
        call assert_equals(2.0_dp, cepa_shift(2,100)) 

        cepa_method = 'acpf'
        call init_cepa_shifts()

        call assert_equals(2.0_dp/real(nel,dp), cepa_shift_single(1))
        call assert_equals(4.0_dp/real(nel,dp), cepa_shift_single(2))

        call assert_equals(2.0_dp/real(nel,dp), cepa_shift(1,1))
        call assert_equals(4.0_dp/real(nel,dp), cepa_shift(2,1))

        call assert_equals(2.0_dp/real(nel,dp), cepa_shift_double(1))
        call assert_equals(4.0_dp/real(nel,dp), cepa_shift_double(2))

        call assert_equals(2.0_dp/real(nel,dp), cepa_shift(1,2))
        call assert_equals(4.0_dp/real(nel,dp), cepa_shift(2,2))

        call assert_equals(1.0_dp, cepa_shift(1,0))
        call assert_equals(1.0_dp, cepa_shift(1,3))
        call assert_equals(2.0_dp, cepa_shift(2,0))
        call assert_equals(2.0_dp, cepa_shift(2,100))

        nel = 4 
        cepa_method = 'aqcc'

        call init_cepa_shifts()

        call assert_equals(10.0_dp/12.0, cepa_shift_single(1))
        call assert_equals(20.0_dp/12.0, cepa_shift_single(2))

        call assert_equals(10.0_dp/12.0, cepa_shift(1,1))
        call assert_equals(20.0_dp/12.0, cepa_shift(2,1))

        call assert_equals(1.0_dp, cepa_shift(1,0))
        call assert_equals(1.0_dp, cepa_shift(1,3))

        call assert_equals(10.0_dp/12.0, cepa_shift_double(1))
        call assert_equals(20.0_dp/12.0, cepa_shift_double(2))

        call assert_equals(10.0_dp/12.0, cepa_shift(1,2))
        call assert_equals(20.0_dp/12.0, cepa_shift(2,2))

        call assert_equals(2.0_dp, cepa_shift(2,0))
        call assert_equals(2.0_dp, cepa_shift(2,100))

        nel = -1 
        deallocate(diagsft)
        niftot = -1
        deallocate(projedet)

    end subroutine init_cepa_shifts_test

end program test_cepa_shifts

