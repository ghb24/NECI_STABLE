#include "macros.h" 

program test_tJ_model 

    use tJ_model 
    use fruit 

    implicit none 

    integer :: failed_count 

    call init_fruit() 
    call tJ_model_test_driver() 
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1 

contains 

    subroutine tJ_model_test_driver() 

        call run_test_case(init_tJ_model_test, 'init_tJ_model_test')
        call run_test_case(gen_excit_tJ_model_test, "gen_excit_tJ_model_test")
        call run_test_case(create_cum_list_tJ_model_test, "create_cum_list_tJ_model_test")
        call run_test_case(find_elec_in_ni_test, "find_elec_in_ni_test")
        call run_test_case(calc_pgen_tJ_model_test, "calc_pgen_tJ_model_test")
        call run_test_case(init_heisenberg_model_test, "init_heisenberg_model_test")
        call run_test_case(gen_excit_heisenberg_model_test, "gen_excit_heisenberg_model_test")
        call run_test_case(create_cum_list_heisenberg_test, "create_cum_list_heisenberg_test")
        call run_test_case(calc_pgen_heisenberg_model_test, "calc_pgen_heisenberg_model_test")
        call run_test_case(setup_exchange_matrix_test, "setup_exchange_matrix_test")
        call run_test_case(get_helement_tJ_test, "get_helement_tJ_test")
        call run_test_case(get_helement_heisenberg_test, "get_helement_heisenberg_test")
        call run_test_case(get_diag_helement_heisenberg_test, "get_diag_helement_heisenberg_test")
        call run_test_case(get_offdiag_helement_heisenberg_test, "get_offdiag_helement_heisenberg_test")
        call run_test_case(determine_optimal_time_step_tJ_test, "determine_optimal_time_step_tJ_test")
        call run_test_case(determine_optimal_time_step_heisenberg_test, "determine_optimal_time_step_heisenberg_test")
        call run_test_case(get_umat_el_heisenberg_test, "get_umat_el_heisenberg_test")

    end subroutine tJ_model_test_driver

    subroutine init_tJ_model_test

        print *, "" 
        print *, "testing: init_tJ_model: "
        call assert_true(.false.)

    end subroutine init_tJ_model_test

    subroutine init_heisenberg_model_test

        print *, "" 
        print *, "testing: init_heisenberg_model: "
        call assert_true(.false.)

    end subroutine init_heisenberg_model_test

    subroutine gen_excit_tJ_model_test

        print *, "" 
        print *, "testing: gen_excit_tJ_model"
        call assert_true(.false.)

    end subroutine gen_excit_tJ_model_test

    subroutine calc_pgen_tJ_model_test

        print *, ""
        print *, "testing: calc_pgen_tJ_model"

        call assert_true(.false.)

    end subroutine calc_pgen_tJ_model_test

    subroutine create_cum_list_tJ_model_test

        use SystemData, only: nel, bhub
        use lattice_mod, only: lattice
        use real_space_hubbard, only: lat 
        use bit_rep_data, only: NIfTot
        use constants, only: dp, n_int 

        integer(n_int), allocatable :: ilut(:) 
        integer, allocatable :: ic_list(:)
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum 

        nel = 2 
        NIfTot = 0
        nbasis = 8 

        bhub = -1.0
        exchange_j = 1.0 

        allocate(ilut(0:NIfTot)) 
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)

        ! i also have to setup the matrix element calculation.. duh.. 

        print *, ""
        print *, "testing: create_cum_list_tJ_model"

        call create_cum_list_tJ_model(ilut, 1, [2,4], & 
            cum_arr, cum_sum, ic_list)

        call assert_true(.false.)
        nel = -1 
        NIfTot = -1 
        nbasis = -1 

    end subroutine create_cum_list_tJ_model_test

    subroutine gen_excit_heisenberg_model_test

        print *, ""
        print *, "testing: gen_excit_heisenberg_model"
        call assert_true(.false.)

    end subroutine gen_excit_heisenberg_model_test

    subroutine create_cum_list_heisenberg_test

        print *, ""
        print *, "testing: create_cum_list_heisenberg"
        call assert_true(.false.)

    end subroutine create_cum_list_heisenberg_test

    subroutine calc_pgen_heisenberg_model_test

        print *, ""
        print *, "testing: calc_pgen_heisenberg_model"
        call assert_true(.false.)

    end subroutine calc_pgen_heisenberg_model_test

    subroutine find_elec_in_ni_test
        use SystemData, only: nel, nbasis
        
        print *, ""
        print *, "testing: find_elec_in_ni" 

        nel = 3
        nbasis = 8 
        call assert_equals(3, find_elec_in_ni([1,2,3],3))
        call assert_equals(2, find_elec_in_ni([1,2,3],2))
        call assert_equals(1, find_elec_in_ni([1,2,3],1))

        call assert_equals(-1, find_elec_in_ni([1,2,3],4))

        call assert_equals(2, find_elec_in_ni([3,7,8],7))
        call assert_equals(1, find_elec_in_ni([3,7,8],3))
        call assert_equals(3, find_elec_in_ni([3,7,8],8))

        call assert_equals(-1, find_elec_in_ni([3,7,8],1))
        call assert_equals(-1, find_elec_in_ni([3,7,8],4))
        call assert_equals(-1, find_elec_in_ni([3,7,8],5))

        nel = -1
        nbasis = -1

    end subroutine find_elec_in_ni_test

    subroutine setup_exchange_matrix_test

        print *, "" 
        print *, "testing: setup_exchange_matrix"
        call assert_true(.false.)

    end subroutine setup_exchange_matrix_test

    subroutine get_helement_tJ_test

        print *, "" 
        print *, "testing: get_helement_tJ" 

        call assert_true(.false.)
        print *, ""
        print *, " for both module procdeures"
        call assert_true(.false.)
    end subroutine get_helement_tJ_test

    subroutine get_helement_heisenberg_test

        print *, "" 
        print *, "testing: get_helement_heisenberg" 

        call assert_true(.false.)
        print *, ""
        print *, " for both module procdeures"
        call assert_true(.false.)
    end subroutine get_helement_heisenberg_test

    subroutine get_diag_helement_heisenberg_test

        print *, ""
        print *, "testing: get_diag_helement_heisenberg" 
        call assert_true(.false.)

    end subroutine get_diag_helement_heisenberg_test

    subroutine get_offdiag_helement_heisenberg_test

        print *, "" 
        print *, "testing: get_offdiag_helement_heisenberg"
        call assert_true(.false.)

    end subroutine get_offdiag_helement_heisenberg_test
    
    subroutine determine_optimal_time_step_tJ_test 

        print *, ""
        print *, "testing: determine_optimal_time_step_tJ"
        call assert_true(.false.) 

    end subroutine determine_optimal_time_step_tJ_test 

    subroutine determine_optimal_time_step_heisenberg_test 

        print *, "" 
        print *, "testing: determine_optimal_time_step_heisenberg" 
        call assert_true(.false.) 

    end subroutine determine_optimal_time_step_heisenberg_test 

    subroutine get_umat_el_heisenberg_test

        print *, "" 
        print *, "testing: get_umat_el_heisenberg"
        call assert_true(.false.)
    end subroutine get_umat_el_heisenberg_test

end program test_tJ_model
