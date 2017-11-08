#include "macros.h" 

program test_tJ_model 

    use tJ_model 
    use fruit 

    implicit none 

    integer :: failed_count 

    t_tJ_model = .true. 

    call init_fruit() 
    call tJ_model_test_driver() 
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1 

contains 

    subroutine tJ_model_test_driver() 

        call run_test_case(init_tJ_model_test, 'init_tJ_model_test')
        call run_test_case(init_heisenberg_model_test, "init_heisenberg_model_test")
        call run_test_case(setup_exchange_matrix_test, "setup_exchange_matrix_test")
        call run_test_case(find_elec_in_ni_test, "find_elec_in_ni_test")
        call run_test_case(get_umat_el_heisenberg_test, "get_umat_el_heisenberg_test")
        call run_test_case(get_helement_tJ_test, "get_helement_tJ_test")
        call run_test_case(get_helement_heisenberg_test, "get_helement_heisenberg_test")
        call run_test_case(get_diag_helement_heisenberg_test, "get_diag_helement_heisenberg_test")
        call run_test_case(get_offdiag_helement_heisenberg_test, "get_offdiag_helement_heisenberg_test")
        call run_test_case(determine_optimal_time_step_tJ_test, "determine_optimal_time_step_tJ_test")
        call run_test_case(determine_optimal_time_step_heisenberg_test, "determine_optimal_time_step_heisenberg_test")
        call run_test_case(create_cum_list_tJ_model_test, "create_cum_list_tJ_model_test")
        call run_test_case(create_cum_list_heisenberg_test, "create_cum_list_heisenberg_test")
        call run_test_case(gen_excit_tJ_model_test, "gen_excit_tJ_model_test")
        call run_test_case(calc_pgen_tJ_model_test, "calc_pgen_tJ_model_test")
        call run_test_case(gen_excit_heisenberg_model_test, "gen_excit_heisenberg_model_test")
        call run_test_case(calc_pgen_heisenberg_model_test, "calc_pgen_heisenberg_model_test")

    end subroutine tJ_model_test_driver

    subroutine init_tJ_model_test
        use SystemData, only: lattice_type, length_x, length_y, nbasis, nel, nbasis
        use OneEInts, only: tmat2d 
        use FciMCData, only: tsearchtau, tsearchtauoption, ilutref
        use CalcData, only: tau 
        use tau_search, only: max_death_cpt
        use procedure_pointers, only: get_umat_el, generate_excitation
        use real_space_hubbard, only: lat_tau_factor
        use bit_rep_data, only: nifd, NIfTot

        print *, "" 
        print *, "testing: init_tJ_model: "
        lattice_type = 'chain' 
        length_x = 100
        length_y = 1 
        nel = 2 
        exchange_j = 1
        nifd = 0
        NIfTot = 0
        allocate(ilutref(0:NIfTot,1))
        ilutref = 9

        call init_tJ_model() 

        call assert_equals(lat%get_ndim(), 1) 
        call assert_equals(lat%get_nsites(), 100)
        call assert_equals(lat%get_length(), 100) 
        call assert_true(lat%is_periodic()) 
        call assert_equals(2, lat%get_nconnect_max())

        call assert_equals(lat%get_site_index(1), 1) 
        call assert_equals(lat%get_site_index(100), 100) 
        call assert_equals(lat%get_site_index(50), 50) 
        call assert_equals(lat%get_neighbors(1), [100, 2], size(lat%get_neighbors(1)))
        call assert_equals(lat%get_neighbors(2), [1,3], size(lat%get_neighbors(2)))
        call assert_equals(lat%get_neighbors(100), [99,1], size(lat%get_neighbors(2)))

        call assert_equals(2, lat%get_num_neighbors(1))
        call assert_equals(2, lat%get_num_neighbors(2))
        call assert_equals(2, lat%get_num_neighbors(100))
        ! i actually do not need to have the common lattice type or? 
        ! when i want to test a chain i could just use a chain or? 

        call assert_equals(200, nbasis) 
        call assert_true(associated(tmat2d))
        call assert_true(associated(g1))
        call assert_equals(0.0_dp, ecore)
        call assert_true(.not. tsearchtau)
        call assert_true(tsearchtauoption)
        call assert_equals(0.0_dp, max_death_cpt)
        call assert_true(associated(get_umat_el))
        call assert_true(associated(generate_excitation))
        call assert_equals(0.25 * lat_tau_factor, tau)

        length_x = -1 
        length_y = -1 
        nel = -1 
        nifd = -1
        NIfTot = -1 
        deallocate(ilutref)
        nbasis = 0

    end subroutine init_tJ_model_test

    subroutine init_heisenberg_model_test
        use SystemData, only: lattice_type, length_x, length_y, nbasis, nel, & 
                              ecore
        use OneEInts, only: tmat2d 
        use FciMCData, only: tsearchtau, tsearchtauoption 
        use CalcData, only: tau 
        use tau_search, only: max_death_cpt
        use procedure_pointers, only: get_umat_el, generate_excitation
        use real_space_hubbard, only: lat_tau_factor
        use bit_rep_data, only: nifd, NIfTot 

        lattice_type = 'chain' 
        length_x = 3
        length_y = 1 
        nel = 3 
        exchange_j = 1
        nifd = 0
        NIfTot = 0
        allocate(ilutref(0:NIfTot,1))
        ilutref = 25

        print *, "" 
        print *, "testing: init_heisenberg_model: "
        call init_heisenberg_model()

        call assert_equals(lat%get_ndim(), 1) 
        call assert_equals(lat%get_nsites(), 3)
        call assert_equals(lat%get_length(), 3) 
        call assert_true(lat%is_periodic()) 
        call assert_equals(2, lat%get_nconnect_max())

        call assert_equals(lat%get_site_index(1), 1) 
        call assert_equals(lat%get_site_index(2), 2) 
        call assert_equals(lat%get_neighbors(1), [3, 2], size(lat%get_neighbors(1)))
        call assert_equals(lat%get_neighbors(2), [1,3], size(lat%get_neighbors(2)))
        call assert_equals(lat%get_neighbors(3), [1,2], size(lat%get_neighbors(2)))

        call assert_equals(2, lat%get_num_neighbors(1))
        call assert_equals(2, lat%get_num_neighbors(2))
        ! i actually do not need to have the common lattice type or? 
        ! when i want to test a chain i could just use a chain or? 

        call assert_equals(6, nbasis) 
        call assert_true(associated(tmat2d))
        call assert_true(associated(g1))
        call assert_equals(0.0_dp, ecore)
        call assert_true(.not. tsearchtau)
        call assert_true(tsearchtauoption)
        call assert_equals(0.0_dp, max_death_cpt)
        call assert_true(associated(get_umat_el))
        call assert_true(associated(generate_excitation))
        call assert_equals(0.25 * lat_tau_factor, tau)

        length_x = -1 
        length_y = -1 
        nel = -1 
        nifd = -1
        NIfTot = -1 
        deallocate(ilutref)
        nbasis = 0

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
        use SystemData, only: nbasis 
        use real_space_hubbard, only: lat
        use lattice_mod, only: lattice

        print *, "" 
        print *, "testing: setup_exchange_matrix"
        exchange_j = 1.0 
        lat => lattice('chain', 2, 1, 1, .true.,.true.,.true.)
        nbasis = 4 
        call setup_exchange_matrix(lat)
        call assert_equals([h_cast(0.0),h_cast(0.0),h_cast(0.0),h_cast(1.0)],exchange_matrix(1,:),4)
        call assert_equals([h_cast(0.0),h_cast(0.0),h_cast(1.0),h_cast(0.0)],exchange_matrix(2,:),4)
        call assert_equals([h_cast(0.0),h_cast(1.0),h_cast(0.0),h_cast(0.0)],exchange_matrix(3,:),4)
        call assert_equals([h_cast(1.0),h_cast(0.0),h_cast(0.0),h_cast(0.0)],exchange_matrix(4,:),4)

        exchange_j = 0 
        nbasis = -1 

    end subroutine setup_exchange_matrix_test

    subroutine get_helement_tJ_test
        use SystemData, only: nel, bhub, nbasis 
        use Determinants, only: get_helement 
        use procedure_pointers, only: get_umat_el
        use lattice_mod, only: lattice 
        use real_space_hubbard, only: lat 
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: ilutI(:), ilutJ(:)

        NIfTot = 0
        nifd = 0
        nbasis = 4 
        get_umat_el => get_umat_el_heisenberg

        nel = 2 
        bhub = 2.0 
        exchange_j = 1.0 


        lat => lattice('chain', 2, 1, 1, .false.,.false.,.false.) 
        call init_tmat(lat) 
        call setup_exchange_matrix(lat) 

        print *, "" 
        print *, "testing: get_helement_tJ" 
        print *, "first test hubbard like single excitations: " 

        call assert_equals(h_cast(2.0_dp), get_helement([1,2],[1,4],1))
        call assert_equals(h_cast(-2.0_dp), get_helement([1,2],[2,3],1))

        call assert_equals(h_cast(2.0_dp), get_helement([1,2],[1,4]))
        call assert_equals(h_cast(-2.0_dp), get_helement([1,2],[2,3]))

        ! what about the diagonal elements? although this is a heisenberg 
        ! like setup.. 
        print *, "" 
        print *, "for diagonal elements: "
        call assert_equals(h_cast(0.0), get_helement([1,2],[1,2],0))
        call assert_equals(h_cast(0.0), get_helement([1,2],[1,2]))

        ! but i am not yet sure about the double counting.. 
        call assert_equals(h_cast(-0.5), get_helement([1,4],[1,4],0))
        call assert_equals(h_cast(-0.5), get_helement([2,3],[2,3]))

        call assert_equals(h_cast(0.0), get_helement([1,3],[1,3]))
        call assert_equals(h_cast(0.0), get_helement([2,4],[2,4]))

        print *, "" 
        print *, "for exchange contributions: "
        call assert_equals(h_cast(1.0), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0), get_helement([2,3],[1,4]))

        print *, "" 
        print *, "and for bigger systems: "
        nbasis = 6 
        lat => lattice('chain', 3, 1, 1, .true., .true., .true.)
        call init_tmat(lat) 
        call setup_exchange_matrix(lat) 

        call assert_equals(h_cast(0.0), get_helement([1,3],[1,3],0))
        call assert_equals(h_cast(0.0), get_helement([1,3],[1,3]))

        call assert_equals(h_cast(2.0), get_helement([1,3],[1,5],1))
        call assert_equals(h_cast(-2.0), get_helement([1,3],[3,5],1))

        call assert_equals(h_cast(2.0), get_helement([1,3],[1,5]))
        call assert_equals(h_cast(-2.0), get_helement([1,3],[3,5]))

        niftot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet([1,2],ilutI)
        call encodebitdet([1,2],ilutJ)

        call assert_equals(h_cast(0.0), get_helement([1,2],[1,2],ilutI,ilutJ))

        call encodebitdet([2,4],ilutJ)
        call assert_equals(h_cast(0.0), get_helement([1,2],[2,4],ilutI,ilutJ))

        call encodebitdet([2,3],ilutJ)
        call assert_equals(h_cast(-2.0), get_helement([1,2],[2,3],ilutI,ilutJ))

        call encodebitdet([1,4],ilutJ)
        call assert_equals(h_cast(2.0), get_helement([1,2],[1,4],ilutJ,ilutJ))

        call assert_equals(h_cast(1.0), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0), get_helement([1,6],[2,5]))

        call assert_equals(h_cast(-2.0), get_helement([1,4],[4,5]))
        call assert_equals(h_cast(2.0), get_helement([1,4],[1,6]))

        call assert_equals(h_cast(0.0), get_helement([1,4],[3,6]))

        call assert_equals(h_cast(0.0), get_helement([1,3],[2,4]))
        call assert_equals(h_cast(0.0), get_helement([1,3],[2,3]))

        lat => lattice('square', 2,2,1,.true.,.true.,.true.)
        nel = 4 
        nbasis = 8

        call init_tmat(lat) 
        call setup_exchange_matrix(lat) 

        call assert_equals(h_cast(0.0), get_helement([1,3,5,7],[1,3,5,7]))
        call assert_equals(h_cast(0.0), get_helement([2,4,6,8],[2,4,6,8]))

        call assert_equals(h_cast(-1.0), get_helement([1,3,6,8],[1,3,6,8]))
        call assert_equals(h_cast(-1.0), get_helement([2,4,5,7],[2,4,5,7]))

        call assert_equals(h_cast(-1.0), get_helement([1,4,5,8],[1,4,5,8]))

        call assert_equals(h_cast(-2.0), get_helement([1,4,6,7],[1,4,6,7]))

        NIfTot = -1 
        nifd = -1
        nbasis = -1 
        nel = -1 

    end subroutine get_helement_tJ_test

    subroutine get_helement_heisenberg_test
        use SystemData, only: nel, bhub, nbasis 
        use Determinants, only: get_helement 
        use procedure_pointers, only: get_umat_el
        use lattice_mod, only: lattice 
        use real_space_hubbard, only: lat 
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: ilutI(:), ilutJ(:)

        NIfTot = 0
        nifd = 0
        nbasis = 4 
        get_umat_el => get_umat_el_heisenberg

        nel = 2 
        bhub = 2.0 
        exchange_j = 1.0 

        lat => lattice('chain', 2, 1, 1, .false.,.false.,.false.) 
        call setup_exchange_matrix(lat) 

        t_tJ_model = .false. 
        t_heisenberg_model = .true. 

        print *, "" 
        print *, "testing: get_helement_heisenberg" 

        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,4],1))
        call assert_equals(h_cast(-0.0_dp), get_helement([1,2],[2,3],1))

        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,4]))
        call assert_equals(h_cast(-0.0_dp), get_helement([1,2],[2,3]))

        print *, "" 
        print *, "for diagonal elements: "
        call assert_equals(h_cast(0.0), get_helement([1,2],[1,2],0))
        call assert_equals(h_cast(0.0), get_helement([1,2],[1,2]))

        ! but i am not yet sure about the double counting.. 
        call assert_equals(h_cast(-0.25), get_helement([1,4],[1,4],0))
        call assert_equals(h_cast(-0.25), get_helement([2,3],[2,3]))

        call assert_equals(h_cast(0.25), get_helement([1,3],[1,3]))
        call assert_equals(h_cast(0.25), get_helement([2,4],[2,4]))

        print *, "" 
        print *, "for exchange contributions: "
        call assert_equals(h_cast(1.0), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0), get_helement([2,3],[1,4]))

        print *, "" 
        print *, "and for bigger systems: "
        nbasis = 6 
        lat => lattice('chain', 3, 1, 1, .true., .true., .true.)
        call setup_exchange_matrix(lat) 

        call assert_equals(h_cast(0.25), get_helement([1,3],[1,3],0))
        call assert_equals(h_cast(0.25), get_helement([1,3],[1,3]))

        call assert_equals(h_cast(0.0), get_helement([1,3],[1,5],1))
        call assert_equals(h_cast(-0.0), get_helement([1,3],[3,5],1))

        call assert_equals(h_cast(0.0), get_helement([1,3],[1,5]))
        call assert_equals(h_cast(-0.0), get_helement([1,3],[3,5]))

        niftot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet([1,2],ilutI)
        call encodebitdet([1,2],ilutJ)

        call assert_equals(h_cast(0.0), get_helement([1,2],[1,2],ilutI,ilutJ))

        call encodebitdet([2,4],ilutJ)
        call assert_equals(h_cast(0.0), get_helement([1,2],[2,4],ilutI,ilutJ))

        call encodebitdet([2,3],ilutJ)
        call assert_equals(h_cast(-0.0), get_helement([1,2],[2,3],ilutI,ilutJ))

        call encodebitdet([1,4],ilutJ)
        call assert_equals(h_cast(0.0), get_helement([1,2],[1,4],ilutJ,ilutJ))

        call assert_equals(h_cast(1.0), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0), get_helement([1,6],[2,5]))

        call assert_equals(h_cast(-0.0), get_helement([1,4],[4,5]))
        call assert_equals(h_cast(0.0), get_helement([1,4],[1,6]))

        call assert_equals(h_cast(0.0), get_helement([1,4],[3,6]))

        call assert_equals(h_cast(0.0), get_helement([1,3],[2,4]))
        call assert_equals(h_cast(0.0), get_helement([1,3],[2,3]))

        NIfTot = -1 
        nifd = -1
        nbasis = -1 
        nel = -1 

    end subroutine get_helement_heisenberg_test

    subroutine get_diag_helement_heisenberg_test
        use SystemData, only: nel, nbasis
        use bit_rep_data, only: NIfTot
        use real_space_hubbard, only: lat 
        
        nel = 4 
        lat => lattice('square', 2, 2, 1, .true., .true., .true.) 
        nbasis = 8
        NIfTot = 0
        exchange_j = 1.0

        call setup_exchange_matrix(lat) 

        print *, ""
        print *, "testing: get_diag_helement_heisenberg" 
        t_tJ_model = .true. 
        t_heisenberg_model = .false. 

        call assert_equals(h_cast(0.0), get_diag_helement_heisenberg([1,3,5,7]))
        call assert_equals(h_cast(0.0), get_diag_helement_heisenberg([2,4,6,8]))

        call assert_equals(h_cast(0.0), get_diag_helement_heisenberg([1,2,3,4]))
        call assert_equals(h_cast(-1.0), get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-2.0), get_diag_helement_heisenberg([1,4,6,7]))

        t_tJ_model = .false. 
        t_heisenberg_model = .true. 

        call assert_equals(h_cast(1.0), get_diag_helement_heisenberg([1,3,5,7]))
        call assert_equals(h_cast(1.0), get_diag_helement_heisenberg([2,4,6,8]))

        call assert_equals(h_cast(-0.0), get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-1.0), get_diag_helement_heisenberg([1,4,6,7]))



        lat => lattice('triangle',2,2,1,.true.,.true.,.true.)
        call setup_exchange_matrix(lat)

        t_tJ_model = .true.
        t_heisenberg_model = .false.

        call assert_equals(h_cast(-2.0),get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-2.0), get_diag_helement_heisenberg([1,4,6,7]))

        call assert_equals(h_cast(0.0), get_diag_helement_heisenberg([1,3,5,7]))

        t_tJ_model = .false. 
        t_heisenberg_model = .true. 

        call assert_equals(h_cast(6.0/4.0), get_diag_helement_heisenberg([1,3,5,7]))

        call assert_equals(h_cast(-0.5),get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-0.5), get_diag_helement_heisenberg([1,4,6,7]))


        nel = -1
        nbasis = -1
        NIfTot = -1

    end subroutine get_diag_helement_heisenberg_test

    subroutine get_offdiag_helement_heisenberg_test
        use SystemData, only: nel, nbasis 
        use real_space_hubbard, only: lat 
        use lattice_mod, only: lattice

        integer :: ex(2,2) 

        nel = 2 
        nbasis = 4 
        lat => lattice('chain', 2, 1, 1, .true., .true., .true.) 
        call setup_exchange_matrix(lat) 

        print *, "" 
        print *, "testing: get_offdiag_helement_heisenberg"
        ex(1,:) = [1,2]
        ex(2,:) = [3,4]

        call assert_equals(h_cast(0.0), get_offdiag_helement_heisenberg([1,2],


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
        use SystemData, only: nbasis
        use lattice_mod, only: lattice 
        use real_space_hubbard, only: lat 

        print *, "" 
        print *, "testing: get_umat_el_heisenberg"

        lat => lattice('chain', 4, 1, 1, .true., .true., .true.) 
        exchange_j = 1.0
        nbasis = 8
        call setup_exchange_matrix(lat)

        call assert_equals(h_cast(0.0), get_umat_el_heisenberg(1,1,1,1))
        call assert_equals(h_cast(0.0), get_umat_el_heisenberg(1,2,3,4))
        call assert_equals(h_cast(0.0), get_umat_el_heisenberg(1,1,2,2))
        call assert_equals(h_cast(0.0), get_umat_el_heisenberg(1,3,1,3))
        call assert_equals(h_cast(0.0), get_umat_el_heisenberg(1,3,3,1))
        call assert_equals(h_cast(0.0), get_umat_el_heisenberg(3,1,1,3))

        ! how do i really encode the exchange in the umat in terms of spatial 
        ! orbitals?! i am not sure i do it right 
        ! since i access the umat only depending on the electrons k,l usually 
        ! it should be fine if i always return the exhange j when the 
        ! orbitals fit 
        call assert_equals(h_cast(1.0), get_umat_el_heisenberg(1,2,1,2))
        call assert_equals(h_cast(1.0), get_umat_el_heisenberg(1,2,2,1))
        call assert_equals(h_cast(1.0), get_umat_el_heisenberg(2,1,1,2))
        call assert_equals(h_cast(1.0), get_umat_el_heisenberg(2,1,2,1))

        call assert_equals(h_cast(1.0), get_umat_el_heisenberg(1,4,1,4))
        call assert_equals(h_cast(1.0), get_umat_el_heisenberg(2,3,3,2))

        nbasis = -1

    end subroutine get_umat_el_heisenberg_test

end program test_tJ_model
