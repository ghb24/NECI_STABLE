#include "macros.h"

! i guess for a new module it would be nice to have one module which tests 
! all the functions and programs of a module .. check that out.. 

! all the global variables in NECI make it a bit tough to do unit-tests

! no.. i guess, atleast in my framework it would be better to write one 
! unit test module which tests all the functionality of a other module 
program test_real_space_hubbard 
    
    use real_space_hubbard
    use constants 
    use fruit 
    use SystemData, only: lattice_type, t_new_real_space_hubbard, t_trans_corr, & 
                          trans_corr_param, t_lattice_model
    use lattice_mod, only: lat

    implicit none 

    integer :: failed_count 

    t_new_real_space_hubbard = .true.
    t_lattice_model = .true.

    call init_fruit()
    ! run the test-driver 
    call real_space_hubbard_test_driver()
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine real_space_hubbard_test_driver() 
        ! this is the main function which calls all the other tests 
        
        ! or try running it with the provided runner of fruit: 
        call run_test_case(get_umat_el_hub_test, "get_umat_el_hub_test")
        call run_test_case(init_tmat_test, "init_tmat_test")
        call run_test_case(get_helement_test, "get_helement_test")
        call run_test_case(gen_excit_rs_hubbard_test, "gen_excit_rs_hubbard_test")
        call run_test_case(init_real_space_hubbard_test, "init_real_space_hubbard_test")
        call run_test_case(trans_corr_fac_test, "trans_corr_fac_test")
        call run_test_case(create_cum_list_rs_hubbard_test, "create_cum_list_rs_hubbard_test")
        call run_test_case(create_avail_neighbors_list_test, "create_avail_neighbors_list_test")
        call run_test_case(calc_pgen_rs_hubbard_test, "calc_pgen_rs_hubbard_test")
        call run_test_case(create_neel_state_chain_test, "create_neel_state_chain_test")
        call run_test_case(create_neel_state_test, "create_neel_state_test")
        call run_test_case(determine_optimal_time_step_test, "determine_optimal_time_step_test")
        call run_test_case(get_optimal_correlation_factor_test, "get_optimal_correlation_factor_test")
        call run_test_case(calc_n_double_test, "calc_n_double_test")
        call run_test_case(right_most_zero_test, "right_most_zero_test")
        call run_test_case(create_one_spin_basis_test, "create_one_spin_basis_test")
        call run_test_case(set_alpha_beta_spins_test, "set_beta_spins_test")
        call run_test_case(combine_spin_basis_test, "combine_spin_basis_test")
        call run_test_case(create_all_open_shell_dets_test, "create_all_open_shell_dets_test")

    end subroutine real_space_hubbard_test_driver

    subroutine get_optimal_correlation_factor_test
        use SystemData, only: uhub, bhub 
        use lattice_mod, only: lattice

        uhub = 1.0
        bhub = 1.0 

        lat => lattice('chain', 2, 1, 1, .true.,.true.,.true.)

        print *, "" 
        print *, "testing: get_optimal_correlation_factor: "

        call assert_equals(-log(1.25_dp), get_optimal_correlation_factor())
        uhub = 4.0
        call assert_equals(-log(2.0_dp), get_optimal_correlation_factor())

        lat => lattice('square', 2, 2, 1, .true.,.true.,.true.)
        call assert_equals(-log(1.5_dp), get_optimal_correlation_factor())
        uhub = 8.0
        call assert_equals(-log(2.0_dp), get_optimal_correlation_factor())

        uhub = 0.0
        bhub = 0.0 


    end subroutine get_optimal_correlation_factor_test

    subroutine create_neel_state_test
        use SystemData, only: nel, lattice_type, nbasis, length_x, length_y
        use lattice_mod, only: lattice


        print *, "" 
        print *, "testing: create_neel_state" 
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)
        length_x = 4
        length_y = 1

        lattice_type = 'chain'
        nel = 4
        nbasis = 8
        call assert_equals([1,4,5,8], create_neel_state(), 4)

        nel = 3
        call assert_equals([1,4,5], create_neel_state(), 3)

        lat => lattice('square', 2, 2, 1, .true., .true., .true.)
        length_x = 2
        length_y = 2
        lattice_type = 'square'
        nel = 4
        call assert_equals([1,4,6,7], create_neel_state(), 4) 

        nel = 3
        call assert_equals([1,4,6], create_neel_state(), 3) 

        lat => lattice('rectangle', 3, 2, 1, .true.,.true.,.true.)
        length_x = 3

        lattice_type = 'rectangle'
        nel = 6
        nbasis = 12
        call assert_equals([1,4,5,8,9,12], create_neel_state(), 6)

        nel = 4
        call assert_equals([1,4,5,8], create_neel_state(), 4)

        nel = 3
        call assert_equals([1,4,5], create_neel_state(), 3)

        lat => lattice('rectangle', 2,4,1,.true.,.true.,.true.)
        length_x = 2
        length_y = 4
        nel = 8 
        nbasis = 16
        call assert_equals([1, 4, 6, 7, 9, 12, 14, 15], create_neel_state(), 8)

        nel = 7 
        call assert_equals([1, 4, 6, 7, 9, 12, 14], create_neel_state(), 7)

        nel = 6
        call assert_equals([1, 4, 6, 7, 9, 12], create_neel_state(), 6)

        lat => lattice('rectangle', 3, 4, 1, .true., .true., .true.)
        length_x = 3
        nel = 12 
        nbasis = 24
        call assert_equals([1, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21, 24], & 
            create_neel_state(), 12)

        lattice_type = 'tilted'
        lat => lattice(lattice_type, 2, 2, 1, .true., .true., .true.)
        length_x = 2
        length_y = 2

        nbasis = 16
        nel = 8 
        call assert_equals([1, 3, 6, 7, 10, 11, 14, 16], create_neel_state(),8)

        nel = 7 
        call assert_equals([1, 3, 6, 7, 10, 11, 14], create_neel_state(),7)

        nel = 6 
        call assert_equals([1, 3, 6, 7, 10, 11], create_neel_state(),6)

        nel = 5
        call assert_equals([1, 3, 6, 7, 10], create_neel_state(),5)

        lat => lattice(lattice_type, 3, 3, 1, .true., .true., .true.)
        length_x = 3
        length_y = 3

        nbasis  = 36
        nel = 18 
        call assert_equals([1,  3, 6, 7,  9, 12, 13, 16, 17,  20, 21, 24, 25, 28, & 
            30, 31, 34, 36], create_neel_state(), 18)

        nel = 17 
        call assert_equals([1,  3, 6, 7,  9, 12, 13, 16, 17,  20, 21, 24, 25, 28, & 
            30, 31, 34], create_neel_state(), 17)

        nel = 16
        call assert_equals([1,  3, 6, 7,  9, 12, 13, 16, 17,  20, 21, 24, 25, 28, & 
            30, 31], create_neel_state(), 16)

        nel = 10
        call assert_equals([1,  3, 6, 7,  9, 12, 13, 16, 17,  20], create_neel_state(), 10)

        length_x = -1 
        length_y = -1 
        nel = -1 
        nbasis = -1

    end subroutine create_neel_state_test

    subroutine create_neel_state_chain_test
        use SystemData, only: nel

        print *, "" 
        print *, "testing: create_neel_state_chain"

        nel = 1
        call assert_equals([1], create_neel_state_chain(), 1)
        nel = 2 
        call assert_equals([1,4], create_neel_state_chain(), 2)
        nel = 3 
        call assert_equals([1,4,5], create_neel_state_chain(), 3)
        nel = 4
        call assert_equals([1,4,5,8], create_neel_state_chain(), 4)

        nel = -1

    end subroutine create_neel_state_chain_test

    subroutine create_cum_list_rs_hubbard_test
        use SystemData, only: nel 
        use bit_rep_data, only: niftot
        use Detbitops, only: encodebitdet
        use constants, only: n_int, dp
        use OneEInts, only: tmat2d

        integer(n_int), allocatable :: ilut(:)
        real(dp) :: cum_sum
        real(dp), allocatable :: cum_arr(:)
        integer, allocatable :: neighbors(:)

        trans_corr_param = 0.0
        nel = 2
        niftot = 0

        allocate(ilut(0:niftot))

        print *, ""
        print *, "testing: create_cum_list_rs_hubbard" 
        call encodebitdet([1,2], ilut)
!         neighbors = [3,5,7,9]
        allocate(tmat2d(10,10))
        tmat2d = 1.0

        call create_cum_list_rs_hubbard(ilut, 1, [3,5,7,9], cum_arr, cum_sum)
        call assert_equals([1.0,2.0,3.0,4.0], real(cum_arr), 4)
        call assert_equals(4.0_dp, cum_sum)
        
        call create_cum_list_rs_hubbard(ilut, 2, [4,6,8,10], cum_arr, cum_sum)
        call assert_equals([1.0,2.0,3.0,4.0], real(cum_arr), 4)
        call assert_equals(4.0_dp, cum_sum)

        call create_cum_list_rs_hubbard(ilut, 2, [4,6,8], cum_arr, cum_sum)
        call assert_equals([1.0,2.0,3.0], real(cum_arr), 3)
        call assert_equals(3.0_dp, cum_sum)

        call create_cum_list_rs_hubbard(ilut, 1, [1,3,5], cum_arr, cum_sum)
        call assert_equals([0.0,1.0,2.0], real(cum_arr), 3)
        call assert_equals(2.0_dp, cum_sum)

        call create_cum_list_rs_hubbard(ilut, 2, [4,2,6], cum_arr, cum_sum)
        call assert_equals([1.0,1.0,2.0], real(cum_arr), 3)
        call assert_equals(2.0_dp, cum_sum)

        call encodebitdet([1,3], ilut)
        call create_cum_list_rs_hubbard(ilut, 1, [9,3,5,1], cum_arr, cum_sum)
        call assert_equals([1.0,1.0,2.0,2.0], real(cum_arr), 4)
        call assert_equals(2.0_dp, cum_sum)

        print *, ""
        print *, "and now with a transcorrelated hamiltonian with K = 1.0"
        trans_corr_param = 1.0 
        t_trans_corr = .true.
        call encodebitdet([1,2], ilut)

        call create_cum_list_rs_hubbard(ilut, 1, [3,5,7,9], cum_arr, cum_sum)
        call assert_equals([exp(1.0),2*exp(1.0),3*exp(1.0),4*exp(1.0)], real(cum_arr), 4)
        call assert_equals(4*exp(1.0_dp), cum_sum)

        call create_cum_list_rs_hubbard(ilut, 2, [4,6,8,10], cum_arr, cum_sum)
        call assert_equals([exp(1.0),2*exp(1.0),3*exp(1.0),4*exp(1.0)], real(cum_arr), 4)
        call assert_equals(4*exp(1.0_dp), cum_sum)

        call encodebitdet([1,3], ilut)
        call create_cum_list_rs_hubbard(ilut, 1, [9,3,5,1], cum_arr, cum_sum)
        call assert_equals([1.0,1.0,2.0,2.0], real(cum_arr), 4)
        call assert_equals(2.0_dp, cum_sum)

        nel = 4 
        call encodebitdet([1,2,3,6],ilut)
        call create_cum_list_rs_hubbard(ilut, 1, [3,5,7,9], cum_arr, cum_sum)
        call assert_equals([0.0,1.0,1.0+exp(1.0),1.0+2*exp(1.0)], real(cum_arr), 4)
        call assert_equals(1.0+2*exp(1.0_dp), cum_sum)

        call create_cum_list_rs_hubbard(ilut, 3, [1,5,7,9], cum_arr, cum_sum)
        call assert_equals([0.0,exp(-1.0),exp(-1.0)+1.0,exp(-1.0)+2.0],&
            real(cum_arr), 4)
        call assert_equals(exp(-1.0_dp)+2.0_dp, cum_sum)

        nel = 3 
        call encodebitdet([1,2,3], ilut)
        call create_cum_list_rs_hubbard(ilut, 2, [4,8], cum_arr, cum_sum) 
        call assert_equals([1.0, 1.0+exp(1.0)], real(cum_arr), 2)
        call assert_equals(1.0+exp(1.0_dp), cum_sum)

        t_trans_corr = .false.
        deallocate(tmat2d)
        nel = -1 
        niftot = -1
        trans_corr_param = 0.0

    end subroutine create_cum_list_rs_hubbard_test

    subroutine create_avail_neighbors_list_test
        use SystemData, only: nel
        use constants, only: n_int 
        use Detbitops, only: encodebitdet 
        use bit_rep_data, only: niftot

        integer(n_int), allocatable :: ilut(:)
        integer, allocatable :: orbs(:)
        integer :: n_orbs

        nel = 2 
        niftot = 0

        allocate(ilut(0:niftot))

        print *, "" 
        print *, "testing: create_avail_neighbors_list"

        call encodebitdet([1,2], ilut)

        call create_avail_neighbors_list(ilut, [3,5,7,9], orbs, n_orbs)
        call assert_equals(4, n_orbs)
        call assert_equals([3,5,7,9], orbs, 4)

        call create_avail_neighbors_list(ilut, [2,4,6,8], orbs, n_orbs)
        call assert_equals(3, n_orbs)
        call assert_equals([4,6,8], orbs, 3)

        call create_avail_neighbors_list(ilut, [1,3,5,7,9], orbs, n_orbs)
        call assert_equals(4, n_orbs)
        call assert_equals([3,5,7,9], orbs, 4)

        call create_avail_neighbors_list(ilut, [1,2,3,5,7,9], orbs, n_orbs)
        call assert_equals(4, n_orbs)
        call assert_equals([3,5,7,9], orbs, 4)

        call create_avail_neighbors_list(ilut, [1,2], orbs, n_orbs)
        call assert_equals(0, n_orbs)

        nel = -1 
        niftot = -1

    end subroutine create_avail_neighbors_list_test

    subroutine calc_pgen_rs_hubbard_test
        use SystemData, only: nel, nbasis, bhub
        use constants, only: n_int, dp
        use lattice_mod, only: lattice
        use bit_rep_data, only: niftot 
        use Detbitops, only: encodebitdet 
        use OneEInts, only: tmat2d

        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2)
        integer, allocatable :: nI(:)

        ex = 0

        nel = 2
        niftot = 0
        bhub = 1.0
        allocate(ilut(0:niftot))

        print *, "" 
        print *, "testing: calc_pgen_rs_hubbard" 
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)
        nbasis = 8
        call init_get_helement_hubbard()

        trans_corr_param = 0.0 
        t_trans_corr = .false. 
        allocate(ni(nel))
        nI = [1,2]
        call encodebitdet(nI, ilut) 
        ex(1,1) = 1 
        ex(2,1) = 3

        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ni, ilut, ex, 2))
        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ni, ilut, ex, 0))
        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ni, ilut, ex, -1))
        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ni, ilut, ex, 3))

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ni,ilut,ex,1))

        ex(1,1) = 2
        ex(2,1) = 8

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ni,ilut,ex,1))

        t_trans_corr = .true. 

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ni,ilut,ex,1))

        ex(1,1) = 1 
        ex(2,1) = 3

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ni,ilut,ex,1))

        ex(1,1) = 2
        ex(2,1) = 8
        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ni,ilut,ex,1))

        trans_corr_param = 1.0 
        call assert_equals(1.0/(4.0_dp), calc_pgen_rs_hubbard(ni,ilut,ex,1))

        nel = 3 
        nI = [1,2,3]
        call encodebitdet(nI, ilut) 
        ex(1,1) = 2
        ex(2,1) = 4
        ! think about the 
!         call assert_equals(1.0
        ! some rounding errors below, otherwise correct
        call assert_equals(real(1.0/3.0_dp*1.0/(1.0+exp(1.0_dp))), real(calc_pgen_rs_hubbard(ni, ilut, ex,1)),1e-12)

        ex(2,1) = 8 
        call assert_equals(exp(1.0_dp)/(3.0*(1.0+exp(1.0_dp))), calc_pgen_rs_hubbard(ni, ilut, ex,1))

        ex(1,1) = 1
        ex(2,1) = 7
        call assert_equals(1.0/3.0_dp, calc_pgen_rs_hubbard(ni, ilut, ex,1))

        nel = -1
        niftot = -1
        nbasis = -1
        trans_corr_param = 0.0
        deallocate(tmat2d)

    end subroutine calc_pgen_rs_hubbard_test

    subroutine trans_corr_fac_test
        use Detbitops, only: encodebitdet
        use constants, only: n_int
        use SystemData, only: nel
        use bit_rep_data, only: niftot

        integer(n_int), allocatable :: ilut(:)

        print *, "" 
        print *, "Testing transcorrelation factor" 

        niftot = 0 
        trans_corr_param = 1.0_dp

        allocate(ilut(0:niftot))

        ilut = 0

        nel = 2
        call encodebitdet([1,2],ilut)

        call assert_equals(exp(1.0_dp), trans_corr_fac(ilut,1,3))
        call assert_equals(exp(1.0_dp), trans_corr_fac(ilut,2,4))
        call assert_equals(exp(-1.0_dp), trans_corr_fac(ilut,3,1))
        call assert_equals(exp(-1.0_dp), trans_corr_fac(ilut,4,2))

        ilut = 0
        call encodebitdet([1,4],ilut)

        call assert_equals(exp(-1.0_dp), trans_corr_fac(ilut,1,3))
        call assert_equals(exp(-1.0_dp), trans_corr_fac(ilut,4,2))
        call assert_equals(1.0_dp, trans_corr_fac(ilut,1,5))
        call assert_equals(1.0_dp, trans_corr_fac(ilut,4,6))

        ilut = 0
        nel = 4
        call encodebitdet([1,2,3,6], ilut)
        call assert_equals(1.0_dp, trans_corr_fac(ilut,2,4))
        call assert_equals(1.0_dp, trans_corr_fac(ilut,1,5))

        niftot = -1
        trans_corr_param = 0.0

    end subroutine trans_corr_fac_test

    subroutine init_real_space_hubbard_test
        use SystemData, only: lattice_type, length_x, length_y, nbasis, nel, &
                              bhub, uhub, ecore, t_new_real_space_hubbard
        use OneEInts, only: tmat2d
        use fcimcdata, only: pSingles, pDoubles, tsearchtau, tsearchtauoption
        use CalcData, only: tau
        use tau_search, only: max_death_cpt
        use procedure_pointers, only: get_umat_el, generate_excitation
        use lattice_mod, only: lattice_deconstructor

        print *, "" 
        print *, "testing: init_real_space_hubbard"

        print *, "for a 100-site periodic chain: "
        lattice_type = 'chain' 
        length_x = 100
        length_y = 1
        nel = 2 
        bhub = 1

        call init_real_space_hubbard() 

        ! test the lattice strucure 

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
        call assert_equals(1.0_dp, pSingles)
        call assert_equals(0.0_dp, pDoubles)
        call assert_true(.not. tsearchtau)
        call assert_true(tsearchtauoption)
        call assert_equals(0.0_dp, max_death_cpt)
        call assert_true(associated(get_umat_el))
        call assert_true(associated(generate_excitation))
        call assert_equals(0.25 * lat_tau_factor, tau)

        call lattice_deconstructor(lat)

        lattice_type = '' 
        length_x = -1 
        length_y = -1 
        nel = -1 
        bhub = 0
        nbasis = -1
        tau = 0.0_dp
        deallocate(tmat2d)
        nullify(get_umat_el)
        nullify(generate_excitation)

        print *, "" 
        print *, "testing: 3x3 square lattice with 2 electrons" 
        lattice_type = 'square'
        length_x = 3
        length_y = 3 
        nel = 2
        bhub = 1 

        call init_real_space_hubbard() 
        call assert_equals(2, lat%get_ndim()) 
        call assert_equals(9, lat%get_nsites())
        call assert_true(lat%is_periodic())
        call assert_equals(4, lat%get_nconnect_max())

        call assert_equals(18, nbasis)
        call assert_true(associated(tmat2d))
        call assert_true(associated(g1))
        call assert_equals(0.0_dp, ecore)
        call assert_equals(1.0_dp, pSingles)
        call assert_equals(0.0_dp, pDoubles)
        call assert_true(.not. tsearchtau)
        call assert_true(tsearchtauoption)
        call assert_equals(0.0_dp, max_death_cpt)
        call assert_true(associated(get_umat_el))
        call assert_true(associated(generate_excitation))
        call assert_equals(1.0/8.0_dp * lat_tau_factor, tau)


        call lattice_deconstructor(lat)

        lattice_type = '' 
        length_x = -1 
        length_y = -1 
        nel = -1 
        bhub = 0
        nbasis = -1
        tau = 0.0
        deallocate(tmat2d)
        nullify(get_umat_el)
        nullify(generate_excitation)

        print *, "" 
        print *, "testing: 3x3 square lattice with 2 electrons" 
        lattice_type = 'triangle'
        length_x = 3
        length_y = 3 
        nel = 2
        bhub = 1 

        call init_real_space_hubbard() 
        call assert_equals(2, lat%get_ndim()) 
        call assert_equals(9, lat%get_nsites())
        call assert_true(lat%is_periodic())
        call assert_equals(6, lat%get_nconnect_max())

        call assert_equals(18, nbasis)
        call assert_true(associated(tmat2d))
        call assert_true(associated(g1))
        call assert_equals(0.0_dp, ecore)
        call assert_equals(1.0_dp, pSingles)
        call assert_equals(0.0_dp, pDoubles)
        call assert_true(.not. tsearchtau)
        call assert_true(tsearchtauoption)
        call assert_equals(0.0_dp, max_death_cpt)
        call assert_true(associated(get_umat_el))
        call assert_true(associated(generate_excitation))
        call assert_equals(1.0/12.0_dp * lat_tau_factor, tau)

        call lattice_deconstructor(lat)

        lattice_type = '' 
        length_x = -1 
        length_y = -1 
        nel = -1 
        bhub = 0
        nbasis = -1
        deallocate(tmat2d)
        nullify(get_umat_el)
        nullify(generate_excitation)


    end subroutine init_real_space_hubbard_test

    subroutine gen_excit_rs_hubbard_test
        use dsfmt_interface, only: dsfmt_init 
        use lattice_mod, only: lattice, lattice_deconstructor
        use SystemData, only: nel, nbasis
        use bit_rep_data, only: niftot
        use Detbitops, only: encodebitdet
        use constants, only: n_int, dp
        use fcimcdata, only: excit_gen_store_type
        
        integer, allocatable :: nI(:), nJ(:)
        integer(n_int), allocatable :: ilutI(:), ilutJ(:)
        integer :: ex(2,2), ic 
        logical :: tpar 
        real(dp) :: pgen 
        HElement_t(dp) :: hel
        type(excit_gen_store_type) :: store
        logical :: found_all, t_found(6)
        
        print *, "" 
        print *, "testing gen_excit_rs_hubbard"
        lat => lattice('chain', 4, 1, 1, .false., .false., .false.)
        nel = 2 
        call dsfmt_init(1)
        allocate(nI(nel))
        allocate(nJ(nel))

        nI = [1,2]
        niftot = 0 
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet(nI, ilutI)
        
        found_all = .false.
        t_found = .false.

        nbasis = 8
        call init_get_helement_hubbard()

        ! try to get all the excitations here. 
        do while(.not. found_all)
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_found(1)) then 
                t_found(1) = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(3, ex(2,1))
                call assert_equals(6, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.5_dp, pgen)
            end if
            if (all(nJ == [1,4]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(4, ex(2,1))
                call assert_equals(9, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.5_dp, pgen)
            end if
            found_all = all(t_found(1:2))
        end do

        nI = [3,4]
        call encodebitdet(nI, ilutI)

        found_all = .false. 
        t_found = .false. 

        ! try to get all the excitations here. 
        do while(.not. found_all)
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_found(1)) then 
                t_found(1) = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(4, ex(1,1))
                call assert_equals(2, ex(2,1))
                call assert_equals(6, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [1,4]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(3, ex(1,1))
                call assert_equals(1, ex(2,1))
                call assert_equals(9, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            if (all(nJ == [3,6]) .and. .not. t_found(3)) then 
                t_found(3) = .true.
                call assert_equals([3,6], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(4, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(36, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [4,5]) .and. .not. t_found(4)) then 
                t_found(4) = .true. 
                call assert_equals([4,5], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(3, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(24, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            found_all = all(t_found(1:4))
        end do

        ! and also try it on a periodic chain 
        lat => lattice('chain', 3, 1, 1, .true., .true., .true.)
        nbasis = 6
        call init_tmat(lat)

        nI = [1,2]
        call encodebitdet(nI, ilutI) 

        found_all = .false. 
        t_found = .false. 

        ! try to get all the excitations here. 
        do while(.not. found_all)
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_found(1)) then 
                t_found(1) = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(3, ex(2,1))
                call assert_equals(6, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [1,4]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(4, ex(2,1))
                call assert_equals(9, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            if (all(nJ == [1,6]) .and. .not. t_found(3)) then 
                t_found(3) = .true.
                call assert_equals([1,6], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(33, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [2,5]) .and. .not. t_found(4)) then 
                t_found(4) = .true. 
                call assert_equals([2,5], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(18, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            found_all = all(t_found(1:4))
        end do

        nel = -1 
        call lattice_deconstructor(lat)

        print *, "" 
        print *, "for a periodic 2x2 square lattice with 2 electrons" 
        lat => lattice('square', 2, 2, 1, .true., .true., .true.) 
        nel = 2 

        nbasis = 8
        call init_tmat(lat)

        nI = [1,2] 
        found_all = .false. 
        t_found = .false. 

        do while(.not. found_all) 
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [1,4]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(4, ex(2,1))
                call assert_equals(9, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [2,3]) .and. .not. t_found(2)) then
                t_found(2) = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(3, ex(2,1))
                call assert_equals(6, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [1,6]) .and. .not. t_found(3)) then 
                t_found(3) = .true. 
                call assert_equals([1,6], nJ, 2)
                call assert_equals(1, ic)
                call assert_equals(2, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(33, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [2,5]) .and. .not. t_found(4)) then
                t_found(4) = .true.
                call assert_equals([2,5], nJ, 2)
                call assert_equals(1, ic)
                call assert_equals(1, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(18, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            found_all = all(t_found(1:4))
        end do

        print *, "" 
        print *, "for a periodic 2x2 triangular lattice: "
        t_found = .false. 
        found_all = .false.
        lat => lattice('triangle', 2,2,1,.true.,.true.,.true.)

        call init_tmat(lat)

        do while(.not. found_all) 
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [1,4]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(4, ex(2,1))
                call assert_equals(9, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if
            if (all(nJ == [2,3]) .and. .not. t_found(2)) then
                t_found(2) = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(3, ex(2,1))
                call assert_equals(6, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if
            if (all(nJ == [1,6]) .and. .not. t_found(3)) then 
                t_found(3) = .true. 
                call assert_equals([1,6], nJ, 2)
                call assert_equals(1, ic)
                call assert_equals(2, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(33, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if
            if (all(nJ == [2,5]) .and. .not. t_found(4)) then
                t_found(4) = .true.
                call assert_equals([2,5], nJ, 2)
                call assert_equals(1, ic)
                call assert_equals(1, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(18, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if
            if (all(nJ == [1,8]) .and. .not. t_found(5)) then 
                t_found(5) = .true. 
                call assert_equals([1,8], nJ, 2)
                call assert_equals(1, ic)
                call assert_equals(2, ex(1,1))
                call assert_equals(8, ex(2,1))
                call assert_equals(129, int(ilutJ(0)))
                call assert_true(.not.tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if
            if (all(nJ == [2,7]) .and. .not. t_found(6)) then 
                t_found(6) = .true. 
                call assert_equals([2,7], nJ, 2)
                call assert_equals(1, ic)
                call assert_equals(1, ex(1,1))
                call assert_equals(7, ex(2,1))
                call assert_equals(66, int(ilutJ(0)))
                call assert_true(tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if

            found_all = all(t_found)
        end do

        nel = -1 
        nbasis = -1
        call lattice_deconstructor(lat)

        print *, "" 
        print *, "and now with the transcorrelated Hamiltonian with K = 1"
        ! i definetly have to do the test for the transcorrelated hamiltonian 
        ! now.. 


    end subroutine gen_excit_rs_hubbard_test

    subroutine get_helement_test
        use SystemData, only: nel, tCSF, bhub, uhub, nbasis, G1
        use bit_rep_data, only: nifd, niftot
        use Determinants, only: get_helement
        use constants, only: dp, n_int
        use OneEInts, only: tmat2d
        use procedure_pointers, only: get_umat_el
        use lattice_mod, only: lattice
        use detbitops, only: encodebitdet
        
        ! do i have to set all the flags.. 
        integer(n_int), allocatable :: ilutI(:), ilutJ(:)

        print *, ""
        print *, "testing get_helement for the real-space hubbard model"

        nbasis = 4 
        get_umat_el => get_umat_el_hub

        nel = 2
        bhub = 1.0

        allocate(G1(nbasis))
        G1([1,3])%ms = -1 
        G1([2,4])%ms = 1

        print *, "here?"
        lat => lattice('chain', 2, 1, 1, .false., .false., .false.)
        call init_tmat(lat)
        call init_get_helement_hubbard()

        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2],0))
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2]))

        uhub = 1.0
        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,2],0))
        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,2]))

        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,4],1))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,2],[2,3],1))

        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,4]))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,2],[2,3]))

        nbasis = 6 
        deallocate(g1)
        allocate(g1(nbasis))
        g1([1,3,5])%ms = -1 
        g1([2,4,6])%ms = 1 

        lat => lattice('chain', 3, 1, 1, .true., .true., .true.)
        call init_tmat(lat)

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,3],0))
        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,3]))

        call assert_equals(h_cast(1.0_dp), get_helement([1,3],[1,5],1))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,3],[3,5],1))

        call assert_equals(h_cast(1.0_dp), get_helement([1,3],[1,5]))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,3],[3,5]))

        niftot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet([1,2],ilutI)
        call encodebitdet([1,2],ilutJ)

        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,2],ilutI,ilutJ))

        call encodebitdet([2,4],ilutJ)
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[2,4],ilutI,ilutJ))

        call encodebitdet([2,3],ilutJ)
        call assert_equals(h_cast(-1.0_dp), get_helement([1,2],[2,3],ilutI,ilutJ))

        call encodebitdet([1,4],ilutJ)
        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,4],ilutJ,ilutJ))

        nel = 4 
        call assert_equals(h_cast(2.0_dp), get_helement([1,2,3,4],[1,2,3,4]))
        call assert_equals(h_cast(2.0_dp), get_helement([1,2,3,4],[1,2,3,4],0))

        call assert_equals(h_cast(1.0_dp), get_helement([1,2,3,4],[1,2,3,6]))
        call assert_equals(h_cast(1.0_dp), get_helement([1,2,3,4],[1,2,3,6],1))

        call assert_equals(h_cast(1.0_dp), get_helement([1,2,3,4],[1,3,4,6],1))
        call assert_equals(h_cast(1.0_dp), get_helement([1,2,3,4],[1,3,4,6]))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,2,3,4],[2,3,4,5],1))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,2,3,4],[2,3,4,5]))

        ! for a 2x2 square lattice 
        lat => lattice('square', 2,2,1,.true.,.true.,.true.)
        nel = 2
        nbasis = 8 
        call init_tmat(lat) 

        deallocate(g1)
        allocate(g1(nbasis))
        g1(1:7:2)%ms = -1 
        g1(2:8:2)%ms = 1 

        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,2]))
        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,4]))
        call assert_equals(h_cast(1.0_dp), get_helement([1,2],[1,6]))
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,8]))
        call assert_equals(h_cast(-1.0_dp), get_helement([1,2],[2,3]))

        print *, "" 
        print *, "and now for transcorrelated hamiltonian with K = 1" 
        t_trans_corr = .true.
        trans_corr_param = 1.0
        call assert_equals(h_cast(1.0 * exp(1.0_dp)), get_helement([1,2],[1,4]))
        call assert_equals(h_cast(1.0 * exp(-1.0_dp)), get_helement([1,4],[1,2]))
        call assert_equals(h_cast(-1.0 * exp(1.0_dp)), get_helement([1,2],[2,3]))
        call assert_equals(h_cast(-1.0 * exp(-1.0_dp)), get_helement([2,3],[1,2]))
        call assert_equals(h_cast(0.0_dp), get_helement([2,3],[2,5]))
        call assert_equals(h_cast(1.0_dp), get_helement([2,3],[2,7]))
        call assert_equals(h_cast(0.0_dp), get_helement([2,3],[3,8]))
        call assert_equals(h_cast(-1.0_dp), get_helement([2,3],[3,6]))

        print *, ""

        nel = -1 
        nbasis = -1 
        nullify(get_umat_el)
        deallocate(g1) 
        deallocate(tmat2d)
        t_trans_corr = .false.

    end subroutine get_helement_test

    subroutine init_tmat_test
        use lattice_mod, only: lattice
        use SystemData, only: nbasis, bhub
        use OneEInts, only: tmat2d


        print *, "" 
        print *, "testing: init_tmat" 

        print *, "for a 4 site, non-periodic chain geometry"
        nbasis = 8 
        bhub = 1.0
        lat => lattice('chain', 4, 1, 1, .false., .false., .false.)

        call init_tmat(lat)

        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), &
            h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(1,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0), &
            h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(2,:))), 8)
        call assert_equals([h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), &
            h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(3,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), & 
            h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(4,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), &
            h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(7,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), &
            h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(8,:))), 8)

        print *, ""
        print *, "for a 4 site periodic chain geometry: "
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)

        call init_tmat(lat)

        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), &
            h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0)], h_cast(real(tmat2d(1,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0),&
            h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0)], h_cast(real(tmat2d(2,:))), 8)
        call assert_equals([h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), &
            h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(3,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), &
            h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(4,:))), 8)
        call assert_equals([h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), &
            h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(7,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), &
            h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(8,:))), 8)

        ! todo: more tests for other lattices later!

        print *, "" 
        print *, "for a 2x2 periodic square lattice" 
        lat => lattice('square', 2,2, 1, .true., .true., .true.) 

        call init_tmat(lat) 

        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(1,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(2,:))), 8)
        call assert_equals([h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0)], h_cast(real(tmat2d(3,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0)], h_cast(real(tmat2d(4,:))), 8)
        call assert_equals([h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0)], h_cast(real(tmat2d(5,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0)], h_cast(real(tmat2d(6,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(7,:))), 8)
        call assert_equals([h_cast(0.0), h_cast(0.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(1.0), h_cast(0.0), h_cast(0.0)], h_cast(real(tmat2d(8,:))), 8)

        nbasis = -1 
        bhub = 0.0
        deallocate(tmat2d)

    end subroutine init_tmat_test


    subroutine test_init_lattice() 
 
    end subroutine test_init_lattice

    subroutine get_umat_el_hub_test
        use SystemData, only: uhub, nbasis

        print *, ""
        print *, "testing get_umat_el_hub" 
! 
        ! for the mateles then..
!         ecore = 0.0
!         tcsf = .false. 
!         texch = .false. 
!         treltvy = .false. 
!         niftot = 0
!         ! bits_n_int .. 
!         ! n_int .. 

        uhub = 1.0
        nbasis = 8
        
        call assert_equals(h_cast(1.0_dp), get_umat_el_hub(1,1,1,1))
        call assert_equals(h_cast(1.0_dp), get_umat_el_hub(2,2,2,2))
        call assert_equals(h_cast(1.0_dp), get_umat_el_hub(3,3,3,3))
        call assert_equals(h_cast(1.0_dp), get_umat_el_hub(4,4,4,4))
        call assert_equals(h_cast(0.0_dp), get_umat_el_hub(1,2,3,4))
        call assert_equals(h_cast(0.0_dp), get_umat_el_hub(1,1,1,4))
        call assert_equals(h_cast(0.0_dp), get_umat_el_hub(2,2,3,4))
        call assert_equals(h_cast(0.0_dp), get_umat_el_hub(3,2,3,4))

        uhub = 0.0
        nbasis = -1

    end subroutine

    subroutine determine_optimal_time_step_test 
        use SystemData, only: nel, nOccAlpha, nOccBeta, uhub, bhub, & 
                              t_new_real_space_hubbard, t_tJ_model, t_heisenberg_model
        use lattice_mod, only: lattice, determine_optimal_time_step

        real(dp) :: time, time_death

        t_new_real_space_hubbard = .true.
        t_trans_corr = .false.

        print *, ""
        print *, "testing determine_optimal_time_step"
        print *, "for a 4x1 periodic chain: "
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)

        nel = 6
        nOccAlpha = 3
        nOccBeta = 3 
        uhub = 1 
        bhub = 1

        call assert_equals(1.0_dp/real(abs(bhub) * nel * 2, dp), determine_optimal_time_step())

        bhub = 2 
        call assert_equals(1.0_dp/real(abs(bhub) * nel * 2, dp), determine_optimal_time_step())

        nel = 8 
        call assert_equals(1.0_dp/real(abs(bhub) * nel * 2, dp), determine_optimal_time_step())

        bhub = -2 
        call assert_equals(1.0_dp/real(abs(bhub) * nel * 2, dp), determine_optimal_time_step())

        time = determine_optimal_time_step(time_death)
        call assert_equals(1.0_dp/real(abs(uhub * nOccAlpha),dp), time_death)

        uhub = 4
        time = determine_optimal_time_step(time_death)
        call assert_equals(1.0_dp/real(abs(uhub * nOccAlpha),dp), time_death)

        nOccAlpha = 4 
        time = determine_optimal_time_step(time_death)
        call assert_equals(1.0_dp/real(abs(uhub * nOccBeta),dp), time_death)
        
        print *, ""
        print *, "for a 3x3 periodic square: "
        lat => lattice('square', 3,3,1,.true.,.true.,.true.)

        nel = 6
        nOccAlpha = 3
        nOccBeta = 3 
        uhub = 1 
        bhub = 1

        call assert_equals(1.0_dp/real(abs(bhub) * nel * 4, dp), determine_optimal_time_step())

        bhub = 2 
        call assert_equals(1.0_dp/real(abs(bhub) * nel * 4, dp), determine_optimal_time_step())

        nel = 8 
        call assert_equals(1.0_dp/real(abs(bhub) * nel * 4, dp), determine_optimal_time_step())

        bhub = -2 
        call assert_equals(1.0_dp/real(abs(bhub) * nel * 4, dp), determine_optimal_time_step())

        time = determine_optimal_time_step(time_death)
        call assert_equals(1.0_dp/real(abs(uhub * nOccAlpha),dp), time_death)

        uhub = 4
        time = determine_optimal_time_step(time_death)
        call assert_equals(1.0_dp/real(abs(uhub * nOccAlpha),dp), time_death)

        nOccAlpha = 4 
        time = determine_optimal_time_step(time_death)
        call assert_equals(1.0_dp/real(abs(uhub * nOccBeta),dp), time_death)

        nel = -1
        nOccAlpha = -1
        nOccBeta = -1
        uhub = 0
        bhub = 0
        
    end subroutine determine_optimal_time_step_test 

    subroutine calc_n_double_test

        real(dp), allocatable :: n_double(:)

        print *, ""
        print *, "testing: calc_n_double" 

        n_double = calc_n_double(4,2,2) 

        call assert_equals(3, size(n_double))
        call assert_equals(6.0_dp, n_double(1))
        call assert_equals(24.0_dp, n_double(2))
        call assert_equals(6.0_dp, n_double(3))
        call assert_equals(36.0_dp, sum(n_double))

        n_double = calc_n_double(4,1,2) 
        call assert_equals(2, size(n_double)) 
        call assert_equals(12.0_dp, n_double(1)) 
        call assert_equals(12.0_dp, n_double(2)) 
        call assert_equals(24.0_dp, sum(n_double))

        n_double = calc_n_double(4,1,1) 
        call assert_equals(2, size(n_double)) 
        call assert_equals(12.0_dp, n_double(1))
        call assert_equals(4.0_dp, n_double(2)) 
        call assert_equals(16.0_dp, sum(n_double))

        n_double = calc_n_double(6,3,3) 
        call assert_equals(4, size(n_double))
        call assert_equals(20.0_dp, n_double(1))
        call assert_equals(180.0_dp, n_double(2))
        call assert_equals(180.0_dp, n_double(3))
        call assert_equals(20.0_dp, n_double(4))
        call assert_equals(400.0_dp, sum(n_double))

        n_double = calc_n_double(6,3,2)
        call assert_equals(3, size(n_double))
        call assert_equals(60.0_dp, n_double(1))
        call assert_equals(180.0_dp, n_double(2))
        call assert_equals(60.0_dp, n_double(3))

        n_double = calc_n_double(6,2,2) 
        call assert_equals(90.0_dp, n_double(1))
        call assert_equals(120.0_dp, n_double(2))
        call assert_equals(15.0_dp, n_double(3))

    end subroutine calc_n_double_test

    subroutine right_most_zero_test
        use bit_rep_data, only: niftot, nifd
        
        integer(n_int) :: i 

        nifd = 0
        niftot = 0

        print *, "" 
        print *, "testing: right_most_zero:" 
        i = b'1001'
        call assert_equals(2, right_most_zero(i, 4))

        i = b'0011'
        call assert_equals(3, right_most_zero(i, 4))

        i = b'1010'
        call assert_equals(3, right_most_zero(i, 4))

        i = b'1100' 
        call assert_equals(5, right_most_zero(i, 4)) 

        call assert_equals(4, right_most_zero(int(b'1100',n_int), 3))

        call assert_equals(-1, right_most_zero(0_n_int, 3))

        call assert_equals(5, right_most_zero(huge(0_n_int),4))
        
        niftot = -1
        nifd = -1
    end subroutine right_most_zero_test

    subroutine create_one_spin_basis_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: alpha(:)
        integer :: i 

        niftot = 0
        nifd = 0
        print *, ""
        print *, "testing: create_one_spin_basis"

        alpha = create_one_spin_basis(4,2) 
        
        call assert_equals(6, size(alpha)) 
        call assert_equals(b'0011',int(alpha(1)))
        call assert_equals((b'0101'),int(alpha(2)))
        call assert_equals((b'0110'),int(alpha(3)))
        call assert_equals((b'1001'),int(alpha(4)))
        call assert_equals((b'1010'),int(alpha(5)))
        call assert_equals((b'1100'),int(alpha(6)))

        alpha = create_one_spin_basis(5,2) 

        call assert_equals(10, size(alpha))
        call assert_equals((b'00011'), int(alpha(1)))
        call assert_equals((b'11000'), int(alpha(10)))

        niftot = -1
        nifd = -1

    end subroutine create_one_spin_basis_test

    subroutine set_alpha_beta_spins_test
        use bit_rep_data, only: niftot, nifd

        niftot = 0
        nifd = 0

        print *, ""
        print *, "testing: set_alpha_beta_spins "

        call assert_equals((b'101'),int(set_alpha_beta_spins(int(b'11',n_int),4,.true.)))
        call assert_equals((b'00010001'), int(set_alpha_beta_spins(int(b'0101',n_int),4,.true.)))
        call assert_equals((b'01000001'), int(set_alpha_beta_spins(int(b'1001',n_int),4,.true.)))

        call assert_equals((b'1010'),int(set_alpha_beta_spins(int(b'11',n_int),4,.false.)))
        call assert_equals((b'00100010'), int(set_alpha_beta_spins(int(b'0101',n_int),4,.false.)))
        call assert_equals((b'10000010'), int(set_alpha_beta_spins(int(b'1001',n_int),4,.false.)))

        call assert_equals((b'10100000'), int(set_alpha_beta_spins(int(b'1100',n_int),4,.false.)))
        call assert_equals((b'01010000'), int(set_alpha_beta_spins(int(b'1100',n_int),4,.true.)))

        niftot = -1
        nifd = -1
    end subroutine set_alpha_beta_spins_test

    subroutine combine_spin_basis_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: basis(:), spin_basis(:)
        integer :: i 

        niftot = 0
        nifd = 0

        print *, ""
        print *, "testing: combine_spin_basis" 

        basis = combine_spin_basis(4,2,2,6,int([3,5,6,9,10,12],n_int),.false.)

        call assert_equals(6, size(basis)) 
        call assert_equals((b'01011010'),int(basis(1)))
        call assert_equals((b'01100110'),int(basis(2)))
        call assert_equals((b'01101001'),int(basis(3)))
        call assert_equals((b'10100101'),int(basis(6)))

        spin_basis = create_one_spin_basis(6,3) 

        basis = combine_spin_basis(6,3,3,20,spin_basis,.false.)

        call assert_equals(20, size(basis))
        call assert_equals((b'010101101010'), int(basis(1)))
        call assert_equals((b'101010010101'), int(basis(20)))

        basis = combine_spin_basis(4,2,1,12,int([3,5,6,9,10,12],n_int),.false.)

        call assert_equals(12, size(basis)) 
        call assert_equals((b'00011010'), int(basis(1)))
        call assert_equals((b'01001010'), int(basis(2)))
        call assert_equals((b'00100110'), int(basis(3)))
        call assert_equals((b'10100001'), int(basis(11)))
        call assert_equals((b'10100100'), int(basis(12)))

        basis = combine_spin_basis(4,1,1,12,int([1,2,4,8],n_int),.false.)
        
        call assert_equals(12, size(basis))
        call assert_equals((b'00000110'), int(basis(1)))
        call assert_equals((b'00010010'), int(basis(2)))
        call assert_equals((b'01000010'), int(basis(3)))
        call assert_equals((b'10010000'), int(basis(12)))

        niftot = -1
        nifd = -1
    end subroutine combine_spin_basis_test

    subroutine create_all_open_shell_dets_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: basis(:)

        niftot = 0
        nifd = 0

        print *, "" 
        print *, "testing: create_all_open_shell_dets"

        basis = create_all_open_shell_dets(4,2,2)

        call assert_equals(6, size(basis)) 
        call assert_equals(int(b'01011010',n_int),basis(1))
        call assert_equals(int(b'01100110',n_int),basis(2))
        call assert_equals(int(b'01101001',n_int),basis(3))
        call assert_equals(int(b'10100101',n_int),basis(6))

        niftot = -1
        nifd = -1

    end subroutine create_all_open_shell_dets_test

end program test_real_space_hubbard

