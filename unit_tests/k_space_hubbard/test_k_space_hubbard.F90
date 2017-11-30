#include "macros.h"

! i guess for a new module it would be nice to have one module which tests 
! all the functions and programs of a module .. check that out.. 

! all the global variables in NECI make it a bit tough to do unit-tests

! no.. i guess, atleast in my framework it would be better to write one 
! unit test module which tests all the functionality of a other module 
program test_k_space_hubbard 
    
    use k_space_hubbard
    use constants 
    use fruit 
    use SystemData, only: t_k_space_hubbard, t_lattice_model, nel, nbasis, & 
                          t_trans_corr, G1, nBasisMax

    implicit none 

    integer :: failed_count 


    call init_fruit()
    ! run the test-driver 
    call k_space_hubbard_test_driver()
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine init_k_space_unit_tests()
        ! since there is so much annoying outside dependency, mainly due to 
        ! momentum-conservation symmetry, lets just initialize all 
        ! necessary stuff here in front so testing is not so annoying.. 

        print *, ""
        print *, "initializing k_space unit tests"
        print *, "for simplicity do the unit tests on a 4 electron system" 
        print *, "in 4 spatial orbitals" 
        nel = 4 
        nbasis = 8 
        NIfTot = 0
        nifd = 0 

        t_k_space_hubbard = .true. 
        t_lattice_model = .true. 

        ! setup G1 properly
        !todo 
        ! do a function like: which depending on the lattice sets up everything
        ! necessary for this type of lattice! yes!
        call setup_g1(lat) 

        ! setup nBasisMax and also the same for nBasisMax
        call setup_nbasismax(lat)
        !todo

    end subroutine init_k_space_unit_tests

    subroutine k_space_hubbard_test_driver() 
        ! with all the annying symmetry stuff to set up, testing the 
        ! k-space hubbard is really annoying.. 
        ! this is the main function which calls all the other tests 
       
        call run_test_case(setup_g1_test, "setup_g1_test")
        call run_test_case(setup_nbasismax_test, "setup_nbasismax_test")

        call init_k_space_unit_tests()
        call run_test_case(get_diag_helement_k_sp_hub_test, "get_diag_helement_k_sp_hub_test")
        call run_test_case(get_offdiag_helement_k_sp_hub_test, "get_offdiag_helement_k_sp_hub_test")
        call run_test_case(get_helement_k_space_hub_test, "get_helement_k_space_hub_test")
        call run_test_case(pick_spin_opp_elecs_test, "pick_spin_opp_elecs_test")
        call run_test_case(pick_from_cum_list_test, "pick_from_cum_list_test")
        call run_test_case(create_ab_list_hubbard_test, "create_ab_list_hubbard_test")
        call run_test_case(pick_ab_orbitals_hubbard_test, "pick_ab_orbitals_hubbard_test")
        call run_test_case(calc_pgen_k_space_hubbard_test, "calc_pgen_k_space_hubbard_test")
        call run_test_case(gen_excit_k_space_hub_test, "gen_excit_k_space_hub_test")
        call run_test_case(pick_three_opp_elecs_test, "pick_three_opp_elecs_test")
        call run_test_case(pick_spin_par_elecs_test, "pick_spin_par_elecs_test")
        call run_test_case(pick_a_orbital_hubbard_test, "pick_a_orbital_hubbard_test")
        call run_test_case(pick_bc_orbitals_hubbard_test,"pick_bc_orbitals_hubbard_test")
        call run_test_case(create_ab_list_par_hubbard_test, "create_ab_list_par_hubbard_test")
        call run_test_case(pick_ab_orbitals_par_hubbard_test, "pick_ab_orbitals_par_hubbard_test")
        call run_test_case(get_transferred_momentum_test, "get_transferred_momentum_test")
        call run_test_case(get_orb_from_kpoints_three_test, "get_orb_from_kpoints_three_test")
        call run_test_case(create_bc_list_hubbard_test, "create_bc_list_hubbard_test")
        call run_test_case(get_3_body_helement_ks_hub_test, "get_3_body_helement_ks_hub_test")
        call run_test_case(check_momentum_sym_test, "check_momentum_sym_test")
        call run_test_case(find_minority_spin_test, "find_minority_spin_test")
        call run_test_case(calc_pgen_k_space_hubbard_transcorr_test, "calc_pgen_k_space_hubbard_transcorr_test")
        call run_test_case(calc_pgen_k_space_hubbard_par_test, "calc_pgen_k_space_hubbard_transcorr")
        call run_test_case(calc_pgen_k_space_hubbard_triples_test, "calc_pgen_k_space_hubbard_triples_test")

    end subroutine k_space_hubbard_test_driver

    subroutine get_diag_helement_k_sp_hub_test

        print *, ""
        print *, "testing: get_diag_helement_k_sp_hub" 
        call assert_true(.false.)

        print *, "" 
        print *, "and now for 2-body transcorrelation: "
        t_trans_corr_2body = .true.
        call assert_true(.false.)

    end subroutine get_diag_helement_k_sp_hub_test

    subroutine get_offdiag_helement_k_sp_hub_test

        print *, ""
        print *, "testing: get_offdiag_helement_k_sp_hub" 
        call assert_true(.false.)

    end subroutine get_offdiag_helement_k_sp_hub_test

    subroutine get_helement_k_space_hub_test

        print *, ""
        print *, "testing: get_helement_k_space_hub_test" 
        call assert_true(.false.)
    end subroutine get_helement_k_space_hub_test

    subroutine pick_spin_opp_elecs_test
        print *, ""
        print *, "testing: pick_spin_opp_elecs"
        call assert_true(.false.)

    end subroutine pick_spin_opp_elecs_test

    subroutine pick_from_cum_list_test

        print *, ""
        print *, "testing: pick_from_cum_list"
        call assert_true(.false.)

    end subroutine pick_from_cum_list_test

    subroutine create_ab_list_hubbard_test

        print *, "" 
        print *, "testing: create_ab_list_hubbard"
        call assert_true(.false.)
        
        ! todo: i also have to do that for 2-body-transcorrelation, which 
        ! leads to parallel double excitations in the k-space hubbard 
        ! case -> check here if the get_orb_from_kpoints() function, works 
        ! correctly for ispn /= 2 and thub!

    end subroutine create_ab_list_hubbard_test

    subroutine calc_pgen_k_space_hubbard_test

        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard"
        call assert_true(.false.)

    end subroutine calc_pgen_k_space_hubbard_test

    subroutine gen_excit_k_space_hub_test

        print *, ""
        print *, "testing: gen_excit_k_space_hub" 
        call assert_true(.false.)

    end subroutine gen_excit_k_space_hub_test

    subroutine pick_three_opp_elecs_test

        print *, ""
        print *, "testing: pick_spin_opp_elecs"
        call assert_true(.false.)

    end subroutine pick_three_opp_elecs_test

    subroutine pick_spin_par_elecs_test

        print *, ""
        print *, "testing: pick_spin_par_elecs"
        call assert_true(.false.)
    end subroutine pick_spin_par_elecs_test

    subroutine pick_a_orbital_hubbard_test

        print *, ""
        print *, "testing: pick_a_orbital_hubbard "
        call assert_true(.false.) 

    end subroutine pick_a_orbital_hubbard_test

    subroutine pick_ab_orbitals_hubbard_test

        print *, "" 
        print *, "testing: pick_ab_orbitals_hubbard"
        call assert_true(.false.)

    end subroutine pick_ab_orbitals_hubbard_test

    subroutine pick_bc_orbitals_hubbard_test

        print *, ""
        print *, "testing: pick_bc_orbitals_hubbard"
        call assert_true(.false.)

    end subroutine pick_bc_orbitals_hubbard_test

    subroutine create_ab_list_par_hubbard_test

        print *, ""
        print *, "testing: create_ab_list_par_hubbard"

        call assert_true(.false.)

    end subroutine create_ab_list_par_hubbard_test

    subroutine pick_ab_orbitals_par_hubbard_test

        print *, ""
        print *, "testing: pick_ab_orbitals_par_hubbard"
        call assert_true(.false.)

    end subroutine pick_ab_orbitals_par_hubbard_test

    subroutine get_transferred_momentum_test

        print *, ""
        print *, "testing: get_transferred_momentum"
        call assert_true(.false.)

    end subroutine get_transferred_momentum_test

    subroutine get_orb_from_kpoints_three_test

        print *, ""
        print *, "testing: get_orb_from_kpoints_three: "
        call assert_true(.false.)

    end subroutine get_orb_from_kpoints_three_test

    subroutine create_bc_list_hubbard_test

        print *, ""
        print *, "testing: create_bc_list_hubbard"
        call assert_true(.false.)

    end subroutine create_bc_list_hubbard_test

    subroutine get_3_body_helement_ks_hub_test

        print *, "" 
        print *, "testing: get_3_body_helement_ks_hub" 
        call assert_true(.false.)

    end subroutine get_3_body_helement_ks_hub_test

    subroutine check_momentum_sym_test

        print *, ""
        print *, "testing: check_momentum_sym"
        call assert_true(.false.)

    end subroutine check_momentum_sym_test

    subroutine find_minority_spin_test

        print *, ""
        print *, "testing: find_minority_spin"
        call assert_true(.false.)

    end subroutine find_minority_spin_test
end program test_k_space_hubbard
