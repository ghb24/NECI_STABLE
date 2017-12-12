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
                          t_trans_corr, G1, nBasisMax, nOccBeta, nOccAlpha, & 
                          bhub, uhub, omega, trans_corr_param_2body, & 
                          t_trans_corr, t_trans_corr_2body, trans_corr_param, & 
                          thub, tpbc, treal, ttilt, TSPINPOLAR, & 
                          tCPMD, tVASP, tExch, tHphf
    use bit_rep_data, only: niftot, nifd
    use lattice_mod, only: lat, lattice
    use dsfmt_interface, only: dsfmt_init
    use OneEInts, only: GetTMatEl, tOneElecDiag, tCPMDSymTMat
    use procedure_pointers, only: get_umat_el
    use gen_coul_ueg_mod, only: get_hub_umat_el
    use IntegralsData, only: umat

    implicit none 

    integer :: failed_count 


    call init_fruit()
    call dsfmt_init(0)

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

        ! set the appropriate flags: 
        thub = .true.
        tpbc = .true. 
        treal = .false. 
        ! todo: also test for tilted lattices! since Ali likes them so much.. 
        TSPINPOLAR = .false. 
        tOneElecDiag = .false.

        ttilt =  .false. 
        tCPMDSymTMat = .false.
        tvasp = .false.
        tCPMD = .false. 

        tExch = .true. 

        thphf = .false.

        bhub = -1.0

        lat => lattice('chain', 4, 1, 1,.true.,.true.,.true.,'k-space')

        t_k_space_hubbard = .true. 
        t_lattice_model = .true. 

        ! setup nBasisMax and also the same for nBasisMax
        call setup_nbasismax(lat)

        ! setup G1 properly
        ! do a function like: which depending on the lattice sets up everything
        ! necessary for this type of lattice! yes!
        call setup_g1(lat) 

        ! also need the tmat ready.. 
        call setup_tmat_k_space(lat)

        ! and i also have to setup the symmetry table... damn.. 
        ! i have to setup umat also or
        uhub = 1.0
        omega = 8.0

        ! and i have to allocate umat.. 
        allocate(umat(1))
        umat = h_cast(real(uhub,dp)/real(omega,dp))

        get_umat_el => get_hub_umat_el

        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
    end subroutine init_k_space_unit_tests

    subroutine k_space_hubbard_test_driver() 
        ! with all the annying symmetry stuff to set up, testing the 
        ! k-space hubbard is really annoying.. 
        ! this is the main function which calls all the other tests 
       
        call run_test_case(setup_g1_test, "setup_g1_test")
        call run_test_case(setup_nbasismax_test, "setup_nbasismax_test")
        call run_test_case(setup_tmat_k_space_test, "setup_tmat_k_space_test")

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
        call run_test_case(calc_pgen_k_space_hubbard_par_test, "calc_pgen_k_space_hubbard_par_test")
        call run_test_case(calc_pgen_k_space_hubbard_triples_test, "calc_pgen_k_space_hubbard_triples_test")
        call run_test_case(make_triple_test, "make_triple_test")
        call run_test_case(make_double_test, "make_double_test")
        call run_test_case(make_single_test, "make_single_test")
        call run_test_case(three_body_transcorr_fac_test, "three_body_transcorr_fac_test")
        call run_test_case(two_body_transcorr_factor_test, "two_body_transcorr_factor_test")
        call run_test_case(epsilon_kvec_test, "epsilon_kvec_test")
        call run_test_case(same_spin_transcorr_factor_test, "same_spin_transcorr_factor_test")
        call run_test_case(get_one_body_diag_test, "get_one_body_diag_test")

    end subroutine k_space_hubbard_test_driver

    subroutine setup_g1_test

        class(lattice), pointer :: ptr 
        integer :: i 
        
        print *, "" 
        print *, "testing: setup_g1" 
        ptr => lattice('chain', 4, 1, 1,.true.,.true.,.true.,'k-space')
        nBasis = 8

        tpbc = .true. 
        treal = .false. 
        ttilt = .false.

        call setup_g1(ptr) 

        ! check ms: 
        do i = 1, 7, 2 
            call assert_equals(-1, G1(i)%ms)
        end do
        do i = 2, 8, 2
            call assert_equals(1, G1(i)%ms)
        end do

        call assert_equals([(0,i=1,8)], G1(1:8)%k(2),8)
        call assert_equals([(0,i=1,8)], G1(1:8)%k(3),8)
        call assert_equals([-1,-1,0,0,1,1,2,2],G1(1:8)%k(1),8)

        deallocate(G1)
        nBasisMax = 0
        nbasis = -1

    end subroutine setup_g1_test

    subroutine setup_nbasismax_test

        use SystemData, only: nBasisMax
        class(lattice), pointer :: ptr 
        print *, "" 
        print *, "testing: setup_nbasismax"
        ptr => lattice('chain', 4, 1, 1,.true.,.true.,.true.,'k-space')

        nBasis = 8
        tpbc = .true. 
        treal = .false. 

        call setup_nbasismax(ptr)

        call assert_equals(0, nBasisMax(1,3))
        call assert_equals(0, nBasisMax(1,4))
        call assert_equals(1, nBasisMax(2,4))
        call assert_equals(4, nBasisMax(1,5))
        call assert_equals(1, nBasisMax(2,5))
        call assert_equals(1, nBasisMax(3,5))

        nBasisMax = 0
        nbasis = -1

    end subroutine setup_nbasismax_test

    subroutine setup_tmat_k_space_test

        class(lattice), pointer :: ptr 
        integer :: i 

        ptr => lattice('chain', 4, 1, 1, .true., .true., .true.,'k-space')
        print *, "" 
        print *, "testing: setup_tmat_k_space"
        tpbc = .true. 
        treal = .false. 
        tOneElecDiag = .false. 
        ttilt = .false. 
        bhub = 1.0
        nBasis = 8 

        call setup_tmat_k_space(ptr)

        call assert_equals(h_cast(0.0_dp), GetTMatEl(1,2))

        call assert_equals(h_cast(2.0_dp*cos(-PI/2.0_dp)), GetTMatEl(1,1))
        call assert_equals(h_cast(2.0_dp*cos(-PI/2.0_dp)), GetTMatEl(2,2))
        call assert_equals(h_cast(2.0_dp*cos(PI/2.0_dp)), GetTMatEl(5,5))
        call assert_equals(h_cast(2.0_dp*cos(PI/2.0_dp)), GetTMatEl(6,6))
        call assert_equals(h_cast(2.0_dp), GetTMatEl(3,3))
        call assert_equals(h_cast(2.0_dp), GetTMatEl(4,4))
        call assert_equals(h_cast(-2.0_dp), GetTMatEl(7,7))
        call assert_equals(h_cast(-2.0_dp), GetTMatEl(8,8))

        deallocate(G1) 
        nBasisMax = 0
        deallocate(tmat2d)
        nbasis = -1

    end subroutine setup_tmat_k_space_test

    subroutine get_diag_helement_k_sp_hub_test

        print *, ""
        print *, "testing: get_diag_helement_k_sp_hub" 
        umat = 0.0_dp
        call assert_equals(h_cast(-4.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

        umat = 2*uhub/omega
        call assert_equals(h_cast(-3.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

        print *, "" 
        print *, "and now for 2-body transcorrelation: "
        t_trans_corr_2body = .true.
        ! test it for 0 transcorrelation
        trans_corr_param_2body = 0.0_dp
        call assert_equals(h_cast(-3.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

        trans_corr_param_2body = 1.0_dp 

        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
        print *, get_diag_helement_k_sp_hub([1,2,3,4])

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

        integer, allocatable :: nI(:)
        integer :: elecs(2)
        real(dp) :: p_elec

        nel = 2 
        nOccBeta = 1
        nOccAlpha = 1
        allocate(nI(nel))

        print *, ""
        print *, "testing: pick_spin_opp_elecs"
        nI = [1,2] 
        call pick_spin_opp_elecs(nI, elecs, p_elec)

        call assert_equals(1.0_dp, p_elec) 
        if (elecs(1) == 1) then 
            call assert_equals(2, elecs(2)) 
        else if (elecs(1) == 2) then
            call assert_equals(1, elecs(2))
        end if

        nel = -1
        nOccBeta = -1
        nOccAlpha = -1

    end subroutine pick_spin_opp_elecs_test

    subroutine pick_from_cum_list_test

        integer :: ind
        real(dp) :: pgen
        print *, ""
        print *, "testing: pick_from_cum_list"
        call pick_from_cum_list([0.0_dp,1.0_dp],1.0_dp, ind, pgen)

        call assert_equals(2, ind) 
        call assert_equals(1.0_dp, pgen) 

        call pick_from_cum_list([1.0_dp,1.0_dp],1.0_dp, ind, pgen)

        call assert_equals(1, ind) 
        call assert_equals(1.0_dp, pgen) 

        call pick_from_cum_list([1.0_dp,2.0_dp],2.0_dp, ind, pgen)
        call assert_equals(0.5_dp, pgen) 

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

        integer :: elecs(3), sum_ms
        real(dp) :: p_elec 

        nel = 3 
        nOccBeta = 2 
        nOccAlpha = 1 

        print *, ""
        print *, "testing: pick_three_opp_elecs"
        call pick_three_opp_elecs([1,2,3], elecs, p_elec, sum_ms)
        call assert_equals(1.0_dp, p_elec)
        call assert_equals(-1, sum_ms) 
        call assert_equals(6, sum(elecs))

        nOccAlpha = 2
        nOccBeta = 1 

        call pick_three_opp_elecs([1,2,4], elecs, p_elec, sum_ms)
        call assert_equals(1.0_dp, p_elec)
        call assert_equals(1, sum_ms) 
        call assert_equals(6, sum(elecs))

        nel = 5
        nOccAlpha = 4 
        call pick_three_opp_elecs([1,2,4,6,8], elecs, p_elec, sum_ms) 
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(1, sum_ms) 
        call assert_true(any(elecs == 1))

        nOccBeta = 4 
        nOccAlpha = 1
        call pick_three_opp_elecs([1,3,5,7,8], elecs, p_elec, sum_ms) 
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(-1, sum_ms) 
        call assert_true(any(elecs == 8))

        nel = -1
        nOccBeta = -1 
        nOccAlpha = -1

    end subroutine pick_three_opp_elecs_test

    subroutine pick_spin_par_elecs_test

        integer :: elecs(2), ispn
        real(dp) :: p_elec
        
        print *, ""
        print *, "testing: pick_spin_par_elecs"
        nel = 2 
        nOccBeta = 2
        call pick_spin_par_elecs([1,3],elecs,p_elec, ispn)
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(1, ispn) 
        call assert_equals(3, sum(elecs))

        nOccAlpha = 2 
        call pick_spin_par_elecs([2,4],elecs,p_elec, ispn)
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(3, ispn) 
        call assert_equals(3, sum(elecs))

        nel = 4 
        call pick_spin_par_elecs([1,2,3,4], elecs, p_elec, ispn) 
        call assert_equals(0.5_dp, p_elec) 
        if (ispn == 1) then 
            call assert_equals(4, sum(elecs))
        else if (ispn == 3) then 
            call assert_equals(6, sum(elecs))
        end if
        
        nel = -1 
        nOccBeta = -1 
        nOccAlpha = -1

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

        integer :: nel, ex(2,3)
        integer, allocatable :: nI(:) 
        logical :: tpar

        tpar = .false.

        nel = 4 
        allocate(ni(nel))

        ni = [1,2,3,4]

        print *, "" 
        print *, "testing: get_3_body_helement_ks_hub" 
        ex(1,:) = [1,3,5] 
        ex(2,:) = [2,4,6]

        ! here spin does not fit:
        call assert_equals(h_cast(0.0_dp), get_3_body_helement_ks_hub(ni,ex,tpar))
        ex(1,:) = [1,2,3]
        ex(2,:) = [4,5,6]
        call assert_equals(h_cast(0.0_dp), get_3_body_helement_ks_hub(ni,ex,tpar))
        ex(2,:) = [5,6,7]

        ! and here momentum does not fit
        call assert_equals(h_cast(0.0_dp), get_3_body_helement_ks_hub(ni,ex,tpar))

        ! and here it should fit.
        ex(1,:) = [3,6,7]
        ex(2,:) = [1,2,5]
        call assert_equals(h_cast(-4*three_body_prefac), get_3_body_helement_ks_hub(ni,ex,tpar))

        ! maybe do more tests..

    end subroutine get_3_body_helement_ks_hub_test

    subroutine check_momentum_sym_test

        print *, ""
        print *, "testing: check_momentum_sym"
        ! use the already setup up 4 site chain.. the input to this is 
        ! with spin-orbitals.. or no.. it is with spatial orbs! no it is 
        ! spin-orbital! but the spin is also checked for symmetry! 
        ! although it is not only momentum symmetry! it is also 
        ! spin symmetry!! 
        call assert_true(check_momentum_sym([1],[1]))
        call assert_true(.not.check_momentum_sym([1],[2]))
        call assert_true(check_momentum_sym([2],[2]))
        call assert_true(.not.check_momentum_sym([3],[1]))

        call assert_true(check_momentum_sym([1,2],[2,1]))
        call assert_true(check_momentum_sym([2,2],[2,2]))

        ! and it takes variable sizes of input..
        call assert_true(check_momentum_sym([1,5],[3,3]))

        call assert_true(check_momentum_sym([6,6],[4,8]))
        call assert_true(.not.check_momentum_sym([5,5],[4,8]))

        call assert_true(check_momentum_sym([5,8],[1,4]))
    end subroutine check_momentum_sym_test

    subroutine find_minority_spin_test

        print *, ""
        print *, "testing: find_minority_spin"
        call assert_equals(1, find_minority_spin([1,2,4]))
        call assert_equals(3, find_minority_spin([3,2,4]))

        call assert_equals(2, find_minority_spin([1,2,3]))
        call assert_equals(2, find_minority_spin([3,2,5]))

    end subroutine find_minority_spin_test

    subroutine calc_pgen_k_space_hubbard_transcorr_test

        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_transcorr"
        call assert_true(.false.)

    end subroutine calc_pgen_k_space_hubbard_transcorr_test

    subroutine calc_pgen_k_space_hubbard_par_test

        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_par"
        call assert_true(.false.)

    end subroutine calc_pgen_k_space_hubbard_par_test

    subroutine calc_pgen_k_space_hubbard_triples_test

        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_triples"
        call assert_true(.false.)

    end subroutine calc_pgen_k_space_hubbard_triples_test

    subroutine make_triple_test


        integer, allocatable :: nI(:), nJ(:)
        integer :: ex(2,3)
        logical :: tpar, tpar_2

        nel = 3 

        allocate(nI(nel))
        allocate(nJ(nel))

        
        print *, "" 
        print *, "testing: make_triple" 
        print *, "testing implicitly: FindExcitDet!"

        nI = [1,2,3] 
        call make_triple(nI,nJ,[1,2,3],[4,5,7],ex,tpar) 
        call assert_equals([4,5,7], nJ, 3)
        call assert_equals([1,2,3], ex(1,:),3)
        call assert_equals([4,5,7], ex(2,:),3)
        call assert_true(.not.tpar)

        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        ! and now more complicated stuff:
        nI = [1,2,3]
        call make_triple(nI,nJ,[1,3,2],[7,5,4],ex,tpar) 
        call assert_equals([4,5,7], nJ, 3)
        call assert_equals([1,2,3], ex(1,:),3)
        call assert_equals([4,5,7], ex(2,:),3)
        call assert_true(.not.tpar)

        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        nI = [1,2,5]
        call make_triple(nI,nJ,[3,1,2],[3,4,7],ex,tpar) 
        call assert_equals([3,4,7], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([3,4,7], ex(2,:),3)
        call assert_true(.not.tpar)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,5]

        call make_triple(nI,nJ,[3,2,1],[8,7,3],ex,tpar) 
        call assert_equals([3,7,8], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([3,7,8], ex(2,:),3)
        call assert_true(.not.tpar)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,5]

        call make_triple(nI,nJ,[3,2,1],[4,7,9],ex,tpar) 
        call assert_equals([4,7,9], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([4,7,9], ex(2,:),3)
        call assert_true(.not.tpar)
        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,5]

        call make_triple(nI, nJ, [1,2,3], [3,4,7], ex, tpar) 
        call assert_true(.not.tpar)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,5]

        nel = 4 

        deallocate(nJ); allocate(NJ(nel))
        deallocate(nI); allocate(nI(nel))
        
        nI = [1,2,5,7]
        call make_triple(nI,nJ,[1,2,3],[3,6,9],ex,tpar)
        call assert_true(tpar)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        nel = -1

    end subroutine make_triple_test

    subroutine make_double_test

        use get_excit, only: make_double
        integer, allocatable :: nJ(:),ni(:)
        integer :: ex(2,2) 
        logical :: tpar, tpar_2

        print *, "" 
        print *, "testing: make_double" 
        print *, "to be consistent with the sign conventions! "

        ! to test this really strange sign convention also call all the other 
        ! routines here, which test sign.. 

        nel = 2
        allocate(nJ(nel))
        allocate(ni(nel)); 
        ni = [1,2]

        call make_double([1,2],nJ, 1,2, 3,4, ex, tpar)
        call assert_equals([3,4], nJ, 2) 
        call assert_equals([1,2],ex(1,:),2)
        call assert_equals(reshape([1,3,2,4],[2,2]),ex, 2,2)
        call assert_true(.not. tpar)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        call make_double([1,2],nJ, 1,2, 5,4, ex, tpar)
        call assert_true(.not. tpar)
        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        call make_double([1,2],nJ, 1,2, 3,6, ex, tpar) 
        call assert_true(.not. tpar)
        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        nel = 3 
        deallocate(nJ); allocate(nJ(nel))
        nI = [1,2,4]

        call make_double([1,2,4],nJ,1,2,5,6,ex,tpar)
        call assert_equals([4,5,6], nJ, 3)
        call assert_true(.not. tpar)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,4]

        call make_double([1,2,4],nj,1,2,3,6,ex,tpar)
        call assert_true(tpar)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,4]

        call make_double([1,2,4],nJ, 1, 2, 6, 7, ex, tpar)
        call assert_true(.not. tpar)
        nI = [1,2,4]

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        call make_double([1,2,3],nJ, 1, 2, 4, 7, ex, tpar)
        call assert_true(.not. tpar)
        nI = [1,2,3]

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        nel = -1


    end subroutine make_double_test

    subroutine make_single_test

        print *, "" 
        print *, "testing: make_single" 
        print *, "for consistency in the sign"

        call assert_true(.false.) 

    end subroutine make_single_test

    subroutine three_body_transcorr_fac_test

        integer :: p(3), q(3), k(3)
        real(dp) :: test 

        nel = 4
        nOccBeta = 2
        nOccAlpha = 2

        print *, "" 
        print *, "testing: three_body_transcorr_fac"
        p = 0 
        q = 0
        k = 0

        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,1))
        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,-1))

        q(1) = 1 
        k(1) = 2 
        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,1),1e-12)

        p(1) = 2 
        q(1) = 0 
        k(1) = 1 
        call assert_equals(h_cast(0.0_dp), &
            three_body_transcorr_fac([1,2,3,4], p,q,k,1))

        ! is there a non-zero combination? 
        p(1) = 0
        q(1) = 2
        k(1) = 1
        call assert_equals(h_cast(three_body_prefac*8), & 
            three_body_transcorr_fac([1,2,3,4],p,q,k,1))
        call assert_equals(h_cast(three_body_prefac*8), & 
            three_body_transcorr_fac([1,2,3,4],p,q,k,-1))
        
        nOccBeta = 4 
        nOccAlpha = 0
        call assert_equals(h_cast(0.0_dp), & 
            three_body_transcorr_fac([1,3,5,7],p,q,k,-1))

        call assert_equals(h_cast(three_body_prefac*8), & 
            three_body_transcorr_fac([1,3,5,7],p,q,k,1))

        nOccBeta = 2
        nOccAlpha = 2

    end subroutine three_body_transcorr_fac_test

    subroutine two_body_transcorr_factor_test

        integer :: p(3), k(3)

        print *, "" 
        print *, "testing: two_body_transcorr_factor" 
        p = 0
        k = 0

        call assert_equals(h_cast(-4.0_dp/omega * (cosh(trans_corr_param_2body) - 1)), & 
            two_body_transcorr_factor(p,k))

        p(1) = 1
        call assert_equals(h_cast(0.0_dp), two_body_transcorr_factor(p,k),1.e-12)
        call assert_equals(h_cast(-4.0_dp/omega * sinh(trans_corr_param_2body)), & 
            two_body_transcorr_factor([2,0,0],[2,0,0]))

        call assert_equals(h_cast(4.0_dp/omega * sinh(trans_corr_param_2body)), & 
            two_body_transcorr_factor([0,0,0],[2,0,0]))

        call assert_equals(h_cast(-2.0_dp/omega *(exp(-trans_corr_param_2body)-1)), & 
            two_body_transcorr_factor([0,0,0],[1,0,0]),1.e-12)


        call assert_equals(h_cast(4.0_dp/omega * (cosh(trans_corr_param_2body) - 1)), & 
            two_body_transcorr_factor([2,0,0],[0,0,0]))

    end subroutine two_body_transcorr_factor_test

    subroutine epsilon_kvec_test

        integer :: k(3)

        ! depending on the lattice dimension.. 
        ! and also the tilted has other values or?? 
        print *, "" 
        print *, "testing: epsilon_kvec"
        k = 0 
        call assert_equals(h_cast(2.0_dp), epsilon_kvec(k))
        call assert_equals(h_cast(0.0_dp), epsilon_kvec([1,0,0]),1e-12)
        call assert_equals(h_cast(0.0_dp), epsilon_kvec([-1,0,0]),1e-12)
        call assert_equals(h_cast(-2.0_dp), epsilon_kvec([2,0,0]))

    end subroutine epsilon_kvec_test

    subroutine same_spin_transcorr_factor_test

        integer :: k(3) 
        integer, allocatable :: nI(:) 

        nel = 4
        allocate(nI(nel))

        nI = [1,2,3,4] 

        print *, "" 
        print *, "testing: same_spin_transcorr_factor" 

        call assert_equals(h_cast(three_body_prefac*4), same_spin_transcorr_factor([1,2,3,4],[0,0,0],1),1.e-12)
        call assert_equals(h_cast(three_body_prefac*4), same_spin_transcorr_factor([1,2,3,4],[0,0,0],-1),1.e-12)
        call assert_equals(h_cast(0.0_dp), same_spin_transcorr_factor([1,2,3,4],[1,0,0],-1),1.e-12)
        call assert_equals(h_cast(-three_body_prefac*4), same_spin_transcorr_factor([1,2,3,4],[2,0,0],-1),1.e-12)

        call assert_equals(h_cast(0.0_dp), same_spin_transcorr_factor([1,3,5,7],[0,0,0],1))
        call assert_equals(h_cast(0.0_dp), same_spin_transcorr_factor([1,3,5,7],[0,0,0],-1))

    end subroutine same_spin_transcorr_factor_test

    subroutine get_one_body_diag_test

        integer, allocatable :: nI(:)
        class(lattice), pointer :: ptr
        integer :: i 

        nel = 4

        allocate(nI(nel))
        nI = [1,2,3,4]
        
        ! i think i also want to use the lattice functionality in the 
        ! k-space hubbard model.. but i still have to think about how to do that!
        ! allow an additional input flag! 

!         ptr => lattice('chain', 4, 1, 1, .true., .true., .true.,'k-space')

!         call setup_tmat_k_space(ptr) 
        print *, "" 
        print *, "testing: get_one_body_diag" 
        ! the spin = 1 means i want the diagonal contribution of the alpha 
        ! electrons! 
        call assert_equals(h_cast(2.0_dp), get_one_body_diag(nI,1))
        call assert_equals(h_cast(2.0_dp), get_one_body_diag(nI,-1))
        call assert_equals(h_cast(4.0_dp), get_one_body_diag(nI))

        nI = [1,2,5,6]
        call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,1),1.e-8)
        call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,-1),1e-8)
        call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI),1.e-8)

        nI = [1,3,5,7]
        call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,1))
        call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI,-1))
        call assert_equals(h_cast(0.0_dp), get_one_body_diag(nI))


    end subroutine get_one_body_diag_test

end program test_k_space_hubbard
