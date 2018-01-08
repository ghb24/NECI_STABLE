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
                          tCPMD, tVASP, tExch, tHphf, tNoSymGenRandExcits, tKPntSym
    use bit_rep_data, only: niftot, nifd
    use lattice_mod, only: lat, lattice
    use dsfmt_interface, only: dsfmt_init
    use OneEInts, only: GetTMatEl, tOneElecDiag, tCPMDSymTMat
    use procedure_pointers, only: get_umat_el
    use gen_coul_ueg_mod, only: get_hub_umat_el
    use IntegralsData, only: umat
    use DetBitOps, only: EncodeBitDet
    use fcimcdata, only: pDoubles, pParallel
    use DetBitOps, only: ilut_lt, ilut_gt
    use sort_mod, only: sort

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

        tKPntSym = .true.
        ttilt =  .false. 
        tCPMDSymTMat = .false.
        tvasp = .false.
        tCPMD = .false. 

        tExch = .true. 

        thphf = .false.

        bhub = -1.0

        tNoSymGenRandExcits = .true. 

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

!         call setup_kPointToBasisFn(lat)

        call setup_k_space_hub_sym(lat) 

        ! and i also have to setup the symmetry table... damn.. 
        ! i have to setup umat also or
        uhub = 1.0
        omega = 4.0

        ! and i have to allocate umat.. 
        allocate(umat(1))
        umat = h_cast(real(uhub,dp)/real(omega,dp))

        get_umat_el => get_hub_umat_el

        trans_corr_param_2body = 0.1
        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)

    end subroutine init_k_space_unit_tests

    subroutine k_space_hubbard_test_driver() 
        ! with all the annying symmetry stuff to set up, testing the 
        ! k-space hubbard is really annoying.. 
        ! this is the main function which calls all the other tests 
       
        call run_test_case(setup_g1_test, "setup_g1_test")
        call run_test_case(setup_nbasismax_test, "setup_nbasismax_test")
        call run_test_case(setup_tmat_k_space_test, "setup_tmat_k_space_test")
        call run_test_case(setup_kPointToBasisFn_test, "setup_kPointToBasisFn_test")

        call init_k_space_unit_tests()
        call run_test_case(test_3e_4orbs_par, "test_3e_4orbs_par")
        call run_test_case(test_3e_4orbs_trip, "test_3e_4orbs_trip")
        call run_test_case(test_4e_ms1, "test_4e_ms1")
        call run_test_case(test_4e_ms0_mom_1, "test_4e_ms0_mom_1")
        call run_test_case(test_3e_ms1, "test_3e_ms1")

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
        call run_test_case(get_transferred_momenta_test, "get_transferred_momenta_test")
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
        call run_test_case(three_body_transcorr_fac_test, "three_body_transcorr_fac_test")
        call run_test_case(two_body_transcorr_factor_test, "two_body_transcorr_factor_test")
        call run_test_case(epsilon_kvec_test, "epsilon_kvec_test")
        call run_test_case(same_spin_transcorr_factor_test, "same_spin_transcorr_factor_test")
        call run_test_case(get_one_body_diag_test, "get_one_body_diag_test")
        call run_test_case(gen_parallel_double_hubbard_test, "gen_parallel_double_hubbard_test")
        call run_test_case(gen_triple_hubbard_test, "gen_triple_hubbard_test")
        call run_test_case(gen_excit_k_space_hub_test_stochastic, "gen_excit_k_space_hub_test_stochastic") 
        call run_test_case(gen_excit_k_space_hub_transcorr_test_stoch, "gen_excit_k_space_hub_transcorr_test_stoch")

    end subroutine k_space_hubbard_test_driver

    subroutine test_3e_4orbs_par

        integer :: hilbert_nI(3,6), i, j, work(18), info
        HElement_t(dp) :: hamil(6,6), hamil_trancorr(6,6), tmp_hamil(6,6)
        real(dp) :: ev_real(6), ev_cmpl(6), left_ev(1,6), right_ev(1,6)
        real(dp) :: t_mat(6,6), trans_mat(6,6)

        nOccBeta = 2 
        nOccAlpha = 1
        nel = 3

        hilbert_nI(:,1) = [1,4,5]
        hilbert_nI(:,2) = [2,3,5]
        hilbert_nI(:,3) = [1,3,6]
        hilbert_nI(:,4) = [3,7,8]
        hilbert_nI(:,5) = [1,2,7]
        hilbert_nI(:,6) = [5,6,7]

        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "

        hamil = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil)

        t_mat = get_tranformation_matrix(hamil,2)

        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))
        
        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.
        
        hamil_trancorr = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        print *, "transformed hamiltonian" 
        call print_matrix(trans_mat)

        print *, "eigen-values: ", calc_eigenvalues(trans_mat)
        
    end subroutine test_3e_4orbs_par

    subroutine test_3e_4orbs_trip

        integer :: hilbert_nI(3,6), i, j, work(18), info
        HElement_t(dp) :: hamil(6,6), hamil_trancorr(6,6), tmp_hamil(6,6)
        real(dp) :: ev_real(6), ev_cmpl(6), left_ev(1,6), right_ev(1,6)
        real(dp) :: t_mat(6,6), trans_mat(6,6)

        nOccBeta = 1 
        nOccAlpha = 2
        nel = 3

        hilbert_nI(:,1) = [2,3,4]
        hilbert_nI(:,2) = [1,2,6]
        hilbert_nI(:,3) = [4,5,8]
        hilbert_nI(:,4) = [4,6,7]
        hilbert_nI(:,5) = [2,7,8]
        hilbert_nI(:,6) = [3,6,8]

        ! first create the non-transcorrelated Hamiltonian
        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "
        hamil = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil)
        
        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.
        hamil_trancorr = create_hamiltonian(hilbert_nI)

        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        t_mat = get_tranformation_matrix(hamil,2) 
        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        print *, "transformed hamiltonian" 
        call print_matrix(trans_mat)

        print *, "eigen-values: ", calc_eigenvalues(trans_mat)
        
    end subroutine test_3e_4orbs_trip

    subroutine test_4e_ms0_mom_1

        integer :: hilbert_nI(4,8)
        HElement_t(dp) :: hamil(8,8),hamil_trancorr(8,8), tmp_hamil(8,8) 
        real(dp) :: ev_real(8), ev_cmpl(8), left_ev(1,8), right_ev(1,8) 
        integer :: work(24), info, n
        real(dp) :: t_mat(8,8), trans_mat(8,8) 

        nOccBeta = 2 
        nOccAlpha = 2
        nel = 4 

        hilbert_nI(:,1) = [3,4,5,8] 
        hilbert_nI(:,2) = [1,2,4,5] 
        hilbert_nI(:,3) = [3,4,6,7]
        hilbert_nI(:,4) = [1,4,7,8]
        hilbert_nI(:,5) = [1,5,6,8]
        hilbert_nI(:,6) = [2,3,7,8]
        hilbert_nI(:,7) = [1,2,3,6]
        hilbert_nI(:,8) = [2,5,6,7]

        t_trans_corr_2body = .false. 

        hamil = create_hamiltonian(hilbert_nI)

        t_mat = get_tranformation_matrix(hamil,4) 

        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        t_trans_corr_2body = .true. 

        hamil_trancorr = create_hamiltonian(hilbert_nI) 

        print *, "un-correlated hamiltonian: "
        call print_matrix(hamil) 

        print *, "eigenvalues: ", calc_eigenvalues(hamil)

        print *, "neci-correlated hamiltonian: "
        call print_matrix(hamil_trancorr)
        print *, "eigenvalues: ", calc_eigenvalues(hamil_trancorr)

        print *, "transformed hamiltonian" 
        call print_matrix(trans_mat)
        print *, "eigenvalues: ", calc_eigenvalues(trans_mat)

    end subroutine test_4e_ms0_mom_1

    function calc_eigenvalues(matrix) result(e_values)
        real(dp) :: matrix(:,:) 
        real(dp) :: e_values(size(matrix,1))

        integer :: n, info, work(3*size(matrix,1))
        real(dp) :: tmp_matrix(size(matrix,1),size(matrix,2)),dummy_val(size(matrix,1))
        real(dp) :: dummy_vec_1(1,size(matrix,1)), dummy_vec_2(1,size(matrix,1))

        n = size(matrix,1)

        tmp_matrix = matrix 

        call dgeev('N','N', n, tmp_matrix, n, e_values, &
            dummy_val, dummy_vec_1,1,dummy_vec_2,1,work,3*n,info)

    end function calc_eigenvalues

    function create_hamiltonian(list_nI) result(hamil) 
        ! quite specific hamiltonian creation for my tests.. 
        integer, intent(in) :: list_nI(:,:) 
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j 

        do i = 1, size(list_nI,2)
            do j = 1, size(list_nI,2)
                hamil(i,j) = get_helement_k_space_hub(list_nI(:,j),list_nI(:,i))
            end do
        end do

    end function create_hamiltonian

    subroutine print_matrix(matrix) 
        ! print a 2-D real matrix 
        real(dp), intent(in) :: matrix(:,:)
        
        integer :: i 

        do i = 1, size(matrix,1)
            print *, matrix(i,:)
        end do

    end subroutine print_matrix

    subroutine test_4e_ms1

        integer :: hilbert_nI(4,4), i, j, three_e(3,3) 
        integer(n_int) :: hilbert_ilut(0:niftot,4)
        HElement_t(dp) :: hamil(4,4), hamil_trancorr(4,4), tmp_hamil(4,4)
        HElement_t(dp) :: tmp_3(3,3), hamil_3(3,3), hamil_3_trans(3,3)
        real(dp) :: ev_real(4), ev_cmpl(4)
        real(dp) :: left_ev(1,4), right_ev(1,4)
        real(dp) :: shift_12, shift_34, two_body, three_body, three_body_1, three_body_2
        real(dp) :: two_body_1, two_body_2
        integer :: work(12), info
        real(dp) :: t_mat(4,4), trans_mat(4,4)

        nOccBeta = 3 
        nOccAlpha = 1 
        nel = 4

        hilbert_nI(:,1) = [1,2,3,5]
        hilbert_nI(:,2) = [3,4,5,7]
        hilbert_nI(:,3) = [1,3,7,8] 
        hilbert_nI(:,4) = [1,5,6,7]

        ! first create the non-transcorrelated Hamiltonian
        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "
        
        hamil = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil)

        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.

        hamil_trancorr = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        t_mat = get_tranformation_matrix(hamil,3) 

        trans_mat = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        print *, "transformed hamiltonian: "
        call print_matrix(trans_mat)

        print *, "eigen-values: ", calc_eigenvalues(trans_mat)

        ! i know now, where there is an error-- between those 2 matrix 
        ! elements: 
!         print *, "============== Excitation:  ======================"
!         print *, "nI:", hilbert_nI(:,1)
!         print *, "nJ:", hilbert_nI(:,2)
!         print *, "H_ij: ", hamil_trancorr(1,2), hamil_trancorr(2,1)
! 
!         ! the excitation is from (1,2) -> (4,7)
!         print *, "excitation: (1,2) -> (4,7)" 
!         ! ex: should be + 
!         ! 1 2 
!         ! 4 7 
!         ! total 1 2 3 5 should be - 
! 
!         two_body_1 = two_body_transcorr_factor(G1(1)%k, G1(7)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(2)%k, G1(4)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,1), G1(1)%k,G1(2)%k,G1(7)%k,-1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,1),G1(2)%k,G1(1)%k,G1(4)%k,1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
!         print *, "excitation: (2,5) -> (7,8)" 
!         ! ex: should be + 
!         ! 2 5 
!         ! 7 8 
!         ! total 1 2 3 5 should be - 
! 
!         two_body_1 = two_body_transcorr_factor(G1(2)%k, G1(8)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(5)%k, G1(7)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,1), G1(2)%k,G1(5)%k,G1(8)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,1),G1(5)%k,G1(2)%k,G1(7)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
!         print *, "excitation: (2,3) -> (6,7)" 
!         ! ex: should be - 
!         ! 2 3 
!         ! 6 7 
!         ! total 1 2 3 5 should be + 
! 
!         two_body_1 = two_body_transcorr_factor(G1(2)%k, G1(6)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(3)%k, G1(7)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,1), G1(2)%k,G1(3)%k,G1(6)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,1),G1(3)%k,G1(2)%k,G1(7)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
!         print *, "excitation: (4,5) -> (1,8)" 
!         ! ex: should be - 
!         ! 4 5 
!         ! 1 8 
!         ! total 3 4 5 7 should be + 
! 
!         two_body_1 = two_body_transcorr_factor(G1(4)%k, G1(8)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(5)%k, G1(1)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,2), G1(4)%k,G1(5)%k,G1(8)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,2),G1(5)%k,G1(4)%k,G1(1)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
! 
!         print *, "excitation: (3,4) -> (1,6)" 
!         ! ex: should be - 
!         ! 3 4 
!         ! 1 6 
!         ! total 3 4 5 7 should be - 
! 
!         two_body_1 = two_body_transcorr_factor(G1(4)%k, G1(6)%k)! + & 
!         two_body_2 = two_body_transcorr_factor(G1(3)%k, G1(1)%k)
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,2), G1(4)%k,G1(3)%k,G1(6)%k,1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,2),G1(3)%k,G1(4)%k,G1(1)%k,-1)
! 
!         print *, "two-body: ", two_body_1, two_body_2
!         print *, "three-body: ", three_body_1, three_body_2
! 
! 
!         print *, "============== Excitation:  ======================"
!         print *, "nI: ", hilbert_nI(:,3)
!         print *, "nJ: ", hilbert_nI(:,4)
! 
!         print *, "H_ij: ", hamil_trancorr(3,4), hamil_trancorr(4,3)
! 
!         ! the excitation is from (3,8) -> (5,6)
!         print *, "excitation: (3,8) -> (5,6)"
!         ! ex: should have - sign 
!         ! 3 8 
!         ! 5 6
!         ! sign in total: 1 3 7 8 should be -
!         two_body = two_body_transcorr_factor(G1(3)%k, G1(5)%k) + & 
!                 two_body_transcorr_factor(G1(8)%k, G1(6)%k)
! 
!         print *, "two-body: ", two_body
! 
!         three_body_1 = three_body_transcorr_fac(hilbert_nI(:,3), G1(3)%k,G1(8)%k,G1(5)%k,-1)! + & 
!         three_body_2 = three_body_transcorr_fac(hilbert_nI(:,3),G1(8)%k,G1(3)%k,G1(6)%k,1)
! 
!         print *, "three-body: ", three_body_1, three_body_2

! 
!         print *, "and for k=0"
!         hilbert_nI(:,1) = [1,3,4,5]
!         hilbert_nI(:,2) = [1,2,3,7]
!         hilbert_nI(:,3) = [3,5,6,7]
!         hilbert_nI(:,4) = [1,5,7,8]
! 
!         ! first create the non-transcorrelated Hamiltonian
!         t_trans_corr_2body = .false. 
!         print *, "un-transcorrelated Hamiltonian: "
!         do i = 1, 4 
!             do j = 1, 4 
!                 hamil(i,j) = get_helement_k_space_hub(hilbert_nI(:,i),hilbert_nI(:,j))
!             end do
!             print *, hamil(i,:)
!         end do
!         
!         ! use the lapack routines to solve these quickly.. 
!         tmp_hamil = hamil
!         call dgeev('N','N',4,tmp_hamil,4,ev_real,ev_cmpl,left_ev,1,right_ev,1,work,12,info)
!         print *, "eigen-values: ", ev_real
! 
!         print *, "transcorrelated Hamiltonian: "
!         t_trans_corr_2body = .true.
!         do i = 1, 4 
!             do j = 1, 4 
!                 hamil_trancorr(i,j) = get_helement_k_space_hub(hilbert_nI(:,i),hilbert_nI(:,j))
!             end do
!             print *, hamil_trancorr(i,:)
!         end do
! 
!         tmp_hamil = hamil_trancorr
!         call dgeev('N','N',4,tmp_hamil,4,ev_real,ev_cmpl,left_ev,1,right_ev,1,work,12,info)
!         print *, "eigen-values: ", ev_real
! 
! !         call stop_all("here","for now")

    end subroutine test_4e_ms1

    function matrix_exponential(matrix) result(exp_matrix)
        ! calculate the matrix exponential of a real, symmetric 2-D matrix with lapack 
        ! routines
        ! i need A = UDU^-1
        ! e^A = Ue^DU^-1
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: exp_matrix(size(matrix,1),size(matrix,2))

        ! maybe i need to allocate this stuff:
        real(dp) :: vectors(size(matrix,1),size(matrix,2)), values(size(matrix,1))
        real(dp) :: work(3*size(matrix,1)-1), inverse(size(matrix,1),size(matrix,2))
        real(dp) :: exp_diag(size(matrix,1),size(matrix,2))
        integer :: info, n

        n = size(matrix,1)

        ! first i need to diagonalise the matrix and calculate the 
        ! eigenvectors 
        vectors = matrix

        call dsyev('V', 'U', n, vectors, n, values, work, 3*n-1,info)

        ! now i have the eigenvectors, which i need the inverse of 
        ! it is rotation only or? so i would just need a transpose or?
!         inverse = matrix_inverse(vectors) 
        inverse = transpose(vectors)

        ! i need to construct exp(eigenvalues) as a diagonal matrix! 
        exp_diag = matrix_diag(exp(values))

        exp_matrix = matmul(matmul(vectors,exp_diag),inverse)
        
    end function matrix_exponential

    function matrix_diag(vector) result(diag) 
        ! constructs a diagonal matrix with the vector on the diagonal 
        real(dp), intent(in) :: vector(:)
        real(dp) :: diag(size(vector),size(vector))

        integer :: i 

        diag = 0.0_dp

        do i = 1, size(vector)
            diag(i,i) = vector(i)
        end do

    end function matrix_diag

    function matrix_inverse(matrix) result(inverse)
        ! from fortran-wiki! search there for "matrix+inversion"
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: inverse(size(matrix,1),size(matrix,2))
        character(*), parameter :: this_routine = "matrix_inverse"

        real(dp) :: work(size(matrix,1))
        integer :: ipiv(size(matrix,1))
        integer :: n, info

        inverse = matrix 
        n = size(matrix,1)

        call dgetrf(n,n,inverse,n,ipiv,info)

        if (info /= 0) call stop_all(this_routine, "matrix singular!")

        call dgetri(n, inverse, n, ipiv, work, n, info)

        if (info /= 0) call stop_all(this_routine, "matrix inversion failed!")

    end function matrix_inverse

    subroutine setup_all(ptr)
        class(lattice), intent(in) :: ptr


        nBasisMax = 0
        nullify(G1)
        nullify(tmat2d)
        deallocate(kPointToBasisFn)

        call setup_nbasismax(ptr)
        call setup_g1(ptr)
        call setup_tmat_k_space(ptr)
        call setup_kPointToBasisFn(ptr)
!         call setup_k_space_hub_sym(ptr)

        omega = real(ptr%get_nsites(),dp)
        nBasis = 2*ptr%get_nsites()

!         allocate(umat(1))
        uhub = 1.0
        bhub = -1.0
        umat = h_cast(real(uhub,dp)/real(omega,dp))

        trans_corr_param_2body = 0.1
        three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)

    end subroutine setup_all

    subroutine test_3e_ms1

        integer :: hilbert_nI(3,3), i, j, work(9), info
        HElement_t(dp) :: hamil(3,3), hamil_trancorr(3,3), tmp_hamil(3,3)
        real(dp) :: ev_real(3), ev_cmpl(3), left_ev(1,3), right_ev(1,3)
        real(dp) :: test(3,3), t_mat(3,3), trans_hamil(3,3)

        integer :: n_pairs

        nOccAlpha = 2
        nOccBeta = 1
        nel = 3

        ! i have to setup the whole system
        lat => lattice('chain',3,1,1,.true.,.true.,.true.,'k-space')

        call setup_all(lat)

        print *, "also test 3 electron system for consistency!" 

        hilbert_nI(:,1) = [2,3,4]
        hilbert_nI(:,2) = [4,5,6] 
        hilbert_nI(:,3) = [1,2,6]

        ! first create the non-transcorrelated Hamiltonian
        t_trans_corr_2body = .false. 
        print *, "un-transcorrelated Hamiltonian: "

        hamil = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil)

        ! the originial hamiltonian also gives me the transformation matrix 
        ! do it the stupid way 
        t_mat = get_tranformation_matrix(hamil, 2) 
        
        ! use the lapack routines to solve these quickly.. 
        print *, "eigen-values: ", calc_eigenvalues(hamil)

        print *, "transcorrelated Hamiltonian: "
        t_trans_corr_2body = .true.

        hamil_trancorr = create_hamiltonian(hilbert_nI)
        call print_matrix(hamil_trancorr)

        print *, "eigen-values: ", calc_eigenvalues(hamil_trancorr)

        trans_hamil = matmul(matmul(matrix_exponential(-t_mat),hamil),matrix_exponential(t_mat))

        print *, "transformed hamiltonian: " 
        call print_matrix(trans_hamil)

        print *, "eigen-values: ", calc_eigenvalues(trans_hamil)

!         ! i know now, where there is an error-- between those 2 matrix 
!         ! elements: 
!         print *, "nI:", hilbert_nI(:,1)
!         print *, "nJ:", hilbert_nI(:,2)
!         print *, "H_ij: ", hamil_trancorr(1,2), hamil_trancorr(2,1)
! 
!         ! check the individual contribs here! 
!         ! excitation: (2,3) -> (5,6)
!         ! ex: should have a + sign
!         ! 2 3 
!         ! 5 6 
!         ! 2 3 4 -> should have an overall + sign
!         print *, "Excitation (2,3) -> (5,6)"
!         print *, "two-body: ", two_body_transcorr_factor(G1(2)%k,G1(6)%k), & 
!                                two_body_transcorr_factor(G1(3)%k,G1(5)%k)
! 
!         print *, "three_body: ", three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(2)%k,G1(3)%k,G1(6)%k,1), &
!                                  three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(3)%k,G1(2)%k,G1(5)%k,-1)
! 
!         print *, "excitation: (3,4) -> (1,6): "
!         ! ex: should have a - sign!
!         ! 3 4
!         ! 1 6
!         ! 2 3 4 -> should have an overal - sign
!         print *, "two-body: ", two_body_transcorr_factor(G1(3)%k,G1(1)%k), & 
!                                two_body_transcorr_factor(G1(4)%k,G1(6)%k)
! 
!         print *, "three_body: ", three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(3)%k,G1(4)%k,G1(1)%k,-1), &
!                                  three_body_transcorr_fac(hilbert_nI(:,1), &
!                                     G1(4)%k,G1(3)%k,G1(6)%k,1)
! 
!         print *, "excitation: (4,5) -> (1,2): "
!         ! ex: should have a + sign
!         ! 4 5
!         ! 1 2
!         ! 4 5 6 -> should have an overall + sign
!         print *, "two-body: ", two_body_transcorr_factor(G1(4)%k,G1(2)%k), & 
!                                two_body_transcorr_factor(G1(5)%k,G1(1)%k)
! 
!         print *, "three_body: ", three_body_transcorr_fac(hilbert_nI(:,2), &
!                                     G1(4)%k,G1(5)%k,G1(2)%k,1), &
!                                  three_body_transcorr_fac(hilbert_nI(:,2), &
!                                     G1(5)%k,G1(4)%k,G1(1)%k,-1)


        ! both sign conventions agree here! maybe it has to do with this?!


        call stop_all("here","for now")
    end subroutine test_3e_ms1

    function get_tranformation_matrix(hamil, n_pairs) result(t_matrix)
        ! n_pairs is actually also a global system dependent quantitiy.. 
        ! which actually might be helpful.. but input it here! 
        HElement_t(dp), intent(in) :: hamil(:,:)
        integer, intent(in) :: n_pairs
        real(dp) :: t_matrix(size(hamil,1),size(hamil,2))

        integer :: i, j

        t_matrix = 0.0_dp 

        do i = 1, size(hamil,1)
            do j = 1, size(hamil,1)
                if (i == j) then 
                    t_matrix(i,i) = n_pairs 
                else 
                    if (abs(hamil(i,j)) > EPS) then 
                        t_matrix(i,j) = sign(1.0_dp, hamil(i,j))
                    end if
                end if
            end do
        end do

        t_matrix = trans_corr_param_2body/omega * t_matrix

    end function get_tranformation_matrix


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

    subroutine setup_kPointToBasisFn_test

        use SymExcitDataMod, only: kPointToBasisFn
        class(lattice), pointer :: ptr
        integer :: i
        print *, "" 
        print *, "testing: setup_kPointToBasisFn" 
        
        ptr => lattice('chain', 4, 1, 1, .true., .true., .true.,'k-space')

        call setup_kPointToBasisFn(ptr)

        call assert_equals(1, kPointToBasisFn(-1,0,0,1))
        call assert_equals(2, kPointToBasisFn(-1,0,0,2))
        call assert_equals(3, kPointToBasisFn(0,0,0,1))
        call assert_equals(4, kPointToBasisFn(0,0,0,2))
        call assert_equals(5, kPointToBasisFn(1,0,0,1))
        call assert_equals(6, kPointToBasisFn(1,0,0,2))
        call assert_equals(7, kPointToBasisFn(2,0,0,1))

        deallocate(kPointToBasisFn)

    end subroutine setup_kPointToBasisFn_test

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

        umat = uhub/omega

        ! todo: tests for actual transcorrelation! 
        call assert_equals(h_cast(1.0_dp), get_diag_helement_k_sp_hub([1,2,3,4]))

    end subroutine get_diag_helement_k_sp_hub_test

    subroutine get_offdiag_helement_k_sp_hub_test

        integer :: ex(2,2)
        integer, allocatable :: nI(:)

        nel = 2 
        allocate(nI(nel)); nI = [3,4] 

        print *, ""
        print *, "testing: get_offdiag_helement_k_sp_hub" 
        t_trans_corr_2body = .false. 

        ! 0 due to spin-symmetry:
        ex(1,:) = [1,3]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ex(1,:) = [1,2]
        ex(2,:) = [6,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! 0 due to momentum symmetry: 
        ex(2,:) = [3,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! this should contribute:
        ex(1,:) = [1,6]
        ex(2,:) = [3,4]
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ex(1,:) = [6,1]
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(2,:) = [1,6]
        ex(1,:) = [3,4]
        call assert_equals(h_cast(-uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.true.))

        t_trans_corr = .true. 
        trans_corr_param = 2.0_dp 

        call assert_equals(h_cast(uhub/real(omega,dp)*exp(4.0_dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(1,:) = [1,6]
        ex(2,:) = [3,4]

        call assert_equals(h_cast(uhub/real(omega,dp)*exp(-4.0_dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ex(1,:) = [6,1]
        call assert_equals(h_cast(uhub/real(omega,dp)*exp(-4.0_dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        t_trans_corr = .false. 
        trans_corr_param = 0.0_dp 

        t_trans_corr_2body = .true. 
        trans_corr_param_2body = 1.0_dp

        ex(1,:) = [1,3]
        ex(2,:) = [2,4] 
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(1,:) = [1,2]
        ex(2,:) = [3,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! this should contribute:
        ex(1,:) = [1,6]
        ex(2,:) = [3,4]
        ! and have a triples contribution now!
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ! this should be independent of order of electrons
        ex(1,:) = [6,1]
        call assert_equals(h_cast(uhub/real(omega,dp)), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

!         three_body_prefac = 1.0_dp
        ex(1,:) = [1,6]
        ex(2,:) = [3,4]
        print *, "---------------------------------" 
        print *, "test order for matrix elements triples contrib to 'normal'"
        print *, "(1,6) -> (3,4): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        ex(1,:) = [6,1]
        print *, "(6,1) -> (3,4): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        ex(2,:) = [4,3]
        print *, "(6,1) -> (4,3): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        ex(1,:) = [1,6]
        print *, "(1,6) -> (4,3): ", get_offdiag_helement_k_sp_hub(nI,ex,.false.)
        print *, "---------------------------------" 

        ! 0 due to momentum symmetry
        ex(1,:) = [1,5]
        ex(2,:) = [3,7]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ex(1,:) = [2,6]
        ex(2,:) = [4,8]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! this should contribute: 
        ex(1,:) = [1,3]
        ex(2,:) = [5,7] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! the order in EX should not matter.. figure that out! 
        ex(1,:) = [3,1]
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        print *," --------------------------" 
        print *, "test order of matrix elements for parallel excitations: " 
        ex(1,:) = [1,3]
        ex(2,:) = [5,7] 
        print *, "(1,3) -> (5,7): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        ex(1,:) = [3,1]
        print *, "(3,1) -> (5,7): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        ex(2,:) = [7,5]
        print *, "(3,1) -> (7,5): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        ex(1,:) = [1,3]
        print *, "(1,3) -> (7,5): ", get_offdiag_helement_k_sp_hub(nI, ex,.false.)
        print *," --------------------------" 

        
        ex(1,:) = [2,4]
        ex(2,:) = [6,8] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))

        ! order of orbitals should also not matter! 
        ex(2,:) = [8,6] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.))
        ! and this should have opposite sign.... 
        ! why should this have the opposite sign?  the spin should not matter! 
        ex(2,:) = [1,3]
        ex(1,:) = [5,7] 
        call assert_equals(h_cast(-4.0_dp*three_body_prefac), get_offdiag_helement_k_sp_hub(nI,ex,.false.),1.0e-12)

        nel = 4 
        nI = [3,4,6,7]
        ex(1,:) = [4,6]
        ex(2,:) = [2,8]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI, ex,.false.))

        ex(1,:) = [3,7] ! k = 2 
        ex(2,:) = [1,5] ! k = 0 
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_k_sp_hub(nI, ex,.false.))

        t_trans_corr_2body = .false. 

    end subroutine get_offdiag_helement_k_sp_hub_test

    subroutine get_helement_k_space_hub_test

        integer, allocatable :: ni(:), nJ(:)
        integer :: ex(2,3), ic_ret

        nel = 2
        allocate(nI(nel)); allocate(nJ(nel)); 

        print *, ""
        print *, "testing: get_helement_k_space_hub_test" 
        nI = [1,2]
        nJ = [3,4] 

        ic_ret = -1
        call assert_equals(h_cast(0.0_dp), get_helement_k_space_hub(nI,nJ,ic_ret))
        call assert_equals(2, ic_ret)

        call assert_equals(h_cast(uhub/omega), get_helement_k_space_hub([1,6],[3,4]))

        nel = 4
        deallocate(nI); deallocate(nJ); allocate(nJ(nel)); allocate(nI(nel)); 
        nI = [3,6,7,8]
        nJ = [1,2,5,8]
        ic_ret = -1
        call assert_equals(h_cast(0.0_dp), get_helement_k_space_hub(nI,nJ,ic_ret))
        call assert_equals(3, ic_ret) 

        t_trans_corr_2body = .true. 
        call assert_equals(h_cast(-4.0*three_body_prefac), get_helement_k_space_hub(nI,nJ))
        ! todo: more tests! 

        t_trans_corr_2body = .false. 

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

        nel = 4 
        nOccBeta = 2
        nOccAlpha = 2 
        deallocate(nI); allocate(nI(nel)); nI = [1,2,3,4]

        call pick_spin_opp_elecs(nI, elecs, p_elec) 
        call assert_equals(0.25_dp, p_elec) 
        call assert_true(.not. same_spin(nI(elecs(1)),nI(elecs(2))))

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

        integer, allocatable :: nI(:), orb_list(:,:)
        integer(n_int), allocatable :: ilutI(:) 
        real(dp), allocatable :: cum_arr(:) 
        real(dp) :: cum_sum, cpt 
        integer :: tgt 

        nel = 4 
        allocate(nI(nel)); allocate(ilutI(0:niftot)); allocate(orb_list(8,2)); 
        allocate(cum_arr(8))

        print *, "" 
        print *, "testing: create_ab_list_hubbard"
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        call create_ab_list_hubbard(nI, ilutI,[1,2], orb_list, cum_arr, cum_sum)

        call assert_equals(0.25_dp, cum_sum) 
        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum)
        call assert_equals(0.25_dp, cum_sum) 

        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum, 1, cpt)
        call assert_equals(0.25_dp, cum_sum) 
        call assert_equals(0.0_dp, cpt) 

        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum, 7, cpt)
        call assert_equals(0.5_dp, cpt)
        call create_ab_list_hubbard(nI, ilutI,[3,4], orb_list, cum_arr, cum_sum, 8, cpt)
        call assert_equals(0.5_dp, cpt)


        ! todo: i also have to do that for 2-body-transcorrelation, which 
        ! leads to parallel double excitations in the k-space hubbard 
        ! case -> check here if the get_orb_from_kpoints() function, works 
        ! correctly for ispn /= 2 and thub!

        ! todo: more tests! 

    end subroutine create_ab_list_hubbard_test

    subroutine calc_pgen_k_space_hubbard_test

        integer :: nI(4), ex(2,2)
        integer(n_int) :: ilutI(0:0)

        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard"

        ni = [1,2,3,4] 
        call EncodeBitDet(nI, ilutI) 

        ex(1,:) = [1,2]
        ex(2,:) = [5,6]

        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 3))
        call assert_equals(0.25_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 2))

        ex(1,:) = [3,4]
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 2))

        ex(2,:) = [7,8] 
        call assert_equals(0.25_dp, calc_pgen_k_space_hubbard(nI, ilutI, ex, 2))

    end subroutine calc_pgen_k_space_hubbard_test

    subroutine gen_excit_k_space_hub_test

        integer :: nI(4), ex(2,2), nJ(4)
        integer(n_int) :: ilutI(0:0), ilutJ(0:0)
        HElement_t(dp) :: hel
        real(dp) :: pgen 
        type(excit_gen_store_type) :: store
        logical :: tpar, t_found(6), found_all
        integer :: ic

        real(dp) :: p_elec = 0.25_dp

        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        print *, ""
        print *, "testing: gen_excit_k_space_hub" 

        t_found = .false. 
        found_all = .false. 

        do while (.not. found_all) 
            call gen_excit_k_space_hub(nI, ilutI, nJ, ilutJ, 0, ic, ex, & 
                    tpar, pgen, hel, store) 

            if (all(nJ == [3,4,5,6]) .and. .not. t_found(1)) then 
                t_found(1) = .true.
                ! do asserts: 
                call assert_equals(p_elec, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [2,3,5,8]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [2,3,6,7]) .and. .not. t_found(3)) then 
                t_found(3) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [1,4,5,8]) .and. .not. t_found(4)) then 
                t_found(4) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if

            if (all(nJ == [1,4,6,7]) .and. .not. t_found(5)) then 
                t_found(5) = .true. 
                call assert_equals(p_elec/2.0, pgen)
                call assert_true(.not. tpar)
            end if
            
            if (all(nJ == [1,2,7,8]) .and. .not. t_found(6))  then 
                t_found(6) = .true. 
                call assert_equals(p_elec, pgen)
                call assert_true(.not. tpar)
            end if

            found_all = all(t_found)

        end do

    end subroutine gen_excit_k_space_hub_test

    subroutine gen_parallel_double_hubbard_test

        integer :: nI(4), nJ(4), ex(2,2) 
        integer(n_int) :: ilutI(0:0), ilutJ(0:0) 
        real(dp) :: pgen 
        logical :: tpar, found_all, t_found(2)

        print *, "" 
        print *, "testing: gen_parallel_double_hubbard "

        t_trans_corr_2body = .true. 
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        found_all = .false. 
        t_found = .false. 

        do while (.not. found_all)
            call gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tPar, pgen) 

            if (all(nJ == [2,4,5,7]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals(0.5_dp, pgen)
                call assert_true(tpar)
            end if

            if (all(nJ == [1,3,6,8]) .and. .not. t_found(2)) then 
                t_found(2) = .true. 
                call assert_equals(0.5_dp, pgen)
                call assert_true(tpar)
            end if 

            found_all = all(t_found)
        end do

        t_trans_corr_2body = .false. 

    end subroutine gen_parallel_double_hubbard_test

    subroutine gen_triple_hubbard_test

        integer :: nI(4), nJ(4), ex(2,3) 
        integer(n_int) :: ilutI(0:0), ilutJ(0:0) 
        real(dp) :: pgen 
        logical :: tpar, found_all, t_found(2)

        print *, "" 
        print *, "testing: gen_triple_hubbard "

        found_all = .false. 
        t_found = .false. 
        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true. 

        do while (.not. found_all) 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tpar, pgen) 

            if (all(nJ == [1,2,4,5]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals(0.125_dp, pgen)
            end if
            found_all = t_found(1)

        end do

        nI = [3,4,5,8] 
        call EncodeBitDet(nI, ilutI) 

        found_all = .false. 
        t_found = .false. 

        do while (.not. found_all) 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tpar, pgen) 

            if (all(nJ == [1,2,3,6]) .and. .not. t_found(1)) then 
                t_found(1) = .true. 
                call assert_equals(0.125_dp, pgen)
            end if
            found_all = t_found(1)

        end do

        t_trans_corr_2body = .false.

    end subroutine gen_triple_hubbard_test

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
        call assert_equals(1.0_dp/6.0_dp, p_elec) 
        call assert_equals(1, sum_ms) 
        call assert_true(any(elecs == 1))

        nOccBeta = 4 
        nOccAlpha = 1
        call pick_three_opp_elecs([1,3,5,7,8], elecs, p_elec, sum_ms) 
        call assert_equals(1.0_dp/6.0_dp, p_elec) 
        call assert_equals(-1, sum_ms) 
        call assert_true(any(elecs == 5))

        nel = 5
        nOccBeta = 3
        nOccAlpha = 2 
        call pick_three_opp_elecs([1,2,3,4,5], elecs, p_elec, sum_ms) 
        if (sum_ms == 1) then 
            call assert_equals(1.0_dp/10.0_dp, p_elec)
        else 
            call assert_equals(7.0_dp/60.0_dp, p_elec,1.0e-12)
        end if

        call pick_three_opp_elecs([1,2,3,4,5], elecs, p_elec, sum_ms) 
        if (sum_ms == 1) then 
            call assert_equals(1.0_dp/10.0_dp, p_elec)
        else 
            call assert_equals(7.0_dp/60.0_dp, p_elec, 1.0e-12)
        end if

        nel = 4
        nOccBeta = 2
        nOccAlpha = 2 
        call pick_three_opp_elecs([1,2,3,4], elecs, p_elec) 
        call assert_equals(0.25_dp, p_elec)

    end subroutine pick_three_opp_elecs_test

    subroutine pick_spin_par_elecs_test

        integer :: elecs(2), ispn
        real(dp) :: p_elec
        integer :: nI(6)
        
        print *, ""
        print *, "testing: pick_spin_par_elecs"
        nel = 2 
        nOccBeta = 2
        nOccAlpha = 0
        call pick_spin_par_elecs([1,3],elecs,p_elec, ispn)
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(1, ispn) 
        call assert_equals(3, sum(elecs))

        nOccAlpha = 2 
        nOccBeta = 0
        call pick_spin_par_elecs([2,4],elecs,p_elec, ispn)
        call assert_equals(1.0_dp, p_elec) 
        call assert_equals(3, ispn) 
        call assert_equals(3, sum(elecs))

        nel = 4 
        nOccBeta = 2
        call pick_spin_par_elecs([1,2,3,4], elecs, p_elec, ispn) 
        call assert_equals(0.5_dp, p_elec) 
        if (ispn == 1) then 
            call assert_equals(4, sum(elecs))
        else if (ispn == 3) then 
            call assert_equals(6, sum(elecs))
        end if
        
        nel = 6 
        nOccBeta = 3 
        nOccAlpha = 3

        call pick_spin_par_elecs([1,2,3,4,5,6], elecs, p_elec) 
        call assert_equals(1.0_dp/6.0_dp, p_elec) 
        nI = [1,2,3,4,5,6]
        call assert_true(same_spin(nI(elecs(1)),nI(elecs(2))))

        nel = 4
        nOccBeta = 2
        nOccAlpha = 2
    end subroutine pick_spin_par_elecs_test

    subroutine pick_a_orbital_hubbard_test

        integer(n_int) :: ilutI(0:0)
        integer :: orb, sum_ms
        real(dp) :: p_orb

        print *, ""
        print *, "testing: pick_a_orbital_hubbard "
        call EncodeBitDet([1,2,3,4], ilutI)

        call pick_a_orbital_hubbard(ilutI, orb, p_orb, -1) 
        call assert_true(orb == 6 .or. orb == 8) 
        call assert_equals(0.5_dp, p_orb)

        call pick_a_orbital_hubbard(ilutI, orb, p_orb, 1) 
        call assert_true(orb == 5 .or. orb == 7) 
        call assert_equals(0.5_dp, p_orb)

        call pick_a_orbital_hubbard(ilutI, orb, p_orb) 
        call assert_equals(0.25_dp, p_orb)

    end subroutine pick_a_orbital_hubbard_test

    subroutine pick_ab_orbitals_hubbard_test

        integer, allocatable :: nI(:)
        integer(n_int), allocatable :: ilutI(:) 
        integer :: orbs(2)
        real(dp) :: p_orb 

        allocate(nI(nel)); allocate(ilutI(0:niftot))

        print *, "" 
        print *, "testing: pick_ab_orbitals_hubbard"
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)

        call pick_ab_orbitals_hubbard(nI, ilutI, [1,2], orbs, p_orb)

        call assert_equals(1.0_dp, p_orb) 
        call assert_equals(11, sum(orbs))

        call pick_ab_orbitals_hubbard(nI, ilutI, [3,4], orbs, p_orb)
        call assert_equals(1.0_dp, p_orb) 
        call assert_equals(15, sum(orbs))

    end subroutine pick_ab_orbitals_hubbard_test

    subroutine pick_bc_orbitals_hubbard_test

        integer:: nI(4)
        integer(n_int) :: ilutI(0:0) 
        integer :: orbs(2)
        real(dp) :: p_orb 


        t_trans_corr_2body = .true. 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)

        print *, ""
        print *, "testing: pick_bc_orbitals_hubbard"
        call pick_bc_orbitals_hubbard(nI, ilutI,[3,6,7],2,orbs,p_orb)
        call assert_equals(6, sum(orbs))
        call assert_equals(1.0_dp, p_orb) 

        nI = [3,4,5,8]
        call EncodeBitDet(nI, ilutI)
        call pick_bc_orbitals_hubbard(nI, ilutI,[4,5,8],1,orbs,p_orb)
        call assert_equals(8, sum(orbs))
        call assert_equals(1.0_dp, p_orb) 

    end subroutine pick_bc_orbitals_hubbard_test

    subroutine create_ab_list_par_hubbard_test

        integer:: nI(4), orb_list(4,2), tgt
        integer(n_int) :: ilutI(0:0) 
        real(dp) :: cum_sum, cpt, cum_arr(4)

        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true.
        print *, ""
        print *, "testing: create_ab_list_par_hubbard"
        call create_ab_list_par_hubbard(nI, ilutI, [1,3], orb_list, cum_arr, cum_sum)
        call assert_true(cum_sum > 0.0_dp)

        call create_ab_list_par_hubbard(nI, ilutI, [1,3], orb_list, cum_arr, cum_sum, 5, cpt)
        call assert_equals(0.5_dp, cpt) 
        call create_ab_list_par_hubbard(nI, ilutI, [1,3], orb_list, cum_arr, cum_sum, 7, cpt)
        call assert_equals(0.5_dp, cpt) 

        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum)
        call assert_true(cum_sum > 0.0_dp)

        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum, 1, cpt)
        call assert_equals(0.0_dp, cpt)
        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum, 5, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_ab_list_par_hubbard(nI, ilutI, [2,4], orb_list, cum_arr, cum_sum, 6, cpt)
        call assert_equals(0.5_dp, cpt)

        nI = [1,2,4,5] 
        call EncodeBitDet(nI, ilutI) 
        call create_ab_list_par_hubbard(nI, ilutI, [1,5], orb_list, cum_arr, cum_sum)
        call assert_equals(0.0_dp, cum_sum)

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)
        call create_ab_list_par_hubbard(nI, ilutI, [4,6], orb_list, cum_arr, cum_sum)

        print *, "" 
        print *, "cum_sum: ", cum_sum
        print *, "cum_arr: ", cum_arr
        print *, "orb_list(:,1):", orb_list(:,1)
        print *, "orb_list(:,2):", orb_list(:,2)

        t_trans_corr_2body = .false.

    end subroutine create_ab_list_par_hubbard_test

    subroutine pick_ab_orbitals_par_hubbard_test

        integer:: nI(4), orbs(2)
        integer(n_int) :: ilutI(0:0) 
        real(dp) :: p_orb 

        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true.

        print *, ""
        print *, "testing: pick_ab_orbitals_par_hubbard"
        call pick_ab_orbitals_par_hubbard(nI, ilutI, [1,3], orbs, p_orb)
        call assert_equals(12, sum(orbs))
        call assert_equals(1.0_dp, p_orb)

        call pick_ab_orbitals_par_hubbard(nI, ilutI, [2,4], orbs, p_orb)
        call assert_equals(14, sum(orbs))
        call assert_equals(1.0_dp, p_orb)

        nI = [1,2,4,5] 
        call EncodeBitDet(nI, ilutI) 
        call pick_ab_orbitals_par_hubbard(nI, ilutI, [1,5], orbs, p_orb)
        call assert_equals(0.0_dp, p_orb)

        t_trans_corr_2body = .false.

    end subroutine pick_ab_orbitals_par_hubbard_test

    subroutine get_transferred_momenta_test

        integer :: ex(2,2), k_vec_a(3), k_vec_b(3)

        print *, ""
        print *, "testing: get_transferred_momenta"
        ex(1,:) = [1,2]
        ex(2,:) = [3,4] 
        call get_transferred_momenta(ex, k_vec_a, k_vec_b) 

        call assert_equals(1, k_vec_a(1)) 
        call assert_equals(-1, k_vec_b(1))

        ex(1,:) = [1,3] 
        ex(2,:) = [5,7]
        call get_transferred_momenta(ex, k_vec_a, k_vec_b) 

        call assert_equals(2, k_vec_a(1))
        call assert_equals(-1, k_vec_b(1))

    end subroutine get_transferred_momenta_test

    subroutine get_orb_from_kpoints_three_test

        print *, ""
        print *, "testing: get_orb_from_kpoints_three: "
        call assert_equals(3, get_orb_from_kpoints_three([1,2,3],1,2))
        call assert_equals(5, get_orb_from_kpoints_three([1,2,3],4,5))
        call assert_equals(5, get_orb_from_kpoints_three([1,3,5],1,3))
        call assert_equals(2, get_orb_from_kpoints_three([2,4,6],4,6))

        call assert_equals(7, get_orb_from_kpoints_three([1,2,3],7,8))

        call assert_equals(7, get_orb_from_kpoints_three([3,4,5],6,8))

    end subroutine get_orb_from_kpoints_three_test

    subroutine create_bc_list_hubbard_test

        integer :: nI(4), orb_list(4,2), tgt 
        integer(n_int) :: ilutI(0:0) 
        real(dp) :: cum_arr(4), cum_sum, cpt 

        t_trans_corr_2body = .true. 
        print *, ""
        print *, "testing: create_bc_list_hubbard"
        nI = [1,2,3,4]
        call EncodeBitDet(nI, ilutI)
        call create_bc_list_hubbard(nI, ilutI, [1,2,3],6,orb_list, cum_arr, cum_sum)
        call assert_equals(0.0_dp, cum_sum) 

        call create_bc_list_hubbard(nI, ilutI, [1,2,4],5,orb_list, cum_arr, cum_sum)
        call assert_equals(0.0_dp, cum_sum) 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)
        call create_bc_list_hubbard(nI, ilutI, [3,6,7],2,orb_list, cum_arr, cum_sum)

        call assert_true(cum_sum > 0.0_dp) 

        nI = [3,4,5,8]
        call EncodeBitDet(nI, ilutI)
        call create_bc_list_hubbard(nI, ilutI, [4,5,8],1,orb_list, cum_arr, cum_sum)
        call assert_true(cum_sum > 0.0_dp) 
        call create_bc_list_hubbard(nI, ilutI, [4,5,8],1,orb_list, cum_arr, cum_sum, 4, cpt)
        call assert_equals(0.0_dp, cpt) 
        call create_bc_list_hubbard(nI, ilutI, [4,5,8],1,orb_list, cum_arr, cum_sum, 2, cpt)
        call assert_equals(0.5_dp, cpt)

        t_trans_corr_2body = .false. 

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

        ! and the order of the involved electrons should not change the 
        ! matrix element! ... damn.. it does.. i need to have some 
        ! convention i think.. as it is in the spin-opposite excitations 
        ! in the "normal" method.. 
        print *, "--------------------------" 
        print *, "testing order influence on sign: "
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [3,7,6]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [6,3,7]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [6,7,3] 
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(1,:) = [7,6,3]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [1,5,2]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [2,1,5]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [5,1,2]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)
        ex(2,:) = [5,2,1]
        print *, "(",ex(1,:),") -> (",ex(2,:),"): ", get_3_body_helement_ks_hub(nI, ex, .false.)

        ! fix this sign incoherence above! 
        call assert_true(.false.)

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

        integer :: nI(4), ex(2,3)
        integer(n_int) :: ilutI(0:0) 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)

        pDoubles = 0.8 
        pParallel = 0.2 

        t_trans_corr_2body = .true. 
        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_transcorr"
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,4))
        
        ex(1,:) = [3,4,0] ! k = 0
        ex(2,:) = [2,5,0] ! k = 0
        ! this should contribute! 
        call assert_equals(0.8*0.8/4.0_dp, calc_pgen_k_space_hubbard_transcorr(nI,iluti,ex,2))

        ex(1,:) = [3,6,0] ! k = 1
        ex(2,:) = [1,8,0] ! k = 1
        call assert_equals(0.8*0.8/4.0_dp, calc_pgen_k_space_hubbard_transcorr(nI,iluti,ex,2))

        ex(1,:) = [3,7,0] ! k=2
        ex(2,:) = [1,5,0] ! k=0
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,2))

        ex(1,:) = [4,6,0] ! k = 1
        ex(2,:) = [2,8,0] ! k = 1
        ! should contribute! no, becuase diagonal part is 0!
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,2))

        nI = [1,4,6,7]
        call EncodeBitDet(nI, ilutI)
        call assert_equals(0.8*0.2/2.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,2))

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI)

        ! the triple should be: 
        ex(1,:) = [3,6,7] 
        ex(2,:) = [1,2,5] 
        call assert_equals(0.2/16.0_dp, calc_pgen_k_space_hubbard_transcorr(ni,iluti,ex,3),1.0e-12)

        t_trans_corr_2body = .false.

    end subroutine calc_pgen_k_space_hubbard_transcorr_test

    subroutine calc_pgen_k_space_hubbard_par_test

        integer :: nI(4), ex(2,2)
        integer(n_int) :: ilutI(0:0) 

        nI = [1,2,3,4] 
        call EncodeBitDet(nI, ilutI)

        t_trans_corr_2body = .true.
        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_par"
        ex(1,:) = [1,3]
        ex(2,:) = [5,7] 

        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,3))

        call assert_equals(0.5_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))
        ex(1,:) = [2,4]
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))
        ex(2,:) = [6,8]
        call assert_equals(0.5_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))
        ex(1,:) = [1,4]
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_par(nI,ilutI,ex,2))

        t_trans_corr_2body = .false.

    end subroutine calc_pgen_k_space_hubbard_par_test

    subroutine calc_pgen_k_space_hubbard_triples_test

        integer :: nI(4), ex(2,3)
        integer(n_int) :: ilutI(0:0) 

        nI = [3,4,6,7]
        call EncodeBitDet(nI, ilutI) 

        t_trans_corr_2body = .true.
        print *, "" 
        print *, "testing: calc_pgen_k_space_hubbard_triples"
        ex(1,:) = [3,6,7]
        ex(2,:) = [1,2,5]

        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,0))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,1))
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,2))

        call assert_equals(1.0_dp/(2*8.0_dp), calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,3))

        ex(1,:) = [4,5,8] 
        call assert_equals(0.0_dp, calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,3))
        ex(2,:) = [1,2,6]
        nI = [3,4,5,8] 
        call EncodeBitDet(nI, ilutI)
        call assert_equals(1.0_dp/(2*8.0_dp), calc_pgen_k_space_hubbard_triples(nI, ilutI, ex,3))

        t_trans_corr_2body = .false.

    end subroutine calc_pgen_k_space_hubbard_triples_test

    subroutine make_triple_test

        integer, allocatable :: nI(:), nJ(:)
        integer :: ex(2,3), ex2(2,3)
        logical :: tpar, tpar_2, tpar_3, tpar_4
        integer(n_int) :: ilutI(0:nifd), ilutJ(0:nifd)

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

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        ! and now more complicated stuff:
        nI = [1,2,3]
        call make_triple(nI,nJ,[1,3,2],[7,5,4],ex,tpar) 
        call assert_equals([4,5,7], nJ, 3)
        call assert_equals([1,2,3], ex(1,:),3)
        call assert_equals([4,5,7], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        nI = [1,2,5]
        call make_triple(nI,nJ,[3,1,2],[3,4,7],ex,tpar) 
        call assert_equals([3,4,7], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([3,4,7], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        call make_triple(nI,nJ,[3,2,1],[8,7,3],ex,tpar) 
        call assert_equals([3,7,8], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([3,7,8], ex(2,:),3)
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        call make_triple(nI,nJ,[3,2,1],[4,7,9],ex,tpar) 
        call assert_equals([4,7,9], nJ, 3)
        call assert_equals([1,2,5], ex(1,:),3)
        call assert_equals([4,7,9], ex(2,:),3)
        call assert_true(.not.tpar)
        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        call make_triple(nI, nJ, [1,2,3], [3,4,7], ex, tpar) 
        call assert_true(.not.tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)

        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)
        nI = [1,2,5]

        nel = 4 

        deallocate(nJ); allocate(NJ(nel))
        deallocate(nI); allocate(nI(nel))
        
        nI = [1,2,5,7]
        call make_triple(nI,nJ,[1,2,3],[3,6,9],ex,tpar)
        call assert_true(tpar)

        ex2(1,1) = 3 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 3 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2, 3)


        ex(1,:) = [1,2,3]
        call FindExcitDet(ex, nI, 3, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        nel = -1

    end subroutine make_triple_test

    subroutine make_double_test

        use get_excit, only: make_double
        integer, allocatable :: nJ(:),ni(:)
        integer :: ex(2,2), ex2(2,2)
        logical :: tpar, tpar_2, tpar_3, tpar_4
        integer(n_int) :: ilutI(0:nifd), ilutJ(0:nifd)

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

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        call assert_true(tpar .eqv. tpar_3)

        ni = [1,2]

        call make_double([1,2],nJ, 1,2, 5,4, ex, tpar)
        call assert_true(.not. tpar)
        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2],nJ, 1,2, 3,6, ex, tpar) 
        call assert_true(.not. tpar)
        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        ni = [1,2]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_3)

        nel = 3 
        deallocate(nJ); allocate(nJ(nel))
        nI = [1,2,4]

        call make_double([1,2,4],nJ,1,2,5,6,ex,tpar)
        call assert_equals([4,5,6], nJ, 3)
        call assert_true(.not. tpar)

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,4]

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2,4],nj,1,2,3,6,ex,tpar)
        call assert_true(tpar)

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)
        nI = [1,2,4]

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2,4],nJ, 1, 2, 6, 7, ex, tpar)
        call assert_true(.not. tpar)
        nI = [1,2,4]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call make_double([1,2,3],nJ, 1, 2, 4, 7, ex, tpar)
        call assert_true(.not. tpar)
        nI = [1,2,3]

        ex2(1,1) = 2 
        call GetExcitation(nI,nJ,nel,ex2,tpar_4)
        call assert_equals(ex, ex2, 2,2)
        call assert_true(tpar .eqv. tpar_4)

        call EncodeBitDet(nI, ilutI)
        call EncodeBitDet(nJ, ilutJ)
        ex2(1,1) = 2 
        call GetBitExcitation(ilutI, ilutJ, ex2, tpar_3)
        call assert_equals(ex, ex2, 2,2)

        call FindExcitDet(ex, nI, 2, tpar_2)
        call assert_true(tpar .eqv. tpar_2)

        nel = -1


    end subroutine make_double_test

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

    subroutine gen_excit_k_space_hub_test_stochastic

        use GenRandSymExcitNUMod, only: TestGenRandSymExcitNU
        integer :: nI(nel)

        print *, "" 
        print *, "try to implement a stochastic version to check the" 
        print *, "calculated pgen and the actual one.. " 
        ! i should write a general test-runner for this.. which takes 
        ! an excitation generator as an input.. that would be nice! 
        ! something like: " 
        ni = [1,2,3,4]
        
        call setup_k_total(nI) 

!         call TestGenRandSymExcitNU(nI, 1000, 1.0, 2)
        call run_excit_gen_tester(gen_excit_k_space_hub, "gen_excit_k_space_hub", & 
            gen_all_excits=gen_all_excits_k_space_hubbard) 

    end subroutine gen_excit_k_space_hub_test_stochastic

    subroutine run_excit_gen_tester(excit_gen, excit_gen_name, opt_nI, opt_n_iters, & 
            gen_all_excits) 
        use procedure_pointers, only: generate_excitation_t, generate_all_excits_t
        use util_mod, only: binary_search

        procedure(generate_excitation_t) :: excit_gen 
        integer, intent(in), optional :: opt_nI(nel), opt_n_iters
        procedure(generate_all_excits_t), optional :: gen_all_excits
        character(*), intent(in) :: excit_gen_name
        character(*), parameter :: this_routine = "run_excit_gen_tester"

        integer :: i, nI(nel), n_iters 
        integer :: default_n_iters = 100000
        integer :: default_n_dets = 10 

        integer(n_int) :: ilut(0:niftot), tgt_ilut(0:niftot)
        integer :: nJ(nel), n_excits, ex(2,2), ic, ex_flag, i_unused = 0
        type(excit_gen_store_type) :: store 
        logical :: tPar, found_all
        real(dp) :: pgen, contrib
        HElement_t(dp) :: hel 
        integer(n_int), allocatable :: det_list(:,:)
        real(dp), allocatable :: contrib_list(:)
        logical, allocatable :: generated_list(:) 
        integer :: n_dets, n_generated, pos

        ASSERT(nel > 0)
        ! and also nbasis and stuff.. 
        ASSERT(nbasis > 0) 
        ASSERT(nel <= nbasis) 

        ! use some default values if not provided: 
        ! nel must be set! 
        if (present(opt_nI)) then 
            nI = opt_nI
        else 
            ! use HF-det as default
            nI = [(i, i = 1, nel)]
        end if

        if (present(opt_n_iters)) then 
            n_iters = opt_n_iters
        else 
            n_iters = default_n_iters
        end if

        ! i have to rewrite this routine, to a part which 
        ! creates all excitations and another who runs on possibly 
        ! multiple excitations! 
!         if (present(opt_n_dets)) then 
!             n_dets = opt_n_dets
!         else 
!             n_dets = default_n_dets
!         end if

        ! the problem here is now.. we want to calulate all the possible 
        ! excitations.. i would need a general routine, which does that 
        ! with only the hamiltonian knowledge.. this is not really there.. 
        ! i should have a routine: 
        ! for this special setup, which is tested.. 

        ! for some special systems we should provide a routine to 
        ! calculate all the excitations (hubbard, UEG eg.) 
        if (present(gen_all_excits)) then 
            call gen_all_excits(nI, n_excits, det_list)
        else 
            call gen_all_excits_default(nI, n_excits, det_list)
        end if

        print *, "total possible excitations: ", n_excits
        do i = 1, n_excits 
            call writebitdet(6, det_list(:,i),.true.)
        end do

        ! call this below now for the number of specified determinants 
        ! (also use excitations of the first inputted, to be really 
        !   consistent) 

        n_dets = min(n_dets, n_excits) 

        print *, "---------------------------------"
        print *, "testing: ", excit_gen_name
        print *, "for ", n_dets, " determinants" 
        print *, " and ", n_iters, " iterations "

        call EncodeBitDet(nI, ilut) 

        print *, "for starting determinant: ", nI 

        ! Lists to keep track of things
        allocate(generated_list(n_excits))
        allocate(contrib_list(n_excits))
        generated_list = .false.
        contrib_list = 0

        n_generated = 0
        contrib = 0.0_dp

        do i = 1, n_iters
            if (mod(i, 1000) == 0) then 
                print *, i, "/" ,n_iters, " - ", contrib / real(n_excits*i,dp) 
            end if
            call excit_gen(nI, ilut, nJ, tgt_ilut, ex_flag, ic, ex, tpar, pgen, & 
                        hel, store) 

            if (nJ(1) == 0) cycle 
            call EncodeBitDet(nJ, tgt_ilut) 
            pos = binary_search(det_list, tgt_ilut, nifd+1)
            if (pos < 0) then 
                print *, "nJ: ", nJ 
                print *, "ilutJ:", tgt_ilut
                call stop_all(this_routine, 'Unexpected determinant generated')
            else 
                generated_list(pos) = .true. 
                n_generated = n_generated + 1 

                contrib = contrib + 1.0_dp / pgen 
                contrib_list(pos) = contrib_list(pos) + 1.0_dp / pgen 
            end if 
        end do

        print *, n_generated, " dets generated in ", n_iters, " iterations " 
        print *, 100.0_dp * (n_iters - n_generated) / real(n_iters,dp), "% abortion rate" 
        print *, "Averaged contribution: ", contrib / real(n_excits*n_iters,dp)

        ! check all dets are generated: 
        call assert_true(all(generated_list))

        print *, "=================================="
        print *, "Contribution List: "
        do i = 1, n_excits 
            call writebitdet(6, det_list(:,i), .false.)
            print *, contrib_list(i)/real(n_iters,dp)
        end do
        ! and check the uniformity of the excitation generation
        call assert_true(all(abs(contrib_list / n_iters - 1.0_dp) < 0.01_dp))

    end subroutine run_excit_gen_tester

    subroutine gen_all_excits_k_space_hubbard(nI, n_excits, det_list) 

        use SystemData, only: tNoBrillouin, tUseBrillouin
        use neci_intfce, only: GenSymExcitIt2
        use GenRandSymExcitNUMod, only: IsMomentumAllowed

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits 
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        character(*), parameter :: this_routine = "gen_all_excits_k_space_hubbard"

        logical :: brillouin_tmp(2), tpar
        integer :: iMaxExcit, nStore(6),nExcitMemLen(1), excitcount, ex(2,2)
        integer :: nJ(nel), ierr, exFlag, iExcit
        integer(n_int) :: iLutnJ(0:niftot)
        integer, allocatable :: EXCITGEN(:)
        integer(n_int), allocatable :: triple_dets(:,:), temp_dets(:,:)
        integer :: n_triples, save_excits


!         exFlag = 2
!          
!         ! resuse the really weirdly implemented stuff in symrandexcit2... 
! 
!         ! The old excitation generator will not generate singles from the HF
!         ! unless tNoBrillouin is set
!         brillouin_tmp(1) = tNoBrillouin
!         brillouin_tmp(2) = tUseBrillouin
!         tNoBrillouin = .true.
!         tUseBrillouin = .false.
! 
!         !Find the number of symmetry allowed excitations there should be by looking at the full excitation generator.
!         !Setup excit generators for this determinant
!         iMaxExcit=0
!         nStore(1:6)=0
!         CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
!         ALLOCATE(EXCITGEN(nExcitMemLen(1)),stat=ierr)
!         IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating excitation generator")
!         EXCITGEN(:)=0
!         CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
!     !    CALL GetSymExcitCount(EXCITGEN,DetConn)
!         excitcount=0
! 
!     lp2: do while(.true.)
!             CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
!             IF(nJ(1).eq.0) exit lp2
!             IF (IsMomentumAllowed(nJ)) THEN
!                 ex(1,1) = 2
!                 call GetExcitation(nI,nJ,nel,ex,tpar)
!                 ! also test for the spin in the hubbard model! 
!                 if (t_trans_corr_2body .or. .not. same_spin(ex(1,1),ex(1,2))) then 
!                     excitcount=excitcount+1
!                     CALL EncodeBitDet(nJ,iLutnJ)
!                 end if
!             ENDIF
!         end do lp2
! 
!         ! and now do it again and fill the determinants.. 
!         allocate(det_list(0:niftot,excitcount)) 
!         EXCITGEN = 0
!         CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
!         excitcount = 0
! 
!     lp3: do while(.true.)
!             CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
!             IF(nJ(1).eq.0) exit lp3
!             IF (IsMomentumAllowed(nJ)) THEN
!                 ex(1,1) = 2
!                 call GetExcitation(nI,nJ,nel,ex,tpar)
!                 ! also test for the spin in the hubbard model! 
!                 if (t_trans_corr_2body .or. .not. same_spin(ex(1,1),ex(1,2))) then 
!                     excitcount=excitcount+1
!                     CALL EncodeBitDet(nJ,iLutnJ)
!                     det_list(:,excitcount) = ilutnJ
!                 end if
!             ENDIF
!         end do lp3
! 
!         n_excits = excitcount

        call gen_all_doubles_k_space(nI, n_excits, det_list)

        if (t_trans_corr_2body) then 
            save_excits = n_excits
            ! also account for triple excitations
            call gen_all_triples_k_space(nI, n_triples, triple_dets)

            n_excits = n_excits + n_triples

            allocate(temp_dets(0:niftot, save_excits), source = det_list(:,1:save_excits))

            deallocate(det_list) 

            allocate(det_list(0:niftot,n_excits)) 

            det_list(:,1:save_excits) = temp_dets 

            det_list(:,save_excits+1:n_excits) = triple_dets

        end if

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_excits_k_space_hubbard

    subroutine gen_all_excits_default(nI, n_excits, det_list) 
        use SymExcit3, only: CountExcitations3, GenExcitations3
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits 
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        character(*), parameter :: this_routine = "gen_all_excits_default"

        integer :: n_singles, n_doubles, n_dets, ex(2,2), ex_flag
        integer :: nJ(nel) 
        logical :: tpar, found_all 
        integer(n_int) :: ilut(0:niftot)

        n_excits = -1

        call EncodeBitDet(nI, ilut) 

        ! for reference in the "normal" case it looks like that: 
        call CountExcitations3(nI, 2, n_singles, n_doubles) 

        n_excits = n_singles + n_doubles 

        print *, "n_singles: ", n_singles
        print *, "n_doubles: ", n_doubles
        
        allocate(det_list(0:niftot,n_excits)) 
        n_dets = 0
        found_all = .false. 
        ex = 0
        ex_flag = 2 
        call GenExcitations3 (nI, ilut, nJ, ex_flag, ex, tpar, found_all, &
                              .false.)

        do while (.not. found_all)
            n_dets = n_dets + 1
            call EncodeBitDet (nJ, det_list(:,n_dets))

            call GenExcitations3 (nI, ilut, nJ, ex_flag, ex, tpar, &
                                  found_all, .false.)
        end do

        if (n_dets /= n_excits) then
            print *, "expected number of excitations: ", n_excits
            print *, "actual calculated ones: ", n_dets
            call stop_all(this_routine,"Incorrect number of excitations found")
        end if

        ! Sort the dets, so they are easy to find by binary searching
        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_excits_default

    subroutine gen_excit_k_space_hub_transcorr_test_stoch

        use bit_reps, only: decode_bit_det
        integer :: nI(nel), n_excits, i, nJ(nel), n_triples
        integer(n_int), allocatable :: det_list(:,:)

        print *, "" 
        print *, "testing: gen_excit_k_space_hub_transcorr" 
        print *, "first for a system with no possible triples, due to momentum conservation"

        pDoubles = 0.8 
        pParallel = 0.2 
        t_trans_corr_2body = .true. 

        nI = [1,2,3,4] 
        call setup_k_total(nI)
        call gen_all_triples_k_space(nI, n_excits, det_list)

        print *, "number of triple excitations: ", n_excits
        do i = 1, n_excits
            call writebitdet(6, det_list(:,i),.true.)
        end do

        call gen_all_excits_k_space_hubbard(nI, n_excits, det_list)

        ! for this momentum sector there are no, triple excitations valid.. 
        ! so test that for now! 
        call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
            "gen_excit_k_space_hub_transcorr",opt_ni = nI, & 
            gen_all_excits = gen_all_excits_k_space_hubbard) 

        do i = 1, n_excits 
            call decode_bit_det(nJ, det_list(:,i))
            call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
                "gen_excit_k_space_hub_transcorr",opt_ni = nJ, & 
                gen_all_excits = gen_all_excits_k_space_hubbard) 
        end do


        print *, "" 
        print *, "and now for a system with triples: "
        nI = [3,4,6,7]
        call setup_k_total(nI) 

        call gen_all_triples_k_space(nI, n_triples, det_list)

        print *, "number of triple excitations: ", n_triples
        do i = 1, n_triples
            call writebitdet(6, det_list(:,i),.true.)
        end do

        call gen_all_excits_k_space_hubbard(nI, n_excits, det_list)

        ! for this momentum sector there are no, triple excitations valid.. 
        ! so test that for now! 
        call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
            "gen_excit_k_space_hub_transcorr",opt_ni = nI, & 
            gen_all_excits = gen_all_excits_k_space_hubbard) 

        do i = 1, n_excits 
            call decode_bit_det(nJ, det_list(:,i))
            call run_excit_gen_tester(gen_excit_k_space_hub_transcorr_test, &
                "gen_excit_k_space_hub_transcorr",opt_ni = nJ, & 
                gen_all_excits = gen_all_excits_k_space_hubbard) 
        end do

    end subroutine gen_excit_k_space_hub_transcorr_test_stoch

end program test_k_space_hubbard
