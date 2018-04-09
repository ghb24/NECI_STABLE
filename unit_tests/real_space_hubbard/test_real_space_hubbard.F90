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
                          trans_corr_param, t_lattice_model, t_trans_corr_hop, brr, & 
                          t_trans_corr_2body, trans_corr_param_2body, &
                          t_trans_corr_new, t_uniform_excits, tHPHF

    use lattice_mod, only: lat, init_dispersion_rel_cache

    use sort_mod, only: sort
    
    use util_mod, only: get_free_unit

    use unit_test_helpers

    use dsfmt_interface, only: dsfmt_init 

    use fcimcdata, only: pSingles, pDoubles

    use lattice_models_utils, only: gen_all_excits_r_space_hubbard, &
                                    create_hilbert_space_realspace

    use HPHFRandexcitmod, only: gen_hphf_excit

    use k_space_hubbard, only: setup_k_space_hub_sym, setup_symmetry_table
    implicit none 

    integer :: failed_count 

    t_new_real_space_hubbard = .true.
    t_lattice_model = .true.

    call dsfmt_init(1)
    call init_fruit()
    ! run the test-driver 
    call exact_test()
    call real_space_hubbard_test_driver()
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine real_space_hubbard_test_driver() 
        ! this is the main function which calls all the other tests 
        
        ! or try running it with the provided runner of fruit: 
!         call run_test_case(gen_excit_rs_hubbard_hphf_test_stoch, "gen_excit_rs_hubbard_hphf_test_stoch")
!         call run_test_case(gen_excit_rs_hubbard_transcorr_hphf_test_stoch, "gen_excit_rs_hubbard_transcorr_hphf_test_stoch")
        call run_test_case(gen_excit_rs_hubbard_transcorr_uniform_hphf_test_stoch, "gen_excit_rs_hubbard_transcorr_uniform_hphf_test_stoch")
!         call run_test_case(gen_excit_rs_hubbard_test_stoch, "gen_excit_rs_hubbard_test_stoch")
!         call run_test_case(gen_excit_rs_hubbard_transcorr_test_stoch, "gen_excit_rs_hubbard_transcorr_test_stoch")
        call run_test_case(gen_excit_rs_hubbard_transcorr_uniform_test_stoch, "gen_excit_rs_hubbard_transcorr_uniform_test_stoch")
        call run_test_case(get_umat_el_hub_test, "get_umat_el_hub_test")
        call run_test_case(init_tmat_test, "init_tmat_test")
        call run_test_case(get_helement_test, "get_helement_test")
        call run_test_case(gen_excit_rs_hubbard_test, "gen_excit_rs_hubbard_test")
        call run_test_case(init_real_space_hubbard_test, "init_real_space_hubbard_test")
        call run_test_case(trans_corr_fac_test, "trans_corr_fac_test")
        call run_test_case(create_cum_list_rs_hubbard_test, "create_cum_list_rs_hubbard_test")
        call run_test_case(create_avail_neighbors_list_test, "create_avail_neighbors_list_test")
        call run_test_case(calc_pgen_rs_hubbard_test, "calc_pgen_rs_hubbard_test")
        call run_test_case(determine_optimal_time_step_test, "determine_optimal_time_step_test")
        call run_test_case(get_optimal_correlation_factor_test, "get_optimal_correlation_factor_test")
        call run_test_case(get_offdiag_helement_rs_hub_test, "get_offdiag_helement_rs_hub_test")

    end subroutine real_space_hubbard_test_driver

    subroutine exact_test

        use DetCalcData, only: nkry, nblk, b2l, ncycle
        use lanczos_wrapper, only: frsblk_wrapper

        integer :: i, n_eig, n_orbs, n_states
        integer, allocatable :: ni(:), hilbert_space(:,:)
        real(dp), allocatable :: e_values(:), e_vecs(:,:)
        integer(n_int), allocatable :: dummy(:,:)
        real(dp) :: j
        real(dp), allocatable :: j_vec(:)


        lat => lattice('tilted', 3, 3, 1,.true.,.true.,.true.)
        uhub = 8
        bhub = -1

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs

        call init_realspace_tests

        nel = 18
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,3,6,7,9,12,13,16,17,20,21,24,25,28,30,31,34,36]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call setup_arr_brr(lat)
        call init_hopping_transcorr()
!         call setup_symmetry_table()
!         call setup_k_space_hub_sym(lat) 
!         call init_dispersion_rel_cache()
!         call init_umat_rs_hub_transcorr()

        j_vec = linspace(-2.0,2.0,100)
        print *, "H diag: "
        t_trans_corr_hop = .true.
        do i = 1, size(j_vec) 
            trans_corr_param= j_vec(i)
            print *, J_vec(i), get_diag_helemen_rs_hub_transcorr_hop(nI)
        end do


        t_trans_corr_hop = .false.
        call stop_all("here","now")

        call create_hilbert_space_realspace(n_orbs, nOccAlpha, nOccBeta, & 
            n_states, hilbert_space, dummy)
        n_eig = 1

        allocate(e_values(n_eig))
        allocate(e_vecs(n_eig, size(hilbert_space,2)))

        print *, "size hilbert: ", size(hilbert_space, 2)
        nblk = 4
        nkry = 8 
        ncycle = 200
        b2l = 1.0e-13_dp

        print *, "nkry: ", nkry
        print *, "nblk: ", nblk
        print *, "b2l: ", b2l
        print *, "ncycle: ", ncycle

        ! try too big systems here: 
        call frsblk_wrapper(hilbert_space, size(hilbert_space, 2), n_eig, e_values, e_vecs)
        print *, "e_value lanczos:", e_values(1)

        j = 0.1_dp
        call exact_transcorrelation(lat, nI, j_vec, real(uhub,dp), hilbert_space)

!         call stop_all("here","now")

    end subroutine exact_test

    subroutine exact_transcorrelation(lat, nI, J, U, hilbert_space)
        class(lattice), intent(in) :: lat 
        integer, intent(in) :: nI(nel) 
        real(dp) :: J(:), U 
        integer, intent(in) :: hilbert_space(:,:)
        character(*), parameter :: this_routine = "exact_transcorrelation" 

        integer :: n_states, iunit, ind, i, k, l
        real(dp), allocatable :: e_values(:), e_vec(:,:), gs_vec(:)
        real(dp) :: gs_energy_orig, gs_energy, hf_coeff_hop(size(J))
        HElement_t(dp), allocatable :: hamil(:,:), hamil_hop(:,:), hamil_hop_neci(:,:), &
                                       diff(:,:)
        real(dp), allocatable :: t_mat(:,:)
        real(dp), allocatable :: neci_eval(:), e_vec_hop(:,:)
        character(30) :: filename, J_str

        ! initialize correctly for transcorrelation tests
        ! just do it for the hopping transcorrelation now! 
        n_states = size(hilbert_space,2)
        print *, "creating original hamiltonian: "

        t_trans_corr_2body = .false.
        t_trans_corr = .false.
        t_trans_corr_hop = .false.

        hamil = create_hamiltonian(hilbert_space)

        print *, "diagonalizing original hamiltonian: " 
        allocate(e_values(n_states));        e_values = 0.0_dp
        allocate(e_vec(n_states, n_states)); e_vec = 0.0_dp
        allocate(gs_vec(n_states));          gs_vec = 0.0_dp
        do i = 1, size(hamil,1)
            do k = 1, size(hamil,2)
                if (isnan(hamil(i,k)) .or. is_inf(hamil(i,k))) print *, i,k,hamil(i,k)
            end do
        end do
        call eig(hamil, e_values, e_vec) 

!         print *, "real-space hamiltonian: "
!         call print_matrix(hamil)

        ! find the ground-state
        ind = minloc(e_values,1) 
        gs_energy_orig = e_values(ind) 

        ! how do i need to access the vectors to get the energy? 
        ! eigenvectors are stored in the columns!!
        gs_vec = abs(e_vec(:,ind))

        call sort(gs_vec)

        ! and flip order.. 
        gs_vec = gs_vec(n_states:1:-1)

        ! and i think i want the sorted by maximum of the GS
        print *, "original ground-state energy: ", gs_energy_orig
        ! and write the ground-state-vector to a file 
        iunit = get_free_unit()
        open(iunit, file = 'gs_vec_orig') 
        do i = 1, n_states
            write(iunit, *) gs_vec(i)
        end do
        close(iunit) 

        t_trans_corr_hop = .true. 

        if (associated(brr)) deallocate(brr)
        allocate(brr(2*lat%get_nsites()))
        brr = [(i, i = 1, 2*lat%get_nsites())]

        ! first create the exact similarity transformation 
        allocate(t_mat(size(hamil,1),size(hamil,2)), source = hamil) 
        ! and set the diagonal elements to 0
        do i = 1, size(t_mat,1)
            t_mat(i,i) = 0.0_dp
        end do

        allocate(e_vec_hop(n_states, size(J))) 
        e_vec_hop = 0.0_dp

        do i = 1, size(J) 
            print *, "J = ", J(i)

            write(J_str, *) J(i) 
            filename = 'gs_vec_trans_J_' // trim(adjustl((J_str)))

            hamil_hop = similarity_transform(hamil, J(i) * t_mat)

            ! for the neci hopping hamiltonian: 
            trans_corr_param = J(i)
            ! i need to deallocate umat to recompute for a new value of J! 
            if (allocated(umat_rs_hub_trancorr_hop)) deallocate(umat_rs_hub_trancorr_hop)
            call init_realspace_tests()

            hamil_hop_neci = create_hamiltonian(hilbert_space)

            print *, "diagonalizing the transformed hamiltonian: " 
            call eig(hamil_hop, e_values, e_vec) 

            ! find the ground-state
            ind = minloc(e_values,1) 
            gs_energy = e_values(ind) 
            print *, "transformed ground-state energy: ", gs_energy 

            if (abs(gs_energy - gs_energy_orig) > 1.e-12) then 
                call stop_all("HERE!", "energy incorrect!")
            end if
            ! how do i need to access the vectors to get the energy? 
            gs_vec = abs(e_vec(:,ind))
            call sort(gs_vec)

            gs_vec = gs_vec(n_states:1:-1)

            hf_coeff_hop(i) = gs_vec(1)
            e_vec_hop(:,i) = gs_vec

            neci_eval = calc_eigenvalues(hamil_hop_neci)

            print *, "neci ground-state energy: ", minval(neci_eval) 

            if (abs(gs_energy_orig - minval(neci_eval)) > 1.0e-12) then 
                if (n_states < 20) then 
                    print *, "basis: " 
                    call print_matrix(hilbert_space)
                    print *, "original hamiltonian: "
                    call print_matrix(hamil)
                    print *, "hopping transcorr hamiltonian neci: "
                    call print_matrix(hamil_hop_neci)
                    print *, "hopping transcorr transformed: "

                    print *, "difference: " 
                    allocate(diff(size(hamil_hop,1),size(hamil_hop,2)))
                    diff = hamil_hop - hamil_hop_neci
                    where (abs(diff) < EPS) diff = 0.0_dp

                    call print_matrix(diff)
                end if
                print *, "orig E0:    ", gs_energy_orig
                print *, "hopping E0: ", minval(neci_eval)
!                 print *, "diagonal similarity transformed: "
!                 do l = 1, size(hamil_hop,1)
!                     print *, hamil_hop(l,l)
!                 end do
!                 print *, "diagonal neci hamil: " 
!                 do l = 1, size(hamil_hop_neci,1)
!                     print *, hamil_hop_neci(l,l)
!                 end do

                call stop_all("here", "hopping transcorrelated energy not correct!")
            end if
        end do

        iunit = get_free_unit() 
        open(iunit, file = "hf_coeff_hop")
        do i = 1, size(J)
            write(iunit, *) J(i), hf_coeff_hop(i)
        end do

        ! maybe plot all transformed into one file.. 
        iunit = get_free_unit() 
        open(iunit, file = "gs_vec_hop")
        ! the important quantitiy is J over U i guess or? 
        ! i am not sure.. 
        write(iunit, *) "# J: ", J
        do i = 1, n_states 
            write(iunit, *) e_vec_hop(i,:)
        end do
        close(iunit)

        t_trans_corr_hop = .false.

    end subroutine exact_transcorrelation

    subroutine init_realspace_tests

        get_umat_el => get_umat_el_hub
        call init_tmat(lat) 

        allocate(G1(nbasis)) 
        G1(1:nbasis-1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1
        call init_get_helement_hubbard()
        t_lattice_model = .true.

    end subroutine init_realspace_tests

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

    subroutine gen_excit_rs_hubbard_transcorr_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .true. 
        trans_corr_param = 0.1_dp

        uhub = 10
        bhub = -1
        pSingles = 0.9_dp
        pDoubles = 1.0_dp - pSingles

        lat => lattice('rectangle', 2, 3, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs
        if (associated(brr)) deallocate(brr)
        allocate(brr(nbasis))
        brr = [(i,i=1,nBasis)]

        call init_realspace_tests()

        nel = 6
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,4,5,8,9,12]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr, & 
            "gen_excit_rs_hubbard_transcorr", nI, n_iters, gen_all_excits_r_space_hubbard)

        t_trans_corr_hop = .false.

    end subroutine gen_excit_rs_hubbard_transcorr_test_stoch

    subroutine gen_excit_rs_hubbard_transcorr_uniform_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .true. 
        trans_corr_param = 0.1_dp
        t_uniform_excits = .true.

        uhub = 10
        bhub = -1
        pSingles = 0.9_dp
        pDoubles = 1.0_dp - pSingles

        lat => lattice('rectangle', 2, 3, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs
        if (associated(brr)) deallocate(brr)
        allocate(brr(nbasis))
        brr = [(i,i=1,nBasis)]

        call init_realspace_tests()

        nel = 6
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,4,5,8,9,12]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr_uniform, & 
            "gen_excit_rs_hubbard_transcorr_uniform", nI, n_iters, gen_all_excits_r_space_hubbard)

        t_trans_corr_hop = .false.
        t_uniform_excits = .false.

    end subroutine gen_excit_rs_hubbard_transcorr_uniform_test_stoch

    subroutine gen_excit_rs_hubbard_transcorr_hphf_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .true. 
        trans_corr_param = 0.1_dp
        tHPHF = .true.

        uhub = 10
        bhub = -1
        pSingles = 0.9_dp
        pDoubles = 1.0_dp - pSingles

        lat => lattice('rectangle', 2, 3, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs
        if (associated(brr)) deallocate(brr)
        allocate(brr(nbasis))
        brr = [(i,i=1,nBasis)]

        call init_realspace_tests()

        nel = 6
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,4,5,8,9,12]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        tHPHF = .false.
        t_trans_corr_hop = .false.

    end subroutine gen_excit_rs_hubbard_transcorr_hphf_test_stoch

    subroutine gen_excit_rs_hubbard_transcorr_uniform_hphf_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .true. 
        trans_corr_param = 0.1_dp
        t_uniform_excits = .true.
        tHPHF = .true.

        uhub = 10
        bhub = -1
        pSingles = 0.9_dp
        pDoubles = 1.0_dp - pSingles

        lat => lattice('rectangle', 2, 3, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs
        if (associated(brr)) deallocate(brr)
        allocate(brr(nbasis))
        brr = [(i,i=1,nBasis)]

        call init_realspace_tests()

        nel = 6
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,4,5,8,9,12]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        tHPHF = .false.
        t_trans_corr_hop = .false.
        t_uniform_excits = .false.

    end subroutine gen_excit_rs_hubbard_transcorr_uniform_hphf_test_stoch

    subroutine gen_excit_rs_hubbard_hphf_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .false. 
        tHPHF = .true. 

        uhub = 10
        bhub = -1
        lat => lattice('rectangle', 2, 3, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs

        call init_realspace_tests()

        nel = 6
        allocate(nI(nel))
        nI = [(i, i = 1, nel)]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        print *, ""
        print *, "first for the original model "

        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "" 
        print *, "now for the original transcorrelation: "
        t_trans_corr = .true. 
        trans_corr_param = 0.1_dp 
        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "" 
        print *, "now for the 'new' transcorr: "
        t_trans_corr_new = .true. 
        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "" 
        print *, "now for the sum of the 'new' with the 2-body: "
        t_trans_corr_2body = .true. 
        trans_corr_param_2body = 0.1_dp
        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "now the original and the 2-body: "
        t_trans_corr_new = .false. 
        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, ""
        print *, "now for the 2-body: "
        t_trans_corr = .false. 
        t_trans_corr_new = .false. 
        call run_excit_gen_tester(gen_hphf_excit, "gen_hphf_excit", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        tHPHF = .false.
        t_trans_corr_2body = .false. 

    end subroutine gen_excit_rs_hubbard_hphf_test_stoch

    subroutine gen_excit_rs_hubbard_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .false. 

        uhub = 10
        bhub = -1
        lat => lattice('rectangle', 2, 3, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs

        call init_realspace_tests()

        nel = 4
        allocate(nI(nel))
        nI = [(i, i = 1, nel)]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel 
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        print *, ""
        print *, "first for the original model "

        call run_excit_gen_tester(gen_excit_rs_hubbard, "gen_excit_rs_hubbard", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "" 
        print *, "now for the original transcorrelation: "
        t_trans_corr = .true. 
        trans_corr_param = 0.1_dp 

        call run_excit_gen_tester(gen_excit_rs_hubbard, "gen_excit_rs_hubbard", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "" 
        print *, "now for the 'new' transcorr: "
        t_trans_corr_new = .true. 
        call run_excit_gen_tester(gen_excit_rs_hubbard, "gen_excit_rs_hubbard", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "" 
        print *, "now for the sum of the 'new' with the 2-body: "
        t_trans_corr_2body = .true. 
        trans_corr_param_2body = 0.1_dp
        call run_excit_gen_tester(gen_excit_rs_hubbard, "gen_excit_rs_hubbard", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, "now the original and the 2-body: "
        t_trans_corr_new = .false. 
        call run_excit_gen_tester(gen_excit_rs_hubbard, "gen_excit_rs_hubbard", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        print *, ""
        print *, "now for the 2-body: "
        t_trans_corr = .false. 
        t_trans_corr_new = .false. 
        call run_excit_gen_tester(gen_excit_rs_hubbard, "gen_excit_rs_hubbard", & 
            nI, n_iters, gen_all_excits_r_space_hubbard)

        t_trans_corr_2body = .false. 


    end subroutine gen_excit_rs_hubbard_test_stoch

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

    subroutine get_offdiag_helement_rs_hub_test
        use SystemData, only: nel, nbasis , bhub, t_trans_corr, trans_corr_param_2body, &
                              t_trans_corr_2body, trans_corr_param
        use bit_rep_data, only: niftot, nifd 
        use lattice_mod, only: lattice

        nel = 2 
        nbasis = 8 
        niftot = 0
        nifd = 0
        bhub = -1.0 

        lat => lattice('chain', 4, 1,1, .true., .true., .true.)
        call init_tmat(lat)

        t_trans_corr = .false. 
        t_trans_corr_2body = .false. 

        trans_corr_param_2body = 1.0_dp

        print *, ""
        print *, "testing: get_offdiag_helement_rs_hub" 

        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_rs_hub([1,2],[2,4],.false.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_rs_hub([1,2],[1,3],.true.))
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_rs_hub([1,2],[1,5],.true.))

        t_trans_corr_2body = .true.

        call assert_equals(h_cast(exp(-1.0_dp)), get_offdiag_helement_rs_hub([1,2],[1,3],.true.))
        call assert_equals(h_cast(exp(-1.0_dp)), get_offdiag_helement_rs_hub([1,2],[2,4],.true.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_rs_hub([1,3],[1,3],.true.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_rs_hub([1,3],[2,4],.true.))

        call assert_equals(h_cast(exp(1.0_dp)), get_offdiag_helement_rs_hub([1,4],[1,3],.true.))

        t_trans_corr = .true. 
        trans_corr_param = 1.0

        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_rs_hub([1,4],[1,3],.true.))

        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_rs_hub([1,2],[1,3],.false.))
        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_rs_hub([1,2],[2,4],.false.))

        t_trans_corr = .false. 
        t_trans_corr_2body = .false. 

        nel = -1
        nbasis = -1
        niftot = -1
        nifd = -1

    end subroutine get_offdiag_helement_rs_hub_test

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
end program test_real_space_hubbard

