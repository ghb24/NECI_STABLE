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
                          t_trans_corr_new, t_uniform_excits, tHPHF, &
                          length_x, length_y

    use lattice_mod, only: lat, init_dispersion_rel_cache

    use sort_mod, only: sort

    use util_mod, only: get_free_unit

    use unit_test_helpers

    use dsfmt_interface, only: dsfmt_init

    use fcimcdata, only: pSingles, pDoubles

    use lattice_models_utils, only: gen_all_excits_r_space_hubbard, &
                                    create_hilbert_space_realspace, &
                                    gen_all_doubles_k_space, &
                                    gen_all_singles_rs_hub_default

    use HPHFRandexcitmod, only: gen_hphf_excit, finddetspinsym

    use k_space_hubbard, only: setup_k_space_hub_sym, setup_symmetry_table

    use double_occ_mod, only: count_double_orbs

    use Detbitops, only: count_open_orbs


    implicit none

    integer :: failed_count

    t_new_real_space_hubbard = .true.
    t_lattice_model = .true.

    call dsfmt_init(1)
    call init_fruit()
    ! run the test-driver
    call exact_test()
    call stop_all("here", "for now")
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
!         call run_test_case(gen_excit_rs_hubbard_transcorr_uniform_hphf_test_stoch, "gen_excit_rs_hubbard_transcorr_uniform_hphf_test_stoch")
!         call run_test_case(gen_excit_rs_hubbard_test_stoch, "gen_excit_rs_hubbard_test_stoch")
        call run_test_case(gen_excit_rs_hubbard_transcorr_test_stoch, "gen_excit_rs_hubbard_transcorr_test_stoch")
        call run_test_case(gen_excit_rs_hubbard_transcorr_uniform_test_stoch, "gen_excit_rs_hubbard_transcorr_uniform_test_stoch")
        call stop_all("here", "for now")
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

        integer :: i, n_eig, n_orbs, n_states, iunit
        integer, allocatable :: ni(:), hilbert_space(:,:), nJ(:), flip(:)
        real(dp), allocatable :: e_values(:), e_vecs(:,:), e_vecs_right(:,:), e_vecs_left(:,:)
        integer(n_int), allocatable :: dummy(:,:)
        real(dp) :: j
        real(dp), allocatable :: j_vec(:), e_orig(:), test_evec(:,:)
        real(dp) :: exact_double_occ, e_pot_orig, e_kin_orig, overlap, &
                    double_occ_t, e_pot_t, e_kin_t, e_kin_sim
        HElement_t(dp) :: H_hop, H_spin
        logical :: t_optimize_corr_param, t_do_diag_elements, t_do_exact_transcorr, &
                   t_do_exact_double_occ, t_j_vec, t_input_U, t_calc_singles
        HElement_t(dp), allocatable :: hamil(:,:), t_mat(:,:), hamil_hop(:,:), &
            gutzwiller(:,:), hamil_onsite(:,:), u_mat(:,:), t_mat_t(:,:)
        integer :: n_excits, k, flip_excits
        integer(n_int), allocatable :: singles(:,:), flip_singles(:,:)
        real(dp) :: sum_singles, sum_singles_t, phase
        real(dp), allocatable :: sign_list(:), flip_sign(:)
        logical :: t_start_neel, t_flip, t_input_nel, t_input_lattice

        t_optimize_corr_param  = .false.
        t_do_diag_elements = .false.
        t_do_exact_transcorr = .false.
        t_do_exact_double_occ = .false.
        t_j_vec = .true.
        t_input_U = .true.
        t_calc_singles = .true.
        t_start_neel = .true.
        t_flip = .false.
        phase = 1.0_dp
        t_input_nel = .true.
        t_input_lattice = .true.
        t_twisted_bc = .true.
        twisted_bc(1) = 0.1_dp

        if (t_input_lattice) then
            print *, "input lattice type: (chain,square,rectangle,tilted)"
            read(*,*) lattice_type
            print *, "input x-dim: "
            read(*,*) length_x
            print *, "input y-dim: "
            read(*,*) length_y
        else
            lattice_type = 'square'
            length_x = 2
            length_y = 2
        end if

        if (t_input_nel) then
            print *, "input number of electrons: "
            read(*,*) nel
            print *, "neel-state will be used!"
            t_start_neel = .true.
        else
            nel = 3
            t_start_neel = .true.
        end if

        t_trans_corr_hop = .false.
        lat => lattice(lattice_type, length_x, length_y, 1,.true.,.true.,.true.)
        t_trans_corr_hop = .false.

        if (t_input_U) then
            print *, "input U:"
            read(*,*) uhub
        else
            uhub = 12
        end if
        bhub = -1

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs

        call init_realspace_tests

!         nel = 3
!         allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
!         nI = [1,3,6,7,9,12,13,16,17,20,21,24,25,28,30,31,34,36]

        if (t_start_neel) then
            nI = create_neel_state()
            print *, "neel-state: ", nI
        else
            nI = [1,4,5,8,9,12]
        end if

        if (t_flip) then
            allocate(flip(nel), source = 0)
            call finddetspinsym(nI,flip,nel)
        end if

!         nI = [1,4]
!         nI = [1,2,3,4,5,6,7]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call setup_arr_brr(lat)

!         call init_hopping_transcorr()
!         call setup_symmetry_table()
!         call setup_k_space_hub_sym(lat)
!         call init_dispersion_rel_cache()
!         call init_umat_rs_hub_transcorr()

        if (t_j_vec) then
            j_vec = linspace(-0.5,0.5,100)
        else
            allocate(j_vec(1), source = -0.17_dp)
        end if

        if (t_do_diag_elements) then
            iunit = get_free_unit()
            open(iunit, file = 'diag_elements')
            if (t_calc_singles) then
                call gen_all_singles_rs_hub_default(nI, n_excits, singles, sign_list)
                allocate(nJ(nel), source = 0)
                print *, "number of singles: ", n_excits
                write(iunit, *) "# J, Hd hop, H_ij H_ji*"
                if (t_flip) then
                    call gen_all_singles_rs_hub_default(flip, flip_excits, flip_singles, flip_sign)
                end if
            else
                write(iunit, *) "# J, Hd hop, Hd spin"
            end if
            print *, "H diag hopping: "
            t_recalc_umat = .true.
            t_recalc_tmat = .true.
            do i = 1, size(j_vec)
                t_trans_corr_hop = .true.
                trans_corr_param = j_vec(i)
                ! i need to deallocate umat every time..
                if (allocated(umat_rs_hub_trancorr_hop)) deallocate(umat_rs_hub_trancorr_hop)
                call init_umat_rs_hub_transcorr()
                H_hop = get_diag_helemen_rs_hub_transcorr_hop(nI)

                if (t_calc_singles) then
                    sum_singles = 0.0_dp
                    sum_singles_t = 0.0_dp

                    do k = 1, n_excits
                        call decode_bit_det(nJ, singles(:,k))
                        sum_singles = sum_singles + sign_list(k)*get_helement_lattice(nI,nJ)
                        sum_singles_t = sum_singles_t + sign_list(k)*get_helement_lattice(nJ,nI)
                    end do

                    if (t_flip) then
                        do k = 1, flip_excits
                            call decode_bit_det(nJ, flip_singles(:,k))
                            sum_singles = sum_singles + &
                                phase * flip_sign(k) * get_helement_lattice(flip,nJ)
                            sum_singles_t = sum_singles_t + &
                                phase * flip_sign(k) * get_helement_lattice(nJ,flip)
                        end do
                    end if

                    write(iunit,*) J_vec(i), H_hop, sum_singles, sum_singles_t
                else
                    t_trans_corr_hop = .false.
                    t_spin_dependent_transcorr = .true.
                    call init_tmat_rs_hub_spin_transcorr()
                    H_spin = get_diag_helemen_rs_hub_transcorr_spin(nI)
                    t_spin_dependent_transcorr = .false.
                    write(iunit,*) J_vec(i), H_hop, H_spin
                end if
            end do
            close(iunit)
!             call stop_all("here","now")
        end if

        t_trans_corr_hop = .false.

        call create_hilbert_space_realspace(n_orbs, nOccAlpha, nOccBeta, &
            n_states, hilbert_space, dummy)

        n_eig = size(hilbert_space,2)

        allocate(e_values(n_eig))
        allocate(e_orig(n_eig))
        allocate(e_vecs(n_eig, size(hilbert_space,2)))
        allocate(e_vecs_right(n_eig, size(hilbert_space,2)))
        allocate(e_vecs_left(n_eig, size(hilbert_space,2)))

        print *, "size hilbert: ", size(hilbert_space, 2)
        hamil = create_hamiltonian(hilbert_space)

        print *, "hamil:"
        call print_matrix(hamil)
!         print *, "orig: e-values:"
!         do i = 1, size(hilbert_space,2)
!             print *, e_orig(i)
!         end do

#ifndef __CMPLX
        call eig(hamil, e_orig, e_vecs)

        ! first create the exact similarity transformation
        allocate(t_mat(size(hamil,1),size(hamil,2)), source = hamil)
        ! and set the diagonal elements to 0
        do i = 1, size(t_mat,1)
            t_mat(i,i) = 0.0_dp
        end do
        allocate(gutzwiller(size(hamil,1),size(hamil,2)), source = 0.0_dp)
        allocate(u_mat(size(hamil,1),size(hamil,2)), source = 0.0_dp)
        do i = 1, size(hamil,1)
            gutzwiller(i,i) =  hamil(i,i) / real(uhub,dp)
            u_mat(i,i) = hamil(i,i)
        end do

        hamil_onsite = similarity_transform(hamil, j_vec(1) * gutzwiller)

        hamil_hop = similarity_transform(hamil, j_vec(1) * t_mat)

        t_mat_t = similarity_transform(t_mat, j_vec(1) * gutzwiller)

        call eig(hamil_onsite, e_values, e_vecs_right)

        print *, "onsite e_values correct?: "
        do i = 1, size(hilbert_space,2)
            if (abs(e_orig(i) - e_values(i)) > 1e-7) then
                print *, e_values(i)
            end if
        end do

        call eig(hamil_onsite, e_values, e_vecs_left, .true.)

        nblk = 4
        nkry = 8
        ncycle = 200
        b2l = 1.0e-13_dp

        print *, "nkry: ", nkry
        print *, "nblk: ", nblk
        print *, "b2l: ", b2l
        print *, "ncycle: ", ncycle

        allocate(test_evec(size(hilbert_space,2),1), source = 0.0_dp)
        test_evec(:,1) = matmul(matrix_exponential(-2.0*j_vec(1)*gutzwiller), e_vecs_left(:,1))
        test_evec(:,1) = test_evec(:,1) / norm(test_evec(:,1))

        print *, "test right e_vec: ", dot_product(test_evec(:,1), e_vecs_right(:,1))

        ! try too big systems here:
!         call frsblk_wrapper(hilbert_space, size(hilbert_space, 2), n_eig, e_values, e_vecs)
!         print *, "e_value lanczos:", e_values(1)

        if (t_do_exact_double_occ) then
            exact_double_occ = calc_exact_double_occ(n_states, dummy, e_vecs(:,1))
            print *, "--- orig ---"
            print *, "exact double occupancy: ", exact_double_occ
            e_pot_orig = real(uhub,dp) * exact_double_occ
            print *, "potential energy: ", real(uhub,dp) * exact_double_occ
            print *, "<U>:", dot_product(e_vecs(:,1),matmul(u_mat, e_vecs(:,1)))
            e_kin_orig = dot_product(e_vecs(:,1),matmul(t_mat,e_vecs(:,1)))
            print *, "kinetic energy: " , e_kin_orig
            print *, "E tot: ", e_kin_orig + e_pot_orig
            print *, "E0: ", e_orig(1)
            print *, "transcorr j = ", j_vec(1)
            overlap = dot_product(e_vecs_left(:,1),e_vecs_right(:,1))
            double_occ_t = calc_exact_double_occ(n_states, dummy, e_vecs_right(:,1), &
                e_vecs_left(:,1)) / overlap
            print *, "exact double occupancy transcorr: ", double_occ_t
            e_pot_t = double_occ_t * real(uhub,dp)
            e_kin_t = dot_product(e_vecs_left(:,1), matmul(t_mat, e_vecs_right(:,1))) !/ &
!                 overlap
            e_kin_sim = dot_product(e_vecs_left(:,1), matmul(t_mat_t, e_vecs_right(:,1))) / &
                overlap
            print *, "E pot t: ", e_pot_t
            print *, "E kin t: ", e_kin_t
            print *, "E kin sim: ", e_kin_sim
            print *, "E tot: ", e_pot_t + e_kin_sim
            print *, "E0 t: ", e_values(1)

        end if

        print *, "overlap: ", dot_product(e_vecs_left(:,1),e_vecs_right(:,1))
!         call stop_all("here","now")

        call eig(hamil_hop, e_values, e_vecs_right)

        print *, "hop e_values correct?: "
        do i = 1, size(hilbert_space,2)
            if (abs(e_orig(i) - e_values(i)) > 1e-7) then
                print *, e_values(i)
            end if
        end do

        if (t_do_exact_transcorr) then
            call exact_transcorrelation(lat, nI, j_vec, real(uhub,dp), hilbert_space)
        end if

        if (t_optimize_corr_param) then
            call optimize_correlation_parameters()
        end if

#endif
        call stop_all("here","now")

    end subroutine exact_test

    subroutine optimize_correlation_parameters()
        ! routine to optimize the correlation parameter in real-space
        ! based on hongjuns projection formula
        real(dp), allocatable :: g(:), hf_coeff(:), phi_coeff(:), energy(:)
        integer :: n_j = 100, i, iunit
        integer(n_int), allocatable :: hf_states(:,:), phi_states(:,:)
        real(dp) :: J
        logical :: t_full_ed = .true.

        ! first we need the HF solution in real-space
        ! is still need to decide, if i want to store the basis in ilut or
        ! nI representation..
        call get_real_space_hf(hf_states, hf_coeff)

        ! for the single parameter jastrow ansatz i just need to loop over
        ! J, or do a bisection search to find the optimal J, if it is a
        ! convex function
        g = linspace(0.0,1.0,n_j)
        allocate(energy(n_j), source = 0.0_dp)

        iunit = get_free_unit()
        open(iunit, file = 'corr_param')

        do i = 1, n_j
            ! i need to set the hamiltonian with the correct j
            J = -log(g(i))
            call set_H_transcorr(J)

            ! then i need to apply H to the real-space HF solution
            call apply_H(hf_states, hf_coeff, phi_states, phi_coeff)

            ! then i need to calc the 'expectation' values of the
            ! correlation operator
            ! <hf|t|phi> - <hf|t|hf><hf|phi>
            energy(i) = get_corr_overlap(hf_states,hf_coeff,phi_states,phi_coeff) &
                - get_corr_overlap(hf_states,hf_coeff,hf_states,hf_coeff) &
                * get_overlap(hf_states,hf_coeff,phi_states,phi_coeff)

            ! and then i just output it to post-process with python
            write(iunit,*) J, energy(i)
        end do
        close(iunit)


    end subroutine optimize_correlation_parameters

    real(dp) function get_corr_overlap(phi_L_states, phi_L_coeff, phi_R_states, phi_R_coeff)
        integer(n_int), intent(in) :: phi_L_states(:,:), phi_R_states(:,:)
        real(dp), intent(in) :: phi_L_coeff(:), phi_R_coeff(:)

    end function get_corr_overlap

    real(dp) function get_overlap(phi_L_states, phi_L_coeff, phi_R_states, phi_R_coeff)
        integer(n_int), intent(in) :: phi_L_states(:,:), phi_R_states(:,:)
        real(dp), intent(in) :: phi_L_coeff(:), phi_R_coeff(:)

    end function get_overlap


    subroutine apply_H(hf_states, hf_coeff, phi_states, phi_coeff)
        integer(n_int), intent(in) :: hf_states(:,:)
        real(dp), intent(in) :: hf_coeff(:)
        integer(n_int), intent(out), allocatable :: phi_states(:,:)
        real(dp), intent(out) :: phi_coeff(:)

    end subroutine apply_H

    subroutine set_H_transcorr(j)
        real(dp), intent(in) :: j

    end subroutine set_H_transcorr

    subroutine get_real_space_hf(hf_states, hf_coeff)
        integer(n_int), intent(out), allocatable :: hf_states(:,:)
        real(dp), intent(out), allocatable :: hf_coeff(:)
        ! routine for a given system and filling, which gives me the real-space
        ! HF solution. essentially i only need to do a fourier transformation
        ! of the k-space HF solution, or if i have the eigenvectors of the
        ! t_ij matrix, it is just a unitary transformation of HF

    end subroutine get_real_space_hf

    function calc_exact_double_occ(n_states, ilut_list, e_vec, e_vec_left) result(double_occ)
        integer, intent(in) :: n_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot,n_states)
        real(dp), intent(in) :: e_vec(1,n_states)
        real(dp), intent(in), optional :: e_vec_left(1,n_states)

        real(dp) :: double_occ

        integer :: i

        double_occ = 0.0_dp

        if (present(e_vec_left)) then
            do i = 1, n_states
                double_occ = double_occ + e_vec(1,i) * e_vec_left(1,i) &
                    * real(count_double_orbs(ilut_list(0:nifd,i)),dp)
            end do

        else
            do i = 1, n_states
                double_occ = double_occ + e_vec(1,i)**2 &
                    * real(count_double_orbs(ilut_list(0:nifd,i)),dp)
            end do
        end if

    end function calc_exact_double_occ

#ifndef __CMPLX
    subroutine exact_transcorrelation(lat, nI, J, U, hilbert_space)
        class(lattice), intent(in) :: lat
        integer, intent(in) :: nI(nel)
        real(dp) :: J(:), U
        integer, intent(in) :: hilbert_space(:,:)
        character(*), parameter :: this_routine = "exact_transcorrelation"

        integer :: n_states, iunit, ind, i, k, l, flip(nel)
        real(dp), allocatable :: e_values(:), e_vec(:,:), gs_vec(:)
        real(dp) :: gs_energy_orig, gs_energy, hf_coeff_hop(size(J)), gs_energy_spin, &
                    hf_coeff_spin(size(J))
        HElement_t(dp), allocatable :: hamil(:,:), hamil_hop(:,:), hamil_hop_neci(:,:), &
                                       diff(:,:), hamil_spin(:,:), hamil_spin_neci(:,:)
        HElement_t(dp), allocatable :: t_mat(:,:), t_mat_spin(:,:)
        real(dp), allocatable :: neci_eval(:), e_vec_hop(:,:), e_vec_spin(:,:), neci_spin_eval(:)
        real(dp), allocatable :: e_vec_hop_left(:,:)
        character(30) :: filename, J_str
        logical :: t_calc_singles, t_flip, t_norm_inside, t_norm_inside_sen
        real(dp), allocatable :: neel_states(:), singles(:), j_opt(:), &
            norm_inside(:), norm_inside_left(:), norm_inside_sen(:), &
            norm_inside_sen_left(:)
        integer :: neel_ind, flip_ind, ic_inside, ic, sen, sen_inside
        real(dp) :: phase
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)

        t_calc_singles = .true.
        ! also consider the spin-flipped of the neel state
        t_flip = .false.
        phase = 1.0_dp

        t_norm_inside = .true.
        ic_inside = 2

        t_norm_inside_sen = .true.
        sen_inside = 6

        ! initialize correctly for transcorrelation tests
        ! just do it for the hopping transcorrelation now!
        n_states = size(hilbert_space,2)
        print *, "creating original hamiltonian: "

        t_trans_corr_2body = .false.
        t_trans_corr = .false.
        t_trans_corr_hop = .false.
        t_spin_dependent_transcorr = .false.

        hamil = create_hamiltonian(hilbert_space)
        t_mat_spin = create_spin_dependent_hopping(hilbert_space)

        print *, "diagonalizing original hamiltonian: "
        allocate(e_values(n_states));        e_values = 0.0_dp
        allocate(e_vec(n_states, n_states)); e_vec = 0.0_dp
        allocate(gs_vec(n_states));          gs_vec = 0.0_dp
        do i = 1, size(hamil,1)
            do k = 1, size(hamil,2)
                if (isnan(real(hamil(i,k))) .or. is_inf(real(hamil(i,k)))) print *, i,k,hamil(i,k)
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
        allocate(e_vec_hop_left(n_states,size(J)), source = 0.0_dp)
        allocate(e_vec_spin(n_states, size(J)))
        e_vec_spin = 0.0_dp

        t_recalc_umat = .true.
        t_recalc_tmat = .true.

        if (n_states <= 16) then
            print *, "hilbert-space:"
            call print_matrix(hilbert_space)
            print *, "original hamiltonian: "
            call print_matrix(hamil)

        end if

        if (t_calc_singles) then
            allocate(neel_states(n_states), source = 0.0_dp)
            allocate(singles(n_states), source = 0.0_dp)
            allocate(j_opt(size(j)), source = 0.0_dp)

            ! find neel-state(only one for now..)
            do i = 1, n_states
                if (all(hilbert_space(:,i) == nI)) then
                    neel_ind = i
                    neel_states(i) = 1.0_dp
                end if
            end do

            if (t_flip) then
                call finddetspinsym(nI,flip,nel)
                do i = 1, n_states
                    if (all(hilbert_space(:,i) == flip)) then
                        flip_ind = i
                        neel_states(i) = phase
                    end if
                end do
                ! and normalize
                neel_states = neel_states/norm(neel_states,2)
            end if

            singles = matmul(hamil,neel_states)
            singles(neel_ind) = 0.0_dp

            if (t_flip) then
                singles(flip_ind) = 0.0_dp
            end if

        end if


        do i = 1, size(J)
            print *, "J = ", J(i)

            write(J_str, *) J(i)
            filename = 'gs_vec_trans_J_' // trim(adjustl((J_str)))

            hamil_hop = similarity_transform(hamil, J(i) * t_mat)

            if (t_calc_singles) then
                j_opt(i) = dot_product(singles, matmul(hamil_hop, neel_states))
            end if

            trans_corr_param = J(i)

            t_trans_corr_hop = .true.
            ! i need to deallocate umat to recompute for a new value of J!
            if (allocated(umat_rs_hub_trancorr_hop)) deallocate(umat_rs_hub_trancorr_hop)
            call init_realspace_tests()

            hamil_hop_neci = create_hamiltonian(hilbert_space)
            t_trans_corr_hop = .false.

            ! for the neci hopping hamiltonian:
            if (n_states <= 16) then
                if (i == 1) then
                    print *, "J: ", J(i)
                    print *, "exact hop hamil: "
                    call print_matrix(hamil_hop)
                    print *, "neci hop hamil: "
                    call print_matrix(hamil_hop_neci)
                end if
            end if

            t_spin_dependent_transcorr = .true.
            if (allocated(umat_rs_hub_trancorr_hop)) deallocate(umat_rs_hub_trancorr_hop)
            if (allocated(tmat_rs_hub_spin_transcorr)) deallocate(tmat_rs_hub_spin_transcorr)
            call init_realspace_tests()

            hamil_spin_neci = create_hamiltonian(hilbert_space)
            t_spin_dependent_transcorr = .false.

            print *, "diagonalizing the transformed hamiltonian: "
            call eig(hamil_hop, e_values, e_vec)

            ! find the ground-state
            ind = minloc(e_values,1)
            gs_energy = e_values(ind)
            print *, "transformed ground-state energy: ", gs_energy

            if (abs(gs_energy - gs_energy_orig) > 1.e-10) then
                call stop_all("HERE!", "energy incorrect!")
            end if
            ! how do i need to access the vectors to get the energy?
            e_vec_hop(:,i) = e_vec(:,ind)
            gs_vec = abs(e_vec(:,ind))
            call sort(gs_vec)

            gs_vec = gs_vec(n_states:1:-1)

            hf_coeff_hop(i) = gs_vec(1)

            ! also do the left-ev for the norm calcs
            call eig(hamil_hop, e_values, e_vec, .true.)
            ind = minloc(e_values,1)
            e_vec_hop_left(:,i) = e_vec(:,ind)

            hamil_spin = similarity_transform(hamil, J(i) * t_mat_spin)

            call eig(hamil_spin, e_values, e_vec)
            ind = minloc(e_values,1)
            gs_energy_spin = e_values(ind)
            print *, "spin-transformed ground-state energy: ", gs_energy_spin

            if (abs(gs_energy_spin - gs_energy_orig) > 1.e-10) then
                call stop_all("HERE", "spin-transformed energy incorrect!")
            end if

            gs_vec = abs(e_vec(:,ind))
            call sort(gs_vec)
            gs_vec = gs_vec(n_states:1:-1)

            hf_coeff_spin(i) = gs_vec(1)
            e_vec_spin(:,i) = gs_vec

            neci_eval = calc_eigenvalues(hamil_hop_neci)

            print *, "neci ground-state energy: ", minval(neci_eval)

            if (abs(gs_energy_orig - minval(neci_eval)) > 1.0e-8) then
                if (n_states < 20) then
                    print *, "hopping transformed NECI eigenvalue wrong"
                    print *, "basis: "
                    call print_matrix(transpose(hilbert_space))
                    print *, "original hamiltonian: "
                    call print_matrix(hamil)
                    print *, "exactly transformed hamiltonian: "
                    call print_matrix(hamil_hop)
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

            neci_spin_eval = calc_eigenvalues(hamil_spin_neci)

            print *, "neci spin ground-state energy: ", minval(neci_spin_eval)

            if (abs(gs_energy_orig - minval(neci_spin_eval)) > 1.0e-10) then
                if (n_states < 20) then
                    print *, "spin transformed NECI eigenvalue wrong"
                    print *, "basis: "
                    call print_matrix(transpose(hilbert_space))
                    print *, "spin-tmat: "
                    call print_matrix(t_mat_spin)
                    print *, "original hamiltonian: "
                    call print_matrix(hamil)
                    print *, "exactly transformed hamiltonian: "
                    call print_matrix(hamil_spin)
                    print *, "spin transcorr hamiltonian neci: "
                    call print_matrix(hamil_spin_neci)

                    print *, "difference: "
                    allocate(diff(size(hamil_hop,1),size(hamil_hop,2)))
                    diff = hamil_spin - hamil_spin_neci
                    where (abs(diff) < EPS) diff = 0.0_dp

                    call print_matrix(diff)
                end if
                print *, "orig E0:    ", gs_energy_orig
                print *, "spin E0: ", minval(neci_spin_eval)
!                 print *, "diagonal similarity transformed: "
!                 do l = 1, size(hamil_hop,1)
!                     print *, hamil_hop(l,l)
!                 end do
!                 print *, "diagonal neci hamil: "
!                 do l = 1, size(hamil_hop_neci,1)
!                     print *, hamil_hop_neci(l,l)
!                 end do

                call stop_all("here", "spin transcorrelated energy not correct!")
            end if

        end do

!         print *, "hamil-spin-exact:"
!         call print_matrix(hamil_spin)
!         print *, "hamil-spin-neci:"
!         call print_matrix(hamil_spin_neci)
!         print *, "hamil-hop-exact:"
!         call print_matrix(hamil_hop)
!         print *, "hamil-hop-neci:"
!         call print_matrix(hamil_hop_neci)
!

        if (t_norm_inside) then
            allocate(norm_inside(size(j)), source = 0.0_dp)
            allocate(norm_inside_left(size(J)), source = 0.0_dp)

            call EncodeBitDet(nI, ilutI)
            do i = 1, size(j)
                do k = 1, n_states
                    call EncodeBitDet(hilbert_space(:,k), ilutJ)
                    ic  = findbitexcitlevel(ilutI,ilutJ)
                    if (ic <= ic_inside) then
                        norm_inside(i) = norm_inside(i) + &
                            e_vec_hop(k,i)**2
                        norm_inside_left(i) = norm_inside_left(i) + &
                            e_vec_hop_left(k,i)**2
                    end if
                end do
            end do

            iunit = get_free_unit()
            open(iunit, file = 'norm_inside')
            write(iunit,*) "# J left right"
            do i = 1, size(j)
                write(iunit,*) J(i), norm_inside(i), norm_inside_left(i)
            end do
            close(iunit)
        end if

        if (t_norm_inside_sen) then
            allocate(norm_inside_sen(size(j)), source = 0.0_dp)
            allocate(norm_inside_sen_left(size(J)), source = 0.0_dp)

            call EncodeBitDet(nI, ilutI)
            do i = 1, size(j)
                do k = 1, n_states
                    call EncodeBitDet(hilbert_space(:,k), ilutJ)
                    sen = count_open_orbs(ilutJ)
                    if (sen >= sen_inside) then
                        norm_inside_sen(i) = norm_inside_sen(i) + &
                            e_vec_hop(k,i)**2
                        norm_inside_sen_left(i) = norm_inside_sen_left(i) + &
                            e_vec_hop_left(k,i)**2
                    end if
                end do
            end do

            norm_inside_sen = sqrt(norm_inside_sen)
            norm_inside_sen_left = sqrt(norm_inside_sen_left)

            iunit = get_free_unit()
            open(iunit, file = 'norm_inside_sen')
            write(iunit,*) "# J left right"
            do i = 1, size(j)
                write(iunit,*) J(i), norm_inside_sen(i), norm_inside_sen_left(i)
            end do
            close(iunit)
        end if



        if (t_calc_singles) then
            iunit = get_free_unit()
            open(iunit, file = 'exact_j_opt')
            do i = 1, size(j)
                write(iunit,*) j(i), j_opt(i)
            end do
            close(iunit)
        end if

        iunit = get_free_unit()
        open(iunit, file = "hf_coeff_hop")
        do i = 1, size(J)
            write(iunit, *) J(i), hf_coeff_hop(i)
        end do
        close(iunit)

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

        iunit = get_free_unit()
        open(iunit, file = "hf_coeff_spin")
        do i = 1, size(J)
            write(iunit, *) J(i), hf_coeff_spin(i)
        end do
        close(iunit)

        ! maybe plot all transformed into one file..
        iunit = get_free_unit()
        open(iunit, file = "gs_vec_spin")
        ! the important quantitiy is J over U i guess or?
        ! i am not sure..
        write(iunit, *) "# J: ", J
        do i = 1, n_states
            write(iunit, *) e_vec_spin(i,:)
        end do
        close(iunit)

        t_trans_corr_hop = .false.

    end subroutine exact_transcorrelation
#endif

    subroutine init_realspace_tests

        if (t_trans_corr_hop .or. t_spin_dependent_transcorr) then
            get_umat_el => get_umat_rs_hub_trans
        else
            get_umat_el => get_umat_el_hub
        end if
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
        trans_corr_param = 0.05_dp

        uhub = 16
        bhub = -1
        pSingles = 0.9_dp
        pDoubles = 1.0_dp - pSingles

        lat => lattice('chain', 4, 1, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs
        if (associated(brr)) deallocate(brr)
        allocate(brr(nbasis))
        brr = [(i,i=1,nBasis)]

        call init_realspace_tests()

        nel = 2
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,4]

        nOccAlpha = 0
        nOccBeta = 0

        do i = 1, nel
            if (is_beta(nI(i))) nOccBeta = nOccBeta + 1
            if (is_alpha(nI(i))) nOccAlpha = nOccAlpha + 1
        end do

        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr, &
            "gen_excit_rs_hubbard_transcorr", nI, n_iters, gen_all_excits_r_space_hubbard)

        nI = [1,2]
        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr, &
            "gen_excit_rs_hubbard_transcorr", nI, n_iters, gen_all_excits_r_space_hubbard)

        nI = [1,6]
        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr, &
            "gen_excit_rs_hubbard_transcorr", nI, n_iters, gen_all_excits_r_space_hubbard)

        nI = [1,8]
        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr, &
            "gen_excit_rs_hubbard_transcorr", nI, n_iters, gen_all_excits_r_space_hubbard)
        nI = [2,3]
        call run_excit_gen_tester(gen_excit_rs_hubbard_transcorr, &
            "gen_excit_rs_hubbard_transcorr", nI, n_iters, gen_all_excits_r_space_hubbard)
        t_trans_corr_hop = .false.

    end subroutine gen_excit_rs_hubbard_transcorr_test_stoch

    subroutine gen_excit_rs_hubbard_transcorr_uniform_test_stoch

        integer, allocatable :: nI(:)
        integer :: n_iters, n_orbs, i

        n_iters = 1000000

        t_trans_corr_hop = .true.
        trans_corr_param = 0.05_dp
        t_uniform_excits = .true.

        uhub = 16
        bhub = -1
        pSingles = 0.9_dp
        pDoubles = 1.0_dp - pSingles

        lat => lattice('chain', 4, 1, 1,.true.,.true.,.true.)

        n_orbs = lat%get_nsites()
        nBasis = 2 * n_orbs
        if (associated(brr)) deallocate(brr)
        allocate(brr(nbasis))
        brr = [(i,i=1,nBasis)]

        call init_realspace_tests()

        nel = 2
        allocate(nI(nel))
!         nI = [(i, i = 1, nel)]
        nI = [1,4]

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

        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ilut, ex, 2))
        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ilut, ex, 0))
        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ilut, ex, -1))
        call assert_equals(0.0_dp, calc_pgen_rs_hubbard(ilut, ex, 3))

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ilut,ex,1))

        ex(1,1) = 2
        ex(2,1) = 8

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ilut,ex,1))

        t_trans_corr = .true.

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ilut,ex,1))

        ex(1,1) = 1
        ex(2,1) = 3

        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ilut,ex,1))

        ex(1,1) = 2
        ex(2,1) = 8
        call assert_equals(0.25_dp, calc_pgen_rs_hubbard(ilut,ex,1))

        trans_corr_param = 1.0
        call assert_equals(1.0/(4.0_dp), calc_pgen_rs_hubbard(ilut,ex,1))

        nel = 3
        nI = [1,2,3]
        call encodebitdet(nI, ilut)
        ex(1,1) = 2
        ex(2,1) = 4
        ! think about the
!         call assert_equals(1.0
        ! some rounding errors below, otherwise correct
        call assert_equals(real(1.0/3.0_dp*1.0/(1.0+exp(1.0_dp))), real(calc_pgen_rs_hubbard(ilut, ex,1)),1e-12)

        ex(2,1) = 8
        call assert_equals(exp(1.0_dp)/(3.0*(1.0+exp(1.0_dp))), calc_pgen_rs_hubbard(ilut, ex,1))

        ex(1,1) = 1
        ex(2,1) = 7
        call assert_equals(1.0/3.0_dp, calc_pgen_rs_hubbard(ilut, ex,1))

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

