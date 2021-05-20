#include "macros.h"

program test_tJ_model

    use SystemData, only: t_tJ_model, t_heisenberg_model, exchange_j, t_lattice_model, &
                          nSpatOrbs
    use tJ_model
    use fruit
    use lattice_mod, only: lat
    use constants, only: maxExcit, pi
    use lattice_models_utils, only : csf_purify, create_heisenberg_fock_space, &
                                     create_heisenberg_fock_space_guga
    use fcimcdata, only: ilutref
    use Detbitops, only: encodebitdet
    use matrix_util, only: eig, print_matrix, store_hf_coeff, my_minloc, my_minval, &
                           print_vec
    use guga_bitRepOps, only: write_guga_list, calcstepvector, calcB_vector_ilut
    use util_mod, only: near_zero
    use guga_matrixElements, only: calcDiagExchangeGUGA_nI
    use bit_rep_data, only: GugaBits, IlutBits
    use dsfmt_interface, only: dsfmt_init, genrand_real2_dSFMT
    use guga_excitations, only: csf_to_sds_ilut, csf_vector_to_sds


    implicit none

    integer :: failed_count
    logical :: t_exact_study

    t_exact_study = .false.
    t_tJ_model = .true.
    t_lattice_model = .true.

    call init_fruit()
    if (t_exact_study) call exact_study()
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
        call run_test_case(get_offdiag_helement_tJ_test, "get_offdiag_helement_tJ_test")

    end subroutine tJ_model_test_driver

    subroutine exact_study

        use SystemData, only: lattice_type, length_x, length_y, nel, orbital_order, &
                              tGUGA, t_input_order, t_open_bc_x, t_open_bc_y, &
                              t_open_bc_z
        use util_mod, only: get_free_unit
        use lattice_models_utils, only: create_all_open_shell_dets
        use binomial_lookup, only: factorial
        use unit_test_helpers, only: create_lattice_hamil_ilut
        use guga_init, only: init_guga
        use bit_reps, only: init_bit_rep
        use guga_excitations, only: create_hamiltonian_guga

        real(dp), allocatable :: hamil(:,:), bosonic_hamil(:,:), hamil_3(:,:)
        complex(dp), allocatable :: hamil_2(:,:), e_vecs_cmplx(:,:)
        integer(n_int), allocatable :: hilbert_space(:,:), ilutBipart(:), &
                                 ilutZigZag(:), ilutVolcano(:), fock_space(:,:), &
                                 sds(:,:), ilutCompact(:)
        real(dp), allocatable :: hf_coeffs(:), off_diag_sum(:), abs_off_diag_sum(:), &
                                 n_diff_sign(:), sum_pos(:), sum_neg(:), &
                                 sum_3_cycle(:), sum_n_cycle(:), e_vecs(:,:), &
                                 e_values(:), sum_non_zeros(:), non_zeros_in_gs(:), &
                                 sign_coherent_gs(:), gs_vec(:), gs_vec_bipart(:), &
                                 bosonic_diff(:), bosonic_energy(:), bos_e_vecs(:,:), &
                                 localisation(:), ref_energies(:), all_weights(:,:), &
                                 all_ref_e(:,:), bandwith(:), gap(:), exchange(:), &
                                 translation_matrix(:,:), single_csf_error(:), &
                                 single_csf_spin_corr(:,:), exact_spin_corr(:), &
                                 bipart_spin_corr(:), orig_spin_corr(:), &
                                 final_energy(:), opt_energy(:), opt_spin_corr(:), &
                                 weights(:), one_orb_entanglement(:), &
                                 two_orb_entanglement(:,:), mutual_entanglement(:,:), &
                                 running_entanglement(:), spin_corr_from_gs(:), &
                                 density_density(:), spin_corr_from_gs_gij(:)
        logical :: t_input_perm, t_periodic, t_input_lattice, t_input_bc, t_input_spin
        logical :: t_input_state, t_input_j, t_translation, t_single_csf_rdm
        logical :: t_annealing, t_entanglement, t_input_guga, t_perms
        real(dp) :: hf_coeff, gs_orig, hf_orig, hf_bipart, gs_bipart, bos_diff
        integer :: n_beta, n_alpha, n_perms, hf_ind, i, iunit_write, tot_spin
        integer :: gs_ind, n_periodic, iunit_read, j, target_state, most_important_ind
        integer :: bosonic_error_ind, n_guga
        real(dp) :: high_ref
        integer, allocatable :: nI(:), ref(:), ref_bipart(:), step(:), most_important_order(:)
        integer, allocatable :: bos_err_order(:), opt_order(:), ref_compact(:), nJ(:)
        integer(n_int), allocatable :: refs(:,:)
        character(2) :: num
        complex(dp) :: fac

        t_input_perm = .true.
        t_input_order = .true.
        t_input_bc = .true.
        t_input_spin = .true.
        t_input_j = .false.

        t_input_lattice = .true.
        t_input_state = .true.
        t_translation = .false.
        t_single_csf_rdm = .false.
        t_annealing = .false.
        t_entanglement = .true.

        t_input_guga = .true.

        t_perms = .true.

        if (t_input_lattice) then
            print *, "input lattice type: (chain,square,rectangle,tilted)"
            read(*,*) lattice_type
            print *, "input x-dim: "
            read(*,*) length_x
            if (lattice_type == 'chain') then
                length_y = 1
            else
                print *, "input y-dim: "
                read(*,*) length_y
            end if
        else
            lattice_type = 'chain'
            length_x = 8
            length_y = 1
        end if

        if (t_input_j) then
            print *, "input Heisenberg J (+ -> antiferro, - -> ferro)"
            read(*,*) exchange_j
        else
            exchange_j = 1.0_dp
        end if

        if (t_input_bc) then
            if (lattice_type == 'chain') then
                print *, "periodic (1..yes, 0..no)"
                read(*,*) n_periodic

                if (n_periodic == 1) then
                    t_periodic = .true.
                else
                    t_periodic = .false.
                end if
                t_open_bc_x = .not. t_periodic
                t_open_bc_y = .not. t_periodic
                t_open_bc_z = .not. t_periodic
            else
                print *, "x periodic? (1..yes,0..no)"
                read(*,*) n_periodic
                if (n_periodic == 1) then
                    t_open_bc_x = .false.
                else
                    t_open_bc_x = .true.
                end if
                print *, "y periodic? (1..yes,0..no)"
                read(*,*) n_periodic
                if (n_periodic == 1) then
                    t_open_bc_y = .false.
                else
                    t_open_bc_y = .true.
                end if
            end if
        end if

        if (t_single_csf_rdm) then
            t_bipartite_order = .false.
            t_input_order = .false.
        end if

        lat => lattice(lattice_type, length_x, length_y, 1, &
            t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y,&
            t_periodic_z = .true.)

        nSpatOrbs = lat%get_nsites()

        if (t_input_guga) then
            print *, "use guga (1) or not (0)?"
            read(*,*) n_guga

            if (n_guga == 1) then
                tGUGA = .true.
            else
                tGUGA = .false.
            end if
        end if


        if (t_input_spin) then
            print *, "input desired total spin"
            read(*,*) tot_spin

            if (mod(nSpatOrbs,2) /= mod(tot_spin,2) .or. tot_spin > nSpatOrbs/2) then
                print *, "incorrect total spin for chosen nSpatOrbs"

            end if
        else
            if (mod(nSpatOrbs,2) == 0) then
                tot_spin = 0
            else
                tot_spin = 1
            end if
        end if

        if (t_input_state) then
            print *, "which state to target? (0..ground-state, 1..first excited state, and so on"
            read(*,*) target_state
            target_state = target_state + 1
        else
            target_state = 1
        end if

        nBasis = nSpatOrbs * 2

        nel = nSpatOrbs
        t_heisenberg_model = .true.
        nifd = 0
        allocate(nI(nel), source = [(2 * i - mod(i,2), i = 1, nel)])
        allocate(nJ(nel), source = 0)
        allocate(step(nSpatOrbs), source = 0)
        allocate(ref(nel), source = [(2 * i - mod(i,2), i = 1, nel)])
        allocate(ref_bipart(nel), source = 0)
        ref_bipart(1:ceiling(nSpatorbs/2.)) = [(2 * i - 1, i = 1, ceiling(nSpatorbs/2.))]
        ref_bipart(ceiling(nSpatorbs/2.)+1:) = [(2 * i, i = ceiling(nSpatorbs/2.)+1,nSpatorbs)]
        allocate(ref_compact(nel), source = 0)
        ref_compact(1) = 1
        ref_compact(2:nel-1) = [(2 * i -1 + mod(i,2), i = 2, nel-1)]
        ref_compact(nel) = 2 * nel
        print *, "nI: ", nI
        print *, "ref-bipart:", ref_bipart
        print *, "ref-Compact: ", ref_compact
        allocate(ilutRef(0:0,1), source = 0_n_int)
        allocate(ilutBipart(0:0), source = 0_n_int)
        allocate(ilutCompact(0:0), source  = 0_n_int)
        call EncodeBitDet(ref_bipart, ilutBipart)
        call EncodeBitDet(ref_compact, ilutCompact)


        call EncodeBitDet(nI, ilutRef)

        allocate(orbital_order(nSpatOrbs), source = [(i, i = 1, nSpatOrbs)])
        allocate(most_important_order(nSpatOrbs), source = 0)
        allocate(bos_err_order(nSpatOrbs), source = 0)

        if (lattice_type == 'chain') then
            orbital_order(1:nSpatOrbs-1:2) = [(i, i = 1, nSpatOrbs/2)]
            orbital_order(2:nSpatOrbs:2) = [(i, i = nSpatOrbs/2 + 1, nSpatOrbs)]
        else if (lattice_type == 'rect' .or. lattice_type == 'rectangle' .or. &
            lattice_type == 'square') then
            if (length_x == 4) then
                orbital_order = [1,5,2,6,7,3,8,4]
            else if (length_x == 5) then
                orbital_order = [1,5,2,6,7,3,8,4,9,10]
            else if (length_x == 3) then
                if (length_y == 3) then
                    orbital_order = [1,6,2,8,3,7,4,9,5]
                else
                    orbital_order = [1,5,2,6,3,4]
                end if
            end if
        else if (lattice_type == 'tilted') then
            orbital_order = [1,2,5,3,6,4,7,8]
        end if

        call init_bit_rep()

        if (tGUGA) then
            call init_guga()
            call init_guga_heisenberg_model()
            call init_get_helement_heisenberg_guga()
        else
            call init_heisenberg_model()
            call init_get_helement_heisenberg()
        end if

        print *, "orbital order: ", orbital_order
        call print_vec(orbital_order, "bipartite-order")

        n_beta = ceiling(nSpatOrbs / 2.) + tot_spin / 2
        n_alpha = floor(nSpatOrbs / 2.) - tot_spin / 2
        print *, "tot spin: ", tot_spin
        print *, "n-beta: ", n_beta
        print *, "n-alpha: ", n_alpha

        hilbert_space = create_all_open_shell_dets(nSpatOrbs, n_alpha, n_beta)
        ! and then I need to loop over the permutations
        if (tGUGA) then
            hilbert_space = csf_purify(hilbert_space, tot_spin, nel)
            hamil = create_hamiltonian_guga(hilbert_space)
        else
            hamil = create_lattice_hamil_ilut(hilbert_space)
        end if

        iunit_write = get_free_unit()
        open(iunit_write, file = 'hilbert-space', status = 'replace', action = 'write')
        if (tGUGA) then
            do i = 1, size(hilbert_space,2)
                step = calcStepvector(hilbert_space(:,i))
                do j = 1, nSpatOrbs - 1
                    write(iunit_write, '(i2)', advance = 'no') step(j)
                end do
                write(iunit_write, '(i2)', advance = 'yes') step(nSpatOrbs)
            end do
        else
            do i = 1, size(hilbert_space,2)
                call decode_bit_det(nJ, hilbert_space(:,i))
                do j = 1, nel - 1
                    write(iunit_write, '(i2)', advance = 'no') nJ(j)
                end do
                write(iunit_write, '(i2)', advance = 'yes') nJ(nel)
            end do
        end if
        close(iunit_write)

        iunit_read = get_free_unit()

        if (t_single_csf_rdm) then
            nI = [1,3,6,7,10,11,14,15,18,19,22,23,26,27, 30, 32]
            allocate(ilutZigZag(0:0), source = 0_n_int)
            call EncodeBitDet(nI, ilutZigZag)

            nI = [1,3,5,7,9,11,13,16,17, 20, 22, 24, 26, 28, 30, 32]
            allocate(ilutVolcano(0:0), source = 0_n_int)
            call EncodeBitDet(nI, ilutVolcano)
            ! for every csf in the hilbert space I want to calculate the
            ! spin-spin correlation function now.. (via the 2-rdm this
            ! would be possible.. figure it out..
            ! optimally i should do everything in here or?
            ! and the exact spin-spin is given!
            print *, "reading in 'exact-spin-corr'"
            open(iunit_read, file = 'exact-spin-corr', status = 'old', action = 'read')
            allocate(exact_spin_corr(nSpatOrbs), source = 0.0_dp)
            allocate(bipart_spin_corr(nSpatorbs), source = 0.0_dp)
            allocate(orig_spin_corr(nSpatorbs), source = 0.0_dp)

            do i = 1, nSpatOrbs
                read(iunit_read, *) exact_spin_corr(i)
            end do
            close(iunit_read)
            print *, "done"

            print *, "looping over all CSFs: ", size(hilbert_space,2)
            print *, "to produce single CSF spin-correlation functions"
            allocate(single_csf_spin_corr(size(hilbert_space,2), nSpatOrbs), source = 0.0_dp)
            allocate(single_csf_error(size(hilbert_space,2)), source = 0.0_dp)
            do i = 1, size(hilbert_space, 2)
                single_csf_spin_corr(i,:) = &
                    calc_single_csf_spin_corr(hilbert_space(:,i), start = 1)
                single_csf_error(i) = &
                    sqrt(sum((single_csf_spin_corr(i,:) - exact_spin_corr)**2))
            end do

            bipart_spin_corr = calc_single_csf_spin_corr(ilutBipart, start = 1)
            orig_spin_corr = calc_single_csf_spin_corr(ilutRef(:,1), start = 1)
            call print_vec(single_csf_error, 'single-csf-spin-corr-err')
            call print_vec(bipart_spin_corr, 'bipart-spin-corr')
            call print_vec(orig_spin_corr, 'orig_spin_corr')

            i = minloc(single_csf_error,1)


            print *, "csf with smallest spin-corr error:"
            call write_det_guga(6, hilbert_space(:,i))
            call print_vec(single_csf_spin_corr(i,:), 'single-csf-spin-corr')

            ! print *, "20 'best' CSFs:"
            ! do i = 2, 20
            !     j = my_minloc(single_csf_error, i)
            !     call write_det_guga(6, hilbert_space(:,j))
            ! end do

            print *, "exchange matrix corresponding to this CSF:"

            iunit_write = get_free_unit()
            open(iunit_write, file = 'exchange-matrix-opt', &
                status = 'replace', action = 'write')
            call print_matrix(guga_exchange_matrix(hilbert_space(:,i), nSpatOrbs), iunit_write)
            close(iunit_write)

            if (t_annealing) then
                allocate(final_energy(size(hilbert_space,2)), source = 0.0_dp)
                allocate(opt_energy(size(hilbert_space,2)), source = 0.0_dp)
                allocate(opt_spin_corr(nSpatorbs), source = 0.0_dp)
                ! do i = 1, size(hilbert_space,2)
                    i = 1311
                    print *, "i: ", i
                    call write_det_guga(6, hilbert_space(:,1311))
                    call simulated_annealing(guga_exchange_matrix(hilbert_space(:,i),&
                        nSpatorbs), spin_free_exchange, opt_order, final_energy(i), &
                        opt_energy(i))
                ! end do
                ! call print_vec(final_energy, "final-annealing-energies", t_index = .true.)
                ! call print_vec(opt_energy, "optimal-annealing-energies", t_index = .true.)
                nI = [1,14,2,5,10,7,8,3,9,13,6,11,16,12,15,5]
                opt_spin_corr = calc_single_csf_spin_corr(hilbert_space(:,i), start = 1)
                call print_vec(opt_spin_corr(nI), "opt-spin-corr")
                print *, "spin-corr diff: ", sqrt(sum((exact_spin_corr - &
                    opt_spin_corr(nI))**2))
            end if

            ! call stop_all("here", "for now")
        end if
        if (t_perms) then
            open(iunit_read, file = 'permutations', status = 'old', action = 'read')
            n_perms = factorial(nSpatOrbs - 1)
        else
            n_perms = 1
        end if

        allocate(e_values(size(hilbert_space,2)), source = 0.0_dp)
        allocate(e_vecs(size(hilbert_space,2),size(hilbert_space,2)), source = 0.0_dp)
        allocate(hf_coeffs(n_perms), source = 0.0_dp)
        allocate(off_diag_sum(n_perms), source = 0.0_dp)
        allocate(abs_off_diag_sum(n_perms), source = 0.0_dp)
        allocate(n_diff_sign(n_perms), source = 0.0_dp)
        allocate(sum_pos(n_perms), source = 0.0_dp)
        allocate(sum_neg(n_perms), source = 0.0_dp)
        allocate(sum_3_cycle(n_perms), source = 0.0_dp)
        allocate(sum_n_cycle(n_perms), source = 0.0_dp)
        allocate(sum_non_zeros(n_perms), source = 0.0_dp)
        allocate(non_zeros_in_gs(n_perms), source = 0.0_dp)
        allocate(sign_coherent_gs(n_perms), source = 0.0_dp)
        allocate(gs_vec(size(hilbert_space,2)), source = 0.0_dp)
        allocate(gs_vec_bipart(size(hilbert_space,2)), source = 0.0_dp)
        allocate(refs(0:0,n_perms), source = 0_n_int)
        allocate(bosonic_diff(n_perms), source = 0.0_dp)
        allocate(bos_e_vecs(size(hilbert_space,2),size(hilbert_space,2)), source = 0.0_dp)
        allocate(bosonic_energy(size(hilbert_space,2)), source = 0.0_dp)
        allocate(localisation(n_perms), source = 0.0_dp)
        allocate(ref_energies(n_perms), source = 0.0_dp)
        allocate(all_weights(size(hilbert_space,2),n_perms), source = 0.0_dp)
        allocate(all_ref_e(size(hilbert_space,2), n_perms), source = 0.0_dp)
        allocate(bandwith(n_perms), source = 0.0_dp)
        allocate(gap(n_perms), source = 0.0_dp)
        allocate(exchange(size(hilbert_space,2)), source = 0.0_dp)

        if (tGUGA) then
            do i = 1, size(hilbert_space,2)
                ! exchange(i) = full_exchange(hilbert_space(:,i), nSpatOrbs)
                iunit_write = get_free_unit()
                write(num, '(i2)') i
                open(iunit_write, file = 'exchange-matrix.' // trim(adjustl(num)), &
                    status = 'replace', action = 'write')
                call print_matrix(guga_exchange_matrix(hilbert_space(:,i), nSpatOrbs), iunit_write)
                close(iunit_write)
            end do
        end if

        call eig(hamil, e_values, e_vecs)

        print *, "CSF eigenvalues: ", e_values
        call store_hf_coeff(e_values, e_vecs, target_state, hf_orig, hf_ind, gs_ind)
        call decode_bit_det(ref, hilbert_space(:, hf_ind))

        gs_orig = e_values(gs_ind)
        gs_vec = e_vecs(:,gs_ind)

        print *, "original Hamiltonian (CSFs)"
        call print_matrix(hamil)


        allocate(spin_corr_from_gs(nel), source = 0.0_dp)
        allocate(spin_corr_from_gs_gij(nel), source = 0.0_dp)
        if (tGUGA) then
            spin_corr_from_gs = spin_corr_chain_ordered(local_spin(hilbert_space, gs_vec))
        else
            spin_corr_from_gs = spin_corr_sds(gs_vec, hilbert_space)
            spin_corr_from_gs_gij = spin_corr_sds_gij(gs_vec, hilbert_space)
            density_density = density_corr_sds(gs_vec, hilbert_space)
        end if

        call print_vec(spin_corr_from_gs, 'spin-corr-from-gs', &
            t_index = .true.)
        if (.not.tGUGA) then
            call print_vec(density_density, 'density-density', &
                t_index = .true.)
            call print_vec(spin_corr_from_gs_gij, 'normalized-spin-corr', &
                t_index = .true.)
            iunit_write = get_free_unit()
            open(iunit_write, file = 'spin-corr-and-dens', status = 'replace', action = 'write')
            write(iunit_write, *) "#  orb   |   <S_1 S_i>   |   <S_1 S_i> - <S_1><S_i> | ni\s n_j\s"
            do i = 1, size(density_density)
                write(iunit_write, '(i3,3G25.17)') i, spin_corr_from_gs(i), &
                    spin_corr_from_gs_gij(i), density_density(i)
            end do
            close(iunit_write)
        end if

        print *, "original GS vec (CSFs): "
        do i = 1, size(hilbert_space,2)
            call write_det_guga(6, hilbert_space(:,i), .false.)
            print *, gs_vec(i)
        end do

        print *, "original GS vec (SDs):"
        call csf_vector_to_sds(hilbert_space, gs_vec, sds, weights)
        call write_guga_list(6, sds)

        ! call init_get_helement_heisenberg()
        ! hamil = create_lattice_hamil_ilut(sds)
        !
        ! print *, "orig Hamiltonian (SDs)"
        ! call print_matrix(hamil)
        !
        ! if (allocated(e_values)) deallocate(e_values)
        ! if (allocated(e_vecs)) deallocate(e_vecs)
        ! allocate(e_values(size(hamil,1)), source = 0.0_dp)
        ! allocate(e_vecs(size(hamil,1), size(hamil,2)), source = 0.0_dp)
        !
        ! call eig(hamil, e_values, e_vecs)
        ! if (t_open_bc_x) then
        !     print *, "SD eigenvalues: ", e_values - (nSpatorbs-1)/4.
        ! else
        !     print *, "SD eigenvalues: ", e_values - (nSpatorbs)/4.
        ! end if
        !
        ! call print_matrix(e_vecs)
        ! call write_guga_list(6, sds)
        !
        if (t_entanglement) then

            call csf_vector_to_sds(hilbert_space, gs_vec, sds, weights)

            allocate(one_orb_entanglement(nSpatorbs), source = 0.0_dp)
            allocate(two_orb_entanglement(nSpatorbs,nSpatorbs), source = 0.0_dp)
            allocate(mutual_entanglement(nSpatorbs, nSpatorbs), source = 0.0_dp)
            allocate(running_entanglement(nSpatorbs-1), source = 0.0_dp)

            call get_entanglement_measures(sds, weights, &
                one_orb_entanglement, two_orb_entanglement, mutual_entanglement, &
                running_entanglement)

            print *, "============= Ground state: ==============="
            print *, "one_orb_entanglement: ", one_orb_entanglement
            print *, "two_orb_entanglement:"
            call print_matrix(two_orb_entanglement)
            print *, "mutual_entanglement: "
            call print_matrix(mutual_entanglement)
            print *, "running_entanglement: ", running_entanglement

            call csf_to_sds_ilut(IlutRef(:,1), sds, weights)

            call get_entanglement_measures(sds, weights, &
                one_orb_entanglement, two_orb_entanglement, mutual_entanglement, &
                running_entanglement)

            print *, "============= Orig state: ==============="
            call print_entanglement(one_orb_entanglement, two_orb_entanglement, &
                mutual_entanglement, running_entanglement)


            call csf_to_sds_ilut(IlutBipart, sds, weights)

            call get_entanglement_measures(sds, weights, &
                one_orb_entanglement, two_orb_entanglement, mutual_entanglement, &
                running_entanglement)

            print *, "============= Bipart state: ==============="
            call print_entanglement(one_orb_entanglement, two_orb_entanglement, &
                mutual_entanglement, running_entanglement)

            call csf_to_sds_ilut(ilutCompact, sds, weights)

            call get_entanglement_measures(sds, weights, &
                one_orb_entanglement, two_orb_entanglement, mutual_entanglement, &
                running_entanglement)

            print *, "============= Compact state: ==============="
            call print_entanglement(one_orb_entanglement, two_orb_entanglement, &
                mutual_entanglement, running_entanglement)


            call stop_all("for now", "here")
        end if

        call print_vec(local_spin(hilbert_space, gs_vec), "loc-spin-gs-orig", &
            t_index = .true., t_zero = .true.)
        call print_vec(spin_corr_chain_ordered(local_spin(hilbert_space, gs_vec)), &
            "spin-corr-gs-orig", t_index = .true.)

        ! call print_vec(spin_corr_chain_ordered(local_spin(hilbert_space(:,hf_ind:hf_ind+1), &
            ! [1.0_dp,0.0_dp])),'test',t_index = .true.)

        iunit_write = get_free_unit()
        open(iunit_write, file = 'orig-hamil', status = 'replace', action = 'write')
        call print_matrix(hamil, iunit_write)
        close(iunit_write)

        call print_vec(gs_vec, "gs-vec-orig")

        bosonic_hamil = create_bosonic_hamil(hamil)
        call eig(bosonic_hamil, bosonic_energy, bos_e_vecs)

        bos_diff = gs_orig - my_minval(bosonic_energy, target_state)

        iunit_write = get_free_unit()
        open(iunit_write, file = 'orig-and-bipartite-data', status = 'replace', action = 'write')
        write(iunit_write, *) "# hf-coeff     off-diag-sum    abs-off-diag-sum    &
            &n-diff-sign    sum-pos     sum-neg     sum-3-cycle     sum-n-cycle &
            &   sum-non-zeros   non-zeros-in-gs     bosonic-diff    localisation &
            &   ref-energies    gap     bandwith    sign-coherent-gs"
        write(iunit_write, '(G25.17)', advance = 'no') hf_orig
        write(iunit_write, '(G25.17)', advance = 'no') sum_off_diag(hamil)
        write(iunit_write, '(G25.17)', advance = 'no') sum_off_diag(abs(hamil))
        write(iunit_write, '(G25.17)', advance = 'no') num_diff_sign(hamil)
        write(iunit_write, '(G25.17)', advance = 'no') sum_signed_off_diag(hamil, 1.0_dp)
        write(iunit_write, '(G25.17)', advance = 'no') sum_signed_off_diag(hamil, -1.0_dp)
        write(iunit_write, '(G25.17)', advance = 'no') cycle_flow(hamil, hf_ind, 3)
        write(iunit_write, '(G25.17)', advance = 'no') cycle_flow(hamil, hf_ind, size(hamil,1))
        write(iunit_write, '(G25.17)', advance = 'no') sum_non_zero_off_diag(hamil)
        write(iunit_write, '(G25.17)', advance = 'no') get_num_non_zeros(gs_vec)
        write(iunit_write, '(G25.17)', advance = 'no') gs_orig - my_minval(bosonic_energy, target_state)
        write(iunit_write, '(G25.17)', advance = 'no') norm(gs_vec,4)
        write(iunit_write, '(G25.17)', advance = 'no') hamil(hf_ind,hf_ind)
        write(iunit_write, '(G25.17)', advance = 'no') my_minval(diag(hamil),2) - minval(diag(hamil))
        write(iunit_write, '(G25.17)', advance = 'no') maxval(diag(hamil)) - minval(diag(hamil))
        write(iunit_write, '(G25.17)', advance = 'yes') log2real(is_vec_sign_coherent(gs_vec))

        lat => lattice(lattice_type, length_x, length_y, 1, &
            t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y, &
            t_periodic_z = .true., t_bipartite_order = .true.)

        if (tGUGA) then
            ! call init_guga_heisenberg_model()
            call init_get_helement_heisenberg_guga()
            hamil = create_hamiltonian_guga(hilbert_space)
        else
            call init_get_helement_heisenberg()
            hamil = create_lattice_hamil_ilut(hilbert_space)
        end if

        call eig(hamil, e_values, e_vecs)

        print *, "original e-values: "
        print *, e_values

        call store_hf_coeff(e_values, e_vecs, target_state, hf_bipart, hf_ind, gs_ind)

        gs_bipart = e_values(gs_ind)
        gs_vec_bipart = e_vecs(:,gs_ind)
        call decode_bit_det(ref_bipart, hilbert_space(:, hf_ind))

        call print_vec(local_spin(hilbert_space, gs_vec_bipart), "loc-spin-gs-bipart", &
            t_index = .true., t_zero = .true.)
        call print_vec(spin_corr_chain_ordered(local_spin(hilbert_space, gs_vec_bipart)), &
            "spin-corr-gs-bipart", t_index = .true.)

        bosonic_hamil = create_bosonic_hamil(hamil)
        call eig(bosonic_hamil, bosonic_energy, bos_e_vecs)

        write(iunit_write, '(G25.17)', advance = 'no')  hf_bipart
        write(iunit_write, '(G25.17)', advance = 'no')  sum_off_diag(hamil)
        write(iunit_write, '(G25.17)', advance = 'no')  sum_off_diag(abs(hamil))
        write(iunit_write, '(G25.17)', advance = 'no')  num_diff_sign(hamil)
        write(iunit_write, '(G25.17)', advance = 'no')  sum_signed_off_diag(hamil, 1.0_dp)
        write(iunit_write, '(G25.17)', advance = 'no')  sum_signed_off_diag(hamil, -1.0_dp)
        write(iunit_write, '(G25.17)', advance = 'no')  cycle_flow(hamil, hf_ind, 3)
        write(iunit_write, '(G25.17)', advance = 'no')  cycle_flow(hamil, hf_ind, size(hamil,1))
        write(iunit_write, '(G25.17)', advance = 'no')  sum_non_zero_off_diag(hamil)
        write(iunit_write, '(G25.17)', advance = 'no')  get_num_non_zeros(gs_vec_bipart)
        write(iunit_write, '(G25.17)', advance = 'no')  gs_bipart - my_minval(bosonic_energy, target_state)
        write(iunit_write, '(G25.17)', advance = 'no') norm(gs_vec_bipart,4)
        write(iunit_write, '(G25.17)', advance = 'no') hamil(hf_ind,hf_ind)
        write(iunit_write, '(G25.17)', advance = 'no') my_minval(diag(hamil),2) - minval(diag(hamil))
        write(iunit_write, '(G25.17)', advance = 'no') maxval(diag(hamil)) - minval(diag(hamil))
        write(iunit_write, '(G25.17)', advance = 'yes') log2real(is_vec_sign_coherent(gs_vec_bipart))

        close(iunit_write)

        iunit_write = get_free_unit()
        open(iunit_write, file = 'bipartite-hamil', status = 'replace', action = 'write')
        call print_matrix(hamil, iunit_write)
        close(iunit_write)

        call print_vec(gs_vec_bipart, "gs-vec-bipart")

        iunit_write = get_free_unit()
        open(iunit_write, file = 'gs-energy', status = 'replace', action = 'write')
        write(iunit_write, '(G25.17)') gs_orig
        close(iunit_write)


        print *, "size(hilbert_space): ", shape(hilbert_space)
        print *, "n_perms: ", n_perms
        print *, "original gs energy: ", gs_orig
        print *, "original HF coeff: ", hf_orig
        print *, "original gs vec non-zeros: ", get_num_non_zeros(gs_vec)
        print *, "original gs vec sign coherent?:", is_vec_sign_coherent(gs_vec)
        print *, "original reference: ", ref
        print *, "original bosonic diff: ", bos_diff

        print *, "bipartite gs energy:", gs_bipart
        print *, "bipartite hf coeff: ", hf_bipart
        print *, "bipartite gs vec non-zeros: ", get_num_non_zeros(gs_vec_bipart)
        print *, "bipartite gs vec sign coherent?:", is_vec_sign_coherent(gs_vec_bipart)
        print *, "bipartite reference: ", ref_bipart
        print *, "bipartite bosonic diff: ", gs_bipart - my_minval(bosonic_energy, target_state)

        ! call stop_all("here", "now")

        most_important_ind = 0
        bosonic_error_ind = 0

        high_ref = 0.0_dp

        if (t_translation) then

            allocate(hamil_2(size(hamil,1),size(hamil,2)), source = cmplx(0.0_dp,0.0_dp,kind=dp))
            allocate(hamil_3(size(hamil,1),size(hamil,2)), source = 0.0_dp)
            allocate(e_vecs_cmplx(size(hamil,1),size(hamil,2)), source = cmplx(0.0_dp,0.0_dp,kind=dp))

            print *, "orig hamil:"
            call print_matrix(hamil)
            do i = 1, 2
                if (i == 1) orbital_order = [1,3,2,5,4,6]
                if (i == 2) orbital_order = [2,3,1,5,4,6]
                if (i == 3) orbital_order = [4,5,2,3,1,6]
                if (i == 4) orbital_order = [4,5,1,3,2,6]
                ! orbital_order = cshift(orbital_order, 1)

                print *, "orbital_order: ", orbital_order
                lat => lattice(lattice_type, length_x, length_y, 1, &
                    t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y, &
                    t_periodic_z = .true., t_bipartite_order = .true.)

                if (tGUGA) then
                    ! call init_guga_heisenberg_model()
                    call init_get_helement_heisenberg_guga()
                    hamil = create_hamiltonian_guga(hilbert_space)
                else
                    call init_get_helement_heisenberg()
                    hamil = create_lattice_hamil_ilut(hilbert_space)
                end if

                print *, "intermediate:"
                call print_matrix(hamil)
                print *, "eigenvalues:"
                e_values = 0.0_dp
                e_vecs = 0.0_dp
                call eig(hamil, e_values, e_vecs)
                print *, e_values
                print *, "eigenvectors:"
                call print_matrix(e_vecs)

                ! fac = exp(2 * pi * cmplx(0.0,1.0_dp) / nSpatOrbs * i)

                fac = 1.0_dp
                hamil_3 = hamil_3 + fac * hamil

            end do

            print *, "translated hamil:"
            hamil_3 = hamil_3 / 2.
            call print_matrix(hamil_3)

            e_values = 0.0_dp
            e_vecs = 0.0_dp
            ! call eig(hamil_2, e_values, e_vecs_cmplx)
            call eig(hamil_3, e_values, e_vecs)

            print *, "tranlated e-values: "
            print *, e_values
            print *, "tranlated e-vectors: "
            call print_matrix(e_vecs)

        end if

        if (.not. t_perms) then
            call stop_all("here", "now")
        end if
        do i = 1, n_perms

            read(iunit_read,*) orbital_order

            lat => lattice(lattice_type, length_x, length_y, 1, &
                t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y, &
                t_periodic_z = .true., t_bipartite_order = .true.)

            if (tGUGA) then
                ! call init_guga_heisenberg_model()
                call init_get_helement_heisenberg_guga()
                hamil = create_hamiltonian_guga(hilbert_space)
            else
                call init_get_helement_heisenberg()
                hamil = create_lattice_hamil_ilut(hilbert_space)
            end if

            gap(i) = my_minval(diag(hamil), 2) - minval(diag(hamil))
            bandwith(i) = maxval(diag(hamil)) - minval(diag(hamil))
            call eig(hamil, e_values, e_vecs)

            forall (j = 1:size(hilbert_space,2)) all_ref_e(j,i) = hamil(j,j)

            ! n_orbs - 1 because I do not consider cycles

            ! then I need to sort by the eigenvalues to get the GS
            ! and store the largest coeff. from the GS
            call store_hf_coeff(e_values, e_vecs, target_state, hf_coeff, hf_ind, gs_ind)
            gs_vec = e_vecs(:, gs_ind)

            all_weights(:,i) = gs_vec

            if (hf_coeff > high_ref) then
                high_ref = hf_coeff
                most_important_ind = i
                most_important_order = orbital_order
            end if

            ref_energies(i) = hamil(hf_ind,hf_ind)

            bosonic_hamil = create_bosonic_hamil(hamil)
            call eig(bosonic_hamil, bosonic_energy, bos_e_vecs)

            bosonic_diff(i) = e_values(gs_ind) - my_minval(bosonic_energy, target_state)

            if (bosonic_diff(i) < 0.0_dp .and. (.not. near_zero(bosonic_diff(i)))) then
                bosonic_error_ind = i
                bos_err_order = orbital_order
            end if

            refs(:,i) = hilbert_space(:, hf_ind)

            if (abs(gs_orig - e_values(gs_ind)) > 1.e-12_dp) then
                print *, "energies diffier: "
                print *, "orig E: ", gs_orig
                print *, "new E: ", e_values(gs_ind)
            end if

            hf_coeffs(i) = hf_coeff
            ! then I could also look the other quantities.. sum of
            ! off-diagonal Hamiltonian entries?
            off_diag_sum(i) = sum_off_diag(hamil)

            ! and maybe also store the sum of absolute values?
            abs_off_diag_sum(i) = sum_off_diag(abs(hamil))

            ! and maybe also store the difference in number of negative and
            ! positive off-diagonal entries
            n_diff_sign(i) = num_diff_sign(hamil)

            ! and also store 'just' the sum of negative and positive ones
            sum_pos(i) = sum_signed_off_diag(hamil, 1.0_dp)
            sum_neg(i) = sum_signed_off_diag(hamil, -1.0_dp)

            ! and store 'sign flux' to the reference!
            ! https://stackoverflow.com/questions/16436165/detecting-cycles-in-an-adjacency-matrix/16437341#16437341
            ! to get the coefficient of all n cycles which lead again to the
            ! starting state I need the diagonal entry of the reference of
            ! A ^ n
            ! count all 3-cycles as they are the "closest" and most
            ! important and also the sum of all cycles > 2 which go
            ! back to the original
            sum_3_cycle(i) = cycle_flow(hamil, hf_ind, 3)
            sum_n_cycle(i) = cycle_flow(hamil, hf_ind, size(hamil,1))

            sum_non_zeros(i) = sum_non_zero_off_diag(hamil)
            non_zeros_in_gs(i) = get_num_non_zeros(gs_vec)
            sign_coherent_gs(i) = log2real(is_vec_sign_coherent(gs_vec))
            localisation(i) = norm(gs_vec,4)

        end do

        close(iunit_read)

        ! then print out everything, which also closes all files..
        call print_vec(hf_coeffs, "hf-coeffs")
        call print_vec(off_diag_sum, "off-diag-sum")
        call print_vec(abs_off_diag_sum, "abs-off-diag-sum")
        call print_vec(n_diff_sign, "n-diff-sign")
        call print_vec(sum_pos, "sum-pos")
        call print_vec(sum_neg, "sum-neg")
        call print_vec(sum_3_cycle, "sum-3-cycle")
        call print_vec(sum_n_cycle, "sum-n-cycle")
        call print_vec(sum_non_zeros, "sum-non-zeros")
        call print_vec(non_zeros_in_gs, "non-zeros-in-gs")
        call print_vec(sign_coherent_gs, "sign-coherent-gs")
        call print_vec(bosonic_diff, "bosonic-diff")
        call print_vec(localisation, "localisation")
        call print_vec(ref_energies, "ref-energies")
        call print_vec(gap, "gaps")
        call print_vec(bandwith, "bandwith")

        iunit_write = get_free_unit()
        open(iunit_write, file = 'references', status = 'replace', action = 'write')
        do i = 1, n_perms
            step = calcstepvector(refs(:,i))
            do j = 1, nSpatOrbs - 1
                write(iunit_write, '(i2)', advance = 'no') step(j)
            end do
            write(iunit_write, '(i2)', advance = 'yes') step(nSpatOrbs)
        end do
        close(iunit_write)

        iunit_write = get_free_unit()
        open(iunit_write, file = 'all-data', status = 'replace', action = 'write')
        write(iunit_write, *) "# hf-coeff     off-diag-sum    abs-off-diag-sum    &
            &n-diff-sign    sum-pos     sum-neg     sum-3-cycle     sum-n-cycle &
            &   sum-non-zeros   non-zeros-in-gs     bosonic-diff    localisation &
            &   ref-energies    gap     bandwith    sign-coherent-gs"
        do i = 1, n_perms
            write(iunit_write, '(G25.17)',advance = 'no') hf_coeffs(i)
            write(iunit_write, '(G25.17)',advance = 'no') off_diag_sum(i)
            write(iunit_write, '(G25.17)',advance = 'no') abs_off_diag_sum(i)
            write(iunit_write, '(G25.17)',advance = 'no') n_diff_sign(i)
            write(iunit_write, '(G25.17)',advance = 'no') sum_pos(i)
            write(iunit_write, '(G25.17)',advance = 'no') sum_neg(i)
            write(iunit_write, '(G25.17)',advance = 'no') sum_3_cycle(i)
            write(iunit_write, '(G25.17)',advance = 'no') sum_n_cycle(i)
            write(iunit_write, '(G25.17)',advance = 'no') sum_non_zeros(i)
            write(iunit_write, '(G25.17)',advance = 'no') non_zeros_in_gs(i)
            write(iunit_write, '(G25.17)',advance = 'no') bosonic_diff(i)
            write(iunit_write, '(G25.17)',advance = 'no') localisation(i)
            write(iunit_write, '(G25.17)',advance = 'no') ref_energies(i)
            write(iunit_write, '(G25.17)',advance = 'no') gap(i)
            write(iunit_write, '(G25.17)',advance = 'no') bandwith(i)
            write(iunit_write, '(G25.17)',advance = 'yes') sign_coherent_gs(i)
        end do
        close(iunit_write)


        do i = 1, size(hilbert_space,2)
            iunit_write = get_free_unit()
            write(num,'(i2)') i
            open(iunit_write, file = 'weight_and_ref.' // trim(adjustl(num)), status = 'replace', action = 'write')
            do j = 1, n_perms
                write(iunit_write, *) j, all_weights(i,j), all_ref_e(i,j)
            end do
            close(iunit_write)
        end do
        ! also print the vector and hamil of the most important state
        print *, "print vector and hamil for most_important_order: ", most_important_order


        call print_vec(most_important_order, "compact-order")

        orbital_order = most_important_order

        lat => lattice(lattice_type, length_x, length_y, 1, &
            t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y, &
            t_periodic_z = .true., t_bipartite_order = .true.)

        if (tGUGA) then
            ! call init_guga_heisenberg_model()
            call init_get_helement_heisenberg_guga()
            hamil = create_hamiltonian_guga(hilbert_space)
        else
            call init_get_helement_heisenberg()
            hamil = create_lattice_hamil_ilut(hilbert_space)
        end if

        call eig(hamil, e_values, e_vecs)

        ! n_orbs - 1 because I do not consider cycles

        ! then I need to sort by the eigenvalues to get the GS
        ! and store the largest coeff. from the GS
        call store_hf_coeff(e_values, e_vecs, target_state, hf_coeff, hf_ind, gs_ind)
        gs_vec = e_vecs(:, gs_ind)

        print *, "compact Hamiltonian (CSFs)"
        call print_matrix(hamil)

        print *, "compact GS vec (CSFs): "
        do i = 1, size(hilbert_space,2)
            call write_det_guga(6, hilbert_space(:,i), .false.)
            print *, gs_vec(i)
        end do

        print *, "compact GS vec (SDs):"
        call csf_vector_to_sds(hilbert_space, gs_vec, sds, weights)
        ! call csf_to_sds_ilut(hilbert_space(:,2), sds, weights)
        call write_guga_list(6, sds)

        call print_vec(local_spin(hilbert_space, gs_vec), "loc-spin-gs-compact", &
            t_index = .true., t_zero = .true.)
        call print_vec(spin_corr_chain_ordered(local_spin(hilbert_space, gs_vec)), &
            "spin-corr-gs-compact", t_index = .true.)

        iunit_write = get_free_unit()
        open(iunit_write, file = 'most-hamil', status = 'replace', action = 'write')
        call print_matrix(hamil, iunit_write)
        close(iunit_write)

        iunit_write = get_free_unit()
        call print_vec(gs_vec, "gs-vec-most")

        ! and if we have a bosonic energy which is higher, also look at
        ! hamiltonians..
        if (bosonic_error_ind > 0) then
            print *, "we do have a state with higher bosonic energy than fermionic!"

            orbital_order = bos_err_order

            lat => lattice(lattice_type, length_x, length_y, 1, &
                t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y, &
                t_periodic_z = .true., t_bipartite_order = .true.)

            if (tGUGA) then
                ! call init_guga_heisenberg_model()
                call init_get_helement_heisenberg_guga()
                hamil = create_hamiltonian_guga(hilbert_space)
            else
                call init_get_helement_heisenberg()
                hamil = create_lattice_hamil_ilut(hilbert_space)
            end if

            call eig(hamil, e_values, e_vecs)

            ! n_orbs - 1 because I do not consider cycles

            ! then I need to sort by the eigenvalues to get the GS
            ! and store the largest coeff. from the GS
            call store_hf_coeff(e_values, e_vecs, target_state, hf_coeff, hf_ind, gs_ind)
            gs_vec = e_vecs(:, gs_ind)

            iunit_write = get_free_unit()
            open(iunit_write, file = 'bos-error-hamil-ferm', status = 'replace', action = 'write')
            call print_matrix(hamil, iunit_write)
            close(iunit_write)

            bosonic_hamil = create_bosonic_hamil(hamil)
            call eig(bosonic_hamil, bosonic_energy, bos_e_vecs)

            iunit_write = get_free_unit()
            open(iunit_write, file = 'bos-error-hamil-bos', status = 'replace', action = 'write')
            call print_matrix(bosonic_hamil, iunit_write)
            close(iunit_write)

        end if

        call stop_all("here", "for now")

    end subroutine exact_study


    subroutine print_entanglement(one_orb_entanglement, two_orb_entanglement, &
            mutual_entanglement, running_entanglement)
        real(dp), intent(in) :: one_orb_entanglement(:), two_orb_entanglement(:,:), &
                                mutual_entanglement(:,:), running_entanglement(:)

        print *, "one_orb_entanglement: ", one_orb_entanglement
        print *, "two_orb_entanglement:"
        call print_matrix(two_orb_entanglement)
        print *, "mutual_entanglement: "
        call print_matrix(mutual_entanglement)
        print *, "running_entanglement: ", running_entanglement

    end subroutine print_entanglement

    subroutine get_entanglement_measures(hilbert_space, state_vec, &
            one_orb_entanglement, two_orb_entanglement, mutual_entanglement, &
            running_entanglement)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), intent(in) :: state_vec(:)
        real(dp), intent(out) :: one_orb_entanglement(nSpatorbs), &
                                 two_orb_entanglement(nSpatorbs, nSpatorbs), &
                                 mutual_entanglement(nSpatorbs, nSpatorbs), &
                                 running_entanglement(nSpatorbs-1)
        integer :: i, j
        integer, allocatable :: ind(:)
        real(dp) :: test(nSpatorbs-1)


        do i = 1, nSpatorbs
            one_orb_entanglement(i) = get_orbital_entanglement(hilbert_space, &
                state_vec, [i])
        end do

        do i = 1, nSpatorbs
            do j = 1, nSpatorbs
                two_orb_entanglement(i,j) = get_orbital_entanglement(hilbert_space, &
                    state_vec, [i, j])
            end do
        end do

        do i = 1, nSpatorbs
            do j = 1, nSpatorbs
                if (i /= j) then
                    mutual_entanglement(i,j) = one_orb_entanglement(i) + &
                        one_orb_entanglement(j) - two_orb_entanglement(i,j)
                end if
            end do
        end do

        do i = 1, nSpatorbs - 1
            if (allocated(ind)) deallocate(ind)
            allocate(ind(i), source = [(j, j = 1,i)])

            running_entanglement(i) = get_orbital_entanglement(hilbert_space, &
                state_vec, ind)
        end do

        ! do i = 1, nSpatorbs - 1
        !     if (allocated(ind)) deallocate(ind)
        !     allocate(ind(i), source = [(j, j = nSpatorbs - i + 1, nSpatorbs)])
        !     print *, "ind: ", ind
        !
        !     test(i) = get_orbital_entanglement(hilbert_space, &
        !         state_vec, ind)
        !
        ! end do
        !
        ! print *, "test: ", test

    end subroutine get_entanglement_measures

    real(dp) function get_orbital_entanglement(hilbert_space, state_vec, tgt_orbs)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), intent(in) :: state_vec(:)
        integer, intent(in) :: tgt_orbs(:)
        character(*), parameter :: this_routine = "get_orbital_entanglement"

        real(dp), allocatable :: orbital_rdm(:,:), e_values(:), e_vecs(:,:)
        integer :: i

        get_orbital_entanglement = 0.0_dp

        ! for now it only works in the Heisenberg case!
        orbital_rdm = get_orbital_rdm(hilbert_space, state_vec, tgt_orbs)

        ! need to diagonalize
        allocate(e_values(size(orbital_rdm,2)), source = 0.0_dp)
        allocate(e_vecs(size(orbital_rdm,1), size(orbital_rdm,2)), source = 0.0_dp)

        call eig(orbital_rdm, e_values, e_vecs)

        get_orbital_entanglement = -sum(e_values*log(e_values), &
            .not. near_zero(e_values))

    end function get_orbital_entanglement

    function get_orbital_rdm(hilbert_space, state_vec, tgt_orbs) &
            result(orbital_rdm)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), intent(in) :: state_vec(:)
        integer, intent(in) :: tgt_orbs(:)
        real(dp), allocatable :: orbital_rdm(:,:)

        integer(n_int), allocatable :: local_fock_space(:,:), &
                                       bit_mask(:,:)
        integer :: i, j

        local_fock_space = create_heisenberg_fock_space(size(tgt_orbs))
        bit_mask = create_fock_bit_masks(local_fock_space, tgt_orbs)

        allocate(orbital_rdm(size(local_fock_space,2), size(local_fock_space,2)), &
            source = 0.0_dp)

        do i = 1, size(local_fock_space, 2)
            do j = 1, size(local_fock_space, 2)
                orbital_rdm(i,j) = get_orbital_rdm_entry(hilbert_space, &
                    state_vec, bit_mask(:,i), bit_mask(:,j))
            end do
        end do

    end function get_orbital_rdm

    pure real(dp) function get_orbital_rdm_entry(hilbert_space, state_vec, &
                                mask_i, mask_j)
        integer(n_int), intent(in) :: hilbert_space(:,:), mask_i(:), mask_j(:)
        real(dp), intent(in) :: state_vec(:)

        integer :: i, j

        get_orbital_rdm_entry = 0.0_dp

        do i = 1, size(hilbert_space,2)
            if (is_in(hilbert_space(1,i), mask_i(1))) then
                do j = 1, size(hilbert_space, 2)
                    if (is_in(hilbert_space(1,j), mask_j(1))) then
                        if (is_complementary(hilbert_space(1,i), &
                            hilbert_space(1,j), mask_i(1), mask_j(1))) then
                            get_orbital_rdm_entry = &
                                get_orbital_rdm_entry + state_vec(i) * state_vec(j)
                        end if
                    end if
                end do
            end if
        end do

    end function get_orbital_rdm_entry

    pure logical function is_complementary(state_i, state_j, mask_i, mask_j)
        integer(n_int), intent(in) :: state_i, state_j, mask_i, mask_j

        is_complementary = .false.
        if (ieor(state_i, mask_i) == ieor(state_j, mask_j)) is_complementary = .true.


    end function is_complementary

    pure logical function is_in(state, mask)
        integer(n_int), intent(in) :: state, mask

        is_in = .false.
        if (popcnt(iand(state, mask)) == popcnt(mask)) is_in = .true.

    end function is_in

    pure function create_fock_bit_masks(fock_space, orbitals) result(bit_masks)
        integer(n_int), intent(in) :: fock_space(:,:)
        integer, intent(in) :: orbitals(:)
        integer(n_int), allocatable :: bit_masks(:,:)

        integer :: i, j, orb

        allocate(bit_masks(0:0, size(fock_space,2)), source = 0_n_int)

        do i = 1, size(fock_space,2)
            do j = 1, size(orbitals)
                orb = orbitals(j)
                call mvbits(fock_space(:,i), 2 * (j - 1), 2, bit_masks(:,i), 2 * (orb - 1))
            end do
        end do

    end function create_fock_bit_masks

    subroutine simulated_annealing(cost_matrix, connect_matrix, order, final_energy, opt_energy)
        real(dp), intent(in) :: cost_matrix(:,:), connect_matrix(:,:)
        integer, intent(out), allocatable :: order(:)
        real(dp), intent(out) :: final_energy, opt_energy

        integer, allocatable :: old_order(:), new_order(:), all_orders(:,:)
        integer :: i, n_dim, n_tot_cycle, n_update_temp, ind_i, ind_j, j, iunit
        real(dp) :: old_energy, new_energy, diff_energy, temp, temp_change
        real(dp) :: prob, r
        real(dp), allocatable :: all_energies(:)
        logical :: old_mask(size(connect_matrix,1), size(connect_matrix,2)), &
                   new_mask(size(connect_matrix,1), size(connect_matrix,2))

        ASSERT(size(cost_matrix,1) == size(cost_matrix,2))
        ASSERT(size(connect_matrix,1) == size(connect_matrix,2))
        ASSERT(size(cost_matrix,1) == size(connect_matrix,1))

        call dsfmt_init(1)

        n_dim = size(cost_matrix,1)
        allocate(old_order(n_dim), source = [(i, i = 1, n_dim)])
        allocate(new_order(n_dim), source = old_order)

        old_mask = (.not.near_zero(connect_matrix))

        ! print *, "original connectivity matrix:"
        ! do i = 1, size(old_mask,1)
        !     do j = 1, size(old_mask,2) - 1
        !         write(6, '(l3)', advance = 'no') old_mask(i,j)
        !     end do
        !     write(6, '(l3)', advance = 'yes') old_mask(i,ubound(old_mask,2))
        ! end do

        old_energy = sum(cost_matrix, old_mask) / 4.0_dp

        ! print *, "original energy: ", old_energy

        temp = 10.0_dp
        temp_change = 0.5
        n_tot_cycle = 100000
        n_update_temp = 10000
        prob = 0.0_dp

        allocate(all_energies(n_tot_cycle), source = 0.0_dp)
        allocate(all_orders(n_tot_cycle, nSpatorbs), source = 0)

        do i = 1, n_tot_cycle
            ! lower the temperature
            if (mod(i, n_update_temp) == 0) temp = temp * temp_change
            ! if (mod(i, 100) == 0) print *, "i, energy, temp: ", &
            !     i, old_energy, temp

            old_energy = sum(cost_matrix, old_mask) / 4.0_dp
            all_energies(i) = old_energy
            all_orders(i,:) = old_order
            ! draw two random numbers > 1 (i want to keep the order at 1 at the
            ! beginning. this is fine, since this only excludes equivalent
            ! cycles
            call propose_swap(old_mask, old_order, new_mask, new_order)

            new_energy = sum(cost_matrix, new_mask) / 4.0_dp

            diff_energy = old_energy - new_energy
            prob = min(1.0_dp, exp(diff_energy / temp))

            if (genrand_real2_dSFMT() < prob) then
                ! here we accept
                old_mask = new_mask
                old_order = new_order
            end if
        end do

        call print_vec(all_energies, 'annealed-energies', t_index = .true.)
        i = minloc(all_energies,1)

        opt_energy = all_energies(i)
        final_energy = old_energy
        print *, "permutation with lowest energy:", all_orders(i,:)
        call print_vec(all_orders(i,:), 'lowest-annealed-order')

        print *, "final permutation:", all_orders(n_tot_cycle,:)
        call print_vec(all_orders(n_tot_cycle,:), 'final-annealed-order')

        iunit = get_free_unit()
        open(iunit, file = 'final-connection-matrix', status = 'replace', action = 'write')
        print *, "final connectivity matrix:"
        do i = 1, size(old_mask,1)
            write(iunit, '(i3)', advance = 'no') i
            do j = 1, size(old_mask,2) - 1
                write(6, '(l3)', advance = 'no') old_mask(i,j)
                if (old_mask(i,j)) then
                    write(iunit, '(i3)', advance = 'no') j
                end if
            end do
            write(6, '(l3)', advance = 'yes') old_mask(i,ubound(old_mask,2))
            if (old_mask(i,ubound(old_mask,2))) then
                write(iunit, '(i3)', advance = 'yes') size(old_mask,2)
            else
                write(iunit,*)
            end if
        end do
        close(iunit)
        !


    end subroutine simulated_annealing

    subroutine propose_swap(old_mask, old_order, new_mask, new_order)
        logical, intent(in) :: old_mask(:,:)
        integer, intent(in) :: old_order(:)
        logical, intent(out) :: new_mask(:,:)
        integer, intent(out) :: new_order(:)

        logical :: temp_mask(size(old_mask,1),size(old_mask,2))

        integer :: ind_i, ind_j

        ASSERT(size(old_mask,1) == size(old_mask,2))
        ASSERT(size(new_mask,1) == size(new_mask,2))
        ASSERT(size(old_mask,1) == size(new_mask,1))
        ASSERT(size(old_order,1) == size(old_mask,1))
        ASSERT(size(new_order,1) == size(old_mask,1))

        new_order = old_order
        temp_mask = old_mask

        call draw_two_indices(size(old_mask,1), ind_i, ind_j)

        temp_mask(:,ind_i) = old_mask(:,ind_j)
        temp_mask(:,ind_j) = old_mask(:,ind_i)

        new_mask = temp_mask

        new_mask(ind_i,:) = temp_mask(ind_j,:)
        new_mask(ind_j,:) = temp_mask(ind_i,:)

        new_order(ind_i) = old_order(ind_j)
        new_order(ind_j) = old_order(ind_i)

    end subroutine propose_swap

    subroutine draw_two_indices(n_dim, ind_i, ind_j)
        integer, intent(in) :: n_dim
        integer, intent(out) :: ind_i, ind_j

        integer :: i

        ASSERT(n_dim > 2)
        ind_i = 2 + int(genrand_real2_dSFMT() * (n_dim - 1))

        do while (i < 1000)
            ind_j = 2 + int(genrand_real2_dSFMT() * (n_dim - 1))
            if (ind_j /= ind_i) return
        end do
    end subroutine draw_two_indices

    pure function calc_single_csf_spin_corr(ilut, start) result(vec)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: start
        real(dp) :: vec(nSpatorbs)

        integer :: nI(nel), j, start_, end_


        call decode_bit_det(nI, ilut)

        do j = 1, nSpatorbs
            start_ = min(start, j)
            end_ = max(start, j)
            vec(j) = -0.5_dp - calcDiagExchangeGUGA_nI(start_, end_, nI) / 2.0_dp
        end do

    end function calc_single_csf_spin_corr

    pure function translation_matrix(n_orb) result(matrix)
        integer, intent(in) :: n_orb
        real(dp) :: matrix(n_orb, n_orb)

        integer :: i

        matrix = 0.0_dp

        do i = 1, n_orb - 1
            matrix(i, i + 1) = 1.0_dp
        end do

        matrix(n_orb, 1) = 1.0_dp

    end function translation_matrix

    pure function guga_exchange_matrix(state, n_orb) result(exchange)
        integer(n_int), intent(in) :: state(0:)
        real(dp) :: exchange(n_orb, n_orb)
        integer, intent(in) :: n_orb

        integer :: i, j, nI(nel)

        call decode_bit_det(nI, state)

        exchange = 0.0_dp

        do i = 1, n_orb - 1
            do j =  i + 1, n_orb
                exchange(i,j) = calcDiagExchangeGUGA_nI(i, j, nI)
            end do
        end do

    end function guga_exchange_matrix


    pure real(dp) function full_exchange(state, n_orb)
        integer(n_int), intent(in) :: state(0:)
        integer, intent(in) :: n_orb

        integer :: i, j, nI(nel)

        call decode_bit_det(nI, state)

        full_exchange = 0.0_dp


        do i = 1, n_orb - 1
            do j = i + 1, n_orb
                full_exchange = full_exchange + calcDiagExchangeGUGA_nI(i,j,nI)
            end do
        end do

    end function full_exchange

    pure real(dp) function norm(vec, order)
        real(dp), intent(in) :: vec(:)
        integer, intent(in), optional :: order
        integer order_
        def_default(order_, order, 2)

        norm = (sum(abs(vec)**order))**(real(1.0/real(order)))

    end function norm

    pure function create_bosonic_hamil(hamil) result(bosonic)
        real(dp), intent(in) :: hamil(:,:)
        real(dp) :: bosonic(size(hamil,1), size(hamil,2))

        real(dp) :: copy(size(hamil,1), size(hamil,2))
        integer :: i

        copy = hamil

        forall(i = 1:size(copy,1)) copy(i,i) = 0.0_dp

        bosonic = diag_matrix(hamil) - abs(copy)


    end function create_bosonic_hamil

    pure function local_spin(states, weights) result(loc_spin)
        integer(n_int), intent(in) :: states(:,:)
        real(dp), intent(in) :: weights(:)
        real(dp) :: loc_spin(nSpatOrbs)

        integer :: i
        real(dp) :: spin(nSpatOrbs)

        loc_spin = 0.0_dp

        do i = 1, size(weights)

            spin = calcB_vector_ilut(states(:,i)) / 2.0_dp

            loc_spin = loc_spin + weights(i)**2 * spin * (spin + 1.0_dp)

        end do

    end function local_spin

    pure function spin_corr_chain_ordered(loc_spin) result(spin_corr)
        real(dp), intent(in) :: loc_spin(nSpatOrbs)
        real(dp) :: spin_corr(nSpatOrbs)

        integer :: i, j

        spin_corr = 0.0_dp

        ! first site -> same
        spin_corr(1) = loc_spin(1)

        ! check my notes for this derivation
        do i = 2, nSpatOrbs
            spin_corr(i) = (loc_spin(i) - i * loc_spin(1)) / 2.0_dp

            do j = 2, i -1
                spin_corr(i) = spin_corr(i) - (i - j + 1) * spin_corr(j)
            end do
        end do

    end function spin_corr_chain_ordered

    pure function local_sz_vec(vec, hilbert_space) result(loc_sz)
        real(dp), intent(in) :: vec(:)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), allocatable :: loc_sz(:)

        integer :: i

        allocate(loc_sz(nel), source = 0.0_dp)

        do i = 1, nel
            loc_sz(i) = local_sz(vec, hilbert_space, i)
        end do

    end function local_sz_vec

    real(dp) pure function local_sz(vec, hilbert_space, orb)
        real(dp), intent(in) :: vec(:)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        integer, intent(in) :: orb

        integer :: i, nI(nel)

        local_sz = 0.0_dp

        do i = 1, size(vec)
            call decode_bit_det(nI, hilbert_space(:,i))

            if (is_beta(nI(orb))) then
                local_sz = local_sz + vec(i)**2 * 0.5_dp
            else
                local_sz = local_sz - vec(i)**2 * 0.5_dp
            end if
        end do

    end function local_sz

    pure function spin_corr_sds_gij(vec, hilbert_space) result(spin_corr)
        real(dp), intent(in) :: vec(:)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), allocatable :: spin_corr(:)

        integer :: i, j, nJ(nel), ms_0, ms_j
        real(dp) :: loc_spin(nel), loc_0

        allocate(spin_corr(nel), source = 0.0_dp)


        loc_spin = local_sz_vec(vec, hilbert_space)

        spin_corr(1) = 0.25_dp - loc_spin(1)**2

        do i = 2, nel
            do j = 1, size(hilbert_space,2)
                call decode_bit_det(nJ, hilbert_space(:,j))

                ms_0 = get_spin_pn(nJ(1))

                ms_j = get_spin_pn(nJ(i))

                spin_corr(i) = spin_corr(i) + &
                    vec(j)**2 * ms_0 * ms_j / 4.0_dp

            end do

            spin_corr(i) = spin_corr(i) - loc_spin(1) * loc_spin(i)
        end do

    end function spin_corr_sds_gij



    pure function spin_corr_sds(vec, hilbert_space) result(spin_corr)
        real(dp), intent(in) :: vec(:)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), allocatable :: spin_corr(:)

        integer :: i, j, nJ(nel), ms_0, ms_j

        allocate(spin_corr(nel), source = 0.0_dp)

        spin_corr(1) = 0.25_dp

        do i = 2, nel
            do j = 1, size(hilbert_space,2)
                call decode_bit_det(nJ, hilbert_space(:,j))

                ms_0 = get_spin_pn(nJ(1))

                ms_j = get_spin_pn(nJ(i))

                spin_corr(i) = spin_corr(i) + vec(j)**2 * ms_0 * ms_j / 4.0_dp

            end do
        end do

    end function spin_corr_sds

    pure function density_corr_sds(vec, hilbert_space) result(dens)
        integer(n_int), intent(in) :: hilbert_space(:,:)
        real(dp), intent(in) :: vec(:)
        real(dp), allocatable :: dens(:)

        integer :: i, j, nJ(nel)
        real(dp) :: n_up_0, n_do_0, n_up_i, n_do_i
        allocate(dens(nel), source = 0.0_dp)

        dens(1) = 1.0_dp

        do j = 1, size(hilbert_space,2)
            call decode_bit_det(nJ, hilbert_space(:,j))

            if (is_beta(nJ(1))) then
                n_up_0 = 1.0_dp
                n_do_0 = 0.0_dp
            else
                n_up_0 = 0.0_dp
                n_do_0 = 1.0_dp
            end if

            do i = 2, nel

                if (is_beta(nJ(i))) then
                    n_up_i = 1.0_dp
                    n_do_i = 0.0_dp
                else
                    n_up_i = 0.0_dp
                    n_do_i = 1.0_dp
                end if

                dens(i) = dens(i) + vec(j)**2 * (n_up_0 * n_up_i + n_do_0 * n_do_i)

            end do
        end do

    end function density_corr_sds

    pure function diag_matrix(matrix) result(diag)
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: diag(size(matrix,1), size(matrix,2 ))

        integer :: i

        diag = 0.0_dp

        forall(i = 1:size(matrix,1)) diag(i,i) = matrix(i,i)

    end function diag_matrix


    pure function diag(matrix) result(diag_entries)
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: diag_entries(size(matrix,1))

        integer :: i

        diag_entries = 0.0_dp

        forall(i = 1:size(matrix,1)) diag_entries(i) = matrix(i,i)

    end function diag

    elemental pure real(dp) function log2real(a)
        logical, intent(in) :: a

        if (a) then
            log2real = 1.0_dp
        else
            log2real = 0.0_dp
        end if

    end function log2real

    real(dp) function get_num_non_zeros(vector)
        real(dp), intent(in) :: vector(:)

        real(dp) :: copy(size(vector,1))

        copy = vector

        where (abs(copy) <= EPS) copy = 0.0_dp
        where (abs(copy) > EPS) copy = 1.0_dp

        get_num_non_zeros = sum(copy)

    end function get_num_non_zeros

    pure logical function is_vec_sign_coherent(vector)
        real(dp), intent(in) :: vector(:)

        real(dp) :: copy(size(vector,1))

        copy = vector

        where (abs(copy) <= EPS) copy = 0.0_dp
        where (copy < 0.0_dp) copy = -1.0_dp
        where (copy > 0.0_dp) copy = 1.0_dp


        if (abs(sum(copy)) == sum(abs(copy))) then
            is_vec_sign_coherent = .true.
        else
            is_vec_sign_coherent = .false.
        end if

    end function is_vec_sign_coherent


    real(dp) function sum_non_zero_off_diag(matrix)
        real(dp), intent(in) :: matrix(:,:)

        integer :: i
        real(dp) :: copy(size(matrix,1),size(matrix,2))

        copy = matrix

        forall(i = 1:size(copy,1)) copy(i,i) = 0.0_dp

        where (abs(copy) > EPS) copy = 1.0_dp
        where (abs(copy) <= EPS) copy = 0.0_dp

        sum_non_zero_off_diag = sum(copy)


    end function sum_non_zero_off_diag

    subroutine set_diag(matrix, val)
        real(dp), intent(inout) :: matrix(:,:)
        real(dp), intent(in) :: val

        integer :: i

        forall (i = 1:size(matrix,1)) matrix(i,i) = val

    end subroutine set_diag

    real(dp) function sum_off_diag(matrix)
        real(dp), intent(in) :: matrix(:,:)

        integer :: i
        real(dp) :: copy(size(matrix,1),size(matrix,2))

        copy = matrix

        forall (i = 1:size(copy,1)) copy(i,i) = 0.0_dp


        sum_off_diag = sum(copy)

    end function sum_off_diag

    real(dp) function num_diff_sign(matrix)
        real(dp), intent(in) :: matrix(:,:)

        integer :: i
        real(dp) :: copy(size(matrix,1),size(matrix,2))


        copy = matrix

        forall (i = 1:size(copy,1)) copy(i,i) = 0.0_dp

        where (copy < 0.0_dp) copy = -1.0_dp
        where (copy > 0.0_dp) copy =  1.0_dp

        num_diff_sign = sum(copy)

    end function num_diff_sign

    real(dp) function sum_signed_off_diag(matrix, sgn)
        real(dp), intent(in) :: matrix(:,:), sgn

        integer :: i
        real(dp) :: copy(size(matrix,1),size(matrix,2))

        copy = matrix

        forall (i = 1:size(copy,1)) copy(i,i) = 0.0_dp

        where (sgn * matrix < 0.0_dp) copy = 0.0_dp

        sum_signed_off_diag = sum(copy)

    end function sum_signed_off_diag

    real(dp) function cycle_flow(matrix, ind, order)
        real(dp), intent(in) :: matrix(:,:)
        integer, intent(in) :: ind, order

        integer :: i
        real(dp) :: copy(size(matrix,1),size(matrix,2)), &
                    prod(size(matrix,1),size(matrix,2))

        copy = matrix

        forall(i = 1:size(copy,1)) copy(i,i) = 0.0_dp

        prod = matmul(copy, copy)

        do i = 1, order - 2
            prod = matmul(prod, copy)
        end do

        cycle_flow = prod(ind,ind)

    end function cycle_flow


    subroutine init_tJ_model_test

        use SystemData, only: lattice_type, length_x, length_y, nbasis, nel, &
                              nbasis
        use OneEInts, only: tmat2d
        use FciMCData, only: tsearchtau, tsearchtauoption, ilutref
        use CalcData, only: tau
        use procedure_pointers, only: get_umat_el, generate_excitation
        use real_space_hubbard, only: lat_tau_factor
        use bit_rep_data, only: nifd, NIfTot
        use bit_reps, only: init_bit_rep

        print *, ""
        print *, "testing: init_tJ_model: "
        lattice_type = 'chain'
        length_x = 100
        length_y = 1
        nel = 2
        exchange_j = 1
        nbasis = 200
        bhub = -1.0
        tau = 0.0_dp
        call init_bit_rep()
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
        call assert_equals(lat%get_neighbors(1), [3,2], size(lat%get_neighbors(1)))
        call assert_equals(lat%get_neighbors(2), [1,3], size(lat%get_neighbors(2)))
        call assert_equals(lat%get_neighbors(3), [2,1], size(lat%get_neighbors(2)))

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
        use dsfmt_interface, only: dsfmt_init
        use lattice_mod, only: lattice, lattice_deconstructor
        use SystemData, only: nel, bhub, exchange_j, nbasis
        use bit_rep_data, only: niftot, nifd
        use constants, only: n_int, dp
        use fcimcdata, only: excit_gen_store_type

        integer, allocatable :: nI(:), nJ(:)
        integer(n_int), allocatable :: ilutI(:), ilutJ(:)
        integer :: ex(2,maxExcit), ic
        logical :: tpar
        real(dp) :: pgen
        HElement_t(dp) :: hel
        type(excit_gen_store_type) :: store
        logical :: found_all, t_found(6)

        print *, ""
        print *, "testing: gen_excit_tJ_model"

        nel = 2
        lat => lattice('chain', 4, 1, 1, .false., .false., .false.)
        call dsfmt_init(1)
        allocate(nI(nel))
        allocate(nJ(nel))

        nI = [1,4]
        NIfTot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet(nI, ilutI)

        found_all = .false.
        t_found = .false.

        bhub = -1.0
        exchange_j = -4.0
        t_tJ_model = .true.
        t_heisenberg_model = .false.

        nbasis = 8

        call init_tmat(lat)
        call setup_exchange_matrix(lat)

        do while (.not. found_all)
            call gen_excit_tJ_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_found(1)) then
                t_found(1) = .true.
                call assert_equals(2, ic)
                call assert_equals([1,4], ex(1,:),2)
                call assert_equals([2,3], ex(2,:),2)
                call assert_equals(6, int(ilutJ(0)))
                call assert_true(.not. tpar)
                ! no it is not just this probability!!
                ! because it could have also happened that we would have
                ! chosen electron (4) first and then did a flip!
                ! here it gets tricky!
                call assert_equals(real(0.5*(1+2.0/3.0_dp),dp), pgen)
            end if

            if (all(nJ == [1,6]) .and. .not. t_found(2)) then
                t_found(2) = .true.
                call assert_equals(1, ic)
                call assert_equals(4, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(33, int(ilutJ(0)))
                call assert_true(.not. tpar)
                call assert_equals(real(1.0/6.0_dp,dp), pgen)
            end if

            found_all = all(t_found(1:2))
        end do

        nI = [2,3]
        call EncodeBitDet(nI, ilutI)

        found_all = .false.
        t_found = .false.

        do while (.not. found_all)
            call gen_excit_tJ_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
                hel, store)

            if (all(nJ == [1,4]) .and. .not. t_found(1)) then
                t_found(1) = .true.
                call assert_equals(2, ic)
                call assert_equals([2,3],ex(1,:),2)
                call assert_equals([1,4],ex(2,:),2)
                call assert_equals(9, int(ilutJ(0)))
                call assert_true(.not. tpar)
                call assert_equals(real(0.5*(1+2.0/3.0_dp),dp), pgen)
            end if

            if (all(nJ == [2,5]) .and. .not. t_found(2)) then
                t_found(2) = .true.
                call assert_equals(1, ic)
                call assert_equals(3, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(18, int(ilutJ(0)))
                call assert_true(.not. tpar)
                call assert_equals(real(1.0/6.0_dp,dp), pgen)
            end if
            found_all = all(t_found(1:2))

        end do

        nI = [3,6]
        call EncodeBitDet(nI, ilutI)

        found_all = .false.
        t_found = .false.

        do while (.not. found_all)
            call gen_excit_tJ_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
                hel, store)

            if (all(nJ == [1,6]) .and. .not. t_found(1)) then
                t_found(1) = .true.
                call assert_equals(1, ic)
                call assert_true(.not. tpar)
                call assert_equals(real(1.0/6.0_dp,dp), pgen)
            end if

            if (all(nJ == [4,5]) .and. .not. t_found(2)) then
                t_found(2) = .true.
                call assert_equals(2, ic)
                call assert_true(.not. tpar)
                call assert_equals(4.0/6.0_dp, pgen)
            end if

            if (all(nJ == [3,8]) .and. .not. t_found(3)) then
                t_found(3) = .true.
                call assert_equals(1,ic)
                call assert_true(.not. tpar)
                call assert_equals(1.0/6.0_dp, pgen)
            end if
            found_all = all(t_found(1:3))

        end do

        nel = -1
        NIfTot = -1
        nifd = -1
        nbasis = -1

    end subroutine gen_excit_tJ_model_test

    subroutine calc_pgen_tJ_model_test
        use bit_rep_data, only: niftot, nifd
        use SystemData, only: nel, nbasis
        use lattice_mod, only: lattice
        use constants, only: n_int, dp
        use DetBitOps, only: EncodeBitDet

        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2)

        print *, ""
        print *, "testing: calc_pgen_tJ_model"

        NIfTot = 0
        nifd = 0

        nel = 2
        nbasis = 8

        lat => lattice('chain', 4, 1, 1, .false., .false., .false.)
        exchange_j = -4.0
        bhub = 1.0
        call init_tmat(lat)
        call setup_exchange_matrix(lat)

        allocate(ilut(0:niftot))

        call EncodeBitDet([1,4], ilut)

        ex(1,:) = [1,4]
        ex(2,:) = [2,3]

        call assert_equals(0.0_dp, calc_pgen_tJ_model(ilut,ex,0))
        call assert_equals(0.0_dp, calc_pgen_tJ_model(ilut,ex,3))

        call assert_equals(0.5_dp*(1+2.0/3.0_dp), calc_pgen_tJ_model(ilut, ex, 2))

        ex(1,1) = 4
        ex(2,1) = 6

        call assert_equals(1.0/6.0_dp, calc_pgen_tJ_model(ilut,ex,1))

        call EncodeBitDet([2,3], ilut)
        ex(1,:) = [2,3]
        ex(2,:) = [1,4]

        call assert_equals(0.5_dp*(1+2.0/3.0_dp), calc_pgen_tJ_model(ilut, ex, 2))

        ex(1,1) = 3
        ex(2,1) = 5

        call assert_equals(1.0/6.0_dp, calc_pgen_tJ_model(ilut,ex,1))

        call EncodeBitDet([3,6], ilut)

        ex(1,1) = 3
        ex(2,1) = 1

        call assert_equals(1.0/6.0_dp, calc_pgen_tJ_model(ilut,ex,1))

        ex(1,1) = 6
        ex(2,1) = 8

        call assert_equals(1.0/6.0_dp, calc_pgen_tJ_model(ilut,ex,1))

        ex(1,:) = [3,6]
        ex(2,:) = [4,5]
        call assert_equals(4.0/6.0_dp, calc_pgen_tJ_model(ilut, ex, 2))

        nel = -1
        nbasis = -1
        NIfTot = -1
        nifd = -1

    end subroutine calc_pgen_tJ_model_test

    subroutine create_cum_list_tJ_model_test

        use SystemData, only: nel, bhub
        use lattice_mod, only: lattice
        use bit_rep_data, only: NIfTot, nifd
        use constants, only: dp, n_int
        use DetBitOps, only: EncodeBitDet

        integer(n_int), allocatable :: ilut(:)
        integer, allocatable :: ic_list(:)
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum, cpt

        t_tJ_model = .true.
        t_heisenberg_model = .false.

        nel = 2
        NIfTot = 0
        nifd = 0
        nbasis = 8

        bhub = -1.0
        exchange_j = -2.0

        allocate(ilut(0:NIfTot))
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)
        call init_tmat(lat)
        call setup_exchange_matrix(lat)

        call EncodeBitDet([1,2],ilut)
        ! i also have to setup the matrix element calculation.. duh..

        print *, ""
        print *, "testing: create_cum_list_tJ_model"
        call create_cum_list_tJ_model(ilut, 1, [2,4], &
            cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,2.0], real(cum_arr), 2)
        call assert_equals(2.0_dp, cum_sum)
        call assert_equals([1,1], ic_list,2)

        call EncodeBitDet([1,4], ilut)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,2.0], real(cum_arr), 2)
        call assert_equals(2.0_dp, cum_sum)
        call assert_equals([2,1],ic_list,2)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,2.0], real(cum_arr), 2)
        call assert_equals(2.0_dp, cum_sum)
        call assert_equals([2,1],ic_list,2)

        print *, ""
        print *, "and also for a provided target"
        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list,4,cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list,7,cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list,5,cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list,3,cpt)
        call assert_equals(0.0_dp, cpt)

        print *, ""
        print *, "and more electrons: "
        nel = 3
        call EncodeBitDet([1,4,7],ilut)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,1.0], real(cum_arr), 2)
        call assert_equals(1.0_dp, cum_sum)
        call assert_equals([2,0],ic_list,2)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,2.0], real(cum_arr), 2)
        call assert_equals(2.0_dp, cum_sum)
        call assert_equals([2,1],ic_list,2)

        call create_cum_list_tJ_model(ilut, 7, [1,3], cum_arr, cum_sum, ic_list)
        call assert_equals([0.0,1.0], real(cum_arr), 2)
        call assert_equals(1.0_dp, cum_sum)
        call assert_equals([0,1],ic_list,2)

        call EncodeBitDet([1,4,8],ilut)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,2.0], real(cum_arr), 2)
        call assert_equals(2.0_dp, cum_sum)
        call assert_equals([2,2],ic_list,2)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list)
        call assert_equals([1.0,2.0], real(cum_arr), 2)
        call assert_equals(2.0_dp, cum_sum)
        call assert_equals([2,1],ic_list,2)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list, 2, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list, 1, cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list, 6, cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list, 5, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 4, [1,3], cum_arr, cum_sum, ic_list, 7, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list, 7, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list, 8, cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list, 3, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list, 4, cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_tJ_model(ilut, 1, [2,4], cum_arr, cum_sum, ic_list, 5, cpt)
        call assert_equals(0.0_dp, cpt)

        call EncodeBitDet([1,4,7],ilut)
        call create_cum_list_tJ_model(ilut, 7, [1,3], cum_arr, cum_sum, ic_list, 1, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 7, [1,3], cum_arr, cum_sum, ic_list, 2, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 7, [1,3], cum_arr, cum_sum, ic_list, 5, cpt)
        call assert_equals(1.0_dp, cpt)

        call create_cum_list_tJ_model(ilut, 7, [1,3], cum_arr, cum_sum, ic_list, 6, cpt)
        call assert_equals(0.0_dp, cpt)

        nel = -1
        NIfTot = -1
        nbasis = -1

    end subroutine create_cum_list_tJ_model_test

    subroutine gen_excit_heisenberg_model_test
        use dsfmt_interface, only: dsfmt_init
        use lattice_mod, only: lattice, lattice_deconstructor
        use SystemData, only: nel, exchange_j, nbasis
        use bit_rep_data, only: niftot, nifd
        use Detbitops, only: encodebitdet
        use constants, only: n_int, dp
        use fcimcdata, only: excit_gen_store_type

        integer, allocatable :: nI(:), nJ(:)
        integer(n_int), allocatable :: ilutI(:), ilutJ(:)
        integer :: ex(2,maxExcit), ic
        logical :: tpar
        real(dp) :: pgen
        HElement_t(dp) :: hel
        type(excit_gen_store_type) :: store
        logical :: found_all, t_found(6)

        print *, ""
        print *, "testing: gen_excit_heisenberg_model"

        nel = 2
        lat => lattice('chain', 2, 1, 1, .false., .false., .false.)
        call dsfmt_init(1)
        allocate(nI(nel))
        allocate(nJ(nel))

        nI = [1,4]
        NIfTot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet(nI, ilutI)

        found_all = .false.
        t_found = .false.

        exchange_j = 2.0
        t_tJ_model = .false.
        t_heisenberg_model = .true.

        nbasis = 4

        call setup_exchange_matrix(lat)

        call gen_excit_heisenberg_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
            hel, store)

        call assert_equals([2,3], nJ, 2)
        call assert_equals(2, ic)
        call assert_equals([1,4],ex(1,:),2)
        call assert_equals([2,3],ex(2,:),2)
        call assert_true(.not. tpar)
        call assert_equals(1.0_dp, pgen)

        ni = [2,3]
        call EncodeBitDet(nI, ilutI)

        call gen_excit_heisenberg_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
            hel, store)

        call assert_equals([1,4], nJ, 2)
        call assert_equals(2, ic)
        call assert_equals([1,4],ex(2,:),2)
        call assert_equals([2,3],ex(1,:),2)
        call assert_true(.not. tpar)
        call assert_equals(1.0_dp, pgen)

        nel = 4
        nbasis = 8
        lat => lattice('chain', 4, 1, 1, .false., .false., .false.)

        call setup_exchange_matrix(lat)
        deallocate(nI)
        deallocate(nJ)
        allocate(nI(nel))
        allocate(nJ(nel))

        nI = [1, 3, 6, 8]
        call EncodeBitDet(nI, ilutI)

        call gen_excit_heisenberg_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
            hel, store)

        call assert_equals([1,4,5,8], nJ, 4)
        call assert_true(.not. tpar)
        call assert_equals(0.5_dp, pgen)

        nI = [1,4,5,8]
        call EncodeBitDet(ni, iluti)

        do while (.not. found_all)
            call gen_excit_heisenberg_model(nI, ilutI, nJ, ilutJ, 2, ic, ex, tpar, pgen,  &
                hel, store)

            if (all(nJ == [2,3,5,8]) .and. .not. t_found(1)) then
                t_found(1) = .true.
                call assert_true(.not. tpar)
                call assert_equals(3.0/8.0_dp, pgen)
            end if

            if (all(nJ == [1,3,6,8]) .and. .not. t_found(2)) then
                t_found(2) = .true.
                call assert_true(.not. tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            if (all(nJ == [1,4,6,7]) .and. .not. t_found(3)) then
                t_found(3) = .true.
                call assert_true(.not. tpar)
                call assert_equals(3.0/8.0_dp, pgen)
            end if

            found_all = all(t_found(1:3))
        end do

        nel = -1
        NIfTot = -1
        nifd = -1
        nbasis = -1


    end subroutine gen_excit_heisenberg_model_test

    subroutine create_cum_list_heisenberg_test
        use SystemData, only: nel
        use lattice_mod, only: lattice
        use bit_rep_data, only: NIfTot
        use constants, only: dp, n_int
        use DetBitOps, only: EncodeBitDet

        integer(n_int), allocatable :: ilut(:)
        integer, allocatable :: ic_list(:)
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum, cpt

        t_tJ_model = .false.
        t_heisenberg_model = .true.

        nel = 4
        NIfTot = 0
        nbasis = 8

        exchange_j = -2.0

        allocate(ilut(0:NIfTot))
        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)
        call setup_exchange_matrix(lat)

        ! the heisenberg assumes that we have only and all orbitals singly
        ! occupied!
        call EncodeBitDet([1,3,5,7],ilut)
        print *, ""
        print *, "testing: create_cum_list_heisenberg"
        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum)
        call assert_equals([0.0,0.0], real(cum_arr), 2)
        call assert_equals(0.0_dp, cum_sum)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 3, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 7, cpt)
        call assert_equals(0.0_dp, cpt)
        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 5, cpt)
        call assert_equals(0.0_dp, cpt)
        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 4, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 3, [1,5], cum_arr, cum_sum)
        call assert_equals([0.0,0.0], real(cum_arr), 2)
        call assert_equals(0.0_dp, cum_sum)

        call create_cum_list_heisenberg(ilut, 3, [1,5], cum_arr, cum_sum, 1, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 3, [1,5], cum_arr, cum_sum, 5, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 3, [1,5], cum_arr, cum_sum, 2, cpt)
        call assert_equals(0.0_dp, cpt)

        call EncodeBitDet([1,4,6,7], ilut)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum)
        call assert_equals([1.0_dp, 1.0_dp], cum_arr, 2)
        call assert_equals(1.0_dp, cum_sum)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 3, cpt)
        call assert_equals(1.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 4, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 7, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 1, [3,7], cum_arr, cum_sum, 8, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 4, [2,6], cum_arr, cum_sum)
        call assert_equals([1.0_dp, 1.0_dp], cum_arr, 2)
        call assert_equals(1.0_dp, cum_sum)

        call create_cum_list_heisenberg(ilut, 4, [2,6], cum_arr, cum_sum, 2, cpt)
        call assert_equals(1.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 4, [2,6], cum_arr, cum_sum, 6, cpt)
        call assert_equals(0.0_dp, cpt)

        call EncodeBitDet([1,4,5,8], ilut)
        call create_cum_list_heisenberg(ilut, 5, [3,7], cum_arr, cum_sum)
        call assert_equals([1.0_dp,2.0_dp], cum_arr, 2)
        call assert_equals(2.0_dp, cum_sum)

        call create_cum_list_heisenberg(ilut, 5, [3,7], cum_arr, cum_sum, 6, cpt)
        call assert_equals(0.0_dp, cpt)

        call create_cum_list_heisenberg(ilut, 5, [3,7], cum_arr, cum_sum, 3, cpt)
        call assert_equals(0.5_dp, cpt)

        call create_cum_list_heisenberg(ilut, 4, [2,6], cum_arr, cum_sum)
        call assert_equals([1.0_dp,2.0_dp], cum_arr, 2)
        call assert_equals(2.0_dp, cum_sum)

        call create_cum_list_heisenberg(ilut, 4, [2,6], cum_arr, cum_sum, 6, cpt)
        call assert_equals(0.5_dp, cpt)

        nel = -1
        niftot = -1
        nbasis = -1

    end subroutine create_cum_list_heisenberg_test

    subroutine calc_pgen_heisenberg_model_test
        use bit_rep_data, only: niftot, nifd
        use SystemData, only: nel, nbasis
        use lattice_mod, only: lattice
        use constants, only: n_int, dp
        use DetBitOps, only: EncodeBitDet

        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2)

        print *, ""
        print *, "testing: calc_pgen_heisenberg_model"

        NIfTot = 0
        nifd = 0

        nel = 2
        nbasis = 4

        lat => lattice('chain', 2, 1, 1, .false., .false., .false.)
        exchange_j = 1.0
        call setup_exchange_matrix(lat)

        allocate(ilut(0:niftot))

        call EncodeBitDet([1,4], ilut)

        ex(1,:) = [1,4]
        ex(2,:) = [2,3]

        call assert_equals(0.0_dp, calc_pgen_heisenberg_model(ilut, ex, 0))
        call assert_equals(0.0_dp, calc_pgen_heisenberg_model(ilut, ex, 1))
        call assert_equals(0.0_dp, calc_pgen_heisenberg_model(ilut, ex, 3))
        call assert_equals(1.0_dp, calc_pgen_heisenberg_model(ilut, ex, 2))

        call encodebitdet([2,3],ilut)
        ex(1,:) = [2,3]
        ex(2,:) = [1,4]

        call assert_equals(1.0_dp, calc_pgen_heisenberg_model(ilut, ex, 2))

        nel = 4
        nbasis = 8
        lat => lattice('chain', 4, 1, 1, .false., .false., .false.)

        call setup_exchange_matrix(lat)

        call EncodeBitDet([1,3,6,8],ilut)

        ex(1,:) = [3,6]
        ex(2,:) = [4,5]

        call assert_equals(0.5_dp, calc_pgen_heisenberg_model(ilut,ex,2))

        call EncodeBitDet([1,4,5,8], ilut)
        ex(1,:) = [1,4]
        ex(2,:) = [2,3]

        call assert_equals(3.0/8.0_dp, calc_pgen_heisenberg_model(ilut,ex,2))

        ex(1,:) = [4,5]
        ex(2,:) = [3,6]
        call assert_equals(0.25_dp, calc_pgen_heisenberg_model(ilut,ex,2))

        ex(1,:) = [5,8]
        ex(2,:) = [6,7]
        call assert_equals(3.0/8.0_dp, calc_pgen_heisenberg_model(ilut,ex,2))

        nel = -1
        nbasis = -1
        niftot = -1
        nifd = -1

    end subroutine calc_pgen_heisenberg_model_test

    subroutine setup_exchange_matrix_test
        use SystemData, only: nbasis
        use real_space_hubbard, only: lat
        use lattice_mod, only: lattice

        print *, ""
        print *, "testing: setup_exchange_matrix"
        exchange_j = 2.0
        lat => lattice('chain', 2, 1, 1, .true.,.true.,.true.)
        nbasis = 4
        call setup_exchange_matrix(lat)
        call assert_equals([0.0,0.0,0.0,1.0],real(exchange_matrix(1,:)),4)
        call assert_equals([0.0,0.0,1.0,0.0],real(exchange_matrix(2,:)),4)
        call assert_equals([0.0,1.0,0.0,0.0],real(exchange_matrix(3,:)),4)
        call assert_equals([1.0,0.0,0.0,0.0],real(exchange_matrix(4,:)),4)

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
        exchange_j = 2.0

        lat => lattice('chain', 2, 1, 1, .false.,.false.,.false.)
        call init_tmat(lat)
        call setup_exchange_matrix(lat)
        call init_get_helement_tj()

        t_lattice_model = .true.

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
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2],0))
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2]))

        ! but i am not yet sure about the double counting..
        call assert_equals(h_cast(-1.0_dp), get_helement([1,4],[1,4],0))
        call assert_equals(h_cast(-1.0_dp), get_helement([2,3],[2,3]))

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,3]))
        call assert_equals(h_cast(0.0_dp), get_helement([2,4],[2,4]))

        print *, ""
        print *, "for exchange contributions: "
        call assert_equals(h_cast(1.0_dp), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0_dp), get_helement([2,3],[1,4]))

        print *, ""
        print *, "and for bigger systems: "
        nbasis = 6
        lat => lattice('chain', 3, 1, 1, .true., .true., .true.)
        call init_tmat(lat)
        call setup_exchange_matrix(lat)

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,3],0))
        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,3]))

        call assert_equals(h_cast(2.0_dp), get_helement([1,3],[1,5],1))
        call assert_equals(h_cast(-2.0_dp), get_helement([1,3],[3,5],1))

        call assert_equals(h_cast(2.0_dp), get_helement([1,3],[1,5]))
        call assert_equals(h_cast(-2.0_dp), get_helement([1,3],[3,5]))

        niftot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet([1,2],ilutI)
        call encodebitdet([1,2],ilutJ)

        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2],ilutI,ilutJ))

        call encodebitdet([2,4],ilutJ)
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[2,4],ilutI,ilutJ))

        call encodebitdet([2,3],ilutJ)
        call assert_equals(h_cast(-2.0_dp), get_helement([1,2],[2,3],ilutI,ilutJ))

        call encodebitdet([1,4],ilutJ)
        call assert_equals(h_cast(2.0_dp), get_helement([1,2],[1,4],ilutJ,ilutJ))

        call assert_equals(h_cast(1.0_dp), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0_dp), get_helement([1,6],[2,5]))

        call assert_equals(h_cast(-2.0_dp), get_helement([1,4],[4,5]))
        call assert_equals(h_cast(2.0_dp), get_helement([1,4],[1,6]))

        call assert_equals(h_cast(0.0_dp), get_helement([1,4],[3,6]))

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[2,4]))
        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[2,3]))

        lat => lattice('square', 2,2,1,.true.,.true.,.true.)
        nel = 4
        nbasis = 8

        call init_tmat(lat)
        call setup_exchange_matrix(lat)

        call assert_equals(h_cast(0.0_dp), get_helement([1,3,5,7],[1,3,5,7]))
        call assert_equals(h_cast(0.0_dp), get_helement([2,4,6,8],[2,4,6,8]))

        call assert_equals(h_cast(-2.0_dp), get_helement([1,3,6,8],[1,3,6,8]))
        call assert_equals(h_cast(-2.0_dp), get_helement([2,4,5,7],[2,4,5,7]))

        call assert_equals(h_cast(-2.0_dp), get_helement([1,4,5,8],[1,4,5,8]))

        call assert_equals(h_cast(-4.0_dp), get_helement([1,4,6,7],[1,4,6,7]))

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
        exchange_j = 2.0

        lat => lattice('chain', 2, 1, 1, .false.,.false.,.false.)
        call setup_exchange_matrix(lat)
        call init_get_helement_heisenberg()

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
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2],0))
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2]))

        ! but i am not yet sure about the double counting..
        call assert_equals(h_cast(-0.5_dp), get_helement([1,4],[1,4],0))
        call assert_equals(h_cast(-0.5_dp), get_helement([2,3],[2,3]))

        call assert_equals(h_cast(0.5_dp), get_helement([1,3],[1,3]))
        call assert_equals(h_cast(0.5_dp), get_helement([2,4],[2,4]))

        print *, ""
        print *, "for exchange contributions: "
        call assert_equals(h_cast(1.0_dp), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0_dp), get_helement([2,3],[1,4]))

        print *, ""
        print *, "and for bigger systems: "
        nbasis = 6
        lat => lattice('chain', 3, 1, 1, .true., .true., .true.)
        call setup_exchange_matrix(lat)

        call assert_equals(h_cast(0.5_dp), get_helement([1,3],[1,3],0))
        call assert_equals(h_cast(0.5_dp), get_helement([1,3],[1,3]))

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,5],1))
        call assert_equals(h_cast(-0.0_dp), get_helement([1,3],[3,5],1))

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[1,5]))
        call assert_equals(h_cast(-0.0_dp), get_helement([1,3],[3,5]))

        niftot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet([1,2],ilutI)
        call encodebitdet([1,2],ilutJ)

        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,2],ilutI,ilutJ))

        call encodebitdet([2,4],ilutJ)
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[2,4],ilutI,ilutJ))

        call encodebitdet([2,3],ilutJ)
        call assert_equals(h_cast(-0.0_dp), get_helement([1,2],[2,3],ilutI,ilutJ))

        call encodebitdet([1,4],ilutJ)
        call assert_equals(h_cast(0.0_dp), get_helement([1,2],[1,4],ilutJ,ilutJ))

        call assert_equals(h_cast(1.0_dp), get_helement([1,4],[2,3]))
        call assert_equals(h_cast(1.0_dp), get_helement([1,6],[2,5]))

        call assert_equals(h_cast(-0.0_dp), get_helement([1,4],[4,5]))
        call assert_equals(h_cast(0.0_dp), get_helement([1,4],[1,6]))

        call assert_equals(h_cast(0.0_dp), get_helement([1,4],[3,6]))

        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[2,4]))
        call assert_equals(h_cast(0.0_dp), get_helement([1,3],[2,3]))

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

        call assert_equals(h_cast(0.0_dp), get_diag_helement_heisenberg([1,3,5,7]))
        call assert_equals(h_cast(0.0_dp), get_diag_helement_heisenberg([2,4,6,8]))

        call assert_equals(h_cast(0.0_dp), get_diag_helement_heisenberg([1,2,3,4]))
        call assert_equals(h_cast(-1.0_dp), get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-2.0_dp), get_diag_helement_heisenberg([1,4,6,7]))

        t_tJ_model = .false.
        t_heisenberg_model = .true.

        call assert_equals(h_cast(1.0_dp), get_diag_helement_heisenberg([1,3,5,7]))
        call assert_equals(h_cast(1.0_dp), get_diag_helement_heisenberg([2,4,6,8]))

        call assert_equals(h_cast(-0.0_dp), get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-1.0_dp), get_diag_helement_heisenberg([1,4,6,7]))



        lat => lattice('triangle',2,2,1,.true.,.true.,.true.)
        call setup_exchange_matrix(lat)

        t_tJ_model = .true.
        t_heisenberg_model = .false.

        call assert_equals(h_cast(-2.0_dp),get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-2.0_dp), get_diag_helement_heisenberg([1,4,6,7]))

        call assert_equals(h_cast(0.0_dp), get_diag_helement_heisenberg([1,3,5,7]))

        t_tJ_model = .false.
        t_heisenberg_model = .true.

        call assert_equals(h_cast(6.0/4.0_dp), get_diag_helement_heisenberg([1,3,5,7]))

        call assert_equals(h_cast(-0.5_dp),get_diag_helement_heisenberg([1,3,6,8]))
        call assert_equals(h_cast(-0.5_dp), get_diag_helement_heisenberg([1,4,6,7]))


        nel = -1
        nbasis = -1
        NIfTot = -1

    end subroutine get_diag_helement_heisenberg_test

    subroutine get_offdiag_helement_heisenberg_test
        use SystemData, only: nel, nbasis
        use real_space_hubbard, only: lat
        use lattice_mod, only: lattice
        use bit_rep_data, only: NIfTot

        integer :: ex(2,2)

        nel = 2
        nbasis = 4
        NIfTot = 0

        exchange_j = 2.0
        lat => lattice('chain', 2, 1, 1, .true., .true., .true.)
        call setup_exchange_matrix(lat)

        print *, ""
        print *, "testing: get_offdiag_helement_heisenberg"
        ex(1,:) = [1,2]
        ex(2,:) = [3,4]


        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([1,2],ex,.true.))
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([1,2],ex,.false.))

        ex(1,:) = [1,4]
        ex(2,:) = [2,3]

        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_heisenberg([1,4], ex,.false.))
        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_heisenberg([1,4], ex,.true.))

        ex(1,:) = [2,3]
        ex(2,:) = [1,4]

        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_heisenberg([2,3], ex,.false.))

        nel = 4
        nbasis = 8
        lat => lattice('triangle', 2, 2, 1, .true., .true., .true.)
        call setup_exchange_matrix(lat)

        ex(1,:) = [3,6]
        ex(2,:) = [4,5]
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_heisenberg([1,2,3,6],ex,.false.))
        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_heisenberg([1,2,3,6],ex,.true.))

        ex(2,:) = [7,8]

        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([1,2,3,6],ex,.false.))

        ex(1,:) = [1,2]

        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([1,2,3,6],ex,.false.))

        ex(2,:) = [4,5]

        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([1,2,3,6],ex,.false.))

        ex(1,:) = [4,5]
        ex(2,:) = [3,6]

        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_heisenberg([1,2,4,5],ex,.false.))

        ex(1,:) = [5,8]
        ex(2,:) = [2,3]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([1,4,5,8],ex,.false.))

        ex(2,:) = [1,4]
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_heisenberg([2,3,5,8],ex,.false.))

        nel = -1
        nbasis = -1
        NIfTot = -1

    end subroutine get_offdiag_helement_heisenberg_test

    subroutine determine_optimal_time_step_tJ_test
        use SystemData, only: nel
        use lattice_mod, only: lat, determine_optimal_time_step

        t_tJ_model = .true.
        t_heisenberg_model = .false.

        bhub = 1.0
        exchange_j = 1.0
        nel = 2

        lat => lattice('chain',2,1,1,.true.,.true.,.true.)

        print *, ""
        print *, "testing: determine_optimal_time_step_tJ"

        call assert_equals(1.0/real(2*2,dp), determine_optimal_time_step())

        exchange_j = 2.0

        call assert_equals(1.0/real(2*4,dp), determine_optimal_time_step())

        bhub = 4
        call assert_equals(1.0/real(2*8,dp), determine_optimal_time_step())

        nel = -1
        bhub = 0
        exchange_j = 0

    end subroutine determine_optimal_time_step_tJ_test

    subroutine determine_optimal_time_step_heisenberg_test
        use SystemData, only: nel
        use lattice_mod, only: lat, determine_optimal_time_step

        t_tJ_model = .false.
        t_heisenberg_model = .true.

        nel = 2
        exchange_j = 1.0
        lat => lattice('triangle', 3,3,1,.true.,.true.,.true.)

        print *, ""
        print *, "testing: determine_optimal_time_step_heisenberg"
        call assert_equals(1.0/real(2*6,dp), determine_optimal_time_step())

        exchange_j = 2.0

        call assert_equals(1.0/real(2*12,dp), determine_optimal_time_step())

        nel = -1
        exchange_j = 0

    end subroutine determine_optimal_time_step_heisenberg_test

    subroutine get_umat_el_heisenberg_test
        use SystemData, only: nbasis
        use lattice_mod, only: lattice
        use real_space_hubbard, only: lat

        print *, ""
        print *, "testing: get_umat_el_heisenberg"

        lat => lattice('chain', 4, 1, 1, .true., .true., .true.)
        exchange_j = 2.0
        nbasis = 8
        call setup_exchange_matrix(lat)

        call assert_equals(h_cast(0.0_dp), get_umat_el_heisenberg(1,1,1,1))
        call assert_equals(h_cast(0.0_dp), get_umat_el_heisenberg(1,2,3,4))
        call assert_equals(h_cast(0.0_dp), get_umat_el_heisenberg(1,1,2,2))
        call assert_equals(h_cast(0.0_dp), get_umat_el_heisenberg(1,3,1,3))
        call assert_equals(h_cast(0.0_dp), get_umat_el_heisenberg(1,3,3,1))
        call assert_equals(h_cast(0.0_dp), get_umat_el_heisenberg(3,1,1,3))

        ! how do i really encode the exchange in the umat in terms of spatial
        ! orbitals?! i am not sure i do it right
        ! since i access the umat only depending on the electrons k,l usually
        ! it should be fine if i always return the exchange j when the
        ! orbitals fit
        call assert_equals(h_cast(1.0_dp), get_umat_el_heisenberg(1,2,1,2))
        call assert_equals(h_cast(1.0_dp), get_umat_el_heisenberg(1,2,2,1))
        call assert_equals(h_cast(1.0_dp), get_umat_el_heisenberg(2,1,1,2))
        call assert_equals(h_cast(1.0_dp), get_umat_el_heisenberg(2,1,2,1))

        call assert_equals(h_cast(1.0_dp), get_umat_el_heisenberg(1,4,1,4))
        call assert_equals(h_cast(1.0_dp), get_umat_el_heisenberg(2,3,3,2))

        nbasis = -1

    end subroutine get_umat_el_heisenberg_test

    subroutine get_offdiag_helement_tJ_test
        use SystemData, only: nel, t_trans_corr, trans_corr_param, nbasis, bhub, &
                             t_trans_corr_2body
        use bit_rep_data, only: NIfTot
        use lattice_mod, only: lat
        use real_space_hubbard, only: init_tmat

        print *, ""
        print *, "testing: get_offdiag_helement_tJ"

        nbasis = 4
        nel = 2
        NIfTot = 0
        lat => lattice('chain', 2, 1, 1, .false.,.false.,.true.)
        bhub = -1.0
        call init_tmat(lat)

        t_trans_corr = .false.
        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_tJ([1,2],[1,3],.false.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_tJ([1,2],[1,3],.true.))

        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_tJ([1,2],[2,4],.false.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_tJ([1,2],[2,4],.true.))

        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_tJ([1,2],[1,4],.false.))
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_tJ([1,2],[2,3],.true.))

        ! the important test is the transcorrelated influence here
        t_trans_corr = .true.
        trans_corr_param = 1.0

        call assert_equals(h_cast(-1.0*exp(-1.0_dp)),get_offdiag_helement_tJ([1,2],[1,3],.false.))
        call assert_equals(h_cast(1.0*exp(3.0_dp)), get_offdiag_helement_tJ([1,2],[4,2],.true.))

        lat => lattice('square', 2, 2, 1, .true.,.true.,.true.)
        nbasis = 8
        call init_tmat(lat)

        nel = 4

        call assert_equals(h_cast(-exp(1.0_dp)), get_offdiag_helement_tJ([1,4,6,7],[1,3],.false.))
        call assert_equals(h_cast(-exp(1.0_dp)), get_offdiag_helement_tJ([1,4,6,7],[2,4],.false.))
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_tJ([1,4,6,7],[2,3],.false.))
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_tJ([1,4,6,7],[1,7],.false.))

        call assert_equals(h_cast(exp(1.0_dp)), get_offdiag_helement_tJ([1,2,3,6],[1,3],.true.))
        call assert_equals(h_cast(exp(1.0_dp)), get_offdiag_helement_tJ([1,2,3,6],[2,4],.true.))

        call assert_equals(h_cast(exp(1.0_dp)), get_offdiag_helement_tJ([1,2,3,6],[7,5],.true.))
        call assert_equals(h_cast(exp(1.0_dp)), get_offdiag_helement_tJ([1,2,3,6],[8,4],.true.))

        call assert_equals(h_cast(1.0_dp*exp(-3.0_dp)), (get_offdiag_helement_tJ([1,2,7,8],[1,3],.true.)),1e-12_dp)
        call assert_equals(h_cast(1.0*exp(5.0_dp)), get_offdiag_helement_tJ([1,2,7,8],[3,7],.true.))

        t_trans_corr = .false.
        t_trans_corr_2body = .true.
        nbasis = 4
        nel = 2
        NIfTot = 0
        lat => lattice('chain', 2, 1, 1, .false.,.false.,.true.)
        bhub = -1.0
        call init_tmat(lat)

        t_trans_corr = .false.
        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_tJ([1,2],[1,3],.false.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_tJ([1,2],[1,3],.true.))

        call assert_equals(h_cast(-1.0_dp), get_offdiag_helement_tJ([1,2],[2,4],.false.))
        call assert_equals(h_cast(1.0_dp), get_offdiag_helement_tJ([1,2],[2,4],.true.))

        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_tJ([1,2],[1,4],.false.))
        call assert_equals(h_cast(0.0_dp), get_offdiag_helement_tJ([1,2],[2,3],.true.))

        t_trans_corr_2body = .false.

        nel = -1
        nbasis = -1
        NIfTot = -1
        t_trans_corr = .false.

    end subroutine get_offdiag_helement_tJ_test
end program test_tJ_model
