#include "macros.h"

program test_tJ_model

    use SystemData, only: t_tJ_model, t_heisenberg_model, exchange_j, t_lattice_model
    use tJ_model
    use fruit
    use lattice_mod, only: lat
    use constants, only: maxExcit
    use lattice_models_utils, only : csf_purify
    use fcimcdata, only: ilutref
    use Detbitops, only: encodebitdet
    use matrix_util, only: eig, print_matrix

    implicit none

    integer :: failed_count
    logical :: t_exact_study

    t_exact_study = .true.
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

        HElement_t(dp), allocatable :: hamil(:,:)
        integer(n_int), allocatable :: hilbert_space(:,:)
        real(dp), allocatable :: hf_coeffs(:), off_diag_sum(:), abs_off_diag_sum(:), &
                                 n_diff_sign(:), sum_pos(:), sum_neg(:), &
                                 sum_3_cycle(:), sum_n_cycle(:), e_vecs(:,:), &
                                 e_values(:), sum_non_zeros(:), non_zeros_in_gs(:), &
                                 sign_coherent_gs(:), gs_vec(:), gs_vec_bipart(:)
        logical :: t_input_perm, t_periodic, t_input_lattice, t_input_bc
        real(dp) :: hf_coeff, gs_orig, hf_orig, hf_bipart, gs_bipart
        integer :: n_beta, n_alpha, n_orbs, n_perms, hf_ind, i, iunit_write, tot_spin
        integer :: gs_ind, n_periodic, iunit_read
        integer, allocatable :: nI(:)

        t_input_perm = .true.
        t_input_order = .true.
        t_input_bc = .true.

        t_input_lattice = .true.

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

        lat => lattice(lattice_type, length_x, length_y, 1, &
            t_periodic_x = .not. t_open_bc_x, t_periodic_y = .not. t_open_bc_y,&
            t_periodic_z = .true.)

        n_orbs = lat%get_nsites()
        nBasis = n_orbs * 2

        nel = n_orbs
        tGUGA = .true.
        t_heisenberg_model = .true.
        nifd = 0
        exchange_j = 1.0_dp
        allocate(nI(nel), source = [(2 * i - mod(i,2), i = 1, nel)])
        print *, "nI: ", nI
        allocate(ilutRef(0:0,1), source = 0_n_int)

        call EncodeBitDet(nI, ilutRef)

        allocate(orbital_order(n_orbs), source = [(i, i = 1, n_orbs)])
        if (lattice_type == 'chain') then
            orbital_order(1:n_orbs/2) = [(2 * i - 1, i = 1, n_orbs/2)]
            orbital_order(n_orbs/2 + 1:n_orbs) = [(2 * i, i = 1, n_orbs/2)]
        else if (lattice_type == 'rect' .or. lattice_type == 'rectangle') then
            orbital_order = [1,5,2,6,7,3,8,4]
        else if (lattice_type == 'tilted') then
            orbital_order = [1,2,4,6,3,5,7,8]
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


        n_alpha = n_orbs / 2
        n_beta = n_orbs / 2
        tot_spin = 0

        hilbert_space = create_all_open_shell_dets(n_orbs, n_alpha, n_beta)
        ! and then I need to loop over the permutations
        if (tGUGA) then
            hilbert_space = csf_purify(hilbert_space, tot_spin, n_orbs, nel)
            hamil = create_hamiltonian_guga(hilbert_space)
        else
            hamil = create_lattice_hamil_ilut(hilbert_space)
        end if

        iunit_read = get_free_unit()
        open(iunit_read, file = 'permutations', status = 'old', action = 'read')

        n_perms = factorial(n_orbs - 1)

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
        allocate(gs_vec(n_perms), source = 0.0_dp)
        allocate(gs_vec_bipart(n_perms), source = 0.0_dp)

        call eig(hamil, e_values, e_vecs)
        call store_gs_hf_coeff(e_values, e_vecs, hf_orig, hf_ind, gs_ind)

        gs_orig = e_values(gs_ind)
        gs_vec = e_vecs(:,gs_ind)

        iunit_write = get_free_unit()
        open(iunit_write, file = 'orig-and-bipartite-data', status = 'replace', action = 'write')
        write(iunit_write, *) "# hf-coeff     off-diag-sum    abs-off-diag-sum    &
            &n-diff-sign    sum-pos     sum-neg     sum-3-cycle     sum-n-cycle &
            &   sum-non-zeros   non-zeros-in-gs     sign-coherent-gs"
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

        call store_gs_hf_coeff(e_values, e_vecs, hf_bipart, hf_ind, gs_ind)

        gs_bipart = e_values(gs_ind)
        gs_vec_bipart = e_vecs(:,gs_ind)

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
        write(iunit_write, '(G25.17)', advance = 'yes') log2real(is_vec_sign_coherent(gs_vec_bipart))

        close(iunit_write)

        print *, "size(hilbert_space): ", shape(hilbert_space)
        print *, "n_perms: ", n_perms
        print *, "original gs energy: ", gs_orig
        print *, "original HF coeff: ", hf_orig
        print *, "original gs vec non-zeros: ", get_num_non_zeros(gs_vec)
        print *, "original gs vec sign coherent?:", is_vec_sign_coherent(gs_vec)

        print *, "bipartite gs energy:", gs_bipart
        print *, "bipartite hf coeff: ", hf_bipart
        print *, "original gs vec non-zeros: ", get_num_non_zeros(gs_vec_bipart)
        print *, "original gs vec sign coherent?:", is_vec_sign_coherent(gs_vec_bipart)

        call stop_all("here", "now")
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

            call eig(hamil, e_values, e_vecs)

            ! n_orbs - 1 because I do not consider cycles

            ! then I need to sort by the eigenvalues to get the GS
            ! and store the largest coeff. from the GS
            call store_gs_hf_coeff(e_values, e_vecs, hf_coeff, hf_ind, gs_ind)
            gs_vec = e_vecs(:, gs_ind)

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

        iunit_write = get_free_unit()
        open(iunit_write, file = 'all-data', status = 'replace', action = 'write')
        write(iunit_write, *) "# hf-coeff     off-diag-sum    abs-off-diag-sum    &
            &n-diff-sign    sum-pos     sum-neg     sum-3-cycle     sum-n-cycle &
            &   sum-non-zeros   non-zeros-in-gs     sign-coherent-gs"
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
            write(iunit_write, '(G25.17)',advance = 'yes') sign_coherent_gs(i)
        end do
        close(iunit_write)
        call stop_all("here", "for now")

    end subroutine exact_study

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


    subroutine print_vec(vec, filename)
        real(dp), intent(in) :: vec(:)
        character(*), intent(in), optional :: filename

        integer :: iunit, i

        if (present(filename)) then
            iunit = get_free_unit()
            open(iunit, file = filename, status = 'replace', action = 'write')

            do i = 1, size(vec,1)
                write(iunit, *) vec(i)
            end do

            close(iunit)
        else
            do i = 1, size(vec,1)
                print *, vec(i)
            end do
        end if

    end subroutine print_vec

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

    subroutine store_gs_hf_coeff(e_values, e_vecs, hf_coeff, hf_ind, gs_ind)
        real(dp), intent(in) :: e_values(:), e_vecs(:,:)
        real(dp), intent(out) :: hf_coeff
        integer, intent(out) :: hf_ind, gs_ind

        real(dp) :: gs_vec(size(e_values))

        gs_ind = minloc(e_values,1)

        gs_vec = abs(e_vecs(:,gs_ind))

        hf_ind = maxloc(gs_vec,1)
        hf_coeff = gs_vec(hf_ind)

    end subroutine store_gs_hf_coeff

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
        ! it should be fine if i always return the exhange j when the
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
