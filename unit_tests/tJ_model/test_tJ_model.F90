#include "macros.h"

program test_tJ_model

    use SystemData, only: t_tJ_model, t_heisenberg_model, exchange_j, t_lattice_model
    use tJ_model
    use fruit
    use lattice_mod, only: lat

    implicit none

    integer :: failed_count

    t_tJ_model = .true.
    t_lattice_model = .true.

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

    subroutine init_tJ_model_test
        use SystemData, only: lattice_type, length_x, length_y, nbasis, nel, nbasis
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
        integer :: ex(2,2), ic
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
