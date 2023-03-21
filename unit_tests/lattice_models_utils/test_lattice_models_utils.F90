#include "macros.h"
! replace _template by whatever one needs
program test_lattice_models_utils

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case, assert_true, assert_equals
    use lattice_models_utils, only: find_minority_spin, pick_spin_par_elecs, &
        pick_three_opp_elecs, pick_spin_opp_elecs, make_ilutJ, &
        get_orb_from_kpoints, get_ispn, get_occ_neighbors, &
        get_spin_density_neighbors, find_elec_in_ni, &
        get_orb_from_kpoints_three, create_all_open_shell_dets, &
        get_spin_opp_neighbors, create_one_spin_basis, calc_n_double, &
        create_neel_state_chain, create_neel_state, &
        pick_from_cum_list, combine_spin_basis, set_alpha_beta_spins, &
        right_most_zero
    use constants, only: dp, n_int
    use lattice_mod, only: lat
    use dsfmt_interface, only: dsfmt_init
    use SystemData, only: t_k_space_hubbard
    use lattice_mod, only: lattice
    use k_space_hubbard, only: setup_nbasismax, setup_g1, setup_kPointToBasisFn
    use unit_test_helpers, only: setup_arr_brr
    use SystemData, only: bhub, nn_bhub, nEl, nBasis, nOccBeta, nOccAlpha

    implicit none

    integer :: failed_count

    call init_fruit()
    call dsfmt_init(0)
    call lattice_models_utils_test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

    subroutine lattice_models_utils_test_driver
        call run_test_case(create_neel_state_chain_test, "create_neel_state_chain_test")
        call run_test_case(create_neel_state_test, "create_neel_state_test")
        call run_test_case(calc_n_double_test, "calc_n_double_test")
        call run_test_case(right_most_zero_test, "right_most_zero_test")
        call run_test_case(create_one_spin_basis_test, "create_one_spin_basis_test")
        call run_test_case(set_alpha_beta_spins_test, "set_beta_spins_test")
        call run_test_case(combine_spin_basis_test, "combine_spin_basis_test")
        call run_test_case(create_all_open_shell_dets_test, "create_all_open_shell_dets_test")
        call run_test_case(get_spin_opp_neighbors_test, "get_spin_density_neighbors")
        call run_test_case(find_elec_in_ni_test, "find_elec_in_ni_test")
        call run_test_case(get_occ_neighbors_test, "get_occ_neighbors_test")
        call run_test_case(get_spin_density_neighbors_test, "get_spin_density_neighbors_test")
        call run_test_case(test_get_ispn, "test_get_ispn")
        call run_test_case(get_orb_from_kpoints_test, "get_orb_from_kpoints_test")
        call run_test_case(make_ilutJ_test, "make_ilutJ_test")
        call run_test_case(pick_spin_opp_elecs_test, "pick_spin_opp_elecs_test")
        call run_test_case(pick_from_cum_list_test, "pick_from_cum_list_test")
        call run_test_case(pick_three_opp_elecs_test, "pick_three_opp_elecs_test")
        call run_test_case(pick_spin_par_elecs_test, "pick_spin_par_elecs_test")
        call run_test_case(find_minority_spin_test, "find_minority_spin_test")
        call run_test_case(get_orb_from_kpoints_three_test, "get_orb_from_kpoints_three_test")
    end subroutine lattice_models_utils_test_driver

    subroutine find_minority_spin_test

        print *, ""
        print *, "testing: find_minority_spin"
        call assert_equals(1, find_minority_spin([1,2,4]))
        call assert_equals(3, find_minority_spin([3,2,4]))

        call assert_equals(2, find_minority_spin([1,2,3]))
        call assert_equals(2, find_minority_spin([3,2,5]))

    end subroutine find_minority_spin_test

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

    subroutine get_orb_from_kpoints_test
        use SystemData, only: G1, nBasis, nel
        use symexcitdatamod, only: KPointToBasisFn

        nbasis = 4
        nel = 2

        if (associated(g1)) deallocate(g1)
        if (allocated(KPointToBasisFn)) deallocate(KPointToBasisFn)
        allocate(G1(nbasis))
        allocate(KPointToBasisFn(-1:2,-1:2,-1:1,2))

        print *, ""
        print *, "testing: get_orb_from_kpoints: "
        print *, "with necessary global data: "
        print *, "get_ispn"
        print *, "G1 "
        print *, "kpointtobasisfn"
        print *, "nBasis: ", nBasis
        print *, "nel: ", nel

        ! i have to setup G1 and kpointtobasisfn correctly
        G1(1)%k = [1,0,0]
        G1(2)%k = [0,1,0]
        G1(3)%k = [0,0,1]

        KPointToBasisFn(1,1,-1,2) = 3
        KPointToBasisFn(2,0,-1,1) = 4
        kpointtobasisfn(-1,2,0,2) = 5

        call assert_equals(3, get_orb_from_kpoints(1,2,3))
        call assert_equals(4, get_orb_from_kpoints(1,1,3))
        call assert_equals(5, get_orb_from_kpoints(2,2,1))

        nBasis = - 1
        nel = -1

        deallocate(G1)
        deallocate(KPointToBasisFn)

    end subroutine get_orb_from_kpoints_test

    subroutine make_ilutJ_test
        use constants, only: n_int
        use bit_rep_data, only: niftot
        use SystemData, only: nbasis

        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2)

        niftot = 0
        nBasis = 2

        allocate(ilut(0:niftot))

        print *, ""
        print *, "testing: make_ilutJ"
        print *, "with necessary global data:"
        print *, "niftot: ", niftot
        print *, "nBasis: ", nBasis

        ilut = 0_n_int
        ex(1,:) = [1,0]
        ex(2,:) = [1,0]

        ! i guess this test could be dependend if the machine is little or
        ! big endian..
        ! ..01 = 1
        ! for the intel compilers i have to cast this to a "normal" integer
        ! so it finds the correct assert_equals routine..
        call assert_equals([1], int(make_ilutJ(ilut, ex, 1)), 1)
        ex(2,1) = 2
        ! ..10 = 2
        call assert_equals([2], int(make_ilutJ(ilut, ex, 1)), 1)

        ex(1,2) = 1
        ex(2,2) = 1
        ! ..11 = 3
        call assert_equals([3], int(make_ilutJ(ilut, ex, 2)), 1)

        niftot = -1
        nBasis = -1

    end subroutine make_ilutJ_test

    subroutine test_get_ispn

        print *, ""
        print *, "testing: get_ispn()"

        call assert_equals(1, get_ispn([1,1]))
        call assert_equals(2, get_ispn([1,2]))
        call assert_equals(3, get_ispn([2,2]))

    end subroutine test_get_ispn

    subroutine get_occ_neighbors_test
        use bit_rep_data, only: NIfTot
        use Detbitops, only: EncodeBitDet
        use SystemData, only: nel
        use lattice_mod, only: lattice
        use constants, only: n_int

        integer(n_int), allocatable :: ilut(:)

        NIfTot = 0
        allocate(ilut(0:NIfTot))

        nel = 3
        lat => lattice('chain',4, 1, 1, .true., .true., .true.)

        call EncodeBitDet([1,2,5], ilut)

        print *, ""
        print *, "testing: get_occ_neighbors "
        call assert_equals(0.0_dp, get_occ_neighbors(ilut,1))
        call assert_equals(3.0_dp, get_occ_neighbors(ilut,2))
        call assert_equals(0.0_dp, get_occ_neighbors(ilut,3))
        call assert_equals(3.0_dp, get_occ_neighbors(ilut,4))

        call EncodeBitDet([1,2,3], ilut)
        call assert_equals(1.0_dp, get_occ_neighbors(ilut,1))
        call assert_equals(2.0_dp, get_occ_neighbors(ilut,2))
        call assert_equals(1.0_dp, get_occ_neighbors(ilut,3))

        NIfTot = -1
        nel = -1

    end subroutine get_occ_neighbors_test

    subroutine get_spin_density_neighbors_test
        use bit_rep_data, only: niftot
        use Detbitops, only: EncodeBitDet
        use SystemData, only: nel
        use lattice_mod, only: lattice
        use constants, only: n_int

        integer(n_int), allocatable :: ilut(:)

        print *, ""
        print *, " testing: get_spin_density_neighbors "

        NIfTot = 0
        lat => lattice('chain',4,1,1,.true.,.true.,.true.)
        allocate(ilut(0:NIfTot))

        nel = 3
        call encodebitdet([1,2,3], ilut)

        call assert_equals(-0.5_dp, get_spin_density_neighbors(ilut,1))
        call assert_equals(0.0_dp, get_spin_density_neighbors(ilut,2))
        call assert_equals(-0.5_dp, get_spin_density_neighbors(ilut,3))
        call assert_equals(0.0_dp, get_spin_density_neighbors(ilut,4))

        call encodebitdet([1,4,5], ilut)
        call assert_equals(0.5_dp, get_spin_density_neighbors(ilut,1))
        call assert_equals(-1.0_dp, get_spin_density_neighbors(ilut,2))
        call assert_equals(0.5_dp, get_spin_density_neighbors(ilut,3))
        call assert_equals(-1.0_dp, get_spin_density_neighbors(ilut,4))

        nel = -1
        NIfTot = -1

    end subroutine get_spin_density_neighbors_test

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


    subroutine get_orb_from_kpoints_three_test

        print *, ""
        print *, "testing: get_orb_from_kpoints_three: "
        nel = 4
        nbasis = 8
        bhub = -1.0_dp
        nn_bhub = 0.0_dp
        t_k_space_hubbard = .true.
        lat => lattice('chain', 4,1,1,.true.,.true.,.true.,'k-space')

        call setup_nbasismax(lat)
        call setup_arr_brr(lat)
        call setup_g1(lat)
        call setup_kPointToBasisFn(lat)

        call assert_equals(3, get_orb_from_kpoints_three([1,2,3],1,2))
        call assert_equals(5, get_orb_from_kpoints_three([1,2,3],4,5))
        call assert_equals(5, get_orb_from_kpoints_three([1,3,5],1,3))
        call assert_equals(2, get_orb_from_kpoints_three([2,4,6],4,6))

        call assert_equals(7, get_orb_from_kpoints_three([1,2,3],7,8))

        call assert_equals(7, get_orb_from_kpoints_three([3,4,5],6,8))

        t_k_space_hubbard = .false.

    end subroutine get_orb_from_kpoints_three_test


    subroutine create_all_open_shell_dets_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: basis(:,:)

        niftot = 0
        nifd = 0

        print *, ""
        print *, "testing: create_all_open_shell_dets"

        basis = create_all_open_shell_dets(4,2,2)

        call assert_equals(6, size(basis,2))
        call assert_equals(int(b'01011010',n_int),basis(1,1))
        call assert_equals(int(b'01100110',n_int),basis(1,2))
        call assert_equals(int(b'01101001',n_int),basis(1,3))
        call assert_equals(int(b'10100101',n_int),basis(1,6))

        niftot = -1
        nifd = -1

    end subroutine create_all_open_shell_dets_test

    subroutine get_spin_opp_neighbors_test
        use bit_rep_data, only: NIfTot
        use Detbitops, only: EncodeBitDet
        use SystemData, only: nel
        use lattice_mod, only: lattice
        use constants, only: n_int

        integer(n_int), allocatable :: ilut(:)

        print *, ""
        print *, "testing: get_spin_opp_neighbors: "

        NIfTot = 0
        lat => lattice('chain', 4,1,1,.true.,.true.,.true.)

        allocate(ilut(0:NIfTot))
        nel = 3
        call encodebitdet([1,2,3], ilut)

        call assert_equals(1.0_dp, get_spin_opp_neighbors(ilut,2))
        call assert_equals(0.0_dp, get_spin_opp_neighbors(ilut,1))
        call assert_equals(1.0_dp, get_spin_opp_neighbors(ilut,3))
        call assert_equals(1.0_dp, get_spin_opp_neighbors(ilut,4))
        call assert_equals(0.0_dp, get_spin_opp_neighbors(ilut,5))
        call assert_equals(1.0_dp, get_spin_opp_neighbors(ilut,6))
        call assert_equals(1.0_dp, get_spin_opp_neighbors(ilut,7))
        call assert_equals(1.0_dp, get_spin_opp_neighbors(ilut,8))

        nel = -1
        NIfTot = -1
    end subroutine get_spin_opp_neighbors_test


    subroutine combine_spin_basis_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: basis(:,:), spin_basis(:)

        niftot = 0
        nifd = 0

        print *, ""
        print *, "testing: combine_spin_basis"

        basis = combine_spin_basis(4,2,2,6,int([3,5,6,9,10,12],n_int),.false.)

        ! for some really strange reason those b'xx' literal are 128-bit integers
        ! in gfortran..
        call assert_equals(6, size(basis,2))
        call assert_equals(int(b'01011010',n_int),(basis(1,1)))
        call assert_equals(int(b'01100110',n_int),(basis(1,2)))
        call assert_equals(int(b'01101001',n_int),(basis(1,3)))
        call assert_equals(int(b'10100101',n_int),(basis(1,6)))

        spin_basis = create_one_spin_basis(6,3)

        basis = combine_spin_basis(6,3,3,20,spin_basis,.false.)

        call assert_equals(20, size(basis))
        call assert_equals(int(b'010101101010',n_int), (basis(1,1)))
        call assert_equals(int(b'101010010101',n_int), (basis(1,20)))

        basis = combine_spin_basis(4,2,1,12,int([3,5,6,9,10,12],n_int),.false.)

        call assert_equals(12, size(basis))
        call assert_equals(int(b'00011010',n_int), (basis(1,1)))
        call assert_equals(int(b'01001010',n_int), (basis(1,2)))
        call assert_equals(int(b'00100110',n_int), (basis(1,3)))
        call assert_equals(int(b'10100001',n_int), (basis(1,11)))
        call assert_equals(int(b'10100100',n_int), (basis(1,12)))

        basis = combine_spin_basis(4,1,1,12,int([1,2,4,8],n_int),.false.)

        call assert_equals(12, size(basis))
        call assert_equals(int(b'00000110',n_int), (basis(1,1)))
        call assert_equals(int(b'00010010',n_int), (basis(1,2)))
        call assert_equals(int(b'01000010',n_int), (basis(1,3)))
        call assert_equals(int(b'10010000',n_int), (basis(1,12)))

        niftot = -1
        nifd = -1
    end subroutine combine_spin_basis_test


    subroutine set_alpha_beta_spins_test
        use bit_rep_data, only: niftot, nifd

        niftot = 0
        nifd = 0

        print *, ""
        print *, "testing: set_alpha_beta_spins "

        call assert_equals(int(b'101',n_int),(set_alpha_beta_spins(int(b'11',n_int),4,.true.)))
        call assert_equals(int(b'00010001',n_int), (set_alpha_beta_spins(int(b'0101',n_int),4,.true.)))
        call assert_equals(int(b'01000001',n_int), (set_alpha_beta_spins(int(b'1001',n_int),4,.true.)))

        call assert_equals(int(b'1010',n_int),(set_alpha_beta_spins(int(b'11',n_int),4,.false.)))
        call assert_equals(int(b'00100010',n_int), (set_alpha_beta_spins(int(b'0101',n_int),4,.false.)))
        call assert_equals(int(b'10000010',n_int), (set_alpha_beta_spins(int(b'1001',n_int),4,.false.)))

        call assert_equals(int(b'10100000',n_int), (set_alpha_beta_spins(int(b'1100',n_int),4,.false.)))
        call assert_equals(int(b'01010000',n_int), (set_alpha_beta_spins(int(b'1100',n_int),4,.true.)))

        niftot = -1
        nifd = -1
    end subroutine set_alpha_beta_spins_test

    subroutine create_one_spin_basis_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int), allocatable :: alpha(:)

        niftot = 0
        nifd = 0
        print *, ""
        print *, "testing: create_one_spin_basis"

        alpha = create_one_spin_basis(4,2)

        call assert_equals(6, size(alpha))
        call assert_equals(int(b'0011',n_int),(alpha(1)))
        call assert_equals(int(b'0101',n_int),(alpha(2)))
        call assert_equals(int(b'0110',n_int),(alpha(3)))
        call assert_equals(int(b'1001',n_int),(alpha(4)))
        call assert_equals(int(b'1010',n_int),(alpha(5)))
        call assert_equals(int(b'1100',n_int),(alpha(6)))

        alpha = create_one_spin_basis(5,2)

        call assert_equals(10, size(alpha))
        call assert_equals(int(b'00011',n_int), (alpha(1)))
        call assert_equals(int(b'11000',n_int), (alpha(10)))

        niftot = -1
        nifd = -1

    end subroutine create_one_spin_basis_test


    subroutine right_most_zero_test
        use bit_rep_data, only: niftot, nifd

        integer(n_int) :: i

        nifd = 0
        niftot = 0

        print *, ""
        print *, "testing: right_most_zero:"
        i = int(b'1001',n_int)
        call assert_equals(2, right_most_zero(i, 4))

        i = int(b'0011',n_int)
        call assert_equals(3, right_most_zero(i, 4))

        i = int(b'1010',n_int)
        call assert_equals(3, right_most_zero(i, 4))

        i = int(b'1100',n_int)
        call assert_equals(5, right_most_zero(i, 4))

        call assert_equals(4, right_most_zero(int(b'1100',n_int), 3))

        call assert_equals(-1, right_most_zero(0_n_int, 3))

        call assert_equals(5, right_most_zero(huge(0_n_int),4))

        niftot = -1
        nifd = -1
    end subroutine right_most_zero_test


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

end program test_lattice_models_utils
