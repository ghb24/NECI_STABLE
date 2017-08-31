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
    use SystemData, only: lattice_type

    implicit none 

    integer :: failed_count 

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
        call run_test_case(determine_optimal_time_step_test, "determine_optimal_time_step_test")
        call run_test_case(gen_excit_rs_hubbard_test, "gen_excit_rs_hubbard_test")
        call run_test_case(init_real_space_hubbard_test, "init_real_space_hubbard_test")

    end subroutine real_space_hubbard_test_driver

    subroutine init_real_space_hubbard_test
        use SystemData, only: lattice_type, length_x, length_y, nbasis, nel, &
                              bhub, uhub, ecore
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
        call assert_equals(0.25, tau)

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
        use SystemData, only: nel 
        use bit_rep_data, only: niftot
        use Detbitops, only: encodebitdet
        use constants, only: n_int, dp
        use fcimcdata, only: excit_gen_store_type
        
        integer, allocatable :: nI(:), nJ(:)
        integer(n_int), allocatable :: ilutI(:), ilutJ(:)
        integer :: ex(2,2), ic 
        logical :: tpar 
        real(dp) :: pgen 
        real(dp) :: hel
        type(excit_gen_store_type) :: store
        logical :: found_all, t_1, t_2, t_3, t_4
        
        print *, "" 
        print *, "testing gen_excit_rs_hubbard"
        lat => lattice('chain', 4, 1, .false., .false.)
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
        t_1 = .false. 
        t_2 = .false. 

        ! try to get all the excitations here. 
        do while(.not. found_all)
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_1) then 
                t_1 = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(3, ex(2,1))
                call assert_equals(6, ilutJ(0))
                call assert_true(tpar)
                call assert_equals(0.5_dp, pgen)
            end if
            if (all(nJ == [1,4]) .and. .not. t_2) then 
                t_2 = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(4, ex(2,1))
                call assert_equals(9, ilutJ(0))
                call assert_true(.not.tpar)
                call assert_equals(0.5_dp, pgen)
            end if
            found_all = (t_1 .and. t_2)
        end do

        nI = [3,4]
        call encodebitdet(nI, ilutI)

        found_all = .false. 
        t_1 = .false. 
        t_2 = .false. 
        t_3 = .false. 
        t_4 = .false. 

        ! try to get all the excitations here. 
        do while(.not. found_all)
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_1) then 
                t_1 = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(4, ex(1,1))
                call assert_equals(2, ex(2,1))
                call assert_equals(6, ilutJ(0))
                call assert_true(tpar)
                call assert_equals(0.25, pgen)
            end if
            if (all(nJ == [1,4]) .and. .not. t_2) then 
                t_2 = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(3, ex(1,1))
                call assert_equals(1, ex(2,1))
                call assert_equals(9, ilutJ(0))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            if (all(nJ == [3,6]) .and. .not. t_3) then 
                t_3 = .true.
                call assert_equals([3,6], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(4, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(36, ilutJ(0))
                call assert_true(.not.tpar)
                call assert_equals(0.25, pgen)
            end if
            if (all(nJ == [4,5]) .and. .not. t_4) then 
                t_4 = .true. 
                call assert_equals([4,5], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(3, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(24, ilutJ(0))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            found_all = (t_1 .and. t_2 .and. t_3 .and. t_4)
        end do

        ! and also try it on a periodic chain 
        lat => lattice('chain', 3, 1, .true., .true.)

        nI = [1,2]
        call encodebitdet(nI, ilutI) 

        found_all = .false. 
        t_1 = .false. 
        t_2 = .false. 
        t_3 = .false. 
        t_4 = .false.

        ! try to get all the excitations here. 
        do while(.not. found_all)
            call gen_excit_rs_hubbard(nI, ilutI, nJ, ilutJ, 1, ic, ex, tpar, pgen,  & 
                hel, store)

            if (all(nJ == [2,3]) .and. .not. t_1) then 
                t_1 = .true.
                call assert_equals([2,3], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(3, ex(2,1))
                call assert_equals(6, ilutJ(0))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [1,4]) .and. .not. t_2) then 
                t_2 = .true. 
                call assert_equals([1,4], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(4, ex(2,1))
                call assert_equals(9, ilutJ(0))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            if (all(nJ == [1,6]) .and. .not. t_3) then 
                t_3 = .true.
                call assert_equals([1,6], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(2, ex(1,1))
                call assert_equals(6, ex(2,1))
                call assert_equals(33, ilutJ(0))
                call assert_true(.not.tpar)
                call assert_equals(0.25_dp, pgen)
            end if
            if (all(nJ == [2,5]) .and. .not. t_4) then 
                t_4 = .true. 
                call assert_equals([2,5], nJ, 2)
                call assert_equals(1, ic) 
                call assert_equals(1, ex(1,1))
                call assert_equals(5, ex(2,1))
                call assert_equals(18, ilutJ(0))
                call assert_true(tpar)
                call assert_equals(0.25_dp, pgen)
            end if

            found_all = (t_1 .and. t_2 .and. t_3 .and. t_4)
        end do

        nel = -1 
        call lattice_deconstructor(lat)

    end subroutine gen_excit_rs_hubbard_test

    subroutine determine_optimal_time_step_test 
        use SystemData, only: nel, nOccAlpha, nOccBeta, uhub, bhub
        use lattice_mod, only: lattice

        real(dp) :: time, time_death
        class(lattice), pointer :: ptr

        print *, ""
        print *, "testing determine_optimal_time_step"

        lat => lattice('chain', 4, 1, .true., .true.)

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
        
        nel = -1
        nOccAlpha = -1
        nOccBeta = -1
        uhub = 0
        bhub = 0
        
    end subroutine determine_optimal_time_step_test 

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
        class(lattice), pointer :: ptr
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

        ptr => lattice('chain', 2, 1, .false., .false.)
        call init_tmat(ptr)

        call assert_equals(0.0_dp, get_helement([1,2],[1,2],0))
        call assert_equals(0.0_dp, get_helement([1,2],[1,2]))

        uhub = 1.0
        call assert_equals(1.0_dp, get_helement([1,2],[1,2],0))
        call assert_equals(1.0_dp, get_helement([1,2],[1,2]))

        call assert_equals(1.0_dp, get_helement([1,2],[1,4],1))
        call assert_equals(-1.0_dp, get_helement([1,2],[2,3],1))

        call assert_equals(1.0_dp, get_helement([1,2],[1,4]))
        call assert_equals(-1.0_dp, get_helement([1,2],[2,3]))

        nbasis = 6 
        deallocate(g1)
        allocate(g1(nbasis))
        g1([1,3,5])%ms = -1 
        g1([2,4,6])%ms = 1 

        ptr => lattice('chain', 3, 1, .true., .true.)
        call init_tmat(ptr)

        call assert_equals(0.0, get_helement([1,3],[1,3],0))
        call assert_equals(0.0, get_helement([1,3],[1,3]))

        call assert_equals(1.0, get_helement([1,3],[1,5],1))
        call assert_equals(-1.0, get_helement([1,3],[3,5],1))

        call assert_equals(1.0, get_helement([1,3],[1,5]))
        call assert_equals(-1.0, get_helement([1,3],[3,5]))

        niftot = 0
        nifd = 0
        allocate(ilutI(0:niftot))
        allocate(ilutJ(0:niftot))

        call encodebitdet([1,2],ilutI)
        call encodebitdet([1,2],ilutJ)

        call assert_equals(1.0, get_helement([1,2],[1,2],ilutI,ilutJ))

        call encodebitdet([2,4],ilutJ)
        call assert_equals(0.0, get_helement([1,2],[2,4],ilutI,ilutJ))

        call encodebitdet([2,3],ilutJ)
        call assert_equals(-1.0, get_helement([1,2],[2,3],ilutI,ilutJ))

        call encodebitdet([1,4],ilutJ)
        call assert_equals(1.0, get_helement([1,2],[1,4],ilutJ,ilutJ))

        nel = 4 
        call assert_equals(2.0, get_helement([1,2,3,4],[1,2,3,4]))
        call assert_equals(2.0, get_helement([1,2,3,4],[1,2,3,4],0))

        call assert_equals(1.0, get_helement([1,2,3,4],[1,2,3,6]))
        call assert_equals(1.0, get_helement([1,2,3,4],[1,2,3,6],1))

        call assert_equals(1.0, get_helement([1,2,3,4],[1,3,4,6],1))
        call assert_equals(1.0, get_helement([1,2,3,4],[1,3,4,6]))
        call assert_equals(-1.0, get_helement([1,2,3,4],[2,3,4,5],1))
        call assert_equals(-1.0, get_helement([1,2,3,4],[2,3,4,5]))

        nel = -1 
        nbasis = -1 
        nullify(get_umat_el)
        deallocate(g1) 
        deallocate(tmat2d)

    end subroutine get_helement_test

    subroutine init_tmat_test
        use lattice_mod, only: lattice
        use SystemData, only: nbasis, bhub
        use OneEInts, only: tmat2d

        class(lattice), pointer :: ptr

        print *, "" 
        print *, "testing: init_tmat" 

        print *, "for a 4 site, non-periodic chain geometry"
        nbasis = 8 
        bhub = 1.0
        ptr => lattice('chain', 4, 1, .false., .false.)

        call init_tmat(ptr)

        call assert_equals([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], tmat2d(1,:), 8)
        call assert_equals([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], tmat2d(2,:), 8)
        call assert_equals([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], tmat2d(3,:), 8)
        call assert_equals([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], tmat2d(4,:), 8)
        call assert_equals([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], tmat2d(7,:), 8)
        call assert_equals([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], tmat2d(8,:), 8)

        print *, ""
        print *, "for a 4 site periodic chain geometry: "
        ptr => lattice('chain', 4, 1, .true., .true.)

        call init_tmat(ptr)

        call assert_equals([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0], tmat2d(1,:), 8)
        call assert_equals([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], tmat2d(2,:), 8)
        call assert_equals([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], tmat2d(3,:), 8)
        call assert_equals([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], tmat2d(4,:), 8)
        call assert_equals([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0], tmat2d(7,:), 8)
        call assert_equals([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], tmat2d(8,:), 8)

        ! todo: more tests for other lattices later!

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

        call assert_equals(1.0, get_umat_el_hub(1,1,1,1))
        call assert_equals(1.0, get_umat_el_hub(2,2,2,2))
        call assert_equals(1.0, get_umat_el_hub(3,3,3,3))
        call assert_equals(1.0, get_umat_el_hub(4,4,4,4))
        call assert_equals(0.0, get_umat_el_hub(1,2,3,4))
        call assert_equals(0.0, get_umat_el_hub(1,1,1,4))
        call assert_equals(0.0, get_umat_el_hub(2,2,3,4))
        call assert_equals(0.0, get_umat_el_hub(3,2,3,4))

        uhub = 0.0
        nbasis = -1

    end subroutine



end program test_real_space_hubbard

