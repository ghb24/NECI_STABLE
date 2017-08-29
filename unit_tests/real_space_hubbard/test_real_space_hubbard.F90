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
!         call run_test_case(test_init_lattice, "test_init_lattice")

    end subroutine real_space_hubbard_test_driver

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
        nel = -1 


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

