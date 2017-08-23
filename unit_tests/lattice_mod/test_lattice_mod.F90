#include "macros.h" 

! test my lattice module here, if the initializations and getter and 
! setter functions work as intended 

program test_lattice_mod 

    use lattice_mod 
    use fruit 

    implicit none 

    integer :: failed_count 

    call init_fruit() 
    ! run my tests: 
    call lattice_mod_test_driver
    call fruit_summary 
    call fruit_finalize

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine lattice_mod_test_driver

        call run_test_case(test_init_lattice_chain, "test_init_lattice_chain")
        call run_test_case(test_init_lattice_star, "test_init_lattice_star")
!         call run_test_case(test_init_lattice_aim_chain, "test_init_lattice_aim_chain")

    end subroutine lattice_mod_test_driver 

    subroutine test_init_lattice_star
        class(lattice), pointer :: ptr

        integer :: i

        print *, ""
        print *, "initialize a 1 site 'star' geometry lattice: "
        ptr => lattice('star', 1, 1, .false., .false.)

        call assert_equals(1, ptr%get_ndim())
        call assert_equals(1, ptr%get_length())
        call assert_equals(1, ptr%get_nsites())
        call assert_equals(0, ptr%get_nconnect_max())
        call assert_true(.not. ptr%is_periodic())

        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(ptr%get_neighbors(1), [-1], 1)

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize a 2-site 'star' geometry lattice: "
        ptr => lattice('star', 2, 1, .false., .false.)
        call assert_equals(1, ptr%get_ndim())
        call assert_equals(1, ptr%get_length())
        call assert_equals(2, ptr%get_nsites())
        call assert_equals(1, ptr%get_nconnect_max())
        call assert_true(.not. ptr%is_periodic())

        call assert_equals(ptr%get_site_index(1), 1) 
        call assert_equals(ptr%get_site_index(2), 2) 
        call assert_equals(ptr%get_neighbors(1), [2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1], size(ptr%get_neighbors(2)))

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize a 100 site 'star' geometry lattice: " 
        ptr => lattice('star', 100, 1, .false., .false.) 
        call assert_equals(1, ptr%get_ndim())
        call assert_equals(1, ptr%get_length())
        call assert_equals(100, ptr%get_nsites())
        call assert_equals(99, ptr%get_nconnect_max())
        call assert_true(.not. ptr%is_periodic())

        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(100, ptr%get_site_index(100))

        call assert_equals([1], ptr%get_neighbors(2), 1)
        call assert_equals([1], ptr%get_neighbors(100), 1)
        call assert_equals( [(i, i = 2, 100)], ptr%get_neighbors(1),99)

    end subroutine test_init_lattice_star
! 
!     subroutine test_init_lattice_aim_chain
!         class(lattice), pointer :: ptr 
! 
!         print *, "initialize a 'aim-chain' lattice with 1 bath site: "
!         ptr => lattice('aim-chain')
! 
!         call assert_true(.false.)
! 
!     end subroutine test_init_lattice_aim_chain

    subroutine test_init_lattice_chain 
        ! test the specific initializers 
        ! implicitly test the setter and getter functions or? 
        class(lattice), pointer :: ptr 

        print *, "initialize a periodic one-site 'chain' lattice: " 
        ptr => lattice('chain', 1, 1, .true., .true.) 

        call assert_equals(ptr%get_ndim(), 1) 
        call assert_equals(ptr%get_nsites(), 1) 
        call assert_equals(ptr%get_length(), 1)
        call assert_true(ptr%is_periodic()) 
        ! hm.. for this case the n_connect_max is not 2. since it is a special 
        ! case.. for the 1 site and the 2 site non-periodic it is only 1 or 0
        ! but this special cases shouldnt matter.. or? 
        ! leave it for now, so i am reminded that i might have to adjust that!
        call assert_equals(0, ptr%get_nconnect_max())

        ! todo: it would be really fancy to overload the round brackets 
        ! to give the get function of the index value of our lattice!:
        ! for a lattice of length 1 we need some special initialization.. 
        ! thats good to have these edge cases!
        ! i need a public getter for the site indices.. 
        ! do i want to put it into the site type or in the lattice type? 
        ! ptr%get_index(1) or ptr%sites(1)%get_index 
        ! in the first i could check if the index is too high and i would 
        ! not have to make so much public..
        call assert_equals(ptr%get_site_index(1), 1)
        ! and i want to have a get_neighbors routine
        ! apparently there is no assert equal for vectors of ints?? 
        ! thats BS! there is, but one has to give the additional number 
        ! of elements input!
        call assert_equals(ptr%get_neighbors(1), [-1], 1)

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

        print *, ""
        print *, "initialize a non-periodic two-site 'chain' lattice: " 
        ptr => lattice('chain', 2, 1, .false., .false.) 
        call assert_equals(ptr%get_ndim(), 1) 
        call assert_equals(ptr%get_nsites(), 2)
        call assert_equals(ptr%get_length(), 2) 
        call assert_true(.not. ptr%is_periodic()) 
        call assert_equals(1, ptr%get_nconnect_max())

        call assert_equals(ptr%get_site_index(1), 1) 
        call assert_equals(ptr%get_site_index(2), 2) 
        call assert_equals(ptr%get_neighbors(1), [2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1], size(ptr%get_neighbors(2)))

        call lattice_deconstructor(ptr) 

        call assert_true(.not. associated(ptr))

        print *, "" 
        print *, "initialize a periodic 100 site 'chain' lattice: "
        ptr => lattice('chain', 0, 100, .true., .true.)
        call assert_equals(ptr%get_ndim(), 1) 
        call assert_equals(ptr%get_nsites(), 100)
        call assert_equals(ptr%get_length(), 100) 
        call assert_true(ptr%is_periodic()) 
        call assert_equals(2, ptr%get_nconnect_max())

        call assert_equals(ptr%get_site_index(1), 1) 
        call assert_equals(ptr%get_site_index(100), 100) 
        call assert_equals(ptr%get_site_index(50), 50) 
        call assert_equals(ptr%get_neighbors(1), [100, 2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1,3], size(ptr%get_neighbors(2)))
        call assert_equals(ptr%get_neighbors(100), [99,1], size(ptr%get_neighbors(2)))

        ! i actually do not need to have the common lattice type or? 
        ! when i want to test a chain i could just use a chain or? 

    end subroutine test_init_lattice_chain

end program test_lattice_mod
