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

        call run_test_case(test_init_lattice, "test_init_lattice")

    end subroutine lattice_mod_test_driver 

    subroutine test_init_lattice 
        ! test the specific initializers 
        ! implicitly test the setter and getter functions or? 
        class(lattice), pointer :: ptr 

        print *, "initialize a periodic one-site 'chain' lattice: " 
        ptr => lattice('chain', 1, 1, .true., .true.) 

        call assert_equals(ptr%get_ndim(), 1) 
        call assert_equals(ptr%get_nsites(), 1) 
        call assert_equals(ptr%get_length(), 1)
        call assert_true(ptr%is_periodic()) 
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

        call assert_equals(ptr%get_site_index(1), 1) 
        call assert_equals(ptr%get_site_index(100), 100) 
        call assert_equals(ptr%get_site_index(50), 50) 
        call assert_equals(ptr%get_neighbors(1), [100, 2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1,3], size(ptr%get_neighbors(2)))
        call assert_equals(ptr%get_neighbors(100), [99,1], size(ptr%get_neighbors(2)))

        ! i actually do not need to have the common lattice type or? 
        ! when i want to test a chain i could just use a chain or? 

    end subroutine test_init_lattice

end program test_lattice_mod
