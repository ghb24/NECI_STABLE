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
        
!         call test_init_lattice() 

        ! or try running it with the provided runner of fruit: 
        ! unnamed: 
!         call run_test_case(test_init_lattice)
!         call run_test_case(test_init_lattice, "test_init_lattice")

    end subroutine real_space_hubbard_test_driver

    subroutine test_init_lattice() 
        ! the problem in this case is, that it tests global variables.. 
        ! and we have the problem with the names of the function.. 
        ! i guess i need to use simons workaround here too.. 

        ! or i use the provided routine in fruit: 
        ! but there should be a nice way provided by fruit already or?? 

        ! i could output stuff here, which only will be printed if 
!         ! the tests fail with the correct flags to ctest!
!         print *, "Testing: init_lattice():"
! 
!         print *, "testing 'chain' setting: "
!         lattice_type = "chain"
!         call init_lattice
!         ! check the default values 
!         call assert_equals(n_dim, 1)
!         call assert_equals(n_connect_max, 2) 
!         call assert_equals(n_sites, length_x) 
!         call assert_equals(length_y, 1) 
! 
!         print *, ""
!         print *, "testing 'square' setting:"
!         lattice_type = "square"
!         call init_lattice() 
!         call assert_equals(n_dim, 2) 
!         call assert_equals(length_x, 1)
!         call assert_equals(length_y, 1) 
!         call assert_equals(n_sites, 2) 
!         call assert_equals(n_connect_max, 4)
! 
        ! but can i pack more then one unit test into this? 
        ! and how do i find the failing tests then? 
        ! i have to write this here as the to be tested config. so either i 
        ! write individual small test routines for each setup or i am 
        ! lazy enough to do more then one test in here and hope that 
        ! fruit provides enough output

    end subroutine test_init_lattice

!     subroutine test_nearest_neighbors
!         ! routine to check the nearest neighbors list for a specific 
!         ! lattice type 
!         ! for the square lattice 
! 
!         ! to make testing and the overall code more readable and 
!         ! maintainable i should not use so many global variables.. 
!         
!         print *, "testing: nearest_neighbors functionality: " 
! 
!         print *, "'chain' setting (periodic)" 
!         ! for the chain setting we order the orbitals left to write 
!         ! of course. so we only have to take boundary conditions into 
!         ! account and we still have to decide if we want to store this 
!         ! info in spatial or spin-orbitals?! 
!         ! the use of this will be to pick an empty orbital after an 
!         ! electron(in spin-orbitals indication) was picked (for single 
!         ! excitations!) 
!         ! for a single impurity site it would also be beneficial, if 
!         ! we have the spin-orbital neighbors. 
!         ! for double excitations? this only happens for 
!         ! multiple impurity sites.. but there the double excitations 
!         ! will be treated differently anyway.. since we want to 
!         ! distinguish between bath and impurity sites there 
!         ! -> so store it in spin-orbitals 
!         ! and: lets store the neighbor list in an ordered manner! 
! 
!         ! setup some sites and test neighbors
!         ! 2 sites
!         ! and i decided to test both the setup and the get routine in 
!         ! one go, since we do actually need both here anyway
!         t_periodic = .true.
!         lattice_type = 'chain' 
!         length_x = 2 
!         ! i should call as little other functions before i test this 
!         ! unit.. maybe only set_nearest_neighbors?
! 
!         call create_neighbor_list(n_sit = 2, lat_type= 'chain', len_x = 2,  & 
!             len_y = 1, t_periodic = .true.)
! 
!         ! spin-orbital 1 then only has 3 as neighbor
!         call assert_equals(get_nearest_neighbors(1), [3])
!         call assert_equals(get_nearest_neighbors(2), [4])
!         call assert_equals(get_nearest_neighbors(3), [1])
! 
!     end subroutine test_nearest_neighbors

end program test_real_space_hubbard

