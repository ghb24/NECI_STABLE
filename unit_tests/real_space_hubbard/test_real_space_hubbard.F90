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
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine real_space_hubbard_test_driver() 
        ! this is the main function which calls all the other tests 
        
        call test_init_lattice() 

    end subroutine real_space_hubbard_test_driver

    subroutine test_init_lattice() 
        ! the problem in this case is, that it tests global variables.. 
        ! and we have the problem with the names of the function.. 
        ! i guess i need to use simons workaround here too.. 

        ! or i use the provided routine in fruit: 
        ! but there should be a nice way provided by fruit already or?? 
        call set_case_name("test_init_lattice") 


        lattice_type = "chain"
        call init_lattice

        call assert_equals(n_dim, 1)

        call set_case_name('_not_set_')

    end subroutine test_init_lattice

end program test_real_space_hubbard

