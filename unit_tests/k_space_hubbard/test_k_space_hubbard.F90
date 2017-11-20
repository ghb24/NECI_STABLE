#include "macros.h"

! i guess for a new module it would be nice to have one module which tests 
! all the functions and programs of a module .. check that out.. 

! all the global variables in NECI make it a bit tough to do unit-tests

! no.. i guess, atleast in my framework it would be better to write one 
! unit test module which tests all the functionality of a other module 
program test_k_space_hubbard 
    
    use k_space_hubbard
    use constants 
    use fruit 
    use lattice_mod, only: lat

    implicit none 

    integer :: failed_count 


    call init_fruit()
    ! run the test-driver 
    call k_space_hubbard_test_driver()
    call fruit_summary()
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine k_space_hubbard_test_driver() 
        ! this is the main function which calls all the other tests 
       
    end subroutine k_space_hubbard_test_driver

end program test_k_space_hubbard
