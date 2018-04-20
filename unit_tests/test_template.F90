! replace _template by whatever one needs 
program test_template

    use fruit
    ! use template

    implicit none 

    integer :: failed_count 

    call init_fruit() 
    call template_test_driver() 
    call fruit_summary() 
    call fruit_finalize()

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine template_test_driver() 

        call run_test_case(first_test, "first_test")

    end subroutine template_test_driver

    subroutine first_test() 

        call assert_true(.false.)
        call assert_equals(0,1)

    end subroutine first_test

end program test_template

