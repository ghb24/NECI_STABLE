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


