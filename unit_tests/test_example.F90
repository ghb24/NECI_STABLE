#include "test_macros.h"
module example_test
    use fruit
    implicit none
contains

    subroutine test_example
        integer :: res
        res = 13
        call assert_equals(res, 13)
        call assert_equals(res, 42)
    end subroutine

    subroutine example_test_driver
        TEST(test_example)
    end subroutine

end module
