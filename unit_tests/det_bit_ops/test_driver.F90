#include "test_macros.h"
module det_bit_ops_tests
    use fruit
    implicit none
contains

    subroutine det_bit_ops_drive_tests
        use countbits_tests
        TEST(countbits_drive_tests)
    end subroutine

end module
