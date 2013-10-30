interface
    function fn2 (i, j, k, l, fn) result(hel)
        use constants, only: dp
        implicit none
        integer, intent(in) :: i, j, k, l
        HElement_t :: hel
#include "umat_ptr.h"
    end function
end interface
