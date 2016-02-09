#include "macros.h"

! data module of the real-time implementation of FCIQMC

module real_time_data
    
    implicit none

    ! number of copies to start the real-time simulation from
    integer :: n_real_time_copies, cnt_real_time_copies
    logical :: t_prepare_real_time

end module real_time_data
