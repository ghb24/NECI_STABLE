#include "macros.h"

! data module of the real-time implementation of FCIQMC

module real_time_data
    
    use FciMCData, only: perturbation
    implicit none

    ! number of copies to start the real-time simulation from
    integer :: n_real_time_copies, cnt_real_time_copies
    logical :: t_prepare_real_time

    ! for the actual calculation: 
    ! global flag indicating real-time calculation
    logical :: t_real_time_fciqmc

    type(perturbation), allocatable :: gf_type(:)

end module real_time_data
