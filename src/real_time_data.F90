#include "macros.h"

! data module of the real-time implementation of FCIQMC

module real_time_data
    
    use constants, only: dp, int64
    use FciMCData, only: perturbation
    implicit none

    ! number of copies to start the real-time simulation from
    integer :: n_real_time_copies, cnt_real_time_copies
    logical :: t_prepare_real_time

    ! for the actual calculation: 
    ! global flag indicating real-time calculation
    logical :: t_real_time_fciqmc

    ! store the type of greensfunction as an integer, or maybe later on even
    ! create a type like nick did, to store all the general information like 
    ! timestep, which orbitals are manipulated etc. in it 
    integer :: gf_type      ! -1 indicates a the lesser, 1 the greater GF

    ! need a global variable for the overlap <y(0)|y(t)> 
    ! determined by the max. cycle
!     real(dp), allocatable :: gf_overlap(:)
    ! for tests now only make it of length 1
    real(dp) :: gf_overlap(1)

    ! need additional info of the original walker number and the number of 
    ! walkers remaining in the perturbed ground state
    integer(int64) :: TotWalkers_pert, TotWalkers_orig


end module real_time_data
