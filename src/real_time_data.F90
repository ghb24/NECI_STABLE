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

    type real_time_type
        ! number of starting vectors = number of parallel mneci calculations
        ! also necessary to have atleast that many popsfile.n
        integer :: n_starting_vec = 1

        ! use a equidistant time-step -> this turns off automated time-step
        ! optimization. but this way i can define a specific end time 
        logical :: t_equidistant_time = .false.
        ! this input should also be provided with a predefined time-step
        real(dp) :: time_step = -1.0_dp
        ! and a end time to stop the simulation afterwards
        real(dp) :: max_time = -1.0_dp

        ! later also store the type of operators and spinorbitals in this 
        ! type! 

        ! store the type of the greensfunction calculated 
        !  1 ... greater GF: creation operator applied! 
        ! -1 ... lesser GF: annihilation operator applied!
        integer :: gf_type = 0

    end type real_time_type

    type(real_time_type) :: real_time_info




        
end module real_time_data
