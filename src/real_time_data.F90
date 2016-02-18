#include "macros.h"

! data module of the real-time implementation of FCIQMC

module real_time_data
    
    use constants, only: dp, int64, n_int
    use FciMCData, only: perturbation, ll_node
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
    real(dp), allocatable :: gf_overlap(:)
    ! also store the norm of the time-evolved wavefunction
    real(dp), allocatable :: wf_norm(:)

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
        ! also store the number of time-steps to be calculated
        integer :: n_time_steps = -1

        ! later also store the type of operators and spinorbitals in this 
        ! type! 

        ! store the type of the greensfunction calculated 
        !  1 ... greater GF: creation operator applied! 
        ! -1 ... lesser GF: annihilation operator applied!
        integer :: gf_type = 0

        ! to reduce the explosive spread of walkers through the 
        ! Hilbert space a small imaginery energy can be introduced in
        ! the Schroedinger equation id/dt y(t) = (H-E0-ie)y(t)
        real(dp) :: damping = 0.0_dp

    end type real_time_type

    type(real_time_type) :: real_time_info

    ! info if the FCIDUMP integrals are complex 
    logical :: t_complex_ints

    ! need a 2nd list to combine y(n) + k1/2 in the 2nd order RK method
    ! from this list the spawns k2 are created, which can be stored in a 
    ! reinitialized spawnvec arrays
    integer(n_int), allocatable, target :: temp_det_list(:,:)
    integer(n_int), pointer :: temp_det_pointer(:,:)

    ! also use hash table to access those lists 
    type(ll_node), pointer :: temp_det_hash(:)
    integer :: temp_n_hashes
    real(dp) :: temp_hash_frac 

    ! also need the freeslot quantitiy 
    integer, allocatable :: temp_freeslot(:)
    ! also need temporary variables to store the iStart- and iEndFreeSlot 
    ! variables(actually dont need tmp_start var.)
    integer :: temp_iEndFreeslot

end module real_time_data
