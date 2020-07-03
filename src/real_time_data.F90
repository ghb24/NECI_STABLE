#include "macros.h"

! data module of the real-time implementation of FCIQMC

module real_time_data

    use constants
    use MemoryManager, only: TagIntType
    use FciMCData, only: perturbation, ll_node, fcimc_iter_data
    implicit none

    ! number of copies to start the real-time simulation from
    integer :: n_real_time_copies, cnt_real_time_copies
    logical :: t_prepare_real_time

    ! for the actual calculation:
    ! global flag indicating real-time calculation
    ! rotated_time_setup: flag to indicate whether pure real time is
    ! used or not
    logical :: t_real_time_fciqmc, t_new_stats_file, t_rotated_time, tStaticShift, &
         tDynamicCoreSpace, tRealTimePopsfile, tStabilizerShift, tLimitShift, &
         tDynamicAlpha, tDynamicDamping, tInfInit, tStartVariation, tOverpopulate, &
         tNewOverlap, tOnlyPositiveShift, tHFOverlap

    logical :: tLowerThreshold
    ! also use a second iter_data type to keep track of the 2 distinct
    ! spawning events
    type(fcimc_iter_data) :: second_spawn_iter_data

    ! store the type of greensfunction as an integer, or maybe later on even
    ! create a type like nick did, to store all the general information like
    ! timestep, which orbitals are manipulated etc. in it
    integer :: gf_type      ! -1 indicates a the lesser, 1 the greater GF
    integer :: normsize ! number of combinations avaliable for calculating overlaps
    integer :: gf_count ! number of different Green's functions to be evaluated
    integer :: allGfs ! 0 -> only one GF, 1 -> all lesser GFs, 2 -> all greater GFs
    integer :: stepsAlpha ! number of cycles after which the rotation angle is updated
    integer :: rotThresh ! number of walkers on which variation of alpha starts
    ! prefactors for rescaling of the wavefunction, if enabled
    real(dp) :: alphaDamping, etaDamping ! prefactors for adjustment of rotation angle and damping
    real(dp) :: stabilizerThresh

    ! the angle used for rotation of time into complex plane
    ! it is better readable to use new variables for real and imaginary part
    ! of the rotated time
    real(dp) :: tau_imag, tau_real
    ! logging the total elapsed time on both axes
    real(dp) :: elapsedRealTime, elapsedImagTime
    ! value of static shift if enabled
    real(dp) :: asymptoticShift
    ! maximum shift that can be used and counter to track the cycles
    real(dp) :: shiftLimit
    integer, allocatable :: numCycShiftExcess(:)
    real(dp), allocatable :: popSnapshot(:), allPopSnapShot(:)
    integer, allocatable :: snapshotOrbs(:)
    integer :: numSnapshotOrbs
    ! averaged overlaps
    real(dp), allocatable :: overlap_real(:), overlap_imag(:)
    ! and dynamic reduced norm, used only when rescaling the wave function
    ! i.e. outdated feature
    complex(dp), allocatable :: dyn_norm_red(:,:)
    ! norms are complex because they are taken between different replicas
    ! -> not necesserily entirely real
    ! need a global variable for the overlap <y(0)|y(t)>
    ! determined by the max. cycle
!     real(dp), allocatable :: gf_overlap(:)
    ! for tests now only make it of length 1
    complex(dp), allocatable :: gf_overlap(:,:)
    ! also store the norm of the time-evolved wavefunction
    complex(dp), allocatable :: dyn_norm_psi(:)
    ! the ground state energy of the N-particle problem which is added as a global shift
    real(dp), allocatable :: gs_energy(:)
    ! the damping due to the DiagSft contribution in the imaginary time propagation
    real(dp), allocatable :: shift_damping(:)
    real(dp), allocatable :: TotPartsPeak(:)
    ! Buffer to store values of alpha when dynamically updating it to get some
    ! information on trends in alpha
    real(dp), allocatable :: alphaLog(:)
    integer :: alphaLogSize, alphaLogPos

    ! also store the current overlap of the cycle..

    complex(dp), allocatable :: current_overlap(:,:)

    ! flag to indicate usage of verlet scheme
    logical :: tVerletScheme, tVerletSweep

    ! buffers for the verlet scheme
    integer(n_int), allocatable :: spawnBuf(:,:)
    integer :: spawnBufSize

    ! initialization variables for the verlet scheme
    integer :: iterInit

    ! cache for delta psi
    integer(n_int), pointer :: dpsi_cache(:,:)
    integer :: dpsi_size, max_cache_size, backup_size

    ! also store the norm of the perturbed ground state to adjust the overlap
    complex(dp), allocatable :: pert_norm(:,:)

    ! need additional info of the original walker number and the number of
    ! walkers remaining in the perturbed ground state
    integer(int64) :: TotWalkers_orig

    type perturbed_state
       ! type for a state that is used as a reference state for the overlap
       integer :: nDets
       integer(n_int), allocatable :: dets(:,:)
    end type perturbed_state

    type(perturbed_state), allocatable :: overlap_states(:)

    type real_time_type
        ! and a end time to stop the simulation afterwards
        real(dp) :: max_time = -1.0_dp

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
        real(dp) :: time_angle = 0.0_dp

    end type real_time_type

    type(real_time_type) :: real_time_info

    ! for hubbard model: Info if the perturbation operators are in k-space
    logical :: t_kspace_operators
    ! and the corresponding phases for the fourier transformation
    real(dp), allocatable :: phase_factors(:)

    ! need a 2nd list to combine y(n) + k1/2 in the 2nd order RK method
    ! from this list the spawns k2 are created, which can be stored in a
    ! reinitialized spawnvec arrays
    integer(n_int), allocatable, target :: temp_det_list(:,:)
    integer(n_int), pointer :: temp_det_pointer(:,:)

    ! also use hash table to access those lists
    type(ll_node), pointer :: temp_det_hash(:)
    integer :: temp_n_hashes
    real(dp) :: temp_hash_frac

    ! also need to store the original number of determinants(and walkers maybe)
    ! of the y(n) list to reload correctly
    integer(int64) :: temp_totWalkers
    integer :: MaxSpawnedDiag

    ! also start to store the diagonal "spawns" in the second rt-fciqmc loop
    ! in a seperate Spawned Parts array, to better keep track of stats and
    ! do seperate different distinct processes
    integer(n_int), pointer :: DiagParts(:,:)

    ! also need the associated actual arrays to point to
    integer(n_int), allocatable, target :: DiagVec(:,:)

    ! i dont think i need a hash table to go with that..
    ! but i need this valid_spawned list thingy..
    ! which holds the next free slot to spawn to.. for each proc
    integer, allocatable :: temp_freeslot(:)
    integer :: temp_iendfreeslot, valid_diag_spawns

    ! also keep a global var. of the number of diag-spawns in a cycle
    ! and limit the number of spawns per determinant
    integer :: n_diag_spawned, nspawnMax

    ! keep stats track of the two runge kutta steps seperately -> need new
    ! global variables!
    real(dp), allocatable :: TotPartsStorage(:), TotPartsLastAlpha(:)

    integer(TagIntType) :: DiagVecTag = 0

    ! use a global integer to specifiy the current runge-kutta step (1 or 2)
    ! to keep stats correclty
    integer :: runge_kutta_step

    ! for saving/loading the trajectories
    logical :: tLogTrajectory, tReadTrajectory, tLiveTrajectory
    integer :: iunitCycLog
    real(dp), allocatable :: tauCache(:), alphaCache(:)
    character(255) :: trajFile
    ! These are required in the extended semi-stochastic treatment
    logical :: tGenerateCoreSpace, tGZero

    ! For corespace construction
    real(dp) :: wn_threshold
    integer(n_int), pointer :: core_space_buf(:,:)
    integer :: csbuf_size, corespace_log_interval
    type(ll_node), pointer :: ssht(:)

    ! energy-dependent damping (quadratic)
    logical :: t_quad_damp = .false.

end module real_time_data
