#include "macros.h"

module global_det_data
  use SystemData, only: nel
    use CalcData, only: tContTimeFCIMC, tContTimeFull, tStoredDets, tActivateLAS, &
                        tSeniorInitiators, tAutoAdaptiveShift, tPairedReplicas, tReplicaEstimates
    use LoggingData, only: tRDMonFly, tExplicitAllRDM, tTransitionRDMs
    use FciMCData, only: MaxWalkersPart
    use constants
    use util_mod
    implicit none

    ! This is the tag for allocation/deallocation.
    private :: glob_tag, glob_det_tag
    ! This is the data used to find the elements inside the storage array.
    private :: pos_hel, len_hel
    integer :: glob_tag = 0
    integer :: glob_det_tag = 0

    ! The diagonal matrix element is always stored. As it is a real value it
    ! always has a length of 1 (never cplx). Therefore, encode these values
    ! as parameters to assist optimisation.
    integer :: SeniorsNum, AllSeniorsNum
    integer, parameter :: pos_hel = 1, len_hel = 1

    !The initial population of a determinant at spawning time.
    integer :: pos_spawn_pop, len_spawn_pop

    !The integral of imaginary time since the spawing of this determinant.
    integer :: pos_tau_int, len_tau_int
    
    !The integral of shift since the spawning of this determinant.
    integer :: pos_shift_int, len_shift_int

    ! Total spawns and the number of successful ones according the initiator criterion
    integer :: len_tot_spawns, len_acc_spawns
    integer :: pos_tot_spawns, pos_acc_spawns

    ! Average sign and first occupation of iteration.
    private :: pos_av_sgn, len_av_sgn, pos_iter_occ, len_iter_occ
    integer :: pos_av_sgn, len_av_sgn
    integer :: pos_iter_occ, len_iter_occ

    ! Average sign and first occupation of iteration, for transition RDMs.
    private :: pos_av_sgn_transition, len_av_sgn_transition
    private :: pos_iter_occ_transition, len_iter_occ_transition
    integer :: pos_av_sgn_transition, len_av_sgn_transition
    integer :: pos_iter_occ_transition, len_iter_occ_transition

    ! The total length of average sign and first occupation iteration, for
    ! both the standard and transition RDMs.
    integer :: len_av_sgn_tot, len_iter_occ_tot

    integer :: pos_spawn_rate, len_spawn_rate

    ! Legth of arrays storing estimates to be written to the replica_est file
    integer :: replica_est_len
    ! global storage of history of determinants: number of pos/neg spawns and 
    ! time since a determinant died
    integer :: len_pos_spawns, len_neg_spawns, len_death_timer, len_occ_time
    integer :: pos_pos_spawns, pos_neg_spawns, pos_death_timer, pos_occ_time

    ! lenght of the determinant and its position
    integer :: pos_det_orbs, len_det_orbs

    ! And somewhere to store the actual data.
    real(dp), pointer :: global_determinant_data(:,:) => null()
    integer, pointer :: global_determinants(:,:) => null()

    ! Routines for setting the properties of both standard and transition
    ! RDMs in a combined manner.
    interface set_av_sgn_tot
        module procedure set_av_sgn_tot_sgl
        module procedure set_av_sgn_tot_all
    end interface

    interface set_iter_occ_tot
        module procedure set_iter_occ_tot_sgl
        module procedure set_iter_occ_tot_all
    end interface

    ! Routines for extracting the properties for standard RDMs.
    interface get_av_sgn_standard
        module procedure get_av_sgn_standard_sgl
        module procedure get_av_sgn_standard_all
    end interface

    interface get_iter_occ_standard
        module procedure get_iter_occ_standard_sgl
        module procedure get_iter_occ_standard_all
    end interface

    ! Routines for extracting the properties for transition RDMs.
    interface get_av_sgn_transition
        module procedure get_av_sgn_transition_sgl
        module procedure get_av_sgn_transition_all
    end interface

    interface get_iter_occ_transition
        module procedure get_iter_occ_transition_sgl
        module procedure get_iter_occ_transition_all
    end interface

    ! Routines for extracting the properties of both standard and transition RDMs.
    interface get_av_sgn_tot
        module procedure get_av_sgn_tot_sgl
        module procedure get_av_sgn_tot_all
    end interface

    interface get_iter_occ_tot
        module procedure get_iter_occ_tot_sgl
        module procedure get_iter_occ_tot_all
    end interface

contains

    subroutine init_global_det_data(nrdms_standard, nrdms_transition)

        use FciMCData, only: var_e_num, rep_est_overlap
        use FciMCData, only: var_e_num_all, rep_est_overlap_all
        use FciMCData, only: e_squared_num, e_squared_num_all
        use FciMCData, only: en2_pert, en2_pert_all
        use FciMCData, only: en2_new, en2_new_all
        use FciMCData, only: precond_e_num, precond_denom
        use FciMCData, only: precond_e_num_all, precond_denom_all

        ! Initialise the global storage of determinant specific persistent
        ! data
        !
        ! --> This is the data that should not be transmitted with each
        !     particle
        ! --> It is not stored in the bit representation

        integer, intent(in) :: nrdms_standard, nrdms_transition

        integer :: tot_len
        integer :: ierr
        character(*), parameter :: t_r = 'init_global_det_data'

        ! The position and size of diagonal matrix elements in the array.
        ! This is set as a module wide parameter, rather than as runtime, as
        ! it is constant, and this aids optimisation.
        ! pos_hel = 1
        ! len_hel = 1

        len_spawn_pop = lenof_sign
        if(tSeniorInitiators)then
            len_tau_int = inum_runs
            len_shift_int = inum_runs
        else
            len_tau_int = 0
            len_shift_int = 0
        end if
        if(tAutoAdaptiveShift)then
            len_tot_spawns = inum_runs
            len_acc_spawns = inum_runs
        else
            len_tot_spawns = 0
            len_acc_spawns = 0
        end if
        

        ! If we are using calculating RDMs stochastically, need to include the
        ! average sign and the iteration on which it became occupied.
        if (tRDMonFly .and. .not. tExplicitAllRDM) then
            len_av_sgn = 2*nrdms_standard
            len_iter_occ = 2*nrdms_standard
            ! The total lengths, including both standard and transition RDMs.
            len_av_sgn_tot = 2*nrdms_standard
            len_iter_occ_tot = 2*nrdms_standard
            ! If we are calculating transition RDMs, then we also need to
            ! include sign averages over different sets of blocks,
            ! corresponding to the ground state combined with all other
            ! excited states.
            if (tTransitionRDMs) then
                len_av_sgn_transition = 2*nrdms_transition
                len_iter_occ_transition = 2*nrdms_transition
                len_av_sgn_tot = len_av_sgn + len_av_sgn_transition
                len_iter_occ_tot = len_iter_occ + len_iter_occ_transition
            end if
            write(6, '(" The average current signs before death will be stored&
                       & for use in the RDMs.")')
        else
            len_av_sgn = 0
            len_iter_occ = 0
        end if

        ! If we are using continuous time, and storing the spawning rates
        len_spawn_rate = 0
        if (tContTimeFCIMC .and. tContTimeFull) then
            len_spawn_rate = 1
        end if

        if(tStoredDets) then
           len_det_orbs = nel
        else
           len_det_orbs = 0
        endif

        if(tActivateLAS) then
           len_pos_spawns = lenof_sign
           len_neg_spawns = lenof_sign
        else
           len_pos_spawns = 0
           len_neg_spawns = 0
        endif
        len_death_timer = 1

        ! Get the starting positions
        pos_spawn_pop = pos_hel+len_hel
        pos_tau_int = pos_spawn_pop+len_spawn_pop
        pos_shift_int = pos_tau_int+len_tau_int
        pos_tot_spawns = pos_shift_int+len_shift_int
        pos_acc_spawns = pos_tot_spawns+len_tot_spawns
        pos_av_sgn = pos_acc_spawns + len_acc_spawns
        pos_av_sgn_transition = pos_av_sgn + len_av_sgn
        pos_iter_occ = pos_av_sgn_transition + len_av_sgn_transition
        pos_iter_occ_transition = pos_iter_occ + len_iter_occ
        pos_spawn_rate = pos_iter_occ_transition + len_iter_occ_transition
        pos_pos_spawns = pos_spawn_rate + len_spawn_rate
        pos_neg_spawns = pos_pos_spawns + len_pos_spawns
        pos_death_timer = pos_neg_spawns + len_neg_spawns
        pos_occ_time = pos_death_timer + pos_death_timer

        tot_len = len_hel + len_spawn_pop + len_tau_int + len_shift_int + len_tot_spawns + len_acc_spawns + &
             len_av_sgn_tot + len_iter_occ_tot + len_pos_spawns + len_neg_spawns + &
             len_death_timer + len_occ_time

        if (tPairedReplicas) then
            replica_est_len = lenof_sign/2
        else
            replica_est_len = lenof_sign
        end if

        if (tReplicaEstimates) then
            allocate(var_e_num(replica_est_len), stat=ierr)
            allocate(rep_est_overlap(replica_est_len), stat=ierr)
            allocate(var_e_num_all(replica_est_len), stat=ierr)
            allocate(rep_est_overlap_all(replica_est_len), stat=ierr)
            allocate(e_squared_num(replica_est_len), stat=ierr)
            allocate(e_squared_num_all(replica_est_len), stat=ierr)
            allocate(en2_pert(replica_est_len), stat=ierr)
            allocate(en2_pert_all(replica_est_len), stat=ierr)
            allocate(en2_new(replica_est_len), stat=ierr)
            allocate(en2_new_all(replica_est_len), stat=ierr)
            allocate(precond_e_num(replica_est_len), stat=ierr)
            allocate(precond_denom(replica_est_len), stat=ierr)
            allocate(precond_e_num_all(replica_est_len), stat=ierr)
            allocate(precond_denom_all(replica_est_len), stat=ierr)
        end if

        ! Allocate and log the required memory (globally)
        allocate(global_determinant_data(tot_len, MaxWalkersPart), stat=ierr)
        log_alloc(global_determinant_data, glob_tag, ierr)

        if(tStoredDets) then
           allocate(global_determinants(len_det_orbs, MaxWalkersPart), stat=ierr)
           log_alloc(global_determinants, glob_det_tag, ierr)
        endif

        write(6,'(a,f14.6,a)') &
            ' Determinant related persistent storage requires: ', &
            (4.0_dp*real(len_det_orbs* MaxWalkersPart,dp) &
            + 8.0_dp * real(tot_len * MaxWalkersPart,dp)) / 1048576_dp, &
            ' Mb / processor'

        ! As an added safety feature
        global_determinant_data = 0.0_dp
        if(tStoredDets) global_determinants = 0
                         
    end subroutine


    subroutine clean_global_det_data()

        character(*), parameter :: this_routine = 'clean_global_det_data'

        ! Cleanup the global storage for determinant specific data

        if(associated(global_determinants)) then
           deallocate(global_determinants)
           log_dealloc(glob_det_tag)
           glob_det_tag = 0
           nullify(global_determinants)
        endif

        if (associated(global_determinant_data)) then
            deallocate(global_determinant_data)
            log_dealloc(glob_tag)
            glob_tag = 0
            nullify(global_determinant_data)
        end if

    end subroutine


    ! These are accessor functions to currentH --> they are data access/setting
    !
    ! Access the global_determinant_data structure by site index (j).

    subroutine set_det_diagH(j, hel_r)

        ! Diagonal matrix elements are real --> store a real value

        integer, intent(in) :: j
        real(dp), intent(in) :: hel_r

        global_determinant_data(pos_hel, j) = hel_r

    end subroutine

    function det_diagH(j) result(hel_r)
        
        integer, intent(in) :: j
        real(dp) :: hel_r

        hel_r = global_determinant_data(pos_hel, j)

    end function


    subroutine set_spawn_pop(j, part, t)

        integer, intent(in) :: j, part
        real(dp), intent(in) :: t

        global_determinant_data(pos_spawn_pop+part-1, j) = t

    end subroutine

    function get_spawn_pop(j, part) result(t)
        
        integer, intent(in) :: j, part
        real(dp) :: t

        t = global_determinant_data(pos_spawn_pop+part-1, j)

    end function

    subroutine set_all_spawn_pops(j, t)

        integer, intent(in) :: j
        real(dp), dimension(lenof_sign), intent(in) :: t

        global_determinant_data(pos_spawn_pop:pos_spawn_pop+len_spawn_pop-1, j) = t(:)

    end subroutine

    function get_all_spawn_pops(j) result(t)
        
        integer, intent(in) :: j
        real(dp), dimension(lenof_sign) :: t

        t(:) = global_determinant_data(pos_spawn_pop:pos_spawn_pop+len_spawn_pop-1, j)

    end function


    subroutine reset_all_tau_ints(j)

        integer, intent(in) :: j

        global_determinant_data(pos_tau_int:pos_tau_int+len_tau_int-1, j) = 0.0_dp

    end subroutine

    subroutine reset_tau_int(j, run)

        integer, intent(in) :: j, run

        global_determinant_data(pos_tau_int+run-1, j) = 0.0_dp

    end subroutine

    subroutine update_tau_int(j, run, t)

        integer, intent(in) :: j, run
        real(dp), intent(in) :: t

        global_determinant_data(pos_tau_int+run-1, j) = global_determinant_data(pos_tau_int+run-1, j) + t

    end subroutine

    function get_tau_int(j, run) result(t)
        
        integer, intent(in) :: j,run
        real(dp) :: t

        t = global_determinant_data(pos_tau_int+run-1, j)

    end function

    subroutine reset_all_shift_ints(j)

        integer, intent(in) :: j

        global_determinant_data(pos_shift_int:pos_shift_int+len_shift_int-1, j) = 0.0_dp

    end subroutine

    subroutine reset_shift_int(j, run)

        integer, intent(in) :: j, run

        global_determinant_data(pos_shift_int+run-1, j) = 0.0_dp

    end subroutine

    subroutine update_shift_int(j, run, t)

        integer, intent(in) :: j, run
        real(dp), intent(in) :: t

        global_determinant_data(pos_shift_int+run-1, j) = global_determinant_data(pos_shift_int+run-1, j) + t

    end subroutine

    function get_shift_int(j, run) result(t)
        
        integer, intent(in) :: j, run
        real(dp) :: t

        t = global_determinant_data(pos_shift_int+run-1, j)

    end function

    subroutine reset_all_tot_spawns(j)

        integer, intent(in) :: j

        global_determinant_data(pos_tot_spawns:pos_tot_spawns+len_tot_spawns-1, j) = 0.0_dp

    end subroutine

    subroutine reset_tot_spawns(j, run)

        integer, intent(in) :: j, run

        global_determinant_data(pos_tot_spawns+run-1, j) = 0.0_dp

    end subroutine

    subroutine update_tot_spawns(j, run, t)

        integer, intent(in) :: j, run
        real(dp), intent(in) :: t

        global_determinant_data(pos_tot_spawns+run-1, j) = global_determinant_data(pos_tot_spawns+run-1, j) + t

    end subroutine

    function get_tot_spawns(j, run) result(t)
        
        integer, intent(in) :: j, run
        real(dp) :: t

        t = global_determinant_data(pos_tot_spawns+run-1, j)

    end function

    subroutine set_tot_acc_spawns(fvals, ndets, initial) 
      implicit none
      integer, intent(in) :: ndets
      real(dp), intent(in) :: fvals(2*inum_runs, ndets)
      integer, intent(in), optional :: initial

      integer :: j, run, start

      if(present(initial)) then
         start = initial
      else
         start = 1
      endif

      ! set all values of tot/acc spawns using the read-in values from fvals
      ! this is used in popsfile read-in to get the values from the previous calculation
      do j = 1, ndets
         do run=1, inum_runs
            global_determinant_data(pos_acc_spawns+run-1,j + start - 1) = &
                 fvals(run,j)
            global_determinant_data(pos_tot_spawns+run-1,j + start - 1) = &
                 fvals(inum_runs+run,j)
         end do
      end do
    end subroutine set_tot_acc_spawns

#ifdef __USE_HDF
    ! nasty bit of code to cope with hdf5 I/O which is using integer(hsize_t)
    subroutine set_tot_acc_spawn_hdf5Int(fvals, j)
      use hdf5
      implicit none
      integer(hsize_t), intent(in) :: fvals(:)
      integer, intent(in) :: j

      integer :: run
      real(dp) :: realVal = 0.0_dp

      do run = 1, inum_runs
         global_determinant_data(pos_acc_spawns+run-1,j) = transfer(fvals(run),realVal)
         global_determinant_data(pos_tot_spawns+run-1,j) = transfer(fvals(run+inum_runs),realVal)
      end do
    end subroutine set_tot_acc_spawn_hdf5Int
#endif

    subroutine writeFFuncAsInt(ndets, fvals)
      implicit none
      integer, intent(in) :: ndets
      integer(n_int), intent(inout) :: fvals(:,:)

      integer :: j, k

      ! write the acc. and tot. spawns per determinant in a contiguous array
      ! fvals(:,j) = (acc, tot) for determinant j (2*inum_runs in size)
      do j = 1, nDets
         do k = 1, inum_runs
            fvals(k,j) = transfer(get_acc_spawns(j,k), fvals(k,j))
         end do
         do k = 1, inum_runs
            fvals(k+inum_runs,j) = transfer(get_tot_spawns(j,k), fvals(k,j))
         end do
      end do
    end subroutine writeFFuncAsInt

    subroutine writeFFunc(ndets, fvals)
      implicit none
      integer, intent(in) :: ndets
      real(dp), intent(inout) :: fvals(:,:)

      integer :: j, k

      ! write the acc. and tot. spawns per determinant in a contiguous array
      ! fvals(:,j) = (acc, tot) for determinant j (2*inum_runs in size)
      do j = 1, nDets
         do k = 1, inum_runs
            fvals(k,j) = get_acc_spawns(j,k)
         end do
         do k = 1, inum_runs
            fvals(k+inum_runs,j) = get_tot_spawns(j,k)
         end do
      end do
    end subroutine writeFFunc

  !------------------------------------------------------------------------------------------!

    subroutine reset_all_acc_spawns(j)

        integer, intent(in) :: j

        global_determinant_data(pos_acc_spawns:pos_acc_spawns+len_acc_spawns-1, j) = 0.0_dp

    end subroutine

    subroutine reset_acc_spawns(j, run)

        integer, intent(in) :: j, run

        global_determinant_data(pos_acc_spawns+run-1, j) = 0.0_dp

    end subroutine

    subroutine update_acc_spawns(j, run, t)

        integer, intent(in) :: j, run
        real(dp), intent(in) :: t

        global_determinant_data(pos_acc_spawns+run-1, j) = global_determinant_data(pos_acc_spawns+run-1, j) + t

    end subroutine

    function get_acc_spawns(j, run) result(t)
        
        integer, intent(in) :: j, run
        real(dp) :: t

        t = global_determinant_data(pos_acc_spawns+run-1, j)

    end function

    subroutine set_av_sgn_tot_sgl (j, part, av_sgn)

        integer, intent(in) :: j, part
        real(dp), intent(in) :: av_sgn

        global_determinant_data(pos_av_sgn + part - 1, j) = av_sgn

    end subroutine

    subroutine set_av_sgn_tot_all (j, av_sgn)

        integer, intent(in) :: j
        real(dp), intent(in) :: av_sgn(len_av_sgn_tot)

        global_determinant_data(pos_av_sgn:pos_av_sgn + len_av_sgn_tot - 1, j) = &
            av_sgn

    end subroutine

    subroutine set_iter_occ_tot_sgl (j, part, iter_occ)

        ! It is unusual, but all the RDM code uses real(dp) values for the
        ! current iteration. As floats store integers exactly up to a
        ! sensible limit, this is just fine, and simplifies this code.
        !
        ! But it is weird.

        integer, intent(in) :: j, part
        real(dp), intent(in) :: iter_occ

        global_determinant_data(pos_iter_occ + part - 1, j) = iter_occ

    end subroutine

    subroutine set_iter_occ_tot_all (j, iter_occ)

        ! It is unusual, but all the RDM code uses real(dp) values for the
        ! current iteration. As floats store integers exactly up to a
        ! sensible limit, this is just fine, and simplifies this code.
        !
        ! But it is weird.

        integer, intent(in) :: j
        real(dp), intent(in) :: iter_occ(len_iter_occ_tot)

        global_determinant_data(pos_iter_occ: &
                                pos_iter_occ + len_iter_occ_tot - 1, j) = iter_occ

    end subroutine

    ! -----Routines for extracting the properties of standard RDMs------------

    function get_av_sgn_standard_sgl (j, part) result(av_sgn)

        integer, intent(in) :: j, part
        real(dp) :: av_sgn

        av_sgn = global_determinant_data(pos_av_sgn + part - 1, j)

    end function

    function get_av_sgn_standard_all (j) result(av_sgn)

        integer, intent(in) :: j
        real(dp) :: av_sgn(len_av_sgn)

        av_sgn = global_determinant_data(pos_av_sgn:&
                                       pos_av_sgn + len_av_sgn - 1, j)

    end function

    function get_iter_occ_standard_sgl (j, part) result(iter_occ)

        integer, intent(in) :: j, part
        real(dp) :: iter_occ

        iter_occ = global_determinant_data(pos_iter_occ + part - 1, j)

    end function

    function get_iter_occ_standard_all (j) result(iter_occ)

        integer, intent(in) :: j
        real(dp) :: iter_occ(len_iter_occ)

        iter_occ = global_determinant_data(pos_iter_occ: &
                                   pos_iter_occ + len_iter_occ - 1, j)

    end function

    ! -----Routines for extracting the properties of transition RDMs-----------

    function get_av_sgn_transition_sgl (j, part) result(av_sgn)

        integer, intent(in) :: j, part
        real(dp) :: av_sgn

        av_sgn = global_determinant_data(pos_av_sgn_transition + part - 1, j)

    end function

    function get_av_sgn_transition_all (j) result(av_sgn)

        integer, intent(in) :: j
        real(dp) :: av_sgn(len_av_sgn_transition)

        av_sgn = global_determinant_data(pos_av_sgn_transition:&
                                       pos_av_sgn_transition + len_av_sgn_transition - 1, j)

    end function

    function get_iter_occ_transition_sgl (j, part) result(iter_occ)

        integer, intent(in) :: j, part
        real(dp) :: iter_occ

        iter_occ = global_determinant_data(pos_iter_occ_transition + part - 1, j)

    end function

    function get_iter_occ_transition_all (j) result(iter_occ)

        integer, intent(in) :: j
        real(dp) :: iter_occ(len_iter_occ)

        iter_occ = global_determinant_data(pos_iter_occ_transition: &
                                   pos_iter_occ_transition + len_iter_occ_transition - 1, j)

    end function

    ! -----Routines for extracting the properties of both standard and--------
    ! -----transition RDMs together-------------------------------------------

    function get_av_sgn_tot_sgl (j, part) result(av_sgn)

        integer, intent(in) :: j, part
        real(dp) :: av_sgn

        av_sgn = global_determinant_data(pos_av_sgn + part - 1, j)

    end function

    function get_av_sgn_tot_all (j) result(av_sgn)

        integer, intent(in) :: j
        real(dp) :: av_sgn(len_av_sgn_tot)

        av_sgn = global_determinant_data(pos_av_sgn:&
                                       pos_av_sgn + len_av_sgn_tot - 1, j)

    end function

    function get_iter_occ_tot_sgl (j, part) result(iter_occ)

        integer, intent(in) :: j, part
        real(dp) :: iter_occ

        iter_occ = global_determinant_data(pos_iter_occ + part - 1, j)

    end function

    function get_iter_occ_tot_all (j) result(iter_occ)

        integer, intent(in) :: j
        real(dp) :: iter_occ(len_iter_occ_tot)

        iter_occ = global_determinant_data(pos_iter_occ: &
                                   pos_iter_occ + len_iter_occ_tot - 1, j)

    end function

    function get_spawn_rate(j) result(rate)

        integer, intent(in) :: j
        real(dp) :: rate
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'get_spawn_rate'
#endif

        ASSERT(tContTimeFCIMC .and. tContTimeFull)
        rate = global_determinant_data(pos_spawn_rate, j)

    end function
    
    subroutine set_spawn_rate(j, rate) 

        integer, intent(in) :: j
        real(dp), intent(in) :: rate
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'set_spawn_rate'
#endif

        ASSERT(tContTimeFCIMC .and. tContTimeFull)
        global_determinant_data(pos_spawn_rate, j) = rate

    end subroutine

  !------------------------------------------------------------------------------------------!
  ! functions for storing information on the history of spawns onto a determinant
  !------------------------------------------------------------------------------------------!

    subroutine store_spawn(j, spawn_sgn)
      implicit none
      integer, intent(in) :: j
      real(dp), intent(in) :: spawn_sgn(lenof_sign)

      integer :: part

      do part = 1, lenof_sign
         if(spawn_sgn(part) > eps) then
            global_determinant_data(pos_pos_spawns+part-1,j) = &
                 global_determinant_data(pos_pos_spawns+part-1,j) + abs(spawn_sgn(part))
         else if(spawn_sgn(part) < -eps) then
            global_determinant_data(pos_neg_spawns+part-1,j) = &
                 global_determinant_data(pos_neg_spawns+part-1,j) + abs(spawn_sgn(part))
         endif
      end do
    end subroutine store_spawn

  !------------------------------------------------------------------------------------------!

    pure function get_pos_spawns(j) result(avSpawn)
      implicit none
      integer, intent(in) :: j
      real(dp) :: avSpawn(lenof_sign)

      avSpawn = global_determinant_data(pos_pos_spawns:(pos_pos_spawns + lenof_sign - 1),j)
    end function get_pos_spawns

  !------------------------------------------------------------------------------------------!

    pure function get_neg_spawns(j) result(avSpawn)
      implicit none
      integer, intent(in) :: j
      real(dp) :: avSpawn(lenof_sign)

      avSpawn = global_determinant_data(pos_neg_spawns:(pos_neg_spawns + lenof_sign - 1),j)
    end function get_neg_spawns

  !------------------------------------------------------------------------------------------!

    subroutine clock_occ_time(j)
      implicit none
      integer, intent(in) :: j
      
      global_determinant_data(pos_occ_time,j) = global_determinant_data(pos_occ_time,j) + 1
      global_determinant_data(pos_pos_spawns:(pos_death_timer-1),j) = &
           global_determinant_data(pos_pos_spawns:(pos_death_timer-1),j) * &
           (global_determinant_data(pos_occ_time,j)-1)/global_determinant_data(pos_occ_time,j)

    end subroutine clock_occ_time

  !------------------------------------------------------------------------------------------!

    subroutine reset_occ_time(j)
      implicit none
      integer, intent(in) :: j
      
      global_determinant_data(pos_occ_time,j) = 0
    end subroutine reset_occ_time

  !------------------------------------------------------------------------------------------!

    subroutine clock_death_timer(j)
      implicit none
      integer, intent(in) :: j
      
      global_determinant_data(pos_death_timer,j) = global_determinant_data(pos_death_timer,j) + 1.0_dp
    end subroutine clock_death_timer

  !------------------------------------------------------------------------------------------!

    function get_death_timer(j) result(niter)
      implicit none
      integer, intent(in) :: j
      real(dp) :: niter

      niter = global_determinant_data(pos_death_timer,j)

    end function get_death_timer

  !------------------------------------------------------------------------------------------!

    subroutine mark_death(j) 
      implicit none
      integer, intent(in) :: j

      global_determinant_data(pos_death_timer,j) = -1.0_dp
    end subroutine mark_death
      
      
  !------------------------------------------------------------------------------------------!

    subroutine reset_death_timer(j)
      implicit none
      integer, intent(in) :: j
      
      global_determinant_data(pos_death_timer,j) = 0.0_dp
    end subroutine reset_death_timer

  !------------------------------------------------------------------------------------------!
  !    Global storage for storing nI for each occupied determinant to save time for
  !    conversion from ilut to nI
  !------------------------------------------------------------------------------------------!

    subroutine store_decoding(j, nI)
      implicit none
      integer, intent(in) :: j, nI(nel)
      
      if(tStoredDets) then
         global_determinants(:,j) = nI
      endif
    end subroutine store_decoding

    function get_determinant(j) result(nI)
      implicit none
      integer, intent(in) :: j
      integer :: nI(nel)
      
      if(tStoredDets) then
         nI = global_determinants(:,j)
      else
         nI = 0
      endif
    end function get_determinant

end module
