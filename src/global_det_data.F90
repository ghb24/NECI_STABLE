#include "macros.h"

module global_det_data

    use CalcData, only: tContTimeFCIMC, tContTimeFull
    use LoggingData, only: tRDMonFly, tExplicitAllRDM, tTransitionRDMs
    use FciMCData, only: MaxWalkersPart
    use constants
    use util_mod
    implicit none

    ! This is the tag for allocation/deallocation.
    private :: glob_tag
    integer :: glob_tag = 0

    ! This is the data used to find the elements inside the storage array.
    private :: pos_hel, len_hel

    ! The diagonal matrix element is always stored. As it is a real value it
    ! always has a length of 1 (never cplx). Therefore, encode these values
    ! as parameters to assist optimisation.
    integer, parameter :: pos_hel = 1, len_hel = 1

    !The initial population of a determinant at spawning time.
    integer, parameter :: pos_spawn_pop = 2, len_spawn_pop = 1

    !The integral of imaginary time since the spawing of this determinant.
    integer, parameter :: pos_tau_int = 3, len_tau_int = 1
    
    !The integral of shift since the spawning of this determinant.
    integer, parameter :: pos_shift_int = 4, len_shift_int = 1

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

    ! And somewhere to store the actual data.
    real(dp), pointer :: global_determinant_data(:,:) => null()

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

        ! Get the starting positions
        pos_av_sgn = pos_shift_int + len_shift_int
        pos_av_sgn_transition = pos_av_sgn + len_av_sgn
        pos_iter_occ = pos_av_sgn_transition + len_av_sgn_transition
        pos_iter_occ_transition = pos_iter_occ + len_iter_occ
        pos_spawn_rate = pos_iter_occ_transition + len_iter_occ_transition

        tot_len = len_hel + len_spawn_pop + len_tau_int + len_shift_int + len_av_sgn_tot + len_iter_occ_tot + len_spawn_rate

        ! Allocate and log the required memory (globally)
        allocate(global_determinant_data(tot_len, MaxWalkersPart), stat=ierr)
        log_alloc(global_determinant_data, glob_tag, ierr)

        write(6,'(a,f14.6,a)') &
            ' Determinant related persistent storage requires: ', &
            8.0_dp * real(tot_len * MaxWalkersPart,dp) / 1048576_dp, &
            ' Mb / processor'

        ! As an added safety feature
        global_determinant_data = 0.0_dp
                         
    end subroutine


    subroutine clean_global_det_data()

        character(*), parameter :: this_routine = 'clean_global_det_data'

        ! Cleanup the global storage for determinant specific data

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


    subroutine set_spawn_pop(j, t)

        integer, intent(in) :: j
        real(dp), intent(in) :: t

        global_determinant_data(pos_spawn_pop, j) = t

    end subroutine

    function get_spawn_pop(j) result(t)
        
        integer, intent(in) :: j
        real(dp) :: t

        t = global_determinant_data(pos_spawn_pop, j)

    end function

    subroutine reset_tau_int(j)

        integer, intent(in) :: j

        global_determinant_data(pos_tau_int, j) = 0.0_dp

    end subroutine

    subroutine update_tau_int(j, t)

        integer, intent(in) :: j
        real(dp), intent(in) :: t

        global_determinant_data(pos_tau_int, j) = global_determinant_data(pos_tau_int, j) + t

    end subroutine

    function get_tau_int(j) result(t)
        
        integer, intent(in) :: j
        real(dp) :: t

        t = global_determinant_data(pos_tau_int, j)

    end function

    subroutine reset_shift_int(j)

        integer, intent(in) :: j

        global_determinant_data(pos_shift_int, j) = 0.0_dp

    end subroutine

    subroutine update_shift_int(j, t)

        integer, intent(in) :: j
        real(dp), intent(in) :: t

        global_determinant_data(pos_shift_int, j) = global_determinant_data(pos_shift_int, j) + t

    end subroutine

    function get_shift_int(j) result(t)
        
        integer, intent(in) :: j
        real(dp) :: t

        t = global_determinant_data(pos_shift_int, j)

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

end module
