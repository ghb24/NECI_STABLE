#include "macros.h"

module tau_main
    use util_mod, only: EnumBase_t, stop_all, near_zero
    use constants, only: dp, inum_runs, EPS, stdout
    use Parallel_neci, only: MPIAllReduce, MPI_MAX, iProcIndex, root
    use FciMCData, only: iter, tSinglePartPhase, VaryShiftIter
    use CalcData, only: tPrecond
    use util_mod, only: clamp

    better_implicit_none

    private
    public :: TauSearchMethod_t, possible_tau_search_methods, &
              tau_search_method, input_tau_search_method, &
              tau_start_val, possible_tau_start, end_of_search_reached, &
              tau_stop_method, possible_tau_stop_methods, &
              scale_tau_to_death, log_death_magnitude, &
              tau, taufactor, min_tau, max_tau, &
              scale_tau_to_death_triggered, t_scale_tau_to_death, &
              max_death_cpt, assign_value_to_tau, &
              stop_options, find_tau_from_refdet_conn, &
              readpops_but_tau_not_from_popsfile, &
              max_permitted_spawn, &
              init_tau, finalize_tau, stop_tau_search

    protected :: tau


    real(dp) :: tau = 0._dp
        !! The time-step itself

    type, extends(EnumBase_t) :: TauSearchMethod_t
        character(20) :: str
    end type

    type :: PossibleTauSearchMethods_t
        type(TauSearchMethod_t) :: &
            OFF = TauSearchMethod_t(1, 'Off'), &
            CONVENTIONAL = TauSearchMethod_t(2, 'Conventional'), &
            HISTOGRAMMING = TauSearchMethod_t(3, 'Histogramming')
    end type

    type(PossibleTauSearchMethods_t), parameter :: &
        possible_tau_search_methods = PossibleTauSearchMethods_t()

    type(TauSearchMethod_t) :: &
        tau_search_method = possible_tau_search_methods%OFF

    type(TauSearchMethod_t), allocatable :: input_tau_search_method

    type, extends(EnumBase_t) :: StopMethod_t
        character(45) :: str
    end type

    type :: PossibleStopMethods_t
        type(StopMethod_t) :: &
            var_shift = StopMethod_t(1, 'Variable Shift reached'), &
            max_iter = StopMethod_t(2, 'n-th iteration reached'), &
            max_eq_iter = StopMethod_t(3, 'n-th iteration after variable shift reached'), &
            no_change = StopMethod_t(4, 'n iterations without change of tau'), &
            n_opts = StopMethod_t(5, 'n optimizations of tau'), &
            changevars = StopMethod_t(6, 'Manual change via `CHANGEVARS` file'), &
            off = StopMethod_t(7, 'Off')
    end type

    type(PossibleStopMethods_t), parameter :: possible_tau_stop_methods = PossibleStopMethods_t()

    type(StopMethod_t) :: tau_stop_method = possible_tau_stop_methods%var_shift

    type :: TauSearchData_t
        integer :: last_change_of_tau = 0
            !! At which iteration was tau changed last?
        integer :: n_opts = 0
            !! How often was tau changed?
    end type

    type(TauSearchData_t) :: search_data

    type :: StopOptions_t
        integer :: max_iter = huge(0)
            !! Number of iterations, after which we stop searching.
        integer :: max_eq_iter = huge(0)
            !! Number of iterations **after** reaching variable shift mode,
            !!      after which we stop searching.
        integer :: max_iter_without_change = huge(0)
            !! Number of iterations without a change of tau,
            !!      after which we stop searching.
        integer :: max_n_opts = huge(0)
            !! Number of optimizations of tau, after which we stop searching
    end type

    type(StopOptions_t) :: stop_options

    type, extends(EnumBase_t) :: TauStartVal_t
        character(40) :: str
    end type

    type :: PossibleStartValTau_t
        type(TauStartVal_t) :: &
            user_given = TauStartVal_t(1, 'User defined'), &
            tau_factor = TauStartVal_t(2, 'Tau factor'), &
            from_popsfile = TauStartVal_t(3, 'Popsfile'), &
            refdet_connections = TauStartVal_t(4, 'Reference determinant connections'), &
            deterministic = TauStartVal_t(5, 'Deterministic'), &
            not_needed = TauStartVal_t(6, 'Not required')
    end type

    type(PossibleStartValTau_t), parameter :: possible_tau_start = PossibleStartValTau_t()

    type(TauStartVal_t), allocatable :: tau_start_val

    real(dp) :: min_tau = 0._dp, max_tau = huge(max_tau), taufactor = 0._dp

    logical :: scale_tau_to_death_triggered = .false., t_scale_tau_to_death = .false.
    real(dp) :: max_death_cpt = 0._dp, max_permitted_spawn


    logical :: readpops_but_tau_not_from_popsfile = .false.

    interface
        ! This is implemented in a submodule
        module subroutine find_tau_from_refdet_conn()
        end subroutine

        module subroutine init_tau()
        end subroutine

        module subroutine stop_tau_search(stop_method)
            type(StopMethod_t), intent(in) :: stop_method
        end subroutine

        module subroutine finalize_tau()
        end subroutine
    end interface

contains

    elemental function end_of_search_reached(curr_tau_search_method, stop_method) result(res)
        type(TauSearchMethod_t), intent(in) :: curr_tau_search_method
        type(StopMethod_t), intent(in) :: stop_method
        logical :: res
        integer :: run
        character(*), parameter :: this_routine = 'end_of_search_reached'

        if (curr_tau_search_method == possible_tau_search_methods%OFF) then
            res = .true.
        else
            if (stop_method == possible_tau_stop_methods%off) then
                res = .false.
            else if (stop_method == possible_tau_stop_methods%var_shift) then
                res = .false.
                do run = 1, inum_runs
                    res = .not. (tSinglePartPhase(run) .or. (tPrecond .and. iter <= 80))
                    if (res) exit
                end do
            else if (stop_method == possible_tau_stop_methods%max_iter) then
                res = iter >= stop_options%max_iter
            else if (stop_method == possible_tau_stop_methods%max_eq_iter) then
                res = any(.not. tSinglePartPhase) .and. (iter - maxval(VaryShiftIter)) >= stop_options%max_eq_iter
            else if (stop_method == possible_tau_stop_methods%no_change) then
                res = (iter - search_data%last_change_of_tau) >= stop_options%max_iter_without_change
            else if (stop_method == possible_tau_stop_methods%n_opts) then
                res = search_data%n_opts >= stop_options%max_n_opts
            else
                call stop_all(this_routine, "has to be implemented")
            end if
        end if
    end function

    subroutine scale_tau_to_death()
        real(dp) :: mpi_tmp, tau_death
        debug_function_name("scale_tau_to_death")
        ! Check that the override has actually occurred.
        ASSERT(scale_tau_to_death_triggered)
        ASSERT(tau_search_method == possible_tau_search_methods%OFF)

        ! The range of tau is restricted by particle death. It MUST be <=
        ! the value obtained to restrict the maximum death-factor to 1.0.
        call MPIAllReduce(max_death_cpt, MPI_MAX, mpi_tmp)
        max_death_cpt = mpi_tmp
        ! again, this only makes sense if there has been some death
        if (max_death_cpt > EPS) then
            tau_death = 1.0_dp / max_death_cpt

            ! If this actually constrains tau, then adjust it!
            if (tau_death < tau) then
                call assign_value_to_tau(tau_death, 'Update due to particle death magnitude.')
            end if
        end if

        ! Condition met --> no need to do this again next iteration
        scale_tau_to_death_triggered = .false.
    end subroutine

    subroutine log_death_magnitude(mult)
        real(dp), intent(in) :: mult

        if (mult > max_death_cpt) then
            max_death_cpt = mult
            if (t_scale_tau_to_death .and. tau_search_method == possible_tau_search_methods%OFF) then
                scale_tau_to_death_triggered = .true.
            end if
        end if
    end subroutine

    subroutine assign_value_to_tau(new_tau, reason)
        !! Assign `new_tau` to `tau`
        !!
        !! `new_tau` has to be `min_tau <= new_tau <= max_tau`.
        !! If the change of `tau` is sufficiently large (determined by `threshhold`),
        !! then `reason` is printed and the change is registered.
        !! This is relevant for the stop-methods that depend on the last relevant
        !! change of tau.
        real(dp), intent(in) :: new_tau
        character(len=*), intent(in) :: reason
            !! Message that gets printed when change was sufficiently large.
        character(*), parameter :: this_routine = 'assign_value_to_tau'

        if (.not. (min_tau <= new_tau .and. new_tau <= max_tau)) then
            call stop_all(&
                    this_routine, &
                    '.not. (min_tau <= new_tau .and. new_tau <= max_tau)')
        end if


        if (large_change(tau, new_tau)) then
            if (iProcIndex == root) then
                write(stdout, '(A, E13.6, 1x, A, E13.6)') '>>> Changing tau:', tau, '->', new_tau
                write(stdout, '(A, A)') '>>> Reason: ', reason
            end if
            search_data%last_change_of_tau = iter
            search_data%n_opts = search_data%n_opts + 1
        end if
        tau = new_tau
    end subroutine

    elemental function large_change(old_tau, new_tau) result(res)
        !! If the change of `old_tau` to `new_tau` is considered large.
        real(dp), intent(in) :: old_tau, new_tau
        logical :: res
        real(kind(new_tau)), parameter :: threshhold = 0.001_dp
            !! Threshhold for the relative change of tau.
        if (near_zero(old_tau)) then
            ! This is at initialization
            res = .false.
        else if (abs(old_tau - new_tau) / old_tau > threshhold) then
            res = .true.
        else
            res = .false.
        end if
    end function
end module
