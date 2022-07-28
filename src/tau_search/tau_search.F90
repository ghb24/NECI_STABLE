#include "macros.h"

module tau_search
    use util_mod, only: EnumBase_t, stop_all
    use constants, only: dp, inum_runs, EPS
    use Parallel_neci, only: MPIAllReduce, MPI_MAX, iProcIndex
    use FciMCData, only: iter, tSinglePartPhase
    use CalcData, only: tau, tPrecond

    better_implicit_none

    ! private
    public :: TauSearchMethod_t, possible_tau_search_methods, &
              tau_search_method, input_tau_search_method, &
              tau_start_val, possible_tau_start, end_of_search_reached, &
              tau_stop_method, possible_tau_stop_methods, &
              scale_tau_to_death, log_death_magnitude



    type, extends(EnumBase_t) :: TauSearchMethod_t
    end type

    type :: PossibleTauSearchMethods_t
        type(TauSearchMethod_t) :: &
            OFF = TauSearchMethod_t(1), &
            CONVENTIONAL = TauSearchMethod_t(2), &
            HISTOGRAMMING = TauSearchMethod_t(3)
    end type

    type(PossibleTauSearchMethods_t), parameter :: &
        possible_tau_search_methods = PossibleTauSearchMethods_t()

    type(TauSearchMethod_t) :: &
        tau_search_method = possible_tau_search_methods%OFF

    type(TauSearchMethod_t), allocatable :: input_tau_search_method

    type, extends(EnumBase_t) :: StopMethod_t
    end type

    type :: PossibleStopMethods_t
        type(StopMethod_t) :: &
            var_shift = StopMethod_t(1), &
            after_iter = StopMethod_t(2), &
            no_change = StopMethod_t(3), &
            off = StopMethod_t(4)
    end type

    type(PossibleStopMethods_t), parameter :: possible_tau_stop_methods = PossibleStopMethods_t()

    type(StopMethod_t) :: tau_stop_method = possible_tau_stop_methods%var_shift

    type, extends(EnumBase_t) :: TauStartVal_t
    end type

    type :: PossibleStartValTau_t
        type(TauStartVal_t) :: &
            user_given = TauStartVal_t(1), &
            tau_factor = TauStartVal_t(2), &
            from_popsfile = TauStartVal_t(3), &
            deterministic = TauStartVal_t(4)
    end type

    type(PossibleStartValTau_t), parameter :: possible_tau_start = PossibleStartValTau_t()

    type(TauStartVal_t), allocatable :: tau_start_val

    real(dp) :: min_tau = 0._dp, max_tau = huge(max_tau)

    logical :: scale_tau_to_death_triggered = .false., t_scale_tau_to_death = .false.
    real(dp) :: max_death_cpt = 0._dp


    interface
        elemental module function end_of_search_reached(curr_tau_search_method) result(res)
            type(TauSearchMethod_t), intent(in) :: curr_tau_search_method
            logical :: res
        end function
    end interface

contains

    elemental function end_of_search_reached(curr_tau_search_method) result(res)
        type(TauSearchMethod_t), intent(in) :: curr_tau_search_method
        logical :: res
        integer :: run
        character(*), parameter :: this_routine = 'end_of_search_reached'

        if (curr_tau_search_method == possible_tau_search_methods%OFF) then
            res = .true.
        else
            if (tau_stop_method == possible_tau_stop_methods%off) then
                res = .false.
            else if (tau_stop_method == possible_tau_stop_methods%var_shift) then
                res = .false.
                do run = 1, inum_runs
                    res = .not. (tSinglePartPhase(run) .or. (tPrecond .and. iter <= 80))
                    if (res) exit
                end do
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
                tau = tau_death

                root_print "******"
                root_print "WARNING: Updating time step due to particle death &
                     &magnitude"
                root_print "This occurs despite variable shift mode"
                root_print "Updating time-step. New time-step = ", tau
                root_print "******"
            end if
        end if

        ! Condition met --> no need to do this again next iteration
        scale_tau_to_death_triggered = .false.
    end subroutine

    subroutine log_death_magnitude(mult)

        ! The same as above, but for particle death

        real(dp), intent(in) :: mult

        if (mult > max_death_cpt) then
            max_death_cpt = mult
            if (t_scale_tau_to_death .and. tau_search_method == possible_tau_search_methods%OFF) then
                scale_tau_to_death_triggered = .true.
            end if
        end if

    end subroutine



end module
