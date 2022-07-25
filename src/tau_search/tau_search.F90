#include "macros.h"

module tau_search
    use util_mod, only: EnumBase_t, stop_all
    use constants, only: dp, inum_runs
    use FciMCData, only: iter, tSinglePartPhase
    use CalcData, only: tPrecond

    better_implicit_none

    ! private
    public :: TauSearchMethod_t, possible_tau_search_methods, &
              tau_search_method, input_tau_search_method, &
              tau_start_val, possible_tau_start, end_of_search_reached, &
              tau_stop_method, possible_tau_stop_methods


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
        tau_search_method = possible_tau_search_methods%OFF, &
        input_tau_search_method = possible_tau_search_methods%OFF

    type, extends(EnumBase_t) :: StopMethod_t
    end type

    type :: PossibleStopMethods_t
        type(StopMethod_t) :: &
            var_shift = StopMethod_t(1), &
            after_iter = StopMethod_t(2), &
            no_change = StopMethod_t(3)
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

    logical :: tSearchTauDeath = .false., scale_tau_to_death = .false.

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
            if (tau_stop_method == possible_tau_stop_methods%var_shift) then
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

end module
