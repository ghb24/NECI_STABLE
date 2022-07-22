#include "macros.h"

module tau_search
    use util_mod, only: EnumBase_t
    use constants, only: dp

    better_implicit_none

    ! private
    public :: TauSearchMethod_t, possible_tau_search_methods, &
        tau_search_method, input_tau_search_method


    type, extends(EnumBase_t) :: TauSearchMethod_t
    end type

    type :: PossibleTauSearchMethods_t
        type(TauSearchMethod_t) :: &
            CONVENTIONAL = TauSearchMethod_t(1), &
            HISTOGRAMMING = TauSearchMethod_t(2)
    end type

    type(PossibleTauSearchMethods_t), parameter :: &
        possible_tau_search_methods = PossibleTauSearchMethods_t()

    type(TauSearchMethod_t), allocatable :: &
        tau_search_method, input_tau_search_method

    real(dp) :: min_tau = 0._dp, max_tau = huge(max_tau)

    logical :: scale_tau_to_death

contains


end module
