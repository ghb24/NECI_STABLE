#include "macros.h"
#:include "../macros.fpph"
module guga_pchb_singles_weights_mod
    use util_mod, only: EnumBase_t, stop_all, pure_stop_all
    use SystemData, only: nel, nSpatOrbs
    use constants, only: dp
    use util_mod, only: EnumBase_t
    use fortran_strings, only: to_upper
    use integrals_neci, only: h, g

    better_implicit_none
    private
    public :: get_weight_t, PropVec_PCHB_SinglesWeighting_t, PropVec_PCHB_SinglesWeighting_vals_t, &
        set_get_weight_pointer

    abstract interface
        real(dp) function get_weight_t(i, a)
            import :: dp
            integer, intent(in) :: i, a
                !! Spatial orbital index
        end function
    end interface


    type, extends(EnumBase_t) :: PropVec_PCHB_SinglesWeighting_t
        character(20) :: str
    contains
        procedure :: to_str
    end type


    type :: PropVec_PCHB_SinglesWeighting_vals_t
        type(PropVec_PCHB_SinglesWeighting_t) :: &
            UNIFORM = PropVec_PCHB_SinglesWeighting_t(1, 'UNIFORM'), &
            ABS_INTEGRAL = PropVec_PCHB_SinglesWeighting_t(2, 'ABS-INTEGRAL'), &
            RHF_EXPR = PropVec_PCHB_SinglesWeighting_t(3, 'RHF-EXPRESSION')
    contains
        procedure, nopass :: from_str => weight_from_str
    end type

    type(PropVec_PCHB_SinglesWeighting_vals_t), parameter :: weighting_vals = PropVec_PCHB_SinglesWeighting_vals_t()



contains

    subroutine set_get_weight_pointer(weight_choice, ptr)
        type(PropVec_PCHB_SinglesWeighting_t), intent(in) :: weight_choice
        procedure(get_weight_t), pointer, intent(out) :: ptr
        routine_name("set_get_weight_pointer")

        if (weight_choice == weighting_vals%UNIFORM) then
            ptr => get_weight_uniform
        else if (weight_choice == weighting_vals%ABS_INTEGRAL) then
            ptr => get_weight_abs_integral
        else if (weight_choice == weighting_vals%RHF_EXPR) then
            ptr => get_weight_RHF_expr
        else
            call stop_all(this_routine, "Invalid weight choice: "//weight_choice%str)
        end if
    end subroutine

    elemental function weight_from_str(str) result(res)
        character(*), intent(in) :: str
        type(PropVec_PCHB_SinglesWeighting_t) :: res
        routine_name("alg_from_str")

        select case(to_upper(str))
        case(weighting_vals%UNIFORM%str)
            res = weighting_vals%UNIFORM
        case(weighting_vals%ABS_INTEGRAL%str)
            res = weighting_vals%ABS_INTEGRAL
        case(weighting_vals%RHF_EXPR%str)
            res = weighting_vals%RHF_EXPR
        case default
            call pure_stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

    pure function get_weight_uniform(i, a) result(weight)
        integer, intent(in) :: i, a
            !! Spatial orbital index
        real(dp) :: weight
        routine_name("get_weight_uniform")
        @:pure_ASSERT(all(1 <= [i, a] .and. [i, a] <= nSpatOrbs))
        weight = merge(1._dp, 0._dp, i /= a)
    end function


    pure function get_weight_abs_integral(i, a) result(weight)
        integer, intent(in) :: i, a
            !! Spatial orbital index
        real(dp) :: weight
        routine_name("get_weight_abs_integral")

        real(dp) :: two_el_term
        integer :: r
        @:pure_ASSERT(all(1 <= [i, a] .and. [i, a] <= nSpatOrbs))

        two_el_term = 0._dp
        do R = 1, nSpatOrbs
            ! Triangle inequality makes it non-zero, whenever integrals are non-zero:
            ! |a| + |b|  >=  | a + b |
            two_el_term = two_el_term + abs(g(i, a, r, r)) + abs(g(i, r, r, a))
        end do
        weight = abs(h(i, a)) + two_el_term
    end function

    pure function get_weight_RHF_expr(i, a) result(weight)
        integer, intent(in) :: i, a
            !! Spatial orbital index
        real(dp) :: weight
        routine_name("get_weight_uniform")

        real(dp) :: two_el_term
        integer :: r
        @:pure_ASSERT(all(1 <= [i, a] .and. [i, a] <= nSpatOrbs))

        two_el_term = 0._dp
        do R = 1, nSpatOrbs
            two_el_term = two_el_term + abs(2 * g(i, a, r, r) - g(i, r, r, a))
        end do
        weight = abs(h(i, a)) + two_el_term
    end function

    pure function to_str(weighting) result(res)
        class(PropVec_PCHB_SinglesWeighting_t), intent(in) :: weighting
        character(:), allocatable :: res
        routine_name("to_str")
        if (weighting == weighting_vals%UNIFORM) then
            res = weighting_vals%UNIFORM%str
        else if (weighting == weighting_vals%ABS_INTEGRAL) then
            res = weighting_vals%ABS_INTEGRAL%str
        else if (weighting == weighting_vals%RHF_EXPR) then
            res = weighting_vals%RHF_EXPR%str
        else
            call pure_stop_all(this_routine, "Should not be here.")
        end if
    end function

end module
