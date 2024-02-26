#include "macros.h"
#:include "../macros.fpph"
module guga_pchb_singles_weights_mod
    use util_mod, only: EnumBase_t, stop_all
    use SystemData, only: nel, nSpatOrbs
    use constants, only: dp
    use util_mod, only: EnumBase_t
    use error_handling_neci, only: stop_all
    use fortran_strings, only: to_upper
    use procedure_pointers, only: get_umat_el
    use OneEInts, only: GetTMatEl

    better_implicit_none
    private
    public :: get_weight_t, GAS_PCHB_SinglesWeighting_t, GAS_PCHB_SinglesWeighting_vals_t, &
        set_get_weight_pointer

    abstract interface
        real(dp) function get_weight_t(i, a)
            import :: dp
            integer, intent(in) :: i, a
                !! Spatial orbital index
        end function
    end interface


    type, extends(EnumBase_t) :: GAS_PCHB_SinglesWeighting_t
        character(20) :: str
    end type


    type :: GAS_PCHB_SinglesWeighting_vals_t
        type(GAS_PCHB_SinglesWeighting_t) :: &
            UNIFORM = GAS_PCHB_SinglesWeighting_t(1, 'UNIFORM'), &
            ABS_INTEGRAL = GAS_PCHB_SinglesWeighting_t(2, 'ABS-INTEGRAL'), &
            RHF_EXPR = GAS_PCHB_SinglesWeighting_t(3, 'RHF-EXPRESSION')
    contains
        procedure, nopass :: from_str => weight_from_str
    end type

    type(GAS_PCHB_SinglesWeighting_vals_t), parameter :: weighting_vals = GAS_PCHB_SinglesWeighting_vals_t()



contains

    subroutine set_get_weight_pointer(weight_choice, ptr)
        type(GAS_PCHB_SinglesWeighting_t), intent(in) :: weight_choice
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
        type(GAS_PCHB_SinglesWeighting_t) :: res
        routine_name("alg_from_str")

        select case(to_upper(str))
        case(weighting_vals%UNIFORM%str)
            res = weighting_vals%UNIFORM
        case(weighting_vals%ABS_INTEGRAL%str)
            res = weighting_vals%ABS_INTEGRAL
        case(weighting_vals%RHF_EXPR%str)
            res = weighting_vals%RHF_EXPR
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
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


    HElement_t(dp) elemental function h(i, a)
        !! Return the 1el integral \( h_{I, A) \)
        !!
        !! I and A are **spatial** indices.
        !! Follows the definition of the purple book.
        integer, intent(in) :: i, a
        h = GetTMATEl(2 * i, 2 * a)
    end function


    HElement_t(dp) elemental function g(i, a, j, b)
        integer, intent(in) :: i, a, j, b
        !! Return the 2el integral \( g_{I, A, J, B} \)

        !! i, a, j, and b are **spatial ** indices.
        !! Order follows the definition of the purple book (chemist's notation).
        !! \( g_{I, A, J, B} = < I J | 1 / r | A B > \)
        g = get_umat_el(i, j, a, b)
    end function
end module
