#include "macros.h"
module guga_pchb_doubles_weights_mod
    use SystemData, only: nel, nSpatOrbs, t_exchange_pchb
    use constants, only: dp
    use util_mod, only: EnumBase_t
    use guga_data, only: DistinctDouble_t, DistinctDoubleTypeVals_t, &
        distinct_doubles
    use error_handling_neci, only: stop_all
    use fortran_strings, only: to_upper
    use integrals_neci, only: g

    better_implicit_none
    private
    public :: WeightingChoice_t, get_doubles_weight_t, WeightingChoice_vals_t, &
        set_doubles_weight_pointer


    abstract interface
        real(dp) function get_doubles_weight_t(exc)
            import :: DistinctDouble_t, dp
            type(DistinctDouble_t), intent(in) :: exc
        end function
    end interface


    type, extends(EnumBase_t) :: WeightingChoice_t
        character(20) :: str
    contains
        procedure :: to_str => to_str_WeightingChoice_t
    end type

    type :: WeightingChoice_vals_t
        type(WeightingChoice_t) :: &
            INITIAL_DOBRAUTZ = WeightingChoice_t(1, 'INITIAL'), &
            UNIFORM = WeightingChoice_t(2, 'UNIFORM'), &
            ABS_INTEGRAL = WeightingChoice_t(3, 'ABS-INTEGRAL'), &
            N4_TETRAHEDRON = WeightingChoice_t(4, 'N4-TETRAHEDRON'), &
            RANDOM = WeightingChoice_t(5, 'RANDOM'), &
            TWEAKED_WERNER = WeightingChoice_t(6, 'TWEAKED-WERNER'), &
            NEW_WERNER = WeightingChoice_t(7, 'NEW-WERNER')
    contains
        procedure, nopass :: from_str => WeightingChoice_from_str
    end type

    type(WeightingChoice_vals_t), parameter :: weight_choice_vals = WeightingChoice_vals_t()

contains

    subroutine set_doubles_weight_pointer(weight_choice, ptr)
        type(WeightingChoice_t), intent(in) :: weight_choice
        procedure(get_doubles_weight_t), pointer, intent(out) :: ptr
        routine_name("set_get_weight_pointer")

        if (weight_choice == weight_choice_vals%INITIAL_DOBRAUTZ) then
            ptr => initial_werner_get_weight
        else if (weight_choice == weight_choice_vals%TWEAKED_WERNER) then
            ptr => tweaked_werner_get_weight
        else if (weight_choice == weight_choice_vals%NEW_WERNER) then
            ptr => new_werner_get_weight
        else if (weight_choice == weight_choice_vals%UNIFORM) then
            ptr => uniform_get_weight
        else if (weight_choice == weight_choice_vals%ABS_INTEGRAL) then
            ptr => abs_integral_get_weight
        else if (weight_choice == weight_choice_vals%N4_TETRAHEDRON) then
            ptr => N4_tetrahedron_get_weight
        else if (weight_choice == weight_choice_vals%RANDOM) then
            ptr => random_get_weight
        else
            call stop_all(this_routine, "Invalid weight choice: "//weight_choice%str)
        end if
    end subroutine

    pure function WeightingChoice_from_str(str) result(res)
        character(*), intent(in) :: str
        type(WeightingChoice_t) :: res
        routine_name("WeightingChoice_from_str")

        select case(to_upper(str))
        case(weight_choice_vals%INITIAL_DOBRAUTZ%str)
            res = weight_choice_vals%INITIAL_DOBRAUTZ
        case(weight_choice_vals%TWEAKED_WERNER%str)
            res = weight_choice_vals%TWEAKED_WERNER
        case(weight_choice_vals%NEW_WERNER%str)
            res = weight_choice_vals%NEW_WERNER
        case(weight_choice_vals%UNIFORM%str)
            res = weight_choice_vals%UNIFORM
        case(weight_choice_vals%ABS_INTEGRAL%str)
            res = weight_choice_vals%ABS_INTEGRAL
        case(weight_choice_vals%N4_TETRAHEDRON%str)
            res = weight_choice_vals%N4_TETRAHEDRON
        case(weight_choice_vals%RANDOM%str)
            res = weight_choice_vals%RANDOM
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

    pure function initial_werner_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable

        integer :: i, j, k, l
        HElement_t(dp) :: cpt1, cpt2, cpt3, cpt4

        call exc%idx%extract_ijkl(i, j, k, l)

        select case (exc%typ%val)
        ! 0 el transferred
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val)
            weight = abs(g(i, j, k, l) + g(i, l, k, j)) / 2.0_dp
        ! 0 el transferred
        case(dd%ijij%val, dd%jiji%val)
            weight = abs(g(i, j, k, l)) / 2.0_dp
        ! 0 el transferred
        case (dd%iijj%val)
            ! weight = abs(g(i, l, k, i))
            weight = abs(g(i, l, l, i))
        ! 1 el transferred
        case(dd%jikk%val, dd%ijkk%val)
            weight = abs(g(i, k, k, j))
        ! 1 el transferred
        case(dd%iikj%val, dd%iijk%val)
            weight = abs(g(i, l, k, i))
        ! 2 el transferred
        case(dd%kilj%val, dd%ikjl%val)
            cpt1 = 2 * g(i, j, k, l)
            cpt2 = 2 * g(i, l, k, j)

            cpt3 = abs(cpt1 + cpt2) / 2.0_dp
            cpt4 = abs(cpt1 - cpt2) / 2.0_dp

            if (t_exchange_pchb) then
                weight = abs(cpt4)
            else
                weight = maxval(abs([cpt1, cpt2, cpt3, cpt4]))
            end if
        ! 2 el transferred
        case(dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            cpt1 = 2 * g(i, j, k, l)
            cpt2 = 2 * g(i, l, k, j)

            cpt3 = abs(cpt2 - cpt1)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1) / 2.0_dp

            if (t_exchange_pchb) then
                weight = abs(cpt2)
            else
                weight = maxval(abs([cpt1, cpt2 / 2.0_dp, cpt3, cpt4]))
            end if
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    pure function new_werner_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable

        integer :: i, j, k, l

        call exc%idx%extract_ijkl(i, j, k, l)

        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val)
            weight = abs(g(i, j, k, l)) + abs(g(i, l, k, j))
        case(dd%ijij%val, dd%jiji%val)
            weight = abs(g(i, j, k, l)) + 1e-1_dp
        case(dd%jikk%val, dd%ijkk%val)
            weight = abs(g(i, k, k, j))
        case(dd%iikj%val, dd%iijk%val, dd%iijj%val)
            weight = abs(g(i, l, k, i))
        case(dd%kilj%val, dd%ikjl%val, dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = 2 * (abs(g(i, j, k, l)) + abs(g(i, l, k, j)))
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    pure function tweaked_werner_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable

        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = initial_werner_get_weight(exc) + 1e-4_dp
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    function uniform_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("uniform_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = 1._dp
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    function random_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("uniform_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = 1._dp + get_r() * 100
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select

    contains
        ! We cannot use the own random number generator of NECI,
        ! it is not yet initialized
        real(dp) function get_r()
            call random_number(get_r)
        end function
    end function

    function abs_integral_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("abs_integral_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        integer :: i, j, k, l

        call exc%idx%extract_ijkl(i, j, k, l)
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = abs(g(i, j, k, l)) + abs(g(i, l, k, j))
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    function N4_tetrahedron_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("N4_tetrahedron_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        integer :: i, j, k, l

        if (nSpatOrbs /= 12) then
            call stop_all(this_routine, "N4 tetrahedron required")
        end if

        call exc%idx%extract_ijkl(i, j, k, l)
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = merge(3, 1, mod(i - 1, 3) == mod(j - 1, 3)) * merge(3, 1, mod(k - 1, 3) == mod(l - 1, 3))
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function



    pure function to_str_WeightingChoice_t(this) result(res)
        class(WeightingChoice_t), intent(in) :: this
        character(:), allocatable :: res
        res = trim(this%str)
    end function

end module
