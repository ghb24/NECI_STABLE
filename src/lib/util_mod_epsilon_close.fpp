#include "macros.h"
#:set primitive_types = {'real': {'sp', 'dp'}, 'complex': {'sp', 'dp'}}

module util_mod_epsilon_close
    use constants, only : sp, dp, EPS
    implicit none

    private
    public :: isclose, operator(.isclose.), near_zero

    interface operator(.isclose.)
    #:for type, kinds in primitive_types.items()
    #:for kind in kinds
        module procedure isclose_for_operator_${type}$_${kind}$
    #:endfor
    #:endfor
    end interface

    interface isclose
    #:for type, kinds in primitive_types.items()
    #:for kind in kinds
        module procedure isclose_${type}$_${kind}$
    #:endfor
    #:endfor
    end interface

    interface near_zero
    #:for type, kinds in primitive_types.items()
    #:for kind in kinds
        module procedure near_zero_${type}$_${kind}$
    #:endfor
    #:endfor
    end interface

contains
    #:for type, kinds in primitive_types.items()
    #:for kind in kinds

        elemental function isclose_${type}$_${kind}$(a, b, atol, rtol) result(res)
            ${type}$(${kind}$), intent(in) :: a, b
            real(dp), intent(in), optional :: atol, rtol
            logical :: res

            real(dp) :: atol_, rtol_

            def_default(rtol_, rtol, 1d-5)
            def_default(atol_, atol, EPS)

            ! Change to <= max(atol_, rtol_ * max(abs(a), abs(b)))
            res = abs(a - b) <= (atol_ + rtol_ * abs(b))
        end function

    !> Operator functions may only have two arguments.
        elemental function isclose_for_operator_${type}$_${kind}$(a, b) result(res)
            ${type}$(${kind}$), intent(in) :: a, b
            logical :: res

            res = isclose(a, b)
        end function

        elemental function near_zero_${type}$_${kind}$(x, epsilon) result(res)
            ${type}$(${kind}$), intent(in) :: x
            real(dp), intent(in), optional :: epsilon
            logical :: res

            real(dp) :: epsilon_

            def_default(epsilon_, epsilon, EPS)

            res = abs(x) < epsilon_
        end function

    #:endfor
    #:endfor
end module
