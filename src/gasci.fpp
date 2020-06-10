module gasci
    implicit none

    private
    public :: operator(==), operator(/=), possible_GAS_exc_gen, GAS_exc_gen


    type :: GAS_exc_gen_t
        integer :: val
    end type

    type :: possible_GAS_exc_gen_t
        type(GAS_exc_gen_t) :: &
            DISCONNECTED = GAS_exc_gen_t(1), &
            GENERAL = GAS_exc_gen_t(2)
    end type

    type(possible_GAS_exc_gen_t), parameter :: possible_GAS_exc_gen = possible_GAS_exc_gen_t()

    type(GAS_exc_gen_t) :: GAS_exc_gen = possible_GAS_exc_gen%GENERAL

    interface operator(==)
        module procedure eq_GAS_exc_gen_t
    end interface

    interface operator(/=)
        module procedure neq_GAS_exc_gen_t
    end interface

contains

    logical pure function eq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        eq_GAS_exc_gen_t = lhs%val == rhs%val
    end function

    logical pure function neq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        neq_GAS_exc_gen_t = lhs%val /= rhs%val
    end function
end module gasci
