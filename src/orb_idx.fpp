#include "macros.h"
#:include "macros.fpph"
#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']

module orb_idx_mod
    implicit none
    private
    public :: OrbIdx_t, SpinOrbIdx_t, SpatOrbIdx_t, size, &
        SpinProj_t, calc_spin, calc_spin_raw, alpha, beta, &
        operator(==), operator(/=), operator(+), operator(-), &
        sum

    type, abstract :: OrbIdx_t
        integer, allocatable :: idx(:)
    end type

    !> We assume order [beta_1, alpha_1, beta_2, alpha_2, ...]
    type, extends(OrbIdx_t) :: SpinOrbIdx_t
    end type

    interface SpinOrbIdx_t
        module procedure SpinOrbIdx_t_from_SpatOrbIdx_t
        module procedure construction_from_array_SpinOrbIdx_t
    end interface

    type, extends(OrbIdx_t) :: SpatOrbIdx_t
    end type

    type :: SpinProj_t
        integer :: val
    end type

    type(SpinProj_t), parameter :: beta = SpinProj_t(-1), alpha = SpinProj_t(1)

    interface size
#:for type in OrbIdxTypes
        module procedure size_${type}$
#:endfor
    end interface

    interface operator(==)
#:for type in OrbIdxTypes
        module procedure eq_${type}$
#:endfor
        module procedure eq_SpinProj_t_SpinProj_t
    end interface

    interface operator(/=)
#:for type in OrbIdxTypes
        module procedure neq_${type}$
#:endfor
        module procedure neq_SpinProj_t_SpinProj_t
    end interface

    interface operator (+)
        module procedure add_SpinProj_t_SpinProj_t
    end interface

    interface operator (-)
        module procedure sub_SpinProj_t_SpinProj_t
        module procedure neg_SpinProj_t
    end interface

    interface sum
        module procedure sum_SpinProj_t
    end interface

    contains

    DEBUG_IMPURE function construction_from_array_SpinOrbIdx_t(idx, m_s) result(res)
        integer, intent(in) :: idx(:)
        type(SpinProj_t), intent(in), optional :: m_s
        type(SpinOrbIdx_t) :: res
        character(*), parameter :: this_routine = 'construction_from_array_SpinOrbIdx_t'

        if (present(m_s)) then
            @:ASSERT(any(m_s == [alpha, beta]))
            associate(spins => calc_spin_raw(idx))
                res%idx = pack(idx, spins%val == m_s%val)
            end associate
        else
            res%idx = idx
        end if
    end function

#:for orb_idx_type in OrbIdxTypes
    pure function size_${orb_idx_type}$(orbs) result(res)
        type(${orb_idx_type}$), intent(in) :: orbs
        integer :: res
        res = size(orbs%idx)
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    pure function eq_${orb_idx_type}$(lhs, rhs) result(res)
        type(${orb_idx_type}$), intent(in) :: lhs, rhs
        logical :: res(size(lhs))
        res = lhs%idx == rhs%idx
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    pure function neq_${orb_idx_type}$(lhs, rhs) result(res)
        type(${orb_idx_type}$), intent(in) :: lhs, rhs
        logical :: res(size(lhs))
        res = lhs%idx /= rhs%idx
    end function
#:endfor

    DEBUG_IMPURE function SpinOrbIdx_t_from_SpatOrbIdx_t(spat_orbs, m_s) result(res)
        type(SpatOrbIdx_t), intent(in) :: spat_orbs
        type(SpinProj_t), intent(in), optional :: m_s
        type(SpinOrbIdx_t) :: res

        character(*), parameter :: this_routine = 'SpinOrbIdx_t_from_SpatOrbIdx_t'


        if (present(m_s)) then
            @:ASSERT(any(m_s == [alpha, beta]))
            res%idx = f(spat_orbs%idx(:), m_s)
        else
            allocate(res%idx(2 * size(spat_orbs)))
            res%idx(1::2) = f(spat_orbs%idx(:), beta)
            res%idx(2::2) = f(spat_orbs%idx(:), alpha)
        end if
        contains
            elemental function f(spat_orb_idx, m_s) result(spin_orb_idx)
                integer, intent(in) :: spat_orb_idx
                type(SpinProj_t), intent(in) :: m_s
                integer :: spin_orb_idx
                if (m_s == alpha) then
                    spin_orb_idx = 2 * spat_orb_idx
                else
                    spin_orb_idx = 2 * spat_orb_idx - 1
                end if
            end function
    end function

    pure function calc_spin(orbs) result(res)
        type(SpinOrbIdx_t), intent(in) :: orbs
        type(SpinProj_t) :: res(size(orbs))
        res = calc_spin_raw(orbs%idx)
    end function

    elemental function add_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        type(SpinProj_t), intent(in) :: lhs, rhs
        type(SpinProj_t) :: res
        res%val = lhs%val + rhs%val
    end function

    pure function sum_SpinProj_t(V) result(res)
        type(SpinProj_t), intent(in) :: V(:)
        type(SpinProj_t) :: res
        integer :: i
        res = SpinProj_t(0)
        do i = 1, size(V)
            res = res + V(i)
        end do
    end function

    elemental function sub_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        type(SpinProj_t), intent(in) :: lhs, rhs
        type(SpinProj_t) :: res
        res%val = lhs%val - rhs%val
    end function

    elemental function neg_SpinProj_t(m_s) result(res)
        type(SpinProj_t), intent(in) :: m_s
        type(SpinProj_t) :: res
        res%val = -m_s%val
    end function

    elemental function eq_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        type(SpinProj_t), intent(in) :: lhs, rhs
        logical :: res
        res = lhs%val == rhs%val
    end function

    elemental function neq_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        type(SpinProj_t), intent(in) :: lhs, rhs
        logical :: res
        res = lhs%val == rhs%val
    end function

    elemental function calc_spin_raw(orb_idx) result(res)
        integer, intent(in) :: orb_idx
        type(SpinProj_t) :: res
        res = merge(alpha, beta, mod(orb_idx, 2) == 0)
    end function
end module
