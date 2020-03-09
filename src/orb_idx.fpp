#include "macros.h"
#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']

module orb_idx_mod
    implicit none
    private
    public :: OrbIdx_t, SpinOrbIdx_t, SpatOrbIdx_t, size, &
        Spin_t, calc_spin, calc_spin_raw, spin_values, operator(==)

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

    type :: SpinValues_t
        integer :: alpha = 1, beta = -1
    end type

    type(SpinValues_t), parameter :: spin_values = SpinValues_t()

    type :: Spin_t
        integer, allocatable :: m_s(:)
    end type

    interface Spin_t
        module procedure construction_from_number_Spin_t
    end interface

    interface size
#:for type in OrbIdxTypes + ['Spin_t']
        module procedure size_${type}$
#:endfor
    end interface

    interface operator(==)
#:for type in OrbIdxTypes
        module procedure eq_${type}$
#:endfor
    end interface

    contains

    pure function construction_from_array_SpinOrbIdx_t(idx, spin) result(res)
        integer, intent(in) :: idx(:)
        type(Spin_t), intent(in), optional :: spin

        type(SpinOrbIdx_t) :: res

        if (present(spin)) then
            if (size(spin) == 1) then
                res%idx = pack(idx, calc_spin_raw(idx) == spin%m_s(1))
            else
                res%idx = pack(idx, calc_spin_raw(idx) == spin%m_s(:))
            end if
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

    pure function SpinOrbIdx_t_from_SpatOrbIdx_t(spat_orbs, spins) result(res)
        type(SpatOrbIdx_t), intent(in) :: spat_orbs
        type(Spin_t), intent(in), optional :: spins
        type(SpinOrbIdx_t) :: res

        character(*), parameter :: this_routine = 'SpinOrbIdx_t_from_SpatOrbIdx_t'

        if (present(spins)) then
            if (size(spins) == 1) then
                res%idx = f(spat_orbs%idx(:), spins%m_s(1))
            else
                res%idx = f(spat_orbs%idx(:), spins%m_s(:))
            end if
        else
            allocate(res%idx(2 * size(spat_orbs)))
            res%idx(1::2) = f(spat_orbs%idx(:), spin_values%beta)
            res%idx(2::2) = f(spat_orbs%idx(:), spin_values%alpha)
        end if
        contains
            elemental function f(spat_orb_idx, spin_idx) result(spin_orb_idx)
                integer, intent(in) :: spat_orb_idx, spin_idx
                integer :: spin_orb_idx
                if (spin_idx == spin_values%alpha) then
                    spin_orb_idx = 2 * spat_orb_idx
                else
                    spin_orb_idx = 2 * spat_orb_idx - 1
                end if
            end function
    end function

    pure function construction_from_number_Spin_t(spin) result(res)
        integer, intent(in) :: spin
        type(Spin_t) :: res
        res%m_s = [spin]
    end function

    pure function size_Spin_t(spins) result(res)
        type(Spin_t), intent(in) :: spins
        integer :: res
        res = size(spins%m_s)
    end function

    pure function calc_spin(orbs) result(res)
        type(SpinOrbIdx_t), intent(in) :: orbs
        type(Spin_t) :: res
        integer :: i
        res%m_s = calc_spin_raw(orbs%idx)
    end function

    integer elemental function calc_spin_raw(orb_idx)
        integer, intent(in) :: orb_idx
        calc_spin_raw = merge(spin_values%alpha, spin_values%beta, mod(orb_idx, 2) == 0)
    end function
end module
