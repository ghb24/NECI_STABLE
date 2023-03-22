#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']

module orb_idx_mod
    use constants, only: n_int, stdout
    use fortran_strings, only: str
    use bit_rep_data, only: nIfTot, nIfD
    use bit_reps, only: decode_bit_det
    use DetBitOps, only: EncodeBitDet
    use util_mod, only: ilex_leq => lex_leq, ilex_geq => lex_geq, stop_all, &
        operator(.div.)
    implicit none
    private
    public :: OrbIdx_t, SpinOrbIdx_t, SpatOrbIdx_t, size, &
              SpinProj_t, calc_spin, calc_spin_raw, alpha, beta, &
              operator(==), operator(/=), operator(+), operator(-), &
              sum, to_ilut, lex_leq, lex_geq, write_det, get_spat

    type, abstract :: OrbIdx_t
        integer, allocatable :: idx(:)
    end type

    !> We assume order [beta_1, alpha_1, beta_2, alpha_2, ...]
    type, extends(OrbIdx_t) :: SpinOrbIdx_t
    contains
        procedure, nopass :: from_ilut => from_ilut_SpinOrbIdx_t
    end type

    interface SpinOrbIdx_t
        module procedure SpinOrbIdx_t_from_SpatOrbIdx_t
        module procedure construction_from_array_SpinOrbIdx_t
    end interface

    type, extends(OrbIdx_t) :: SpatOrbIdx_t
    end type

    type :: SpinProj_t
        integer :: val
            !! Twice the spin projection as integer. \( S_z = 2 \cdot \text{val} \)
    contains
        private
        procedure :: eq_SpinProj_t_SpinProj_t
        generic, public :: operator(==) => eq_SpinProj_t_SpinProj_t
        procedure :: neq_SpinProj_t_SpinProj_t
        generic, public :: operator(/=) => neq_SpinProj_t_SpinProj_t
        procedure :: add_SpinProj_t_SpinProj_t
        generic, public :: operator(+) => add_SpinProj_t_SpinProj_t
        procedure :: sub_SpinProj_t_SpinProj_t
        procedure :: neg_SpinProj_t
        generic, public :: operator(-) => sub_SpinProj_t_SpinProj_t, neg_SpinProj_t
    end type

    type(SpinProj_t), parameter :: beta = SpinProj_t(-1), alpha = SpinProj_t(1)

    interface size
#:for type in OrbIdxTypes
        module procedure size_${type}$
#:endfor
    end interface

    interface write_det
#:for type in OrbIdxTypes
        module procedure write_det_${type}$
#:endfor
    end interface

    interface operator(==)
#:for type in OrbIdxTypes
        module procedure eq_${type}$
#:endfor
    end interface

    interface operator(/=)
#:for type in OrbIdxTypes
        module procedure neq_${type}$
#:endfor
    end interface

#:for type in OrbIdxTypes
    interface lex_leq
        module procedure lex_leq_${type}$
    end interface

    interface lex_geq
        module procedure lex_geq_${type}$
    end interface
#:endfor

    interface sum
        module procedure sum_SpinProj_t
    end interface

contains

    pure function construction_from_array_SpinOrbIdx_t(idx, m_s) result(res)
        integer, intent(in) :: idx(:)
        type(SpinProj_t), intent(in), optional :: m_s
        type(SpinOrbIdx_t) :: res
        type(SpinProj_t) :: spins(size(idx))
        character(*), parameter :: this_routine = 'construction_from_array_SpinOrbIdx_t'

        if (present(m_s)) then
            @:pure_ASSERT(any(m_s == [alpha, beta]))
            spins = calc_spin_raw(idx)
            res%idx = pack(idx, spins == m_s)
        else
            res%idx = idx
        end if
    end function

#:for orb_idx_type in OrbIdxTypes
    pure function size_${orb_idx_type}$ (orbs) result(res)
        type(${orb_idx_type}$), intent(in) :: orbs
        integer :: res
        res = size(orbs%idx)
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    pure function eq_${orb_idx_type}$ (lhs, rhs) result(res)
        type(${orb_idx_type}$), intent(in) :: lhs, rhs
        logical :: res(size(lhs))
        res = lhs%idx == rhs%idx
    end function
#:endfor

#:for orb_idx_type in OrbIdxTypes
    pure function neq_${orb_idx_type}$ (lhs, rhs) result(res)
        type(${orb_idx_type}$), intent(in) :: lhs, rhs
        logical :: res(size(lhs))
        res = lhs%idx /= rhs%idx
    end function
#:endfor

    pure function SpinOrbIdx_t_from_SpatOrbIdx_t(spat_orbs, m_s) result(res)
        type(SpatOrbIdx_t), intent(in) :: spat_orbs
        type(SpinProj_t), intent(in), optional :: m_s
        type(SpinOrbIdx_t) :: res

        character(*), parameter :: this_routine = 'SpinOrbIdx_t_from_SpatOrbIdx_t'

        if (present(m_s)) then
            @:pure_ASSERT(any(m_s == [alpha, beta]))
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

    pure function to_ilut(det_I) result(ilut)
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(kind=n_int) :: iLut(0:nIfTot)
        call EncodeBitDet(det_I%idx, ilut)
    end function

    pure function calc_spin(orbs) result(res)
        type(SpinOrbIdx_t), intent(in) :: orbs
        type(SpinProj_t) :: res(size(orbs))
        res = calc_spin_raw(orbs%idx)
    end function

    elemental function add_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        class(SpinProj_t), intent(in) :: lhs, rhs
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
        class(SpinProj_t), intent(in) :: lhs, rhs
        type(SpinProj_t) :: res
        res%val = lhs%val - rhs%val
    end function

    elemental function neg_SpinProj_t(m_s) result(res)
        class(SpinProj_t), intent(in) :: m_s
        type(SpinProj_t) :: res
        res%val = -m_s%val
    end function

    elemental function eq_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        class(SpinProj_t), intent(in) :: lhs, rhs
        logical :: res
        res = lhs%val == rhs%val
    end function

    elemental function neq_SpinProj_t_SpinProj_t(lhs, rhs) result(res)
        class(SpinProj_t), intent(in) :: lhs, rhs
        logical :: res
        res = lhs%val /= rhs%val
    end function

    elemental function calc_spin_raw(orb_idx) result(res)
        integer, intent(in) :: orb_idx
        type(SpinProj_t) :: res
        res = merge(alpha, beta, mod(orb_idx, 2) == 0)
    end function

#:for type in OrbIdxTypes
    pure function lex_leq_${type}$ (lhs, rhs) result(res)
        type(${type}$), intent(in) :: lhs, rhs
        logical :: res
        character(*), parameter :: this_routine = 'lex_lt_${type}$'

        @:pure_ASSERT(size(lhs) == size(rhs))
        res = ilex_leq(lhs%idx, rhs%idx)
    end function

    pure function lex_geq_${type}$ (lhs, rhs) result(res)
        type(${type}$), intent(in) :: lhs, rhs
        logical :: res
        character(*), parameter :: this_routine = 'lex_gt_${type}$'

        @:pure_ASSERT(size(lhs) == size(rhs))
        res = ilex_geq(lhs%idx, rhs%idx)
    end function
#:endfor

#:for type in OrbIdxTypes
    subroutine write_det_${type}$ (det_I, i_unit, advance)
        type(${type}$), intent(in) :: det_I
        integer, intent(in), optional :: i_unit
        logical, intent(in), optional :: advance

        integer :: i, i_unit_
        character(:), allocatable :: advance_str, format

        @:def_default(i_unit_, i_unit, stdout)

        if (present(advance)) then
            if (advance) then
                advance_str = 'yes'
            else
                advance_str = 'no'
            end if
        else
            advance_str = 'yes'
        end if

        write(i_unit_, "(a)", advance='no') '${type}$(['

        if (size(det_I) == 0) then
            write(i_unit_, "(a)", advance=advance_str) '])'
        else
            format = "(I"//str(int(log10(real(maxval(det_I%idx)))) + 2)//", a)"

            do i = 1, size(det_I) - 1
                write(i_unit_, format, advance='no') det_I%idx(i), ','
            end do
            write(i_unit_, format, advance=advance_str) det_I%idx(size(det_I)), '])'
        end if
    end subroutine
#:endfor

    pure function from_ilut_SpinOrbIdx_t(ilut) result(res)
        integer(n_int), intent(in) :: ilut(0:nIfD)
        type(SpinOrbIdx_t) :: res
        integer :: n_el
        n_el = sum(popCnt(ilut))
        allocate(res%idx(n_el))
        call decode_bit_det(res%idx, ilut)
    end function


    elemental function get_spat(iorb) result(res)
        !! Return the spatial orbital of iorb
        integer, intent(in) :: iorb
        integer :: res
        res = (iorb + 1) .div. 2
    end function
end module
