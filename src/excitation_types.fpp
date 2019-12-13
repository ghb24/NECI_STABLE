#include "macros.h"
#:include "macros.fpph"
#:set ExcitationTypes = ['NoExc_t', 'SingleExc_t', 'DoubleExc_t', 'TripleExc_t']

module excitation_types
    use constants, only: dp
    use SystemData, only: nEl
    implicit none
    private
    public :: excitation_t, NoExc_t, SingleExc_t, DoubleExc_t, TripleExc_t, &
        UNKNOWN, defined, get_excitation


!> Arbitrary non occuring (?!) orbital index.
    integer, parameter :: UNKNOWN = -10**5

    type :: excitation_t
    end type

!> Represents a No-Op excitation.
    type, extends(excitation_t) :: NoExc_t
    end type

!> Represents the orbital indices of a single excitation.
!> The array is sorted like:
!> [src1, tgt2]
    type, extends(excitation_t) :: SingleExc_t
        integer :: val(2) = UNKNOWN
    end type

!> Represents the orbital indices of a double excitation.
!> The array is sorted like:
!> [[src1, src2],
!>  [tgt1, tgt2]]
    type, extends(excitation_t) :: DoubleExc_t
        integer :: val(2, 2) = UNKNOWN
    end type

!> Represents the orbital indices of a triple excitation.
!> The array is sorted like:
!> [[src1, src2, src3],
!>  [tgt1, tgt2, tgt3]]
    type, extends(excitation_t) :: TripleExc_t
        integer :: val(2, 3) = UNKNOWN
    end type

!> Additional constructors for the excitation types from integers instead
!> of an integer array.
    #:for excitation_t in ExcitationTypes[1:]
    interface ${excitation_t}$
        module procedure from_integer_${excitation_t}$
    end interface
    #:endfor

!> Returns true if all sources and targets are not UNKNOWN.
    interface defined
    #:for excitation_t in ExcitationTypes
        module procedure defined_${excitation_t}$
    #:endfor
    end interface


contains

    #:for excitation_t in ExcitationTypes[1:]
        elemental function defined_${excitation_t}$(exc) result(res)
            type(${excitation_t}$), intent(in) :: exc
            logical :: res

            res = all(exc%val /= UNKNOWN)
        end function
    #:endfor

    elemental function defined_NoExc_t(exc) result(res)
        type(NoExc_t), intent(in) :: exc
        logical :: res

        res = .true.
    end function

    pure function from_integer_SingleExc_t(src, tgt) result(res)
        integer, intent(in), optional :: src, tgt
        type(SingleExc_t) :: res

        ! The values default to UNKNOWN, in the type.
        if (present(src)) res%val(1) = src
        if (present(tgt)) res%val(1) = tgt
    end function

    pure function from_integer_DoubleExc_t(src1, tgt1, src2, tgt2) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2
        type(DoubleExc_t) :: res

        ! The values default to UNKNOWN, in the type.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
    end function

    pure function from_integer_TripleExc_t(src1, tgt1, src2, tgt2, src3, tgt3) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2, src3, tgt3
        type(TripleExc_t) :: res

        ! The values default to UNKNOWN, in the type.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
        if (present(src2)) res%val(1, 3) = src3
        if (present(tgt2)) res%val(2, 3) = tgt3
    end function

    subroutine get_excitation(nI, nJ, IC, exc, tParity)
        integer, intent(in) :: nI(nEl), nJ(nEl), IC
        class(excitation_t), allocatable, intent(out) :: exc
        logical, intent(out) :: tParity

        select case (IC)
        case(0)
            allocate(NoExc_t :: exc)
        case(1)
            allocate(SingleExc_t :: exc)
        case(2)
            allocate(DoubleExc_t :: exc)
        case(3)
            allocate(TripleExc_t :: exc)
        case default
            allocate(exc)
        end select

        select type (exc)
        type is (SingleExc_t)
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (DoubleExc_t)
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (TripleExc_t)
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        end select
    end subroutine get_excitation
end module
