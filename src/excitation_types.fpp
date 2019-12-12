#include "macros.h"
#:include "macros.fpph"
#:set ExcitationTypes = ['NoExc_t', 'SingleExc_t', 'DoubleExc_t', 'TripleExc_t']

module excitation_types
    use constants, only: dp
    implicit none
    private
    public :: NoExc_t, SingleExc_t, DoubleExc_t, TripleExc_t, UNKNOWN, &
        defined


!> Arbitrary non occuring (?!) orbital index.
    integer, parameter :: UNKNOWN = -10**5

!> Represents a No-Op excitation.
    type :: NoExc_t
    end type

!> Represents the orbital indices of a single excitation.
!> The array is sorted like:
!> [src1, tgt2]
    type :: SingleExc_t
        integer :: val(2) = UNKNOWN
    end type

!> Represents the orbital indices of a double excitation.
!> The array is sorted like:
!> [[src1, src2],
!>  [tgt1, tgt2]]
    type :: DoubleExc_t
        integer :: val(2, 2) = UNKNOWN
    end type

!> Represents the orbital indices of a triple excitation.
!> The array is sorted like:
!> [[src1, src2, src3],
!>  [tgt1, tgt2, tgt3]]
    type :: TripleExc_t
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
        type(DoubleExc_t) :: res

        ! The values default to UNKNOWN, in the type.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
        if (present(src2)) res%val(1, 3) = src3
        if (present(tgt2)) res%val(2, 3) = tgt3
    end function

end module
