#include "macros.h"
#:include "macros.fpp"
#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module new_sltcnd_mod
    use constants, only: dp
    use SystemData, only: nEl
    use sltcnd_mod, only: sltcnd_1, sltcnd_2
    implicit none

    private
    public :: SingleExc_t, DoubleExc_t, UNKNOWN, sltcnd_excit, defined

    integer, parameter :: UNKNOWN = -10**5

!> Represents the orbital indices of a single excitation.
!> The array is sorted like:
!> [src1, tgt2]
    type :: SingleExc_t
        sequence
        integer :: val(2) = [UNKNOWN, UNKNOWN]
    end type

!> Represents the orbital indices of a double excitation.
!> The array is sorted like:
!> [[src1, src2],
!>  [tgt1, tgt2]]
    type :: DoubleExc_t
        sequence
        integer :: val(2, 2) = reshape([UNKNOWN, UNKNOWN, UNKNOWN, UNKNOWN], [2, 2])
    end type

    #:for excitation_t in ExcitationTypes
    interface ${excitation_t}$
        module procedure from_integer_${excitation_t}$
    end interface
    #:endfor

    interface defined
    #:for excitation_t in ExcitationTypes
        module procedure defined_${excitation_t}$
    #:endfor
    end interface

    interface sltcnd_excit
    #:for excitation_t in ExcitationTypes
        module procedure sltcnd_excit_${excitation_t}$
    #:endfor
    end interface

contains

    #:for excitation_t in ExcitationTypes
    elemental function defined_${excitation_t}$(exc) result(res)
        type(${excitation_t}$), intent(in) :: exc
        logical :: res
        res = all(exc%val /= UNKNOWN)
    end function
    #:endfor

    HElement_t(dp) function sltcnd_excit_SingleExc_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(SingleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity

        sltcnd_excit_SingleExc_t = sltcnd_1(ref, exc%val, tParity)
    end function

    HElement_t(dp) function sltcnd_excit_DoubleExc_t(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(DoubleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity

        @:unused_var(ref)

        sltcnd_excit_DoubleExc_t = sltcnd_2(exc%val, tParity)
    end function

    pure function from_integer_SingleExc_t(src, tgt) result(res)
        integer, intent(in), optional :: src, tgt
        type(SingleExc_t) :: res
        if (present(src)) res%val(1) = src
        if (present(tgt)) res%val(1) = tgt
    end function

    pure function from_integer_DoubleExc_t(src1, tgt1, src2, tgt2) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2
        type(DoubleExc_t) :: res
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(1, 2) = tgt1
        if (present(src2)) res%val(2, 1) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
    end function

end module new_sltcnd_mod
