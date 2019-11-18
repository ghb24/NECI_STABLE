#include "macros.h"

module new_sltcnd_mod
    use constants, only: dp
    use SystemData, only: nEl
    use sltcnd_mod, only: sltcnd_1, sltcnd_2
    implicit none

    private
    public :: SingleExc_t, DoubleExc_t, sltcnd_excit, UNKNOWN, &
        defined, matrix

    integer, parameter :: UNKNOWN = 0

    type :: SingleExc_t
        sequence
        integer :: src = UNKNOWN, tgt = UNKNOWN
    end type

    type :: DoubleExc_t
        sequence
        integer :: src1 = UNKNOWN, tgt1 = UNKNOWN
        integer :: src2 = UNKNOWN, tgt2 = UNKNOWN
    end type

    interface defined
        module procedure s_defined, d_defined
    end interface

    interface matrix
        module procedure s_matrix, d_matrix
    end interface

    interface sltcnd_excit
        module procedure s_sltcnd_excit, d_sltcnd_excit
    end interface

contains

    logical elemental function s_defined(exc)
        type(SingleExc_t), intent(in) :: exc
        s_defined = exc%src /= UNKNOWN .and. exc%tgt /= UNKNOWN
    end function

    logical elemental function d_defined(exc)
        type(DoubleExc_t), intent(in) :: exc
        d_defined = exc%src1 /= UNKNOWN .and. exc%tgt1 /= UNKNOWN &
                    .and. exc%src2 /= UNKNOWN .and. exc%tgt2 /= UNKNOWN
    end function

    pure function s_matrix(exc) result(res)
        type(SingleExc_t), intent(in) :: exc
        integer :: res(2)
        res(1) = exc%src
        res(2) = exc%tgt
    end function

    pure function d_matrix(exc) result(res)
        type(DoubleExc_t), intent(in) :: exc
        integer :: res(2, 2)
        res(1, 1) = exc%src1; res(1, 2) = exc%src2
        res(2, 1) = exc%tgt1; res(2, 2) = exc%tgt2
    end function

    HElement_t(dp) function s_sltcnd_excit(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(SingleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity

        s_sltcnd_excit = sltcnd_1(ref, matrix(exc), tParity)
    end function

    HElement_t(dp) function d_sltcnd_excit(ref, exc, tParity)
        integer, intent(in) :: ref(nel)
        type(DoubleExc_t), intent(in) :: exc
        logical, intent(in) :: tParity

        unused_var(ref)

        d_sltcnd_excit = sltcnd_2(matrix(exc), tParity)
    end function

end module new_sltcnd_mod
