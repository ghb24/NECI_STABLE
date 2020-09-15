#include "macros.h"
#:include "macros.fpph"
#:set excitations = ['NoExc_t', 'FurtherExc_t', 'SingleExc_t', 'DoubleExc_t', 'TripleExc_t']
#:set trivial_excitations = excitations[:2]
#:set non_trivial_excitations = excitations[2:]

!>  @brief
!>      A module for representing different excitations.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  There is one abstract base class excitation_t that represents
!>  an arbitrary excitation.
!>  Possible excitations are the
!>  non trivial excitations SingleExc_t, DoubleExc_t, and TripleExc_t
!>  and the trivial excitations NoExc_t and FurtherExc_t.
!>
!>  The non trivial excitations can "know" only some indices
!>  and leave the rest UNKNOWN, which is an arbitrary integer constant.
!>
!>  The procedures create_excitation, get_excitation, and get_bit_excitation
!>  can be used, to create excitations from nIs, or iluts at runtime.
module excitation_types
    use constants, only: dp, n_int, bits_n_int
    use bit_rep_data, only: nIfTot
    use util_mod, only: stop_all
    use SystemData, only: nEl
    use orb_idx_mod, only: SpinOrbIdx_t
    use sets_mod, only: disjoint, subset, is_sorted, special_union_complement
    implicit none
    private
    public :: Excitation_t, NoExc_t, SingleExc_t, DoubleExc_t, TripleExc_t, &
        FurtherExc_t, UNKNOWN, defined, dyn_defined, get_last_tgt, set_last_tgt, &
        create_excitation, get_excitation, get_bit_excitation, &
        ilut_excite, excite, dyn_excite


!> Arbitrary non occuring (?!) orbital index.
    integer, parameter :: UNKNOWN = -20

!>  @brief
!>      Abstract base class for excitations.
    type, abstract :: Excitation_t
    end type

!>  @brief
!>      Represents a No-Op excitation.
    type, extends(Excitation_t) :: NoExc_t
    end type

!>  @brief
!>      Represents the orbital indices of a single excitation.
!>      The array is sorted like:
!>      [src1, tgt2]
    type, extends(Excitation_t) :: SingleExc_t
        integer :: val(2) = UNKNOWN
    end type

!>  @brief
!>      Represents the orbital indices of a double excitation.
!>      The array is sorted like:
!>      [[src1, src2],
!>      [tgt1, tgt2]]
    type, extends(Excitation_t) :: DoubleExc_t
        integer :: val(2, 2) = UNKNOWN
    end type

!> Represents the orbital indices of a triple excitation.
!> The array is sorted like:
!> [[src1, src2, src3],
!>  [tgt1, tgt2, tgt3]]
    type, extends(Excitation_t) :: TripleExc_t
        integer :: val(2, 3) = UNKNOWN
    end type

!> Represents an excitation with so many differing indices,
!> that it is for sure zero.
    type, extends(Excitation_t) :: FurtherExc_t
    end type

    #:for Excitation_t in non_trivial_excitations
!>  @brief
!>      Additional constructors for the excitation types from integers instead
!>      of an integer array.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The non trivial excitations SingleExc_t, DoubleExc_t, and TripleExc_t
!>  are initialized by passing the respective integer arrays into
!>  the type.
!>  Alternatively one can use integer arguments to initialize.
!>  Omitted indices are set to UNKNOWN.
!>
!>  \code{.unparsed}
!>  SingleExc_t([1, 2]) == SingleExc_t(src=1, tgt=2)
!>  ! If the target should be UNKNOWN, just omit it
!>  SingleExc_t(src=1)
!>  \endcode
!>
!>  The signature is (src_1, tgt_1, src_2, tgt_2, ...).
!>  depending on the actual type.
!>
!>  @param[in] src_i
!>  @param[in] tgt_i
    interface ${Excitation_t}$
        module procedure from_integer_${Excitation_t}$
    end interface
    #:endfor

!>  @brief
!>     Return true if all sources and targets are not UNKNOWN.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @param[in] exc, A non_trivial_excitation.
    interface defined
    #:for Excitation_t in non_trivial_excitations + ['NoExc_t']
        module procedure defined_${Excitation_t}$
    #:endfor
    end interface

!>  @brief
!>     Get the last target of a non trivial excitation.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @param[in] exc, A non_trivial_excitation.
    interface get_last_tgt
    #:for Excitation_t in non_trivial_excitations
        module procedure get_last_tgt_${Excitation_t}$
    #:endfor
    end interface

!>  @brief
!>     Set the last target of a non trivial excitation.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @param[in] exc, A non_trivial_excitation.
!>  @param[in] tgt, Index of target.
    interface set_last_tgt
    #:for Excitation_t in non_trivial_excitations
        module procedure set_last_tgt_${Excitation_t}$
    #:endfor
    end interface

!>  @brief
!>     Perform the excitation on a given determinant.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is assumed, that the excitations are non trivial.
!>  I.e. for single excitations source /= target
!>  and for double excitations the set of sources and targets has
!>  to be disjoint.
!>
!>  @param[in] det_I, A Slater determinant of SpinOrbIdx_t.
!>  @param[in] exc, NoExc_t, SingleExc_t, or DoubleExc_t.
    interface excite
    #:for det_type in ['nI', 'SpinOrbIdx_t']
        #:for Excitation_t in ['NoExc_t', 'SingleExc_t', 'DoubleExc_t']
            module procedure excite_${det_type}$_${Excitation_t}$
        #:endfor
    #:endfor
    end interface

    interface ilut_excite
    #:for det_type in ['Ilut_t']
        #:for Excitation_t in ['NoExc_t', 'SingleExc_t', 'DoubleExc_t']
        module procedure excite_${det_type}$_${Excitation_t}$
        #:endfor
    #:endfor
    end interface

contains

    #:for Excitation_t in non_trivial_excitations
    elemental function defined_${Excitation_t}$ (exc) result(res)
        type(${Excitation_t}$), intent(in) :: exc
        logical :: res

        res = all(exc%val /= UNKNOWN)
    end function
    #:endfor
    elemental function defined_NoExc_t(exc) result(res)
        type(NoExc_t), intent(in) :: exc
        logical :: res
        @:unused_var(exc)
        res = .true.
    end function

    elemental function dyn_defined(exc) result(res)
        class(Excitation_t), intent(in) :: exc
        logical :: res

        select type (exc)
        type is (NoExc_t)
            res = defined(exc)
        type is (SingleExc_t)
            res = defined(exc)
        type is (DoubleExc_t)
            res = defined(exc)
        type is (TripleExc_t)
            res = defined(exc)
        end select
    end function

    pure function from_integer_SingleExc_t(src, tgt) result(res)
        integer, intent(in), optional :: src, tgt
        type(SingleExc_t) :: res

        ! The values default to UNKNOWN automatically.
        if (present(src)) res%val(1) = src
        if (present(tgt)) res%val(2) = tgt
    end function

    pure function from_integer_DoubleExc_t(src1, tgt1, src2, tgt2) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2
        type(DoubleExc_t) :: res

        ! The values default to UNKNOWN automatically.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
    end function

    pure function from_integer_TripleExc_t(src1, tgt1, src2, tgt2, src3, tgt3) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2, src3, tgt3
        type(TripleExc_t) :: res

        ! The values default to UNKNOWN automatically.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
        if (present(src2)) res%val(1, 3) = src3
        if (present(tgt2)) res%val(2, 3) = tgt3
    end function

!>  @brief
!>      Create an excitation from an excitation matrix and excitation level IC
!>
!>  @param[out] exc, An excitation of type excitation_t.
!>      By using select type(exc) one can select the actual type at runtime
!>      **and** statically dispatch as much as possible at runtime.
!>  @param[in] ic, The excitation level. (1=SingleExc_t, 2=DoubleExc_t, ...)
!>  @param[in] ex, An excitation matrix as in the %val component of
!>      the excitation types.
    subroutine create_excitation(exc, ic, ex)
        integer, intent(in) :: IC
        integer, intent(in), optional :: ex(2, ic)
        class(Excitation_t), allocatable, intent(out) :: exc
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'create_excitation'
#endif

        select case (IC)
        case (0)
            allocate(NoExc_t :: exc)
        case (1)
            allocate(SingleExc_t :: exc)
        case (2)
            allocate(DoubleExc_t :: exc)
        case (3)
            allocate(TripleExc_t :: exc)
        case (4:)
            allocate(FurtherExc_t :: exc)
        case default
#ifdef DEBUG_
            call stop_all(this_routine, 'invalid IC < 0 passed.')
#endif
        end select

        if (present(ex)) then
            select type (exc)
            type is (SingleExc_t)
                exc%val = ex(:, 1)
            type is (DoubleExc_t)
                exc%val = ex
            type is (TripleExc_t)
                exc%val = ex
            end select
        end if
    end subroutine

!>  @brief
!>      Create an excitation from nI to nJ where the excitation level
!>      is already known.
!>
!>  @param[in] nI, An array of occupied orbital indices.
!>  @param[in] nJ, An array of occupied orbital indices.
!>  @param[in] ic, The excitation level. (1=SingleExc_t, 2=DoubleExc_t, ...)
!>  @param[out] exc, An excitation of type excitation_t.
!>      By using select type(exc) one can select the actual type at runtime
!>      **and** statically dispatch as much as possible at compile time.
!>  @param[out] tParity, The parity of the excitation.
    subroutine get_excitation(nI, nJ, IC, exc, tParity)
        integer, intent(in) :: nI(nEl), nJ(nEl), IC
        class(Excitation_t), allocatable, intent(out) :: exc
        logical, intent(out) :: tParity

        call create_excitation(exc, IC)

        ! The compiler has to statically know, what the type of exc is.
        select type (exc)
        type is (SingleExc_t)
            exc%val(1) = 1
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (DoubleExc_t)
            exc%val(1, 1) = 2
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (TripleExc_t)
            exc%val(1, 1) = 3
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        end select
    end subroutine get_excitation

!>  @brief
!>      Create an excitation from ilutI to ilutJ where the excitation level
!>      is already known.
!>
!>  @param[in] ilutI, A bitmask encoding occupation of spin orbitals.
!>  @param[in] ilutJ, A bitmask encoding occupation of spin orbitals.
!>  @param[in] ic, The excitation level. (1=SingleExc_t, 2=DoubleExc_t, ...)
!>  @param[out] exc, An excitation of type excitation_t.
!>      By using select type(exc) one can select the actual type at runtime
!>      **and** statically dispatch as much as possible at runtime.
!>  @param[out] tParity, The parity of the excitation.
    subroutine get_bit_excitation(ilutI, ilutJ, IC, exc, tParity)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: IC
        class(Excitation_t), allocatable, intent(out) :: exc
        logical, intent(out) :: tParity

        call create_excitation(exc, IC)

        ! The compiler has to statically know, what the type of exc is.
        select type (exc)
        type is (SingleExc_t)
            exc%val(1) = IC
            call GetBitExcitation(iLutI, iLutJ, exc%val, tParity)
        type is (DoubleExc_t)
            exc%val(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, exc%val, tParity)
        type is (TripleExc_t)
            exc%val(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, exc%val, tParity)
        end select
    end subroutine get_bit_excitation

    pure subroutine set_last_tgt_SingleExc_t(exc, tgt)
        type(SingleExc_t), intent(inout) :: exc
        integer, intent(in) :: tgt
        exc%val(2) = tgt
    end subroutine set_last_tgt_SingleExc_t

    #:for Excitation_t in ['DoubleExc_t', 'TripleExc_t']
    pure subroutine set_last_tgt_${Excitation_t}$ (exc, tgt)
        type(${Excitation_t}$), intent(inout) :: exc
        integer, intent(in) :: tgt
        exc%val(2, size(exc%val, 2)) = tgt
    end subroutine set_last_tgt_${Excitation_t}$
    #:endfor

    pure function get_last_tgt_SingleExc_t(exc) result(res)
        type(SingleExc_t), intent(in) :: exc
        integer :: res
        res = exc%val(2)
    end function get_last_tgt_SingleExc_t

    #:for Excitation_t in ['DoubleExc_t', 'TripleExc_t']
    pure function get_last_tgt_${Excitation_t}$ (exc) result(res)
        type(${Excitation_t}$), intent(in) :: exc
        integer :: res

        res = exc%val(2, size(exc%val, 2))
    end function get_last_tgt_${Excitation_t}$
    #:endfor

    pure function excite_nI_NoExc_t(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(NoExc_t), intent(in) :: exc
        integer :: res(size(det_I))
        @:unused_var(exc)
        res = det_I
    end function

    pure function excite_nI_SingleExc_t(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(SingleExc_t), intent(in) :: exc
        integer :: res(size(det_I))
        character(*), parameter :: this_routine = 'excite_SingleExc_t'

        @:pure_ASSERT(defined(exc), exc%val)
        associate(src => exc%val(1), tgt => exc%val(2))
            @:pure_ASSERT(src /= tgt)
            @:pure_ASSERT(disjoint([tgt], det_I))
            @:pure_ASSERT(subset([src], det_I))
            res = special_union_complement(det_I, [tgt], [src])
        end associate
    end function


    pure function excite_nI_DoubleExc_t(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(DoubleExc_t), intent(in) :: exc
        integer :: res(size(det_I))
        character(*), parameter :: this_routine = 'excite_DoubleExc_t'

        integer :: src(2), tgt(2), i

        @:pure_ASSERT(defined(exc), exc%val)
        src = exc%val(1, :)
        tgt = exc%val(2, :)
        if (src(1) > src(2)) call swap(src(1), src(2))
        if (tgt(1) > tgt(2)) call swap(tgt(1), tgt(2))
        @:pure_ASSERT(is_sorted(src))
        @:pure_ASSERT(is_sorted(tgt))
        @:pure_ASSERT(disjoint(src, tgt))
        @:pure_ASSERT(disjoint(tgt, det_I))
        @:pure_ASSERT(subset(src, det_I))

        res = special_union_complement(det_I, tgt, src)

    contains
        pure subroutine swap(a, b)
            integer, intent(inout) :: a, b
            integer :: tmp
            tmp = a
            a = b
            b = tmp
        end subroutine
    end function


    pure function excite_SpinOrbIdx_t_NoExc_t(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(NoExc_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        character(*), parameter :: this_routine = 'excite_NoExc_t'
        res%idx = excite(det_I%idx, exc)
    end function

    pure function excite_SpinOrbIdx_t_SingleExc_t(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SingleExc_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        character(*), parameter :: this_routine = 'excite_SingleExc_t'
        res%idx = excite(det_I%idx, exc)
    end function


    pure function excite_SpinOrbIdx_t_DoubleExc_t(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(DoubleExc_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        character(*), parameter :: this_routine = 'excite_DoubleExc_t'
        res%idx = excite(det_I%idx, exc)
    end function


    pure function excite_Ilut_t_NoExc_t(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(NoExc_t), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        @:unused_var(exc)
        res = ilut_I
    end function

    pure function excite_Ilut_t_SingleExc_t(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(SingleExc_t), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        character(*), parameter :: this_routine = 'excite_SingleExc_t'

        associate(src => exc%val(1), tgt => exc%val(2))
            @:pure_ASSERT(defined(exc), exc%val)
            @:pure_ASSERT(src /= tgt, src, tgt)
            res = ilut_I
            clr_orb(res, src)
            set_orb(res, tgt)
        end associate
    end function

    pure function excite_Ilut_t_DoubleExc_t(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(DoubleExc_t), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        character(*), parameter :: this_routine = 'excite_DoubleExc_t'

        integer :: src(2), tgt(2), i

        src = exc%val(1, :)
        tgt = exc%val(2, :)
        @:pure_ASSERT(defined(exc), exc%val)
        do i = 1, 2
            @:pure_ASSERT(all(src(i) /= tgt), src(i), tgt)
        end do
        res = ilut_I
        clr_orb(res, src(1))
        clr_orb(res, src(2))
        set_orb(res, tgt(1))
        set_orb(res, tgt(2))
    end function

    pure function dyn_excite(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        class(Excitation_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'dyn_excite'
        type(SpinOrbIdx_t) :: res

        select type (exc)
        type is (NoExc_t)
            res = excite(det_I, exc)
        type is (SingleExc_t)
            res = excite(det_I, exc)
        type is (DoubleExc_t)
            res = excite(det_I, exc)
        end select
    end function
end module
