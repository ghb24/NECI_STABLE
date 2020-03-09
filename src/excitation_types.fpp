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
    use constants, only: dp, n_int
    use bit_rep_data, only: nIfTot
    use SystemData, only: nEl
    use orb_idx_mod, only: SpinOrbIdx_t
    implicit none
    private
    public :: excitation_t, NoExc_t, SingleExc_t, DoubleExc_t, TripleExc_t, &
        FurtherExc_t, UNKNOWN, defined, last_tgt_unknown, set_last_tgt, &
        create_excitation, get_excitation, get_bit_excitation, excite


!> Arbitrary non occuring (?!) orbital index.
    integer, parameter :: UNKNOWN = -10**5

!>  @brief
!>      Abstract base class for excitations.
    type, abstract :: excitation_t
    end type

!>  @brief
!>      Represents a No-Op excitation.
    type, extends(excitation_t) :: NoExc_t
    end type

!>  @brief
!>      Represents the orbital indices of a single excitation.
!>      The array is sorted like:
!>      [src1, tgt2]
    type, extends(excitation_t) :: SingleExc_t
        integer :: val(2) = UNKNOWN
    end type

!>  @brief
!>      Represents the orbital indices of a double excitation.
!>      The array is sorted like:
!>      [[src1, src2],
!>      [tgt1, tgt2]]
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

!> Represents an excitation with so many differing indices,
!> that it is for sure zero.
    type, extends(excitation_t) :: FurtherExc_t
    end type

    #:for excitation_t in non_trivial_excitations
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
    interface ${excitation_t}$
        module procedure from_integer_${excitation_t}$
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
    #:for excitation_t in non_trivial_excitations
        module procedure defined_${excitation_t}$
    #:endfor
    end interface

!>  @brief
!>     Return true if all sources are known and all targets
!>     except the last are known.
!>
!>  @author Oskar Weser
!>
!>  @details
!>
!>  @param[in] exc, A non_trivial_excitation.
    interface last_tgt_unknown
    #:for excitation_t in non_trivial_excitations
        module procedure last_tgt_unknown_${excitation_t}$
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
    #:for excitation_t in non_trivial_excitations
        module procedure set_last_tgt_${excitation_t}$
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
    #:for excitation_t in ['NoExc_t', 'SingleExc_t', 'DoubleExc_t']
        module procedure excite_${excitation_t}$
    #:endfor
    end interface


contains

    #:for excitation_t in non_trivial_excitations
        elemental function defined_${excitation_t}$(exc) result(res)
            type(${excitation_t}$), intent(in) :: exc
            logical :: res

            res = all(exc%val /= UNKNOWN)
        end function
    #:endfor

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
        class(excitation_t), allocatable, intent(out) :: exc
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'create_excitation'
#endif

        select case (IC)
        case(0)
            allocate(NoExc_t :: exc)
        case(1)
            allocate(SingleExc_t :: exc)
        case(2)
            allocate(DoubleExc_t :: exc)
        case(3)
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
        class(excitation_t), allocatable, intent(out) :: exc
        logical, intent(out) :: tParity

        call create_excitation(exc, IC)

        ! The compiler has to statically know, what the type of exc is.
        select type (exc)
        type is (SingleExc_t)
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (DoubleExc_t)
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (TripleExc_t)
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
        class(excitation_t), allocatable, intent(out) :: exc
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

    #:for excitation_t in ['DoubleExc_t', 'TripleExc_t']
        pure subroutine set_last_tgt_${excitation_t}$(exc, tgt)
            type(${excitation_t}$), intent(inout) :: exc
            integer, intent(in) :: tgt
            exc%val(2, size(exc%val, 2)) = tgt
        end subroutine set_last_tgt_${excitation_t}$
    #:endfor


    #:for excitation_t in non_trivial_excitations
        pure function last_tgt_unknown_${excitation_t}$(exc) result(res)
            type(${excitation_t}$), intent(in) :: exc
            logical :: res
            integer :: flattened(size(exc%val))
            flattened(:) = pack(exc%val, .true.)
            res = (all(flattened(: size(flattened) - 1) /= UNKNOWN) &
                   .and.  flattened(size(flattened)) == UNKNOWN)
        end function last_tgt_unknown_${excitation_t}$
    #:endfor

    pure function excite_NoExc_t(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(NoExc_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        @:unused_var(exc)
        res = det_I
    end function

    DEBUG_IMPURE function excite_SingleExc_t(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(SingleExc_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        character(*), parameter :: this_routine = 'excite_SingleExc_t'

        associate(src => exc%val(1), tgt => exc%val(2))
            ASSERT(defined(exc))
            ASSERT(src /= tgt)
            ASSERT(any(src == det_I%idx))
            ASSERT(all(tgt /= det_I%idx))
            res%idx = insert_delete_sorted(det_I%idx, [src], [tgt])
        end associate
    end function


    DEBUG_IMPURE function excite_DoubleExc_t(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(DoubleExc_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        character(*), parameter :: this_routine = 'excite_DoubleExc_t'

        integer :: src(2), tgt(2), i


        src = exc%val(1, :)
        tgt = exc%val(2, :)
        ASSERT(defined(exc))
        do i = 1, 2
            ASSERT(all(src(i) /= tgt))
            ASSERT(any(src(i) == det_I%idx))
            ASSERT(all(tgt(i) /= det_I%idx))
        end do

        if (src(1) > src(2)) call swap(src(1), src(2))
        if (tgt(1) > tgt(2)) call swap(tgt(1), tgt(2))

        res%idx = insert_delete_sorted(det_I%idx, src, tgt)

        contains
            pure subroutine swap(a, b)
                integer, intent(inout) :: a, b
                integer :: tmp
                tmp = a
                a = b
                b = tmp
            end subroutine
    end function

    !> Merge C into A and remove values of B.
    !> Preconditions (not tested!):
    !>      1. B is a subset of A
    !>      2. A and C are disjoint
    !>      3. B and C are disjoint
    !>      4. B and C have the same size.
    !>      5. A, B, and C are sorted.
    pure function insert_delete_sorted(A, B, C) result(D)
        integer, intent(in) :: A(:), B(:), C(:)
        integer :: D(size(A))

        integer :: i, j, k, l

        i = 1
        j = 1
        k = 1
        l = 1
        do while(l <= size(D))
            ! Only indices from C have to be added to A
            if (i > size(A)) then
                D(l) = C(k)
                k = k + 1
                l = l + 1
            ! No more indices from B have to be deleted in A
            else if (j > size(B)) then
                if (A(i) < C(k)) then
                    D(l) = A(i)
                    i = i + 1
                    l = l + 1
                else if (A(i) > C(k)) then
                    D(l) = C(k)
                    k = k + 1
                    l = l + 1
                end if
            ! No more indices have to be added from C to A
            else if (k > size(C)) then
                if (A(i) /= B(j)) then
                    D(l) = A(i)
                    i = i + 1
                    l = l + 1
                else
                    i = i + 1
                    j = j + 1
                end if
            else if (A(i) < C(k)) then
                if (A(i) /= B(j)) then
                    D(l) = A(i)
                    i = i + 1
                    l = l + 1
                else
                    i = i + 1
                    j = j + 1
                end if
            else if (A(i) > C(k)) then
                if (A(i) /= B(j)) then
                    D(l) = C(k)
                    k = k + 1
                    l = l + 1
                else
                    D(l) = C(k)
                    i = i + 1
                    j = j + 1
                    k = k + 1
                    l = l + 1
                end if
            end if
        end do
    end function

end module
