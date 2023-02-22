#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

#:set max_excit_rank = 3
    ! all excitations with rank higher than max_excit_rank are definitely zero
#:set excit_ranks = list(range(max_excit_rank + 1))
    ! note that this excludes further excitations, which must be handled manually
#:set excitations = [f'Excite_{i}_t' for i in excit_ranks + ['Further']]
    ! Excite_Further_t is for all ranks > max_excit_rank
#:set defined_excitations = excitations[:-1]
#:set trivial_excitations = [excitations[0], excitations[-1]]
#:set non_trivial_excitations = excitations[1:-1]

!>  A module for representing different excitations.
!>
!>  There is one abstract base class excitation_t that represents
!>  an arbitrary excitation.
!>  Possible excitations are the
!>  non trivial excitations Excite_{1,2,3}_t
!>  and the trivial excitations Excite_0_t and Excite_Further_t.
!>
!>  The non trivial excitations can "know" only some indices
!>  and leave the rest UNKNOWN, which is an arbitrary integer constant.
!>
!>  The procedures create_excitation, get_excitation, and get_bit_excitation
!>  can be used, to create excitations from nIs, or iluts at runtime.
module excitation_types
    use orb_idx_mod, only: size
    use constants, only: dp, n_int, bits_n_int, maxExcit
    use bit_rep_data, only: nIfTot
    use util_mod, only: stop_all
    use SystemData, only: nEl
    use orb_idx_mod, only: SpinOrbIdx_t
    use sets_mod, only: disjoint, subset, is_sorted, special_union_complement
    use DetBitOps, only: GetBitExcitation
    use neci_intfce, only: GetExcitation
    implicit none
    private
    public :: Excitation_t, UNKNOWN, defined, dyn_defined, get_last_tgt, set_last_tgt, &
        create_excitation, get_excitation, get_bit_excitation, &
        ilut_excite, excite, dyn_excite, dyn_nI_excite
    #:for excit in excitations
    public :: ${excit}$
    #:endfor

    !> Arbitrary non occuring (?!) orbital index.
    integer, parameter :: UNKNOWN = -20

    !>  Abstract base class for excitations.
    type, abstract :: Excitation_t
    end type

    #:for rank, excite_t in zip(excit_ranks, defined_excitations)
    type, extends(Excitation_t) :: ${excite_t}$
        !! Represents the orbital indices of a ${rank}$-order excitation
        !! The array is sorted like:
        !! [srcs, tgts]
        integer :: val(2, ${rank}$) = UNKNOWN
    end type
    #:endfor

    type, extends(Excitation_t) :: Excite_Further_t
        !! Represents an excitation with so many different indices, it has to be zero
    end type

    #:for Excitation_t in non_trivial_excitations
!>  Additional constructors for the excitation types from integers instead
!>  of an integer array.
!>
!>  The non trivial excitations
!>  are initialized by passing the respective integer arrays into
!>  the type.
!>  Alternatively one can use integer arguments to initialize.
!>  Omitted indices are set to UNKNOWN.
!>
!>  \code{.unparsed}
!>  Excite_t([1, 2]) == Excite_t(src=1, tgt=2)
!>  ! If the target should be UNKNOWN, just omit it
!>  Excite_t(src=1)
!>  \endcode
!>
!>  The signature is `(src_1, tgt_1, src_2, tgt_2, ...)`.
!>  depending on the actual type.
    interface ${Excitation_t}$
        module procedure from_integer_${Excitation_t}$
    end interface
    #:endfor

! workaround for pesky intel compiler errors
! this should be removed if IFORT_ works properly at some pointer in the future
#ifdef IFORT_
    interface Excite_0_t
        module procedure Excite_0_t_ctor
    end interface
#endif


!>  Return true if all sources and targets are not UNKNOWN.
    interface defined
    #:for Excitation_t in defined_excitations
        module procedure defined_${Excitation_t}$
    #:endfor
    end interface

!>  Get the last target of a non trivial excitation.
    interface get_last_tgt
    #:for Excitation_t in non_trivial_excitations
        module procedure get_last_tgt_${Excitation_t}$
    #:endfor
    end interface

!>  Set the last target of a non trivial excitation.
    interface set_last_tgt
    #:for Excitation_t in non_trivial_excitations
        module procedure set_last_tgt_${Excitation_t}$
    #:endfor
    end interface

!>  Perform the excitation on a given determinant.
!>
!>  It is assumed that the excitations are non trivial.
!>  I.e. for single excitations source /= target
!>  and for double excitations the set of sources and targets has
!>  to be disjoint.
    interface excite
    #:for det_type in ['nI', 'SpinOrbIdx_t']
        #:for Excitation_t in defined_excitations
            module procedure excite_${det_type}$_${Excitation_t}$
        #:endfor
    #:endfor
    end interface

    interface ilut_excite
    #:for det_type in ['Ilut_t']
        #:for Excitation_t in defined_excitations
        module procedure excite_${det_type}$_${Excitation_t}$
        #:endfor
    #:endfor
    end interface

    interface get_excitation
        module procedure get_excitation_old, get_excitation_new
    end interface

contains

! workaround for pesky intel compiler errors
! this should be removed if IFORT_ works properly at some pointer in the future
#ifdef IFORT_
    type(Excite_0_t) function Excite_0_t_ctor() result(this)
        integer :: tmpval(2, 0)
        tmpval = UNKNOWN
        this%val = tmpval
    end function
#endif

    #:for Excitation_t in non_trivial_excitations
    elemental function defined_${Excitation_t}$ (exc) result(res)
        type(${Excitation_t}$), intent(in) :: exc
        logical :: res

        res = all(exc%val /= UNKNOWN)
    end function
    #:endfor
    elemental function defined_Excite_0_t(exc) result(res)
        type(Excite_0_t), intent(in) :: exc
        logical :: res
        @:unused_var(exc)
        res = .true.
    end function

    elemental function dyn_defined(exc) result(res)
        class(Excitation_t), intent(in) :: exc
        logical :: res

        select type (exc)
        #:for Excitation_t in defined_excitations
        type is (${Excitation_t}$)
            res = defined(exc)
        #:endfor
        end select
    end function

    pure function from_integer_Excite_1_t(src, tgt) result(res)
        integer, intent(in), optional :: src, tgt
        type(Excite_1_t) :: res

        ! The values default to UNKNOWN automatically.
        if (present(src)) res%val(1, 1) = src
        if (present(tgt)) res%val(2, 1) = tgt
    end function

    pure function from_integer_Excite_2_t(src1, tgt1, src2, tgt2) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2
        type(Excite_2_t) :: res

        ! The values default to UNKNOWN automatically.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
    end function

    pure function from_integer_Excite_3_t(src1, tgt1, src2, tgt2, src3, tgt3) result(res)
        integer, intent(in), optional :: src1, tgt1, src2, tgt2, src3, tgt3
        type(Excite_3_t) :: res

        ! The values default to UNKNOWN automatically.
        if (present(src1)) res%val(1, 1) = src1
        if (present(tgt1)) res%val(2, 1) = tgt1
        if (present(src2)) res%val(1, 2) = src2
        if (present(tgt2)) res%val(2, 2) = tgt2
        if (present(src2)) res%val(1, 3) = src3
        if (present(tgt2)) res%val(2, 3) = tgt3
    end function

!>  Create an excitation from an excitation matrix and excitation level IC
    subroutine create_excitation(exc, ic, ex)
        !>  The excitation level. (1=Excite_1_t, 2=Excite_2_t, ...)
        integer, intent(in) :: IC
        !>  An excitation matrix as in the %val component of
        !>      the excitation types.
        integer, intent(in), optional :: ex(2, ic)
        !>  An excitation of type excitation_t.
        !>      By using select type(exc) one can select the actual type at runtime
        !>      **and** statically dispatch as much as possible at runtime.
        class(Excitation_t), allocatable, intent(out) :: exc
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'create_excitation'
#endif

        select case (IC)
        #:for rank, excite_t in zip(excit_ranks, defined_excitations)
        case (${rank}$)
            allocate(${excite_t}$ :: exc)
        #:endfor
        case (${max_excit_rank + 1}$ :)
            allocate(Excite_Further_t :: exc)
        case (:-1)
#ifdef DEBUG_
            call stop_all(this_routine, 'invalid IC < 0 passed.')
#endif
        end select

        if (present(ex)) then
            select type (exc)
            type is (Excite_1_t)
                exc%val = ex
            type is (Excite_2_t)
                exc%val = ex
            type is (Excite_3_t)
                exc%val = ex
            end select
        end if
    end subroutine

!>  Create an excitation from nI to nJ where the excitation level
!>  is already known.
    subroutine get_excitation_new(nI, nJ, IC, exc, tParity)
        !> Two Slater determinants in nI format.
        integer, intent(in) :: nI(nEl), nJ(nEl)
        !>  The excitation level. (1=Excite_1_t, 2=Excite_2_t, ...)
        integer, intent(in) :: IC
        !>  An excitation of type excitation_t.
        !>      By using select type(exc) one can select the actual type at runtime
        !>      **and** statically dispatch as much as possible at compile time.
        class(Excitation_t), allocatable, intent(out) :: exc
        !>  The parity of the excitation.
        logical, intent(out) :: tParity

        call create_excitation(exc, IC)

        ! The compiler has to statically know, what the type of exc is.
        select type (exc)
        type is (Excite_1_t)
            exc%val(1, 1) = 1
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (Excite_2_t)
            exc%val(1, 1) = 2
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        type is (Excite_3_t)
            exc%val(1, 1) = 3
            call GetExcitation(nI, nJ, nel, exc%val, tParity)
        end select
    end subroutine get_excitation_new

    subroutine get_excitation_old(nI, nJ, ic, exc, tParity)
        integer, intent(in) :: nI(nEl), nJ(nEl), IC
        integer, intent(out) :: exc(2, maxExcit)
        logical, intent(out) :: tParity
        character(*), parameter :: this_routine = 'get_excitation_old'
        @:ASSERT(any(ic == [1, 2, 3]))
        exc(1, 1) = ic
        call GetExcitation(nI, nJ, nel, exc, tParity)
    end subroutine


!>  Create an excitation from ilutI to ilutJ where the excitation level
!>  is already known.
    subroutine get_bit_excitation(ilutI, ilutJ, IC, exc, tParity)
        !> Two Slater determinants in bitmask format.
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        !>  The excitation level. (1=Excite_1_t, 2=Excite_2_t, ...)
        integer, intent(in) :: IC
        !>  The parity of the excitation.
        class(Excitation_t), allocatable, intent(out) :: exc
        logical, intent(out) :: tParity

        call create_excitation(exc, IC)

        ! The compiler has to statically know, what the type of exc is.
        select type (exc)
        type is (Excite_1_t)
            exc%val(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, exc%val, tParity)
        type is (Excite_2_t)
            exc%val(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, exc%val, tParity)
        type is (Excite_3_t)
            exc%val(1, 1) = IC
            call GetBitExcitation(iLutI, iLutJ, exc%val, tParity)
        end select
    end subroutine get_bit_excitation

    #:for Excitation_t in non_trivial_excitations
    pure subroutine set_last_tgt_${Excitation_t}$ (exc, tgt)
        type(${Excitation_t}$), intent(inout) :: exc
        integer, intent(in) :: tgt
        exc%val(2, size(exc%val, 2)) = tgt
    end subroutine set_last_tgt_${Excitation_t}$
    #:endfor

    #:for Excitation_t in non_trivial_excitations
    pure function get_last_tgt_${Excitation_t}$ (exc) result(res)
        type(${Excitation_t}$), intent(in) :: exc
        integer :: res

        res = exc%val(2, size(exc%val, 2))
    end function get_last_tgt_${Excitation_t}$
    #:endfor

    pure function excite_nI_Excite_0_t(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(Excite_0_t), intent(in) :: exc
        integer :: res(size(det_I))
        @:unused_var(exc)
        res = det_I
    end function

    #:for rank, excite_t in zip(excit_ranks[1:], non_trivial_excitations)
    pure function excite_nI_${excite_t}$(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(${excite_t}$), intent(in) :: exc
        integer :: res(size(det_I))
        character(*), parameter :: this_routine = 'excite_nI_${excite_t}$'

        integer :: src(${rank}$), tgt(${rank}$)

        @:pure_ASSERT(defined(exc))
        src = exc%val(1, :)
        tgt = exc%val(2, :)
        @:sort(integer, src)
        @:sort(integer, tgt)
        @:pure_ASSERT(is_sorted(src))
        @:pure_ASSERT(is_sorted(tgt))
        @:pure_ASSERT(disjoint(src, tgt))
        @:pure_ASSERT(disjoint(tgt, det_I))
        @:pure_ASSERT(subset(src, det_I))

        res = special_union_complement(det_I, tgt, src)
    end function
    #:endfor

    #:for Excitation_t in defined_excitations
    pure function excite_SpinOrbIdx_t_${Excitation_t}$(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(${Excitation_t}$), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        res%idx = excite(det_I%idx, exc)
    end function
    #:endfor

    #:for rank, excite_t in zip(excit_ranks, defined_excitations)
    pure function excite_Ilut_t_${excite_t}$(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(${excite_t}$), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        character(*), parameter :: this_routine = 'excite_Ilut_t_${excite_t}$'

        integer :: src(${rank}$), tgt(${rank}$), i

        src = exc%val(1, :)
        tgt = exc%val(2, :)
        @:pure_ASSERT(defined(exc))
        res = ilut_I
        do i = 1, ${rank}$
            @:pure_ASSERT(all(src(i) /= tgt))
            clr_orb(res, src(i))
            set_orb(res, tgt(i))
        end do
    end function
    #:endfor

    pure function dyn_excite(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        class(Excitation_t), intent(in) :: exc
        type(SpinOrbIdx_t) :: res

        select type (exc)
        #:for Excitation_t in non_trivial_excitations
        type is (${Excitation_t}$)
            res = excite(det_I, exc)
        #:endfor
        end select
    end function

    pure function dyn_nI_excite(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        class(Excitation_t), intent(in) :: exc
        integer :: res(size(det_I))

        select type (exc)
        #:for Excitation_t in defined_excitations
        type is (${Excitation_t}$)
            res = excite(det_I, exc)
        #:endfor
        class default
            call stop_all("dyn_nI_excite", "Excitation type invalid.")
        end select
    end function
end module
