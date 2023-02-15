#include "macros.h"
#:include "macros.fpph"
! while in principle further excitations are possible, in all current implementations
! of NECI, they are always zero
#:set excit_ranks = [0, 1, 2, 3]
#:set excitations = [f'Excite_{i}_t' for i in excit_ranks]
#:set trivial_excitations = excitations[:1]
#:set non_trivial_excitations = excitations[1:]
#:set classic_abinit_excitations = excitations[:3]

!>  A module for representing different excitations.
!>
!>  There is one abstract base class excitation_t that represents
!>  an arbitrary excitation.
!>  Possible excitations are the
!>  non trivial excitations Excite_{1,2,3}_t
!>  and the trivial excitation Excite_0_t.
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
        ilut_excite, excite, dyn_excite
    #:for excit in excitations
    public :: ${excit}$
    #:endfor

    !> Arbitrary non occuring (?!) orbital index.
    integer, parameter :: UNKNOWN = -20

    !>  Abstract base class for excitations.
    type, abstract :: Excitation_t
    end type

    #:for rank in excit_ranks
    type, extends(Excitation_t) :: Excite_${rank}$_t
        !! Represents the orbital indices of a ${rank}$-order excitation
        !! The array is sorted like:
        !! [srcs, tgts]
        integer :: val(2, ${rank}$) = UNKNOWN
    end type
    #:endfor

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

#ifdef IFORT_
    interface Excite_0_t
        module procedure Excite_0_t_ctor
    end interface
#endif


!>  Return true if all sources and targets are not UNKNOWN.
    interface defined
    #:for Excitation_t in excitations
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
        #:for Excitation_t in classic_abinit_excitations
            module procedure excite_${det_type}$_${Excitation_t}$
        #:endfor
    #:endfor
    end interface

    interface ilut_excite
    #:for det_type in ['Ilut_t']
        #:for Excitation_t in classic_abinit_excitations
        module procedure excite_${det_type}$_${Excitation_t}$
        #:endfor
    #:endfor
    end interface

    interface get_excitation
        module procedure get_excitation_old, get_excitation_new
    end interface

contains

! workaround for pesky intel compiler errors
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
        #:for Excitation_t in excitations
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
        case (0)
            allocate(Excite_0_t :: exc)
        case (1)
            allocate(Excite_1_t :: exc)
        case (2)
            allocate(Excite_2_t :: exc)
        case (3)
            allocate(Excite_3_t :: exc)
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

    pure function excite_nI_Excite_1_t(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(Excite_1_t), intent(in) :: exc
        integer :: res(size(det_I))
        character(*), parameter :: this_routine = 'excite_nI_Excite_1_t'

        @:pure_ASSERT(defined(exc))
        associate(src => exc%val(1, 1), tgt => exc%val(2, 1))
            @:pure_ASSERT(src /= tgt)
            @:pure_ASSERT(disjoint([tgt], det_I))
            @:pure_ASSERT(subset([src], det_I))
            res = special_union_complement(det_I, [tgt], [src])
        end associate
    end function


    pure function excite_nI_Excite_2_t(det_I, exc) result(res)
        integer, intent(in) :: det_I(:)
        type(Excite_2_t), intent(in) :: exc
        integer :: res(size(det_I))
        character(*), parameter :: this_routine = 'excite_nI_Excite_2_t'
        integer :: src(2), tgt(2)

        @:pure_ASSERT(defined(exc))
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

#:for Excitation_t in classic_abinit_excitations
    pure function excite_SpinOrbIdx_t_${Excitation_t}$(det_I, exc) result(res)
        type(SpinOrbIdx_t), intent(in) :: det_I
        type(${Excitation_t}$), intent(in) :: exc
        type(SpinOrbIdx_t) :: res
        res%idx = excite(det_I%idx, exc)
    end function
#:endfor

    pure function excite_Ilut_t_Excite_0_t(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(Excite_0_t), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        @:unused_var(exc)
        res = ilut_I
    end function

    pure function excite_Ilut_t_Excite_1_t(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(Excite_1_t), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        character(*), parameter :: this_routine = 'excite_Ilut_t_Excite_1_t'

        associate(src => exc%val(1, 1), tgt => exc%val(2, 1))
            @:pure_ASSERT(defined(exc))
            @:pure_ASSERT(src /= tgt)
            res = ilut_I
            clr_orb(res, src)
            set_orb(res, tgt)
        end associate
    end function

    pure function excite_Ilut_t_Excite_2_t(ilut_I, exc) result(res)
        integer(n_int), intent(in) :: ilut_I(:)
        type(Excite_2_t), intent(in) :: exc
        integer(n_int) :: res(0:size(ilut_I) - 1)
        character(*), parameter :: this_routine = 'excite_Ilut_t_Excite_2_t'

        integer :: src(2), tgt(2), i

        src = exc%val(1, :)
        tgt = exc%val(2, :)
        @:pure_ASSERT(defined(exc))
        do i = 1, 2
            @:pure_ASSERT(all(src(i) /= tgt))
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
        type(SpinOrbIdx_t) :: res

        select type (exc)
        type is (Excite_0_t)
            res = excite(det_I, exc)
        type is (Excite_1_t)
            res = excite(det_I, exc)
        type is (Excite_2_t)
            res = excite(det_I, exc)
        end select
    end function
end module
