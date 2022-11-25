#include "macros.h"

module excitation_generators
    use constants, only: dp, n_int, maxExcit
    use util_mod, only: operator(.isclose.)
    use SystemData, only: nel, nBasis
    use bit_rep_data, only: NIfTot
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use dSFMT_interface, only: genrand_real2_dSFMT
    use procedure_pointers, only: generate_excitation_t
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite
    use SymExcitDataMod, only: ScratchSize
    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet
    use sets_mod, only: operator(.complement.)
    use symexcit3, only: gen_excits
    use sort_mod, only: sort
    use procedure_pointers, only: generate_excitation, gen_all_excits
    implicit none
    private
    public :: ExcitationGenerator_t, &
        SingleExcitationGenerator_t, DoubleExcitationGenerator_t, TripleExcitationGenerator_t, &
        gen_exc_sd, get_pgen_sd, gen_all_excits_sd, ClassicAbInitExcitationGenerator_t

    type, abstract :: ExcitationGenerator_t
    contains
        procedure(BoundGenExc_t), public, deferred :: gen_exc
        procedure(BoundGetPgen_t), public, deferred :: get_pgen
        procedure(BoundGenAllExcits_t), public, deferred :: gen_all_excits
        procedure(BoundFinalize_t), public, deferred :: finalize
    end type

    type, abstract, extends(ExcitationGenerator_t) :: SingleExcitationGenerator_t
    contains
        procedure, public :: gen_all_excits => FCI_singles_gen_all_excits
    end type

    type, abstract, extends(ExcitationGenerator_t) :: DoubleExcitationGenerator_t
    contains
        procedure, public :: gen_all_excits => FCI_doubles_gen_all_excits
    end type

    type, abstract, extends(ExcitationGenerator_t) :: TripleExcitationGenerator_t
    end type

    type, abstract, extends(ExcitationGenerator_t) :: ClassicAbInitExcitationGenerator_t
        !! this abstract excitation generator covers all ab initio Hamiltonians
        !! in the typical sense (i.e. up to double excitations)
        private
        class(DoubleExcitationGenerator_t), public, allocatable :: doubles_generator
        ! NOTE: Change into class(SingleExcitationGenerator_t), allocatable
        !   if you want to change singles_generators at runtime.
        class(SingleExcitationGenerator_t), public, allocatable :: singles_generator
    contains
        private
        procedure, public :: gen_exc => abinit_PCHB_gen_exc
        procedure, public :: get_pgen => abinit_PCHB_get_pgen
        procedure, public :: gen_all_excits => abinit_PCHB_gen_all_excits
        procedure, public :: finalize => abinit_PCHB_finalize
    end type

    abstract interface
        subroutine BoundGenExc_t(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                         ex, tParity, pGen, hel, store, part_type)
            import :: ExcitationGenerator_t, n_int, dp, excit_gen_store_type, nEl, NifTot, maxExcit
            implicit none
            class(ExcitationGenerator_t), intent(inout) :: this
            integer, intent(in) :: nI(nel), exFlag
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
            integer(n_int), intent(out) :: ilutJ(0:NifTot)
            real(dp), intent(out) :: pGen
            logical, intent(out) :: tParity
            HElement_t(dp), intent(out) :: hel
            type(excit_gen_store_type), intent(inout), target :: store
            integer, intent(in), optional :: part_type
        end subroutine BoundGenExc_t

        subroutine BoundGenAllExcits_t(this, nI, n_excits, det_list)
            import :: ExcitationGenerator_t, n_int, nEl
            class(ExcitationGenerator_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
            integer, intent(out) :: n_excits
            integer(n_int), allocatable, intent(out) :: det_list(:,:)
        end subroutine BoundGenAllExcits_t

        subroutine BoundFinalize_t(this)
            import :: ExcitationGenerator_t
            class(ExcitationGenerator_t), intent(inout) :: this
        end subroutine BoundFinalize_t

        real(dp) function BoundGetPgen_t(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
            import :: ExcitationGenerator_t, n_int, dp, ScratchSize, maxExcit, NifTot, nEl
            class(ExcitationGenerator_t), intent(inout) :: this
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: ex(2, maxExcit), ic
            integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        end function BoundGetPgen_t
    end interface
contains

    !>  @brief
    !>  The excitation generator subroutine for singles and doubles
    subroutine gen_exc_sd(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                        ex, tParity, pGen, hel, store, part_type, &
                        singles_generator, doubles_generator)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        class(SingleExcitationGenerator_t), intent(inout) :: singles_generator
        class(DoubleExcitationGenerator_t), intent(inout) :: doubles_generator

        if (genrand_real2_dSFMT() < pSingles) then
            ic = 1
            call singles_generator%gen_exc(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                                ex, tParity, pGen, hel, store, part_type)
            pgen = pgen * pSingles
        else
            ic = 2
            call doubles_generator%gen_exc(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                                ex, tParity, pGen, hel, store, part_type)
            pgen = pgen * pDoubles
        end if
    end subroutine


    function get_pgen_sd(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2, &
                         singles_generator, doubles_generator) result(pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        class(SingleExcitationGenerator_t), intent(inout) :: singles_generator
        class(DoubleExcitationGenerator_t), intent(inout) :: doubles_generator
        real(dp) :: pgen

        if (ic == 1) then
            pgen = pSingles * singles_generator%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        else if (ic == 2) then
            pgen = (1.0 - pSingles) * doubles_generator%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        else
            pgen = 0.0_dp
        end if
    end function

    subroutine gen_all_excits_sd(nI, n_excits, det_list, singles_generator, doubles_generator)
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        class(SingleExcitationGenerator_t), intent(in) :: singles_generator
        class(DoubleExcitationGenerator_t), intent(in) :: doubles_generator

        integer :: idummy
        integer(n_int), allocatable :: singles(:, :), doubles(:, :)

        call singles_generator%gen_all_excits(nI, idummy, singles)
        call doubles_generator%gen_all_excits(nI, idummy, doubles)

        n_excits = size(singles, 2) + size(doubles, 2)
        allocate(det_list(0:niftot, n_excits))

        det_list(:, : size(singles, 2)) = singles(:, :)
        det_list(:, size(singles, 2) + 1 :) = doubles(:, :)
        call sort(det_list, ilut_lt, ilut_gt)
    end subroutine

    subroutine FCI_singles_gen_all_excits(this, nI, n_excits, det_list)
        class(SingleExcitationGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        unused_var(this)
        call gen_excits(nI, n_excits, det_list, ex_flag=1)
    end subroutine

    subroutine FCI_doubles_gen_all_excits(this, nI, n_excits, det_list)
        class(DoubleExcitationGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        unused_var(this)
        call gen_excits(nI, n_excits, det_list, ex_flag=2)
    end subroutine

    !!! ClassicAbInitExcitationGenerator_t methods !!!
    subroutine abinit_PCHB_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
        ex, tParity, pGen, hel, store, part_type)
        class(ClassicAbInitExcitationGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        call gen_exc_sd(nI, ilutI, nJ, ilutJ, exFlag, ic, &
        ex, tParity, pGen, hel, store, part_type, &
        this%singles_generator, this%doubles_generator)
    end subroutine abinit_PCHB_gen_exc

    real(dp) function abinit_PCHB_get_pgen(&
            this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) &
                result(pgen)
        class(ClassicAbInitExcitationGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        pgen = get_pgen_sd(&
                        nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2, &
                        this%singles_generator, this%doubles_generator)
    end function abinit_PCHB_get_pgen


    subroutine abinit_PCHB_gen_all_excits(this, nI, n_excits, det_list)
        class(ClassicAbInitExcitationGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        call gen_all_excits_sd(nI, n_excits, det_list, &
                               this%singles_generator, this%doubles_generator)
    end subroutine abinit_PCHB_gen_all_excits

    subroutine abinit_PCHB_finalize(this)
        class(ClassicAbInitExcitationGenerator_t), intent(inout) :: this
        call this%doubles_generator%finalize()
        call this%singles_generator%finalize()
        deallocate(this%singles_generator)
        deallocate(this%doubles_generator)
    end subroutine abinit_PCHB_finalize

end module
