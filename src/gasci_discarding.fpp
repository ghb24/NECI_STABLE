#include "macros.h"
#:include "macros.fpph"

module gasci_discarding
    use constants, only: n_int, dp, maxExcit
    use SystemData, only: nel
    use FciMCData, only: excit_gen_store_type
    use bit_rep_data, only: NIfTot
    use sort_mod, only: sort
    use SymExcitDataMod, only: ScratchSize

    use pchb_excitgen, only: PCHB_FCI_excit_generator_t
    use excitation_generators, only: ExcitationGenerator_t, SingleExcitationGenerator_t, DoubleExcitationGenerator_t
    use gasci, only: GASSpec_t
    use gasci_util, only: GAS_gen_all_excits => gen_all_excits
    implicit none

    private

    public :: GAS_DiscardingGenerator_t

    type, extends(ExcitationGenerator_t) :: GAS_DiscardingGenerator_t
        private
        type(PCHB_FCI_excit_generator_t) :: FCI_generator
        class(GASSpec_t), allocatable :: GAS_spec
    contains
        private
        procedure, public :: init
        procedure, public :: finalize
        procedure, public :: gen_exc
        procedure, public :: get_pgen
        procedure, public :: gen_all_excits
    end type
contains

    !> This GAS excitation generator just uses a FCI excitation generator
    !> and discards excitations which are not in the GAS space.
    subroutine gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                       ex, tParity, pGen, hel, store, part_type)
        class(GAS_DiscardingGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'generate_nGAS_excitation'

        integer :: src_copy(maxExcit)

#ifdef WARNING_WORKAROUND_
        hel = h_cast(0.0_dp)
#endif
        ASSERT(this%GAS_spec%contains_det(nI))

        call this%FCI_generator%gen_exc(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex, tParity, pGen, hel, store, part_type)
        if (nJ(1) /= 0) then
            if (.not. this%GAS_spec%contains_det(nJ)) then
                src_copy(:ic) = ex(1, :ic)
                call sort(src_copy)
                ex(1, :ic) = src_copy(:ic)
                ex(2, :ic) = ex(2, :ic)
                nJ(1) = 0
                ilutJ = 0_n_int
            end if
        end if
    end subroutine

    function get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_DiscardingGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        pgen = this%FCI_generator%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
    end function


    subroutine init(this, GAS_spec)
        class(GAS_DiscardingGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        unused_var(this)
        this%GAS_spec = GAS_spec
        call this%FCI_generator%init()
    end subroutine

    subroutine finalize(this)
        class(GAS_DiscardingGenerator_t), intent(inout) :: this
        unused_var(this)
        call this%FCI_generator%finalize()
    end subroutine


    subroutine gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_DiscardingGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        call GAS_gen_all_excits(this%GAS_spec, nI, n_excits, det_list)
    end subroutine

end module gasci_discarding
