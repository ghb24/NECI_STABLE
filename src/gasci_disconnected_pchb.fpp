#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_disconnected_pchb
    use constants, only: n_int, dp, int64, maxExcit
    use SystemData, only: nel, tGASSpinRecoupling
    use bit_rep_data, only: NIfTot
    use sort_mod, only: sort
    use excitation_types, only: DoubleExc_t
    use orb_idx_mod, only: calc_spin_raw


    use SymExcitDataMod, only: scratchSize
    use GenRandSymExcitNUMod, only: uniform_single_excit_wrapper, calc_pgen_symrandexcit2
    use FciMCData, only: pSingles, pDoubles

    use pchb_factory, only: PCHB_excitation_generator_t


    use gasci, only: GAS_specification, GASSpec_t
    use gasci_disconnected, only: generate_nGAS_single, init_disconnected_GAS, clearGAS
    use FciMCData, only: excit_gen_store_type

    implicit none

    private

    public :: gen_GASCI_pchb, disconnected_GAS_PCHB

    type, extends(PCHB_excitation_generator_t) :: PCHB_GAS_excit_generator_t
        private
    contains
        procedure :: init_GAS_pchb_excitgen
        generic, public :: init => init_GAS_pchb_excitgen

        procedure :: clear => clear_GAS_pchb_excitgen

        procedure, nopass :: is_allowed => GAS_allowed
    end type

    type(PCHB_GAS_excit_generator_t) :: disconnected_GAS_PCHB

contains

    subroutine gen_GASCI_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex_mat, tParity, pGen, hel, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'gen_GASCI_pchb'

        @:unused_var(exFlag, part_type, store)
        @:ASSERT(GAS_specification%contains(nI))

        call disconnected_GAS_PCHB%gen_excit(nI, ilutI, nJ, ilutJ, ic, ex_mat, tParity, pgen, hel, store)
    end subroutine gen_GASCI_pchb

    subroutine wrapper_gen_GAS_single(nI, ilutI, nJ, ilutJ, ex, par, store, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        type(excit_gen_store_type), intent(inout), target :: store
        real(dp), intent(out) :: pgen
        @:unused_var(store)

        call generate_nGAS_single(nI, ilutI, nJ, ilutJ, ex, par, pgen)
    end subroutine

    subroutine init_GAS_pchb_excitgen(this)
        class(PCHB_GAS_excit_generator_t), intent(inout) :: this

        call this%init(wrapper_gen_GAS_single)
        call init_disconnected_GAS(GAS_specification)
    end subroutine

    subroutine clear_GAS_pchb_excitgen(this)
        class(PCHB_GAS_excit_generator_t), intent(inout) :: this

        call this%finalize()
        call clearGAS()
    end subroutine

    logical pure function GAS_allowed(exc)
        type(DoubleExc_t), intent(in) :: exc

        integer :: GAS_spaces(2, 2)

        GAS_spaces = GAS_specification%get_iGAS(exc%val)

        @:sort(integer, GAS_spaces(1, :), <=)
        @:sort(integer, GAS_spaces(2, :), <=)

        GAS_allowed = all(GAS_spaces(1, :) == GAS_spaces(2, :))
    end function
end module gasci_disconnected_pchb
