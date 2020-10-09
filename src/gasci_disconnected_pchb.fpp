#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_disconnected_pchb
    use constants, only: n_int, dp, int64, maxExcit
    use util_mod, only: intswap
    use SystemData, only: nel, tGASSpinRecoupling
    use bit_rep_data, only: NIfTot
    use excitation_types, only: DoubleExc_t
    use orb_idx_mod, only: calc_spin_raw, SpinProj_t, operator(==)

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
        private
        generic, public :: init => init_GAS_pchb_excitgen
        procedure, public :: clear => clear_GAS_pchb_excitgen

        procedure :: init_GAS_pchb_excitgen

        procedure, nopass, public :: is_allowed => GAS_allowed
    end type

    !>  @brief
    !>  The disconnected_GAS_PCHB excitation generator subroutine.
    !>
    !>  @details
    !>  This is a wrapper around `disconnected_GAS_PCHB%gen_excit`
    !>  to match the function pointer interface.
    type(PCHB_GAS_excit_generator_t) :: disconnected_GAS_PCHB

contains

    !>  @brief
    !>  The disconnected_GAS_PCHB excitation generator subroutine.
    !>
    !>  @details
    !>  This is a wrapper around `disconnected_GAS_PCHB%gen_excit`
    !>  to match the function pointer interface.
    !>  The interface is common to all excitation generators, see proc_ptrs.F90
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

        ! It is not possible to define `wrapper_gen_GAS_single` as
        ! internal procedure, because the lifetime of internal procedures
        ! ends with the scope of the defining procedure.
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

        integer :: src_spaces(2), tgt_spaces(2)

        src_spaces = GAS_specification%get_iGAS(exc%val(1, :))
        tgt_spaces = GAS_specification%get_iGAS(exc%val(2, :))

        if (all(src_spaces == tgt_spaces) .and. src_spaces(1) == src_spaces(2)) then
            ! All electrons come from the same space and there are no restrictions
            ! regarding recoupling.
            GAS_allowed = .true.
        else if (tGASSpinRecoupling) then
            ! Ensure that the number of particles per GAS space stays constant.
            if (src_spaces(1) > src_spaces(2)) call intswap(src_spaces(1), src_spaces(2))
            if (tgt_spaces(1) > tgt_spaces(2)) call intswap(tgt_spaces(1), tgt_spaces(2))
            GAS_allowed = all(src_spaces == tgt_spaces)
        else
            ! Ensure that the number of particles and the spin projection
            ! per GAS space stays constant.
            block
                type(SpinProj_t) :: src_spins(2), tgt_spins(2)
                #:set spin_swap = functools.partial(swap, 'SpinProj_t', "", 0)
                src_spins = calc_spin_raw(exc%val(1, :))
                tgt_spins = calc_spin_raw(exc%val(2, :))
                if (src_spaces(1) > src_spaces(2)) then
                    call intswap(src_spaces(1), src_spaces(2))
                    @:spin_swap(src_spins(1), src_spins(2))
                end if
                if (tgt_spaces(1) > tgt_spaces(2)) then
                    call intswap(tgt_spaces(1), tgt_spaces(2))
                    @:spin_swap(tgt_spins(1), tgt_spins(2))
                end if

                GAS_allowed = all((src_spaces == tgt_spaces) .and. (src_spins == tgt_spins))
            end block
        end if
    end function
end module gasci_disconnected_pchb
