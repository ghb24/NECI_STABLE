#include "macros.h"
module pchb_excitgen
    use constants, only: n_int, dp, maxExcit
    use util_mod, only: operator(.div.)
    use SystemData, only: nel, nBasis, t_pchb_weighted_singles
    use bit_rep_data, only: NIfTot
    use GenRandSymExcitNUMod, only: uniform_single_excit_wrapper
    use excit_gens_int_weighted, only: gen_single_4ind_ex, pgen_single_4ind
    use gasci, only: GASSpec_t
    use gasci_general_pchb, only: GAS_PCHB_excit_gen_t
    use FciMCData, only: pSingles, excit_gen_store_type, pDoubles
    use SymExcitDataMod, only: ScratchSize
    use pchb_factory, only: PCHB_excitation_generator_t
    use GenRandSymExcitNUMod, only: calc_pgen_symrandexcit2
    use excitation_types, only: DoubleExc_t
    implicit none

    private

    public :: gen_rand_excit_pchb, PCHB_FCI

    type :: PCHB_FCI_excit_generator_t
        private
        type(GAS_PCHB_excit_gen_t) :: GAS_exc_generator
    contains
        private
        procedure, public :: init => init_FCI_pchb_excitgen
        procedure, public :: finalize
        procedure, public :: gen_excit
        procedure, public :: calc_pgen => calc_pgen_pchb
    end type

    type(PCHB_FCI_excit_generator_t) :: PCHB_FCI


contains


    !>  @brief
    !>  The excitation generator subroutine for Full CI PCHB.
    !>
    !>  @details
    !>  This is a wrapper to match the function pointer interface.
    !>  The interface is common to all excitation generators, see proc_ptrs.F90
    !>
    !>  For singles, use the uniform or weighted excitgen.
    subroutine gen_rand_excit_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tpar, &
                                   pgen, helgen, store, part_type)
        ! The interface is common to all excitation generators, see proc_ptrs.F90
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        unused_var(exFlag); unused_var(part_type)
        helgen = h_cast(0._dp)

        call PCHB_FCI%gen_excit(nI, ilutI, nJ, ilutJ, ic, ex, tpar, store, pgen)
    end subroutine gen_rand_excit_pchb


    subroutine init_FCI_pchb_excitgen(this)
        class(PCHB_FCI_excit_generator_t), intent(out) :: this

        ! It is not possible to define the wrapper functions as
        ! internal procedure, because the lifetime of internal procedures
        ! ends with the scope of the defining procedure.
        if (t_pchb_weighted_singles) then
            call this%GAS_exc_generator%init(CAS(nEl, nBasis .div. 2), weighted_single_excit_wrapper, calc_pgen_weighted_single)
        else
            call this%GAS_exc_generator%init(CAS(nEl, nBasis .div. 2), gen_uniform_single, calc_pgen_uniform_single)
        end if

    contains

        !> Create a GAS specification with one GAS space hence CAS.
        type(GASSpec_t) pure function CAS(n_el, n_spat_orb)
            integer, intent(in) :: n_el, n_spat_orb
            integer :: i
            CAS = GASSpec_t(n_min=[n_el], n_max=[n_el], &
                                 spat_GAS_orbs=[(1, i = 1, n_spat_orb)])
        end function
    end subroutine

    subroutine finalize(this)
        class(PCHB_FCI_excit_generator_t) :: this
        call this%GAS_exc_generator%finalize()
    end subroutine


    subroutine gen_excit(this, nI, ilutI, nJ, ilutJ, ic, ex, tpar, store, pgen)
        class(PCHB_FCI_excit_generator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        type(excit_gen_store_type), intent(inout), target :: store
        real(dp), intent(out) :: pGen
        call this%GAS_exc_generator%gen_excit(nI, ilutI, nJ, ilutJ, ic, ex, tpar, store, pgen)
    end subroutine


    function calc_pgen_pchb(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(PCHB_FCI_excit_generator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        pgen = this%GAS_exc_generator%calc_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
    end function calc_pgen_pchb


    !> Wrapper function to create a weighted single excitation
    subroutine weighted_single_excit_wrapper(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pGen
        type(excit_gen_store_type), intent(inout), target :: store

        unused_var(store)
        ! Call the 4ind-weighted single excitation generation
        call gen_single_4ind_ex(nI, ilutI, nJ, ilutJ, ex, tpar, pgen)
    end subroutine weighted_single_excit_wrapper

    !> Wrapper function to calculate pgen for weighted single excitation
    function calc_pgen_weighted_single(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        unused_var(ClassCount2); unused_var(ClassCountUnocc2); unused_var(ic)
        pgen = pgen_single_4ind(nI, ilutI, ex(1, 1), ex(2, 1))
    end function

    !> Wrapper function for creating a uniform single excitation
    subroutine gen_uniform_single(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pGen
        type(excit_gen_store_type), intent(inout), target :: store

        call uniform_single_excit_wrapper(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
        pgen = pgen / pSingles
    end subroutine

    !> Wrapper function to calculate pgen for uniform single excitation
    function calc_pgen_uniform_single(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        unused_var(ilutI)
        call calc_pgen_symrandexcit2(nI, ex, ic, ClassCount2, ClassCountUnocc2, pDoubles, pGen)
        pgen = pgen / pSingles
    end function
end module pchb_excitgen
