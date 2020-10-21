#include "macros.h"
module pchb_excitgen
    use constants, only: n_int, dp, maxExcit
    use SystemData, only: nel, t_pchb_weighted_singles
    use bit_rep_data, only: NIfTot
    use GenRandSymExcitNUMod, only: uniform_single_excit_wrapper
    use excit_gens_int_weighted, only: gen_single_4ind_ex, pgen_single_4ind
    use FciMCData, only: pSingles, excit_gen_store_type, pDoubles
    use SymExcitDataMod, only: ScratchSize
    use pchb_factory, only: PCHB_excitation_generator_t
    use GenRandSymExcitNUMod, only: calc_pgen_symrandexcit2
    use excitation_types, only: DoubleExc_t
    implicit none

    private

    public :: gen_rand_excit_pchb, PCHB_FCI

    type, extends(PCHB_excitation_generator_t) :: PCHB_FCI_excit_generator_t
        private
    contains
        private
        procedure :: init_FCI_pchb_excitgen
        generic, public :: init => init_FCI_pchb_excitgen

        procedure, nopass, public :: is_allowed => all_allowed
    end type

    type(PCHB_FCI_excit_generator_t) :: PCHB_FCI
contains

    logical pure function all_allowed(exc)
        type(DoubleExc_t), intent(in) :: exc
        unused_var(exc)
        all_allowed = .true.
    end function


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

        call PCHB_FCI%gen_excit(nI, ilutI, nJ, ilutJ, ic, ex, tpar, pgen, helgen, store)
    end subroutine gen_rand_excit_pchb


    subroutine init_FCI_pchb_excitgen(this)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this

        ! It is not possible to define the wrapper functions as
        ! internal procedure, because the lifetime of internal procedures
        ! ends with the scope of the defining procedure.
        if (t_pchb_weighted_singles) then
            call this%init(weighted_single_excit_wrapper, calc_pgen_weighted_single)
        else
            call this%init(gen_uniform_single, calc_pgen_uniform_single)
        end if
    end subroutine

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
