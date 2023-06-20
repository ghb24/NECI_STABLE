#include "macros.h"
module exc_gen_class_wrappers
    use constants, only: dp, n_int, maxExcit
    use excitation_generators, only: FCISingleExcitationGenerator_t
    use FciMCData, only: pSingles, pDoubles, excit_gen_store_type
    use SystemData, only: nel
    use bit_rep_data, only: NIfTot
    use SymExcitDataMod, only: ScratchSize
    use util_mod, only: stop_all

    use GenRandSymExcitNUMod, only: uniform_single_excit_wrapper, calc_pgen_symrandexcit2
    use excit_gens_int_weighted, only: gen_single_4ind_ex, pgen_single_4ind
    implicit none
    private
    public :: UniformSingles_t, WeightedSingles_t

    type, extends(FCISingleExcitationGenerator_t) :: UniformSingles_t
    contains
        private
        procedure, public :: gen_exc => UniformSingles_gen_exc
        procedure, public :: get_pgen => UniformSingles_get_pgen
        procedure, public :: finalize => UniformSingles_do_nothing
    end type

    type, extends(FCISingleExcitationGenerator_t) :: WeightedSingles_t
    contains
        private
        procedure, public :: gen_exc => WeightedSingles_gen_exc
        procedure, public :: get_pgen => WeightedSingles_get_pgen
        procedure, public :: finalize => WeightedSingles_do_nothing
    end type

contains

    subroutine UniformSingles_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                             ex, tParity, pGen, hel, store, part_type)
        class(UniformSingles_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        unused_var(this); unused_var(exFlag); unused_var(part_type)
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0.0_dp)
#endif
        ic = 1
        call uniform_single_excit_wrapper(nI, ilutI, nJ, ilutJ, ex, tParity, store, pgen)
        pgen = pgen / pSingles
    end subroutine

    function UniformSingles_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(UniformSingles_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        unused_var(this); unused_var(ilutI)
        call calc_pgen_symrandexcit2(nI, ex, ic, ClassCount2, ClassCountUnocc2, pDoubles, pgen)
        pgen = pgen / pSingles
    end function

    subroutine UniformSingles_do_nothing(this)
        class(UniformSingles_t), intent(inout) :: this
        unused_var(this)
    end subroutine

    subroutine WeightedSingles_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                             ex, tParity, pGen, hel, store, part_type)
        class(WeightedSingles_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        unused_var(this); unused_var(store); unused_var(exFlag); unused_var(part_type);
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0.0_dp)
#endif
        ic = 1
        call gen_single_4ind_ex(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
    end subroutine

    function WeightedSingles_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(WeightedSingles_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'WeightedSingles_get_pgen'
        unused_var(this); unused_var(ClassCount2); unused_var(ClassCountUnocc2);
        unused_var(ilutI)
        ASSERT(ic == 1)
        pgen = pgen_single_4ind(nI, ilutI, ex(1, 1), ex(2, 1))
    end function

    subroutine WeightedSingles_do_nothing(this)
        class(WeightedSingles_t), intent(inout) :: this
        unused_var(this)
    end subroutine
end module
