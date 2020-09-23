#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_general_pchb
    use constants, only: n_int, dp, int64, maxExcit
!     use util_mod, only: choose, stop_all, cumsum, binary_search_first_ge
    use FciMCData, only: excit_gen_store_type
    use SystemData, only: nEl
    use bit_rep_data, only: NIfTot

    use gasci, only: GAS_specification, GASSpec_t

    implicit none

    private

    public :: gen_general_GASCI_pchb

contains

    !>  @brief
    !>  The disconnected_GAS_PCHB excitation generator subroutine.
    !>
    !>  @details
    !>  This is a wrapper around `disconnected_GAS_PCHB%gen_excit`
    !>  to match the function pointer interface.
    !>  The interface is common to all excitation generators, see proc_ptrs.F90
    subroutine gen_general_GASCI_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
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


        @:unused_var(nI, ilutI, nJ, ilutJ, exFlag, ic, ex_mat, tParity, pGen, hel, store, part_type)

    end subroutine gen_general_GASCI_pchb


end module gasci_general_pchb
