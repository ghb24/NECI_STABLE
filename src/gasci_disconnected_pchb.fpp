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

    use pchb_excitgen, only: init_pchb_excitgen, finalize_pchb_excitgen, gen_rand_excit_pchb

!     use util_mod, only: fuseIndex, linearIndex, intswap, getSpinIndex, near_zero

    use gasci, only: GAS_specification, GASSpec_t
    use FciMCData, only: excit_gen_store_type, projEDet

    implicit none

    private

    public :: gen_GASCI_pchb, init_GASCI_pchb, finalize_GASCI_pchb

contains

    !> This GAS excitation generator just uses a FCI excitation generator
    !> and discards excitations which are not in the GAS space.
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

        integer :: src_copy(maxExcit)



        @:unused_var(exFlag, part_type, store)

#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        @:ASSERT(GAS_specification%contains(nI))

        call gen_rand_excit_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex_mat, tParity, pGen, hel, store, part_type)

    end subroutine gen_GASCI_pchb

    subroutine init_GASCI_pchb(projEDet)
        integer, intent(in) :: projEDet(:)
        call init_pchb_excitgen(projEDet(:), GAS_allowed)
        contains

            logical pure function GAS_allowed(exc)
                type(DoubleExc_t), intent(in) :: exc

                integer :: GAS_spaces(2, 2)

                GAS_spaces = GAS_specification%get_iGAS(exc%val)

                @:sort(integer, GAS_spaces(1, :), <=)
                @:sort(integer, GAS_spaces(2, :), <=)

                GAS_allowed = all(GAS_spaces(1, :) == GAS_spaces(2, :))
            end function
    end subroutine

    subroutine finalize_GASCI_pchb()
        call finalize_pchb_excitgen()
    end subroutine


end module gasci_disconnected_pchb
