#include "macros.h"
#:include "macros.fpph"

module gasci_discarding
    use constants, only: n_int, dp, maxExcit
    use SystemData, only: nel
    use bit_rep_data, only: NIfTot
    use sort_mod, only: sort

!     use pchb_excitgen, only: init_pchb_excitgen, finalize_pchb_excitgen, gen_rand_excit_pchb
    use pchb_excitgen, only: PCHB_FCI
    use gasci, only: GAS_specification, GASSpec_t
    use FciMCData, only: excit_gen_store_type
    implicit none

    private

    public :: gen_GASCI_discarding, init_GASCI_discarding, finalize_GASCI_discarding

contains

    !> This GAS excitation generator just uses a FCI excitation generator
    !> and discards excitations which are not in the GAS space.
    subroutine gen_GASCI_discarding(nI, ilutI, nJ, ilutJ, exFlag, ic, &
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
        character(*), parameter :: this_routine = 'generate_nGAS_excitation'

        integer :: src_copy(maxExcit)



        @:unused_var(exFlag, part_type, store)

#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        @:ASSERT(GAS_specification%contains_det(nI))

        call PCHB_FCI%gen_excit(nI, ilutI, nJ, ilutJ, ic, ex_mat, tParity, pgen, hel, store)

        if (nJ(1) /= 0) then
            if (.not. GAS_specification%contains_det(nJ)) then
                src_copy(:ic) = ex_mat(1, :ic)
                call sort(src_copy)
                ex_mat(1, :ic) = src_copy(:ic)
                ex_mat(2, :ic) = ex_mat(2, :ic)
                nJ(1) = 0
                ilutJ = 0_n_int
            end if
        end if
    end subroutine gen_GASCI_discarding

    subroutine init_GASCI_discarding()
        call PCHB_FCI%init()
    end subroutine

    subroutine finalize_GASCI_discarding()
        call PCHB_FCI%finalize()
    end subroutine

end module gasci_discarding
