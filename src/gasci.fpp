#include "macros.h"
#:include "macros.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module gasci
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
        near_zero
    use sort_mod, only : sort
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        last_tgt_unknown, set_last_tgt
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    implicit none

    private
    public :: isValidExcit, loadGAS, generate_nGAS_excitation, clearGAS, &
      GAS_specification, nGAS

    !> Speficies the GAS spaces.
    !> It is assumed, that the GAS spaces are contigous.
    !> The indices are:
    !>  n_orbs_per_GAS(1:nGAS), n_min(1:nGAS), n_max(1:nGAS)
    !> n_orbs_per_GAS(iGAS, iRep) specifies how many orbitals are in
    !> the `iGAS` GAS space in the `iRep` Irrep.
    !> n_min(iGAS) specifies the cumulated! minimum particle number per GAS space.
    !> n_max(iGAS) specifies the cumulated! maximum particle number per GAS space.
    type :: GAS_specification_t
        integer, allocatable :: n_orbs_per_GAS(:), n_min(:), n_max(:)
    end type

    type(GAS_specification_t) :: GAS_specification



contains

    pure function get_nGAS(GAS_specification) result(res)

    end function

    pure function is_valid(GAS_specification) result(res)
        type(GAS_specification_t), intent(in) :: GAS_specification
        logical :: res

        logical :: shapes_match, nEl_correct

        shapes_match = size(GAS_specification%n_min) == size(GAS_specification%n_max) &
                       .and. size(GAS_specification%n_orbs_per_GAS, 1) == size(GAS_specification%n_min)
                       .and. size(GAS_specification%n_orbs_per_GAS, 1) == nGAS
! TODO: insert
!                        .and. size(GAS_specification%n_orbs_per_GAS, 2) == nSym
        nEl_correct = GAS_specification%n_min(nGAS) == GAS_specification%n_max(nGAS) &
                      .and. GAS_specification%n_min(nGAS) == nEl

        res = shapes_match .and. nEl_correct
    end function

    pure function is_disconnected(GAS_specification) result(res)
        type(GAS_specification_t), intent(in) :: GAS_specification
        logical :: res
        res = all(GAS_specification%n_min(:) == GAS_specification%n_max(:))
    end function


end module gasci
