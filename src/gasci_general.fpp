#include "macros.h"
#:include "macros.fpph"

module gasci_general
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
        near_zero, cumsum, operator(.isclose.), lex_leq, stop_all
    use sort_mod, only: sort
    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet
    use sets_mod, only: is_sorted, complement, union, disjoint
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pSingles, pDoubles, excit_gen_store_type
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        get_last_tgt, set_last_tgt, defined, UNKNOWN, &
        excite, ilut_excite
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, alpha, beta, &
        operator(-), operator(==), operator(/=), sum
    use growing_buffers, only: buffer_int_2D_t

    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit

    use gasci, only: GAS_specification, GASSpec_t

    use gasci_util, only: draw_from_cum_list, get_cumulative_list, get_possible_spaces, get_possible_holes
    implicit none

    private

    public :: gen_exc_single


contains

end module gasci_general
