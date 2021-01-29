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

    !>  @brief
    !>  Generate a single excitation under GAS constraints.
    subroutine gen_exc_single(GAS_spec, det_I, ilutI, nJ, ilutJ, ex_mat, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_exc_SingleExc_t'

        type(SingleExc_t) :: exc
        integer, allocatable :: possible_holes(:)
        real(dp) :: pgen_particle, pgen_hole
        real(dp), allocatable :: c_sum(:)
        integer :: i, elec

        ! Pick any random electron
        elec = int(genrand_real2_dSFMT() * nEl) + 1
        exc%val(1) = det_I(elec)
        pgen_particle = 1.0_dp / real(nEl, kind=dp)

        ! Get a hole with the same spin projection
        ! while fullfilling GAS-constraints.
        block
            integer :: deleted(1)
            deleted(1) = exc%val(1)
            possible_holes = get_possible_holes(&
                GAS_spec, det_I, add_holes=deleted, excess=-calc_spin_raw(exc%val(1)))
        end block


        if (size(possible_holes) == 0) then
            call zero_result()
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt \in possible_holes
        c_sum = get_cumulative_list(det_I, exc, possible_holes)
        call draw_from_cum_list(c_sum, i, pgen_hole)
        @:ASSERT(i == 0 .neqv. (0.0_dp < pgen_hole .and. pgen_hole <= 1.0_dp))
        if (i /= 0) then
            exc%val(2) = possible_holes(i)
            call make_single(det_I, nJ, elec, exc%val(2), ex_mat, par)
            ilutJ = ilut_excite(ilutI, exc)
        else
            call zero_result()
        end if

        pgen = pgen_particle * pgen_hole
        @:ASSERT(nJ(1) == 0 .neqv. 0.0_dp < pgen .and. pgen <= 1.0_dp)
        contains

            subroutine zero_result()
                ex_mat(:, 1) = exc%val(:)
                nJ(1) = 0
                ilutJ = 0_n_int
            end subroutine zero_result
    end subroutine
end module gasci_general
