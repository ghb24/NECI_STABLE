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

    public :: gen_GASCI_general

    public :: gen_exc_single


contains

    !>  @brief
    !>  The actual general GAS excitation generator.
    subroutine gen_GASCI_general(nI, ilutI, nJ, ilutJ, exFlag, ic, &
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


        @:unused_var(exFlag, part_type, store)

#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        if (genrand_real2_dSFMT()  >= pDoubles) then
            ic = 1
            call gen_exc_single(GAS_specification, nI, ilutI, &
                                nJ, ilutJ, ex_mat, tParity, pgen)
            pgen = pgen * (1.0_dp - pDoubles)
        else
            ic = 2
            call gen_exc_double(GAS_specification, nI, ilutI,&
                                nJ, ilutJ, ex_mat, tParity, pgen)
            pgen = pgen * pDoubles
        end if

        @:ASSERT(0.0_dp <= pgen .and. pgen <= 1.0_dp, pgen)
        @:ASSERT(all(ex_mat(:, :ic) /= UNKNOWN) .or. nJ(1) == 0)
    end subroutine gen_GASCI_general


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

    !>  @brief
    !>  Generate a double excitation under GAS constraints.
    subroutine gen_exc_double(GAS_spec, det_I, ilutI, nJ, ilutJ, ex_mat, par, pgen)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'gen_exc_double'

        type(DoubleExc_t) :: exc, reverted_exc
        integer, allocatable :: possible_holes(:)
        integer :: deleted(2)
        ! Spin of second electron
        type(SpinProj_t) :: m_s_1
        real(dp) :: pgen_particles, &
            ! These are arrays, because the pgen might be different if we
            ! pick AB or BA. Let i, j denote the particles and a, b the holes.
            ! pgen_particles == p (i j)
            ! pgen_first_pick == p(a | i j) == p(b | i j)
            ! pgen_second_pick == [p(b | a i j), p(a | b i j)]
            ! pgen == p_double * p (i j) * p(a | i j) * sum([p(b | a i j), p(a | b i j)])
            !      == p_double * p (i j) * p(b | i j) * sum([p(b | a i j), p(a | b i j)])
            pgen_first_pick, pgen_second_pick(2)
        real(dp) :: r
        real(dp), allocatable :: c_sum(:)
        integer :: i

        integer :: elecs(2), sym_product, ispn, sum_ml

        call pick_biased_elecs(det_I, elecs, exc%val(1, :), &
                               sym_product, ispn, sum_ml, pgen_particles)
        @:ASSERT(exc%val(1, 1) /= exc%val(2, 1), exc%val)

        deleted = exc%val(1, :)
        ! Get possible holes for the first particle, while fullfilling and spin- and GAS-constraints.
        ! and knowing that a second particle will be created afterwards!
        possible_holes = get_possible_holes(&
            GAS_spec, det_I, add_holes=deleted, n_total=2, excess=-sum(calc_spin_raw(deleted)))
        @:ASSERT(disjoint(possible_holes, det_I))

        if (size(possible_holes) == 0) then
            pgen = pgen_particles
            call zeroResult()
            return
        end if

        r = genrand_real2_dSFMT()
        ! Pick randomly one hole with arbitrary spin
        exc%val(2, 1) = possible_holes(int(r * real(size(possible_holes), kind=dp)) + 1)
        pgen_first_pick = 1.0_dp / real(size(possible_holes), dp)
        m_s_1 = calc_spin_raw(exc%val(2, 1))
        @:ASSERT(any(m_s_1 == [alpha, beta]))


        ! Pick second hole.
        ! The total spin projection of the created particles has to add up
        ! to the total spin projection of the deleted particles.
        ! Get possible holes for the second particle,
        ! while fullfilling GAS- and Spin-constraints.
        block
            type(SpinProj_t) :: excess

            excess = m_s_1 - sum(calc_spin_raw(deleted))

            @:ASSERT(abs(excess%val) <= 1)

            possible_holes = get_possible_holes(&
                    GAS_spec, det_I, add_holes=deleted, &
                    add_particles=exc%val(2:2, 1), &
                    n_total=1, excess=excess)
        end block

        @:ASSERT(disjoint(possible_holes, exc%val(2:2, 1)))
        @:ASSERT(disjoint(possible_holes, det_I))

        if (size(possible_holes) == 0) then
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt2 in possible_holes
        c_sum = get_cumulative_list(det_I, exc, possible_holes)
        call draw_from_cum_list(c_sum, i, pgen_second_pick(1))

        if (i /= 0) then
            exc%val(2, 2) = possible_holes(i)
        else
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if
        @:ASSERT(defined(exc))

        ! We could have picked the holes the other way round and have to
        ! determine the probability of picking tgt1 with spin m_s_1 upon picking tgt2 first.
        block
            integer :: src1, tgt1, src2, tgt2
            src1 = exc%val(1, 1); tgt1 = exc%val(2, 1)
            src2 = exc%val(1, 2); tgt2 = exc%val(2, 2)

            @:ASSERT(src1 /= tgt2 .and. src2 /= tgt2, src1, tgt2, src2)
            possible_holes = get_possible_holes(&
                    GAS_spec, det_I, add_holes=deleted, &
                    add_particles=exc%val(2, 2:2), &
                    n_total=1, excess=calc_spin_raw(tgt2) - sum(calc_spin_raw(deleted)))
            @:ASSERT(disjoint(possible_holes, exc%val(2, 2:2)))
            @:ASSERT(disjoint(possible_holes, det_I))

            if (size(possible_holes) == 0) then
                pgen = pgen_particles * pgen_first_pick * pgen_second_pick(1)
                call zeroResult()
                return
            end if
            ! Possible_holes has to contain tgt1.
            ! we look up its index with binary search
            i = binary_search_first_ge(possible_holes, tgt1)
            @:ASSERT(i /= -1, tgt1, possible_holes)

            reverted_exc = DoubleExc_t(src1=src1, tgt1=tgt2, src2=src2)
            c_sum = get_cumulative_list(det_I, reverted_exc, possible_holes)
            if (i == 1) then
                pgen_second_pick(2) = c_sum(1)
            else
                pgen_second_pick(2) = (c_sum(i) - c_sum(i - 1))
            end if

            if (i /= 0) then
                call make_double(det_I, nJ, elecs(1), elecs(2), tgt1, tgt2, ex_mat, par)
                ilutJ = ilut_excite(ilutI, exc)
            else
                ilutJ = 0
            end if
        end block

        pgen = pgen_particles * pgen_first_pick * sum(pgen_second_pick)

        @:ASSERT(0.0_dp < pgen .and. pgen <= 1.0_dp)
        @:ASSERT(all(ex_mat(:, :2) /= UNKNOWN))

        contains

            subroutine zeroResult()
                integer :: src_copy(2)

                src_copy(:) = exc%val(1, :)
                call sort(src_copy)
                ex_mat(1, :2) = src_copy
                ex_mat(2, :2) = exc%val(2, :)
                nJ(1) = 0
                ilutJ = 0_n_int
            end subroutine zeroResult
    end subroutine gen_exc_double
end module gasci_general
