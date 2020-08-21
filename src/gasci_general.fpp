#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

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
    implicit none

    private

    public :: gen_GASCI_general, gen_all_excits
    public :: get_possible_spaces, get_possible_holes, &
        get_available_singles, get_available_doubles


    interface get_cumulative_list
        #:for Excitation_t in ExcitationTypes
            module procedure get_cumulative_list_${Excitation_t}$
        #:endfor
    end interface

contains

    !>  @brief
    !>      Return the GAS spaces, where one particle can be created.
    !>
    !>  @details
    !>  It can be proven, that the spaces where a particle can
    !>  be created are contigous, so a two element integer array
    !>  [lower_bound, upper_bound] is returned.
    !>  As long as lower_bound <= iGAS <= upper_bound is true,
    !>  the created particle will lead to a valid Slater-Determinant.
    !>
    !>  It is **assumed** that the input determinant is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  Checks are only performed in DEBUG compilation mode and
    !>  the return value is undefined, if this is not the case!
    !>
    !>  It is possible to delete additional particles with the
    !>  optional argument `add_holes` before checking the
    !>  validity of particle creation.
    !>  On the other hand it is possible to create particles before
    !>  checking the validity of particle creation.
    !>
    !>  If more than one particle should be created, the optional argument
    !>  n_total (default 1) should be used.
    !>
    !>  If no creation is allowed by the GAS constraints,
    !>  the bounds will be returned as integer constant EMPTY_BOUNDS.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] particles_per_GAS, The particles per GAS space.
    !>  @param[in] add_holes, optional, An index of orbitals
    !>      where particles should be deleted before creating the new particle.
    !>  @param[in] add_particles, optional, An index of orbitals
    !>      where particles should be created before creating the new particle.
    !>  @param[in] n_total, optional, The total number of particles
    !>      that will be created. Defaults to one (integer).
    function get_possible_spaces(GAS_spec, particles_per_GAS, add_holes, add_particles, n_total) result(spaces)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: &
            particles_per_GAS(GAS_spec%nGAS())
        integer, intent(in), optional :: add_holes(:), add_particles(:)
        integer, intent(in), optional :: n_total
        character(*), parameter :: this_routine = 'get_possible_spaces'
        integer :: spaces(2)

        !> Lower and upper bound for spaces where a particle can be created.
        !> If no particle can be created, then spaces == 0 .
        integer :: n_total_, i, iGAS, lower_bound, upper_bound

        integer :: &
        !> Cumulated number of particles per iGAS
            cum_n_particle(GAS_spec%nGAS()), &
        !> Cumulated deficit per iGAS
            deficit(GAS_spec%nGAS()), &
        !> Cumulated vacant orbitals per iGAS
            vacant(GAS_spec%nGAS())

        @:def_default(n_total_, n_total, 1)

        block
            integer :: B(GAS_spec%nGAS()), C(GAS_spec%nGAS())
            if (present(add_holes)) then
                B = GAS_spec%count_per_GAS(add_holes)
            else
                B = 0
            end if
            if (present(add_particles)) then
                C = GAS_spec%count_per_GAS(add_particles)
            else
                C = 0
            end if
            cum_n_particle = cumsum(particles_per_GAS - B + C)
        end block

        deficit = GAS_spec%cn_min(:) - cum_n_particle(:)
        vacant = GAS_spec%cn_max(:) - cum_n_particle(:)

        if (any(n_total_ < deficit) .or. all(vacant < n_total_)) then
            spaces = 0
            return
        end if

        ! Find the first index, where a particle has to be created.
        do iGAS = 1, GAS_spec%nGAS()
            if (deficit(iGAS) == n_total_) exit
        end do
        upper_bound = iGAS

        ! We assume that it is possible to create a particle at least in
        ! the last GAS space.
        ! Search from behind the first occurence where it is not possible
        ! anymore to create a particle.
        ! The lower bound is one GAS index above.
        do iGAS = GAS_spec%nGAS(), 1, -1
            if (vacant(iGAS) <= 0) exit
        end do
        lower_bound = iGAS + 1

        if (lower_bound > upper_bound .or. lower_bound > GAS_spec%nGAS()) then
            spaces = 0
        else
            spaces = [lower_bound, upper_bound]
        end if
    end function


    function get_possible_holes(GAS_spec, det_I, add_holes, add_particles, n_total, excess) result(possible_holes)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, intent(in), optional :: add_holes(:)
        integer, intent(in), optional :: add_particles(:)
        integer, intent(in), optional :: n_total
        type(SpinProj_t), intent(in), optional :: excess
        character(*), parameter :: this_routine = 'get_possible_holes'

        integer, allocatable :: possible_holes(:)

        integer :: &
            splitted(GAS_spec%max_GAS_size(), GAS_spec%nGAS()), &
            splitted_sizes(GAS_spec%nGAS())

        integer :: spaces(2)
        integer :: n_total_

        @:def_default(n_total_, n_total, 1)
        @:ASSERT(1 <= n_total_)
        if (present(excess)) then
            @:ASSERT(abs(excess%val) <= n_total_, excess, n_total_)
        end if

        call GAS_spec%split_per_GAS(det_I, splitted, splitted_sizes)

        ! Note, that non-present optional arguments can be passed
        ! into optional arguments without checking!
        spaces = get_possible_spaces(&
             GAS_spec, splitted_sizes, add_holes=add_holes, &
             add_particles=add_particles, n_total=n_total_)

        if (all(spaces == 0)) then
            possible_holes = [integer::]
            return
        end if

        block
            integer :: i, k, iGAS, incr, curr_value, iGAS_min_val
            integer :: L(spaces(1) : spaces(2)), counter(spaces(1) : spaces(2))
            integer, allocatable :: possible_values(:)
            type(SpinProj_t) :: m_s

            m_s = SpinProj_t(0)
            if (present(excess)) then
                if (abs(excess%val) == n_total_) m_s = -SpinProj_t(sign(1, excess%val))
            end if


            L = GAS_spec%GAS_size([(i, i = spaces(1), spaces(2))])
            if (m_s == beta) then
                allocate(possible_values(sum(L) .div. 2))
                counter = 1
                incr = 2
            else if (m_s == alpha) then
                allocate(possible_values(sum(L) .div. 2))
                counter = 2
                incr = 2
            else
                allocate(possible_values(sum(L)))
                counter = 1
                incr = 1
            end if

            ! Here we merge the values from splitted_orbitals sortedly into
            ! possible_values
            i = 1
            do while (any(counter <= L))
                curr_value = huge(curr_value)
                do iGAS = spaces(1), spaces(2)
                    if (counter(iGAS) <= GAS_spec%GAS_size(iGAS)) then
                        if (GAS_spec%get_orb_idx(counter(iGAS), iGAS) < curr_value) then
                            curr_value = GAS_spec%get_orb_idx(counter(iGAS), iGAS)
                            iGAS_min_val = iGAS
                        end if
                    end if
                end do

                counter(iGAS_min_val) = counter(iGAS_min_val) + incr
                possible_values(i) = curr_value
                i = i + 1
            end do

            if (present(add_particles)) then
                possible_holes = complement(possible_values, union(det_I, add_particles))
            else
                possible_holes = complement(possible_values, det_I)
            end if
        end block
    end function



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
            exc%val(2) = 0
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
            nJ = 0
            ilutJ = 0
        end if

        pgen = pgen_particle * pgen_hole
        @:ASSERT(all(nJ == 0) .neqv. 0.0_dp < pgen .and. pgen <= 1.0_dp)
    end subroutine

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
        type(SpinProj_t) :: m_s_1, m_s_2
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
        integer :: i, elec

        integer :: elecs(2), sym_product, ispn, sum_ml, tgt(2)
        integer :: ms, nJBase(nel)
        logical :: tExchange

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

    DEBUG_IMPURE subroutine draw_from_cum_list(c_sum, idx, pgen)
        real(dp), intent(in) :: c_sum(:)
        integer, intent(out) :: idx
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'draw_from_cum_list'

        @:ASSERT((c_sum(size(c_sum)) .isclose. 0.0_dp) &
            .or. (c_sum(size(c_sum)) .isclose. 1.0_dp))
        @:ASSERT(is_sorted(c_sum))

        ! there might not be such an excitation
        if (c_sum(size(c_sum)) > 0) then
            ! find the index of the target hole in the cumulative list
            idx = binary_search_first_ge(c_sum, genrand_real2_dSFMT())
            @:ASSERT(1 <= idx .and. idx <= size(c_sum))

            ! adjust pgen with the probability for picking tgt from the cumulative list
            if (idx == 1) then
                pgen = c_sum(1)
            else
                pgen = c_sum(idx) - c_sum(idx - 1)
            end if
        else
            idx = 0
            pgen = 0.0_dp
        end if
    end subroutine

#:for excitation_t in ExcitationTypes
    function get_cumulative_list_${excitation_t}$(det_I, incomplete_exc, possible_holes) result(cSum)
        integer, intent(in) :: det_I(:)
        type(${excitation_t}$), intent(in) :: incomplete_exc
        integer, intent(in) :: possible_holes(:)
        real(dp) :: cSum(size(possible_holes))
        character(*), parameter :: this_routine = 'get_cumulative_list_${excitation_t}$'

        real(dp) :: previous
        type(${excitation_t}$) :: exc
        integer :: i

        @:ASSERT(get_last_tgt(exc) == UNKNOWN)
        exc = incomplete_exc

        ! build the cumulative list of matrix elements <src|H|tgt>
        previous = 0.0_dp
        do i = 1, size(possible_holes)
            call set_last_tgt(exc, possible_holes(i))
            cSum(i) = abs(sltcnd_excit(det_I, exc, .false.)) + previous
            previous = cSum(i)
        end do

        ! Normalize
        if (near_zero(cSum(size(cSum)))) then
            cSum(:) = 0.0_dp
        else
            cSum(:) = cSum(:) / cSum(size(cSum))
        end if
    end function get_cumulative_list_${excitation_t}$
#:endfor



    DEBUG_IMPURE function get_available_singles(GAS_spec, det_I) result(singles_exc_list)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: singles_exc_list(:, :)
        character(*), parameter :: this_routine = 'get_available_singles'

        integer, allocatable :: possible_holes(:)
        integer, allocatable :: tmp_buffer(:, :)
        integer :: i, j, src1, tgt1
        type(SpinProj_t) :: m_s_1
        type(buffer_int_2D_t) :: buffer


        @:ASSERT(GAS_spec%contains(det_I))

        call buffer%init(size(det_I), 1000_int64)

        do i = 1, size(det_I)
            src1 = det_I(i)
            m_s_1 = calc_spin_raw(src1)
            possible_holes = get_possible_holes(&
                        GAS_spec, det_I, add_holes=det_I(i:i), &
                        excess=-m_s_1)
            do j = 1, size(possible_holes)
                tgt1 = possible_holes(j)
                call buffer%add_val(excite(det_I, SingleExc_t(src1, tgt1)))
            end do
        end do
        call buffer%dump(singles_exc_list)

        @:sort(integer, singles_exc_list, rank=2, along=2, comp=lex_leq)
    end function

    DEBUG_IMPURE function get_available_doubles(GAS_spec, det_I) result(doubles_exc_list)
        type(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: doubles_exc_list(:, :)
        character(*), parameter :: this_routine = 'get_available_doubles'

        integer, allocatable :: first_pick_possible_holes(:), second_pick_possible_holes(:), deleted(:)
        integer :: i, j, k, l, src1, src2, tgt1, tgt2
        type(SpinProj_t) :: m_s_1
        type(buffer_int_2D_t) :: buffer

        @:ASSERT(GAS_spec%contains(det_I))

        call buffer%init(size(det_I), 1000_int64)

        do i = 1, size(det_I)
            do j = i + 1, size(det_I)
                src1 = det_I(i)
                src2 = det_I(j)
                deleted = det_I([i, j])
                first_pick_possible_holes = get_possible_holes(GAS_spec, det_I, &
                                        add_holes=deleted, excess=-sum(calc_spin_raw(deleted)), &
                                        n_total=2)
                @:ASSERT(disjoint(first_pick_possible_holes, det_I))
                do k = 1, size(first_pick_possible_holes)
                    tgt1 = first_pick_possible_holes(k)
                    m_s_1 = calc_spin_raw(tgt1)
                    @:ASSERT(any(m_s_1 == [alpha, beta]))

                    second_pick_possible_holes = get_possible_holes(&
                            GAS_spec, det_I, add_holes=deleted, &
                            add_particles=[tgt1], &
                            n_total=1, excess=m_s_1 - sum(calc_spin_raw(deleted)))

                    @:ASSERT(disjoint(second_pick_possible_holes, [tgt1]))
                    @:ASSERT(disjoint(second_pick_possible_holes, det_I))

                    do l = 1, size(second_pick_possible_holes)
                        tgt2 = second_pick_possible_holes(l)
                        call buffer%add_val(excite(det_I, DoubleExc_t(src1, tgt1, src2, tgt2)))
                    end do
                end do
            end do
        end do

        call buffer%dump(doubles_exc_list)

        @:sort(integer, doubles_exc_list, rank=2, along=2, comp=lex_leq)

        remove_double_appearances : block
            integer, allocatable :: tmp_buffer(:, :)
            allocate(tmp_buffer, mold=doubles_exc_list)
            j = 1
            tmp_buffer(:, j) = doubles_exc_list(:, 1)
            do i = 2, size(doubles_exc_list, 2)
                if (any(doubles_exc_list(:, i - 1) /= doubles_exc_list(:, i))) then
                    j = j + 1
                    tmp_buffer(:, j) = doubles_exc_list(:, i)
                end if
            end do
            doubles_exc_list = tmp_buffer(:, : j)
        end block remove_double_appearances
    end function

    subroutine gen_all_excits(nI, n_excits, det_list)
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)

        integer, allocatable :: singles(:, :), doubles(:, :)
        integer :: i, j, k

        singles = get_available_singles(GAS_specification, nI)
        doubles = get_available_doubles(GAS_specification, nI)

        n_excits = size(singles, 2) + size(doubles, 2)
        allocate(det_list(0:niftot, n_excits))
        j = 1
        do i = 1, size(singles, 2)
            call EncodeBitDet(singles(:, i), det_list(:, j))
            j = j + 1
        end do

        do i = 1, size(doubles, 2)
            call EncodeBitDet(doubles(:, i), det_list(:, j))
            j = j + 1
        end do

        call sort(det_list, ilut_lt, ilut_gt)
    end subroutine gen_all_excits
end module gasci_general
