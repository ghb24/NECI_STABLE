#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module gasci_general
    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    use constants, only: dp, n_int, maxExcit
    use util_mod, only: binary_search_first_ge, &
        near_zero, operator(.isclose.), custom_findloc
    use sort_mod, only: sort
    use DetBitOps, only: EncodeBitDet
    use sets_mod, only: is_sorted, disjoint
    use bit_rep_data, only: NIfTot
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: excit_gen_store_type, pParallel
    use get_excit, only: make_double, make_single
    use SymExcitDataMod, only: ScratchSize
    use excit_gens_int_weighted, only: pick_biased_elecs, get_pgen_pick_biased_elecs
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        get_last_tgt, set_last_tgt, defined, UNKNOWN, &
        excite, ilut_excite
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, &
        operator(-), operator(==), operator(/=), sum
    use util_mod, only: intswap

    use excitation_generators, only: &
        SingleExcitationGenerator_t, DoubleExcitationGenerator_t, ExcitationGenerator_t,  &
        gen_exc_sd, gen_all_excits_sd, get_pgen_sd

    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit

    use gasci, only: GASSpec_t
    use gasci_util, only: &
        get_possible_holes, get_available_singles, get_available_doubles, &
        get_cumulative_list, draw_from_cum_list
    implicit none

    private
    public :: GAS_heat_bath_ExcGenerator_t, GAS_singles_heat_bath_ExcGen_t

    !> The heath bath GAS on-the-fly excitation generator
    type, extends(SingleExcitationGenerator_t) :: GAS_singles_heat_bath_ExcGen_t
        private
        type(GASSpec_t) :: GAS_spec
    contains
        private
        procedure, public :: finalize => GAS_singles_do_nothing

        !> Get the GAS allowed holes for a given determinant and a chosen particle.
        procedure, public :: get_possible_holes => GAS_singles_possible_holes
        procedure, public :: gen_exc => GAS_singles_gen_exc
        procedure, public :: get_pgen => GAS_singles_get_pgen
        procedure, public :: gen_all_excits => GAS_singles_gen_all_excits
    end type

    !> The heath bath GAS on-the-fly excitation generator
    type, extends(DoubleExcitationGenerator_t) :: GAS_doubles_heat_bath_ExcGenerator_t
        private
        type(GASSpec_t) :: GAS_spec
    contains
        procedure, public :: finalize => GAS_doubles_do_nothing

        !> Get the GAS allowed holes for a given determinant and a chosen particle.
        procedure, public :: get_possible_holes => GAS_doubles_possible_holes
        procedure, public :: gen_exc => GAS_doubles_gen_exc
        procedure, public :: get_pgen => GAS_doubles_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_gen_all_excits
    end type


    type, extends(ExcitationGenerator_t) :: GAS_heat_bath_ExcGenerator_t
        private
        type(GAS_singles_heat_bath_ExcGen_t) :: singles_generator
        type(GAS_doubles_heat_bath_ExcGenerator_t) :: doubles_generator
    contains
        procedure, public :: finalize => GAS_heat_bath_finalize
        procedure, public :: gen_exc => GAS_heat_bath_gen_exc
        procedure, public :: get_pgen => GAS_heat_bath_get_pgen
        procedure, public :: gen_all_excits => GAS_heat_bath_gen_all_excits
    end type


    interface GAS_heat_bath_ExcGenerator_t
        module procedure construct_GAS_heat_bath_ExcGenerator_t
    end interface


    interface GAS_singles_heat_bath_ExcGen_t
        module procedure construct_GAS_singles_heat_bath_ExcGen_t
    end interface

contains

    pure function construct_GAS_heat_bath_ExcGenerator_t(GAS_spec) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_heat_bath_ExcGenerator_t) :: res
        res%singles_generator = GAS_singles_heat_bath_ExcGen_t(GAS_spec)
        res%doubles_generator = GAS_doubles_heat_bath_ExcGenerator_t(GAS_spec)
    end function

    pure function construct_GAS_singles_heat_bath_ExcGen_t(GAS_spec) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_singles_heat_bath_ExcGen_t) :: res
        res%GAS_spec = GAS_spec
    end function

    subroutine GAS_singles_do_nothing(this)
        class(GAS_singles_heat_bath_ExcGen_t), intent(inout) :: this
        ! no Initialization => no Finalization
        @:unused_var(this)
    end subroutine

    pure function GAS_singles_possible_holes(this, nI, src) result(unoccupied)
        class(GAS_singles_heat_bath_ExcGen_t), intent(in) :: this
        integer, intent(in) :: nI(:), src
        integer, allocatable :: unoccupied(:)
        unoccupied = get_possible_holes(this%GAS_spec, nI, add_holes=[src], &
                                        n_total=1, excess=-calc_spin_raw(src))
    end function

    !>  @brief
    !>  Generate a single excitation under GAS constraints.
    subroutine GAS_singles_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                   ex, tParity, pGen, hel, store, part_type)
        class(GAS_singles_heat_bath_ExcGen_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nEl), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "GAS_singles_gen_exc"

        type(SingleExc_t) :: exc
        integer, allocatable :: possible_holes(:)
        real(dp) :: pgen_particle, pgen_hole
        real(dp), allocatable :: c_sum(:)
        integer :: i, elec

        @:unused_var(store, exFlag, part_type)
        ic = 1
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0.0_dp)
#endif
        ! Pick any random electron
        elec = int(genrand_real2_dSFMT() * nEl) + 1
        exc%val(1) = nI(elec)
        pgen_particle = 1.0_dp / real(nEl, kind=dp)

        ! Get a hole with the same spin projection
        ! while fullfilling GAS-constraints.
        possible_holes = this%get_possible_holes(nI, exc%val(1))

        if (size(possible_holes) == 0) then
            call zero_result()
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt \in possible_holes
        c_sum = get_cumulative_list(nI, exc, possible_holes)
        call draw_from_cum_list(c_sum, i, pgen_hole)
        @:ASSERT(i == 0 .neqv. (0.0_dp < pgen_hole .and. pgen_hole <= 1.0_dp))
        if (i /= 0) then
            exc%val(2) = possible_holes(i)
            call make_single(nI, nJ, elec, exc%val(2), ex, tParity)
            ilutJ = ilut_excite(ilutI, exc)
        else
            call zero_result()
        end if

        pgen = pgen_particle * pgen_hole
        @:ASSERT(nJ(1) == 0 .neqv. 0.0_dp < pgen .and. pgen <= 1.0_dp)
        contains

            subroutine zero_result()
                ex(:, 1) = exc%val(:)
                nJ(1) = 0
                ilutJ = 0_n_int
            end subroutine zero_result
    end subroutine

    function GAS_singles_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_singles_heat_bath_ExcGen_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_singles_get_pgen'

        real(dp) :: pgen_particle, pgen_hole
        real(dp), allocatable :: c_sum(:)
        integer, allocatable :: possible_holes(:)
        integer :: idx
        @:unused_var(ClassCount2, ClassCountUnocc2, ilutI)
        @:ASSERT(ic == 1)

        pgen_particle = 1.0_dp / real(nEl, kind=dp)
        possible_holes = this%get_possible_holes(nI, ex(1, 1))
        @:ASSERT(any(ex(2, 1) == possible_holes))
        idx = custom_findloc(possible_holes, ex(2, 1))
        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt \in possible_holes
        c_sum = get_cumulative_list(nI, SingleExc_t(src=ex(1, 1)), possible_holes)
        if (idx == 1) then
            pgen_hole = c_sum(1)
        else
            pgen_hole = c_sum(idx) - c_sum(idx - 1)
        end if
        pgen = pgen_particle * pgen_hole
    end function

    subroutine GAS_singles_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_singles_heat_bath_ExcGen_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        integer, allocatable :: singles(:, :)
        integer :: i

        singles = get_available_singles(this%GAS_spec, nI)
        n_excits = size(singles, 2)
        allocate(det_list(0:niftot, n_excits))
        do i = 1, size(singles, 2)
            call EncodeBitDet(singles(:, i), det_list(:, i))
        end do
    end subroutine


    subroutine GAS_doubles_do_nothing(this)
        class(GAS_doubles_heat_bath_ExcGenerator_t), intent(inout) :: this
        ! no Initialization => no Finalization
        @:unused_var(this)
    end subroutine

    pure function GAS_doubles_possible_holes(this, nI, src, tgt) result(unoccupied)
        class(GAS_doubles_heat_bath_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(:), src(2)
        integer, intent(in), optional :: tgt
        integer, allocatable :: unoccupied(:)
        type(SpinProj_t) :: excess
        character(*), parameter :: this_routine = "GAS_doubles_possible_holes"

        if (present(tgt)) then
            excess = calc_spin_raw(tgt) - sum(calc_spin_raw(src))
            unoccupied = get_possible_holes(&
                this%GAS_spec, nI, add_holes=src, add_particles=[tgt], &
                n_total=1, excess=excess)
        else
            unoccupied = get_possible_holes(&
                this%GAS_spec, nI, add_holes=src, n_total=2, excess=-sum(calc_spin_raw(src)))
        end if
        @:pure_ASSERT(disjoint(unoccupied, nI))
    end function

    !>  @brief
    !>  Generate a single excitation under GAS constraints.
    subroutine GAS_doubles_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                   ex, tParity, pGen, hel, store, part_type)
        class(GAS_doubles_heat_bath_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nEl), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'GAS_doubles_gen_exc'

        type(DoubleExc_t) :: exc
        integer, allocatable :: possible_holes(:)
        integer :: deleted(2)
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

        @:unused_var(exFlag, store, part_type)
        ic = 2
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0.0_dp)
#endif

        call pick_biased_elecs(nI, elecs, exc%val(1, :), &
                               sym_product, ispn, sum_ml, pgen_particles)
        @:ASSERT(exc%val(1, 1) /= exc%val(2, 1))

        deleted = exc%val(1, :)

        ! Get possible holes for the first particle, while fullfilling and spin- and GAS-constraints.
        ! and knowing that a second particle will be created afterwards!
        possible_holes = this%get_possible_holes(nI, deleted)

        if (size(possible_holes) == 0) then
            pgen = pgen_particles
            call zeroResult()
            return
        end if

        r = genrand_real2_dSFMT()
        ! Pick randomly one hole with arbitrary spin
        exc%val(2, 1) = possible_holes(int(r * real(size(possible_holes), kind=dp)) + 1)
        pgen_first_pick = 1.0_dp / real(size(possible_holes), dp)

        ! Pick second hole.
        ! The total spin projection of the created particles has to add up
        ! to the total spin projection of the deleted particles.
        ! Get possible holes for the second particle,
        ! while fullfilling GAS- and Spin-constraints.
        possible_holes = this%get_possible_holes(nI, deleted, exc%val(2, 1))
        @:ASSERT(disjoint(possible_holes, exc%val(2:2, 1)))

        if (size(possible_holes) == 0) then
            pgen = pgen_particles * pgen_first_pick
            call zeroResult()
            return
        end if

        ! build the cumulative list of matrix elements <src|H|tgt>
        ! with tgt2 in possible_holes
        c_sum = get_cumulative_list(nI, exc, possible_holes)
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
            type(DoubleExc_t) :: reverted_exc
            src1 = exc%val(1, 1); src2 = exc%val(1, 2);
            tgt1 = exc%val(2, 1); tgt2 = exc%val(2, 2);
            @:ASSERT(src1 /= tgt2 .and. src2 /= tgt2, src1, tgt2, src2)

            possible_holes = this%get_possible_holes(nI, deleted, tgt2)

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
            c_sum = get_cumulative_list(nI, reverted_exc, possible_holes)
            if (i == 1) then
                pgen_second_pick(2) = c_sum(1)
            else
                pgen_second_pick(2) = (c_sum(i) - c_sum(i - 1))
            end if

            if (i /= 0) then
                call make_double(nI, nJ, elecs(1), elecs(2), tgt1, tgt2, ex, tParity)
                ilutJ = ilut_excite(ilutI, exc)
            else
                ilutJ = 0
            end if
        end block

        pgen = pgen_particles * pgen_first_pick * sum(pgen_second_pick)

        @:ASSERT(0.0_dp < pgen .and. pgen <= 1.0_dp)
        @:ASSERT(all(ex(:, :2) /= UNKNOWN))

        contains

            subroutine zeroResult()
                integer :: src_copy(2)

                src_copy(:) = exc%val(1, :)
                call sort(src_copy)
                ex(1, :2) = src_copy
                ex(2, :2) = exc%val(2, :)
                nJ(1) = 0
                ilutJ = 0_n_int
            end subroutine zeroResult
    end subroutine

    function GAS_doubles_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_doubles_heat_bath_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_doubles_get_pgen'

        real(dp) :: pgen_particles, &
            ! These are arrays, because the pgen might be different if we
            ! pick AB or BA. Let i, j denote the particles and a, b the holes.
            ! pgen_particles == p (i j)
            ! pgen_first_pick == p(a | i j) == p(b | i j)
            ! pgen_second_pick == [p(b | a i j), p(a | b i j)]
            ! pgen == p_double * p (i j) * p(a | i j) * sum([p(b | a i j), p(a | b i j)])
            !      == p_double * p (i j) * p(b | i j) * sum([p(b | a i j), p(a | b i j)])
            pgen_first_pick, pgen_second_pick(2)
        type(DoubleExc_t) :: exc, reverted_exc
        real(dp), allocatable :: c_sum(:)
        integer, allocatable :: possible_holes(:)
        integer :: i
        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 2)

        exc = DoubleExc_t(ex)

        pgen_particles = get_pgen_pick_biased_elecs(&
                calc_spin_raw(exc%val(1, 1)) == calc_spin_raw(exc%val(1, 2)), &
                pParallel, par_elec_pairs, AB_elec_pairs)

        possible_holes = this%get_possible_holes(nI, exc%val(1, :))
        pgen_first_pick = 1.0_dp / real(size(possible_holes), dp)

        possible_holes = this%get_possible_holes(nI, exc%val(1, :), exc%val(2, 1))
        i = binary_search_first_ge(possible_holes, exc%val(2, 2))
        c_sum = get_cumulative_list(nI, exc, possible_holes)
        if (i == 1) then
            pgen_second_pick(1) = c_sum(1)
        else
            pgen_second_pick(1) = (c_sum(i) - c_sum(i - 1))
        end if

        reverted_exc = exc
        call intswap(reverted_exc%val(2, 1), reverted_exc%val(2, 2))

        possible_holes = this%get_possible_holes(nI, reverted_exc%val(1, :), reverted_exc%val(2, 1))
        i = binary_search_first_ge(possible_holes, reverted_exc%val(2, 2))
        c_sum = get_cumulative_list(nI, reverted_exc, possible_holes)
        if (i == 1) then
            pgen_second_pick(2) = c_sum(1)
        else
            pgen_second_pick(2) = (c_sum(i) - c_sum(i - 1))
        end if

        pgen = pgen_particles * pgen_first_pick * sum(pgen_second_pick)
    end function

    subroutine GAS_doubles_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_doubles_heat_bath_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        integer, allocatable :: doubles(:, :)
        integer :: i

        doubles = get_available_doubles(this%GAS_spec, nI)
        n_excits = size(doubles, 2)
        allocate(det_list(0:niftot, n_excits))
        do i = 1, size(doubles, 2)
            call EncodeBitDet(doubles(:, i), det_list(:, i))
        end do
    end subroutine


    subroutine GAS_heat_bath_finalize(this)
        class(GAS_heat_bath_ExcGenerator_t), intent(inout) :: this
        call this%doubles_generator%finalize()
        call this%singles_generator%finalize()
    end subroutine


    subroutine GAS_heat_bath_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                ex, tParity, pGen, hel, store, part_type)
        class(GAS_heat_bath_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nEl), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        call gen_exc_sd(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                        ex, tParity, pGen, hel, store, part_type, &
                        this%singles_generator, this%doubles_generator)
    end subroutine


    real(dp) function GAS_heat_bath_get_pgen(&
            this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        class(GAS_heat_bath_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        GAS_heat_bath_get_pgen = get_pgen_sd(&
                        nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2, &
                        this%singles_generator, this%doubles_generator)
    end function


    subroutine GAS_heat_bath_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_heat_bath_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        call gen_all_excits_sd(nI, n_excits, det_list, &
                               this%singles_generator, this%doubles_generator)
    end subroutine


end module gasci_general
