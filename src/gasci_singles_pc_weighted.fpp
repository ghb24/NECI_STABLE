#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_singles_pc_weighted
    use constants, only: dp, int64, stdout, n_int, bits_n_int, maxExcit
    use util_mod, only: operator(.div.), stop_all, EnumBase_t, near_zero
    use bit_rep_data, only: NIfTot, nIfD
    use bit_reps, only: decode_bit_det
    use SymExcitDataMod, only: ScratchSize
    use SystemData, only: nEl, nBasis, G1
    use excitation_generators, only: SingleExcitationGenerator_t
    use FciMCData, only: excit_gen_store_type
    use dSFMT_interface, only: genrand_real2_dSFMT
    use aliasSampling, only: AliasSampler_2D_t
    use get_excit, only: make_single
    use excitation_types, only: SingleExc_t
    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use orb_idx_mod, only: calc_spin_raw, operator(==)
    use OneEInts, only: GetTMatEl
    use UMatCache, only: GTID, get_umat_el
    use CDF_sampling_mod, only: CDF_Sampler_t

    use matrix_util, only: print_matrix

    better_implicit_none
    private
    public :: PC_singles_weighted_t, possible_PC_singles_weighted, PC_weighted_singles, &
        Base_PC_SinglesLocalised_t, PC_UniformSingles_t, PC_SinglesHOnly_t, &
        PC_HAndGTerm_t, PC_HAndGTermBothAbs_t, from_keyword

    type, extends(EnumBase_t) :: PC_singles_weighted_t
    end type

    type :: possible_PC_singles_weighted_t
        type(PC_singles_weighted_t) :: &
            UNIFORM = PC_singles_weighted_t(1), &
            H_ONLY = PC_singles_weighted_t(2), &
                !! \( |h_{I, A}| \)
            H_AND_G_TERM = PC_singles_weighted_t(3), &
                !! \( | h_{I, A} + \sum_{R} g_{I, A, R, R} - g_{I, R, R, A} | \)
            H_AND_G_TERM_BOTH_ABS = PC_singles_weighted_t(4)
                !! \( | h_{I, A} | + | \sum_{R} g_{I, A, R, R} - g_{I, R, R, A} | \)
    end type

    type(possible_PC_singles_weighted_t), parameter :: &
        possible_PC_singles_weighted = possible_PC_singles_weighted_t()

    type(PC_singles_weighted_t) :: PC_weighted_singles = possible_PC_singles_weighted%UNIFORM

    type, abstract, extends(SingleExcitationGenerator_t) :: Base_PC_SinglesLocalised_t
        type(AliasSampler_2D_t) :: sampler
            !! p(A | I, i_sg)
            !! The probability of picking the hole A after having picked particle I
            !! in the supergroup i_sg.
        real(dp), allocatable :: weights(:, :, :)
            !! The weights w_{A, I, i_sg} for the excitation of I -> A.
            !! They are made independent of the determinant by various approximations, which
            !! are implemented in the inheriting classes by overwriting `get_weight`.
            !! For example setting \( w_{A, I, i_sg} = | h_{I, A} + \sum_{R} g_{I, A, R, R} - g_{I, R, R, A} | \)
            !! where \(R\) runs over all orbitals instead of only the occupied.
        class(GASSpec_t), allocatable :: GAS_spec
            !! The GAS specification
        type(SuperGroupIndexer_t), pointer :: indexer => null()
            !! The Supergroup indexer.
            !! This is only a pointer because components cannot be targets
            !! otherwise. :-(
        logical, public :: use_lookup = .false.
            !! Use a lookup for the supergroup index in global_det_data.
        logical, public :: create_lookup = .false.
    contains
        private
        procedure, public :: init
        procedure, public :: gen_exc
        procedure, public :: get_pgen
        procedure, public :: finalize
        procedure(get_weight_t), nopass, deferred :: get_weight
    end type

    abstract interface
        real(dp) elemental function get_weight_t(exc)
            import :: SingleExc_t, dp
            implicit none
            type(SingleExc_t), intent(in) :: exc
        end function
    end interface

    type, extends(Base_PC_SinglesLocalised_t) :: PC_UniformSingles_t
    contains
        private
        procedure, nopass :: get_weight => Uniform_get_weight
    end type

    type, extends(Base_PC_SinglesLocalised_t) :: PC_SinglesHOnly_t
    contains
        private
        procedure, nopass :: get_weight => PC_hOnly_get_weight
    end type

    type, extends(Base_PC_SinglesLocalised_t) :: PC_HAndGTerm_t
    contains
        private
        procedure, nopass :: get_weight => PC_HAndGTerm_get_weight
    end type

    type, extends(Base_PC_SinglesLocalised_t) :: PC_HAndGTermBothAbs_t
    contains
        private
        procedure, nopass :: get_weight => PC_HAndGTermBothAbs_get_weight
    end type


contains

    pure function from_keyword(w) result(res)
        character(*), intent(in) :: w
        type(PC_singles_weighted_t) :: res
        routine_name("from_keyword")
        select case(w)
        case('UNIFORM')
            res = possible_PC_singles_weighted%UNIFORM
        case('H-ONLY')
            res = possible_PC_singles_weighted%H_ONLY
        case('H-AND-G-TERM')
            res = possible_PC_singles_weighted%H_AND_G_TERM
        case('H-AND-G-TERM-BOTH-ABS')
            res = possible_PC_singles_weighted%H_AND_G_TERM_BOTH_ABS
        case default
            call Stop_All(this_routine, trim(w)//" not a valid PC-WEIGHTED singles generator")
        end select
    end function

    subroutine init(this, GAS_spec, use_lookup, create_lookup)
        class(GASSpec_t), intent(in) :: GAS_spec
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this
        logical, intent(in) :: use_lookup, create_lookup
        routine_name("init")
        integer, allocatable :: supergroups(:, :)
        integer :: n_supergroups, nBI

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'PC particles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'PC particles is using the supergroup lookup'

        nBI = this%GAS_spec%n_spin_orbs()
        supergroups = this%indexer%get_supergroups()
        n_supergroups = size(supergroups, 2)

        call this%sampler%shared_alloc([nBi, n_supergroups], nBI, 'PC_singles')
        allocate(this%weights(nBI, nBI, n_supergroups), source=0._dp)
        block
            integer :: i_sg, src, tgt
            type(SingleExc_t) :: exc
            do i_sg = 1, n_supergroups
                do src = 1, nBi
                    do tgt = 1, nBi
                        exc = SingleExc_t(src, tgt)
                        if (this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) &
                                .and. calc_spin_raw(src) == calc_spin_raw(tgt) &
                                .and. src /= tgt &
                                .and. symmetry_allowed(exc) &
                        ) then
                            this%weights(tgt, src, i_sg) = this%get_weight(exc)
                        end if
                    end do
                    call this%sampler%setup_entry(src, i_sg, this%weights(:, src, i_sg))
                end do
            end do
        end block

    contains
            ! For single excitations it is simple
            logical pure function symmetry_allowed(exc)
                use SymExcitDataMod, only: SpinOrbSymLabel
                type(SingleExc_t), intent(in) :: exc
                symmetry_allowed = SpinOrbSymLabel(exc%val(1)) == SpinOrbSymLabel(exc%val(2))
            end function
    end subroutine

    subroutine gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        integer :: elec, src, tgt, i_sg
        real(dp) :: p_src, p_tgt

        @:unused_var(exFlag, part_type)
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        ic = 1

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if



        block
            integer :: unoccupied(this%GAS_spec%n_spin_orbs() - nEl), idx
            type(CDF_Sampler_t) :: tgt_sampler, src_sampler
            real(dp), allocatable :: weights(:, :)

            call decode_bit_det(unoccupied, not(ilutI))

            weights = this%weights(unoccupied, nI, i_sg)

            src_sampler = CDF_Sampler_t(sum(weights, dim=1))

            call src_sampler%sample(elec, p_src)
            if (elec == 0) then
                call make_invalid()
                return
            end if
            src = nI(elec)


            tgt_sampler = CDF_Sampler_t(this%weights(unoccupied, src, i_sg))

            call tgt_sampler%sample(idx, p_tgt)
            if (idx == 0) then
                call make_invalid()
                return
            end if
            tgt = unoccupied(idx)
        end block

        pGen = p_src * p_tgt
        call make_single(nI, nJ, elec, tgt, ex, tParity)
        ilutJ = ilutI
        clr_orb(ilutJ, src)
        set_orb(ilutJ, tgt)
        contains

            subroutine make_invalid()
                nJ(1) = 0
                ilutJ = 0_n_int
            end subroutine
    end subroutine gen_exc

    subroutine finalize(this)
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this

        call this%sampler%finalize()
        deallocate(this%indexer)
        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine

    function get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(p_gen)
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        debug_function_name("get_pgen")
        real(dp) :: p_gen
        integer :: i_sg

        @:ASSERT(ic == 1)
        @:unused_var(ClassCount2, ClassCountUnocc2, ilutI)
        i_sg = this%indexer%idx_nI(nI)
        associate (src => ex(1, 1), tgt => ex(2, 1))
            p_gen = 1._dp / real(nEl, dp) * this%sampler%get_prob(src, i_sg, tgt)
        end associate
    end function


    elemental function Uniform_get_weight(exc) result(w)
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: w
        @:unused_var(exc)
        w = 1._dp
    end function


    elemental function PC_hOnly_get_weight(exc) result(w)
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: w
        w = abs(h(exc%val(1), exc%val(2)))
    end function


    elemental function PC_HAndGTerm_get_weight(exc) result(w)
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: w
        integer :: R
        real(dp) :: two_el_term
        associate(I => exc%val(1), A => exc%val(2))
            two_el_term = 0._dp
            do R = 1, nBasis
                two_el_term = two_el_term + g(I, A, R, R) - G(I, R, R, A)
            end do
            w = abs(h(I, A) + two_el_term)
        end associate
    end function


    elemental function PC_HAndGTermBothAbs_get_weight(exc) result(w)
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: w
        integer :: R
        real(dp) :: two_el_term
        associate(I => exc%val(1), A => exc%val(2))
            two_el_term = 0._dp
            do R = 1, nBasis
                two_el_term = two_el_term + abs(g(I, A, R, R) - G(I, R, R, A))
            end do
            w = abs(h(I, A)) + two_el_term
        end associate
    end function


    real(dp) elemental function h(I, A)
        !! Return the 1el integral \( h_{I, A) \)
        !!
        !! I and A are **spin** indices.
        !! Follows the definition of the purple book.
        integer, intent(in) :: I, A
        h = GetTMATEl(I, A)
    end function


    real(dp) elemental function g(I, A, J, B)
        integer, intent(in) :: I, A, J, B
        !! Return the 2el integral \( g_{I, A, J, B} \)
        !!
        !! I, A, J, and B are **spin** indices.
        !! Order follows the definition of the purple book.
        if (calc_spin_raw(I) == calc_spin_raw(A) .and. calc_spin_raw(J) == calc_spin_raw(B)) then
            g = get_umat_el(gtID(I), gtID(A), gtID(J), gtID(B))
        else
            g = 0._dp
        end if
    end function


end module gasci_singles_pc_weighted
