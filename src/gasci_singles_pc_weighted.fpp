#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_singles_pc_weighted
    use constants, only: dp, int64, stdout, n_int, bits_n_int, maxExcit
    use util_mod, only: operator(.div.), stop_all, EnumBase_t
    use bit_rep_data, only: NIfTot, nIfD
    use SymExcitDataMod, only: ScratchSize
    use SystemData, only: nEl
    use excitation_generators, only: SingleExcitationGenerator_t
    use FciMCData, only: excit_gen_store_type
    use dSFMT_interface, only: genrand_real2_dSFMT
    use aliasSampling, only: AliasSampler_2D_t
    use get_excit, only: make_single
    use excitation_types, only: SingleExc_t
    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use orb_idx_mod, only: calc_spin_raw, operator(==)
    better_implicit_none
    private
    public :: PC_UniformSingles_t, PC_singles_weighted_t, possible_PC_singles_weighted

    type, extends(EnumBase_t) :: PC_singles_weighted_t
    end type

    type :: possible_PC_singles_weighted_t
        type(PC_singles_weighted_t) :: &
            UNIFORM = PC_singles_weighted_t(1)
    end type

    type(possible_PC_singles_weighted_t), parameter :: &
        possible_PC_singles_weighted = possible_PC_singles_weighted_t()

    type, abstract, extends(SingleExcitationGenerator_t) :: Base_PC_SinglesLocalised_t
        type(AliasSampler_2D_t) :: sampler
            !! p(A | I, i_sg)
            !! The probability of picking the hole A after having picked particle I
            !! in the supergroup i_sg.

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
        procedure, public :: init => PC_SinglesLocalised_init
        procedure, public :: gen_exc => PC_SinglesLocalised_gen_exc
        procedure, public :: get_pgen => PC_SinglesLocalised_get_pgen
        procedure, public :: finalize => PC_SinglesLocalised_finalize
        procedure(get_weight_t), deferred :: get_weight
    end type

    abstract interface
        real(dp) pure function get_weight_t(this, exc)
            import :: Base_PC_SinglesLocalised_t, SingleExc_t, dp
            implicit none
            class(Base_PC_SinglesLocalised_t), intent(in) :: this
            type(SingleExc_t), intent(in) :: exc
        end function
    end interface

    type, extends(Base_PC_SinglesLocalised_t) :: PC_UniformSingles_t
    contains
        private
        procedure :: get_weight => Uniform_get_weight
    end type


contains

    subroutine PC_SinglesLocalised_init(this, GAS_spec, use_lookup, create_lookup)
        class(GASSpec_t), intent(in) :: GAS_spec
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this
        logical, intent(in) :: use_lookup, create_lookup
        routine_name("PC_SinglesLocalised_init")
        real(dp), allocatable :: weights(:)
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
        allocate(weights(nBI))
        block
            integer :: i_sg, I, A
            type(SingleExc_t) :: exc
            do i_sg = 1, n_supergroups
                do I = 1, nBi
                    do A = 1, nBi
                        exc = SingleExc_t(I, A)
                        if (this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) &
                                .and. calc_spin_raw(I) == calc_spin_raw(A) &
                                .and. I /= A &
                                .and. symmetry_allowed(exc) &
                        ) then
                            weights(A) = this%get_weight(exc)
                        else
                            weights(A) = 0._dp
                        end if
                    end do
                    call this%sampler%setup_entry(I, i_sg, weights)
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

    subroutine PC_SinglesLocalised_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
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

        real(dp) :: p_src, p_tgt
        integer :: elec, src, tgt, i_sg

        @:unused_var(exFlag, part_type)
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        ic = 1

        elec = int(genrand_real2_dSFMT() * nel) + 1
        src = nI(elec)
        p_src = 1._dp / real(nEl, dp)

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        call this%sampler%sample(src, i_sg, tgt, p_tgt)

        pGen = p_src * p_tgt
        if (IsOcc(ilutI, tgt)) then
            nJ(1) = 0
            ilutJ = 0_n_int
        else
            call make_single(nI, nJ, elec, tgt, ex, tParity)
            ilutJ = ilutI
            clr_orb(ilutJ, src)
            set_orb(ilutJ, tgt)
        end if
    end subroutine PC_SinglesLocalised_gen_exc

    subroutine PC_SinglesLocalised_finalize(this)
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this

        call this%sampler%finalize()
        deallocate(this%indexer)
        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine

    function PC_SinglesLocalised_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(p_gen)
        class(Base_PC_SinglesLocalised_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        debug_function_name("PC_SinglesLocalised_get_pgen")
        real(dp) :: p_gen
        integer :: i_sg

        @:ASSERT(ic == 1)
        @:unused_var(ClassCount2, ClassCountUnocc2, ilutI)
        i_sg = this%indexer%idx_nI(nI)
        associate (src => ex(1, 1), tgt => ex(2, 1))
            p_gen = 1._dp / real(nEl, dp) * this%sampler%get_prob(src, i_sg, tgt)
        end associate
    end function


    pure function Uniform_get_weight(this, exc) result(w)
        class(PC_UniformSingles_t), intent(in) :: this
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: w
        @:unused_var(this, exc)
        w = 1._dp
    end function
end module gasci_singles_pc_weighted
