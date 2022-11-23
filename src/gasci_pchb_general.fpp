!! this module covers the parts of the code shared by both the modules
!! gasci_pchb and gasci_pchb_uhf
#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_general
    ! use constants, only: dp, n_int
    ! use excitation_generators, only: SingleExcitationGenerator_t
    ! use gasci, only: GASSpec_t
    ! use gasci_supergroup_index, only: SuperGroupIndexer_t
    ! use util_mod, only: EnumBase_t
    ! @jph remove unused modules (not sure how to determine?)
    use constants, only: n_int, dp, int64, maxExcit, stdout, bits_n_int, int32
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==), operator(/=), alpha, beta
    use exc_gen_class_wrappers, only: UniformSingles_t
    use util_mod, only: fuseIndex, getSpinIndex, near_zero, swap, &
        operator(.div.), operator(.implies.), EnumBase_t, &
        operator(.isclose.), swap, stop_all
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use SymExcitDataMod, only: pDoubNew, ScratchSize
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite
    use sltcnd_mod, only: sltcnd_excit
    use procedure_pointers, only: generate_single_excit_t
    use aliasSampling, only: AliasSampler_3D_t
    use UMatCache, only: gtID, numBasisIndices
    use FciMCData, only: pSingles, excit_gen_store_type, pParallel, &
        projEDet, GAS_PCHB_init_time
    use excit_gens_int_weighted, only: pick_biased_elecs, get_pgen_pick_biased_elecs
    use shared_ragged_array, only: shared_ragged_array_int32_t
    use parallel_neci, only: iProcIndex_intra
    use get_excit, only: make_single
    use timing_neci, only: timer, set_timer, halt_timer

    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    use bit_rep_data, only: NIfTot, nIfD
    use bit_reps, only: decode_bit_det
    use sort_mod, only: sort
    use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt

    use gasci, only: GASSpec_t
    use gasci_general, only: GAS_singles_heat_bath_ExcGen_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pc_select_particles, only: &
        ParticleSelector_t, PC_WeightedParticlesOcc_t, &
        PC_FastWeightedParticles_t, UniformParticles_t

    use display_matrices, only: write_matrix

    use excitation_generators, only: &
            ExcitationGenerator_t, SingleExcitationGenerator_t, &
            DoubleExcitationGenerator_t, gen_exc_sd, get_pgen_sd, gen_all_excits_sd
    better_implicit_none

    ! @jph
    ! I think this will include:
    ! - [x] possible singles
    ! - [x] PC_UNIFORM implementations (?)
    ! - [x] discarding uniform
    ! - [x] test compilation / testing
    ! - [ ] clean public/private
    ! - [ ] abstract doubles generator
    ! - [ ] move GAS_PCHB_ExcGenerator_t but with abstract doubles exc gentor
    ! - [ ] test compilation / testing
    ! - [ ] UHF PCHB unit tests
    ! - [ ] UHF PCHB implementation, largely copied from RHF implementation
    !           (which now inherits pchb_general)
    ! - [ ] do the same with the pchb_excitgen files
    ! - [ ] clean use statements
    ! NOTE pchb hole selection is different between the two

    private
    public :: possible_GAS_singles, GAS_PCHB_singles_generator, &
            PCHB_particle_selections, PCHB_ParticleSelection_t, &
            GAS_used_singles_t, GAS_singles_PC_uniform_ExcGenerator_t, &
            GAS_singles_DiscardingGenerator_t, GAS_PCHB_particle_selection

    !> The precomputed GAS uniform excitation generator
    type, extends(SingleExcitationGenerator_t) :: GAS_singles_PC_uniform_ExcGenerator_t
        private
        ! allowed_holes(:, src, i_sg)
        ! is a bitmask that returns for a given supergroup `i_sg` and `src`
        ! the GAS allowed holes.
        integer(n_int), allocatable :: allowed_holes(:, :, :)
        class(GASSpec_t), allocatable :: GAS_spec
        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        !> Use a lookup for the supergroup index in global_det_data.
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.
    contains
        private
        procedure, public :: init => GAS_singles_uniform_init
        procedure, public :: finalize => GAS_singles_uniform_finalize

        !> Get the GAS allowed holes for a given determinant and a chosen particle.
        procedure, public :: get_possible_holes => GAS_singles_uniform_possible_holes
        procedure, public :: gen_exc => GAS_singles_uniform_gen_exc
        procedure, public :: get_pgen => GAS_singles_uniform_get_pgen
        procedure, public :: gen_all_excits => GAS_singles_uniform_gen_all_excits
    end type


    type, extends(SingleExcitationGenerator_t) :: GAS_singles_DiscardingGenerator_t
        private
        type(UniformSingles_t) :: FCI_singles_generator
        class(GASSpec_t), allocatable :: GAS_spec
    contains
        private
        procedure, public :: finalize => GAS_discarding_singles_finalize
        procedure, public :: gen_exc => GAS_discarding_singles_gen_exc
        procedure, public :: get_pgen => GAS_discarding_singles_get_pgen
        procedure, public :: gen_all_excits => GAS_discarding_singles_gen_all_excits
    end type


    type, extends(EnumBase_t) :: GAS_used_singles_t
    end type

    type :: possible_GAS_singles_t
        type(GAS_used_singles_t) :: &
            ON_FLY_HEAT_BATH = GAS_used_singles_t(1), &
            DISCARDING_UNIFORM = GAS_used_singles_t(2), &
            PC_UNIFORM = GAS_used_singles_t(3)
    end type

    type(possible_GAS_singles_t), parameter :: possible_GAS_singles = possible_GAS_singles_t()

    type(GAS_used_singles_t) :: GAS_PCHB_singles_generator = possible_GAS_singles%PC_UNIFORM

    type, extends(EnumBase_t) :: PCHB_ParticleSelection_t
    end type

    type :: possible_PCHB_ParticleSelection_t
        type(PCHB_ParticleSelection_t) :: &
            UNIFORM = PCHB_ParticleSelection_t(1), &
            PC_WEIGHTED = PCHB_ParticleSelection_t(2), &
            PC_WEIGHTED_APPROX = PCHB_ParticleSelection_t(3)
    end type

    type(possible_PCHB_ParticleSelection_t), parameter :: &
        PCHB_particle_selections = possible_PCHB_ParticleSelection_t()


    interface GAS_singles_DiscardingGenerator_t
        module procedure construct_GAS_singles_DiscardingGenerator_t
    end interface

    type(PCHB_ParticleSelection_t) :: GAS_PCHB_particle_selection = PCHB_particle_selections%PC_WEIGHTED


contains
    ! @jph TODO implementation

    subroutine GAS_singles_uniform_init(this, GAS_spec, use_lookup, create_lookup)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        integer, allocatable :: supergroups(:, :)
        character(*), parameter :: this_routine = 'GAS_singles_uniform_init'

        integer :: i_sg, src, tgt

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        if (create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'GAS singles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        this%use_lookup = use_lookup
        if (use_lookup) write(stdout, *) 'GAS singles is using the supergroup lookup'
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        allocate(this%allowed_holes(0 : nIfD, this%GAS_spec%n_spin_orbs(), size(supergroups, 2)), source=0_n_int)


        ! Find for each supergroup (i_sg)
        ! the allowed holes `tgt` for a given `src`.
        do i_sg = 1, size(supergroups, 2)
            ! Note that the loop cannot take a symmetric shape.
            ! It may be, that (1 -> 5) is allowed, but (5 -> 1) is not.
            do src = 1, this%GAS_spec%n_spin_orbs()
                do tgt = 1, this%GAS_spec%n_spin_orbs()
                    if (src /= tgt &
                            .and. calc_spin_raw(src) == calc_spin_raw(tgt) &
                            .and. this%GAS_spec%is_allowed(SingleExc_t(src, tgt), supergroups(:, i_sg)) &
                            .and. symmetry_allowed(SingleExc_t(src, tgt))) then
                        call my_set_orb(this%allowed_holes(:, src, i_sg), tgt)
                    end if
                end do
            end do
        end do

    contains

            ! Ugly Fortran syntax rules forces us to not use the macro.
            ! Even associate would not help here
            ! https://stackoverflow.com/questions/65734764/non-one-indexed-array-from-associate
            pure subroutine my_set_orb(ilut, orb)
                integer(n_int), intent(inout) :: ilut(0 : nIfD)
                integer, intent(in) :: orb
                set_orb(ilut, orb)
            end subroutine

            ! For single excitations it is simple
            logical pure function symmetry_allowed(exc)
                use SymExcitDataMod, only: SpinOrbSymLabel
                type(SingleExc_t), intent(in) :: exc
                symmetry_allowed = SpinOrbSymLabel(exc%val(1)) == SpinOrbSymLabel(exc%val(2))
            end function
    end subroutine GAS_singles_uniform_init

    subroutine GAS_singles_uniform_finalize(this)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this

        deallocate(this%allowed_holes)
        deallocate(this%indexer)

        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine GAS_singles_uniform_finalize

    !> @brief
    !> For a determinant nI and a spin orbital src return
    !>  the GAS allowed orbitals with the same spin as src which are not occupied in nI.
    function GAS_singles_uniform_possible_holes(this, nI, ilutI, src, use_lookup, store) result(unoccupied)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel), src
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
        logical, intent(in) :: use_lookup
        type(excit_gen_store_type), optional, intent(in) :: store
        integer, allocatable :: unoccupied(:)
        integer(n_int) :: ilut_unoccupied(0 : nIfD)
        integer :: i_sg
        character(*), parameter :: this_routine = 'GAS_PC_possible_holes'

        @:ASSERT(use_lookup .implies. present(store))
        if (use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        ilut_unoccupied = iand(this%allowed_holes(:, src, i_sg), not(ilutI))
        allocate(unoccupied(sum(popcnt(ilut_unoccupied))))
        call decode_bit_det(unoccupied, ilut_unoccupied)
    end function GAS_singles_uniform_possible_holes


    !> @brief
    !> This is the uniform singles excitation generator which uses precomputed indices
    !> to generate only GAS allowed excitations.
    subroutine GAS_singles_uniform_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        integer :: elec, src, tgt
        integer, allocatable :: unoccupied(:)

        @:unused_var(exFlag, part_type)
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        ic = 1

        elec = int(genrand_real2_dSFMT() * nel) + 1
        src = nI(elec)

        unoccupied = this%get_possible_holes(nI, ilutI, src, this%use_lookup, store)

        ! NOTE: this is actually possible for some systems.
        if (size(unoccupied) == 0) then
            pgen = 1._dp / nEl
            nJ(1) = 0
            ilutJ = 0_n_int
            return
        end if

        ! NOTE: The `tgt` could be drawn from `unoccupied` with weights according
        !  to the matrix elements < nI | H |  E_{src}^{tgt} nI >.
        !  Probably not worth it and I am too lazy to implement it now.
        tgt = unoccupied(int(genrand_real2_dSFMT() * size(unoccupied)) + 1)
        pgen = 1._dp / (nEl * size(unoccupied))

        call make_single(nI, nJ, elec, tgt, ex, tParity)
        ilutJ = ilutI
        clr_orb(ilutJ, src)
        set_orb(ilutJ, tgt)
    end subroutine GAS_singles_uniform_gen_exc

    function GAS_singles_uniform_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_PC_get_pgen'

        integer :: src
        integer, allocatable :: unoccupied(:)

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 1)
        src = ex(1, 1)
        unoccupied = this%get_possible_holes(nI, ilutI, src, use_lookup=.false.)
        @:ASSERT(size(unoccupied) > 0)
        pgen = 1._dp / (nEl * size(unoccupied))
    end function GAS_singles_uniform_get_pgen


    subroutine GAS_singles_uniform_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=1)
    end subroutine GAS_singles_uniform_gen_all_excits



    pure function construct_GAS_singles_DiscardingGenerator_t(GAS_spec) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_singles_DiscardingGenerator_t) :: res
        res%GAS_spec = GAS_spec
        res%FCI_singles_generator = UniformSingles_t()
    end function construct_GAS_singles_DiscardingGenerator_t

    subroutine GAS_discarding_singles_finalize(this)
        class(GAS_singles_DiscardingGenerator_t), intent(inout) :: this
        call this%FCI_singles_generator%finalize()
    end subroutine GAS_discarding_singles_finalize

    subroutine GAS_discarding_singles_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(GAS_singles_DiscardingGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'GAS_discarding_singles_gen_exc'

        integer :: src_copy(maxExcit)

        @:ASSERT(this%GAS_spec%contains_conf(nI))

        call this%FCI_singles_generator%gen_exc(&
                    nI, ilutI, nJ, ilutJ, exFlag, ic, &
                    ex, tParity, pGen, hel, store, part_type)
        if (nJ(1) /= 0) then
            if (.not. this%GAS_spec%contains_conf(nJ)) then
                src_copy(:ic) = ex(1, :ic)
                call sort(src_copy)
                ex(1, :ic) = src_copy(:ic)
                ex(2, :ic) = ex(2, :ic)
                nJ(1) = 0
                ilutJ = 0_n_int
            end if
        end if
    end subroutine GAS_discarding_singles_gen_exc

    function GAS_discarding_singles_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_singles_DiscardingGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        pgen = this%FCI_singles_generator%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
    end function GAS_discarding_singles_get_pgen


    subroutine GAS_discarding_singles_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_singles_DiscardingGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=1)
    end subroutine GAS_discarding_singles_gen_all_excits


end module gasci_pchb_general
