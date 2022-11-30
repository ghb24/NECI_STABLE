#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

! The main idea of the precomputed Heat bath sampling (PCHB) is taken from
!    J. Li, M. Otten, A. A. Holmes, S. Sharma, and C. J. Umrigar, J. Comput. Phys. 149, 214110 (2018).
! and described there.
! The main "ingredient" are precomputed probability distributions p(ab | ij) to draw a, b holes
! when i, j electrons were chosen.
! This requires #{i, j | i < j} probability distributions.
!
! The improved version to use spatial orbital indices to save memory is described in
!    Guther K. et al., J. Chem. Phys. 153, 034107 (2020).
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx.
! This gives #{i, j | i < j} * 3 probability distributions.
! NOTE: This is only relevant for RHF-type calculations (see gasci_pchb_rhf.fpp)
!
! The generalization to GAS spaces is available in a preprint (should be available in JCTC soon as well)
!    https://chemrxiv.org/engage/chemrxiv/article-details/61447e60b1d4a605d589af2e
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and GAS constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in gasci_supergroup_index.f90

module gasci_pchb_general
    !! precomputed heat bath implementation for GASCI. This modules implements
    !! the excitation generator which builds on either gasci_pchb and gasci_pchb_uhf
    !! depending on if we are working in spin or spatial orbitals.

    ! note: these two are used discreetly by macros (bad style)
    ! constants : bits_n_int
    ! orb_idx_mod : operator(==)
    use constants, only: n_int, dp, maxExcit, stdout, int32, bits_n_int
    use orb_idx_mod, only: calc_spin_raw, operator(==)
    use fortran_strings, only: to_upper
    use exc_gen_class_wrappers, only: UniformSingles_t
    use util_mod, only: operator(.implies.), EnumBase_t, stop_all
    use dSFMT_interface, only: genrand_real2_dSFMT
    use SymExcitDataMod, only: ScratchSize
    use excitation_types, only: SingleExc_t
    use FciMCData, only: excit_gen_store_type, GAS_PCHB_init_time
    use excit_gens_int_weighted, only: pick_biased_elecs
    use get_excit, only: make_single
    use timing_neci, only: timer, set_timer, halt_timer

    use SystemData, only: nEl
    use bit_rep_data, only: NIfTot, nIfD
    use bit_reps, only: decode_bit_det
    use sort_mod, only: sort

    use gasci, only: GASSpec_t
    use gasci_general, only: GAS_singles_heat_bath_ExcGen_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pc_select_particles, only: &
        ParticleSelector_t, PC_WeightedParticlesOcc_t, &
        PC_FastWeightedParticles_t, UniformParticles_t, &
        PCHB_particle_selections, PCHB_ParticleSelection_t
    use gasci_singles_pc_weighted, only: PC_SinglesOptions_t, &
        Base_PC_Weighted_t, do_allocation, print_options, &
        possible_PC_singles_weighting, possible_PC_singles_drawing

    use excitation_generators, only: ClassicAbInitExcitationGenerator_t, &
            ExcitationGenerator_t, SingleExcitationGenerator_t


    use gasci_pchb_rhf, only: GAS_doubles_RHF_PCHB_ExcGenerator_t
    better_implicit_none


    private
    public :: GAS_PCHB_ExcGenerator_t, &
        GAS_PCHB_options_t, GAS_PCHB_options, &
        singles_from_keyword, possible_GAS_singles

    type, extends(EnumBase_t) :: GAS_used_singles_t
    end type

    type :: possible_GAS_singles_t
        type(GAS_used_singles_t) :: &
            ON_FLY_HEAT_BATH = GAS_used_singles_t(1), &
            DISCARDING_UNIFORM = GAS_used_singles_t(2), &
            BITMASK_UNIFORM = GAS_used_singles_t(3), &
            PC_WEIGHTED = GAS_used_singles_t(4)
    end type

    type(possible_GAS_singles_t), parameter :: possible_GAS_singles = possible_GAS_singles_t()

    type :: GAS_PCHB_options_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(GAS_used_singles_t) :: singles
        logical :: UHF
            !! Do a spin-projection resolved calculation.
        type(PC_SinglesOptions_t) :: PC_singles_options = PC_SinglesOptions_t(&
            possible_PC_singles_weighting%UNDEFINED, possible_PC_singles_drawing%UNDEFINED)
            !! Only relevant if `singles == possible_PCHB_singles%PC_WEIGHTED`
        logical :: use_lookup= .false., create_lookup= .false.
            !! Use and/or create/manage the supergroup lookup.
    end type

    type(GAS_PCHB_options_t) :: GAS_PCHB_options = GAS_PCHB_options_t( &
        PCHB_particle_selections%PC_WEIGHTED_APPROX, &
        possible_GAS_singles%PC_WEIGHTED, &
        PC_singles_options=PC_SinglesOptions_t(&
            possible_PC_singles_weighting%H_AND_G_TERM_BOTH_ABS, &
            possible_PC_singles_drawing%APPROX &
        ), &
        UHF=.false. &
    )

    type, extends(ClassicAbInitExcitationGenerator_t) :: GAS_PCHB_ExcGenerator_t
    contains
        private
        procedure, public :: init => GAS_PCHB_init
    end type


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


    interface GAS_singles_DiscardingGenerator_t
        module procedure construct_GAS_singles_DiscardingGenerator_t
    end interface



contains


    pure function singles_from_keyword(w) result(res)
        !! Parse a given keyword into the possible weighting schemes
        character(*), intent(in) :: w
        type(GAS_used_singles_t) :: res
        routine_name("singles_from_keyword")
        select case(to_upper(w))
        case('UNIFORM')
            res = possible_GAS_singles%BITMASK_UNIFORM
        case('ON-THE-FLY-HEAT-BATH')
            res = possible_GAS_singles%ON_FLY_HEAT_BATH
        case('DISCARDING-UNIFORM')
            res = possible_GAS_singles%DISCARDING_UNIFORM
        case('PC-WEIGHTED')
            res = possible_GAS_singles%PC_WEIGHTED
        case default
            call stop_all(this_routine, trim(w)//" not a valid singles generator for FCI PCHB.")
        end select
    end function


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


    subroutine GAS_PCHB_init(this, GAS_spec, options)
        !! Initialize the PCHB excitation generator.
        !!
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
            !!  The GAS specifications for the excitation generator.
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_PCHB_options_t), intent(in) :: options
        routine_name("GAS_PCHB_init")

        call set_timer(GAS_PCHB_init_time)

        if (options%singles == possible_GAS_singles%DISCARDING_UNIFORM) then
            write(stdout, *) 'GAS discarding singles activated'
            allocate(this%singles_generator, source=GAS_singles_DiscardingGenerator_t(GAS_spec))
        else if (options%singles == possible_GAS_singles%BITMASK_UNIFORM) then
            write(stdout, *) 'GAS precomputed singles activated'
            allocate(GAS_singles_PC_uniform_ExcGenerator_t :: this%singles_generator)
            select type(generator => this%singles_generator)
            type is(GAS_singles_PC_uniform_ExcGenerator_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, options%use_lookup, create_lookup=.false.)
            end select
        else if (options%singles == possible_GAS_singles%ON_FLY_HEAT_BATH) then
            write(stdout, *) 'GAS heat bath on the fly singles activated'
            allocate(this%singles_generator, source=GAS_singles_heat_bath_ExcGen_t(GAS_spec))
        else if (options%singles == possible_GAS_singles%PC_WEIGHTED) then
            call print_options(options%PC_singles_options, stdout)
            call do_allocation(this%singles_generator, options%PC_singles_options%drawing)
            select type(generator => this%singles_generator)
            class is(Base_PC_Weighted_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, options%PC_singles_options%weighting, &
                                    options%use_lookup, create_lookup=.false.)
            end select
        else
            call stop_all(this_routine, "Invalid choise for singles.")
        end if


        ! @jph at the moment only RHF -- implement UHF
        if (options%UHF) then
            call stop_all('gas init', 'UHF PCHB not yet implemented :(')
        else
            allocate(GAS_doubles_RHF_PCHB_ExcGenerator_t :: this%doubles_generator)
            select type(generator => this%doubles_generator)
            type is(GAS_doubles_RHF_PCHB_ExcGenerator_t)
                call generator%init(&
                    GAS_spec, options%use_lookup, options%create_lookup, &
                    options%particle_selection)
            end select
        end if
        call halt_timer(GAS_PCHB_init_time)
    end subroutine GAS_PCHB_init


end module gasci_pchb_general
