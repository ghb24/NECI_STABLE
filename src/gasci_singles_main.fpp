#:include "macros.fpph"
#include "macros.h"
module gasci_singles_main
    use util_mod, only: EnumBase_t, stop_all, operator(.implies.)
    use SystemData, only: nEl
    use excitation_generators, only: SingleExcitationGenerator_t
    use constants, only: n_int, dp, maxExcit, bits_n_int, stdout
    use dSFMT_interface, only: genrand_real2_dSFMT
    use SymExcitDataMod, only: ScratchSize
    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_on_the_fly_heat_bath, only: GAS_singles_heat_bath_ExcGen_t
    use get_excit, only: make_single
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use excitation_types, only: Excite_1_t, spin_allowed
    use FciMCData, only: excit_gen_store_type, GAS_PCHB_init_time
    use fortran_strings, only: to_upper
    use bit_rep_data, only: NIfTot, nIfD
    use bit_reps, only: decode_bit_det

    use gasci_singles_pc_weighted, only: PC_Weighted_t, &
        do_allocation, print_options, &
        PC_WeightedSinglesOptions_t, PC_WeightedSinglesOptions_vals_t
    better_implicit_none

    private
    public :: GAS_PCHB_SinglesAlgorithm_t, &
        GAS_singles_PC_uniform_ExcGenerator_t, &
        GAS_singles_heat_bath_ExcGen_t, allocate_and_init, &
        GAS_PCHB_SinglesOptions_t, GAS_PCHB_SinglesOptions_vals_t, GAS_PCHB_singles_options_vals
    ! Reexpose the stuff from gasci_singles_pc_weighted
    public :: PC_WeightedSinglesOptions_t, PC_Weighted_t, do_allocation, &
        print_options


    type, extends(EnumBase_t) :: GAS_PCHB_SinglesAlgorithm_t
    end type

    type :: GAS_PCHB_SinglesAlgorithm_vals_t
        type(GAS_PCHB_SinglesAlgorithm_t) :: &
            ON_FLY_HEAT_BATH = GAS_PCHB_SinglesAlgorithm_t(1), &
            BITMASK_UNIFORM = GAS_PCHB_SinglesAlgorithm_t(2), &
            PC_WEIGHTED = GAS_PCHB_SinglesAlgorithm_t(3)
        contains
            procedure, nopass :: from_str => singles_from_keyword
    end type

    type(GAS_PCHB_SinglesAlgorithm_vals_t), parameter :: GAS_used_singles_vals = GAS_PCHB_SinglesAlgorithm_vals_t()

    type :: GAS_PCHB_SinglesOptions_vals_t
        type(GAS_PCHB_SinglesAlgorithm_vals_t) :: algorithm = GAS_PCHB_SinglesAlgorithm_vals_t()
        type(PC_WeightedSinglesOptions_vals_t) :: PC_weighted = PC_WeightedSinglesOptions_vals_t()
    end type

    type(GAS_PCHB_SinglesOptions_vals_t), parameter :: GAS_PCHB_singles_options_vals = GAS_PCHB_SinglesOptions_vals_t()

    type :: GAS_PCHB_SinglesOptions_t
        type(GAS_PCHB_SinglesAlgorithm_t) :: algorithm
        type(PC_WeightedSinglesOptions_t) :: PC_weighted = PC_WeightedSinglesOptions_t(&
            GAS_PCHB_singles_options_vals%PC_weighted%weighting%UNDEFINED, GAS_PCHB_singles_options_vals%PC_weighted%drawing%UNDEFINED)
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

contains


    subroutine allocate_and_init(GAS_spec, options, use_lookup, generator)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_PCHB_SinglesOptions_t), intent(in) :: options
        logical, intent(in) :: use_lookup
            !! Use the supergroup lookup
        class(SingleExcitationGenerator_t), allocatable, intent(inout) :: generator
        routine_name("gasci_singles_main::allocate_and_init")

        if (allocated(generator)) then
            call generator%finalize()
            deallocate(generator)
        end if
        if (options%algorithm == GAS_used_singles_vals%BITMASK_UNIFORM) then
            write(stdout, *) 'GAS precomputed singles activated'
            allocate(GAS_singles_PC_uniform_ExcGenerator_t :: generator)
            select type(generator => generator)
            type is(GAS_singles_PC_uniform_ExcGenerator_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, use_lookup, create_lookup=.false.)
            end select
        else if (options%algorithm == GAS_used_singles_vals%ON_FLY_HEAT_BATH) then
            write(stdout, *) 'GAS heat bath on the fly singles activated'
            allocate(generator, source=GAS_singles_heat_bath_ExcGen_t(GAS_spec))
        else if (options%algorithm == GAS_used_singles_vals%PC_WEIGHTED) then
            call print_options(options%PC_weighted, stdout)
            call do_allocation(generator, options%PC_weighted%drawing)
            select type(generator)
            class is(PC_Weighted_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, options%PC_weighted%weighting, &
                                    use_lookup, create_lookup=.false.)
            end select
        else
            call stop_all(this_routine, "Invalid choice for singles.")
        end if
    end subroutine



    pure function singles_from_keyword(w) result(res)
        !! Parse a given keyword into the possible weighting schemes
        character(*), intent(in) :: w
        type(GAS_PCHB_SinglesAlgorithm_t) :: res
        routine_name("singles_from_keyword")
        select case(to_upper(w))
        case('UNIFORM')
            res = GAS_used_singles_vals%BITMASK_UNIFORM
        case('ON-THE-FLY-HEAT-BATH')
            res = GAS_used_singles_vals%ON_FLY_HEAT_BATH
        case('PC-WEIGHTED')
            res = GAS_used_singles_vals%PC_WEIGHTED
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
                    if (src == tgt) cycle
                    associate(exc => Excite_1_t(src, tgt))
                    if (spin_allowed(exc) &
                            .and. this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) &
                            .and. symmetry_allowed(exc)) then
                        call my_set_orb(this%allowed_holes(:, src, i_sg), tgt)
                    end if
                    end associate
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
                type(Excite_1_t), intent(in) :: exc
                symmetry_allowed = SpinOrbSymLabel(exc%val(1, 1)) == SpinOrbSymLabel(exc%val(2, 1))
            end function
    end subroutine GAS_singles_uniform_init

    subroutine GAS_singles_uniform_finalize(this)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this

        if (allocated(this%allowed_holes)) then
            deallocate(this%allowed_holes, this%indexer)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
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



end module gasci_singles_main
