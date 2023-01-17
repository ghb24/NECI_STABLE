#include "macros.h"
module pchb_excitgen
    use constants, only: n_int, dp, maxExcit, stdout
    use SystemData, only: nel, nBasis
    use fortran_strings, only: to_upper
    use util_mod, only: operator(.div.), EnumBase_t, stop_all, operator(.implies.)
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t
    use SystemData, only: tUHF

    use gasci, only: GASSpec_t, LocalGASSpec_t
    use excitation_generators, only: ClassicAbInitExcitationGenerator_t, &
        SingleExcitationGenerator_t
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t
    use gasci_singles_pc_weighted, only: &
        PC_WeightedSinglesOptions_t,  PC_WeightedSinglesOptions_vals_t, &
        PC_singles_drawing_vals, &
        PC_singles_weighting_vals, PC_Weighted_t, do_allocation, print_options
    use gasci_pchb_doubles_main, only: PCHB_DoublesOptions_t, PCHB_DoublesOptions_vals_t, &
        doubles_allocate_and_init => allocate_and_init, &
        possible_PCHB_hole_selection, PCHB_particle_selection_vals
    use gasci_pchb_doubles_select_particles, only: PCHB_ParticleSelection_t, PCHB_particle_selection_vals
    better_implicit_none

    private

    public :: PCHB_FCI_excit_generator_t, &
        FCI_PCHB_options_t, FCI_PCHB_Options_vals_t, FCI_PCHB_options_vals, &
        FCI_PCHB_options
    public :: FCI_PCHB_SinglesOptions_t
    ! Reexpose
    public :: PCHB_DoublesOptions_t


    type, extends(ClassicAbInitExcitationGenerator_t) :: PCHB_FCI_excit_generator_t
    contains
        private
        procedure, public :: init
    end type

    type, extends(EnumBase_t) :: FCI_PCHB_singles_algorithm_t
    end type

    type :: FCI_PCHB_singles_algorithm_vals_t
        type(FCI_PCHB_singles_algorithm_t) :: &
            ON_FLY_HEAT_BATH = FCI_PCHB_singles_algorithm_t(1), &
            UNIFORM = FCI_PCHB_singles_algorithm_t(2), &
            PC_WEIGHTED = FCI_PCHB_singles_algorithm_t(3)
        contains
            procedure, nopass :: from_str => singles_from_keyword
    end type

    type(FCI_PCHB_singles_algorithm_vals_t), parameter :: &
        possible_PCHB_singles = FCI_PCHB_singles_algorithm_vals_t()

    type :: FCI_PCHB_SinglesOptions_vals_t
        type(FCI_PCHB_singles_algorithm_vals_t) :: &
            algorithm = FCI_PCHB_singles_algorithm_vals_t()
        type(PC_WeightedSinglesOptions_vals_t) :: &
            PC_weighted = PC_WeightedSinglesOptions_vals_t()
    end type

    type(FCI_PCHB_SinglesOptions_vals_t), parameter :: &
        FCI_PCHB_singles_options_vals = FCI_PCHB_SinglesOptions_vals_t()

    type :: FCI_PCHB_SinglesOptions_t
        type(FCI_PCHB_singles_algorithm_t) :: algorithm
        type(PC_WeightedSinglesOptions_t) :: PC_weighted = PC_WeightedSinglesOptions_t(&
            FCI_PCHB_singles_options_vals%PC_weighted%weighting%UNDEFINED, &
            FCI_PCHB_singles_options_vals%PC_weighted%drawing%UNDEFINED)
    end type

    type :: FCI_PCHB_Options_t
        type(FCI_PCHB_SinglesOptions_t) :: singles
        type(PCHB_DoublesOptions_t) :: doubles
        logical :: spinorb_resolved
            !! Do a spin-projection-resolved calculation.
    contains
        procedure :: assert_validity
    end type

    type :: FCI_PCHB_Options_vals_t
        type(FCI_PCHB_SinglesOptions_vals_t) :: singles = FCI_PCHB_SinglesOptions_vals_t()
        type(PCHB_DoublesOptions_vals_t) :: doubles = PCHB_DoublesOptions_vals_t()
    end type

    type(FCI_PCHB_Options_vals_t), parameter :: FCI_PCHB_options_vals = FCI_PCHB_Options_vals_t()

    type(FCI_PCHB_options_t) :: FCI_PCHB_options = FCI_PCHB_options_t(&
        FCI_PCHB_SinglesOptions_t(&
            FCI_PCHB_options_vals%singles%algorithm%PC_WEIGHTED, &
            PC_WeightedSinglesOptions_t(&
                FCI_PCHB_options_vals%singles%PC_weighted%weighting%H_AND_G_TERM_BOTH_ABS, &
                FCI_PCHB_options_vals%singles%PC_weighted%drawing%FAST_WEIGHTED &
            ) &
        ), &
        PCHB_DoublesOptions_t( &
            FCI_PCHB_options_vals%doubles%particle_selection%FAST_WEIGHTED, &
            FCI_PCHB_options_vals%doubles%hole_selection%SPATORB_FAST_WEIGHTED &
        ), &
        spinorb_resolved=.false. &
    )


contains

    pure function singles_from_keyword(w) result(res)
        !! Parse a given keyword into the possible weighting schemes
        character(*), intent(in) :: w
        type(FCI_PCHB_singles_algorithm_t) :: res
        routine_name("from_keyword")
        select case(to_upper(w))
        case('UNIFORM')
            res = possible_PCHB_singles%UNIFORM
        case('ON-THE-FLY-HEAT-BATH')
            res = possible_PCHB_singles%ON_FLY_HEAT_BATH
        case('PC-WEIGHTED')
            res = possible_PCHB_singles%PC_WEIGHTED
        case default
            call stop_all(this_routine, trim(w)//" not a valid singles generator for GAS PCHB.")
        end select
    end function


    subroutine init(this, options)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this
        type(FCI_PCHB_options_t), intent(in) :: options
        character(*), parameter :: this_routine = 'pchb_excitgen::init'

        ! CAS is implemented as a special case of GAS with only one GAS space.
        ! Since a GAS specification with one GAS space is trivially disconnected, there
        ! is no point to use the lookup.
        ! @jph TODO implement UHF
        call doubles_allocate_and_init(&
                CAS_spec(n_el=nEl, n_spat_orbs=nBasis .div. 2), options%doubles, &
                .false., this%doubles_generator)
        call singles_allocate_and_init(options%singles, this%singles_generator)
    end subroutine


    subroutine singles_allocate_and_init(options, generator)
        type(FCI_PCHB_SinglesOptions_t), intent(in) :: options
        class(SingleExcitationGenerator_t), allocatable, intent(inout) :: generator
        routine_name("pchb_excitgen::singles_allocate_and_init")

        if (allocated(generator)) then
            call generator%finalize()
            deallocate(generator)
        end if

        ! luckily many of the singles generators don't require initialization.
        if (options%algorithm == possible_PCHB_singles%ON_FLY_HEAT_BATH) then
            allocate(WeightedSingles_t :: generator)
        else if (options%algorithm == possible_PCHB_singles%UNIFORM) then
            allocate(UniformSingles_t :: generator)
        else if (options%algorithm == possible_PCHB_singles%PC_WEIGHTED) then
            call print_options(options%PC_weighted, stdout)
            call do_allocation(generator, options%PC_weighted%drawing)
            select type(generator)
            class is(PC_Weighted_t)
                call generator%init(CAS_spec(n_el=nEl, n_spat_orbs=nBasis .div. 2), &
                                    options%PC_weighted%weighting, &
                                    use_lookup=.false., create_lookup=.false.)
            end select
        else
            call stop_all(this_routine, "Invalid PCHB_singles in FCI PCHB init.")
        end if
    end subroutine


    type(LocalGASSpec_t) pure function CAS_spec(n_el, n_spat_orbs)
        integer, intent(in) :: n_el, n_spat_orbs
        integer :: i
        ! It does not matter if we use local or cumulative GAS
        ! constraints
        CAS_spec = LocalGASSpec_t(&
            n_min=[n_el], n_max=[n_el], spat_GAS_orbs=[(1, i = 1, n_spat_orbs)])
    end function


    subroutine assert_validity(this)
        class(FCI_PCHB_options_t), intent(in) :: this
        routine_name("assert_validity")

        if (.not. (this%singles%algorithm == possible_PCHB_singles%PC_WEIGHTED &
               .implies. (this%singles%PC_weighted%weighting /= PC_singles_weighting_vals%UNDEFINED &
                        .and. this%singles%PC_weighted%drawing /= PC_singles_drawing_vals%UNDEFINED))) then
            call stop_all(this_routine, "PC-WEIGHTED singles require valid PC_weighted options.")
        end if

        if (.not. (tUHF .implies. this%spinorb_resolved)) then
            call stop_all(this_routine, "UHF requires spin-resolved PCHB")
        end if

        if (.not. (this%spinorb_resolved &
                    .implies. (this%doubles%hole_selection == possible_PCHB_hole_selection%SPINORB_FAST_WEIGHTED &
                                .or. this%doubles%hole_selection == possible_PCHB_hole_selection%SPINORB_FULLY_WEIGHTED))) then
            call stop_all(this_routine, "Spin-resolved excitation generation requires spin-resolved hole generation.")
        end if

    end subroutine

end module pchb_excitgen
