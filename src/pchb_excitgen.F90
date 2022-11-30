#include "macros.h"
module pchb_excitgen
    use constants, only: n_int, dp, maxExcit, stdout
    use SystemData, only: nel, nBasis
    use fortran_strings, only: to_upper
    use util_mod, only: operator(.div.), EnumBase_t, stop_all
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t

    use gasci, only: GASSpec_t, LocalGASSpec_t
    use excitation_generators, only: ClassicAbInitExcitationGenerator_t
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t
    use gasci_singles_pc_weighted, only: PC_WeightedSinglesOptions_t, possible_PC_singles_drawing, &
        possible_PC_singles_weighting, Base_PC_Weighted_t, do_allocation, print_options
    use gasci_pchb_rhf, only: GAS_doubles_RHF_PCHB_ExcGenerator_t
    use gasci_pc_select_particles, only: PCHB_ParticleSelection_t, PCHB_particle_selections
    better_implicit_none

    private

    public :: PCHB_FCI_excit_generator_t, &
        FCI_PCHB_options_t, FCI_PCHB_options, &
        possible_PCHB_singles, singles_from_keyword

    type, extends(ClassicAbInitExcitationGenerator_t) :: PCHB_FCI_excit_generator_t
    contains
        private
        procedure, public :: init
    end type

    type, extends(EnumBase_t) :: PCHB_used_singles_t
    end type

    type :: possible_PCHB_singles_t
        type(PCHB_used_singles_t) :: &
            ON_FLY_HEAT_BATH = PCHB_used_singles_t(1), &
            UNIFORM = PCHB_used_singles_t(2), &
            PC_WEIGHTED = PCHB_used_singles_t(3)
    end type

    type(possible_PCHB_singles_t), parameter :: possible_PCHB_singles = possible_PCHB_singles_t()

    type :: FCI_PCHB_options_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(PCHB_used_singles_t) :: singles
        logical :: UHF
            !! Do a spin-projection resolved calculation.
        type(PC_WeightedSinglesOptions_t) :: PC_singles_options = PC_WeightedSinglesOptions_t(&
            possible_PC_singles_weighting%UNDEFINED, possible_PC_singles_drawing%UNDEFINED)
            !! Only relevant if `singles == possible_PCHB_singles%PC_WEIGHTED`
    end type

    type(FCI_PCHB_options_t) :: FCI_PCHB_options = FCI_PCHB_options_t(&
        PCHB_particle_selections%PC_WEIGHTED_APPROX, &
        possible_PCHB_singles%PC_WEIGHTED, &
        UHF=.false., &
        PC_singles_options=PC_WeightedSinglesOptions_t(&
            possible_PC_singles_weighting%H_AND_G_TERM_BOTH_ABS, &
            possible_PC_singles_drawing%APPROX &
        ) &
    )


contains

    pure function singles_from_keyword(w) result(res)
        !! Parse a given keyword into the possible weighting schemes
        character(*), intent(in) :: w
        type(PCHB_used_singles_t) :: res
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
        if (options%UHF) then
            call stop_all(this_routine, 'UHF PCHB not yet implemented :(')
        else
            allocate(GAS_doubles_RHF_PCHB_ExcGenerator_t :: this%doubles_generator)
            select type(generator => this%doubles_generator)
            type is(GAS_doubles_RHF_PCHB_ExcGenerator_t)
                call generator%init(&
                CAS_spec(n_el=nEl, n_spat_orbs=nBasis .div. 2), &
                use_lookup=.false., create_lookup=.false., &
                PCHB_particle_selection=options%particle_selection)
            end select
        end if

        ! luckily the singles generators don't require initialization.
        if (options%singles == possible_PCHB_singles%ON_FLY_HEAT_BATH) then
            allocate(WeightedSingles_t :: this%singles_generator)
        else if (options%singles == possible_PCHB_singles%UNIFORM) then
            allocate(UniformSingles_t :: this%singles_generator)
        else if (options%singles == possible_PCHB_singles%PC_WEIGHTED) then
            call print_options(options%PC_singles_options, stdout)
            call do_allocation(this%singles_generator, options%PC_singles_options%drawing)
            select type(generator => this%singles_generator)
            class is(Base_PC_Weighted_t)
                call generator%init(CAS_spec(n_el=nEl, n_spat_orbs=nBasis .div. 2), &
                                    options%PC_singles_options%weighting, &
                                    use_lookup=.false., create_lookup=.false.)
            end select
        else
            call stop_all(this_routine, "Invalid PCHB_singles in FCI PCHB init.")
        end if
    contains

        type(LocalGASSpec_t) pure function CAS_spec(n_el, n_spat_orbs)
            integer, intent(in) :: n_el, n_spat_orbs
            integer :: i
            ! It does not matter if we use local or cumulative GAS
            ! constraints
            CAS_spec = LocalGASSpec_t(&
                n_min=[n_el], n_max=[n_el], spat_GAS_orbs=[(1, i = 1, n_spat_orbs)])
        end function
    end subroutine

end module pchb_excitgen
