#include "macros.h"
#:include "macros.fpph"

module gasci_pchb_doubles_main
    use util_mod, only: EnumBase_t, stop_all
    use excitation_generators, only: DoubleExcitationGenerator_t
    use fortran_strings, only: to_upper
    use gasci, only: GASSpec_t
    use gasci_pchb_doubles_rhf_fastweighted, only: GAS_doubles_RHF_PCHB_ExcGenerator_t
    use gasci_pchb_doubles_UHF_fullyweighted, only: GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t
    use gasci_pchb_doubles_select_particles, only: PCHB_ParticleSelection_t,  PCHB_particle_selection_vals, &
        PCHB_ParticleSelection_vals_t
    better_implicit_none

    private

    public :: PCHB_HoleSelection_t, possible_PCHB_hole_selection, &
        PCHB_DoublesOptions_t, allocate_and_init, &
        PCHB_DoublesOptions_vals_t, doubles_options_vals
    ! reexpose stuff from doubles particle selection
    public :: PCHB_ParticleSelection_t, PCHB_particle_selection_vals

    type, extends(EnumBase_t) :: PCHB_HoleSelection_t
    end type

    type :: PCHB_HoleSelection_vals_t
        type(PCHB_HoleSelection_t) :: &
            RHF_FAST_WEIGHTED = PCHB_HoleSelection_t(1), &
            RHF_FULLY_WEIGHTED = PCHB_HoleSelection_t(2), &
            UHF_FAST_WEIGHTED = PCHB_HoleSelection_t(3), &
            UHF_FULLY_WEIGHTED = PCHB_HoleSelection_t(4)
        contains
            procedure, nopass :: from_str => select_holes_from_keyword
    end type

    type(PCHB_HoleSelection_vals_t), parameter :: possible_PCHB_hole_selection = PCHB_HoleSelection_vals_t()

    type :: PCHB_DoublesOptions_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(PCHB_HoleSelection_t) :: hole_selection
    end type

    type :: PCHB_DoublesOptions_vals_t
        type(PCHB_ParticleSelection_vals_t) :: particle_selection = PCHB_ParticleSelection_vals_t()
        type(PCHB_HoleSelection_vals_t) :: hole_selection = PCHB_HoleSelection_vals_t()
    end type

    type(PCHB_DoublesOptions_vals_t), parameter :: doubles_options_vals = PCHB_DoublesOptions_vals_t()

contains

    subroutine allocate_and_init(GAS_spec, options, use_lookup, generator)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(PCHB_DoublesOptions_t), intent(in) :: options
        logical, intent(in) :: use_lookup
            !! Use the supergroup lookup
        class(DoubleExcitationGenerator_t), allocatable, intent(inout) :: generator
        routine_name("gasci_pchb_doubles_main::allocate_and_init")

        if (options%hole_selection == possible_PCHB_hole_selection%RHF_FAST_WEIGHTED) then
            allocate(GAS_doubles_RHF_PCHB_ExcGenerator_t :: generator)
        else if (options%hole_selection == possible_PCHB_hole_selection%RHF_FULLY_WEIGHTED) then
            call stop_all(this_routine, "PCHB RHF fully weighted hole selection not yet implemented.")
        else if (options%hole_selection == possible_PCHB_hole_selection%UHF_FAST_WEIGHTED) then
            call stop_all(this_routine, "PCHB UHF fast weighted hole selection not yet implemented.")
        else if (options%hole_selection == possible_PCHB_hole_selection%UHF_FULLY_WEIGHTED) then
            allocate(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t :: generator)
        else
            call stop_all(this_routine, "Invalid hole selection algorithm for PCHB doubles.")
        end if

        select type(generator)
        type is(GAS_doubles_RHF_PCHB_ExcGenerator_t)
            call generator%init(&
                GAS_spec, use_lookup, use_lookup, options%particle_selection)
        type is(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t)
            call generator%init(&
                GAS_spec, use_lookup, use_lookup, options%particle_selection)
        class default
            call stop_all(this_routine, "Error. Should never be here.")
        end select
    end subroutine

    pure function select_holes_from_keyword(w) result(res)
        !! Parse a given keyword into the possible weighting schemes
        character(*), intent(in) :: w
        type(PCHB_HoleSelection_t) :: res
        routine_name("gasci_pchb_doubles_main::select_holes_from_keyword")
        select case(to_upper(w))
        case('RHF-FAST-WEIGHTED')
            res = possible_PCHB_hole_selection%RHF_FAST_WEIGHTED
        case('UHF-FULLY-WEIGHTED')
            res = possible_PCHB_hole_selection%UHF_FULLY_WEIGHTED
        case default
            call stop_all(this_routine, trim(w)//" not a valid hole selection for GAS PCHB.")
        end select
    end function

end module gasci_pchb_doubles_main
