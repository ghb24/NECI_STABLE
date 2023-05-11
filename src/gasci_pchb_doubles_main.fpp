#include "macros.h"
#:include "macros.fpph"

module gasci_pchb_doubles_main
    use constants, only: stdout
    use SystemData, only: tUHF
    use util_mod, only: EnumBase_t, stop_all
    use excitation_generators, only: DoubleExcitationGenerator_t
    use fortran_strings, only: to_upper
    use fortran_strings, only: split
    use gasci, only: GASSpec_t
    use gasci_pchb_doubles_spatorb_fastweighted, only: GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t
    use gasci_pchb_doubles_spinorb_fullyweighted, only: GAS_PCHB_DoublesSpinorbFullyWeightedExcGenerator_t
    use gasci_pchb_doubles_spinorb_fastweighted, only: GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t
    use gasci_pchb_doubles_select_particles, only: PCHB_ParticleSelection_t, &
        PCHB_ParticleSelection_vals_t
    better_implicit_none

    private

    public :: PCHB_HoleSelection_t, &
        PCHB_DoublesOptions_t, allocate_and_init, &
        PCHB_DoublesOptions_vals_t, doubles_options_vals
    ! reexpose stuff from doubles particle selection
    public :: PCHB_ParticleSelection_t

    type, extends(EnumBase_t) :: PCHB_HoleSelection_t
    end type

    type :: PCHB_HoleSelection_vals_t
        type(PCHB_HoleSelection_t) :: &
            FAST_FAST = PCHB_HoleSelection_t(1), &
            FULL_FULL = PCHB_HoleSelection_t(2)
        contains
            procedure, nopass :: from_str => select_holes_from_keyword
    end type

    type :: PCHB_DoublesOptions_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(PCHB_HoleSelection_t) :: hole_selection
        logical, allocatable :: spin_orb_resolved
            !! We want to distinguish between "not yet a value"
            !!  and the actual value.
    end type

    type :: PCHB_DoublesOptions_vals_t
        type(PCHB_ParticleSelection_vals_t) :: particle_selection = PCHB_ParticleSelection_vals_t()
        type(PCHB_HoleSelection_vals_t) :: hole_selection = PCHB_HoleSelection_vals_t()
        contains
            procedure, nopass :: from_str => select_doubles_option_from_keyword
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

        if (options%spin_orb_resolved) then
            if (options%hole_selection == doubles_options_vals%hole_selection%FAST_FAST) then
                allocate(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t :: generator)
            else if (options%hole_selection == doubles_options_vals%hole_selection%FULL_FULL) then
                allocate(GAS_PCHB_DoublesSpinorbFullyWeightedExcGenerator_t :: generator)
            end if
        else
            if (tUHF) call stop_all(this_routine, "spatial-orbital-resolved PCHB generator not compatible with UHF FCIDUMP.")
            if (options%hole_selection == doubles_options_vals%hole_selection%FAST_FAST) then
                allocate(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t :: generator)
            else if (options%hole_selection == doubles_options_vals%hole_selection%FULL_FULL) then
                call stop_all(this_routine, "PCHB spatorb-resolved fully weighted hole selection not implemented.")
            end if
        end if

        select type(generator)
        type is(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t)
            call generator%init(&
                GAS_spec, use_lookup, use_lookup, options%particle_selection)
        type is(GAS_PCHB_DoublesSpinorbFullyWeightedExcGenerator_t)
            call generator%init(&
                GAS_spec, use_lookup, use_lookup, options%particle_selection)
        type is(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t)
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
        case('FAST-FAST')
            res = doubles_options_vals%hole_selection%FAST_FAST
        case('FULL-FULL')
            res = doubles_options_vals%hole_selection%FULL_FULL
        case default
            call stop_all(this_routine, trim(w)//" not a valid hole selection for GAS PCHB.")
        end select
    end function

    pure function select_doubles_option_from_keyword(w) result(res)
        !! Parse a given keyword into the possible weighting schemes
        character(*), intent(in) :: w
        type(PCHB_DoublesOptions_t) :: res
        associate(tokens => split(w, ':'))
            res%particle_selection = doubles_options_vals%particle_selection%from_str(tokens(1)%str)
            res%hole_selection = doubles_options_vals%hole_selection%from_str(tokens(2)%str)
        end associate
    end function

end module gasci_pchb_doubles_main
