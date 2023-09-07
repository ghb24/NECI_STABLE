#include "macros.h"
#:include "../macros.fpph"

module guga_gas_pchb_doubles_main
    use constants, only: stdout
    use SystemData, only: tUHF
    use util_mod, only: EnumBase_t, stop_all
    use fortran_strings, only: to_upper
    use fortran_strings, only: split
    use gasci, only: GASSpec_t

    use guga_base_class, only: DoublesGUGABase_t
    use guga_gas_pchb_doubles_fast_fast, only: GAS_PCHB_FastFast_t
    use guga_gas_pchb_doubles_select_particles, only: &
        PCHB_ParticleSelection_t, &
        PCHB_ParticleSelection_vals_t

    better_implicit_none

    private

    public :: PCHB_DoublesOptions_t, PCHB_DoublesOptions_vals_t, allocate_and_init

    type, extends(EnumBase_t) :: PCHB_HoleSelection_t
        character(9) :: str
    end type

    type :: PCHB_HoleSelection_vals_t
        type(PCHB_HoleSelection_t) :: &
            FAST_FAST = PCHB_HoleSelection_t(1, 'FAST-FAST')
    end type

    type :: PCHB_DoublesOptions_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(PCHB_HoleSelection_t) :: hole_selection
    contains
        procedure :: is_valid
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
        class(DoublesGUGABase_t), allocatable, intent(inout) :: generator
        routine_name("gasci_pchb_doubles_main::allocate_and_init")

        if (options%hole_selection == doubles_options_vals%hole_selection%FAST_FAST) then
            allocate(GAS_PCHB_FastFast_t :: generator)
        end if

        select type(generator)
        type is(GAS_PCHB_FastFast_t)
            call generator%init(&
                GAS_spec, options%particle_selection, use_lookup, use_lookup)
        class default
            call stop_all(this_routine, "Error. Should never be here.")
        end select
    end subroutine

    logical elemental function is_valid(this)
        class(PCHB_DoublesOptions_t), intent(in) :: this
        @:unused_var(this)
        is_valid = .true.
    end function

end module
