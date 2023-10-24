#include "macros.h"
#:include "../macros.fpph"

module guga_gas_pchb_doubles_main
    use constants, only: stdout
    use SystemData, only: tUHF
    use util_mod, only: EnumBase_t, stop_all
    use fortran_strings, only: to_upper, split
    use gasci, only: GASSpec_t

    use guga_base_class, only: DoublesGUGABase_t
    use guga_gas_pchb_doubles_select_holes, only: GAS_PCHB_FastFast_t, &
        PCHB_HoleAlgorithm_t, PCHB_HoleAlgorithm_vals_t, &
        WeightingChoice_t, WeightingChoice_vals_t
    use guga_gas_pchb_doubles_select_particles, only: &
        PCHB_ParticleSelection_t, &
        PCHB_ParticleSelection_vals_t

    better_implicit_none

    private

    public :: PCHB_DoublesOptions_t, PCHB_DoublesOptions_vals_t, allocate_and_init

    type :: PCHB_DoublesOptions_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(PCHB_HoleAlgorithm_t)  :: hole_selection
        type(WeightingChoice_t)  :: weighting
    contains
        procedure :: is_valid
        procedure :: to_str
    end type

    type :: PCHB_DoublesOptions_vals_t
        type(PCHB_ParticleSelection_vals_t) :: &
            particle_selection = PCHB_ParticleSelection_vals_t()
        type(PCHB_HoleAlgorithm_vals_t) :: &
            hole_selection = PCHB_HoleAlgorithm_vals_t()
        type(WeightingChoice_vals_t) :: &
            weighting = WeightingChoice_vals_t()
    end type

    type(PCHB_DoublesOptions_vals_t), parameter :: &
        doubles_options_vals = PCHB_DoublesOptions_vals_t()

contains

    subroutine allocate_and_init(GAS_spec, options, use_lookup, create_lookup, generator)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(PCHB_DoublesOptions_t), intent(in) :: options
        logical, intent(in) :: use_lookup, create_lookup
            !! Use the supergroup lookup
        class(DoublesGUGABase_t), allocatable, intent(inout) :: generator
        routine_name("gasci_pchb_doubles_main::allocate_and_init")

        if (options%hole_selection == doubles_options_vals%hole_selection%FAST_FAST) then
            allocate(GAS_PCHB_FastFast_t :: generator)
        end if

        select type(generator)
        type is(GAS_PCHB_FastFast_t)
            call generator%init(&
                GAS_spec, options%particle_selection, options%weighting, use_lookup, create_lookup)
        class default
            call stop_all(this_routine, "Error. Should never be here.")
        end select
    end subroutine

    logical elemental function is_valid(this)
        class(PCHB_DoublesOptions_t), intent(in) :: this
        @:unused_var(this)
        is_valid = .true.
    end function

    pure function to_str(this) result(res)
        class(PCHB_DoublesOptions_t), intent(in) :: this
        character(:), allocatable :: res
        res = this%particle_selection%to_str()//':'//this%hole_selection%to_str()//' '//this%weighting%to_str()
    end function

end module
