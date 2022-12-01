#:include "macros.fpph"
#include "macros.h"
module gasci_pchb_doubles_main
    use util_mod, only: EnumBase_t
    use gasci_pchb_doubles_select_particles, only: PCHB_ParticleSelection_t, PCHB_particle_selections
    better_implicit_none

    private

    public :: PCHB_HoleSelection_t, possible_PCHB_hole_selection, &
        PCHB_DoublesOptions_t
    ! Reexpose stuff from particle_selection
    public :: PCHB_ParticleSelection_t, PCHB_particle_selections

    type, extends(EnumBase_t) :: PCHB_HoleSelection_t
    end type

    type :: possible_PCHB_HoleSelection_t
        type(PCHB_HoleSelection_t) :: &
            RHF_FAST_WEIGHTED = PCHB_HoleSelection_t(1), &
            UHF_FAST_WEIGHTED = PCHB_HoleSelection_t(2), &
            UHF_FULLY_WEIGHTED = PCHB_HoleSelection_t(3)
    end type

    type(possible_PCHB_HoleSelection_t), parameter :: possible_PCHB_hole_selection = possible_PCHB_HoleSelection_t()

    type :: PCHB_DoublesOptions_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(PCHB_HoleSelection_t) :: hole_selection
    end type

contains

end module gasci_pchb_doubles_main
