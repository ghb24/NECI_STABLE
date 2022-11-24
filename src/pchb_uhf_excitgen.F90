#include "macros.h"

module pchb_uhf_excit
    !! precomputed heat bath excitation generator, open shell ("UHF")
    ! @jph
    use excitation_generators, only: ExcitationGenerator_t, SingleExcitationGenerator_t
    ! @jph not completely confident about importing this one vs writing it myself
    ! use gasci_pchb_rhf, only: PCHB_ParticleSelection_t
    better_implicit_none
    private
    public :: PCHB_UHF_FCI_excit_generator_t

    ! @jph
    ! FCI_PCHB_particle_selection, &
    !     FCI_PCHB_singles, possible_PCHB_singles

    type, extends(ExcitationGenerator_t) :: PCHB_UHF_FCI_excit_generator_t
        private
        type(GAS_doubles_PCHB_UHF_ExcGenerator_t) :: doubles_generator
        class(SingleExcitationGenerator_t), allocatable :: singles_generator
    contains
    ! deferred from parent class
    ! gen_exc
    ! get_pgen
    ! gen_all_excits
    ! finalize
        private
        procedure, public :: init
        procedure, public :: finalize
        procedure, public :: gen_exc
        procedure, public :: get_pgen
        procedure, public :: gen_all_excits
    end type PCHB_UHF_FCI_excit_generator_t

    ! type(PCHB_ParticleSelection_t) :: FCI_PCHB_UHF_particle_selection = PCHB_particle_selections

    type, extends(EnumBase_t) :: PCHB_UHF_used_singles_t
    end type PCHB_UHF_used_singles_t

    type :: possible_PCHB_UHF_singles_t
        type(PCHB_used_singles_t) :: &
            ON_FLY_HEAT_BATH = PCHB_UHF_used_singles_t(1), &
            UNIFORM = PCHB_UHF_used_singles_t(2)
    end type

    type(possible_PCHB_UHF_singles_t), parameter :: possible_PCHB_UHF_singles = possible_PCHB_UHF_singles

    type(PCHB_UHF_used_singles_t) :: FCI_PCHB_UHF_singles = possible_PCHB_UHF_singles%UNIFORM

contains

end module pchb_uhf_excit