#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_uhf
    !! same as gasci_pchb but with spin-orbitals
    !! useful for UHF-format FCIDUMP files

    ! @jph reformulate as notes instead of self-ramblings when completed and I
    !    know what's going on
    !==========================================================================!
    ! key point is I can get rid of this line in gasci_pchb
    ! integer, parameter :: SAME_SPIN = 1, OPP_SPIN_NO_EXCH = 2, OPP_SPIN_EXCH = 3
    !   and any loops over `samplerIndex`
    ! which simplifies the code, purportedly...
    ! heavy lifting will be in this file, but I will need a pchb_uhf_excitgen too
    ! not sure where I will put it but I will need to also specify to use this
    ! excit gen when both PCHB and UHF are selected as input
    ! also very important: unit tests
    ! unit_tests/gasci/test_gasci_general_pchb.F90
    ! unit_tests/pcpp_excitgen/test_pchb_excitgen.F90
    ! unit_tests/pcpp_excitgen/test_aliasTables.F90 (don't need to change)
    !==========================================================================!
    use excitation_generators, only: ExcitationGenerator_t, SingleExcitationGenerator_t, &
                                doubleExcitationGenerator_t, get_pgen_sd
    use constants, only: n_int
    use dSFMT_interface, only: genrand_real2_dSFMT
    better_implicit_none

    private
    ! @jph
    ! public :: GAS_PCHB_ExcGenerator_t, use_supergroup_lookup, &
    !     GAS_doubles_PCHB_ExcGenerator_t, &
    !     possible_GAS_singles, GAS_PCHB_singles_generator, &
    !     GAS_PCHB_particle_selection, PCHB_particle_selections, &
    !     PCHB_ParticleSelection_t

    type, extends(ExcitationGenerator_t) :: GAS_PCHB_UHF_ExcGenerator_t
        !! full GAS PCHB excitation generator, working on spin orbitals
        private
        ! @jph not sure why GAS excitgen is different
        type(GAS_PCHB_UHF_DoublesExcGenerator_t) :: doubles_generator
        ! @jph make allocatable
        ! allocatable to allow choice of singles generator at runtime
        ! class(SingleExcitationGenerator_t), allocatable :: singles_generator
        class(SingleExcitationGenerator_t) :: singles_generator

    contains
        private
        procedure, public :: init => init_GAS_PCHB_UHF
        procedure, public :: finalize => finalize_GAS_PCHB_UHF
        procedure, public :: gen_exc => gen_exc_GAS_PCHB_UHF
        procedure, public :: get_pgen => get_pgen_GAS_PCHB_UHF
        procedure, public :: gen_all_excits => gen_all_excits_GAS_PCHB_UHF
    end type GAS_PCHB_UHF_ExcGenerator_t

    type, extends(doubleExcitationGenerator_t) :: GAS_PCHB_UHF_DoublesExcGenerator_t
        private
        ! @jph TODO
    contains
        private
        ! @jph TODO
    end type GAS_PCHB_UHF_DoublesExcGenerator_t




end module gasci_pchb_uhf