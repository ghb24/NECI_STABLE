#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_uhf
    !! same as gasci_pchb but with spin-orbitals
    !! useful for UHF-format FCIDUMP files
    !! only implements generators that work on spin-orbitals

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
    use excitation_generators, only: doubleExcitationGenerator_t
    use aliasSampling, only: AliasSampler_3D_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pc_select_particles, only: ParticleSelector_t
    use gasci_general, only: GASSpec_t
    better_implicit_none

    private
    ! @jph
    public :: GAS_PCHB_UHF_DoublesExcGenerator_t
    ! public :: GAS_PCHB_ExcGenerator_t, use_supergroup_lookup, &
    !     GAS_doubles_PCHB_ExcGenerator_t, &
    !     possible_GAS_singles, GAS_PCHB_singles_generator, &
    !     GAS_PCHB_particle_selection, PCHB_particle_selections, &
    !     PCHB_ParticleSelection_t


    type, extends(doubleExcitationGenerator_t) :: GAS_doubles_UHF_PCHB_ExcGenerator_t
        private
        logical, public :: use_lookup = .false.
            !! use a lookup for the supergroup index
        logical, public :: create_lookup
            !! create adn manage the supergroup index
        type(AliasSampler_3D_t) :: pchb_samplers
            !! the shape is (fused_number_of_double_excitations, 3, n_supergroup)

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(ParticleSelector_t), allocatable :: particle_selector
        class(GASSpec_t), allocatable :: GAS_spec
        real(dp), allocatable :: pExch(:, :)
        integer, allocatable :: tgtOrbs(:, :)
    contains
        private
        procedure, public :: init => GAS_doubles_PCHB_UHF_init
        procedure, public :: finalize => GAS_doubles_PCHB_UHF_finalize
        procedure, public :: gen_exc => GAS_doubles_PCHB_uhf_gen_exc
        procedure, public :: get_pgen => GAS_doubles_PCHB_uhf_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_PCHB_uhf_gen_all_excits

        procedure :: compute_samplers => GAS_doubles_PCHB_uhf_compute_samplers
    end type GAS_doubles_UHF_PCHB_ExcGenerator_t

contains

! @jph TODO implementation for doubles

    subroutine GAS_doubles_PCHB_UHF_init(this)
        !! initalises the UHF PCHB doubles excitation generator
        !!
        !! more specifically, sets up a lookup table for ab -> (a,b) and
        !! sets up the alias table for picking ab given ij with prob~Hijab
        ! @jph stub
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t) :: this

    end subroutine GAS_doubles_PCHB_UHF_init

    subroutine GAS_doubles_PCHB_UHF_finalize(this)
        !! deallocates the sampler and mapper
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t) :: this
        call this%pchb_samplers%finalize()
        call this%particle_selector%finalize()
        deallocate(this%particle_selector)
        deallocate(this%tgtOrbs)
        deallocate(this%pExch)
        deallocate(this%indexer)

        if (this%create_lookup) nullify(lookup_supergroup_indexer)

    end subroutine GAS_doubles_PCHB_UHF_finalize

    subroutine GAS_doubles_PCHB_uhf_gen_exc()
        ! @jph stub

    end subroutine GAS_doubles_PCHB_uhf_gen_exc

    subroutine GAS_doubles_PCHB_uhf_get_pgen()
        ! @jph stub

    end subroutine GAS_doubles_PCHB_uhf_get_pgen

    subroutine GAS_doubles_PCHB_uhf_gen_all_excits()
        ! @jph stub

    end subroutine GAS_doubles_PCHB_uhf_gen_all_excits

    subroutine GAS_doubles_PCHB_uhf_compute_samplers()
        ! @jph stub

    end subroutine GAS_doubles_PCHB_uhf_compute_samplers


end module gasci_pchb_uhf