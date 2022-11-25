#include "macros.h"
module pchb_excitgen
    use constants, only: n_int, dp, maxExcit
    use SystemData, only: nel, nBasis
    use util_mod, only: operator(.div.), EnumBase_t, stop_all
    use bit_rep_data, only: NIfTot
    use FciMCData, only: pSingles, excit_gen_store_type, pDoubles
    use SymExcitDataMod, only: ScratchSize
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t

    use gasci, only: GASSpec_t, LocalGASSpec_t
    use excitation_generators, only: ClassicAbInitExcitationGenerator_t, SingleExcitationGenerator_t
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t
    use gasci_pchb_rhf, only: GAS_doubles_RHF_PCHB_ExcGenerator_t
        ! @jph TODO should no longer need this once this is abstracted(!)
    use gasci_pc_select_particles, only: PCHB_ParticleSelection_t, PCHB_particle_selections
    better_implicit_none

    private

    public :: PCHB_FCI_excit_generator_t, FCI_PCHB_particle_selection, &
        FCI_PCHB_singles, possible_PCHB_singles

    type, extends(ClassicAbInitExcitationGenerator_t) :: PCHB_FCI_excit_generator_t
    contains
        private
        procedure, public :: init
    end type

    type(PCHB_ParticleSelection_t) :: FCI_PCHB_particle_selection = PCHB_particle_selections%PC_WEIGHTED


    type, extends(EnumBase_t) :: PCHB_used_singles_t
    end type

    type :: possible_PCHB_singles_t
        type(PCHB_used_singles_t) :: &
            ON_FLY_HEAT_BATH = PCHB_used_singles_t(1), &
            UNIFORM = PCHB_used_singles_t(2)
    end type

    type(possible_PCHB_singles_t), parameter :: possible_PCHB_singles = possible_PCHB_singles_t()

    type(PCHB_used_singles_t) :: FCI_PCHB_singles = possible_PCHB_singles%UNIFORM

contains

    subroutine init(this, PCHB_particle_selection, PCHB_singles, is_uhf)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        type(PCHB_used_singles_t), intent(in) :: PCHB_singles
        logical, intent(in), optional :: is_uhf
        logical :: is_uhf_
        character(*), parameter :: this_routine = 'pchb_excitgen::init'

        def_default(is_uhf_, is_uhf, .false.)
        ! CAS is implemented as a special case of GAS with only one GAS space.
        ! Since a GAS specification with one GAS space is trivially disconnected, there
        ! is no point to use the lookup.
        ! @jph TODO implement RHF
        if (is_uhf_) then
            call stop_all(this_routine, 'UHF PCHB not yet implemented :(')
        else
            allocate(GAS_doubles_RHF_PCHB_ExcGenerator_t :: this%doubles_generator)
            select type(generator => this%doubles_generator)
            type is(GAS_doubles_RHF_PCHB_ExcGenerator_t)
                call generator%init(&
                CAS_spec(n_el=nEl, n_spat_orbs=nBasis .div. 2), &
                use_lookup=.false., create_lookup=.false., &
                PCHB_particle_selection=PCHB_particle_selection)
            end select
        end if

        ! luckily the singles generators don't require initialization.
        if (PCHB_singles == possible_PCHB_singles%ON_FLY_HEAT_BATH) then
            allocate(WeightedSingles_t :: this%singles_generator)
        else if (PCHB_singles == possible_PCHB_singles%UNIFORM) then
            allocate(UniformSingles_t :: this%singles_generator)
        else
            call stop_all(this_routine, "Invalid PCHB_singles in FCI PCHB init.")
        end if
    contains

        type(LocalGASSpec_t) pure function CAS_spec(n_el, n_spat_orbs)
            integer, intent(in) :: n_el, n_spat_orbs
            integer :: i
            ! It does not matter if we use local or cumulative GAS
            ! constraints
            CAS_spec = LocalGASSpec_t(n_min=[n_el], n_max=[n_el], spat_GAS_orbs=[(1, i = 1, n_spat_orbs)])
        end function
    end subroutine

end module pchb_excitgen
