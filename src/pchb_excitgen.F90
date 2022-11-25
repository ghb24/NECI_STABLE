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
    use excitation_generators, only: ExcitationGenerator_t, SingleExcitationGenerator_t, get_pgen_sd, gen_exc_sd, gen_all_excits_sd
    use exc_gen_class_wrappers, only: UniformSingles_t, WeightedSingles_t
    use gasci_pchb_rhf, only: GAS_doubles_RHF_PCHB_ExcGenerator_t
        ! @jph TODO should no longer need this once this is abstracted(!)
    use gasci_pc_select_particles, only: PCHB_ParticleSelection_t, PCHB_particle_selections
    better_implicit_none

    private

    public :: PCHB_FCI_excit_generator_t, FCI_PCHB_particle_selection, &
        FCI_PCHB_singles, possible_PCHB_singles

    type, extends(ExcitationGenerator_t) :: PCHB_FCI_excit_generator_t
        private
        type(GAS_doubles_RHF_PCHB_ExcGenerator_t) :: doubles_generator
        class(SingleExcitationGenerator_t), allocatable :: singles_generator
    contains
        private
        procedure, public :: init
        procedure, public :: finalize
        procedure, public :: gen_exc
        procedure, public :: get_pgen
        procedure, public :: gen_all_excits
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

    subroutine init(this, PCHB_particle_selection, PCHB_singles)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        type(PCHB_used_singles_t), intent(in) :: PCHB_singles
        character(*), parameter :: this_routine = 'pchb_excitgen::init'
        ! CAS is implemented as a special case of GAS with only one GAS space.
        ! Since a GAS specification with one GAS space is trivially disconnected, there
        ! is no point to use the lookup.
        call this%doubles_generator%init(&
                CAS_spec(n_el=nEl, n_spat_orbs=nBasis .div. 2), &
                use_lookup=.false., create_lookup=.false., &
                PCHB_particle_selection=PCHB_particle_selection)

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


    subroutine finalize(this)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this
        call this%singles_generator%finalize()
        deallocate(this%singles_generator)
        call this%doubles_generator%finalize()
    end subroutine

    function get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        pgen = get_pgen_sd(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2, &
                           this%singles_generator, this%doubles_generator)
    end function

    subroutine gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(PCHB_FCI_excit_generator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        call gen_exc_sd(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                        ex, tParity, pGen, hel, store, part_type, &
                        this%singles_generator, this%doubles_generator)
    end subroutine

    subroutine gen_all_excits(this, nI, n_excits, det_list)
        class(PCHB_FCI_excit_generator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        call gen_all_excits_sd(nI, n_excits, det_list, &
                               this%singles_generator, this%doubles_generator)
    end subroutine

end module pchb_excitgen
