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
    use constants, only: dp, n_int, maxExcit
    use util_mod, only: operator(.isclose.)
    use FciMCData, only: excit_gen_store_type
    use SystemData, only: nEl
    use SymExcitDataMod, only: ScratchSize
    use bit_rep_data, only: nIfTot
    use excitation_generators, only: doubleExcitationGenerator_t
    use FciMCData, only:
    use aliasSampling, only: AliasSampler_3D_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pc_select_particles, only: ParticleSelector_t
    use gasci, only: GASSpec_t
    better_implicit_none

    private
    ! @jph
    public :: GAS_doubles_UHF_PCHB_ExcGenerator_t
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
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t), intent(inout) :: this
        call this%pchb_samplers%finalize()
        call this%particle_selector%finalize()
        @:safe_deallocate(this%particle_selector)
        @:safe_deallocate(this%tgtOrbs)
        @:safe_deallocate(this%pExch)
        deallocate(this%indexer) ! pointer, so allocated(...) does not work

        if (this%create_lookup) nullify(lookup_supergroup_indexer)
    end subroutine GAS_doubles_PCHB_UHF_finalize

    subroutine GAS_doubles_PCHB_uhf_gen_exc(&
                this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                ex, tParity, pGen, hel, store, part_type)
        !! given the initial determinant (both as nI and ilut), create a random
        !! doubles excitation using the Hamiltonian matrix elements as weights
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t), intent(inout) :: this
            !! the exctitation generator
        integer, intent(in) :: nI(nel), exFlag
            !! determinant to excite from
            !! unused in this generator
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
            !! determint from which to excite, ilus format
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
            !! the excited determinant upon return
            !! excitation order (for doubles generator, always == 2)
            !! excitation matrix nI -> nJ
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
            !! excited determinant, ilut format
        real(dp), intent(out) :: pGen
            !! probability of generating the excitation
        logical, intent(out) :: tParity
            !! the parity of the excitation nI -> nJ
        HElement_t(dp), intent(out) :: hel
            !! matrix element Hijab
        type(excit_gen_store_type), intent(inout), target :: store
            !!
        integer, intent(in), optional :: part_type
            !! unused in this generator
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_uhf_gen_exc'

        integer :: i_sg ! supergroup index
        integer :: src ! particles (I, J)
        integer :: elecs(2) ! particle indeices in nI


        @:unused_var(exFlag, part_type)
        ic = 2
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0._dp) ! macro
#endif
        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI
        end if

        ! pick two random electrons
        call this%particle_selector%draw(nI, i_sg, elecs, src, pgen)
        if (src(1) == 0) then
            call invalidate()
            return
        end if

        ! @jph stub

    contains

        subroutine invalidate()
            nJ = 0
            ilutJ = 0_n_int
            ex(1, 1 : 2) = src
            ex(2, 1 : 2) = orbs
        end subroutine invalidate

    end subroutine GAS_doubles_PCHB_uhf_gen_exc

    real(dp) function GAS_doubles_PCHB_uhf_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        ! @jph docs
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_uhf_get_pgen'

        ! @jph stub

    end function GAS_doubles_PCHB_uhf_get_pgen

    subroutine GAS_doubles_PCHB_uhf_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        ! @jph stub

    end subroutine GAS_doubles_PCHB_uhf_gen_all_excits

    subroutine GAS_doubles_PCHB_uhf_compute_samplers(this)
        class(GAS_doubles_UHF_PCHB_ExcGenerator_t), intent(in) :: this
        ! @jph stub

    end subroutine GAS_doubles_PCHB_uhf_compute_samplers


end module gasci_pchb_uhf