#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pc_singles_localised
    use constants, only: dp, int64, stdout, n_int
    use SystemData, only: nEl
    use excitation_generators, only: SingleExcitationGenerator_t
    use FciMCData, only: excit_gen_store_type
    use aliasSampling, only: AliasSampler_2D_t
    use gasci, only: GASSpec_t
    better_implicit_none
    public :: PC_SinglesLocalised_t


    type, extends(SingleExcitationGenerator_t) :: PC_SinglesLocalised_t
        type(AliasSampler_2D_t) :: sampler
            !! p(A | I, i_sg)
            !! The probability of picking the hole A after having picked particle I
            !! in the supergroup i_sg.
    contains
        private
        procedure, public :: init => PC_SinglesLocalised_init
        procedure, public :: gen_exc => PC_SinglesLocalised_gen_exc
        procedure, public :: get_pgen => PC_SinglesLocalised_get_pgen
        procedure, public :: finalize => PC_SinglesLocalised_finalize
    end type

contains

    subroutine PC_SinglesLocalised_init(this, GAS_spec, weights)
        class(GASSpec_t), intent(in) :: GAS_spec
        class(PC_SinglesLocalised_t), intent(inout) :: this
        real(dp), intent(in) :: weights(:, :, :)
        integer :: n_supergroups

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'PC particles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'PC particles is using the supergroup lookup'

        nBI = this%GAS_spec%n_spin_orbs()
        @:ASSERT(nBI == size(weights, 1) .and. nBI == size(weights, 2))

        n_supergroups = this%indexer%n_supergroups()
        @:ASSERT(n_supergroups == size(weights, 3))

        call this%sampler%shared_alloc([size(weights, dim=2), size(weights, dim=3)], size(weights, dim=1), 'PC_singles')
        block
            integer :: i_sg, I
            do i_sg = 1, size(weights, 3)
                do I = 1, size(weights, 2)
                    call this%J_sampler%setup_entry(I, i_sg, weights(:, I, i_sg))
                end do
            end do
        end block
    end subroutine

    subroutine PC_SinglesLocalised_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(PC_SinglesLocalised_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        elec = int(genrand_real2_dSFMT() * nel) + 1
        src = nI(elec)
    end subroutine BoundGenExc_t

    subroutine PC_SinglesLocalised_finalize(this)
        class(PC_SinglesLocalised_t), intent(inout) :: this
    end subroutine

    real(dp) function PC_SinglesLocalised_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        class(PC_SinglesLocalised_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
    end function

end module
