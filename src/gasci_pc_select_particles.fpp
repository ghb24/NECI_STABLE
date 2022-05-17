#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

! The main idea of the precomputed Heat bath sampling (PCHB) is taken from
!    J. Li, M. Otten, A. A. Holmes, S. Sharma, and C. J. Umrigar, J. Comput. Phys. 149, 214110 (2018).
! and described there.
! The main "ingredient" are precomputed probability distributions p(ab | ij) to draw a, b holes
! when i, j electrons were chosen.
! This requires #{i, j | i < j} probability distributions.
!
! The improved version to use spatial orbital indices to save memory is described in
!    Guther K. et al., J. Chem. Phys. 153, 034107 (2020).
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx.
! This gives #{i, j | i < j} * 3 probability distributions.
!
! The generalization to GAS spaces is available in a preprint (should be available in JCTC soon as well)
!    https://chemrxiv.org/engage/chemrxiv/article-details/61447e60b1d4a605d589af2e
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and GAS constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in gasci_supergroup_index.f90

module gasci_pc_select_particles
    use constants, only: dp, int64, stdout
    use aliasSampling, only: AliasSampler_1D_t, AliasSampler_2D_t
    use SystemData, only: nEl
    use sets_mod, only: is_set, operator(.in.)
    use util_mod, only: stop_all
    use UMatCache, only: numBasisIndices
    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use sets_mod, only: empty_int
    better_implicit_none
    public :: PC_WeightedParticles_t

    !> The precomputed GAS uniform excitation generator
    type :: PC_WeightedParticles_t
        private
        type(AliasSampler_1D_t) :: I_sampler
        !> The shape is (number_of_spin_orbs, n_supergroup)
        type(AliasSampler_2D_t) :: J_sampler

        class(GASSpec_t), allocatable :: GAS_spec
        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        !> Use a lookup for the supergroup index in global_det_data.
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.
    contains
        private
        procedure, public :: init
        procedure, public :: finalize

        procedure, public :: draw
        procedure, public :: get_pgen
    end type
contains

    subroutine init(this, GAS_spec, weights, use_lookup, create_lookup)
        class(PC_WeightedParticles_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        real(dp), intent(in) :: weights(:, :, :)
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = 'init'

        real(dp), allocatable :: w(:)
        integer, allocatable :: supergroups(:, :)
        integer :: nBI
        integer :: i_sg

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

        supergroups = this%indexer%get_supergroups()
        @:ASSERT(size(supergroups, 2) == size(weights, 3))

        call this%I_sampler%shared_alloc(size(supergroups, 2), nBI, 'PC_particles_i')
        call this%J_sampler%shared_alloc([nBI, size(supergroups, 2)], nBI, 'PC_particles_j')
        allocate(w(nBI))

        block
            integer :: I
            do i_sg = 1, size(supergroups, 2)
                do I = 1, nBI
                    call this%J_sampler%setup_entry(I, i_sg, weights(:, I, i_sg))
                end do
                call this%I_sampler%setup_entry(i_sg, sum(weights(:, :, i_sg), dim=1))
            end do
        end block
    end subroutine

    subroutine draw(this, nI, elecs, srcs, p)
        class(PC_WeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: srcs(2), elecs(2)
        real(dp), intent(out) :: p
        character(*), parameter :: this_routine = 'draw'

        integer :: i_sg
        real(dp) :: renorm_I, p_I, renorm_J, p_J

        ! TODO(@Oskar): Optimize
        i_sg = this%indexer%idx_nI(nI)

        ! TODO(@Oskar): Optimize
        renorm_I = sum(this%i_sampler%get_prob(i_sg, nI))
        call this%i_sampler%constrained_sample(i_sg, nI, renorm_I, elecs(1), srcs(1), p_I)
        @:ASSERT(srcs(1) .in. nI)

        ! TODO(@Oskar): Optimize
        renorm_J = sum(abs(this%J_sampler%get_prob(srcs(1), i_sg, nI)))

        ! Note that p(I | I) is automatically zero and cannot be drawn
        call this%j_sampler%constrained_sample(srcs(1), i_sg, nI, renorm_J, elecs(2), srcs(2), p_J)
        @:ASSERT(srcs(2) .in. nI)
        @:ASSERT(srcs(1) /= srcs(2))
        @:ASSERT(elecs(1) /= elecs(2))
        @:ASSERT(all(nI(elecs) == srcs))

        p = p_I * p_J
    end subroutine

    pure function get_pgen(this, nI, I, J) result(p)
        class(PC_WeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(in) :: I, J
        real(dp) :: p

        integer :: i_sg
        real(dp) :: renorm_I, p_I, renorm_J, p_J

        renorm_I = sum(this%I_sampler%get_prob(i_sg, nI))
        p_I = this%i_sampler%constrained_getProb(i_sg, nI, renorm_I, I)

        ! TODO(@Oskar): Optimize
        renorm_J = sum(this%J_sampler%get_prob(I, i_sg, nI))
        p_J = this%j_sampler%constrained_getProb(I, i_sg, nI, renorm_J, J)

        p = p_I * p_J
    end function


    subroutine finalize(this)
        class(PC_WeightedParticles_t), intent(inout) :: this

        call this%I_sampler%finalize()
        deallocate(this%indexer)
        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine
end module
