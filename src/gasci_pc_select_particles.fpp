#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pc_select_particles
    use constants, only: dp, int64, stdout
    use util_mod, only: EnumBase_t
    use aliasSampling, only: AliasSampler_1D_t, AliasSampler_2D_t
    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pParallel
    use sets_mod, only: is_set, operator(.in.)
    use excit_gens_int_weighted, only: pick_biased_elecs, get_pgen_pick_biased_elecs
    use util_mod, only: stop_all, operator(.isclose.), swap, binary_search_int
    use UMatCache, only: numBasisIndices
    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use sets_mod, only: empty_int
    better_implicit_none
    public :: ParticleSelector_t, PC_WeightedParticlesOcc_t, &
        UniformParticles_t, PC_FastWeightedParticles_t, &
        PCHB_ParticleSelection_t, PCHB_particle_selections, &
        GAS_PCHB_particle_selection

    type, abstract :: ParticleSelector_t
    contains
        procedure(Draw_t), public, deferred :: draw
        procedure(GetPgen_t), public, deferred :: get_pgen
        procedure(Finalize_t), public, deferred :: finalize
    end type

    abstract interface
        subroutine Finalize_t(this)
            import :: ParticleSelector_t
            implicit none
            class(ParticleSelector_t), intent(inout) :: this
        end subroutine

        real(dp) pure function GetPgen_t(this, nI, i_sg, I, J)
            import :: dp, ParticleSelector_t, nEl
            implicit none
            class(ParticleSelector_t), intent(in) :: this
            integer, intent(in) :: nI(nEl), i_sg
                !! The determinant in nI-format and the supergroup index
            integer, intent(in) :: I, J
                !! The particles.
        end function

        subroutine Draw_t(this, nI, i_sg, elecs, srcs, p)
            import :: dp, ParticleSelector_t, nEl
            class(ParticleSelector_t), intent(in) :: this
            integer, intent(in) :: nI(nEl), i_sg
                !! The determinant in nI-format and the supergroup index
            integer, intent(out) :: srcs(2), elecs(2)
                !! The chosen particles \(I, J\) and their index in `nI`.
                !! It is guaranteed that `scrs(1) < srcs(2)`.
            real(dp), intent(out) :: p
                !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
                !! This is the probability of drawing two particles from
                !! a given determinant \(D_i\) regardless of order.
        end subroutine
    end interface

    type, extends(ParticleSelector_t) :: UniformParticles_t
    contains
        private
        procedure, public :: draw => draw_UniformParticles_t
        procedure, public :: get_pgen => get_pgen_UniformParticles_t
        procedure, public :: finalize => finalize_UniformParticles_t
    end type

    type, extends(ParticleSelector_t), abstract :: PC_WeightedParticles_t
        private
        type(AliasSampler_1D_t) :: I_sampler
            !! The shape is (n_supergroup) -> number_of_spin_orbs
        type(AliasSampler_2D_t) :: J_sampler
            !! The shape is (number_of_spin_orbs, n_supergroup) -> number_of_spin_orbs

        class(GASSpec_t), allocatable :: GAS_spec
        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        logical, public :: use_lookup = .false.
            !! Use a lookup for the supergroup index in global_det_data.
        logical, public :: create_lookup = .false.
            !! Create **and** manage! the supergroup index lookup in global_det_data.
    contains
        private
        procedure, public :: init => init_PC_WeightedParticles_t
        procedure, public :: finalize => finalize_PC_WeightedParticles_t
    end type

    type, extends(PC_WeightedParticles_t) :: PC_WeightedParticlesOcc_t
    contains
        private
        procedure, public :: draw => draw_PC_WeightedParticlesOcc_t
        procedure, public :: get_pgen => get_pgen_PC_WeightedParticlesOcc_t
    end type

    type, extends(PC_WeightedParticles_t) :: PC_FastWeightedParticles_t
    contains
        private
        procedure, public :: draw => draw_PC_FastWeightedParticles_t
        procedure, public :: get_pgen => get_pgen_PC_FastWeightedParticles_t
    end type

    ! PCHB particle selections

    type, extends(EnumBase_t) :: PCHB_ParticleSelection_t
    end type

    type :: possible_PCHB_ParticleSelection_t
        type(PCHB_ParticleSelection_t) :: &
            UNIFORM = PCHB_ParticleSelection_t(1), &
            PC_WEIGHTED = PCHB_ParticleSelection_t(2), &
            PC_WEIGHTED_APPROX = PCHB_ParticleSelection_t(3)
    end type

    type(possible_PCHB_ParticleSelection_t), parameter :: &
        PCHB_particle_selections = possible_PCHB_ParticleSelection_t()

    type(PCHB_ParticleSelection_t) :: GAS_PCHB_particle_selection = PCHB_particle_selections%PC_WEIGHTED


contains

    subroutine draw_UniformParticles_t(this, nI, i_sg, elecs, srcs, p)
        class(UniformParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(out) :: srcs(2), elecs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        integer :: sym_prod, ispn, sum_ml
        @:unused_var(this, i_sg)
        call pick_biased_elecs(nI, elecs, srcs, sym_prod, ispn, sum_ml, p)
    end subroutine

    pure function get_pgen_UniformParticles_t(this, nI, i_sg, I, J) result(p)
        !! Calculates \( p(\{I, J\}) \Big |_{D_i} \)
        !!
        !! This is the probability of drawing two particles from
        !! a given determinant \(D_i\) regardless of order.
        !!
        !! Note that the unordered probability is given by the ordered
        !! probability as:
        !! $$ p(\{I, J\}) \Big |_{D_i} = p((I, J)) \Big |_{D_i}
        !!              + p((J, I)) \Big |_{D_i} \quad.$$
        !! In addition we have
        !! $$ p((I, J)) \Big |_{D_i} \neq p((J, I)) \Big |_{D_i} $$
        !! so we have to actually calculate the probability of drawing
        !! two given particles in different order.
        class(UniformParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: I, J
            !! The particles.
        real(dp) :: p
        @:unused_var(this, nI, i_sg)

        p = get_pgen_pick_biased_elecs(&
            is_beta(I) .eqv. is_beta(J), pParallel, &
            par_elec_pairs, AB_elec_pairs)
    end function

    subroutine finalize_UniformParticles_t(this)
        class(UniformParticles_t), intent(inout) :: this
        @:unused_var(this)
        ! Nothing to do
    end subroutine


    subroutine init_PC_WeightedParticles_t(this, GAS_spec, weights, use_lookup, create_lookup)
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

    subroutine draw_PC_WeightedParticlesOcc_t(this, nI, i_sg, elecs, srcs, p)
        class(PC_WeightedParticlesOcc_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(out) :: srcs(2), elecs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        character(*), parameter :: this_routine = 'draw'

        real(dp) :: renorm_first, p_first(2)
        real(dp) :: renorm_second(2), p_second(2)

        renorm_first = sum(this%i_sampler%get_prob(i_sg, nI))
        call this%i_sampler%constrained_sample(&
            i_sg, nI, renorm_first, elecs(1), srcs(1), p_first(1))
        @:ASSERT(srcs(1) .in. nI)

        renorm_second(1) = sum(this%J_sampler%get_prob(srcs(1), i_sg, nI))

        ! Note that p(I | I) is automatically zero and cannot be drawn
        call this%j_sampler%constrained_sample(&
            srcs(1), i_sg, nI, renorm_second(1), elecs(2), srcs(2), p_second(1))
        if (srcs(2) == 0) then
            elecs(:) = 0; srcs(:) = 0; p = 1._dp
            return
        end if
        @:ASSERT(srcs(2) .in. nI)
        @:ASSERT(srcs(1) /= srcs(2))
        @:ASSERT(elecs(1) /= elecs(2))
        @:ASSERT(all(nI(elecs) == srcs))

        ! We could have picked them the other way round.
        ! Account for that.
        ! The renormalization for the first electron is the same
        p_first(2) = this%i_sampler%constrained_getProb(&
            i_sg, nI, renorm_first, srcs(2))
        renorm_second(2) = sum(this%J_sampler%get_prob(srcs(2), i_sg, nI))
        p_second(2) = this%J_sampler%constrained_getProb(&
            srcs(2), i_sg, nI, renorm_second(2), srcs(1))

        p = sum(p_first * p_second)

        if (srcs(1) > srcs(2)) then
            call swap(srcs(1), srcs(2))
            call swap(elecs(1), elecs(2))
        end if

        @:ASSERT(p .isclose. this%get_pgen(nI, i_sg, srcs(1), srcs(2)))
    end subroutine

    pure function get_pgen_PC_WeightedParticlesOcc_t(this, nI, i_sg, I, J) result(p)
        !! Calculates \( p(\{I, J\}) \Big |_{D_i} \)
        !!
        !! This is the probability of drawing two particles from
        !! a given determinant \(D_i\) regardless of order.
        !!
        !! Note that the unordered probability is given by the ordered
        !! probability as:
        !! $$ p(\{I, J\}) \Big |_{D_i} = p((I, J)) \Big |_{D_i}
        !!              + p((J, I)) \Big |_{D_i} \quad.$$
        !! In addition we have
        !! $$ p((I, J)) \Big |_{D_i} \neq p((J, I)) \Big |_{D_i} $$
        !! so we have to actually calculate the probability of drawing
        !! two given particles in different order.
        class(PC_WeightedParticlesOcc_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: I, J
            !! The particles.
        real(dp) :: p

        real(dp) :: renorm_first, p_first(2), renorm_second(2), p_second(2)

        ! The renormalization for the first electron is the same,
        ! regardless of order.
        renorm_first = sum(this%I_sampler%get_prob(i_sg, nI))

        p_first(1) = this%i_sampler%constrained_getProb(i_sg, nI, renorm_first, I)
        p_first(2) = this%i_sampler%constrained_getProb(i_sg, nI, renorm_first, J)

        renorm_second(1) = sum(this%J_sampler%get_prob(I, i_sg, nI))
        p_second(1) = this%j_sampler%constrained_getProb(&
            I, i_sg, nI, renorm_second(1), J)

        renorm_second(2) = sum(this%J_sampler%get_prob(J, i_sg, nI))
        p_second(2) = this%j_sampler%constrained_getProb(&
            J, i_sg, nI, renorm_second(2), I)

        p = sum(p_first * p_second)
    end function


    subroutine draw_PC_FastWeightedParticles_t(this, nI, i_sg, elecs, srcs, p)
        class(PC_FastWeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(out) :: srcs(2), elecs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        character(*), parameter :: this_routine = 'draw'
        real(dp) :: p_J_I, p_I_J

        elecs(1) = int(genrand_real2_dSFMT() * nEl) + 1
        srcs(1) = nI(elecs(1))

        call this%j_sampler%sample(srcs(1), i_sg, srcs(2), p_J_I)

        elecs(2) = int(binary_search_int(nI, srcs(2)))
        if (elecs(2) == -1) then
            elecs(:) = 0; srcs(:) = 0; p = 1._dp
            return
        end if
        @:ASSERT((srcs(1) .in. nI) .and. (srcs(2) .in. nI))
        @:ASSERT(srcs(1) /= srcs(2))
        @:ASSERT(elecs(1) /= elecs(2))
        @:ASSERT(all(nI(elecs) == srcs))


        p_I_J = this%J_sampler%get_prob(srcs(2), i_sg, srcs(1))
        p = (p_J_I + p_I_J) / nEl

        if (srcs(1) > srcs(2)) then
            call swap(srcs(1), srcs(2))
            call swap(elecs(1), elecs(2))
        end if

        @:ASSERT(p .isclose. this%get_pgen(nI, i_sg, srcs(1), srcs(2)))
    end subroutine

    pure function get_pgen_PC_FastWeightedParticles_t(this, nI, i_sg, I, J) result(p)
        !! Calculates \( p(\{I, J\}) \Big |_{D_i} \)
        !!
        !! This is the probability of drawing two particles from
        !! a given determinant \(D_i\) regardless of order.
        !!
        !! Note that the unordered probability is given by the ordered
        !! probability as:
        !! $$ p(\{I, J\}) \Big |_{D_i} = p((I, J)) \Big |_{D_i}
        !!              + p((J, I)) \Big |_{D_i} \quad.$$
        !! In addition we have
        !! $$ p((I, J)) \Big |_{D_i} \neq p((J, I)) \Big |_{D_i} $$
        !! so we have to actually calculate the probability of drawing
        !! two given particles in different order.
        class(PC_FastWeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: I, J
            !! The particles.
        real(dp) :: p

        @:unused_var(nI)
        p = (this%J_sampler%get_prob(I, i_sg, J) &
              + this%J_sampler%get_prob(J, i_sg, I)) / nEl
    end function


    subroutine finalize_PC_WeightedParticles_t(this)
        class(PC_WeightedParticles_t), intent(inout) :: this

        call this%I_sampler%finalize()
        call this%J_sampler%finalize()
        deallocate(this%indexer)
        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine
end module
