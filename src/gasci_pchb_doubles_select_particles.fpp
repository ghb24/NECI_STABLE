#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_doubles_select_particles
    use constants, only: dp, int64, stdout, n_int, bits_n_int
    use fortran_strings, only: to_upper
    use bit_rep_data, only: nIfD
    use util_mod, only: EnumBase_t
    use aliasSampling, only: AliasSampler_1D_t, AliasSampler_2D_t
    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pParallel
    use sets_mod, only: is_set, operator(.in.)
    use excit_gens_int_weighted, only: pick_biased_elecs, get_pgen_pick_biased_elecs
    use util_mod, only: stop_all, operator(.isclose.), swap, &
        binary_search_int, EnumBase_t, operator(.div.)
    use UMatCache, only: numBasisIndices
    use MPI_wrapper, only: root
    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use sets_mod, only: empty_int
    better_implicit_none
    private
    public :: ParticleSelector_t, PC_FullyWeightedParticles_t, &
        UniformParticles_t, PC_FastWeightedParticles_t, &
        PCHB_particle_selection_vals, PCHB_ParticleSelection_t, allocate_and_init, &
        PCHB_ParticleSelection_vals_t



    type, extends(EnumBase_t) :: PCHB_ParticleSelection_t
    end type

    type :: PCHB_ParticleSelection_vals_t
        type(PCHB_ParticleSelection_t) :: &
            UNIFORM = PCHB_ParticleSelection_t(1), &
            FULLY_WEIGHTED = PCHB_ParticleSelection_t(2), &
                !! We draw from \( p(I)|_{D_i} \) and then \( p(J | I)_{J \in D_i} \)
                !! and both probabilites come from the PCHB weighting scheme.
                !! We guarantee that \(I\) and \(J\) are occupied.
            WEIGHTED = PCHB_ParticleSelection_t(3), &
                !! We draw \( \tilde{p}(I)|_{D_i} \) uniformly and then \( p(J | I)_{J \in D_i} \)
                !! The second distribution comes from the PCHB weighting scheme.
                !! We guarantee that \(I\) and \(J\) are occupied.
            FAST_WEIGHTED = PCHB_ParticleSelection_t(4)
                !! We draw \( \tilde{p}(I)|_{D_i} \) uniformly and then \( p(J | I)_{J} \).
                !! The second distribution comes from the PCHB weighting scheme.
                !! We guarantee that \(I\) is occupied.
    contains
        procedure, nopass :: from_str => from_keyword
    end type

    type(PCHB_ParticleSelection_vals_t), parameter :: &
        PCHB_particle_selection_vals = PCHB_ParticleSelection_vals_t()


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

        subroutine Draw_t(this, nI, ilutI, i_sg, elecs, srcs, p)
            import :: dp, ParticleSelector_t, nEl, n_int, nIfD
            class(ParticleSelector_t), intent(in) :: this
            integer, intent(in) :: nI(nEl), i_sg
                !! The determinant in nI-format and the supergroup index
            integer(n_int), intent(in) :: ilutI(0 : nIfD)
                !! The determinant in bitmask format
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

    type, extends(ParticleSelector_t), abstract :: PC_Particles_t
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

    type, extends(PC_Particles_t) :: PC_FullyWeightedParticles_t
    contains
        private
        procedure, public :: draw => draw_PC_FullyWeightedParticles_t
        procedure, public :: get_pgen => get_pgen_PC_FullyWeightedParticles_t
    end type

    type, extends(PC_Particles_t) :: PC_WeightedParticles_t
    contains
        private
        procedure, public :: draw => draw_PC_WeightedParticles_t
        procedure, public :: get_pgen => get_pgen_PC_WeightedParticles_t
    end type

    type, extends(PC_Particles_t) :: PC_FastWeightedParticles_t
    contains
        private
        procedure, public :: draw => draw_PC_FastWeightedParticles_t
        procedure, public :: get_pgen => get_pgen_PC_FastWeightedParticles_t
    end type

contains

    pure function from_keyword(w) result(res)
        !! Parse a given keyword into the possible particle selection schemes
        character(*), intent(in) :: w
        type(PCHB_ParticleSelection_t) :: res
        routine_name("from_keyword")
        select case(to_upper(w))
        case('UNIFORM')
            res = PCHB_particle_selection_vals%UNIFORM
        case('FULLY-WEIGHTED')
            res = PCHB_particle_selection_vals%FULLY_WEIGHTED
        case('WEIGHTED')
            res = PCHB_particle_selection_vals%WEIGHTED
        case('FAST-WEIGHTED')
            res = PCHB_particle_selection_vals%FAST_WEIGHTED
        case default
            call stop_all(this_routine, trim(w)//" not a valid doubles particle selection scheme.")
        end select
    end function

    subroutine allocate_and_init(PCHB_particle_selection, GAS_spec, IJ_weights, use_lookup, particle_selector)
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        class(GASSpec_t), intent(in) :: GAS_spec
        real(dp), intent(in) :: IJ_weights(:, :, :)
        logical, intent(in) :: use_lookup
        class(ParticleSelector_t), allocatable, intent(inout) :: particle_selector
        routine_name("gasci_pchb_doubles_select_particles::allocate_and_init")
        if (PCHB_particle_selection == PCHB_particle_selection_vals%FULLY_WEIGHTED) then
            allocate(PC_FullyWeightedParticles_t :: particle_selector)
            select type(particle_selector)
            type is(PC_FullyWeightedParticles_t)
                call particle_selector%init(GAS_spec, IJ_weights, use_lookup, .false.)
            end select
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%WEIGHTED) then
            allocate(PC_WeightedParticles_t :: particle_selector)
            select type(particle_selector)
            type is(PC_WeightedParticles_t)
                call particle_selector%init(GAS_spec, IJ_weights, use_lookup, .false.)
            end select
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%FAST_WEIGHTED) then
            allocate(PC_FastWeightedParticles_t :: particle_selector)
            select type(particle_selector)
            type is(PC_FastWeightedParticles_t)
                call particle_selector%init(GAS_spec, IJ_weights, use_lookup, .false.)
            end select
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%UNIFORM) then
            allocate(UniformParticles_t :: particle_selector)
        else
            call stop_all(this_routine, 'Invalid particle selection.')
        end if
    end subroutine

    subroutine draw_UniformParticles_t(this, nI, ilutI, i_sg, elecs, srcs, p)
        class(UniformParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
        integer, intent(out) :: srcs(2), elecs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        integer :: sym_prod, ispn, sum_ml
        @:unused_var(this, ilutI, i_sg)
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
        class(PC_Particles_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        real(dp), intent(in) :: weights(:, :, :)
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = 'init'

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
        block
            integer :: I
            do i_sg = 1, size(supergroups, 2)
                do I = 1, nBI
                    call this%J_sampler%setup_entry(I, i_sg, root, weights(:, I, i_sg))
                end do
                call this%I_sampler%setup_entry(i_sg, root, sum(weights(:, :, i_sg), dim=1))
            end do
        end block
    end subroutine

    subroutine draw_PC_FullyWeightedParticles_t(this, nI, ilutI, i_sg, elecs, srcs, p)
        class(PC_FullyWeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
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

        renorm_first = sum(this%I_sampler%get_prob(i_sg, nI))
        call this%I_sampler%constrained_sample(&
            i_sg, nI, ilutI, renorm_first, elecs(1), srcs(1), p_first(1))
        @:ASSERT(nI(elecs(1)) == srcs(1))

        renorm_second(1) = sum(this%J_sampler%get_prob(srcs(1), i_sg, nI))

        ! Note that p(I | I) is automatically zero and cannot be drawn
        call this%J_sampler%constrained_sample(&
            srcs(1), i_sg, nI, ilutI, renorm_second(1), elecs(2), srcs(2), p_second(1))
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
        p_first(2) = this%I_sampler%constrained_getProb(&
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

    pure function get_pgen_PC_FullyWeightedParticles_t(this, nI, i_sg, I, J) result(p)
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
        class(PC_FullyWeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: I, J
            !! The particles.
        real(dp) :: p

        real(dp) :: renorm_first, p_first(2), renorm_second(2), p_second(2)

        ! The renormalization for the first electron is the same,
        ! regardless of order.
        renorm_first = sum(this%I_sampler%get_prob(i_sg, nI))

        p_first(1) = this%I_sampler%constrained_getProb(i_sg, nI, renorm_first, I)
        p_first(2) = this%I_sampler%constrained_getProb(i_sg, nI, renorm_first, J)

        renorm_second(1) = sum(this%J_sampler%get_prob(I, i_sg, nI))
        p_second(1) = this%J_sampler%constrained_getProb(&
            I, i_sg, nI, renorm_second(1), J)

        renorm_second(2) = sum(this%J_sampler%get_prob(J, i_sg, nI))
        p_second(2) = this%J_sampler%constrained_getProb(&
            J, i_sg, nI, renorm_second(2), I)

        p = sum(p_first * p_second)
    end function

    subroutine draw_PC_WeightedParticles_t(this, nI, ilutI, i_sg, elecs, srcs, p)
        class(PC_WeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
        integer, intent(out) :: srcs(2), elecs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        character(*), parameter :: this_routine = 'draw'
        real(dp) :: renorm_second(2), p_second(2)

        elecs(1) = int(genrand_real2_dSFMT() * nEl) + 1
        srcs(1) = nI(elecs(1))

        renorm_second(1) = sum(this%J_sampler%get_prob(srcs(1), i_sg, nI))
        ! Note that p(I | I) is automatically zero and cannot be drawn
        call this%J_sampler%constrained_sample(&
            srcs(1), i_sg, nI, ilutI, renorm_second(1), elecs(2), srcs(2), p_second(1))

        @:ASSERT((srcs(1) .in. nI) .and. (srcs(2) .in. nI))
        @:ASSERT(srcs(1) /= srcs(2))
        @:ASSERT(elecs(1) /= elecs(2))
        @:ASSERT(all(nI(elecs) == srcs))

        ! We could have picked them the other way round.
        ! Account for that.
        ! The renormalization for the first electron is the same
        renorm_second(2) = sum(this%J_sampler%get_prob(srcs(2), i_sg, nI))
        p_second(2) = this%J_sampler%constrained_getProb(&
            srcs(2), i_sg, nI, renorm_second(2), srcs(1))

        p = sum(p_second) / nEl

        if (srcs(1) > srcs(2)) then
            call swap(srcs(1), srcs(2))
            call swap(elecs(1), elecs(2))
        end if

        @:ASSERT(p .isclose. this%get_pgen(nI, i_sg, srcs(1), srcs(2)))
    end subroutine

    pure function get_pgen_PC_WeightedParticles_t(this, nI, i_sg, I, J) result(p)
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
        class(PC_WeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: I, J
            !! The particles.
        real(dp) :: p

        real(dp) :: renorm_second(2), p_second(2)

        ! The probability for the first electron is always 1 / nEl

        renorm_second(1) = sum(this%J_sampler%get_prob(I, i_sg, nI))
        p_second(1) = this%J_sampler%constrained_getProb(&
            I, i_sg, nI, renorm_second(1), J)

        renorm_second(2) = sum(this%J_sampler%get_prob(J, i_sg, nI))
        p_second(2) = this%J_sampler%constrained_getProb(&
            J, i_sg, nI, renorm_second(2), I)

        p = sum(p_second) / nEl
    end function



    subroutine draw_PC_FastWeightedParticles_t(this, nI, ilutI, i_sg, elecs, srcs, p)
        class(PC_FastWeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
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

        call this%J_sampler%sample(srcs(1), i_sg, srcs(2), p_J_I)

        if (.not. IsOcc(ilutI, srcs(2))) then
            elecs(:) = 0; srcs(:) = 0; p = 1._dp
        else
            elecs(2) = int(binary_search_int(nI, srcs(2)))
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
        end if
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
        class(PC_Particles_t), intent(inout) :: this

        if (allocated(this%GAS_spec)) then
            call this%I_sampler%finalize()
            call this%J_sampler%finalize()
            ! Yes, we assume, that either all or none are allocated
            deallocate(this%indexer, this%GAS_spec)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
        end if
    end subroutine
end module
