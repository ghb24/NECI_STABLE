#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module guga_gas_pchb_doubles_select_particles
    use constants, only: dp, n_int, stdout
    use fortran_strings, only: to_upper
    use bit_rep_data, only: nIfD
    use aliasSampling, only: AliasSampler_1D_t, AliasSampler_2D_t
    use SystemData, only: nEl, ElecPairs, nSpatOrbs
    use sets_mod, only: operator(.in.)
    use util_mod, only: stop_all, EnumBase_t, swap, operator(.isclose.), near_zero
    use UMatCache, only: gtID
    use generalized_supergroup, only: AlsoGUGA_SuperGroupIndexer_t, lookup_supergroup_indexer
    use MPI_wrapper, only: iProcIndex_intra
    use dSFMT_interface, only: genrand_real2_dSFMT

    use guga_bitrepops, only: CSF_Info_t
    use guga_excitations, only: pick_elec_pair_uniform_guga

    better_implicit_none
    private
    public :: ParticleSelector_t, &
        UniformParticles_t, PCHB_ParticleSelection_t, allocate_and_init, &
        PCHB_ParticleSelection_vals_t, PCHB_particle_selection_vals



    type, extends(EnumBase_t) :: PCHB_ParticleSelection_t
        private
        character(20) :: str
    contains
        procedure :: to_str
    end type

    type :: PCHB_ParticleSelection_vals_t
        type(PCHB_ParticleSelection_t) :: &
            UNIF_UNIF = PCHB_ParticleSelection_t(1, 'UNIF-UNIF'), &
                !! Both particles are drawn uniformly.
            FULL_FULL = PCHB_ParticleSelection_t(2, 'FULL-FULL'), &
                !! We draw from \( p(I)|_{D_i} \) and then \( p(J | I)_{J \in D_i} \)
                !! and both probabilites come from the PCHB weighting scheme.
                !! We guarantee that \(I\) and \(J\) are occupied.
            UNIF_FULL = PCHB_ParticleSelection_t(3, 'UNIF-FULL')
                !! We draw \( \tilde{p}(I)|_{D_i} \) uniformly and then \( p(J | I)_{J \in D_i} \)
                !! The second distribution comes from the PCHB weighting scheme.
                !! We guarantee that \(I\) and \(J\) are occupied.
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

        real(dp) pure function GetPgen_t(this, nI, csf_I, i_sg, I, J)
            import :: dp, CSF_Info_t, ParticleSelector_t, nEl
            implicit none
            integer, intent(in) :: nI(nEl), i_sg
                !! The CSF in nI-format and the supergroup
            class(ParticleSelector_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_I
            integer, intent(in) :: I, J
                !! The particles.
        end function

        subroutine Draw_t(this, nI, ilutI, csf_i, i_sg, srcs, p)
            import :: dp, CSF_Info_t, ParticleSelector_t, nEl, n_int, nIfD
            class(ParticleSelector_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
                !! The determinant in nI-format and the supergroup index
            integer(n_int), intent(in) :: ilutI(0 : nIfD)
                !! The determinant in bitmask format
            type(CSF_Info_t), intent(in) :: csf_I
            integer, intent(in) :: i_sg
                !! The supergroup index
            integer, intent(out) :: srcs(2)
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
            !! The shape is (n_supergroup) -> 2 * number_of_spat_orbs
            !!
            !! @note
            !! To make doubly occupied orbitals more probable we have doubly the
            !! amount of spatial orbitals and treat it as if there were spin-orbitals.
            !! @endnote
            !!
        type(AliasSampler_2D_t) :: j_sampler
            !! The shape is (number_of_spat_orbs, n_supergroup) -> number_of_spat_orbs

        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        class(AlsoGUGA_SuperGroupIndexer_t), pointer :: indexer => null()
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

contains

    pure function to_str(options) result(res)
        class(PCHB_ParticleSelection_t), intent(in) :: options
        character(:), allocatable :: res
        res = trim(options%str)
    end function

    pure function from_keyword(w) result(res)
        !! Parse a given keyword into the possible particle selection schemes
        character(*), intent(in) :: w
        type(PCHB_ParticleSelection_t) :: res
        routine_name("from_keyword")
        select case(to_upper(w))
        case(PCHB_particle_selection_vals%UNIF_UNIF%str)
            res = PCHB_particle_selection_vals%UNIF_UNIF
        case(PCHB_particle_selection_vals%FULL_FULL%str)
            res = PCHB_particle_selection_vals%FULL_FULL
        case(PCHB_particle_selection_vals%UNIF_FULL%str)
            res = PCHB_particle_selection_vals%UNIF_FULL
        case default
            call stop_all(this_routine, trim(w)//" not a valid doubles particle selection scheme.")
        end select
    end function

    subroutine allocate_and_init(PCHB_particle_selection, indexer, ij_weights, rank_with_info, use_lookup, particle_selector)
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        class(AlsoGUGA_SuperGroupIndexer_t), intent(in) :: indexer
        real(dp), intent(in) :: ij_weights(:, :, :)
            !! The weights for the spatial orbitals ij
            !! nSpatOrb * n_SpatOrb * n_supergroup
        integer, intent(in) :: rank_with_info
            !! The **intra-node** rank that contains the weights to be used
            !! all other arr of all other ranks are ignored (and can be allocated with size 0).
        logical, intent(in) :: use_lookup
        class(ParticleSelector_t), allocatable, intent(inout) :: particle_selector

        routine_name("gasci_pchb_doubles_select_particles::allocate_and_init")
        if (PCHB_particle_selection == PCHB_particle_selection_vals%UNIF_UNIF) then
            allocate(UniformParticles_t :: particle_selector)
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%FULL_FULL) then
            allocate(PC_FullyWeightedParticles_t :: particle_selector)
            select type(particle_selector)
            type is(PC_FullyWeightedParticles_t)
                call particle_selector%init(indexer, IJ_weights, rank_with_info, use_lookup, .false.)
            end select
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%UNIF_FULL) then
            allocate(PC_WeightedParticles_t :: particle_selector)
            select type(particle_selector)
            type is(PC_WeightedParticles_t)
                call particle_selector%init(indexer, IJ_weights, rank_with_info, use_lookup, .false.)
            end select
        else
            call stop_all(this_routine, 'Invalid particle selection.')
        end if
    end subroutine

    subroutine draw_UniformParticles_t(this, nI, ilutI, csf_i, i_sg, srcs, p)
        class(UniformParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i_sg
            !! The supergroup index
        integer, intent(out) :: srcs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        integer :: sym_prod, sum_ml
        routine_name("draw_UniformParticles_t")
        @:unused_var(this, ilutI, csf_i, i_sg)
        call pick_elec_pair_uniform_guga(nI, srcs, sym_prod, sum_ml, p)

        @:ASSERT(srcs(1) < srcs(2))
        srcs = gtID(srcs)
        associate(i => srcs(1), j => srcs(2))
            if (i /= j) then
                p = csf_i%occ_real(i) * csf_i%occ_real(j) * p
            end if
        end associate
    end subroutine

    pure function get_pgen_UniformParticles_t(this, nI, csf_i, i_sg, i, j) result(p)
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
            !! The CSF in nI-format
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j
            !! The particles.
        real(dp) :: p
        @:unused_var(this, nI, csf_i, i, j, i_sg)
        if (i /= j) then
            p = csf_i%occ_real(i) * csf_i%occ_real(j) / real(ElecPairs, dp)
        else
            p = 1.0_dp / real(ElecPairs, dp)
        end if
    end function

    subroutine finalize_UniformParticles_t(this)
        class(UniformParticles_t), intent(inout) :: this
        @:unused_var(this)
        ! Nothing to do
    end subroutine


    subroutine init_PC_WeightedParticles_t(this, indexer, weights, rank_with_info, use_lookup, create_lookup)
        class(PC_Particles_t), intent(inout) :: this
        class(AlsoGUGA_SuperGroupIndexer_t), intent(in) :: indexer
        real(dp), intent(in) :: weights(:, :, :)
            !! `w(j, i, i_sg)`, The weight of picking `j` after picking `i` under supergroup `i_sg`
        integer, intent(in) :: rank_with_info
            !! The **intra-node** rank that contains the weights to be used
            !! all other arr of all other ranks are ignored (and can be allocated with size 0).
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = 'init_PC_WeightedParticles_t'

        integer, allocatable :: supergroups(:, :)
        integer :: i_sg, nBi

        nBi = this%indexer%n_spin_orbs()
        @:ASSERT(2 * nSpatOrbs == nBi)

        allocate(this%indexer, source=indexer)
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

        supergroups = this%indexer%get_supergroups()

        if (iProcIndex_intra == rank_with_info) then
            @:ASSERT(nSpatOrbs == size(weights, 1) .and. nSpatOrbs == size(weights, 2))
            @:ASSERT(size(supergroups, 2) == size(weights, 3))
        end if

        call this%I_sampler%shared_alloc(size(supergroups, 2), nBi, 'PC_particles_i')
        call this%J_sampler%shared_alloc([nBi, size(supergroups, 2)], nBi, 'PC_particles_j')
        block
            integer :: I
            real(dp) :: dummy(2)
            real(dp), allocatable :: i_weight(:), j_weight(:)
            dummy = 10000._dp
            if (iProcIndex_intra == rank_with_info) then
                allocate(i_weight(nBi), j_weight(nBi), source=0._dp)
            else
                allocate(i_weight(0), j_weight(0))
            end if
            do i_sg = 1, size(supergroups, 2)
                do I = 1, nBI
                    if (iProcIndex_intra == rank_with_info) then
                        J_weight(1::2) = weights(:, gtID(I), i_sg)
                        J_weight(2::2) = J_weight(1::2)
                        J_weight(I) = 0._dp
                        call this%J_sampler%setup_entry(I, i_sg, rank_with_info, J_weight)
                    else
                        call this%J_sampler%setup_entry(I, i_sg, rank_with_info, dummy)
                    end if
                end do

                if (iProcIndex_intra == rank_with_info) then
                    i_weight(1::2) = sum(weights(:, :, i_sg), dim=1)
                    i_weight(2::2) = i_weight(1::2)
                    call this%I_sampler%setup_entry(i_sg, rank_with_info, i_weight)
                else
                    call this%I_sampler%setup_entry(i_sg, rank_with_info, dummy)
                end if
            end do
        end block
    end subroutine

    subroutine finalize_PC_WeightedParticles_t(this)
        class(PC_Particles_t), intent(inout) :: this

        if (associated(this%indexer)) then
            ! Yes, we assume, that either all or none are allocated
            call this%I_sampler%finalize()
            call this%j_sampler%finalize()
            deallocate(this%indexer)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
        end if
    end subroutine


    subroutine draw_PC_FullyWeightedParticles_t(this, nI, ilutI, csf_i, i_sg, srcs, p)
        class(PC_FullyWeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i_sg
            !! The supergroup index
        integer, intent(out) :: srcs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) <= srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        character(*), parameter :: this_routine = 'draw'

        real(dp) :: renorm_first, p_first(2)
        real(dp) :: renorm_second(2), p_second(2)
        integer :: spin_srcs(2), dummy

        renorm_first = sum(this%I_sampler%get_prob(i_sg, nI))
        call this%I_sampler%constrained_sample(&
            i_sg, nI, ilutI, renorm_first, dummy, spin_srcs(1), p_first(1))
        srcs(1) = gtID(spin_srcs(1))

        renorm_second(1) = sum(this%J_sampler%get_prob(spin_srcs(1), i_sg, nI))
        ! Note that p(I | I) is automatically zero and cannot be drawn
        call this%J_sampler%constrained_sample(&
            spin_srcs(1), i_sg, nI, ilutI, renorm_second(1), dummy, spin_srcs(2), p_second(1))
        if (spin_srcs(2) == 0) then
            p = 1._dp; srcs(:) = 0
            return
        end if
        srcs(2) = gtID(spin_srcs(2))

        if (srcs(1) == srcs(2)) then
            @:ASSERT(csf_i%occ_int(srcs(1)) == 2)
            p = 2._dp * p_first(1) * p_second(1)
        else
            ! From now on we assume srcs(1) /= srcs(2)!

            ! We could have picked them the other way round.
            ! Account for that.
            ! The renormalization for the first electron is the same
            p_first(2) = this%I_sampler%constrained_getProb(&
                i_sg, nI, renorm_first, spin_srcs(2))
            renorm_second(2) = sum(this%J_sampler%get_prob(spin_srcs(2), i_sg, nI))
            p_second(2) = this%J_sampler%constrained_getProb(&
                spin_srcs(2), i_sg, nI, renorm_second(2), spin_srcs(1))

            p = sum(csf_i%occ_int(srcs) * p_first(:) * csf_i%occ_int(srcs(2:1:-1)) * p_second(:))

            if (srcs(1) > srcs(2)) call swap(srcs(1), srcs(2))
        end if

        @:ASSERT(p .isclose. this%get_pgen(nI, csf_i, i_sg, srcs(1), srcs(2)), srcs, nI, csf_i%stepvector)
    end subroutine

    pure function get_pgen_PC_FullyWeightedParticles_t(this, nI, csf_i, i_sg, i, j) result(p)
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
        type(CSF_Info_t), intent(in) :: csf_I
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: i, j
            !! The particles in spatial format
        real(dp) :: p
        routine_name("get_pgen_PC_FullyWeightedParticles_t")

        @:pure_ASSERT(csf_I%stepvector(i) /= 0 .and. csf_I%stepvector(j) /= 0)

        associate(i_u => 2 * i - 1, i_d => 2 * i, &
                  j_u => 2 * j - 1, j_d => 2 * j)

            if (i == j) then
                if (csf_i%occ_int(i) /= 2) then
                    p = 0._dp
                    return
                else
                    p = 2 * p_first(i_u) * p_second(i_u, i_d)
                end if
            else
                p = &
                        p_first(i_u) * p_second(i_u, j_u) &
                        + p_first(i_u) * p_second(i_u, j_d) &
                        + p_first(i_d) * p_second(i_d, j_u) &
                        + p_first(i_d) * p_second(i_d, j_d) &
                        + p_first(j_u) * p_second(j_u, i_u) &
                        + p_first(j_u) * p_second(j_u, i_d) &
                        + p_first(j_d) * p_second(j_d, i_u) &
                        + p_first(j_d) * p_second(j_d, i_d)
            end if
        end associate

        contains

            real(dp) pure function p_first(I)
                integer, intent(in) :: I
                    !! Spin orbital index, i.e. u-d-index
                ! TODO replace with Ilut lookup
                if (I .in. nI) then
                    p_first = this%i_sampler%constrained_getProb(&
                        i_sg, nI, sum(this%I_sampler%get_prob(i_sg, nI)), I)
                else
                    p_first = 0._dp
                end if
            end function

            real(dp) pure function p_second(I, J)
                integer, intent(in) :: I, J
                    !! Spin orbital index, i.e. u-d-index
                ! TODO replace with Ilut lookup
                if ((I .in. nI) .and. (J .in. nI)) then
                    p_second = this%j_sampler%constrained_getProb(&
                                I, i_sg, nI, sum(this%j_sampler%get_prob(I, i_sg, nI)), J)
                else
                    p_second = 0._dp
                end if
            end function
    end function



    subroutine draw_PC_WeightedParticles_t(this, nI, ilutI, csf_i, i_sg, srcs, p)
        class(PC_WeightedParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i_sg
            !! The supergroup index
        integer, intent(out) :: srcs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) <= srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.

        integer :: spin_srcs(2), dummy
        real(dp) :: r_dummy, renorm


        spin_srcs(1) = nI(int(genrand_real2_dSFMT() * nEl) + 1)
        srcs(1) = gtID(spin_srcs(1))

        renorm = sum(this%J_sampler%get_prob(spin_srcs(1), i_sg, nI))
        ! Note that p(I | I) is automatically zero and cannot be drawn
        call this%J_sampler%constrained_sample(&
            spin_srcs(1), i_sg, nI, ilutI, renorm, dummy, spin_srcs(2), r_dummy)
        if (spin_srcs(2) == 0) then
            srcs(:) = 0; spin_srcs(:) = 0; p = 1._dp; srcs(:) = 0
            return
        end if
        srcs(2) = gtID(spin_srcs(2))

        if (srcs(1) > srcs(2)) call swap(srcs(1), srcs(2))
        p = this%get_pgen(nI, csf_i, i_sg, srcs(1), srcs(2))
    end subroutine


    pure function get_pgen_PC_WeightedParticles_t(this, nI, csf_i, i_sg, i, j) result(p)
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
        type(CSF_Info_t), intent(in) :: csf_I
        integer, intent(in) :: nI(nEl), i_sg
            !! The determinant in nI-format and the supergroup index
        integer, intent(in) :: i, j
            !! The particles in spatial format
        real(dp) :: p
        routine_name("get_pgen_PC_WeightedParticles_t")

        @:pure_ASSERT(csf_I%stepvector(i) /= 0 .and. csf_I%stepvector(j) /= 0)

        associate(i_u => 2 * i - 1, i_d => 2 * i, &
                  j_u => 2 * j - 1, j_d => 2 * j)

            if (i == j) then
                if (csf_i%occ_int(i) /= 2) then
                    p = 0_dp
                    return
                else
                    p = 2 * p_first(i_u) * p_second(i_u, i_d)

                    @:pure_ASSERT(near_zero(p_second(i_u, i_u)) &
                            .and. near_zero(p_second(i_d, i_d)))

                    @:pure_ASSERT(p_first(i_u) * p_second(i_u, i_d) &
                            .isclose. p_first(i_d) * p_second(i_d, i_u))
                end if
            else
                p = &
                        p_first(i_u) * p_second(i_u, j_u) &
                        + p_first(i_u) * p_second(i_u, j_d) &
                        + p_first(i_d) * p_second(i_d, j_u) &
                        + p_first(i_d) * p_second(i_d, j_d) &
                        + p_first(j_u) * p_second(j_u, i_u) &
                        + p_first(j_u) * p_second(j_u, i_d) &
                        + p_first(j_d) * p_second(j_d, i_u) &
                        + p_first(j_d) * p_second(j_d, i_d)
            end if
        end associate

        contains

            real(dp) pure function p_first(I)
                integer, intent(in) :: I
                    !! Spin orbital index, i.e. u-d-index
                if (I .in. nI) then
                    p_first = 1._dp / nEl
                else
                    p_first = 0._dp
                end if
            end function

            real(dp) pure function p_second(I, J)
                integer, intent(in) :: I, J
                    !! Spin orbital index, i.e. u-d-index
                if ((I .in. nI) .and. (J .in. nI)) then
                    p_second = this%j_sampler%constrained_getProb(&
                                I, i_sg, nI, sum(this%j_sampler%get_prob(I, i_sg, nI)), J)
                else
                    p_second = 0._dp
                end if
            end function
    end function
end module
