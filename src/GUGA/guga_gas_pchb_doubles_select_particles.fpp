#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module guga_gas_pchb_doubles_select_particles
    use constants, only: dp, n_int
    use fortran_strings, only: to_upper
    use bit_rep_data, only: nIfD
    use aliasSampling, only: AliasSampler_1D_t, AliasSampler_2D_t
    use SystemData, only: nEl, ElecPairs
    use sets_mod, only: operator(.in.)
    use util_mod, only: stop_all, EnumBase_t
    use UMatCache, only: gtID
    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer

    use guga_bitrepops, only: CSF_Info_t
    use guga_excitations, only: pick_elec_pair_uniform_guga

    better_implicit_none
    private
    public :: ParticleSelector_t, &
        UniformParticles_t, PCHB_ParticleSelection_t, allocate_and_init, &
        PCHB_ParticleSelection_vals_t, PCHB_particle_selection_vals



    type, extends(EnumBase_t) :: PCHB_ParticleSelection_t
        private
        character(9) :: str
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
            UNIF_FULL = PCHB_ParticleSelection_t(3, 'UNIF-FULL'), &
                !! We draw \( \tilde{p}(I)|_{D_i} \) uniformly and then \( p(J | I)_{J \in D_i} \)
                !! The second distribution comes from the PCHB weighting scheme.
                !! We guarantee that \(I\) and \(J\) are occupied.
            UNIF_FAST = PCHB_ParticleSelection_t(4, 'UNIF-FAST')
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

        real(dp) pure function GetPgen_t(this, nI, csf_I, I, J)
            import :: dp, CSF_Info_t, ParticleSelector_t, nEl
            implicit none
            integer, intent(in) :: nI(nEl)
                !! The CSF in nI-format
            class(ParticleSelector_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_I
            integer, intent(in) :: I, J
                !! The particles.
        end function

        subroutine Draw_t(this, nI, ilutI, csf_i, srcs, p)
            import :: dp, CSF_Info_t, ParticleSelector_t, nEl, n_int, nIfD
            class(ParticleSelector_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
                !! The determinant in nI-format and the supergroup index
            integer(n_int), intent(in) :: ilutI(0 : nIfD)
                !! The determinant in bitmask format
            type(CSF_Info_t), intent(in) :: csf_I
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

    ! type, extends(ParticleSelector_t), abstract :: PC_Particles_t
    !     private
    !     type(AliasSampler_1D_t) :: I_sampler
    !         !! The shape is (n_supergroup *) -> number_of_spat_orbs * 2
    !     type(AliasSampler_2D_t) :: J_sampler
    !         !! The shape is (number_of_spat_orbs * 2, n_supergroup) -> number_of_spat_orbs * 2
    !
    !     class(GASSpec_t), allocatable :: GAS_spec
    !     ! This is only a pointer because components cannot be targets
    !     ! otherwise. :-(
    !     type(SuperGroupIndexer_t), pointer :: indexer => null()
    !     logical, public :: use_lookup = .false.
    !         !! Use a lookup for the supergroup index in global_det_data.
    !     logical, public :: create_lookup = .false.
    !         !! Create **and** manage! the supergroup index lookup in global_det_data.
    ! contains
    !     private
    !     procedure, public :: init => init_PC_WeightedParticles_t
    !     procedure, public :: finalize => finalize_PC_WeightedParticles_t
    ! end type
    !
    ! type, extends(PC_Particles_t) :: PC_FullyWeightedParticles_t
    ! contains
    !     private
    !     procedure, public :: draw => draw_PC_FullyWeightedParticles_t
    !     procedure, public :: get_pgen => get_pgen_PC_FullyWeightedParticles_t
    ! end type
    !
    ! type, extends(PC_Particles_t) :: PC_WeightedParticles_t
    ! contains
    !     private
    !     procedure, public :: draw => draw_PC_WeightedParticles_t
    !     procedure, public :: get_pgen => get_pgen_PC_WeightedParticles_t
    ! end type
    !
    ! type, extends(PC_Particles_t) :: PC_FastWeightedParticles_t
    ! contains
    !     private
    !     procedure, public :: draw => draw_PC_FastWeightedParticles_t
    !     procedure, public :: get_pgen => get_pgen_PC_FastWeightedParticles_t
    ! end type

contains

    pure function to_str(options) result(res)
        class(PCHB_ParticleSelection_t), intent(in) :: options
        character(9) :: res
        res = options%str
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
        case(PCHB_particle_selection_vals%UNIF_FAST%str)
            res = PCHB_particle_selection_vals%UNIF_FAST
        case default
            call stop_all(this_routine, trim(w)//" not a valid doubles particle selection scheme.")
        end select
    end function

    subroutine allocate_and_init(PCHB_particle_selection, particle_selector)
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        class(ParticleSelector_t), allocatable, intent(inout) :: particle_selector
        routine_name("gasci_pchb_doubles_select_particles::allocate_and_init")
        if (PCHB_particle_selection == PCHB_particle_selection_vals%UNIF_UNIF) then
            allocate(UniformParticles_t :: particle_selector)
        else
            call stop_all(this_routine, 'Invalid particle selection.')
        end if
    end subroutine

    subroutine draw_UniformParticles_t(this, nI, ilutI, csf_i, srcs, p)
        class(UniformParticles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
            !! The determinant in nI-format and the supergroup index
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
            !! The determinant in bitmask format
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(out) :: srcs(2)
            !! The chosen particles \(I, J\) and their index in `nI`.
            !! It is guaranteed that `scrs(1) < srcs(2)`.
        real(dp), intent(out) :: p
            !! The probability of drawing \( p(\{I, J\}) \Big |_{D_i} \).
            !! This is the probability of drawing two particles from
            !! a given determinant \(D_i\) regardless of order.
        integer :: sym_prod, sum_ml
        routine_name("draw_UniformParticles_t")
        @:unused_var(this, ilutI, csf_i)
        call pick_elec_pair_uniform_guga(nI, srcs, sym_prod, sum_ml, p)
        ! make a different picker, which does not bias towards doubly
        ! occupied orbitals here too!

        @:ASSERT(srcs(1) < srcs(2))
        srcs = gtID(srcs)
        associate(i => srcs(1), j => srcs(2))
            if (i /= j) then
                if (csf_i%stepvector(i) == 3) p = 2.0_dp * p
                if (csf_i%stepvector(j) == 3) p = 2.0_dp * p
            end if
        end associate
    end subroutine

    pure function get_pgen_UniformParticles_t(this, nI, csf_i, i, j) result(p)
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
        integer, intent(in) :: nI(nEl)
            !! The CSF in nI-format
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j
            !! The particles.
        real(dp) :: p
        @:unused_var(this, nI, csf_i, i, j)
        p = 1.0_dp / real(ElecPairs, dp)
        if (i /= j) then
            if (csf_i%stepvector(i) == 3) p = 2.0_dp * p
            if (csf_i%stepvector(j) == 3) p = 2.0_dp * p
        end if
    end function

    subroutine finalize_UniformParticles_t(this)
        class(UniformParticles_t), intent(inout) :: this
        @:unused_var(this)
        ! Nothing to do
    end subroutine
end module
