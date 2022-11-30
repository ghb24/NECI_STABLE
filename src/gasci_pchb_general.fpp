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
! NOTE: This is only relevant for RHF-type calculations (see gasci_pchb_rhf.fpp)
!
! The generalization to GAS spaces is available in a preprint (should be available in JCTC soon as well)
!    https://chemrxiv.org/engage/chemrxiv/article-details/61447e60b1d4a605d589af2e
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and GAS constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in gasci_supergroup_index.f90

module gasci_pchb_general
    !! Precomputed Heat Bath Implementation for GASCI. This modules implements
    !! the excitation generator which builds on either gasci_pchb and gasci_pchb_uhf
    !! depending on if we are working in spin or spatial orbitals.

    use constants, only: stdout
    use util_mod, only: stop_all
    use timing_neci, only: set_timer, halt_timer
    use FciMCData, only: GAS_PCHB_init_time

    use gasci, only: GASSpec_t
    use gasci_pc_select_particles, only: PCHB_particle_selections, PCHB_ParticleSelection_t
    use gasci_singles_main_mod, only: GAS_used_singles_t, possible_GAS_singles, &
        GAS_singles_PC_uniform_ExcGenerator_t, GAS_singles_DiscardingGenerator_t, &
        GAS_singles_heat_bath_ExcGen_t
    use gasci_singles_pc_weighted, only: PC_SinglesOptions_t, &
        Base_PC_Weighted_t, do_allocation, print_options, &
        possible_PC_singles_weighting, possible_PC_singles_drawing

    use excitation_generators, only: ClassicAbInitExcitationGenerator_t


    use gasci_pchb_rhf, only: GAS_doubles_RHF_PCHB_ExcGenerator_t
    better_implicit_none


    private
    public :: GAS_PCHB_ExcGenerator_t, GAS_PCHB_options_t, GAS_PCHB_options

    type :: GAS_PCHB_options_t
        type(PCHB_ParticleSelection_t) :: particle_selection
        type(GAS_used_singles_t) :: singles
        logical :: UHF
            !! Do a spin-projection resolved calculation.
        type(PC_SinglesOptions_t) :: PC_singles_options = PC_SinglesOptions_t(&
            possible_PC_singles_weighting%UNDEFINED, possible_PC_singles_drawing%UNDEFINED)
            !! Only relevant if `singles == possible_PCHB_singles%PC_WEIGHTED`
        logical :: use_lookup= .false.
            !! Use and/or create/manage the supergroup lookup.
    end type

    type(GAS_PCHB_options_t) :: GAS_PCHB_options = GAS_PCHB_options_t( &
        PCHB_particle_selections%PC_WEIGHTED_APPROX, &
        possible_GAS_singles%PC_WEIGHTED, &
        PC_singles_options=PC_SinglesOptions_t(&
            possible_PC_singles_weighting%H_AND_G_TERM_BOTH_ABS, &
            possible_PC_singles_drawing%APPROX &
        ), &
        UHF=.false., &
        use_lookup=.true. &
    )

    type, extends(ClassicAbInitExcitationGenerator_t) :: GAS_PCHB_ExcGenerator_t
    contains
        private
        procedure, public :: init => GAS_PCHB_init
    end type


contains



    subroutine GAS_PCHB_init(this, GAS_spec, options)
        !! Initialize the PCHB excitation generator.
        !!
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
            !!  The GAS specifications for the excitation generator.
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_PCHB_options_t), intent(in) :: options
        routine_name("GAS_PCHB_init")

        call set_timer(GAS_PCHB_init_time)

        if (options%singles == possible_GAS_singles%DISCARDING_UNIFORM) then
            write(stdout, *) 'GAS discarding singles activated'
            allocate(this%singles_generator, source=GAS_singles_DiscardingGenerator_t(GAS_spec))
        else if (options%singles == possible_GAS_singles%BITMASK_UNIFORM) then
            write(stdout, *) 'GAS precomputed singles activated'
            allocate(GAS_singles_PC_uniform_ExcGenerator_t :: this%singles_generator)
            select type(generator => this%singles_generator)
            type is(GAS_singles_PC_uniform_ExcGenerator_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, options%use_lookup, create_lookup=.false.)
            end select
        else if (options%singles == possible_GAS_singles%ON_FLY_HEAT_BATH) then
            write(stdout, *) 'GAS heat bath on the fly singles activated'
            allocate(this%singles_generator, source=GAS_singles_heat_bath_ExcGen_t(GAS_spec))
        else if (options%singles == possible_GAS_singles%PC_WEIGHTED) then
            call print_options(options%PC_singles_options, stdout)
            call do_allocation(this%singles_generator, options%PC_singles_options%drawing)
            select type(generator => this%singles_generator)
            class is(Base_PC_Weighted_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, options%PC_singles_options%weighting, &
                                    options%use_lookup, create_lookup=.false.)
            end select
        else
            call stop_all(this_routine, "Invalid choise for singles.")
        end if


        ! @jph at the moment only RHF -- implement UHF
        if (options%UHF) then
            call stop_all('gas init', 'UHF PCHB not yet implemented :(')
        else
            allocate(GAS_doubles_RHF_PCHB_ExcGenerator_t :: this%doubles_generator)
            select type(generator => this%doubles_generator)
            type is(GAS_doubles_RHF_PCHB_ExcGenerator_t)
                call generator%init(&
                    GAS_spec, options%use_lookup, options%use_lookup, &
                    options%particle_selection)
            end select
        end if
        call halt_timer(GAS_PCHB_init_time)
    end subroutine GAS_PCHB_init


end module gasci_pchb_general
