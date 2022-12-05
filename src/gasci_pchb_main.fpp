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

module gasci_pchb_main
    !! Precomputed Heat Bath Implementation for GASCI. This modules implements
    !! the excitation generator which builds on either gasci_pchb and gasci_pchb_uhf
    !! depending on if we are working in spin or spatial orbitals.

    use constants, only: stdout
    use util_mod, only: stop_all, EnumBase_t, operator(.implies.)
    use timing_neci, only: set_timer, halt_timer
    use FciMCData, only: GAS_PCHB_init_time
    use SystemData, only: tUHF

    use gasci, only: GASSpec_t
    use gasci_singles_main, only: &
        GAS_PCHB_SinglesOptions_t, GAS_PCHB_SinglesOptions_vals_t, GAS_PCHB_singles_options_vals, &
        PC_WeightedSinglesOptions_t, &
        singles_allocate_and_init => allocate_and_init
    use gasci_pchb_doubles_main, only: PCHB_DoublesOptions_t, &
        PCHB_DoublesOptions_vals_t, &
        doubles_allocate_and_init => allocate_and_init

    use excitation_generators, only: ClassicAbInitExcitationGenerator_t

    better_implicit_none


    private
    public :: GAS_PCHB_ExcGenerator_t, &
        GAS_PCHB_options_t, GAS_PCHB_SinglesOptions_vals_t, GAS_PCHB_options_vals, &
        GAS_PCHB_options


    type :: GAS_PCHB_options_t
        type(GAS_PCHB_SinglesOptions_t) :: singles
        type(PCHB_DoublesOptions_t) :: doubles
        logical :: UHF
            !! Do a spin-projection resolved calculation.
        logical :: use_lookup = .false.
            !! Use and/or create/manage the supergroup lookup.
    contains
        procedure :: assert_validity
    end type

    type :: GAS_PCHB_options_vals_t
        type(GAS_PCHB_SinglesOptions_vals_t) :: singles = GAS_PCHB_SinglesOptions_vals_t()
        type(PCHB_DoublesOptions_vals_t) :: doubles = PCHB_DoublesOptions_vals_t()
    end type

    type(GAS_PCHB_options_vals_t), parameter :: GAS_PCHB_options_vals = GAS_PCHB_options_vals_t()

    type(GAS_PCHB_options_t) :: GAS_PCHB_options = GAS_PCHB_options_t( &
        GAS_PCHB_SinglesOptions_t(&
            GAS_PCHB_options_vals%singles%algorithm%PC_weighted, &
            PC_WeightedSinglesOptions_t(&
                GAS_PCHB_options_vals%singles%PC_weighted%weighting%H_AND_G_TERM_BOTH_ABS, &
                GAS_PCHB_options_vals%singles%PC_weighted%drawing%APPROX &
            ) &
        ), &
        PCHB_DoublesOptions_t( &
            GAS_PCHB_options_vals%doubles%particle_selection%FAST_WEIGHTED, &
            GAS_PCHB_options_vals%doubles%hole_selection%RHF_FAST_WEIGHTED &
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

        call set_timer(GAS_PCHB_init_time)

        call options%assert_validity()

        call singles_allocate_and_init(GAS_spec, options%singles, options%use_lookup, this%singles_generator)
        call doubles_allocate_and_init(GAS_spec, options%doubles, options%use_lookup, this%doubles_generator)

        call halt_timer(GAS_PCHB_init_time)
    end subroutine GAS_PCHB_init

    subroutine assert_validity(this)
        class(GAS_PCHB_options_t), intent(in) :: this
        routine_name("assert_validity")

        associate(singles => GAS_PCHB_options_vals%singles)
        if (.not. (this%singles%algorithm == singles%algorithm%PC_WEIGHTED &
               .implies. (this%singles%PC_weighted%weighting /= singles%PC_weighted%weighting%UNDEFINED &
                        .and. this%singles%PC_weighted%drawing /= singles%PC_weighted%drawing%UNDEFINED))) then
            call stop_all(this_routine, "PC-WEIGHTED singles require valid PC_weighted options.")
        end if
        end associate

        if (.not. (tUHF .implies. this%UHF)) then
            call stop_all(this_routine, "UHF requires spin-resolved PCHB")
        end if

        associate(doubles => GAS_PCHB_options_vals%doubles)
        if (.not. (this%UHF &
                    .implies. (this%doubles%hole_selection == doubles%hole_selection%UHF_FAST_WEIGHTED &
                                .or. this%doubles%hole_selection == doubles%hole_selection%UHF_FULLY_WEIGHTED))) then
            call stop_all(this_routine, "Spin resolved excitation generation requires spin resolved hole generation.")
        end if
        end associate

    end subroutine


end module gasci_pchb_main
