#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

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
! The generalization to PropVec spaces is available in a preprint (should be available in JCTC soon as well)
!    https://chemrxiv.org/engage/chemrxiv/article-details/61447e60b1d4a605d589af2e
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and PropVec constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in gasci_supergroup_index.f90

module property_vector_pchb_main
    !! Precomputed Heat Bath Implementation for PropVecCI. This modules implements
    !! the excitation generator PropVecCI PCHB either resolve in spin- or spatial-
    !! orbitals.

    use constants, only: stdout
    use util_mod, only: stop_all, EnumBase_t, operator(.implies.)
    use timing_neci, only: set_timer, halt_timer
    use FciMCData, only: PropVec_PCHB_init_time
    use SystemData, only: tUHF, nBasis
    use fortran_strings, only: Token_t, join

    use property_vector_index, only: PropertyIndexer_t
    use gasci_singles_main, only: &
        PropVec_PCHB_SinglesOptions_t, PropVec_PCHB_SinglesOptions_vals_t, PropVec_PCHB_singles_options_vals, &
        PC_WeightedSinglesOptions_t, &
        singles_allocate_and_init => allocate_and_init
    use property_vector_pchb_doubles_main, only: PCHB_DoublesOptions_t, &
        PCHB_DoublesOptions_vals_t, &
        doubles_allocate_and_init => allocate_and_init

    use excitation_generators, only: ClassicAbInitExcitationGenerator_t

    better_implicit_none


    private
    public :: PropVec_PCHB_ExcGenerator_t, &
        PropVec_PCHB_options_t, PropVec_PCHB_SinglesOptions_vals_t, PropVec_PCHB_options_vals, &
        PropVec_PCHB_OptionsUserInput_t, PropVec_PCHB_user_input_vals, PropVec_PCHB_user_input, &
        decide_on_PCHB_options, PCHB_OptionSelection_t, PCHB_OptionSelection_vals_t


    type :: PropVec_PCHB_options_t
        type(PropVec_PCHB_SinglesOptions_t) :: singles
        type(PCHB_DoublesOptions_t) :: doubles
        logical :: use_lookup = .false.
            !! Use and/or create/manage the supergroup lookup.
    contains
        procedure :: assert_validity
        procedure :: to_str
    end type

    type :: PropVec_PCHB_options_vals_t
        type(PropVec_PCHB_SinglesOptions_vals_t) :: singles = PropVec_PCHB_SinglesOptions_vals_t()
        type(PCHB_DoublesOptions_vals_t) :: doubles = PCHB_DoublesOptions_vals_t()
    end type

    type(PropVec_PCHB_options_vals_t), parameter :: PropVec_PCHB_options_vals = PropVec_PCHB_options_vals_t()

    type, extends(EnumBase_t) :: PCHB_OptionSelection_t
    end type

    type :: PCHB_OptionSelection_vals_t
        type(PCHB_OptionSelection_t) :: &
            LOCALISED = PCHB_OptionSelection_t(1), &
            DELOCALISED = PCHB_OptionSelection_t(2), &
            MANUAL = PCHB_OptionSelection_t(3)
    end type

    type :: PropVec_PCHB_OptionsUserInput_t
        type(PCHB_OptionSelection_t) :: option_selection
        type(PropVec_PCHB_options_t), allocatable :: options
    end type

    type :: PropVec_PCHB_OptionsUserInput_vals_t
        type(PCHB_OptionSelection_vals_t) :: option_selection = PCHB_OptionSelection_vals_t()
        type(PropVec_PCHB_options_vals_t) :: options = PropVec_PCHB_options_vals_t()
    end type

    type(PropVec_PCHB_OptionsUserInput_vals_t), parameter :: PropVec_PCHB_user_input_vals = PropVec_PCHB_OptionsUserInput_vals_t()

    type(PropVec_PCHB_OptionsUserInput_t), allocatable :: PropVec_PCHB_user_input

    type, extends(ClassicAbInitExcitationGenerator_t) :: PropVec_PCHB_ExcGenerator_t
    contains
        private
        procedure, public :: init => PropVec_PCHB_init
    end type


contains


    subroutine PropVec_PCHB_init(this, indexer, options)
        !! Initialize the PCHB excitation generator.
        !!
        class(PropVec_PCHB_ExcGenerator_t), intent(inout) :: this
        class(PropertyIndexer_t), intent(in) :: indexer
            !!  The PropVec specifications for the excitation generator.
        type(PropVec_PCHB_options_t), intent(in) :: options

        call set_timer(PropVec_PCHB_init_time)

        call options%assert_validity()

        write(stdout, '(A)') 'PropVec PCHB with ' // options%to_str() // ' is used'
        call singles_allocate_and_init(indexer, options%singles, options%use_lookup, this%singles_generator)
        call doubles_allocate_and_init(indexer, options%doubles, options%use_lookup, this%doubles_generator)

        call halt_timer(PropVec_PCHB_init_time)
    end subroutine PropVec_PCHB_init

    subroutine assert_validity(this)
        class(PropVec_PCHB_options_t), intent(in) :: this
        routine_name("assert_validity")

        associate(singles => PropVec_PCHB_options_vals%singles)
        if (.not. (this%singles%algorithm == singles%algorithm%PC_WEIGHTED &
               .implies. (this%singles%PC_weighted%drawing /= singles%PC_weighted%drawing%UNDEFINED))) then
            call stop_all(this_routine, "PC-WEIGHTED singles require valid PC_weighted options.")
        end if
        end associate

        if (.not. (tUHF .implies. this%doubles%spin_orb_resolved)) then
            call stop_all(this_routine, "Spin-resolved excitation generation requires spin-resolved hole generation.")
        end if
    end subroutine


    pure function decide_on_PCHB_options(PropVec_PCHB_user_input, loc_nBasis, loc_nEl, loc_tUHF) result(res)
        type(PropVec_PCHB_OptionsUserInput_t), intent(in) :: PropVec_PCHB_user_input
        integer, intent(in) :: loc_nBasis, loc_nEl
        logical, intent(in) :: loc_tUHF
        type(PropVec_PCHB_options_t) :: res
        routine_name("decide_on_PCHB_options")
        associate(hole_selection_2 => merge(PropVec_PCHB_options_vals%doubles%hole_selection%FAST_FAST, &
                                            PropVec_PCHB_options_vals%doubles%hole_selection%FULL_FULL, &
                                            loc_nBasis > 4 * loc_nEl))
        if (PropVec_PCHB_user_input%option_selection == PropVec_PCHB_user_input_vals%option_selection%LOCALISED) then
            associate(single_selection => merge(PropVec_PCHB_options_vals%singles%PC_WEIGHTED%drawing%UNIF_FAST, &
                                                PropVec_PCHB_options_vals%singles%PC_WEIGHTED%drawing%UNIF_FULL, &
                                                loc_nBasis > 4 * loc_nEl))
                res = PropVec_PCHB_options_t(&
                        PropVec_PCHB_SinglesOptions_t(&
                            PropVec_PCHB_options_vals%singles%algorithm%PC_weighted, &
                            PC_WeightedSinglesOptions_t(single_selection) &
                        ), &
                        PCHB_DoublesOptions_t( &
                            PropVec_PCHB_options_vals%doubles%particle_selection%UNIF_FULL, &
                            hole_selection_2 &
                        ), &
                        use_lookup=.true. &
                )
            end associate
        else if (PropVec_PCHB_user_input%option_selection == PropVec_PCHB_user_input_vals%option_selection%DELOCALISED) then
            res = PropVec_PCHB_options_t(&
                    PropVec_PCHB_SinglesOptions_t(&
                        PropVec_PCHB_options_vals%singles%algorithm%BITMASK_UNIFORM &
                    ), &
                    PCHB_DoublesOptions_t( &
                        PropVec_PCHB_options_vals%doubles%particle_selection%UNIF_UNIF, &
                        hole_selection_2 &
                    ), &
                    use_lookup=.true. &
            )
        else if (PropVec_PCHB_user_input%option_selection == PropVec_PCHB_user_input_vals%option_selection%MANUAL) then
            res = PropVec_PCHB_user_input%options
        else
            call stop_all(this_routine, "Should not be here.")
        end if
        end associate

        if (.not. allocated(res%doubles%spin_orb_resolved)) then
            res%doubles%spin_orb_resolved = &
                loc_tUHF &
                .or. (res%doubles%hole_selection == PropVec_PCHB_options_vals%doubles%hole_selection%FULL_FULL)
        end if
    end function

    function to_str(options) result(res)
        class(PropVec_PCHB_options_t), intent(in) :: options
        character(:), allocatable :: res
        res = join([Token_t(options%singles%to_str()), Token_t(options%doubles%to_str())], ' ')
    end function

end module property_vector_pchb_main
