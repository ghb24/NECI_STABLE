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
! NOTE: This is only relevant for RHF-type calculations
!
! The generalization to PropVec spaces is published in
!    https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00936
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and PropVec constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in GAS/gasci_supergroup_index.fpp

module guga_prop_vec_pchb_main
    !! Precomputed Heat Bath Implementation for PropVecCI. This modules implements
    !! the excitation generator PropVecCI PCHB either resolve in spin- or spatial-
    !! orbitals.

    use constants, only: stdout, stderr
    use util_mod, only: stop_all, EnumBase_t, operator(.implies.)
    use SystemData, only: nBasis
    use fortran_strings, only: split, to_upper
    use input_parser_mod, only: TokenIterator_t

    use property_vector_index, only: AlsoGUGA_PropertyIndexer_t
    use guga_base_class, only: AbInitGUGABase_t

    use guga_prop_vec_pchb_singles_main, only: &
        PropVec_PCHB_SinglesOptions_t, PropVec_PCHB_SinglesOptions_vals_t, &
        singles_allocate_and_init => allocate_and_init

    use guga_prop_vec_pchb_doubles_main, only: &
        PCHB_DoublesOptions_t, PCHB_DoublesOptions_vals_t, &
        doubles_allocate_and_init => allocate_and_init

    better_implicit_none
    private
    public :: PropVec_PCHB_ExcGenerator_t, PropVec_PCHB_options_t, options_vals, &
        PropVec_PCHB_options_vals_t, from_tokens, &
        PropVec_PCHB_OptionsUserInput_t, PropVec_PCHB_OptionsUserInput_vals_t, &
        PropVec_PCHB_user_input, PropVec_PCHB_user_input_vals,  decide_on_PCHB_options

    type, extends(AbInitGUGABase_t) :: PropVec_PCHB_ExcGenerator_t
    contains
        procedure :: init
        procedure :: finalize
    end type

    type :: PropVec_PCHB_options_t
        type(PropVec_PCHB_SinglesOptions_t) :: singles
        type(PCHB_DoublesOptions_t) :: doubles
        logical :: use_lookup
    contains
        procedure :: assert_validity
    end type

    type :: PropVec_PCHB_options_vals_t
        type(PropVec_PCHB_SinglesOptions_vals_t) :: &
            singles = PropVec_PCHB_SinglesOptions_vals_t()
        type(PCHB_DoublesOptions_vals_t) :: &
            doubles = PCHB_DoublesOptions_vals_t()
    end type

    type(PropVec_PCHB_options_vals_t), parameter :: options_vals = PropVec_PCHB_options_vals_t()

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

contains

    subroutine init(this, indexer, options)
        !! Initialize the PCHB excitation generator.
        !!
        class(PropVec_PCHB_ExcGenerator_t), intent(inout) :: this
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: indexer
            !!  The PropVec specifications for the excitation generator.
        type(PropVec_PCHB_options_t), intent(in) :: options
        call options%assert_validity()

        call singles_allocate_and_init(indexer, options%singles, options%use_lookup, .false., this%singles_generator)
        call doubles_allocate_and_init(indexer, options%doubles, options%use_lookup, options%use_lookup, this%doubles_generator)
    end subroutine

    subroutine finalize(this)
        class(PropVec_PCHB_ExcGenerator_t), intent(inout) :: this

        if (allocated(this%singles_generator)) then
            ! Yes I assume that either all or none are allocated
            call this%singles_generator%finalize()
            call this%doubles_generator%finalize()
            deallocate(this%singles_generator, this%doubles_generator)
        end if
    end subroutine

    subroutine assert_validity(this)
        class(PropVec_PCHB_options_t), intent(in) :: this
        routine_name("assert_validity")
        if (.not. this%singles%is_valid()) then
            call stop_all(this_routine, 'singles options are not valid')
        end if
        if (.not. this%doubles%is_valid()) then
            call stop_all(this_routine, 'doubles options are not valid')
        end if
    end subroutine

    function from_tokens(tokens, use_lookup) result(res)
        type(TokenIterator_t), intent(inout) :: tokens
        logical, intent(in) :: use_lookup
        type(PropVec_PCHB_OptionsUserInput_t) :: res
        routine_name("from_tokens")

        if (tokens%remaining_items() == 0) then
            call stop_all(this_routine, "No item specified. Please specify one of {LOCALISED, DELOCALISED, MANUAL}.")
        end if
        select case (to_upper(tokens%next()))
        case ("LOCALISED")
            res%option_selection &
                = PropVec_PCHB_user_input_vals%option_selection%LOCALISED
        case ("DELOCALISED")
            res%option_selection &
                = PropVec_PCHB_user_input_vals%option_selection%DELOCALISED
        case ("MANUAL")
            res%option_selection &
                = PropVec_PCHB_user_input_vals%option_selection%MANUAL
            allocate(res%options)
            if (to_upper(tokens%next()) /= "DRAWING") then
                call stop_all(this_routine, "Speciy Drawing")
            end if
            res%options%singles%algorithm = options_vals%singles%algorithm%from_str(tokens%next())
            associate(splitted => split(tokens%next(), ':'))
                res%options%doubles%particle_selection = options_vals%doubles%particle_selection%from_str(splitted(1)%str)
                res%options%doubles%hole_selection = options_vals%doubles%hole_selection%from_str(splitted(2)%str)
            end associate
            if (to_upper(tokens%next()) /= "WEIGHTING") then
                call stop_all(this_routine, "Speciy WEIGHTING")
            end if
            res%options%singles%weighting = options_vals%singles%weighting%from_str(tokens%next())
            res%options%doubles%weighting = options_vals%doubles%weighting%from_str(tokens%next())

            ! To make sure we did not forget a component
            res%options = PropVec_PCHB_options_t(res%options%singles, res%options%doubles, use_lookup)
        case default
            write(stderr, *) "Please specify one of {LOCALISED, DELOCALISED, MANUAL}."
            write(stderr, *) "With localised and delocalised the fastest PCHB sampler is chosen automatically."
            write(stderr, *) "Specify delocalised for Hartree-Fock like orbitals and"
            write(stderr, *) "localised for localised orbitals."
            write(stderr, *) "With manual one can select the PCHB sampler themselves."
            call stop_all(this_routine, "Wrong input to PCHB.")
        end select


    end function


    pure function decide_on_PCHB_options(PropVec_PCHB_user_input, loc_nBasis, loc_nEl) result(res)
        type(PropVec_PCHB_OptionsUserInput_t), intent(in) :: PropVec_PCHB_user_input
        integer, intent(in) :: loc_nBasis, loc_nEl
        type(PropVec_PCHB_options_t) :: res
        routine_name("decide_on_PCHB_options")

        associate(hole_selection_2 => merge(options_vals%doubles%hole_selection%FAST_FAST, &
                                            options_vals%doubles%hole_selection%FULL_FULL, &
                                            loc_nBasis > 6 * loc_nEl))
        if (PropVec_PCHB_user_input%option_selection == PropVec_PCHB_user_input_vals%option_selection%LOCALISED) then
            res = PropVec_PCHB_options_t(&
                    PropVec_PCHB_SinglesOptions_t(&
                        options_vals%singles%algorithm%UNIF_FULL, &
                        options_vals%singles%weighting%RHF_EXPR &
                    ), &
                    PCHB_DoublesOptions_t( &
                        options_vals%doubles%particle_selection%UNIF_FULL, &
                        hole_selection_2, &
                        options_vals%doubles%weighting%INITIAL_DOBRAUTZ &
                    ), &
                    use_lookup=.true. &
            )
        else if (PropVec_PCHB_user_input%option_selection == PropVec_PCHB_user_input_vals%option_selection%DELOCALISED) then
            res = PropVec_PCHB_options_t(&
                    PropVec_PCHB_SinglesOptions_t(&
                        options_vals%singles%algorithm%UNIF_UNIF, &
                        options_vals%singles%weighting%UNIFORM &
                    ), &
                    PCHB_DoublesOptions_t( &
                        options_vals%doubles%particle_selection%UNIF_FULL, &
                        hole_selection_2, &
                        options_vals%doubles%weighting%INITIAL_DOBRAUTZ &
                    ), &
                    use_lookup=.true. &
            )
        else if (PropVec_PCHB_user_input%option_selection == PropVec_PCHB_user_input_vals%option_selection%MANUAL) then
            res = PropVec_PCHB_user_input%options
        else
            call stop_all(this_routine, "Should not be here.")
        end if
        end associate
    end function

end module
