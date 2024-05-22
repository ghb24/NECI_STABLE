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
! The generalization to GAS spaces is published in
!    https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00936
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and GAS constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in GAS/gasci_supergroup_index.fpp

module guga_gas_pchb_main
    !! Precomputed Heat Bath Implementation for GASCI. This modules implements
    !! the excitation generator GASCI PCHB either resolve in spin- or spatial-
    !! orbitals.

    use constants, only: stdout
    use util_mod, only: stop_all, EnumBase_t, operator(.implies.)
    use SystemData, only: nBasis
    use fortran_strings, only: split, to_upper
    use input_parser_mod, only: TokenIterator_t

    use generalized_supergroup, only: AlsoGUGA_PropertyIndexer_t
    use guga_base_class, only: AbInitGUGABase_t

    use guga_gas_pchb_singles_main, only: &
        GAS_PCHB_SinglesOptions_t, GAS_PCHB_SinglesOptions_vals_t, &
        singles_allocate_and_init => allocate_and_init

    use guga_gas_pchb_doubles_main, only: &
        PCHB_DoublesOptions_t, PCHB_DoublesOptions_vals_t, &
        doubles_allocate_and_init => allocate_and_init

    better_implicit_none
    private
    public :: GAS_PCHB_ExcGenerator_t, GAS_PCHB_options_t, options_vals, &
        GAS_PCHB_options_vals_t, from_tokens, options

    type, extends(AbInitGUGABase_t) :: GAS_PCHB_ExcGenerator_t
    contains
        procedure :: init
        procedure :: finalize
    end type

    type :: GAS_PCHB_options_t
        type(GAS_PCHB_SinglesOptions_t) :: singles
        type(PCHB_DoublesOptions_t) :: doubles
        logical :: use_lookup
    contains
        procedure :: assert_validity
    end type

    type :: GAS_PCHB_options_vals_t
        type(GAS_PCHB_SinglesOptions_vals_t) :: &
            singles = GAS_PCHB_SinglesOptions_vals_t()
        type(PCHB_DoublesOptions_vals_t) :: &
            doubles = PCHB_DoublesOptions_vals_t()
    end type

    type(GAS_PCHB_options_vals_t), parameter :: options_vals = GAS_PCHB_options_vals_t()
    type(GAS_PCHB_options_t) :: options

contains

    subroutine init(this, indexer, options)
        !! Initialize the PCHB excitation generator.
        !!
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: indexer
            !!  The GAS specifications for the excitation generator.
        type(GAS_PCHB_options_t), intent(in) :: options
        call options%assert_validity()

        call singles_allocate_and_init(indexer, options%singles, options%use_lookup, .false., this%singles_generator)
        call doubles_allocate_and_init(indexer, options%doubles, options%use_lookup, options%use_lookup, this%doubles_generator)
    end subroutine

    subroutine finalize(this)
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this

        if (allocated(this%singles_generator)) then
            ! Yes I assume that either all or none are allocated
            call this%singles_generator%finalize()
            call this%doubles_generator%finalize()
            deallocate(this%singles_generator, this%doubles_generator)
        end if
    end subroutine

    subroutine assert_validity(this)
        class(GAS_PCHB_options_t), intent(in) :: this
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
        type(GAS_PCHB_options_t) :: res
        routine_name("from_tokens")

        if (to_upper(tokens%next()) /= "DRAWING") then
            call stop_all(this_routine, "Speciy Drawing")
        end if
        res%singles%algorithm = options_vals%singles%algorithm%from_str(tokens%next())
        associate(splitted => split(tokens%next(), ':'))
            res%doubles%particle_selection = options_vals%doubles%particle_selection%from_str(splitted(1)%str)
            res%doubles%hole_selection = options_vals%doubles%hole_selection%from_str(splitted(2)%str)
        end associate
        if (to_upper(tokens%next()) /= "WEIGHTING") then
            call stop_all(this_routine, "Speciy WEIGHTING")
        end if
        res%singles%weighting = options_vals%singles%weighting%from_str(tokens%next())
        res%doubles%weighting = options_vals%doubles%weighting%from_str(tokens%next())

        ! To make sure we did not forget a component
        res = GAS_PCHB_options_t(res%singles, res%doubles, use_lookup)
    end function

end module
