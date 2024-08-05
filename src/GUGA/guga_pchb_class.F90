#include "macros.h"
module guga_pchb_class

    use constants, only: n_int, dp
    use bit_rep_data, only: GugaBits
    use SystemData, only: nel, nSpatOrbs
    use guga_data, only: ExcitationInformation_t
    use guga_bitrepops, only: CSF_Info_t
    use guga_base_class, only: GUGABase_t
    use gasci_supergroup_index, only: get_CAS_indexer
    use input_parser_mod, only: TokenIterator_t
    use guga_prop_vec_pchb_main, only: &
        PropVec_PCHB_ExcGenerator_t, PropVec_PCHB_options_t, &
        PropVec_PCHB_options_vals_t, PropVec_options_from_tokens => from_tokens, &
        PropVec_PCHB_OptionsUserInput_t, PropVec_PCHB_OptionsUserInput_vals_t, &
        PropVec_decide_on_PCHB_options => decide_on_PCHB_options


    better_implicit_none

    private
    public :: PCHB_t, FCI_PCHB_options_t, FCI_PCHB_options_vals_t, options_vals, &
        from_tokens, decide_on_PCHB_options, FCI_PCHB_user_input

    type, extends(GUGABase_t) :: PCHB_t
        private

        type(PropVec_PCHB_ExcGenerator_t) :: PropVec_sampler
    contains
        private
        procedure, public :: init
        procedure, public :: finalize

        procedure, public :: pickorbitals_single
        procedure, public :: pickorbitals_double
        procedure, public :: calc_pgen
        procedure, public :: calc_orbital_pgen_contr
        procedure, public :: calc_orbital_pgen_contr_start
        procedure, public :: calc_orbital_pgen_contr_end
    end type PCHB_t

    type :: FCI_PCHB_options_t
        type(PropVec_PCHB_options_t) :: prop_vec_pchb
    end type

    type :: FCI_PCHB_options_vals_t
        type(PropVec_PCHB_options_vals_t) :: prop_vec_pchb = PropVec_PCHB_options_vals_t()
    end type

    type :: FCI_PCHB_OptionsUserInput_t
        type(PropVec_PCHB_OptionsUserInput_t) :: prop_vec_pchb
    end type

    type :: FCI_PCHB_OptionsUserInput_vals_t
        type(PropVec_PCHB_OptionsUserInput_vals_t) :: &
            prop_vec_pchb = PropVec_PCHB_OptionsUserInput_vals_t()
    end type

    type(FCI_PCHB_OptionsUserInput_t), allocatable :: FCI_PCHB_user_input


    type(FCI_PCHB_options_vals_t), parameter :: options_vals = FCI_PCHB_options_vals_t()

contains

    subroutine init(this, options)
        class(PCHB_t), intent(inout) :: this
        type(FCI_PCHB_options_t), intent(in) :: options
        call this%PropVec_sampler%init(get_CAS_indexer(nEl, nSpatOrbs), &
            options%prop_vec_pchb)
    end subroutine

    subroutine finalize(this)
        class(PCHB_t), intent(inout) :: this
        call this%PropVec_sampler%finalize()
    end subroutine


    subroutine pickorbitals_single(this, nI, ilut, csf_i, excitInfo, pgen)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%PropVec_sampler%pickOrbitals_single(nI, ilut, csf_i, excitInfo, pgen)
    end subroutine

    subroutine pickorbitals_double(this, nI, ilut, csf_i, excitInfo, pgen)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%PropVec_sampler%pickOrbitals_double(nI, ilut, csf_i, excitInfo, pgen)
    end subroutine

    subroutine calc_orbital_pgen_contr_start(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%PropVec_sampler%calc_orbital_pgen_contr_start(nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine calc_orbital_pgen_contr_end(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%PropVec_sampler%calc_orbital_pgen_contr_end(nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine calc_orbital_pgen_contr(this, nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt
        call this%PropVec_sampler%calc_orbital_pgen_contr(nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
    end subroutine


    function calc_pgen(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        pgen = this%PropVec_sampler%calc_pgen(nI, ilutI, csf_i, excitInfo)
    end function

    function from_tokens(tokens, use_lookup) result(res)
        type(TokenIterator_t), intent(inout) :: tokens
        logical, intent(in) :: use_lookup
            !! Do a lookup of supergroups
        type(FCI_PCHB_OptionsUserInput_t) :: res
        res = FCI_PCHB_OptionsUserInput_t(PropVec_options_from_tokens(tokens, use_lookup))
    end function


    pure function decide_on_PCHB_options(FCI_PCHB_user_input, loc_nBasis, loc_nEl) result(res)
        type(FCI_PCHB_OptionsUserInput_t), intent(in) :: FCI_PCHB_user_input
        integer, intent(in) :: loc_nBasis, loc_nEl
        type(FCI_PCHB_options_t) :: res
        res = FCI_PCHB_options_t(PropVec_decide_on_PCHB_options(&
            FCI_PCHB_user_input%prop_vec_pchb, loc_nBasis, loc_nEl))
    end function

end module guga_pchb_class
