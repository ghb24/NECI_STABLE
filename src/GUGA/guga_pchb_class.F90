#include "macros.h"
module guga_pchb_class

    use aliasSampling, only: AliasSampler_1D_t
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm, bits_n_int
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel, G1, &
                          nBasis, nSpatOrbs, ElecPairs, &
                          t_analyze_pchb, t_old_pchb, t_exchange_pchb
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut, CSF_Info_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: near_zero, fuseIndex, binary_search_first_ge, &
                        get_free_unit, stop_all, operator(.isclose.)
    use guga_excitations, only: assign_excitinfo_values_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, get_guga_integral_contrib_spat, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo, init_singleWeight, &
                                checkCompatibility_single
    use guga_bitrepops, only: identify_excitation, encode_excit_info, extract_excit_info, &
                              contract_2_rdm_ind
    use bit_reps, only: decode_bit_det
    use shared_array, only: shared_array_int64_t, shared_array_real_t
    use MPI_wrapper, only: iProcIndex_intra, iprocindex
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use procedure_pointers, only: get_umat_el

    use guga_gas_pchb_mod, only: GAS_PCHB_t

    use gasci, only: LocalGASSpec_t

    use guga_base_class, only: GUGABaseClass_t


    better_implicit_none

    private
    public :: PCHB_t

    type, extends(GUGABaseClass_t) :: PCHB_t
        private

        type(GAS_PCHB_t) :: GAS_sampler
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

        procedure, public :: get_valid_orbs
    end type PCHB_t

contains

    subroutine init(this)
        class(PCHB_t), intent(inout) :: this
        call this%GAS_sampler%init(CAS_spec(nEl, nSpatOrbs), &
                                   use_lookup=.false., create_lookup=.false.)

    contains

        type(LocalGASSpec_t) pure function CAS_spec(n_el, n_spat_orbs)
            integer, intent(in) :: n_el, n_spat_orbs
            integer :: i
            ! It does not matter if we use local or cumulative GAS
            ! constraints
            CAS_spec = LocalGASSpec_t(n_min=[n_el], n_max=[n_el], spat_GAS_orbs=[(1, i = 1, n_spat_orbs)])
        end function
    end subroutine

    subroutine finalize(this)
        class(PCHB_t), intent(inout) :: this
        call this%GAS_sampler%finalize()
    end subroutine


    subroutine pickorbitals_single(this, ilut, nI, csf_i, excitInfo, pgen)
        class(PCHB_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%GAS_sampler%pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine pickorbitals_double(this, ilut, nI, csf_i, excitInfo, pgen)
        class(PCHB_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%GAS_sampler%pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine calc_orbital_pgen_contr_start(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%GAS_sampler%calc_orbital_pgen_contr_start(csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine calc_orbital_pgen_contr_end(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%GAS_sampler%calc_orbital_pgen_contr_end(csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine calc_orbital_pgen_contr(this, csf_i, occ_orbs, above_cpt, below_cpt)
        class(PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt
        call this%GAS_sampler%calc_orbital_pgen_contr(csf_i, occ_orbs, above_cpt, below_cpt)
    end subroutine


    function calc_pgen(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        pgen = this%GAS_sampler%calc_pgen(nI, ilutI, csf_i, excitInfo)
    end function

    function get_valid_orbs(this, csf_i, src_orb) result(valid_orbs)
        class(PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: src_orb
        integer, allocatable :: valid_orbs(:)
        valid_orbs = this%GAS_sampler%get_valid_orbs(csf_i, src_orb)
    end function get_valid_orbs

end module guga_pchb_class
