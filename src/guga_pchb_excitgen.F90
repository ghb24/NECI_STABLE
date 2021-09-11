#include "macros.h"
module guga_pchb_excitgen
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel
    use guga_data, only: tNewDet, ExcitationInformation_t, gen_type, excit_type
    use guga_bitRepOps, only: CSF_Info_t
    use guga_pchb_class, only: GugaAliasSampler_t

    implicit none

    private
    public :: init_guga_pchb_excitgen, finalize_pchb_excitgen_guga
    public :: pick_orbitals_pure_uniform_singles, pick_orbitals_double_pchb
    public :: calc_orbital_pgen_contr_pchb, calc_orbital_pgen_contr_start_pchb, calc_orbital_pgen_contr_end_pchb
    public :: calc_pgen_guga_pchb

    type(GugaAliasSampler_t), allocatable :: pchb_sampler

contains

    subroutine init_guga_pchb_excitgen()
        allocate(pchb_sampler)
        call pchb_sampler%init()
    end subroutine

    subroutine finalize_pchb_excitgen_guga()
        call pchb_sampler%finalize()
        deallocate(pchb_sampler)
    end subroutine

    subroutine pick_orbitals_pure_uniform_singles(ilut, nI, csf_info, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_info
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call pchb_sampler%pick_orbitals_pure_uniform_singles(ilut, nI, csf_info, excitInfo, pgen)
    end subroutine

    subroutine pick_orbitals_double_pchb(ilut, nI, csf_info, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_info
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call pchb_sampler%pick_orbitals_double_pchb(ilut, nI, csf_info, excitInfo, pgen)
    end subroutine

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    pure subroutine calc_orbital_pgen_contr_pchb(csf_info, occ_orbs, cpt_a, cpt_b)
        type(CSF_Info_t), intent(in) :: csf_info
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b
        call pchb_sampler%calc_orbital_pgen_contr_pchb(csf_info, occ_orbs, cpt_a, cpt_b)
    end subroutine calc_orbital_pgen_contr_pchb


    ! i think it would be better if i 'just' reimplement:
    pure subroutine calc_orbital_pgen_contr_start_pchb(csf_info, occ_orbs, a, orb_pgen)
        type(CSF_Info_t), intent(in) :: csf_info
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call pchb_sampler%calc_orbital_pgen_contr_start_pchb(csf_info, occ_orbs, a, orb_pgen)
    end subroutine

    pure subroutine calc_orbital_pgen_contr_end_pchb(csf_info, occ_orbs, a, orb_pgen)
        type(CSF_Info_t), intent(in) :: csf_info
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call pchb_sampler%calc_orbital_pgen_contr_end_pchb(csf_info, occ_orbs, a, orb_pgen)
    end subroutine calc_orbital_pgen_contr_end_pchb

    function calc_pgen_guga_pchb(ilutI, csf_info, ilutJ, excitInfo) result(pgen)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_info
        type(ExcitationInformation_t), intent(in), optional :: excitInfo
        real(dp) :: pgen
        pgen = pchb_sampler%calc_pgen_guga_pchb(ilutI, csf_info, ilutJ, excitInfo)
    end function

end module guga_pchb_excitgen
