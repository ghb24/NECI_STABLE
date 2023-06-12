#include "macros.h"
module guga_base_class
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type
    use guga_bitRepOps, only: CSF_Info_t

    better_implicit_none

    private
    public :: GUGABaseClass_t

    public :: init_guga_class, finalize_guga_class
    public :: class_pick_singles, class_pick_doubles
    public :: class_calc_orb_pgen_contr, class_calc_orb_pgen_contr_str, class_calc_orb_pgen_contr_end
    public :: class_calc_pgen
    public :: guga_exc_generator

    type, abstract :: GUGABaseClass_t
    contains
        procedure(PickOrbitals_t), deferred :: pickOrbitals_single
        procedure(PickOrbitals_t), deferred :: pickOrbitals_double

        procedure(CalcOrbitalPgenContr_t), deferred :: calc_orbital_pgen_contr_start
        procedure(CalcOrbitalPgenContr_t), deferred :: calc_orbital_pgen_contr_end
        procedure(calc_orbital_pgen_contr_t), deferred :: calc_orbital_pgen_contr

        procedure(CalcPgen_t), deferred :: calc_pgen

        procedure(BoundFinalize_t), deferred :: finalize
    end type

    interface
        ! This is implemented in a submodule
        module subroutine init_guga_class()
        end subroutine
    end interface

    class(GUGABaseClass_t), allocatable :: guga_exc_generator

    abstract interface
        subroutine PickOrbitals_t(this, ilut, nI, csf_i, excitInfo, pgen)
            import :: GUGABaseClass_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(GUGABaseClass_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine PickOrbitals_t

        subroutine CalcOrbitalPgenContr_t(this, csf_i, occ_orbs, orb_a, orb_pgen)
            import :: GUGABaseClass_t, dp, CSF_Info_t
            implicit none
            class(GUGABaseClass_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen
        end subroutine CalcOrbitalPgenContr_t

        subroutine calc_orbital_pgen_contr_t(this, csf_i, occ_orbs, above_cpt, below_cpt)
            import :: GUGABaseClass_t, dp, CSF_Info_t
            implicit none
            class(GUGABaseClass_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine calc_orbital_pgen_contr_t


        function CalcPgen_t(this, nI, ilutI, csf_i, excitInfo) result(pgen)
            import :: GUGABaseClass_t, dp, CSF_Info_t, n_int, GUGABits, ExcitationInformation_t, nEl
            implicit none
            class(GUGABaseClass_t), intent(in) :: this
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(in) :: excitInfo
            real(dp) :: pgen
        end function

        subroutine BoundFinalize_t(this)
            import :: GUGABaseClass_t
            implicit none
            class(GUGABaseClass_t), intent(inout) :: this
        end subroutine
    end interface

contains

    subroutine finalize_guga_class()
        call guga_exc_generator%finalize()
        deallocate(guga_exc_generator)
    end subroutine

    subroutine class_pick_singles(ilut, nI, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call guga_exc_generator%pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine class_pick_doubles(ilut, nI, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call guga_exc_generator%pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine class_calc_orb_pgen_contr(csf_i, occ_orbs, cpt_a, cpt_b)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b
        call guga_exc_generator%calc_orbital_pgen_contr(csf_i, occ_orbs, cpt_a, cpt_b)
    end subroutine class_calc_orb_pgen_contr


    ! i think it would be better if i 'just' reimplement:
    subroutine class_calc_orb_pgen_contr_str(csf_i, occ_orbs, a, orb_pgen)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call guga_exc_generator%calc_orbital_pgen_contr_start(csf_i, occ_orbs, a, orb_pgen)
    end subroutine

    subroutine class_calc_orb_pgen_contr_end(csf_i, occ_orbs, a, orb_pgen)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call guga_exc_generator%calc_orbital_pgen_contr_end(csf_i, occ_orbs, a, orb_pgen)
    end subroutine class_calc_orb_pgen_contr_end

    function class_calc_pgen(nI, ilutI, csf_i, excitInfo) result(pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        pgen = guga_exc_generator%calc_pgen(nI, ilutI, csf_i, excitInfo)
    end function


end module
