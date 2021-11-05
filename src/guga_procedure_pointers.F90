module guga_procedure_pointers

    use constants, only: dp, n_int, int_rdm

    use guga_data, only: ExcitationInformation_t, WeightData_t

    use guga_bitRepOps, only: CSF_Info_t

    use bit_rep_data, only: GugaBits, IlutBits

    use SystemData, only: nel


    implicit none

    private

    public :: general_weight_dummy, general_weight_zero, pickorbitals_single, &
              pickorbitals_double, calc_orbital_pgen_contr, &
              pick_first_orbital, orb_pgen_contrib_type_3, orb_pgen_contrib_type_2, &
              calc_off_diag_guga_ref, &
              calc_orbital_pgen_contrib_start, calc_orbital_pgen_contrib_end

    abstract interface
        subroutine PickOrbitals_t(ilut, nI, csf_i, excitInfo, pgen)
            import :: dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine PickOrbitals_t

        subroutine CalcOrbitalPgenContr_t(csf_i, occ_orbs, orb_a, orb_pgen)
            import :: dp, CSF_Info_t
            implicit none
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen
        end subroutine CalcOrbitalPgenContr_t

        subroutine calc_orbital_pgen_contr_t(csf_i, occ_orbs, above_cpt, below_cpt)
            import :: dp, CSF_Info_t
            implicit none
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine calc_orbital_pgen_contr_t

        subroutine calc_mixed_contr_t(ilut, t, csf_i, excitInfo, pgen, integral)
            import :: dp, n_int, CSF_Info_t, GugaBits, ExcitationInformation_t
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot), &
                                             t(0:GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(inout) :: excitInfo
            real(dp), intent(out) :: pgen
            HElement_t(dp), intent(out) :: integral
        end subroutine calc_mixed_contr_t

        subroutine calc_mixed_start_contr_t(ilut, csf_i, t, excitInfo, branch_pgen, pgen, &
                                            integral, rdm_ind, rdm_mat)
            import :: dp, n_int, int_rdm, GugaBits, CSF_Info_t, ExcitationInformation_t
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot), &
                                             t(0:GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            real(dp), intent(inout) :: branch_pgen
            type(ExcitationInformation_t), intent(inout) :: excitInfo
            real(dp), intent(out) :: pgen
            HElement_t(dp), intent(out) :: integral
            integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
            real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        end subroutine calc_mixed_start_contr_t

        ! maybe scrap all the below and only to one general one.
        ! for minus and plus functions:
        function general_weight_dummy(nSwitches, bVal, dat) result(weight)
            import :: WeightData_t, dp
            implicit none
            real(dp), intent(in) :: nSwitches, bval
            type(WeightData_t), intent(in) :: dat
            real(dp) :: weight
        end function general_weight_dummy

        ! for zero function:
        function general_weight_zero(negSwitches, posSwitches, bVal, dat) &
            result(weight)
            import :: WeightData_t, dp
            implicit none
            real(dp), intent(in) :: negSwitches, posSwitches, bVal
            type(WeightData_t), intent(in) :: dat
            real(dp) :: weight
        end function general_weight_zero

        subroutine pick_first_orbital_t(i, pgen, excit_lvl, excit_typ)
            import :: dp
            integer, intent(out) :: i, excit_lvl, excit_typ
            real(dp), intent(out) :: pgen
        end subroutine pick_first_orbital_t

        function orb_pgen_contrib_type_t() result(orb_pgen)
            import :: dp
            real(dp) :: orb_pgen
        end function orb_pgen_contrib_type_t

        function calc_off_diag_guga_t(ilut, csf_i, run, exlevel) result(hel)
            import :: n_int, dp, IlutBits, CSF_Info_t
            implicit none
            integer(n_int), intent(in) :: ilut(0:IlutBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in), optional :: run
            integer, intent(out), optional :: exlevel
            HElement_t(dp) :: hel
        end function calc_off_diag_guga_t

    end interface

    procedure(PickOrbitals_t), pointer :: pickOrbitals_single => null()
    procedure(PickOrbitals_t), pointer :: pickOrbitals_double => null()
    procedure(calc_orbital_pgen_contr_t), pointer :: calc_orbital_pgen_contr => null()

    procedure(pick_first_orbital_t), pointer :: pick_first_orbital => null()
    procedure(orb_pgen_contrib_type_t), pointer :: orb_pgen_contrib_type_3 => null()
    procedure(orb_pgen_contrib_type_t), pointer :: orb_pgen_contrib_type_2 => null()

    procedure(calc_off_diag_guga_t), pointer :: calc_off_diag_guga_ref => null()

    procedure(CalcOrbitalPgenContr_t), pointer :: calc_orbital_pgen_contrib_start => null()
    procedure(CalcOrbitalPgenContr_t), pointer :: calc_orbital_pgen_contrib_end => null()

end module
