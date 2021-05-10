module guga_procedure_pointers

    use constants, only: dp, n_int, int_rdm

    use guga_data, only: ExcitationInformation_t, WeightData_t

    use bit_rep_data, only: GugaBits, IlutBits

    use SystemData, only: nel


    implicit none

    private

    public :: general_weight_dummy, general_weight_zero, pickorbitals_single, &
              pickorbitals_double, calc_orbital_pgen_contr, calc_mixed_contr, &
              calc_mixed_start_r2l_contr, calc_mixed_end_l2r_contr, &
              calc_mixed_start_l2r_contr, calc_mixed_end_r2l_contr, &
              pick_first_orbital, orb_pgen_contrib_type_3, orb_pgen_contrib_type_2, &
              calc_off_diag_guga_ref, gen_single_excit_guga, gen_double_excit_guga, &
              calc_orbital_pgen_contrib_start, calc_orbital_pgen_contrib_end

    abstract interface
        subroutine PickOrbitals_t(ilut, nI, excitInfo, pgen)
            import :: dp, n_int, GugaBits, ExcitationInformation_t, nel
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine PickOrbitals_t

        subroutine CalcOrbitalPgenContr_t(occ_orbs, orb_a, orb_pgen)
            import :: dp
            implicit none
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen

        end subroutine CalcOrbitalPgenContr_t

        subroutine calc_orbital_pgen_contr_t(ilut, occ_orbs, above_cpt, below_cpt)
            import :: dp, n_int, GugaBits
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine calc_orbital_pgen_contr_t

        subroutine calc_mixed_contr_t(ilut, t, excitInfo, pgen, integral)
            import :: dp, n_int, GugaBits, ExcitationInformation_t
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot), &
                                             t(0:GugaBits%len_tot)
            type(ExcitationInformation_t), intent(inout) :: excitInfo
            real(dp), intent(out) :: pgen
            HElement_t(dp), intent(out) :: integral
        end subroutine calc_mixed_contr_t

        subroutine calc_mixed_start_contr_t(ilut, t, excitInfo, branch_pgen, pgen, &
                                            integral, rdm_ind, rdm_mat)
            import :: dp, n_int, int_rdm, GugaBits, ExcitationInformation_t
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot), &
                                             t(0:GugaBits%len_tot)
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

        function calc_off_diag_guga_t(ilut, run, exlevel) result(hel)
            import :: n_int, dp, IlutBits
            implicit none
            integer(n_int), intent(in) :: ilut(0:IlutBits%len_tot)
            integer, intent(in), optional :: run
            integer, intent(out), optional :: exlevel
            HElement_t(dp) :: hel
        end function calc_off_diag_guga_t

        subroutine CreateSingleExcitGUGA_t(ilut, nI, excitation, pgen)
            import :: n_int, dp, GugaBits, nel
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(out) :: excitation(0:GugaBits%len_tot)
            real(dp), intent(out) :: pgen
        end subroutine CreateSingleExcitGUGA_t

        subroutine CreateDoubleExcitGUGA_t(ilut, nI, excitation, pgen, excit_typ)
            import :: n_int, dp, GugaBits, nel
            implicit none
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(out) :: excitation(0:GugaBits%len_tot)
            real(dp), intent(out) :: pgen
            integer, intent(out) :: excit_typ(2)
        end subroutine CreateDoubleExcitGUGA_t

    end interface

    procedure(PickOrbitals_t), pointer :: pickOrbitals_single
    procedure(PickOrbitals_t), pointer :: pickOrbitals_double
    procedure(calc_orbital_pgen_contr_t), pointer :: calc_orbital_pgen_contr
    procedure(calc_mixed_contr_t), pointer :: calc_mixed_contr
    procedure(calc_mixed_start_contr_t), pointer :: calc_mixed_start_l2r_contr
    procedure(calc_mixed_start_contr_t), pointer :: calc_mixed_start_r2l_contr
    procedure(calc_mixed_start_contr_t), pointer :: calc_mixed_end_r2l_contr
    procedure(calc_mixed_start_contr_t), pointer :: calc_mixed_end_l2r_contr

    procedure(pick_first_orbital_t), pointer :: pick_first_orbital
    procedure(orb_pgen_contrib_type_t), pointer :: orb_pgen_contrib_type_3
    procedure(orb_pgen_contrib_type_t), pointer :: orb_pgen_contrib_type_2

    procedure(calc_off_diag_guga_t), pointer :: calc_off_diag_guga_ref

    procedure(CreateSingleExcitGUGA_t), pointer :: gen_single_excit_guga
    procedure(CreateDoubleExcitGUGA_t), pointer :: gen_double_excit_guga

    procedure(CalcOrbitalPgenContr_t), pointer :: calc_orbital_pgen_contrib_start
    procedure(CalcOrbitalPgenContr_t), pointer :: calc_orbital_pgen_contrib_end

end module
