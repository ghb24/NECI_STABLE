module guga_procedure_pointers

    implicit none

    abstract interface
        subroutine pickOrbitals_t(ilut, nI, excitInfo, pgen)
            use constants, only: dp, n_int
            use bit_reps, only: nifguga
            use guga_data, only: excitationInformation
            use SystemData, only: nel
            implicit none
            integer(n_int), intent(in) :: ilut(0:nifguga)
            integer, intent(in) :: nI(nel)
            type(excitationInformation), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine pickOrbitals_t

        subroutine calc_orbital_pgen_contr_t(ilut, occ_orbs, above_cpt, below_cpt)
            use constants, only: dp, n_int
            use bit_reps, only: nifguga
            implicit none
            integer(n_int), intent(in) :: ilut(0:nifguga)
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine calc_orbital_pgen_contr_t

        subroutine calc_mixed_contr_t(ilut, t, excitInfo, pgen, integral)
            use constants, only: dp, n_int
            use bit_reps, only: nifguga
            use guga_data, only: excitationInformation
            implicit none
            integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
            type(excitationInformation), intent(inout) :: excitInfo
            real(dp), intent(out) :: pgen
            HElement_t(dp), intent(out) :: integral
        end subroutine calc_mixed_contr_t

        subroutine calc_mixed_start_contr_t(ilut, t, excitInfo, branch_pgen, pgen, &
                integral)
            use constants, only: dp, n_int
            use bit_reps, only: nifguga
            use guga_data, only: excitationInformation
            implicit none
            integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
            real(dp), intent(inout) :: branch_pgen
            type(excitationInformation), intent(inout) :: excitInfo
            real(dp), intent(out) :: pgen
            HElement_t(dp), intent(out) :: integral
        end subroutine calc_mixed_start_contr_t

        ! maybe scrap all the below and only to one general one.
        ! for minus and plus functions:
        function general_weight_dummy(nSwitches, bVal, dat) result(weight)
            use guga_data, only: weight_data
            use constants, only: dp
            implicit none
            real(dp), intent(in) :: nSwitches, bval
            type(weight_data), intent(in) :: dat
            real(dp) :: weight
        end function general_weight_dummy

        ! for zero function:
        function general_weight_zero(negSwitches, posSwitches, bVal, dat) &
                result(weight)
            use guga_data, only: weight_data
            use constants, only: dp
            implicit none
            real(dp), intent(in) :: negSwitches, posSwitches, bVal
            type(weight_data), intent(in) :: dat
            real(dp) :: weight
        end function general_weight_zero

        subroutine pick_first_orbital_t(i, pgen, excit_lvl, excit_typ)
            use constants, only: dp
            integer, intent(out) :: i, excit_lvl, excit_typ
            real(dp), intent(out) :: pgen
        end subroutine pick_first_orbital_t

        function orb_pgen_contrib_type_t() result(orb_pgen)
            use constants, only: dp
            real(dp) :: orb_pgen
        end function orb_pgen_contrib_type_t

        function calc_off_diag_guga_t(ilut, run, exlevel) result(hel)
            use constants, only: n_int, dp
            use bit_reps, only: niftot
            implicit none
            integer(n_int), intent(in) :: ilut(0:niftot)
            integer, intent(in), optional :: run
            integer, intent(out), optional :: exlevel
            HElement_t(dp) :: hel
        end function calc_off_diag_guga_t
    end interface

    procedure(pickOrbitals_t), pointer :: pickOrbitals_single
    procedure(pickOrbitals_t), pointer :: pickOrbitals_double
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


end module
