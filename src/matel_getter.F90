#include "macros.h"

module matel_getter

    ! This module contains getters for matrix elements <D|H|D> and <D0|H|D>
    ! as stored in global_determinant_data.  Moved here from load_balance so
    ! to avoid circular dependencies.

    use constants, only: dp, n_int
    use FciMCData, only: ilutRef, GetOffDiagMatel_Time
    use SystemData, only: nel, tHPHF, tNoBrillouin, tGUGA, t_3_body_excits, &
        t_ueg_3_body, t_mol_3_body
    use DetBitOps, only: FindBitExcitLevel
    use bit_reps, only: decode_bit_det
    use bit_rep_data, only: NIfTot
    use determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use guga_procedure_pointers, only: calc_off_diag_guga_ref
    use guga_bitRepOps, only: CSF_Info_t
    use timing_neci, only: set_timer, halt_timer

    implicit none
    private
    public get_diagonal_matel, get_off_diagonal_matel

contains

    function get_diagonal_matel(nI, ilut) result(diagH)
        ! Get the diagonal element for a determinant nI with ilut representation ilut

        ! In:  nI        - The determinant to evaluate
        !      ilut      - Bit representation (only used with HPHF)
        ! Ret: diagH     - The diagonal matrix element
        implicit none
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp) :: diagH

        if (tHPHF) then
            diagH = hphf_diag_helement(nI, ilut)
        else
            diagH = get_helement(nI, nI, 0)
        end if

    end function get_diagonal_matel

    function get_off_diagonal_matel(det, ilut) result(offdiagH)
        ! Get the off-diagonal element between the reference determinant
        ! a determinant det with ilut representation ilut, but only if it
        ! would contribute to the projected energy.

        ! In:  det       - The determinant to evaluate
        !      ilut      - Bit representation
        ! Ret: offdiagH  - The off-diagonal matrix element
        implicit none
        integer, intent(in) :: det(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        HElement_t(dp) :: offdiagH
        integer(n_int) :: ilut0(0:NIfTot)
        integer :: exlevel, det0(nel)

        call set_timer(GetOffDiagMatel_Time)

        ! We operate with the reference for run=1 here, which works for
        ! inum_runs==1 and inum_runs>1 with a common reference.  The projected
        ! energy for inum_runs>1 with different references is calculated
        ! separately anyway, so the return value does not get used.
        ilut0 = ilutRef(:,1)

        ! Get reference and excitation level with respect to it.
        call decode_bit_det(det0, ilut0)
        exlevel = FindBitExcitLevel(ilut0, ilut, t_hphf_ic=.true.)

        ! Compute matrix element depending on the calculation type.
        offdiagH = h_cast(0_dp)
        if (tGUGA) then
            ! TODO(@Oskar): Perhaps keep csf_i calculated?
            if (exlevel/=0) offdiagH = calc_off_diag_guga_ref(ilut, &
                CSF_Info_t(ilut), exlevel=exlevel)
        else
            if (  exlevel==2 .or. &
                ( exlevel==1 .and. tNoBrillouin ) .or. &
                ( exlevel==3 .and. (t_3_body_excits .or. &
                                    t_ueg_3_body .or. &
                                    t_mol_3_body) ) ) then
                if (tHPHF) then
                    offdiagH = hphf_off_diag_helement(det0, det, ilut0, ilut)
                else
                    offdiagH = get_helement(det, det0, exlevel, ilut, ilut0)
                end if
            end if
        end if

        call halt_timer(GetOffDiagMatel_Time)

    end function get_off_diagonal_matel

end module matel_getter
