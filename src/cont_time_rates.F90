#include "macros.h"
module cont_time_rates

    ! This module is part of the continuous time implementation in
    ! cont_time_fcimc.F90. It is split out for circuit breaking purposes
    ! (due to mutual dependencies of Annihilation.F90)

    use dSFMT_interface, only: genrand_real2_dSFMT
    use CalcData, only: tContTimeFull, DiagSft
    use Determinants, only: get_helement
    use SymExcit3, only: GenExcitations3
    use DetBitOps, only: EncodeBitDet
    use bit_rep_data, only: NIfTot
    use SystemData, only: nel
    use constants
    implicit none

contains

    function spawn_rate_full(det, ilut) result(rate)

        ! Calculate the spawning rate for continuous time spawning:
        !
        ! R_i = sum_j R_ij = sum_j abs(helement(i, j))
        !
        ! This in principle includes the self element, but we EXCLUDE that
        ! here, and add that term in on the fly in the continuous time loop
        ! (as it needs to be adjusted for the shift).

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: det(nel)
        real(dp) :: rate

        integer :: ex(2,2), flag, det_spwn(nel)
        logical :: par, found_all

        ASSERT(tContTimeFull)

        rate = 0
        found_all = .false.
        det_spwn = 0
        flag = 3
        ex = 0
        do while (.true.)

            ! Generate the next connected det
            call GenExcitations3(det, ilut, det_spwn, flag, ex, par, &
                                 found_all, .false.)
            if (found_all) &
                exit

            rate = rate + abs(get_helement(det, det_spwn))
        end do

    end function

    subroutine cont_time_gen_excit_full(det, ilut, rate, hdiag, det_spwn, &
                                        ilut_spwn, hoffdiag, ic, part_type)

        integer, intent(in) :: det(nel), part_type
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: rate, hdiag
        integer, intent(out) :: det_spwn(nel)
        integer(n_int), intent(out) :: ilut_spwn(0:NIfTot)
        HElement_t, intent(out) :: hoffdiag
        integer, intent(out) :: ic
        character(*), parameter :: this_routine = 'cont_time_gen_excit_full'

        integer :: ex(2,2), flag
        logical :: found_all, par
        real(dp) :: r, cum_rate

        ! Generate a random number to use throughout
        ! We wish to select the first det such that the cumulative sum of
        ! the rates becomes larger than this.
        r = genrand_real2_dSFMT() * rate

        ! Test the diagonal element first (as it will be by far the largest
        ! term in the cumulative sum). If we are generating this, only worry
        ! about the IC value, which is tested by the outer routine.
        cum_rate = abs(hdiag - DiagSft(part_type))
        if (cum_rate >= r) then
            det_spwn = det
            ic = 0
            return
        end if

        ! Use the slow, complete, excitation generator to generate all allowed
        ! excits and sum them in.
        found_all = .false.
        det_spwn = 0
        flag = 3
        ex = 0
        do while (.true.)

            ! Generate the next connected det
            call GenExcitations3(det, ilut, det_spwn, flag, ex, par, &
                                 found_all, .false.)

            ! If we have got to this point without making a selection, then
            ! something has gone very wrong...
            ! (Don't stop all, as it could occur (Just) due to rounding errors
            ! when it would still be valid)
            if (found_all) then
                write(iout, '("Potential cumulative sum error: ")')
                write(iout, *) cum_rate, rate, r
                exit
            end if

            hoffdiag = get_helement(det, det_spwn)
            cum_rate = cum_rate + abs(hoffdiag)
            if (cum_rate >= r) exit
        end do

        ! Now we want to use this determinant
        ! n.b. We are returning a non-zero (but invalid) IC value. The rest
        !      of the code for the full scheme doesn't care.
        ilut_spwn = 0
        call EncodeBitDet(det_spwn, ilut_spwn)
        ic = 99

    end subroutine


end module
