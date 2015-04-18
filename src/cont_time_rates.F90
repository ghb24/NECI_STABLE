#include "macros.h"
module cont_time_rates

    ! This module is part of the continuous time implementation in
    ! cont_time_fcimc.F90. It is split out for circuit breaking purposes
    ! (due to mutual dependencies of Annihilation.F90)

    use procedure_pointers, only: generate_excitation, encode_child, &
                                  get_spawn_helement
    use CalcData, only: tContTimeFull, DiagSft, cont_time_max_overspawn
    use Determinants, only: get_helement, write_det
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: excit_gen_store_type
    use SymExcit3, only: GenExcitations3
    use MemoryManager, only: TagIntType
    use LoggingData, only: FCIMCDebug
    use DetBitOps, only: EncodeBitDet
    use bit_rep_data, only: NIfTot
    use SystemData, only: nel, LMS
    use constants
    use util_mod
    implicit none
    save

    type(excit_gen_store_type) :: secondary_gen_store
    real(dp), allocatable :: oversample_factors(:,:)
    integer(TagIntType) :: ostag

    integer :: cont_spawn_attempts
    integer :: cont_spawn_success

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
        character(*), parameter :: this_routine = 'spawn_rate_full'

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

    subroutine cont_time_gen_excit(det, ilut, rate, hdiag, det_spwn, &
                                   ilut_spwn, hoffdiag, ic, part_type, nopen, &
                                   nspawn, store)

        integer, intent(in) :: det(nel), part_type, nopen
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: rate, hdiag
        integer, intent(out) :: det_spwn(nel)
        integer(n_int), intent(out) :: ilut_spwn(0:NIfTot)
        HElement_t, intent(out) :: hoffdiag
        integer, intent(out) :: ic, nspawn
        type(excit_gen_store_type), intent(inout) :: store
        character(*), parameter :: this_routine = 'cont_time_gen_excit_full'

        real(dp) :: probs(2), pgen, old, new_fac
        real(dp) :: rate_diag, rate_offdiag, pdiag, r, pneeded, pkeep
        integer :: ex(2,2), i, j
        logical :: tParity
        HElement_t :: helgen

        ! By default we spawn one particle
        nspawn = 1

        ! Test the diagonal element first (as it will be by far the largest
        ! term in the cumulative sum). If we are generating this, only worry
        ! about the IC value, which is tested by the outer routine.
        rate_diag = abs(hdiag - DiagSft(part_type))
        pdiag = rate_diag / rate
        r = genrand_real2_dSFMT()
        if (pdiag >= r) then
            det_spwn = det
            ic = 0
            return
        end if

        ! Adjust the random number so that it is accross [0, 1)
        r = (r - pdiag) / (1.0_dp - pdiag)
        rate_offdiag = rate - rate_diag

        ! Obtain the singles/doubles probability from the oversampling
        pSingles = oversample_factors(1, nopen) / rate_offdiag
        pDoubles = oversample_factors(2, nopen) / rate_offdiag
        probs = (/ pSingles, pDoubles /)

        cont_spawn_attempts = cont_spawn_attempts + 1
        call generate_excitation(det, ilut, det_spwn, ilut_spwn, 3, ic, ex, &
                                 tParity, pgen, helgen, store)
        IFDEBUG(FCIMCDebug,3) then
            write(iout, '("SP att: ",f12.5)', advance='no') pgen
            call write_det(iout, det_spwn, .true.)
            call neci_flush(iout)
        end if

        if (.not. IsNullDet(det_spwn)) then

            ! Get the diagonal matrix element
            hoffdiag = get_spawn_helement(det, det_spwn, ilut, ilut_spwn, &
                                          ic, ex, tParity, helgen)

            ! What is the required acceptance rate?
            pneeded = abs(hoffdiag) / rate_offdiag

            ! Is our generation rate high enough to give the correct
            ! distribution, or do we need to increase our oversampling?
            pkeep = pneeded / pgen
            if (pkeep > cont_time_max_overspawn) then

                old = oversample_factors(ic, nopen)
                new_fac = abs(hoffdiag) * probs(ic) &
                        / (pgen * cont_time_max_overspawn)
                if (new_fac > old) then
                    oversample_factors(ic, nopen) = new_fac

                    write(iout, '("New oversampling factors")')
                    do j = 1, 2
                        write(iout, '("ic",i1," ")', advance='no') j
                        write(iout, *) (oversample_factors(j,i),i=LMS,nel,2)
                    end do
                end if
                nspawn = stochastic_round_r(cont_time_max_overspawn, r)

            else if (pneeded == 0) then

                ! We shouldn't be spawning this particle at all
                det_spwn(1) = 0
                nspawn = 0

            else

                ! We are sampling enough --> This is the normal path
                nspawn = stochastic_round_r(pkeep, r)

            end if


            ! If this is going to survive, then encode it!
            if (nspawn /= 0) &
                call encode_child(ilut, ilut_spwn, ic, ex)

        end if
        
        ! Keep track of spawning successes
        cont_spawn_success = cont_spawn_success + nspawn

    end subroutine


end module
