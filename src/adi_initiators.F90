#include "macros.h"
module adi_initiators
    use adi_data, only: ilutRefAdi, nRefs, tAllDoubsInitiators, tAdiActive, &
                        nExChecks, nExCheckFails
    use bit_rep_data, only: test_flag, NIfTot, extract_sign
    use bit_reps, only: set_flag, clr_flag, decode_bit_det
    use constants, only: lenof_sign, n_int, dp, inum_runs, eps
    use SystemData, only: nel, max_ex_level

contains

    !------------------------------------------------------------------------------------------!

    function check_static_init(ilut, nI, sgn, ex, run) result(staticInit)
        use bit_rep_data, only: flag_adi_checked, flag_static_init
        use adi_data, only: tReferenceChanged
        implicit none
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        real(dp), intent(in) :: sgn(lenof_sign)
        integer, intent(in) ::  nI(nel), ex, run
        logical :: staticInit

        staticInit = .false.
        ! Doubles are always initiators if the corresponding flag is set
        if (tAdiActive) then
            ! Only check for the adi flags if this was not done before, or
            ! if the superinitiator pool changed
            if (.not. test_flag(ilut, flag_adi_checked) .or. tReferenceChanged) &
                call set_adi_flags(ilut, nI, sgn, ex)
            ! Check if the sign of the determinant changed, if yes, re-evaluate the adi flags
            ! only for this run
            ! suspend the sign change check for now
!       else
!          if(check_sign_changed(ilut, sgn, run)) &
!               call set_adi_flags_run(ilut, nI, sgn, run)
!       end if
            ! If the adi check was not done yet, do it now for all runs and then
            ! set the corresponding flags
            ! By doing so, we ensure the check is done exactly one per run
            ! return the flag of the requested run
            staticInit = test_flag(ilut, flag_static_init(run))
        end if
    end function check_static_init

    !------------------------------------------------------------------------------------------!

    subroutine set_adi_flags(ilut, nI, sgn, ex)
        use bit_rep_data, only: flag_adi_checked
        ! This sets the adi flags (flag_adi_checked and flag_static_init)
        ! for a given ilut
        implicit none
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        real(dp), intent(in) :: sgn(lenof_sign)
        integer, intent(in) :: nI(nel), ex
        integer :: ir

        ! this part has to be done independently for each run, since the coherence
        ! criterium depends on the exact population
        do ir = 1, inum_runs
            call set_adi_flags_run(ilut, nI, sgn, ex, ir)
        end do

        ! Now, the adi criterium was checked
        call set_flag(ilut, flag_adi_checked)

    end subroutine set_adi_flags

    !------------------------------------------------------------------------------------------!

    subroutine set_adi_flags_run(ilut, nI, sgn, ex, ir)
        use bit_rep_data, only: flag_static_init
        use adi_data, only: nCoherentDoubles
        ! This sets the adi flags for a given run ir
        implicit none
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        real(dp), intent(in) :: sgn(lenof_sign)
        integer, intent(in) :: ir, nI(nel), ex
        logical :: staticInit

        ! Check if its a static initiator on run ir
        staticInit = adi_criterium(ilut, nI, sgn, ex, ir)
        if (staticInit) then
            ! if so, set the flag
            call set_flag(ilut, flag_static_init(ir))
            nCoherentDoubles = nCoherentDoubles + 1
        else
            ! else, clear it
            call clr_flag(ilut, flag_static_init(ir))
        end if

    end subroutine set_adi_flags_run

    !------------------------------------------------------------------------------------------!

    function adi_criterium(ilut, nI, sgn, ex, run) result(staticInit)
        ! This is the adi-initiator criterium expansion
        ! I expect it to grow further
        use adi_references, only: initialize_c_caches, update_coherence_check, &
                                  eval_coherence
        use adi_data, only: tWeakCoherentDoubles, tAvCoherentDoubles, &
                            tUseCaches, exLvlRef
        use DetBitOps, only: FindBitExcitLevel
        implicit none
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: sgn(lenof_sign)
        integer, intent(in) :: run, nI(nel), ex
        integer :: exLevel, i
        logical :: staticInit, tCCache, tCouplingPossible, tSI
        ! cache for the weak coherence check
        HElement_t(dp) :: signedCache
        real(dp) :: unsignedCache
        integer :: connections

        tSI = .false.
        staticInit = .false.
        tCCache = tWeakCoherentDoubles .or. tAvCoherentDoubles

        exLevel = 0
        if (tCCache) call initialize_c_caches(signedCache, unsignedCache, connections)
        ! Important : Only compare to the already initialized reference
        do i = 1, nRefs
            ! First, check if the excitation level differs by more than 2
            tCouplingPossible = .true.
            if (tUseCaches) then
                if (abs(ex - exLvlRef(i)) > max_ex_level) tCouplingPossible = .false.
            end if
            nExChecks = nExChecks + 1
            if (tCouplingPossible) then
                exLevel = FindBitExcitLevel(ilutRefAdi(:, i), ilut)
                ! We only need to do this if the excitation level is below 3
                if (exLevel < max_ex_level + 1) then
                    ! Check if we are sign-coherent if this is desired
                    if (unset_incoherent_initiator(exLevel, ilut, nI, sgn, i, staticInit, run)) &
                        ! If we find the determinant to be incoherent, we do not apply
                        ! the ADI rules and instead
                        return

                    if (tCCache) then
                        call update_coherence_check(ilut, nI, i, &
                                                    signedCache, unsignedCache, connections)
                        if (exLevel == 0) tSI = .true.
                    end if

                    ! Set the doubles to initiators
                    call set_si_initiator(exLevel, staticInit)
                    ! We have to keep looping over the SIs, even if we already have a double,
                    ! if we want to compute Xi. Else, we can exit if we found a valid double
                    if (staticInit .and. .not. tCCache) exit
                end if
            else
                nExCheckFails = nExCheckFails + 1
            end if
        end do

        ! superinitiators themselves are always static initiators -> no coherence check then
        if (tCCache .and. staticInit .and. .not. tSI) &
            call eval_coherence(signedCache, unsignedCache, sgn(run), connections, staticInit)

    end function adi_criterium

!------------------------------------------------------------------------------------------!

    function unset_incoherent_initiator(exLevel, ilut, nI, sgn, iRef, &
                                        staticInit, run) result(tSuspendADI)
        use adi_data, only: tStrictCoherentDoubles, nIncoherentDets
        use adi_references, only: check_sign_coherence
        ! This version of the coherence check now only prevents us from setting
        ! the initiator flag via ADI, but not by regular means. So we can't remove 'normal'
        ! initiators anymore and also not make deterministic space determinants non-initiators
        implicit none
        integer, intent(in) :: exLevel, iRef
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: sgn(lenof_sign)
        logical, intent(inout) :: staticInit
        integer, intent(in) :: run, nI(nel)
        logical :: tSuspendADI

        tSuspendADI = .false.
        if (tStrictCoherentDoubles .and. (exLevel > 0)) then
            if (.not. check_sign_coherence(ilut, nI, sgn, iRef, run)) then
                ! If not, do not let the determinant be an initiator
                ! Note that this strikes anytime, even if it is coherent with
                ! the reference to which it is the double, but not with some other one
                ! (if we have e.g. another reference that is a single)
                tSuspendADI = .true.
                staticInit = .false.
                nIncoherentDets = nIncoherentDets + 1
            end if
        end if
    end function unset_incoherent_initiator

!------------------------------------------------------------------------------------------!

    subroutine set_si_initiator(exLevel, staticInit)
        implicit none
        integer, intent(in) :: exLevel
        logical, intent(inout) :: staticInit

        if (exLevel <= 2 .and. tAllDoubsInitiators) then
            staticInit = .true.
        end if
    end subroutine set_si_initiator

    !------------------------------------------------------------------------------------------!
end module adi_initiators
