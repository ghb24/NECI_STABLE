module adi_initiators
  use adi_data, only: ilutRefAdi, nRefs, tAllDoubsInitiators, tAllSingsInitiators, tAdiActive
  use bit_rep_data, only: test_flag, NIfTot
  use bit_reps, only: set_flag, clr_flag
  use constants, only: lenof_sign, n_int, dp, inum_runs
  use SystemData, only: nel

contains

  !------------------------------------------------------------------------------------------!

  function check_static_init(ilut, nI, sgn, run) result(staticInit)
    use bit_rep_data, only: flag_adi_checked, flag_static_init
    use adi_data, only: tReferenceChanged
    implicit none
    integer(n_int), intent(inout) :: ilut(0:NIfTot)
    real(dp), intent(in) :: sgn(lenof_sign)
    integer, intent(in) :: run, nI(nel)
    logical :: staticInit

    staticInit = .false.
    ! Doubles are always initiators if the corresponding flag is set
    if(tAdiActive) then
       if(test_flag(ilut, flag_adi_checked)) then
          ! Check if the sign of the determinant changed, if yes, re-evaluate the adi flags
          ! only for this run
          if(check_sign_changed(ilut, sgn, run)) call set_adi_flags_run(ilut, nI, sgn, run)
       else
          ! else, evaluate all adi flags and then set flag_adi_checked
          call set_adi_flags(ilut, nI, sgn)
       endif
       ! If the adi check was not done yet, do it now for all runs and then
       ! set the corresponding flags
       ! By doing so, we ensure the check is done exactly one per run
       ! return the flag of the requested run
       staticInit = test_flag(ilut, flag_static_init(run))
    endif
  end function check_static_init

  !------------------------------------------------------------------------------------------!

  subroutine set_adi_flags(ilut, nI, sgn)
    use bit_rep_data, only: flag_adi_checked
    ! This sets the adi flags (flag_positive, flag_adi_checked and flag_static_init) 
    ! for a given ilut
    implicit none
    integer(n_int), intent(inout) :: ilut(0:NIfTot)
    real(dp), intent(in) :: sgn(lenof_sign)
    integer, intent(in) :: nI(nel)
    integer :: ir

    do ir = 1, inum_runs
       call set_adi_flags_run(ilut, nI, sgn, ir)
    end do

    ! Now, the adi criterium was checked
    call set_flag(ilut, flag_adi_checked)

  end subroutine set_adi_flags

  !------------------------------------------------------------------------------------------!

  subroutine set_adi_flags_run(ilut, nI, sgn, ir)
    use bit_rep_data, only: flag_static_init
    ! This sets the adi flags for a given run ir
    implicit none
    integer(n_int), intent(inout) :: ilut(0:NIfTot)
    real(dp), intent(in) :: sgn(lenof_sign)
    integer, intent(in) :: ir, nI(nel)
    logical :: staticInit

    ! Check if its a static initiator on run ir
    staticInit = adi_criterium(ilut, nI, sgn, ir)
    if(staticInit) then
       ! if so, set the flag
       call set_flag(ilut, flag_static_init(ir))
    else
       ! else, clear it
       call clr_flag(ilut, flag_static_init(ir))
    endif
    ! and set the flag_positive
    call assign_flag_positive(ilut, sgn, ir)

  end subroutine set_adi_flags_run

  !------------------------------------------------------------------------------------------!

  function adi_criterium(ilut, nI, sgn, run) result(staticInit)
    ! This is the adi-initiator criterium expansion
    ! I expect it to grow further
    use adi_references, only: giovannis_check
    use adi_data, only: tInitiatorsSubspace
    use DetBitOps, only: FindBitExcitLevel
    implicit none
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    real(dp), intent(in) :: sgn(lenof_sign)
    integer, intent(in) :: run, nI(nel)
    integer :: exLevel, i
    logical :: staticInit

    staticInit = .false.
    ! This is Giovanni's CAS-initiator criterium
    if(tInitiatorsSubspace) then
       if(giovannis_check(ilut)) then
          staticInit = .true.
       endif
    endif
    
    exLevel = 0
    ! Important : Only compare to the already initialized reference
    do i = 1, nRefs
       ! TODO: Check if i marks a superinitiator for the current run
       exLevel = FindBitExcitLevel(ilutRefAdi(:,i),ilut)
       ! We only need to do this if the excitation level is below 3
       if(exLevel < 3) then
          ! Check if we are sign-coherent if this is desired
          if(unset_incoherent_initiator(exLevel, ilut, nI, sgn, i, staticInit, run)) &
                                ! If we find the determinant to be incoherent, we do not apply
                                ! the ADI rules and instead 
               return

          ! Set the doubles to initiators
          call set_double_initiator(exLevel, staticInit)

          ! If desired, also set singles as initiators
          call set_single_initiator(exLevel, staticInit)
       endif
    enddo

  end function adi_criterium

!------------------------------------------------------------------------------------------!

  function unset_incoherent_initiator(exLevel, ilut, nI, sgn, iRef, &
       staticInit, run) result(tSuspendADI)
    use adi_data, only: tCoherentDoubles, nIncoherentDets
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
    if(tCoherentDoubles .and. (exLevel > 0)) then
       if(.not. check_sign_coherence(ilut,nI,sgn,iRef,run)) then
          ! If not, do not let the determinant be an initiator
          ! Note that this strikes anytime, even if it is coherent with 
          ! the reference to which it is the double, but not with some other one
          ! (if we have e.g. another reference that is a single)
          tSuspendADI = .true.
          staticInit = .false.
          nIncoherentDets = nIncoherentDets + 1
       endif
    endif
  end function unset_incoherent_initiator

!------------------------------------------------------------------------------------------!

  subroutine set_double_initiator(exLevel, staticInit)
    use adi_data, only:  nCoherentDoubles
    implicit none
    integer, intent(in) :: exLevel
    logical, intent(inout) :: staticInit

    if(exLevel == 2 .and. tAllDoubsInitiators) then
       staticInit = .true.
       ! also, log this event
       nCoherentDoubles = nCoherentDoubles + 1
    endif
  end subroutine set_double_initiator

  !------------------------------------------------------------------------------------------!

  subroutine set_single_initiator(exLevel, staticInit)
    use adi_data, only: nCoherentSingles
    implicit none
    integer, intent(in) :: exLevel
    logical, intent(out) :: staticInit

    if(exLevel == 1 .and. tAllSingsInitiators) then
       staticInit = .true.
       nCoherentSingles = nCoherentSingles + 1
    endif
  end subroutine set_single_initiator

  !------------------------------------------------------------------------------------------!

  function check_sign_changed(ilut, sgn, run) result(changed)
    use bit_rep_data, only: flag_positive
    implicit none
    integer(n_int), intent(inout) :: ilut(0:NIfTot)
    real(dp), intent(in) :: sgn(lenof_sign)
    integer, intent(in) :: run
    logical :: changed

#ifdef __CMPLX
    ! in the complex code, the "sign" i.e. the phase always changed as 
    ! soon as some walker is created
    ! I think the complex coherence check has to be approximated, else basically
    ! nothing will ever be coherent
    changed = .true.
#else
    changed = .true.
    if((test_flag(ilut, flag_positive) .eq. sgn(run) > 0)) changed = .false.
#endif
    ! If the sign changed, adjust the sign flag
    if(changed) call assign_flag_positive(ilut, sgn, run)
  end function check_sign_changed

!------------------------------------------------------------------------------------------!


  subroutine assign_flag_positive(ilut, sgn, run)
    use bit_rep_data, only: flag_positive
    implicit none
    integer(n_int), intent(inout) :: ilut(0:NIfTot)
    real(dp), intent(in) :: sgn(lenof_sign)
    integer, intent(in) :: run

    if(sgn(run) > 0) then 
       call set_flag(ilut, flag_positive)
    else
       call clr_flag(ilut, flag_positive)
    endif
  end subroutine assign_flag_positive

  !------------------------------------------------------------------------------------------!
end module adi_initiators
