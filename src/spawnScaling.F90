module spawnScaling
  ! This module allows to dynamically rescale the number of spawning attempts
  ! to prevent blooms -> increase the number of spawns to unbias
  ! This is done only for determinants attempting less than maxSpawnRescale spawns
  ! and has to be done on all spawns simultaneously: All spawns, including the ones
  ! already performed, have to be re-scaled
  use constants
  use bit_rep_data, only: NIfTot, test_flag, flag_rescale
  use fcimc_helper, only: decide_num_to_spawn
  use CalcData, only: initiatorWalkNo
  use FciMCData, only: tEScaleWalkers
  use procedure_pointers, only: scaleFunction
  use util_mod, only: near_zero
  implicit none
  private
  public :: resetScale, currentSpawnScale, rescaleIfFlagged, maxSpawnRescale, &
         scaleCondition, initSpawnScaling, finalizeSpawnScaling
  ! Fraction of an initiator that is ok to spawn in a single attempt
  real(dp) :: perSpawnInitRatio = 2.0_dp
  ! Maximum scale factor (everything above is deemed unstable)
  integer :: maxScaleFactor = 20
  ! scale factor
  integer :: spawnScale
  integer :: maxSpawnRescale = 1
  
contains

  !------------------------------------------------------------------------------------------!
  ! Interface function to update scaling factors
  !------------------------------------------------------------------------------------------!

  subroutine resetScale()
    implicit none

    ! and a factor of 1
    spawnScale = 1
  end subroutine resetScale

  !------------------------------------------------------------------------------------------!

  subroutine rescaleIfFlagged(ilut, WalkersToSpawn)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(inout) :: WalkersToSpawn
    integer, parameter :: defaultScale = 2
    ! check if this determinant has been problematic
    if(test_flag(ilut, flag_rescale) .and. (WalkersToSpawn <= maxSpawnRescale)) then
       ! rescale by default
       spawnScale = defaultScale
       ! add the extra spawn
       WalkersToSpawn = WalkersToSpawn + spawnScale - 1
    endif
  end subroutine rescaleIfFlagged

  !------------------------------------------------------------------------------------------!

  function currentSpawnScale() result(cScale)
    integer :: cScale

    cScale = spawnScale
  end function currentSpawnScale
  
  !------------------------------------------------------------------------------------------!

  ! determine whether a scaled bloom occurred on a given spawn
  subroutine scaleCondition(nSpawn, hdiag, factor)
    ! Determine if a new scaling factor for spawns shall be introduced
    ! If yes, also set the new scale factor and return the ratio
    ! Input: nSpawn - spawned weight
    !        hdiag - diagonal matrix element of the parent determinant
    !        spawnScale - array containing the scale factors per level
    !        scaleLevel - current scaling level
    !        factor - on return, the ratio between old and new scaling factor
    implicit none
    real(dp), intent(in) :: nSpawn
    real(dp), intent(in) :: hdiag
    integer, intent(out) :: factor
    real(dp) :: threshold
    logical :: tRescale

    ! maximum weight to be spawned in a single attempt
    threshold = real(InitiatorWalkNo,dp) / perSpawnInitRatio
    if(tEScaleWalkers) then
       threshold = threshold * scaleFunction(hdiag)
    endif
    ! check if we want to do a rescaling here
    ! do so if the scaled spawn would exceed the threshold and we can still scale
    ! no recursive rescaling!
    tRescale = (nSpawn > threshold) .and. (spawnScale == 1)

    if(tRescale .and. .not. near_zero(threshold)) then
       ! rescale this spawn down to the threshold
       ! this is essentially a stochastic round with some unused arguments
       call decide_num_to_spawn(nSpawn / threshold, hdiag, 1.0_dp, factor)
       factor = max(factor, 1)
       ! if the factor turned out to be one, no action is required
       if(factor > 1 .and. factor < maxScaleFactor) then
          spawnScale = factor * spawnScale
          ! increase the scaleLevel by 1 IF the new factor is not excessive
       else
          factor = 1
       endif
    else
       ! no scaling required
       factor = 1
    endif

  end subroutine scaleCondition


  !------------------------------------------------------------------------------------------!

  subroutine initSpawnScaling()
    implicit none

  end subroutine initSpawnScaling

  !------------------------------------------------------------------------------------------!

  subroutine finalizeSpawnScaling()
    implicit none

  end subroutine finalizeSpawnScaling

  
end module spawnScaling
