module spawnScaling
  ! This module allows to dynamically rescale the number of spawning attempts
  ! to prevent blooms -> increase the number of spawns to unbias
  ! This is done recursively as follows: We start at a scaleLevel of 1, with a
  ! factor of 1, this is the normal mode. When the need to scale
  ! arises, we set the factor of level 2 and increase the scaleLevel to 2.
  ! After a number of spawning attempts at this level equal to the factor, we
  ! decrease the level again.
  use constants
  use fcimc_helper, only: decide_num_to_spawn
  use CalcData, only: initiatorWalkNo
  use FciMCData, only: tEScaleWalkers
  use procedure_pointers, only: scaleFunction
  use util_mod, only: near_zero
  implicit none
  private
  public :: setLevelTimer, resetScale, updateLevelTimer, currentSpawnScale, &
         scaleCondition, initSpawnScaling, finalizeSpawnScaling
  ! Fraction of an initiator that is ok to spawn in a single attempt
  real(dp) :: perSpawnInitRatio = 2.0_dp
  ! Maximum recursion depth for rescaling
  integer :: maxScaleLevel = 15
  ! Maximum scale factor (everything above is deemed unstable)
  integer :: maxScaleFactor = 20
  ! Timer to count the number of spawns per level
  integer, allocatable :: levelTimer(:)
  ! Target number of spawns per level
  integer, allocatable :: levelLifeTime(:)
  ! scale factor per level
  integer, allocatable :: spawnScale(:)
  ! The current scaleLevel
  integer :: scaleLevel
  
contains

  !------------------------------------------------------------------------------------------!
  ! Interface function to update scaling factors
  !------------------------------------------------------------------------------------------!

  subroutine resetScale()
    implicit none

    ! start at level 1
    scaleLevel = 1
    ! and a factor of 1
    spawnScale = 1
  end subroutine resetScale

  !------------------------------------------------------------------------------------------!

  subroutine setLevelTimer()
    implicit none

    ! this is the number of spawns attempted at this level
    levelTimer(scaleLevel) = 0
    ! and this is the total number of spawns we will attempt
    levelLifetime(scaleLevel) = spawnScale(scaleLevel)
  end subroutine setLevelTimer

  
  !------------------------------------------------------------------------------------------!

  subroutine updateLevelTimer()
    implicit none

    ! Do the scaleLevel bookkeeping:
    ! Only for levels larger than one, an action is required
    if(scaleLevel > 1) then
       ! add one spawn at the current level
       levelTimer(scaleLevel) = levelTimer(scaleLevel) + 1
       ! if we are at the max number of spawns at this level, decrease it
       if(levelTimer(scaleLevel) == levelLifeTime(scaleLevel)) then
          scaleLevel = scaleLevel - 1
       endif
    endif
    
  end subroutine updateLevelTimer

  !------------------------------------------------------------------------------------------!

  function currentSpawnScale() result(cScale)
    integer :: cScale

    cScale = spawnScale(scaleLeveL)
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
    tRescale = (nSpawn > threshold * spawnScale(scaleLevel)) &
         .and. (scaleLevel < maxScaleLevel)

    if(tRescale .and. .not. near_zero(threshold)) then
       ! rescale this spawn down to the threshold
       ! this is essentially a stochastic round with some unused arguments
       call decide_num_to_spawn(nSpawn / threshold, hdiag, 1.0_dp, factor)
       factor = max(factor, 1)
       ! if the factor turned out to be one, no action is required
       if(factor > 1) then
          spawnScale(scaleLevel+1) = spawnScale(scaleLevel) * factor
          ! increase the scaleLevel by 1 IF the new factor is not excessive
          if(spawnScale(scaleLevel+1) < maxScaleFactor) then
             scaleLevel = scaleLevel + 1
             if(scaleLevel > 2) then
                print *, "WARNING: Recursive rescaling required"
             endif
          else
             ! this would be excessive, do not scale
             factor = 1
          end if
       endif
    else
       ! no scaling required
       factor = 1
    endif

  end subroutine scaleCondition


  !------------------------------------------------------------------------------------------!

  subroutine initSpawnScaling()
    implicit none

    allocate(levelTimer(maxScaleLevel))
    allocate(levelLifeTime(maxScaleLevel))
    allocate(spawnScale(maxScaleLevel))
  end subroutine initSpawnScaling

  !------------------------------------------------------------------------------------------!

  subroutine finalizeSpawnScaling()
    implicit none

    if(allocated(levelTimer)) deallocate(levelTimer)
    if(allocated(levelLifeTime)) deallocate(levelLifeTime)
    if(allocated(spawnScale)) deallocate(spawnScale)
  end subroutine finalizeSpawnScaling

  
end module spawnScaling
