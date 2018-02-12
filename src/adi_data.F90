module adi_data

  use iso_c_hack
  use constants
  use FciMCData, only: ll_node
  implicit none
  save
  
  ! Number of references for all-doubs-initiators and (important) number of references 
  ! currently to check
  integer :: nRefs, nTZero, maxNRefs, nRefsSings, nRefsDoubs
  ! References for the purpose of the ADI scheme
  integer(n_int), allocatable :: ilutRefAdi(:,:)
  ! Store the signs and determinants separately, so they dont need to be 
  ! reconstructed on each coherence check
  integer, allocatable :: nIRef(:,:), exLvlRef(:)
  real(dp), allocatable :: signsRef(:,:)
  integer :: nIncoherentDets, nCoherentDoubles, nCoherentSingles, &
       AllCoherentSingles, AllCoherentDoubles, AllIncoherentDets, htBlock
  type(ll_node), pointer :: SIHash(:)
  logical :: tReferenceChanged, tSetupSIs, tUseCaches
  ! Data for the update of nrefs
  logical :: tSingleSteps, tVariableNRef
  real(dp) :: lastAllNoatHF
  integer :: lastNRefs
  ! desired reference population and tolerance
  integer :: targetRefPop, targetRefPopTol, nRefUpdateInterval

  ! Flags for the alldoublesinitiators feature
  logical :: tAllDoubsInitiators, tDelayAllDoubsInits, tAllSingsInitiators, tDelayAllSingsInits
  logical :: tSetDelayAllSingsInits, tSetDelayAllDoubsInits, tDelayGetRefs
  integer :: allDoubsInitsDelay, nExProd, superInitiatorLevel, SIUpdateInterval
  logical :: tAdiActive, tStrictCoherentDoubles, tWeakCoherentDoubles, tAvCoherentDoubles
  ! Thresholds for xi and populations
  real(dp) :: NoTypeN, coherenceThreshold, SIThreshold
  logical :: tReadRefs, tProductReferences, tAccessibleDoubles, tAccessibleSingles, tSuppressSIOutput

  ! Additional variables for giovannis check
  integer(n_int), allocatable :: g_markers(:)
  integer :: g_markers_num
  logical :: tInitiatorsSubspace

  integer :: nExCheckFails, nExChecks, allNExCheckFails, allNExChecks

  !Minimum number of connections to SI in order for the sign-coherence parameter to be valid
  integer :: minSIConnect

end module adi_data
