module CCMCData
   use HElem
   implicit none
   save
   real*8   dT1SqCuml
   logical  tExactCluster  ! Go through all combinations of excitors to make all clusters
   logical  tExactSpawn    ! For each cluster, go through all connected dets, and spawn there
   integer  nSpawnings     ! The number of spawning events per cluster if not tExactSpawn
   logical  tCCMCFCI       ! Run CCMC code without excitation clusters, recovering the FCIMC result
   logical  tAmplitudes    ! Use real numbers to indicate the amplitudes rather than stochastically sampling
   real*8   dInitAmplitude ! Specify the initial amplitude for use in CCMC amplitude calculations.
   real*8   dProbSelNewExcitor !The probability that the cluster selection algorithm terminates after each addition of an excitor.
   LOGICAL  tSpawnProp     ! Set if we use spawning proportional to the cluster amplitude rather than equally

   LOGICAL  tCCBuffer      ! Buffer the CC Amplitudes - this is useful when there are many cluster selections which lead to the same collapsed det. It creates a combined amplitude of the det first and spawns from that.
   

!This contains information as to a chosen Cluster
TYPE Cluster 
   INTEGER, allocatable :: SelectedExcitors(:,:)      !(0:NIfTot,nEl)  !The excitors which make up this cluster
   INTEGER, allocatable :: SelectedExcitorIndices(:)  !(nEl)     !The indices in a list of excitors, of the excitors which make up this cluster
   INTEGER, allocatable :: iLutDetCurr(:)             !(0:NIfTot) The determinant made from collapsing this cluster in bit representation
   INTEGER, allocatable :: DetCurr(:)                 !(nEl) The determinant made from collapsing this cluster.
   INTEGER  iSize
   INTEGER  iSgn                                      !The sign of the determinant after collapsing the cluster
   INTEGER iExcitLevel                                !The excitation level of the resultant det

   REAL*8   dAbsAmplitude
! dAbsAmplitude is the product of the coefficients of the excitors with the relevant normalizations.
! i.e. abs ( N0  (tI/N0) (tJ/N0) ... )
   REAL*8   dSelectionProb
! dSelectionProb is the probability that the cluster was selected



!The following are old
   REAL*8   dProbNorm
! dProbNorm is the prob that a cluster in this level would've been chosen had they been equally weighted
   REAL*8   dClusterProb
!dClusterProb is the Probability of having chosen this cluster excitor (normalized such that <1/dClusterProb> = 1)
   REAL*8   dClusterNorm                           
!  dClusterNorm is the probability that this cluster was chosen, given the level had already been selected.  This includes multiple selections of the same excitor as well as combinations of excitors which produce a 0 sign.
END TYPE Cluster

TYPE ClustSelector 
   INTEGER iIndex
   LOGICAL tFull     !Set if we are to generate all possible clusters
   INTEGER iMaxSize  !The maximum size of a cluster
   INTEGER nSelects  !If we're stochastically sampling the cluster space, this is the number of samples we take
   REAL*8 dProbSelNewExcitor  !The probability that we quit at every stage of selecting a new excitor for a cluster  
   TYPE(Cluster) C

END TYPE ClustSelector

TYPE Spawner 
   LOGICAL tFull        !Set if we go through all possible spawnees sequentially
   INTEGER nSpawnings   !The number of spawning events to attempt if we are randomly spawning
   INTEGER iIndex       !The index of the current spawning event
   INTEGER iMaxExcitLevel !The Max level away from the reference that an excit canbe
!   INTEGER, allocatable :: nI(:)           !The det from which to spawn
!   INTEGER, allocatable :: iLutnI(:)       !The det from which to spawn
   LOGICAL bValid       !Set if a valid det was found
   Type(Cluster), pointer:: C
   INTEGER, allocatable :: Scratch1(:)
   INTEGER, allocatable :: Scratch2(:)
   INTEGER, allocatable :: nJ(:)           !The det which is spawned to.
   INTEGER, allocatable :: iLutnJ(:)       !The det from which to spawn
   INTEGER iExcitLevel                     !The excitation level of the resultant det from the composite cluster
   TYPE(HElement)       :: HIJ
   REAL*8               :: dProbSpawn      !Prob that we spawned here (including the number of spawning events)
   INTEGER              :: ExcitMat(2,2)   !Internal data corresponding to the excitation matrix of the last generated det.
   INTEGER              :: ExFlag
END TYPE Spawner 

TYPE CCTransitionLog
   REAL*8, allocatable :: dProbTransition(:,:,:,:) !(2,2,nClust,nClust)
   REAL*8, allocatable :: dProbClust(:,:)  !(2,nClust)
   REAL*8, allocatable :: dProbUniqClust(:,:)  !(2,-1:nClust)
   INTEGER MaxIndex
END TYPE CCTransitionLog

end module CCMCData
