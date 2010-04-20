MODULE FciMCData
      use constants, only: dp,int64,n_int
      USE global_utilities
      IMPLICIT NONE
      SAVE

      INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor

      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: WalkVecDets(:,:)                !Contains determinant list
      INTEGER , ALLOCATABLE , TARGET :: WalkVecSign(:)                    !Contains sign list (1 = positive, -1 = negative)
      REAL(KIND=dp) , ALLOCATABLE , TARGET :: WalkVecH(:)                    !Diagonal hamiltonian element
      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: SpawnVec(:,:),SpawnVec2(:,:)
      INTEGER , ALLOCATABLE , TARGET :: SpawnSignVec(:),SpawnSignVec2(:)
    
      INTEGER :: WalkVecDetsTag=0,WalkVecSignTag=0
      INTEGER :: WalkVecHTag=0
      INTEGER :: SpawnVecTag=0,SpawnVec2Tag=0,SpawnSignVecTag=0,SpawnSignVec2Tag=0

!Pointers to point at the correct arrays for use
      INTEGER(KIND=n_int) , POINTER :: CurrentDets(:,:)
      INTEGER , POINTER :: CurrentSign(:)
      REAL*8 , POINTER :: CurrentH(:)
      INTEGER(KIND=n_int) , POINTER :: SpawnedParts(:,:),SpawnedParts2(:,:)
      INTEGER , POINTER :: SpawnedSign(:),SpawnedSign2(:)

      INTEGER :: ParentInitiator                                !This is a variable for the CASSTAR approximation - keeps track of where spawned walkers have come from.
      INTEGER :: NoAbortedInCAS,NoAbortedOutCAS,NoInCAS,NoOutCAS,HighPopNeg,HighPopPos,MaxInitPopNeg,MaxInitPopPos
      REAL*8 :: AllGrowRateAbort,NoAborted,AllNoAborted,NoAddedInitiators,AllNoAddedInitiators,NoInitDets,AllNoInitDets
      REAL(KIND=dp) :: AllNoAbortedOld 
      REAL*8 :: NoNonInitDets,NoInitWalk,NoNonInitWalk,NoDoubSpawns,InitRemoved,AllInitRemoved
      REAL*8 :: AllNoNonInitDets,AllNoInitWalk,AllNoNonInitWalk,AllNoDoubSpawns
      REAL*8 :: NoExtraInitDoubs,AllNoExtraInitDoubs
      LOGICAL :: tHFInitiator,tPrintHighPop
 
      REAL*8 :: AvDiagSftAbort,SumDiagSftAbort,DiagSftAbort     !This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying, and
                                                                !the instantaneous shift, including the number of aborted as though they had lived.
    
      INTEGER , ALLOCATABLE :: HFDet(:)       !This will store the HF determinant
      INTEGER :: HFDetTag=0

      INTEGER :: MaxWalkersPart,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps,MaxWalkersAnnihil
      INTEGER :: TotParts,TotPartsOld
      INTEGER :: exFlag=3

!The following variables are calculated as per processor, but at the end of each update cycle, are combined to the root processor
      REAL*8 :: GrowRate,DieRat,ProjectionE,SumENum
      INTEGER(KIND=int64) :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
      REAL*8 :: AvSign           !This is the average sign of the particles on each node
      REAL*8 :: AvSignHFD        !This is the average sign of the particles at HF or Double excitations on each node
      INTEGER(KIND=int64) :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
      INTEGER :: Annihilated      !This is the number annihilated on one processor
      INTEGER :: NoatHF           !This is the instantaneous number of particles at the HF determinant
      INTEGER :: NoatDoubs
      INTEGER :: Acceptances      !This is the number of accepted spawns - this is only calculated per node.
      REAL*8 :: AccRat            !Acceptance ratio for each node over the update cycle
      INTEGER :: PreviousCycles   !This is just for the head node, so that it can store the number of previous cycles when reading from POPSFILE
      INTEGER :: NoBorn,NoDied
      INTEGER :: SpawnFromSing  !These will output the number of particles in the last update cycle which have been spawned by a single excitation.
      INTEGER :: AllSpawnFromSing
      INTEGER :: HFPopCyc         !This is the number of update cycles which have a HF particle at some point
      INTEGER :: HFCyc            !This is the number of HF*sign particles on a given processor over the course of the update cycle
      REAL*8 :: AllHFCyc          !This is the sum of HF*sign particles over all processors over the course of the update cycle
      REAL*8 :: ENumCyc           !This is the sum of doubles*sign*Hij on a given processor over the course of the update cycle
      REAL*8 :: AllENumCyc        !This is the sum of double*sign*Hij over all processors over the course of the update cycle
      REAL*8 :: ProjEIter,ProjEIterSum    !This is the energy estimator where each update cycle contributes an energy and each is given equal weighting.
      integer :: iPartBloom   ! The maximum number of children spawned from a
                              ! single excitation. Used to calculate blooms.

!These are the global variables, calculated on the root processor, from the values above
      REAL*8 :: AllGrowRate
      REAL(KIND=dp) :: AllTotWalkers,AllTotWalkersOld,AllTotParts,AllTotPartsOld
      INTEGER(KIND=int64) :: AllSumWalkersCyc
      INTEGER :: AllAnnihilated,AllNoatHF,AllNoatDoubs
      REAL*8 :: AllSumNoatHF,AllSumENum,AllAvSign,AllAvSignHFD
      INTEGER :: AllNoBorn,AllNoDied,MaxSpawned
  
      HElement_t :: rhii
      REAL*8 :: Hii,Fii

      LOGICAL :: TSinglePartPhase                 !This is true if TStartSinglePart is true, and we are still in the phase where the shift is fixed and particle numbers are growing

!      INTEGER :: mpilongintegertype               !This is used to create an MPI derived type to cope with 8 byte integers

      LOGICAL :: TDebug                           !Debugging flag
      INTEGER :: MaxIndex

      LOGICAL :: tErrorBlocking=.false.           !This becomes true when the blocking error analysis begins, and initiates the calling of the blocking routine.
      LOGICAL :: tShiftBlocking=.false.

      LOGICAL :: TTruncSpace=.false.              !This is a flag set as to whether the excitation space should be truncated or not.
      LOGICAL :: TFlippedSign=.false.             !This is to indicate when the sign of the particles have been flipped. This is needed for the calculation of the ACF

      type(timer) :: Walker_Time, Annihil_Time,ACF_Time, Sort_Time, &
                           Comms_Time, AnnSpawned_time, AnnMain_time, &
                           BinSearch_time

!These are variables needed for the FixCASshift option in which an active space is chosen and the shift fixed only for determinants within this space
!The SpinInvBRR vector stores the energy ordering for each spatial orbital, which is the inverse of the BRR vector
      INTEGER, ALLOCATABLE :: SpinInvBRR(:)
      INTEGER :: SpinInvBRRTag=0
      INTEGER :: CASmin=0,CASmax=0

      ! The approximate fraction of singles and doubles. This is calculated
      ! using the HF determinant, if using non-uniform random excitations.
      real*8 :: pDoubles, pSingles
      
      ! Bit representation of the HF determinant
      integer(kind=n_int), allocatable :: iLutHF(:)
    
      REAL(4) :: IterTime
    
      REAL(KIND=dp) , ALLOCATABLE :: Histogram(:),AllHistogram(:),InstHist(:),AllInstHist(:),AttemptHist(:),AllAttemptHist(:),SpawnHist(:),AllSpawnHist(:)
      REAL(KIND=dp) , ALLOCATABLE :: AvAnnihil(:),AllAvAnnihil(:),InstAnnihil(:),AllInstAnnihil(:)
      REAL(KIND=dp) , ALLOCATABLE :: SinglesAttemptHist(:),AllSinglesAttemptHist(:),SinglesHist(:),AllSinglesHist(:),DoublesHist(:),AllDoublesHist(:),DoublesAttemptHist(:),AllDoublesAttemptHist(:)
      REAL(KIND=dp) , ALLOCATABLE :: SinglesHistOccOcc(:),SinglesHistOccVirt(:),SinglesHistVirtOcc(:),SinglesHistVirtVirt(:)
      REAL(KIND=dp) , ALLOCATABLE :: AllSinglesHistOccOcc(:),AllSinglesHistVirtOcc(:),AllSinglesHistOccVirt(:),AllSinglesHistVirtVirt(:)

      INTEGER :: MaxDet,iOffDiagNoBins
      INTEGER , ALLOCATABLE :: HistMinInd(:),HistMinInd2(:)

      INTEGER , ALLOCATABLE :: DoublesDets(:,:)
      INTEGER :: DoublesDetsTag,NoDoubs

      INTEGER , ALLOCATABLE :: ValidSpawnedList(:) !This is used for the direct annihilation, and ValidSpawnedList(i) indicates the next free slot in the processor iProcIndex ( 0 -> nProcessors-1 )
      INTEGER , ALLOCATABLE :: InitialSpawnedSlots(:) !This is set up as the initial ValidSpawnedList elements, so that it does not need to be reevaluated each time.

      INTEGER :: WalkersDiffProc

      LOGICAL :: tGenMatHEl=.true.      !This is whether to generate matrix elements as generating excitations for the HPHF option

      INTEGER :: VaryShiftCycles                    !This is the number of update cycles that the shift has allowed to vary for.
      INTEGER :: VaryShiftIter                     !This is the iteration that the shift can vary.
      REAL*8 :: AvDiagSft,SumDiagSft                !This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying.

      REAL*8 , ALLOCATABLE :: HistHamil(:,:),AllHistHamil(:,:),AvHistHamil(:,:),AllAvHistHamil(:,:) !These arrays are for histogramming the hamiltonian when tHistHamil is set.
      REAL*8 :: TotImagTime
            
      INTEGER(KIND=n_int) , ALLOCATABLE :: CASMask(:)        !These are masking arrays for the core and external orbitals in the cas space
      INTEGER(KIND=n_int) , ALLOCATABLE :: CoreMask(:)       !These are masking arrays for the Core orbitals in the cas space

      INTEGER , ALLOCATABLE :: RandomHash(:)    !This is a random indexing scheme by which the orbital indices are randomised to attempt to provide a better hashing performance

      INTEGER :: HFIter
      REAL*8 :: ENumIter,IterEnergy

      REAL*8 :: HFShift     !A 'shift'-like value for the total energy which is taken from the growth of walkers on the HF determinant.
      REAL*8 :: InstShift   !An instantaneous value for the shift from the growth of walkers.
      INTEGER :: OldAllNoatHF

      INTEGER :: iHFProc    !Processor index for HF determinant

      !This data is for calculating the highest population determinant, and potentially restarting the calculation based on this determinant, or changing the determiant which the energy is calculated from.
      INTEGER :: iHighestPop
      INTEGER , ALLOCATABLE :: ProjEDet(:)
      INTEGER(KIND=n_int) , ALLOCATABLE :: HighestPopDet(:),iLutRef(:)

      ! These are variables used to control the behaviour of PerformFciMCycPar
      ! without passing them directly to it.
      character(150) :: bloom_warn_string
      integer :: max_calc_ex_level
      
      
      
      !*****************  Redundant variables ************************
    
      TYPE ExcitGenerator
          INTEGER , POINTER :: ExcitData(:)=>null()   !This stores the excitation generator
          INTEGER :: nExcitMemLen                     !This is the length of the excitation generator
          INTEGER :: nPointed                         !This indicates the number of elements in the excitation pointer arrays which are pointing to this position
      END TYPE

      TYPE ExcitPointer
          INTEGER , POINTER :: PointToExcit(:)        !This is a pointer to the excitation generator in ExcitGens
          INTEGER :: IndexinExArr                     !This is the index in ExcitGens which we are pointing at
      END TYPE

      TYPE(ExcitGenerator) , ALLOCATABLE :: ExcitGens(:)  !This is the array to store all the excitation generators
      INTEGER , ALLOCATABLE :: FreeIndArray(:)            !This is a circular list of the free positions to put excitation generators in ExcitGens.
      INTEGER :: BackofList,FrontOfList                   !This indicates where in the list we are.
                                                        !We add indices which become available in ExcitGens to the front of the list, and when we use them up,
                                                        !take off the back of the list. It is circular and so will repeat indefinitly.

      INTEGER , ALLOCATABLE , TARGET :: WalkVec2Dets(:,:),WalkVec2Sign(:)
      REAL(KIND=dp) , ALLOCATABLE , TARGET :: WalkVec2H(:)
      INTEGER , ALLOCATABLE :: IndexTable(:),Index2Table(:)                               !Indexing for the annihilation
      INTEGER , ALLOCATABLE :: ProcessVec(:),Process2Vec(:)                               !Index for process rank of original walker
      INTEGER(KIND=int64) , ALLOCATABLE :: HashArray(:),Hash2Array(:)                         !Hashes for the walkers when annihilating
      INTEGER :: HashArrayTag=0,Hash2ArrayTag=0,IndexTableTag=0,Index2TableTag=0,ProcessVecTag=0,Process2VecTag=0
      INTEGER :: WalkVe2HTag=0,WalkVec2DetsTag=0,WalkVec2SignTag=0
      INTEGER , POINTER :: NewDets(:,:)
      INTEGER , POINTER :: NewSign(:)
      REAL*8 , POINTER :: NewH(:)
      TYPE(ExcitPointer) , POINTER :: CurrentExcits(:), NewExcits(:)

      TYPE(ExcitGenerator) :: HFExcit         !This is the excitation generator for the HF determinant
      INTEGER(KIND=int64) :: HFHash               !This is the hash for the HF determinant
!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle for each variable
      INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
      INTEGER :: CullInfo(10,3)

!      INTEGER , ALLOCATABLE :: AllowedDetList(:,:) !If tListDets is on, this array will fill with allowed determinants to spawn at
!      INTEGER :: NAllowedDetList   !This is the number of allowed determinants to spawn at in the AllowedDetList.


END MODULE FciMCData
