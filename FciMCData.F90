MODULE FciMCData
      USE HElem
      USE global_utilities
      IMPLICIT NONE
      SAVE

      INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
      INTEGER , PARAMETER :: r2=kind(0.d0)
      INTEGER , PARAMETER :: i2=SELECTED_INT_KIND(18)
!      REAL*8 , PARAMETER :: HEpsilon=0.D+0
    
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

      TYPE(ExcitPointer) , ALLOCATABLE , TARGET :: WalkVecExcits(:),WalkVec2Excits(:)   !This will store the excitation generators for the particles on each node

      INTEGER , ALLOCATABLE , TARGET :: WalkVecDets(:,:),WalkVec2Dets(:,:)                !Contains determinant list
      INTEGER , ALLOCATABLE , TARGET :: WalkVecSign(:),WalkVec2Sign(:)                    !Contains sign list (1 = positive, -1 = negative)
!      INTEGER , ALLOCATABLE , TARGET :: WalkVecIC(:),WalkVec2IC(:)                        !Contains excit level list
      REAL(KIND=r2) , ALLOCATABLE , TARGET :: WalkVecH(:),WalkVec2H(:)                    !Diagonal hamiltonian element
      INTEGER , ALLOCATABLE :: IndexTable(:),Index2Table(:)                               !Indexing for the annihilation
      INTEGER , ALLOCATABLE :: ProcessVec(:),Process2Vec(:)                               !Index for process rank of original walker
      INTEGER(KIND=i2) , ALLOCATABLE :: HashArray(:),Hash2Array(:)                         !Hashes for the walkers when annihilating
      INTEGER , ALLOCATABLE , TARGET :: SpawnVec(:,:),SpawnVec2(:,:)
      INTEGER , ALLOCATABLE , TARGET :: SpawnSignVec(:),SpawnSignVec2(:)
    
      INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0
      INTEGER :: WalkVecHTag=0,WalkVec2HTag=0
      INTEGER :: HashArrayTag=0,Hash2ArrayTag=0,IndexTableTag=0,Index2TableTag=0,ProcessVecTag=0,Process2VecTag=0
      INTEGER :: SpawnVecTag=0,SpawnVec2Tag=0,SpawnSignVecTag=0,SpawnSignVec2Tag=0

!Pointers to point at the correct arrays for use
      INTEGER , POINTER :: CurrentDets(:,:), NewDets(:,:)
      INTEGER , POINTER :: CurrentSign(:), NewSign(:)
!      INTEGER , POINTER :: CurrentIC(:), NewIC(:)
      REAL*8 , POINTER :: CurrentH(:), NewH(:)
      INTEGER , POINTER :: SpawnedParts(:,:),SpawnedParts2(:,:)
      INTEGER , POINTER :: SpawnedSign(:),SpawnedSign2(:)
      TYPE(ExcitPointer) , POINTER :: CurrentExcits(:), NewExcits(:)

      INTEGER :: ParentInitiator,NoAborted,AllNoAborted                     !This is a variable for the CASSTAR approximation - keeps track of where spawned walkers have come from.
      INTEGER :: NoAbortedInCAS,NoAbortedOutCAS,NoInCAS,NoOutCAS
      REAL*8 :: AllGrowRateAbort
      REAL(KIND=r2) :: AllNoAbortedOld 
      INTEGER :: NoAddedInitiators,NoInitDets,NoNonInitDets,NoInitWalk,NoNonInitWalk,NoDoubSpawns
      INTEGER :: AllNoAddedInitiators,AllNoInitDets,AllNoNonInitDets,AllNoInitWalk,AllNoNonInitWalk,AllNoDoubSpawns
      INTEGER :: NoExtraInitDoubs,AllNoExtraInitDoubs
 
      REAL*8 :: AvDiagSftAbort,SumDiagSftAbort,DiagSftAbort     !This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying, and
                                                                !the instantaneous shift, including the number of aborted as though they had lived.
    
      INTEGER , ALLOCATABLE :: HFDet(:)       !This will store the HF determinant
      INTEGER :: HFDetTag=0
      TYPE(ExcitGenerator) :: HFExcit         !This is the excitation generator for the HF determinant
      INTEGER(KIND=i2) :: HFHash               !This is the hash for the HF determinant

      INTEGER :: MaxWalkersPart,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps,MaxWalkersAnnihil
      INTEGER :: TotParts,TotPartsOld
      INTEGER :: exFlag=3

!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle for each variable
      INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
      INTEGER :: CullInfo(10,3)

!The following variables are calculated as per processor, but at the end of each update cycle, are combined to the root processor
      REAL*8 :: GrowRate,DieRat,ProjectionE,SumENum
      INTEGER*8 :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
      REAL*8 :: AvSign           !This is the average sign of the particles on each node
      REAL*8 :: AvSignHFD        !This is the average sign of the particles at HF or Double excitations on each node
      INTEGER(KIND=i2) :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
!      REAL*8 :: MeanExcitLevel    
!      INTEGER :: MinExcitLevel
!      INTEGER :: MaxExcitLevel
      INTEGER :: Annihilated      !This is the number annihilated on one processor
      INTEGER :: LocalAnn         !This is the number of locally annihilated particles
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
!      REAL*8 :: DetsNorm          !This is the sum of the square of the particles on each determinant. It will be meaningless unless in serial, or using rotoannihil, or better - direct annihil.

!These are the global variables, calculated on the root processor, from the values above
      REAL*8 :: AllGrowRate
!      REAL*8 :: AllMeanExcitLevel
!      INTEGER :: AllMinExcitLevel
      REAL(KIND=r2) :: AllTotWalkers,AllTotWalkersOld,AllTotParts,AllTotPartsOld!,AllDetsNorm
!      INTEGER :: AllMaxExcitLevel
      INTEGER(KIND=i2) :: AllSumWalkersCyc
      INTEGER :: AllAnnihilated,AllNoatHF,AllNoatDoubs,AllLocalAnn
      REAL*8 :: AllSumNoatHF,AllSumENum,AllAvSign,AllAvSignHFD
      INTEGER :: AllNoBorn,AllNoDied,MaxSpawned
  
      REAL*8 :: MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector

      LOGICAL :: tCleanRun    !This will be true if the options are there to use the clean code (this should be faster.)

      TYPE(HElement) :: rhii
      REAL*8 :: Hii,Fii

      REAL*8 , ALLOCATABLE :: GraphRhoMat(:,:)    !This stores the rho matrix for the graphs in resumFCIMC
      INTEGER :: GraphRhoMatTag=0

      REAL*8 , ALLOCATABLE :: GraphVec(:)         !This stores the final components for the propagated graph in ResumFCIMC
      INTEGER :: GraphVecTag=0

      REAL*8 , ALLOCATABLE :: GraphKii(:)         !This stores the diagonal Kii matrix elements for the determinants in the graph
      INTEGER :: GraphKiiTag=0

      INTEGER , ALLOCATABLE :: DetsinGraph(:,:)   !This stores the determinants in the graph created for ResumFCIMC
      INTEGER :: DetsinGraphTag=0

      LOGICAL :: TSinglePartPhase                 !This is true if TStartSinglePart is true, and we are still in the phase where the shift is fixed and particle numbers are growing

!      INTEGER :: mpilongintegertype               !This is used to create an MPI derived type to cope with 8 byte integers

      LOGICAL :: TBalanceNodes                    !This is true when the nodes need to be balanced

      LOGICAL :: TDebug                           !Debugging flag
      INTEGER :: MaxIndex

      LOGICAL :: tErrorBlocking=.false.           !This becomes true when the blocking error analysis begins, and initiates the calling of the blocking routine.
      LOGICAL :: tShiftBlocking=.false.

      LOGICAL :: tNoSpinSymExcitgens=.false.      !Whether to use the new full excitation generators which work with non-spin-symmetric symmetries.

      LOGICAL :: TTruncSpace=.false.              !This is a flag set as to whether the excitation space should be truncated or not.
      LOGICAL :: TFlippedSign=.false.             !This is to indicate when the sign of the particles have been flipped. This is needed for the calculation of the ACF

      TYPE(timer), save :: Walker_Time,Annihil_Time,ACF_Time,Sort_Time,Comms_Time,AnnSpawned_time,AnnMain_time,BinSearch_time,CSF_H_Time,timer_A,timer_B,timer_C,timer_D,timer_E, timer_F

!These are arrays used to store the autocorrelation function
      INTEGER , ALLOCATABLE :: WeightatDets(:)                   !First index - det which is stored, second - weight on proc at that iteration
      INTEGER , ALLOCATABLE :: AutoCorrDets(:,:)                  !(NEl,NoAutoDets)
      INTEGER :: WeightatDetsTag=0,AutoCorrDetsTag=0
      INTEGER :: NoAutoDets

!These are variables relating to calculating the population density of an excitation level for use with TLocalAnnihilation
!Local annihilation is currently commented out
      REAL*8 , ALLOCATABLE :: ApproxExcitDets(:)    !This is the approximate size of each excitation level
      INTEGER , ALLOCATABLE :: PartsInExcitLevel(:) !This is the number of walkers for a given iteration in that excitation level
    
!Variables for magnetisation
      INTEGER , ALLOCATABLE :: MagDets(:,:)       !This is to hold the NoMagDets-1 magnetic determinant paths (HF is always magnetic)
      INTEGER , ALLOCATABLE :: MagDetsSign(:)      !This is to hold the sign of the NoMagDets-1 magnetic determinants (HF is always positive)
      INTEGER :: MagDetsTag=0,MagDetsSignTag=0

!These are variables needed for the FixCASshift option in which an active space is chosen and the shift fixed only for determinants within this space
!The SpinInvBRR vector stores the energy ordering for each spatial orbital, which is the inverse of the BRR vector
      INTEGER, ALLOCATABLE :: SpinInvBRR(:)
      INTEGER :: SpinInvBRRTag=0
      INTEGER :: CASmin=0,CASmax=0

      REAL*8 :: pDoubles                          !This is the approximate fraction of excitations which are doubles. This is calculated if we are using non-uniform
                                                !random excitations.
      INTEGER , ALLOCATABLE :: iLutHF(:)          !This is the bit representation of the HF determinant.
    
      REAL*8 , ALLOCATABLE :: OneRDM(:,:)         !This is the 1 electron reduced density matrix.
      INTEGER :: OneRDMTag=0                      !It is calculated as an FCIMC run progresses.  As the run tends towards the correct wavefunction, diagonalisation 
      INTEGER :: Orbs(2)                          !of the 1RDM gives linear combinations of the HF orbitals which tend towards the natural orbitals of the system.

      REAL(4) :: IterTime
    
      REAL(KIND=r2) , ALLOCATABLE :: Histogram(:),AllHistogram(:),InstHist(:),AllInstHist(:),AttemptHist(:),AllAttemptHist(:),SpawnHist(:),AllSpawnHist(:)
      REAL(KIND=r2) , ALLOCATABLE :: AvAnnihil(:),AllAvAnnihil(:),InstAnnihil(:),AllInstAnnihil(:)
      REAL(KIND=r2) , ALLOCATABLE :: SinglesAttemptHist(:),AllSinglesAttemptHist(:),SinglesHist(:),AllSinglesHist(:),DoublesHist(:),AllDoublesHist(:),DoublesAttemptHist(:),AllDoublesAttemptHist(:)
      REAL(KIND=r2) , ALLOCATABLE :: SinglesHistOccOcc(:),SinglesHistOccVirt(:),SinglesHistVirtOcc(:),SinglesHistVirtVirt(:)
      REAL(KIND=r2) , ALLOCATABLE :: AllSinglesHistOccOcc(:),AllSinglesHistVirtOcc(:),AllSinglesHistOccVirt(:),AllSinglesHistVirtVirt(:)

      INTEGER :: MaxDet,iOffDiagNoBins
      INTEGER , ALLOCATABLE :: HistMinInd(:),HistMinInd2(:)

      INTEGER :: GuideFuncDetsTag,GuideFuncSignTag,DetstoRotateTag,SigntoRotateTag,DetstoRotate2Tag,SigntoRotate2Tag,AlliInitGuideParts,InitGuideFuncSignTag
      INTEGER :: GuideFuncHFIndex,GuideFuncHF
      REAL*8 :: GuideFuncDoub
      INTEGER , ALLOCATABLE :: GuideFuncDets(:,:),GuideFuncSign(:),DetstoRotate(:,:),SigntoRotate(:),DetstoRotate2(:,:),SigntoRotate2(:),InitGuideFuncSign(:)

      INTEGER , ALLOCATABLE :: DomExcIndex(:),DomDets(:,:)
      INTEGER :: DomExcIndexTag,DomDetsTag,iMinDomLev,iMaxDomLev,iNoDomDets
      INTEGER , ALLOCATABLE :: MinorStarDets(:,:),MinorSpawnDets(:,:),MinorStarParent(:,:),MinorSpawnParent(:,:),MinorStarSign(:),MinorSpawnSign(:)
      INTEGER , ALLOCATABLE :: MinorSpawnDets2(:,:),MinorSpawnSign2(:),MinorSpawnParent2(:,:)
      INTEGER :: MinorStarDetsTag,MinorSpawnDetsTag,MinorStarParentTag,MinorSpawnParentTag,MinorStarSignTag,MinorSpawnSignTag,MinorStarHiiTag,MinorStarHijTag,NoMinorWalkers
      INTEGER :: MinorSpawnDets2Tag,MinorSpawnSign2Tag,MinorSpawnParent2Tag,MinorAnnihilated,AllMinorAnnihilated
      TYPE(HElement), ALLOCATABLE :: MinorStarHii(:),MinorStarHij(:)
      REAL*8 :: AllNoMinorWalkers

      INTEGER , ALLOCATABLE :: DoublesDets(:,:)
      INTEGER :: DoublesDetsTag,NoDoubs

      INTEGER , ALLOCATABLE :: ValidSpawnedList(:) !This is used for the direct annihilation, and ValidSpawnedList(i) indicates the next free slot in the processor iProcIndex ( 0 -> nProcessors-1 )
      INTEGER , ALLOCATABLE :: InitialSpawnedSlots(:) !This is set up as the initial ValidSpawnedList elements, so that it does not need to be reevaluated each time.

      INTEGER :: WalkersDiffProc

      INTEGER , ALLOCATABLE :: AllowedDetList(:,:) !If tListDets is on, this array will fill with allowed determinants to spawn at
      INTEGER :: NAllowedDetList   !This is the number of allowed determinants to spawn at in the AllowedDetList.

      LOGICAL :: tGenMatHEl=.true.      !This is whether to generate matrix elements as generating excitations for the HPHF option

      INTEGER :: VaryShiftCycles                    !This is the number of update cycles that the shift has allowed to vary for.
      INTEGER :: VaryShiftIter                     !This is the iteration that the shift can vary.
      REAL*8 :: AvDiagSft,SumDiagSft                !This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying.

      REAL*8 , ALLOCATABLE :: HistHamil(:,:),AllHistHamil(:,:),AvHistHamil(:,:),AllAvHistHamil(:,:) !These arrays are for histogramming the hamiltonian when tHistHamil is set.
      REAL*8 :: TotImagTime


END MODULE FciMCData
