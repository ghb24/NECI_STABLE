MODULE FciMCData
      use, intrinsic :: iso_c_binding
	  use SystemData, only: BasisFN
      use constants, only: dp, int64, n_int, lenof_sign
      use global_utilities

      implicit none
      save

      INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor

      ! Units used to write to files
      integer :: fcimcstats_unit ! FCIMCStats
      integer :: initiatorstats_unit ! INITIATORStats
      integer :: ComplexStats_unit ! COMPLEXStats

      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: WalkVecDets(:,:)                !Contains determinant list
      REAL(KIND=dp) , ALLOCATABLE , TARGET :: WalkVecH(:)                    !Diagonal hamiltonian element
      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: SpawnVec(:,:),SpawnVec2(:,:)
    
      INTEGER :: WalkVecDetsTag=0
      INTEGER :: WalkVecHTag=0
      INTEGER :: SpawnVecTag=0,SpawnVec2Tag=0

!Pointers to point at the correct arrays for use
      INTEGER(KIND=n_int) , POINTER :: CurrentDets(:,:)
      REAL*8 , POINTER :: CurrentH(:)
      INTEGER(KIND=n_int) , POINTER :: SpawnedParts(:,:),SpawnedParts2(:,:)

      INTEGER :: ParentInitiator                                !This is a variable for the CASSTAR approximation - keeps track of where spawned walkers have come from.
      INTEGER :: NoAbortedInCAS,NoAbortedOutCAS,NoInCAS,NoOutCAS,HighPopNeg,HighPopPos,MaxInitPopNeg,MaxInitPopPos

    integer(int64) :: NoAborted, NoAddedInitiators, NoInitDets, NoNonInitDets
    integer(int64) :: NoInitWalk, NoNonInitWalk, NoDoubspawns
    integer(int64) :: NoExtraInitDoubs, InitRemoved

    integer(int64) :: AllNoAborted, AllNoAddedInitiators, AllNoInitDets
    integer(int64) :: AllNoNonInitDets, AllNoInitWalk, AllNoNonInitWalk
    integer(int64) :: AllNodoubSpawns, AllNoExtraInitDoubs, AllInitRemoved
    integer(int64) :: AllNoAbortedOld, AllGrowRateAbort

      LOGICAL :: tHFInitiator,tPrintHighPop
 
      REAL*8 :: AvDiagSftAbort,SumDiagSftAbort,DiagSftAbort     !This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying, and
                                                                !the instantaneous shift, including the number of aborted as though they had lived.

      REAL*8 :: DiagSftRe,DiagSftIm     !For complex walkers - this is just for info - not used for population control.
    
      INTEGER , ALLOCATABLE :: HFDet(:)       !This will store the HF determinant
      INTEGER :: HFDetTag=0

      INTEGER :: MaxWalkersPart,PreviousNMCyc,Iter,NoComps,MaxWalkersAnnihil
      integer(int64) :: TotWalkers, TotWalkersOld
      integer(int64), dimension(lenof_sign) :: TotParts, TotPartsOld
      INTEGER :: exFlag=3

!The following variables are calculated as per processor, but at the end of each update cycle, are combined to the root processor
      REAL*8 :: GrowRate,DieRat,ProjectionE,SumENum
      integer(int64) :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
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
      integer(int64) :: AllTotWalkers, AllTotWalkersOld
      integer(int64), dimension(lenof_sign) :: AllTotParts, AllTotPartsOld
      integer(int64) :: AllSumNoatHF
      INTEGER(KIND=int64) :: AllSumWalkersCyc
      INTEGER :: AllAnnihilated,AllNoatHF,AllNoatDoubs
      REAL*8 :: AllSumENum,AllAvSign,AllAvSignHFD
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

      LOGICAL , PARAMETER :: tGenMatHEl=.true.      !This is whether to generate matrix elements as generating excitations for the HPHF option

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

      ! Store data about all processors for calculating load balancing
      integer(int64) :: MaxWalkersProc, MinWalkersProc

	  TYPE(BasisFN) :: HFSym


      ! ********************** FCIMCPar control variables *****************
      ! Store data from one fcimc iteration
      !  --> We can deal with different types of iteration separately
      type fcimc_iter_data
          integer, dimension(lenof_sign) :: nborn
          integer, dimension(lenof_sign) :: ndied
          integer, dimension(lenof_sign) :: nannihil
          integer, dimension(lenof_sign) :: naborted
          integer, dimension(lenof_sign) :: update_growth, update_growth_tot
          integer(int64), dimension(lenof_sign) :: tot_parts_old
          integer :: update_iters
      end type
      
      ! These are variables used to control the behaviour of PerformFciMCycPar
      ! without passing them directly to it.
      character(150) :: bloom_warn_string
      integer :: max_calc_ex_level
      type(fcimc_iter_data), target :: iter_data_fciqmc
      type(fcimc_iter_data), target :: iter_data_ccmc



      ! Here are the FUNCTION POINTERS for use with PerformFciMCycPar
      ! Use with extreme care, and keep your interfaces up to date or bad
      ! things (namely segfaults) will happen

      type(c_ptr) :: ptr_excit_generator
      type(c_ptr) :: ptr_attempt_create
      type(c_ptr) :: ptr_get_spawn_helement
      type(c_ptr) :: ptr_new_child_stats
      type(c_ptr) :: ptr_encode_child
      type(c_ptr) :: ptr_attempt_die
      type(c_ptr) :: ptr_iter_data

      integer :: yama_global (4)
      
      !*****************  Redundant variables ************************
    

      INTEGER , ALLOCATABLE , TARGET :: WalkVec2Dets(:,:),WalkVec2Sign(:)
      REAL(KIND=dp) , ALLOCATABLE , TARGET :: WalkVec2H(:)
      INTEGER , ALLOCATABLE :: IndexTable(:),Index2Table(:)                               !Indexing for the annihilation
      INTEGER , ALLOCATABLE :: ProcessVec(:),Process2Vec(:)                               !Index for process rank of original walker
      INTEGER(KIND=int64) , ALLOCATABLE :: HashArray(:),Hash2Array(:)                         !Hashes for the walkers when annihilating
      INTEGER :: HashArrayTag=0,Hash2ArrayTag=0,IndexTableTag=0,Index2TableTag=0,ProcessVecTag=0,Process2VecTag=0
      INTEGER :: WalkVe2HTag=0,WalkVec2DetsTag=0,WalkVec2SignTag=0
      INTEGER , POINTER :: NewDets(:,:)
      INTEGER , POINTER :: NewSign(:)

      ! Only used in FciMC, but put here to allow access to a data module for
      ! the sorting routines etc.
      type excitGenerator
          integer, pointer :: excitdata(:) ! The excitation generator
          integer :: nExcitMemLen = 0           ! Length of excitation gen.
          logical :: excitGenForDet = .false.   ! true if excitation generator
                                               ! stored corresponds to the
                                               ! determinant.
      end type

      type(ExcitGenerator) :: HFExcit         !This is the excitation generator for the HF determinant
      integer(int64) :: HFHash               !This is the hash for the HF determinant
!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle for each variable
      INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
      INTEGER :: CullInfo(10,3)

      ! Used for modifying the ReadPops procedures, so that we can call 
      ! InitFCIMCCalcPar again without reading the popsfile.
      logical :: tPopsAlreadyRead


      interface assignment(=)
          module procedure excitgenerator_assign
      end interface

contains
    
    pure subroutine excitgenerator_init (egen)
        type(excitGenerator), intent(inout) :: egen
        nullify(egen%excitdata)
        egen%nExcitMemLen = 0
        egen%excitGenForDet = .false.
    end subroutine

    pure subroutine excitgenerator_destroy (egen)
        type(excitGenerator), intent(inout) :: egen
        if (associated(egen%excitdata)) deallocate(egen%excitdata)
        nullify(egen%excitdata)
        egen%nExcitMemLen = 0
        egen%excitGenForDet = .false.
    end subroutine

    elemental subroutine excitgenerator_assign (lhs, rhs)
        type(excitGenerator), intent(inout) :: lhs
        type(excitGenerator), intent(in) :: rhs

        if (associated(lhs%excitData)) deallocate(lhs%excitData)

        ! Do we actually want the excitation generator? Is it for the correct
        ! determinant?
        if (rhs%excitGenForDet) then
            ! Now copy the excitation generator.
            allocate (lhs%excitData(rhs%nExcitMemLen))

            lhs%excitData = rhs%excitData
        endif
        lhs%nExcitMemLen = rhs%nExcitMemLen
        lhs%excitGenForDet = rhs%excitGenForDet
    end subroutine

END MODULE FciMCData
