#include "macros.h"

MODULE FciMCData
      use iso_c_hack
      use SystemData, only: BasisFN
      use constants, only: dp, int64, n_int, lenof_sign, sp, inum_runs
      use SymExcitDataMod, only: excit_gen_store_type
      use MemoryManager, only: TagIntType
      use global_utilities         
      use Parallel_neci, only: MPIArg
      use ras_data

      implicit none
      save

      ! Type for creating linked lists for the linear scaling algorithm.
      type ll_node
          integer :: ind
          type(ll_node), pointer :: next => null()
      end type

      integer :: iPopsTimers    !Number of timed popsfiles written out (initiatlised to 1)

      real(dp) :: MaxTimeExit   !Max time before exiting out of MC
      logical :: tTimeExit      !Whether to exit out of MC after an amount of runtime

      ! Units used to write to files
      integer :: fcimcstats_unit ! FCIMCStats
      integer :: fcimcstats_unit2 ! FCIMCStats
      integer :: initiatorstats_unit ! INITIATORStats
      integer :: ComplexStats_unit ! COMPLEXStats
      integer :: Tot_Unique_Dets_Unit 

      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: WalkVecDets(:,:)                !Contains determinant list
      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: SpawnVec(:,:),SpawnVec2(:,:)
      INTEGER(KIND=n_int) , ALLOCATABLE , TARGET :: SpawnVecKP(:,:), SpawnVecKP2(:,:)

      INTEGER(TagIntType) :: WalkVecDetsTag=0
      INTEGER(TagIntType) :: SpawnVecTag=0,SpawnVec2Tag=0

!Pointers to point at the correct arrays for use
      INTEGER(KIND=n_int) , POINTER :: CurrentDets(:,:)
      INTEGER(KIND=n_int) , POINTER :: SpawnedParts(:,:),SpawnedParts2(:,:)
      INTEGER(KIND=n_int) , POINTER :: SpawnedPartsKP(:,:), SpawnedPartsKP2(:,:)

      ! In some instances (such as when applying a perturbation operator) it is
      ! useful to store the vector read in from the popsfile in a separate
      ! array. This is what popsfile_dets is used for.
      integer(n_int), allocatable :: popsfile_dets(:,:)
      ! If true then the above array will be allocated to be the same size as
      ! WalkVecDets.
      logical :: alloc_popsfile_dets

      INTEGER(KIND=n_int) , ALLOCATABLE :: Spawned_Parents(:,:)
      INTEGER , ALLOCATABLE :: Spawned_Parents_Index(:,:)
      INTEGER :: Spawned_ParentsTag, Spawned_Parents_IndexTag
      REAL(dp) :: SumSigns, SumSpawns
      REAL(dp), dimension(lenof_sign) :: AvNoatHF
      LOGICAL :: tFillingStochRDMonFly, tFillingExplicRDMonFly
      logical :: tFill_RDM
      integer :: IterLastRDMFill
      integer :: Spawned_Parts_Zero, HFInd
      integer :: IterRDMStart
      integer, dimension(inum_runs) :: IterRDM_HF
      real(dp), dimension(lenof_sign) :: InstNoatHf
      logical :: tGhostChild, tFinalRDMEnergy

      INTEGER(KIND=n_int) , ALLOCATABLE :: TempSpawnedParts(:,:)
      INTEGER :: TempSpawnedPartsTag, TempSpawnedPartsInd

      ! Be able to store a list of the current initiators
      integer(n_int), allocatable :: CurrentInits(:,:)
      integer :: max_inits
      integer(TagIntType) :: CurrentInitTag=0

      INTEGER :: NoAbortedInCAS,NoAbortedOutCAS,NoInCAS,NoOutCAS, HighPopNeg, HighPopPos
      REAL(dp) :: MaxInitPopNeg,MaxInitPopPos

    real(dp), dimension(lenof_sign) :: NoAborted, AllNoAborted, AllNoAbortedOld
    real(dp), dimension(lenof_sign) :: NoRemoved, AllNoRemoved, AllNoRemovedOld
    integer(int64), dimension(lenof_sign) :: NoAddedInitiators, NoInitDets, NoNonInitDets
    real(dp), dimension(lenof_sign) :: NoInitWalk, NoNonInitWalk
    integer(int64), dimension(lenof_sign) :: NoExtraInitDoubs, InitRemoved
    integer :: no_spatial_init_dets

    integer(int64), dimension(lenof_sign) :: AllNoAddedInitiators, AllNoInitDets
    integer(int64), dimension(lenof_sign) :: AllNoNonInitDets
    real(dp), dimension(inum_runs) :: AllNoInitWalk, AllNoNonInitWalk
    integer(int64), dimension(lenof_sign) :: AllNoExtraInitDoubs, AllInitRemoved
    integer(int64), dimension(lenof_sign) :: AllGrowRateAbort

      LOGICAL :: tHFInitiator,tPrintHighPop, tcurr_initiator
      logical :: tHashWalkerList    !Option to store occupied determinant in a hash table
      integer, allocatable :: FreeSlot(:)   !List of the free slots in the main list
      integer :: iStartFreeSlot     !=1 at the beginning of an iteration, will increment
      !as free slots are used up for *newly spawned* walkers onto previously unoccupied determinants only
      integer :: iEndFreeSlot   !Position of last free slot, so after we exceed this, just add the the end of the main list.
 
!      real(dp) :: AvDiagSftAbort,SumDiagSftAbort,DiagSftAbort     
!This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying, and
                                     !the instantaneous shift, including the number of aborted as though they had lived.

      real(dp) :: DiagSftRe,DiagSftIm     !For complex walkers - this is just for info - not used for population control.
    
      INTEGER , ALLOCATABLE :: HFDet(:), HFDet_True(:)       !This will store the HF determinant
      INTEGER(TagIntType) :: HFDetTag=0

      INTEGER :: MaxWalkersPart,PreviousNMCyc,Iter,NoComps,MaxWalkersAnnihil
      integer :: MaxWalkersUncorrected
      integer(int64) :: TotWalkers, TotWalkersOld
      real(dp), dimension(lenof_sign) :: TotParts, TotPartsOld
      real(dp), dimension(inum_runs) :: norm_psi_squared
      real(dp), dimension(inum_runs) :: norm_semistoch_squared
      real(dp), dimension(inum_runs) :: all_norm_psi_squared
#ifdef __CMPLX
      real(dp) :: norm_psi
      ! The norm of the wavefunction in just the semi-stochastic space.
      real(dp) :: norm_semistoch
#else
      real(dp), dimension(inum_runs) :: norm_psi
      ! The norm of the wavefunction in just the semi-stochastic space.
      real(dp),dimension(inum_runs) :: norm_semistoch
#endif

      INTEGER :: exFlag=3
      real(dp) :: AccumRDMNorm, AccumRDMNorm_Inst, AllAccumRDMNorm
      
      !Hash tables to point to the correct determinants in CurrentDets
      type(ll_node), pointer :: HashIndex(:) 
      integer :: nWalkerHashes    ! The length of hash table.
      real(dp) :: HashLengthFrac

!The following variables are calculated as per processor, but at the end of each update cycle, 
!are combined to the root processor
      real(dp) :: GrowRate,DieRat
      HElement_t, dimension(inum_runs) :: SumENum

      ! The averaged projected energy - calculated from accumulated values.
      HElement_t, dimension(inum_runs) :: ProjectionE

      ! The averaged projected energy - calculated over the last update cycle
      HElement_t, dimension(inum_runs) :: proje_iter

      ! The averaged 'absolute' projected energy - calculated over the last update cycle
      ! The magnitude of each contribution is taken before it is summed in
      HElement_t, dimension(inum_runs) :: AbsProjE

      real(dp), dimension(inum_runs) :: trial_numerator, tot_trial_numerator
      real(dp), dimension(inum_runs) :: trial_denom, tot_trial_denom

      ! This is the sum over all previous cycles of the number of particles at the HF determinant
      real(dp), dimension(lenof_sign) :: SumNoatHF 
      real(dp) :: AvSign           !This is the average sign of the particles on each node
      real(dp) :: AvSignHFD        !This is the average sign of the particles at HF or Double excitations on each node
      real(dp), dimension(inum_runs) :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
      Real(dp), dimension(lenof_sign) :: Annihilated      !This is the number annihilated on one processor
      REAL(dp), DIMENSION(lenof_sign) :: NoatHF           !This is the instantaneous number of particles at the HF determinant
      REAL(dp), dimension(inum_runs) :: NoatDoubs
      REAL(dp), dimension(inum_runs) :: Acceptances      !This is the number of accepted spawns - this is only calculated per node.
      real(dp), dimension(inum_runs) :: AccRat            !Acceptance ratio for each node over the update cycle
      INTEGER :: PreviousCycles
      !This is just for the head node, so that it can store the number of previous cycles when reading from POPSFILE
      REAL(dp),dimension(inum_runs) :: NoBorn,NoDied
      REAL(dp),dimension(inum_runs) :: SpawnFromSing  
      !These will output the number of particles in the last update cycle which have been spawned by a single excitation.
      REAL(dp), dimension(inum_runs) :: AllSpawnFromSing
      REAL(dp), DIMENSION(lenof_sign) :: HFCyc           
      !This is the number of HF*sign particles on a given processor over the course of the update cycle
      HElement_t  :: AllHFCyc(inum_runs) 
      !This is the sum of HF*sign particles over all processors over the course of the update cycle
      HElement_t :: OldAllHFCyc(inum_runs) 
      !This is the old *average* (not sum) of HF*sign over all procs over previous update cycle
      HElement_t :: ENumCyc(inum_runs)
      !This is the sum of doubles*sign*Hij on a given processor over the course of the update c
      HElement_t :: AllENumCyc(inum_runs)
      !This is the sum of double*sign*Hij over all processors over the course of the update cyc
      HElement_t :: ENumCycAbs(inum_runs)
      !This is the sum of abs(doubles*sign*Hij) on a given processor "" "" "" 
      HElement_t :: AllENumCycAbs(inum_runs)
      !This is the sum of abs(double*sign*Hij) over all processors over the course of the updat

      ! The projected energy over the current update cycle.
      HElement_t :: ProjECyc(inum_runs)
      
      real(dp) :: bloom_sizes(0:2), bloom_max(0:2)
      integer :: bloom_count(0:2), all_bloom_count(0:2)

!These are the global variables, calculated on the root processor, from the values above
      real(dp), dimension(inum_runs) :: AllGrowRate
      integer(int64) :: AllTotWalkers, AllTotWalkersOld
      real(dp), dimension(lenof_sign) :: AllTotParts, AllTotPartsOld
      real(dp), dimension(lenof_sign) :: AllSumNoatHF
      real(dp), dimension(inum_runs) :: AllSumWalkersCyc
      real(dp), dimension(inum_runs) :: OldAllAvWalkersCyc    !This is the average number of walkers each iteration over the previous update cycle
      REAL(dp), dimension(inum_runs) :: AllAnnihilated
      REAL(dp), dimension(inum_runs) :: AllNoAtDoubs
      REAl(dp), DIMENSION(lenof_sign) :: AllNoatHF
      HElement_t, dimension(inum_runs) :: sum_proje_denominator, &
                        cyc_proje_denominator, all_cyc_proje_denominator, &
                        all_sum_proje_denominator
      real(dp) :: AllAvSign,AllAvSignHFD
      INTEGER :: MaxSpawned
      REAL(dp), dimension(inum_runs) :: AllNoBorn,AllNoDied

      HElement_t, dimension(inum_runs) :: AllSumENum
  
      HElement_t :: rhii
      real(dp) :: Hii,Fii
!This is true if TStartSinglePart is true, and we are still in the phase where the 
!shift is fixed and particle numbers are growing
      LOGICAL, dimension(inum_runs) :: TSinglePartPhase                 

!      INTEGER :: mpilongintegertype               !This is used to create an MPI derived type to cope with 8 byte integers

      LOGICAL :: TDebug                           !Debugging flag
      INTEGER :: MaxIndex

      integer, dimension(inum_runs) :: iBlockingIter                    !The iteration to begin the automatic blocking from 
 !This becomes true when the blocking error analysis begins, and initiates the calling of the blocking routine.
      LOGICAL :: tErrorBlocking=.false.           
      LOGICAL :: tShiftBlocking=.false.

      LOGICAL :: TTruncSpace=.false.      !This is a flag set as to whether the excitation space should be truncated or not.
!This is to indicate when the sign of the particles have been flipped. This is needed for the calculation of the ACF
      LOGICAL :: TFlippedSign=.false.             

      type(timer) :: Walker_Time, Annihil_Time,ACF_Time, Sort_Time, &
                           Comms_Time, AnnSpawned_time, AnnMain_time, &
                           BinSearch_time, SemiStoch_Comms_Time, &
                           SemiStoch_Multiply_Time, Trial_Search_Time, &
                           SemiStoch_Init_Time, Trial_Init_Time, &
                           kp_generate_time
      
      ! Store the current value of S^2 between update cycles
      real(dp), dimension(inum_runs) :: curr_S2, curr_S2_init

      integer :: HolesInList    !This is for tHashWalkerList and indicates the number of holes in the main list this iter

!These are variables needed for the FixCASshift option in which an active space is chosen and the 
!shift fixed only for determinants within this space
!The SpinInvBRR vector stores the energy ordering for each spatial orbital, which is the inverse of the BRR vector
      INTEGER, ALLOCATABLE :: SpinInvBRR(:)
      INTEGER(TagIntType) :: SpinInvBRRTag=0
      INTEGER :: CASmin=0,CASmax=0

      ! The approximate fraction of singles and doubles. This is calculated
      ! using the HF determinant, if using non-uniform random excitations.
      real(dp) :: pDoubles, pSingles, pParallel
      integer :: nSingles, nDoubles
      
      ! Bit representation of the HF determinant
      integer(kind=n_int), allocatable :: iLutHF(:), iLutHF_True(:)
    
      REAL(KIND=sp) :: IterTime
    
      REAL(KIND=dp) , ALLOCATABLE :: AttemptHist(:),AllAttemptHist(:),SpawnHist(:),AllSpawnHist(:)
      REAL(KIND=dp) , ALLOCATABLE :: AvAnnihil(:,:),AllAvAnnihil(:,:),InstAnnihil(:,:),AllInstAnnihil(:,:)
      REAL(KIND=dp) , ALLOCATABLE :: SinglesAttemptHist(:),AllSinglesAttemptHist(:),SinglesHist(:),AllSinglesHist(:)
      REAL(KIND=dp) , ALLOCATABLE :: DoublesHist(:),AllDoublesHist(:),DoublesAttemptHist(:),AllDoublesAttemptHist(:)
      REAL(KIND=dp) , ALLOCATABLE :: SinglesHistOccOcc(:),SinglesHistOccVirt(:),SinglesHistVirtOcc(:),SinglesHistVirtVirt(:)
      REAL(KIND=dp) , ALLOCATABLE :: AllSinglesHistOccOcc(:),AllSinglesHistVirtOcc(:),AllSinglesHistOccVirt(:)
      REAL(KIND=dp) , ALLOCATABLE :: AllSinglesHistVirtVirt(:)

      real(dp), allocatable :: spin_det_hist(:,:)

      INTEGER :: MaxDet,iOffDiagNoBins

      INTEGER , ALLOCATABLE :: DoublesDets(:,:)
      INTEGER(TagIntType) :: DoublesDetsTag
      INTEGER :: NoDoubs
!This is used for the direct annihilation, and ValidSpawnedList(i) indicates the next 
!free slot in the processor iProcIndex ( 0 -> nProcessors-1 )
      INTEGER , ALLOCATABLE :: ValidSpawnedList(:) 
 !This is set up as the initial ValidSpawnedList elements, so that it does not need to be reevaluated each time.
      INTEGER , ALLOCATABLE :: InitialSpawnedSlots(:) 

      ! For use in kp-fciqmc code. How many different determinants are spawned onto this processor?
      integer :: max_spawned_ind

      integer :: WalkersDiffProc, PartsDiffProc

      !This is whether to generate matrix elements as generating excitations for the HPHF/MI/ISK options
      LOGICAL , PARAMETER :: tGenMatHEl=.true.      

      INTEGER, dimension(inum_runs) :: VaryShiftCycles  !number of update cycles that the shift has allowed to vary for.
      INTEGER, dimension(inum_runs) :: VaryShiftIter                     !This is the iteration that the shift can vary.
      !This is the average diagonal shift value since it started varying, and the sum of the shifts since it started varying.
      real(dp), dimension(inum_runs) :: AvDiagSft,SumDiagSft                

!These arrays are for histogramming the hamiltonian when tHistHamil is set.
      real(dp) , ALLOCATABLE :: HistHamil(:,:),AllHistHamil(:,:),AvHistHamil(:,:),AllAvHistHamil(:,:) 
      real(dp) :: TotImagTime
            
      INTEGER(KIND=n_int) , ALLOCATABLE :: CASMask(:)        !These are masking arrays for the core 
                                                             !and external orbitals in the cas space
      INTEGER(KIND=n_int) , ALLOCATABLE :: CoreMask(:)       !These are masking arrays for the Core orbitals in the cas space

      INTEGER , ALLOCATABLE :: RandomHash(:)    !This is a random indexing scheme by which the orbital indices 
                                                !are randomised to attempt to provide a better hashing performance
      integer, allocatable :: RandomHash2(:)    !Another random index scheme for the hashing used by tHashWalkerList

      real(dp), dimension(inum_runs) :: HFShift     !A 'shift'-like value for the total energy which is taken
                              !from the growth of walkers on the HF determinant.
      real(dp), dimension(inum_runs) :: InstShift   !An instantaneous value for the shift from the growth of walkers.
      REAL(dp), DIMENSION(lenof_sign) :: OldAllNoatHF

      INTEGER :: iHFProc    !Processor index for HF determinant

      !This data is for calculating the highest population determinant, 
      !and potentially restarting the calculation based on this determinant, 
      !or changing the determiant which the energy is calculated from.
      INTEGER :: iHighestPop
      INTEGER :: QuadDetsEst !Estimate of the number of symmetry allowed determinants at excit level 4
      INTEGER :: DoubDetsEst !Estimate of the number of symmetry allowed determinants at excit level 2
      INTEGER , ALLOCATABLE :: ProjEDet(:)
      INTEGER(KIND=n_int) , ALLOCATABLE :: HighestPopDet(:),iLutRef(:)
      INTEGER(n_int) , ALLOCATABLE :: iLutRefFlip(:)     !If we are using HPHF and projecting onto 
                                                        !an open-shell determinant, then it is useful
                                                        !to store the spin-coupled determinant, 
                                                        !so we can calculate projection onto both.
      INTEGER , ALLOCATABLE :: RefDetFlip(:)
      LOGICAL :: tSpinCoupProjE
      
      !Extra data recorded for using RealCoefficients
      INTEGER :: WalkersToSpawn
      LOGICAL :: blank_det

      ! Store data about all processors for calculating load balancing
      integer(int64) :: MaxWalkersProc, MinWalkersProc
      real(dp) :: MaxPartsProc, MinPartsProc

      TYPE(BasisFN) :: HFSym
      integer :: iMaxBloom !If tMaxBloom is on, this stores the largest bloom to date.

      real(dp) :: proje_denominator_cyc(lenof_sign)
      real(dp) :: proje_denominator_sum(lenof_sign)
      logical :: tRestart   !Whether to restart a calculation
      real(dp), dimension(inum_runs) :: InputDiagSft  !Diag shift from the input file if needed to be reset after a restart
      

      ! ********************** FCIMCPar control variables *****************
      ! Store data from one fcimc iteration
      !  --> We can deal with different types of iteration separately
      type fcimc_iter_data
          real(dp), dimension(lenof_sign) :: nborn
          real(dp), dimension(lenof_sign) :: ndied
          real(dp), dimension(lenof_sign) :: nannihil
          real(dp), dimension(lenof_sign) :: naborted
          real(dp), dimension(lenof_sign) :: nremoved
          real(dp), dimension(lenof_sign) :: update_growth, update_growth_tot
          real(dp), dimension(lenof_sign) :: tot_parts_old
          integer :: update_iters
      end type
      
      ! These are variables used to control the behaviour of PerformFciMCycPar
      ! without passing them directly to it.
      character(150) :: bloom_warn_string
      integer :: max_calc_ex_level
      type(fcimc_iter_data), target :: iter_data_fciqmc
      type(fcimc_iter_data), target :: iter_data_ccmc

      integer :: yama_global (4)

      ! Used for modifying the ReadPops procedures, so that we can call 
      ! InitFCIMCCalcPar again without reading the popsfile.
      logical :: tPopsAlreadyRead

!      ! Excitation generation storage 
      type(excit_gen_store_type) :: fcimc_excit_gen_store

      !Tau searching variables
      logical :: tSearchTau
      real(dp) :: MaxTau

      !Variables for diagonalisation of the walker subspace
      integer :: unitWalkerDiag

      !*****************  Yucky globals for AJWT iter-dependent hashes ***********
      integer :: hash_iter       ! An iteration number added to make iteration-dependent hashes
! -Ln_2 (Cycletime), where CycleTime is the average number of cycles until a det returns to its processor
      integer :: hash_shift      

      ! This array stores the Hamiltonian matrix, or part of it, when performing a diagonalisation. It is currently
      ! only used for the code for the Davidson method and semi-stochastic method.
      real(dp), allocatable, dimension(:,:) :: hamiltonian

      integer(TagIntType) :: HamTag, DavidsonTag

      ! Semi-stochastic data.

      ! The diagonal elements of the core-space Hamiltonian (with Hii taken away).
      real(dp), allocatable, dimension(:) :: core_ham_diag
            
      ! This stores the entire core space from all processes, on each process.
      integer(n_int), allocatable, dimension(:,:) :: core_space

      ! This stores all the amplitudes of the walkers in the deterministic space. This vector has the size of the part
      ! of the deterministic space stored on *this* processor only. It is therefore used to store the deterministic vector
      ! on this processor, before it is combined to give the whole vector, which is stored in full_determ_vector.
      ! Later in the iteration, it is also used to store the result of the multiplication by the core Hamiltonian on
      ! full_determ_vector.
      real(dp), allocatable, dimension(:,:) :: partial_determ_vector
      real(dp), allocatable, dimension(:,:) :: full_determ_vector
      real(dp), allocatable, dimension(:,:) :: full_determ_vector_av

      real(dp), allocatable, dimension(:,:) :: partial_determ_vecs_kp
      real(dp), allocatable, dimension(:,:) :: full_determ_vecs_kp

      ! determ_proc_sizes(i) holds the core space size on processor i.
      integer(MPIArg), allocatable, dimension(:) :: determ_proc_sizes
      ! determ_proc_indices(i) holds sum(determ_proc_sizes(i-1)), that is, the
      ! total number of core states on all processors up to processor i.
      ! (determ_proc_indices(1) == 0).
      integer(MPIArg), allocatable, dimension(:) :: determ_proc_indices
      ! The total size of the core space on all processors.
      integer(MPIArg) :: determ_space_size
      ! determ_space_size_int is identical to determ_space_size, but converted
      ! to the default integer kind.
      integer :: determ_space_size_int

      ! This vector will store the indicies of the deterministic states in CurrentDets. This is worked out in the main loop.
      integer, allocatable, dimension(:) :: indices_of_determ_states

      ! If true (as is the case by default) then semi-stochastic calculations
      ! will start from the ground state of the core space
      logical :: tStartCoreGroundState

      ! Trial wavefunction data.

      ! This list stores the iluts from which the trial wavefunction is formed, but only those that reside on this processor.
      ! Not that this is deallocated after initialisation when using the tTrialHash option (turned on by default).
      integer(n_int), allocatable, dimension(:,:) :: trial_space
      ! The number of states in the trial vector space.
      integer :: trial_space_size = 0
      ! This list stores the iluts from which the trial wavefunction is formed, but only those that reside on this processor.
      ! Not that this is deallocated after initialisation when using the tTrialHash option (turned on by default).
      integer(n_int), allocatable, dimension(:,:) :: con_space
      ! The number of states in the space connected to (but not including) the trial vector space.
      integer :: con_space_size = 0

      ! This vector stores the trial wavefunction itself.
      ! Not that this is deallocated after initialisation when using the tTrialHash option (turned on by default).
      real(dp), allocatable, dimension(:) :: trial_wf
      ! This vector stores the values of trial_wf for the occupied trial state in CurrentDets, in the same order as these
      ! states in CurrentDets. If not all trial states are occupied then the final elements store junk and aren't used.
      real(dp), allocatable, dimension(:) :: occ_trial_amps
      ! Holds the number of occupied trial states in CurrentDets.
      integer :: ntrial_occ
      ! Holds the index of the element in occ_trial_amps to access when the next trial space state is found.
      integer :: trial_ind
      ! When new states are about to inserted into CurrentDets, a search is performed to see if any of them are in the
      ! trial space. If they are then the corresponding amplitudes are stored in this vector until they are merged
      ! into occ_trial_amps, to prevent overwriting during the merging process.
      real(dp), allocatable, dimension(:) :: trial_temp

      real(dp) :: trial_energy
      ! This vector's elements store the quantity
      ! \sum_j H_{ij} \psi^T_j,
      ! where \psi is the trial wavefunction. These elements are stored only in the space of states which are connected
      ! to *but not included in* the trial vector space.
      ! Not that this is deallocated after initialisation when using the tTrialHash option (turned on by default).
      real(dp), allocatable, dimension(:) :: con_space_vector
      ! This vector stores the values of con_space_vector for the occupied connected state in CurrentDets, in the same order as
      ! these states in CurrentDets. If not all connected states are occupied then the final elements store junk and aren't used.
      real(dp), allocatable, dimension(:) :: occ_con_amps
      ! Holds the number of occupied connected states in CurrentDets.
      integer :: ncon_occ
      ! Holds the index of the element in occ_con_amps to access when the next connected space state is found.
      integer :: con_ind
      ! When new states are about to inserted into CurrentDets, a search is performed to see if any of them are in the
      ! connected space. If they are then the corresponding amplitudes are stored in this vector until they are merged
      ! into occ_con_amps, to prevent overwriting during the merging process.
      real(dp), allocatable, dimension(:) :: con_temp

      ! If index i in CurrentDets is a trial state then index i of this array stores the corresponding amplitude of
      ! trial_wf. If index i in the CurrentDets is a connected state then index i of this array stores the corresponding
      ! amplitude of con_space_vector. Else, it will be zero.
      real(dp), allocatable, dimension(:) :: current_trial_amps
      ! Only used with the linscalefcimcalgo option: Because in AnnihilateSpawendParts trial and connected states are
      ! sorted in the same order, a smaller section of the trial and connected space can be searched for each state.
      ! These indices hold the indices to be searched from next time.
      integer :: min_trial_ind, min_conn_ind

      ! If true (which it is by default) then the trial and connected space states are stored in a trial_ht and
      ! con_ht and are accessed by a hash lookup.
      logical :: tTrialHash
      ! If true, include the walkers that get cancelled by the initiator criterion in the trial energy estimate.
      ! tTrialHash needs to be used for this option.
      logical :: tIncCancelledInitEnergy

      ! Semi-stochastic tags:
      integer(TagIntType) :: CoreTag, FDetermTag, FDetermAvTag, PDetermTag, IDetermTag, CoreSpaceTag
      logical :: tInDetermSpace

      ! Trial wavefunction tags:
      integer(TagIntType) :: TrialTag, ConTag, ConVecTag, TrialWFTag, TempTag, CurrentTrialTag
      integer(TagIntType) :: TrialTempTag, ConTempTag, OccTrialTag, OccConTag

      ! Data for performing the direct-ci Davidson algorithm.
      type(ras_parameters) :: davidson_ras
      type(ras_class_data), allocatable, dimension(:) :: davidson_classes
      integer, allocatable, dimension(:,:) :: davidson_strings
      integer(n_int), allocatable, dimension(:,:) :: davidson_iluts
      type(direct_ci_excit), allocatable, dimension(:) :: davidson_excits

      real(dp) :: max_cyc_spawn, all_max_cyc_spawn

      ! Type containing information on a perturbation operator, constructed from
      ! a string of creation and annihilation operators.
      type perturbation
          ! The number of annihilation operators in the perturbation.
          integer :: nannihilate = 0
          ! The orbitals to be annihilated.
          integer, allocatable :: ann_orbs(:)
          ! The elements in the ilut representation where the occupation of the above orbs are encoded.
          integer, allocatable :: ann_elems(:)
          ! The positions of the bits in the bitstring representation where the above orbs are encoded.
          integer, allocatable :: ann_bits(:)

          ! The number of creation operators in the perturbation.
          integer :: ncreate = 0
          ! The orbitals to be created.
          integer, allocatable :: crtn_orbs(:)
          ! The elements in the ilut representation where the occupation of the above orbs are encoded.
          integer, allocatable :: crtn_elems(:)
          ! The positions of the bits in the bitstring representation where the above orbs are encoded.
          integer, allocatable :: crtn_bits(:)
      end type perturbation

      type(perturbation), allocatable :: pops_pert(:)

contains

    subroutine set_initial_global_data(ndets, ilut_list)

        use bit_rep_data, only: NIfTot, NIfDBO, extract_sign
        use Parallel_neci, only: iProcIndex, MPISumAll

        ! Take in a list of determinants and calculate and set all of the
        ! global data needed for the start of a FCIQMC calculation.

        integer(int64), intent(in) :: ndets
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot,ndets)

        integer :: i, run
        real(dp) :: real_sign(lenof_sign)
        character(*), parameter :: t_r = 'set_initial_global_data'

        TotParts = 0.0_dp
        NoAtHF = 0.0_dp

        ! First, find the population of the walkers in walker_list.
        do i = 1, ndets
            call extract_sign(ilut_list(:,i), real_sign)
            TotParts = TotParts + abs(real_sign)
            if ( all(ilut_list(0:NIfDBO,i) == iLutRef(0:NIfDBO)) ) NoAtHF = real_sign
        end do

        TotWalkers = ndets
        TotWalkersOld = TotWalkers
        call MPISumAll(TotWalkers,AllTotWalkers)
        AllTotWalkersOld = AllTotWalkers

        TotPartsOld = TotParts
        call MPISumAll(TotParts, AllTotParts)
        AllTotPartsOld = AllTotParts

        call MPISumAll(NoatHF, AllNoatHF)
        OldAllNoatHF = AllNoatHF

#ifdef __CMPLX
        OldAllAvWalkersCyc = sum(AllTotParts)
#else
        OldAllAvWalkersCyc = AllTotParts
#endif

        do run = 1, inum_runs
            OldAllHFCyc(run) = ARR_RE_OR_CPLX(AllNoatHF, run)
        end do

        AllNoAbortedOld(:) = 0.0_dp

        iter_data_fciqmc%tot_parts_old = AllTotParts
        
        ! Calculate the projected energy for this iteration.
        do run = 1, inum_runs
            if (ARR_RE_OR_CPLX(AllSumNoAtHF,run) /= 0) &
                ProjectionE(run) = AllSumENum(run) / ARR_RE_OR_CPLX(AllSumNoatHF,run)
        enddo 

        if (iProcIndex == iHFProc) then
            SumNoatHF(:) = AllSumNoatHF(:)
            SumENum(:) = AllSumENum(:)
            InstNoatHF(:) = NoatHF(:)

            if ( any(AllNoatHF /= NoatHF) ) then
                call stop_all(t_r, "HF particles spread across different processors.")
            endif
        endif

    end subroutine set_initial_global_data

END MODULE FciMCData
