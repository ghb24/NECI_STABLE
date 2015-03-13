module CalcData

    use constants
    use MemoryManager, only: TagIntType
    implicit none

    save

! Type containing logicals, which will tell various routines which
! routines to call when generating a subspace. This is used, for example,
! by semi-stochastic and trial wave function initialisation.
type subspace_in
  ! Just the Hartree-Fock determinant.
  logical :: tHF = .false.
  ! Use the most populated states in CurrentDets.
  logical :: tPops = .false.
  ! Read states from a file.
  logical :: tRead = .false.
  ! Use the space of all single and double excitations from the
  ! Hartree-Fock determinant (and also include the HF determinant).
  logical :: tDoubles = .false.
  ! Use all connections to the Hartree-Fock.
  logical :: tHFConn = .false.
  ! Use a CAS space.
  logical :: tCAS = .false.
  ! Use a RAS space.
  logical :: tRAS = .false.
  ! Use the determinants with the largest amplitudes in the MP1
  ! wave function.
  logical :: tMP1 = .false.
  ! Use the iterative approach described by Petruzielo et. al. (PRL 109, 230201).
  logical :: tOptimised = .false.
  ! Like the optimised space, but instead of diagonalising the space
  ! each iteration to find which states to keep, we keep the states
  ! with the lowest energies.
  logical :: tLowE = .false.
  ! Use the entire FCI space.
  logical :: tFCI = .false.
  ! Use the entire FCI space in Heisenberg model calculations.
  logical :: tHeisenbergFCI = .false.
end type subspace_in

LOGICAL :: TSTAR,TTROT,TGrowInitGraph
LOGICAL :: TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
LOGICAL :: TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER,tMP2Standalone
LOGICAL :: EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET
LOGICAL :: TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
LOGICAL :: TLADDER,TMC,TREADRHO,TRHOIJ,TBiasing,TMoveDets
LOGICAL :: TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph
LOGICAL :: TInitStar,TNoSameExcit,TLanczos,TStarTrips
LOGICAL :: TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
LOGICAL ::THDiag,TMCStar,TReadPops,TBinCancel,TFCIMC,TMCDets,tDirectAnnihil
LOGICAL :: tDetermProj, tFTLM, tSpecLanc, tExactSpec, tExactDiagAllSym
LOGICAL :: TFullUnbias,TNoAnnihil,tStartMP1
LOGICAL :: TRhoElems,TReturnPathMC,TSignShift
LOGICAL :: THFRetBias,TProjEMP2,TFixParticleSign
LOGICAL :: TStartSinglePart,TRegenExcitgens
LOGICAL :: TUnbiasPGeninProjE, tCheckHighestPopOnce
LOGICAL :: tCheckHighestPop,tRestartHighPop,tChangeProjEDet
LOGICAL :: tRotoAnnihil,tSpawnAsDet
LOGICAL :: tTruncCAS,tTruncInitiator,tAddtoInitiator    !Truncation the FCIMC excitation space by CAS
LOGICAL :: tInitIncDoubs,tWalkContGrow,tAnnihilatebyRange
logical :: tReadPopsRestart, tReadPopsChangeRef, tInstGrowthRate
logical :: tAllRealCoeff, tUseRealCoeffs
logical :: tRealSpawnCutoff
logical :: tRealCoeffByExcitLevel
integer :: RealCoeffExcitThresh
real(dp) :: RealSpawnCutoff, OccupiedThresh, InitiatorOccupiedThresh
logical :: tEnhanceRemainder
logical :: tInitOccThresh
logical :: tRPA_QBA     !RPA calculation with QB approximation
logical :: tStartCAS    !Start FCIMC dynamic with walkers distributed according to CAS diag.
logical :: tShiftonHFPop    !Adjust shift in order to keep the population on HF constant, rather than total pop.

! Base hash values only on spatial orbitals
! --> All dets with same spatial structure on the same processor.
logical :: tSpatialOnlyHash

! Do we truncate spawning based on the number of unpaired electrons
logical :: tTruncNOpen
integer :: trunc_nopen_max

logical :: tMaxBloom    !If this is on, then we only print out a bloom warning if it is the biggest to date.

INTEGER :: NWHTAY(3,10),NPATHS,NoMoveDets,NoMCExcits,NShiftEquilSteps
INTEGER :: NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED,HApp,iFullSpaceIter
INTEGER :: IMCSTEPS,IEQSTEPS,MDK(5),Iters,NDets,iDetGroup
INTEGER :: CUR_VERT,NHISTBOXES,I_P,LinePoints,iMaxExcitLevel
INTEGER :: NMCyc,StepsSft,CLMax
INTEGER :: NEquilSteps,iSampleRDMIters
real(dp) :: InitialPart
real(dp), allocatable :: InitialPartVec(:)
INTEGER :: OccCASorbs,VirtCASorbs,iAnnInterval
integer :: iPopsFileNoRead, iPopsFileNoWrite,iRestartWalkNum
real(dp) :: iWeightPopRead
integer :: MaxWalkerBloom   !Max number of walkers allowed in one bloom before reducing tau
INTEGER(int64) :: HFPopThresh
real(dp) :: InitWalkers, maxnoathf, InitiatorWalkNo

! The average number of excitations to be performed from each walker.
real(dp) :: AvMCExcits

integer :: iReadWalkersRoot !The number of walkers to read in on the head node in each batch during a popsread

real(dp) :: g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
real(dp) :: G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10)
real(dp) :: BETAP,RHOEPSILON,DBETA,STARCONV,GraphBias
real(dp) :: GrowGraphsExpo,Tau,SftDamp,ScaleWalkers
real(dp) :: PRet,FracLargerDet
real(dp) :: MemoryFacPart
real(dp) :: MemoryFacSpawn,SinglesBias,TauFactor,StepsSftImag

real(dp) :: MemoryFacInit

real(dp), allocatable, target :: DiagSft(:)

real(dp) :: GraphEpsilon
real(dp) :: PGenEpsilon
real(dp) :: InputTargetGrowRate
integer(int64) :: InputTargetGrowRateWalk
real(dp), allocatable :: TargetGrowRate(:)
! Number of particles before targetgrowrate kicks in
integer(int64), allocatable :: TargetGrowRateWalk(:)
integer(int64) :: iExitWalkers  !Exit criterion, based on total walker number


!// additional from NECI.F
INTEGER, Allocatable :: MCDet(:)
INTEGER(TagIntType) :: tagMCDet=0
real(dp) :: RHOEPS ! calculated from RHOEPSILON

!// set if we include no triple-excitations as the 3rd vertex in 3+ vertex graphs.
LOGICAL :: lNoTriples

LOGICAL tUseProcsAsNodes  !Set if we treat each processor as its own node.
INTEGER iLogicalNodeSize  !An alternative to the above, create logical nodes of at most this size.
                          ! 0 means use physical nodes.

    logical :: tJumpShift

! Perform a Davidson calculation if true.
logical :: tDavidson
! Should the HF determinant be put on its own processor?
logical :: tUniqueHFNode

! Options relating to the semi-stochastic code.
logical :: tSemiStochastic ! Performing a semi-stochastic simulation if true.
! Input type describing which space(s) type to use.
type(subspace_in) :: ss_space_in
! Options regarding splitting the space into core and non-core elements. Needed, for example when performing a
! semi-stochastic simulation, to specify the deterministic space.
logical :: tCSFCore ! Use CSFs for the core states.
logical :: tSparseCoreHamil ! Use a sparse representation of the core Hamiltonian.
! When using a CAS deterministic space, these integers store the number of orbitals above and below the Fermi energy to
! include in the CAS active space (the occupied and virtual orbitals).
integer :: OccDetermCASOrbs
integer :: VirtDetermCASOrbs

! This type refers to data that specifies optimised core spaces. It is passed to
! the routine generate_optimised_core. See this routine for an explanation.
type opt_space_data
    ! The number of generation loops for the algorithm to perform. 
    integer :: ngen_loops
    ! If true then perform cutoff of determinants each iteration by using
    ! amplitudes les than the values in cutoff_amps, else perform cutoff using the
    ! number of determinants in cutoff_nums.
    logical :: tAmpCutoff
    ! If tAmpCutoff is .true. then this will be allocated and used to perform the
    ! cutoff of determinants each iteration.
    real(dp), allocatable :: cutoff_amps(:)
    ! If tAmpCutoff is .false. then this will be allocated and used to perform the
    ! cutoff of determinants each iteration.
    integer, allocatable :: cutoff_nums(:)
end type

! Optimised space data for generating a semi-stohastic space.
type(opt_space_data) :: determ_opt_data
! Optimised space data for generating a trial wave function space.
type(opt_space_data) :: trial_opt_data

! If this is true then set a limit on the maximum deterministic space size.
logical :: tLimitDetermSpace
! This is maximum number of elements in the deterministic space, if tLimitDetermSpace is true.
integer :: max_determ_size
! The number of states to use from a POPSFILE for the core space.
integer :: n_core_pops
! This option gives the maximum excitation level to go up to when generating the low energy deterministic space.
integer :: low_e_core_excit
! This integer specifies the number of states to keep for each iteration of the low energy core generation.
integer :: low_e_core_num_keep
! When using tLowECore, if this option is true then all doubles will be kept.
logical :: tLowECoreAllDoubles
! When using the tMP1Core option, this specifies how many determinants to keep.
integer :: semistoch_mp1_ndets

! If this is non-zero then we turn semi-stochastic semistoch_shift_iter
! iterations after the shift starts to vary.
integer :: semistoch_shift_iter

! If true then, if using a deterministic space of all singles and doubles, no
! stochastic spawning will be attempted from the Hartree-Fock. This is allowed
! because all 'spawnings' from the Hartree-Fock in this case will be
! deterministic.
logical :: tDetermHFSpawning

! Options relating to the trial wavefunction.
logical :: tTrialWavefunction ! Use a trial wavefunction-based energy estimator.
! Input type describing which space(s) type to use.
type(subspace_in) :: trial_space_in
! When using a CAS trial space, these integers store the number of orbitals above and below the Fermi energy to
! include in the CAS active space (the occupied and virtual orbitals).
integer :: OccTrialCASOrbs
integer :: VirtTrialCASOrbs
! If this is true then set a limit on the maximum trial space size.
logical :: tLimitTrialSpace
! This is maximum number of elements in the trial space, if tLimitDetermSpace is true.
integer :: max_trial_size
! The number of states to use from a POPSFILE for the trial space.
integer :: n_trial_pops
! This option gives the maximum excitation level to go up to when generating the low energy trial space.
integer :: low_e_trial_excit
! This integer specifies the number of states to keep for each iteration of the low energy trial generation.
integer :: low_e_trial_num_keep
! When using tLowETrial, if this option is true then all doubles will be kept.
logical :: tLowETrialAllDoubles
! When using the tMP1Trial option, this specifies how many determinants to keep.
integer :: trial_mp1_ndets

! True if running a kp-fciqmc calculation.
logical :: tKP_FCIQMC

! If this is true then we only start varying shift when we get *below* the target population, rather than above it.
! This is useful when we want to start from a large and high-energy population and let many walkers quickly die with
! a constant shift (i.e. finite-temperature calculations).
logical :: tLetInitialPopDie

! Calculate the norms of the *unperturbed* POPSFILE wave functions and output them to a file.
logical :: tWritePopsNorm
real(dp) :: pops_norm
integer :: pops_norm_unit

! What is the maximum energy, above which all particles are treated as
! initiators
real(dp) :: InitiatorCutoffEnergy, InitiatorCutoffWalkNo

! Should we aggregate particle counts across all of the simulations in
! determining which sites are initiators (in system-replica mode).
logical :: tMultiReplicaInitiators = .false.

! Do we make sites into initiators if they have survived more than a certain
! period of time?
logical :: tSurvivalInitiatorThreshold, tSurvivalInitMultThresh
real(dp) :: im_time_init_thresh, init_survival_mult

! Do we make sites into initiators depending on the number of spawns to them?
logical :: tSpawnCountInitiatorThreshold
integer :: init_spawn_thresh

! Are we orthogonalising replicas?
logical :: tOrthogonaliseReplicas, tReplicaSingleDetStart
integer :: orthogonalise_iter

! If true then a hash table is kept for the spawning array and is used when
! new spawnings are added to the spawned list, to prevent adding the same
! determinant multiple times.
logical :: use_spawn_hash_table

end module CalcData
