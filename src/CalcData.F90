module CalcData

    use constants
    use ras_data, only: ras_parameters
    use MemoryManager, only: TagIntType
    implicit none

    save

! This type refers to data that specifies 'optimised' spaces. It is passed to
! the routine generate_optimised_core. See this routine for an explanation.
type opt_space_data
    ! The number of generation loops for the algorithm to perform. 
    integer :: ngen_loops = 1
    ! If true then perform cutoff of determinants each iteration by using
    ! amplitudes les than the values in cutoff_amps, else perform cutoff using the
    ! number of determinants in cutoff_nums.
    logical :: tAmpCutoff = .false.
    ! If tAmpCutoff is .true. then this will be allocated and used to perform the
    ! cutoff of determinants each iteration.
    real(dp), allocatable :: cutoff_amps(:)
    ! If tAmpCutoff is .false. then this will be allocated and used to perform the
    ! cutoff of determinants each iteration.
    integer, allocatable :: cutoff_nums(:)
end type

! Type containing information for routines which generate subspaces. This is
! used, for example, by semi-stochastic and trial wave function initialisation.
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

    ! When using a CAS deterministic space, these integers store the number of orbitals above and below the Fermi energy to
    ! include in the CAS active space (the occupied and virtual orbitals).
    integer :: occ_cas
    integer :: virt_cas

    ! Optimised space data for generating a semi-stohastic space.
    type(opt_space_data) :: opt_data

    ! Paremeters for RAS space calculations.
    type(ras_parameters) :: ras
    
    ! If this is true then set a limit on the maximum deterministic space size.
    logical :: tLimitSpace = .false.
    ! This is maximum number of elements in the deterministic space, if tLimitDetermSpace is true.
    integer :: max_size = 0
    ! The number of states to use from a POPSFILE for the core space.
    integer :: npops = 0
    ! When using the tMP1Core option, this specifies how many determinants to keep.
    integer :: mp1_ndets = 0
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
real(dp) :: PRet,FracLargerDet,pop_change_min
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

    logical :: tJumpShift, tPopsJumpShift

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

! If this is non-zero then we turn semi-stochastic semistoch_shift_iter
! iterations after the shift starts to vary.
integer :: semistoch_shift_iter

! If true then, if using a deterministic space of all singles and doubles, no
! stochastic spawning will be attempted from the Hartree-Fock. This is allowed
! because all 'spawnings' from the Hartree-Fock in this case will be
! deterministic.
logical :: tDetermHFSpawning

! Options relating to the trial wavefunction.

! If true at a given point during a simulation then we are currently
! calculating trial wave function-based energy estimates.
logical :: tTrialWavefunction
! How many excited states to calculate in the trial space, for the
! trial wave functions estimates
integer :: ntrial_ex_calc = 0
! Input type describing which space(s) type to use.
type(subspace_in) :: trial_space_in

! If true then start using a trial estimator later on in the calculation.
logical :: tStartTrialLater = .false.
! How many iterations after the shift starts to vary should be turn on the use
! of trial estimators?
integer :: trial_shift_iter

! If false then create the trial wave function by diagonalising the
! Hamiltonian in the trial subspace.
! If true then create the trial wave function by taking the weights from the
! QMC simulation in the trial subspace, at the point the that the trial
! wave function is turned on.
logical :: qmc_trial_wf = .false.

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
logical :: tOrthogonaliseSymmetric
integer :: orthogonalise_iter
! Information on a trial space to create trial excited states with.
type(subspace_in) :: init_trial_in

! If true then a hash table is kept for the spawning array and is used when
! new spawnings are added to the spawned list, to prevent adding the same
! determinant multiple times.
logical :: use_spawn_hash_table

logical :: tMultipleInitialRefs = .false.
integer, allocatable :: initial_refs(:,:)

! Array to specify how to reorder the trial states (which are by default
! ordered by the energy in the trial space).
! First the trial states for the energy estimates:
integer, allocatable :: trial_est_reorder(:)
! And also the trial states used for the intial states:
integer, allocatable :: trial_init_reorder(:)

! If true then, when using the orthogonalise-replicas option, print out the
! overlaps between replicas in a separate file.
logical :: tPrintReplicaOverlaps = .true.

! Keep track of when the calculation began (globally)
real(sp) :: s_global_start

! Use continuous time FCIQMC
logical :: tContTimeFCIMC, tContTimeFull
real(dp) :: cont_time_max_overspawn

! Are we doing an mneci run where each state is represented by two FCIQMC
! replicas?
logical :: tPairedReplicas = .false.

! If true then swap the sign of the FCIQMC wave function if the sign of the
! Hartree-Fock population becomes negative.
logical :: tPositiveHFSign = .false.

end module CalcData
