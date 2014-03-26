module CalcData

    use constants, only: dp,int64, lenof_sign, inum_runs, n_int
    use MemoryManager, only: TagIntType
    implicit none

    save

LOGICAL :: TSTAR,TTROT,TGrowInitGraph
LOGICAL :: TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
LOGICAL :: TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER,tMP2Standalone
LOGICAL :: EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET
LOGICAL :: TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
LOGICAL :: TLADDER,TMC,TREADRHO,TRHOIJ,TBiasing,TMoveDets
LOGICAL :: TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph
LOGICAL :: TInitStar,TNoSameExcit,TLanczos,TStarTrips
LOGICAL :: TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
LOGICAL ::THDiag,TMCStar,TReadPops,TBinCancel,TFCIMC,TMCDets,tDirectAnnihil,tCCMC
LOGICAL :: tDetermProj, tFTLM
LOGICAL :: TFullUnbias,TNoAnnihil,tStartMP1
LOGICAL :: TRhoElems,TReturnPathMC,TSignShift
LOGICAL :: THFRetBias,TProjEMP2,TFixParticleSign
LOGICAL :: TStartSinglePart,TRegenExcitgens
LOGICAL :: TUnbiasPGeninProjE, tCheckHighestPopOnce
LOGICAL :: tCheckHighestPop,tRestartHighPop,tChangeProjEDet
LOGICAL :: tRotoAnnihil,tRegenDiagHEls,tSpawnAsDet,tFindGroundDet
LOGICAL :: tTruncCAS,tTruncInitiator,tDelayTruncInit,tAddtoInitiator    !Truncation the FCIMC excitation space by CAS
LOGICAL :: tInitIncDoubs,tWalkContGrow,tAnnihilatebyRange,tRetestAddtoInit
logical :: tReadPopsRestart, tReadPopsChangeRef, tInstGrowthRate
logical :: tAllRealCoeff, tUseRealCoeffs
logical :: tRealSpawnCutoff
logical :: tRealCoeffByExcitLevel
integer :: RealCoeffExcitThresh
real(dp) :: RealSpawnCutoff, OccupiedThresh
logical :: tEnhanceRemainder
logical :: tRPA_QBA     !RPA calculation with QB approximation
logical :: tStartCAS    !Start FCIMC dynamic with walkers distributed according to CAS diag.
logical :: tPopsMapping !Map popsfile from smaller basis onto larger basis
logical :: tShiftonHFPop    !Adjust shift in order to keep the population on HF constant, rather than total pop.

! Base hash values only on spatial orbitals
! --> All dets with same spatial structure on the same processor.
logical :: tSpatialOnlyHash

! Do we allow walkers to survive (in the initiator approx.) if a determinant
! with the same spatial configuration is an initiator?
logical :: tSpawnSpatialInit

!These options mean that only initiators can spawn walkers.
!tSpawn_Only_Init_Grow means that this option is removed once variable shift is entered.
logical :: tSpawn_Only_Init,tSpawn_Only_Init_Grow

! Do we truncate spawning based on the number of unpaired electrons
logical :: tTruncNOpen
integer :: trunc_nopen_max

logical :: tMaxBloom    !If this is on, then we only print out a bloom warning if it is the biggest to date.

INTEGER :: NWHTAY(3,10),NPATHS,NoMoveDets,NoMCExcits,IterTruncInit,NShiftEquilSteps
INTEGER :: NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED,HApp,iFullSpaceIter
INTEGER :: IMCSTEPS,IEQSTEPS,MDK(5),Iters,NDets,iDetGroup
INTEGER :: CUR_VERT,NHISTBOXES,I_P,LinePoints,iMaxExcitLevel
INTEGER :: NMCyc,StepsSft,CLMax
INTEGER :: NEquilSteps
real(dp) :: InitialPart
real(dp), dimension(lenof_sign) :: InitialPartVec
INTEGER :: OccCASorbs,VirtCASorbs,iAnnInterval
integer :: iPopsFileNoRead, iPopsFileNoWrite,iWeightPopRead,iRestartWalkNum
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
real(dp) :: GrowMaxFactor,CullFactor,PRet,FracLargerDet
real(dp) :: MemoryFacPart,MemoryFacAnnihil
real(dp) :: MemoryFacSpawn,SinglesBias,TauFactor,StepsSftImag

real(dp) :: MemoryFacInit

real(dp), dimension(inum_runs), target :: DiagSft

real(dp) :: GraphEpsilon
real(dp) :: PGenEpsilon
real(dp), dimension(inum_runs) :: TargetGrowRate
integer(int64), dimension(inum_runs) :: TargetGrowRateWalk    !Number of walkers before targetgrowrate kicks in
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

logical :: tContinueAfterMP2 ! UEG option only
    logical :: tJumpShift

! If true, the initiator threshold for a given state will be equal to the (minimum) number of spawning
! attempts away from the core space it took to create the state.
logical :: tVaryInitThresh

! Perform a Davidson calculation if true.
logical :: tDavidson
! Should the HF determinant be put on its own processor?
logical :: tUniqueHFNode

! Options relating to the semi-stochastic code.
logical :: tSemiStochastic ! Performing a semi-stochastic simulation if true.
! Options regarding splitting the space into core and non-core elements. Needed, for example when performing a
! semi-stochastic simulation, to specify the deterministic space.
logical :: tCSFCore ! Use CSFs for the core states.
logical :: tOptimisedCore ! Generate an optimised deterministic space by diagonalising part of the space.
logical :: tDoublesCore ! Use single and double excitations for the core states.
logical :: tCASCore ! Use Determinants where orbitals within an active space can differ from the Hartree-Fock for core states.
logical :: tRASCore ! Use a RAS space for the core space (see ras.F90 for definition).
logical :: tPopsCore ! Use the most populated states from a POPSFILE for the core space.
logical :: tReadCore ! Read in the core space from the CORESPACE file.
! Like the optimised core space, but instead of diagonalising the space each iteration to find which states to keep, we
! keep the states with the lowest energies.
logical :: tLowECore
! If a number of determinants smaller than the number of singles and doubles to the HF is requested, then the amplitude
! of the determinants in the MP1 wave function will be used to determine which to keep. Otherwise all singles and
! doubles are kept.
logical :: tMP1Core 
! Use the entire Hilbert space as the core space.
logical :: tFCICore
logical :: tHeisenbergFCICore
logical :: tSparseCoreHamil ! Use a sparse representation of the core Hamiltonian.
! cas_determ_bitmask has all bits that refer to the active space set, and all other bits unset.
! cas_not_determ_bitmask is simply the result after the not operation is applied to cas_determ_bitmask.
integer(n_int), allocatable, dimension(:) :: cas_determ_bitmask
integer(n_int), allocatable, dimension(:) :: cas_determ_not_bitmask
! Bitmasks with all bits corresponding to orbitals in RAS1 and RAS3, repectively, set.
integer(n_int), allocatable, dimension(:) :: core_ras1_bitmask
integer(n_int), allocatable, dimension(:) :: core_ras3_bitmask
! When using a CAS deterministic space, these integers store the number of orbitals above and below the Fermi energy to
! include in the CAS active space (the occupied and virtual orbitals).
integer :: OccDetermCASOrbs
integer :: VirtDetermCASOrbs
! When using tOptimisedCore, this specifies the minimum amplitude that a basis state should have in the ground state
! (of the deterministic space generated) in order to be kept in the next deterministic space. The first component
! refers to the first iteration, the second component to the second iteration...
real(dp), allocatable, dimension(:) :: determ_space_cutoff_amp
! When using tOptimisedCore, this specifies how many basis states (of the deterministic space generated) should be kept
! in the next deterministic space. The first component refers to the first iteration, the second component to the second
! iteration... i.e. if determ_space_cutoff_num(1) = 5 then in the first iteration, the 5 basis states with the largest
! amplitudes in the ground state are kept.
integer, allocatable, dimension(:) :: determ_space_cutoff_num
! If this logical is true, then the cutoff criterion is done using the amplitude. If false, it is done using a fixed number.
logical :: tDetermAmplitudeCutoff
! When using the optimised core option for semi-stochastic simulations, this option specifies how many iterations of
! the generation procedure should be performed.
integer :: num_det_generation_loops
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

! Options relating to the trial wavefunction.
logical :: tTrialWavefunction ! Use a trial wavefunction-based energy estimator.
logical :: tDoublesTrial ! Use single and double exciations for the trial space.
logical :: tCASTrial ! Use a CAS space for the trial space.
logical :: tOptimisedTrial ! Generate an optimised trial space by diagonalisaing part of the space.
! As for determ_space_cutoff_amp and determ_space_cutoff_num above, but the following two quantities refer to the trial space
! generation rather than the deterministic space generation.
logical :: tReadTrial ! Read in the trial space from the TRIALSPACE file.
logical :: tPopsTrial ! Use the most populated states from a POPSFILE for the trial space.
! Like the optimised trial space, but instead of diagonalising the space each iteration to find which states to keep, we
! keep the states with the lowest energies.
logical :: tLowETrial
! If a number of determinants smaller than the number of singles and doubles to the HF is requested, then the amplitude
! of the determinants in the MP1 wave function will be used to determine which to keep. Otherwise all singles and
! doubles are kept.
logical :: tMP1Trial
! Use the entire Hilbert space as the trial space.
logical :: tFCITrial
logical :: tHeisenbergFCITrial
real(dp), allocatable, dimension(:) :: trial_space_cutoff_amp
integer, allocatable, dimension(:) :: trial_space_cutoff_num
! When using a CAS trial space, these integers store the number of orbitals above and below the Fermi energy to
! include in the CAS active space (the occupied and virtual orbitals).
integer :: OccTrialCASOrbs
integer :: VirtTrialCASOrbs
! If this logical is true, then the cutoff criterion is done using the amplitude. If false, it is done using a fixed number.
logical :: tTrialAmplitudeCutoff
! As for num_determ_generation_loops above, but here for the trial wavefunction generation.
integer :: num_trial_generation_loops
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

! If true, apply creation and annihilation operators (using the orbitals in the arrays below) to the wave function
! in the POPSFILE before starting from it.
logical :: tPerturbPops
! The number of annihilation operators in the perturbation.
integer :: n_pops_annihilate
! The orbitals to be annihilated.
integer, allocatable :: annihilate_orbs(:)
! The elements in the ilut representation where the occupation of the above orbs are encoded.
integer, allocatable :: annihilate_elems(:)
! The positions of the bits in the bitstring representation where the above orbs are encoded.
integer, allocatable :: annihilate_bits(:)
! The number of creation operators in the perturbation.
integer :: n_pops_creation
! The orbitals to be created.
integer, allocatable :: creation_orbs(:)
! The elements in the ilut representation where the occupation of the above orbs are encoded.
integer, allocatable :: creation_elems(:)
! The positions of the bits in the bitstring representation where the above orbs are encoded.
integer, allocatable :: creation_bits(:)

! Calculate the norms of the *unperturbed* POPSFILE wave functions and output them to a file.
logical :: tWritePopsNorm
real(dp) :: pops_norm(lenof_sign)
integer :: pops_norm_unit

end module CalcData
