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
    ! Automatically choosing 10% of the total initiator space, if this number
    ! is larger than 50000, then use npops = 50000
    logical :: tPopsAuto = .false.
    ! Read states from a file.
    logical :: tRead = .false.
    ! Use the space of all single and double (+triple) excitations from the
    ! Hartree-Fock determinant (and also include the HF determinant).
    logical :: tDoubles = .false.
    logical :: tTriples = .false.
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
    ! Actually use the space of all connections to the chosen space
    logical :: tAllConnCore = .false.
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
    ! If true then the space generated from the pops-core option may be very
    ! slightly unoptimal, but should be very good, and sometimes exactly the
    ! same as when this logical is false. Only applies to the pops-* options.
    logical :: tApproxSpace = .false.
    ! if tApproxSpace is true, then nApproxSpace times target pops-core space is the size of the array
    ! kept on each processor, 1 =< nApproxSpace =< nProcessors. The larger nApproxSpace, the more
    ! memory is consumed and the slower (but more accurate) is the semi-stochastic initialisation
    ! (see subroutine generate_space_most_populated).
    integer :: nApproxSpace = 10
    ! When using the tMP1Core option, this specifies how many determinants to keep.
    integer :: mp1_ndets = 0

    character(255) :: read_filename
end type subspace_in

LOGICAL :: TSTAR,TTROT,TGrowInitGraph
LOGICAL :: TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
LOGICAL :: TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER,tMP2Standalone
LOGICAL :: EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET
LOGICAL :: TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
LOGICAL :: TLADDER,TMC,TREADRHO,TRHOIJ,TBiasing,TMoveDets
LOGICAL :: TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph
LOGICAL :: TInitStar,TNoSameExcit,TLanczos,TStarTrips, tFCIDavidson
LOGICAL :: TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
LOGICAL ::THDiag,TMCStar,TReadPops,TBinCancel,TFCIMC,TMCDets,tDirectAnnihil
LOGICAL :: tDetermProj, tFTLM, tSpecLanc, tExactSpec, tExactDiagAllSym
LOGICAL :: tDetermProjApproxHamil
LOGICAL :: TFullUnbias,TNoAnnihil,tStartMP1
LOGICAL :: TRhoElems,TReturnPathMC,TSignShift
LOGICAL :: THFRetBias,TProjEMP2,TFixParticleSign
LOGICAL :: TStartSinglePart,TRegenExcitgens
LOGICAL :: TUnbiasPGeninProjE, tCheckHighestPopOnce
LOGICAL :: tCheckHighestPop,tRestartHighPop,tChangeProjEDet
LOGICAL :: tRotoAnnihil,tSpawnAsDet
LOGICAL :: tTruncCAS ! Truncation of the FCIMC excitation space by a CAS
logical :: tTruncInitiator, tAddtoInitiator, tInitCoherentRule, tGlobalInitFlag
logical :: tEN2, tEN2Init, tEN2Truncated, tEN2Started, tEN2Rigorous

LOGICAL :: tSeniorInitiators !If a det. has lived long enough (called a senior det.), it is added to the initiator space.
LOGICAL :: tWalkContGrow,tAnnihilatebyRange
logical :: tReadPopsRestart, tReadPopsChangeRef, tInstGrowthRate
logical :: tL2GrowRate
logical :: tAllRealCoeff, tUseRealCoeffs
logical :: tRealSpawnCutoff
logical :: tRealCoeffByExcitLevel
integer :: RealCoeffExcitThresh
real(dp) :: RealSpawnCutoff, OccupiedThresh
logical :: tRPA_QBA     !RPA calculation with QB approximation
logical :: tStartCAS    !Start FCIMC dynamic with walkers distributed according to CAS diag.
logical :: tShiftonHFPop    !Adjust shift in order to keep the population on HF constant, rather than total pop.
logical :: tSpecifiedTau
logical :: tInitializeCSF
real(dp) :: S2Init
logical :: tFixedN0 !Fix the reference population by using projected energy as shift.
logical :: tTrialShift !Fix the overlap with trial wavefunction by using trial energy as shift.
logical :: tSkipRef(1:inum_runs_max) !Skip spawing onto reference det and death/birth on it. One flag for each run.
logical :: tFixTrial(1:inum_runs_max) !Fix trial overlap by determinstically updating one det. One flag for each run.
integer :: N0_Target !The target reference population in fixed-N0 mode
real(dp) :: TrialTarget !The target for trial overlap in trial-shift mode
logical :: tAdaptiveShift !Whether any of the adaptive shift schemes is used
logical :: tCoreAdaptiveShift = .false. ! Whether the adaptive shift is also applied to the corespace
logical :: tLinearAdaptiveShift !Make shift depends on the population linearly
real(dp) :: LAS_Sigma !Population which below the shift is set to zero
real(dp) :: LAS_F1 !Shift modification factor at AdaptiveShiftSigma
real(dp) :: LAS_F2 !Shift modification factor at InitiatorWalkNo
logical :: tAutoAdaptiveShift !Let the modification factor of adaptive shift be computed autmatically
real(dp) :: AAS_Thresh !Number of spawn under which below the shift is set to zero
real(dp) :: AAS_Expo !Exponent of the modification factor, value 1 is default. values 0 means going back to full shift.
real(dp) :: AAS_Cut  !The modification factor should never go below this.
logical :: tAAS_MatEle !Use the magnitude of |Hij| in the modifcation factor i.e. sum_{accepted} |H_ij| / sum_{all attempts} |H_ij|
logical :: tAAS_MatEle2 !Use the weight |Hij|/(Hjj-E) in the modifcation factor
logical :: tAAS_MatEle3 !Same as MatEle2 but use weight of one for accepted moves.
logical :: tAAS_MatEle4 !Same as MatEle2 but use E_0 in the weight of accepted moves.
real(dp) :: AAS_DenCut !Threshold on the denominators of MatEles
real(dp) :: AAS_Const
logical :: tExpAdaptiveShift !Make the shift depends on the population exponentialy
real(dp) :: EAS_Scale !Scale parameter of exponentail adaptive shift
logical :: tAS_TrialOffset !Whether the adaptive shift scheme should be applied with respect to trial-wf energy not HF energy
real(dp) :: ShiftOffset ! An offest for the adaptive shift
! Giovannis option for using only initiators for the RDMs (off by default)
logical :: tOutputInitsRDM = .false.
logical :: tNonInitsForRDMs = .true.
logical :: tNonVariationalRDMs = .false.
! Adaptive shift RDM correction using initiators as reference
logical :: tInitsRDMRef, tInitsRDM
! Base hash values only on spatial orbitals
! --> All dets with same spatial structure on the same processor.
logical :: tSpatialOnlyHash

! if all determinants are stored to prevent the need for conversion each iteration
logical :: tStoredDets

! Do we truncate spawning based on the number of unpaired electrons
logical :: tTruncNOpen
integer :: trunc_nopen_max

! introduce a new truncation scheme based on the difference of seniority
! compared to the reference determinant
logical :: t_trunc_nopen_diff = .false.
integer :: trunc_nopen_diff = 0

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
real(dp) :: MaxWalkerBloom   !Max number of walkers allowed in one bloom before reducing tau
INTEGER(int64) :: HFPopThresh
real(dp) :: InitWalkers, maxnoathf, InitiatorWalkNo, ErrThresh
! Options for dynamic rescaling of spawn attempts + blooms
logical :: tScaleBlooms = .false.
real(dp) :: max_allowed_spawn

real(dp) :: SeniorityAge !A threshold on the life time of a determinat (measured in its halftime) to become a senior determinant.

! The average number of excitations to be performed from each walker.
real(dp) :: AvMCExcits
! Optionally: allow this number to change during runtime
logical :: tDynamicAvMCEx
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
! for consistency with forgetting the walkcontgrow keyword and hdf5 read-in
! use a temporary storage of the read-in diags-shift
real(dp), allocatable :: hdf5_diagsft(:)

real(dp) :: GraphEpsilon
real(dp) :: PGenEpsilon
real(dp) :: InputTargetGrowRate
integer(int64) :: InputTargetGrowRateWalk
real(dp), allocatable :: TargetGrowRate(:)
! Number of particles before targetgrowrate kicks in
integer(int64), allocatable :: TargetGrowRateWalk(:)
integer(int64) :: iExitWalkers  !Exit criterion, based on total walker number

! Lanczos initialisation of the wavefunctions
logical :: t_lanczos_init
logical :: t_lanczos_store_vecs
logical :: t_lanczos_orthogonalise
logical :: t_force_lanczos
integer :: lanczos_max_restarts
integer :: lanczos_max_vecs
integer :: lanczos_energy_precision
integer :: lanczos_ritz_overlap_precision


!// additional from NECI.F
INTEGER, Allocatable :: MCDet(:)
INTEGER(TagIntType) :: tagMCDet=0
real(dp) :: RHOEPS ! calculated from RHOEPSILON

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
logical :: tDynamicCoreSpace, tStaticCore, tIntervalSet ! update the corespace
integer :: coreSpaceUpdateCycle, semistochStartIter
! Input type describing which space(s) type to use.
type(subspace_in) :: ss_space_in

! Options regarding splitting the space into core and non-core elements. Needed, for example when performing a
! semi-stochastic simulation, to specify the deterministic space.
logical :: tSparseCoreHamil ! Use a sparse representation of the core Hamiltonian.

! If this is non-zero then we turn semi-stochastic on semistoch_shift_iter
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

! if we want to choose a specific excited states as the trial wf, if we
! have a reasonable estimate. this must be done for all replicas if
! multiple are used
logical :: t_choose_trial_state = .false.
integer, allocatable :: trial_excit_choice(:)

! Input type describing which space(s) type to use.
type(subspace_in) :: trial_space_in

! If true then start using a trial estimator later on in the calculation.
logical :: tStartTrialLater = .false.
! How many iterations after the shift starts to vary should be turn on the use
! of trial estimators?
integer :: trial_shift_iter

! Update the trial wf?
logical :: tDynamicTrial
integer :: trialSpaceUpdateCycle

! If false then create the trial wave function by diagonalising the
! Hamiltonian in the trial subspace.
! If true then create the trial wave function by taking the weights from the
! QMC simulation in the trial subspace, at the point the that the trial
! wave function is turned on.
logical :: qmc_trial_wf = .false.

! Define a space in which all determinants are initiators
logical :: tInitiatorSpace
type(subspace_in) :: i_space_in

! If true then initiators can only be those determinants in the defined fixed space.
logical :: tPureInitiatorSpace

! Run FCIQMC in the truncated space of all connections to the initiator space
logical :: tAllConnsPureInit

! Allow all spawns with (no) sign change
! The modi here are:  0, no changes to initiator approx are made
!                     >0 (commonly 1), same-sign spawns are always allowed
!                     <0 (commonly -1), opp. sign spawns are always allowed
integer :: allowedSpawnSign = 0

! If this is true, don't allow non-initiators to spawn to another non-initiator,
! even if it is occupied.
logical :: tSimpleInit

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

! Are we orthogonalising replicas?
logical :: tOrthogonaliseReplicas, tReplicaSingleDetStart
logical :: tOrthogonaliseSymmetric
integer :: orthogonalise_iter

! test reintroducing an overlap.
logical :: t_test_overlap = .false.
real(dp) :: overlap_eps = 1.0e-5_dp
integer :: n_stop_ortho = -1
! Information on a trial space to create trial excited states with.

type(subspace_in) :: init_trial_in

! Start wave function from solutions with a trial space
logical :: tTrialInit

! If true then a hash table is kept for the spawning array and is used when
! new spawnings are added to the spawned list, to prevent adding the same
! determinant multiple times.
logical :: use_spawn_hash_table

! Used when the user specifies multiple reference states manually, using the
! multiple-initial-refs option, similar to the definedet option but for more
! than one states, when using multiple replicas.
logical :: tMultipleInitialRefs = .false.
integer, allocatable :: initial_refs(:,:)

! As for the mutliple-initial-refs options above, but the
! multiple-initial-states option allows the user to manually specify the
! actual starting states, rather than the starting reference functions.
logical :: tMultipleInitialStates = .false.
integer, allocatable :: initial_states(:,:)

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
real(dp) :: s_global_start

! Use continuous time FCIQMC
logical :: tContTimeFCIMC, tContTimeFull
real(dp) :: cont_time_max_overspawn

! Are we doing an mneci run where each state is represented by two FCIQMC
! replicas?
logical :: tPairedReplicas = .false.

! Calculate and print estimates which use the replica approach to a file
logical :: tReplicaEstimates

logical :: tSetInitFlagsBeforeDeath

! If true then swap the sign of the FCIQMC wave function if the sign of the
! Hartree-Fock population becomes negative.
logical :: tPositiveHFSign = .false.

! If true, then set the initial shift for each replica (in jobs with multiple
! different references) based on the corresponding references that have
! been assigned.
logical :: tMultiRefShift = .false.

! Keep track of where in the calculation sequence we are.
integer :: calc_seq_no

! introduce a min_tau value to set a minimum of tau for the automated tau
! search
logical :: t_min_tau = .false.
real(dp) :: min_tau_global = 1.0e-7_dp

! alis suggestion: have an option after restarting to keep the time-step
! fixed to the values obtained from the POPSFILE
logical :: t_keep_tau_fixed = .false.

logical :: tPopsAlias = .false.
character(255) :: aliasStem
! new tau-search using HISTOGRAMS:
logical :: t_hist_tau_search = .false., t_hist_tau_search_option = .false.
logical :: t_fill_frequency_hists = .false.

! also use a logical, read-in in the case of a continued run, which turns
! off the tau-search independent of the input and uses the time-step
! pSingles and pDoubles values from the previous calculation.
logical :: t_previous_hist_tau = .false.

! it can be forced to do a tau-search again, if one provides an additional
! input restart-hist-tau-search in addition to the the hist-tau-search
! keyword in case the tau-search is not converged enough
logical :: t_restart_hist_tau = .false.

logical :: t_consider_par_bias = .false.

! quickly implement a control parameter to test the order of matrix element
! calculation in the transcorrelated approach
logical :: t_test_order = .false.
! for a set amount of iterations
! (in the restart case for now!)
integer :: hist_search_delay = 0

! maybe also introduce a mixing between the old and new quantities in the
! histogramming tau-search, since it is a stochastic process now
logical :: t_mix_ratios = .false.
! and choose a mixing ration p_new = (1-mix_ratio)*p_old + mix_ratio * p_new
! for now default it to 1.0_dp, meaning if it is not inputted, i only
! take the new contribution, like it is already done, and if it is
! inputted, without an additional argument default it to 0.7_dp
real(dp) :: mix_ratio = 1.0_dp

logical :: t_test_hist_tau = .false.
! real(dp) :: frq_step_size = 0.1_dp
!
!
! ! i need bin arrays for all types of possible spawns:
! integer, allocatable :: frequency_bins_singles(:), frequency_bins_para(:), &
!                         frequency_bins_anti(:), frequency_bins_doubles(:), &
!                         frequency_bins(:)
!
! ! for the rest of the tau-search, reuse the quantities from the "standard"
! ! tau search, like enough_sing, etc. although they are not global yet..
! ! so maybe define new ones to not get confused

! and i also need to truncate the spawns maybe:
logical :: t_truncate_spawns = .false.
logical :: t_truncate_unocc = .false., t_truncate_multi = .false.

logical :: t_prone_walkers, t_activate_decay
real(dp) :: n_truncate_spawns = 2.0_dp


! flags for global storage
logical :: tLogAverageSpawns, tActivateLAS
! threshold value to make something an initiator based on spawn coherence
real(dp) :: spawnSgnThresh
integer :: minInitSpawns


! introduce a cutoff for the matrix elements, to be more consistent with
! UMATEPS (let the default be zero, so no matrix elements are ignored!)
logical :: t_matele_cutoff = .false.
real(dp) :: matele_cutoff = EPS

! alis new idea to increase the chance of non-initiators to spawn to
! already occupied determinant
logical :: t_back_spawn = .false., t_back_spawn_option = .false.
! logical to control where first orbital is chosen from
logical :: t_back_spawn_occ_virt = .false.
! also use an integer to maybe start the backspawning later, or otherwise
! it may never correctly grow
integer :: back_spawn_delay = 0

! new more flexible implementation:
logical :: t_back_spawn_flex = .false., t_back_spawn_flex_option = .false.
! also make an combination of the flexible with occ-virt with an additional
! integer which manages the degree of how much you want to de-excite
! change now: we also want to enable to increase the excitation by possibly
! 1 -> maybe I should rename this than so that minus indicates de-excitation?!
integer :: occ_virt_level = 0

! make variables for automated tau determination, globally available
! 4ind-weighted variables:
real(dp) :: gamma_sing, gamma_doub, gamma_opp, gamma_par, max_death_cpt, &
            max_permitted_spawn
real(dp) :: gamma_trip
logical :: enough_sing, enough_doub, enough_opp, enough_par, consider_par_bias
logical :: enough_trip
real(dp) :: gamma_sum

real(dp) :: gamma_sing_spindiff1, gamma_doub_spindiff1, gamma_doub_spindiff2
integer :: cnt_sing, cnt_doub, cnt_opp, cnt_par, cnt_trip
! guga-specific:
integer :: cnt_four, cnt_three_same, cnt_three_mixed, cnt_two_same, cnt_two_mixed
integer :: n_opp, n_par
integer :: cnt_sing_hist, cnt_doub_hist, cnt_opp_hist, cnt_par_hist


! guga non-weighted excitation generator tau-update variables
real(dp) :: gamma_two_same, gamma_two_mixed, gamma_three_same, gamma_three_mixed, &
            gamma_four
logical :: enough_two, enough_two_same, enough_two_mixed, enough_three, &
           enough_three_same, enough_three_mixed, enough_four

! introducing an new way to adapt the time-step through H_ij/pgen frequency
! analysis: for this we need to store a histogram of the H_ij/pgens
! across all processors which are accumulated during a FCIQMC run
! the bins and boundaries need to be able to be adjusted during run-time
! to store the number of elements
integer, allocatable :: frequency_bins(:)!, all_frequency_bins(:)
! to store the boundaries of bins
real(dp), allocatable :: frequency_bounds(:)!, all_frequency_bounds(:)
logical :: t_frequency_analysis = .false. ! flag to initiate the new analysis
! change how this is approached to avoid MPI communication issues
! fix the size of the bins to 10M and the bound to 1M and the step-size to
! 0.1 and ignore all the frequency ratios above that.. and assume these
! happen really seldomly and would be cut-off anyway with the integration
! technique..
integer :: n_frequency_bins = 100000 ! optional input to adjust the number of bins
real(dp) :: max_frequency_bound = 10000.0_dp
! and also store data for the MPI communication
real(dp) :: all_max_bound = 0.0_dp
integer :: all_n_bins = 0
! use a global step-size, so no numericall error creeps in ..
real(dp) :: frq_step_size = 0.1_dp

! for automated tau-search i need more histograms to distinguish between
! the different types of excitations..
integer, allocatable :: frequency_bins_singles(:), frequency_bins_para(:), &
                        frequency_bins_anti(:), frequency_bins_doubles(:)
real(dp), allocatable :: frequency_bounds_singles(:), frequency_bounds_para(:), &
                         frequency_bounds_anti(:), frequency_bounds_doubles(:)

! for the nosym guga implementation also use different bins for the mixed
! and alike types of excitations
integer, allocatable :: frequency_bins_type2(:), frequency_bins_type2_diff(:), &
                        frequency_bins_type3(:), frequency_bins_type3_diff(:), &
                        frequency_bins_type4(:)
real(dp), allocatable :: frequency_bounds_type2(:), frequency_bounds_type2_diff(:), &
                         frequency_bounds_type3(:), frequency_bounds_type3_diff(:), &
                         frequency_bounds_type4(:)

! use also an input dependent ratio cutoff for the time-step adaptation
real(dp) :: frq_ratio_cutoff = 0.999_dp

! also use an additional flag to turn the new tau-search off but keep some
! of its functionality anyway..

! do that for the nosym guga now too deliberately! because there seems to be
! some matrix element or dynamics problem..
integer :: cnt_type2_same, cnt_type2_diff, cnt_type3_same, cnt_type3_diff, &
           cnt_type4
! and also logicals if i have enough of the excitations
logical :: enough_sing_hist, enough_doub_hist, enough_par_hist, enough_opp_hist

! keep also track of the H_ij/pgen value at the integration threshold to
! determine what we should do with the spawning events above that ..
real(dp) :: int_ratio_singles, int_ratio_para, int_ratio_anti, int_ratio_doubles


! introduce a new logical to decide if we want to calculate matrix elements
! by applying the full hamiltonian(the old way) or use the new guga matrix
! element calculation routines. to compare the influence on the time per
! iteration
logical :: t_guga_mat_eles = .true.

! introduce a flag to read the pSingles/pDoubles quantity even though the
! tau-search may be turned off
! do i want to change this to the default behavior? and indicate it
! with "no-read_probs?" to do otherwise? i think so, because why wouldn't
! i always use that, since we are always using tau-search, and even if
! we dont use it, the pSingles etc. are stored anyway, and if the are not
! stored they are 0, and in this case they are not read in anyway!
logical :: t_read_probs = .false.
! also need multiple new specific excitation type probabilites, but they are
! defined in FciMCdata module!
! move tSpinProject here to avoid circular dependencies
logical :: tSpinProject

! Use a Jacobi preconditioner in evolution equation
logical :: tPreCond

! Do we perform the death step before the communication of spawnings, or after?
! If tReplicaEstimates is is true, then it is essential that tDeathBeforeComms
! is false. To calculate these replica-based estimates, we use the summed-together
! spawnings (i.e. after communication), and the current walkers *before* death
! has been performed. So if tReplicaEstimates = .true., then we *must* have
! tDeathBeforeComms = .false.
logical :: tDeathBeforeComms

! Allow the user to input the following values for the excitation generator
real(dp) :: pSinglesIn, pParallelIn, pDoublesIn

! If true then allow set_initial_run_references to be called
logical :: tSetInitialRunRef
! If true, then when using the pops-core option, don't allow the
! pops-core-approx option to take over, even in the default case
! where it has been decided that it is efficient and appropriate.
logical :: tForceFullPops

! work on a new approximation for GUGA where I truncate based on the
! pgens or matrix elements during the excitation generation
logical :: t_trunc_guga_pgen = .false.
logical :: t_trunc_guga_matel = .false.
real(dp) :: trunc_guga_pgen = 1.0e-4_dp
real(dp) :: trunc_guga_matel = 1.0e-4_dp

! try this truncation also only for noninits
logical :: t_trunc_guga_pgen_noninits = .false.

! try a back-spawn like idea for the guga approximations
logical :: t_guga_back_spawn = .false.
logical :: t_guga_back_spawn_trunc = .false.

! this integer indicates if we want to
! -2    only treat double excitations, decreasing the excit-lvl by 2 fully
! -1    treat single and doubly excits decreasing excit-lvl by 1 or 1 fully
!  0    treat all excitations leaving the excit-lvl unchanged or lowering fully
!  1    also treat excits increasing excit-lvl up to 1 full
integer :: n_guga_back_spawn_lvl = 0


logical :: tLogGreensfunction

! Whether global_determinant_data should be moved alongside determinants
! during load balancing
logical :: tMoveGlobalDetData

! Whether we should allow non-initiators to spawn to empty dets if these
! already exist in CurrentDets.
logical :: tAllowSpawnEmpty

end module CalcData
