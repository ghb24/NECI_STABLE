#include "macros.h"

module LoggingData

    use constants, only: dp,int64
    use MemoryManager, only: TagIntType

    implicit none

    save

    INTEGER :: iDiagSubspaceIter
    LOGICAL :: tDiagWalkerSubspace
    INTEGER ILOGGING,ILOGGINGDef
    INTEGER :: iGlobalTimerLevel=40
    INTEGER nPrintTimer,G_VMC_LOGCOUNT
    INTEGER HFLOGLEVEL,iWritePopsEvery,StartPrintOrbOcc
    INTEGER PreVarLogging,WavevectorPrint,NoHistBins,HistInitPopsIter
    real(dp) MaxHistE,OffDiagMax,OffDiagBinRange,PopsfileTimer
    LOGICAL TDistrib,TPopsFile,TCalcWavevector,TDetPops,tROFciDump,tROHistOffDiag,tROHistDoubExc,tROHistOnePartOrbEn
    LOGICAL tPrintPopsDefault
    LOGICAL TZeroProjE,TWriteDetE,TAutoCorr,tBinPops,tIncrementPops,tROHistogramAll,tROHistER,tROHistSingExc
    LOGICAL tRoHistOneElInts
    LOGICAL tROHistVirtCoulomb,tPrintInts,tHistEnergies,tTruncRODump,tRDMonFly,tDiagRDM,tDo_Not_Calc_2RDM_est
    LOGICAL tPrintFCIMCPsi,tCalcFCIMCPsi,tPrintSpinCoupHEl,tIterStartBlock,tHFPopStartBlock,tInitShiftBlocking
    LOGICAL tTruncDumpbyVal, tPrintRODump, tno_RDMs_to_read, tReadRDMs, tNoNewRDMContrib
    LOGICAL tWriteTransMat,tPrintOrbOcc,tHistInitPops,tPrintOrbOccInit, tWriteMultRDMs
    LOGICAL twrite_normalised_RDMs, tWriteSpinFreeRDM, twrite_RDMs_to_read
    LOGICAL tNoNOTransform, tPrint1RDM, tPrintInitiators
    INTEGER NoACDets(2:4),iPopsPartEvery,iWriteHistEvery,NHistEquilSteps
    INTEGER IterRDMonFly, RDMExcitLevel, RDMEnergyIter, IterWriteRDMs
    INTEGER FCIMCDebug !FciMC Debugging Level 0-6.  Default 0

    ! Logical(4) datatypes for compilation with builds of openmpi that don't
    ! have support for logical(8). Gah.
    logical :: tExplicitAllRDM, tChangeVarsRDM
    logical :: tPopAutoAdaptiveShift
    LOGICAL tSaveBlocking !Do not overwrite blocking files
    INTEGER iWriteBlockingEvery !How often to write out blocking files
    INTEGER IterStartBlocking,HFPopStartBlocking,NoDumpTruncs
    INTEGER(TagIntType)  OrbOccsTag,HistInitPopsTag,AllHistInitPopsTag,NoTruncOrbsTag,TruncEvaluesTag
    INTEGER , ALLOCATABLE :: NoTruncOrbs(:),HistInitPops(:,:),AllHistInitPops(:,:)
    real(dp) , ALLOCATABLE :: TruncEvalues(:),OrbOccs(:),DoubsUEG(:,:,:,:),DoubsUEGLookup(:)
    LOGICAL, ALLOCATABLE :: DoubsUEGStore(:,:,:)
    LOGICAL :: tBlockEveryIteration
    LOGICAL tLogDets       ! Write out the DETS and SymDETS files.
    LOGICAL tLogComplexPops     ! Write out complex walker information
    LOGICAL tLogEXLEVELStats    ! Write L_{0,1,2} norms of weights by exlevel
    LOGICAL tMCOutput
    logical :: tDumpForcesInfo
    logical :: tPrintLagrangian  !Print out the 1RDM,2RDM and Lagrangian to file at the end of a run as long as 2RDM is calculated
    real(dp) :: ThreshOccRDM
    logical :: tFullHFAv, tThreshOccRDMDiag
    logical :: tRDMInstEnergy

    logical :: tCalcInstantS2, tCalcInstSCpts, tCalcInstantS2Init
    integer :: instant_s2_multiplier, instant_s2_multiplier_init
    integer :: iHighPopWrite

    !Just do a blocking analysis on previous data
    logical :: tJustBlocking
    integer :: iBlockEquilShift,iBlockEquilProjE
    logical :: tDiagAllSpaceEver,tCalcVariationalEnergy

    ! Do we want to split the popsfile up into multiple bits?
    logical :: tSplitPops

    ! What is the mininum weight on a determinant required for it to be
    ! included in a binary pops file?
    real(dp) :: binarypops_min_weight

    ! If true, output the core/trial spaces to a file.
    logical :: tWriteCore
    logical :: tWriteTrial

    ! If true then, at the end of a calculation, find the write_end_core_size
    ! most populated determinants and write them to a CORESPACE file.
    logical :: tWriteCoreEnd
    integer :: write_end_core_size

    ! If true, output to a file the FCIQMC amplitudes in the trial space against the amplitudes of the trial wavefunction.
    logical :: tCompareTrialAmps
    ! Output the above data to a file every compare_amps_period iterations.
    integer :: compare_amps_period
    logical :: tDipoles !Do we want to calculate the dipole moments
    logical :: tHistExcitToFrom

    logical :: log_cont_time_survivals, tNoWarnIC0Bloom, tDumpHamilBinary, &
               tDumpHamilOverlap
    logical :: tFCIMCStats2

    !If we want to force the Cauchy--Schwarz inequality (e.g. if we know the 1RDM is undersampled)
    logical :: tForceCauchySchwarz
    ! If we'd like to rotate the NOs again so as to obtain broken symmetry NOs
    logical :: tBrokenSymNOs,tBreakSymNOs
    real(dp) :: occ_numb_diff
    integer :: rottwo,rotthree,rotfour,local_cutoff
    integer, allocatable :: RotNOs(:)
    integer(TagIntType) :: tagRotNOs

    ! If not true then don't output data tables to FCIMCStats, INITIATORStats or standard output.
    logical :: tPrintDataTables

    ! Should we output the load-balanced distribution?
    logical :: tOutputLoadDistribution

    logical :: tHDF5PopsRead, tHDF5PopsWrite

    ! Whether to write another HDF5 popsfile with dets restricted to a maximum
    ! exitation level
    logical :: tHDF5TruncPopsWrite
    ! The maximum excitation level of dets writen to truncated HDF5 popsfile
    integer :: iHDF5TruncPopsEx

    ! Whether to calculate and print the instanenous project energy of
    ! wavefunction printed to popsfile
    logical :: tPopsInstProjE

    ! Whether to accumulate the population of determinants and write them
    ! to the popsfile
    logical :: tAccumPops
    logical :: tAccumPopsActive

    ! When to start accumlating the population
    integer :: iAccumPopsIter

    ! Number of iterations the empty dets are kept during accumlation
    integer :: iAccumPopsExpire

    ! Maximum excitation level of empty dets to keep during accumlation
    integer :: iAccumPopsMaxEx

    logical :: tOldRDMs = .false.

    logical :: tTransitionRDMs = .false.

    ! If true, then read in 2-RDM popsfiles and then output 1-RDMs calculated
    ! directly from these.
    logical :: tPrint1RDMsFrom2RDMPops = .false.
    ! If true, then read in spinfree 2-RDM files and then output 1-RDMs
    ! calculated directly from these.
    logical :: tPrint1RDMsFromSpinfree = .false.

    ! Whether we are calculating the estimates of properties
    logical ::  tCalcPropEst
    integer :: iNumPropToEst
    ! The name of the integral file for each of the property to be estimated
    character(100), allocatable :: EstPropFile(:)

    ! like rdms, as it is a bit similar, access the double occupancy
    ! measurement in the logging section!
    logical :: t_calc_double_occ = .false.
    ! also use a optional input parameter to start averaging the
    ! double occupancy only after a certain number of steps after the
    ! shift changes
    integer :: equi_iter_double_occ = 0
    logical :: t_calc_double_occ_av = .false.
    ! I essentially only need a local and a global storage for the
    ! the expectation vaulue <n_u n_d>
    ! and also some storage for the instantaneous, averaged, summed over
    ! stuff etc..
!     real(dp) :: n_double_occ_loc, n_double_occ_all
    ! [Werner Dobrautz 4.4.2017]
    ! changes belonging to the histogram tau-search
    ! for now always print out the histograms at the end, maybe change that
    ! behavior in the future
    logical :: t_print_frq_histograms = .true.

    ! up to which exLvl we track the number of initiators per exLvl
    integer :: maxInitExLvlWrite
    integer, allocatable :: initsPerExLvl(:), allInitsPerExLvl(:)

    ! if this is true, force moving fcimcstats and initiatorstats files, and accumulate stats in new files
    logical :: t_no_append_stats = .false.

    ! for giovanni introduce a new keyword to print out the spin-resolved
    ! 2-rdms, especially for systems with ms /= 0
    logical :: t_spin_resolved_rdms = .false.

    ! Alis idea to histogram bad p(a|ij) with the triplet (aij)
    logical :: t_log_ija = .false.
    ! i need bins of (nbasis,nbasis,nbasis) size (although i guess spatial
    ! orbitals would be enough)
    ! store more information for these.. but how do i MPI communicate that?
    ! not now!
    ! do it seperately for singles and parallel and opposite spin excitations
    ! and also store the number of symmetry allowed orbitals! just to check
    integer, allocatable :: ija_bins_sing(:), all_ija_bins_sing(:), &
                            ija_bins_para(:,:,:), all_ija_bins(:,:,:), &
                            ija_bins_anti(:,:,:), &
                            ija_orbs_sing(:), all_ija_orbs_sing(:), &
                            ija_orbs_para(:,:,:), all_ija_orbs(:,:,:), &
                            ija_orbs_anti(:,:,:)


    ! also use a use defined threshold
    real(dp) :: ija_thresh = 1.0e-7_dp

    ! for the output of the references in the adi-mode
    logical :: tWriteRefs
    character(255) :: ref_filename
    ! for the histogramming of the acceptance rates used in the adaptive shift mode
    logical :: t_hist_fvals
    integer :: enGrid, arGrid

end module LoggingData
