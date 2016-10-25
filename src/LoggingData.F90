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
    LOGICAL tROHistVirtCoulomb,tPrintInts,tHistEnergies,tTruncRODump,tRDMonFly,tDiagRDM,tDo_Not_Calc_RDMEnergy
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

    logical :: t_print_frq_histograms = .false.
end module LoggingData
