#include "macros.h"

MODULE Logging
            
    use constants, only: dp,int64
    use input_neci
    use MemoryManager, only: LogMemAlloc, LogMemDealloc,TagIntType
    use SystemData, only: nel, LMS, nbasis, tHistSpinDist, nI_spindist, &
                          hist_spin_dist_iter
    use constants, only: n_int, size_n_int, bits_n_int
    use bit_rep_data, only: NIfTot, NIfD
    use DetBitOps, only: EncodeBitDet
    use hist_data, only: iNoBins, tHistSpawn, BinRange
    use errors, only: Errordebug 
    use LoggingData

    IMPLICIT NONE

    contains

    subroutine SetLogDefaults()
      != Set defaults for Logging data items.

      use default_sets
      implicit none

      tDipoles = .false.
      tPrintInitiators = .false.
      tDiagAllSpaceEver = .false.
      tCalcVariationalEnergy = .false.
      tJustBlocking = .false.
      iBlockEquilShift = 0
      iBlockEquilProjE = 0
      ErrorDebug = 0
      iHighPopWrite = 15    !How many highest weighted determinants to write out at the end of an FCIQMC calc.
      tDiagWalkerSubspace = .false.
      iDiagSubspaceIter = 1
      tSplitProjEHist = .false.
      tSplitProjEHistG = .false.
      tSplitProjEHistK3 = .false.
      PopsfileTimer=0.0_dp
      tMCOutput=.true.
      tLogComplexPops=.false.
      iWriteBlockingEvery=1000
      tSaveBlocking=.false.
      OffDiagBinRange=0.001
      OffDiagMax=1.0_dp
      BinRange=0.001
      iNoBins=100000
      tHistEnergies=.false.
      tHistHamil=.false.
      iWriteHamilEvery=-1
      tHistSpawn=.false.
      iWriteHistEvery=-1
      NoACDets(:)=0
      TAutoCorr=.false.
      MaxHistE=50.0_dp
      NoHistBins=200
      iWritePopsEvery=100000
      TCalcWavevector=.false.
      WavevectorPrint=100
      TPopsFile=.true.
      tIncrementPops = .false.
      tPrintPopsDefault=.true.
      TDistrib=.false.
      ILOGGINGDef=0
      iGlobalTimerLevel=40
      nPrintTimer=10
      HFLOGLEVEL=0
      PreVarLogging=0
      TDetPops=.false.
      TZeroProjE=.false.
      TWriteDetE=.false.
      iPopsPartEvery=1
      tBinPops=.false.
      tROHistogramAll=.false.
      tROFciDump=.true.
      tTruncRODump=.false.
      tTruncDumpbyVal=.false.
      tROHistER=.false.
      tROHistDoubExc=.false.
      tROHistOffDiag=.false.
      tROHistSingExc=.false.
      tROHistOnePartOrbEn=.false.
      tROHistOneElInts=.false.
      tPrintInts=.false.
      tPrintSpinCoupHEl=.false.
      tPrintFCIMCPsi=.false.
      tCalcFCIMCPsi=.false.
      NHistEquilSteps=0
      tPrintDoubsUEG=.false.
      StartPrintDoubsUEG=0
      tPrintOrbOcc=.false.
      StartPrintOrbOcc=0
      tPrintOrbOccInit=.false.
      CCMCDebug=0
      FCIMCDebug=0
      tHFPopStartBlock=.false.
      tIterStartBlock=.false.
      IterStartBlocking=0
      HFPopStartBlocking=100
      tInitShiftBlocking=.false.
      IterShiftBlock=0
      NoDumpTruncs=0
      tWriteTransMat=.false.
      tCCMCLogTransitions=.false.
      tCCMCLogUniq=.true.
      tHistInitPops=.false.
      tHistSpinDist = .false.
      HistInitPopsIter=100000
      hist_spin_dist_iter = 1000
      tLogDets=.false.
      tCalcInstantS2 = .false.
      tCalcInstantS2Init = .false.
      tCalcInstSCpts = .false.
      instant_s2_multiplier = 1
      tRDMonFly=.false.
      tChangeVarsRDM = .false.
      RDMEnergyIter=100
      tDiagRDM=.false.
      tPrint1RDM = .false.
      tNoNOTransform = .false.
      tPrintRODump=.false.
      IterRDMonFly=0
      RDMExcitLevel=1
      tDo_Not_Calc_RDMEnergy = .false.
      tExplicitAllRDM = .false.
      tHF_S_D_Ref = .false.
      tHF_S_D = .false.
      tHF_Ref_Explicit = .false.
      twrite_normalised_RDMs = .true. 
      twrite_RDMs_to_read = .false.
      tno_RDMs_to_read = .false.
      tReadRDMs = .false.
      tNoNewRDMContrib=.false.
      IterWriteRDMs = 10000
      tWriteMultRDMs = .false.
      tInitiatorRDM = .false.
      tThreshOccRDMDiag=.false.
      ThreshOccRDM=2.0_dp
      tDumpForcesInfo = .false.
      tPrintLagrangian = .false.
      instant_s2_multiplier_init = 1
      binarypops_min_weight = 0
      tSplitPops = .false.
      tWriteCore = .false.
      tWriteTrial = .false.
      tCompareTrialAmps = .false.
      compare_amps_period = 0
      tForceCauchySchwarz = .false.
      tBrokenSymNOs = .false.
      occ_numb_diff = 0.001_dp

! Feb08 defaults
      IF(Feb08) THEN
          !Mcpaths set
          ILOGGINGDef=2
      ENDIF

    end subroutine SetLogDefaults



    subroutine LogReadInput()

        ! Read the logging section from the input file

        logical eof
        integer :: i, ierr
        character(100) :: w
        character(*), parameter :: t_r = 'LogReadInput'

      ILogging=iLoggingDef

      logging: do
        call read_line(eof)
        if (eof) then
            exit
        end if
        call readu(w)
        select case(w)

        case("REBLOCKSHIFT")
            !Abort all other calculations, and just block data again with given equilibration time (in iterations)
            tJustBlocking = .true.
            call readi(iBlockEquilShift)
        case("REBLOCKPROJE")
            !Abort all other calculations, and just block data again with given equilibration time (in iterations)
            tJustBlocking = .true.
            call readi(iBlockEquilProjE)
        case("HIGHLYPOPWRITE")
            !At the end of an FCIMC calculation, how many highly populated determinants should we write out?
            call readi(iHighPopWrite)
        case("DIAGWALKERSUBSPACE")
            !Diagonalise walker subspaces every iDiagSubspaceIter iterations
            tDiagWalkerSubspace = .true.
            call readi(iDiagSubspaceIter)
        case("DIAGALLSPACEEVER")
            !Diagonalise all space ever visited in the fciqmc dynamic. This will be written out each time HistSpawn is
            tDiagAllSpaceEver=.true.
        case("CALCVARIATIONALENERGY")
            !Calculate the variational energy of the FCIQMC dynamic each time Histspawn is calculated
            tCalcVariationalEnergy=.true.
        case("SPLITPROJE")
            !Partition contribution from doubles, and write them out
            tSplitProjEHist=.true.
        case("SPLITPROJE-G")
            !Partition contribution from doubles, and write them out; bin according to g
            tSplitProjEHistG=.true.
        case("SPLITPROJE-K3")
            !Partition contribution from doubles, and write them out; bin according to k3
            tSplitProjEHistK3=.true.
        case("NOMCOUTPUT")
            !No output to stdout from the fcimc or ccmc iterations
            tMCOutput=.false.
        case("LOGCOMPLEXWALKERS")
            !This means that the complex walker populations are now logged.
            tLogComplexPops=.true.

        case("PRINTNEWBLOCKING")
!This is the iteration interval period to write out the blocking files.
            call readi(iWriteBlockingEvery)
        case("SAVEBLOCKING")
!In this case, blocking files are not overwritten each time they are printed out, but 
            tSaveBlocking=.true.
        case("ERRORBLOCKING")
!Performs blocking analysis on the errors in the instantaneous projected energy to get the error involved.
!This is default on, but can be turned off with this keyword followed by OFF.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tHFPopStartBlock=.false.
                end select
            ELSE
                tHFPopStartBlock=.true.
            ENDIF

        case("SHIFTERRORBLOCKING")
!Performs blocking analysis on the errors in the instantaneous projected energy to get the error involved.
!This is default on, but can be turned off with this keyword followed by OFF.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tInitShiftBlocking=.false.
                end select
            ELSE
                tInitShiftBlocking=.true.
            ENDIF

        case("BLOCKINGSTARTITER")
!This keyword can be used if we want to start the blocking error analysis at a particular iteration.            
!If it is a negative integer, then this means that the blocking will start when we come out of fixed shift mode.
            tIterStartBlock=.true.
            tHFPopStartBlock=.false.
            call readi(IterStartBlocking)

        case("SHIFTBLOCKINGSTARTITER")
!This keyword can be used if we want to start the blocking error analysis of the shift at a particular 
!iteration after the shift begins to change.            
            call readi(IterShiftBlock)

        case("BLOCKINGSTARTHFPOP")            
!This keyword can be used if we want to start the blocking error analysis at a particular HF population.
!The current default is 100.
            tHFPopStartBlock=.true.
            call readi(HFPopStartBlocking)
 
        case("ROFCIDUMP")
!Turning this option on prints out a new FCIDUMP file at the end of the orbital rotation.  At the moment, the rotation is very slow
!so this will prevent us having to do the transformation every time we run a calculation on a particular system
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROFCIDUMP=.false.
                end select
            ELSE
                tROFCIDUMP=.true.
            ENDIF

        case("TRUNCROFCIDUMP")
!This options truncates the rotated FCIDUMP file by removing the specified number of virtual orbitals, based on the occupation
!numbers given by diagonalisation of the MP2 variational density matrix.
            tTruncRODump=.true.
            NoDumpTruncs=1
            ALLOCATE(NoTruncOrbs(NoDumpTruncs),stat=ierr)
            CALL LogMemAlloc('NoTruncOrbs',NoDumpTruncs,4,'Logging',NoTruncOrbsTag,ierr)
            NoTruncOrbs(:)=0
            do i=1,NoDumpTruncs
                call readi(NoTruncOrbs(i))
            enddo

        case("MULTTRUNCROFCIDUMP")
!This option allows us to specify multiple truncations, so that one calculation will print out multiple 
!ROFCIDUMP files with different
!levels of truncation - prevents us from having to do multiple identical CISD calculations to get the different truncations.
            tTruncRODump=.true.
            call readi(NoDumpTruncs)
            ALLOCATE(NoTruncOrbs(NoDumpTruncs),stat=ierr)
            CALL LogMemAlloc('NoTruncOrbs',NoDumpTruncs,4,'Logging',NoTruncOrbsTag,ierr)
            NoTruncOrbs(:)=0
            do i=1,NoDumpTruncs
                call readi(NoTruncOrbs(i))
            enddo

        case("MULTTRUNCVALROFCIDUMP")     
!This option allows us to specify particular cutoffs values for the eigenvalues - and print out multiply 
!ROFCIDUMP files with orbitals
!with eigenvalues below these removed.
            tTruncRODump=.true.
            tTruncDumpbyVal=.true.
            call readi(NoDumpTruncs)
            ALLOCATE(TruncEvalues(NoDumpTruncs),stat=ierr)
            CALL LogMemAlloc('TruncEvalues',NoDumpTruncs,8,'Logging',TruncEvaluesTag,ierr)
            TruncEvalues(:)=0.0_dp
            do i=1,NoDumpTruncs
                call readf(TruncEvalues(i))
            enddo

        case("WRITETRANSFORMMAT")
!This option writes out the transformation matrix used to convert the HF orbitals into the natural orbitals.  
!This can then be read into
!QChem to produce the natural orbital cube files and then visualise them using VMD.  Note : Currently, 
!because of Fortran 90's weird dealings
!with writing and reading binary - this option is only compatible with QChem if the code is compiled using 
!PGI - this will be fixed at 
!some stage.  Also - QChem INTDUMP files must be used to be compatible.  
            tWriteTransMat=.true.


        case("HIST-SPIN-DIST")
            ! Histogram the distribution of walkers within determinants of the
            ! given spatial configuration
            ! --> The determinant is specified using SPIN orbitals, but these
            !     are converted to a spatial structure for use.

            tHistSpinDist = .true.
            call readi(hist_spin_dist_iter)
            if (.not. allocated(nI_spindist)) &
                allocate(nI_spindist(nel))
            nI_spindist = 0
            i = 1
            do while (item < nitems .and. i <= nel)
                call geti(nI_spindist(i))
                i = i+1
            enddo


        case("ROHISTOGRAMALL")
!This option goes with the orbital rotation routine.  If this keyword is included, all possible histograms are included.
!These maybe also turned off/on with individual keywords.
!As it stands, the bins run from -1 to 1 with increments of 0.05. These parameters may be made options in the future.
            tROHistogramAll=.true.
            tROHistOffDiag=.true.
            tROHistDoubExc=.true.
            tROHistSingExc=.true.
            tROHistOneElInts=.true.
            tROHistOnePartOrbEn=.true.
            tROHistER=.true.
            tROHistVirtCoulomb=.true.

       case("ROHISTOFFDIAG")
!This option creates a histogram of the <ij|kl> terms where i<k and j<l.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistOffDiag=.false.
                end select
            ELSE
                tROHistOffDiag=.true.
            ENDIF

       case("ROHISTDOUBEXC")
!This option creates a histogram of the 2<ij|kl>-<ij|lk> terms, the off diagonal hamiltonian elements for double excitations.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistDoubExc=.false.
                end select
            ELSE
                tROHistDoubExc=.true.
            ENDIF
 
       case("ROHISTSINGEXC")
!This option creates a histogram of the single excitation hamiltonian elements.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistSingExc=.false.
                end select
            ELSE
                tROHistSingExc=.true.
            ENDIF

       case("ROHISTER")
!This option creates a histogram of the <ii|ii> terms, the ones that are maximised in the edmiston-reudenberg localisation.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistER=.false.
                end select
            ELSE
                tROHistER=.true.
            ENDIF
  
       case("ROHISTONEElINTS")
!This option creates a histogram of the one electron integrals, the <i|h|i> terms.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistOneElInts=.false.
                end select
            ELSE
                tROHistOneElInts=.true.
            ENDIF
         

       case("ROHISTONEPARTORBEN")
!This option creates a histogram of the one particle orbital energies, epsilon_i = <i|h|i> + sum_j [<ij||ij>].
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistOnePartOrbEn=.false.
                end select
            ELSE
                tROHistOnePartOrbEn=.true.
            ENDIF
 
       case("ROHISTVIRTCOULOMB")
!This option creates a histogram of the coulomb integrals, <ij|ij>, where i and j are both virtual and i<j.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tROHistVirtCoulomb=.false.
                end select
            ELSE
                tROHistVirtCoulomb=.true.
            ENDIF
        
        case("PRINTINTEGRALS")
!This option prints 2 files containing the values of certain integrals at each rotation iteration.  This is so that we can see the
!effect the rotation is having on all values, other than just the one we are max/minimising.
            tPrintInts=.true.

        case("PRINTTRICONNECTIONS")
!This option takes each generated pair of determinant and excitation and finds 3rd determinant to make up a triangular connection.
!The product of the three connecting elements are then histogrammed in two separate files. 
!In one, the triangular connections that combine 
!to be sign coherent are recorded, and in the other, those which are sign incoherent.
            CALL Stop_All(t_r,"PRINTTRICONNECTIONS option depreciated")
!            call readf(TriConMax)
!            call readi(NoTriConBins)
!            tPrintTriConnections=.true.

        case("HISTTRICONNELEMENTS")
!This keyword takes the above triangles of connected determinants and histograms each connecting 
!element that contributes to the triangle.
!It then prints these according to whether they are single or double connecting elements.
!It also prints a histogram and the average size of the Hjk elements (regardless of whether or not they are zero).
            CALL Stop_All(t_r,"HISTTRICONNELEMENTS option depreciated")
!            call readf(TriConHElSingMax)
!            call readf(TriConHElDoubMax)
!            call readi(NoTriConHElBins)
!            tHistTriConHEls=.true.

        case("PRINTHELACCEPTSTATS")
!This keyword prints out an extra file that keeps track of the H elements involved in spawning attempts 
!that are accepted or not accepted.
!It prints out the average H elements where spawning is accepted and the average where it is not accepted.
            CALL Stop_All(t_r,"PRINTHELACCEPTSTATS option depreciated")
!            tPrintHElAccept=.true.

        case("PRINTSPINCOUPHELS")
!This option prints out the number of positive and negative (and their sums) H elements connecting two spin coupled determinants.
            tPrintSpinCoupHEl=.true.

        case("HISTINITIATORPOPS")
!This option prints out a file (at every HistInitPopsIter iteration) containing the 
!natural log of the populations of the initiator determinants 
!and the number with this population. The range of populations histogrammed goes from ln(N_add) -> ln(1,000,000) with 50,000 bins.
            tHistInitPops=.true.
            call readi(HistInitPopsIter)

        case("CALCRDMONFLY")
!This keyword sets the calculation to calculate the reduced density matrix on the fly.  
!This starts at IterRDMonFly iterations after the shift changes.
!If RDMExcitLevel = 1, only the 1 electron RDM is found, if RDMExcitLevel = 2, 
! only the 2 electron RDM is found and if RDMExcitLevel = 3, both are found. 
            tRDMonFly=.true.
            call readi(RDMExcitLevel)
            call readi(IterRDMonFly)
            call readi(RDMEnergyIter)

        case("DIAGFLYONERDM")
!This sets the calculation to diagonalise the *1* electron reduced density matrix.   
!The eigenvalues give the occupation numbers of the natural orbitals (eigenfunctions).
            tDiagRDM=.true.

        case("NONOTRANSFORM")
! This tells the calc that we don't want to print the NO_TRANSFORM matrix.            
! i.e. the diagonalisation is just done to get the correlation entropy.
            tNoNOTransform = .true.

        case("PRINTONERDM")
! This prints the OneRDM, regardless of whether or not we are calculating just the 1-RDM, or the 2-RDM.            
            tPrint1RDM = .true.

        case("PRINTRODUMP")
            tPrintRODump=.true.
            tROFciDump = .true.
! This is to do with the calculation of the MP2 or CI natural orbitals.  
!This should be used if we want the transformation matrix of the 
! natural orbitals to be found, but no ROFCIDUMP file to be printed (i.e. 
!the integrals don't need to be transformed).  This is so that at the end 
! of a calculation, we may get the one body reduced density matrix from the 
!wavefunction we've found, and then use the MOTRANSFORM file printed to 
! visualise the natural orbitals with large occupation numbers.

        case("FORCECAUCHYSCHWARZ")
            tForceCauchySchwarz=.true.
!This forces the inequality gamma_pq <= sqrt(gamma_pp * gamma_qq) is obeyed.
!we choose min(gamma_pq, sqrt(gamma_pp * gamma_qq) to ensure a positive-definite matrix
!May be of use for getting orbitals from an approximate initial FCIQMC calc.

        case("BROKENSYMNOS")
            tBrokenSymNOs = .true.
            call readf(occ_numb_diff)
! This is to rotate the obtained natural orbitals (NOs) again in order to obtain
! symmetry broken NOs: pairs of NOs whose occupation numbers differ by less 
! than the specified threshold occ_numb_diff (relative difference, i.e. difference
! divided by absolute value) will be rotated so as to 
! maximally localise them using and Edminston Ruedenberg type localisation
! A new FCIDUMP file (BSFCIDUMP) with the rotated NOs is printed out

        case("DIPOLE_MOMENTS")
            !Calculate the dipole moments if we are in molpro
            tDipoles = .true.

        case("CALCRDMENERGY")
!This takes the 1 and 2 electron RDM and calculates the energy using the RDM expression.            
!For this to be calculated, RDMExcitLevel must be = 3, so there is a check to make sure this 
!is so if the CALCRDMENERGY keyword is present.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tDo_Not_Calc_RDMEnergy=.true.
                end select
            ELSE
                tDo_Not_Calc_RDMEnergy=.false.
            ENDIF

        case("EXPLICITALLRDM")
!Explicitly calculates all the elements of the RDM.            
            tExplicitAllRDM = .true.

        case("HFREFRDMEXPLICIT")
!Uses the HF as a reference and explicitly calculates the RDM to find the energy - should be same as projected energy, 
!when printing out every shift update.
            tHF_Ref_Explicit = .true.

        case("HFSDRDM")
!Calculate the RDM for the HF, singles and doubles only - symmetrically.            
            tHF_S_D = .true.

        case("HFSDREFRDM")
!Uses the HF, singles and doubles as a multiconfigurational reference and calculates the RDM to find the energy.            
            tHF_S_D_Ref = .true.

        case("WRITEINITIATORS")
! Requires a popsfile to be written out.  Writes out the initiator populations. 
            tPrintInitiators = .true.
        
        case("WRITERDMSTOREAD")
! Writes out the unnormalised RDMs (in binary), so they can be read back in, and the calculations restarted at a later point 
! This is also tied to the POPSFILE/BINARYPOPS keyword - so if we're writing a normal POPSFILE, we'll write this too, 
! unless **WRITERDMSTOREAD** OFF is used.
            IF(item.lt.nitems) THEN
                call readu(w)
                select case(w)
                    case("OFF")
                        tno_RDMs_to_read = .true. 
                        twrite_RDMs_to_read = .false. 
                end select
            ELSE
                twrite_RDMs_to_read = .true. 
                tno_RDMs_to_read = .false. 
            ENDIF

        case("NONORMRDMS")            
! Does not print out the normalised (final) RDMs - to be used if you know the calculation will not be converged, and don't  
! want to take up disk space.
            twrite_normalised_RDMs = .false.

        case("READRDMS")
! Read in the RDMs from a previous calculation, and continue accumulating the RDMs from the very beginning of this restart. 
            tReadRDMs = .true.

        case("NONEWRDMCONTRIB")
            ! To be used with READRDMs.  This option makes sure that we don't add in any new contributions to the RDM if filling stochastically
            !This is useful if we want to read in an RDM from another calculation and then just print out the analysis, without adding in any more information.
      tNoNewRDMContrib=.true.

        case("WRITERDMSEVERY")
! Write out the normalised, hermitian RDMs every IterWriteRDMs iterations.  
            tWriteMultRDMs = .true.
            call readi(IterWriteRDMs)

        case("INITIATORRDM")
! Use only the determinants that are (on average) initiators to calculate the RDMs.
            tInitiatorRDM = .true.
        
        case("THRESHOCCONLYRDMDIAG")
            !Only add in a contribution to the diagonal elements of the RDM if the average sign of the determinant is greater than [ThreshOccRDM]
            tThreshOccRDMDiag=.true.
            call Getf(ThreshOccRDM)

        case("DUMPFORCESINFO")
! Using the finalised 2RDM, calculate the Lagrangian X used for the calculation of the forces, 
!and dump all these in Molpro-friendly format
            tDumpForcesInfo = .true.
        
        case("PRINTLAGRANGIAN")
            ! Print out the Lagrangian X to file (Only works in conjuction with DUMPFORCESINFO: otherwise, this option does nothing)
            tPrintLagrangian = .true.
        
        case("AUTOCORR")
!This is a Parallel FCIMC option - it will calculate the largest weight MP1 determinants and histogramm them
!HF Determinant is always histogrammed. NoACDets(2) is number of doubles. NoACDets(3) is number of triples and NoACDets(4) is 
!number of quads to histogram.
            TAutoCorr=.true.
            CALL Stop_All("LogReadInput","The ACF code has been commented out in the FCIMCPar module")
            do i=2,4
                IF(item.lt.nitems) call readi(NoACDets(i))
            enddo
        case("DETPOPS")
!This option no longer works...
            TDetPops=.true.
        case("DISTRIBS")
            TDistrib=.true.
        case("HISTPARTENERGIES")
!This will histogram the hamiltonian matrix elements of the particles in the parallel FCIMC algorithm.
            tHistEnergies=.true.
            call readf(BinRange)
            call readi(iNoBins)
            call readf(OffDiagBinRange)
            call readf(OffDiagMax)
            IF(OffDiagMax.lt.0.0_dp) THEN
                OffDiagMax=-OffDiagMax
            ENDIF
        case("HISTSPAWN")
!This option will histogram the spawned wavevector, averaged over all previous iterations. 
!It scales horrifically and can only be done for small systems
!which can be diagonalized. It requires a diagonalization initially to work. 
!It can write out the average wavevector every iWriteHistEvery.
            tHistSpawn=.true.
            IF(item.lt.nitems) call readi(iWriteHistEvery)
        case("HISTHAMIL")
!This option will histogram the spawned hamiltonian, averaged over all previous iterations. It scales horrifically 
!and can only be done for small systems
!which can be diagonalized. It will write out the hamiltonian every iWriteHamilEvery.
            tHistHamil=.true.
            IF(item.lt.nitems) call readi(iWriteHamilEvery)
        case("BLOCKEVERYITER")
!This will block the projected energy every iteration with the aim of achieving accurate error estimates. 
!However, this does require a small amount of additional communication.
            tBlockEveryIteration=.true.
        case("PRINTFCIMCPSI")
            tPrintFCIMCPsi=.true.
            tCalcFCIMCPsi=.true.
        case("HISTEQUILSTEPS")
!This option sets the histogramming to only be done after the specified number of iterations.            
            call readi(NHistEquilSteps)
        case("PRINTORBOCCS")
!This option initiates the above histogramming of determinant populations and then at the end of the 
!spawning uses these to find the normalised  
!contribution of each orbital to the total wavefunction.  
            tPrintOrbOcc=.true.
            IF(item.lt.nitems) call readi(StartPrintOrbOcc)
        case("PRINTDOUBSUEG")
!This option initiates the above histogramming of doubles for the UEG
!            if (.not.tUEG) call stop_all("Logging","Printdoubs doesn't work with systems other than UEG")
            tPrintDoubsUEG=.true.
            IF(item.lt.nitems) call readi(StartPrintDoubsUEG)
        case("PRINTORBOCCSINIT")
!This option initiates the above histogramming of determinant populations and then 
!at the end of the spawning uses these to find the normalised  
!contribution of each orbital to the total wavefunction.  
            tPrintOrbOcc=.true.
            tPrintOrbOccInit=.true.
            IF(item.lt.nitems) call readi(StartPrintOrbOcc)
        case("POPSFILE")
! This is so that the determinants at the end of the MC run are written
! out, to enable them to be read back in using READPOPS in the Calc section,
! if you want to restart the simulation at a later date.  !iWritePopsEvery
! will write the configuration of particles out each time the iteration
! passes that many.
            TPopsFile=.true.
            IF(item.lt.nitems) THEN
                call readi(iWritePopsEvery)
                IF(iWritePopsEvery.lt.0) THEN
!If a negative argument is supplied to iWritePopsEvery, then the POPSFILE will 
!never be written out, even at the end of a simulation.
!If it is exactly zero, this will be the same as without any argument, and a
!popsfile will only be written out in the instance of a clean exit
                    TPopsFile=.false.
                    tPrintPopsDefault=.false.
                ELSEIF(iWritePopsEvery.gt.0) THEN
                    tPrintPopsDefault=.false.
                ENDIF
            ENDIF
        case("REDUCEDPOPSFILE")
!A reduced popsfile works in exactly the same way as a normal popsfile, but only every iPopsPartEvery particle is printed out.
            TPopsFile=.true.
            call readi(iWritePopsEvery)
            call readi(iPopsPartEvery)
        case("POPSFILETIMER")
            call readf(PopsfileTimer)   !Write out a POPSFILE every "PopsfileTimer" hours.

        case("BINARYPOPS")
            ! This means that the popsfile (full or reduced) will now be 
            ! written out in binary format. This should now take up less 
            ! space, and be written quicker.
            !
            ! By default, all particles are written into the popsfile. If 
            ! a minimum weight is proveded, only those particles with at ]east
            ! that weight are included.
            tBinPops= .true.
            if (item < nitems) then
                call readf(binarypops_min_weight)
            end if

        case("INCREMENTPOPS")
! Don't overwrite existing POPSFILES.
            tIncrementPops = .true.
        case("CCMCDEBUG")
!CCMC debugging level. Takes an integer 0-6
            call readi(CCMCDebug)
        case("FCIMCDEBUG")
!FCIQMC debugging level. Takes an integer 0-6
            call readi(FCIMCDebug)
        case("ERRORDEBUG")
!Error analysus debugging level. Takes an integer 0-6
            call readi(ErrorDebug)
        case("CCMCLOGTRANSITIONS")
            tCCMCLogTransitions=.true.
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
               case("NONUNIQUE")
                  tCCMCLogUniq=.false.
               case("UNIQUE")
                  tCCMCLogUniq=.true.
               case default
                  CALL report("Logging keyword CCMCLOGTRANSITIONS "//trim(w)       &
     &               //" not recognised",.true.)
               end select
            enddo
        case("WRITEDETE")
!This logging option will write out the energies of all determinants which have been spawned at in the simulation
! The two input options are the number of bins, and the maximum determinant energy to be histogrammed.
            TWriteDetE=.true.
            IF(item.lt.nitems) call readi(NoHistBins)
            IF(item.lt.nitems) call readf(MaxHistE)
        case("ZEROPROJE")
! This is for FCIMC when reading in from a POPSFILE. If this is on, then the energy 
! estimator will be restarted.
            TZeroProjE=.true.
        case("WAVEVECTORPRINT")
! This is for FCIMC - if on, it will calculate the exact eigenvector and
! values initially, and then print out the running wavevector every
! WavevectorPrint MC steps. However, this is slower.
            TCalcWavevector=.true.
            call readi(WavevectorPrint)
        case("MCPATHS")
            ILOGGING = IOR(ILOGGING,2**1)
        case("BLOCKING")
            ILOGGING = IOR(ILOGGING,2**13)
        case("PREVAR")
            ILOGGING = IOR(ILOGGING,2**14)
        case("FMCPR")
!  We log the value
            ILOGGING = IOR(ILOGGING,2**0)
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
               case("LABEL")
                   ILOGGING = IOR(ILOGGING,2**2)
               case("RHO")
                   ILOGGING = IOR(ILOGGING,2**3)
               case("1000")
                   ILOGGING = IOR(ILOGGING,2**9)
               case("EXCITATION")
                   ILOGGING = IOR(ILOGGING,2**12)
               case("XIJ")
                   ILOGGING = IOR(ILOGGING,2**6)
               case("")
                   ILOGGING = IOR(ILOGGING,2**2)
               case default
                  CALL report("Logging keyword FMCPR "//trim(w)       &
     &               //" not recognised",.true.)
               end select
            enddo
        case("CALCPATH")
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
               case("LABEL")
                   ILOGGING = IOR(ILOGGING,2**4)
               case("RHO")
                   ILOGGING = IOR(ILOGGING,2**5)
               case("")
                   ILOGGING = IOR(ILOGGING,2**4)
               case default
                  CALL report("Logging keyword CALCPATH "//trim(w)    &
     &               //" not recognised",.true.)
               end select
            enddo
        case("XIJ")
            ILOGGING = IOR(ILOGGING,2**6)
        case("HAMILTONIAN")
            ILOGGING = IOR(ILOGGING,2**7)
        case("PSI")
            ILOGGING = IOR(ILOGGING,2**8)
        case("TIMING")
            do while(item.lt.nitems)
                call readu(w)
                select case(w)
                case("LEVEL")
                    call readi(iGlobalTimerLevel)
                case("PRINT")
                    call readi(nPrintTimer)
                case default
                    call reread(-1)
                    call readi(iGlobalTimerLevel)
                end select
            end do
        case("VERTEX")
            do while(item.lt.nitems)
               call readu(w)
               select case(w)
              ! case("1000")
              !     ILOGGING = IOR(ILOGGING,2**9)
               case("EVERY")
                   ILOGGING = IOR(ILOGGING,2**10)
               case default
                  call reread(-1)
                  call geti(G_VMC_LOGCOUNT)
                  ILOGGING = IOR(ILOGGING,2**9)
!                  CALL report("Logging keyword VERTEX "//trim(w)    &
!     &               //" not recognised",.true.)
               end select
            end do
        case("HFBASIS")
            ILOGGING = IOR(ILOGGING,2**11)
        case("HFLOGLEVEL")
            call geti(HFLOGLEVEL)
        case("SAVEPREVARLOGGING")
             PreVarLogging=iLogging
             iLogging=iLoggingDef
        case("DETS")
            tLogDets=.true.
        case("DETERMINANTS")
            tLogDets=.true.

        case ("INSTANT-S2-FULL")
            ! Calculate an instantaneous value for S^2, and output it to the
            ! relevant column in the FCIMCStats file.
            !
            ! The second parameter is a multiplier such that we only calculate
            ! S^2 once for every n update cycles (it must be on an update
            ! cycle such that norm_psi_squared is correct)
            tCalcInstantS2 = .true.
            if (item < nitems) &
                call readi (instant_s2_multiplier)

        case ("INSTANT-S2-INIT")
            ! Calculate an instantaneous value ofr S^2 considering only the
            ! initiators, and output it to the relevant column in the
            ! FCIMCStats file.
            !
            ! The second parameter is a multiplier such that we only calculate
            ! S^2 once for every n update cycles (it must be an update
            ! cycle such that norm_psi_squared is correct).
            tCalcInstantS2Init = .true.
            if (item < nitems) &
                call readi (instant_s2_multiplier_init)

        case ("INSTANT-S-CPTS")
            ! Calculate components of the wavefunction with each value of S.
            ! n.b. This is NOT quantitatively correct.
            !      --> Only of QUALITATIVE utility.
            tCalcInstSCpts = .true.

        case ("SPLIT-POPS")
            ! Do we want to split a popfile up into multiple parts which are
            ! output on each of the nodes, and need to be combined/split-up and
            ! distributed to the nodes on our (sequential) time rather than
            ! on multi-processor time?
            tSplitPops = .true.
            tBinPops = .true.

        case("WRITE-CORE")
            ! Output the semi-stochastic core space to a file.
            tWriteCore = .true.

        case("WRITE-TRIAL")
            ! Output the trial wavefunction space to a file.
            tWriteTrial = .true.

        case("COMPARE-TRIAL-AND-FCIQMC-AMPS")
            tCompareTrialAmps = .true.
            call readi(compare_amps_period)

        case("ENDLOG")
            exit logging
        case default
           CALL report("Logging keyword "//trim(w)//" not recognised",.true.)
        end select
      end do logging
    END SUBROUTINE LogReadInput


END MODULE Logging
