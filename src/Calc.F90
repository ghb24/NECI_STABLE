#include "macros.h"

MODULE Calc
        
    use CalcData
    use SystemData, only: beta, nel, STOT, tCSF, LMS, tSpn, AA_elec_pairs, &
                          BB_elec_pairs, par_elec_pairs, AB_elec_pairs, &
                          AA_hole_pairs, BB_hole_pairs, AB_hole_pairs, &
                          par_hole_pairs, hole_pairs, nholes_a, nholes_b, &
                          nholes, UMATEPS, tHub
    use Determinants, only: write_det
    use spin_project, only: spin_proj_interval, tSpinProject, &
                            spin_proj_gamma, spin_proj_shift, &
                            spin_proj_cutoff, spin_proj_stochastic_yama, &
                            spin_proj_spawn_initiators, spin_proj_no_death, &
                            spin_proj_iter_count, spin_proj_nopen_max, &
                            disable_spin_proj_varyshift
    use default_sets
    use Determinants, only: iActiveBasis, SpecDet, tSpecDet, nActiveSpace, &
                            tDefineDet
    use DetCalc, only: iObs, jObs, kObs, DETINV, &
                       icilevel, tBlock, tCalcHMat, tEnergy, tRead, &
                       tFindDets
    use DetCalcData, only: B2L, nKry, nEval, nBlk, nCycle
    use IntegralsData, only: tNeedsVirts
    use FciMCData, only: tTimeExit,MaxTimeExit, InputDiagSft, tSearchTau, &
                         nWalkerHashes, HashLengthFrac, tSearchTauDeath, &
                         tTrialHash, tIncCancelledInitEnergy, MaxTau, &
                         tStartCoreGroundState, pParallel, pops_pert, &
                         alloc_popsfile_dets, tSearchTauOption, &
                         sFAlpha, tEScaleWalkers, sFBeta, sFTag, tLogNumSpawns
    use adi_data, only: maxNRefs, nRefs, tAllDoubsInitiators, tDelayGetRefs, &
         tDelayAllDoubsInits, tAllSingsInitiators, tDelayAllSingsInits, tSetDelayAllDoubsInits, &
         tSetDelayAllSingsInits, nExProd, NoTypeN, tAdiActive, tReadRefs, SIUpdateInterval, &
         tProductReferences, tAccessibleDoubles, tAccessibleSingles, &
         tReferenceChanged, superInitiatorLevel, allDoubsInitsDelay, tStrictCoherentDoubles, &
         tWeakCoherentDoubles, tAvCoherentDoubles, coherenceThreshold, SIThreshold, &
         tSuppressSIOutput, targetRefPop, targetRefPopTol, tSingleSteps, tVariableNRef, &
         nRefsSings, nRefsDoubs, minSIConnect, tWeightedConnections, tSignedRepAv
    use ras_data, only: core_ras, trial_ras
    use load_balance, only: tLoadBalanceBlocks
    use ftlm_neci
    use spectral_data
    use spectral_lanczos, only: n_lanc_vecs_sl
    use exact_spectrum
    use perturbations, only: init_perturbation_creation, init_perturbation_annihilation

    implicit none

contains

    subroutine SetCalcDefaults()
          use FciMCData, only: hash_shift
        
        ! Set defaults for Calc data items.

        ! Values for old parameters.
        ! These have no input options to change the defaults, but are used in
        ! the code.
          InputTargetGrowRateWalk = 500000
          InputTargetGrowRate = 0.0_dp
          InitialPart = 1.0_dp
          B2L = 1.0e-13_dp
          TMC = .false.
          NHISTBOXES = 0
          TREADRHO = .false.
          TRHOIJ = .false.
          TBEGRAPH = .false.

          if (Nov11) then
              tInstGrowthRate=.true.
          else
              tInstGrowthRate = .false.
          end if

!       Calc defaults 
          iSampleRDMIters = -1
          tStartCoreGroundState = .true.
          HashLengthFrac = 0.7_dp
          nWalkerHashes=0
          tTrialHash=.true.
          tIncCancelledInitEnergy = .false.
          iExitWalkers=-1
          FracLargerDet=1.2_dp
          iReadWalkersRoot=0 
          tShiftonHFPop=.false.
          MaxWalkerBloom=2
          tSearchTau=.true.
          tSearchTauOption = .true.
          tSearchTauDeath = .false.

          t_lanczos_init = .false.
          t_lanczos_store_vecs = .true.
          t_lanczos_orthogonalise = .false.
          t_force_lanczos = .false.
          lanczos_max_restarts = 10
          lanczos_max_vecs = 40
          lanczos_energy_precision = 8
          lanczos_ritz_overlap_precision = 4

          tTimeExit=.false.
          MaxTimeExit=0.0_dp
          tMaxBloom=.false.
          iRestartWalkNum=0
          iWeightPopRead=0.0_dp
          tCheckHighestPop = .true.
          tChangeProjEDet = .true.
          StepsSftImag=0.0_dp
          TauFactor=0.0_dp
          tStartMP1=.false.
          tStartCAS=.false.
          iAnnInterval=1
          tTruncCAS=.false.
          iFullSpaceIter=0
          iDetGroup=2
          tFindDets=.false.
          SinglesBias=1.0_dp
          tSpawnAsDet=.false.
          tDirectAnnihil=.true.
          tRotoAnnihil=.false.
          OccCASorbs=0
          VirtCASorbs=0
          TUnbiasPGeninProjE=.false.
          TRegenExcitgens=.false.
          MemoryFacPart=10.0_dp
          MemoryFacSpawn=3.0_dp
          MemoryFacInit = 0.3_dp
          TStartSinglePart=.true.
          TFixParticleSign=.false.
          TProjEMP2=.false.
          THFRetBias=.false.
          TSignShift=.false.
          tFixedN0 = .false.
          tSkipRef(:) = .false.
          tTrialShift = .false.
          tFixTrial(:) = .false.
          TrialTarget = 0.0
          tAdaptiveShift = .false.
          AdaptiveShiftSigma = 1.0
          AdaptiveShiftF1 = 0.0
          AdaptiveShiftF2 = 1.0
          tAutoAdaptiveShift = .false.
          AdaptiveShiftThresh = 10
          AdaptiveShiftExpo = 1
          AdaptiveShiftCut = -1 !If the user does not specify a value, this will be set to 1.0/HFConn later
          tAAS_MatEle = .false. 
          tInitsEnergy = .false.
          tInitsRDMRef = .false.
          NEquilSteps=0
          NShiftEquilSteps=1000
          TRhoElems=.false.
          TReturnPathMC=.false.
          CLMax=NEl
          PRet=1.0_dp
          TNoAnnihil=.false.
          TFullUnbias=.false.
          TFCIMC=.false.
          tRPA_QBA=.false.
          TMCDets=.false.
          TBinCancel=.false.  
          ScaleWalkers=1.0_dp
          TReadPops=.false.
          iPopsFileNoRead = 0
          iPopsFileNoWrite = 0
          tWalkContGrow=.false.
          StepsSft=100
          SftDamp=10.0_dp
          Tau=0.0_dp
          InitWalkers=3000.0_dp
          NMCyc = -1
          HApp=1
          TMCStar=.false.
          THDiag=.false.
          GrowGraphsExpo=2.0_dp
          TGrowInitGraph=.false.
          AvMCExcits=1.0_dp
          TMaxExcit=.false.
          TFullDiag=.false.
          TSinglesExcitSpace=.false.
          TOneExcitConn=.false.
          TStarTrips=.false.
          TLanczos=.false.
          tDavidson=.false.
          TNoSameExcit=.false.
          TInitStar=.false.
          NoMoveDets=1
          TMoveDets=.false.
          GraphBias=0.99_dp
          TBiasing=.false.
          NDets=400
          Iters=10
          TGraphMorph=.false.
          LinePoints=10
          TSTARSTARS=.false.
          TDIAGNODES=.false.
          STARPROD=.false.
          TCALCHMAT = .false.
          TStar=.false.
          TENERGY = .false.
          NEVAL = 0
          TREAD = .false.
          NBLK = 4
          NKRY = 8
          TBLOCK = .false.
          ICILEVEL = 0
          TNEWEXCITATIONS=.FALSE.
          NWHTAY(:,:) = 0
          NCYCLE = 200
          IOBS=16
          JOBS=16
          KOBS=16
          NDETWORK = 50000
          I_HMAX=0
          I_VMAX=0
          g_MultiWeight(:)=0.0_dp
!This is whether to calculate the expected variance for a MC run when doing full sum (seperate denominator and numerator at present
          TVARCALC(:)=.false.              
          TBIN=.false.
          TVVDISALLOW=.false.
          tMCDirectSum=.FALSE.
          TMPTHEORY=.FALSE.
          tMP2Standalone=.FALSE.
          TMODMPTHEORY=.FALSE.
          G_VMC_PI = 0.95_dp
          G_VMC_SEED = -7
          G_VMC_FAC = 16.0_dp
          TUPOWER=.false.
          G_VMC_EXCITWEIGHT(:)=0.0_dp
          G_VMC_EXCITWEIGHTS(:,:)=0.0_dp
          EXCITFUNCS(:)=.false.
          EXCITFUNCS(10)=.true.
          NPaths = 1
          iActiveBasis=0
          nActiveSpace(:)=0 
          TNPDERIV = .false.
          TMONTE = .false.
          IMCSTEPS = 0
          IEQSTEPS = 0
          BETAEQ = 0.0_dp
          TMCDET = .false.
          MDK(:) = 0
          DETINV = 0
          TSPECDET = .false.
          TTROT=.true.
          BETA = 1000.0_dp
          BETAP=1.0e-4_dp
          TBETAP=.false.
          RHOEPSILON=1.0e-6_dp
          DBETA=-1.0_dp
          GraphEpsilon=0.0_dp
          PGenEpsilon=0.0_dp
          StarConv=1.0e-3_dp
          calcp_sub2vstar=.false.
          calcp_logweight=.false.
          TENPT=.false.
          TLADDER=.false. 
          tDefineDet=.false.
          tTruncInitiator=.false.
          tAddtoInitiator=.false.
          tSTDInits = .false.
          tActivateLAS = .false.
          tLogAverageSpawns = .false.
          spawnSgnThresh = 3.0_dp
          minInitSpawns = 20
          lingerTime = 300
          tTimedDeaths = .false.
          tAVReps = .false.
          tGlobalInitFlag = .false.
          tInitCoherentRule=.true.
          InitiatorWalkNo=3.0_dp
          ErrThresh = 0.3
          tSeniorInitiators =.false.
          SeniorityAge=1.0_dp
          tInitIncDoubs=.false.
          MaxNoatHF=0.0_dp
          HFPopThresh=0
          tSpatialOnlyHash = .false.
          tStoredDets = .false.
          tNeedsVirts=.true.! Set if we need virtual orbitals  (usually set).  Will be unset 
          !(by Calc readinput) if I_VMAX=1 and TENERGY is false

          lNoTriples=.false.
          tReadPopsChangeRef = .false.
          tReadPopsRestart = .false.
          iLogicalNodeSize = 0 !Meaning use the physical node size
          tAllRealCoeff=.false.
          tUseRealCoeffs = .false.
          tRealCoeffByExcitLevel=.false.
          RealCoeffExcitThresh=2
          tRealSpawnCutoff=.false.
          RealSpawnCutoff=1.0e-5_dp
          OccupiedThresh=1.0_dp
          tJumpShift = .true.
!Feb 08 default set.
          IF(Feb08) THEN
              RhoEpsilon=1.0e-8_dp
          ENDIF

          ! Spin Projection defaults
          spin_proj_gamma = 0.1_dp
          tSpinProject  = .false.
          spin_proj_stochastic_yama = .false.
          spin_proj_spawn_initiators = .true.
          spin_proj_no_death = .false.
          spin_proj_interval = 5
          spin_proj_shift = 0.0_dp
          spin_proj_cutoff = 0
          spin_proj_iter_count = 1
          spin_proj_nopen_max = -1
          disable_spin_proj_varyshift = .false.
          tUseProcsAsNodes=.false.

          ! Truncation based on number of unpaired electrons
          tTruncNOpen = .false.

          ! initiators based on number of open orbs
          tSeniorityInits = .false.
          initMaxSenior = 0

          ! keep spawns up to a given seniority + excitation level
          tSpawnSeniorityBased = .false.
          numMaxExLvlsSet = 0
          allocate(maxKeepExLvl(0))
          tLargeMatelSurvive = .false.
          spawnMatelThresh = 1.0_dp

          ! trunaction for spawns/based on spawns
          t_truncate_unocc = .false.
          t_prone_walkers = .false.
          t_activate_decay = .false.

          hash_shift=0
          tUniqueHFNode = .false.

          ! Semi-stochastic and trial wavefunction options.
          tSemiStochastic = .false.
          tCSFCore = .false.
          tDynamicCoreSpace = .false.
          tIntervalSet = .false.
          tStaticCore = .true.
          coreSpaceUpdateCycle = 400
          semistoch_shift_iter = 0
          tTrialWavefunction = .false.
          tDynamicTrial = .false.
          trialSpaceUpdateCycle = 400
          tKP_FCIQMC = .false.
          tLetInitialPopDie = .false.
          tWritePopsNorm = .false.
          pops_norm_unit = 0
          n_init_vecs_ftlm = 20
          n_lanc_vecs_ftlm = 20
          nbeta_ftlm = 100
          delta_beta_ftlm = 0.1_dp
          n_lanc_vecs_sl = 20
          nomega_spectral = 100
          delta_omega_spectral = 0.01_dp
          min_omega_spectral = 0.0_dp
          spectral_broadening = 0.05_dp
          spectral_ground_energy = 0.0_dp
          tIncludeGroundSpectral = .false.
          alloc_popsfile_dets = .false.
          tDetermHFSpawning = .true.

          pParallel = 0.5_dp

          MaxTau = 1.0_dp
          pop_change_min = 50
          tOrthogonaliseReplicas = .false.
          tOrthogonaliseSymmetric = .false.
          orthogonalise_iter = 0
          tReplicaSingleDetStart = .false.

          use_spawn_hash_table = .false.

          ! Continuous time FCIQMC control
          tContTimeFCIMC = .false.
          tContTimeFull = .false.
          cont_time_max_overspawn = 4.0

          tLoadBalanceBlocks = .true.
          tPopsJumpShift = .false.
          calc_seq_no = 1

          ! Superinitiator flags and thresholds
          tAllDoubsInitiators = .false.
          tDelayAllDoubsInits = .false.
          allDoubsInitsDelay = 0
          tAllSingsInitiators = .false.
          tDelayAllSingsInits = .false.
          tSetDelayAllSingsInits = .false.
          tSetDelayAllDoubsInits = .false.
          ! By default, we have one reference for the purpose of all-doubs-initiators      
          nRefs = 1
          maxNRefs = 400
          targetRefPop = 1000
          targetRefPopTol = 80
          tVariableNref = .false.
          tSingleSteps = .true.
          tReadRefs = .false.
          tDelayGetRefs = .false.
          tProductReferences = .false.
          tAccessibleSingles = .false.
          tAccessibleDoubles = .false.
          tSuppressSIOutput = .true.
          nExProd = 2
          NoTypeN = InitiatorWalkNo
          tStrictCoherentDoubles = .false.
          tWeakCoherentDoubles = .true.
          tAvCoherentDoubles = .true.
          superInitiatorLevel = 0
          coherenceThreshold = 0.5
          SIThreshold = 0.95
          SIUpdateInterval = 100
          tAdiActive = .false.
          minSIConnect = 1
          
          ! Walker scaling with energy
          ! do not use scaled walkers
          tEScaleWalkers = .false.
          tLogNumSpawns = .false.
          sFAlpha = 1.0_dp
          sFBeta = 1.0_dp
          sFTag = 0

          ! Epstein-Nesbet second-order correction logicals.
          tEN2 = .false.
          tEN2Init = .false.
          tEN2Truncated = .false.
          tEN2Started = .false.

          ! Giovannis option for RDMs without non-initiators
          tNonInitsForRDMs = .true.

        end subroutine SetCalcDefaults

        SUBROUTINE CalcReadInput()
          USE input_neci
          Use Determinants, only : iActiveBasis, SpecDet, tagSpecDet, tSpecDet, nActiveSpace
          Use Determinants, only : tDefineDet, DefDet, tagDefDet
          use SystemData, only : Beta,nEl
          Use DetCalc, only: iObs, jObs, kObs, B2L, DETINV
          Use DetCalc, only: icilevel, nBlk, nCycle, nEval, nKry, tBlock, tCalcHMat
          Use DetCalc, only: tEnergy, tRead,tFindDets
          use IntegralsData, only: tNeedsVirts,NFROZEN
          use UMatCache, only: gen2CPMDInts
          use FciMCData, only: hash_shift, davidson_ras
          use ras_data
          use global_utilities
          use Parallel_neci, only : nProcessors
          use util_mod, only: addToIntArray
          use LoggingData, only: tLogDets, tWriteUnocc
          IMPLICIT NONE
          LOGICAL eof
          CHARACTER (LEN=100) w
          CHARACTER (LEN=100) input_string
          CHARACTER(*),PARAMETER :: t_r='CalcReadInput'
          integer :: l, i, j, line, ierr
          integer :: tempMaxNoatHF,tempHFPopThresh
          logical :: tExitNow
          integer :: ras_size_1, ras_size_2, ras_size_3, ras_min_1, ras_max_3
          integer :: npops_pert, npert_spectral_left, npert_spectral_right
          integer :: maxKeepNOpenBuf, maxKeepExLvlBuf
          real(dp) :: InputDiagSftSingle

          ! Allocate and set this default here, because we don't have inum_runs
          ! set when the other defaults are set.
          if(.not.allocated(InputDiagSft)) allocate(InputDiagSft(inum_runs))
          InputDiagSft=0.0_dp

          calc: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)

            case("HAMILTONIAN")
                TCALCHMAT = .true.
                IF(item.lt.nitems) THEN
                   call readu(w)
                   select case(w)
                   case("STAR")
                      TSTAR=.TRUE. 
                   case default
                   call report ("Keyword "//trim(w)//                 &
     &                " not recognised",.true.)
                   end select
                ENDIF
            case("ENERGY")
                TENERGY = .true.
                TCALCHMAT = .true.
                tLogDets=.true.
            case("LANCZOS-STORE-VECTORS")
                ! default
                t_lanczos_init = .true.
            case("LANCZOS-STORE-VECTORS-ORTHOGONALISE")
                t_lanczos_init = .true.
                t_lanczos_orthogonalise = .true.
            case("LANCZOS-NO-STORE-VECTORS")
                t_lanczos_init = .true.
                t_lanczos_store_vecs = .true.
            case("LANCZOS-FORCE")
                t_force_lanczos = .true.
            case("LANCZOS-MAX-SUBSPACE-SIZE")
                call readi(lanczos_max_vecs)
            case("LANCZOS-MAX-RESTARTS")
                call readi(lanczos_max_restarts)
            case("LANCZOS-ENERGY-PRECISION")
                call readi(lanczos_energy_precision)
            case("LANCZOS-RITZ-OVERLAP-PRECISION")
                call readi(lanczos_ritz_overlap_precision)


            case("LANCZOS")
!Sets the diagonaliser for the GraphMorph algorithm to be Lanczos
                tLanczos=.true.
            case("EIGENVALUES")
                call readi(NEVAL)
            case("READ")
                TREAD = .true.
            case("COMPLETE")
                NBLK = 0
            case("BLOCKS")
                call geti(NBLK)
            case("KRYLOV")
                call geti(NKRY)
            case("ACCURACY")
                call getf(B2L)
            case("BLOCK")
                call readu(w)
                select case(w)
                case("OFF")
                    TBLOCK = .false.
                case("ON")
                    TBLOCK = .true.
                case default
                    TBLOCK = .true.
                end select
            case("EXCITE")
                call geti(ICILEVEL)
            case("EXCITATIONS")
                call readu(w)
                select case(w)
                case("NEW")
                   TNEWEXCITATIONS=.TRUE.
                case("OLD")
                   TNEWEXCITATIONS=.FALSE.
                case default
                   call inpgetexcitations(NWHTAY(1,1),w)
                end select
            case("STEPS")
                call geti(NCYCLE)
            case("POSITION")
                call geti(IOBS)
                call geti(JOBS)
                call geti(KOBS)
            case("WORKOUT")
                call geti(NDETWORK)
! Using the keyword CONSTRUCTNATORBS includes a calculation of the 1 electron reduced 
! density matrix (1-RDM) as the FCIMC calculation progresses.
! Diagonalisation of this matrix gives linear combinations of the HF orbitals which 
! tend towards the natural orbitals.
! The EQUILSTEPS keyword specifies the number of iterations which must pass before the 
! population of the singles is counted towards the projection energy. 
            case("CONSTRUCTNATORBS")
                CALL Stop_All(t_r,"CONSTRUCTNATORBS option depreciated")
!                tConstructNOs = .true.
            case("ENDCALC")
                exit calc
            case("METHODS")
                if(I_HMAX.ne.0) call report("METHOD already set",.true.)
                I_HMAX=-10
                I_VMAX=1
                tExitNow = .false.
                do while (.not. tExitNow)
                   call read_line(eof)
                   if (eof) then
                      call report("Incomplete input file",.true.)
                   end if
                   call readu(w)
                   select case(trim(w))
                   case("METHOD")
                     I_VMAX=I_VMAX+1
                      NWHTAY(3,I_VMAX)=I_VMAX
                     call inpgetmethod(NWHTAY(1,I_VMAX),NWHTAY(2,I_VMAX),&
     &                I_VMAX)
                    case("EXCITATIONS")
                       call readu(w)
                       call inpgetexcitations(NWHTAY(2,I_VMAX),w)
                     case("CYCLES")
                        call readi(NWHTAY(2,I_VMAX))
                        if ( NWHTAY(1,I_VMAX).ne. -7.and.                  &
       &                     NWHTAY(1,I_VMAX).ne.-19 ) then
                           call report(trim(w)//" only valid for MC "      &
       &                    //"method",.true.)
                       end if
                     case("VERTICES")
                        call geti(NWHTAY(3,I_VMAX))
                     case("MULTIMCWEIGHT")
                        call getf(g_MultiWeight(I_VMAX))
                   case("CALCVAR")
                       if ( NWHTAY(1,I_VMAX).NE.-20 ) then
                           call report("Keyword "//trim(w)//"            &
     &                      only valid for HDIAG routine",.true.)
                       else
                      TVARCALC(I_VMAX)=.true.
                       end if
                    case("DAVIDSON")
                        I_VMAX = I_VMAX + 1
                        tDavidson = .true.
                        call geti(ras_size_1)  ! Number of spatial orbitals in RAS1.
                        call geti(ras_size_2)  ! Number of spatial orbitals in RAS2.
                        call geti(ras_size_3)  ! Number of spatial orbitals in RAS3.
                        call geti(ras_min_1)  ! Min number of electrons (alpha and beta) in RAS1 orbs. 
                        call geti(ras_max_3)  ! Max number of electrons (alpha and beta) in RAS3 orbs.
                        davidson_ras%size_1 = int(ras_size_1,sp)
                        davidson_ras%size_2 = int(ras_size_2,sp)
                        davidson_ras%size_3 = int(ras_size_3,sp)
                        davidson_ras%min_1 = int(ras_min_1,sp)
                        davidson_ras%max_3 = int(ras_max_3,sp)
                    case("ENDMETHODS")
                       tExitNow = .true.

                    case default
                        write(6,*) 'REPORT' // trim(w)
                       !call report ("Keyword "//trim(w)//" not recognized",.true.)
                    end select
 
                end do

                       
            case("METHOD")

                if(I_HMAX.ne.0) call report("METHOD already set",.true.)
                call inpgetmethod(I_HMAX,NWHTAY(1,1),0)

            case("RDMSAMPLINGITERS")
                !How many iterations do we want to sample the RDM for?
                call readi(iSampleRDMIters)
            case("CYCLES")
                call readi(NWHTAY(1,1))
                if ( I_HMAX .ne. -7.and.                              &
     &               I_HMAX .ne. -19) then
                    call report(trim(w)//" only valid for MC "        &
     &                 //"method",.true.)
                end if

             case("INITS-RDM")
                ! only take into account initiators when calculating RDMs
                tNonInitsForRDMs = .false.
            case("VVDISALLOW")
                TVVDISALLOW=.TRUE.
            case("MCDIRECTSUM")
                TMCDIRECTSUM=.TRUE.
            case("EPSTEIN-NESBET")
!True if Epstein-Nesbet PT rather than Moller-Plesset
                tENPT=.TRUE.
            case("LADDER")
                tLadder=.TRUE.
            case("MODMPTHEORY")
                TMODMPTHEORY=.TRUE.
            case("MPTHEORY")
                TMPTHEORY=.TRUE.
                if (item.lt.nitems) then
                    ! Something else remains on the line.
                    call readu(w)
                    select case(w)
                    case("ONLY")
                        tMP2Standalone=.true.
                    end select
                end if
!                if ( I_HMAX .ne. -7.and.
!     &               I_HMAX .ne. -19) then
!                    call report(trim(w)//" only valid for MC " 
!     &                 //"method",.true.)
!                end if
            case("MAXVERTICES")
                if ( I_VMAX .ne. 0 ) then
                   call report("Cannot reset MAXVERTICES",.true.)
                endif 
                call readi(I_VMAX)
            case("IMPORTANCE")
                call readf(G_VMC_PI)
!                if ( I_HMAX .ne. -7 ) then
!                    call report(trim(w)//" only valid for MC "  
!     &                //"method",.true.)
!                end if
            case("SEED")
                call readi(G_VMC_SEED)
!                if ( I_HMAX .ne. -7 ) then
!                    call report(trim(w)//" only valid for MC " 
!     &                //"method",.true.)
!                end if
            case("BIAS")
                call readf(G_VMC_FAC)

!                if ( I_HMAX .ne. -7 ) then
!                    call report(trim(w)//" only valid for MC " 
!     &              //" method",.true.)
!                end if
            case("STARCONVERGE")
                call readf(STARCONV)
                if((NWHTAY(1,I_VMAX).ne.0).and.(NWHTAY(1,I_VMAX).ne.-21)&
     &               .and.(NWHTAY(1,I_VMAX).ne.-9)) then
                    call report(trim(w)//" only valid for STAR method",.true.)
                end if
            case("UFORM-POWER")
                TUPOWER=.true.
            case("CHEMPOTWEIGHTING")
                call readf(g_VMC_ExcitWeights(1,1))
                call readf(g_VMC_ExcitWeights(2,1))
                call readf(G_VMC_EXCITWEIGHT(1))
                DO l=1,6
                  IF(EXCITFUNCS(l)) THEN
                      call report(trim(w)//" only valid if another weighting scheme not specified",.true.)
                  ENDIF
                ENDDO
                EXCITFUNCS(4)=.true.
            case("CHEMPOT-TWOFROM")
                call readf(g_VMC_ExcitWeights(1,1))
                call readf(g_VMC_ExcitWeights(2,1))
                call readf(g_VMC_ExcitWeights(3,1))
                call readf(G_VMC_EXCITWEIGHT(1))
                DO l=1,6
                  IF(EXCITFUNCS(l)) THEN
                      call report(trim(w)//" only valid if "          &
     &             //" another weighting scheme not specified",.true.)
                  ENDIF
                ENDDO
                EXCITFUNCS(5)=.true.
            case("POLYEXCITWEIGHT")
                call readf(g_VMC_ExcitWeights(1,1))
                call readf(g_VMC_ExcitWeights(2,1))
                call readf(g_VMC_ExcitWeights(3,1))
                call readf(G_VMC_EXCITWEIGHT(1))
                DO l=1,6
                  IF(EXCITFUNCS(l)) THEN
                      call report(trim(w)//" only valid if "          &
     &             //" another weighting scheme not specified",.true.)
                  ENDIF
                ENDDO
                EXCITFUNCS(2)=.true.
            case("POLYEXCITBOTH")
                call readf(g_VMC_ExcitWeights(1,1))
                call readf(g_VMC_ExcitWeights(2,1))
                call readf(g_VMC_ExcitWeights(3,1))
                call readf(g_VMC_ExcitWeights(4,1))
                call readf(G_VMC_EXCITWEIGHT(1))
                DO l=1,6
                  IF(EXCITFUNCS(l)) THEN
                      call report(trim(w)//" only valid if "          &
     &             //" another weighting scheme not specified",.true.)
                  ENDIF
                ENDDO
                EXCITFUNCS(3)=.true.
            case("EXCITWEIGHTING")
                           write(6,*) '---------------->excitweighting'
                             call neci_flush(6)  
                call readf(g_VMC_ExcitWeights(1,1))
                call readf(g_VMC_ExcitWeights(2,1))
                call readf(G_VMC_EXCITWEIGHT(1))
                IF(item.lt.nitems) call readf(g_VMC_ExcitWeights(3,1))
                DO l=1,6
                  IF(EXCITFUNCS(l)) THEN
                      call report(trim(w)//" only valid if "          &
     &             //" another weighting scheme not specified",.true.)
                  ENDIF
                ENDDO
                EXCITFUNCS(1)=.true.

             case("STEPEXCITWEIGHTING")
!This excitation weighting involves a step function between the virtual and occupied electon 
!manifold (i.e. step is at the chemical potential)
!When choosing an electron to move, the probability of selecting it is 1 if the electron is in the virtual manifold
!and (g_VMC_ExcitWeights(1,1) if in the virtual manifold. When choosing where to excite to, the situation is reversed, 
!and the probability of selecting it is
!1 if the electron is in the occupied manifold and g_VMC_ExcitWeights(2,1) if in the occupied manifold. 
!U-weighting is the third parameter as before.
                call readf(g_VMC_ExcitWeights(1,1))
                call readf(g_VMC_ExcitWeights(2,1))
                call readf(G_VMC_EXCITWEIGHT(1))
                DO l=1,6
                    IF(EXCITFUNCS(l)) THEN
                      call report(trim(w)//" only valid if "          &
     &             //" another weighting scheme not specified",.true.)
                  ENDIF
                ENDDO
                EXCITFUNCS(6)=.true.
            case("PATHS")
                call readu(w)
                select case(w)
                case("ALL")
                   NPATHS=-1
                case("ACTIVE")
                   if(item.lt.nitems) then
                      call readu(w)
                      select case(w)
                      case("ORBITALS")
                         NPATHS=-3
                         call readi(nActiveSpace(1))
                         call readi(nActiveSpace(2))
                      case("SETS")
                         NPATHS=-2
                         call readi(nActiveSpace(1))
                         call readi(nActiveSpace(2))
                      case default
                          call report(trim(w)//" unknown",.true.)
                      end select
                   else
                     NPATHS=-2
                     nActiveSpace(:)=0 
                   endif
                case default
                   call reread(-1)
                   call geti(NPATHS)
                end select
                iActiveBasis=nPaths
            case("ALLPATHS")
                NPATHS = -1
            case("DERIV")
                TNPDERIV = .true.
               if (DBETA .lt. 0 ) then
                  call report("Only calculate energy with derivatives"&
     &            //" if delta_beta positive",.true.)
                   TNPDERIV = .false.
               end if
            case("CIMC")
                TMONTE = .true.
            case("MCSTEPS")
                call readi(IMCSTEPS)
                if ( .not. TMONTE ) then
                    call report(trim(w)//" only relevant if CI space" &
     &              //" monte carlo is performed.",.true.)
                end if
            case("EQSTEPS")
                call readi(IEQSTEPS)
                if ( .not. TMONTE ) then
                    call report(trim(w)//" only relevant if CI space" &
     &              //" monte carlo is performed.",.true.)
                end if
            case("BETAEQ")
                call readf(BETAEQ)
                if ( .not. TMONTE ) then
                    call report(trim(w)//" only relevant if CI space" &
     &              //" monte carlo is performed.",.true.)
                end if
            case("DETSYM")
                TMCDET = .true.
                do I = 1,5
                  call readi(MDK(I))
                end do
                if ( .not. TMONTE ) then
                    call report(trim(w)//" only relevant if CI space" &
     &               //" monte carlo is performed.",.true.)
                end if
            case("DETINV")
                call readi(DETINV)
            case("INSPECT")
                TSPECDET = .true.
                ALLOCATE(SPECDET(NEL-NFROZEN),STAT=ierr)
                CALL LogMemAlloc('SPECDET',NEL-NFROZEN,4,t_r,tagSPECDET,ierr)
                SPECDET(1)=0
                if(item.lt.nitems) then
!Cannot specify frozen core orbitals if want to specify a determinant?
!This is because LOGREAD has not been called yet, and NFROZEN not specified.
                   do I = 1,NEL-NFROZEN
                     call geti(SPECDET(I))
                   end do
                endif
            case("LOGICALNODESIZE")
!Sets the Logical node size to this value, rather than using the physical node size.
!Use to simulate a multi-node process on a single node.
               call geti(iLogicalNodeSize)
            case("DEFINEDET")
!This defines the reference determinant to be that specified in the input here, rather than the determinant 
!chosen from the lowest energy orbitals.
!The 'HF' energy calculated should the be that of the defined determinant.
                tDefineDet=.true.
                if(.not.allocated(DefDet)) then
                  ALLOCATE(DefDet(NEl),stat=ierr)
                  CALL LogMemAlloc('DefDet',NEl,4,t_r,tagDefDet,ierr)
                end if
                DefDet(:)=0
                do i=1,NEl
                    call geti(DefDet(i))
                enddo

            case("MULTIPLE-INITIAL-REFS")
                tMultipleInitialRefs = .true.
                allocate(initial_refs(nel, inum_runs), stat=ierr)
                initial_refs = 0

                do line = 1, inum_runs
                    call read_line(eof)
                    do i = 1, nel
                        call geti(initial_refs(i, line))
                    end do
                end do

            case("MULTIPLE-INITIAL-STATES")
                tMultipleInitialStates = .true.
                allocate(initial_states(nel, inum_runs), stat=ierr)
                initial_states = 0

                do line = 1, inum_runs
                    call read_line(eof)
                    do i = 1, nel
                        call geti(initial_states(i, line))
                    end do
                end do

            case("FINDGUIDINGFUNCTION")
! At the end of a calculation, this keyword sets the spawning calculation to print out the iGuideDets
! most populated determinants, to be read in as a guiding (or annihilating) function in a following calculation.
                CALL Stop_All(t_r,"FINDGUIDINGFUNCTION option depreciated")
!                tFindGuide=.true.
!                call geti(iGuideDets)

            case("USEGUIDINGFUNCTION")
! This keyword sets the calculationg to read in a guiding function from a file GUIDINGFUNC.  This function then sits
! in the back of a calculation, able to annihilate particle, but not allowed to spawn or die.
! iInitGuideParts specifies how many walkers start on the HF determinant, and the remaining determinants are populated
! based on their populations from the previous calculation relative to the HF.
                CALL Stop_All(t_r,"USEGUIDINGFUNCTION option depreciated")
!                tUseGuide=.true.
!                call geti(iInitGuideParts)

            case("SPAWNDOMINANTONLY")
! This option sets the calculation to read in from a file DOMINANTDETS.  The determinants from this file make up a list of 
! determinants on which spawning is allowed for the excitation levels included.
! Spawning onto determinants that have the listed excitation level, but are not read in from this file is forbidden.
                CALL Stop_All(t_r,"SPAWNDOMINANTONLY option depreciated")
                
!                tSpawnDominant=.true.

            case("PRINTDOMINANTDETS")
! This option finds the iNoDominantDets most populated determinants with excitation level 
!between MinExcDom and MaxExcDom and
! prints them to a file named DOMINANTDETS.  This can be later read in as the allowed determinants 
!for spawing in a restricted calc.
                CALL Stop_All(t_r,"PRINTDOMINANTDETS option depreciated")
!                tPrintDominant=.true.
!                call geti(iNoDominantDets)
!                call geti(MinExcDom)
!                call geti(MaxExcDom)

            case("PRINTDOMSPINCOUPLED")
! This option finds the iNoDominantDets most populated determinants with excitation level between 
!MinExcDom and MaxExcDom and
! prints them to a file named DOMINANTDETS.  This can be later read in as the allowed determinants 
!for spawing in a restricted calc.
                CALL Stop_All(t_r,"PRINTDOMSPINCOUPLED option depreciated")
!                if(item.lt.nitems) then
!                    call readu(w)
!                    select case(w)
!                    case("OFF")
!                        tNoDomSpinCoup=.true.
!                    end select
!                else
!                    tNoDomSpinCoup=.false.
!                end if

            case("STARMINORDETERMINANTS")
                CALL Stop_All(t_r,"STARMINORDETERMINANTS option depreciated")
!                tMinorDetsStar=.true.
! This option goes along with the SPAWNDOMINANTONLY option.  However, if this keyword is present, 
!spawning onto determinants that are not in the 
! dominant determinants list is allowed, however once spawned into this "insignificant" region, 
!walkers may only spawn back onto the determinant 
! from which they came.  In the mean time, walkers on "insignificant" determinants may live/die 
!and annihilate like any others.
! This is a second order perturbation to the SPAWNDOMINANTONLY approximation.

            case("TROTTER")
                TTROT = .true.
            case("BETA")
                call getf(BETA)
            case("BETAOVERP")
                call getf(BETAP)
                TBETAP = .true.
            case("TIMESTEPS")
                BETAP = 0
                call geti(I_P)
                if ( TBETAP ) then
                    call report("Warning - declared beta/p and p. Using p.",.true.)
                end if
            case("DELTABETA")
                call getf(DBETA)
            case("RHOEPSILON")
                call getf(RHOEPSILON)
            case("GRAPHEPSILON")
                call getf(GraphEpsilon)
            case("PGENEPSILON")
                call getf(PGenEpsilon)
!This indicates the number of times the eigenvalues of the star matrix should be evaluated to 
!achieve the linear approximation when STARSTARS set,
            case("LINEPOINTSSTAR")
                call geti(LinePoints)
!This is the number of vertices in the Graph Morph graph. Alternativly, it is used by ResumFCIMC, as the 
!size of their graphs. Then, if it is negative, the graph is all possible connections
            case("GRAPHSIZE")
                call geti(NDets)
!This is the number of times to systematically improve the Graph using the morphing algorithm
            case("ITERATIONS")
                call geti(Iters)
            case("GRAPHBIAS")
                call getf(GraphBias)
                TBiasing=.true.
            case("MOVEDETS")
                call geti(NoMoveDets)
                TMoveDets=.true.
            case("INITSTAR")
                TInitStar=.true.
            case("NOSAMEEXCIT")
                TNoSameExcit=.true.
            case("NOTRIPLES")
                lNoTriples=.true.
            case("ONEEXCITCONN")
!This means that determinants can only be attached to each other if they differ by one excitation level from HF
                TOneExcitConn=.true.
            case("SINGLESEXCITSPACE")
!This means that the space accessible to the morphing algorithm is the space of single 
!excitations of the determinants in the graph.
                TSinglesExcitSpace=.true.
            case("FULLDIAGTRIPS")
!When constructing a star of triples from each double star, then this tag results in a full 
!diagonalisation of this matrix.
                TFullDiag=.true.
            case("AVERAGEMCEXCITS")
! This sets the average number of spawning attempts from each walker.
                call getf(AvMCExcits)
            case("GROWINITGRAPH")
!In GraphMorph, this means that the initial graph is grown non-stochastically from the excitations 
!of consecutive determinants
                TGrowInitGraph=.true.
            case("GROWGRAPHSEXPO")
!In GraphMorph, this is the exponent to which the components of the excitation vector and eigenvector 
!will be raised to turn them into probabilities.
                call getf(GrowGraphsExpo)
            case("HAPP")
!For graph MC, this indicates the number of local applications of the hamiltonian to random determinants 
!before the trial eigenvector is updated
                call geti(HApp)
            case("MAXEXCIT")
!This imposes a maximum excitation level to the space that GraphMorph can explore. Note: FCIMC uses EXCIT 
!to indicate a maximum excit level.
                TMaxExcit=.true.
                call geti(iMaxExcitLevel)
            case("INITWALKERS")
!For FCIMC, this is the number of walkers to start with
                call getf(InitWalkers)
            case("TOTALWALKERS")
!This is now input as the total number, rather than the number per processor, and it is changed to the number per processor here.
                call getf(InitWalkers)
                InitWalkers=NINT(REAL(InitWalkers,dp)/REAL(nProcessors,dp),int64)
            case("TIME")
                !Input the desired runtime (in MINUTES) before exiting out of the MC.
                call getf(MaxTimeExit)
                MaxTimeExit=MaxTimeExit*60.0_dp    !Change straightaway so that MaxTimeExit corresponds to SECONDS!
                tTimeExit=.true.
            case("MAXNOATHF")
!If the number of walkers at the HF determinant reaches this number, the shift is allowed to change. 
!(This is the total number across all processors).                
!If a second integer is present, this determinants the threshhold for the HF population.  If the HF 
!population drops below MaxNoatHF-HFPopThresh, the
!number of walkers is allowed to grow again until MaxNoatHF is reachieved.
!Without the second integer, MaxNoatHF-HFPopThresh=0, and the HF population can drop to 0 without any consequences.
                call geti(tempMaxNoatHF)
                MaxNoatHF=tempMaxNoatHF
                if(item.lt.nitems) then
                    call geti(tempHFPopThresh)
                    HFPopThresh=tempHFPopThresh
                else
                    HFPopThresh=int(MaxNoatHF,int64) 
                end if
            case("HASH_SHIFT")
                call geti(hash_shift)
            case("NMCYC")
!For FCIMC, this is the number of MC cycles to perform
                call geti(NMCyc)
            case("DIAGSHIFT")
!For FCIMC, this is the amount extra the diagonal elements will be shifted. This is proportional to the deathrate of 
!walkers on the determinant
                if (nitems == 2) then
                    call getf(InputDiagSftSingle)
                    InputDiagSft = InputDiagSftSingle
                else
                    if (inum_runs /= nitems-1) call stop_all(t_r, "The number of initial shifts input is not equal to &
                                                                  &the number of replicas being used.")
                    do i = 1, inum_runs
                        call getf(InputDiagSft(i))
                    end do
                end if

            case("TAUFACTOR")
!For FCIMC, this is the factor by which 1/(HF connectivity) will be multiplied by to give the timestep for the calculation.
                tSearchTau=.false.  !Tau is set, so don't search for it.
                tSearchTauOption = .false.
                call getf(TauFactor)
            case("TAU")
                ! For FCIMC, this can be considered the timestep of the
                ! simulation. It is a constant which will increase/decrease
                ! the rate of spawning/death for a given iteration.
                call getf(Tau)

                ! If SEARCH is provided, use this value as the starting value
                ! for tau searching
                if (item < nitems) then
                    call readu(w)
                    select case(w)
                    case("SEARCH")
                        tSearchTau = .true.
                        tSearchTauOption = .true.
                    case default
                        tSearchTau = .false.
                        tSearchTauOption = .false.
                    end select
                else
                    tSearchTau = .false.
                    tSearchTauOption = .false.
                end if

            case("MAX-TAU")
                ! For tau searching, set a maximum value of tau. This places
                ! a limit to prevent craziness at the start of a calculation
                call getf(MaxTau)

            case("MIN-TAU")
                ! introduce a minimum value of tau for the automated 
                ! tau-search to limit the lower bound of the automated 
                ! tau-search 
                t_min_tau = .true. 

                if (item < nitems) then
                    call getf(min_tau_global)
                end if

                ! also only use that if the automated tau-search is used 
                ! so enable the automated tau search here 
                tSearchTau = .true. 
                tSearchTauOption = .true.

            case("KEEPTAUFIXED")
                ! option for a restarted run to keep the tau, read in from the 
                ! POPSFILE and other parameters, as pSingles, pParallel 
                ! fixed for the remainder of the run, even if we keep 
                ! growing the walkers 
                t_keep_tau_fixed = .true.

                ! here i need to turn off the tau-search option
                tSearchTau = .false.
                tSearchTauOption = .false.

            case("HIST-TAU-SEARCH","NEW-TAU-SEARCH")
                ! [Werner Dobrautz, 4.4.2017:]
                ! the new tau search method using histograms of the 
                ! H_ij / pgen ratio and integrating the histograms up to 
                ! a certain value, to obtain the time-step and not using 
                ! only the worst case H_ij / pgen ration
                
                ! this option has 3 possible input parameters: 
                ! 1) the integration cutoff in percentage [0.999 default]
                ! 2) the number of bins used [100000 default]
                ! 3) the upper bound of the bins [10000.0 default]
                t_hist_tau_search = .true.
                t_hist_tau_search_option = .true.
                t_fill_frequency_hists = .true.

                ! turn off the other tau-search, if by mistake both were 
                ! chosen! 
                if (tSearchTau .or. tSearchTauOption) then 
                   write(iout, &
                       '("(WARNING: both the histogramming and standard tau&
                       &-search option were chosen! TURNING STANDARD VERSION OFF!")')
                   tSearchTau = .false.
                   tSearchTauOption = .false.
               end if

               if (item < nitems) then
                   call getf(frq_ratio_cutoff)
               end if

               if (item < nitems) then 
                   call geti(n_frequency_bins)
                   
                   ! check that not too many bins are used which may crash 
                   ! the MPI communication of the histograms! 
                   if (n_frequency_bins > 1000000) then 
                       write(iout, &
                           '("WARNING: maybe too many bins used for the &
                           &histograms! This might cause MPI problems!")')
                   end if
               end if

               if (item < nitems) then 
                   call getf(max_frequency_bound) 
               end if

            case("RESTART-HIST-TAU-SEARCH", "RESTART-NEW-TAU-SEARCH")
                ! [Werner Dobrautz 5.5.2017:]
                ! a keyword, which in case of a continued run from a 
                ! previous hist-tau-search run restarts the histogramming 
                ! tau-search anyway, in case the tau-search is not yet 
                ! converged enough
                t_restart_hist_tau = .true.

                if (item < nitems) then 
                    call geti(hist_search_delay) 
                end if

            case ("TEST-HIST-TAU", "LESS-MPI-HEAVY")
                ! test a change to the tau search to avoid those nasty 
                ! MPI communications each iteration
                t_test_hist_tau = .true. 

            case("TRUNCATE-SPAWNS")
                ! [Werner Dobrautz, 4.4.2017:]
                ! in combination with the above HIST-TAU-SEARCH option I 
                ! also introduced a truncation keyword for spawning events
                ! which are missed by the integrated time-step. 
                ! to limit the effect of these possible large blooms I 
                ! implemented a truncation of those. But this might be an 
                ! uncontrolled approximation, so be careful! 
                t_truncate_spawns = .true. 
                if (item < nitems) then 
                    call getf(n_truncate_spawns)
                    if(item < nitems) then
                       call readu(w)
                       select case(w)
                       case("UNOCC")
                          t_truncate_unocc = .true.
                       case("MULTI")
                          t_truncate_multi = .true.
                       case default
                          t_truncate_unocc = .false.
                       end select
                    endif
                end if

             case("PRONE-DETERMINANTS")
                ! when close to running out of memory, start culling 
                ! the population by removing lonely spawns
                t_prone_walkers = .true.
                
            case("MIX-RATIOS")
                ! pablos idea: mix the old and new contributions and not 
                ! only take the new ones, since we are doing a stochastic 
                ! process now, maybe make that the default behavior..
                t_mix_ratios = .true.

                if (item < nitems) then
                    call getf(mix_ratio)
                else
                    ! if no additional input default it to 0.7
                    mix_ratio = 0.7_dp
                end if

            case("MATRIX-CUTOFF")
                ! [Werner Dobrautz 26.4.2017:]
                ! introduce a matrix element cutoff similar to the 
                ! UMATEPS quantity when ignoring 2-body integrals
                t_matele_cutoff = .true.
                if (item < nitems) then
                    call getf(matele_cutoff)
                else
                    ! does this work? is umateps already defined properly?
                    matele_cutoff = UMATEPS
                    print *, "TEST cutoff: ", matele_cutoff
                end if


            case("MAXWALKERBLOOM")
                !Set the maximum allowed walkers to create in one go, before reducing tau to compensate.
                call getf(MaxWalkerBloom)
            case("SHIFTDAMP")
!For FCIMC, this is the damping parameter with respect to the update in the DiagSft value for a given number of MC cycles.
                call getf(SftDamp)

            case("LINSCALEFCIMCALGO")
                ! Use the linear scaling FCIMC algorithm
                ! This option is now deprecated, as it is default.
                write(iout, '("WARNING: LINSCALEFCIMCALGO option has been &
                              &deprecated, and now does nothing")')
                !call stop_all(t_r, "Option LINSCALEFCIMCALGO deprecated")

            case("PARTICLE-HASH-MULTIPLIER")
                ! Determine the absolute length of the hash table relative to
                ! the target number of walkers (InitWalkers)
                ! 
                ! By default this value is 0.7 (see above)
                call getf(HashLengthFrac)

            case("SEMI-STOCHASTIC")
                tSemiStochastic = .true.
                ! If there is ane extra item, it should specify that we turn
                ! semi-stochastic on later.
                if (item < nitems) then
                   if(item < nitems) &
                      call geti(semistoch_shift_iter)
                    tSemiStochastic = .false.
                    tStartCoreGroundState = .false.
                end if
            case("CSF-CORE")
                if(item.lt.nitems) then
                   call geti(STOT)
                else
                   STOT=0
                endif
                tCSFCore = .true.
                tCSF = .true.
                LMS = STOT
            case("DOUBLES-CORE")
                ss_space_in%tDoubles = .true.
            case("HF-CONN-CORE")
                ss_space_in%tDoubles = .true.
                ss_space_in%tHFConn = .true.
            case("CAS-CORE")
                ss_space_in%tCAS = .true.
                tSpn = .true.
                call geti(ss_space_in%occ_cas)  !Number of electrons in CAS 
                call geti(ss_space_in%virt_cas)  !Number of virtual spin-orbitals in CAS
            case("RAS-CORE")
                ss_space_in%tRAS = .true.
                call geti(ras_size_1)  ! Number of spatial orbitals in RAS1.
                call geti(ras_size_2)  ! Number of spatial orbitals in RAS2.
                call geti(ras_size_3)  ! Number of spatial orbitals in RAS3.
                call geti(ras_min_1)  ! Min number of electrons (alpha and beta) in RAS1 orbs. 
                call geti(ras_max_3)  ! Max number of electrons (alpha and beta) in RAS3 orbs.
                ss_space_in%ras%size_1 = int(ras_size_1,sp)
                ss_space_in%ras%size_2 = int(ras_size_2,sp)
                ss_space_in%ras%size_3 = int(ras_size_3,sp)
                ss_space_in%ras%min_1 = int(ras_min_1,sp)
                ss_space_in%ras%max_3 = int(ras_max_3,sp)
            case("OPTIMISED-CORE")
                ss_space_in%tOptimised = .true.
            case("OPTIMISED-CORE-CUTOFF-AMP")
                ss_space_in%opt_data%tAmpCutoff = .true.
                ss_space_in%opt_data%ngen_loops = nitems - 1
                allocate(ss_space_in%opt_data%cutoff_amps(ss_space_in%opt_data%ngen_loops))
                do I = 1, ss_space_in%opt_data%ngen_loops
                    call getf(ss_space_in%opt_data%cutoff_amps(I))
                end do
            case("OPTIMISED-CORE-CUTOFF-NUM")
                ss_space_in%opt_data%tAmpCutoff = .false.
                ss_space_in%opt_data%ngen_loops = nitems - 1
                allocate(ss_space_in%opt_data%cutoff_nums(ss_space_in%opt_data%ngen_loops))
                do I = 1, ss_space_in%opt_data%ngen_loops
                    call geti(ss_space_in%opt_data%cutoff_nums(I))
                end do
            case("FCI-CORE")
                ss_space_in%tFCI = .true.
            case("HEISENBERG-FCI-CORE")
                ss_space_in%tHeisenbergFCI = .true.
            case("HF-CORE")
                ss_space_in%tHF = .true.
            case("POPS-CORE")
                ss_space_in%tPops = .true.
                call geti(ss_space_in%npops)
            case("POPS-CORE-APPROX")
                ss_space_in%tPops = .true.
                ss_space_in%tApproxSpace = .true.
                call geti(ss_space_in%npops)
                if(item.lt.nitems) then
                   call geti(ss_space_in%nApproxSpace)
                endif
            case("MP1-CORE")
                ss_space_in%tMP1 = .true.
                call geti(ss_space_in%mp1_ndets)
            case("READ-CORE")
                ss_space_in%tRead = .true.
                ss_space_in%read_filename = 'CORESPACE'
            case("MAX-CORE-SIZE")
                ss_space_in%tLimitSpace = .true.
                call geti(ss_space_in%max_size)
            case("DYNAMIC-CORE")
                tDynamicCoreSpace = .true.
                tIntervalSet = .true.
                if( item < nitems) call geti(coreSpaceUpdateCycle)
            case("STATIC-CORE")
                tDynamicCoreSpace = .false.
                tIntervalSet = .true.
            case("STOCHASTIC-HF-SPAWNING")
                tDetermHFSpawning = .false.

            case("TRIAL-WAVEFUNCTION")
                if (item == nitems) then
                    tTrialWavefunction = .true.
                else if (item < nitems) then
                    tStartTrialLater = .true.
                    call geti(trial_shift_iter)
                end if

            case("NUM-TRIAL-STATES-CALC")
                call geti(ntrial_ex_calc)
            case("QMC-TRIAL-WF")
               qmc_trial_wf = .true.
            case("MAX-TRIAL-SIZE")
                trial_space_in%tLimitSpace = .true.
                call geti(trial_space_in%max_size)
            case("MP1-TRIAL")
                trial_space_in%tMP1 = .true.
                call geti(trial_space_in%mp1_ndets)
            case("DOUBLES-TRIAL")
                trial_space_in%tDoubles = .true.
             case("DYNAMIC-TRIAL")
                ! Update the trial wavefunction periodically
                tDynamicTrial = .true.
                if(item < nitems) call geti(trialSpaceUpdateCycle)
            case("CAS-TRIAL")
                trial_space_in%tCAS = .true.
                tSpn = .true.
                call geti(trial_space_in%occ_cas) ! Number of electrons in CAS 
                call geti(trial_space_in%virt_cas) ! Number of virtual spin-orbitals in CAS
            case("RAS-TRIAL")
                trial_space_in%tRAS = .true.
                call geti(ras_size_1)  ! Number of spatial orbitals in RAS1.
                call geti(ras_size_2)  ! Number of spatial orbitals in RAS2.
                call geti(ras_size_3)  ! Number of spatial orbitals in RAS3.
                call geti(ras_min_1)  ! Min number of electrons (alpha and beta) in RAS1 orbs. 
                call geti(ras_max_3)  ! Max number of electrons (alpha and beta) in RAS3 orbs.
                trial_space_in%ras%size_1 = int(ras_size_1,sp)
                trial_space_in%ras%size_2 = int(ras_size_2,sp)
                trial_space_in%ras%size_3 = int(ras_size_3,sp)
                trial_space_in%ras%min_1 = int(ras_min_1,sp)
                trial_space_in%ras%max_3 = int(ras_max_3,sp)
            case("OPTIMISED-TRIAL")
                trial_space_in%tOptimised = .true.
            case("OPTIMISED-TRIAL-CUTOFF-AMP")
                trial_space_in%opt_data%tAmpCutoff = .true.
                trial_space_in%opt_data%ngen_loops = nitems - 1
                allocate(trial_space_in%opt_data%cutoff_amps(trial_space_in%opt_data%ngen_loops))
                do I = 1, trial_space_in%opt_data%ngen_loops
                    call getf(trial_space_in%opt_data%cutoff_amps(I))
                end do
            case("OPTIMISED-TRIAL-CUTOFF-NUM")
                trial_space_in%opt_data%tAmpCutoff = .false.
                trial_space_in%opt_data%ngen_loops = nitems - 1
                allocate(trial_space_in%opt_data%cutoff_nums(trial_space_in%opt_data%ngen_loops))
                do I = 1, trial_space_in%opt_data%ngen_loops
                    call geti(trial_space_in%opt_data%cutoff_nums(I))
                end do
            case("HF-TRIAL")
                trial_space_in%tHF = .true.
            case("POPS-TRIAL")
                trial_space_in%tPops = .true.
                call geti(trial_space_in%npops)
            case("POPS-TRIAL-APPROX")
                trial_space_in%tPops = .true.
                trial_space_in%tApproxSpace = .true.
                call geti(trial_space_in%npops)
            case("READ-TRIAL")
                trial_space_in%tRead = .true.
                trial_space_in%read_filename = 'TRIALSPACE'
            case("FCI-TRIAL")
                trial_space_in%tFCI = .true.
            case("HEISENBERG-FCI-TRIAL")
                trial_space_in%tHeisenbergFCI = .true.
            case("TRIAL-BIN-SEARCH")
                tTrialHash = .false.
                write(iout,*) "WARNING: Disabled trial hashtable. Load balancing "//&
                     "is not supported in this mode and might break the trial energy"
            case("TRIAL-ESTIMATE-REORDER")
                allocate(trial_est_reorder(inum_runs))
                trial_est_reorder = 0
                do i = 1, inum_runs
                    call geti(trial_est_reorder(i))
                end do
            case("TRIAL-INIT-REORDER")
                allocate(trial_init_reorder(inum_runs))
                trial_init_reorder = 0
                do i = 1, inum_runs
                    call geti(trial_init_reorder(i))
                end do
            case("MP1-INIT")
                init_trial_in%tMP1 = .true.
                call geti(init_trial_in%mp1_ndets)
            case("DOUBLES-INIT")
                init_trial_in%tDoubles = .true.
            case("CAS-INIT")
                init_trial_in%tCAS = .true.
                tSpn = .true.
                call geti(init_trial_in%occ_cas) ! Number of electrons in CAS
                call geti(init_trial_in%virt_cas) ! Number of virtual spin-orbitals in CAS
            case("RAS-INIT")
                init_trial_in%tRAS = .true.
                call geti(ras_size_1)  ! Number of spatial orbitals in RAS1.
                call geti(ras_size_2)  ! Number of spatial orbitals in RAS2.
                call geti(ras_size_3)  ! Number of spatial orbitals in RAS3.
                call geti(ras_min_1)  ! Min number of electrons (alpha and beta) in RAS1 orbs.
                call geti(ras_max_3)  ! Max number of electrons (alpha and beta) in RAS3 orbs.
                init_trial_in%ras%size_1 = int(ras_size_1,sp)
                init_trial_in%ras%size_2 = int(ras_size_2,sp)
                init_trial_in%ras%size_3 = int(ras_size_3,sp)
                init_trial_in%ras%min_1 = int(ras_min_1,sp)
                init_trial_in%ras%max_3 = int(ras_max_3,sp)
            case("OPTIMISED-INIT")
                init_trial_in%tOptimised = .true.
            case("OPTIMISED-INIT-CUTOFF-AMP")
                init_trial_in%opt_data%tAmpCutoff = .true.
                init_trial_in%opt_data%ngen_loops = nitems - 1
                allocate(init_trial_in%opt_data%cutoff_amps(init_trial_in%opt_data%ngen_loops))
                do I = 1, init_trial_in%opt_data%ngen_loops
                    call getf(init_trial_in%opt_data%cutoff_amps(I))
                end do
            case("OPTIMISED-INIT-CUTOFF-NUM")
                init_trial_in%opt_data%tAmpCutoff = .false.
                init_trial_in%opt_data%ngen_loops = nitems - 1
                allocate(init_trial_in%opt_data%cutoff_nums(init_trial_in%opt_data%ngen_loops))
                do I = 1, init_trial_in%opt_data%ngen_loops
                    call geti(init_trial_in%opt_data%cutoff_nums(I))
                end do
            case("HF-INIT")
                init_trial_in%tHF = .true.
            case("POPS-INIT")
                init_trial_in%tPops = .true.
                call geti(init_trial_in%npops)
            case("READ-INIT")
                init_trial_in%tRead = .true.
            case("FCI-INIT")
                init_trial_in%tFCI = .true.
                tStartSinglePart = .false.
            case("HEISENBERG-FCI-INIT")
                init_trial_in%tHeisenbergFCI = .true.
            case("START-FROM-HF")
                tStartCoreGroundState = .false.
            case("INC-CANCELLED-INIT-ENERGY")
!If true, include the spawnings cancelled due the the initiator criterion in the trial energy.
                tIncCancelledInitEnergy = .true.
            case("STEPSSHIFTIMAG")
!For FCIMC, this is the amount of imaginary time which will elapse between updates of the shift.
                call getf(StepsSftImag)
            case("STEPSSHIFT")
!For FCIMC, this is the number of steps taken before the Diag shift is updated
                if(tFixedN0 .or. tTrialShift)then
                    write(6,*) "WARNING: 'STEPSSHIFT' cannot be changed. &
                    & 'FIXED-N0' or 'TRIAL-SHIFT' is already specified and sets this parameter to 1."
                else
                    call geti(StepsSft)
                end if
            case("FIXED-N0")
#ifdef __CMPL
                call stop_all(t_r, 'FIXED-N0 currently not implemented for complex')
#endif
                tFixedN0 = .true.
                call geti(N0_Target)
                !In this mode, the shift should be updated every iteration.
                !Otherwise, the dynamics is biased.
                StepsSft = 1
                !Also avoid changing the reference determinant.
                tReadPopsChangeRef = .false.
                tChangeProjEDet = .false.
            case("TRIAL-SHIFT")
#ifdef __CMPL
                call stop_all(t_r, 'TRIAL-SHIFT currently not implemented for complex')
#endif
                if (item.lt.nitems) then
                    call getf(TrialTarget)
                end if
                tTrialShift = .true.
                StepsSft = 1
            case("ADAPTIVE-SHIFT")
                if (item.lt.nitems) then
                    call getf(AdaptiveShiftSigma)
                end if
                if (item.lt.nitems) then
                    call getf(AdaptiveShiftF1)
                    if(AdaptiveShiftF1<0.0 .or. AdaptiveShiftF1>1.0)then
                        call stop_all(t_r, 'AdaptiveShiftF1 is a scaling parameter and should be between 0.0 and 1.0')
                    end if
                end if
                if (item.lt.nitems) then
                    call getf(AdaptiveShiftF2)
                    if(AdaptiveShiftF2<0.0 .or. AdaptiveShiftF2>1.0)then
                        call stop_all(t_r, 'AdaptiveShiftF2 is a scaling parameter and should be between 0.0 and 1.0')
                    end if
                end if
                tAdaptiveShift = .true.
            case("AUTO-ADAPTIVE-SHIFT")
                tAutoAdaptiveShift = .true.
                tAdaptiveShift = .true.
                if (item.lt.nitems) then
                    call getf(AdaptiveShiftThresh)
                end if

                if (item.lt.nitems) then
                    call getf(AdaptiveShiftExpo)
                end if

                if (item.lt.nitems) then
                    call getf(AdaptiveShiftCut)
                end if
            case("AAS-MATELE")
                tAAS_MatEle = .true.
                !When using the MatEle, the default value of 10 becomes meaningless
                AdaptiveShiftThresh = 0.0
             case("INITS-PROJE")
                ! only take initiators into account for the projected energy
                tInitsEnergy = .true.
             case("INITS-GAMMA0")
                ! use the density matrix obtained from the initiator space to 
                ! correct for the adaptive shift
                tInitsRDMRef = .true.
            case("EXITWALKERS")
!For FCIMC, this is an exit criterion based on the total number of walkers in the system.
                call getiLong(iExitWalkers)

            case("TARGETGROWRATE")
                ! For FCIMC, this is the target growth rate once in vary shift mode.
                call getf(InputTargetGrowRate)
                call getiLong(InputTargetGrowRateWalk)

            case("READPOPS")
!For FCIMC, this indicates that the initial walker configuration will be read in from the file POPSFILE, which must be present.
!DiagSft and InitWalkers will be overwritten with the values in that file.
                TReadPops=.true.
                tStartSinglePart=.false.
                if (item.lt.nitems) then
                    call readi(iPopsFileNoRead)
                    iPopsFileNoWrite = iPopsFileNoRead
                    iPopsFileNoRead = -iPopsFileNoRead-1
                end if
            case("WALKERREADBATCH")
                !The number of walkers to read in on the head node in each batch during a popsread
                call readi(iReadWalkersRoot)
            case("POPSFILEMAPPING")
!This indicates that we will be mapping a popsfile from a smaller basis calculation, into a bigger basis calculation.
!Requires a "mapping" file.
                call stop_all(t_r,'POPSFILEMAPPING deprecated')
            case("READPOPSTHRESH")
!When reading in a popsfile, this will only save the determinant, if the number of particles on this 
!determinant is greater than iWeightPopRead.
                tReadPops=.true.
                call readf(iWeightPopRead)
                if (item.lt.nitems) then
                    call readi(iPopsFileNoRead)
                    iPopsFileNoWrite = iPopsFileNoRead
                    iPopsFileNoRead = -iPopsFileNoRead-1
                end if
            case("READPOPS-CHANGEREF")
                ! When reading in a pops file, use the most highly weighted
                ! determinant as the reference determinant for calculating
                ! the projected energy.
                ! Equivalent to PROJE-CHANGEREF at this point.
                tReadPopsChangeRef = .true.
            case("READPOPS-RESTARTNEWREFDET")
                ! When reading in a popsfile, restart the calculation
                ! according to the other parameters in the input file, but
                ! using the most highly weighted determinant as the reference
                ! determinant.
                tReadPopsRestart = .true.
            case("WALKCONTGROW")
!This option goes with the above READPOPS option.  If this is present - the INITWALKERS value is not 
!overwritten, and the walkers are continued to be allowed to grow before reaching                
!this value.  Without this keyword, when a popsfile is read in, the number of walkers is kept at the number 
!in the POPSFILE regardless of whether the shift had been allowed to change in the previous calc.
                tWalkContGrow=.true.
            case("SCALEWALKERS")
!For FCIMC, if this is a way to scale up the number of walkers, after having read in a POPSFILE
                call getf(ScaleWalkers)
            case("BINCANCEL")
!This is a seperate method to cancel down to find the residual walkers from a list, involving binning the walkers 
!into their determinants. This has to refer to the whole space, and so is much slower.
                TBinCancel=.true.
            case("REFSHIFT")
!With this, the shift is changed in order to keep the population on the reference determinant fixed, rather 
!than the total population.
                tShiftonHFPop=.true.
            case("STARTMP1")
!For FCIMC, this has an initial configuration of walkers which is proportional to the MP1 wavefunction
!                CALL Stop_All(t_r,"STARTMP1 option depreciated")
                TStartMP1=.true.
                TStartSinglePart=.false.
                if(item.lt.nitems) then
                    !Allow us to specify a desired number of particles to start with, so that the shift doesn't
                    !change dramatically to start with.
                    call getf(InitialPart)
                endif
            case("STARTCAS")
!For FCIMC, this has an initial configuration of walkers which is proportional to the MP1 wavefunction
!                CALL Stop_All(t_r,"STARTMP1 option depreciated")
                TStartCAS=.true.
                TStartSinglePart=.false.
                call geti(OccCASOrbs)  !Number of electrons in CAS 
                call geti(VirtCASOrbs)  !Number of virtual spin-orbitals in CAS
                if(item.lt.nitems) then
                    !Allow us to specify a desired number of particles to start with, so that the shift doesn't
                    !change dramatically to start with.
                    call getf(InitialPart)
                endif
            case("EQUILSTEPS")
!For FCIMC, this indicates the number of cycles which have to
!pass before the energy of the system from the doubles (HF)
!or singles (natural orbitals) population is counted.
                call geti(NEquilSteps)
            case("SHIFTEQUILSTEPS")
!This is the number of steps after the shift is allowed to change that it begins averaging the shift value.                
                call geti(NShiftEquilSteps)
            case("NOBIRTH")
!For FCIMC, this means that the off-diagonal matrix elements become zero, and so all we get is an exponential 
!decay of the initial populations on the determinants, at a rate which can be exactly calculated and compared against.
                CALL Stop_All(t_r,"NOBIRTH option depreciated")
!                TNoBirth=.true.
            case("MCDIFFUSE")
                CALL Stop_All(t_r,"MCDIFFUSE option depreciated")
!                TDiffuse=.true.
!Lambda indicates the amount of diffusion compared to spawning in the FCIMC algorithm.
!                call getf(Lambda)
            case("FLIPTAU")
!This indicates that time is to be reversed every FlipTauCyc cycles in the FCIMC algorithm. This might 
!help with undersampling problems.
                CALL Stop_All(t_r,"FLIPTAU option depreciated")
!                TFlipTau=.true.
!                call geti(FlipTauCyc)
            case("NON-PARTCONSDIFF")
!This is a seperate partitioning of the diffusion matrices in FCIMC in which the antidiffusion matrix (+ve connections) 
!create a net increase of two particles.
                CALL Stop_All(t_r,"NON-PARTCONSDIFF option depreciated")
!                TExtraPartDiff=.true.
            case("FULLUNBIASDIFF")
!This is for FCIMC, and fully unbiases for the diffusion process by summing over all connections
                CALL Stop_All(t_r,"FULLUNBIASDIFF option depreciated")
!                TFullUnbias=.true.
            case("NODALCUTOFF")
!This is for all types of FCIMC, and constrains a determinant to be of the same sign as the MP1 wavefunction at 
!that determinant, if the normalised component of the MP1 wavefunction is greater than the NodalCutoff value.
                CALL Stop_All(t_r,"NODALCUTOFF option depreciated")
!                TNodalCutoff=.true.
!                call getf(NodalCutoff)
            case("NOANNIHIL")
!For FCIMC, this removes the annihilation of particles on the same determinant step.
                TNoAnnihil=.true.
            case("MAXCHAINLENGTH")
!For closed path MC, this is the maximum allowed chain length before it is forced to come back
                call geti(CLMAX)
            case("RETURNBIAS")
!For closed path MC, this is the return bias at any point to spawn at the parent determinant
                call getf(PRet)
            case("RHOAPP")
!This is for resummed FCIMC, it indicates the number of propagation steps around each subgraph before 
!particles are assigned to the nodes
                CALL Stop_All(t_r,"RHOAPP option depreciated")
!                call geti(RhoApp)
            case("SIGNSHIFT")
!This is for FCIMC and involves calculating the change in shift depending on the absolute value of the 
!sum of the signs of the walkers.
!This should hopefully mean that annihilation is implicitly taken into account.
                TSignShift=.true.
            case("HFRETBIAS")
!This is a simple guiding function for FCIMC - if we are at a double excitation, then we return to the HF 
!determinant with a probability PRet.
!This is unbiased by the acceptance probability of returning to HF.
                THFRetBias=.true.
                call getf(PRet)
            case("EXCLUDERANDGUIDE")
!This is an alternative method to unbias for the HFRetBias. It invloves disallowing random 
!excitations back to the guiding function (HF Determinant)
                CALL Stop_All(t_r,"EXCLUDERANDGUIDE option depreciated")
!                TExcludeRandGuide=.true.
            case("PROJECTE-MP2")
!This will find the energy by projection of the configuration of walkers onto the MP2 wavefunction.
                TProjEMP2=.true.
            case("PROJE-CHANGEREF")

                ! If there is a determinant larger than the current reference,
                ! then swap references on the fly
                ! 
                ! The first parameter specifies the relative weights to trigger
                ! the change.

                ! The second parameter specifies an absolute minimum weight.

                tCheckHighestPop=.true.
                tChangeProjEDet=.true.
                IF(item.lt.nitems) then
                    call Getf(FracLargerDet)
                ENDIF
                if (item < nitems) then
                    call getf(pop_change_min)
                endif

            case("NO-CHANGEREF")

                ! Now that changing the reference determinant is default
                ! behaviour, we want a way to turn that off!

                tReadPopsChangeRef = .false.
                tChangeProjEDet = .false.
                
            case("AVGROWTHRATE")

                ! This option will average the growth rate over the update 
                ! cycle when updating the shift.

                if (item < nitems) then
                    call readu(w)
                    select case(w)
                    case("OFF")
                        tInstGrowthRate = .true.
                    case default
                        tInstGrowthRate = .false.
                    end select
                else
                    tInstGrowthRate = .false.
                end if

             case("L2-GROWRATE")
                ! use the L2-norm instead of the L1 norm to get the shift
                tL2GrowRate = .true.

            case("RESTARTLARGEPOP")
                tCheckHighestPop=.true.
                tRestartHighPop=.true.
                IF(item.lt.nitems) then
                    call Getf(FracLargerDet)
                ENDIF
                IF(item.lt.nitems) then
                    call Geti(iRestartWalkNum)
                ENDIF
            case("FIXPARTICLESIGN")
!This uses a modified hamiltonian, whereby all the positive off-diagonal hamiltonian matrix elements are zero. 
!Instead, their diagonals are modified to change the
!on-site death rate. Particles now have a fixed (positive) sign which cannot be changed and so no annihilation occurs.
                TFixParticleSign=.true.
            case("MAXBLOOMWARNONLY")
                !This means that we only get a particle bloom warning if the bloom is larger than any previous blooming event.
                tMaxBloom=.true.
            case("STARTSINGLEPART")
!A FCIMC option - this will start the simulation with a single positive particle at the HF, and fix the 
!shift at its initial value, until the number of particles gets to the INITPARTICLES value.
                TStartSinglePart=.true.
                IF(item.lt.nitems) THEN
                    !If an optional integer keyword is added, then InitialPart will indicate the number of
                    !particles to start at the HF determinant.
                    call readf(InitialPart)
                    if (InitialPart < 0) then
                        ! Turn StartSinglePart off.
                        tStartSinglePart = .false.
                        InitialPart = 1
                    end if
                ENDIF
            case("MEMORYFACPART")
!An FCIMC option - MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
                CALL Getf(MemoryFacPart)
            case("MEMORYFACANNIHIL")
!!An FCIMC option - MemoryFac is the factor by which space will be made available for particles sent to 
!the processor during annihilation compared to InitWalkers. This will generally want to be larger than
!memoryfacPart, because the parallel annihilation may not be exactly load-balanced because of differences 
!in the wavevector and uniformity of the hashing algorithm.
                call stop_all(t_r,'MEMORYFACANNIHIL should not be needed any more')
            case("MEMORYFACSPAWN")
!A parallel FCIMC option for use with ROTOANNIHILATION. This is the factor by which space will be made 
!available for spawned particles each iteration. 
!Several of these arrays are needed for the annihilation process. With ROTOANNIHILATION, MEMORYFACANNIHIL 
!is redundant, but MEMORYFACPART still need to be specified.
                CALL Getf(MemoryFacSpawn)
            case("MEMORYFACINIT")
                ! If we are maintaining a list of initiators on each
                ! processor, this is the factor of InitWalkers which will be
                ! used for the size
                call getf(MemoryFacInit)
            case("REGENEXCITGENS")
!An FCIMC option. With this, the excitation generators for the walkers will NOT be stored, and regenerated 
!each time. This will be slower, but save on memory.
                TRegenExcitGens=.true.

            case("REGENDIAGHELS")
                ! A parallel FCIMC option. With this, the diagonal elements of
                ! the hamiltonian matrix will not be stored with each particle.
                ! This will generally be slower, but save on memory.
                call stop_all(t_r, 'This option (REGENDIAGHELS) has been &
                                   &deprecated')

            case("FIXSHELLSHIFT")
!An FCIMC option. With this, the shift is fixed at a value given here, but only for excitations which are less than 
!<ShellFix>. This will almost definitly give the wrong answers for both the energy
!and the shift, but may be of use in equilibration steps to maintain particle density at low excitations, before writing 
!out the data and letting the shift change.
                CALL Stop_All(t_r,"FIXSHELLSHIFT option depreciated")
!                TFixShiftShell=.true.
!                CALL Geti(ShellFix)
!                CALL Getf(FixShift)
            case("FIXKIISHIFT")
!A Parallel FCIMC option. Similar to FixShellShift option, but will fix the shifts of the particles which have a diagonal
!matrix element Kii of less than the cutoff, FixedKiiCutOff.
                CALL Stop_All(t_r,"FIXKIISHIFT option depreciated")
!                TFixShiftKii=.true.
!                CALL Getf(FixedKiiCutoff)
!                CALL Getf(FixShift)
            
            case("FIXCASSHIFT")                
!A Parallel FCIMC option similar to the FixShellShift and FixShiftKii options.
!In this option, an active space is chosen containing a certain number of highest occupied spin orbitals (OccCASorbs) and
!lowest unoccupied spin orbitals (VirtCASorbs).  The shift is then fixed only for determinants 
!which have completely occupied spin orbitals for those lower in energy than the active space, 
!and completely unoccupied spin orbitals above the active space.  i.e. the electrons are only excited within the active space.  
                CALL Stop_All(t_r,"FIXKIISHIFT option depreciated")
!                TFixCASShift=.true.
!                call Geti(OccCASorbs)
!                call Geti(VirtCASorbs)
!                call Getf(FixShift)

            case("TRUNCATECAS")
!A Parallel FCIMC option. With this, the excitation space of the determinants will only include 
!the determinants accessible to the CAS
!space as specified here. 
!In this option, an active space is chosen containing a certain number of highest occupied spin orbitals (OccCASorbs) and
!lowest unoccupied spin orbitals (VirtCASorbs).  The determinant only for determinants 
!which have completely occupied spin orbitals for those lower in energy than the active space, 
!and completely unoccupied spin orbitals above the active space.  i.e. the electrons are only excited within the active space.  
                tTruncCAS=.true.
                call Geti(OccCASOrbs)
                call Geti(VirtCASOrbs)

            case("TRUNCINITIATOR")
!This option goes along with the above TRUNCATECAS option.  This means that walkers are allowed to spawn on 
!determinants outside the active space, however if this is done, they
!can only spawn back on to the determinant from which they came.  This is the star approximation from the CAS space. 
                tTruncInitiator=.true.
             case("AVSPAWN-INITIATORS")
! Create initiators based on the average spawn onto some determinant                
                tActivateLAS = .true.
                if(item < nitems) call getf(spawnSgnThresh)
                if(item < nitems) call geti(minInitSpawns)
                
             case("DELAY-DEATHS")
                ! have determinants persits for a number of iterations after their
                ! occupation reached 0
                tTimedDeaths = .true.
                tWriteUnocc = .true.
                ! read in said number of iterations
                if(item < nitems) call geti(lingerTime)

             case("REPLICA-GLOBAL-INITIATORS")
! with this option, all replicas will use the same initiator flag, which is then set 
! depending on the avereage population, else, the initiator flag is set for each replica
! using the population of that replica
                tGlobalInitFlag = .true.

             case("AVERAGE-REPLICAS")
                ! average the replica populations if they are not sign coherent
                tAVReps = .true.
                
             case("REPLICA-COHERENT-INITS")
                ! require initiators to be coherent across replcias
                tReplicaCoherentInits = .true.

            case("NO-COHERENT-INIT-RULE")
                tInitCoherentRule=.false.

             case("ALL-SENIORITY-INITS")
                ! make all determinants with at most initMaxSenior open orbitals initiators
                tSeniorityInits = .true.
                ! the maximum number of open orbs, default is 0
                if(item < nitems) call getI(initMaxSenior)

             case("ALL-SENIORITY-SURVIVE")
                ! keep all spawns, regardless of initiator criterium, onto 
                ! determinants up to a given Seniority level and excitation level
                tSpawnSeniorityBased = .true.
                do while(item < nitems) 
                   ! Default: Max Seniority level 0
                   call geti(maxKeepNOpenBuf)
                   ! Default: Max excit level 8
                   call geti(maxKeepExLvlBuf)
                   call addToIntArray(maxKeepExLvl,maxKeepNOpenBuf+1,maxKeepExLvlBuf)
                   numMaxExLvlsSet = maxKeepNOpenBuf+1
                end do

             case("LARGE-MATEL-SURVIVE")
                ! keep all spawns with a matrix element larger than a given threshold
                tLargeMatelSurvive = .true.
                call getf(spawnMatelThresh)

! Epstein-Nesbet second-order perturbation using the stochastic spawnings to correct initiator error.
            case("EN2-INITIATOR")
                tEN2 = .true.
                tEN2Init = .true.

! Epstein-Nesbet second-order perturbation using stochastic spawnings. However, this is not used to
! correct initiator error. Currently, it is only used for the full non-initiator scheme when applied
! to a truncated space. Then, an EN2 correction is applied to the space outside that truncation.
            case("EN2-TRUNCATED")
                tEN2 = .true.
                tEN2Truncated = .true.

            case("KEEPDOUBSPAWNS")
!This means that two sets of walkers spawned on the same determinant with the same sign will live, 
!whether they've come from inside or outside the CAS space.  Before, if both of these
!were from outside the space, they would've been aborted.
!                tKeepDoubleSpawns=.true.
!This option is now on permanently by default and cannot be turned off.

            case("ADDTOINITIATOR")
!This option means that if a determinant outside the initiator space becomes significantly populated - 
!it is essentially added to the initiator space and is allowed to spawn where it likes.
!The minimum walker population for a determinant to be added to the initiator space is InitiatorWalkNo.
                tAddtoInitiator=.true.
                call getf(InitiatorWalkNo)

            case("SENIOR-INITIATORS")
!This option means that if a determinant has lived  long enough (called a 'senior determinant'),
!it is added to the initiaor space. A determinant is considered 'senior' if its life time (measured in its halftime) exceeds SeniortyAge.
                tSeniorInitiators =.true.
                if(item.lt.nitems) then
                    call getf(SeniorityAge)
                end if
            case("INITIATOR-ENERGY-CUTOFF")
                !
                ! Specify both a threshold an an addtoinitiator value for
                ! varying the thresholds
                call stop_all(t_r,'Deprecated Option')

            case("SPAWNONLYINIT", "SPAWNONLYINITGROWTH")
                call stop_all(t_r, 'Option (SPAWNONLYINIT) deprecated')

            case("RETESTINITPOP")                
!This keyword is on by default.  It corresponds to the original initiator algorithm whereby a determinant may 
!be added to the initiator space if its population becomes higher 
!than InitiatorWalkNo (above), but if the pop then drops below this, the determinant is removed again from the initiator space.
!Having this on means the population is tested at every iteration, turning it off means that once a determinant 
!becomes an initiator by virtue of its population, it remains an initiator 
!for the rest of the simulation.
 
            case("INCLDOUBSINITIATOR")
!This keyword includes any doubly excited determinant in the 'initiator' space so that it may spawn as usual
!without any restrictions.
                tInitIncDoubs=.true.

            case("UNBIASPGENINPROJE")
!A FCIMC serial option. With this, walkers will be accepted with probability tau*hij. i.e. they will not unbias 
!for PGen in the acceptance criteria, but in the term for the projected energy.
                TUnbiasPGeninProjE=.true.
            case("ANNIHILATEONPROCS")
!A parallel FCIMC option. With this, walkers will only be annihilated with other walkers on the same processor. 
                CALL Stop_All(t_r,"ANNIHILATEONPROCS option depreciated")
!                TAnnihilonproc=.true.
            case("ANNIHILATDISTANCE")
!A Serial FCIMC experimental option. With this, walkers have the ability to annihilate each other as long as 
!they are connected, which they will do with probability = Lambda*Hij
                CALL Stop_All(t_r,"ANNIHILATEONPROCS option depreciated")
!                TDistAnnihil=.true.
!                call Getf(Lambda)
            case("ANNIHILATERANGE")
!This option should give identical results whether on or off. It means that hashes are histogrammed and sent 
!to processors, rather than sent due to the value of mod(hash,nprocs).
!This removes the need for a second full sorting of the list of hashes, but may have load-balancing issues for the algorithm.
!This now is on by default, and can only be turned off by specifying OFF after the input.
                CALL Stop_All(t_r,"ANNIHILATEONPROCS option depreciated")
!                IF(item.lt.nitems) then
!                    call readu(w)
!                    select case(w)
!                    case("OFF")
!                        tAnnihilatebyrange=.false.
!                    end select
!                ELSE
!                    tAnnihilatebyrange=.true.
!                ENDIF
            case("ROTOANNIHILATION")
!A parallel FCIMC option which is a different - and hopefully better scaling - algorithm. This is substantially 
!different to previously. It should involve much less memory.
!MEMORYFACANNIHIL is no longer needed (MEMORYFACPART still is), and you will need to specift a MEMORYFACSPAWN 
!since newly spawned walkers are held on a different array each iteration.
!Since the newly-spawned particles are annihilated initially among themselves, you can still specift 
!ANNIHILATEATRANGE as a keyword, which will change things.
                CALL Stop_All(t_r,"ROTOANNIHILATION option depreciated")
!                tRotoAnnihil=.true.
            case("DIRECTANNIHILATION")
!A parallel FCIMC option which is a different annihilation algorithm. It has elements in common with both 
!rotoannihilation and the hashing annihilation, but hopefully will be quicker and
!better scaling with number of processors. It has no explicit loop over processors.
                tDirectAnnihil=.true.
            case("LOCALANNIHIL")
!A parallel FCIMC experimental option. This will attempt to compensate for undersampled systems, by 
!including extra annihilation for walkers which are the sole occupier of determiants
!This annihilation is governed by the parameter Lambda, which is also used in other circumstances 
!as a variable, but should not be used at the same time.
                CALL Stop_All(t_r,"LOCALANNIHIL option depreciated")
!                TLocalAnnihilation=.true.
!                call Getf(Lambda)
            case("ANNIHILATEEVERY")
!In FCIMC, this will result in annihilation only every iAnnInterval iterations
                call Geti(iAnnInterval)
            case("GLOBALSHIFT")
                ! Parallel FCIMC option which has been removed.
                call stop_all (t_r, "GLOBALSHIFT - option removed")

            case("RANDOMISEHASHORBS")
                ! This will create a random 1-to-1 mapping between the 
                ! orbitals, which should hopefully improve load balancing.
                ! (now on always - sds)
                call stop_all (t_r, "RANDOMISEHASHORBS - option removed &
                                    &(now default)")

            case("SPATIAL-ONLY-HASH")
                ! Base hash values only on spatial orbitals
                ! --> All determinants with the same spatial structure will
                !     end up on the same processor
                tSpatialOnlyHash = .true.

             case("STORE-DETS")
                ! store all determinants in their decoded form in memory
                ! this gives a speed-up at the cost of the memory required for storing 
                ! all of them
                tStoredDets = .true.

            case("SPAWNASDETS")
!This is a parallel FCIMC option, which means that the particles at the same determinant on each processor, 
!will choose the same determinant to attempt spawning to and the 
!probability of a successful spawn will be multiplied by the number of particles on the determinant.
                tSpawnAsDet=.true.
            case("MAGNETIZE")
!This is a parallel FCIMC option. It chooses the largest weighted MP1 components and records their sign. 
!If then a particle occupies this determinant and is of the opposite sign, it energy,
!i.e. diagonal matrix element is raised by an energy given by BField.
                CALL Stop_All(t_r,"MAGNETIZE option depreciated")
!                tMagnetize=.true.
!                tSymmetricField=.false.
!                call Geti(NoMagDets)
!                call Getf(BField)

            case("FINDGROUNDDET")
                call stop_all(t_r, 'Option (FINDGROUNDDET) deprecated')

            case("STARORBS")
!A parallel FCIMC option. Star orbs means that determinants which contain these orbitals can only be spawned 
!at from the HF determinant, and conversly, can only spawn back at the HF determinant.
                CALL Stop_All(t_r,"STARORBS option depreciated")
!                call geti(iStarOrbs)
!                if(item.lt.nitems) then
!                    call readu(w)
!                    select case(w)
!                    case("NORETURN")
!!This option will mean that particles spawned at these high energy determinants will not be allowed to 
!spawn back at HF, but will be left to die.
!                        tNoReturnStarDets=.true.
!                    case("ALLSPAWNSTARDETS")
!!This option will mean that all particles can spawn at the star determinants and annihilation will take place 
!there. Once there however, they are left to die, and cannot spawn anywhere else.
!                        tAllSpawnStarDets=.true.
!                    end select
!                else
!                    tNoReturnStarDets=.false.
!                endif
!                tStarOrbs=.true.
            case("EXCITETRUNCSING")
!This is a parallel FCIMC option, where excitations between determinants where at least one of the determinants 
!is above iHighExcitsSing will be restricted to be single excitations.
                CALL Stop_All(t_r,"EXCITETRUNCSING option depreciated")
!                tHighExcitsSing=.true.
!                call readi(iHighExcitsSing)
            case("MAGNETIZESYM")
!A parallel FCIMC option. Similar to the MAGNETIZE option, but in addition to the energy being raised for 
!particles of the opposite sign, the energy is lowered by the same amount for particles
!of 'parallel' sign.
                CALL Stop_All(t_r,"MAGNETIZESYM option depreciated")
!                call Geti(NoMagDets)
!                call Getf(BField)
!                tSymmetricField=.true.
!                tMagnetize=.true.
            case("SINGLESBIAS")
!This is a parallel FCIMC option, where the single excitations from any determinant will be favoured compared 
!to the simple ratio of number of doubles to singles from HF by multiplying the number of singles by this factor.
                call Getf(SinglesBias)
            case("JUSTFINDDETS")
!This option is to be used in conjunction with the diagonalization methods. With this, all the determinants 
!will be enumerated, but the hamiltonian will not be calculated,
!and the energies not calculated. This is needed when the full list of determinants is needed for later on.
                tFindDets=.true.
            case("EXPANDSPACE")
                call report(" "//trim(w)//" is a depreciated option - look at EXPANDFULLSPACE",.true.)
            case("EXPANDFULLSPACE")
!Read in a value of the iteration to expand to the full space.
                call geti(iFullSpaceIter)
            case("MULTIPLEDETSSPAWN")
!This option creates connections from iDetGroup randomly chosen determinants and attempts to spawn from them 
!all at once. This should hopefully mean that annihilations are implicitly done.
                CALL Stop_All(t_r,"MULTIPLEDETSSPAWN option depreciated")
!                tMultipleDetsSpawn=.true.
!                call Geti(iDetGroup)

            case("SPIN-PROJECT")
                ! Enable spin projection (spin_project.F90).
                ! Optional argument specifies no. of iterations between
                ! each application of stochastic spin projection.
                tSpinProject = .true.
                if (item < nitems) call geti (spin_proj_interval)

            case("SPIN-PROJECT-GAMMA")
                ! Change the value of delta-gamma used by the spin projection
                ! routines. Similar to modifying tau for normal FCIQMC.
                call getf (spin_proj_gamma)

            case("SPIN-PROJECT-SHIFT")
                ! Change the value of delta-gamma used by the spin projection
                ! routines. Similar to modifying tau for normal FCIQMC.
                call getf (spin_proj_shift)

            case("SPIN-PROJECT-CUTOFF")
                ! Change the minimum number of walkers required for spin
                ! projection to be applied to a determinant
                call geti (spin_proj_cutoff)

            case("SPIN-PROJECT-STOCHASTIC-YAMA")
                ! Only project via one Yamanouchi symbol on each iteration, 
                ! selecting that symbol stochastically.
                spin_proj_stochastic_yama = .true.

            case("SPIN-PROJECT-NOPEN-LIMIT")
                ! Determine the largest number of unpaired electrons a
                ! determinant may have for us to apply spin projectino to it.
                !
                ! --> Attempt to reduce the exponential scaling of the
                !     projection sum.
                call geti (spin_proj_nopen_max)

            case("SPIN-PROJECT-SPAWN-INITIATORS")
                ! If TRUNCINITIATOR is set, then ensure that all children of
                ! initiators created by spin projection are made into
                ! initiators.
                spin_proj_spawn_initiators = readt_default (.true.)
                if (spin_proj_spawn_initiators) &
                    write(6,*) 'Disabling spin projected progeny of &
                               &initiators automatically being initiators'

            case("SPIN-PROJECT-NO-DEATH")
                ! Only spawn, don't die, particles in spin projection
                spin_proj_no_death = readt_default (.true.)
                if (spin_proj_no_death) &
                    write(6,*) 'Disabling death for spin projection'

            case("SPIN-PROJECT-ITER-COUNT")
                ! How many times should the spin projection step be applied 
                ! on each occasion it gets called? (default 1)
                call geti (spin_proj_iter_count)

            case("SPIN-PROJECT-VARYSHIFT-OFF")
                ! When VARYSHIFT is enabled, turn spin projection off.
                ! TODO: Should this be made default?
                if (item < nitems) then
                    call readu(w)
                    select case(w)
                    case("OFF")
                        disable_spin_proj_varyshift = .false.
                    case default
                        disable_spin_proj_varyshift = .true.
                    end select
                else
                    disable_spin_proj_varyshift = .true.
                endif

            case("TRUNC-NOPEN")
                ! Truncate determinant spawning at a specified number of
                ! unpaired electrons.
                tTruncNOpen = .true.
                call geti (trunc_nopen_max)

            case("WEAKINITIATORS")
                !Additionally allow the children of initiators to spawn freely
                !This adaptation is applied stochastically with probability weakthresh
                !Hence weakthresh = 1 --> Always on where applicable.
                !weakthresh = 0 --> The original initiator scheme is maintained.
                call stop_all(t_r,'Deprecated option')

            case("ALLREALCOEFF")
                tAllRealCoeff=.true.
                tUseRealCoeffs = .true.
                !Turn on continuous spawning/death
                !Kill populations n<1 with probability 1-n
            case("REALCOEFFBYEXCITLEVEL")
                tRealCoeffByExcitLevel=.true.
                tUseRealCoeffs = .true.
                call readi(RealCoeffExcitThresh)
            case("KEEPWALKSMALL")
                call stop_all(t_r,'Deprecated Option')
            case("REALSPAWNCUTOFF")
                tRealSpawnCutoff=.true.
                call Getf(RealSpawnCutoff)
            case("SETOCCUPIEDTHRESH")
                call Getf(OccupiedThresh)
            case("SETINITOCCUPIEDTHRESH")
                call stop_all(t_r,'Deprecated option')

            case("JUMP-SHIFT")
                ! When variable shift is enabled, jump the shift to the value
                ! predicted by the projected energy!
                ! --> Reduce the waiting time while the number of particles is
                !     growing.
                !
                ! This is now the default behaviour. Use JUMP-SHIFT OFF to
                ! disable it (likely only useful in some of the tests).
                tJumpShift = .true.
                if (item < nitems) then
                    call readu(w)
                    select case(w)
                    case("OFF", "FALSE")
                        tJumpShift = .false.
                    case default
                    end select
                end if

            case("UNIQUE-HF-NODE")
                ! Assign the HF processor to a unique node.
                ! TODO: Set a default cutoff criterion for this
                tUniqueHFNode = .true.

            case("LET-INIT-POP-DIE")
                tLetInitialPopDie = .true.

            case("POPS-ANNIHILATE")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npops_pert)
                if (.not. allocated(pops_pert)) then
                    allocate(pops_pert(npops_pert))
                else
                    if (npops_pert /= size(pops_pert)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npops_pert
                    call read_line(eof)
                    pops_pert(i)%nannihilate = nitems
                    allocate(pops_pert(i)%ann_orbs(nitems))
                    do j = 1, nitems
                        call readi(pops_pert(i)%ann_orbs(j))
                    end do
                    ! Create the rest of the annihilation-related
                    ! components of the pops_pert object.
                    call init_perturbation_annihilation(pops_pert(i))
                end do

            case("POPS-CREATION")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npops_pert)
                if (.not. allocated(pops_pert)) then
                    allocate(pops_pert(npops_pert))
                else
                    if (npops_pert /= size(pops_pert)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npops_pert
                    call read_line(eof)
                    pops_pert(i)%ncreate = nitems
                    allocate(pops_pert(i)%crtn_orbs(nitems))
                    do j = 1, nitems
                        call readi(pops_pert(i)%crtn_orbs(j))
                    end do
                    ! Create the rest of the creation-related
                    ! components of the pops_pert object.
                    call init_perturbation_creation(pops_pert(i))
                end do

            case("WRITE-POPS-NORM")
                tWritePopsNorm = .true.

            ! Options relating to finite-temperature Lanczos calculations.
            case("NUM-INIT-VECS-FTLM")
                call geti(n_init_vecs_ftlm)
            case("NUM-LANC-VECS-FTLM")
                call geti(n_lanc_vecs_ftlm)
            case("NUM-BETA-FTLM")
                call geti(nbeta_ftlm)
            case("BETA-FTLM")
                call getf(delta_beta_ftlm)

            ! Options relating to exact spectral calculations.
            case("NUM-LANC-VECS-SPECTRAL")
                call geti(n_lanc_vecs_sl)
            case("NUM-OMEGA-SPECTRAL")
                call geti(nomega_spectral)
            case("OMEGA-SPECTRAL")
                call getf(delta_omega_spectral)
            case("MIN-OMEGA-SPECTRAL")
                call getf(min_omega_spectral)
            case("BROADENING_SPECTRAL")
                call getf(spectral_broadening)
            case("INCLUDE-GROUND-SPECTRAL")
                tIncludeGroundSpectral = .true.
            case("GROUND-ENERGY-SPECTRAL")
                call getf(spectral_ground_energy)

            case("LEFT-ANNIHILATE-SPECTRAL")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npert_spectral_left)
                if (.not. allocated(left_perturb_spectral)) then
                    allocate(left_perturb_spectral(npert_spectral_left))
                else
                    if (npert_spectral_left /= size(left_perturb_spectral)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert_spectral_left
                    call read_line(eof)
                    left_perturb_spectral(i)%nannihilate = nitems
                    allocate(left_perturb_spectral(i)%ann_orbs(nitems))
                    do j = 1, nitems
                        call readi(left_perturb_spectral(i)%ann_orbs(j))
                    end do
                    ! Create the rest of the annihilation-related
                    ! components of the left_perturb_spectral object.
                    call init_perturbation_annihilation(left_perturb_spectral(i))
                end do
            case("LEFT-CREATION-SPECTRAL")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npert_spectral_left)
                if (.not. allocated(left_perturb_spectral)) then
                    allocate(left_perturb_spectral(npert_spectral_left))
                else
                    if (npert_spectral_left /= size(left_perturb_spectral)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert_spectral_left
                    call read_line(eof)
                    left_perturb_spectral(i)%ncreate = nitems
                    allocate(left_perturb_spectral(i)%crtn_orbs(nitems))
                    do j = 1, nitems
                        call readi(left_perturb_spectral(i)%crtn_orbs(j))
                    end do
                    ! Create the rest of the creation-related
                    ! components of the left_perturb_spectral object.
                    call init_perturbation_creation(left_perturb_spectral(i))
                end do

            case("RIGHT-ANNIHILATE-SPECTRAL")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npert_spectral_right)
                if (.not. allocated(right_perturb_spectral)) then
                    allocate(right_perturb_spectral(npert_spectral_right))
                else
                    if (npert_spectral_right /= size(right_perturb_spectral)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert_spectral_right
                    call read_line(eof)
                    right_perturb_spectral(i)%nannihilate = nitems
                    allocate(right_perturb_spectral(i)%ann_orbs(nitems))
                    do j = 1, nitems
                        call readi(right_perturb_spectral(i)%ann_orbs(j))
                    end do
                    ! Create the rest of the annihilation-related
                    ! components of the right_perturb_spectral object.
                    call init_perturbation_annihilation(right_perturb_spectral(i))
                end do
            case("RIGHT-CREATION-SPECTRAL")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npert_spectral_right)
                if (.not. allocated(right_perturb_spectral)) then
                    allocate(right_perturb_spectral(npert_spectral_right))
                else
                    if (npert_spectral_right /= size(right_perturb_spectral)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert_spectral_right
                    call read_line(eof)
                    right_perturb_spectral(i)%ncreate = nitems
                    allocate(right_perturb_spectral(i)%crtn_orbs(nitems))
                    do j = 1, nitems
                        call readi(right_perturb_spectral(i)%crtn_orbs(j))
                    end do
                    ! Create the rest of the creation-related
                    ! components of the right_perturb_spectral object.
                    call init_perturbation_creation(right_perturb_spectral(i))
                end do

            case("TAU-CNT-THRESHOLD")
                write(6,*) 'WARNING: This option is unused in this branch'

            case("INITIATOR-SURVIVAL-CRITERION")
                ! If a site survives for at least a certain number of
                ! iterations, it should be treated as an initiator.
                ! --> Soft expand the range of the initiators in the Hilbert
                !     space
                call stop_all(t_r,'Deprecated option')

            case("INITIATOR-SURVIVAL-MULTIPLIER")
                ! If a site survives for a certain multiple of how long it
                ! would _expect_ to have survived, then it should be treated
                ! as an initiator
                ! --> A more flexible version of INITIATOR-SURVIVAL-CRITERION
                call stop_all(t_r,'Deprecated option')

            case("INITIATOR-SPAWN-CRITERION")
                ! A site becomes an initiator once a certain number of
                ! spawns have occurred to it (these must be independent
                ! spawns, rather than a certain magnitude of spawning)
                call stop_all(t_r,'Deprecated option')

            case("MULTI-REPLICA-INITIATORS")
                ! Aggregate particle counts across all of the simulation
                ! replicas to determine which sites are considered to be
                ! initiators.
                ! Obviously, this only does anything with system-replicas
                ! set...
                call stop_all(t_r,'Option Deprecated')

            case("ORTHOGONALISE-REPLICAS")
                ! Apply Gram Schmidt ortgogonalisation to replicas, starting
                ! with replica 1, so that we will collect excited states of
                ! a given symmetry
                tOrthogonaliseReplicas = .true.
                if (item < nitems) then
                    call readi(orthogonalise_iter)
                endif

                ! Don't start all replicas from the deterministic ground state
                ! when using this option.
                tStartCoreGroundState = .false.

            case("ORTHOGONALISE-REPLICAS-SYMMETRIC")
                ! Use the Lowdin (symmetric) orthogonaliser instead of the 
                ! Gram Schmidt one from the ORTHOGONALISE-REPLICAS option
                tOrthogonaliseReplicas = .true.
                tOrthogonaliseSymmetric = .true.
                if (item < nitems) then
                    call readi(orthogonalise_iter)
                end if

                ! Don't start all replicas from the deterministic ground state
                ! when using this option.
                tStartCoreGroundState = .false.

            case("REPLICA-SINGLE-DET-START")
                ! If we want to start off multiple replicas from single dets
                ! chosen fairly naively as excited states of the HF, then use
                ! this option
                tReplicaSingleDetStart = .true.

            case("DONT-PRINT-OVERLAPS")
                ! Don't print overlaps between replicas when using the
                ! orthogonalise-replicas option.
                tPrintReplicaOverlaps = .false.

            case("USE-SPAWN-HASH-TABLE")
                use_spawn_hash_table = .true.

            case("CONT-TIME-FULL")
                ! Use the full continuous time scheme, not the approximated
                ! oversampled scheme
                ! --> Needs to calculate the spawning rate for each det as it
                !     appears, so is slow
                tContTimeFull = .true.

            case("CONT-TIME-MAX-OVERSPAWN")
                ! Efficient continuous time propagation requires a fine
                ! interplay between the oversampling rate, and the maximum
                ! spawn allowed
                call readf(cont_time_max_overspawn)

            case("POSITIVE-HF-SIGN")
                tPositiveHFSign = .true.

            case("LOAD-BALANCE-BLOCKS")
                ! When load balancing, have more blocks to distribute than
                ! there are processors to do the redistributing. This allows
                ! us to shuffle walkers around in the system
                tLoadBalanceBlocks = .true.
                if (item < nitems) then
                    call readu(w)
                    select case(w)
                    case("OFF", "NO", "DISABLE")
                        tLoadBalanceBlocks = .false.
                    case default
                        tLoadBalanceBlocks = .true.
                    end select

                    if (tLoadBalanceBlocks) then
                        write(iout, '("WARNING: LOAD-BALANCE-BLOCKS option is &
                                    &now enabled by default.")')
                    end if
                end if
            
            case("POPS-JUMP-SHIFT")
                ! Use the same logic as JUMP-SHIFT, but reset the shift value
                ! after restarting with a POPSFILE
                !
                ! --> This prevents undesirable behaviour if the simulation is
                !     restarted with a different FCIDUMP file (i.e. during
                !     CASSCF calculations).
                tPopsJumpShift = .true.

            case("MULTI-REF-SHIFT")
                tMultiRefShift = .true.

            case("MP2-FIXED-NODE")
                call stop_all(t_r,'Deprecated option')

            case("INTERPOLATE-INITIATOR")
                ! Implement interpolation between aborting particles
                ! due to the initiator criterion, and accepting them, based
                ! on the ratio of the parents coefficient and the value of
                ! InitiatorWalkNo
                !
                ! This modifies the acceptance criterion such that
                !
                ! alpha = alpha_min + (((n_parent - OccupiedThresh) / (InitiatorWalkNo - OccupiedThresh)) ** gamma) * (alpha_max - alpha_min)
                !
                ! Additional optional parameters (with default):
                !
                ! i)   alpha_min (0.0)
                ! ii)  alpha_max (1.0)
                ! iii) gamma     (1.0)
                call stop_all(t_r,'Deprecated option')
            case("SHIFT-PROJECT-GROWTH")
                ! Extrapolate the expected number of walkers at the end of the
                ! _next_ update cycle for calculating the shift. i.e. use
                !
                ! log((N_t + (N_t - N_(t-1))) / N_t)
                call stop_all(t_r,'Option deprecated')

            case ("BACK-SPAWN")
                ! Alis idea to increase the chance of non-initiators to spawn
                ! to occupied determinants
                ! and out of laziness this is only introduced for 
                ! 4ind-weighted-2 and above excitation generators! 
                ! maybe for hubbard model too, but lets see..
                t_back_spawn = .true.
                t_back_spawn_option = .true.

                if (item < nitems) then 
                    t_back_spawn = .false.
                    call geti(back_spawn_delay)
                end if

            case ("BACK-SPAWN-OCC-VIRT")
                t_back_spawn = .true.
                t_back_spawn_occ_virt = .true.
                
                t_back_spawn_option = .true.

                if (item < nitems) then 
                    t_back_spawn = .false.

                    call geti(back_spawn_delay)
                end if

            case("BACK-SPAWN-FLEX")
                t_back_spawn_flex = .true. 
                t_back_spawn_flex_option = .true.

                if (item < nitems) then 
                    t_back_spawn_flex = .false.

                    call geti(back_spawn_delay)
                end if

                ! can be value: -1, 0(default), 1, 2) 
                ! to indicate (de-)excitation
                if (item < nitems) then 
                    call geti(occ_virt_level)
                end if

             case("ALL-DOUBS-INITIATORS")
                ! Set all doubles to be treated as initiators
                ! If truncinitiator is not set, this does nothing
                tAllDoubsInitiators = .true.   
                ! If given, take the number of references for doubles
                if(item < nitems) call geti(nRefsDoubs)

             case("ALL-DOUBS-INITIATORS-DELAY")
                ! Only start after this number of steps in variable shift mode with 
                ! the all-doubs-initiators
                if(item < nitems) call geti(allDoubsInitsDelay)
                tSetDelayAllDoubsInits = .true.
                tSetDelayAllSingsInits = .true.

             case("ALL-SINGS-INITIATORS")
                ! Make the singles of given references initiators
                tAllSingsInitiators = .true.
                ! If given, take the number of references for singles
                if(item < nitems) call geti(nRefsSings)
                
             case("READ-REFERENCES")
                ! Instead of generating new references, read in existing ones
                tReadRefs = .true.
                
             case("EXCITATION-PRODUCT-REFERENCES")
                ! Also add all excitation products of references to the reference space
                tProductReferences = .true.

             case("COHERENT-REFERENCES")
                ! Only make those doubles/singles initiators that are sign coherent
                ! with their reference(s)
                if(item < nitems) then
                   call readu(w)
                   select case(w)
                   case("STRICT")
                      tStrictCoherentDoubles = .true.
                   case("WEAK")
                      ! This is recommended, we first check if there is a sign 
                      ! tendency and then if it agrees with the sign on the det
                      tAvCoherentDoubles = .true.
                      tWeakCoherentDoubles = .true.
                   case("XI")
                      ! This is a minimalistic version that should not
                      ! be used unless you know what you're doing
                      tWeakCoherentDoubles = .true.
                   case("AV")
                      ! only using av ignores sign tendency and can overestimate
                      ! the correctness of a sign
                      tAvCoherentDoubles = .true.
		   case("OFF")
		      ! do not perform a coherence check
		      tAvCoherentDoubles = .false.
		      tWeakCoherentDoubles = .false.
                   case default
                      ! default is WEAK
                      tAvCoherentDoubles = .true.
                      tWeakCoherentDoubles = .true.
                   end select
                else
                   tWeakCoherentDoubles = .true.
                   tAvCoherentDoubles = .true.
                endif                   

             case("SECONDARY-SUPERINITIATORS")
                ! Enable superinitiators by coherence criteria
                superInitiatorLevel = 1
                ! As the secondary SIs are now self-consistently determined, it is
                ! highly unlikely that a level above 1 will do anything more
                if(item < nItems) call readi(superInitiatorLevel)
                
             case("DYNAMIC-SUPERINITIATORS")
                ! Re-evaluate the superinitiators every SIUpdateInterval steps
                ! Beware, this can be very expensive
		! By default, it is 100, to turn it off, use 0
                call readi(SIUpdateInterval)
		
       	     case("STATIC-SUPERINITIATORS")
	        ! Do not re-evaluate the superinitiators
		SIUpdateInterval = 0

             case("INITIATOR-COHERENCE-THRESHOLD")
                ! Set the minimal coherence parameter for superinitiator-related
                ! initiators
                call readf(coherenceThreshold)

             case("SUPERINITIATOR-COHERENCE-THRESHOLD")
                ! set the minimal coherence parameter for superinitiators
                call readf(SIThreshold)

             case("MIN-SI-CONNECTIONS")
                ! set the minimal number of connections with superinititators for 
                ! superinitiators-related initiators
                call readi(minSIConnect)
                ! optionally, allow to weight the connections with the population
                if(item < nItems) then
                   call readu(w)
                   select case(w)
                   case("WEIGHTED")
                      tWeightedConnections = .true.
                   case("UNWEIGHTED")
                      tWeightedConnections = .false.
                   case default
                      tWeightedConnections = .false.
                   end select
                endif

             case("SIGNED-REPLICA-AVERAGE")
                tSignedRepAv = .true.
                if(item < nitems) then
                   call readu(w)
                   select case(w)
                   case("OFF")
                      tSignedRepAv = .false.
                   case default
                      tSignedRepAv = .true.
                   end select
                endif

             case("ENERGY-SCALED-WALKERS")
                ! the amplitude unit of a walker shall be scaled with energy
                tEScaleWalkers = .true.
                ! and the number of spawns shall be logged
                tLogNumSpawns = .true.
                sfTag = 0
                if(item < nItems) then
                   call readu(w)
                   select case(w)
                   case("EXPONENTIAL")
                      sfTag = 1
                      sFAlpha = 0.1
                   case("POWER")
                      sfTag = 0
                   case("EXP-BOUND")
                      sfTag = 3
                      sFAlpha = 0.1
                      sFBeta = 0.01
                   case("NEGATIVE")
                      sfTag = 2
                   case default
                      sfTag = 0
                      call stop_all(t_r, "Invalid argument 1 of ENERGY-SCALED-WALKERS")
                   end select
                endif
                if(item < nitems) &
                     ! an optional prefactor for scaling 
                   call readf(sFAlpha)
                if(item < nitems) &
                     ! an optional exponent for scaling
                     call readf(sFBeta)
                ! set the cutoff to the minimal value
                RealSpawnCutoff = sFBeta

             case("SUPERINITIATOR-POPULATION-THRESHOLD")
                ! set the minimum value for superinitiator population
                call readf(NoTypeN)

	     case("SUPPRESS-SUPERINITIATOR-OUTPUT")	
	        ! just for backwards-compatibility

             case("WRITE-SUPERINITIATOR-OUTPUT")
                ! Do not output the newly generated superinitiators upon generation
                tSuppressSIOutput = .false.
                
             case("TARGET-REFERENCE-POP")
                tVariableNRef = .true.
                if(item < nItems) call readi(targetRefPop)

            case default
                call report("Keyword "                                &
     &            //trim(w)//" not recognized in CALC block",.true.)
            end select

          end do calc

          IF((.not.TReadPops).and.(ScaleWalkers.ne.1.0_dp)) THEN
              call report("Can only specify to scale walkers if READPOPS is set",.true.)
          ENDIF

          ! Set if we need virtual orbitals  (usually set).  Will be unset (by
          ! Calc readinput) if I_VMAX=1 and TENERGY is false
          if(.not.tEnergy.and.I_VMAX.eq.1)  tNeedsVirts=.false.

          ! If the max vertex level is 2 or less, then we just need to calculate
          ! <ij|ab> and never need <ib|aj> for double excitations.  We do need
          ! them if we're doing a complete diagonalisation.
          gen2CPMDInts=MAXVAL(NWHTAY(3,:)).ge.3.or.TEnergy


        END SUBROUTINE CalcReadInput



        Subroutine CalcInit()
          use constants, only: dp
          use SystemData, only: G1, Alat, Beta, BRR, ECore, LMS, nBasis, nBasisMax, STot,tCSF,nMsh,nEl
          use SystemData, only: tUEG,nOccAlpha,nOccBeta,ElecPairs,tExactSizeSpace,tMCSizeSpace,MaxABPairs,tMCSizeTruncSpace
          use IntegralsData, only: FCK, CST, nMax, UMat
          use IntegralsData, only: HFEDelta, HFMix, NHFIt, tHFCalc
          Use Determinants, only: FDet, tSpecDet, SpecDet, get_helement
          Use DetCalc, only: DetInv, nDet, tRead
          Use DetCalcData, only:  ICILevel
          use hilbert_space_size, only: FindSymSizeofSpace, FindSymSizeofTruncSpace 
          use hilbert_space_size, only: FindSymMCSizeofSpace, FindSymMCSizeExcitLevel
          use global_utilities
          
          real(dp) CalcT, CalcT2, GetRhoEps
          
          
          INTEGER I, IC,J, norb
          INTEGER nList
          HElement_t(dp) HDiagTemp
          character(*), parameter :: this_routine='CalcInit'

          Allocate(MCDet(nEl))
          call LogMemAlloc('MCDet',nEl,4,this_routine,tagMCDet)

          IF(NPATHS.EQ.-1) THEN
             WRITE(6,*) 'NPATHS=-1.  SETTING NPATHS to NDET'
             NPATHS=NDET
          ENDIF
          IF(NDET.GT.0.AND.ABS(DETINV)+NPATHS.GT.NDET) THEN
             WRITE(6,*) 'DETINV+NPATHS=',ABS(DETINV)+NPATHS,'>NDET=',NDET
             WRITE(6,*) 'Setting DETINV and NPATHS to 0'
             DETINV=0
             NPATHS=0
          ENDIF

          IF(THFCALC) THEN
             WRITE(6,*) "Calculating Hartree-Fock Basis"
             WRITE(6,*) "Max Iterations:",NHFIT
             WRITE(6,*) "FMIX,EDELTA",HFMIX,HFEDELTA
          ENDIF
          IF(TMONTE) THEN 
             WRITE(6,*) 'MC Determinant Symmetry:'
             WRITE(6,*) (MDK(I),I=1,4)
          ENDIF
! Thus would appear to be obsolete

!          IF(G_VMC_FAC.LE.0) THEN
!             WRITE(6,*) "G_VMC_FAC=",G_VMC_FAC
!             call stop_all(this_routine, "G_VNC_FAC LE 0")
!          ENDIF

          IF(BETAP.NE.0.0_dp) THEN 
             I_P=NINT(BETA/BETAP)
             IF(.not.tFCIMC) THEN
                 WRITE(6,*) 'BETAP=',BETAP
                 WRITE(6,*) 'RESETTING P '
                 IF(I_P.GT.100000) WRITE(6,*) '*** WARNING I_P=',I_P
             ENDIF
          ENDIF

          IF(.not.tFCIMC) WRITE(6,*) 'BETA, P :',BETA,I_P
       
!C         DBRAT=0.001
!C         DBETA=DBRAT*BETA
          IF(.not.tFCIMC) WRITE(6,*) "DBETA=",DBETA

          IF(.NOT.TREAD) THEN
!             CALL WRITETMAT(NBASIS)
             IC=0
             HDiagTemp = get_helement(fDet, fDet, 0)
             WRITE(6,*) '<D0|H|D0>=',real(HDiagTemp,dp)
             WRITE(6,*) '<D0|T|D0>=',CALCT(FDET,NEL)
          
             IF(TUEG) THEN
!  The actual KE rather than the one-electron part of the Hamiltonian
                WRITE(6,*) 'Kinetic=',CALCT2(FDET,NEL,G1,ALAT,CST)
             ENDIF
          ENDIF

! Find out the number of alpha and beta electrons. For restricted calculations, these should be the same.
          if (tCSF) then
              nOccAlpha = (nel / 2) + LMS 
              nOccBeta =  (nel / 2) - LMS
          else
              nOccAlpha=0
              nOccBeta=0
              do i=1,NEl
                  CALL GETUNCSFELEC(FDET(I),J,IC)
                  IF(G1(J)%Ms.eq.1) THEN
                      ! Orbital is an alpha orbital
                     nOccAlpha=nOccAlpha+1
                  ELSE
                     nOccBeta=nOccBeta+1
                  ENDIF
              enddo
          end if
          WRITE(6,"(A,I5,A,I5,A)") " FDet has ",nOccAlpha," alpha electrons, and ",nOccBeta," beta electrons."
          ElecPairs=(NEl*(NEl-1))/2
          MaxABPairs=(nBasis*(nBasis-1)/2)

          ! And stats on the number of different types of electron pairs
          ! that can be found
          AA_elec_pairs = nOccAlpha * (nOccAlpha - 1) / 2
          BB_elec_pairs = nOccBeta * (nOccBeta - 1) / 2
          par_elec_pairs = AA_elec_pairs + BB_elec_pairs
          AB_elec_pairs = nOccAlpha * nOccBeta
          if (AA_elec_pairs + BB_elec_pairs + AB_elec_pairs /= ElecPairs) &
              call stop_all(this_routine, "Calculation of electron pairs failed")

          write(6,*) '    ', AA_elec_pairs, &
              ' alpha-alpha occupied electron pairs'
          write(6,*) '    ', BB_elec_pairs, &
              ' beta-beta occupied electron pairs'
          write(6,*) '    ', AB_elec_pairs, &
              ' alpha-beta occupied electron pairs'

          ! Get some stats about available numbers of holes, etc.
          ASSERT(.not. btest(nbasis, 0))
          norb = nbasis / 2
          nholes = nbasis - nel
          nholes_a = norb - nOccAlpha
          nholes_b = norb - nOccBeta

          ! And count the available hole pairs!
          hole_pairs = nholes * (nholes - 1) / 2
          AA_hole_pairs = nholes_a * (nholes_a - 1) / 2
          BB_hole_pairs = nholes_b * (nholes_b - 1) / 2
          AB_hole_pairs = nholes_a * nholes_b
          par_hole_pairs = AA_hole_pairs + BB_hole_pairs
          if (par_hole_pairs + AB_hole_pairs /= hole_pairs) &
              call stop_all(this_routine, "Calculation of hole pairs failed")

          write(6,*) '    ', AA_hole_pairs, 'alpha-alpha (vacant) hole pairs'
          write(6,*) '    ', BB_hole_pairs, 'beta-beta (vacant) hole pairs'
          write(6,*) '    ', AB_hole_pairs, 'alpha-beta (vacant) hole pairs'

          IF(tExactSizeSpace) THEN
              IF(ICILevel.eq.0) THEN
                  CALL FindSymSizeofSpace(6)
              ELSE
                  CALL FindSymSizeofTruncSpace(6)
              ENDIF
          endif
          IF(tMCSizeSpace) THEN
              CALL FindSymMCSizeofSpace(6) 
          ENDIF
          if(tMCSizeTruncSpace) then
              CALL FindSymMCSizeExcitLevel(6)
          endif

          IF(TMCDET) THEN
!C.. Generate the determinant from which we start the MC
             NLIST=1
             CALL GENSYMDETSS(MDK,NEL,G1,BRR,NBASIS,MCDET,NLIST,NBASISMAX)
             IF(NLIST.EQ.0) THEN
!C.. we couldn't find a det of that symmetry
                call stop_all(this_routine, 'Cannot find MC start determinant of correct symmetry')
             ENDIF
          ELSE
!C             CALL GENRANDOMDET(NEL,NBASIS,MCDET)
             DO I=1,NEL
                MCDET(I)=FDET(I)
             ENDDO
          ENDIF
          IF(TMONTE) THEN
             WRITE(6,"(A)",advance='no') 'MC Start Det: '
             call write_det (6, mcDet, .true.)
          ENDIF
!C.. we need to calculate a value for RHOEPS, so we approximate that
!C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this
!C.. we have put TMAT instead of ZIA
          IF(I_HMAX.NE.-20) THEN
!C.. If we're using rhos,
             RHOEPS=GETRHOEPS(RHOEPSILON,BETA,NEL,BRR,I_P)

!             WRITE(6,*) "RHOEPS:",RHOEPS
          ELSE
!C.. we're acutally diagonalizing H's, so we just leave RHOEPS as RHOEPSILON
             RHOEPS=RHOEPSILON
          ENDIF

        End Subroutine CalcInit
    
    
    
        subroutine CalcDoCalc(kp)
          use SystemData, only: Alat, Arr,Brr, Beta, ECore, G1, LMS, LMS2, nBasis,NMSH, nBasisMax
          use SystemData, only: SymRestrict, tCSFOLD, tParity, tSpn, ALat, Beta,tMolpro,tMolproMimic
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB,nEl
          Use DetCalcData, only : nDet, nEval, nmrks, w
          USE FciMCParMod , only : FciMCPar
          use RPA_Mod, only: RunRPA_QBA
          use DetCalc, only: CK, DetInv, tEnergy, tRead
          Use Determinants, only: FDet, nActiveBasis, SpecDet, tSpecDet
          use IntegralsData, only: FCK, NMAX, UMat, FCK
          use IntegralsData, only: HFEDelta, HFMix,nTay
          Use LoggingData, only: iLogging
          use Parallel_Calc
          use util_mod, only: get_free_unit, NECI_ICOPY
          use sym_mod
          use davidson_neci, only: DavidsonCalcType, DestroyDavidsonCalc
          use davidson_neci, only: davidson_direct_ci_init, davidson_direct_ci_end, perform_davidson
          use hamiltonian_linalg, only: direct_ci_type
          use kp_fciqmc, only: perform_kp_fciqmc, perform_subspace_fciqmc
          use kp_fciqmc_data_mod, only: tExcitedStateKP
          use kp_fciqmc_procs, only: kp_fciqmc_data
          use util_mod, only: int_fmt

          real(dp) :: EN,WeightDum, EnerDum
          real(dp), allocatable :: final_energy(:)
          integer :: iSeed, iunit, i
          type(kp_fciqmc_data), intent(inout) :: kp
          character(*), parameter :: this_routine = 'CalcDoCalc'
          type(DavidsonCalcType) :: davidsonCalc

          iSeed = 7

          IF (tMP2Standalone) then
              call ParMP2(FDet)
! Parallal 2v sum currently for testing only.
!          call Par2vSum(FDet)
          ELSE IF(tDavidson) then
              davidsonCalc = davidson_direct_ci_init(.true.)
              call perform_davidson(davidsonCalc, direct_ci_type, .true.)
              call davidson_direct_ci_end(davidsonCalc)
              call DestroyDavidsonCalc(davidsonCalc)

          ELSE IF(NPATHS.NE.0.OR.DETINV.GT.0) THEN
!Old and obsiolecte
!             IF(TRHOIJND) THEN
!C.. We're calculating the RHOs for interest's sake, and writing them,
!C.. but not keeping them in memory
!                  WRITE(6,*) "Calculating RHOS..."
!                  WRITE(6,*) "Using approx NTAY=",NTAY
!                  CALL CALCRHOSD(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,        &
!     &               NBASISMAX,G1,nBasis,BRR,NMSH,FCK,NMAX,ALAT,UMAT,             &
!     &               NTAY,RHOEPS,NWHTAY,ECORE)
!             ENDIF

             if(tFCIMC) then
                 call FciMCPar(final_energy)
                 if ((.not.tMolpro) .and. (.not.tMolproMimic)) then
                     if (allocated(final_energy)) then
                         do i = 1, size(final_energy)
                             write(6,'(1X,"Final energy estimate for state",1X,'//int_fmt(i)//',":",g25.14)') &
                                 i, final_energy(i)
                         end do
                     end if
                 end if
             elseif(tRPA_QBA) then
                call RunRPA_QBA(WeightDum,EnerDum)
                WRITE(6,*) "Summed approx E(Beta)=",EnerDum
             elseif(tKP_FCIQMC) then
                 if (tExcitedStateKP) then
                     call perform_subspace_fciqmc(kp)
                 else
                     call perform_kp_fciqmc(kp)
                 end if
              ENDIF
          endif
          IF(TMONTE.and..not.tMP2Standalone) THEN
!             DBRAT=0.01
!             DBETA=DBRAT*BETA
             WRITE(6,*) "I_HMAX:",I_HMAX
             WRITE(6,*) "Calculating MC Energy..."
             CALL neci_flush(6)
             IF(NTAY(1).GT.0) THEN
                WRITE(6,*) "Using approx RHOs generated on the fly, NTAY=",NTAY(1)
!C.. NMAX is now ARR
                call stop_all(this_routine, "DMONTECARLO2 is now non-functional.")
             ELSEIF(NTAY(1).EQ.0) THEN
                IF(TENERGY) THEN
                   WRITE(6,*) "Using exact RHOs generated on the fly"
!C.. NTAY=0 signifying we're going to calculate the RHO values when we
!C.. need them from the list of eigenvalues.   
!C.. Hide NMSH=NEVAL
!C..         FCK=W
!C..         ZIA=CK
!C..         UMAT=NDET
!C..         ALAT=NMRKS
!C..         NMAX=ARR
                call stop_all(this_routine, "DMONTECARLO2 is now non-functional.")
!                   EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,IMCSTEPS,             &
!     &                G1,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS,                                 &
!     &                NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,RHOEPS,NWHTAY,ILOGGING,ECORE,BETAEQ)
                ELSE
                   call stop_all(this_routine, "TENERGY not set, but NTAY=0" )
                ENDIF
             ENDIF
             WRITE(6,*) "MC Energy:",EN
!CC           WRITE(12,*) DBRAT,EN
          ENDIF
         
!C.. /AJWT
        End Subroutine CalcDoCalc

        Subroutine CalcCleanup()
           != Clean up (e.g. via deallocation) mess from Calc routines.
           use global_utilities
           character(*), parameter :: this_routine='CalcCleanup'

           deallocate(MCDet)
           call LogMemDealloc(this_routine,tagMCDet)

        End Subroutine CalcCleanup

      END MODULE Calc
      


      subroutine inpgetmethod(I_HMAX,NWHTAY,I_V)
         use constants
         use input_neci
         use CalcData, only: calcp_sub2vstar, calcp_logWeight, tMCDirectSum, &
                             g_multiweight, g_vmc_fac, tMPTheory, StarProd, &
                             tDiagNodes, tStarStars, tGraphMorph, tStarTrips, &
                             tHDiag, tMCStar, tFCIMC, tMCDets, tRhoElems, &
                             tReturnPathMC, tUseProcsAsNodes, tRPA_QBA, &
                             tDetermProj, tFTLM, TSpecLanc, tContTimeFCIMC, &
                             tExactSpec, tExactDiagAllSym
         use RPA_Mod, only: tDirectRPA
         use LoggingData, only: tCalcFCIMCPsi
         implicit none
         integer I_HMAX,NWHTAY,I_V
         CHARACTER(LEN=16) w
         do while ( item .lt. nitems )
           call readu(w)
           select case(w)
           case("VERTEX")
               call readu(w)
               select case(w)
               case("FCIMC")
                   I_HMAX=-21
                   TFCIMC=.true.
                   tUseProcsAsNodes=.true.
                   do while(item.lt.nitems)
                      call readu(w)
                      select case(w)
                      case("CONT-TIME")
                          tContTimeFCIMC = .true.
                      case("MCDIFFUSION")
!                          TMCDiffusion=.true.
                          CALL Stop_All("inpgetmethod","MCDIFFUSION option depreciated")
                      case("RESUMFCIMC")
!                          TResumFCIMC=.true.
                          CALL Stop_All("inpgetmethod","MCDIFFUSION option depreciated")
                      case default
                          call report("Keyword error with "//trim(w),.true.)
                      endselect
                   enddo
               case("RPA")
                  tRPA_QBA=.true.
                  tDirectRPA=.false.
                  do while(item.lt.nitems)
                      call readu(w)
                      select case(w)
                      case("DIRECT")
                          tDirectRPA=.true.
                      endselect
                  enddo
               case("RETURNPATHMC")
                   I_HMAX=-21
                   TReturnPathMC=.true.
                   call readu(w)
                   select case(w)
                   case("RHOELEMS")
                       TRhoElems=.true.
                   endselect
               case("MCDets")
                   I_HMAX=-21
                   TMCDets=.true.
               case("SUM")
                  do while(item.lt.nitems)
                   call readu(w)
                   select case(w)
                   case("OLD")
                       I_HMAX = -1
                   case("NEW")
                       I_HMAX = -8
                   case("HDIAG")
                       I_HMAX = -20
                   case("READ")
                       I_HMAX=-14
                   case("SUB2VSTAR")
                       CALCP_SUB2VSTAR=.TRUE.
                   case("LOGWEIGHT")
                       CALCP_LOGWEIGHT=.TRUE.
                   case default
                       call report("Error - must specify OLD or NEW vertex sum method",.true.)
                   end select
                  enddo
               case("MC","MCMETROPOLIS")
                  I_HMAX = -7
                   call readu(w)
                   select case(w)
                   case("HDIAG")
                       I_HMAX = -19
                   end select
                  tMCDirectSum=.FALSE.
                  IF(I_V.GT.0) g_MultiWeight(I_V)=1.0_dp
               case("MCDIRECT")
                  I_HMAX = -7
                  tMCDirectSum=.TRUE.
                   call readu(w)
                   select case(w)
                   case("HDIAG")
                       I_HMAX = -19
                   end select
                  G_VMC_FAC=0.0_dp
               case("MCMP")
                  tMCDirectSum=.TRUE.
                  I_HMAX = -19
                  G_VMC_FAC=0.0_dp
                  TMPTHEORY=.TRUE.
               case("GRAPHMORPH")
                   TGraphMorph=.true.
                   I_HMAX=-21
                   call readu(w)
                   select case(w)
                   case("HDIAG")
                       !If this is true, then it uses the hamiltonian matrix to determinant coupling to excitations, 
                       !and to diagonalise to calculate the energy
                       THDiag=.true.
                   endselect
               case("STAR")
                  I_HMAX=0
                  do while(item.lt.nitems)
                     call readu(w)
                     select case(w)
                     case("NEW")
                        I_HMAX=-21
                     case("OLD")
                        I_HMAX=-9
                     case("NODAL")
                        TDIAGNODES=.TRUE.
                     case("STARSTARS")
                         TSTARSTARS=.true.
                     case("MCSTAR")
                         NWHTAY=IBSET(NWHTAY,0)
                         TMCSTAR=.true.
                     case("STARPROD")
                        STARPROD=.TRUE.
                     case("TRIPLES")
                         TStarTrips=.TRUE.
                     case("COUNTEXCITS")
                        NWHTAY=IBSET(NWHTAY,8)
                     case("ADDSINGLES")
                        NWHTAY=IBSET(NWHTAY,7)
                        IF(I_HMAX.NE.-21)  call report(        &
     &                     "Error - cannot use ADDSINGLES"     &
     &                     //" without STAR NEW",.true.)
                     case("DIAG")
                         NWHTAY=IBCLR(NWHTAY,0)
                     case("POLY")
                         NWHTAY=IBSET(NWHTAY,0)
                     case("POLYMAX")
                         NWHTAY=IBSET(NWHTAY,0)
                         NWHTAY=IBSET(NWHTAY,1)
                     case("POLYCONVERGE")
                         NWHTAY=IBSET(NWHTAY,0)
                         NWHTAY=IBSET(NWHTAY,2)
                     case("POLYCONVERGE2")
                         NWHTAY=IBSET(NWHTAY,0)
                         NWHTAY=IBSET(NWHTAY,6)
                     case("H0")
                         NWHTAY=IBSET(NWHTAY,5)
                         if(I_HMAX.ne.-21) call report ("H0 "  &
     &              //"can only be specified with POLY... NEW")
                     case default
                       call report("Error - must specify DIAG" &
     &               //" or POLY vertex star method",.true.)
                      end select
                  enddo
!                  IF(TSTARSTARS.and..not.BTEST(NWHTAY,0)) THEN 
!                      call report("STARSTARS must be used with " &
!     &                 //"a poly option",.true.)
!                  ENDIF
                  IF(STARPROD.and.BTEST(NWHTAY,0)) THEN
                      call report("STARPROD can only be "      &
     &               //"specified with DIAG option",.true.)
                   ENDIF
                  if(i_hmax.eq.0)                              &
     &          call report("OLD/NEW not specified for STAR",  &
     &                 .true.)
               case("DETERM-PROJ")
                   tDetermProj = .true.
                   I_HMAX=-21
                   TFCIMC=.true.
                   tUseProcsAsNodes=.true.
               case("FTLM")
                   tFTLM = .true.
                   I_HMAX=-21
                   TFCIMC=.true.
                   tUseProcsAsNodes=.true.
               case("EXACT-SPECTRUM")
                   tExactSpec = .true.
                   I_HMAX=-21
                   TFCIMC=.true.
                   tUseProcsAsNodes=.true.
               case("EXACT-DIAG")
                   tExactDiagAllSym = .true.
                   I_HMAX=-21
                   TFCIMC=.true.
                   tUseProcsAsNodes=.true.
               case("SPECTRAL-LANCZOS")
                   tSpecLanc = .true.
                   I_HMAX=-21
                   TFCIMC=.true.
                   tUseProcsAsNodes=.true.
               case default
                   call report("Keyword error with "//trim(w),     &
         &                 .true.)
               end select
           case default
               call report("Error.  Method not specified."     &
     &           //" Stopping.",.true.)
           end select
         end do

      end subroutine inpgetmethod



      subroutine inpgetexcitations(NWHTAY,w)
         use input_neci
         IMPLICIT NONE
         INTEGER NWHTAY
         CHARACTER(LEN=16) w
!         call readu(w)
         select case(w)
         case("FORCEROOT")
            NWHTAY=IOR(NWHTAY,1)
         case("FORCETREE")
            NWHTAY=IOR(NWHTAY,2)
         case("SINGLES")
            NWHTAY=IOR(NWHTAY,8)
         case("DOUBLES")
            NWHTAY=IOR(NWHTAY,16)
         case("ALL")
            NWHTAY=0
         case default
            call report("Keyword error with EXCITATIONS "//trim(w),.true.)
         end select
      end subroutine inpgetexcitations




! Given an input RHOEPSILON, create Fermi det D out of lowest orbitals and get RHOEPS (which is rhoepsilon * exp(-(beta/P)<D|H|D>
      FUNCTION GETRHOEPS(RHOEPSILON,BETA,NEL,BRR,I_P)
         Use Determinants, only: get_helement, write_det
         use constants, only: dp
         use SystemData, only: BasisFN
         use sort_mod
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I,I_P
         INTEGER BRR(*)
         real(dp) RHOEPSILON,BETA,GETRHOEPS
         HElement_t(dp) BP, tmp
         DO I=1,NEL
            NI(I)=BRR(I)
         ENDDO
         call sort (nI)
         BP=-BETA/I_P
         tmp = RHOEPSILON * exp(BP*get_helement(nI, nI, 0))
         GETRHOEPS = sqrt(tmp * tmp)
         RETURN
      END FUNCTION GetRhoEps



! Calculate the kinetic energy of the UEG (this differs from CALCT by including the constant CST
      FUNCTION CALCT2(NI,NEL,G1,ALAT,CST)
         use constants, only: dp
         use SystemData, only: BasisFN, kvec, k_lattice_constant, TUEG2
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I,J
         TYPE(BasisFN) G1(*)
         real(dp) ALAT(4),CST,TMAT,CALCT2
         LOGICAL ISCSF_old

         CALCT2=0.0_dp
         IF(iscsf_old(NI,NEL)) RETURN

         !===============================
         if (TUEG2) then
            DO J=1,NEL
                 I=NI(J)
                 TMAT=real(kvec(I, 1)**2+kvec(I, 2)**2+kvec(I, 3)**2, dp)
                 TMAT=0.5_dp*TMAT*k_lattice_constant**2
                 CALCT2=CALCT2+TMAT
            ENDDO
            return
        end if ! TUEG2
         !===============================
         DO J=1,NEL
            I=NI(J)
           TMAT=((ALAT(1)**2)*((G1(I)%K(1)**2)/(ALAT(1)**2)+   &
               (G1(I)%K(2)**2)/(ALAT(2)**2)+                   &
               (G1(I)%K(3)**2)/(ALAT(3)**2)))
           TMAT=TMAT*CST
           CALCT2=CALCT2+TMAT
         ENDDO
         RETURN
      END FUNCTION CALCT2
