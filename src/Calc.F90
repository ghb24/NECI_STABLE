#include "macros.h"

MODULE Calc
        
    use CalcData
    use SystemData, only: beta, nel, STOT, tCSF, LMS, tSpn, AA_elec_pairs, &
                          BB_elec_pairs, par_elec_pairs, AB_elec_pairs, &
                          AA_hole_pairs, BB_hole_pairs, AB_hole_pairs, &
                          par_hole_pairs, hole_pairs, nholes_a, nholes_b, &
                          nholes
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
    use CCMCData, only: dInitAmplitude, dProbSelNewExcitor, nSpawnings, &
                        tSpawnProp, nClustSelections, tExactEnergy,     &
                        dClustSelectionRatio,tSharedExcitors
    use FciMCData, only: tTimeExit,MaxTimeExit, InputDiagSft, tSearchTau, &
                         nWalkerHashes, tHashWalkerList, HashLengthFrac, &
                         tTrialHash, tIncCancelledInitEnergy, MaxTau, &
                         tStartCoreGroundState, pParallel, pops_pert, &
                         alloc_popsfile_dets
    use semi_stoch_gen, only: core_ras
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
          InitialPart=1
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
          tStartCoreGroundState = .true.
          tHashWalkerList=.false.
          HashLengthFrac=0.0_dp
          nWalkerHashes=0
          tTrialHash=.true.
          tIncCancelledInitEnergy = .false.
          iExitWalkers=-1
          FracLargerDet=1.2_dp
          iReadWalkersRoot=0 
          tShiftonHFPop=.false.
          MaxWalkerBloom=-1
          tSearchTau=.true.
          InputDiagSft=0.0_dp
          tTimeExit=.false.
          MaxTimeExit=0.0_dp
          tMaxBloom=.false.
          iRestartWalkNum=0
          iWeightPopRead=1.0e-12
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
          MemoryFacSpawn=0.5_dp
          MemoryFacInit = 0.3_dp
          TStartSinglePart=.true.
          TFixParticleSign=.false.
          TProjEMP2=.false.
          THFRetBias=.false.
          TSignShift=.false.
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
          TCCMC=.false.
          TMCDets=.false.
          TBinCancel=.false.  
          ScaleWalkers=1.0_dp
          TReadPops=.false.
          iPopsFileNoRead = 0
          iPopsFileNoWrite = 0
          tWalkContGrow=.false.
          StepsSft=100
          SftDamp=10.0
          Tau=0.0_dp
          InitWalkers=3000
          dInitAmplitude=1.0_dp
          dProbSelNewExcitor=0.7_dp
          nSpawnings=1
          nClustSelections=1
          dClustSelectionRatio=1
          tExactEnergy=.false.
          tSharedExcitors=.false.
          tSpawnProp=.false.
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
          g_MultiWeight(:)=0
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
          G_VMC_FAC = 16
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
          BETAEQ = 0
          TMCDET = .false.
          MDK(:) = 0
          DETINV = 0
          TSPECDET = .false.
          TTROT=.true.
          BETA = 1000
          BETAP=1.0e-4_dp
          TBETAP=.false.
          RHOEPSILON=1.0e-6_dp
          DBETA=-1.0_dp
          GraphEpsilon=0
          PGenEpsilon=0
          StarConv=1.0e-3_dp
          calcp_sub2vstar=.false.
          calcp_logweight=.false.
          TENPT=.false.
          TLADDER=.false. 
          tDefineDet=.false.
          tTruncInitiator=.false.
          tAddtoInitiator=.false.
          InitiatorWalkNo=10.0_dp
          tInitIncDoubs=.false.
          MaxNoatHF=0
          HFPopThresh=0
          tSpatialOnlyHash = .false.
          tNeedsVirts=.true.! Set if we need virtual orbitals  (usually set).  Will be unset 
          !(by Calc readinput) if I_VMAX=1 and TENERGY is false

          lNoTriples=.false.
          tReadPopsChangeRef = .false.
          tReadPopsRestart = .false.
          iLogicalNodeSize = 0 !Meaning use the physical node size
          tAllRealCoeff=.false.
          tEnhanceRemainder=.true.
          tUseRealCoeffs = .false.
          tRealCoeffByExcitLevel=.false.
          RealCoeffExcitThresh=2
          tRealSpawnCutoff=.false.
          RealSpawnCutoff=1.0e-5_dp
          OccupiedThresh=1.0_dp
          tInitOccThresh=.false.
          InitiatorOccupiedThresh=0.1_dp
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
          spin_proj_shift = 0
          spin_proj_cutoff = 0
          spin_proj_iter_count = 1
          spin_proj_nopen_max = -1
          disable_spin_proj_varyshift = .false.
          tUseProcsAsNodes=.false.

          ! Truncation based on number of unpaired electrons
          tTruncNOpen = .false.

          hash_shift=0
          tUniqueHFNode = .false.

          ! Semi-stochastic and trial wavefunction options.
          tSemiStochastic = .false.
          tCSFCore = .false.
          tDoublesCore = .false.
          tCASCore = .false.
          tRASCore = .false.
          tOptimisedCore = .false.
          tFCICore = .false.
          tHeisenbergFCICore = .false.
          tPopsCore = .false.
          tReadCore = .false.
          tLowECore = .false.
          tMP1Core = .false.
          num_det_generation_loops = 1
          n_core_pops = 0
          low_e_core_excit = 0
          low_e_core_num_keep = 0
          semistoch_mp1_ndets = 0
          tLowECoreAllDoubles = .false.
          tLimitDetermSpace = .false.
          tLimitTrialSpace = .false.
          max_determ_size = 0
          max_trial_size = 0
          tDetermAmplitudeCutoff = .false.
          tTrialWavefunction = .false.
          tDoublesTrial = .false.
          tCASTrial = .false.
          tOptimisedTrial =.false.
          tPopsTrial = .false.
          tReadTrial = .false.
          tLowETrial = .false.
          tMP1Trial = .false.
          tFCITrial = .false.
          tHeisenbergFCITrial = .false.
          num_trial_generation_loops = 1
          n_trial_pops = 0
          low_e_trial_excit = 0
          low_e_trial_num_keep = 0
          trial_mp1_ndets = 0
          tLowETrialAllDoubles = .false.
          tTrialAmplitudeCutoff = .false.
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

          pParallel = 0.5

          InitiatorCutoffEnergy = 99.99e99_dp
          InitiatorCutoffWalkNo = 99.0_dp

          tSurvivalInitiatorThreshold = .false.
          tSurvivalInitMultThresh = .false.
          im_time_init_thresh = 0.1_dp
          init_survival_mult = 3.0_dp
          MaxTau = 1.0_dp

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
          use CCMCData, only: dInitAmplitude,dProbSelNewExcitor,nSpawnings,tSpawnProp,nClustSelections
          use CCMCData, only: tExactEnergy,tSharedExcitors
          use FciMCData, only: hash_shift, davidson_ras
          use ras_data
          use global_utilities
          use Parallel_neci, only : nProcessors
          use LoggingData, only: tLogDets
          IMPLICIT NONE
          LOGICAL eof
          CHARACTER (LEN=100) w
          CHARACTER (LEN=100) input_string
          CHARACTER(*),PARAMETER :: t_r='CalcReadInput'
          integer :: l, i, j, ierr
          integer :: tempMaxNoatHF,tempHFPopThresh
          logical :: tExitNow
          integer :: ras_size_1, ras_size_2, ras_size_3, ras_min_1, ras_max_3
          integer :: npops_pert, npert_spectral_left, npert_spectral_right

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

            case("CYCLES")
                call readi(NWHTAY(1,1))
                if ( I_HMAX .ne. -7.and.                              &
     &               I_HMAX .ne. -19) then
                    call report(trim(w)//" only valid for MC "        &
     &                 //"method",.true.)
                end if
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
                ALLOCATE(DefDet(NEl),stat=ierr)
                CALL LogMemAlloc('DefDet',NEl,4,t_r,tagDefDet,ierr)
                DefDet(:)=0
                do i=1,NEl
                    call geti(DefDet(i))
                enddo

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
            case("INITAMPLITUDE")
!For Amplitude CCMC the initial amplitude.
                call getf(dInitAmplitude)
            case("CLUSTERSIZEBIAS")
                call getf(dProbSelNewExcitor)
            case("NSPAWNINGS")
!For Amplitude CCMC the number of spawnings for each cluster.
                call geti(nSpawnings)
            case("HASH_SHIFT")
                call geti(hash_shift)
            case("NCLUSTSELECTIONS")
!For Particle CCMC the number of  cluster.
                call geti(nClustSelections)
            case("CLUSTSELECTIONRATIO")
!For Particle CCMC the number of  cluster.
                call getf(dClustSelectionRatio)
            case("CCMCEXACTENERGY")
               tExactEnergy=.true.
            case("CCMCSHAREDEXCITORS")
               tSharedExcitors=.true.
            case("SPAWNPROP")
!For Amplitude CCMC use NSPAWNINGS as a total number of spawnings, and distribute them according to the Amplitudes of clusters.
               tSpawnProp=.true.
            case("NMCYC")
!For FCIMC, this is the number of MC cycles to perform
                call geti(NMCyc)
            case("DIAGSHIFT")
!For FCIMC, this is the amount extra the diagonal elements will be shifted. This is proportional to the deathrate of 
!walkers on the determinant
                call getf(InputDiagSft)

            case("TAUFACTOR")
!For FCIMC, this is the factor by which 1/(HF connectivity) will be multiplied by to give the timestep for the calculation.
                tSearchTau=.false.  !Tau is set, so don't search for it.
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
                    case default
                        tSearchTau = .false.
                    end select
                else
                    tSearchTau = .false.
                end if

            case("MAX-TAU")
                ! For tau searching, set a maximum value of tau. This places
                ! a limit to prevent craziness at the start of a calculation
                call getf(MaxTau)

            case("MAXWALKERBLOOM")
                !Set the maximum allowed walkers to create in one go, before reducing tau to compensate.
                call geti(MaxWalkerBloom)
            case("SHIFTDAMP")
!For FCIMC, this is the damping parameter with respect to the update in the DiagSft value for a given number of MC cycles.
                call getf(SftDamp)
            case("LINSCALEFCIMCALGO")
                !Use the linear scaling FCIMC algorithm
                !Instead of the absolute length of the hash table, read in the fraction of initwalkers that it wants to be.
!                call geti(nWalkerHashes)
                tHashWalkerList=.true.
                if(item.lt.nitems) then
                    call getf(HashLengthFrac)
                else
                    HashLengthFrac=0.7_dp
                endif
            case("SEMI-STOCHASTIC")
                tSemiStochastic = .true.
                if (item < nitems) then
                    call geti(semistoch_mp1_ndets)
                    tMP1Core = .true.
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
                tDoublesCore = .true.
            case("CAS-CORE")
                tCASCore = .true.
                tSpn = .true.
                call geti(OccDetermCASOrbs)  !Number of electrons in CAS 
                call geti(VirtDetermCASOrbs)  !Number of virtual spin-orbitals in CAS
            case("RAS-CORE")
                tRASCore = .true.
                call geti(ras_size_1)  ! Number of spatial orbitals in RAS1.
                call geti(ras_size_2)  ! Number of spatial orbitals in RAS2.
                call geti(ras_size_3)  ! Number of spatial orbitals in RAS3.
                call geti(ras_min_1)  ! Min number of electrons (alpha and beta) in RAS1 orbs. 
                call geti(ras_max_3)  ! Max number of electrons (alpha and beta) in RAS3 orbs.
                core_ras%size_1 = int(ras_size_1,sp)
                core_ras%size_2 = int(ras_size_2,sp)
                core_ras%size_3 = int(ras_size_3,sp)
                core_ras%min_1 = int(ras_min_1,sp)
                core_ras%max_3 = int(ras_max_3,sp)
            case("OPTIMISED-CORE")
                tOptimisedCore = .true.
            case("OPTIMISED-CORE-CUTOFF-AMP")
                tDetermAmplitudeCutoff = .true.
                num_det_generation_loops = nitems - 1
                allocate(determ_space_cutoff_amp(num_det_generation_loops))
                do I = 1, num_det_generation_loops
                    call getf(determ_space_cutoff_amp(I))
                end do
            case("OPTIMISED-CORE-CUTOFF-NUM")
                tDetermAmplitudeCutoff = .false.
                num_det_generation_loops = nitems - 1
                allocate(determ_space_cutoff_num(num_det_generation_loops))
                do I = 1, num_det_generation_loops
                    call geti(determ_space_cutoff_num(I))
                end do
            case("FCI-CORE")
                tFCICore = .true.
            case("HEISENBERG-FCI-CORE")
                tHeisenbergFCICore = .true.
            case("POPS-CORE")
                tPopsCore = .true.
                call geti(n_core_pops)
            case("READ-CORE")
                tReadCore = .true.
            case("LOW-ENERGY-CORE")
    ! Input values: The first integer is the maximum excitation level to go up to.
    !               The second integer is the maximum number of states to keep for a subsequent iteration.
    !               If desired, you can put "All-Doubles" after these two integers to keep all singles and doubles.
    !               If max-core-size is specified then this value will be used to select the number of states kept
    !               after the *final* iteration.
                tLowECore = .true.
                call geti(low_e_core_excit)
                call geti(low_e_core_num_keep)
                if (nitems > 3) then
                    call geta(input_string)
                    if (trim(input_string) == "All-Doubles") then
                        tLowECoreAllDoubles = .true.
                    else
                        call stop_all("SysReadInput","Input string is not recognised.")
                    end if
                end if
            case("MAX-CORE-SIZE")
                tLimitDetermSpace = .true.
                call geti(max_determ_size)
            case("MAX-TRIAL-SIZE")
                tLimitTrialSpace = .true.
                call geti(max_trial_size)
            case("TRIAL-WAVEFUNCTION")
                tTrialWavefunction = .true.
                if (item < nitems) then
                    call geti(trial_mp1_ndets)
                    tMP1Trial = .true.
                end if
            case("DOUBLES-TRIAL")
                tDoublesTrial = .true.
            case("CAS-TRIAL")
                tCASTrial = .true.
                tSpn = .true.
                call geti(OccTrialCASOrbs)  !Number of electrons in CAS 
                call geti(VirtTrialCASOrbs)  !Number of virtual spin-orbitals in CAS
            case("OPTIMISED-TRIAL")
                tOptimisedTrial = .true.
            case("OPTIMISED-TRIAL-CUTOFF-AMP")
                tTrialAmplitudeCutoff = .true.
                num_trial_generation_loops = nitems - 1
                allocate(trial_space_cutoff_amp(num_trial_generation_loops))
                do I = 1, num_trial_generation_loops
                    call getf(trial_space_cutoff_amp(I))
                end do
            case("OPTIMISED-TRIAL-CUTOFF-NUM")
                tTrialAmplitudeCutoff = .false.
                num_trial_generation_loops = nitems - 1
                allocate(trial_space_cutoff_num(num_trial_generation_loops))
                do I = 1, num_trial_generation_loops
                    call geti(trial_space_cutoff_num(I))
                end do
            case("POPS-TRIAL")
                tPopsTrial = .true.
                call geti(n_trial_pops)
            case("READ-TRIAL")
                tReadTrial = .true.
            case("LOW-ENERGY-TRIAL")
    ! Input values: The first integer is the maximum excitation level to go up to.
    !               The second integer is the maximum number of states to keep for a subsequent iteration.
    !               If desired, you can put "All-Doubles" after these two integers to keep all singles and doubles.
    !               If max-trial-size is specified then this value will be used to select the number of states kept
    !               after the *final* iteration.
                tLowETrial = .true.
                call geti(low_e_trial_excit)
                call geti(low_e_trial_num_keep)
                if (nitems > 3) then
                    call geta(input_string)
                    if (trim(input_string) == "All-Doubles") then
                        tLowETrialAllDoubles = .true.
                    else
                        call stop_all("SysReadInput","Input string is not recognised.")
                    end if
                end if
            case("FCI-TRIAL")
                tFCITrial = .false.
            case("HEISENBERG-FCI-TRIAL")
                tHeisenbergFCITrial = .false.
            case("TRIAL-BIN-SEARCH")
                tTrialHash = .false.
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
                call geti(StepsSft)
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
                tCheckHighestPop=.true.
                tChangeProjEDet=.true.
                IF(item.lt.nitems) then
                    call Getf(FracLargerDet)
                ENDIF
                
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

            case("INITIATOR-ENERGY-CUTOFF")
                !
                ! Specify both a threshold an an addtoinitiator value for
                ! varying the thresholds
                call getf(InitiatorCutoffEnergy)
                call getf(InitiatorCutoffWalkNo)

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
                tEnhanceRemainder=.false.
                !When we do the removal step with AllRealCoeff, on the occasions where these pops are *not* removed,
                !Keep their population the same, rather than resetting as a value of 1 (which is technically correct)
                !This "bug" produced initiator-like (no plateau) behaviour, so may be of interest
            case("REALSPAWNCUTOFF")
                tRealSpawnCutoff=.true.
                call Getf(RealSpawnCutoff)
            case("SETOCCUPIEDTHRESH")
                call Getf(OccupiedThresh)
            case("SETINITOCCUPIEDTHRESH")
                tInitOccThresh=.true.
                tAllRealCoeff=.true.
                call Getf(InitiatorOccupiedThresh)

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
                tSurvivalInitiatorThreshold = .true.
                if (item < nitems) then
                    call readf(im_time_init_thresh)
                end if

            case("INITIATOR-SURVIVAL-MULTIPLIER")
                ! If a site survives for a certain multiple of how long it
                ! would _expect_ to have survived, then it should be treated
                ! as an initiator
                ! --> A more flexible version of INITIATOR-SURVIVAL-CRITERION
                tSurvivalInitMultThresh = .true.
                if (item < nitems) then
                    call readf(init_survival_mult)
                end if

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
          HElement_t HDiagTemp
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
!             STOP "G_VNC_FAC LE 0"
!          ENDIF

          IF(BETAP.NE.0) THEN 
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
             WRITE(6,*) '<D0|H|D0>=',HDiagTemp
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
                STOP 'Cannot find MC start determinant of correct symmetry'
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
    
    
    
        Subroutine CalcDoCalc()
          use SystemData, only: Alat, Arr,Brr, Beta, ECore, G1, LMS, LMS2, nBasis,NMSH, nBasisMax
          use SystemData, only: SymRestrict, tCSFOLD, tParity, tSpn, ALat, Beta,tMolpro,tMolproMimic
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB,nEl
          Use DetCalcData, only : nDet, nEval, nmrks, w
          USE FciMCParMod , only : FciMCPar
          use RPA_Mod, only: RunRPA_QBA
          USE CCMC, only: CCMCStandalone,CCMCStandaloneParticle
          use CCMCData, only: tAmplitudes
          use DetCalc, only: CK, DetInv, tEnergy, tRead
          Use Determinants, only: FDet, nActiveBasis, SpecDet, tSpecDet
          use IntegralsData, only: FCK, NMAX, UMat, FCK
          use IntegralsData, only: HFEDelta, HFMix,nTay
          Use LoggingData, only: iLogging
          use Parallel_Calc
          use util_mod, only: get_free_unit, NECI_ICOPY
          use sym_mod
          use davidson_neci, only: davidson_direct_ci_init, davidson_direct_ci_end, perform_davidson
          use davidson_neci, only: direct_ci_type
          use kp_fciqmc, only: perform_kp_fciqmc
          use kp_fciqmc_procs, only: kp

!Calls
!          real(dp) DMonteCarlo2
!Local Vars
          real(dp) EN,WeightDum,EnerDum
          integer iSeed,iunit
          iSeed=7 

!C.. we need to calculate a value for RHOEPS, so we approximate that
!C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this 
!C.. If we haven't already calced RHOEPS, do it now
          Call DoExactVertexCalc()

          IF (tMP2Standalone) then
              call ParMP2(FDet)
! Parallal 2v sum currently for testing only.
!          call Par2vSum(FDet)
          ELSE IF(tDavidson) then
              call davidson_direct_ci_init()
              call perform_davidson(direct_ci_type, .true.)
              call davidson_direct_ci_end()
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
                 call FciMCPar(WeightDum,EnerDum)

                 if((.not.tMolpro).and.(.not.tMolproMimic)) WRITE(6,*) "Summed approx E(Beta)=",EnerDum
             elseif(tCCMC) then
                  if(tAmplitudes) THEN
                     CALL CCMCStandAlone(WeightDum,EnerDum)
                  else
                     CALL CCMCStandaloneParticle(WeightDum,EnerDum)
                  endif
                  WRITE(6,*) "Summed approx E(Beta)=",EnerDum
             elseif(tRPA_QBA) then
                call RunRPA_QBA(WeightDum,EnerDum)
                WRITE(6,*) "Summed approx E(Beta)=",EnerDum
             elseif(tKP_FCIQMC) then
                 call perform_kp_fciqmc(kp)
             else


                 WRITE(6,*) "Calculating ",NPATHS," W_Is..."
                 iunit =get_free_unit()
                 IF(BTEST(ILOGGING,1)) THEN
                    IF(I_HMAX.EQ.-10) THEN
                       OPEN(iunit,FILE="MCSUMMARY",STATUS="UNKNOWN")
                       WRITE(iunit,*) "Calculating ",NPATHS," W_Is..."
                       CLOSE(iunit)
                    ELSE
                       OPEN(iunit,FILE="MCPATHS",STATUS="UNKNOWN")
                       WRITE(iunit,*) "Calculating ",NPATHS," W_Is..."
                       CLOSE(iunit)
                    ENDIF
                 ENDIF
                 IF(NTAY(1).GT.0) THEN
                    WRITE(6,*) "Using list of determinants."
                    WRITE(6,*) "Using approx RHOs generated on the fly,NTAY=",NTAY(1)
    !C.. we haven't calculated the energy, so we're calculating the weights
    !C.. with approx RHOs
                    IF(TENERGY) THEN
    !C.. If we've generated a list of dets
    !C.. Instead of NMAX, we put ARR
                         CALL CALCRHOPII2(BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,               &
         &                 NBASISMAX,G1,nBasis,BRR,NMSH,FCK,ARR,ALAT,UMAT,NTAY,          &
         &                 RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,           &
         &                 DETINV,TSPECDET,SPECDET)
    !                      WRITE(6,*) "Out Here 2"
    !                      CALL neci_flush(6)
                    ELSE
                       IF(TCSFOLD) THEN
                          IF(.NOT.TSPECDET) THEN
                             WRITE(6,*) "SPECDET not specified. Using Fermi determinant ONLY"
                             TSPECDET=.TRUE.
                             CALL NECI_ICOPY(NEL,FDET,1,SPECDET,1)
                          ENDIF
                       ENDIF
                      IF(TPARITY) THEN
                          WRITE(6,*) "Using symmetry restriction:"
                          CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
                      ENDIF
                      IF(TSPN) THEN
                          WRITE(6,*) "Using spin restriction:",LMS
                      ENDIF
    !C.. Instead of NMAX we have ARR 
                      CALL CALCRHOPII3(BETA,I_P,I_HMAX,I_VMAX,NEL,                         &
         &               NBASISMAX,G1,nBasis,BRR,NMSH,FCK,ARR,ALAT,UMAT,NTAY,              &
         &               RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,               &
         &               DETINV,TSPN,LMS2,TPARITY,SymRestrict,TSPECDET,SPECDET,            &
         &               nActiveBasis)
                    ENDIF
                 ELSE
                     WRITE(6,*) "Invalid combination of NTAY and TENERGY.  No NPATHS calculated"
                     WRITE(6,*) "NTAY: ",NTAY(1)," TENERGY: ",TENERGY
                 ENDIF
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
                STOP "DMONTECARLO2 is now non-functional."
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
                STOP "DMONTECARLO2 is now non-functional."
!                   EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,IMCSTEPS,             &
!     &                G1,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS,                                 &
!     &                NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,RHOEPS,NWHTAY,ILOGGING,ECORE,BETAEQ)
                ELSE
                   STOP "TENERGY not set, but NTAY=0" 
                ENDIF
             ENDIF
             WRITE(6,*) "MC Energy:",EN
!CC           WRITE(12,*) DBRAT,EN
          ENDIF
         
!C.. /AJWT
        End Subroutine CalcDoCalc

        Subroutine DoExactVertexCalc()
          use SystemData, only: Alat, Beta, Brr, ECORE, G1, nBasis, nBasisMax,nMsh, Arr,nEl
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB
          use IntegralsData, only: fck, nMax, UMat,nTay
          Use DetCalc, only: cK, nDet, nEval, tEnergy, tRead, W, NMRKS, DetInv
          Use Determinants, only: specdet, tSpecDet
          Use LoggingData, only: iLogging
          Use util_mod, only: get_free_unit
          Use DetCalc, only: tFindDets
          use sym_mod
          real(dp) flri, flsi
          real(dp) En, ExEn, GSEn
          real(dp) RH
          INTEGER iDeg, III, iunit
          Type(BasisFN) iSym
          LOGICAL tWarn
          
          real(dp) CalcMCEn, CalcDLWDB, DoExMC
            
          IF(TENERGY.and.(.not.tFindDets)) THEN
             RHOEPS=RHOEPSILON*EXP(-BETA*(W(1))/I_P)
!            WRITE(6,*) "RHOEPS:",RHOEPS
             IF(TREAD) THEN
                EXEN=CALCMCEN(NEVAL,W,BETA)
                WRITE(6,"(A,F19.5)") "EXACT E(BETA)=",EXEN
                GSEN=CALCDLWDB(1,NDET,NEVAL,CK,W,BETA)
                WRITE(6,"(A,F19.5)") "EXACT DLWDB(D0)=",GSEN
             ENDIF
             iunit = get_free_unit()
             OPEN(iunit,FILE='RHOPIIex',STATUS='UNKNOWN')
             IF(NDETWORK.EQ.0.OR.NDETWORK.GT.NDET) NDETWORK=NDET
             DO III=1,NDETWORK
             
                CALL CALCRHOPII(III,NDET,NEVAL,CK,W,BETA,I_P,FLRI,FLSI,TWARN)
                IF(TWARN) THEN
                   IF(III.EQ.1) THEN
                      WRITE(6,*) "Warning received from CALCRHOPII."
                      IF(TREAD) THEN
                         WRITE(6,*) "TREAD set. Cannot calculate RHOII."
                      ELSE
                      WRITE(6,*) "Calculating RHOII using 1st order Taylor."
                      ENDIF
                   ENDIF
                   IF(.NOT.TREAD) THEN
                   CALL CALCRHO2(NMRKS(1:NEl,III),NMRKS(1:NEl,III),BETA,I_P,NEL, G1,NBASIS,NMSH, &
                    FCK,NMAX,ALAT,UMAT,RH,1,0,ECORE)
!C                   WRITE(6,*) RH
                   FLRI=LOG(RH)
                   FLSI=FLSI-I_P*FLRI
                   ENDIF
                ENDIF
                call write_det (iunit, NMRKS(:,III), .false.)
                GSEN=CALCDLWDB(III,NDET,NEVAL,CK,W,BETA)
                CALL GETSYM(NMRKS(:,III),NEL,G1,NBASISMAX,ISYM)
                CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
                WRITE(iunit,"(4G25.16,I5)") EXP(FLSI+I_P*FLRI),FLRI*I_P,FLSI,GSEN,IDEG
             ENDDO
!C             CLOSE(17)
             CLOSE(iunit)
          ENDIF
        
          IF(TMONTE.AND.TENERGY.AND.NTAY(1).EQ.-1) THEN
             WRITE(6,*) "Calculating Exact MC Energy..."
             EN=DOEXMC(NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,0.0_dp,IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS)
          ENDIF
          IF(TBEGRAPH) THEN
             STOP 'BEGRAPH not implemented'
             IF(TENERGY) THEN
                IF(NTAY(1).NE.0) THEN
!                   CALL DOBEGRAPH(NDET,NEVAL,CK,W,I_P,ILOGGING,G1,NMRKS,nEl,NBASISMAX,nBasis,BRR)
                ELSE
!C..     NTAY=0 signifying we're going to calculate the RHO values when we
!C..     need them from the list of eigenvalues.   
!C..     Hide NMSH=NEVAL
!C..          FCK=W
!C..          ZIA=CK
!C..          UMAT=NDET
!C..          ALAT=NMRKS        
!                   CALL DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,             &
!     &                NBASISMAX,G1,nBasis,BRR,NEVAL,W,CK,NMAX,NMRKS,NDET,      &
!     &                NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
                ENDIF
             ELSE
!                CALL DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,                &
!     &                NBASISMAX,G1,nBasis,BRR,NMSH,FCK,NMAX,ALAT,UMAT,         &
!     &                NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
             ENDIF
          ENDIF
          IF(NTAY(1).EQ.0.AND.TENERGY) THEN
              WRITE(6,*) "Using exact RHOs generated on the fly"
!C.. we've calculated energies, and we're passing them through to
!C.. calculate the exact RHOS
!C.. NTAY=0 signifying we're going to calculate the RHO values when we
!C.. need them from the list of eigenvalues.   
!C.. Hide NMSH=NEVAL
!C..          FCK=W
!C..          ZIA=CK
!C..          UMAT=NDET
!C..          ALAT=NMRKS
!C..          NMAX=ARR
!              CALL CALCRHOPII2(BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,                        &
!     &               NBASISMAX,G1,nBasis,BRR,NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,           &
!     &                RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,           &
!     &                DETINV,TSPECDET,SPECDET)
!a
              call stop_all("DoExactVertexCalc","DoExactVertexCalc non-functional.")
          endif
            
        End Subroutine DoExactVertexCalc

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
         use UMatCache , only : TSTARSTORE
         use CalcData , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TMCDIRECTSUM,g_Multiweight,G_VMC_FAC,TMPTHEORY
         use CalcData, only : STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph,TStarTrips,THDiag,TMCStar,TFCIMC,TMCDets,tCCMC
         use CalcData , only : TRhoElems,TReturnPathMC, tUseProcsAsNodes,tRPA_QBA, tDetermProj, tFTLM, tSpecLanc
         use CalcData, only: tExactSpec, tExactDiagAllSym
         use RPA_Mod, only : tDirectRPA
         use CCMCData, only: tExactCluster,tCCMCFCI,tAmplitudes,tExactSpawn,tCCBuffer,tCCNoCuml
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
               case("CCMC")
                  !Piggy-back on the FCIMC code
                  I_HMAX=-21
                  TFCIMC=.true.
                  TCCMC=.true.
                  tExactCluster=.false.
                  tExactSpawn=.false.
                  tCCMCFCI=.false.
                  tAmplitudes=.false.
                  tCCBuffer=.false.
                  tCCNoCuml=.false.
                  do while(item.lt.nitems)
                    call readu(w)
                    select case(w)
                    case("PARTICLE")
                       tFCIMC=.false.  !We don't use the FCIMC code, but use a standalone.
                       tAmplitudes=.false.
                    case("AMPLITUDE")
                       tAmplitudes=.true.
                       tFCIMC=.false.  !We don't use the FCIMC code, but use a standalone.
                       tCalcFCIMCPsi=.true. !We want a full list of dets
                    case("EXACTCLUSTER")
                       tExactCluster=.true.
                    case("EXACTSPAWN")
                       tExactSpawn=.true.
                    case("FCI")
                       tCCMCFCI=.true.
                    case("BUFFER")
                        tCCBuffer=.true.
                    case("NOCUML")
                       tCCNoCuml=.true.
                    case default
                       call report("Keyword error with "//trim(w),.true.)
                    end select
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
                        IF(TSTARSTORE) call report("Error - "  &
     &                   //"can only use STARSTOREREAD with "  &
     &                   //"double excitations of HF",.true.)
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
         HElement_t BP, tmp
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
