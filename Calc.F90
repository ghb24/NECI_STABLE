MODULE Calc
        
        use CalcData

        IMPLICIT NONE

        contains

        subroutine SetCalcDefaults()
          != Set defaults for Calc data items.

          Use Determinants, only : iActiveBasis, SpecDet, tSpecDet, nActiveSpace
          Use DetCalc, only: iObs, jObs, kObs, tCorr, B2L, tRhoOfR, tFodM, DETINV
          Use DetCalc, only: icilevel, nBlk, nCycle, nEval, nKry, tBlock, tCalcHMat
          Use DetCalc, only: tEnergy, tRead
          use IntegralsData, only: tNeedsVirts
          use SystemData, only : Beta,nEl
          use default_sets
          implicit none

!       Values for old parameters.
!       These have no input options to change the defaults, but are used in the code.
          TRHOOFR = .false.
          TCORR = .false.
          TFODM = .false.
          B2L = 1.D-13
          TMC = .false.
          NHISTBOXES = 0
          TREADRHO = .false.
          TRHOIJ = .false.
          TBEGRAPH = .false.


!       Calc defaults 
          tSymmetricField=.false.
          NoMagDets=1
          BField=0.D0
          tMagnetize=.false.
          tFixShiftKii=.false.
          FixedKiiCutoff=0.D0
          tFixCASShift=.false.
          OccCASorbs=0
          VirtCASorbs=0
          tGlobalSftCng=.false.
          TLocalAnnihilation=.false.
          TDistAnnihil=.false.
          TAnnihilonproc=.false.
          TUnbiasPGeninProjE=.false.
          TFixShiftShell=.false.
          FixShift=0.D0
          ShellFix=0
          TRegenExcitgens=.false.
          MemoryFac=30.D0
          MemoryFacExcit=30.D0
          TStartSinglePart=.false.
          TFixParticleSign=.false.
          TProjEMP2=.false.
          TExcludeRandGuide=.false.
          THFRetBias=.false.
          TSignShift=.false.
          NEquilSteps=0
          RhoApp=10
          TResumFCIMC=.false.
          TRhoElems=.false.
          TReturnPathMC=.false.
          CLMax=NEl
          PRet=1.D0
          TMCDiffusion=.false.
          TNoAnnihil=.false.
          TNodalCutoff=.false.
          NodalCutoff=0.75
          TFullUnbias=.false.
          TExtraPartDiff=.false.
          TFlipTau=.false.
          FlipTauCyc=73     !A prime
          Lambda=0.D0
          TDiffuse=.false.
          TNoBirth=.false.
          GrowMaxFactor=20.D0
          CullFactor=2.D0
          TStartMP1=.false.
          TFCIMC=.false.
          TMCDets=.false.
          TBinCancel=.false.  
          ScaleWalkers=1.D0
          TReadPops=.false.
          StepsSft=100
          SftDamp=10.0
          Tau=0.D0
          InitWalkers=3000
          NMCyc=2000
          DiagSft=0.D0
          HApp=1
          TMCStar=.false.
          THDiag=.false.
          GrowGraphsExpo=2.D0
          TGrowInitGraph=.false.
          NoMCExcits=5000
          TMCExcitSpace=.false.
          TMaxExcit=.false.
          TFullDiag=.false.
          TSinglesExcitSpace=.false.
          TOneExcitConn=.false.
          TStarTrips=.false.
          TLanczos=.false.
          TNoSameExcit=.false.
          TInitStar=.false.
          NoMoveDets=1
          TMoveDets=.false.
          GraphBias=0.99
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
          TVARCALC(:)=.false.              !This is whether to calculate the expected variance for a MC run when doing full sum (seperate denominator and numerator at present)
          TBIN=.false.
          TVVDISALLOW=.false.
          tMCDirectSum=.FALSE.
          TMPTHEORY=.FALSE.
          tMP2Standalone=.FALSE.
          TMODMPTHEORY=.FALSE.
          G_VMC_PI = 0.95
          G_VMC_SEED = -7
          G_VMC_FAC = 16
          TUPOWER=.false.
          G_VMC_EXCITWEIGHT(:)=0.D0
          G_VMC_EXCITWEIGHTS(:,:)=0.D0
          EXCITFUNCS(:)=.false.
          EXCITFUNCS(10)=.true.
          NPaths=0
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
          BETA=0.D0
          BETAP=1.D-4
          TBETAP=.false.
          RHOEPSILON=1.D-6
!DBETA now has  three elements
!          DBETA(1) is DBETA
!          DBETA(2) is GRAPHEPSILON
!          DBETA(3) is PGENEPSILON
          DBETA(1)=-1.D0
          DBETA(2)=0.D0
          DBETA(3)=0.D0
          StarConv=1.d-3
          calcp_sub2vstar=.false.
          calcp_logweight=.false.
          TENPT=.false.
          TLADDER=.false. 

          tNeedsVirts=.true.! Set if we need virtual orbitals  (usually set).  Will be unset (by Calc readinput) if I_VMAX=1 and TENERGY is false

          lNoTriples=.false.

!Feb 08 default set.
          IF(Feb08) THEN
              RhoEpsilon=1.D-08
          ENDIF
      
        end subroutine SetCalcDefaults



        SUBROUTINE CalcReadInput()
          USE input
          Use Determinants, only : iActiveBasis, SpecDet, tagSpecDet, tSpecDet, nActiveSpace
          use SystemData, only : Beta,nEl
          Use DetCalc, only: iObs, jObs, kObs, tCorr, B2L, tRhoOfR, tFodM, DETINV
          Use DetCalc, only: icilevel, nBlk, nCycle, nEval, nKry, tBlock, tCalcHMat
          Use DetCalc, only: tEnergy, tRead
          use IntegralsData, only: tNeedsVirts,NFROZEN
          use UMatCache, only: gen2CPMDInts
          use global_utilities
          IMPLICIT NONE
          LOGICAL eof
          CHARACTER (LEN=100) w
          CHARACTER(*),PARAMETER :: t_r='CalcReadInput'
          INTEGER :: l,i,ierr

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
            case("LANCZOS")
!Sets the diagonaliser for the GraphMorph algorithm to be Lanczos
                TLanczos=.true.
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
            case("ENDCALC")
                exit calc
            case("METHODS")
                if(I_HMAX.ne.0) call report("METHOD already set",.true.)
                I_HMAX=-10
                I_VMAX=1
                methods: do
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
                      
                   case("ENDMETHODS")
                      exit methods
                   case default
                      call report ("Keyword "//trim(w)//" not recognized",.true.)
                   end select
                end do methods
            
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
!This excitation weighting involves a step function between the virtual and occupied electon manifold (i.e. step is at the chemical potential)
!When choosing an electron to move, the probability of selecting it is 1 if the electron is in the virtual manifold
!and (g_VMC_ExcitWeights(1,1) if in the virtual manifold. When choosing where to excite to, the situation is reversed, and the probability of selecting it is
!1 if the electron is in the occupied manifold and g_VMC_ExcitWeights(2,1) if in the occupied manifold. U-weighting is the third parameter as before.
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
               if (DBETA(1) .lt. 0 ) then
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
                call getf(DBETA(1))
            case("RHOEPSILON")
                call getf(RHOEPSILON)
            case("GRAPHEPSILON")
                call getf(DBETA(2))
            case("PGENEPSILON")
                call getf(DBETA(3))
!This indicates the number of times the eigenvalues of the star matrix should be evaluated to achieve the linear approximation when STARSTARS set,
            case("LINEPOINTSSTAR")
                call geti(LinePoints)
!This is the number of vertices in the Graph Morph graph. Alternativly, it is used by ResumFCIMC, as the size of their graphs. Then, if it is negative, the graph is all possible connections
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
!This means that the space accessible to the morphing algorithm is the space of single excitations of the determinants in the graph.
                TSinglesExcitSpace=.true.
            case("FULLDIAGTRIPS")
!When constructing a star of triples from each double star, then this tag results in a full diagonalisation of this matrix.
                TFullDiag=.true.
            case("MCEXCITSPACE")
!In GraphMorph, this means that the space of excitations is chosen randomly
!It is also an option in FCIMC, where it indicates the number of excitations to be chosen randomly from each chosen walker
                TMCExcitSpace=.true.
                call geti(NoMCExcits)
            case("GROWINITGRAPH")
!In GraphMorph, this means that the initial graph is grown non-stochastically from the excitations of consecutive determinants
                TGrowInitGraph=.true.
            case("GROWGRAPHSEXPO")
!In GraphMorph, this is the exponent to which the components of the excitation vector and eigenvector will be raised to turn them into probabilities.
                call getf(GrowGraphsExpo)
            case("HAPP")
!For graph MC, this indicates the number of local applications of the hamiltonian to random determinants before the trial eigenvector is updated
                call geti(HApp)
            case("MAXEXCIT")
!This imposes a maximum excitation level to the space that GraphMorph can explore. Note: FCIMC uses EXCIT to indicate a maximum excit level.
                TMaxExcit=.true.
                call geti(iMaxExcitLevel)
            case("INITWALKERS")
!For FCIMC, this is the number of walkers to start with
                call geti(InitWalkers)
            case("NMCYC")
!For FCIMC, this is the number of MC cycles to perform
                call geti(NMCyc)
            case("DIAGSHIFT")
!For FCIMC, this is the amount extra the diagonal elements will be shifted. This is proportional to the deathrate of walkers on the determinant
                call getf(DiagSft)
            case("TAU")
!For FCIMC, this can be considered the timestep of the simulation. It is a constant which will increase/decrease the rate of spawning/death for a given iteration.
                call getf(Tau)
            case("SHIFTDAMP")
!For FCIMC, this is the damping parameter with respect to the update in the DiagSft value for a given number of MC cycles.
                call getf(SftDamp)
            case("STEPSSHIFT")
!For FCIMC, this is the number of steps taken before the Diag shift is updated
                call geti(StepsSft)
            case("READPOPS")
!For FCIMC, this indicates that the initial walker configuration will be read in from the file POPSFILE, which must be present.
!DiagSft and InitWalkers will be overwritten with the values in that file.
                TReadPops=.true.
            case("SCALEWALKERS")
!For FCIMC, if this is a way to scale up the number of walkers, after having read in a POPSFILE
                call getf(ScaleWalkers)
            case("BINCANCEL")
!This is a seperate method to cancel down to find the residual walkers from a list, involving binning the walkers into their determinants. This has to refer to the whole space, and so is much slower.
                TBinCancel=.true.
            case("STARTMP1")
!For FCIMC, this has an initial configuration of walkers which is proportional to the MP1 wavefunction
                TStartMP1=.true.
            case("GROWMAXFACTOR")
!For FCIMC, this is the factor to which the initial number of particles is allowed to go before it is culled
                call getf(GrowMaxFactor)
            case("CULLFACTOR")
!For FCIMC, this is the factor to which the total number of particles is reduced once it reaches the GrowMaxFactor limit
                call getf(CullFactor)
            case("EQUILSTEPS")
!For FCIMC, this indicates the number of cycles which have to
!pass before the energy of the system from the doubles
!population is counted
                call geti(NEquilSteps)
            case("NOBIRTH")
!For FCIMC, this means that the off-diagonal matrix elements become zero, and so all we get is an exponential decay of the initial populations on the determinants, at a rate which can be exactly calculated and compared against.
                TNoBirth=.true.
            case("MCDIFFUSE")
                TDiffuse=.true.
!Lambda indicates the amount of diffusion compared to spawning in the FCIMC algorithm.
                call getf(Lambda)
            case("FLIPTAU")
!This indicates that time is to be reversed every FlipTauCyc cycles in the FCIMC algorithm. This might help with undersampling problems.
                TFlipTau=.true.
                call geti(FlipTauCyc)
            case("NON-PARTCONSDIFF")
!This is a seperate partitioning of the diffusion matrices in FCIMC in which the antidiffusion matrix (+ve connections) create a net increase of two particles.
                TExtraPartDiff=.true.
            case("FULLUNBIASDIFF")
!This is for FCIMC, and fully unbiases for the diffusion process by summing over all connections
                TFullUnbias=.true.
            case("NODALCUTOFF")
!This is for all types of FCIMC, and constrains a determinant to be of the same sign as the MP1 wavefunction at that determinant, if the normalised component of the MP1 wavefunction is greater than the NodalCutoff value.
                TNodalCutoff=.true.
                call getf(NodalCutoff)
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
!This is for resummed FCIMC, it indicates the number of propagation steps around each subgraph before particles are assigned to the nodes
                call geti(RhoApp)
            case("SIGNSHIFT")
!This is for FCIMC and involves calculating the change in shift depending on the absolute value of the sum of the signs of the walkers.
!This should hopefully mean that annihilation is implicitly taken into account.
                TSignShift=.true.
            case("HFRETBIAS")
!This is a simple guiding function for FCIMC - if we are at a double excitation, then we return to the HF determinant with a probability PRet.
!This is unbiased by the acceptance probability of returning to HF.
                THFRetBias=.true.
                call getf(PRet)
            case("EXCLUDERANDGUIDE")
!This is an alternative method to unbias for the HFRetBias. It invloves disallowing random excitations back to the guiding function (HF Determinant)
                TExcludeRandGuide=.true.
            case("PROJECTE-MP2")
!This will find the energy by projection of the configuration of walkers onto the MP2 wavefunction.
                TProjEMP2=.true.
            case("FIXPARTICLESIGN")
!This uses a modified hamiltonian, whereby all the positive off-diagonal hamiltonian matrix elements are zero. Instead, their diagonals are modified to change the
!on-site death rate. Particles now have a fixed (positive) sign which cannot be changed and so no annihilation occurs.
                TFixParticleSign=.true.
            case("STARTSINGLEPART")
!A FCIMC option - this will start the simulation with a single positive particle at the HF, and fix the shift at its initial value, until the number of particles gets to the INITPARTICLES value.
                TStartSinglePart=.true.
            case("MEMORYFAC")
!An FCIMC option - MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
                CALL Getf(MemoryFac)
            case("MEMORYFACEXCIT")
!!An FCIMC option - MemoryFac is the factor by which space will be made available for excitation generators compared to InitWalkers. This can be smaller than
!memoryfac, because the excitation generator array is not used in the parallel annihilation, which may not be exactly load-balanced because of differences in the wavevector.
                CALL Getf(MemoryFacExcit)
            case("REGENEXCITGENS")
!An FCIMC option. With this, the excitation generators for the walkers will NOT be stored, and regenerated each time. This will be slower, but save on memory.
                TRegenExcitGens=.true.
            case("FIXSHELLSHIFT")
!An FCIMC option. With this, the shift is fixed at a value given here, but only for excitations which are less than <ShellFix>. This will almost definitly give the wrong answers for both the energy
!and the shift, but may be of use in equilibration steps to maintain particle density at low excitations, before writing out the data and letting the shift change.
                TFixShiftShell=.true.
                CALL Geti(ShellFix)
                CALL Getf(FixShift)
            case("FIXKIISHIFT")
!A Parallel FCIMC option. Similar to FixShellShift option, but will fix the shifts of the particles which have a diagonal
!matrix element Kii of less than the cutoff, FixedKiiCutOff.
                TFixShiftKii=.true.
                CALL Getf(FixedKiiCutoff)
                CALL Getf(FixShift)
            
            case("FIXCASSHIFT")                
!A Parallel FCIMC option similar to the FixShellShift and FixShiftKii options.  In this option, an active space is chosen containing a certain number of highest occupied spin orbitals (OccCASorbs) and
!lowest unoccupied spin orbitals (VirtCASorbs).  The shift is then fixed only for determinants which have completely occupied spin orbitals for those lower in energy than the active space, and completely unoccupied spin orbitals above the active space.  i.e. the electrons are only excited within the active space.  
                TFixCASShift=.true.
                call Geti(OccCASorbs)
                call Geti(VirtCASorbs)
                call Getf(FixShift)

            case("UNBIASPGENINPROJE")
!A FCIMC serial option. With this, walkers will be accepted with probability tau*hij. i.e. they will not unbias for PGen in the acceptance criteria, but in the term for the projected energy.
                TUnbiasPGeninProjE=.true.
            case("ANNIHILATEONPROCS")
!A parallel FCIMC option. With this, walkers will only be annihilated with other walkers on the same processor. 
                TAnnihilonproc=.true.
            case("ANNIHILATDISTANCE")
!A Serial FCIMC experimental option. With this, walkers have the ability to annihilate each other as long as they are connected, which they will do with probability = Lambda*Hij
                TDistAnnihil=.true.
                call Getf(Lambda)
            case("LOCALANNIHIL")
!A parallel FCIMC experimental option. This will attempt to compensate for undersampled systems, by including extra annihilation for walkers which are the sole occupier of determiants
!This annihilation is governed by the parameter Lambda, which is also used in other circumstances as a variable, but should not be used at the same time.
                TLocalAnnihilation=.true.
                call Getf(Lambda)
            case("GLOBALSHIFT")
!A parallel FCIMC option. It is generally recommended to have this option on. This will calculate the growth rate of the system as a simple ratio of the total walkers on all processors
!before and after update cycle. This however is incompatable with culling, and so is removed for update cycles with this in. 
                tGlobalSftCng=.true.
            case("MAGNETIZE")
!This is a parallel FCIMC option. It chooses the largest weighted MP1 components and records their sign. If then a particle occupies this determinant and is of the opposite sign, it energy,
!i.e. diagonal matrix element is raised by an energy given by BField.
                tMagnetize=.true.
                tSymmetricField=.false.
                call Geti(NoMagDets)
                call Getf(BField)
            case("MAGNETIZESYM")
!A parallel FCIMC option. Similar to the MAGNETIZE option, but in addition to the energy being raised for particles of the opposite sign, the energy is lowered by the same amount for particles
!of 'parallel' sign.
                call Geti(NoMagDets)
                call Getf(BField)
                tSymmetricField=.true.
                tMagnetize=.true.
            case default
                call report("Keyword "                                &
     &            //trim(w)//" not recognized in CALC block",.true.)
            end select
          end do calc
          IF((.not.TReadPops).and.(ScaleWalkers.ne.1.D0)) THEN
              call report("Can only specify to scale walkers if READPOPS is set",.true.)
          ENDIF
          IF(tFixShiftKii.and.tFixShiftShell) THEN
              call Stop_All("ReadCalc","Cannot have both fixshellshift and fixshiftkii options")
          ENDIF
          IF(tFixShiftKii.and.tFixCASShift) THEN
              call Stop_All("ReadCalc","Cannot have both fixshiftkii and fixCASshift options")
          ENDIF
          IF(tFixCASShift.and.tFixShiftShell) THEN
              call Stop_All("ReadCalc","Cannot have both fixCASshift and fixshiftshell options")
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
          use SystemData, only: G1, Alat, Beta, BRR, ECore, LMS, nBasis, nBasisMax, STot,tCSF,nMsh,nEl
          use SystemData, only: tUEG
          use IntegralsData, only: FCK, CST, nMax, UMat
          use IntegralsData, only: HFEDelta, HFMix, NHFIt, tHFCalc
          Use Determinants, only: FDet, tSpecDet, SpecDet, GetHElement2
          Use DetCalc, only: DetInv, nDet, tRead
          use global_utilities
          
          REAL*8 CalcT, CalcT2, GetRhoEps
          
          
          INTEGER I, IC
          INTEGER nList
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

          IF(TCSF.AND.TSPECDET) THEN
             WRITE(6,*) "TSPECDET set.  SPECDET is"
             CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
             CALL NECI_ICOPY(NEL,SPECDET,1,FDET,1)
             CALL GETCSFFROMDET(FDET,SPECDET,NEL,STOT,LMS)
             WRITE(6,*) "CSF with 2S=",STOT," and 2Sz=",LMS," now in SPECDET is"
             CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
          ENDIF
          IF(TSPECDET.AND.SPECDET(1).EQ.0) THEN
             WRITE(6,*) "TSPECDET set, but invalid.  using FDET"
             CALL NECI_ICOPY(NEL,FDET,1,SPECDET,1)
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
             WRITE(6,*) ' BETAP=',BETAP
             WRITE(6,*) ' RESETTING P '
             IF(I_P.GT.100000) WRITE(6,*) ' *** WARNING I_P=',I_P
          ENDIF

          WRITE(6,*) ' BETA, P :',BETA,I_P
       
!C         DBRAT=0.001
!C         DBETA=DBRAT*BETA
          WRITE(6,*) "DBETA=",DBETA

          IF(.NOT.TREAD) THEN
!             CALL WRITETMAT(NBASIS)
             IC=0
             WRITE(6,*) '<D0|H|D0>=',GETHELEMENT2(FDET,FDET,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
             WRITE(6,*) '<D0|T|D0>=',CALCT(FDET,NEL,G1,NBASIS)
             IF(TUEG) THEN
!  The actual KE rather than the one-electron part of the Hamiltonian
                WRITE(6,*) 'Kinetic=',CALCT2(FDET,NEL,G1,ALAT,NBASIS,CST)
             ENDIF
          ENDIF
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
             CALL WRITEDET(6,MCDET,NEL,.TRUE.)
          ENDIF
!C.. we need to calculate a value for RHOEPS, so we approximate that
!C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this
!C.. we have put TMAT instead of ZIA
          IF(I_HMAX.NE.-20) THEN
!C.. If we're using rhos,
             RHOEPS=GETRHOEPS(RHOEPSILON,BETA,NEL,NBASISMAX,G1,nBasis,BRR, NMSH,FCK,NMAX,ALAT,UMAT,I_P,ECORE)

             WRITE(6,*) "RHOEPS:",RHOEPS
          ELSE
!C.. we're acutally diagonalizing H's, so we just leave RHOEPS as RHOEPSILON
             RHOEPS=RHOEPSILON
          ENDIF

        End Subroutine CalcInit
    
    
    
        Subroutine CalcDoCalc()
          use SystemData, only: Alat, Arr,Brr, Beta, ECore, G1, LMS, LMS2, nBasis,NMSH, nBasisMax
          use SystemData, only: SymRestrict, tCSF, tParity, tSpn, ALat, Beta
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB,nEl
          Use DetCalc, only : CK, DetInv, nDet, nEval, tEnergy, tRead, nmrks, w
          Use Determinants, only: FDet, nActiveBasis, SpecDet, tSpecDet
          use IntegralsData, only: FCK, NMAX, UMat, FCK
          use IntegralsData, only: HFEDelta, HFMix,nTay
          Use Logging, only: iLogging
          use Parallel_Calc
!          Use MCDets, only: MCDetsCalc
!Calls
          REAL*8 DMonteCarlo2
!Local Vars
          REAL*8 EN, ExEn, GsEN
          REAL*8 FLRI, FLSI
          REAL*8 RH
          LOGICAL tWarn
          integer iSeed
          iSeed=7 

!          call MCDetsCalc(FDet, iSeed, nwhtay(1))
!          stop
    
!C.. we need to calculate a value for RHOEPS, so we approximate that
!C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this 
!C.. If we haven't already calced RHOEPS, do it now
          Call DoExactVertexCalc()

          IF (tMP2Standalone) then
              call ParMP2(FDet)
! Parallal 2v sum currently for testing only.
!          call Par2vSum(FDet)
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
             WRITE(6,*) "Calculating ",NPATHS," W_Is..."
             IF(BTEST(ILOGGING,1)) THEN
                IF(I_HMAX.EQ.-10) THEN
                   OPEN(11,FILE="MCSUMMARY",STATUS="UNKNOWN")
                   WRITE(11,*) "Calculating ",NPATHS," W_Is..."
                   CLOSE(11)
                ELSE
                   OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
                   WRITE(11,*) "Calculating ",NPATHS," W_Is..."
                   CLOSE(11)
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
                     CALL CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,         &
     &                 NBASISMAX,G1,nBasis,BRR,NMSH,FCK,ARR,ALAT,UMAT,NTAY,          &
     &                 RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,           &
     &                 DETINV,TSPECDET,SPECDET)
                ELSE
                   IF(TCSF) THEN
                      IF(.NOT.TSPECDET) THEN
                         WRITE(6,*) "SPECDET not specified. Using Fermi determinant ONLY"
                         TSPECDET=.TRUE.
                         CALL NECI_ICOPY(NEL,FDET,1,SPECDET,1)
                      ENDIF
                   ENDIF
!C.. Instead of NMAX we have ARR
                  IF(TPARITY) THEN
                      WRITE(6,*) "Using symmetry restriction:"
                      CALL WRITEALLSYM(6,SymRestrict,.TRUE.)
                  ENDIF
                  IF(TSPN) THEN
                      WRITE(6,*) "Using spin restriction:",LMS
                  ENDIF
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
          IF(TMONTE.and..not.tMP2Standalone) THEN
!             DBRAT=0.01
!             DBETA=DBRAT*BETA
             WRITE(6,*) "I_HMAX:",I_HMAX
             WRITE(6,*) "Calculating MC Energy..."
             CALL FLUSH(6)
             IF(NTAY(1).GT.0) THEN
                WRITE(6,*) "Using approx RHOs generated on the fly, NTAY=",NTAY(1)
!C.. NMAX is now ARR
                EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,IMCSTEPS,G1,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS,NMSH,FCK,ARR,ALAT,UMAT,NTAY,RHOEPS,NWHTAY,ILOGGING,ECORE,BETAEQ) 
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
                   EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,IMCSTEPS,             &
     &                G1,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS,                                 &
     &                NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,RHOEPS,NWHTAY,ILOGGING,ECORE,BETAEQ)
                ELSE
                   STOP "TENERGY not set, but NTAY=0" 
                ENDIF
             ENDIF
             WRITE(6,*) "MC Energy:",EN
!CC           WRITE(12,*) DBRAT,EN
          ENDIF
         
!C.. /AJWT
        End Subroutine

        Subroutine DoExactVertexCalc()
          use SystemData, only: Alat, Beta, Brr, ECORE, G1, nBasis, nBasisMax,nMsh, Arr,nEl
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB
          use IntegralsData, only: fck, nMax, UMat,nTay
          Use DetCalc, only: cK, nDet, nEval, tEnergy, tRead, W, NMRKS, DetInv
          Use Determinants, only: specdet, tSpecDet
          Use Logging, only: iLogging
          real*8 flri, flsi
          REAL*8 En, ExEn, GSEn
          REAL*8 RH
          INTEGER iDeg, III
          Type(BasisFN) iSym
          LOGICAL tWarn
          
          REAL*8 CalcMCEn, CalcDLWDB, DoExMC
            
          IF(TENERGY) THEN
             RHOEPS=RHOEPSILON*EXP(-BETA*(W(1))/I_P)
            WRITE(6,*) "RHOEPS:",RHOEPS
             IF(TREAD) THEN
                EXEN=CALCMCEN(NDET,NEVAL,CK,W,BETA,0.D0)
                WRITE(6,"(A,F19.5)") "EXACT E(BETA)=",EXEN
                GSEN=CALCDLWDB(1,NDET,NEVAL,CK,W,BETA,0.D0)
                WRITE(6,"(A,F19.5)") "EXACT DLWDB(D0)=",GSEN
             ENDIF
             OPEN(14,FILE='RHOPIIex',STATUS='UNKNOWN')
             IF(NDETWORK.EQ.0.OR.NDETWORK.GT.NDET) NDETWORK=NDET
             DO III=1,NDETWORK
             
                CALL CALCRHOPII(III,NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,0.D0,FLRI,FLSI,TWARN)
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
                   CALL CALCRHO2(NMRKS(1,III),NMRKS(1,III),BETA,I_P,NEL, NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,1,0,ECORE)
!C                   WRITE(6,*) RH
                   FLRI=LOG(RH)
                   FLSI=FLSI-I_P*FLRI
                   ENDIF
                ENDIF
                CALL WRITEDET(14,NMRKS(1,III),NEL,.FALSE.)
                GSEN=CALCDLWDB(III,NDET,NEVAL,CK,W,BETA,0.D0)
                CALL GETSYM(NMRKS(1,III),NEL,G1,NBASISMAX,ISYM)
                CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
                WRITE(14,"(4G25.16,I5)") EXP(FLSI+I_P*FLRI),FLRI*I_P,FLSI,GSEN,IDEG
             ENDDO
!C             CLOSE(17)
             CLOSE(14)
          ENDIF
        
          IF(TMONTE.AND.TENERGY.AND.NTAY(1).EQ.-1) THEN
             WRITE(6,*) "Calculating Exact MC Energy..."
             EN=DOEXMC(NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,0.D0,IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS)
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
              CALL CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,           &
     &               NBASISMAX,G1,nBasis,BRR,NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,           &
     &                RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,           &
     &                DETINV,TSPECDET,SPECDET)

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
         use input
         use UMatCache , only : TSTARSTORE
         use CalcData , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TMCDIRECTSUM,g_Multiweight,G_VMC_FAC,TMPTHEORY
         use CalcData, only : STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph,TStarTrips,THDiag,TMCStar,TFCIMC,TMCDets,TMCDiffusion
         use CalcData , only : TRhoElems,TReturnPathMC, TResumFCIMC
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
                   call readu(w)
                   select case(w)
                   case("MCDIFFUSION")
                       TMCDiffusion=.true.
                   case("RESUMFCIMC")
                       TResumFCIMC=.true.
                   endselect
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
                  IF(I_V.GT.0) g_MultiWeight(I_V)=1.D0
               case("MCDIRECT")
                  I_HMAX = -7
                  tMCDirectSum=.TRUE.
                   call readu(w)
                   select case(w)
                   case("HDIAG")
                       I_HMAX = -19
                   end select
                  G_VMC_FAC=0.D0
               case("MCMP")
                  tMCDirectSum=.TRUE.
                  G_VMC_FAC=0.D0
                  TMPTHEORY=.TRUE.
               case("GRAPHMORPH")
                   TGraphMorph=.true.
                   I_HMAX=-21
                   call readu(w)
                   select case(w)
                   case("HDIAG")
                       !If this is true, then it uses the hamiltonian matrix to determinant coupling to excitations, and to diagonalise to calculate the energy
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
         use input
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

! Calculate 'PATHS' using on-the-fly generated determinants
!  Determinants are generated by GENNEXTDET according to the various symmetry specifications given by the use
!  nActiveBasis will restrict the space of generated determinants to excitations from and to a given set of
!  basis functions
      SUBROUTINE CALCRHOPII3(BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,DETINV,TSPN,LMS,TPARITY,SymRestrict,     &
     &            TSPECDET,SPECDET,nActiveBasis)
         USE HElem
         use global_utilities
         use SystemData, only: BasisFN,BasisFNSize
         IMPLICIT NONE
         include 'irat.inc'
         INTEGER I_HMAX,NEL,NBASIS,I_VMAX
         INTEGER,ALLOCATABLE :: LSTE(:) !(NEL,NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)??!!
         INTEGER,ALLOCATABLE :: ICE(:)  !(NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)??!!
         TYPE(HElement)  UMAT(*)
         TYPE(HElement),allocatable  :: RIJLIST(:)
         integer,save :: tagRIJList=0,tagLSTE=0,tagICE=0
         REAL*8 BETA,FCK(*),ALAT(*),RHOEPS
         INTEGER NPATHS,NI(NEL),I_P,nBasisMax(5,*)
         INCLUDE 'gndwork.inc'
         INTEGER Work(GNDWorkSize+2*NEL)
         TYPE(BASISFN) G1(NBASIS)
         INTEGER BRR(NBASIS),NMSH,NMAX(*),NTAY,ILOGGING
         INTEGER III,NWHTAY,I,IMAX,ILMAX,LMS
         TYPE(BasisFN) ISYM,SymRestrict
         LOGICAL TSPN,TPARITY,TSYM
         REAL*8 DBETA,ECORE
         TYPE(HDElement) WLRI,WLSI,WLRI1,WLRI2,WLSI1,WLSI2,WI,DLWDB
         TYPE(HDElement) TOT,WLRI0,WLSI0,WINORM,HElP,NORM
         LOGICAL TNPDERIV,TDONE,TFIRST
         INTEGER DETINV
         INTEGER ISTART,IEND,IDEG
         LOGICAL TSPECDET,TLOG
         INTEGER SPECDET(NEL)
         TYPE(HDElement) DLWDB2,DLWDB3,DLWDB4,TOT2
         INTEGER nActiveBasis(2)
         type(timer), save :: proc_timer
         character(len=*), parameter :: thisroutine='CALCRHOPII3'
         TLOG=BTEST(ILOGGING,1)
         HElP=HDElement(I_P)
         TSYM=.TRUE.
         TOT=0.D0
         TOT2=0.D0
         NORM=0.D0
         IMAX=I_HMAX
         IF(I_VMAX.GT.IMAX) IMAX=I_VMAX
         proc_timer%timer_name=thisroutine
         call set_timer(proc_timer)
         WRITE(6,*) "Entering CALCRHOPII3..."
!         ILMAX=NDET
!.. We don't need to store lists for I_HMAX=-8
         ILMAX=(NBASIS-NEL)**2*NEL*NEL/4
         IF((I_HMAX.GE.-10.AND.I_HMAX.LE.-7).OR.I_HMAX.LE.-12) ILMAX=1
         allocate(LSTE((1+ILMAX)*NEL*IMAX))
         call LogMemAlloc('LSTE',size(LSTE),4/IRAT,thisroutine,tagLSTE)
         allocate(ICE((1+ILMAX)*IMAX))
         call LogMemAlloc('ICE',size(ICE),4/IRAT,thisroutine,tagICE)
         allocate(RIJList((1+ILMAX)*IMAX*2))
         call LogMemAlloc('RIJList',(1+ILMAX)*IMAX*2,8,thisroutine, tagRIJList)
!:         CALL PRINT_MEMORY()
         IF(I_VMAX.NE.0) THEN
            WRITE(6,*) "Using Vertex approximation.  I_VMAX=",I_VMAX
            IF(I_HMAX.EQ.0) WRITE(6,*) "I_HMAX=0.  Summing all I_HMAX up to P using contour"
            IF(I_HMAX.GT.0) WRITE(6,*) "I_HMAX=",I_HMAX
         ELSEIF(I_HMAX.NE.0) THEN
            WRITE(6,*) "Using hop-restricted paths. I_HMAX:",I_HMAX
         ELSE
            WRITE(6,*) "I_HMAX=I_VMAX=0. Using rho diagonalisation."
         ENDIF
         IF(TLOG) THEN
            IF(I_HMAX.EQ.-10) THEN
               OPEN(11,FILE="MCSUMMARY",STATUS="UNKNOWN")
               WRITE(11,*) "Calculating ",NPATHS," W_Is..."
               CLOSE(11)
            ELSE
               OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
               WRITE(11,*) "Calculating ",NPATHS," W_Is..."
               CLOSE(11)
            ENDIF
            OPEN(15,FILE='RHOPII',STATUS='UNKNOWN')
         ENDIF
         IF(DETINV.NE.0) THEN
            ISTART=ABS(DETINV)
            IEND=ABS(DETINV)
         ELSEIF(TSPECDET) THEN
            WRITE(6,*) "Calculating vertex series for specific det:"
            CALL WRITEDET(6,SPECDET,NEL,.TRUE.) 
            ISTART=-1
            IEND=1
         ELSE
            ISTART=1
            IEND=NPATHS
         ENDIF

         III=0
         TDONE=.FALSE.
         TFIRST=.TRUE.
         IF(.NOT.TSPECDET) THEN
            CALL GENNEXTDET(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,TPARITY,SymRestrict,ISYM,NI,.TRUE.,TDONE,WORK,nActiveBasis)
         ENDIF
         DO WHILE(III.NE.IEND.AND..NOT.TDONE)
          III=III+1
          IF(TSPECDET) THEN
             TDONE=.FALSE.
             CALL NECI_ICOPY(NEL,SPECDET,1,NI,1)
             IDEG=1
          ELSE
             CALL GENNEXTDET(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,TPARITY,SymRestrict,ISYM,NI,.FALSE.,TDONE,WORK,nActiveBasis)
             CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
          ENDIF
          IF(III.GE.ISTART.AND..NOT.TDONE) THEN
            IF(NPATHS.EQ.1.AND..NOT.TSPECDET) CALL WRITEDET(6,NI,NEL,.TRUE.) 
            CALL MCPATHSR3(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
     &         RHOEPS,LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,WLRI,WLSI,DBETA,DLWDB2)
            IF(TLOG) THEN
               WRITE(15,"(I12)",advance='no') III
               CALL WRITEDET(15,NI,NEL,.FALSE.)
               WRITE(15,"(3G25.16)",advance='no') EXP(WLSI+HElP*WLRI),WLRI*HElP,WLSI
            ENDIF
            IF(TFIRST) THEN
               TFIRST=.FALSE.
               WLRI0=WLRI
               WLSI0=WLSI
            ENDIF  
            IF(TNPDERIV) THEN
!.. if we're calculating the derivatives too
               CALL MCPATHSR3(NI,BETA+DBETA,I_P,I_HMAX, I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,  &
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY, ILOGGING,ECORE,ILMAX,WLRI1,WLSI1,DBETA,DLWDB3)
               CALL MCPATHSR3(NI,BETA-DBETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,   &
     &            NMAX,ALAT,UMAT,NTAY,RHOEPS,LSTE,ICE,RIJLIST,NWHTAY, ILOGGING,ECORE,ILMAX,WLRI2,WLSI2,DBETA,DLWDB4)
               DLWDB=-(HElP*(WLRI1-WLRI2)+(WLSI1-WLSI2))/HDElement(2*DBETA)
            ELSE
               DLWDB=DLWDB2
            ENDIF
!.. we calculate the energy with weightings normalized to the weight of
!.. the Fermi determinant, otherwise the numbers blow up
            WINORM=EXP(HElP*(WLRI-WLRI0)+(WLSI-WLSI0))
            NORM=NORM+HDElement(IDEG)*WINORM
            TOT=TOT+HDElement(IDEG)*WINORM*DLWDB
            IF(TLOG) WRITE(15,"(G25.16,I5)") DLWDB,IDEG
            IF(DETINV.EQ.III) THEN
               IF(TLOG) CALL FLUSH(15)
               WRITE(6,*) "Investigating det ",DETINV
               CALL FLUSH(6)
               CALL WIRD_SUBSET(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
     &            RHOEPS,ILOGGING,TSYM,ECORE)
            ENDIF
           ELSE
! Correct for overcounting
               III=III-1
           ENDIF
          ENDDO
         IF(TLOG) CLOSE(15)
         IF(TFIRST) THEN
            WRITE(6,*) "*** NO determinants found to calculate***"
         ELSE
            WRITE(6,*) "Total ",III," determinants summed."
         ENDIF
         WRITE(6,*) "Summed approx E(Beta)=",TOT/NORM
         deallocate(RIJList,LSTE,ICE)
         call LogMemDealloc(thisroutine,tagRIJList)
         call LogMemDealloc(thisroutine,tagLSTE)
         call LogMemDealloc(thisroutine,tagICE)
         call halt_timer(proc_timer)
         RETURN
      END SUBROUTINE CALCRHOPII3  



! Given an input RHOEPSILON, create Fermi det D out of lowest orbitals and get RHOEPS (which is rhoepsilon * exp(-(beta/P)<D|H|D>
      REAL*8 FUNCTION GETRHOEPS(RHOEPSILON,BETA,NEL,NBASISMAX,G1,NHG, BRR,NMSH,FCK,NMAX,ALAT,UMAT,I_P,ECORE)
         Use Determinants, only: GetHElement2
         USE HElem
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I,nBasisMax(5,*),I_P
         INTEGER BRR(*),NMSH,NMAX,NHG
         COMPLEX*16 FCK(*)
         REAL*8 RHOEPSILON,BETA,ECORE,ALAT(*)
         TYPE(HElement) BP,UMat(*)
         TYPE(BasisFN) G1(*)
         DO I=1,NEL
            NI(I)=BRR(I)
         ENDDO
         CALL NECI_SORTI(NEL,NI)
         BP=HElement(-BETA/I_P)
         GETRHOEPS=DSQRT(SQ(HElement(RHOEPSILON)*EXP(BP*GETHELEMENT2(NI, &
     &      NI,NEL,NBASISMAX,G1,NHG,BRR,NMSH,FCK,NMAX,ALAT,UMAT,        &
     &         0,ECORE))))
         RETURN
      END FUNCTION GetRhoEps



! Calculate the kinetic energy of the UEG (this differs from CALCT by including the constant CST
      REAL*8 FUNCTION CALCT2(NI,NEL,G1,ALAT,NBASIS,CST)
         USE HElem
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,I,J
         TYPE(BasisFN) G1(*)
         REAL*8 ALAT(4),CST,TMAT
         LOGICAL ISCSF
         CALCT2=0.D0
         IF(ISCSF(NI,NEL)) RETURN
         DO J=1,NEL
            I=NI(J)
           TMAT=((ALAT(1)**2)*((G1(I)%K(1)**2)/(ALAT(1)**2)+   &
     &         (G1(I)%K(2)**2)/(ALAT(2)**2)+                   &
     &         (G1(I)%K(3)**2)/(ALAT(3)**2)))
           TMAT=TMAT*CST
           CALCT2=CALCT2+TMAT
         ENDDO
         RETURN
      END FUNCTION CALCT2
