
MODULE Calc
        USE input
        
        USE System , only : NEL,Feb08,defaults
        USE Integrals , only : NFROZEN
        Use Determinants, only :nActiveSpace
        
        IMPLICIT NONE
        save

        LOGICAL TSTAR,TTROT,TMCExcitSpace
        LOGICAL TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
        LOGICAL TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER
        LOGICAL EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET
        LOGICAL TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
        LOGICAL TLADDER,TMC,TREADRHO,TRHOIJ,TBiasing,TMoveDets
        LOGICAL TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph
        LOGICAL TInitStar,TNoCross,TNoSameExcit,TLanczos,TStarTrips
        LOGICAL TMaxExcit,TOneExcitConn,TSinglesExcitSpace,TFullDiag
        
        INTEGER NWHTAY(3,10),NPATHS,NoMoveDets,NoMCExcits
        INTEGER NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED
        INTEGER IMCSTEPS,IEQSTEPS,MDK(5),Iters,NDets
        INTEGER CUR_VERT,NHISTBOXES,I_P,LinePoints,iMaxExcitLevel
        
        
        REAL*8 g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
        REAL*8 G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10)
        REAL*8 BETAP,RHOEPSILON,DBETA(3),STARCONV,GraphBias



!// additional from NECI.F
        INTEGER, Allocatable :: MCDet(:)
        REAL*8 RHOEPS ! calculated from RHOEPSILON

!// set if we include no triple-excitations as the 3rd vertex in 3+ vertex graphs.
        LOGICAL lNoTriples
        contains

        SUBROUTINE CalcReadInput()
            Use Determinants, only : iActiveBasis, SpecDet, tSpecDet
            Use System, only : Beta
            Use DetCalc, only: iObs, jObs, kObs, tCorr, B2L, tRhoOfR, tFodM, DETINV
            Use DetCalc, only: icilevel, nBlk, nCycle, nEval, nKry, tBlock, tCalcHMat
            Use DetCalc, only: tEnergy, tRead
            Use Integrals, only: tNeedsVirts
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
        INTEGER :: l,i,ierr

!     Values for old parameters.
!     These have no input options to change the defaults, but are used in the code.
      NoMCExcits=5000
      TMCExcitSpace=.false.
      TMaxExcit=.false.
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
          TFullDiag=.false.
          TSinglesExcitSpace=.false.
          TOneExcitConn=.false.
          TStarTrips=.false.
          TLanczos=.false.
          TNoSameExcit=.false.
          TNoCross=.false.
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
!          IP_SPECDET=0
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

!Feb 08 defaults
          IF(Feb08) THEN
              RhoEpsilon=1.D-08
          ENDIF
         

      
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
     &                  " not recognised",.true.)
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
     &               I_VMAX)
                  case("EXCITATIONS")
                     call readu(w)
                     call inpgetexcitations(NWHTAY(2,I_VMAX),w)
                  case("CYCLES")
                     call readi(NWHTAY(2,I_VMAX))
                     if ( NWHTAY(1,I_VMAX).ne. -7.and.                  &
     &                    NWHTAY(1,I_VMAX).ne.-19 ) then
                        call report(trim(w)//" only valid for MC "      &
     &                   //"method",.true.)
                     end if
                  case("VERTICES")
                     call geti(NWHTAY(3,I_VMAX))
                  case("MULTIMCWEIGHT")
                     call getf(g_MultiWeight(I_VMAX))
                  case("CALCVAR")
                      if ( NWHTAY(1,I_VMAX).NE.-20 ) then
                          call report("Keyword "//trim(w)//"            &
     &                     only valid for HDIAG routine",.true.)
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
     &                 I_HMAX .ne. -19) then
                      call report(trim(w)//" only valid for MC "        &
     &                   //"method",.true.)
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
!                  if ( I_HMAX .ne. -7.and.
!     &                 I_HMAX .ne. -19) then
!                      call report(trim(w)//" only valid for MC " 
!     &                   //"method",.true.)
!                  end if
              case("MAXVERTICES")
                  if ( I_VMAX .ne. 0 ) then
                     call report("Cannot reset MAXVERTICES",.true.)
                  endif 
                  call readi(I_VMAX)
              case("IMPORTANCE")
                  call readf(G_VMC_PI)
!                  if ( I_HMAX .ne. -7 ) then
!                      call report(trim(w)//" only valid for MC "  
!     &                  //"method",.true.)
!                  end if
              case("SEED")
                  call readi(G_VMC_SEED)
!                  if ( I_HMAX .ne. -7 ) then
!                      call report(trim(w)//" only valid for MC " 
!     &                  //"method",.true.)
!                  end if
              case("BIAS")
                  call readf(G_VMC_FAC)

!                  if ( I_HMAX .ne. -7 ) then
!                      call report(trim(w)//" only valid for MC " 
!     &                //" method",.true.)
!                  end if
              case("STARCONVERGE")
                  call readf(STARCONV)
                if((NWHTAY(1,I_VMAX).ne.0).and.(NWHTAY(1,I_VMAX).ne.-21)&
     &                 .and.(NWHTAY(1,I_VMAX).ne.-9)) then
                      call report(trim(w)//" only valid for STAR "      &
     &                //" method",.true.)
                  end if
              case("UFORM-POWER")
                  TUPOWER=.true.
              case("CHEMPOTWEIGHTING")
                  call readf(g_VMC_ExcitWeights(1,1))
                  call readf(g_VMC_ExcitWeights(2,1))
                  call readf(G_VMC_EXCITWEIGHT(1))
                  DO l=1,5
                    IF(EXCITFUNCS(l)) THEN
                        call report(trim(w)//" only valid if "          &
     &               //" another weighting scheme not specified",.true.)
                    ENDIF
                  ENDDO
                  EXCITFUNCS(4)=.true.
              case("CHEMPOT-TWOFROM")
                  call readf(g_VMC_ExcitWeights(1,1))
                  call readf(g_VMC_ExcitWeights(2,1))
                  call readf(g_VMC_ExcitWeights(3,1))
                  call readf(G_VMC_EXCITWEIGHT(1))
                  DO l=1,5
                    IF(EXCITFUNCS(l)) THEN
                        call report(trim(w)//" only valid if "          &
     &               //" another weighting scheme not specified",.true.)
                    ENDIF
                  ENDDO
                  EXCITFUNCS(5)=.true.
              case("POLYEXCITWEIGHT")
                  call readf(g_VMC_ExcitWeights(1,1))
                  call readf(g_VMC_ExcitWeights(2,1))
                  call readf(g_VMC_ExcitWeights(3,1))
                  call readf(G_VMC_EXCITWEIGHT(1))
                  DO l=1,5
                    IF(EXCITFUNCS(l)) THEN
                        call report(trim(w)//" only valid if "          &
     &               //" another weighting scheme not specified",.true.)
                    ENDIF
                  ENDDO
                  EXCITFUNCS(2)=.true.
              case("POLYEXCITBOTH")
                  call readf(g_VMC_ExcitWeights(1,1))
                  call readf(g_VMC_ExcitWeights(2,1))
                  call readf(g_VMC_ExcitWeights(3,1))
                  call readf(g_VMC_ExcitWeights(4,1))
                  call readf(G_VMC_EXCITWEIGHT(1))
                  DO l=1,5
                    IF(EXCITFUNCS(l)) THEN
                        call report(trim(w)//" only valid if "          &
     &               //" another weighting scheme not specified",.true.)
                    ENDIF
                  ENDDO
                  EXCITFUNCS(3)=.true.
              case("EXCITWEIGHTING")
                  call readf(g_VMC_ExcitWeights(1,1))
                  call readf(g_VMC_ExcitWeights(2,1))
                  call readf(G_VMC_EXCITWEIGHT(1))
                  IF(item.lt.nitems) call readf(g_VMC_ExcitWeights(3,1))
                  DO l=1,5
                    IF(EXCITFUNCS(l)) THEN
                        call report(trim(w)//" only valid if "          &
     &               //" another weighting scheme not specified",.true.)
                    ENDIF
                  ENDDO
                  EXCITFUNCS(1)=.true.
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
     &              //" if delta_beta positive",.true.)
                     TNPDERIV = .false.
                 end if
              case("CIMC")
                  TMONTE = .true.
              case("MCSTEPS")
                  call readi(IMCSTEPS)
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space" &
     &                //" monte carlo is performed.",.true.)
                  end if
              case("EQSTEPS")
                  call readi(IEQSTEPS)
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space" &
     &                //" monte carlo is performed.",.true.)
                  end if
              case("BETAEQ")
                  call readf(BETAEQ)
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space" &
     &                //" monte carlo is performed.",.true.)
                  end if
              case("DETSYM")
                  TMCDET = .true.
                  do I = 1,5
                    call readi(MDK(I))
                  end do
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space" &
     &                 //" monte carlo is performed.",.true.)
                  end if
              case("DETINV")
                  call readi(DETINV)
              case("INSPECT")
                  TSPECDET = .true.
                  ALLOCATE(SPECDET(NEL-NFROZEN),STAT=ierr)
                  CALL MemAlloc(ierr,SPECDET,NEL-NFROZEN,'SPECDET')
!                  call MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')
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
                      call report("Warning - declared beta/p and p."    &
     &                //"Using p.",.true.)
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
!This is the number of vertices in the Graph Morph graph.
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
              case("NOCROSSING")
                  TNoCross=.true.
                  call report("NOCROSSING option not yet working",.true.)
              case("NOSAMEEXCIT")
                  TNoSameExcit=.true.
              case("NOTRIPLES")
                  lNoTriples=.true.
              case("MAXEXCIT")
!This imposes a maximum excitation level to the space that GraphMorph can explore
                  TMaxExcit=.true.
                  call geti(iMaxExcitLevel)
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
                  TMCExcitSpace=.true.
                  call geti(NoMCExcits)
              case default
                  call report("Keyword "                                &
     &              //trim(w)//" not recognized in CALC block",.true.)
              end select
            end do calc
          if(.not.tEnergy.and.I_VMAX.eq.1)  tNeedsVirts=.false.! Set if we need virtual orbitals  (usually set).  Will be unset (by Calc readinput) if I_VMAX=1 and TENERGY is false
        END SUBROUTINE CalcReadInput
        Subroutine CalcInit()
            Use System, only: G1, Alat, Beta, BRR, ECore, LMS, nBasis, nBasisMax, STot,tCSF
            Use System, only: tUEG
            Use Integrals, only: FCK, CST, nMax, nMsh, UMat
            Use Integrals, only: HFEDelta, HFMix, NHFIt, tHFCalc
            Use Determinants, only: FDet, tSpecDet, SpecDet, GetHElement2
            Use DetCalc, only: DetInv, nDet, tRead
            
            REAL*8 CalcT, CalcT2, GetRhoEps
            
            
            INTEGER I, IC
            INTEGER nList

        Allocate(MCDet(nEl))

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
         CALL ICOPY(NEL,SPECDET,1,FDET,1)
         CALL GETCSFFROMDET(FDET,SPECDET,NEL,STOT,LMS)
         WRITE(6,*) "CSF with 2S=",STOT," and 2Sz=",LMS," now in SPECDET is"
         CALL WRITEDET(6,SPECDET,NEL,.TRUE.)
      ENDIF
      IF(TSPECDET.AND.SPECDET(1).EQ.0) THEN
         WRITE(6,*) "TSPECDET set, but invalid.  using FDET"
         CALL ICOPY(NEL,FDET,1,SPECDET,1)
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

!      IF(G_VMC_FAC.LE.0) THEN
!         WRITE(6,*) "G_VMC_FAC=",G_VMC_FAC
!         STOP "G_VNC_FAC LE 0"
!      ENDIF

       IF(BETAP.NE.0) THEN 
         I_P=NINT(BETA/BETAP)
         WRITE(6,*) ' BETAP=',BETAP
         WRITE(6,*) ' RESETTING P '
         IF(I_P.GT.100000) WRITE(6,*) ' *** WARNING I_P=',I_P
      ENDIF

      WRITE(6,*) ' BETA, P :',BETA,I_P
       
!C      DBRAT=0.001
!C      DBETA=DBRAT*BETA
      WRITE(6,*) "DBETA=",DBETA

      IF(.NOT.TREAD) THEN
!         CALL WRITETMAT(NBASIS)
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
!C         CALL GENRANDOMDET(NEL,NBASIS,MCDET)
         DO I=1,NEL
            MCDET(I)=FDET(I)
         ENDDO
      ENDIF
      IF(TMONTE) THEN
         WRITE(6,"(A,$)") 'MC Start Det: '
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
        Use System, only: Alat, Arr,Brr, Beta, ECore, G1, LMS, LMS2, nBasis, nBasisMax
        Use System, only: SymRestrict, tCSF, tParity, tSpn, ALat, Beta
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB
        Use DetCalc, only : CK, DetInv, nDet, nEval, tEnergy, tRead, nmrks, w
        Use Determinants, only: FDet, nActiveBasis, SpecDet, tSpecDet
        Use Integrals, only: FCK, NMAX, NMSH, UMat, FCK
        Use Integrals, only: HFEDelta, HFMix,nTay
        Use Logging, only: iLogging
!Calls
        REAL*8 DMonteCarlo2
!Local Vars
        REAL*8 EN, ExEn, GsEN
        REAL*8 FLRI, FLSI
        REAL*8 RH
        LOGICAL tWarn
    
    !C.. we need to calculate a value for RHOEPS, so we approximate that
!C.. RHO_II~=exp(-BETA*H_II/p).  RHOEPS is a %ge of this 
!C.. If we haven't already calced RHOEPS, do it now
        Call DoExactVertexCalc()
        IF(NPATHS.NE.0.OR.DETINV.GT.0) THEN
!Old and obsiolecte
!            IF(TRHOIJND) THEN
!C.. We're calculating the RHOs for interest's sake, and writing them,
!C.. but not keeping them in memory
               !WRITE(6,*) "Calculating RHOS..."
!               WRITE(6,*) "Using approx NTAY=",NTAY
!               CALL CALCRHOSD(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,        &
!     &            NBASISMAX,G1,nBasis,BRR,NMSH,FCK,NMAX,ALAT,UMAT,             &
!     &            NTAY,RHOEPS,NWHTAY,ECORE)
!            ENDIF
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
                    CALL CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,           &
     &           NBASISMAX,G1,nBasis,BRR,NMSH,FCK,ARR,ALAT,UMAT,NTAY,              &
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,           &
     &            DETINV,TSPECDET,SPECDET)
               ELSE
                  IF(TCSF) THEN
                     IF(.NOT.TSPECDET) THEN
                        WRITE(6,*) "SPECDET not specified. Using Fermi determinant ONLY"
                        TSPECDET=.TRUE.
                        CALL ICOPY(NEL,FDET,1,SPECDET,1)
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
                CALL CALCRHOPII3(BETA,I_P,I_HMAX,I_VMAX,NEL,                          &
     &           NBASISMAX,G1,nBasis,BRR,NMSH,FCK,ARR,ALAT,UMAT,NTAY,                  &
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,               &
     &            DETINV,TSPN,LMS2,TPARITY,SymRestrict,TSPECDET,SPECDET,            &
     &            nActiveBasis)
            ENDIF
        ELSE
            WRITE(6,*) "Invalid combination of NTAY and TENERGY.  No NPATHS calculated"
            WRITE(6,*) "NTAY: ",NTAY(1)," TENERGY: ",TENERGY
        ENDIF
      ENDIF
      IF(TMONTE) THEN
!C            DBRAT=0.01
!C            DBETA=DBRAT*BETA
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
!C..      FCK=W
!C..      ZIA=CK
!C..      UMAT=NDET
!C..      ALAT=NMRKS
!C..      NMAX=ARR
                  EN=DMONTECARLO2(MCDET,I_P,BETA,DBETA,I_HMAX,I_VMAX,IMCSTEPS,          &
     &               G1,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS,                                 &
     &         NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,RHOEPS,NWHTAY,ILOGGING,ECORE,BETAEQ)
               ELSE
                  
                  STOP "TENERGY not set, but NTAY=0" 
               ENDIF
            ENDIF
            WRITE(6,*) "MC Energy:",EN
!CC            WRITE(12,*) DBRAT,EN
       ENDIF
         
!C.. /AJWT
        End Subroutine
    Subroutine DoExactVertexCalc()
        Use System, only: Alat, Beta, Brr, ECORE, G1, nBasis, nBasisMax, Arr
        use System, only: Symmetry,SymmetrySize,SymmetrySizeB,BasisFN,BasisFNSize,BasisFNSizeB
        Use Integrals, only: fck, nMax, nMsh, UMat,nTay
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
!C               WRITE(6,*) RH
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
!C         CLOSE(17)
         CLOSE(14)
      ENDIF
    
      IF(TMONTE.AND.TENERGY.AND.NTAY(1).EQ.-1) THEN
         WRITE(6,*) "Calculating Exact MC Energy..."
         EN=DOEXMC(NDET,NEVAL,CK,W,BETA,I_P,ILOGGING,0.D0,IMCSTEPS,G1,NMRKS,NEL,NBASISMAX,nBasis,BRR,IEQSTEPS)
      ENDIF
      IF(TBEGRAPH) THEN
         IF(TENERGY) THEN
            IF(NTAY(1).NE.0) THEN
               CALL DOBEGRAPH(NDET,NEVAL,CK,W,I_P,ILOGGING,G1,NMRKS,nEl,NBASISMAX,nBasis,BRR)
            ELSE
!C.. NTAY=0 signifying we're going to calculate the RHO values when we
!C.. need them from the list of eigenvalues.   
!C.. Hide NMSH=NEVAL
!C..      FCK=W
!C..      ZIA=CK
!C..      UMAT=NDET
!C..      ALAT=NMRKS        
               CALL DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,             &
     &            NBASISMAX,G1,nBasis,BRR,NEVAL,W,CK,NMAX,NMRKS,NDET,      &
     &            NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
            ENDIF
         ELSE
            CALL DOBEGRAPHAP(I_P,I_HMAX,I_VMAX,NEL,NDET,                &
     &            NBASISMAX,G1,nBasis,BRR,NMSH,FCK,NMAX,ALAT,UMAT,         &
     &            NTAY,RHOEPS,NWHTAY,NPATHS,ILOGGING)
         ENDIF
      ENDIF
            IF(NTAY(1).EQ.0.AND.TENERGY) THEN
               WRITE(6,*) "Using exact RHOs generated on the fly"
!C.. we've calculated energies, and we're passing them through to
!C.. calculate the exact RHOS
!C.. NTAY=0 signifying we're going to calculate the RHO values when we
!C.. need them from the list of eigenvalues.   
!C.. Hide NMSH=NEVAL
!C..      FCK=W
!C..      ZIA=CK
!C..      UMAT=NDET
!C..      ALAT=NMRKS
!C..      NMAX=ARR
              CALL CALCRHOPII2(NMRKS,BETA,I_P,I_HMAX,I_VMAX,NEL,NDET,           &
     &           NBASISMAX,G1,nBasis,BRR,NEVAL,W,CK,ARR,NMRKS,NDET,NTAY,           &
     &            RHOEPS,NWHTAY,NPATHS,ILOGGING,ECORE,TNPDERIV,DBETA,           &
     &            DETINV,TSPECDET,SPECDET)

            endif
            
                End Subroutine DoExactVertexCalc
         Subroutine CalcCleanup()
         End Subroutine CalcCleanup
      END MODULE Calc
      subroutine inpgetmethod(I_HMAX,NWHTAY,I_V)
         use input
         use UMatCache , only : TSTARSTORE
         USE Calc , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TMCDIRECTSUM,g_Multiweight,G_VMC_FAC,TMPTHEORY
         USE Calc, only : STARPROD,TDIAGNODES,TSTARSTARS,TGraphMorph,TStarTrips
         implicit none
         integer I_HMAX,NWHTAY,I_V
         CHARACTER(LEN=16) w
                  do while ( item .lt. nitems )
                    call readu(w)
                    select case(w)
                    case("VERTEX")
                        call readu(w)
                        select case(w)
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
                              case("STARPROD")
                                 STARPROD=.TRUE.
                              case("TRIPLES")
                                  TStarTrips=.TRUE.
                              case("COUNTEXCITS")
                                 NWHTAY=IBSET(NWHTAY,8)
                              case("ADDSINGLES")
                                 NWHTAY=IBSET(NWHTAY,7)
                                 IF(I_HMAX.NE.-21)  call report(        &
     &                              "Error - cannot use ADDSINGLES"     &
     &                              //" without STAR NEW",.true.)
                                 IF(TSTARSTORE) call report("Error - "  &
     &                            //"can only use STARSTOREREAD with "  &
     &                            //"double excitations of HF",.true.)
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
     &                       //"can only be specified with POLY... NEW")
                              case default
                                call report("Error - must specify DIAG" &
     &                        //" or POLY vertex star method",.true.)
                               end select
                           enddo
!                           IF(TSTARSTARS.and..not.BTEST(NWHTAY,0)) THEN 
!                               call report("STARSTARS must be used with " &
!     &                          //"a poly option",.true.)
!                           ENDIF
                           IF(STARPROD.and.BTEST(NWHTAY,0)) THEN
                               call report("STARPROD can only be "      &
     &                        //"specified with DIAG option",.true.)
                            ENDIF
                           if(i_hmax.eq.0)                              &
     &                   call report("OLD/NEW not specified for STAR",  &
     &                          .true.)
                        case default
                        call report("Keyword error with "//trim(w),     &
     &                          .true.)
                        end select
                    case default
                        call report("Error.  Method not specified."     &
     &                    //" Stopping.",.true.)
                    end select
               end do
      end
      subroutine inpgetexcitations(NWHTAY,w)
         use input
         IMPLICIT NONE
         INTEGER NWHTAY
         CHARACTER(LEN=16) w
!                  call readu(w)
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
                        call report("Keyword error with EXCITATIONS "   &
     &                     //trim(w),                                   &
     &                          .true.)
                  end select
      end

