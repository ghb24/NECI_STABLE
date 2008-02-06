
      MODULE CALCREAD
        USE input
        USE SYSREAD , only : NEL,Feb08,defaults
        USE INTREAD , only : NFROZEN
        IMPLICIT NONE

        LOGICAL TCALCHMAT,TSTAR,TENERGY,TREAD,TBLOCK,TTROT
        LOGICAL TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
        LOGICAL TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER
        LOGICAL EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET,TSPECDET
        LOGICAL TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
        LOGICAL TLADDER,TRHOOFR,TCORR,TFODM,TMC,TREADRHO,TRHOIJ
        LOGICAL TBEGRAPH,STARPROD,TDIAGNODES
        
        INTEGER NEVAL,NBLK,NKRY,ICILEVEL,NWHTAY(3,10),NCYCLE,NPATHS
        INTEGER IOBS,JOBS,KOBS,NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED
        INTEGER nActiveSpace(2),IMCSTEPS,IEQSTEPS,MDK(5),DETINV
        INTEGER CUR_VERT,NHISTBOXES,I_P
        
        REAL*8 B2L,g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
        REAL*8 G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10),BETA
        REAL*8 BETAP,RHOEPSILON,DBETA(3),STARCONV,CHEMPOT

        INTEGER, DIMENSION(:), POINTER :: SPECDET

        contains

        SUBROUTINE readinputcalc()
        IMPLICIT NONE
        LOGICAL eof
        CHARACTER (LEN=100) w
        INTEGER :: l,i,ierr

!     Values for old parameters.
!     These have no input options to change the defaults, but are used in the code.
      TRHOOFR = .false.
      TCORR = .false.
      TFODM = .false.
      TMC = .false.
      NHISTBOXES = 0
      TREADRHO = .false.
      TRHOIJ = .false.
      TBEGRAPH = .false.


!       Calc defaults      
          TDIAGNODES=.false.
          STARPROD=.false.
          TCALCHMAT = .false.
          TStar=.false.
          TENERGY = .false.
          NEVAL = 0
          TREAD = .false.
          NBLK = 4
          NKRY = 8
          B2L = 1.D-13
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
                     read(w,"(I)") NPATHS
                  end select
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
              case default
                  call report("Keyword "                                &
     &              //trim(w)//" not recognized in CALC block",.true.)
              end select
            end do calc
        END SUBROUTINE

      END MODULE CALCREAD
