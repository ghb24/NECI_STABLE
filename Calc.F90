
      MODULE Calc
        USE input
        USE System , only : NEL,Feb08,defaults
        USE Integrals , only : NFROZEN
        IMPLICIT NONE

        LOGICAL TCALCHMAT,TSTAR,TENERGY,TREAD,TBLOCK,TTROT
        LOGICAL TNEWEXCITATIONS,TVARCALC(0:10),TBIN,TVVDISALLOW
        LOGICAL TMCDIRECTSUM,TMPTHEORY,TMODMPTHEORY,TUPOWER
        LOGICAL EXCITFUNCS(10),TNPDERIV,TMONTE,TMCDET,TSPECDET
        LOGICAL TBETAP,CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,TENPT
        LOGICAL TLADDER,TRHOOFR,TCORR,TFODM,TMC,TREADRHO,TRHOIJ
        LOGICAL TBEGRAPH,STARPROD,TDIAGNODES,TSTARSTARS
        
        INTEGER NEVAL,NBLK,NKRY,ICILEVEL,NWHTAY(3,10),NCYCLE,NPATHS
        INTEGER IOBS,JOBS,KOBS,NDETWORK,I_HMAX,I_VMAX,G_VMC_SEED
        INTEGER nActiveSpace(2),IMCSTEPS,IEQSTEPS,MDK(5),DETINV
        INTEGER CUR_VERT,NHISTBOXES,I_P,LinePoints
        
        REAL*8 B2L,g_MultiWeight(0:10),G_VMC_PI,G_VMC_FAC,BETAEQ
        REAL*8 G_VMC_EXCITWEIGHT(10),G_VMC_EXCITWEIGHTS(6,10),BETA
        REAL*8 BETAP,RHOEPSILON,DBETA(3),STARCONV,CHEMPOT

        INTEGER, DIMENSION(:), POINTER :: SPECDET

        contains

        SUBROUTINE CalcReadInput()
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
                     call reread(-1)
                     call geti(NPATHS)
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
!This indicates the number of times the eigenvalues of the star matrix should be evaluated to achieve the linear approximation when STARSTARS set,
              case("LINEPOINTSSTAR")
                  call geti(LinePoints)
              case default
                  call report("Keyword "                                &
     &              //trim(w)//" not recognized in CALC block",.true.)
              end select
            end do calc
        END SUBROUTINE CalcReadInput
      END MODULE Calc
      subroutine inpgetmethod(I_HMAX,NWHTAY,I_V)
         use input
         use UMatCache , only : TSTARSTORE
         USE Calc , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT,         &
     &          TMCDIRECTSUM,g_Multiweight,G_VMC_FAC,TMPTHEORY,         &
     &          STARPROD,TDIAGNODES,TSTARSTARS
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
                                call report("Error - must specify OLD"  &
     &                         //" or NEW vertex sum method",.true.)
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

