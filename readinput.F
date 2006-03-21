      SUBROUTINE READINPUT(NEL,TSPN,LMS,TPARITY,IPARITY,NMAXX,NMAXY,
     &   NMAXZ,
     &	NMSH,BOX,BOA,COA,TALPHA,ALPHA,ISTATE,NEVAL,B2L,NCYCLE,NBLK,
     &	NKRY,TRHOOFR,TREAD,TCORR,TFODM,IOBS,JOBS,KOBS,TENERGY,TMC,
     &	NHISTBOXES,BETA,I_P,TRHOIJ,RHOEPSILON,TREADRHO,NPATHS,TMONTE,
     &  IMCSTEPS,IEQSTEPS,TTROT,I_HMAX,TCALCHMAT,NTAY,NWHTAY,I_VMAX,
     &  ILOGGING,TUEG,TBEGRAPH,THFBASIS,THUB,UHUB,BHUB,ETRIAL,TNPDERIV,
     & BETAP,BETAEQ,TMCDET,MDK,DETINV,THFCALC,NHFIT,TTILT,ITILTX,ITILTY,
     & TBLOCK,HFMIX,HFEDELTA,TREADINT,TPBC,TREAL,NDETWORK,ORBORDER,
     &   NFROZEN,NTFROZEN,TRHF,HFCDELTA,IHFMETHOD,TREADHF,TSPECDET,
     &  IP_SPECDET,DBETA,TREADTUMAT,THFORDER,TCSF,STOT,ICILEVEL,FUEGRS)
      USE input
      IMPLICIT REAL*8    (A-H,O-Z)
      DIMENSION IPARITY(3)
      LOGICAL TSPN,TALPHA,TPARITY,TCORR,TFODM,TCALCHMAT,TUEG,     
     &   TREAD,TRHOOFR,TENERGY,TMC,TRHOIJ,TREADRHO,TMONTE,TTROT,  
     &   TBEGRAPH,THFBASIS,THUB,TNPDERIV,TMCDET,THFCALC,TTILT,    
     &   TBLOCK,TREADINT,TPBC,TREAL,TRHF,TREADHF,TSPECDET,TREADTUMAT,  
     &   THFORDER,TCSF,TBETAP, eof
      CHARACTER*255 INP
      INTEGER SPECDET(*)
!
      CHARACTER(LEN=16) title  
      CHARACTER(LEN=16) w
      INTEGER NPATHS,STOT
      INTEGER MDK(5),DETINV,ORBORDER(8,2)
      INTEGER NWHTAY(2,10)
      REAL*8  ORBORDER2(8)
      POINTER (IP_SPECDET,SPECDET)
      INCLUDE 'vmc.inc'
      INCLUDE 'fcoul.inc'
!     --------------------------------------------
!     Set some sensible defaults here
      title = ""
      ISTATE = 1
      NMAXX = 0
      NMAXY = 0
      NMAXZ = 0
      TUEG = .false.
      THUB = .false.
      UEG = .false.
      UHUB = 4
      BHUB = -1
      TPBC = .false.
      TREAL = .false.
      TTILT = .false.
      TALPHA = .false.
      NEVAL = 0
      TBLOCK = .false.
      ICILEVEL = 0
      B2L = 1.D-13
      NCYCLE = 200
      NBLK = 4
      NKRY = 8
      TREAD = .false.
      TCALCHMAT = .false.
      TENERGY = .false.
      NWHTAY(1,1) = 0
      G_VMC_PI = 0.95
      G_VMC_SEED = -7
      G_VMC_FAC = 16

      LMS=0
      STOT=0
      TSPN=.false.
      TCSF=.false.

      BETAP=1.D-4
      DBETA=-1.D0
      RHOEPSILON=1.D-6
      TTROT=.true.

      TNPDERIV = .false.
      TMONTE = .false.
      IMCSTEPS = 0
      IEQSTEPS = 0
      BETAEQ = 0
      IHFMETHOD = 1
      TREADBASIS = .false.
      TREADUMAT = .false.
      TMCDET = .false.
      FCOUL=1.D0
      FUEGRS=0.D0
      do I = 1,5
        MDK(I) = 0
      end do
      DETINV = 0
      NDETWORK = 50000
      TSPECDET = .false.
!     Values for old parameters.
      TRHOOFR = .false.
      TCORR = .false.
      TFODM = .false.
      TMC = .false.
      NHISTBOXES = 0
      NTAY = 1
      TREADRHO = .false.
      TRHOIJ = .false.
      TBEGRAPH = .false.

      IOBS=16
      JOBS=16
      KOBS=16
      NMSH=32

      ir=1
      IF(IARGC().GT.0) THEN
         CALL GETARG(1,INP)
         OPEN(1,FILE=INP,STATUS="OLD")
      ELSE
         ir=5
      ENDIF

!     --------------------------------------------
      call input_options(echo_lines=.false.,skip_blank_lines=.true.)
      main: do
        call read_line(eof)
        if (eof) exit          ! Exit if at end of input file.
        call readu(w)
        select case(w)
        case("TITLE")
            do while ( item .lt. nitems )
              call reada(w)
              title = trim(title)//" "//trim(w)
            end do 
            write (6,*) title
        case("SYSTEM")
            call readu(w)
            select case(w)
            case("READ","GENERIC")
                TREADINT = .true.
                call readu(w)
                select case(w)
                case("ORDER")
                    THFORDER = .true.
                end select
            case("HUBBARD")
                THUB = .true.
                TPBC=.true.
            case("UEG")
                TUEG = .true.
            case default
                call report ("System type "//trim(w)    
     &                       //" not valid",.true.)
            end select
            system: do
              call read_line(eof)
              if (eof) then
                  call report("Incomplete input file",.true.)
              end if
              call readu(w)
              select case(w)
              case("ELECTRONS","NEL")
                  call geti(NEL)
              case("SPIN-RESTRICT")
                  if(item.lt.nitems) then
                     call geti(LMS)
                  else
                     LMS=0
                  endif
                  TSPN = .true.
              case("CSF")
                  if(item.lt.nitems) then
                     call geti(STOT)
                  else
                     STOT=0
                  endif
                  TCSF = .true.
              case("SYM")
                  TPARITY = .true.
                  do I = 1,3
                    call geti(IPARITY(I))
                  end do
              case("CELL")
                  call geti(NMAXX)
                  call geti(NMAXY)
                  call geti(NMAXZ)
              case("MESH")
                  call geti(NMSH)
              case("BOXSIZE")
                  call getf(BOX)
                  if(item.lt.nitems) then
                     call getf(BOA)
                     call getf(COA)
                  else
                     BOA=1.D0
                     COA=1.D0
                  endif
              case("RS")
                  call getf(FUEGRS)
              case("COULOMB")
                  call getf(FCOUL)
              case("U")
                  call getf(UHUB)
              case("B")
                  call getf(BHUB)
              case("REAL")
                  TREAL = .true.
              case("APERIODIC")
                  TPBC = .false.
              case("TILT")
                  TTILT = .true.
                  call geti(ITILTX)
                  call geti(ITILTY)
              case("ALPHA")
                  TALPHA = .true.
                  call getf(ALPHA)
              case("STATE")
                  call geti(ISTATE)
                  if ( ISTATE /= 1 ) then
                      call report("Require ISTATE to be left set as"  
     &                //" 1",.true.)
                  end if
              case("ENDSYS") 
                  exit system
              case default
                  call report("Keyword "         
     &              //trim(w)//" not recognized in SYSTEM block",.true.)
              end select
            end do system
            if(NEL.eq.0)
     &         call report("Number of electrons cannot be zero.",.true.)
            if(THUB.OR.TUEG.OR..NOT.TREADINT) then
               if(NMAXX.EQ.0)
     &            call report("Must specify CELL "
     &            //"- the number of basis functions in each dim.",
     &            .true.)
               if(.NOT.THUB.AND.BOX.EQ.0.D0)
     &            call report("Must specify BOX size.",.true.)
               if(TTILT.AND..NOT.THUB)
     &            call report("TILT can only be specified with "
     &         //"HUBBARD.",.true.)
            endif
        case("CALC")
            calc: do
              call read_line(eof)
              if (eof) then
                  exit
              end if
              call readu(w)
              select case(w)
              case("HAMILTONIAN")
                  TCALCHMAT = .true.
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
              case("METHOD")
                  do while ( item .lt. nitems )
                    call readu(w)
                    select case(w)
                    case("VERTEX")
                        call readu(w)
                        select case(w)
                        case("SUM")
                            call readu(w)
                            select case(w)
                            case("OLD")
                                I_HMAX = -1
                            case("NEW")
                                I_HMAX = -8
                            case default
                                call report("Error - must specify OLD" 
     &                         //" or NEW vertex sum method",.true.)
                            end select
                        case default
                        call report("Keyword error with "//trim(w),  
     &                          .true.)
                        end select
                    case("MC")
                        I_HMAX = -7
                    case default
                        call report("Error.  Method not specified."    
     &                    //" Stopping.",.true.)
                    end select
                  end do
              case("CYCLES")
                  call readi(NWHTAY(1,1))
                  if ( I_HMAX .ne. -7 ) then
                      call report(trim(w)//" only valid for MC " 
     &                   //"method",.true.)
                  end if
              case("MAXVERTICES")
                  call readi(I_VMAX)
              case("IMPORTANCE")
                  call readf(G_VMC_PI)
                  if ( I_HMAX .ne. -7 ) then
                      call report(trim(w)//" only valid for MC "  
     &                  //"method",.true.)
                  end if
              case("SEED")
                  call readi(G_VMC_SEED)
                  if ( I_HMAX .ne. -7 ) then
                      call report(trim(w)//" only valid for MC " 
     &                  //"method",.true.)
                  end if
              case("BIAS")
                  call readf(G_VMC_FAC)
                  if ( I_HMAX .ne. -7 ) then
                      call report(trim(w)//" only valid for MC " 
     &                //" method",.true.)
                  end if
              case("PATHS")
                  call readi(NPATHS)
              case("ALLPATHS")
                  NPATHS = -1
              case("DERIV")
                  TNPDERIV = .true.
                 if (DBETA .lt. 0 ) then
                    call report("Only calculate energy with derivatives"
     &              //" if delta_beta positive",.true.)
                     TNPDERIV = .false.
                 end if
              case("CIMC")
                  TMONTE = .true.
              case("MCSTEPS")
                  call readi(IMCSTEPS)
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space"
     &                //" monte carlo is performed.",.true.)
                  end if
              case("EQSTEPS")
                  call readi(IEQSTEPS)
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space"
     &                //" monte carlo is performed.",.true.)
                  end if
              case("BETAEQ")
                  call readf(BETAEQ)
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space"
     &                //" monte carlo is performed.",.true.)
                  end if
              case("DETSYM")
                  TMCDET = .true.
                  do I = 1,5
                    call readi(MDK(I))
                  end do
                  if ( .not. TMONTE ) then
                      call report(trim(w)//" only relevant if CI space"
     &                 //" monte carlo is performed.",.true.)
                  end if
              case("DETINV")
                  call readi(DETINV)
              case("INSPECT")
                  TSPECDET = .true.
                  call MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')
                  SPECDET(1)=0
                  if(item.lt.nitems) then
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
                      call report("Warning - declared beta/p and p."  
     &                //"Using p.",.true.)
                  end if
              case("DELTABETA")
                  call getf(DBETA)
              case("RHOEPSILON")
                  call getf(RHOEPSILON)
              case default
                  call report("Keyword "         
     &              //trim(w)//" not recognized in CALC block",.true.)
              end select
            end do calc
        case("INTEGRAL")
            integral: do
              call read_line(eof)
              if (eof) then
                  exit
              end if
              call readu(w)
              select case(w)
              case("HF")
                  THFBASIS = .true.
              case("CALCULATE")
                  THFCALC = .true.
              case("MAXITERATIONS")
                  call geti(NHFIT)
              case("MIX")
                  call getf(HFMIX)
              case("THRESHOLD")
                  do while ( item .lt. nitems )
                    call readu(w)
                    select case(w)
                    case("ENERGY")
                        call readf(HFEDELTA)
                    case("ORBITAL")
                        call readf(HFCDELTA)
                    case default
                        call report(trim(w)//" not valid THRESHOLD" 
     &               //"OPTION.  Specify ENERGY or ORBITAL convergence"
     &                 //" threshold.",.true.)
                    end select
                  end do
              case("RHF")
                  TRHF = .true.
              case("UHF")
                  TRHF = .false.
              case("HFMETHOD")
                  call readu(w)
                  select case(w)
                  case("DESCENT")
                      call readu(w)
                      select case(w)
                      case("OTHER")
                          IHFMETHOD = 2
                      case("SINGLES")
                          IHFMETHOD = 1
                      case default
                          call report(trim(w)//" not valid DESCENT"
     &                    //" option",.true.)
                      end select
                  case("STANDARD")
                      IHFMETHOD = 0
                  case default
                      call report(trim(w)//" not valid HF method", 
     &                 .true.)
                  end select
              case("READ")
                  do while ( item .lt. nitems )
                    call readu(w)
                    select case(w)
                    case("MATRIX")
                        TREADTUMAT = .true.
                    case("BASIS")
                        TREADHF = .true.
                    case default
                        call report(trim(w)//" is an invalid HF read"
     &                  //" option.",.true.)
                    end select
                  end do
              case("FREEZE")
                  call readi(NFROZEN)
                  call readi(NTFROZEN)
                  if ( NPATHS .ne. 1 ) then
                      call report("Can only freeze orbitals when " 
     &                 //"NPATHS=1",.true.)
                  end if
                  if ( mod(NFROZEN,2).ne.0 .or. mod(NTFROZEN,2).ne.0 ) 
     &             then
                      call report("NFROZEN and NTFROZEN must be"  
     &                //"multiples of 2",.true.)
                  end if
              case("ORDER")
                  I = 1
                  do while ( item .lt. nitems )
                    call readf(ORBORDER2(I))
                    I = I + 1
                  end do
                  DO I=1,8
C.. two ways of specifying open orbitals
C.. if orborder2(I,1) is integral, then if it's odd, we have a single
C.. open orbital
                     IF(ORBORDER2(I).EQ.INT(ORBORDER2(I))) THEN
                        ORBORDER(I,1)=IAND(INT(ORBORDER2(I)),65534)
                        IF((INT(ORBORDER2(I))-ORBORDER(I,1)).GT.0) THEN
C.. we have an open orbital
                           ORBORDER(I,2)=2
                        ELSE
                           ORBORDER(I,2)=0
                        ENDIF
                     ELSE
C.. non-integral.  The integral part is the number of closed oribtals,
C.. and the fractional*1000 is the number of open orbitals.
C.. e.g. 6.002 would mean 6 closed and 2 open
C.. which would have orborder(I,1)=6, orborder(I,2)=4
C.. but say 5.002 would be meaningless as the integral part must be a
C.. multiple of 2
                        ORBORDER(I,1)=INT(ORBORDER2(I)+0.000001)
                        ORBORDER(I,2)=INT((ORBORDER2(I)-ORBORDER(I,1)+
     &                                 0.000001)*1000)*2
                     ENDIF
                  ENDDO
              case("ENDINT")
                   exit integral
              case default
                  call report("Keyword "  
     &          //trim(w)//" not recognized in INTEGRAL block",.true.)
              end select
            end do integral
        case("LOGGING")
            ILOGGING = 0
            logging: do
              call read_line(eof)
              if (eof) then
                  exit
              end if
              call readu(w)
              select case(w)
              case("MCPATHS")
                  ILOGGING = ILOGGING + 2**0
              case("FMCPR")
                  call readu(w)
                  select case(w)
                  case("LABEL")
                      ILOGGING = ILOGGING + 2**1
                  case("RHO")
                      ILOGGING = ILOGGING + 2**2
                  case default
                      ILOGGING = ILOGGING + 2**1
                  end select
              case("CALPATH")
                  call readu(w)
                  select case(w)
                  case("LABEL")
                      ILOGGING = ILOGGING + 2**3
                  case("RHO")
                      ILOGGING = ILOGGING + 2**4
                  case default
                      ILOGGING = ILOGGING + 2**3
                  end select
              case("XIJ")
                  ILOGGING = ILOGGING + 2**5
              case("HAMILTONIAN")
                  ILOGGING = ILOGGING + 2**6
              case("PSI")
                  ILOGGING = ILOGGING + 2**7
              case("VERTEX")
                  call readu(w)
                  select case(w)
                  case("1000")
                      ILOGGING = ILOGGING + 2**8
                  case("EVERY")
                      ILOGGING = ILOGGING + 2**9
                  end select
              case("HFBASIS")
                  ILOGGING = ILOGGING + 2**10
              case("ENDLOG")
                  exit logging
              end select
            end do logging
        case("END")
            exit
!       Error message if unknown keyword
        case default
            call report ("Keyword "//trim(w)//" not recognized",.true.)
        end select
      end do main
C..   Not sure why this statement is necessary, even if TSPECDET =
C..   false.
      CALL MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')
C..   Testing ILOGGING
C     ILOGGING = 0771
      WRITE (6,*) 'ILOGGING after input routine', ILOGGING
      IF(IR.EQ.1) CLOSE(1)
      END
C ==------------------------------------------------------------------==

      SUBROUTINE READOLDINPUT(NEL,TSPN,LMS,TPARITY,IPARITY,NMAXX,NMAXY,
     &   NMAXZ,
     &  NMSH,BOX,BOA,COA,TALPHA,ALPHA,ISTATE,NEVAL,B2L,NCYCLE,NBLK,
     &  NKRY,TRHOOFR,TREAD,TCORR,TFODM,IOBS,JOBS,KOBS,TENERGY,TMC,
     &  NHISTBOXES,BETA,I_P,TRHOIJ,RHOEPSILON,TREADRHO,NPATHS,TMONTE,
     &  IMCSTEPS,IEQSTEPS,TTROT,I_HMAX,TCALCHMAT,NTAY,NWHTAY,I_VMAX,
     &  ILOGGING,TUEG,TBEGRAPH,THFBASIS,THUB,UHUB,BHUB,ETRIAL,TNPDERIV,
     & BETAP,BETAEQ,TMCDET,MDK,DETINV,THFCALC,NHFIT,TTILT,ITILTX,ITILTY,
     & TBLOCK,HFMIX,HFEDELTA,TREADINT,TPBC,TREAL,NDETWORK,ORBORDER,
     &   NFROZEN,NTFROZEN,TRHF,HFCDELTA,IHFMETHOD,TREADHF,TSPECDET,
     &  IP_SPECDET,DBETA,TREADTUMAT,THFORDER,TCSF,STOT,ICILEVEL,FUEGRS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPARITY(3)
      LOGICAL TSPN,TALPHA,TPARITY,TCORR,TFODM,TCALCHMAT,TUEG,
     &   TREAD,TRHOOFR,TENERGY,TMC,TRHOIJ,TREADRHO,TMONTE,TTROT,
     &   TBEGRAPH,THFBASIS,THUB,TNPDERIV,TMCDET,THFCALC,TTILT,
     &   TBLOCK,TREADINT,TPBC,TREAL,TRHF,TREADHF,TSPECDET,TREADTUMAT,
     &   THFORDER,TCSF,TSPECDET2
      CHARACTER*255 INP
      INTEGER NPATHS,STOT
      INTEGER MDK(5),DETINV,ORBORDER(8,2)
      INTEGER SPECDET(*)
      INTEGER NWHTAY(2,10)
      REAL*8 ORBORDER2(8)
      POINTER (IP_SPECDET,SPECDET)
      REAL*8 FUEGRS
      INCLUDE 'vmc.inc'
      INCLUDE 'fcoul.inc'

      IFILE=1
      IF(IARGC().GT.0) THEN
         CALL GETARG(1,INP)
         OPEN(1,FILE=INP,STATUS="OLD")
      ELSE
         IFILE=5
      ENDIF
C..Read in input parameters
      READ(IFILE,*) NEL
      READ(IFILE,*) TSPN,TCSF
      READ(IFILE,*) LMS,STOT
      READ(IFILE,*) TREADINT,THFORDER
      READ(IFILE,*) TPARITY,(IPARITY(I),I=1,3)
      READ(IFILE,*) NMAXX,NMAXY,NMAXZ
      READ(IFILE,*) NMSH
      READ(IFILE,*) BOX,BOA,COA
      READ(IFILE,*) TUEG,FUEGRS,FCOUL
      READ(IFILE,*) THUB,UHUB,BHUB,TPBC,TREAL
      READ(IFILE,*) TTILT,ITILTX,ITILTY
      READ(IFILE,*) TALPHA,ALPHA
      READ(IFILE,*) ISTATE
      READ(IFILE,*) NEVAL
      READ(IFILE,*) TBLOCK,ICILEVEL
      READ(IFILE,*) B2L
      READ(IFILE,*) NCYCLE
      READ(IFILE,*) NBLK
      READ(IFILE,*) NKRY
      READ(IFILE,*) TRHOOFR
      READ(IFILE,*) TREAD
      READ(IFILE,*) TCORR
      READ(IFILE,*) TFODM
      READ(IFILE,*) IOBS,JOBS,KOBS 
      READ(IFILE,*) TCALCHMAT
      READ(IFILE,*) TENERGY
      READ(IFILE,*) TMC
      READ(IFILE,*) TTROT
      READ(IFILE,*) NHISTBOXES
      READ(IFILE,*) BETA,BETAP
      READ(IFILE,*) I_P
      READ(IFILE,*) DBETA
      READ(IFILE,*) NTAY
      READ(IFILE,*) TRHOIJ
      READ(IFILE,*) RHOEPSILON
      READ(IFILE,*) TREADRHO
      READ(IFILE,*) I_HMAX
      READ(IFILE,*) I_VMAX
      IF(I_HMAX.EQ.-10) THEN
         IF(I_VMAX.GT.10) STOP 'I_HMAX=-10 and I_VMAX>10'
         READ(IFILE,*) ((NWHTAY(J,I),J=1,2),I=2,I_VMAX),G_VMC_PI,
     &      G_VMC_SEED,G_VMC_FAC
      ELSE
         READ(IFILE,*) NWHTAY(1,1),G_VMC_PI,G_VMC_SEED,G_VMC_FAC
      ENDIF
      WRITE(6,*) "GVMCFAC:",G_VMC_FAC
      IF(G_VMC_SEED.EQ.0) G_VMC_SEED=-7
      READ(IFILE,*) NPATHS,TNPDERIV
      READ(IFILE,*) TMONTE
      READ(IFILE,*) IMCSTEPS
      READ(IFILE,*) IEQSTEPS,BETAEQ   
      READ(IFILE,*) ILOGGING
      READ(IFILE,*) TBEGRAPH
      READ(IFILE,*) THFBASIS, THFCALC,NHFIT,HFMIX,HFEDELTA,
     &               TRHF,HFCDELTA
      READ(IFILE,*) IHFMETHOD,TREADHF,TREADTUMAT
      READ(IFILE,*) ETRIAL
      READ(IFILE,*) TMCDET,MDK(1),MDK(2),MDK(3),MDK(4),MDK(5)
      READ(IFILE,*) DETINV
      READ(IFILE,*) NDETWORK
      READ(IFILE,*) (ORBORDER2(I),I=1,8)
      DO I=1,8
C.. two ways of specifying open orbitals
C.. if orborder2(I,1) is integral, then if it's odd, we have a single
C.. open orbital
         IF(ORBORDER2(I).EQ.INT(ORBORDER2(I))) THEN
            ORBORDER(I,1)=IAND(INT(ORBORDER2(I)),65534)
            IF((INT(ORBORDER2(I))-ORBORDER(I,1)).GT.0) THEN
C.. we have an open orbital
               ORBORDER(I,2)=2
            ELSE
               ORBORDER(I,2)=0
            ENDIF
         ELSE
C.. non-integral.  The integral part is the number of closed oribtals,
C.. and the fractional*1000 is the number of open orbitals.
C.. e.g. 6.002 would mean 6 closed and 2 open
C.. which would have orborder(I,1)=6, orborder(I,2)=4
C.. but say 5.002 would be meaningless as the integral part must be a
C.. multiple of 2
            ORBORDER(I,1)=INT(ORBORDER2(I)+0.000001)
            ORBORDER(I,2)=INT((ORBORDER2(I)-ORBORDER(I,1)+0.000001)
     &                        *1000)*2
         ENDIF
      ENDDO
      READ(IFILE,*) NFROZEN,NTFROZEN
      READ(IFILE,*) TSPECDET,TSPECDET2
      CALL MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')
      IF(TSPECDET) THEN
         IF(TSPECDET2) THEN
            READ(IFILE,*) (SPECDET(I),I=1,NEL-NFROZEN)
         ELSE
            SPECDET(1)=0
         ENDIF
      ELSE
         READ(IFILE,*)
      ENDIF
      IF(IFILE.EQ.1) CLOSE(1)
      END

