      SUBROUTINE READINPUT(FILENAME,
     &   NEL,TSPN,LMS,TPARITY,IPARITY,NMAXX,NMAXY,
     &   NMAXZ,
     &	NMSH,BOX,BOA,COA,TALPHA,ALPHA,ISTATE,NEVAL,B2L,NCYCLE,NBLK,
     &	NKRY,TRHOOFR,TREAD,TCORR,TFODM,IOBS,JOBS,KOBS,TENERGY,TMC,
     &	NHISTBOXES,BETA,I_P,TRHOIJ,RHOEPSILON,TREADRHO,NPATHS,TMONTE,
     &  IMCSTEPS,IEQSTEPS,TTROT,I_HMAX,TCALCHMAT,NTAY,NWHTAY,I_VMAX,
     &  ILOGGING,TUEG,TBEGRAPH,THFBASIS,THUB,UHUB,BHUB,ETRIAL,TNPDERIV,
     & BETAP,BETAEQ,TMCDET,MDK,DETINV,THFCALC,NHFIT,TTILT,ITILTX,ITILTY,
     & TBLOCK,HFMIX,HFEDELTA,TREADINT,TPBC,TREAL,NDETWORK,ORBORDER,
     &   NFROZEN,NTFROZEN,TRHF,HFCDELTA,IHFMETHOD,TREADHF,TSPECDET,
     &  IP_SPECDET,DBETA,TREADTUMAT,THFORDER,TCSF,STOT,ICILEVEL,FUEGRS,
     &   TCPMD,ios,TDFREAD,fRc)
      USE input
      USE HElement
      USE UMatCache
      IMPLICIT REAL*8    (A-H,O-Z)
      INTEGER IPARITY(5)
      LOGICAL TSPN,TALPHA,TPARITY,TCORR,TFODM,TCALCHMAT,TUEG,     
     &   TREAD,TRHOOFR,TENERGY,TMC,TRHOIJ,TREADRHO,TMONTE,TTROT,  
     &   TBEGRAPH,THFBASIS,THUB,TNPDERIV,TMCDET,THFCALC,TTILT,    
     &   TBLOCK,TREADINT,TPBC,TREAL,TRHF,TREADHF,TSPECDET,TREADTUMAT,  
     &   THFORDER,TCSF,TBETAP, eof,TCPMD,TDFREAD,check
      CHARACTER*255 INP
      INTEGER SPECDET(*)
!
      CHARACTER(LEN=32) title  
      CHARACTER(LEN=100) w
      INTEGER NPATHS,STOT,z,zz,k,l
      INTEGER MDK(5),DETINV,ORBORDER(8,2)
      INTEGER NWHTAY(3,10)
      INTEGER vv,cc,kk
      INTEGER PREIH
      REAL*8  ORBORDER2(8)
      POINTER (IP_SPECDET,SPECDET)
      CHARACTER*255 FILENAME
      INTEGER ios
      REAL*8 DBETA(3)
      REAL*8 fRc
      INTEGER NTAY(2)
      INCLUDE 'vmc.inc'
      INCLUDE 'fcoul.inc'
      INCLUDE 'sym.inc'
      INCLUDE 'uhfdet.inc'
      INCLUDE 'time.inc'
      include 'calcp.inc'
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
      G_VMC_EXCITWEIGHT(:)=0.D0
      TPREVAR=.false.
      TVARCALC=.false.              !This is whether to calculate the expected variance for a MC run when doing full sum (seperate denominator and numerator at present)
      USEVAR=0
      g_VMC_ExcitWeights=0.D0
      EXCITFUNCS=.false.
      EXCITFUNCS(10)=.true.
      TGRIDVAR=.false.
      TLINEVAR=.false.
      MEMSAV=.false.
      TUPOWER=.false.
      IP_SPECDET=0
      LMS=0
      STOT=0
      TSPN=.false.
      TCSF=.false.

      BETAP=1.D-4
!DBETA now has  two elements
! DBETA(1) is DBETA
! DBETA(2) is GRAPHEPSILON
! DBETA(3) is PGENEPSILON
      DBETA(1)=-1.D0
      DBETA(2)=0.D0
      DBETA(3)=0.D0
      RHOEPSILON=1.D-6
      TTROT=.true.
      PREWEIGHTEPS=1.D-8
      
      TMPTHEORY=.FALSE.
      TNEWEXCITATIONS=.FALSE.

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
      fRc=0.D0
      do I = 1,5
        MDK(I) = 0
      end do
      DETINV = 0
      NDETWORK = 50000
      TSPECDET = .false.

      NSLOTSINIT=1024
      I_VMAX=0
      FCOULDAMPBETA=-1.D0
      COULDAMPORB=0

      !Precalc defaults
      PRE_IVMAX=0      
      PRE_TAY=0
      PRE_TAYLOG=.false.
      PRE_TAYREAL=0.D0
      PRE_TAYREAL(2,:)=0.1      !TOLERANCE   
      TRUECYCLES=0
      TOTALERROR=0.D0
      calcp_sub2vstar=.false.
!     Values for old parameters.
      TVVDISALLOW=.false.
      TRHOOFR = .false.
      TCORR = .false.
      TFODM = .false.
      TMC = .false.
      NHISTBOXES = 0
      NTAY(1) = 1
      NTAY(2) = 1
      TREADRHO = .false.
      TRHOIJ = .false.
      TBEGRAPH = .false.
      TDFREAD=.false.

      IOBS=16
      JOBS=16
      KOBS=16
      NMSH=32
      OrbECutoff=-1e50
      iGlobalTimerLevel=40
      ir=1
      IF(FILENAME.NE.'') THEN
         WRITE(6,*) "Reading from file: ",TRIM(FILENAME)
         OPEN(1,FILE=FILENAME,STATUS="OLD",err=99,iostat=ios)
      ELSEIF(IARGC().GT.0) THEN
         CALL GETARG(1,INP)
         WRITE(6,*) "Reading from file: ",INP
         OPEN(1,FILE=INP,STATUS="OLD")
      ELSE
         ir=5
         WRITE(6,*) "Reading from STDIN"
      ENDIF

!     --------------------------------------------
      write (6,'(/,64("*"),/)')
      call input_options(echo_lines=.true.,skip_blank_lines=.true.)
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
        case("SYSTEM")
            call readu(w)
            select case(w)
            case("DFREAD")
                TREADINT = .true.
                TDFREAD = .true.
                call readu(w)
                select case(w)
                case("ORDER")
                    THFORDER = .true.
                end select
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
            case("CPMD")
                TCPMD = .true.
                call readu(w)
                select case(w)
                case("ORDER")
                    THFORDER = .true.
                end select
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
                  do I = 1,4
                    call geti(IPARITY(I))
                  end do
C.. the last number is the symmetry specification - and is placed in position 5
                  IPARITY(5)=IPARITY(4)
                  IPARITY(4)=0
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
              case("USEBRILLOUINTHEOREM")
                TUSEBRILLOUIN=.TRUE. 
              case("RS")
                  call getf(FUEGRS)
              case("EXCHANGE-CUTOFF")
                  iPeriodicDampingType=2
                  if(item.lt.nitems) then
                     call getf(fRc)
                  endif
              case("EXCHANGE-ATTENUATE")
                  iPeriodicDampingType=1
                  if(item.lt.nitems) then
                     call getf(fRc)
                  endif
              case("COULOMB")
                  call getf(FCOUL)
              case("COULOMB-DAMPING")
                      call report("Coulomb damping feature removed"  
     &                ,.true.)
                  
C                  call readu(w)
C                  select case(w)
C                  case("ENERGY")
C                     call getf(FCOULDAMPMU)
C                     call getf(FCOULDAMPBETA)
C                  case("ORBITAL")
C                     call geti(COULDAMPORB)
C                     call getf(FCOULDAMPBETA)
C                  end select
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
              case("ENERGY-CUTOFF")
                  call getf(OrbECutoff)
              case("ENDSYS") 
                  exit system
              case default
                  call report("Keyword "         
     &              //trim(w)//" not recognized in SYSTEM block",.true.)
              end select
            end do system
            if(NEL.eq.0)
     &         call report("Number of electrons cannot be zero.",.true.)
            if(THUB.OR.TUEG.OR..NOT.(TREADINT.OR.TCPMD)) then
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
        case("PRECALC")
            preIV_MAX=1
            preIH=-10
            precalc: do
               call read_line(eof)
               if (eof) then
                   call report("No arguments for PRECALC block",.true.)
               end if
               call readu(w)
               select case(trim(w))
               case("VERTEX")
                   preIV_MAX=preIV_MAX+1
                   pre_TAY(3,preIV_MAX)=preIV_MAX
                   call inpgetprecalc(pre_TAY(1,preIV_MAX))
                
               case("TOLERANCE")
                   call readf(pre_TAYREAL(2,preIV_MAX))
                   if ( pre_TAYREAL(2,preIV_MAX) .le. 0.D0 ) then
                       call report(trim(w)//" must take a value "
     &                  //"more than 0.D0",.true.)
                   endif
               case("LINEVAR")
                   TLINEVAR(preIV_MAX)=.true.
                   IF(pre_TAYLOG(4,preIV_MAX)) THEN
                       call report(trim(w)//" can not be applied"
     &                 //" in conjunction with USEIMPORT",.true.)
                   ENDIF
                   do zz=1,3
                       call readf(LINEVARPAR(preIV_MAX,zz))
                   enddo
               case("GRIDVAR")
                   TGRIDVAR(preIV_MAX)=.true.
                   IF(pre_TAYLOG(1,preIV_MAX)) THEN
                       call report(trim(w)//" can not be applied"
     &                      //" in conjunction with USE",.true.)
                   ENDIF
                   do zz=1,6
                       call readf(GRIDVARPAR(preIV_MAX,zz))
                   enddo
               case("TOTALERROR")
                   call readf(TOTALERROR)
               case("TRUECYCLES")
                   call readi(TRUECYCLES)
               case("PREGRAPHEPSILON")
                   call readf(PREWEIGHTEPS)
               case("USE")
                    pre_TAYLOG(1,preIV_MAX)=.true.    
                    IF(TGRIDVAR(preIV_MAX)) THEN
                        call report(trim(w)//" can not be applied"
     &                      //" in conjunction with GRIDVAR",.true.)
                    ENDIF
                    !USE is specified on its own - all vertex levels take this value
                    IF(nitems.eq.1) THEN
                        do k=2,10
                            IF(k.ne.preIV_MAX) USEVAR(k,:)=0
                            IF(k.eq.preIV_MAX) THEN
                                do z=1,8
                                    USEVAR(preIV_MAX,z)=z+1
                                ENDDO
                            ENDIF
                        enddo
                    ELSE
                    call readi(USEVAR(preIV_MAX,1))
                    do k=2,8
                        if(item.lt.nitems) THEN
                            call readi(USEVAR(preIV_MAX,k))
                    
                        do z=1,(preIV_MAX-1)
                            do zz=1,8
                           if (USEVAR(z,zz).eq.USEVAR(preIV_MAX,k)) then
                           call report(trim(w)//" can only be applied"
     &                      //" to a vertex level once",.true.)
                            endif
                            enddo
                        enddo
                    
                        endif
                    enddo
                    ENDIF
               case("MEMORISE")
                   MEMSAV(preIV_MAX)=.true.
                   if ((pre_TAY(1,preIV_MAX).ne.-7).and.
     &              (pre_TAY(1,preIV_MAX).ne.-19)) then
                        call report(trim(w)//" not yet valid for full"
     &                  //" sum precalc",.true.)
                   endif
               case("NONE")
                    pre_TAY(3,preIV_MAX)=0
               case("CYCLES")
                   call readi(pre_TAY(2,preIV_MAX))
                   if ((pre_TAY(1,preIV_MAX).ne. -7).and.
     &                  (pre_TAY(1,preIV_MAX).ne. -19)) then
                      call report(trim(w)//" only valid for MC"
     &                 //" method",.true.)
                   end if
               case("UEPSILON")
                   call readf(pre_TAYREAL(1,preIV_MAX))
               case("FINDC")
                   pre_TAYLOG(2,preIV_MAX)=.true.
               case("FINDD")
                   pre_TAYLOG(5,preIV_MAX)=.true.
               case("USED")
                   pre_TAYLOG(6,preIV_MAX)=.true.
               case("FINDIMPORT")
                   if ( preIV_MAX .lt. 3 ) then
                       call report(trim(w)//" only valid for vertex "
     &                  //"levels of 3 or higher",.true.)
                   endif
                   pre_TAYLOG(3,preIV_MAX)=.true.
               case("USEIMPORT")
                   if ( preIV_MAX .lt. 3 ) then
                       call report(trim(w)//" only valid for vertex "
     &                  //"levels of 3 or higher",.true.)
                   endif
                   pre_TAYLOG(4,preIV_MAX)=.true.
                       do zz=1,(preIV_MAX-1)
                          if ( pre_TAYLOG(4,zz) ) then
                             call report(trim(w)//" can only be applied"
     &                       //" to a single vertex level",.true.)
                         endif
                       enddo
               case("ENDPRECALC")
                   exit precalc
               case default
                   call report ("Keyword "//trim(w)//
     &              " not recognised",.true.)
               end select
               end do precalc
                   
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
                    call inpgetmethod(NWHTAY(1,I_VMAX),NWHTAY(2,I_VMAX),
     &               I_VMAX)
                  case("EXCITATIONS")
                     call readu(w)
                     call inpgetexcitations(NWHTAY(2,I_VMAX),w,I_VMAX)
                  case("CYCLES")
                     call readi(NWHTAY(2,I_VMAX))
                     if ( NWHTAY(1,I_VMAX).ne. -7.and.
     &                    NWHTAY(1,I_VMAX).ne.-19 ) then
                        call report(trim(w)//" only valid for MC " 
     &                   //"method",.true.)
                     end if
                  case("VERTICES")
                     call geti(NWHTAY(3,I_VMAX))
                  case("MULTIMCWEIGHT")
                     call getf(g_MultiWeight(I_VMAX))
                  case("CALCVAR")
                      if ( NWHTAY(1,I_VMAX).NE.-20 ) then
                          call report("Keyword "//trim(w)//" 
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
                  if ( I_HMAX .ne. -7.and.
     &                 I_HMAX .ne. -19) then
                      call report(trim(w)//" only valid for MC " 
     &                   //"method",.true.)
                  end if
              case("VVDISALLOW")
                  TVVDISALLOW=.TRUE.
              case("MCDIRECTSUM")
                  TMCDIRECTSUM=.TRUE.
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
              case("UFORM-POWER")
                  TUPOWER=.true.
              case("CHEMPOTWEIGHTING")
                  call readf(g_VMC_ExcitWeights(1,1))
                  call readf(g_VMC_ExcitWeights(2,1))
                  call readf(G_VMC_EXCITWEIGHT(1))
                  DO l=1,5
                    IF(EXCITFUNCS(l)) THEN
                        call report(trim(w)//" only valid if "
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
                        call report(trim(w)//" only valid if "
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
                        call report(trim(w)//" only valid if "
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
                        call report(trim(w)//" only valid if "
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
                        call report(trim(w)//" only valid if "
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
                     NPATHS=-2
                  case default
                     read(w,"(I)") NPATHS
                  end select
              case("ALLPATHS")
                  NPATHS = -1
              case("DERIV")
                  TNPDERIV = .true.
                 if (DBETA(1) .lt. 0 ) then
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
                  call getf(DBETA(1))
              case("RHOEPSILON")
                  call getf(RHOEPSILON)
              case("GRAPHEPSILON")
                  call getf(DBETA(2))
              case("PGENEPSILON")
                  call getf(DBETA(3))
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
              case("RHO-1STORDER")
                  NTAY(2)=4
              case("FOCK-PARTITION")
                  NTAY(2)=2
              case("FOCK-PARTITION-LOWDIAG")
                  NTAY(2)=3
              case("FOCK-PARTITION-DCCORRECT-LOWDIAG")
                  NTAY(2)=5
              case("DIAG-PARTITION")
                  NTAY(2)=1
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
                  case("MODIFIED")
                     IHFMETHOD=3
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
                  if ( mod(NFROZEN,2).ne.0 .or.
     &             (NTFROZEN.GT.0 .and. mod(NTFROZEN,2).ne.0) ) then
                      call report("NFROZEN and (+ve) NTFROZEN must be"  
     &                //"multiples of 2",.true.)
                  end if
                  if (
     &             (NTFROZEN.LT.0 .and. mod(NEL-NTFROZEN,2).ne.0) ) then
                      call report("-ve NTFROZEN must be same parity  "
     &                //"as NEL",.true.)
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
              case("UMATCACHE")
                  call geti(NSLOTSINIT)
              case("NOUMATCACHE")
                  NSLOTSINIT=-1
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
                        CALL report("Logging keyword FMCPR "//trim(w)
     &                     //" not recognised",.true.)
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
                        CALL report("Logging keyword CALCPATH "//trim(w)
     &                     //" not recognised",.true.)
                     end select
                  enddo
              case("XIJ")
                  ILOGGING = IOR(ILOGGING,2**6)
              case("HAMILTONIAN")
                  ILOGGING = IOR(ILOGGING,2**7)
              case("PSI")
                  ILOGGING = IOR(ILOGGING,2**8)
              case("TIMING")
                  call readi(iGlobalTimerLevel)
              case("VERTEX")
                  do while(item.lt.nitems)
                     call readu(w)
                     select case(w)
                    ! case("1000")
                    !     ILOGGING = IOR(ILOGGING,2**9)
                     case("EVERY")
                         ILOGGING = IOR(ILOGGING,2**10)
                     case default
                        READ(w,"I") G_VMC_LOGCOUNT
                        ILOGGING = IOR(ILOGGING,2**9)
                      !  CALL report("Logging keyword VERTEX "//trim(w)
     &                !     //" not recognised",.true.)
                     end select
                  end do
              case("HFBASIS")
                  ILOGGING = IOR(ILOGGING,2**11)
              case("ENDLOG")
                  exit logging
              case default
                 CALL report("Logging keyword "//trim(w)
     &                     //" not recognised",.true.)
              end select
            end do logging
        case("END")
            exit
!       Error message if unknown keyword
        case default
            call report ("Keyword "//trim(w)//" not recognized",.true.)
        end select
      end do main
      write (6,'(/,64("*"),/)')
C.. We still need a specdet space even if we don't have a specdet.
      IF(IP_SPECDET.EQ.0) call MEMORY(IP_SPECDET,NEL-NFROZEN,'SPECDET')
C..   Testing ILOGGING
C     ILOGGING = 0771
      IF(I_VMAX.EQ.0.AND.NPATHS.NE.0) 
     &   STOP 'NPATHS!=0 and I_VMAX=0.  VERTEX SUM max level not set'
      WRITE (6,"(A,Z4)") 'ILOGGING after input routine', ILOGGING
      IF(IR.EQ.1) CLOSE(1)
   99 IF (ios.gt.0) THEN
          WRITE (6,*) 'Problem reading input file ',TRIM(FILENAME)
      END IF
      
      !Make sure there aren't more precalc levels than true vertex levels - problems otherwise
      IF(preIV_MAX.gt.I_VMAX) THEN
          CALL report("There cannot be more precalc vertex levels "
     &     //"than vertex levels in the main program",.true.)
      ENDIF
      
      !We make sure that in precalc, a use statement is not specified more than once for any vertex level
      IF(preIV_MAX.ne.0) THEN
          do vv=2,I_VMAX
            check=.false.
            do kk=2,preIV_MAX
                do cc=1,8
                    IF((USEVAR(kk,cc).eq.vv).and.(check)) THEN
                     CALL report("Can only specify to use precalc "
     &               //"parameters on a given vertex level once",.true.)
                    ENDIF
                    IF(USEVAR(kk,cc).eq.vv) check=.true.
                enddo
            enddo
            !If use isn't specified for a vertex level, use the values given in the input file
!           Done later now
!            IF(.not.check) THEN
!                g_VMC_ExcitWeights(:,vv)=g_VMC_ExcitWeights(:,1)
!                G_VMC_EXCITWEIGHT(vv)=G_VMC_EXCITWEIGHT(1)
!            ENDIF
        enddo
      ENDIF
      
      !If not doing precalc, set all weighting parameters to the ones in the input file
      IF(preIV_MAX.eq.0) THEN
        do vv=2,I_VMAX
            g_VMC_ExcitWeights(:,vv)=g_VMC_ExcitWeights(:,1)
            G_VMC_EXCITWEIGHT(vv)=G_VMC_EXCITWEIGHT(1)
        enddo
      ENDIF

      !IF THERE IS NO WEIGHTING FUNCTION
      do vv=1,9
          IF(EXCITFUNCS(vv)) EXCITFUNCS(10)=.false.
      enddo
      !IF FINDD or USED specified without using Excitweighting option
      do vv=2,preIV_MAX
         IF((pre_TAY(1,vv).eq.-20).and.((NWHTAY(1,vv).eq.-7).or.        &
     &    (NWHTAY(1,vv).eq.-19))) THEN
         CALL report("Full precalc cannot be used on a vertex level"
     &   //" which is only sampled using MC in the main program",.true.)
          ENDIF
          IF((pre_TAYLOG(5,vv).or.pre_TAYLOG(6,vv)).and.                &
     &         (.not.EXCITFUNCS(1))) THEN
               CALL report("Logging keyword FINDD and USED"
     &         //" can only be used with EXCITWEIGHTING",.true.)
          ENDIF
          IF(TGRIDVAR(vv).and.((.not.EXCITFUNCS(1)).and.(.not.          &
     &         EXCITFUNCS(4)))) THEN
               CALL report("GRIDVAR option only available with" 
     &          //" excitation functions with two variables",.true.)
          ENDIF
          IF(TLINEVAR(vv).and.((.not.pre_TAYLOG(3,vv)).and.             &
     &              (.not.pre_TAYLOG(2,vv)))) THEN
              CALL report("LINEVAR option only available with"
     &          //" FINDIMPORT or FINDC",.true.)
          ENDIF
      ENDDO
      IF((preIV_MAX.ne.0).AND.(.NOT.TMCDIRECTSUM)) THEN
          CALL report("Precalculation can only work with the"
     &          //" MCDIRECTSUM option enabled",.true.)
      ENDIF
      IF((TOTALERROR.ne.0.D0).AND.(TRUECYCLES.ne.0)) THEN
          CALL report("Only TRUECYCLES or TOTALERROR can be"            &
     &      //" specified in precalc block",.true.)
      ENDIF
      IF((TOTALERROR.ne.0.D0).AND.(preIV_MAX.ne.I_VMAX)) THEN
          CALL report("TOTALERROR can only be used if the precalc"      &
     &    //" levels are equal to the main block vertex levels",.true.)
      ENDIF
      
      RETURN
      END

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
                        call report("Keyword error with EXCITATIONS "
     &                     //trim(w),  
     &                          .true.)
                  end select
      end
      subroutine inpgetprecalc(preIH)
         use input
         implicit none
         include 'uhfdet.inc'
         integer preIH
         CHARACTER(LEN=16) w
                do while ( item .lt. nitems )
                  call readu(w)
                  select case(w)
                  case("HDIAG")
                      call readu(w)
                      select case(w)
                      case("FULL")
                          preIH=-20
                      case("MC")
                          preIH=-19
                          IF(.NOT.TUSEBRILLOUIN) THEN
                              write(6,*) "Warning  USEBRILLOUINTHEOREM"
     &                       //" might need to be specified in system "
     &                       //"block to use MC-PRECALC"
                          ENDIF
                      case default
                         call report("Error - must specify FULL"
     &                   //" or MC after HDIAG in PRECALC block",.true.)
                      end select
                   case("RHODIAG")
                         call readu(w)
                         select case(w)
                         case("FULL")
                             preIH=-8
                         case("MC")
                             preIH=-7
                         case default
                           call report("Error - must specify FULL or "
     &                     //"MC after RHODIAG in PRECALC block",.true.)
                         end select
                     case default
                           call report("Keyword error with "//trim(w),
     &                     .true.)
                     end select
                 end do
        end
                   
                  
         
      subroutine inpgetmethod(I_HMAX,NWHTAY,I_V)
         use input
         implicit none
         include 'vmc.inc'
         include 'calcp.inc'
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
                            case default
                                call report("Error - must specify OLD" 
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
                              case("H0")
                                  NWHTAY=IBSET(NWHTAY,5)
                                  if(I_HMAX.ne.-21) call report ("H0 "
     &                       //"can only be specified with POLY... NEW")
                              case default
                                call report("Error - must specify DIAG" 
     &                        //" or POLY vertex star method",.true.)
                               end select
                           enddo
                           if(i_hmax.eq.0)
     &                   call report("OLD/NEW not specified for STAR",
     &                          .true.)
                        case default
                        call report("Keyword error with "//trim(w),  
     &                          .true.)
                        end select
                    case default
                        call report("Error.  Method not specified."    
     &                    //" Stopping.",.true.)
                    end select
               end do
      end

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
      INTEGER NWHTAY(3,10)
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

