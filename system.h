C     ==================================================================
C     == MAXSP  : MAXIMUM NUMBER OF ATOMIC SPECIES                    ==
C     == NHX    : MAXIMUM NUMBER OF GAUSS-HERMIT POINTS               ==
C     ==        : MAXIMUM NUMBER OF BETA FUNCTIONS (VANDERBILT)       ==
C     ==        : MAXIMUM NUMBER OF PROJECTORS FOR NON-LOCAL PP       ==
C     == NACX   : SIZE OF THE ACCUMULATOR ARRAY                       ==
C     == MAXDIS : MAX. NUMBER OF VECTORS IN DIIS                      ==
C     == MXGDIS : MAX. NUMBER OF VECTORS IN GDIIS                     ==
C     == NBRX   : NUMBER OF DISTINCT RADIAL BETA FUNCTIONS            ==
C     == LX     : 2*LMAX-1 WHERE LMAX=3 IF S,P AND D FUNCTIONS        ==
C     == NLX    : COMBINED ANGULAR MOMENTUM FOR LLi (VANDERBILT)      ==
C     == MX     : 2*LX-1                            (VANDERBILT)      ==
C     == LIX    : MAX. LLi                          (VANDERBILT)      ==
C     == MIX    : 2*LIX-1                           (VANDERBILT)      ==
C     == LMAXX  : MAXIMUM NUMBER OF ANGULAR MOMENTUM                  ==
C     == LM1X   : LMAXX-1                                             ==
C     ==================================================================
      INTEGER   MAXSP,NHX,NACX,MAXDIS,MXGDIS,NBRX,LX,NLX,MX,LIX,MIX,
     &          LMAXX,LM1X
      PARAMETER (MAXSP=20)
      PARAMETER (NHX=100)
      PARAMETER (NACX=50)
      PARAMETER (MAXDIS=20)
      PARAMETER (MXGDIS=20)
      PARAMETER (NBRX=6)
      PARAMETER (LX=5)
      PARAMETER (NLX=9,MX=2*LX-1,LIX=3,MIX=LIX*2-1)
      PARAMETER (LMAXX=4,LM1X=LMAXX-1)
C     ==================================================================
C     == MAXSYS : Maximum value for system                            ==
C     ==          NSX   Number of atomic species                      ==
C     ==          NAX   Maximum number of atoms per species           ==
C     ==          NCORX (NSX*(NSX+1))/2 (detsp.F)                     ==
C     ==          NHXS  Maximum number for all species as NHX         ==
C     ==          LPMAX Maximum value of angular momentum             ==
C     ==          MMAXX Maximum number of spline points               ==
C     ==--------------------------------------------------------------==
      INTEGER        NSX,NAX,NSAX,NCORX,NHXS,LPMAX,MMAXX
      COMMON/MAXSYS/ NSX,NAX,NSAX,NCORX,NHXS,LPMAX,MMAXX
C     ==================================================================
C     == SPAR   : Global array dimension (parallel work)              ==
C     ==          NHGS = sum of NHG for all processors                ==
C     ==          NHGLS= sum of NHGL     --                           ==
C     ==          NHGKS= sum of NHGK     --                           ==
C     ==          NGWS = sum of NGW      --                           ==
C     ==          NGWKS= sum of NGWKS    --                           ==
C     ==          NGWLS= sum of NGWL     --                           ==
C     ==          NR1S,NR2S,NR3S total mesh dimension                 ==
C     ==================================================================
      INTEGER        NHGS,NHGLS,NHGKS,NGWS,NGWKS,NGWLS,
     &               NLDXS,NR1S,NR2S,NR3S
      COMMON/SPAR/   NHGS,NHGLS,NHGKS,NGWS,NGWKS,NGWLS,
     &               NLDXS,NR1S,NR2S,NR3S
C     ==--------------------------------------------------------------==
C     == NCPW   : Number of plane waves                               ==
C     ==          NHG  number of G-components for electronic density  ==
C     ==          NHGL number of G-shells for electronic density      ==
C     ==          NGW  number of G-components for wavefunctions       ==
C     ==          NGWL number of G-shells for wavefunctions           ==
C     ==--------------------------------------------------------------==
      INTEGER        NHG,NHGL,NGW,NGWL
      COMMON/NCPW/   NHG,NHGL,NGW,NGWL
C     ==--------------------------------------------------------------==
C     == NKPT   : Used with K-points options                          ==
C     ==          NHGK = 2 * NHG (no inversion symmetry)              ==
C     ==          NGWK = 2 * NGW (no inversion symmetry)              ==
C     ==          NKPNT - number of k points  (in memory)             ==
C     ==          NKPTS - number of total k points >= NKPNT           ==
C     ==          NBLKP - number of k points block                    ==
C     ==          NKPBL(NBLKP) - number of k points for each block    ==
C     ==          KPBEG(NBLKP) - number - 1 of first k-pt in block    ==
C     ==--------------------------------------------------------------==
      INTEGER        NHGK,NGWK,NKPNT,NKPTS,NBLKP
      COMMON/NKPT/   NHGK,NGWK,NKPNT,NKPTS,NBLKP
      INTEGER        NKPBL(*),KPBEG(*)
      POINTER        (IP_NKPBL,NKPBL),(IP_KPBEG,KPBEG)
      COMMON/NKIP/   IP_NKPBL,IP_KPBEG
C     ==--------------------------------------------------------------==
C     == PARM   : Information about supercell and mesh                ==
C     ==          ALAT lattice parameter                              ==
C     ==          A1(3), A2(3), A3(3) crystal vector basis            ==
C     ==          OMEGA Supercell volume                              ==
C     ==          TPIBA=2.D0*PI/ALAT                                  ==
C     ==          TPIBA2=TPIBA*TPIBA                                  ==
C     ==          APBC(4) used by PBC routine                         ==
C     ==                  for periodic boundary condition             ==
C     ==          IBRAV Determines shape and symmetry of the supercell==
C     ==              =0 Isolated system (cubic cell)                 ==
C     ==               1 Simple cubic                                 ==
C     ==               2 Face centred cubic                           ==
C     ==               3 Body centred cubic                           ==
C     ==               4 Hexagonal cell                               ==
C     ==               5 Trigonal or rhombohedral cell                ==
C     ==               6 Tetragonal                                   ==
C     ==               7 Body tetragonal                              ==
C     ==               8 Orthorhombic                                 ==
C     ==              12 Monoclinic                                   ==
C     ==              14 Triclinic                                    ==
C     ==          NR1, NR2, NR3 mesh dimension for each processor     ==
C     ==--------------------------------------------------------------==
      REAL*8         ALAT,A1(3),A2(3),A3(3),OMEGA,TPIBA,TPIBA2,APBC(4)
      INTEGER        IBRAV,NR1,NR2,NR3
      COMMON/PARM/   ALAT,A1,A2,A3,OMEGA,TPIBA,TPIBA2,APBC,
     *               IBRAV,NR1,NR2,NR3
C     ==--------------------------------------------------------------==
C     == FPAR   : Leading dimensions for some arrays,especially FFT   ==
C     ==          KR1, KR2, KR3 FFT mesh dimension for each proc.     ==
C     ==          KR1S, KR2S, KR3S FFT mesh for all proc.             ==
C     ==          KRX Used in FFTPRP routine with groups              ==
C     ==          KRY, KRZ UNUSED                                     ==
C     ==          NNR1 Number of real components for density          ==
C     ==               for each processor                             ==
C     ==--------------------------------------------------------------==
      INTEGER        KR1,KR2,KR3,KR1S,KR2S,KR3S,KRX,KRY,KRZ,NNR1
      COMMON/FPAR/   KR1,KR2,KR3,KR1S,KR2S,KR3S,KRX,KRY,KRZ,NNR1
C     ==--------------------------------------------------------------==
C     == ACCU   : Accumulator quantities                              ==
C     ==          NACC number of quantities calculated per self-c iter==
C     ==          ACC(1:NACC) values of calculated quantities         ==
C     ==--------------------------------------------------------------==
      INTEGER        NACC
      REAL*8         ACC
      COMMON/ACCU/   ACC(NACX),NACC
C     ==================================================================
C     == ALL THE INFO FOR PARALLEL                                    ==
C     ==================================================================
C     ==--------------------------------------------------------------==
C     == MAXCPU  : Maximum number of CPU
C     ==--------------------------------------------------------------==
      INTEGER       MAXCPU
      PARAMETER     (MAXCPU=512)
C     ==--------------------------------------------------------------==
C     == PARENT  : .TRUE if IP = SOURCE                               ==
C     ==--------------------------------------------------------------==
      LOGICAL       PARENT
      COMMON/PARAL/ PARENT
C     ==--------------------------------------------------------------==
C     == NCPUS   : value of NCPUS environment variable                ==
C     == NPROC   : TOTAL NUMBER OF PROCESSORS                         ==
C     == ME      : CPU index                                          ==
C     == MEPOS   : equal to ME                                        ==
C     == SOURCE  : MASTER CPU index                                   ==
C     == IGEQ0   : PROCESSOR INDEX OF G=0 COMPONENT                   ==
C     == ALLGRP  : GIVEN VALUE FOR BROADCAST COMMUNICATION            ==
C     == NHRAYS  : number of rays for the processor for the density   ==
C     == NGRAYS  : number of rays for the processor for the wavefunc. ==
C     ==--------------------------------------------------------------==
      INTEGER       NCPUS,NPROC,ME,MEPOS,SOURCE,IGEQ0,ALLGRP,
     &              NHRAYS,NGRAYS
      COMMON/PARAI/ NCPUS,NPROC,ME,MEPOS,SOURCE,IGEQ0,ALLGRP,
     &              NHRAYS,NGRAYS
C     ==--------------------------------------------------------------==
C     == NRXPL   : NR1 for each processor                             ==
C     == NRZPL   : Number of z-planes for each processor              ==
C     == SPARM   : Information about dimension of distributed data    ==
C     == NST12   : index of the first and last distributed states     ==
C     ==         : per proc                                           ==
C     == PGROUP  : Pointer position --> PE number                     ==
C     == NLINK   : Pointer PE number --> position - 1                 ==
C     ==--------------------------------------------------------------==
      INTEGER       NRXPL(0:MAXCPU,2),NRZPL(0:MAXCPU,2),
     &              SPARM(9,0:MAXCPU),NST12(0:MAXCPU,2),
     &              PGROUP(0:MAXCPU),NLINK(0:MAXCPU)
      COMMON/PARAP/ NRXPL,NRZPL,SPARM,NST12,PGROUP,NLINK
C     ==================================================================
C     == MAPPING G VECTORS-->PROCESSES                                ==
C     ==================================================================
C     == SUBDIVIDE THE PROCESSORS IN GROUPS (see GROUPS routine)      ==
C     == MAPGP(1:NHG) : Use for distribution of wavefunction          ==
C     == MAXGRP=MAXCPU : Maximum taskgroups (for distributed FFT)     ==
C     == NOGRP  : Number of groups                                    ==
C     == NPGRP  : Number of processors per groups                     ==
C     == MEOGRP : Index of the group for the processor                ==
C     == MEPGRP : Index of the processor in the group                 ==
C     == MPEN   : Dimensions used for groups (see FFTPRP routine)     ==
C     == MPENM  : ''                  ''          ''                  ==
C     == NOLIST(1:MAXGRP) : List of processors in the orbital group   ==
C     == NPLIST(1:MAXGRP) : List of processors in the plane wave group==
C     ==--------------------------------------------------------------==
      INTEGER       MAPGP(*)
      POINTER       (IP_MAPGP,MAPGP)
      COMMON/MAPP/  IP_MAPGP
      INTEGER       MAXGRP
      PARAMETER     (MAXGRP=MAXCPU)
      INTEGER       NOGRP,NPGRP,MEOGRP,MEPGRP,MPEN,MPENM,
     &              NOLIST(MAXGRP),NPLIST(MAXGRP)
      COMMON/GROUP/ NOGRP,NPGRP,MEOGRP,MEPGRP,MPEN,MPENM,NOLIST,NPLIST
C     ==================================================================
C     == MAPPING ATOMS-->PROCESSES                                    ==
C     ==================================================================
      INTEGER NATPE,NORBPE,IATPT,IPEPT,IATPE
      COMMON/ATPT/NATPE,NORBPE,IP_IATPT,IP_IPEPT,IP_IATPE
      DIMENSION IATPT(2,*),IPEPT(2,0:NPROC),IATPE(*)
      POINTER(IP_IATPT,IATPT),(IP_IPEPT,IPEPT),(IP_IATPE,IATPE)
C     ==================================================================
C     == THE CONTROL PARAMETERS                                       ==
C     ==================================================================
C     == MLOGO : NUMBER OF LOGICAL VARIABLES IN CNTL COMMON           ==
C     == MINTE : NUMBER OF INTEGER VARIABLES IN CNTI COMMON           ==
C     == MREAL : NUMBER OF REAL*8  VARIABLES IN CNTR COMMON           ==
C     ==--------------------------------------------------------------==
      INTEGER   MLOGO,MINTE,MREAL
      PARAMETER (MLOGO=84)
      PARAMETER (MINTE=37)
      PARAMETER (MREAL=33)
C     ==================================================================
C     == A LOT OF ITEMS IN THE COMMON BLOCKS IS NEEDED FOR PARALLEL   ==
C     ==================================================================
C     == CNTL1  : Use for the broadcast (dummy variable)              ==
C     == MD     : Molecular dynamics option                           ==
C     == GEOPT  : Geometry optimization                               ==
C     == WFOPT  : Wavefunction optimization                           ==
C     == TLOWD  : Lowdin orthogonalization instead of Gram-Schmidt    ==
C     == TRANE  : Randomize wavefunctions                             ==
C     == TRANP  : Randomize ionic coordinates                         ==
C     == TC     : Temperature control for electrons                   ==
C     == TCP    : Temperature control for ions                        ==
C     == ANNEI  : Simulated annealing for ions                        ==
C     == ANNEE  : Simulated annealing for electrons                   ==
C     == ANNEC  : Simulated annealing for cell                        ==
C     == QUENCHP: Quench ionic velocities                             ==
C     == QUENCHE: Quench electronic velocities (Car-Parrinello)       ==
C     == QUENCHB: Quench system to the Born-Oppenheimer surface       ==
C     == CNSTMD : Iterative orthogonalization (Car-Parrinello)        ==
C     == TSDE   : Steepest descent (electronic coordinates)           ==
C     == TSDP   : Steepest descent (ionic coordinates)                ==
C     == DIIS   : Wavefunction optimization by DIIS                   ==
C     == PREC   : Preconditioning for Steepest descent (electronic)   ==
C     == TSCALE : Use Scaled ionic coordinates                        ==
C     == BOHR   : Use atomic units for ionic coordinates (bohr)       ==
C     == PCG    : Wavefunction optimization by Precond. Conj. Gradient==
C     == GDIIS  : Geometry optimization by GDIIS                      ==
C     == RFO    : Geometry optimization by Rational Function Approx.  ==
C     == BFGS   : Geometry optimization by quasi-Newton update        ==
C     == TGC    : Gradient correction                                 ==
C     == TGCX   : Gradient correction for exchange part               ==
C     == TGCC   : Gradient correction for correlation part            ==
C     == TSMOOTH: Smoothing of the density                            ==
C     == CALDIP : Calculate Dipole dynamics                           ==
C     == TNOSEE : Nose-Hoover thermostats for electrons dynamics      ==
C     == TNOSEP : Nose-Hoover thermostats for ions dynamics           ==
C     == TNOSES : One thermostat chain per species                    ==
C     == KSENER : Calculate Kohn-Sham energies for the final potential==
C     == WCOMP  : Wavefunction in compressed form                     ==
C     == NONORT : Use nonorthogonal orbital                           ==
C     == BIGMEM : Big memory option                                   ==
C     == THARM  : Harmonic reference system integration               ==
C     == TMASS  : Scaled electron masses                              ==
C     == KRWFN  : Store real space representation of wavefunctions    ==
C     == PCGMIN : Perform quadratic line search                       ==
C     == PROPER : Properties calculation                              ==
C     == VIBRAT : Vibrational analysis                                ==
C     == TLSD   : Local Spin Density Approximation                    ==
C     == TPRES  : Stress tensor calculation                           ==
C     == SIMUL  : Combined geometry/wavefunction scheme               ==
C     == TPRCP  : Parrinello-Rahman-Car-Parrinello                    ==
C     == TNOSEC : Nose-Hoover thermostats for cell dynamics           ==
C     == QUENCHC: Quench electronic velocities (Parrinello-Rahman)    ==
C     == TSDC   : Cell dynamics by Steepest descent                   ==
C     == TRANC  : Randomize cell                                      ==
C     == TRFNL  : Calculate nonlocal pseudopotential in real space    ==
C     == TMEMCHK: Memory checking                                     ==
C     == TPENN  : NPT dynamics                                        ==
C     == TFINT  : Free Energy Functional (Alavi's method)             ==
C     == TDIAG  : Diagonalisation scheme (Lanczos or Davidson)        ==
C     == TDAVI  : Davidson diagonalisation                            ==
C     == TLANC  : Lanczos diagonalisation                             ==
C     == TFRSBLK: New Block Lanczos diagonalisation                   ==
C     == TDIPD  : Dipole Dynamics                                     ==
C     == TINTER : Interface to classical MD program                   ==
C     == TCC    : Temperature Control for the cell                    ==
C     == TEPOT  : Electrostatic Potential                             ==
C     == TEXPOT : External Potential                                  ==
C     == TPATH  : Path Integrals                                      ==
C     == TRESCALE : Re-adjust ionic velocities only after restart to  ==
C     ==          desired temperature TEMPW                           ==
C     == TSYMRHO: SYMMTRIZE DENSITY                                   ==
C     ==--------------------------------------------------------------==
      LOGICAL      CNTL1,MD,TMDBO,GEOPT,WFOPT,TLOWD,
     *             TRANE,TRANP,THEAD,TC,TCP,ANNEI,ANNEE,ANNEC,
     *             QUENCHP,QUENCHE,QUENCHB,CNSTMD,
     *             STOPOT,TCHNGW,TSDE,TSDP,DIIS,PREC,TSCALE,
     *             BOHR,TCHNGR,PCG,GDIIS,RFO,BFGS,TGC,TGCX,TGCC,TSMOOTH,
     *             CALDIP,TNOSEE,TNOSEP,TNOSES,TIMING,
     *             KSENER,ORBROT,RCOMP,WCOMP,NONORT,BIGMEM,GAUDYN,
     *             THARM,TMASS,FCNSTR,KRWFN,PCGMIN,POSVER,
     *             PROPER,VIBRAT,INITCM,FINALCM,TLSD,TPRES,SIMUL,
     *             TPRCP,TNOSEC,QUENCHC,TSDC,TRANC,TRFNL,
     *             TMEMCHK,TPENN,TFINT,TDIAG,TDAVI,TLANC,TFRSBLK,
     *             TDIPD,TINTER,TCC,TEPOT,TEXPOT,TPATH,TRESCALE,
     *             TSAMPL,TPREC,TFDIST,TSYMRHO
      COMMON/CNTL/ CNTL1,MD,TMDBO,GEOPT,WFOPT,TLOWD,
     *             TRANE,TRANP,THEAD,TC,TCP,ANNEI,ANNEE,ANNEC,
     *             QUENCHP,QUENCHE,QUENCHB,CNSTMD,
     *             STOPOT,TCHNGW,TSDE,TSDP,DIIS,PREC,TSCALE,
     *             BOHR,TCHNGR,PCG,GDIIS,RFO,BFGS,TGC,TGCX,TGCC,TSMOOTH,
     *             CALDIP,TNOSEE,TNOSEP,TNOSES,TIMING,
     *             KSENER,ORBROT,RCOMP,WCOMP,NONORT,BIGMEM,GAUDYN,
     *             THARM,TMASS,FCNSTR,KRWFN,PCGMIN,POSVER,
     *             PROPER,VIBRAT,INITCM,FINALCM,TLSD,TPRES,SIMUL,
     *             TPRCP,TNOSEC,QUENCHC,TSDC,TRANC,TRFNL,
     *             TMEMCHK,TPENN,TFINT,TDIAG,TDAVI,TLANC,TFRSBLK,
     *             TDIPD,TINTER,TCC,TEPOT,TEXPOT,TPATH,TRESCALE,
     *             TSAMPL,TPREC,TFDIST,TSYMRHO
C     ==================================================================
C     == CNTI1  : Use for the broadcast (dummy variable)              ==
C     == NOMORE : Maximum number of steps                             ==
C     == MAXIT  : Parameters for the RATTLE step                      ==
C     == MDIIS  : Maximum number of vectors retained for DIIS (Wave.) ==
C     == INSYS  : Logic unit for properties (by default=5)            ==
C     == MGDIIS : Size of GDIIS matrix (geometry optimization)        ==
C     == NPARA  : Initial Hessian (Unit, Disco or Schlegel)           ==
C     == NCHP   : Number of Nose thermostats (ions)                   ==
C     == NCHB   : Number of Nose thermostats (cell)                   ==
C     == NCHS   : Number of Nose thermostats (electrons)              ==
C     == NCALLS0: Number of Yoshida-Suzuki steps (Nose)               ==
C     == NIT0   : Number of integration cycles (Nose)                 ==
C     == NKSSTA : Number of Kohn-Sham eigenvalues                     ==
C     == IMOVIE : Write movie file every imovie steps                 ==
C     == IPROJ  : Electronic gradient projection to be used           ==
C     == INWFUN : Wavefunction initialization (1=random, 2=atomic)    ==
C     == NRESTF : Number of different restart files                   ==
C     == NTRANS : Use in Lanczos scheme                               ==
C     == NTRAJ  : Trajectories are saved on file every NTRAJ          ==
C     == NPRES  : Stress tensor calculation each NPRES                ==
C     == NSPLP  : Number of points for spline                         ==
C     == NSORDER: Saddle point order for geo. opt. by rational func.  ==
C     == RCOMPB : Give index of read wavefunction compression scheme  ==
C     == WCOMPB : Give index of written wave. compression scheme      ==
C     == NDAVV  : Davidson parameter (n. of iteratation for one diag. ==
C     == N_FRIES: Lanczos para. (n. of iteratation for one diag.)     ==
C     == NKRY_MAX: Lanczos para.(Krylov subspace dimension)           ==
C     == NKRY_BLOCK: Lanczos para.(Block dimension)                   ==
C     == NPDIP:   Store dipole moments each NPDIP steps               ==
C     == IFTYPE : Interface to classical MD program (EGO)             ==
C     == ICMET  : Interface to classical MD program                   ==
C     ==--------------------------------------------------------------==
      INTEGER      CNTI1,NOMORE,MAXIT,NPS,MDIIS,INSYS,
     *             MGDIIS,NPARA,NCHP,NCHB,NCHS,NCALLS0,NIT0,
     *             KSSTA,NKSSTA,IMOVIE,IPROJ,IVDBRD,IVDBWR,
     *             NTGAUS,INWFUN,NRESTF,NTRANS,NTRAJ,NPRES,NSPLP,
     *             NSORDER,RCOMPB,WCOMPB,NDAVV,N_FRIES,NKRY_MAX,
     *             NKRY_BLOCK,NPDIP,IFTYPE,ICMET,NPREC
      COMMON/CNTI/ CNTI1,NOMORE,MAXIT,NPS,MDIIS,INSYS,
     *             MGDIIS,NPARA,NCHP,NCHB,NCHS,NCALLS0,NIT0,
     *             KSSTA,NKSSTA,IMOVIE,IPROJ,IVDBRD,IVDBWR,
     *             NTGAUS,INWFUN,NRESTF,NTRANS,NTRAJ,NPRES,NSPLP,
     *             NSORDER,RCOMPB,WCOMPB,NDAVV,N_FRIES,NKRY_MAX,
     *             NKRY_BLOCK,NPDIP,IFTYPE,ICMET,NPREC
C     ==================================================================
C     == CNTR1  : Use for the broadcast (dummy variable)              ==
C     == DELT   : Time step                                           ==
C     == EMASS  : Electronic mass (Car-Parrinello)                    ==
C     == EPSOG  : Parameters for the RATTLE step                      ==
C     == EKINW  : Electron dynamics with rescaling of velocities      ==
C     ==          Average kinetic energy                              ==
C     == TOLL   : Tollerance (electron dynamics)                      ==
C     == AMPRE  : Randomize wavefunction parameter                    ==
C     == TEMPW  : Ion dynamics with rescaling of velocities           ==
C     ==          desired temperature                                 ==
C     == TOLP   : Tollerance (ion dynamics)                           ==
C     == ANNERI : Simulated annealing parameter (ions)                ==
C     == ANNERE : Simulated annealing parameter (electrons)           ==
C     == ANNERC : Simulated annealing parameter (cell)                ==
C     == HTHRS  : Hamiltonian cutoff                                  ==
C     == TOLOG  : Convergence orbital tollerance                      ==
C     == TOLNG  : Convergence geometry tollerance                     ==
C     == ECUT   : Cutoff energy                                       ==
C     == SDELTA : Smooth density parameter                            ==
C     == SMF    : ''      ''      ''                                  ==
C     == GCEPS  : Gradient correction cutoff                          ==
C     == WNOSE0 : Characteristic frequency (electrons -- Nose)        ==
C     == WNOSP0 : Characteristic frequency (ions -- Nose)             ==
C     == EPSDAV : Davidson parameter                                  ==
C     == AMPRP  : Randomize ionic positions parameter                 ==
C     == FDIFF  : Finite difference step                              ==
C     == CMASS  : Fictitious MD cell mass                             ==
C     == AMPRC  : Randomize initial cell parameters, amplitude        ==
C     == TEMPC  : Target cell temperature(Kelvin -- Nose)             ==
C     == WNOSC0 : Characteristic frequency (cell -- Nose)             ==
C     == B2LIMIT: Lanczos parameter                                   ==
C     == EKINHR : Cell dynamics with rescaling of velocities          ==
C     ==          Average kinetic energy                              ==
C     == TOLC   : Tollerance (cell dynamics)                          ==
C     == TOLCG  : Convergence cell tollerance                         ==
C     == NEDOF0 : Scaling for elec. DOF (Nose)                        ==
C     ==--------------------------------------------------------------==
      REAL*8       CNTR1,DELT,EMASS,EPSOG,EKINW,TOLL,AMPRE,TEMPW,TOLP,
     *             ANNERI,HTHRS,TOLOG,TOLNG,ECUT,SDELTA,SMF,GCEPS,
     *             WNOSE0,WNOSP0,EPSDAV,AMPRP,FDIFF,CMASS,AMPRC,TEMPC,
     *             WNOSC0,B2LIMIT,EKINHR,TOLC,ANNERC,TOLCG,ANNERE,
     *             NEDOF0
      COMMON/CNTR/ CNTR1,DELT,EMASS,EPSOG,EKINW,TOLL,AMPRE,TEMPW,TOLP,
     *             ANNERI,HTHRS,TOLOG,TOLNG,ECUT,SDELTA,SMF,GCEPS,
     *             WNOSE0,WNOSP0,EPSDAV,AMPRP,FDIFF,CMASS,AMPRC,TEMPC,
     *             WNOSC0,B2LIMIT,EKINHR,TOLC,ANNERC,TOLCG,ANNERE,
     *             NEDOF0
C     ==================================================================
C     == FILEPATH  OPTION (OTHERWISE="./")                            ==
C     ==   FPATH   : PATHNAME OF THE DIRECTORY TO COPY FILES          ==
C     ==   IAPATH  : FIRST NON-BLANK CHARACTER OF FPATH               ==
C     ==   IEPATH  : LAST  NON-BLANK CHARACTER OF FPATH               ==
C     ==--------------------------------------------------------------==
      CHARACTER*80 FPATH
      INTEGER      IAPATH,IEPATH
      COMMON/PATH/ FPATH,IAPATH,IEPATH
C     ==================================================================
C     == DUAL OPTION: SIZE OF THE DENSITY MESH VERSUS                 ==
C     ==              THE WAVEFUNCTION CUTOFF (DEFAULT=4)             ==
C     ==      DUAL=.TRUE. IF USE THIS OPTION                          ==
C     ==--------------------------------------------------------------==
      REAL*8        CDUAL
      LOGICAL       DUAL
      COMMON/DUAL00/CDUAL,DUAL
C     ==================================================================
      INTEGER MAXRF,NSTEPWR,NFNOW,NRCOPY
      PARAMETER (MAXRF=10)
      COMMON/RESTF/NSTEPWR(MAXRF),NFNOW,NRCOPY
C     ==================================================================

