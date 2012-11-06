module mcpathsismc
    use constants, only: dp,int64,sp,sizeof_int
    use util_mod, only: NECI_ICOPY,isnan_neci
   contains
!C.. Calculate RHO^(P)_II without having a stored H matrix
!C.. SAMPLE over distinct nodes, e.g. IJKLI, with paths up to I_HMAX
!C.. generated from these, and summed (e.g IJILKJI), up to max H
!C.. In theory more efficient because RHO_IJ,RHO_JK, etc are calculated
!C.. just once for all these paths.
!C.. I_VMAX is the max number of distinct vertices in a path.
!C.. I_HMAX is the max number of hops in a path.
!C.. NWHTAY contains the number of samples to take for each level
!C.. This is a simple importance sampling method
!C.. NMAX has ARR hidden in it.
      SUBROUTINE  MCPATHSR4(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1, &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,     &
     &               NWHTAY,ILOGGING,ECORE,WLRI,WLSI,DBETA,DLWDB,I_VMIN)
         Use Determinants, only: get_helement, write_det
         USE MCStat
         use CalcData , only : G_VMC_PI,G_VMC_FAC,G_VMC_SEED,CUR_VERT,  &
     &          g_MultiWeight,TMPTHEORY,TVVDISALLOW,TMCDIRECTSUM
         USE Logging , only : G_VMC_LOGCOUNT
         use SystemData, only: BasisFN
         use sym_mod, only: getsym
         use global_utilities
         use mcpathsdata, only: egp
         IMPLICIT NONE
         INTEGER I_VMAX,NEL,NBASIS
         TYPE(MCStats) MCSt
         INTEGER IPATH(NEL,0:I_VMAX)
         INTEGER NI(NEL),I_P,I_HMAX,NMSH,NMAX,NTAY(2),NWHTAY
         INTEGER ILOGGING,I,I_VMIN
         type(timer), save :: proc_timer
         type(timer), save :: proc_timer2
         CHARACTER(20) STR
         real(dp) TOTAL,RHOII(0:I_VMAX)
         HElement_t RHOIJ(0:I_VMAX,0:I_VMAX),RH
         real(dp) ALAT(3),RHOEPS,BETA
         HElement_t UMAT(*)
         complex(dp) FCK(*)
         real(dp) ECORE
         real(dp) FF
         real(dp) WLRI,WLSI
         TYPE(BasisFN) G1(*),KSYM
         INTEGER nBasisMax(5,*),BRR(*)
!CNEL,0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
!C0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER BTABLE(0:I_VMAX)
         LOGICAL TLOG
         INTEGER ICOUNT,ISEED,I_V,L,LT
         real(dp) DBETA
         INTEGER IEXCITS
         type(egp) EXCITGEN(0:I_VMAX) 
         real(dp) DLWDB,DLWDB2,DLWDB3
         HElement_t HIJS(0:I_VMAX)
         INTEGER INODE2(NEL)
         INTEGER, allocatable, target :: NMEM(:)
         INTEGER, target :: NMEMLEN(1)
         real(dp) ODLWDB,OWEIGHT,DLWDBSQ
         real(dp) OPROB,R
         INTEGER ITREE
         INTEGER I_VCUR,I_VM1,I_VM2,I_OVCUR,IOCLS,NTR
         INTEGER ICHANGED
         real(dp) VWEIGHTS(2,0:I_VMAX)
         INTEGER TST
         real(dp) PFAC
         real(dp) WMIN,SUMDLWDB
         real(dp) RAN2
!C.. do some blocking analysis
         INTEGER STORE(6)
         real(dp) DLWDBCORE,WCORE
         INTEGER IOV,IGV,IACC
         LOGICAL TSEQ,TBLOCKING
         real(dp) PREJ,PGR
         REAL(sp) OTIME,NTIME,tarr(2),neci_etime
         integer(int64) LP
         HElement_t :: hel

         ISEED=0  !Init the seed
         OTIME=neci_etime(tarr)
         TST=0
         proc_timer%timer_name='MCPATHSR4 '
         call set_timer(proc_timer)
         TLOG=BTEST(ILOGGING,1)
         TBLOCKING=BTEST(ILOGGING,13)
         IF(TLOG.AND.I_VMIN.EQ.0) THEN
            OPEN(11,FILE="MCPATHS",STATUS="OLD",POSITION='APPEND')
!C.. go to end of file
!            I=FSEEK(11,0,2)
            call write_det (11, NI, .true.)
         ENDIF
         IF(BTEST(ILOGGING,2)) OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,9)) OPEN(12,FILE="VERTEXMC",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,10)) OPEN(22,FILE="VMC",STATUS="UNKNOWN")
!         IF(tMCDirectSum)
!     &      OPEN(122,FILE="MCDIRECTSUM",STATUS="UNKNOWN")
!C.. Set the first node to I_I
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,G1,NBASIS,              &
     &            NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(0)=RH
         RHOIJ(0,0)=RHOII(0)
         WLRI=LOG(RHOII(0))
         hel = get_helement (nI, nI, 0)
         HIJS(0) = hel
         TOTAL=1.0_dp
!C.. These variables have been passed in and represent the precalculated values for 1..I_VMIN-1 vertices
         DLWDBCORE=DLWDB
         WCORE=WLSI
         IF(TMPTHEORY) THEN
            DLWDBCORE=0.0_dp
            WCORE=0.0_dp
            CALL Create(MCSt,I_VMAX,INT(LOG(0.0_dp+NWHTAY)/LOG(2.0_dp)+1),DLWDBCORE,WCore)
            WRITE(6,*) "MCPATHSR4 MP Theory"
         ElseIf(tMCDirectSum) Then
            IF(I_VMIN.NE.0) THEN
                CALL Create(MCSt,I_VMIN,                    &
     &          INT(LOG(0.0_dp+NWHTAY)/LOG(2.0_dp)+1),          &
     &          DLWDBCORE/WCore,WCore)
            ELSE
                CALL Create(MCSt,I_VMAX,                    &
     &          INT(LOG(0.0_dp+NWHTAY)/LOG(2.0_dp)+1),          &
     &          DLWDBCORE/WCore,WCore)
            ENDIF
            WRITE(6,*) "MCPATHSR4 Direct Sum"
            !Print out Eref
!            write(6,*) DLWDBCORE/WCore
            !call neci_flush(6)
         ELSEIF(I_VMIN.NE.0) THEN
            CALL Create(MCSt,I_VMAX,                        &
     &         INT(LOG(NWHTAY/(G_VMC_FAC**4))/LOG(2.0_dp)+1), &
     &      DLWDBCORE/WCORE,WCore)  
            IF(G_VMC_FAC.GT.1.0_dp.OR.G_VMC_FAC.LT.-1.0_dp)     &
     &       STOP "Invalid MULTI MC BIAS"
            WRITE(6,*) "MCPATHSR4 MultiMC.  Bias=",G_VMC_FAC
         ELSE
            DLWDBCORE=HIJS(0)
            CALL Create(MCSt,I_VMAX,                        &
     &         INT(LOG(0.0_dp+NWHTAY)/LOG(2.0_dp)+1),DLWDBCORE,WCore)
         ENDIF
         IF(DBETA.LT.0.0_dp.AND.I_VMIN.EQ.0) THEN
            DLWDB=HIJS(0)
         ENDIF
         IF(TLOG.AND.I_VMIN.EQ.0)                           &
     &         WRITE(11,"(I12,2G25.16,F19.7,2I12,G25.12)")  &
     &         1,TOTAL,TOTAL,0.0_dp,1,1,DLWDB
!C.. we're working in block space, so we work out our current symmetry.
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
         IF(I_HMAX.EQ.-3.OR.I_HMAX.EQ.-4) THEN
            IEXCITS=0
            ISEED=G_VMC_SEED
            CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IEXCITS,   &
     &      ISEED,INODE2)
         ELSE
            TOTAL=0.0_dp
            IF(DBETA.LT.0.0_dp.AND.I_VMIN.EQ.0) THEN
               DLWDB=0.0_dp
            ENDIF
         ENDIF
         I_VM1=2
         I_VM2=I_VMAX
         IF(I_HMAX.EQ.-7.OR.I_HMAX.LE.-12) THEN
!C.. Setup the spin excit generator
            STORE(1)=0
            CALL GENSYMEXCITIT2(NI,NEL,G1,NBASIS,            &
     &         .TRUE.,NMEMLEN,INODE2,I,STORE,3)
            allocate(NMEM(NMEMLEN(1)))
            NMEM(1)=0
            CALL GENSYMEXCITIT2(NI,NEL,G1,NBASIS,            &
     &         .TRUE.,NMEM,INODE2,I,STORE,3)
!C.. Count the excitations (and generate a random one which we throw)
            CALL GENRANDSYMEXCITIT2(NI,NEL,                  &
     &         NMEM,INODE2,ISEED,IEXCITS,PGR)
            IF(I_VMIN.NE.0) THEN
               I_VM1=I_VMIN
               I_VM2=I_VMIN
               I_OVCUR=1
            ELSE
!This needs to be initialized somewhere, and here seems as good a place as any
                 I_OVCUR=0
            ENDIF
         ENDIF
         DO I_V=I_VM1,I_VM2
            CUR_VERT=I_V
            WRITE(STR,"(A,I4)") "FMCPR_",I_VMIN
            proc_timer2%timer_name=STR
            call set_timer(proc_timer2)
            L=0
            LT=0
            BTABLE(0)=0
            ISEED=G_VMC_SEED
            SUMDLWDB=0.0_dp
            IF((I_HMAX.NE.-3.AND.I_HMAX.NE.-4).OR.I_V.LE.IEXCITS+1) THEN
!C.. Do importance sampling with probabilities.  No norm needed
               DLWDB3=0.0_dp
               IOCLS=0
               I_VCUR=1
               NTR=0
               ITREE=1
               DLWDBSQ=0.0_dp
               PFAC=G_VMC_PI
               WMIN=1.0e-8_dp_dp
               IF(I_HMAX.LT.-12.AND.I_HMAX.GT.-19) THEN
!C.. Setup the first graph for memory MC
                  OWEIGHT=1.0_dp
                  ODLWDB=HIJS(0)
!C.. The VWEIGHTS hold the weights of the current vertices.
!C.. These are currently set to (for 1) (rho_jj/rho_ii)**P
!C.. and for 2 (rho_jj/rhoii)**-P
                  VWEIGHTS(1,0)=1.0_dp
                  VWEIGHTS(2,0)=1.0_dp
                  VWEIGHTS(1,1)=VWEIGHTS(1,0)
                  VWEIGHTS(2,1)=VWEIGHTS(2,0)
                  CALL NECI_ICOPY(NEL,IPATH(1,0),1,IPATH(1,I_VCUR),1)
               ELSE
                  OWEIGHT=0.0_dp
                  ODLWDB=0.0_dp
               ENDIF
               IF(I_VMIN.NE.0) THEN
                  ODLWDB=DLWDBCORE
                  OWEIGHT=WCORE
                  OPROB=1   
               ENDIF
               IF(NWHTAY.LT.1) WRITE(6,*) "Warning: MC CYCLE count <0: "  &
     &               ,NWHTAY
               WRITE(6,*) NWHTAY, " MC Cycles"
!               WRITE(6,*) I_V
!               WRITE(6,*) G_VMC_EXCITWEIGHTS(:,CUR_VERT)
!               WRITE(6,*) G_VMC_EXCITWEIGHT(CUR_VERT)
!               CALL neci_flush(6)
               DO ICOUNT=1,NWHTAY
!C                  IF(ICOUNT.EQ.135) TST=1
!C                  DLWDB2=0.0_dp
                  IF(I_HMAX.EQ.-3) THEN
                     FF=FMCPR4B(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,       &
     &                  G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,            &
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,         &
     &                  ECORE,ISEED,DBETA,DLWDB2,HIJS)
                  ELSEIF(I_HMAX.EQ.-4) THEN
                     FF=FMCPR4C(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,       &
     &                  G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,        &
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,         &
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ElseIf(tMCDirectSum) Then
                    oWeight=0.0_dp !Ensure generated graph is "accepted"
                    IF(g_MultiWeight(0).NE.0) THEN
                     R=RAN2(ISEED)*g_MultiWeight(0)
                     I_VCUR=I_VMIN-1
                     DO WHILE (R.GE.0.AND.I_VCUR.LT.I_VMAX)
                        I_VCUR=I_VCUR+1
                        R=R-g_MultiWeight(I_VCUR)
                     ENDDO
                     pFac=g_MultiWeight(I_VCUR)/g_MultiWeight(0)
                    ELSE
                     pFac=1
                     I_VCUR=I_V
                    ENDIF
                  
                  !  write(6,*) ODLWDB,OPROB
                    !Generate an I_V vertex graph:
!      WRITE(43,"(I3,4G25.16)") I_VCUR,ECORE,DLWDB2, ODLWDB,OWEIGHT
                            
                    FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,NBASISMAX,   &
     &                  G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,           &
     &                   RHOEPS,RHOII,RHOIJ,I_HMAX,ILOGGING,              &
     &                   ECORE,ISEED,DBETA,DLWDB2,HIJS,NMEM,              &
     &                   ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,   &
     &                   IACC,0.0_dp,I_VMAX,EXCITGEN)
                    
!     WRITE(43,"(I3,4G25.16)") I_VCUR,ECORE,DLWDB2, ODLWDB,OWEIGHT
!                  CALL neci_flush(43)
                    
                    !Set up values for the AddGraph routine later:
                    ioV=i_vcur
                    igV=i_vcur
!                    i_VCur=i_v
!MPTheory has weights = 1, and does not need to be further divided.
                    IF(.NOT.tMPTheory) THEN
                     FF=oWeight/(oProb) 
                    ENDIF
                    DLWDB2=oWeight*DLWDB2/(oProb)
                  ELSE
!C.. see if we want to increase the number of vertices
                    I_OVCUR=I_VCUR
                    R=RAN2(ISEED)
                    ICHANGED=0
                    IF(I_HMAX.EQ.-13) THEN
                        STOP "I_HMAX=-13 no longer supported"
!                   !
                    ELSEIF(I_HMAX.EQ.-7.OR.I_HMAX.EQ.-19) THEN
!C.. We decide on a graph size.  See 22/8/05 #1,2
!C.. Probability of graph size v is prop to G_VMC_FAC**(v)
                        I_VCUR=0
!                   !   .
                        IF(I_VMIN.EQ.0) THEN
!C.. we're doing MC with lots of different vertex levels
                         IF(G_VMC_FAC.GT.0) THEN
                           DO WHILE(I_VCUR.LT.1.OR.I_VCUR.GT.I_VMAX)
                              R=RAN2(ISEED)
                              I_VCUR=I_VMAX+1+int(LOG(R)/LOG(G_VMC_FAC),sizeof_int)
                           ENDDO
                           PFAC=G_VMC_FAC**(I_VCUR-1)
                         ELSEIF(G_VMC_FAC.EQ.0) THEN
                           R=RAN2(ISEED)
                           I_VCUR=int(R*real(I_VMAX-1,dp)+2.0_dp,sizeof_int)
                           PFAC=1.0_dp/(I_VMAX-1)
                         ELSE
                           DO WHILE(I_VCUR.LT.1.OR.I_VCUR.GT.I_VMAX)
                              R=RAN2(ISEED)
                              I_VCUR=-FLOOR(LOG(R)/LOG(-G_VMC_FAC))
                           ENDDO
                           PFAC=ABS(G_VMC_FAC)**(-I_VCUR)
                         ENDIF
                         IOV=I_OVCUR
                         IGV=I_VCUR
!C.. Markov chain with random generation.
                         FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,       &
     &                     NBASISMAX,G1,                                 &
     &                     NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,          &
     &                   RHOEPS,RHOII,RHOIJ,I_HMAX,ILOGGING,             &
     &                     ECORE,ISEED,DBETA,DLWDB2,HIJS,NMEM,           &
     &                    ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC, &
     &                     IACC,0.0_dp,I_VMAX,EXCITGEN)
!                   !   .
                        ELSE
!..   We're doing MC where we have either the gestalt 1..I_VMIN-1 vertex object
!.. or graphs of size I_VMIN=I_VMAX

                         IF(G_VMC_FAC.LE.0.0_dp) THEN
!.. We're doing non-stochastic time MC but with a renormalized composite gestalt.
!.. G_VMC_FAC is the prob of generating a non-gestalt.
                           IOV=I_OVCUR
                           IF(RAN2(ISEED).GT.ABS(G_VMC_FAC)) THEN
!.. We generate a gestalt
                              I_VCUR=1
                              IGV=1
                              PFAC=1-ABS(G_VMC_FAC)
                              FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,    &
     &                        NBASISMAX,G1,                                &
     &                     NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,            &
     &                   RHOEPS,RHOII,RHOIJ,I_HMAX,ILOGGING,               &
     &                     ECORE,ISEED,DBETA,DLWDBCORE,HIJS,NMEM,          &
     &                    ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,   &
     &                        IACC,WCORE,I_VMAX,EXCITGEN)
                           ELSE
                              R=RAN2(ISEED)*g_MultiWeight(0)
                              I_VCUR=I_VMIN-1
                              DO WHILE (R.GE.0.AND.I_VCUR.LT.I_VMAX)
                                 I_VCUR=I_VCUR+1
                                 R=R-g_MultiWeight(I_VCUR)
                              ENDDO
                              IGV=I_VCUR
                              PFAC=ABS(G_VMC_FAC)*g_MultiWeight(I_VCUR)    &
     &                                       /g_MultiWeight(0)
!C.. Markov chain with random generation.
                              FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,    &
     &                           NBASISMAX,G1,                             &
     &                     NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,            &
     &                   RHOEPS,RHOII,RHOIJ,I_HMAX,ILOGGING,               &
     &                     ECORE,ISEED,DBETA,DLWDB2,HIJS,NMEM,             &
     &                    ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,   &
     &                           IACC,0.0_dp,I_VMAX,EXCITGEN)
                           ENDIF
!                   !   .
                         ELSE
!.. Stochastic time MC
!.. We choose the gestalt with probability G_VMC_FAC
                         TSEQ=.FALSE.
                         IF(I_OVCUR.EQ.1) TSEQ=.TRUE.
                         IF(TSEQ) THEN
                           I_VCUR=1
!.. The last graph was a 1-vertex gestalt.  See how many 1->1 transitions we are likely to generate in a row
                           LP=int(LOG(1.0_dp-RAN2(ISEED))/LOG(1.0_dp-G_VMC_FAC),sizeof_int)
                           IF(LP.GT.0) THEN
                              IGV=1
                              I_OVCUR=1
                              OWEIGHT=WCORE
                              ODLWDB=DLWDBCORE
                              DLWDB2=DLWDBCORE/WCORE
                              OPROB=1.0_dp-ABS(G_VMC_FAC)
                              IOCLS=0
                              FF=1.0_dp
                              ITREE=1
                              CALL AddGraph(MCSt,LP,1,FF,DLWDB2/FF,     &
     &                              WCORE,0,ITREE,1,1,1,                &
     &                              BTEST(ILOGGING,10),OPROB,           &
     &                              TBLOCKING)
!                              WRITE(39,*) 1,1,LP,1
                           ENDIF
                         ENDIF
!                   !
!C                         WRITE(6,*) TSEQ,I_OVCUR,I_VCUR
                         IF(.NOT.TSEQ) THEN
!.. We're at a V-graph. Consider all possible strings of rejected V->0 processes.
!.. First work out the acceptance prob of V->0, pa, and thus the rejectance prob, pr=1-pa, which is high.
!.. The probability of a string of N of these rejectances is (p0 pr)^N, so work out how many we're going to
!.. get in a row
                           IF(TVVDISALLOW) THEN
                              PREJ=((WCORE)/ABS((OWEIGHT)))*OPROB
                              PREJ=1.0_dp-MIN(PREJ,1.0_dp)
                              PGR=PREJ*(1.0_dp) !-G_VMC_FAC)
                           ELSE
                              PREJ=((WCORE)/ABS((OWEIGHT)))*OPROB/(1.0_dp-G_VMC_FAC)
                              PREJ=1.0_dp-MIN(PREJ,1.0_dp)
                              PGR=PREJ*(1.0_dp-G_VMC_FAC)
                           ENDIF
!C.. L is the number of strings of rejections (i.e. of V's we add
                           LP=int(LOG((1.0_dp-RAN2(ISEED)))/LOG(PGR),sizeof_int)
!C                           WRITE(6,*) "V->V",PREJ,L,I_OVCUR,DLWDB2
                           IF(LP.GT.0) THEN
                              CALL AddGraph(MCSt,LP,I_OVCUR,FF,              &
     &                          DLWDB2/FF,OWEIGHT,IOCLS,ITREE,0,I_OVCUR,     &
     &                              1,BTEST(ILOGGING,10),OPROB,              &
     &                          BTEST(ILOGGING,13))
!                              WRITE(39,*) I_OVCUR,1,LP,0
                           ENDIF
!C.. Now if we generate a 1, we accept with prob p1 pa (out of p1 pa + pV=1-(1-pV)*pr), but we need to consider any V's separately, with prob p1
!                   !
!Either we've disallowed V->V' transitions or we generate a 1-v graph
                           IF(TVVDISALLOW.OR.                                &
     &   RAN2(ISEED)*(1-(1-G_VMC_FAC)*PREJ).GE.ABS(G_VMC_FAC)) THEN
                              FF=1.0_dp
                              IOV=I_OVCUR
                              I_OVCUR=1
                              IGV=1
                              DLWDB2=DLWDBCORE/WCORE
                              ODLWDB=DLWDBCORE
                              OWEIGHT=WCORE
                              ITREE=1
!                              LP=1
!  This addgraph occurs later.  We can just set it up
                              IACC=1
!                              CALL AddGraph(MCSt,LP,1,FF,DLWDB2/FF,
!     &                           OWEIGHT,ITREE,1,I_VMIN,IGV)
                              I_VCUR=1
                              IF(TVVDISALLOW) THEN
                                 OPROB=1
                              ELSE
                                 OPROB=1.0_dp-ABS(G_VMC_FAC)
                              ENDIF
                              IOCLS=0
                              TSEQ=.FALSE.
!                              WRITE(39,*) IOV,IGV,1,IACC
                           ELSE
                              TSEQ=.TRUE.
!                              I_OVCUR=I_VMIN
                           ENDIF
                         ENDIF
!                   !
                         IF(TSEQ) THEN
!.. Generate a graph between I_VMIN and I_VMAX
                              R=RAN2(ISEED)*g_MultiWeight(0)
                              I_VCUR=I_VMIN-1
                              DO WHILE (R.GE.0.AND.I_VCUR.LT.I_VMAX)
                                 I_VCUR=I_VCUR+1
                                 R=R-g_MultiWeight(I_VCUR)
                              ENDDO
                              IOV=I_OVCUR
                              IGV=I_VCUR
                              PFAC=ABS(G_VMC_FAC)*g_MultiWeight(I_VCUR)     &
     &                                       /g_MultiWeight(0)
!C                           WRITE(6,*) IOV,I_VCUR,I_OVCUR
                              IF(TVVDISALLOW) OPROB=1
                              FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,     &
     &                         NBASISMAX,G1,NBASIS,                         &
     &                         NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,         &
     &                         RHOII,RHOIJ,I_HMAX,ILOGGING,                 &
     &                         ECORE,ISEED,DBETA,DLWDB2,HIJS,NMEM,          &
     &                         ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,    &
     &                         PFAC,IACC,0.0_dp,I_VMAX,EXCITGEN)
!                              WRITE(39,*) IOV,IGV,1,IACC
!C                              WRITE(6,*) IACC
                         ENDIF
                        ENDIF
                       ENDIF
!                   !
                    ELSEIF(I_HMAX.EQ.-12) THEN
!C.. Markov chain with special generation - does not work
                        STOP "FMCPR4D3 has been removed."
!                   !
                    ENDIF
                  ENDIF
                  LP=1
                  If (FF.eq.0) Then
         CALL AddGraph(MCSt,LP,I_VCUR,(0.0_dp),(0.0_dp)               &
     &               ,OWEIGHT,IOCLS, ITREE,                       &
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),OPROB,       &
     &               TBLOCKING)
                Else

              CALL AddGraph(MCSt,LP,I_VCUR,FF,DLWDB2/FF,OWEIGHT,IOCLS,   &
     &               ITREE,                                              &
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),OPROB,              &
     &                  TBLOCKING)
                  EndIf
                  IF(BTEST(ILOGGING,9).AND.                              &
     &              MOD(ICOUNT,G_VMC_LOGCOUNT).EQ.0) THEN
                 !    CALL WriteStats(MCSt,12)
                     If(tMCDirectSum) Then
                       Write(12, "(I9, 6G25.16)") iCount,             &
     &                 MCSt%wWeightedDelta(0)                         &
     &                 /(MCSt%wWeighting(0)+(iCount+1)*WCore),        &
     &                 MCSt%wWeightedDelta(0),                        &
     &                 MCSt%wWeighting(0), (iCount+1)*WCore,          &
     &      MCSt%wWeightedDeltaSq(0)/iCount                           &
     &             -(MCSt%wWeightedDelta(0)/iCount)**2,               &
     &      MCSt%wWeightingSq(0)/iCount                               &
     &             -(MCSt%wWeighting(0)/iCount)**2
                     Else
                       Write(12, "(I9, 3G25.16)") iCount,             &
     &                  MCSt%wWeightedDelta(0)/MCSt%wWeighting(0),    &
     &                  MCSt%wWeightedDelta(0),MCSt%wWeighting(0)
                     EndIf
                     CALL neci_flush(12)
                  ENDIF
               ENDDO   ! end loop over mc cycles  nwhtay
            ENDIF
            call halt_timer(proc_timer2)
!            F(I_V)=F(I_V)/NWHTAY
!            FSQ(I_V)=FSQ(I_V)/NWHTAY
!            DLWDB=DLWDB+DLWDB3/NWHTAY
!            STD=SQRT(ABS((FSQ(I_V)-F(I_V)*F(I_V))))
!            TOTAL=TOTAL+F(I_V)
!            WRITE(6,*) I_V,F(I_V),TOTAL,get_total_time(proc_timer2),L,LT
!            IF(TLOG)
!     &         WRITE(11,"(I12,2G25.16,G19.7,I12,2G19.7)") I_V,F(I_V),
!     &            TOTAL,get_total_time(proc_timer2),L,STD,DLWDB3/NWHTAY
!C ,F(I_V),FSQ(I_V)
            NTIME=neci_etime(tarr)
            IF(TLOG.AND.(I_HMAX.EQ.-7.OR.I_HMAX.LE.-12)) THEN
              IF(I_VM1.EQ.I_VM2) THEN
!.. just give the additional components for this vertex level
                CALL WriteLongStats2(MCSt,11,WCORE, real(NTIME-OTIME,dp))
              ELSE
!.. rejig the sums so the result is Sum w Delta / Sum w 
                CALL WriteLongStats(MCSt,11,WCORE,DLWDBCORE,         &
     &                                            real(NTIME-OTIME,dp))
              ENDIF
            ENDIF
            CALL neci_flush(11)
!            IF(ISNAN_neci(F(I_V))) THEN
!C.. save all log files
!               ITIME=neci_etime(tarr)
!               CALL neci_flush(11)
!caa            CALL neci_flush(10)
!               CALL LOGNAN(NI,NEL,BETA,ITIME)
!               WRITE(6,*) "WARNING: nan found at time",ITIME
!               WRITE(6,"(A)",advance='no') "  nan det=("
!               DO K=1,NEL
!                  WRITE(6,"(I3,A)",advance='no') NI(K),","
!               ENDDO
!               WRITE(6,"(A)",advance='no') "),"
!            ENDIF
         ENDDO
! Finish off any sequences
         LP=0
         IOV=I_VCUR
         I_VCUR=0
         IGV=0
         IACC=1
         If (FF.eq.0) Then
         CALL AddGraph(MCSt,LP,I_VCUR,(0.0_dp),(0.0_dp)                &
     &               ,OWEIGHT,IOCLS, ITREE,                        &
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),0.0_dp,         &
     &               BTEST(ILOGGING,13))             
         Else
         CALL AddGraph(MCSt,LP,I_VCUR,FF,DLWDB2/FF,OWEIGHT,IOCLS,  &
     &               ITREE,                                        &
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),0.0_dp,         &
     &               BTEST(ILOGGING,13))
         EndIf
         IF(I_HMAX.EQ.-7.OR.I_HMAX.LE.-12) THEN
!C.. setup spin excitation generator
            deallocate(NMEM)
         ENDIF
         IF(TLOG.AND.I_VMIN.EQ.0) CLOSE(11)
         IF(BTEST(ILOGGING,2)) CLOSE(10)
         IF(BTEST(ILOGGING,9)) CLOSE(12)
         IF(BTEST(ILOGGING,10)) CLOSE(22)
         IF(BTEST(ILOGGING,13)) THEN
            Call CalcStDev(MCSt,R)
            WRITE(6,"(A,I4,A,G25.16)") "Error in MC ",I_VM2,":",R
            CALL neci_flush(6)
        ENDIF
         IF(TMPTHEORY) THEN
            Call GetStats(MCSt,0,TOTAL,WLSI,DLWDB2)
            WLSI=LOG(TOTAL)
            DLWDB=DLWDB2
            FF=HIJS(0)
            FF=FF+DLWDB/TOTAL
!  If we're doing a full sum with all levels, then we print data out here.
!  Otherwise it will be done in MCPATHSR10 which calls us
            IF(I_VMIN.EQ.0) THEN
               WRITE(STR,"(A2,I1,A10)") "MP",I_VMAX," ENERGY = "
               WRITE(6,"(A14,G25.16)") STR,FF
            ENDIF
!            DO I=2,I_VMAX
!               FF=FF+MP2E(I)
!               WRITE(STR,"(A2,I1,A10)") "MP",I," ENERGY = "
!               WRITE(6,"(A14,G)") STR,FF
!            ENDDO
         ELSEIF(I_VMIN.EQ.0.OR.tMCDirectSum) THEN
            IF(I_VMIN.NE.0) THEN  !We're just doing this vertex level, not the complete sum
! GetStats(MCS,iV,wAvgWeighting,wAvgWeightedValue,wAvgDelta)
               Call GetStats(MCSt,0,TOTAL,DLWDB2,WLSI)
!               WRITE(6,"(A,I3,4G)") "MCSt:",I_VMIN,TOTAL,WLSI,
!     &            DLWDB2,DLWDB
               TOTAL=TOTAL
               WLSI=TOTAL
!WLSI gets the weight of this level
!DLWDB has the  w E~ sum of past levels, and we add the w E~ of this level to it.
               DLWDB=DLWDB+DLWDB2
            ELSE
               WLSI=1.0_dp !weight of 1v graph
               TOTAL=TOTAL+WLSI
               DLWDB=HIJS(0)
               DO I_V=I_VM1,I_VM2
                  Call GetStats(MCSt,I_V,FF,DLWDB2,WLSI)
                  TOTAL=TOTAL+FF
                  DLWDB=DLWDB+DLWDB2 
               ENDDO
               WLSI=LOG(TOTAL)
               DLWDB=DLWDB/TOTAL   !DLWDBCORE is Etilde of 1v graph
            ENDIF
         ELSE
            Call GetStats(MCSt,0,TOTAL,DLWDB2,WLSI)
            DLWDB=DLWDB2
            IF(TOTAL.LT.0.0_dp) THEN
               TOTAL=-TOTAL
               DLWDB=-DLWDB
            ENDIF
            DLWDB=DLWDB
            WLSI=TOTAL-WCORE
         ENDIF
         call halt_timer(proc_timer)
         CALL Delete(MCSt)
         RETURN
      END SUBROUTINE mcpathsr4


!C.. A function which chooses a random set of I_V connected dets, working out
!C.. loop contribution for that set.
!C.. All nodes are distinct.  Paths IJIKJI etc.
!C.. are generated by permutation from IJKI, and summed up to length I_HMAX
!C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
!C.. This function assumes that there are enough available excitations to 
!C.. form a loop of length I_V.  If not it will probably hang.
      real(dp) FUNCTION FMCPR4(NI,BETA,I_P,IPATH,I_V,NEL,   &
     &   G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,          &
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,       &
     &   ECORE,ISEED)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,I_V,NI(NEL),I_P,IPATH(NEL,0:I_V)
         Type(BasisFN) G1(*)
         INTEGER NBASIS,NMAX
         INTEGER NTAY(2),NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,ALAT(*),UMAT(*),ECORE
         complex(dp) FCK(*)
         real(dp) RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER I,J,K,INODE(NEL)
         INTEGER IADJ(0:I_V-1,0:I_V-1),ICE,IPATH2(0:I_V-1)
         real(dp) RH,CALCPATHS
         INTEGER ICOUNT,IGETEXCITLEVEL,ICMPDETS
         LOGICAL BR
         LOGICAL TLOG,TLOG2
         real(dp) RAN2
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
!C.. Pick random elements of the path and generate excitations
         DO I=1,I_V-1
            K=INT(RAN2(ISEED)*I)
            BR=.TRUE.
            DO WHILE(BR)
               BR=.FALSE.
               CALL GENRANDOMEXCIT(IPATH(1,K),NEL,NBASIS,      &
     &            ABS(NTAY(1)*2),ISEED,INODE)
               DO J=0,I-1
                  IF(ICMPDETS(INODE,IPATH(1,J),NEL).EQ.0) BR=.TRUE.
               ENDDO
            ENDDO
            CALL NECI_ICOPY(NEL,INODE,1,IPATH(1,I),1)
            DO J=0,I
               ICE=IGETEXCITLEVEL(INODE,IPATH(1,J),NEL)
               IF(ICE.LE.2) THEN
                  IADJ(I,J)=1
                  IADJ(J,I)=1
               ELSE
                  IADJ(I,J)=0
                  IADJ(J,I)=0
               ENDIF
               CALL CALCRHO2(INODE,IPATH(1,J),BETA,I_P,NEL,    &
     &            G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,           &
     &            RH,NTAY,ICE,ECORE)
               IF(ABS(RH).GT.RHOEPS) THEN
                  RHOIJ(I,J)=RH
                  RHOIJ(J,I)=RH
               ELSE
                  RHOIJ(I,J)=0.0_dp
                  RHOIJ(J,I)=0.0_dp
               ENDIF
            ENDDO   
            RHOII(I)=RHOIJ(I,I)
         ENDDO
!C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
!C.. matrix elements.  We now call CALCPATHS to get the path weight
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))                        &
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,.TRUE.)
         ENDIF   
         FMCPR4=CALCPATHS(IPATH,RHOII,RHOIJ,I_V,I_HMAX,  &
     &         I_P,NWHTAY,NEL,I_V,ILOGGING)
!C.. Now we need to call the recursive IMCPR4N which will count all the
!C.. ways we could possibly generate this path.
         IPATH2(0)=0
!C.. As we traverse the paths, we set the diagonal elements to 0 to
!C.. indicate we have traversed a particular vertex
         IADJ(0,0)=0
         ICOUNT=IMCPR4N(IADJ,I_V,IPATH2,1)
         IF(TLOG) WRITE(10,*) FMCPR4,ICOUNT 
         FMCPR4=FMCPR4/ICOUNT
         RETURN
      END function fmcpr4

      RECURSIVE INTEGER FUNCTION IMCPR4N(IADJ,I_V,IPATH,IND)      &
     &   RESULT (IMCPR4NRES)
         use sort_mod
         IMPLICIT NONE
         INTEGER I_V,IADJ(0:I_V-1,0:I_V-1),IPATH(0:I_V-1),IND
         INTEGER I,J,INODE,ITOT
!C.. Go through all the nodes currently in the path, and through each
!C.. possible attachment for each node
         IF(IND.EQ.I_V) THEN
            IMCPR4NRES=1
            RETURN
         ENDIF
         ITOT=0
         DO I=0,IND-1
            INODE=IPATH(I)
            DO J=0,I_V-1
!C.. If there's a connection and we haven't been to that node before
               IF(IADJ(INODE,J).NE.0.AND.IADJ(J,J).NE.0) THEN
                  IPATH(IND)=J
                  IADJ(J,J)=0
                  IADJ(INODE,J)=0
                  IADJ(J,INODE)=0
!C.. add that node to the path, and recurse
                  ITOT=ITOT+IMCPR4N(IADJ,I_V,IPATH,IND+1)
                  IADJ(J,J)=1
                  IADJ(INODE,J)=1
                  IADJ(J,INODE)=1
               ENDIF
            ENDDO
         ENDDO
         IMCPR4NRES=ITOT
         RETURN
      END function imcpr4n

      SUBROUTINE GENRANDOMEXCIT(NI,NEL,NBASIS,IEXLEVEL,ISEED,NJ)
         use sort_mod
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,IEXLEVEL,ISEED,NJ(NEL)
         INTEGER I,J,K,IEX,IEXL2
         LOGICAL BR
         real(dp) RAN2
         IF(IEXLEVEL.GT.2)                                     &
     &    STOP "Cannot handle more than double excitations."
         IF(IEXLEVEL.LT.2) IEXL2=IEXLEVEL
         IF(IEXLEVEL.EQ.2) THEN
            IEX=int(RAN2(ISEED)*real((NBASIS-NEL)*NEL)*                     &
     &         (1.0_dp+real((NBASIS-NEL-1)*(NEL-1),dp)/4.0_dp),sizeof_int)
            IF(IEX.LT.(NBASIS-NEL)*NEL) THEN
               IEXL2=1
            ELSE
               IEXL2=2
            ENDIF
!C            WRITE(41,"(I)",advance='no') IEXL2
         ENDIF
         CALL NECI_ICOPY(NEL,NI,1,NJ,1)
         DO IEX=1,IEXL2
            BR=.TRUE.
!C.. Find an electron we haven't excited before
            DO WHILE (BR)
               I=INT(RAN2(ISEED)*NEL)+1
               IF(NJ(I).EQ.NI(I)) BR=.FALSE.
            ENDDO
            BR=.TRUE.
            DO WHILE (BR)
               BR=.FALSE.
               J=INT(RAN2(ISEED)*NBASIS)+1
               DO K=1,NEL
!C.. If the new basis fn's in our original det or our new one, we loop again
                  IF(NI(K).EQ.J.OR.NJ(K).EQ.J) BR=.TRUE.
               ENDDO
            ENDDO
            NJ(I)=J
         ENDDO
         call sort (nJ)
!C         WRITE(41,*) IGETEXCITLEVEL(NI,NJ,NEL),NI(1),NI(2),NJ(1),NJ(2)
         RETURN
      END subroutine genrandomexcit

!C.. A function which chooses a random set of I_V connected dets, working out
!C.. loop contribution for that set, and weighting with the appropriate
!C.. probability - not importance sampling with the RHOs
!C.. All nodes are distinct.  Paths IJIKJI etc.
!C.. are generated by permutation from IJKI, and summed up to length I_HMAX
!C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
!C.. This function assumes that there are enough available excitations to 
!C.. form a loop of length I_V.  If not it will probably hang.
!C..
!C.. All singles and doubles are allowed (even if RHO_IJ=0) so this is
!C.. not very efficient.
      real(dp) FUNCTION FMCPR4B(NI,BETA,I_P,IPATH,I_V,NEL,             &
     &   NBASISMAX,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,           &
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,                  &
     &   ECORE,ISEED,DBETA,DLWDB,HIJS)
         use SystemData , only : BasisFN
         use Determinants, only: get_helement
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_V,I_P,IPATH(NEL,0:I_V)
         INTEGER nBasisMax(5,*),NBASIS,NMAX
         Type(BasisFn) G1(*)
         INTEGER NTAY(2),NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,ALAT(*),ECORE
         HElement_t UMAT(*)
         complex(dp) FCK(*)
         real(dp) RHOEPS,RHOII(0:I_V)
         HElement_t RHOIJ(0:I_V,0:I_V)
         INTEGER NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT
         INTEGER I,ICE,IC
         real(dp) XIJ(0:I_V-1,0:I_V-1)         
         real(dp) RH,CALCPATHS_N
         LOGICAL LISINPATH
         INTEGER IGETEXCITLEVEL
         LOGICAL TLOG,TLOG2,TLOG3
         real(dp) RP,DBETA,DLWDB
         HElement_t HIJS(0:I_V)
         LOGICAL ISUHFDET
         INTEGER ICLS
         HElement_t :: hel
         RP=0
         I_VNEXT=1
         CALL NECI_ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0
!C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
!C.. doesn't bother calculating the RHO_JJ.
!C         CALL GENSYMDETSSDN(NI,KSYM,NEL,G1,BRR,NBASIS,0,
!C     &         IC,NBASISMAX,
!C     &      BETA,I_P,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,X,RP,0.0_dp,
!C     &      RHOEPS)
         CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IC,         &
     &      ISEED,INODE2)
!C.. Diagonal X elements contain the normalization of that node
!C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         NLIST=1
         ICOUNT=0
         DO WHILE(I_VNEXT.LT.I_V)
!C.. pick a random excitation of where we are (INODE) of the appropriate
!C.. symmetry.  We don't worry about weighting this at the moment
!C            WRITE(57,"(A1)",advance='no') "N"
            CALL GENRANDOMSPINEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,IC,   &
     &      ISEED,INODE2)
!C            CALL WRITEDET(35,INODE2,NEL,.FALSE.)
            CALL NECI_ICOPY(NEL,INODE2,1,INODE,1)
!C            CALL GENRANDOMEXCITSYM(INODE,NEL,NBASIS,2,ISEED,
!C     &         INODE2)
!C.. If the new node is not in the path, add it.
            IF(.NOT.LISINPATH(INODE2,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL NECI_ICOPY(NEL,INODE2,1,IPATH(1,I_VNEXT),1)
!C.. Count the number of connections from this node (RP=0)
!C               NLIST=0
!C               CALL GENSYMDETSSDN(INODE2,KSYM,NEL,G1,BRR,NBASIS,LSTE,
!C     &            NLIST,NBASISMAX,BETA,I_P,NMSH,FCK,NMAX,
!C     &            ALAT,UMAT,NTAY,ECORE,X,RP,0.0_dp,RHOEPS)
               CALL GENRANDOMSPINEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,IC,      &
     &            ISEED,INODE2)

               XIJ(I_VNEXT,I_VNEXT)=IC
               NLIST=IC
!CNLIST
!C.. Update the rho and X (probability) matrices with this new node
               CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,                        &
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,                         &
     &             RH,NTAY,0,ECORE)
               RHOII(I_VNEXT)=RH
               RHOIJ(I_VNEXT,I_VNEXT)=RH
               DO I=0,I_VNEXT-1
                  ICE=IGETEXCITLEVEL(INODE,IPATH(1,I),NEL)
                  IF(ICE.LE.2) THEN
                     XIJ(I_VNEXT,I)=1.0_dp/NLIST
                     IF(ISUHFDET(IPATH(1,I),NEL)) THEN
                        IF(ICE.EQ.2) THEN
                           XIJ(I,I_VNEXT)=1.0_dp/XIJ(I,I)
                        ELSE
                           XIJ(I,I_VNEXT)=0.0_dp
                        ENDIF
                     ELSE
                        XIJ(I,I_VNEXT)=1.0_dp/XIJ(I,I)
                     ENDIF   
                  ELSE
                     XIJ(I_VNEXT,I)=0.0_dp
                     XIJ(I,I_VNEXT)=0.0_dp
                  ENDIF
                  IF(I.EQ.0) THEN
                     hel = get_helement (nI, iNode)
                     HIJS(I_vNext) = hel
                  ENDIF
                  CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,                &
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,                         &
     &             RH,NTAY,ICE,ECORE)
                  IF(ABS(RH).GT.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.0_dp
                     RHOIJ(I_VNEXT,I)=0.0_dp
                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.500) THEN
!C.. give up trying to find something to attach, and go home
               FMCPR4B=0.0_dp
               RETURN
            ENDIF
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
!C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
!C.. matrix elements.
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))                                    &
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,.TRUE.)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
!C.. GETPATHPROB gives us the probability of generating the path
         RH=GETPATHPROB(XIJ,I_V)
!C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         FMCPR4B=CALCPATHS_N(RHOII,RHOIJ,I_V,I_HMAX,                       &
     &         I_P,RH*NWHTAY,DBETA,DLWDB,HIJS,ICLS)
         IF(TLOG) WRITE(10,"(3E25.16, I7)") FMCPR4B,RH,DLWDB,ICLS
         IF(RH.GT.0.0_dp) THEN
!C.. Unbias the sum 
            FMCPR4B=FMCPR4B/RH
            DLWDB=DLWDB/RH
         ELSE
            FMCPR4B=0.0_dp
            DLWDB=0.0_dp
         ENDIF
!C         WRITE(35,*)
         RETURN
      END function fmcpr4b

!C.. A function which chooses a random set of I_V connected dets, working out
!C.. loop contribution for that set, and weighting with the appropriate
!C.. probability - not importance sampling with the RHOs
!C.. All nodes are distinct.  Paths IJIKJI etc.
!C.. are generated by permutation from IJKI, and summed up to length I_HMAX
!C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
!C.. This function assumes that there are enough available excitations to 
!C.. form a loop of length I_V.  If not it will probably hang.
!C.. This chooses vertices with weights according to the RHO_IJ values
!C.. and RHO_II values.
!C.. (26/09/04)
!C.. If a RHOIJ is non-zero, then its weight XIJ will be non-zero.
!C.. XIJ is proportional to RHO_JJ**ABS(RP). 
      real(dp) FUNCTION FMCPR4C(NI,BETA,I_P,IPATH,I_V,NEL,                   &
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,             &
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,                        &
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS)
         use SystemData , only : BasisFN
         use Determinants, only: get_helement
         IMPLICIT NONE
         TYPE(BasisFN) :: KSYM,G1(*)
         INTEGER NEL,I_V,NI(NEL),I_P,IPATH(NEL,0:I_V)
         INTEGER nBasisMax(5,*),NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY(2),NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,ALAT(*),ECORE
         HElement_t UMAT(*)
         complex(dp) FCK(*)
         real(dp) RHOEPS,RHOII(0:I_V)
         HElement_t RHOIJ(0:I_V,0:I_V)
         INTEGER LSTE(NEL),NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT,ICURNODE
         INTEGER I,ICE,IC,IONODE
         real(dp) XIJ(0:I_V-1,0:I_V-1)         
         real(dp) RH,CALCPATHS_N,X
         INTEGER IGETEXCITLEVEL,IISINPATH
         LOGICAL TLOG,TLOG2,TLOG3
         real(dp) RP,PP
         real(dp) DBETA,DLWDB
         HElement_t HIJS(0:I_V)
         INTEGER ICLS
         real(dp) RAN2
         HElement_t :: hel
!C.. we hard code RP as P/50, although this should be an empirical
!C.. parameter.
         RP=-I_P/50.0_dp
!C..I_P
         I_VNEXT=1
         CALL NECI_ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0
!C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
!C.. doesn't bother calculating the RHO_JJ.
         CALL GENSYMDETSSDN(NI,KSYM,NEL,G1,BRR,NBASIS,(/0/),             &
     &         IC,NBASISMAX,BETA,I_P,NMSH,FCK,NMAX,ALAT,             &
     &         UMAT,NTAY,ECORE,X,RP,0.0_dp,RHOEPS)
!C.. Diagonal X elements contain the normalization of that node
!C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=X
         IF(X.EQ.0.0_dp) THEN
            FMCPR4C=0.0_dp
            RETURN
         ENDIF       
!C
         NLIST=1
         ICOUNT=0
         ICURNODE=0
         DO WHILE(I_VNEXT.LT.I_V)
!C.. XIJ(ICURNODE,ICURNDOE) is the norm const for the current node
            PP=RAN2(ISEED)*XIJ(ICURNODE,ICURNODE)
            IC=1
            CALL GENSYMDETSSDN(IPATH(1,ICURNODE),KSYM,NEL,G1,BRR,NBASIS,      &
     &         INODE2,IC,NBASISMAX,BETA,I_P,NMSH,FCK,NMAX,ALAT,               &
     &         UMAT,NTAY,ECORE,X,RP,PP,RHOEPS)
            CALL NECI_ICOPY(NEL,INODE2,1,INODE,1)
!C            CALL GENRANDOMEXCITSYM(INODE,NEL,NBASIS,2,ISEED,
!C     &         INODE2)
!C.. If the new node is not in the path, add it.
            IONODE=ICURNODE
            ICURNODE=IISINPATH(INODE2,IPATH,NEL,I_VNEXT,-1) 
            IF(ICURNODE.EQ.-1) THEN
               ICURNODE=I_VNEXT
               CALL NECI_ICOPY(NEL,INODE2,1,IPATH(1,I_VNEXT),1)
!C.. Count the number of connections from this node (RP=0)
               NLIST=0
               CALL GENSYMDETSSDN(INODE,KSYM,NEL,G1,BRR,NBASIS,LSTE,          &
     &            NLIST,NBASISMAX,BETA,I_P,NMSH,FCK,NMAX,                     &
     &            ALAT,UMAT,NTAY,ECORE,X,RP,0.0_dp,RHOEPS)
               XIJ(I_VNEXT,I_VNEXT)=X
!C.. Update the rho and X (probability) matrices with this new node
               CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,                        &
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,                         &
     &             RH,NTAY,0,ECORE)
               RHOII(I_VNEXT)=RH
               RHOIJ(I_VNEXT,I_VNEXT)=RH
               DO I=0,I_VNEXT-1
                  ICE=IGETEXCITLEVEL(INODE,IPATH(1,I),NEL)
                  CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,                &
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,                         &
     &             RH,NTAY,ICE,ECORE)
                  IF(ABS(RH).GT.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
!C                     XIJ(I_VNEXT,I)=1.0_dp/XIJ(I_VNEXT,I_VNEXT)
!C                     XIJ(I,I_VNEXT)=1.0_dp/XIJ(I,I)
                     XIJ(I_VNEXT,I)=                                          &
     &                  RHOII(I)**ABS(RP)/XIJ(I_VNEXT,I_VNEXT)
                     XIJ(I,I_VNEXT)=                                          &
     &                  RHOII(I_VNEXT)**ABS(RP)/XIJ(I,I)
                     if(I == 0) then
                        hel = get_helement(nI, iNode, ICE)
                        HIJS(I_vNext) = hel
                     endif
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.0_dp
                     RHOIJ(I_VNEXT,I)=0.0_dp
                     XIJ(I_VNEXT,I)=0.0_dp
                     XIJ(I,I_VNEXT)=0.0_dp
                     IF(I.EQ.0) HIJS(I_VNEXT)=0.0_dp
                  ENDIF
!C                  IF(ICE.LE.2) THEN
!C                     XIJ(I_VNEXT,I)=RHOII(I)**RP/X
!C                     XIJ(I,I_VNEXT)=RHOII(I_VNEXT)**RP/XIJ(I,I)
!C                     XIJ(I_VNEXT,I)=RHOIJ(I,I_VNEXT)**RP/X
!C                     XIJ(I,I_VNEXT)=RHOIJ(I,I_VNEXT)**RP/XIJ(I,I)
!C                     XIJ(I_VNEXT,I)=RHOIJ(I,I_VNEXT)**RP/X
!C                     XIJ(I,I_VNEXT)=RHOIJ(I,I_VNEXT)**RP/XIJ(I,I)
!C                  ELSE
!C                     XIJ(I_VNEXT,I)=0.0_dp
!C                     XIJ(I,I_VNEXT)=0.0_dp
!C                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.500) THEN
!C.. give up trying to find something to attach, and go home
               FMCPR4C=0.0_dp
               RETURN
            ENDIF
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
!C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
!C.. matrix elements.
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))                                          &
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,.TRUE.)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
!C.. GETPATHPROB gives us the probability of generating the path
         RH=GETPATHPROB(XIJ,I_V)
!C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         FMCPR4C=CALCPATHS_N(RHOII,RHOIJ,I_V,I_HMAX,                       &
     &         I_P,RH*NWHTAY,DBETA,DLWDB,HIJS,ICLS)
         IF(TLOG) WRITE(10,"(3E25.16, I7)") FMCPR4C,RH,DLWDB,ICLS
         IF(RH.GT.0.0_dp) THEN
!C.. Unbias the sum 
            FMCPR4C=FMCPR4C/RH
            DLWDB=DLWDB/RH
         ELSE
            FMCPR4C=0.0_dp
         ENDIF
         RETURN
      END function fmcpr4c


!C.. Generate a path probability for the modified path generation
!C.. algorithm in FMCPR4D

      real(dp) FUNCTION GETPATHPROB2(XIJ,I_V)
         use CalcData , only : G_VMC_PI
         IMPLICIT NONE
         INTEGER I_V
         real(dp) XIJ(I_V,I_V)
         INTEGER IPATH(I_V)
         real(dp) M(I_V,I_V)
         real(dp) RET
         real(dp) PI,PJ
         real(dp) INV(I_V,I_V)
         IF(I_V.EQ.2) THEN
            GETPATHPROB2=XIJ(1,2)
         ELSEIF(I_V.EQ.3) THEN
            IF(G_VMC_PI.EQ.-1.0_dp) THEN
                PI=0.5
                PJ=0.5
            ELSE
                PI=G_VMC_PI
                PJ=1.0_dp-PI
            ENDIF
            GETPATHPROB2= XIJ(1,2)*(PI*XIJ(1,3)+PJ*XIJ(2,3))/     &
     &                     (1-PI*XIJ(1,2)-PJ*XIJ(2,1))            &
     &                  +XIJ(1,3)*(PI*XIJ(1,2)+PJ*XIJ(3,2))/      &
     &                     (1-PI*XIJ(1,3)-PJ*XIJ(3,1))
!C            IPATH(1)=1
!C            RET=0.0_dp
!C            CALL GETPP2_R(IPATH,XIJ,M,I_V,2,RET,1.0_dp,INV)
!C            WRITE(61,*) GETPATHPROB2,RET
!C         ENDIF
         ELSE
!C.. We require a recursive function summing over all possible
!C.. permutations of jkl... and calculating appropriate inverse matrix
!C.. elements
            IPATH(1)=1
            M(1,1)=1.0_dp
            RET=0.0_dp
!C            STOP 'Cannot handle new path gen with IV_MAX>3'
            CALL GETPP2_R(IPATH,XIJ,M,I_V,2,RET,1.0_dp,INV)
            GETPATHPROB2=RET
         ENDIF
         RETURN
      END function getpathprob2
   
      real(dp) FUNCTION GETPATHPROB(XIJ,I_V)
         IMPLICIT NONE
         INTEGER I_V
         real(dp) XIJ(I_V,I_V)
         INTEGER IPATH(I_V)
         real(dp) M(I_V,I_V)
         real(dp) RET
         IF(I_V.EQ.2) THEN
            GETPATHPROB=XIJ(1,2)
         ELSEIF(I_V.EQ.3) THEN
            GETPATHPROB= XIJ(1,2)*(XIJ(2,3)+XIJ(2,1)*XIJ(1,3))/      &
     &                     (1-XIJ(1,2)*XIJ(2,1))                     &
     &                  +XIJ(1,3)*(XIJ(3,2)+XIJ(3,1)*XIJ(1,2))/      &
     &                     (1-XIJ(1,3)*XIJ(3,1))
!C         ENDIF
         ELSE
!C.. We require a recursive function summing over all possible
!C.. permutations of jkl... and calculating appropriate inverse matrix
!C.. elements
            IPATH(1)=1
            M(1,1)=1.0_dp
            RET=0.0_dp
            CALL GETPP_R(IPATH,XIJ,M,I_V,2,RET,1.0_dp)
            GETPATHPROB=RET
         ENDIF
         RETURN
      END function getpathprob

      RECURSIVE SUBROUTINE GETPP2_R(IPATH,XIJ,M,I_VMAX,I_V,RET,TMP,INV)
         use CalcData , only : G_VMC_PI
         IMPLICIT NONE
         INTEGER I_VMAX,I_V,IPATH(I_VMAX)
         real(dp) XIJ(I_VMAX,I_VMAX),M(I_VMAX,I_VMAX)
         real(dp) INV(I_VMAX,I_VMAX)
         real(dp) WORK(I_VMAX,I_VMAX),RET,TMP
         INTEGER I,J,K,IPIVOT(I_VMAX),INFO,L
         real(dp) P(I_VMAX)
         LOGICAL T
!C.. M contains the I-X matrix being constructed
!C.. Add this node
         IF(I_V.GT.I_VMAX) THEN
            RET=RET+TMP
            RETURN
         ENDIF
         DO J=2,I_VMAX
            T=.TRUE.
            DO K=2,I_V-1
               IF(IPATH(K).EQ.J) T=.FALSE.
            ENDDO
!C.. If we're not already in the path
            IF(T) THEN
               IPATH(I_V)=J
!C.. Create the prob matrix
               IF(G_VMC_PI.eq.-1.0_dp) THEN
                   DO I=1,I_V-1
                    P(I)=1.0_dp/(I_V-1.0_dp)
                   ENDDO
               ELSE
                IF(I_V.GT.2) THEN
                  P(1)=G_VMC_PI
                ELSE
                  P(1)=1.0_dp
                ENDIF
                DO I=2,I_V-1
                     P(I)=(1.0_dp-G_VMC_PI)/(I_V-2.0_dp)
                ENDDO
               ENDIF
               P(I_V)=0.0_dp
!C.. Create matrix M 
               DO I=1,I_V
                  DO K=1,I_V
                     INV(I,K)=0.0_dp
                     DO L=1,I_V
                        IF(L.NE.K) INV(I,K)=INV(I,K)              &
     &                     -P(L)*XIJ(IPATH(L),IPATH(K))
                     ENDDO
                     IF(I.EQ.I_V) INV(I,K)=0.0_dp
                     IF(I.EQ.K) THEN
                        INV(I,K)=1.0_dp+INV(I,K)
                     ENDIF
                  ENDDO
!C                  WRITE(61,*) (INV(I,K),K=1,I_V)
               ENDDO
!C.. Invert the matrix, and pass the multiplicative factor on.
!C.. LU decomp
               CALL DGETRF(I_V,I_V,INV,I_VMAX,IPIVOT,INFO)
!C.. Inverse after LU decomp
               CALL DGETRI(I_V,INV,I_VMAX,IPIVOT,WORK,I_VMAX,INFO)
               CALL GETPP2_R(IPATH,XIJ,M,I_VMAX,I_V+1,RET,        &
     &            TMP*INV(I_V-1,I_V),INV)
            ENDIF
         ENDDO
         RETURN
      END subroutine getpp2_r

      RECURSIVE SUBROUTINE GETPP_R(IPATH,XIJ,M,I_VMAX,I_V,RET,TMP)
         IMPLICIT NONE
         INTEGER I_VMAX,I_V,IPATH(I_VMAX)
         real(dp) XIJ(I_VMAX,I_VMAX),M(I_VMAX,I_VMAX)
         real(dp) INV(I_VMAX,I_VMAX)
         real(dp) WORK(I_VMAX,I_VMAX),RET,TMP
         INTEGER I,J,K,IPIVOT(I_VMAX),INFO
         LOGICAL T
!C.. M contains the I-X matrix being constructed
!C.. Add this node
         IF(I_V.GT.I_VMAX) THEN
            RET=RET+TMP
            RETURN
         ENDIF
         DO J=2,I_VMAX
            T=.TRUE.
            DO K=2,I_V-1
               IF(IPATH(K).EQ.J) T=.FALSE.
            ENDDO
!C.. If we're not already in the path
            IF(T) THEN
               IPATH(I_V)=J
               DO I=1,I_V-1
                  M(I,I_V)=-XIJ(IPATH(I),J)
                  M(I_V,I)=-XIJ(J,IPATH(I))
               ENDDO
               M(I_V,I_V)=1.0_dp
               CALL DCOPY(I_VMAX**2,M,1,INV,1)
               DO I=1,I_V-1
                  INV(I_V,I)=0.0_dp
               ENDDO
!C.. Invert the matrix, and pass the multiplicative factor on.
               CALL DGETRF(I_V,I_V,INV,I_VMAX,IPIVOT,INFO)
               CALL DGETRI(I_V,INV,I_VMAX,IPIVOT,WORK,I_VMAX,INFO)
               CALL GETPP_R(IPATH,XIJ,M,I_VMAX,I_V+1,RET,               &
     &            TMP*INV(I_V-1,I_V))
            ENDIF
         ENDDO
         RETURN
      END subroutine getpp_r

      SUBROUTINE GENRANDOMEXCITSYM(NI,NEL,NBASIS,IEXLEVEL,ISEED,        &
     &            NJ)
         use SystemData, only: BasisFN
         use sort_mod
         use sym_mod, only: lChkSym
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,IEXLEVEL,ISEED,NJ(NEL)
         INTEGER I,J,K,IEX,IEXL2
         LOGICAL BR,BR2
         real(dp) RAN2
         IF(IEXLEVEL.GT.2)                                              &
     &    STOP "Cannot handle more than double excitations."
         IF(IEXLEVEL.LT.2) THEN
            STOP "No sym excitations for IEXLEVEL<2"
         ENDIF
         IEXL2=IEXLEVEL
         BR2=.TRUE.
         DO WHILE(BR2)
            CALL NECI_ICOPY(NEL,NI,1,NJ,1)
            DO IEX=1,IEXL2
               BR=.TRUE.
!C.. Find an electron we haven't excited before
               DO WHILE (BR)
                  I=INT(RAN2(ISEED)*NEL)+1
                  IF(NJ(I).EQ.NI(I)) BR=.FALSE.
               ENDDO
               BR=.TRUE.
               DO WHILE (BR)
                  BR=.FALSE.
                  J=INT(RAN2(ISEED)*NBASIS)+1
                  DO K=1,NEL
!C.. If the new basis fn's in our original det or our new one, we loop again
                     IF(NI(K).EQ.J.OR.NJ(K).EQ.J) BR=.TRUE.
                  ENDDO
               ENDDO
               NJ(I)=J
            ENDDO
            call sort (nJ)
            call stop_all('GENRANDOMEXCITSYM',                             &
     &                            'BROKEN CODE! Fix LCKSYM call.')
            !BR=LCHKSYM(NI,NJ,NEL,G1,NBASISMAX)
            BR2=.NOT.BR
         ENDDO
         RETURN
      END subroutine genrandomexcitsym


      SUBROUTINE GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,NEXCITS,       &
     &      ISEED,NJ)
         use SystemData, only: BasisFN
         use sort_mod
         use sym_mod, only: getsym
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,NEXCITS,ISEED,NJ(NEL)
         INTEGER NB2,NA,NB,NEX1,NEX2
         INTEGER IEX1,IEX2,IIEX1,IIEX2
         INTEGER nBasisMax(5,*)
         INTEGER NEX1A,NEX1B,NEX2A,NEX2B
         type(BasisFN) G1(*)
         INTEGER NEXAB
         real(dp) R,IEX
         LOGICAL ISUHFDET
         real(dp) RAN2
         TYPE(BasisFN) KSym
!C         CALL WRITEDET(56,NI,NEL,.TRUE.)
!C         CALL neci_flush(56)
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
            CALL NECI_ICOPY(NEL,NI,1,NJ,1)
         NB2=NBASIS/2
         NB=(NEL-KSYM%Ms)/2
         NA=NEL-NB
         NEX1A=NA*(NB2-NA)
         NEX1B=NB*(NB2-NB)
         IF(ISUHFDET(NI,NEL)) THEN
            NEX1A=0
            NEX1B=0
         ENDIF
         NEX2A=NA*(NB2-NA)*(NA-1)*(NB2-NA-1)/4
         NEX2B=NB*(NB2-NB)*(NB-1)*(NB2-NB-1)/4
         NEXAB=NA*NB*(NB2-NA)*(NB2-NB)
         NEXCITS=NEX1A+NEX1B+NEX2A+NEX2B+NEXAB
!C         WRITE(56,*) NEX1A,NEX1B,NEX2A,NEX2B,NEXAB
         NEX1B=NEX1B+NEX1A
         NEX2A=NEX2A+NEX1B
         NEX2B=NEX2B+NEX2A
         NEXAB=NEXAB+NEX2B
         NEX2=0
         R=RAN2(ISEED)
         IEX=R*NEXCITS
         IF(IEX.LT.NEX1A) THEN
!C            WRITE(57,"(A2)",advance='no') "A1"
            NEX1=1
         ELSEIF(IEX.LT.NEX1B) THEN
!C            WRITE(57,"(A2)",advance='no') "B1"
            NEX1=-1
         ELSEIF(IEX.LT.NEX2A) THEN
!C            WRITE(57,"(A2)",advance='no') "A2"
            NEX1=1
            NEX2=1
         ELSEIF(IEX.LT.NEX2B) THEN
!C            WRITE(57,"(A2)",advance='no') "B2"
            NEX1=-1
            NEX2=-1
         ELSE
!C            WRITE(57,"(A2)",advance='no') "AB"
            NEX1=1
            NEX2=-1
         ENDIF
         IEX1=0
         CALL FINDELECSPIN(NI,NEL,NEX1,G1,ISEED,IEX1)
         IF(NEX2.NE.0) THEN
            IEX2=IEX1
            CALL FINDELECSPIN(NI,NEL,NEX2,G1,ISEED,IEX2)
         ELSE
            IEX2=0
         ENDIF
         IIEX1=0
         IIEX2=0
         CALL FINDNEWELECSPIN(NI,NEL,NEX1,G1,ISEED,NBASIS,IIEX1,IIEX1)
         IF(IEX2.NE.0) THEN
            CALL FINDNEWELECSPIN(NI,NEL,NEX2,G1,ISEED,NBASIS,IIEX2,           &
     &         IIEX1)
            NJ(IEX2)=IIEX2
         ENDIF
         NJ(IEX1)=IIEX1 
         call sort (nJ)
!C         CALL WRITEDET(57,NJ,NEL,.TRUE.)
!C         CALL neci_flush(57)
         RETURN
      END subroutine genrandomspinexcit
      
      SUBROUTINE FINDELECSPIN(NI,NEL,NSPIN,G1,ISEED,IEX)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NSPIN,IEX,ISEED,IEL
         TYPE(BasisFN) G1(*)
         LOGICAL BR
         real(dp) RAN2
         BR=.TRUE.
         DO WHILE(BR)
            IEL=int(RAN2(ISEED)*real(NEL+1,dp),sizeof_int)
            IF(G1(NI(IEL))%Ms.EQ.NSPIN.AND.IEL.NE.IEX) BR=.FALSE.
         ENDDO
         IEX=IEL
         RETURN
      END subroutine findelecspin
      SUBROUTINE FINDNEWELECSPIN(NI,NEL,NSPIN,G1,ISEED,NBASIS,IEX,IEX1)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NSPIN,IEX,ISEED,IEL,IEX1,NBASIS
         TYPE(BasisFN) G1(*)
         LOGICAL BR
         INTEGER I
         real(dp) RAN2
         BR=.TRUE.
         DO WHILE(BR)
            IEL=int(RAN2(ISEED)*real(NBASIS+1,dp),sizeof_int)
            IF(G1(IEL)%Ms.EQ.NSPIN) THEN
               BR=.FALSE.
               IF(IEX1.EQ.IEL) BR=.TRUE.
               DO I=1,NEL
                  IF(IEL.EQ.NI(I)) BR=.TRUE.
               ENDDO
            ENDIF
         ENDDO
         IEX=IEL
         RETURN
      END subroutine findnewelecspin 


!C.. Like MCPATHSR4C, this selects a n-vertex graph.  it chooses which
!Cvertex to excite from randomly however.
!C.. 
!C.. This does a markov chain monte carlo also - the probabilities are 
!C.. weighted with the values of S' (the double-counting corrected
!C.. weight.
!.. NMAX has ARR hidden in it
      FUNCTION FMCPR4D2(NI,BETA,I_P,IPATH,I_V,NEL,                      &
     &   NBASISMAX,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,              &
     &   RHOEPS,RHOII,RHOIJ,I_HMAX,ILOGGING,                            &
     &   ECORE,ISEED,DBETA,DLWDB,HIJS,NMEM,OETILDE,OPROB,               &
     &   I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,IACC,INWI,I_VMAX,EXCITGEN)
         USE SystemData , only : BasisFN
         use CalcData, only: tMPTheory, tMCDirectSum, pGenEpsilon
         use CalcData , only : TMPTHEORY,TMCDIRECTSUM,PGenEpsilon
         use mcpathsdata, only: egp
         IMPLICIT NONE
         real(dp) FMCPR4D2
         INTEGER NEL,NI(NEL),I_V,I_P,IPATH(NEL,0:I_V),I_VMAX
         INTEGER nBasisMax(5,*),NBASIS,NMAX
         Type(BasisFn) G1(*)
         INTEGER NTAY(2),I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,ALAT(*),ECORE
         complex(dp) FCK(*)
         INTEGER I_OVCUR,IOCLS,ITREE
         HElement_t UMat(*) 
         real(dp) RHOEPS
         real(dp) RHOII(0:I_V),INWI
         HElement_t RHOIJ(0:I_V,0:I_V)
         INTEGER I,IC
         real(dp) XIJ(0:I_V-1,0:I_V-1)         
         real(dp) CALCPATHS_N
         real(dp) RHX
         LOGICAL TLOG,TLOG2,TLOG3
         real(dp) DBETA
         real(dp) DLWDB
         HElement_t HIJS(0:I_V)
         INTEGER NMEM(:)
         type(egp) ExcitGen(0:I_V)
         INTEGER ICLS
         real(dp) ORHOII(0:I_OVCUR)
         HElement_t ORHOIJ(0:I_OVCUR,0:I_OVCUR)
         INTEGER OIPATH(NEL,0:I_OVCUR)
         HElement_t OHIJS(0:I_OVCUR)
         real(dp) OXIJ(0:I_OVCUR-1,0:I_OVCUR-1),OPROB,PR,R2
         real(dp) OETILDE,WEIGHT,OWEIGHT,ETILDE
         INTEGER IACC
         real(dp) PFAC
         real(dp) RAN2
         real(dp) MPEs(2:i_VMax)
         
!C.. Take a copy of the old path and rho matrix etc.
         IF(OWEIGHT.ne.0.0_dp) THEN
            CALL DCOPY(I_OVCUR+1,RHOII,1,ORHOII,1)
            CALL DCOPY((I_OVCUR+1)**2,RHOIJ,1,ORHOIJ,1)
            CALL DCOPY(I_OVCUR+1,HIJS,1,OHIJS,1)
            CALL DCOPY(I_OVCUR*I_OVCUR,XIJ,1,OXIJ,1)
            CALL NECI_ICOPY(NEL*(1+I_OVCUR),IPATH,1,OIPATH,1)
        ENDIF
         IF(.not.abs(INWI).gt.0.0_dp) THEN
!.. Generate a new graph
            RHX=-1
            IC=0
            TLOG=BTEST(ILOGGING,2)
            TLOG2=BTEST(ILOGGING,3)
            TLOG3=BTEST(ILOGGING,6)
            DO WHILE (RHX.LE.PGenEpsilon)
               CALL FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,I_V,XIJ,        &
     &           NBASISMAX,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,              &
     &           NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,NMEM,I_HMAX,     &
     &           EXCITGEN)
!C.. Now get the parameters for the new path 
            
!C.. GETPATHPROB gives us the probability of generating the path
               RHX=GETPATHPROB2(XIJ,I_V)*PFAC
               IC=IC+1
               IF(IC.GT.10000) THEN
                  WRITE(6,*) "Have thrown away 10000 graphs with pGen<",   &
     &               PGenEpsilon 
              STOP "Have thrown away 10000 graphs with pGen<pGenEpsilon"
               ENDIF
            ENDDO
            ICLS=0
            IF(TMPTHEORY) THEN
! MP Theory requires HDIAG
               MPEs=(0.0_dp)
               CALL AddMPEnergy(RHOIJ,i_V,i_vmax,NMAX,nBasis,              &
     &            iPath,nEl,.false.,ECORE,MPEs)
               ETILDE=0.0_dp
               DO I=2,I_VMAX
                  ETILDE=ETILDE+MPEs(I)
               ENDDO
               WEIGHT=1.0_dp
            ELSE
!C.. CALCPATHS gives us the contribution of the path
               WEIGHT=CALCPATHS_N(RHOII,RHOIJ,I_V,I_HMAX,                  &
     &            I_P,0.0_dp,DBETA,ETILDE,HIJS,ICLS)
            ENDIF
         ELSE
!C.. we use the gestalt
            RHX=PFAC
            ICLS=1
            WEIGHT=INWI
            ETILDE=DLWDB
         ENDIF
!C.. Now do an accept/reject
         IF(.NOT.abs(OWEIGHT).gt.0.0_dp) THEN
!C.. this is the first time round, so we automatically accept the new
!C.. configuration
            OETILDE=ETILDE
            OPROB=RHX
            I_OVCUR=I_V
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            IACC=1
         ELSE
!C         CALL WRITEPATH(6,IPATH,I_V,NEL,.FALSE.)
!C.. Accept the new configuration with prob (|S'(new)|/|S'(old)|)*(OPROB/RH)  

!C.. acc(A->B)/acc(B->A) = p(B) gen(B->A) / (p(A) gen(A->B))
!C..                     = p(B) gen(A) / (p(A) gen(B))

            PR=(OPROB/RHX)*ABS((WEIGHT/OWEIGHT))
!            WRITE(40,"(5G)",advance='no') OPROB,RH,(WEIGHT), (OWEIGHT),PR
            IF(RHX.EQ.0.0_dp.OR..NOT.abs(OWEIGHT).gt.0.0_dp) PR=0.0_dp
            IF(isnan_neci(WEIGHT)) THEN
                WRITE(60,*) WEIGHT,ETILDE,RHX
                CALL WRITEPATH(60,IPATH,I_V,NEL,.FALSE.)
                CALL WRITERHOMAT(60,RHOIJ,I_V,.TRUE.)
                CALL WRITE_XMATRIX(60,XIJ,I_V)
            ENDIF
            R2=RAN2(ISEED)
            IF(PR.LT.R2) THEN
                IACC=0
            ELSE
                IACC=1
            ENDIF
!            CALL WRITEPATHEX(40,OIPATH,I_OVCUR,NEL,.FALSE.)
!            WRITE(40,"(A)",advance='no'),"--->"
!            CALL WRITEPATHEX(40,IPATH,I_V,NEL,.FALSE.)
!            WRITE(40,*) IACC
         ENDIF
         CALL CLASSPATHS(FMCPR4D2,DLWDB,1.0_dp,RHOIJ,                              &
     &            I_V,ICLS)
! If we're writing graphs, we write here before acceptance/rejectance
         IF(BTEST(ILOGGING,2)) THEN
            IF(BTEST(ILOGGING,12)) THEN
! We log excitations
               CALL WRITEPATHEX(10,IPATH,I_V,NEL,.FALSE.)
            ELSE
! We log dets
               CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            ENDIF
            IF(BTEST(ILOGGING,3))                                             &
     &            CALL WRITERHOMAT(10,RHOIJ,I_V,.TRUE.)
            IF(BTEST(ILOGGING,6)) CALL WRITE_XMATRIX(10,XIJ,I_V)
            WRITE(10,"(3E25.16, 2I7)") WEIGHT,RHX,ETILDE,ICLS,IACC
         ENDIF
          
         IF(IACC.EQ.0) THEN
!C.. reject the new, and copy in the old
            ETILDE=OETILDE
            RHX=OPROB
            I_V=I_OVCUR
            ICLS=IOCLS
            WEIGHT=OWEIGHT
            CALL NECI_ICOPY(NEL*(1+I_OVCUR),OIPATH,1,IPATH,1)
         ELSE
            OETILDE=ETILDE
            OPROB=RHX
            I_OVCUR=I_V
!               ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            IACC=1
         ENDIF
! If we're not writing graphs, we write here after
         IF(BTEST(ILOGGING,0).AND..NOT.BTEST(ILOGGING,2)) THEN
            WRITE(10,"(3E25.16, 2I7)") WEIGHT,RHX,ETILDE,ICLS,IACC
         ENDIF
         IF(TMPTHEORY.AND..NOT.tMCDirectSum) THEN
            DLWDB=ETILDE/(OPROB)
            FMCPR4D2=OPROB
         ELSEIF(tMCDirectSum.AND.RHX.GT.0.0_dp.AND.                          &
     &                                   (abs(WEIGHT).gt.0.0_dp)) THEN
!C.. Unbias the sum 
            DLWDB=ETILDE/WEIGHT
            FMCPR4D2=WEIGHT/WEIGHT
         ELSEIF(RHX.GT.0.0_dp.AND.abs(WEIGHT).gt.0.0_dp) THEN
!C.. Unbias the sum 
            DLWDB=ETILDE/(ABS((WEIGHT)))
            FMCPR4D2=WEIGHT/(ABS((WEIGHT)))
         ELSE
            FMCPR4D2=0.0_dp
            DLWDB=0.0_dp
         ENDIF
         CALL GETTREENESS(ICLS,ITREE,FMCPR4D2,I_V)
      END function fmcpr4d2











!.. Generate a new graph.  Called from FMCPR4D2
!.. NMAX has ARR  hidden in it
      SUBROUTINE FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,I_V,XIJ,        &
     &   NBASISMAX,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMat,                   &
     &   NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,NMEMNI,I_HMAX,        &
     &   PVERTMEMS2)
         use CalcData , only : G_VMC_PI
         use Determinants, only: get_helement, write_det
         use SystemData, only: BasisFN,Arr
         use mcpathsdata, only: egp
         IMPLICIT NONE
         INTEGER NEL,I_V,IPATH(NEL,0:I_V),NI(NEL)
         real(dp) BETA,ALAT(*),ECORE
         HElement_t Umat(*)
         TYPE(BasisFN) G1(*)
         INTEGER I_P,nBasisMax(5,*),NBASIS,NMSH
         INTEGER NMAX,NTAY(2)
         type(egp) PVERTMEMS2(0:I_V)
         complex(dp) FCK(*)
         INTEGER I_HMAX
         real(dp) RHOEPS
         real(dp) RHOII(0:I_V)
         HElement_t RHOIJ(0:I_V,0:I_V)
         real(dp) XIJ(0:I_V-1,0:I_V-1)         
         HElement_t HIJS(0:I_V)
         INTEGER ISEED
         real(dp) PEXCIT(I_V)
         type(egp) PVERTMEMS(0:I_V)
         INTEGER I_VNEXT
         INTEGER INODE(NEL),INODE2(NEL)
         INTEGER, target :: NMEMNI(:)
         INTEGER, target :: NEWEXLEN(1)
         INTEGER NEXNODE
         integer, pointer :: curex(:)
         integer, pointer :: newex(:) => null()
         INTEGER STORE(6)
         real(dp) R
         HElement_t RH
         INTEGER ICE,I,ICOUNT,IC
         INTEGER IGETEXCITLEVEL
         LOGICAL LISINPATH
         LOGICAL ISVALIDDET
         real(dp) RAN2
         real(dp) pGen,pGen2
         HElement_t :: hel

         !Deallocate Excitation Generators after precalc
         IF(I_HMAX.EQ.154) THEN
             DO I=1,I_V-1
                 deallocate(PVERTMEMS2(I)%p)
             ENDDO
             RETURN
!          STOP 'I_HMAX of 154 should not be called here (mcpathsismc.F)'
          ENDIF
!C.. PEXCIT(NODE) is the probability of selecting NODE and
!C.. those before it.
         IF(I_HMAX.EQ.-19) THEN
!.. We're actually storing H matrices in RHOIJ rather than rho matrices
            hel = get_helement (nI, nI, 0)
            RHOIJ(0,0) = hel
            RHOII(1)=BETA
         ENDIF
         PEXCIT(1)=1.0_dp
         !LOC(NMEMNI) is a pointer to NMEM
         PVERTMEMS(0)%p=>NMEMNI
         I_VNEXT=1
         CALL NECI_ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0

!C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
!C.. doesn't bother calculating the RHO_JJ.

!24/04/07 - ghb24 - Removing (hopefully) redundant routine calls - now
!       that there is a separate weighting routine, it is not needed...
!          CALL GENRANDSYMEXCITIT2(NI,NEL,
!     &         NMEMNI,INODE2,ISEED,IC,pGen)
!C.. Diagonal X elements contain the normalization of that node
!C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         ICOUNT=0
         DO WHILE(I_VNEXT.LT.I_V)
!C.. pick a random excitation of where we are (INODE) of the appropriate
!C.. symmetry.  We don't worry about weighting this at the moment
!C.. We pick which node we'regoing to excite from
            IF(I_VNEXT.GT.1) THEN
                IF(G_VMC_PI.eq.-1.0_dp) THEN
                    
                    DO I=1,I_VNEXT
                        PEXCIT(I)=1.0_dp/(I_VNEXT+0.0_dp)
                    ENDDO
                ELSE
                    DO I=2,I_VNEXT
                        PEXCIT(I)=(1.0_dp-G_VMC_PI)/(I_VNEXT-1.0_dp)
                    ENDDO
                    PEXCIT(1)=G_VMC_PI
                ENDIF
            ENDIF
!C.. Just make really sure the last one catches everytihng
            PEXCIT(I_VNEXT)=2.0_dp
            R=RAN2(ISEED)
            NEXNODE=0
            DO WHILE(R.GE.0)
               NEXNODE=NEXNODE+1
               R=R-PEXCIT(NEXNODE)
            ENDDO
! Vertex labels start at 0
            NEXNODE=NEXNODE-1
!C.. set this chosen node to be the pivot
            CUREX=>PVERTMEMS(NEXNODE)%p

!C.. generate a random excitation
            CALL GENRANDSYMEXCITIT2(IPATH(1,NEXNODE),NEL,         &
     &            CUREX,INODE2,ISEED,IC,pGen)


            CALL NECI_ICOPY(NEL,INODE2,1,INODE,1)

!C.. If the new node is not in the path, add it.
            IF(.NOT.ISVALIDDET(INODE,NEL)) THEN
               WRITE(6,*) "INVALID DET"
               call write_det (6, INODE, .true.)
               STOP "INVALID DET"
            ENDIF
            IF(.NOT.LISINPATH(INODE,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL NECI_ICOPY(NEL,INODE,1,IPATH(1,I_VNEXT),1)

!C.. Setup the spin excit generator
               STORE(1)=0
               !determine how much memory needed
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,           &
     &            .TRUE.,NEWEXLEN(1),INODE2,IC,STORE,3)
               !allocate memory
               nullify(newex)
               allocate(NewEx(NewExLen(1)))
               NEWEX(1)=0
               PVERTMEMS(I_VNEXT)%p=>NEWEX
               !Generate Excitation generator
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,           &
     &            .TRUE.,NEWEX,INODE2,IC,STORE,3)
!C.. Count the excitations (and generate a random one which we throw)
 
!24/04/07 - ghb24 - Removing (hopefully) redundant routine calls - now
!       that there is a separate weighting routine, it is not needed..
!               CALL GENRANDSYMEXCITIT2(INODE,NEL,
!     &            NEWEX,INODE2,ISEED,IC,pGen2)

               XIJ(I_VNEXT,I_VNEXT)=IC
!C.. Update the rho and X (probability) matrices with this new node
               IF(I_HMAX.EQ.-19) THEN
!  using H-elements in the rho matrix
                 RH = get_helement (iNode, iNode, 0)
                  RHOIJ(I_VNEXT,I_VNEXT)=RH
               ELSE
                  CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,         &
     &            G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,              &
     &             RH,NTAY,0,ECORE)
                  RHOII(I_VNEXT)=RH
                  RHOIJ(I_VNEXT,I_VNEXT)=RH
               ENDIF
               DO I=0,I_VNEXT-1
                  ICE=IGETEXCITLEVEL(INODE,IPATH(1,I),NEL)
                  IF(ICE.LE.2) THEN
!  Find the prob that we would've generated this node INODE from any of the previous.
                     IF(I.EQ.NEXNODE) THEN
!  If this is the node we actually generated from, we have the prob in pGen, otherwise we calculate it
                        pGen2=pGen
                     ELSE
                        CUREX=>PVERTMEMS(I)%p
!  NMAX is really Arr
                        CALL GenExcitProb(IPATH(1,I),INODE,nEl,      &
     &                   CUREX,G1,nBasisMax,Arr,nBasis,pGen2)
                     ENDIF
                     XIJ(I,I_VNEXT)=pGen2
!  Now work out the reverse prob of generation.
!  NMAX is really Arr
                     CALL GenExcitProb(INODE,IPATH(1,I),nEl,      &
     &                NEWEX,G1,nBasisMax,Arr,nBasis,pGen2)
                     XIJ(I_VNEXT,I)=pGen2
      
! We supercede the old unbiased random choice with a new one.
!                     XIJ(I_VNEXT,I)=1.0_dp/IC
!                     IF(ISUHFDET(IPATH(1,I),NEL,NBASISMAX)) THEN
!                        IF(ICE.EQ.2) THEN
!                           XIJ(I,I_VNEXT)=1.0_dp/XIJ(I,I)
!                        ELSE
!                           XIJ(I,I_VNEXT)=0.0_dp
!                        ENDIF
!                     ELSE
!                        XIJ(I,I_VNEXT)=1.0_dp/XIJ(I,I)
!                     ENDIF   
                  ELSE
                     XIJ(I_VNEXT,I)=0.0_dp
                     XIJ(I,I_VNEXT)=0.0_dp
                  ENDIF
                  IF(I_HMAX.EQ.-19) THEN
!  using H-elements in the rho matrix
                     RH = get_helement (iNode, iPath(:,I))
                     IF(I.EQ.0) HIJS(I_VNEXT)=RH
                  ELSE
                     IF(I.EQ.0) THEN
                         hel = get_helement (nI, iNode)
                         HIJS(I_VNEXT) = hel
                     ENDIF
                     CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,    &
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,                &
     &             RH,NTAY,ICE,ECORE)
                  ENDIF
                  IF(abs(RH).gt.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.0_dp
                     RHOIJ(I_VNEXT,I)=0.0_dp
                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.50000) THEN
!C.. give up trying to find something to attach, and go home
               WRITE(6,*) "WARNING: Unable to find attachee to vertex"
               call write_det (6, IPATH(1,NEXNODE), .true.)
               RHOII(1)=0.0_dp
               RETURN
            ENDIF
         ENDDO
         DO I=1,I_V-1
             deallocate(PVERTMEMS(I)%p)
         ENDDO
!C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
!C.. matrix elements.
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1,I_V),1)
      END subroutine fmcpr4d2gengraph


      SUBROUTINE GETTREENESS(ICLASS,ITREE,WEIGHT,I_V)
         IMPLICIT NONE
         INTEGER ICLASS,ITREE,I_V,N,ICL
         real(dp) WEIGHT
         ITREE=0
         IF(WEIGHT.EQ.0.0_dp) RETURN
         N=0
         ICL=ICLASS
         DO WHILE (ICL.NE.0)
            IF(IAND(ICL,1).EQ.1) N=N+1
            ICL=ICL/2
         ENDDO
         IF(N.EQ.I_V-1) ITREE=1
         RETURN
      END subroutine gettreeness

end module

!  Calculate the generation probability of a given graph for debugging purposes.
!  if iUnit<>0, write this to this unit
      SUBROUTINE CalcWriteGraphPGen(iUnit,iGraph,iV,nEl,G1,       &
     &            nBasisMax,Arr,nBasis,pGenGraph,EXCITGEN)
         USE symexcit2
         use constants, only: dp
         use SystemData, only: BasisFN
         use mcpathsdata, only: egp
         use mcpathsismc, only: getpathprob2
         IMPLICIT NONE
         INTEGER iUnit
         INTEGER iV,nEl       ! Vertices and number of electrons
         INTEGER iGraph(nEl,0:iV)  ! The graph
         INTEGER STORE(6)
         INTEGER nBasis,IC,INODE2(NEL)
         INTEGER, allocatable, target :: NEWEX(:)
         integer, target :: NEWEXLEN(1)
         !Array of pointers
         type(EGP) EXCITGEN(0:iV)

         TYPE(BasisFN) G1(nBasis)
         INTEGER nBasisMax(*)
         real(dp) Arr(nBasis,2)
         real(dp) pGen, pGenGraph


         real(dp) XIJ(0:iV-1,0:iV-1)  ! XIJ(I,J) is the pgen of J from I
         INTEGER i,j,ISEED
         ISEED=0
         
        DO i=0,iV-1
!  First work out which exictor we need to work out the excitation i->?
!C            DO j=0,iV-1
!C!  See if vertex j was excited from vertex i.  If so, then LOCTAB(1,j) is the excitor for vertex i
!C               IF(LOCTAB(3,j).EQ.i) THEN
!C                  IP_NEWEX=LOCTAB(1,j)
!C                  EXIT
!C               ENDIF
!C            ENDDO
!C            IF(j.EQ.iV-1) THEN
!We've not found an excitor, so we create one
            STORE(1)=0
!C.. Setup the excit generator for the last vertex
            CALL GENSYMEXCITIT2(iGraph(1,i),NEL,G1,NBASIS,           &
     &              .TRUE.,NEWEXLEN,INODE2,IC,STORE,3)
            allocate(NEWEX(NEWEXLEN(1)))
            NEWEX(1)=0
!C.. Count the excitations (and generate a random one which we throw)
            CALL GENSYMEXCITIT2(iGraph(1,i),NEL,G1,NBASIS,              &
     &              .TRUE.,NEWEX,INODE2,IC,STORE,3)
            CALL GENRANDSYMEXCITIT2(iGraph(1,i),NEL,                    &
     &               NEWEX,INODE2,ISEED,IC,PGEN)
!
            DO j=0,iV-1
                IF(i.EQ.j) THEN
                    XIJ(i,i)=IC
                ELSE
                    CALL GenExcitProb(iGraph(1,i),iGraph(1,j),nEl,      &
     &                      NEWEX,G1,nBasisMax,Arr,nBasis,pGen)
                    XIJ(i,j)=pGen
                ENDIF
            ENDDO
            deallocate(NEWEX)
        ENDDO
        pGenGraph=GETPATHPROB2(XIJ,iV)
         
         IF(iUnit.NE.0) CALL WRITE_XMATRIX(iUnit,XIJ,iV)

      END SUBROUTINE CALCWRITEGRAPHPGEN
