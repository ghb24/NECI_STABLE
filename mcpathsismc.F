
C.. Calculate RHO^(P)_II without having a stored H matrix
C.. SAMPLE over distinct nodes, e.g. IJKLI, with paths up to I_HMAX
C.. generated from these, and summed (e.g IJILKJI), up to max H
C.. In theory more efficient because RHO_IJ,RHO_JK, etc are calculated
C.. just once for all these paths.
C.. I_VMAX is the max number of distinct vertices in a path.
C.. I_HMAX is the max number of hops in a path.
C.. NWHTAY contains the number of samples to take for each level
C.. This is a simple importance sampling method
C.. NMAX has ARR hidden in it.
      SUBROUTINE  MCPATHSR4(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,
     &              NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &               LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,
     &               WLRI,WLSI,DBETA,DLWDB,I_VMIN)
         USE HElement
         USE MCStats
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(MCStats) MCStat
         INTEGER IPATH(NEL,0:I_VMAX)
         INTEGER NI(NEL),I_P,I_HMAX,NMSH,NMAX,NTAY,NWHTAY
         INTEGER ILOGGING,ISUB,I,J,K,I_VMIN
         CHARACTER*20 STR,STR2
         TYPE(HDElement) TOTAL,RHOII(0:I_VMAX)
         TYPE(HElement) RHOIJ(0:I_VMAX,0:I_VMAX),RH
         TYPE(HDElement) FLNSUM,FLNRII
         REAL*8 ALAT(3),RHOEPS,BETA
         TYPE(HElement) UMAT(*),TMat(*)
         COMPLEX*16 FCK(*)
         REAL*8 ECORE
         TYPE(HDElement) FF
         TYPE(HDElement) WLRI,WLSI
         TYPE(BasisFN) G1(*)
         INTEGER NBASISMAX(5,2),BRR(*)
         INTEGER LSTE(*),ILMAX,I_VMAX,NEL,NBASIS
CNEL,0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER ICE(*)
         TYPE(HElement) RIJLIST(*)
C0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER NLIST(0:I_VMAX-1),LSTP(0:I_VMAX-1),BTABLE(0:I_VMAX)
         LOGICAL TLOG,TSYM,LISNAN
         INTEGER ICOUNT,ISEED,I_V,ISUB2,L,LT
         REAL*8 FNCR,GETNLOOPS,DBETA
         INTEGER KSYM(5),IEXCITS,ITIME
         TYPE(HDElement) FMCPR4
         INTEGER FSEEK
         INTEGER EXCITGEN(:) 
         TYPE(HDElement) FMCPR4B,FMCPR4C,FMCPR4D2
         TYPE(HDElement) FMCPR4D3
         TYPE(HDElement) DLWDB,DLWDB2,DLWDB3
         TYPE(HElement) HIJS(0:I_VMAX)
         INTEGER INODE2(NEL)
         INTEGER NMEM(*),NMEMLEN
         POINTER (IP_NMEM,NMEM)
         TYPE(HDElement) ODLWDB,OWEIGHT,DLWDBSQ
         REAL*8 OPROB,R
         INTEGER ITREE
         INTEGER I_VCUR,I_VM1,I_VM2,I_OVCUR,IOCLS,NTR
         INTEGER ICHANGED
         REAL*8 VWEIGHTS(2,0:I_VMAX)
         INTEGER VEXCITS(0:I_VMAX),TST
         REAL*8 PFAC
         TYPE(HDElement) WMIN,SUMDLWDB
         REAL*8 RAN2
C.. do some blocking analysis
         INTEGER NN,BMAX,OBMAX
         INTEGER STORE(6)
         TYPE(HDElement) DLWDBCORE,WCORE
         INTEGER IOV,IGV,IACC
         INTEGER TIGET
         LOGICAL TSEQ,TBLOCKING
         REAL*8 PREJ,PGR
         REAL*8 TIMEC,OTIME,NTIME
         INTEGER*8 LP
         INCLUDE 'vmc.inc'
         OTIME=TIMEC()*0.001
         TST=0
         CALL TISET('MCPATHSR4 ',ISUB)
         TLOG=BTEST(ILOGGING,1)
         TBLOCKING=BTEST(ILOGGING,13)
         IF(TLOG.AND.I_VMIN.EQ.0) THEN
            OPEN(11,FILE="MCPATHS",STATUS="OLD")
C.. go to end of file
            I=FSEEK(11,0,2)
            CALL WRITEDET(11,NI,NEL,.TRUE.)
         ENDIF
         IF(BTEST(ILOGGING,2))
     &      OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,9))
     &      OPEN(12,FILE="VERTEXMC",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,10))
     &      OPEN(22,FILE="VMC",STATUS="UNKNOWN")
!         IF(tMCDirectSum)
!     &      OPEN(122,FILE="MCDIRECTSUM",STATUS="UNKNOWN")
C.. Set the first node to I_I
         CALL ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,
     &            BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(0)=RH
         RHOIJ(0,0)=RHOII(0)
         WLRI=LOG(RHOII(0))
         HIJS(0)=GETHELEMENT2(NI,NI,NEL,NBASISMAX,
     &      G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,0,ECORE)
         TOTAL=1.D0
C.. These variables have been passed in and represent the precalculated values for 1..I_VMIN-1 vertices
         DLWDBCORE=DLWDB
         WCORE=WLSI
         IF(TMPTHEORY) THEN
            DLWDBCORE=0.D0
            WCORE=0.D0
            CALL Create(MCStat,I_VMAX,
     &         INT(LOG(0.D0+NWHTAY)/LOG(2.D0)+1),DLWDBCORE,WCore)
            WRITE(6,*) "MCPATHSR4 MP Theory"
         ElseIf(tMCDirectSum) Then
            IF(I_VMIN.NE.0) THEN
                CALL Create(MCStat,I_VMIN,
     &          INT(LOG(0.D0+NWHTAY)/LOG(2.D0)+1),
     &          DLWDBCORE/WCore,WCore)
            ELSE
                CALL Create(MCStat,I_VMAX,
     &          INT(LOG(0.D0+NWHTAY)/LOG(2.D0)+1),
     &          DLWDBCORE/WCore,WCore)
            ENDIF
            WRITE(6,*) "MCPATHSR4 Direct Sum"
            !Print out Eref
!            write(6,*) DLWDBCORE/WCore
            !call flush(6)
         ELSEIF(I_VMIN.NE.0) THEN
            CALL Create(MCStat,I_VMAX,
     &         INT(LOG(NWHTAY/(G_VMC_FAC**4))/LOG(2.D0)+1),
     &      DLWDBCORE/WCORE,WCore)
            IF(G_VMC_FAC.GT.1.D0.OR.G_VMC_FAC.LT.-1.D0)
     &       STOP "Invalid MULTI MC BIAS"
            WRITE(6,*) "MCPATHSR4 MultiMC.  Bias=",G_VMC_FAC
         ELSE
            DLWDBCORE=HIJS(0)
            CALL Create(MCStat,I_VMAX,
     &         INT(LOG(0.D0+NWHTAY)/LOG(2.D0)+1),DLWDBCORE,WCore)
         ENDIF
         IF(DBETA.LT.0.D0.AND.I_VMIN.EQ.0) THEN
            DLWDB=HIJS(0)
         ENDIF
         IF(TLOG.AND.I_VMIN.EQ.0) 
     &         WRITE(11,"(I12,2G25.16,F19.7,2I12,G25.12)") 
     &         1,TOTAL,TOTAL,0.D0,1,1,DLWDB
C.. we're working in block space, so we work out our current symmetry.
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
         IF(I_HMAX.EQ.-3.OR.I_HMAX.EQ.-4) THEN
            IEXCITS=0
            ISEED=G_VMC_SEED
            CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IEXCITS,
     &      ISEED,INODE2)
         ELSE
            TOTAL=0.D0
            IF(DBETA.LT.0.D0.AND.I_VMIN.EQ.0) THEN
               DLWDB=0.D0
            ENDIF
         ENDIF
         I_VM1=2
         I_VM2=I_VMAX
         IF(I_HMAX.EQ.-7.OR.I_HMAX.LE.-12) THEN
C.. setup spin excitation generator
C            CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
C     &               .TRUE.,NMEMLEN,INODE2,ISEED,IEXCITS,0)
C            CALL MEMORY(IP_NMEM,NMEMLEN,'NMEM')
C            NMEM(1)=0
C            CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
C     &               .FALSE.,NMEM,INODE2,ISEED,IEXCITS,0)
C.. Setup the spin excit generator
            STORE(1)=0
            CALL GENSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         .TRUE.,NMEMLEN,INODE2,I,0,STORE,3)
            CALL MEMORY(IP_NMEM,NMEMLEN,'NMEM')
            NMEM(1)=0
            CALL GENSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         .TRUE.,NMEM,INODE2,I,0,STORE,3)
C.. Count the excitations (and generate a random one which we throw)
            CALL GENRANDSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         NMEM,INODE2,ISEED,IEXCITS,0,UMAT,NMAX,PGR)
            IF(I_VMIN.NE.0) THEN
               I_VM1=I_VMIN
               I_VM2=I_VMIN
               I_OVCUR=1
            ENDIF
            VEXCITS(0)=IP_NMEM
         ENDIF
         DO I_V=I_VM1,I_VM2
            CUR_VERT=I_V
            WRITE(STR,"(A,I4)") "FMCPR_",I_VMIN
            CALL TISET(STR,ISUB2)
            L=0
            LT=0
            BTABLE(0)=0
            ISEED=G_VMC_SEED
            SUMDLWDB=0.D0
            IF((I_HMAX.NE.-3.AND.I_HMAX.NE.-4).OR.I_V.LE.IEXCITS+1) THEN
C.. Do importance sampling with probabilities.  No norm needed
               DLWDB3=0.D0
               IOCLS=0
               I_VCUR=1
               NTR=0
               ITREE=1
               DLWDBSQ=0.D0
               PFAC=G_VMC_PI
               WMIN=1.D-8
               IF(I_HMAX.LT.-12.AND.I_HMAX.GT.-19) THEN
C.. Setup the first graph for memory MC
                  OWEIGHT=1.D0
                  ODLWDB=HIJS(0)
C.. The VWEIGHTS hold the weights of the current vertices.
C.. These are currently set to (for 1) (rho_jj/rho_ii)**P
C.. and for 2 (rho_jj/rhoii)**-P
                  VWEIGHTS(1,0)=1.D0
                  VWEIGHTS(2,0)=1.D0
                  VWEIGHTS(1,1)=VWEIGHTS(1,0)
                  VWEIGHTS(2,1)=VWEIGHTS(2,0)
                  CALL ICOPY(NEL,IPATH(1,0),1,IPATH(1,I_VCUR),1)
               ELSE
                  OWEIGHT=0.D0
                  ODLWDB=0.D0
               ENDIF
               IF(I_VMIN.NE.0) THEN
                  ODLWDB=DLWDBCORE
                  OWEIGHT=WCORE
                  OPROB=1   
               ENDIF
               IF(NWHTAY.LT.1) WRITE(6,*) "Warning: MC CYCLE count <0: "
     &               ,NWHTAY
               WRITE(6,*) NWHTAY, " MC Cycles"
!               WRITE(6,*) I_V
!               WRITE(6,*) G_VMC_EXCITWEIGHTS(:,CUR_VERT)
!               WRITE(6,*) G_VMC_EXCITWEIGHT(CUR_VERT)
!               CALL FLUSH(6)
               DO ICOUNT=1,NWHTAY
C                  IF(ICOUNT.EQ.135) TST=1
C                  DLWDB2=0.D0
                  IF(I_HMAX.EQ.-3) THEN
                     FF=FMCPR4B(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ELSEIF(I_HMAX.EQ.-4) THEN
                     FF=FMCPR4C(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ElseIf(tMCDirectSum) Then
                    oWeight=0.D0 !Ensure generated graph is "accepted"
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
!      WRITE(43,"(I3,4G25.16)") I_VCUR,ECORE,DLWDB2%v, ODLWDB%v,OWEIGHT%v
                            
                    FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                   ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,
     &                   ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,
     &                   IACC,0.D0,I_VMAX,EXCITGEN)
                    
!     WRITE(43,"(I3,4G25.16)") I_VCUR,ECORE,DLWDB2%v, ODLWDB%v,OWEIGHT%v
!                  CALL FLUSH(43)
                    
                    !Set up values for the AddGraph routine later:
                    ioV=i_vcur
                    igV=i_vcur
!                    i_VCur=i_v
!MPTheory has weights = 1, and does not need to be further divided.
                    IF(.NOT.tMPTheory) THEN
                     FF=oWeight/HDElement(oProb) 
                    ENDIF
                    DLWDB2=oWeight*DLWDB2/HDElement(oProb)
                  ELSE
C.. see if we want to increase the number of vertices
                    I_OVCUR=I_VCUR
                    R=RAN2(ISEED)
                    ICHANGED=0
                    IF(I_HMAX.EQ.-13) THEN
!                   !
                       IF(R.LT.0.40D0) THEN
                        IF(I_HMAX.EQ.-13.AND.I_VCUR.GT.1) THEN
                         CALL MCPR4D4DEC(IPATH,NEL,I_VCUR,RHOIJ,I_HMAX,
     &                        I_P,HIJS,VWEIGHTS,VEXCITS,G1,NBASIS,
     &                        NBASISMAX,OWEIGHT,ILOGGING,DBETA,
     &                        ODLWDB,DLWDB2,FF,IOCLS,
     &                        ITREE,IACC,ISEED,PFAC,WMIN)
                           ICHANGED=1
!                   ! 
                        ELSE
                           IF(I_VCUR.GT.1) I_VCUR=I_VCUR-1
                        ENDIF
                       ELSEIF(R.LT.0.8D0) THEN
                        IF(I_HMAX.EQ.-13.AND.I_VCUR.LT.I_VMAX) THEN
                         CALL MCPR4D4INC(IPATH,NEL,I_VCUR,RHOIJ,I_HMAX,
     &                        I_P,BETA,HIJS,VWEIGHTS,VEXCITS,G1,NBASIS,
     &                        NBASISMAX,OWEIGHT,ILOGGING,DBETA,RHOEPS,
     &                        ODLWDB,DLWDB2,FF,IOCLS,
     &                        ITREE,IACC,ISEED,BRR,NMSH,FCK,TMat,NTAY,
     &                        NMAX,ALAT,UMAT,ECORE,PFAC,WMIN)
                           ICHANGED=1
                        ELSE
                           IF(I_VCUR.LT.I_VMAX) I_VCUR=I_VCUR+1
                        ENDIF
                       ENDIF
!                   !
                    ELSEIF(I_HMAX.EQ.-7.OR.I_HMAX.EQ.-19) THEN
C.. We decide on a graph size.  See 22/8/05 #1,2
C.. Probability of graph size v is prop to G_VMC_FAC**(v)
                        I_VCUR=0
!                   !   .
                        IF(I_VMIN.EQ.0) THEN
C.. we're doing MC with lots of different vertex levels
                         IF(G_VMC_FAC.GT.0) THEN
                           DO WHILE(I_VCUR.LT.1.OR.I_VCUR.GT.I_VMAX)
                              R=RAN2(ISEED)
                              I_VCUR=I_VMAX+1+LOG(R)/LOG(G_VMC_FAC)
                           ENDDO
                           PFAC=G_VMC_FAC**(I_VCUR-1)
                         ELSEIF(G_VMC_FAC.EQ.0) THEN
                           R=RAN2(ISEED)
                           I_VCUR=R*(I_VMAX-1)+2
                           PFAC=1.D0/(I_VMAX-1)
                         ELSE
                           DO WHILE(I_VCUR.LT.1.OR.I_VCUR.GT.I_VMAX)
                              R=RAN2(ISEED)
                              I_VCUR=-FLOOR(LOG(R)/LOG(-G_VMC_FAC))
                           ENDDO
                           PFAC=ABS(G_VMC_FAC)**(-I_VCUR)
                         ENDIF
                         IOV=I_OVCUR
                         IGV=I_VCUR
C.. Markov chain with random generation.
                         FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,
     &                     NBASISMAX,G1,
     &                     NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                     ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,
     &                    ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,
     &                     IACC,0.D0,I_VMAX,EXCITGEN)
!                   !   .
                        ELSE
!..   We're doing MC where we have either the gestalt 1..I_VMIN-1 vertex object
!.. or graphs of size I_VMIN=I_VMAX

                         IF(G_VMC_FAC.LE.0.D0) THEN
!.. We're doing non-stochastic time MC but with a renormalized composite gestalt.
!.. G_VMC_FAC is the prob of generating a non-gestalt.
                           IOV=I_OVCUR
                           IF(RAN2(ISEED).GT.ABS(G_VMC_FAC)) THEN
!.. We generate a gestalt
                              I_VCUR=1
                              IGV=1
                              PFAC=1-ABS(G_VMC_FAC)
                              FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,
     &                        NBASISMAX,G1,
     &                     NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                     ECORE,ISEED,KSYM,DBETA,DLWDBCORE,HIJS,NMEM,
     &                    ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,
     &                        IACC,WCORE,I_VMAX,EXCITGEN)
                           ELSE
                              R=RAN2(ISEED)*g_MultiWeight(0)
                              I_VCUR=I_VMIN-1
                              DO WHILE (R.GE.0.AND.I_VCUR.LT.I_VMAX)
                                 I_VCUR=I_VCUR+1
                                 R=R-g_MultiWeight(I_VCUR)
                              ENDDO
                              IGV=I_VCUR
                              PFAC=ABS(G_VMC_FAC)*g_MultiWeight(I_VCUR)
     &                                       /g_MultiWeight(0)
C.. Markov chain with random generation.
                              FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,
     &                           NBASISMAX,G1,
     &                     NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                     ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,
     &                    ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,
     &                           IACC,0.D0,I_VMAX,EXCITGEN)
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
                           LP=LOG(1-RAN2(ISEED))/LOG(1-G_VMC_FAC)
                           IF(LP.GT.0) THEN
                              IGV=1
                              I_OVCUR=1
                              OWEIGHT=WCORE
                              ODLWDB=DLWDBCORE
                              DLWDB2=DLWDBCORE/WCORE
                              OPROB=1.D0-ABS(G_VMC_FAC)
                              IOCLS=0
                              FF=1.D0
                              ITREE=1
                              CALL AddGraph(MCStat,LP,1,FF,DLWDB2/FF,
     &                              WCORE,0,ITREE,1,1,1,
     &                              BTEST(ILOGGING,10),OPROB,
     &                              TBLOCKING)
!                              WRITE(39,*) 1,1,LP,1
                           ENDIF
                         ENDIF
!                   !
C                         WRITE(6,*) TSEQ,I_OVCUR,I_VCUR
                         IF(.NOT.TSEQ) THEN
!.. We're at a V-graph. Consider all possible strings of rejected V->0 processes.
!.. First work out the acceptance prob of V->0, pa, and thus the rejectance prob, pr=1-pa, which is high.
!.. The probability of a string of N of these rejectances is (p0 pr)^N, so work out how many we're going to
!.. get in a row
                           IF(TVVDISALLOW) THEN
                              PREJ=(DREAL(WCORE)/ABS(DREAL(OWEIGHT)))*
     &                           OPROB
                              PREJ=1.D0-MIN(PREJ,1.D0)
                              PGR=PREJ*(1.D0) !-G_VMC_FAC)
                           ELSE
                              PREJ=(DREAL(WCORE)/ABS(DREAL(OWEIGHT)))*
     &                           OPROB/(1.D0-G_VMC_FAC)
                              PREJ=1.D0-MIN(PREJ,1.D0)
                              PGR=PREJ*(1.D0-G_VMC_FAC)
                           ENDIF
C.. L is the number of strings of rejections (i.e. of V's we add
                           LP=LOG((1-RAN2(ISEED)))/LOG(PGR)
C                           WRITE(6,*) "V->V",PREJ,L,I_OVCUR,DLWDB2
                           IF(LP.GT.0) THEN
                              CALL AddGraph(MCStat,LP,I_OVCUR,FF,
     &                          DLWDB2/FF,OWEIGHT,IOCLS,ITREE,0,I_OVCUR,
     &                              1,BTEST(ILOGGING,10),OPROB,
     &                          BTEST(ILOGGING,13))
!                              WRITE(39,*) I_OVCUR,1,LP,0
                           ENDIF
C.. Now if we generate a 1, we accept with prob p1 pa (out of p1 pa + pV=1-(1-pV)*pr), but we need to consider any V's separately, with prob p1
!                   !
!Either we've disallowed V->V' transitions or we generate a 1-v graph
                           IF(TVVDISALLOW.OR.
     &   RAN2(ISEED)*(1-(1-G_VMC_FAC)*PREJ).GE.ABS(G_VMC_FAC)) THEN
                              FF=1.D0
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
!                              CALL AddGraph(MCStat,LP,1,FF,DLWDB2/FF,
!     &                           OWEIGHT,ITREE,1,I_VMIN,IGV)
                              I_VCUR=1
                              IF(TVVDISALLOW) THEN
                                 OPROB=1
                              ELSE
                                 OPROB=1.D0-ABS(G_VMC_FAC)
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
                              PFAC=ABS(G_VMC_FAC)*g_MultiWeight(I_VCUR)
     &                                       /g_MultiWeight(0)
C                           WRITE(6,*) IOV,I_VCUR,I_OVCUR
                              IF(TVVDISALLOW) OPROB=1
                              FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,
     &                         NBASISMAX,G1,NBASIS,BRR,
     &                         NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &                         RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                         ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,
     &                         ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,
     &                         PFAC,IACC,0.D0,I_VMAX,EXCITGEN)
!                              WRITE(39,*) IOV,IGV,1,IACC
C                              WRITE(6,*) IACC
                         ENDIF
                        ENDIF
                       ENDIF
!                   !
                    ELSEIF(I_HMAX.EQ.-12) THEN
C.. Markov chain with special generation - does not work
                        FF=FMCPR4D3(NI,BETA,I_P,IPATH,I_VCUR,NEL,
     &                  NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,
     &                  ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT,
     &                  IACC)
!                   !
                    ELSEIF(I_HMAX.EQ.-13.AND.ICHANGED.EQ.0) THEN
C.. Real Markov Chain MC with memory - i.e. move one node from the last
C.. graph to generate the new one.
                        CALL MCPR4D4MOVE(NI,BETA,I_P,I_VCUR,NEL,
     &                     NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,
     &                     ALAT,UMAT,NTAY,RHOEPS,NWHTAY,I_HMAX,ILOGGING,
     &                     ECORE,ISEED,KSYM,DBETA,DLWDB2,IPATH,RHOIJ,
     &                     HIJS,FF,ODLWDB,IOCLS,ITREE,OWEIGHT,IACC,
     &                     VWEIGHTS,VEXCITS,PFAC,TST,WMIN)
C.. FF for ETILDE, ODLWDB for OETILDE, 
!                   !
                    ENDIF
                  ENDIF
                  LP=1
                  If (FF%v.eq.0) Then
         CALL AddGraph(MCStat,LP,I_VCUR,HDElement(0.D0),HDElement(0.D0)
     &               ,OWEIGHT,IOCLS, ITREE,
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),OPROB,
     &               TBLOCKING)
                Else

              CALL AddGraph(MCStat,LP,I_VCUR,FF,DLWDB2/FF,OWEIGHT,IOCLS,
     &               ITREE,
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),OPROB,
     &                  TBLOCKING)
                  EndIf
                  IF(BTEST(ILOGGING,9).AND.
     &              MOD(ICOUNT,G_VMC_LOGCOUNT).EQ.0) THEN
                 !    CALL WriteStats(MCStat,12)
                     If(tMCDirectSum) Then
                       Write(12, "I9, 6G25.16") iCount,
     &                 MCStat%wWeightedDelta(0)%v
     &                 /(MCStat%wWeighting(0)%v+(iCount+1)*WCore%v),
     &                 MCStat%wWeightedDelta(0)%v,
     &                 MCStat%wWeighting(0)%v, (iCount+1)*WCore%v,
     &      MCStat%wWeightedDeltaSq(0)%v/iCount
     &             -(MCStat%wWeightedDelta(0)%v/iCount)**2,
     &      MCStat%wWeightingSq(0)%v/iCount
     &             -(MCStat%wWeighting(0)%v/iCount)**2
                     Else
                       Write(12, "I9, 3G25.16") iCount, 
     &                  MCStat%wWeightedDelta(0)/MCStat%wWeighting(0),
     &                  MCStat%wWeightedDelta(0),MCStat%wWeighting(0)
                     EndIf
                     CALL FLUSH(12)
                  ENDIF
               ENDDO   ! end loop over mc cycles  nwhtay
            ENDIF
            CALL TIHALT(STR,ISUB2)
!            F(I_V)=F(I_V)/NWHTAY
!            FSQ(I_V)=FSQ(I_V)/NWHTAY
!            DLWDB=DLWDB+DLWDB3/NWHTAY
!            STD=SQRT(ABS(DREAL(FSQ(I_V)-F(I_V)*F(I_V))))
!            TOTAL=TOTAL+F(I_V)
c            WRITE(6,*) I_V,F(I_V),TOTAL,TIGET(ISUB2),L,LT
!            IF(TLOG)
!     &         WRITE(11,"(I12,2G25.16,G19.7,I12,2G19.7)") I_V,F(I_V),
!     &            TOTAL,TIGET(ISUB2),L,STD,DLWDB3/NWHTAY
C ,F(I_V),FSQ(I_V)
            NTIME=TIMEC()*0.001
            IF(TLOG.AND.(I_HMAX.EQ.-7.OR.I_HMAX.LE.-12)) THEN
              IF(I_VM1.EQ.I_VM2) THEN
!.. just give the additional components for this vertex level
             CALL WriteLongStats2(MCStat,11,WCORE,DLWDBCORE,NTIME-OTIME)
               ELSE
!.. rejig the sums so the result is Sum w Delta / Sum w
              CALL WriteLongStats(MCStat,11,WCORE,DLWDBCORE,NTIME-OTIME)
               ENDIF
            ENDIF
            CALL FLUSH(11)
!            IF(LISNAN(F(I_V))) THEN
C.. save all log files
!               ITIME=TIMEC()
!               CALL FLUSH(11)
caa            CALL FLUSH(10)
!               CALL LOGNAN(NI,NEL,BETA,ITIME)
!               WRITE(6,*) "WARNING: nan found at time",ITIME
!               WRITE(6,"(A,$)") "  nan det=("
!               DO K=1,NEL
!                  WRITE(6,"(I3,A,$)") NI(K),","
!               ENDDO
!               WRITE(6,"(A,$)") "),"
!            ENDIF
         ENDDO
! Finish off any sequences
         LP=0
         IOV=I_VCUR
         I_VCUR=0
         IGV=0
         IACC=1
         If (FF%v.eq.0) Then
         CALL AddGraph(MCStat,LP,I_VCUR,HDElement(0.D0),HDElement(0.D0)
     &               ,OWEIGHT,IOCLS, ITREE,
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),0.D0,
     &               BTEST(ILOGGING,13))
         Else
         CALL AddGraph(MCStat,LP,I_VCUR,FF,DLWDB2/FF,OWEIGHT,IOCLS,
     &               ITREE,
     &               IACC,IOV,IGV,BTEST(ILOGGING,10),0.D0,
     &               BTEST(ILOGGING,13))
         EndIf
         IF(I_HMAX.EQ.-7.OR.I_HMAX.LE.-12) THEN
C.. setup spin excitation generator
            CALL FREEM(IP_NMEM)
         ENDIF
         IF(TLOG.AND.I_VMIN.EQ.0) CLOSE(11)
         IF(BTEST(ILOGGING,2)) CLOSE(10)
         IF(BTEST(ILOGGING,9)) CLOSE(12)
         IF(BTEST(ILOGGING,10)) CLOSE(22)
         IF(BTEST(ILOGGING,13)) THEN
            Call CalcStDev(MCStat,R)
            WRITE(6,"(A,I4,A,G)") "Error in MC ",I_VM2,":",R
            CALL FLUSH(6)
        ENDIF
         IF(TMPTHEORY) THEN
            Call GetStats(MCStat,0,TOTAL,WLSI,DLWDB2)
            WLSI=LOG(TOTAL)
            DLWDB=DLWDB2
            FF=HIJS(0)
            FF=FF+DLWDB/TOTAL
!  If we're doing a full sum with all levels, then we print data out here.
!  Otherwise it will be done in MCPATHSR10 which calls us
            IF(I_VMIN.EQ.0) THEN
               WRITE(STR,"(A2,I1,A10)") "MP",I_VMAX," ENERGY = "
               WRITE(6,"(A14,G)") STR,FF
            ENDIF
!            DO I=2,I_VMAX
!               FF=FF+MP2E(I)
!               WRITE(STR,"(A2,I1,A10)") "MP",I," ENERGY = "
!               WRITE(6,"(A14,G)") STR,FF
!            ENDDO
         ELSEIF(I_VMIN.EQ.0.OR.tMCDirectSum) THEN
            IF(I_VMIN.NE.0) THEN  !We're just doing this vertex level, not the complete sum
! GetStats(MCS,iV,wAvgWeighting,wAvgWeightedValue,wAvgDelta)
               Call GetStats(MCStat,0,TOTAL,DLWDB2,WLSI)
!               WRITE(6,"(A,I3,4G)") "MCStat:",I_VMIN,TOTAL,WLSI,
!     &            DLWDB2,DLWDB
               TOTAL=TOTAL
               WLSI=TOTAL
!WLSI gets the weight of this level
!DLWDB has the  w E~ sum of past levels, and we add the w E~ of this level to it.
               DLWDB=DLWDB+DLWDB2
            ELSE
               WLSI=1.D0 !weight of 1v graph
               TOTAL=TOTAL+WLSI
               DLWDB=HIJS(0)
               DO I_V=I_VM1,I_VM2
                  Call GetStats(MCStat,I_V,FF,DLWDB2,WLSI)
                  TOTAL=TOTAL+FF
                  DLWDB=DLWDB+DLWDB2 
               ENDDO
               WLSI=LOG(TOTAL)
               DLWDB=DLWDB/TOTAL   !DLWDBCORE is Etilde of 1v graph
            ENDIF
         ELSE
            Call GetStats(MCStat,0,TOTAL,DLWDB2,WLSI)
            DLWDB=DLWDB2
            IF(TOTAL%v.LT.0.D0) THEN
               TOTAL%v=-TOTAL%v
               DLWDB%v=-DLWDB%v
            ENDIF
            DLWDB=DLWDB
            WLSI=TOTAL-WCORE
         ENDIF
         CALL TIHALT('MCPATHSR4 ',ISUB)
         CALL Delete(MCStat)
         RETURN
      END


C.. A function which chooses a random set of I_V connected dets, working out
C.. loop contribution for that set.
C.. All nodes are distinct.  Paths IJIKJI etc.
C.. are generated by permutation from IJKI, and summed up to length I_HMAX
C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
C.. This function assumes that there are enough available excitations to 
C.. form a loop of length I_V.  If not it will probably hang.
      REAL*8 FUNCTION FMCPR4(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         REAL*8 RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER I,J,K,INODE(NEL)
         INTEGER IADJ(0:I_V-1,0:I_V-1),ICE,IPATH2(0:I_V-1)
         REAL*8 RH,CALCPATHS
         INTEGER ICOUNT,IMCPR4N,IGETEXCITLEVEL,ICMPDETS
         LOGICAL BR
         LOGICAL TLOG,TLOG2
         REAL*8 RAN2
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
C.. Pick random elements of the path and generate excitations
         DO I=1,I_V-1
            K=INT(RAN2(ISEED)*I)
            BR=.TRUE.
            DO WHILE(BR)
               BR=.FALSE.
               CALL GENRANDOMEXCIT(IPATH(1,K),NEL,NBASIS,
     &            ABS(NTAY*2),ISEED,INODE)
               DO J=0,I-1
                  IF(ICMPDETS(INODE,IPATH(1,J),NEL).EQ.0) BR=.TRUE.
               ENDDO
            ENDDO
            CALL ICOPY(NEL,INODE,1,IPATH(1,I),1)
            DO J=0,I
               ICE=IGETEXCITLEVEL(INODE,IPATH(1,J),NEL)
               IF(ICE.LE.2) THEN
                  IADJ(I,J)=1
                  IADJ(J,I)=1
               ELSE
                  IADJ(I,J)=0
                  IADJ(J,I)=0
               ENDIF
               CALL CALCRHO2(INODE,IPATH(1,J),BETA,I_P,NEL,
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &            RH,NTAY,ICE,ECORE)
               IF(ABS(RH).GT.RHOEPS) THEN
                  RHOIJ(I,J)=RH
                  RHOIJ(J,I)=RH
               ELSE
                  RHOIJ(I,J)=0.D0
                  RHOIJ(J,I)=0.D0
               ENDIF
            ENDDO   
            RHOII(I)=RHOIJ(I,I)
         ENDDO
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.  We now call CALCPATHS to get the path weight
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
         ENDIF   
         FMCPR4=CALCPATHS(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,I_V,ILOGGING)
C.. Now we need to call the recursive IMCPR4N which will count all the
C.. ways we could possibly generate this path.
         IPATH2(0)=0
C.. As we traverse the paths, we set the diagonal elements to 0 to
C.. indicate we have traversed a particular vertex
         IADJ(0,0)=0
         ICOUNT=IMCPR4N(IADJ,I_V,IPATH2,1)
         IF(TLOG) WRITE(10,*) FMCPR4,ICOUNT 
         FMCPR4=FMCPR4/ICOUNT
         RETURN
      END

      RECURSIVE INTEGER FUNCTION IMCPR4N(IADJ,I_V,IPATH,IND)
         IMPLICIT NONE
         INTEGER IADJ(0:I_V-1,0:I_V-1),I_V,IPATH(0:I_V-1),IND
         INTEGER I,J,INODE,ITOT
C.. Go through all the nodes currently in the path, and through each
C.. possible attachment for each node
         IF(IND.EQ.I_V) THEN
            IMCPR4N=1
            RETURN
         ENDIF
         ITOT=0
         DO I=0,IND-1
            INODE=IPATH(I)
            DO J=0,I_V-1
C.. If there's a connection and we haven't been to that node before
               IF(IADJ(INODE,J).NE.0.AND.IADJ(J,J).NE.0) THEN
                  IPATH(IND)=J
                  IADJ(J,J)=0
                  IADJ(INODE,J)=0
                  IADJ(J,INODE)=0
C.. add that node to the path, and recurse
                  ITOT=ITOT+IMCPR4N(IADJ,I_V,IPATH,IND+1)
                  IADJ(J,J)=1
                  IADJ(INODE,J)=1
                  IADJ(J,INODE)=1
               ENDIF
            ENDDO
         ENDDO
         IMCPR4N=ITOT
         RETURN
      END

      SUBROUTINE GENRANDOMEXCIT(NI,NEL,NBASIS,IEXLEVEL,ISEED,NJ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,IEXLEVEL,ISEED,NJ(NEL)
         INTEGER I,J,K,IEX,IEXL2
         LOGICAL BR
         INTEGER IGETEXCITLEVEL
         REAL*8 RAN2
         IF(IEXLEVEL.GT.2)
     &    STOP "Cannot handle more than double excitations."
         IF(IEXLEVEL.LT.2) IEXL2=IEXLEVEL
         IF(IEXLEVEL.EQ.2) THEN
            IEX=RAN2(ISEED)*(NBASIS-NEL)*NEL*
     &         (1+(NBASIS-NEL-1)*(NEL-1)/4.D0)
            IF(IEX.LT.(NBASIS-NEL)*NEL) THEN
               IEXL2=1
            ELSE
               IEXL2=2
            ENDIF
C            WRITE(41,"(I,$)") IEXL2
         ENDIF
         CALL ICOPY(NEL,NI,1,NJ,1)
         DO IEX=1,IEXL2
            BR=.TRUE.
C.. Find an electron we haven't excited before
            DO WHILE (BR)
               I=INT(RAN2(ISEED)*NEL)+1
               IF(NJ(I).EQ.NI(I)) BR=.FALSE.
            ENDDO
            BR=.TRUE.
            DO WHILE (BR)
               BR=.FALSE.
               J=INT(RAN2(ISEED)*NBASIS)+1
               DO K=1,NEL
C.. If the new basis fn's in our original det or our new one, we loop again
                  IF(NI(K).EQ.J.OR.NJ(K).EQ.J) BR=.TRUE.
               ENDDO
            ENDDO
            NJ(I)=J
         ENDDO
         CALL SORTI(NEL,NJ)
C         WRITE(41,*) IGETEXCITLEVEL(NI,NJ,NEL),NI(1),NI(2),NJ(1),NJ(2)
         RETURN
      END

C.. A function which chooses a random set of I_V connected dets, working out
C.. loop contribution for that set, and weighting with the appropriate
C.. probability - not importance sampling with the RHOs
C.. All nodes are distinct.  Paths IJIKJI etc.
C.. are generated by permutation from IJKI, and summed up to length I_HMAX
C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
C.. This function assumes that there are enough available excitations to 
C.. form a loop of length I_V.  If not it will probably hang.
C..
C.. All singles and doubles are allowed (even if RHO_IJ=0) so this is
C.. not very efficient.
      REAL*8 FUNCTION FMCPR4B(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         REAL*8 RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER KSYM(5)
         INTEGER LSTE(NEL),NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT
         INTEGER I,J,K,ICE,IC
         REAL*8 XIJ(0:I_V-1,0:I_V-1)         
         REAL*8 RH,CALCPATHS_N,X,GETPATHPROB
         LOGICAL LISINPATH
         INTEGER IGETEXCITLEVEL
         LOGICAL TLOG,TLOG2,TLOG3
         REAL*8 RP,DBETA,DLWDB,DLWDB2
         REAL*8 GETHELEMENT2,HIJS(0:I_V)
         LOGICAL ISUHFDET
         INTEGER ICLS
         RP=0
         I_VNEXT=1
         CALL ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0
C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
C.. doesn't bother calculating the RHO_JJ.
C         CALL GENSYMDETSSDN(NI,KSYM,NEL,G1,BRR,NBASIS,0,
C     &         IC,NBASISMAX,
C     &      BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,
C     &      RHOEPS)
         CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IC,
     &      ISEED,INODE2)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         NLIST=1
         ICOUNT=0
         DO WHILE(I_VNEXT.LT.I_V)
C.. pick a random excitation of where we are (INODE) of the appropriate
C.. symmetry.  We don't worry about weighting this at the moment
C            WRITE(57,"(A1,$)") "N"
            CALL GENRANDOMSPINEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,IC,
     &      ISEED,INODE2)
C            CALL WRITEDET(35,INODE2,NEL,.FALSE.)
            CALL ICOPY(NEL,INODE2,1,INODE,1)
C            CALL GENRANDOMEXCITSYM(INODE,NEL,NBASIS,G1,2,ISEED,
C     &         NBASISMAX,KSYM,INODE2)
C.. If the new node is not in the path, add it.
            IF(.NOT.LISINPATH(INODE2,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL ICOPY(NEL,INODE2,1,IPATH(1,I_VNEXT),1)
C.. Count the number of connections from this node (RP=0)
C               NLIST=0
C               CALL GENSYMDETSSDN(INODE2,KSYM,NEL,G1,BRR,NBASIS,LSTE,
C     &            NLIST,NBASISMAX,BETA,I_P,NMSH,FCK,ZIA,NMAX,
C     &            ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,RHOEPS)
               CALL GENRANDOMSPINEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,IC,
     &            ISEED,INODE2)

               XIJ(I_VNEXT,I_VNEXT)=IC
               NLIST=IC
CNLIST
C.. Update the rho and X (probability) matrices with this new node
               CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,0,ECORE)
               RHOII(I_VNEXT)=RH
               RHOIJ(I_VNEXT,I_VNEXT)=RH
               DO I=0,I_VNEXT-1
                  ICE=IGETEXCITLEVEL(INODE,IPATH(1,I),NEL)
                  IF(ICE.LE.2) THEN
                     XIJ(I_VNEXT,I)=1.D0/NLIST
                     IF(ISUHFDET(IPATH(1,I),NEL)) THEN
                        IF(ICE.EQ.2) THEN
                           XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
                        ELSE
                           XIJ(I,I_VNEXT)=0.D0
                        ENDIF
                     ELSE
                        XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
                     ENDIF   
                  ELSE
                     XIJ(I_VNEXT,I)=0.D0
                     XIJ(I,I_VNEXT)=0.D0
                  ENDIF
                  IF(I.EQ.0) THEN
                     HIJS(I_VNEXT)=GETHELEMENT2(NI,INODE,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &                  ICE,ECORE)
                  ENDIF
                  CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,ICE,ECORE)
                  IF(ABS(RH).GT.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.D0
                     RHOIJ(I_VNEXT,I)=0.D0
                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.500) THEN
C.. give up trying to find something to attach, and go home
               FMCPR4B=0.D0
               RETURN
            ENDIF
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
C.. GETPATHPROB gives us the probability of generating the path
         RH=GETPATHPROB(XIJ,I_V)
C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         FMCPR4B=CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,RH*NWHTAY,NEL,I_V,ILOGGING,DBETA,DLWDB,HIJS,ICLS)
         IF(TLOG) WRITE(10,"(3E25.16, I7)") FMCPR4B,RH,DLWDB,ICLS
         IF(RH.GT.0.D0) THEN
C.. Unbias the sum 
            FMCPR4B=FMCPR4B/RH
            DLWDB=DLWDB/RH
         ELSE
            FMCPR4B=0.D0
            DLWDB=0.D0
         ENDIF
C         WRITE(35,*)
         RETURN
      END

C.. A function which chooses a random set of I_V connected dets, working out
C.. loop contribution for that set, and weighting with the appropriate
C.. probability - not importance sampling with the RHOs
C.. All nodes are distinct.  Paths IJIKJI etc.
C.. are generated by permutation from IJKI, and summed up to length I_HMAX
C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
C.. This function assumes that there are enough available excitations to 
C.. form a loop of length I_V.  If not it will probably hang.
C.. This chooses vertices with weights according to the RHO_IJ values
C.. and RHO_II values.
C.. (26/09/04)
C.. If a RHOIJ is non-zero, then its weight XIJ will be non-zero.
C.. XIJ is proportional to RHO_JJ**ABS(RP). 
      REAL*8 FUNCTION FMCPR4C(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         REAL*8 RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER KSYM(5)
         INTEGER LSTE(NEL),NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT,ICURNODE
         INTEGER I,J,K,ICE,IC,IONODE
         REAL*8 XIJ(0:I_V-1,0:I_V-1)         
         REAL*8 RH,CALCPATHS_N,X,GETPATHPROB
         LOGICAL LISINPATH
         INTEGER IGETEXCITLEVEL,IISINPATH
         LOGICAL TLOG,TLOG2,TLOG3
         REAL*8 RP,PP
         REAL*8 DBETA,DLWDB,HIJS(0:I_V)
         REAL*8 GETHELEMENT2
         INTEGER ICLS
         REAL*8 RAN2
C.. we hard code RP as P/50, although this should be an empirical
C.. parameter.
         RP=-I_P/50.D0
C..I_P
         I_VNEXT=1
         CALL ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0
C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
C.. doesn't bother calculating the RHO_JJ.
         CALL GENSYMDETSSDN(NI,KSYM,NEL,G1,BRR,NBASIS,0,
     &         IC,NBASISMAX,BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,
     &         UMAT,NTAY,ECORE,X,RP,0.D0,RHOEPS)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=X
         IF(X.EQ.0.D0) THEN
            FMCPR4C=0.D0
            RETURN
         ENDIF       
C
         NLIST=1
         ICOUNT=0
         ICURNODE=0
         DO WHILE(I_VNEXT.LT.I_V)
C.. XIJ(ICURNODE,ICURNDOE) is the norm const for the current node
            PP=RAN2(ISEED)*XIJ(ICURNODE,ICURNODE)
            IC=1
            CALL GENSYMDETSSDN(IPATH(1,ICURNODE),KSYM,NEL,G1,BRR,NBASIS,
     &         INODE2,IC,NBASISMAX,BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,
     &         UMAT,NTAY,ECORE,X,RP,PP,RHOEPS)
            CALL ICOPY(NEL,INODE2,1,INODE,1)
C            CALL GENRANDOMEXCITSYM(INODE,NEL,NBASIS,G1,2,ISEED,
C     &         NBASISMAX,KSYM,INODE2)
C.. If the new node is not in the path, add it.
            IONODE=ICURNODE
            ICURNODE=IISINPATH(INODE2,IPATH,NEL,I_VNEXT,-1) 
            IF(ICURNODE.EQ.-1) THEN
               ICURNODE=I_VNEXT
               CALL ICOPY(NEL,INODE2,1,IPATH(1,I_VNEXT),1)
C.. Count the number of connections from this node (RP=0)
               NLIST=0
               CALL GENSYMDETSSDN(INODE,KSYM,NEL,G1,BRR,NBASIS,LSTE,
     &            NLIST,NBASISMAX,BETA,I_P,NMSH,FCK,ZIA,NMAX,
     &            ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,RHOEPS)
               XIJ(I_VNEXT,I_VNEXT)=X
C.. Update the rho and X (probability) matrices with this new node
               CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,0,ECORE)
               RHOII(I_VNEXT)=RH
               RHOIJ(I_VNEXT,I_VNEXT)=RH
               DO I=0,I_VNEXT-1
                  ICE=IGETEXCITLEVEL(INODE,IPATH(1,I),NEL)
                  CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,ICE,ECORE)
                  IF(ABS(RH).GT.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
C                     XIJ(I_VNEXT,I)=1.D0/XIJ(I_VNEXT,I_VNEXT)
C                     XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
                     XIJ(I_VNEXT,I)=
     &                  RHOII(I)**ABS(RP)/XIJ(I_VNEXT,I_VNEXT)
                     XIJ(I,I_VNEXT)=
     &                  RHOII(I_VNEXT)**ABS(RP)/XIJ(I,I)
                  IF(I.EQ.0)
     &               HIJS(I_VNEXT)=GETHELEMENT2(NI,INODE,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &                  ICE,ECORE)
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.D0
                     RHOIJ(I_VNEXT,I)=0.D0
                     XIJ(I_VNEXT,I)=0.D0
                     XIJ(I,I_VNEXT)=0.D0
                     IF(I.EQ.0) HIJS(I_VNEXT)=0.D0
                  ENDIF
C                  IF(ICE.LE.2) THEN
C                     XIJ(I_VNEXT,I)=RHOII(I)**RP/X
C                     XIJ(I,I_VNEXT)=RHOII(I_VNEXT)**RP/XIJ(I,I)
C                     XIJ(I_VNEXT,I)=RHOIJ(I,I_VNEXT)**RP/X
C                     XIJ(I,I_VNEXT)=RHOIJ(I,I_VNEXT)**RP/XIJ(I,I)
C                     XIJ(I_VNEXT,I)=RHOIJ(I,I_VNEXT)**RP/X
C                     XIJ(I,I_VNEXT)=RHOIJ(I,I_VNEXT)**RP/XIJ(I,I)
C                  ELSE
C                     XIJ(I_VNEXT,I)=0.D0
C                     XIJ(I,I_VNEXT)=0.D0
C                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.500) THEN
C.. give up trying to find something to attach, and go home
               FMCPR4C=0.D0
               RETURN
            ENDIF
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
C.. GETPATHPROB gives us the probability of generating the path
         RH=GETPATHPROB(XIJ,I_V)
C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         FMCPR4C=CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,RH*NWHTAY,NEL,I_V,ILOGGING,DBETA,DLWDB,HIJS,ICLS)
         IF(TLOG) WRITE(10,"(3E25.16, I7)") FMCPR4C,RH,DLWDB,ICLS
         IF(RH.GT.0.D0) THEN
C.. Unbias the sum 
            FMCPR4C=FMCPR4C/RH
            DLWDB=DLWDB/RH
         ELSE
            FMCPR4C=0.D0
         ENDIF
         RETURN
      END


C.. Like MCPATHSR4C, this selects a n-vertex graph.  it chooses which
Cvertex to excite from randomly however.
C.. 
      REAL*8 FUNCTION FMCPR4D(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS,NMEM)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         REAL*8 RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER KSYM(5)
         INTEGER LSTE(NEL),NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT
         INTEGER I,J,K,ICE,IC
         REAL*8 XIJ(0:I_V-1,0:I_V-1)         
         REAL*8 RH,CALCPATHS_N,X,GETPATHPROB2
         LOGICAL LISINPATH
         INTEGER IGETEXCITLEVEL
         LOGICAL TLOG,TLOG2,TLOG3
         REAL*8 RP,DBETA,DLWDB,DLWDB2
         REAL*8 GETHELEMENT2,HIJS(0:I_V)
         LOGICAL ISUHFDET
         REAL*8 PEXCIT(I_V),R
         INTEGER NEXNODE
         INTEGER NMEM(*)
         INTEGER CUREX(*),NEWEX(*),NEWEXLEN
         INTEGER PVERTMEMS(0:I_V)
         POINTER (IP_CUREX,CUREX),(IP_NEWEX,NEWEX)
         INTEGER ICLS
         INTEGER STORE(6)
         REAL*8 RAN2
         REAL*8 pGen,pGen2
         INCLUDE 'vmc.inc'
C.. PEXCIT(NODE) is the probability of selecting NODE and
C.. those before it.
         PEXCIT(1)=1.D0
         PVERTMEMS(0)=LOC(NMEM)
C.. we hard code RP as P/50, although this should be an empirical
C.. parameter.
         RP=0
         I_VNEXT=1
         CALL ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0
C.. Count the number of adjacent nodes to us in IC.
          CALL GENRANDSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         NMEM,INODE2,ISEED,IC,0,UMAT,NMAX,pGen)
C         CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
C     &         .FALSE.,NMEM,INODE2,ISEED,IC,0)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         NLIST=1
         ICOUNT=0
         DO WHILE(I_VNEXT.LT.I_V)
C.. pick a random excitation of where we are (INODE) of the appropriate
C.. symmetry.  The prob of picking the root is G_VMC_PI, and the rest are split equally.
C.. We pick which node we'regoing to excite from
            IF(I_VNEXT.GT.1) THEN
               DO I=2,I_VNEXT
                  PEXCIT(I)=(1.D0-G_VMC_PI)/(I_VNEXT-1.D0)
               ENDDO
               PEXCIT(1)=G_VMC_PI
            ENDIF
C.. Just make really sure the last one catches everytihng
            PEXCIT(I_VNEXT)=2.D0
            R=RAN2(ISEED)
            NEXNODE=0
            DO WHILE(R.GE.0)
               NEXNODE=NEXNODE+1
               R=R-PEXCIT(NEXNODE)
            ENDDO
! 0 is the first node
            NEXNODE=NEXNODE-1

C.. set the current to be the pivot
            IP_CUREX=PVERTMEMS(NEXNODE)


C.. generate a random excitation
            CALL GENRANDSYMEXCITIT2(IPATH(1,NEXNODE),NEL,G1,NBASIS,
     &         NBASISMAX,CUREX,INODE2,ISEED,IC,0,UMAT,NMAX,pGen)

C            CALL GENRANDSYMEXCITIT(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
C     &         NBASISMAX,.FALSE.,CUREX,INODE2,ISEED,IC,0)
            CALL ICOPY(NEL,INODE2,1,INODE,1)

C.. If the new node is not in the path, add it.
            IF(.NOT.LISINPATH(INODE,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL ICOPY(NEL,INODE,1,IPATH(1,I_VNEXT),1)





C.. Setup the spin excit generator
               STORE(1)=0
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            .TRUE.,NEWEXLEN,INODE2,IC,0,STORE,3)
               CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
               NEWEX(1)=0
               PVERTMEMS(I_VNEXT)=IP_NEWEX
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            .TRUE.,NEWEX,INODE2,IC,0,STORE,3)
C.. Count the excitations (and generate a random one which we throw)
               CALL GENRANDSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            NEWEX,INODE2,ISEED,IC,0,UMAT,NMAX,pGen2)





C.. Count the number of connections from this node (RP=0)
C               CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
C     &               .TRUE.,NEWEXLEN,INODE2,ISEED,IC,0)
C               CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
C               NEWEX(1)=0
C               PVERTMEMS(I_VNEXT)=IP_NEWEX
C               WRITE(6,*) I_VNEXT,IP_NEWEX
C               CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
C     &               .FALSE.,NEWEX,INODE2,ISEED,IC,0)

               XIJ(I_VNEXT,I_VNEXT)=IC
               NLIST=IC
C.. Update the rho and X (probability) matrices with this new node
               CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,0,ECORE)
               RHOII(I_VNEXT)=RH
               RHOIJ(I_VNEXT,I_VNEXT)=RH
               DO I=0,I_VNEXT-1
                  ICE=IGETEXCITLEVEL(INODE,IPATH(1,I),NEL)
                  IF(ICE.LE.2) THEN
!  Find the prob that we would've generated this node INODE from any of the previous.
                     IF(I.EQ.NEXNODE) THEN
!  If this is the node we actually generated from, we have the prob in pGen, otherwise we calculate it
                        pGen2=pGen
                     ELSE
                        IP_CUREX=PVERTMEMS(I)
!  NMAX is really Arr
                        CALL GenExcitProb(IPATH(1,I),INODE,nEl,
     &                   CUREX,G1,nBasisMax,UMat,NMAX,nBasis,pGen2)
                     ENDIF
                     XIJ(I,I_VNEXT)=pGen2
!  Now work out the reverse prob of generation.
!  NMAX is really Arr
                     CALL GenExcitProb(INODE,IPATH(1,I),nEl,
     &                NEWEX,G1,nBasisMax,UMat,NMAX,nBasis,pGen2)
                     XIJ(I_VNEXT,I)=pGen2
      
! We supercede the unbiased generation probability with a weighted one.
!                     XIJ(I_VNEXT,I)=1.D0/NLIST
!                     IF(ISUHFDET(IPATH(1,I),NEL)) THEN
!                        IF(ICE.EQ.2) THEN
!                           XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
!                        ELSE
!                           XIJ(I,I_VNEXT)=0.D0
!                        ENDIF
!                     ELSE
!                        XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
!                     ENDIF   
                  ELSE
                     XIJ(I_VNEXT,I)=0.D0
                     XIJ(I,I_VNEXT)=0.D0
                  ENDIF
                  IF(I.EQ.0) THEN
                     HIJS(I_VNEXT)=GETHELEMENT2(NI,INODE,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &                  ICE,ECORE)
                  ENDIF
                  CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,ICE,ECORE)
                  IF(ABS(RH).GT.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.D0
                     RHOIJ(I_VNEXT,I)=0.D0
                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.500) THEN
C.. give up trying to find something to attach, and go home
               FMCPR4D=0.D0
               RETURN
            ENDIF
         ENDDO
         DO I=1,I_V-1
            CALL FREEM(PVERTMEMS(I))
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
C.. GETPATHPROB gives us the probability of generating the path
         RH=GETPATHPROB2(XIJ,I_V)
C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         FMCPR4D=CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,RH*NWHTAY,NEL,I_V,ILOGGING,DBETA,DLWDB,HIJS,ICLS)
         IF(TLOG) WRITE(10,"(3E25.16, I7)") FMCPR4D,RH,DLWDB,ICLS
         IF(RH.GT.0.D0) THEN
C.. Unbias the sum 
            FMCPR4D=FMCPR4D/RH
            DLWDB=DLWDB/RH
         ELSE
            FMCPR4D=0.D0
            DLWDB=0.D0
         ENDIF
C         WRITE(35,*)
C         CALL MEMORY_CHECK()
         RETURN
      END

C.. Generate a path probability for the modified path generation
C.. algorithm in FMCPR4D

      REAL*8 FUNCTION GETPATHPROB2(XIJ,I_V)
         IMPLICIT NONE
         INTEGER I_V
         REAL*8 XIJ(I_V,I_V)
         INTEGER IPATH(I_V)
         REAL*8 M(I_V,I_V)
         REAL*8 RET
         REAL*8 PI,PJ
         REAL*8 INV(I_V,I_V)
         INCLUDE 'vmc.inc'
         IF(I_V.EQ.2) THEN
            GETPATHPROB2=XIJ(1,2)
         ELSEIF(I_V.EQ.3) THEN
            PI=G_VMC_PI
            PJ=1.D0-PI
            GETPATHPROB2= XIJ(1,2)*(PI*XIJ(1,3)+PJ*XIJ(2,3))/
     &                     (1-PI*XIJ(1,2)-PJ*XIJ(2,1))
     &                  +XIJ(1,3)*(PI*XIJ(1,2)+PJ*XIJ(3,2))/
     &                     (1-PI*XIJ(1,3)-PJ*XIJ(3,1))
C            IPATH(1)=1
C            RET=0.D0
C            CALL GETPP2_R(IPATH,XIJ,M,I_V,2,RET,1.D0,INV)
C            WRITE(61,*) GETPATHPROB2,RET
C         ENDIF
         ELSE
C.. We require a recursive function summing over all possible
C.. permutations of jkl... and calculating appropriate inverse matrix
C.. elements
            IPATH(1)=1
            M(1,1)=1.D0
            RET=0.D0
C            STOP 'Cannot handle new path gen with IV_MAX>3'
            CALL GETPP2_R(IPATH,XIJ,M,I_V,2,RET,1.D0,INV)
            GETPATHPROB2=RET
         ENDIF
         RETURN
      END
   
      REAL*8 FUNCTION GETPATHPROB(XIJ,I_V)
         IMPLICIT NONE
         INTEGER I_V
         REAL*8 XIJ(I_V,I_V)
         INTEGER IPATH(I_V)
         REAL*8 M(I_V,I_V)
         REAL*8 RET
         IF(I_V.EQ.2) THEN
            GETPATHPROB=XIJ(1,2)
         ELSEIF(I_V.EQ.3) THEN
            GETPATHPROB= XIJ(1,2)*(XIJ(2,3)+XIJ(2,1)*XIJ(1,3))/
     &                     (1-XIJ(1,2)*XIJ(2,1))
     &                  +XIJ(1,3)*(XIJ(3,2)+XIJ(3,1)*XIJ(1,2))/
     &                     (1-XIJ(1,3)*XIJ(3,1))
C         ENDIF
         ELSE
C.. We require a recursive function summing over all possible
C.. permutations of jkl... and calculating appropriate inverse matrix
C.. elements
            IPATH(1)=1
            M(1,1)=1.D0
            RET=0.D0
            CALL GETPP_R(IPATH,XIJ,M,I_V,2,RET,1.D0)
            GETPATHPROB=RET
         ENDIF
         RETURN
      END

      RECURSIVE SUBROUTINE GETPP2_R(IPATH,XIJ,M,I_VMAX,I_V,RET,TMP,INV)
         IMPLICIT NONE
         INTEGER I_VMAX,I_V,IPATH(I_VMAX)
         REAL*8 XIJ(I_VMAX,I_VMAX),M(I_VMAX,I_VMAX)
         REAL*8 INV(I_VMAX,I_VMAX)
         REAL*8 WORK(I_VMAX,I_VMAX),RET,TMP
         INTEGER I,J,K,IPIVOT(I_VMAX),INFO,L
         REAL*8 P(I_VMAX)
         LOGICAL T
         INCLUDE 'vmc.inc'
C.. M contains the I-X matrix being constructed
C.. Add this node
         IF(I_V.GT.I_VMAX) THEN
            RET=RET+TMP
            RETURN
         ENDIF
         DO J=2,I_VMAX
            T=.TRUE.
            DO K=2,I_V-1
               IF(IPATH(K).EQ.J) T=.FALSE.
            ENDDO
C.. If we're not already in the path
            IF(T) THEN
               IPATH(I_V)=J
C.. Create the prob matrix
               IF(I_V.GT.2) THEN
                  P(1)=G_VMC_PI
               ELSE
                  P(1)=1.D0
               ENDIF
               DO I=2,I_V-1
                     P(I)=(1.D0-G_VMC_PI)/(I_V-2.D0)
               ENDDO
               P(I_V)=0.D0
C.. Create matrix M 
               DO I=1,I_V
                  DO K=1,I_V
                     INV(I,K)=0.D0
                     DO L=1,I_V
                        IF(L.NE.K) INV(I,K)=INV(I,K)
     &                     -P(L)*XIJ(IPATH(L),IPATH(K))
                     ENDDO
                     IF(I.EQ.I_V) INV(I,K)=0.D0
                     IF(I.EQ.K) THEN
                        INV(I,K)=1.D0+INV(I,K)
                     ENDIF
                  ENDDO
C                  WRITE(61,*) (INV(I,K),K=1,I_V)
               ENDDO
C.. Invert the matrix, and pass the multiplicative factor on.
C.. LU decomp
               CALL DGETRF(I_V,I_V,INV,I_VMAX,IPIVOT,INFO)
C.. Inverse after LU decomp
               CALL DGETRI(I_V,INV,I_VMAX,IPIVOT,WORK,I_VMAX,INFO)
               CALL GETPP2_R(IPATH,XIJ,M,I_VMAX,I_V+1,RET,
     &            TMP*INV(I_V-1,I_V),INV)
            ENDIF
         ENDDO
         RETURN
      END

      RECURSIVE SUBROUTINE GETPP_R(IPATH,XIJ,M,I_VMAX,I_V,RET,TMP)
         IMPLICIT NONE
         INTEGER I_VMAX,I_V,IPATH(I_VMAX)
         REAL*8 XIJ(I_VMAX,I_VMAX),M(I_VMAX,I_VMAX)
         REAL*8 INV(I_VMAX,I_VMAX)
         REAL*8 WORK(I_VMAX,I_VMAX),RET,TMP
         INTEGER I,J,K,IPIVOT(I_VMAX),INFO
         LOGICAL T
C.. M contains the I-X matrix being constructed
C.. Add this node
         IF(I_V.GT.I_VMAX) THEN
            RET=RET+TMP
            RETURN
         ENDIF
         DO J=2,I_VMAX
            T=.TRUE.
            DO K=2,I_V-1
               IF(IPATH(K).EQ.J) T=.FALSE.
            ENDDO
C.. If we're not already in the path
            IF(T) THEN
               IPATH(I_V)=J
               DO I=1,I_V-1
                  M(I,I_V)=-XIJ(IPATH(I),J)
                  M(I_V,I)=-XIJ(J,IPATH(I))
               ENDDO
               M(I_V,I_V)=1.D0
               CALL DCOPY(I_VMAX**2,M,1,INV,1)
               DO I=1,I_V-1
                  INV(I_V,I)=0.D0
               ENDDO
C.. Invert the matrix, and pass the multiplicative factor on.
               CALL DGETRF(I_V,I_V,INV,I_VMAX,IPIVOT,INFO)
               CALL DGETRI(I_V,INV,I_VMAX,IPIVOT,WORK,I_VMAX,INFO)
               CALL GETPP_R(IPATH,XIJ,M,I_VMAX,I_V+1,RET,
     &            TMP*INV(I_V-1,I_V))
            ENDIF
         ENDDO
         RETURN
      END

      SUBROUTINE GENRANDOMEXCITSYM(NI,NEL,NBASIS,G1,IEXLEVEL,ISEED,
     &            NBASISMAX,KSYM,NJ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,IEXLEVEL,ISEED,NJ(NEL)
         INTEGER I,J,K,IEX,IEXL2
         LOGICAL BR,BR2
         INTEGER IGETEXCITLEVEL
         INCLUDE 'basis.inc'
         TYPE(BasisFN) KSYM
         INTEGER G1(*),NBASISMAX(5,2)
         LOGICAL LCHKSYM
         REAL*8 RAN2
         IF(IEXLEVEL.GT.2)
     &    STOP "Cannot handle more than double excitations."
         IF(IEXLEVEL.LT.2) THEN
            STOP "No sym excitations for IEXLEVEL<2"
         ENDIF
         IEXL2=IEXLEVEL
         BR2=.TRUE.
         DO WHILE(BR2)
            CALL ICOPY(NEL,NI,1,NJ,1)
            DO IEX=1,IEXL2
               BR=.TRUE.
C.. Find an electron we haven't excited before
               DO WHILE (BR)
                  I=INT(RAN2(ISEED)*NEL)+1
                  IF(NJ(I).EQ.NI(I)) BR=.FALSE.
               ENDDO
               BR=.TRUE.
               DO WHILE (BR)
                  BR=.FALSE.
                  J=INT(RAN2(ISEED)*NBASIS)+1
                  DO K=1,NEL
C.. If the new basis fn's in our original det or our new one, we loop again
                     IF(NI(K).EQ.J.OR.NJ(K).EQ.J) BR=.TRUE.
                  ENDDO
               ENDDO
               NJ(I)=J
            ENDDO
            CALL SORTI(NEL,NJ)
            BR=LCHKSYM(NI,NJ,NEL,G1,NBASISMAX)
            BR2=.NOT.BR
         ENDDO
         RETURN
      END

      SUBROUTINE GETEXCITCOUNTSYM(INODE,NEL,NBASIS,G1,NBASISMAX,
     &            KSYM,BRR,ICOUNT)
         IMPLICIT NONE
         INTEGER NEL,INODE(NEL),NBASIS,G1(*),NBASISMAX(5,2)
         INTEGER ICOUNT,BRR(NBASIS)
         INCLUDE 'basis.inc'
         TYPE(BasisFN) KSYM
         ICOUNT=0
         CALL GENSYMDETSSD(INODE,KSYM,NEL,G1,BRR,NBASIS,0,
     &        ICOUNT,NBASISMAX,0.D0,0,0,0,0,0,0,0,0,0.D0,0.D0,0.D0,0.D0,
     &        0.D0,0.D0)
         RETURN
      END

      SUBROUTINE GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,NEXCITS,
     &      ISEED,NJ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,NEXCITS,ISEED,NJ(NEL)
         INTEGER I,J,K,IEXL2,NB2,NA,NB,NEX1,NEX2
         INTEGER IEX1,IEX2,IIEX1,IIEX2
         LOGICAL BR
         INTEGER NBASISMAX(5,6)
         INTEGER NEX1A,NEX1B,NEX2A,NEX2B,G1(*)
         INTEGER NEXAB
         REAL*8 R,IEX
         LOGICAL ISUHFDET
         REAL*8 RAN2
         INCLUDE 'basis.inc'
         TYPE(BasisFN) KSym
C         CALL WRITEDET(56,NI,NEL,.TRUE.)
C         CALL FLUSH(56)
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
            CALL ICOPY(NEL,NI,1,NJ,1)
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
C         WRITE(56,*) NEX1A,NEX1B,NEX2A,NEX2B,NEXAB
         NEX1B=NEX1B+NEX1A
         NEX2A=NEX2A+NEX1B
         NEX2B=NEX2B+NEX2A
         NEXAB=NEXAB+NEX2B
         NEX2=0
         R=RAN2(ISEED)
         IEX=R*NEXCITS
         IF(IEX.LT.NEX1A) THEN
C            WRITE(57,"(A2,$)") "A1"
            NEX1=1
         ELSEIF(IEX.LT.NEX1B) THEN
C            WRITE(57,"(A2,$)") "B1"
            NEX1=-1
         ELSEIF(IEX.LT.NEX2A) THEN
C            WRITE(57,"(A2,$)") "A2"
            NEX1=1
            NEX2=1
         ELSEIF(IEX.LT.NEX2B) THEN
C            WRITE(57,"(A2,$)") "B2"
            NEX1=-1
            NEX2=-1
         ELSE
C            WRITE(57,"(A2,$)") "AB"
            NEX1=1
            NEX2=-1
         ENDIF
         IEX1=0
         CALL FINDELECSPIN(NI,NEL,NEX1,G1,ISEED,NBASIS,IEX1)
         IF(NEX2.NE.0) THEN
            IEX2=IEX1
            CALL FINDELECSPIN(NI,NEL,NEX2,G1,ISEED,NBASIS,IEX2)
         ELSE
            IEX2=0
         ENDIF
         IIEX1=0
         IIEX2=0
         CALL FINDNEWELECSPIN(NI,NEL,NEX1,G1,ISEED,NBASIS,IIEX1,IIEX1)
         IF(IEX2.NE.0) THEN
            CALL FINDNEWELECSPIN(NI,NEL,NEX2,G1,ISEED,NBASIS,IIEX2,
     &         IIEX1)
            NJ(IEX2)=IIEX2
         ENDIF
         NJ(IEX1)=IIEX1 
         CALL SORTI(NEL,NJ)
C         CALL WRITEDET(57,NJ,NEL,.TRUE.)
C         CALL FLUSH(57)
         RETURN
      END
      
      SUBROUTINE FINDELECSPIN(NI,NEL,NSPIN,G1,ISEED,NBASIS,IEX)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,NSPIN,IEX,ISEED,IEL
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS
         LOGICAL BR
         REAL*8 RAN2
         BR=.TRUE.
         DO WHILE(BR)
            IEL=RAN2(ISEED)*NEL+1
            IF(G1(NI(IEL))%Ms.EQ.NSPIN.AND.IEL.NE.IEX) BR=.FALSE.
         ENDDO
         IEX=IEL
         RETURN
      END
      SUBROUTINE FINDNEWELECSPIN(NI,NEL,NSPIN,G1,ISEED,NBASIS,IEX,IEX1)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,NSPIN,IEX,ISEED,IEL,IEX1,NBASIS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         LOGICAL BR
         INTEGER I
         REAL*8 RAN2
         BR=.TRUE.
         DO WHILE(BR)
            IEL=RAN2(ISEED)*NBASIS+1
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
      END
      LOGICAL FUNCTION ISUHFDET(NI,NEL)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         INTEGER I
         INCLUDE 'uhfdet.inc'
         ISUHFDET=.FALSE.
         IF(.NOT.TUSEBRILLOUIN) RETURN
            ISUHFDET=.TRUE.
            DO I=1,NEL
               IF(NI(I).NE.NUHFDET(I)) ISUHFDET=.FALSE.
            ENDDO
C         ISUHFDET=.FALSE.
         RETURN
      END


C.. Like MCPATHSR4C, this selects a n-vertex graph.  it chooses which
Cvertex to excite from randomly however.
C.. 
C.. This does a markov chain monte carlo also - the probabilities are 
C.. weighted with the values of S' (the double-counting corrected
C.. weight.
!.. NMAX has ARR hidden in it
      FUNCTION FMCPR4D2(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS,NMEM,OETILDE,OPROB,
     &   I_OVCUR,IOCLS,ITREE,OWEIGHT,PFAC,IACC,INWI,I_VMAX,EXCITGEN)
         USE HElement
         IMPLICIT NONE
         TYPE(HDElement) FMCPR4D2
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V,I_VMAX
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,ALAT(*),ECORE
         COMPLEX*16 FCK
         TYPE(HElement) UMat(*), TMat(*)
         REAL*8 RHOEPS
         TYPE(HDElement) RHOII(0:I_V),INWI
         TYPE(HElement) RHOIJ(0:I_V,0:I_V)
         INTEGER KSYM(5)
         INTEGER LSTE(NEL),NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT
         INTEGER I,J,K,ICE,IC
         REAL*8 XIJ(0:I_V-1,0:I_V-1)         
         TYPE(HDElement) CALCPATHS_N
         REAL*8 X,GETPATHPROB2,RHX
         LOGICAL LISINPATH
         INTEGER IGETEXCITLEVEL,EXCITGEN(0:I_V)
         LOGICAL TLOG,TLOG2,TLOG3
         REAL*8 RP,DBETA(3)
         TYPE(HDElement) DLWDB,DLWDB2
         TYPE(HElement)HIJS(0:I_V)
         LOGICAL ISUHFDET
         REAL*8 PEXCIT(I_V),R
         INTEGER NEXNODE
         INTEGER NMEM(*)
         INTEGER CUREX(*),NEWEX(*),NEWEXLEN
         INTEGER PVERTMEMS(0:I_V)
         POINTER (IP_CUREX,CUREX),(IP_NEWEX,NEWEX)
         INTEGER ICLS
         TYPE(HDElement) ORHOII(0:I_OVCUR)
         TYPE(HElement) ORHOIJ(0:I_OVCUR,0:I_OVCUR)
         INTEGER OIPATH(NEL,0:I_OVCUR)
         TYPE(HElement) OHIJS(0:I_OVCUR)
         REAL*8 OXIJ(0:I_OVCUR-1,0:I_OVCUR-1),OPROB,PR,R2
         TYPE(HDElement) OETILDE,WEIGHT,OWEIGHT,ETILDE
         INTEGER I_OVCUR,IOCLS,ITREE
         LOGICAL LISNAN,ISVALIDDET
         INTEGER IACC
         REAL*8 PFAC
         REAL*8 RAN2
         INTEGER STORE(6)
         TYPE(HDElement) MPEs(2:i_VMax)
         INCLUDE 'vmc.inc'
         
C.. Take a copy of the old path and rho matrix etc.
         IF(OWEIGHT%v.ne.0.D0) THEN
            CALL DCOPY(I_OVCUR+1,RHOII,1,ORHOII,1)
            CALL DCOPY((I_OVCUR+1)**2,RHOIJ,1,ORHOIJ,1)
            CALL DCOPY(I_OVCUR+1,HIJS,1,OHIJS,1)
            CALL DCOPY(I_OVCUR*I_OVCUR,XIJ,1,OXIJ,1)
            CALL ICOPY(NEL*(1+I_OVCUR),IPATH,1,OIPATH,1)
        ENDIF
         IF(.NOT.(INWI.AGT.0.D0)) THEN
!.. Generate a new graph
!  DBETA(3) has PGENEPSILON hidden in it.
            RHX=-1
            IC=0
            TLOG=BTEST(ILOGGING,2)
            TLOG2=BTEST(ILOGGING,3)
            TLOG3=BTEST(ILOGGING,6)
            DO WHILE (RHX.LE.DBETA(3))
               CALL FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,I_V,XIJ,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,
     &           NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,NMEM,I_HMAX,
     &           EXCITGEN)
C.. Now get the parameters for the new path 
            
C.. GETPATHPROB gives us the probability of generating the path
               RHX=GETPATHPROB2(XIJ,I_V)*PFAC
               IC=IC+1
               IF(IC.GT.10000) THEN
                  WRITE(6,*) "Have thrown away 10000 graphs with pGen<",
     &               DBETA(3) 
                  STOP "Have thrown away 10000 graphs with pGen<"
     &               //"pGenEpsilon"
               ENDIF
            ENDDO
            ICLS=0
            IF(TMPTHEORY) THEN
! MP Theory requires HDIAG
               CALL AZZERO(MPEs,I_VMAX-1)
               CALL AddMPEnergy(RHOIJ,i_V,i_vmax,NMAX,nBasis,
     &            iPath,nEl,.false.,ECORE,MPEs)
               ETILDE=0.D0
               DO I=2,I_VMAX
                  ETILDE=ETILDE+MPEs(I)
               ENDDO
               WEIGHT=1.D0
            ELSE
C.. CALCPATHS gives us the contribution of the path
               WEIGHT=CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &            I_P,0.D0,NEL,I_V,ILOGGING,DBETA,ETILDE,HIJS,ICLS)
            ENDIF
         ELSE
C.. we use the gestalt
            RHX=PFAC
            ICLS=1
            WEIGHT=INWI
            ETILDE=DLWDB
         ENDIF
C.. Now do an accept/reject
         IF(.NOT.(OWEIGHT.AGT.0.D0)) THEN
C.. this is the first time round, so we automatically accept the new
C.. configuration
            OETILDE=ETILDE
            OPROB=RHX
            I_OVCUR=I_V
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            IACC=1
         ELSE
C         CALL WRITEPATH(6,IPATH,I_V,NEL,.FALSE.)
C.. Accept the new configuration with prob (|S'(new)|/|S'(old)|)*(OPROB/RH)  

C.. acc(A->B)/acc(B->A) = p(B) gen(B->A) / (p(A) gen(A->B))
C..                     = p(B) gen(A) / (p(A) gen(B))

            PR=(OPROB/RHX)*ABS(DREAL(WEIGHT/OWEIGHT))
!            WRITE(40,"(5G,$)") OPROB,RH,DREAL(WEIGHT), DREAL(OWEIGHT),PR
            IF(RHX.EQ.0.D0.OR..NOT.(OWEIGHT.AGT.0.D0)) PR=0.D0
            IF(.not.TPREVAR) THEN
                IF(LISNAN(WEIGHT)) THEN
                    WRITE(60,*) WEIGHT,ETILDE,RHX
                    CALL WRITEPATH(60,IPATH,I_V,NEL,.FALSE.)
                    CALL WRITERHOMAT(60,RHOIJ,I_V,NEL,.TRUE.)
                    CALL WRITE_XMATRIX(60,XIJ,I_V)
                ENDIF
                R2=RAN2(ISEED)
                IF(PR.LT.R2) THEN
                    IACC=0
                ELSE
                    IACC=1
                ENDIF
            ELSE
                IACC=1
            ENDIF
!            CALL WRITEPATHEX(40,OIPATH,I_OVCUR,NEL,.FALSE.)
!            WRITE(40,"(A,$)"),"--->"
!            CALL WRITEPATHEX(40,IPATH,I_V,NEL,.FALSE.)
!            WRITE(40,*) IACC
         ENDIF
         CALL CLASSPATHS(FMCPR4D2,DLWDB,1.D0,RHOIJ,
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
            IF(BTEST(ILOGGING,3))
     &            CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            IF(BTEST(ILOGGING,6)) CALL WRITE_XMATRIX(10,XIJ,I_V)
            WRITE(10,"(3E25.16, 2I7)") WEIGHT,RHX,ETILDE,ICLS,IACC
         ENDIF
          
         IF(IACC.EQ.0) THEN
C.. reject the new, and copy in the old
            ETILDE=OETILDE
            RHX=OPROB
            I_V=I_OVCUR
            ICLS=IOCLS
            WEIGHT=OWEIGHT
            CALL ICOPY(NEL*(1+I_OVCUR),OIPATH,1,IPATH,1)
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
            DLWDB=ETILDE/HDElement(OPROB)
            FMCPR4D2=OPROB
         ELSEIF(tMCDirectSum.AND.RHX.GT.0.D0.AND.(WEIGHT.AGT.0.D0)) THEN
C.. Unbias the sum 
            DLWDB=ETILDE/WEIGHT
            FMCPR4D2=WEIGHT/WEIGHT
         ELSEIF(RHX.GT.0.D0.AND.(WEIGHT.AGT.0.D0)) THEN
C.. Unbias the sum 
            DLWDB=ETILDE/HDElement(ABS(DREAL(WEIGHT)))
            FMCPR4D2=WEIGHT/HDElement(ABS(DREAL(WEIGHT)))
         ELSE
            FMCPR4D2=0.D0
            DLWDB=0.D0
         ENDIF
         CALL GETTREENESS(ICLS,ITREE,FMCPR4D2,I_V)
      END











!.. Generate a new graph.  Called from FMCPR4D2
!.. NMAX has ARR  hidden in it
      SUBROUTINE FMCPR4D2GENGRAPH(NI,NEL,BETA,I_P,IPATH,I_V,XIJ,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMat,
     &   NTAY,RHOEPS,RHOII,RHOIJ,ECORE,ISEED,HIJS,NMEMNI,I_HMAX,
     &   PVERTMEMS2)
         USE HElement
         IMPLICIT NONE
         INCLUDE 'vmc.inc'
         INCLUDE 'basis.inc'
         INTEGER NEL,I_V,IPATH(NEL,0:I_V),NI(NEL)
         REAL*8 BETA,ALAT(*),ECORE
         TYPE(HElement) TMat(*),Umat(*)
         TYPE(BasisFN) G1(*)
         INTEGER I_P,NBASISMAX(*),NBASIS,BRR(*),NMSH
         INTEGER NMAX,NTAY
         INTEGER PVERTMEMS2(0:I_V)
         COMPLEX*16 FCK(*)
         INTEGER I_HMAX

         REAL*8 RHOEPS
         TYPE(HDElement) RHOII(0:I_V)
         TYPE(HElement) RHOIJ(0:I_V,0:I_V)
         REAL*8 XIJ(0:I_V-1,0:I_V-1)         
         TYPE(HElement) HIJS(0:I_V)
         INTEGER ISEED

         REAL*8 PEXCIT(I_V)
         INTEGER PVERTMEMS(0:I_V)
         INTEGER I_VNEXT
         INTEGER INODE(NEL),INODE2(NEL)
         INTEGER NMEMNI(*)
         INTEGER CUREX(*),NEWEX(*),NEWEXLEN
         INTEGER NEXNODE
         POINTER (IP_CUREX,CUREX),(IP_NEWEX,NEWEX)
         INTEGER STORE(6)
         REAL*8 R
         TYPE(HElement) RH
         INTEGER ICE,I,ICOUNT,IC

         INTEGER IGETEXCITLEVEL
         LOGICAL ISUHFDET,LISINPATH
         LOGICAL LISNAN,ISVALIDDET
         REAL*8 RAN2
         REAL*8 pGen,pGen2

         !Deallocate Excitation Generators after precalc
         IF(I_HMAX.EQ.154) THEN
             DO I=1,I_V-1
                 CALL FREEM(PVERTMEMS2(I))
             ENDDO
             RETURN
!          STOP 'I_HMAX of 154 should not be called here (mcpathsismc.F)'
          ENDIF
C.. PEXCIT(NODE) is the probability of selecting NODE and
C.. those before it.
         IF(I_HMAX.EQ.-19) THEN
!.. We're actually storing H matrices in RHOIJ rather than rho matrices
            RHOIJ(0,0)=GETHELEMENT2(NI,NI,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,
     &                  0,ECORE)
            RHOII(1)=BETA
         ENDIF
         PEXCIT(1)=1.D0
         !LOC(NMEMNI) is a pointer to NMEM
         PVERTMEMS(0)=LOC(NMEMNI)
         I_VNEXT=1
         CALL ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0

C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
C.. doesn't bother calculating the RHO_JJ.

!24/04/07 - ghb24 - Removing (hopefully) redundant routine calls - now
!       that there is a separate weighting routine, it is not needed...
!          CALL GENRANDSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
!     &         NMEMNI,INODE2,ISEED,IC,0,UMAT,NMAX,pGen)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         ICOUNT=0
         DO WHILE(I_VNEXT.LT.I_V)
C.. pick a random excitation of where we are (INODE) of the appropriate
C.. symmetry.  We don't worry about weighting this at the moment
C.. We pick which node we'regoing to excite from
            IF(I_VNEXT.GT.1) THEN
               DO I=2,I_VNEXT
                  PEXCIT(I)=(1.D0-G_VMC_PI)/(I_VNEXT-1.D0)
               ENDDO
               PEXCIT(1)=G_VMC_PI
            ENDIF
C.. Just make really sure the last one catches everytihng
            PEXCIT(I_VNEXT)=2.D0
            R=RAN2(ISEED)
            NEXNODE=0
            DO WHILE(R.GE.0)
               NEXNODE=NEXNODE+1
               R=R-PEXCIT(NEXNODE)
            ENDDO
! Vertex labels start at 0
            NEXNODE=NEXNODE-1
C.. set this chosen node to be the pivot
            IP_CUREX=PVERTMEMS(NEXNODE)

C.. generate a random excitation
            CALL GENRANDSYMEXCITIT2(IPATH(1,NEXNODE),NEL,G1,NBASIS,
     &            NBASISMAX,CUREX,INODE2,ISEED,IC,0,UMAT,NMAX,pGen)


            CALL ICOPY(NEL,INODE2,1,INODE,1)

C.. If the new node is not in the path, add it.
            IF(.NOT.ISVALIDDET(INODE,NEL)) THEN
               WRITE(6,*) "INVALID DET"
               CALL WRITEDET(6,INODE,NEL,.TRUE.)
               STOP "INVALID DET"
            ENDIF
            IF(.NOT.LISINPATH(INODE,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL ICOPY(NEL,INODE,1,IPATH(1,I_VNEXT),1)

C.. Setup the spin excit generator
               STORE(1)=0
               !determine how much memory needed
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            .TRUE.,NEWEXLEN,INODE2,IC,0,STORE,3)
               !allocate memory
               CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
               NEWEX(1)=0
               PVERTMEMS(I_VNEXT)=IP_NEWEX
               !Generate Excitation generator
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            .TRUE.,NEWEX,INODE2,IC,0,STORE,3)
               !Save the pointers to the excitation generator
               IF(TPREVAR) THEN
                    PVERTMEMS2=PVERTMEMS
               ENDIF
C.. Count the excitations (and generate a random one which we throw)
 
!24/04/07 - ghb24 - Removing (hopefully) redundant routine calls - now
!       that there is a separate weighting routine, it is not needed..
!               CALL GENRANDSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
!     &            NEWEX,INODE2,ISEED,IC,0,UMAT,NMAX,pGen2)

               XIJ(I_VNEXT,I_VNEXT)=IC
C.. Update the rho and X (probability) matrices with this new node
               IF(I_HMAX.EQ.-19) THEN
!  using H-elements in the rho matrix
                 RH=GETHELEMENT2(INODE,INODE,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,
     &                  UMAT,0,ECORE)
                  RHOIJ(I_VNEXT,I_VNEXT)=RH
               ELSE
                  CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,TMat,TMat,NMAX,ALAT,UMAT,
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
                        IP_CUREX=PVERTMEMS(I)
!  NMAX is really Arr
                        CALL GenExcitProb(IPATH(1,I),INODE,nEl,
     &                   CUREX,G1,nBasisMax,UMat,NMAX,nBasis,pGen2)
                     ENDIF
                     XIJ(I,I_VNEXT)=pGen2
!  Now work out the reverse prob of generation.
!  NMAX is really Arr
                     CALL GenExcitProb(INODE,IPATH(1,I),nEl,
     &                NEWEX,G1,nBasisMax,UMat,NMAX,nBasis,pGen2)
                     XIJ(I_VNEXT,I)=pGen2
      
! We supercede the old unbiased random choice with a new one.
!                     XIJ(I_VNEXT,I)=1.D0/IC
!                     IF(ISUHFDET(IPATH(1,I),NEL,NBASISMAX)) THEN
!                        IF(ICE.EQ.2) THEN
!                           XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
!                        ELSE
!                           XIJ(I,I_VNEXT)=0.D0
!                        ENDIF
!                     ELSE
!                        XIJ(I,I_VNEXT)=1.D0/XIJ(I,I)
!                     ENDIF   
                  ELSE
                     XIJ(I_VNEXT,I)=0.D0
                     XIJ(I,I_VNEXT)=0.D0
                  ENDIF
                  IF(I_HMAX.EQ.-19) THEN
!  using H-elements in the rho matrix
                     RH=GETHELEMENT2(INODE,IPATH(1,I),NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,
     &                  UMAT,ICE,ECORE)
                     IF(I.EQ.0) HIJS(I_VNEXT)=RH
                  ELSE
                     IF(I.EQ.0) THEN
                      HIJS(I_VNEXT)=GETHELEMENT2(NI,INODE,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,
     &                  ICE,ECORE)
                     ENDIF
                     CALL CALCRHO2(INODE,IPATH(1,I),BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,
     &             RH,NTAY,ICE,ECORE)
                  ENDIF
                  IF(RH.AGT.RHOEPS) THEN
                     RHOIJ(I_VNEXT,I)=RH
                     RHOIJ(I,I_VNEXT)=RH
                  ELSE
                     RHOIJ(I,I_VNEXT)=0.D0
                     RHOIJ(I_VNEXT,I)=0.D0
                  ENDIF
               ENDDO
               I_VNEXT=I_VNEXT+1
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.500) THEN
C.. give up trying to find something to attach, and go home
               WRITE(6,*) "WARNING: Unable to find attachee to vertex"
               CALL WRITEDET(6,IPATH(1,NEXNODE-1),NEL,.TRUE.)
               RHOII(1)=0.D0
               RETURN
            ENDIF
         ENDDO
         !Don't deallocate if want to use the graphs later for MC precalc
         IF(.NOT.TPREVAR) THEN
            DO I=1,I_V-1
                CALL FREEM(PVERTMEMS(I))
            ENDDO
         ENDIF
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)
      END


      SUBROUTINE GETTREENESS(ICLASS,ITREE,WEIGHT,I_V)
         IMPLICIT NONE
         INTEGER ICLASS,ITREE,I_V,N,ICL
         REAL*8 WEIGHT
         ITREE=0
         IF(WEIGHT.EQ.0.D0) RETURN
         N=0
         ICL=ICLASS
         DO WHILE (ICL.NE.0)
            IF(IAND(ICL,1).EQ.1) N=N+1
            ICL=ICL/2
         ENDDO
         IF(N.EQ.I_V-1) ITREE=1
         RETURN
      END


C.. 11/7/05
C.. Like MCPATHSR4C, this selects a n-vertex graph.  it chooses which
C.. vertex to excite from randomly however, and it throws out a few
C.. tentacles to find excitations, and weights them according to their
C.. primitive probability, building up the graph this way
C.. 
C.. This does a markov chain monte carlo also - the probabilities are 
C.. weighted with the values of S' (the double-counting corrected
C.. weight.

C.. 14/7/05 - This does not work - in order to calculate the 
C.. generation probability successfully, one needs to sample all
C.. possible ways of having thrown out tentacles.  This is not an 
C.. easily computable quantity.  What is used here is merely the 
C.. generation probability using this set of tentacles, which will give
C.. a wrongly biased MC
      REAL*8 FUNCTION FMCPR4D3(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS,NMEM,OETILDE,OPROB,
     &   I_OVCUR,IOCLS,ITREE,OWEIGHT,IACC)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         REAL*8 RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER KSYM(5)
         INTEGER LSTE(NEL),NLIST
         INTEGER INODE(NEL),I_VNEXT,INODE2(NEL),ICOUNT
         INTEGER I,J,K,ICE,IC
         REAL*8 XIJ(0:I_V-1,0:I_V-1)         
         REAL*8 RH,CALCPATHS_N,X,GETPATHPROB2
         LOGICAL LISINPATH
         INTEGER IGETEXCITLEVEL
         LOGICAL TLOG,TLOG2,TLOG3
         REAL*8 RP,DBETA,DLWDB,DLWDB2
         REAL*8 GETHELEMENT2,HIJS(0:I_V)
         LOGICAL ISUHFDET
         REAL*8 PEXCIT(I_V),R,ETILDE
         INTEGER NEXNODE
         INTEGER NMEM(*)
         INTEGER CUREX(*),NEWEX(*),NEWEXLEN
         INTEGER PVERTMEMS(0:I_V)
         POINTER (IP_CUREX,CUREX),(IP_NEWEX,NEWEX)
         INTEGER ICLS
         REAL*8 ORHOII(0:I_V),ORHOIJ(0:I_V,0:I_V)
         INTEGER OIPATH(NEL,0:I_V)
         REAL*8 OHIJS(0:I_V),OXIJ(0:I_V-1,0:I_V-1)
         REAL*8 OETILDE,OPROB,PR,R2,WEIGHT,OWEIGHT
         INTEGER I_OVCUR,IOCLS,ITREE
         LOGICAL LISNAN,ISVALIDDET
         INTEGER IACC
         REAL*8 NRHOII(0:I_V),NRHOIJ(0:I_V,0:I_V)
         INTEGER NIPATH(NEL,0:I_V),IACCNODE
         REAL*8 NHIJS(0:I_V),NXIJ(0:I_V-1,0:I_V-1)
         REAL*8 CURWEIGHT,WEIGHTTOT,PGEN,CURDLWDB,CURICLS
         INTEGER IEXNUM,NTENT
         REAL*8 RAN2
         INTEGER STORE(6)
         INCLUDE 'vmc.inc'
C         WRITE(6,*) "START:",LOC(NMEM)

C.. number of tentacles
         NTENT=4
C.. PEXCIT(NODE) is the probability of selecting NODE and
C.. those before it.
         PEXCIT(1)=1.D0
         PVERTMEMS(0)=LOC(NMEM)
         RP=0
         I_VNEXT=1
         CALL ICOPY(NEL,IPATH(1,0),1,INODE,1)
         IC=0
C.. Take a copy of the old path and rho matrix etc.
         CALL DCOPY(I_V+1,RHOII,1,ORHOII,1)
         CALL DCOPY((I_V+1)**2,RHOIJ,1,ORHOIJ,1)
         CALL DCOPY(I_V+1,HIJS,1,OHIJS,1)
         CALL DCOPY(I_V*I_V,XIJ,1,OXIJ,1)
         CALL ICOPY(NEL*(1+I_V),IPATH,1,OIPATH,1)













C.. Get the excitation count
            CALL GENRANDSYMEXCITIT2(NI,NEL,G1,NBASIS,
     &         NBASISMAX,NMEM,INODE2,ISEED,IC,0,UMAT,NMAX,PGEN)




CC.. setup excitation generator
C         CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
C     &         .FALSE.,NMEM,INODE2,ISEED,IC,0)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         NLIST=1
         ICOUNT=0
         PGEN=1.D0
         IF(I_V.EQ.1) THEN
            CURWEIGHT=1.D0
            CURDLWDB=HIJS(0)
            CURICLS=0
         ENDIF
         DO WHILE(I_VNEXT.LT.I_V)
C.. pick a random excitation of where we are (INODE) of the appropriate
C.. symmetry.  We don't worry about weighting this at the moment
C.. We pick which node we'regoing to excite from
            IF(I_VNEXT.GT.1) THEN
               DO I=2,I_VNEXT
                  PEXCIT(I)=(1.D0-G_VMC_PI)/(I_VNEXT-1.D0)
               ENDDO
               PEXCIT(1)=G_VMC_PI
            ENDIF
C.. Just make really sure the last one catches everytihng
            PEXCIT(I_VNEXT)=2.D0
            R=RAN2(ISEED)
            NEXNODE=0
            DO WHILE(R.GE.0)
               NEXNODE=NEXNODE+1
               R=R-PEXCIT(NEXNODE)
            ENDDO

C.. NEXNODE-1 is the excitation node

C.. set the current to be the pivot
            IP_CUREX=PVERTMEMS(NEXNODE-1)

C.. we now generate a number of different excitations, accepting and
C.. rejecting them on the basis of their weights.
            CURWEIGHT=0.D0
            WEIGHTTOT=0.D0
            CALL DCOPY(I_V+1,RHOII,1,NRHOII,1)
            CALL DCOPY((I_V+1)**2,RHOIJ,1,NRHOIJ,1)
            CALL DCOPY(I_V+1,HIJS,1,NHIJS,1)
            CALL DCOPY(I_V*I_V,XIJ,1,NXIJ,1)
            CALL ICOPY(NEL*(1+I_V),IPATH,1,NIPATH,1)
C.. Hard code at 4 excitations.
            IEXNUM=0
            DO WHILE (IEXNUM.LT.NTENT)

C.. generate a random excitation
            CALL GENRANDSYMEXCITIT2(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
     &         NBASISMAX,CUREX,INODE2,ISEED,IC,0,UMAT,NMAX,PGEN)


CC.. Generate an excitation
C               CALL GENRANDSYMEXCITIT(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
C     &            NBASISMAX,.FALSE.,CUREX,INODE2,ISEED,IC,0)
               CALL ICOPY(NEL,INODE2,1,INODE,1)
               IF(.NOT.ISVALIDDET(INODE,NEL)) THEN
                  WRITE(6,*) "INVALID DET"
                  CALL WRITEDET(6,INODE,NEL,.TRUE.)
                  STOP "INVALID DET"
               ENDIF

C.. If the new node is not in the path, add it.
               IF(.NOT.LISINPATH(INODE,IPATH,NEL,I_VNEXT,-1)) THEN
                  IEXNUM=IEXNUM+1
                  CALL ICOPY(NEL,INODE,1,NIPATH(1,I_VNEXT),1)
C.. Update the rho and X (probability) matrices with this new node
                  CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,
     &              NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &              RH,NTAY,0,ECORE)
                  NRHOII(I_VNEXT)=RH
                  NRHOIJ(I_VNEXT,I_VNEXT)=RH
C.. We don't know how many connections there are from the new node, but
C.. we will later.  Just put a dummy value in for now
                  NXIJ(I_VNEXT,I_VNEXT)=1.D0
                  DO I=0,I_VNEXT-1
                     ICE=IGETEXCITLEVEL(INODE,NIPATH(1,I),NEL)
                     IF(ICE.LE.2) THEN
C.. We don't know how many connections there are from the new node, but
C.. we will later.  Just put a dummy value in for now
                        NXIJ(I_VNEXT,I)=1.D0
                        IF(ISUHFDET(NIPATH(1,I),NEL)) THEN
                           IF(ICE.EQ.2) THEN
                              NXIJ(I,I_VNEXT)=1.D0/NXIJ(I,I)
                           ELSE
                              NXIJ(I,I_VNEXT)=0.D0
                           ENDIF
                        ELSE
                           NXIJ(I,I_VNEXT)=1.D0/NXIJ(I,I)
                        ENDIF   
                     ELSE
                        NXIJ(I_VNEXT,I)=0.D0
                        NXIJ(I,I_VNEXT)=0.D0
                     ENDIF
                     IF(I.EQ.0) THEN
                        NHIJS(I_VNEXT)=GETHELEMENT2(NI,INODE,NEL,
     &                  NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &                  ICE,ECORE)
                     ENDIF
                     CALL CALCRHO2(INODE,NIPATH(1,I),BETA,I_P,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,
     &                  UMAT,RH,NTAY,ICE,ECORE)
                     IF(ABS(RH).GT.RHOEPS) THEN
                        NRHOIJ(I_VNEXT,I)=RH
                        NRHOIJ(I,I_VNEXT)=RH
                     ELSE
                        NRHOIJ(I,I_VNEXT)=0.D0
                        NRHOIJ(I_VNEXT,I)=0.D0
                     ENDIF
                  ENDDO
C.. by giving CALCPATHS_N 0.D0 as the weighting, it doesn't call
C.. classpaths, or make use of ICLS - we do this later
                  WEIGHT=CALCPATHS_N(NIPATH,NRHOII,NRHOIJ,I_V,I_HMAX,
     &               I_P,0.D0,NEL,I_V,ILOGGING,DBETA,ETILDE,NHIJS,ICLS)
                  IF(CURWEIGHT.EQ.0.D0.OR.
     &               ABS(WEIGHT/CURWEIGHT).GT.RAN2(ISEED)) THEN
C.. Accept the new node
C.. Take a copy of the old path and rho matrix etc.
                     CALL DCOPY(I_V+1,NRHOII,1,RHOII,1)
                     CALL DCOPY((I_V+1)**2,NRHOIJ,1,RHOIJ,1)
                     CALL DCOPY(I_V+1,NHIJS,1,HIJS,1)
                     CALL DCOPY(I_V*I_V,NXIJ,1,XIJ,1)
                     CALL ICOPY(NEL*(1+I_V),NIPATH,1,IPATH,1)
                     CURWEIGHT=WEIGHT
                     CURDLWDB=DLWDB
                  ELSE
C.. keep the old node
                  ENDIF
                  WEIGHTTOT=WEIGHTTOT+ABS(WEIGHT)
               ENDIF
               ICOUNT=ICOUNT+1
               IF(ICOUNT.GT.500) THEN
C.. give up trying to find something to attach, and go home
                  FMCPR4D3=0.D0
                  RETURN
               ENDIF
            ENDDO



C.. Create an excitation generator for this node
C.. Setup the spin excit generator
               STORE(1)=0
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            .TRUE.,NEWEXLEN,INODE2,IC,0,STORE,3)
               CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
               NEWEX(1)=0
               PVERTMEMS(I_VNEXT)=IP_NEWEX
               CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            .TRUE.,NEWEX,INODE2,IC,0,STORE,3)
C.. Count the excitations (and generate a random one which we throw)
               CALL GENRANDSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,
     &            NEWEX,INODE2,ISEED,IC,0,UMAT,NMAX,PGEN)







CC.. Create an excitation generator for this node
C            CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
C     &               .TRUE.,NEWEXLEN,INODE2,ISEED,IC,0)
C            CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
C            NEWEX(1)=0
C            PVERTMEMS(I_VNEXT)=IP_NEWEX
CC.. Count the number of connections to this node
C            CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
C     &               .FALSE.,NEWEX,INODE2,ISEED,IC,0)

            XIJ(I_VNEXT,I_VNEXT)=IC
            DO I=0,I_VNEXT-1
               NXIJ(I_VNEXT,I)=1.D0/IC
            ENDDO
            I_VNEXT=I_VNEXT+1
            PGEN=PGEN*CURWEIGHT*NTENT/(WEIGHTTOT*IC)
         ENDDO
         DO I=1,I_V-1
            CALL FREEM(PVERTMEMS(I))
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)

C.. Now get the parameters for the new path 

C.. GETPATHPROB gives us the probability of generating the path
C         RH=GETPATHPROB2(XIJ,I_V)
         RH=PGEN
C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         WEIGHT=CURWEIGHT
         ETILDE=CURDLWDB
         ICLS=CURICLS
         IF(OWEIGHT.EQ.0.D0) THEN
C.. this is the first time round, so we automatically accept the new
C.. configuration
            OETILDE=ETILDE
            OPROB=RH
            I_OVCUR=I_V
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            IACC=1
         ELSE
C.. Accept the new configuration with prob (|S'(new)|/|S'(old)|)*(OPROB/RH)  

C.. acc(A->B)/acc(B->A) = p(B) gen(B->A) / (p(A) gen(A->B))
C..                     = p(B) gen(A) / (p(A) gen(B))

            PR=(OPROB/RH)*ABS(WEIGHT/OWEIGHT)
            IF(RH.EQ.0.D0.OR.WEIGHT.EQ.0.D0) PR=0.D0
            R2=RAN2(ISEED)
            IF(LISNAN(WEIGHT)) THEN
               WRITE(60,*) WEIGHT,ETILDE,RH
               CALL WRITEPATH(60,IPATH,I_V,NEL,.FALSE.)
               CALL WRITERHOMAT(60,RHOIJ,I_V,NEL,.TRUE.)
               CALL WRITE_XMATRIX(60,XIJ,I_V)
            ENDIF
            IF(PR.LT.R2) THEN
C.. reject the new, and copy in the old
               ETILDE=OETILDE
               RH=OPROB
               I_V=I_OVCUR
               ICLS=IOCLS
               WEIGHT=OWEIGHT
               IACC=0
            ELSE
               OETILDE=ETILDE
               OPROB=RH
               I_OVCUR=I_V
               ICLS=0
               IOCLS=ICLS
               OWEIGHT=WEIGHT
               IACC=1
            ENDIF
         ENDIF
         IF(RH.GT.0.D0.AND.FMCPR4D3.NE.0.D0) THEN
C.. Unbias the sum 
            DLWDB=ETILDE/(ABS(WEIGHT))
            FMCPR4D3=WEIGHT/(ABS(WEIGHT))
         ELSE
            FMCPR4D3=0.D0
            DLWDB=0.D0
         ENDIF
C.. If we've a new graph, IOCLS=0, and classpaths regenerates the class.
C.. otherwise, we just reuse the old version
         CALL CLASSPATHS(FMCPR4D3,DLWDB,1.D0,RHOIJ,
     &            I_V,IOCLS)
         CALL GETTREENESS(IOCLS,ITREE,FMCPR4D3,I_V)
          IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
         IF(TLOG) WRITE(10,"(3E25.16, 2I7)") WEIGHT,RH,ETILDE,ICLS,IACC
         RETURN
      END

C.. 12/7/05
C.. Markov Chain MC.  Generate a new graph each time from the old one by
C.. either
C.. a) adding or removing a vertex
C.. b) Moving one vertex.  THIS SUBROUTINE


C.. To move a vertex in the graph, we first select the vertex to move, k.
C.. This is done according to the prob 
C..   Psel(k|G)=RHO_kk^-P/Sum_(i' in G) Rho_i'i' ^ -P

C.. Next, we throw out NTENT tentacles to excitations k' from the graph
C.. having removed k, each weighted with w(k') = rho_k'k' ^ P
C.. Once G' has been created (from G -k +k'), we work out its primitive
C.. weight w'(G'), and accept G'  with prob

C.. Pacc(G'|G)=min(1,|w'(G')| Psel(k'|G') w(k)/ (|w(G)| Psel(k|G) w(k')))

C.. VWEIGHTS(1,I_V) holds the sum(rhoii^P) of the graph.  Lower elements
C.. have the rhoii^P for the vertices. (excluding the root (position 0))
C.. and VWEIGHTS(2,I_V) holds sum(rhoii^-P) etc.

C.. NB TMAT has sometimes been hidden in ZIA
      SUBROUTINE MCPR4D4MOVE(NI,BETA,I_P,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,NWHTAY,I_HMAX,ILOGGING,
     &   ECORE,ISEED,KSYM,DBETA,ETILDERET,IPATH,RHOIJ,HIJS,WEIGHTRET,
     &   OETILDE,IOCLS,ITREE,OWEIGHT,IACC,VWEIGHTS,VEXCITS,PFAC,TST,
     &   WMIN)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,KSYM(*)
         REAL*8 BETA,DBETA
         INTEGER NBASISMAX(*),G1(*),NBASIS,BRR(*),NMSH,NMAX
         REAL*8 FCK(*),TMAT(*),ALAT(*),UMAT(*),ECORE
         INTEGER NWHTAY,NTAY,I_HMAX,ILOGGING
         REAL*8 RHOEPS
         INTEGER ISEED,I_V
         INTEGER IPATH(NEL,0:I_V)
         REAL*8 RHOIJ(0:I_V,0:I_V)
         REAL*8 HIJS(0:I_V)
         REAL*8 ETILDE,WEIGHT
         REAL*8 OETILDE,OPROB
         INTEGER IOCLS,ITREE,IACC
         REAL*8 VWEIGHTS(2,0:I_V)
         REAL*8 R,RH,PR
         INTEGER IVERT,IVERTMOVE,ITENTNEW,ICOUNT,ITENT
         PARAMETER NTENT=1
         INTEGER TENTDETS(NEL,NTENT)
         REAL*8 WTENTS(NTENT),WTTOT,WT1,WT2
         POINTER (IP_CUREX,CUREX)
         INTEGER CUREX(*),ICE,CUREXLEN
         INTEGER VEXCITS(0:I_V)
         REAL*8 NRHOIJ(0:I_V,0:I_V)
         REAL*8 NHIJS(0:I_V)
         INTEGER NIPATH(NEL,0:I_V)
         REAL*8 OWEIGHT
         INTEGER ICONN,I,ICLS,IC
         REAL*8 ETILDERET,WEIGHTRET,PFAC,WMIN

         LOGICAL ISVALIDDET,ISUHFDET,LISINPATH,AREDETSEXCITS
         INTEGER IGETEXCITLEVEL
         REAL*8 GETHELEMENT2,CALCPATHS_N

         INTEGER TST
         REAL*8 RAN2
         IF(I_V.EQ.1) THEN
C.. Single vertex graph always rejected
            WEIGHT=OWEIGHT
            ETILDE=OETILDE
            IACC=0
            PR=0.D0
            CALL ICOPY(NEL*2,IPATH,1,NIPATH,1)
        ELSE
C.. First pick a vertex to move according to prob in VWEIGHTS(2)
C.. (excluding vertex 0, the pivot)
         R=RAN2(ISEED)*(VWEIGHTS(2,I_V)-VWEIGHTS(2,0))
         IVERT=1
         DO WHILE(VWEIGHTS(2,IVERT).LT.R)
            R=R-VWEIGHTS(2,IVERT)
            IVERT=IVERT+1
         ENDDO
         IF(IVERT.GE.I_V) IVERT=I_V-1
         IVERTMOVE=IVERT
         ICOUNT=0
         ITENT=0
         WTTOT=0.D0
C.. Now throw out some tentacles
         DO WHILE(ITENT.LT.NTENT)
C.. Pick a vertex randomly in the graph to throw from
            IVERT=RAN2(ISEED)*I_V
            IF(IVERT.NE.IVERTMOVE) THEN
C.. generate an excitation from this vertex
               ITENT=ITENT+1
               IP_CUREX=VEXCITS(IVERT)
               CALL GENRANDSYMEXCITIT(IPATH(1,IVERT),NEL,G1,NBASIS,
     &           NBASISMAX,.FALSE.,CUREX,TENTDETS(1,ITENT),ISEED,ICE,0)
               IF(.NOT.ISVALIDDET(TENTDETS(1,ITENT),NEL)) THEN
                  WRITE(6,*) "INVALID DET"
                  CALL WRITEDET(6,TENTDETS(1,ITENT),NEL,.TRUE.)
                  STOP "INVALID DET"
               ENDIF
C.. We cannot deal with the case that the tentacle is already in the graph
C.. later as this mucks up generation probabilities - so we deal with it
C.. now
               IF(.NOT.LISINPATH(TENTDETS(1,ITENT),IPATH,
     &            NEL,I_V,IVERTMOVE)) THEN
                  CALL CALCRHO2(TENTDETS(1,ITENT),TENTDETS(1,ITENT),
     &                  BETA,I_P,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,
     &                  UMAT,RH,NTAY,0,ECORE)
C                  WRITE(6,*)
C     &               (RH/RHOIJ(0,0)),(RH/RHOIJ(0,0))**(I_P*PFAC),WTTOT
                  RH=(RH/RHOIJ(0,0))**(I_P*PFAC)
                  WTENTS(ITENT)=RH
                  WTTOT=WTTOT+RH
               ELSE
                  ITENT=ITENT-1
               ENDIF
            ENDIF
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.NTENT*100) THEN
C.. we've tried 10 times as many tentacles as we need, and we can't find
Cenough.  Something's wrong
               STOP "Not enough tentacles found"
            ENDIF
         ENDDO
C.. Now pick a tentacle vertex to move according to prob in VWEIGHTS(2)
         R=RAN2(ISEED)*WTTOT
         IVERT=1
         DO WHILE(WTENTS(IVERT).LT.R)
            R=R-WTENTS(IVERT)
            IVERT=IVERT+1
         ENDDO
         IF(IVERT.GT.NTENT) IVERT=NTENT
         ITENTNEW=IVERT
C.. Now setup the Rho and HIJ matrix for the new graph.
C.. First copy the old one         
C.. Take a copy of the old path and rho matrix etc.
         CALL DCOPY((I_V+1)**2,RHOIJ,1,NRHOIJ,1)
         CALL DCOPY(I_V+1,HIJS,1,NHIJS,1)
         CALL ICOPY(NEL*(1+I_V),IPATH,1,NIPATH,1)
         CALL ICOPY(NEL,TENTDETS(1,ITENTNEW),1,NIPATH(1,IVERTMOVE),1)
C.. Fill out the parts which have changed
         ICONN=0
         IACC=1
         DO I=0,I_V-1
            IF(I.EQ.IVERTMOVE) THEN
               ICE=-1
               NHIJS(I)=GETHELEMENT2(NI,NIPATH(1,I),NEL,
     &                  NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,UMAT,
     &                  ICE,ECORE)
               ICE=0
               CALL CALCRHO2(NIPATH(1,I),NIPATH(1,I),BETA,
     &                  I_P,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,
     &                  UMAT,RH,NTAY,0,ECORE)
               NRHOIJ(IVERTMOVE,I)=RH
            ELSE
               IF(AREDETSEXCITS(NIPATH(1,IVERTMOVE),
     &                       NIPATH(1,I),NEL,NBASISMAX,ICE)) THEN
                  ICONN=ICONN+1
                  CALL CALCRHO2(NIPATH(1,IVERTMOVE),NIPATH(1,I),BETA,
     &                  I_P,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,
     &                  UMAT,RH,NTAY,ICE,ECORE)
                  IF(ABS(RH).GT.RHOEPS) THEN
                     NRHOIJ(IVERTMOVE,I)=RH
                     NRHOIJ(I,IVERTMOVE)=RH
                  ELSE
                     NRHOIJ(I,IVERTMOVE)=0.D0
                     NRHOIJ(IVERTMOVE,I)=0.D0
                  ENDIF
                  IF(I.EQ.0) THEN
C.. fill out the last rows too
                     NRHOIJ(I_V,IVERTMOVE)=NRHOIJ(0,IVERTMOVE)
                     NRHOIJ(IVERTMOVE,I_V)=NRHOIJ(0,IVERTMOVE)
                  ENDIF
                  IF(ICE.EQ.0) THEN
C.. We've picked a tentacle which is already in the graph - a disaster,
C.. creating a diallowed graph. We set to reject automatically

C.. We're not allowed to do this, as it mucks up the generation
C.. probability for graphs. Must reject the tentacles earlier

C.. If we've reached here, it's an error therefore
                     STOP "Tentacle already in Path"
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
C.. Now we have all the information we need to work out the weight of
C.. the new graph
         WEIGHT=0.D0
         IF(IACC)
     &       WEIGHT=CALCPATHS_N(NIPATH,0.D0,NRHOIJ,I_V,I_HMAX,
     &               I_P,0.D0,NEL,I_V,ILOGGING,DBETA,ETILDE,NHIJS,ICLS)
C.. we now accept or reject with the correct prob
C.. Pacc(G'|G)=min(1,|w'(G')| Psel(k'|G') w(k)/ (|w(G)| Psel(k|G) w(k')))
         IF(ABS(WEIGHT).LT.WMIN) WEIGHT=0.D0
         WT1=(NRHOIJ(IVERTMOVE,IVERTMOVE)/RHOIJ(0,0))**(I_P*PFAC)
         WT2=1.D0/WT1
         PR=ABS(WEIGHT)*(WT2/(WTTOT-VWEIGHTS(2,IVERTMOVE)+WT2))*WT1
     &      /(ABS(OWEIGHT)*(VWEIGHTS(2,IVERTMOVE)/WTTOT)
     &         *VWEIGHTS(1,IVERTMOVE))
         R=RAN2(ISEED)
         IF(IACC.EQ.1.AND.R.GT.PR) IACC=0
        ENDIF
         IF(BTEST(ILOGGING,10)) THEN
            CALL WRITEPATH(24,IPATH,I_V,NEL,.FALSE.)
            WRITE(24,"(A,$)") "->"
            CALL WRITEPATH(24,NIPATH,I_V,NEL,.FALSE.)
            WRITE(24,*) IACC
            WRITE(24,"(7E25.16)") PR,WEIGHT,
     &         (WT2/(WTTOT-VWEIGHTS(2,IVERTMOVE)+WT2)),WT1,
     &         OWEIGHT,VWEIGHTS(2,IVERTMOVE)/WTTOT,
     &         VWEIGHTS(1,IVERTMOVE)
C  ,VWEIGHTS(2,IVERTMOVE)/WTTOT
         ENDIF 
         IF(IACC) THEN
C.. Accept the new
            OETILDE=ETILDE
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            CALL DCOPY((I_V+1)**2,NRHOIJ,1,RHOIJ,1)
            CALL DCOPY(I_V+1,NHIJS,1,HIJS,1)
            CALL ICOPY(NEL*(1+I_V),NIPATH,1,IPATH,1)
            VWEIGHTS(1,I_V)=VWEIGHTS(1,I_V)-VWEIGHTS(1,IVERTMOVE)+WT1
            VWEIGHTS(2,I_V)=VWEIGHTS(2,I_V)-VWEIGHTS(2,IVERTMOVE)+WT2
            VWEIGHTS(1,IVERTMOVE)=WT1
            VWEIGHTS(2,IVERTMOVE)=WT2
C.. Free the old excitor
            CALL FREEM(VEXCITS(IVERTMOVE))
C.. Now create a new excitor for this new det.
            CALL GENRANDSYMEXCITIT(IPATH(1,IVERTMOVE),NEL,G1,NBASIS,
     &               NBASISMAX,.TRUE.,CUREXLEN,0,ISEED,IC,0)
C            WRITE(24,"(2I,$)") IP_CUREX,CUREXLEN
            CALL MEMORY(IP_CUREX,CUREXLEN,'CUREX')
C            WRITE(24,"(I)") IP_CUREX
            CUREX(1)=0
            VEXCITS(IVERTMOVE)=IP_CUREX
         ELSE
            ETILDE=OETILDE
            WEIGHT=OWEIGHT
            ICLS=IOCLS
         ENDIF         
         IF(WEIGHT.NE.0.D0) THEN
C.. Unbias the sum 
            ETILDERET=ETILDE/(ABS(WEIGHT))
            WEIGHTRET=WEIGHT/(ABS(WEIGHT))
         ELSE
            ETILDERET=0.D0
            WEIGHTRET=0.D0
         ENDIF
C.. If we've a new graph, IOCLS=0, and classpaths regenerates the class.
C.. otherwise, we just reuse the old version
C         IF(TST.EQ.1) THEN
C            TST=1
C         ENDIF
         CALL CLASSPATHS(WEIGHT,ETILDE,1.D0,RHOIJ,
     &            I_V,IOCLS)
         CALL GETTREENESS(IOCLS,ITREE,WEIGHTRET,I_V)
         IF(BTEST(ILOGGING,2)) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            WRITE(10,"(3E25.16, 3I5,A2)") WEIGHT,PR,ETILDE,IOCLS,IACC,
     &         ITREE,"M"
            WRITE(10,"(3E25.16)") (HIJS(I),I=0,2)
         ENDIF
         RETURN
      END

C.. Add a vertex to a graph.  We do this with prob
C.. Pacc(G   ->G')
C.. Pacc(ijk->ijkl) =
C..   |w'(ijkl)| w2(l)/Sum_(j' in G') w2(j')
C.. -------------------------------------------------------------------
C..  (|w'(ijk)| Sum_(j' in G&j'-l) w1(j')/N(j'))/Sum_(j''in G) w1(j'')
      SUBROUTINE MCPR4D4INC(IPATH,NEL,I_V,RHOIJ,I_HMAX,I_P,BETA,
     &                        HIJS,VWEIGHTS,VEXCITS,G1,NBASIS,
     &                        NBASISMAX,OWEIGHT,ILOGGING,DBETA,RHOEPS,
     &                        OETILDE,ETILDERET,WEIGHTRET,IOCLS,
     &                        ITREE,IACC,ISEED,BRR,NMSH,FCK,TMAT,NTAY,
     &                        NMAX,ALAT,UMAT,ECORE,PFAC,WMIN)
         IMPLICIT NONE
         INTEGER I_V,NEL,IPATH(NEL,0:I_V)
         REAL*8 RHOIJ(0:I_V,0:I_V)
         REAL*8 HIJS(0:I_V),VWEIGHTS(2,0:I_V)
         INTEGER VEXCITS(0:I_V+1)
         REAL*8 RHOIJ2(0:I_V+1,0:I_V+1)
         REAL*8 HIJS2(0:I_V+1)
         INTEGER IPATH2(NEL,0:I_V+1)
         INTEGER ISEED,IC
         REAL*8 R,WTL,WTJ,PR,WT1,WT2
         INTEGER IVERT,IVERTMOVE,I,J,I2,J2
         REAL*8 WEIGHT,OWEIGHT,ETILDE,OETILDE
         POINTER (IP_CUREX,CUREX)
         INTEGER CUREX(*),ICE,NL(NEL),CUREXLEN
   
         INTEGER BRR(*),NMSH,I_HMAX,NMAX,I_P,NTAY
         REAL*8 FCK,TMAT(*),ALAT(*),UMAT(*),ECORE
         REAL*8 DBETA,RHOEPS,BETA

         INTEGER G1(*),NBASIS,NBASISMAX(*),ILOGGING

         REAL*8 ETILDERET,WEIGHTRET
         INTEGER ITREE,IACC,ICLS,IOCLS
         REAL*8 RH,PFAC,WMIN
 
         LOGICAL AREDETSEXCITS,ISVALIDDET,LISINPATH
         REAL*8 GETHELEMENT2,CALCPATHS_N
         LOGICAL LBR
         INTEGER ICOUNT
         REAL*8 RAN2
C.. First pick a vertex to excite from 
         R=RAN2(ISEED)*VWEIGHTS(1,I_V)
         IVERT=0
         DO WHILE(VWEIGHTS(1,IVERT).LT.R)
            R=R-VWEIGHTS(1,IVERT)
            IVERT=IVERT+1
         ENDDO
         IF(IVERT.GE.I_V) IVERT=I_V-1
         IVERTMOVE=IVERT
C.. Now use the excitor for this vertex to generate a new vertex
         IP_CUREX=VEXCITS(IVERTMOVE)
         ICOUNT=0
         LBR=.TRUE.
         IACC=1
         DO WHILE(LBR)
            ICOUNT=ICOUNT+1
            CALL GENRANDSYMEXCITIT(IPATH(1,IVERTMOVE),NEL,G1,NBASIS,
     &         NBASISMAX,.FALSE.,CUREX,NL,ISEED,ICE,0)
            IF(.NOT.ISVALIDDET(NL,NEL)) THEN
               WRITE(6,*) "INVALID DET"
               CALL WRITEDET(6,NL,NEL,.TRUE.)
               STOP "INVALID DET"
            ENDIF
            IF(.NOT.LISINPATH(NL,IPATH,NEL,I_V,-1)) THEN
               LBR=.FALSE.
            ENDIF
            
            IF(ICOUNT.GT.100) THEN
C               WRITE(6,*) "Vertices:",I_V
C               CALL WRITEDET(6,IPATH(1,0),NEL,.TRUE.)
C               STOP "Unable to extend graph after 100 tries"
               IACC=0
               LBR=.FALSE.
            ENDIF
         ENDDO
         WTJ=0.D0
C.. We've now generated our new node.  We need to find out the
C.. probability of having generated it.  Go through all nodes, and
C.. work out how connected it is.
         IF(IACC.EQ.1) THEN
         DO I=0,I_V-1
C.. fill up the new rho matrix
            CALL DCOPY(I_V,RHOIJ(0,I),1,RHOIJ2(0,I),1)
            RHOIJ2(I_V+1,I)=RHOIJ(0,I)
            RHOIJ2(I,I_V+1)=RHOIJ(0,I)
            IF(AREDETSEXCITS(NL,
     &                       IPATH(1,I),NEL,NBASISMAX,ICE)) THEN
C.. we're joined.  Calculate the total number of dets connected to this det
               IP_CUREX=VEXCITS(I)
               CALL GETSYMEXCITCOUNT(IPATH(1,I),NEL,G1,NBASIS,NBASISMAX,
     &            CUREX,IC,0)
C.. Remove the number that are in this graph, G
               DO J=0,I_V-1
                  IF(I.NE.J.AND.RHOIJ(I,J).NE.0) IC=IC-1
               ENDDO
C.. discount the vertex we're about to add IVERTMOVE, as we are working
C.. out the probability of having generated the new node from the old
C.. Graph, G
C.. now accumulate the weight of this graph
               WTJ=WTJ+VWEIGHTS(1,IVERTMOVE)/IC
               
C.. now get the rho matrix element
               CALL CALCRHO2(NL,IPATH(1,I),BETA,I_P,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,
     &                  UMAT,RH,NTAY,ICE,ECORE)
               IF(ABS(RH).GE.RHOEPS) THEN
                  RHOIJ2(I,I_V)=RH
                  RHOIJ2(I_V,I)=RH
                  IF(I.EQ.0)
     &               HIJS2(I_V)=GETHELEMENT2(IPATH(1,0),NL,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,
     &                  NMAX,ALAT,UMAT,ICE,ECORE)
               ELSE
                  RHOIJ2(I,I_V)=0.D0
                  RHOIJ2(I_V,I)=0.D0
                  IF(I.EQ.0) HIJS2(I_V)=0.D0
               ENDIF
            ELSE
               RHOIJ2(I,I_V)=0.D0
               RHOIJ2(I_V,I)=0.D0
               IF(I.EQ.0) HIJS2(I_V)=0.D0
            ENDIF 
         ENDDO
C.. Normalize WTJ
         WTJ=WTJ/VWEIGHTS(1,I_V)
C.. Finish off the rho matrix, and copy the H mat etc.
         CALL CALCRHO2(NL,NL,BETA,I_P,NEL,
     &                  NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMAT,NMAX,ALAT,
     &                  UMAT,RH,NTAY,0,ECORE)
         RHOIJ2(I_V,I_V)=RH
         RHOIJ2(I_V+1,I_V+1)=RHOIJ(0,0)
         RHOIJ2(I_V,I_V+1)=RHOIJ2(0,I_V)
         RHOIJ2(I_V+1,I_V)=RHOIJ2(0,I_V)
         CALL DCOPY(I_V,HIJS,1,HIJS2,1)
         HIJS2(I_V+1)=HIJS(0)
         HIJS2(0)=HIJS(0)
         CALL ICOPY(NEL*I_V,IPATH(1,0),1,IPATH2(1,0),1)
         CALL ICOPY(NEL,IPATH(1,0),1,IPATH2(1,I_V+1),1)
         CALL ICOPY(NEL,NL,1,IPATH2(1,I_V),1)

C.. Now work out the weight of the new det
         WT1=(RH/RHOIJ(0,0))**(I_P*PFAC)
         WT2=1.D0/WT1
C.. and the weight of the new graph
         WEIGHT=CALCPATHS_N(IPATH2,0.D0,RHOIJ2,I_V+1,I_HMAX,
     &              I_P,0.D0,NEL,I_V+1,ILOGGING,DBETA,ETILDE,HIJS2,ICLS)
         IF(ABS(WEIGHT).LT.WMIN) WEIGHT=0.D0
         PR=ABS((WEIGHT*WT2/(VWEIGHTS(2,I_V)+WT2-VWEIGHTS(2,0)))
     &         /(OWEIGHT*WTJ))
         R=RAN2(ISEED)
        ELSE
C.. fake a ipath2
         CALL ICOPY(NEL*(I_V+1),IPATH,1,IPATH2,1)
         CALL ICOPY(NEL,IPATH,1,IPATH2(1,I_V+1),1)
         PR=0.D0
         WEIGHT=0.D0
         ETILDE=0.D0
        ENDIF
         IF(IACC.EQ.1.AND.R.GE.PR) IACC=0
         IF(BTEST(ILOGGING,10)) THEN
            CALL WRITEPATH(24,IPATH,I_V,NEL,.FALSE.)
            WRITE(24,"(A,$)") "->"
            CALL WRITEPATH(24,IPATH2,I_V+1,NEL,.FALSE.)
            WRITE(24,*) IACC 
            WRITE(24,"(5E25.16)") PR,WEIGHT,
     &         WT2/(VWEIGHTS(2,I_V)+WT2-VWEIGHTS(2,0)),
     &         OWEIGHT,WTJ
         ENDIF 
         IF(IACC) THEN
C.. Accept the new
            OETILDE=ETILDE
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            CALL DCOPY((I_V+2)**2,RHOIJ2,1,RHOIJ,1)
            CALL DCOPY(I_V+2,HIJS2,1,HIJS,1)
            CALL ICOPY(NEL*(I_V+2),IPATH2,1,IPATH,1)
            VWEIGHTS(1,I_V+1)=VWEIGHTS(1,I_V)+WT1
            VWEIGHTS(2,I_V+1)=VWEIGHTS(2,I_V)+WT2
            VWEIGHTS(1,I_V)=WT1
            VWEIGHTS(2,I_V)=WT2
            
C.. Now create a new excitor for this new det.
             
            CALL GENRANDSYMEXCITIT(NL,NEL,G1,NBASIS,NBASISMAX,
     &               .TRUE.,CUREXLEN,0,ISEED,IC,0)
            CALL MEMORY(IP_CUREX,CUREXLEN,'CUREX')
            CUREX(1)=0
            VEXCITS(I_V)=IP_CUREX
            I_V=I_V+1
         ELSE
            ETILDE=OETILDE
            WEIGHT=OWEIGHT
            ICLS=IOCLS
         ENDIF         
         IF(WEIGHT.NE.0.D0) THEN
C.. Unbias the sum 
            ETILDERET=ETILDE/(ABS(WEIGHT))
            WEIGHTRET=WEIGHT/(ABS(WEIGHT))
         ELSE
            ETILDERET=0.D0
            WEIGHTRET=0.D0
         ENDIF
C.. If we've a new graph, IOCLS=0, and classpaths regenerates the class.
C.. otherwise, we just reuse the old version
         CALL CLASSPATHS(WEIGHT,ETILDE,1.D0,RHOIJ,
     &            I_V,IOCLS)
         CALL GETTREENESS(IOCLS,ITREE,WEIGHTRET,I_V)
         IF(BTEST(ILOGGING,2)) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            WRITE(10,"(3E25.16, 3I5,A2)") WEIGHT,PR,ETILDE,IOCLS,IACC,
     &         ITREE,"A"
            WRITE(10,"(3E25.16)") (HIJS(I),I=0,2)
         ENDIF   
         RETURN
      END


C.. Remove a vertex from a graph.  We do this with prob
C.. Pacc(G   ->G')
C.. Pacc(ijkl->ijk) =
C..  (|w'(ijk)| Sum_(j' in G'&j'-l) w1(j')/N(j'))/Sum_(j''in G') w1(j'')
C.. -------------------------------------------------------------------
C..   |w'(ijkl)| w2(l)/Sum_(j' in G) w2(j')

C.. N(j) is the number of dets connected to j which are not in G'

C.  We choose the vertex to remove with prob rho_ii^-P (w2)
      SUBROUTINE MCPR4D4DEC(IPATH,NEL,I_V,RHOIJ,I_HMAX,I_P,
     &                        HIJS,VWEIGHTS,VEXCITS,G1,NBASIS,
     &                        NBASISMAX,OWEIGHT,ILOGGING,DBETA,
     &                        OETILDE,ETILDERET,WEIGHTRET,IOCLS,
     &                        ITREE,IACC,ISEED,PFAC,WMIN)
         IMPLICIT NONE
         INTEGER I_V,NEL,IPATH(NEL,0:I_V)
         REAL*8 RHOIJ(0:I_V,0:I_V)
         REAL*8 HIJS(0:I_V),VWEIGHTS(2,0:I_V)
         INTEGER VEXCITS(0:I_V)
         REAL*8  RHOIJ2(0:I_V-1,0:I_V-1)
         REAL*8  HIJS2(0:I_V-1)
         INTEGER IPATH2(NEL,0:I_V-1)
         INTEGER ISEED
         REAL*8 R,WTL,WTJ,PR
         INTEGER IVERT,IVERTMOVE,I,J,I2,J2
         REAL*8 WEIGHT,OWEIGHT,ETILDE,OETILDE
         POINTER (IP_CUREX,CUREX)
         INTEGER CUREX(*),ICE,IC


         INTEGER G1(*),NBASIS,NBASISMAX(*),ILOGGING
         INTEGER I_P,I_HMAX

         REAL*8 ETILDERET,WEIGHTRET,DBETA,BETA
         INTEGER ITREE,IACC,ICLS,IOCLS
         REAL*8 RH,PFAC,WMIN
 
         LOGICAL AREDETSEXCITS,ISVALIDDET
         REAL*8 CALCPATHS_N
         REAL*8 RAN2
C.. Choose the vertex to remove
C.. (we cannot remove the pivot, so we exclude that)
         R=RAN2(ISEED)*(VWEIGHTS(2,I_V)-VWEIGHTS(2,0))
         IVERT=1
         DO WHILE(VWEIGHTS(2,IVERT).LT.R)
            R=R-VWEIGHTS(2,IVERT)
            IVERT=IVERT+1
         ENDDO
         IF(IVERT.GE.I_V) IVERT=I_V-1
         IVERTMOVE=IVERT
C.. Normalize WTL
         WTL=VWEIGHTS(2,IVERTMOVE)/(VWEIGHTS(2,I_V)-VWEIGHTS(2,0))
C.. Now we need to traverse the RHO matrix to find to which dets
C.. we are joined
         WTJ=0.D0
         DO I=0,I_V-1
C..
C.. An interesting question at this point is should we use dets which 
C.. are actually connected, or merely could be connected.  The latter
C.. seems more sensible.
            IF(I.NE.IVERTMOVE.AND.
     &         AREDETSEXCITS(IPATH(1,IVERTMOVE),
     &                       IPATH(1,I),NEL,NBASISMAX,ICE)) THEN
C.. we're joined.  Calculate the total number of dets connected to this det
               IP_CUREX=VEXCITS(I)
               CALL GETSYMEXCITCOUNT(IPATH(1,I),NEL,G1,NBASIS,NBASISMAX,
     &            CUREX,IC,0)
C.. Remove the number that are in this graph
               DO J=0,I_V-1
                  IF(I.NE.J.AND.RHOIJ(I,J).NE.0) IC=IC-1
               ENDDO
C.. discount the vertex we're about to remove IVERTMOVE, as we need 
C.. to remove dets connected in G' not G
               IC=IC+1
C.. now accumulate the weight of this graph
               WTJ=WTJ+VWEIGHTS(1,I)/IC
            ENDIF 
         ENDDO
C.. Normalize WTJ
         WTJ=WTJ/(VWEIGHTS(1,I_V)-VWEIGHTS(1,IVERTMOVE))
C.. Now we have all the relevant generating probabilities, we can create
C.. a new matrix without the removed det, and generate its weight.
         I2=0
         DO I=0,I_V
            IF(I.NE.IVERTMOVE) THEN
               HIJS2(I2)=HIJS(I)
               CALL ICOPY(NEL,IPATH(1,I),1,IPATH2(1,I2),1)
               J2=0
               DO J=0,I_V
                  IF(J.NE.IVERTMOVE) THEN         
                     RHOIJ2(I2,J2)=RHOIJ(I,J)
                     J2=J2+1
                  ENDIF
               ENDDO
               I2=I2+1
            ENDIF
         ENDDO      
         WEIGHT=CALCPATHS_N(IPATH2,0.D0,RHOIJ2,I_V-1,I_HMAX,
     &              I_P,0.D0,NEL,I_V-1,ILOGGING,DBETA,ETILDE,HIJS2,ICLS)
         IF(ABS(WEIGHT).LT.WMIN) WEIGHT=0.D0
         PR=ABS(WEIGHT*WTJ/(OWEIGHT*WTL))
         R=RAN2(ISEED)
         IACC=0
         IF(R.LT.PR) IACC=1
         IF(BTEST(ILOGGING,10)) THEN
            CALL WRITEPATH(24,IPATH,I_V,NEL,.FALSE.)
            WRITE(24,"(A,$)") "->"
            CALL WRITEPATH(24,IPATH2,I_V-1,NEL,.FALSE.)
            WRITE(24,*) IACC 
            WRITE(24,"5E25.16") PR,WEIGHT,WTJ,OWEIGHT,WTL
         ENDIF 
         IF(IACC) THEN
C.. Accept the new
            OETILDE=ETILDE
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
            CALL DCOPY(I_V**2,RHOIJ2,1,RHOIJ,1)
            CALL DCOPY(I_V,HIJS2,1,HIJS,1)
            CALL ICOPY(NEL*I_V,IPATH2,1,IPATH,1)
            VWEIGHTS(1,I_V)=VWEIGHTS(1,I_V)-VWEIGHTS(1,IVERTMOVE)
            VWEIGHTS(2,I_V)=VWEIGHTS(2,I_V)-VWEIGHTS(2,IVERTMOVE)
C.. Free the excitor
            CALL FREEM(VEXCITS(IVERTMOVE))
            DO I=IVERTMOVE,I_V-1
               VWEIGHTS(1,I)=VWEIGHTS(1,I+1)
               VWEIGHTS(2,I)=VWEIGHTS(2,I+1)
               VEXCITS(I)=VEXCITS(I+1)
            ENDDO
            I_V=I_V-1
         ELSE
            ETILDE=OETILDE
            WEIGHT=OWEIGHT
            ICLS=IOCLS
         ENDIF         
         IF(WEIGHT.NE.0.D0) THEN
C.. Unbias the sum 
            ETILDERET=ETILDE/(ABS(WEIGHT))
            WEIGHTRET=WEIGHT/(ABS(WEIGHT))
         ELSE
            ETILDERET=0.D0
            WEIGHTRET=0.D0
         ENDIF
C.. If we've a new graph, IOCLS=0, and classpaths regenerates the class.
C.. otherwise, we just reuse the old version
         CALL CLASSPATHS(WEIGHT,ETILDE,1.D0,RHOIJ,
     &            I_V,IOCLS)
         CALL GETTREENESS(IOCLS,ITREE,WEIGHTRET,I_V)
         IF(BTEST(ILOGGING,2)) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(BTEST(ILOGGING,3))
     &         CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            WRITE(10,"(3E25.16, 3I5,A2)") WEIGHT,PR,ETILDE,IOCLS,IACC,
     &         ITREE,"D"
            WRITE(10,"(3E25.16)") (HIJS(I),I=0,2)
         ENDIF   
         RETURN
      END



