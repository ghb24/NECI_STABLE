module mcpaths
    use util_mod, only: isnan
contains

!C.. Calculate RHO^(P)_II without having a stored H matrix
!C.. Sum over distinct nodes, e.g. IJKLI, with paths up to I_HMAX
!C.. generated from these, and summed (e.g IJILKJI), up to max H
!C.. In theory more efficient because RHO_IJ,RHO_JK, etc are calculated
!C.. just once for all these paths.
!C.. I_VMAX is the max number of distinct vertices in a path.
!C.. I_HMAX is the max number of hops in a path.
!C.. If we must calculate DLWDB along with WI, we set its value.
!C.. If we leave it as it is, we will be called again to calculated it.
!C..  MCPATHSR10 allows different methods at each vertex level
!C.. DLWDB = d ln (wi) / d Beta.

!..  wi is split into ln wi = P ln rii + ln si
!..  rii = rho_ii = <Di|exp(-beta H/P)|Di>


      SUBROUTINE  MCPATHSR10(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1, &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,NTAY,RHOEPS,      &
     &               LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,       &
     &               WLRI,WLSI,DBETA,DLWDB)
         Use Determinants, only: get_helement, write_det
         use constants, only: dp
         use SystemData, only: BasisFN
!         USE PreCalc ,  only: GETVARS
         USE STARDIAGMOD , only: fMCPR3StarNewExcit
         USE GraphMorph , only : MorphGraph
         USE StarDiagTripMod , only : StarDiagTrips
#ifdef PARALLEL
         USE FciMCParMod , only : FciMCPar
#endif
         USE FciMCMod , only : FciMC
         USE ReturnPathMCMod , only : ReturnPathMC
         USE NODEDIAG , only : fMCPR3StarNodes
         use CalcData , only : G_VMC_FAC,CUR_VERT,g_MultiWeight,         &
     &          TMPTHEORY,TMCDIRECTSUM,TDIAGNODES,TGraphMorph,           &
     &          calcp_logweight,TFCIMC,TReturnPathMC,tFCIMCSerial 
         use CalcData, only: tCCMC
         use CalcData, only: TStarTrips
         USE Logging , only : G_VMC_LOGCOUNT
         USE PrecalcData , only : PREIV_MAX,TPREVAR
         USE CCMC, only: CCMCStandalone,CCMCStandaloneParticle
         use CCMCData, only:  tAmplitudes
         use global_utilities
         use mcpathsdata, only: EGP
         use mcpathshdiag, only: fmcpr3b2
         use util_mod, only: isnan, NECI_ICOPY
         IMPLICIT NONE
         TYPE(BasisFN) :: G1(*),KSYM
         INTEGER I_VMAX,NEL,NBASIS
         INTEGER IPATH(NEL,0:I_VMAX)
         INTEGER NI(NEL)
         COMPLEX*16 FCK(*)
         INTEGER I_P,I_HMAX,BRR(*),NMSH,NMAX
         INTEGER NTAY(2),NWHTAY(3,I_VMAX),ILOGGING,I,I_V
         INTEGER L,LT,K,J
         real*4 otime,itime,etime,tarr(2)
         REAL*8 BETA,ECORE
         real(dp) WLRI,WLSI
         HElement_t UMat(*),RH
         REAL*8 NTOTAL

         real(dp) F(2:I_VMAX)
         CHARACTER*40 STR
         real(dp) TOTAL,RHOII(0:I_VMAX)
         HElement_t RHOIJ(0:I_VMAX,0:I_VMAX)
         real(dp) FLNSUM,FLNRII
         REAL*8 ALAT(3),RHOEPS
         INTEGER  nBasisMax(5,*)
         INTEGER LSTE(*),ILMAX
!CNEL,0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER ICE(*)
         HElement_t RIJLIST(:,:)
!C0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER NLIST(0:I_VMAX-1),LSTP(0:I_VMAX-1),BTABLE(0:I_VMAX)
         LOGICAL TLOG,TSYM
         REAL*8 DBETA
         real(dp) DLWDB,DLWDB2,EREF,ODLWDB2
         HElement_t  HIJS(0:I_VMAX)
         TYPE(EGP) LOCTAB(I_VMAX)
         real(dp) FMCPR4B,FMCPR4C,FMCPR4D,FMCPR4D2
         real(dp) FMCPR3STAR,FMCPR3NVSTAR
         INTEGER I_CHMAX,CNWHTAY
         INTEGER ISEED,ICOUNT
         real(dp) FF,DLWDB3,TEMPTOT
         REAL*8 FSQ(2:I_VMAX),STD
         INTEGER NMEMLEN,INODE2(NEL)
         REAL*8 NTIME,VARSUM
         INTEGER IFRZ(0:NBASIS,I_VMAX),I_V1,I_V2,T
         REAL*8 OSI,ODLWDB,OPROB
         real(dp) MP2E(2:I_VMAX)
         LOGICAL TLOGP
         type(timer), save :: proc_timer
         type(timer), save :: proc_timer2
         type(timer), save :: proc_timerpre
         proc_timer%timer_name='MCPATHSR10'
         call set_timer(proc_timer)
         WRITE(6,*) "Running Multi-level graph calculation."
!C.. This will store the MP2 energy
         DO I=2,I_VMAX
            MP2E(I)=0.D0
         ENDDO
         ISEED=7
         TLOGP=BTEST(ILOGGING,14)
         TLOG=BTEST(ILOGGING,1)
         IF(NTAY(1).LE.0) THEN
            TSYM=.FALSE.
         ELSE
            TSYM=.TRUE.
         ENDIF
         IF(TLOG) THEN
            OPEN(11,FILE="MCSUMMARY",STATUS="OLD",POSITION='APPEND')
!C.. go to end of file
!            I=FSEEK(11,0,2)
            call write_det (11, NI, .true.)
         ENDIF
         IF(BTEST(ILOGGING,0)) THEN
            OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
            IF(TMPTHEORY) OPEN(13,FILE="MP2PATHS",STATUS="UNKNOWN")
         ENDIF
         IF(BTEST(ILOGGING,9)) OPEN(12,FILE="VERTEXMC",STATUS="UNKNOWN")
!C.. Set the first node to I_I
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         HIJS(0) = get_helement (nI, nI, 0)
         RHOII(0)=RH
         RHOIJ(0,0)=RH
         WLRI=LOG(RHOII(0))
         TOTAL=1.D0
         DLWDB=HIJS(0)
         IF(TLOG) WRITE(11,"(I12,2G25.16,F19.7,2I12,2F19.7)") 1,TOTAL,TOTAL,0.D0,1,1,HIJS(0)

         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)


         NTOTAL=1.D0
         I_V1=0

         !PRECALC block called
         if(PREIV_MAX.ne.0) then
             
           TPREVAR=.TRUE. 
            
            IF(.not.TMCDIRECTSUM) THEN
                WRITE(6,*) "PRECALC only valid if MCDIRECTSUM is being used"
                CALL FLUSH(6)
                STOP
            ENDIF
            
            proc_timerPRE%timer_name='PRECALC '
            call set_timer(proc_timerPRE)
            CALL GETVARS(NI,BETA,I_P,IPATH,2,                           &
     &       G1,NMSH,FCK,NMAX,UMAT,                                     &
     &       NTAY,RHOEPS,RHOII,RHOIJ,LOCTAB,TSYM,ECORE,DBETA,DLWDB2,    &
     &       HIJS,L,LT,IFRZ,MP2E,NTOTAL,DLWDB,TOTAL,TLOGP,KSYM,NWHTAY,  &
     &       I_VMAX)
            call halt_timer(proc_timerPRE)
         end if

!C.. I_V is the number of vertices in the path
         
         DO I=2,I_VMAX
            I_V=NWHTAY(3,I)
!CI_P
!c            WRITE(6,"I2") I_V
            CUR_VERT=I_V
            WRITE(STR,"(A,I5)") "FMCPR",I_V
            proc_timer2%timer_name=STR
            call set_timer(proc_timer2)
            OTIME=etime(tarr)
            L=0
            LT=0
            BTABLE(0)=-1
            IF(TPREVAR) THEN
                DO T=2,I_VMAX
                    MP2E(T)=0.D0
                ENDDO
                NTOTAL=1
                TOTAL=1.D0
                DLWDB=HIJS(0)
                TPREVAR=.FALSE.
            ENDIF

            DLWDB2=0.D0
            I_CHMAX=NWHTAY(1,I)
            CNWHTAY=NWHTAY(2,I)
            IF(I_CHMAX.EQ.-1) THEN
!C.. This code generates list of all excitations at each level
            WRITE(6,"(A6,I3,A,I3,A)") "Level ",I," Full Sum Old ",I_V," Vertex"
            F(I_V)=FMCPR3(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,G1,NBASIS,   &
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,                 &
     &            0,RHOII,RHOIJ,LSTE,ICE,RIJLIST,                          &
     &            CNWHTAY,L,LT,I_CHMAX,NLIST,LSTP,                         &
     &            BTABLE,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB2,HIJS,      &
     &            MP2E,NTOTAL)
            ELSEIF(I_CHMAX.EQ.-8) THEN
            WRITE(6,"(A6,I3,A,I3,A)") "Level ",I," Full Sum New ",I_V," Vertex"

               EREF=DLWDB/TOTAL
            
!C.. This code generates excitations on the fly
               F(I_V)=FMCPR3B(NI,BETA,I_P,IPATH,I_V,NEL,                   &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,        &
     &         RHOEPS,0,RHOII,RHOIJ,CNWHTAY,I_CHMAX,LOCTAB,                &
     &         ILOGGING,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,       &
     &         MP2E,NTOTAL,EREF,VARSUM,TOTAL)
            ELSEIF(I_CHMAX.EQ.-14) THEN
!.. We read in this vertex level from the MCPATHS file
            WRITE(6,"(A6,I3,A,I3,A)") "Level ",I," MCPATHS read ",I_V," Vertex"
               CALL ReadMCPaths(I_V,F(I_V),DLWDB2)
            ELSEIF(I_CHMAX.EQ.-20) THEN
            WRITE(6,"(A6,I3,A,I3,A)") "Level ",I," Full Sum HDiag ",I_V," Vertex"
!C.. This code generates excitations on the fly, and diagonalizes a
!C.. matrix of hij instead of rhoij

        EREF=DLWDB/TOTAL
!        write (6,*) "from mcpaths.F", EREF 
            
               F(I_V)=FMCPR3B2(NI,BETA,I_P,IPATH,I_V,NEL,                  &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,        &
     &         RHOEPS,0,RHOIJ,CNWHTAY,I_CHMAX,LOCTAB,                      &
     &         ILOGGING,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,       &
     &        MP2E,NTOTAL,I_VMAX,EREF,VARSUM,TOTAL)
            ELSEIF(I_CHMAX.EQ.-9) THEN
!C.. This code generates a star consisting of all the two-vertex terms,
!C.. forcing them to be disconnected.
               F(I_V)=FMCPR3STAR(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,      &
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,                 &
     &            LSTE,ICE,RIJLIST,L,LT,                                   &
     &            CNWHTAY,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB2)
            ELSEIF(I_CHMAX.EQ.-21) THEN
               IF(TDIAGNODES) THEN
! This code prediagonalises connected nodes of virtual orbitals, before solving as a star.
                F(I_V)=fMCPR3StarNodes(NI,BETA,I_P,NEL,NBASISMAX,G1,       &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,        &
     &              L,LT,CNWHTAY,ILOGGING,TSYM,ECORE,DBETA,DLWDB2)         
               ELSEIF(TGraphMorph) THEN
!This involves iterativly improving a graph in order to maximise the connections between determinants.
                    IF(I_VMAX.ne.2) STOP 'Error in vertex number'
                    CALL MorphGraph(F(I_V),DLWDB2)
               ELSEIF(TStarTrips) THEN
!Triple excitations are prediagonalised from each double
                    CALL StarDiagTrips(DLWDB2,F(I_V))
               ELSEIF(TFCIMC) THEN
!A MC simulation involving replicating particles is run
#ifdef PARALLEL
!                    WRITE(6,*) "Get Here!: ",I_V,F(I_V),DLWDB2
                  IF(tFCIMCSerial) THEN
                    CALL FciMC(F(I_V),DLWDB2)
                  ELSE
                    CALL FciMCPar(F(I_V),DLWDB2)
                  ENDIF
!                    WRITE(6,*) "Get Here!: ",I_V,F(I_V),DLWDB2
#else
                    CALL FciMC(F(I_V),DLWDB2)
#endif
               ELSEIF(tCCMC) THEN
                  if(tAmplitudes) THEN
                     CALL CCMCStandAlone(F(I_V),DLWDB2)
                  else
                     CALL CCMCStandaloneParticle(F(I_V),DLWDB2)
                  endif
               ELSEIF(TReturnPathMC) THEN
!A MC simulation involving replicating particles, constrained to returning paths is run
                    CALL ReturnPathMC(F(I_V),DLWDB2)
               ELSE
!C.. This code generates a star consisting of all the two-vertex terms,
!C.. forcing them to be disconnected. - this uses new excitation generators
               F(I_V)=FMCPR3STARNewExcit(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS, &
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,L,LT,               &
     &            CNWHTAY,ILOGGING,TSYM,ECORE,DBETA,DLWDB2,MP2E)
               ENDIF
            ELSEIF(I_CHMAX.EQ.-11) THEN
               IF(I_V.EQ.I_VMAX) THEN
               F(I_V)=FMCPR3NVSTAR(NI,BETA,I_P,IPATH,I_V,NEL,                 &
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,                &
     &   RHOEPS,0,RHOII,RHOIJ,LSTE,ICE,RIJLIST,CNWHTAY,L,LT,I_HMAX,           &
     &   NLIST,LSTP,BTABLE,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB,HIJS)
               ELSE
                  F(I_V)=0.D0
                  L=0
                  LT=0
               ENDIF
            ELSEIF(I_CHMAX.EQ.-7.OR.I_CHMAX.EQ.-19) THEN
             IF(G_VMC_FAC.NE.0.D0) THEN
!C..   We couple a vertex monte carlo at this vertex level (and subsequent is g_MultiWeight has non-zero values) to a known result at previous
!C.. levels by regarding all previous levels as a composite object with the Et and wi
!C.. we have already worked out.g_MultiWeight(0) contains the sum of the weightings of the levels 
                g_MultiWeight(0)=0.D0
                DO J=1,I_VMAX
                  IF(g_MultiWeight(J).NE.0.D0) THEN
                      IF(g_MultiWeight(0).EQ.0.D0) THEN
!  The first of the multi - get its vertex level
                          I_V1=NWHTAY(3,J)
                      ENDIF
                      g_MultiWeight(0)=g_MultiWeight(0)+g_MultiWeight(J)
                  ENDIF
                ENDDO
              IF(g_MultiWeight(0).EQ.0.D0) THEN
! In fact we're not doing a multi
                 I_V1=I_V
                 I_V2=I_V
              ELSE
!  The last of the multi
                 I_V2=NWHTAY(3,I_VMAX)
              ENDIF
             ELSE  !If G_VMC_FAC.EQ.0
                I_V1=I_V
                I_V2=I_VMAX
             ENDIF
             IF(I_V1.EQ.I_V2) THEN
                IF(tMCDirectSum) THEN
                  WRITE(6,"(A,I3,A,I3,A)") "Level ",I," Direct MC Sum HDiag ",I_V1," Vertex"
               ELSE
                  WRITE(6,"(A,I3,A,I3,A)") "Level ",I," Bias MC Sum HDiag ",I_V1," Vertex"
               ENDIF
             ELSE
               IF(tMCDirectSum) THEN
                  WRITE(6,"(A,I3,A,I3,A,I3,A)") "Level ",I," Direct MC Sum HDiag ",I_V1," to ",I_V2," Vertex"
               ELSE
                  WRITE(6,"(A,I3,A,I3,A,I3,A)") "Level ",I," Bias MC Sum HDiag ",I_V1," to ",I_V2," Vertex"
               ENDIF
            ENDIF
   
      
              WLSI=TOTAL
              CALL MCPATHSR4(NI,BETA,I_P,I_CHMAX,I_V2,NEL,NBASISMAX,G1,       &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,           &
     &               LSTE,ICE,RIJLIST,CNWHTAY,ILOGGING,ECORE,ILMAX,           &
     &               WLRI,WLSI,DBETA,DLWDB,I_V1)
               F(I_V)=WLSI
               IF(TMPTHEORY) THEN
!C  Add the cumulative change
                  MP2E(I_V)=MP2E(I_V)+DLWDB
               ENDIF
!If the routine does many levels of MC, it rebiases them itself, taking into account 
!the DLWDB2 we give it, and adds them into DLWDB.
!If it just does a single level of MC, it does the same.
               DLWDB2=0.D0
              CALL WRITECLASSPATHS()
            ELSE
               F(I_V)=0.D0
               FSQ(I_V)=0.D0
               DLWDB3=0.D0
               OSI=0.D0
               DO ICOUNT=1,CNWHTAY
                  DLWDB2=0.D0
                  IF(I_CHMAX.EQ.-3) THEN
                     FF=FMCPR4B(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,          &
     &                  G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,           &
     &                 RHOEPS,RHOII,RHOIJ,CNWHTAY,I_CHMAX,ILOGGING,TSYM,      &
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ELSEIF(I_CHMAX.EQ.-4) THEN
                     FF=FMCPR4C(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,          &
     &                  G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,           &
     &                 RHOEPS,RHOII,RHOIJ,CNWHTAY,I_CHMAX,ILOGGING,TSYM,      &
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ENDIF
                  F(I_V)=F(I_V)+FF
                  FSQ(I_V)=FSQ(I_V)+abs(FF)**2
                  DLWDB3=DLWDB3+DLWDB2
                  IF(BTEST(ILOGGING,9).AND.MOD(ICOUNT,G_VMC_LOGCOUNT).EQ.0) THEN
                    WRITE(12,"(I3,I10,4F19.7)") I_V,ICOUNT/1000,F(I_V)/ICOUNT,&
     &                   SQRT(FSQ(I_V)/ICOUNT-abs(F(I_V)**2/ICOUNT)),         &   
     &                   DLWDB3/ICOUNT,(DLWDB+DLWDB3/ICOUNT)/(F(I_V)/ICOUNT+TOTAL)
                     CALL FLUSH(12)
                  ENDIF
               ENDDO
               F(I_V)=F(I_V)/CNWHTAY
               FSQ(I_V)=FSQ(I_V)/CNWHTAY
               DLWDB2=DLWDB3/CNWHTAY
               STD=SQRT(ABS(FSQ(I_V)-abs(F(I_V))**2))
            ENDIF

            call halt_timer(proc_timer2)
            NTIME=etime(tarr)
            TOTAL=TOTAL+F(I_V)
            DLWDB=DLWDB+DLWDB2
!            WRITE(6,*) "Get Here2: ",TOTAL,DLWDB
!            WRITE(6,*) "MCP",I_V,TOTAL,DLWDB
            IF(TLOG.AND.I_V1.EQ.0) THEN
                WRITE(11,"(I12,2G25.16,F19.7,I12,2G25.12)") I_V,F(I_V),TOTAL,NTIME-OTIME,L,STD,DLWDB2
               CALL FLUSH(11)
            ENDIF
            IF(ISNAN(TOTAL)) THEN
!C.. save all log files
               ITIME=etime(tarr)
               CALL FLUSH(11)
!               CALL LOGNAN(NI,NEL,BETA,ITIME)
               WRITE(6,*) "WARNING: nan found at time",ITIME
               WRITE(6,"(A)",advance='no') "  nan det="
               call write_det (6, NI, .true.)
            ENDIF
            CALL WRITECLASSPATHS()
! If we've done a multi, we don't need to do any more
            IF(I_V1.NE.0.AND.G_VMC_FAC.NE.0.D0) EXIT
!  Write out partial MP energy info for every cycle
         IF(abs(MP2E(2)).gt.0.D0) THEN
            TEMPTOT=HIJS(0)
            DO J=2,I_VMAX
               TEMPTOT=TEMPTOT+MP2E(J)
               WRITE(STR,"(A10,I1,A10)") "Partial MP",J," ENERGY = "
               WRITE(6,"(A14,G25.16)") STR,TEMPTOT
            ENDDO
         ENDIF
         ENDDO
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,0)) THEN
            CLOSE(10)
            IF(TMPTHEORY) CLOSE(13)
         ENDIF
         IF(abs(MP2E(2)).gt.0.D0) THEN
            FF=HIJS(0)
            DO I=2,I_VMAX
               FF=FF+MP2E(I)
               WRITE(STR,"(A2,I1,A10)") "MP",I," ENERGY = "
               WRITE(6,"(A14,G25.16)") STR,FF
            ENDDO
         ENDIF
         IF(BTEST(ILOGGING,9)) CLOSE(12)
         IF(TOTAL.LT.0) THEN
            DLWDB=-DLWDB
            TOTAL=-TOTAL
         ENDIF
         if(calcp_logweight) then
!  For LogWeight, we are passing around and summing log x[G] in the weights.  The total weight, w_i=prod_G x[G]
!  We return the log of the weight, which is simply the sum of log x[G] which is in total
            wlsi=total
         else
            WLSI=LOG(TOTAL)
!C.. DLWDB contains <D|H exp(-b H)|D>/rhoii^P and we want 
!C.<D|H exp(-b H)|D> / <D|exp(-b H)|D>
!C.. <D|exp(-b H)|D>=exp(P WLRI + WLSI)
!C         WRITE(6,*) DLWDB,TOTAL,EXP(WLRI), HIJS(0)
            DLWDB=DLWDB/TOTAL
         endif
         call halt_timer(proc_timer)
         CALL N_MEMORY_CHECK()
!         WRITE(6,*) "Get Here 3",WLSI,TOTAL,DLWDB
         RETURN
      END SUBROUTINE

!C.. Calculate RHO^(P)_II without having a stored H matrix
!C.. Sum over distinct nodes, e.g. IJKLI, with paths up to I_HMAX
!C.. generated from these, and summed (e.g IJILKJI), up to max H
!C.. In theory more efficient because RHO_IJ,RHO_JK, etc are calculated
!C.. just once for all these paths.
!C.. I_VMAX is the max number of distinct vertices in a path.
!C.. I_HMAX is the max number of hops in a path.
!C.. If we must calculate DLWDB along with WI, we set its value.
!C.. If we leave it as it is, we will be called again to calculated it.
!C.. MCPATHSR3 has the same method at each vertex level
      SUBROUTINE  MCPATHSR3(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,       &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,           &
     &               LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,            &
     &               WLRI,WLSI,DBETA,DLWDB)
         use constants, only: dp
         use SystemData, only: BasisFN
         Use Determinants, only: get_helement, write_det
         USE STARDIAGMOD , only : fMCPR3StarNewExcit
         USE NODEDIAG , only : fMCPR3StarNodes
         USE GraphMorph , only : MorphGraph
         USE StarDiagTripMod , only : StarDiagTrips
#ifdef PARALLEL
         USE FciMCParMod , only : FciMCPar
#else         
         USE FciMCMod , only : FciMC
#endif
         USE ReturnPathMCMod , only : ReturnPathMC
         USE PrecalcData , only : TPREVAR
         use CalcData , only : TMPTHEORY,TDIAGNODES,TGraphMorph
         use CalcData, only : calcp_logweight,TFCIMC,TReturnPathMC
         use CalcData, only: TStarTrips
         use global_utilities
         use mcpathsdata, only: EGP
         use mcpathshdiag, only: fmcpr3b2
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         TYPE(BasisFN) G1(*)
         INTEGER I_VMAX,NEL,NBASIS
         INTEGER IPATH(NEL,0:I_VMAX)
         INTEGER NI(NEL)
         COMPLEX*16 FCK(*)
         HElement_t UMat(*), RH
         INTEGER I_P,I_HMAX,BRR(*),NMSH,NMAX
         INTEGER NTAY(2),NWHTAY(3,I_VMAX),ILOGGING,I,I_V
         type(timer), save :: proc_timer,proc_timer2
         INTEGER L,LT,ITIME,K
         REAL*8 BETA,ECORE
         real*4 etime,tarr(2)
         real(dp) WLRI,WLSI
         real(dp) F(2:I_VMAX)
         CHARACTER*40 STR
         HElement_t RHOIJ(0:I_VMAX,0:I_VMAX)
         real(dp) RHOII(0:I_VMAX),TOTAL,FLNSUM,FLNRII
         REAL*8 ALAT(3),RHOEPS
         INTEGER nBasisMax(5,*)
         INTEGER, allocatable :: LSTE(:,:,:)
         INTEGER ILMAX
!CNEL,0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER, allocatable :: ICE(:,:)
         HElement_t RIJLIST(:,:)
!C0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER NLIST(0:I_VMAX-1),LSTP(0:I_VMAX-1),BTABLE(0:I_VMAX)
         LOGICAL TLOG,TSYM
         REAL*8 DBETA
         HElement_t HIJS(0:I_VMAX)
         type(EGP), target ::  LOCTAB(I_VMAX)
         real(dp) FMCPR4B,FMCPR4C,FMCPR4D,FMCPR4D2
         real(dp) DLWDB,DLWDB2,EREF
         REAL*8 NTIME,OTIME,VARSUM
         INTEGER IFRZ(0:NBASIS,I_VMAX)
         real(dp) MP2E(2:I_VMAX),H00,FF
         REAL*8 NTOTAL
         real(dp) FMCPR3STAR
         write(6,*) "MCPATHSR3:  I_HMAX=",I_HMAX
! Init the weight of the 1-v graph
         WLSI=1.D0

         TPREVAR=.FALSE.
!C.. This is where we will store the MP2 energy
         DO I=2,I_VMAX
            MP2E(I)=0.D0
            
         ENDDO
!C         CALL N_MEMORY_CHECK()
!C         CALL WRITEDET(6,NI,NEL,.TRUE.)
         CALL CLEARCLASSPATHS()
         IF(NTAY(1).LE.0) THEN
            TSYM=.FALSE.
         ELSE
            TSYM=.TRUE.
         ENDIF
         IF(I_VMAX.EQ.0.AND.I_HMAX.EQ.0) THEN
            CALL WIRHODIAG(NI,BETA,I_P,NEL,                             &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,     &
     &         RHOEPS,ILOGGING,TSYM,ECORE,WLRI,WLSI)                    
            RETURN
         ENDIF
!C         IF(I_HMAX.EQ.-2) THEN
!C            CALL WIRHODIAG_SUB(NI,BETA,I_P,NEL,
!C     &         NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
!C     &         RHOEPS,ILOGGING,TSYM,ECORE,WLRI,WLSI,I_VMAX)
!C            RETURN
!C         ENDIF
         IF(I_HMAX.EQ.-3.OR.I_HMAX.EQ.-4.OR.I_HMAX.EQ.-7.OR.(I_HMAX.LE.-12.AND.I_HMAX.GT.-20)) THEN
!C.. Importance sampling monte-carlo
!C.. H=-3 is unbiased, and H=-4 is biased
            CALL MCPATHSR4(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,  &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,     &
     &               LSTE,ICE,RIJLIST,NWHTAY(1,1),ILOGGING,ECORE,ILMAX,      &
     &               WLRI,WLSI,DBETA,DLWDB,0)
            CALL WRITECLASSPATHS()
            RETURN
         ELSEIF(I_HMAX.EQ.-5) THEN
!C.. Markov Chain Monte Carlo
            CALL MCPATHSR5(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,  &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,     &
     &               LSTE,ICE,RIJLIST,NWHTAY(1,1),ILOGGING,ECORE,ILMAX,      &
     &               WLRI,WLSI)
            CALL WRITECLASSPATHS()
            RETURN
         ELSEIF(I_HMAX.EQ.-6) THEN
!C.. Pick the largest cluster
            CALL MCPATHSR6(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,  &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,     &
     &               LSTE,ICE,RIJLIST,NWHTAY(1,1),ILOGGING,ECORE,ILMAX,      &
     &               WLRI,WLSI,DBETA,DLWDB)
            RETURN
         ELSEIF(I_HMAX.EQ.-10) THEN
!C.. Use a different method at each vertex level
!            WRITE(6,*) "Do we go in here?"
            CALL MCPATHSR10(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1, &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,     &
     &               LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,      &
     &               WLRI,WLSI,DBETA,DLWDB)
            RETURN
         ENDIF
         IF(I_VMAX.EQ.0) THEN
!C.. Pick graphs which may not contain distinct vertices 
            CALL MCPATHSR2(NI,BETA,I_P,I_HMAX,NEL,NBASISMAX,G1,         &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,     &
     &               LSTE,ICE,NWHTAY(1,1),ILOGGING,ECORE,ILMAX,WLRI,WLSI,    &
     &               TSYM)
            RETURN
         ENDIF
         proc_timer%timer_name='MCPATHSR3 '
         call set_timer(proc_timer)
         OTIME=etime(tarr)
         TLOG=BTEST(ILOGGING,1)
         IF(TLOG) THEN
            OPEN(11,FILE="MCPATHS",STATUS="OLD",POSITION='APPEND')
!C.. go to end of file
!            I=FSEEK(11,0,2)
            call write_det (11, NI, .true.)
         ENDIF
         IF(BTEST(ILOGGING,0)) THEN
            OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
            IF(TMPTHEORY) OPEN(13,FILE="MP2PATHS",STATUS="UNKNOWN")
         ENDIF
!C.. Set the first node to I_I
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         HIJS(0) = get_helement (nI, nI, 0)
!C     &            HIJS(0))
         RHOII(0)=RH
         RHOIJ(0,0)=RH
         WLRI=LOG(RHOII(0))
         TOTAL=1.D0
         NTOTAL=1.D0
         DLWDB=HIJS(0)
         H00=HIJS(0)
         IF(TLOG) WRITE(11,"(I12,2G25.16,F19.7,2I12,F19.7)") 1,TOTAL,TOTAL,0.D0,1,1,H00
!c         WRITE(6,*) 0,TOTAL,TOTAL,0
!C.. I_V is the number of vertices in the path
         DO I_V=2,I_VMAX
!CI_P
!c            WRITE(6,"I2") I_V
            WRITE(STR,"(A,I5)") "FMCPR",I_V
            proc_timer2%timer_name=STR
            call set_timer(proc_timer2)
            L=0
            LT=0
            BTABLE(0)=-1
            DLWDB2=0.D0
            IF(I_HMAX.EQ.-1) THEN
!C.. This code generates list of all excitations at each level
            WRITE(6,"(A,I3,A)") " Full Sum Old ",I_V," Vertex"
            F(I_V)=FMCPR3(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,G1,NBASIS,   &
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,                 &
     &            0,RHOII,RHOIJ,LSTE,ICE,RIJLIST,                          &
     &            NWHTAY(1,1),L,LT,I_HMAX,NLIST,LSTP,                           &
     &            BTABLE,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB2,HIJS,      &
     &            MP2E,NTOTAL)
            ELSEIF(I_HMAX.EQ.-8) THEN
!C.. This code generates excitations on the fly
            WRITE(6,"(A,I3,A)") " Full Sum New ",I_V," Vertex"

               EREF=DLWDB/TOTAL
            
               F(I_V)=FMCPR3B(NI,BETA,I_P,IPATH,I_V,NEL,                   &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,        &
     &        RHOEPS,0,RHOII,RHOIJ,NWHTAY(1,1),I_HMAX,LOCTAB,                   &
     &         ILOGGING,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,       &
     &         MP2E,NTOTAL,EREF,VARSUM,TOTAL)
            ELSEIF(I_HMAX.EQ.-20) THEN 
!C.. This code generates excitations on the fly, and diagonalizes a matrix of
!C.. HIJ rather than RHOIJ
            WRITE(6,"(A,I3,A)") " Full Sum HDiag ",I_V," Vertex"

            EREF=DLWDB/TOTAL
            
               F(I_V)=FMCPR3B2(NI,BETA,I_P,IPATH,I_V,NEL,                  &
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,        &
     &        RHOEPS,0,RHOIJ,NWHTAY(1,1),I_HMAX,LOCTAB,                         &
     &         ILOGGING,TSYM,ECORE,DBETA,DLWDB2,HIJS,L,LT,IFRZ,1.D0,       &
     &        MP2E,NTOTAL,I_VMAX,EREF,VARSUM,TOTAL)
            ELSEIF(I_HMAX.EQ.-9) THEN
!C.. This code generates a star consisting of all the two-vertex terms,
!C.. forcing them to be disconnected.
               F(I_V)=FMCPR3STAR(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,      &
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,                 &
     &            LSTE,ICE,RIJLIST,L,LT,                                   &
     &            NWHTAY(1,1),ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB2)
            ELSEIF(I_HMAX.EQ.-21) THEN
               IF(TDIAGNODES) THEN
! This code prediagonalises connected nodes of virtual orbitals, before solving as a star.
                F(I_V)=fMCPR3StarNodes(NI,BETA,I_P,NEL,NBASISMAX,G1,       &
     &              NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,        &
     &              L,LT,NWHTAY(1,1),ILOGGING,TSYM,ECORE,DBETA,DLWDB2)       
               ELSEIF(TGraphMorph) THEN
!This involves iterativly improving a graph in order to maximise the connections between determinants.
                    IF(I_VMAX.ne.2) STOP 'Error in vertex number'
                    CALL MorphGraph(F(I_V),DLWDB2)
               ELSEIF(TStarTrips) THEN
!Excited stars of triples from all doubles are prediagonalised
                    CALL StarDiagTrips(DLWDB2,F(I_V))
               ELSEIF(TFCIMC) THEN
!A MC simulation involving replicating particles is run
#ifdef PARALLEL
!                    WRITE(6,*) "Get Here?!"
                    CALL FciMCPar(F(I_V),DLWDB2)
!                    WRITE(6,*) "Get Here??!"
#else
                    CALL FciMC(F(I_V),DLWDB2)
#endif
               ELSEIF(TReturnPathMC) THEN
!A MC simulation involving replicating particles, constrained to returning paths is run
                    CALL ReturnPathMC(F(I_V),DLWDB2)
               ELSE
!C.. This code generates a star consisting of all the two-vertex terms,
!C.. forcing them to be disconnected. - this uses new excitation generators
               F(I_V)=FMCPR3STARNewExcit(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS, &
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,L,LT,               &
     &            NWHTAY(1,1),ILOGGING,TSYM,ECORE,DBETA,DLWDB2,MP2E)               
               ENDIF
            ENDIF
            call halt_timer(proc_timer2)
            NTIME=etime(tarr)
            TOTAL=TOTAL+F(I_V)
            DLWDB=DLWDB+DLWDB2
!c            WRITE(6,*) I_V,F(I_V),TOTAL,get_total_time(proc_timer2),L,LT
            IF(TLOG) THEN
               WRITE(11,"(I12,3G25.16,2I12,G25.16)",advance='no')             &
     &            I_V,F(I_V),TOTAL,NTIME-OTIME,L,LT,DLWDB2
               IF(I_HMAX.EQ.-20) THEN
                  DO I=2,I_VMAX
                     WRITE(11,"(G25.16)",advance='no') MP2E(I)
                  ENDDO
               ENDIF
               WRITE(11,*)
               CALL FLUSH(11)
            ENDIF
            IF(ISNAN(TOTAL)) THEN
!C.. save all log files
               ITIME=etime(tarr)
               CALL FLUSH(11)
!               CALL LOGNAN(NI,NEL,BETA,ITIME)
               WRITE(6,*) "WARNING: nan found at time",ITIME
               WRITE(6,"(A)",advance='no') "  nan det="
               call write_det (6, NI, .true.)
            ENDIF
            CALL WRITECLASSPATHS()
         ENDDO
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,0)) THEN
            CLOSE(10)
            IF(TMPTHEORY) CLOSE(13)
         ENDIF
         if(calcp_logweight) then
!  For LogWeight, we are passing around and summing log x[G] in the weights.  The total weight, w_i=prod_G x[G]
!  We return the log of the weight, which is simply the sum of log x[G] which is in total
            wlsi=total
         else
            WLSI=LOG(TOTAL)
!C.. DLWDB contains <D|H exp(-b H)|D>/rhoii^P and we want 
!C.<D|H exp(-b H)|D> / <D|exp(-b H)|D>
!C.. <D|exp(-b H)|D>=exp(P WLRI + WLSI)
!C         WRITE(6,*) DLWDB,TOTAL,EXP(WLRI), HIJS(0)
            DLWDB=DLWDB/TOTAL
         endif
         TOTAL=HIJS(0)
! Only print out MP energies if I_VMAX>1 [MP2E(2:I_VMAX)]
         IF (I_VMAX.gt.1) THEN 
             IF(abs(MP2E(2)).gt.0.D0) THEN
                FF=HIJS(0)
                DO I=2,I_VMAX
                   FF=FF+MP2E(I)
                   WRITE(STR,"(A2,I1,A10)") "MP",I," ENERGY = "
                   WRITE(6,"(A14,G25.16)") STR,FF
                ENDDO
             ENDIF
         ENDIF
         call halt_timer(proc_timer)
         CALL N_MEMORY_CHECK()
         RETURN
      END SUBROUTINE

!C.. A function to loop recursively over each node set choosing a different
!C.. node for each set.  All nodes are distinct.  Paths IJIKJI etc.
!C.. are generated by permutation from IJKI, and summed up to length I_HMAX
!C.. using the appropriate weightings (Z-sums) from CALCPATH7.(26/01/04). 
      RECURSIVE FUNCTION FMCPR3(NI,BETA,I_P,IPATH,I_V,NEL,                 &
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,             &
     &   RHOEPS,I_VIND,RHOII,RHOIJ,LSTE,ICE,RIJLIST,NWHTAY,L,LT,I_HMAX,    &
     &   NLIST,LSTP,BTABLE,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB,HIJS,     &
     &   MP2E,NTOTAL)  RESULT (FMCPR3RES)
         use constants, only: dp
         use SystemData, only: BasisFN
         Use Determinants, only: get_helement
         USE Logging , only : G_VMC_LOGCOUNT
         use CalcData , only : TMPTHEORY
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         TYPE(BasisFN) G1(*),ISYM
         real(dp) FMCPR3RES
         INTEGER I_V,NEL,I_P,nBasisMax(5,*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY(2),I_VIND,NWHTAY,L,LT,ILOGGING,J,K,I_VMAX
         INTEGER I,IVLMAX,II
         REAL*8 ALAT(*),ECORE
         COMPLEX*16 FCK(*)
         real(dp) CALCPATHS_N
         HElement_t UMat(*),R
         INTEGER IPATH(NEL,0:I_V)
         real(dp) RHOCUML,TOTAL,RHOII(0:I_V)
         HElement_t RHOIJ(0:I_V,0:I_V)
         INTEGER INODE(NEL),ILMAX
!C.. LSTE is a list of excitations (which we will generate)
!C.. ICE is the IC of each excitation (i.e. how much it differs from us (INODE)
!C.. LSTP is an index into LSTE of what has been processed in it so far
!C..      (i.e. it is the first un-processed node)
!C.. NLIST contains the legnth of LSTE for each node
         INTEGER LSTE(NEL,0:ILMAX,0:I_V-1),NI(NEL)
         INTEGER LSTP(0:I_V),LSTP2(0:I_V),NLIST(0:I_V)
         INTEGER ICE(0:ILMAX,0:I_V-1)
         HElement_t RIJLIST(0:ILMAX,0:I_V-1)
         INTEGER IVLEVEL,I_HMAX,IC,BTABLE(0:I_V)
         REAL*8 BETA,RHOEPS,DBETA
         LOGICAL LSAME,TSYM
         LOGICAL TLOG,TLOG2,TLOG3,TLOG4,TLOG5
         real(dp) DLWDB,DLWDB2
         HElement_t HIJS(0:I_V),RH
         INTEGER ICLS
         real(dp) MP2E(:),NTOTAL
         INTEGER EXFLAG
!C.. Allow both single and double excitations
         EXFLAG=3
         TLOG=BTEST(ILOGGING,0)
         TLOG5=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,12)
         TLOG4=BTEST(ILOGGING,9)
         TLOG=TLOG.AND..NOT.TLOG4
         LT=LT+1
 
         TOTAL=0.D0
!C.. This is the current node (set by our parent)         
         CALL NECI_ICOPY(NEL,IPATH(1:NEL,I_VIND),1,INODE,1)
!C.. Set the Zeroth node in our excitation list to be ourselves
!C.. as we'll always want to exclude ourselves from being counted again
!C.. (unused)
         CALL NECI_ICOPY(NEL,INODE,1,LSTE(1:NEL,0,I_VIND),1)
         ICE(0,I_VIND)=0
!C.. Find RHO_II for this node
         IF(INODE(1).EQ.0) THEN
            R=1.D0
         ENDIF
         CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(I_VIND)=RH
         RHOIJ(I_VIND,I_VIND)=RH
!C.. we note that if this node has rho=0 (i.e. RH=0) then we just return
         IF(.not.abs(RH).gt.0.D0) THEN
            FMCPR3RES=0.D0
            RETURN
         ENDIF
         IF(I_VIND.EQ.(I_V-1)) THEN
!C.. If we're at the last node we call CALCPATHS to generate all 
!C.. the paths for it
            CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,I_V),1) 
            RHOII(I_V)=RHOII(0)
            IF(TLOG5) THEN
               IF(.NOT.TLOG3) THEN
                  CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
               ELSE
                  CALL WRITEPATHEX(10,IPATH,I_V,NEL,.FALSE.)
               ENDIF
            ENDIF
            
            IF(TLOG2) CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
!C..
            ICLS=0 
            TOTAL=TOTAL+CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,             &
     &         I_P,1.D0,NEL,I_VMAX,ILOGGING,DBETA,DLWDB2,HIJS,ICLS)  
            NTOTAL=NTOTAL+TOTAL
!C.. Sum up the components of <D|H exp(-b H)|D>
            DLWDB=DLWDB+DLWDB2
               IF(TLOG) WRITE(10,"(2E25.16, I7)") TOTAL,DLWDB2,ICLS
!caa     call flush(10) 
            L=L+1
            FMCPR3RES=TOTAL
            IF(I_V.EQ.2.AND.TMPTHEORY) THEN
!C.. Calculate the MP2 Energy
!C.. NMAX has ARR hidden in it.
               CALL ADDMP2E(HIJS,NMAX,NBASIS,IPATH,NEL,TLOG,MP2E)
            ENDIF
            IF(TLOG4.AND.MOD(L,G_VMC_LOGCOUNT).EQ.0) THEN
!C.. log every 1000
               IF(TMPTHEORY) THEN
                  WRITE(10,*) L,NTOTAL,DLWDB,MP2E
               ELSE
                  WRITE(10,*) L,NTOTAL,DLWDB
               ENDIF
               CALL FLUSH(10)
            ENDIF
            RETURN
         ENDIF
!C.. Find all nodes connected to this node
!C.. This is of maximum order NTAY*2
         NLIST(I_VIND)=ILMAX
!C         CALL GETSYM(INODE,NEL,G1,NBASISMAX,ISYM)
!C         CALL GENSYMEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,.TRUE.,ISYM,
!C     &            LSTE(1,1,I_VIND),ICE(1,I_VIND),NLIST(I_VIND))
         CALL GENEXCIT(INODE,ABS(NTAY(1)*2),NBASIS,NEL,                 &
     &         LSTE(1,1,I_VIND),ICE(1,I_VIND),NLIST(I_VIND),1,G1,TSYM,  &
     &         NBASISMAX,.FALSE.)
!C         WRITE(37,*) I_VIND,NLIST(I_VIND)
         IF(NLIST(I_VIND).GT.(NBASIS-NEL)**2*NEL*NEL) THEN
            WRITE(6,*) "WARNING on excitations"
         ENDIF
!C.. First we need to check that all the nodes in our list are not
!C.. members of the excitation lists of nodes further back in the line
!C.. which have already been processed, or which are yet to be processed
!C.. i.e. if they're already in a list somewhere, then we don't need to
!C.. include them, as they will be processed later.
!C.. As the lists are ordered, this is quite a small job.

!C.. 20050307 - it seems that a better way to do this (which does not
!C.. require that the lists be sorted, and thus rather frees up the 
!C.. excitation generation routine), is to see whether each excitation
!C.. is connected to a previous vertex in the current list.  If it is,
!C.. and we're attempting to connect it to a new vertex, then we should
!C.. remove it from our list.
!C.. This procedure need not be done here, but can be done as we recurse
!C.. through the excitations.
         CALL ELIMDUPS(LSTE,LSTP,I_VIND,NEL,NLIST,ILMAX,NI)

!C.. We now find all the nodes which actually have a connection to us
!C.. eliminating any others in the list.  As this requires that we calc
!C.. RHOIJ, so we store it too.
!C.. II points to the next free position in the list, as we move along,
!C.. rewriting the list to remove the nodes which are zero, and to which
!C.. we are not connected.  
         II=1
         DO I=1,NLIST(I_VIND)
            IF(LSTE(1,I,I_VIND).NE.0) THEN
               CALL CALCRHO2(INODE,LSTE(1:NEL,I,I_VIND),BETA,I_P,NEL,      &
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,         &
     &            RH,NTAY,ICE(I,I_VIND),ECORE)  
               IF(abs(RH).GE.RHOEPS) THEN
                  IF(II.NE.I) CALL NECI_ICOPY(NEL,LSTE(1:NEL,I,I_VIND),1,  &
     &                    LSTE(1:NEL,II,I_VIND),1)
                  ICE(II,I_VIND)=ICE(I,I_VIND)
                  RIJLIST(II,I_VIND)=RH
                  II=II+1
               ENDIF
            ENDIF
         ENDDO
         NLIST(I_VIND)=II-1            
            
!C         WRITE(6,*) NLIST(I_VIND)
!C         DO I=1,NLIST(I_VIND)
!C            CALL WRITEDET(6,LSTE(1,I,I_VIND),NEL,.TRUE.)
!C         ENDDO
!C.. We recurse over all possibilities for the next node:
!C..  1) Nodes connected to us (excluding those connected to previous nodes
!C..        (whether we have processed them or not))
!C..        Once this has been done, our list of connected nodes is marked as
!C..        entirely processed.
!C..  2) Nodes connected to previous nodes which have not yet been processed
!C..        (excluding those which are connected to further previous nodes. 
!C..         NB: These may be connected to us, but will have been excluded
!C..         in 1 as they are connected to previous nodes to 1
!C..  3) Further recursions down the (direct) line from us back to I.  Previous 
!C..        completed lines will have been entirely processed, so we don't
!C..        need to attempt to add nodes to then (but we do need to ensure
!C..        that new nodes don't connect to them)
!C.. 
!C..   All of this boils down to looking at all previous nodes (irrespective
!C..   of which line they're in), and attempting to attach nodes to the
!C..   positions after the parts which have been processed.  This requires
!C..   a COPY of the processed index table upon which to work, as we will 
!C..   be modifying that (i.e. processing nodes after this one or our
!C..   predecessors), GIVEN we have connected up this node.  We will then
!C..   wish to return to past nodes, and connect up further ones (past
!C..   and instead of this one)
        LSTP(I_VIND)=1
        IF(NWHTAY.EQ.0) THEN
!C.. allow chain graphs
         IVLMAX=I_VIND
        ELSE
         IVLMAX=0
        ENDIF
         IVLEVEL=IVLMAX
        DO WHILE(IVLEVEL.GE.0)
!C.. we first take a copy of the pointer table
         DO II=0,I_VIND
            LSTP2(II)=LSTP(II)
         ENDDO
!C.. We must recurse over all the nodes connected to the node at IVLEVEL
         I=LSTP(IVLEVEL)
!C.. remind us which node we're doing
         IF(IVLEVEL.LT.I_VIND) THEN
            CALL NECI_ICOPY(NEL,IPATH(1:NEL,IVLEVEL),1,INODE,1)
         ENDIF
         DO WHILE (I.LE.NLIST(IVLEVEL))
!C.. Tell future recursions that we've processed this node
            LSTP2(IVLEVEL)=I+1
!C.. We also need to update this node in the table we use as the master table
!C.. This is ok as we're not modifying another vertex level's space, and 
!C.. necessary to tell future copies of this master how far we've processed
            IF(IVLEVEL.EQ.I_VIND) LSTP(I_VIND)=I+1
!C.. Because of the elimination of duplicates, we need to ensure that this
!C.. node in the list hasn't been set to 0 (i.e. it was a duplicate, and has
!C.. been removed)
            IF(LSTE(1,I,IVLEVEL).NE.0) THEN
!C.. For VLEVEL=I_VIND, the list at this level has had removed all nodes
!C.. which are connected to previous nodes.
!C.. For VLEVEL<I_VIND, there may be some unprocessed nodes connected to the 
!C.. node at VLEVEL, which are also connected to a node at a higher VLEVEL
!C.. (e.g. this node).  We DO need to process these, because we specifically
!C.  removed them from the list of nodes connected to (e.g.) this level 
!C.. for processing later (i.e. now).
!C.. Calculate RHO_IJ - the connection between the new node and this node
!C.. LSTE(1:NEL,I,IVLEVEL) is the node we're attempting to connect to
!               WRITE(6,"(A)",advance='no') "MCPR"
!               CALL WRITEDET(6,LSTE(1,I,IVLEVEL),NEL,.FALSE.)
               RH=RIJLIST(I,IVLEVEL)
!               WRITE(6,*) RH
!C.. see if it's worth going further
               RHOIJ(IVLEVEL,I_VIND+1)=RH
#ifdef __CMPLX
               RHOIJ(I_VIND+1,IVLEVEL)=DCONJG(RH)
#else
               RHOIJ(I_VIND+1,IVLEVEL)=(RH)
#endif
               IF(.NOT.(abs(RH).GE.RHOEPS)) THEN
!C.. we're not connected to this node, so we remove it from our list
                  DO II=1,NEL
                     LSTE(II,I,IVLEVEL)=0
                  ENDDO
               ELSE
!C.. we now calculate the RHOIJ for the previous nodes to connect to the new one
                  DO II=0,I_VIND           
!C.. we've already calculated the RHO to the node to which we are connected.
                        CALL CALCRHO2(IPATH(1:NEL,II), LSTE(1:NEL,I,IVLEVEL),BETA,I_P, &
     &                   NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,                    &
     &                   ALAT,UMAT,RH,NTAY,-1,ECORE)
                        IF(.NOT.(abs(RH).GE.RHOEPS)) THEN
!C.. These two vertices are not connected, which is as we want
                           RH=0.D0
                        ENDIF
!C                        ELSE
!C.. These two vertices are connected

!C***********
!C      STOP 'NB - LOOK AT RED GRAPH PIC'


#ifdef __CMPLX
                        RHOIJ(I_VIND+1,II)=DCONJG(RH)
#else
                        RHOIJ(I_VIND+1,II)=(RH)
#endif
                        RHOIJ(II,I_VIND+1)=RH
                  ENDDO
!C.. 20050307 - At this point we note that if a vertex is connected to
!C.. any vertices before it in the path (excluding the one we excited 
!C.. it from), we need to ignore this vertex, as it will have already been
!C.. counted before.

!C.. Add the new node to the path, and then recurse from it
                  CALL NECI_ICOPY(NEL,LSTE(1:NEL,I,IVLEVEL),1,IPATH(1:NEL,I_VIND+1),1)
!C.. Get the H element so we can calculate the energy
                  HIJS(I_VIND+1) = get_helement (iPath(:,0), iPath(:,I_VIND+1))
                  BTABLE(I_VIND+1)=IVLEVEL
                  TOTAL=TOTAL+                                                      &
     &                     FMCPR3(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,              &
     &                     G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,                   &
     &                     NTAY,RHOEPS,I_VIND+1,RHOII,RHOIJ,LSTE,                   &
     &                     ICE,RIJLIST,NWHTAY,L,LT,I_HMAX,NLIST,                    &
     &                     LSTP2,BTABLE,ILOGGING,TSYM,ECORE,ILMAX,                  &
     &                     DBETA,DLWDB,HIJS,MP2E,NTOTAL)
               ENDIF
            ENDIF
            I=I+1
          ENDDO
!C.. Recurse back through the direct line of descendents
          IVLEVEL=BTABLE(IVLEVEL)
         ENDDO
         FMCPR3RES=TOTAL
         RETURN
      END FUNCTION



!C.. A function to loop recursively over each node set choosing a different
!C.. node for each set.  All nodes are distinct.  Paths IJIKJI etc.
!C.. are generated by permutation from IJKI, and summed up to length I_HMAX
!C.. using the appropriate weightings (Z-sums) from CALCPATH7.(26/01/04).

!C.. This version doesn't need to generate excitation lists, so 
!C.. calculates excitations on the fly.  Overheads are bigger, but
!C.. scaling should be better. (9/3/05)

!C.. Various flags are available in NWHTAY
!C.. Bit 8 (512).  If set the lowest 8 bits of NWHTAY correspond to freezing data.
!C..                  (which is no longer implemented)
!C.. Bit 0 (1).  If 0, then allow connections to anywhere in a graph
!C..             If 1, then only allow connections to the root
!C.. Bit 1 (2).  If 0 allow all topologies.
!C..             If 1, disallow a det if it's already connected to one in the graph
!C.. Bit 3 (8).  If 0 allow all types of excitations
!C..             If 1 only allow singles.
!C.. Bit 4 (16)  If 0 allow all types of excitations
!C.,             If 1, only allow doubles.
         RECURSIVE FUNCTION FMCPR3B(NI,BETA,I_P,IPATH,I_V,NEL,          &
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,                   &
     &   RHOEPS,I_VIND,RHOII,RHOIJ,NWHTAY,I_HMAX,LOCTAB,                         &
     &   ILOGGING,TSYM,ECORE,DBETA,DLWDB,HIJS,L,LT,IFRZ,FSCALE,MP2E,             &
     &   NTOTAL,EREF,VARSUM,WREF) RESULT (FMCPR3BRES)
         use constants, only: dp
         use SystemData, only: BasisFN
         Use Determinants, only: get_helement
         use CalcData , only : TVARCALC,TMPTHEORY
         USE Logging , only : G_VMC_LOGCOUNT
         USE PrecalcData , only : TPREVAR,PREWEIGHTEPS
         use mcpathsdata, only: EGP
         use legacy_data, only: irat
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         TYPE(BasisFN) G1(*),ISYM
         INTEGER I_V,NEL,I_P,nBasisMax(5,*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY(2),I_VIND,NWHTAY,ILOGGING,J,K,II
         INTEGER I
         COMPLEX*16 FCK(*)
         HElement_t UMAT(*),R
         REAL*8 ALAT(*),ECORE,PR
         real(dp) TOTAL,FMCPR3BRES,WREF,EREF,PROB
         real(dp) CALCPATHS_N
         INTEGER IPATH(NEL,0:I_V)
         real(dp) RHOII(0:I_V),DLWDB,DLWDB2
         HElement_t RHOIJ(0:I_V,0:I_V),RH,HIJS(0:I_V)
         INTEGER INODE(NEL)
         INTEGER NI(NEL),NJ(NEL)
         INTEGER I_HMAX
         REAL*8 BETA,RHOEPS
         LOGICAL TSYM
         LOGICAL TLOG,TLOG2,TLOG3,TLOG4,TLOG5,TLOG6
         REAL*8 DBETA
         INTEGER ICLS
         INTEGER,pointer :: NMEM(:)
         INTEGER NMEMLEN
         INTEGER, pointer :: OGEN(:)
         INTEGER, pointer :: CURGEN(:)
!C.. LOCTAB(1)%p is the address of the generator used to create node 1 in
!C.. the path (i.e. J).  LOCTAB(1)%l is the length of the generator (i.e. 
!C.. the amount of memory used to store it)
         TYPE(EGP) LOCTAB(:)
         TYPE(EGP) LOCTAB2(I_V)
         LOGICAL TFAIL,TNEXT,T
         INTEGER L,LT,IVLEVEL,IEXFROM,IVLMAX,IVLMIN
         INTEGER ICMPDETS
         INTEGER IC
         INTEGER DUMMY(0:I_V)
         INTEGER IFRZ(0:NBASIS,I_V),IFRZ2(0:NBASIS)
         INTEGER EX(2,2),ICIL,ICILMAX
         INTEGER STORE(6)
         REAL*8 FSCALE,FSC2
         real(dp) MP2E(:),NTOTAL
         INTEGER EXFLAG
         LOGICAL ISCONNECTEDDET
         REAL*8 VARSUM,SumX,SumY,SumXY,SumXsq,SumYsq,SumP
         DATA SumP/0.D0/
         SAVE SumX,SumY,SumXY,SumXsq,SumYsq,SumP
!         write(6, *) "FMCPR3B I_VIND:",I_VIND
         IF(BTEST(NWHTAY,0).AND.I_VIND.EQ.0) WRITE(6,*) "FMCPR3 ForceRoot."
         IF(BTEST(NWHTAY,1).AND.I_VIND.EQ.0) WRITE(6,*) "FMCPR3 ForceStar."
         SELECT CASE (IAND(NWHTAY,24))
         CASE(0)
!C.. Allow both singles and doubles
            EXFLAG=3
         CASE(8)
!C.. only singles
            IF(I_VIND.EQ.0) WRITE(6,*) "FMCPR3 Only singles."
            EXFLAG=1
         CASE(16)
            IF(I_VIND.EQ.0) WRITE(6,*) "FMCPR3 Only doubles."
            EXFLAG=2
         CASE(24)
            STOP "Invalid combination of flags in NWHTAY"
         END SELECT
         IF(I_VIND.EQ.0) THEN
!C.. 1st time in
            IFRZ(0:NBASIS,1:I_V) =0
            IF(IAND(NWHTAY,8).NE.0) THEN
!C.. Force freexing
               IFRZ(0,1)=NWHTAY
            ELSE
               IFRZ(0,1)=0
            ENDIF
         ENDIF
         IFRZ2(0:NBASIS) =0
         CALL NECI_ICOPY(NBASIS+1,IFRZ(0,I_VIND+1),1,IFRZ2,1)
          
         
!C            DO I=0,NBASIS
!C               WRITE(10,"(I2)",advance='no'),IFRZ2(I)
!C            ENDDO
!C         WRITE(10,*) "V_",I_VIND
!C.. LOCTAB(1)%p is the address of the generator used to create node 1 in
!C.. the path (i.e. J).  LOCTAB(1)%l is the length of the generator (i.e. 
!C.. the amount of memory used to store it)
         TLOG=BTEST(ILOGGING,0)
         TLOG6=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,12)
         TLOG4=BTEST(ILOGGING,9)
         TLOG5=BTEST(ILOGGING,6)
         TLOG=TLOG.AND..NOT.TLOG4
         TOTAL=0.D0
         LT=LT+1
!C.. This is the current node (set by our parent)         
         CALL NECI_ICOPY(NEL,IPATH(1:NEL,I_VIND),1,INODE,1)
         CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,NBASISMAX,                       &
     &   G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(I_VIND)=RH
         RHOIJ(I_VIND,I_VIND)=RH
!C.. we note that if this node has rho_II=0 (i.e. RH=0) then we just return
         IF(.not.abs(RH).gt.0.D0) THEN
            FMCPR3BRES=0.D0
            RETURN
         ENDIF
         IF(I_VIND.EQ.(I_V-1)) THEN
!C.. If we're at the last node we call CALCPATHS to generate all 
!C.. the paths for it
            CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,I_V),1) 
            RHOII(I_V)=RHOII(0)
!            CALL WRITEPATH(6,IPATH,I_VIND+1,NEL,.TRUE.)
!            DO J=1,I_VIND
!C            DO I=0,NBASIS
!C               WRITE(10,"(I2)",advance='no'),IFRZ(I,J)
!C            ENDDO
!              WRITE(6,*) J,":",loc(loctab(J)%p),LOCTAB(J)%v
!            ENDDO
            IF(TLOG6) THEN
               IF(.NOT.TLOG3) THEN
                  CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
               ELSE
                  CALL WRITEPATHEX(10,IPATH,I_V,NEL,.FALSE.)
               ENDIF
            ENDIF
            IF(TLOG2) CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
!C.. 
            ICLS=0
            TOTAL=TOTAL+ CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,                      &
     &         I_P,FSCALE,NEL,I_V,ILOGGING,DBETA,DLWDB2,HIJS,ICLS)      
!C.. Sum up the components of <D|H exp(-b H)|D>
            DLWDB=DLWDB+DLWDB2
            NTOTAL=NTOTAL+TOTAL
            IF(TLOG5) THEN
!  Log XIJS (usually for debugging), and the pgen
!  NMAX has Arr hidden in it
              CALL CalcWriteGraphPGen(10,IPATH,I_V,nEl,LOCTAB,G1,             &
     &            nBasisMax,UMat,NMAX,nBasis,PR,DUMMY)
                  WRITE(10,"(3E25.16, I7)") TOTAL,PR,DLWDB2,ICLS
            ELSE
               IF(TLOG) WRITE(10,"(2E25.16, I7)") TOTAL,DLWDB2,ICLS
            ENDIF
            L=L+1
            FMCPR3BRES=TOTAL
            IF(I_V.EQ.2.AND.TMPTHEORY) THEN
!C.. Calculate the MP2 Energy
!C.. NMAX has ARR hidden in it.
               CALL ADDMP2E(HIJS,NMAX,NBASIS,IPATH,NEL,TLOG,MP2E)
            ENDIF
            IF(TLOG4.AND.MOD(L,G_VMC_LOGCOUNT).EQ.0) THEN
!C.. log every 1000
               IF(TMPTHEORY) THEN
                  WRITE(10,"(I10,3E25.16)") L,NTOTAL,DLWDB,MP2E
               ELSE
                  WRITE(10,"(I10,2E25.16)") L,NTOTAL,DLWDB
               ENDIF
               CALL FLUSH(10)
            ENDIF
!C            WRITE(10,*) (LOCTAB(I)%v,I=1,I_V-1)
!C            WRITE(6,*) "X"
        
            IF (TVARCALC(I_V).or.TPREVAR) THEN
                J=0
                IF(TOTAL.ge.PREWEIGHTEPS) THEN

                   CALL  CalcWriteGraphPGen(J,IPATH,I_V,nEl,LOCTAB,G1,           &
     &                       nBasisMax,UMat,NMAX,nBasis,Prob,DUMMY)
                   SumX  =SumX   + DLWDB2-(EREF*TOTAL)
                   SumY  =SumY   + TOTAL
                   SumXsq=SumXsq + (DLWDB2-(EREF*TOTAL))**2/Prob
                   SumYsq=SumYsq + ((TOTAL)**2)/Prob
                   SumXY =SumXY  + (DLWDB2-(EREF*TOTAL))*TOTAL/Prob
                   SumP=SumP+Prob
     
                 ENDIF
             ENDIF
            RETURN
        ENDIF
!C.. We recurse over all possibilities for the next node:
!C..  1) Nodes connected to us (excluding those connected to previous nodes)
!C..     We perform this exclusion by looking at whether a newly
!C..     generated node connected to us is connected to a previous node
!C..     which existed before the one we were generated from.  (If we are
!C..     connected to a node created after the one we were created from,
!C..     we have been specifically excluded from that node's excitation
!C..     list because we were connected to a prior node, so we, now being
!C...    generated from that prior node, now need to be counted)
!C..  2) Nodes connected to previous (i.e. before the one we're connected to )
!C..         nodes which have not yet been processed
!C..        (excluding those which are connected to further previous nodes. 
!C..         NB: These may be connected to us, but will have been excluded
!C..         in 1 as they are connected to previous nodes to 1
!C..  3) Further recursions down the (direct) line from us back to I.  Previous 
!C..        completed lines will have been entirely processed, so we don't
!C..        need to attempt to add nodes to them (but we do need to ensure
!C..        that new nodes don't connect to them)
!C.. 
!C..   All of this boils down to looking at all previous nodes (irrespective
!C..   of which line they're in), and attempting to attach nodes to the
!C..   positions after the parts which have been processed.

!C.. We start by creating an excitation generator connected to the
!C.. current node

!C.. Initialiaze the excitation generators
         CALL GETSYM(INODE,NEL,G1,NBASISMAX,ISYM)
         STORE(1)=0
         CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,                      &
     &         .TRUE.,NMEMLEN,NJ,IC,IFRZ(0,I_VIND+1),STORE,EXFLAG)
!C         CALL GENSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,.TRUE.,ISYM,
!C     &         .TRUE.,NMEMLEN,NJ,IC,IFRZ(0,I_VIND+1))
         allocate(NMEM(NMEMLEN))
!C         WRITE(6,"(A,I)") "NMEM",NMEMLEN
!         write(6,*) "alloc NMEM", loc(NMEM)
         NMEM(1)=0
         CALL GENSYMEXCITIT2(INODE,NEL,G1,NBASIS,NBASISMAX,                      &
     &         .TRUE.,NMEM,NJ,IC,IFRZ(0,I_VIND+1),STORE,EXFLAG)
!C.. I_VIND is the node that has just been chosen (so LOCTAB(I_VIND) is a
!C.. generator for that node.  We need to generate from that node, so we
!C.. store at I_VIND+1

!C.. We now choose the next node, attempting to generate from all the
!C.. nodes in the direct line back to I

!C.. We do this by taking a copy of LOCTAB (excluding the generator
!C.. from this node, which we haven't put in yet)
         IF(I_VIND>0) LOCTAB2(1:I_VIND)=LOCTAB(1:I_VIND)
!!#if defined(POINTER8)
!!         CALL NECI_ICOPY(I_VIND*3*2,LOCTAB,1,LOCTAB2,1)
!!#else
!!         CALL NECI_ICOPY(I_VIND*3,LOCTAB,1,LOCTAB2,1)
!!#endif
!C.. Now start at this node, and work backwards along the path to I.
!C.. (This is equivalent to working through all previous nodes in the
!C.. path, as those which are not on our direct route back to I will have 
!C.. their generators exhausted, and thus a zero-entry in the LOCTAB.
         
         IF(IAND(NWHTAY,1).EQ.0) THEN
!C.. Either a normal graph sum or a chain/star sum

!C.. start the iterator at "us" (i.e. get all excitations connected to
!C.. us) and then work back down the direct tree to the root.
!C.. allow chain graphs
            IVLMAX=I_VIND+1
            IVLMIN=0
         ELSE
!C.. only stars start at what remains from the iterator the last one
!C.. was generated from (or just the root if there is none) 
            IVLMAX=I_VIND
            IF(IVLMAX.LT.1) IVLMAX=1
            IVLMIN=IVLMAX-1
         ENDIF
         IVLEVEL=IVLMAX

!LOCTAB contains
!C.. LOCTAB(1,V) is the address of the generator used to create node V in
!C.. the path (the first node is node 0).
!C.. LOCTAB(2,V) is the length of the generator (i.e. the amount of memory used to store it)
!C.. LOCTAB(3,V) is the vertex from which vertex V was excited
!C.   The last node never has an excitation generator, as it never needs one in the full graph generation scheme.


!C.. Set these just in case
         CURGEN=>NMEM
         LOCTAB(I_VIND+1)%l=NMEMLEN
         LOCTAB(I_VIND+1)%v=IVLEVEL-1
         
         DO WHILE (IVLEVEL.GT.IVLMIN)
!C.. If we're at I_VIND+1, then we don't need to create a new copy of the 
!C.. generator, as no other sub-recursion is using it yet
            TNEXT=.FALSE.
            CALL NECI_ICOPY(NBASIS+1,IFRZ(0,IVLEVEL),1,IFRZ2,1)
            IF(IVLEVEL.LE.I_VIND) THEN
!C.. create a copy
               nullify(CURGEN)
               allocate(CURGEN(LOCTAB2(IVLEVEL)%l))
!               write(6,*) "alloc CURGEN for level ",IVLEVEL, loc(CURGEN)
               OGEN=>LOCTAB2(IVLEVEL)%p
               CURGEN(1:LOCTAB2(IVLEVEL)%l)=OGEN(1:LOCTAB2(IVLEVEL)%l)
!!               CALL NECI_ICOPY(LOCTAB2(2,IVLEVEL),OGEN,1,CURGEN,1)
!To stop optimization problem with pgi6.1-6
               I_VIND=I_VIND+1
               I_VIND=I_VIND-1
               IF(IAND(IFRZ(0,I_VIND+1),8).NE.0) THEN
!C.. We need to reset the generator if we're only generating stars.
!C.. Because we now have some frozen orbitals, the original excitation
!C.. will not be re-generated.
                  CALL RESETEXIT2(IPATH(1,LOCTAB(IVLEVEL)%v),NEL,G1,          &
     &               NBASIS,NBASISMAX,CURGEN,IFRZ2)
               ENDIF
               LOCTAB2(I_VIND+1)%l=LOCTAB(IVLEVEL)%l
               LOCTAB2(I_VIND+1)%v=LOCTAB(IVLEVEL)%v
            ELSE
               CURGEN=>NMEM
               LOCTAB2(I_VIND+1)%l=NMEMLEN
               LOCTAB2(I_VIND+1)%v=IVLEVEL-1
            ENDIF
            IEXFROM=LOCTAB2(I_VIND+1)%v
            LOCTAB2(I_VIND+1)%p=>CURGEN
!            write(6,*) "goes to LOCTAB2(",I_VIND+1,")%p", loc(LOCTAB2(I_VIND+1)%p)
            DO WHILE(.NOT.TNEXT)
!C               CALL EXCIT_DUMP(10,IPATH(1,IEXFROM),NEL,G1,NBASIS,
!C     &         NBASISMAX,CURGEN,IFRZ2)
!C.. Now use the generator to make the next node,NJ
               CALL GENSYMEXCITIT2(IPATH(1,IEXFROM),NEL,G1,NBASIS,NBASISMAX,     &
     &         .FALSE.,CURGEN,NJ,IC,IFRZ2,STORE,EXFLAG)
!C               CALL GENSYMEXCITIT(IPATH(1,IEXFROM),NEL,G1,
!C     &            NBASIS,NBASISMAX,.TRUE.,
!C     &            ISYM,.FALSE.,CURGEN,NJ,IC,IFRZ2)
!C               CALL WRITEDET(6,IPATH(1,IEXFROM),NEL,.FALSE.)
!C               CALL WRITEDET(6,NJ,NEL,.FALSE.)
!                WRITE(17,*) LT,NJ(:)
!C               WRITE(6,"(3I)") I_VIND,LOC(CURGEN),CURGEN(6)
!C.. Check to see it's actually been generated
               IF(NJ(1).EQ.0) THEN
                  TFAIL=.TRUE.
                  TNEXT=.TRUE.
               ELSE
                  TFAIL=.FALSE.
               ENDIF
!C.. Now make sure it's not already in the path
               IF(.NOT.TFAIL) THEN
                  DO J=0,I_VIND
                     IF(ICMPDETS(NJ,IPATH(1,J),NEL).EQ.0) THEN
                        TFAIL=.TRUE.
                     ENDIF
                  ENDDO
               ENDIF
               IF(.NOT.TFAIL) THEN
!C.. see if we're connected to what we were excited from
                  CALL CALCRHO2(IPATH(1,IEXFROM),NJ,BETA,I_P,                    &
     &                   NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,              &
     &                   ALAT,UMAT,RH,NTAY,IC,ECORE)
                  IF(.NOT.(abs(RH).GE.RHOEPS)) THEN
!C.. if we're not connected to this node, we fail now
                     RH=0.D0
                     TFAIL=.TRUE.
                  ENDIF
#ifdef __CMPLX
                  RHOIJ(I_VIND+1,IEXFROM)=DCONJG(RH)
#else
                  RHOIJ(I_VIND+1,IEXFROM)=(RH)
#endif
                  RHOIJ(IEXFROM,I_VIND+1)=RH
               ENDIF
            
!C.. Work out the connectivity of the new node
               DO II=0,I_VIND
!C.. we don't need to calc the connection to the node we excited from
!C.. because we've just done that
                  IF(.NOT.TFAIL.AND.II.NE.IEXFROM) THEN
                     CALL CALCRHO2(IPATH(1,II),NJ,BETA,I_P,                      &
     &                   NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,              &
     &                   ALAT,UMAT,RH,NTAY,-1,ECORE)
                     IF(.NOT.(abs(RH).GE.RHOEPS)) RH=0.D0
#ifdef __CMPLX
                     RHOIJ(I_VIND+1,II)=DCONJG(RH)
#else
                     RHOIJ(I_VIND+1,II)=(RH)
#endif
                     RHOIJ(II,I_VIND+1)=RH
!LOCTAB contains
!C.. LOCTAB(1,V) is the address of the generator used to create node V in
!C.. the path (the first node is node 0).
!C.. LOCTAB(2,V) is the length of the generator (i.e. the amount of memory used to store it)
!C.. LOCTAB(3,V) is the vertex from which vertex V was excited
!C.   The last node never has an excitation generator, as it never needs one in the full graph generation scheme.
!                     IF(abs(RH ).gt. 0.D0) THEN
                     IF(II<I_VIND) then
!                      WRITE(6,*) "ICD", II
!                      call flush(6)
                      IF(IsConnectedDet(IPATH(1,II),NJ,nEl,LOCTAB(II+1)%p,             &
     &                  G1,nBasisMax,nBasis)) THEN
!if rhoeps is zero then always set TFAIL.  if rhoeps isn't zero, then set TFAIL if Rh isn't zero (i.e. we've included it before)
                       IF(RH.NE.0.D0.OR.RHOEPS.EQ.0.D0) THEN
!C.. If the node to which this node (NJ) is attached was known about at
!C.. the time of node II, we are allowed to have a connection between
!C.. node NJ and II.  Otherwise, this node must not be attached to II, as
!C.. it will've already been counted as one of those attached to II
!    If we excluded it through a rhoeps cutoff, we ARE allowed to connect to it
!                        CALL WRITEDET(6,NJ,NEL,.FALSE.)
!                        CALL WRITEDET(6,IPATH(1,II),NEL,.FALSE.)

!C.. We disallow this connection if IAND(IFRZ2(0),2) to require that there are
!C.. no loops
                        IF(IAND(NWHTAY,2).NE.0.OR.IEXFROM.GT.II) THEN
!C.. we're not allowed to count this node
                           TFAIL=.TRUE.                  
                        ENDIF
                       ENDIF
                      ENDIF
                     endif
!                     WRITE(6,*) RH,TFAIL
                  ENDIF
               ENDDO
               IF(.NOT.TFAIL) THEN
!C.. If we've got a node we're allowed to count
!C.. Deal with freezing orbitals if we have to
!C..
                  ICIL=0
                  ICILMAX=1
                  FSC2=FSCALE
                  DO WHILE(ICIL.LT.ICILMAX)
                   IF(IAND(NWHTAY,8).NE.0.AND.I_VIND.LT.(I_V-2)) THEN
                     IF(ICIL.EQ.0) THEN
                        EX(1,1)=2
                        CALL GETEXCITATION(IPATH(1,IEXFROM),NJ,NEL,EX,T)
!C                        WRITE(10,*) EX(1,1),EX(1,2),EX(2,1),EX(2,2)
                        ICILMAX=2
                        IF(EX(1,2).NE.0) ICILMAX=ICILMAX+1
                     ENDIF
!C.. We freeze the orbitals we're exciting from and to
                     CALL NECI_ICOPY(NBASIS+1,IFRZ2,1,IFRZ(0,I_VIND+1),1)
!C.. freeze full orbitals
!C.. 0) second excited full (if there is one)
!C.. 1) first excited full 
!C.. 2) both excited full (and set to subtract this contrib)
                     IF(ICIL.EQ.ICILMAX-3) THEN
                        IFRZ(EX(1,2),I_VIND+1)=1
                     ELSEIF(ICIL.EQ.ICILMAX-2) THEN
                        IFRZ(EX(1,1),I_VIND+1)=1
                     ELSEIF(ICIL.EQ.ICILMAX-1) THEN
                        IFRZ(EX(1,1),I_VIND+1)=1
                        IF(EX(1,2).NE.0) IFRZ(EX(1,2),I_VIND+1)=1
                        FSC2=-FSC2
                     ENDIF
!C.. freeze the empties
                     IFRZ(EX(2,1),I_VIND+1)=1
                     IF(EX(2,2).NE.0) IFRZ(EX(2,2),I_VIND+1)=1
                     CALL NECI_ICOPY(NBASIS+1,IFRZ(0,I_VIND+1),1,IFRZ(0,I_VIND+2),1)
                   ELSE
                     CALL NECI_ICOPY(NBASIS+1,IFRZ(0,I_VIND+1),1,IFRZ(0,I_VIND+2),1)
                   ENDIF
!C.. Add it to the path, and recurse
                   CALL NECI_ICOPY(NEL,NJ,1,IPATH(1,I_VIND+1),1)
!C.. Get the H element so we can calculate the energy
                   HIJS(I_VIND+1) = get_helement (iPath(:,0), iPath(:,I_VIND+1))
                   
                   TOTAL=TOTAL+FMCPR3B(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,                &
     &                     G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,                   &
     &                     NTAY,RHOEPS,I_VIND+1,RHOII,RHOIJ,NWHTAY,                 &
     &                     I_HMAX,LOCTAB2,ILOGGING,TSYM,ECORE,                      &
     &                     DBETA,DLWDB,HIJS,L,LT,IFRZ,FSC2,MP2E,NTOTAL,             &
     &                     EREF,VARSUM,WREF)
                   ICIL=ICIL+1
                  ENDDO
               ENDIF
            ENDDO
            IF(IVLEVEL.LE.I_VIND) THEN
!C.. free the copy
!               Write(6,*) "Level ",IVLEVEL," of ",I_VIND
!               write(6,*) "dealloc CURGEN", loc(CURGEN)
               deallocate(CURGEN)
               nullify(CURGEN)
!               write(6,*) "nullify LOCTAB2(",IVLEVEL,")%p: ", loc(LOCTAB2(IVLEVEL)%p)
               nullify(LOCTAB2(IVLEVEL)%p)
            ENDIF
!C.. for the next iterator, we want to find the node to which we are
!C.. connected, and continue from after where this node was excited
            IVLEVEL=LOCTAB2(IVLEVEL)%v
         ENDDO
!         Write(6,*) "End of ",I_VIND
!         write(6,*) "dealloc NMEM", loc(NMEM)
         deallocate(NMEM)
!         write(6,*) "nullify LOCTAB(",I_VIND+1,")%p: ", loc(LOCTAB(I_VIND+1)%p)
!         nullify(LOCTAB(I_VIND+1)%p)
         FMCPR3BRES=TOTAL

         If (TVARCALC(I_V).and.(I_VIND.eq.0).and.(.not.TPREVAR)) Then


            SumYsq=SumYsq+(2*WREF*SumY)+(((WREF)**2)*SumP)
            SumY=SumY+(WREF*SumP)
            SumXY=SumXY+(WREF*SumX)
!           write(6,"A, 5G26.15") "Terms in estimated variance:",
!     &        SumX, SumY, SumXsq, SumYsq, SumXY
         write(6,"(2A,I3,A,G23.16)") "Expected MC Variance for ratio ",       &
     &       "vertex level", I_V," is ", (SumX/SumY)**2*(SumXsq/SumX**2       &
     &                      +SumYsq/SumY**2-2*SumXY/(SumX*SumY))           
            write(6,"(2A,G24.16)") "Sum of Pgens of graphs included",         &
     &         " in variance is ", SumP
            SumP=0.D0
            VARSUM=0.D0
            SumX=0.D0
            SumY=0.D0
            SumXsq=0.D0
            SumYsq=0.D0
            SumXY=0.D0
        End If
                          
        IF (TPREVAR.and.I_VIND.eq.0) THEN

           SumYsq=SumYsq+(2*WREF*SumY)+(((WREF)**2)*SumP)
           SumY=SumY+(WREF*SumP)
           SumXY=SumXY+(WREF*SumX)
!           OPEN(44,FILE="HDIAGVARTERMS",STATUS="UNKNOWN")
!           WRITE(43,*) EREF,WREF
          VARSUM=((SumX/SumY)**2)*((SumXsq/(SumX**2))+(SumYsq/(SumY**2))-2*SumXY/(SumX*SumY))
!           WRITE(44,*) g_VMC_ExcitWeights
!            WRITE(44,("6G25.16")) SumX,SumY,(SumXsq-(SumX**2)),
!     &               (SumYsq-(SumY**2)),(SumXY-(SumX*SumY)),VARSUM
           SumP=0.D0
           SumX=0.D0
           SumY=0.D0
           SumXsq=0.D0
           SumYsq=0.D0
           SumXY=0.D0
       ENDIF
!       WRITE(6,*) "Leaving ",I_VIND
!       call flush(6) 
       RETURN
      END FUNCTION

end module mcpaths

      SUBROUTINE WRITEPATH(NUNIT,IPATH,I_V,NEL,LTERM)
         use Determinants, only: write_det
         IMPLICIT NONE
         INTEGER NUNIT,IPATH(NEL,0:I_V),I_V,NEL
         LOGICAL LTERM
         INTEGER J,K
         
         WRITE(NUNIT,"(A)",advance='no') "["
         DO J=0,I_V-1
            call write_det (NUNIT, IPATH(1,J), .false.)
            WRITE(NUNIT,"(A)",advance='no') ","
         ENDDO
         WRITE(NUNIT,"(A)",advance='no') "] "
         IF(LTERM) WRITE(NUNIT,*)
         RETURN
      END

      SUBROUTINE WRITEPATHEX(NUNIT,IPATH,I_V,NEL,LTERM)
         IMPLICIT NONE
         INTEGER NUNIT,IPATH(NEL,0:I_V),I_V,NEL
         LOGICAL LTERM,T
         INTEGER J,K,EX(2,2)
!C.. First determine the excitation
         WRITE(NUNIT,"(A)",advance='no') "["
         DO J=1,I_V-1
            DO K=0,I_V-2
               EX(1,1)=2
!C               WRITE(6,*) J,K
               CALL GETEXCITATION(IPATH(1,K),IPATH(1,J),NEL,EX,T)
!C               WRITE(6,*) J,K
               IF(EX(1,1).GE.0) EXIT
            ENDDO
            IF(EX(1,2).EQ.0) THEN
               WRITE(NUNIT,"(Z1,A,I5,A,I5,A)",advance='no') K,"(",EX(1,1),")->(",EX(2,1),"),"
            ELSE
               WRITE(NUNIT,"(Z1,A,I5,A,I5,A,I5,A,I5,A)",advance='no') K, &
     &         "(",EX(1,1),",",EX(1,2),")->(",EX(2,1),",",EX(2,2),"),"
            ENDIF
         ENDDO
         WRITE(NUNIT,"(A)",advance='no') "]"
         IF(LTERM) WRITE(NUNIT,*)
         CALL FLUSH(NUNIT)
         RETURN
      END

      SUBROUTINE WRITERHOMAT(NUNIT,RHOIJ,I_V,NEL,LTERM)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NUNIT,I_V,NEL
         LOGICAL LTERM
         HElement_t RHOIJ(0:I_V,0:I_V)
         INTEGER J,K
         WRITE(NUNIT,"(A)",advance='no') "("
         DO J=0,I_V-1
            DO K=0,I_V-1
#ifdef __CMPLX
               WRITE(NUNIT,"(2E25.16,A)",advance='no') RHOIJ(J,K),","
#else
               WRITE(NUNIT,"(E25.16,A)",advance='no') RHOIJ(J,K),","
#endif
            ENDDO
            WRITE(NUNIT,"(A)",advance='no') "|"
         ENDDO
         WRITE(NUNIT,"(A)",advance='no') ")"
         IF(LTERM) WRITE(NUNIT,*)
         RETURN
      END
      SUBROUTINE ELIMDUPS(LSTE,LSTP,I_V,NEL,NLIST,NLISTMAX,NI)
         INTEGER NEL,NLISTMAX
         INTEGER LSTE(1:NEL,0:NLISTMAX,0:I_V)
         INTEGER LSTP(0:I_V),NLIST(0:I_V),NI(1:NEL)
         INTEGER I_V,I,J,IND,NMAX,IC,NMAX2,K
!C.. First we check none of our list are NI
         NMAX=NLIST(I_V)
         J=1
         DO WHILE (J.LE.NMAX)
            IF(ICMPDETS(NI,LSTE(1:NEL,J,I_V),NEL).EQ.0) THEN
               DO K=1,NEL
                  LSTE(K,J,I_V)=0
               ENDDO
               J=NMAX
            ENDIF
            J=J+1
         ENDDO
!C.. Work through previous lists, eliminating all elements from our list which
!C.. are in other lists.
         DO I=0,I_V-1
            IND=1
            J=1
            NMAX2=NLIST(I)
!CLSTP(I)-1
            DO WHILE(IND.LE.NMAX.AND.J.LE.NMAX2)
!C.. check to see if we've got a valid element in our new list.  If not, we
!C.. need to move on to the next element
               DO WHILE(IND.LE.NMAX.AND.LSTE(1,IND,I_V).EQ.0)
                  IND=IND+1
               ENDDO
!C.. Now we compare the old list elements with this one, until we find one
!C.. that's >= it, or we run out of elements
               IF(IND.LE.NMAX) THEN
                  IC=ICMPDETS(LSTE(1:NEL,J,I),LSTE(1:NEL,IND,I_V),NEL)
                  DO WHILE(IC.LT.0.AND.J.LT.NMAX2)
                     J=J+1
                    IC=ICMPDETS(LSTE(1:NEL,J,I),LSTE(1:NEL,IND,I_V),NEL)
                  ENDDO
                  IF(IC.EQ.0) THEN
!C.. We've found one in our new list that's the same as the one in an old list
!C.. so we set it to zero
                     DO K=1,NEL
                        LSTE(K,IND,I_V)=0
                     ENDDO
                  ENDIF
!C.. either we've found an old element that's greater than our current new one
!C.. or we've eliminated one, or we've run out of old elements.  Either way, 
!C.. we can go to the next new element
                  IND=IND+1
                  IF(IC.LE.0) J=J+1
               ENDIF
            ENDDO
         ENDDO
         RETURN
      END


!C.. Generate a star of all n-vertex subgraphs - each of which has been
!C.. diagonalised. 
!C.. This calls ADDDIAGSTAR which writes the relevant roots to a file,
!C.. and also includes the removal of subgraphs
      RECURSIVE REAL*8 FUNCTION FMCPR3NVSTAR(NI,BETA,I_P,IPATH,I_V,NEL,    &
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,             &
     &   RHOEPS,I_VIND,RHOII,RHOIJ,LSTE,ICE,RIJLIST,NWHTAY,L,LT,I_HMAX,    &
     &   NLIST,LSTP,BTABLE,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB,HIJS      &
     &   ) RESULT (FMCPR3NVSTARRES)
         
         use SystemData, only: BasisFN
         use Determinants, only: get_helement
         use constants, only: dp
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         type(BasisFN) ISYM
         INTEGER I_V,NEL,I_P,NBASISMAX(*),G1(*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY(2),I_VIND,NWHTAY,L,LT,ILOGGING,J,K,I_VMAX
         INTEGER I,IVLMAX,II
         REAL*8 FCK(*), ALAT(*),UMAT(*),ECORE,TOTAL,R
         REAL*8 CALCPATHS_N
         REAL*8 RHOCUML
         INTEGER IPATH(NEL,0:I_V)
         REAL*8 RHOII(0:I_V)
         REAL*8 RHOIJ(0:I_V,0:I_V)
         INTEGER INODE(NEL),ILMAX
!C.. LSTE is a list of excitations (which we will generate)
!C.. ICE is the IC of each excitation (i.e. how much it differs from us (INODE)
!C.. LSTP is an index into LSTE of what has been processed in it so far
!C..      (i.e. it is the first un-processed node)
!C.. NLIST contains the length of LSTE for each node
         INTEGER LSTE(NEL,0:ILMAX,0:I_V-1),NI(NEL)
         INTEGER LSTP(0:I_V),LSTP2(0:I_V),NLIST(0:I_V)
         INTEGER ICE(0:ILMAX,0:I_V-1)
         HElement_t RIJLIST(0:ILMAX,0:I_V-1)
         INTEGER IVLEVEL,I_HMAX,IC,BTABLE(0:I_V)
         REAL*8 RH,BETA,RHOEPS
         LOGICAL LSAME,TSYM
         LOGICAL TLOG,TLOG2
         REAL*8 DBETA,DLWDB,DLWDB2
         REAL*8 HIJS(0:I_V)
         INTEGER ICLS
         HElement_t :: hel
         IF(I_VIND.EQ.0) THEN
!C.. setup
            OPEN(19,FILE="STARDATA",STATUS="UNKNOWN")
         ENDIF

         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         LT=LT+1
 
         TOTAL=0.D0
!C.. This is the current node (set by our parent)         
         CALL NECI_ICOPY(NEL,IPATH(1:NEL,I_VIND),1,INODE,1)
!C.. Set the Zeroth node in our excitation list to be ourselves
!C.. as we'll always want to exclude ourselves from being counted again
!C.. (unused)
         CALL NECI_ICOPY(NEL,INODE,1,LSTE(1:NEL,0,I_VIND),1)
         ICE(0,I_VIND)=0
!C.. Find RHO_II for this node
         IF(INODE(1).EQ.0) THEN
            R=1.D0
         ENDIF
         CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,NBASISMAX,           &
     &   G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(I_VIND)=RH
         RHOIJ(I_VIND,I_VIND)=RH
!C.. we note that if this node has rho=0 (i.e. RH=0) then we just return
         IF(RH.EQ.0.D0) THEN
            FMCPR3NVSTARRES=0.D0
            RETURN
         ENDIF
         IF(I_VIND.GT.0) THEN
!C.. if we've got a conencted graph of more than one vertex, we send
!C.. it off to ADDDIAGSTAR which writes it to a file.
            CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,I_VIND+1),1) 
            RHOII(I_V)=RHOII(0)
            IF(TLOG)  CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)
            IF(TLOG2) CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
            CALL ADDDIAGSTAR(IPATH,NEL,RHOIJ,I_VIND+1,I_V,ILOGGING,HIJS,1,L)
         ENDIF
         IF(I_VIND.EQ.(I_V-1)) THEN
!C.. If we've made a graph of the right size, we return
            FMCPR3NVSTARRES=0.D0
            RETURN
         ENDIF


!C.. Find all nodes connected to this node
!C.. This is of maximum order NTAY*2
         NLIST(I_VIND)=ILMAX
!C         CALL GETSYM(INODE,NEL,G1,NBASISMAX,ISYM)
!C         CALL GENSYMEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,.TRUE.,ISYM,
!C     &            LSTE(1,1,I_VIND),ICE(1,I_VIND),NLIST(I_VIND))
         CALL GENEXCIT(INODE,ABS(NTAY(1)*2),NBASIS,NEL,                    &
     &         LSTE(1,1,I_VIND),ICE(1,I_VIND),NLIST(I_VIND),1,G1,TSYM,     &
     &         NBASISMAX,.FALSE.)
!C         WRITE(37,*) I_VIND,NLIST(I_VIND)
         IF(NLIST(I_VIND).GT.(NBASIS-NEL)**2*NEL*NEL) THEN
            WRITE(6,*) "WARNING on excitations"
         ENDIF
!C.. First we need to check that all the nodes in our list are not
!C.. members of the excitation lists of nodes further back in the line
!C.. which have already been processed, or which are yet to be processed
!C.. i.e. if they're already in a list somewhere, then we don't need to
!C.. include them, as they will be processed later.
!C.. As the lists are ordered, this is quite a small job.

!C.. 20050307 - it seems that a better way to do this (which does not
!C.. require that the lists be sorted, and thus rather frees up the 
!C.. excitation generation routine), is to see whether each excitation
!C.. is connected to a previous vertex in the current list.  If it is,
!C.. and we're attempting to connect it to a new vertex, then we should
!C.. remove it from our list.
!C.. This procedure need not be done here, but can be done as we recurse
!C.. through the excitations.
         CALL ELIMDUPS(LSTE,LSTP,I_VIND,NEL,NLIST,ILMAX,NI)



!C.. We now find all the nodes which actually have a connection to us
!C.. eliminating any others in the list.  As this requires that we calc
!C.. RHOIJ, so we store it too.
!C.. II points to the next free position in the list, as we move along,
!C.. rewriting the list to remove the nodes which are zero, and to which
!C.. we are not connected.  
         II=1
         DO I=1,NLIST(I_VIND)
            IF(LSTE(1,I,I_VIND).NE.0) THEN
               CALL CALCRHO2(INODE,LSTE(1:NEL,I,I_VIND),BETA,I_P,NEL,         &
     &            NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,            &
     &            RH,NTAY,ICE(I,I_VIND),ECORE)
               IF(ABS(RH).GT.RHOEPS) THEN
                  IF(II.NE.I)  CALL NECI_ICOPY(NEL,LSTE(1:NEL,I,I_VIND),1,               &
     &                    LSTE(1:NEL,II,I_VIND),1)          
                  ICE(II,I_VIND)=ICE(I,I_VIND)
                  RIJLIST(II,I_VIND)=RH
                  II=II+1
               ENDIF
            ENDIF
         ENDDO
         NLIST(I_VIND)=II-1            
            
!C.. We recurse over all possibilities for the next node:
!C..  1) Nodes connected to us (excluding those connected to previous nodes
!C..        (whether we have processed them or not))
!C..        Once this has been done, our list of connected nodes is marked as
!C..        entirely processed.
!C..  2) Nodes connected to previous nodes which have not yet been processed
!C..        (excluding those which are connected to further previous nodes. 
!C..         NB: These may be connected to us, but will have been excluded
!C..         in 1 as they are connected to previous nodes to 1
!C..  3) Further recursions down the (direct) line from us back to I.  Previous 
!C..        completed lines will have been entirely processed, so we don't
!C..        need to attempt to add nodes to then (but we do need to ensure
!C..        that new nodes don't connect to them)
!C.. 
!C..   All of this boils down to looking at all previous nodes (irrespective
!C..   of which line they're in), and attempting to attach nodes to the
!C..   positions after the parts which have been processed.  This requires
!C..   a COPY of the processed index table upon which to work, as we will 
!C..   be modifying that (i.e. processing nodes after this one or our
!C..   predecessors), GIVEN we have connected up this node.  We will then
!C..   wish to return to past nodes, and connect up further ones (past
!C..   and instead of this one)
        LSTP(I_VIND)=1
        IF(NWHTAY.EQ.0) THEN
!C.. allow chain graphs
         IVLMAX=I_VIND
        ELSE
         IVLMAX=0
        ENDIF
         IVLEVEL=IVLMAX
        DO WHILE(IVLEVEL.GE.0)
!C.. we first take a copy of the pointer table
         DO II=0,I_VIND
            LSTP2(II)=LSTP(II)
         ENDDO
!C.. We must recurse over all the nodes connected to the node at IVLEVEL
         I=LSTP(IVLEVEL)
!C.. remind us which node we're doing
         IF(IVLEVEL.LT.I_VIND) THEN
            CALL NECI_ICOPY(NEL,IPATH(1:NEL,IVLEVEL),1,INODE,1)
         ENDIF
         DO WHILE (I.LE.NLIST(IVLEVEL))
!C.. Tell future recursions that we've processed this node
            LSTP2(IVLEVEL)=I+1
!C.. We also need to update this node in the table we use as the master table
!C.. This is ok as we're not modifying another vertex level's space, and 
!C.. necessary to tell future copies of this master how far we've processed
            IF(IVLEVEL.EQ.I_VIND) LSTP(I_VIND)=I+1
!C.. Because of the elimination of duplicates, we need to ensure that this
!C.. node in the list hasn't been set to 0 (i.e. it was a duplicate, and has
!C.. been removed)
            IF(LSTE(1,I,IVLEVEL).NE.0) THEN
!C.. For VLEVEL=I_VIND, the list at this level has had removed all nodes
!C.. which are connected to previous nodes.
!C.. For VLEVEL<I_VIND, there may be some unprocessed nodes connected to the 
!C.. node at VLEVEL, which are also connected to a node at a higher VLEVEL
!C.. (e.g. this node).  We DO need to process these, because we specifically
!C.  removed them from the list of nodes connected to (e.g.) this level 
!C.. for processing later (i.e. now).
!C.. Calculate RHO_IJ - the connection between the new node and this node
!C.. LSTE(1:NEL,I,IVLEVEL) is the node we're attempting to connect to
               RH=RIJLIST(I,IVLEVEL)
!C.. see if it's worth going further
               RHOIJ(IVLEVEL,I_VIND+1)=RH
               RHOIJ(I_VIND+1,IVLEVEL)=RH
               IF(ABS(RH).LE.RHOEPS) THEN
!C.. we're not connected to this node, so we remove it from our list
                  DO II=1,NEL
                     LSTE(II,I,IVLEVEL)=0
                  ENDDO
               ELSE
!C.. we now calculate the RHOIJ for the previous nodes to connect to the new one
                  DO II=0,I_VIND           
!C.. we've already calculated the RHO to the node to which we are connected.
                        CALL CALCRHO2(IPATH(1:NEL,II), LSTE(1:NEL,I,IVLEVEL),BETA,I_P, &
     &                   NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,                    &
     &                   ALAT,UMAT,RH,NTAY,-1,ECORE)
                        IF(ABS(RH).LE.RHOEPS) THEN
!C.. These two vertices are not connected, which is as we want
                           RH=0.D0
                        ENDIF
!C                        ELSE
!C.. These two vertices are connected

!C***********
!C      STOP 'NB - LOOK AT RED GRAPH PIC'


                        RHOIJ(I_VIND+1,II)=RH
                        RHOIJ(II,I_VIND+1)=RH
                  ENDDO
!C.. 20050307 - At this point we note that if a vertex is connected to
!C.. any vertices before it in the path (excluding the one we excited 
!C.. it from), we need to ignore this vertex, as it will have already been
!C.. counted before.

!C.. Add the new node to the path, and then recurse from it
                  CALL NECI_ICOPY(NEL,LSTE(1:NEL,I,IVLEVEL),1,IPATH(1:NEL,I_VIND+1),1)
!C.. Get the H element so we can calculate the energy
                  
                  hel = get_helement (iPath(:,0), iPath(:,I_VIND+1))
                  HIJS(I_VIND+1) = hel
                  BTABLE(I_VIND+1)=IVLEVEL
                  TOTAL=TOTAL+ FMCPR3NVSTAR(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,         &
     &                     G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,                &
     &                     NTAY,RHOEPS,I_VIND+1,RHOII,RHOIJ,LSTE,                &
     &                     ICE,RIJLIST,NWHTAY,L,LT,I_HMAX,NLIST,                 &
     &                     LSTP2,BTABLE,ILOGGING,TSYM,ECORE,ILMAX,               &
     &                     DBETA,DLWDB,HIJS)
               ENDIF
            ENDIF
            I=I+1
          ENDDO
!C.. Recurse back through the direct line of descendents
          IVLEVEL=BTABLE(IVLEVEL)
         ENDDO
         IF(I_VIND.EQ.0) THEN
!C.. if we're right at the end, we call the star diag
            CLOSE(19)
            CALL READSTARDIAG(L,I_P,RHOII(0),HIJS(0),FMCPR3NVSTAR,DLWDB)
         ELSE
            FMCPR3NVSTARRES=TOTAL
         ENDIF
         RETURN
      END

      RECURSIVE SUBROUTINE ADDDIAGSTAR(IPATH,NEL,RHOIJ,I_V,I_VMAX,ILOGGING,HIJS,N,L)
         IMPLICIT NONE
         INTEGER NEL,I_V,IPATH(NEL,0:I_V),I_VMAX
         INTEGER ILOGGING
         REAL*8 RHOII(0:I_V)
         REAL*8 RHOIJ(0:I_VMAX,0:I_VMAX)
         REAL*8 HIJS(0:I_VMAX)
         
         REAL*8 MAT(I_V,I_V)
         REAL*8 HIJS2(1:I_V)
         
         REAL*8 WLIST(I_V-1)
         INTEGER INFO
         REAL*8 WORK(3*(I_V-1))
         INTEGER ILEN
         INTEGER I,J,K,II,JJ,N,L
         REAL*8 HIA,RIA
         INTEGER IPATH2(NEL,0:I_V-1)
!C         CALL WRITEPATH(19,IPATH,I_V,NEL,.FALSE.)
!C         CALL WRITERHOMAT(19,RHOIJ,I_V,NEL,.FALSE.)
!C         WRITE(19,*) N

         ILEN=I_V-1
         IF(ILEN.EQ.1) THEN
            WRITE(19,"(3G21.12,I4)") RHOIJ(1,1),RHOIJ(0,1),HIJS(1),N
            L=L+1
            RETURN
         ENDIF
            
         MAT=0.d0
!C.. Now we fill the RIJ array
         II=0
         DO I=1,ILEN
            DO J=I,ILEN
               MAT(I,J)=RHOIJ(I,J)
            ENDDO
         ENDDO 

!C.. Diagonalize
         CALL DSYEV('V','U',ILEN,MAT,I_V,WLIST,WORK,3*ILEN,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
!C.. MAT now contains the eigenvectors, and WLIST the eigenvalues
!C.. we now need to transform the Hij and Rij elements into this new
!C.. Corresponding to eigenvalue lambda_a, Rho_ia = Sum_j rho_ij v_ja
         DO I=1,ILEN
            RIA=0.D0
            HIA=0.D0
            DO J=1,ILEN
               RIA=RIA+RHOIJ(0,J)*MAT(J,I)
               HIA=HIA+HIJS(J)*MAT(J,I)
            ENDDO
            WRITE(19,"(3G21.12,I4)") WLIST(I),RIA,HIA,N
         ENDDO
         L=L+ILEN
!C.. Now deal with the subgraphs we have to remove
         DO K=1,ILEN
            II=0
            DO I=0,ILEN
               IF(I.NE.K) THEN
                  II=II+1
                  JJ=0
                  DO J=0,ILEN
                     IF(J.NE.K) THEN
                        JJ=JJ+1
                        MAT(II,JJ)=RHOIJ(I,J)
                     ENDIF
                  ENDDO
                  HIJS2(II)=HIJS(I)
!C                  CALL NECI_ICOPY(NEL,IPATH(1,I),1,IPATH2(1,II-1),1)
               ENDIF
            ENDDO
!C            CALL NECI_ICOPY(NEL,IPATH(1,0),1,IPATH2(1,I_V-1),1)
            CALL ADDDIAGSTAR(IPATH2,NEL,MAT,I_V-1,I_V-1,ILOGGING,HIJS2,-N,L)
         ENDDO
         RETURN
      END
!C.. Use an iterative Order(N) root-finding method to diagonalize the
!C.. star matrix read in from STARDATA
      SUBROUTINE READSTARDIAG(NROOTS1,I_P,RII,HII,SI,DLWDB)
         use global_utilities
         use sort_mod
         IMPLICIT NONE
         INTEGER NROOTS,I_P,NROOTS1
         REAL*8 SI,DLWDB,RII,HII
         REAL*8 ROOTS(0:*),HIJS(0:*),RIJS(0:*),RIIS(0:*)
         POINTER (IP_ROOTS,ROOTS),(IP_HIJS,HIJS),(IP_RIJS,RIJS)
         POINTER (IP_RIIS,RIIS),(IP_NR,NR)
         INTEGER NR(0:*)
         type(timer), save :: proc_timer
         INTEGER I,J,N,II,OII
         REAL*8 NORM,A,B
         NROOTS=NROOTS1
         proc_timer%timer_name='RDSTARDIAG'
         call set_timer(proc_timer)
         CALL N_MEMORY(IP_RIIS,NROOTS+1,'RIIS')
         CALL N_MEMORY(IP_ROOTS,NROOTS+1,'ROOTS')
         CALL N_MEMORY(IP_HIJS,NROOTS+1,'HIJS')
         CALL N_MEMORY(IP_RIJS,NROOTS+1,'RIJS')
         CALL N_MEMORY(IP_NR,NROOTS+1,'NR')
         NR(0)=1
         RIJS(0)=1.D0
!CRII
         RIIS(0)=1.D0
!CRII
         HIJS(0)=HII
         OPEN(19,FILE="STARDATA",STATUS="UNKNOWN")
         DO I=1,NROOTS
            READ(19,*) RIIS(I),RIJS(I),HIJS(I),NR(I)
            RIIS(I)=RIIS(I)/RII
            RIJS(I)=RIJS(I)/RII
         ENDDO
         CLOSE(19)
!C.. we need to sort A and B (and the list of hamil values) into ascending A order
!C         DO I=1,NROOTS
!C            WRITE(6,*) RIIS(I),RIJS(I),NR(I)
!C         ENDDO
           call sort (RIIS(1:nroots), rijs(1:nroots), hijs(1:nroots), &
                      nr(1:nroots))
!C         WRITE(6,*)
!C         DO I=1,NROOTS
!C            WRITE(6,*) RIIS(I),RIJS(I),NR(I)
!C         ENDDO
         II=1
         OII=1
         A=RIJS(1)*RIJS(1)*NR(1)
         DO I=2,NROOTS
            IF(RIIS(I).EQ.RIIS(I-1)) THEN
               II=II+1
               RIIS(II)=RIIS(I)
               RIJS(II)=RIJS(I)
               HIJS(II)=HIJS(I)
               NR(II)=NR(I)
               A=A+RIJS(I)*RIJS(I)*NR(I)
            ELSE
               IF(OII.NE.II) THEN
!C.. if we've just finished a string of nodes with the same diagonal
!C.. element, we see if the total was zero.  If it was, we remove them
!C.. all
                  IF(SQRT(ABS(A)).LT.1.D-9) II=OII
               ENDIF
               OII=II
               IF(II.EQ.0.OR.NR(II).NE.0) II=II+1
               RIIS(II)=RIIS(I)
               RIJS(II)=RIJS(I)
               HIJS(II)=HIJS(I)
               NR(II)=NR(I)
               A=RIJS(I)*RIJS(I)*NR(I)
            ENDIF
         ENDDO
         NROOTS=II
!C         WRITE(6,*)
!C         DO I=1,NROOTS
!C            WRITE(6,*) RIIS(I),RIJS(I),NR(I)
!C         ENDDO
!C.. Find the eigenvalues
         CALL FINDROOTSTAR2(NROOTS,RIIS,RIJS,NR,ROOTS)
!C         WRITE(6,*)
!C         DO I=1,NROOTS
!C            WRITE(6,*) RIIS(I)*RII,ROOTS(I)*RII
!C         ENDDO
         SI=0.D0
         DO I=NROOTS,0,-1
!C            WRITE(6,*) I,"ROOT:",RII*ROOTS(I),RII*RIIS(I)
            IF(ROOTS(I).EQ.RIIS(I).OR.ROOTS(I).EQ.0.D0) THEN
!C.. If we're in a degenerate set of eigenvectors, we calculate things a
!C.. little differently
!C.. k is the vertex which the degeneracies couple to
!C.. and j is the current degenerate element
!C.. <Di|H|Dj><Dj|Psi> L**P <Psi|Di>
!C.. <Dk|Psi>=rho_ij/NORM
!C.. <Dj|Psi>=-rho_ik/NORM
!C.. NORM=rho_ij**2+rho_ik**2
!C.. <Di|Psi>=0 so we have no contributions at all!


            ELSE
!C.. We need to calculate the normalization of each eigenvector
               NORM=1.D0
               DO J=1,NROOTS
                  NORM=NORM+NR(J)*(RIJS(J)/(ROOTS(I)-RIIS(J)))**2
               ENDDO
!C               WRITE(6,*) "NORM**2:",NORM
!C               WRITE(6,*) 1,1.D0/SQRT(NORM)
!C               DO J=1,NROOTS
!C                  WRITE(6,*) NR(J),
!C     &               RII*(RIJS(J)/(ROOTS(I)-RIIS(J)))/SQRT(NORM)
!C               ENDDO 
!C.. We add in the first element of the eigenvector * lambda**P
               SI=SI+(ROOTS(I)**I_P)*1.D0/NORM
!C               WRITE(6,*) "SI",SI
!C.. We do not currently calculate the DLWDB

!C               IF(DBETA.NE.0.D0) THEN
!C                  DLWDB2=HII
!C                  WRITE(6,*) LIST(1,2),SQRT(1/NORM)
!C                  DO J=1,NROOTS
!C.. NOTE - this is not correct - we might need imaginaries here.
!C                     DLWDB2=DLWDB2+
!C     &                  HIJS(J)*(LIST(J,1)/(ROOTS(I+1)-LIST(J,0)))
!C                WRITE(6,*) LIST(J,2),
!C     &            LIST(J,1)/((ROOTS(I+1)-LIST(J,0))*SQRT(NORM))
!C                  ENDDO
!C                  DLWDB=DLWDB+DLWDB2*(ROOTS(I+1)**I_P)/NORM
!C               ENDIF
!C               WRITE(6,*) ROOTS(I+1)**I_P,DLWDB2*(ROOTS(I+1)**I_P)/NORM
!C               WRITE(6,*)
            ENDIF
         ENDDO
         SI=SI-1.D0
!C         DLWDB=DLWDB-LIST(1,2)
         CALL N_FREEM(IP_RIIS)
         CALL N_FREEM(IP_ROOTS)
         CALL N_FREEM(IP_HIJS)
         CALL N_FREEM(IP_RIJS)
         CALL N_FREEM(IP_NR)
         call halt_timer(proc_timer)
         RETURN
      END

!.. Read in values for a weight and ETilde from an MCPATHS file
      SUBROUTINE ReadMCPaths(iV,wWeight,wETilde)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iV,i1,i2,i3
         CHARACTER*1 c
         REAL*8 r1,r2,r3,r4
         real(dp) wWeight,wETilde
         
         OPEN(43,FILE="MCPATHS",STATUS="OLD")
         READ(43,*) c
         READ(43,*) c
         i1=0
         DO WHILE(i1.NE.iV)
            READ(43,*) i1,r1,r2,r3,i2,i3,r4
         ENDDO
         wWeight=r1
         wETilde=r4
         CLOSE(43)
      END
      
!.. Function to get Arr out of NMAX - passed in two integers      
      FUNCTION ORBENERGY(ARR,orbnum)
          use constants, only: dp
          IMPLICIT NONE
          REAL*8 ARR(*)
          INTEGER orbnum
          REAL*8 ORBENERGY
          ORBENERGY=ARR(orbnum)
          RETURN
      END FUNCTION ORBENERGY
