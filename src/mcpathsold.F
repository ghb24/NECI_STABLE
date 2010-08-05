C.. Calculate RHO^(P)_II without having a stored H matrix
C.. Sum over gradually increasing H.  Paths generated are e.g. IJIKLI
C.. which do not necessarily have distinct nodes
      SUBROUTINE  MCPATHSR2(NI,BETA,I_P,I_HMAX,NEL,NBASISMAX,G1,
     &               NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &               LSTE,ICE,NWHTAY,ILOGGING,ECORE,ILMAX,WLRI,WLSI,
     &               TSYM)
         use global_utilities
         use util_mod, only: NECI_ICOPY
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION IPATH(NEL,0:I_HMAX)
         INTEGER NI(NEL)
         REAL*8 F(2:I_HMAX)
         CHARACTER*20 STR
         REAL*8 TOTAL,RHOII(0:I_HMAX)
         REAL*8 FLNSUM,FLNRII
         REAL*8 ALAT(3),RHOEPS
         INTEGER  G1(5,NBASIS),nBasisMax(5,*),ILMAX
         INTEGER LSTE(NEL,ILMAX,0:I_HMAX-1)
         INTEGER ICE(ILMAX,0:I_HMAX-1),NLIST
         LOGICAL TLOG
         type(timer), save :: proc_timer,proc_timer2
         proc_timer%timer_name='MCPATHSR2 '
         call set_timer(proc_timer)
         TLOG=BTEST(ILOGGING,1)
         IF(TLOG)
     &      OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,2))
     &      OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
C.. Allocate some mem
C.. We'll start with I_H=2
C.. Set the first node to I_I
c         OPEN(10,FILE="HPATHS",STATUS="UNKNOWN")
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(0)=RH
         WLRI=LOG(RHOII(0))
         TOTAL=1.D0
         IF(TLOG) WRITE(11,"(I12,2G25.16,F19.7,2I12)")
     &          0,TOTAL,TOTAL,0.D0,0,0
c         WRITE(6,*) 0,TOTAL,TOTAL,0
         DO I_H=2,I_HMAX
CI_P
c            WRITE(6,"I2") I_H
            WRITE(STR,"(A,I5)") "FMCPR",I_H
            proc_timer2%timer_name=STR
            call set_timer(proc_timer2)
            L=0
            LT=0
            F(I_H)=FMCPR2(NI,BETA,I_P,IPATH,I_H,NEL,NBASISMAX,G1,NBASIS,
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &            0,1.D0,RHOII,LSTE,ICE,NWHTAY,L,LT,ILOGGING,TSYM,
     &            ECORE,ILMAX)
            call halt_timer(proc_timer2)
            TOTAL=TOTAL+F(I_H)
c            WRITE(6,*) I_H,F(I_H),TOTAL,get_total_time(proc_timer2),L,LT
            IF(TLOG)
     &            WRITE(11,"(I12,2G25.16,F19.7,2I12)")
     &               I_H,F(I_H),TOTAL,get_total_time(proc_timer2),L,LT
         ENDDO
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,2)) CLOSE(10)
         WLSI=LOG(TOTAL)
C         CALL N_FREEM(IP_LSTE)
C         CALL N_FREEM(IP_ICE)
         call halt_timer(proc_timer)
         RETURN
      END
C.. A function to loop recursively over each node set choosing a different
C.. node for each set. (7/11/03)
C..  It then weights using CALCPATHn (n=5 for rec/taylor2 hybrid) (18/11/03)
C.. This does not use a stored RHO or H matrix
      RECURSIVE REAL*8 FUNCTION FMCPR2(NI,BETA,I_P,IPATH,I_H,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,I_HIND,RHOCUML,RHOII,LSTE,ICE,NWHTAY,L,LT,ILOGGING,
     &   TSYM,ECORE,ILMAX) RESULT (FMCPR2RES)
         
         use util_mod, only: NECI_ICOPY
         IMPLICIT REAL*8 (A-H,O-Z)
         REAL*8 RHOCUML
         INTEGER IPATH(NEL,0:I_H)
         REAL*8 RHOII(0:I_H)
         INTEGER INODE(NEL),ILMAX
C.. LSTE is a list of excitations (which we will generate)
C.. ICE is the IC of each excitation (i.e. how much it differs from us (INODE)

         INTEGER LSTE(NEL,ILMAX,0:I_H-1),NI(NEL)
         INTEGER ICE(ILMAX,0:I_H-1),NLIST
         REAL*8 RH,BETA,RHOEPS
         LOGICAL LSAME,TSYM
         LOGICAL TLOG1,TLOG2
         TLOG1=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TOTAL=0.D0
C.. This is the current node (set by our parent)         
         CALL NECI_ICOPY(NEL,IPATH(1:NEL,I_HIND),1,INODE,1)
C.. Find all nodes connected to this node
C.. This is of maximum order NTAY*2
         IF(I_HIND.EQ.I_H-1) LT=LT+1
         IF(I_HIND.GT.INT(I_H/2)) THEN
C.. if we're at the last but one node, the only allowable excitation is I
C.. although it might not be allowed.
            I=IGETEXCITLEVEL(NI,INODE,NEL)
C.. if the current node is more than the allowable excitations from I, then
C.. give up and return
            IF(I.GT.2*ABS(NTAY)*(I_H-I_HIND)) THEN
               FMCPR2RES=0.D0
               RETURN
            ENDIF
         ENDIF
         IF(I_HIND.EQ.I_H-1) THEN
            L=L+1
C.. otherwise just set the only possible allowed excitation to I
            CALL NECI_ICOPY(NEL,NI,1,LSTE(1:NEL,1,I_HIND),1)
            I=IGETEXCITLEVEL(NI,INODE,NEL)
            ICE(1,I_HIND)=I
            NLIST=1
         ELSE
C.. 1 at end means we need at least 1 excitation
            CALL GENEXCIT(INODE,ABS(NTAY*2),NBASIS,NEL,
     &         LSTE(1,1,I_HIND),ICE(1,I_HIND),NLIST,1,G1,TSYM,
     &         NBASISMAX,.FALSE.)

         ENDIF
C.. Find RHO_II for this node
         CALL CALCRHO2(INODE,INODE,BETA,I_P,NEL,NBASISMAX,
     &   G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(I_HIND)=RH
         I=1
C.. We must recurse over all the nodes connected to this node
         DO WHILE (I.LE.NLIST)
          CALL CALCRHO2(INODE,LSTE(1:NEL,I,I_HIND),BETA,I_P,NEL,
     &      NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,
     &            ICE(I,I_HIND),ECORE)
C.. see if it's worth going further
          IF(ABS(RH).GT.RHOEPS) THEN
            IF(I_HIND.EQ.(I_H-1)) THEN
C.. If we're at the last but one node, we can only connect to I_I
C.. (this has been accounted for at the list stage
C.. Bingo        
               CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,I_H),1) 
               RHOII(I_H)=RHOII(0)
               IF(TLOG1) THEN
                  WRITE(10,"(A)",advance='no') "["
                  DO J=0,I_H
                     WRITE(10,"(A)",advance='no')
                     DO K=1,NEL
                        WRITE(10,"(I3,A)",advance='no') IPATH(K,J),","
                     ENDDO
                     WRITE(10,"(A)",advance='no') "),"
                  ENDDO
                  WRITE(10,*) "]"
                  IF(TLOG2) THEN
                     WRITE(10,"(A)",advance='no') "("
                     DO J=1,I_H-1
                        WRITE(10,"(E25.16,A)",advance='no')
     &                                            RHOII(J)/RHOII(0),","
                     ENDDO
                     WRITE(10,*) ")"
                  ENDIF
               ENDIF
               TOTAL=TOTAL+
     &               CALCPATH7(IPATH,RHOII,RHOCUML*RH,I_H,
     &                  I_P,NWHTAY)
               IF(TLOG1) WRITE(10,*) TOTAL
            ELSE
C.. Recurse to the next node
C.. We're not allowed to connect to ourself again (accounted for in the
C..  excitations)
C.. If we're the antepenultimate node, we're not allowed to connect
C.. to I_I as the penultimate node, because that HAS to connect to 
C,, I_I as the ultimate node
               IF(I_HIND.NE.(I_H-2)) THEN
                  LSAME=.FALSE.
               ELSE
                  LSAME=.TRUE.
                  DO J=1,NEL
                     IF(LSTE(J,I,I_HIND).NE.NI(J)) LSAME=.FALSE.
                  ENDDO
               ENDIF
               IF(.NOT.LSAME) THEN
                  CALL NECI_ICOPY(NEL,LSTE(1:NEL,I,I_HIND),1,
     &                  IPATH(1:NEL,I_HIND+1),1)
                  TOTAL=TOTAL+
     &                  FMCPR2(NI,BETA,I_P,IPATH,I_H,NEL,NBASISMAX,
     &                   G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,
     &                   NTAY,RHOEPS,I_HIND+1,RHOCUML*RH,RHOII,LSTE,ICE,
     &                   NWHTAY,L,LT,ILOGGING,TSYM,ECORE,ILMAX)
               ENDIF
            ENDIF
          ENDIF
          I=I+1
         ENDDO
         FMCPR2RES=TOTAL
         RETURN
      END

C.. Go systematically through paths with I_H=0,2,3,... non-self-loops
C.. For each path, generate the weight, and weight appropriately. (7/11/03)
C.. The weight can only be calculated up to I_H=4 (18/11/03) using Taylor
C.. expansion.
C.. This still requires the whole Rho matrix.
C.. We return ln (RHOII^(P)), decomposing it into
C.. P ln(RHOII) +ln(1+RHOII^-2 Sum_J (RHOIJ^2 w2(RHOJJ,RHOII,P)) + ...)
C.. returning 
C..  FLNRII = ln(RHOII) and
C..  FLNSUM =ln(1+RHOII^-2 Sum_J (RHOIJ^2 w2(RHOJJ,RHOII,P)) + ...)

      SUBROUTINE MCPATHSR(I_I,I_P,NDET,RHOIJROWS,RHOIJ,RHOIJLAB,
     &               NRHOMAXWIDTH,I_HMAX,
     &                   NWHTAY,ILOGGING,FLNRII,FLNSUM)
         use global_utilities
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER NDET,NRHO
         INTEGER RHOIJLAB(*),RHOIJROWS(NDET+1)
         REAL*8 RHOIJ(*)
         DIMENSION IPATH(0:I_HMAX)
         DIMENSION RHOII(0:I_HMAX)
         INTEGER ISEED
         REAL*8 F(2:I_HMAX)
         CHARACTER*10 STR
         REAL*8 TOTAL
         REAL*8 FLNSUM,FLNRII
         LOGICAL TLOG
         type(timer), save :: proc_timer,proc_timer2
         proc_timer%timer_name='MCPATHSR  '
         call set_timer(proc_timer)
         TLOG=BTEST(ILOGGING,1)
         IF(TLOG) OPEN(11,FILE="MCPATHS",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,2)) OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
C.. We'll start with I_H=2
C.. Set the first node to I_I
         IPATH(0)=I_I
         RHOII(0)=GETRHOELEMENT(I_I,I_I,RHOIJ,RHOIJLAB,RHOIJROWS,NDET)
         FLNRII=LOG(RHOII(0))
         TOTAL=1.D0
C         WRITE(10,*) "I=",I_I
         IF(TLOG) WRITE(11,"(I12,2G25.16,F19.7,2I12)")
     &       0,TOTAL,TOTAL,0.D0,0,0
c         WRITE(6,*) 0,TOTAL,TOTAL,0
         DO I_H=2,I_HMAX
CI_P
c            WRITE(6,"I2") I_H
            WRITE(STR,"(A,I5)") "FMCPR",I_H
            proc_timer2%timer_name=STR
            call set_timer(proc_timer2)
            L=0
            LT=0
            F(I_H)=FMCPR(I_I,I_P,NDET,RHOIJROWS,RHOIJ,RHOIJLAB,
     &          NRHOMAXWIDTH,IPATH,ISEED,I_H,NWHTAY,0,1.D0,RHOII,L,LT,
     &          ILOGGING)
            call halt_timer(proc_timer2)
            TOTAL=TOTAL+F(I_H)
c            WRITE(6,*) I_H,F(I_H),TOTAL,get_total_time(proc_timer2),L,LT
c     &                  ,(L+0.D0)/LT
            IF(TLOG)
     &         WRITE(11,"(I12,2G25.16,F19.7,2I12)")
     &         I_H,F(I_H),TOTAL,get_total_time(proc_timer2),L,LT
         ENDDO
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,3)) CLOSE(10)
         FLNSUM=LOG(TOTAL)
         call halt_timer(proc_timer)
         RETURN
      END

C.. A function to loop recursively over each node set choosing a different
C.. node for each set. (7/11/03)
C..  It then weights using CALCPATHn (n=5 for rec/taylor2 hybrid) (18/11/03)
      RECURSIVE REAL*8 FUNCTION FMCPR(I_I,I_P,NDET,RHOIJROWS,RHOIJ,
     &      RHOIJLAB,NRHOMAXWIDTH,IPATH,ISEED,I_H,NWHTAY,
     &      I_HIND,RHOCUML,RHOII,L,LT,ILOGGING) RESULT (FMCPRRES)
         
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER NDET,NRHO,I_H
         INTEGER RHOIJLAB(*),RHOIJROWS(NDET+1)
         REAL*8 RHOIJ(*)
         REAL*8 RHOCUML
         INTEGER IPATH(0:I_H)
         REAL*8 RHOII(0:I_H)
         INTEGER ISEED
         INTEGER INTERSECT(NRHOMAXWIDTH)
         LOGICAL TLOG1, TLOG2
         TLOG1=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)

         TOTAL=0.D0
C.. This is the current node (set by our parent)
         INODE=IPATH(I_HIND)
C.. Find all nodes connected to this node
         IMIN=RHOIJROWS(INODE)
         IMAX=RHOIJROWS(INODE+1)        
C.. Find RHO_II for this node
         I=IMIN
         RHOII(I_HIND)=0.D0
         DO WHILE(I<IMAX)
            IF(RHOIJLAB(I).EQ.INODE) THEN
               RHOII(I_HIND)=RHOIJ(I)
               I=IMAX
            ENDIF
            I=I+1
         ENDDO
         I=IMIN
C.. We must recurse over all the nodes connected to this node
         DO WHILE (I<IMAX)
            IF(I_HIND.EQ.(I_H-1)) THEN
C.. If we're at the last but one node, we can only connect to I_I
               IF(RHOIJLAB(I).EQ.I_I) THEN
C.. Bingo         
                  IPATH(I_H)=I_I
                  RHOII(I_H)=RHOII(0)
 
                 IF(TLOG1) THEN
                     WRITE(10,"(A)",advance='no') "["
                     DO J=0,I_H
                        WRITE(10,"(I3,A)",advance='no') IPATH(J),","
                     ENDDO
                     WRITE(10,*) "]"
                     IF(TLOG2) THEN
                        WRITE(10,"(A)",advance='no') "("
                        DO J=1,I_H-1
                          WRITE(10,"(E25.16,A)",advance='no')
     &                                            RHOII(J)/RHOII(0),","
                        ENDDO
                        WRITE(10,*) ")"
                     ENDIF
                  ENDIF
                  L=L+1
                  TOTAL=TOTAL+
     &               CALCPATH7(IPATH,RHOII,RHOCUML*RHOIJ(I),I_H,
     &                  I_P,NWHTAY)
                  IF(TLOG1) WRITE(10,*) TOTAL
                  I=IMAX
               ENDIF
            ELSE
C.. Recurse to the next node
               LT=LT+1
               IF(RHOIJLAB(I).NE.INODE) THEN
C.. We're not allowed to connect to ourself again
C.. If we're the antepenultimate node, we're not allowed to connect
C.. to I_I as the penultimate node, because that HAS to connect to 
C,, I_I as the ultimate node
                  IF(I_HIND.NE.(I_H-2).OR.RHOIJLAB(I).NE.I_I) THEN
                     IPATH(I_HIND+1)=RHOIJLAB(I)
                     TOTAL=TOTAL+
     &                     FMCPR(I_I,I_P,NDET,RHOIJROWS,RHOIJ,RHOIJLAB,
     &                     NRHOMAXWIDTH,IPATH,ISEED,I_H,NWHTAY,I_HIND+1,
     &                     RHOCUML*RHOIJ(I),RHOII,L,LT,ILOGGING) 
                  ENDIF
               ENDIF
            ENDIF
            I=I+1
         ENDDO
         FMCPRRES=TOTAL
         RETURN
      END


C.. Randomly walk through paths length I_P in CI space. v inefficient (5/11/03)
C.. (Requires whole rho matrix)
      SUBROUTINE MCPATHS0(I_I,I_P,NDET,RHOIJROWS,RHOIJ,RHOIJLAB,
     &               NRHOMAXWIDTH)
         use global_utilities
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER NDET,NRHO
         INTEGER RHOIJLAB(*),RHOIJROWS(NDET+1)
         REAL*8 RHOIJ(*)
         INTEGER IPATH(I_P+1)
         INTEGER ISEED
         INTEGER INTERSECT(NRHOMAXWIDTH)
         REAL*8 PROB(NRHOMAXWIDTH)
         REAL*8 RAN2
         REAL(KIND=4) :: SecondSeed=0.0
         type(timer), save :: proc_timer
         proc_timer%timer_name='MCPATHS0  '
         call set_timer(proc_timer)
C.. Randomize the seed:
         call system_clock(iseed)
         IF(MOD(ISEED,2).EQ.0) ISEED=ISEED-1
C.. THe seed must be odd
         
         INIT=5
         RTOT=0.D0               
C.. Start with all nodes at I
         DO N=1,I_P+1
            IPATH(N)=I_I
         ENDDO    

C.. Open our output data file
         OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
C.. Do 1000 iterations
         WRITE(10,"(I12,E25.16)") I_P,GETRHOELEMENT(I_I,I_I,RHOIJ,
     &                RHOIJLAB,RHOIJROWS,NDET)**I_P
         DO I=1,100
C.. pick a random element N in IPATH excluding start and end
            N=INT(RAN2(ISEED)*(I_P-1)+2)
            N=INT(RAN2(ISEED)*(I_P-1)+2)
C.. The next free position in the list of intersections
            INTPOS=1
C.. Now get the intersections of IPATH(N-1) and IPATH(N+1)
            IROWJ=RHOIJROWS(IPATH(N-1))
            IROWK=RHOIJROWS(IPATH(N+1))
            ILJ=RHOIJROWS(IPATH(N-1)+1)
            ILK=RHOIJROWS(IPATH(N+1)+1)
            K=IROWK
            J=IROWJ
            DO WHILE ((J.LE.ILJ).AND.(K.LE.ILK))
               ILABJ=RHOIJLAB(J)
               DO WHILE ((RHOIJLAB(K).LT.ILABJ).AND.(K.LE.ILK))
                  K=K+1
               ENDDO
               IF((ILABJ.EQ.RHOIJLAB(K)).AND.(K.LE.ILK)) THEN
                  INTERSECT(INTPOS)=ILABJ
                  INTPOS=INTPOS+1
               ENDIF
               J=J+1
            ENDDO
C.. In order to Monte carlo, we need to calculate the sum of all possible
C.. RHO_(N-1),N RHO_N,(N+1) to weight our calculations
C.. save old W
            OW=W
            SUM=0.D0
            DO J=1,INTPOS
               P=GETRHOELEMENT(IPATH(J-1),INTERSECT(J),RHOIJ,
     &                        RHOIJLAB,RHOIJROWS,NDET)
     &             *GETRHOELEMENT(INTERSECT(J),IPATH(J+1),RHOIJ,
     &                        RHOIJLAB,RHOIJROWS,NDET)
               SUM=SUM+ABS(P)
               PROB(J)=P
            ENDDO
            W=SUM
C.. We now have INTPOS intersections.  we generate with appropriate probablity
C.. If we're initializing, we just choose randomly, otherwise we weight
            IF(INIT.EQ.0) THEN
               POS=RAN2(ISEED)*W
               SUM=0.D0
               IPOS=0
               DO WHILE (SUM.LT.POS)
                  IPOS=IPOS+1
                  SUM=SUM+ABS(PROB(IPOS))
               ENDDO
            ELSE
               IPOS=INT(RAN2(ISEED)*(INTPOS-1))+1
               INIT=INIT-1
            ENDIF
C.. Now check acceptance
C            IF(INIT.EQ.0) THEN
C               ACC=W
C.. Now update the path
            IPATH(N)=INTERSECT(IPOS)
C.. Now calc rho
            RH=1.D0
            DO J=2,I_P+1
               RH=RH*GETRHOELEMENT(IPATH(J-1),IPATH(J),RHOIJ,
     &                        RHOIJLAB,RHOIJROWS,NDET)
            ENDDO
               
C.. Now count the number of self-loops, and output the data
            ISELF=0
            DO J=2,I_P+1
               IF(IPATH(J).EQ.IPATH(J-1)) ISELF=ISELF+1
            ENDDO
            RTOT=RTOT+RH
           WRITE(6,"(I3,E25.16,2I3,A)",advance='no')
     &                                   ISELF,RH,INIT,INTPOS," int.  "
C.. Let's print out the whole path
            DO J=1,I_P+1
               WRITE(6,"(I4,A)",advance='no') IPATH(J),"," 
            ENDDO
            WRITE(6,*) W
            WRITE(10,"(I12,2E25.16)") ISELF,RH,RTOT
            
         ENDDO       
         CLOSE(10)
         call halt_timer(proc_timer)
         RETURN
      END

      SUBROUTINE LOGNAN(NI,NEL,BETA,ITIME)
         INTEGER NEL,ITIME,NI(NEL)
         REAL*8 BETA
         CHARACTER*40 STR
         WRITE(STR,*) 'mkdir "nan',ITIME,'"'
         I=SYSTEM(STR)
         WRITE(STR,*) 'cp * "nan',ITIME,'"'
         I=SYSTEM(STR)
         WRITE(STR,*) 'nan',ITIME,"/LOG"
         OPEN(12,FILE=STR,STATUS="UNKNOWN")
         WRITE(12,*) "WARNING: nan found at time",ITIME
         WRITE(12,*) "Beta=",BETA
         WRITE(12,"(A)",advance='no') "  nan det=("
         DO K=1,NEL
            WRITE(12,"(I3,A)",advance='no') NI(K),","
         ENDDO
         WRITE(12,"(A)",advance='no') "),"
         CLOSE(12)
         RETURN

      END

C      SUBROUTINE WRITE_PATH(IUNIT,IPATH,NEL,I_V)
C         IMPLICIT NONE
C         INTEGER NEL,I_V,IPATH(NEL,0:I_V)
C         INTEGER K,J,IUNIT
C         WRITE(IUNIT,"(A)",advance='no') "["
C         DO J=0,I_V
C            WRITE(IUNIT,"(A)",advance='no') "("
C            DO K=1,NEL
C               WRITE(IUNIT,"(I3,A)",advance='no') IPATH(K,J),","
C            ENDDO
C            WRITE(IUNIT,"(A)",advance='no') "),"
C         ENDDO
C         WRITE(IUNIT,"(A)",advance='no') "] "
C         RETURN
C      END
      
C      SUBROUTINE WRITE_PATH2(IUNIT,IPATH,NEL,I_V,LTERM)
C         IMPLICIT NONE
C         INTEGER NEL,I_V,IPATH(NEL,0:I_V)
C         INTEGER K,J,IUNIT
C         LOGICAL LTERM
C         WRITE(IUNIT,"(A)",advance='no') "["
C         DO J=0,I_V
C            WRITE(IUNIT,"(A)",advance='no') "("
C            DO K=1,NEL
C               WRITE(IUNIT,"(I3,A)",advance='no') IPATH(K,J),","
C            ENDDO
C            WRITE(IUNIT,"(A)",advance='no') "),"
C         ENDDO
C         WRITE(IUNIT,"(A)",advance='no') "] "
C         IF(LTERM) WRITE(IUNIT,*)
C         RETURN
C      END
C      SUBROUTINE WRITE_RHOMATRIX(IUNIT,RHOIJ,NEL,I_V)
C         IMPLICIT NONE
C         INTEGER NEL,I_V,J,K,IUNIT
C         REAL*8 RHOIJ(0:I_V,0:I_V)      
C         WRITE(IUNIT,"(A)",advance='no') "("
C         DO J=0,I_V-1
C            DO K=0,I_V-1
C               WRITE(IUNIT,"(E25.16,A)",advance='no') RHOIJ(J,K),","
C            ENDDO
C            WRITE(IUNIT,"(A)",advance='no') "|"
C         ENDDO
C         WRITE(IUNIT,*) ")"
C         RETURN
C      END

      SUBROUTINE WRITE_XMATRIX(IUNIT,RHOIJ,I_V)
         IMPLICIT NONE
         INTEGER I_V,J,K,IUNIT
         REAL*8 RHOIJ(0:I_V-1,0:I_V-1)      
         WRITE(IUNIT,"(A)",advance='no') "{"
         DO J=0,I_V-1
            DO K=0,I_V-1
               WRITE(IUNIT,"(E25.16,A)",advance='no') RHOIJ(J,K),","
            ENDDO
            WRITE(IUNIT,"(A)",advance='no') "|"
         ENDDO
         WRITE(IUNIT,*) "}"
         RETURN
      END

