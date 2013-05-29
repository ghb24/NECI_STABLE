! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
!Molpro requires compilation with the NAG compiler.
!This is not going to be possible without extensive rewriting to
!get around ajwt's horrific hiding of variables inside other types.
!Therefore, I am slashing all this file. GHB(24/8/11)
#ifndef MOLPRO

C.. Calculate RHO^(P)_II without having a stored H matrix
C.. SAMPLE over distinct nodes, e.g. IJKLI, with paths up to I_HMAX
C.. generated from these, and summed (e.g IJILKJI), up to max H
C.. In theory more efficient because RHO_IJ,RHO_JK, etc are calculated
C.. just once for all these paths.
C.. I_VMAX is the max number of distinct vertices in a path.
C.. I_HMAX is the max number of hops in a path.
C.. NWHTAY contains the number of samples to take for each level
C.. This is an approach based on Markov Chains
      SUBROUTINE  MCPATHSR5(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,
     &               NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &               NWHTAY,ILOGGING,ECORE, WLRI,WLSI)
         use global_utilities
         use util_mod, only: isnan_neci, neci_etime
         use util_mod, only: NECI_ICOPY
         use constants, only: dp,sp
!         IMPLICIT real(dp) (A-H,O-Z)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_VMAX
         integer :: IPATH(NEL,0:I_VMAX)
         real(dp) F(2:I_VMAX),FSQ(2:I_VMAX),UMAT(*)
         CHARACTER(20) STR
         real(dp) TOTAL
         real(dp) RHOII(0:I_VMAX),RHOIJ(0:I_VMAX,0:I_VMAX)
         real(dp) ALAT(3),RHOEPS,RH,BETA
         INTEGER  G1(*),nBasisMax(5,*),NMAX
CNEL,0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
C0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER BTABLE(0:I_VMAX),ILOGGING,K,NBASIS,I_P
         LOGICAL TLOG
         INTEGER ICOUNT,ISEED,I_VMOVE,I_HMAX,NWHTAY,LT,L,I_V
         INTEGER IADJ(0:I_VMAX-1,0:I_VMAX-1),ADJ,NTAY,NMSH
         real(dp) OVALUE,WEIGHT,NVALUE,WLSI,fck,BRR(*)
         LOGICAL LCANMOVE
         EXTERNAL FMCPR4WFUNC,FMCPR4VFUNC,FMCPR4UFUNC
         real(dp) RNORM,ND,STD,GETNLOOPS,FMCPR5,WLRI,ECORE
         real(sp) itime,tarr(2)
         type(timer), save :: proc_timer
         type(timer), save :: proc_timer2
         proc_timer%timer_name='MCPATHSR5 '
         call set_timer(proc_timer)
         TLOG=BTEST(ILOGGING,1)
         IF(TLOG) THEN
            OPEN(11,FILE="MCPATHS",STATUS="OLD",POSITION='APPEND')
C.. go to end of file
!            I=FSEEK(11,0,2)
            WRITE(11,"(A)",advance='no') "("
            DO K=1,NEL
               WRITE(11,"(I3,A)",advance='no') NI(K),","
            ENDDO
            WRITE(11,"(A)") ")"
         ENDIF
         IF(BTEST(ILOGGING,2))
     &      OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
C.. Set the first node to I_I
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(0)=RH
         RHOIJ(0,0)=RH
         IADJ(0,0)=1
         WLRI=LOG(RHOII(0))
         TOTAL=1.0_dp
         IF(TLOG) WRITE(11,"(I12,2G25.16,F19.7,2I12)") 
     &         1,TOTAL,TOTAL,0.0_dp,1,1
c         WRITE(6,*) 0,TOTAL,TOTAL,0
C.. I_V is the number of vertices in the path
         DO I_V=2,I_VMAX
            WRITE(STR,"(A,I5)") "FMCPR",I_V
            proc_timer2%timer_name=STR
            call set_timer(proc_timer2)
            L=0
            LT=0
            BTABLE(0)=0
            F(I_V)=0.0_dp
            FSQ(I_V)=0.0_dp
            ISEED=7
C.. I_VMOVE=0 indicates we need to generate a valid path.
            CALL GENVALIDPATH(BETA,I_P,IPATH,I_V,NEL,
     &         NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &         ECORE,ISEED,I_VMOVE,NVALUE,IADJ,ADJ,LCANMOVE,
     &         FMCPR4WFUNC,FMCPR4VFUNC,WEIGHT)
            IF(WEIGHT.GT.0.0_dp) THEN
            DO ICOUNT=1,NWHTAY
               OVALUE=NVALUE
               NVALUE=
     &            FMCPR5(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &              RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &             ECORE,ISEED,I_VMOVE,OVALUE,WEIGHT,IADJ,ADJ,LCANMOVE,
     &               FMCPR4WFUNC,FMCPR4VFUNC)
               F(I_V)=F(I_V)+NVALUE
               FSQ(I_V)=FSQ(I_V)+NVALUE*NVALUE
            ENDDO
            F(I_V)=F(I_V)/NWHTAY
            FSQ(I_V)=FSQ(I_V)/NWHTAY
            ISEED=7
            CALL GENVALIDPATH(BETA,I_P,IPATH,I_V,NEL,
     &         NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &         ECORE,ISEED,I_VMOVE,NVALUE,IADJ,ADJ,LCANMOVE,
     &         FMCPR4UFUNC,FMCPR4WFUNC,WEIGHT)
            RNORM=0.0_dp
            DO ICOUNT=1,NWHTAY
               OVALUE=NVALUE
               NVALUE=
     &            FMCPR5(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &              RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &             ECORE,ISEED,I_VMOVE,OVALUE,WEIGHT,IADJ,ADJ,LCANMOVE,
     &               FMCPR4UFUNC,FMCPR4WFUNC)
               RNORM=RNORM+NVALUE
            ENDDO
            RNORM=RNORM/NWHTAY
            ND=GETNLOOPS(I_V,NBASIS,NEL)
C            WRITE(6,*) F(I_V),F(I_V)*RNORM*ND,RNORM,ND,RNORM*ND
            RNORM=RNORM*ND
            F(I_V)=F(I_V)*RNORM
            FSQ(I_V)=FSQ(I_V)*RNORM*RNORM
            ELSE
               F(I_V)=0.0_dp
            ENDIF
            call halt_timer(proc_timer2)
            TOTAL=TOTAL+F(I_V)
c            WRITE(6,*) I_V,F(I_V),TOTAL,get_total_time(proc_timer2),L,LT
            STD=SQRT(ABS(FSQ(I_V)-F(I_V)*F(I_V)))
            IF(TLOG)
     &         WRITE(11,"(I12,2G25.16,F19.7,I12,F19.7)")
     &            I_V,F(I_V),TOTAL,get_total_time(proc_timer2),L,STD
            IF(ISNAN_neci(F(I_V))) THEN
C.. save all log files
               ITIME=neci_etime(tarr)
               CALL neci_flush(11)
               CALL LOGNAN(NI,NEL,BETA,ITIME)
               WRITE(6,*) "WARNING: nan found at time",ITIME
               WRITE(6,"(A)",advance='no') "  nan det=("
               DO K=1,NEL
                  WRITE(6,"(I3,A)",advance='no') NI(K),","
               ENDDO
               WRITE(6,"(A)",advance='no') "),"
            ENDIF
         ENDDO
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,2)) CLOSE(10)
         WLSI=LOG(TOTAL)
         call halt_timer(proc_timer)
         RETURN
      END



C.. A function which chooses a random set of I_V connected dets, working out
C.. loop contribution for that set.
C.. All nodes are distinct.  Paths IJIKJI etc.
C.. are generated by permutation from IJKI, and summed up to length I_HMAX
C.. using the appropriate weightings (Z-sums) from CALCPATHS.(03/07/04). 
C.. This function assumes that there are enough available excitations to 
C.. form a loop of length I_V.  If not it will probably hang.
      FUNCTION FMCPR5(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &   ECORE,ISEED,I_VMOVE,VALUE,WEIGHT,IADJ,ADJ,LCANMOVE,
     &   WFUNC,VFUNC)
         use constants, only: dp
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,I_V,IPATH(NEL,0:I_V)
         INTEGER nBasisMax(5,*),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE,FMCPR5
         real(dp) RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         real(dp) NRHOII(0:I_V),NRHOIJ(0:I_V,0:I_V)
         INTEGER I_NVMOVE,INPATH(NEL,0:I_V)
         INTEGER INODE(NEL),I_VMOVE
         INTEGER IADJ(0:I_V-1,0:I_V-1),ADJ
         INTEGER NIADJ(0:I_V-1,0:I_V-1),NADJ
         real(dp) VALUE,RAT
         INTEGER ICOUNT
         LOGICAL LCANMOVE
         LOGICAL LACC
         LOGICAL TLOG,TLOG2
         real(dp) WEIGHT,NWEIGHT
         real(dp) WFUNC,VFUNC
         real(dp) RAN2
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
C.. First decide what MC move - either change the moving vertex or
C.. move the moving vertex.
         IF(RAN2(ISEED).LT.0.25_dp) THEN
C.. Our move just involves changing which vertex we move 1 in 4 times
            CALL GENNEWVMOVE(I_V,I_VMOVE,IADJ,ADJ,LCANMOVE,ISEED)
            FMCPR5=VALUE
            RETURN
         ENDIF
         LACC=.FALSE.
         IF(LCANMOVE) THEN
C.. Generate a node to replace the one at I_VMOVE, and put it in INODE
            CALL GENPATHEXCIT(IPATH,I_VMOVE,NEL,NBASIS,NTAY,ISEED,INODE,
     &         IADJ,I_V)
C.. Now insert the node in the path at I_VMOVE.  The path must then be
C.. lexically reordered and I_VMOVE,IADJ and the RHOIJ matrix 
C.. modified appropriately
C.. The relevant elements of IADJ and RHOIJ are copied to the new path's
            CALL MODIFYPATH(IPATH,I_VMOVE,RHOII,RHOIJ,IADJ,INODE,
     &               NRHOII,NRHOIJ,NIADJ,I_NVMOVE,NEL,INPATH,I_V)
            CALL MODIFYADJCONMAT(INPATH,I_NVMOVE,NRHOIJ,NRHOII,NIADJ,
     &               INODE,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,RHOEPS,
     &           NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,I_V,NADJ,I_V-1)
C.. NIADJ now contains the adjacency matrix, and NRHOIJ the connectivity
C.. of the new path
C.. We now must calculate the weight of the new move
            NWEIGHT=WFUNC(INPATH,NRHOII,NRHOIJ,I_V,I_HMAX,
     &            I_P,NWHTAY,NEL,ILOGGING,I_NVMOVE)
            RAT=NWEIGHT*NADJ/(WEIGHT*ADJ)
C.. Check for acceptance
            IF(RAT.GE.1.0_dp) THEN
               LACC=.TRUE.
            ELSE
               IF(RAN2(ISEED).LE.RAT) LACC=.TRUE.
            ENDIF
         ENDIF
         IF(LACC) THEN
C.. We accept the new move.  Copy NRHOII to RHOII, NRHOIJ to RHOIJ etc.
            ADJ=NADJ
            I_VMOVE=I_NVMOVE
            WEIGHT=NWEIGHT
            CALL NECI_ICOPY(NEL*(I_V+1),INPATH,1,IPATH,1)
            CALL NECI_ICOPY(I_V*I_V,NIADJ,1,IADJ,1)
            CALL DCOPY((I_V+1)**2,NRHOIJ,1,RHOIJ,1)
            CALL DCOPY(I_V+1,NRHOII,1,RHOII,1)
         ENDIF        
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.  We now call CALCPATHS to get the path weight
         CALL NECI_ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)   
            IF(TLOG2) CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
         ENDIF   
C.. Now calculate a new value if required
         IF(LACC) THEN
            VALUE=VFUNC(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,ILOGGING,I_VMOVE)/WEIGHT
         ENDIF
         FMCPR5=VALUE
         IF(TLOG) WRITE(10,*) FMCPR5,ICOUNT 
         RETURN
      END

      SUBROUTINE GENNEWVMOVE(I_V,I_VMOVE,IADJ,ADJ,LCANMOVE,ISEED)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I_V,I_VMOVE,IADJ(0:I_V-1,0:I_V-1),ISEED,ADJ
         LOGICAL LCANMOVE
         INTEGER I,J,K,L
         real(dp) ADJ2(I_V-1,I_V-1),ADJ3(I_V-1,I_V-1)
         real(dp) RAN2
         I_VMOVE=INT(RAN2(ISEED)*(I_V-1))+1
         ADJ=0
         DO I=0,I_V-1
            IF(I.NE.I_VMOVE) ADJ=ADJ+IADJ(I,I_VMOVE)
         ENDDO
C.. If the new moving vertex is only connected to one other, it must be
C.. allowed to move, otherwise we need to check it's movable
         IF(ADJ.EQ.1) THEN
            LCANMOVE=.TRUE.
         ELSE
C            IPTEMP(0)=I_VMOVE
C            LCANMOVE=IMCPR4NC(IADJ,I_V,IPTEMP,1).EQ.ADJ 
C.. we form the (I_V-1)th power of the adjacency matrix without the moving node
C.. and check the first row has all non-zeroes (i.e. we can get from the
C.. root to any other node)
            K=0
            DO I=0,I_V-1
               IF(I.NE.I_VMOVE) THEN
                  K=K+1
                  L=0
                  DO J=0,I_V-1
                     IF(J.NE.I_VMOVE) THEN
                        L=L+1
                        ADJ2(K,L)=IADJ(I,J)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            CALL MATPOWER(ADJ2,I_V-1,ADJ3,I_V-2)
            LCANMOVE=.TRUE.
            DO I=1,I_V-1
               IF(ADJ3(1,I).EQ.0) LCANMOVE=.FALSE.
            ENDDO
         ENDIF
         RETURN
      END

      SUBROUTINE GENPATHEXCIT(IPATH,I_VMOVE,NEL,NBASIS,NTAY,ISEED,INODE,
     &         IADJ,I_V)
         use mcpathsismc, only: genrandomexcit
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I_V,NEL,I_VMOVE,NBASIS,NTAY,ISEED,IPATH(NEL,0:I_V)
         INTEGER IADJ(0:I_V-1,0:I_V-1),INODE(NEL)
         INTEGER IEXCITS(0:I_V-1),I,J,IEX,IEXNODE,E
         LOGICAL BR,LISINPATH
         real(dp) RAN2
C.. We first calculate the number of allowed excitations for each vertex
C.. (excluding I_VMOVE) in our path.
         E=NBASIS-NEL
         DO I=0,I_V-1
            IF(I.EQ.I_VMOVE) THEN
               IEXCITS(I)=IEXCITS(I-1)
            ELSE
               IEX=E*NEL+E*NEL*(E-1)*(NEL-1)/4
C.. Subtract out the number of adjacent vertices which are already in
C.. the path
               DO J=0,I_V-1
                  IF(J.NE.I) IEX=IEX-IADJ(J,I)
               ENDDO
               IF(I.GT.0) IEX=IEX+IEXCITS(I-1)
               IEXCITS(I)=IEX
            ENDIF
         ENDDO
C.. IEXCITS(I_V-1) now contains the total number of excitations
C.. available from all the vertices in the path.  Some of these will be
C.. double-counted, but that bias is dealt with elsewhere.
         BR=.TRUE.
         IEX=RAN2(ISEED)*IEXCITS(I_V-1)
         DO I=I_V-1,0,-1
            IF(IEX.LT.IEXCITS(I)) IEXNODE=I
         ENDDO
         DO WHILE(BR)
            CALL GENRANDOMEXCIT(IPATH(1,IEXNODE),NEL,NBASIS,
     &            ABS(NTAY*2),ISEED,INODE)
            BR=LISINPATH(INODE,IPATH,NEL,I_V,I_VMOVE)
         ENDDO
         RETURN
      END

      SUBROUTINE MODIFYADJCONMAT(IPATH,I_NVMOVE,NRHOIJ,NRHOII,NIADJ,
     &               INODE,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,RHOEPS,
     &               NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,I_V,NADJ,
     &               I_VMAX)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,I_V,IPATH(NEL,0:I_V),NIADJ(0:I_V-1,0:I_V-1),NADJ
         INTEGER INODE(NEL),I_P,NBASISMAX(*),G1(*),NBASIS,BRR(*)
         INTEGER NMSH,NMAX,NTAY,I_NVMOVE,I_VMAX
         real(dp) NRHOIJ(0:I_V,0:I_V),NRHOII(0:I_V),RHOEPS
         real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE
         INTEGER J,ICE
         real(dp) RH
         INTEGER IGETEXCITLEVEL
C.. Now generate all the RHOIJ elements for this new node
         NADJ=0
         DO J=0,I_VMAX
            ICE=IGETEXCITLEVEL(INODE,IPATH(1,J),NEL)
            CALL CALCRHO2(INODE,IPATH(1,J),BETA,I_P,NEL,
     &        NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,
     &        RH,NTAY,ICE,ECORE)
            IF(ICE.LE.2) THEN
               ICE=1
            ELSE
               ICE=0
            ENDIF
            NADJ=NADJ+ICE
            IF(ABS(RH).LE.RHOEPS) RH=0.0_dp
            NIADJ(I_NVMOVE,J)=ICE
            NIADJ(J,I_NVMOVE)=ICE
            NRHOIJ(I_NVMOVE,J)=RH
            NRHOIJ(J,I_NVMOVE)=RH
            IF(I_NVMOVE.EQ.J) NRHOII(J)=RH
         ENDDO   
         RETURN
      END

      

      SUBROUTINE GENVALIDPATH(BETA,I_P,IPATH,I_V,NEL,
     &         NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &         ECORE,ISEED,I_VMOVE,OVALUE,IADJ,ADJ,LCANMOVE,
     &         WFUNC,VFUNC,WEIGHT)
         use util_mod, only: NECI_ICOPY
         use constants, only: dp
         use mcpathsismc, only: genrandomexcit
         IMPLICIT NONE
         INTEGER I_V,NEL,I_P,IPATH(NEL,0:I_V)
         INTEGER nBasisMax(5,*),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE
         real(dp) RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER I,K,INODE(NEL),I_VMOVE,ADJ
         INTEGER IADJ(0:I_V-1,0:I_V-1)
         real(dp) RH,OVALUE,WEIGHT
         INTEGER ICOUNT
         LOGICAL BR
         LOGICAL LCANMOVE
         real(dp) VFUNC,WFUNC
         LOGICAL LISINPATH
         real(dp) RAN2
C.. we need to generate a valid path
         DO I=1,I_V-1
            BR=.TRUE.
            K=INT(RAN2(ISEED)*I)
            ICOUNT=0
            DO WHILE(BR)
               CALL GENRANDOMEXCIT(IPATH(1,K),NEL,NBASIS,
     &            ABS(NTAY*2),ISEED,INODE)
               IF(.NOT.LISINPATH(INODE,IPATH,NEL,I,-1)) THEN
                  CALL CALCRHO2(IPATH(1,K),INODE,BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,
     &             RH,NTAY,-1,ECORE)
                  IF(ABS(RH).GT.RHOEPS) BR=.FALSE.
                  ICOUNT=ICOUNT+1
C.. If we cannot find a det that's connected, we give up
                  IF(ICOUNT.GT.10000) BR=.FALSE.
               ENDIF
            ENDDO
            CALL NECI_ICOPY(NEL,INODE,1,IPATH(1,I),1)
            CALL MODIFYADJCONMAT(IPATH,I,RHOIJ,RHOII,IADJ,
     &               INODE,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,RHOEPS,
     &               NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,I_V,ADJ,I)
         ENDDO
         I_VMOVE=I_V-1
         OVALUE=VFUNC(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,ILOGGING,I_VMOVE)
         WEIGHT=WFUNC(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,ILOGGING,I_VMOVE)
         LCANMOVE=.TRUE.
         RETURN
      END
#else
      SUBROUTINE  MCPATHSR5(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,
     &               NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &               NWHTAY,ILOGGING,ECORE, WLRI,WLSI)
        use constants, only: dp
        IMPLICIT NONE
        INTEGER NEL,NI(NEl),I_P,I_VMAX,I_HMAX
        INTEGER nBasisMax(5,*),NMAX,NBASIS,NMSH,NTAY,NWHTAY,ILOGGING
        real(dp) Beta,Brr(*),fck,ALAT(3),UMAT(*),RHOEPS,ECORE,WLRI,WLSI
        integer G1(*)
        character(*), parameter :: t_r="MCPATHSR5"
        call stop_all(t_r,"Routine depreciated. See mcpathsmcmc.F head")
      end subroutine MCPATHSR5


      SUBROUTINE GENNEWVMOVE(I_V,I_VMOVE,IADJ,ADJ,LCANMOVE,ISEED)
        IMPLICIT NONE
        INTEGER :: I_V,I_VMOVE,IADJ(0:I_V-1,0:I_V-1),ISEED,ADJ
        LOGICAL :: LCANMOVE
        character(*), parameter :: t_r="GENNEWVMOVE"

        call stop_all(t_r,"Routine depreciated. See mcpathsmcmc.F head")

      end subroutine GenNewVMove

        
      SUBROUTINE GENPATHEXCIT(IPATH,I_VMOVE,NEL,NBASIS,NTAY,ISEED,INODE,
     &         IADJ,I_V)
        IMPLICIT NONE
        INTEGER :: I_V,NEL,I_VMOVE,NBASIS,NTAY,ISEED,IPATH(NEL,0:I_V)
        INTEGER :: IADJ(0:I_V-1,0:I_V-1),INODE(NEL)
        character(*), parameter :: t_r="GENPATHEXCIT"
        call stop_all(t_r,"Routine depreciated. See mcpathsmcmc.F head")

      END SUBROUTINE GENPATHEXCIT



      SUBROUTINE MODIFYADJCONMAT(IPATH,I_NVMOVE,NRHOIJ,NRHOII,NIADJ,
     &               INODE,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,RHOEPS,
     &               NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,I_V,NADJ,
     &               I_VMAX)
        use constants, only: dp
        IMPLICIT NONE
        INTEGER NEL,I_V,IPATH(NEL,0:I_V),NIADJ(0:I_V-1,0:I_V-1),NADJ
        INTEGER INODE(NEL),I_P,NBASISMAX(*),G1(*),NBASIS,BRR(*)
        INTEGER NMSH,NMAX,NTAY,I_NVMOVE,I_VMAX
        real(dp) NRHOIJ(0:I_V,0:I_V),NRHOII(0:I_V),RHOEPS
        real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE
        character(*), parameter :: t_r="MODIFYADJCONMAT"
        call stop_all(t_r,"Routine depreciated. See mcpathsmcmc.F head")

      END SUBROUTINE MODIFYADJCONMAT


      SUBROUTINE GENVALIDPATH(BETA,I_P,IPATH,I_V,NEL,
     &         NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &         RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &         ECORE,ISEED,I_VMOVE,OVALUE,IADJ,ADJ,LCANMOVE,
     &         WFUNC,VFUNC,WEIGHT)
        use constants, only: dp
        IMPLICIT NONE
        INTEGER I_V,NEL,I_P,IPATH(NEL,0:I_V)
        INTEGER nBasisMax(5,*),G1,NBASIS,BRR(NBASIS),NMAX
        INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
        real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE
        real(dp) RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
        real(dp) VFUNC,WFUNC
        LOGICAL LCANMOVE
        real(dp) OVALUE,WEIGHT
        INTEGER IADJ(0:I_V-1,0:I_V-1),ADJ
        INTEGER I_VMOVE
        character(*), parameter :: t_r="GENVALIDPATH"
        call stop_all(t_r,"Routine depreciated. See mcpathsmcmc.F head")
      END SUBROUTINE GENVALIDPATH
      
      FUNCTION FMCPR5(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,
     &   ECORE,ISEED,I_VMOVE,VALUE,WEIGHT,IADJ,ADJ,LCANMOVE,
     &   WFUNC,VFUNC)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,I_V,IPATH(NEL,0:I_V)
         INTEGER nBasisMax(5,*),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE,FMCPR5
         real(dp) RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER I_VMOVE
         INTEGER IADJ(0:I_V-1,0:I_V-1),ADJ
         LOGICAL LCANMOVE
         real(dp) WEIGHT,VALUE
         real(dp) WFUNC,VFUNC
         character(*), parameter :: t_r="FMCPR5"
         FMCPR5=0.0_dp
        call stop_all(t_r,"Routine depreciated. See mcpathsmcmc.F head")

       END FUNCTION FMCPR5
#endif
      SUBROUTINE MODIFYPATH(IPATH,I_VMOVE,RHOII,RHOIJ,IADJ,INODE,
     &               NRHOII,NRHOIJ,NIADJ,I_NVMOVE,NEL,INPATH,I_V)
         use util_mod, only: NECI_ICOPY
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,I_V,IPATH(NEL,0:I_V),IADJ(0:I_V-1,0:I_V-1)
         INTEGER INPATH(NEL,0:I_V),NIADJ(0:I_V-1,0:I_V-1)
         INTEGER I_VMOVE,I_NVMOVE,INODE(NEL)
         real(dp) RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         real(dp) NRHOII(0:I_V),NRHOIJ(0:I_V,0:I_V)

C.. In fact, the path does not need to be lexically ordered, so we just
C.. copy
         CALL NECI_ICOPY(NEL*(I_V+1),IPATH,1,INPATH,1)
         CALL NECI_ICOPY(NEL,INODE,1,INPATH(1,I_VMOVE),1)
         CALL NECI_ICOPY(I_V*I_V,IADJ,1,NIADJ,1)
         CALL DCOPY((I_V+1)**2,RHOIJ,1,NRHOIJ,1)
         CALl DCOPY(I_V+1,RHOII,1,NRHOII,1)
         I_NVMOVE=I_VMOVE
         RETURN
      END

C.. calculates B=A**N where B and A are LxL matrices
      SUBROUTINE MATPOWER(A,L,B,N)
         use constants, only: dp
         real(dp) A(L,L),B(L,L)
         real(dp) C(L,L,0:1)
         INTEGER N,N2
         INTEGER I,S
         LOGICAL LM
         N2=N
         S=0
         LM=BTEST(N2,0)
         N2=N2/2
         IF(LM) THEN
            CALL DCOPY(L*L,A,1,B,1)
         ELSE
            B(1:L,1:L)=0
            DO I=1,L
               B(I,I)=1.0_dp
            ENDDO
         ENDIF
         DO WHILE(N2.NE.0)
            LM=BTEST(N2,0)
            N2=N2/2
            CALL DGEMM('N','N',L,L,L,1.0_dp,C(1,1,S),L,C(1,1,S),L,
     &        0.0_dp, C(1,1,1-S),L)
            S=1-S
            IF(LM) THEN
               CALL DGEMM('N','N',L,L,L,1.0_dp,C(1,1,S),L,B,L,0.0_dp,
     &            C(1,1,1-S),L)
               CALL DCOPY(L*L,C(1,1,1-S),1,B,1)
            ENDIF
         ENDDO
      END
      
      FUNCTION GETNLOOPS(I_V,NBASIS,NEL)
        use util_mod, only: factrl
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I_V,NBASIS,NEL,I
         real(dp) NL,FNCR,NL2
         real(dp) N,E,GETNLOOPS
         N=NEL
         E=NBASIS-NEL
         IF(I_V.EQ.2) THEN
            NL=N*E*(4+(N-1)*(E-1))/4
         ELSEIF(I_V.EQ.3) THEN
            NL=N*E*(8*(N+E-2)+(N-1)*(E-1)*(12*(N+E-2)+
     &         (N-2)*(N-3)+(E-2)*(E-3)+
     &         (N-2)*(E-2)*(16 +4*(N+E-6)+(N-3)*(E-3))))/16.0_dp
         ELSEIF(NBASIS.GT.30) THEN
            NL=(FNCR(NBASIS,NEL)**I_V)/FACTRL(I_V)
         ELSE
            NL=1.0_dp
            NL2=NEL*(NBASIS-NEL)
     &         +(NBASIS-NEL)*(NBASIS-NEL-1)*NEL*(NEL-1)/4
            DO I=2,I_V
               NL=NL*NL2
               NL2=NL2-1
            ENDDO
         ENDIF
         GETNLOOPS=NL
         RETURN
      END
      
      LOGICAL FUNCTION LISINPATH(NI,IPATH,NEL,I_V,I_NV)
         IMPLICIT NONE
         INTEGER I,I_NV
         INTEGER NEL,I_V,NI(NEL),IPATH(NEL,0:I_V-1)
         INTEGER ICMPDETS
         LISINPATH=.FALSE.
         DO I=0,I_V-1
            IF(I.NE.I_NV.AND.ICMPDETS(NI,IPATH(1,I),NEL).EQ.0)
     &          LISINPATH=.TRUE.
         ENDDO
         RETURN
      END

      INTEGER FUNCTION IISINPATH(NI,IPATH,NEL,I_V,I_NV)
         IMPLICIT NONE
         INTEGER I,I_NV
         INTEGER NEL,I_V,NI(NEL),IPATH(NEL,0:I_V-1)
         INTEGER ICMPDETS
         IISINPATH=-1
         DO I=0,I_V-1
            IF(I.NE.I_NV.AND.ICMPDETS(NI,IPATH(1,I),NEL).EQ.0)
     &          IISINPATH=I
         ENDDO
         RETURN
      END

      FUNCTION FMCPR4WFUNC(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,ILOGGING,I_VMOVE)
         use constants, only: dp
         INTEGER NEL,I_V,I_HMAX,I_P,IPATH(NEL,0:I_V),ILOGGING
         INTEGER NWHTAY,I_VMOVE
         real(dp) RHOII(0:I_V),RHOIJ(0:I_V,0:I_V),FMCPR4WFUNC
         INTEGER I,J
         LOGICAL L
         FMCPR4WFUNC=1.0_dp
         DO I=0,I_V-1
            L=.FALSE.
            DO J=0,I_V-1
               IF(I.NE.J) THEN
                  IF(ABS(RHOIJ(I,J)).GT.0.0_dp) L=.TRUE.
               ENDIF
            ENDDO
            IF(.NOT.L) FMCPR4WFUNC=0.0_dp
         ENDDO

         ! Avoid warnings
         ipath=ipath; rhoii(0)=rhoii(0); i_p=i_p; nwhtay=nwhtay
         ilogging=ilogging; i_vmove=i_vmove; i_hmax=i_hmax
      END
      
      FUNCTION FMCPR4VFUNC(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,ILOGGING,I_VMOVE)
         use constants, only: dp
         INTEGER NEL,I_V,I_HMAX,I_P,IPATH(NEL,0:I_V),ILOGGING
         INTEGER NWHTAY,I_VMOVE
         real(dp) RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         real(dp) CALCPATHS,FMCPR4VFUNC
         FMCPR4VFUNC=CALCPATHS(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,I_V,ILOGGING)

         ! Avoid warnings
         i_vmove=i_vmove
      END
      
      FUNCTION FMCPR4UFUNC(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,NWHTAY,NEL,ILOGGING,I_VMOVE)
         use constants, only: dp
         INTEGER NEL,I_V,I_HMAX,I_P,IPATH(NEL,0:I_V),ILOGGING
         INTEGER NWHTAY
         real(dp) RHOII(0:I_V),RHOIJ(0:I_V,0:I_V),FMCPR4UFUNC
         FMCPR4UFUNC=1.0_dp

         ! Avoid warnings
         ipath=ipath; rhoii(0)=rhoii(0); i_v=i_v; i_hmax=i_hmax
         i_p=i_p; nwhtay=nwhtay;nel=nel; ilogging=ilogging
         i_vmove=i_vmove; rhoij=rhoij
      END
      RECURSIVE INTEGER FUNCTION IMCPR4NC(IADJ,I_V,IPATH,IND)
     &    RESULT (IMCPR4NCRES)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I,J,INODE,ITOT
         INTEGER I_V,IADJ(0:I_V-1,0:I_V-1),IPATH(0:I_V-1),IND
         LOGICAL BR
C.. Go through all the nodes currently in the path, and through each
C.. possible attachment for each node
         IF(IND.EQ.1) THEN
C.. We need to recurse through nodes connected to us, removing them
            BR=.TRUE.
            I=0
            DO WHILE(BR)
               IF(IADJ(I,IPATH(0)).EQ.1.AND.IPATH(0).NE.I) THEN
                  BR=.FALSE.
                  IADJ(I,IPATH(0))=0
                  IADJ(IPATH(0),I)=0
                  IMCPR4NCRES=IMCPR4NC(IADJ,I_V,IPATH,1)
                  IPATH(1)=I
                  IADJ(I,I)=0
                  IMCPR4NCRES=IMCPR4NCRES+IMCPR4NC(IADJ,I_V,IPATH,2)
                  IADJ(I,IPATH(0))=1
                  IADJ(IPATH(0),I)=1
                  IADJ(I,I)=1
               ENDIF
               I=I+1
               IF(I.GE.I_V) BR=.FALSE.
            ENDDO
         ELSE
            IF(IPATH(IND).EQ.0) THEN
               IMCPR4NCRES=1
               RETURN
            ENDIF
            ITOT=0
            INODE=IPATH(IND)
            J=0
            BR=.TRUE.
            DO WHILE(BR)
C.. If there's a connection and we haven't been to that node before
               IF(IADJ(INODE,J).NE.0.AND.IADJ(J,J).NE.0) THEN
                  IPATH(IND)=J
                  IADJ(J,J)=0
                  IADJ(INODE,J)=0
                  IADJ(J,INODE)=0
C.. add that node to the path, and recurse
                  ITOT=ITOT+IMCPR4NC(IADJ,I_V,IPATH,IND+1)
                  IF(ITOT.EQ.1) BR=.FALSE.
                  IADJ(J,J)=1
                  IADJ(INODE,J)=1
                  IADJ(J,INODE)=1
               ENDIF
               J=J+1
               IF(J.GE.I_V) BR=.FALSE.
            ENDDO
            IMCPR4NCRES=ITOT
         ENDIF
         RETURN
      END
