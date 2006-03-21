
C.. A stochastic diagonaliser - generate an optimal set of determinants
C.. (a graph) using the weight of the graph to choose the dets.
      SUBROUTINE  STOCHDIAG(NI,BETA,I_P,I_VMAX,NEL,NBASISMAX,G1,
     &               NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &               NWHTAY,ILOGGING,ECORE,WLRI,WLSI,DBETA,DLWDB)
         IMPLICIT NONE
         INTEGER NI(NEL),I_P,NMSH,NMAX,NTAY,NWHTAY
         INTEGER ILOGGING,ISUB,I,J,K
         REAL*8 ALAT(3),RHOEPS,BETA,UMAT(*),ZIA(*),FCK(*)
         REAL*8 ECORE,WLRI,WLSI
         INTEGER G1(5,NBASIS),NBASISMAX(5,2),BRR(*)
         INTEGER I_VMAX,NEL,NBASIS
         INTEGER ISEED
         REAL*8 DBETA,DLWDB
         REAL*8 GETHELEMENT2
         INTEGER NMEM(*),NMEMLEN
         POINTER (IP_NMEM,NMEM)

         INTEGER FSEEK
         INTEGER IPATH(NEL,0:I_VMAX)
         REAL*8 RHOII(0:I_VMAX),RHOIJ(0:I_VMAX,0:I_VMAX)
         CALL TISET('STOCHDIAG ',ISUB)
         IF(BTEST(ILOGGING,2))
     &      OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
         IF(BTEST(ILOGGING,9))
     &      OPEN(12,FILE="VERTEXMC",STATUS="UNKNOWN")
C.. Set the first node to I_I
         CALL ICOPY(NEL,NI,1,IPATH(1,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,
     &            BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(0)=RH
         RHOIJ(0,0)=RHOII(0)
         WLRI=LOG(RHOII(0))
         HIJS(0)=GETHELEMENT2(NI,NI,NEL,NBASISMAX,
     &      G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,0,ECORE)












         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
         IF(I_HMAX.NE.-7) THEN
C.. we need to count the excitations
C            CALL GETEXCITCOUNTSYM(NI,NEL,NBASIS,G1,NBASISMAX,
C     &            KSYM,BRR,IEXCITS)
         ELSE
            IEXCITS=0
            ISEED=7
            CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IEXCITS,
     &      ISEED,INODE2)
            TOTAL=0.D0
            IF(DBETA.LT.0.D0) THEN
               DLWDB=0.D0
            ENDIF
         ENDIF
         I_VM2=I_VMAX
         IF(I_HMAX.EQ.-7) THEN
C.. setup spin excitation generator
            CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
     &               .TRUE.,NMEMLEN,INODE2,ISEED,IEXCITS,0)
            CALL MEMORY(IP_NMEM,NMEMLEN,'NMEM')
            NMEM(1)=0
            CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
     &               .TRUE.,NMEMLEN,INODE2,ISEED,IEXCITS,0)
            I_VM2=2
         ENDIF
         DO I_V=0,I_VMAX
            IVL(I_V)=0
            WGHTT(I_V)=0.D0
            WGHTP(I_V)=0.D0
            WGHTM(I_V)=0.D0
            NTREES(I_V)=0
            NNTP(I_V)=0
            NNTM(I_V)=0
         ENDDO
         DO I_V=2,I_VM2
            WRITE(STR,"(A,I5)") "FMCPR",I_V
            CALL TISET(STR,ISUB2)
            L=0
            LT=0
            BTABLE(0)=0
            F(I_V)=0.D0
            FSQ(I_V)=0.D0
            ISEED=7
            IF(I_V.LE.IEXCITS+1) THEN
C.. Do importance sampling with probabilities.  No norm needed
               DLWDB3=0.D0
               OSI=0.D0
               OWEIGHT=0.D0
               IOCLS=0
               I_VCUR=1
               NTR=0
               ITREE=0
               DO ICOUNT=1,NWHTAY
                  DLWDB2=0.D0
                  IF(I_HMAX.EQ.-3) THEN
                  FF=FMCPR4B(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ELSEIF(I_HMAX.EQ.-4) THEN
                     FF=FMCPR4C(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS)
                  ELSE
C.. see if we want to increase the number of vertices
                     I_OVCUR=I_VCUR
                     R=RAN(ISEED)
                     IF(R.LT.0.1D0) THEN
                        IF(I_VCUR.GT.1) I_VCUR=I_VCUR-1
                     ELSEIF(R.LT.0.2D0) THEN
                        IF(I_VCUR.LT.I_VMAX) I_VCUR=I_VCUR+1
                     ENDIF         
                     FF=FMCPR4D2(NI,BETA,I_P,IPATH,I_VCUR,NEL,NBASISMAX,
     &                  G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &                  RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &                  ECORE,ISEED,KSYM,DBETA,DLWDB2,HIJS,NMEM,
     &                  OSI,ODLWDB,OPROB,I_OVCUR,IOCLS,ITREE,OWEIGHT)
                     IVL(I_VCUR)=IVL(I_VCUR)+1
                     IVL(0)=IVL(0)+1
                  ENDIF
                  F(I_V)=F(I_V)+FF
                  FSQ(I_V)=FSQ(I_V)+FF*FF
                  DLWDB3=DLWDB3+DLWDB2
C.. some sort of compiler bug needs this
                  F(I_V)=F(I_V)+0.D0
                  DLWDB3=DLWDB3+0.D0
                  NTR=NTR+ITREE
                  IF(ITREE.EQ.1) THEN
                     WGHTT(I_VCUR)=WGHTT(I_VCUR)+OWEIGHT
                     NTREES(I_VCUR)=NTREES(I_VCUR)+1
                     WGHTT(0)=WGHTT(0)+OWEIGHT
                     NTREES(0)=NTREES(0)+1
                  ELSE
                     IF(OWEIGHT.GT.0.D0) THEN
                        WGHTP(I_VCUR)=WGHTP(I_VCUR)+OWEIGHT
                        NNTP(I_VCUR)=NNTP(I_VCUR)+1
                        WGHTP(0)=WGHTP(0)+OWEIGHT
                        NNTP(0)=NNTP(0)+1
                     ELSE
                        WGHTM(I_VCUR)=WGHTM(I_VCUR)-OWEIGHT
                        NNTM(I_VCUR)=NNTM(I_VCUR)+1
                        WGHTM(0)=WGHTM(0)-OWEIGHT
                        NNTM(0)=NNTM(0)+1
                     ENDIF
                  ENDIF
C                  WRITE(32,*) FF,DLWDB2,F(I_V),DLWDB3
                  IF(BTEST(ILOGGING,9).AND.MOD(ICOUNT,1000).EQ.0) THEN
                    WRITE(12,"(I3,I10,8G19.7)") I_V,ICOUNT/1000,
     &                   F(I_V)/ICOUNT,
     &                   SQRT(FSQ(I_V)/ICOUNT-(F(I_V)/ICOUNT)**2),
     &                   DLWDB3/ICOUNT,
     &                   (DLWDB+DLWDB3/ICOUNT)/(F(I_V)/ICOUNT+TOTAL),
     &                   (NTR+0.D0)/ICOUNT,
     &                   (NTREES(0)+0.D0)/IVL(0),
     &                   (NNTP(0)+0.D0)/IVL(0),
     &                   (NNTM(0)+0.D0)/IVL(0)
                     CALL FLUSH(12)
                  ENDIF
               ENDDO
            ENDIF

            F(I_V)=F(I_V)/NWHTAY
            FSQ(I_V)=FSQ(I_V)/NWHTAY
            DLWDB=DLWDB+DLWDB3/NWHTAY
            STD=SQRT(ABS(FSQ(I_V)-F(I_V)*F(I_V)))
            CALL TIHALT(STR,ISUB2)
            TOTAL=TOTAL+F(I_V)
c            WRITE(6,*) I_V,F(I_V),TOTAL,TIGET(STR,ISUB2),L,LT
            IF(TLOG)
     &         WRITE(11,"(I12,2G25.16,F19.7,I12,2F19.7)") I_V,F(I_V),
     &            TOTAL,TIGET(STR,ISUB2),L,STD,DLWDB3/NWHTAY
C ,F(I_V),FSQ(I_V)
            IF(TLOG.AND.I_HMAX.EQ.-7) THEN
               WRITE(STR2,"(A,I5,A)") "(A,",I_VMAX+1,"I)"
               WRITE(11,STR2) "GRAPHS(V)",(IVL(I),I=0,I_VMAX)
               WRITE(11,STR2) "TREES(V)",(NTREES(I),I=0,I_VMAX)
               WRITE(11,STR2) "NON-TR+(V)",(NNTP(I),I=0,I_VMAX)
               WRITE(11,STR2) "NON-TR-(V)",(NNTM(I),I=0,I_VMAX)
               WRITE(STR2,"(A,I5,A)") "(A,",I_VMAX+1,"G)"
               WRITE(11,STR2) "WGHTT(V)",(WGHTT(I),I=0,I_VMAX)
               WRITE(11,STR2) "WGHT+(V)",(WGHTP(I),I=0,I_VMAX)
               WRITE(11,STR2) "WGHT-(V)",(WGHTM(I),I=0,I_VMAX)
            ENDIF
            CALL FLUSH(11)
            IF(LISNAN(F(I_V))) THEN
C.. save all log files
               ITIME=TIMEC()
               CALL FLUSH(11)
               CALL FLUSH(10)
               CALL LOGNAN(NI,NEL,BETA,ITIME)
               WRITE(6,*) "WARNING: nan found at time",ITIME
               WRITE(6,"(A,$)") "  nan det=("
               DO K=1,NEL
                  WRITE(6,"(I3,A,$)") NI(K),","
               ENDDO
               WRITE(6,"(A,$)") "),"
            ENDIF
         ENDDO
         IF(I_HMAX.EQ.-7) THEN
C.. setup spin excitation generator
            CALL FREEM(IP_NMEM)
         ENDIF
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,2)) CLOSE(10)
         IF(BTEST(ILOGGING,9)) CLOSE(12)
         WLSI=LOG(TOTAL)
         DLWDB=DLWDB/TOTAL
         CALL TIHALT('MCPATHSR4 ',ISUB)
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
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
C.. Pick random elements of the path and generate excitations
         DO I=1,I_V-1
            K=INT(RAN(ISEED)*I)
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
            CALL WRITE_PATH(10,IPATH,NEL,I_V)
            IF(TLOG2) CALL WRITE_RHOMATRIX(10,RHOIJ,NEL,I_V)
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
         IF(IEXLEVEL.GT.2)
     &    STOP "Cannot handle more than double excitations."
         IF(IEXLEVEL.LT.2) IEXL2=IEXLEVEL
         IF(IEXLEVEL.EQ.2) THEN
            IEX=RAN(ISEED)*(NBASIS-NEL)*NEL*
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
               I=INT(RAN(ISEED)*NEL)+1
               IF(NJ(I).EQ.NI(I)) BR=.FALSE.
            ENDDO
            BR=.TRUE.
            DO WHILE (BR)
               BR=.FALSE.
               J=INT(RAN(ISEED)*NBASIS)+1
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
                     IF(ISUHFDET(IPATH(1,I),NEL,NBASISMAX)) THEN
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
            CALL WRITE_PATH(10,IPATH,NEL,I_V)
            IF(TLOG2) CALL WRITE_RHOMATRIX(10,RHOIJ,NEL,I_V)
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
            PP=RAN(ISEED)*XIJ(ICURNODE,ICURNODE)
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
            CALL WRITE_PATH(10,IPATH,NEL,I_V)
            IF(TLOG2) CALL WRITE_RHOMATRIX(10,RHOIJ,NEL,I_V)
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
C         WRITE(56,*) "N"
C         WRITE(57,*) "N"
C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
C.. doesn't bother calculating the RHO_JJ.
C         CALL GENSYMDETSSDN(NI,KSYM,NEL,G1,BRR,NBASIS,0,
C     &         IC,NBASISMAX,
C     &      BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,
C     &      RHOEPS)
CC         CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IC,
CC     &      ISEED,INODE2)
         CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
     &         .FALSE.,NMEM,INODE2,ISEED,IC,0)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         NLIST=1
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
            R=RAN(ISEED)
            NEXNODE=0
            DO WHILE(R.GE.0)
               NEXNODE=NEXNODE+1
               R=R-PEXCIT(NEXNODE)
            ENDDO
C            WRITE(56,*) NEXNODE
CC            CALL GENRANDOMSPINEXCIT(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
CC     &      NBASISMAX,IC,ISEED,INODE2)

C.. set the current to be the pivot
            IP_CUREX=PVERTMEMS(NEXNODE-1)


C.. setup spin excitation generator
            CALL GENRANDSYMEXCITIT(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
     &         NBASISMAX,.FALSE.,CUREX,INODE2,ISEED,IC,0)
            CALL ICOPY(NEL,INODE2,1,INODE,1)
C            CALL WRITEDET(6,INODE,NEL,.TRUE.)

C            CALL GENRANDOMEXCITSYM(INODE,NEL,NBASIS,G1,2,ISEED,
C     &         NBASISMAX,KSYM,INODE2)
C.. If the new node is not in the path, add it.
            IF(.NOT.LISINPATH(INODE,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL ICOPY(NEL,INODE,1,IPATH(1,I_VNEXT),1)
C.. Count the number of connections from this node (RP=0)
C               NLIST=0
C               CALL GENSYMDETSSDN(INODE2,KSYM,NEL,G1,BRR,NBASIS,LSTE,
C     &            NLIST,NBASISMAX,BETA,I_P,NMSH,FCK,ZIA,NMAX,
C     &            ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,RHOEPS)
               CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
     &               .TRUE.,NEWEXLEN,INODE2,ISEED,IC,0)
               CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
               NEWEX(1)=0
               PVERTMEMS(I_VNEXT)=IP_NEWEX
C               WRITE(6,*) I_VNEXT,IP_NEWEX
               CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
     &               .FALSE.,NEWEX,INODE2,ISEED,IC,0)
CC               CALL GENRANDOMSPINEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,IC,
CC     &            ISEED,INODE2)

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
                     IF(ISUHFDET(IPATH(1,I),NEL,NBASISMAX)) THEN
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
               FMCPR4D=0.D0
               RETURN
            ENDIF
         ENDDO
         DO I=1,I_V-1
C            WRITE(6,*) PVERTMEMS(I)
            CALL FREEM(PVERTMEMS(I))
C            WRITE(6,*) "N"
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)
         IF(TLOG) THEN
            CALL WRITE_PATH(10,IPATH,NEL,I_V)
            IF(TLOG2) CALL WRITE_RHOMATRIX(10,RHOIJ,NEL,I_V)
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
         INTEGER KSYM(5),G1(5,NBASIS),NBASISMAX(5,2)
         LOGICAL LCHKSYM
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
                  I=INT(RAN(ISEED)*NEL)+1
                  IF(NJ(I).EQ.NI(I)) BR=.FALSE.
               ENDDO
               BR=.TRUE.
               DO WHILE (BR)
                  BR=.FALSE.
                  J=INT(RAN(ISEED)*NBASIS)+1
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
         INTEGER NEL,INODE(NEL),NBASIS,G1(5,NBASIS),NBASISMAX(5,2)
         INTEGER KSYM(5),ICOUNT,BRR(NBASIS)
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
         INTEGER KSYM(5),NEX1A,NEX1B,NEX2A,NEX2B,G1(5,*)
         INTEGER NEXAB
         REAL*8 R,IEX
         LOGICAL ISUHFDET
C         CALL WRITEDET(56,NI,NEL,.TRUE.)
C         CALL FLUSH(56)
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
            CALL ICOPY(NEL,NI,1,NJ,1)
         NB2=NBASIS/2
         NB=(NEL-KSYM(4))/2
         NA=NEL-NB
         NEX1A=NA*(NB2-NA)
         NEX1B=NB*(NB2-NB)
         IF(ISUHFDET(NI,NEL,NBASISMAX)) THEN
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
         R=RAN(ISEED)
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
         INTEGER NI(NEL),NEL,NSPIN,IEX,ISEED,IEL,G1(5,*)
         INTEGER NBASIS
         LOGICAL BR
         BR=.TRUE.
         DO WHILE(BR)
            IEL=RAN(ISEED)*NEL+1
            IF(G1(4,NI(IEL)).EQ.NSPIN.AND.IEL.NE.IEX) BR=.FALSE.
         ENDDO
         IEX=IEL
         RETURN
      END
      SUBROUTINE FINDNEWELECSPIN(NI,NEL,NSPIN,G1,ISEED,NBASIS,IEX,IEX1)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,NSPIN,IEX,ISEED,IEL,IEX1,G1(5,*),NBASIS
         LOGICAL BR
         INTEGER I
         BR=.TRUE.
         DO WHILE(BR)
            IEL=RAN(ISEED)*NBASIS+1
            IF(G1(4,IEL).EQ.NSPIN) THEN
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
      LOGICAL FUNCTION ISUHFDET(NI,NEL,NBASISMAX)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASISMAX(5,6)
         INTEGER I
         ISUHFDET=.FALSE.
         IF(NBASISMAX(4,5).EQ.1) THEN
            ISUHFDET=.TRUE.
            DO I=1,NEL
               IF(NI(I).NE.I) ISUHFDET=.FALSE.
            ENDDO
         ENDIF
C         ISUHFDET=.FALSE.
         RETURN
      END


C.. Like MCPATHSR4C, this selects a n-vertex graph.  it chooses which
Cvertex to excite from randomly however.
C.. 
C.. This does a markov chain monte carlo also - the probabilities are 
C.. weighted with the values of S' (the double-counting corrected
C.. weight.
      REAL*8 FUNCTION FMCPR4D2(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ISEED,KSYM,DBETA,DLWDB,HIJS,NMEM,OSI,ODLWDB,OPROB,
     &   I_OVCUR,IOCLS,ITREE,OWEIGHT)
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
         REAL*8 ORHOII(0:I_V),ORHOIJ(0:I_V,0:I_V)
         INTEGER OIPATH(NEL,0:I_V)
         REAL*8 OHIJS(0:I_V),OXIJ(0:I_V-1,0:I_V-1)
         REAL*8 OSI,ODLWDB,OPROB,PR,R2,WEIGHT,OWEIGHT
         INTEGER I_OVCUR,IOCLS,ITREE
         LOGICAL LISNAN,ISVALIDDET
         INCLUDE 'vmc.inc'
C         WRITE(6,*) "START:",LOC(NMEM)
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
C.. Take a copy of the old path and rho matrix etc.
         CALL DCOPY(I_V+1,RHOII,1,ORHOII,1)
         CALL DCOPY((I_V+1)**2,RHOIJ,1,ORHOIJ,1)
         CALL DCOPY(I_V+1,HIJS,1,OHIJS,1)
         CALL DCOPY(I_V*I_V,XIJ,1,OXIJ,1)
         CALL ICOPY(NEL*(1+I_V),IPATH,1,OIPATH,1)


C         WRITE(56,*) "N"
C         WRITE(57,*) "N"
C.. Count the number of adjacent nodes to us in IC.  If RP==0, it
C.. doesn't bother calculating the RHO_JJ.
C         CALL GENSYMDETSSDN(NI,KSYM,NEL,G1,BRR,NBASIS,0,
C     &         IC,NBASISMAX,
C     &      BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,
C     &      RHOEPS)
CC         CALL GENRANDOMSPINEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,IC,
CC     &      ISEED,INODE2)
         CALL GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
     &         .FALSE.,NMEM,INODE2,ISEED,IC,0)
C.. Diagonal X elements contain the normalization of that node
C.. In this case that is the number of symmetric adjacent nodes

         XIJ(0,0)=IC
         NLIST=1
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
            R=RAN(ISEED)
            NEXNODE=0
            DO WHILE(R.GE.0)
               NEXNODE=NEXNODE+1
               R=R-PEXCIT(NEXNODE)
            ENDDO
C            WRITE(56,*) NEXNODE
CC            CALL GENRANDOMSPINEXCIT(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
CC     &      NBASISMAX,IC,ISEED,INODE2)

C.. set the current to be the pivot
            IP_CUREX=PVERTMEMS(NEXNODE-1)


C.. setup spin excitation generator
            CALL GENRANDSYMEXCITIT(IPATH(1,NEXNODE-1),NEL,G1,NBASIS,
     &         NBASISMAX,.FALSE.,CUREX,INODE2,ISEED,IC,0)
            CALL ICOPY(NEL,INODE2,1,INODE,1)
C            CALL WRITEDET(6,INODE,NEL,.TRUE.)

C            CALL GENRANDOMEXCITSYM(INODE,NEL,NBASIS,G1,2,ISEED,
C     &         NBASISMAX,KSYM,INODE2)
C.. If the new node is not in the path, add it.
            IF(.NOT.ISVALIDDET(INODE,NEL)) THEN
               WRITE(6,*) "INVALID DET"
               STOP "INVALID DET"
            ENDIF
            IF(.NOT.LISINPATH(INODE,IPATH,NEL,I_VNEXT,-1)) THEN
               CALL ICOPY(NEL,INODE,1,IPATH(1,I_VNEXT),1)
C.. Count the number of connections from this node (RP=0)
C               NLIST=0
C               CALL GENSYMDETSSDN(INODE2,KSYM,NEL,G1,BRR,NBASIS,LSTE,
C     &            NLIST,NBASISMAX,BETA,I_P,NMSH,FCK,ZIA,NMAX,
C     &            ALAT,UMAT,NTAY,ECORE,X,RP,0.D0,RHOEPS)
               CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
     &               .TRUE.,NEWEXLEN,INODE2,ISEED,IC,0)
               CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
C               WRITE(6,*) "M",IP_NEWEX
               NEWEX(1)=0
               PVERTMEMS(I_VNEXT)=IP_NEWEX
C               WRITE(6,*) I_VNEXT,IP_NEWEX
               CALL GENRANDSYMEXCITIT(INODE,NEL,G1,NBASIS,NBASISMAX,
     &               .FALSE.,NEWEX,INODE2,ISEED,IC,0)
CC               CALL GENRANDOMSPINEXCIT(INODE,NEL,G1,NBASIS,NBASISMAX,IC,
CC     &            ISEED,INODE2)

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
                     IF(ISUHFDET(IPATH(1,I),NEL,NBASISMAX)) THEN
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
               FMCPR4D2=0.D0
               RETURN
            ENDIF
         ENDDO
         DO I=1,I_V-1
C            WRITE(6,*) "D",PVERTMEMS(I)
            CALL FREEM(PVERTMEMS(I))
C            WRITE(6,*) "N"
         ENDDO
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
         TLOG3=BTEST(ILOGGING,6)
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.
         CALL ICOPY(NEL,NI,1,IPATH(1,I_V),1)

C.. Now get the parameters for the new path 

C.. GETPATHPROB gives us the probability of generating the path
         RH=GETPATHPROB2(XIJ,I_V)
C.. CALCPATHS gives us the contribution of the path
         ICLS=0
         WEIGHT=CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &         I_P,0.D0,NEL,I_V,ILOGGING,DBETA,DLWDB,HIJS,ICLS)
         FMCPR4D2=WEIGHT
         IF(OSI.EQ.0.D0) THEN
C.. this is the first time round, so we automatically accept the new
C.. configuration
            ODLWDB=DLWDB
            OSI=FMCPR4D2
            OPROB=RH
            I_OVCUR=I_V
            ICLS=0
            IOCLS=ICLS
            OWEIGHT=WEIGHT
         ELSE
C.. Accept the new configuration with prob (|S'(new)|/|S'(old)|)*(OPROB/RH)  

C.. acc(A->B)/acc(B->A) = p(B) gen(B->A) / (p(A) gen(A->B))
C..                     = p(B) gen(A) / (p(A) gen(B))

            PR=(OPROB/RH)*ABS(FMCPR4D2/OSI)
            IF(RH.EQ.0.D0.OR.OSI.EQ.0.D0) PR=0.D0
            R2=RAN(ISEED)
            IF(LISNAN(FMCPR4D2)) THEN
               WRITE(60,*) FMCPR4D2,DLWDB,RH
               CALL WRITE_PATH(60,IPATH,NEL,I_V)
               CALL WRITE_RHOMATRIX(60,RHOIJ,NEL,I_V)
               CALL WRITE_XMATRIX(60,XIJ,I_V)
            ENDIF
            IF(PR.LT.R2) THEN
C.. reject the new, and copy in the old
               DLWDB=ODLWDB
               FMCPR4D2=OSI
               RH=OPROB
               I_V=I_OVCUR
               ICLS=IOCLS
               WEIGHT=OWEIGHT
C               CALL DCOPY(I_V+1,ORHOII,1,RHOII,1)
C               CALL DCOPY((I_V+1)**2,ORHOIJ,1,RHOIJ,1)
C               CALL DCOPY(I_V+1,OHIJS,1,HIJS,1)
C               CALL DCOPY(I_V*I_V,OXIJ,1,XIJ,1)
C               CALL ICOPY(NEL*(1+I_V),OIPATH,1,IPATH,1)
            ELSE
C               WRITE(6,*) "ACCEPT",PR
               ODLWDB=DLWDB
               OSI=FMCPR4D2
               OPROB=RH
               I_OVCUR=I_V
               ICLS=0
               IOCLS=ICLS
               OWEIGHT=WEIGHT
            ENDIF
         ENDIF
         IF(TLOG) WRITE(10,"(3E25.16, I7)") FMCPR4D2,RH,DLWDB,ICLS
         IF(RH.GT.0.D0.AND.FMCPR4D2.NE.0.D0) THEN
C.. Unbias the sum 
            DLWDB=DLWDB/(ABS(FMCPR4D2))
            FMCPR4D2=FMCPR4D2/(ABS(FMCPR4D2))
         ELSE
            FMCPR4D2=0.D0
            DLWDB=0.D0
         ENDIF
         CALL CLASSPATHS(FMCPR4D2,DLWDB,1.D0,RHOIJ,
     &            I_V,IOCLS)
         CALL GETTREENESS(IOCLS,ITREE,FMCPR4D2,I_V)
          IF(TLOG) THEN
            CALL WRITE_PATH(10,IPATH,NEL,I_V)
            IF(TLOG2) CALL WRITE_RHOMATRIX(10,RHOIJ,NEL,I_V)
            IF(TLOG3) CALL WRITE_XMATRIX(10,XIJ,I_V)
         ENDIF   
C         WRITE(6,*) "END"
C         CALL MEMORY_CHECK()
         RETURN
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
