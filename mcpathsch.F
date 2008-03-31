
C.. Calculate RHO^(P)_II without having a stored H matrix
C.. SAMPLE over distinct nodes, e.g. IJKLI, with paths up to I_HMAX
C.. generated from these, and summed (e.g IJILKJI), up to max H
C.. In theory more efficient because RHO_IJ,RHO_JK, etc are calculated
C.. just once for all these paths.
C.. I_VMAX is the max number of distinct vertices in a path.
C.. I_HMAX is the max number of hops in a path.
C.. NWHTAY contains the number of samples to take for each level

C.. This approach creates a cluster of increasing size by adding to it
C.. the best vertex at each level.
      SUBROUTINE  MCPATHSR6(NI,BETA,I_P,I_HMAX,I_VMAX,NEL,NBASISMAX,G1,
     &               NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,RHOEPS,
     &               LSTE,ICE,RIJLIST,NWHTAY,ILOGGING,ECORE,ILMAX,
     &               WLRI,WLSI)
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION IPATH(NEL,0:I_VMAX)
         INTEGER NI(NEL)
         REAL*8 F(2:I_VMAX),FSQ(2:I_VMAX)
         CHARACTER*20 STR
         REAL*8 TOTAL
         REAL*8 RHOII(0:I_VMAX),RHOIJ(0:I_VMAX,0:I_VMAX)
         REAL*8 FLNSUM,FLNRII
         REAL*8 ALAT(3),RHOEPS
         INTEGER  G1(*),NBASISMAX(5,2)
         INTEGER LSTE(*),ILMAX
CNEL,0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER ICE(*)
         REAL*8 RIJLIST(*)
C0:NBASIS*NBASIS*NEL*NEL,0:I_VMAX-1)
         INTEGER NLIST(0:I_VMAX-1),LSTP(0:I_VMAX-1),BTABLE(0:I_VMAX)
         LOGICAL TLOG,TSYM,LISNAN
         INTEGER ICOUNT,ISEED,I_VMOVE
         INTEGER IADJ(0:I_VMAX-1,0:I_VMAX-1),ADJ
         REAL*8 FNCR
         REAL*8 OVALUE,WEIGHT,NVALUE
         LOGICAL LCANMOVE,BR
         REAL*8 RNORM,ND,STD
         CALL TISET('MCPATHSR6 ',ISUB)
         TLOG=BTEST(ILOGGING,1)
         IF(TLOG) THEN
            OPEN(11,FILE="MCPATHS",STATUS="OLD")
C.. go to end of file
            Call ENDSeek(11)
!            I=FSEEK(11,0,2)
            WRITE(11,"(A,$)") "("
            DO K=1,NEL
               WRITE(11,"(I3,A,$)") NI(K),","
            ENDDO
            WRITE(11,"(A)") ")"
         ENDIF
         IF(BTEST(ILOGGING,2))
     &      OPEN(10,FILE="PATHS",STATUS="UNKNOWN")
C.. Set the first node to I_I
         CALL ICOPY(NEL,NI,1,IPATH(1:NEL,0),1)
         CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,
     &            BRR,NMSH,FCK,NMAX,ALAT,UMAT,RH,NTAY,0,ECORE)
         RHOII(0)=RH
         RHOIJ(0,0)=RH
         IADJ(0,0)=1
         WLRI=LOG(RHOII(0))
         TOTAL=1.D0
         IF(TLOG) WRITE(11,"(I12,2G25.16,F19.7,2I12)") 
     &         1,TOTAL,TOTAL,0.D0,1,1
c         WRITE(6,*) 0,TOTAL,TOTAL,0
C.. I_V is the number of vertices in the path
         I_V=2
         BR=.TRUE.
         DO WHILE (BR)
            WRITE(STR,"(A,I5)") "FMCPR",I_V
            CALL TISET(STR,ISUB2)
            L=0
            LT=0
            BTABLE(0)=0
            ISEED=7
C.. I_VMOVE=0 indicates we need to generate a valid path.
            CALL WIDENRHOMAT(RHOIJ,I_V,I_VMAX)
            NVALUE=
     &            FMCPR6(NI,BETA,I_P,IPATH,I_V,NEL,NBASISMAX,
     &              G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &              RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &             ECORE,ICOUNT)
            F(I_V)=NVALUE
            CALL TIHALT(STR,ISUB2)
c            WRITE(6,*) I_V,F(I_V),TOTAL,TIGET(STR,ISUB2),L,LT
            IF(TLOG)
     &         WRITE(11,"(I12,2G25.16,F19.7,I12)")
     &            I_V,F(I_V),TOTAL,TIGET(STR,ISUB2),ICOUNT
            IF(LISNAN(F(I_V))) THEN
C.. save all log files
               ITIME=TIMEC()
               CALL FLUSH(11)
               CALL LOGNAN(NI,NEL,BETA,ITIME)
               WRITE(6,*) "WARNING: nan found at time",ITIME
               WRITE(6,"(A,$)") "  nan det=("
               DO K=1,NEL
                  WRITE(6,"(I3,A,$)") NI(K),","
               ENDDO
               WRITE(6,"(A,$)") "),"
            ENDIF
            IF(ICOUNT.GT.0) I_V=I_V+1
            IF(I_V.GT.I_VMAX.OR.ICOUNT.EQ.0) BR=.FALSE.
         ENDDO
         IF(TLOG) CLOSE(11)
         IF(BTEST(ILOGGING,2)) CLOSE(10)
         WLSI=LOG(F(I_V-1))
         CALL TIHALT('MCPATHSR6 ',ISUB)
         RETURN
      END

      REAL*8 FUNCTION FMCPR6(NI,BETA,I_P,IPATH,I_V,NEL,
     &   NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &   RHOEPS,RHOII,RHOIJ,NWHTAY,I_HMAX,ILOGGING,TSYM,
     &   ECORE,ICOUNT)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),I_P,IPATH(NEL,0:I_V),I_V
         INTEGER NBASISMAX(5,2),G1,NBASIS,BRR(NBASIS),NMAX
         INTEGER NTAY,NWHTAY,I_HMAX,ILOGGING,ISEED,NMSH
         LOGICAL TSYM
         REAL*8 BETA,FCK(*),ALAT(*),UMAT(*),ECORE
         REAL*8 RHOEPS,RHOII(0:I_V),RHOIJ(0:I_V,0:I_V)
         INTEGER I,J,K,NELEC,NADJ,ICE,NBEST
         REAL*8 RH,CALCPATHS,VALUE,NVALUE,RAT
         LOGICAL LCANMOVE,LISINPATH
         LOGICAL TLOG,TLOG2,LDONE
         INTEGER INDJ(NEL),NJ(NEL),KI(5),KJ(5)
         INTEGER IGETEXCITLEVEL,ICOUNT
         TLOG=BTEST(ILOGGING,2)
         TLOG2=BTEST(ILOGGING,3)
C.. First decide what MC move - either change the moving vertex or
C.. move the moving vertex.
C.. enumerate all excitations of this graph working out which is best
         LDONE=.FALSE.
         INDJ(1)=0
         NELEC=1
         CALL ICOPY(NEL,IPATH(1,0),1,IPATH(1,I_V),1)
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KI)
         ICOUNT=0
         DO I=1,4
            KJ(I)=0
         ENDDO
         DO WHILE (.NOT.LDONE)
            CALL GENNEXTSYMDET(KI,NEL,G1,BRR,NBASIS,
     &               LDONE,NJ,IPATH(1,I_V-1),INDJ,KJ,NELEC,NBASISMAX)
C.. Now insert the node in the path at I_VMOVE.  The path must then be
C.. lexically reordered and I_VMOVE,IADJ and the RHOIJ matrix 
C.. modified appropriately
C.. The relevant elements of IADJ and RHOIJ are copied to the new path's
            IF(.NOT.LDONE
     &      .AND..NOT.LISINPATH(IPATH(1,I_V-1),IPATH,NEL,I_V-1,-1)) THEN
               NADJ=0
               DO J=0,I_V-1
                  ICE=IGETEXCITLEVEL(IPATH(1,I_V-1),IPATH(1,J),NEL)
                  CALL CALCRHO2(IPATH(1,I_V-1),IPATH(1,J),BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,
     &             RH,NTAY,ICE,ECORE)
                  IF(ICE.LE.2) THEN
                     ICE=1
                  ELSE
                     ICE=0
                  ENDIF
                  IF(ABS(RH).LE.RHOEPS) THEN
                     RH=0.D0
                     ICE=0
                  ENDIF
                  NADJ=NADJ+ICE
                  RHOIJ(I_V-1,J)=RH
                  RHOIJ(J,I_V-1)=RH
                  IF(J.EQ.I_V-1) RHOII(I_V-1)=RH
               ENDDO
               NADJ=2
C.. IPATH now contains the path, and RHOII and RHOIJ the appropriate
C.. matrix elements.  We now call CALCPATHS to get the path weight
C.. Now calculate a new value if required
               IF(NADJ.GT.1) THEN
                  ICOUNT=ICOUNT+1
                  CALL ICOPY(NEL,IPATH(1,0),1,IPATH(1,I_V),1)
                  RHOII(I_V)=RHOII(0)
                  IF(TLOG) THEN
                     CALL WRITEPATH(10,IPATH,I_V,NEL,.FALSE.)   
                     IF(TLOG2) CALL WRITERHOMAT(10,RHOIJ,I_V,NEL,.TRUE.)
                  ENDIF   
                  VALUE=CALCPATHS(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &               I_P,NWHTAY,NEL,I_V,ILOGGING)
                  IF(VALUE.GT.FMCPR6) THEN
                     CALL ICOPY(NEL,IPATH(1,I_V-1),1,NBEST,1)
                     FMCPR6=VALUE
                  ENDIF
                  IF(TLOG) WRITE(10,*) VALUE 
               ENDIF
            ENDIF
         ENDDO
C.. Restore the best one
         IF(ICOUNT.GT.0) THEN
            CALL ICOPY(NEL,NBEST,1,IPATH(1,I_V-1),1)         
            DO J=0,I_V-1
               ICE=IGETEXCITLEVEL(IPATH(1,I_V-1),IPATH(1,J),NEL)
               CALL CALCRHO2(IPATH(1,I_V-1),IPATH(1,J),BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,
     &           RH,NTAY,ICE,ECORE)
               IF(ICE.LE.2) THEN
                  ICE=1
               ELSE
                  ICE=0
               ENDIF
               IF(ABS(RH).LE.RHOEPS) THEN
                  RH=0.D0
                  ICE=0
               ENDIF
               NADJ=NADJ+ICE
               RHOIJ(I_V-1,J)=RH
               RHOIJ(J,I_V-1)=RH
               IF(J.EQ.I_V-1) RHOII(I_V-1)=RH
            ENDDO
         ELSE
            FMCPR6=0.D0
         ENDIF
         RETURN
      END


C.. move elements in RHOIJ from a 0:I_V-1 **2 to a 0:I_V **2 matrix
      SUBROUTINE WIDENRHOMAT(RHOIJ,I_V,I_VMAX)
         IMPLICIT NONE
         REAL*8 RHOIJ(0:(I_VMAX+1)**2-1)
         INTEGER I_VMAX,I_V,I,J
         DO J=I_V-1,0,-1
            DO I=I_V-2,0,-1
               RHOIJ(J*(I_V+1)+I)=RHOIJ(J*I_V+I)
            ENDDO
            RHOIJ(J*(I_V+1)+I_V)=0.D0
            RHOIJ(J*(I_V+1)+I_V-1)=0.D0
         ENDDO
         RETURN
      END
