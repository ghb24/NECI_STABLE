      SUBROUTINE WIRHODIAG(NI,BETA,I_P,NEL,
     &      NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &      ECORE,WLRI,WLSI)
         use constants, only: dp
         use global_utilities
         use legacy_data, only: irat
         use sym_mod, only: lchksymd
         use SystemData, only: BasisFn
         use util_mod, only: NECI_ICOPY
         use MemoryManager, only: TagIntType
         IMPLICIT NONE
         real(dp) BETA,FCK(*),ALAT(3),UMAT(*),ECORE
         INTEGER NEL,NBASIS,nBasisMax(5,*)
         INTEGER I_P,BRR(NBASIS),NTAY,NMSH,NMAX
         type(BasisFn) G1(*)
         INTEGER NI(NEL),K
         LOGICAL LACC
         INTEGER NLIST,NLMAX
         real(dp) RH,RHII
         real(dp) HDMAX,RHRATMAX
         INTEGER ILFREE,I,J
         real(dp), ALLOCATABLE :: RIJMAT(:),WLIST(:),RJJLIST(:),WORK(:)
         INTEGER, ALLOCATABLE :: LSTE(:,:)
         INTEGER err
         INTEGER WORKL
         INTEGER*4 INFO
         real(dp) WLRI,WLSI,SI
         INTEGER(TagIntType) ,SAVE :: tagRIJMAT=0,tagWLIST=0
         INTEGER(TagIntType),SAVE :: tagLSTE=0,tagRJJLIST=0,tagWORK=0 
         type(timer), save :: proc_timer
         character(*),parameter :: thisroutine='WIRHODIAG'
         proc_timer%timer_name=thisroutine
         call set_timer(proc_timer)
         HDMAX=5000.D0
         RHRATMAX=EXP(-BETA*HDMAX/I_P)

C.. Get how much mem we need
         NLIST=0
         CALL GENSYMDETS(NI,NEL,G1,BRR,NBASIS,LSTE,NLIST,NBASISMAX)
C.. account for LSTE(0)
         NLMAX=NLIST+1
         allocate(LSTE(NEL,0:NLIST),stat=err)
         CALL LogMemAlloc('LSTE',NLMAX*NEL,4,thisroutine,tagLSTE,err)
         allocate(RJJLIST(0:NLIST),stat=err)
         CALL LogMemAlloc('RJJLIST',NLMAX,8,thisroutine,tagRJJLIST,err)

         CALL NECI_ICOPY(NEL,NI,1,LSTE(1,0),1)
         CALL GENSYMDETS(NI,NEL,G1,BRR,NBASIS,LSTE(1,1),NLIST,NBASISMAX)
         
         IF(NLIST.GT.NLMAX) THEN
            WRITE(6,*) "WARNING on excitations"
         ENDIF

         CALL CALCRHO2(LSTE(1,0),LSTE(1,0),BETA,I_P,NEL,
     &         G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &         RH,NTAY,0,ECORE)
         WRITE(35,"(A)",advance='no') "("
         DO K=1,NEL
            WRITE(35,"(I3,A)",advance='no') NI(K),","
         ENDDO
         WRITE(35,"(A,G25.16,I4)",advance='no')
     &                                   ")", -(I_P/BETA)*LOG(RH),NLIST
         RJJLIST(0)=RH
         WLRI=LOG(RH)
         RHII=RH
         ILFREE=1
         DO I=1,NLIST
            IF(LSTE(1,I).NE.0) THEN
               LACC=LCHKSYMD(LSTE(1,0),LSTE(1,I),
     &               NEL,G1,NBASISMAX)
               IF(LACC) THEN
                  CALL CALCRHO2(LSTE(1,I),LSTE(1,I),BETA,I_P,NEL,
     &              G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &               RH,NTAY,0,ECORE)
C                  WRITE(35,"(G25.16)",advance='no') -(I_P/BETA)*LOG(RH)
                  IF(RH/RHII.GT.RHRATMAX) THEN
                     IF(ILFREE.LT.I)
     &                 CALL NECI_ICOPY(NEL,LSTE(1,I),1,LSTE(1,ILFREE),1)
                     RJJLIST(ILFREE)=RH
                     ILFREE=ILFREE+1
C                     WRITE(35,*) 1
                  ELSE
C                     WRITE(35,*) 0
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         NLIST=ILFREE
C         WRITE(6,*) NLIST," dets in energy range."
         allocate(RIJMAT(nLIST**2),stat=err)
         CALL LogMemAlloc('RIJMAT',NLIST**2,8,thisroutine,tagRIJMAT,err)
         allocate(WLIST(NLIST),stat=err)
         CALL LogMemAlloc('WLIST',NLIST,8,thisroutine,tagWLIST,err)
         WORKL=3*NLIST
         allocate(WORK(WORKL),stat=err)
         CALL LogMemAlloc('WORK',WORKL,8,thisroutine,tagWORK,err)

C.. Now we fill the RIJ array
         DO I=0,NLIST-1
            RIJMAT(I*NLIST+I+1)=RJJLIST(I)
            DO J=I+1,NLIST-1
               CALL CALCRHO2(LSTE(1,I),LSTE(1,J),BETA,I_P,NEL,
     &            G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &            RH,NTAY,-1,ECORE)
               RIJMAT(I*NLIST+J+1)=RH
            ENDDO
         ENDDO  

C.. Diagonalize
         CALL DSYEV('V','L',NLIST,RIJMAT,NLIST,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         deallocate(WORK)
         call LogMemDealloc(thisroutine,tagWORK)
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
         SI=0.D0
C.. divide through by highest eigenvalue to stop things blowing up
         DO I=0,NLIST-1
            SI=SI+RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)
     &               *((WLIST(I+1)/WLIST(NLIST))**I_P)
C            WRITE(6,"(I5,3G25.16)") I,WLIST(I+1),RIJMAT(I*NLIST+1),WI
         ENDDO
C         DO I=0,NLIST-1
C            DO J=0,NLIST-1
C               WRITE(34,"(2I5)",advance='no') I,J
C               WRITE(34,"(A)",advance='no') "  ("
C               DO K=1,NEL
C                  WRITE(34,"(I3,A)",advance='no') LSTE(K,J),","
C               ENDDO
C            WRITE(34,"(A,G25.16") ")",RIJMAT(I*NLIST+J+1)
C            ENDDO
C         ENDDO
         WLSI=LOG(SI)+I_P*LOG(WLIST(NLIST)/RHII)
         WRITE(35,"(I5,3G25.16)") NLIST+1,WLRI,WLSI,EXP(I_P*WLRI+WLSI)
         deallocate(WLIST,RIJMAT,LSTE,RJJLIST)
         call LogMemDealloc(thisroutine,tagWLIST)
         call LogMemDealloc(thisroutine,tagRIJMAT)
         call LogMemDealloc(thisroutine,tagLSTE)
         call LogMemDealloc(thisroutine,tagRJJLIST)
         call halt_timer(proc_timer)
         RETURN
      END
 
         SUBROUTINE WIRD_SUBSET(NI,BETA,I_P,NEL,
     &      NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,
     &      ECORE)
         use global_utilities
         use legacy_data, only: irat
         use util_mod, only: NECI_ICOPY
         use MemoryManager, only: TagIntType
         use SystemData, only: BasisFN
         use constants, only: dp
         IMPLICIT NONE
         real(dp) BETA,FCK(*),ALAT(3),UMAT(*),ECORE
         INTEGER NEL,NBASIS,nBasisMax(5,*)
         INTEGER I_P,BRR(NBASIS),NTAY,NMSH,NMAX
         INTEGER K,ID,IDI,IDJ
         TYPE(BasisFN) G1(*)
         INTEGER NI(NEL)
         LOGICAL LACC
         INTEGER NLIST,NLMAX
         real(dp) RH,RHII
         real(dp) HDMAX,RHRATMAX
         INTEGER ILFREE,I,J
         real(dp), ALLOCATABLE :: RIJMAT(:),WLIST(:),RJJLIST(:)
         real(dp), ALLOCATABLE :: VARRAY(:)
         INTEGER, ALLOCATABLE :: LSTE(:,:),LSTP(:)
         real(dp), allocatable :: WORK(:)
         INTEGER(TagIntType), SAVE :: tagRIJMAT,tagWLIST,tagRJJLIST
         INTEGER(TagIntType), SAVE :: tagVARRAY
         INTEGER(TagIntType), SAVE :: tagLSTE,tagLSTP,tagWORK
         type(timer), save :: proc_timer
         INTEGER err
         INTEGER WORKL
         INTEGER*4 INFO
         real(dp) WLRI,WLSI,SI,SN,SNO,SNC,SNN
         INTEGER IPLEN,NCOUNT
         INTEGER BINOMIAL
         LOGICAL LSUB
         character(*),parameter :: thisroutine='WIRHODIAG'
         proc_timer%timer_name=thisroutine
         call set_timer(proc_timer)
         HDMAX=5000.D0
         RHRATMAX=EXP(-BETA*HDMAX/I_P)

C.. Get how much mem we need
         NLIST=0
         CALL GENSYMDETS(NI,NEL,G1,BRR,NBASIS,LSTE,NLIST,NBASISMAX)

C.. account for LSTE(0)
         NLMAX=NLIST+1
         allocate(LSTE(NEL,0:NLIST),stat=err)
         CALL LogMemAlloc('LSTE',NLMAX*NEL,4,thisroutine,tagLSTE,err)
         allocate(RJJLIST(0:NLIST),stat=err)
         CALL LogMemAlloc('RJJLIST',NLMAX,8,thisroutine,tagRJJLIST,err)
         allocate(RIJMAT(nLMAX**2),stat=err)
         CALL LogMemAlloc('RIJMAT',NLMAX**2,8,thisroutine,tagRIJMAT,err)
         allocate(WLIST(NLIST),stat=err)
         CALL LogMemAlloc('WLIST',NLIST,8,thisroutine,tagWLIST,err)
         WORKL=3*NLMAX
         allocate(WORK(WORKL),stat=err)
         CALL LogMemAlloc('WORK',WORKL,4,thisroutine,tagWORK,err)
         allocate(LSTP(0:nList),stat=err)
         CALL LogMemAlloc('LSTP',NLMAX,4,thisroutine,tagLSTP,err)
         allocate(VARRAY(0:nList),stat=err)
         CALL LogMemAlloc('VARRAY',NLMAX,8,thisroutine,tagVARRAY,err)

         CALL NECI_ICOPY(NEL,NI,1,LSTE(1,0),1)
         CALL GENSYMDETS(NI,NEL,G1,BRR,NBASIS,LSTE(1,1),NLIST,NBASISMAX)

         CALL CALCRHO2(LSTE(1,0),LSTE(1,0),BETA,I_P,NEL,
     &         G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &         RH,NTAY,0,ECORE)
         LSTP(0)=0

         LSUB=.FALSE.
C         LSUB=.TRUE.
C         OPEN(36,FILE="/dev/null",STATUS="UNKNOWN")
         IF(LSUB) THEN
            WRITE(36,"(A)",advance='no') "("
            DO K=1,NEL
               WRITE(36,"(I3,A)",advance='no') NI(K),","
            ENDDO
            WRITE(36,"(A,G25.16,I5)") ")", -(I_P/BETA)*LOG(RH),NLIST
         ENDIF
         WRITE(37,"(A)",advance='no') "("
         DO K=1,NEL
            WRITE(37,"(I3,A)",advance='no') NI(K),","
         ENDDO
         WRITE(37,"(A,I5)") ")",NLIST
         RJJLIST(0)=RH
         WLRI=LOG(RH)
         RHII=RH
C.. List all the dets
         IF(LSUB) THEN
            DO I=1,NLIST
               WRITE(36,"(A)",advance='no') "("
               DO K=1,NEL
                  WRITE(36,"(I3,A)",advance='no') LSTE(K,I),","
               ENDDO
               WRITE(36,"(A)",advance='no') ") | "
            ENDDO
            WRITE(36,*)
         ENDIF
         SN=0.D0
         SNC=0.D0
         DO I=1,NLIST
            CALL CALCRHO2(LSTE(1,I),LSTE(1,I),BETA,I_P,NEL,
     &              G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &               RH,NTAY,0,ECORE)
            IF(LSUB)
     &         WRITE(36,"(I5,G25.16,A)",advance='no')
     &                                            I,-(I_P/BETA)*LOG(RH)
     &                  ," | "
         ENDDO
         IF(LSUB) WRITE(36,*)
C.. Recurse over all lengths of permutation         
         DO IPLEN=0,NLIST
C.. Initiate with the first IPLEN dets
            NCOUNT=0
            DO I=1,IPLEN
               LSTP(I)=I
            ENDDO
            NLIST=IPLEN+1
            LACC=.TRUE.
            SNO=SN
C            IF(IPLEN.EQ.0) THEN
C               SN=1.D0
C            ELSE
               SN=0.D0
C            ENDIF
            DO WHILE(LACC)
               NCOUNT=NCOUNT+1
               ILFREE=1
               DO I=1,NLIST-1
                  ID=LSTP(I)
                  IF(LSUB) THEN
                      WRITE(36,"(I5,A)",advance='no') ID,","
                      WRITE(38,"(I5,A)",advance='no') ID,","
                  ENDIF
                  CALL CALCRHO2(LSTE(1,ID),LSTE(1,ID),BETA,I_P,NEL,
     &              G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &               RH,NTAY,0,ECORE)
                  RJJLIST(ILFREE)=RH
                  ILFREE=ILFREE+1
               ENDDO
               IF(LSUB) WRITE(36,"(A)",advance='no') " | "
C.. Now we fill the RIJ array
               DO I=0,NLIST-1
                  RIJMAT(I*NLIST+I+1)=RJJLIST(I)
                  IDI=LSTP(I)
                  DO J=I+1,NLIST-1
                     IDJ=LSTP(J)
                     CALL CALCRHO2(LSTE(1,IDI),LSTE(1,IDJ),BETA,I_P,NEL,
     &            G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &            RH,NTAY,-1,ECORE)
                     RIJMAT(I*NLIST+J+1)=RH
                  ENDDO
               ENDDO  
               IF(LSUB) THEN
                  WRITE(38,*)
                  DO I=0,NLIST-1
                     DO J=0,NLIST-1
                        IF(J.LT.I) THEN
                           RH=RIJMAT(J*NLIST+I+1)
                        ELSE
                           RH=RIJMAT(I*NLIST+J+1)
                        ENDIF
                        WRITE(38,"(G20.10)",advance='no') RH
                     ENDDO
                     WRITE(38,*)
                  ENDDO
                  WRITE(38,*)
               ENDIF
C.. Diagonalize
               CALL DSYEV('V','L',NLIST,
     &            RIJMAT,NLIST,WLIST,WORK,WORKL,INFO)
               IF(INFO.NE.0) THEN
                  WRITE(6,*) 'DYSEV error: ',INFO
                  STOP
               ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
               SI=0.D0
C.. divide through by highest eigenvalue to stop things blowing up
               DO I=0,NLIST-1
                  SI=SI+RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)
     &                     *((WLIST(I+1)/WLIST(NLIST))**I_P)
C            WRITE(6,"(I5,3G25.16)") I,WLIST(I+1),RIJMAT(I*NLIST+1),WI
               ENDDO
               IF(LSUB) THEN
                  DO I=0,NLIST-1
                     WRITE(38,"(G20.10)",advance='no') WLIST(I+1)
                     DO J=0,NLIST-1
                        WRITE(38,"(G20.10)",advance='no')
     &                                              RIJMAT(I*NLIST+J+1)
                     ENDDO
                     WRITE(38,*)
                  ENDDO
                  WRITE(38,*)
                  WRITE(38,*)
               ENDIF
C         DO I=0,NLIST-1
C            DO J=0,NLIST-1
C               WRITE(34,"(2I5)",advance='no') I,J
C               WRITE(34,"(A)",advance='no') "  ("
C               DO K=1,NEL
C                  WRITE(34,"(I3,A)",advance='no') LSTE(K,J),","
C               ENDDO
C            WRITE(34,"(A,G25.16") ")",RIJMAT(I*NLIST+J+1)
C            ENDDO
C         ENDDO
               WLSI=LOG(SI)+I_P*LOG(WLIST(NLIST)/RHII)
              IF(LSUB) 
     &         WRITE(36,"(I5,4G25.16)") NLIST,WLRI,WLSI,EXP(WLSI),
     &                  EXP(I_P*WLRI+WLSI)
               SN=SN+EXP(WLSI)
              CALL GETNEXTPERM(LSTP(1),NLMAX-1,IPLEN,LACC,IPLEN,NLMAX-1)
            ENDDO
            SNN=SN
C            WRITE(37,*) SN
            DO I=0,IPLEN-1
               SN=SN-BINOMIAL(NLMAX-1-I,IPLEN-I)*VARRAY(I)
            ENDDO
            VARRAY(IPLEN)=SN
            SNC=SNC+SN
            WRITE(37,"(I4,3G25.16,I5)") IPLEN+1,SN,SNC,SNN,NCOUNT
            IF(NLMAX.GT.10) CALL neci_flush(37)
C            IF(IPLEN.EQ.NLMAX-1) WRITE(37,*) IPLEN+1,0.D0, SNN
         ENDDO
         deallocate(WLIST,RIJMAT,LSTE,RJJLIST,WORK)
         call LogMemDealloc(thisroutine,tagWORK)
         call LogMemDealloc(thisroutine,tagWLIST)
         call LogMemDealloc(thisroutine,tagRIJMAT)
         call LogMemDealloc(thisroutine,tagLSTE)
         call LogMemDealloc(thisroutine,tagRJJLIST)
         call halt_timer(proc_timer)
C         CLOSE(36)
         RETURN
      END


      RECURSIVE SUBROUTINE GETNEXTPERM(LSTP,NL,IPLEN,LACC,IPOS,IMAX)
         IMPLICIT NONE
         LOGICAL LACC
         INTEGER NL,IPOS,IMAX,IPLEN,LSTP(IPLEN)
         IF(IPOS.LT.1.OR.IPLEN.LT.1) THEN
            LACC=.FALSE.
            RETURN
         ENDIF
         IF(LSTP(IPOS).EQ.IMAX) THEN
            IF(IPOS.GT.1) THEN
               CALL GETNEXTPERM(LSTP,NL,IPLEN,LACC,IPOS-1,IMAX-1)
               LSTP(IPOS)=LSTP(IPOS-1)+1
            ELSE
               LACC=.FALSE.
            ENDIF
         ELSE
            LSTP(IPOS)=LSTP(IPOS)+1
         ENDIF
         RETURN
      END 

      INTEGER FUNCTION BINOMIAL(N,R)
         INTEGER N,R,I
         IF(R.GT.N.OR.N.LT.0.OR.R.LT.0) THEN
            BINOMIAL=0
         ELSE
            BINOMIAL=1
            DO I=1,R
               BINOMIAL=(BINOMIAL*(N-I+1))/I
            ENDDO
         ENDIF
         RETURN
      END

      FUNCTION RHODIAG_CP(RHOIJ,I_P,I_V)
         use global_utilities
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I_P,I_V
         type(timer), save :: proc_timer
         real(dp) RHOIJ(0:I_V,0:I_V)
         real(dp) RIJMAT(I_V,I_V),RHODIAG_CP
         real(dp) WLIST(I_V),WORK(3*I_V)         
         INTEGER*4 INFO
         INTEGER I,J
         real(dp) SI
         real(dp) RII
         RII=RHOIJ(0,0)
C.. Diagonalize
         proc_timer%timer_name='RHODIAG_CP'
         call set_timer(proc_timer)
         RIJMAT(1:I_V,1:I_V)=0.d0
         DO I=1,I_V
            DO J=I,I_V
               RIJMAT(I,J)=RHOIJ(I-1,J-1)+0.D0
               RIJMAT(I,J)=RIJMAT(I,J)+0.D0
               SI=RHOIJ(I-1,J-1)
            ENDDO
         ENDDO
C         RIJMAT(1,1)=1.D0
C         RIJMAT(1,2)=1.D0
C         RIJMAT(2,1)=1.D0
C         RIJMAT(2,2)=1.D0
C         WRITE(6,*) ((RIJMAT(I,J),J=1,I_V),I=1,I_V)
         CALL DSYEV('V','U',I_V,RIJMAT(1,1),I_V,
     &            WLIST(1),WORK(1),3*I_V,INFO)
C         WRITE(6,*) ((RIJMAT(I,J),J=1,I_V),I=1,I_V)
C         WRITE(6,*) (WLIST(I),I=1,I_V)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
         SI=0.D0
C.. divide through by RHOII^P
         DO I=1,I_V
            SI=SI+RIJMAT(1,I)*RIJMAT(1,I)
     &               *((WLIST(I)/RII)**I_P)
C           WRITE(6,"(I5,3G25.16)") I,WLIST(I+1),RIJMAT(I*NLIST+1),WI
         ENDDO
         RHODIAG_CP=SI
C         DO I=0,NLIST-1
C            DO J=0,NLIST-1
C               WRITE(34,"(2I5)",advance='no') I,J
C               WRITE(34,"(A)",advance='no') "  ("
C               DO K=1,NEL
C                  WRITE(34,"(I3,A)",advance='no') LSTE(K,J),","
C               ENDDO
C            WRITE(34,"(A,G25.16") ")",RIJMAT(I*NLIST+J+1)
C            ENDDO
C         ENDDO
         
         call halt_timer(proc_timer)
         RETURN
      END


C.. 29/6/06 Based on RHODIAG_CPP, instead of diagonalizing
C.. a matrix of RHOIJ elements, and raising the result ^P,
C.. we diagonalize the HIJ elements and work out e^-beta lambda.

C.. See 1/7/04
      RECURSIVE FUNCTION HDIAG_CPP(HIJ,
     &       I_P,I_V,IMISS,TSUB,BETA,DLWDB,HIJS) result (HDIAG_CPPRES)
         use constants, only: dp
         use global_utilities
         use HElem
         IMPLICIT NONE
         INTEGER I_P,I_V
         type(timer), save :: proc_timer
         HElement_t HIJ(I_V+1,I_V+1),RIJMAT(I_V,I_V)
         real(dp) WLIST(I_V),WORK(3*I_V)         
         HElement_t NWORK(4*I_V)
         INTEGER*4 INFO
         INTEGER I,J,IMISS,II,IJ
         real(dp) SI,SI2
C.. do we subtract out lower vertices here or later?
         LOGICAL TSUB
         real(dp) BETA
         HElement_t HIJS(I_V+1),HIJS2(I_V),DLWT,T,U
         real(dp) HDiag_CPPRES
         real(dp) R,S,DLWDB,DD2
         real(dp) RII
! Optimise the 1V case
         IF(I_V.EQ.1) THEN
            DLWDB=HIJ(1,1)
            HDIAG_CPPRES=1.D0
            RETURN
         ENDIF
         R=HIJ(1,1)
         S=-BETA
         R=R*S
         RII=EXP(R)
C.. Diagonalize
C         WRITE(6,*) "...",I_V,IMISS
         proc_timer%timer_name='HDIAG_CPP '
         call set_timer(proc_timer,55)
         RIJMAT(1:I_V,1:I_V)=(0.d0)
         DLWDB=0.D0
         II=0
         DO I=1,I_V+1
            IF(I.NE.IMISS) THEN
               IJ=II
               II=II+1
               DO J=I,I_V+1
                  IF(J.NE.IMISS) THEN
                     IJ=IJ+1
                     RIJMAT(II,IJ)=HIJ(I,J)
                  ENDIF
               ENDDO
               HIJS2(II)=HIJS(I)
            ENDIF
         ENDDO
         SI=0.D0
C.. Now subtract out the smaller submatrices first
C.. In order to count the subsets only once, we need to only
C.. remove up to IMISS 
         IF(TSUB) THEN
            DO I=2,IMISS-1
CIMISS-1
               DD2=0.D0
               SI=SI-HDIAG_CPP(RIJMAT,I_P,I_V-1,I,TSUB,
     &            BETA,DD2,HIJS2)
               DLWDB=DLWDB-DD2
            ENDDO
         ENDIF
         IF(HElement_t_size.EQ.1) THEN
            CALL DSYEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,WORK,3*I_V,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. now calculate exp(-beta lambda) for each eigenvalue, with the 
C.. appropriate projection onto the root
            SI2=0.D0
            DLWT=0.D0
            DO I=1,I_V
C            WRITE(6,*) WLIST(I),RIJMAT(1,I)
               R=HIJ(1,1)
               R=EXP(-BETA*(WLIST(I)-R))
               S=RIJMAT(1,I)*RIJMAT(1,I)
               SI2=SI2+S*R
C/RII
               T=R
C/RII
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
               DO J=1,I_V
                  U=HIJS2(J)*RIJMAT(J,I)*RIJMAT(1,I)
                  DLWT=DLWT+U*T
!                 WRITE(6,*) I,J,HIJS(J),RIJMAT(J,I),RIJMAT(1,I),U,T,DLWT
               ENDDO
            ENDDO
         ELSE
!.. The complex case
            CALL ZHEEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,NWORK,4*I_V,WORK,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'ZHEEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. now calculate exp(-beta lambda) for each eigenvalue, with the 
C.. appropriate projection onto the root
            SI2=0.D0
            DLWT=0.D0
            DO I=1,I_V
C            WRITE(6,*) WLIST(I),RIJMAT(1,I)
               S=abs(RIJMAT(1,I))**2
               R=HIJ(1,1)
               R=EXP(-BETA*(WLIST(I)-R))
               SI2=SI2+S*R
C%/RII
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
               U=R
C/RII
               DO J=1,I_V
#ifdef __CMPLX
                  T=HIJS2(J)*RIJMAT(J,I)*DCONJG(RIJMAT(1,I))
#else
                  T=HIJS2(J)*RIJMAT(J,I)*(RIJMAT(1,I))
#endif
                  DLWT=DLWT+T*U
!                 WRITE(6,*) I,J,HIJS(J),RIJMAT(J,I),RIJMAT(1,I),T,U,DLWT
               ENDDO
            ENDDO
         ENDIF
         S=DLWT
         DLWDB=DLWDB+S
         HDIAG_CPPRES=SI+SI2
         call halt_timer(proc_timer)
         RETURN
      END
C.. See 1/7/04
      RECURSIVE FUNCTION RHODIAG_CPP(RHOIJ,
     &    I_P,I_V,IMISS,TSUB,DBETA,DLWDB,HIJS,tLogWeight) 
     &      RESULT (RHODIAG_CPPRES)
         use constants, only: dp
         use global_utilities
         use HElem
         IMPLICIT NONE
         INTEGER I_P,I_V
         type(timer), save :: proc_timer
         HElement_t RHOIJ(I_V+1,I_V+1),HIJS(I_V+1),HIJS2(I_V)
         HElement_t RIJMAT(I_V,I_V),NWORK(4*I_V),S,DLWT,S2
         real(dp) WORK(3*I_V),R
         real(dp) WLIST(I_V)         
         INTEGER*4 INFO
         INTEGER I,J,IMISS,II,IJ
         real(dp) SI,SI2,SS,SS2
C.. do we subtract out lower vertices here or later?
         LOGICAL TSUB
         real(dp) DLWDB,RII,DD2
         real(dp) RhoDiag_CPPRES
         real(dp) DBETA
         LOGICAL tLogWeight
! Optimise the 1V case
         IF(I_V.EQ.1) THEN
            DLWDB=HIJS(1)
            if(tLogWeight) then
               RHODIAG_CPPRES=0.d0
            else
               RHODIAG_CPPRES=1.D0
            ENDIF
            RETURN
         ENDIF
         RII=RHOIJ(1,1)
C.. Diagonalize
C         WRITE(6,*) "...",I_V,IMISS
         proc_timer%timer_name='RHODIAG_C2'
         call set_timer(proc_timer,55)
         RIJMAT(1:I_V,1:I_V)=(0.d0)
         IF(DBETA.NE.0.D0) DLWDB=0.D0
         II=0
         DO I=1,I_V+1
            IF(I.NE.IMISS) THEN
               IJ=II
               II=II+1
               DO J=I,I_V+1
                  IF(J.NE.IMISS) THEN
                     IJ=IJ+1
                     RIJMAT(II,IJ)=RHOIJ(I,J)
                  ENDIF
               ENDDO
               HIJS2(II)=HIJS(I)
            ENDIF
         ENDDO
         SI=0.D0
C.. Now subtract out the smaller submatrices first
C.. In order to count the subsets only once, we need to only
C.. remove up to IMISS 
         IF(TSUB) THEN
            DO I=2,IMISS-1
               DD2=0.D0
               IF(tLogWeight) THEN
! we return e~' instead of w' E~'
! e~'[G]=e~[G]-sum_{g in G} e~'[g]
!  we ignore w' returned in SI
!  DD2 returns e~'[g]
!  where g is this graph missing out vertex I.
                SI=SI-RHODIAG_CPP(RIJMAT,I_P,I_V-1,I,TSUB,
     &            DBETA,DD2,HIJS2,tLogWeight)
                IF(DBETA.NE.0.D0) DLWDB=DLWDB-DD2
               ELSE
! we return w' E~'
! w'[G]E~'[G]=w[G]E~[G]-sum_{g in G} w'[g] E~'[g]
!SI is w' and DLWDB is E~'  
!  The next line subtracts out the contribution from this graph, but missing out vertex I.
!   done recursively, this makes w'[G]E~'[G]
                SI=SI-RHODIAG_CPP(RIJMAT,I_P,I_V-1,I,TSUB,
     &            DBETA,DD2,HIJS2,tLogWeight)
                IF(DBETA.NE.0.D0) DLWDB=DLWDB-DD2
               ENDIF
C            WRITE(6,*) "SI=",SI
            ENDDO
         ENDIF
C         WRITE(6,*) I_V
         IF(HElement_t_size.EQ.1) THEN
            CALL DSYEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,WORK,3*I_V,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. divide through by RHOII^P
            DLWT=0.D0
            SI2=0.D0
            DO I=1,I_V
C               WRITE(6,*) WLIST(I),RIJMAT(1,I)
               R=(RII)
               R=((WLIST(I)/R)**I_P)
               S=R
               SS=S*RIJMAT(1,I)*RIJMAT(1,I)
               SI2=SI2+SS
               IF(DBETA.NE.0.D0) THEN
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                  DO J=1,I_V
                     S2=S*RIJMAT(J,I)*RIJMAT(1,I)
                     DLWT=DLWT+S2*HIJS2(J)
                  ENDDO
               ENDIF
            ENDDO
! DLWT is the value of w[G] E~[G]
! SI2 is w[G]
!DLWDB contains the subtracted out subgraphs
           if(tLogWeight) then
! we return e~' instead of w' E~'
! e~'[G]=e~[G]-sum_{g in G} e~'[g]
            S=DLWT
            S2=SI2
            S=S/S2
            SS=S
            DLWDB=SS+DLWDB
           else
            SS=DLWT
            DLWDB=DLWDB+SS
           endif
         ELSE
!.. The complex case
            CALL ZHEEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,NWORK,4*I_V,WORK,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'ZHEEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. divide through by RHOII^P
            SI2=0.D0
            DLWT=0.D0
            DO I=1,I_V
               R=(RII)
               SS=((WLIST(I)/R)**I_P)
               SS2=abs(RIJMAT(1,I))**2
               SI2=SI2+SS*SS2
               IF(DBETA.NE.0.D0) THEN
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                  DO J=1,I_V
                     S=SS
#ifdef __CMPLX
                     S=S*RIJMAT(J,I)*DCONJG(RIJMAT(1,I))
#else
                     S=S*RIJMAT(J,I)*(RIJMAT(1,I))
#endif
                     DLWT=DLWT+S*HIJS2(J)
                  ENDDO
               ENDIF
            ENDDO
! DLWT is the value of w[G] E~[G]
! SI2 is w[G]
!DLWDB contains the subtracted out subgraphs
          if(tLogWeight) then
! we return e~' instead of w' E~'
! e~'[G]=e~[G]-sum_{g in G} e~'[g]

            S=DLWT
            S2=SI2
            S=S/S2
            SS=S
            DLWDB=SS+DLWDB
           else
            SS=DLWT
            DLWDB=DLWDB+SS
           endif
         ENDIF



         if(tLogWeight) then
!  We pass the log x[G] around
            SI2=LOG(SI2)
            RHODIAG_CPPRES=SI+SI2
!  We return zero as the weight of the graph
         ELSE
            RHODIAG_CPPRES=SI+SI2
         ENDIF

         call halt_timer(proc_timer)
         RETURN
      END
!  As rhodiag_cpp, but only deal with this vertex level, and subtract out the two-vertex star contribtion from this vertex level.
      FUNCTION RHODIAG_CPPS2VS(RHOIJ,
     &            I_P,I_V,DBETA,DLWDB,HIJS)
         use constants, only: dp
         use global_utilities
         use HElem
         IMPLICIT NONE
         INTEGER I_P,I_V
         type(timer), save :: proc_timer
         HElement_t RHOIJ(I_V+1,I_V+1),HIJS(I_V+1),HIJS2(I_V)
         HElement_t RIJMAT(I_V,I_V),NWORK(4*I_V),S,DLWT,S2
         real(dp) WORK(3*I_V),R
         real(dp) WLIST(I_V)         
         INTEGER*4 INFO
         INTEGER I,J
         real(dp) SI,SI2,SS,SS2
C.. do we subtract out lower vertices here or later?
         real(dp) DLWDB,RII
         real(dp) RhoDiag_CPPS2VS
         real(dp) DBETA
! Optimise the 1V case
         IF(I_V.EQ.1) THEN
            RII=HIJS(1)
            DLWDB=DLWDB-RII
            RHODIAG_CPPS2VS=-1.D0
            RETURN
         ENDIF
         RII=RHOIJ(1,1)
C.. Diagonalize
C         WRITE(6,*) "...",I_V,IMISS
         proc_timer%timer_name='RHODIAG_C2'
         call set_timer(proc_timer,55)
         RIJMAT(1:I_V,1:I_V)=(0.d0)
         DO J=1,I_V
            RIJMAT(1,J)=RHOIJ(1,J)
            RIJMAT(J,J)=RHOIJ(J,J)
            HIJS2(J)=HIJS(J)
         ENDDO
         SI=0.D0
!  No need to deal with small submatrices
         IF(HElement_t_size.EQ.1) THEN
            CALL DSYEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,WORK,3*I_V,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. divide through by RHOII^P
            DLWT=0.D0
            SI2=0.D0
            DO I=1,I_V
C               WRITE(6,*) WLIST(I),RIJMAT(1,I)
               R=(RII)
               R=((WLIST(I)/R)**I_P)
               S=R
               SS=S*RIJMAT(1,I)*RIJMAT(1,I)
               SI2=SI2+SS
               IF(DBETA.NE.0.D0) THEN
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                  DO J=1,I_V
                     S2=S*RIJMAT(J,I)*RIJMAT(1,I)
                     DLWT=DLWT+S2*HIJS2(J)
                  ENDDO
               ENDIF
            ENDDO
C            S=DLWDB
C            S=S+DLWT
            SS=DLWT
!We're doing a subtraction
            DLWDB=DLWDB-SS
C            WRITE(6,*) I_V,DLWDB
         ELSE
!.. The complex case
            CALL ZHEEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,NWORK,4*I_V,WORK,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'ZHEEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. divide through by RHOII^P
            SI2=0.D0
            DLWT=0.D0
            DO I=1,I_V
               R=(RII)
               SS=((WLIST(I)/R)**I_P)
               SS2=abs(RIJMAT(1,I))**2
               SI2=SI2+SS*SS2
               IF(DBETA.NE.0.D0) THEN
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                  DO J=1,I_V
                     S=SS
#ifdef __CMPLX
                     S=S*RIJMAT(J,I)*DCONJG(RIJMAT(1,I))
#else
                     S=S*RIJMAT(J,I)*(RIJMAT(1,I))
#endif
                     DLWT=DLWT+S*HIJS2(J)
                  ENDDO
               ENDIF
            ENDDO
            S=DLWDB
! In subtraction mode
            S=S-DLWT
            DLWDB=S
         ENDIF
         RHODIAG_CPPS2VS=SI-SI2
         call halt_timer(proc_timer)
         RETURN
      END
! As HDIAG_CPP, but return the subtraction out of the 2-vertex star part.
      FUNCTION HDIAG_CPPS2VS(HIJ,
     &            I_V,BETA,DLWDB,HIJS)
         use constants, only: dp
         use global_utilities
         use HElem
         IMPLICIT NONE
         INTEGER I_V
         type(timer), save :: proc_timer
         HElement_t HIJ(I_V+1,I_V+1),RIJMAT(I_V,I_V)
         real(dp) WLIST(I_V),WORK(3*I_V)         
         HElement_t NWORK(4*I_V)
         INTEGER*4 INFO
         INTEGER I,J,II
         real(dp) SI,SI2
C.. do we subtract out lower vertices here or later?
         real(dp) BETA
         HElement_t HIJS(I_V+1),HIJS2(I_V),DLWT,T,U
         real(dp) HDiag_CPPS2VS
         real(dp) R,S,DLWDB
         real(dp) RII
         call stop_all('HDIAG_CPPS2VS','Routine depreciated')
         j=1
! Optimise the 1V case
         IF(I_V.EQ.1) THEN
            RII=HIJ(1,1)
            DLWDB=DLWDB-RII
            HDIAG_CPPS2VS=-1.D0
            RETURN
         ENDIF
         R=HIJ(1,1)
         S=-BETA
         R=R*S
         RII=EXP(R)
C.. Diagonalize
C         WRITE(6,*) "...",I_V,IMISS
         proc_timer%timer_name='HDIAG_CPP '
         call set_timer(proc_timer,55)
         RIJMAT(1:I_V,1:I_V)=(0.d0)
         II=0
         DO I=1,I_V
            RIJMAT(I,J)=HIJ(I,J)
            RIJMAT(J,J)=HIJ(J,J)
            HIJS2(I)=HIJS(I)
         ENDDO
         SI=0.D0
C.. Now subtract out the smaller submatrices first
C.. In order to count the subsets only once, we need to only
C.. remove up to IMISS 
         IF(HElement_t_size.EQ.1) THEN
            CALL DSYEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,WORK,3*I_V,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. now calculate exp(-beta lambda) for each eigenvalue, with the 
C.. appropriate projection onto the root
            SI2=0.D0
            DLWT=0.D0
            DO I=1,I_V
C            WRITE(6,*) WLIST(I),RIJMAT(1,I)
               R=HIJ(1,1)
               R=EXP(-BETA*(WLIST(I)-R))
               S=RIJMAT(1,I)*RIJMAT(1,I)
               SI2=SI2+S*R
C/RII
               T=R
C/RII
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
               DO J=1,I_V
                  U=HIJS2(J)*RIJMAT(J,I)*RIJMAT(1,I)
                  DLWT=DLWT+U*T
!                 WRITE(6,*) I,J,HIJS(J),RIJMAT(J,I),RIJMAT(1,I),U,T,DLWT
               ENDDO
            ENDDO
         ELSE
!.. The complex case
            CALL ZHEEV('V','U',I_V,RIJMAT,I_V,
     &            WLIST,NWORK,4*I_V,WORK,INFO)
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'ZHEEV error: ',INFO
               STOP
            ENDIF
C.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
C.. now calculate exp(-beta lambda) for each eigenvalue, with the 
C.. appropriate projection onto the root
            SI2=0.D0
            DLWT=0.D0
            DO I=1,I_V
C            WRITE(6,*) WLIST(I),RIJMAT(1,I)
               S=abs(RIJMAT(1,I))**2
               R=HIJ(1,1)
               R=EXP(-BETA*(WLIST(I)-R))
               SI2=SI2+S*R
C%/RII
C.. calculate <D|H exp(-b H)|D>/RHO_ii^P
               U=R
C/RII
               DO J=1,I_V
#ifdef __CMPLX
                  T=HIJS2(J)*RIJMAT(J,I)*DCONJG(RIJMAT(1,I))
#else
                  T=HIJS2(J)*RIJMAT(J,I)*(RIJMAT(1,I))
#endif
                  DLWT=DLWT+T*U
!                 WRITE(6,*) I,J,HIJS(J),RIJMAT(J,I),RIJMAT(1,I),T,U,DLWT
               ENDDO
            ENDDO
         ENDIF
         S=DLWT
!Subtraction mode
         DLWDB=DLWDB-S
         HDIAG_CPPS2VS=SI-SI2
         call halt_timer(proc_timer)
         RETURN
      END
