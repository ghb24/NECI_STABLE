C.. Ti=Sum_a,b(c_ia* c_ib <u_a(1)|h(1)|u_b(1)>)
C..   =Sum_a(|c_ia|^2 <u_a|(p^2)/2|u_a>)
      SUBROUTINE CALCHFTMAT(TMAT,TMAT2,NBASIS,HFBASIS,NORBUSED)
         IMPLICIT NONE
         INTEGER NBASIS
         REAL*8 TMAT(NBASIS,NBASIS),HFBASIS(NBASIS,NBASIS)
         REAL*8 TMAT2(NBASIS,NBASIS)
         INTEGER I,J,A,B,NORBUSED
         REAL*8 SUM,S1
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         OPEN(10,FILE='TMAT2',STATUS='UNKNOWN')
         DO I=1,NORBUSED
            DO J=1,NORBUSED
               SUM=0.D0
               DO A=1,NBASIS
                  DO B=1,NBASIS
                    SUM=SUM+HFBASIS(I,A)*HFBASIS(J,B)*TMAT(A,B)
                  ENDDO
               ENDDO
               TMAT2(I,J)=SUM
               IF(ABS(SUM).GT.1.D-10) WRITE(10,*) I,J,TMAT2(I,J)
            ENDDO
         ENDDO
         CLOSE(10)
         RETURN
      END
C.. Ti=Sum_a,b(c_ia* c_ib <u_a(1)|h(1)|u_b(1)>)
C..   =Sum_a(|c_ia|^2 <u_a|(p^2)/2|u_a>)
      SUBROUTINE READHFTMAT(TMAT,TMAT2,NBASIS,HFBASIS)
         IMPLICIT NONE
         INTEGER NBASIS
         REAL*8 TMAT(NBASIS,NBASIS),HFBASIS(NBASIS,NBASIS)
         REAL*8 TMAT2(NBASIS,NBASIS),R
         INTEGER I,J,A,B
         REAL*8 SUM,S1
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         OPEN(10,FILE='TMAT2',STATUS='OLD')
         I=0
         DO WHILE(.NOT.(I.EQ.NBASIS.AND.J.EQ.NBASIS))
            READ(10,*,END=11) I,J,R
            TMAT2(I,J)=R
            TMAT2(J,I)=R
         ENDDO
11       CLOSE(10)
         RETURN
      END
      SUBROUTINE READHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,ISS,
     &            HFE,ARR,BRR)
         IMPLICIT NONE
         INTEGER NBASIS,ISS,NBASISMAX(5,2),G1(*)
         REAL*8 UMAT(*)
         REAL*8 UMAT2(*),ARR(NBASIS)
         REAL*8 HFBASIS(NBASIS,NBASIS),SUM,OS,HFE(NBASIS)
         INTEGER I,J,K,L,A,B,C,D,ISUB
         INTEGER ID1,ID2,ID3,ID4,NHG,ISSM,BRR(NBASIS)
         LOGICAL LSPN
         INTEGER UMatInd
         NHG=NBASIS
         WRITE(6,*) 'READING HF UMAT'
C ==--------------------------------------------------------------------=
         OPEN(10,FILE='UMAT2',STATUS='OLD')
         DO WHILE (.TRUE.)
            
            READ(10,*,END=12) A,B,C,D,SUM
              UMAT2(UMatInd(A,B,C,D))=SUM
C..Symmetries not needed by UMatInd
!              UMAT2(C,D,A,B)=SUM
!              UMAT2(B,A,D,C)=SUM
!              UMAT2(D,C,B,A)=SUM
         ENDDO
12       CLOSE(10)
      END

      SUBROUTINE SETUPHFBASIS(NBASISMAX,G1,NBASIS,HFE,ARR,BRR)
         IMPLICIT NONE
         INTEGER NBASIS,NBASISMAX(5,2)
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(nBasis)
         REAL*8 ARR(NBASIS,2)
         REAL*8 HFE(NBASIS)
         INTEGER BRR(NBASIS),ORBORDER(8)
         INTEGER I
 
C.. We now need to modify G1.  Pretend that all basis functions have
C.. a different X quantum number, and all have the same spin.
         DO I=1,NBASIS
            G1(I)%k(1)=0
            G1(I)%Ms=-2*(MOD(I,2))+1
            G1(I)%Sym%s=1
            ARR(I,1)=HFE(I)
            ARR(I,2)=HFE(I)
            BRR(I)=I
         ENDDO
C.. Now modify NBASISMAX
         NBASISMAX(4,2)=1
         NBASISMAX(4,1)=-1
         NBASISMAX(3,2)=0
         NBASISMAX(3,1)=0
         NBASISMAX(2,1)=0
         NBASISMAX(2,2)=0
         NBASISMAX(1,1)=0
C         NBASISMAX(1,2)=NBASIS-1
         NBASISMAX(5,2)=0

C.. Generic spatial symmetry
         NBASISMAX(3,3)=1
         CALL GENMOLPSYMTABLE(1,G1,NBASIS,ARR,BRR)
         CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
         CALL GENMOLPSYMREPS(1,G1,NBASIS,ARR,BRR)

         RETURN
      END
      SUBROUTINE CALCHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,ISS,
     &            HFE,ARR,BRR,NORBUSED)
         IMPLICIT NONE
         INTEGER NBASIS,G1(*),ISS,NBASISMAX(5,2)
         REAL*8 UMAT(*)
         REAL*8 UMAT2(*),ARR(NBASIS)
         REAL*8 UMATT((((NHG/ISS)*(NHG/ISS-1))/2)**2)
         REAL*8 HFBASIS(NBASIS,NBASIS),SUM,OS,HFE(NBASIS)
         INTEGER I,J,K,L,A,B,C,D,ISUB
         INTEGER ID1,ID2,ID3,ID4,NHG,ISSM,BRR(NBASIS),NORBUSED
         LOGICAL LSPN
         INTEGER UMatInd
         NHG=NBASIS
         CALL TISET('CALCHFUMAT',ISUB)
         WRITE(6,*) 'CALCULATING HF UMAT'
         STOP "HF UMAT calculation broken through UMAT reindexing."
     &      //"  Please fix"
C ==--------------------------------------------------------------------==
         OPEN(10,FILE='UMAT2',STATUS='UNKNOWN')
C.. A, B, C, D denote basis fns in the HF basis.  These basis fns
C.. alternate in spin i.e. fn 1 has alpha, fn 2 beta, fn 3 alpha etc.
C.. We need to take into account the spin when writing out the new U
C.. matrix.  For <A(1) B(2) |U| C(1) D(2)> to be non-zero, A and C must
C.. have the same spin, and B and D must have the same spin.
C.. Thus the sum A+C must be even, as must B+D
         WRITE(6,*) "Index 1..."
         DO A=1,NORBUSED
          DO J=1,NHG/ISS
           DO K=1,NHG/ISS
            DO L=1,NHG/ISS
             SUM=0.D0
             DO I=1,NHG
              CALL GTID(NBASISMAX,I,ID1)
              SUM=SUM+HFBASIS(A,I)
     &               *UMAT(UMatInd(ID1,J,K,L))
             ENDDO
             UMATT(UMatInd(J,K,L,A))=SUM
!             IF(ABS(SUM).GT.1.D-9) WRITE(6,*) J,K,L,A,SUM
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         WRITE(6,*) "Index 2..."
         DO A=1,NORBUSED
          DO B=1,NORBUSED
           DO K=1,NHG/ISS
            DO L=1,NHG/ISS
             SUM=0.D0
             DO J=1,NHG
              CALL GTID(NBASISMAX,J,ID2)
              SUM=SUM+HFBASIS(B,J)
     &               *UMATT(UMatInd(ID2,K,L,A))
             ENDDO
             UMAT2(UMatInd(K,L,A,B))=SUM
!             IF(ABS(SUM).GT.1.D-9) WRITE(6,*) K,L,A,B,SUM
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         WRITE(6,*) "Index 3..."
         DO A=1,NORBUSED
          DO B=1,NORBUSED
           DO C=1,NORBUSED
            DO L=1,NHG/ISS
             SUM=0.D0
             DO K=1,NHG
              CALL GTID(NBASISMAX,K,ID3)
              SUM=SUM+HFBASIS(C,K)
     &               *UMAT2(UMatInd(ID3,L,A,B))
             ENDDO
             UMATT(UMatInd(L,A,B,C))=SUM
!             IF(ABS(SUM).GT.1.D-9) WRITE(6,*) L,A,B,C,SUM
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         CALL AZZERO(UMAT2,NBASIS**4)
         WRITE(6,*) "Index 4..."
         DO A=1,NORBUSED
          DO B=1,NORBUSED
           DO C=1,NORBUSED
            DO D=1,NORBUSED
             LSPN=(MOD(A+C,2).EQ.0).AND.(MOD(B+D,2).EQ.0)
             IF(LSPN.AND.UMAT2(UMatInd(A,B,C,D)).EQ.0.D0) THEN
              SUM=0.D0
              DO L=1,NHG
               CALL GTID(NBASISMAX,L,ID4)
               SUM=SUM+HFBASIS(D,L)
     &                *UMATT(UMatInd(ID4,A,B,C))
              ENDDO
              UMAT2(UMatInd(A,B,C,D))=SUM
C..Symmetries not needed by UMatInd
!              UMAT2(C,D,A,B)=SUM
!              UMAT2(B,A,D,C)=SUM
!              UMAT2(D,C,B,A)=SUM
               IF(ABS(SUM).GT.1.D-10)
     &          WRITE(10,'(4I7,F19.9)') A,B,C,D,SUM
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         CLOSE(10)
         WRITE(6,*) ' !!! FINISHED CALCULATING HF UMAT !!! '
         CALL TIHALT('CALCHFUMAT',ISUB)
         RETURN
      END
      SUBROUTINE READHFBASIS(HFBASIS,HFE,G1,NBASIS)
         IMPLICIT NONE
         INTEGER NBASIS,NQNS(5),NN
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(nBasis)
         REAL*8 HFBASIS(NBASIS,NBASIS),HFE(NBASIS)
         INTEGER I,L,J,NB,NE,IG,N
         REAL*8 VAL
         INTEGER IFINDBASISFN
         WRITE(6,*) "Loading HF BASIS"
         OPEN(10,FILE='HFBASIS',STATUS='OLD')
         READ(10,*)
         READ(10,*) NB,NE
C.. NE is NEVAL, and NB is NBASIS/2
C.. NBASIS is the number of orbitals, so *2 to get # spinorbitals
C         IF(NE.NE.NEL) STOP 'NEL in HFBASIS <> NEL'

         IF(NE.NE.NB) STOP 'NEVAL <> NBASIS in HFBASIS not supported'
         IF(NB*2.NE.NBASIS) STOP 'NBASIS in HFBASIS <> NHG'
         DO I=1,NB
            DO L=-1,1,2
               READ(10,*)
               READ(10,*) HFE(I*2+(L-1)/2)
               NQNS(4)=L
               DO J=1,NB
                  READ(10,*) NN,NQNS(1),NQNS(2),NQNS(3),VAL
                  IG=IFINDBASISFN(NQNS,G1,NBASIS)
                  HFBASIS(I*2+(L-1)/2,J*2+(L-1)/2)=VAL
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
                ENDDO
            ENDDO
         ENDDO
         CLOSE(10)
         RETURN
      END



      SUBROUTINE CALCHFBASIS(NBASIS,ISS,NBASISMAX,G1,ARR,BRR,ECORE,
     &      TMAT,UMAT,HFE,HFBASIS,NHFIT,NEL,MS,HFMIX,EDELTA,CDELTA,TRHF,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         IMPLICIT NONE
         INTEGER NBASIS,G1(*),ISS,NBASISMAX(5,2)
         REAL*8 UMAT(*)
         REAL*8 TMAT(NBASIS,NBASIS),ECORE
         REAL*8 ARR(NBASIS)
         INTEGER BRR(NBASIS)
         REAL*8 HFBASIS(NBASIS,NBASIS),HFE(NBASIS)
         REAL*8 FMAT(*),OFMAT(*)
         POINTER (IP_FMAT,FMAT),(IP_OFMAT,OFMAT)
         REAL*8 DMAT(*),ODMAT(*)
         POINTER (IP_DMAT,DMAT),(IP_ODMAT,ODMAT)
         REAL*8 WORK(*)
         POINTER (IP_WORK,WORK)
         REAL*8 HFES(*)
         POINTER (IP_HFES,HFES)
         REAL*8 R1(*),R2(*)
         POINTER (IP_R1,R1),(IP_R2,R2)
         INTEGER NHFIT,NEL,NSTART(NEL)
         INTEGER I,J,K,L,ISPN,ISPS,IS1,IS2
         REAL*8 TOT,HFMIX,EDELTA,CDELTA
         INTEGER ID1,ID2,ID3,ID4,INFO,IHFIT
         INTEGER MS,IHFMETHOD
         INTEGER NSPINS,NSBASIS
         LOGICAL TRHF,TREADHF
         REAL*8 FRAND
         INTEGER HFDET(NEL)
         INTEGER ILOGGING
C.. Allocate mem etc.
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
C         IF(TRHF) STOP 'RHF not supported'
         NSPINS=1+(NBASISMAX(4,2)-NBASISMAX(4,1))/2
         NSBASIS=NBASIS/NSPINS
C         IF(TRHF) NSPINS=1
         CALL MEMORY(IP_FMAT,NSBASIS*NSBASIS*NSPINS,"FMAT")
         CALL MEMORY(IP_OFMAT,NSBASIS*NSBASIS*NSPINS,"OFMAT")
         CALL MEMORY(IP_DMAT,NSBASIS*NSBASIS*NSPINS,"DMAT")
         CALL MEMORY(IP_ODMAT,NSBASIS*NSBASIS*NSPINS,"ODMAT")
         CALL MEMORY(IP_WORK,NSBASIS*3,"WORK")
         CALL MEMORY(IP_HFES,NSBASIS*NSPINS,"HFES")
         CALL MEMORY(IP_R1,NSBASIS*NSBASIS,"R1")
         CALL MEMORY(IP_R2,NSBASIS*NSBASIS,"R2")
C.. Generate initial HFBASIS vectors as the energy ordered single
C.. particle basis fns, separated into up and down blocks
         CALL AZZERO(FMAT,NSBASIS*NSBASIS*NSPINS)
         CALL AZZERO(DMAT,NSBASIS*NSBASIS*NSPINS)
         IF(IHFMETHOD.EQ.0.OR.IHFMETHOD.EQ.-1) THEN
            CALL UHFSCF(NBASIS,ISS,NBASISMAX,G1,ARR,BRR,ECORE,
     &         TMAT,UMAT,HFE,HFBASIS,NHFIT,NEL,MS,FMAT,DMAT,ODMAT,WORK,
     &         NSPINS,NSBASIS,HFES,OFMAT,HFMIX,EDELTA,CDELTA,TRHF,R1,R2,
     &         IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         ELSE
            CALL UHFGRADDESC(NBASIS,ISS,NBASISMAX,G1,ARR,BRR,ECORE,
     &         TMAT,UMAT,HFE,HFBASIS,NHFIT,NEL,MS,NSPINS,NSBASIS,HFES,
     &      HFMIX,FMAT,OFMAT,DMAT,ODMAT,EDELTA,CDELTA,R1,R2,WORK,TRHF,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         ENDIF

         CALL FREEM(IP_FMAT)
         CALL FREEM(IP_OFMAT)
         CALL FREEM(IP_DMAT)
         CALL FREEM(IP_ODMAT)
         CALL FREEM(IP_WORK)
         CALL FREEM(IP_R1)
         CALL FREEM(IP_R2)
         CALL FREEM(IP_HFES)
         
         RETURN
      END
C.. Unrestricted HF SCF code
      SUBROUTINE UHFSCF(NBASIS,ISS,NBASISMAX,G1,ARR,BRR,ECORE,
     &      TMAT,UMAT,HFE,HFBASIS,NHFIT,NEL,MS,FMAT,DMAT,ODMAT,WORK,
     &      NSPINS,NSBASIS,HFES,OFMAT,HFMIX,EDELTA,CDELTA,TRHF,R1,R2,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS
         INTEGER NBASIS,G1(*),ISS,NBASISMAX(5,2)
         REAL*8 UMAT(NSBASIS,NSBASIS,NSBASIS,NSBASIS)
         REAL*8 TMAT(NBASIS,NBASIS),ECORE
         REAL*8 ARR(NBASIS)
         INTEGER BRR(NBASIS)
         REAL*8 HFBASIS(NBASIS,NBASIS),HFE(NBASIS)
         REAL*8 HFES(NSBASIS,NSPINS)
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 OFMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 ODMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 WORK(NBASIS*3),HFMIX
         REAL*8 R1(NSBASIS,NSBASIS),R2(NSBASIS,NSBASIS)
         INTEGER NHFIT,NEL,NSTART(NEL)
         INTEGER I,J,K,L,ISPN,ISPN2,NELS(NSPINS)
         INTEGER INDS(NBASIS)
         REAL*8 TOT,ELAST,EDELTA,ECUR,TOT2,CDELTA
         INTEGER ID1,ID2,ID3,ID4,INFO,IHFIT
         INTEGER MS,NELEX,IHFMETHOD,IRHFB
         REAL*8 F
         REAL*8 RMSD,FRAND
         LOGICAL BR
         LOGICAL TRHF,TREADHF
         INTEGER HFDET(*)
         INTEGER ILOGGING
         F=HFMIX
C         EDELTA=1.D-8
         ELAST=1.D20
         WRITE(6,*) "Performing Hartree-Fock SCF diagonalisation..."
         IF(IHFMETHOD.EQ.-1) THEN
            WRITE(6,*) "Method -1:Rotational Mixing "
         ELSEIF(IHFMETHOD.EQ.0) THEN
            WRITE(6,*) "Method 0:Linear Mixing"
         ENDIF
         IF(NSPINS.EQ.2) THEN
            NELS(2)=(MS+NEL)/2
            NELS(1)=NEL-NELS(2)
            WRITE(6,*) " Beta, Alpha: ",NELS(1),NELS(2)
         ELSE
            NELS(1)=NEL
         ENDIF
         IF(TREADHF) THEN
            CALL  READHFFMAT(NBASIS,FMAT,HFES,G1,NSPINS,NSBASIS,.FALSE.)
         ELSE
          CALL GENHFGUESS(FMAT,NSPINS,NSBASIS,R1,R2,WORK,BRR,G1,.FALSE.,
     &            MS,FRAND,NELS,HFDET)
         ENDIF
      WRITE(6,*) "Iteration   Energy     MSD"
      BR=.TRUE.
      IHFIT=1
      IRHFB=0
      IF(TRHF.AND.NSPINS.GT.1) THEN
         IF(NELS(2).GT.NELS(1)) THEN
            IRHFB=2
         ELSE
            IRHFB=1
         ENDIF
      ENDIF
      DO WHILE(BR)
         IF(IRHFB.GT.0) THEN
            CALL DCOPY(NSBASIS*NSBASIS,FMAT(1,1,IRHFB),1,
     &         FMAT(1,1,3-IRHFB),1)
            CALL DCOPY(NSBASIS*NSBASIS,DMAT(1,1,IRHFB),1,
     &         DMAT(1,1,3-IRHFB),1)
         ENDIF 
         CALL DCOPY(NSBASIS*NSBASIS*NSPINS,DMAT,1,ODMAT,1)
         CALL DCOPY(NSBASIS*NSBASIS*NSPINS,FMAT,1,OFMAT,1)
C.. Construct the Density Matrix
         CALL GENDMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,.FALSE.)
C.. See how much our density matrix has changed from last time
         RMSD=0.D0
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS 
               DO J=1,NSBASIS
                  IF(.NOT.TRHF.OR.TRHF.AND.ISPN.EQ.IRHFB) 
     &               RMSD=RMSD+(DMAT(I,J,ISPN)-ODMAT(I,J,ISPN))**2
               ENDDO
            ENDDO
         ENDDO
         RMSD=SQRT(RMSD/(NSBASIS*NSBASIS*NSPINS))
C.. replace our HF orbitals in FMAT with the Fock matrix
C.. FMAT just stores the results and the HF orbitals (currently) in it
C.. are not used in calculating the F matrix
         CALL GENFMAT(FMAT,DMAT,TMAT,UMAT,NSBASIS,NSPINS,NBASISMAX,TRHF)
         CALL DIAGFMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,HFES,WORK,
     &         TMAT,ECORE,ECUR)
         IF(IRHFB.GT.0) THEN
            CALL DCOPY(NSBASIS*NSBASIS,FMAT(1,1,IRHFB),1,
     &         FMAT(1,1,3-IRHFB),1)
            CALL DCOPY(NSBASIS,HFES(1,IRHFB),1,
     &         HFES(1,3-IRHFB),1)
         ENDIF 
C.. FMAT now contains HF orbitals again, and ECUR the Fock Energy
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
         WRITE(6,*) IHFIT,ECUR,RMSD
C.. Now add back in some of our original F matrix
         IF(IHFMETHOD.EQ.-1) THEN
            CALL HFROTMIX(FMAT,OFMAT,NSPINS,NSBASIS,F,R1,R2,WORK)
         ELSE
            CALL HFLINMIX(FMAT,OFMAT,NSPINS,NSBASIS,F,R1,R2,WORK)
         ENDIF
         IHFIT=IHFIT+1
         IF(IHFIT.GT.NHFIT) THEN
            WRITE(6,*) "*** WARNING Hartree-Fock did not converge ***"
            BR=.FALSE.
         ENDIF
         IF(ABS(ECUR-ELAST).LT.EDELTA.AND.RMSD.LT.CDELTA
     &      .AND.IHFIT.GT.5) THEN
            WRITE(6,*) "*** Hartree-Fock converged in ",
     &         IHFIT," iterations."
            WRITE(6,*) "*** HF ENERGY=",ECUR
            BR=.FALSE.
         ENDIF
         ELAST=ECUR
         
      ENDDO
      IF(BTEST(ILOGGING,11))
     &   CALL WRITEHFPSIALL(NBASIS,FMAT,HFES,G1,NSPINS,NSBASIS,.FALSE.) 
C.. We write out HFMAT  
      CALL AZZERO(HFBASIS,NBASIS*NBASIS)
         DO I=1,NSBASIS
            DO ISPN=1,NSPINS
              K=(I-1)*NSPINS+ISPN
               DO J=1,NSBASIS
                  L=(J-1)*NSPINS+ISPN
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
                  HFBASIS(K,L)=FMAT(J,I,ISPN)
               ENDDO
               HFE(K)=HFES(I,ISPN)
            ENDDO
         ENDDO
         RETURN
      END
      SUBROUTINE GENDMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,LTRANS)
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS,NELS(NSPINS)
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DMAT(NSBASIS,NSBASIS,NSPINS)
         LOGICAL LTRANS

         INTEGER ISPN,I,J,K
         REAL*8 TOT
C.. Construct the Density Matrix
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS
               DO J=I,NSBASIS
                  TOT=0.D0
C.. Sum over occupied orbitals for our spin
                  DO K=1,NELS(ISPN)
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
                     IF(LTRANS) THEN
                        TOT=TOT+FMAT(K,I,ISPN)*FMAT(K,J,ISPN)
                     ELSE
                        TOT=TOT+FMAT(I,K,ISPN)*FMAT(J,K,ISPN)
                     ENDIF
                  ENDDO
                  DMAT(I,J,ISPN)=TOT
                  DMAT(J,I,ISPN)=TOT
               ENDDO
            ENDDO 
         ENDDO
         RETURN
      END 
      SUBROUTINE GENFMAT(FMAT,DMAT,TMAT,UMAT,NSBASIS,NSPINS,NBASISMAX,
     &                  TRHF)
         IMPLICIT NONE
         INTEGER NSBASIS,NSPINS,NBASISMAX(5,2)
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 TMAT(NSBASIS*NSPINS,NSBASIS*NSPINS)
         REAL*8 UMAT(*)
         INTEGER I,J,K,L,ISPN,ISPN2,ID1,ID2,ID3,ID4
         REAL*8 TOT
         LOGICAL TRHF,TTRANS
         REAL*8 RHFMUL
         INTEGER UMatInd
C         IF(TRHF) THEN
C            RHFMUL=2.D0
C         ELSE
            RHFMUL=1.D0
C         ENDIF
C.. Construct the Fock Matrix
         CALL AZZERO(FMAT,NSBASIS*NSBASIS*NSPINS)
C. <ui|F|uj>=tij+Sum_kl(Dkl <ui uk|U|uj ul>-s<ui uk|U|ul uj>) where s=1 if i,k same spin
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS
               DO J=I,NSBASIS
                  TOT=TMAT((I-1)*NSPINS+ISPN,(J-1)*NSPINS+ISPN)
C.. Now sum in the alpha and beta block of u matrix
                  CALL GTID(NBASISMAX,(I-1)*NSPINS+ISPN,ID1)
                  CALL GTID(NBASISMAX,(J-1)*NSPINS+ISPN,ID2)
                  DO ISPN2=1,NSPINS
                     DO K=1,NSBASIS
                        DO L=1,NSBASIS
                           CALL GTID(NBASISMAX,(K-1)*NSPINS+ISPN2,ID3)
                           CALL GTID(NBASISMAX,(L-1)*NSPINS+ISPN2,ID4)
                           TOT=TOT+RHFMUL*
     &                    DMAT(K,L,ISPN2)*UMAT(UMatInd(ID1,ID3,ID2,ID4))
                           IF(ISPN2.EQ.ISPN)
     &                        TOT=TOT-
     &                   DMAT(K,L,ISPN2)*UMAT(UMatInd(ID1,ID3,ID4,ID2))
                        ENDDO
                     ENDDO
                  ENDDO
                  FMAT(I,J,ISPN)=TOT
                  FMAT(J,I,ISPN)=TOT
               ENDDO
            ENDDO
         ENDDO
         RETURN
      END
      SUBROUTINE DIAGFMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,HFES,WORK,
     &         TMAT,ECORE,ECUR)
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS,NELS(NSPINS)
         REAL*8 HFES(NSBASIS,NSPINS)
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 TMAT(NSBASIS*NSPINS,NSBASIS*NSPINS)
         REAL*8 ECUR,ECORE
         REAL*8 WORK(3*NSBASIS)

         INTEGER INFO
         REAL*8 TOT,TOT2
         INTEGER ISPN,I,J,K         

C.. First calculate the HF energy double counting contrib
         TOT=0.D0
         DO ISPN=1,NSPINS
            DO J=1,NSBASIS
               DO K=1,NSBASIS
                  TOT=TOT+DMAT(J,K,ISPN)
     &     *(FMAT(J,K,ISPN)-TMAT((J-1)*NSPINS+ISPN,(K-1)*NSPINS+ISPN))
               ENDDO
            ENDDO
         ENDDO

C.. Now diagonalize the Fock matrix 
         DO ISPN=1,NSPINS
            CALL DSYEV('V','U',NSBASIS,FMAT(1,1,ISPN),NSBASIS,
     &         HFES(1,ISPN),WORK,3*NSBASIS, INFO)
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
         ENDDO
C.. now calculate the sum of the occupied Fock orbitals
         TOT2=0.D0
         DO ISPN=1,NSPINS
            DO I=1,NELS(ISPN)
               TOT2=TOT2+HFES(I,ISPN)
            ENDDO
         ENDDO
C.. subtract out the double counting term, and add in the core energy
         ECUR=TOT2-TOT/2.D0+ECORE
         RETURN
      END
      SUBROUTINE WRITEHFPSIALL(NBASIS,FMAT,HFES,G1,NSPINS,NSBASIS,
     &         TRANSP)
         IMPLICIT NONE
         INTEGER NBASIS,NSBASIS,NSPINS,I,J,K
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(nBasis)
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS),HFES(NSBASIS,NSPINS)
         LOGICAL TRANSP
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector unless transposed
         OPEN(10,FILE='HFBASIS',STATUS='UNKNOWN')
         WRITE(10,*) 'NBASIS,NEVAL'
         WRITE(10,*) NSBASIS,NSBASIS
         DO K=1,NSBASIS
            IF(NSPINS.GT.1) THEN
               WRITE(10,*) ' BETA ELECTRON NO: ' , K
               WRITE(10,*) HFES(K,1)
               DO I=1,NSBASIS
                  IF(TRANSP) THEN
                     WRITE(10,'(I7,1X,3I7,F19.9)') 
     &               I,(G1((I-1)*NSPINS+1)%k(J),J=1,3),FMAT(K,I,1)
                  ELSE
                     WRITE(10,'(I7,1X,3I7,F19.9)') 
     &               I,(G1((I-1)*NSPINS+1)%k(J),J=1,3),FMAT(I,K,1)
                  ENDIF
               ENDDO
            ENDIF
            WRITE(10,*) 'ALPHA ELECTRON NO: ' , K
            WRITE(10,*) HFES(K,2)
            DO I=1,NSBASIS
               IF(TRANSP) THEN
                  WRITE(10,'(I7,1X,3I7,F19.9)') 
     &            I,(G1((I-1)*NSPINS+2)%k(J),J=1,3),FMAT(K,I,2)
               ELSE
                  WRITE(10,'(I7,1X,3I7,F19.9)') 
     &            I,(G1((I-1)*NSPINS+2)%k(J),J=1,3),FMAT(I,K,2)
               ENDIF
            ENDDO
         ENDDO
         CLOSE(10)
         RETURN
      END
C.. Generate initial density matrix, as well as a guess at the HF DET
      SUBROUTINE GENHFGUESS(FMAT,NSPINS,NSBASIS,R1,R2,WORK,BRR,G1,TRANS,
     &            LMS,FRAND,NELS,HFDET)
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS),PI,R
         INTEGER ISPN,I,J,ISEED
         REAL*8 RAND
         INTEGER BRR(NSBASIS*NSPINS),NELR,IREAL,IS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER LMS
C.. Working space
         REAL*8 R1(*),R2(*),WORK(*),FRAND
         LOGICAL TRANS
         INTEGER NELS(NSPINS)
         INTEGER HFDET(*),NEL
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector unless transposed
         PI=3.141592653589793D0
         ISEED=1007
         R=RAND()
         CALL AZZERO(FMAT,NSBASIS*NSBASIS*NSPINS)
         WRITE(6,*) "Generating HF Guess..."
         NEL=0
         DO IS=1,NSPINS
            IF(LMS.GE.0) THEN
               ISPN=IS
            ELSE
               ISPN=NSPINS+1-IS
            ENDIF
            IREAL=1
            WRITE(6,"(A,I2,A,$)") "Spin ",IS,":"
            
            DO I=1,NSBASIS
               DO WHILE(G1(BRR(IREAL))%Ms.NE.(-3+2*ISPN))
                  IREAL=IREAL+1
               ENDDO
               NELR=(BRR(IREAL)-1)/NSPINS+1
               IF(TRANS) THEN
                  FMAT(I,NELR,ISPN)=1
               ELSE
                  FMAT(NELR,I,ISPN)=1
               ENDIF
               IF(I.LE.NELS(IS)) THEN
                  WRITE(6,"(I4,A,$)"),BRR(IREAL),","
                  NEL=NEL+1
                  HFDET(NEL)=BRR(IREAL)
               ENDIF
               IREAL=IREAL+1
            ENDDO
            WRITE(6,*)
            DO I=1,NSBASIS
               DO J=1,NSBASIS
                  FMAT(I,J,ISPN)=FMAT(I,J,ISPN)+FRAND*RAND()
               ENDDO
            ENDDO
            CALL GRAMSCHMIDT(FMAT(1,1,ISPN),NSBASIS)
C            CALL LOWDIN_ORTH(FMAT(1,1,ISPN),NSBASIS,R1,R2,WORK)
         ENDDO
         CALL SORTI(NEL,HFDET)
         CALL WRITEDET(6,HFDET,NEL,.TRUE.)
      END
 
      SUBROUTINE GRAMSCHMIDT(MAT,LEN)
C.. MAT(IELEMENT,IVECTOR)
         IMPLICIT NONE
         INTEGER LEN
         REAL*8 MAT(LEN,LEN)
         REAL*8 NORM,DOT,SNORM
         INTEGER I,J,K
         DO I=1,LEN
C.. First dot with all lower vectors, and remove their components
            DO J=1,I-1
               DOT=0.D0
               NORM=0.D0
               DO K=1,LEN
                  DOT=DOT+MAT(K,J)*MAT(K,I)
               ENDDO
               DO K=1,LEN
                  MAT(K,I)=MAT(K,I)-MAT(K,J)*DOT
               ENDDO 
            ENDDO
            NORM=0.D0
            DO K=1,LEN
               NORM=NORM+MAT(K,I)**2
            ENDDO        
            SNORM=SQRT(NORM)    
C            WRITE(6,*) NORM
            DO K=1,LEN
C               WRITE(6,*) MAT(K,I),MAT(K,I)/SNORM
               MAT(K,I)=MAT(K,I)/SNORM
            ENDDO
         ENDDO
         RETURN
      END

      SUBROUTINE HFLINMIX(FMAT,OFMAT,NSPINS,N,FMIX,R1,R2,WORK)
         IMPLICIT NONE
         INTEGER N,NSPINS
         REAL*8 FMIX,FMAT(N,N,NSPINS),OFMAT(N,N,NSPINS)
         REAL*8 R1(*),R2(*),WORK(*)
         INTEGER I,J,ISPN
         DO ISPN=1,NSPINS
            DO I=1,N
               DO J=1,N
                  FMAT(I,J,ISPN)=FMIX*FMAT(I,J,ISPN)
     &               +(1.D0-FMIX)*OFMAT(I,J,ISPN)
               ENDDO
            ENDDO
            IF(N.LT.3) THEN
C.. For N<2, Lowdin will return in FMAT exactly what we started with in OFMAT, which is
C.. rather pointless in mixing, so we use Gram-Schmidt to mix things up a
C.. bit.
               CALL GRAMSCHMIDT(FMAT(1,1,ISPN),N)
            ELSE
               CALL LOWDIN_ORTH(FMAT(1,1,ISPN),N,R1,R2,WORK)
            ENDIF
         ENDDO
      END

      SUBROUTINE HFROTMIX(FMAT,OFMAT,NSPINS,N,FMIX,R1,R2,WORK)
         IMPLICIT NONE
         INTEGER NSPINS,N
         REAL*8 FMIX,FMAT(N,N,NSPINS)
         REAL*8 OFMAT(N,N,NSPINS)
         REAL*8 R1(N,N),R2(N,N),WORK(3*N),ALPH
         INTEGER I,J,ISPN
C.. OFMAT is the old HF orbitals
C.. FMAT is the new HF orbitals
C.. as HF orbitals make an orthoganal set, then FMAT will merely be a
C.. rotation of OFMAT around some axis
C.. If we want to only include part of the new HF orbitals, we should
C.. only rotate OFMAT by part of that amount.

C.. F' = OFMAT.  F=FMAT.  R is the rotation
C.. F=R F'.  We want F'' (the new FMAT we're to generate) to be
C.. F'' = R^a F'  (where a=FMIX is between 0 and 1)
C.. R^a = (I +(R-I))^a = I+a(R-I)+a(a-1)/2 (R-I)^2 + ... (Taylor)

C.. R = F F'T as F and F' are orthogonal matrices
         DO ISPN=1,NSPINS


C.. Work out R = R1=1.D0 * F * F'T + 0.D0*R1
            CALL DGEMM('N','T',N,N,N,1.D0,FMAT(1,1,ISPN),N,
     &            OFMAT(1,1,ISPN),N,0.D0,R1,N)
C.. now let P=R1=R-I
            DO I=1,N
               R1(I,I)=R1(I,I)-1.D0
            ENDDO
C.. R^a = I+aP(I+(a-1)/2 P(I+ (a-2)/3 P (I+...) ) ) F'
C.. Let FMAT be the accumulator
            CALL DCOPY(N*N,OFMAT(1,1,ISPN),1,FMAT(1,1,ISPN),1)
C.. Go up to 2nd order taylor
            DO J=1,0,-1
C.. Set R2=I
               CALL UNITMX(R2,N)
C.. Work out R2=(a-J)/(J+1.D0) * R1 * FMAT + 1.D0*R2.  a=FMIX
               CALL DGEMM('N','N',N,N,N,(FMIX-J)/(J+1.D0),
     &            R1,N,FMAT(1,1,ISPN),N,1.D0,R2,N)
               CALL DCOPY(N*N,R2,1,FMAT(1,1,ISPN),1)
            ENDDO
C.. Now reorthoganalise, as this is just an approximation    
            CALL LOWDIN_ORTH(FMAT(1,1,ISPN),N,R1,R2,WORK)
         ENDDO
      END


C.. Lowdin Orthoganalize
C.. for any non-singular R, let S=R RT
C.. P = S^(-1/2) R is orthogonal.
C.. MAT is NxN and is returned as an orthogal matrix
C.. R1 and R2 are NxN workspaces
      SUBROUTINE LOWDIN_ORTH(MAT,N,R1,R2,WORK)
         IMPLICIT NONE
         INTEGER N
         REAL*8 MAT(N,N),R1(N,N),R2(N,N)
         REAL*8 L(N),LL,WORK(3*N)
         INTEGER INFO,I,J
C.. R=MAT
C.. S= R1=1.D0 * R * RT + 0.D0*R1
         CALL DGEMM('N','T',N,N,N,1.D0,MAT,N,MAT,N,0.D0,R1,N)
C.. Diagonalize S=R1 into eigenvectors U=R1 and eigenvalues L
         CALL DSYEV('V','U',N,R1,N,L,WORK,N*3,INFO)
C.. eigenvector 1 is given by R1(I,1)
         IF(INFO.NE.0) THEN
            WRITE(6,*) "INFO=",INFO," on diag in LOWDIN_ORTH. Stopping"
            STOP 'Error in LOWDIN_ORTH.'
         ENDIF
C.. Calculate P = S^(-1/2) R = U L^(-1/2) UT R.  U=R1
C.. First let R2=U R. U=R1.  R=MAT
C.. R2=1.D0 * R1 * MAT + 0.D0*R1
         CALL DGEMM('T','N',N,N,N,1.D0,R1,N,MAT,N,0.D0,R2,N)
C.. Now let R2=L^(-1/2) (U R) = L^(-1/2) R2
C.. row I is multiplied by (L(I))^(-1/2)
         DO I=1,N
            LL=L(I)**(-0.5D0)
            DO J=1,N
               R2(I,J)=R2(I,J)*LL
            ENDDO
         ENDDO
C.. Now let MAT = P = U (L^(-1/2) UT R) = U R2.  U=R1
C.. MAT=1.D0 * U * R2 + 0.D0*MAT
         CALL DGEMM('N','N',N,N,N,1.D0,R1,N,R2,N,0.D0,MAT,N)
C.. and we should be done, with an orthoganal matrix in MAT
      END                       

C.. Unrestricted Hartree Fock with a gradient descent method.
      SUBROUTINE UHFGRADDESC(NBASIS,ISS,NBASISMAX,G1,ARR,BRR,ECORE,
     &      TMAT,UMAT,HFE,HFBASIS,NHFIT,NEL,MS,NSPINS,NSBASIS,HFES,
     &      HFMIX,CMAT,OCMAT,DEDCIJ,DMAT,EDELTA,CDELTA,R1,R2,WORK,TRHF,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS
         INTEGER NBASIS,G1(*),ISS,NBASISMAX(5,2)
         REAL*8 UMAT(*)
         REAL*8 TMAT(NBASIS,NBASIS),ECORE
         REAL*8 ARR(NBASIS)
         INTEGER BRR(NBASIS)
         REAL*8 HFBASIS(NBASIS,NBASIS),HFE(NBASIS)
         REAL*8 HFES(NSBASIS,NSPINS)
         REAL*8 CMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 OCMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DEDCIJ(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 WORK(NBASIS*3)
         INTEGER INORDER(100,2),ILOGGING
         REAL*8 EORDER(100,2)
C,HFMIX
         REAL*8 R1(NSBASIS,NSBASIS),R2(NSBASIS,NSBASIS)
         INTEGER NHFIT,NEL,IHFMETHOD,NELEX2
C,NSTART(NEL)
         INTEGER I,J,K,L,ISPN,NELS(NSPINS)
         REAL*8 ELAST,HFMIX,ECUR,FRAND
C         INTEGER INDS(NBASIS)
C         REAL*8 TOT,ELAST,EDELTA,ECUR,TOT2,CDELTA
C         INTEGER ID1,ID2,ID3,ID4,INFO,
         INTEGER IHFIT
         INTEGER MS
C         REAL*8 F
C         REAL*8 RMSD
         LOGICAL BR
         LOGICAL TRHF,TREADHF
         REAL*8 SCRRES,TOT,MIX,TOT2
         INTEGER NDET1(0:NEL+1),NDET2(NEL),NSPN(NSPINS),NELEX
         INTEGER NELNEW,M,ISGNCH
         INTEGER IPOSO, IPOSN,IRHFB,JSPN
         REAL*8 RMSD, EDELTA,CDELTA,EN,ECUR2
         INTEGER HFDET(*)
         
         REAL*8 GETHELEMENT2T
         IF(NSBASIS.GT.99) STOP 'ERROR - hardcoded NSBASIS limit of 100'
         ELAST=1.D20
         WRITE(6,*) "Performing Hartree-Fock Gradient Descent..."
         IF(IHFMETHOD.EQ.1) THEN
            WRITE(6,*) "Method 1:Singles replacement "
         ELSEIF(IHFMETHOD.EQ.2) THEN
            WRITE(6,*) "Method 2:Explicit differential"
         ENDIF
         IF(NSPINS.EQ.2) THEN
            NELS(2)=(MS+NEL)/2
            NELS(1)=NEL-NELS(2)
            WRITE(6,*) " Beta, Alpha: ",NELS(1),NELS(2)
         ELSE
            NELS(1)=NEL
         ENDIF
C.. Cij - i corresponds to rows and new basis functions, phi_i
C..       j corresponds to columns and old basis functions, u_j
C.. phi_i=sum_j=1,M cij u_j
         
C         CALL AZZERO(VMAT,NSPINS*NSBASIS*NSBASIS)
         IF(TREADHF) THEN
           CALL  READHFFMAT(NBASIS,CMAT,HFES,G1,NSPINS,NSBASIS,.TRUE.)
         ELSE
           CALL GENHFGUESS(CMAT,NSPINS,NSBASIS,R1,R2,WORK,BRR,G1,.TRUE.,
     &           MS,FRAND,NELS,HFDET)
         ENDIF
C         DO ISPN=1,NSPINS
C         DO I=1,NSBASIS
C            WRITE(6,*) (FMAT(I,J,ISPN),J=1,NSBASIS)
C         ENDDO
C         ENDDO

C.. Initialize our HF det in NDET1
         DO I=1,NSPINS
            NSPN(I)=0
         ENDDO
         I=1
         NDET1(0)=0
         NDET1(NEL+1)=NBASIS+1
         DO WHILE (I.LE.NEL)
            DO ISPN=1,NSPINS
               IF(NSPN(ISPN).LT.NELS(ISPN)) THEN
                  NDET1(I)=NSPN(ISPN)*NSPINS+ISPN
                  NSPN(ISPN)=NSPN(ISPN)+1
                  I=I+1
               ENDIF
            ENDDO
         ENDDO 
         CALL WRITEDET(6,NDET1(1),NEL,.TRUE.)
         BR=.TRUE.
         WRITE(6,*) "Iteration   Energy     MSD   Fock Energy"
         IHFIT=0
         IF(TRHF.AND.NSPINS.GT.1) THEN
            IF(NELS(2).GT.NELS(1)) THEN
               IRHFB=2
            ELSE
               IRHFB=1
            ENDIF
         ENDIF
         DO WHILE(BR)
            IF(IRHFB.GT.0) THEN
               CALL DCOPY(NSBASIS*NSBASIS,CMAT(1,1,IRHFB),1,
     &            CMAT(1,1,3-IRHFB),1)
            ENDIF 
            CALL DCOPY(NSBASIS*NSBASIS*NSPINS,CMAT,1,OCMAT,1)
            IHFIT=IHFIT+1
C.. First Calculate dE/dcij         
C.. dE/dcij is automatically 0 if i>N as the HF det only depends on
C.. phi_1 to phi_N.  All values of j must be iterated as each phi_i is
C.. dependent on all u_j
            ECUR=GETHELEMENT2T(NDET1(1),NDET1(1),NEL,NBASISMAX,
     &               G1,NBASIS,TMAT,UMAT,ECORE,0,CMAT,NSBASIS,NSPINS)
C.. Calculate the Gradient
            IF(IHFMETHOD.EQ.1) THEN
             CALL CALCDEDCIJ(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,TMAT,UMAT,
     &         ECORE, NBASIS,G1,NBASISMAX,NELS,NEL,ECUR)            
            ELSEIF(IHFMETHOD.EQ.2) THEN
              CALL CALCDEDCIJ2(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,TMAT,
     &         UMAT,ECORE, NBASIS,G1,NBASISMAX,NELS,NEL,ECUR,ISS)
            ENDIF
C.. DEDCIJ now comtains all elements of dE/dcij
C.. To move down the slope, we subtract a small amount of this from cij,
C.. and re-orthogonalise
C.. HFMIX is -ve
C              WRITE(6,*) ((CMAT(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
C              WRITE(6,*)
C              WRITE(6,*) ((DEDCIJ(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
C.. modify the velocity.
C              CALL DGEMA('N','N',NSBASIS,NSBASIS,1.D0,VMAT(1,1,ISPN),
C     &         NSBASIS,MIX,DEDCIJ(1,1,ISPN),NSBASIS,R1,NSBASIS)
C              CALL DCOPY(NSBASIS*NSBASIS,R1,1,VMAT(1,1,ISPN),1)
C              CALL DGEMA('N','N',NSBASIS,NSBASIS,1.D0,CMAT(1,1,ISPN),
C     &         NSBASIS,01.D0,VMAT(1,1,ISPN),NSBASIS,R1,NSBASIS)

C.. Remove the projection of the "force" already in the direction of
C.. the coefficients
            DO ISPN=1,NSPINS
              CALL AZZERO(R1,NSBASIS*NSBASIS)
              DO I=1,NELS(ISPN)
               TOT=0.D0
               DO J=1,NSBASIS
                  TOT=TOT+DEDCIJ(I,J,ISPN)*CMAT(I,J,ISPN)
               ENDDO
               TOT2=0.D0
               DO J=1,NSBASIS
                  R1(I,J)=DEDCIJ(I,J,ISPN)-TOT*CMAT(I,J,ISPN)
                  TOT2=TOT2+R1(I,J)**2
               ENDDO
               TOT2=SQRT(TOT2)
               DO J=1,NSBASIS
C                  R1(I,J)=R1(I,J)/TOT2
               ENDDO
              ENDDO
              MIX=-HFMIX
C/ABS(ECUR)
C              IF(ABS(ECUR).LT.1.D-4) MIX=HFMIX
              IF(IHFIT.GT.NHFIT) BR=.FALSE.
              CALL DGEMA('N','N',NSBASIS,NSBASIS,1.D0,CMAT(1,1,ISPN),
     &         NSBASIS,MIX,R1,NSBASIS,R2,NSBASIS)
              CALL DCOPY(NSBASIS*NSBASIS,R2,1,CMAT(1,1,ISPN),1)
C              WRITE(6,*)
C              WRITE(6,*) ((CMAT(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
C              WRITE(6,*)
C              WRITE(6,*)
              CALL LOWDIN_ORTH(CMAT(1,1,ISPN),NSBASIS,R1,R2,WORK)
            ENDDO
            RMSD=0.D0
            DO ISPN=1,NSPINS
               DO I=1,NSBASIS 
                  DO J=1,NSBASIS
                   IF(.NOT.TRHF.OR.TRHF.AND.ISPN.EQ.IRHFB) 
     &                RMSD=RMSD+(CMAT(I,J,ISPN)-OCMAT(I,J,ISPN))**2
                  ENDDO
               ENDDO
            ENDDO
            RMSD=SQRT(RMSD/(NSBASIS*NSBASIS*NSPINS))
            IF(IHFIT.GT.NHFIT) THEN
               WRITE(6,*) "** WARNING Hartree-Fock did not converge **"
               BR=.FALSE.
            ENDIF
            IF(ABS(ECUR-ELAST).LT.EDELTA.AND.RMSD.LT.CDELTA
     &         .AND.IHFIT.GT.5) THEN
               WRITE(6,*) "*** Hartree-Fock converged in ",
     &            IHFIT," iterations."
               WRITE(6,*) "*** HF ENERGY=",ECUR
               BR=.FALSE.
            ENDIF
            ELAST=ECUR

C.. Construct the Density Matrix
         CALL GENDMAT(NSPINS,NSBASIS,NELS,CMAT,DMAT,.TRUE.)
C.. Use the Density Matrix to generate the Fock matrix (in DEDCIJ)
         CALL GENFMAT(DEDCIJ,DMAT,TMAT,UMAT,NSBASIS,NSPINS,NBASISMAX,
     &      TRHF)
         CALL DIAGFMAT(NSPINS,NSBASIS,NELS,DEDCIJ,DMAT,HFES,WORK,
     &         TMAT,ECORE,ECUR2)

            WRITE(6,"(I6,$)") IHFIT
            WRITE(6,*) ECUR,RMSD,ECUR2
C.. DEDCIJ now contains HF orbitals, and ECUR the Fock Energy
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector





C.. calculate the orbital energies every time
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS
               NELEX=(I-1)*NSPINS+1+ISPN-1
               EN=0.D0
               CALL GETTRTMATEL(NELEX,NELEX,TMAT,NBASIS,CMAT,NSBASIS,
     &               NSPINS,EN)
               HFES(I,ISPN)=EN
               DO JSPN=1,NSPINS
                  DO J=1,NELS(JSPN)
                     NELEX2=(J-1)*NSPINS+1+JSPN-1
                     IF(NELEX.NE.NELEX2) THEN
C.. we're not allowed to count the current electron again
                        CALL GETTRUMATEL(NELEX,NELEX2,NELEX,NELEX2,
     &                     UMAT,NBASISMAX,NBASIS,
     &                     CMAT,NSBASIS,NSPINS,ISS,EN)
                        
C                        HFES(I,ISPN)=HFES(I,ISPN)+EN
                        EN=0.D0
                        IF(ISPN.EQ.JSPN)
     &                     CALL GETTRUMATEL(NELEX,NELEX2,NELEX2,NELEX,
     &                        UMAT,NBASISMAX,NBASIS,
     &                        CMAT,NSBASIS,NSPINS,ISS,EN)
C                        HFES(I,ISPN)=HFES(I,ISPN)-EN
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            DO I=1,NSBASIS
               INORDER(I,ISPN)=I
            ENDDO
            CALL DCOPY(NSBASIS,HFES(1,ISPN),1,EORDER(1,ISPN),1)
            CALL SORT2(NSBASIS,EORDER(1,ISPN),INORDER(1,ISPN))
         ENDDO
      IF(BR) THEN
C.. Now re-order the orbitals
C.. although we don't actually use this info
         DO ISPN=1,NSPINS
            CALL AZZERO(R1,NSBASIS*NSBASIS)
            DO I=1,NSBASIS
C.. If R1(2,1) is occupied, then postmultiplying C by R1 moves
C.. column 2 in C to column 1.
C.. INORDER(1,ISPN) is the old index of the lowest energy orb
               R1(I,INORDER(I,ISPN))=1
            ENDDO
C.. Work out R2=1.D0 CMAT*R1+ 0.D0*R2
            CALL DGEMM('N','N',NSBASIS,NSBASIS,NSBASIS,1.D0,
     &            R1,NSBASIS,CMAT(1,1,ISPN),NSBASIS,
     &            0.D0,R2,NSBASIS)
C            CALL DCOPY(NSBASIS*NSBASIS,R2,1,CMAT(1,1,ISPN),1)
            CALL DGEMM('N','N',NSBASIS,NSBASIS,NSBASIS,1.D0,
     &            R1,NSBASIS,OCMAT(1,1,ISPN),NSBASIS,
     &            0.D0,R2,NSBASIS)
C            CALL DCOPY(NSBASIS*NSBASIS,R2,1,OCMAT(1,1,ISPN),1)
         ENDDO            
 
            DO I=1,NSBASIS
               DO ISPN=1,NSPINS
C                  WRITE(6,*) I,ISPN*2-3,HFES(I,ISPN)
C                  WRITE(6,*) I,ISPN*2-3,HFES(INORDER(I,ISPN),ISPN)
C,
C     &               INORDER(I,ISPN)
               ENDDO
            ENDDO
            DO ISPN=1,NSPINS
               DO I=1,NSBASIS
C                  WRITE(6,*) I,ISPN*2-3,EORDER(I,ISPN),INORDER(I,ISPN)
CHFES(INORDER(I,ISPN),ISPN),
C     &               INORDER(I,ISPN)
               ENDDO
            ENDDO
      ENDIF
         IF(MOD(IHFIT,10).EQ.0.AND.BTEST(ILOGGING,11))
     &    CALL WRITEHFPSIALL(NBASIS,CMAT,HFES,G1,NSPINS,NSBASIS,.TRUE.) 
      ENDDO
      IF(BTEST(ILOGGING,11))
     &   CALL WRITEHFPSIALL(NBASIS,CMAT,HFES,G1,NSPINS,NSBASIS,.TRUE.) 
C.. We write out HFMAT  
      CALL AZZERO(HFBASIS,NBASIS*NBASIS)
         DO I=1,NSBASIS
            DO ISPN=1,NSPINS
              K=(I-1)*NSPINS+ISPN
               DO J=1,NSBASIS
                  L=(J-1)*NSPINS+ISPN
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
                  HFBASIS(K,L)=CMAT(I,J,ISPN)
               ENDDO
               HFE(K)=HFES(I,ISPN)
            ENDDO
         ENDDO
         RETURN
         
      END 



      SUBROUTINE CALCDEDCIJ(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,TMAT,UMAT,
     &         ECORE, NBASIS,G1,NBASISMAX,NELS,NEL,ECUR)            
         IMPLICIT NONE
         REAL*8 CMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DEDCIJ(NSBASIS,NSBASIS,NSPINS)
         INTEGER NSBASIS,NSPINS,NBASIS
         INTEGER NDET1(0:NEL+1),NDET2(NEL),NELS(NSPINS),NEL
         REAL*8 UMAT(*),TMAT(*),ECORE
         INTEGER NBASISMAX(*),G1(*)

         INTEGER I,J,K,L,M,ISPN
         INTEGER NELEX,NELNEW,IPOSO,IPOSN,ISGNCH
         REAL*8 ECUR,SCRRES,GETHELEMENT2T,TOT
            CALL AZZERO(DEDCIJ,NSBASIS*NSBASIS*NSPINS)
            DO ISPN=1,NSPINS
             DO I=1,NELS(ISPN)
C.. First find out which new basis function index this is
               NELEX=(I-1)*NSPINS+1+ISPN-1
               DO J=1,NSBASIS
C.. first part is d_ji <Psi|H|Psi>.  d=cT
                  TOT=ECUR*CMAT(I,J,ISPN)
C.. now add in factor for all single replacements of NELEX
C.. We only need to look at the same spin, as we multiply by dki, which
C.. is zero if i and k have different spins
                  DO K=NELS(ISPN)+1,NSBASIS
                     NELNEW=(K-1)*NSPINS+1+ISPN-1
                     M=1
                     DO L=1,NEL
                        IF(NDET1(L).EQ.NELEX) THEN
                           IPOSO=L
                        ELSEIF(NDET1(L-1).LT.NELNEW
     &                     .AND.NDET1(L).GT.NELNEW) THEN
C.. We slot in the new det here
                           NDET2(M)=NELNEW
                           IF(M.LT.NEL) NDET2(M+1)=NDET1(L)
                           IPOSN=M
                           M=M+2
                        ELSE
                           NDET2(M)=NDET1(L)
                           M=M+1
                        ENDIF
                     ENDDO
C.. If we haven't yet put it the new electron, we put it at the end
                     IF(M.EQ.NEL) NDET2(NEL)=NELNEW

C.. calculate how many positions move from Old to new electron, and work
C.. out appropriate sign change
                     ISGNCH=(-1)**(IPOSN-IPOSO+NEL)
C.. At this point, NDET2 contains a det which is Psi with  phi_i
C.. replaced by phi_k and rerdered.  
C.. Now calculate <NDET2 | H | PSI>
                     SCRRES=GETHELEMENT2T(NDET1(1),NDET2,NEL,NBASISMAX,
     &                  G1,NBASIS,TMAT,UMAT,ECORE,1,CMAT,NSBASIS,NSPINS)
C.. We multiply by the amount of phi_k in u_j, as well as the sign
                     TOT=TOT+ISGNCH*SCRRES*CMAT(K,J,ISPN)
                  ENDDO
C.. TOT now contains dE/dcij, so we store this
                  DEDCIJ(I,J,ISPN)=2*TOT
               ENDDO
              ENDDO
            ENDDO
         RETURN
      END

      SUBROUTINE CALCDEDCIJ2(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,TMAT,UMAT,
     &         ECORE, NBASIS,G1,NBASISMAX,NELS,NEL,ECUR,ISS)            
         IMPLICIT NONE
         REAL*8 CMAT(NSBASIS,NSBASIS,NSPINS)
         REAL*8 DEDCIJ(NSBASIS,NSBASIS,NSPINS)
         INTEGER NSBASIS,NSPINS,NBASIS
         INTEGER NDET1(0:NEL+1),NELS(NSPINS),NEL
         REAL*8 UMAT(*)
         REAL*8 TMAT(NBASIS,NBASIS),ECORE
         INTEGER NBASISMAX(*),G1(*),ISS
         INTEGER ISPN,I,J,K,A,B,C,F,JSPN,KSPN,JJ
         INTEGER IDA,IDB,IDC,IDF
         REAL*8 TOT,TOT1,TOT2,ECUR,TOT1B
         INTEGER UMatInd
C.. We calculate dE/dcij as
C.. dE/dc_kf = 2(Sum_a c_ka <a|h|f> 
C.. +Sum_j Sum_abc c_ja c_kb c_jc (<af|U|cb>-<af|U|bc>+<ab|U|cf>-<ab|U|fc>)) 
         CALL AZZERO(DEDCIJ,NSBASIS*NSBASIS*NSPINS)
         DO KSPN=1,NSPINS
C.. K is an HF orbital
            DO K=1,NELS(KSPN)
C.. F is the basis orbital
               DO F=1,NSBASIS
C.. TMAT ID
                  IDF=(F-1)*NSPINS+1+KSPN-1
                  TOT=0.D0
C.. deal with the one-electron integrals first
                  DO A=1,NSBASIS
                     IDA=(A-1)*NSPINS+1+KSPN-1
                     TOT=TOT+CMAT(K,A,KSPN)*TMAT(IDA,IDF)
                  ENDDO
C.. UMAT ID
                  CALL GTID(NBASISMAX,(F-1)*NSPINS+1+KSPN-1,IDF)
                  DO I=1,NEL
                     JJ=NDET1(I)
                     JSPN=MOD(JJ-1,NSPINS)+1
                     J=(JJ-1)/NSPINS+1
C.. Here A and C correspond to J and B corresponds to K
                     DO A=1,NSBASIS
                        CALL GTID(NBASISMAX,(A-1)*NSPINS+1+JSPN-1,IDA)
                        TOT1=0.D0
                        DO B=1,NSBASIS
                           CALL GTID(NBASISMAX,
     &                        (B-1)*NSPINS+1+KSPN-1,IDB)
                           TOT1B=0.D0
                           DO C=1,NSBASIS
                              TOT2=0.D0
                              CALL GTID(NBASISMAX,
     &                           (C-1)*NSPINS+1+JSPN-1,IDC)
                              TOT2=TOT2+UMAT(UMatInd(IDA,IDF,IDC,IDB))
                              TOT2=TOT2+UMAT(UMatInd(IDA,IDB,IDC,IDF))
                              IF(KSPN.EQ.JSPN) THEN                  
                                TOT2=TOT2-UMAT(UmatInd(IDA,IDF,IDB,IDC))
                                TOT2=TOT2-UMAT(UMatInd(IDA,IDB,IDF,IDC))
                              ENDIF
                              TOT1B=TOT1B+TOT2*CMAT(J,C,JSPN)
                           ENDDO
                           TOT1=TOT1+TOT1B*CMAT(K,B,KSPN)
                        ENDDO
                        TOT=TOT+TOT1*CMAT(J,A,JSPN)
                     ENDDO
                  ENDDO
                  DEDCIJ(K,F,KSPN)=2*TOT
               ENDDO
            ENDDO
         ENDDO
         RETURN
      END

      SUBROUTINE READHFFMAT(NBASIS,FMAT,HFES,G1,NSPINS,NSBASIS,TRANSP)
         IMPLICIT NONE
         INTEGER NBASIS,G1(*),NQNS(5),NN,NSPINS,NSBASIS
         REAL*8 FMAT(NSBASIS,NSBASIS,NSPINS),HFES(NSBASIS,NSPINS)
         INTEGER I,L,J,NB,NE,IG,N
         REAL*8 VAL
         INTEGER IFINDBASISFN
         LOGICAL TRANSP
         WRITE(6,*) "Loading HF BASIS"
         OPEN(10,FILE='HFBASIS',STATUS='OLD')
         READ(10,*)
         READ(10,*) NB,NE
C.. NE is NEVAL, and NB is NBASIS/2
C.. NBASIS is the number of orbitals, so *2 to get # spinorbitals
C         IF(NE.NE.NEL) STOP 'NEL in HFBASIS <> NEL'

         IF(NE.NE.NB) STOP 'NEVAL <> NBASIS in HFBASIS not supported'
         IF(NB*2.NE.NBASIS) STOP 'NBASIS in HFBASIS <> NHG'
         DO I=1,NB
            DO L=-1,1,2
               READ(10,*)
               READ(10,*) HFES(I,(L+3)/2)
               NQNS(4)=L
               DO J=1,NB
                  READ(10,*) NN,NQNS(1),NQNS(2),NQNS(3),VAL
                  IG=IFINDBASISFN(NQNS,G1,NBASIS)
                  IF(TRANSP) THEN
                     FMAT(I,J,(L+3)/2)=VAL
                  ELSE
                     FMAT(J,I,(L+3)/2)=VAL
                  ENDIF
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
                ENDDO
            ENDDO
         ENDDO
         CLOSE(10)
         RETURN
      END
      SUBROUTINE ORDERBASISHF(ARR,BRR,HFE,HFBASIS,G1,NBASIS,FDET,NEL)
         IMPLICIT NONE
         INTEGER NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS,2),HFE(NBASIS),HFBASIS(NBASIS,NBASIS)
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         INTEGER I,J,IBF,G1(*),ITOT,NEL,FDET(NEL),ICUR
         REAL*8 MX,OEN
         ICUR=1
         DO I=1,NBASIS
            IF(ICUR.LT.NEL.AND.FDET(ICUR+1).EQ.I) ICUR=ICUR+1
            MX=0.D0
            IBF=0
            DO J=1,NBASIS
               IF(ABS(HFBASIS(I,J)).GT.MX) THEN
                  MX=ABS(HFBASIS(I,J))
                  IBF=J
               ENDIF
            ENDDO
            IF(I.EQ.FDET(ICUR).AND.MX.LT.0.95D0) THEN
C.. The largest element isn't big enough, so we abort
               WRITE(6,*) "Largest coeff of HF basis fn ",I," is ",MX
               WRITE(6,*) "Aborting ORDERBASISHF"
               CALL FLUSH(6)
               STOP "ORDERBASISHF failed - HF Basis not converged"
            ENDIF
            ARR(I,1)=HFE(I)
            ARR(IBF,2)=HFE(I)
            BRR(I)=IBF
         ENDDO
C.. We need to now go through each set of degenerate orbitals, and make
C.. the correct ones are paired together in BRR otherwise bad things
C.. happen in FREEZEBASIS
C.. We do this by ensuring that within a degenerate set, the BRR are in
C.. ascending order
         OEN=ARR(1,1)
         J=1
C         G1(3,BRR(1))=J
         ITOT=1
         DO I=2,NBASIS
            IF(ABS(ARR(I,1)-OEN).GT.1.D-4) THEN
C.. We don't have degenerate orbitals
C.. First deal with the last set of degenerate orbitals
C.. We sort them into order of BRR
               CALL SORT2_(ITOT,BRR(I-ITOT),ARR(I-ITOT,1))
C.. now setup the new degenerate set.
               J=J+1
               ITOT=1
            ELSE
               ITOT=ITOT+1
            ENDIF
            OEN=ARR(I,1)
C.. If we've got a generic spatial sym or hf we mark degeneracies
C               G1(3,BRR(I))=J
         ENDDO
         RETURN
      END
