C.. Ti=Sum_a,b(c_ia* c_ib <u_a(1)|h(1)|u_b(1)>)
C..   =Sum_a(|c_ia|^2 <u_a|(p^2)/2|u_a>)
      SUBROUTINE CALCHFTMAT(NBASIS,HFBASIS,NORBUSED)
        use constants, only: dp
        USE OneEInts, only : GetTMATEl,TMAT2D,TMAT2D2,TSTARSTORE
         IMPLICIT NONE
         INTEGER :: NBASIS
         real(dp) HFBASIS(NBASIS,NBASIS)
         INTEGER I,J,A,B,NORBUSED
         real(dp) SUM1
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         OPEN(10,FILE='TMAT2',STATUS='UNKNOWN')
         DO I=1,NORBUSED
            DO J=1,NORBUSED
               SUM1=0.0_dp
               DO A=1,NBASIS
                  DO B=1,NBASIS
                    SUM1=SUM1+HFBASIS(I,A)*HFBASIS(J,B)*
     &                   (GetTMATEl(A,B))
                  ENDDO
               ENDDO
               TMAT2D2(I,J)=SUM1
               IF(ABS(SUM1).GT.1.0e-10_dp) WRITE(10,*) I,J,TMAT2D2(I,J)
            ENDDO
         ENDDO
         CLOSE(10)
         RETURN
      END
C.. Ti=Sum_a,b(c_ia* c_ib <u_a(1)|h(1)|u_b(1)>)
C..   =Sum_a(|c_ia|^2 <u_a|(p^2)/2|u_a>)
      SUBROUTINE READHFTMAT(NBASIS)
         use constants, only: dp
         USE OneEInts, only : TMAT2D2
         IMPLICIT NONE
         INTEGER NBASIS
         real(dp) R
         INTEGER I,J
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         OPEN(10,FILE='TMAT2',STATUS='OLD')
         I=0
         DO WHILE(.NOT.(I.EQ.NBASIS.AND.J.EQ.NBASIS))
            READ(10,*,END=11) I,J,R
            TMAT2D2(I,J)=R
            TMAT2D2(J,I)=R
         ENDDO
11       CLOSE(10)
         RETURN
      END
      SUBROUTINE READHFUMAT(UMAT2,NBASIS)
         USE UMatCache , only : UMatInd
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NBASIS
         real(dp) UMAT2(*)
         real(dp) SUM1
         INTEGER A,B,C,D
         INTEGER NHG
         NHG=NBASIS
         WRITE(6,*) 'READING HF UMAT'
C ==--------------------------------------------------------------------=
         OPEN(10,FILE='UMAT2',STATUS='OLD')
         DO WHILE (.TRUE.)
            
            READ(10,*,END=12) A,B,C,D,SUM1
              UMAT2(UMatInd(A,B,C,D,0,0))=SUM1
C..Symmetries not needed by UMatInd
!              UMAT2(C,D,A,B)=SUM1
!              UMAT2(B,A,D,C)=SUM1
!              UMAT2(D,C,B,A)=SUM1
         ENDDO
12       CLOSE(10)
      END

      SUBROUTINE SETUPHFBASIS(NBASISMAX,G1,NBASIS,HFE,ARR,BRR)
         use SystemData, only: Symmetry,BasisFN
         use constants, only: dp
         use sym_mod
         IMPLICIT NONE
         INTEGER NBASIS,nBasisMax(5,*)
         TYPE(BasisFN) G1(nBasis)
         real(dp) ARR(NBASIS,2)
         real(dp) HFE(NBASIS)
         INTEGER BRR(NBASIS),ORBORDER(8)
         INTEGER I
#ifdef __CMPLX
         call stop_all('SETUPHFBASIS', 
     &                 'HF not implemented for complex orbitals.')
#endif
 
C.. We now need to modify G1.  Pretend that all basis functions have
C.. a different X quantum number, and all have the same spin.
         DO I=1,NBASIS
            G1(I)%k(1)=0
            G1(I)%Ms=-2*(MOD(I,2))+1
            G1(I)%Sym=TotSymRep()
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
         CALL GENMOLPSYMTABLE(1,G1,NBASIS)
         CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
         CALL GENMOLPSYMREPS()

         RETURN
      END
      SUBROUTINE CALCHFUMAT(UMAT,UMAT2,NBASIS,HFBASIS,ISS,NORBUSED)
         USE UMatCache , only : UMatInd,GTID
         use constants, only: dp
         use global_utilities
         IMPLICIT NONE
         INTEGER NBASIS,ISS
         real(dp) UMAT(*)
         real(dp) UMAT2(*)
         real(dp) UMATT((((NBASIS/ISS)*(NBASIS/ISS-1))/2)**2)
         real(dp) HFBASIS(NBASIS,NBASIS),SUM1
         INTEGER I,J,K,L,A,B,C,D
         type(timer), save :: proc_timer
         INTEGER ID1,ID2,ID3,ID4,NHG,NORBUSED
         LOGICAL LSPN
         NHG=NBASIS
         proc_timer%timer_name='CALCHFUMAT'
         call set_timer(proc_timer)
         WRITE(6,*) 'CALCULATING HF UMAT'
       STOP "HF UMAT calc broken through U/TMAT reindexing. Please fix"
#ifdef __CMPLX
         call stop_all('CALCHFUMAT', 
     &                 'HF not implemented for complex orbitals.')
#endif
C ==--------------------------------------------------------------------==
         OPEN(10,FILE='UMAT2',STATUS='UNKNOWN')
C.. A, B, C, D denote basis fns in the HF basis.  These basis fns
C.. alternate in spin i.e. fn 1 has alpha, fn 2 beta, fn 3 alpha etc.
C.. We need to take into account the spin when writing out the new U
C.. matrix.  For <A(1) B(2) |U| C(1) D(2)> to be non-zero, A and C must
C.. have the same spin, and B and D must have the same spin.
C.. Thus the sum1 A+C must be even, as must B+D
         WRITE(6,*) "Index 1..."
         DO A=1,NORBUSED
          DO J=1,NHG/ISS
           DO K=1,NHG/ISS
            DO L=1,NHG/ISS
             SUM1=0.0_dp
             DO I=1,NHG
              ID1 = GTID(I)
              SUM1=SUM1+HFBASIS(A,I)
     &               *UMAT(UMatInd(ID1,J,K,L,0,0))
             ENDDO
             UMATT(UMatInd(J,K,L,A,0,0))=SUM1
!             IF(ABS(SUM1).GT.1.0e-9_dp) WRITE(6,*) J,K,L,A,SUM1
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         WRITE(6,*) "Index 2..."
         DO A=1,NORBUSED
          DO B=1,NORBUSED
           DO K=1,NHG/ISS
            DO L=1,NHG/ISS
             SUM1=0.0_dp
             DO J=1,NHG
              ID2 = GTID(J)
              SUM1=SUM1+HFBASIS(B,J)
     &               *UMATT(UMatInd(ID2,K,L,A,0,0))
             ENDDO
             UMAT2(UMatInd(K,L,A,B,0,0))=SUM1
!             IF(ABS(SUM1).GT.1.0e-9_dp) WRITE(6,*) K,L,A,B,SUM1
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         WRITE(6,*) "Index 3..."
         DO A=1,NORBUSED
          DO B=1,NORBUSED
           DO C=1,NORBUSED
            DO L=1,NHG/ISS
             SUM1=0.0_dp
             DO K=1,NHG
              ID3 = GTID(K)
              SUM1=SUM1+HFBASIS(C,K)
     &               *UMAT2(UMatInd(ID3,L,A,B,0,0))
             ENDDO
             UMATT(UMatInd(L,A,B,C,0,0))=SUM1
!             IF(ABS(SUM1).GT.1.0e-9_dp) WRITE(6,*) L,A,B,C,SUM1
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         UMAT2(1:nBasis**4)=0.0_dp
         WRITE(6,*) "Index 4..."
         DO A=1,NORBUSED
          DO B=1,NORBUSED
           DO C=1,NORBUSED
            DO D=1,NORBUSED
             LSPN=(MOD(A+C,2).EQ.0).AND.(MOD(B+D,2).EQ.0)
             IF(LSPN.AND.UMAT2(UMatInd(A,B,C,D,0,0)).EQ.0.0_dp) THEN
              SUM1=0.0_dp
              DO L=1,NHG
               ID4 = GTID(L)
               SUM1=SUM1+HFBASIS(D,L)
     &                *UMATT(UMatInd(ID4,A,B,C,0,0))
              ENDDO
              UMAT2(UMatInd(A,B,C,D,0,0))=SUM1
C..Symmetries not needed by UMatInd
!              UMAT2(C,D,A,B)=SUM1
!              UMAT2(B,A,D,C)=SUM1
!              UMAT2(D,C,B,A)=SUM1
               IF(ABS(SUM1).GT.1.0e-10_dp)
     &          WRITE(10,'(4I7,F19.9)') A,B,C,D,SUM1
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         CLOSE(10)
         WRITE(6,*) ' !!! FINISHED CALCULATING HF UMAT !!! '
         call halt_timer(proc_timer)
         RETURN
      END
      SUBROUTINE READHFBASIS(HFBASIS,HFE,G1,NBASIS)
         use SystemData, only: BasisFN
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NBASIS,NQNS(5),NN
         TYPE(BasisFN) G1(nBasis)
         real(dp) HFBASIS(NBASIS,NBASIS),HFE(NBASIS)
         INTEGER I,L,J,NB,NE,IG
         real(dp) VAL
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



      SUBROUTINE CALCHFBASIS(NBASIS,NBASISMAX,G1,BRR,ECORE,
     &      UMAT,HFE,HFBASIS,NHFIT,NEL,MS,HFMIX,EDELTA,CDELTA,TRHF,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         use constants, only: dp
         use HElem, only: HElement_t_size
         use global_utilities
         use SystemData, only: BasisFN
         use MemoryManager, only: TagIntType
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS
         INTEGER NBASIS,nBasisMax(5,*)
         TYPE(BasisFN) G1(nBasis)
         real(dp) UMAT(*)
         real(dp) ECORE
         INTEGER BRR(NBASIS)
         HElement_t HFBASIS(NBASIS,NBASIS)
         real(dp) HFE(NBASIS)
         HElement_t,allocatable :: FMAT(:),OFMAT(:)
         HElement_t,allocatable :: DMAT(:),ODMAT(:)
         HElement_t,allocatable :: WORK(:)
         real(dp),allocatable :: HFES(:)
         HElement_t,allocatable :: R1(:),R2(:)
         integer(TagIntType), save :: tagR1=0,tagR2=0,tagHFES=0
         integer(TagIntType), save :: tagFMAT=0,tagOFMat=0
         integer(TagIntType), save :: tagWork=0,tagDMAT=0,tagODMat=0
         INTEGER NHFIT,NEL
         real(dp) HFMIX,EDELTA,CDELTA
         INTEGER MS,IHFMETHOD
         LOGICAL TRHF,TREADHF
         real(dp) FRAND
         INTEGER HFDET(NEL)
         INTEGER ILOGGING
         character(*),parameter :: this_routine='CALCHFBASIS'
C.. Allocate mem etc.
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
C         IF(TRHF) STOP 'RHF not supported'
         NSPINS=1+(NBASISMAX(4,2)-NBASISMAX(4,1))/2
         NSBASIS=NBASIS/NSPINS
C         IF(TRHF) NSPINS=1
         allocate(FMAT(NSBASIS*NSBASIS*NSPINS))
         call LogMemAlloc('FMAT',NSBASIS*NSBASIS*NSPINS,
     &                           HElement_t_size*8,this_routine,tagFMAT)
         allocate(OFMAT(NSBASIS*NSBASIS*NSPINS))
         call LogMemAlloc('OFMAT',NSBASIS*NSBASIS*NSPINS,
     &                          HElement_t_size*8,this_routine,tagOFMAT)
         allocate(DMAT(NSBASIS*NSBASIS*NSPINS))
         call LogMemAlloc('DMAT',NSBASIS*NSBASIS*NSPINS,
     &                           HElement_t_size*8,this_routine,tagDMAT)
         allocate(ODMAT(NSBASIS*NSBASIS*NSPINS))
         call LogMemAlloc('ODMAT',NSBASIS*NSBASIS*NSPINS,
     &                          HElement_t_size*8,this_routine,tagODMAT)
         allocate(WORK(NSBASIS*3))
         call LogMemAlloc('WORK',NSBASIS*3,HElement_t_size*8,
     &                                             this_routine,tagWORK)
         allocate(R1(NSBASIS*NSBASIS))
         call LogMemAlloc('R1',NSBASIS*NSBASIS,HElement_t_size*8,
     &                                             this_routine,tagR1)
         allocate(R2(NSBASIS*NSBASIS))
         call LogMemAlloc('R2',NSBASIS*NSBASIS,HElement_t_size*8,
     &                                             this_routine,tagR2)
         allocate(HFES(NSBASIS*NSPINS))
         call LogMemAlloc('HFES',NSBASIS*NSPINS,8,this_routine,tagHFES)
C.. Generate initial HFBASIS vectors as the energy ordered single
C.. particle basis fns, separated into up and down blocks
         FMAT=(0.0_dp)
         DMAT=(0.0_dp)
         IF(IHFMETHOD.EQ.0.OR.IHFMETHOD.EQ.-1) THEN
            CALL UHFSCF(NBASIS,G1,BRR,ECORE,
     &         HFE,HFBASIS,NHFIT,NEL,MS,FMAT,DMAT,ODMAT,WORK,
     &         NSPINS,NSBASIS,HFES,OFMAT,HFMIX,EDELTA,CDELTA,TRHF,R1,R2,
     &         IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         ELSE
            CALL UHFGRADDESC(NBASIS,NBASISMAX,G1,BRR,ECORE,
     &         UMAT,HFE,HFBASIS,NHFIT,NEL,MS,NSPINS,NSBASIS,HFES,
     &      HFMIX,FMAT,OFMAT,DMAT,ODMAT,EDELTA,CDELTA,R1,R2,WORK,TRHF,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         ENDIF

         deallocate(FMAT,OFMAT,DMAT,ODMAT,WORK,R1,R2,HFES)
         call LogMemDealloc(this_routine,tagFMAT)
         call LogMemDealloc(this_routine,tagOFMAT)
         call LogMemDealloc(this_routine,tagDMAT)
         call LogMemDealloc(this_routine,tagODMAT)
         call LogMemDealloc(this_routine,tagWORK)
         call LogMemDealloc(this_routine,tagR1)
         call LogMemDealloc(this_routine,tagR2)
         call LogMemDealloc(this_routine,tagHFES)
         
         RETURN
      END
C.. Unrestricted HF SCF code
      SUBROUTINE UHFSCF(NBASIS,G1,BRR,ECORE,
     &      HFE,HFBASIS,NHFIT,NEL,MS,FMAT,DMAT,ODMAT,WORK,
     &      NSPINS,NSBASIS,HFES,OFMAT,HFMIX,EDELTA,CDELTA,TRHF,R1,R2,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         use SystemData, only: BasisFn
         use constants, only: dp
         use HElem, only: HElement_t_size
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS
         INTEGER NBASIS
         TYPE(BasisFn) G1(*)
         real(dp) ECORE
         INTEGER BRR(NBASIS)
         HElement_t HFBASIS(NBASIS,NBASIS)
         real(dp) HFES(NSBASIS,NSPINS),HFE(NBASIS)
         HElement_t FMAT(NSBASIS,NSBASIS,NSPINS)
         HElement_t OFMAT(NSBASIS,NSBASIS,NSPINS)
         HElement_t DMAT(NSBASIS,NSBASIS,NSPINS)
         HElement_t ODMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) WORK(NBASIS*3),HFMIX
         real(dp) R1(NSBASIS,NSBASIS),R2(NSBASIS,NSBASIS)
         INTEGER NHFIT,NEL
         INTEGER I,J,K,L,ISPN,NELS(NSPINS)
         real(dp) ELAST,EDELTA,ECUR,CDELTA
         INTEGER IHFIT
         INTEGER MS,IHFMETHOD,IRHFB
         real(dp) F
         real(dp) RMSD,FRAND
         LOGICAL BR
         LOGICAL TRHF,TREADHF
         INTEGER HFDET(*)
         INTEGER ILOGGING
         F=HFMIX
C         EDELTA=1.0e-8_dp
         ELAST=1.D20
         WRITE(6,*) "Performing Hartree-Fock SCF diagonalisation..."
         IF(IHFMETHOD.EQ.-1) THEN
            WRITE(6,*) "Method -1:Rotational Mixing "
         ELSEIF(IHFMETHOD.EQ.0) THEN
            WRITE(6,*) "Method 0:Linear Mixing"
         ENDIF
         WRITE(6,*) "HF Mixing",HFMIX
         WRITE(6,*) "E Thresh:",EDELTA
         WRITE(6,*) "RMSD Thresh:",CDELTA
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
          CALL GENHFGUESS(FMAT,NSPINS,NSBASIS,BRR,G1,.FALSE.,
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
            CALL DCOPY(NSBASIS*NSBASIS*HElement_t_size,
     &           FMAT(1,1,IRHFB),1,FMAT(1,1,3-IRHFB),1)
            CALL DCOPY(NSBASIS*NSBASIS*HElement_t_size,
     &           DMAT(1,1,IRHFB),1,DMAT(1,1,3-IRHFB),1)
         ENDIF 
       CALL DCOPY(NSBASIS*NSBASIS*NSPINS*HElement_t_size,
     &          DMAT,1,ODMAT,1)
       CALL DCOPY(NSBASIS*NSBASIS*NSPINS*HElement_t_size,
     &          FMAT,1,OFMAT,1)
C.. Construct the Density Matrix
         CALL GENDMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,.FALSE.)
C.. See how much our density matrix has changed from last time
         RMSD=0.0_dp
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS 
               DO J=1,NSBASIS
                  IF(.NOT.TRHF.OR.TRHF.AND.ISPN.EQ.IRHFB) 
     &               RMSD=RMSD+abs(DMAT(I,J,ISPN)-ODMAT(I,J,ISPN))**2
               ENDDO
            ENDDO
         ENDDO
         RMSD=SQRT(RMSD/(NSBASIS*NSBASIS*NSPINS))
C.. replace our HF orbitals in FMAT with the Fock matrix
C.. FMAT just stores the results and the HF orbitals (currently) in it
C.. are not used in calculating the F matrix
         CALL GENFMAT(FMAT,DMAT,NSBASIS,NSPINS)
         CALL DIAGFMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,HFES,WORK,
     &         ECORE,ECUR)
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
         IF(ABS(ECUR-ELAST).LE.EDELTA.AND.RMSD.LE.CDELTA
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
      HFBASIS=(0.0_dp)
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
         use Integrals_neci, only: GetUMatEl
         USE LoggingData, only: HFLogLevel
         use SystemData, only: BasisFN
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS,NELS(NSPINS)
         HElement_t FMAT(NSBASIS,NSBASIS,NSPINS)
         HElement_t DMAT(NSBASIS,NSBASIS,NSPINS)
         LOGICAL LTRANS

         INTEGER ISPN,I,J,K
         HElement_t TOT
C.. Construct the Density Matrix
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS
               DO J=I,NSBASIS
                  TOT=0.0_dp
C.. Sum over occupied orbitals for our spin
                  DO K=1,NELS(ISPN)
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
#ifdef __CMPLX
                     IF(LTRANS) THEN
                        TOT=TOT+CONJG(FMAT(K,I,ISPN))*FMAT(K,J,ISPN)
                     ELSE
                        TOT=TOT+FMAT(I,K,ISPN)*CONJG(FMAT(J,K,ISPN))
                     ENDIF
#else
                     IF(LTRANS) THEN
                        TOT=TOT+(FMAT(K,I,ISPN))*FMAT(K,J,ISPN)
                     ELSE
                        TOT=TOT+FMAT(I,K,ISPN)*(FMAT(J,K,ISPN))
                     ENDIF
#endif
                  ENDDO
                  DMAT(I,J,ISPN)=TOT
                  DMAT(J,I,ISPN)=TOT
               ENDDO
            ENDDO 
            IF(HFLogLevel.gt.0) 
     &      CALL WRITE_HEMATRIX("DMAT",NSBASIS,NSBASIS,DMAT(1,1,ISPN))
         ENDDO
         RETURN
      END 
      SUBROUTINE GENFMAT(FMAT,DMAT,NSBASIS,NSPINS)
         USE OneEInts, only: GetTMATEl
         USE UMatCache, only: GTID
         use Integrals_neci, only: GetUMatEl
         USE LoggingData, only: HFLogLevel
         use constants, only: dp
         use SystemData, only: BasisFN
         use IntegralsData, only: DMatEpsilon
         IMPLICIT NONE
         INTEGER NSBASIS,NSPINS
         HElement_t FMAT(NSBASIS,NSBASIS,NSPINS)
         HElement_t DMAT(NSBASIS,NSBASIS,NSPINS)
!         real(dp) TMAT(NSBASIS*NSPINS,NSBASIS*NSPINS)
         HElement_t TOT
         INTEGER I,J,K,L,ISPN,ISPN2,ID1,ID2,ID3,ID4
         real(dp) RHFMUL
C         IF(TRHF) THEN
C            RHFMUL=2.0_dp
C         ELSE
            RHFMUL=1.0_dp
C         ENDIF
C.. Construct the Fock Matrix
         FMAT=(0.0_dp)
C. <ui|F|uj>=tij+Sum_kl(Dkl <ui uk|U|uj ul>-s<ui uk|U|ul uj>) where s=1 if i,k same spin
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS
               DO J=I,NSBASIS
                  TOT=GetTMATEl((I-1)*NSPINS+ISPN,
     &                     (J-1)*NSPINS+ISPN)
!                  WRITE(6,*) I,J,TOT
C.. Now sum in the alpha and beta block of u matrix
                  ID1 = GTID((I-1)*NSPINS+ISPN)
                  ID2 = GTID((J-1)*NSPINS+ISPN)
                  DO ISPN2=1,NSPINS
                     DO K=1,NSBASIS
                        DO L=1,NSBASIS
                           ID3 = GTID((K-1)*NSPINS+ISPN2)
                           ID4 = GTID((L-1)*NSPINS+ISPN2)
                           IF(abs(DMAT(K,L,ISPN2)).gt.DmatEpsilon)
     &                      TOT=TOT+(RHFMUL)*
     &                  DMAT(K,L,ISPN2)*GetUMatEl(ID1,ID3,ID2,ID4)
!                           WRITE(6,"(4I3)", advance='no') I,J,K,L
!                           WRITE(6,*) (RHFMUL)*
!     &                  DMAT(K,L,ISPN2)*GetUMatEl(NBasisMax,UMat,ALAT,
!     &                     NSBASIS*NSPINS,ISS,G1,ID1,ID3,ID2,ID4)
                           IF(ISPN2.EQ.ISPN
     &                       .AND.(abs(DMAT(K,L,ISPN2)).gt.DmatEpsilon))
     &                        TOT=TOT-
     &                  DMAT(K,L,ISPN2)*GetUMatEl(ID1,ID3,ID4,ID2)
!                           IF(ISPN2.EQ.ISPN)
!     &                        WRITE(6,*)
!     &                  DMAT(K,L,ISPN2)*GetUMatEl(NBasisMax,UMat,ALat,
!     &                     NSBASIS*NSPINS,ISS,G1,ID1,ID3,ID4,ID2)
                        ENDDO
                     ENDDO
                  ENDDO
                  FMAT(I,J,ISPN)=TOT
#ifdef __CMPLX
                  FMAT(J,I,ISPN)=CONJG(TOT)
#endif
               ENDDO
            ENDDO
            IF(HFLogLevel.gt.0) 
     &      CALL WRITE_HEMATRIX("FMAT",NSBASIS,NSBASIS,FMAT(1,1,ISPN))
         ENDDO
         RETURN
      END
      SUBROUTINE DIAGFMAT(NSPINS,NSBASIS,NELS,FMAT,DMAT,HFES,WORK,
     &         ECORE,ECUR)
         USE OneEInts, only : GetTMatEl
         USE LoggingData, only :HFLogLevel
         use constants, only: dp
         use HElem, only: HElement_t_size
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS,NELS(NSPINS)
         real(dp) HFES(NSBASIS,NSPINS)
         HElement_t FMAT(NSBASIS,NSBASIS,NSPINS)
         HElement_t DMAT(NSBASIS,NSBASIS,NSPINS)
!         real(dp) TMAT(NSBASIS*NSPINS,NSBASIS*NSPINS)
         real(dp) ECUR,ECORE
         HElement_t WORK(3*NSBASIS),RWORK(3*NSBASIS)

         INTEGER*4 INFO
         HElement_t TOT,TOT2
         INTEGER ISPN,I,J,K         

C.. First calculate the HF energy double counting contrib
         TOT=0.0_dp
         DO ISPN=1,NSPINS
            DO J=1,NSBASIS
               DO K=1,NSBASIS
                  TOT=TOT+DMAT(J,K,ISPN)
     &     *(FMAT(J,K,ISPN)-
     &      GetTMATEl((J-1)*NSPINS+ISPN,(K-1)*NSPINS+ISPN))
               ENDDO
            ENDDO
         ENDDO

C.. Now diagonalize the Fock matrix 
         DO ISPN=1,NSPINS
            IF(HElement_t_size.EQ.1) THEN
               CALL DSYEV('V','U',NSBASIS,FMAT(1,1,ISPN),NSBASIS,
     &            HFES(1,ISPN),WORK,3*NSBASIS, INFO)
            ELSE
               CALL ZHEEV('V','U',NSBASIS,FMAT(1,1,ISPN),NSBASIS,
     &            HFES(1,ISPN),WORK,3*NSBASIS,RWORK, INFO)
            ENDIF
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               STOP
            ENDIF
            IF(HFLogLevel.gt.0) 
     &      CALL WRITE_HEMATRIX("EIGVEC",NSBASIS,NSBASIS,FMAT(1,1,ISPN))
            IF(HFLogLevel.gt.0) 
     &      CALL NECI_WRITE_MATRIX("EIGVAL",1,NSBASIS,HFES(1,ISPN))
         ENDDO
C.. now calculate the sum of the occupied Fock orbitals
         TOT2=0.0_dp
         DO ISPN=1,NSPINS
            DO I=1,NELS(ISPN)
               TOT2=TOT2+(HFES(I,ISPN))
            ENDDO
         ENDDO
C.. subtract out the double counting term, and add in the core energy
!         WRITE(6,*) "ECUR=",TOT2-TOT/(2.0_dp)+(ECORE)
         ECUR=(TOT2)-(TOT)/2.0_dp+ECORE
         RETURN
      END
      SUBROUTINE WRITEHFPSIALL(NBASIS,FMAT,HFES,G1,NSPINS,NSBASIS,
     &         TRANSP)
         use SystemData, only: BasisFN
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NBASIS,NSBASIS,NSPINS,I,J,K
         TYPE(BasisFN) G1(nBasis)
         real(dp) FMAT(NSBASIS,NSBASIS,NSPINS),HFES(NSBASIS,NSPINS)
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
      SUBROUTINE GENHFGUESS(FMAT,NSPINS,NSBASIS,BRR,G1,TRANS,
     &            LMS,FRAND,NELS,HFDET)
         use SystemData, only: BasisFN
         use Determinants, only: write_det
         use sort_mod
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NSPINS,NSBASIS
         HElement_t FMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) PI,R
         INTEGER ISPN,I,J,ISEED,irand
         real(dp) RAN2
         INTEGER BRR(NSBASIS*NSPINS),NELR,IREAL,IS
         TYPE(BasisFN) G1(*)
         INTEGER LMS
C.. Working space
         real(dp) FRAND
         LOGICAL TRANS
         INTEGER NELS(NSPINS)
         INTEGER HFDET(*),NEL
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector unless transposed
         PI=3.141592653589793_dp
         ISEED=1007
         irand=1
         R=RAN2(irand)
         FMAT=(0.0_dp)
         WRITE(6,*) "Generating HF Guess..."
         NEL=0
         DO IS=1,NSPINS
            IF(LMS.LT.0) THEN
               ISPN=IS
            ELSE
               ISPN=NSPINS+1-IS
            ENDIF
            IREAL=1
            WRITE(6,"(A,I2,A)",advance='no') "Spin ",IS,":"
            
            DO I=1,NSBASIS
               DO WHILE(G1(BRR(IREAL))%Ms.NE.(-3+2*ISPN))
                  IREAL=IREAL+1
               ENDDO
               NELR=(BRR(IREAL)-1)/NSPINS+1
               IF(TRANS) THEN
                  FMAT(I,NELR,ISPN)=(1.0_dp)
               ELSE
                  FMAT(NELR,I,ISPN)=(1.0_dp)
               ENDIF
               IF(I.LE.NELS(IS)) THEN
                  WRITE(6,"(I4,A)",advance='no') BRR(IREAL),","
                  NEL=NEL+1
                  HFDET(NEL)=BRR(IREAL)
               ENDIF
               IREAL=IREAL+1
            ENDDO
            WRITE(6,*)
            DO I=1,NSBASIS
               DO J=1,NSBASIS
                  FMAT(I,J,ISPN)=FMAT(I,J,ISPN)+
     &               (FRAND*RAN2(irand))
               ENDDO
            ENDDO
            CALL GRAMSCHMIDT(FMAT(1,1,ISPN),NSBASIS)
C            CALL LOWDIN_ORTH(FMAT(1,1,ISPN),NSBASIS,R1,R2,WORK)
         ENDDO
         call sort (HFDet(1:nel))
         call write_det (6, HFDET, .true.)
      END
C.. Generate an initial determinant of specified spin
      SUBROUTINE GENFDET(BRR,G1,NBASIS,LMS,NEL,FDET)
         use SystemData, only: BasisFN,tMolpro
         use util_mod, only: int_fmt
         use sort_mod
         IMPLICIT NONE
         INTEGER BRR(*)
         TYPE(BasisFN) G1(*)
         INTEGER LMS
C.. Working space
         INTEGER NELS(2)
         INTEGER FDET(*),NEL,NELR,IREAL,IS,ISPN,NBASIS,NSBASIS,I
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector unless transposed
         NSBASIS=NBASIS/2
         NELS(2)=(LMS+NEL)/2
         NELS(1)=NEL-NELS(2)
         WRITE(6,'(4(A,'//int_fmt(nel)//'))') " N_alpha:",NELS(1),      &
     &   " ; N_beta:",  NELS(2)," ; LMS:",LMS," ; NEl:",NEL
         if(tMolpro) then
             !Assume that orbitals are ordered by occupation number
             !Assuming LMS positive, NELS(2) > NELS(1)
             do i=1,NELS(1)*2
                 FDET(i)=i  !Fill up closed shells
             enddo
             do i=1,LMS !Fill up open shells in odd orbitals
                 FDET(NELS(1)*2+i)=NELS(1)*2+(i*2-1)
             enddo
             return
         endif
         NEL=0
         DO IS=1,2
            IF(LMS.GE.0) THEN
               ISPN=3-IS
            ELSE
               ISPN=IS
            ENDIF
            IREAL=1
            DO I=1,NSBASIS
               DO WHILE(G1(BRR(IREAL))%Ms.NE.(-3+2*ISPN))
                  IREAL=IREAL+1
               ENDDO
               NELR=(BRR(IREAL)-1)/2+1
               IF(I.LE.NELS(IS)) THEN
!                  WRITE(6,"(I4,A)",advance='no'),BRR(IREAL),","
                  NEL=NEL+1
                  FDET(NEL)=BRR(IREAL)
               ENDIF
               IREAL=IREAL+1
            ENDDO
         ENDDO
         call sort (fDet(1:nel))
      END

      SUBROUTINE HFLINMIX(FMAT,OFMAT,NSPINS,N,FMIX,R1,R2,WORK)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N,NSPINS
         real(dp) FMIX
         HElement_t FMAT(N,N,NSPINS),OFMAT(N,N,NSPINS)
         real(dp) R1(*),R2(*),WORK(*)
         INTEGER I,J,ISPN
         DO ISPN=1,NSPINS
            DO I=1,N
               DO J=1,N
                  FMAT(I,J,ISPN)=(FMIX)*FMAT(I,J,ISPN)
     &               +(1.0_dp-FMIX)*OFMAT(I,J,ISPN)
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
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NSPINS,N
         real(dp) FMIX,FMAT(N,N,NSPINS)
         real(dp) OFMAT(N,N,NSPINS)
         real(dp) R1(N,N),R2(N,N),WORK(3*N)
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


C.. Work out R = R1=1.0_dp * F * F'T + 0.0_dp*R1
            CALL DGEMM('N','T',N,N,N,1.0_dp,FMAT(1,1,ISPN),N,
     &            OFMAT(1,1,ISPN),N,0.0_dp,R1,N)
C.. now let P=R1=R-I
            DO I=1,N
               R1(I,I)=R1(I,I)-1.0_dp
            ENDDO
C.. R^a = I+aP(I+(a-1)/2 P(I+ (a-2)/3 P (I+...) ) ) F'
C.. Let FMAT be the accumulator
            CALL DCOPY(N*N,OFMAT(1,1,ISPN),1,FMAT(1,1,ISPN),1)
C.. Go up to 2nd order taylor
            DO J=1,0,-1
C.. Set R2=I
               R2=0.0_dp
               DO I=1,N
                   R2(I,I)=1.0_dp
               END DO
C.. Work out R2=(a-J)/(J+1.0_dp) * R1 * FMAT + 1.0_dp*R2.  a=FMIX
               CALL DGEMM('N','N',N,N,N,(FMIX-J)/(J+1.0_dp),
     &            R1,N,FMAT(1,1,ISPN),N,1.0_dp,R2,N)
               CALL DCOPY(N*N,R2,1,FMAT(1,1,ISPN),1)
            ENDDO
C.. Now reorthoganalise, as this is just an approximation    
            CALL LOWDIN_ORTH(FMAT(1,1,ISPN),N,R1,R2,WORK)
         ENDDO
      END


C.. Unrestricted Hartree Fock with a gradient descent method.
      SUBROUTINE UHFGRADDESC(NBASIS,NBASISMAX,G1,BRR,ECORE,
     &      UMAT,HFE,HFBASIS,NHFIT,NEL,MS,NSPINS,NSBASIS,HFES,
     &      HFMIX,CMAT,OCMAT,DEDCIJ,DMAT,EDELTA,CDELTA,R1,R2,WORK,TRHF,
     &      IHFMETHOD,TREADHF,FRAND,HFDET,ILOGGING)
         use SystemData, only: BasisFn
         use constants, only: dp
         use Determinants, only: write_det
         use sort_mod
         IMPLICIT NONE
         TYPE(BasisFn) G1(*)
         INTEGER NSPINS,NSBASIS
         INTEGER NBASIS,nBasisMax(5,*)
         real(dp) UMAT(*)
         real(dp) ECORE
         INTEGER BRR(NBASIS)
         real(dp) HFBASIS(NBASIS,NBASIS),HFE(NBASIS)
         real(dp) HFES(NSBASIS,NSPINS)
         real(dp) CMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) OCMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) DEDCIJ(NSBASIS,NSBASIS,NSPINS)
         real(dp) DMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) WORK(NBASIS*3)
         INTEGER INORDER(100,2),ILOGGING
         real(dp) EORDER(100,2)
C,HFMIX
         real(dp) R1(NSBASIS,NSBASIS),R2(NSBASIS,NSBASIS)
         INTEGER NHFIT,NEL,IHFMETHOD,NELEX2
C,NSTART(NEL)
         INTEGER I,J,K,L,ISPN,NELS(NSPINS)
         real(dp) ELAST,HFMIX,ECUR,FRAND
C         INTEGER INDS(NBASIS)
C         real(dp) TOT,ELAST,EDELTA,ECUR,TOT2,CDELTA
C         INTEGER ID1,ID2,ID3,ID4,INFO,
         INTEGER IHFIT
         INTEGER MS
C         real(dp) F
C         real(dp) RMSD
         LOGICAL BR
         LOGICAL TRHF,TREADHF
         real(dp) TOT,MIX,TOT2
         INTEGER NDET1(0:NEL+1),NSPN(NSPINS),NELEX
         INTEGER IRHFB,JSPN
         real(dp) RMSD, EDELTA,CDELTA,EN,ECUR2
         INTEGER HFDET(*)
         
         real(dp) GETHELEMENT2T
         irhfb = 0
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
         
C         VMAT=0.0_dp
         IF(TREADHF) THEN
           CALL  READHFFMAT(NBASIS,CMAT,HFES,G1,NSPINS,NSBASIS,.TRUE.)
         ELSE
           CALL GENHFGUESS(CMAT,NSPINS,NSBASIS,BRR,G1,.TRUE.,
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
         call write_det (6, NDET1(1), .true.)
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
     &               NBASIS,ECORE,0,CMAT,NSBASIS,NSPINS)
C.. Calculate the Gradient
            IF(IHFMETHOD.EQ.1) THEN
             CALL CALCDEDCIJ(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,
     &         ECORE, NBASIS,NBASISMAX,NELS,NEL,ECUR)            
            ELSEIF(IHFMETHOD.EQ.2) THEN
              CALL CALCDEDCIJ2(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,
     &         UMAT,NELS,NEL)
            ENDIF
C.. DEDCIJ now comtains all elements of dE/dcij
C.. To move down the slope, we subtract a small amount of this from cij,
C.. and re-orthogonalise
C.. HFMIX is -ve
C              WRITE(6,*) ((CMAT(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
C              WRITE(6,*)
C              WRITE(6,*) ((DEDCIJ(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
C.. modify the velocity.
C              R1 = VMAT(:,:,ISPN) + MIX*DEDCIJ(:,:,ISPN)
C              CALL DCOPY(NSBASIS*NSBASIS,R1,1,VMAT(1,1,ISPN),1)
C              R1 = CMAT(:,:,ISPN) + 0.1_dp*VMAT(:,:,ISPN)

C.. Remove the projection of the "force" already in the direction of
C.. the coefficients
            DO ISPN=1,NSPINS
              R1=0.0_dp
              DO I=1,NELS(ISPN)
               TOT=0.0_dp
               DO J=1,NSBASIS
                  TOT=TOT+DEDCIJ(I,J,ISPN)*CMAT(I,J,ISPN)
               ENDDO
               TOT2=0.0_dp
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
C              IF(ABS(ECUR).LT.1.0e-4_dp) MIX=HFMIX
              IF(IHFIT.GT.NHFIT) BR=.FALSE.
              R2 = CMAT(:,:,ISPN) + R1
              CALL DCOPY(NSBASIS*NSBASIS,R2,1,CMAT(1,1,ISPN),1)
C              WRITE(6,*)
C              WRITE(6,*) ((CMAT(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
C              WRITE(6,*)
C              WRITE(6,*)
              CALL LOWDIN_ORTH(CMAT(1,1,ISPN),NSBASIS,R1,R2,WORK)
            ENDDO
            RMSD=0.0_dp
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
         CALL GENFMAT(DEDCIJ,DMAT,NSBASIS,NSPINS)
         CALL DIAGFMAT(NSPINS,NSBASIS,NELS,DEDCIJ,DMAT,HFES,WORK,
     &         ECORE,ECUR2)

            WRITE(6,"(I6)",advance='no') IHFIT
            WRITE(6,*) ECUR,RMSD,ECUR2
C.. DEDCIJ now contains HF orbitals, and ECUR the Fock Energy
C.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
C.. vector





C.. calculate the orbital energies every time
         DO ISPN=1,NSPINS
            DO I=1,NSBASIS
               NELEX=(I-1)*NSPINS+1+ISPN-1
               EN=0.0_dp
               CALL GETTRTMATEL(NELEX,NELEX,CMAT,NSBASIS,
     &               NSPINS,EN)
               HFES(I,ISPN)=EN
               DO JSPN=1,NSPINS
                  DO J=1,NELS(JSPN)
                     NELEX2=(J-1)*NSPINS+1+JSPN-1
                     IF(NELEX.NE.NELEX2) THEN
C.. we're not allowed to count the current electron again
                        CALL GETTRUMATEL(NELEX,NELEX2,NELEX,NELEX2,
     &                     CMAT,NSBASIS,NSPINS,EN)
                        
C                        HFES(I,ISPN)=HFES(I,ISPN)+EN
                        EN=0.0_dp
                        IF(ISPN.EQ.JSPN)
     &                     CALL GETTRUMATEL(NELEX,NELEX2,NELEX2,NELEX,
     &                        CMAT,NSBASIS,NSPINS,EN)
C                        HFES(I,ISPN)=HFES(I,ISPN)-EN
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            DO I=1,NSBASIS
               INORDER(I,ISPN)=I
            ENDDO
            CALL DCOPY(NSBASIS,HFES(1,ISPN),1,EORDER(1,ISPN),1)
            call sort (eorder(1:nsBasis, iSpn), 
     &                 inOrder(1:nsBasis, iSpn))
         ENDDO
      IF(BR) THEN
C.. Now re-order the orbitals
C.. although we don't actually use this info
         DO ISPN=1,NSPINS
            R1=0.0_dp
            DO I=1,NSBASIS
C.. If R1(2,1) is occupied, then postmultiplying C by R1 moves
C.. column 2 in C to column 1.
C.. INORDER(1,ISPN) is the old index of the lowest energy orb
               R1(I,INORDER(I,ISPN))=1
            ENDDO
C.. Work out R2=1.0_dp CMAT*R1+ 0.0_dp*R2
            CALL DGEMM('N','N',NSBASIS,NSBASIS,NSBASIS,1.0_dp,
     &            R1,NSBASIS,CMAT(1,1,ISPN),NSBASIS,
     &            0.0_dp,R2,NSBASIS)
C            CALL DCOPY(NSBASIS*NSBASIS,R2,1,CMAT(1,1,ISPN),1)
            CALL DGEMM('N','N',NSBASIS,NSBASIS,NSBASIS,1.0_dp,
     &            R1,NSBASIS,OCMAT(1,1,ISPN),NSBASIS,
     &            0.0_dp,R2,NSBASIS)
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
      HFBASIS=0.0_dp
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



      SUBROUTINE CALCDEDCIJ(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,
     &         ECORE, NBASIS,NBASISMAX,NELS,NEL,ECUR)            
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NSBASIS,NSPINS,NBASIS,NEL
         real(dp) CMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) DEDCIJ(NSBASIS,NSBASIS,NSPINS)
         INTEGER NDET1(0:NEL+1),NDET2(NEL),NELS(NSPINS)
         real(dp) ECORE
         INTEGER NBASISMAX(*)

         INTEGER I,J,K,L,M,ISPN
         INTEGER NELEX,NELNEW,IPOSO,IPOSN,ISGNCH
         real(dp) ECUR,SCRRES,GETHELEMENT2T,TOT
            DEDCIJ=0.0_dp
            iposo = 0
            iposn = 0
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
     &                  NBASIS,ECORE,1,CMAT,NSBASIS,NSPINS)
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

      SUBROUTINE CALCDEDCIJ2(CMAT,DEDCIJ,NDET1,NSPINS,NSBASIS,UMAT,
     &         NELS,NEL)
         USE UMatCache, only : UMatInd,GTID
         use OneEInts, only: GetTMatEl
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NSBASIS,NSPINS,NEL
         real(dp) CMAT(NSBASIS,NSBASIS,NSPINS)
         real(dp) DEDCIJ(NSBASIS,NSBASIS,NSPINS)
         INTEGER NDET1(0:NEL+1),NELS(NSPINS)
         real(dp) UMAT(*)
         INTEGER I,J,K,A,B,C,F,JSPN,KSPN,JJ
         INTEGER IDA,IDB,IDC,IDF
         real(dp) TOT,TOT1,TOT2,TOT1B
#ifdef __CMPLX
         call stop_all('CALCDEDCIJ2', 
     &                 'HF not implemented for complex orbitals.')
#endif
C.. We calculate dE/dcij as
C.. dE/dc_kf = 2(Sum_a c_ka <a|h|f> 
C.. +Sum_j Sum_abc c_ja c_kb c_jc (<af|U|cb>-<af|U|bc>+<ab|U|cf>-<ab|U|fc>)) 
         DEDCIJ=0.0_dp
         DO KSPN=1,NSPINS
C.. K is an HF orbital
            DO K=1,NELS(KSPN)
C.. F is the basis orbital
               DO F=1,NSBASIS
C.. TMAT ID
                  IDF=(F-1)*NSPINS+1+KSPN-1
                  TOT=0.0_dp
C.. deal with the one-electron integrals first
                  DO A=1,NSBASIS
                     IDA=(A-1)*NSPINS+1+KSPN-1
                     TOT=TOT+CMAT(K,A,KSPN)*(GetTMATEl(IDA,IDF))
                  ENDDO
C.. UMAT ID
                  IDF = GTID((F-1)*NSPINS+1+KSPN-1)
                  DO I=1,NEL
                     JJ=NDET1(I)
                     JSPN=MOD(JJ-1,NSPINS)+1
                     J=(JJ-1)/NSPINS+1
C.. Here A and C correspond to J and B corresponds to K
                     DO A=1,NSBASIS
                        IDA = GTID((A-1)*NSPINS+1+JSPN-1)
                        TOT1=0.0_dp
                        DO B=1,NSBASIS
                           IDB = GTID((B-1)*NSPINS+1+KSPN-1)
                           TOT1B=0.0_dp
                           DO C=1,NSBASIS
                              TOT2=0.0_dp
                              IDC = GTID((C-1)*NSPINS+1+JSPN-1)
                              TOT2=TOT2+UMAT(UMatInd(IDA,IDF,IDC,IDB
     &                                       ,0,0))
                              TOT2=TOT2+UMAT(UMatInd(IDA,IDB,IDC,IDF
     &                                       ,0,0))
                              IF(KSPN.EQ.JSPN) THEN                  
                              TOT2=TOT2-UMAT(UmatInd(IDA,IDF,IDB,IDC
     &                                       ,0,0))
                              TOT2=TOT2-UMAT(UMatInd(IDA,IDB,IDF,IDC
     &                                       ,0,0))
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
         use SystemData, only: BasisFn
         use constants, only: dp
         IMPLICIT NONE
         TYPE(BasisFn) G1(*)
         INTEGER NBASIS,NQNS(5),NN,NSPINS,NSBASIS
         real(dp) FMAT(NSBASIS,NSBASIS,NSPINS),HFES(NSBASIS,NSPINS)
         INTEGER I,L,J,NB,NE,IG
         real(dp) VAL
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
      SUBROUTINE ORDERBASISHF(ARR,BRR,HFE,HFBASIS,NBASIS,FDET,NEL)
         use sort_mod
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NBASIS,BRR(NBASIS)
         real(dp) ARR(NBASIS,2),HFE(NBASIS),HFBASIS(NBASIS,NBASIS)
C.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
         INTEGER I,J,IBF,ITOT,NEL,FDET(NEL),ICUR
         real(dp) MX,OEN
         ICUR=1
         DO I=1,NBASIS
             if (ICUR < NEL) then
                if (FDET(ICUR+1) == I) ICUR = ICUR+1
            endif
            MX=0.0_dp
            IBF=0
            DO J=1,NBASIS
               IF(ABS(HFBASIS(I,J)).GT.MX) THEN
                  MX=ABS(HFBASIS(I,J))
                  IBF=J
               ENDIF
            ENDDO
            IF(I.EQ.FDET(ICUR).AND.MX.LT.0.95_dp) THEN
C.. The largest element isn't big enough, so we abort
               WRITE(6,*) "Largest coeff of HF basis fn ",I," is ",MX
               WRITE(6,*) "Aborting ORDERBASISHF"
               CALL neci_flush(6)
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
            IF(ABS(ARR(I,1)-OEN).GT.1.0e-4_dp) THEN
C.. We don't have degenerate orbitals
C.. First deal with the last set of degenerate orbitals
C.. We sort them into order of BRR
               call sort (brr(i-itot:i-1), arr(i-itot:i-1,1))
!               CALL SORT2_(ITOT,BRR(I-ITOT),ARR(I-ITOT,1))
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

      SUBROUTINE Write_HEMatrix(CHAR,M,N,A)
      use constants, only: dp 
      use HElem, only: HElement_t_size
      IMPLICIT NONE
      CHARACTER(*) CHAR
      Integer I,M,N,J
      HElement_t A(M,N)
      WRITE(6,*) CHAR
      DO I=1,M
         IF(HElement_t_size.EQ.1) THEN
           WRITE(6,"(12E15.6)") (A(I,J),J=1,N)
         ELSE
           WRITE(6,"(6('(',E15.6,',',E15.6,')'))") (A(I,J),J=1,N)
         ENDIF
      ENDDO
c 1000 FORMAT(12E15.6)
      RETURN
      END
