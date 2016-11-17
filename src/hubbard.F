C ===================================================================
      SUBROUTINE CALCUMATHUBREAL(NBASIS,UHUB, UMAT)
         use UMatCache , only : UMatInd
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NBASIS
         real(dp) UMAT(*)
         real(dp) UHUB
         INTEGER I,IC
         OPEN(10,FILE='UMAT',STATUS='UNKNOWN')
         IC=1
!Only do this for each spatial orbital.
         DO I=1,NBASIS/2
            UMAT(UMatInd(I,I,I,I))=UHUB
!            UMAT(UMatInd(I,I+1,I,I+1))=UHUB
!            UMAT(UMatInd(I+1,I,I+1,I))=UHUB
!            UMAT(I,I+1,I,I+1)=UHUB
!            UMAT(I+1,I,I+1,I)=UHUB
            WRITE(10,'(4I7,F19.9)') I,I,I,I,UHUB
         ENDDO
         CLOSE(10)
         RETURN
      END


      SUBROUTINE HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,TOTSUM,TREAL)
C returns the non-interacting energy of state with 
C..quatum numbers (i,j,k) for Hubbard model of 
C..lengths (LX,LY,LZ) PBC.
C..Returned in sum
      use constants, only: Pi, Pi2, THIRD,dp
      use SystemData, only: twisted_bc
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER nBasisMax(5,*),AX,AY
      INTEGER LX,LY,LZ,K,J,I
      LOGICAL TTILT,TREAL
      real(dp) :: TOTSUM,BHUB
      real(dp) :: tb_x, tb_y
      tb_x = twisted_bc(1)
      tb_y = twisted_bc(2)
      LX=NBASISMAX(1,2)-NBASISMAX(1,1)+1
      LY=NBASISMAX(2,2)-NBASISMAX(2,1)+1
      LZ=NBASISMAX(3,2)-NBASISMAX(3,1)+1
      AX=NBASISMAX(1,4)
      AY=NBASISMAX(2,4)
      IF(TREAL) THEN
         TOTSUM=0.0_dp
      ELSE
       IF(TTILT) THEN
CCCC.. NBASISMAX goes from -NMAXX+1 to MAXX so LX=2MAXX
         LX=NBASISMAX(1,5)
         LX=(LX*(AX*AX+AY*AY))
         TOTSUM=COS(2*PI*((I+tb_x)*AX+(J+tb_y)*AY)/LX)
         TOTSUM=TOTSUM+COS(2*PI*((I+tb_x)*AY-(J+tb_y)*AX)/LX)
       ELSE
         IF(LX.GT.1) TOTSUM=COS(2*PI*(I + twisted_bc(1))/(LX))
         IF(LY.GT.1) TOTSUM=TOTSUM+COS(2*PI*(J + twisted_bc(2))/(LY))
         IF(LZ.GT.1) TOTSUM=TOTSUM+COS(2*PI*K/(LZ))
       ENDIF
       TOTSUM=TOTSUM*2.0_dp*BHUB
      ENDIF
      RETURN
      END

      SUBROUTINE HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,TOTSUM,TREAL)
C returns the non-interacting energy of state with 
C..quatum numbers (i,j,k) for Hubbard model of 
C..lengths (LX,LY,LZ).  NON-PBC
C..Returned in sum
      use constants, only: dp
      use constants, only: Pi, Pi2, THIRD
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER nBasisMax(5,*),AX,AY
      LOGICAL TTILT,TREAL
      INTEGER II,JJ,KK,K,J,I,LZ,LX,LY
      real(dp) :: BHUB,TOTSUM
      LX=NBASISMAX(1,2)-NBASISMAX(1,1)+1
      LY=NBASISMAX(2,2)-NBASISMAX(2,1)+1
      LZ=NBASISMAX(3,2)-NBASISMAX(3,1)+1
      AX=NBASISMAX(1,4)
      AY=NBASISMAX(2,4)
      IF(TTILT) THEN
CCCC.. NBASISMAX goes from -NMAXX+1 to MAXX so LX=2MAXX
         LX=NBASISMAX(1,5)
         LX=(LX*(AX*AX+AY*AY))
         TOTSUM=COS(PI*(I*AX-J*AY)/(LX+1))+COS(PI*(I*AY+J*AX)/(LX+1))
      ELSE
         II=I-NBASISMAX(1,1)+1
         JJ=J-NBASISMAX(2,1)+1
         KK=K-NBASISMAX(3,1)+1
         IF(LX.GT.1) TOTSUM=COS(PI*II/(LX+1))
         IF(LY.GT.1) TOTSUM=TOTSUM+COS(PI*JJ/(LY+1))
         IF(LZ.GT.1) TOTSUM=TOTSUM+COS(PI*KK/(LZ+1))
C         IF(LY.GT.1) TOTSUM=TOTSUM+COS(PI*J/(LY+1))
C         IF(LZ.GT.1) TOTSUM=TOTSUM+COS(PI*K/(LZ+1))
      ENDIF
      TOTSUM=TOTSUM*2.0_dp*BHUB
      IF(TREAL) TOTSUM=0.0_dp
      RETURN
      END

      SUBROUTINE CALCTMATHUB(NBASIS,NBASISMAX,BHUB,TTILT,G1,TREAL,TPBC)
         use constants, only: dp
         use SystemData, only: BasisFN
         USE OneEInts, only: TMat2D,TMATSYM,SetupTMAT
         USE Parallel_neci, only : iProcIndex
         IMPLICIT NONE
         INTEGER NBASIS,nBasisMax(5,*)
         TYPE(BasisFN) G1(nBasis)
         real(dp) BHUB
         INTEGER iSize
         INTEGER I,J
         INTEGER DX,DY,DZ,LX,LY,LZ
         LOGICAL TTILT,TREAL,TPBC
         real(dp) TOTSUM
         character(*), parameter :: this_routine = 'CALCTMATHUB'
         IF(iProcIndex.eq.0) OPEN(10,FILE='TMAT',STATUS='UNKNOWN')
         CALL SetupTMAT(NBASIS,2,iSize)
         TOTSUM=0.0_dp
         LX=NBASISMAX(1,2)-NBASISMAX(1,1)
         LY=NBASISMAX(2,2)-NBASISMAX(2,1)
         LZ=NBASISMAX(3,2)-NBASISMAX(3,1)
         IF(LY.EQ.0) LY=-1
         IF(LZ.EQ.0) LZ=-1
         IF(TREAL) THEN
            DO I=1,NBASIS
               DO J=1,NBASIS
                  DX=ABS(G1(I)%k(1)-G1(J)%k(1))
                  DY=ABS(G1(I)%k(2)-G1(J)%k(2))
                  DZ=ABS(G1(I)%k(3)-G1(J)%k(3))
                  TOTSUM=BHUB
                  IF(TPBC) THEN
!This bit is only for if the hubbard lattice only has one site in a certain dimension
                     IF(DX.EQ.LX) THEN
                        DX=1
                        IF(LX.EQ.1) TOTSUM=TOTSUM+BHUB
                     ENDIF
                     IF(DY.EQ.LY) THEN
                        DY=1
                        IF(LY.EQ.1) TOTSUM=TOTSUM+BHUB
                     ENDIF
                     IF(DZ.EQ.LZ) THEN
                        DZ=1
                        IF(LZ.EQ.1) TOTSUM=TOTSUM+BHUB
                     ENDIF
                  ENDIF
                  IF(DX+DY+DZ.EQ.1.AND.G1(I)%Ms.EQ.G1(J)%MS) THEN
!This is for if the site can interact with a periodic image.
                    TMAT2D(I,J)=TOTSUM
                    IF(iProcIndex.eq.0) WRITE(10,*) I,J,TOTSUM
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            IF(TPBC) THEN
               DO I=1,NBASIS
                  CALL HUBKIN(G1(I)%k(1),G1(I)%k(2),G1(I)%k(3),
     &               NBASISMAX,BHUB,
     &               TTILT,TOTSUM,TREAL)
                  TMAT2D(I,I)=TOTSUM
                  IF(iProcIndex.eq.0) WRITE(10,*) I,I,TMAT2D(I,I)
               ENDDO
            ELSE
               DO I=1,NBASIS
                  CALL HUBKINN(G1(I)%k(1),G1(I)%k(2),G1(I)%k(3),
     &               NBASISMAX,BHUB,
     &               TTILT,TOTSUM,TREAL)
                  TMAT2D(I,I)=TOTSUM
                  IF(iProcIndex.eq.0) WRITE(10,*) I,I,TMAT2D(I,I)
               ENDDO
            ENDIF
         ENDIF
         IF(iProcIndex.eq.0) CLOSE(10)
         RETURN
      END
C.. NBASISMAX descriptor (1,3)
C
C HUBBARD:
C 0 Non-Tilted Lattice - pbc
C 1 Tilted Lattice - pbc
C 2 Non-Tilted lattice - no pbc
C 3 Tilted Lattice - no pbc
C.. four following are REAL
C 4 Non-Tilted Lattice - pbc
C 5 Tilted Lattice - pbc
C 6 Non-Tilted lattice - no pbc
C 7 Tilted Lattice - no pbc
C
      SUBROUTINE SETBASISLIM_HUB(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,
     &            TPBC,TREAL)
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),NMAXX,NMAXY,NMAXZ,LEN
         LOGICAL TPBC,TREAL
         IF(TPBC) THEN
            NBASISMAX(1,3)=0
         ELSE
C.. Non-tilted, not pbc
            NBASISMAX(1,3)=2
         ENDIF
         IF(TREAL) NBASISMAX(1,3)=NBASISMAX(1,3)+4
C         IF(.NOT.TPBC.AND..NOT.TREAL) THEN
C.. non-pbc has Huckel MOs starting k from 1
C            NBASISMAX(1,1)=1 
C            NBASISMAX(2,1)=1 
C            NBASISMAX(3,1)=1
C            NBASISMAX(1,2)=NMAXX
C            NBASISMAX(2,2)=NMAXY
C            NBASISMAX(3,2)=NMAXZ
C         ELSE
         IF(MOD(NMAXX,2).EQ.0) THEN
            NBASISMAX(1,2)=NMAXX/2
            NBASISMAX(1,1)=-NMAXX/2+1
         ELSE
            NBASISMAX(1,2)=NMAXX/2
            NBASISMAX(1,1)=-NMAXX/2
         ENDIF
         IF(MOD(NMAXY,2).EQ.0) THEN
            NBASISMAX(2,2)=NMAXY/2
            NBASISMAX(2,1)=-NMAXY/2+1
         ELSE
            NBASISMAX(2,2)=NMAXY/2
            NBASISMAX(2,1)=-NMAXY/2
         ENDIF
         IF(MOD(NMAXZ,2).EQ.0) THEN
            NBASISMAX(3,2)=NMAXZ/2
            NBASISMAX(3,1)=-NMAXZ/2+1
         ELSE
            NBASISMAX(3,2)=NMAXZ/2
            NBASISMAX(3,1)=-NMAXZ/2
         ENDIF
C         ENDIF
         NBASISMAX(1,4)=0
         NBASISMAX(2,4)=1
         NBASISMAX(1,5)=NMAXX
         NBASISMAX(2,5)=NMAXY
         NBASISMAX(3,5)=NMAXZ
         LEN=NMAXX*NMAXY*NMAXZ*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
      END

      SUBROUTINE SETBASISLIM_HUBTILT(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,
     &         TPBC,ITILTX,ITILTY)
         use constants, only : sizeof_int,dp
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),NMAXX,NMAXY,NMAXZ,LEN
         LOGICAL TPBC
         INTEGER ITILTX,ITILTY
         character(*), parameter :: this_routine = 'SETBASISLIM_HUBTILT'
         IF(TPBC) THEN
C.. Indicate tilted
            NBASISMAX(1,3)=1
         ELSE
C.. Not periodic boundaries
            NBASISMAX(1,3)=3
         ENDIF
         NBASISMAX(1,2)=int(real(NMAXX,dp)*                          
     &          (real(ITILTX+ITILTY,dp)/2.0_dp),sizeof_int)
         NBASISMAX(1,1)=-NBASISMAX(1,2)
C+1-MOD(NMAXX,2)
         IF(NMAXY.NE.NMAXX)
     &      call stop_all(this_routine,
     &                   'CANNOT HANDLE NON-SQUARE TILTED HUBBARD')
         NBASISMAX(2,2)=int(real(NMAXX,dp)*        
     &          (real(ITILTY+ITILTX,dp)/2.0_dp),sizeof_int)
         NBASISMAX(2,1)=-NBASISMAX(2,2)
C+1-MOD(NMAXX,2)
         IF(NMAXZ.GT.1) call stop_all(this_routine,
     &                          'CANNOT HANDLE TILTED 3D HUBBARD')
         NBASISMAX(3,2)=0
         NBASISMAX(3,1)=0
         NBASISMAX(1,4)=ITILTX
         NBASISMAX(2,4)=ITILTY
         NBASISMAX(1,5)=NMAXX
         NBASISMAX(2,5)=NMAXY
         NBASISMAX(3,5)=NMAXZ
!Len is number of basis functions
         LEN=NMAXX*NMAXY*(ITILTX*ITILTX+ITILTY*ITILTY)
     &         *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
C         LEN=NMAXX*NMAXY*2
C     &         *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
      END

!Generate the Sym table for a Hubbard Lattice
      Subroutine GenHubMomIrrepsSymTable(G1,nBasis,nBasisMax)
         use constants, only: int64
         use global_utilities
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nSym,SymConjTab,SymClasses,SymLabels
         use SymData, only: nSymLabels,tAbelian,SymTable
         use SymData, only: tagSymConjTab,tagSymClasses,tagSymLabels
         use SymData, only: tagSymTable
         use sym_mod, only: RoundSym, AddElecSym
         IMPLICIT NONE
         INTEGER nBasis,nBasisMax(*)
         TYPE(BasisFN) G1(nBasis)
         INTEGER I,J,nSyms,K
         TYPE(BasisFN) NQNS,S
         INTEGER iFindBasisFn
         character(*),parameter:: this_routine='GenHubMomIrrepsSymTable'
         nSym=nBasis/2
         WRITE(6,"(A,I3,A)") "Generating abelian symmetry table with",
     &    nBasis/2 , " generators for Hubbard momentum" 

C.. Now generate a list of sym labels.
         write (6,*) 'SIZES',nSymLabels, nBasis, allocated(symlabels),
     &     associated(symclasses), allocated(symconjtab), 
     &     allocated(symtable)
         if (allocated(SymLabels)) then
             write (6,'(a/a)') 
     &                'Warning: symmetry info already allocated.',
     &                'Deallocating and reallocating.'
             deallocate(SymLabels)
             call LogMemDealloc(this_routine,tagSymLabels)
         end if
         allocate(SymLabels(nSym))
         call LogMemAlloc('SymLabels',nSym,SymmetrySize,this_routine,
     &                                                     tagSymLabels)
         if (associated(SymClasses)) then
             deallocate(SymClasses)
             call LogMemDealloc(this_routine,tagSymClasses)
         end if
         allocate(SymClasses(nBasis))
         call LogMemAlloc('SymClasses',nBasis,4,this_routine,
     &                                                    tagSymClasses)
         NSYMLABELS=NSYM
         nSyms=1
         tAbelian=.false.
         DO I=1,NBASIS/2
C.. place the sym label of each state in SymClasses(ISTATE).
            IF(ALL(G1(I*2)%k.EQ.0)) THEN
               SymClasses(I)=1
               NQNS=G1(I*2)
            ELSE
               nSyms=nSyms+1
               SymClasses(I)=nSyms
            ENDIF
C.. list the symmetry string of each sym label
C..            SymLabels(I)%s=2_int64**(I-1)
            SymLabels(i)%s = i
         ENDDO
C.. Setup the symmetry product table
         if (allocated(SymTable)) then
             deallocate(SymTable)
             call LogMemDealloc(this_routine,tagSymTable)
         end if
         allocate(SymTable(nSym,nSym))
         call LogMemAlloc('SymTable',nSym**2,SymmetrySize,this_routine,
     &                                                      tagSymTable)
         if (allocated(SymConjTab)) then
             deallocate(SymConjTab)
             call LogMemDealloc(this_routine,tagSymConjTab)
         end if
         allocate(SymConjTab(nSym))
         call LogMemAlloc('SymConjTable',nSym,4,this_routine,
     &                                                    tagSymConjTab)
         SYMTABLE(1:NSYM,1:NSYM)=Symmetry(0)
         DO I=1,nBasis/2
            NQNS%k=-G1(I*2)%k
            CALL RoundSym(NQNS,nBasisMax)
            J=iFindBasisFn(NQNS,G1,nBasis)
            IF(J.EQ.0) THEN 
               WRITE(6,*) "Cannot Find symmetry conjugate to basis fn ",
     &            I*2
               call stop_all(this_routine,
     &                     "Cannot find symmetry conjugate.")
            ENDIF
            SymConjTab(SymClasses(I))=SymClasses(J/2)
            DO J=1,nBasis/2
               S=G1(I*2)
               CALL AddElecSym(J*2,G1,nBasisMax,S)
               CALL RoundSym(S,nBasisMax)
               NQNs%k=S%k
               K=iFindBasisFn(NQNS,G1,nBasis)
               IF(K.EQ.0) THEN 
                WRITE(6,*) "Cannot find symmetry product of basis fns ",
     &               I*2, J*2
                  call stop_all(this_routine,
     &                        "Cannot find symmetry product.")
               ENDIF
               SymTable(SymClasses(I),SymClasses(J))=
     &            SymLabels(SymClasses(K/2))
            ENDDO
         ENDDO
         DO I=1,nBasis/2
            G1(I*2-1)%Sym=SymLabels(SymClasses(I))
            G1(I*2)%Sym=SymLabels(SymClasses(I))
         ENDDO
         WRITE(6,*) "Symmetry, Symmetry Conjugate"
         DO I=1,NSYM
            WRITE(6,*) I,SymConjTab(I)
         ENDDO
      End

!Hubbard Sym Reps are different from normal ones.  In hubbard, we have split degenerate sets into 1D subcomponents,
!  but a complete degenerate set will count as a sym rep (i.e. if it is filled, its sym can be discounted),
!  not just one of the subcomponents
      SUBROUTINE GENHUBSYMREPS(NBASIS,ARR,BRR)
         use global_utilities
         use constants, only: dp
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymReps,tagSymReps
         IMPLICIT NONE
         INTEGER I,J
         INTEGER NBASIS,BRR(NBASIS)
         real(dp) ARR(NBASIS)
         character(*), parameter :: this_routine='GenHubSymReps'
         
C.. now work out which reps are degenerate and label them
         allocate(SymReps(2,nBasis))
         call LogMemAlloc('SymReps',2*nBasis,4,this_routine,tagSymReps)
         J=0
         
         DO I=1,NBASIS
            IF(I.EQ.1) THEN
C.. we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ELSE
                IF(ABS(ARR(I)-ARR(I-1)).LT.1.0e-5_dp) THEN
C.. we have the same degenerate rep as the previous entry
                   SYMREPS(2,J)=SYMREPS(2,J)+1
                ELSE
C.. we have a new rep
                   J=J+1
                   SYMREPS(2,J)=1
                ENDIF
            ENDIF
            SYMREPS(1,BRR(I))=J
         ENDDO
!         DO I=1,NBASIS
!            WRITE(6,*) SYMREPS(1,I),SYMREPS(2,I)
!         ENDDO
      END

      LOGICAL FUNCTION KALLOWED(G,NBASISMAX)
        ! See if a given G vector is within a (possibly tilted) unit cell.
        ! Used to generate the basis functions for the hubbard model (or perhaps electrons in boxes)
        use constants, only: dp
        use SystemData, only: BasisFN
        IMPLICIT NONE
        Type(BasisFN), intent(in) :: G
        INTEGER nBasisMax(5,*),NMAXX,I,J,AX,AY
        INTEGER KX,KY
        real(dp) MX,MY,XX,YY
        LOGICAL TALLOW
        TALLOW=.TRUE.
        IF(NBASISMAX(3,3).EQ.1) THEN
      !.. spatial symmetries
            IF(G%k(1).NE.0) TALLOW=.FALSE.
        ELSEIF(NBASISMAX(3,3).EQ.0) THEN
      !.. Hubbard
            IF(NBASISMAX(1,3).EQ.1) THEN
      !.. Tilted hubbard
                NMAXX=NBASISMAX(1,5)
      !         NMAXY=
                AX=NBASISMAX(1,4)
                AY=NBASISMAX(2,4)
      !.. (XX,YY) is the position of the bottom right corner of the unit cell
                XX=((AX+AY)/2.0_dp)*NMAXX
                YY=((AY-AX)/2.0_dp)*NMAXX
                MX=XX*AX+YY*AY
                MY=XX*AY-YY*AX
                I=G%k(1)
                J=G%k(2)
                KX=I*AX+J*AY
                KY=I*AY-J*AX
                IF(KX.GT.MX) TALLOW=.FALSE.
                IF(KY.GT.MY) TALLOW=.FALSE.
                IF(KX.LE.-MX) TALLOW=.FALSE.
                IF(KY.LE.-MY) TALLOW=.FALSE.
            ELSEIF(NBASISMAX(1,3).GE.4.OR.NBASISMAX(1,3).EQ.2) THEN
      !.. Real space Hubbard
                IF(G%k(1).EQ.0.AND.G%k(2).EQ.0.AND.G%k(3).EQ.0) THEN
                   TALLOW=.TRUE.
                ELSE
                    TALLOW=.FALSE.
                ENDIF
      !      ELSEIF(NBASISMAX(1,3).EQ.2) THEN
      !.. mom space non-pbc non-tilt hub - parity sym
      !          IF(  (G(1).EQ.0.OR.G(1).EQ.1)
      !     &       .AND.(G(2).EQ.0.OR.G(2).EQ.1)
      !     &       .AND.(G(3).EQ.0.OR.G(3).EQ.1)) THEN
      !              TALLOW=.TRUE.
      !          ELSE
      !              TALLOW=.FALSE.
      !          ENDIF
      !      ELSEIF(NBASISMAX(1,3).EQ.2) THEN
      !.. non-pbc hubbard
      !          TALLOW=.TRUE.
      !          IF(G(1).GT.NBASISMAX(1,2).OR.G(1).LT.NBASISMAX(1,1))
      !     &       TALLOW=.FALSE.
      !          IF(G(2).GT.NBASISMAX(2,2).OR.G(2).LT.NBASISMAX(2,1))
      !     &       TALLOW=.FALSE.
      !          IF(G(3).GT.NBASISMAX(3,2).OR.G(3).LT.NBASISMAX(3,1))
      !     &       TALLOW=.FALSE.
            ELSE
      !.. Normal Hubbard
                TALLOW=.TRUE.
                IF(G%k(1).GT.NBASISMAX(1,2).OR.G%k(1).LT.NBASISMAX(1,1))
     &               TALLOW=.FALSE.
                IF(G%k(2).GT.NBASISMAX(2,2).OR.G%k(2).LT.NBASISMAX(2,1))
     &               TALLOW=.FALSE.
                IF(G%k(3).GT.NBASISMAX(3,2).OR.G%k(3).LT.NBASISMAX(3,1))
     &              TALLOW=.FALSE.
            ENDIF
        ENDIF         
        KALLOWED=TALLOW
        RETURN
      END FUNCTION KALLOWED

      INTEGER FUNCTION IFINDBASISFN(NQNS,G1,NBASIS)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER NBASIS,I,J
         TYPE(BasisFN) G1(NBASIS),NQNS
         LOGICAL L
         DO I=1,NBASIS
            L=.TRUE.
            DO J=1,3
               IF(G1(I)%k(J).NE.NQNS%k(J)) L=.FALSE.
            ENDDO
            IF(G1(I)%Sym%s.NE.NQNS%Sym%s) L=.FALSE.
            IF(G1(I)%Ms.NE.NQNS%Ms) L=.FALSE.
            IF(L) THEN
               IFINDBASISFN=I
               RETURN
            ENDIF
         ENDDO
         IFINDBASISFN=0
         RETURN
      END FUNCTION IFINDBASISFN
