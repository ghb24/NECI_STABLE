      SUBROUTINE CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER NBASISMAX(5,2), BRR(LEN),LEN
         TYPE(BASISFN) G1(LEN)
         REAL*8 ARR(LEN,2)
         INTEGER I,J
         TYPE(Symmetry) iDecomp
         LOGICAL IMPROP(NKPS)
         LOGICAL ALL1D
         integer*8  ComposeAbelianSym
C         NBASISMAX(1,1)=1
C         NBASISMAX(1,2)=NSTATES
C         NBASISMAX(1,3)=2
C         NBASISMAX(4,1)=-1
C         NBASISMAX(4,2)=1
         NBASISMAX(1,1)=0
         NBASISMAX(2,1)=0
         NBASISMAX(3,1)=0
         IF(BIGCELL.AND.PROJECT) THEN
            NBASISMAX(1,2)=CELLSIZES(1)-1
            NBASISMAX(2,2)=CELLSIZES(2)-1
            NBASISMAX(3,2)=CELLSIZES(3)-1
         ELSE
            NBASISMAX(1,2)=0
            NBASISMAX(2,2)=0
            NBASISMAX(3,2)=0
         ENDIF
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1

C.. set ISPINSKIP=0, to tell the SCRs that there's no UMAT
         NBASISMAX(2,3)=0
         IF(TKP) THEN
!  We've got a k-point code, so we're not currently using spatial symmetry.  We use translational symmetry according to the k-points
!            JSS: Use Abelian symmetry formulation.
             CALL GenKPtIrreps(KPNT_CHAR,NKPTRANS,NKPS,KPNTIND,NSTATES)
         ELSE
             CALL ExtractSymLabels
             CALL GENIRREPS(TKP,IMPROPER_OP,NROTOP)
             CALL GENSYMTABLE
         END IF

C         CALL FAKEEIGENVALUES

         CALL IAZZERO(G1,LEN*BasisFNSize)
         ALL1D=.TRUE.
         DO I=1,NSTATES
           if (TAbelian) then
               IDECOMP%s=ComposeAbelianSym(KpntSym(:,KPntInd(I)))
           else
               CALL DECOMPOSEREP(SYMLABELCHARS(1,SymClasses(I)),IDECOMP)
           end if
            IF(ABS(SYMLABELCHARS(1,SymClasses(I))-1.D0).GE.1.D-6)
     &          ALL1D=.FALSE.
            G1(I*2-1)%SYM=IDECOMP
            G1(I*2)%SYM=IDECOMP
            G1(I*2-1)%MS=-1
            G1(I*2)%MS=1
            IF(BIGCELL.AND.PROJECT) THEN
               DO J=1,3
                  G1(I*2-1)%K(J)=K_VECTORS(J,TRANS_LABEL(I))
                  G1(I*2)%K(J)=K_VECTORS(J,TRANS_LABEL(I))
               ENDDO
            ELSE
               DO J=1,3
                  G1(I*2-1)%K(J)=0
                  G1(I*2)%K(J)=0
               ENDDO
            ENDIF
C.. Fake the eigenvalues
            ARR(2*I-1,1)=EIGENVALUES(I)
            ARR(2*I,1)=EIGENVALUES(I)
            ARR(2*I-1,2)=EIGENVALUES(I)
            ARR(2*I,2)=EIGENVALUES(I)
            BRR(2*I-1)=2*I-1
            BRR(2*I)=2*I
         ENDDO
         IF (TAbelian) THEN
            WRITE (6,*) 'Using Abelian symmetry formulation.'
            NBASISMAX(5,1)=0
            NBASISMAX(5,2)=NSYM-1
         ELSE IF (ALL1D) THEN
C.. If all reps are 1D, we can block determinants in different symmetries
            NBASISMAX(5,1)=0
            NBASISMAX(5,2)=NSYM-1
            WRITE(6,*)
     &   "All orbitals have 1D symmetry. Using blocking symmetry."
         ELSE
            WRITE(6,*)
     &   "Multidimensional symmetries detected. ",
     &   "Not using symmetry blocking."
         ENDIF
C.. show it's a generic spatial basis
         NBASISMAX(3,3)=1
      END
C.. Called to retrieve number of basis functions from CPMD
      SUBROUTINE CPMDSYSTEMINIT(LEN)
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER LEN
C         CALL FAKEREADCHARTABLE(NSTATES)
         LEN=NSTATES*2
         IF(TSPACEGP) THEN
C.. We're dealing with a complete space group rather than just a point
C.. group
            NROT=NSpGpOp
         ELSE
            NROT=NROTOP
         ENDIF
         IF(TKP) THEN
            NROT=NKPTRANS
         ENDIF
      END
      
      SUBROUTINE GENCPMDSYMREPS(G1,NBASIS,ARR,BRR)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INTEGER BRR(NBASIS),NBASIS,G1(*)
         REAL*8 ARR(NBASIS,2)
         CALL GENSYMREPS(G1,NBASIS,ARR,BRR,DEGENTOL)
      END

      SUBROUTINE EXTRACTSYMLABELS
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER I,J,K
         INTEGER LABELS(2,NSTATES)
C         IF(TSPACEGP) THEN
C            NSYMLABELS=NSTATES
C         ELSE
            NSYMLABELS=0
            DO I=1,NSTATES
               IF(ROT_LABEL(I).GT.NSYMLABELS) NSYMLABELS=ROT_LABEL(I)
            ENDDO
C         ENDIF
         CALL MEMORY(IP_SYMLABELCHARS,NSYMLABELS*NROT*2,'SYMLABELCH')
         CALL AZZERO(SYMLABELCHARS,NSYMLABELS*NROT*2)
         CALL MEMORY(IP_SYMLABELS,nSymLabels,'SYMLABELS')
         CALL MEMORY(IP_SymClasses,NSTATES,'SymClasses')
C.. Two versions - one for space groups and one for point groups
         IF(TSPACEGP) THEN
C.. space group is more complex as we have separate translational and
C.. rotational labels.  We will need to compare each pair of t&r labels
C.. with a stored list.
            WRITE(6,*) "Using space group representation."
C            NSYMLABELS=0
C            DO I=1,NSTATES
C               DO K=1,NSYMLABELS
C                  IF(LABELS(1,K).EQ.ROT_LABEL(I)
C     &               .AND.LABELS(2,K).EQ.TRANS_LABEL(I)) THEN
CC.. We can re-use the old label
C                     SymClasses(I)=K
C                     EXIT
C                  ENDIF
C               END DO
C               IF(K.GT.NSYMLABELS) THEN
C                  NSYMLABELS=K
CC.. we've found a new symlabel - fill it
C                  DO J=1,NROT
CC.. Indices need to be inverted so cannot straight ICOPY
C                     SYMLABELCHARS(J,K)=GROUP_CHAR(I,J)
C                  ENDDO
C                  SymClasses(I)=K
C                  LABELS(1,K)=ROT_LABEL(I)
C                  LABELS(2,K)=TRANS_LABEL(I)
C               ENDIF
C            ENDDO
            DO I=1,NSTATES
               IF(SYMLABELCHARS(1,ROT_LABEL(I)).EQ.(0.D0,0.D0)) THEN
C.. we've found an uninitialized symlabel - fill it
                  DO J=1,NROT
C.. Indices need to be inverted so cannot straight ICOPY
                     SYMLABELCHARS(J,ROT_LABEL(I))=GROUP_CHAR(J,I)
                  ENDDO
               ENDIF
               SymClasses(I)=ROT_LABEL(I)
            ENDDO
         ELSE
C.. point group
            WRITE(6,*) "Using point group representation."
            DO I=1,NSTATES
               IF(SYMLABELCHARS(1,ROT_LABEL(I)).EQ.(0.D0,0.D0)) THEN
C.. we've found an uninitialized symlabel - fill it
                  DO J=1,NROT
C.. Indices need to be inverted so cannot straight ICOPY
                     SYMLABELCHARS(J,ROT_LABEL(I))=NINT(ROT_CHAR(I,J))
                  ENDDO
               ENDIF
               SymClasses(I)=ROT_LABEL(I)
            ENDDO
         ENDIF
         WRITE(6,*) "SYMMETRY CLASSES"
         CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
C.. Allocate memory gor irreps.
C.. Assume there will be no more than 64 irreps
         CALL MEMORY(IP_IRREPCHARS,NROT*64*2,"IRREPCH")
      END

      SUBROUTINE FAKEEIGENVALUES
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INTEGER I
         CALL MEMORY(IP_EIGENVALUES,NSTATES,"EIGENVALUES")
         DO I=1,NSTATES
            EIGENVALUES(I)=I
         ENDDO
      END 
      SUBROUTINE FAKEREADCHARTABLE(NSTATES)
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INCLUDE 'sym.inc'
C.. First set up the labels and chars
         INTEGER CHARS(16,7)
         DATA ((CHARS(J,I),J=1,16),I=1,7)
     &               /1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &               2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2,
     &               2, 0, 0,-2, 0, 0, 0, 0,-2, 0, 0, 2, 0, 0, 0, 0,
     &               1,-1,-1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1, 1,-1,-1,
     &               4, 0, 0,-4, 0, 0, 0, 0,-4, 0, 0, 4, 0, 0, 0, 0,
     &               1,-1,-1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,
     &               2,-2,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2/
         INTEGER LABS(19)
         DATA LABS /1,2,2,1,3,3,4,2,2,1,5,5,5,5,
     &              3,3,6,7,7/
         INTEGER NSTATES,I,J
         INTEGER ARRAY(*)
         POINTER (IP_ARRAY,ARRAY)
         NSYMLABELS=7
         NROT=16
         NSTATES=19
         CALL MEMORY(IP_SYMLABELS,NSYMLABELS*SymmetrySize,"SYMLABELS")
         CALL MEMORY(IP_SymClasses,NSTATES,"SymClasses")
         CALL MEMORY(IP_SYMLABELCHARS,NROT*NSYMLABELS*2,"SYMLABELCH")
         IP_ARRAY=IP_SYMLABELCHARS
         CALL MEMORY(IP_IRREPCHARS,NROT*NSYMLABELS*2,"IRREPCH")
         DO I=1,7
            DO J=1,16
               ARRAY((I-1)*16+J)=CHARS(J,I)
            ENDDO
         ENDDO
         DO I=1,NSTATES
            SymClasses(I)=LABS(I)
         ENDDO
         WRITE(6,*) "SYMMETRY CLASSES", IP_SYMLABELCHARS,
     &      LOC(SYMLABELCHARS)
         CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
      END

C.. Calculate the h_ij elements for CPMD orbitals.
C.. E_i is the eigenvalue for state i.
C.. f_i is the occupation number for state i.
C.. delta_ij is the Dirac delta.
C.. (ab|cd) is the 4-index Coulomb integral between states, i,j,k,l
C.. h_ij=(E_i+xi) delta_ij - Sum_k f_k (ik|jk)
      SUBROUTINE CPMDINIT2INDINT(NHG,NORBUSED,NBASISMAX,ISS,G1,NEL,
     &      ECORE,TORDER,ARR,BRR,iCacheFlag)
         USE HElem
         use System, only: BasisFN
         USE UMatCache, only: nstatesumat=>nstates,tUMAT2D,
     &                        TMat2D,tCPMDSymTMat,GetUMatEl,TMatSym
         USE UMatCache, only: SetUmatTrans,SetUMatCacheFlag,SetupUMat2D,
     &                        SetupUMatCache,TMatInd    
         IMPLICIT NONE
         INTEGER NHG,NTFROZEN
         INCLUDE 'cpmddata.inc'
         TYPE(HElement) CPMD1EINT(NSTATES,NSTATES)
         INTEGER I,J,K,NSTATESUSED,ISUB,NORBUSED,II,JJ
         REAL*8 F,ALAT(3)  !ALAT just dummy for GetUMatEl.
         TYPE(HElement) U,TOT,sum1
         INTEGER ISS,NBASISMAX(*)
         TYPE(HElement) UMAT(0)  ! UMAT dummy for GetUMatEl.

         TYPE(BASISFN) G1(NHG)
         REAL*8 ECORE
         INTEGER NEL
         LOGICAL LCHKSYM
         REAL*8 DIAGTMAT(NSTATES)
         INTEGER STATEORDER(NSTATES)
         LOGICAL TORDER
         INTEGER BRR(NSTATES*2),NLOCK
         REAL*8 ARR(NSTATES*2,2)
         complex*16,allocatable,save :: HarInt(:,:),FockMat(:,:)
         INTEGER iCacheFlag
         ! For ZHEEV (eigenvalues of FockMatrix).
         integer :: info
         real*8,allocatable :: rwork(:),FockEigv(:)
         complex*16,allocatable :: work(:)
         CALL TISET('CPMDINIT2I',ISUB)

         ECORE=EIONION
         WRITE(6,*) "Core Energy:", ECORE 
         NLOCK=(NEL+1)/2
         IF(TORDER) THEN
C.. We have to cache all states, as we don't, as yet, know which ones we're going to freeze because we have to reorder.
            NSTATESUSED=NSTATES
         ELSE
            NSTATESUSED=NORBUSED/2
         ENDIF
         IF(.NOT.BTEST(iCacheFlag,1)) THEN
!Don't initialize the cache if it's already there
            allocate(HarInt(NStatesUsed,NStatesUsed))
!            allocate(FockMat(nStatesUsed,nStatesUsed))
            call Memory(IP_PIInt,nStates,'PIInt')
C.. If we're freezing we allocate a small cache first to store the <ij|kj> integrals, which will be deallocated and a larger cache made later.
            CALL SETUPUMATCACHE(NSTATESUSED,NSTATESUSED.NE.NSTATES)
         ENDIF
         CALL CPMDGEN1EINTS(NStatesUsed,G1,CPMD1EINT)
         CALL MEMORY_CHECK()
         OPEN(10,FILE='TMAT',STATUS='UNKNOWN')
         IF(.NOT.BTEST(iCacheFlag,1)) THEN
!Don't initialize the cache if it's already there
C.. Set the UMAT cache to cache things in the lowest possible position
            CALL SETUMATCACHEFLAG(1)
C   JSS  Calculate the elements of the type <ij|ij> and <ij|ji> efficiently.
            CALL SETUPUMAT2D(G1,HarInt)
            ! Must now calculate the hartree integrals and periodic integrals 
            ! (otherwise done in CPMDAntiSymIntEl) if not using UMAT2D.
            if(.not.TUMAT2D) call CalcHarPIInts(HarInt,NStatesUsed)
            CALL SETUMATCACHEFLAG(0)
         ENDIF

! JSS. 25/01/08.  Calculate the Fock matrix.
! F_ij = < i | -1/2 \nabla^2 + v_ext | j > + \sum_k [ <ik|jk> - <ik|kj> ]
!      = - h^{KS}_ij - <i|v_xc|j> - \sum_k <ik|kj>
!      = \epsilon_i \delta_ij - \sum_k <ik|kj>
! where h^{KS}_ij = < i | -1/2 \nabla^2 + v_ext + v_har + v_xc | j >
! and h^{KS}_ii = \epsilon_i
!         FockMat(:,:)=dcmplx(0.d0,0.d0)

         WRITE(6,*) "Calculating TMAT"
         DO I=1,NSTATESUSED
            II=I*2-1
            DO J=I,NSTATESUSED
               JJ=J*2-1

               IF(LCHKSYM(G1(II),G1(JJ))) THEN

                  TOT=0.D0
                  IF(I.EQ.J) THEN
                      TOT=TOT+HElement(EIGENVALUES(I)+PIInt(I)/2.d0)
!                      FockMat(I,J)=FockMat(I,J)+Eigenvalues(I)
                  END IF
                  TOT=TOT-CPMD1EINT(I,J)
                  TOT=TOT-HElement(HarInt(i,j))
!                  FockMat(i,j)=FockMat(i,j)-dcmplx(CPMD1EInt(I,J))
                  IF (tCPMDSymTMat) THEN
                      TMATSYM(TMatInd(II+1,JJ+1))=TOT ! Indexing needs to be compatiable with TMatInd.
                  ELSE
                      TMAT2D(II,JJ)=TOT
                      TMAT2D(II+1,JJ+1)=TOT
                      TMAT2D(JJ,II)=dconjg(TOT)
                      TMAT2D(JJ+1,II+1)=dconjg(TOT)
                  ENDIF
                  WRITE(10,*) I,J,TOT
                  IF(I.EQ.J) THEN
                     DIAGTMAT(I)=SQRT(SQ(TOT))
                     STATEORDER(I)=I
                  ENDIF

                  ! Subtract exchange integrals from Fock matrix.
!                 do K=1,nStates
!                     call CPMDGetOcc(K,F) ! F is the occupation number of state K.
!                      if (F.gt.1.d-5) then
!                          ! Get exchange integral <ik|kj>
!                          if (K.gt.nStatesUsed) then
!                              write (6,"(a,i4,a,f10.6)") 
!     &                             "Top frozen state ",K,
!     &                             " has occupation number",F
!                              stop "Top Frozen states occupied."
!                          end if
!                          U=GetUMatEl(nBasisMax,UMat,ALAT,NHG,ISS,G1,
!     &                                I,K,K,J)
!                          FockMat(i,j)=FockMat(i,j)-F*dcmplx(U)/2
!                      end if
!                  end do
!                 FockMat(j,i)=dconjg(FockMat(i,j))

               ENDIF

            END DO
         END DO

C.. Now we order it if we need to
         IF(TORDER) THEN
            WRITE(6,*) "Re-ordering CPMD orbitals according ",
     &        "to one-electron energies."
           CALL SORT2(NSTATESUSED-NLOCK,DIAGTMAT(NLOCK+1),
     &         STATEORDER(NLOCK+1))
C.. Now copy to BRR
            DO I=1,NSTATES
               BRR(2*I-1)=STATEORDER(I)*2-1
               BRR(2*I)=STATEORDER(I)*2
               ARR(2*I-1,1)=DIAGTMAT(I)
               ARR(2*I,1)=DIAGTMAT(I)
            ENDDO
C.. we only need a translation table if we've reordered the states.
C.. This should save some time in the UMAT lookup.
            CALL SETUMATTRANS(STATEORDER)
            CALL WRITEBASIS(6,G1,NSTATES*2,ARR,BRR)
         ENDIF
C.. Set the UMAT cache to cache normally again
         WRITE(6,*) "Finished TMAT"
         CLOSE(10)

!        Obtain the fock eigenvalues (based upon the KS orbitals).
!         allocate(FockEigv(nStatesUsed))
!         allocate(work(2*nStatesUsed))
!         allocate(rwork(3*nStatesUsed-2))
!         call zheev('N','U',nStatesUsed,FockMat,nStatesUsed,FockEigv,
!     &               work,2*nStatesUsed,rwork,info)
!         if (info.eq.0) then
!             ! Success! Output eigenvalues.
!             write (6,'(/a)') 'i  KS eigenvalues   Fock eigenvalues'
!             do i=1,nStatesUsed
!                 write (6,'(i4,2f12.8)') i,Eigenvalues(i),FockEigv(i)
!             end do
!             write (6,'(/)')
!
!         else
!             write (6,'(a,i)') 'Failure in ZHEEV. INFO=',INFO
!             write (6,'(a)') 'Fock eigenvalues not computed.'
!         end if
!         deallocate(FockMat)
!         deallocate(FockEigv)
!         deallocate(work)
!         deallocate(rwork)

         IF(.not.BTEST(iCacheFlag,0)) THEN
            deallocate(HarInt)
            call FreeM(IP_PIInt)
         endif
         CALL TIHALT('CPMDINIT2I',ISUB)
      END   
