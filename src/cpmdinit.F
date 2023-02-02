      SUBROUTINE CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
         use constants, only: dp
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,
     &                         NullBasisFn
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: KPntSym,tAbelian,nSym,SymClasses,
     &                      SymLabelChars
         use CPMDData, only: EigenValues,Project,BigCell,Trans_Label
         use CPMDData, only: K_Vectors,KPntInd,CellSizes,TKP,Improper_OP
         use CPMDData, only: nStates,NKPTrans,nKPS,nRotOp
         use sym_mod, only: ComposeAbelianSym, GenKPtIrreps
         use sym_mod, only: DECOMPOSEREP, GENIRREPS, GENSYMTABLE
         use util_mod, only: near_zero
         IMPLICIT NONE
         INTEGER LEN,nBasisMax(5,*), BRR(LEN)
         TYPE(BASISFN) G1(LEN)
         real(dp) ARR(LEN,2)
         INTEGER I,J
         TYPE(Symmetry) iDecomp
         LOGICAL ALL1D
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
             CALL GenKPtIrreps(NKPTRANS,NKPS,KPNTIND,NSTATES)
         ELSE
             CALL ExtractSymLabels
             CALL GENIRREPS(TKP,IMPROPER_OP,NROTOP)
             CALL GENSYMTABLE
         END IF

         G1(1:LEN)=NullBasisFn
         ALL1D=.TRUE.
         DO I=1,NSTATES
           if (TAbelian) then
               IDECOMP%s=ComposeAbelianSym(KpntSym(:,KPntInd(I)))
           else
               CALL DECOMPOSEREP(SYMLABELCHARS(1,SymClasses(I)),IDECOMP)
           end if
            IF(ABS(SYMLABELCHARS(1,SymClasses(I))-1.0_dp).GE.1.0e-6_dp)
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
C.. the eigenvalues
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
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nRot
         use CPMDData, only:tSpaceGp,nSpGpOp,nRotOp,nKPTrans,tKP,nStates
         IMPLICIT NONE
         INTEGER LEN
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

      SUBROUTINE GENCPMDSYMREPS(G1,NBASIS,ARR)
         use constants, only: dp
         use CPMDData, only: DegenTol
         use SystemData, only: BasisFN
         use sym_mod, only: GENSYMREPS
         IMPLICIT NONE
         INTEGER NBASIS
         type(BasisFn) G1(*)
         real(dp) ARR(NBASIS,2)
         CALL GENSYMREPS(G1,NBASIS,ARR,DEGENTOL)
      END

      SUBROUTINE EXTRACTSYMLABELS
         use global_utilities
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nSymLabels,SymLabelChars,nRot,SymClasses
         use SymData, only: SymLabels,SymLabelChars,SymClasses
         use SymData, only: IRREPCHARS
         use SymData, only: tagSymLabels,tagSymLabelChars,tagSymClasses
         use SymData, only: tagIRREPCHARS
         use CPMDData, only: nStates,Rot_Label,tSpaceGp,Group_Char
         use CPMDData, only: Rot_Char
         use sym_mod, only: WriteIrrepTab
         use util_mod, only: near_zero
         use constants
         IMPLICIT NONE
         INTEGER I,J
         character(*), parameter :: this_routine='ExtractSymLabels'
         NSYMLABELS=0
         DO I=1,NSTATES
           IF(ROT_LABEL(I).GT.NSYMLABELS) NSYMLABELS=ROT_LABEL(I)
         ENDDO
         allocate(SymLabelChars(nRot,nSymLabels))
         call LogMemAlloc('SymLabelChars',nSymLabels*nRot*2,16,
     &                                    this_routine,tagSymLabelChars)
         SYMLABELCHARS=0.0_dp
         allocate(SymLabels(nSymLabels))
         call LogMemAlloc('SymLabels',nSymLabels,4,this_routine,
     &                                                     tagSymLabels)
         allocate(SymClasses(nStates))
         call LogMemAlloc('SymClasses',nStates,4,this_routine,
     &                                                    tagSymClasses)
C.. Two versions - one for space groups and one for point groups
         IF(TSPACEGP) THEN
C.. space group is more complex as we have separate translational and
C.. rotational labels.  We will need to compare each pair of t&r labels
C.. with a stored list.
            WRITE(6,*) "Using space group representation."
            DO I=1,NSTATES
               IF(near_zero(SYMLABELCHARS(1,ROT_LABEL(I)))) THEN
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
               IF(near_zero(SYMLABELCHARS(1,ROT_LABEL(I)))) THEN
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
         allocate(IRREPCHARS(NROT,NROT)) ! nSym==nRot
         call LogMemAlloc('IRREPCHARS',nROT**2,16,this_routine,
     &                                                    tagIRREPCHARS)
         !CALL N_MEMORY(IP_IRREPCHARS,NROT*64*2,"IRREPCH")
      END

C.. Calculate the h_ij elements for CPMD orbitals.
C.. E_i is the eigenvalue for state i.
C.. f_i is the occupation number for state i.
C.. delta_ij is the Dirac delta.
C.. (ab|cd) is the 4-index Coulomb integral between states, i,j,k,l
C.. h_ij=(E_i+xi) delta_ij - Sum_k f_k (ik|jk)
      SUBROUTINE CPMDINIT2INDINT(NHG,NORBUSED,G1,NEL,
     &      ECORE,TORDER,ARR,BRR,iCacheFlag)
         use constants, only: dp
         use CPMDData, only: nStates,PIInt,EIonIon,EigenValues,tagPIInt
         use SystemData, only: BasisFN
         USE UMatCache, only: nstatesumat=>nstates,tUMAT2D
         USE UMatCache, only: SetUmatTrans,SetUMatCacheFlag,SetupUMat2D,
     &                        SetupUMatCache
         use OneEInts, only: TMat2D,tCPMDSymTMat,TMatSym
         use OneEInts, only: TMatInd
         use global_utilities
         use sort_mod
         use sym_mod, only: LCHKSYM
         use MemoryManager, only: TagIntType
         use util_mod, only: stop_all
         IMPLICIT NONE
         INTEGER NHG
         HElement_t(dp) CPMD1EINT(NSTATES,NSTATES)
         INTEGER I,J,NSTATESUSED,NORBUSED,II,JJ
         type(timer), save :: proc_timer
         HElement_t(dp) TOT

         TYPE(BASISFN) G1(NHG)
         real(dp) ECORE
         INTEGER NEL
         real(dp) DIAGTMAT(NSTATES)
         INTEGER STATEORDER(NSTATES)
         LOGICAL TORDER
         INTEGER BRR(NSTATES*2),NLOCK
         real(dp) ARR(NSTATES*2,2)
         complex(dp),allocatable,save :: HarInt(:,:)
         INTEGER iCacheFlag
         ! For ZHEEV (eigenvalues of FockMatrix).
         integer(TagIntType) :: tagHarInt=0
         character(*), parameter :: this_routine='CPMDINIT2INDINT'
         proc_timer%timer_name='CPMDINIT2I'
         call set_timer(proc_timer)

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
            call LogMemAlloc('HarInt',nStatesUsed**2,16,this_routine,
     &                                                        tagHarInt)
            allocate(PIInt(nStates))
            call LogMemAlloc('PIInt',nStates,8,this_routine,tagPIInt)
!            allocate(FockMat(nStatesUsed,nStatesUsed))
C.. If we're freezing we allocate a small cache first to store the <ij|kj> integrals, which will be deallocated and a larger cache made later.
            CALL SETUPUMATCACHE(NSTATESUSED,NSTATESUSED.NE.NSTATES)
         ENDIF
         call stop_all(this_routine,'Code deprecated')
!         CALL N_MEMORY_CHECK()
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
!         FockMat(:,:)=dcmplx(0.0_dp,0.0_dp)

         WRITE(6,*) "Calculating TMAT"
         DO I=1,NSTATESUSED
            II=I*2-1
            DO J=I,NSTATESUSED
               JJ=J*2-1

               IF(LCHKSYM(G1(II),G1(JJ))) THEN

                  TOT=0.0_dp
                  IF(I.EQ.J) THEN
                      TOT=TOT+(EIGENVALUES(I)+PIInt(I)/2.0_dp)
!                      FockMat(I,J)=FockMat(I,J)+Eigenvalues(I)
                  END IF
                  TOT=TOT-CPMD1EINT(I,J)
                  TOT=TOT-real(HarInt(i,j),dp)
!                  FockMat(i,j)=FockMat(i,j)-dcmplx(CPMD1EInt(I,J))
                  IF (tCPMDSymTMat) THEN
                      TMATSYM(TMatInd(II+1,JJ+1))=TOT ! Indexing needs to be compatiable with TMatInd.
                  ELSE
                      TMAT2D(II,JJ)=TOT
                      TMAT2D(II+1,JJ+1)=TOT
#ifdef CMPLX_
                      TMAT2D(JJ,II)=conjg(TOT)
                      TMAT2D(JJ+1,II+1)=conjg(TOT)
#else
                      TMAT2D(JJ,II)=(TOT)
                      TMAT2D(JJ+1,II+1)=(TOT)
#endif
                  ENDIF
                  WRITE(10,*) I,J,TOT
                  IF(I.EQ.J) THEN
                     DIAGTMAT(I)=ABS(TOT)
                     STATEORDER(I)=I
                  ENDIF

                  ! Subtract exchange integrals from Fock matrix.
!                 do K=1,nStates
!                     call CPMDGetOcc(K,F) ! F is the occupation number of state K.
!                      if (F.gt.1.0e-5_dp) then
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
!                 FockMat(j,i)=conjg(FockMat(i,j))

               ENDIF

            END DO
         END DO

C.. Now we order it if we need to
         IF(TORDER) THEN
            WRITE(6,*) "Re-ordering CPMD orbitals according ",
     &        "to one-electron energies."
           call sort (diagTMAT(nlock+1:nStatesUsed),
     &                stateOrder(nlock+1:nStatesUsed))
!           CALL NECI_SORT2(NSTATESUSED-NLOCK,DIAGTMAT(NLOCK+1),
!     &         STATEORDER(NLOCK+1))
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
            call LogMemDealloc(this_routine,tagHarInt)
            deallocate(PIInt)
            call LogMemDealloc(this_routine,tagPIInt)
         endif
         call halt_timer(proc_timer)
      END
