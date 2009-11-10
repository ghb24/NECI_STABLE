C.. Symmetries can be unset(=0), or have bits set for the different irreps included
C.. bit 0 corresponds to totally symmetric.
      
! JSS: if Abelian symmetry, then
! symmetry=0 corresponds to the symmetric representation.  It is not
! possible to unset symmetries for k-point jobs.
      
C.. To multiply symmetries, for each bit set in each of the two symmetries, we look up the 
C.. product in the symmetry table, and OR that with the product.

! JSS: In Abelian (k-point) symmetry, a representation can be described by "quantum" numbers.  
! The multiplication of two irreps is equal to the sum of such vectors
! describing the irreps, subject to a modulo of the periodic conditions
! due to the size of the symmetry cell (which corresponds to the k-point
! mesh).
! Note that many symmetry parameters for CPMD-NECI jobs are set in
! kpntrep.F in CPMD source.

         FUNCTION SYMPROD(ISYM1,ISYM2)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymTable,nProp,tAbelian,TwoCycleSymGens
         IMPLICIT NONE
         TYPE(Symmetry) ISYM1, ISYM2
         TYPE(Symmetry) SYMPROD
         TYPE(Symmetry) IS1,IS2
         INTEGER I,J,Abel1(3),Abel2(3)
         integer*8  ComposeAbelianSym
         if (TAbelian) then

             IF(TwoCycleSymGens) THEN
!For molecular systems, we can only have a maximum of 8 irreps, and so can do a simple
!xor to get the symmetry.

                 SymProd%s=IEOR(ISym1%s,ISym2%s)

             ELSE
             
                 call DecomposeAbelianSym(ISym1%s,Abel1)
                 call DecomposeAbelianSym(ISym2%s,Abel2)
!Slightly faster when calling a lot to do it in an array operation
                 Abel1(:)=modulo(Abel1(:)+Abel2(:),NProp(:))
!                do i=1,3
!                    Abel1(1)=modulo(Abel1(1)+Abel2(1),NProp(i))
!                    Abel1(2)=modulo(Abel1(2)+Abel2(2),NProp(i))
!                    Abel1(3)=modulo(Abel1(3)+Abel2(3),NProp(i))
!                end do
                 SymProd%s=ComposeAbelianSym(Abel1)
             ENDIF
         else
             IF(ISYM1%s.EQ.0.OR.ISYM2%s.EQ.0) THEN
                SYMPROD%s=0
                RETURN
             ENDIF
             IF (.not.allocated(SYMTABLE)) 
     &                               STOP 'SYMMETRY TABLE NOT ALLOCATED'
             IS1=ISYM1
             I=1
             SYMPROD%s=0
             DO WHILE(IS1%s.NE.0)
                IF(BTEST(IS1%s,0)) THEN
                   IS2=ISYM2
                   J=1
                   DO WHILE(IS2%s.NE.0)
                      IF(BTEST(IS2%s,0)) THEN
                         SYMPROD%s=IOR(SYMPROD%s,SYMTABLE(I,J)%s)
                      ENDIF
!  RSHIFT(,1)
                      IS2%s=ISHFT(IS2%s,-1)
                      J=J+1
                   ENDDO
                ENDIF
!RSHIFT(,1)
                IS1%s=ISHFT(IS1%s,-1)
                I=I+1
             ENDDO
         end if
         RETURN
      END

      FUNCTION SymConj(s2)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: tAbelian,nProp,SymConjTab,TwoCycleSymGens
         IMPLICIT NONE
         TYPE(Symmetry) s,SymConj,s2
         INTEGER i,AbelConj(3)
         integer*8 ComposeAbelianSym
         if (TAbelian) then
             IF(TwoCycleSymGens) THEN
!For molecular systems, we only have symmetry generators which are two cycles, and so the inverse of
!a symmetry is simply itself.

                 SymConj%s=s2%s
                 RETURN
             ELSE
! K-point symmetry has k_-i=k_i.  We store k-vectors from 0
! to N rather than -N/2 to N.  Hence k_-i=mod(-k_i+N,N) for
! each component of the vector.
! This also works for abelian groups made out of symmetry 
! generators which are 2-cycles
                 call DecomposeAbelianSym(s2%s,AbelConj)
                 do i=1,3
                   AbelConj(i)=modulo(-AbelConj(i)+NProp(i),NProp(i))
                 end do
                 SymConj%s=ComposeAbelianSym(AbelConj)

             ENDIF
         else
             SymConj%s=0
             s=s2
             i=1
             DO WHILE(s%s.NE.0)
!LSHIFT(,)
                ! 1_8 is 1 in integer*8: we need to have consistent
                ! kinds for bit-wise operations.
                IF(BTEST(s%s,0)) SymConj%s=IOR(SymConj%s,
     &         ISHFT(1_8,SymConjTab(I)-1))
!RSHIFT(,1)
                s%s=ISHFT(s%s,-1)
                i=i+1
             ENDDO
         end if
      END


      SUBROUTINE WRITESYMTABLE(IUNIT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymTable,nSym
         IMPLICIT NONE
         INTEGER IUNIT,I,J
         
         DO I=0,NSYM-1
            DO J=0,NSYM-1
               CALL WRITESYM(IUNIT,SYMTABLE(I+1,J+1),.FALSE.)
            ENDDO
            WRITE(IUNIT,*)
         ENDDO
         IF(NSYM.EQ.0) THEN
            WRITE(6,*) "No Symmetry table found."
         ENDIF 
      END

      LOGICAL FUNCTION LSYMSYM(SYM)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: tAbelian
         implicit none
         Type(Symmetry) SYM
         if (TAbelian) then
             LSymSym=Sym%s.eq.0
         else
             LSYMSYM=(SYM%s.EQ.0.OR.BTEST(SYM%s,0))
         end if
         RETURN
      END


C.. Generate a symmetry table with molpro symmetries.
C.. irreps are simply made out of up to three generators.
C.. MOLPRO classifies irreps with bits 0-3 corresponding to those generators.
C.. symmetry products are merely exclusive ors of the molpro irrep numbers
C.. We set each of the MOLPRO irreps to a bit in our symmetry specifier.
C.. A1 corresponds to bit 0 (i.e. irrep 1)
      SUBROUTINE GENMOLPSYMTABLE(NSYMMAX,G1,NBASIS,ARR,BRR)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nProp,PropBitLen,SymClasses,nSymLabels
         use SymData, only: tAbelian,SymLabels
         use SymData, only: tagSymLabels,tagSymClasses
         use global_utilities
         IMPLICIT NONE
         INTEGER NSYMMAX,nSymGen
         INTEGER I,J
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS)
         character(*), parameter :: this_routine='GenMolPSymTable'
        
         TAbelian=.true.
         PropBitLen=1
         nSymGen=INT(DLOG(NSYMMAX+0.D0)/DLOG(2.D0)+.4)
         WRITE(6,"(A,I3,A)") "Generating abelian symmetry table with",  &
     &      nSymGen, " generators" 
         WRITE(6,*) "Number of symmetry classes: ",nSymMax

         nprop=1
         nprop(1:nSymGen)=2
C.. Now generate a list of sym labels.
         NSYMLABELS=NSYMMAX
         allocate(SymLabels(nSymLabels))
         call LogMemAlloc('SymLabels',nSymLabels,SymmetrySize,
     &                                        this_routine,tagSymLabels)
         allocate(SymClasses(nBasis))
         call LogMemAlloc('SymClasses',nBasis,4,this_routine,
     &                                                    tagSymClasses)
         DO I=1,NBASIS,2
C.. place the sym label of each state in SymClasses(ISTATE).  For molp sym, this is 
C.. the log_2 of the symmetry bit string
            IF(G1(I)%Sym%s.EQ.0) THEN
C.. we don't have symmetry, so fake it.
               SymClasses((I+1)/2)=1
            ELSE
!          SymClasses((I+1)/2)=INT(DLOG(G1(I)%Sym%s+0.D0)/DLOG(2.D0)+1.4)
          SymClasses((I+1)/2)=G1(I)%Sym%s+1
            ENDIF
         ENDDO
C.. list the symmetry string of each sym label
         DO I=1,NSYMLABELS
!            SYMLABELS(I)%s=2**(I-1)
            SYMLABELS(I)%s=I-1
         ENDDO
      END

C.. Freeze the SYM LABELS and reps
C.. NHG is the old number of orbs
C.. NBASIS is the new number of orbs
C.. GG(I) is the new index of the (old) orb I

      SUBROUTINE FREEZESYMLABELS(NHG,NBASIS,GG,FRZ)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only:SymReps,SymClasses,SymClasses2,tagSymClasses2
         use global_utilities
         IMPLICIT NONE
         INTEGER NHG,NBASIS,GG(NHG)
         INTEGER I
         INTEGER NSL(NBASIS)
         LOGICAL FRZ
         character(*), parameter :: this_routine='FreezeSymLabels'
C.. SYMLABELS is used to classify all states which transform with the same symmetry
C.. for the excitation generation routines
C.. Each state's symmetry falls into a class SymClasses(ISTATE).
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SymLabels(ISYMLABEL)
C.. The characters of this class are stored in SYMLABELCHARS(1:NROT, SymClasses(ISTATE))
C.. The total number of symmetry labels is NSYMLABELS
!.. SYMREPS(1,IBASISFN) contains the numnber of the representation
!.. of which IBASISFN is a part.
         IF(.NOT.FRZ) THEN
            DO I=1,NHG,2
                IF(GG(I).NE.0) THEN
                    NSL((GG(I)+1)/2)=SymClasses((I+1)/2)
                ENDIF
            ENDDO
            DO I=1,NBASIS/2
                SymClasses(I)=NSL(I)
C               WRITE(6,*) "SL",I,SymClasses(I)
            ENDDO
            DO i=1,nhg
                IF(GG(i).ne.0) NSL(GG(i))=Symreps(1,i)
            enddo
            DO i=1,nbasis
                Symreps(1,i)=NSL(i)
            enddo
        ELSE
            IF(associated(SYMCLASSES2)) STOP 'Problem in freezing'
            allocate(SymClasses2(nBasis/2))
            call LogMemAlloc('SymClasses2',nBasis/2,4,this_routine,
     &                                                   tagSymClasses2)
            DO I=1,NHG,2
                IF(GG(I).NE.0) THEN
                    NSL((GG(I)+1)/2)=SymClasses((I+1)/2)
                ENDIF
            ENDDO
            DO I=1,NBASIS/2
                SymClasses2(I)=NSL(I)
C               WRITE(6,*) "SL",I,SymClasses(I)
            ENDDO
            DO i=1,nhg
                IF(GG(i).ne.0) NSL(GG(i))=Symreps(1,i)
            enddo
            DO i=1,nbasis
                Symreps(1,i)=NSL(i)
            enddo
        ENDIF
!        WRITE(6,*) "Sym Reps after Freezing"
!         DO i=1,nbasis
!             WRITE(6,*) i,Symreps(1,i),Symreps(2,i)
!         enddo
         
      END

      SUBROUTINE GENMOLPSYMREPS(NSYMMAX)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,Arr
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,Brr
         use SystemData, only: tSymIgnoreEnergies,nBasis,G1
         use SymData, only: SymReps,tagSymReps
         use global_utilities
         IMPLICIT NONE
         INTEGER NSYMMAX
         INTEGER I,J
!         TYPE(BasisFN) G1(NBASIS)
!         INTEGER NBASIS,BRR(NBASIS)
!         REAL*8 ARR(NBASIS)
         character(*), parameter :: this_routine='GENMOLPSYMREPS'
         LOGICAL tNew
         
C.. now work out which reps are degenerate and label them
         allocate(SymReps(2,nBasis))
         call LogMemAlloc('SymReps',2*nBasis,4,this_routine,tagSymReps)
         SymReps(:,:)=0
         J=0
         DO I=1,NBASIS
!             WRITE(6,*) I,nbasis
            tNew=.true.
            IF(tSymIgnoreEnergies.AND.MOD(I,2).EQ.0) THEN
!Pair even orbs up with the odd ones.
               SYMREPS(2,J)=SYMREPS(2,J)+1
               tNew=.false.
            ELSEIF(I.gt.1) THEN
                IF((ABS(ARR(I,1)-ARR(I-1,1)).LT.1.D-5)
     &       .AND.(G1(BRR(I))%Sym%s.EQ.G1(BRR(I-1))%Sym%s)) THEN
C.. we have the same degenerate rep as the previous entry
                    SYMREPS(2,J)=SYMREPS(2,J)+1
                  tNew=.false.
                ENDIF
            ENDIF
            IF(tNew) THEN
C.. we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ENDIF
            SYMREPS(1,BRR(I))=J
         ENDDO
!         WRITE(6,*) "Sym Reps MOLPRO"
!         DO i=1,nbasis
!             WRITE(6,*) i,Symreps(1,i),Symreps(2,i)
!         enddo
      END

C.. delete a symmetry table if one existed.
      SUBROUTINE ENDSYM()
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use global_utilities
         use SymData, only: SymTable,SymConjTab,SymReps,SymClasses
         use SymData, only: SymClasses2,SymLabels,SymLabelChars
         use SymData, only: IRREPChars,SymStatePairs,SymLabelList
         use SymData, only: SymLabelCounts
         use SymData, only: tagSymTable,tagSymConjTab,tagSymReps
         use SymData, only: tagSymClasses2,tagSymLabels,tagSymLabelChars
         use SymData, only: tagIRREPChars,tagSymStatePairs,tagSymClasses
         use SymData, only: tagSymLabelCounts,tagSymLabelList
         use SymData, only: SymPairProds,tagSymPairProds
         use SymExcitDataMod , only : SymLabelList2,SymLabelCounts2
         use SymExcitDataMod , only : OrbClassCount
         IMPLICIT NONE
         character(*), parameter :: this_routine='EndSym'
         if (allocated(SymTable)) then
             deallocate(SymTable)
             call LogMemDealloc(this_routine,tagSymTable)
         end if
         if (allocated(SymConjTab)) then
             deallocate(SymConjTab)
             call LogMemDealloc(this_routine,tagSymConjTab)
         end if
         if (allocated(SymReps)) then
             deallocate(SymReps)
             call LogMemDealloc(this_routine,tagSymReps)
         end if
         if (associated(SymClasses)) then
             deallocate(SymClasses)
             call LogMemDealloc(this_routine,tagSymClasses)
         end if
         nullify(SymClasses)
         if (associated(SymClasses2)) then
             deallocate(SymClasses2)
             call LogMemDealloc(this_routine,tagSymClasses2)
         end if
         nullify(SymClasses2)
         if (allocated(SymLabels)) then
             deallocate(SymLabels)
             call LogMemDealloc(this_routine,tagSymLabels)
         end if
         if (allocated(SymLabelChars)) then
             deallocate(SymLabelChars)
             call LogMemDealloc(this_routine,tagSymLabelChars)
         end if
         if (allocated(IRREPCHARS)) then
             deallocate(IRREPCHARS)
             call LogMemDealloc(this_routine,tagIRREPCHARS)
         end if
         if (allocated(SymStatePairs)) then
             deallocate(SymStatePairs)
             call LogMemDealloc(this_routine,tagSymStatePairs)
         end if
         if (allocated(SymLabelList)) then
             deallocate(SymLabelList)
             call LogMemDealloc(this_routine,tagSymLabelList)
         end if
         if (allocated(SymLabelCounts)) then
             deallocate(SymLabelCounts)
             call LogMemDealloc(this_routine,tagSymLabelCounts)
         end if
         if (allocated(SymLabelList2)) then
             deallocate(SymLabelList2)
         end if
         if (allocated(SymLabelCounts2)) then
             deallocate(SymLabelCounts2)
         end if
         if (allocated(OrbClassCount)) then
             deallocate(OrbClassCount)
         end if
         if (allocated(SymPairProds)) then
             deallocate(SymPairProds)
             call LogMemDealloc(this_routine,tagSymPairProds)
         end if
      END

C.. Precompute a list of the symmetry product of all pairs of symmetry labels
      SUBROUTINE GENSymStatePairs(NSTATES,FRZ)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymLabelList,SymClasses
         use SymData, only: SymClasses2,SymPairProds,SymPairProdSize
         use SymData, only: SymStatePairs,nSymPairProds,nSymLabels
         use SymData, only: tagSymPairProds,tagSymLabelList
         use SymData, only: tagSymLabelCounts,tagSymStatePairs
         use SymData, only: SymPairProd
         use SymData, only: tAbelianFastExcitGen,tAbelian
         use SymData, only: tStoreStateList
         use global_utilities
         IMPLICIT NONE
         INTEGER I,IDECOMP,TOT,NPRODS,NSTATES
         INTEGER TEMPLIST(NSTATES)
         LOGICAL FRZ
         character(*), parameter :: this_routine='GenSymStatePairs'

         if(tAbelianFastExcitGen.AND..NOT.tAbelian) THEN
            WRITE(6,*) "Fast Abelian excitation generators specified,"
     &      ,"but abelian symmetry not in use.  Using slow generators."
            tAbelianFastExcitGen=.false.
         endif

         if(.not.tAbelianFastExcitGen) THEN
!We are not in abelian fast excitgen - we are always storing state pairs, whether specified or not.
             tStoreStateList=.true.
         endif

         !May need to deallocate, since this info is allocated in storage of UMAT before freezing
         if (allocated(SymPairProds)) then
             deallocate(SymPairProds)
             call LogMemDealloc(this_routine,tagSymPairProds)
         end if
         if (allocated(SymLabelList)) then
             deallocate(SymLabelList)
             call LogMemDealloc(this_routine,tagSymLabelList)
         end if
         if (allocated(SymLabelCounts)) then
             deallocate(SymLabelCounts)
             call LogMemDealloc(this_routine,tagSymLabelCounts)
         end if
         if (allocated(SymStatePairs)) then
             deallocate(SymStatePairs)
             call LogMemDealloc(this_routine,tagSymStatePairs)
         end if

         allocate(SymLabelList(nStates))
         call LogMemAlloc('SymLabelList',nStates,4,this_routine,
     &                                                  tagSymLabelList)
         allocate(SymLabelCounts(2,nSymLabels))
         call LogMemAlloc('SymLabelCounts',2*nSymLabels,4,this_routine,
     &                                                tagSymLabelCounts)
C.. First deal with listing single states
         DO I=1,NSTATES
            SYMLABELLIST(I)=I
            IF(FRZ) THEN
               TEMPLIST(I)=SymClasses2(I)
            ELSE
               TEMPLIST(I)=SymClasses(I)
            ENDIF
         ENDDO
C.. order according to sym label, so SYMLABELLIST gets a list of states grouped under SYMLABEL
         CALL NECI_SORT2I(NSTATES,TEMPLIST,SYMLABELLIST)
         SYMLABELCOUNTS(:,:)=0
         SYMLABELCOUNTS(1,TEMPLIST(1))=1
         SYMLABELCOUNTS(2,TEMPLIST(1))=1
         DO I=2,NSTATES
            IF(TEMPLIST(I).NE.TEMPLIST(I-1)) THEN
C.. add a new sym label
               SYMLABELCOUNTS(2,TEMPLIST(I))=1
               SYMLABELCOUNTS(1,TEMPLIST(I))=I
C.. sort the symlabellist
               CALL NECI_SORTI(SYMLABELCOUNTS(2,TEMPLIST(I-1)),
     &            SYMLABELLIST(SYMLABELCOUNTS(1,TEMPLIST(I-1))))
            ELSE
               SYMLABELCOUNTS(2,TEMPLIST(I))=
     &            SYMLABELCOUNTS(2,TEMPLIST(I))+1
            ENDIF
         ENDDO
         CALL NECI_SORTI(SYMLABELCOUNTS(2,TEMPLIST(I-1)),
     &         SYMLABELLIST(SYMLABELCOUNTS(1,TEMPLIST(I-1))))
!         DO I=1,NSYMLABELS
!            WRITE(6,*) "NSL",I,SYMLABELCOUNTS(1,I),SYMLABELCOUNTS(2,I)
!     &         SymLabels(I)
!         ENDDO

C.. Now deal with pairs of states

         if(.not.tStoreStateList) then
!.. We don't bother listing all pairs of orbs, because we can calculate the number
!.. and they're easy to generate. 
!.. Instead of listing all pairs of states, we can list all pairs of sym classes (labels), 
!..  ordered according to their sym prod.
            allocate(SymPairProds(nSymLabels**2))
            call LogMemAlloc('SymPairProds',nSymLabels**2,
     &                   SymPairProdSize*8,this_routine,tagSymPairProds)
            SymPairProds(:)=SymPairProd(Symmetry(0),0,0,0,0)
C.. Now enumerate all pairs, and classify their product, but don't store them.
            nSymPairProds=0
            CALL GenSymPairs(nSymLabels,0)

C.. Now sort the SymPairProds into order
            CALL SORTSYMARR(nSymPairProds,SymPairProds,SymPairProdSize)
            TOT=0
            DO I=1,nSymPairProds
!               WRITE(6,"(I4,Z8,4I4)")
!     &             I,SymPairProds(I)%Sym,SymPairProds(I)%nPairs,
!     &           SymPairProds(I)%nIndex, SymPairProds(I)%nPairsStateSS,
!     &           SymPairProds(I)%nPairsStateOS
               SymPairProds(I)%nIndex=TOT
               TOT=TOT+SymPairProds(I)%nPairs
               SymPairProds(I)%nPairs=0
               SymPairProds(I)%nPairsStateSS=0
               SymPairProds(I)%nPairsStateOS=0
            ENDDO
            WRITE(6,*) TOT," Symmetry PAIRS"
            WRITE(6,*) NSYMPAIRPRODS,
     &         " DISTINCT ORBITAL PAIR PRODUCT SYMS"
            allocate(SymStatePairs(2,0:TOT-1))
            call LogMemAlloc('SymStatePairs',2*TOT,4,this_routine,
     &                                                 tagSymStatePairs)
            SymStatePairs(:,:)=0
            CALL GenSymPairs(nSymLabels,1)
!            WRITE(6,*) "Sym State Pairs"
!            DO I=0,TOT-1
!               WRITE(6,*) I,SymStatePairs(1:2,I)
!            ENDDO

!            WRITE(6,*) "SymLabelList: ",SymLabelList(:)
!            WRITE(6,*) "***","SymLabelCounts..."
!            WRITE(6,*) SymLabelCounts(1,:)
!            WRITE(6,*) "***"
!            WRITE(6,*) SymLabelCounts(2,:)

         else
!.. Non-abelian symmetry requires us to go through and work out all the possible pairs of orbs.
            allocate(SymPairProds(nSymLabels**2))
            call LogMemAlloc('SymPairProds',nSymLabels**2,
     &                   SymPairProdSize*8,this_routine,tagSymPairProds)
            SymPairProds=SymPairProd(Symmetry(0),0,0,0,0)
C.. Now enumerate all pairs, and classify their product, but don't store them.
            nSymPairProds=0
            CALL GENALLSymStatePairs(NSTATES,.FALSE.,FRZ)

C.. Now sort the SymPairProds into order
            CALL SORTSYMARR(nSymPairProds,SymPairProds,SymPairProdSize)
            TOT=0
!            write(6,*) "SymPairs",nSymPairProds
            DO I=1,nSymPairProds
               SymPairProds(I)%nIndex=TOT
               TOT=TOT+SymPairProds(I)%nPairs
               SymPairProds(I)%nPairs=0
            ENDDO
            WRITE(6,*) TOT," STATE PAIRS"
            WRITE(6,*) NSYMPAIRPRODS,
     &         " DISTINCT ORBITAL PAIR PRODUCT SYMS"
            allocate(SymStatePairs(2,0:TOT-1))
            call LogMemAlloc('SymStatePairs',2*TOT,4,this_routine,
     &                                                 tagSymStatePairs)
            SymStatePairs(:,:)=0
            CALL GENALLSymStatePairs(NSTATES,.TRUE.,FRZ)
         endif
      END

!= Generates Symmetry pairs in three passes:
!=  iPass
!=    0   Count number of pairs of symmetries for each possible symmetry product
!=    1   Store each pair of symmetries for each symmetry product
!=    2   Count the number of pairs of STATES for each pair of symmetries.  
!=
!=  NB This differs from GENAllSymStatePairs which goes through every possible 
!=     pair of states (and so is O(N^2)), and eventually stores them all.
!=     Here we only store pairs of symmetries (but calculate the number of pairs of states)
!=      This will only work for Abelian symmetries.

      SUBROUTINE GenSymPairs(nSymLabels,iPass)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabels,SymClasses,SymClasses2
         use SymData, only: SymPairProds,SymStatePairs,nSymPairProds
         use SymData, only: SymLabelCounts
         IMPLICIT NONE
         INTEGER iPass
         INTEGER I,J
         TYPE(Symmetry) PROD
         TYPE(Symmetry) SYMPROD
         INTEGER nSymLabels,iProd
         INTEGER iSS,iOS
         DO I=1,nSymLabels
            DO J=I,nSymLabels
!               WRITE(6,*) I,J
               PROD=SYMPROD(SymLabels(I),
     &               SymLabels(J))
               CALL FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
               IF(iProd.EQ.nSymPairProds+1) THEN
                  nSymPairProds=nSymPairProds+1
                  SymPairProds(iProd)%Sym=Prod
                  SymPairProds(iProd)%nIndex=0
                  SymPairProds(iProd)%nPairs=0
                  SymPairProds(iProd)%nPairsStateSS=0
                  SymPairProds(iProd)%nPairsStateOS=0
               ENDIF
                  
!.. iOS counts the number of pairs of spin-orbitals with the opposite spin, which give rise to the
!.. given symmetry product. iSS is for same spin orbital pairs.
               iOS=SymLabelCounts(2,I)*SymLabelCounts(2,J)
               if(i.ne.j) iOS=iOS*2
!.. Same spin has n(n-1)/2 if same state
               if(i.ne.J) then
                 iSS=SymLabelCounts(2,I)*SymLabelCounts(2,J)
               else
                 iSS=(SymLabelCounts(2,I)*(SymLabelCounts(2,J)-1))/2
               endif
               if(iOS.gt.0.or.iSS.gt.0) THEN
                  IF(iPass.eq.1) THEN
C.. put the pair into the list of pairs.
                     SymStatePairs(1,SymPairProds(iProd)%nIndex
     &                            +SymPairProds(iProd)%nPairs) =I
                     SymStatePairs(2,SymPairProds(iProd)%nIndex
     &                            +SymPairProds(iProd)%nPairs) =J
!                     WRITE(6,"(3I5,Z10,3I5)") 
!     &                  iProd,I,J,PROD,SymPairProds(iProd)%nIndex
!     &                            +SymPairProds(iProd)%nPairs,
!     &                           SymPairProds(iProd)%nIndex,
!     &                            SymPairProds(iProd)%nPairs

                  ENDIF
C.. increment the counter in the pairlist
                 SymPairProds(iProd)%nPairs=SymPairProds(iProd)%nPairs+1
                 SymPairProds(iProd)%nPairsStateOS=
     &              SymPairProds(iProd)%nPairsStateOS+iOS
                 SymPairProds(iProd)%nPairsStateSS=
     &              SymPairProds(iProd)%nPairsStateSS+iSS
!                write(6,*) "NN",SymLabelCounts(2,I),SymLabelCounts(2,J),
!     &           SymPairProds(iProd)%nPairsStateSS,
!     &           SymPairProds(iProd)%nPairsStateOS
               endif
            ENDDO
         ENDDO
      END

      SUBROUTINE GENALLSymStatePairs(NSTATES,TSTORE,FRZ)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabels,SymClasses,SymClasses2
         use SymData, only: SymPairProds,SymStatePairs,nSymPairProds
         IMPLICIT NONE
         LOGICAL TSTORE,FRZ
         INTEGER I,J
         TYPE(Symmetry) PROD
         TYPE(Symmetry) SYMPROD
         INTEGER NSTATES,iProd
         DO I=1,NSTATES
            DO J=I,NSTATES
C               WRITE(6,*) I,J,SymClasses(I),SymClasses(J)
                IF(FRZ) THEN
               PROD=SYMPROD(SymLabels(SymClasses2(I)),
     &               SymLabels(SymClasses2(J)))
                ELSE
               PROD=SYMPROD(SymLabels(SymClasses(I)),
     &               SymLabels(SymClasses(J)))
                ENDIF
               CALL FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
               IF(TSTORE) THEN
C.. put the pair into the list of pairs.
                  SymStatePairs(1,SymPairProds(iProd)%nIndex
     &                            +SymPairProds(iProd)%nPairs) =I
                  SymStatePairs(2,SymPairProds(iProd)%nIndex
     &                            +SymPairProds(iProd)%nPairs) =J
               ENDIF
               IF(iProd.EQ.nSymPairProds+1) THEN
                  nSymPairProds=nSymPairProds+1
                  SymPairProds(iProd)%Sym=Prod
                  SymPairProds(iProd)%nIndex=0
                  SymPairProds(iProd)%nPairs=0
               ENDIF
C.. incrememnt the counter in the pairlist
               SymPairProds(iProd)%nPairs=SymPairProds(iProd)%nPairs+1
            ENDDO
         ENDDO
      END

      SUBROUTINE FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymPairProd
         implicit none
         TYPE(SymPairProd) SymPairProds(nSymPairProds)
         INTEGER nSymPairProds,iProd
         TYPE(Symmetry) Prod
         LOGICAL SYMEQ
         DO iProd=1,nSymPairProds
            IF(SYMEQ(SymPairProds(iProd)%Sym,Prod)) EXIT
         ENDDO
      END
!.. SYMREPS is used to group together degenerate sets of orbitals of the same sym
!.. (e.g. the six orbitals which might make up a T2g set), and is used for working 
!.. out the symmetry of a determinant in GETSYM
!.. It uses that fact that even for non-abelian groups a completely filled degenerate symmetry set is totally symmetric.
!..  Thus each member of a set of states which when completely filled gives a totally symmetric det should be labelled with the same symrep
C..   SYMREPS(2,*) has two sets of data:
C.. SYMREPS(1,IBASISFN) contains the numnber of the representation
C.. of which IBASISFN is a part.
C.. SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
C.. The new method does the following:
C.. Identify all the completely filled reps.
C..   use ADDELECSYM to add together the momenta of these.
C..   These together are totally symmetric
C.. Identify all part-filled degenerate non-reduced representations
C..   Use SYMPROD and ADDELECSYM to generate the resultant symmetry of these


      SUBROUTINE GETSYM(NI2,NEL,G1,NBASISMAX,ISYM)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,tFixLz
         use SymData, only: SymReps
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),nBasisMax(5,*)
         TYPE(BasisFn) G1(*),ISym
         INTEGER I,J,NI2(NEL)
         INTEGER NREPS(NEL),NELECS(NEL),SSYM
         TYPE(Symmetry) SYMPROD
         LOGICAL ISCSF,ISC
         I=1
         NREPS(1:NEL)=0
         CALL SETUPSYM(NBASISMAX,ISYM)
         ISC=ISCSF(NI2,NEL)
         IF(tFixLz) CALL GetLz(NI2,NEL,ISYM%Ml)
         IF(ISC) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(NI2(I),NI(I),SSYM)
            ENDDO
         ELSE
!            CALL NECI_ICOPY(NEL,NI2,1,NI,1)
            NI(1:NEL)=NI2(1:NEL)
         ENDIF
         DO I=1,NEL
C.. Count all electrons in each rep
C.. NREPS(J) is the rep, and NELECS(J) is the number of electrons in that rep

            J=1
            DO WHILE(J.LT.NEL.AND.NREPS(J).NE.0)
               IF(NREPS(J).EQ.SYMREPS(1,NI(I))) THEN
C.. We've found the slot for the rep.  increment it and leave.
                  NELECS(J)=NELECS(J)+1
                  J=NEL
               ENDIF
               J=J+1
            ENDDO
            IF(J.LE.NEL) THEN
C.. need to put the new rep in a new space
               NREPS(J)=SYMREPS(1,NI(I))
               NELECS(J)=1
            ENDIF
         END DO
C.. now go through and see which are closed and which open
         DO I=1,NEL
            J=1
            DO WHILE(NREPS(J).NE.SYMREPS(1,NI(I)))
               J=J+1
            ENDDO
C.. electron NI(I) is in rep NREPS(J)
            IF(NELECS(J).NE.SYMREPS(2,NREPS(J))) THEN
C.. we don't have a closed shell
C.. add the sym product
               ISYM%Sym=SYMPROD(ISYM%Sym,G1(NI(I))%Sym)
            ENDIF
C.. add the momentum
            CALL ADDELECSYM(NI(I),G1,NBASISMAX,ISYM)
         ENDDO
C.. round the momentum
         CALL ROUNDSYM(ISYM,NBASISMAX)
         IF(ISC) CALL CSFGETSPIN(NI2,NEL,ISYM%Ms) 
         RETURN
      END

      SUBROUTINE GetLz(nI,NElec,Lz)
        use SystemData , only : G1
        INTEGER :: nI(NElec),NElec,Lz,i
        Lz=0
        do i=1,NElec
            Lz=Lz+G1(nI(i))%Ml
        enddo
      END SUBROUTINE GetLz

C.. Given a set of characters of states, generate all relevant irreps which span the set of characters.
      SUBROUTINE GENIRREPS(TKP,IMPROPER_OP,NROTOP)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: IRREPCHARS,nRot,SymLabelChars,nSym,tAbelian
         use SymData, only: SymLabels,nSymLabels
         IMPLICIT NONE
         INTEGER I,J,K
         LOGICAL LDO,LDO2
         TYPE(Symmetry) iDecomp
         INTEGER NEXTSYMLAB
         LOGICAL GETIRREPDECOMP
         COMPLEX*16 REPCHARS(NROT,NSYMLABELS*10)
         INTEGER NREPS,NROTOP
         REAL*8 NORM
         LOGICAL TKP,INV,IMPROPER_OP(NROTOP)
         NREPS=0
C.. Initialize the table with the totally symmetric rep.
         INV=.FALSE.
         DO I=1,NROT
            IRREPCHARS(I,1)=1
            IF(IMPROPER_OP(MOD(I-1,NROTOP)+1).and..not.TKP) INV=.TRUE.
         ENDDO
         NSYM=1
         IF(INV) THEN
            WRITE(6,*) "Inversion centre detected"
            NSYM=NSYM+1
C.. There's an inversion centre, so we can immediately create an A1u irrep
            DO I=1,NROT
               IF(IMPROPER_OP(MOD(I-1,NROTOP)+1)) THEN
                  IRREPCHARS(I,NSYM)=-1
               ELSE
                  IRREPCHARS(I,NSYM)=1
               ENDIF
            ENDDO
C            CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         ENDIF
         LDO=.TRUE.
         NEXTSYMLAB=1
         LDO2=.TRUE.
         DO WHILE(LDO.OR.LDO2)
C            CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
C            WRITE(6,*) NREPS," non-reducible"
C            CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
C.. First see if all the products of chars are decomposable
            LDO=.FALSE.
            NREPS=0
         lp1:DO I=1,NSYM
               DO J=I,NSYM
                  NREPS=NREPS+1
                  IF(NREPS.GT.NSYMLABELS*10) STOP 'TOO MANY REPS'
                  DO K=1,NROT
                     REPCHARS(K,NREPS)=DCONJG(IRREPCHARS(K,I))
     &                                    *IRREPCHARS(K,J)
                  ENDDO
                  
!                  WRITE(6,*) NREPS,"PROD",I,J
                  CALL N_MEMORY_CHECK
!                  CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDPRD")
                  IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM,TAbelian)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                     IF(ABS(NORM-NROT).LE.1.D-2) THEN
C.. if it's an irrep
                        NSYM=NSYM+1
                        IF(NSYM.GT.64) STOP "MORE than 64 irreps"
                        DO K=1,NROT
                           IRREPCHARS(K,NSYM)=REPCHARS(K,NREPS)
                        ENDDO
C                        CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
                        NREPS=NREPS-1
                        LDO=.TRUE.
                        EXIT lp1
                     ELSE
C                        WRITE(6,*) "IDECOMP:", IDECOMP,NORM,"SYMS:",NSYM
C                      CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"REMAIN")
C.. It's not an irrep, but we cannot reduce it.  Store only if we think we've got all the irreps.
!                        WRITE(6,*) "NR",NREPS,LDO2
                        IF(LDO2) NREPS=NREPS-1
!                        NREPS=NREPS-1
                     ENDIF
                  ELSE
C                     WRITE(6,*) "IDECOMP:", IDECOMP
                     NREPS=NREPS-1
                  ENDIF
               END DO
            END DO lp1
C            WRITE(6,*) LDO,NEXTSYMLAB,NSYMLABELS
            IF(LDO) CYCLE
C.. Check to see if the next symlabel's char is decomposable
        lp2: DO WHILE (NEXTSYMLAB.LE.NSYMLABELS)
               NREPS=NREPS+1
               IF(NREPS.GT.NSYMLABELS*10) STOP 'TOO MANY REPS'
               DO I=1,NROT
                  REPCHARS(I,NREPS)=SYMLABELCHARS(I,NEXTSYMLAB)
               ENDDO
C               CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDST ")
               IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM,TAbelian)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                  IF(ABS(NORM-NROT).LE.1.D-2) THEN
C.. if it's an irrep
                     NSYM=NSYM+1
                     IF(NSYM.GT.64) STOP "MORE than 64 irreps"
                     DO I=1,NROT
                        IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
                     ENDDO
C                     CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
                     NREPS=NREPS-1
                     LDO=.TRUE.
                     EXIT lp2
                  ELSE
C                     WRITE(6,*) "IDECOMP:", IDECOMP,NORM,"SYMS:",NSYM
C                     CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"REMAIN")
C.. It's not an irrep, but we cannot reduce it.  Store only if we think we've got all the irreps.
                     IF(LDO2) NREPS=NREPS-1
                  ENDIF
               ELSE
C                  WRITE(6,*) "IDECOMP:", IDECOMP
                  NREPS=NREPS-1
               ENDIF
               NEXTSYMLAB=NEXTSYMLAB+1
               IF(.NOT.LDO) THEN
C.. We've not manage to add any more irreps, so we have achieved self-consistency.  Do one more pass to check, saving all C.. non-reducible reps
                  LDO=.TRUE.
                  LDO2=.FALSE.
                  NREPS=0
               ENDIF
            END DO lp2
         ENDDO
C.. 
         WRITE(6,*) "IRREP TABLE"
         CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         IF(NREPS.GT.0) THEN
            WRITE(6,*) NREPS," non-reducible"
               CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
C            IF(NREPS.GT.1) THEN
               STOP "More than 1 non-reducible reps found."
C            ENDIF
C.. we can cope with a single reducible rep.
C            NSYM=NSYM+1
C            DO I=1,NROT
C               IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
C            ENDDO
         ENDIF
C..   Classify each of the symlabels with its decomposition into irreps
         DO I=1,NSYMLABELS
            CALL DECOMPOSEREP(SYMLABELCHARS(1,I),IDECOMP)
            SymLabels(I)=IDECOMP
         ENDDO
      END


C.. Display irrep table      
      SUBROUTINE WRITEIRREPTAB(IUNIT,CHARS,NROT,NSYM)
         IMPLICIT NONE
         INTEGER IUNIT,NROT,NSYM
         COMPLEX*16 CHARS(NROT,NSYM)
         CHARACTER*6 STR
         INTEGER I,J
         LOGICAL LCOMP,LREAL
         LCOMP=.FALSE.
         LREAL=.FALSE.
         DO I=1,NSYM
            DO J=1,NROT
               IF(ABS(REAL(CHARS(J,I))).GT.1.D-2.AND.
     &            ABS(DIMAG(CHARS(J,I))).GT.1.D-2) LCOMP=.TRUE.
              IF(ABS(REAL(CHARS(J,I))-NINT(REAL(CHARS(J,I))))
     &          .GT.1.D-2)
     &            LREAL=.TRUE.
            IF(ABS(DIMAG(CHARS(J,I))-NINT(DIMAG(CHARS(J,I))))
     &          .GT.1.D-2)
     &            LREAL=.TRUE.
            ENDDO
         ENDDO
         DO I=1,NSYM
            WRITE(STR,"(A,I3)"), "SYM", I
            CALL WRITECHARSF(IUNIT,CHARS(1,I),NROT,STR,LCOMP,LREAL)
         ENDDO
         WRITE(IUNIT,*)
      END 
C.. Display a line of characters
      SUBROUTINE WRITECHARS(IUNIT,CHARS,NROT,STR)
         IMPLICIT NONE
         INTEGER IUNIT,NROT
         COMPLEX*16 CHARS(NROT)
         INTEGER I,J
         CHARACTER*6 STR
         LOGICAL LCOMP,LREAL
C.. First do a check for the format
            LCOMP=.FALSE.
            LREAL=.FALSE.
            DO J=1,NROT
               IF(ABS(REAL(CHARS(J))).GT.1.D-2.AND.
     &            ABS(DIMAG(CHARS(J))).GT.1.D-2) LCOMP=.TRUE.
               IF(ABS(REAL(CHARS(J))-NINT(REAL(CHARS(J)))).GT.1.D-2)
     &            LREAL=.TRUE.
              IF(ABS(DIMAG(CHARS(J))-NINT(DIMAG(CHARS(J)))).GT.1.D-2)
     &            LREAL=.TRUE.
            ENDDO
            CALL WRITECHARSF(IUNIT,CHARS,NROT,STR,LCOMP,LREAL)
      END
      SUBROUTINE WRITECHARSF(IUNIT,CHARS,NROT,STR,LCOMP,LREAL)
         IMPLICIT NONE
         INTEGER IUNIT,NROT
         COMPLEX*16 CHARS(NROT)
         INTEGER I,J
         CHARACTER*6 STR
         LOGICAL LCOMP,LREAL
            WRITE(IUNIT,"(A6,A)",advance='no'), STR,":   "
            DO J=1,NROT
               IF(LCOMP) THEN
                  IF(LREAL) THEN
                     WRITE(IUNIT,"(A,2G16.9,A)",advance='no') "(",
     &                  NINT(REAL(CHARS(J))*1000)/1000.D0,
     &                  NINT(DIMAG(CHARS(J))*1000)/1000.D0
     &                  ,")"
                  ELSE
                     WRITE(IUNIT,"(A,2F6.3,A)",advance='no')
     &                                                 "(",CHARS(J),")"
                  ENDIF
               ELSE
                  IF(ABS(DIMAG(CHARS(J))).GT.1.D-2) THEN
C.. write in terms of I.
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G14.9,A)",advance='no')
     &                                                     CHARS(J)," "
                     ELSE                        
                        IF(ABS(DIMAG(CHARS(J))+1.D0).LT.1.D-2) THEN
                           WRITE(IUNIT,"(A)",advance='no') " -I "
                        ELSEIF(ABS(DIMAG(CHARS(J))-1.D0).LT.1.D-2) THEN
                           WRITE(IUNIT,"(A)",advance='no') "  I "
                        ELSE 
                         WRITE(IUNIT,"(I2,A)",advance='no')
     &                                           NINT(DIMAG(CHARS(J))),
     &                     "I "
                        ENDIF
                     ENDIF
                  ELSE
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G21.9)",advance='no')
     &                                           REAL(CHARS(J)), "    "
                     ELSE
                        WRITE(IUNIT,"(I3)",advance='no')
     &                                             NINT(REAL(CHARS(J)))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            WRITE(IUNIT,*)
      END

C.. Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
C.. CHARS at the end contains the remainder after the decomposition.
      SUBROUTINE DECOMPOSEREP(CHARSIN,IDECOMP)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nRot,nSym,tAbelian,IRREPCHARS
         IMPLICIT NONE
         TYPE(Symmetry) IDECOMP
         COMPLEX*16 CHARS(NROT),CHARSIN(NROT),TOT
         REAL*8 CNORM
         INTEGER I,J
         REAL*8 NORM,DIFF
         if (TAbelian) then
             ! We shouldn't be here!  Using symmetry "quantum" numbers
             ! rather than irreps.
             stop "Should not be decomposing irreps with Abelian sym"
         end if
         IDECOMP%s=0
         CALL DCOPY(NROT*2,CHARSIN,1,CHARS,1)
C         WRITE(6,*) "Decompose Rep"
C         CALL WRITECHARS(6,CHARS,NROT,"REP   ")
C,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+CHARS(J)*CHARS(J)
         ENDDO
         DO I=1,NSYM
            TOT=0
!            CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IR")
!            CALL WRITECHARS(6,CHARS(1),NROT,"CH")
            DO J=1,NROT
               TOT=TOT+DCONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
!            WRITE(6,*) I,TOT
            IF(TOT.NE.0) THEN
C.. Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+DCONJG(IRREPCHARS(J,I))*IRREPCHARS(J,I)
               ENDDO
C               WRITE(6,*) "IRREP ",I,(TOT+0.D0)/NORM
               DIFF=ABS(TOT-NINT(ABS(TOT/NORM))*NORM)
               IF(DIFF.GT.1.D-2) THEN
                  WRITE(6,*) 'Symmetry decomposition not complete'
                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IRREP ")
                  CALL WRITECHARS(6,CHARS,NROT,"CHARS ")
                  WRITE(6,*) "Dot product: ",(TOT+0.D0)/NORM,TOT,NORM
                  STOP 'Incomplete symmetry decomposition'
C.. The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
C.. Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
               ELSEIF(ABS(TOT).GT.1.D-2) THEN
C.. We've found an (ir)rep which is wholly in CHARS
                  IDECOMP%s=IBSET(IDECOMP%s,I-1)
                  CNORM=0
C                  WRITE(6,*) I,DIFF,TOT,TOT/NORM
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+DCONJG(CHARS(J))*CHARS(J)
                  ENDDO
C                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"DIRREP")
C                  CALL WRITECHARS(6,CHARS,NROT,"DCHARS")
               ENDIF
            ENDIF
         ENDDO
      END
   
 
C.. Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
C.. CHARS at the end contains the remainder after the decomposition.
C.. Return .FALSE. if the decomposition is complete and CHARS contains only 0.
C.. This is used internally in the symmetry routine and destroys CHARS.  For general decomposition,
C,, use DECOMPOSEREP
      LOGICAL FUNCTION GETIRREPDECOMP(CHARS,IRREPCHARS,NIRREPS,NROT,
     &         IDECOMP,CNORM,TAbelian)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         INTEGER NIRREPS, NROT
         TYPE(Symmetry) IDECOMP
         COMPLEX*16 IRREPCHARS(NROT,NIRREPS),CHARS(NROT)
         REAL*8 CNORM, NORM,DIFF
         COMPLEX*16 TOT
         INTEGER I,J
         logical TAbelian
         if (TAbelian) then
             ! We shouldn't be here!  Using symmetry "quantum" numbers
             ! rather than irreps.
             stop "Should not be decomposing irreps with Abelian sym"
         end if
         IDECOMP%s=0
C,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+DCONJG(CHARS(J))*CHARS(J)
         ENDDO
         DO I=1,NIRREPS
            TOT=0
            DO J=1,NROT
               TOT=TOT+DCONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
            IF(ABS(TOT).GE.1.D-2) THEN
C.. Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+DCONJG(IRREPCHARS(J,I))*IRREPCHARS(J,I)
               ENDDO
C               WRITE(6,*) "IRREP ",I,(TOT+0.D0)/NORM
C                CALL WRITECHARS(6,CHARS,NROT,"REP   ")
C                CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IRREP ")
               DIFF=ABS(TOT-NINT(ABS(TOT/NORM))*NORM)
               IF(DIFF.GE.1.D-2
     &            .AND.CNORM.EQ.NROT) THEN
C.. The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
C.. Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
C                  DO J=1,NROT
C                    IRREPCHARS(J,I)=IRREPCHARS(J,I)-CHARS(J)*TOT/CNORM
C                  ENDDO
C                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"NOW   ")
               ELSEIF(DIFF.LT.1.D-2) THEN
C.. We've found an (ir)rep which is wholly in CHARS
                  IDECOMP%s=IBSET(IDECOMP%s,I-1)
                  CNORM=0
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+DCONJG(CHARS(J))*CHARS(J)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         GETIRREPDECOMP=.FALSE.
         DO J=1,NROT
            IF(ABS(CHARS(J)).GT.1.D-2) GETIRREPDECOMP=.TRUE.
         ENDDO
      END


      SUBROUTINE GENSYMTABLE
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: IRREPCHARS,SymConjTab,tAbelian,nSym,SymTable
         use SymData, only: nRot,tagSymTable,tagSymConjTab
         use global_utilities
         IMPLICIT NONE
         INTEGER I,J,K
         LOGICAL GETIRREPDECOMP
         COMPLEX*16 CHARS(NROT)
         TYPE(Symmetry) IDECOMP
         REAL*8 CNORM
         character(*), parameter :: this_routine='GENSYMTABLE'
         allocate(SymTable(nSym,nSym))
         call LogMemAlloc('SymTable',nSym**2,SymmetrySize,this_routine,
     &                                                      tagSymTable)
         allocate(SymConjTab(nSym))
         call LogMemAlloc('SymConjTab',nSym,4,this_routine,
     &                                                    tagSymConjTab)
         DO I=1,NSYM
            DO K=1,NROT
               CHARS(K)=DCONJG(IRREPCHARS(K,I))
            ENDDO
            IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,
     &         IDECOMP,CNORM,TAbelian)) THEN
               WRITE(6,*) "Conjugate of SYM ",I,
     &            " not reducible,"
               CALL WRITECHARS(6,CHARS,NROT,"REMAIN")
               STOP "Symmetry table element not conjugable"
            ENDIF
            K=0
            DO WHILE(.NOT.BTEST(IDECOMP%s,0))
               K=K+1
!RSHIFT(,1)
               IDECOMP%s=ISHFT(IDECOMP%s,-1)
            ENDDO
            IF(IDECOMP%s.NE.1) THEN
               WRITE(6,*) "Conjugate of SYM ",I,
     &            " not a single SYM,"
               STOP
            ENDIF
            SymConjTab(I)=K+1
            DO J=I,NSYM
               DO K=1,NROT
                  CHARS(K)=IRREPCHARS(K,I)*IRREPCHARS(K,J)
               ENDDO
               IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,
     &            IDECOMP,CNORM,TAbelian)) THEN
                  WRITE(6,*) "Multiplication of SYMS ",I,J,
     &               " not reducible,"
                  CALL WRITECHARS(6,CHARS,NROT,"REMAIN")
                  STOP "Symmetry table element not reducible"
               ENDIF
               SYMTABLE(I,J)=IDECOMP
               SYMTABLE(J,I)=IDECOMP
C               WRITE(6,"(2I3,B12)") I,J,IDECOMP
            ENDDO
         ENDDO
         WRITE(6,*) "Symmetry, Symmetry Conjugate"
         DO I=1,NSYM
            WRITE(6,*) I,SymConjTab(I)
         ENDDO
      END

      SUBROUTINE GENSYMREPS(G1,NBASIS,ARR,BRR,DEGENTOL)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymReps,tAbelian,tagSymReps
         use global_utilities
         IMPLICIT NONE
         INTEGER I,J
         TYPE(BasisFN) G1(nBasis)
         INTEGER NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS,2)
         REAL*8 DEGENTOL
         LOGICAL SYMEQ
         character(*), parameter :: this_routine='GenSymReps'

C.. now work out which reps are degenerate and label them
         allocate(SymReps(2,nBasis))
         call LogMemAlloc('SymReps',2*nBasis,4,this_routine,tagSymReps)
         J=0
         DO I=1,NBASIS
C            WRITE(6,*) "SR2",I
            IF(I.GT.1.AND.ABS(ARR(I,2)-ARR(I-1,2)).LT.DEGENTOL
     &         .AND.(TAbelian.OR.SYMEQ(G1(I)%Sym,G1(I-1)%Sym))) THEN
C.. we have the same degenerate rep as the previous entry
               SYMREPS(2,J)=SYMREPS(2,J)+1
            ELSE
C.. we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ENDIF
            SYMREPS(1,I)=J
         ENDDO
C         DO I=1,NBASIS
C            WRITE(6,*) "SR1",SYMREPS(1,I),SYMREPS(2,I)
C         ENDDO   
      END

C.  Irrep symmetries are specified in SYM(5).
C.. if SYM(5)=0, we assume it's totally symmetric
C.. Other irreps contributing to the symmetry have bits set in 
C.. SYM.
C.. e.g. if irreps are a1,a2,b1,b2
      LOGICAL FUNCTION LCHKSYM(ISYM,JSYM)
         use SystemData, only: BasisFN,Symmetry
         IMPLICIT NONE
         TYPE(BASISFN) ISYM,JSYM
         INTEGER I,IS,JS
         TYPE(SYMMETRY) SYMPROD,SymConj
         LOGICAL LSYMSYM
         LCHKSYM=.TRUE.
         DO I=1,3
            IF(ISYM%K(I).NE.JSYM%K(I)) LCHKSYM=.FALSE.
         ENDDO
         IF(ISYM%Ms.NE.JSYM%Ms) LCHKSYM=.FALSE.
C.. if the symmetry product of I and J doesn't contain the totally
C.. symmetric irrep, we set sym to .FALSE.
        LCHKSYM=LCHKSYM.AND.LSYMSYM(SYMPROD(SymConj(ISYM%SYM),JSYM%SYM))
      RETURN
      END
      
      LOGICAL FUNCTION LCHKSYMD(NI,NJ,NEL,G1,NBASISMAX)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BASISFN) ISYM,JSYM,G1(*)
         INTEGER NEL,NI(NEL),NJ(NEL),nBasisMax(5,*)
         LOGICAL LCHKSYM
         CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
         CALL GETSYM(NJ,NEL,G1,NBASISMAX,JSYM)
         LCHKSYMD=LCHKSYM(ISYM,JSYM)
         RETURN
      END
C.. NBASISMAX descriptor (1,3)
C
C HUBBARD:
C  BITS
C. 0 Tilted
C. 1 non-pbc
C. 2 real-space
C.. which effects to values
C.. MOM SPACE
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
C (3,3)
C -2 Particle in a box
C -1 UEG
C 0 Hubbard
C 1 Generic spatial


C.. This only works for momentum variables - 1-4
      SUBROUTINE ADDELECSYM(IEL,G1,NBASISMAX,ISYM)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BASISFN) ISYM,G1(*)
         INTEGER IELEC,nBasisMax(5,*)
         INTEGER I,IEL,SSYM
         CALL GETUNCSFELEC(IEL,IELEC,SSYM)
        IF(NBASISMAX(1,3).LT.4) THEN
C.. Momentum space
            DO I=1,3
               ISYM%K(I)=ISYM%K(I)+G1(IELEC)%K(I)
            ENDDO
C.. Symmetry space
         ELSEIF(NBASISMAX(3,3).EQ.0.AND.NBASISMAX(1,3).GE.4) THEN
C.. We have no symmetries, so do nothing. (we're in real space)
C.. except Ms
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
C.. deal with momentum
            DO I=1,3
               ISYM%k(I)=ISYM%k(I)+G1(IELEC)%k(I)
            ENDDO
         ENDIF
         ISYM%MS=ISYM%MS+G1(IELEC)%MS
C.. SSYM keeps track of the total S change on adding this electron
C.. (it is +/-CSF_NSBASIS)
         I=ISYM%MS+0
         ISYM%Ms=I+SSYM
         RETURN
      END 
      
      SUBROUTINE ROUNDSYM(ISYM,NBASISMAX)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BasisFN) ISYM
         INTEGER nBasisMax(5,*)
         INTEGER I
         IF(NBASISMAX(3,3).EQ.-2) THEN
C.. particle in a box
C.. parity symmetries
            DO I=1,3
               ISYM%k(I)=MOD(ISYM%k(I),2)
            ENDDO
         ELSEIF(NBASISMAX(3,3).EQ.-1) THEN
C.. UEG (can't remember the symmetries of that
C.. probably momentum  conservation)
         ELSEIF(NBASISMAX(3,3).EQ.0) THEN
C.. Hubbard model
            IF(NBASISMAX(1,3).LT.2) THEN
C.. momentum conservation - various PBC issues
               CALL MOMPBCSYM(ISYM,NBASISMAX)
C            ELSEIF(NBASISMAX(1,3).EQ.2) THEN
C.. non-pbc mom space has parity symmetry
C               DO I=1,3
C                  ISYM(I)=MOD(ISYM(I),2)
C               ENDDO
            ELSEIF(NBASISMAX(1,3).GE.2) THEN
C.. we're in real space so no sym
               DO I=1,3
                  ISYM%k(I)=0
               ENDDO
            ENDIF
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
C.. Generic spatial symmetries
C..         We need do nothing.
C.. However, there is still momentum conservation - various PBC issues
               CALL MOMPBCSYM(ISYM,NBASISMAX)
         ENDIF
         RETURN 
      END

C.. NBASISMAX descriptor (1,3)
C
C HUBBARD:
C  BITS
C. 0 Tilted
C. 1 non-pbc
C. 2 real-space
C.. which effects to values
C.. MOM SPACE
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
C (3,3)
C -2 Particle in a box
C -1 UEG
C 0 Hubbard
C 1 Generic spatial


      SUBROUTINE MOMPBCSYM(K1,NBASISMAX)
C.. NB the third column of NBASISMAX tells us whether it is tilted
         IMPLICIT NONE
         INTEGER K1(3),nBasisMax(5,*)
         INTEGER J,LDIM,AX,AY,LENX,LENY,KK2,T1,T2
         REAL*8 R1,R2,NORM
         AX=NBASISMAX(1,4)
         AY=NBASISMAX(2,4)
         LENX=NBASISMAX(1,5)
         LENY=NBASISMAX(2,5)
         IF(NBASISMAX(1,3).EQ.0.OR.NBASISMAX(1,3).EQ.0) THEN
C.. A non-tilted lattice with PBC
            DO J=1,3
C..  non-tilted
               KK2=K1(J)
               LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
               KK2=MOD(KK2,LDIM)
               IF(KK2.LT.NBASISMAX(J,1)) KK2=KK2+LDIM
               IF(KK2.GT.NBASISMAX(J,2)) KK2=KK2-LDIM
               K1(J)=KK2 
            ENDDO
         ELSEIF(NBASISMAX(1,3).EQ.1) THEN
C.. we have a tilted lattice with PBC
C.. we want the a1,a2 components of k
            NORM=AX*AX+AY*AY
            R1=(AX*K1(1)+AY*K1(2))/NORM
            R2=(AX*K1(2)-AY*K1(1))/NORM
            R1=R1/LENX+0.5D0
            R2=R2/LENY+0.5D0
C.. T1= highest integer less than R1
            T1=INT(ABS(R1))
            IF(R1.LT.0.D0) THEN
               T1=-T1
               IF(T1.NE.R1) T1=T1-1
            ENDIF
            T2=INT(ABS(R2))
            IF(R2.LT.0.D0) THEN
               T2=-T2
               IF(T2.NE.R2) T2=T2-1
            ENDIF
            IF(R1.EQ.T1) T1=T1-1
            IF(R1.GT.1.D0.OR.R1.LE.0.D0) R1=R1-T1
            IF(R2.GE.1.D0.OR.R2.LT.0.D0) R2=R2-T2
            R1=(R1-0.5D0)*LENX
            R2=(R2-0.5D0)*LENY
            K1(1)=NINT(R1*AX-R2*AY)
            K1(2)=NINT(R1*AY+R2*AX)
         ENDIF
         RETURN
      END

      LOGICAL FUNCTION SYMLT(A,B)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) A,B
         IF(A%s.GE.0) THEN
            IF(B%s.GE.0) THEN
               SYMLT=A%s.LT.B%s
            ELSE
               SYMLT=.TRUE.
            ENDIF
         ELSE
            IF(B%s.GE.0) THEN
               SYMLT=.FALSE.
            ELSE
               SYMLT=A%s.LT.B%s
            ENDIF
         ENDIF    
         RETURN
      END
      LOGICAL FUNCTION SYMNE(A,B)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) A,B
         SYMNE=A%s.NE.B%s
         RETURN
      END
      LOGICAL FUNCTION SYMEQ(A,B)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) A,B
         SYMEQ=A%s.EQ.B%s
         RETURN
!Need to cope with 'unsigned integers'
      END
      LOGICAL FUNCTION SYMGT(A,B)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) A,B
         IF(A%s.GE.0) THEN
            IF(B%s.GE.0) THEN
               SYMGT=A%s.GT.B%s
            ELSE
               SYMGT=.FALSE.
            ENDIF
         ELSE
            IF(B%s.GE.0) THEN
               SYMGT=.TRUE.
            ELSE
               SYMGT=A%s.GT.B%s
            ENDIF
         ENDIF    
         RETURN
      END

      integer Function FindSymLabel(s)
         use SystemData, only: Symmetry
         use SymData, only: SymLabels,nSymLabels

         IMPLICIT NONE
         Type(Symmetry) s
         logical symeq
         integer i
         do i=1,nSymLabels
            if(symeq(SymLabels(i),s)) exit
         enddo
         if(i.gt.nSymLabels) i=0
         FindSymLabel=i
         return
      end
C.. A binary search to find VAL in TAB.  TAB is sorted
      SUBROUTINE BINARYSEARCHSYM(VAL,TAB,ROWLEN,LEN,LOC)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) VAL
         INTEGER LOC,LEN,ROWLEN
         INTEGER TAB(ROWLEN,LEN)
         INTEGER I,J,IFIRST,N,ILAST
         LOGICAL SYMNE,SYMLT,SYMGT,SYMEQ
         I=1
         J=LEN
         IFIRST=I
         ILAST=J
         DO WHILE(J-I.GE.1)
            N=(I+J)/2
C            WRITE(6,"(3I4)",advance='no') I,J,N
C            CALL WRITESYM(6,TAB(1,I),.FALSE.)
C            CALL WRITESYM(6,TAB(1,J),.FALSE.)
C            CALL WRITESYM(6,TAB(1,N),.TRUE.)
            IF(SYMLT(TAB(1,N),VAL).AND.I.NE.N) THEN
               IF(SYMNE(TAB(1,N),TAB(1,IFIRST))) IFIRST=N
C.. reset the lower limit
               I=N
            ELSEIF(SYMGT(TAB(1,N),VAL)) THEN
               IF(SYMNE(TAB(1,N),TAB(1,ILAST))) ILAST=N
C.. reset the upper limit
               J=N
            ELSEIF(SYMEQ(TAB(1,N),VAL)) THEN
C.. bingo, we've got it!
               LOC=N
               RETURN
            ELSE
C.. we've reached a situation where I and J's entries have the same value, and it's
C.. not the one we want.  Leave the loop.
               I=J
            ENDIF
         ENDDO
         IF(SYMEQ(TAB(1,I),VAL)) THEN
            LOC=I
         ELSEIF(SYMEQ(TAB(1,J),VAL)) THEN
            LOC=J
         ELSE
C.. Failure
            LOC=0
         ENDIF
      END
      SUBROUTINE GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &      TPARITY,IPARITY,TSETUP,TDONE,IMAX,ISYM)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: tAbelian
         IMPLICIT NONE
         INTEGER NEL,nBasisMax(5,*)
         INTEGER LMS
         TYPE(BasisFN) IPARITY,ISYM,IMax(2)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE,KALLOWED,TMORE2
         INTEGER ILEV
         integer*8 TotSymRep,MinSymRep,MaxSymRep
         IF(TSETUP) THEN
            DO ILEV=1,3
               IF(TPARITY) THEN
                  IMAX(1)%k(iLev)=IPARITY%k(ILEV)
                  IMAX(2)%k(ILEV)=IPARITY%k(ILEV)
               ELSE
                  IMAX(1)%k(iLev)=NBASISMAX(ILEV,1)
                  IMAX(2)%k(iLev)=NBASISMAX(ILEV,2)
C                  IF(NBASISMAX(1,3).EQ.2) THEN
C.. hubbard non-pbc mom space
C                     IMAX(ILEV,1)=IMAX(ILEV,1)*NEL
C                     IMAX(ILEV,2)=IMAX(ILEV,2)*NEL
C                  ENDIF
               ENDIF
            ENDDO
            IF(TSPN) THEN
               IMAX(1)%Ms=LMS
               IMAX(2)%Ms=LMS
            ELSE
               IMAX(1)%Ms=NBASISMAX(4,1)*NEL
               IMAX(2)%Ms=NBASISMAX(4,2)*NEL
            ENDIF
C.. If we're specifying a sym (TPARITY) in IPARITY(5), and
C.. we have a system with all 1D reducible orbs, then we put
C.. that into IMAX
            IF(NBASISMAX(5,2).NE.0.OR.TAbelian) THEN
               IF(TPARITY) THEN
                  IMAX(1)%Sym%s=IPARITY%Sym%s
                  IMAX(2)%Sym%s=IMAX(1)%Sym%s
               ELSE
                  IMAX(1)%Sym%s=MinSymRep(nBasisMax)
                  IMAX(2)%Sym%s=MaxSymRep(nBasisMax)
               ENDIF
            ELSE
C.. we've got a sym system with polydimensional irreps, which leads to
C.. dets with combinations of irreps, so we cannot put sym into blocks
                
!  JSS: if only 1D symmetries, then a determinant can only interact with
!  other determinants of the same symmetry.  This applies to Abelian
!  groups.  If there are multi-dimensional irreps, then this is no
!  longer the case, so we set the symmetries to be 0 (i.e. ignore
!  symmetry when generating determinants which interact).  This is not
!  equivalent to setting %s=0 if the Abelian case (which corresponds to
!  the totally symmetric irrep).
               IMAX(1)%Sym%s=0
               IMAX(2)%Sym%s=0
            ENDIF
            TDONE=.FALSE.
            CALL DOSYMLIMDEGEN(IMAX,NBASISMAX)
            ISym=IMax(1)
         ENDIF
         IF(TSETUP.AND.KALLOWED(ISYM,NBASISMAX)) RETURN
C.. Go to the next sym.
         TMORE2=.TRUE.
         TMORE=.TRUE.
         ILEV=5
         DO WHILE(TMORE2)
            DO WHILE (ILEV.GT.0)
               IF(ILEV.EQ.5) THEN
                  IF(IMAX(1)%Sym%s.NE.0) THEN
C.. symmetry specifiers are incremented by multiplying*2 (unless there are no syms counted)
                     ISYM%Sym%s=ISYM%Sym%s*2
                  ELSE
                     Call IncrSym(ISym%Sym)
                  ENDIF
                  IF(ISYM%Sym%s.EQ.IMAX(1)%Sym%s) THEN
                     ILEV=ILEV-1
                     IF(ILEV.EQ.0) THEN
                        TMORE2=.FALSE.
C.. If we've run out of syms, we give up
                        TMORE=.FALSE.
                     ENDIF
                  ELSEIF(KALLOWED(ISYM,NBASISMAX)) THEN
                     TMORE2=.FALSE.
                     ILEV=0
                  ENDIF
               ELSE
                  ISYM%k(ILEV)=ISYM%k(ILEV)+1
                  IF(ISYM%k(ILEV).GT.IMAX(2)%k(ILEV)) THEN
                     ISYM%k(ILEV)=IMAX(1)%k(ILEV)
                     ILEV=ILEV-1
                     IF(ILEV.EQ.0) THEN
                        TMORE2=.FALSE.
C.. If we've run out of syms, we give up
                        TMORE=.FALSE.
                     ENDIF
                  ELSEIF(ILEV.LT.4) THEN
C.. We've just incremented one of the higher columns, now go down to the
C.. lower ones.
                     ILEV=ILEV+1
                     ISYM%k(ILEV)=IMAX(1)%k(ILEV)-1
                     
                  ELSEIF(KALLOWED(ISYM,NBASISMAX)) THEN
                     TMORE2=.FALSE.
                     ILEV=0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         TDONE=.NOT.TMORE
      END
      SUBROUTINE DOSYMLIMDEGEN(IMAX,NBASISMAX)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BasisFN) IMax(2)
         INTEGER nBasisMax(5,*),I
         IF(NBASISMAX(3,3).EQ.0) THEN
            DO I=1,3
               IF(IMax(2)%k(I).NE.IMAX(1)%k(I)) IMAX(1)%k(I)=0
            ENDDO
         ENDIF
C.. always a spin symmetry
         IF(IMAX(1)%Ms.NE.IMAX(2)%Ms) IMAX(1)%Ms=0
      END
      SUBROUTINE GETSYMDEGEN(ISYM,NBASISMAX,IDEGEN)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BasisFN) ISym,ISym2
         INTEGER nBasisMax(5,*),IDEGEN,I,J
         LOGICAL KALLOWED,TDO
         IDEGEN=0
         IF(NBASISMAX(3,3).EQ.0) THEN
C.. Hubbard
            DO I=0,7
               TDO=.TRUE.
               DO J=1,3
                  IF(.NOT.BTEST(I,J-1)) THEN
                     ISYM2%k(J)=ISYM%k(J)
                  ELSE
                     ISYM2%k(J)=-ISYM%k(J)
                     IF(ISYM%k(J).EQ.0) TDO=.FALSE.
                  ENDIF
               ENDDO
               IF(TDO.AND.KALLOWED(ISYM2,NBASISMAX)) IDEGEN=IDEGEN+1
            ENDDO
         ELSE
            IDEGEN=1
         ENDIF
C.. Spin
         IF(ISYM%Ms.NE.0) IDEGEN=IDEGEN*2
      END

C.. Initialize symmetry to take into account the core electrons
      SUBROUTINE SETUPSYM(NBASISMAX,ISYM)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: FrozenSym
         IMPLICIT NONE
         INTEGER nBasisMax(5,*)
         TYPE(BasisFN) ISym
         ISym=FrozenSym
         RETURN
      END

      SUBROUTINE WRITEALLSYM(IUNIT,SYM,LTERM)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER IUNIT
         TYPE(BASISFN) SYM
         LOGICAL LTERM
         INTEGER J
         WRITE(IUNIT,"(4I5)",advance='no')
     &                                SYM%K(1),SYM%K(2),SYM%K(3),SYM%MS
         CALL WRITESYM(IUNIT,SYM%SYM,LTERM)
      END
      SUBROUTINE WRITESYM(IUNIT,SYM,LTERM)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,tCPMD
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nSym,tAbelian
         IMPLICIT NONE
         INTEGER IUNIT
         TYPE(SYMMETRY) SYM
         LOGICAL LTERM
         INTEGER Abel(3)
         INTEGER J
         IF(TAbelian) THEN
            CALL DecomposeAbelianSym(SYM,Abel)
            if (tCPMD) then
              WRITE(IUNIT,'(" (",I2,",",I2,",",I2,")",I2)',advance='no')
     &                                                  Abel(1:3)
            else
              WRITE(IUNIT,'(" (",I2,",",I2,",",I2,")",I2)',advance='no')
     &                                                  Abel(1:3),SYM%s
            end if
         ELSEIF(NSYM.LE.16) THEN
            WRITE(IUNIT,"(Z5)",advance='no') SYM
         ELSEIF(NSYM.LE.24) THEN
            WRITE(IUNIT,"(Z7)",advance='no') SYM
         ELSEIF(NSYM.LE.32) THEN
            WRITE(IUNIT,"(Z9)",advance='no') SYM
         ELSEIF(NSYM.LE.40) THEN
            WRITE(IUNIT,"(Z11)",advance='no') SYM
         ELSEIF(NSYM.LE.48) THEN
            WRITE(IUNIT,"(Z13)",advance='no') SYM
         ELSEIF(NSYM.LE.56) THEN
            WRITE(IUNIT,"(Z15)",advance='no') SYM
         ELSE
            WRITE(IUNIT,"(Z17)",advance='no') SYM
         ENDIF
         IF(LTERM) WRITE(IUNIT,*)
      END
      SUBROUTINE SetupFreezeAllSym(Sym)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: FrozenSym
         IMPLICIT NONE
         TYPE(BasisFN) Sym
C.. Set to be totally symmetric
         FrozenSym=Sym
      END
      SUBROUTINE SetupFreezeSym(Sym)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: FrozenSym
         IMPLICIT NONE
         TYPE(BasisFN) Sym
!         FrozenSym=NullBasisFn
C.. Set to be totally symmetric
         FrozenSym=Sym
      END
 
!Deal with K-point symmetries, using translational symmetry operations.
      ! JSS: use Abelian symmetry formulation (allows us to go beyond 64
      ! symmetry operations, and hence deal with larger k-point meshes).
      SUBROUTINE GenKPtIrreps(nTranslat,nKps,KpntInd,nStates)
      use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
      use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
      use SymData, only: SymLabels,SymLabelChars,nRot,nSymLabels,KPntSym
      use SymData, only: SymClasses,tagSymLabelChars,tagSymLabels
      use SymData, only: tagSymClasses
      use global_utilities
      IMPLICIT NONE
      INTEGER nTranslat,nKps,KpntInd(nStates)
      INTEGER I,J,nStates
      integer*8  ComposeAbelianSym
      type(Symmetry) SymConj,SymProd
      character(*), parameter :: this_routine='GenKPtIrreps'
      nSymLabels=nKps
      nRot=nTranslat
      allocate(SymLabelChars(nRot,nSymLabels))
      call LogMemAlloc('SymLabelChars',nSymLabels*nRot,16,this_routine,
     &                                                 tagSymLabelChars)
      allocate(SymLabels(nSymLabels))
      call LogMemAlloc('SymLabels',nSymLabels,4,this_routine,
     &                                                     tagSymLabels)
      allocate(SymClasses(nStates))
      call LogMemAlloc('SymClasses',nStates,4,this_routine,
     &                                                    tagSymClasses)
      SYMLABELCHARS=0.d0
      DO I=1,nStates
        SymClasses(I)=KpntInd(I)
        SymLabels(KPntInd(I))%s=ComposeAbelianSym(KpntSym(:,KPntInd(I)))
      END DO
      write (6,*) 
      write(6,'(a11," |",a13,"|",a10)')' K-vector',' Label ','Conjugate'
      write (6,'(39("-"))')
      do i=1,nSymLabels
        write (6,'("(",3i3,")"," | ")',advance='no')
!,z8," | ",z8)') 
     &                   KpntSym(:,I)
         call writesym(6,SymLabels(I),.false.)
         write(6,'(A)',advance='no') " | "
         call writesym(6,SymConj(SymLabels(I)),.true.)
      end do
!      write (6,'(/,a)') 'Symmetry Multiplication Table'
!      do i=1,nSymLabels
!        do j=1,nSymLabels
!          write (6,'(z12)',advance='no') SymProd(SymLabels(I),SymLabels(J))
!        end do
!        write (6,*) 
!      end do
!      write (6,'(/)') 



!        WRITE(6,*) "SYMMETRY CLASSES"
!        CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
!.. Allocate memory gor irreps.
!.. Assume there will be no more than 64 irreps
!        CALL N_MEMORY(IP_IRREPCHARS,NROT*64*2,"IRREPCH")
      END

      subroutine  DecomposeAbelianSym(ISym,AbelSym)
      ! Store the symmetry index as integer*8.  For Abelian symmetry
      ! we need to have 3 numbers stored in this index.  We store
      ! according to isym=\sum_i AbelSym(i)*32768**(i-1).
      ! This allows point groups with more than 64 irreps to be used in
      ! the point group is Abelian (as all translational/k-point
      ! symmetries are).
      ! Decompose the symmetry label back into the appropriate "quantum"
      ! numbers.
      ! Store the symmetry index as integer*8.  For Abelian symmetry
      ! we need to have 3 numbers stored in this index.  We store
      ! according to isym=1+\sum_i AbelSym(i)*32768**(i-1).
      ! Note that many symmetry parameters for CPMD-NECI jobs are set in
      ! kpntrep.F in CPMD source.
      ! Decompose the symmetry label back into the appropriate
      ! numbers...
      use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
      use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
      use SymData, only: PropBitLen
      implicit none
      integer*8 Isym
      integer  AbelSym(3)
!RShift
      AbelSym(3)=IShft(Isym,-(PropBitLen*2))
!RShift
      AbelSym(2)=Iand(IShft(ISym,-PropBitLen),2_8**PropBitLen-1)
      AbelSym(1)=Iand(Isym,2_8**PropBitLen-1)
      return
      end subroutine DecomposeAbelianSym

      integer*8 function ComposeAbelianSym(AbelSym)
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
          use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
          use SymData, only: PropBitLen
          implicit none
          integer  AbelSym(3)
          integer*8 TempVar
          TempVar=AbelSym(3)
!LShift
          ComposeAbelianSym=IShft(Tempvar,PropBitLen*2)
     &                     +IShft(AbelSym(2),PropBitLen)
     &                      +AbelSym(1)
      end function ComposeAbelianSym


      function TotSymRep()
          ! Our definition of the totally symmetric representation
          ! changes according to whether we're using Abelian/k-point
          ! symmetry or the standard symmetry.  It's just a matter of
          ! convenience, rather than some deep theoretical insight!
          use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
          use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,tUEG
          use SymData, only: tAbelian
          implicit none
          Type(Symmetry) TotSymRep 
          if (TAbelian.or.tUEG) then
              TotSymRep%s=0
          else
              TotSymRep%s=1
          end if
      end function TotSymRep

      ! nBasisMax might well be needed in the future in these functions.
      integer*8 function MinSymRep(nBasisMax)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: tAbelian
         implicit none
         integer nBasisMax(5,*)
         if(TAbelian) then
            MinSymRep=0
         else
            MinSymRep=0
         endif
      end function MinSymRep
      integer*8 function MaxSymRep(nBasisMax)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: tAbelian,nProp
         implicit none
         integer nBasisMax(5,*)
         integer abel(3)
         integer*8 ComposeAbelianSym
         abel(:)=nprop(:)-1
         if(TAbelian) then
            MaxSymRep=ComposeAbelianSym(abel)
         else
            MaxSymRep=0
         endif
      end function MaxSymRep
      subroutine IncrSym(Sym)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: tAbelian,nProp
         implicit none
         type(Symmetry) Sym
         integer abel(3),i
         logical lcont
         integer*8 ComposeAbelianSym
         if(TAbelian) then
            call DecomposeAbelianSym(Sym%s,abel)
            i=1
            lcont=.true.
            do while(i.lt.4.and.lcont)
               abel(i)=mod(abel(i)+1,nprop(i))
               lcont=abel(i).eq.0
               i=i+1
            enddo
            Sym%s=ComposeAbelianSym(abel)
         else
            Sym%s=0
         endif
      end subroutine IncrSym

      SUBROUTINE GETSYMTMATSIZE(Nirrep,nBasis,iSS,iSize)
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: SymLabelCounts,SymLabelCountsCum
        use SymData, only: SymLabelIntsCum
        use SymData, only: tagSymLabelIntsCum,tagSymLabelCountsCum
        use global_utilities
        implicit none
        integer Nirrep,nBasis,iSS,nBi,i,basirrep,t
        integer*8 iSize
        character, parameter :: this_routine='GetSymTMatSize'
        nBi=nBasis/iSS
        iSize=0
        allocate(SymLabelIntsCum(nIrrep))
        call LogMemAlloc('SymLabelIntsCum',nIrrep,4,this_routine,
     &                                               tagSymLabelIntsCum)
        allocate(SymLabelCountsCum(nIrrep))
        call LogMemAlloc('SymLabelCountsCum',nIrrep,4,this_routine,
     &                                             tagSymLabelCountsCum)
        SYMLABELINTSCUM(1:Nirrep)=0
        SYMLABELCOUNTSCUM(1:Nirrep)=0
        do i=1,Nirrep
            basirrep=SYMLABELCOUNTS(2,i)
            iSize=iSize+(basirrep*(basirrep+1))/2
            SYMLABELINTSCUM(i)=iSize
            IF(i.eq.1) THEN
                SYMLABELCOUNTSCUM(i)=0
            ELSE
                DO t=1,(i-1)
                    SYMLABELCOUNTSCUM(i)=SYMLABELCOUNTSCUM(i)+
     &                   SYMLABELCOUNTS(2,t)
                ENDDO
            ENDIF
            write(6,*) basirrep,SYMLABELINTSCUM(i),SYMLABELCOUNTSCUM(i)
            call flush(6)
        enddo
        iSize=iSize+2
        !This is to allow the index of '-1' in the array to give a zero value
      END SUBROUTINE GETSYMTMATSIZE


      SUBROUTINE WRITESYMCLASSES(NBASIS)
        USE HElem
        USE UMatCache
        use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
        use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
        use SymData, only: nSymLabels,SymClasses,SymLabelCounts
        IMPLICIT NONE
        INTEGER I,NBASIS
        
        DO I=1,NBASIS/2
            WRITE(13,*) I,SYMCLASSES(I)
            CALL FLUSH(13)
        ENDDO
        DO I=1,NSYMLABELS
            WRITE(13,*) I,SYMLABELCOUNTS(2,I)
        ENDDO
        WRITE(13,*) "******************"
        CALL FLUSH(13)
      END

!This routine *stochastically* finds the size of the determinant space. For certain symmetries, its hard to find the
!allowed size of the determinant space. However, it can be simply found using a MC technique.
      SUBROUTINE FindSymMCSizeofSpace(IUNIT)
         use SymData, only : TwoCycleSymGens
         use SystemData, only: nEl,G1,nBasis,nOccAlpha,nOccBeta,G1
         use SystemData, only: tUEG, NIfTot
         use Determinants, only : FDet
         use DetCalc, only : ICILevel
         use SystemData, only : CalcDetCycles, CalcDetPrint,tFixLz
         use mt95, only : genrand_real2 
         use soft_exit, only : ChangeVars
         use Parallel
         use DetBitops, only: EncodeBitDet
         IMPLICIT NONE
         INTEGER :: IUNIT,OverallSym,j,SpatOrbs,FDetMom
         INTEGER :: FDetiLut(0:NIfTot),ExcitLev,iLut(0:NIfTot)
         INTEGER :: FDetSym,TotalSym,TotalMom,alpha,beta,ierr,Momx,Momy
         INTEGER :: Momz
         INTEGER*8 :: Accept,AcceptAll,i
         INTEGER*8 :: ExcitBin(0:NEl),ExcitBinAll(0:NEl)
         REAL*8 :: FullSpace,r,Choose,Frac
         REAL*8 :: SizeLevel(0:NEl) 
         LOGICAL :: tTruncSpace,tDummy,tDummy2,tSoftExitFound
         LOGICAL :: tNotAllowed,tAcc

         IF((.not.TwoCycleSymGens).and.(.not.tUEG)) THEN
             WRITE(IUNIT,*) "Only for molecular abelian symmetry "      &
     &      //" calculations can the exact size of the determinant "    &
     &      //" space be calculated currently..."
             WRITE(IUNIT,*) "Skipping size of space calculation..."
             RETURN
         ENDIF

         WRITE(IUNIT,*) "Calculating exact size of symmetry-allowed "   &
     &       //"determinant space using MC..."
         WRITE(IUNIT,*) CalcDetCycles, " MC cycles will be used, and "  &
     &       //"statistics printed out every ",CalcDetPrint," cycles."
         FDetSym=0
         FDetMom=0
         ExcitBin(:)=0
         ExcitBinAll(:)=0

         do i=1,NEl
            FDetSym=IEOR(FDetSym,INT(G1(FDet(i))%Sym%S,4))
            IF(tFixLz) FDetMom=FDetMom+G1(FDet(i))%Ml
         enddo

         IF(ICILevel.gt.0) THEN
             tTruncSpace=.true.
         ELSE
             tTruncSpace=.false.
         ENDIF

         CALL EncodeBitDet(FDet,FDetiLut)

         WRITE(IUNIT,*) "Symmetry of HF determinant is: ",FDetSym
         IF(tFixLz) THEN
             WRITE(IUNIT,*) "Momentum of HF determinant is: ",FDetMom
         ENDIF

         SpatOrbs=nBasis/2

         Accept=0
         
         FullSpace=Choose(SpatOrbs,nOccAlpha)
         FullSpace=FullSpace*Choose(SpatOrbs,nOccBeta)

         WRITE(IUNIT,*) "Size of space neglecting all but Sz symmetry: "&
     &      ,FullSpace

         CALL FLUSH(IUNIT)

         IF(iProcIndex.eq.0) THEN
             OPEN(14,file="SpaceMCStats",status='unknown',              &
     &           form='formatted')
         ENDIF

         do i=1,CalcDetCycles

             TotalSym=0
             TotalMom=0
             ExcitLev=0
             Momx=0
             Momy=0
             Momz=0
             iLut(:)=0

             !Create random determinant (Correct Sz symmetry)
             !Loop over alpha electrons
             do j=1,nOccAlpha

                 tNotAllowed=.true.
                 do while(tNotAllowed)

                     CALL genrand_real2(r)
                     alpha=2*(INT(SpatOrbs*r)+1)
                     IF(.not.BTEST(iLut((alpha-1)/32)
     &                      ,mod((alpha-1),32))) THEN
                         !Has *not* been picked before
                         iLut((alpha-1)/32)=
     &                       IBSET(iLut((alpha-1)/32),mod(alpha-1,32))
                         tNotAllowed=.false.
                     ENDIF
                 enddo

                 TotalSym=IEOR(TotalSym,INT((G1(alpha)%Sym%S),4))
                 IF(tFixLz) THEN
                     TotalMom=TotalMom+G1(alpha)%Ml
                 ENDIF
                 IF(tUEG) THEN
                     Momx=Momx+G1(alpha)%k(1)
                     Momy=Momy+G1(alpha)%k(2)
                     Momz=Momz+G1(alpha)%k(3)
                 ENDIF
           IF(.not.BTEST(FDetiLut((alpha-1)/32),mod((alpha-1),32))) THEN
                     !orbital chosen is *not* in the reference determinant
                     ExcitLev=ExcitLev+1
                 ENDIF

                 !Test
!                 IF((alpha.lt.2).or.(alpha.gt.nBasis)) THEN
!                     CALL Stop_All("FindSymMCSizeofSpace","Error "      &
!     &                   //"calculating whether determinant is allowed")
!                 ENDIF

             enddo

             !Loop over beta electrons
             do j=1,nOccBeta

                 tNotAllowed=.true.
                 do while(tNotAllowed)
                     CALL genrand_real2(r)
                     beta=2*(INT(SpatOrbs*r)+1)-1
                     IF(.not.BTEST(iLut((beta-1)/32)
     &                      ,mod((beta-1),32))) THEN
                         !Has *not* been picked before
                         iLut((beta-1)/32)=
     &                       IBSET(iLut((beta-1)/32),mod(beta-1,32))
                         tNotAllowed=.false.
                     ENDIF
                 enddo

                 TotalSym=IEOR(TotalSym,INT((G1(beta)%Sym%S),4))
                 IF(tFixLz) THEN
                     TotalMom=TotalMom+G1(beta)%Ml
                 ENDIF
                 IF(tUEG) THEN
                     Momx=Momx+G1(beta)%k(1)
                     Momy=Momy+G1(beta)%k(2)
                     Momz=Momz+G1(beta)%k(3)
                 ENDIF
            IF(.not.BTEST(FDetiLut((beta-1)/32),mod((beta-1),32))) THEN
                     !orbital chosen is *not* in the reference determinant
                     ExcitLev=ExcitLev+1
                 ENDIF
                 
                 !Test
!                 IF((beta.lt.1).or.(beta.gt.(nBasis-1))) THEN
!                     CALL Stop_All("FindSymMCSizeofSpace","Error "      &
!     &                   //"calculating whether determinant is allowed")
!                 ENDIF

             enddo

             tAcc=.false.
             IF(TotalSym.eq.FDetSym) THEN
             !Allow/disallow the determinant
                 IF(tFixLz) THEN
                     IF(TotalMom.eq.FDetMom) THEN
                         IF(tTruncSpace) THEN
                             IF(ExcitLev.le.ICILevel) THEN
                                 Accept=Accept+1
                                 tAcc=.true.
                             ENDIF
                         ELSE
                             Accept=Accept+1
                             tAcc=.true.
                         ENDIF
                     ENDIF
                 ELSE
                     IF(tTruncSpace) THEN
                         IF(ExcitLev.le.ICILevel) THEN
                             IF(tUEG) THEN
                                 IF((Momx.eq.0).and.(Momy.eq.0)         &
     &                                 .and.(Momz.eq.0)) THEN
                                    Accept=Accept+1
                                    tAcc=.true.
                                 ENDIF
                             ELSE
                                 Accept=Accept+1
                                 tAcc=.true.
                             ENDIF
                         ENDIF
                     ELSE
                         IF(tUEG) THEN
                             IF((Momx.eq.0).and.(Momy.eq.0).and.        &
     &                               (Momz.eq.0)) THEN
                                Accept=Accept+1
                                tAcc=.true.
                             ENDIF
                         ELSE
                             Accept=Accept+1
                             tAcc=.true.
                         ENDIF
                     ENDIF
                 ENDIF
             ENDIF

             IF(tAcc) THEN
!Add to correct bin for the excitation level
                 ExcitBin(ExcitLev)=ExcitBin(ExcitLev)+1
             ENDIF

             
             IF(mod(i,CalcDetPrint).eq.0) THEN
                 !Write out statistics
#ifdef PARALLEL
!                 WRITE(6,*) Accept,AcceptAll
                 CALL MPI_Reduce(Accept,AcceptAll,1,                    &
     &            MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
                 CALL MPI_Reduce(ExcitBin(0:NEl),ExcitBinAll(0:NEl),    &
     &            NEl+1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
                 AcceptAll=Accept
                 ExcitBinAll(0:NEl)=ExcitBin(0:NEl)
#endif
                 Frac=REAL(AcceptAll,8)/REAL(i*nProcessors,8)
                 do j=0,NEl
                     SizeLevel(j)=(REAL(ExcitBinAll(j),8)/              &
     &                 REAL(AcceptAll,8))*Frac*FullSpace
                 enddo
                 IF(iProcIndex.eq.0) THEN
                     WRITE(14,"(2I16,G25.10,F30.5)",advance='no') i,    &
     &                  AcceptAll,Frac,Frac*FullSpace
                     do j=0,NEl
                         WRITE(14,"(F30.5)",advance='no') SizeLevel(j)
                     enddo
                     WRITE(14,"(A)") ""
                 ENDIF

                 AcceptAll=0
                 ExcitBinAll(0:NEl)=0

                 CALL ChangeVars(tDummy,tSoftExitFound,tDummy2)
                 IF(tSoftExitFound) EXIT

             ENDIF

         enddo

#ifdef PARALLEL
         CALL MPI_Reduce(Accept,AcceptAll,1,                            &
     &     MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_Reduce(ExcitBin(0:NEl),ExcitBinAll(0:NEl),            &
     &     NEl+1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
         AcceptAll=Accept
         ExcitBinAll(0:NEl)=ExcitBin(0:NEl)
#endif
         Frac=REAL(AcceptAll,8)/REAL(i*nProcessors,8)
         do j=0,NEl
             SizeLevel(j)=(REAL(ExcitBinAll(j),8)/                      &
     &         REAL(AcceptAll,8))*Frac*FullSpace
         enddo

         IF(iProcIndex.eq.0) THEN
             WRITE(14,"(2I16,G25.10,F30.5)",advance='no') i,AcceptAll,  &
     &           Frac,Frac*FullSpace
                 do j=0,NEl
                     WRITE(14,"(F30.5)",advance='no') SizeLevel(j)
                 enddo
                 WRITE(14,"(A)") ""
             CLOSE(14)
         ENDIF

         WRITE(IUNIT,*) "*** MC size of space: ",Frac*FullSpace
         WRITE(IUNIT,*) "Individual excitation level contributions: "
         do j=0,NEl
             WRITE(IUNIT,"(I5,F30.5)") j,SizeLevel(j)
         enddo
         CALL FLUSH(IUNIT)

      END SUBROUTINE FindSymMCSizeofSpace

!This routine finds the size of the determinant space in terms, including all symmetry allowed determinants.
!This is written to IUNIT. This is only available for molecular (i.e. abelian) systems with a maximum of eigth irreps.
!This is done in a very crude way. Feel free to optimise it!
      SUBROUTINE FindSymSizeofSpace(IUNIT)
         use SymData , only : TwoCycleSymGens
         use SystemData , only : nEl,G1,nBasis,nOccAlpha,nOccBeta
         use Determinants , only : FDet
         IMPLICIT NONE
         INTEGER :: ClassCounts(2,0:7),Lim0a,Lim0b,Lim1a,Lim1b,Lim2a
         INTEGER :: Lima(0:7),Limb(0:7),a0,a1,a2,a3,a4,a5,a6,a7,NAlph
         INTEGER :: b0,b1,b2,b3,b4,b5,b6,b7,NBet,i,IUNIT,OverallSym
         INTEGER :: FDetSym,Lim2b
         REAL*8 :: Choose,Space,SpaceGrow
         LOGICAL :: Sym(0:7)

         IF(.not.TwoCycleSymGens) THEN
             WRITE(IUNIT,*) "Only for molecular abelian symmetry "      &
     &      //" calculations can the exact size of the determinant "    &
     &      //" space be calculated currently..."
             WRITE(IUNIT,*) "Skipping size of space calculation..."
             RETURN
         ENDIF

         WRITE(IUNIT,*) "Calculating exact size of symmetry-allowed "   &
     &       //"determinant space..."
         FDetSym=0
         do i=1,NEl
            FDetSym=IEOR(FDetSym,INT(G1(FDet(i))%Sym%S,4))
         enddo
         WRITE(6,*) "Symmetry of HF determinant is: ",FDetSym
         CALL FLUSH(IUNIT)
         ClassCounts(:,:)=0
!First, we need to find the number of spatial orbitals in each symmetry irrep.
         do i=1,nBasis,1
             IF(G1(i)%Ms.eq.1) THEN
                 ClassCounts(1,INT(G1(i)%Sym%S,4))=                     &
     &              ClassCounts(1,INT(G1(i)%Sym%S,4))+1
             ELSE

                 ClassCounts(2,INT(G1(i)%Sym%S,4))=                     &
     &              ClassCounts(2,INT(G1(i)%Sym%S,4))+1
             ENDIF
         enddo
         do i=0,7
             IF(mod((ClassCounts(1,i)+ClassCounts(2,i)),2).ne.0) THEN
!                 STOP 'Error counting determinants'
                 WRITE(6,*) 'WARNING: Different number of symmetries    &
     &           between the alpha and beta orbitals.'
             ENDIF
!             ClassCounts(i)=ClassCounts(i)/2
         enddo

         Lima(0)=min(nOccAlpha,ClassCounts(1,0))
         Limb(0)=min(nOccBeta,ClassCounts(2,0))
         Lima(1)=min(nOccAlpha,ClassCounts(1,1))
         Limb(1)=min(nOccBeta,ClassCounts(2,1))
         Lima(2)=min(nOccAlpha,ClassCounts(1,2))
         Limb(2)=min(nOccBeta,ClassCounts(2,2))
         Lima(3)=min(nOccAlpha,ClassCounts(1,3))
         Limb(3)=min(nOccBeta,ClassCounts(2,3))
         Lima(4)=min(nOccAlpha,ClassCounts(1,4))
         Limb(4)=min(nOccBeta,ClassCounts(2,4))
         Lima(5)=min(nOccAlpha,ClassCounts(1,5))
         Limb(5)=min(nOccBeta,ClassCounts(2,5))
         Lima(6)=min(nOccAlpha,ClassCounts(1,6))
         Limb(6)=min(nOccBeta,ClassCounts(2,6))
         Lima(7)=min(nOccAlpha,ClassCounts(1,7))
         Limb(7)=min(nOccBeta,ClassCounts(2,7))
         Space=0.D0

!         WRITE(6,*) ClassCounts(:)
!         WRITE(6,*) "***"
!         WRITE(6,*) Lima(:),Limb(:)

!Loop over each irrep twice, once for alpha electrons and once for beta.
         do a0=0,Lima(0)
         do b0=0,Limb(0)
             IF(mod(a0+b0,2).eq.1) THEN
                 Sym(0)=.true.
             ELSE
                 Sym(0)=.false.
             ENDIF
         do a1=0,Lima(1)
         do b1=0,Limb(1)
             IF(mod(a1+b1,2).eq.1) THEN
                 Sym(1)=.true.
             ELSE
                 Sym(1)=.false.
             ENDIF
         do a2=0,Lima(2)
         do b2=0,Limb(2)
             IF(mod(a2+b2,2).eq.1) THEN
                 Sym(2)=.true.
             ELSE
                 Sym(2)=.false.
             ENDIF
         do a3=0,Lima(3)
         do b3=0,Limb(3)
             IF(mod(a3+b3,2).eq.1) THEN
                 Sym(3)=.true.
             ELSE
                 Sym(3)=.false.
             ENDIF
         do a4=0,Lima(4)
         do b4=0,Limb(4)
             IF(mod(a4+b4,2).eq.1) THEN
                 Sym(4)=.true.
             ELSE
                 Sym(4)=.false.
             ENDIF
         do a5=0,Lima(5)
         do b5=0,Limb(5)
             IF(mod(a5+b5,2).eq.1) THEN
                 Sym(5)=.true.
             ELSE
                 Sym(5)=.false.
             ENDIF
         do a6=0,Lima(6)
         do b6=0,Limb(6)
             IF(mod(a6+b6,2).eq.1) THEN
                 Sym(6)=.true.
             ELSE
                 Sym(6)=.false.
             ENDIF
         do a7=0,Lima(7)
         do b7=0,Limb(7)
             IF(mod(a7+b7,2).eq.1) THEN
                 Sym(7)=.true.
             ELSE
                 Sym(7)=.false.
             ENDIF

             OverallSym=0
             do i=0,7
                IF(Sym(i)) THEN
                    OverallSym=IEOR(OverallSym,i)
                ENDIF
            enddo
            IF(OverallSym.eq.FDetSym) THEN
                NAlph=a0+a1+a2+a3+a4+a5+a6+a7
                NBet=b0+b1+b2+b3+b4+b5+b6+b7

                IF((NAlph.eq.NOccAlpha).and.(NBet.eq.NOccBeta)) THEN

                    SpaceGrow=1.D0
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,0),a0)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,0),b0)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,1),a1)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,1),b1)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,2),a2)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,2),b2)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,3),a3)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,3),b3)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,4),a4)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,4),b4)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,5),a5)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,5),b5)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,6),a6)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,6),b6)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(1,7),a7)
                    SpaceGrow=SpaceGrow*Choose(ClassCounts(2,7),b7)
                    Space=Space+SpaceGrow
                ENDIF
            ENDIF

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         WRITE(IUNIT,"(A,G25.16)") " *EXACT* size of symmetry allowed " &
     &       //"space of determinants is: ",Space
         CALL FLUSH(IUNIT)

      END SUBROUTINE FindSymSizeofSpace


!This routine finds the size of the determinant space in terms, including all symmetry allowed determinants.
!This is written to IUNIT. This is only available for molecular (i.e. abelian) systems with a maximum of eigth irreps.
!This is done in a very crude way. Feel free to optimise it!
      SUBROUTINE FindSymSizeofTruncSpace(IUNIT)
         use SymData , only : TwoCycleSymGens
         use SystemData , only : nEl,G1,nBasis,nOccAlpha,nOccBeta,Brr
         use Determinants , only : FDet
         use DetCalc , only : ICILevel
         IMPLICIT NONE
         INTEGER :: ClassCountsOcc(0:7),Lim0a,Lim0b,Lim1a,Lim1b,Lim2a
         INTEGER :: ClassCountsVirt(0:7),Lim2b,NAlphOcc,NAlphVirt
         INTEGER :: ClassCountsOccMax(0:7),ClassCountsVirtMax(0:7)
         INTEGER :: LimaOcc(0:7),LimbOcc(0:7),LimaVirt(0:7)
         INTEGER :: LimbVirt(0:7)
         INTEGER :: a0o,a1o,a2o,a3o,a4o,a5o,a6o,a7o
         INTEGER :: a0v,a1v,a2v,a3v,a4v,a5v,a6v,a7v
         INTEGER :: b0o,b1o,b2o,b3o,b4o,b5o,b6o,b7o,OverallSym
         INTEGER :: b0v,b1v,b2v,b3v,b4v,b5v,b6v,b7v,NBetOcc,i,IUNIT
         INTEGER :: FDetSym,NBetVirt
         REAL*8 :: Choose,Space,SpaceGrow,SpaceGrow2
         LOGICAL :: Sym(0:7)

         IF(.not.TwoCycleSymGens) THEN
             WRITE(IUNIT,*) "Only for molecular abelian symmetry "      &
     &      //" calculations can the exact size of the determinant "    &
     &      //" space be calculated currently..."
             WRITE(IUNIT,*) "Skipping size of space calculation..."
             RETURN
         ENDIF

         WRITE(IUNIT,*) "Calculating exact size of symmetry-allowed "   &
     &       //"determinant space..."
         FDetSym=0
         do i=1,NEl
            FDetSym=IEOR(FDetSym,INT(G1(FDet(i))%Sym%S,4))
         enddo
         WRITE(6,*) "Symmetry of HF determinant is: ",FDetSym
         CALL FLUSH(IUNIT)
         ClassCountsOcc(:)=0
         ClassCountsVirt(:)=0
!First, we need to find the number of spatial orbitals in each symmetry irrep.
!We need to separate this into occupied and virtual. 
         do i=1,NEl,1
             ClassCountsOcc(INT(G1(BRR(i))%Sym%S,4))=                   &
     &           ClassCountsOcc(INT(G1(BRR(i))%Sym%S,4))+1
         enddo
 
         do i=NEL+1,nBasis,1
             ClassCountsVirt(INT(G1(BRR(i))%Sym%S,4))=                  &
     &           ClassCountsVirt(INT(G1(BRR(i))%Sym%S,4))+1
         enddo

!These are still in spin orbitals, so check there are multiple of 2 values in 
!each symmetry irrep and then divide by two because we deal with alpha and beta separately.         
         do i=0,7
         IF(mod((ClassCountsOcc(i)+ClassCountsVirt(i)),2).ne.0) THEN
             STOP 'Error counting determinants'
         ENDIF
         ClassCountsOccMax(i)=CEILING(REAL(ClassCountsOcc(i))/2.D0)
         ClassCountsVirtMax(i)=CEILING(REAL(ClassCountsVirt(i))/2.D0)
         ClassCountsOcc(i)=FLOOR(REAL(ClassCountsOcc(i))/2.D0)
         ClassCountsVirt(i)=FLOOR(REAL(ClassCountsVirt(i))/2.D0)
         
!         ClassCounts(i)=ClassCounts(i)/2
         enddo

         IF(nOccAlpha.gt.nOccBeta) THEN
             do i=0,7
                 LimaOcc(i)=min(nOccAlpha,ClassCountsOccMax(i))
                 LimbOcc(i)=min(nOccBeta,ClassCountsOcc(i))
                 LimaVirt(i)=min(ICILevel,ClassCountsVirtMax(i))
                 LimbVirt(i)=min(ICILevel,ClassCountsVirt(i))
             enddo
         ELSE
             do i=0,7
                 LimaOcc(i)=min(nOccAlpha,ClassCountsOcc(i))
                 LimbOcc(i)=min(nOccBeta,ClassCountsOccMax(i))
                 LimaVirt(i)=min(ICILevel,ClassCountsVirt(i))
                 LimbVirt(i)=min(ICILevel,ClassCountsVirtMax(i))
             enddo
         ENDIF
 
         Space=0.D0

!Loop over each irrep twice, once for alpha electrons and once for beta.
!a0 is the number of alpha electrons in symmetry 0.
!b0 is the number of beta electrons in symmetry 0.         
         do a0o=0,LimaOcc(0)
         do b0o=0,LimbOcc(0)
         do a0v=0,LimaVirt(0)
         do b0v=0,LimbVirt(0)
             IF(mod(a0o+b0o+a0v+b0v,2).eq.1) THEN
                 Sym(0)=.true.
             ELSE
                 Sym(0)=.false.
             ENDIF
         do a1o=0,LimaOcc(1)
         do b1o=0,LimbOcc(1)
         do a1v=0,LimaVirt(1)
         do b1v=0,LimbVirt(1)
             IF(mod(a1o+b1o+a1v+b1v,2).eq.1) THEN
                 Sym(1)=.true.
             ELSE
                 Sym(1)=.false.
             ENDIF
         do a2o=0,LimaOcc(2)
         do b2o=0,LimbOcc(2)
         do a2v=0,LimaVirt(2)
         do b2v=0,LimbVirt(2)
             IF(mod(a2o+b2o+a2v+b2v,2).eq.1) THEN
                 Sym(2)=.true.
             ELSE
                 Sym(2)=.false.
             ENDIF

         do a3o=0,LimaOcc(3)
         do b3o=0,LimbOcc(3)
         do a3v=0,LimaVirt(3)
         do b3v=0,LimbVirt(3)
             IF(mod(a3o+b3o+a3v+b3v,2).eq.1) THEN
                 Sym(3)=.true.
             ELSE
                 Sym(3)=.false.
             ENDIF
         do a4o=0,LimaOcc(4)
         do b4o=0,LimbOcc(4)
         do a4v=0,LimaVirt(4)
         do b4v=0,LimbVirt(4)
             IF(mod(a4o+b4o+a4v+b4v,2).eq.1) THEN
                 Sym(4)=.true.
             ELSE
                 Sym(4)=.false.
             ENDIF
         do a5o=0,LimaOcc(5)
         do b5o=0,LimbOcc(5)
         do a5v=0,LimaVirt(5)
         do b5v=0,LimbVirt(5)
             IF(mod(a5o+b5o+a5v+b5v,2).eq.1) THEN
                 Sym(5)=.true.
             ELSE
                 Sym(5)=.false.
             ENDIF
         do a6o=0,LimaOcc(6)
         do b6o=0,LimbOcc(6)
         do a6v=0,LimaVirt(6)
         do b6v=0,LimbVirt(6)
             IF(mod(a6o+b6o+a6v+b6v,2).eq.1) THEN
                 Sym(6)=.true.
             ELSE
                 Sym(6)=.false.
             ENDIF
         do a7o=0,LimaOcc(7)
         do b7o=0,LimbOcc(7)
         do a7v=0,LimaVirt(7)
         do b7v=0,LimbVirt(7)
             IF(mod(a7o+b7o+a7v+b7v,2).eq.1) THEN
                 Sym(7)=.true.
             ELSE
                 Sym(7)=.false.
             ENDIF

             OverallSym=0
             do i=0,7
                IF(Sym(i)) THEN
                    OverallSym=IEOR(OverallSym,i)
                ENDIF
             enddo
             IF(OverallSym.eq.FDetSym) THEN
                NAlphOcc=a0o+a1o+a2o+a3o+a4o+a5o+a6o+a7o
                NBetOcc=b0o+b1o+b2o+b3o+b4o+b5o+b6o+b7o
                NAlphVirt=a0v+a1v+a2v+a3v+a4v+a5v+a6v+a7v
                NBetVirt=b0v+b1v+b2v+b3v+b4v+b5v+b6v+b7v


                IF(((NAlphOcc+NAlphVirt).eq.NOccAlpha)                  &
     &          .and.((NBetOcc+NBetVirt).eq.NOccBeta)) THEN
                IF((NAlphVirt+NBetVirt).le.ICILevel) THEN

                IF(nOccAlpha.gt.nOccBeta) THEN

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(0),a0o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(0),b0o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(0),a0v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(0),b0v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(1),a1o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(1),b1o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(1),a1v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(1),b1v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(2),a2o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(2),b2o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(2),a2v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(2),b2v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(3),a3o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(3),b3o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(3),a3v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(3),b3v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(4),a4o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(4),b4o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(4),a4v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(4),b4v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(5),a5o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(5),b5o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(5),a5v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(5),b5v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(6),a6o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(6),b6o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(6),a6v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(6),b6v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(7),a7o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(7),b7o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(7),a7v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(7),b7v)

                Space=Space+SpaceGrow

                ELSE

                SpaceGrow=1.D0
                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(0),a0o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(0),b0o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(0),a0v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(0),b0v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(1),a1o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(1),b1o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(1),a1v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(1),b1v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(2),a2o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(2),b2o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(2),a2v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(2),b2v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(3),a3o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(3),b3o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(3),a3v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(3),b3v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(4),a4o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(4),b4o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(4),a4v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(4),b4v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(5),a5o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(5),b5o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(5),a5v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(5),b5v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(6),a6o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(6),b6o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(6),a6v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(6),b6v)

                SpaceGrow=SpaceGrow*Choose(ClassCountsOcc(7),a7o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsOccMax(7),b7o)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirt(7),a7v)
                SpaceGrow=SpaceGrow*Choose(ClassCountsVirtMax(7),b7v)

                Space=Space+SpaceGrow
                ENDIF
                ENDIF
                ENDIF
            ENDIF

         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo

         WRITE(IUNIT,"(A,G25.16)") " *EXACT* size of symmetry allowed " &
     &       //"space of determinants is: ",Space
         CALL FLUSH(IUNIT)

      END SUBROUTINE FindSymSizeofTruncSpace

