module sym_mod

use constants, only: dp,int64,sizeof_int
use SymExcitDataMod, only: SymTableLabels
use SystemData, only: tKpntSym, tNoSymGenRandExcits, thub
implicit none

contains

!   Symmetries can be unset(=0), or have bits set for the different irreps included
!   bit 0 corresponds to totally symmetric.
      
! JSS: if Abelian symmetry, then
! symmetry=0 corresponds to the symmetric representation.  It is not
! possible to unset symmetries for k-point jobs.
      
!   To multiply symmetries, for each bit set in each of the two symmetries, we look up the 
!   product in the symmetry table, and OR that with the product.

! JSS: In Abelian (k-point) symmetry, a representation can be described by "quantum" numbers.  
! The multiplication of two irreps is equal to the sum of such vectors
! describing the irreps, subject to a modulo of the periodic conditions
! due to the size of the symmetry cell (which corresponds to the k-point
! mesh).
! Note that many symmetry parameters for CPMD-NECI jobs are set in
! kpntrep.F in CPMD source.

         FUNCTION SYMPROD(ISYM1,ISYM2)
         use SystemData, only: Symmetry
         use SystemData, only: BasisFN
         use SymData, only: SymTable,nProp,tAbelian,TwoCycleSymGens, nSymLabels
         IMPLICIT NONE
         TYPE(Symmetry) ISYM1, ISYM2
         TYPE(Symmetry) SYMPROD
         TYPE(Symmetry) IS1,IS2
         INTEGER I,J,Abel1(3),Abel2(3)
         character(*), parameter :: this_routine = 'SYMPROD'
         if (TAbelian) then

             IF(TwoCycleSymGens) THEN
!For molecular systems, we can only have a maximum of 8 irreps, and so can do a simple
!xor to get the symmetry.

                 SymProd%s=IEOR(ISym1%s,ISym2%s)

             ELSE
             
                 call DecomposeAbelianSym(ISym1%s,Abel1)
                 call DecomposeAbelianSym(ISym2%s,Abel2)
!Slightly faster when calling a lot to do it in an array operation
                 Abel1=modulo(Abel1+Abel2,NProp)
                 SymProd%s=ComposeAbelianSym(Abel1)
             ENDIF
         else
             IF(ISYM1%s.EQ.0.OR.ISYM2%s.EQ.0) THEN
                SYMPROD%s=0
                RETURN
             ENDIF
             ! can i do a quick fix, if all the symmetries are always 1 ? 
             ! so the product should also be 1 always...
             if (isym1%s == isym2%s .and. isym1%s == 1) then 
                 symprod%s = 1
                 return
             end if
!              IF (.not.allocated(SYMTABLE)) call stop_all(this_routine, 'SYMMETRY TABLE NOT ALLOCATED')
             IS1=ISYM1
             I=1
             SYMPROD%s=0
             if (tHub .or. tKPntSym) then
                 ! try new hubbard kpoint symmetry storage..
                 if (isym1%s > nSymLabels .or. isym2%s > nSymLabels) then
                     return
                 endif
                 symprod = SYMTABLE(isym1%s,isym2%s)
             else
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
         end if
         RETURN
      END FUNCTION SYMPROD

      FUNCTION SymConj(s2)
         use SystemData, only: Symmetry
         use SystemData, only: BasisFN
         use SymData, only: tAbelian,nProp,SymConjTab,TwoCycleSymGens
         IMPLICIT NONE
         TYPE(Symmetry) s,SymConj,s2
         INTEGER i,AbelConj(3)
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
                ! 1_8 is 1 in integer(int64): we need to have consistent
                ! kinds for bit-wise operations.
                IF(BTEST(s%s,0)) SymConj%s=IOR(SymConj%s,ISHFT(1_int64,int(SymConjTab(I)-1,int64)))
!RSHIFT(,1)
                s%s=ISHFT(s%s,-1)
                i=i+1
             ENDDO
         end if
      END FUNCTION SYMCONJ


      SUBROUTINE WRITESYMTABLE(IUNIT)
         use SystemData, only: Symmetry
         use SystemData, only: BasisFN
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
      END SUBROUTINE WRITESYMTABLE

      LOGICAL FUNCTION LSYMSYM(SYM)
         use SystemData, only: Symmetry
         use SystemData, only: BasisFN
         use SymData, only: tAbelian
         implicit none
         Type(Symmetry) SYM
         if (TAbelian) then
             LSymSym=Sym%s.eq.0
         else
             LSYMSYM=(SYM%s.EQ.0.OR.BTEST(SYM%s,0))
         end if
         RETURN
      END FUNCTION LSYMSYM


!   Generate a symmetry table with molpro symmetries.
!   irreps are simply made out of up to three generators.
!   MOLPRO classifies irreps with bits 0-3 corresponding to those generators.
!   symmetry products are merely exclusive ors of the molpro irrep numbers
!   We set each of the MOLPRO irreps to a bit in our symmetry specifier.
!   A1 corresponds to bit 0 (i.e. irrep 1)
      SUBROUTINE GENMOLPSYMTABLE(NSYMMAX,G1,NBASIS)
         use SystemData, only: Symmetry,SymmetrySize
         use SystemData, only: BasisFN, tUEG, tHUB, treal
         use SymData, only: nProp,SymClasses,nSymLabels
         use SymData, only: tAbelian,SymLabels, TwoCycleSymGens
         use SymData, only: tagSymLabels,tagSymClasses
         use SymData, only: SymConjTab, tagSymConjTab
         use util_mod, only: int_fmt
         use global_utilities
         IMPLICIT NONE
         INTEGER NSYMMAX,nSymGen
         INTEGER I,ILABEL
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS
         character(*), parameter :: this_routine='GenMolPSymTable'

         TAbelian=.true.
         nSymGen=INT(log(NSYMMAX+0.0_dp)/log(2.0_dp)+.4_dp)
         WRITE(6,"(A,I3,A)") "  Generating abelian symmetry table with",&
            nSymGen, " generators" 
         WRITE(6,'(A,'//int_fmt(nSymMax)//')')                          &
                              "  Number of symmetry classes: ",nSymMax

         ! We actually use momentum conservation directly for the UEG and
         ! Hubbard mode so just fake the symmetry information here.
         ! WARNING: do *not* use SymConj etc for these systems without fixing
         ! this---functions which rely upon the wavevectors being encoded into
         ! a symmetry integer will not work.
         if (TwoCycleSymGens .or. tUEG .or. tHUB) then
             ! Set propogation information.
             ! If not TwoCycleSymGens we assume the user has already
             ! done so...
             nprop=1
             nprop(1:min(nSymGen,3))=2
         end if

!   Now generate a list of sym labels.
         NSYMLABELS=NSYMMAX
         allocate(SymLabels(nSymLabels))
         call LogMemAlloc('SymLabels',nSymLabels,SymmetrySize,this_routine,tagSymLabels)
         allocate(SymClasses(nBasis))
         call LogMemAlloc('SymClasses',nBasis,4,this_routine,tagSymClasses)
         allocate(SymConjTab(nSymlabels))
         call LogMemAlloc('SymConjTab',nSymlabels,4,this_routine,tagSymConjTab)
         if (TwoCycleSymGens .or. tUEG) then
             DO I=1,NBASIS,2
!   place the sym label of each state in SymClasses(ISTATE).  For molp sym, this is 
!   the log_2 of the symmetry bit string
                IF(G1(I)%Sym%s.EQ.0) THEN
!   we don't have symmetry, so fake it.
                   SymClasses((I+1)/2)=1
                ELSE
              SymClasses((I+1)/2)=int(G1(I)%Sym%s,sizeof_int)+1
                ENDIF
             ENDDO
!   list the symmetry string of each sym label
             DO I=1,NSYMLABELS
                SYMLABELS(I)%s=I-1
                ! Abelian representations are self-inverses if the group is
                ! real.
                SymConjTab(I) = I
             ENDDO
         else if (.not.tHUB .or. treal) then
             ! Hubbard symmetry info set up in GenHubMomIrrepsSymTable.
             ! except for the real-space lattice!
             symlabels(:)%s = -1
             do i = 1, nbasis, 2
                 do ilabel = 1, nsymlabels
                     if (symlabels(ilabel)%s == g1(i)%sym%s) then
                         ! Already found this symmetry label.
                         symclasses((i+1)/2) = ilabel
                         exit
                     else if (symlabels(ilabel)%s == -1) then
                         ! Have not found this label...
                         symclasses((i+1)/2) = ilabel
                         symlabels(ilabel)%s = g1(i)%sym%s
                         exit
                     end if
                 end do
             end do
             ! Find inverses.
             do ilabel = 1, nsymlabels
                 do i = 1, nsymlabels
                     if (SymEq(symlabels(i), SymConj(symlabels(ilabel)))) then
                         SymConjTab(ilabel) = i
                         exit
                     end if
                 end do
             end do
!             WRITE(6,*) "Label, Sym, SymConjLabel, SymConj, SymProd"
!             do i=1,nsymlabels
!                 WRITE(6,"(5I12)") i,symlabels(i),SymConjTab(i),symlabels(SymConjTab(i)),
!SYMPROD(symlabels(i),symlabels(SymConjTab(i)))
!             enddo
         end if
      END SUBROUTINE GENMOLPSYMTABLE

!   Freeze the SYM LABELS and reps
!   NHG is the old number of orbs
!   NBASIS is the new number of orbs
!   GG(I) is the new index of the (old) orb I

      SUBROUTINE FREEZESYMLABELS(NHG,NBASIS,GG,FRZ)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only:SymReps,SymClasses,SymClasses2,tagSymClasses2
         use global_utilities
         IMPLICIT NONE
         INTEGER NHG,NBASIS,GG(NHG)
         INTEGER I
         INTEGER NSL(NBASIS)
         LOGICAL FRZ
         character(*), parameter :: this_routine='FreezeSymLabels'
!   SYMLABELS is used to classify all states which transform with the same symmetry
!   for the excitation generation routines
!   Each state's symmetry falls into a class SymClasses(ISTATE).
!   The symmetry bit string, decomposing the sym label into its component irreps is in 
!   SymLabels(ISYMLABEL)
!   The characters of this class are stored in SYMLABELCHARS(1:NROT, SymClasses(ISTATE))
!   The total number of symmetry labels is NSYMLABELS
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
!               WRITE(6,*) "SL",I,SymClasses(I)
            ENDDO
            DO i=1,nhg
                IF(GG(i).ne.0) NSL(GG(i))=Symreps(1,i)
            enddo
            DO i=1,nbasis
                Symreps(1,i)=NSL(i)
            enddo
        ELSE
            IF(associated(SYMCLASSES2)) call stop_all(this_routine, 'Problem in freezing')
            allocate(SymClasses2(nBasis/2))
            call LogMemAlloc('SymClasses2',nBasis/2,4,this_routine,tagSymClasses2)
            DO I=1,NHG,2
                IF(GG(I).NE.0) THEN
                    NSL((GG(I)+1)/2)=SymClasses((I+1)/2)
                ENDIF
            ENDDO
            DO I=1,NBASIS/2
                SymClasses2(I)=NSL(I)
!               WRITE(6,*) "SL",I,SymClasses(I)
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
         
      END SUBROUTINE FREEZESYMLABELS

      SUBROUTINE GENMOLPSYMREPS()
         use SystemData, only: Symmetry,Arr
         use SystemData, only: BasisFN,Brr
         use SystemData, only: tSymIgnoreEnergies,nBasis,G1,tKPntSym
         use SymData, only: SymReps,tagSymReps
         use global_utilities
         IMPLICIT NONE
         INTEGER I,J
!         TYPE(BasisFN) G1(NBASIS)
!         INTEGER NBASIS,BRR(NBASIS)
!         real(dp) ARR(NBASIS)
         character(*), parameter :: this_routine='GENMOLPSYMREPS'
         LOGICAL tNew

         if(tKPntSym) then
             !These symmetry routines only work for cases where all irreps are their
             !own inverse. In systems with multiple kpoints, this will not be the
             !case. Setup the symreps for non-abelian symmetries.
             CALL GENSYMREPS(G1,NBASIS,ARR,1.e-6_dp)
             return
         endif
         
!   now work out which reps are degenerate and label them
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
                IF((ABS(ARR(I,1)-ARR(I-1,1)).LT.1.0e-5_dp) &
                   .AND.(G1(BRR(I))%Sym%s.EQ.G1(BRR(I-1))%Sym%s)) THEN
!   we have the same degenerate rep as the previous entry
                    SYMREPS(2,J)=SYMREPS(2,J)+1
                  tNew=.false.
                ENDIF
            ENDIF
            IF(tNew) THEN
!   we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ENDIF
            SYMREPS(1,BRR(I))=J
         ENDDO
!         WRITE(6,*) "Sym Reps MOLPRO"
!         DO i=1,nbasis
!             WRITE(6,*) i,Symreps(1,i),Symreps(2,i)
!         enddo
      END SUBROUTINE GENMOLPSYMREPS

!   delete a symmetry table if one existed.
      SUBROUTINE ENDSYM()
         use SystemData, only: Symmetry, BasisFN
         use global_utilities
         use SymData
         use SymExcitDataMod , only : SymLabelList2,SymLabelCounts2, sym_label_list_spat
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
         if (associated(SymIndex)) then
             deallocate(SymIndex)
             call LogMemDealloc(this_routine,tagSymIndex)
         end if
         if (associated(SymIndex2)) then
             deallocate(SymIndex2)
             call LogMemDealloc(this_routine,tagSymIndex2)
         end if
         if (allocated(SymLabelList2)) then
             deallocate(SymLabelList2)
         end if
         if (allocated(sym_label_list_spat)) then
             deallocate(sym_label_list_spat)
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
      END SUBROUTINE ENDSYM

!   Precompute a list of the symmetry product of all pairs of symmetry labels
      SUBROUTINE GENSymStatePairs(NSTATES,FRZ)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: SymLabelCounts,SymLabelList,SymClasses
         use SymData, only: SymClasses2,SymPairProds,SymPairProdSize
         use SymData, only: SymStatePairs,nSymPairProds,nSymLabels
         use SymData, only: tagSymPairProds,tagSymLabelList
         use SymData, only: tagSymLabelCounts,tagSymStatePairs
         use SymData, only: SymPairProd
         use SymData, only: tAbelianFastExcitGen,tAbelian
         use SymData, only: tStoreStateList
         use global_utilities
         use sort_mod
         IMPLICIT NONE
         INTEGER I,TOT,NSTATES
         INTEGER TEMPLIST(NSTATES), lo, hi
         LOGICAL FRZ
         character(*), parameter :: this_routine='GenSymStatePairs'

         if(tAbelianFastExcitGen.AND..NOT.tAbelian) THEN
            WRITE(6,*) "Fast Abelian excitation generators specified,", &
             "but abelian symmetry not in use.  Using slow generators."
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
         call LogMemAlloc('SymLabelList',nStates,4,this_routine,tagSymLabelList)
         allocate(SymLabelCounts(2,nSymLabels))
         call LogMemAlloc('SymLabelCounts',2*nSymLabels,4,this_routine,tagSymLabelCounts)
!   First deal with listing single states
         DO I=1,NSTATES
            SYMLABELLIST(I)=I
            IF(FRZ) THEN
               TEMPLIST(I)=SymClasses2(I)
            ELSE
               TEMPLIST(I)=SymClasses(I)
            ENDIF
         ENDDO
!   order according to sym label, so SYMLABELLIST gets a list of states grouped under SYMLABEL
         call sort (tempList, symLabelList) ! 1:nStates
         SYMLABELCOUNTS(:,:)=0
         SYMLABELCOUNTS(1,TEMPLIST(1))=1
         SYMLABELCOUNTS(2,TEMPLIST(1))=1
         DO I=2,NSTATES
            IF(TEMPLIST(I).NE.TEMPLIST(I-1)) THEN
!   add a new sym label
               SYMLABELCOUNTS(2,TEMPLIST(I))=1
               SYMLABELCOUNTS(1,TEMPLIST(I))=I
!   sort the symlabellist
               lo = symLabelcounts(1, tempList(i-1))
               hi = lo + symLabelCounts(2, tempList(i-1)) - 1
               call sort (symLabelList (lo:hi))
            ELSE
               SYMLABELCOUNTS(2,TEMPLIST(I))=SYMLABELCOUNTS(2,TEMPLIST(I))+1
            ENDIF
         ENDDO
         lo = symLabelcounts(1,tempList(i-1))
         hi = lo + symLabelCounts(2,tempList(i-1)) - 1
         call sort (symLabelList(lo:hi))
!         DO I=1,NSYMLABELS
!            WRITE(6,*) "NSL",I,SYMLABELCOUNTS(1,I),SYMLABELCOUNTS(2,I)
!     &         SymLabels(I)
!         ENDDO

!   Now deal with pairs of states

         if(.not.tStoreStateList) then
!.. We don't bother listing all pairs of orbs, because we can calculate the number
!.. and they're easy to generate. 
!.. Instead of listing all pairs of states, we can list all pairs of sym classes (labels), 
!..  ordered according to their sym prod.
            allocate(SymPairProds(nSymLabels**2))
            call LogMemAlloc('SymPairProds',nSymLabels**2, &
                         SymPairProdSize*8,this_routine,tagSymPairProds)
            SymPairProds(:)=SymPairProd(Symmetry(0),0,0,0,0)
!   Now enumerate all pairs, and classify their product, but don't store them.
            nSymPairProds=0
            CALL GenSymPairs(nSymLabels,0)

!   Now sort the SymPairProds into order
            call sort (symPairProds(1:nSymPairProds))
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
            WRITE(6,*) NSYMPAIRPRODS, " DISTINCT ORBITAL PAIR PRODUCT SYMS"
            allocate(SymStatePairs(2,0:TOT-1))
            call LogMemAlloc('SymStatePairs',2*TOT,4,this_routine,tagSymStatePairs)
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
            call LogMemAlloc('SymPairProds',nSymLabels**2,SymPairProdSize*8,this_routine,tagSymPairProds)
            SymPairProds=SymPairProd(Symmetry(0),0,0,0,0)
!   Now enumerate all pairs, and classify their product, but don't store them.
            nSymPairProds=0
            CALL GENALLSymStatePairs(NSTATES,.FALSE.,FRZ)

!   Now sort the SymPairProds into order
            call sort (symPairProds(1:nSymPairProds))
            TOT=0
!            write(6,*) "SymPairs",nSymPairProds
            DO I=1,nSymPairProds
               SymPairProds(I)%nIndex=TOT
               TOT=TOT+SymPairProds(I)%nPairs
               SymPairProds(I)%nPairs=0
            ENDDO
            WRITE(6,*) TOT," STATE PAIRS"
            WRITE(6,*) NSYMPAIRPRODS," DISTINCT ORBITAL PAIR PRODUCT SYMS"
            allocate(SymStatePairs(2,0:TOT-1))
            call LogMemAlloc('SymStatePairs',2*TOT,4,this_routine,tagSymStatePairs)
            SymStatePairs(:,:)=0
            CALL GENALLSymStatePairs(NSTATES,.TRUE.,FRZ)
         endif
      END SUBROUTINE GENSymStatePairs

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
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: SymLabels
         use SymData, only: SymPairProds,SymStatePairs,nSymPairProds
         use SymData, only: SymLabelCounts
         IMPLICIT NONE
         INTEGER iPass
         INTEGER I,J
         TYPE(Symmetry) PROD
         INTEGER nSymLabels,iProd
         INTEGER iSS,iOS
         DO I=1,nSymLabels
            DO J=I,nSymLabels
!               WRITE(6,*) I,J
               PROD=SYMPROD(SymLabels(I),SymLabels(J))
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
!   put the pair into the list of pairs.
                     SymStatePairs(1,SymPairProds(iProd)%nIndex+SymPairProds(iProd)%nPairs) = I
                     SymStatePairs(2,SymPairProds(iProd)%nIndex+SymPairProds(iProd)%nPairs) = J
!                     WRITE(6,"(3I5,Z10,3I5)") 
!     &                  iProd,I,J,PROD,SymPairProds(iProd)%nIndex
!     &                            +SymPairProds(iProd)%nPairs,
!     &                           SymPairProds(iProd)%nIndex,
!     &                            SymPairProds(iProd)%nPairs

                  ENDIF
!   increment the counter in the pairlist
                 SymPairProds(iProd)%nPairs=SymPairProds(iProd)%nPairs+1
                 SymPairProds(iProd)%nPairsStateOS=SymPairProds(iProd)%nPairsStateOS+iOS
                 SymPairProds(iProd)%nPairsStateSS=SymPairProds(iProd)%nPairsStateSS+iSS
!                write(6,*) "NN",SymLabelCounts(2,I),SymLabelCounts(2,J),
!     &           SymPairProds(iProd)%nPairsStateSS,
!     &           SymPairProds(iProd)%nPairsStateOS
               endif
            ENDDO
         ENDDO
      END SUBROUTINE GenSymPairs

      SUBROUTINE GENALLSymStatePairs(NSTATES,TSTORE,FRZ)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: SymLabels,SymClasses,SymClasses2
         use SymData, only: SymPairProds,SymStatePairs,nSymPairProds
         IMPLICIT NONE
         LOGICAL TSTORE,FRZ
         INTEGER I,J
         TYPE(Symmetry) PROD
         INTEGER NSTATES,iProd
         DO I=1,NSTATES
            DO J=I,NSTATES
!               WRITE(6,*) I,J,SymClasses(I),SymClasses(J)
                IF(FRZ) THEN
               PROD=SYMPROD(SymLabels(SymClasses2(I)),SymLabels(SymClasses2(J)))
                ELSE
               PROD=SYMPROD(SymLabels(SymClasses(I)),SymLabels(SymClasses(J)))
                ENDIF
               CALL FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
               IF(TSTORE) THEN
!   put the pair into the list of pairs.
                  SymStatePairs(1,SymPairProds(iProd)%nIndex+SymPairProds(iProd)%nPairs) =I
                  SymStatePairs(2,SymPairProds(iProd)%nIndex+SymPairProds(iProd)%nPairs) =J
               ENDIF
               IF(iProd.EQ.nSymPairProds+1) THEN
                  nSymPairProds=nSymPairProds+1
                  SymPairProds(iProd)%Sym=Prod
                  SymPairProds(iProd)%nIndex=0
                  SymPairProds(iProd)%nPairs=0
               ENDIF
!   incrememnt the counter in the pairlist
               SymPairProds(iProd)%nPairs=SymPairProds(iProd)%nPairs+1
            ENDDO
         ENDDO
      END SUBROUTINE GENALLSymStatePairs

      SUBROUTINE FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: SymPairProd
         implicit none
         INTEGER nSymPairProds,iProd
         TYPE(SymPairProd) SymPairProds(nSymPairProds)
         TYPE(Symmetry) Prod
         DO iProd=1,nSymPairProds
            IF(SYMEQ(SymPairProds(iProd)%Sym,Prod)) EXIT
         ENDDO
      END SUBROUTINE FindSymProd
!.. SYMREPS is used to group together degenerate sets of orbitals of the same sym
!.. (e.g. the six orbitals which might make up a T2g set), and is used for working 
!.. out the symmetry of a determinant in GETSYM
!.. It uses that fact that even for non-abelian groups a completely filled degenerate symmetry set is totally symmetric.
!..  Thus each member of a set of states which when completely filled gives a totally 
!symmetric det should be labelled with the same symrep
!     SYMREPS(2,*) has two sets of data:
!   SYMREPS(1,IBASISFN) contains the numnber of the representation
!   of which IBASISFN is a part.
!   SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
!   The new method does the following:
!   Identify all the completely filled reps.
!     use ADDELECSYM to add together the momenta of these.
!     These together are totally symmetric
!   Identify all part-filled degenerate non-reduced representations
!     Use SYMPROD and ADDELECSYM to generate the resultant symmetry of these

      subroutine getsym_wrapper(det, sym)

          use SystemData, only: G1, nel, nBasisMax, BasisFn

          integer, intent(in) :: det(nel)
          type(basisfn), intent(out) :: sym

          call getsym(det, nel, G1, nBasisMax, sym)

      end subroutine


      SUBROUTINE GETSYM(NI2,NEL,G1,NBASISMAX,ISYM)
         use SystemData, only: Symmetry, BasisFN, tFixLz, lnosymmetry
         use SymData, only: SymReps,tAbelian
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),nBasisMax(5,*)
         TYPE(BasisFn) G1(*),ISym
         INTEGER I,J,NI2(NEL)
         INTEGER NREPS(NEL),NELECS(NEL),SSYM
         LOGICAL iscsf_old,ISC
         if (lnosymmetry) then
             isym%sym%s = 0
             return
         end if
         I=1
         NREPS(1:NEL)=0
         CALL SETUPSYM(ISYM)
         ISC=iscsf_old(NI2,NEL)
         IF(tFixLz) THEN
            CALL GetLz(NI2,NEL,ISYM%Ml)
         ELSE
            ISYM%Ml=0
         ENDIF
         IF(ISC) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(NI2(I),NI(I),SSYM)
            ENDDO
         ELSE
!            CALL NECI_ICOPY(NEL,NI2,1,NI,1)
            NI(1:NEL)=NI2(1:NEL)
         ENDIF
         IF(tAbelian) THEN !For Abelian symmetry we don't need the symreps malarky.
            DO I=1,NEL 
               ISYM%Sym=SYMPROD(ISYM%Sym,G1(NI(I))%Sym)
!   add the momentum
               CALL ADDELECSYM(NI(I),G1,NBASISMAX,ISYM)
            ENDDO
         ELSE
               DO I=1,NEL
      !   Count all electrons in each rep
      !   NREPS(J) is the rep, and NELECS(J) is the number of electrons in that rep

                  J=1
                  DO WHILE(J.LT.NEL)
                     IF(NREPS(J).EQ.0) exit
                     IF(NREPS(J).EQ.SYMREPS(1,NI(I))) THEN
      !   We've found the slot for the rep.  increment it and leave.
                        NELECS(J)=NELECS(J)+1
                        J = NEL
                     ENDIF
                     J=J+1
                  ENDDO
                  IF(J.LE.NEL) THEN
      !   need to put the new rep in a new space
                     NREPS(J)=SYMREPS(1,NI(I))
                     NELECS(J)=1
                  ENDIF
               END DO
      !   now go through and see which are closed and which open
               DO I=1,NEL
                  J=1
                  DO WHILE(NREPS(J).NE.SYMREPS(1,NI(I)))
                     J=J+1
                  ENDDO
      !   electron NI(I) is in rep NREPS(J)
                  IF(NELECS(J).NE.SYMREPS(2,NREPS(J))) THEN
      !   we don't have a closed shell
      !   add the sym product
                     ISYM%Sym=SYMPROD(ISYM%Sym,G1(NI(I))%Sym)
                  ENDIF
      !   add the momentum
                  CALL ADDELECSYM(NI(I),G1,NBASISMAX,ISYM)
               ENDDO
         ENDIF
!   round the momentum
         CALL ROUNDSYM(ISYM,NBASISMAX)
         IF(ISC) CALL CSFGETSPIN(NI2,NEL,ISYM%Ms) 
         RETURN
      END SUBROUTINE GETSYM

      SUBROUTINE GetLz(nI,NElec,Lz)
        use SystemData , only : G1
        INTEGER :: NElec
        INTEGER :: nI(NElec),Lz,i
        Lz=0
        do i=1,NElec
            Lz=Lz+G1(nI(i))%Ml
        enddo
      END SUBROUTINE GetLz

!   Given a set of characters of states, generate all relevant irreps which span the set of characters.
      SUBROUTINE GENIRREPS(TKP,IMPROPER_OP,NROTOP)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: IRREPCHARS,nRot,SymLabelChars,nSym,tAbelian
         use SymData, only: SymLabels,nSymLabels
         IMPLICIT NONE
         INTEGER I,J,K
         LOGICAL LDO,LDO2
         TYPE(Symmetry) iDecomp
         INTEGER NEXTSYMLAB
         complex(dp) REPCHARS(NROT,NSYMLABELS*10)
         INTEGER NREPS,NROTOP
         real(dp) NORM
         LOGICAL TKP,INV,IMPROPER_OP(NROTOP)
         character(*), parameter :: this_routine = 'GENIRREPS'
         NREPS=0
!   Initialize the table with the totally symmetric rep.
         INV=.FALSE.
         DO I=1,NROT
            IRREPCHARS(I,1)=1
            IF(IMPROPER_OP(MOD(I-1,NROTOP)+1).and..not.TKP) INV=.TRUE.
         ENDDO
         NSYM=1
         IF(INV) THEN
            WRITE(6,*) "Inversion centre detected"
            NSYM=NSYM+1
!   There's an inversion centre, so we can immediately create an A1u irrep
            DO I=1,NROT
               IF(IMPROPER_OP(MOD(I-1,NROTOP)+1)) THEN
                  IRREPCHARS(I,NSYM)=-1
               ELSE
                  IRREPCHARS(I,NSYM)=1
               ENDIF
            ENDDO
!            CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         ENDIF
         LDO=.TRUE.
         NEXTSYMLAB=1
         LDO2=.TRUE.
         DO WHILE(LDO.OR.LDO2)
!            CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
!            WRITE(6,*) NREPS," non-reducible"
!            CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
!   First see if all the products of chars are decomposable
            LDO=.FALSE.
            NREPS=0
         lp1:DO I=1,NSYM
               DO J=I,NSYM
                  NREPS=NREPS+1
                  IF(NREPS.GT.NSYMLABELS*10) call stop_all(this_routine, 'TOO MANY REPS')
                  DO K=1,NROT
                     REPCHARS(K,NREPS)=CONJG(IRREPCHARS(K,I))*IRREPCHARS(K,J)
                  ENDDO
                  
!                  WRITE(6,*) NREPS,"PROD",I,J
!                  CALL N_MEMORY_CHECK
!                  CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDPRD")
                  IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,NSYM,NROT,IDECOMP,NORM,TAbelian)) THEN
!   CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                     IF(ABS(NORM-NROT).LE.1.0e-2_dp) THEN
!   if it's an irrep
                        NSYM=NSYM+1
                        IF(NSYM.GT.64) call stop_all(this_routine, "MORE than 64 irreps")
                        DO K=1,NROT
                           IRREPCHARS(K,NSYM)=REPCHARS(K,NREPS)
                        ENDDO
!                        CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
                        NREPS=NREPS-1
                        LDO=.TRUE.
                        EXIT lp1
                     ELSE
!                        WRITE(6,*) "IDECOMP:", IDECOMP,NORM,"SYMS:",NSYM
!                      CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"REMAIN")
!   It's not an irrep, but we cannot reduce it.  Store only if we think we've got all the irreps.
!                        WRITE(6,*) "NR",NREPS,LDO2
                        IF(LDO2) NREPS=NREPS-1
!                        NREPS=NREPS-1
                     ENDIF
                  ELSE
!                     WRITE(6,*) "IDECOMP:", IDECOMP
                     NREPS=NREPS-1
                  ENDIF
               END DO
            END DO lp1
!            WRITE(6,*) LDO,NEXTSYMLAB,NSYMLABELS
            IF(LDO) CYCLE
!   Check to see if the next symlabel's char is decomposable
        lp2: DO WHILE (NEXTSYMLAB.LE.NSYMLABELS)
               NREPS=NREPS+1
               IF(NREPS.GT.NSYMLABELS*10) call stop_all(this_routine, 'TOO MANY REPS')
               DO I=1,NROT
                  REPCHARS(I,NREPS)=SYMLABELCHARS(I,NEXTSYMLAB)
               ENDDO
!               CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDST ")
               IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,NSYM,NROT,IDECOMP,NORM,TAbelian)) THEN
!   CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                  IF(ABS(NORM-NROT).LE.1.0e-2_dp) THEN
!   if it's an irrep
                     NSYM=NSYM+1
                     IF(NSYM.GT.64) call stop_all(this_routine, "MORE than 64 irreps")
                     DO I=1,NROT
                        IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
                     ENDDO
!                     CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
                     NREPS=NREPS-1
                     LDO=.TRUE.
                     EXIT lp2
                  ELSE
!                     WRITE(6,*) "IDECOMP:", IDECOMP,NORM,"SYMS:",NSYM
!                     CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"REMAIN")
!   It's not an irrep, but we cannot reduce it.  Store only if we think we've got all the irreps.
                     IF(LDO2) NREPS=NREPS-1
                  ENDIF
               ELSE
!                  WRITE(6,*) "IDECOMP:", IDECOMP
                  NREPS=NREPS-1
               ENDIF
               NEXTSYMLAB=NEXTSYMLAB+1
               IF(.NOT.LDO) THEN
!   We've not manage to add any more irreps, so we have achieved self-consistency.  
!Do one more pass to check, saving all C.. non-reducible reps
                  LDO=.TRUE.
                  LDO2=.FALSE.
                  NREPS=0
               ENDIF
            END DO lp2
         ENDDO
!   
         WRITE(6,*) "IRREP TABLE"
         CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         IF(NREPS.GT.0) THEN
            WRITE(6,*) NREPS," non-reducible"
               CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
!            IF(NREPS.GT.1) THEN
               call stop_all(this_routine, "More than 1 non-reducible reps found.")
!            ENDIF
!   we can cope with a single reducible rep.
!            NSYM=NSYM+1
!            DO I=1,NROT
!               IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
!            ENDDO
         ENDIF
!     Classify each of the symlabels with its decomposition into irreps
         DO I=1,NSYMLABELS
            CALL DECOMPOSEREP(SYMLABELCHARS(1,I),IDECOMP)
            SymLabels(I)=IDECOMP
         ENDDO
      END SUBROUTINE GENIRREPS


!   Display irrep table      
      SUBROUTINE WRITEIRREPTAB(IUNIT,CHARS,NROT,NSYM)
         IMPLICIT NONE
         INTEGER IUNIT,NROT,NSYM
         complex(dp) CHARS(NROT,NSYM)
         CHARACTER(6) STR
         INTEGER I,J
         LOGICAL LCOMP,LREAL
         LCOMP=.FALSE.
         LREAL=.FALSE.
         DO I=1,NSYM
            DO J=1,NROT
               IF(ABS(REAL(CHARS(J,I))).GT.1.0e-2_dp.AND.ABS(AIMAG(CHARS(J,I))).GT.1.0e-2_dp) LCOMP=.TRUE.
               IF(ABS(REAL(CHARS(J,I))-NINT(REAL(CHARS(J,I)))).GT.1.0e-2_dp) LREAL=.TRUE.
               IF(ABS(AIMAG(CHARS(J,I))-NINT(AIMAG(CHARS(J,I)))).GT.1.0e-2_dp) LREAL=.TRUE.
            ENDDO
         ENDDO
         DO I=1,NSYM
            WRITE(STR,"(A,I3)") "SYM", I
            CALL WRITECHARSF(IUNIT,CHARS(1,I),NROT,STR,LCOMP,LREAL)
         ENDDO
         WRITE(IUNIT,*)
      END  SUBROUTINE WRITEIRREPTAB
!   Display a line of characters
      SUBROUTINE WRITECHARS(IUNIT,CHARS,NROT,STR)
         IMPLICIT NONE
         INTEGER IUNIT,NROT
         complex(dp) CHARS(NROT)
         INTEGER J
         CHARACTER(6) STR
         LOGICAL LCOMP,LREAL
!   First do a check for the format
            LCOMP=.FALSE.
            LREAL=.FALSE.
            DO J=1,NROT
               IF(ABS(REAL(CHARS(J))).GT.1.0e-2_dp.AND.ABS(AIMAG(CHARS(J))).GT.1.0e-2_dp) LCOMP=.TRUE.
               IF(ABS(REAL(CHARS(J))-NINT(REAL(CHARS(J)))).GT.1.0e-2_dp) LREAL=.TRUE.
              IF(ABS(AIMAG(CHARS(J))-NINT(AIMAG(CHARS(J)))).GT.1.0e-2_dp) LREAL=.TRUE.
            ENDDO
            CALL WRITECHARSF(IUNIT,CHARS,NROT,STR,LCOMP,LREAL)
      END SUBROUTINE WRITECHARS
      SUBROUTINE WRITECHARSF(IUNIT,CHARS,NROT,STR,LCOMP,LREAL)
         IMPLICIT NONE
         INTEGER IUNIT,NROT
         complex(dp) CHARS(NROT)
         INTEGER J
         CHARACTER(6) STR
         LOGICAL LCOMP,LREAL
            WRITE(IUNIT,"(A6,A)",advance='no') STR,":   "
            DO J=1,NROT
               IF(LCOMP) THEN
                  IF(LREAL) THEN
                     WRITE(IUNIT,"(A,2G16.9,A)",advance='no') "(",  &
                       NINT(REAL(CHARS(J))*1000)/1000.0_dp,           &
                       NINT(AIMAG(CHARS(J))*1000)/1000.0_dp           &
                       ,")"
                  ELSE
                     WRITE(IUNIT,"(A,2F6.3,A)",advance='no') "(",CHARS(J),")"
                  ENDIF
               ELSE
                  IF(ABS(AIMAG(CHARS(J))).GT.1.0e-2_dp) THEN
!   write in terms of I.
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G14.9,A)",advance='no') CHARS(J)," "
                     ELSE                        
                        IF(ABS(AIMAG(CHARS(J))+1.0_dp).LT.1.0e-2_dp) THEN
                           WRITE(IUNIT,"(A)",advance='no') " -I "
                        ELSEIF(ABS(AIMAG(CHARS(J))-1.0_dp).LT.1.0e-2_dp) THEN
                           WRITE(IUNIT,"(A)",advance='no') "  I "
                        ELSE 
                         WRITE(IUNIT,"(I2,A)",advance='no') NINT(AIMAG(CHARS(J))), "I "
                        ENDIF
                     ENDIF
                  ELSE
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G21.9)",advance='no') REAL(CHARS(J)), "    "
                     ELSE
                        WRITE(IUNIT,"(I3)",advance='no') NINT(REAL(CHARS(J)))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            WRITE(IUNIT,*)
      END SUBROUTINE WRITECHARSF

!   Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
!   CHARS at the end contains the remainder after the decomposition.
      SUBROUTINE DECOMPOSEREP(CHARSIN,IDECOMP)
         use SystemData, only: Symmetry, BasisFn
         use SymData, only: nRot,nSym,tAbelian,IRREPCHARS
         IMPLICIT NONE
         TYPE(Symmetry) IDECOMP
         complex(dp) CHARS(NROT),CHARSIN(NROT),TOT
         real(dp) CNORM
         INTEGER I,J
         real(dp) NORM,DIFF
         character(*), parameter :: this_routine = 'DECOMPOSEREP'
         if (TAbelian) then
             ! We shouldn't be here!  Using symmetry "quantum" numbers
             ! rather than irreps.
             call stop_all(this_routine, "Should not be decomposing irreps with Abelian sym")
         end if
         IDECOMP%s=0
         CALL DCOPY(NROT*2,CHARSIN,1,CHARS,1)
!         WRITE(6,*) "Decompose Rep"
!         CALL WRITECHARS(6,CHARS,NROT,"REP   ")
!,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+real(CHARS(J)*CHARS(J),dp)
         ENDDO
         DO I=1,NSYM
            TOT=0
!            CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IR")
!            CALL WRITECHARS(6,CHARS(1),NROT,"CH")
            DO J=1,NROT
               TOT=TOT+CONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
!            WRITE(6,*) I,TOT
            IF (abs(TOT) > 1.0e-12_dp) THEN
!   Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+real(CONJG(IRREPCHARS(J,I))*IRREPCHARS(J,I),dp)
               ENDDO
!               WRITE(6,*) "IRREP ",I,(TOT+0.0_dp)/NORM
               DIFF=ABS(TOT-NINT(ABS(TOT/NORM))*NORM)
               IF(DIFF.GT.1.0e-2_dp) THEN
                  WRITE(6,*) 'Symmetry decomposition not complete'
                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IRREP ")
                  CALL WRITECHARS(6,CHARS,NROT,"CHARS ")
                  WRITE(6,*) "Dot product: ",(TOT+0.0_dp)/NORM,TOT,NORM
                  call stop_all(this_routine, 'Incomplete symmetry decomposition')
!   The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
!   Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
               ELSEIF(ABS(TOT).GT.1.0e-2_dp) THEN
!   We've found an (ir)rep which is wholly in CHARS
                  IDECOMP%s=IBSET(IDECOMP%s,I-1)
                  CNORM=0
!                  WRITE(6,*) I,DIFF,TOT,TOT/NORM
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+real(CONJG(CHARS(J))*CHARS(J),dp)
                  ENDDO
!                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"DIRREP")
!                  CALL WRITECHARS(6,CHARS,NROT,"DCHARS")
               ENDIF
            ENDIF
         ENDDO
      END SUBROUTINE DECOMPOSEREP
   
 
!   Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
!   CHARS at the end contains the remainder after the decomposition.
!   Return .FALSE. if the decomposition is complete and CHARS contains only 0.
!   This is used internally in the symmetry routine and destroys CHARS.  For general decomposition,
!,, use DECOMPOSEREP
      LOGICAL FUNCTION GETIRREPDECOMP(CHARS,IRREPCHARS,NIRREPS,NROT,IDECOMP,CNORM,TAbelian)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         INTEGER NIRREPS, NROT
         TYPE(Symmetry) IDECOMP
         complex(dp) IRREPCHARS(NROT,NIRREPS),CHARS(NROT)
         real(dp) CNORM, NORM,DIFF
         complex(dp) TOT
         INTEGER I,J
         logical TAbelian
         character(*), parameter :: this_routine = 'GETIRREPDECOMP'
         if (TAbelian) then
             ! We shouldn't be here!  Using symmetry "quantum" numbers
             ! rather than irreps.
             call stop_all(this_routine, "Should not be decomposing irreps with Abelian sym")
         end if
         IDECOMP%s=0
!,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+real(CONJG(CHARS(J))*CHARS(J),dp)
         ENDDO
         DO I=1,NIRREPS
            TOT=0
            DO J=1,NROT
               TOT=TOT+CONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
            IF(ABS(TOT).GE.1.0e-2_dp) THEN
!   Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+real(CONJG(IRREPCHARS(J,I))*IRREPCHARS(J,I),dp)
               ENDDO
!               WRITE(6,*) "IRREP ",I,(TOT+0.0_dp)/NORM
!                CALL WRITECHARS(6,CHARS,NROT,"REP   ")
!                CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IRREP ")
               DIFF=ABS(TOT-NINT(ABS(TOT/NORM))*NORM)
               IF(DIFF.GE.1.0e-2_dp.AND. abs(CNORM - NROT) < 1.0e-12_dp) THEN
!   The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
!   Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
!                  DO J=1,NROT
!                    IRREPCHARS(J,I)=IRREPCHARS(J,I)-CHARS(J)*TOT/CNORM
!                  ENDDO
!                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"NOW   ")
               ELSEIF(DIFF.LT.1.0e-2_dp) THEN
!   We've found an (ir)rep which is wholly in CHARS
                  IDECOMP%s=IBSET(IDECOMP%s,I-1)
                  CNORM=0
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+real(CONJG(CHARS(J))*CHARS(J),dp)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         GETIRREPDECOMP=.FALSE.
         DO J=1,NROT
            IF(ABS(CHARS(J)).GT.1.0e-2_dp) GETIRREPDECOMP=.TRUE.
         ENDDO
      END FUNCTION GETIRREPDECOMP


      SUBROUTINE GENSYMTABLE
         use SystemData, only: Symmetry, BasisFN, SymmetrySize
         use SymData, only: IRREPCHARS,SymConjTab,tAbelian,nSym,SymTable
         use SymData, only: nRot,tagSymTable,tagSymConjTab
         use global_utilities
         IMPLICIT NONE
         INTEGER I,J,K
         complex(dp) CHARS(NROT)
         TYPE(Symmetry) IDECOMP
         real(dp) CNORM
         character(*), parameter :: this_routine='GENSYMTABLE'
         allocate(SymTable(nSym,nSym))
         call LogMemAlloc('SymTable',nSym**2,SymmetrySize,this_routine,tagSymTable)
         allocate(SymConjTab(nSym))
         call LogMemAlloc('SymConjTab',nSym,4,this_routine,tagSymConjTab)
         DO I=1,NSYM
            DO K=1,NROT
               CHARS(K)=CONJG(IRREPCHARS(K,I))
            ENDDO
            IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,IDECOMP,CNORM,TAbelian)) THEN
               WRITE(6,*) "Conjugate of SYM ",I," not reducible,"
               CALL WRITECHARS(6,CHARS,NROT,"REMAIN")
               call stop_all(this_routine, "Symmetry table element not conjugable")
            ENDIF
            K=0
            DO WHILE(.NOT.BTEST(IDECOMP%s,0))
               K=K+1
!RSHIFT(,1)
               IDECOMP%s=ISHFT(IDECOMP%s,-1)
            ENDDO
            IF(IDECOMP%s.NE.1) THEN
               WRITE(6,*) "Conjugate of SYM ",I," not a single SYM,"
               call stop_all(this_routine, 'Incorrect sym conjugate')
            ENDIF
            SymConjTab(I)=K+1
            DO J=I,NSYM
               DO K=1,NROT
                  CHARS(K)=IRREPCHARS(K,I)*IRREPCHARS(K,J)
               ENDDO
               IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,IDECOMP,CNORM,TAbelian)) THEN
                  WRITE(6,*) "Multiplication of SYMS ",I,J," not reducible,"
                  CALL WRITECHARS(6,CHARS,NROT,"REMAIN")
                  call stop_all(this_routine, "Symmetry table element not reducible")
               ENDIF
               SYMTABLE(I,J)=IDECOMP
               SYMTABLE(J,I)=IDECOMP
!               WRITE(6,"(2I3,B12)") I,J,IDECOMP
            ENDDO
         ENDDO
         WRITE(6,*) "Symmetry, Symmetry Conjugate"
         DO I=1,NSYM
            WRITE(6,*) I,SymConjTab(I)
         ENDDO
      END SUBROUTINE GENSYMTABLE

      SUBROUTINE GENSYMREPS(G1,NBASIS,ARR,DEGENTOL)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: SymReps,tAbelian,tagSymReps
         use global_utilities
         IMPLICIT NONE
         INTEGER I,J
         INTEGER NBASIS
         TYPE(BasisFN) G1(nBasis)
         real(dp) ARR(NBASIS,2)
         real(dp) DEGENTOL
         logical lTmp
         character(*), parameter :: this_routine='GenSymReps'

!   now work out which reps are degenerate and label them
         allocate(SymReps(2,nBasis))
         call LogMemAlloc('SymReps',2*nBasis,4,this_routine,tagSymReps)
         J=0
         DO I=1,NBASIS
!            WRITE(6,*) "SR2",I
            ltmp = .false.
            if (i > 1) then
                if (abs(arr(i,2) - arr(i-1, 2)) < degentol .and. &
                    (tAbelian .or. symeq(G1(i)%sym, G1(i-1)%sym))) then
                    ! We have the same degenerate rep as the previous entry
                    symreps(2, J) = symreps(2, J) + 1
                    lTmp = .true.
                endif
            endif
            if (.not. lTmp) then
                ! We have a new rep
                J = J + 1
                symreps(2, J) = 1
            endif
            SYMREPS(1,I)=J
         ENDDO
!         DO I=1,NBASIS
!            WRITE(6,*) "SR1",SYMREPS(1,I),SYMREPS(2,I)
!         ENDDO   
      END SUBROUTINE GENSYMREPS

!.  Irrep symmetries are specified in SYM(5).
!   if SYM(5)=0, we assume it's totally symmetric
!   Other irreps contributing to the symmetry have bits set in 
!   SYM.
!   e.g. if irreps are a1,a2,b1,b2
      LOGICAL FUNCTION LCHKSYM(ISYM,JSYM)
         use SystemData, only: BasisFN,Symmetry
         IMPLICIT NONE
         TYPE(BASISFN) ISYM,JSYM
         INTEGER I
         LCHKSYM=.TRUE.
         DO I=1,3
            IF(ISYM%K(I).NE.JSYM%K(I)) LCHKSYM=.FALSE.
         ENDDO
         IF(ISYM%Ms.NE.JSYM%Ms) LCHKSYM=.FALSE.
         IF(ISYM%Ml.NE.JSYM%Ml) LCHKSYM=.FALSE.
!   if the symmetry product of I and J doesn't contain the totally
!   symmetric irrep, we set sym to .FALSE.
        LCHKSYM=LCHKSYM.AND.LSYMSYM(SYMPROD(SymConj(ISYM%SYM),JSYM%SYM))
      RETURN
      END FUNCTION LCHKSYM
      
      LOGICAL FUNCTION LCHKSYMD(NI,NJ,NEL,G1,NBASISMAX)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BASISFN) ISYM,JSYM,G1(*)
         INTEGER NEL,NI(NEL),NJ(NEL),nBasisMax(5,*)
         CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
         CALL GETSYM(NJ,NEL,G1,NBASISMAX,JSYM)
         LCHKSYMD=LCHKSYM(ISYM,JSYM)
         RETURN
      END FUNCTION LCHKSYMD
!   NBASISMAX descriptor (1,3)
!
! HUBBARD:
!  BITS
!. 0 Tilted
!. 1 non-pbc
!. 2 real-space
!   which effects to values
!   MOM SPACE
! 0 Non-Tilted Lattice - pbc
! 1 Tilted Lattice - pbc
! 2 Non-Tilted lattice - no pbc
! 3 Tilted Lattice - no pbc
!   four following are REAL
! 4 Non-Tilted Lattice - pbc
! 5 Tilted Lattice - pbc
! 6 Non-Tilted lattice - no pbc
! 7 Tilted Lattice - no pbc
!
! (3,3)
! -2 Particle in a box
! -1 UEG
! 0 Hubbard
! 1 Generic spatial


!   This only works for momentum variables - 1-4
      SUBROUTINE ADDELECSYM(IEL,G1,NBASISMAX,ISYM)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BASISFN) ISYM,G1(*)
         INTEGER IELEC,nBasisMax(5,*)
         INTEGER I,IEL,SSYM
         CALL GETUNCSFELEC(IEL,IELEC,SSYM)
        IF(NBASISMAX(1,3).LT.4) THEN
!   Momentum space
            DO I=1,3
               ISYM%K(I)=ISYM%K(I)+G1(IELEC)%K(I)
            ENDDO
!   Symmetry space
         ELSEIF(NBASISMAX(3,3).EQ.0.AND.NBASISMAX(1,3).GE.4) THEN
!   We have no symmetries, so do nothing. (we're in real space)
!   except Ms
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
!   deal with momentum
            DO I=1,3
               ISYM%k(I)=ISYM%k(I)+G1(IELEC)%k(I)
            ENDDO
         ENDIF
         ISYM%MS=ISYM%MS+G1(IELEC)%MS
         ISYM%Ml=ISYM%Ml+G1(IELEC)%Ml
!   SSYM keeps track of the total S change on adding this electron
!   (it is +/-CSF_NSBASIS)
         I=ISYM%MS+0
         ISYM%Ms=I+SSYM
         RETURN
      END SUBROUTINE ADDELECSYM
      
      SUBROUTINE ROUNDSYM(ISYM,NBASISMAX)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BasisFN) ISYM
         INTEGER nBasisMax(5,*)
         INTEGER I
         IF(NBASISMAX(3,3).EQ.-2) THEN
!   particle in a box
!   parity symmetries
            DO I=1,3
               ISYM%k(I)=MOD(ISYM%k(I),2)
            ENDDO
         ELSEIF(NBASISMAX(3,3).EQ.-1) THEN
!   UEG (can't remember the symmetries of that
!   probably momentum  conservation)
         ELSEIF(NBASISMAX(3,3).EQ.0) THEN
!   Hubbard model
            IF(NBASISMAX(1,3).LT.2) THEN
!   momentum conservation - various PBC issues

!ALEX PLEASE CHECK.

                CALL MOMPBCSYM(ISYM%k,NBASISMAX)
!            ELSEIF(NBASISMAX(1,3).EQ.2) THEN
!   non-pbc mom space has parity symmetry
!               DO I=1,3
!                  ISYM(I)=MOD(ISYM(I),2)
!               ENDDO
            ELSEIF(NBASISMAX(1,3).GE.2) THEN
!   we're in real space so no sym
               DO I=1,3
                  ISYM%k(I)=0
               ENDDO
            ENDIF
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
!   Generic spatial symmetries
!           We need do nothing.
!   However, there is still momentum conservation - various PBC issues

!ALEX PLEASE CHECK.

                CALL MOMPBCSYM(ISYM%k,NBASISMAX)
         ENDIF
         RETURN 
      END SUBROUTINE ROUNDSYM

!   NBASISMAX descriptor (1,3)
!
! HUBBARD:
!  BITS
!. 0 Tilted
!. 1 non-pbc
!. 2 real-space
!   which effects to values
!   MOM SPACE
! 0 Non-Tilted Lattice - pbc
! 1 Tilted Lattice - pbc
! 2 Non-Tilted lattice - no pbc
! 3 Tilted Lattice - no pbc
!   four following are REAL
! 4 Non-Tilted Lattice - pbc
! 5 Tilted Lattice - pbc
! 6 Non-Tilted lattice - no pbc
! 7 Tilted Lattice - no pbc
!
! (3,3)
! -2 Particle in a box
! -1 UEG
! 0 Hubbard
! 1 Generic spatial

!   This function imposes periodic boundary conditions for the Hubbard Model
!   K1(3) is a K vector which is then mapped back into the Brillouin Zone of the Full cell.
!   The PBC are such that, a point, (KX,KY) (in terms of the cell lattice vectors), has bounds
!   -LENX/2 < KX <= LENX/2  and -LENY/2 <= KY < LENY/2

!   My ASCII art isn't up to drawing a picture unfortunately.

      SUBROUTINE MOMPBCSYM(K1,NBASISMAX)
!   NB the third column of NBASISMAX tells us whether it is tilted
         IMPLICIT NONE
         INTEGER K1(3),nBasisMax(5,*)
         INTEGER J,LDIM,AX,AY,LENX,LENY,KK2,T1,T2
         real(dp) R1,R2,NORM
!   (AX,AY) is the tilt of the lattice, and corresponds to the lattice vector of the cell.  The other lattice vector is (-AY,AX).
!These are expressed in terms of the primitive Hubbard lattice vectors
         AX=NBASISMAX(1,4)
         AY=NBASISMAX(2,4)
!   LENX is the length (i.e. number of lattice vectors) in direction (AX,AY).  LENY is in the other lattice vector direction
         LENX=NBASISMAX(1,5)
         LENY=NBASISMAX(2,5)
         IF(NBASISMAX(1,3).EQ.0.OR.NBASISMAX(1,3).EQ.0) THEN
!   A non-tilted lattice with PBC
            DO J=1,3
!    non-tilted
               KK2=K1(J)
               LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
               KK2=MOD(KK2,LDIM)
               IF(KK2.LT.NBASISMAX(J,1)) KK2=KK2+LDIM
               IF(KK2.GT.NBASISMAX(J,2)) KK2=KK2-LDIM
               K1(J)=KK2 
            ENDDO
         ELSEIF(NBASISMAX(1,3).EQ.1) THEN
!   we have a tilted lattice with PBC
!   we want the a1,a2 components of k
            NORM=AX*AX+AY*AY
            R1=(AX*K1(1)+AY*K1(2))/NORM
            R2=(AX*K1(2)-AY*K1(1))/NORM
            R1=R1/LENX+0.5_dp
            R2=R2/LENY+0.5_dp
!   R1 is now in terms of the shifted extended unit cell.  We want 0<R1<=1
!   R2 is now in terms of the shifted extended unit cell.  We want 0<=R2<1
!   T1 and T2 will be the vectors to translate the point K1 back into our cell
            T1=INT(ABS(R1))
            IF(R1.LT.0.0_dp) THEN
               T1=-T1
               IF (abs(t1 - r1) > 1.0e-12_dp) T1=T1-1
            ENDIF
!   Now T1= highest integer less than R1  (i.e. FLOOR)
!   We want to include R1=1, so we have one fewer translations in that case.
            IF(abs(r1 - t1) < 1.0e-12_dp) T1=T1-1
            T2=INT(ABS(R2))
            IF(R2.LT.0.0_dp) THEN
               T2=-T2
               IF(abs(t2 - r2) > 1.0e-12_dp) T2=T2-1
            ENDIF
!   The following line conserves the top-right rule for edge effects
            IF(abs(r1 - t1) < 1.0e-12_dp) T1=T1-1
!   Now T2= highest integer less than R2  (i.e. FLOOR)
!   Do the translation
            IF(R1.GT.1.0_dp.OR.R1.LE.0.0_dp) R1=R1-T1
            IF(R2.GE.1.0_dp.OR.R2.LT.0.0_dp) R2=R2-T2
!   Now convert back to the the Coords in the extended brillouin zone (ish)
            R1=(R1-0.5_dp)*LENX
            R2=(R2-0.5_dp)*LENY
!   Now K1 is defined as the re-scaled (R1,R2) vector
            K1(1)=NINT(R1*AX-R2*AY)
            K1(2)=NINT(R1*AY+R2*AX)
         ENDIF
         RETURN
      END SUBROUTINE MOMPBCSYM

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
      END FUNCTION SYMLT
      LOGICAL FUNCTION SYMNE(A,B)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) A,B
         SYMNE=A%s.NE.B%s
         RETURN
      END FUNCTION SYMNE
      LOGICAL FUNCTION SYMEQ(A,B)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) A,B
         SYMEQ=A%s.EQ.B%s
         RETURN
!Need to cope with 'unsigned integers'
      END FUNCTION SYMEQ
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
      END FUNCTION SYMGT

      integer Function FindSymLabel(s)
         use SystemData, only: Symmetry
         use SymData, only: SymLabels,nSymLabels

         IMPLICIT NONE
         Type(Symmetry) s
         integer i
         do i=1,nSymLabels
            if(symeq(SymLabels(i),s)) exit
         enddo
         if(i.gt.nSymLabels) i=0
         FindSymLabel=i
         return
      end Function FindSymLabel
!   A binary search to find VAL in TAB.  TAB is sorted
      SUBROUTINE BINARYSEARCHSYM(VAL,TAB,LEN,LOC)
         use SystemData, only: Symmetry
         IMPLICIT NONE
         TYPE(Symmetry) VAL
         INTEGER LOC,LEN
         type(Symmetry) TAB(LEN)
         INTEGER I,J,IFIRST,N,ILAST
         I=1
         J=LEN
         IFIRST=I
         ILAST=J
         DO WHILE(J-I.GE.1)
            N=(I+J)/2
!            WRITE(6,"(3I4)",advance='no') I,J,N
!            CALL WRITESYM(6,TAB(1,I),.FALSE.)
!            CALL WRITESYM(6,TAB(1,J),.FALSE.)
!            CALL WRITESYM(6,TAB(1,N),.TRUE.)
            IF(SYMLT(TAB(N),VAL).AND.I.NE.N) THEN
               IF(SYMNE(TAB(N),TAB(IFIRST))) IFIRST=N
!   reset the lower limit
               I=N
            ELSEIF(SYMGT(TAB(N),VAL)) THEN
               IF(SYMNE(TAB(N),TAB(ILAST))) ILAST=N
!   reset the upper limit
               J=N
            ELSEIF(SYMEQ(TAB(N),VAL)) THEN
!   bingo, we've got it!
               LOC=N
               RETURN
            ELSE
!   we've reached a situation where I and J's entries have the same value, and it's
!   not the one we want.  Leave the loop.
               I=J
            ENDIF
         ENDDO
         IF(SYMEQ(TAB(I),VAL)) THEN
            LOC=I
         ELSEIF(SYMEQ(TAB(J),VAL)) THEN
            LOC=J
         ELSE
!   Failure
            LOC=0
         ENDIF
      END SUBROUTINE BINARYSEARCHSYM

!This function when called repeatedly, generates all allowed symmmetries.
!It appears to not have been modified for Lz symmetry, so will only generate Lz=0
! AJWT 20110121
      SUBROUTINE GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,TPARITY,IPARITY,TSETUP,TDONE,IMAX,ISYM)
         use SystemData, only: Symmetry, BasisFN, NullBasisFn
         use SymData, only: tAbelian
         IMPLICIT NONE
         INTEGER NEL,nBasisMax(5,*)
         INTEGER LMS
         TYPE(BasisFN) IPARITY,ISYM,IMax(2)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE,KALLOWED,TMORE2
         INTEGER ILEV
         IF(TSETUP) THEN
            IMAX(1)=NullBasisFn
            IMAX(2)=NullBasisFn
            DO ILEV=1,3
               IF(TPARITY) THEN
                  IMAX(1)%k(iLev)=IPARITY%k(ILEV)
                  IMAX(2)%k(ILEV)=IPARITY%k(ILEV)
               ELSE
                  IMAX(1)%k(iLev)=NBASISMAX(ILEV,1)
                  IMAX(2)%k(iLev)=NBASISMAX(ILEV,2)
!                  IF(NBASISMAX(1,3).EQ.2) THEN
!   hubbard non-pbc mom space
!                     IMAX(ILEV,1)=IMAX(ILEV,1)*NEL
!                     IMAX(ILEV,2)=IMAX(ILEV,2)*NEL
!                  ENDIF
               ENDIF
            ENDDO
            IF(TSPN) THEN
               IMAX(1)%Ms=LMS
               IMAX(2)%Ms=LMS
            ELSE
               IMAX(1)%Ms=NBASISMAX(4,1)*NEL
               IMAX(2)%Ms=NBASISMAX(4,2)*NEL
            ENDIF
!   If we're specifying a sym (TPARITY) in IPARITY(5), and
!   we have a system with all 1D reducible orbs, then we put
!   that into IMAX
            IF(NBASISMAX(5,2).NE.0.OR.TAbelian) THEN
               IF(TPARITY) THEN
                  IMAX(1)%Sym%s=IPARITY%Sym%s
                  IMAX(2)%Sym%s=IMAX(1)%Sym%s
               ELSE
                  IMAX(1)%Sym%s=MinSymRep()
                  IMAX(2)%Sym%s=MaxSymRep()
               ENDIF
            ELSE
!   we've got a sym system with polydimensional irreps, which leads to
!   dets with combinations of irreps, so we cannot put sym into blocks
                
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
!   Go to the next sym.
         TMORE2=.TRUE.
         TMORE=.TRUE.
         ILEV=5
         DO WHILE(TMORE2)
            DO WHILE (ILEV.GT.0)
               IF(ILEV.EQ.5) THEN
                  IF(IMAX(1)%Sym%s.NE.0) THEN
!   symmetry specifiers are incremented by multiplying*2 (unless there are no syms counted)
                     ISYM%Sym%s=ISYM%Sym%s*2
                  ELSE
                     Call IncrSym(ISym%Sym)
                  ENDIF
                  IF(ISYM%Sym%s.EQ.IMAX(1)%Sym%s) THEN
                     ILEV=ILEV-1
                     IF(ILEV.EQ.0) THEN
                        TMORE2=.FALSE.
!   If we've run out of syms, we give up
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
!   If we've run out of syms, we give up
                        TMORE=.FALSE.
                     ENDIF
                  ELSEIF(ILEV.LT.4) THEN
!   We've just incremented one of the higher columns, now go down to the
!   lower ones.
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
      END SUBROUTINE GENNEXTSYM
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
!   always a spin symmetry
         IF(IMAX(1)%Ms.NE.IMAX(2)%Ms) IMAX(1)%Ms=0
      END SUBROUTINE DOSYMLIMDEGEN
      SUBROUTINE GETSYMDEGEN(ISYM,NBASISMAX,IDEGEN)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BasisFN) ISym,ISym2
         INTEGER nBasisMax(5,*),IDEGEN,I,J
         LOGICAL KALLOWED,TDO
         IDEGEN=0
         IF(NBASISMAX(3,3).EQ.0) THEN
!   Hubbard
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
!   Spin
         IF(ISYM%Ms.NE.0) IDEGEN=IDEGEN*2
      END SUBROUTINE GETSYMDEGEN

!   Initialize symmetry to take into account the core electrons
      SUBROUTINE SETUPSYM(ISYM)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: FrozenSym
         IMPLICIT NONE
         TYPE(BasisFN) ISym
         ISym=FrozenSym
         RETURN
      END SUBROUTINE SETUPSYM

      SUBROUTINE WRITEALLSYM(IUNIT,SYM,LTERM)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER IUNIT
         TYPE(BASISFN) SYM
         LOGICAL LTERM
         WRITE(IUNIT,"(4I5)",advance='no')SYM%K(1),SYM%K(2),SYM%K(3),SYM%MS
         CALL WRITESYM(IUNIT,SYM%SYM,LTERM)
      END SUBROUTINE WRITEALLSYM
      SUBROUTINE WRITESYM(IUNIT,SYM,LTERM)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: nSym,tAbelian, TwoCycleSymGens
         IMPLICIT NONE
         INTEGER IUNIT
         TYPE(SYMMETRY) SYM
         LOGICAL LTERM
         INTEGER Abel(3)
         IF(TAbelian) THEN
            CALL DecomposeAbelianSym(SYM%s,Abel)
            if (TwoCycleSymGens) then
              WRITE(IUNIT,'(" (",I2,",",I2,",",I2,")",I2)',advance='no') Abel(1:3),SYM%s
            else
              WRITE(IUNIT,'(" (",I2,",",I2,",",I2,")",I2)',advance='no') Abel(1:3)
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
      END SUBROUTINE WRITESYM
      SUBROUTINE SetupFreezeAllSym(Sym)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: FrozenSym
         IMPLICIT NONE
         TYPE(BasisFN) Sym
!   Set to be totally symmetric
         FrozenSym=Sym
      END SUBROUTINE SetupFreezeAllSym
      SUBROUTINE SetupFreezeSym(Sym)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: FrozenSym
         IMPLICIT NONE
         TYPE(BasisFN) Sym
!         FrozenSym=NullBasisFn
!   Set to be totally symmetric
         FrozenSym=Sym
      END SUBROUTINE SetupFreezeSym
 
!Deal with K-point symmetries, using translational symmetry operations.
      ! JSS: use Abelian symmetry formulation (allows us to go beyond 64
      ! symmetry operations, and hence deal with larger k-point meshes).
      SUBROUTINE GenKPtIrreps(nTranslat,nKps,KpntInd,nStates)
      use SystemData, only: Symmetry, BasisFN
      use SymData, only: SymLabels,SymLabelChars,nRot,nSymLabels,KPntSym
      use SymData, only: SymClasses,tagSymLabelChars,tagSymLabels
      use SymData, only: tagSymClasses
      use global_utilities
      IMPLICIT NONE
      INTEGER I,nStates
      INTEGER nTranslat,nKps,KpntInd(nStates)
      character(*), parameter :: this_routine='GenKPtIrreps'
      nSymLabels=nKps
      nRot=nTranslat
      allocate(SymLabelChars(nRot,nSymLabels))
      call LogMemAlloc('SymLabelChars',nSymLabels*nRot,16,this_routine,tagSymLabelChars)
      allocate(SymLabels(nSymLabels))
      call LogMemAlloc('SymLabels',nSymLabels,4,this_routine,tagSymLabels)
      allocate(SymClasses(nStates))
      call LogMemAlloc('SymClasses',nStates,4,this_routine,tagSymClasses)
      SYMLABELCHARS=0.0_dp
      DO I=1,nStates
        SymClasses(I)=KpntInd(I)
        SymLabels(KPntInd(I))%s=ComposeAbelianSym(KpntSym(:,KPntInd(I)))
      END DO
      write (6,*) 
      write(6,'(a11," |",a13,"|",a10)')' K-vector',' Label ','Conjugate'
      write (6,'(39("-"))')
      do i=1,nSymLabels
        write (6,'("(",3i3,")"," | ")',advance='no') KpntSym(:,I)
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
      END SUBROUTINE GenKPtIrreps

      subroutine  DecomposeAbelianSym(ISym,AbelSym)
      ! Store the symmetry index as integer(int64).  For Abelian symmetry
      ! we need to have 3 numbers stored in this index.  We store
      ! according to isym=\sum_i AbelSym(i)*32768**(i-1).
      ! This allows point groups with more than 64 irreps to be used in
      ! the point group is Abelian (as all translational/k-point
      ! symmetries are).
      ! Decompose the symmetry label back into the appropriate "quantum"
      ! numbers.
      ! Store the symmetry index as integer(int64).  For Abelian symmetry
      ! we need to have 3 numbers stored in this index.  We store
      ! according to isym=1+\sum_i AbelSym(i)*32768**(i-1).
      ! Note that many symmetry parameters for CPMD-NECI jobs are set in
      ! kpntrep.F in CPMD source.
      ! Decompose the symmetry label back into the appropriate
      ! numbers...
      use SystemData, only: Symmetry, BasisFN
      use SymData, only: PropBitLen
      implicit none
      integer(int64) Isym
      integer  AbelSym(3)
!RShift
      AbelSym(3)=int(IShft(Isym,-(PropBitLen*2)),sizeof_int)
!RShift
      AbelSym(2)=int(Iand(IShft(ISym,-PropBitLen),2_int64**PropBitLen-1),sizeof_int)
      AbelSym(1)=int(Iand(Isym,2_int64**PropBitLen-1),sizeof_int)
      return
      end subroutine DecomposeAbelianSym

      integer(int64) function ComposeAbelianSym(AbelSym)
          use SystemData, only: Symmetry, BasisFN
          use SymData, only: PropBitLen
          implicit none
          integer  AbelSym(3)
          integer(int64) TempVar
          TempVar=AbelSym(3)
!LShift
          ComposeAbelianSym=IShft(Tempvar,PropBitLen*2)    &
                            +IShft(AbelSym(2),PropBitLen)  &
                            +AbelSym(1)
      end function ComposeAbelianSym


      function TotSymRep()
          ! Our definition of the totally symmetric representation
          ! changes according to whether we're using Abelian/k-point
          ! symmetry or the standard symmetry.  It's just a matter of
          ! convenience, rather than some deep theoretical insight!
          use SystemData, only: Symmetry, BasisFN,tUEG,treal
          use SymData, only: tAbelian
          implicit none
          Type(Symmetry) TotSymRep 
          if (TAbelian.or.tUEG.or.treal) then
              TotSymRep%s=0
          else
              TotSymRep%s=1
          end if
      end function TotSymRep

      ! nBasisMax might well be needed in the future in these functions.
      integer(int64) function MinSymRep()
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: tAbelian
         implicit none
         if(TAbelian) then
            MinSymRep=0
         else
            MinSymRep=0
         endif
      end function MinSymRep
      integer(int64) function MaxSymRep()
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: tAbelian,nProp
         implicit none
         integer abel(3)
         abel(:)=nprop(:)-1
         if(TAbelian) then
            MaxSymRep=ComposeAbelianSym(abel)
         else
            MaxSymRep=0
         endif
      end function MaxSymRep
      subroutine IncrSym(Sym)
         use SystemData, only: Symmetry, BasisFN
         use SymData, only: tAbelian,nProp
         implicit none
         type(Symmetry) Sym
         integer abel(3),i
         logical lcont
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
        use SystemData, only: Symmetry, BasisFN
        use SymData, only: SymLabelCounts,SymLabelCountsCum
        use SymData, only: SymLabelIntsCum
        use SymData, only: tagSymLabelIntsCum,tagSymLabelCountsCum
        use global_utilities
        implicit none
        integer Nirrep,nBasis,iSS,nBi,i,basirrep,t
        integer(int64) iSize
        character(*), parameter :: this_routine='GetSymTMATSize'
        nBi=nBasis/iSS
        iSize=0
        allocate(SymLabelIntsCum(nIrrep))
        call LogMemAlloc('SymLabelIntsCum',nIrrep,4,this_routine,tagSymLabelIntsCum)
        allocate(SymLabelCountsCum(nIrrep))
        call LogMemAlloc('SymLabelCountsCum',nIrrep,4,this_routine,tagSymLabelCountsCum)
        SYMLABELINTSCUM(1:Nirrep)=0
        SYMLABELCOUNTSCUM(1:Nirrep)=0
        do i=1,Nirrep
            basirrep=SYMLABELCOUNTS(2,i)
            iSize=iSize+(basirrep*(basirrep+1))/2
            SYMLABELINTSCUM(i)=int(iSize,sizeof_int)
            IF(i.eq.1) THEN
                SYMLABELCOUNTSCUM(i)=0
            ELSE
                DO t=1,(i-1)
                    SYMLABELCOUNTSCUM(i)=SYMLABELCOUNTSCUM(i)+SYMLABELCOUNTS(2,t)
                ENDDO
            ENDIF
            write(6,*) basirrep,SYMLABELINTSCUM(i),SYMLABELCOUNTSCUM(i)
            call neci_flush(6)
        enddo
        iSize=iSize+2
        !This is to allow the index of '-1' in the array to give a zero value
      END SUBROUTINE GETSYMTMATSIZE


    ! This function returns the label (0 -> nSymlabels-1) of the symmetry
    ! product of two symmetry labels.
    PURE INTEGER FUNCTION RandExcitSymLabelProd(SymLabel1,SymLabel2)
        IMPLICIT NONE
        INTEGER , INTENT(IN) :: SymLabel1,SymLabel2

        IF(tNoSymGenRandExcits) THEN
            RandExcitSymLabelProd=0
        ELSEIF(tKPntSym) THEN
            !Look up the symmetry in the product table for labels (returning labels, not syms)
            RandExcitSymLabelProd=SymTableLabels(SymLabel1,SymLabel2)
        ELSE
            RandExcitSymLabelProd=IEOR(SymLabel1,SymLabel2)
!            WRITE(6,*) "***",SymLabel1,SymLabel2,RandExcitSymLabelProd
        ENDIF

    END FUNCTION RandExcitSymLabelProd

end module sym_mod
