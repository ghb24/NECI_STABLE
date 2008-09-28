
!! NB SYMSETUPEXCITS uses a store which has INTEGER*8 to cope with 64-bit machines, and passes this directly to SYMSETUPEXCITS2
!!  Iterators which use GENSYMEXCITIT2 use an INTEGER STORE, and create an INTEGER*8 STORE2 internally, and pass this to SYMSETUPEXCITS2

!! SYMSETUPEXCITS takes an INTEGER*8 STORE(6) which stores the pointers used by SYMSETUPEXCITS2
!!   This STORE can then be passed to SYMGENEXCITS.
C.. Setup for symmetry routine below.
      SUBROUTINE SYMSETUPEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,TCOUNT,
     &      ICOUNT,ILEVEL,iMinElec1,iMaxElec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER G1(nBasis),NBASISMAX(*)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER*8 STORE(6)
         INTEGER ICOUNT
         LOGICAL TCOUNT
         INTEGER ILEVEL
         INTEGER iMinElec1,iMaxElec1
         CALL IAZZERO(STORE,6*2)
         CALL MEMORY(IP_DSTORE,SymClassSize*NEL+NBASIS/32+1
     &      +SymmetrySize*(NEL*NEL+1),
     &      'DSTORE')
         STORE(1)=IP_DSTORE
         CALL SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,TCOUNT,
     &      ICOUNT,DSTORE(1),
     &      DSTORE(SymClassSize*NEL+1),
     &      DSTORE(SymClassSize*NEL+1+NBASIS/32+1),ILEVEL,iMinElec1,
     &      iMaxElec1)
C.. If we're just counting, we don't need to keep DSTORE
         IF(TCOUNT) CALL FREEM(IP_DSTORE)
      END


!=  ILEVEL is 1 for singles, 2 for just doubles, and 3 for both.
!=  A a new symmetry excitation generation algorithm.
!=  Algorithm: First list all the different symmetry classes in NI.
!=  For each pair of symmetry classes, determine its symmetry product.
!=  Classify and store each possible pair of orbitals under its symmetry product, [].
!=  With each of the symmetry products [], calculate []' such that []x[]'
!=  contains the totally symmetric rep.
!=  We do this by checking whether any of the []' in the 
!=  global symprods table multiply by [] to give the symmetric rep.

!= Excitation Generator Scaling   Timing                    Memory
!= 
!= Create Class List              O(N)                      4 N
!=   Gives NCL =O(NSYM)
!= Create Class Remainder info    O(NCL)                    
!= Create Class Sym Product info  O(NCL^2)                  4 NPR 
!=   Gives NPR ~= NCL^2
!==   for abelian  NPR<=NSYM
!= Save list of occ orbital pairs O(N^2) /NPROC             2 O(N^2) /NPROC
!=   This is sorted according
!=    to the symmetry product
!=                                                                   For abelian, this isn't
!=                                  [ O(NCL^2) ]          [ NSYM ]   needed as we can very easily
!=                                                                   determine this from a list of 
!=                                                                   occs sorted by sym.  We should
!=                                                                   store counts however
!= Count number of allowed virt
!=    pairs of each sym prod      O(NPR M^2 / NCL^2)        O(NCL^2)
!=    Also counts num of doubles
!=                                                                   For abelian this is simply
!=                                  [ NSYM ]              [ NSYM ]   the difference between the num
!=                                                                   of occ pairs of a given sym and
!=                                                                   the total num of pairs.
!= Calculate the number of
!=    single excitations          O(NCL^3)
!=                                  [ NSYM ]
!= --- Exit here for init
!= 
!= Store singles classes          O(NCL^3)                  O(NCL^2)
!=                                  [ NSYM ]                 [ NSYM ]
!= Store doubles classes          O(NPR NSYM)               O(NPR^2)
!=                                  [ NSYM^2 ]               [ NSYM^2 ] 

      SUBROUTINE SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &   TCOUNT,ICOUNT,CLASSES,ILUT,SYMPRODS,ILEVEL,iMinElec1,iMaxElec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClass,SymClasses,SymLabelCounts
         use SymData, only: SymStatePairs,SymLabels,SymPairProds
         use SymData, only: nSymPairProds,nSymLabels
         use global_utilities
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL),NBASIS
         TYPE(BasisFN) G1(nBasis)
         INTEGER nBasisMax(5,*)
         LOGICAL TCOUNT
         INTEGER*8 STORE(*)


         TYPE(SymClass) CLASSES(*)
         INTEGER CLASSCOUNT(2,NEL)
         INTEGER THISCLASSCOUNT(2,NEL)
!  ThisClassCount is used to list only electrons which this processor deals with
         INTEGER PREVCLASSCOUNT(2,NEL)
!  PrevClassCount is used to list electrons which lower indexed processors deal with
         TYPE(Symmetry) SYMPRODS(0:NEL*NEL)
         TYPE(Symmetry) PR,PR2
         INTEGER NPR
         INTEGER SYMPRODCOUNT(3,0:NEL*NEL)
         INTEGER NPAIRS,ICC1,ICC2,ICC3,ICC,ICC4
         INTEGER I,J,NCL,K,L
         LOGICAL TNEW
         INTEGER ORBPAIRS(2,*)
         POINTER (IP_ORBPAIRS,ORBPAIRS)
         INTEGER SYMPRODIND(2,3,0:*)
         POINTER (IP_SYMPRODIND,SYMPRODIND)
         INTEGER NALLOWPPS(3,*)
         TYPE(Symmetry) SPP
         POINTER(IP_NALLOWPPS,NALLOWPPS)
         LOGICAL L1B,L1A,L2B,L2A
         INTEGER NEXCITS,NSING,NDOUB,NEXCITTYPES

         INTEGER EXCITTYPES(5,*)
         POINTER (IP_EXCITTYPES,EXCITTYPES)

         INTEGER ICOUNT
         TYPE(Symmetry) SYMPROD,SymConj
         LOGICAL LSYMSYM,SYMEQ
         INTEGER ISPN
         INTEGER ILUT(0:NBASIS/32)
         INTEGER,SAVE :: ISUB=0
         INTEGER ILEVEL
         Type(Symmetry) TotSymRep

         LOGICAL ISUHFDET

         INTEGER iMinElec1,iMaxElec1
         type(timer), save :: proc_timer
         proc_timer%timer_name='SYMSUEXCIT'
         
         call set_timer(proc_timer,65)

C.. List the symmetry classes in NI.
C.. For each orb, just see if its class is in the list, and if not add it
C.. Each state's symmetry falls into a class SymClasses(ISTATE).

C.. We divide the occupied electrons into classes.  A class correpsonds to a unique value of SymClasses(iState) for state (i.e. spatial orbital) iState.
C.. The CLASSCOUNT arrays contain the number of electrons of each class.
C.. We keep three lists:   CLASSCOUNT  which takes into account all electrons.  
C..                     THISCLASSCOUNT which takes into account electrons allocated to this processor
C..  and                PREVCLASSCOUNT which takes into account electrons allocated to the previous processor.
C.. These both have the format (Spn,iCl), where Spin=1,2 corresponding to alpha and beta and
C..                                             iCl is the index of the Class.
C.. CLASSES(iCl)%SymLab is the symmetry label (i.e. entry in SymClasses) of this class (NB not a symmetry itself)

         NCL=0
         DO I=1,NEL
            TNEW=.TRUE.
!            WRITE(6,*) I,SymClasses((NI(I)+1)/2)
            DO J=1,NCL
               IF(SymClasses((NI(I)+1)/2).EQ.CLASSES(J)%SymLab) THEN
                  ISPN=(G1(NI(I))%Ms+3)/2
                  CLASSCOUNT(ISPN,J)=CLASSCOUNT(ISPN,J)+1
                  IF(I.GE.iMinElec1.AND.i.LE.iMaxElec1) THEN
                     THISCLASSCOUNT(ISPN,J)=THISCLASSCOUNT(ISPN,J)+1
                  ENDIF
                  IF(I.LT.iMinElec1) THEN
                     PREVCLASSCOUNT(ISPN,J)=PREVCLASSCOUNT(ISPN,J)+1
                  ENDIF
                  TNEW=.FALSE.
                  EXIT
               ENDIF
            ENDDO
            IF(TNEW) THEN
C.. add the new class
C.. ISPN=1 is beta, ISPN=2 is alpha
               ISPN=(G1(NI(I))%Ms+3)/2
               NCL=NCL+1
               CLASSES(NCL)%SymLab=SymClasses((NI(I)+1)/2)
               CLASSCOUNT(ISPN,NCL)=1
               CLASSCOUNT(3-ISPN,NCL)=0
               IF(I.GE.iMinElec1.AND.i.LE.iMaxElec1) THEN
                  THISCLASSCOUNT(ISPN,NCL)=1
                  THISCLASSCOUNT(3-ISPN,NCL)=0
               ELSE
                  THISCLASSCOUNT(1:2,NCL)=0
               ENDIF
               IF(I.LT.iMinElec1) THEN
                  PREVCLASSCOUNT(ISPN,NCL)=1
                  PREVCLASSCOUNT(3-ISPN,NCL)=0
               ELSE
                  PREVCLASSCOUNT(1:2,NCL)=0
               ENDIF
            ENDIF
         ENDDO
C.. Now calculate the resultant symmetry from removing one of each of the classes (for singles)
C.. i.e. if we remove j, what's the symmetry of SYMCONJ(<ikl... |) (X) SYM( |ikl...>)
C..   NB to cope with non-abelians, this actually takes the prod sym of the classes of each orb, 
C..   rather than the orb itself, so may end up classifying a single as allowed when it's not.
C.. we'll use this in the future to determina the allowed j->a
C.. This is stored in CLASSES(iCl)%SymRem
         DO I=1,NCL
            PR=TotSymRep()
            DO J=1,NCL
               L=CLASSCOUNT(1,I)+CLASSCOUNT(2,I)
               IF(J.EQ.I) L=L-1
               DO K=1,L
                  PR=SYMPROD( PR,SymConj(SymLabels(CLASSES(J)%SymLab)))
                  PR=SYMPROD( PR,SymLabels(CLASSES(J)%SymLab))
               ENDDO
            ENDDO
            CLASSES(I)%SymRem=PR
!            WRITE(6,"(A,I4,2Z4,I)") "CLASS ",NCL,CLASSES(I)%SymLab,
!     &            CLASSES(I)%SymRem,I
         ENDDO
!            WRITE(6,*) NCL," Symmetry Classes"

C.. Now determine all possible symmetry products of pairs of orbs, and classify them.
C.. SYMPRODS(iPr) is the a symmetry product.
C.. SymProdCount(iType,iPr) is the number of prods of a given type.  iType=0 (beta/beta), 1 (alpha/beta), 2 (alpha/alpha)
C.. NPairs totals all SymProdCounts.

C.. For each pair of symmetry classes, determine its symmetry product.
C.. We list them and then sort them
         CALL IAZZERO(SYMPRODCOUNT,(NEL*NEL+1)*3)
         NPR=0
         NPAIRS=0
         DO I=1,NCL
            DO J=I,NCL
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SymLabels(ISYMLABEL)
               PR=SYMPROD( SymLabels(CLASSES(I)%SymLab),
     &                     SymLabels(CLASSES(J)%SymLab))
               TNEW=.TRUE.
               DO K=1,NPR
                  IF(SYMEQ(SYMPRODS(K),PR)) THEN
                     TNEW=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF(TNEW) THEN
                  NPR=NPR+1
                  K=NPR
                  SYMPRODS(NPR)=PR
               ENDIF
C.. 1 corresponds to beta beta, 2 to alpha/beta and beta/alpha, and 3 to alpha, alpha
C.. 1 corresponds to beta beta, 2 to alpha/beta and beta/alpha, and 3 to alpha, alpha
!  The number of pairs is somewhat more complex when we split electron 1 between processors
!  For systems where we require the upper-triangle (excluding diagonal) we have
!  NThis*(NAll+(NAll-NThis)-1)/2  pairs
!  For cases where we allow all pairs, it's just NThis*NAll
               IF(I.EQ.J) THEN
                  ICC1=THISCLASSCOUNT(1,I)*
     &               (2*CLASSCOUNT(1,I)-THISCLASSCOUNT(1,I)
     &                  -2*PREVCLASSCOUNT(1,I)-1)/2
!If we have all electrons being considered here then the alpha beta lot is simple, and just Nalpha*Nbeta
!                  ICC2=CLASSCOUNT(1,I)*CLASSCOUNT(2,I)
!  However, if we are splitting electron 1 between processors, alpha electrons here can only interact with beta
!    electrons with higher index, and the same witb beta on alpha
!  Therefore we have NThisA*(NAllB-NPrevB)+ NThisB*(NAllA-NPrevA-NThisA)
                  ICC2=THISCLASSCOUNT(1,I)*
     &               (CLASSCOUNT(2,I)-PREVCLASSCOUNT(2,I))
     &               + THISCLASSCOUNT(2,I)*
     &         (CLASSCOUNT(1,I)-PREVCLASSCOUNT(1,I)-THISCLASSCOUNT(1,I))

                  ICC3=THISCLASSCOUNT(2,I)*
     &               (2*CLASSCOUNT(2,I)-THISCLASSCOUNT(2,I)
     &               -2*PREVCLASSCOUNT(2,I)-1)/2
               ELSE
                  ICC1=THISCLASSCOUNT(1,I)*CLASSCOUNT(1,J)
                  ICC2=THISCLASSCOUNT(1,I)*CLASSCOUNT(2,J)
     &                +THISCLASSCOUNT(2,I)*CLASSCOUNT(1,J)
                  ICC3=THISCLASSCOUNT(2,I)*CLASSCOUNT(2,J)
               ENDIF
               SYMPRODCOUNT(1,K)=SYMPRODCOUNT(1,K)+ICC1
               SYMPRODCOUNT(2,K)=SYMPRODCOUNT(2,K)+ICC2
               SYMPRODCOUNT(3,K)=SYMPRODCOUNT(3,K)+ICC3
     
               NPAIRS=NPAIRS+ICC1+ICC2+ICC3
!                  WRITE(6,"(6I6)") CLASSCOUNT(1:2,I),
!     &  THISCLASSCOUNT(1:2,I), PREVCLASSCOUNT(1:2,I)
!                  WRITE(6,"(6I6)") CLASSCOUNT(1:2,J),
!     &  THISCLASSCOUNT(1:2,J), PREVCLASSCOUNT(1:2,J)
!                  WRITE(6,"(A,I6,Z4,5I6)") "SYMPROD ",K,SYMPRODS(K),I,J,
!     &               ICC1,ICC2,ICC3
            ENDDO
         ENDDO
         SYMPRODS(0)%s=0 ! For historical reasons.  Should never be needed.
         SYMPRODCOUNT(1,0)=0
         SYMPRODCOUNT(2,0)=0
         SYMPRODCOUNT(3,0)=0
C         DO I=1,NPR
C            WRITE(6,"(Z4,3I4)") SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
C         ENDDO

C.  We sort the array of Symmetry products (and its associated SymProdCount info), into ascending order of product, for easy searching later.
         CALL SORTSYMIARR(NPR,SymProds(1),SymmetrySize,
     &                                             SymProdCount(:,1),3)
!         WRITE(6,*) NPR," Symmetry Products"
!         WRITE(6,*) NPAIRS," Orbital Pairs"
C.. Create a cumulative sum.
         DO I=1,NPR
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)+SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)+SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)+SYMPRODCOUNT(3,I)
C            WRITE(6,*) SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
         ENDDO

C..  SymProdCount(iType,iPr) is the number of pairs of orbs up to and including this (iType, iPr) set.

C.. Allocate enough memory to store the index
         IF(STORE(5).EQ.0) THEN
            CALL MEMORY(IP_SYMPRODIND,2*3*(NPR+1),'SYMPRODIND')   
         ELSE
C.. we've acually got space to store it somewhere
            IP_SYMPRODIND=STORE(5)
         ENDIF
         CALL IAZZERO(SYMPRODIND,2*3*(NPR+1))

C.. Now shift this such that SYMPRODCOUNT(ISPN,I) is the index of
C.. the first orbital pair of SYMPROD(I) with spin ISPN in ORBPAIRS
C.. Store in SYMPRODIND(1,ISPN,I) too
         DO I=NPR,1,-1
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)
            SYMPRODIND(1,3,I)=SYMPRODCOUNT(3,I)
            SYMPRODIND(1,2,I)=SYMPRODCOUNT(2,I)
            SYMPRODIND(1,1,I)=SYMPRODCOUNT(1,I)
         ENDDO
C         DO I=1,NPR
C            WRITE(6,"(Z4,3I4)") SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
C         ENDDO
C.. We allocate enough memory to store all the pairs.
C.. Each pair consists of (ORB1,ORB2) where ORB1<ORB2
         IF(STORE(4).EQ.0) THEN
            CALL MEMORY(IP_ORBPAIRS,2*NPAIRS/IRAT+1,'ORBPAIRS')
         ELSE
            IP_ORBPAIRS=STORE(4)
         ENDIF

C.. Classify and store each possible pair of orbitals under its symmetry product.
C.. With each symmetry product (denoted []), calculate []'' which is the 
C.. remainder after having removed the orbitals.


C.. Now classify each pair of orbs under its sym prod
C.. SYMPRODS(I) specifies a sym prod.  Orbitals for this are stored in
C.. ORBPAIRS starting at index SYMPRODIND(1,ISPN,I)+1.
C.. The first free space is at SYMPRODIND(1,ISPN,I)+SYMPRODIND(2,ISPN,I)+1.
         CALL IAZZERO(ORBPAIRS,2*NPAIRS)
         L=0
         DO I=iMinElec1,iMaxElec1
            
            DO J=I+1,NEL
               L=L+1
               PR=SYMPROD(G1(NI(I))%Sym,G1(NI(J))%Sym)
C.. Find the product in the sorted symprod list
               CALL BINARYSEARCHSYM(PR,SYMPRODS(1),SymmetrySize,NPR,K)
               IF(K.EQ.0) THEN
                  WRITE(6,*) "Occupied Symmetry Products"
                  DO L=0,NPR
                     WRITE(6,"(I4)",advance='no') L
                     CALL WRITESYM(6,SYMPRODS(L),.TRUE.)
                  ENDDO
                  WRITE(6,"(A)",advance='no') "Illegal Symmetry"
                  CALL WRITESYM(6,PR,.TRUE.)
                  STOP 'Illegal Symmetry Found'
               ENDIF
               ISPN=(G1(NI(I))%Ms+G1(NI(J))%Ms)/2+2

C.. Increment the offset for the next one
               SYMPRODIND(2,ISPN,K)=SYMPRODIND(2,ISPN,K)+1
               ICC=SYMPRODIND(1,ISPN,K)+SYMPRODIND(2,ISPN,K)
C..                 The index          +  offset of the element

               ORBPAIRS(1,ICC)=NI(I)
               ORBPAIRS(2,ICC)=NI(J)
               
!               WRITE(6,"(A,2I4,2Z4,2I4)") 
!     &            "SP",ICC,ISPN,PR,SYMPRODS(K),NI(I),NI(J)
C               WRITE(6,*) L,NI(I),NI(J),PR,K
            ENDDO
         ENDDO

C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements

C.. With each of the symmetry products, calculate []' such that []x[]'
C.. contains the totally symmetric rep.
C.. We do this by checking whether any of the []' in the 
C.. global symprods table multiplied by [] give the sym rep.


C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) is a global list of all SYMPRODs available, the number of pairs of
C.. symlabels (listed in SymStatePairs), and the index of the start of this list
C.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
C.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of symlabels whose prod is of that symmetry.

C.. Now go through the list of all pairs, finding out how many pairs are to be excluded as they contain some of the
C.. orbitals in this determinant.
C.. We first create a quick lookup table to enable us to quickly check whether a given orbital is in this determinant
C.. in order 1, not order NEL.
         CALL IAZZERO(ILUT,NBASIS/32+1)
         DO I=1,NEL
            ILUT((NI(I)-1)/32)=IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32))
         ENDDO
!         DO I=0,NBASIS/32
!            WRITE(6,"(A,Z10)") "LUT: ",ILUT(I)
!         ENDDO
C.. Now look through the list of our pairs.  For each pair sym of the complete list which has a 
C.. symmetric product with any of our pair syms, we work out how many allowed pairs there are in
C.. the complete list, and store that value in NALLOWPPS
         IF(STORE(3).EQ.0) THEN
            CALL MEMORY(IP_NALLOWPPS,3*NSYMPAIRPRODS/IRAT+1,'NALLOWPPS')
         ELSE
            IP_NALLOWPPS=STORE(3)
         ENDIF
         NDOUB=0
         NEXCITTYPES=0
         DO I=1,NPR
            DO J=1,NSYMPAIRPRODS
C.. If it's the first time around, we flag NALLOWPPS as invalid
               IF(I.EQ.1) THEN
                  NALLOWPPS(1,J)=-1
                  NALLOWPPS(2,J)=-1
                  NALLOWPPS(3,J)=-1
               ENDIF
C.. [] is the SymProd(I).  For each sympairprod, check if [] x []' contains the sym rep
C               WRITE(6,*) SYMPRODS(I),SYMPAIRPRODS(J)
               SPP=SymPairProds(J)%Sym
               IF(LSYMSYM(SYMPROD(SymConj(SYMPRODS(I)), SPP))) THEN
C.. Check if we've worked out the allowed number of excitation pairs, otherwise work it out
                  IF(NALLOWPPS(1,J).EQ.-1) THEN
C.. Zero the counters
                     NALLOWPPS(1,J)=0
                     NALLOWPPS(2,J)=0
                     NALLOWPPS(3,J)=0

                     DO K=SymPairProds(J)%nIndex,
     &                 SymPairProds(J)%nIndex+SymPairProds(J)%nPairs-1

C.. Now check according to ISPN
C.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
C.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
                        ICC1=SymStatePairs(1,K)*2-1
                        ICC2=ICC1+1
                        ICC3=SymStatePairs(2,K)*2-1
                        ICC4=ICC3+1
                        L1B=.NOT.BTEST(ILUT((ICC1-1)/32),
     &                        MOD(ICC1-1,32))
                        L1A=.NOT.BTEST(ILUT((ICC2-1)/32),
     &                        MOD(ICC2-1,32))
                        L2B=.NOT.BTEST(ILUT((ICC3-1)/32),
     &                        MOD(ICC3-1,32))
                        L2A=.NOT.BTEST(ILUT((ICC4-1)/32),
     &                        MOD(ICC4-1,32))
C.. L1B is set if the beta of the first virtual is not in NI, i.e. is allowed
!  NALLOWPPS is the number of ALLOWED virtuals
                        IF(ICC1.EQ.ICC3) THEN
!                           NALLOWPPS(1,J)=NALLOWPPS(1,J)+1
!                           NALLOWPPS(3,J)=NALLOWPPS(3,J)+1
C.. The virtuals we are exciting to are the same state.  We're only allowed an AB excitation
                         IF(L1B.AND.L2A) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                        ELSE
                         IF(L1B.AND.L2B) NALLOWPPS(1,J)=NALLOWPPS(1,J)+1
                         IF(L1B.AND.L2A) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                         IF(L1A.AND.L2B) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                         IF(L1A.AND.L2A) NALLOWPPS(3,J)=NALLOWPPS(3,J)+1
                        ENDIF
C                        WRITE(6,*) SPP,SymStatePairs(1,K)*2-1,
C     &                     SymStatePairs(2,K)*2-1
C     &                     ,L1B,L1A,L2B,L2A
                     ENDDO
                  ENDIF
!.. Now count all the excitations allowed for this set of []->[]'
                  DO ISPN=1,3
                     K=NALLOWPPS(ISPN,J)
                     L=SYMPRODIND(2,ISPN,I)
                     ICC=L*K
C                     WRITE(6,"(2I5,2Z6,4I5)") I,J,SYMPRODS(I),
C     &                  SYMPAIRPRODS(J),ISPN,L,K,ICC
                     NDOUB=NDOUB+ICC
                     IF(ICC.GT.0) NEXCITTYPES=NEXCITTYPES+1
                  ENDDO
               ENDIF
            ENDDO
         ENDDO       

!         WRITE(6,*) "Number of double excitations: ",NDOUB


C.. Now calculate the number of single excitations

C.. SYMLABELLIST holds a list of states grouped under symmlabel
C.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
C.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I

         NSING=0
         IF(.NOT.ISUHFDET(NI,NEL)) THEN
          DO I=1,NCL
C.. For each class, see what sym prods it can interact with
            DO J=1,NSYMLABELS
               IF(LSYMSYM(SYMPROD(CLASSES(I)%SymRem,
     &                  SYMPROD(SymConj(SymLabels(J)),
     &                         SymLabels(CLASSES(I)%SymLab))))) THEN
                  DO ISPN=1,2
                     ICC=0
! Work out which virtuals are available
                     DO K=1,NCL
                       IF(CLASSES(K)%SymLab.EQ.J)
     &                  ICC=CLASSCOUNT(ISPN,K)
                     ENDDO
                     ICC=THISCLASSCOUNT(ISPN,I)*
     &                  (SYMLABELCOUNTS(2,J)-ICC)
                     NSING=NSING+ICC
                     IF(ICC.NE.0) NEXCITTYPES=NEXCITTYPES+1
C                     WRITE(6,"(2Z4,4I)") SymLabels(CLASSES(I)),
C     &                  SymLabels(J),ISPN,ICC,
C     &                  CLASSCOUNT(ISPN,I),
C     &                  (SYMLABELCOUNTS(2,J)-CLASSCOUNT(ISPN,I))
                  ENDDO
               ENDIF
            ENDDO
          ENDDO         
         ENDIF

!         WRITE(6,*) "Number of single excitations: ",NSING
         NEXCITS=0
         IF(BTEST(ILEVEL,0)) NEXCITS=NEXCITS+NSING
         IF(BTEST(ILEVEL,1)) NEXCITS=NEXCITS+NDOUB
         ICOUNT=NEXCITS
C         WRITE(6,*) "Total number of excitations: ",NEXCITS
C         WRITE(6,*) "Total number of excitation types: ",NEXCITTYPES

         IF(TCOUNT) THEN
C.. If we're just counting, we're done, so we get rid of some pointers.
            IP_EXCITTYPES=0
C.. However, we do save the length of the memory required.
C.. EXCITTYPES
            STORE(2)=(NEXCITTYPES*5)
C.. NALLOWPPS
            STORE(3)=3*NSYMPAIRPRODS
C.. ORBPAIRS
            STORE(4)=2*NPAIRS
C.. SYMPRODIND
            STORE(5)=2*3*(NPR+1)
C.. indicate that these are lengths
            STORE(6)=0
            
            CALL FREEM(IP_NALLOWPPS)
            CALL FREEM(IP_ORBPAIRS)
            CALL FREEM(IP_SYMPRODIND)
         ELSE
C.. Now allocate memory to store all the excitation types if there hasn't been one already allocated.
C.. This will have to be manually deallocated later.
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
         IF(STORE(2).EQ.0) THEN
          CALL MEMORY(IP_EXCITTYPES,(NEXCITTYPES*5)/IRAT+1,'EXCITTYPES')
         ELSE
            IP_EXCITTYPES=STORE(2)
         ENDIF

C.. Now store the singles:
         NEXCITTYPES=0
C         CALL WRITEDET(6,NI,NEL,.TRUE.)
         IF(.NOT.ISUHFDET(NI,NEL)) THEN
          IF(BTEST(ILEVEL,0)) THEN
            DO I=1,NCL
C.. For each class, see what sym prods it can interact with
               DO J=1,NSYMLABELS

                  PR=SYMPROD(SymConj(SymLabels(J)),
     &                       SymLabels(CLASSES(I)%SymLab))
                  PR2=SYMPROD(PR,Classes(I)%SymRem)
C                  WRITE(6,"(2I4,5Z4)") I,J,
C     &                        CLASSES(I)%SymRem,
C     &                         SymLabels(J),
C     &                         SymLabels(CLASSES(I)%SymLab),
C     &                         PR,PR2
                  IF(LSYMSYM(PR2)) THEN
                     DO ISPN=1,2
                        ICC=0
                        DO K=1,NCL
                           IF(CLASSES(K)%SymLab.EQ.J) 
     &                        ICC=THISCLASSCOUNT(ISPN,K)
                        ENDDO
                        ICC=THISCLASSCOUNT(ISPN,I)*
     &                     (SYMLABELCOUNTS(2,J)-ICC)
                        IF(ICC.NE.0) THEN
                           NEXCITTYPES=NEXCITTYPES+1
                           EXCITTYPES(1,NEXCITTYPES)=1
                           EXCITTYPES(2,NEXCITTYPES)=ISPN
                           EXCITTYPES(3,NEXCITTYPES)=I
                           EXCITTYPES(4,NEXCITTYPES)=J
                           EXCITTYPES(5,NEXCITTYPES)=ICC
C                           WRITE(6,"(A,I,Z4,A,Z4,A,Z4)") "SINGLE",
C     &                        NEXCITTYPES,
C     &                        SymLabels(2,Classes(I)%SymLab),"(",
C     &                        Classes(I)%SymRem,
C     &                        ")->",SymLabels(2,J)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
          ENDIF
         ENDIF
C.. And the doubles
         IF(BTEST(ILEVEL,1)) THEN
            DO I=1,NPR
               DO J=1,nSymPairProds
                  SPP=SymPairProds(J)%Sym
                  IF(LSYMSYM(SYMPROD(SymConj(SYMPRODS(I)), SPP))) THEN
                     DO ISPN=1,3
                        K=NALLOWPPS(ISPN,J)
                        L=SYMPRODIND(2,ISPN,I)
                        ICC=L*K
                        IF(ICC.NE.0) THEN
                           NEXCITTYPES=NEXCITTYPES+1
                           EXCITTYPES(1,NEXCITTYPES)=2
                           EXCITTYPES(2,NEXCITTYPES)=ISPN
                           EXCITTYPES(3,NEXCITTYPES)=I
                           EXCITTYPES(4,NEXCITTYPES)=J
                           EXCITTYPES(5,NEXCITTYPES)=ICC
C                        WRITE(6,"(2I5,2Z6,4I5)") I,J,SYMPRODS(I),
C     &                     SYMPAIRPRODS(J),ISPN,L,K,ICC
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
C.. Store all the pointers we need
         STORE(2)=IP_EXCITTYPES
         STORE(3)=IP_NALLOWPPS
         STORE(4)=IP_ORBPAIRS
         STORE(5)=IP_SYMPRODIND
         STORE(6)=NEXCITTYPES
         ENDIF
C..      ENDIF(.NOT.TCOUNT)
         call halt_timer(proc_timer)
      END

      SUBROUTINE SYMSETUPEXCITS3(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &   TCOUNT,ICOUNT,CLASSES,ILUT,SYMPRODS,ILEVEL,iMinElec1,iMaxElec1)
         use global_utilities
         use SystemData, only: Symmetry,BasisFN
         use SymData, only: SymClass
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL),NBASIS
         TYPE(BasisFN) G1(nBasis)
         INTEGER nBasisMax(5,*)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER*8 STORE(6)
         INTEGER ICOUNT
         LOGICAL TCOUNT
         INTEGER ILEVEL
         INTEGER iMinElec1,iMaxElec1

         TYPE(SymClass) CLASSES(*)
         TYPE(Symmetry) SYMPRODS(0:NEL*NEL)
         INTEGER CLASSCOUNT(2,NEL)
         INTEGER THISCLASSCOUNT(2,NEL)
!  ThisClassCount is used to list only electrons which this processor deals with
         INTEGER PREVCLASSCOUNT(2,NEL)
!  PrevClassCount is used to list electrons which lower indexed processors deal with
         INTEGER SYMPRODCOUNT(3,0:NEL*NEL)
         INTEGER SYMPRODIND(2,3,0:*)
         POINTER (IP_SYMPRODIND,SYMPRODIND)
         INTEGER ILUT(0:NBASIS/32)
         INTEGER EXCITTYPES(5,*)
         POINTER (IP_EXCITTYPES,EXCITTYPES)
         INTEGER nPairs, nSing, nDoub, nExcits
         INTEGER NALLOWPPS(3,*)
         POINTER(IP_NALLOWPPS,NALLOWPPS)
         INTEGER ORBPAIRS(2,*)
         POINTER (IP_ORBPAIRS,ORBPAIRS)
         INTEGER nCl,nExcitTypes,nSymPairProds,nPr

         LOGICAL ISUHFDET
 
         INTEGER I         
         type(timer), save :: proc_timer
         proc_timer%timer_name='SYMSUEXCIT'
         
         call set_timer(proc_timer,65)
         Call SymSetupExcits_CreateClassList(nI,nEl,Classes,
     & iMinElec1, iMaxElec1, ThisClassCount, PrevClassCount,ClassCount,
     & G1, nCl)
         CALL IAZZERO(SYMPRODCOUNT,(NEL*NEL+1)*3)
         Call SymSetupExcits_CreateClassSymProds(nPr,nPairs,nCl,
     &   SymProds, ThisClassCount, PrevClassCount, ClassCount,Classes,
     &   SymProdCount)
C.. Allocate enough memory to store the index
         IF(STORE(5).EQ.0) THEN
            CALL MEMORY(IP_SYMPRODIND,2*3*(NPR+1),'SYMPRODIND')   
         ELSE
C.. we've acually got space to store it somewhere
            IP_SYMPRODIND=STORE(5)
         ENDIF
         CALL IAZZERO(SYMPRODIND,2*3*(NPR+1))

C.. Now shift this such that SYMPRODCOUNT(ISPN,I) is the index of
C.. the first orbital pair of SYMPROD(I) with spin ISPN in ORBPAIRS
C.. Store in SYMPRODIND(1,ISPN,I) too
         DO I=NPR,1,-1
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)
            SYMPRODIND(1,3,I)=SYMPRODCOUNT(3,I)
            SYMPRODIND(1,2,I)=SYMPRODCOUNT(2,I)
            SYMPRODIND(1,1,I)=SYMPRODCOUNT(1,I)
         ENDDO
C         DO I=1,NPR
C            WRITE(6,"(Z4,3I4)") SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
C         ENDDO
C.. We allocate enough memory to store all the pairs.
C.. Each pair consists of (ORB1,ORB2) where ORB1<ORB2
         IF(STORE(4).EQ.0) THEN
            CALL MEMORY(IP_ORBPAIRS,2*NPAIRS/IRAT+1,'ORBPAIRS')
         ELSE
            IP_ORBPAIRS=STORE(4)
         ENDIF

         Call SymSetupExcits_StoreOccPairs(OrbPairs, nPairs, nPr,
     &    iMinElec1, iMaxElec1,G1,SymProdInd,SymProds,nI,nEl)
C.. Now go through the list of all pairs, finding out how many pairs are to be excluded as they contain some of the
C.. orbitals in this determinant.
C.. We first create a quick lookup table to enable us to quickly check whether a given orbital is in this determinant
C.. in order 1, not order NEL.
         CALL IAZZERO(ILUT,NBASIS/32+1)
         DO I=1,NEL
            ILUT((NI(I)-1)/32)=IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32))
         ENDDO
!         DO I=0,NBASIS/32
!            WRITE(6,"(A,Z10)") "LUT: ",ILUT(I)
!         ENDDO
C.. Now look through the list of our pairs.  For each pair sym of the complete list which has a 
C.. symmetric product with any of our pair syms, we work out how many allowed pairs there are in
C.. the complete list, and store that value in NALLOWPPS
         IF(STORE(3).EQ.0) THEN
            CALL MEMORY(IP_NALLOWPPS,3*NSYMPAIRPRODS/IRAT+1,'NALLOWPPS')
         ELSE
            IP_NALLOWPPS=STORE(3)
         ENDIF
         Call SymSetupExcits_CountVirtProds(nDoub, nExcitTypes,nPr,
     &   nSymPairProds, nAllowPPS,iLUT)
         NEXCITS=0
         IF(BTEST(ILEVEL,0)) NEXCITS=NEXCITS+NSING
         IF(BTEST(ILEVEL,1)) NEXCITS=NEXCITS+NDOUB
         ICOUNT=NEXCITS
C         WRITE(6,*) "Total number of excitations: ",NEXCITS
C         WRITE(6,*) "Total number of excitation types: ",NEXCITTYPES

         IF(TCOUNT) THEN
C.. If we're just counting, we're done, so we get rid of some pointers.
            IP_EXCITTYPES=0
C.. However, we do save the length of the memory required.
C.. EXCITTYPES
            STORE(2)=(NEXCITTYPES*5)
C.. NALLOWPPS
            STORE(3)=3*NSYMPAIRPRODS
C.. ORBPAIRS
            STORE(4)=2*NPAIRS
C.. SYMPRODIND
            STORE(5)=2*3*(NPR+1)
C.. indicate that these are lengths
            STORE(6)=0
            
            CALL FREEM(IP_NALLOWPPS)
            CALL FREEM(IP_ORBPAIRS)
            CALL FREEM(IP_SYMPRODIND)
         ELSE
C.. Now allocate memory to store all the excitation types if there hasn't been one already allocated.
C.. This will have to be manually deallocated later.
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
         IF(STORE(2).EQ.0) THEN
          CALL MEMORY(IP_EXCITTYPES,(NEXCITTYPES*5)/IRAT+1,'EXCITTYPES')
         ELSE
            IP_EXCITTYPES=STORE(2)
         ENDIF

         IF(.NOT.ISUHFDET(NI,NEL)) THEN
          Call SymSetupExcits_CountSingles(nSing,nCl,nExcitTypes,
     &   ThisClassCount, ClassCount,Classes)
         ELSE
          nSing=0
         ENDIF
         nExcitTypes=0
         IF(.NOT.ISUHFDET(NI,NEL)) THEN
            IF(BTEST(ILEVEL,0)) THEN
               Call SymSetupExcits_StoreSingles(nExcitTypes,nCl,
     &   Classes,ThisClassCount,ExcitTypes)
            ENDIF
         ENDIF
         IF(BTEST(ILEVEL,1)) THEN
            Call SymSetupExcits_StoreDoubles(nPr,nSymPairProds,
     &   nAllowPPS,ExcitTypes,nExcitTypes,SymProds,SymProdInd)
         ENDIF
C.. Store all the pointers we need
         STORE(2)=IP_EXCITTYPES
         STORE(3)=IP_NALLOWPPS
         STORE(4)=IP_ORBPAIRS
         STORE(5)=IP_SYMPRODIND
         STORE(6)=NEXCITTYPES
         ENDIF
C..      ENDIF(.NOT.TCOUNT)
         call halt_timer(proc_timer)

      END Subroutine SymSetupExcits3

C.. List the symmetry classes in NI.
C.. For each orb, just see if its class is in the list, and if not add it
C.. Each state's symmetry falls into a class SymClasses(ISTATE).

C.. We divide the occupied electrons into classes.  A class correpsonds to a unique value of SymClasses(iState) for state (i.e. spatial orbital) iState.
C.. The CLASSCOUNT arrays contain the number of electrons of each class.
C.. We keep three lists:   CLASSCOUNT  which takes into account all electrons.  
C..                     THISCLASSCOUNT which takes into account electrons allocated to this processor
C..  and                PREVCLASSCOUNT which takes into account electrons allocated to the previous processor.
C.. These both have the format (Spn,iCl), where Spin=1,2 corresponding to alpha and beta and
C..                                             iCl is the index of the Class.
C.. CLASSES(iCl)%SymLab is the symmetry label (i.e. entry in SymClasses) of this class (NB not a symmetry itself)

      Subroutine SymSetupExcits_CreateClassList(nI,nEl,Classes,
     & iMinElec1, iMaxElec1, ThisClassCount, PrevClassCount,ClassCount,
     &  G1, nCl)
         use SystemData, only: BasisFN
         use SystemData, only: Symmetry
         use SymData, only: SymClass,SymClasses,SymLabels
         IMPLICIT NONE
         INTEGER nEl, nI(nEl)
         Type(BasisFN) G1(*)
         TYPE(SymClass) CLASSES(*)     ! Data about each class
         INTEGER nCl                   ! Will be the number of classes
         INTEGER CLASSCOUNT(2,nEl)     ! Number of alpha and beta electrons in each class
         INTEGER THISCLASSCOUNT(2,nEl) 
!  ThisClassCount is used to list only electrons which this processor deals with
         INTEGER PREVCLASSCOUNT(2,nEl)
!  PrevClassCount is used to list electrons which lower indexed processors deal with
         INTEGER iMinElec1,iMaxElec1

         LOGICAL tNew
         INTEGER I,J,K,L
         INTEGER iSpn
         Type(Symmetry) Pr
         
         Type(Symmetry) TotSymRep
         TYPE(Symmetry) SYMPROD,SymConj
         
         NCL=0
         DO I=1,NEL
            TNEW=.TRUE.
!            WRITE(6,*) I,SymClasses((NI(I)+1)/2)
            DO J=1,NCL
               IF(SymClasses((NI(I)+1)/2).EQ.CLASSES(J)%SymLab) THEN
                  ISPN=(G1(NI(I))%Ms+3)/2
                  CLASSCOUNT(ISPN,J)=CLASSCOUNT(ISPN,J)+1
                  IF(I.GE.iMinElec1.AND.i.LE.iMaxElec1) THEN
                     THISCLASSCOUNT(ISPN,J)=THISCLASSCOUNT(ISPN,J)+1
                  ENDIF
                  IF(I.LT.iMinElec1) THEN
                     PREVCLASSCOUNT(ISPN,J)=PREVCLASSCOUNT(ISPN,J)+1
                  ENDIF
                  TNEW=.FALSE.
                  EXIT
               ENDIF
            ENDDO
            IF(TNEW) THEN
C.. add the new class
C.. ISPN=1 is beta, ISPN=2 is alpha
               ISPN=(G1(NI(I))%Ms+3)/2
               NCL=NCL+1
               CLASSES(NCL)%SymLab=SymClasses((NI(I)+1)/2)
               CLASSCOUNT(ISPN,NCL)=1
               CLASSCOUNT(3-ISPN,NCL)=0
               IF(I.GE.iMinElec1.AND.i.LE.iMaxElec1) THEN
                  THISCLASSCOUNT(ISPN,NCL)=1
                  THISCLASSCOUNT(3-ISPN,NCL)=0
               ELSE
                  THISCLASSCOUNT(1:2,NCL)=0
               ENDIF
               IF(I.LT.iMinElec1) THEN
                  PREVCLASSCOUNT(ISPN,NCL)=1
                  PREVCLASSCOUNT(3-ISPN,NCL)=0
               ELSE
                  PREVCLASSCOUNT(1:2,NCL)=0
               ENDIF
            ENDIF
         ENDDO
C.. Now calculate the resultant symmetry from removing one of each of the classes (for singles)
C.. i.e. if we remove j, what's the symmetry of SYMCONJ(<ikl... |) (X) SYM( |ikl...>)
C..   NB to cope with non-abelians, this actually takes the prod sym of the classes of each orb, 
C..   rather than the orb itself, so may end up classifying a single as allowed when it's not.
C.. we'll use this in the future to determina the allowed j->a
C.. This is stored in CLASSES(iCl)%SymRem
         DO I=1,NCL
            PR=TotSymRep()
            DO J=1,NCL
               L=CLASSCOUNT(1,I)+CLASSCOUNT(2,I)
               IF(J.EQ.I) L=L-1
               DO K=1,L
                  PR=SYMPROD( PR,SymConj(SymLabels(CLASSES(J)%SymLab)))
                  PR=SYMPROD( PR,SymLabels(CLASSES(J)%SymLab))
               ENDDO
            ENDDO
            CLASSES(I)%SymRem=PR
!            WRITE(6,"(A,I4,2Z4,I)") "CLASS ",NCL,CLASSES(I)%SymLab,
!     &            CLASSES(I)%SymRem,I
         ENDDO
!            WRITE(6,*) NCL," Symmetry Classes"
      End Subroutine

C.. Now determine all possible symmetry products of pairs of orbs, and classify them.
C.. SYMPRODS(iPr) is the a symmetry product.
C.. SymProdCount(iType,iPr) is the number of prods of a given type.  iType=0 (beta/beta), 1 (alpha/beta), 2 (alpha/alpha)
C.. NPairs totals all SymProdCounts.

C.. For each pair of symmetry classes, determine its symmetry product.
C.. We list them and then sort them
      Subroutine SymSetupExcits_CreateClassSymProds(nPr,nPairs,nCl,
     &   SymProds, ThisClassCount, PrevClassCount, ClassCount,Classes,
     &   SymProdCount)
         use SystemData, only: Symmetry,SymmetrySize
         use SymData, only: SymClass,SymClasses,SymLabels
         IMPLICIT NONE
         INTEGER nPr,nPairs,nCl
         TYPE(Symmetry) SymProds(0:*)
         INTEGER CLASSCOUNT(2,nCl)     ! Number of alpha and beta electrons in each class
         INTEGER THISCLASSCOUNT(2,nCl) 
!  ThisClassCount is used to list only electrons which this processor deals with
         INTEGER PREVCLASSCOUNT(2,nCl)
!  PrevClassCount is used to list electrons which lower indexed processors deal with
         INTEGER SYMPRODCOUNT(3,0:*)
         TYPE(SymClass) CLASSES(*)
        
         INTEGER I,J,K
         LOGICAL tNew
         INTEGER ICC1,ICC2,ICC3
          
         Type(Symmetry) TotSymRep
         TYPE(Symmetry) SYMPROD,SymConj
         LOGICAL SYMEQ

         Type(Symmetry) Pr
         NPR=0
         NPAIRS=0
         DO I=1,NCL
            DO J=I,NCL
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SymLabels(ISYMLABEL)
               PR=SYMPROD( SymLabels(CLASSES(I)%SymLab),
     &                     SymLabels(CLASSES(J)%SymLab))
               TNEW=.TRUE.
               DO K=1,NPR
                  IF(SYMEQ(SYMPRODS(K),PR)) THEN
                     TNEW=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF(TNEW) THEN
                  NPR=NPR+1
                  K=NPR
                  SYMPRODS(NPR)=PR
               ENDIF
C.. 1 corresponds to beta beta, 2 to alpha/beta and beta/alpha, and 3 to alpha, alpha
C.. 1 corresponds to beta beta, 2 to alpha/beta and beta/alpha, and 3 to alpha, alpha
!  The number of pairs is somewhat more complex when we split electron 1 between processors
!  For systems where we require the upper-triangle (excluding diagonal) we have
!  NThis*(NAll+(NAll-NThis)-1)/2  pairs
!  For cases where we allow all pairs, it's just NThis*NAll
               IF(I.EQ.J) THEN
                  ICC1=THISCLASSCOUNT(1,I)*
     &               (2*CLASSCOUNT(1,I)-THISCLASSCOUNT(1,I)
     &                  -2*PREVCLASSCOUNT(1,I)-1)/2
!If we have all electrons being considered here then the alpha beta lot is simple, and just Nalpha*Nbeta
!                  ICC2=CLASSCOUNT(1,I)*CLASSCOUNT(2,I)
!  However, if we are splitting electron 1 between processors, alpha electrons here can only interact with beta
!    electrons with higher index, and the same witb beta on alpha
!  Therefore we have NThisA*(NAllB-NPrevB)+ NThisB*(NAllA-NPrevA-NThisA)
                  ICC2=THISCLASSCOUNT(1,I)*
     &               (CLASSCOUNT(2,I)-PREVCLASSCOUNT(2,I))
     &               + THISCLASSCOUNT(2,I)*
     &         (CLASSCOUNT(1,I)-PREVCLASSCOUNT(1,I)-THISCLASSCOUNT(1,I))

                  ICC3=THISCLASSCOUNT(2,I)*
     &               (2*CLASSCOUNT(2,I)-THISCLASSCOUNT(2,I)
     &               -2*PREVCLASSCOUNT(2,I)-1)/2
               ELSE
                  ICC1=THISCLASSCOUNT(1,I)*CLASSCOUNT(1,J)
                  ICC2=THISCLASSCOUNT(1,I)*CLASSCOUNT(2,J)
     &                +THISCLASSCOUNT(2,I)*CLASSCOUNT(1,J)
                  ICC3=THISCLASSCOUNT(2,I)*CLASSCOUNT(2,J)
               ENDIF
               SYMPRODCOUNT(1,K)=SYMPRODCOUNT(1,K)+ICC1
               SYMPRODCOUNT(2,K)=SYMPRODCOUNT(2,K)+ICC2
               SYMPRODCOUNT(3,K)=SYMPRODCOUNT(3,K)+ICC3
     
               NPAIRS=NPAIRS+ICC1+ICC2+ICC3
!                  WRITE(6,"(6I6)") CLASSCOUNT(1:2,I),
!     &  THISCLASSCOUNT(1:2,I), PREVCLASSCOUNT(1:2,I)
!                  WRITE(6,"(6I6)") CLASSCOUNT(1:2,J),
!     &  THISCLASSCOUNT(1:2,J), PREVCLASSCOUNT(1:2,J)
!                  WRITE(6,"(A,I6,Z4,5I6)") "SYMPROD ",K,SYMPRODS(K),I,J,
!     &               ICC1,ICC2,ICC3
            ENDDO
         ENDDO
         SYMPRODS(0)%s=0 ! For historical reasons.  Should never be needed.
         SYMPRODCOUNT(1,0)=0
         SYMPRODCOUNT(2,0)=0
         SYMPRODCOUNT(3,0)=0
C         DO I=1,NPR
C            WRITE(6,"(Z4,3I4)") SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
C         ENDDO

C.  We sort the array of Symmetry products (and its associated SymProdCount info), into ascending order of product, for easy searching later.
         CALL SORTSYMIARR(NPR,SymProds(1),SymmetrySize,
     &                                             SymProdCount(:,1),3)
!         WRITE(6,*) NPR," Symmetry Products"
!         WRITE(6,*) NPAIRS," Orbital Pairs"
C.. Create a cumulative sum.
         DO I=1,NPR
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)+SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)+SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)+SYMPRODCOUNT(3,I)
C            WRITE(6,*) SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
         ENDDO
      End Subroutine

C.. Classify and store each possible pair of orbitals under its symmetry product.
C.. With each symmetry product (denoted []), calculate []'' which is the 
C.. remainder after having removed the orbitals.


C.. Now classify each pair of orbs under its sym prod
C.. SYMPRODS(I) specifies a sym prod.  Orbitals for this are stored in
C.. ORBPAIRS starting at index SYMPRODIND(1,ISPN,I)+1.
C.. The first free space is at SYMPRODIND(1,ISPN,I)+SYMPRODIND(2,ISPN,I)+1.
      Subroutine SymSetupExcits_StoreOccPairs(OrbPairs, nPairs, nPr,
     &    iMinElec1, iMaxElec1,G1,SymProdInd,SymProds,nI,nEl)
         use SystemData, only: Symmetry,BasisFN,SymmetrySize
         IMPLICIT NONE
         INTEGER ORBPAIRS(2,nPairs)
         INTEGER nPairs
         INTEGER nPr
         INTEGER iMinElec1, iMaxElec1,nEl,nI(nEl)
         TYPE(Symmetry) SymProds(0:*)

         INTEGER I,J,K,L
         Type(Symmetry) Pr
         INTEGER iSpn
         INTEGER iCC
         TYPE(BasisFN) G1(*)
         INTEGER SYMPRODIND(2,3,0:*)

         Type(Symmetry) SymProd

         CALL IAZZERO(ORBPAIRS,2*NPAIRS)
         L=0
         DO I=iMinElec1,iMaxElec1
            
            DO J=I+1,NEL
               L=L+1
               PR=SYMPROD(G1(NI(I))%Sym,G1(NI(J))%Sym)
C.. Find the product in the sorted symprod list
               CALL BINARYSEARCHSYM(PR,SYMPRODS(1),SymmetrySize,NPR,K)
               IF(K.EQ.0) THEN
                  WRITE(6,*) "Occupied Symmetry Products"
                  DO L=0,NPR
                     WRITE(6,"(I4)",advance='no') L
                     CALL WRITESYM(6,SYMPRODS(L),.TRUE.)
                  ENDDO
                  WRITE(6,"(A)",advance='no') "Illegal Symmetry"
                  CALL WRITESYM(6,PR,.TRUE.)
                  STOP 'Illegal Symmetry Found'
               ENDIF
               ISPN=(G1(NI(I))%Ms+G1(NI(J))%Ms)/2+2

C.. Increment the offset for the next one
               SYMPRODIND(2,ISPN,K)=SYMPRODIND(2,ISPN,K)+1
               ICC=SYMPRODIND(1,ISPN,K)+SYMPRODIND(2,ISPN,K)
C..                 The index          +  offset of the element

               ORBPAIRS(1,ICC)=NI(I)
               ORBPAIRS(2,ICC)=NI(J)
               
!               WRITE(6,"(A,2I4,2Z4,2I4)") 
!     &            "SP",ICC,ISPN,PR,SYMPRODS(K),NI(I),NI(J)
C               WRITE(6,*) L,NI(I),NI(J),PR,K
            ENDDO
         ENDDO
      End Subroutine

C.. With each of the symmetry products, calculate []' such that []x[]'
C.. contains the totally symmetric rep.
C.. We do this by checking whether any of the []' in the 
C.. global symprods table multiplied by [] give the sym rep.

      Subroutine SymSetupExcits_CountVirtProds(nDoub, nExcitTypes,nPr,
     &   nSymPairProds, SymProdInd, SymProds, nAllowPPS,iLUT)
         use SystemData, only: Symmetry,BasisFN
         use SymData, only: SymPairProds,SymStatePairs
         IMPLICIT NONE
         INTEGER nDoub
         INTEGER nExcitTypes
         INTEGER nPr
         INTEGER nSymPairProds
         INTEGER nAllowPPS(3,*)
         INTEGER iLUT(*)
         INTEGER SYMPRODIND(2,3,0:*)
         TYPE(Symmetry) SYMPRODS(0:*)

         INTEGER I,J,K,L
         Type(Symmetry) SPP
         INTEGER ICC1,ICC2,ICC3,ICC4
         LOGICAL L1B,L1A,L2B,L2A
         INTEGER ICC,iSpn
   
         Type(Symmetry) SymConj,SymProd
         LOGICAL LSYMSYM

C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements


C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) is a global list of all SYMPRODs available, the number of pairs of
C.. symlabels (listed in SymStatePairs), and the index of the start of this list
C.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
C.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of symlabels whose prod is of that symmetry.

         NDOUB=0
         NEXCITTYPES=0
         DO I=1,NPR
            DO J=1,NSYMPAIRPRODS
C.. If it's the first time around, we flag NALLOWPPS as invalid
               IF(I.EQ.1) THEN
                  NALLOWPPS(1,J)=-1
                  NALLOWPPS(2,J)=-1
                  NALLOWPPS(3,J)=-1
               ENDIF
C.. [] is the SymProd(I).  For each sympairprod, check if [] x []' contains the sym rep
C               WRITE(6,*) SYMPRODS(I),SYMPAIRPRODS(J)
               SPP=SymPairProds(J)%Sym
               IF(LSYMSYM(SYMPROD(SymConj(SYMPRODS(I)), SPP))) THEN
C.. Check if we've worked out the allowed number of excitation pairs, otherwise work it out
                  IF(NALLOWPPS(1,J).EQ.-1) THEN
C.. Zero the counters
                     NALLOWPPS(1,J)=0
                     NALLOWPPS(2,J)=0
                     NALLOWPPS(3,J)=0

                     DO K=SymPairProds(J)%nIndex,
     &                 SymPairProds(J)%nIndex+SymPairProds(J)%nPairs-1

C.. Now check according to ISPN
C.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
C.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
                        ICC1=SymStatePairs(1,K)*2-1
                        ICC2=ICC1+1
                        ICC3=SymStatePairs(2,K)*2-1
                        ICC4=ICC3+1
                        L1B=.NOT.BTEST(ILUT((ICC1-1)/32),
     &                        MOD(ICC1-1,32))
                        L1A=.NOT.BTEST(ILUT((ICC2-1)/32),
     &                        MOD(ICC2-1,32))
                        L2B=.NOT.BTEST(ILUT((ICC3-1)/32),
     &                        MOD(ICC3-1,32))
                        L2A=.NOT.BTEST(ILUT((ICC4-1)/32),
     &                        MOD(ICC4-1,32))
C.. L1B is set if the beta of the first virtual is not in NI, i.e. is allowed
!  NALLOWPPS is the number of ALLOWED virtuals
                        IF(ICC1.EQ.ICC3) THEN
!                           NALLOWPPS(1,J)=NALLOWPPS(1,J)+1
!                           NALLOWPPS(3,J)=NALLOWPPS(3,J)+1
C.. The virtuals we are exciting to are the same state.  We're only allowed an AB excitation
                         IF(L1B.AND.L2A) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                        ELSE
                         IF(L1B.AND.L2B) NALLOWPPS(1,J)=NALLOWPPS(1,J)+1
                         IF(L1B.AND.L2A) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                         IF(L1A.AND.L2B) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                         IF(L1A.AND.L2A) NALLOWPPS(3,J)=NALLOWPPS(3,J)+1
                        ENDIF
C                        WRITE(6,*) SPP,SymStatePairs(1,K)*2-1,
C     &                     SymStatePairs(2,K)*2-1
C     &                     ,L1B,L1A,L2B,L2A
                     ENDDO
                  ENDIF
!.. Now count all the excitations allowed for this set of []->[]'
                  DO ISPN=1,3
                     K=NALLOWPPS(ISPN,J)
                     L=SYMPRODIND(2,ISPN,I)
                     ICC=L*K
C                     WRITE(6,"(2I5,2Z6,4I5)") I,J,SYMPRODS(I),
C     &                  SYMPAIRPRODS(J),ISPN,L,K,ICC
                     NDOUB=NDOUB+ICC
                     IF(ICC.GT.0) NEXCITTYPES=NEXCITTYPES+1
                  ENDDO
               ENDIF
            ENDDO
         ENDDO       

!         WRITE(6,*) "Number of double excitations: ",NDOUB
      End Subroutine

      Subroutine SymSetupExcits_CountSingles(nSing,nCl,nExcitTypes,
     &   ThisClassCount, ClassCount,Classes)
         use SystemData, only: Symmetry
         use SymData, only: SymClass,SymLabelCounts,nSymLabels,SymLabels
         IMPLICIT NONE
         INTEGER nSing
         INTEGER nCl
         INTEGER nExcitTypes
         INTEGER ThisClassCount(2,*)
         INTEGER ClassCount(2,*)
         Type(SymClass) Classes(nCl)
         
         Type(Symmetry) SymProd,SymConj
         LOGICAL LSymSym

         INTEGER I,J,K
         INTEGER iSpn
         INTEGER iCC
         NSING=0
         DO I=1,NCL
C.. For each class, see what sym prods it can interact with
            DO J=1,NSYMLABELS
               IF(LSYMSYM(SYMPROD(CLASSES(I)%SymRem,
     &                  SYMPROD(SymConj(SymLabels(J)),
     &                         SymLabels(CLASSES(I)%SymLab))))) THEN
                  DO ISPN=1,2
                     ICC=0
! Work out which virtuals are available
                     DO K=1,NCL
                       IF(CLASSES(K)%SymLab.EQ.J)
     &                  ICC=CLASSCOUNT(ISPN,K)
                     ENDDO
                     ICC=THISCLASSCOUNT(ISPN,I)*
     &                  (SYMLABELCOUNTS(2,J)-ICC)
                     NSING=NSING+ICC
                     IF(ICC.NE.0) NEXCITTYPES=NEXCITTYPES+1
C                     WRITE(6,"(2Z4,4I)") SymLabels(CLASSES(I)),
C     &                  SymLabels(J),ISPN,ICC,
C     &                  CLASSCOUNT(ISPN,I),
C     &                  (SYMLABELCOUNTS(2,J)-CLASSCOUNT(ISPN,I))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO         

!         WRITE(6,*) "Number of single excitations: ",NSING
      End Subroutine
      
      Subroutine SymSetupExcits_StoreSingles(nExcitTypes,nCl,
     &   Classes,ThisClassCount,ExcitTypes)
         Use SystemData, only: Symmetry
         use SymData, only: SymClass,SymPairProds,nSymLabels,SymLabels
         use SymData, only: SymLabelCounts
         IMPLICIT NONE
         INTEGER nExcitTypes
         INTEGER nCl
         Type(SymClass) Classes(*)
         INTEGER ThisClassCount(2,*)
         INTEGER ExcitTypes(5,*)

         INTEGER I,J,K
         Type(Symmetry) Pr,Pr2
         INTEGER iSpn
         INTEGER iCC

         Type(Symmetry) SymConj,SymProd
         LOGICAL LSYMSYM

C.. Now store the singles:
         NEXCITTYPES=0
C         CALL WRITEDET(6,NI,NEL,.TRUE.)
         DO I=1,NCL
C.. For each class, see what sym prods it can interact with
            DO J=1,NSYMLABELS

               PR=SYMPROD(SymConj(SymLabels(J)),
     &                       SymLabels(CLASSES(I)%SymLab))
               PR2=SYMPROD(PR,Classes(I)%SymRem)
C                  WRITE(6,"(2I4,5Z4)") I,J,
C     &                        CLASSES(I)%SymRem,
C     &                         SymLabels(J),
C     &                         SymLabels(CLASSES(I)%SymLab),
C     &                         PR,PR2
               IF(LSYMSYM(PR2)) THEN
                  DO ISPN=1,2
                     ICC=0
                     DO K=1,NCL
                        IF(CLASSES(K)%SymLab.EQ.J) 
     &                        ICC=THISCLASSCOUNT(ISPN,K)
                     ENDDO
                     ICC=THISCLASSCOUNT(ISPN,I)*
     &                     (SYMLABELCOUNTS(2,J)-ICC)
                     IF(ICC.NE.0) THEN
                        NEXCITTYPES=NEXCITTYPES+1
                        EXCITTYPES(1,NEXCITTYPES)=1
                        EXCITTYPES(2,NEXCITTYPES)=ISPN
                        EXCITTYPES(3,NEXCITTYPES)=I
                        EXCITTYPES(4,NEXCITTYPES)=J
                        EXCITTYPES(5,NEXCITTYPES)=ICC
C                           WRITE(6,"(A,I,Z4,A,Z4,A,Z4)") "SINGLE",
C     &                        NEXCITTYPES,
C     &                        SymLabels(2,Classes(I)%SymLab),"(",
C     &                        Classes(I)%SymRem,
C     &                        ")->",SymLabels(2,J)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      End Subroutine

      Subroutine SymSetupExcits_StoreDoubles(nPr,nSymPairProds,
     &   nAllowPPS,ExcitTypes,nExcitTypes,SymProds,SymProdInd)
         Use SystemData, only: Symmetry
         Use SymData, only: SymPairProds
         IMPLICIT NONE
         INTEGER nPr
         INTEGER nSymPairProds
         INTEGER nAllowPPS(3,*)
         INTEGER ExcitTypes(5,*)
         INTEGER nExcitTypes
         TYPE(Symmetry) SYMPRODS(0:*)
         INTEGER SYMPRODIND(2,3,0:*)
         
         Type(Symmetry) SymConj,SymProd
         LOGICAL LSYMSYM

         INTEGER I,J,K,L
         INTEGER ICC
         INTEGER iSpn
         Type(Symmetry) SPP
         
C.. And the doubles
         DO I=1,NPR
            DO J=1,nSymPairProds
               SPP=SymPairProds(J)%Sym
               IF(LSYMSYM(SYMPROD(SymConj(SYMPRODS(I)), SPP))) THEN
                  DO ISPN=1,3
                     K=NALLOWPPS(ISPN,J)
                     L=SYMPRODIND(2,ISPN,I)
                     ICC=L*K
                     IF(ICC.NE.0) THEN
                        NEXCITTYPES=NEXCITTYPES+1
                        EXCITTYPES(1,NEXCITTYPES)=2
                        EXCITTYPES(2,NEXCITTYPES)=ISPN
                        EXCITTYPES(3,NEXCITTYPES)=I
                        EXCITTYPES(4,NEXCITTYPES)=J
                        EXCITTYPES(5,NEXCITTYPES)=ICC
C                        WRITE(6,"(2I5,2Z6,4I5)") I,J,SYMPRODS(I),
C     &                     SYMPAIRPRODS(J),ISPN,L,K,ICC
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      End Subroutine


!  Must be called after SYMSETUPEXCITS.  Needs INTEGER*8 STORE and frees memory therein
      SUBROUTINE SYMGENEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,LSTE,
     &   ICLIST,ICOUNT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         use global_utilities
         IMPLICIT NONE
         INTEGER NBASIS,NEL,NI(NEL),G1(nBasis),NBASISMAX(*)
         INTEGER*8 STORE(6)
         INTEGER ORBPAIRS(2,*)
         POINTER (IP_ORBPAIRS,ORBPAIRS)
         INTEGER SYMPRODIND(2,3,0:*)
         POINTER (IP_SYMPRODIND,SYMPRODIND)
         INTEGER NALLOWPPS(3,*)
         TYPE(Symmetry) SPP
         POINTER(IP_NALLOWPPS,NALLOWPPS)

         INTEGER EXCITTYPES(5,*)
         POINTER (IP_EXCITTYPES,EXCITTYPES)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         type(timer), save :: proc_timer
         INTEGER I
         INTEGER LSTE(NEL,*),ICLIST(*),ICOUNT
         proc_timer%timer_name='SYMGNEXCIT'
         call set_timer(proc_timer,65)

         IP_DSTORE=STORE(1)
         IP_EXCITTYPES=STORE(2)
         IP_NALLOWPPS=STORE(3)
         IP_ORBPAIRS=STORE(4)
         IP_SYMPRODIND=STORE(5)
         
         CALL SYMGENALLEXCITS(NI,NEL,EXCITTYPES,STORE(6),DSTORE(1),
     &            SYMPRODIND,DSTORE(SymClassSize*NEL+1),ORBPAIRS,G1,
     &            NBASISMAX,LSTE,ICLIST,ICOUNT)
   
         CALL FREEM(IP_EXCITTYPES)
         CALL FREEM(IP_NALLOWPPS)
         CALL FREEM(IP_ORBPAIRS)
         CALL FREEM(IP_SYMPRODIND)
         CALL FREEM(IP_DSTORE)

         call halt_timer(proc_timer)
 
      END
      SUBROUTINE SYMGENALLEXCITS(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
     &               SYMPRODIND,ILUT,ORBPAIRS,G1,NBASISMAX,
     &               LSTE,ICLIST,ICOUNT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES),NEXCITTYPES
         
         INTEGER I,J,K,L
         INTEGER NJ(NEL),NK(NEL)
         TYPE(SymClass) CLASSES(*)
         INTEGER IEXCIT
         INTEGER IFROM,ITO,ISPN
         INTEGER IFROMSL,ITOSL
         TYPE(Symmetry) SPP
         LOGICAL L1,L2
         INTEGER ICC1,ICC2,ICC3,ICC4
         LOGICAL L1B,L1A,L2B,L2A

         INTEGER ORBPAIRS(2,*)
         INTEGER ILUT(0:*)
         INTEGER SYMPRODIND(2,3,0:*)
         INTEGER G1(*),NBASISMAX(*)
         INTEGER LSTE(NEL,*),ICLIST(*),ICOUNT
         type(BasisFN) ISYM

         ICOUNT=0

C.. Now work out all the excitations
         DO IEXCIT=1,NEXCITTYPES
            IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
C.. a single excitation
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
               ISPN=EXCITTYPES(2,IEXCIT)-2
               IFROM=CLASSES(EXCITTYPES(3,IEXCIT))%SymLab
               ITO=EXCITTYPES(4,IEXCIT)
               J=0
               I=1
C               WRITE(6,*) ISPN,IFROM,ITO,EXCITTYPES(5,IEXCIT)
C.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
C.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I
               DO WHILE (I.LE.NEL)
                  L1=J.LT.SYMLABELCOUNTS(2,IFROM)
                  DO WHILE(L1)
                     IFROMSL=(SYMLABELLIST(SYMLABELCOUNTS(1,IFROM)+J)*2
     &                  +ISPN)
C                     WRITE(6,*) I,NI(I),J,IFROM,IFROMSL
                     IF(IFROMSL.LT.NI(I)) J=J+1
                     L1=IFROMSL.LT.NI(I)
     &                  .AND.J.LT.SYMLABELCOUNTS(2,IFROM)
                  ENDDO
C                  WRITE(6,*) I,J,SYMLABELCOUNTS(2,IFROM)
                  IF(J.LT.SYMLABELCOUNTS(2,IFROM).AND.IFROMSL.EQ.NI(I))
     &               THEN
C.. We've found an orb in NI with the correct sym.  Now go through the list of possible excitations of it
                     K=0
                     L=0
                     DO WHILE(K.LT.SYMLABELCOUNTS(2,ITO))
                        L2=.TRUE.
                        ITOSL=2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)
     &                              +ISPN
                        DO WHILE (L2)
                           L=L+1
                           L2=L.LE.NEL.AND.NI(L).LT.ITOSL
                        ENDDO
                        IF(L.GT.NEL.OR.NI(L).NE.ITOSL) THEN
C.. We've found a virtual into which we can excite our occupied orb.
                           CALL ICOPY(NEL,NI,1,NK,1)
                           NK(I)=ITOSL
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=1
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM 
                        ENDIF
                        L=L-1
                        K=K+1
                     ENDDO
                  ENDIF
                  I=I+1
               ENDDO
            ELSE
C.. a double excitation
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
               ISPN=EXCITTYPES(2,IEXCIT)
               IFROM=EXCITTYPES(3,IEXCIT)
               ITO=EXCITTYPES(4,IEXCIT)
C               WRITE(6,*) "EXC",IEXCIT,ISPN,IFROM,ITO,SYMPRODIND(2,ISPN,
C     &            IFROM)
C.. Go through the list of pairs with a given symprod.
C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements
               DO I=1,SYMPRODIND(2,ISPN,IFROM)
C.. Now go through the list of virtual pairs, excluding those with orbitals in NI 
                  SPP=SymPairProds(ITO)%Sym
C                  WRITE(6,*) "SPP:", SPP
C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
C.. states (listed in SymStatePairs), and the index of the start of this list
C.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
C.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose prod is of that symmetry.
                  DO K=SymPairProds(ITO)%nIndex,
     &              SymPairProds(ITO)%nIndex+SymPairProds(ITO)%nPairs-1
C.. Now check according to ISPN
C.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
C.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
                     ICC1=SymStatePairs(1,K)*2-1
                     ICC2=ICC1+1
                     ICC3=SymStatePairs(2,K)*2-1
                     ICC4=ICC3+1
C                     WRITE(6,*) ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
C     &                      ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)
C                     WRITE(6,*) ISPN,ICC1,ICC3
                     L1B=BTEST(ILUT((ICC1-1)/32),
     &                     MOD(ICC1-1,32))
                     L1A=BTEST(ILUT((ICC2-1)/32),
     &                     MOD(ICC2-1,32))
                     L2B=BTEST(ILUT((ICC3-1)/32),
     &                     MOD(ICC3-1,32))
                     L2A=BTEST(ILUT((ICC4-1)/32),
     &                     MOD(ICC4-1,32))
C                     WRITE(6,*) L1B,L1A,L2B,L2A
C.. L1B is set if the beta of the first virtual is in NI, i.e. is disallowed
                     IF(ISPN.EQ.1) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                        IF(ICC1.NE.ICC3.AND..NOT.(L1B.OR.L2B)) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC1
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC3
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ELSEIF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
                        IF(.NOT.(L1B.OR.L2A)) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC1
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC4
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
                        IF(.NOT.(L1A.OR.L2B).AND.ICC1.NE.ICC3) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC2
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC3
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ELSEIF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                        IF(ICC1.NE.ICC3.AND..NOT.(L1A.OR.L2A)) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC2
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC4
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO 
C         WRITE(6,*) "COUNT:",ICOUNT
      END

! For a symmetry excitation generator stored in NMEM, return the total number of excitations from it in iCount.
      SUBROUTINE GetSymExcitCount(NMEM,iCount)
         IMPLICIT NONE
         INTEGER NMEM(*)
         INTEGER iCount
         iCount=NMEM(23)
         RETURN
      END SUBROUTINE GetSymExcitCount

C.. IF(TSETUP) Generate an iterator which allows up to double excitations to be generated
C.. one at a time (in an unordered fashion) from a given det.  THis needs to be called twice,
C.. first with STORE(1)=0, and then again.  Finally it can be called with TSETUP=.FALSE. 
C.. to actually generate the excitations
!! First run.  TSETUP=.TRUE.  STORE(1)=0.  NMEM can just be an INTEGER
!!  requires INTEGER STORE(6) to hold initialization data.  NMEM is returned containing the memory required to hold the excitation generator (call this nLength)
!! Second run,  TSETUP=.TRUE.  STORE(:) contains data from first run.  NMEM now is an array of nLenght INTEGERs.  The excitation generator data is stored in this.
!! Third run,   TSETUP=.FALSE.  will generate an excitation of NI, and put it in NEL.  

      SUBROUTINE GENSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER G1(nBasis),NBASISMAX(*)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER ICOUNT
         INTEGER NMEM(*),NJ(NEL),IC,IFRZ
         LOGICAL TSETUP
         INTEGER ILEVEL
         CALL GenSymExcitIt2Par(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL,1,NEL)
      END

!GenSymExcitIt3Par and GenSymExcitIt3 are versions GenSymExcitIt2 which do not require information taken from System and Integrals modules to be specified on the command line
      SUBROUTINE GENSYMEXCITIT3(NI,TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         IMPLICIT NONE 
         INCLUDE 'irat.inc'
         INTEGER NI(NEL)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER*8 STORE2(6)
         INTEGER ICOUNT
         INTEGER NMEM(*),NJ(NEL),IC,IFRZ,I
         LOGICAL TSETUP
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         INTEGER ILEVEL

         INTEGER iMinElec1, iMaxElec1
         CALL GENSYMEXCITIT2Par(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL,1,nEl)
      END

      SUBROUTINE GENSYMEXCITIT3Par(NI,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL,iMinElec1,iMaxElec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         IMPLICIT NONE 
         INCLUDE 'irat.inc'
         INTEGER NI(NEL)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER*8 STORE2(6)
         INTEGER ICOUNT
         INTEGER NMEM(*),NJ(NEL),IC,IFRZ,I
         LOGICAL TSETUP
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         INTEGER ILEVEL

         INTEGER iMinElec1, iMaxElec1
         CALL GENSYMEXCITIT2Par(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL,iMinElec1,iMaxElec1)
      END


C.. A modified version of GENSYMEXCITIT2.  This has min and max electron 1 specifiers used for parallelization
      SUBROUTINE GENSYMEXCITIT2Par(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE,ILEVEL,iMinElec1,iMaxElec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER G1(nBasis),NBASISMAX(*)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER*8 STORE2(6)
         INTEGER ICOUNT
         INTEGER NMEM(*),NJ(NEL),IC,IFRZ,I
         LOGICAL TSETUP
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         INTEGER ILEVEL

         INTEGER iMinElec1, iMaxElec1
         IF(TSETUP) THEN
C.. This is the first time we've been called for setup.
!  We pass information back in STORE, but actually hold it in STORE2 during the SYMSETUPEXCITS2
            IF(STORE(1).EQ.0) THEN
               CALL IAZZERO(STORE2,6*2)
               CALL MEMORY(IP_DSTORE,SymClassSize*NEL+NBASIS/32+1
     &                     +SymmetrySize*(NEL*NEL+1),'DSTORE')
               STORE2(1)=IP_DSTORE
C.. Just count.
               CALL SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE2,
     &            .TRUE.,ICOUNT,DSTORE(1), DSTORE(SymClassSize*NEL+1),
     &          DSTORE(SymClassSize*NEL+1+NBASIS/32+1),ILEVEL,iMinElec1,
     &           iMaxElec1)

               CALL FREEM(IP_DSTORE)
               DO i=1,6
                  STORE(i)=STORE2(i)
               ENDDO
               STORE(1)=SymClassSize*NEL+NBASIS/32+1+
     &            SymmetrySize*(NEL*NEL+1)
               NMEM(1)=23+STORE(1)+STORE(2)+STORE(3)+STORE(4)+STORE(5)
C               WRITE(6,"(A6,6I5)"), "SIZE:",NMEM(1),(STORE(I),I=1,5)
            ELSE
C.. The second setup.  Now NMEM is allocated, we store all the info
C.. NMEM is as follows:
C..   1           -  5           STORE(1:5)
C..   6           -  6           STORE(6) = NEXCITTYPES
C.. Data for the iterators
C..   7  IEXCIT
C..   8  ISPN
C..   9  IFROM
C..   10 ITO
C..   11 I
C..   12 J
C..   13 K
C..   14 L
C..   15 ICC(1:4)
C..   19 LS(1:2,1:2)  ([1 or 2], [A or B])
C..   23 NEXCIT      total number of excitations
C..  (STORE(1)=24)-  STORE(2)-1  DSTORE
C..   STORE(2)    -  STORE(3)-1  EXCITTYPES
C..   STORE(3)    -  STORE(4)-1  NALLOWPPS
C..   STORE(4)    -  STORE(5)-1  ORBPAIRS
C..   STORE(5)    -  ...         SYMPRODIND

C..   DSTORE(1)   -  DSTORE(NEL*SymClassSize) CLASSES
C..   DSTORE(NEL*SymClassSize+1) - ,,,        ILUT

               CALL IAZZERO(NMEM,23)
               ICOUNT=24
C.. Put the addresses in store
               DO I=1,5
                  STORE2(I)=LOC(NMEM(ICOUNT))
                  ICOUNT=ICOUNT+STORE(I)
               ENDDO
               CALL IAZZERO(NMEM(6),ICOUNT-6)
               STORE2(6)=0
               NMEM(11)=-1
               NMEM(7)=0
               IP_DSTORE=STORE2(1)
C               CALL DUMPIMEMORY(6,NMEM,ICOUNT-1)
!!      SUBROUTINE SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
!!     &   TCOUNT,ICOUNT,CLASSES,ILUT,SYMPRODS,ILEVEL)
               CALL SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE2,
     &            .FALSE.,ICOUNT,DSTORE(1), DSTORE(SymClassSize*NEL+1),
     &          DSTORE(SymClassSize*NEL+1+NBASIS/32+1),ILEVEL,iMinElec1,
     &          iMaxElec1)
               NMEM(6)=STORE2(6)
               NMEM(23)=ICOUNT
               IC=ICOUNT
C.. and now instead the indices within NMEM
               ICOUNT=24
               DO I=1,5
                  NMEM(I)=ICOUNT
                  ICOUNT=ICOUNT+STORE(I)
               ENDDO
C.. Second setup finally complete.
C               CALL DUMPIMEMORY(6,NMEM,ICOUNT-1)
            ENDIF
         ELSE
C.. Actually generate a det
C            WRITE(6,"(A,Z10,8I4)",advance='no') "GET",LOC(NMEM(1)),
C     &         (NMEM(I),I=7,14)
            IP_DSTORE=LOC(NMEM(NMEM(1)))
!!      SUBROUTINE SYMGENEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
!!     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,
!!     &               I,J,K,L,ICC,LS,
!!     &               NK,IC)
            CALL SYMGENEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NJ,IC,iMinElec1,
     &         iMaxElec1)
C            CALL WRITEDET(6,NI,NEL,.FALSE.)
C            WRITE(6,"(A)",advance='no'), "->"
C            IF(NJ(1).NE.0) THEN
C               CALL WRITEDET(6,NJ,NEL,.TRUE.)
C            ELSE 
C               WRITE(6,*) "(    0,)"
C            ENDIF
         ENDIF
      END
      SUBROUTINE RESETEXIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &            NMEM,IFRZ)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,NBASIS,G1(nBasis),NBASISMAX(*)
         INTEGER NMEM(*),IFRZ(*)
!         CALL IAZZERO(NMEM,22)
         NMEM(11)=-1  ! Set I=-1 (i.e. start from the first electron again)
         NMEM(7)=0  !set iExcit=0 (i.e. start from the first excit class again
      END
C                  CALL RESETEXIT(IPATH(1,LOCTAB(3,IVLEVEL)),NEL,G1,
C     &               NBASIS,NBASISMAX,CURGEN,IFRZ2)

      SUBROUTINE SYMGENEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,
     &               I,J,K,L,ICC,LS,
     &               NK,IC,iMinElec1,iMaxElec1)
      use SystemData , only : TSTOREASEXCITATIONS
      use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
      use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
      use SymData, only: SymLabelCounts,SymStatePairs,SymClass
      use SymData, only: SymLabelList,SymPairProds
      IMPLICIT NONE
      INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES),NEXCITTYPES
      
      INTEGER I,J,K,L
      INTEGER NK(NEL)
      TYPE(SymClass) CLASSES(*)
      INTEGER IEXCIT
      INTEGER IFROM,ITO,ISPN
      INTEGER IFROMSL,ITOSL
      TYPE(Symmetry) SPP
      LOGICAL L1,L2
      INTEGER ICC(4)
C.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.
      LOGICAL LS(2,2)

      INTEGER ORBPAIRS(2,*)
      INTEGER ILUT(0:*)
      INTEGER SYMPRODIND(2,3,0:*)
      INTEGER IC
      LOGICAL tDebugPrint
      INTEGER iMinElec1, iMaxElec1
      tDebugPrint=.false.

      DO WHILE(.TRUE.)
C.. see if we need a new EXCIT
!         WRITE(6,*) I,iExcit,nExcitTypes
         IF(I.LT.0) THEN
C.. move to the next excitation
            IEXCIT=IEXCIT+1
            IF(IEXCIT.GT.NEXCITTYPES) THEN
C.. We're done
               NK(1)=0
               RETURN
            ENDIF
            IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
C.. a single
               ISPN=EXCITTYPES(2,IEXCIT)-2
               IFROM=CLASSES(EXCITTYPES(3,IEXCIT))%SymLab
               ITO=EXCITTYPES(4,IEXCIT)
               L=0
               K=SYMLABELCOUNTS(2,ITO)
!  J will run along symblabels
               J=0
!  I will run along the determinant
               I=0
!  Take into account the we may only be doing a subset of the electrons
               DO WHILE(I.LT.NEL)
                  IF(NI(I+1).GE.iMinElec1) EXIT
                  I=I+1
               ENDDO
!  I now is the index of the first electron allowed in our subset
               IF(tDebugPrint)
     &         WRITE(6,*) "EXC",IEXCIT,ISPN,IFROM,ITO,SYMPRODIND(2,ISPN,
     &            IFROM)
            ELSE
C.. a double
               ISPN=EXCITTYPES(2,IEXCIT)
               IFROM=EXCITTYPES(3,IEXCIT)
               ITO=EXCITTYPES(4,IEXCIT)
               I=0
               L=5
               K=-2
               IF(tDebugPrint)
     &         WRITE(6,*) "EXC",IEXCIT,ISPN,IFROM,ITO,SYMPRODIND(2,ISPN,
     &            IFROM)
            ENDIF
         ENDIF
         IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
C.. singles
C.. We always need a new K.  K and L run in parallel
            K=K+1
C.. We've stored IFROMSL in ICC1
            IFROMSL=ICC(1)
C.. For each K, we check to see if it's valid
            IF(K.GE.SYMLABELCOUNTS(2,ITO)) THEN
C.. no more possible K, so we get a new I
               K=-1
               L=0
C.. I and J run in parallel, and for each J, we check whether it's in 
C.. the det.  If it is, then we carry on.
C.. SYMLABELLIST holds a list of states grouped under symmlabel, and ordered within that symlabel
C.. SYMLABELCOUNTS(1,J) is the index within SYMLABELLIST of the first state of symlabel J
C.. SYMLABELCOUNTS(2,J) is the number of states with symlabel J
!..  J runs along the symlabel looking for states of that symmetry
!..  I runs along the determinant comparing against electrons of that symmetry.
!..  We're looking for electrons in the det which are also in the SYMLABELLIST, and we run I and J
!..   in parallel along each to keep the search efficient.
!  instead of NEL, we have iMaxElec1
               L2=I.LT.iMaxElec1
               L1=.TRUE.
!  This loop, dependent on L2 gets a new I
               DO WHILE(L2)
                  I=I+1
                  L1=J.LT.SYMLABELCOUNTS(2,IFROM)
!  This loop dependent on L1 gets a new J
                  DO WHILE(L1)
                     IFROMSL=(SYMLABELLIST(SYMLABELCOUNTS(1,IFROM)+J)*2
     &                     +ISPN)
                     ICC(1)=IFROMSL
C                     WRITE(6,*) I,NI(I),J,IFROM,IFROMSL
!  If J points to an electron lower than the current in the det, inc J
                     IF(IFROMSL.LT.NI(I)) J=J+1
!  If J was a lower elec, and there are more J to get, then carry on in this loop
                     L1=IFROMSL.LT.NI(I)
     &                     .AND.J.LT.SYMLABELCOUNTS(2,IFROM)
                  ENDDO
C.. see if we need a new I because NI(I)<IFROMSL which is the current electron corresponging to J.
!  Instead of NEL we have iMaxElec1
                  L2=IFROMSL.NE.NI(I).AND.I.LT.iMaxElec1
!  If I.EQ.NEL, then we're run out of electrons, because the next electron we get will be in pos NEL+1
               ENDDO
C.. If we've gone too far, signal a new excit
!  If L1 is true then we never made it into the get new J loop because L2 was FALSE because there were no more electrons to get in the det. Alternatively we've hit the last electron in the det, and not found a correspondence in the sym excit list.  Either way we signal to give up on this set of excitations, as we cannot find a valid I.
!  Instead of NEL we use iMaxElec1
               IF(L1.OR.(I.EQ.iMaxElec1.AND.IFROMSL.NE.NI(I))) I=-1

!  Cycle back to the beginning with our new valid I (or invalid I signalling find a new excit type), and look for a K
               CYCLE
            ENDIF
C.. it's a valid K, but is it in the det already?
C.. ITOSL is the orb it corresponds to
            ITOSL=2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)
     &                              +ISPN
C.. Check if it's in the det
            L2=.TRUE.
C            WRITE(6,*) SYMLABELCOUNTS(2,ITO),ITOSL            
            DO WHILE (L2)
               L=L+1
               L2=L.LE.NEL.AND.NI(L).LT.ITOSL
            ENDDO
            IF(L.LE.NEL.AND.NI(L).EQ.ITOSL) THEN
C.. We've found an L in the det which is the same as ITOSL. 
C.. we go round again, getting another K
               CYCLE
            ENDIF
            L=L-1
C.. hoorah!  We've got an I in the det, and a K not in the det.  Create an excitation
            IF(tStoreAsExcitations) THEN
!The excitation storage starts with -1.  The next number is the excitation level,L .  Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
               NK(1)=-1
               NK(2)=1
               NK(3)=0
               NK(4)=NI(I)
               NK(5)=ITOSL
               STOP 'tStoreAsExcitations not tested for singles.'
            ELSE
               CALL ICOPY(NEL,NI,1,NK,1)
               IF(tDebugPrint)
     &        WRITE(6,*) "[",NK(I),"->",ITOSL,"]"
               NK(I)=ITOSL
               CALL SORTI(NEL,NK)
               IC=1
            ENDIF
C.. quit the do loop
            EXIT
         ELSE
C.. doubles
            SPP=SymPairProds(ITO)%Sym
C.. See if we need a new K.  L is the spin label, and goes from 1..4.
            IF(L.GT.4) THEN
               L=1
               K=K+1
               IF(K.LT.0.OR.
     &           K.GE.SymPairProds(ITO)%nIndex+SymPairProds(ITO)%nPairs)
     &            THEN
C.. K is invalid.  Get a new I
                  K=1
                  I=I+1
                  IF(I.GT.SYMPRODIND(2,ISPN,IFROM)) THEN
C.. I is now invalid
                     I=-1
                     CYCLE
                  ELSE
C.. reset K
                     K=SymPairProds(ITO)%nIndex
                  ENDIF
               ENDIF
C.. We've got a new K, so we need to reset some variables
               ICC(1)=SymStatePairs(1,K)*2-1
               ICC(2)=ICC(1)+1
               ICC(3)=SymStatePairs(2,K)*2-1
               ICC(4)=ICC(3)+1
C                     WRITE(6,*) ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
C     &                      ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)
C                     WRITE(6,*) ISPN,ICC(1),ICC(3)
               LS(1,1)=BTEST(ILUT((ICC(1)-1)/32),
     &                     MOD(ICC(1)-1,32))
               LS(1,2)=BTEST(ILUT((ICC(2)-1)/32),
     &                     MOD(ICC(2)-1,32))
               LS(2,1)=BTEST(ILUT((ICC(3)-1)/32),
     &                     MOD(ICC(3)-1,32))
               LS(2,2)=BTEST(ILUT((ICC(4)-1)/32),
     &                     MOD(ICC(4)-1,32))
            ENDIF
C.. Now check for an excitation
            IF(L.EQ.1) THEN
               L=2
               IF(ISPN.EQ.1.AND.ICC(1).NE.ICC(3).AND.
     &            .NOT.(LS(1,1).OR.LS(2,1))) THEN
                  IF(tStoreAsExcitations) THEN
!The excitation storage starts with -1.  The next number is the excitation level,L .  Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
                     NK(1)=-1
                     NK(2)=2
                     NK(3)=1
                     NK(4)=ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(5)=ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(6)=ICC(1)
                     NK(7)=ICC(3)
                  ELSE
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(1),","
                           NK(J)=ICC(1)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(3),"]"
                           NK(J)=ICC(3)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
                     CALL SORTI(NEL,NK)
                  ENDIF
                  IC=2
                  EXIT
               ENDIF
            ENDIF
            IF(L.EQ.2) THEN
               L=3
               IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
                  IF(.NOT.(LS(1,1).OR.LS(2,2))) THEN
                   IF(tStoreAsExcitations) THEN
!The excitation storage starts with -1.  The next number is the excitation level,L .  Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
                     NK(1)=-1
                     NK(2)=2
                     NK(3)=1
                     NK(4)=ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(5)=ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(6)=ICC(1)
                     NK(7)=ICC(4)
                   ELSE
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(1),","
                           NK(J)=ICC(1)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(4),"]"
                           NK(J)=ICC(4)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
                     CALL SORTI(NEL,NK)
                   ENDIF
                   IC=2
                   EXIT
                  ENDIF
               ENDIF
            ENDIF
            IF(L.EQ.3) THEN
               L=4
               IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
                 IF(.NOT.(LS(1,2).OR.LS(2,1)).AND.ICC(1).NE.ICC(3)) THEN
                   IF(tStoreAsExcitations) THEN
!The excitation storage starts with -1.  The next number is the excitation level,L .  Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
                     NK(1)=-1
                     NK(2)=2
                     NK(3)=1
                     NK(4)=ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(5)=ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(6)=ICC(2)
                     NK(7)=ICC(3)
                   ELSE
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(2),","
                           NK(J)=ICC(2)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(3),"]"
                           NK(J)=ICC(3)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
                     CALL SORTI(NEL,NK)
                   ENDIF
                   IC=2
                   EXIT
                  ENDIF
               ENDIF
            ENDIF
            IF(L.EQ.4) THEN
               L=5
               IF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                 IF(ICC(1).NE.ICC(3).AND..NOT.(LS(1,2).OR.LS(2,2))) THEN
                   IF(tStoreAsExcitations) THEN
!The excitation storage starts with -1.  The next number is the excitation level,L .  Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
                     NK(1)=-1
                     NK(2)=2
                     NK(3)=1
                     NK(4)=ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(5)=ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)     
                     NK(6)=ICC(2)
                     NK(7)=ICC(4)
                   ELSE
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(2),","
                           NK(J)=ICC(2)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(4),"]"
                           NK(J)=ICC(4)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
                     CALL SORTI(NEL,NK)
                   ENDIF
                   IC=2
                   EXIT
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO 
      END
C.. Generate a random symmetry allowed excitation.  GENSYMEXCITIT2 must be called twice
C.. before this.
C.. This version is based on GenRandSymExcitIt2, but does not need values it extracts from modules
!  IC is the excitation level, and iCount the total number of allowed excitations
      SUBROUTINE GenRandSymExcitIt3(NI,NMEM,NJ,ISEED,IC,IFRZ,PGEN,
     &      ICOUNT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         use SystemData, only: Arr,nMax
         use IntegralsData, only: UMat
         use SymData, only: SymClassSize
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NI(NEL)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER STORE(6)
         INTEGER NMEM(*),NJ(NEL),IC,ICOUNT,IFRZ,I
         INTEGER ISEED
         REAL*8 PGEN
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
C.. The second setup.  Now NMEM is allocated, we store all the info
C.. NMEM is as follows:
C..   1           -  5           STORE(1:5)
C..   6           -  6           STORE(6) = NEXCITTYPES
C.. Data for the iterators
C..   23 NEXCIT      total number of excitations
C..  (STORE(1)=24)-  STORE(2)-1  DSTORE
C..   STORE(2)    -  STORE(3)-1  EXCITTYPES
C..   STORE(3)    -  STORE(4)-1  NALLOWPPS
C..   STORE(4)    -  STORE(5)-1  ORBPAIRS
C..   STORE(5)    -  ...         SYMPRODIND

C..   DSTORE(1)   -  DSTORE(NEL*SymClassSize) CLASSES
C..   DSTORE(NEL*SymClassSize+1) - ...       ILUT
C.. Actually generate a det
         IP_DSTORE=LOC(NMEM(NMEM(1)))
         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
         ICOUNT=NMEM(23)
      END
C.. Generate a random symmetry allowed excitation.  GENSYMEXCITIT2 must be called twice
C.. before this.
C.. ICOUNT returned is the total number of possible excitations.  To get the excitation level,
C..   use GenRandSymExcitIt3
      SUBROUTINE GENRANDSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         NMEM,NJ,ISEED,ICOUNT,IFRZ,UMAT,ARR,PGEN)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER G1(nBasis),NBASISMAX(*)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER STORE(6)
         INTEGER NMEM(*),NJ(NEL),ICOUNT,IFRZ,I,IC
         LOGICAL TSETUP
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         INTEGER ISEED
         REAL*8 UMAT(*),PGEN
         REAL*8 ARR(NBASIS,2)
C.. The second setup.  Now NMEM is allocated, we store all the info
C.. NMEM is as follows:
C..   1           -  5           STORE(1:5)
C..   6           -  6           STORE(6) = NEXCITTYPES
C.. Data for the iterators
C..   23 NEXCIT      total number of excitations
C..  (STORE(1)=24)-  STORE(2)-1  DSTORE
C..   STORE(2)    -  STORE(3)-1  EXCITTYPES
C..   STORE(3)    -  STORE(4)-1  NALLOWPPS
C..   STORE(4)    -  STORE(5)-1  ORBPAIRS
C..   STORE(5)    -  ...         SYMPRODIND

C..   DSTORE(1)   -  DSTORE(NEL*SymClassSize) CLASSES
C..   DSTORE(NEL*SymClassSize+1) - ...       ILUT
C.. Actually generate a det
         IP_DSTORE=LOC(NMEM(NMEM(1)))
         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
         ICOUNT=NMEM(23)
      END
      SUBROUTINE SYMGENRANDEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,       &
     &               CLASSES,                                           &
     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,    &
     &               I,J,K,L,ICC,LS,ITOTAL,ISEED,                       &
     &               NK,IC,G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)      
      USE HElem
      USE symexcit2
      use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
      use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
      use SymData, only: SymClass,SymLabelCounts,SymLabelList
      INTEGER NBASIS
      INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES),NEXCITTYPES
      
      INTEGER I,J,K,L
      INTEGER NK(NEL)
      TYPE(SymClass) CLASSES(*)
      INTEGER IEXCIT
      INTEGER IFROM,ITO,ISPN
      INTEGER IFROMSL,ITOSL
      TYPE(Symmetry) SPP
      LOGICAL L1,L2
      INTEGER ICC(4)
!.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.
      LOGICAL LS(2,2)

      INTEGER ORBPAIRS(2,*)
      INTEGER ILUT(0:*),ITOTAL,ISEED
      INTEGER SYMPRODIND(2,3,0:*)
      INTEGER IC,ICOUNT
      REAL*8 RAN2
         INTEGER ICC1,ICC2,ICC3,ICC4
         LOGICAL L1B,L1A,L2B,L2A
      LOGICAL ISVALIDDET
      TYPE(BasisFN) G1(nBasis)
      INTEGER nBasisMax(5,*)
      TYPE(HElement) UMAT(*)
      REAL*8 RCOUNT,R
      REAL*8 PGEN
      TYPE(ExcitWeight) ews(*)
      POINTER (IP_ExcitWeight,ews)
      REAL*8 ARR(NBASIS,2)
!.. First generate a random excitation 
      PGEN=1.D0/ITOTAL
      I=ITOTAL*RAN2(ISEED)
      IEXCIT=1
      DO WHILE(EXCITTYPES(5,IEXCIT).LE.I)
         I=I-EXCITTYPES(5,IEXCIT)
         IEXCIT=IEXCIT+1
      ENDDO
!.. We've chosen the excitation.  Now pick a random from and to
      IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
!.. a single
         ISPN=EXCITTYPES(2,IEXCIT)-2
         IFROM=CLASSES(EXCITTYPES(3,IEXCIT))%SymLab
         ITO=EXCITTYPES(4,IEXCIT)
         ICOUNT=0
!  K is a loop variable.  We go through the loop once to count the number of possible excitations
!  We then pick a random number up to that number, and choose that excitation
         ot:DO K=0,1
            J=0
            I=1
!               WRITE(6,*) ISPN,IFROM,ITO,EXCITTYPES(5,IEXCIT)
!.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
!.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I
            DO WHILE (I.LE.NEL)
               L1=J.LT.SYMLABELCOUNTS(2,IFROM)
               DO WHILE(L1)
                  IFROMSL=(SYMLABELLIST(SYMLABELCOUNTS(1,IFROM)+J)*2
     &                        +ISPN)
!                     WRITE(6,*) I,NI(I),J,IFROM,IFROMSL
                  IF(IFROMSL.LT.NI(I)) J=J+1
                  L1=IFROMSL.LT.NI(I).AND.J.LT.SYMLABELCOUNTS(2,IFROM)
               ENDDO
!                  WRITE(6,*) I,J,SYMLABELCOUNTS(2,IFROM)
               IF(J.LT.SYMLABELCOUNTS(2,IFROM).AND.IFROMSL.EQ.NI(I)) 
     &            THEN
!.. We've found an orb in NI with the correct sym. 
                  IF(K.EQ.0) THEN
                     ICOUNT=ICOUNT+1
                  ELSEIF(ICOUNT.EQ.0) THEN
                     EXIT ot
                  ELSE
                     ICOUNT=ICOUNT-1
                  ENDIF
               ENDIF
               I=I+1
            ENDDO
            ICOUNT=ICOUNT*RAN2(ISEED)
         ENDDO ot
!.. NI(I) is now the orbital to excite
! Now go through the list of possible excitations of it
         ICOUNT=0
!  J is a loop variable.  We go through the loop once to count the number of possible excitations
!  We then pick a random number up to that number, and choose that excitation
         DO J=0,1
            K=0
            L=0
            DO WHILE(K.LT.SYMLABELCOUNTS(2,ITO))
               ITOSL=2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)+ISPN
               L2=L.EQ.0.OR.NI(L).LT.ITOSL
               DO WHILE (L2)
                  L=L+1
                  L2=L.LE.NEL.AND.NI(L).LT.ITOSL
               ENDDO
               IF(L.GT.NEL.OR.NI(L).NE.ITOSL) THEN
                  IF(J.EQ.0) THEN
                     ICOUNT=ICOUNT+1
                  ELSEIF(ICOUNT.EQ.0) THEN
!.. We've found a virtual into which we can excite our occupied orb.
                     CALL ICOPY(NEL,NI,1,NK,1)
                     NK(I)=ITOSL
                     CALL SORTI(NEL,NK)
                     IF(.NOT.ISVALIDDET(NK,NEL)) THEN
            
                       WRITE(6,*) "INVALID DET SYMGENRANDEXCITIT SINGLE"
                       CALL WRITEDET(6,NK,NEL,.TRUE.)
                       CALL WRITEDET(6,NI,NEL,.TRUE.)
                       WRITE(6,*) I,ITOSL,K,L
                       WRITE(6,*)                                       &
     &                  (2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)+ISPN,  &
     &                     K=0,SYMLABELCOUNTS(2,ITO)-1)
                        STOP "INVALID DET SYMGENRANDEXCITIT SINGLE"
                     ENDIF
                     IC=1
                     RETURN
                  ELSE
                     ICOUNT=ICOUNT-1
                  ENDIF
               ENDIF
               K=K+1
            ENDDO
            ICOUNT=ICOUNT*RAN2(ISEED)
         ENDDO
      ELSE
!.. a double excitation
!.. We store each excitation type as:
!.. 1   TYPE (single=1, double=2)
!.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
!.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
!.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
!.. 5  COUNT (Total number of excitations in this category)
         ISPN=EXCITTYPES(2,IEXCIT)
         IFROM=EXCITTYPES(3,IEXCIT)
         ITO=EXCITTYPES(4,IEXCIT)
!.. Go through the list of pairs with a given symprod.
!.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
!.. SYMPRODS(I) in ORBPAIRS
!.. SYMPRODIND(2,ISPN,I) contains the number of such elements
         pGen=(0.D0+ExcitTypes(5,iExcit))/iTotal !pGen is the prob of choosing iExcit
!  Allocate memory to store all weights and excitations in this category
         CALL MEMORY(IP_ExcitWeight,ExcitWeightSize*
     &      SYMPRODIND(2,ISPN,IFROM),'ExcitWeigh')
         CALL EnumExcitFromWeights(EXCITTYPES(1,IEXCIT),ews,OrbPairs,
     &      SymProdInd,
     &      RCOUNT,ICOUNT,G1,NBASISMAX,UMAT,ARR,NBASIS)
!  We choose a symprod of this symmetry to excite FROM with probs in ews
         pGen=pGen/rCount  ! The norm factor for the FROMs
         rCount=rCount*RAN2(ISEED)
         I=0
!  This is a little inefficient, and could perhaps be binary searched if cumulative weights were also stored
         DO WHILE (rCount.GT.0.D0)
            I=I+1
            rCount=rCount-ews(I)%Weight
         ENDDO   
!  I is the index of the specific SYM PAIR FROM
         pGen=pGen*ews(I)%Weight
         
!  pGen is the prob of choosing a specific FROM (given having chosen iExcit proportional to the number of excitations in each iExcit)
!           times the prob of choosing iExcit
         CALL FREEM(IP_ExcitWeight)


!  Allocate memory to store all weights and excitations in this category
         CALL MEMORY(IP_ExcitWeight,ExcitWeightSize*
     &      EXCITTYPES(5,IEXCIT)/SYMPRODIND(2,ISPN,IFROM),'ExcitWeigh')

!.. Now go through the list of virtual pairs, excluding those with orbitals in NI 
!.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
!.. states (listed in SymStatePairs), and the index of the start of this list
!.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
!.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose prod is of that symmetry.
         RCOUNT=0.D0
         ICOUNT=0
!  L is a loop variable.  We go through the loop once to count the number of possible excitations
!  We then pick a random number up to that number, and choose that excitation
!  This is modified to weight each excitation according to its U-matrix element, so each
!  Excitation's weight is added to RCOUNT to give a normalization.  RCOUNT*RAND then corresponds to 
!  a point in the sum of weights which is found by subtracting weights until we reach 0.
!  This is only done for the excitations TO because we need to know the excitations FROM to generate the weight
         CALL EnumExcitWeights(EXCITTYPES(1,IEXCIT),I,ILUT,ews,OrbPairs,
     &      SymProdInd,RCOUNT,ICOUNT,G1,NBASISMAX,UMAT,ARR,NBASIS)

         pGen=pGen/rCount           ! Normalize for the sum of TOs
         
!  We've now got a list of iCount weights in ews, and total norm rCount.
!  Select a random one
         rCount=rCount*RAN2(ISEED)
         K=0
!  This is a little inefficient, and could perhaps be binary searched if cumulative weights were also stored
         DO WHILE (rCount.GT.0.D0)
            K=K+1
            rCount=rCount-ews(K)%Weight
         ENDDO   
!  ews(K) is the excitation we've chosen
         PGEN=PGEN*ews(K)%Weight !weight with the weight of this TO
!  pGen is the prob of choosing a specific TO (given the FROM, and the iExcit)
!           times prob of choosing a specific FROM (given having chosen iExcit proportional to the number of excitations in each iExcit)
!           times the prob of choosing iExcit
         CALL ICOPY(NEL,NI,1,NK,1)
         DO J=1,NEL
            IF(NI(J).EQ.ews(K)%I) THEN
               NK(J)=ews(K)%A
            ENDIF
            IF(NI(J).EQ.ews(K)%J) THEN
               NK(J)=ews(K)%B
            ENDIF
         ENDDO
         CALL SORTI(NEL,NK)
         IC=2
         CALL FREEM(IP_ExcitWeight)
         RETURN
      ENDIF
      END 
      LOGICAL FUNCTION IsConnectedDet(nI,nJ,nEl,nIExcitor,G1,nBasisMax,
     &      nBasis)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize 
         IMPLICIT NONE
         INTEGER nEl,nI(nEl),nJ(nEl),nIExcitor(*),nBasis,nBasisMax(*)
         TYPE(BasisFN) G1(nBasis) 
         INTEGER DStore(*)
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         POINTER(IP_DSTORE,DSTORE)
         IP_DSTORE=LOC(nIExcitor(nIExcitor(1)))
         CALL IsConnectedDet2(nI,nJ,nEl,G1,nBasisMax,nBasis,
     &      nIExcitor(nIExcitor(4)),nIExcitor(nIExcitor(5)),
     &      DStore(SymClassSize*nEl+1),
     &      DStore(SymClassSize*nEl+1+nBasis/32+1),
     &      nIExcitor(nIExcitor(2)),nIExcitor(23),IsConnectedDet)
         RETURN
      END

!  We are required to enshroud GenExcitProbInternal as F90 doesn't seem to want to allow any casting.
      SUBROUTINE IsConnectedDet2(nI,nJ,nEl,G1,nBasisMax,
     &      nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,
     &      tIsConnected)
         USE HElem
         USE symexcit2
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER nEl,nI(nEl),nJ(nEl),nBasis,nBasisMax(5,*)
         INTEGER ExcitTypes(5,*)
         TYPE(BasisFN) G1(nBasis) 
         INTEGER iLUT(*),iTotal
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,0:*)
         TYPE(Symmetry) SymProds(0:*)
         LOGICAL tIsConnected
         CALL IsConnectedDetInternal(nI,nJ,nEl,G1,nBasisMax,nBasis,
     & OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,tIsConnected)
      END

!We wish to calculate what excitation class the excitation NI->NJ falls into, with the appropriate 
! IFROM class and index within that, IFROMINDEX, and ITO class, and index within that. ITOINDEX
! After that, we generate the probability that nJ would be an excitation from nI.
! We require the excitation generator of I, which is used to pass info to GenExcitProbInternal (in symexcit2.F90)
      SUBROUTINE GenExcitProb(nI,nJ,nEl,nIExcitor,G1,nBasisMax,UMat,Arr,
     &      nBasis,pGen)
         USE HElem
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         IMPLICIT NONE
         INTEGER nEl,nI(nEl),nJ(nEl),nIExcitor(*),nBasis,nBasisMax(*)
         TYPE(BasisFN) G1(nBasis) 
         TYPE(HElement) UMat(*)
         INTEGER DStore(*)
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         POINTER(IP_DSTORE,DSTORE)
         REAL*8 pGen
         REAL*8 Arr(nBasis,2)
         IP_DSTORE=LOC(nIExcitor(nIExcitor(1)))
         CALL GenExcitProb2(nI,nJ,nEl,G1,nBasisMax,UMat,Arr,nBasis,
     &      nIExcitor(nIExcitor(4)),nIExcitor(nIExcitor(5)),
     &      DStore(SymClassSize*nEl+1),
     &      DStore(SymClassSize*nEl+1+nBasis/32+1),
     &      nIExcitor(nIExcitor(2)),nIExcitor(23),pGen)
      END
! For reference below is the call and declaration of SYMGENRANDEXCITIT to extract the locations of teh variables.
!  NMEM is nIExcitor
!         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
!     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
!     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
!     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
!     &         G1,NBASISMAX,UMAT,NBASIS,PGEN)
!         IC=NMEM(23)
!      END
!     SymProds(0) = DSTORE(SymClassSize*NEL+1+nBasis/32+1))
!      SUBROUTINE SYMGENRANDEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,       &
!     &               CLASSES,                                           &
!     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,    &
!     &               I,J,K,L,ICC,LS,ITOTAL,ISEED,                       &
!     &               NK,IC,G1,NBASISMAX,UMAT,NBASIS,PGEN)               


!  We are required to enshroud GenExcitProbInternal as F90 doesn't seem to want to allow any casting.
      SUBROUTINE GenExcitProb2(nI,nJ,nEl,G1,nBasisMax,UMat,Arr,
     &      nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,
     &      pGen)
         USE HElem
         USE symexcit2
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER nEl,nI(nEl),nJ(nEl),nBasis,nBasisMax(5,*)
         INTEGER ExcitTypes(5,*)
         TYPE(BasisFN) G1(nBasis) 
         TYPE(HElement) UMat(*)
         INTEGER iLUT(*),iTotal
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,0:*)
         TYPE(Symmetry) SymProds(0:*)
         REAL*8 pGen
         REAL*8 Arr(nBasis,2)
         CALL GenExcitProbInternal(nI,nJ,nEl,G1,nBasisMax,UMat,Arr,
     &  nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,pGen)
      END
!  Calculate the generation probability of a given graph for debugging purposes.
!  if iUnit<>0, write this to this unit
      SUBROUTINE CalcWriteGraphPGen(iUnit,iGraph,iV,nEl,LOCTAB,G1,
     &            nBasisMax,UMat,Arr,nBasis,pGenGraph,EXCITGEN)
         USE symexcit2
         USE HElem
         USE PrecalcData , only : PRE_TAY,TPREVAR
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER iUnit
         INTEGER iV,nEl       ! Vertices and number of electrons
         INTEGER iGraph(nEl,0:iV)  ! The graph
         INTEGER LOCTAB(3,iV)    ! Three components for each vertex in the graph: 
C.. LOCTAB(1,V) is the address of the generator used to create node V in
C.. the path (the first node is node 0).
C.. LOCTAB(2,V) is the length of the generator (i.e. the amount of memory used to store it)
C.. LOCTAB(3,V) is the vertex from which vertex V was were excited
C.   The last node never has an excitation generator, as it never needs one in the full graph generation scheme.
C.    We must create one.
         INTEGER STORE(6)
         POINTER (IP_NEWEX,NEWEX)
         INTEGER NEWEX(*),NEWEXLEN
         !Array of pointers
#if defined(POINTER8)
         INTEGER*8 EXCITGEN(0:iV)
#else
         INTEGER EXCITGEN(0:iV)
#endif

         TYPE(BasisFN) G1(nBasis)
         INTEGER nBasisMax(*)
         TYPE(HElement) UMat(*)
         REAL*8 Arr(nBasis,2)
         INTEGER nBasis,IC,INODE2(NEL)
         REAL*8 pGen, pGenGraph

         REAL*8 XIJ(0:iV-1,0:iV-1)  ! XIJ(I,J) is the pgen of J from I
         REAL*8 GETPATHPROB2
         INTEGER i,j,ISEED
         ISEED=0
         
         IF(TPREVAR.AND.((PRE_TAY(1,iV).EQ.-7).OR.                      
     &                       (PRE_TAY(1,iV).EQ.-19))) THEN

            DO i=0,iV-1
                IP_NEWEX=EXCITGEN(i)
                DO j=0,iV-1
                    IF(i.EQ.j) THEN
                        !Don't think this is used
                        XIJ(i,i)=0
                    ELSE
                        CALL GenExcitProb(iGraph(1,i),iGraph(1,j),nEl,
     &                      NEWEX,G1,nBasisMax,UMat,Arr,nBasis,pGen)
                        XIJ(i,j)=pGen
                    ENDIF
                ENDDO
            ENDDO
            pGenGraph=GETPATHPROB2(XIJ,iV)
            
         ELSE
         
            DO i=0,iV-1
!  First work out which exictor we need to work out the excitation i->?
C            DO j=0,iV-1
C!  See if vertex j was excited from vertex i.  If so, then LOCTAB(1,j) is the excitor for vertex i
C               IF(LOCTAB(3,j).EQ.i) THEN
C                  IP_NEWEX=LOCTAB(1,j)
C                  EXIT
C               ENDIF
C            ENDDO
C            IF(j.EQ.iV-1) THEN
!We've not found an excitor, so we create one
                STORE(1)=0
C.. Setup the excit generator for the last vertex
                CALL GENSYMEXCITIT2(iGraph(1,i),NEL,G1,NBASIS,NBASISMAX,
     &              .TRUE.,NEWEXLEN,INODE2,IC,0,STORE,3)
                CALL MEMORY(IP_NEWEX,NEWEXLEN,'NEWEX')
                NEWEX(1)=0
C.. Count the excitations (and generate a random one which we throw)
                CALL GENSYMEXCITIT2(iGraph(1,i),NEL,G1,NBASIS,NBASISMAX,
     &              .TRUE.,NEWEX,INODE2,IC,0,STORE,3)
                CALL GENRANDSYMEXCITIT2(iGraph(1,i),NEL,G1,NBASIS,
     &               NBASISMAX,NEWEX,INODE2,ISEED,IC,0,UMAT,Arr,PGEN)
!
                DO j=0,iV-1
                    IF(i.EQ.j) THEN
                        XIJ(i,i)=IC
                    ELSE
                        CALL GenExcitProb(iGraph(1,i),iGraph(1,j),nEl,
     &                      NEWEX,G1,nBasisMax,UMat,Arr,nBasis,pGen)
                        XIJ(i,j)=pGen
                    ENDIF
                ENDDO
                CALL FREEM(IP_NEWEX)
            ENDDO
            pGenGraph=GETPATHPROB2(XIJ,iV)
         ENDIF
         
         IF(iUnit.NE.0) CALL WRITE_XMATRIX(iUnit,XIJ,iV)

      END


      SUBROUTINE DUMPIMEMORY(IUNIT,MEM,INLEN)
         IMPLICIT NONE
         INTEGER IUNIT,I,J
         INTEGER*1 MEM(*),CH
         INTEGER ILEN,IPOST,INLEN,LEN,ILEAD
         INTEGER ADDR,LEAD,POST,LINES
#ifdef POINTER8
         INTEGER*8 LOC
#else
         INTEGER LOC
#endif
         INTEGER,PARAMETER :: LINELEN=16
         ILEN=INLEN*4
         LEN=ILEN
         ADDR=16*(LOC(MEM)/LINELEN)
         ILEAD=LOC(MEM)-ADDR
         IPOST=MOD(LEAD+LEN,LINELEN)
         IF(IPOST.NE.0) IPOST=16-IPOST
         LINES=(ILEAD+ILEN+IPOST)/LINELEN
         DO I=0,LINES-1
            POST=IPOST
            LEN=ILEN
            LEAD=ILEAD
            WRITE(IUNIT,"(Z8.8)",advance='no') ADDR
            DO J=0,LINELEN-1
               IF(J.EQ.LINELEN/2) WRITE(IUNIT,'(A2)',advance='no') "  "
               IF(LEAD.EQ.0) THEN
                  IF(LEN.EQ.0) THEN
                     WRITE(IUNIT,"(A3)",advance='no') "   "
                  ELSE
                     WRITE(IUNIT,"(Z3.2)",advance='no')
     &                                              ,MEM(I*LINELEN+J+1)
                     LEN=LEN-1
                  ENDIF
               ELSE
                  WRITE(IUNIT,"(A3)",advance='no') "   "
                  LEAD=LEAD-1
               ENDIF 
            ENDDO
            WRITE(IUNIT,"(A3)",advance='no'),"   "
            DO J=0,LINELEN-1
               IF(J.EQ.LINELEN/2) WRITE(IUNIT,'(A1)',advance='no') " "
               IF(ILEAD.EQ.0) THEN
                  IF(ILEN.NE.0) THEN
                     CH=MEM(I*LINELEN+J+1)
                     IF(CH.GE.32) THEN
                        WRITE(IUNIT,"(A1)",advance='no'),CHAR(CH)
                     ELSE
                        WRITE(IUNIT,"(A1)",advance='no'),"."
                     ENDIF
                     ILEN=ILEN-1
                  ELSE
                     WRITE(IUNIT,"(A1)",advance='no')," "
                  ENDIF
               ELSE
                  WRITE(IUNIT,"(A1)",advance='no')," "
                  ILEAD=ILEAD-1
               ENDIF
            ENDDO
            ILEAD=LEAD
            IPOST=POST
            ILEN=LEN
            WRITE(IUNIT,*)
            ADDR=ADDR+LINELEN
         ENDDO
      END

