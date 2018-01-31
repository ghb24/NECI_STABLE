      ! A wrapper routine to allow symsetupexcit3_worker to be called
      ! with incorrect arguments --> Allows Alex's hacking!!
      SUBROUTINE SYMSETUPEXCITS3(NI,NEL,G1,NBASIS,STORE,SPIin,ETin,
     &   NAPin, OPin,
     &   TCOUNT,ICOUNT,CLASSES,ILUT,SYMPRODS,ILEVEL,iMinElec1,iMaxElec1)

          use systemdata, only: Basisfn, Symmetry
          use SymData, only: SymClass
          use SymExcit2, only: symsetupexcits3_worker

          integer :: nel, nI(nel), nBasis, store(6), iCount, iLevel
          integer :: ilut(0:nBasis/32), iMinElec1, iMaxElec1
          type(basisfn) G1(nBasis)
          logical :: tCount
          type(symclass) :: classes(*)
          type(symmetry) :: symprods(0:nel*nel)  !This nel*nel is the maximum it could possibly be

          integer :: SPIin(1:2,1:3,1:*), ETin(1:5,*), NAPin(1:3,*)
          integer :: OPin(1:2,*)

          call symsetupexcits3_worker (nI, nel, g1, nbasis, store,
     &             SPIin, ETin, NAPin, OPin, tCount, iCount, classes,
     &             ilut, symprods, iLevel, iMinElec1, iMaxElec1)

      end subroutine


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
         use SymData, only: tAbelianFastExcitgen
         use sym_mod, only: SYMPROD, TotSymRep, SYMCONJ, SYMEQ
         use sym_mod, only: FindSymLabel
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
         
         NCL=0  !This is our running count of how many symmetry classes we have
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

         IF(.not.tAbelianFastExcitgen) THEN
!This is not needed for abelian symmetry, since orb i and a must be the same symmetry class
C.. Now calculate the resultant symmetry from removing one of each of the classes (for singles)
C.. i.e. if we remove j, what's the symmetry of SYMCONJ(<ikl... |) (X) SYM( |ikl...>)
C..   NB to cope with non-abelians, this actually takes the prod sym of the classes of each orb, 
C..   rather than the orb itself, so may end up classifying a single as allowed when it's not.
C.. we'll use this in the future to determine the allowed j->a
C.. This is stored in CLASSES(iCl)%SymRem
             DO I=1,NCL
                PR=TotSymRep()
                DO J=1,NCL
                   L=CLASSCOUNT(1,I)+CLASSCOUNT(2,I)
                   IF(J.EQ.I) L=L-1
                   DO K=1,L
                      PR=SYMPROD( PR,SymConj(
     &                   SymLabels(CLASSES(J)%SymLab)))
                      PR=SYMPROD( PR,SymLabels(CLASSES(J)%SymLab))
                   ENDDO
                ENDDO
                CLASSES(I)%SymRem=PR
!            WRITE(6,"(A,I4,2Z4,I)") "CLASS ",NCL,CLASSES(I)%SymLab,
!     &            CLASSES(I)%SymRem,I
             ENDDO
!            WRITE(6,*) NCL," Symmetry Classes"
         ENDIF
      End Subroutine

C.. Now determine all possible symmetry products of pairs of orbs, and classify them.
C.. SYMPRODS(iPr) is the a symmetry product.
C.. SymProdCount(iType,iPr) is the number of prods of a given type.  iType=0 (beta/beta), 1 (alpha/beta), 2 (alpha/alpha)
C.. NPairs totals all SymProdCounts.

C.. For each pair of symmetry classes, determine its symmetry product.
C.. We list them and then sort them
      Subroutine SymSetupExcits_CreateCSProds(nPr,nPairs,nCl,
     &   SymProds, ThisClassCount, PrevClassCount, ClassCount,Classes,
     &   SymProdCount)
         use SystemData, only: Symmetry,SymmetrySize, nel
         use SymData, only: SymClass,SymClasses,SymLabels
         use sort_mod
         use sym_mod, only: SymProd, symeq
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
!K will then give the index for the unique symmetry product.
                  K=NPR
                  SYMPRODS(NPR)=PR  !Add the new symmetry product label to SymProds
               ENDIF
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
     
!Find total number of spin-orbital pairs of all spin-combos
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
         call sort (symProds(1:npr), symProdCount(:,1:npr))
!         WRITE(6,*) NPR," Symmetry Products"
!         WRITE(6,*) NPAIRS," Orbital Pairs"
C.. Create a cumulative sum.
         DO I=1,NPR
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)+SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)+SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)+SYMPRODCOUNT(3,I)
!            WRITE(6,*) SYMPRODS(I),SYMPRODCOUNT(1,I),
!     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
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
         use sym_mod, only: SymProd, BINARYSEARCHSYM, WRITESYM
         IMPLICIT NONE
         INTEGER nPairs
         INTEGER ORBPAIRS(2,nPairs)
         INTEGER nPr
         INTEGER iMinElec1, iMaxElec1,nEl,nI(nEl)
         TYPE(Symmetry) SymProds(0:*)
         INTEGER I,J,K,L
         Type(Symmetry) Pr
         INTEGER iSpn
         INTEGER iCC
         TYPE(BasisFN) G1(*)
         INTEGER SYMPRODIND(1:2,1:3,1:nPr)
         character(*), parameter :: t_r = 'SymSetupExcits_StoreOccPairs'

         ORBPAIRS(1:2,1:NPAIRS)=0
         L=0
         DO I=iMinElec1,iMaxElec1
            
            DO J=I+1,NEL
               L=L+1
               PR=SYMPROD(G1(NI(I))%Sym,G1(NI(J))%Sym)
C.. Find the product in the sorted symprod list
               CALL BINARYSEARCHSYM(PR,SYMPRODS(1),NPR,K)
               IF(K.EQ.0) THEN
                  WRITE(6,*) "Occupied Symmetry Products"
                  DO L=0,NPR
                     WRITE(6,"(I4)",advance='no') L
                     CALL WRITESYM(6,SYMPRODS(L),.TRUE.)
                  ENDDO
                  WRITE(6,"(A)",advance='no') "Illegal Symmetry"
                  CALL WRITESYM(6,PR,.TRUE.)
                  call stop_all(t_r, 'Illegal Symmetry Found')
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
!               WRITE(6,*) L,NI(I),NI(J),PR,K
            ENDDO
         ENDDO
      End Subroutine

C.. With each of the symmetry products, calculate []' such that []x[]'
C.. contains the totally symmetric rep.
C.. We do this by checking whether any of the []' in the 
C.. global symprods table multiplied by [] give the sym rep.

      Subroutine SymSetupExcits_CountVirtProds(nDoub, nExcitTypes,nPr,
     &   SymProdInd, SymProds, nAllowPPS,iLUT)
         use SystemData, only: Symmetry,BasisFN
         use SymData, only: SymPairProds,SymStatePairs,nSymPairProds
         use sym_mod, only: LSYMSYM, SYMPROD, SymConj
         IMPLICIT NONE
         INTEGER nDoub
         INTEGER nExcitTypes
         INTEGER nPr
         INTEGER nAllowPPS(3,*)
         INTEGER iLUT(0:*)
         INTEGER SYMPRODIND(2,3,1:*)
         TYPE(Symmetry) SYMPRODS(0:*)

         INTEGER I,J,K,L
         Type(Symmetry) SPP
         INTEGER ICC1,ICC2,ICC3,ICC4
         LOGICAL L1B,L1A,L2B,L2A
         INTEGER ICC,iSpn
   
C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements


C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) is a global list of all SYMPRODs available, the number of pairs of
C.. symlabels (listed in SymStatePairs), and the index of the start of this list
C.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
C.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of symlabels whose prod is of that symmetry.

         NDOUB=0
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
                     L=NALLOWPPS(ISPN,J)
                     K=SYMPRODIND(2,ISPN,I)
                     ICC=L*K
!                     WRITE(6,*) "ET",NEXCITTYPES
!                     WRITE(6,"(14I5)") I,J,SYMPRODS(I),
!     &                  SYMPAIRPRODS(J),ISPN,L,K,ICC
                     NDOUB=NDOUB+ICC
                     IF(ICC.GT.0) NEXCITTYPES=NEXCITTYPES+1
                  ENDDO
               ENDIF
            ENDDO
         ENDDO       

!         WRITE(6,*) "Number of double excitations: ",NDOUB
      End Subroutine

      Subroutine SymSetupExcitsAb_CountVProds(nDoub,nExcitTypes,
     &            nCl,nPr,SymProds,SymProdInd,Classes,ClassCount,
     &            nAllowPPS) ! Calculated not enumerated
         use SymData, only: SymPairProds,nSymPairProds,SymmetrySize
         use SymData, only: Symmetry,SymLabelCounts,SymLabels,SymClass
         use sym_mod, only: SymConj, SymProd, FindSymLabel
         use sym_mod, only: BINARYSEARCHSYM
         IMPLICIT NONE
         Type(Symmetry) SymProds(0:*)
         INTEGER nPr,nCl
         INTEGER ClassCount(2,*)
         Type(SymClass) Classes(nCl)
         INTEGER nAllowPPS(1:3,nSymPairProds)
C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements
         INTEGER SymProdInd(1:2,1:3,1:nPr)
         INTEGER nDoub
         INTEGER nExcitTypes
         INTEGER I,L,K,KK,ICC,ISPN,lab,iOccs,iVirts,iClass,inCl
         Type(Symmetry) s
         Logical tDebugPrint
         tDebugPrint=.false.
         nDoub=0 
         nAllowPPS=0
! For each Symmetry product count the number of allowed doubles.
! Loop over the complete list of pair types
!  Each pair prod in the full list can only interact with prods of the same sym in the list of occ pairs.
         DO I=1,nSymPairProds
! Now search for that symmetry within our occupied pair counts
           CALL BINARYSEARCHSYM(SymPairProds(I)%Sym,SYMPRODS(1),
     &         NPR,K)
!K now has the index of that sym in the occ sym prod list.

            IF(tDebugPrint)
     &     WRITE(6,*) "Sym:",SymPairProds(I)%Sym
           IF(K.NE.0) THEN
!.. Now count all the excitations allowed for this set of []->[]'
               DO ISPN=1,3
! First get the virts
                  If(iSpn.eq.2) THEN
                     L=SymPairProds(I)%nPairsStateOS
                  else
                     L=SymPairProds(I)%nPairsStateSS
                  ENDIF
! L now contains total number (including occ and single excits) of double excits of this sym.
! We must remove the single excitations.
! The number of pairs of doubles in the total list with sym prod SymProds(K)
!  which contain just one of the occupied orbitals can be worked out by
!  multiplying SymProds(K) by the inverse of the symmetry of each occupied class.
!   This will give the sym of the virtual which would make up the pair.  
!   Since we know how many occ there are in each occ class, and can calculated
!   How many virts there are in its complement, that gives us the number of 
!  'single excitations' masquerading as doubles.

!  First loop over occupied classes
                  DO iClass=1,nCl
!  SymLabels(Classes(iClass)%SymLab) is the sym of the class
!  WARNING - here we use SymConj to mean 'inverse'.
!   This seems right, although it would be nice to have it really proven.
                     if(tDebugPrint)
     &                  WRITE(6,*) "Occ Class,sym:",iClass,
     &                  SymLabels(Classes(iClass)%SymLab)
                     s=SymProd(
     &                   SymConj(SymLabels(Classes(iClass)%SymLab)),
     &                   SymPairProds(I)%Sym)
                     lab=FindSymLabel(s)
!   SymLabelCounts(2,lab) is the total number of states with the inverse sym.
                     iVirts=SymLabelCounts(2,lab)
                     if(tDebugPrint)
     &                  Write(6,*) "Virts,lab",iVirts,lab
                     do iNCl=1,nCl
                        if(Classes(iNCl)%SymLab.eq.lab) exit
                     enddo
                     if(tDebugPrint)
     &                  WRITE(6,*) "Virt Class,sym:",iNCl,s
                     if(iSpn.eq.2) then
!The mixed spin case
                        if(tDebugPrint)
     &                  WRITE(6,*) "OS Tot:",L
                        if(iNCL.le.nCl) then
!  We've found the class corresponding to the inverse, so we know how many occs it has
                           iOccs=ClassCount(1,iClass)
                           L=L-iOccs*(iVirts-ClassCount(2,iNCL))
                           if(tDebugPrint)
     &                        WRITE(6,*) "Aocc Bvirt removed:", 
     &         iOccs*(iVirts-ClassCount(2,iNCL)),ClassCount(2,iNCL)
                           iOccs=ClassCount(2,iClass)
                           L=L-iOccs*(iVirts-ClassCount(1,iNCL))
                           if(tDebugPrint)
     &                        WRITE(6,*) "Bocc Avirt removed:", 
     &         iOccs*(iVirts-ClassCount(1,iNCL)),ClassCount(1,iNCL)
                        else
                           iOccs=ClassCount(1,iClass)
                           L=L-iOccs*iVirts
                           if(tDebugPrint)
     &                     WRITE(6,*) "Aocc Bvirt removed:",iOccs*iVirts
                           iOccs=ClassCount(2,iClass)
                           L=L-iOccs*iVirts
                           if(tDebugPrint)
     &                     WRITE(6,*) "Bocc Avirt removed:",iOccs*iVirts
                        endif
                     else
! same spin case
                        iOccs=ClassCount((iSpn+1)/2,iClass)
                        if(tDebugPrint)
     &                  WRITE(6,*) "SS Tot:",L
                        if(iNCL.le.nCl) then
                          L=L-iOccs*(iVirts-ClassCount((iSpn+1)/2,iNCL))
                        else
                          L=L-iOccs*iVirts
                        endif
                        if(tDebugPrint)
     &                  WRITE(6,*) "SS singles removed:",L
                     endif
                  enddo
! The number of occupied pairs with sym prod SymProds(K)
                  KK=SymProdInd(2,iSpn,K)
                  if(tDebugPrint)
     &                  WRITE(6,*) "Occ pairs", KK
! Also remove the occupied pairs from the count of virtuals
                  L=L-KK
! Save this for later
                  nAllowPPS(iSpn,I)=L
                  ICC=L*KK
                  if(tDebugPrint)
     &                  WRITE(6,*) "Occ removed", L
                  if(tDebugPrint) WRITE(6,*) "ET",NEXCITTYPES
!                  if(tDebugPrint)
!     &                  WRITE(6,"(14I5)") K,I,
!     &                  SYMPAIRPRODS(I),ISPN,L,KK,ICC
                  NDOUB=NDOUB+ICC
                  IF(ICC.GT.0) NEXCITTYPES=NEXCITTYPES+1
               ENDDO
             
           ENDIF 
         ENDDO                
      End Subroutine

      Subroutine SymSetupExcits_CountSingles(nSing,nCl,nExcitTypes,
     &   ThisClassCount, ClassCount,Classes)
         use SystemData, only: Symmetry
         use SymData, only: SymClass,SymLabelCounts,nSymLabels,SymLabels
         use sym_mod, only: LSYMSYM, SYMPROD, SymConj
         IMPLICIT NONE
         INTEGER nSing
         INTEGER nCl
         INTEGER nExcitTypes
         INTEGER ThisClassCount(2,*)
         INTEGER ClassCount(2,*)
         Type(SymClass) Classes(nCl)
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

!A simple routine for counting the number of single excitations in abelian systems.
!orbital 'a' must have the same symmetry and spin as orbital 'b'
!Simply run through symmetry classes in NI to count these.
!NExcitTypes counts the number of different ways of getting an excitation, i.e. Different spin or symmetry
      Subroutine SymSetupExcitsAb_CountSing(nSing,nCl,
     &      nExcitTypes,ThisClassCount,Classes)
         use SystemData, only: Symmetry
         use SymData, only: SymClass,SymLabelCounts,nSymLabels,SymLabels
         IMPLICIT NONE
         INTEGER nSing
         INTEGER nCl
         INTEGER nExcitTypes
         INTEGER ThisClassCount(2,*)
         Type(SymClass) Classes(nCl)
         INTEGER I,ICC
         nSing=0
!Loop over all symmetry classes of your electrons
         DO I=1,nCl

!First deal with alpha spins - there are ThisClassCount(1,I) alpha spins in this determinant.
!We need to know how many unoccupied alpha spins with the same symmetry there are.
!There are SymLabelCounts(1,Classes(I)%SymLab) total alpha spin-orbitals with the same symmetry.
!SymLabelCounts(2,Classes(I)%SymLab)-ThisClassCount(1,I) indicates the number of unoccupied alpha orbitals with the required symmetry.
!Each electron can excite to each unoccupied orbital, giving us the number of single excitations for this symmetry and spin

             IF(ThisClassCount(1,I).ne.0) THEN
!Aren't going to be any singles if there are no electrons of that spin in the symmetry...
                ICC=(SymLabelCounts(2,Classes(I)%SymLab)-
     &                      ThisClassCount(1,I))*ThisClassCount(1,I)
                IF(ICC.ne.0) THEN
!We have found some alpha->alpha single excitations
                    nSing=nSing+ICC
                    nExcitTypes=nExcitTypes+1
                ENDIF

             ENDIF

!Now do the same for beta spins
             IF(ThisClassCount(2,I).ne.0) THEN
                ICC=(SymLabelCounts(2,Classes(I)%SymLab)-
     &                      ThisClassCount(2,I))*ThisClassCount(2,I)
                IF(ICC.ne.0) THEN
!We have found some single excitations
                    nSing=nSing+ICC
                    nExcitTypes=nExcitTypes+1
                ENDIF
             ENDIF

         ENDDO

      End Subroutine SymSetupExcitsAb_CountSing

!Based on the Count singles analogue, this fills the ExcitTypes array with
!symmetry and excitation type information.
      Subroutine SymSetupExcitsAb_StoreSing(nExcitTypes,nCl,
     &      Classes,ThisClassCount,ExcitTypes)
         use SystemData, only: Symmetry
         use SymData, only: SymClass,SymLabelCounts,nSymLabels,SymLabels
         IMPLICIT NONE
         INTEGER nCl
         INTEGER nExcitTypes
         INTEGER ThisClassCount(2,*)
         Type(SymClass) Classes(nCl)
         INTEGER I,ICC
         INTEGER ExcitTypes(5,*)
!nExcitTypes is reset to zero before this routine - the nExcitTypes will be incremented for doubles later.
!Loop over all symmetry classes of your electrons
         DO I=1,nCl

!First deal with alpha spins - there are ThisClassCount(1,I) alpha spins in this determinant.
!We need to know how many unoccupied alpha spins with the same symmetry there are.
!There are SymLabelCounts(1,Classes(I)%SymLab) total alpha spin-orbitals with the same symmetry.
!SymLabelCounts(2,Classes(I)%SymLab)-ThisClassCount(1,I) indicates the number of unoccupied alpha orbitals with the required symmetry.
!Each electron can excite to each unoccupied orbital, giving us the number of single excitations for this symmetry and spin
             IF(ThisClassCount(1,I).ne.0) THEN
!Aren't going to be any singles if there are no electrons of that spin in the symmetry...
                ICC=(SymLabelCounts(2,Classes(I)%SymLab)-
     &                      ThisClassCount(1,I))*ThisClassCount(1,I)
                IF(ICC.ne.0) THEN
!We have found some alpha->alpha single excitations of a new 'type'
                    nExcitTypes=nExcitTypes+1
!Put information on these type of excitation into the excittypes array
                    ExcitTypes(1,nExcitTypes)=1     !Its a single
                    ExcitTypes(2,nExcitTypes)=1     !Its an alpha->alpha excitation
                    ExcitTypes(3,nExcitTypes)=I     !Symmetry class of occupied orbital
                    ExcitTypes(4,nExcitTypes)=Classes(I)%SymLab     
                                                    !Symmetry class of virtual orbital (the same!)
                                                    !However, it wants to be the symmetry class from all
                                                    !symmetry classes
                    ExcitTypes(5,nExcitTypes)=ICC   !Number of this type of excitation
                ENDIF

             ENDIF

!Now do the same for beta spins
             IF(ThisClassCount(2,I).ne.0) THEN
                ICC=(SymLabelCounts(2,Classes(I)%SymLab)-
     &                      ThisClassCount(2,I))*ThisClassCount(2,I)
                IF(ICC.ne.0) THEN
!We have found some beta single excitations of a new 'type'
                    nExcitTypes=nExcitTypes+1
!Put information on these type of excitation into the excittypes array
                    ExcitTypes(1,nExcitTypes)=1     !Its a single
                    ExcitTypes(2,nExcitTypes)=2     !Its an beta->beta excitation
                    ExcitTypes(3,nExcitTypes)=I     !Symmetry class of occupied orbital
                    ExcitTypes(4,nExcitTypes)=Classes(I)%SymLab     
                                                    !Symmetry class of virtual orbital (the same!)
                                                    !However, it wants to be the symmetry class from all
                                                    !symmetry classes
                    ExcitTypes(5,nExcitTypes)=ICC   !Number of this type of excitation
                ENDIF
             ENDIF

         ENDDO

      End Subroutine SymSetupExcitsAb_StoreSing

      
      Subroutine SymSetupExcits_StoreSingles(nExcitTypes,nCl,
     &   Classes,ThisClassCount,ExcitTypes)
         Use SystemData, only: Symmetry
         use SymData, only: SymClass,SymPairProds,nSymLabels,SymLabels
         use SymData, only: SymLabelCounts
         use sym_mod, only: SYMPROD, SymConj, LSYMSYM
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
         use sym_mod, only: LSYMSYM, SYMPROD, SymConj
         IMPLICIT NONE
         INTEGER nPr
         INTEGER nSymPairProds
         INTEGER nAllowPPS(3,*)
         INTEGER ExcitTypes(5,*)
         INTEGER nExcitTypes
         TYPE(Symmetry) SYMPRODS(0:*)
         INTEGER SYMPRODIND(2,3,1:*)
         
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
!                        WRITE(6,"(A,3I5,4I5)") "SD",NEXCITTYPES,I,J,
!     &                     ISPN,L,K,ICC
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      End Subroutine


      SUBROUTINE SYMGENALLEXCITS(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
     &               SYMPRODIND,ILUT,ORBPAIRS,LSTE,ICLIST,ICOUNT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         use SymData, only: tAbelianFastExcitGen
         IMPLICIT NONE
         INTEGER NEXCITTYPES
         INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES)
         INTEGER I,J,K,L
         INTEGER NK(NEL)
         TYPE(SymClass) CLASSES(*)
         INTEGER IEXCIT
         INTEGER IFROM,ITO,ISPN
         INTEGER IFROMSL,ITOSL
         INTEGER ExcitMat(2,2)
         TYPE(Symmetry) SPP
         LOGICAL L1,L2
         INTEGER ICC1,ICC2,ICC3,ICC4
         LOGICAL L1B,L1A,L2B,L2A,TParity

         INTEGER ORBPAIRS(2,*)
         INTEGER ILUT(0:*)
         INTEGER SYMPRODIND(2,3,1:*)
         INTEGER LSTE(NEL,*),ICLIST(*),ICOUNT
         character(66) error_msg
         character(*), parameter :: this_routine = 'SYMGENALLEXCITS'

         ICOUNT=0
         if(tAbelianFastExcitGen) then
             error_msg='Fast Abelian excitation generators do not work '
     &       //'for SymGenAllExcits'
             call stop_all(this_routine, "Wrong excitation generator")
         end if

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
!                           CALL NECI_ICOPY(NEL,NI,1,NK,1)
                           NK(:)=NI(:)
!Find the ordered excited determinant
                           ExcitMat(1,1)=I  !This is the position of the orbital to excite
                           ExcitMat(2,1)=ITOSL
                           CALL FindExcitDet(ExcitMat,NK,1,TParity)
!                           NK(I)=ITOSL
!                           CALL NECI_SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=1
!                           CALL NECI_ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
                           LSTE(1:NEL,ICOUNT)=NK(:)
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
!                           CALL NECI_ICOPY(NEL,NI,1,NK,1)
                           NK(1:NEL)=NI(1:NEL)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,1)=J
                                 ExcitMat(2,1)=ICC1
!                                 NK(J)=ICC1
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(2,1)=J
                                 ExcitMat(2,2)=ICC3
!                                 NK(J)=ICC3
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL FindExcitDet(ExcitMat,NK,2,TParity)
!                           CALL NECI_SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
!                           CALL NECI_ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
                           LSTE(1:NEL,ICOUNT)=NK(1:NEL)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ELSEIF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
                        IF(.NOT.(L1B.OR.L2A)) THEN
!                           CALL NECI_ICOPY(NEL,NI,1,NK,1)
                           NK(1:NEL)=NI(1:NEL)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,1)=J
                                 ExcitMat(2,1)=ICC1
!                                 NK(J)=ICC1
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,2)=J
                                 ExcitMat(2,2)=ICC4
!                                 NK(J)=ICC4
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL FindExcitDet(ExcitMat,NK,2,TParity)
!                           CALL NECI_SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
!                           CALL NECI_ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
                           LSTE(1:NEL,ICOUNT)=NK(1:NEL)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
                        IF(.NOT.(L1A.OR.L2B).AND.ICC1.NE.ICC3) THEN
!                           CALL NECI_ICOPY(NEL,NI,1,NK,1)
                           NK(1:NEL)=NI(1:NEL)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,1)=J
                                 ExcitMat(2,1)=ICC2
!                                 NK(J)=ICC2
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,2)=J
                                 ExcitMat(2,2)=ICC3
!                                 NK(J)=ICC3
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
!                           CALL NECI_SORTI(NEL,NK)
                           CALL FindExcitDet(ExcitMat,NK,2,TParity)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
!                           CALL NECI_ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
                           LSTE(1:NEL,ICOUNT)=NK(1:NEL)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ELSEIF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                        IF(ICC1.NE.ICC3.AND..NOT.(L1A.OR.L2A)) THEN
!                           CALL NECI_ICOPY(NEL,NI,1,NK,1)
                           NK(1:NEL)=NI(1:NEL)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,1)=J
                                 ExcitMat(2,1)=ICC2
!                                 NK(J)=ICC2
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 ExcitMat(1,2)=J
                                 ExcitMat(2,2)=ICC4
!                                 NK(J)=ICC4
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
!                           CALL NECI_SORTI(NEL,NK)
                           CALL FindExcitDet(ExcitMat,NK,2,TParity)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
!                           CALL NECI_ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
                           LSTE(1:NEL,ICOUNT)=NK(1:NEL)
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
         use SystemData , only : tAssumeSizeExcitgen
         IMPLICIT NONE
         INTEGER NMEM(*)
         INTEGER iCount
         IF(tAssumeSizeExcitgen) THEN
             iCount=NMEM(1)
         ELSE
            iCount=NMEM(23)
         ENDIF
         RETURN
      END SUBROUTINE GetSymExcitCount


!GenSymExcitIt3Par and GenSymExcitIt3 are versions GenSymExcitIt2 which do not require information taken from System and Integrals modules to be specified on the command line
      SUBROUTINE GENSYMEXCITIT3(NI,TSETUP,NMEM,NJ,IC,STORE,ILEVEL)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         use legacy_data, only: irat
         IMPLICIT NONE 
         INTEGER NI(NEL)
         INTEGER, pointer :: DSTORE(:)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER NMEM(*),NJ(NEL),IC
         LOGICAL TSETUP
         INTEGER ILEVEL

         CALL GENSYMEXCITIT2Par(NI,NEL,G1,NBASIS,
     &         TSETUP,NMEM,NJ,IC,STORE,ILEVEL,1,nEl)
      END

      SUBROUTINE GENSYMEXCITIT3Par(NI,
     &         TSETUP,NMEM,NJ,IC,STORE,ILEVEL,iMinElec1,iMaxElec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         use legacy_data, only: irat
         IMPLICIT NONE 
         INTEGER NI(NEL)
         INTEGER, pointer :: DSTORE(:)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER, target ::  NMEM(*)
         INTEGER NJ(NEL),IC
         LOGICAL TSETUP
         INTEGER ILEVEL

         INTEGER iMinElec1, iMaxElec1
         CALL GENSYMEXCITIT2Par(NI,NEL,G1,NBASIS,
     &         TSETUP,NMEM,NJ,IC,STORE,ILEVEL,iMinElec1,iMaxElec1)
      END

C.. A modified version of GENSYMEXCITIT2.  This has min and max electron 1 specifiers used for parallelization
      SUBROUTINE GENSYMEXCITIT2Par(NI,NEL,G1,NBASIS,
     &         TSETUP,NMEM,NJ,IC,STORE,ILEVEL,iMinElec1,iMaxElec1)

         use SystemData, only: basisfn
         use SymExcit2, only: gensymexcitit2par_worker

         integer :: nel, nI(nel), nbasis, nmem(*), nJ(nel), ic, store(6)
         integer :: ilevel, iMinElec1, iMaxElec1
         type(basisfn) :: g1(nbasis)
         logical :: tsetup
         
         call gensymexcitit2par_worker (nI, nel, g1, nbasis, tSetup,
     &                      nmem, nJ, ic, store, ilevel, iminelec1,
     &                      imaxelec1)


      end subroutine

      SUBROUTINE RESETEXIT2(NMEM)
         IMPLICIT NONE
         INTEGER NMEM(*)
!         NMEM(1:22)=0
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
      INTEGER NEXCITTYPES
      INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES)
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
      LOGICAL LS(2,2),TParity

      INTEGER ORBPAIRS(2,*)
      INTEGER ILUT(0:*)
      INTEGER SYMPRODIND(2,3,1:*)
      INTEGER IC,ExcitMat(2,2)
      LOGICAL tDebugPrint
      INTEGER iMinElec1, iMaxElec1
      character(*), parameter :: this_routine = 'SYMGENEXCITIT'
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
               call stop_all(this_routine,
     &                'tStoreAsExcitations not tested for singles.')
            ELSE
!               CALL NECI_ICOPY(NEL,NI,1,NK,1)
               NK(1:NEL)=NI(1:NEL)
               IF(tDebugPrint)
     &        WRITE(6,*) "[",NK(I),"->",ITOSL,"]"
               ExcitMat(1,1)=I
               ExcitMat(2,1)=ITOSL
               CALL FindExcitDet(ExcitMat,NK,1,TParity)
!               NK(I)=ITOSL
!               CALL NECI_SORTI(NEL,NK)
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
!                     CALL NECI_ICOPY(NEL,NI,1,NK,1)
                     NK(1:NEL)=NI(1:NEL)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(1),","
                           ExcitMat(1,1)=J
                           ExcitMat(2,1)=ICC(1)
!                           NK(J)=ICC(1)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(3),"]"
!                           NK(J)=ICC(3)
                           ExcitMat(1,2)=J
                           ExcitMat(2,2)=ICC(3)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
!                     CALL NECI_SORTI(NEL,NK)
                     CALL FindExcitDet(ExcitMat,NK,2,TParity)
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
!                     CALL NECI_ICOPY(NEL,NI,1,NK,1)
                     NK(1:NEL)=NI(1:NEL)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(1),","
                           ExcitMat(1,1)=J
                           ExcitMat(2,1)=ICC(1)
!                           NK(J)=ICC(1)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(4),"]"
                           ExcitMat(1,2)=J
                           ExcitMat(2,2)=ICC(4)
!                           NK(J)=ICC(4)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
                     CALL FindExcitDet(ExcitMat,NK,2,TParity)
!                     CALL NECI_SORTI(NEL,NK)
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
!                     CALL NECI_ICOPY(NEL,NI,1,NK,1)
                     NK(1:NEL)=NI(1:NEL)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(2),","
                           ExcitMat(1,1)=J
                           ExcitMat(2,1)=ICC(2)
!                           NK(J)=ICC(2)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(3),"]"
!                           NK(J)=ICC(3)
                           ExcitMat(1,2)=J
                           ExcitMat(2,2)=ICC(3)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
!                     CALL NECI_SORTI(NEL,NK)
                     CALL FindExcitDet(ExcitMat,NK,2,TParity)
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
!                     CALL NECI_ICOPY(NEL,NI,1,NK,1)
                     NK(1:NEL)=NI(1:NEL)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,"(A,I3,A,I3,A)",advance='no')
     &                                        "[",NK(J),"->",ICC(2),","
                           ExcitMat(1,1)=J
                           ExcitMat(2,1)=ICC(2)
!                           NK(J)=ICC(2)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           IF(tDebugPrint)
     &                       WRITE(6,*) NK(J),"->",ICC(4),"]"
!                           NK(J)=ICC(4)
                           ExcitMat(1,2)=J
                           ExcitMat(2,2)=ICC(4)
                        ENDIF
                     ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
!                     CALL NECI_SORTI(NEL,NK)
                     CALL FindExcitDet(ExcitMat,NK,2,TParity)
                   ENDIF
                   IC=2
                   EXIT
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO 
      END

      Subroutine SymGenExcitIt2(nI,nEl,ExcitTypes,nExcitTypes,
     &   Classes, SymProdInd, iLUT, OrbPairs, iExcit, iSpn,
     &   iFrom, iTo, I, J, K, L, iCC, LS,
     &   nK, iC, iMinElec1, iMaxElec1)
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER NEXCITTYPES
         INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES)
         INTEGER I,J,K,L
         INTEGER NK(NEL)
         TYPE(SymClass) CLASSES(*)
         INTEGER IEXCIT
         INTEGER IFROM,ITO,ISPN
         TYPE(Symmetry) SPP
         INTEGER ICC(4)
C.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.
         LOGICAL LS(2,2)

         INTEGER ORBPAIRS(2,*)
         INTEGER ILUT(0:*)
         INTEGER SYMPRODIND(2,3,1:*)
         INTEGER IC
         LOGICAL tDebugPrint
         INTEGER iMinElec1, iMaxElec1

         INTEGER iRet
         tDebugPrint=.false.
         if(tDebugPrint) write(6,*) "Entering SymGenExcitIt2"

         DO WHILE(.TRUE.)
! Indicate that we need to keep on cycling.  If it is changed
!  by the later routines we stop the loop
            iRet=0
C.. see if we need a new EXCIT
            if(tDebugPrint) WRITE(6,*) "I,iExcit,nExcitTypes",
     &         I,iExcit,nExcitTypes
            IF(I.LT.0) THEN
C.. move to the next excitation
               IEXCIT=IEXCIT+1
               IF(IEXCIT.GT.NEXCITTYPES) THEN
C.. We're done
                  NK(1)=0
                  RETURN
               ENDIF
               IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
                  Call SymGenExcitIt_SetupSingle(iSpn,iFrom,iTo,
     &      iExcit, ExcitTypes, Classes, SymProdInd,
     &      I,J,K,L,nI,nEl,iMinElec1,tDebugPrint)
               ELSE
                  Call SymGenExcitIt_SetupDouble(iSpn,iFrom,iTo,
     &      iExcit, ExcitTypes, SymProdInd,
     &      I,K,L,tDebugPrint)
               ENDIF 
            ENDIF
            IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
               Call SymGenExcitIt_GenSingle(iRet,iSpn,iFrom,iTo,iCC,
     &   I,J,K,L,tDebugPrint,iMaxElec1,iC,
     &   nI,nK,nEl)
            ELSE
               Call SymGenExcitIt_GenDouble(iRet,SPP,
     &  I,K,L,iTo,iSpn,iFrom,iCC,iLUT,LS,OrbPairs,SymProdInd,
     &   nI,nK,nEl,iC,tDebugPrint)
            ENDIF 
            if(iRet.eq.1) cycle
            if(iRet.eq.2) return
         ENDDO 
      End Subroutine SymGenExcitIt2

      Subroutine SymGenExcitIt_SetupSingle(iSpn,iFrom,iTo,
     &      iExcit, ExcitTypes, Classes, SymProdInd,
     &      I,J,K,L,nI,nEl,iMinElec1,tDebugPrint)
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER nEl,nI(nEl)
         INTEGER ExcitTypes(5,*)
         
         INTEGER I,J,K,L
         TYPE(SymClass) Classes(*)
         INTEGER iExcit
         INTEGER iFrom,iTo,iSpn

         INTEGER SymProdInd(2,3,1:*)
         LOGICAL tDebugPrint
         INTEGER iMinElec1


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
     &      WRITE(6,*) "EXC",IEXCIT,ISPN,IFROM,ITO,SYMPRODIND(2,ISPN,
     &         IFROM)
      End Subroutine SymGenExcitIt_SetupSingle

      Subroutine SymGenExcitIt_SetupDouble(iSpn,iFrom,iTo,
     &      iExcit, ExcitTypes, SymProdInd,
     &      I,K,L,tDebugPrint)
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER ExcitTypes(5,*)
         
         INTEGER I,K,L
         INTEGER iExcit
         INTEGER iFrom,iTo,iSpn

         INTEGER SymProdInd(2,3,1:*)
         LOGICAL tDebugPrint

C.. a double
         ISPN=EXCITTYPES(2,IEXCIT)
         IFROM=EXCITTYPES(3,IEXCIT)
         ITO=EXCITTYPES(4,IEXCIT)
         I=0
         L=5 !An invalid spin
         K=-2 !When we +1 this is still invalid so it resets the pair index
         IF(tDebugPrint)
     &      WRITE(6,*) "EXC",IEXCIT,ISPN,IFROM,ITO,SYMPRODIND(2,ISPN,
     &            IFROM)
      End Subroutine SymGenExcitIt_SetupDouble

      Subroutine SymGenExcitIt_GenSingle(iRet,iSpn,iFrom,iTo,iCC,
     &   I,J,K,L,tDebugPrint,iMaxElec1,iC,
     &   nI,nK,nEl)
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER nEl,nI(nEl)
         
         INTEGER I,J,K,L
         INTEGER nK(nEl)
         INTEGER iFrom,iTo,iSpn
         INTEGER ICC(4)


         INTEGER IFROMSL,ITOSL
         LOGICAL L1,L2,TParity
C.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.

         INTEGER iC,ExcitMat(2,2)
         LOGICAL tDebugPrint
         INTEGER iMaxElec1

         INTEGER iRet
         character(*),parameter:: this_routine='SymGenExcitIt_GenSingle'
         
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
            iRet=1
            RETURN               
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
            L2 = .false.
            if (L <= nel) then
                if (nI(L) < iToSL) L2 = .true.
            endif
         ENDDO
         IF(L.LE.NEL) then
             if (NI(L).EQ.ITOSL) THEN
C.. We've found an L in the det which is the same as ITOSL. 
C.. we go round again, getting another K
                iRet=1
                RETURN               
            endif
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
            call stop_all(this_routine,
     &              'tStoreAsExcitations not tested for singles.')
         ELSE
!            CALL NECI_ICOPY(NEL,NI,1,NK,1)
            NK(1:NEL)=NI(1:NEL)
            IF(tDebugPrint)
     &        WRITE(6,*) "[",NK(I),"->",ITOSL,"]"
            ExcitMat(1,1)=I
            ExcitMat(2,1)=ITOSL
!            NK(I)=ITOSL
!            CALL NECI_SORTI(NEL,NK)
            CALL FindExcitDet(ExcitMat,NK,1,TParity)
            IC=1
         ENDIF
C.. quit the do loop
         iRet=2
         RETURN
      End Subroutine SymGenExcitIt_GenSingle

!Because the fast abelian excitation generators don't keep a stored list of virtual pairs, they need a different excitation generator routine.  We generate the list on the fly.
      Subroutine SymGenExcitItOld_GenDouble(iRet,SPP,
     &   I,K,L,iTo,iSpn,iFrom,iCC,iLUT,LS,OrbPairs,SymProdInd,
     &   nI,nK,nEl,iC,tDebugPrint)
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         
         INTEGER I,K,L
         INTEGER NK(NEL)
         INTEGER IFROM,ITO,ISPN
         TYPE(Symmetry) SPP
         INTEGER ICC(4)
C.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.
         LOGICAL LS(2,2)

         INTEGER ORBPAIRS(2,*)
         INTEGER ILUT(0:*)
         INTEGER SYMPRODIND(2,3,1:*)
         INTEGER IC
         LOGICAL tDebugPrint

         INTEGER iRet
         iRet=0
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
                  iRet=1
                  RETURN
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
               Call SymGenExcitIt_MakeDouble(
     &            ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &            ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &            ICC(1),ICC(3),nI,nK,nEl,tDebugPrint)
               IC=2
               iRet=2
               RETURN
            ENDIF
         ENDIF
         IF(L.EQ.2) THEN
            L=3
            IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
               IF(.NOT.(LS(1,1).OR.LS(2,2))) THEN
                  Call SymGenExcitIt_MakeDouble(
     &             ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ICC(1),ICC(4),nI,nK,nEl,tDebugPrint)
                  IC=2
                  iRet=2
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         IF(L.EQ.3) THEN
            L=4
            IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
              IF(.NOT.(LS(1,2).OR.LS(2,1)).AND.ICC(1).NE.ICC(3)) THEN
                Call SymGenExcitIt_MakeDouble(
     &             ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ICC(2),ICC(3),nI,nK,nEl,tDebugPrint)
                IC=2
                iRet=2
                RETURN
               ENDIF
            ENDIF
         ENDIF
         IF(L.EQ.4) THEN
            L=5
            IF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the same and neither are in NI, then allow
              IF(ICC(1).NE.ICC(3).AND..NOT.(LS(1,2).OR.LS(2,2))) THEN
                Call SymGenExcitIt_MakeDouble(
     &             ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ICC(2),ICC(4),nI,nK,nEl,tDebugPrint)
                IC=2
                iRet=2
                RETURN
               ENDIF
            ENDIF
         ENDIF
      End Subroutine SymGenExcitItOld_GenDouble

!Get a new pair from the global list (if there is one) or generate one afresh.  All according to the pair type iTo which is the index into SymPairProds.
!   iLooped will be set to 1 by this if we've run out of pairs and had to go back to the beginning.
!   The pair is retured in iTo1 and iTo2.
!  K is a counter for us to keep track of where we are.  If K<-1 on entry, we
!  start afresh and set iLooped to -1.
      Subroutine SymGenExcitIt_GetNextPair(K,iTo,iLooped,iTo1,iTo2,
     &      tDebugPrint)
         use SymData, only: tStoreStateList
         use SymData, only: SymLabelCounts,SymStatePairs
         use SymData, only: SymLabelList,SymPairProds
         use SymData, only: SymClasses
         IMPLICIT NONE
         INTEGER K,iTo
         INTEGER iLooped
         INTEGER iTo1,iTo2

         INTEGER iLabelPairIndexOff
         INTEGER iLabelPairBase
         INTEGER iLabel1Index,iLabel2Index
         INTEGER iLabel1,iLabel2
         LOGICAL tDebugPrint
         iLooped=0
         if(.not.tStoreStateList) then
! We must generate the pair as we don't have a stored list
! SymStatePairs holds a list of SymLabels.
! SymPairProds(iTo) tells is which set of pairs we're using
! iTo is the symmetry of the 'to' spin-orb pair
! we split K into 3 integers as we need to store 3 variables rather than 1
!   (LabelPairIndexOff, Label1Index, Label2Index)
!bits      30-24         23-12           11-0
!The index of the current state of the 1st symmetry of the pair (Label1Index)             
!The index of the current state of the 2nd symmetry of the pair (Label2Index)

!we loop over these with Label2Index changing fastest, then Label1Index, then iLabelPairIndexOff
! LabelPairIndexOff:  Which pair of symmetries (in SymStatePairs)
!                       This is the index for a given iTo.  The base is iLabelPairBase
! Label1Index:        The index of the current first state within its symmetry class
! Label2Index:        The index of the current second state within it symmetry class

!  the symmetries of the states are determined by the LAbelPairIndexOff.
!   We iterate through all possible [LabelPair,(state1,state2)], in order state2, state1, LabelPair.


            if(K.lt.0) then
               iLabelPairIndexOff=0
               iLabel1Index=0
               iLabel2Index=-1
               iLooped=-1 !indicate we've reset from nothing
            else
               iLabelPairIndexOff = ishft(iand(K, 2130706432), -24)
               iLabel1Index = ishft(iand(K, 16773120), -12)
               iLabel2Index = iand(K, 4095)
            endif
            iLabelPairBase=SymPairProds(iTo)%nIndex
            iLabel1=SymStatePairs(1,iLabelPairBase+iLabelPairIndexOff)
            iLabel2=SymStatePairs(2,iLabelPairBase+iLabelPairIndexOff)
! We cycle through these with Label2Index fastest
            iLabel2Index=iLabel2Index+1
! This is a do rather than an if because if iLabel1==iLabel2, it's possible
!  that there are no valid combinations after an increment, so we have to go round again.
            if(iLabel2Index.ge.SymLabelCounts(2,iLabel2)) then
! Increment label 1
               iLabel1Index=iLabel1Index+1
               if(iLabel1Index.ge.SymLabelCounts(2,iLabel1)) then
! Increment the Pair Index
                  iLabelPairIndexOff=iLabelPairIndexOff+1
                  if(iLabelPairIndexOff.ge.SymPairProds(iTo)%nPairs)
     &                then
                     iLooped=1 ! indicate we're looping
                     iLabelPairIndexOff=0
                  endif
                  iLabel1=SymStatePairs(1,
     &               iLabelPairBase+iLabelPairIndexOff)
                  iLabel2=SymStatePairs(2,
     &               iLabelPairBase+iLabelPairIndexOff)
                  iLabel1Index=0
               endif
!  If the two are from the same symlabel, we enforce 2's index>1
               if(iLabel1.eq.iLabel2) then
                  iLabel2Index=iLabel1Index
               else
                  iLabel2Index=0
               endif
            endif
! Now extract the states and turn them into orbitals.
            iTo1=(SymLabelList(SymLabelCounts(1,iLabel1)+iLabel1Index))
     &            *2-1 
            iTo2=(SymLabelList(SymLabelCounts(1,iLabel2)+iLabel2Index))
     &            *2-1 
! Put the K index back together
            K=IShft(iLabelPairIndexOff,24)
     &            +IShft(iLabel1Index,12)+iLabel2Index
            if(tDebugPrint) THEN
               Write(6,*) "To: ",iTo, "LabelPair:", iLabelPairBase,"+",
     &            iLabelPairIndexOff,
     &            "/",SymPairProds(iTo)%nPairs
               Write(6,"(A10,I4,A1,I4,A1,I4,I4)") "Label1: (",iLabel1,
     &            ")",iLabel1Index, "/",
     &            SymLabelCounts(2 ,iLabel1), iTo1
               Write(6,"(A10,I4,A1,I4,A1,I4,I4)") "Label2: (",iLabel2,
     &            ")",iLabel2Index, "/",
     &            SymLabelCounts(2 ,iLabel2), iTo2
               WRITE(6,"(A,Z10,A,I3)") "K:",K," Looped:",iLooped
            endif
         else
! We have a stored list, so just get the next from there.
            K=K+1
            IF(K.lt.0.or.
     &           K.GE.SymPairProds(ITO)%nIndex+SymPairProds(ITO)%nPairs)
     &            THEN
               iLooped=1 !we've looped
               if(K.lt.0) iLooped=-1 ! we've reset
               K=SymPairProds(ITO)%nIndex
            ENDIF
            iTo1=SymStatePairs(1,K)*2-1
            iTo2=SymStatePairs(2,K)*2-1
            if(tDebugPrint) THEN
               Write(6,*) "To: ",iTo
               WRITE(6,"(A10,I4,A1,I4)") "Label1: (",
     &            SymClasses((iTo1+1)/2),")",iTo1
               WRITE(6,"(A10,I4,A1,I4)") "Label2: (",
     &            SymClasses((iTo2+1)/2),")",iTo2
               WRITE(6,"(A,I10,A,I3)") "K:",K," Looped:",iLooped
            endif
         ENDIF
         RETURN
      End Subroutine SymGenExcitIt_GetNextPair

      Subroutine SymGenExcitIt_GenDouble(iRet,SPP,
     &   I,K,L,iTo,iSpn,iFrom,iCC,iLUT,LS,OrbPairs,SymProdInd,
     &   nI,nK,nEl,iC,tDebugPrint)
         use SystemData , only : TSTOREASEXCITATIONS
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymLabelCounts,SymStatePairs,SymClass
         use SymData, only: SymLabelList,SymPairProds
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         
         INTEGER I,K,L, iTo, iFrom, iSpn
         INTEGER NK(NEL)
         TYPE(Symmetry) SPP
         INTEGER ICC(4)
C.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.
         LOGICAL LS(2,2)

         INTEGER ORBPAIRS(2,*)
         INTEGER ILUT(0:*)
         INTEGER SYMPRODIND(2,3,1:*)
         INTEGER IC
         LOGICAL tDebugPrint

         INTEGER iRet
         INTEGER iLooped
         INTEGER iTo1,iTo2
         iRet=0
C.. doubles
         SPP=SymPairProds(ITO)%Sym
C.. See if we need a new K.  L is the spin label, and goes from 1..4.
         IF(L.GT.4) THEN
            L=1
            Call SymGenExcitIt_GetNextPair(K,iTo,iLooped,iTo1,iTo2,
     &         tDebugPrint)
            IF(iLooped.ne.0) THEN
C.. K is invalid.  Get a new I
               I=I+1
               IF(I.GT.SYMPRODIND(2,ISPN,IFROM)) THEN
C.. I is now invalid
                  I=-1
                  iRet=1
                  RETURN
               ENDIF
            ENDIF
C.. We've got a new K, so we need to reset some variables
            ICC(1)=iTo1
            ICC(2)=iTo1+1
            ICC(3)=iTo2
            ICC(4)=iTo2+1
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
               Call SymGenExcitIt_MakeDouble(
     &            ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &            ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &            ICC(1),ICC(3),nI,nK,nEl,tDebugPrint)
               IC=2
               iRet=2
               RETURN
            ENDIF
         ENDIF
         IF(L.EQ.2) THEN
            L=3
            IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
               IF(.NOT.(LS(1,1).OR.LS(2,2))) THEN
                  Call SymGenExcitIt_MakeDouble(
     &             ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ICC(1),ICC(4),nI,nK,nEl,tDebugPrint)
                  IC=2
                  iRet=2
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         IF(L.EQ.3) THEN
            L=4
            IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
              IF(.NOT.(LS(1,2).OR.LS(2,1)).AND.ICC(1).NE.ICC(3)) THEN
                Call SymGenExcitIt_MakeDouble(
     &             ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ICC(2),ICC(3),nI,nK,nEl,tDebugPrint)
                IC=2
                iRet=2
                RETURN
               ENDIF
            ENDIF
         ENDIF
         IF(L.EQ.4) THEN
            L=5
            IF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
              IF(ICC(1).NE.ICC(3).AND..NOT.(LS(1,2).OR.LS(2,2))) THEN
                Call SymGenExcitIt_MakeDouble(
     &             ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I),
     &             ICC(2),ICC(4),nI,nK,nEl,tDebugPrint)
                IC=2
                iRet=2
                RETURN
               ENDIF
            ENDIF
         ENDIF
      End Subroutine SymGenExcitIt_GenDouble
! Actually go through and generate the excitation given a det and from and to.
!  nI is the reference, nK has the new det returned in it and nEl is the # electrons.
      Subroutine SymGenExcitIt_MakeDouble(iFrom1,iFrom2,
     &      iTo1,iTo2,nI,nK,nEl,tDebugPrint)
         use SystemData , only : TSTOREASEXCITATIONS
         IMPLICIT NONE
         INTEGER iFrom1,iFrom2
         INTEGER iTo1,iTo2,ExcitMat(2,2)
         INTEGER nEl,nI(nEl),nK(nEl)
         LOGICAL tDebugPrint,TParity
         INTEGER J
         IF(tStoreAsExcitations) THEN
!The excitation storage starts with -1.  The next number is the excitation level,L .  Next is the parity of the permutation required to lineup occupied->excited.  Then follows a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
            NK(1)=-1
            NK(2)=2
            NK(3)=1
            NK(4)=iFrom1
            NK(5)=iFrom2
            NK(6)=iTo1
            NK(7)=iTo2
         ELSE
!            CALL NECI_ICOPY(NEL,NI,1,NK,1)
            NK(1:NEL)=NI(1:NEL)
            DO J=1,NEL
               IF(NI(J).EQ.iFrom1) THEN
                  IF(tDebugPrint)
     &               WRITE(6,"(A,I3,A,I3,A)",advance='no') 
     &                  "[",iFrom1,"->",
     &                   iTo1,","
                  ExcitMat(1,1)=J
                  ExcitMat(2,1)=iTo1
!                  NK(J)=iTo1
               ENDIF
               IF(NI(J).EQ.iFrom2) THEN
                  IF(tDebugPrint)
     &               WRITE(6,*) iFrom2,"->",iTo2,"]"
                  ExcitMat(1,2)=J
                  ExcitMat(2,2)=iTo2
!                  NK(J)=iTo2
               ENDIF
            ENDDO
C                     CALL WRITEDET(6,NK,NEL,.TRUE.)
!            CALL NECI_SORTI(NEL,NK)
            CALL FindExcitDet(ExcitMat,NK,2,TParity)
         ENDIF
      End Subroutine SymGenExcitIt_MakeDouble
      

      LOGICAL FUNCTION IsConnectedDet(nI,nJ)
         use SystemData, only: nel
         IMPLICIT NONE
         INTEGER nI(nEl), nJ(nEl)
         CALL IsConnectedDet2(nI,nJ,IsConnectedDet)
      END

!  We are required to enshroud GenExcitProbInternal as F90 doesn't seem to want to allow any casting.
      SUBROUTINE IsConnectedDet2(nI,nJ,tIsConnected)
         use SystemData, only: nel
         use symexcit2
         IMPLICIT NONE
         INTEGER nI(nEl),nJ(nEl)
         LOGICAL tIsConnected
         CALL IsConnectedDetInternal(nI,nJ,tIsConnected)
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
!     SymProds(0) = DSTORE(SymClassSize*NEL+1+NIfTot+1))
!      SUBROUTINE SYMGENRANDEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,       &
!     &               CLASSES,                                           &
!     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,    &
!     &               I,J,K,L,ICC,LS,ITOTAL,ISEED,                       &
!     &               NK,IC,G1,NBASISMAX,UMAT,NBASIS,PGEN)               


!  We are required to enshroud GenExcitProbInternal as F90 doesn't seem to want to allow any casting.
      SUBROUTINE GenExcitProb2(nI,nJ,nEl,G1,nBasisMax,Arr,
     &      nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,
     &      pGen)
         use constants, only: dp
         USE symexcit2
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER nEl,nI(nEl),nJ(nEl),nBasis,nBasisMax(5,*)
         INTEGER ExcitTypes(5,*)
         TYPE(BasisFN) G1(nBasis) 
         INTEGER iLUT(*),iTotal
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,1:*)
         TYPE(Symmetry) SymProds(0:*)
         real(dp) pGen
         real(dp) Arr(nBasis,2)
         CALL GenExcitProbInternal(nI,nJ,nEl,G1,nBasisMax,Arr,
     &  nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,
     &   iTotal,pGen)
      END


      SUBROUTINE DUMPIMEMORY(IUNIT,MEM,INLEN)
         IMPLICIT NONE
         INTEGER IUNIT, INLEN!,I,J
         INTEGER*1 MEM(*)!,CH
!         INTEGER ILEN,IPOST,LEN,ILEAD
!         INTEGER ADDR,LEAD,POST,LINES
!#ifdef POINTER8
!         integer(int64) LOC
!#else
!         INTEGER LOC
!#endif
!         INTEGER,PARAMETER :: LINELEN=16
         call stop_all('DUMPIMEMORY',
     &   'Not functional with CRAY compilers.  Fix if you want to use!')

         ! Stop warnings
         iunit=iunit; inlen=inlen; mem(1)=mem(1)
!         ILEN=INLEN*4
!         LEN=ILEN
!!         ADDR=16*(LOC(MEM)/LINELEN)
!         ILEAD=LOC(MEM)-ADDR
!         IPOST=MOD(LEAD+LEN,LINELEN)
!         IF(IPOST.NE.0) IPOST=16-IPOST
!         LINES=(ILEAD+ILEN+IPOST)/LINELEN
!         DO I=0,LINES-1
!            POST=IPOST
!            LEN=ILEN
!            LEAD=ILEAD
!            WRITE(IUNIT,"(Z8.8)",advance='no') ADDR
!            DO J=0,LINELEN-1
!               IF(J.EQ.LINELEN/2) WRITE(IUNIT,'(A2)',advance='no') "  "
!               IF(LEAD.EQ.0) THEN
!                  IF(LEN.EQ.0) THEN
!                     WRITE(IUNIT,"(A3)",advance='no') "   "
!                  ELSE
!                     WRITE(IUNIT,"(Z3.2)",advance='no')
!     &                                              ,MEM(I*LINELEN+J+1)
!                     LEN=LEN-1
!                  ENDIF
!               ELSE
!                  WRITE(IUNIT,"(A3)",advance='no') "   "
!                  LEAD=LEAD-1
!               ENDIF 
!            ENDDO
!            WRITE(IUNIT,"(A3)",advance='no'),"   "
!            DO J=0,LINELEN-1
!               IF(J.EQ.LINELEN/2) WRITE(IUNIT,'(A1)',advance='no') " "
!               IF(ILEAD.EQ.0) THEN
!                  IF(ILEN.NE.0) THEN
!                     CH=MEM(I*LINELEN+J+1)
!                     IF(CH.GE.32) THEN
!                        WRITE(IUNIT,"(A1)",advance='no'),CHAR(CH)
!                     ELSE
!                        WRITE(IUNIT,"(A1)",advance='no'),"."
!                     ENDIF
!                     ILEN=ILEN-1
!                  ELSE
!                     WRITE(IUNIT,"(A1)",advance='no')," "
!                  ENDIF
!               ELSE
!                  WRITE(IUNIT,"(A1)",advance='no')," "
!                  ILEAD=ILEAD-1
!               ENDIF
!            ENDDO
!            ILEAD=LEAD
!            IPOST=POST
!            ILEN=LEN
!            WRITE(IUNIT,*)
!            ADDR=ADDR+LINELEN
!         ENDDO
      END

