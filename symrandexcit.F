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
!         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
!     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
!     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
!     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
!     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
         CALL SymGenRandExcitIt2(NI,NEL,NMEM(NMEM(2)),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         nMem(23),ISEED,NJ,IC,PGEN)
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
!         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
!     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
!     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
!     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
!     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
         CALL SymGenRandExcitIt2(NI,NEL,NMEM(NMEM(2)),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         nMem(23),ISEED,NJ,IC,PGEN)
         ICOUNT=NMEM(23)
!         call writedet(6,nI,nEl,.true.)
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
      IMPLICIT NONE
      INTEGER NBASIS
      INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES),NEXCITTYPES
      
      INTEGER I,J,K,L
      INTEGER NK(NEL)
      TYPE(SymClass) CLASSES(*)
      INTEGER IEXCIT
      INTEGER IFROM,ITO,ISPN
      INTEGER IFROMSL,ITOSL,ExcitMat(2,2)
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
!                     CALL NECI_ICOPY(NEL,NI,1,NK,1)
                     NK(1:NEL)=NI(1:NEL)
                     ExcitMat(1,1)=I
                     ExcitMat(2,1)=ITOSL
                     CALL FindExcitDet(ExcitMat,NK,1)
!                     NK(I)=ITOSL
!                     CALL NECI_SORTI(NEL,NK)
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
         CALL N_MEMORY(IP_ExcitWeight,ExcitWeightSize*
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
         CALL N_FREEM(IP_ExcitWeight)


!  Allocate memory to store all weights and excitations in this category
         CALL N_MEMORY(IP_ExcitWeight,ExcitWeightSize*
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
!         CALL NECI_ICOPY(NEL,NI,1,NK,1)
         NK(1:NEL)=NI(1:NEL)
         DO J=1,NEL
            IF(NI(J).EQ.ews(K)%I) THEN
               ExcitMat(1,1)=J
               ExcitMat(2,1)=ews(K)%A
!               NK(J)=ews(K)%A
            ENDIF
            IF(NI(J).EQ.ews(K)%J) THEN
               ExcitMat(1,2)=J
               ExcitMat(2,2)=ews(K)%B
!               NK(J)=ews(K)%B
            ENDIF
         ENDDO
!         CALL NECI_SORTI(NEL,NK)
         CALL FindExcitDet(ExcitMat,NK,2)
         IC=2
         CALL N_FREEM(IP_ExcitWeight)
         RETURN
      ENDIF
      END 

!== SymGenRandExcitIt2
!== Based on SymGenRandExcitIt, but modular, and capable of acceleration
!== for abelian symmetry groups, where lists of symmetry pairs needn't be
!== stored.

      SUBROUTINE SymGenRandExcitIt2(NI,NEL,EXCITTYPES,                  &
     &               CLASSES,                                           &
     &               SYMPRODIND,ILUT,ORBPAIRS,iTotal,ISEED,             &
     &               NK,IC,PGEN)      
         USE symexcit2
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClass,SymLabelCounts,SymLabelList
         IMPLICIT NONE
         INTEGER nEl,nI(nEl)
         INTEGER I
         INTEGER iSpn,iFrom,iTo
         INTEGER iC
         REAL*8 pGen
         INTEGER ILUT(0:*)
         TYPE(SymClass) CLASSES(*)
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit
         INTEGER OrbPairs(2,*)
         INTEGER iTotal
         INTEGER SymProdInd(2,3,0:*)
         INTEGER nK(nEl)
         INTEGER iSeed
         
         REAL*8 RAN2   
!.. First generate a random excitation 
         PGEN=1.D0/ITOTAL
         I=ITOTAL*RAN2(ISEED)
!.. Now work out which Excitation Type it's in.
         IEXCIT=1
         DO WHILE(EXCITTYPES(5,IEXCIT).LE.I)
            I=I-EXCITTYPES(5,IEXCIT)
            IEXCIT=IEXCIT+1
         ENDDO
!..   iExcit is the Excitation Type
         IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
!.. It's a single
            ISPN=EXCITTYPES(2,IEXCIT)-2
            IFROM=CLASSES(EXCITTYPES(3,IEXCIT))%SymLab
            ITO=EXCITTYPES(4,IEXCIT)
            Call SymGenRandExcitIt2_Single_GetFrom(iSeed,nI,nEl,iSpn,
     &          iFrom,I)
            Call SymGenRandExcitIt2_Single_GetTo(iSeed,nI,nEl,I,iSpn,
     &         iTo,nK)
!.. Indicate we're a single
            iC=1
         ELSE
!..  A double
!.. We store each excitation type as:
!.. 1   TYPE (single=1, double=2)
!.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
!.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
!.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
!.. 5  COUNT (Total number of excitations in this category)
            ISPN=EXCITTYPES(2,IEXCIT)
            IFROM=EXCITTYPES(3,IEXCIT)
            ITO=EXCITTYPES(4,IEXCIT)
            Call SymGenRandExcitIt2_Double_GetFrom(iSeed,ExcitTypes,
     &       iExcit,iSpn, iFrom, OrbPairs, SymProdInd, iTotal, I,pGen)

            Call SymGenRandExcitIt2_Double_GetTo(iSeed,nI,ExcitTypes,
     &         iExcit,OrbPairs, SymProdInd,iSpn,iFrom,iLUT, I,nK,pGen)
!.. Indicate we're a double.
            iC=2
         ENDIF
      End Subroutine SymGenRandExcitIt2

!Return the index of orbital to excit within nI into I. i.e. nI(I) will be excited.
      Subroutine SymGenRandExcitIt2_Single_GetFrom(iSeed,nI,nEl,iSpn,
     &       iFrom,I)
         use SymData, only: SymClass,SymLabelCounts,SymLabelList
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nEl,nI(nEl)
         INTEGER iSpn, iFrom
         INTEGER I

         INTEGER iCount
         INTEGER iFromSL
         INTEGER J,K
         LOGICAL L1
         REAL*8 RAN2   
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
      End Subroutine SymGenRandExcitIt2_Single_GetFrom

!.. Locate an orbital to excite to, and create the excitation, copying it back to nI.
!..  I is the index in I to excite from.
      Subroutine SymGenRandExcitIt2_Single_GetTo(iSeed,nI,nEl,I,iSpn,
     &   iTo,nK)
      use SymData, only: SymClass,SymLabelCounts,SymLabelList
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nEl,nI(nEl)
         INTEGER I,ExcitMat(2,2)
         INTEGER iSpn,iTo
         INTEGER iC

         INTEGER iCount
         INTEGER iToSL
         INTEGER J,K,L
         INTEGER NK(nEl)
         LOGICAL L2
      
         LOGICAL IsValidDet
         REAL*8 RAN2   
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
!                     CALL NECI_ICOPY(NEL,NI,1,NK,1)
                     NK(1:NEL)=NI(1:NEL)
                     ExcitMat(1,1)=I
                     ExcitMat(2,1)=ITOSL
!                     NK(I)=ITOSL
!                     CALL NECI_SORTI(NEL,NK)
                     CALL FindExcitDet(ExcitMat,NK,1)
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
                     RETURN
                  ELSE
                     ICOUNT=ICOUNT-1
                  ENDIF
               ENDIF
               K=K+1
            ENDDO
            ICOUNT=ICOUNT*RAN2(ISEED)
         ENDDO
      End Subroutine SymGenRandExcitIt2_Single_GetTo

!.. Get the FROM pair.  Return as index I into the part of OrbPairs pointed to
!.. by the SymProdInd correspnding to (iSpn,iFrom) 
      Subroutine SymGenRandExcitIt2_Double_GetFrom(iSeed,ExcitTypes,
     &       iExcit,iSpn, iFrom,
     &      OrbPairs, SymProdInd, iTotal, 
     &      I, pGen)
         USE symexcit2
         use SystemData, only: G1, nBasis, nBasisMax
         use SystemData, only: Arr
         use IntegralsData, only: UMat
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER iSpn
         INTEGER iFrom
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit
         INTEGER OrbPairs(2,*)
         INTEGER iTotal
         INTEGER SymProdInd(2,3,0:*)
         INTEGER I
         REAL*8 pGen
         
         REAL*8 rCount
         TYPE(ExcitWeight) ews(*)
         POINTER (IP_ExcitWeight,ews)
         INTEGER iCount
         REAL*8 RAN2
!.. Go through the list of pairs with a given symprod.
!.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
!.. SYMPRODS(I) in ORBPAIRS
!.. SYMPRODIND(2,ISPN,I) contains the number of such elements
         pGen=(0.D0+ExcitTypes(5,iExcit))/iTotal !pGen is the prob of choosing iExcit
!  Allocate memory to store all weights and excitations in this category
         CALL N_MEMORY(IP_ExcitWeight,ExcitWeightSize*
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
         CALL  N_FREEM(IP_ExcitWeight)
      End Subroutine SymGenRandExcitIt2_Double_GetFrom

      Subroutine SymGenRandExcitIt2_Double_GetTo(iSeed,nI,ExcitTypes,
     &      iExcit,OrbPairs, SymProdInd,iSpn,iFrom,iLUT,
     &      I,nK, pGen)
         USE symexcit2
         use SystemData, only: G1, nBasis, nBasisMax
         use SystemData, only: Arr,nEl
         use IntegralsData, only: UMat
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nI(nEl)
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit,ExcitMat(2,2)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,0:*)
         INTEGER iSpn,iFrom
         INTEGER iLUT(0:*)
         INTEGER I
         REAL*8 pGen
        
         INTEGER NK(nEl) 
         TYPE(ExcitWeight) ews(*)
         POINTER (IP_ExcitWeight,ews)
         REAL*8 rCount
         INTEGER iCount
         INTEGER K,J
         REAL*8 RAN2
!  Allocate memory to store all weights and excitations in this category
         CALL N_MEMORY(IP_ExcitWeight,ExcitWeightSize*
     &      EXCITTYPES(5,IEXCIT)/SYMPRODIND(2,ISPN,IFROM),'ExcitWeigh')

!.. Now go through the list of virtual pairs, ng those with orbitals in NI 
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
!         CALL NECI_ICOPY(NEL,NI,1,NK,1)
         NK(1:NEL)=NI(1:NEL)
         DO J=1,NEL
            IF(NI(J).EQ.ews(K)%I) THEN
               ExcitMat(1,1)=J
               ExcitMat(2,1)=ews(K)%A
!               NK(J)=ews(K)%A
            ENDIF
            IF(NI(J).EQ.ews(K)%J) THEN
               ExcitMat(1,2)=J
               ExcitMat(2,2)=ews(K)%B
!               NK(J)=ews(K)%B
            ENDIF
         ENDDO
!         CALL NECI_SORTI(NEL,NK)
         CALL FindExcitDet(ExcitMat,NK,2)
         CALL  N_FREEM(IP_ExcitWeight)
         RETURN
      End Subroutine SymGenRandExcitIt2_Double_GetTo
