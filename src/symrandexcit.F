!This wrapper also passes back the excitation matrix and parity of excitation permutation, so that matrix elements can be worked out simply between the original determinant
!and the excitation.
      SUBROUTINE GenRandSymExcitIt4(NI,NMEM,NJ,ISEED,IC,PGEN,
     &        ICOUNT,EXCITMAT,PARITY)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         use constants, only: dp,sizeof_int
         use SystemData, only: Arr,nMax,tAssumeSizeExcitgen
         use IntegralsData, only: UMat
         use SymData, only: SymClassSize,nSymPairProds,nSymLabels
         use legacy_data, only: irat
         IMPLICIT NONE
         INTEGER NI(NEL)
         INTEGER,pointer :: DSTORE(:)
         INTEGER EXCITMAT(2,2), parity
         INTEGER,target :: NMEM(*)
         integer NJ(NEL),IC,ICOUNT
         INTEGER ISEED,Pos1,Pos2,Pos3
         real(dp) PGEN

         IF(tAssumeSizeExcitgen) THEN
!Excitgen is set up differently...
! 1         NEXCIT
! 2         NEXCITTYPES
! 3                                 -    SymClassSize*NEL+nBasis/32+3                                 DSTORE
! SymClassSize*NEL+nBasis/32+4      -    SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+2  EXCITTYPES
! SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+3      -    SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+2  ORBPAIRS
! SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+3   -   SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+(nSymLabels+1)*6+2  SYMPRODIND

! DSTORE(1)      -   SymClassSize*NEL         CLASSES
! DSTORE(NEL*SymClassSize+1)                     ILUT

            DSTORE=>NMEM(3:NEL*SymClassSize+4)
             Pos1=SymClassSize*NEL+nBasis/32+4
             Pos2=SymClassSize*NEL+nBasis/32+4+15*nSymPairProds+
     &             10*nSymLabels
             Pos3=SymClassSize*NEL+nBasis/32+4+15*nSymPairProds+
     &             10*nSymLabels+NEL*(NEL-1)
            
            CALL SymGenRandExcitIt2(NI,NEL,NMEM(Pos1),DSTORE(1),
     &         NMEM(Pos3),DSTORE(SymClassSize*NEL+1),NMEM(Pos2),
     &         nMem(1),ISEED,NJ,IC,PGEN,EXCITMAT,PARITY)
            ICOUNT=NMEM(1)

         ELSE
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
            DSTORE=>NMEM(NMEM(1):NMEM(1)+
     &                SymClassSize*NEL+(nBasis/32)+1
     &                     +SymmetrySize*(NEL*NEL+1))  !Point to DSTORE
!         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
!     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
!     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
!     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
!     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
            CALL SymGenRandExcitIt2(NI,NEL,NMEM(NMEM(2)),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         nMem(23),ISEED,NJ,IC,PGEN,EXCITMAT,PARITY)
            ICOUNT=NMEM(23)
         ENDIF
      END SUBROUTINE GenRandSymExcitIt4



C.. Generate a random symmetry allowed excitation.  GENSYMEXCITIT2 must be called twice
C.. before this.
C.. This version is based on GenRandSymExcitIt2, but does not need values it extracts from modules
!  IC is the excitation level, and iCount the total number of allowed excitations
      SUBROUTINE GenRandSymExcitIt3(NI,NMEM,NJ,ISEED,IC,PGEN,
     &      ICOUNT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nEl,G1, nBasis,nBasisMax
         use SystemData, only: Arr,nMax,tAssumeSizeExcitgen
         use IntegralsData, only: UMat
         use SymData, only: SymClassSize,nSymPairProds,nSymLabels
         use legacy_data, only: irat
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NI(NEL)
         INTEGER, pointer :: DSTORE(:)
         INTEGER EXCITMAT(2,2), parity
         INTEGER,target :: NMEM(*)
         integer  NJ(NEL),IC,ICOUNT
         INTEGER ISEED,Pos1,Pos2,Pos3
         real(dp) PGEN
         IF(tAssumeSizeExcitgen) THEN
!Excitgen is set up differently...
! 1         NEXCIT
! 2         NEXCITTYPES
! 3                                 -    SymClassSize*NEL+nBasis/32+3                                 DSTORE
! SymClassSize*NEL+nBasis/32+4      -    SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+2  EXCITTYPES
! SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+3      -    SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+2  ORBPAIRS
! SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+3   -   SymClassSize*NEL+nBasis/32+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+(nSymLabels+1)*6+2  SYMPRODIND

! DSTORE(1)      -   SymClassSize*NEL         CLASSES
! DSTORE(NEL*SymClassSize+1)                     ILUT

             DSTORE=>NMEM(3:NEL*SymClassSize+4)
             Pos1=SymClassSize*NEL+nBasis/32+4
             Pos2=SymClassSize*NEL+nBasis/32+4+15*nSymPairProds+
     &             10*nSymLabels
             Pos3=SymClassSize*NEL+nBasis/32+4+15*nSymPairProds+
     &             10*nSymLabels+NEL*(NEL-1)
            
            CALL SymGenRandExcitIt2(NI,NEL,NMEM(Pos1),DSTORE(1),
     &         NMEM(Pos3),DSTORE(SymClassSize*NEL+1),NMEM(Pos2),
     &         nMem(1),ISEED,NJ,IC,PGEN,EXCITMAT,PARITY)
            ICOUNT=NMEM(1)

         ELSE
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
            DSTORE=>NMEM(NMEM(1):NMEM(1)+
     &                SymClassSize*NEL+(nBasis/32)+1
     &                     +SymmetrySize*(NEL*NEL+1))  !Point to DSTORE
!         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
!     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
!     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
!     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
!     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
            CALL SymGenRandExcitIt2(NI,NEL,NMEM(NMEM(2)),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         nMem(23),ISEED,NJ,IC,PGEN,EXCITMAT,PARITY)
            ICOUNT=NMEM(23)
         ENDIF
      END



C.. Generate a random symmetry allowed excitation.  GENSYMEXCITIT2 must be called twice
C.. before this.
C.. ICOUNT returned is the total number of possible excitations.  To get the excitation level,
C..   use GenRandSymExcitIt3
      SUBROUTINE GENRANDSYMEXCITIT2(NI,NEL,
     &         NMEM,NJ,ISEED,ICOUNT,PGEN)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,nBasis
         use SymData, only: SymClassSize
         use legacy_data, only: irat
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         INTEGER,pointer :: DSTORE(:)
         INTEGER EXCITMAT(2,2), parity
         INTEGER,target :: NMEM(*)
         integer NJ(NEL),ICOUNT,IC
         INTEGER ISEED
         real(dp) PGEN
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
            DSTORE=>NMEM(NMEM(1):NMEM(1)+
     &                SymClassSize*NEL+(nBasis/32)+1
     &                     +SymmetrySize*(NEL*NEL+1))  !Point to DSTORE
!         CALL SYMGENRANDEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
!     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
!     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
!     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NMEM(23),ISEED,NJ,IC,
!     &         G1,NBASISMAX,UMAT,ARR,NBASIS,PGEN)
         CALL SymGenRandExcitIt2(NI,NEL,NMEM(NMEM(2)),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),
     &         nMem(23),ISEED,NJ,IC,PGEN,EXCITMAT,PARITY)
         ICOUNT=NMEM(23)
!         call writedet(6,nI,nEl,.true.)
      END
      SUBROUTINE SYMGENRANDEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,       &
     &               CLASSES,                                           &
     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,    &
     &               I,J,K,L,ITOTAL,ISEED,                              &
     &               NK,IC,NBASISMAX,ARR,NBASIS,PGEN)      
      use constants, only: dp,sizeof_int
      USE symexcit2
      use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
      use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
      use SymData, only: SymClass,SymLabelCounts,SymLabelList
      use Determinants, only: write_det
      IMPLICIT NONE
      INTEGER NBASIS
      INTEGER NEXCITTYPES
      INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES)
      
      INTEGER I,J,K,L
      INTEGER NK(NEL)
      TYPE(SymClass) CLASSES(*)
      INTEGER IEXCIT
      INTEGER IFROM,ITO,ISPN, parity
      INTEGER IFROMSL,ITOSL,ExcitMat(2,2)
      LOGICAL L1,L2
!.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.

      INTEGER ORBPAIRS(2,*)
      INTEGER ILUT(0:*),ITOTAL,ISEED
      INTEGER SYMPRODIND(2,3,1:*)
      INTEGER IC,ICOUNT
      real(dp) RAN2
      LOGICAL ISVALIDDET
      INTEGER nBasisMax(5,*)
      real(dp) RCOUNT
      real(dp) PGEN
      TYPE(ExcitWeight),pointer :: ews(:)
      real(dp) ARR(NBASIS,2)
!.. First generate a random excitation 
      PGEN=1.0_dp/ITOTAL
      I=int(ITOTAL*RAN2(ISEED),sizeof_int)
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
            ICOUNT=int(ICOUNT*RAN2(ISEED),sizeof_int)
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
                     CALL FindExcitDet(ExcitMat,NK,1,Parity)
!                     NK(I)=ITOSL
!                     CALL NECI_SORTI(NEL,NK)
                     IF(.NOT.ISVALIDDET(NK,NEL)) THEN
            
                       WRITE(6,*) "INVALID DET SYMGENRANDEXCITIT SINGLE"
                       call write_det (6, NK, .true.)
                       call write_det (6, NI, .true.)
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
            ICOUNT=int(ICOUNT*RAN2(ISEED),sizeof_int)
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
         pGen=(0.0_dp+ExcitTypes(5,iExcit))/iTotal !pGen is the prob of choosing iExcit
!  Allocate memory to store all weights and excitations in this category
         allocate(ews(SYMPRODIND(2,ISPN,IFROM)))
         CALL EnumExcitFromWeights(EXCITTYPES(1,IEXCIT),ews,OrbPairs,
     &      SymProdInd,
     &      RCOUNT,ICOUNT,ARR,NBASIS)
!  We choose a symprod of this symmetry to excite FROM with probs in ews
         pGen=pGen/rCount  ! The norm factor for the FROMs
         rCount=rCount*RAN2(ISEED)
         I=0
!  This is a little inefficient, and could perhaps be binary searched if cumulative weights were also stored
         DO WHILE (rCount.GT.0.0_dp)
            I=I+1
            rCount=rCount-ews(I)%Weight
         ENDDO   
!  I is the index of the specific SYM PAIR FROM
         pGen=pGen*ews(I)%Weight
         
!  pGen is the prob of choosing a specific FROM (given having chosen iExcit proportional to the number of excitations in each iExcit)
!           times the prob of choosing iExcit
         deallocate(ews)


!  Allocate memory to store all weights and excitations in this category
         allocate(ews(EXCITTYPES(5,IEXCIT)/SYMPRODIND(2,ISPN,IFROM)))

!.. Now go through the list of virtual pairs, excluding those with orbitals in NI 
!.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
!.. states (listed in SymStatePairs), and the index of the start of this list
!.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
!.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose prod is of that symmetry.
         RCOUNT=0.0_dp
         ICOUNT=0
!  L is a loop variable.  We go through the loop once to count the number of possible excitations
!  We then pick a random number up to that number, and choose that excitation
!  This is modified to weight each excitation according to its U-matrix element, so each
!  Excitation's weight is added to RCOUNT to give a normalization.  RCOUNT*RAND then corresponds to 
!  a point in the sum of weights which is found by subtracting weights until we reach 0.
!  This is only done for the excitations TO because we need to know the excitations FROM to generate the weight
         CALL EnumExcitWeights(EXCITTYPES(1,IEXCIT),I,ILUT,ews,OrbPairs,
     &      SymProdInd,RCOUNT,ICOUNT,NBASISMAX,ARR,NBASIS)

         pGen=pGen/rCount           ! Normalize for the sum of TOs
         
!  We've now got a list of iCount weights in ews, and total norm rCount.
!  Select a random one
         rCount=rCount*RAN2(ISEED)
         K=0
!  This is a little inefficient, and could perhaps be binary searched if cumulative weights were also stored
         DO WHILE (rCount.GT.0.0_dp)
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
         CALL FindExcitDet(ExcitMat,NK,2,Parity)
         IC=2
         deallocate(ews)
         RETURN
      ENDIF
      END 


!Return the index of orbital to excit within nI into I. i.e. nI(I) will be excited.
!In the case where we have no weighting, the can be done very simply
!iSpn is -1 for alpha and 0 for beta.
!iFrom is the symmetry label of the electron to pick.
!I is the returned orbital to excite
      Subroutine SymGenRandExcitIt2_Sing_GF_NoW(iSeed,
     &       nEl,iSpn,iFrom,I,nI)
         use SymData, only: SymClass,SymLabelCounts,SymLabelList
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nEl,nI(nEl)
         INTEGER iSpn, iFrom
         INTEGER I
         INTEGER State
         INTEGER SpinOrb
         INTEGER J,OrbChosen,nOrbs
         LOGICAL InvalidOrb
         real(dp) RAN2   

!Algorithm is to simply keep drawing orbitals of the required symmetry at random, until we
!pick one which is in NI.
         InvalidOrb=.true.

!Number of states with allowed symmetry is given by SymLabelCounts(2,Sym)
!Draw an orbital at random from this number, and find its related orbital, given at index SymLabelCounts(1,Sym)+(DrawnOrb-1)
!See if the spin orbital relating to this state is in NI. If not, redraw.
         nOrbs=SymLabelCounts(2,iFrom)
!         WRITE(6,*) "SYM: ",iFrom
!         WRITE(6,*) "NORBS: ",nOrbs
!         WRITE(6,*) "iSpn: ",iSpn

         do while(InvalidOrb)

            OrbChosen=INT(Ran2(iSeed)*nOrbs)+1      !This is the orbital picked at random from the number with desired symmetry
!            WRITE(6,*) "OrbChosen: ",OrbChosen
            State=SymLabelList(SymLabelCounts(1,iFrom)+OrbChosen-1) !This is the index of the state chosen
!            WRITE(6,*) "State: ",State
            SpinOrb=(2*State)+iSpn      !This is the index of the spin-orbital we've chosen
!            WRITE(6,*) "SpinOrb: ",SpinOrb
!            WRITE(6,*) "ILUT: ", ILUT(0:1)

!Now we need to check if the spin-orbital relating to the state is in NI
!This should be binary searched for speed. Even better, we could just to a bit-test on ILUT
            DO j=1,NEl
                IF(nI(j).eq.SpinOrb) THEN
!We've found an allowed orbital
                    InvalidOrb=.false.
                    I=j
!                    IF(.not.BTEST(ILUT((SpinOrb-1)/32)
!     &                       ,mod(SpinOrb-1,32))) THEN
!                        STOP 'ILUT failure'
!                    ENDIF
                    EXIT
                ENDIF
            ENDDO

        enddo

      End Subroutine SymGenRandExcitIt2_Sing_GF_NoW

!Return the index of orbital to excit within nI into I. i.e. nI(I) will be excited.
      Subroutine SymGenRandExcitIt2_Sing_GF(iSeed,nI,nEl,iSpn,
     &       iFrom,I)
         use SymData, only: SymClass,SymLabelCounts,SymLabelList
         use constants, only: dp,sizeof_int
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nEl,nI(nEl)
         INTEGER iSpn, iFrom
         INTEGER I

         INTEGER iCount
         INTEGER iFromSL
         INTEGER J,K
         LOGICAL L1
         real(dp) RAN2   
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
            ICOUNT=int(ICOUNT*RAN2(ISEED),sizeof_int)
         ENDDO ot
!.. NI(I) is now the orbital to excite
      End Subroutine SymGenRandExcitIt2_Sing_GF

!.. Locate an orbital to excite to, and create the excitation, copying it back to nI.
!..  I is the index in NI to excite from.
! This is for the case that there is no excitation weighting
      Subroutine SymGenRandExcitIt2_Sing_GT_NoW(iSeed,nI,nEl,
     &              I,iSpn,iTo,nK,ExcitMat,Parity)
      use SymData, only: SymClass,SymLabelCounts,SymLabelList
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nEl,nI(nEl)
         INTEGER I,ExcitMat(2,2)
         INTEGER iSpn,iTo
         INTEGER OrbChosen,norbs,spinorb
         INTEGER J,State, parity
         INTEGER NK(nEl)
         LOGICAL InDi
         real(dp) RAN2   

!Number of states with allowed symmetry is given by SymLabelCounts(2,Sym)
!Draw an orbital at random from this number, and find its related orbital, given at index SymLabelCounts(1,Sym)+(DrawnOrb-1)
!See if the spin orbital relating to this state is in NI. If it is, redraw.
         nOrbs=SymLabelCounts(2,iTo)

         do while(.true.)

            OrbChosen=INT(Ran2(iSeed)*nOrbs)+1      !This is the orbital picked at random from the number with desired symmetry
            State=SymLabelList(SymLabelCounts(1,iTo)+OrbChosen-1) !This is the index of the state chosen
            SpinOrb=(2*State)+iSpn      !This is the index of the spin-orbital we've chosen

!Now we need to check if the spin-orbital relating to the state is in NI
!This should be binary searched for speed. Even better, we could just to a bit-test on ILUT
            InDi=.false.
            DO j=1,NEl
                IF(nI(j).eq.SpinOrb) THEN
!We've found the same orbital in Di - this is not a valid virtual to excite to.
                    InDi=.true.
                    EXIT
                ENDIF
            ENDDO
!            WRITE(6,*) "NI: ",NI(:)
!            WRITE(6,*) "Orb to vacate: ",NI(I)
!            WRITE(6,*) "Orb to occupy: ",SpinOrb
            IF(.not.InDi) THEN
                NK(1:NEl)=NI(1:NEl)
                ExcitMat(1,1)=I
                ExcitMat(2,1)=SpinOrb
                CALL FindExcitDet(ExcitMat,NK,1,Parity)
!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
!                Excitlevel=iGetExcitLevel(NI,NK,NEL)
!                IF(Excitlevel.ne.1) THEN
!                    WRITE(6,*) "Have not created a single excitation"
!                   CALL WRITEDET(6,NI,NEL,.TRUE.)
!                   CALL WRITEDET(6,NK,NEL,.TRUE.)
!                   WRITE(6,*) I,SpinOrb
!                    STOP "Have not created a single excitation"
!                ENDIF
!                IF(.NOT.ISVALIDDET(NK,NEL)) THEN
!                   WRITE(6,*) "INVALID DET SYMGENRANDEXCITIT SINGLE"
!                   CALL WRITEDET(6,NI,NEL,.TRUE.)
!                   CALL WRITEDET(6,NK,NEL,.TRUE.)
!                   WRITE(6,*) I,SpinOrb
!                   STOP "INVALID DET SYMGENRANDEXCITIT SINGLE"
!                ENDIF
                EXIT
            ENDIF

        enddo

      End Subroutine SymGenRandExcitIt2_Sing_GT_NoW


!.. Locate an orbital to excite to, and create the excitation, copying it back to nI.
!..  I is the index in I to excite from.
      Subroutine SymGenRandExcitIt2_Single_GetTo(iSeed,nI,nEl,I,iSpn,
     &   iTo,nK,ExcitMat,Parity)
      use SymData, only: SymClass,SymLabelCounts,SymLabelList
        use Determinants, only: write_det
         use constants, only: dp,sizeof_int
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nEl,nI(nEl)
         INTEGER I,ExcitMat(2,2)
         INTEGER iSpn,iTo, parity
         INTEGER iCount
         INTEGER iToSL
         INTEGER J,K,L
         INTEGER NK(nEl)
         LOGICAL L2, lTmp
      
         LOGICAL IsValidDet
         real(dp) RAN2   
! Now go through the list of possible excitations of it
         ICOUNT=0
!  J is a loop variable.  We go through the loop once to count the number of possible excitations
!  We then pick a random number up to that number, and choose that excitation
         DO J=0,1
            K=0
            L=0
            DO WHILE(K.LT.SYMLABELCOUNTS(2,ITO))
               ITOSL=2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)+ISPN
               L2 = .true.
               if (L /= 0) then
                   if (nI(L) >= iToSL) L2 = .false.
               endif
               DO WHILE (L2)
                  L=L+1
                  !L2=L.LE.NEL.AND.NI(L).LT.ITOSL
                  L2 = .false.
                  if (L <= NEL) then
                      if (nI(L) < itoSL) L2 = .true.
                  endif
               ENDDO
               lTmp = .true.
               if (L <= nel) then
                   if (nI(L) == iToSL) lTmp = .false.
               endif
               IF(lTmp) THEN
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
                     CALL FindExcitDet(ExcitMat,NK,1,Parity)
                     IF(.NOT.ISVALIDDET(NK,NEL)) THEN
            
                       WRITE(6,*) "INVALID DET SYMGENRANDEXCITIT SINGLE"
                       call write_det (6, NK, .true.)
                       call write_det (6, NI, .true.)
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
            ICOUNT=int(ICOUNT*RAN2(ISEED),sizeof_int)
         ENDDO
      End Subroutine SymGenRandExcitIt2_Single_GetTo

!.. Get the FROM pair.  Return as index I into the part of OrbPairs pointed to
!.. by the SymProdInd correspnding to (iSpn,iFrom) 
!.. ExcitTypes is the whole excittypes array.iExcit the label the the interested type.
!.. iSpn is the spin (1, 2 or 3)
!.. iFrom is the symmetry label of the orbital pair to choose
!.. OrbPairs is the list of all occupied pairs, ordered by symmetry and spin.
!.. iTotal is the total number of excitations.
      Subroutine SymGenRandExcitIt2_Doub_GF_NoW(iSeed,
     &       iSpn, iFrom,
     &       SymProdInd,I,nI,OrbPairs)
         use SystemData , only: NEl
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iSeed,nI(NEl)
         INTEGER iSpn
         INTEGER iFrom,OrbPairs(2,*)
         INTEGER SymProdInd(2,3,1:*)
         INTEGER I,k,j
         LOGICAL found
         real(dp) RAN2
!!.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
!!.. SYMPRODS(I) in ORBPAIRS
!!.. SYMPRODIND(2,ISPN,I) contains the number of such elements
!.. Get the FROM pair.  Return as index I into the part of OrbPairs pointed to
!.. by the SymProdInd correspnding to (iSpn,iFrom) 
!There are SymProdInd(2,iSpn,iFrom) excitations with the desired spin and symmetry.
!There should be no reason why one should be preferred over another, therefore, just pick at random

          I=INT(RAN2(iSeed)*SymProdInd(2,iSpn,iFrom))+1
!         WRITE(6,*) "OrbPair: ", OrbPairs(1,I+SymProdInd(1,iSpn,iFrom)),
!     &            OrbPairs(2,I+SymProdInd(1,iSpn,iFrom))
          do j=1,2
              found=.false.
              do k=1,NEl
               IF(nI(k).eq.OrbPairs(j,I+SymProdInd(1,iSpn,iFrom))) THEN
                      found=.true.
                      exit
                  ENDIF
              enddo
              IF(.not.found) THEN
                  WRITE(6,*) "Cannot locate found occupied orbital pair"
                  WRITE(6,*) iSpn,iFrom,
     &               OrbPairs(j,I+SymProdInd(1,iSpn,iFrom))
                  WRITE(6,*) nI(:)
                  STOP 'Cannot locate found occupied orbital pair'
              ENDIF
          enddo
!          WRITE(6,*) "NI= ",NI(:)
!          WRITE(6,*) "i,j = ",OrbPairs(1,I+SymProdInd(1,iSpn,iFrom)),
!     &          OrbPairs(2,I+SymProdInd(1,iSpn,iFrom))

      End Subroutine SymGenRandExcitIt2_Doub_GF_NoW


      Subroutine SymGenRandExcitIt2_Doub_GF(iSeed,ExcitTypes,
     &       iExcit,iSpn, iFrom,
     &      OrbPairs, SymProdInd, iTotal, 
     &      I, pGen)
         USE symexcit2
         use SystemData, only: G1, nBasis, nBasisMax
         use SystemData, only: Arr
         use IntegralsData, only: UMat
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER iSpn
         INTEGER iFrom
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit
         INTEGER OrbPairs(2,*)
         INTEGER iTotal
         INTEGER SymProdInd(2,3,1:*)
         INTEGER I
         real(dp) pGen
         
         real(dp) rCount
         TYPE(ExcitWeight),pointer :: ews(:)
         INTEGER iCount
         real(dp) RAN2

!         WRITE(6,*) "IEXCIT: ",IExcit
!         WRITE(6,*) "ExcitTypes2",ExcitTypes(2,1:4)
!         WRITE(6,*) "ExcitTypes3",ExcitTypes(3,1:4)
!         WRITE(6,*) "ExcitTypes4",ExcitTypes(4,1:4)
!         WRITE(6,*) "ExcitTypes5",ExcitTypes(5,1:4)
!         WRITE(6,*) "SymProdInd(1,1,:): ",SymProdInd(1,1,1:8)
!         WRITE(6,*) "SymProdInd(1,2,:): ",SymProdInd(1,2,1:8)
!         WRITE(6,*) "SymProdInd(1,3,:): ",SymProdInd(1,3,1:8)
!         WRITE(6,*) "SymProdInd(2,1,:): ",SymProdInd(2,1,1:8)
!         WRITE(6,*) "SymProdInd(2,2,:): ",SymProdInd(2,2,1:8)
!         WRITE(6,*) "SymProdInd(2,3,:): ",SymProdInd(2,3,1:8)
!         WRITE(6,*) "***"
!         WRITE(6,*) "OrbPairs: ",OrbPairs(:,0:26)

!.. Go through the list of pairs with a given symprod.
!.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
!.. SYMPRODS(I) in ORBPAIRS
!.. SYMPRODIND(2,ISPN,I) contains the number of such elements
         pGen=(0.0_dp+ExcitTypes(5,iExcit))/iTotal !pGen is the prob of choosing iExcit
!  Allocate memory to store all weights and excitations in this category
         allocate(ews(SYMPRODIND(2,ISPN,IFROM)))
         CALL EnumExcitFromWeights(EXCITTYPES(1,IEXCIT),ews,OrbPairs,
     &      SymProdInd,
     &      RCOUNT,ICOUNT,ARR,NBASIS)
!  We choose a symprod of this symmetry to excite FROM with probs in ews
         pGen=pGen/rCount  ! The norm factor for the FROMs
         rCount=rCount*RAN2(ISEED)
         I=0
!  This is a little inefficient, and could perhaps be binary searched if cumulative weights were also stored
         DO WHILE (rCount.GT.0.0_dp)
            I=I+1
            rCount=rCount-ews(I)%Weight
         ENDDO   
!  I is the index of the specific SYM PAIR FROM
         pGen=pGen*ews(I)%Weight
         
!         WRITE(6,*) "OrbPair: ", OrbPairs(1,I+SymProdInd(1,iSpn,iFrom)),
!     &            OrbPairs(2,I+SymProdInd(1,iSpn,iFrom))
         
         
!  pGen is the prob of choosing a specific FROM (given having chosen iExcit proportional to the number of excitations in each iExcit)
!           times the prob of choosing iExcit
         deallocate(ews)
      End Subroutine SymGenRandExcitIt2_Doub_GF
      
!This is a routine to pick an allowed virtual pair of orbitals to excite to.
!There is no weighting function, so it simply 
      Subroutine SymGenRandExcitIt2_Doub_GT_NoW(iSeed,nI,
     &      ExcitTypes,iExcit,OrbPairs, SymProdInd,iSpn,iFrom,iLUT,
     &      I,nK,ExcitMat,Parity)
         use SystemData, only: nEl
         use SymData, only: tAbelianFastExcitgen
         use SymData, only: SymPairProds,SymStatePairs
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iSeed
         INTEGER nI(nEl)
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit,ExcitMat(2,2)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,1:*)
         INTEGER iSpn,iFrom,iTo,FinalA,FinalB
         INTEGER iLUT(0:*)
         INTEGER I,StatePair,Index1,Orb1,Orb2,Orb1Alpha,Orb1Beta
         INTEGER NK(nEl),Orb2Alpha,Orb2Beta,OrbAlpha,OrbBeta,parity
         LOGICAL :: NotSelectedPair,Check1,Check2
         INTEGER J
         real(dp) RAN2

         iTo=ExcitTypes(4,iExcit)   !This is the symmetry of the virtual pair type

         do while(.true.)   !Cycle until we find an acceptable virtual pair to excite to.
             
!We have an explicit list of state pairs and their symmetry, so we can just pick one at random.
!The number of pairs with desired symmetry is given by SymPairsProds(ExcitTypes(4,iExcit))%nPairs.
             IF((iSpn.eq.1).or.(iSpn.eq.3)) THEN
!Choose a state pair from the list at random
                StatePair=INT(Ran2(iSeed)*
     &              SymPairProds(iTo)%nPairs)+1
!Orb1/2 are the is the beta orbital corresponding to the a state in the state pair with the desired symmetry product             
                Index1=SymPairProds(iTo)%nIndex+StatePair-1

                Orb1=SymStatePairs(1,Index1)*2-1
                Orb2=SymStatePairs(2,Index1)*2-1

                Orb1Beta=Orb1
                Orb2Beta=Orb2
                Orb1Alpha=Orb1+1
                Orb2Alpha=Orb2+1

             ELSE
!We want to pick an alpha-beta spin orbtial combo without biasing the pick.
!The difficulty is due to the fact that some state pairs will have two possible ab combinations (different states)
!but some will have only one, if they are from the same state. To get around this, reject the choice with Prob=0.5 if
!we pick a state pair with two same states. Otherwise, just choose an ab state at random from the states chosen.
                NotSelectedPair=.true.
                do while(NotSelectedPair)

                    StatePair=INT(Ran2(iSeed)*
     &                 SymPairProds(iTo)%nPairs)+1

                    Index1=SymPairProds(iTo)%nIndex+StatePair-1

                    Orb1=SymStatePairs(1,Index1)*2-1
                    Orb2=SymStatePairs(2,Index1)*2-1

                    IF(Orb1.eq.Orb2) THEN
!We have a same-state pair - reject with prob=0.5
                        IF(Ran2(iSeed).le.0.5) THEN
!Accept the state pair
                            NotSelectedPair=.false.
                            OrbAlpha=Orb1+1
                            OrbBeta=Orb1
                        ENDIF
                    ELSE
!We have different states in our state pair. Choose whether to take the ab pair, or ba pair.
                        IF(Ran2(iSeed).le.0.5) THEN
                            NotSelectedPair=.false.
                            OrbAlpha=Orb1+1
                            OrbBeta=Orb2
                        ELSE
                            NotSelectedPair=.false.
                            OrbAlpha=Orb2+1
                            OrbBeta=Orb1
                        ENDIF
                    ENDIF

                enddo

             ENDIF

!             WRITE(6,*) "Chosen orbs: ",iSpn,OrbAlpha,OrbBeta

             IF(iSpn.eq.1) THEN
!We want a beta beta pair.
                 IF(Orb1Beta.ne.Orb2Beta) THEN
!The two states are different - we now just need to see if one of the orbitals is in NI
!The checks will return true if the orbital is in NI.
                     Check1=BTEST(ILUT((Orb1Beta-1)/32),
     &                   MOD(Orb1Beta-1,32))
                     Check2=BTEST(ILUT((Orb2Beta-1)/32),
     &                   MOD(Orb2Beta-1,32))
                     IF(.not.(Check1.or.Check2)) THEN
!We have found an acceptable pair of virtual orbitals. Excite into these.
                         FinalA=Orb1Beta
                         FinalB=Orb2Beta
                         EXIT
                     ENDIF
                 ENDIF

             ELSEIF(iSpn.eq.2) THEN
!We want an alpha beta pair
                Check1=BTEST(ILUT((OrbAlpha-1)/32),
     &              MOD(OrbAlpha-1,32))
                Check2=BTEST(ILUT((OrbBeta-1)/32),
     &              MOD(OrbBeta-1,32))
                IF(.not.(Check1.or.Check2)) THEN
!We have found an acceptable pair of virtual orbitals. Excite into these.
                     FinalA=OrbAlpha
                     FinalB=OrbBeta
                     EXIT
                 ENDIF
             ELSEIF(iSpn.eq.3) THEN
!We want a alpha alpha pair
                 IF(Orb1Alpha.ne.Orb2Alpha) THEN
!The two states are different - we now just need to see if one of the orbitals is in NI
!The checks will return true if the orbital is in NI.
                     Check1=BTEST(ILUT((Orb1Alpha-1)/32),
     &                   MOD(Orb1Alpha-1,32))
                     Check2=BTEST(ILUT((Orb2Alpha-1)/32),
     &                   MOD(Orb2Alpha-1,32))
                     IF(.not.(Check1.or.Check2)) THEN
!We have found an acceptable pair of virtual orbitals. Excite into these.
                         FinalA=Orb1Alpha
                         FinalB=Orb2Alpha
                         EXIT
                     ENDIF
                 ENDIF

             ENDIF

         enddo  !End over loop until we find an allowed virtual pair

!         WRITE(6,*) "OrbPair to vacate: ", 
!     &            OrbPairs(1,I+SymProdInd(1,iSpn,iFrom)),
!     &            OrbPairs(2,I+SymProdInd(1,iSpn,iFrom))
!         WRITE(6,*) "Orbs to fill: ",
!     &           iSpn,OrbAlpha,OrbBeta,ILUT(0)

         NK(1:NEl)=NI(1:NEL)
!ExcitMat(1,*) are the indices in NI to excite from and 2,* are the spin orbs to excite to.
         DO J=1,NEL
             IF(NI(J).EQ.OrbPairs(1,I+SymProdInd(1,iSpn,iFrom))) THEN
               ExcitMat(1,1)=J
               ExcitMat(2,1)=FinalA
            ELSEIF(NI(J).EQ.OrbPairs(2,I+SymProdInd(1,iSpn,iFrom))) THEN
               ExcitMat(1,2)=J
               ExcitMat(2,2)=FinalB
            ENDIF
         ENDDO
         CALL FindExcitDet(ExcitMat,NK,2,Parity)
!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
!         Excitlevel=iGetExcitLevel(NI,NK,NEL)
!         IF(Excitlevel.ne.2) THEN
!             WRITE(6,*) "Have not created a double excitation"
!            CALL WRITEDET(6,NK,NEL,.TRUE.)
!            CALL WRITEDET(6,NI,NEL,.TRUE.)
!            WRITE(6,*) I
!            STOP "Have not created a double excitation"
!         ENDIF
!         IF(.NOT.ISVALIDDET(NK,NEL)) THEN
!             WRITE(6,*) "INVALID DET SYMGENRANDEXCITIT DOUBLE"
!             CALL WRITEDET(6,NK,NEL,.TRUE.)
!             CALL WRITEDET(6,NI,NEL,.TRUE.)
!             STOP "INVALID DET SYMGENRANDEXCITIT DOUBLE"
!         ENDIF

      End Subroutine SymGenRandExcitIt2_Doub_GT_NoW

      Subroutine SymGenRandExcitIt2_Double_GetTo(iSeed,nI,ExcitTypes,
     &      iExcit,OrbPairs, SymProdInd,iSpn,iFrom,iLUT,
     &      I,nK, pGen,ExcitMat,Parity)
         USE symexcit2
         use SystemData, only: G1, nBasis, nBasisMax
         use SystemData, only: Arr,nEl
         use IntegralsData, only: UMat
         use constants, only: dp
         IMPLICIT NONE
         INTEGER iSeed, parity
         INTEGER nI(nEl)
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit,ExcitMat(2,2)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,1:*)
         INTEGER iSpn,iFrom
         INTEGER iLUT(0:*)
         INTEGER I
         real(dp) pGen
         INTEGER NK(nEl) 
         TYPE(ExcitWeight),pointer :: ews(:)
         real(dp) rCount
         INTEGER iCount
         INTEGER K,J
         real(dp) RAN2
!  Allocate memory to store all weights and excitations in this category
         allocate(ews(EXCITTYPES(5,IEXCIT)/SYMPRODIND(2,ISPN,IFROM)))
!         WRITE(6,*)EXCITTYPES(5,IEXCIT)/SYMPRODIND(2,ISPN,IFROM),
!     &         EXCITTYPES(5,IEXCIT),SYMPRODIND(2,ISPN,IFROM)
!.. Now go through the list of virtual pairs, ng those with orbitals in NI 
!.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
!.. states (listed in SymStatePairs), and the index of the start of this list
!.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
!.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose prod is of that symmetry.
         RCOUNT=0.0_dp
         ICOUNT=0
!  L is a loop variable.  We go through the loop once to count the number of possible excitations
!  We then pick a random number up to that number, and choose that excitation
!  This is modified to weight each excitation according to its U-matrix element, so each
!  Excitation's weight is added to RCOUNT to give a normalization.  RCOUNT*RAND then corresponds to 
!  a point in the sum of weights which is found by subtracting weights until we reach 0.
!  This is only done for the excitations TO because we need to know the excitations FROM to generate the weight
         CALL EnumExcitWeights(EXCITTYPES(1,IEXCIT),I,ILUT,ews,OrbPairs,
     &      SymProdInd,RCOUNT,ICOUNT,NBASISMAX,ARR,NBASIS)

!         Write(6,*) "Doub To", iCount,rCount
         pGen=pGen/rCount           ! Normalize for the sum of TOs
         
!  We've now got a list of iCount weights in ews, and total norm rCount.
!  Select a random one
         rCount=rCount*RAN2(ISEED)
         K=0
!  This is a little inefficient, and could perhaps be binary searched if cumulative weights were also stored
         DO WHILE (rCount.GT.0.0_dp)
            K=K+1
            rCount=rCount-ews(K)%Weight
         ENDDO 
!         WRITE(6,*) K,rCount  
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
!         Call WriteDet(6,nK,nEl,.true.)
!         WRITE(6,*) ExcitMat
         CALL FindExcitDet(ExcitMat,NK,2,Parity)
!         Call WriteDet(6,nK,nEl,.true.)
         deallocate(ews)
!         WRiTE(6,*) "ENd GetTo"
         RETURN
      End Subroutine SymGenRandExcitIt2_Double_GetTo
