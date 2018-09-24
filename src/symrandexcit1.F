
!== SymGenRandExcitIt2
!== Based on SymGenRandExcitIt, but modular, and capable of acceleration
!== for abelian symmetry groups, where lists of symmetry pairs needn't be
!== stored.

      SUBROUTINE SymGenRandExcitIt2(NI,NEL,EXCITTYPES,
     &               CLASSES,
     &               SYMPRODIND,ILUT,ORBPAIRS,iTotal,ISEED,
     &               NK,IC,PGEN,EXCITMAT,TPARITY)      
         use CalcData, only: ExcitFuncs
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         USE symexcit2
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: nBasis
         use SymData, only: SymClass,SymLabelCounts,SymLabelList
         use SymData, only: tAbelianFastExcitgen
         use SymData, only: tStoreStateList
         use SystemData, only: tNoRenormRandExcits
         use constants, only: dp,sizeof_int, maxExcit
         IMPLICIT NONE
         INTEGER nEl,nI(nEl)
         INTEGER I
         INTEGER iSpn,iFrom,iTo
         INTEGER iC
         real(dp) pGen
         INTEGER ILUT(0:*),ExcitMat(2,maxExcit)
         LOGICAL TPARITY
         TYPE(SymClass) CLASSES(*)
         INTEGER ExcitTypes(5,*)
         INTEGER iExcit
         INTEGER OrbPairs(2,*)
         INTEGER iTotal
         INTEGER SymProdInd(2,3,1:*)
         INTEGER nK(nEl)
         INTEGER iSeed
         real(dp) RAN2   
!.. irst generate a random excitation 

         PGEN=1.0_dp/ITOTAL
!         I=ITOTAL
!         IEXCIT=1
!         DO WHILE(I.GT.0.AND.EXCITTYPES(5,IEXCIT).LE.I)
!            I=I-EXCITTYPES(5,IEXCIT)
!            WRITE(6,"(A,7I5)") "SE",IEXCIT,
!     &                   EXCITTYPES(1,IEXCIT),
!     &                   EXCITTYPES(3,IEXCIT),
!     &                   EXCITTYPES(4,IEXCIT),
!     &                   EXCITTYPES(2,IEXCIT),
!     &                   EXCITTYPES(5,IEXCIT),I
!            IEXCIT=IEXCIT+1
!         ENDDO
         I=int(real(ITOTAL,dp)*RAN2(ISEED),sizeof_int)
!.. Now work out which Excitation Type it's in.
         IEXCIT=1
         DO WHILE(I.GT.0.AND.EXCITTYPES(5,IEXCIT).LE.I)
            I=I-EXCITTYPES(5,IEXCIT)
            IEXCIT=IEXCIT+1
         ENDDO
!         WRITE(6,*) "Excit:",iExcit,EXCITTYPES(1,IEXCIT)
!..   iExcit is the Excitation Type
         IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
!.. It's a single
            ISPN=EXCITTYPES(2,IEXCIT)-2
            IFROM=CLASSES(EXCITTYPES(3,IEXCIT))%SymLab
            ITO=EXCITTYPES(4,IEXCIT)
            IF(tNoRenormRandExcits) THEN
!There is no excitation weighting - we do not need to enumerate all excitations
                Call SymGenRandExcitIt2_Sing_GF_NoW(iSeed,
     &                  nEl,iSpn,iFrom,I,nI)
                Call SymGenRandExcitIt2_Sing_GT_NoW(iSeed,nI,
     &                  nEl,I,iSpn,iTo,nK,ExcitMat,TParity)
            ELSE
!We have an excitation weighting - enumerate all excitations to get normalisation
                Call SymGenRandExcitIt2_Sing_GF(iSeed,nI,nEl,
     &                  iSpn,iFrom,I)
                Call SymGenRandExcitIt2_Single_GetTo(iSeed,nI,nEl,I,
     &                  iSpn,iTo,nK,ExcitMat,TParity)
            ENDIF
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
            IF(tStoreStateList.and.tNoRenormRandExcits) THEN
!There is no excitation weighting - we do not need to enumerate all excitations in the class
!However, currently, we do need a list of state pairs.
                Call SymGenRandExcitIt2_Doub_GF_NoW(iSeed,
     &               iSpn,iFrom,SymProdInd,I,nI,OrbPairs)

                Call SymGenRandExcitIt2_Doub_GT_NoW(iSeed,
     &               nI,ExcitTypes,iExcit,OrbPairs,SymProdInd,iSpn,
     &               iFrom,iLUT,I,nK,ExcitMat,TParity)
            ELSE
                Call SymGenRandExcitIt2_Doub_GF(iSeed,ExcitTypes,
     &         iExcit,iSpn, iFrom, OrbPairs, SymProdInd, iTotal, I,pGen)

               Call SymGenRandExcitIt2_Double_GetTo(iSeed,nI,ExcitTypes,
     &          iExcit,OrbPairs, SymProdInd,iSpn,iFrom,iLUT, I,nK,pGen,
     &          ExcitMat,TParity)
            ENDIF
!.. Indicate we're a double.
            iC=2
         ENDIF
      End Subroutine
