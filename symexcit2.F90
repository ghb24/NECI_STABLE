MODULE SymExcit2
      
      USE Calc , only : G_VMC_EXCITWEIGHT,G_VMC_EXCITWEIGHTS,CUR_VERT,EXCITFUNCS
      USE Calc , only : TUPOWER
      Use Integrals, only: ChemPot
      IMPLICIT NONE

      TYPE ExcitWeight
! The orbitals excited from
         INTEGER I,J
! The orbitals excited to
         INTEGER A,B
         REAL*8 WEIGHT
      END TYPE ExcitWeight
! Size in terms of reals.
      integer, PARAMETER :: ExcitWeightSize=3
      CONTAINS 

!  Enumerate the weights of all possible determinants to excite from in a given excittype.
      SUBROUTINE EnumExcitFromWeights(ExcitType, ews,OrbPairs, SymProdInd,Norm,iCount,G1,nBasisMax,UMat,Arr,nBasis)
         USE HElem
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER ExcitType(5)
         INTEGER nBasis,nBasisMax(5,*)
         INTEGER OrbPairs(2,*)
         REAL*8 Norm
         INTEGER iCount
         TYPE(HElement) UMat(*)
         REAL*8 Arr(nBasis,2)
         TYPE(ExcitWeight) ews(*)
         TYPE(BasisFN) G1(*)
         INTEGER SymProdInd(2,3,0:*)
         INTEGER iSpn,iFrom,iFromIndex
         INTEGER I
!.. We store each excitation type as:
!.. 1   TYPE (single=1, double=2)
!.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
!.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
!.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
!.. 5  COUNT (Total number of excitations in this category)
         iSpn=EXCITTYPE(2)
         iFrom=EXCITTYPE(3)
!.. Go through the list of pairs with a given symprod.
!.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
!.. SYMPRODS(I) in ORBPAIRS
!.. SYMPRODIND(2,ISPN,I) contains the number of such elements
         iCount=0
         Norm=0.D0
         DO I=1,SYMPRODIND(2,iSpn,iFrom)
            iFromIndex=I+SymProdInd(1,iSpn,iFrom)
            CALL AddExcitFromWeight(                                    &
     &                  OrbPairs(1,iFromIndex),                         &
     &                  OrbPairs(2,iFromIndex),                         &
     &                  ews,Norm,iCount,G1,nBasisMax,UMat,Arr,nBasis)
         ENDDO
      END subroutine
!  Enumerate the excitations and weights of excitations of a given ExcitType.
      SUBROUTINE EnumExcitWeights(ExcitType,iFromIndex,iLUT,ews,OrbPairs,SymProdInd,Norm,iCount,G1,NBASISMAX,UMAT,Arr,NBASIS)
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         USE HElem
         use SymData, only: SymPairProds,SymStatePairs
         INTEGER ExcitType(5)
         INTEGER NBASIS
         INTEGER IFROM,ITO,ISPN
         INTEGER OrbPairs(2,*)
         INTEGER iLUT(0:*)
         INTEGER SymProdInd(2,3,0:*)
         INTEGER ICC1,ICC2,ICC3,ICC4
         INTEGER iCount
         INTEGER iFromIndex
         LOGICAL L1B,L1A,L2B,L2A
         INTEGER nBasisMax(5,*)
         TYPE(BasisFN) G1(*)
         TYPE(HElement) UMAT(*)
         REAL*8 Norm
         TYPE(ExcitWeight) ews(*)
         INTEGER K
         REAL*8 Arr(nBasis,2)
         ISPN=ExcitType(2)
         IFROM=ExcitType(3)
         ITO=ExcitType(4)
! K loops over the possible pairs of states which we are allowed to excite to.
         DO K=SymPairProds(ITO)%nIndex,SymPairProds(ITO)%nIndex+SymPairProds(ITO)%nPairs-1
!.. Now check according to ISPN
!.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
!.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
            ICC1=SymStatePairs(1,K)*2-1
            ICC2=ICC1+1
            ICC3=SymStatePairs(2,K)*2-1
            ICC4=ICC3+1
            L1B=BTEST(ILUT((ICC1-1)/32),MOD(ICC1-1,32))
            L1A=BTEST(ILUT((ICC2-1)/32),MOD(ICC2-1,32))
            L2B=BTEST(ILUT((ICC3-1)/32),MOD(ICC3-1,32))
            L2A=BTEST(ILUT((ICC4-1)/32),MOD(ICC4-1,32))
!.. L1B is set if the beta of the first virtual is in NI, i.e. is disallowed
            IF(ISPN.EQ.1) THEN
!.. If both virtuals aren't the samem and neither are in NI, then allow
               IF(ICC1.NE.ICC3.AND..NOT.(L1B.OR.L2B)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC1,                                           &
     &                  ICC3,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,Arr,NBASIS)
               ENDIF
            ELSEIF(ISPN.EQ.2) THEN
!.. If neither virtuals are in NI, then allow
               IF(.NOT.(L1B.OR.L2A)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC1,                                           &
     &                  ICC4,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,Arr,NBASIS)
               ENDIF
!.. If neither virtuals are in NI, and they're not the same(which would give
!.. us the same excitation as previously), then allow
               IF(.NOT.(L1A.OR.L2B).AND.ICC1.NE.ICC3) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC2,                                           &
     &                  ICC3,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,Arr,NBASIS)
               ENDIF
            ELSEIF(ISPN.EQ.3) THEN
!.. If both virtuals aren't the samem and neither are in NI, then allow
               IF(ICC1.NE.ICC3.AND..NOT.(L1A.OR.L2A)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC2,                                           &
     &                  ICC4,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,Arr,NBASIS)
               ENDIF
            ENDIF
         ENDDO
      END subroutine
! Add the weight of the excitation to the list in ExWeights
! I,J are from, K,L are to
      SUBROUTINE AddExcitWeight(I,J,A,B,ExWeights,Norm,iCount,G1,NBASISMAX,UMAT,Arr,NBASIS)
         USE HElem
         use System, only: BasisFN
         INTEGER I,J,A,B
         REAL*8 R,Norm
         INTEGER nBasisMax(5,*),NBASIS
         TYPE(BasisFN) G1(*)
         TYPE(HElement) UMAT(*)
         TYPE(ExcitWeight) ExWeights(iCount+1)
         INTEGER iCount
         REAL*8 Arr(nBasis,2)
         CALL ExcitWeighting(I,J,A,B,R,G1,NBASISMAX,UMAT,Arr,NBASIS)
         iCount=iCount+1
         ExWeights(iCount)%I=I
         ExWeights(iCount)%J=J
         ExWeights(iCount)%A=A
         ExWeights(iCount)%B=B
         ExWeights(iCount)%Weight=R
         Norm=Norm+R
      END subroutine
         
! Add the weight of the 'from' excitation to the list in ExWeights
! I,J are from
      SUBROUTINE AddExcitFromWeight(I,J,ExWeights,Norm,iCount,G1,nBasisMax,UMat,Arr,nBasis)
         USE HElem
         use System, only: BasisFN
         INTEGER I,J,A,B
         REAL*8 R,Norm
         INTEGER nBasisMax(5,*),nBasis
         TYPE(BasisFN) G1(*)
         TYPE(HElement) UMat(*)
         TYPE(ExcitWeight) ExWeights(iCount+1)
         INTEGER iCount
         REAL*8 Arr(nBasis,2)
         CALL ExcitFromWeighting(I,J,R,G1,nBasisMax,UMat,Arr,nBasis)
         iCount=iCount+1
         ExWeights(iCount)%I=I
         ExWeights(iCount)%J=J
         ExWeights(iCount)%Weight=R
         Norm=Norm+R
      END subroutine
!        A sub called to generate an unnormalised weight for an ij->?? excitation
!          We return a function of the energies of the orbitals, exp(-(ei+ej)/a)
      SUBROUTINE ExcitFromWeighting(I,J,Weight,G1,nBasisMax,UMat,Arr,nBasis)
         USE HElem
         use System, only: BasisFN
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),nBasis
         TYPE(BasisFN) G1(nBasis)
!  We fake ISS
         INTEGER I,J,K,L
         REAL*8 WEIGHT
         TYPE(HElement) UMAT(*),W
         REAL*8 Arr(nBasis,2)
         !No weighting
         IF(EXCITFUNCS(10)) THEN
            Weight=1.D0
         !Exponential weighting
         ELSEIF(EXCITFUNCS(1)) THEN
            Weight=EXP((Arr(I,2)+Arr(J,2))*g_VMC_ExcitWeights(1,CUR_VERT))
         !Polynomial weighting with a cut-off at the chemical potential
         ELSEIF(EXCITFUNCS(4)) THEN
            IF((Arr(I,2)+Arr(J,2)).GT.(2.D0*CHEMPOT)) Weight=1.D0
            IF((Arr(I,2)+Arr(J,2)).LE.(2.D0*CHEMPOT)) THEN
                Weight=(1.D0/((-(Arr(I,2)+Arr(J,2))+(2.D0*CHEMPOT)+1.D0)**g_VMC_ExcitWeights(1,CUR_VERT)))
            ENDIF
         !CHEMPOT-TWOFROM
         ELSEIF(EXCITFUNCS(5)) THEN
            IF((Arr(I,2)+Arr(J,2)).GT.(2.D0*CHEMPOT)) THEN
                Weight=(1.D0/(((Arr(I,2)+Arr(J,2))-(2.D0*CHEMPOT)+1.D0)**g_VMC_ExcitWeights(2,CUR_VERT)))
            ENDIF
            IF((Arr(I,2)+Arr(J,2)).LE.(2.D0*CHEMPOT)) THEN
                Weight=(1.D0/((-(Arr(I,2)+Arr(J,2))+(2.D0*CHEMPOT)+1.D0)**g_VMC_ExcitWeights(1,CUR_VERT)))
            ENDIF
         !PolyExcitWeighting (for both real & virtual orbs)
         ELSEIF(EXCITFUNCS(3)) THEN
            IF((Arr(I,2)+Arr(J,2)).GT.g_VMC_ExcitWeights(1,CUR_VERT)) Weight=1.D0
            IF((Arr(I,2)+Arr(J,2)).LE.g_VMC_ExcitWeights(1,CUR_VERT)) THEN
                Weight=(1.D0/((-(Arr(I,2)+Arr(J,2))+g_VMC_ExcitWeights(1,CUR_VERT)+1.D0)**g_VMC_ExcitWeights(2,CUR_VERT)))
            ENDIF
         ENDIF
!         write(83,"(4G25.16)") (Arr(I,2)+Arr(J,2)), g_VMC_ExcitWeights(1,CUR_VERT), g_VMC_ExcitWeights(2,CUR_VERT), Weight
         RETURN
      END subroutine
!        A sub called to generate an unnormalised weight for a given ij->kl excitation
!          We return a function of the U matrix element (|<ij|u|kl>|^2)^G_VMC_EXCITWEIGHT
      SUBROUTINE EXCITWEIGHTING(I,J,K,L,WEIGHT,G1,NBASISMAX,UMAT,Arr,NBASIS)
         USE HElem
         USE UMatCache , only : GTID
         USE Integrals, only : GetUMatEl
         use System, only: BasisFN
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),NBASIS
         TYPE(BasisFN) G1(NBASIS)
!  We fake ISS
         INTEGER ISS
         INTEGER IDI,IDJ,IDK,IDL
         INTEGER I,J,K,L,ISUB
         REAL*8 WEIGHT,W2
         TYPE(HElement) UMAT(*),W
         REAL*8 Arr(nBasis,2),Alat(3)
         IF(G_VMC_EXCITWEIGHT(CUR_VERT).EQ.0.D0) THEN
            WEIGHT=1.D0
         ELSE
!            CALL TISET('UMATELWT',ISUB)
            ISS=NBASISMAX(2,3)
            CALL GTID(NBASISMAX,I,IDI)
            CALL GTID(NBASISMAX,J,IDJ)
            CALL GTID(NBASISMAX,K,IDK)
            CALL GTID(NBASISMAX,L,IDL)
            W=GetUMatEl(NBASISMAX,UMAT,Alat,NBASIS,ISS,G1,IDI,IDJ,IDK,IDL)
            IF(TUPOWER) THEN
                WEIGHT=1.D0+(SQRT(SQ(W)))**(G_VMC_EXCITWEIGHT(CUR_VERT))
            ELSE
                WEIGHT=EXP(SQRT(SQ(W))*G_VMC_EXCITWEIGHT(CUR_VERT))
            ENDIF
!            CALL TIHALT('UMATELWT',ISUB)
         ENDIF
         IF(.not.EXCITFUNCS(10)) THEN
             IF((EXCITFUNCS(1)).and.(g_VMC_ExcitWeights(3,CUR_VERT).NE.0.D0)) THEN
                W2=ABS(((Arr(I,2)+Arr(J,2))-(Arr(K,2)+Arr(L,2))))
                IF(ABS(W2).LT.1.D-2) W2=1.D-2
                Weight=Weight*W2**g_VMC_ExcitWeights(3,CUR_VERT)
             ENDIF
             IF(EXCITFUNCS(1)) Weight=Weight*EXP(-(Arr(K,2)+Arr(L,2))*g_VMC_ExcitWeights(2,CUR_VERT))
!chempotweighting - using a chemical potential cut-off
             IF(EXCITFUNCS(4)) THEN
                 IF((Arr(K,2)+Arr(L,2)).LT.(CHEMPOT*2.D0)) Weight=Weight
                 IF((Arr(K,2)+Arr(L,2)).GE.(CHEMPOT*2.D0)) THEN
                    Weight=Weight*(1.D0/(((Arr(K,2)+Arr(L,2))-(2.D0*CHEMPOT)+1.D0)**g_VMC_ExcitWeights(2,CUR_VERT)))
                 ENDIF
             ENDIF
!CHEMPOT-TWOFROM
             IF(EXCITFUNCS(5)) THEN
                 IF((Arr(K,2)+Arr(L,2)).LT.(CHEMPOT*2.D0)) Weight=Weight
                 IF((Arr(K,2)+Arr(L,2)).GE.(CHEMPOT*2.D0)) THEN
                     Weight=Weight*(1.D0/(((Arr(K,2)+Arr(L,2))-(2.D0*CHEMPOT)+1.D0)**g_VMC_ExcitWeights(3,CUR_VERT)))
                 ENDIF
             ENDIF
!PolyExcitWeighting
             IF(EXCITFUNCS(2)) THEN
                 IF((Arr(K,2)+Arr(L,2)).LT.g_VMC_ExcitWeights(2,CUR_VERT)) Weight=Weight
                 IF((Arr(K,2)+Arr(L,2)).GE.g_VMC_ExcitWeights(2,CUR_VERT)) THEN
                     Weight=Weight*(1.D0/(((Arr(K,2)+Arr(L,2))-g_VMC_ExcitWeights(2,CUR_VERT)+1.D0)**g_VMC_ExcitWeights(3,CUR_VERT)))
                 ENDIF
             ENDIF
         ENDIF
         RETURN
      END subroutine

!We wish to calculate what excitation class the excitation NI->NJ falls into, with the appropriate 
! IFROM class and index within that, IFROMINDEX, and ITO class, and index within that. ITOINDEX
! After that, we generate the probability that nJ would be an excitation from nI.

!  WARNING - this currently only works for abelian symmetry groups
      SUBROUTINE GenExcitProbInternal(nI,nJ,nEl,G1,nBasisMax,UMat,Arr,nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,pGen)
         USE HElem
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nSymPairProds,SymPairProds
         IMPLICIT NONE
         INTEGER iExcit(2,2)
         LOGICAL L
         INTEGER nI(nEl),nJ(nEl),nEl,nBasis
         INTEGER iFrom,iFromIndex,iTo,iToIndex
         TYPE(BasisFn) G1(nBasis)
         TYPE(Symmetry) Sym,Prod
         TYPE(Symmetry) SYMPROD
         TYPE(Symmetry) SymProds(0:*)
         LOGICAL SYMEQ
         INTEGER nBasisMax(5,*)
         TYPE(HElement) UMat(*)
         TYPE(ExcitWeight) ews(*)
         POINTER(IP_ews,ews)
         INTEGER iLUT(*)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,0:*)
         INTEGER ExcitTypes(5,*)
         REAL*8 Norm
         INTEGER K,NPR,I,iExcitType
         INTEGER iSpn,iCount,nToPairs,iTotal,nFromPairs
         REAL*8 pGen
         REAL*8 Arr(nBasis,2)
         LOGICAL IsUHFDet
         iExcit(1,1)=2
         CALL GETEXCITATION(nI,nJ,nEl,iExcit,L)
!EXCIT(1,*) are the ij... in NI, and EXCIT(2,*) the ab... in NJ
         IF(iExcit(1,1).LE.0) THEN
! More than a double excitation, so not allowed
            IFROM=0
            pGen=0
            RETURN
         ENDIF
!  See if we have a single
         IF(iExcit(1,2).EQ.0) THEN
            if(isUHFDet(nI,nEl)) then
               pGen=0.D0      !HF -> single has 0 prob
            ELSE
!  The UHF det isn't set if Brillouin's Theorem is disabled, and we end up here.
!  NB we can still generate the HF det from a single excitation, and that would end up here.
               pGen=1.D0/iTotal
            ENDIF
            RETURN
         ENDIF
!  We have a double
! Now check the sym prod of i and j to work out the FROM
         Sym=SYMPROD( G1(iExcit(1,1))%Sym,G1(iExcit(1,2))%Sym)
         iFrom=1
         DO WHILE(.NOT.SYMEQ(SymProds(iFrom),Sym))
            iFrom=iFrom+1
         ENDDO
!.. Now to find the appropriate iTo index:
         Sym=SYMPROD( G1(iExcit(2,1))%Sym, G1(iExcit(2,2))%Sym)
!..  Find which SymPairProd it corresponds to 
!.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
!.. symlabels (listed in SymStatePairs), and the index of the start of this list
!.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
!.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of symlabels whose prod is of that symmetry.
         CALL FindSymProd(Sym,SymPairProds,nSymPairProds,iTo)
!.. There are three values of ISPN.  ISPN=1 is beta beta, ISPN=2 is alpha, beta and beta, alpha.  ISPN=3 is alpha alpha
         iSpn=(G1(iExcit(1,1))%Ms+G1(iExcit(1,2))%Ms)/2+2
!.. Now find which excittype we correspons to 
         K=1
         DO WHILE(ExcitTypes(1,K).NE.2.OR.ExcitTypes(2,K).NE.iSpn.OR.ExcitTypes(3,K).NE.iFrom.OR.ExcitTypes(4,K).NE.iTo)
            K=K+1
         ENDDO
!.. K is now the excitation
         iExcitType=K
         pGen=(ExcitTypes(5,iExcitType)+0.D0)/iTotal  ! pGen is the prob of choosing excit iExcitType


!.. We've worked out what class the IFROM was.  Now work out which member of the class it is
!  Now enumerate all the FROMs, and their weightings
         nFromPairs=SymProdInd(2,iSpn,iFrom)
         CALL MEMORY(IP_ews,ExcitWeightSize*nFromPairs,'ExcitWGEPI')
         iCount=0
         Norm=0.D0
         CALL EnumExcitFromWeights(ExcitTypes(1,iExcitType),ews,OrbPairs,SymProdInd,Norm,iCount,G1,nBasisMax,UMat,Arr,nBasis)
         pGen=pGen/Norm    ! divide through by the Norm of the FROMs
         DO I=1,nFromPairs
            IF(ews(I)%I.EQ.iExcit(1,1).AND.ews(I)%J.EQ.iExcit(1,2)) THEN
               pGen=pGen*ews(I)%Weight !Multiply by the weight of the one chosen
               EXIT
            ENDIF
         ENDDO
         iFromIndex=I
         CALL FREEM(IP_ews)
!  pGen is the prob of choosing a specific FROM (given having chosen iExcitType proportional to the number of excitations in each iExcitType)
!           times the prob of choosing iExcitType

!.. Now work out the index of the (a,b) pair within the prod
         nToPairs=ExcitTypes(5,iExcitType)/SymProdInd(2,iSpn,iFrom)
         CALL MEMORY(IP_ews,ExcitWeightSize*nToPairs,'ExcitWGEPI')
         iCount=0
         Norm=0.D0
         CALL EnumExcitWeights(ExcitTypes(1,iExcitType),iFromIndex,iLUT,ews,OrbPairs,SymProdInd,Norm,iCount,G1,nBasisMax,UMat,Arr,nBasis)
!.. Find the (a,b) pair
!.. The prob of all possible excitations in this iTo 
         DO I=1,nToPairs
            IF(ews(I)%A.EQ.iExcit(2,1).AND.ews(I)%B.EQ.iExcit(2,2)) THEN
               pGen=pGen*ews(I)%Weight
               EXIT
            ENDIF
         ENDDO
         pGen=pGen/Norm    ! The norm of the TOs
!  pGen is the prob of choosing a specific TO (given the FROM, and the iExcitType)
!           times prob of choosing a specific FROM (given having chosen iExcitType proportional to the number of excitations in each iExcitType)
!           times the prob of choosing iExcitType
         CALL FREEM(IP_ews)
      END subroutine

!We wish to calculate whether NJ is an excitation of NI.
!WARNING - this currently only works for abelian symmetry groups
      SUBROUTINE IsConnectedDetInternal(nI,nJ,nEl,G1,nBasisMax,nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,tIsConnectedDet)
         USE HElem
         use System, only: Symmetry,SymmetrySize,SymmetrySizeB
         use System, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER iExcit(2,2)
         LOGICAL L
         INTEGER nI(nEl),nJ(nEl),nEl,nBasis
         INTEGER iFrom,iFromIndex,iTo,iToIndex
         TYPE(BasisFn) G1(nBasis)
         TYPE(Symmetry) Sym,Prod
         TYPE(Symmetry) SYMPROD
         TYPE(Symmetry) SymProds(0:*)
         LOGICAL SYMEQ
         INTEGER nBasisMax(5,*)
         TYPE(ExcitWeight) ews(*)
         POINTER(IP_ews,ews)
         INTEGER iLUT(*)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,0:*)
         INTEGER ExcitTypes(5,*)
         REAL*8 Norm
         INTEGER K,NPR,I,iExcitType
         INTEGER iSpn,iCount,nToPairs,iTotal,nFromPairs
         LOGICAL tIsConnectedDet
         LOGICAL IsUHFDet
         iExcit(1,1)=2
         CALL GETEXCITATION(nI,nJ,nEl,iExcit,L)
!EXCIT(1,*) are the ij... in NI, and EXCIT(2,*) the ab... in NJ
         IF(iExcit(1,1).LE.0) THEN
! More than a double excitation, so not allowed
            tIsConnectedDet=.FALSE.
            RETURN
         ENDIF
!  See if we have a single
         IF(iExcit(1,2).EQ.0) THEN
!Warning - this is not necessarily the case, but will do for abelian symmetry groups
            IF(.NOT.(ISUHFDET(NI,NEL).OR.ISUHFDET(NJ,NEL))) THEN
!  The UHF det isn't set if Brillouin's Theorem is disabled, and we end up here.
               tIsConnectedDet=.TRUE.
            ELSE
               tIsConnectedDet=.FALSE.
            ENDIF
            RETURN
         ENDIF
         tIsConnectedDet=.TRUE.
         RETURN
      END subroutine
END MODULE SymExcit2
