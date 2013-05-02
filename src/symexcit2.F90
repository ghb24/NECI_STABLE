MODULE SymExcit2
      
      use CalcData , only : G_VMC_EXCITWEIGHT,G_VMC_EXCITWEIGHTS,CUR_VERT,EXCITFUNCS
      use CalcData , only : TUPOWER
      use IntegralsData, only: ChemPot
      use constants, only: dp
      IMPLICIT NONE

      TYPE ExcitWeight
          SEQUENCE
! The orbitals excited from
         INTEGER I,J
! The orbitals excited to
         INTEGER A,B
         real(dp) WEIGHT
      END TYPE ExcitWeight
! Size in terms of reals.
      integer, PARAMETER :: ExcitWeightSize=3
      CONTAINS 

!  Enumerate the weights of all possible determinants to excite from in a given excittype.
      SUBROUTINE EnumExcitFromWeights(ExcitType, ews,OrbPairs, SymProdInd,Norm,iCount,Arr,nBasis)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER ExcitType(5)
         INTEGER nBasis
         INTEGER OrbPairs(2,*)
         real(dp) Norm
         INTEGER iCount
         real(dp) Arr(nBasis,2)
         TYPE(ExcitWeight) ews(*)
         INTEGER SymProdInd(2,3,1:*)
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
         Norm=0.0_dp
         DO I=1,SYMPRODIND(2,iSpn,iFrom)
            iFromIndex=I+SymProdInd(1,iSpn,iFrom)
            CALL AddExcitFromWeight(                                    &
                       OrbPairs(1,iFromIndex),                         &
                       OrbPairs(2,iFromIndex),                         &
                       ews,Norm,iCount,Arr,nBasis)
         ENDDO
      END subroutine
!  Enumerate the excitations and weights of excitations of a given ExcitType.  
      SUBROUTINE EnumExcitWeights(ExcitType,iFromIndex,iLUT,ews,OrbPairs,SymProdInd,Norm,iCount,NBASISMAX,Arr,NBASIS)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymPairProds,SymStatePairs
         INTEGER ExcitType(5)
         INTEGER NBASIS
         INTEGER IFROM,ITO,ISPN
         INTEGER OrbPairs(2,*)
         INTEGER iLUT(0:*)
         INTEGER SymProdInd(2,3,1:*)
         INTEGER ICC1,ICC2,ICC3,ICC4
         INTEGER iCount
         INTEGER iFromIndex
         LOGICAL L1B,L1A,L2B,L2A
         INTEGER nBasisMax(5,*)
         real(dp) Norm
         TYPE(ExcitWeight) ews(*)
         INTEGER K
         real(dp) Arr(nBasis,2)
         INTEGER iLooped,iTo1,iTo2
         LOGICAL tDebugPrint
         tDebugPrint=.false.
         ISPN=ExcitType(2)
         IFROM=ExcitType(3)
         ITO=ExcitType(4)
! K loops over the possible pairs of states which we are allowed to excite to.
         K=-2
         Call SymGenExcitIt_GetNextPair(K,iTo,iLooped,iTo1,iTo2,  &
     &         tDebugPrint)
         DO WHILE(iLooped.lt.1)
!         DO K=SymPairProds(ITO)%nIndex,SymPairProds(ITO)%nIndex+SymPairProds(ITO)%nPairs-1
!.. Now check according to ISPN
!.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
!.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
            ICC1=iTo1
            ICC2=iTo1+1
            ICC3=iTo2
            ICC4=iTo2+1
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
     &                  ews,Norm,ICOUNT,NBASISMAX,Arr,NBASIS)
               ENDIF
            ELSEIF(ISPN.EQ.2) THEN
!.. If neither virtuals are in NI, then allow
               IF(.NOT.(L1B.OR.L2A)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC1,                                           &
     &                  ICC4,                                           &
     &                  ews,Norm,ICOUNT,NBASISMAX,Arr,NBASIS)
               ENDIF
!.. If neither virtuals are in NI, and they're not the same(which would give
!.. us the same excitation as previously), then allow
               IF(.NOT.(L1A.OR.L2B).AND.ICC1.NE.ICC3) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC2,                                           &
     &                  ICC3,                                           &
     &                  ews,Norm,ICOUNT,NBASISMAX,Arr,NBASIS)
               ENDIF
            ELSEIF(ISPN.EQ.3) THEN
!.. If both virtuals aren't the samem and neither are in NI, then allow
               IF(ICC1.NE.ICC3.AND..NOT.(L1A.OR.L2A)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC2,                                           &
     &                  ICC4,                                           &
     &                  ews,Norm,ICOUNT,NBASISMAX,Arr,NBASIS)
               ENDIF
            ENDIF
         Call SymGenExcitIt_GetNextPair(K,iTo,iLooped,iTo1,iTo2,  &
     &         tDebugPrint)
         ENDDO
      END subroutine
! Add the weight of the excitation to the list in ExWeights
! I,J are from, K,L are to
      SUBROUTINE AddExcitWeight(I,J,A,B,ExWeights,Norm,iCount,NBASISMAX,Arr,NBASIS)
         use SystemData, only: BasisFN
         INTEGER I,J,A,B
         real(dp) R,Norm
         INTEGER nBasisMax(5,*),NBASIS
         INTEGER iCount
         TYPE(ExcitWeight) ExWeights(iCount+1)
         real(dp) Arr(nBasis,2)
         CALL ExcitWeighting(I,J,A,B,R,NBASISMAX,Arr,NBASIS)
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
      SUBROUTINE AddExcitFromWeight(I,J,ExWeights,Norm,iCount,Arr,nBasis)
         use SystemData, only: BasisFN
         INTEGER I,J
         real(dp) R,Norm
         INTEGER nBasis
         INTEGER iCount
         TYPE(ExcitWeight) ExWeights(iCount+1)
         real(dp) Arr(nBasis,2)
         CALL ExcitFromWeighting(I,J,R,Arr,nBasis)
         iCount=iCount+1
         ExWeights(iCount)%I=I
         ExWeights(iCount)%J=J
         ExWeights(iCount)%Weight=R
         Norm=Norm+R
      END subroutine

!        A sub called to generate an unnormalised weight for an ij->?? excitation
!          We return a function of the energies of the orbitals, exp(-(ei+ej)/a)
      SUBROUTINE ExcitFromWeighting(I,J,Weight,Arr,nBasis)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER nBasis
!  We fake ISS
         INTEGER I,J
         real(dp) WEIGHT
         real(dp) Arr(nBasis,2)
         !No weighting
         IF(EXCITFUNCS(10)) THEN
            Weight=1.0_dp
         !Exponential weighting
         ELSEIF(EXCITFUNCS(1)) THEN
            Weight=EXP((Arr(I,2)+Arr(J,2))*g_VMC_ExcitWeights(1,CUR_VERT))
         !Polynomial weighting with a cut-off at the chemical potential
         ELSEIF(EXCITFUNCS(6)) THEN
!Step-function weighting
!First deal with electron I
            IF(Arr(I,2).lt.CHEMPOT) THEN
                Weight=g_VMC_ExcitWeights(1,CUR_VERT)
            ELSE
                Weight=1.0_dp
            ENDIF
!Then J...
            IF(Arr(J,2).lt.CHEMPOT) THEN
                Weight=Weight+g_VMC_ExcitWeights(1,CUR_VERT)
            ELSE
                Weight=Weight+1.0_dp
            ENDIF
         ELSEIF(EXCITFUNCS(4)) THEN
            IF((Arr(I,2)+Arr(J,2)).GT.(2.0_dp*CHEMPOT)) Weight=1.0_dp
            IF((Arr(I,2)+Arr(J,2)).LE.(2.0_dp*CHEMPOT)) THEN
                Weight=(1.0_dp/((-(Arr(I,2)+Arr(J,2))+(2.0_dp*CHEMPOT)+1.0_dp)**g_VMC_ExcitWeights(1,CUR_VERT)))
            ENDIF
         !CHEMPOT-TWOFROM
         ELSEIF(EXCITFUNCS(5)) THEN
            IF((Arr(I,2)+Arr(J,2)).GT.(2.0_dp*CHEMPOT)) THEN
                Weight=(1.0_dp/(((Arr(I,2)+Arr(J,2))-(2.0_dp*CHEMPOT)+1.0_dp)**g_VMC_ExcitWeights(2,CUR_VERT)))
            ENDIF
            IF((Arr(I,2)+Arr(J,2)).LE.(2.0_dp*CHEMPOT)) THEN
                Weight=(1.0_dp/((-(Arr(I,2)+Arr(J,2))+(2.0_dp*CHEMPOT)+1.0_dp)**g_VMC_ExcitWeights(1,CUR_VERT)))
            ENDIF
         !PolyExcitWeighting (for both real & virtual orbs)
         ELSEIF(EXCITFUNCS(3)) THEN
            IF((Arr(I,2)+Arr(J,2)).GT.g_VMC_ExcitWeights(1,CUR_VERT)) Weight=1.0_dp
            IF((Arr(I,2)+Arr(J,2)).LE.g_VMC_ExcitWeights(1,CUR_VERT)) THEN
                Weight=(1.0_dp/((-(Arr(I,2)+Arr(J,2))+g_VMC_ExcitWeights(1,CUR_VERT)+1.0_dp)**g_VMC_ExcitWeights(2,CUR_VERT)))
            ENDIF
         ENDIF
!         write(83,"(4G25.16)") (Arr(I,2)+Arr(J,2)), g_VMC_ExcitWeights(1,CUR_VERT), g_VMC_ExcitWeights(2,CUR_VERT), Weight
         RETURN
      END subroutine
!        A sub called to generate an unnormalised weight for a given ij->kl excitation
!          We return a function of the U matrix element (|<ij|u|kl>|^2)^G_VMC_EXCITWEIGHT
      SUBROUTINE EXCITWEIGHTING(I,J,K,L,WEIGHT,NBASISMAX,Arr,NBASIS)
         USE UMatCache , only : GTID
         use Integrals_neci, only : GetUMatEl
         use SystemData, only: BasisFN
         use global_utilities
         IMPLICIT NONE
         INTEGER nBasisMax(5,*),NBASIS
!  We fake ISS
         INTEGER ISS
         INTEGER IDI,IDJ,IDK,IDL
         INTEGER I,J,K,L
         !type(timer), save :: proc_timer
         real(dp) WEIGHT,W2
         HElement_t W
         real(dp) Arr(nBasis,2)
         IF(G_VMC_EXCITWEIGHT(CUR_VERT).EQ.0.0_dp) THEN
            WEIGHT=1.0_dp
         ELSE
!            proc_timer%timer_name='UMATELWT'
!            call set_timer(proc_timer)
            ISS=NBASISMAX(2,3)
            IDI = GTID(I)
            IDJ = GTID(J)
            IDK = GTID(K)
            IDL = GTID(L)
            W=GetUMatEl(IDI,IDJ,IDK,IDL)
            IF(TUPOWER) THEN
                WEIGHT=1.0_dp+(abs(W))**(G_VMC_EXCITWEIGHT(CUR_VERT))
            ELSE
                WEIGHT=EXP(abs(W)*G_VMC_EXCITWEIGHT(CUR_VERT))
            ENDIF
!            call halt_timer(proc_timer)
         ENDIF
         IF(.not.EXCITFUNCS(10)) THEN
             IF((EXCITFUNCS(1)).and.(g_VMC_ExcitWeights(3,CUR_VERT).NE.0.0_dp)) THEN
                W2=ABS(((Arr(I,2)+Arr(J,2))-(Arr(K,2)+Arr(L,2))))
                IF(ABS(W2).LT.1.0e-2_dp) W2=1.0e-2_dp
                Weight=Weight*W2**g_VMC_ExcitWeights(3,CUR_VERT)
             ELSEIF(EXCITFUNCS(1)) THEN
                Weight=Weight*EXP(-(Arr(K,2)+Arr(L,2))*g_VMC_ExcitWeights(2,CUR_VERT))
             ELSEIF(EXCITFUNCS(6)) THEN
!Step-function weighting using chemical potential cut-off
                 IF(Arr(K,2).gt.CHEMPOT) THEN
                     W2=g_VMC_ExcitWeights(2,CUR_VERT)
                 ELSE
                     W2=1.0_dp
                 ENDIF
                 IF(Arr(L,2).gt.CHEMPOT) THEN
                     W2=W2+g_VMC_ExcitWeights(2,CUR_VERT)
                 ELSE
                     W2=W2+1.0_dp
                 ENDIF
                 Weight=Weight*W2
!chempotweighting - using a chemical potential cut-off
             ELSEIF(EXCITFUNCS(4)) THEN
                 IF((Arr(K,2)+Arr(L,2)).LT.(CHEMPOT*2.0_dp)) Weight=Weight
                 IF((Arr(K,2)+Arr(L,2)).GE.(CHEMPOT*2.0_dp)) THEN
                    Weight=Weight*(1.0_dp/(((Arr(K,2)+Arr(L,2))-(2.0_dp*CHEMPOT)+1.0_dp)**g_VMC_ExcitWeights(2,CUR_VERT)))
                 ENDIF
!CHEMPOT-TWOFROM
             ELSEIF(EXCITFUNCS(5)) THEN
                 IF((Arr(K,2)+Arr(L,2)).LT.(CHEMPOT*2.0_dp)) Weight=Weight
                 IF((Arr(K,2)+Arr(L,2)).GE.(CHEMPOT*2.0_dp)) THEN
                     Weight=Weight*(1.0_dp/(((Arr(K,2)+Arr(L,2))-(2.0_dp*CHEMPOT)+1.0_dp)**g_VMC_ExcitWeights(3,CUR_VERT)))
                 ENDIF
!PolyExcitWeighting
             ELSEIF(EXCITFUNCS(2)) THEN
                 IF((Arr(K,2)+Arr(L,2)).LT.g_VMC_ExcitWeights(2,CUR_VERT)) Weight=Weight
                 IF((Arr(K,2)+Arr(L,2)).GE.g_VMC_ExcitWeights(2,CUR_VERT)) THEN
                     Weight=Weight*(1.0_dp/(((Arr(K,2)+Arr(L,2))-                         &
                     & g_VMC_ExcitWeights(2,CUR_VERT)+1.0_dp)**g_VMC_ExcitWeights(3,CUR_VERT)))
                 ENDIF
             ENDIF
         ENDIF
         RETURN
      END subroutine

!We wish to calculate what excitation class the excitation NI->NJ falls into, with the appropriate 
! IFROM class and index within that, IFROMINDEX, and ITO class, and index within that. ITOINDEX
! After that, we generate the probability that nJ would be an excitation from nI.

!  WARNING - this currently only works for abelian symmetry groups
      SUBROUTINE GenExcitProbInternal(nI,nJ,nEl,G1,nBasisMax,Arr,nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,pGen)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: nSymPairProds,SymPairProds
         use MemoryManager, only: TagIntType
         use sym_mod
         use global_utilities
         IMPLICIT NONE
         INTEGER iExcit(2,2)
         LOGICAL L
         INTEGER nEl,nI(nEl),nJ(nEl),nBasis
         INTEGER iFrom,iFromIndex,iTo
         TYPE(BasisFn) G1(nBasis)
         TYPE(Symmetry) Sym
         TYPE(Symmetry) SymProds(0:*)
         INTEGER nBasisMax(5,*)
         TYPE(ExcitWeight), allocatable :: ews(:)
         integer(TagIntType), save :: tagews=0
         INTEGER iLUT(*)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,1:*)
         INTEGER ExcitTypes(5,*)
         real(dp) Norm
         INTEGER K,I,iExcitType
         INTEGER iSpn,iCount,nToPairs,iTotal,nFromPairs
         real(dp) pGen
         real(dp) Arr(nBasis,2)
         LOGICAL IsUHFDet
         character(*), parameter :: thisroutine='GenExcitProbInternal'
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
               pGen=0.0_dp      !HF -> single has 0 prob
            ELSE
!  The UHF det isn't set if Brillouin's Theorem is disabled, and we end up here.
!  NB we can still generate the HF det from a single excitation, and that would end up here.
               pGen=1.0_dp/iTotal
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
         pGen=(ExcitTypes(5,iExcitType)+0.0_dp)/iTotal  ! pGen is the prob of choosing excit iExcitType


!.. We've worked out what class the IFROM was.  Now work out which member of the class it is
!  Now enumerate all the FROMs, and their weightings
         nFromPairs=SymProdInd(2,iSpn,iFrom)
         allocate(ews(nFromPairs))
         call LogMemAlloc('ExcitWGEPI',nFromPairs,8*ExcitWeightSize,thisroutine,tagEWS)
         iCount=0
         Norm=0.0_dp
         CALL EnumExcitFromWeights(ExcitTypes(1,iExcitType),ews,OrbPairs,SymProdInd,Norm,iCount,Arr,nBasis)
         pGen=pGen/Norm    ! divide through by the Norm of the FROMs
         DO I=1,nFromPairs
            IF(ews(I)%I.EQ.iExcit(1,1).AND.ews(I)%J.EQ.iExcit(1,2)) THEN
               pGen=pGen*ews(I)%Weight !Multiply by the weight of the one chosen
               EXIT
            ENDIF
         ENDDO
         iFromIndex=I
         deallocate(ews)
         call LogMemDealloc(thisroutine,tagEWS)
!  pGen is the prob of choosing a specific FROM (given having chosen iExcitType proportional 
!to the number of excitations in each iExcitType)
!           times the prob of choosing iExcitType

!.. Now work out the index of the (a,b) pair within the prod
         nToPairs=ExcitTypes(5,iExcitType)/SymProdInd(2,iSpn,iFrom)
         allocate(ews(nToPairs))
         call LogMemAlloc('ExcitWGEPI',nToPairs,8*ExcitWeightSize,thisroutine,tagEWS)
         iCount=0
         Norm=0.0_dp
         CALL EnumExcitWeights(ExcitTypes(1,iExcitType),iFromIndex,iLUT,ews,OrbPairs,SymProdInd,Norm,iCount,nBasisMax,Arr,nBasis)
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
!           times prob of choosing a specific FROM (given having chosen iExcitType 
!proportional to the number of excitations in each iExcitType)
!           times the prob of choosing iExcitType
         deallocate(ews)
         call LogMemDealloc(thisroutine,tagEWS)
      END subroutine

!We wish to calculate whether NJ is an excitation of NI.
!WARNING - this currently only works for abelian symmetry groups
      SUBROUTINE IsConnectedDetInternal(nI,nJ,tIsConnectedDet)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB, nel
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         IMPLICIT NONE
         INTEGER iExcit(2,2)
         LOGICAL L
         INTEGER nI(nEl),nJ(nEl)
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
!=---- 
!=SSE_CreateClassList(Classes, *ClassCount)
!=
!= Create Class List              O(N)                      4 N
!=   Gives NCL =O(NSYM)
!=
!= Create Class Remainder info    O(NCL)                    
!=----
!=SSE_CreateClassSymProds(nPr,nPairs,nCl,SymProds,SymProdCount,*ClassCount,Classes)
!=SSEA_CreateClassSymProds(nPr,nCl,SymProds,SymProdCount,*ClassCount,Classes)  !  Better scaling
!=
!= Create Class Sym Product info  O(NCL^2)                  4 NPR 
!=   Gives NPR ~= NCL^2
!=    for abelian  NPR<=NSYM
!=
!=----
!=Form SymProdInd from SymProdCount
!=----
!=SSE_StoreOccPairs(OrbPairs,nPairs,nPr,SymProdInd,SymProds)
!= !Abelian - do nothing
!=
!= Save list of occ orbital pairs O(N^2) /NPROC             2 O(N^2) /NPROC
!=   This is sorted according
!=    to the symmetry product
!=                                                                   For abelian, this isn't
!=                                  [ O(NCL^2) ]          [ NSYM ]   needed as we can very easily
!=                                                                   determine this from a list of 
!=                                                                   occs sorted by sym.  We should
!=                                                                   store counts however
!=
!=----
!=SSE_CountVirtProds(nDoub,nExcitTypes,nPr,nSymPairProds,nAllowPPS,iLUT)
!=SSEA_CountVirtProds(nDoub,nExcitTypes,nPr,nSymPairProds,SymProdCount,nAllowPPS) ! Calculated not enumerated
!=
!= Count number of allowed virt
!=    pairs of each sym prod      O(NPR M^2 / NCL^2)        O(NCL^2)
!=    Also counts num of doubles
!=                                                                   For abelian this is simply
!=                                  [ NSYM ]              [ NSYM ]   the difference between the num
!=                                                                   of occ pairs of a given sym and
!=                                                                   the total num of pairs.
!=
!=----
!=SSE_CountSingles(nSing,nCl,nExcitTypes,*ClassCount,Classes)
!=
!= Calculate the number of
!=    single excitations          O(NCL^3)
!=                                  [ NSYM ]
!= --- Exit here for init
!= 
!=----
!=SSE_StoreSingles(nExcitTypes,nCl,Classes,ThisClassCount,ExcitTupes)
!=
!= Store singles classes          O(NCL^3)                  O(NCL^2)
!=                                  [ NSYM ]                 [ NSYM ]
!=----
!=SSE_StoreDoubles(nPr,nSymPairProds,nAllowPPS,ExcitTypes,nExcitTypes,SymProds,SymProdInd)
!= Store doubles classes          O(NPR NSYM)               O(NPR^2)
!=                                  [ NSYM^2 ]               [ NSYM^2 ] 
      subroutine symsetupexcits3_worker (nI, nel, g1, nbasis, store, &
                                 SPIin, ETin, NAPin, OPin, tCount, iCount, &
                                 classes, ilut, symprods, iLevel, iMinElec1, &
                                 iMaxElec1)
         use global_utilities
         use SystemData, only: Symmetry,BasisFN,tAssumeSizeExcitgen
         use SymData, only: SymClass,nSymPairProds,nSymLabels
         use SymData, only: tAbelianFastExcitGen
         use SymData, only: tStoreStateList
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         TYPE(BasisFN) G1(nBasis)
         INTEGER,pointer :: DSTORE(:)
         INTEGER STORE(6)
         INTEGER ICOUNT
         LOGICAL TCOUNT
         INTEGER ILEVEL
         INTEGER iMinElec1,iMaxElec1

         TYPE(SymClass) CLASSES(*)
         TYPE(Symmetry) SYMPRODS(0:NEL*NEL) !nel*nel is the max it could be
         INTEGER CLASSCOUNT(2,NEL)
         INTEGER THISCLASSCOUNT(2,NEL)
!  ThisClassCount is used to list only electrons which this processor deals with
         INTEGER PREVCLASSCOUNT(2,NEL)
!  PrevClassCount is used to list electrons which lower indexed processors deal with
         INTEGER SYMPRODCOUNT(3,0:NEL*NEL)
         INTEGER, target :: SPIin(1:2,1:3,1:*)
         INTEGER,pointer :: SYMPRODIND(:,:,:)
         INTEGER ILUT(0:nBasis/32)
         INTEGER, target :: ETin(1:5,*)
         INTEGER,pointer :: EXCITTYPES(:,:)
         INTEGER nPairs, nSing, nDoub, nExcits
         INTEGER, target :: NAPin(1:3,*)
         INTEGER,pointer :: NALLOWPPS(:,:)
         INTEGER, target :: OPin(1:2,*)
         INTEGER,pointer :: ORBPAIRS(:,:)
         INTEGER nCl,nExcitTypes,nPr

         LOGICAL ISUHFDET
 
         INTEGER I         
         type(timer), save :: proc_timer
         proc_timer%timer_name='SYMSUEXCIT'
         
         call set_timer(proc_timer,65)
         Call SymSetupExcits_CreateClassList(nI,nEl,Classes, &
      iMinElec1, iMaxElec1, ThisClassCount, PrevClassCount,ClassCount, &
      G1, nCl)
         SYMPRODCOUNT(:,:)=0
         Call SymSetupExcits_CreateCSProds(nPr,nPairs,nCl, &
        SymProds, ThisClassCount, PrevClassCount, ClassCount,Classes, &
        SymProdCount)
!nPr now contains the number of items in SymProds
!.. Allocate enough memory to store the index
         IF(STORE(5).EQ.0) THEN
            allocate(SYMPRODIND(2,3,1:NPR))
         ELSE
            SYMPRODIND=>SPIin(1:2,1:3,1:NPR)
         ENDIF
         SYMPRODIND(1:2,1:3,1:NPR)=0

!.. Now shift this such that SYMPRODCOUNT(ISPN,I) is the index of
!.. the first orbital pair of SYMPROD(I) with spin ISPN in ORBPAIRS
!.. Store in SYMPRODIND(1,ISPN,I) too
         DO I=NPR,1,-1
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)
            SYMPRODIND(1,3,I)=SYMPRODCOUNT(3,I)
            SYMPRODIND(1,2,I)=SYMPRODCOUNT(2,I)
            SYMPRODIND(1,1,I)=SYMPRODCOUNT(1,I)
         ENDDO
!.. We allocate enough memory to store all the pairs.
!.. Each pair consists of (ORB1,ORB2) where ORB1<ORB2
         IF(STORE(4).EQ.0) THEN
            allocate(ORBPAIRS(2,NPAIRS))
         ELSE
            ORBPAIRS=>OPin(1:2,1:NPAIRS)
         ENDIF
         Call SymSetupExcits_StoreOccPairs(OrbPairs, nPairs, nPr, &
         iMinElec1, iMaxElec1,G1,SymProdInd(1:2,1:3,1:NPR),SymProds, &
          nI,nEl)
!.. Now go through the list of all pairs, finding out how many pairs are to be excluded as they contain some of the
!.. orbitals in this determinant.
!.. We first create a quick lookup table to enable us to quickly check whether a given orbital is in this determinant
!.. in order 1, not order NEL.
         IF(((.not.tStoreStateList).and.(.not.TCOUNT)).or.&
               (tStoreStateList)) THEN
!ILUT is not needed in the setup of the excitation generators for abelian symmetry.
!tStoreStateList will be false for abelian symmetry, unless specified otherwise.
            ILUT(:)=0
            DO I=1,NEL
                ILUT((NI(I)-1)/32)= IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32))
!                WRITE(6,*) (NI(I)-1)/32,
!     &               IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32)),
!     &               ILUT(0:NIfTot)
            ENDDO
!..            DO I=0,NIfTot
!..                WRITE(6,*) "ILUT: ",ILUT(i)
!..                WRITE(6,"(A,Z10)") "LUT: ",ILUT(I)
!..            ENDDO
         ENDIF
!.. Now look through the list of our pairs.  For each pair sym of the complete list which has a 
!.. symmetric product with any of our pair syms, we work out how many allowed pairs there are in
!.. the complete list, and store that value in NALLOWPPS
         IF(STORE(3).EQ.0) THEN
            allocate(NALLOWPPS(3,NSYMPAIRPRODS))
!          WRITE(6,*) "Allocating memory for nallowpps"
         ELSE
            nAllowPPS=>NAPin(1:3,1:nSymPairProds)
         ENDIF
         nExcitTypes=0
         if(.not.tStoreStateList) then
!We can calculate the virtual pairs more easily in abelian symmetry.
            Call SymSetupExcitsAb_CountVProds(nDoub,nExcitTypes,&
                nCl,nPr,SymProds,SymProdInd(1:2,1:3,1:NPR),Classes,&
                ClassCount,  &
                nAllowPPS) ! Calculated not enumerated
         else
             Call SymSetupExcits_CountVirtProds(nDoub, nExcitTypes,nPr,&
        SymProdInd(1:2,1:3,1:NPR),SymProds,nAllowPPS,iLUT)
         end if
         IF(.NOT.ISUHFDET(NI,NEL)) THEN
            if(tAbelianFastExcitGen) then
!This can be done quicker for abelian symmetry, whether or not the state pairs are stored or not.
                CALL SymSetupExcitsAb_CountSing(nSing,nCl,&
                   nExcitTypes,ThisClassCount,Classes)
            else
                Call SymSetupExcits_CountSingles(nSing,nCl,nExcitTypes,&
                   ThisClassCount, ClassCount,Classes)
            endif
         ELSE
          nSing=0
         ENDIF
         NEXCITS=0
         IF(BTEST(ILEVEL,0)) NEXCITS=NEXCITS+NSING
         IF(BTEST(ILEVEL,1)) NEXCITS=NEXCITS+NDOUB
         ICOUNT=NEXCITS
!         WRITE(6,*) "Total number of singles: ",NSING
!         WRITE(6,*) "Total number of doubles: ",NDOUB
!         WRITE(6,*) "Total number of excitations: ",NEXCITS
!         WRITE(6,*) "Total number of excitation types: ",NEXCITTYPES

         IF(TCOUNT) THEN
!.. If we're just counting, we're done, so we get rid of some pointers.
!.. However, we do save the length of the memory required.
!.. EXCITTYPES - Number of excitations which can be created from this 'type' (spin, symmetry, number)
!NExcitTypes will always be =< nsympairprods*3 for doubles and nSymLabels*2 for singles.
            STORE(2)=(NEXCITTYPES*5)
!.. NALLOWPPS
            STORE(3)=3*NSYMPAIRPRODS
!.. ORBPAIRS - Storage of all allowed pairs of orbitals. This will always be less than N*(N+1)/2.
            STORE(4)=2*NPAIRS
!.. SYMPRODIND - Indexing system for ORBPAIRS. NPR is the number of symmetry classes, which for abelian
! symmetry will always be less than or equal to nSymLabels.
            STORE(5)=2*3*(NPR+1)
!.. indicate that these are lengths
            STORE(6)=0
           
            deallocate(nAllowPPS)
            deallocate(OrbPairs)
            deallocate(SymProdInd)
         ELSE
!.. Now allocate memory to store all the excitation types if there hasn't been one already allocated.
!.. This will have to be manually deallocated later.
!.. We store each excitation type as:
!.. 1   TYPE (single=1, double=2)
!.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
!.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
!.. 4   TO   (for single, J in SymLabels(J); for double, J in SYMPAIRPRODS(J) )
!.. 5  COUNT (Total number of excitations in this category)
             IF(STORE(2).EQ.0) THEN
                allocate(EXCITTYPES(5,NEXCITTYPES*5))
             ELSE
                ExcitTypes=>ETin(:5,:nExcitTypes)
             ENDIF

             nExcitTypes=0
             IF(.NOT.ISUHFDET(NI,NEL)) THEN
                IF(BTEST(ILEVEL,0)) THEN
                    IF(tAbelianFastExcitGen) THEN
                        Call SymSetupExcitsAb_StoreSing(&
                                nExcitTypes,nCl,Classes,ThisClassCount,&
                               ExcitTypes)
                    ELSE
                        Call SymSetupExcits_StoreSingles(nExcitTypes,&
                        nCl,Classes,ThisClassCount,ExcitTypes)
                    ENDIF
                ENDIF
             ENDIF
             IF(BTEST(ILEVEL,1)) THEN
                Call SymSetupExcits_StoreDoubles(nPr,nSymPairProds,&
        nAllowPPS,ExcitTypes,nExcitTypes,SymProds,&
           SymProdInd(1:2,1:3,1:NPR))
             ENDIF
!.. Store all the pointers we need
!             STORE(2)=IP_EXCITTYPES
!             STORE(3)=IP_NALLOWPPS
!             STORE(4)=IP_ORBPAIRS
!             STORE(5)=IP_SYMPRODIND
             STORE(6)=NEXCITTYPES

             IF(tAssumeSizeExcitgen) THEN
!If we have an assumed size excitation generator, we don't store nAllowPPS, so we have to deallocate this now.
                 deallocate(nAllowPPS)
             ENDIF

         ENDIF
!..      ENDIF(.NOT.TCOUNT)
!         WRITE(6,*) "Total number of excitation types: ",NEXCITTYPES
         call halt_timer(proc_timer)

      END Subroutine

      subroutine gensymexcitit2par_worker (nI, nel, g1, nbasis, tsetup, &
                    nmem, nJ, ic, store, ilevel, iminelec1, imaxelec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SystemData, only: tAssumeSizeExcitgen
         use SymData, only: SymClassSize,nSymPairProds,nSymLabels
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,ierr
         type(BASISfn) G1(nBasis)
         INTEGER , ALLOCATABLE :: SymProdsTemp(:)
         INTEGER, pointer :: DSTORE(:)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE2 will contain the indices of various components of the excitation 
!generator within the memory NMEM, and is passed to SYMSETUPEXCITS2
         INTEGER STORE2(6)
         INTEGER ICOUNT
         INTEGER, target :: NMEM(*)
         INTEGER NJ(NEL),IC,I
         LOGICAL TSETUP
         INTEGER ILEVEL,Pos1,Pos2,Pos3
         INTEGER iMinElec1, iMaxElec1
         
         IF(TSETUP) THEN
!.. This is the first time we've been called for setup.
!  We pass information back in STORE, but actually hold it in STORE2 during the SYMSETUPEXCITS2
            IF(STORE(1).EQ.0) THEN
               IF(tAssumeSizeExcitgen) THEN
!We are assuming the maximum possible size for the blocks of memory which are required for creation of random
!unweighted excitations. This means that we only need to store three arrays.
!EXCITTYPES - Number of excitations which can be created from this 'type' (spin, symmetry, number)
!NExcitTypes will always be =< nsympairprods*3 for doubles and nSymLabels*2 for singles.
!Total memory needed for the excittypes array is less than 5*(3*nSymPairProds+2*nSymLabels)

!ORBPAIRS - Storage of all allowed pairs of occupied orbitals. This will always be less than (N*(N-1)/2)*2.

!SYMPRODIND - Indexing system for ORBPAIRS. NPR is the number of symmetry classes, which for abelian
! symmetry will always be less than or equal to nSymLabels, so we need a total of (nSymLabels+1)*6

!For 'extras', we want space for ILUT (NBASIS/32+1+NIfY) and CLASSES (SymClassSize*NEL) 
                   STORE(1)=SymClassSize*NEl+(nBasis/32)+1    !ILUT and CLASSES
                   STORE(2)=5*(3*nSymPairProds+2*nSymLabels)    !EXCITTYPES (Could be compressed)
                   STORE(3)=0   !Do not need to store nAllowPPS
                   STORE(4)=NEl*(NEl-1) !Maximum memory for ORBPAIRS
                   STORE(5)=(nSymLabels+1)*6    !SYMPRODIND
                   STORE(6)=0   !Will store nExcitTypes, but zero to indicate length

!We also want some memory to store lengths, but do not need iterator information.
                   NMEM(1)=2+STORE(1)+STORE(2)+STORE(4)+STORE(5)

               ELSE


                   STORE2(1:6)=0
                  allocate(DSTORE(SymClassSize*NEL+(nBasis/32)+1&
                          +SymmetrySize*(NEL*NEL+1)))
                  STORE2(1)=1 ! i.e. the 1st index IP_DSTORE
!.. Just count.
                 ! These store(*) are fakes
                 CALL SYMSETUPEXCITS3(NI,NEL,G1,NBASIS,STORE2,&
                 STORE2(1),STORE2(1),STORE2(1),STORE2(1),&
                 .TRUE.,ICOUNT,DSTORE(1), DSTORE(SymClassSize*NEL+1),&
               DSTORE(SymClassSize*NEL+1+(nBasis/32)+1:),&
                   ILEVEL,iMinElec1,iMaxElec1)

                 deallocate(DSTORE)
                   DO i=2,6
                      STORE(i)=STORE2(i)
                   ENDDO
                   STORE(1)=SymClassSize*NEL+(nBasis/32)+1+&
                     SymmetrySize*(NEL*NEL+1)
                 NMEM(1)=23+STORE(1)+STORE(2)+STORE(3)+STORE(4)+STORE(5)
!               WRITE(6,"(A6,6I5)"), "SIZE:",NMEM(1),(STORE(I),I=1,5)
               ENDIF

            ELSE    !Second setup excitgen run - fill memory

               IF(tAssumeSizeExcitgen) THEN
!If we have an assumed size excitgen, things are a little different. We can only use these 
!assumed sized excitation generators, if we are only
! using them to create random excitations as we are not storing iterator information to save space.

! 1         NEXCIT
! 2         NEXCITTYPES
! 3                                 -    SymClassSize*NEL+NIfTot+3                                 DSTORE
! SymClassSize*NEL+NIfTot+4      -    SymClassSize*NEL+NIfTot+1+15*nSymPairProds+10*nSymLabels+2  EXCITTYPES
! SymClassSize*NEL+NIfTot+1+15*nSymPairProds+10*nSymLabels+3      -    
!SymClassSize*NEL+NIfTot+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+2  ORBPAIRS
! SymClassSize*NEL+NIfTot+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+3   -   
!SymClassSize*NEL+NIfTot+1+15*nSymPairProds+10*nSymLabels+NEL*(NEL-1)+(nSymLabels+1)*6+2  SYMPRODIND

! DSTORE(1)      -   SymClassSize*NEL         CLASSES
! DSTORE(NEL*SymClassSize+1)                     ILUT
                   NMEM(1:2)=0
                   Pos1=SymClassSize*NEL+(nBasis/32)+4    !Beginning of EXCITTYPES
                   Pos2=SymClassSize*NEL+(nBasis/32)+4+15*nSymPairProds+&
                        10*nSymLabels                  !Beginning of ORBPAIRS
                   Pos3=SymClassSize*NEL+(nBasis/32)+4+15*nSymPairProds+&
                        10*nSymLabels+NEL*(NEL-1)      !Beginning of SYMPRODIND

!STORE is now used to hold pointers to the arrays wanted. Fill this.
                   STORE2(1)=3   !This is the beginning of the DSTORE array
                   STORE2(2)=Pos1
                   STORE2(3)=0  !nAllowPPS not stored
                   STORE2(4)=Pos2
                   STORE2(5)=Pos3

                  DSTORE=>NMEM(STORE2(1):STORE2(1)+&
                     SymClassSize*NEL+(nBasis/32)+1&
                          +SymmetrySize*(NEL*NEL+1))  !Point to DSTORE

!Since we do not store symprods, we have to allocate more memory to store it temporarily.
                   ALLOCATE(SymProdsTemp(SymmetrySize*(NEL*NEL+1)), stat=ierr)
                   IF(ierr.ne.0) THEN
                       CALL Stop_All("GenSymExcitIt2Par","Cannot &
                            &allocate SymProdsTemp Memory")
                   ENDIF

!Call SymSetupExcits3. Since we are not storing symprods, DStore is shorter than normal, 
!and we just pass through a temporary array to hold it.
                    ! n.b. nAllowPPS not extant
                    CALL SYMSETUPEXCITS3(NI,NEL,G1,NBASIS,STORE2,&
                       NMEM(STORE2(5)),NMEM(STORE2(2)),NMEM(1),&
                       NMEM(STORE2(4)),&
                       .FALSE.,ICOUNT,DSTORE(1), &
                       DSTORE(SymClassSize*NEL+1),&
                       SymProdsTemp,ILEVEL,iMinElec1,iMaxElec1)

                    DEALLOCATE(SymProdsTemp)

                    NMEM(1)=iCount      !Keep record of number of excitations...
                    NMEM(2)=STORE2(6)   !...and number of types of excitation.

                ELSE
!.. The second setup.  Now NMEM is allocated, we store all the info
!.. NMEM is as follows:
!..   1           -  5           STORE(1:5)
!..   6           -  6           STORE(6) = NEXCITTYPES
!.. Data for the iterators
!..   7  IEXCIT
!..   8  ISPN
!..   9  IFROM
!..   10 ITO
!..   11 I
!..   12 J
!..   13 K
!..   14 L
!..   15 ICC(1:4)
!..   19 LS(1:2,1:2)  ([1 or 2], [A or B])
!..   23 NEXCIT      total number of excitations
!..  (STORE(1)=24)-  STORE(2)-1  DSTORE
!..   STORE(2)    -  STORE(3)-1  EXCITTYPES
!..   STORE(3)    -  STORE(4)-1  NALLOWPPS
!..   STORE(4)    -  STORE(5)-1  ORBPAIRS
!..   STORE(5)    -  ...         SYMPRODIND

!..   DSTORE(1)   -  DSTORE(NEL*SymClassSize)      CLASSES
!..   DSTORE(NEL*SymClassSize+1)                   ILUT
!..                            - DSTORE(NEL*SymClassSize+1 +nBasis/32+1)

!..   DSTORE(NEL*SymClassSize+1 +nBasis/32+2)      SymProds
!..                            - DSTORE(NEL*SymClassSize+1 +nBasis/32+2+ SymmetrySize*(1+nPr)
!..                              This last is the end of DSTORE i.e. MEM(STORE(2)-1)

                   NMEM(1:23)=0
                   ICOUNT=24  
!.. Put the indices in store
                   DO I=1,5
                      STORE2(I)=ICOUNT
                      ICOUNT=ICOUNT+STORE(I)
                   ENDDO
                   NMEM(6:ICOUNT-1)=0
                   STORE2(6)=0
                   NMEM(11)=-1
                   NMEM(7)=0
!                  DSTORE=>NMEM(STORE2(1):STORE2(1)+&
!                     SymClassSize*NEL+(nBasis/32)+1&
!                          +SymmetrySize*(NEL*NEL+1))  !Point to DSTORE
                  DSTORE=>NMEM(STORE2(1):STORE2(2)-1) !point to DSTORE
!                   CALL DUMPIMEMORY(6,NMEM,ICOUNT-1)
!!      SUBROUTINE SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
!!     &   TCOUNT,ICOUNT,CLASSES,ILUT,SYMPRODS,ILEVEL)
                 CALL SYMSETUPEXCITS3(NI,NEL,G1,NBASIS,STORE2,&
                       NMEM(STORE2(5)),NMEM(STORE2(2)),NMEM(STORE2(3)),&
                       NMEM(STORE2(4)),&
                  .FALSE.,ICOUNT,DSTORE(1), DSTORE(SymClassSize*NEL+1),&
             DSTORE(SymClassSize*NEL+1+(nBasis/32)+1:),&
     &       ILEVEL,iMinElec1,&
                   iMaxElec1)
                   NMEM(6)=STORE2(6)
                   NMEM(23)=ICOUNT
                   IC=ICOUNT
!.. and now instead the indices within NMEM
                   ICOUNT=24
                   DO I=1,5
                      NMEM(I)=ICOUNT
                      ICOUNT=ICOUNT+STORE(I)
                   ENDDO
!.. Second setup finally complete.
!               CALL DUMPIMEMORY(6,NMEM,ICOUNT-1)
!               WRITE(6,*) "DST",LOC(NMEM(NMEM(1)))
                ENDIF

            ENDIF
         ELSE

             IF(tAssumeSizeExcitGen) THEN
                 CALL Stop_All("GenSymExcitIt2","Cannot use assumed size &
                  &ExcitGens when enumerating all determinants.")
             ENDIF
!.. Actually generate a det
!            WRITE(6,"(A,Z10,8I4)",advance='no') "GET",LOC(NMEM(1)),
!     &         (NMEM(I),I=7,14)
!                  DSTORE=>NMEM(NMEM(1):NMEM(1)+&
!                     SymClassSize*NEL+(nBasis/32)+1&
!                          +SymmetrySize*(NEL*NEL+1))  !Point to DSTORE
                  DSTORE=>NMEM(NMEM(1):NMEM(2)-1) !Point to DStore
!!      SUBROUTINE SYMGENEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
!!     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,
!!     &               I,J,K,L,ICC,LS,
!!     &               NK,IC)
            CALL SYMGENEXCITIT2(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),&
              NMEM(NMEM(5)),DSTORE(SymClassSize*NEL+1),NMEM(NMEM(4)),&
              NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),&
              NMEM(13),NMEM(14),NMEM(15),NMEM(19),NJ,IC,iMinElec1,&
              iMaxElec1)
!            CALL WRITEDET(6,NI,NEL,.FALSE.)
!            WRITE(6,"(A)",advance='no'), "->"
!            IF(NJ(1).NE.0) THEN
!               CALL WRITEDET(6,NJ,NEL,.TRUE.)
!            ELSE 
!               WRITE(6,*) "(    0,)"
!            ENDIF
         ENDIF
      END subroutine

END MODULE SymExcit2
