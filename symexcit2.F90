MODULE SymExcit2
      IMPLICIT NONE

      TYPE ExcitWeight
! The orbitals excited from
         INTEGER I,J
! The orbitals excited to
         INTEGER A,B
         REAL*8 WEIGHT
      END TYPE ExcitWeight
! Size in terms of reals.
      PARAMETER ExcitWeightSize=3
      CONTAINS 

!  Enumerate the excitations and weights of excitations of a given ExcitType.
      SUBROUTINE EnumExcitWeights(ExcitType,iFromIndex,iLUT,ews,OrbPairs,SymProdInd,Norm,iCount,G1,NBASISMAX,UMAT,NBASIS)
         USE HElement
         INCLUDE 'sym.inc'
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
         INTEGER NBASISMAX(*)
         TYPE(BasisFN) G1(*)
         TYPE(HElement) UMAT(*)
         REAL*8 Norm
         TYPE(ExcitWeight) ews(*)
         INTEGER K
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
            L1B=BTEST(ILUT(ICC1/32),MOD(ICC1,32)-1)
            L1A=BTEST(ILUT(ICC2/32),MOD(ICC2,32)-1)
            L2B=BTEST(ILUT(ICC3/32),MOD(ICC3,32)-1)
            L2A=BTEST(ILUT(ICC4/32),MOD(ICC4,32)-1)
!.. L1B is set if the beta of the first virtual is in NI, i.e. is disallowed
            IF(ISPN.EQ.1) THEN
!.. If both virtuals aren't the samem and neither are in NI, then allow
               IF(ICC1.NE.ICC3.AND..NOT.(L1B.OR.L2B)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC1,                                           &
     &                  ICC3,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,NBASIS)
               ENDIF
            ELSEIF(ISPN.EQ.2) THEN
!.. If neither virtuals are in NI, then allow
               IF(.NOT.(L1B.OR.L2A)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC1,                                           &
     &                  ICC4,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,NBASIS)
               ENDIF
!.. If neither virtuals are in NI, and they're not the same(which would give
!.. us the same excitation as previously), then allow
               IF(.NOT.(L1A.OR.L2B).AND.ICC1.NE.ICC3) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC2,                                           &
     &                  ICC3,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,NBASIS)
               ENDIF
            ELSEIF(ISPN.EQ.3) THEN
!.. If both virtuals aren't the samem and neither are in NI, then allow
               IF(ICC1.NE.ICC3.AND..NOT.(L1A.OR.L2A)) THEN
                  CALL AddExcitWeight(                                  &
     &                  ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+iFromIndex),&
     &                  ICC2,                                           &
     &                  ICC4,                                           &
     &                  ews,Norm,ICOUNT,G1,NBASISMAX,UMAT,NBASIS)
               ENDIF
            ENDIF
         ENDDO
      END
! Add the weight of the excitation to the list in ExWeights
! I,J are from, K,L are to
      SUBROUTINE AddExcitWeight(I,J,A,B,ExWeights,Norm,iCount,G1,NBASISMAX,UMAT,NBASIS)
         USE HElement
         INTEGER I,J,A,B
         REAL*8 R,Norm
         INTEGER NBASISMAX(*),NBASIS
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         TYPE(HElement) UMAT(*)
         TYPE(ExcitWeight) ExWeights(iCount+1)
         INTEGER iCount
         CALL ExcitWeighting(I,J,A,B,R,G1,NBASISMAX,UMAT,NBASIS)
         iCount=iCount+1
         ExWeights(iCount)%I=I
         ExWeights(iCount)%J=J
         ExWeights(iCount)%A=A
         ExWeights(iCount)%B=B
         ExWeights(iCount)%Weight=R
         Norm=Norm+R
      END
         
!        A sub called to generate an unnormalised weight for a given ij->kl excitation
!          We return a function of the U matrix element (|<ij|u|kl>|^2)^G_VMC_EXCITWEIGHT
      SUBROUTINE EXCITWEIGHTING(I,J,K,L,WEIGHT,G1,NBASISMAX,UMAT,NBASIS)
         USE HElement
         IMPLICIT NONE
         INTEGER NBASISMAX(5,3),NBASIS
         INCLUDE 'basis.inc'
! We get G_VMC_EXCITWEIGHT from vmc.inc
         INCLUDE 'vmc.inc'
         TYPE(BasisFN) G1(NBASIS)
!  We fake ISS
         INTEGER ISS
         INTEGER IDI,IDJ,IDK,IDL
         INTEGER I,J,K,L
         REAL*8 WEIGHT
         TYPE(HElement) GetUMatEl
         TYPE(HElement) UMAT(*),W
         IF(G_VMC_EXCITWEIGHT.EQ.0.D0) THEN
            WEIGHT=1.D0
            RETURN
         ENDIF
         ISS=NBASISMAX(2,3)
         CALL GTID(NBASISMAX,I,IDI)
         CALL GTID(NBASISMAX,J,IDJ)
         CALL GTID(NBASISMAX,K,IDK)
         CALL GTID(NBASISMAX,L,IDL)
         W=GetUMatEl(NBASISMAX,UMAT,0.0,NBASIS,ISS,G1,IDI,IDJ,IDK,IDL)
         WEIGHT=SQ(W)**G_VMC_EXCITWEIGHT
         RETURN
      END

!We wish to calculate what excitation class the excitation NI->NJ falls into, with the appropriate 
! IFROM class and index within that, IFROMINDEX, and ITO class, and index within that. ITOINDEX
! After that, we generate the probability that nJ would be an excitation from nI.
      SUBROUTINE GenExcitProbInternal(nI,nJ,nEl,G1,nBasisMax,UMat,nBasis,OrbPairs,SymProdInd,iLUT,SymProds,ExcitTypes,iTotal,pGen)
         USE HElement
         IMPLICIT NONE
         INTEGER iExcit(2,2)
         LOGICAL L
         INTEGER nI(nEl),nJ(nEl),nEl,nBasis
         INTEGER iFrom,iFromIndex,iTo,iToIndex
         INCLUDE 'sym.inc'
         TYPE(BasisFn) G1(nBasis)
         TYPE(Symmetry) Sym,Prod
         TYPE(Symmetry) SYMPROD
         TYPE(Symmetry) SymProds(0:*)
         LOGICAL SYMEQ
         INTEGER nBasisMax(*)
         TYPE(HElement) UMat(*)
         TYPE(ExcitWeight) ews(*)
         POINTER(IP_ews,ews)
         INTEGER iLUT(*)
         INTEGER OrbPairs(2,*)
         INTEGER SymProdInd(2,3,0:*)
         INTEGER ExcitTypes(5,*)
         REAL*8 Norm
         INTEGER K,NPR,I
         INTEGER iSpn,iCount,nToPairs,iTotal
         REAL*8 pGen
         iExcit(1,1)=2
         CALL GETEXCITATION(nI,nJ,nEl,iExcit,L)
!EXCIT(1,*) are the ij... in NI, and EXCIT(2,*) the ab... in NJ
         IF(iExcit(1,1).LE.0) THEN
! More than a double excitation, so not allowed
            IFROM=0
            RETURN
         ENDIF
!  See if we have a single
         IF(iExcit(1,2).EQ.0) THEN
            pGen=1.D0/iTotal
            RETURN
         ENDIF
!  We have a double
! Now check the sym prod of i and j to work out the FROM
         Sym=SYMPROD( G1(iExcit(1,1))%Sym,G1(iExcit(1,2))%Sym)
         iFrom=1
         DO WHILE(.NOT.SYMEQ(SymProds(iFrom),Sym))
            iFrom=iFrom+1
         ENDDO
!.. We've worked out what class the IFROM was.  Now work out which member of the class it is
!.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym SYMPRODS(I) in ORBPAIRS
!.. SYMPRODIND(2,ISPN,I) contains the number of such elements
!.. There are three values of ISPN.  ISPN=1 is beta beta, ISPN=2 is alpha, beta and beta, alpha.  ISPN=3 is alpha alpha
         iSpn=(G1(iExcit(1,1))%Ms+G1(iExcit(1,2))%Ms)/2+2
         DO K=1,SymProdInd(2,iSpn,iFrom)
            IF(iExcit(1,1).EQ.OrbPairs(1,SymProdInd(1,iSpn,iFrom)+K).AND. &
     &         iExcit(1,2).EQ.OrbPairs(2,SymProdInd(1,iSpn,iFrom)+K)) THEN
!.. We've found the index.  Save and leave
               iFromIndex=K
               EXIT
            ENDIF
         ENDDO
!.. Now to find the appropriate iTo index:
         Sym=SYMPROD( G1(iExcit(2,1))%Sym, G1(iExcit(2,2))%Sym)
!..  Find which SymPairProd it corresponds to 
!.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
!.. symlabels (listed in SymStatePairs), and the index of the start of this list
!.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
!.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of symlabels whose prod is of that symmetry.
         CALL FindSymProd(Sym,SymPairProds,nSymPairProds,iTo)
!.. Now work out the index of the (a,b) pair within the prod
         K=1
         DO WHILE(ExcitTypes(1,K).NE.2.OR.ExcitTypes(2,K).NE.iSpn.OR.ExcitTypes(3,K).NE.iFrom.OR.ExcitTypes(4,K).NE.iTo)
            K=K+1
         ENDDO
!.. K is now the excitation
         nToPairs=ExcitTypes(5,K)/SymProdInd(2,iSpn,iFrom)
         CALL MEMORY(IP_ews,ExcitWeightSize*nToPairs,'ExcitWGEPI')
         iCount=0
         Norm=0.D0
         CALL EnumExcitWeights(ExcitTypes(1,K),iFromIndex,iLUT,ews,OrbPairs,SymProdInd,Norm,iCount,G1,nBasisMax,UMat,nBasis)
!.. Find the (a,b) pair
!.. The prob of all possible excitations in this iTo 
         pGen=(iCount+0.D0)/iTotal
         DO I=1,nToPairs
            IF(ews(I)%A.EQ.iExcit(2,1).AND.ews(I)%B.EQ.iExcit(2,2)) THEN
               pGen=pGen*ews(I)%Weight
               EXIT
            ENDIF
         ENDDO
         pGen=pGen/Norm
         CALL FREEM(IP_ews)
      END
END MODULE SymExcit2
