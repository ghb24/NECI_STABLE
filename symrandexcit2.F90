MODULE GenRandSymExcitNUMod
!This is a new random excitation generator. It still creates excitations with a normalised and calculable probability,
!but these probabilities are not uniform. They should also be quicker to generate, and have small or no excitation generators
!to store. nI is the root determinant and nJ is the excitation.
!This currently is only available for abelian symmetry. This is because I label the symmetries as integers, but this
!algorithm could theoretically be applied relatively easily to non-abelian symmetry. The symmetries would have to be stored
!in a better way. For molecular systems, symmetry irrep is just symmetry_label-1. Symlabels will need to be used for more
!complicated symmetry tables.

!THEORY

      !  P[i,j,a,b] = P_doub x ( P[i,j] x P[a,b|i,j] ) x (M-N)/(M-N-Q)
      !  P[a,b|i,j] = P[b|i,j] x P[a|b,i,j] + P[a|i,j] x P[b|a,i,j]
      !  P[i,j] = 2/N(N-1) always - others are more difficult to calculate
      !  Difficulties arise if there is ever no allowed final unoccupied orbital to pick.

      !  Spin and spatial symmetry needs to be considered.
      !  The fourth orbital chosen (i.e. one of the unoccupied ones) will have its
      !  symmetry class specified by the classes of the other ones.

      !  i and j are picked randomly. They can be alpha,beta | beta,alpha | alpha,alpha | beta,beta
      !  a is then picked from all orbitals with the allowed spin. This means that a beta,beta i,j pair
      !  will force a to also be beta. There are also a number of forbidden a orbitals which are counted.
      !  These are forbidden since they have no possible b orbital which will give rise to a symmetry and
      !  spin allowed unoccupied a,b pair. The number of these orbitals, Q, is needed to calculate the
      !  normalised probability of generating the excitation.
    use SystemData, only: ALAT,iSpinSkip,tFixLz,iMaxLz,NIfTot
    use SystemData, only: nEl,G1, nBasis,nBasisMax,tNoSymGenRandExcits,tMerTwist
    use SystemData, only: Arr,nMax,tCycleOrbs,nOccAlpha,nOccBeta,ElecPairs,MaxABPairs
    use IntegralsData, only: UMat
    use Determinants, only: GetHElement4
    use SymData, only: nSymLabels,TwoCycleSymGens
    use SymData, only: SymLabelList,SymLabelCounts
    use mt95 , only : genrand_real2
    use SymExcitDataMod 
    use HElem
    IMPLICIT NONE

    INTEGER , PARAMETER :: r2=kind(0.d0)

    contains

!This routine is the same as GenRandSymExcitNU, but you can pass in the ClassCount arrays, so they do not have to be recalculated each time
!we want an excitation. If tFilled is false, then it wil assume that they are unfilled and calculate them. It will then return the arrays
!with tFilled = .true. for use in the next excitation.
!The two arrays want to be integers, both of size (2,1:nSymLabels)
    SUBROUTINE GenRandSymExcitScratchNU(nI,iLut,nJ,pDoub,IC,ExcitMat,tParity,exFlag,pGen,ClassCount2,ClassCountUnocc2,tFilled)
        INTEGER :: nI(NEl),nJ(NEl),IC,ExcitMat(2,2),exFlag
        INTEGER :: ClassCount2(ScratchSize),ElecsWNoExcits
        INTEGER :: ClassCountUnocc2(ScratchSize)
!        INTEGER , SAVE :: Iter=0
        INTEGER :: ILUT(0:NIfTot),i!,DetSym
        LOGICAL :: tNoSuccess,tParity,tFilled
        REAL*8 :: pDoub,pGen,r
        CHARACTER , PARAMETER :: this_routine='GenRandSymExcitNU'

!        Iter=Iter+1
!        WRITE(6,*) Iter,tFilled,nSymLabels

        IF(.not.tFilled) THEN
            IF(.not.TwoCycleSymGens) THEN
!Currently only available for molecular systems, or without using symmetry.
                IF(.not.tNoSymGenRandExcits) THEN
                    WRITE(6,*) "GenRandSymExcitNU can only be used for molecular systems"
                    WRITE(6,*) "This is because of difficulties with other symmetries setup."
                    WRITE(6,*) "If you want to use these excitation generators, then add NOSYMGEN to the input to ignore symmetry while generating excitations."
                    CALL FLUSH(6)
                    CALL Stop_All(this_routine,"GenRandSymExcitNU can only be used for molecular systems using symmetry")
                ENDIF
            ENDIF

!First, we need to do an O[N] operation to find the number of occupied alpha electrons, number of occupied beta electrons
!and number of occupied electrons of each symmetry class and spin. This is similar to the ClassCount array.
!This has the format (Spn,sym), where Spin=1,2 corresponding to alpha and beta.
!For molecular systems, sym runs from 0 to 7. This is NOT general and should be made so using SymLabels.
!This could be stored to save doing this multiple times, but shouldn't be too costly an operation.
            CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)

            tFilled=.true.
        ENDIF
!        WRITE(6,*) ClassCount2(:)
!        WRITE(6,*) "***" 
!        WRITE(6,*) ClassCountUnocc2(:)
!        WRITE(6,*) "***" 

!ExFlag is 1 for singles, 2 for just doubles, and 3 for both.
        IF(ExFlag.eq.3) THEN
!Choose whether to generate a double or single excitation. Prob of generating a double is given by pDoub.
            pDoubNew=pDoub
            IF(tNoSingsPossible) THEN
!This will check if there are any possible single excitations from this determinant
!If there are not, then this will change pDoubNew so that it = 1 and only doubles will be generated.
                CALL CheckIfSingleExcits(ElecsWNoExcits,ClassCount2,ClassCountUnocc2)
            ENDIF

            IF(pDoubNew.gt.1.D0) THEN
                CALL Stop_All(this_routine,"pDoub is greater than 1")
!            ELSEIF(pDoubNew.eq.1.D0) THEN
!                WRITE(6,*) "*****",nI,iLut
!                STOP
            ENDIF
                

            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            IF(r.lt.pDoubNew) THEN
!A double excitation has been chosen to be created.
                IC=2
            ELSE
                IC=1
            ENDIF
        ELSEIF(ExFlag.eq.2) THEN
            IC=2
            pDoubNew=1.D0
        ELSEIF(ExFlag.eq.1) THEN
            IC=1
            tNoSingsPossible=.false.
            pDoubNew=0.D0
        ELSE
            CALL Stop_All(this_routine,"Error in choosing excitations to create.")
        ENDIF

        IF(IC.eq.2) THEN
            CALL CreateDoubExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
        ELSE
            CALL CreateSingleExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen,ElecsWNoExcits)
!            IF(pGen.eq.-1.D0) THEN
!NOTE: ghb24 5/6/09 Cannot choose to create double instead, since you could have chosen a double first and it would have a different pGen.
!                IF((ExFlag.ne.3).or.tHPHF) THEN
!If using HPHF wavefunctions, then we do not want to do this, since it will affect the generation probabilities calculated using CalcNonUniPGens.
!                    CALL Stop_All("GenRandSymExcitNU","Found determinant with no singles, but can only have got here from single. Should never be in this position! (or HPHF is on and this will screw with pGens)")
!                ENDIF
!                pDoubNew=1.D0
!                IC=2
!                CALL CreateDoubExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
!            ENDIF

        ENDIF
!        DetSym=0
!        do i=1,NEl
!            DetSym=IEOR(DetSym,INT(G1(nJ(i))%Sym%S,4))
!        enddo
!        IF(DetSym.ne.6) THEN
!            CALL Stop_All("GenRand","WrongSym")
!        ENDIF

    END SUBROUTINE GenRandSymExcitScratchNU

    SUBROUTINE GenRandSymExcitNU(nI,iLut,nJ,pDoub,IC,ExcitMat,TParity,exFlag,pGen)
        INTEGER :: nI(NEl),nJ(NEl),IC,ExcitMat(2,2),exFlag
        INTEGER :: ClassCount2(ScratchSize),ElecsWNoExcits
        INTEGER :: ClassCountUnocc2(ScratchSize)
        INTEGER :: ILUT(0:NIfTot),i
        LOGICAL :: tNoSuccess,tParity
        REAL*8 :: pDoub,pGen,r
        CHARACTER , PARAMETER :: this_routine='GenRandSymExcitNU'

!        WRITE(6,*) "nSymlabels:", nSymLabels
!        WRITE(6,*) "SymLabelList: "
!        do i=1,nbasis/2
!            WRITE(6,*) SymLabelList(i)
!        enddo
!        WRITE(6,*) "SymLabelCounts:"
!        do i=1,nSymLabels
!            WRITE(6,*) SymLabelCounts(1,i),SymLabelCounts(2,i)
!        enddo
!        WRITE(6,*) "G1:"
!        do i=1,nBasis
!            WRITE(6,*) G1(i)%Sym%S
!        enddo
!        CALL FLUSH(6)
!        STOP

        IF(.not.TwoCycleSymGens) THEN
!Currently only available for molecular systems, or without using symmetry.
            IF(.not.tNoSymGenRandExcits) THEN
                WRITE(6,*) "GenRandSymExcitNU can only be used for molecular systems"
                WRITE(6,*) "This is because of difficulties with other symmetries setup."
                WRITE(6,*) "If you want to use these excitation generators, then add NOSYMGEN to the input to ignore symmetry while generating excitations."
                CALL FLUSH(6)
                CALL Stop_All(this_routine,"GenRandSymExcitNU can only be used for molecular systems using symmetry")
            ENDIF
        ENDIF

!First, we need to do an O[N] operation to find the number of occupied alpha electrons, number of occupied beta electrons
!and number of occupied electrons of each symmetry class and spin. This is similar to the ClassCount array.
!This has the format (Spn,sym), where Spin=1,2 corresponding to alpha and beta.
!For molecular systems, sym runs from 0 to 7. This is NOT general and should be made so using SymLabels.
!This could be stored to save doing this multiple times, but shouldn't be too costly an operation.
        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)

!ExFlag is 1 for singles, 2 for just doubles, and 3 for both.
        IF(ExFlag.eq.3) THEN
!Choose whether to generate a double or single excitation. Prob of generating a double is given by pDoub.
            pDoubNew=pDoub
            IF(tNoSingsPossible) THEN
!This will check if there are any possible single excitations from this determinant
!If there are not, then this will change pDoubNew so that it = 1 and only doubles will be generated.
                CALL CheckIfSingleExcits(ElecsWNoExcits,ClassCount2,ClassCountUnocc2)
            ENDIF
            IF(pDoubNew.gt.1.D0) CALL Stop_All(this_routine,"pDoub is greater than 1")

            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            IF(r.lt.pDoubNew) THEN
!A double excitation has been chosen to be created.
                IC=2
            ELSE
                IC=1
            ENDIF
        ELSEIF(ExFlag.eq.2) THEN
            IC=2
            pDoubNew=1.D0
        ELSEIF(ExFlag.eq.1) THEN
            IC=1
            tNoSingsPossible=.false.
            pDoubNew=0.D0
        ELSE
            CALL Stop_All(this_routine,"Error in choosing excitations to create.")
        ENDIF

        IF(IC.eq.2) THEN
            CALL CreateDoubExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
        ELSE
            CALL CreateSingleExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen,ElecsWNoExcits)
!            IF(pGen.eq.-1.D0) THEN
!NOTE: ghb24 5/6/09 Cannot choose to create double instead, since you could have chosen a double first and it would have a different pGen.
!                IF(ExFlag.ne.3) THEN
!                    CALL Stop_All("GenRandSymExcitNU","Found determinant with no singles, but can only have got here from single. Should never be in this position!")
!                ENDIF
!                pDoubNew=1.D0
!                IC=2
!                CALL CreateDoubExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
!            ENDIF

        ENDIF


    END SUBROUTINE GenRandSymExcitNU

    SUBROUTINE CreateDoubExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
        INTEGER :: nI(NEl),nJ(NEl),ExcitMat(2,2),NExcitOtherWay,OrbB
        INTEGER :: ClassCount2(ScratchSize)
        INTEGER :: ClassCountUnocc2(ScratchSize),SumMl,MlA,MlB
        INTEGER :: ILUT(0:NIfTot),NExcitB,SpinOrbA,OrbA,SymB,NExcitA
        INTEGER :: Elec1Ind,Elec2Ind,SymProduct,iSpn,ForbiddenOrbs,SymA
        REAL*8 :: pGen
        LOGICAL :: tParity,tAOrbFail

!First, we need to pick an unbiased distinct electron pair.
!These have symmetry product SymProduct, and spin pair iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha
        CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,SumMl,-1)

!This routine finds the number of orbitals which are allowed by spin, but not part of any spatial symmetry allowed unoccupied pairs.
!This number is needed for the correct normalisation of the probability of drawing any given A orbital since these can be chucked and redrawn.
        IF(tNoSymGenRandExcits) THEN
            CALL FindNumForbiddenOrbsNoSym(ForbiddenOrbs,ClassCountUnocc2,iSpn)
        ELSE
            CALL FindNumForbiddenOrbs(ForbiddenOrbs,ClassCountUnocc2,SymProduct,iSpn,SumMl)
        ENDIF

!Now we have to pick the first unoccupied orbital. If an orbital is not present in any allowed pairs, it is chucked and a new one drawn.
!The number NExcit is the number of unoccupied orbitals that the orbital was chosen from (including symmetry-forbidden orbital pairs)
!Arguments:     NExcit = Number of possible spin-allowed unoccupied spinorbitals, including forbidden orbs (these will be chucked)
!               SpinOrbA = Spin of the chosen spin-orbital. 1 is alpha, -1 is beta.
!               OrbA = Index of the chosen spinorbital
!               SymB = Symmetry required of the second unoccupied spinorbital, so that sym(A) x sym(B) = SymProduct
!               SymProduct = (intent in), the symmetry of electron pair chosen
!               SumMl = Sum of the Ml values for the two picked electrons
!               MlA = Ml of the a orbital chosen
!               MlB = Required Ml of the b orbital still to be chosen
        CALL PickAOrb(nI,iSpn,ILUT,ClassCountUnocc2,NExcitA,Elec1Ind,Elec2Ind,SpinOrbA,OrbA,SymA,SymB,SymProduct,SumMl,MlA,MlB,ForbiddenOrbs,tAOrbFail)
        IF(tAOrbFail) THEN
!            WRITE(6,*) "A ORB FAIL"
            nJ(1)=0
            pGen=HUGE(0.D0)
            RETURN
        ENDIF

!This routine will pick an unoccupied orbital at random from a specified spin and symmetry class.
!There should definitely be a possible spinorbital, since A was chosen so that there was one.
!We have to make sure with alpha/alpha or beta/beta pairs and when SymProduct=0, that we don't choose the same unoccupied orbital.
!If we do this, then we should chuck and redraw, since there should definitely be another allowed spinorbital in the class.
!We return the number of allowed B's for the A we picked in NExcitB, however we also need to know the number of allowed A's if we
!had picked B first. This will be returned in NExcitOtherWay.
        CALL PickBOrb(nI,iSpn,ILUT,ClassCountUnocc2,SpinOrbA,OrbA,SymA,OrbB,SymB,NExcitB,SymProduct,SumMl,MlA,MlB,NExcitOtherWay)

        CALL FindNewDet(nI,nJ,Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat,tParity)

!        IF(Iter.eq.13) THEN
!            WRITE(6,*) "***",ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay,SymProduct 
!        ENDIF
        CALL FindDoubleProb(ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay,pGen)

    END SUBROUTINE CreateDoubExcit

!This routine creates the final determinant.
    SUBROUTINE FindNewDet(nI,nJ,Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat,tParity)
        INTEGER :: nI(NEl),nJ(NEl),Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat(2,2)
        INTEGER :: iGetExcitLevel_2,ExcitLevel
        LOGICAL :: tParity,IsValidDet,SymAllowed

!First construct ExcitMat
        ExcitMat(1,1)=Elec1Ind
        ExcitMat(2,1)=OrbA
        ExcitMat(1,2)=Elec2Ind
        ExcitMat(2,2)=OrbB
        nJ(:)=nI(:)
        CALL FindExcitDet(ExcitMat,nJ,2,tParity)
!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
!        CALL IsSymAllowedExcit(nI,nJ,2,ExcitMat,SymAllowed)

    END SUBROUTINE FindNewDet

!This routine finds the probability of creating the excitation. See the header of the file for more information on how this works.
    SUBROUTINE FindDoubleProb(ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay,pGen)
        INTEGER :: ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay
        REAL*8 :: pGen!,PabGivenij

!        PabGivenij=(1.D0/real((NExcitA-ForbiddenOrbs),r2))*((1.D0/real(NExcitB,r2))+(1.D0/real(NExcitOtherWay,r2)))
!        pGen=pDoubNew*(1.D0/real(ElecPairs,r2))*PabGivenij
        pGen=pDoubNew*((1.D0/real(NExcitB,r2))+(1.D0/real(NExcitOtherWay,r2)))/(REAL((ElecPairs*(NExcitA-ForbiddenOrbs)),r2))

    END SUBROUTINE FindDoubleProb

    SUBROUTINE PickBOrb(nI,iSpn,ILUT,ClassCountUnocc2,SpinOrbA,OrbA,SymA,OrbB,SymB,NExcit,SymProduct,SumMl,MlA,MlB,NExcitOtherWay)
        INTEGER :: nI(NEl),iSpn,SpinOrbA,OrbA,SymB,NExcit,SymProduct,NExcitOtherWay
        INTEGER :: OrbB,Attempts,SpinOrbB,ChosenUnocc,Ind,SumMl,MlA,MlB
        INTEGER :: ILUT(0:NIfTot),SymA,nOrbs,z,i
        INTEGER :: ClassCountUnocc2(ScratchSize)
        REAL*8 :: r

!We want to calculate the number of possible B's given the symmetry and spin it has to be since we have already picked A.
!We have calculated in NExcit the number of orbitals available for B given A, but we also need to know the number of orbitals to choose from for A IF
!we had picked B first.
        IF(iSpn.eq.2) THEN
!If iSpn=2, then we want to find a spinorbital of the opposite spin of SpinOrbA
            IF(SpinOrbA.eq.-1) THEN
!We have already picked a beta orbital, so now we want to pick an alpha orbital. Find out how many of these there are.
                NExcit=ClassCountUnocc2(ClassCountInd(1,SymB,MlB))
!                WRITE(6,*) "NExcit",NExcit
                NExcitOtherWay=ClassCountUnocc2(ClassCountInd(2,SymA,MlA))
!                WRITE(6,*) "NExcitOtherWay:", NExcitOtherWay
!                SpinOrbB=0  !This is defined differently to SpinOrbA. 0=Alpha, -1=Beta.
                SpinOrbB=1  !This is defined differently to SpinOrbA. 1=Alpha, 2=Beta.
            ELSE
!Want to pick an beta orbital.
                NExcit=ClassCountUnocc2(ClassCountInd(2,SymB,MlB))
                NExcitOtherWay=ClassCountUnocc2(ClassCountInd(1,SymA,MlA))
                SpinOrbB=2
            ENDIF
        ELSEIF(iSpn.eq.1) THEN
!Definitely want a beta orbital
            NExcit=ClassCountUnocc2(ClassCountInd(2,SymB,MlB))
            NExcitOtherWay=ClassCountUnocc2(ClassCountInd(2,SymA,MlA))
            SpinOrbB=2
        ELSE
!Definitely want an alpha orbital
            NExcit=ClassCountUnocc2(ClassCountInd(1,SymB,MlB))
            NExcitOtherWay=ClassCountUnocc2(ClassCountInd(1,SymA,MlA))
            SpinOrbB=1
        ENDIF

        IF((iSpn.ne.2).and.(SymProduct.eq.0).and.(MlA.eq.MlB)) THEN
!In this case, we need to check that we do not pick the same orbital as OrbA. If we do this, then we need to redraw.
!Only when SymProduct=0 will the classes of a and b be the same, and the spins will be different if iSpn=2, so this is the only possibility of a clash.
            NExcit=NExcit-1     !Subtract 1 from the number of possible orbitals since we cannot choose orbital A.
            NExcitOtherWay=NExcitOtherWay-1     !The same goes for the probabilities the other way round.
        ENDIF

!All orbitals with the specified symmetry and spin should be allowed unless it is OrbA. There will be NExcit of these. Pick one at random.
!Check that orbital is not in ILUT and is not = OrbA (Although this can only happen in the circumstance indicated earlier).
!Now we need to choose the final unoccupied orbital.
!There are two ways to do this. We can either choose the orbital we want out of the NExcit possible unoccupied orbitals.
!It would then be necessary to cycle through all orbitals of that symmetry and spin, only counting the unoccupied ones to find the correct determinant.
!Alternatively, we could pick orbitals at random and redraw until we find an allowed one. This would probably be preferable for larger systems.
        IF(tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals in symmetry class with needed spin to find allowed one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            ChosenUnocc=INT(NExcit*r)+1

!Now run through all allowed orbitals until we find this one.
            IF(tNoSymGenRandExcits) THEN
                nOrbs=nBasis/2      !No symmetry, therefore all orbitals of allowed spin possible to generate.
            ELSE
                Ind=ClassCountInd(SpinOrbB,SymB,MlB)
                nOrbs=SymLabelCounts2(2,Ind)
            ENDIF
            z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
            do i=0,nOrbs-1
                IF(tNoSymGenRandExcits) THEN
                    OrbB=(2*(i+1))-(SpinOrbB-1)
                ELSE
!Find the spin orbital index. SymLabelCounts has the index of the state for the given symmetry.
!                    OrbB=(2*SymLabelList(SymLabelCounts(1,SymB+1)+i))+SpinOrbB
!                    OrbB=SymLabelList2(SpinOrbB,(SymLabelCounts2(SpinOrbB,1,SymB+1)+i))
                    OrbB=SymLabelList2(SymLabelCounts2(1,Ind)+i)
                ENDIF

!Find out if the orbital is in the determinant, or is the other unocc picked
                IF((.not.(BTEST(ILUT((OrbB-1)/32),MOD((OrbB-1),32)))).and.(OrbB.ne.OrbA)) THEN
!The orbital is not found in the original determinant - increment z
                    z=z+1
                    IF(z.eq.ChosenUnocc) THEN
!We have got to the determinant that we want to pick.
                        EXIT
                    ENDIF
                ENDIF

            enddo

!We now have our final orbitals.
            IF(z.ne.ChosenUnocc) THEN
                WRITE(6,*) "This is a problem, since there should definitly be an allowed beta orbital once alpha is chosen..."
                CALL Stop_All("PickBOrb","Could not find allowed unoccupied orbital to excite to.")
            ENDIF

        ELSE
! METHOD 2 (Keep drawing orbitals from the desired symmetry and spin until we find one unoccupied)
! =========================

            IF(tNoSymGenRandExcits) THEN
                nOrbs=nBasis/2
            ELSE
!                nOrbs=SymLabelCounts2(SpinOrbB,2,SymB+1)
                Ind=ClassCountInd(SpinOrbB,SymB,MlB)
                nOrbs=SymLabelCounts2(2,Ind)
            ENDIF

            Attempts=0
            do while(.true.)
                
!Draw randomly from the set of orbitals
                IF(tMerTwist) THEN
                    CALL genrand_real2(r)
                ELSE
                    CALL RANLUX(r,1)
                ENDIF
                ChosenUnocc=INT(nOrbs*r)
                IF(tNoSymGenRandExcits) THEN
                    OrbB=(2*(ChosenUnocc+1))-(SpinOrbB-1)
                ELSE
!                    OrbB=(2*SymLabelList(SymLabelCounts(1,SymB+1)+ChosenUnocc))+SpinOrbB
!                    OrbB=SymLabelList2(SpinOrbB,SymLabelCounts2(SpinOrbB,1,SymB+1)+ChosenUnocc)
                     OrbB=SymLabelList2(SymLabelCounts2(1,Ind)+ChosenUnocc)
                ENDIF
!                WRITE(6,*) "B: ",OrbB, nOrbs, SpinOrbB
!                WRITE(6,*) "SymLabelList(1:nBasis): ",SymLabelList(1:nBasis)

!Find out if orbital is in nI or not. Accept if it isn't in it.
                IF((.not.(BTEST(ILUT((OrbB-1)/32),MOD((OrbB-1),32)))).and.(OrbB.ne.OrbA)) THEN
!Orbital not in nI. Accept.
                    EXIT
                ENDIF
                
                IF(Attempts.gt.250) THEN
                    WRITE(6,*) "Cannot find double excitation unoccupied orbital after 250 attempts..."
                    WRITE(6,*) "This is a problem, since there should definitly be an allowed beta orbital once alpha is chosen..."
                    CALL WRITEDET(6,nI,NEL,.true.)
                    CALL Stop_All("PickBOrb","Cannot find double excitation unoccupied orbital after 250 attempts...")
                ENDIF
                Attempts=Attempts+1

            enddo

        ENDIF

    END SUBROUTINE PickBOrb

!This routine does the same as the FindNumForbiddenOrbs routine, but is optimised for when there are no spatial symmetry considerations.    
    SUBROUTINE FindNumForbiddenOrbsNoSym(ForbiddenOrbs,ClassCountUnocc2,iSpn)
        INTEGER :: ClassCountUnocc2(ScratchSize)
        INTEGER :: ForbiddenOrbs,iSpn

!We know that all orbitals are totally symmetric, and that the symproduct=0

        ForbiddenOrbs=0
        IF(iSpn.eq.2) THEN
            IF(ClassCountUnocc2(1).eq.0) THEN
!There are no unoccupied alpha orbitals - are there any beta spins which are now forbidden?
                ForbiddenOrbs=ClassCountUnocc2(2)
            ENDIF
            IF(ClassCountUnocc2(2).eq.0) THEN
!There are no unoccupied beta orbitals - are there any alpha spins which are now forbidden?
                ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(1)
            ENDIF

        ELSEIF(iSpn.eq.1) THEN
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
            IF(ClassCountUnocc2(2).eq.1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                ForbiddenOrbs=1
            ENDIF
        ELSEIF(iSpn.eq.3) THEN
            IF(ClassCountUnocc2(1).eq.1) THEN
                ForbiddenOrbs=1
            ENDIF
        ENDIF

    END SUBROUTINE FindNumForbiddenOrbsNoSym

!This routine finds the number of orbitals which are allowed by spin, but not part of any spatial symmetry allowed unoccupied pairs.
!This number is needed for the correct normalisation of the probability of drawing any given A orbital since these can be chucked and redrawn.
    SUBROUTINE FindNumForbiddenOrbs(ForbiddenOrbs,ClassCountUnocc2,SymProduct,iSpn,SumMl)
        INTEGER :: ClassCountUnocc2(ScratchSize),OrbAMl,SumMl,j,k
        INTEGER :: ForbiddenOrbs,SymProduct,iSpn,i,ConjSym,Ind,l

        ForbiddenOrbs=0

        IF(tFixLz) THEN

            Ind=0
!Run over all possible b symmetries, and count the a orbitals which would be disallow due to the unavailability of a corresponding b orbital.
!            WRITE(6,*) "iSpn: ",iSpn,SymProduct,SumMl
            do k=-iMaxLz,iMaxLz
                OrbAMl=SumMl-k
!                WRITE(6,*) "OrbAMl: ",OrbAMl
                IF(abs(OrbAMl).le.iMaxLz) THEN
                    !If the OrbAMl which would be needed to require this B-symmetry is out of range, then there is no need to consider it - we cannot pick an A orb which would require this symmetry from the B orbital.
                    do i=0,nSymLabels-1
                        ConjSym=IEOR(SymProduct,i)
!                        WRITE(6,*) "ConjSym: ",ConjSym,i
                        do j=1,2
                            Ind=Ind+1
!                            WRITE(6,*) "Alpha/Beta: ",j
!                            WRITE(6,*) "***i***"
                            IF(ClassCountUnocc2(Ind).eq.0) THEN
!There is no point going in here if SymProduct=0
                                !Ignore if already spin-forbidden
                                IF(iSpn.eq.1) THEN
                                    IF(j.eq.1) THEN
                                        CYCLE  !We are only interested in beta orbitals
                                    ELSE
                                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(2,ConjSym,OrbAMl))
                                    ENDIF
                                ELSEIF(iSpn.eq.3) THEN
                                    IF(j.eq.2) THEN
                                        CYCLE   !We are only interested in alpha orbitals
                                    ELSE
                                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(1,ConjSym,OrbAMl))
!                                        IF(Iter.eq.13) WRITE(6,*) "Adding forbidden orb for sym ",i,ConjSym,ForbiddenOrbs
                                    ENDIF
                                ELSEIF(iSpn.eq.2) THEN  !alpha/beta pair - can forbid orbitals both ways.
                                    IF(j.eq.1) THEN
                                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(2,ConjSym,OrbAMl))
!                                        WRITE(6,*) "***",ForbiddenOrbs
                                    ELSE
                                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(1,ConjSym,OrbAMl))
!                                        WRITE(6,*) "***",ForbiddenOrbs
                                    ENDIF
                                ENDIF
                            ELSEIF((ClassCountUnocc2(Ind).eq.1).and.(iSpn.ne.2).and.(SymProduct.eq.0).and.(OrbAMl.eq.k)) THEN
                                !This is the situation where you actually need two orbitals of the given symmetry to allow this a orbital to be chosen.
                                IF(iSpn.eq.1.and.j.eq.1) CYCLE
                                IF(iSpn.eq.3.and.j.eq.2) CYCLE
                                ForbiddenOrbs=ForbiddenOrbs+1
!                                IF(Iter.eq.13) WRITE(6,*) "Extra forbidden orb for symmetry ",i,ForbiddenOrbs
                            ENDIF
                        enddo
                    enddo
                ELSE
                    !All unoccupied orbitals in this Ml block are forbbidden
                    IF(iSpn.eq.2) THEN
                        do l=Ind+1,Ind+nSymLabels*2
!                            WRITE(6,*) l
                            IF(ClassCountUnocc2(l).ne.0) THEN
                                ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(l)
                            ENDIF
                        enddo
                    ELSEIF(iSpn.eq.1) THEN
                    !Forbid all beta unoccupied orbitals (these are the second of the pair
                        do l=Ind+2,Ind+nSymLabels*2,2
                            IF(ClassCountUnocc2(l).ne.0) THEN
                                ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(l)
                            ENDIF
                        enddo
                    ELSE
                        do l=Ind+1,Ind+nSymLabels*2-1,2
                            IF(ClassCountUnocc2(l).ne.0) THEN
                                ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(l)
                            ENDIF
                        enddo
                    ENDIF

                    !Move onto the next k-block of B orbitals.
                    Ind=Ind+nSymLabels*2
                ENDIF
            enddo

        ELSE
            IF(iSpn.eq.2) THEN
!i,j are an alpha/beta pair. The number of forbidden orbitals includes all alphas and betas.

                do i=0,nSymLabels-1
!Run though all symmetries
                    IF(ClassCountUnocc2(ClassCountInd(1,i,0)).eq.0) THEN
!This symmetry has no unoccupied alpha orbitals - does its symmetry conjugate have any unoccupied beta orbitals which are now forbidden?
!If there are no unoccupied orbitals in this conjugate symmetry, then it won't increase the forbidden orbital number, since it can never be chosen.
                        ConjSym=IEOR(SymProduct,i)
                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(2,ConjSym,0)) !No unocc alphas in i, therefore all betas in ConjSym are forbidden
!                        WRITE(6,*) ClassCountUnocc2(2,ConjSym),i,ConjSym
                    ENDIF
                    IF(ClassCountUnocc2(ClassCountInd(2,i,0)).eq.0) THEN
!This symmetry has no unoccupied beta orbitals - does its symmetry conjugate have any unoccupied alpha orbitals which are now forbidden?
!If there are no unoccupied orbitals in this conjugate symmetry, then it won't increase the forbidden orbital number, since it can never be chosen.
                        ConjSym=IEOR(SymProduct,i)
                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(1,ConjSym,0))
!                        WRITE(6,*) ClassCountUnocc2(2,ConjSym),i,ConjSym
                    ENDIF
                enddo

            ELSEIF(iSpn.eq.1) THEN
                IF(SymProduct.ne.0) THEN
!i,j are a beta/beta pair. The number of forbidden orbitals is just betas
                    do i=0,nSymLabels-1
                        IF(ClassCountUnocc2(ClassCountInd(2,i,0)).eq.0) THEN
                            ConjSym=IEOR(SymProduct,i)
                            ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(2,ConjSym,0))
                        ENDIF
                    enddo
                ELSE
!There is a subtle point here, which could change the probabilities.
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
                    do i=0,nSymLabels-1
                        IF(ClassCountUnocc2(ClassCountInd(2,i,0)).eq.1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                            ForbiddenOrbs=ForbiddenOrbs+1
                        ENDIF
                    enddo
                ENDIF
            ELSEIF(iSpn.eq.3) THEN
                IF(SymProduct.ne.0) THEN
!i,j are a alpha/alpha pair. The number of forbidden orbitals is just alphas
                    do i=0,nSymLabels-1
                        IF(ClassCountUnocc2(ClassCountInd(1,i,0)).eq.0) THEN
                            ConjSym=IEOR(SymProduct,i)
                            ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(ClassCountInd(1,ConjSym,0))
                        ENDIF
                    enddo
                ELSE
!There is a subtle point here, which could change the probabilities.
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
                    do i=0,nSymLabels-1
                        IF(ClassCountUnocc2(ClassCountInd(1,i,0)).eq.1) THEN
!The one alpha orbital in this class is forbidden, since it cannot form a pair.
                            ForbiddenOrbs=ForbiddenOrbs+1
                        ENDIF
                    enddo
                ENDIF
            ENDIF
        ENDIF

    END SUBROUTINE FindNumForbiddenOrbs

!Choose the first unoccupied orbital to excite to. Spatial symmetry does not need to be considered, but spin symmetry does.
!After one has been chosen, we have to make sure that an allowed symmetry pair can be formed, otherwise the orbital is
!chucked and a new one drawn.
!Arguments:     NExcit = Number of possible spin-allowed unoccupied spinorbitals, including forbidden orbs (these will be chucked)
!               SpinOrbA = Spin of the chosen spin-orbital. 1 is alpha, -1 is beta.
!               OrbA = Index of the chosen spinorbital
!               SymB = Symmetry required of the second unoccupied spinorbital, so that sym(A) x sym(B) = SymProduct
!               SymProduct = (intent in), the symmetry of electron pair chosen
!               SumMl = Sum of the Ml values for the two picked electrons
!               MlA = Ml of the a orbital chosen
!               MlB = Required Ml of the b orbital still to be chosen
    SUBROUTINE PickAOrb(nI,iSpn,ILUT,ClassCountUnocc2,NExcit,Elec1Ind,Elec2Ind,SpinOrbA,OrbA,SymA,SymB,SymProduct,SumMl,MlA,MlB,ForbiddenOrbs,tAOrbFail)
        INTEGER :: nI(NEl),iSpn,Elec1Ind,Elec2Ind,SpinOrbA,AttemptsOverall,SymA,ForbiddenOrbs
        INTEGER :: NExcit,ChosenUnocc,z,i,OrbA,Attempts,SymB,SymProduct,SumMl,MlA,MlB
        INTEGER :: ILUT(0:NIfTot)
        INTEGER :: ClassCountUnocc2(ScratchSize)
        LOGICAL :: tAOrbFail
        REAL*8 :: r

        IF(iSpn.eq.2) THEN
!There is no restriction on whether we choose an alpha or beta spin, so there are nBasis-NEl possible spinorbitals to choose from.
!Therefore the spinOrbA variable has to be set after we have chosen one. For the other iSpn types, we can set it now.
            NExcit=nBasis-NEl
        ELSEIF(iSpn.eq.1) THEN 
!SpinOrbA indicates the spin of the chosen A orb. This is only really useful for iSpn=2
            SpinOrbA=-1  !Going to pick a beta orb
            NExcit=(nBasis/2)-nOccBeta  !This is the number of unoccupied beta spinorbitals there are.
        ELSE    !iSpn is 3 - alpha/alpha pair.
            SpinOrbA=1 !Going to pick an alpha orb
            NExcit=(nBasis/2)-nOccAlpha !This is the number of unoccupied alpha spinorbitals there are.
        ENDIF

        IF(NExcit.eq.ForbiddenOrbs) THEN
            tAOrbFail=.true.
            RETURN
        ELSE
            tAOrbFail=.false.
        ENDIF

        AttemptsOverall=0
        do while(.true.)
!Keep drawing unoccupied orbitals, until we find one which has an allowed partner to form a symmetry-allowed unoccupied pair.

            IF(iSpn.eq.2) THEN
!Electrons chosen were an alpha/beta pair, therefore first randomly chosen orbital can be an alpha OR beta orbital - no restriction.
            
                IF(tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals to find desired one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
                    IF(tMerTwist) THEN
                        CALL genrand_real2(r)
                    ELSE
                        CALL RANLUX(r,1)
                    ENDIF
                    ChosenUnocc=INT(NExcit*r)+1

!Now run through all allowed orbitals until we find this one.
                    z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
!This could be changed to an O[N] operation, rather than O[M] with a little thought...
                    do i=0,nBasis-1

!Find out if the orbital is in the determinant.
                        IF(.not.(BTEST(ILUT(i/32),MOD(i,32)))) THEN
!The orbital is not found in the original determinant - increment z
                            z=z+1
                            IF(z.eq.ChosenUnocc) THEN
!We have got to the determinant that we want to pick.
                                EXIT
                            ENDIF
                        ENDIF

                    enddo

!We now have our final orbital. The chosen unoccupied orbital is i+1.
                    IF(z.ne.ChosenUnocc) THEN
                        CALL Stop_All("PickAOrb","Could not find allowed unoccupied orbital to excite to.")
                    ENDIF

                    OrbA=i+1    !This is the allowed orbital

                ELSE
! METHOD 2 (Keep drawing orbitals randomly until we find one unoccupied). This should be more efficient, unless we have v. small basis sets.
! =========================

                    Attempts=0
                    do while(.true.)
                        
!Draw randomly from the set of orbitals
                        IF(tMerTwist) THEN
                            CALL genrand_real2(r)
                        ELSE
                            CALL RANLUX(r,1)
                        ENDIF
                        ChosenUnocc=INT(nBasis*r)

!Find out if orbital is in nI or not. Accept if it isn't in it.
                        IF(.not.(BTEST(ILUT(ChosenUnocc/32),MOD(ChosenUnocc,32)))) THEN
!Orbital not in nI. Accept.
                            EXIT
                        ENDIF
                        
                        IF(Attempts.gt.250) THEN
                            WRITE(6,*) "Cannot find A unoccupied orbital after 250 attempts..."
                            CALL WRITEDET(6,nI,NEL,.true.)
                            CALL Stop_All("PickAOrb","Cannot find A unoccupied orbital after 250 attempts...")
                        ENDIF
                        Attempts=Attempts+1

                    enddo

                    OrbA=ChosenUnocc+1  !This is the allowed orbital

                ENDIF

                SpinOrbA=G1(OrbA)%Ms

            ELSE
!We are either constrained to choose a beta orbital, or an alpha orbital - we know that there are NExcit of these.

                IF(tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals to find desired one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
                    IF(tMerTwist) THEN
                        CALL genrand_real2(r)
                    ELSE
                        CALL RANLUX(r,1)
                    ENDIF
                    ChosenUnocc=INT(NExcit*r)+1

!Now run through all allowed orbitals until we find this one.
!Alpha orbitals are the ODD numbered orbitals. Beta orbitals are EVEN numbered orbitals.

                    z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
!This could be changed to an O[N] operation, rather than O[M] with a little thought...
                    do i=1,nBasis/2
                        IF(iSpn.eq.1) THEN
!We want to run through all beta orbitals (odd numbered basis function)
                            OrbA=(2*i)-1
                        ELSE
!We want to run through all alpha orbitals (odd numbered basis functions)
                            OrbA=(2*i)
                        ENDIF

!Find out if the orbital is in the determinant.
                        IF(.not.(BTEST(ILUT((OrbA-1)/32),MOD((OrbA-1),32)))) THEN
!The orbital is not found in the original determinant - increment z
                            z=z+1
                            IF(z.eq.ChosenUnocc) THEN
!We have got to the determinant that we want to pick.
                                EXIT
                            ENDIF
                        ENDIF

                    enddo

!We now have our final orbital. The chosen unoccupied orbital is i+1.
                    IF(z.ne.ChosenUnocc) THEN
                        CALL Stop_All("PickAOrb","Could not find allowed unoccupied orbital to excite to.")
                    ENDIF

!OrbA is the allowed orbital

                ELSE
! METHOD 2 (Keep drawing orbitals randomly until we find one unoccupied). This should be more efficient, unless we have small basis sets.
! =========================

                    Attempts=0
                    do while(.true.)
                        
!Draw randomly from the set of orbitals
                        IF(tMerTwist) THEN
                            CALL genrand_real2(r)
                        ELSE
                            CALL RANLUX(r,1)
                        ENDIF
                        ChosenUnocc=INT((nBasis/2)*r)+1
                        IF(iSpn.eq.1) THEN
!We only want to choose beta orbitals(i.e. odd).
                            ChosenUnocc=(2*ChosenUnocc)-1
                        ELSE
!We only want to choose alpha orbitals(i.e. even).
                            ChosenUnocc=2*ChosenUnocc
                        ENDIF

!Find out if orbital is in nI or not. Accept if it isn't in it.
                        IF(.not.(BTEST(ILUT((ChosenUnocc-1)/32),MOD((ChosenUnocc-1),32)))) THEN
!Orbital not in nI. Accept.
                            EXIT
                        ENDIF
                        
                        IF(Attempts.gt.250) THEN
                            WRITE(6,*) "Cannot find A unoccupied orbital after 250 attempts..."
                            CALL WRITEDET(6,nI,NEL,.true.)
                            CALL Stop_All("PickAOrb","Cannot find A unoccupied orbital after 250 attempts...")
                        ENDIF
                        Attempts=Attempts+1

                    enddo

                    OrbA=ChosenUnocc  !This is the allowed orbital

                ENDIF

            ENDIF

!We now need to test whether this orbital has any symmetry-allowed unoccupied orbitals to form a pair with.
!To test this, we need to find the needed symmetry of B, in order that Sym(A) x Sym(B) = SymProduct
            IF(tNoSymGenRandExcits) THEN
                SymA=0
                SymB=0
                MlA=0
                MlB=0
            ELSE
                SymA=INT(G1(OrbA)%Sym%S,4)
                SymB=IEOR(SymA,SymProduct)
                MlB=0
                MlA=0
                IF(tFixLz) THEN
                    MlA=G1(OrbA)%Ml
                    MlB=SumMl-MlA
                ENDIF
            ENDIF

            IF(abs(MlB).le.iMaxLz) THEN
!Make sure that the B orbital that we would need to pick to conserve momentum is actually in the available range of Ml values.
                IF(iSpn.eq.2) THEN
!We want an alpha/beta unocc pair. 
                    IF(SpinOrbA.eq.1) THEN
!We have picked an alpha orbital - check to see if there are allowed beta unoccupied orbitals from the SymB Class.
                        IF(ClassCountUnocc2(ClassCountInd(2,SymB,MlB)).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                            EXIT
                        ENDIF
                    ELSE
!We have picked a beta orbital - check to see if there are allowed alpha unoccupied orbitals from the SymB Class.
                        IF(ClassCountUnocc2(ClassCountInd(1,SymB,MlB)).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                            EXIT
                        ENDIF
                    ENDIF
                ELSEIF(iSpn.eq.1) THEN
!We want a beta/beta pair.
                    IF((SymProduct.ne.0).or.(MlA.ne.MlB)) THEN
!Check to see if there are any unoccupied beta orbitals in the SymB Class.
                        IF(ClassCountUnocc2(ClassCountInd(2,SymB,MlB)).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                            EXIT
                        ENDIF
                    ELSE
!We want an orbital from the same class. Check that this isn't the only unoccupied beta orbital in the class.
                        IF(ClassCountUnocc2(ClassCountInd(2,SymB,MlB)).ne.1) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                            EXIT
                        ENDIF
                    ENDIF
                ELSE
!We want an alpha/alpha pair.
                    IF((SymProduct.ne.0).or.(MlA.ne.MlB)) THEN
!Check to see if there are any unoccupied alpha orbitals in the SymB Class.
                        IF(ClassCountUnocc2(ClassCountInd(1,SymB,MlB)).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                            EXIT
                        ENDIF
                    ELSE
!We want an orbital from the same class. Check that this isn't the only unoccupied alpha orbital in the class.
                        IF(ClassCountUnocc2(ClassCountInd(1,SymB,MlB)).ne.1) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                            EXIT
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF

            IF(AttemptsOverall.gt.1000) THEN
                WRITE(6,*) "***",NExcit,ForbiddenOrbs
                WRITE(6,*) "Cannot find first allowed unoccupied orbital for given i,j pair after 1000 attempts."
                WRITE(6,*) "It may be that there are no possible excitations from this i,j pair, in which case "
                WRITE(6,*) "the given algorithm is inadequate to describe excitations from such a small space."
                WRITE(6,*) "Try reverting to old excitation generators."
                CALL WRITEDET(6,nI,NEl,.true.)
                WRITE(6,*) "I,J pair; sym_i, sym_j: ",nI(Elec1Ind),nI(Elec2Ind),G1(nI(Elec1Ind))%Sym%S,G1(nI(Elec2Ind))%Sym%S
                CALL Stop_All("PickAOrb","Cannot find first allowed unocc orb for double excitation")
            ENDIF
            AttemptsOverall=AttemptsOverall+1

        enddo

    END SUBROUTINE PickAOrb

!This routine takes determinant nI and returns two randomly chosen electrons, whose index in nI is Elec1Ind and Elec2Ind.
!These electrons have symmetry product SymProduct and spin pairing iSpn, where iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha.
!If IndInp = -1, the pair is picked randomly with prob = 1/ElecPairs. Otherwise, it will choose electron pair given by index IndInp.
    SUBROUTINE PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,SumMl,IndInp)
        INTEGER :: Ind,X,K,Elec1Ind,Elec2Ind,SymProduct,IndInp
        INTEGER :: nI(NEl),iSpn,SumMl
        REAL*8 :: r
!Triangular indexing system.
!This is used for picking two distinct electrons out of all N(N-1)/2 pairs.
!
!   12  13  14  15          1   2   3   4
!       23  24  25              5   6   7
!           34  35      =>          8   9
!               45                      10

!For a given index, ind, there are [N(N-1)/2 - ind] elements at positions larger than ind. Call this number X.
!We want to find out how many rows there are after the row containing the element ind. t rows has t(t+1)/2 elements in it.
!Therefore, to find out the number of rows after ind, we want to find the largest K, such that K(K+1)/2 =< X
!The solution to this is that K =< (SQRT(8*X+1)-1)/2, therefore we can find K (the largest integer which satisfies this).
!We then know the number of rows after the element ind. Therefore, since there are N-1 rows in total, we know we are on row N-1-K.
!This gives us the index of the first electron.
!To find the second electron (i.e. the column), we know that out of the X elements in positions larger than ind, K(K+1)/2 are in the next rows.
!This means that X - K(K+1)/2 are in the same row. There are N-(N-1-K) = 1+K elements in the row chosen, and so the number of elements into the 
!row it is, is given by (1+K) - (X-K(K+1)/2). However, in row z, the first column index is z+1. Therefore the index of the second electron is
!(1+K) - (X-K(K+1)/2) + N-K-1 = N-X+(K(K+1)/2).

!        ElecPairs=(NEl*(NEl-1))/2

! If we want to find an index randomly, IndInp will be -1.
        IF(IndInp.eq.-1) THEN
!Find an index randomly.
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            Ind=INT(ElecPairs*r)+1
        ELSE
!Otherwise we are looking for a specific electron pair specified by IndInp            
            Ind=IndInp
        ENDIF

!X is number of elements at positions larger than ind
        X=ElecPairs-Ind
!K is the number of complete rows after the element ind
        K=INT((SQRT(8.D0*REAL(X,r2)+1.D0)-1.D0)/2.D0)
        Elec1Ind=NEl-1-K
        Elec2Ind=NEl-X+((K*(K+1))/2)

!We now want to find the symmetry product of the two electrons, and the spin product of the two electrons.
        IF(tNoSymGenRandExcits) THEN
            SymProduct=0
            SumMl=0
        ELSE
            SymProduct=INT(IEOR(G1(nI(Elec1Ind))%Sym%S,G1(nI(Elec2Ind))%Sym%S),4)
            IF(tFixLz) THEN
                SumMl=G1(nI(Elec1Ind))%Ml+G1(nI(Elec2Ind))%Ml
            ENDIF
        ENDIF

        IF((G1(nI(Elec1Ind))%Ms)*(G1(nI(Elec2Ind))%Ms).eq.-1) THEN
!We have an alpha beta pair of electrons.
            iSpn=2
        ELSE
            IF(G1(nI(Elec1Ind))%Ms.eq.1) THEN
!We have an alpha alpha pair of electrons.
                iSpn=3
            ELSE
!We have a beta beta pair of electrons.
                iSpn=1
            ENDIF
        ENDIF

    END SUBROUTINE PickElecPair

    SUBROUTINE CheckIfSingleExcits(ElecsWNoExcits,ClassCount2,ClassCountUnocc2)
        INTEGER :: ElecsWNoExcits,ClassCount2(ScratchSize),ClassCountUnocc2(ScratchSize),i,k
        INTEGER :: Ind1,Ind2


!First, we need to find out if there are any electrons which have no possible excitations. This is because these will need to be redrawn and so 
!will affect the probabilities.
        ElecsWNoExcits=0

        IF(.not.tNoSymGenRandExcits) THEN
            IF(tFixLz) THEN
!Here, we also have to check that the electron is momentum allowed.
                do i=1,ScratchSize
!Run through all labels
                    IF((ClassCount2(i).ne.0).and.(ClassCountUnocc2(i).eq.0)) THEN
!If there are electrons in this class with no possible unoccupied orbitals in the same class, these electrons have no single excitations.
                        ElecsWNoExcits=ElecsWNoExcits+ClassCount2(i)
                    ENDIF
                enddo
                
!                do k=-iMaxLz,iMaxLz,1
!                    Ind1=ClassCountInd(1,0,k)
!                    Ind2=ClassCountInd(1,0,-k)
!                    do i=0,nSymLabels*2-1
!                        IF((ClassCount2(i+Ind1).ne.0).and.(ClassCountUnocc2(i+Ind2).eq.0)) THEN
!                            ElecsWNoExcits=ElecsWNoExcits+ClassCount2(i+Ind1)
!                        ENDIF
!                    enddo
!                enddo
            ELSE

                do i=1,ScratchSize
!Run through all labels
                    IF((ClassCount2(i).ne.0).and.(ClassCountUnocc2(i).eq.0)) THEN
!If there are electrons in this class with no possible unoccupied orbitals in the same class, these electrons have no single excitations.
                        ElecsWNoExcits=ElecsWNoExcits+ClassCount2(i)
                    ENDIF
                enddo
            ENDIF

!Need to look for forbidden electrons through all the irreps.

            IF(ElecsWNoExcits.eq.NEl) THEN
!There are no single excitations from this determinant at all. This means the probability to create a double excitation = 1
!Then we will create a double excitation instead.
                pDoubNew=1.D0
                RETURN
            ENDIF

        ENDIF

    END SUBROUTINE CheckIfSingleExcits
        

    SUBROUTINE CreateSingleExcit(nI,nJ,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen,ElecsWNoExcits)
        INTEGER :: ElecsWNoExcits,i,Attempts,nOrbs,z,Orb
        INTEGER :: Eleci,ElecSym,nI(NEl),nJ(NEl),NExcit,iSpn,ChosenUnocc
        INTEGER :: ExcitMat(2,2),ExcitLevel,iGetExcitLevel
        INTEGER :: ClassCount2(ScratchSize)
        INTEGER :: ClassCountUnocc2(ScratchSize),k,ElecK,Ind,SymIndex
        INTEGER :: ILUT(0:NIfTot),Ind1,Ind2
        REAL*8 :: r,pGen
        LOGICAL :: tParity,IsValidDet,SymAllowed


!This will not normally be called.
        IF((.not.tNoSingsPossible).and.(.not.tNoSymGenRandExcits)) THEN
!First, we need to find out if there are any electrons which have no possible excitations. This is because these will need to be redrawn and so 
!will affect the probabilities.
            ElecsWNoExcits=0
!Need to look for forbidden electrons through all the irreps.

            do i=1,ScratchSize
!Run through all labels
                IF((ClassCount2(i).ne.0).and.(ClassCountUnocc2(i).eq.0)) THEN
!If there are electrons in this class with no possible unoccupied orbitals in the same class, these electrons have no single excitations.
                    ElecsWNoExcits=ElecsWNoExcits+ClassCount2(i)
                ENDIF
            enddo

        ENDIF

!        IF(ElecsWNoExcits.eq.NEl) THEN
!!There are no single excitations from this determinant at all. Indicate this by putting pGen=-1.D0.
!!Then we will create a double excitation instead.
!            pGen=-1.D0
!            RETURN
!        ENDIF


        Attempts=0
        do while(.true.)


!Choose an electron randomly...
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            Eleci=INT(NEl*r)+1

!Find symmetry of chosen electron
            IF(tNoSymGenRandExcits) THEN
                ElecSym=0
            ELSE
                ElecSym=INT((G1(nI(Eleci))%Sym%S),4)
                ElecK=G1(nI(Eleci))%Ml
            ENDIF

            IF(G1(nI(Eleci))%Ms.eq.1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                iSpn=1
            ELSE
!Beta orbital
                iSpn=2
            ENDIF
            SymIndex=ClassCountInd(iSpn,ElecSym,ElecK)
            NExcit=ClassCountUnocc2(SymIndex)

            IF(NExcit.ne.0) EXIT    !Have found electron with allowed excitations

            IF(Attempts.gt.250) THEN
                WRITE(6,*) "Cannot find single excitation from electrons after 250 attempts..."
                CALL WRITEDET(6,nI,NEL,.true.)
!                WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
!                WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
                WRITE(6,*) "***"
!                WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
!                WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
                CALL Stop_All("CreateSingleExcit","Cannot find single excitation from electrons after 250 attempts...")
            ENDIF
            Attempts=Attempts+1

        enddo



!Now we need to choose the unoccupied orbital for the chosen electron.
!There are two ways to do this. We can either choose the orbital we want out of the NExcit possible unoccupied orbitals.
!It would then be necessary to cycle through all orbitals of that symmetry and spin, only counting the unoccupied ones to find the correct determinant.
!Alternatively, we could pick orbitals at random and redraw until we find an allowed one. This would probably be preferable for larger systems.
        IF(tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals in symmetry class with needed spin to find allowed one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            ChosenUnocc=INT(NExcit*r)+1

!Now run through all allowed orbitals until we find this one.
            IF(tNoSymGenRandExcits) THEN
                nOrbs=nBasis/2
            ELSE
!                nOrbs=SymLabelCounts(2,ElecSym+1)
!                nOrbs=SymLabelCounts2(iSpn,2,ElecSym+1)
                nOrbs=OrbClassCount(SymIndex)

!                !!REMOVE THIS TEST ONCE WORKING!!
!                IF(nOrbs.ne.SymLabelCounts2(2,SymIndex)) THEN
!                    CALL Stop_All("GetSingleExcit","Error in symmetry arrays")
!                ENDIF
            ENDIF

            z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
            do i=0,nOrbs-1
!Find the spin orbital index. SymLabelCounts has the index of the state for the given symmetry.
                IF(tNoSymGenRandExcits) THEN
                    Orb=(2*(i+1))-(iSpn-1)
                ELSE
!                    Orb=(2*SymLabelList(SymLabelCounts(1,ElecSym+1)+i))-(iSpn-1)
!                    Orb=SymLabelList2(iSpn,SymLabelCounts2(iSpn,1,ElecSym+1)+i)
                    Orb=SymLabelList2(SymLabelCounts2(1,SymIndex)+i)
                ENDIF

!Find out if the orbital is in the determinant.
                IF(.not.(BTEST(ILUT((Orb-1)/32),MOD(Orb-1,32)))) THEN
!The orbital is not found in the original determinant - increment z
                    z=z+1
                    IF(z.eq.ChosenUnocc) THEN
!We have got to the determinant that we want to pick.
                        EXIT
                    ENDIF
                ENDIF

            enddo

!We now have our final orbitals. i=nI(Eleci). a=Orb.
            IF(z.ne.ChosenUnocc) THEN
                CALL Stop_All("CreateSingleExcit","Could not find allowed unoccupied orbital to excite to.")
            ENDIF

        ELSE
! METHOD 2 (Keep drawing orbitals from the desired symmetry and spin until we find one unoccupied)
! =========================

            IF(tNoSymGenRandExcits) THEN
                nOrbs=nBasis/2
            ELSE
                nOrbs=OrbClassCount(SymIndex)

!                !!REMOVE THIS TEST ONCE WORKING!!
!                IF(nOrbs.ne.SymLabelCounts2(2,SymIndex)) THEN
!                    CALL Stop_All("GetSingleExcit","Error in symmetry arrays")
!                ENDIF
!                nOrbs=SymLabelCounts2(iSpn,2,ElecSym+1)
            ENDIF
            Attempts=0
            do while(.true.)
                
!Draw randomly from the set of orbitals
                IF(tMerTwist) THEN
                    CALL genrand_real2(r)
                ELSE
                    CALL RANLUX(r,1)
                ENDIF
                ChosenUnocc=INT(nOrbs*r)
                IF(tNoSymGenRandExcits) THEN
                    Orb=(2*(ChosenUnocc+1))-(iSpn-1)
                ELSE
!                    Orb=(2*SymLabelList(SymLabelCounts(1,ElecSym+1)+ChosenUnocc))-(iSpn-1)
                    Orb=SymLabelList2(SymLabelCounts2(1,SymIndex)+ChosenUnocc)
                ENDIF

!Find out if orbital is in nI or not. Accept if it isn't in it.
                IF(.not.(BTEST(ILUT((Orb-1)/32),MOD(Orb-1,32)))) THEN
!Orbital not in nI. Accept.
                    EXIT
                ENDIF
                
                IF(Attempts.gt.250) THEN
                    WRITE(6,*) "Cannot find single excitation unoccupied orbital after 250 attempts..."
                    WRITE(6,*) "Desired symmetry of unoccupied orbital = ",ElecSym
                    WRITE(6,*) "Number of orbitals (of correct spin) in symmetry = ",nOrbs
                    WRITE(6,*) "Number of orbitals to legitimatly pick = ",NExcit
                    CALL WRITEDET(6,nI,NEL,.true.)
!                    WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
!                    WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
!                    WRITE(6,*) "***"
!                    WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
!                    WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
                    CALL Stop_All("CreateSingleExcit","Cannot find single excitation unoccupied orbital after 250 attempts...")
                ENDIF
                Attempts=Attempts+1

            enddo

        ENDIF

        nJ(:)=nI(:)
!ExcitMat wants to be the index in nI of the orbital to excite from, but returns the actual orbitals.
        ExcitMat(1,1)=Eleci
        ExcitMat(2,1)=Orb
        CALL FindExcitDet(ExcitMat,nJ,1,TParity)

!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
!        CALL IsSymAllowedExcit(nI,nJ,1,ExcitMat,SymAllowed)

!Now we need to find the probability of creating this excitation.
!This is: P_single x P(i) x P(a|i) x N/(N-ElecsWNoExcits)
!        pGen=(1.D0-pDoubNew)*(1.D0/real(NEl,r2))*(1.D0/real(NExcit,r2))*((real(NEl,r2))/(real((NEl-ElecsWNoExcits),r2)))
        pGen=(1.D0-pDoubNew)/(REAL((NExcit*(NEl-ElecsWNoExcits)),r2))

!        WRITE(6,*) "ElecsWNoExcits: ",ElecsWNoExcits

    END SUBROUTINE CreateSingleExcit

!This routine returns two arrays of length ScratchSize, which have information about the number of orbitals, occupied and unoccupied respectively, of each symmetry.
    SUBROUTINE ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)
        INTEGER :: i,nI(NEl)
        INTEGER :: ClassCount2(ScratchSize)
        INTEGER :: ClassCountUnocc2(ScratchSize)
!        INTEGER :: Alph,Bet

!        Alph=0
!        Bet=0
        ClassCount2(:)=0
        ClassCountUnocc2(:)=OrbClassCount(:)
!nOccAlpha and nOccBeta now set in the system block. Since we conserve Sz, these will not change.
!        NOccAlpha=0
!        NOccBeta=0

        IF(tNoSymGenRandExcits) THEN

!            do i=1,NEl
!Create ILUT for O[1] comparison of orbitals in root determinant - This is now read in
!                ILUT((nI(I)-1)/32)=IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32))
!            enddo
!All orbitals are in irrep 0
            ClassCount2(ClassCountInd(1,0,0))=nOccAlpha
            ClassCount2(ClassCountInd(2,0,0))=nOccBeta
            ClassCountUnocc2(ClassCountInd(1,0,0))=ClassCountUnocc2(ClassCountInd(1,0,0))-nOccAlpha
            ClassCountUnocc2(ClassCountInd(2,0,0))=ClassCountUnocc2(ClassCountInd(2,0,0))-nOccBeta
        
        ELSE

            do i=1,NEl
                
!Create ILUT for O[1] comparison of orbitals in root determinant - This is now read in
!                ILUT((nI(I)-1)/32)=IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32))

                IF(G1(nI(i))%Ms.eq.1) THEN
!orbital is an alpha orbital and symmetry of the orbital can be found in G1
!                    WRITE(6,*) G1(nI(i))%Ms,G1(nI(i))%Sym%S
                    ClassCount2(ClassCountInd(1,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))=ClassCount2(ClassCountInd(1,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))+1
                    ClassCountUnocc2(ClassCountInd(1,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))=ClassCountUnocc2(ClassCountInd(1,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))-1
!                    Alph=Alph+1
!                    NOccAlpha=NOccAlpha+1

                ELSE
!orbital is a beta orbital
!                    WRITE(6,*) G1(nI(i))%Ms,G1(nI(i))%Sym%S
                    ClassCount2(ClassCountInd(2,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))=ClassCount2(ClassCountInd(2,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))+1
                    ClassCountUnocc2(ClassCountInd(2,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))=ClassCountUnocc2(ClassCountInd(2,INT(G1(nI(i))%Sym%S,4),G1(nI(i))%Ml))-1
!                    Bet=Bet+1
!                    NOccBeta=NOccBeta+1
                ENDIF
            enddo

!We now want to find ClassCountUnocc2 - the unoccupied version of the array
!SymLabelCounts(2,Sym) gives the number of *spin-orbitals* in each symmetry class.

!We don't need to do this any more, since we store the whole space classcounts in OrbClassCount initially.
!            do i=1,nSymLabels
!                ClassCountUnocc2(ClassCountInd(1,i-1,0))=SymLabelCounts2(1,2,i)-ClassCount2(ClassCountInd(1,i-1,0))
!                ClassCountUnocc2(ClassCountInd(2,i-1,0))=SymLabelCounts2(2,2,i)-ClassCount2(ClassCountInd(2,i-1,0))
!            enddo

!            WRITE(6,*) "Alph=",alph,"Bet=",Bet

        ENDIF

!        WRITE(6,*) nI(:)
!        WRITE(6,*) "***"
!        WRITE(6,*) ScratchSize
!        WRITE(6,*) "***"
!        WRITE(6,*) ClassCount2(:)

    END SUBROUTINE ConstructClassCounts



!This routine will calculate the PGen between two connected determinants, nI and nJ which are IC excitations of each other, using the unbiased scheme.
!Only the excitation matrix is needed (1,*) are the i,j orbs, and (2,*) are the a,b orbs.
!This is the prob of generating nJ FROM nI, not the other way round.
!Passed in is also the ClassCount2 arrays for nI, and the probability of picking a double.
!A word of warning: The routine does not check that the determinants are indeed connected, and may well return a non-zero probability even if they arent.
!Therefore, make sure that they are at most double excitations of each other.
    SUBROUTINE CalcNonUniPGen(Ex,IC,ClassCount2,ClassCountUnocc2,pDoub,pGen)
        REAL*8 :: pDoub,pGen!,PabGivenij
        INTEGER :: ClassCount2(ScratchSize),ForbiddenOrbs,SymA,SymB,SumMl,MlA,MlB,k,Elec1Ml
        INTEGER :: ClassCountUnocc2(ScratchSize),ElecsWNoExcits,i,NExcitOtherWay,Ind1,Ind2
        INTEGER :: SymProduct,OrbI,OrbJ,iSpn,NExcitA,NExcitB,IC,ElecSym,OrbA,OrbB,Ex(2,2)
            
        pDoubNew=pDoub
        IF(tNoSingsPossible) THEN
!This will check if there are any possible single excitations from this determinant
!If there are not, then this will change pDoubNew so that it = 1 and only doubles will be generated.
            CALL CheckIfSingleExcits(ElecsWNoExcits,ClassCount2,ClassCountUnocc2)
        ENDIF

        IF(IC.eq.1) THEN

!First, we need to find out if there are any electrons which have no possible excitations. This is because these will need to be redrawn and so 
!will affect the probabilities.

            IF(.not.tNoSingsPossible) THEN

                ElecsWNoExcits=0
                do i=1,ScratchSize
!Run through all labels
                    IF((ClassCount2(i).ne.0).and.(ClassCountUnocc2(i).eq.0)) THEN
!If there are electrons in this class with no possible unoccupied orbitals in the same class, these electrons have no single excitations.
                        ElecsWNoExcits=ElecsWNoExcits+ClassCount2(i)
                    ENDIF
                enddo
            ENDIF

            IF(tNoSymGenRandExcits) THEN
!Find symmetry of chosen electron
                ElecSym=0
                Elec1Ml=0
            ELSE
                ElecSym=INT((G1(Ex(1,1))%Sym%S),4)
                IF(tFixLz) THEN
                    Elec1Ml=G1(Ex(1,1))%Ml
                ELSE
                    Elec1Ml=0
                ENDIF
            ENDIF


            IF(G1(Ex(1,1))%Ms.eq.1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                NExcitA=ClassCountUnocc2(ClassCountInd(1,ElecSym,Elec1Ml))
            ELSE
!Beta orbital
                NExcitA=ClassCountUnocc2(ClassCountInd(2,ElecSym,Elec1Ml))
            ENDIF

!Now we need to find the probability of creating this excitation.
!This is: P_single x P(i) x P(a|i) x N/(N-ElecsWNoExcits)
!Prob of generating a single is 1-pDoub
!            pGen=(1.D0-pDoub)*(1.D0/real(NEl,r2))*(1.D0/real(NExcitA,r2))*((real(NEl,r2))/(real((NEl-ElecsWNoExcits),r2)))
            pGen=(1.D0-pDoubNew)/(REAL((NExcitA*(NEl-ElecsWNoExcits)),r2))

        ELSE
!Prob of generating a double excitation.
!Find the I and J orbitals
            OrbI=Ex(1,1)
            OrbJ=Ex(1,2)
            OrbA=Ex(2,1)
            OrbB=Ex(2,2)
!            WRITE(6,*) OrbI,OrbJ,OrbA,OrbB

!Find the spin-product of the occupied pair
            IF((G1(OrbI)%Ms)*(G1(OrbJ)%Ms).eq.-1) THEN
!We have an alpha beta pair of electrons.
                iSpn=2
!NExcit is the number of allowed unoccupied orbitals to choose
                NExcitA=nBasis-NEl
            ELSE
                IF(G1(OrbI)%Ms.eq.1) THEN
!We have an alpha alpha pair of electrons.
                    iSpn=3
                    NExcitA=(nBasis/2)-nOccAlpha     !This is the number of unocc alpha spinorbs
                ELSE
!We have a beta beta pair of electrons.
                    iSpn=1
                    NExcitA=(nBasis/2)-nOccBeta      !This is the number of unocc beta spinorbs
                ENDIF
            ENDIF
        
            IF(tFixLz) THEN
                SumMl=G1(OrbI)%Ml+G1(OrbJ)%Ml
                MlA=G1(OrbA)%Ml
                MlB=SumMl-MlA
            ELSE
                SumMl=0
                MlA=0
                MlB=0
            ENDIF
            IF(tNoSymGenRandExcits) THEN
                ElecSym=0
                CALL FindNumForbiddenOrbsNoSym(ForbiddenOrbs,ClassCountUnocc2,iSpn)
                SymA=0
                SymB=0
            ELSE
!Calculate the symmetry product of the occupied orbital pair
                ElecSym=INT(IEOR(G1(OrbI)%Sym%S,G1(OrbJ)%Sym%S),4)
!This will calculate the A orbitals which will have no B pair
                CALL FindNumForbiddenOrbs(ForbiddenOrbs,ClassCountUnocc2,ElecSym,iSpn,SumMl)
!Need to find the symmetries of the unoccupied A and B orbitals.
                SymA=INT(G1(OrbA)%Sym%S,4)
                SymB=IEOR(SymA,ElecSym)
            ENDIF
!            WRITE(6,*) "Check: ",ForbiddenOrbs,ElecSym,iSpn,SymA,SymB,OrbA,OrbB


!We want to calculate the number of possible B's given the symmetry and spin it has to be since we have already picked A.
!We have calculated in NExcit the number of orbitals available for B given A, but we also need to know the number of orbitals to choose from for A IF
!we had picked B first.
            IF(iSpn.eq.2) THEN
!If iSpn=2, then we want to find a spinorbital of the opposite spin of SpinOrbA
                IF((G1(OrbA)%Ms).eq.-1) THEN
!We have already picked a beta orbital, so now we want to pick an alpha orbital. Find out how many of these there are.
                    NExcitB=ClassCountUnocc2(ClassCountInd(1,SymB,MlB))
                    NExcitOtherWay=ClassCountUnocc2(ClassCountInd(2,SymA,MlA))
                ELSE
!Want to pick an beta orbital.
                    NExcitB=ClassCountUnocc2(ClassCountInd(2,SymB,MlB))
                    NExcitOtherWay=ClassCountUnocc2(ClassCountInd(1,SymA,MlA))
                ENDIF
            ELSEIF(iSpn.eq.1) THEN
!Definitely want a beta orbital
                NExcitB=ClassCountUnocc2(ClassCountInd(2,SymB,MlB))
                NExcitOtherWay=ClassCountUnocc2(ClassCountInd(2,SymA,MlA))
            ELSE
!Definitely want an alpha orbital
                NExcitB=ClassCountUnocc2(ClassCountInd(1,SymB,MlB))
                NExcitOtherWay=ClassCountUnocc2(ClassCountInd(1,SymA,MlA))
            ENDIF

            IF((iSpn.ne.2).and.(ElecSym.eq.0).and.(MlB.eq.MlA)) THEN
!In this case, we need to check that we do not pick the same orbital as OrbA. If we do this, then we need to redraw.
!Only when ElecSym=0 will the classes of a and b be the same, and the spins will be different if iSpn=2, so this is the only possibility of a clash.
                NExcitB=NExcitB-1     !Subtract 1 from the number of possible orbitals since we cannot choose orbital A.
                NExcitOtherWay=NExcitOtherWay-1     !The same goes for the probabilities the other way round.
            ENDIF

!            PabGivenij=(1.D0/real((NExcitA-ForbiddenOrbs),r2))*((1.D0/real(NExcitB,r2))+(1.D0/real(NExcitOtherWay,r2)))
!            pGen=pDoub*(1.D0/real(ElecPairs,r2))*PabGivenij
            pGen=pDoubNew*((1.D0/real(NExcitB,r2))+(1.D0/real(NExcitOtherWay,r2)))/(REAL((ElecPairs*(NExcitA-ForbiddenOrbs)),r2))
!            WRITE(6,*) "***",pDoubNew,NExcitB,NExcitOtherWay,ElecPairs,NExcitA,ForbiddenOrbs,ElecSym,iSpn,SymA,SymB

        ENDIF

    END SUBROUTINE CalcNonUniPGen

!ClassCount arrays are not going to be 1D arrays, so that new symmetries can be added more easily.
!This means that an indexing system is needed for the array.
!For spin, alpha=1, beta=2. Sym goes from 0:nSymLabels-1 and Mom goes from -Lmax to Lmax.
    INTEGER FUNCTION ClassCountInd(Spin,Sym,Mom)
        INTEGER :: Spin,Sym,Mom

        IF(tFixLz) THEN
            ClassCountInd=2*nSymLabels*(Mom+iMaxLz)+2*Sym+Spin
        ELSE
            ClassCountInd=2*Sym+Spin
        ENDIF

        IF(tNoSymGenRandExcits) THEN
            IF(Spin.eq.1) THEN
                ClassCountInd=1
            ELSE
                ClassCountInd=2
            ENDIF
        ENDIF

    END FUNCTION ClassCountInd




    !*********  BIASED EXCITATION GENERATION ROUTINES **************!
    


    
!This routine is an importance sampled excitation generator. However, it is currently set up to work with the
!spawning algorithm, since a stochastic choice as to whether the particle is accepted or not is also done within the routine.
!Because of this, tau is needed for the timestep of the simulation, and iCreate is returned as the number of children to create
!on the determinant. If this is zero, then no childred are to be created.
    SUBROUTINE GenRandSymExcitBiased(nI,iLut,nJ,pDoub,IC,ExcitMat,TParity,exFlag,nParts,WSign,tau,iCreate)
        INTEGER :: nI(NEl),nJ(NEl),IC,ExcitMat(2,2),Attempts,exFlag
        INTEGER :: ILUT(0:NIfTot),i,iCreate,nParts,WSign,ElecsWNoExcits
        LOGICAL :: tNoSuccess,tParity
        REAL*8 :: pDoub,pGen,r,tau
        CHARACTER , PARAMETER :: this_routine='GenRandSymExcitBiased'

        IF(.not.TwoCycleSymGens) THEN
!Currently only available for molecular systems, or without using symmetry.
            IF(.not.tNoSymGenRandExcits) THEN
                WRITE(6,*) "GenRandSymExcitBiased can only be used for molecular systems"
                WRITE(6,*) "This is because of difficulties with other symmetries setup."
                WRITE(6,*) "If you want to use these excitation generators, then add NOSYMGEN to the input to ignore symmetry while generating excitations."
                CALL FLUSH(6)
                CALL Stop_All(this_routine,"GenRandSymExcitBiased can only be used for molecular systems using symmetry")
            ENDIF
        ELSEIF(nBasisMax(2,3).eq.1) THEN
            CALL Stop_All(this_routine,"GenRandSymExcitBiased can not be used with UHF systems currently")
        ELSEIF(tFixLz) THEN
            CALL Stop_All(this_routine,"GenRandSymExcitBiased can not be used with tFixLz currently")
        ENDIF
!        MaxABPairs=(nBasis*(nBasis-1)/2)

!ExFlag is 1 for singles, 2 for just doubles, and 3 for both.
        IF(ExFlag.eq.3) THEN
!Choose whether to generate a double or single excitation. Prob of generating a double is given by pDoub.
            pDoubNew=pDoub
            IF(pDoubNew.gt.1.D0) CALL Stop_All(this_routine,"pDoub is greater than 1")

            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            IF(r.lt.pDoubNew) THEN
!A double excitation has been chosen to be created.
                IC=2
            ELSE
                IC=1
            ENDIF
        ELSEIF(ExFlag.eq.2) THEN
            IC=2
            pDoubNew=1.D0
        ELSEIF(ExFlag.eq.1) THEN
            IC=1
            pDoubNew=0.D0
        ELSE
            CALL Stop_All(this_routine,"Error in choosing excitations to create.")
        ENDIF

        IF(IC.eq.2) THEN
            CALL CreateDoubExcitBiased(nI,nJ,ILUT,ExcitMat,tParity,nParts,WSign,Tau,iCreate)
        ELSE
            CALL CreateSingleExcitBiased(nI,nJ,iLut,ExcitMat,tParity,ElecsWNoExcits,nParts,WSign,Tau,iCreate)
            IF(ElecsWNoExcits.eq.NEl) THEN
!                IF(ExFlag.ne.3) THEN
!Should not be changing pDoub since it won't affect all doubles generation equally.
                CALL Stop_All(this_routine,"Found determinant with no singles. Turn symmetry off, or fix code (see non-biased code).")
!                ENDIF
!                pDoubNew=1.D0
!                IC=2
!                CALL CreateDoubExcitBiased(nI,nJ,ILUT,ExcitMat,tParity,nParts,WSign,Tau,iCreate)
            ENDIF

        ENDIF

    END SUBROUTINE GenRandSymExcitBiased

    SUBROUTINE CreateSingleExcitBiased(nI,nJ,iLut,ExcitMat,tParity,ElecsWNoExcits,nParts,WSign,Tau,iCreate)
        Use SystemData, only: FCoul
        INTEGER :: ClassCount2(ScratchSize),i,Attempts,OrbA
        INTEGER :: ClassCountUnocc2(ScratchSize)
        INTEGER :: ElecsWNoExcits,nParts,WSign,iCreate,nI(NEl),nJ(NEl),iLut(0:NIfTot),Ind
        INTEGER :: ExcitMat(2,2),SpawnOrb(nBasis),Eleci,ElecSym,NExcit,VecInd,ispn,EndSymState,j
        REAL*8 :: Tau,SpawnProb(nBasis),NormProb,r,rat
        LOGICAL :: tParity,SymAllowed
        TYPE(HElement) :: rh

!First, we need to do an O[N] operation to find the number of occupied alpha electrons, number of occupied beta electrons
!and number of occupied electrons of each symmetry class and spin. This is similar to the ClassCount array.
!This has the format (Spn,sym), where Spin=1,2 corresponding to alpha and beta.
!For molecular systems, sym runs from 0 to 7. This is NOT general and should be made so using SymLabels.
!This could be stored to save doing this multiple times, but shouldn't be too costly an operation.
        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)

!We need to find out if there are any electrons which have no possible excitations. This is because these will need to be redrawn and so 
!will affect the probabilities.
        ElecsWNoExcits=0

!Need to look for forbidden electrons through all the irreps.
        do i=0,nSymLabels-1
!Run through all labels
            IF((ClassCount2(ClassCountInd(1,i,0)).ne.0).and.(ClassCountUnocc2(ClassCountInd(1,i,0)).eq.0)) THEN
!If there are alpha electrons in this class with no possible unoccupied alpha orbitals in the same class, these alpha electrons have no single excitations.
                ElecsWNoExcits=ElecsWNoExcits+ClassCount2(ClassCountInd(1,i,0))
            ENDIF
            IF((ClassCount2(ClassCountInd(2,i,0)).ne.0).and.(ClassCountUnocc2(ClassCountInd(2,i,0)).eq.0)) THEN
                ElecsWNoExcits=ElecsWNoExcits+ClassCount2(ClassCountInd(2,i,0))
            ENDIF
        enddo

        IF(ElecsWNoExcits.eq.NEl) THEN
!There are no single excitations from this determinant at all.
!Then we will create a double excitation instead.
            RETURN
        ENDIF

!We want to pick an occupied electron - Prob = 1/(N-ElecsWNoExcits)
        Attempts=0
        do while(.true.)

            Attempts=Attempts+1

!Choose an electron randomly...
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            Eleci=INT(NEl*r)+1

!Find symmetry of chosen electron
            ElecSym=INT((G1(nI(Eleci))%Sym%S),4)

            IF(G1(nI(Eleci))%Ms.eq.1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                NExcit=ClassCountUnocc2(ClassCountInd(1,ElecSym,0))
                ispn=1  !Alpha ispn=1, Beta ispn=2
            ELSE
!Beta orbital
                NExcit=ClassCountUnocc2(ClassCountInd(2,ElecSym,0))
                ispn=2
            ENDIF

            IF(NExcit.ne.0) EXIT    !Have found electron with allowed excitations

            IF(Attempts.gt.250) THEN
                WRITE(6,*) "Cannot find single excitation from electrons after 250 attempts..."
                CALL WRITEDET(6,nI,NEL,.true.)
!                WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
!                WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
                WRITE(6,*) "***"
!                WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
!                WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
                CALL Stop_All("CreateSingleExcit","Cannot find single excitation from electrons after 250 attempts...")
            ENDIF

        enddo
        ExcitMat(1,1)=nI(Eleci)

!Now we want to run through all sym+spin allowed excitations of the chosen electron, and determine their matrix elements
!To run just through the states of the required symmetry we want to use SymLabelCounts.

!We also want to take into account spin. We want the spin of the chosen unoccupied orbital to be the same as the chosen occupied orbital.
!Run over all possible a orbitals
        Ind=ClassCountInd(ispn,ElecSym,0)
        EndSymState=SymLabelCounts2(1,Ind)+SymLabelCounts2(2,Ind)-1

        VecInd=1
        NormProb=0.D0

       do j=SymLabelCounts2(1,Ind),EndSymState

!            IF(ispn.eq.-1) THEN
!We want to look through all beta orbitals
            OrbA=SymLabelList2(j)     !This is the spin orbital chosen for a
!                OrbA=(2*SymLabelList(j))-1     !This is the spin orbital chosen for a
!            ELSE
!                OrbA=SymLabelList(1,j)
!            ENDIF

            IF(BTEST(ILUT((OrbA-1)/32),MOD((OrbA-1),32))) THEN
!Orbital is in nI...not an unoccupied orbital
                CYCLE
            ENDIF

!Now we want to find the information about this excitation
            ExcitMat(2,1)=OrbA
            CALL Scr1Excit2(NEl,nBasisMax,nI,nI,G1,nBasis,UMat,ALat,iSpinSkip,FCoul,.true.,rh,ExcitMat,.false.)
        
            SpawnProb(VecInd)=abs(REAL(rh%v,r2))
            SpawnOrb(VecInd)=OrbA
            NormProb=NormProb+SpawnProb(VecInd)
            VecInd=VecInd+1
            IF(VecInd.ge.nBasis) THEN
                CALL Stop_All("CreateSingleExcitBiased","Finding too many virtual pairs...")
            ENDIF

        enddo
        IF((VecInd-1).ne.NExcit) THEN
            CALL Stop_All("CreateSingleExcitBiased","Wrong number of excitations found from chosen electron")
        ENDIF
        IF(VecInd.eq.1) THEN
            CALL Stop_All("CreateSingleExcitBiased","No excitations found from electron, but some should exist")
        ENDIF

!We now have to find out how many children to spawn, based on the value of normprob.
        rat=Tau*NormProb*REAL((NEl-ElecsWNoExcits)*nParts,r2)/(1.D0-PDoubNew)
        iCreate=INT(rat)
        rat=rat-REAL(iCreate)
        IF(tMerTwist) THEN
            CALL genrand_real2(r)
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(rat.gt.r) THEN
!Child is created
            iCreate=iCreate+1
        ENDIF

        IF(iCreate.gt.0) THEN
!We want to spawn particles. This only question now is where. Run through the ab pairs again and choose based on the SpawnProb element.
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            r=r*NormProb

            i=0
            do while(r.gt.0.D0)
                i=i+1
                r=r-SpawnProb(i)
            enddo
            IF(i.gt.VecInd-1) THEN
                CALL Stop_All("CreateSingleExcitBiased","Chosen virtual does not correspond to allowed orbital")
            ENDIF

            OrbA=SpawnOrb(i)

!We now know that we want to create iCreate particles, from orbitals nI(Eleci) -> OrbA.
            nJ(:)=nI(:)
!ExcitMat wants to be the index in nI of the orbital to excite from, but returns the actual orbitals.
            ExcitMat(1,1)=Eleci
            ExcitMat(2,1)=OrbA
            CALL FindExcitDet(ExcitMat,nJ,1,TParity)

!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
!            CALL IsSymAllowedExcit(nI,nJ,1,ExcitMat,SymAllowed)

!Once we have the definitive determinant, we also want to find out what sign the particles we want to create are.
!iCreate is initially positive, so its sign can change depending on the sign of the connection and of the parent particle(s)
            rh=GetHElement4(nI,nJ,1,ExcitMat,tParity)

            IF(WSign.gt.0) THEN
                !Parent particle is positive
                IF(real(rh%v).gt.0.D0) THEN
                    iCreate=-iCreate     !-ve walker created
                ENDIF
            ELSE
                IF(real(rh%v).lt.0.D0) THEN
                    iCreate=-iCreate    !-ve walkers created
                ENDIF
            ENDIF

        ENDIF

    END SUBROUTINE CreateSingleExcitBiased
        

    SUBROUTINE CreateDoubExcitBiased(nI,nJ,iLut,ExcitMat,tParity,nParts,WSign,Tau,iCreate)
        INTEGER :: nI(NEl),nJ(NEl),iLut(0:NIfTot),ExcitMat(2,2),iCreate,iSpn,OrbA,OrbB,SymProduct
        INTEGER :: Elec1Ind,Elec2Ind,nParts,WSign,SumMl
        TYPE(HElement) :: rh
        LOGICAL :: tParity
        REAL*8 :: Tau

!First, we need to pick an unbiased distinct electron pair.
!These have symmetry product SymProduct, and spin pair iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha
!The probability for doing this is 1/ElecPairs.
        CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,SumMl,-1)

!This routine runs through all distinct ab pairs for the chosen ij and stochastically chooses how many particles to create.
!If spawning wants to occur, then it runs through the list again and chooses a pair, which it returns.
        CALL CalcAllab(nI,iLut,Elec1Ind,Elec2Ind,SymProduct,iSpn,OrbA,OrbB,nParts,iCreate,Tau)

!We now know that we want to create iCreate particles, from orbitals nI(Elec1/2Ind) -> OrbA + OrbB.
        IF(iCreate.gt.0) THEN
            CALL FindNewDet(nI,nJ,Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat,tParity)

!Once we have the definitive determinant, we also want to find out what sign the particles we want to create are.
!iCreate is initially positive, so its sign can change depending on the sign of the connection and of the parent particle(s)
            rh=GetHElement4(nI,nJ,2,ExcitMat,tParity)

            IF(WSign.gt.0) THEN
                !Parent particle is positive
                IF(real(rh%v).gt.0.D0) THEN
                    iCreate=-iCreate     !-ve walker created
                ENDIF
            ELSE
                IF(real(rh%v).lt.0.D0) THEN
                    iCreate=-iCreate    !-ve walkers created
                ENDIF
            ENDIF

        ENDIF


    END SUBROUTINE CreateDoubExcitBiased

    SUBROUTINE CalcAllab(nI,ILUT,Elec1Ind,Elec2Ind,SymProduct,iSpn,OrbA,OrbB,nParts,iCreate,Tau)
        use Integrals , only : GetUMatEl
        INTEGER :: nI(NEl),iLut(0:NIfTot),Elec1Ind,Elec2Ind,SymProduct,iSpn,OrbA,OrbB,iCreate
        INTEGER :: SpatOrbi,SpatOrbj,Spini,Spinj,i,aspn,bspn,SymA,SymB,SpatOrba,EndSymState,VecInd
        REAL*8 :: Tau,SpawnProb(MaxABPairs),NormProb,rat,r
        INTEGER :: SpawnOrbs(2,MaxABPairs),j,nParts,SpinIndex,Ind
        TYPE(HElement) :: HEl

!We want the spatial orbital number for the ij pair (Elec1Ind is the index in nI).
!Later, we'll have to use GTID for UHF.
        SpatOrbi=((nI(Elec1Ind)-1)/2)+1
        SpatOrbj=((nI(Elec2Ind)-1)/2)+1
        Spini=G1(nI(Elec1Ind))%Ms
        Spinj=G1(nI(Elec2Ind))%Ms
        VecInd=1
        NormProb=0.D0

        do i=1,nBasis
!Run through all a orbitals
            IF(mod(i,2).eq.0) THEN
!We have an alpha spin...
                aspn=1
                IF(iSpn.eq.1) THEN
!ij is an beta/beta pair, so we only want to pick beta a orbitals.
                    CYCLE
                ENDIF
            ELSE
!We have a beta spin...
                aspn=-1
                IF(iSpn.eq.3) THEN
!ij is a alpha/alpha pair, so we only want to pick alpha a orbitals.
                    CYCLE
                ENDIF
            ENDIF
            
!We also want to check that the a orbital we have picked is not already in the determinant we are exciting from.
            IF(BTEST(ILUT((i-1)/32),MOD((i-1),32))) THEN
!Orbital is in nI...do not make a pair from this.
                CYCLE
            ENDIF

!Now we have to run over all b's which can go with a. We only want unique pairs, so we constrain b to be at least a+1.
!We only want to run over symmetry and spin allowed b's though.
!Find the required symmetry of b.
            SymA=INT(G1(i)%Sym%S,4)
            SymB=IEOR(SymA,SymProduct)
            SpatOrba=((i-1)/2)+1

!We also want to take into account spin.
            IF(ispn.eq.1) THEN
                bspn=-1  !Want beta spin b orbitals
                SpinIndex=2
            ELSEIF(ispn.eq.3) THEN
                bspn=1  !Want alpha spin b orbitals
                SpinIndex=1
            ELSE
!ij pair is an alpha/beta spin pair, therefore b wants to be of opposite spin to a.
                IF(aspn.eq.-1) THEN
!a orbital is a beta orbital, therefore we want b to be an alpha orbital.
                    bspn=1
                    SpinIndex=1
                ELSE
                    bspn=-1
                    SpinIndex=2
                ENDIF
            ENDIF

!To run just through the states of the required symmetry we want to use SymLabelCounts.
!            StartSymState=SymLabelCounts(1,SymB+1)
            Ind=ClassCountInd(SpinIndex,SymB,0)
            EndSymState=SymLabelCounts2(1,Ind)+SymLabelCounts2(2,Ind)-1


!Run over all possible b orbitals
            do j=SymLabelCounts2(1,Ind),EndSymState

!                IF(bspn.eq.-1) THEN
!                    OrbB=(2*SymLabelList(j))-1     !This is the spin orbital chosen for b
!                ELSE
!                    OrbB=(2*SymLabelList(j))
!                ENDIF

                OrbB=SymLabelList2(j)     !This is the spin orbital chosen for b

                IF(OrbB.le.i) THEN
!Since we only want unique ab pairs, ensure that b > a.
                    CYCLE
                ENDIF

                IF(BTEST(ILUT((OrbB-1)/32),MOD((OrbB-1),32))) THEN
!Orbital is in nI...do not make a pair from this.
                    CYCLE
                ENDIF

!We have now found an allowed ab pair to go with the ij pair chosen previously - record its stats.
                IF( Spini.EQ.aspn.and.Spinj.eq.bspn) THEN
                    Hel=GETUMATEL(NBASISMAX,UMAT,ALAT,nBasis,iSpinSkip,G1,SpatOrbi,SpatOrbj,SpatOrba,j)
                ELSE
                    Hel=HElement(0.D0)
                ENDIF
                IF(Spini.EQ.bspn.and.Spinj.EQ.aspn) THEN
                    Hel=Hel-GETUMATEL(NBASISMAX,UMAT,ALAT,nBasis,iSpinSkip,G1,SpatOrbi,SpatOrbj,j,SpatOrba)
                ENDIF

                SpawnProb(VecInd)=abs(REAL(Hel%v,r2))
                SpawnOrbs(1,VecInd)=i
                SpawnOrbs(2,VecInd)=OrbB
                NormProb=NormProb+SpawnProb(VecInd)
                VecInd=VecInd+1
                IF(VecInd.ge.MaxABPairs) THEN
                    CALL Stop_All("CalcAllab","Finding too many ab pairs...")
                ENDIF

            enddo
        enddo

        VecInd=VecInd-1     !This now indicates the total number of ab pairs we have found for the chosen ij.

        IF(VecInd.eq.0) THEN
            CALL Stop_All("CalcAllab","No ab pairs found for the chosen ij")
        ENDIF

!We now have to find out how many children to spawn, based on the value of normprob.
        rat=Tau*NormProb*REAL(ElecPairs*nParts,r2)/PDoubNew
        iCreate=INT(rat)
        rat=rat-REAL(iCreate)
        IF(tMerTwist) THEN
            CALL genrand_real2(r)
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(rat.gt.r) THEN
!Child is created
            iCreate=iCreate+1
        ENDIF

        IF(iCreate.gt.0) THEN
!We want to spawn particles. This only question now is where. Run through the ab pairs again and choose based on the SpawnProb element.
            IF(tMerTwist) THEN
                CALL genrand_real2(r)
            ELSE
                CALL RANLUX(r,1)
            ENDIF
            r=r*NormProb

            i=0
            do while(r.gt.0.D0)
                i=i+1
                r=r-SpawnProb(i)
            enddo
            IF(i.gt.VecInd) THEN
                CALL Stop_All("CalcAllab","Chosen ab pair does not correspond to allowed pair")
            ENDIF

            OrbA=SpawnOrbs(1,i)
            OrbB=SpawnOrbs(2,i)
            
        ENDIF

    END SUBROUTINE CalcAllab


END MODULE GenRandSymExcitNUMod


!Sometimes (especially UHF orbitals), the symmetry routines will not set up the orbitals correctly. Therefore, this routine will set up symlabellist and symlabelcounts
!to be cast in terms of spin orbitals, and the symrandexcit2 routines will use these arrays.
SUBROUTINE SpinOrbSymSetup(tRedoSym)
    use SymExcitDataMod , only : ScratchSize,SymLabelList2,SymLabelCounts2,OrbClassCount
    use GenRandSymExcitNUMod , only : ClassCountInd
    use SymData, only: nSymLabels,TwoCycleSymGens
    use SymData, only: SymLabelList,SymLabelCounts
    use SystemData , only : G1,tFixLz,tNoSymGenRandExcits,nBasis,iMaxLz,tUEG
    IMPLICIT NONE
    INTEGER :: i,j,SymInd
    INTEGER :: Spin
    LOGICAL :: tRedoSym
    INTEGER , ALLOCATABLE :: Temp(:)
    
    IF(tFixLz) THEN
!Calculate the upper array bound for the ClassCount2 arrays. This will be dependant on the number of symmetries needed.
        ScratchSize=2*nSymLabels*(2*iMaxLz+1)
    ELSE
        ScratchSize=2*nSymLabels
    ENDIF

    IF(tNoSymGenRandExcits.or.tUEG.or.(.not.TwoCycleSymGens)) THEN
        ScratchSize=2
    ENDIF



!SymLabelList2 and SymLabelCounts2 are now organised differently, so that it is more efficient, and easier to add new symmetries.
!SymLabelCounts is of size (2,ScratchSize), where 1,x gives the index in SymlabelList2 where the orbitals of symmetry x start.
!SymLabelCounts(2,x) tells us the number of orbitals of spin & symmetry x there are.

!Therefore, if you want to run over all orbitals of a specific symmetry, you want to run over SymLabelList from SymLabelCounts(1,sym) to SymLabelCounts(1,sym)+SymLabelCounts(2,sym)-1

    Allocate(SymLabelList2(nBasis))
    Allocate(SymLabelCounts2(2,ScratchSize))
    SymLabelList2(:)=0          !Indices:   spin-orbital number
    SymLabelCounts2(:,:)=0      !Indices:   index/Number , symmetry(inc. spin)
    Allocate(Temp(ScratchSize))
    
    do j=1,nBasis
        IF(G1(j)%Ms.eq.1) THEN
            Spin=1
        ELSE
            Spin=2
        ENDIF
        IF(tUEG.or.tNoSymGenRandExcits.or.(.not.TwoCycleSymGens)) THEN
            SymInd=ClassCountInd(Spin,0,G1(j)%Ml)
        ELSE
            SymInd=ClassCountInd(Spin,INT(G1(j)%Sym%S,4),G1(j)%Ml)
        ENDIF
        SymLabelCounts2(2,SymInd)=SymLabelCounts2(2,SymInd)+1
    enddo
    SymLabelCounts2(1,1)=1
    do j=2,ScratchSize
        SymLabelCounts2(1,j)=SymLabelCounts2(1,j-1)+SymLabelCounts2(2,j-1)
    enddo
    Temp(:)=SymLabelCounts2(1,:)
    do j=1,nBasis
        IF(G1(j)%Ms.eq.1) THEN
            Spin=1
        ELSE
            Spin=2
        ENDIF
        IF(tUEG.or.tNoSymGenRandExcits.or.(.not.TwoCycleSymGens)) THEN
            SymInd=ClassCountInd(Spin,0,G1(j)%Ml)
        ELSE
            SymInd=ClassCountInd(Spin,INT(G1(j)%Sym%S,4),G1(j)%Ml)
        ENDIF
        SymLabelList2(Temp(SymInd))=j
        Temp(SymInd)=Temp(SymInd)+1
    enddo

    Deallocate(Temp)


!
!    Allocate(SymLabelList2(2,(nBasis/2)))   !This will seperate the alpha and beta states.
!    SymLabelList2(:,:)=0
!    Allocate(SymLabelCounts2(2,2,nSymLabels))     !Indices: (Alpha:Beta, Index:Number , Symmetry)
!    SymLabelCounts2(:,:,:)=0
!    IF(tRedoSym) THEN
!        AlphaCounter=1
!        BetaCounter=1
!        IF(tNoSymGenRandExcits) THEN
!            LoopVar=0
!        ELSE
!            LoopVar=nSymLabels-1
!        ENDIF
!        do i=0,LoopVar
!            tFirstSymAlpha=.true.
!            tFirstSymBeta=.true.
!            CountSymAlpha=0
!            CountSymBeta=0
!            do j=1,nBasis
!                IF(tNoSymGenRandExcits) THEN
!                    Sym=0
!                ELSE
!                    Sym=INT(G1(j)%Sym%S,4)
!                ENDIF
!                IF(G1(j)%Ms.eq.1) THEN
!!Alpha orbital
!                    IF(Sym.eq.i) THEN
!                        IF(tFirstSymAlpha) THEN
!                            SymLabelCounts2(1,1,i+1)=AlphaCounter
!                            tFirstSymAlpha=.false.
!                        ENDIF
!                        SymLabelList2(1,AlphaCounter)=j
!                        AlphaCounter=AlphaCounter+1
!                        CountSymAlpha=CountSymAlpha+1
!                    ENDIF
!                ELSE
!!Beta orbital
!                    IF(Sym.eq.i) THEN
!                        IF(tFirstSymBeta) THEN
!                            SymLabelCounts2(2,1,i+1)=BetaCounter
!                            tFirstSymBeta=.false.
!                        ENDIF
!                        SymLabelList2(2,BetaCounter)=j
!                        BetaCounter=BetaCounter+1
!                        CountSymBeta=CountSymBeta+1
!                    ENDIF
!                ENDIF
!            enddo
!            SymLabelCounts2(1,2,i+1)=CountSymAlpha
!            SymLabelCounts2(2,2,i+1)=CountSymBeta
!        enddo
!    ELSE
!        IF(tNoSymGenRandExcits) THEN
!            do i=1,nBasis,2
!                SymLabelList2(2,(i+1)/2)=i
!                SymLabelList2(1,(i+1)/2)=i+1
!            enddo
!            SymLabelCounts2(1,1,1)=1
!            SymLabelCounts2(2,1,1)=1
!            SymLabelCounts2(1,2,1)=nBasis/2
!            SymLabelCounts2(2,2,1)=nBasis/2
!        ELSE
!            do i=1,nBasis/2
!                SymLabelList2(1,i)=2*SymLabelList(i)
!                SymLabelList2(2,i)=(2*SymLabelList(i))-1
!            enddo
!            do i=1,nSymLabels
!                SymLabelCounts2(1,1,i)=SymLabelCounts(1,i)
!                SymLabelCounts2(2,1,i)=SymLabelCounts(1,i)
!                SymLabelCounts2(1,2,i)=SymLabelCounts(2,i)
!                SymLabelCounts2(2,2,i)=SymLabelCounts(2,i)
!            enddo
!        ENDIF
!    ENDIF

!        WRITE(6,*) 'Symmetries of orbitals 1:nBasis'
!        WRITE(6,*) G1(1:nBasis)%Sym%S

!        do x=1,2
!            IF(x.eq.1) WRITE(6,*) '******* ALPHA ********'
!            IF(x.eq.2) WRITE(6,*) '******* BETA ********'
!            WRITE(6,*) "***"
!            WRITE(6,*) 'SymLabelList2(1:nBasis) - the orbitals in the correct symmetry order, followed by their symmetry'
!            do i=1,(nBasis/2)
!                WRITE(6,*) SymLabelList2(x,i),G1(SymLabelList2(x,i))%Sym%S
!            enddo
!            WRITE(6,*) 'SymLabelCounts2(1,1,S) - the index in SymLabelList2 where symmetry S starts'
!            do i=1,nSymLabels
!                WRITE(6,*) SymLabelCounts2(x,1,i)
!            enddo
!            WRITE(6,*) "***"
!            WRITE(6,*) 'SymLabelCounts2(1,2,S) - the number of orbitals in SymLabelList2 with symmetry S'
!            do i=1,nSymLabels
!                WRITE(6,*) SymLabelCounts2(x,2,i)
!            enddo
!        enddo


    ALLOCATE(OrbClassCount(ScratchSize))
    OrbClassCount(:)=0
    
    
    IF(tNoSymGenRandExcits.or.tUEG.or.(.not.TwoCycleSymGens)) THEN
!All orbitals are in irrep 0
        OrbClassCount(ClassCountInd(1,0,0))=(nBasis/2)
        OrbClassCount(ClassCountInd(2,0,0))=(nBasis/2)
    ELSE
        do i=1,nBasis
            IF(G1(i)%Ms.eq.1) THEN
                OrbClassCount(ClassCountInd(1,INT(G1(i)%Sym%S,4),G1(i)%Ml))=OrbClassCount(ClassCountInd(1,INT(G1(i)%Sym%S,4),G1(i)%Ml))+1
            ELSE
                OrbClassCount(ClassCountInd(2,INT(G1(i)%Sym%S,4),G1(i)%Ml))=OrbClassCount(ClassCountInd(2,INT(G1(i)%Sym%S,4),G1(i)%Ml))+1
            ENDIF
        enddo
    ENDIF


!        ELSE
!!SymLabelCounts(2,1:nSymLabels) gives the number of *states* in each symmetry class.
!!There are therefore equal number of alpha and beta orbitals in each state from which to calculate the unoccupied classcount.
!            do i=1,nSymLabels
!                OrbClassCount(ClassCountInd(1,i-1,0))=SymLabelCounts2(1,2,i)
!                OrbClassCount(ClassCountInd(2,i-1,0))=SymLabelCounts2(2,2,i)
!            enddo
!        ENDIF

END SUBROUTINE SpinOrbSymSetup


!This routine will take a determinant nI, and find Iterations number of excitations of it. It will then histogram these, summing in 1/pGen for every occurance of
!the excitation. This means that all excitations should be 0 or 1 after enough iterations. It will then count the excitations and compare the number to the
!number of excitations generated using the full enumeration excitation generation. This can be done for both doubles and singles, or one of them.
SUBROUTINE TestGenRandSymExcitNU(nI,Iterations,pDoub,exFlag,iWriteEvery)
    Use SystemData , only : NEl,nBasis,G1,nBasisMax,LzTot,NIfTot
    Use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,ConstructClassCounts,ScratchSize
    Use SymData , only : nSymLabels
    use Parallel
!    use soft_exit , only : ChangeVars 
    use DetBitOps , only : EncodeBitDet
    IMPLICIT NONE
    INTEGER :: i,Iterations,exFlag,nI(NEl),nJ(NEl),IC,ExcitMat(2,2),DetConn
    REAL*8 :: pDoub,pGen,AverageContrib,AllAverageContrib
    INTEGER :: ClassCount2(ScratchSize),iLut(0:NIfTot),Scratch1(ScratchSize),Scratch2(ScratchSize),iLutnJ(0:NIfTot)
    INTEGER :: ClassCountUnocc2(ScratchSize),iExcit,iWriteEvery
    LOGICAL :: tParity,SymAllowed,tFilled,tSoftExitFound,tDummy,tDummy2
    REAL*8 , ALLOCATABLE :: DoublesHist(:,:,:,:),SinglesHist(:,:),AllDoublesHist(:,:,:,:),AllSinglesHist(:,:)
    INTEGER , ALLOCATABLE :: EXCITGEN(:)
    INTEGER :: ierr,Ind1,Ind2,Ind3,Ind4,iMaxExcit,nStore(6),nExcitMemLen,j,k,l,DetNum,DetNumS,Lz,excitcount,ForbiddenIter,error

    WRITE(6,*) nI(:)
    WRITE(6,*) Iterations,pDoub,exFlag
    WRITE(6,*) "nSymLabels: ",nSymLabels
    CALL FLUSH(6)

!Find the number of symmetry allowed excitations there should be by looking at the full excitation generator.
!Setup excit generators for this determinant
    iMaxExcit=0
    nStore(1:6)=0
    CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
    ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
    IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
    EXCITGEN(:)=0
    CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGEN,nJ,iMaxExcit,0,nStore,exFlag)
!    CALL GetSymExcitCount(EXCITGEN,DetConn)
    excitcount=0

lp2: do while(.true.)
        CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.false.,EXCITGEN,nJ,iExcit,0,nStore,exFlag)
        IF(nJ(1).eq.0) exit lp2
        CALL GetLz(nJ,NEl,Lz)
        IF(Lz.eq.LzTot) THEN
            excitcount=excitcount+1
            CALL EncodeBitDet(nJ,iLutnJ)
            WRITE(25,*) excitcount,iExcit,iLutnJ(0)
        ENDIF
    enddo lp2

    WRITE(6,*) "Determinant has ",excitcount," total excitations from it."
    CALL FLUSH(6)

!Allocate memory for histogramming determinants
    ALLOCATE(DoublesHist(nBasis,nBasis,nBasis,nBasis),stat=ierr)
    ALLOCATE(AllDoublesHist(nBasis,nBasis,nBasis,nBasis),stat=ierr)
    IF(ierr.ne.0) THEN
        CALL Stop_All("TestGenRandSymExcitNU","Not possible to allocate memory to do histogramming")
    ENDIF
    ALLOCATE(SinglesHist(nBasis,nBasis),stat=ierr)
    ALLOCATE(AllSinglesHist(nBasis,nBasis),stat=ierr)
    IF(ierr.ne.0) THEN
        CALL Stop_All("TestGenRandSymExcitNU","Not possible to allocate memory to do histogramming")
    ENDIF
    DoublesHist(:,:,:,:)=0.D0
    SinglesHist(:,:)=0.D0
    AllDoublesHist(:,:,:,:)=0.D0
    AllSinglesHist(:,:)=0.D0

!    do i=1,NEl
!!Create ILUT for O[1] comparison of orbitals in root determinant - This is now read in
!        ILUT((nI(i)-1)/32)=IBSET(ILUT((NI(i)-1)/32),MOD(NI(i)-1,32))
!    enddo
    CALL EncodeBitDet(nI,iLut)

    tFilled=.false.
    Scratch1(:)=0
    Scratch2(:)=0

    AverageContrib=0.D0
    AllAverageContrib=0.D0
    ForbiddenIter=0
!    pDoub=1.D0
    OPEN(9,FILE="AvContrib",STATUS="UNKNOWN")

    do i=1,Iterations

        IF((mod(i,10000).eq.0).and.(iProcIndex.eq.0)) THEN
            WRITE(6,"(A,I10)") "Iteration: ",i
        ENDIF

        CALL GenRandSymExcitScratchNU(nI,iLut,nJ,pDoub,IC,ExcitMat,TParity,exFlag,pGen,Scratch1,Scratch2,tFilled)
        IF(nJ(1).eq.0) THEN
!            ForbiddenIter=ForbiddenIter+1
            CYCLE
        ENDIF
        AverageContrib=AverageContrib+1.D0/pGen

!        CALL EncodeBitDet(nJ,iLutnJ)
!        IF(IC.eq.1) THEN
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(2,1)
!        ELSE
!            WRITE(6,*) "Double Created"
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(1,2),ExcitMat(2,1),ExcitMat(2,2)
!        ENDIF

        IF(IC.eq.1) THEN
            SinglesHist(ExcitMat(1,1),ExcitMat(2,1))=SinglesHist(ExcitMat(1,1),ExcitMat(2,1))+(1.D0/pGen)
!            SinglesNum(ExcitMat(1,1),ExcitMat(2,1))=SinglesNum(ExcitMat(1,1),ExcitMat(2,1))+1
        ELSE
!Have to make sure that orbitals are in the same order...
            IF(ExcitMat(1,1).gt.ExcitMat(1,2)) THEN
                Ind1=ExcitMat(1,2)
                Ind2=ExcitMat(1,1)
            ELSE
                Ind1=ExcitMat(1,1)
                Ind2=ExcitMat(1,2)
            ENDIF
            IF(ExcitMat(2,1).gt.ExcitMat(2,2)) THEN
                Ind3=ExcitMat(2,2)
                Ind4=ExcitMat(2,1)
            ELSE
                Ind3=ExcitMat(2,1)
                Ind4=ExcitMat(2,2)
            ENDIF
            DoublesHist(Ind1,Ind2,Ind3,Ind4)=DoublesHist(Ind1,Ind2,Ind3,Ind4)+(1.D0/pGen)
        ENDIF
        IF(mod(i,iWriteEvery).eq.0) THEN
            AllAverageContrib=0.D0
#ifdef PARALLEL
            CALL MPI_Reduce(AverageContrib,AllAverageContrib,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
#else            
            AllAverageContrib=AverageContrib
#endif
            IF(iProcIndex.eq.0) THEN
                WRITE(9,*) i,AllAverageContrib/(REAL(i,8)*excitcount*nProcessors)
            ENDIF
!            CALL ChangeVars(tDummy,tSoftExitFound,tDummy2)
!            IF(tSoftExitFound) EXIT
        ENDIF

!Check excitation
!        CALL IsSymAllowedExcit(nI,nJ,IC,ExcitMat,SymAllowed)

    enddo

    CLOSE(9)

#ifdef PARALLEL
    CALL MPI_Reduce(DoublesHist,AllDoublesHist,nBasis**4,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    CALL MPI_Reduce(SinglesHist,AllSinglesHist,nBasis**2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
#else
    AllDoublesHist=DoublesHist
    AllSinglesHist=SinglesHist
#endif

!Now run through arrays normalising them so that numbers are more managable.
    OPEN(8,FILE="DoublesHist",STATUS="UNKNOWN")
    DetNum=0
    do i=1,nBasis-1
        do j=i+1,nBasis
            do k=1,nBasis-1
                do l=k+1,nBasis
                    IF(AllDoublesHist(i,j,k,l).gt.0.D0) THEN
!                        DoublesHist(i,j,k,l)=DoublesHist(i,j,k,l)/real(Iterations,8)
                        DetNum=DetNum+1
                        ExcitMat(1,1)=i
                        ExcitMat(1,2)=j
                        ExcitMat(2,1)=k
                        ExcitMat(2,2)=l
                        CALL FindExcitBitDet(iLut,iLutnJ,2,ExcitMat,NIfTot)
                        WRITE(8,"(I12,F20.12,4I5,I15)") DetNum,AllDoublesHist(i,j,k,l)/(real(Iterations,8)*nProcessors),i,j,k,l,iLutnJ(0)
                    ENDIF
                enddo
            enddo
        enddo
    enddo
    CLOSE(8)
    WRITE(6,*) DetNum," Double excitations found from nI"
    OPEN(9,FILE="SinglesHist",STATUS="UNKNOWN")
    DetNumS=0
    do i=1,nBasis
        do j=1,nBasis
            IF(AllSinglesHist(i,j).gt.0.D0) THEN
                DetNumS=DetNumS+1
                ExcitMat(1,1)=i
                ExcitMat(2,1)=j
                CALL FindExcitBitDet(iLut,iLutnJ,1,ExcitMat,NIfTot)
                WRITE(9,*) DetNumS,AllSinglesHist(i,j)/(real(Iterations,8)*nProcessors),iLutnJ(0)
            ENDIF
        enddo
    enddo
    CLOSE(9)
    WRITE(6,*) DetNumS," Single excitations found from nI"
    IF((DetNum+DetNumS).ne.ExcitCount) THEN
        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)
        WRITE(6,*) "Total determinants = ", ExcitCount
!        WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
!        WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
        WRITE(6,*) "***"
!        WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
!        WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
        CALL Stop_All("TestGenRandSymExcitNU","Not all excitations accounted for...")
    ENDIF

END SUBROUTINE TestGenRandSymExcitNU

SUBROUTINE IsSymAllowedExcit(nI,nJ,IC,ExcitMat,SymAllowed)
    Use SystemData , only : G1,NEl,tFixLz
    Use SystemData , only : Symmetry,tNoSymGenRandExcits
    IMPLICIT NONE
    Type(Symmetry) :: SymProduct,SymProduct2,SYMPROD
    LOGICAL :: SYMEQ,ISVALIDDET,SymAllowed
    INTEGER :: IC,ExcitMat(2,2),nI(NEl),nJ(NEl),ExcitLevel,iGetExcitLevel
    INTEGER :: KOcc,KUnocc

     SymAllowed=.true.
     Excitlevel=iGetExcitLevel(nI,nJ,NEl)
     IF(Excitlevel.ne.IC) THEN
         WRITE(6,*) "Have not created a correct excitation"
        CALL WRITEDET(6,nI,NEL,.TRUE.)
        CALL WRITEDET(6,nJ,NEL,.TRUE.)
        STOP "Have not created a correct excitation"
     ENDIF
     IF(.NOT.ISVALIDDET(nJ,NEL)) THEN
         WRITE(6,*) "INVALID DET"
         CALL WRITEDET(6,nI,NEL,.TRUE.)
         CALL WRITEDET(6,nJ,NEL,.TRUE.)
         STOP "INVALID DET"
     ENDIF
     
     IF(.not.tNoSymGenRandExcits) THEN
         IF(IC.eq.2) THEN
            SymProduct=SYMPROD(G1(ExcitMat(1,1))%Sym,G1(ExcitMat(1,2))%Sym)
            SymProduct2=SYMPROD(G1(ExcitMat(2,1))%Sym,G1(ExcitMat(2,2))%Sym)
            IF(.not.SYMEQ(SymProduct,SymProduct2)) THEN
                SymAllowed=.false.
                CALL Stop_All("IsSymAllowedExcit","Excitation not a valid symmetry allowed double excitation")
            ENDIF
        ELSE
            IF(.not.SYMEQ(G1(ExcitMat(1,1))%Sym,G1(ExcitMat(2,1))%Sym)) THEN
                SymAllowed=.false.
                CALL Stop_All("IsSymAllowedExcit","Excitation not a valid symmetry allowed single excitation")
            ENDIF
            IF(G1(ExcitMat(1,1))%Ms.ne.G1(ExcitMat(2,1))%Ms) THEN
                SymAllowed=.false.
                CALL Stop_All("IsSymAllowedExcit","Excitation not a valid spin-symmetry allowed single excitation")
            ENDIF
        ENDIF
    ENDIF
    IF(tFixLz) THEN
        IF(IC.eq.1) THEN
            IF(G1(ExcitMat(1,1))%Ml.ne.G1(ExcitMat(2,1))%Ml) THEN
                CALL Stop_All("IsSymAllowedExcit","Excitation not a valid momentum allowed single excitation")
            ENDIF
        ELSE
            kOcc=G1(ExcitMat(1,1))%Ml+G1(ExcitMat(1,2))%Ml
            KUnocc=G1(ExcitMat(2,1))%Ml+G1(ExcitMat(2,2))%Ml
            IF(Kocc.ne.KUnocc) THEN
                CALL Stop_All("IsSymAllowedExcit","Excitation not a valid momentum allowed double excitation")
            ENDIF
        ENDIF
    ENDIF
        

END SUBROUTINE IsSymAllowedExcit
