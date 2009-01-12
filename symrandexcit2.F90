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

    use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
    use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
    use SystemData, only: nEl,G1, nBasis,nBasisMax
    use SystemData, only: Arr,nMax,tAssumeSizeExcitgen,tCycleOrbs
    use IntegralsData, only: UMat
    use SymData, only: SymClassSize,nSymPairProds,nSymLabels,TwoCycleSymGens
    use SymData, only: SymLabelList,SymLabelCounts
    IMPLICIT NONE
    SAVE

    INTEGER :: NOccAlpha,NOccBeta,ElecPairs
    REAL*8 :: pDoubNew
    INTEGER , PARAMETER :: r2=kind(0.d0)

    contains

    SUBROUTINE GenRandSymExcitNU(nI,nJ,iSeed,pDoub,IC,ExcitMat,TParity,exFlag,pGen)
        INTEGER :: nI(NEl),nJ(NEl),IC,ExcitMat(2,2),Attempts,exFlag,iSeed
        INTEGER :: ClassCount2(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        INTEGER :: ILUT(0:nBasis/32)
        LOGICAL :: tNoSuccess,tParity
        REAL*8 :: pDoub,pGen,Ran2
        CHARACTER , PARAMETER :: this_routine='GenRandSymExcitNU'

        IF(.not.TwoCycleSymGens) THEN
!Currently only available for molecular systems
            CALL Stop_All(this_routine,"GenRandSymExcitNU can only be used for molecular systems")
        ENDIF

!First, we need to do an O[N] operation to find the number of occupied alpha electrons, number of occupied beta electrons
!and number of occupied electrons of each symmetry class and spin. This is similar to the ClassCount array.
!This has the format (Spn,sym), where Spin=1,2 corresponding to alpha and beta.
!For molecular systems, sym runs from 0 to 7. This is NOT general and should be made so using SymLabels.
!This could be stored to save doing this multiple times, but shouldn't be too costly an operation.
        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2,ILUT)

!ExFlag is 1 for singles, 2 for just doubles, and 3 for both.
        IF(ExFlag.eq.3) THEN
!Choose whether to generate a double or single excitation. Prob of generating a double is given by pDoub.
            pDoubNew=pDoub
            IF(pDoubNew.gt.1.D0) CALL Stop_All(this_routine,"pDoub is greater than 1")

            IF(Ran2(iSeed).lt.pDoubNew) THEN
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
            CALL CreateDoubExcit(nI,nJ,iSeed,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
        ELSE
            CALL CreateSingleExcit(nI,nJ,iSeed,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
            IF(pGen.eq.-1.D0) THEN
                IF(ExFlag.ne.3) THEN
                    CALL Stop_All("GenRandSymExcitNU","Found determinant with no singles, but can only have got here from single. Should never be in this position!")
                ENDIF
                pDoubNew=1.D0
                IC=2
                CALL CreateDoubExcit(nI,nJ,iSeed,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
            ENDIF

        ENDIF

    END SUBROUTINE GenRandSymExcitNU

    SUBROUTINE CreateDoubExcit(nI,nJ,iSeed,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
        INTEGER :: nI(NEl),nJ(NEl),iSeed,ExcitMat(2,2),NExcitOtherWay,OrbB
        INTEGER :: ClassCount2(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        INTEGER :: ILUT(0:nBasis/32),NExcitB,SpinOrbA,OrbA,SymB,NExcitA
        INTEGER :: Elec1Ind,Elec2Ind,SymProduct,iSpn,ForbiddenOrbs,SymA
        REAL*8 :: pGen
        LOGICAL :: tParity

!First, we need to pick an unbiased distinct electron pair.
!These have symmetry product SymProduct, and spin pair iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha
        CALL PickElecPair(nI,Elec1Ind,Elec2Ind,iSeed,SymProduct,iSpn)

!This routine finds the number of orbitals which are allowed by spin, but not part of any spatial symmetry allowed unoccupied pairs.
!This number is needed for the correct normalisation of the probability of drawing any given A orbital since these can be chucked and redrawn.
        CALL FindNumForbiddenOrbs(ForbiddenOrbs,ClassCountUnocc2,SymProduct,iSpn)

!Now we have to pick the first unoccupied orbital. If an orbital is not present in any allowed pairs, it is chucked and a new one drawn.
!The number NExcit is the number of unoccupied orbitals that the orbital was chosen from (including symmetry-forbidden orbital pairs)
!Arguments:     NExcit = Number of possible spin-allowed unoccupied spinorbitals, including forbidden orbs (these will be chucked)
!               SpinOrbA = Spin of the chosen spin-orbital. -1 is alpha, 1 is beta.
!               OrbA = Index of the chosen spinorbital
!               SymB = Symmetry required of the second unoccupied spinorbital, so that sym(A) x sym(B) = SymProduct
!               SymProduct = (intent in), the symmetry of electron pair chosen
        CALL PickAOrb(nI,iSeed,iSpn,ILUT,ClassCountUnocc2,NExcitA,Elec1Ind,Elec2Ind,SpinOrbA,OrbA,SymA,SymB,SymProduct)

!This routine will pick an unoccupied orbital at random from a specified spin and symmetry class.
!There should definitely be a possible spinorbital, since A was chosen so that there was one.
!We have to make sure with alpha/alpha or beta/beta pairs and when SymProduct=0, that we don't choose the same unoccupied orbital.
!If we do this, then we should chuck and redraw, since there should definitely be another allowed spinorbital in the class.
!We return the number of allowed B's for the A we picked in NExcitB, however we also need to know the number of allowed A's if we
!had picked B first. This will be returned in NExcitOtherWay.
        CALL PickBOrb(nI,iSeed,iSpn,ILUT,ClassCountUnocc2,SpinOrbA,OrbA,SymA,OrbB,SymB,NExcitB,SymProduct,NExcitOtherWay)

        CALL FindNewDet(nI,nJ,Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat,tParity)

        CALL FindDoubleProb(ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay,pGen)

    END SUBROUTINE CreateDoubExcit

!This routine creates the final determinant.
    SUBROUTINE FindNewDet(nI,nJ,Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat,tParity)
        INTEGER :: nI(NEl),nJ(NEl),Elec1Ind,Elec2Ind,OrbA,OrbB,ExcitMat(2,2)
        INTEGER :: iGetExcitLevel_2,ExcitLevel
        LOGICAL :: tParity,IsValidDet

!First construct ExcitMat
        ExcitMat(1,1)=Elec1Ind
        ExcitMat(2,1)=OrbA
        ExcitMat(1,2)=Elec2Ind
        ExcitMat(2,2)=OrbB
        nJ(:)=nI(:)
        CALL FindExcitDet(ExcitMat,nJ,2,tParity)

!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
        CALL IsSymAllowedExcit(nI,nJ,2,ExcitMat)

    END SUBROUTINE FindNewDet

!This routine finds the probability of creating the excitation. See the header of the file for more information on how this works.
    SUBROUTINE FindDoubleProb(ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay,pGen)
        INTEGER :: ForbiddenOrbs,NExcitA,NExcitB,NExcitOtherWay
        REAL*8 :: pGen,PabGivenij

        PabGivenij=(1.D0/real((NExcitA-ForbiddenOrbs),r2))*((1.D0/real(NExcitB,r2))+(1.D0/real(NExcitOtherWay,r2)))
        pGen=pDoubNew*(1.D0/real(ElecPairs,r2))*PabGivenij

    END SUBROUTINE FindDoubleProb

    SUBROUTINE PickBOrb(nI,iSeed,iSpn,ILUT,ClassCountUnocc2,SpinOrbA,OrbA,SymA,OrbB,SymB,NExcit,SymProduct,NExcitOtherWay)
        INTEGER :: nI(NEl),iSeed,iSpn,SpinOrbA,OrbA,SymB,NExcit,SymProduct,NExcitOtherWay
        INTEGER :: OrbB,Attempts,SpinOrbB,ChosenUnocc
        INTEGER :: ILUT(0:nBasis/32),SymA,nOrbs,z,i
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        REAL*8 :: Ran2

!We want to calculate the number of possible B's given the symmetry and spin it has to be since we have already picked A.
!We have calculated in NExcit the number of orbitals available for B given A, but we also need to know the number of orbitals to choose from for A IF
!we had picked B first.
        IF(iSpn.eq.2) THEN
!If iSpn=2, then we want to find a spinorbital of the opposite spin of SpinOrbA
            IF(SpinOrbA.eq.-1) THEN
!We have already picked an alpha orbital, so now we want to pick a beta orbital. Find out how many of these there are.
                NExcit=ClassCountUnocc2(2,SymB)
                NExcitOtherWay=ClassCountUnocc2(1,SymA)
                SpinOrbB=0  !This is defined differently to SpinOrbA. 0=Beta, -1=alpha.
            ELSE
!Want to pick an alpha orbital.
                NExcit=ClassCountUnocc2(1,SymB)
                NExcitOtherWay=ClassCountUnocc2(2,SymA)
                SpinOrbB=-1
            ENDIF
        ELSEIF(iSpn.eq.1) THEN
!Definitely want a beta orbital
            NExcit=ClassCountUnocc2(2,SymB)
            NExcitOtherWay=ClassCountUnocc2(2,SymA)
            SpinOrbB=0
        ELSE
!Definitely want an alpha orbital
            NExcit=ClassCountUnocc2(1,SymB)
            NExcitOtherWay=ClassCountUnocc2(1,SymA)
            SpinOrbB=-1
        ENDIF

        IF((iSpn.ne.2).and.(SymProduct.eq.0)) THEN
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
            ChosenUnocc=INT(NExcit*Ran2(iSeed))+1

!Now run through all allowed orbitals until we find this one.
            nOrbs=SymLabelCounts(2,SymB+1)
            z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
            do i=0,nOrbs-1
!Find the spin orbital index. SymLabelCounts has the index of the state for the given symmetry.
                OrbB=(2*SymLabelList(SymLabelCounts(1,SymB+1)+i))+SpinOrbB

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
                CALL Stop_All("PickBOrb","Could not find allowed unoccupied orbital to excite to.")
            ENDIF

        ELSE
! METHOD 2 (Keep drawing orbitals from the desired symmetry and spin until we find one unoccupied)
! =========================

            nOrbs=SymLabelCounts(2,SymB+1)
            Attempts=0
            do while(.true.)
                Attempts=Attempts+1
                
!Draw randomly from the set of orbitals
                ChosenUnocc=INT(nOrbs*Ran2(iSeed))
                OrbB=(2*SymLabelList(SymLabelCounts(1,SymB+1)+ChosenUnocc))+SpinOrbB

!Find out if orbital is in nI or not. Accept if it isn't in it.
                IF((.not.(BTEST(ILUT((OrbB-1)/32),MOD((OrbB-1),32)))).and.(OrbB.ne.OrbA)) THEN
!Orbital not in nI. Accept.
                    EXIT
                ENDIF
                
                IF(Attempts.gt.250) THEN
                    WRITE(6,*) "Cannot find double excitation unoccupied orbital after 250 attempts..."
                    CALL WRITEDET(6,nI,NEL,.true.)
                    CALL Stop_All("PickBOrb","Cannot find double excitation unoccupied orbital after 250 attempts...")
                ENDIF

            enddo

        ENDIF

    END SUBROUTINE PickBOrb

!This routine finds the number of orbitals which are allowed by spin, but not part of any spatial symmetry allowed unoccupied pairs.
!This number is needed for the correct normalisation of the probability of drawing any given A orbital since these can be chucked and redrawn.
    SUBROUTINE FindNumForbiddenOrbs(ForbiddenOrbs,ClassCountUnocc2,SymProduct,iSpn)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        INTEGER :: ForbiddenOrbs,SymProduct,iSpn,i,ConjSym

        ForbiddenOrbs=0
        IF(iSpn.eq.2) THEN
!i,j are an alpha/beta pair. The number of forbidden orbitals includes all alphas and betas.

            do i=0,nSymLabels-1
!Run though all symmetries
                IF(ClassCountUnocc2(1,i).eq.0) THEN
!This symmetry has no unoccupied alpha orbitals - does its symmetry conjugate have any unoccupied beta orbitals which are now forbidden?
!If there are no unoccupied orbitals in this conjugate symmetry, then it won't increase the forbidden orbital number, since it can never be chosen.
                    ConjSym=IEOR(SymProduct,i)
                    ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(2,ConjSym) !No unocc alphas in i, therefore all betas in ConjSym are forbidden
                ENDIF
                IF(ClassCountUnocc2(2,i).eq.0) THEN
!This symmetry has no unoccupied beta orbitals - does its symmetry conjugate have any unoccupied alpha orbitals which are now forbidden?
!If there are no unoccupied orbitals in this conjugate symmetry, then it won't increase the forbidden orbital number, since it can never be chosen.
                    ConjSym=IEOR(SymProduct,i)
                    ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(1,ConjSym)
                ENDIF
            enddo

        ELSEIF(iSpn.eq.1) THEN
            IF(SymProduct.ne.0) THEN
!i,j are a beta/beta pair. The number of forbidden orbitals is just betas
                do i=0,nSymLabels-1
                    IF(ClassCountUnocc2(2,i).eq.0) THEN
                        ConjSym=IEOR(SymProduct,i)
                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(2,ConjSym)
                    ENDIF
                enddo
            ELSE
!There is a subtle point here, which could change the probabilities.
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
                do i=0,nSymLabels-1
                    IF(ClassCountUnocc2(2,i).eq.1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                        ForbiddenOrbs=ForbiddenOrbs+1
                    ENDIF
                enddo
            ENDIF
        ELSEIF(iSpn.eq.3) THEN
            IF(SymProduct.ne.0) THEN
!i,j are a alpha/alpha pair. The number of forbidden orbitals is just alphas
                do i=0,nSymLabels-1
                    IF(ClassCountUnocc2(1,i).eq.0) THEN
                        ConjSym=IEOR(SymProduct,i)
                        ForbiddenOrbs=ForbiddenOrbs+ClassCountUnocc2(1,ConjSym)
                    ENDIF
                enddo
            ELSE
!There is a subtle point here, which could change the probabilities.
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
                do i=0,nSymLabels-1
                    IF(ClassCountUnocc2(1,i).eq.1) THEN
!The one alpha orbital in this class is forbidden, since it cannot form a pair.
                        ForbiddenOrbs=ForbiddenOrbs+1
                    ENDIF
                enddo
            ENDIF
        ENDIF

    END SUBROUTINE FindNumForbiddenOrbs

!Choose the first unoccupied orbital to excite to. Spatial symmetry does not need to be considered, but spin symmetry does.
!After one has been chosen, we have to make sure that an allowed symmetry pair can be formed, otherwise the orbital is
!chucked and a new one drawn.
!Arguments:     NExcit = Number of possible spin-allowed unoccupied spinorbitals, including forbidden orbs (these will be chucked)
!               SpinOrbA = Spin of the chosen spin-orbital. -1 is alpha, 1 is beta.
!               OrbA = Index of the chosen spinorbital
!               SymB = Symmetry required of the second unoccupied spinorbital, so that sym(A) x sym(B) = SymProduct
!               SymProduct = (intent in), the symmetry of electron pair chosen
    SUBROUTINE PickAOrb(nI,iSeed,iSpn,ILUT,ClassCountUnocc2,NExcit,Elec1Ind,Elec2Ind,SpinOrbA,OrbA,SymA,SymB,SymProduct)
        INTEGER :: nI(NEl),iSeed,iSpn,Elec1Ind,Elec2Ind,SpinOrbA,AttemptsOverall,SymA
        INTEGER :: NExcit,ChosenUnocc,z,i,OrbA,Attempts,SymB,SymProduct
        INTEGER :: ILUT(0:nBasis/32)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        REAL*8 :: Ran2

        IF(iSpn.eq.2) THEN
!There is no restriction on whether we choose an alpha or beta spin, so there are nBasis-NEl possible spinorbitals to choose from.
!Therefore the spinOrbA variable has to be set after we have chosen one. For the other iSpn types, we can set it now.
            NExcit=nBasis-NEl
        ELSEIF(iSpn.eq.1) THEN 
!SpinOrbA indicates the spin of the chosen A orb. This is only really useful for iSpn=2
            SpinOrbA=1  !Going to pick a beta orb
            NExcit=(nBasis/2)-nOccBeta  !This is the number of unoccupied beta spinorbitals there are.
        ELSE    !iSpn is 3 - alpha/alpha pair.
            SpinOrbA=-1 !Going to pick an alpha orb
            NExcit=(nBasis/2)-nOccAlpha !This is the number of unoccupied alpha spinorbitals there are.
        ENDIF

        AttemptsOverall=0
        do while(.true.)
!Keep drawing unoccupied orbitals, until we find one which has an allowed partner to form a symmetry-allowed unoccupied pair.
            AttemptsOverall=AttemptsOverall+1

            IF(iSpn.eq.2) THEN
!Electrons chosen were an alpha/beta pair, therefore first randomly chosen orbital can be an alpha OR beta orbital - no restriction.
            
                IF(tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals to find desired one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
                    ChosenUnocc=INT(NExcit*Ran2(iSeed))+1

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
                        Attempts=Attempts+1
                        
!Draw randomly from the set of orbitals
                        ChosenUnocc=INT(nBasis*Ran2(iSeed))

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
                    ChosenUnocc=INT(NExcit*Ran2(iSeed))+1

!Now run through all allowed orbitals until we find this one.
!Alpha orbitals are the ODD numbered orbitals. Beta orbitals are EVEN numbered orbitals.

                    z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
!This could be changed to an O[N] operation, rather than O[M] with a little thought...
                    do i=1,nBasis/2
                        IF(iSpn.eq.1) THEN
!We want to run through all beta orbitals (even numbered basis function)
                            OrbA=(2*i)
                        ELSE
!We want to run through all alpha orbitals (odd numbered basis functions)
                            OrbA=(2*i)-1
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
! METHOD 2 (Keep drawing orbitals randomly until we find one unoccupied). This should be more efficient, unless we have v. small basis sets.
! =========================

                    Attempts=0
                    do while(.true.)
                        Attempts=Attempts+1
                        
!Draw randomly from the set of orbitals
                        ChosenUnocc=INT((nBasis/2)*Ran2(iSeed))+1
                        IF(iSpn.eq.1) THEN
!We only want to choose beta orbitals(i.e. even).
                            ChosenUnocc=2*ChosenUnocc
                        ELSE
!We only want to choose alpha orbitals(i.e. odd).
                            ChosenUnocc=(2*ChosenUnocc)-1
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

                    enddo

                    OrbA=ChosenUnocc  !This is the allowed orbital

                ENDIF

            ENDIF

!We now need to test whether this orbital has any symmetry-allowed unoccupied orbitals to form a pair with.
!To test this, we need to find the needed symmetry of B, in order that Sym(A) x Sym(B) = SymProduct
            SymA=INT(G1(OrbA)%Sym%S,4)
            SymB=IEOR(SymA,SymProduct)
            
            IF(iSpn.eq.2) THEN
!We want an alpha/beta unocc pair. 
                IF(SpinOrbA.eq.-1) THEN
!We have picked an alpha orbital - check to see if there are allowed beta unoccupied orbitals from the SymB Class.
                    IF(ClassCountUnocc2(2,SymB).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                        EXIT
                    ENDIF
                ELSE
!We have picked a beta orbital - check to see if there are allowed alpha unoccupied orbitals from the SymB Class.
                    IF(ClassCountUnocc2(1,SymB).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                        EXIT
                    ENDIF
                ENDIF
            ELSEIF(iSpn.eq.1) THEN
!We want a beta/beta pair.
                IF(SymProduct.ne.0) THEN
!Check to see if there are any unoccupied beta orbitals in the SymB Class.
                    IF(ClassCountUnocc2(2,SymB).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                        EXIT
                    ENDIF
                ELSE
!We want an orbital from the same class. Check that this isn't the only unoccupied beta orbital in the class.
                    IF(ClassCountUnocc2(2,SymB).ne.1) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                        EXIT
                    ENDIF
                ENDIF
            ELSE
!We want an alpha/alpha pair.
                IF(SymProduct.ne.0) THEN
!Check to see if there are any unoccupied alpha orbitals in the SymB Class.
                    IF(ClassCountUnocc2(1,SymB).ne.0) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                        EXIT
                    ENDIF
                ELSE
!We want an orbital from the same class. Check that this isn't the only unoccupied alpha orbital in the class.
                    IF(ClassCountUnocc2(1,SymB).ne.1) THEN
!Success! We have found an allowed A orbital! Exit from loop.
                        EXIT
                    ENDIF
                ENDIF
            ENDIF

            IF(AttemptsOverall.gt.300) THEN
                WRITE(6,*) "Cannot find first allowed unoccupied orbital for given i,j pair after 300 attempts."
                WRITE(6,*) "It may be that there are no possible excitations from this i,j pair, in which case "
                WRITE(6,*) "the given algorithm is inadequate to describe excitations from such a small space."
                WRITE(6,*) "Try reverting to old excitation generators."
                CALL WRITEDET(6,nI,NEl,.true.)
                WRITE(6,*) "I,J pair; sym_i, sym_j: ",nI(Elec1Ind),nI(Elec2Ind),G1(nI(Elec1Ind))%Sym%S,G1(nI(Elec2Ind))%Sym%S
                CALL Stop_All("PickAOrb","Cannot find first allowed unocc orb for double excitation")
            ENDIF

        enddo

    END SUBROUTINE PickAOrb

!This routine takes determinant nI and returns two randomly chosen electrons, whose index in nI is Elec1Ind and Elec2Ind.
!These electrons have symmetry product SymProduct and spin pairing iSpn.
    SUBROUTINE PickElecPair(nI,Elec1Ind,Elec2Ind,iSeed,SymProduct,iSpn)
        INTEGER :: Ind,iSeed,X,K,Elec1Ind,Elec2Ind,SymProduct
        INTEGER :: nI(NEl),iSpn
        REAL*8 :: Ran2
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

        ElecPairs=(NEl*(NEl-1))/2

!Find an index randomly.
        Ind=INT(ElecPairs*Ran2(iSeed))+1

!X is number of elements at positions larger than ind
        X=ElecPairs-Ind
!K is the number of complete rows after the element ind
        K=INT((SQRT(8.D0*REAL(X,r2)+1.D0)-1.D0)/2.D0)
        Elec1Ind=NEl-1-K
        Elec2Ind=NEl-X+((K*(K+1))/2)

!We now want to find the symmetry product of the two electrons, and the spin product of the two electrons.
        SymProduct=INT(IEOR(G1(nI(Elec1Ind))%Sym%S,G1(nI(Elec2Ind))%Sym%S),4)

        IF((G1(nI(Elec1Ind))%Ms)*(G1(nI(Elec2Ind))%Ms).eq.-1) THEN
!We have an alpha beta pair of electrons.
            iSpn=2
        ELSE
            IF(G1(nI(Elec1Ind))%Ms.eq.-1) THEN
!We have an alpha alpha pair of electrons.
                iSpn=3
            ELSE
!We have a beta beta pair of electrons.
                iSpn=1
            ENDIF
        ENDIF

    END SUBROUTINE PickElecPair


    SUBROUTINE CreateSingleExcit(nI,nJ,iSeed,ClassCount2,ClassCountUnocc2,ILUT,ExcitMat,tParity,pGen)
        INTEGER :: ElecsWNoExcits,i,Attempts,nOrbs,z,Orb
        INTEGER :: Eleci,ElecSym,nI(NEl),nJ(NEl),NExcit,iSpn,ChosenUnocc
        INTEGER :: ExcitMat(2,2),ExcitLevel,iGetExcitLevel,iSeed
        INTEGER :: ClassCount2(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        INTEGER :: ILUT(0:nBasis/32)
        REAL*8 :: Ran2,pGen
        LOGICAL :: tParity,IsValidDet

!First, we need to find out if there are any electrons which have no possible excitations. This is because these will need to be redrawn and so 
!will affect the probabilities.
        ElecsWNoExcits=0
        do i=0,nSymLabels-1
!Run through all labels
            IF((ClassCount2(1,i).ne.0).and.(ClassCountUnocc2(1,i).eq.0)) THEN
!If there are alpha electrons in this class with no possible unoccupied alpha orbitals in the same class, these alpha electrons have no single excitations.
                ElecsWNoExcits=ElecsWNoExcits+ClassCount2(1,i)
            ENDIF
            IF((ClassCount2(2,i).ne.0).and.(ClassCountUnocc2(2,i).eq.0)) THEN
                ElecsWNoExcits=ElecsWNoExcits+ClassCount2(2,i)
            ENDIF
        enddo

        IF(ElecsWNoExcits.eq.NEl) THEN
!There are no single excitations from this determinant at all. Indicate this by putting pGen=-1.D0.
!Then we will create a double excitation instead.
            pGen=-1.D0
            RETURN
        ENDIF


        Attempts=0
        do while(.true.)

            Attempts=Attempts+1

!Choose an electron randomly...
            Eleci=INT(NEl*Ran2(iSeed))+1

!Find symmetry of chosen electron
            ElecSym=INT((G1(nI(Eleci))%Sym%S),4)

            IF(G1(nI(Eleci))%Ms.eq.-1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                NExcit=ClassCountUnocc2(1,ElecSym)
                iSpn=2
            ELSE
!Beta orbital
                NExcit=ClassCountUnocc2(2,ElecSym)
                iSpn=1
            ENDIF

            IF(NExcit.ne.0) EXIT    !Have found electron with allowed excitations

            IF(Attempts.gt.250) THEN
                WRITE(6,*) "Cannot find single excitation from electrons after 250 attempts..."
                CALL WRITEDET(6,nI,NEL,.true.)
                WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
                WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
                WRITE(6,*) "***"
                WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
                WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
                CALL Stop_All("CreateSingleExcit","Cannot find single excitation from following electron after 250 attempts...")
            ENDIF

        enddo

!Now we need to choose the unoccupied orbital for the chosen electron.
!There are two ways to do this. We can either choose the orbital we want out of the NExcit possible unoccupied orbitals.
!It would then be necessary to cycle through all orbitals of that symmetry and spin, only counting the unoccupied ones to find the correct determinant.
!Alternatively, we could pick orbitals at random and redraw until we find an allowed one. This would probably be preferable for larger systems.
        IF(tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals in symmetry class with needed spin to find allowed one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
            ChosenUnocc=INT(NExcit*Ran2(iSeed))+1

!Now run through all allowed orbitals until we find this one.
            nOrbs=SymLabelCounts(2,ElecSym+1)
            z=0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
            do i=0,nOrbs-1
!Find the spin orbital index. SymLabelCounts has the index of the state for the given symmetry.
                Orb=(2*SymLabelList(SymLabelCounts(1,ElecSym+1)+i))-(iSpn-1)

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

            nOrbs=SymLabelCounts(2,ElecSym+1)
            Attempts=0
            do while(.true.)
                Attempts=Attempts+1
                
!Draw randomly from the set of orbitals
                ChosenUnocc=INT(nOrbs*Ran2(iSeed))
                Orb=(2*SymLabelList(SymLabelCounts(1,ElecSym+1)+ChosenUnocc))-(iSpn-1)

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
                    WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
                    WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
                    WRITE(6,*) "***"
                    WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
                    WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
                    CALL Stop_All("CreateSingleExcit","Cannot find single excitation unoccupied orbital after 250 attempts...")
                ENDIF

            enddo

        ENDIF

        nJ(:)=nI(:)
!ExcitMat wants to be the index in nI of the orbital to excite from, but returns the actual orbitals.
        ExcitMat(1,1)=Eleci
        ExcitMat(2,1)=Orb
        CALL FindExcitDet(ExcitMat,nJ,1,TParity)

!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
        CALL IsSymAllowedExcit(nI,nJ,1,ExcitMat)

!Now we need to find the probability of creating this excitation.
!This is: P_single x P(i) x P(a|i) x N/(N-ElecsWNoExcits)
        pGen=(1.D0-pDoubNew)*(1.D0/real(NEl,r2))*(1.D0/real(NExcit,r2))*((real(NEl,r2))/(real((NEl-ElecsWNoExcits),r2)))

    END SUBROUTINE CreateSingleExcit


    SUBROUTINE ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2,ILUT)
        INTEGER :: i,nI(NEl)
        INTEGER :: ClassCount2(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        INTEGER :: ILUT(0:nBasis/32)

        ClassCount2(:,:)=0
        ClassCountUnocc2(:,:)=0
        ILUT(:)=0
        NOccAlpha=0
        NOccBeta=0
        do i=1,NEl
            
!Create ILUT for O[1] comparison of orbitals in root determinant
            ILUT((nI(I)-1)/32)=IBSET(ILUT((NI(I)-1)/32),MOD(NI(I)-1,32))

            IF(G1(nI(i))%Ms.eq.-1) THEN
!orbital is an alpha orbital and symmetry of the orbital can be found in G1
!                WRITE(6,*) G1(nI(i))%Ms,G1(nI(i))%Sym%S
                ClassCount2(1,INT(G1(nI(i))%Sym%S,4))=ClassCount2(1,INT(G1(nI(i))%Sym%S,4))+1
                NOccAlpha=NOccAlpha+1

            ELSE
!orbital is a beta orbital
!                WRITE(6,*) G1(nI(i))%Ms,G1(nI(i))%Sym%S
                ClassCount2(2,INT(G1(nI(i))%Sym%S,4))=ClassCount2(2,INT(G1(nI(i))%Sym%S,4))+1
                NOccBeta=NOccBeta+1
            ENDIF
        enddo

!We now want to find ClassCountUnocc2 - the unoccupied version of the array
!SymLabelCounts(2,1:nSymLabels) gives the number of *states* in each symmetry class.
!There are therefore equal number of alpha and beta orbitals in each state from which to calculate the unoccupied classcount.
!Again, we store as 
        do i=1,nSymLabels
            ClassCountUnocc2(1,i-1)=SymLabelCounts(2,i)-ClassCount2(1,i-1)
            ClassCountUnocc2(2,i-1)=SymLabelCounts(2,i)-ClassCount2(2,i-1)
        enddo

    END SUBROUTINE ConstructClassCounts

END MODULE GenRandSymExcitNUMod

!This routine will take a determinant nI, and find Iterations number of excitations of it. It will then histogram these, summing in 1/pGen for every occurance of
!the excitation. This means that all excitations should be 0 or 1 after enough iterations. It will then count the excitations and compare the number to the
!number of excitations generated using the full enumeration excitation generation. This can be done for both doubles and singles, or one of them.
SUBROUTINE TestGenRandSymExcitNU(nI,Iterations,pDoub,iSeed,exFlag)
    Use SystemData , only : NEl,nBasis,G1,nBasisMax
    Use GenRandSymExcitNUMod , only : GenRandSymExcitNU,ConstructClassCounts
    Use SymData , only : nSymLabels
    IMPLICIT NONE
    INTEGER :: i,Iterations,exFlag,nI(NEl),nJ(NEl),iSeed,IC,ExcitMat(2,2),DetConn
    REAL*8 :: pDoub,pGen
    INTEGER :: ClassCount2(2,0:nSymLabels-1)
    INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
    INTEGER :: ILUT(0:nBasis/32)
    LOGICAL :: tParity
    REAL*8 , ALLOCATABLE :: DoublesHist(:,:,:,:),SinglesHist(:,:)
    INTEGER , ALLOCATABLE :: EXCITGEN(:)
    INTEGER :: ierr,Ind1,Ind2,Ind3,Ind4,iMaxExcit,nStore(6),nExcitMemLen,j,k,l,DetNum,DetNumS

    WRITE(6,*) nI(:)
    WRITE(6,*) Iterations,pDoub,iSeed,exFlag
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
    CALL GetSymExcitCount(EXCITGEN,DetConn)
    WRITE(6,*) "Determinant has ",DetConn," total excitations from it."
    CALL FLUSH(6)

!Allocate memory for histogramming determinants
    ALLOCATE(DoublesHist(nBasis,nBasis,nBasis,nBasis),stat=ierr)
    IF(ierr.ne.0) THEN
        CALL Stop_All("TestGenRandSymExcitNU","Not possible to allocate memory to do histogramming")
    ENDIF
    ALLOCATE(SinglesHist(nBasis,nBasis),stat=ierr)
    IF(ierr.ne.0) THEN
        CALL Stop_All("TestGenRandSymExcitNU","Not possible to allocate memory to do histogramming")
    ENDIF
    DoublesHist(:,:,:,:)=0.D0
    SinglesHist(:,:)=0.D0
    
    do i=1,Iterations

        WRITE(6,"(A,I10)") "Iteration: ",i

        CALL GenRandSymExcitNU(nI,nJ,iSeed,pDoub,IC,ExcitMat,TParity,exFlag,pGen)
        IF(IC.eq.1) THEN
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(2,1)
        ELSE
!            WRITE(6,*) "Double Created"
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(1,2),ExcitMat(2,1),ExcitMat(2,2)
        ENDIF

        IF(IC.eq.1) THEN
            SinglesHist(ExcitMat(1,1),ExcitMat(2,1))=SinglesHist(ExcitMat(1,1),ExcitMat(2,1))+(1.D0/pGen)
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

!Check excitation
        CALL IsSymAllowedExcit(nI,nJ,IC,ExcitMat)

    enddo

!Now run through arrays normalising them so that numbers are more managable.
    OPEN(8,FILE="DoublesHist",STATUS="UNKNOWN")
    DetNum=0
    do i=1,nBasis-1
        do j=i+1,nBasis
            do k=1,nBasis-1
                do l=k+1,nBasis
                    IF(DoublesHist(i,j,k,l).gt.0.D0) THEN
!                        DoublesHist(i,j,k,l)=DoublesHist(i,j,k,l)/real(Iterations,8)
                        DetNum=DetNum+1
                        WRITE(8,"(I12,F20.12,4I5)") DetNum,DoublesHist(i,j,k,l)/real(Iterations,8),i,j,k,l
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
            IF(SinglesHist(i,j).gt.0.D0) THEN
                DetNumS=DetNumS+1
                WRITE(9,*) DetNumS,SinglesHist(i,j)/real(Iterations,8)
            ENDIF
        enddo
    enddo
    CLOSE(9)
    WRITE(6,*) DetNumS," Single excitations found from nI"
    IF((DetNum+DetNumS).ne.DetConn) THEN
        CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2,ILUT)
        WRITE(6,*) "Total determinants = ", DetConn
        WRITE(6,*) "ClassCount2(1,:)= ",ClassCount2(1,:)
        WRITE(6,*) "ClassCount2(2,:)= ",ClassCount2(2,:)
        WRITE(6,*) "***"
        WRITE(6,*) "ClassCountUnocc2(1,:)= ",ClassCountUnocc2(1,:)
        WRITE(6,*) "ClassCountUnocc2(2,:)= ",ClassCountUnocc2(2,:)
        CALL Stop_All("TestGenRandSymExcitNU","Not all excitations accounted for...")
    ENDIF

END SUBROUTINE TestGenRandSymExcitNU

SUBROUTINE IsSymAllowedExcit(nI,nJ,IC,ExcitMat)
    Use SystemData , only : G1,NEl
    Use SystemData , only : Symmetry
    IMPLICIT NONE
    Type(Symmetry) :: SymProduct,SymProduct2,SYMPROD
    LOGICAL :: SYMEQ,ISVALIDDET
    INTEGER :: IC,ExcitMat(2,2),nI(NEl),nJ(NEl),ExcitLevel,iGetExcitLevel

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
     
     IF(IC.eq.2) THEN
        SymProduct=SYMPROD(G1(ExcitMat(1,1))%Sym,G1(ExcitMat(1,2))%Sym)
        SymProduct2=SYMPROD(G1(ExcitMat(2,1))%Sym,G1(ExcitMat(2,2))%Sym)
        IF(.not.SYMEQ(SymProduct,SymProduct2)) THEN
            CALL Stop_All("IsSymAllowedExcit","Excitation not a valid symmetry allowed double excitation")
        ENDIF
    ELSE
        IF(.not.SYMEQ(G1(ExcitMat(1,1))%Sym,G1(ExcitMat(2,1))%Sym)) THEN
            CALL Stop_All("IsSymAllowedExcit","Excitation not a valid symmetry allowed single excitation")
        ENDIF
        IF(G1(ExcitMat(1,1))%Ms.ne.G1(ExcitMat(2,1))%Ms) THEN
            CALL Stop_All("IsSymAllowedExcit","Excitation not a valid spin-symmetry allowed single excitation")
        ENDIF
    ENDIF
        

END SUBROUTINE IsSymAllowedExcit
