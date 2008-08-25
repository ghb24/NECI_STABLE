MODULE FciMCMod
    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr
    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    USE Calc , only : TReadPops,ScaleWalkers,TMCExcitSpace,NoMCExcits,TStartMP1
    USE Calc , only : GrowMaxFactor,CullFactor,TMCDets,TNoBirth,Lambda,TDiffuse,FlipTauCyc,TFlipTau
    USE Calc , only : TExtraPartDiff,TFullUnbias,TNodalCutoff,NodalCutoff,TNoAnnihil,TMCDiffusion
    USE Calc , only : NDets,RhoApp,TResumFCIMC,NEquilSteps,TSignShift,THFRetBias,PRet,TExcludeRandGuide
    USE Calc , only : TProjEMP2,TFixParticleSign,TStartSinglePart
    USE Determinants , only : FDet,GetHElement2,GetH0Element3
    USE DetCalc , only : NMRKS,ICILevel
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE Logging , only : iWritePopsEvery,TPopsFile
    USE HElem
    IMPLICIT NONE
    SAVE

    INTEGER, PARAMETER :: r2=kind(0.d0)

    TYPE ExcitGenerator
        INTEGER , ALLOCATABLE :: ExcitData(:)      !This stores the excitation generator
        INTEGER :: nExcitMemLen                    !This is the length of the excitation generator
        LOGICAL :: ExitGenForDet=.false.           !This is true when the excitation generator stored corresponds to the determinant
    END TYPE

    TYPE(ExcitGenerator) , ALLOCATABLE , TARGET :: WalkVecExcits(:),WalkVec2Excits(:)   !This will store the excitation generators for the particles on each node
    INTEGER , ALLOCATABLE , TARGET :: WalkVecDets(:,:),WalkVec2Dets(:,:)    !Contains determinant list
    LOGICAL , ALLOCATABLE , TARGET :: WalkVecSign(:),WalkVec2Sign(:)        !Contains sign list
    INTEGER , ALLOCATABLE , TARGET :: WalkVecIC(:),WalkVec2IC(:)            !Contains excit level list
    REAL*8 , ALLOCATABLE , TARGET :: WalkVecH(:,:),WalkVec2H(:,:)       !First element is diagonal hamiltonian element - second is the connection to HF determinant
    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0
    INTEGER :: WalkVecICTag=0,WalkVec2ICTag=0,WalkVecHTag=0,WalkVec2HTag=0
!Pointers to point at the correct arrays for use
    INTEGER , POINTER :: CurrentDets(:,:), NewDets(:,:)
    LOGICAL , POINTER :: CurrentSign(:), NewSign(:)
    INTEGER , POINTER :: CurrentIC(:), NewIC(:)
    REAL*8 , POINTER :: CurrentH(:,:), NewH(:,:)
    TYPE(ExcitGenerator) , POINTER :: CurrentExcits(:), NewExcits(:)

    INTEGER , ALLOCATABLE :: HFDet(:)       !This will store the HF determinant
    INTEGER :: HFDetTag=0
    TYPE(ExcitGenerator) :: HFExcit         !This is the excitation generator for the HF determinant

!MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
    INTEGER :: MemoryFac=200

    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,TotSign,TotSignOld,PreviousNMCyc,Iter,NoComps
    INTEGER :: exFlag=3

!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle for each variable
    INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
    INTEGER :: CullInfo(10,3)

!The following variables are calculated at the end of each update cycle, are combined to the root processor
    REAL*8 :: GrowRate,DieRat,ProjectionE,SumENum
    INTEGER*8 :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
    REAL*8 :: PosFrac           !This is the fraction of positive particles on each node
    INTEGER :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
    REAL*8 :: MeanExcitLevel
    INTEGER :: MinExcitLevel
    INTEGER :: MaxExcitLevel
    INTEGER :: NoatDoubs,Annihilated,Acceptances
    REAL*8 :: AccRat
    INTEGER :: PreviousCycles   !The number of previous cycles performed before the POPSFILE is read in
    INTEGER :: NoatHF           !This is the number at HF for a given Iteration

!These values are used to calculate the energy when TProjEMP2 is set
    REAL*8 :: SumOverlapMP2     !This is the overlap of the walker distribution with the MP2 wavefunction, summed over all iterations
    REAL*8 :: SumHOverlapMP2    !This is the overlap of the excitations of the walkers with MP2, summed over all iterations
    REAL*8 :: ProjectionMP2     !This is the energy calculated by <Psi_MP1|H|Psi>/<Psi_MP1|Psi>

    TYPE(HElement) :: rhii,FZero
    REAL*8 :: Hii,Fii
    
    REAL*8 , ALLOCATABLE :: GraphRhoMat(:,:)    !This stores the rho matrix for the graphs in resumFCIMC
    INTEGER :: GraphRhoMatTag=0

    REAL*8 , ALLOCATABLE :: GraphVec(:)         !This stores the final components for the propagated graph in ResumFCIMC
    INTEGER :: GraphVecTag=0

    REAL*8 , ALLOCATABLE :: GraphKii(:)         !This stores the diagonal Kii matrix elements for the determinants in the graph
    INTEGER :: GraphKiiTag=0

    INTEGER , ALLOCATABLE :: DetsinGraph(:,:)   !This stores the determinants in the graph created for ResumFCIMC
    INTEGER :: DetsinGraphTag=0

    LOGICAL :: TSinglePartPhase                 !This is true if TStartSinglePart is true, and we are still in the phase where the shift is fixed and particle numbers are growing

    contains

    SUBROUTINE FciMC(Weight,Energyxw)
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j
        LOGICAL :: exists,TestSoftExit
        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'
        TYPE(HElement) :: Hamii

        CALL InitFCIMCCalc()

        WRITE(6,*) ""
        WRITE(6,*) "       Step     Shift   WalkerCng   GrowRate    TotWalkers  Annihil   Proj.E        Proj.MP2    SumNoatHF NoatDoubs  +veWalk        AccRat     MeanExcit  MinExcit MaxEx"
        WRITE(15,*) "#       Step     Shift   WalkerCng   GrowRate    TotWalkers  Annihil   Proj.E        Proj.MP2    SumNoatHF NoatDoubs  +veWalk        AccRat     MeanExcit  MinExcit MaxEx"

!TotWalkersOld is the number of walkers last time the shift was changed
        WRITE(15,"(I12,G15.6,I7,G15.6,I10,I6,2G15.6,I11,I9,3G14.6,2I6)") PreviousCycles+Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,Annihilated,ProjectionE,ProjectionMP2,SumNoatHF,NoatDoubs,1.D0,AccRat,MeanExcitLevel,MaxExcitLevel,MinExcitLevel
        WRITE(6,"(I12,G15.6,I7,G15.6,I10,I6,2G15.6,I11,I9,3G14.6,2I6)") PreviousCycles+Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,Annihilated,ProjectionE,ProjectionMP2,SumNoatHF,NoatDoubs,1.D0,AccRat,MeanExcitLevel,MaxExcitLevel,MinExcitLevel

!Start MC simulation...
        do Iter=1,NMCyc

            CALL PerformFCIMCyc()

!Test that the file SOFTEXIT is still present. If not, then exit cleanly.
            exists=TestSoftExit()
            IF(.not.exists) EXIT

            IF(mod(Iter,StepsSft).eq.0) THEN
!This will find the new shift (and other parameters), and print out the needed values.
                CALL UpdateDiagSft()
            ENDIF

            IF(TPopsFile.and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
!This will write out the configuration of walkers
                CALL WriteToPopsFile()
            ENDIF

!End of MC cycle
        enddo

!Write out popsfile
        Iter=Iter-1
        IF(TPopsFile) CALL WriteToPopsFile()

        Weight=HDElement(0.D0)
        Energyxw=HDElement(ProjectionE)

!Deallocate memory
        CALL DeallocFCIMCMem()

        CLOSE(15)

        RETURN

    END SUBROUTINE FciMC


!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCyc()
        INTEGER :: VecSlot,i,j,k,l
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount,TotWalkersNew
        REAL*8 :: Prob,rat,HDiag,Ran2,TotProb
        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,iGetExcitLevel_2,MaxExcits
        LOGICAL :: WSign,SpawnBias,SameDet,TGenGuideDet
        TYPE(HElement) :: HDiagTemp,HOffDiag

!VecSlot indicates the next free position in NewDets
        VecSlot=1
        TotSign=0

        do j=1,TotWalkers
!j runs through all current walkers

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(CurrentDets(:,j),CurrentIC(j),CurrentH(2,j),CurrentSign(j),j)

            IF(TResumFCIMC) THEN
!Setup excit generators for this determinant
                CALL SetupExitgen(CurrentDets(:,j),CurrentExcits(j))
                CALL ResumGraph(CurrentDets(:,j),CurrentSign(j),VecSlot,CurrentExcits(j),j)

            ELSEIF(TFixParticleSign) THEN

                CALL FixParticleSignIter(j,VecSlot)

            ELSEIF(THFRetBias) THEN

                CALL HFRetBiasIter(j,VecSlot)


            ELSE    !Not using HFRetBias - normal spawn
!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
                CALL SetupExitgen(CurrentDets(:,j),CurrentExcits(j))
                CALL GenRandSymExcitIt3(CurrentDets(:,j),CurrentExcits(j)%ExcitData,nJ,Seed,IC,0,Prob,iCount)

                IF(ICILevel.ne.0) THEN
!We are performing run in a truncated space
                    IF(CurrentIC(j).eq.(ICILevel-1)) THEN
!The current walker is one below the excitation cutoff - if IC is a double, then could go over
                        IF(IC.eq.2) THEN
!Need to check excitation level of excitation
                            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                            IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                                Child=0
                            ELSE
                                Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC)
                            ENDIF
                        ELSE
                            Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC)
                        ENDIF
                    ELSEIF(CurrentIC(j).eq.ICILevel) THEN
!Walker is at the excitation cutoff level - all possible excitations could be disallowed
                        ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                        IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                            Child=0
                        ELSE
                            Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC)
                        ENDIF
                    ELSE
!Excitation cannot be in a dissallowed excitation level - allow it as normal
                        Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC)
                    ENDIF

                ELSE

!Calculate number of children to spawn
                    Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC)
                ENDIF

                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store
                    IF(Child.gt.0) THEN
                        WSign=.true.    !+ve child created
                    ELSE
                        WSign=.false.   !-ve child created
                    ENDIF
!Calculate excitation level, connection to HF and diagonal ham element
                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
                    IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
                        HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                    ENDIF
                    IF(ExcitLevel.eq.0) THEN
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                    ELSE
                        HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                        HDiag=(REAL(HDiagTemp%v,r2))-Hii
                    ENDIF

                    do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them
                        NewDets(:,VecSlot)=nJ(:)
                        NewSign(VecSlot)=WSign
                        NewExcits(VecSlot)%ExitGenForDet=.false.
                        NewIC(VecSlot)=ExcitLevel
                        NewH(1,VecSlot)=HDiag      !Diagonal H-element-Hii
                        NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                        VecSlot=VecSlot+1
                    enddo
                            
                    Acceptances=Acceptances+ABS(Child)  !Sum the number of created children to use in acceptance ratio

                ENDIF   !End if child created

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
                iDie=AttemptDie(CurrentDets(:,j),CurrentH(1,j))
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births

                IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle

                    do l=1,abs(iDie)+1    !We need to copy accross one more, since we need to include the original spared particle
                        NewDets(:,VecSlot)=CurrentDets(:,j)
                        NewSign(VecSlot)=CurrentSign(j)
!Copy excitation generator accross
                        CALL CopyExitgen(CurrentExcits(j),NewExcits(VecSlot))
                        NewIC(VecSlot)=CurrentIC(j)
                        NewH(:,VecSlot)=CurrentH(:,j)
                        VecSlot=VecSlot+1
                    enddo

                ELSEIF(iDie.gt.0) THEN
!This indicates that particles want to be killed. The first kill will simply be performed by not copying accross the original particle.
!Therefore, if iDie = 1, then we can simply ignore it.
!However, after that anti-particles will need to be created on the same determinant.

                    do l=1,iDie-1
                        NewDets(:,VecSlot)=CurrentDets(:,j)
                        IF(CurrentSign(j)) THEN
!Copy accross new anti-particles
                            NewSign(VecSlot)=.FALSE.
                        ELSE
                            NewSign(VecSlot)=.TRUE.
                        ENDIF
!Copy excitation generator accross
                        CALL CopyExitgen(CurrentExcits(j),NewExcits(VecSlot))
                        NewIC(VecSlot)=CurrentIC(j)
                        NewH(:,VecSlot)=CurrentH(:,j)
                        VecSlot=VecSlot+1
                    enddo

                ENDIF   !To kill if

            ENDIF   !Choice of method (Resum/FixParticleSign/HFRetBias) if

!Finish cycling over walkers
        enddo

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+TotWalkers

!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1
        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
        ENDIF

        IF(TNoAnnihil) THEN

            TotWalkers=TotWalkersNew    !Since there is no further cancellation of walkers

!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
            IF(associated(CurrentDets,target=WalkVecDets)) THEN
                CurrentDets=>WalkVec2Dets
                CurrentSign=>WalkVec2Sign
                CurrentIC=>WalkVec2IC
                CurrentH=>WalkVec2H
                CurrentExcits=>WalkVec2Excits
                NewDets=>WalkVecDets
                NewSign=>WalkVecSign
                NewIC=>WalkVecIC
                NewH=>WalkVecH
                NewExcits=>WalkVecExcits
            ELSE
                CurrentDets=>WalkVecDets
                CurrentSign=>WalkVecSign
                CurrentIC=>WalkVecIC
                CurrentH=>WalkVecH
                CurrentExcits=>WalkVecExcits
                NewDets=>WalkVec2Dets
                NewSign=>WalkVec2Sign
                NewIC=>WalkVec2IC
                NewH=>WalkVec2H
                NewExcits=>WalkVec2Excits
            ENDIF

        ELSE

!This routine now cancels down the particles with opposing sign on each determinant
!This routine does not necessarily need to be called every Iter
            CALL AnnihilatePairs(TotWalkersNew)
            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)
!            WRITE(6,*) "Number of annihilated particles= ",TotWalkersNew-TotWalkers,Iter,TotWalkers
        ENDIF

        IF(NoatHF.lt.0) THEN
!Flip the sign if we're beginning to get a negative population on the HF
            WRITE(6,*) "No. at HF < 0 - flipping sign of entire ensemble of particles..."
            CALL FlipSign()
        ENDIF
!Reset the number at HF per iteration
        NoatHF=0


        IF(TSinglePartPhase) THEN
!Do not allow culling if we are still in the single particle phase.
            IF(TotWalkers.gt.InitWalkers) THEN
                WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
                TSinglePartPhase=.false.
            ENDIF
        ELSE
!Check whether we need to cull particles
            IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!Particle number is too large - kill them randomly

!Log the fact that we have made a cull
                NoCulls=NoCulls+1
                IF(NoCulls.gt.10) THEN
                    WRITE(6,*) "Too Many Culls"
                    CALL FLUSH(6)
                    CALL Stop_All("FCIMC","Too Many Culls")
                ENDIF

                IF(TSignShift) THEN
!CullInfo(:,1) is walkers/residualsign before cull
                    CullInfo(NoCulls,1)=TotSign
                ELSE
                    CullInfo(NoCulls,1)=TotWalkers
                ENDIF
                IF(mod(Iter,StepsSft).eq.0) THEN
!CullInfo(:,3) is MC Steps into shift cycle before cull
!This is just before the calculation of the shift - we want the value to be equal to StepsSft
                    CullInfo(NoCulls,3)=StepsSft
                ELSE
                    CullInfo(NoCulls,3)=mod(Iter,StepsSft)
                ENDIF

                WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number..."
                WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number..."
                WRITE(6,"(A,F8.2)") "Population will reduce by a factor of ",CullFactor
                CALL ThermostatParticles(.true.)

            ELSEIF(TotWalkers.lt.(InitWalkers/2)) THEN
!Particle number is too small - double every particle in its current position

!Log the fact that we have made a cull
                NoCulls=NoCulls+1
                IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
                IF(TSignShift) THEN
!CullInfo(:,1) is walkers/residualsign before cull
                    CullInfo(NoCulls,1)=TotSign
                ELSE
                    CullInfo(NoCulls,1)=TotWalkers
                ENDIF
                IF(mod(Iter,StepsSft).eq.0) THEN
!CullInfo(:,3) is MC Steps into shift cycle before cull
!This is just before the calculation of the shift - we want the value to be equal to StepsSft
                    CullInfo(NoCulls,3)=StepsSft
                ELSE
                    CullInfo(NoCulls,3)=mod(Iter,StepsSft)
                ENDIF

                WRITE(6,*) "Doubling particle population to increase total number..."
                CALL ThermostatParticles(.false.)

            ENDIF

        ENDIF

        RETURN

    END SUBROUTINE PerformFCIMCyc
                
    
    SUBROUTINE HFRetBiasIter(j,VecSlot)
        LOGICAL :: SpawnBias,WSign,SameDet,TGenGuideDet
        REAL*8 :: TotProb,Ran2,Prob,HDiag
        INTEGER :: j,VecSlot,Child,iCount,IC,nJ(NEl),MaxExcits,ExcitLevel,iGetExcitLevel_2,iDie,l
        TYPE(HElement) :: HOffDiag,HDiagTemp

!We want the simplest guiding function - if we're at a double, attempt to spawn at HF with PRet probability
        SpawnBias=.false.
        TotProb=1.D0  
            
        IF(CurrentIC(j).eq.2) THEN
!We are at a double - see if we are forced to return to HF
            IF(Ran2(Seed).lt.PRet) THEN
!Ensure that we try to spawn children at HF -   Modify prob of doing this to equal PRet
                SpawnBias=.true.
                IF(TExcludeRandGuide) THEN
!In this method of unbiasing the Guiding function, we unbias completely at this stage
                    Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),HFDet,PRet,2,CurrentH(2,j))
                ELSE
!No need to unbias here - divide by 1, not PRet
                    Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),HFDet,1.D0,2,CurrentH(2,j))
                ENDIF
                IF(Child.ne.0) THEN
                    IF(Child.gt.0) THEN
                        WSign=.true. !We have successfully created at least one positive child at HF
                    ELSE
                        WSign=.false. !We have successfully created at least one negative child at HF
                    ENDIF

                    do l=1,abs(Child)
!                                IF(abs(Child).gt.1) WRITE(6,*) "Multiple children created (returned)"
!Copy across children to HF - can also copy across known excitation generator
                        NewDets(:,VecSlot)=HFDet(:)
                        NewSign(VecSlot)=WSign
                        CALL CopyExitgen(HFExcit,NewExcits(VecSlot))
                        NewIC(VecSlot)=0
                        NewH(1,VecSlot)=0.D0
                        NewH(2,VecSlot)=0.D0
                        VecSlot=VecSlot+1
                    enddo

                    Acceptances=Acceptances+ABS(Child)  !Sum the number of created children to use in acceptance ratio

                ENDIF   !End if child created

            ELSE
                    
                TotProb=1.D0-PRet
                    
            ENDIF   !End if Spawning back at HF due to bias

        ENDIF

        IF(.not.SpawnBias) THEN
!We are either at a double but have decided not to try to spawn back to HF, or we are not at a double

!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry).
            CALL SetupExitgen(CurrentDets(:,j),CurrentExcits(j))
            CALL GenRandSymExcitIt3(CurrentDets(:,j),CurrentExcits(j)%ExcitData,nJ,Seed,IC,0,Prob,iCount)

            !TESTS
            CALL GetSymExcitCount(CurrentExcits(j)%ExcitData,MaxExcits)
            IF(ABS((1.D0/MaxExcits)-Prob).gt.1.D-07) THEN
                WRITE(6,*) "PROBLEM WITH PGENS!"
                WRITE(6,*) MaxExcits,1.D0/MaxExcits,Prob
                CALL Stop_All("PerformFCIMCyc","Problem with PGens")
            ENDIF
            IF((CurrentIC(j).eq.2).and.(TotProb.ne.(1.D0-PRet))) WRITE(6,*) "PROBLEM HERE!!"
            IF((CurrentIC(j).ne.2).and.(TotProb.ne.1.D0)) WRITE(6,*) "PROBLEM HERE!"


            IF(TExcludeRandGuide) THEN
!We must not be allowed to generate a guiding determinant from a double excitation
        
                IF(CurrentIC(j).eq.2) THEN 

                    IF(SameDet(HFDet,nJ,NEl)) THEN
                        TGenGuideDet=.true.    !Have not generated a guiding determinant (HF)
                    ELSE
                        TGenGuideDet=.false.     !Have generated a guiding determinant (HF)
                    ENDIF
                    do while(TGenGuideDet)
!Generate another randomly connected determinant in order to not generate a guiding determinant

                        CALL GenRandSymExcitIt3(CurrentDets(:,j),CurrentExcits(j)%ExcitData,nJ,Seed,IC,0,Prob,iCount)
                        IF(.not.(SameDet(HFDet,nJ,NEl))) TGenGuideDet=.false.
                            
                    enddo   !Have now definitely generated a non-HF determinant
!Need to change the probabilities of generating these excitations, since there is now a determinant which is forbidden

                    Prob=1.D0/(MaxExcits-1)      !Prob is now 1/(N-1) - do not allow excit weighting
                    TotProb=TotProb*Prob
                    
                ELSE
!We are not at a double - want to unbias with normal prob
                    TotProb=Prob
                        
                ENDIF

            ELSE
!The other unbiasing method allows the guiding determinant to be generated, but we unbias it differently.

                IF(CurrentIC(j).eq.2) THEN
                    IF(SameDet(HFDet,nJ,NEl)) THEN
!We are at a double, and have decided to randomly attempt a return to the guiding function HF determinant, so we need to change 
!the Pgen - it wants to be unbiased by dividing by just PGen, not PGen*(1-PRet).
                        TotProb=Prob
                    ELSE
                        TotProb=TotProb*Prob    !TotProb should initially by 1, or 1-PRet if we are at a double
                    ENDIF

                ELSE

                    TotProb=Prob
                    
                ENDIF
                
            ENDIF !Choice of unbiasing methods

!Calculate number of children to spawn
            Child=AttemptCreate(CurrentDets(:,j),CurrentSign(j),nJ,TotProb,IC)

            IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store
                IF(Child.gt.0) THEN
                    WSign=.true.    !+ve child created
                ELSE
                    WSign=.false.   !-ve child created
                ENDIF
!Calculate excitation level, connection to HF and diagonal ham element
                ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
                IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
                    HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                ENDIF
                IF(ExcitLevel.eq.0) THEN
!We know we are at HF - HDiag=0
                    HDiag=0.D0
                ELSE
                    HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                    HDiag=(REAL(HDiagTemp%v,r2))-Hii
                ENDIF

                do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them
                    NewDets(:,VecSlot)=nJ(:)
                    NewSign(VecSlot)=WSign
                    NewExcits(VecSlot)%ExitGenForDet=.false.
                    NewIC(VecSlot)=ExcitLevel
                    NewH(1,VecSlot)=HDiag      !Diagonal H-element-Hii
                    NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                    VecSlot=VecSlot+1
                enddo
                    
                Acceptances=Acceptances+ABS(Child)  !Sum the number of created children to use in acceptance ratio

            ENDIF   !End if child created

        ENDIF   !End if we are trying to create a child which isn't a returning spawn
!We now have to decide whether the parent particle (j) wants to self-destruct or not...
        iDie=AttemptDie(CurrentDets(:,j),CurrentH(1,j))
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births

        IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle

            do l=1,abs(iDie)+1    !We need to copy accross one more, since we need to include the original spared particle
                NewDets(:,VecSlot)=CurrentDets(:,j)
                NewSign(VecSlot)=CurrentSign(j)
!Copy excitation generator accross
                CALL CopyExitgen(CurrentExcits(j),NewExcits(VecSlot))
                NewIC(VecSlot)=CurrentIC(j)
                NewH(:,VecSlot)=CurrentH(:,j)
                VecSlot=VecSlot+1
            enddo

        ELSEIF(iDie.gt.0) THEN
!This indicates that particles want to be killed. The first kill will simply be performed by not copying accross the original particle.
!Therefore, if iDie = 1, then we can simply ignore it.
!However, after that anti-particles will need to be created on the same determinant.

            do l=1,iDie-1
                NewDets(:,VecSlot)=CurrentDets(:,j)
                IF(CurrentSign(j)) THEN
!Copy accross new anti-particles
                    NewSign(VecSlot)=.FALSE.
                ELSE
                    NewSign(VecSlot)=.TRUE.
                ENDIF
!Copy excitation generator accross
                CALL CopyExitgen(CurrentExcits(j),NewExcits(VecSlot))
                NewIC(VecSlot)=CurrentIC(j)
                NewH(:,VecSlot)=CurrentH(:,j)
                VecSlot=VecSlot+1
            enddo

        ENDIF   !To kill if
        RETURN
    END SUBROUTINE HFRetBiasIter

!This impliments a fixed sign approximation, whereby connections in the system which would flip the sign of the 
!particle are resummed into the diagonal matrix element contributing to the on-site death rate. 
!For a real-space analogue, see: D.F.B. ten Haaf, H.J.M. van Bemmel, J.M.J. van Leeuwen, W. van Saarloos, and D.M.
!Ceperley, Phys. Rev. B 51, 13039 (1995)
    SUBROUTINE FixParticleSignIter(Walker,VecSlot)
        INTEGER :: Walker,VecSlot,nStore(6),nJ(NEl),IC,Child,ExcitLevel,iGetExcitLevel_2
        INTEGER :: iDie,l
        TYPE(HElement) :: HConn,HOffDiag,HDiagTemp
        LOGICAL :: WSign
        REAL*8 :: SpinFlipContrib,HConnReal,HDiag,DiagDeath

!Setup excitgen
        CALL SetupExitgen(CurrentDets(:,Walker),CurrentExcits(Walker))
        CALL ResetExit2(CurrentDets(:,Walker),NEl,G1,nBasis,nBasisMax,CurrentExcits(Walker)%ExcitData,0) !Reset excitgen

!Run through all excits
        do while(.true.)
            CALL GenSymExcitIt2(CurrentDets(:,Walker),NEl,G1,nBasis,nBasisMax,.false.,CurrentExcits(Walker)%ExcitData,nJ,IC,0,nStore,exFlag)
            IF(nJ(1).eq.0) EXIT
            HConn=GetHElement2(CurrentDets(:,Walker),nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
            HConnReal=REAL(HConn%v,r2)
            IF(HConnReal.lt.0.D0) THEN
!Connection is negative, therefore attempt to create a particle there - we are running through all walkers, so Prob=1.D0
                Child=AttemptCreate(CurrentDets(:,Walker),CurrentSign(Walker),nJ,1.D0,IC,HConnReal)

                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store
                    IF(Child.gt.0) THEN
                        WSign=.true.    !+ve child created
                    ELSE
                        WSign=.false.   !-ve child created
                        IF(CurrentSign(Walker)) THEN
                            CALL Stop_All("FixParticleSignIter","Particle of opposite sign created - this should not happen with FixParticleSign")
                        ENDIF
                    ENDIF
!Calculate excitation level, connection to HF and diagonal ham element
                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
                    IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
                        HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                    ENDIF
                    IF(ExcitLevel.eq.0) THEN
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                    ELSE
                        HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                        HDiag=(REAL(HDiagTemp%v,r2))-Hii
                    ENDIF

                    do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them
                        NewDets(:,VecSlot)=nJ(:)
                        NewSign(VecSlot)=WSign
                        NewExcits(VecSlot)%ExitGenForDet=.false.
                        NewIC(VecSlot)=ExcitLevel
                        NewH(1,VecSlot)=HDiag      !Diagonal H-element-Hii
                        NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                        VecSlot=VecSlot+1
                    enddo
                                
                    Acceptances=Acceptances+ABS(Child)  !Sum the number of created children to use in acceptance ratio

                ENDIF   !End if child created

            ELSE
!Sum in positive connections to a spin-flip term which modifies the diagonal death-rate
                SpinFlipContrib=SpinFlipContrib+HConnReal

            ENDIF

        enddo   !End of searching through excits

        CALL ResetExit2(CurrentDets(:,Walker),NEl,G1,nBasis,nBasisMax,CurrentExcits(Walker)%ExcitData,0) !Reset excitgen

!Add the spin-flip contribution to the death-rate
        DiagDeath=CurrentH(1,Walker)+SpinFlipContrib

!We now have to decide whether the parent particle (Walker) wants to self-destruct or not...
        iDie=AttemptDie(CurrentDets(:,Walker),DiagDeath)
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births

        IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle

            do l=1,abs(iDie)+1    !We need to copy accross one more, since we need to include the original spared particle
                NewDets(:,VecSlot)=CurrentDets(:,Walker)
                NewSign(VecSlot)=CurrentSign(Walker)
!Copy excitation generator accross
                CALL CopyExitgen(CurrentExcits(Walker),NewExcits(VecSlot))
                NewIC(VecSlot)=CurrentIC(Walker)
                NewH(:,VecSlot)=CurrentH(:,Walker)
                VecSlot=VecSlot+1
            enddo

        ELSEIF(iDie.gt.0) THEN
!This indicates that particles want to be killed. The first kill will simply be performed by not copying accross the original particle.
!Therefore, if iDie = 1, then we can simply ignore it.
!However, after that anti-particles will need to be created on the same determinant.

            do l=1,iDie-1
                NewDets(:,VecSlot)=CurrentDets(:,Walker)
                IF(CurrentSign(Walker)) THEN
!Copy accross new anti-particles
                    NewSign(VecSlot)=.FALSE.
                ELSE
                    NewSign(VecSlot)=.TRUE.
                ENDIF
!Copy excitation generator accross
                CALL CopyExitgen(CurrentExcits(Walker),NewExcits(VecSlot))
                NewIC(VecSlot)=CurrentIC(Walker)
                NewH(:,VecSlot)=CurrentH(:,Walker)
                VecSlot=VecSlot+1
            enddo

        ENDIF   !To kill if

        RETURN

    END SUBROUTINE FixParticleSignIter


    SUBROUTINE ResumGraph(nI,WSign,VecSlot,nIExcitGen,VecInd)
        INTEGER :: nI(NEl),VecSlot,VecInd,ExcitLevel,iGetExcitLevel_2,Create,i,j
        TYPE(ExcitGenerator) :: nIExcitGen
        TYPE(HElement) :: HOffDiag
        LOGICAL :: WSign,ChildSign
        REAL*8 :: Prob,rat,Ran2

        CALL CreateGraph(nI,nIExcitGen,Prob,VecInd)      !Create graph with NDets distinct determinants
        CALL ApplyRhoMat()   !Apply the rho matrix successive times

!First find how many to create at the root determinant
        Create=INT(abs(GraphVec(1)))
        rat=abs(GraphVec(1))-REAL(Create,r2)
        IF(rat.gt.Ran2(Seed)) Create=Create+1
        IF(.not.WSign) Create=-Create
        IF(GraphVec(1).lt.0.D0) Create=-Create
        do j=1,abs(Create)
            NewDets(:,VecSlot)=nI(:)
            IF(Create.lt.0) THEN
                NewSign(VecSlot)=.false.
            ELSE
                NewSign(VecSlot)=.true.
            ENDIF
            NewIC(VecSlot)=CurrentIC(VecInd)
            NewH(:,VecSlot)=CurrentH(:,VecInd)
            CALL CopyExitgen(CurrentExcits(VecInd),NewExcits(VecSlot))
            VecSlot=VecSlot+1
        enddo

        do i=2,NDets
!Now create the new particles according the the final vector GraphVec

            GraphVec(i)=GraphVec(i)/((NDets-1)*Prob)    !Augment the component by the chances of picking that determinant

            Create=INT(abs(GraphVec(i)))
            rat=abs(GraphVec(i))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
            IF(rat.gt.Ran2(Seed)) Create=Create+1
            IF(abs(Create).gt.0) THEN
                IF(.not.WSign) Create=-Create
                IF(GraphVec(i).lt.0.D0) Create=-Create
!Find needed information out about the new particles
!Calculate excitation level, connection to HF. Diagonal ham element info is already stored

                ExcitLevel=iGetExcitLevel_2(HFDet,DetsinGraph(:,i),NEl,NEl)
                IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
                    HOffDiag=GetHElement2(HFDet,DetsinGraph(:,i),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                ENDIF
                IF(Create.lt.0) THEN
                    ChildSign=.false.
                ELSE
                    ChildSign=.true.
                ENDIF

!Now actually create the particles in NewDets and NewSign
                do j=1,abs(Create)
                    NewDets(:,VecSlot)=DetsInGraph(:,i)
                    NewSign(VecSlot)=ChildSign
                    NewIC(VecSlot)=ExcitLevel
                    NewH(1,VecSlot)=GraphKii(i)       !Diagonal H El previously stored
                    NewH(2,VecSlot)=REAL(HOffDiag%v,r2)
                    NewExcits(VecSlot)%ExitGenForDet=.false.
                    VecSlot=VecSlot+1
                enddo

            ENDIF

        enddo

        RETURN

    END SUBROUTINE ResumGraph

    SUBROUTINE CreateGraph(nI,nIExcitGen,Prob,VecInd)
        INTEGER :: nI(NEl),VecInd,nJ(NEl),iCount,IC,i,j,Attempts
        TYPE(ExcitGenerator) :: nIExcitGen
        REAL*8 :: Prob,Kii,ExcitProb
        LOGICAL :: SameDet,CompiPath
        TYPE(HElement) :: Hamij,Hamii

        CALL AZZERO(GraphRhoMat,NDets*NDets)

!Do not need to put the root determinant in the first column of DetsinGraph -
!just assume its there.

        Kii=CurrentH(1,VecInd)      !This is now the Kii element of the root
        GraphRhoMat(1,1)=1.D0-Tau*(Kii-DiagSft)

        i=2
        do while(i.lt.NDets)    !Loop until all determinants found

            CALL GenRandSymExcitIt3(nI,nIExcitGen%ExcitData,nJ,Seed,IC,0,Prob,iCount)

            SameDet=.false.
            do j=2,(i-1)
                IF(CompiPath(nJ,DetsinGraph(:,j),NEl)) THEN
!Determinants are the same as already created determinant - ignore it

                    SameDet=.true.
                    Attempts=Attempts+1
                    IF(Attempts.gt.100) CALL Stop_All("CreateGraphPar","More than 100 attempts needed to grow graph")
                    EXIT
                ENDIF
            enddo

            IF(.not.SameDet) THEN
!Store the unbiased probability of generating excitations from this root - check that it is the same as other excits generated
                IF(i.eq.2) THEN
                    ExcitProb=Prob
                ELSE
                    IF(abs(Prob-ExcitProb).gt.1.D-07) THEN
                        CALL Stop_All("CreateGraph","Excitation probabilities are not uniform - problem here...")
                    ENDIF
                ENDIF

!Determinant is distinct - add it
                DetsinGraph(:,i)=nJ(:)
!First find connection to root
                Hamij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
                GraphRhoMat(1,i)=-Tau*REAL(Hamij%v,r2)
                GraphRhoMat(i,1)=GraphRhoMat(1,i)

!Then find connection to other determinants
                do j=2,(i-1)
                    Hamij=GetHElement2(nJ,DetsInGraph(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
                    GraphRhoMat(i,j)=-Tau*REAL(Hamij%v,r2)
                    GraphRhoMat(j,i)=GraphRhoMat(i,j)
                enddo

!Find diagonal element - and store it for later on...
                Hamii=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                GraphKii(i)=REAL(Hamii%v,r2)-Hii                !Again, the root value is not stored
                GraphRhoMat(i,i)=1.D0-Tau*(GraphKii(i)-DiagSft)

                i=i+1

            ENDIF
        enddo

        RETURN

    END SUBROUTINE CreateGraph

!This applies the rho matrix successive times to a root determinant. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
    SUBROUTINE ApplyRhoMat()
        REAL*8 :: TempVec(NDets)
        INTEGER :: i,j,k

        CALL AZZERO(GraphVec,NDets)
        GraphVec(1)=1.D0        !Set the initial vector to be 1 at the root (i.e. for one walker initially)

        do i=1,RhoApp

            CALL DGEMV('n',NDets,NDets,1.D0,GraphRhoMat,NDets,GraphVec,1,0.D0,TempVec,1)
            CALL DCOPY(NDets,TempVec,1,GraphVec,1)
            CALL AZZERO(TempVec,NDets)

!            do j=1,NDets
!                TempVec(j)=0.D0
!                do k=1,NDets
!                    TempVec(j)=TempVec(j)+GraphRhoMat(j,k)*GraphVec(k)
!                enddo
!            enddo
!            GraphVec(:)=TempVec(:)

        enddo

        RETURN

    END SUBROUTINE ApplyRhoMat

!This will flip the sign of all the particles.
    SUBROUTINE FlipSign()
        INTEGER :: i

        do i=1,TotWalkers
            CurrentSign(i)=.not.CurrentSign(i)
        enddo
        RETURN
    END SUBROUTINE FlipSign


!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!value of the diagonal shift in the hamiltonian in order to compensate for this
    SUBROUTINE UpdateDiagSft()
        IMPLICIT NONE
        INTEGER :: j,k,GrowthSteps

        IF(NoCulls.eq.0) THEN
            IF(TSignShift) THEN
!                WRITE(6,*) TotSign,TotSignOld,TotWalkers,TotWalkersOld
                GrowRate=(ABS(TotSign)+0.D0)/(ABS(TotSignOld)+0.D0)
            ELSE
                GrowRate=(TotWalkers+0.D0)/(TotWalkersOld+0.D0)
            ENDIF
        ELSEIF(NoCulls.eq.1) THEN
!GrowRate is the sum of the individual grow rates for each uninterrupted growth sequence, multiplied by the fraction of the cycle which was spent on it
            IF(TSignShift) THEN
                GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*(ABS(CullInfo(1,1)+0.D0))/(ABS(TotSignOld+0.D0))
                GrowRate=GrowRate+(((StepsSft-CullInfo(1,3))+0.D0)/(StepsSft+0.D0))*((ABS(TotSign)+0.D0)/(ABS(CullInfo(1,2))+0.D0))

            ELSE
                GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((CullInfo(1,1)+0.D0)/(TotWalkersOld+0.D0))
                GrowRate=GrowRate+(((StepsSft-CullInfo(1,3))+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(1,2)+0.D0))
            ENDIF

            NoCulls=0
            CALL IAZZERO(CullInfo,30)
        ELSE
!More than one cull in this update cycle
            IF(TSignShift) THEN
                GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((ABS(CullInfo(1,1))+0.D0)/(ABS(TotSignOld)+0.D0))
                do j=2,NoCulls

!This is needed since the steps between culling is stored cumulatively
                    GrowthSteps=CullInfo(j,3)-CullInfo(j-1,3)
                    GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((ABS(CullInfo(j,1))+0.D0)/(ABS(CullInfo(j-1,2))+0.D0))

                enddo

                GrowthSteps=StepsSft-CullInfo(NoCulls,3)
                GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((ABS(TotSign)+0.D0)/(ABS(CullInfo(NoCulls,2))+0.D0))

            ELSE
                GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((CullInfo(1,1)+0.D0)/(TotWalkersOld+0.D0))
                do j=2,NoCulls

!This is needed since the steps between culling is stored cumulatively
                    GrowthSteps=CullInfo(j,3)-CullInfo(j-1,3)
                    GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((CullInfo(j,1)+0.D0)/(CullInfo(j-1,2)+0.D0))

                enddo

                GrowthSteps=StepsSft-CullInfo(NoCulls,3)
                GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(NoCulls,2)+0.D0))
            ENDIF

            NoCulls=0
            CALL IAZZERO(CullInfo,30)

        ENDIF

        IF(.not.TSinglePartPhase) THEN
            DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
        ENDIF
!        IF((DiagSft).gt.0.D0) THEN
!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!            STOP
!        ENDIF

        MeanExcitLevel=(MeanExcitLevel/(real(SumWalkersCyc,r2)))
        PosFrac=PosFrac/real(SumWalkersCyc,r2)
        ProjectionE=SumENum/(REAL(SumNoatHF,r2))
        AccRat=(REAL(Acceptances,r2))/(REAL(SumWalkersCyc,r2))
        ProjectionMP2=((SumHOverlapMP2/SumOverlapMP2)-Hii)/2.D0

!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers...
        WRITE(15,"(I12,G15.6,I7,G15.6,I10,I6,2G15.6,I11,I9,3G14.6,2I6)") PreviousCycles+Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,Annihilated,ProjectionE,ProjectionMP2,SumNoatHF,NoatDoubs,PosFrac,AccRat,MeanExcitLevel,MinExcitLevel,MaxExcitLevel
        WRITE(6,"(I12,G15.6,I7,G15.6,I10,I6,2G15.6,I11,I9,3G14.6,2I6)") PreviousCycles+Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,Annihilated,ProjectionE,ProjectionMP2,SumNoatHF,NoatDoubs,PosFrac,AccRat,MeanExcitLevel,MinExcitLevel,MaxExcitLevel
!        WRITE(6,*) SumHOverlapMP2,SumOverlapMP2
        CALL FLUSH(15)
        CALL FLUSH(6)

!Now need to reinitialise all variables
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        MeanExcitLevel=0.D0
        SumWalkersCyc=0
        PosFrac=0.D0
        NoatDoubs=0
        Annihilated=0
        Acceptances=0
!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers
        TotSignOld=TotSign

        RETURN

    END SUBROUTINE UpdateDiagSft


!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,Hij0,WSign,j)
        INTEGER :: DetCurr(NEl),ExcitLevel,j
        LOGICAL :: WSign
        REAL*8 :: Hij0      !This is the hamiltonian matrix element between DetCurr and HF

        MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,r2)
        IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
        IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
        IF(ExcitLevel.eq.0) THEN
            IF(WSign) THEN
                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+1
                NoatHF=NoatHF+1
                PosFrac=PosFrac+1.D0
                TotSign=TotSign+1
            ELSE
                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF-1
                NoatHF=NoatHF-1
                TotSign=TotSign-1
            ENDIF
        ELSEIF(ExcitLevel.eq.2) THEN
            NoatDoubs=NoatDoubs+1   !Count number at double excitations in a cycle
!At double excit - sum in energy
            IF(WSign) THEN
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+Hij0
                PosFrac=PosFrac+1.D0
                TotSign=TotSign+1
            ELSE
                IF(Iter.gt.NEquilSteps) SumENum=SumENum-Hij0
                TotSign=TotSign-1
            ENDIF
        ELSE
            IF(WSign) THEN
                PosFrac=PosFrac+1.D0
                TotSign=TotSign+1
            ELSE
                TotSign=TotSign-1
            ENDIF
        ENDIF

        IF(TProjEMP2.and.(Iter.gt.NEquilSteps)) THEN
!Calculate the overlap of the wavefunction with the MP2 wavefunction in order to calculate the energy.
            CALL AddProjMP2Energy(DetCurr,ExcitLevel,Hij0,WSign,j)
        ENDIF

        RETURN

    END SUBROUTINE SumEContrib

!This routine will add the energy contribution from the overlap of the walker with the MP2 wavefunction
!We want to calculate the energy as: sum_n <Psi_MP2 | H | Psi>/ sum_n <Psi_MP2 | Psi> where n is all iterations
    SUBROUTINE AddProjMP2Energy(DetCurr,ExcitLevel,Hij0,WSign,j)
        INTEGER :: DetCurr(NEl),ExcitLevel,j,OrbI,OrbJ,OrbA,OrbB,t
        INTEGER :: ExcitForm(2,2),nJ(NEl),nStore(6),iExcit,nJ_IC
        INTEGER :: iGetExcitLevel_2
        REAL*8 :: Hij0,MP1Compt,Hij,Hjj
        LOGICAL :: WSign,tSign
        TYPE(HElement) :: Hijtemp,TempFjj,Hij2temp

!First calculate the overlap of the wavefunction with the MP2 wavefunction - only doubles/HF will contribute
        IF(ExcitLevel.eq.2) THEN

            TempFjj=GetH0Element3(DetCurr)
            MP1Compt=-Hij0/(Fii-(REAL(TempFjj%v,r2)))
            IF(WSign) THEN
                SumOverlapMP2=SumOverlapMP2+MP1Compt
            ELSE
                SumOverlapMP2=SumOverlapMP2-MP1Compt
            ENDIF
            
!!We need to find the ijab orbitals of the double excitation
!            ExcitForm(1,1)=ExcitLevel
!            
!            CALL GetExcitation(HFDet,DetCurr,NEl,ExcitForm,tSign)
!!The i,j orbitals are ExcitForm(1,1) and (1,2). a/b are (2,1) and (2,2)
!!Find the ordering of the orbitals in terms of energy
!            OrbI=INVBRRSpinOrb(ExcitForm(1,1))
!            OrbJ=INVBRRSpinOrb(ExcitForm(1,2))
!            OrbA=INVBRRSpinOrb(ExcitForm(2,1))
!            OrbB=INVBRRSpinOrb(ExcitForm(2,2))
!
!!Ensure i > j and a > b.
!            IF(OrbI.lt.OrbJ) THEN
!                t=OrbI
!                OrbI=OrbJ
!                OrbJ=t
!            ELSEIF(OrbI.eq.OrbJ) THEN
!                CALL Stop_All("AddProjMP2Energy","Cannot have I.eq.J")
!            ENDIF
!            IF(OrbA.lt.OrbB) THEN
!                t=OrbA
!                OrbA=OrbB
!                OrbB=t
!            ELSEIF(OrbA.eq.OrbB) THEN
!                CALL Stop_All("AddProjMP2Energy","Cannot have A.eq.B")
!            ENDIF
!            IF(OrbI.ge.OrbB) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!            IF((OrbI.gt.NEl).or.(OrbJ.gt.NEl)) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!            IF((OrbA.gt.nBasis).or.(OrbB.gt.nBasis)) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!            IF((OrbA.le.NEl).or.(OrbB.le.NEl)) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!            
!            OrbA=OrbA-NEl     !Virtual orbital indexing cannot be > NEl, so remove to save space.
!            OrbB=OrbB-NEl
!            
!            SumOverlapMP2=SumOverlapMP2+MP2ExcitComps(OrbI,OrbJ,OrbA,OrbB)    !Add component of MP1 Wavefunction

            IF(WSign) THEN
                Hjj=CurrentH(1,j)
            ELSE
                Hjj=-CurrentH(1,j)
            ENDIF
            SumHOverlapMP2=SumHOverlapMP2+(MP1Compt*Hjj)

        ELSEIF(ExcitLevel.eq.0) THEN
!HF MP2 wavefunction component is simply 1 since wavefunction is unnormalised.

            IF(WSign) THEN
                SumOverlapMP2=SumOverlapMP2+1.D0
                SumHOverlapMP2=SumHOverlapMP2+Hii
            ELSE
                SumOverlapMP2=SumOverlapMP2-1.D0
                SumHOverlapMP2=SumHOverlapMP2-Hii
            ENDIF

        ENDIF

!Now, for the more difficult task of finding the overlap of H|Psi> with the MP1 wavefunction
        IF(ExcitLevel.le.4) THEN
!Any walker which is at a quadruple or less excitation level will have excitations which overlap with MP1.
!Setup excit generators for this determinant if needed
            
            CALL SetupExitgen(DetCurr,CurrentExcits(j))
            CALL ResetExit2(DetCurr,NEl,G1,nBasis,nBasisMax,CurrentExcits(j)%ExcitData,0) !Reset excitgen
            
            do while(.true.)
!Run through all excitations - nStore has nothing in it (hopefully not needed)
                CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.false.,CurrentExcits(j)%ExcitData,nJ,iExcit,0,nStore,3)
                IF(nJ(1).eq.0) exit
                nJ_IC=iGetExcitLevel_2(HFDet,nJ,NEl,2)  !Find out if excitation is a double

                IF(nJ_IC.eq.2) THEN
!Excitation is less than or equal to a double - find connection to original walker
                    Hijtemp=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,iExcit,ECore)
                    IF(WSign) THEN
                        Hij=REAL(Hijtemp%v,r2)
                    ELSE
                        Hij=-REAL(Hijtemp%v,r2)
                    ENDIF
!Find MP1 excitation compt
                    TempFjj=GetH0Element3(nJ)
                    Hij2temp=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,nJ_IC,ECore)
                    MP1Compt=-(REAL(Hij2temp%v,r2))/(Fii-(REAL(TempFjj%v,r2)))
                    SumHOverlapMP2=SumHOverlapMP2+(MP1Compt*Hij)

!...and then find overlap of this excitation with the MP1 wavefunction
!                    ExcitForm(1,1)=nJ_IC
!                    CALL GetExcitation(HFDet,nJ,NEl,ExcitForm,tSign)
!!The i,j orbitals are ExcitForm(1,1) and (1,2). a/b are (2,1) and (2,2)
!!Find the ordering of the orbitals in terms of energy
!                    OrbI=INVBRRSpinOrb(ExcitForm(1,1))
!                    OrbJ=INVBRRSpinOrb(ExcitForm(1,2))
!                    OrbA=INVBRRSpinOrb(ExcitForm(2,1))
!                    OrbB=INVBRRSpinOrb(ExcitForm(2,2))
!
!!Ensure i > j and a > b.
!                    IF(OrbI.lt.OrbJ) THEN
!                        t=OrbI
!                        OrbI=OrbJ
!                        OrbJ=t
!                    ELSEIF(OrbI.eq.OrbJ) THEN
!                        CALL Stop_All("AddProjMP2Energy","Cannot have I.eq.J")
!                    ENDIF
!                    IF(OrbA.lt.OrbB) THEN
!                        t=OrbA
!                        OrbA=OrbB
!                        OrbB=t
!                    ELSEIF(OrbA.eq.OrbB) THEN
!                        CALL Stop_All("AddProjMP2Energy","Cannot have A.eq.B")
!                    ENDIF
!                    IF(OrbI.ge.OrbB) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!                    IF((OrbI.gt.NEl).or.(OrbJ.gt.NEl)) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!                    IF((OrbA.gt.nBasis).or.(OrbB.gt.nBasis)) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!                    IF((OrbA.le.NEl).or.(OrbB.le.NEl)) CALL Stop_All("AddProjMP2Energy","Indexing scheme incorrect")
!                    
!                    OrbA=OrbA-NEl     !Virtual orbital indexing cannot be > NEl, so remove to save space.
!                    OrbB=OrbB-NEl
!                    
!!Add <MP1|H|Psi>
!                    SumHOverlapMP2=SumHOverlapMP2+(MP2ExcitComps(OrbI,OrbJ,OrbA,OrbB)*(real(Hijtemp%v,r2)))

                ELSEIF(nJ_IC.eq.0) THEN
!Excitation is the HF determinant - find connection to it

                    Hijtemp=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,iExcit,ECore)
                    IF(WSign) THEN
                        Hij=REAL(Hijtemp%v,r2)
                    ELSE
                        Hij=-REAL(Hijtemp%v,r2)
                    ENDIF

                    SumHOverlapMP2=SumHOverlapMP2+Hij

                ENDIF

            enddo

            CALL ResetExit2(DetCurr,NEl,G1,nBasis,nBasisMax,CurrentExcits(j)%ExcitData,0) !Reset excitgen

        ENDIF

    END SUBROUTINE AddProjMP2Energy


!This routine acts as a thermostat for the simulation - killing random particles if the population becomes too large, or 
!Doubling them if it gets too low...
    SUBROUTINE ThermostatParticles(HighLow)
        IMPLICIT NONE
        LOGICAL :: HighLow
        INTEGER :: VecSlot,i,j,ToCull,Culled,OrigWalkers,Chosen
        REAL*8 :: Ran2

        IF(HighLow) THEN
!The population is too large - cull TotWalkers/CullFactor randomly selected particles

            OrigWalkers=TotWalkers
            ToCull=TotWalkers-nint((TotWalkers+0.D0)/CullFactor)
            Culled=0

            do while (Culled.lt.ToCull)

!Pick a random walker between 1 and TotWalkers
                Chosen=int((Ran2(Seed)*TotWalkers)+1.D0)

                IF(CurrentSign(Chosen)) THEN
!Walker that is removed is positive - therefore totsign is one less than it was before
                    TotSign=TotSign-1
                ELSE
                    TotSign=TotSign+1
                ENDIF

!Move the Walker at the end of the list to the position of the walker we have chosen to destroy
                CurrentDets(:,Chosen)=CurrentDets(:,TotWalkers)
                CurrentSign(Chosen)=CurrentSign(TotWalkers)
                CurrentH(:,Chosen)=CurrentH(:,TotWalkers)
                CurrentIC(Chosen)=CurrentIC(TotWalkers)
                CALL CopyExitgen(CurrentExcits(TotWalkers),CurrentExcits(Chosen))
                CurrentExcits(TotWalkers)%ExitGenForDet=.false.

                TotWalkers=TotWalkers-1
                Culled=Culled+1

            enddo

            IF(TotWalkers.ne.(OrigWalkers-ToCull)) THEN
                WRITE(6,*) "Error in culling walkers..."
                STOP "Error in culling walkers..."
            ENDIF

            IF(TSignShift) THEN
!CullInfo(:,2) is the new number of total walkers/residualsign
                CullInfo(NoCulls,2)=TotSign
            ELSE
                CullInfo(NoCulls,2)=TotWalkers
            ENDIF

        ELSE
!The population is too low - give it a boost by doubling every particle

            VecSlot=TotWalkers+1
            do i=1,TotWalkers
                
                IF(CurrentSign(Chosen)) THEN
!Walker that is doubled is positive - therefore totsign is one more than it was before
                    TotSign=TotSign+1
                ELSE
                    TotSign=TotSign-1
                ENDIF

!Add clone of walker, at the same determinant, to the end of the list
                CurrentDets(:,VecSlot)=CurrentDets(:,i)
                CurrentSign(VecSlot)=CurrentSign(i)
                CurrentH(:,VecSlot)=CurrentH(:,i)
                CurrentIC(VecSlot)=CurrentIC(i)
                CALL CopyExitgen(CurrentExcits(i),CurrentExcits(VecSlot))

                VecSlot=VecSlot+1

            enddo

            TotWalkers=TotWalkers*2

            IF((VecSlot-1).ne.TotWalkers) THEN
                WRITE(6,*) "Problem in doubling all particles..."
                STOP "Problem in doubling all particles..."
            ENDIF

            IF(TSignShift) THEN
!CullInfo(:,2) is the new number of total walkers/residualsign
                CullInfo(NoCulls,2)=TotSign
            ELSE
                CullInfo(NoCulls,2)=TotWalkers
            ENDIF

        ENDIF

        RETURN

    END SUBROUTINE ThermostatParticles


!This routine cancels out particles of opposing sign on the same determinant.
    SUBROUTINE AnnihilatePairs(TotWalkersNew)
        TYPE(ExcitGenerator) :: TempExcit
        REAL*8 :: TempH(2)
        INTEGER :: TempIC
        INTEGER :: TotWalkersNew,j,k,l,DetCurr(NEl),VecSlot,TotWalkersDet
        INTEGER :: DetLT

!First, it is necessary to sort the list of determinants
        CALL SortParts(TotWalkersNew,NewDets(:,1:TotWalkersNew),NEl)

!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to CurrentDets
        j=1
!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
!DetCurr is the current determinant
        CALL ICOPY(NEl,NewDets(:,j),1,DetCurr(:),1)
        TempIC=NewIC(j)
        TempH(:)=NewH(:,j)
        CALL CopyExitGen(NewExcits(j),TempExcit)
        
        VecSlot=1

        do while(j.le.TotWalkersNew)
!Loop over all walkers
            TotWalkersDet=0
            do while ((DetLT(NewDets(:,j),DetCurr,NEl).eq.0).and.(j.le.TotWalkersNew))
!Loop over all walkers on DetCurr and count residual number after cancelling
                IF(NewSign(j)) THEN
                    TotWalkersDet=TotWalkersDet+1
                ELSE
                    TotWalkersDet=TotWalkersDet-1
                ENDIF
                j=j+1
            enddo
!Transfer residual population into VecSlot, along with residual sign
            IF(TotWalkersDet.gt.0) THEN
!Positive sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    CALL ICOPY(NEl,DetCurr(:),1,CurrentDets(:,VecSlot),1)
                    CurrentSign(VecSlot)=.true.
                    CurrentIC(VecSlot)=TempIC
                    CurrentH(:,VecSlot)=TempH(:)
                    CALL CopyExitgen(TempExcit,CurrentExcits(VecSlot))
                    VecSlot=VecSlot+1
                enddo
            ELSE
!Negative sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    CALL ICOPY(NEl,DetCurr(:),1,CurrentDets(:,VecSlot),1)
                    CurrentSign(VecSlot)=.false.
                    CurrentIC(VecSlot)=TempIC
                    CurrentH(:,VecSlot)=TempH(:)
                    CALL CopyExitgen(TempExcit,CurrentExcits(VecSlot))
                    VecSlot=VecSlot+1
                enddo
            ENDIF
!Now update the current determinant
            CALL ICOPY(NEl,NewDets(:,j),1,DetCurr(:),1)
            TempIC=NewIC(j)
            TempH(:)=NewH(:,j)
            CALL CopyExitGen(NewExcits(j),TempExcit)
        enddo
!The new number of residual cancelled walkers is given by one less that VecSlot again.
        TotWalkers=VecSlot-1

        RETURN

    END SUBROUTINE AnnihilatePairs

!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalc()
        USE Calc, only : EXCITFUNCS
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot,error,HFConn
        REAL*8 :: Ran2
        TYPE(HElement) :: rh,TempHii
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMC'

        IF(HElementSize.gt.1) THEN
            CALL Stop_All("FCIMCPar","FciMCPar cannot function with complex orbitals.")
        ENDIF

!Open a file to store output
        OPEN(15,file='FCIMCStats',status='unknown')

!Store information specifically for the HF determinant
        ALLOCATE(HFDet(NEl),stat=ierr)
        CALL LogMemAlloc('HFDet',NEl,4,this_routine,HFDetTag)
        do i=1,NEl
            HFDet(i)=FDet(i)
        enddo

!Setup excitation generator for the HF determinant
        CALL SetupExitgen(HFDet,HFExcit)
        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=G_VMC_Seed

!Calculate Hii
        TempHii=GetHElement2(HFDet,HFDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
        Hii=REAL(TempHii%v,r2)          !Diagonal Hamiltonian element for the HF determinant

        TempHii=GetH0Element3(HFDet)
        Fii=REAL(TempHii%v,r2)          !Fock-energy of the HF determinant

!Initialise variables for calculation on each node
        ProjectionE=0.D0
        ProjectionMP2=0.D0
        PosFrac=0.D0
        SumENum=0.D0
        SumNoatHF=0
        SumOverlapMP2=0.D0
        SumHOverlapMP2=0.D0
        MeanExcitLevel=0.D0
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        NoatDoubs=0
        Acceptances=0
        Annihilated=0
        AccRat=0.D0
        PreviousCycles=0
        IF(TStartSinglePart) THEN
            TSinglePartPhase=.true.
            IF(TReadPops) THEN
                CALL Stop_All("InitFciMCCalc","Cannot read in POPSFILE as well as starting with a single particle")
            ENDIF
        ELSE
            TSinglePartPhase=.false.
        ENDIF

        IF(TResumFciMC) THEN
            IF(NDets.gt.2) THEN
                IF(.not.EXCITFUNCS(10)) THEN
                    WRITE(6,*) "Cannot have an excitation bias with multiple determinant graphs...exiting."
                    CALL Stop_All("InitFCIMCCalc","Cannot have biasing with Graphsizes > 2")
                ENDIF

!Allocate memory for graphs...
                ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
                CALL LogMemAlloc('GraphRhoMat',NDets**2,8,this_routine,GraphRhoMatTag,ierr)
                ALLOCATE(GraphVec(NDets),stat=ierr)
                CALL LogMemAlloc('GraphVec',NDets,8,this_routine,GraphVecTag,ierr)
                ALLOCATE(GraphKii(NDets),stat=ierr)
                CALL LogMemAlloc('GraphKii',NDets,8,this_routine,GraphKiiTag,ierr)
                ALLOCATE(DetsinGraph(NEl,NDets),stat=ierr)
                CALL LogMemAlloc('DetsinGraph',NDets*NEl,4,this_routine,DetsinGraphTag,ierr)

            ELSEIF(NDets.lt.2) THEN
                WRITE(6,*) "Graphs cannot be smaller than two vertices. Exiting."
                CALL Stop_All("InitFCIMCCalc","Graphs cannot be smaller than two vertices")
            ENDIF
            WRITE(6,*) "Resumming in multiple transitions to/from each excitation"
            WRITE(6,"(A,I5,A)") "Graphs to resum will consist of ",NDets," determinants."
        ENDIF
        WRITE(6,*) ""
        WRITE(6,*) "Performing FCIMC...."
        WRITE(6,*) "Maximum connectivity of HF determinant is: ",HFConn

        IF(ICILevel.ne.0) THEN
            IF(TResumFCIMC.or.TFixParticleSign.or.THFRetBias) THEN
                CALL Stop_All("InitFCIMCCalc","ResumFCIMC, FixParticleSign and HFRetBias FCIMC methods do not currently work in a truncated space.")
            ELSE
                WRITE(6,*) "Excitation levels w.r.t. HF restricted to level: ",ICILevel
            ENDIF
        ENDIF

        IF(TStartMP1) THEN
!Start the initial distribution off at the distribution of the MP1 eigenvector

            WRITE(6,*) "This feature not ready yet"
            CALL FLUSH(6)
            CALL Stop_All("InitFCIMC","StartMP1 feature no longer operational")
            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."

        ELSEIF(TReadPops) THEN
!We are reading in an initial distribution of walkers

            WRITE(6,*) "Reading in initial particle configuration from POPSFILE..." 
            CALL ReadFromPopsFile()

        ELSE

!initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkers=MemoryFac*InitWalkers
        
            IF(TStartSinglePart) THEN
                WRITE(6,"(A,I9)") "Initial number of particles set to 1, and shift will be held at 0.D0 until particle number gets to ",InitWalkers
            ELSE
                WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
            ENDIF
            WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft

!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
            CALL IAZZERO(WalkVecDets,NEl*MaxWalkers)
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            CALL IAZZERO(WalkVec2Dets,NEl*MaxWalkers)
            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

            ALLOCATE(WalkVecIC(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecIC',MaxWalkers,4,this_routine,WalkVecICTag,ierr)
            ALLOCATE(WalkVec2IC(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2IC',MaxWalkers,4,this_routine,WalkVec2ICTag,ierr)
            ALLOCATE(WalkVecH(2,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkers*2,8,this_routine,WalkVecHTag,ierr)
            CALL AZZERO(WalkVecH,2*MaxWalkers)
            ALLOCATE(WalkVec2H(2,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2H',MaxWalkers*2,8,this_routine,WalkVec2HTag,ierr)
            CALL AZZERO(WalkVec2H,2*MaxWalkers)

!Allocate pointers to the correct walker arrays
            CurrentDets=>WalkVecDets
            CurrentSign=>WalkVecSign
            CurrentIC=>WalkVecIC
            CurrentH=>WalkVecH
            NewDets=>WalkVec2Dets
            NewSign=>WalkVec2Sign
            NewIC=>WalkVec2IC
            NewH=>WalkVec2H

            IF(TStartSinglePart) THEN
                CurrentDets(:,1)=HFDet(:)
                CurrentIC(1)=0
                CurrentSign(1)=.true.
            ELSE
                do j=1,InitWalkers
                    CurrentDets(:,j)=HFDet(:)
                    CurrentIC(j)=0
                    IF(TFixParticleSign) THEN
                        IF(Ran2(Seed).lt.0.40) THEN
                            CurrentSign(j)=.false.
                        ELSE
                            CurrentSign(j)=.true.
                        ENDIF
                    ELSE
                        CurrentSign(j)=.true.
                    ENDIF
                enddo
            ENDIF

            WRITE(6,*) "Initial memory allocation sucessful..."

            ALLOCATE(WalkVecExcits(MaxWalkers),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("InitFCIMCCalc","Error in allocating walker excitation generators")
            ALLOCATE(WalkVec2Excits(MaxWalkers),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("InitFCIMCCalc","Error in allocating walker excitation generators")

!Allocate pointers to the correct excitation arrays
            CurrentExcits=>WalkVecExcits
            NewExcits=>WalkVec2Excits

            IF(TStartSinglePart) THEN
                CALL CopyExitGen(HFExcit,CurrentExcits(1))
            ELSE
                do j=1,InitWalkers
    !Copy the HF excitation generator accross to each initial particle
                    CALL CopyExitGen(HFExcit,CurrentExcits(j))
                enddo
            ENDIF

            WRITE(6,*) "Initial allocation of excitation generators successful..."
            CALL FLUSH(6)

        ENDIF

        IF(THFRetBias) THEN
            IF((PRet.gt.1.D0).or.(PRet.lt.0.D0)) THEN
                CALL Stop_All("InitFciMCCalc","HFRetBias set, but return bias is not a normalised probability")
            ELSE
                WRITE(6,*) "Return bias to HF determinant set, with PRet = ", PRet
            ENDIF
            IF(TExcludeRandGuide.and.(.not.EXCITFUNCS(10))) THEN
                CALL Stop_All("InitFCIMCCalc","Cannot have excitation weighting if using ExcludeRandGuide unbiasing with a guiding function")
            ENDIF
            IF(TExcludeRandGuide.and.(PRet.lt.0.D0)) THEN
                CALL Stop_All("InitFCIMCCalc","Cannont have PRet=0 with EXCLUDERANDGUIDE unbiasing as cannot return to HF")
            ENDIF
        ENDIF

        IF(TStartSinglePart) THEN
            TotWalkers=1
            TotWalkersOld=1
            TotSign=1
            TotSignOld=1
        ELSE
!TotWalkers contains the number of current walkers at each step
            TotWalkers=InitWalkers
            TotSign=InitWalkers
            TotWalkersOld=InitWalkers
            TotSignOld=InitWalkers
        ENDIF
!TotSign is the sum of the walkers x sign
        ProjectionE=SumENum/(REAL(SumNoatHF,r2))

        CALL IAZZERO(CullInfo,30)
        NoCulls=0

!Routine to initialise the blocking analysis
        CALL InitBlocking()

        IF(TProjEMP2) THEN
!This will perform an MP2 calculation, so that the projection to it can be used to calculate the energy
            CALL InitMP2Calc()
        ENDIF

        RETURN

    END SUBROUTINE InitFCIMCCalc

    SUBROUTINE InitMP2Calc()
        INTEGER :: ierr,i,j,nJ(NEl),IC,nStore(6),ExcitForm(2,2),iExcit
        INTEGER :: A,B,t,MaxComptDet(NEl),Compts
        REAL*8 :: Denom,MaxCompt,MP2Energy,MP1Comp
        LOGICAL :: tSign
        TYPE(HElement) :: Hij
        CHARACTER(Len=*) , PARAMETER :: this_routine='InitMP2Calc'
        
!The MP1 wavefunction components are no longer precalculated, but rather calculated on the fly now.
!Need to calculate all the components of the MP2 wavefunction, and store them in MP2ExcitComps
!MP2 ExcitComps stores the eigenvector for ij->ab, with each dimension representing an index
!i > j and a > b. An indexing scheme and a 1D array could be used.
!The HF component is 1, since we are working with an unnormalised eigenvector.
!        ALLOCATE(MP2ExcitComps(NEl,NEl,nBasis-NEl,nBasis-NEl),stat=ierr)
!        CALL LogMemAlloc("MP2ExcitComps",NEl*NEl*(nBasis-NEl)**2,8,this_routine,MP2ExcitCompsTag,ierr)
!        CALL AZZERO(MP2ExcitComps,NEl*NEl*(nBasis-NEl)**2)
!
!!To store in excitation form - need to be able to order orbitals purely in terms of energy, so setup INVBRR
!!This is different to the INVBRR in UMatCache, since it works with spin-orbitals
!        ALLOCATE(INVBRRSpinOrb(nBasis),stat=ierr)
!        CALL LogMemAlloc("InvBrrSpinOrb",nBasis,4,this_routine,InvBrrSpinOrbTag,ierr)
!        CALL IAZZERO(INVBRRSpinOrb,nBasis)
!        t=0
!        do i=1,nBasis
!            t=t+1
!            INVBRRSpinOrb(BRR(i))=t
!        enddo

!        WRITE(6,*) "INVBRRSpinOrb is: ",INVBRRSpinOrb(:)

        MP2Energy=Hii  !From HF Determinant which has a wavevector component of 1
        MaxCompt=1.D0
        MaxComptDet(:)=HFDet(:)
        Compts=0

!Reset the HF Excitation generator
        CALL ResetExIt2(HFDet,NEl,G1,nBasis,nBasisMax,HFExcit%ExcitData,0)

        do while(.true.)
!Generate double excitations
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.false.,HFExcit%Excitdata,nJ,iExcit,0,nStore,2)
            IF(nJ(1).eq.0) EXIT
            Compts=Compts+1     !Calculate total number of MP1 excitations
!            ExcitForm(1,1)=2    !signify that we are only dealing with double excitations
            IF(iExcit.ne.2) CALL Stop_All("InitMP2Calc","Excitations other than double being generated")
!            CALL GetExcitation(HFDet,nJ,NEl,ExcitForm,tSign)
!!The i,j orbitals are ExcitForm(1,1) and (1,2). a/b are (2,1) and (2,2)
!!Find the ordering of the orbitals in terms of energy
!            I=INVBRRSpinOrb(ExcitForm(1,1))
!            J=INVBRRSpinOrb(ExcitForm(1,2))
!            A=INVBRRSpinOrb(ExcitForm(2,1))
!            B=INVBRRSpinOrb(ExcitForm(2,2))
!
!!Ensure i > j and a > b.
!            IF(I.lt.J) THEN
!                t=I
!                I=J
!                J=t
!            ELSEIF(I.eq.J) THEN
!                CALL Stop_All("InitMP2Calc","Cannot have I.eq.J")
!            ENDIF
!            IF(A.lt.B) THEN
!                t=A
!                A=B
!                B=t
!            ELSEIF(A.eq.B) THEN
!                CALL Stop_All("InitMP2Calc","Cannot have A.eq.B")
!            ENDIF
!            IF(I.ge.B) CALL Stop_All("InitMP2Calc","Indexing scheme incorrect")
!            IF((I.gt.NEl).or.(J.gt.NEl)) CALL Stop_All("InitMP2Calc","Indexing scheme incorrect")
!            IF((A.gt.nBasis).or.(B.gt.nBasis)) CALL Stop_All("InitMP2Calc","Indexing scheme incorrect")
!            IF((A.le.NEl).or.(B.le.NEl)) CALL Stop_All("InitMP2Calc","Indexing scheme incorrect")
!            
!            A=A-NEl     !Virtual orbital indexing cannot be > NEl, so remove to save space.
!            B=B-NEl

!Now put component of MP2 wavefunction into MP2ExcitComps(I,J,A,B)
            Denom=Arr(ExcitForm(2,1),2)-Arr(ExcitForm(1,1),2)+Arr(ExcitForm(2,2),2)-Arr(ExcitForm(1,2),2)
            Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,2,ECore)

!            MP2ExcitComps(I,J,A,B)=-REAL(Hij%v,r2)/Denom    !Store MP2 Wavefunction
            MP1Comp=-REAL(Hij%v,r2)/Denom
            IF((ABS(MP1Comp)).gt.MaxCompt) THEN 
!Find the maximum component of the MP2 wavefunction
                MaxCompt=ABS(MP1Comp)
                MaxComptDet(:)=nJ(:)
            ENDIF
            MP2Energy=MP2Energy-(REAL(Hij%v,r2)**2)/Denom   !Calculate MP2 Energy

        enddo

!Reset the HF Excitation generator
        CALL ResetExIt2(HFDet,NEl,G1,nBasis,nBasisMax,HFExcit%ExcitData,0)

        WRITE(6,"(I7,A)") Compts," MP2 components calculated. Maximum MP2 wavevector component is determinant: "
        CALL WRITEDET(6,MaxComptDet,NEl,.TRUE.)
        WRITE(6,*) "MP2 ENERGY = ",MP2Energy
        CALL FLUSH(6)

        RETURN

    END SUBROUTINE InitMP2Calc

    SUBROUTINE InitBlocking()
        INTEGER :: TotSamples,TotBlocks,i

!First, we want to find the theoretical maximum number of different block sizes there could be in the simulation.
!This is the total number of shiftchange samples theoretically possible if the simulation goes to completion.
        TotSamples=INT(real(NMCyc)/real(StepsSft)) 

!I'm sure theres a more elegant way to do this, but this will calculate the maximum number of different block sizes possible
!This is calculated by seeing if the run will give at least two complete blocks of data.
        TotBlocks=0
        do i=0,100
            IF((TotSamples/(2**i)).ge.2) THEN
                TotBlocks=TotBlocks+1
            ELSE
                EXIT
            ENDIF
        enddo

        WRITE(6,"(A,I6)") "Total number of different sizes blocks possible is: ",TotBlocks

!We will have Totblocks different block sizes, of length: 1,2,4,8,...,2**TotBlocks.
!Allocate memory to hold the stats for these different block lengths.


        RETURN
    END SUBROUTINE InitBlocking

    SUBROUTINE WriteToPopsFile()
        INTEGER :: j,k

        WRITE(6,*) "Writing to POPSFILE..."
        OPEN(17,FILE='POPSFILE',Status='unknown')
        WRITE(17,*) TotWalkers, "   TOTWALKERS"
        WRITE(17,*) DiagSft, "   DIAGSHIFT"
        WRITE(17,*) SumNoatHF, "   SUMNOATHF"
        WRITE(17,*) SumENum,"   SUMENUM ( \sum<D0|H|Psi> )"
        WRITE(17,*) Iter+PreviousCycles, "   PREVIOUS CYCLES"
        do j=1,TotWalkers
            do k=1,NEl
                WRITE(17,"(I5)",advance='no') CurrentDets(k,j)
            enddo
            WRITE(17,*) CurrentSign(j)
        enddo
        CLOSE(17)

        RETURN

    END SUBROUTINE WriteToPopsFile
    
    SUBROUTINE ReadFromPopsFile()
        INTEGER :: ierr,l,j,k,VecSlot,IntegerPart,iGetExcitLevel_2
        REAL*8 :: FracPart,Ran2
        TYPE(HElement) :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsFile'
        LOGICAL :: exists

        INQUIRE(FILE='POPSFILE',EXIST=exists)
        IF(.not.exists) THEN
            CALL Stop_All("ReadFromPopsFile","POPSFILE not present - cannot read in particle configuration")
        ENDIF
        OPEN(17,FILE='POPSFILE',Status='old')
!Read in initial data
        READ(17,*) InitWalkers
        READ(17,*) DiagSft
        READ(17,*) SumNoatHF
        READ(17,*) SumENum
        READ(17,*) PreviousCycles

        WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
        IF(NEquilSteps.gt.0) THEN
            WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
            NEquilSteps=0
        ENDIF

        MaxWalkers=MemoryFac*InitWalkers
!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
        CALL IAZZERO(WalkVecDets,NEl*MaxWalkers)
        ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)

        do l=1,InitWalkers
            READ(17,*) WalkVecDets(1:NEl,l),WalkVecSign(l)
        enddo

        WRITE(6,*) InitWalkers," particles read in from POPSFILE..."

        IF(ScaleWalkers.ne.1.D0) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers
            MaxWalkers=MemoryFac*(nint(InitWalkers*ScaleWalkers))
            
!Allocate more memory for WalkVec2
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            CALL IAZZERO(WalkVec2Dets,NEl*MaxWalkers)
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            VecSlot=1
            do l=1,InitWalkers
                do k=1,IntegerPart
                    WalkVec2Dets(1:NEl,VecSlot)=WalkVecDets(1:NEl,l)
                    WalkVec2Sign(VecSlot)=WalkVecSign(l)
                    VecSlot=VecSlot+1
                enddo
                IF(Ran2(Seed).lt.FracPart) THEN
!Create extra stochastically created particle
                    WalkVec2Dets(1:NEl,VecSlot)=WalkVecDets(1:NEl,l)
                    WalkVec2Sign(VecSlot)=WalkVecSign(l)
                    VecSlot=VecSlot+1
                ENDIF
            enddo

            WRITE(6,*) "Total number of initial walkers is now: ",VecSlot-1
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is now: ", DiagSft
            InitWalkers=nint(InitWalkers*ScaleWalkers)

!Deallocate the arrays used to hold the original walkers, and reallocate with right size
            DEALLOCATE(WalkVecDets)
            CALL LogMemDealloc(this_routine,WalkVecDetsTag)
            DEALLOCATE(WalkVecSign)
            CALL LogMemDealloc(this_routine,WalkVecSignTag)
            ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
            CALL IAZZERO(WalkVecDets,NEl*MaxWalkers)
            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)

!Transfer initial data accross to WalkVecDets
            do l=1,InitWalkers
                WalkVecDets(1:NEl,l)=WalkVec2Dets(1:NEl,l)
                WalkVecSign(l)=WalkVec2Sign(l)
            enddo

!Zero the second arrays
            CALL IAZZERO(WalkVec2Dets,NEl*MaxWalkers)

        ELSE
!If not scaling, second array has not been allocated, allocate it now
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            CALL IAZZERO(WalkVec2Dets,NEl*MaxWalkers)
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)
            
            WRITE(6,*) "Total number of initial walkers is now: ",InitWalkers
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is now: ", DiagSft
        
        ENDIF

        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp

!Now need to allocate other arrays
        ALLOCATE(WalkVecIC(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecIC',MaxWalkers,4,this_routine,WalkVecICTag,ierr)
        ALLOCATE(WalkVec2IC(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2IC',MaxWalkers,4,this_routine,WalkVec2ICTag,ierr)
        ALLOCATE(WalkVecH(2,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecH',MaxWalkers*2,8,this_routine,WalkVecHTag,ierr)
        CALL AZZERO(WalkVecH,2*MaxWalkers)
        ALLOCATE(WalkVec2H(2,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2H',MaxWalkers*2,8,this_routine,WalkVec2HTag,ierr)
        CALL AZZERO(WalkVec2H,2*MaxWalkers)

!Allocate pointers to the correct walker arrays
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign
        CurrentIC=>WalkVecIC
        CurrentH=>WalkVecH
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign
        NewIC=>WalkVec2IC
        NewH=>WalkVec2H
            
        ALLOCATE(WalkVecExcits(MaxWalkers),stat=ierr)
        ALLOCATE(WalkVec2Excits(MaxWalkers),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("InitFCIMCCalc","Error in allocating walker excitation generators")

!Allocate pointers to the correct excitation arrays
        CurrentExcits=>WalkVecExcits
        NewExcits=>WalkVec2Excits
            
!Initialise data for new walkers
        do j=1,InitWalkers
            CurrentIC(j)=iGetExcitLevel_2(HFDet,CurrentDets(:,j),NEl,NEl)
            IF(CurrentIC(j).eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
                HElemTemp=GetHElement2(HFDet,CurrentDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,CurrentIC(j),ECore)
                CurrentH(2,j)=REAL(HElemTemp%v,r2)
            ELSEIF(CurrentIC(j).eq.0) THEN
!We know we are at HF - HDiag=0, and can use HF excitgen
                CALL CopyExitGen(HFExcit,CurrentExcits(j))
                CurrentH(1,j)=0.D0
                CurrentH(2,j)=0.D0
                
            ELSE
                CurrentExcits(j)%ExitGenForDet=.false.
                HElemTemp=GetHElement2(CurrentDets(:,j),CurrentDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                CurrentH(1,j)=(REAL(HElemTemp%v,r2))-Hii
                CurrentH(2,j)=0.D0
            ENDIF

        enddo
            
        WRITE(6,*) "Initial allocation of excitation generators successful..."
        CALL FLUSH(6)
        CLOSE(17)
        
        RETURN

    END SUBROUTINE ReadFromPopsFile



! Based on SORTI, SORTPARTS sorts arrays of integers, representing the determinant the walkers are on
! It then takes all the corresponding info with it
! Dets is the array (length N) of integers to sort
! NElecs is the length (in numbers of integers) of each element of Dets
! Vectors of NewXXX will be sorted correspondingly
    SUBROUTINE SortParts(N,Dets,NElecs)
        TYPE(ExcitGenerator) :: ExcitTemp
        REAL*8 :: HTemp(2)
        INTEGER :: ICTemp
        INTEGER :: TempDet(NElecs)     !This stores a single element of the vector temporarily     
        LOGICAL :: WSignTemp
        INTEGER N,I,L,IR,J,NElecs
        INTEGER Dets(NElecs,N)
        INTEGER DETLT

        IF(N.LE.1) RETURN
        L=N/2+1 
        IR=N
10      CONTINUE
        IF(L.GT.1)THEN
            L=L-1
            CALL ICOPY(NElecs,Dets(1,L),1,TempDet,1)          !Copy Lth element to temp
            HTemp(:)=NewH(:,L)
            ICTemp=NewIC(L)
            WSignTemp=NewSign(L)
            CALL CopyExitgen(NewExcits(L),ExcitTemp)
        ELSE
            CALL ICOPY(NElecs,Dets(1,IR),1,TempDet,1)         !Copy IRth elements to temp
            HTemp(:)=NewH(:,IR)
            ICTemp=NewIC(IR)
            WSignTemp=NewSign(IR)
            CALL CopyExitgen(NewExcits(IR),ExcitTemp)

            CALL ICOPY(NElecs,Dets(1,1),1,Dets(1,IR),1)       !Copy 1st element to IRth element
            NewH(:,IR)=NewH(:,1)
            NewIC(IR)=NewIC(1)
            NewSign(IR)=NewSign(1)
            CALL CopyExitgen(NewExcits(1),NewExcits(IR))
            IR=IR-1
            IF(IR.EQ.1)THEN
                CALL ICOPY(NElecs,TempDet,1,Dets(1,1),1)        !Copy temp element to 1st element
                NewH(:,1)=HTemp(:)
                NewIC(1)=ICTemp
                NewSign(1)=WSignTemp
                CALL CopyExitgen(ExcitTemp,NewExcits(1))
                RETURN
            ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
                IF((DETLT(Dets(1,J),Dets(1,J+1),NElecs)).eq.-1) J=J+1
            ENDIF
            IF((DETLT(TempDet,Dets(1,J),NElecs)).eq.-1)THEN
                CALL ICOPY(NElecs,Dets(1,J),1,Dets(1,I),1)      !Copy Jth element to Ith element
                NewH(:,I)=NewH(:,J)
                NewIC(I)=NewIC(J)
                NewSign(I)=NewSign(J)
                CALL CopyExitgen(NewExcits(J),NewExcits(I))
                I=J
                J=J+J
            ELSE
                J=IR+1
            ENDIF
            GO TO 20
        ENDIF
        CALL ICOPY(NElecs,TempDet,1,Dets(1,I),1)                !Copy from temp element to Ith element
        NewH(:,I)=HTemp(:)
        NewIC(I)=ICTemp
        NewSign(I)=WSignTemp
        CALL CopyExitgen(ExcitTemp,NewExcits(I))
        GO TO 10
    END SUBROUTINE SortParts


!This routine copies an excitation generator from origExcit to NewExit, if the original claims that it is for the correct determinant
    SUBROUTINE CopyExitgen(OrigExit,NewExit)
        TYPE(ExcitGenerator) :: OrigExit,NewExit
        INTEGER :: ierr

        IF(Allocated(NewExit%ExcitData)) THEN
            DEALLOCATE(NewExit%ExcitData)
        ENDIF
        IF(.not.OrigExit%ExitGenForDet) THEN
!See if we actually want the excitation generator - is it for the correct determinant
            NewExit%ExitGenForDet=.false.
            RETURN
        ELSE
!We want to copy the excitation generator
            ALLOCATE(NewExit%ExcitData(OrigExit%nExcitMemLen),stat=ierr)
!            IF(OrigExit%nExcitMemLen.eq.0) THEN
!                CALL Stop_All("CopyExitgen","Problem allocating memory for new excit")
!            ENDIF
            IF(ierr.ne.0) CALL Stop_All("CopyExitgen","Problem with allocating memory for new excitation generator")
            CALL ICOPY(OrigExit%nExcitMemLen,OrigExit%ExcitData,1,NewExit%ExcitData,1)
!            NewExit%ExcitData(:)=OrigExit%ExcitData(:)
            NewExit%nExcitMemLen=OrigExit%nExcitMemLen
            NewExit%ExitGenForDet=.true.
        ENDIF

        RETURN

    END SUBROUTINE CopyExitgen

    SUBROUTINE SetupExitgen(nI,ExcitGen)
        TYPE(ExcitGenerator) :: ExcitGen
        INTEGER :: ierr,iMaxExcit,nExcitMemLen,nJ(NEl)
        INTEGER :: nI(NEl),nStore(6),iCount,IC!,Exitlevel,iGetExcitLevel
!        REAL*8 :: Prob

        IF(ExcitGen%ExitGenForDet) THEN
!The excitation generator is already allocated for the determinant in question - no need to recreate it
            IF(.not.Allocated(ExcitGen%ExcitData)) THEN
                CALL Stop_All("SetupExitgen","Excitation generator meant to already be set up")
            ENDIF

        ELSE

            IF(Allocated(ExcitGen%ExcitData)) THEN
                DEALLOCATE(ExcitGen%ExcitData)
            ENDIF

!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
            iMaxExcit=0
            CALL IAZZERO(nStore,6)
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGen%nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
            ALLOCATE(ExcitGen%ExcitData(ExcitGen%nExcitMemLen),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
            ExcitGen%ExcitData(1)=0
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGen%ExcitData,nJ,iMaxExcit,0,nStore,3)

!Check generation probabilities
!            CALL GenRandSymExcitIt3(nI,ExcitGen%ExcitData,nJ,Seed,IC,0,Prob,iCount)
!            Exitlevel=iGetExcitLevel(nI,HFDet,NEl)
!            IF(ABS((1.D0/iMaxExcit)-Prob).gt.1.D-07) WRITE(6,"I5,I5,I9,2G20.10") ExitLevel,IC,iMaxExcit,1.D0/iMaxExcit,Prob
            
!Indicate that the excitation generator is now correctly allocated.
            ExcitGen%ExitGenForDet=.true.

        ENDIF

    END SUBROUTINE SetupExitgen

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
!You can optionally specify the connection between the determinants if it is already known.
    INTEGER FUNCTION AttemptCreate(DetCurr,WSign,nJ,Prob,IC,Hij)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate
        LOGICAL :: WSign
        REAL*8 :: Prob,Ran2,rat
        REAL*8 , OPTIONAL :: Hij
        TYPE(HElement) :: rh

        IF(present(Hij)) THEN
!Connection between determinants is specified in Hij argument
            rh=HElement(Hij)
        ELSE
!Calculate off diagonal hamiltonian matrix element between determinants
            rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        ENDIF

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
        rat=Tau*abs(rh%v)/Prob

!If probability is > 1, then we can just create multiple children at the chosen determinant
        ExtraCreate=INT(rat)
        rat=rat-REAL(ExtraCreate)


!Stochastically choose whether to create or not according to Ran2
        IF(rat.gt.Ran2(Seed)) THEN
!Child is created - what sign is it?
            IF(WSign) THEN
!Parent particle is positive
                IF(real(rh%v).gt.0.D0) THEN
                    AttemptCreate=-1     !-ve walker created
                ELSE
                    AttemptCreate=1      !+ve walker created
                ENDIF

            ELSE
!Parent particle is negative
                IF(real(rh%v).gt.0.D0) THEN
                    AttemptCreate=1      !+ve walker created
                ELSE
                    AttemptCreate=-1     !-ve walker created
                ENDIF
            ENDIF

        ELSE
!No child particle created
            AttemptCreate=0
        ENDIF

        IF(ExtraCreate.ne.0) THEN
!Need to include the definitely create additional particles from a initial probability > 1

            IF(AttemptCreate.lt.0) THEN
!In this case particles are negative
                AttemptCreate=AttemptCreate-ExtraCreate
            ELSEIF(AttemptCreate.gt.0) THEN
!Include extra positive particles
                AttemptCreate=AttemptCreate+ExtraCreate
            ELSEIF(AttemptCreate.eq.0) THEN
!No particles were stochastically created, but some particles are still definatly created - we need to determinant their sign...
                IF(WSign) THEN
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreate=-1*ExtraCreate    !Additional particles are negative
                    ELSE
                        AttemptCreate=ExtraCreate       !Additional particles are positive
                    ENDIF
                ELSE
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreate=ExtraCreate
                    ELSE
                        AttemptCreate=-1*ExtraCreate
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        RETURN

    END FUNCTION AttemptCreate

!This function tells us whether we should kill the particle at determinant DetCurr
!If also diffusing, then we need to know the probability with which we have spawned. This will reduce the death probability.
!The function allows multiple births(if +ve shift) or deaths from the same particle.
!The returned number is the number of deaths if positive, and the number of births if negative.
    INTEGER FUNCTION AttemptDie(DetCurr,Kii)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),iKill
!        TYPE(HElement) :: rh,rhij
        REAL*8 :: Ran2,rat,Kii

!Calculate the diagonal hamiltonian matrix element for the determinant
!        rh=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!Subtract from the diagonal the value of the lowest hamiltonian matrix element
!        rh=rh-Hii

!Subtract the current value of the shift and multiply by tau
        rat=Tau*(Kii-DiagSft)

        iKill=INT(rat)
        rat=rat-REAL(iKill)

!Stochastically choose whether to die or not
        IF(abs(rat).gt.Ran2(Seed)) THEN
            IF(rat.ge.0.D0) THEN
!Depends whether we are trying to kill or give birth to particles.
                iKill=iKill+1
            ELSE
                iKill=iKill-1
            ENDIF
        ENDIF

        AttemptDie=iKill

        RETURN

    END FUNCTION AttemptDie

!Deallocate memory needed for the simulation
    SUBROUTINE DeallocFCIMCMem()
        INTEGER :: i
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocFCIMCMem'
        
        DEALLOCATE(WalkVecDets)
        CALL LogMemDealloc(this_routine,WalkVecDetsTag)
        DEALLOCATE(WalkVec2Dets)
        CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
        DEALLOCATE(WalkVecSign)
        CALL LogMemDealloc(this_routine,WalkVecSignTag)
        DEALLOCATE(WalkVec2Sign)
        CALL LogMemDealloc(this_routine,WalkVec2SignTag)
        DEALLOCATE(WalkVecIC)
        CALL LogMemDealloc(this_routine,WalkVecICTag)
        DEALLOCATE(WalkVec2IC)
        CALL LogMemDealloc(this_routine,WalkVec2ICTag)
        DEALLOCATE(WalkVecH)
        CALL LogMemDealloc(this_routine,WalkVecHTag)
        DEALLOCATE(WalkVec2H)
        CALL LogMemDealloc(this_routine,WalkVec2HTag)

        IF(TResumFCIMC) THEN
            DEALLOCATE(GraphRhoMat)
            CALL LogMemDealloc(this_routine,GraphRhoMatTag)
            DEALLOCATE(GraphVec)
            CALL LogMemDealloc(this_routine,GraphVecTag)
            DEALLOCATE(GraphKii)
            CALL LogMemDealloc(this_routine,GraphKiiTag)
            DEALLOCATE(DetsinGraph)
            CALL LogMemDealloc(this_routine,DetsinGraphTag)
        ENDIF

        DEALLOCATE(HFDet)
        CALL LogMemDealloc(this_routine,HFDetTag)
        DEALLOCATE(HFExcit%ExcitData)
        do i=1,MaxWalkers
            IF(Allocated(WalkVecExcits(i)%ExcitData)) THEN
                DEALLOCATE(WalkVecExcits(i)%ExcitData)
            ENDIF
            IF(Allocated(WalkVec2Excits(i)%ExcitData)) THEN
                DEALLOCATE(WalkVec2Excits(i)%ExcitData)
            ENDIF
        enddo
        DEALLOCATE(WalkVecExcits)
        DEALLOCATE(WalkVec2Excits)

    END SUBROUTINE DeallocFCIMCMem
        

END MODULE FciMCMod


!This is a determinant comparison function
!DetA and DetB are two determinants. Function will return 0 if they are the same, -1 if A.lt.B and +1 if A.gt.B when ordered.
INTEGER FUNCTION DetLT(DetA,DetB,NEl)
    IMPLICIT NONE
    INTEGER :: DetA(NEl),DetB(NEl),NEl,i

    do i=1,NEl
        IF(DetA(i).lt.DetB(i)) THEN
            DetLT=-1
            RETURN
        ELSEIF(DetA(i).gt.DetB(i)) THEN
            DetLT=+1
            RETURN
        ENDIF
    enddo
    DetLT=0
    RETURN

END FUNCTION DetLT

!A function to get the factorial of the number Num
INTEGER*8 FUNCTION Fact(Num)
    IMPLICIT NONE
    INTEGER :: Num,i

    Fact=Num
    do i=Num-1,2,-1
        Fact=Fact*i
    enddo
    RETURN

END FUNCTION Fact

!Below is the commented out old developmental serial version of the code. This contains more features, though has turned into a bit of a mess.
!It is not deleted as may want salvaging later, but it is now easier to rewrite a seperate more efficient code!

!MODULE FciMCMod
!    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr
!    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
!    USE Calc , only : TReadPops,ScaleWalkers,TMCExcitSpace,NoMCExcits,TStartMP1
!    USE Calc , only : GrowMaxFactor,CullFactor,TMCDets,TNoBirth,Lambda,TDiffuse,FlipTauCyc,TFlipTau
!    USE Calc , only : TExtraPartDiff,TFullUnbias,TNodalCutoff,NodalCutoff,TNoAnnihil,TMCDiffusion
!    USE Calc , only : NDets,RhoApp,TResumFCIMC,NEquilSteps
!    USE Determinants , only : FDet,GetHElement2
!    USE DetCalc , only : NMRKS
!    USE Integrals , only : fck,NMax,nMsh,UMat
!    USE Logging , only : TPopsFile,TCalcWavevector,WavevectorPrint,TDetPops
!    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
!    USE HElem
!    IMPLICIT NONE
!    SAVE
!
!    INTEGER, PARAMETER :: r2=kind(0.d0)
!
!    INTEGER , ALLOCATABLE , TARGET :: WalkVecDets(:,:),WalkVec2Dets(:,:)
!    LOGICAL , ALLOCATABLE , TARGET :: WalkVecSign(:),WalkVec2Sign(:)
!    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0
!
!!Pointers to point at the correct arrays for use
!    INTEGER , POINTER :: CurrentDets(:,:)
!    LOGICAL , POINTER :: CurrentSign(:)
!    INTEGER , POINTER :: NewDets(:,:)
!    LOGICAL , POINTER :: NewSign(:)
!
!!InitPopsVector is used to store the initial populations of the determinants. This can be used for comparison when TNoBirth is true.
!    INTEGER , ALLOCATABLE :: InitPops(:)
!    INTEGER :: InitPopsTag=0
!
!    INTEGER , ALLOCATABLE :: ExcitStore(:,:)
!    INTEGER :: ExcitStoreTag=0
!
!    REAL*8 , ALLOCATABLE :: TransMat(:,:),PopsVec(:)
!    INTEGER :: TransMatTag=0,PopsVecTag=0,SizeOfSpace
!
!!MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
!    INTEGER :: MemoryFac=1000
!
!    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps
!    INTEGER :: exFlag=3
!
!!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!!NoCulls is the number of culls in a given shift update cycle
!    INTEGER :: NoCulls=0
!!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!!Only 10 culls/growth increases are allowed in a given shift cycle
!    INTEGER :: CullInfo(10,3)
!
!    REAL*8 :: GrowRate,DieRat,MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector
!    REAL*8 :: ProjectionE,SumENum,SumE,ProjectionEInst
!
!    INTEGER :: CycwNoHF     !Count the number of iterations which don't have any walkers on HF - in this case, ignore running average of energy, but count it so it is not included in demonimator
!
!    INTEGER*8 :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
!    REAL*8 :: MeanExcitLevel    
!    INTEGER :: MinExcitLevel
!    INTEGER :: MaxExcitLevel
!    INTEGER :: SumWalkersCyc    !Sum of all walkers in an update cycle
!
!    INTEGER :: NetPositive
!
!    TYPE(HElement) :: Hii,rhii,FZero
!
!    TYPE ExcitGenerator
!        INTEGER , ALLOCATABLE :: ExcitData(:)      !This stores the excitation generator
!        INTEGER :: nStore(6)
!    END TYPE
!    TYPE(ExcitGenerator) , ALLOCATABLE :: ExcitGens(:)   !This will store the excitation generators for the walkers in MCDiffusion
!
!    REAL*8 , ALLOCATABLE :: GraphRhoMat(:,:)    !This stores the rho matrix for the graphs in resumFCIMC
!    INTEGER :: GraphRhoMatTag=0
!
!    REAL*8 , ALLOCATABLE :: GraphVec(:)         !This stores the final components for the propagated graph in ResumFCIMC
!    INTEGER :: GraphVecTag=0
!
!    INTEGER , ALLOCATABLE :: DetsinGraph(:,:)   !This stores the determinants in the graph created for ResumFCIMC
!    INTEGER :: DetsinGraphTag=0
!
!    REAL*8 :: RootExcitProb     !This is the probability of generating an excitation from the current root in the ResumFCIMC current graph.
!
!    LOGICAL :: TResumAllConns=.false.    !If set, this means that all possible connections to a walker are created for the graph.
!
!!    REAL*8 :: SumCreateProb     !This is the culmulative probability of creating particles. It can be used as a more accurate determination of shift change for small walker numbers.
!
!    contains
!
!    SUBROUTINE FciMC(Weight,Energyxw)
!        Use MCDets, only : MCDetsCalc
!        IMPLICIT NONE
!        TYPE(HDElement) :: Weight,Energyxw
!        INTEGER :: i,j,iSub,WalkOnDet,DetLT,DetCurr(NEl),ExpectedDets
!        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'
!        TYPE(HElement) :: Hamii
!
!        if(NMCyc.lt.0.or.TMCDets) then
!
!           CALL MCDetsCalc(FDet,G_VMC_Seed,-NMCyc,Tau,SftDamp,10000*InitWalkers,InitWalkers,StepsSft,DiagSft,GrowMaxFactor,CullFactor)
!            Weight=1.d0
!            Energyxw=1.d0
!            return
!        endif
!        CALL TISET('FCIMC',iSub)
!
!        IF(TDiffuse) THEN
!            IF((.NOT.TMCExcitSpace).or.(NoMCExcits.ne.1)) THEN
!                CALL Stop_All("FCIMC","Diffusion can only work with one MCExcitSpace")
!            ENDIF
!            IF((Lambda.gt.1.D0).or.(Lambda.lt.0.D0)) THEN
!                CALL Stop_All("FCIMC","Diffusion coefficient must be between zero and one")
!            ENDIF
!        ENDIF
!
!        IF(HElementSize.gt.1) THEN
!            CALL Stop_All("FCIMC","StarDiagMC cannot function with complex orbitals.")
!        ENDIF
!        
!        OPEN(15,file='FCIMCStats',status='unknown')
!
!!Initialise random number seed
!        Seed=G_VMC_Seed
!!Calculate Hii
!        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
!!Initialise variables for calculation of the running average
!        ProjectionE=0.D0
!        ProjectionEInst=0.D0
!        SumE=0.D0
!        SumENum=0.D0
!        SumNoatHF=0
!        CycwNoHF=0
!        MeanExcitLevel=0.D0
!        MinExcitLevel=Nel+10
!        MaxExcitLevel=0
!        SumWalkersCyc=0
!!        SumCreateProb=0.D0
!
!        IF(TMCDiffusion) THEN
!            CALL MCDiffusion()
!            CLOSE(15)
!            RETURN
!        ENDIF
!
!        IF(TResumFciMC) THEN
!            IF(NDets.lt.0) TResumAllConns=.true.
!            IF(NDets.eq.1.or.NDets.eq.0) CALL Stop_All("CreateGraph","Graphsize too small...")
!            WRITE(6,*) ""
!            WRITE(6,*) "Performing FCIMC...."
!        ELSE
!            WRITE(6,*) ""
!            WRITE(6,*) "Performing FCIMC...."
!        ENDIF
!
!
!        IF(TCalcWaveVector) THEN
!            WRITE(6,*) "Wavevector calculation is only available in star MC..."
!        ENDIF
!
!        IF(TReadPops) THEN
!            
!            IF(TNoBirth) THEN
!                WRITE(6,*) "Cannot use NOBIRTHS with READPOPS"
!                STOP "Cannot use NOBIRTHS with READPOPS"
!            ENDIF
!
!            WRITE(6,*) "Reading in POPSFILE to restart calculation..."
!            OPEN(17,FILE='POPSFILE',Status='old')
!            READ(17,*) InitWalkers
!            READ(17,*) DiagSft
!            READ(17,*) PreviousNMCyc
!            WRITE(6,*) "Initial number of walkers read to be: ", InitWalkers
!        ELSE
!            WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
!        ENDIF
!
!        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
!        IF(TFlipTau) THEN
!            WRITE(6,*) "Flipping the sign of Tau once every :", FlipTauCyc
!        ENDIF
!
!!        IF(DiagSft.gt.0.D0) THEN
!!            CALL Stop_All("StarDiagMC","Intial value of DiagSft should be negative.")
!!        ELSE
!            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
!!        ENDIF
!
!        CALL InitFCIMCCalc()
!
!        IF(.NOT.TNoBirth) THEN
!!Print out initial starting configurations
!            WRITE(6,*) ""
!            WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers        Proj.E      Net+veWalk     Proj.E-Inst   MeanExcitLevel   MinExcitLevel   MaxExcitLevel"
!            WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers         Proj.E      Net+veWalk     Proj.E-Inst   MeanExcitLevel   MinExcitLevel   MaxExcitLevel"
!!TotWalkersOld is the number of walkers last time the shift was changed
!            IF(TReadPops) THEN
!                WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
!                WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
!            ELSE
!                WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7,F16.7,2I5)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst,MeanExcitLevel,MinExcitLevel,MaxExcitLevel
!                WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7,F16.7,2I5)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst,MeanExcitLevel,MinExcitLevel,MaxExcitLevel
!            ENDIF
!        ENDIF
!        
!!Reset TotWalkersOld so that it is the number of walkers now
!        TotWalkersOld=TotWalkers
!
!!Start MC simulation...
!        do Iter=1,NMCyc
!            
!            IF(TResumFciMC) THEN
!                CALL PerformResumFCIMCyc()
!            ELSE
!                CALL PerformFCIMCyc()
!            ENDIF
!
!!            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
!!            CALL FLUSH(6)
!
!            IF(TFlipTau) THEN
!!If we are flipping the sign of tau every FlipTauCyc cycles
!                IF(mod(Iter,FlipTauCyc).eq.0) THEN
!                    Tau=Tau*(-1.D0)
!                    DiagSft=DiagSft*(-1.D0)
!                ENDIF
!            ENDIF
!
!            IF(mod(Iter,StepsSft).eq.0) THEN
!
!                IF(TNoBirth) THEN
!!If we have TNoBirth flag, the populations on each determinant simply exponentially decay, and so at each StepsSft, we want to know the population on each determinant
!!When creating wavevector, store the initial number of particles on each determinant (order them).
!                    do i=1,NoComps
!                        do j=1,NEl
!                            DetCurr(j)=ExcitStore(j,i)
!                        enddo
!                        WalkOnDet=0
!                        do j=1,TotWalkers
!                            IF(DetLT(DetCurr,CurrentDets(:,j),NEl).eq.0) THEN
!!Walker is on chosen determinant
!                                WalkOnDet=WalkOnDet+1
!                            ENDIF
!                        enddo
!
!                        Hamii=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                        Hamii=Hamii-Hii
!!The expected number of determinants now is the original number on each determinant * exp(-tau*Iter*(Hii-DiagSft)
!                        ExpectedDets=exp(-Tau*Iter*((Hamii%v)-DiagSft))*InitPops(i)
!
!                        WRITE(6,*) Iter,i,WalkOnDet,ExpectedDets
!                        WRITE(15,*) Iter,i,WalkOnDet,ExpectedDets
!                        CALL FLUSH(15)
!                        CALL FLUSH(6)
!
!                    enddo
!                    WRITE(6,*) ""
!                    WRITE(15,*) ""
!
!!This can be compared to the actual number of particles on each of the determinants.
!!This will only work when we start with the MP1 wavefunction.
!
!                ELSE
!!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
!                    CALL UpdateDiagSft()
!
!                    IF(TResumFCIMC) THEN
!                        ProjectionE=SumENum/REAL(SumNoatHF,r2)
!                        MeanExcitLevel=MeanExcitLevel/real(SumWalkersCyc,r2)
!                    ENDIF
!
!!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
!                    IF(TReadPops) THEN
!                        WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7,I)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
!                        WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
!                    ELSE
!                        IF(Tau.gt.0.D0) THEN
!                            WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7,F16.7,2I5)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst,MeanExcitLevel,MinExcitLevel,MaxExcitLevel
!                            WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7,F16.7,2I5)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst,MeanExcitLevel,MinExcitLevel,MaxExcitLevel
!                        ELSE
!                            WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
!                            WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
!                        ENDIF
!                    ENDIF
!
!                    MinExcitLevel=Nel+10
!                    MaxExcitLevel=0
!                    MeanExcitLevel=0.D0
!                    SumWalkersCyc=0
!
!                    IF(TDetPops) THEN
!
!                        do i=1,SizeOfSpace
!                            WRITE(67,"(F10.5,$)") TransMat(1,i)
!                        enddo
!                        WRITE(67,*) ""
!
!                    ENDIF
!
!
!                    CALL FLUSH(15)
!                    CALL FLUSH(6)
!
!!Reset TotWalkersOld so that it is the number of walkers now
!                    TotWalkersOld=TotWalkers
!                ENDIF
!
!!                WRITE(6,*) DieRat
!
!            ENDIF
!
!!End of MC cycle
!        enddo
!
!        IF(TPopsFile) THEN
!!Print out current state of simulation, so it can be restarted if desired...
!            OPEN(16,file='POPSFILE',status='unknown')
!            WRITE(16,*) TotWalkers, "   TOTWALKERS"
!            WRITE(16,*) DiagSft, "   DIAG SHIFT"
!            IF(TReadPops) THEN
!                WRITE(16,*) NMCyc+PreviousNMCyc, "   NO. CYCLES"
!            ELSE
!                WRITE(16,*) NMCyc, "   MC CYCLES"
!            ENDIF
!            do i=1,TotWalkers
!                WRITE(16,*) CurrentDets(:,i),CurrentSign(i)
!            enddo
!            CLOSE(16)
!        ENDIF
!
!        Weight=HDElement(0.D0)
!        Energyxw=HDElement(DiagSft)
!
!!Deallocate memory
!        IF(TNoBirth) THEN
!            DEALLOCATE(ExcitStore)
!            CALL LogMemDealloc(this_routine,ExcitStoreTag)
!            DEALLOCATE(InitPops)
!            CALL LogMemDealloc(this_routine,InitPopsTag)
!        ENDIF
!        IF(TDetPops) THEN
!            DEALLOCATE(ExcitStore)
!            CALL LogMemDealloc(this_routine,ExcitStoreTag)
!            DEALLOCATE(PopsVec)
!            CALL LogMemDealloc(this_routine,PopsVecTag)
!            DEALLOCATE(TransMat)
!            CALL LogMemDealloc(this_routine,TransMatTag)
!        ENDIF
!        IF(TResumFciMC) THEN
!            DEALLOCATE(GraphRhoMat)
!            CALL LogMemDealloc(this_routine,GraphRhoMatTag)
!            DEALLOCATE(GraphVec)
!            CALL LogMemDealloc(this_routine,GraphVecTag)
!            DEALLOCATE(DetsinGraph)
!            CALL LogMemDealloc(this_routine,DetsInGraphTag)
!        ENDIF
!        DEALLOCATE(WalkVecDets)
!        CALL LogMemDealloc(this_routine,WalkVecDetsTag)
!        DEALLOCATE(WalkVec2Dets)
!        CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
!        DEALLOCATE(WalkVecSign)
!        CALL LogMemDealloc(this_routine,WalkVecSignTag)
!        DEALLOCATE(WalkVec2Sign)
!        CALL LogMemDealloc(this_routine,WalkVec2SignTag)
!
!        CLOSE(15)
!
!        CALL TIHALT('FCIMC',iSub)
!
!        RETURN
!
!    END SUBROUTINE FciMC
!
!!This performs a resummed FCIMC calculation, where small graphs are created from each walker at the space, and the true matrix propagated around
!    SUBROUTINE PerformResumFCIMCyc()
!        IMPLICIT NONE
!        TYPE(ExcitGenerator) :: ExcitGen
!        INTEGER :: TotWalkersNew,ExcitLevel,nExcitMemLen,iGetExcitLevel,VecSlot,i,j
!        INTEGER :: TotExcits,TotComps
!        TYPE(HElement) :: Hij0
!        REAL*8 :: rat
!
!
!!VecSlot indicates the next free position in NewDets
!        VecSlot=1
!
!        do j=1,TotWalkers
!
!            ExcitLevel=iGetExcitLevel(CurrentDets(:,j),FDet(:),NEl)
!            MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,r2)
!            IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
!            IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
!            IF(ExcitLevel.eq.0) THEN
!                IF(CurrentSign(j)) THEN
!                    SumNoatHF=SumNoatHF+1
!                ELSE
!                    SumNoatHF=SumNoatHF-1
!                ENDIF
!            ELSEIF(ExcitLevel.eq.2) THEN
!!At double excit - sum in energy
!                Hij0=GetHElement2(FDet(:),CurrentDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
!                IF(CurrentSign(j)) THEN
!                    SumENum=SumENum+REAL(Hij0%v,r2)
!                ELSE
!                    SumENum=SumENum-REAL(Hij0%v,r2)
!                ENDIF
!            ENDIF
!
!!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.)
!            CALL SetupExitgen(CurrentDets(:,j),ExcitGen,nExcitMemLen,TotExcits)
!            CALL CreateGraph(CurrentDets(:,j),ExcitGen,TotExcits,TotComps)     !Construct a graph of size NDets with the current particle at the root
!            CALL ApplyRhoMat(TotComps)  !Successivly apply the rho matrix to the particle RhoApp times
!            CALL CreateNewParts(CurrentSign(j),VecSlot,TotComps)   !Create new particles according to the components of GraphVec, and put them into NewVec
!
!!Destroy excitation generators for current walker
!            DEALLOCATE(ExcitGen%ExcitData)
!        
!!Finish cycling over walkers
!        enddo
!
!        SumWalkersCyc=SumWalkersCyc+TotWalkers
!
!!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
!        TotWalkersNew=VecSlot-1
!        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
!        IF(rat.gt.0.9) THEN
!            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
!        ENDIF
!        
!        IF(TNodalCutoff.and.(NodalCutoff.lt.0.D0)) THEN
!!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff
!!then the net number of walkers must be the same sign, or it is set to zero walkers, and they are killed.
!            CALL TestWavevectorNodes(TotWalkersNew,2)
!        ENDIF
!
!        IF(TNoAnnihil) THEN
!!We are not annihilating particles - this will make things much quicker.
!
!!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
!            IF(associated(CurrentDets,target=WalkVecDets)) THEN
!                CurrentDets=>WalkVec2Dets
!                CurrentSign=>WalkVec2Sign
!                NewDets=>WalkVecDets
!                NewSign=>WalkVecSign
!            ELSE
!                CurrentDets=>WalkVecDets
!                CurrentSign=>WalkVecSign
!                NewDets=>WalkVec2Dets
!                NewSign=>WalkVec2Sign
!            ENDIF
!
!            TotWalkers=TotWalkersNew
!
!        ELSE
!!This routine now cancels down the particles with opposing sign on each determinant
!!This routine does not necessarily need to be called every Iter, but it does at the moment, since it is the only way to 
!!transfer the residual particles back onto CurrentDets
!            CALL AnnihilatePairs(TotWalkersNew)
!!            WRITE(6,*) "Number of annihilated particles= ",TotWalkersNew-TotWalkers
!        ENDIF
!
!
!        IF(TNodalCutoff.and.(NodalCutoff.ge.0.D0)) THEN
!!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff
!!then the net number of walkers must be the same sign, or it is set to zero walkers, and they are killed.
!            CALL TestWavevectorNodes(TotWalkers,1)
!        ENDIF
!
!        
!        IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!!Particle number is too large - kill them randomly
!
!!Log the fact that we have made a cull
!            NoCulls=NoCulls+1
!            IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!            CullInfo(NoCulls,1)=TotWalkers
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!               THIS IS DONE INCORRECTLY NOW
!            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!
!            WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number..."
!            WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number..."
!            WRITE(6,"(A,F8.2)") "Population will reduce by a factor of ",CullFactor
!            CALL ThermostatParticles(.true.)
!
!!Need to reduce totwalkersOld, so that the change in shift is also reflected by this
!!The Shift is no longer calculated like this...
!!            TotWalkersOld=nint((TotWalkersOld+0.D0)/CullFactor)
!
!        ELSEIF((TotWalkers.lt.(InitWalkers/2)).and.(.NOT.TNoBirth)) THEN
!!Particle number is too small - double every particle in its current position
!
!!Log the fact that we have made a cull
!            NoCulls=NoCulls+1
!            IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!            CullInfo(NoCulls,1)=TotWalkers
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!               THIS IS DONE INCORRECTLY NOW!!
!            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!            
!            WRITE(6,*) "Doubling particle population to increase total number..."
!            CALL ThermostatParticles(.false.)
!
!!Need to increase TotWalkersOld, so that the change in shift is also reflected by this
!!The shift is no longer calculated like this...
!!            TotWalkersOld=TotWalkersOld*2
!
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE PerformResumFCIMCyc
!
!    SUBROUTINE CreateNewParts(WSign,VecSlot,TotComps)
!        IMPLICIT NONE
!        LOGICAL :: WSign
!        INTEGER :: i,j,VecSlot,Create,ExtraCreate,StochCreate,iKill,TotComps,Components
!        REAL*8 :: Ran2,rat
!
!        IF(TResumAllConns) THEN
!            Components=TotComps
!        ELSE
!            Components=NDets
!        ENDIF
!
!        do i=1,Components
!
!            IF((i.ne.1).and.(.not.TResumAllConns)) THEN
!!Augment the list of creation probabilities by dividing through by the probability of creating a graph with that excitation in it.
!                GraphVec(i)=GraphVec(i)/((NDets-1)*RootExcitProb)
!            ENDIF
!
!            Create=INT(abs(GraphVec(i)))
!
!            rat=abs(GraphVec(i))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
!            IF(rat.gt.Ran2(Seed)) Create=Create+1
!            IF(.not.WSign) Create=-Create
!            IF(GraphVec(i).lt.0.D0) Create=-Create
!!            IF((GraphVec(i).lt.0.D0).and.(i.eq.1)) THEN
!!                CALL Stop_All("CreateParts", "trying to create opposite signed particle on root")
!!            ENDIF
!!            IF(i.eq.1) WRITE(6,*) GraphVec(i),Create, WSign
!
!!Now actually create the particles in NewDets and NewSign
!            do j=1,abs(Create)
!                NewDets(:,VecSlot)=DetsInGraph(:,i)
!                IF(Create.lt.0) THEN
!                    NewSign(VecSlot)=.false.
!                ELSE
!                    NewSign(VecSlot)=.true.
!                ENDIF
!                VecSlot=VecSlot+1
!            enddo
!        enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   T E S T - to reduce to original spawning when RhoApp=1 and NDets=2    !!!!!!!!!!!!!!!!!!!!!!!!
!
!!        do i=1,2
!!            GraphVec(i)=GraphVec(i)/GraphVecProbs(i)
!!        enddo
!!
!!!First deal with spawning prob.....GraphVec(2) has -tau*Hij/Pj    spawning prob is tau*abs(Hij)/Pj so take abs value
!!        rat=abs(GraphVec(2))
!!        ExtraCreate=INT(rat)
!!        rat=rat-REAL(ExtraCreate)
!!        IF(rat.gt.Ran2(Seed)) THEN
!!            StochCreate=1
!!        ELSE
!!            StochCreate=0
!!        ENDIF
!!        IF(WSign) THEN
!!!Parent particle is positive
!!            IF(GraphVec(2).lt.0.D0) THEN
!!                Create=-StochCreate-ExtraCreate     !-ve walker created
!!            ELSE
!!                Create=StochCreate+ExtraCreate      !+ve walker created
!!            ENDIF
!!
!!        ELSE
!!!Parent particle is negative
!!            IF(GraphVec(2).lt.0.D0) THEN
!!                Create=StochCreate+ExtraCreate      !+ve walker created
!!            ELSE
!!                Create=-StochCreate-ExtraCreate     !-ve walker created
!!            ENDIF
!!        ENDIF
!!
!!        do j=1,abs(Create)
!!            NewDets(:,VecSlot)=DetsInGraph(:,2)
!!            IF(Create.lt.0) THEN
!!                NewSign(VecSlot)=.false.
!!            ELSE
!!                NewSign(VecSlot)=.true.
!!            ENDIF
!!            VecSlot=VecSlot+1
!!        enddo
!!
!!!Now deal with death......GraphVec(1) has 1-tau*((Hii-H00)-Sft)    Death probability is tau*((Hii-H00)-Sft)
!!        rat=1.D0-GraphVec(1)        !This is now death prob
!!        iKill=INT(rat)
!!        rat=rat-REAL(iKill)
!!        IF(abs(rat).gt.Ran2(Seed)) THEN
!!            IF(rat.ge.0.D0) THEN
!!!Depends whether we are trying to kill or give birth to particles.
!!                iKill=iKill+1
!!            ELSE
!!                iKill=iKill-1
!!            ENDIF
!!        ENDIF
!!
!!        IF(iKill.eq.0) THEN
!!!Want to keep particle - copy accross
!!            NewDets(:,VecSlot)=DetsInGraph(:,1)
!!            IF(WSign) THEN
!!                NewSign(VecSlot)=.true.
!!            ELSE
!!                NewSign(VecSlot)=.false.
!!            ENDIF
!!            VecSlot=VecSlot+1
!!        ELSEIF(iKill.eq.1) THEN
!!!Don't copy accross
!!        ELSE
!!            WRITE(6,*) Iter,iKill,-(1.D0-GraphVec(1)),GraphVec(1)
!!            CALL Stop_All("CreateParts","Trying to kill/survive > 1")
!!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE CreateNewParts
!
!!This applies the rho matrix successive times to a root determinant. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
!    SUBROUTINE ApplyRhoMat(TotComps)
!        IMPLICIT NONE
!        REAL*8 , ALLOCATABLE :: TempVec(:)
!        INTEGER :: i,j,k,TotComps,ierr,Components
!        
!        IF(TResumAllConns) THEN
!!TotComps is the total dimension of the matrix if all excitations included (i.e. iMaxExcit+1)
!            Components=TotComps
!        ELSE
!            Components=NDets
!        ENDIF
!
!        ALLOCATE(TempVec(Components),stat=ierr)
!        CALL AZZERO(GraphVec,Components)
!        CALL AZZERO(TempVec,Components)
!
!        GraphVec(1)=1.D0    !Set the initial vector to be 1 at the root (i.e. for one walker initially)
!
!        do i=1,RhoApp   !Cycle over the number of times we want to apply the rho matrix
!
!            CALL DGEMV('n',Components,Components,1.D0,GraphRhoMat,Components,GraphVec,1,0.D0,TempVec,1)
!            CALL DCOPY(Components,TempVec,1,GraphVec,1)
!            CALL AZZERO(TempVec,Components)
!
!!            do j=1,Components
!!                TempVec(j)=0.D0
!!                do k=1,Components
!!                    TempVec(j)=TempVec(j)+GraphRhoMat(j,k)*GraphVec(k)
!!                enddo
!!            enddo
!!            do j=1,Components
!!                GraphVec(j)=TempVec(j)
!!            enddo
!            
!        enddo
!
!        DEALLOCATE(TempVec)
!
!        RETURN
!
!    END SUBROUTINE ApplyRhoMat
!            
!
!
!!This creates a graph of size NDets with all determinants attached to nI. It also forms the matrix for it, and puts it in GraphRhoMat, with the dets used in DetsInGraph(NEl,NDets)
!    SUBROUTINE CreateGraph(nI,ExcitGen,TotExcits,TotComps)
!        IMPLICIT NONE
!        INTEGER :: nI(NEl),IC,iCount,nJ(NEl),i,j,Attempts,TotExcits,ierr,TotComps
!        REAL*8 :: Prob,RemovedProb
!        TYPE(ExcitGenerator) :: ExcitGen
!        TYPE(HElement) :: Hamii,Hamij
!        LOGICAL :: SameDet,CompiPath
!
!        IF(TResumAllConns) THEN
!            TotComps=TotExcits+1
!            IF(ALLOCATED(GraphRhoMat)) DEALLOCATE(GraphRhoMat)
!            IF(ALLOCATED(GraphVec)) DEALLOCATE(GraphVec)
!            IF(ALLOCATED(DetsinGraph)) DEALLOCATE(DetsinGraph)
!            ALLOCATE(GraphRhoMat(TotComps,TotComps),stat=ierr)
!            ALLOCATE(GraphVec(TotComps),stat=ierr)
!            ALLOCATE(DetsinGraph(NEl,TotComps),stat=ierr)
!            CALL AZZERO(GraphRhoMat,TotComps**2)
!            CALL IAZZERO(DetsinGraph,NEl*TotComps)
!        ELSE
!            CALL AZZERO(GraphRhoMat,NDets**2)
!        ENDIF
!
!        DetsInGraph(:,1)=nI(:)
!        RemovedProb=0.D0        !RemovedProb is the sum of the probabilities of the excitations already picked. This allows for unbiasing when we only select distinct determinants.
!
!!Find diagonal element for root determinant
!        Hamii=GetHElement2(nI,nI,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!!        GraphRhoMat(1,1)=1.D0-Tau*((REAL(Hamii%v,r2)-REAL(Hii%v,r2))-(DiagSft/REAL(RhoApp,r2)))
!        GraphRhoMat(1,1)=1.D0-Tau*((REAL(Hamii%v,r2)-REAL(Hii%v,r2))-DiagSft)
!
!        IF(TResumAllConns) THEN
!!We want to run through all possible connections to nI...
!            i=2
!            do while(.true.)
!                CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.false.,ExcitGen%ExcitData,nJ,IC,0,ExcitGen%nStore,exFlag)
!                IF(nJ(1).eq.0) EXIT
!!Determinant is distinct - add it
!                IF(i.gt.TotComps) CALL Stop_All("CreateGraph","Problem creating graph")
!                DetsInGraph(:,i)=nJ(:)
!
!!First find connection to root
!                Hamij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!                GraphRhoMat(1,i)=-Tau*REAL(Hamij%v,r2)
!                GraphRhoMat(i,1)=GraphRhoMat(1,i)
!
!!Then find connection to other determinants
!                do j=2,(i-1)
!                    Hamij=GetHElement2(nJ,DetsInGraph(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
!                    GraphRhoMat(i,j)=-Tau*REAL(Hamij%v,r2)
!                    GraphRhoMat(j,i)=GraphRhoMat(i,j)
!                enddo
!
!!Find diagonal element
!                Hamii=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!!                GraphRhoMat(i,i)=1.D0-(Tau*((REAL(Hamii%v,r2)-REAL(Hii%v,r2))-(DiagSft/REAL(RhoApp,r2))))
!                GraphRhoMat(i,i)=1.D0-Tau*((REAL(Hamii%v,r2)-REAL(Hii%v,r2))-DiagSft)
!                i=i+1   !Increment counter
!            
!            enddo   !loop over excitations
!
!        ELSE
!!Here we are choosing a stochastic graph with a smaller number of excitations chosen at random
!        
!            i=2
!            do while(i.le.NDets)    !loop until all connections found
!
!                CALL GenRandSymExcitIt3(nI,ExcitGen%ExcitData,nJ,Seed,IC,0,Prob,iCount)
!
!                SameDet=.false.
!                do j=2,(i-1)
!                    IF(CompiPath(nJ,DetsInGraph(:,j),NEl)) THEN
!!determinants are the same - ignore them
!                        SameDet=.true.
!                        Attempts=Attempts+1     !Increment the attempts counter
!                        IF(Attempts.gt.100) CALL Stop_All("CreateGraph","More than 100 attempts needed to grow graph")
!                        EXIT
!                    ENDIF
!                enddo
!
!                IF(.not.SameDet) THEN
!!Store the unbiased probability of generating excitations from this root - check that it is the same as other excits generated
!                    IF(i.eq.2) THEN
!                        RootExcitProb=Prob
!                    ELSE
!                        IF(abs(Prob-RootExcitProb).gt.1.D-07) THEN
!                            CALL Stop_All("CreateGraph","Excitation probabilities are not uniform - problem here...")
!                        ENDIF
!                    ENDIF
!
!!Determinant is distinct - add it
!                    DetsInGraph(:,i)=nJ(:)
!
!!First find connection to root
!                    Hamij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!                    GraphRhoMat(1,i)=-Tau*REAL(Hamij%v,r2)
!                    GraphRhoMat(i,1)=GraphRhoMat(1,i)
!
!!Then find connection to other determinants
!                    do j=2,(i-1)
!                        Hamij=GetHElement2(nJ,DetsInGraph(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
!                        GraphRhoMat(i,j)=-Tau*REAL(Hamij%v,r2)
!                        GraphRhoMat(j,i)=GraphRhoMat(i,j)
!                    enddo
!
!!Find diagonal element
!                    Hamii=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!!                GraphRhoMat(i,i)=1.D0-(Tau*((REAL(Hamii%v,r2)-REAL(Hii%v,r2))-(DiagSft/REAL(RhoApp,r2))))
!                    GraphRhoMat(i,i)=1.D0-Tau*((REAL(Hamii%v,r2)-REAL(Hii%v,r2))-DiagSft)
!
!                    i=i+1   !increment the excit counter
!                    Attempts=0      !Reset the attempts counter
!
!                ENDIF
!
!            enddo
!
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE CreateGraph
!                
!
!!This is the heart of FCIMC, where the MC Cycles are performed
!    SUBROUTINE PerformFCIMCyc()
!        IMPLICIT NONE
!        INTEGER :: VecSlot,i,j,k,l,DetCurr(NEl),iMaxExcit,nExcitMemLen,nStore(6)
!        INTEGER :: nJ(NEl),ierr,nExcitTag=0,IC,Child,iSubCyc,TotWalkersNew,iCount
!        REAL*8 :: Prob,rat,Kik
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
!        LOGICAL :: WSign
!        LOGICAL :: KeepOrig
!        INTEGER :: CreateAtI,CreateAtJ,tocopy
!        CHARACTER(len=*), PARAMETER :: this_routine='PerformFCIMCyc'
!        
!        CALL TISET('MCyc',iSubCyc)
!        
!!VecSlot indicates the next free position in NewDets
!        VecSlot=1
!
!        do j=1,TotWalkers
!!j runs through all current walkers
!            do k=1,NEl
!                DetCurr(k)=CurrentDets(k,j)
!            enddo
!
!!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
!            iMaxExcit=0
!            CALL IAZZERO(nStore,6)
!            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!            nExcit(1)=0
!            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!
!            IF(TMCExcitSpace) THEN
!!Excitations are picked stochastically
!
!                do i=1,NoMCExcits
!
!                    CALL GenRandSymExcitIt3(DetCurr,nExcit,nJ,Seed,IC,0,Prob,iCount)
!
!                    Child=AttemptCreate(DetCurr,CurrentSign(j),nJ,Prob,IC,Kik)
!!Kik is the off-diagonal hamiltonian matrix element for the walker. This is used for the augmentation of the death term if TDiffuse is on.
!                    IF(Child.gt.0) THEN
!!We have successfully created at least one positive child at nJ
!                        WSign=.true.
!                    ELSE
!!We have successfully created at least one negative child at nJ
!                        WSign=.false.
!                    ENDIF
!                    do l=1,abs(Child)
!                        do k=1,NEl
!                            NewDets(k,VecSlot)=nJ(k)
!                        enddo
!                        NewSign(VecSlot)=WSign
!                        VecSlot=VecSlot+1
!                    enddo
!
!                enddo
!
!            ELSE
!!Run through all possible excitations of each walker
!                
!                CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
!
!                do while(.true.)
!                    CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
!                    IF(nJ(1).eq.0) EXIT
!
!                    Child=AttemptCreate(DetCurr,CurrentSign(j),nJ,1.D0,IC,Kik)
!!Kik is the off-diagonal hamiltonian matrix element for the walker. This is used for the augmentation of the death term if TDiffuse is on.
!                    IF(Child.gt.0) THEN
!!We have successfully created at least one positive child at nJ
!                        WSign=.true.
!                    ELSE
!!We have successfully created at least one negative child at nJ
!                        WSign=.false.
!                    ENDIF
!                    do l=1,abs(Child)
!                        do k=1,NEl
!                            NewDets(k,VecSlot)=nJ(k)
!                        enddo
!                        NewSign(VecSlot)=WSign
!                        VecSlot=VecSlot+1
!                    enddo
!
!                enddo
!
!            ENDIF
!
!            KeepOrig=.true.
!            IF(TDiffuse) THEN
!!Next look at possibility of diffusion to another determinant
!                CALL GenRandSymExcitIt3(DetCurr,nExcit,nJ,Seed,IC,0,Prob,iCount)
!                CALL AttemptDiffuse(DetCurr,nJ,Prob,IC,CurrentSign(j),KeepOrig,CreateAtI,CreateAtJ)
!                !If we want to keep the original walker, then KeepOrig is true, However, we do not want to copy it accross yet, because we want to see if it is killed first in the birth/death process
!!                IF(KeepOrig) THEN
!!                    do k=1,NEl
!!                        NewDets(k,VecSlot)=DetCurr(k)
!!                    enddo
!!                    NewSign(VecSlot)=CurrentSign(j)
!!                    VecSlot=VecSlot+1
!!                ENDIF
!                do l=1,abs(CreateAtI)       !Sum in the number of walkers to create at the original determinant
!                    do k=1,NEl
!                        NewDets(k,VecSlot)=DetCurr(k)
!                    enddo
!                    IF(CreateAtI.gt.0) THEN
!                        NewSign(VecSlot)=.true.
!                    ELSE
!                        NewSign(VecSlot)=.false.
!                    ENDIF
!                    VecSlot=VecSlot+1
!                enddo
!                do l=1,abs(CreateAtJ)       !Add the number of walkers to create at nJ
!                    do k=1,NEl
!                        NewDets(k,VecSlot)=nJ(k)
!                    enddo
!                    IF(CreateAtJ.gt.0) THEN
!                        NewSign(VecSlot)=.true.
!                    ELSE
!                        NewSign(VecSlot)=.false.
!                    ENDIF
!                    VecSlot=VecSlot+1
!                enddo
!            ENDIF
!
!!We now have to decide whether the parent particle (j) wants to self-destruct or not...
!            iDie=AttemptDie(DetCurr,Kik,nExcitMemLen,nExcit,nStore)
!!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
!
!            IF(iDie.le.0) THEN
!!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle
!        
!                IF(KeepOrig) THEN
!                    ToCopy=abs(iDie)+1  !This is because we need to copy accross the original particle too
!                ELSE
!                    IF(iDie.eq.0) THEN
!                        ToCopy=0        !Indicates that want to spare original particle, which has already been copied accross previously, or previously been annihilated
!                    ELSE
!                        ToCopy=abs(iDie)
!                    ENDIF
!                ENDIF
!
!                do l=1,ToCopy    !We need to copy accross one more, since we need to include the original spared particle
!                    do k=1,NEl
!                        NewDets(k,VecSlot)=DetCurr(k)
!                    enddo
!                    NewSign(VecSlot)=CurrentSign(j)
!                    VecSlot=VecSlot+1
!                enddo
!
!            ELSEIF(iDie.gt.0) THEN
!!This indicates that particles on DetCurr want to be killed. The first kill will simply be performed by not copying accross the original particle.
!!Therefore, if iDie = 1, then we can simply ignore it.
!!However, after that anti-particles will need to be created on the same determinant.
!
!                IF(KeepOrig) iDie=iDie-1    
!!This is because we can already kill one particle by not copying accross the particle which was originally there (and still is even after diffusion). 
!!If KeepOrig is false, then the particle has already diffused somewhere else, and so antiparticles need to be created in its place.
!
!                do l=1,iDie
!                    do k=1,NEl
!                        NewDets(k,VecSlot)=DetCurr(k)
!                    enddo
!                    IF(CurrentSign(j)) THEN
!!Copy accross new anti-particles
!                        NewSign(VecSlot)=.FALSE.
!                    ELSE
!                        NewSign(VecSlot)=.TRUE.
!                    ENDIF
!                    VecSlot=VecSlot+1
!                enddo
!            
!            ENDIF
!
!!Destroy excitation generators for current walker
!            DEALLOCATE(nExcit)
!            CALL LogMemDealloc(this_routine,nExcitTag)
!        
!!            rat=(VecSlot+0.D0)/(MaxWalkers+0.D0)
!!            IF(rat.gt.0.9) THEN
!!                WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
!!            ENDIF
!
!!Finish cycling over walkers
!        enddo
!
!!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
!        TotWalkersNew=VecSlot-1
!        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
!        IF(rat.gt.0.9) THEN
!            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
!        ENDIF
!        
!        IF(TNodalCutoff.and.(NodalCutoff.lt.0.D0)) THEN
!!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff
!!then the net number of walkers must be the same sign, or it is set to zero walkers, and they are killed.
!            CALL TestWavevectorNodes(TotWalkersNew,2)
!        ENDIF
!
!        IF(TNoAnnihil) THEN
!!We are not annihilating particles - this will make things much quicker.
!
!!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
!            IF(associated(CurrentDets,target=WalkVecDets)) THEN
!                CurrentDets=>WalkVec2Dets
!                CurrentSign=>WalkVec2Sign
!                NewDets=>WalkVecDets
!                NewSign=>WalkVecSign
!            ELSE
!                CurrentDets=>WalkVecDets
!                CurrentSign=>WalkVecSign
!                NewDets=>WalkVec2Dets
!                NewSign=>WalkVec2Sign
!            ENDIF
!
!            TotWalkers=TotWalkersNew
!
!        ELSE
!!This routine now cancels down the particles with opposing sign on each determinant
!!This routine does not necessarily need to be called every Iter, but it does at the moment, since it is the only way to 
!!transfer the residual particles back onto CurrentDets
!            CALL AnnihilatePairs(TotWalkersNew)
!!            WRITE(6,*) "Number of annihilated particles= ",TotWalkersNew-TotWalkers
!        ENDIF
!
!
!        IF(TNodalCutoff.and.(NodalCutoff.ge.0.D0)) THEN
!!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff
!!then the net number of walkers must be the same sign, or it is set to zero walkers, and they are killed.
!            CALL TestWavevectorNodes(TotWalkers,1)
!        ENDIF
!
!        
!        IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!!Particle number is too large - kill them randomly
!
!!Log the fact that we have made a cull
!            NoCulls=NoCulls+1
!            IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!            CullInfo(NoCulls,1)=TotWalkers
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!
!            WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number..."
!            WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number..."
!            WRITE(6,"(A,F8.2)") "Population will reduce by a factor of ",CullFactor
!            CALL ThermostatParticles(.true.)
!
!!Need to reduce totwalkersOld, so that the change in shift is also reflected by this
!!The Shift is no longer calculated like this...
!!            TotWalkersOld=nint((TotWalkersOld+0.D0)/CullFactor)
!
!        ELSEIF((TotWalkers.lt.(InitWalkers/2)).and.(.NOT.TNoBirth)) THEN
!!Particle number is too small - double every particle in its current position
!
!!Log the fact that we have made a cull
!            NoCulls=NoCulls+1
!            IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!            CullInfo(NoCulls,1)=TotWalkers
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!            
!            WRITE(6,*) "Doubling particle population to increase total number..."
!            CALL ThermostatParticles(.false.)
!
!!Need to increase TotWalkersOld, so that the change in shift is also reflected by this
!!The shift is no longer calculated like this...
!!            TotWalkersOld=TotWalkersOld*2
!
!        ENDIF
!
!        CALL TIHALT('MCyc',iSubCyc)
!
!        RETURN
!
!    END SUBROUTINE PerformFCIMCyc
!
!!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff,
!!then the walkers at that determinant must be of the same sign, or they are killed.
!!Particles indicates the number of particles to look through, and iArray is 1 if the WalkVecDets array is active, and 2 if the WalkVec2Dets array is active.
!    SUBROUTINE TestWavevectorNodes(Particles,iArray)
!        USE Calc , only : i_P
!        USE System , only : Beta
!        USE Integrals , only : nTay
!        IMPLICIT NONE
!        INTEGER :: ParticlesOrig,VecSlot,IC,iGetExcitLevel,j,k,NoatHF,Particles,iArray,NoPositive,NoNegative
!        INTEGER , POINTER :: ActiveVecDets(:,:)
!        LOGICAL , POINTER :: ActiveVecSign(:)
!        REAL*8 :: EigenComp,EnergyNum
!        LOGICAL :: Component
!        TYPE(HElement) :: Hamij,Fj!,rhjj,rhij
!
!!We first need to point to the active array
!        IF(iArray.eq.1) THEN
!            ActiveVecDets => CurrentDets
!            ActiveVecSign => CurrentSign
!        ELSEIF(iArray.eq.2) THEN
!            ActiveVecDets => NewDets
!            ActiveVecSign => NewSign
!        ELSE
!            CALL Stop_All("TestWavevectorNodes","Error with iArray")
!        ENDIF
!
!        EnergyNum=0.D0  !EnergyNum indicates the sum over all the hamiltonian matrix elements between the double excitations and HF
!        NoatHF=0
!        NoPositive=0    !Total number of positive particles
!        NoNegative=0    !Total number of negative particles
!
!        ParticlesOrig=Particles
!!VecSlot indicates the next free position in ActiveVec
!        VecSlot=1
!
!        do j=1,Particles
!!j runs through all current walkers
!!IC labels the excitation level away from the HF determinant
!            IC=iGetExcitLevel(FDet,ActiveVecDets(:,j),NEl)
!
!            IF((IC.gt.2).or.(IC.eq.1)) THEN
!!More than double excitations/singles are not included in the MP1 wavefunction, therefore we can copy these straight across
!                IF(VecSlot.lt.j) THEN
!!Only copy accross to VecSlot if VecSlot<j, otherwise, if VecSlot=j, particle already in right place.
!                    do k=1,NEl
!                        ActiveVecDets(k,VecSlot)=ActiveVecDets(k,j)
!                    enddo
!                    ActiveVecSign(VecSlot)=ActiveVecSign(j)
!                ENDIF
!                VecSlot=VecSlot+1
!
!                IF(ActiveVecSign(j)) THEN
!                    NoPositive=NoPositive+1
!                ELSE
!                    NoNegative=NoNegative+1
!                ENDIF
!
!            ELSEIF((IC.eq.2).or.(IC.eq.0)) THEN
!                IF(IC.eq.2) THEN
!!We are at a double excitation - first find the desired sign of the determinant.
!                    Hamij=GetHElement2(FDet,ActiveVecDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!                    IF((real(Hamij%v)).lt.0) THEN
!!Negative Hamiltonian connection indicates that the component of the MP1 wavevector is positive
!                        Component=.true.
!                    ELSE
!                        Component=.false.
!                    ENDIF
!                ELSE
!!We are at the HF determinant
!                    Hamij=Hii
!                    Component=.true.    !HF component of MP1 wavevector always wants to be positive
!                ENDIF
!
!                IF((ActiveVecSign(j).and.Component).or.((.NOT.ActiveVecSign(j)).and.(.NOT.Component))) THEN
!!The particle is of the same sign as the MP1 wavevector component, so whether or not the determinant is of fixed sign, the particle does not want to be destroyed
!                    IF(VecSlot.lt.j) THEN
!!Only copy accross to VecSlot if VecSlot<j, otherwise, if VecSlot=j, particle already in right place.
!                        do k=1,NEl
!                            ActiveVecDets(k,VecSlot)=ActiveVecDets(k,j)
!                        enddo
!                        ActiveVecSign(VecSlot)=ActiveVecSign(j)
!                    ENDIF
!                    VecSlot=VecSlot+1
!
!!Add to the estimate for the energy if we want to keep the particle
!                    IF(ActiveVecSign(j)) THEN
!                        IF(Iter.gt.NEquilSteps) EnergyNum=EnergyNum+(REAL(Hamij%v,r2))
!                        NoPositive=NoPositive+1
!                    ELSE
!                        IF(Iter.gt.NEquilSteps) EnergyNum=EnergyNum-(REAL(Hamij%v,r2))
!                        NoNegative=NoNegative+1
!                    ENDIF
!                    IF(IC.eq.0) THEN
!                        IF(ActiveVecSign(j)) THEN
!                            IF(Iter.gt.NEquilSteps) NoatHF=NoatHF+1
!                        ELSE
!                            IF(Iter.gt.NEquilSteps) NoatHF=NoatHF-1
!                        ENDIF
!                    ENDIF
!
!                ELSE
!!Particle is of a different sign to the component of the MP1 wavefunction - if the determinant is of fixed sign, we don't want to copy it across. Find if fixed sign...
!                    IF(IC.eq.2) THEN
!                        CALL GetH0Element(ActiveVecDets(:,j),NEl,Arr,nBasis,ECore,Fj)
!!EigenComp is now the component of the normalised MP1 wavefunction for the double excitation WalkVecDets(:,j)
!                        EigenComp=abs(((Hamij%v)/(Fj%v-FZero%v))/MPNorm)
!                    ELSE
!!Here, EigenComp is the component of the normalised MP1 wavefunction at the HF determinant (should always be +ve)
!                        EigenComp=1.D0/MPNorm
!                    ENDIF
!
!                    IF(EigenComp.lt.abs(NodalCutoff)) THEN
!!Determinant does not have a fixed sign - keep the particle
!                        IF(VecSlot.lt.j) THEN
!!Only copy accross to VecSlot if VecSlot<j, otherwise, if VecSlot=j, particle already in right place.
!                            do k=1,NEl
!                                ActiveVecDets(k,VecSlot)=ActiveVecDets(k,j)
!                            enddo
!                            ActiveVecSign(VecSlot)=ActiveVecSign(j)
!                        ENDIF
!                        VecSlot=VecSlot+1
!
!!Add to the estimate for the energy if we want to keep the particle
!                        IF(ActiveVecSign(j)) THEN
!                            IF(Iter.gt.NEquilSteps) EnergyNum=EnergyNum+(REAL(Hamij%v,r2))
!                            NoPositive=NoPositive+1
!                        ELSE
!                            IF(Iter.gt.NEquilSteps) EnergyNum=EnergyNum-(REAL(Hamij%v,r2))
!                            NoNegative=NoNegative+1
!                        ENDIF
!                        IF(IC.eq.0) THEN
!                            IF(ActiveVecSign(j)) THEN
!                                IF(Iter.gt.NEquilSteps) NoatHF=NoatHF+1
!                            ELSE
!                                IF(Iter.gt.NEquilSteps) NoatHF=NoatHF-1
!                            ENDIF
!                        ENDIF
!
!                    ENDIF
!                ENDIF
!            ELSE
!                CALL Stop_All("TestWavevectorNodes","Should not be here - wrong IC calculated")
!            ENDIF
!        enddo   !end loop over all walkers
!
!!Total number of particles now modified by killing of wrong signed particles
!        Particles=VecSlot-1
!    
!        IF((ParticlesOrig-Particles).ne.0) THEN
!            WRITE(6,*) "Some particles killed by Nodal approximation: ", ParticlesOrig-Particles
!        ENDIF
!
!!Calculate the time average of the numerator and demonimator to calculate the running average of the energy - this will then not be affected by the case when there aren't any particles at the HF determinant
!        SumNoatHF=SumNoatHF+NoatHF
!        SumENum=SumENum+EnergyNum
!        ProjectionE=(SumENum/(SumNoatHF+0.D0))-REAL(Hii%v,r2)
!
!        IF(NoatHF.ne.0) THEN
!!The energy cannot be calculated via the projection back onto the HF if there are no particles at HF
!            SumE=SumE+((EnergyNum/(NoatHF+0.D0))-(REAL(Hii%v,r2)))
!            ProjectionEInst=SumE/((Iter-CycwNoHF)+0.D0)
!        ELSE
!            CycwNoHF=CycwNoHF+1         !Record the fact that there are no particles at HF in this run, so we do not bias the average
!!            WRITE(6,*) "No positive particles at reference determinant during iteration: ",Iter
!        ENDIF
!
!        NetPositive=NoPositive-NoNegative
!!        WRITE(6,*) ParticlesOrig,Particles,NetPositive,NoPositive,NoNegative
!
!        RETURN
!
!    END SUBROUTINE TestWavevectorNodes
!
!!This routine acts as a thermostat for the simulation - killing random particles if the population becomes too large, or 
!!Doubling them if it gets too low...
!    SUBROUTINE ThermostatParticles(HighLow)
!        IMPLICIT NONE
!        LOGICAL :: HighLow
!        INTEGER :: VecSlot,i,j,ToCull,Culled,OrigWalkers,Chosen
!        REAL*8 :: Ran2
!
!        IF(HighLow) THEN
!!The population is too large - cull TotWalkers/CullFactor randomly selected particles
!
!            OrigWalkers=TotWalkers
!            ToCull=TotWalkers-nint((TotWalkers+0.D0)/CullFactor)
!            Culled=0
!
!            do while (Culled.lt.ToCull)
!
!!Pick a random walker between 1 and TotWalkers
!                Chosen=int((Ran2(Seed)*TotWalkers)+1.D0)
!
!!Move the Walker at the end of the list to the position of the walker we have chosen to destroy
!                do i=1,NEl
!                    CurrentDets(i,Chosen)=CurrentDets(i,TotWalkers)
!                enddo
!                CurrentSign(Chosen)=CurrentSign(TotWalkers)
!
!                TotWalkers=TotWalkers-1
!                Culled=Culled+1
!
!            enddo
!
!            IF(TotWalkers.ne.(OrigWalkers-ToCull)) THEN
!                WRITE(6,*) "Error in culling walkers..."
!                STOP "Error in culling walkers..."
!            ENDIF
!
!!CullInfo(:,2) is the new number of total walkers
!            CullInfo(NoCulls,2)=TotWalkers
!
!        ELSE
!!The population is too low - give it a boost by doubling every particle
!
!            VecSlot=TotWalkers+1
!            do i=1,TotWalkers
!
!!Add clone of walker, at the same determinant, to the end of the list
!                do j=1,NEl
!                    CurrentDets(j,VecSlot)=CurrentDets(j,i)
!                enddo
!                CurrentSign(VecSlot)=CurrentSign(i)
!
!                VecSlot=VecSlot+1
!
!            enddo
!
!            TotWalkers=TotWalkers*2
!
!            IF((VecSlot-1).ne.TotWalkers) THEN
!                WRITE(6,*) "Problem in doubling all particles..."
!                STOP "Problem in doubling all particles..."
!            ENDIF
!
!!CullInfo(:,2) is the new number of total walkers
!            CullInfo(NoCulls,2)=TotWalkers
!
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE ThermostatParticles
!
!
!!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!!value of the diagonal shift in the hamiltonian in order to compensate for this
!    SUBROUTINE UpdateDiagSft()
!        IMPLICIT NONE
!        INTEGER :: j,k,GrowthSteps
!
!        IF(NoCulls.eq.0) THEN
!            GrowRate=(TotWalkers+0.D0)/(TotWalkersOld+0.D0)
!        ELSEIF(NoCulls.eq.1) THEN
!!GrowRate is the sum of the individual grow rates for each uninterrupted growth sequence, multiplied by the fraction of the cycle which was spent on it
!            GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((CullInfo(1,1)+0.D0)/(TotWalkersOld+0.D0))
!            GrowRate=GrowRate+(((StepsSft-CullInfo(1,3))+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(1,2)+0.D0))
!
!            NoCulls=0
!            CALL IAZZERO(CullInfo,30)
!        ELSE
!            GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((CullInfo(1,1)+0.D0)/(TotWalkersOld+0.D0))
!            do j=2,NoCulls
!    
!!This is needed since the steps between culling is stored cumulatively
!                GrowthSteps=CullInfo(j,3)-CullInfo(j-1,3)
!                GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((CullInfo(j,1)+0.D0)/(CullInfo(j-1,2)+0.D0))
!
!            enddo
!
!            GrowthSteps=StepsSft-CullInfo(NoCulls,3)
!            GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(NoCulls,2)+0.D0))
!
!            NoCulls=0
!            CALL IAZZERO(CullInfo,30)
!
!        ENDIF
!        DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
!!        IF((DiagSft).gt.0.D0) THEN
!!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!!            STOP
!!        ENDIF
!
!    END SUBROUTINE UpdateDiagSft
!
!
!!This routine cancels out particles of opposing sign on the same determinant.
!    SUBROUTINE AnnihilatePairs(TotWalkersNew)
!        IMPLICIT NONE
!        INTEGER :: TotWalkersNew,j,k,l,DetCurr(NEl),VecSlot,TotWalkersDet
!        INTEGER :: DetLT
!
!!First, it is necessary to sort the list of determinants
!        CALL SortDets(TotWalkersNew,NewDets(:,1:TotWalkersNew),NEl,NewSign(1:TotWalkersNew),1)
!
!!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to CurrentDets
!        j=1
!!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
!        do k=1,NEl
!!DetCurr is the current determinant
!            DetCurr(k)=NewDets(k,j)
!        enddo
!        VecSlot=1
!
!        do while(j.le.TotWalkersNew)
!!Loop over all walkers
!            TotWalkersDet=0
!            do while ((DetLT(NewDets(:,j),DetCurr,NEl).eq.0).and.(j.le.TotWalkersNew))
!!Loop over all walkers on DetCurr and count residual number after cancelling
!                IF(NewSign(j)) THEN
!                    TotWalkersDet=TotWalkersDet+1
!                ELSE
!                    TotWalkersDet=TotWalkersDet-1
!                ENDIF
!                j=j+1
!            enddo
!!Transfer residual population into VecSlot, along with residual sign
!            IF(TotWalkersDet.gt.0) THEN
!!Positive sign particles want to populate this determinant
!                do l=1,abs(TotWalkersDet)
!                    do k=1,NEl
!                        CurrentDets(k,VecSlot)=DetCurr(k)
!                    enddo
!                    CurrentSign(VecSlot)=.true.
!                    VecSlot=VecSlot+1
!                enddo
!            ELSE
!!Negative sign particles want to populate this determinant
!                do l=1,abs(TotWalkersDet)
!                    do k=1,NEl
!                        CurrentDets(k,VecSlot)=DetCurr(k)
!                    enddo
!                    CurrentSign(VecSlot)=.false.
!                    VecSlot=VecSlot+1
!                enddo
!            ENDIF
!!Now update the current determinant
!            do k=1,NEl
!                DetCurr(k)=NewDets(k,j)
!            enddo
!        enddo
!!The new number of residual cancelled walkers is given by one less that VecSlot again.
!        TotWalkers=VecSlot-1
!
!        RETURN
!
!    END SUBROUTINE AnnihilatePairs
!
!!This routine calculates the MP1 eigenvector, and uses it as a guide for setting the initial walker configuration
!    SUBROUTINE StartWavevector(WaveType)
!        USE Calc , only : i_P
!        USE System , only : Beta
!        USE Integrals , only : nTay
!        IMPLICIT NONE
!        INTEGER :: ierr,i,j,WaveType,EigenvectorTag=0,k,VecSlot,NoDoublesWalk
!        CHARACTER(len=*), PARAMETER :: this_routine='StartWavevector'
!        REAL*8 :: TypeChange,SumComp,GrowFactor
!        INTEGER :: nStore(6),nExcitMemLen,nJ(NEl),iMaxExcit,nExcitTag=0,iExcit,WalkersOnDet
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        REAL*8 , ALLOCATABLE :: Eigenvector(:)
!        TYPE(HElement) :: rhij,rhjj,Hamij,Fj,MP2E
!
!        IF((WaveType.eq.1).or.(WaveType.eq.2)) THEN
!!If WaveType=1, we want to calculate the MP1 wavevector as our initial configuration, and WaveType=2 is the star wavevector
!        
!            IF(WaveType.eq.1) THEN
!!For MP1 wavefunction, we want TypeChange=1.D0, but star will require the star correlation energy
!                TypeChange=1.D0
!            ELSEIF(WaveType.eq.2) THEN
!!Star energy not yet proparly coded & tested
!                STOP 'Star initial wavefunction not yet working'
!!                StarWeight=fMCPR3StarNewExcit(FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,nTay,RhoEps,iExcit,iMaxExcit,nWHTay,iLogging,TSym,ECore,DBeta,DLWDB,MP2E)
!!                TypeChange=DLWDB
!            ENDIF
!        
!            MP2E=0.D0
!
!!First, generate all excitations, and store their determianants, and rho matrix elements 
!            nStore(1)=0
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!            nExcit(1)=0
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!
!            ALLOCATE(Eigenvector(iMaxExcit+1),stat=ierr)
!            CALL LogMemAlloc('Eigenvector',iMaxExcit+1,8,this_routine,EigenvectorTag,ierr)
!            CALL AZZERO(Eigenvector,iMaxExcit+1)
!
!!Also need to store the determinants which each component of the eigenvector refers to...
!            ALLOCATE(ExcitStore(NEl,iMaxExcit+1),stat=ierr)
!            CALL LogMemAlloc('ExcitStore',(iMaxExcit+1)*NEl,4,this_routine,ExcitStoreTag,ierr)
!            CALL IAZZERO(ExcitStore,(iMaxExcit+1)*NEl)
!            
!!            CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
!            CALL GetH0Element(FDet,NEl,Arr,nBasis,ECore,FZero)
!
!            i=1
!            do j=1,NEl
!                ExcitStore(j,i)=FDet(j)
!            enddo
!            Eigenvector(i)=1.D0
!            SumComp=1.D0
!            
!            do while (.true.)
!                CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
!                IF(nJ(1).eq.0) EXIT
!                i=i+1
!!                CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
!!                CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)
!
!!We want the value of rho_jj/rho_ii
!!                rhjj=rhjj/rhii
!                Hamij=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
!                CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fj)
!
!                do j=1,NEl
!                    ExcitStore(j,i)=nJ(j)
!                enddo
!!                Eigenvector(i)=(rhij%v)/((rhjj%v)-TypeChange)
!                MP2E=MP2E%v-((Hamij%v)**2)/((Fj%v)-(FZero%v))
!                Eigenvector(i)=(Hamij%v)/(Fj%v-FZero%v)
!                SumComp=SumComp+Eigenvector(i)
!            enddo
!
!            NoComps=i
!
!            WRITE(6,*) "Number of components found for new starting wavevector: ",NoComps
!
!            DEALLOCATE(nExcit)
!            CALL LogMemDealloc(this_routine,nExcitTag)
!
!            DiagSft=real(MP2E%v,KIND(0.D0))
!
!        ENDIF
!
!        GrowFactor=(InitWalkers+0.D0)/SumComp
!
!        WRITE(6,*) "Growth factor of initial wavevector is: ",GrowFactor
!
!!Find actual number of new initial walkers
!        NoDoublesWalk=0
!        InitWalkers=0
!        do i=1,NoComps
!            WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
!            InitWalkers=InitWalkers+WalkersOnDet
!            IF(i.ne.1) THEN
!                NoDoublesWalk=NoDoublesWalk+WalkersOnDet
!            ENDIF
!        enddo
!
!        WRITE(6,*) "New number of initial walkers is: ",InitWalkers
!        WRITE(6,*) "Number of walkers on double excitations: ",NoDoublesWalk
!
!!Set the maximum number of walkers allowed
!        MaxWalkers=MemoryFac*InitWalkers
!
!!Allocate memory to hold walkers
!        ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
!        CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
!        ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
!        CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
!        ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
!        CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
!        ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
!        CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)
!    
!        CurrentDets=>WalkVecDets
!        CurrentSign=>WalkVecSign
!        NewDets=>WalkVec2Dets
!        NewSign=>WalkVec2Sign
!
!        VecSlot=1
!!Cycle over components
!        do i=1,NoComps
!!Cycle over number of initial walkers wanted in each component
!            IF(Eigenvector(i).gt.0.D0) THEN
!!Walkers in this determinant are positive
!                WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
!                do j=1,WalkersOnDet
!                    do k=1,NEl
!                        CurrentDets(k,VecSlot)=ExcitStore(k,i)
!                    enddo
!                    CurrentSign(VecSlot)=.true.
!                    VecSlot=VecSlot+1
!                enddo
!
!            ELSE
!                WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
!                do j=1,WalkersOnDet
!                    do k=1,NEl
!                        CurrentDets(k,VecSlot)=ExcitStore(k,i)
!                    enddo
!                    CurrentSign(VecSlot)=.false.
!                    VecSlot=VecSlot+1
!                enddo
!            ENDIF
!
!        enddo
!
!        IF((VecSlot-1).ne.InitWalkers) THEN
!            WRITE(6,*) "Problem in assigning particles proportionally to given wavevector..."
!            STOP "Problem in assigning particles proportionally to given wavevector..."
!        ENDIF
!
!        IF(TNoBirth) THEN
!!Allocate memory to hold all initial determinants and their starting populations
!            ALLOCATE(InitPops(NoComps),stat=ierr)
!            CALL LogMemAlloc('InitPops',NoComps,4,this_routine,InitPopsTag)
!                
!            WRITE(6,*) ""
!            WRITE(6,*) "    Step   Components   WalkersOnDet   Expected"
!            WRITE(15,*) "#    Step   Components   WalkersOnDet   Expected"
!
!            do i=1,NoComps
!                InitPops(i)=abs(nint(Eigenvector(i)*GrowFactor))
!                WRITE(6,"(2I9,2G16.7)") 0,i,InitPops(i),InitPops(i)
!                WRITE(15,"(2I9,2G16.7)") 0,1,InitPops(i),InitPops(i)
!            enddo
!            WRITE(6,*) ""
!            WRITE(15,*) ""
!
!!Can deallocate the eigenvector, but we want to keep the ExcitStore in order to compare the populations at a later date
!            DEALLOCATE(Eigenvector)
!            CALL LogMemDealloc(this_routine,EigenvectorTag)
!
!        ELSE
!
!            DEALLOCATE(Eigenvector)
!            CALL LogMemDealloc(this_routine,EigenvectorTag)
!            DEALLOCATE(ExcitStore)
!            CALL LogMemDealloc(this_routine,ExcitStoreTag)
!
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE StartWavevector
!
!
!!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
!    SUBROUTINE InitFCIMCCalc()
!        USE DetCalc , only : NDet
!        IMPLICIT NONE
!        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
!        INTEGER :: DetLT,VecSlot
!        TYPE(HElement) :: rh
!        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMC'
!
!
!        IF(TStartMP1) THEN
!!Start the initial distribution off at the distribution of the MP1 eigenvector
!
!            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."
!            CALL StartWavevector(1)
!
!        ELSEIF(TReadPops) THEN
!
!!Set the maximum number of walkers allowed
!            MaxWalkers=MemoryFac*InitWalkers
!
!!Allocate memory to hold walkers
!            ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
!            CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
!            ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
!            CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
!            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
!            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
!            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
!            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)
!
!            CurrentDets=>WalkVecDets
!            CurrentSign=>WalkVecSign
!            NewDets=>WalkVec2Dets
!            NewSign=>WalkVec2Sign
!
!            IF((ABS(ScaleWalkers-1.D0)).lt.1.D-8) THEN
!!Read in walker positions
!                do i=1,InitWalkers
!                    READ(17,*) CurrentDets(:,i),CurrentSign(i)
!                enddo
!            ELSE
!!Read in walker positions - we will scale these later...
!                do i=1,InitWalkers
!                    READ(17,*) NewDets(:,i),NewSign(i)
!                enddo
!                WRITE(6,*) "Scaling number of walkers by: ",ScaleWalkers
!                ReadWalkers=InitWalkers
!                InitWalkers=0
!!First, count the total number of initial walkers on each determinant - sort into list
!                CALL SortDets(ReadWalkers,NewDets(:,1:ReadWalkers),NEl,NewSign(1:ReadWalkers),1)
!
!                j=1
!!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
!                do k=1,NEl
!!DetCurr is the current determinant
!                    DetCurr(k)=NewDets(k,j)
!                enddo
!
!                do while(j.le.ReadWalkers)
!!Loop over all walkers
!                    TotWalkersDet=0
!                    do while ((DetLT(NewDets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!!Loop over all walkers on DetCurr and count residual number after cancelling
!                        IF(NewSign(j)) THEN
!                            TotWalkersDet=TotWalkersDet+1
!                        ELSE
!                            TotWalkersDet=TotWalkersDet-1
!                        ENDIF
!                        j=j+1
!                    enddo
!!Now update the current determinant
!                    do k=1,NEl
!                        DetCurr(k)=NewDets(k,j)
!                    enddo
!!Count total number of initial walkers
!                    InitWalkers=InitWalkers+abs(nint((TotWalkersDet+0.D0)*ScaleWalkers))
!                enddo
!                WRITE(6,*) "Total number of walkers is now: ",InitWalkers
!!Set the new maximum number of walkers allowed
!                MaxWalkers=MemoryFac*InitWalkers
!
!!Deallocate old memory block for WalkVec
!                DEALLOCATE(WalkVecDets)
!                CALL LogMemDealloc(this_routine,WalkVecDetsTag)
!                DEALLOCATE(WalkVecSign)
!                CALL LogMemDealloc(this_routine,WalkVecSignTag)
!
!!Allocate memory to hold new maximum number of walkers
!                ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
!                CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
!                ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
!                CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
!                
!                CurrentDets=>WalkVecDets
!                CurrentSign=>WalkVecSign
!
!!Now multiply them up...
!                j=1
!                VecSlot=1
!!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
!                do k=1,NEl
!                    DetCurr(k)=NewDets(k,j)
!                enddo
!!DetCurr is the current determinant
!                do while(j.le.ReadWalkers)
!!Loop over all walkers
!                    TotWalkersDet=0
!                    do while ((DetLT(NewDets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!!Loop over all walkers on DetCurr and count residual number after cancelling
!                        IF(NewSign(j)) THEN
!                            TotWalkersDet=TotWalkersDet+1
!                        ELSE
!                            TotWalkersDet=TotWalkersDet-1
!                        ENDIF
!                        j=j+1
!                    enddo
!!Now multiply up the number of walkers, and insert into CurrentDets
!                    TotWalkersDet=nint((TotWalkersDet+0.D0)*ScaleWalkers)
!                    IF(TotWalkersDet.gt.0) THEN
!                        do l=1,abs(TotWalkersDet)
!                            do k=1,NEl
!                                CurrentDets(k,VecSlot)=DetCurr(k)
!                            enddo
!                            CurrentSign(VecSlot)=.true.
!                            VecSlot=VecSlot+1
!                        enddo
!                    ELSE
!                        do l=1,abs(TotWalkersDet)
!                            do k=1,NEl
!                                CurrentDets(k,VecSlot)=DetCurr(k)
!                            enddo
!                            CurrentSign(VecSlot)=.false.
!                            VecSlot=VecSlot+1
!                        enddo
!                    ENDIF
!                    do k=1,NEl
!                        DetCurr(k)=NewDets(k,j)
!                    enddo
!                enddo
!                IF((VecSlot-1).ne.InitWalkers) THEN
!                    WRITE(6,*) "Problem scaling up walker number - exiting..."
!                    STOP 'Problem scaling up walker number - exiting...'
!                ENDIF
!
!!Now deallocate and reallocate WalkVec2 with correct number of total walkers
!                DEALLOCATE(WalkVec2Dets)
!                CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
!                DEALLOCATE(WalkVec2Sign)
!                CALL LogMemDealloc(this_routine,WalkVec2SignTag)
!                ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
!                CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
!                ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
!                CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)
!
!                NewDets=>WalkVec2Dets
!                NewSign=>WalkVec2Sign
!
!            ENDIF
!
!!End of reading in POPSFILE
!            CLOSE(17)
!
!        ELSE
!!If not reading in from POPSFILE, then we need to initialise the particle positions - start at HF with positive sign
!
!!Set the maximum number of walkers allowed
!            IF(TMCDiffusion) THEN
!                MaxWalkers=InitWalkers
!            ELSE
!                MaxWalkers=MemoryFac*InitWalkers
!            ENDIF
!
!!Allocate memory to hold walkers
!            ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
!            CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
!            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
!            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
!
!            CurrentDets=>WalkVecDets
!            CurrentSign=>WalkVecSign
!
!            IF(.NOT.TMCDiffusion) THEN
!                ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
!                CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)
!                ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
!                CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
!                NewDets=>WalkVec2Dets
!                NewSign=>WalkVec2Sign
!            ENDIF
!
!            do j=1,InitWalkers
!                do k=1,NEl
!                    CurrentDets(k,j)=FDet(k)
!                enddo
!                CurrentSign(j)=.true.
!            enddo
!
!            IF(TNoBirth) THEN
!                
!                WRITE(6,*) ""
!                WRITE(6,*) "    Step   Components   WalkersOnDet   Expected"
!                WRITE(15,*) "#    Step   Components   WalkersOnDet   Expected"
!
!                ALLOCATE(InitPops(1),stat=ierr)
!                CALL LogMemAlloc('InitPops',1,4,this_routine,InitPopsTag)
!                ALLOCATE(ExcitStore(NEl,1),stat=ierr)
!                CALL LogMemAlloc('ExcitStore',NEl,4,this_routine,ExcitStoreTag)
!
!                NoComps=1
!                do j=1,NEl
!                    ExcitStore(j,1)=FDet(j)
!                enddo
!                InitPops(1)=InitWalkers
!                WRITE(6,"(2I9,2G16.7)") 0,1,InitWalkers,InitWalkers
!                WRITE(15,"(2I9,2G16.7)") 0,1,InitWalkers,InitWalkers
!
!            ELSEIF(TDetPops) THEN
!                
!                SizeofSpace=NDET
!                WRITE(6,*) "Size of space is: ", SizeOfSpace
!
!                ALLOCATE(PopsVec(SizeofSpace),stat=ierr)
!                CALL LogMemAlloc('PopsVec',SizeofSpace,8,this_routine,PopsVecTag)
!                CALL AZZERO(PopsVec,SizeofSpace)
!                ALLOCATE(TransMat(SizeofSpace,SizeofSpace),stat=ierr)
!                CALL LogMemAlloc('TransMat',SizeofSpace**2,8,this_routine,TransMatTag)
!                CALL AZZERO(TransMat,SizeofSpace**2)
!
!                IF(DetLT(NMRKS(:,1),FDet,NEl).ne.0) THEN
!                    WRITE(6,*) "Problem with NMRKS"
!                    STOP "Problem with NMRKS"
!                ENDIF
!                
!                do j=1,SizeOfSpace
!                    rh=GetHElement2(FDet,NMRKS(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
!                    WRITE(67,"(F10.5,$)") rh%v
!                enddo
!                WRITE(67,*) ""
!                WRITE(67,*) ""
!
!            ENDIF
!
!        ENDIF
!
!!TotWalkers contains the number of current walkers at each step
!        TotWalkers=InitWalkers
!        TotWalkersOld=InitWalkers
!
!        CALL IAZZERO(CullInfo,30)
!        NoCulls=0
!
!        IF(TNodalCutoff.and.(.not.TMCDiffusion)) CALL CalcNodalSurface()
!
!        IF(TResumFciMC) THEN
!!Allocate memory to hold graphs and corresponding vectors for ResumFciMC.
!            ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
!            CALL LogMemAlloc('GraphRhoMat',NDets**2,8,this_routine,GraphRhoMatTag)
!            ALLOCATE(GraphVec(NDets),stat=ierr)
!            CALL LogMemAlloc('GraphVec',NDets,8,this_routine,GraphVecTag)
!            ALLOCATE(DetsinGraph(NEl,NDets),stat=ierr)
!            CALL LogMemAlloc('DetsinGraph',NEl*NDets,4,this_routine,DetsinGraphTag)
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE InitFCIMCCalc
!
!!This routine calculates the normalisation for the MP1 wavefunction. This is needed if a nodal structure is being applied, and can also calculate the number of determinants which are being constrained.
!    SUBROUTINE CalcNodalSurface()
!        USE Calc , only : i_P
!        USE System , only : Beta
!        USE Integrals , only : nTay
!        IMPLICIT NONE
!        INTEGER :: ierr,i,j,k,Doubs,FixedSign
!        REAL*8 :: EigenComp
!        CHARACTER(len=*), PARAMETER :: this_routine='CalcNodalSurface'
!        INTEGER :: nStore(6),nExcitMemLen,nJ(NEl),iMaxExcit,nExcitTag=0,iExcit
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        TYPE(HElement) :: Fj,Hamij!,rhij,rhjj
!
!        WRITE(6,"(A,F19.9)") "Calculating the nodal structure of the MP1 wavefunction with a normalised cutoff of ",NodalCutoff
!
!!        IF(nTay(2).ne.3) THEN
!!This is no longer needed since the MP1 components are calculated exactly
!!Fock-partition-lowdiag is not set - it must be in order to use the given formulation for the MP1 wavefunction
!!            WRITE(6,*) "FOCK-PARTITION-LOWDIAG is not specified. It must be to use NODALCUTOFF."
!!            WRITE(6,*) "Resetting all rho integrals to use FOCK-PARTITION-LOWDIAG"
!!            nTay(2)=3
!!        ENDIF
!
!!First, generate all excitations, and store their determianants, and rho matrix elements 
!        nStore(1)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,2)
!        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!        nExcit(1)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,2)
!
!!Calculate F0
!        CALL GetH0Element(FDet,NEl,Arr,nBasis,ECore,FZero)
!!        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
!
!        Doubs=0
!!The HF determinant has a component 1         
!        MPNorm=1.D0
!        do while (.true.)
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,2)
!            IF(nJ(1).eq.0) EXIT
!            Doubs=Doubs+1
!            Hamij=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
!            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fj)
!!            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
!!            CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)
!
!!We want the value of rho_jj/rho_ii
!!            rhjj=rhjj/rhii
!!EigenComp is now the component of the MP1 wavefunction for nJ, but unnormalised
!!            EigenComp=(rhij%v)/((rhjj%v)-1.D0)
!            EigenComp=(Hamij%v)/(Fj%v-FZero%v)
!            MPNorm=MPNorm+(EigenComp**2)
!        enddo
!
!        WRITE(6,"(A,I15)") "Number of double excitations found in MP Wavevector: ",Doubs
!
!!Find the normalisation for the MP1 wavevector
!        MPNorm=SQRT(MPNorm)
!
!!Reset the excitation generator, so that we can run through the excitations again and calculate the number of components of the MP1 wavefunction which will be sign-constrained
!        CALL ResetExit2(FDet,NEl,G1,nBasis,nBasisMax,nExcit,0)
!
!        FixedSign=0     !FixedSign is the number of determinants which are fixed in sign with the given value of the NodalCutoff
!
!        IF((1.D0/MPNorm).gt.abs(NodalCutoff)) THEN
!!This indicates that the HF is included in the fixed sign approximation, and so will always be constrained to have net positive particles.
!            FixedSign=FixedSign+1
!            WRITE(6,*) "Hartree-Fock determinant constrained to always have net positive sign"
!        ELSE
!            WRITE(6,"(A,F17.9)") "Hartree-Fock determinant NOT constrained to have net positive sign, since normalised component is only: ", 1.D0/MPNorm
!        ENDIF
!
!        do while(.true.)
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,2)
!            IF(nJ(1).eq.0) EXIT
!            Hamij=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
!            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fj)
!!            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
!!            CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)
!
!!We want the value of rho_jj/rho_ii
!!            rhjj=rhjj/rhii
!!EigenComp is now the component of the MP1 wavefunction for nJ, but now normalised
!!            EigenComp=abs(((rhij%v)/((rhjj%v)-1.D0))/MPNorm)
!            EigenComp=abs(((Hamij%v)/(Fj%v-FZero%v))/MPNorm)
!            IF(EigenComp.gt.abs(NodalCutoff)) THEN
!!Increase the counter of number of determinants with fixed sign
!                FixedSign=FixedSign+1
!            ENDIF
!
!        enddo
!
!        WRITE(6,*) "Total number of determinants constrained by fixed sign: ", FixedSign
!
!        DEALLOCATE(nExcit)
!        CALL LogMemDealloc(this_routine,nExcitTag)
!
!        RETURN
!
!    END SUBROUTINE CalcNodalSurface
!
!!This function tells us whether we want to diffuse from DetCurr to nJ
!    SUBROUTINE AttemptDiffuse(DetCurr,nJ,Prob,IC,WSign,KeepOrig,CreateAtI,CreateAtJ)
!        IMPLICIT NONE
!        INTEGER :: DetCurr(NEl),nJ(NEl),IC,i,CreateAtI,CreateAtJ
!        REAL*8 :: Prob,Ran2,rat
!        LOGICAL :: WSign,KeepOrig
!        TYPE(HElement) :: rh
!
!!Calculate off-diagonal hamiltonian matrix element between determinants
!        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!        
!!rat is the probability of diffusing to nJ
!        rat=Tau*Lambda*abs(rh%v)/Prob
!
!        IF(rat.gt.1.D0) CALL Stop_All("AttemptDiffuse","*** Probability > 1 to diffuse.")
!
!        IF(rat.gt.Ran2(Seed)) THEN
!            IF(TExtraPartDiff) THEN
!!We want to perform the non-total number conserving diffusion matrix - anti-diffusion creates 2 new particles
!                IF(real(rh%v).gt.0.D0) THEN
!!Perform anti-diffusion - particle number will increase
!                    IF(WSign) THEN
!!We have a positive walker
!                        KeepOrig=.true.
!                        CreateAtI=1
!                        CreateAtJ=-1
!                    ELSE
!!We have a negative walker
!                        KeepOrig=.true.
!                        CreateAtI=-1
!                        CreateAtJ=1
!                    ENDIF
!                ELSE
!!Perform diffusion - particle number will remain constant
!                    IF(WSign) THEN
!                        KeepOrig=.false.    !Particle is annihilated by newly created opposing signed particle
!                        CreateAtI=0
!                        CreateAtJ=1
!                    ELSE
!                        KeepOrig=.false.
!                        CreateAtI=0
!                        CreateAtJ=-1
!                    ENDIF
!                ENDIF
!            ELSE
!!We are performing a number-conserving diffusion process - this will have a different diagonal birth/death unbiasing probability
!                IF(real(rh%v).gt.0.D0) THEN
!!Perform anti-diffusion, but conserve total particle number
!                    IF(WSign) THEN
!                        KeepOrig=.false.
!                        CreateAtI=0
!                        CreateAtJ=-1
!                    ELSE
!                        KeepOrig=.false.
!                        CreateAtI=0
!                        CreateAtJ=1
!                    ENDIF
!                ELSE
!!Perform diffusion - this should conserve particle number, and be the same as the other diffusion process
!                    IF(WSign) THEN
!                        KeepOrig=.false.
!                        CreateAtI=0
!                        CreateAtJ=1
!                    ELSE
!                        KeepOrig=.false.
!                        CreateAtI=0
!                        CreateAtJ=-1
!                    ENDIF
!                ENDIF
!            ENDIF
!        ELSE
!!No diffusion will occur...
!            KeepOrig=.true.        !We don't want to copy it accross, because it still can die - wait to see if it dies before copying accross
!            CreateAtI=0
!            CreateAtJ=0
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE AttemptDiffuse
!
!!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
!!Kik returns the a value which can be used to unbias the diffusion at the birth/death stage
!    INTEGER FUNCTION AttemptCreate(DetCurr,WSign,nJ,Prob,IC,Kik)
!        IMPLICIT NONE
!        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate
!        LOGICAL :: WSign
!        REAL*8 :: Prob,Ran2,rat,Kik
!        TYPE(HElement) :: rh
!
!!Calculate off diagonal hamiltonian matrix element between determinants
!        IF(TNoBirth) THEN
!            rh=HElement(0.D0)
!        ELSE
!            rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!        ENDIF
!
!!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
!        IF(TDiffuse) THEN
!            IF(TExtraPartDiff) THEN
!                Kik=Tau*rh%v/Prob
!            ELSE
!                Kik=Tau*abs(rh%v)/Prob
!            ENDIF
!            rat=(1.D0-Lambda)*abs(Kik)
!        ELSE
!            rat=Tau*abs(rh%v)/Prob
!        ENDIF
!
!!        IF(rat.lt.0.D0) THEN
!!            CALL Stop_All("AttemptCreate","*** Probability < 0 to create child.")
!!        ENDIF
!
!!If probability is > 1, then we can just create multiple children at the chosen determinant
!        ExtraCreate=INT(rat)
!        rat=rat-REAL(ExtraCreate)
!
!
!!Stochastically choose whether to create or not according to Ran2
!        IF(rat.gt.Ran2(Seed)) THEN
!!Child is created - what sign is it?
!            IF(WSign) THEN
!!Parent particle is positive
!                IF(real(rh%v).gt.0.D0) THEN
!                    AttemptCreate=-1     !-ve walker created
!                ELSE
!                    AttemptCreate=1      !+ve walker created
!                ENDIF
!
!            ELSE
!!Parent particle is negative
!                IF(real(rh%v).gt.0.D0) THEN
!                    AttemptCreate=1      !+ve walker created
!                ELSE
!                    AttemptCreate=-1     !-ve walker created
!                ENDIF
!            ENDIF
!
!        ELSE
!!No child particle created
!            AttemptCreate=0
!        ENDIF
!
!        IF(ExtraCreate.ne.0) THEN
!!Need to include the definitely create additional particles from a initial probability > 1
!
!            IF(AttemptCreate.lt.0) THEN
!!In this case particles are negative
!                AttemptCreate=AttemptCreate-ExtraCreate
!            ELSEIF(AttemptCreate.gt.0) THEN
!!Include extra positive particles
!                AttemptCreate=AttemptCreate+ExtraCreate
!            ELSEIF(AttemptCreate.eq.0) THEN
!!No particles were stochastically created, but some particles are still definatly created - we need to determinant their sign...
!                IF(WSign) THEN
!                    IF(real(rh%v).gt.0.D0) THEN
!                        AttemptCreate=-1*ExtraCreate    !Additional particles are negative
!                    ELSE
!                        AttemptCreate=ExtraCreate       !Additional particles are positive
!                    ENDIF
!                ELSE
!                    IF(real(rh%v).gt.0.D0) THEN
!                        AttemptCreate=ExtraCreate
!                    ELSE
!                        AttemptCreate=-1*ExtraCreate
!                    ENDIF
!                ENDIF
!            ENDIF
!        ENDIF
!
!        IF(Tau.lt.0.D0) THEN
!!If tau is negative, we are going back in time, and so will actually create antiparticles - flip sign again...
!            AttemptCreate=-AttemptCreate
!        ENDIF
!
!        IF(TDetPops) THEN
!!Here, we want to record the details of every spawning connection
!
!            StoreNumTo=0
!            StoreNumFrom=0
!            do i=1,SizeOfSpace
!                IF((DetLT(nJ,NMRKS(:,i),NEl).eq.0).and.(StoreNumTo.eq.0)) THEN
!!Found position of determinant
!                    StoreNumTo=i
!                ENDIF
!                IF((DetLT(DetCurr,NMRKS(:,i),NEl).eq.0).and.(StoreNumFrom.eq.0)) THEN
!                    StoreNumFrom=i
!                ENDIF
!                IF((StoreNumTo.ne.0).and.(StoreNumFrom.ne.0)) EXIT
!            enddo
!
!            IF(AttemptCreate.lt.0) THEN
!!If creating a negative particle, reduce the connection
!                TransMat(StoreNumFrom,StoreNumTo)=TransMat(StoreNumFrom,StoreNumTo)-(AttemptCreate*0.00001)
!                TransMat(StoreNumTo,StoreNumFrom)=TransMat(StoreNumTo,StoreNumFrom)-(AttemptCreate*0.00001)
!            ELSE
!!If creating a positive particle, increase connection
!                TransMat(StoreNumFrom,StoreNumTo)=TransMat(StoreNumFrom,StoreNumTo)+(AttemptCreate*0.00001)
!                TransMat(StoreNumTo,StoreNumFrom)=TransMat(StoreNumTo,StoreNumFrom)+(AttemptCreate*0.00001)
!            ENDIF
!
!        ENDIF
!
!        RETURN
!
!    END FUNCTION AttemptCreate
!
!!This function tells us whether we should kill the particle at determinant DetCurr
!!If also diffusing, then we need to know the probability with which we have spawned. This will reduce the death probability.
!!The function allows multiple births(if +ve shift) or deaths from the same particle.
!!The returned number is the number of deaths if positive, and the number of births if negative.
!    INTEGER FUNCTION AttemptDie(DetCurr,Kik,nExcitMemLen,nExcit,nStore)
!        IMPLICIT NONE
!        INTEGER :: DetCurr(NEl),DetLT,iKill,nExcitMemLen,nStore(6),nExcit(nExcitMemLen)
!        INTEGER :: nJ(NEl),IC,iMaxExcit,ierr,nExcitTag=0
!        CHARACTER(len=*) , PARAMETER :: this_routine='AttemptDie'
!        TYPE(HElement) :: rh,rhij
!        REAL*8 :: Ran2,rat,Kik
!
!!Test if determinant is FDet - in a strongly single-configuration problem, this will save time
!        IF(DetLT(DetCurr,FDet,NEl).eq.0) THEN
!            rh=0.D0
!        ELSE
!!Calculate the diagonal hamiltonian matrix element for the determinant
!            rh=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!!Subtract from the diagonal the value of the lowest hamiltonian matrix element
!            rh=rh-Hii
!        ENDIF
!
!        IF(TDiffuse) THEN
!!If also diffusing, then the probability of dying must be modified, since the diagonal elements have been altered
!            IF(TExtraPartDiff) THEN
!
!                IF(TFullUnbias) THEN
!
!
!!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
!                    
!                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
!                    
!                    Kik=0.D0    !If we are fully unbiasing, then the unbiasing factor is reset and recalculated from all excitations of DetCurr
!                    do while(.true.)
!                        CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
!                        IF(nJ(1).eq.0) EXIT
!                        rhij=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!                        Kik=Kik+(rhij%v)
!                    enddo
!
!                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
!
!                    Kik=Kik*Tau    !Unbias with the tau*sum of connected elements
!                ENDIF
!
!                rat=(Tau*((rh%v)-DiagSft))+(Lambda*Kik)     !This is now the probability with the correct unbiasing
!
!            ELSE
!                IF(TFullUnbias) THEN
!
!                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
!                    
!                    Kik=0.D0    !If we are fully unbiasing, then the unbiasing factor is reset and recalculated from all excitations of DetCurr
!                    do while(.true.)
!                        CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
!!                        CALL WRITEDET(6,nJ,NEl,.true.)
!                        IF(nJ(1).eq.0) EXIT
!                        rhij=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!                        Kik=Kik+abs(rhij%v)
!                    enddo
!
!                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
!
!                    Kik=Kik*Tau    !Unbias with the tau*sum of connected elements
!                ENDIF
!
!                rat=(Tau*((rh%v)-DiagSft))-(Lambda*Kik)     !This is now the probability with the correct unbiasing
!
!            ENDIF
!
!        ELSE
!!Subtract the current value of the shift and multiply by tau
!            rat=Tau*((rh%v)-DiagSft)
!        ENDIF
!
!!        IF(rat.gt.1.D0) THEN
!!If probs of dying is greater than one, reduce tau
!!            CALL Stop_All("AttemptDie","*** Death probability > 1. *** Tau too large")
!!        ENDIF
!
!        iKill=INT(rat)
!        rat=rat-REAL(iKill)
!
!!Stochastically choose whether to die or not
!        IF(abs(rat).gt.Ran2(Seed)) THEN
!            IF(rat.ge.0.D0) THEN
!!Depends whether we are trying to kill or give birth to particles.
!                iKill=iKill+1
!            ELSE
!                iKill=iKill-1
!            ENDIF
!        ENDIF
!
!        AttemptDie=iKill
!!        IF(AttemptDie.le.-1) WRITE(6,*) Iter,AttemptDie,rat
!
!        RETURN
!
!    END FUNCTION AttemptDie
!
!    SUBROUTINE MCDiffusion()
!        IMPLICIT NONE
!        INTEGER :: ierr,nExcitMemLen,TotExcits
!        CHARACTER(len=*) , PARAMETER :: this_routine='MCDiffusion'
!        INTEGER :: nJ(NEl),nK(NEl),ICJ,ICK,iCountJ,iCountK,i,j,iGetExcitLevel
!        REAL*8 :: ProbJ,ProbK,Ran2,rat,NewHii,SumDeathProb,ExpectedWalkers
!        TYPE(HElement) :: Hij,Hik,tempHii
!        TYPE(HElement) , ALLOCATABLE :: Hi0Array(:)
!        INTEGER , ALLOCATABLE :: ICWalk(:)
!        REAL*8 , ALLOCATABLE :: HiiArray(:)
!
!        CALL InitFCIMCCalc()    !Initialise walkers to be all positive and on the HF
!
!        WRITE(6,*) "Performing MC Diffusion with ",InitWalkers," walkers."
!        WRITE(6,*) " Iteration   Shift    ExpectedGrowthRate   ProjectionE"
!        WRITE(15,*) "# Iteration   Shift    ExpectedGrowthRate   ProjectionE"
!        
!        WRITE(6,"(I9,3G16.7)") 0,DiagSft,0.D0,0.D0
!        WRITE(15,"(I9,3G16.7)") 0,DiagSft,0.D0,0.D0
!
!!Need to store excitation generators for each of the particles
!        ALLOCATE(ExcitGens(InitWalkers),stat=ierr)     !Array to hold excitation generators for each walker
!        ALLOCATE(ICWalk(InitWalkers),stat=ierr)     !Array to hold excitation level for each walker
!        ALLOCATE(Hi0Array(InitWalkers),stat=ierr)   !Array to hold connection back to HF for each walker
!        ALLOCATE(HiiArray(InitWalkers),stat=ierr)   !Array to hold diagonal hamiltonian element for each walker
!        IF(ierr.ne.0) CALL Stop_All("MCDiffusion","Error in allocation")
!
!        CALL SetupExitgen(FDet,ExcitGens(1),nExcitMemLen,TotExcits)
!        ICWalk(1)=0
!        Hi0Array(1)=Hii
!        HiiArray(1)=real(Hii%v,r2)
!
!        do i=2,InitWalkers
!!Copy the excitation generator for FDet to all the other initial walkers
!            ALLOCATE(ExcitGens(i)%ExcitData(nExcitMemLen),stat=ierr)
!            IF(ierr.ne.0) CALL Stop_All("MCDiffusion","Error in allocation")
!            do j=1,nExcitMemLen
!                ExcitGens(i)%ExcitData(j)=ExcitGens(1)%ExcitData(j)
!            enddo
!            ICWalk(i)=0
!            Hi0Array(i)=Hii
!            HiiArray(i)=real(Hii%v,r2)
!        enddo
!
!        SumENum=real(Hii%v,r2)*InitWalkers
!        SumNoatHF=InitWalkers
!
!        SumDeathProb=0.D0
!
!        do Iter=1,NMCyc
!!Perform MC cycles
!            do j=1,InitWalkers
!!Cycle over all walkers
!
!!First we have to see if we're going to perform a self-hop, or allow an attempt at a diffusive move.
!!                rat=Tau/(HiiArray(j)-(real(Hii%v,r2))-DiagSft)      !This is the probability of self-hopping, rather than attempting a diffusive move
!!                rat=Tau*((HiiArray(j)/(real(Hii%v,r2)))-DiagSft)      !This is the probability of self-hopping, rather than attempting a diffusive move
!                rat=0.8*(1.D0-Tau*((HiiArray(j)-(real(Hii%v,r2)))-DiagSft))      !This is the probability of self-hopping, rather than attempting a diffusive move
!                                                        !Higher excitations have smaller prob, so resampled fewer times and lower excitations sampled for longer.
!
!                IF((rat.lt.0.D0).or.(rat.gt.1.D0)) CALL Stop_All("DiffusionMC","Incorrect self-hop probability")
!                IF(rat.gt.Ran2(Seed)) THEN
!!We want to self-hop. Resum in energy, but to not allow an attempted move away from excit. We also update the effect on the shift later.
!                        
!                        SumENum=SumENum+(real(Hi0Array(j)%v,r2))
!                        IF(ICWalk(j).eq.0) SumNoatHF=SumNoatHF+1
!                        NewHii=HiiArray(j)
!
!                ELSE
!!We do not have to self-hop. Attempt to diffuse away from determinant.
!                
!                    CALL GenRandSymExcitIt3(CurrentDets(:,j),ExcitGens(j)%ExcitData,nJ,Seed,ICJ,0,ProbJ,iCountJ)  !First random excitation is to attempt to move to
!                    CALL GenRandSymExcitIt3(CurrentDets(:,j),ExcitGens(j)%ExcitData,nK,Seed,ICK,0,ProbK,iCountK)  !Second is to unbias the diffusion in the birth/death prob
!                    Hij=GetHElement2(CurrentDets(:,j),nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,ICJ,ECore)
!                    Hik=GetHElement2(CurrentDets(:,j),nK,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,ICK,ECore)
!
!!Attempt diffusion away to nJ
!                    rat=Tau*abs(real(Hij%v,r2))/ProbJ
!                
!                    IF(rat.gt.1.D0) CALL Stop_All("AttemptDiffuse","*** Probability > 1 to diffuse.")
!
!                    IF(rat.gt.Ran2(Seed)) THEN
!!Diffusion successful - need to update all the information
!                   
!                        CurrentDets(:,j)=nJ(:)
!
!                        IF(real(Hij%v,r2).gt.0.D0) THEN
!!This is the anti-diffusion
!                            IF(CurrentSign(j)) THEN
!!Walker is positive
!                                CurrentSign(j)=.false.  !Want negative walker
!                            ELSE
!                                CurrentSign(j)=.true.   !Positive walker
!                            ENDIF
!                        ELSE
!!This is diffusion - no need to flip sign, since we want to keep the same signs.
!                        ENDIF
!
!!Sum in new energy
!                        ICWalk(j)=iGetExcitLevel(FDet,nJ,NEl)
!                        IF(ICWalk(j).eq.2) THEN
!                            Hi0Array(j)=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,2,ECore)
!                            SumENum=SumENum+(real(Hi0Array(j)%v,r2))
!                        ELSEIF(ICWalk(j).eq.0) THEN
!                            Hi0Array(j)=Hii
!                            SumENum=SumENum+(real(Hii%v,r2))
!                            SumNoatHF=SumNoatHF+1
!                        ELSE
!                            Hi0Array(j)=HElement(0.D0)
!                        ENDIF
!
!!Create and save new excitation generator
!                        CALL SetupExitgen(nJ,ExcitGens(j),nExcitMemLen,TotExcits)
!
!                        tempHii=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
!                        NewHii=real(tempHii%v,r2)  !This is the new diagonal matrix element
!                                                                                                       
!                    ELSE
!!Attempted diffusion away is not successful - still need to update the energy
!                    
!                        SumENum=SumENum+(real(Hi0Array(j)%v,r2))
!                        IF(ICWalk(j).eq.0) SumNoatHF=SumNoatHF+1
!                        NewHii=HiiArray(j)
!
!                    ENDIF
!
!                ENDIF
!
!                rat=Tau*(((HiiArray(j)-(real(Hii%v,r2)))-DiagSft)-(abs(real(Hik%v,r2))/ProbK))      !This is the prob of death, adjusted to unbias for the diffusion
!
!                HiiArray(j)=NewHii      !Now we can update the HiiArray to take into account if the walker has moved to a new position
!
!                SumDeathProb=SumDeathProb+rat   !Sum the death probabilities
!
!            enddo   !Finsh cycling over all walkers
!            
!            ProjectionE=(SumENum/(SumNoatHF+0.D0))-REAL(Hii%v,r2)
!
!            IF(mod(Iter,StepsSft).eq.0) THEN
!        
!                SumDeathProb=SumDeathProb/(InitWalkers*StepsSft+0.D0)       !This now indicates the average death probability per walker per cycle
!!                WRITE(6,*) SumDeathProb
!
!!                ExpectedWalkers=InitWalkers+0.D0    !This is the number of walkers we start off with at the beginning of the update shift cycle.
!!                do i=1,StepsSft
!!                    ExpectedWalkers=ExpectedWalkers-ExpectedWalkers*SumDeathProb    
!!                enddo                                                               
!
!                IF(SumDeathProb.gt.1.D0) CALL Stop_All("MCDIFFUSION","Average Death prob. > 1")
!
!                ExpectedWalkers=InitWalkers*((1.D0-SumDeathProb)**StepsSft)
!
!                GrowRate=(ExpectedWalkers)/(InitWalkers+0.D0)  !This is the expected growth rate over the previous StepsSft
!
!                DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
!
!                SumDeathProb=0.D0
!                            
!                WRITE(6,"(I9,3G16.7)") Iter,DiagSft,GrowRate,ProjectionE
!                WRITE(15,"(I9,3G16.7)") Iter,DiagSft,GrowRate,ProjectionE
!
!                CALL FLUSH(6)
!                CALL FLUSH(15)
!
!            ENDIF
!
!        enddo
!
!        RETURN
!
!    END SUBROUTINE MCDiffusion
!   
!    SUBROUTINE SetupExitgen(nI,ExcitGen,nExcitMemLen,iMaxExcit)
!        IMPLICIT NONE
!        TYPE(ExcitGenerator) :: ExcitGen
!        INTEGER :: ierr,iMaxExcit,nExcitMemLen,nJ(NEl)
!        INTEGER :: nI(NEl)
!
!        IF(Allocated(ExcitGen%ExcitData)) THEN
!            DEALLOCATE(ExcitGen%ExcitData)
!        ENDIF
!
!!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
!        iMaxExcit=0
!        CALL IAZZERO(ExcitGen%nStore,6)
!        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,ExcitGen%nStore,3)
!        ALLOCATE(ExcitGen%ExcitData(nExcitMemLen),stat=ierr)
!        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
!        ExcitGen%ExcitData(1)=0
!        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGen%ExcitData,nJ,iMaxExcit,0,ExcitGen%nStore,3)
!
!    END SUBROUTINE SetupExitgen
!
!END MODULE FciMCMod 
!
