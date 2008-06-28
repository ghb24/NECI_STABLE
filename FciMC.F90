!!!!!!!!!
! TO DO !
!!!!!!!!!
!Test Dynamic Tau Changing
!Initial Graph-Morph Graph
!Test of LD with MP1
!Test POPSFILE read/write
!Setup Init Star

MODULE FciMCMod
    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr
    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    USE Calc , only : TReadPops,ScaleWalkers,TMCExcitSpace,NoMCExcits,TStartMP1
    USE Calc , only : GrowMaxFactor,CullFactor,TMCDets,TNoBirth,Lambda,TDiffuse,FlipTauCyc,TFlipTau
    USE Calc , only : TExtraPartDiff,TFullUnbias,TNodalCutoff,NodalCutoff,TNoAnnihil
    USE Determinants , only : FDet,GetHElement2
    USE DetCalc , only : NMRKS
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE Logging , only : TPopsFile,TCalcWavevector,WavevectorPrint,TDetPops
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    IMPLICIT NONE
    SAVE

    INTEGER, PARAMETER :: r2=kind(0.d0)

    INTEGER , ALLOCATABLE , TARGET :: WalkVecDets(:,:),WalkVec2Dets(:,:)
    LOGICAL , ALLOCATABLE , TARGET :: WalkVecSign(:),WalkVec2Sign(:)
    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0

!Pointers to point at the correct arrays for use
    INTEGER , POINTER :: CurrentDets(:,:)
    LOGICAL , POINTER :: CurrentSign(:)
    INTEGER , POINTER :: NewDets(:,:)
    LOGICAL , POINTER :: NewSign(:)

!InitPopsVector is used to store the initial populations of the determinants. This can be used for comparison when TNoBirth is true.
    INTEGER , ALLOCATABLE :: InitPops(:)
    INTEGER :: InitPopsTag=0

    INTEGER , ALLOCATABLE :: ExcitStore(:,:)
    INTEGER :: ExcitStoreTag=0

    REAL*8 , ALLOCATABLE :: TransMat(:,:),PopsVec(:)
    INTEGER :: TransMatTag=0,PopsVecTag=0,SizeOfSpace

!MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
    INTEGER :: MemoryFac=10000

    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps
    INTEGER :: exFlag=3

!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle
    INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
    INTEGER :: CullInfo(10,3)

    REAL*8 :: GrowRate,DieRat,MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector
    REAL*8 :: ProjectionE,SumENum,SumE,ProjectionEInst

    INTEGER :: CycwNoHF     !Count the number of iterations which don't have any walkers on HF - in this case, ignore running average of energy, but count it so it is not included in demonimator

    INTEGER*8 :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant

    INTEGER :: NetPositive

    TYPE(HElement) :: Hii,rhii,FZero

    contains

    SUBROUTINE FciMC(Weight,Energyxw)
        Use MCDets, only : MCDetsCalc
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j,iSub,WalkOnDet,DetLT,DetCurr(NEl),ExpectedDets
        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'
        TYPE(HElement) :: Hamii

        if(NMCyc.lt.0.or.TMCDets) then

           CALL MCDetsCalc(FDet,G_VMC_Seed,-NMCyc,Tau,SftDamp,10000*InitWalkers,InitWalkers,StepsSft,DiagSft,GrowMaxFactor,CullFactor)
            Weight=1.d0
            Energyxw=1.d0
            return
        endif
        CALL TISET('FCIMC',iSub)

        IF(TDiffuse) THEN
            IF((.NOT.TMCExcitSpace).or.(NoMCExcits.ne.1)) THEN
                CALL STOPGM("FCIMC","Diffusion can only work with one MCExcitSpace")
            ENDIF
            IF((Lambda.gt.1.D0).or.(Lambda.lt.0.D0)) THEN
                CALL STOPGM("FCIMC","Diffusion coefficient must be between zero and one")
            ENDIF
        ENDIF

        IF(HElementSize.gt.1) THEN
            CALL STOPGM("FCIMC","StarDiagMC cannot function with complex orbitals.")
        ENDIF

        WRITE(6,*) ""
        WRITE(6,*) "Performing FCIMC...."

        IF(TCalcWaveVector) THEN
            WRITE(6,*) "Wavevector calculation is only available in star MC..."
        ENDIF

        OPEN(15,file='FCIMCStats',status='unknown')

        IF(TReadPops) THEN
            
            IF(TNoBirth) THEN
                WRITE(6,*) "Cannot use NOBIRTHS with READPOPS"
                STOP "Cannot use NOBIRTHS with READPOPS"
            ENDIF

            WRITE(6,*) "Reading in POPSFILE to restart calculation..."
            OPEN(17,FILE='POPSFILE',Status='old')
            READ(17,*) InitWalkers
            READ(17,*) DiagSft
            READ(17,*) PreviousNMCyc
            WRITE(6,*) "Initial number of walkers read to be: ", InitWalkers
        ELSE
            WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
        ENDIF

        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
        IF(TFlipTau) THEN
            WRITE(6,*) "Flipping the sign of Tau once every :", FlipTauCyc
        ENDIF

!Initialise random number seed
        Seed=G_VMC_Seed
!Initialise variables for calculation of the running average
        ProjectionE=0.D0
        ProjectionEInst=0.D0
        SumE=0.D0
        SumENum=0.D0
        SumNoatHF=0
        CycwNoHF=0

!Calculate Hii
        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)

!        IF(DiagSft.gt.0.D0) THEN
!            CALL StopGM("StarDiagMC","Intial value of DiagSft should be negative.")
!        ELSE
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
!        ENDIF

        CALL InitFCIMCCalc()

        IF(.NOT.TNoBirth) THEN
!Print out initial starting configurations
            WRITE(6,*) ""
            WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers        Proj.E      Net+veWalk     Proj.E-Inst"
            WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers         Proj.E      Net+veWalk     Proj.E-Inst"
!TotWalkersOld is the number of walkers last time the shift was changed
            IF(TReadPops) THEN
                WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
                WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
            ELSE
                WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
                WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
            ENDIF
        ENDIF
        
!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers

!Start MC simulation...
        do Iter=1,NMCyc
            
            CALL PerformFCIMCyc()

!            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
!            CALL FLUSH(6)

            IF(TFlipTau) THEN
!If we are flipping the sign of tau every FlipTauCyc cycles
                IF(mod(Iter,FlipTauCyc).eq.0) THEN
                    Tau=Tau*(-1.D0)
                    DiagSft=DiagSft*(-1.D0)
                ENDIF
            ENDIF

            IF(mod(Iter,StepsSft).eq.0) THEN

                IF(TNoBirth) THEN
!If we have TNoBirth flag, the populations on each determinant simply exponentially decay, and so at each StepsSft, we want to know the population on each determinant
!When creating wavevector, store the initial number of particles on each determinant (order them).
                    do i=1,NoComps
                        do j=1,NEl
                            DetCurr(j)=ExcitStore(j,i)
                        enddo
                        WalkOnDet=0
                        do j=1,TotWalkers
                            IF(DetLT(DetCurr,CurrentDets(:,j),NEl).eq.0) THEN
!Walker is on chosen determinant
                                WalkOnDet=WalkOnDet+1
                            ENDIF
                        enddo

                        Hamii=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                        Hamii=Hamii-Hii
!The expected number of determinants now is the original number on each determinant * exp(-tau*Iter*(Hii-DiagSft)
                        ExpectedDets=exp(-Tau*Iter*((Hamii%v)-DiagSft))*InitPops(i)

                        WRITE(6,*) Iter,i,WalkOnDet,ExpectedDets
                        WRITE(15,*) Iter,i,WalkOnDet,ExpectedDets
                        CALL FLUSH(15)
                        CALL FLUSH(6)

                    enddo
                    WRITE(6,*) ""
                    WRITE(15,*) ""

!This can be compared to the actual number of particles on each of the determinants.
!This will only work when we start with the MP1 wavefunction.

                ELSE
!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
                    CALL UpdateDiagSft()

!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
                    IF(TReadPops) THEN
                        WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
                        WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
                    ELSE
                        IF(Tau.gt.0.D0) THEN
                            WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
                            WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
                        ELSE
                            WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
                            WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,NetPositive,ProjectionEInst
                        ENDIF
                    ENDIF

                    IF(TDetPops) THEN

                        do i=1,SizeOfSpace
                            WRITE(67,"(F10.5,$)") TransMat(1,i)
                        enddo
                        WRITE(67,*) ""

                    ENDIF


                    CALL FLUSH(15)
                    CALL FLUSH(6)

!Reset TotWalkersOld so that it is the number of walkers now
                    TotWalkersOld=TotWalkers
                ENDIF

!                WRITE(6,*) DieRat

            ENDIF

!End of MC cycle
        enddo

        IF(TPopsFile) THEN
!Print out current state of simulation, so it can be restarted if desired...
            OPEN(16,file='POPSFILE',status='unknown')
            WRITE(16,*) TotWalkers, "   TOTWALKERS"
            WRITE(16,*) DiagSft, "   DIAG SHIFT"
            IF(TReadPops) THEN
                WRITE(16,*) NMCyc+PreviousNMCyc, "   NO. CYCLES"
            ELSE
                WRITE(16,*) NMCyc, "   MC CYCLES"
            ENDIF
            do i=1,TotWalkers
                WRITE(16,*) CurrentDets(:,i),CurrentSign(i)
            enddo
            CLOSE(16)
        ENDIF

        Weight=HDElement(1.D0)
        Energyxw=HDElement(2.D0*DiagSft)

!Deallocate memory
        IF(TNoBirth) THEN
            DEALLOCATE(ExcitStore)
            CALL LogMemDealloc(this_routine,ExcitStoreTag)
            DEALLOCATE(InitPops)
            CALL LogMemDealloc(this_routine,InitPopsTag)
        ENDIF
        IF(TDetPops) THEN
            DEALLOCATE(ExcitStore)
            CALL LogMemDealloc(this_routine,ExcitStoreTag)
            DEALLOCATE(PopsVec)
            CALL LogMemDealloc(this_routine,PopsVecTag)
            DEALLOCATE(TransMat)
            CALL LogMemDealloc(this_routine,TransMatTag)
        ENDIF
        DEALLOCATE(WalkVecDets)
        CALL LogMemDealloc(this_routine,WalkVecDetsTag)
        DEALLOCATE(WalkVec2Dets)
        CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
        DEALLOCATE(WalkVecSign)
        CALL LogMemDealloc(this_routine,WalkVecSignTag)
        DEALLOCATE(WalkVec2Sign)
        CALL LogMemDealloc(this_routine,WalkVec2SignTag)

        CLOSE(15)

        CALL TIHALT('FCIMC',iSub)

        RETURN

    END SUBROUTINE FciMC

!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCyc()
        IMPLICIT NONE
        INTEGER :: VecSlot,i,j,k,l,DetCurr(NEl),iMaxExcit,nExcitMemLen,nStore(6)
        INTEGER :: nJ(NEl),ierr,nExcitTag=0,IC,Child,iSubCyc,TotWalkersNew,iCount
        REAL*8 :: Prob,rat,Kik
        INTEGER , ALLOCATABLE :: nExcit(:)
        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
        LOGICAL :: WSign
        LOGICAL :: KeepOrig
        INTEGER :: CreateAtI,CreateAtJ,tocopy
        CHARACTER(len=*), PARAMETER :: this_routine='PerformFCIMCyc'
        
        CALL TISET('MCyc',iSubCyc)
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1

        do j=1,TotWalkers
!j runs through all current walkers
            do k=1,NEl
                DetCurr(k)=CurrentDets(k,j)
            enddo

!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
            iMaxExcit=0
            CALL IAZZERO(nStore,6)
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
            nExcit(1)=0
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

            IF(TMCExcitSpace) THEN
!Excitations are picked stochastically

                do i=1,NoMCExcits

                    CALL GenRandSymExcitIt3(DetCurr,nExcit,nJ,Seed,IC,0,Prob,iCount)

                    Child=AttemptCreate(DetCurr,CurrentSign(j),nJ,Prob,IC,Kik)
!Kik is the off-diagonal hamiltonian matrix element for the walker. This is used for the augmentation of the death term if TDiffuse is on.
                    IF(Child.gt.0) THEN
!We have successfully created at least one positive child at nJ
                        WSign=.true.
                    ELSE
!We have successfully created at least one negative child at nJ
                        WSign=.false.
                    ENDIF
                    do l=1,abs(Child)
                        do k=1,NEl
                            NewDets(k,VecSlot)=nJ(k)
                        enddo
                        NewSign(VecSlot)=WSign
                        VecSlot=VecSlot+1
                    enddo

                enddo

            ELSE
!Run through all possible excitations of each walker
                
                CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)

                do while(.true.)
                    CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
                    IF(nJ(1).eq.0) EXIT

                    Child=AttemptCreate(DetCurr,CurrentSign(j),nJ,1.D0,IC,Kik)
!Kik is the off-diagonal hamiltonian matrix element for the walker. This is used for the augmentation of the death term if TDiffuse is on.
                    IF(Child.gt.0) THEN
!We have successfully created at least one positive child at nJ
                        WSign=.true.
                    ELSE
!We have successfully created at least one negative child at nJ
                        WSign=.false.
                    ENDIF
                    do l=1,abs(Child)
                        do k=1,NEl
                            NewDets(k,VecSlot)=nJ(k)
                        enddo
                        NewSign(VecSlot)=WSign
                        VecSlot=VecSlot+1
                    enddo

                enddo

            ENDIF

            KeepOrig=.true.
            IF(TDiffuse) THEN
!Next look at possibility of diffusion to another determinant
                CALL GenRandSymExcitIt3(DetCurr,nExcit,nJ,Seed,IC,0,Prob,iCount)
                CALL AttemptDiffuse(DetCurr,nJ,Prob,IC,CurrentSign(j),KeepOrig,CreateAtI,CreateAtJ)
                !If we want to keep the original walker, then KeepOrig is true, However, we do not want to copy it accross yet, because we want to see if it is killed first in the birth/death process
!                IF(KeepOrig) THEN
!                    do k=1,NEl
!                        NewDets(k,VecSlot)=DetCurr(k)
!                    enddo
!                    NewSign(VecSlot)=CurrentSign(j)
!                    VecSlot=VecSlot+1
!                ENDIF
                do l=1,abs(CreateAtI)       !Sum in the number of walkers to create at the original determinant
                    do k=1,NEl
                        NewDets(k,VecSlot)=DetCurr(k)
                    enddo
                    IF(CreateAtI.gt.0) THEN
                        NewSign(VecSlot)=.true.
                    ELSE
                        NewSign(VecSlot)=.false.
                    ENDIF
                    VecSlot=VecSlot+1
                enddo
                do l=1,abs(CreateAtJ)       !Add the number of walkers to create at nJ
                    do k=1,NEl
                        NewDets(k,VecSlot)=nJ(k)
                    enddo
                    IF(CreateAtJ.gt.0) THEN
                        NewSign(VecSlot)=.true.
                    ELSE
                        NewSign(VecSlot)=.false.
                    ENDIF
                    VecSlot=VecSlot+1
                enddo
            ENDIF

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
            iDie=AttemptDie(DetCurr,Kik,nExcitMemLen,nExcit,nStore)
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births

            IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle
        
                IF(KeepOrig) THEN
                    ToCopy=abs(iDie)+1  !This is because we need to copy accross the original particle too
                ELSE
                    IF(iDie.eq.0) THEN
                        ToCopy=0        !Indicates that want to spare original particle, which has already been copied accross previously, or previously been annihilated
                    ELSE
                        ToCopy=abs(iDie)
                    ENDIF
                ENDIF

                do l=1,ToCopy    !We need to copy accross one more, since we need to include the original spared particle
                    do k=1,NEl
                        NewDets(k,VecSlot)=DetCurr(k)
                    enddo
                    NewSign(VecSlot)=CurrentSign(j)
                    VecSlot=VecSlot+1
                enddo

            ELSEIF(iDie.gt.0) THEN
!This indicates that particles on DetCurr want to be killed. The first kill will simply be performed by not copying accross the original particle.
!Therefore, if iDie = 1, then we can simply ignore it.
!However, after that anti-particles will need to be created on the same determinant.

                IF(KeepOrig) iDie=iDie-1    
!This is because we can already kill one particle by not copying accross the particle which was originally there (and still is even after diffusion). 
!If KeepOrig is false, then the particle has already diffused somewhere else, and so antiparticles need to be created in its place.

                do l=1,iDie
                    do k=1,NEl
                        NewDets(k,VecSlot)=DetCurr(k)
                    enddo
                    IF(CurrentSign(j)) THEN
!Copy accross new anti-particles
                        NewSign(VecSlot)=.FALSE.
                    ELSE
                        NewSign(VecSlot)=.TRUE.
                    ENDIF
                    VecSlot=VecSlot+1
                enddo
            
            ENDIF

!Destroy excitation generators for current walker
            DEALLOCATE(nExcit)
            CALL LogMemDealloc(this_routine,nExcitTag)
        
!            rat=(VecSlot+0.D0)/(MaxWalkers+0.D0)
!            IF(rat.gt.0.9) THEN
!                WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
!            ENDIF

!Finish cycling over walkers
        enddo

!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
        TotWalkersNew=VecSlot-1
        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
        ENDIF
        
        IF(TNodalCutoff.and.(NodalCutoff.lt.0.D0)) THEN
!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff
!then the net number of walkers must be the same sign, or it is set to zero walkers, and they are killed.
            CALL TestWavevectorNodes(TotWalkersNew,2)
        ENDIF

        IF(TNoAnnihil) THEN
!We are not annihilating particles - this will make things much quicker.

!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
            CurrentDets=>WalkVec2Dets
            CurrentSign=>WalkVec2Sign
            NewDets=>WalkVecDets
            NewSign=>WalkVecSign

            TotWalkers=TotWalkersNew

        ELSE
!This routine now cancels down the particles with opposing sign on each determinant
!This routine does not necessarily need to be called every Iter, but it does at the moment, since it is the only way to 
!transfer the residual particles back onto CurrentDets
            CALL AnnihilatePairs(TotWalkersNew)
!            WRITE(6,*) "Number of annihilated particles= ",TotWalkersNew-TotWalkers
        ENDIF


        IF(TNodalCutoff.and.(NodalCutoff.ge.0.D0)) THEN
!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff
!then the net number of walkers must be the same sign, or it is set to zero walkers, and they are killed.
            CALL TestWavevectorNodes(TotWalkers,1)
        ENDIF

        
        IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!Particle number is too large - kill them randomly

!Log the fact that we have made a cull
            NoCulls=NoCulls+1
            IF(NoCulls.gt.10) CALL STOPGM("PerformFCIMCyc","Too Many Culls")
!CullInfo(:,1) is walkers before cull
            CullInfo(NoCulls,1)=TotWalkers
!CullInfo(:,3) is MC Steps into shift cycle before cull
            CullInfo(NoCulls,3)=mod(Iter,StepsSft)

            WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number..."
            WRITE(6,*) "Killing randomly selected particles in order to reduce total number..."
            WRITE(6,"(A,F8.2)") "Population will reduce by a factor of ",CullFactor
            CALL ThermostatParticles(.true.)

!Need to reduce totwalkersOld, so that the change in shift is also reflected by this
!The Shift is no longer calculated like this...
!            TotWalkersOld=nint((TotWalkersOld+0.D0)/CullFactor)

        ELSEIF((TotWalkers.lt.(InitWalkers/2)).and.(.NOT.TNoBirth)) THEN
!Particle number is too small - double every particle in its current position

!Log the fact that we have made a cull
            NoCulls=NoCulls+1
            IF(NoCulls.gt.10) CALL STOPGM("PerformFCIMCyc","Too Many Culls")
!CullInfo(:,1) is walkers before cull
            CullInfo(NoCulls,1)=TotWalkers
!CullInfo(:,3) is MC Steps into shift cycle before cull
            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
            
            WRITE(6,*) "Doubling particle population to increase total number..."
            CALL ThermostatParticles(.false.)

!Need to increase TotWalkersOld, so that the change in shift is also reflected by this
!The shift is no longer calculated like this...
!            TotWalkersOld=TotWalkersOld*2

        ENDIF

        CALL TIHALT('MCyc',iSubCyc)

        RETURN

    END SUBROUTINE PerformFCIMCyc

!If TNodalCutoff is set, then we are imposing a nodal boundary on the wavevector - if the MP1 wavefunction has a component which is larger than a NodalCutoff,
!then the walkers at that determinant must be of the same sign, or they are killed.
!Particles indicates the number of particles to look through, and iArray is 1 if the WalkVecDets array is active, and 2 if the WalkVec2Dets array is active.
    SUBROUTINE TestWavevectorNodes(Particles,iArray)
        USE Calc , only : i_P
        USE System , only : Beta
        USE Integrals , only : nTay
        IMPLICIT NONE
        INTEGER :: ParticlesOrig,VecSlot,IC,iGetExcitLevel,j,k,NoatHF,Particles,iArray,NoPositive,NoNegative
        INTEGER , POINTER :: ActiveVecDets(:,:)
        LOGICAL , POINTER :: ActiveVecSign(:)
        REAL*8 :: EigenComp,EnergyNum
        LOGICAL :: Component
        TYPE(HElement) :: Hamij,Fj!,rhjj,rhij

!We first need to point to the active array
        IF(iArray.eq.1) THEN
            ActiveVecDets => CurrentDets
            ActiveVecSign => CurrentSign
        ELSEIF(iArray.eq.2) THEN
            ActiveVecDets => NewDets
            ActiveVecSign => NewSign
        ELSE
            CALL STOPGM("TestWavevectorNodes","Error with iArray")
        ENDIF

        EnergyNum=0.D0  !EnergyNum indicates the sum over all the hamiltonian matrix elements between the double excitations and HF
        NoatHF=0
        NoPositive=0    !Total number of positive particles
        NoNegative=0    !Total number of negative particles

        ParticlesOrig=Particles
!VecSlot indicates the next free position in ActiveVec
        VecSlot=1

        do j=1,Particles
!j runs through all current walkers
!IC labels the excitation level away from the HF determinant
            IC=iGetExcitLevel(FDet,ActiveVecDets(:,j),NEl)

            IF((IC.gt.2).or.(IC.eq.1)) THEN
!More than double excitations/singles are not included in the MP1 wavefunction, therefore we can copy these straight across
                IF(VecSlot.lt.j) THEN
!Only copy accross to VecSlot if VecSlot<j, otherwise, if VecSlot=j, particle already in right place.
                    do k=1,NEl
                        ActiveVecDets(k,VecSlot)=ActiveVecDets(k,j)
                    enddo
                    ActiveVecSign(VecSlot)=ActiveVecSign(j)
                ENDIF
                VecSlot=VecSlot+1

                IF(ActiveVecSign(j)) THEN
                    NoPositive=NoPositive+1
                ELSE
                    NoNegative=NoNegative+1
                ENDIF

            ELSEIF((IC.eq.2).or.(IC.eq.0)) THEN
                IF(IC.eq.2) THEN
!We are at a double excitation - first find the desired sign of the determinant.
                    Hamij=GetHElement2(FDet,ActiveVecDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
                    IF((real(Hamij%v)).lt.0) THEN
!Negative Hamiltonian connection indicates that the component of the MP1 wavevector is positive
                        Component=.true.
                    ELSE
                        Component=.false.
                    ENDIF
                ELSE
!We are at the HF determinant
                    Hamij=Hii
                    Component=.true.    !HF component of MP1 wavevector always wants to be positive
                ENDIF

                IF((ActiveVecSign(j).and.Component).or.((.NOT.ActiveVecSign(j)).and.(.NOT.Component))) THEN
!The particle is of the same sign as the MP1 wavevector component, so whether or not the determinant is of fixed sign, the particle does not want to be destroyed
                    IF(VecSlot.lt.j) THEN
!Only copy accross to VecSlot if VecSlot<j, otherwise, if VecSlot=j, particle already in right place.
                        do k=1,NEl
                            ActiveVecDets(k,VecSlot)=ActiveVecDets(k,j)
                        enddo
                        ActiveVecSign(VecSlot)=ActiveVecSign(j)
                    ENDIF
                    VecSlot=VecSlot+1

!Add to the estimate for the energy if we want to keep the particle
                    IF(ActiveVecSign(j)) THEN
                        EnergyNum=EnergyNum+(REAL(Hamij%v,r2))
                        NoPositive=NoPositive+1
                    ELSE
                        EnergyNum=EnergyNum-(REAL(Hamij%v,r2))
                        NoNegative=NoNegative+1
                    ENDIF
                    IF(IC.eq.0) THEN
                        IF(ActiveVecSign(j)) THEN
                            NoatHF=NoatHF+1
                        ELSE
                            NoatHF=NoatHF-1
                        ENDIF
                    ENDIF

                ELSE
!Particle is of a different sign to the component of the MP1 wavefunction - if the determinant is of fixed sign, we don't want to copy it across. Find if fixed sign...
                    IF(IC.eq.2) THEN
                        CALL GetH0Element(ActiveVecDets(:,j),NEl,Arr,nBasis,ECore,Fj)
!EigenComp is now the component of the normalised MP1 wavefunction for the double excitation WalkVecDets(:,j)
                        EigenComp=abs(((Hamij%v)/(Fj%v-FZero%v))/MPNorm)
                    ELSE
!Here, EigenComp is the component of the normalised MP1 wavefunction at the HF determinant (should always be +ve)
                        EigenComp=1.D0/MPNorm
                    ENDIF

                    IF(EigenComp.lt.abs(NodalCutoff)) THEN
!Determinant does not have a fixed sign - keep the particle
                        IF(VecSlot.lt.j) THEN
!Only copy accross to VecSlot if VecSlot<j, otherwise, if VecSlot=j, particle already in right place.
                            do k=1,NEl
                                ActiveVecDets(k,VecSlot)=ActiveVecDets(k,j)
                            enddo
                            ActiveVecSign(VecSlot)=ActiveVecSign(j)
                        ENDIF
                        VecSlot=VecSlot+1

!Add to the estimate for the energy if we want to keep the particle
                        IF(ActiveVecSign(j)) THEN
                            EnergyNum=EnergyNum+(REAL(Hamij%v,r2))
                            NoPositive=NoPositive+1
                        ELSE
                            EnergyNum=EnergyNum-(REAL(Hamij%v,r2))
                            NoNegative=NoNegative+1
                        ENDIF
                        IF(IC.eq.0) THEN
                            IF(ActiveVecSign(j)) THEN
                                NoatHF=NoatHF+1
                            ELSE
                                NoatHF=NoatHF-1
                            ENDIF
                        ENDIF

                    ENDIF
                ENDIF
            ELSE
                CALL STOPGM("TestWavevectorNodes","Should not be here - wrong IC calculated")
            ENDIF
        enddo   !end loop over all walkers

!Total number of particles now modified by killing of wrong signed particles
        Particles=VecSlot-1
    
        IF((ParticlesOrig-Particles).ne.0) THEN
            WRITE(6,*) "Some particles killed by Nodal approximation: ", ParticlesOrig-Particles
        ENDIF

!Calculate the time average of the numerator and demonimator to calculate the running average of the energy - this will then not be affected by the case when there aren't any particles at the HF determinant
        SumNoatHF=SumNoatHF+NoatHF
        SumENum=SumENum+EnergyNum
        ProjectionE=(SumENum/(SumNoatHF+0.D0))-REAL(Hii%v,r2)

        IF(NoatHF.ne.0) THEN
!The energy cannot be calculated via the projection back onto the HF if there are no particles at HF
            SumE=SumE+((EnergyNum/(NoatHF+0.D0))-(REAL(Hii%v,r2)))
            ProjectionEInst=SumE/((Iter-CycwNoHF)+0.D0)
        ELSE
            CycwNoHF=CycwNoHF+1         !Record the fact that there are no particles at HF in this run, so we do not bias the average
!            WRITE(6,*) "No positive particles at reference determinant during iteration: ",Iter
        ENDIF

        NetPositive=NoPositive-NoNegative
!        WRITE(6,*) ParticlesOrig,Particles,NetPositive,NoPositive,NoNegative

        RETURN

    END SUBROUTINE TestWavevectorNodes

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

!Move the Walker at the end of the list to the position of the walker we have chosen to destroy
                do i=1,NEl
                    CurrentDets(i,Chosen)=CurrentDets(i,TotWalkers)
                enddo
                CurrentSign(Chosen)=CurrentSign(TotWalkers)

                TotWalkers=TotWalkers-1
                Culled=Culled+1

            enddo

            IF(TotWalkers.ne.(OrigWalkers-ToCull)) THEN
                WRITE(6,*) "Error in culling walkers..."
                STOP "Error in culling walkers..."
            ENDIF

!CullInfo(:,2) is the new number of total walkers
            CullInfo(NoCulls,2)=TotWalkers

        ELSE
!The population is too low - give it a boost by doubling every particle

            VecSlot=TotWalkers+1
            do i=1,TotWalkers

!Add clone of walker, at the same determinant, to the end of the list
                do j=1,NEl
                    CurrentDets(j,VecSlot)=CurrentDets(j,i)
                enddo
                CurrentSign(VecSlot)=CurrentSign(i)

                VecSlot=VecSlot+1

            enddo

            TotWalkers=TotWalkers*2

            IF((VecSlot-1).ne.TotWalkers) THEN
                WRITE(6,*) "Problem in doubling all particles..."
                STOP "Problem in doubling all particles..."
            ENDIF

!CullInfo(:,2) is the new number of total walkers
            CullInfo(NoCulls,2)=TotWalkers

        ENDIF

        RETURN

    END SUBROUTINE ThermostatParticles


!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!value of the diagonal shift in the hamiltonian in order to compensate for this
    SUBROUTINE UpdateDiagSft()
        IMPLICIT NONE
        INTEGER :: j,k,GrowthSteps

        IF(NoCulls.eq.0) THEN
            GrowRate=(TotWalkers+0.D0)/(TotWalkersOld+0.D0)
        ELSEIF(NoCulls.eq.1) THEN
!GrowRate is the sum of the individual grow rates for each uninterrupted growth sequence, multiplied by the fraction of the cycle which was spent on it
            GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((CullInfo(1,1)+0.D0)/(TotWalkersOld+0.D0))
            GrowRate=GrowRate+(((StepsSft-CullInfo(1,3))+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(1,2)+0.D0))

            NoCulls=0
            CALL IAZZERO(CullInfo,30)
        ELSE
            GrowRate=((CullInfo(1,3)+0.D0)/(StepsSft+0.D0))*((CullInfo(1,1)+0.D0)/(TotWalkersOld+0.D0))
            do j=2,NoCulls
    
                GrowthSteps=CullInfo(j,3)
                do k=1,j-1
                    GrowthSteps=GrowthSteps-CullInfo(k,3)
                enddo

                GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((CullInfo(j,1)+0.D0)/(CullInfo(j-1,2)+0.D0))

            enddo

            GrowthSteps=StepsSft
            do k=1,NoCulls
                GrowthSteps=GrowthSteps-CullInfo(k,3)
            enddo
            GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(NoCulls,2)+0.D0))

            NoCulls=0
            CALL IAZZERO(CullInfo,30)

        ENDIF
        DiagSft=DiagSft-(log(GrowRate))/(SftDamp*Tau*(StepsSft+0.D0))
!        IF((DiagSft).gt.0.D0) THEN
!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!            STOP
!        ENDIF

    END SUBROUTINE UpdateDiagSft


!This routine cancels out particles of opposing sign on the same determinant.
    SUBROUTINE AnnihilatePairs(TotWalkersNew)
        IMPLICIT NONE
        INTEGER :: TotWalkersNew,j,k,l,DetCurr(NEl),VecSlot,TotWalkersDet
        INTEGER :: DetLT

!First, it is necessary to sort the list of determinants
        CALL SortDets(TotWalkersNew,NewDets(:,1:TotWalkersNew),NEl,NewSign(1:TotWalkersNew),1)

!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to CurrentDets
        j=1
!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
        do k=1,NEl
!DetCurr is the current determinant
            DetCurr(k)=NewDets(k,j)
        enddo
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
                    do k=1,NEl
                        CurrentDets(k,VecSlot)=DetCurr(k)
                    enddo
                    CurrentSign(VecSlot)=.true.
                    VecSlot=VecSlot+1
                enddo
            ELSE
!Negative sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    do k=1,NEl
                        CurrentDets(k,VecSlot)=DetCurr(k)
                    enddo
                    CurrentSign(VecSlot)=.false.
                    VecSlot=VecSlot+1
                enddo
            ENDIF
!Now update the current determinant
            do k=1,NEl
                DetCurr(k)=NewDets(k,j)
            enddo
        enddo
!The new number of residual cancelled walkers is given by one less that VecSlot again.
        TotWalkers=VecSlot-1

        RETURN

    END SUBROUTINE AnnihilatePairs

!This routine calculates the MP1 eigenvector, and uses it as a guide for setting the initial walker configuration
    SUBROUTINE StartWavevector(WaveType)
        USE Calc , only : i_P
        USE System , only : Beta
        USE Integrals , only : nTay
        IMPLICIT NONE
        INTEGER :: ierr,i,j,WaveType,EigenvectorTag=0,k,VecSlot,NoDoublesWalk
        CHARACTER(len=*), PARAMETER :: this_routine='StartWavevector'
        REAL*8 :: TypeChange,SumComp,GrowFactor
        INTEGER :: nStore(6),nExcitMemLen,nJ(NEl),iMaxExcit,nExcitTag=0,iExcit,WalkersOnDet
        INTEGER , ALLOCATABLE :: nExcit(:)
        REAL*8 , ALLOCATABLE :: Eigenvector(:)
        TYPE(HElement) :: rhij,rhjj

        IF((WaveType.eq.1).or.(WaveType.eq.2)) THEN
!If WaveType=1, we want to calculate the MP1 wavevector as our initial configuration, and WaveType=2 is the star wavevector
        
            IF(WaveType.eq.1) THEN
!For MP1 wavefunction, we want TypeChange=1.D0, but star will require the star correlation energy
                TypeChange=1.D0
            ELSEIF(WaveType.eq.2) THEN
!Star energy not yet proparly coded & tested
                STOP 'Star initial wavefunction not yet working'
!                StarWeight=fMCPR3StarNewExcit(FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,nTay,RhoEps,iExcit,iMaxExcit,nWHTay,iLogging,TSym,ECore,DBeta,DLWDB,MP2E)
!                TypeChange=DLWDB
            ENDIF
            

!First, generate all excitations, and store their determianants, and rho matrix elements 
            nStore(1)=0
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
            nExcit(1)=0
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

            ALLOCATE(Eigenvector(iMaxExcit+1),stat=ierr)
            CALL LogMemAlloc('Eigenvector',iMaxExcit+1,8,this_routine,EigenvectorTag,ierr)
            CALL AZZERO(Eigenvector,iMaxExcit+1)

!Also need to store the determinants which each component of the eigenvector refers to...
            ALLOCATE(ExcitStore(NEl,iMaxExcit+1),stat=ierr)
            CALL LogMemAlloc('ExcitStore',(iMaxExcit+1)*NEl,4,this_routine,ExcitStoreTag,ierr)
            CALL IAZZERO(ExcitStore,(iMaxExcit+1)*NEl)
            
            CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)

            i=1
            do j=1,NEl
                ExcitStore(j,i)=FDet(j)
            enddo
            Eigenvector(i)=1.D0
            SumComp=1.D0
            
            IF(nTay(2).ne.3) THEN
!Need to ensure that Low-Diag is being used to recover MP1 wavefunction, where c_j=rh_ij/(rh_jj - 1)
!Fock-partition-lowdiag is not set - it must be in order to use the given formulation for the MP1 wavefunction
                WRITE(6,*) "FOCK-PARTITION-LOWDIAG is not specified. It must be to use NODALCUTOFF."
                WRITE(6,*) "Resetting all rho integrals to use FOCK-PARTITION-LOWDIAG"
                nTay(2)=3
            ENDIF
            
            do while (.true.)
                CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
                IF(nJ(1).eq.0) EXIT
                i=i+1
                CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
                CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)

!We want the value of rho_jj/rho_ii
                rhjj=rhjj/rhii
                do j=1,NEl
                    ExcitStore(j,i)=nJ(j)
                enddo
                Eigenvector(i)=(rhij%v)/((rhjj%v)-TypeChange)
                SumComp=SumComp+Eigenvector(i)
            enddo

            NoComps=i

            WRITE(6,*) "Number of components found for new starting wavevector: ",NoComps

            DEALLOCATE(nExcit)
            CALL LogMemDealloc(this_routine,nExcitTag)

        ENDIF

        GrowFactor=(InitWalkers+0.D0)/SumComp

        WRITE(6,*) "Growth factor of initial wavevector is: ",GrowFactor

!Find actual number of new initial walkers
        NoDoublesWalk=0
        InitWalkers=0
        do i=1,NoComps
            WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
            InitWalkers=InitWalkers+WalkersOnDet
            IF(i.ne.1) THEN
                NoDoublesWalk=NoDoublesWalk+WalkersOnDet
            ENDIF
        enddo

        WRITE(6,*) "New number of initial walkers is: ",InitWalkers
        WRITE(6,*) "Number of walkers on double excitations: ",NoDoublesWalk

!Set the maximum number of walkers allowed
        MaxWalkers=MemoryFac*InitWalkers

!Allocate memory to hold walkers
        ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
        ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
        ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
        ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)
    
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign

        VecSlot=1
!Cycle over components
        do i=1,NoComps
!Cycle over number of initial walkers wanted in each component
            IF(Eigenvector(i).gt.0.D0) THEN
!Walkers in this determinant are positive
                WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
                do j=1,WalkersOnDet
                    do k=1,NEl
                        CurrentDets(k,VecSlot)=ExcitStore(k,i)
                    enddo
                    CurrentSign(VecSlot)=.true.
                    VecSlot=VecSlot+1
                enddo

            ELSE
                WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
                do j=1,WalkersOnDet
                    do k=1,NEl
                        CurrentDets(k,VecSlot)=ExcitStore(k,i)
                    enddo
                    CurrentSign(VecSlot)=.false.
                    VecSlot=VecSlot+1
                enddo
            ENDIF

        enddo

        IF((VecSlot-1).ne.InitWalkers) THEN
            WRITE(6,*) "Problem in assigning particles proportionally to given wavevector..."
            STOP "Problem in assigning particles proportionally to given wavevector..."
        ENDIF

        IF(TNoBirth) THEN
!Allocate memory to hold all initial determinants and their starting populations
            ALLOCATE(InitPops(NoComps),stat=ierr)
            CALL LogMemAlloc('InitPops',NoComps,4,this_routine,InitPopsTag)
                
            WRITE(6,*) ""
            WRITE(6,*) "    Step   Components   WalkersOnDet   Expected"
            WRITE(15,*) "#    Step   Components   WalkersOnDet   Expected"

            do i=1,NoComps
                InitPops(i)=abs(nint(Eigenvector(i)*GrowFactor))
                WRITE(6,"(2I9,2G16.7)") 0,i,InitPops(i),InitPops(i)
                WRITE(15,"(2I9,2G16.7)") 0,1,InitPops(i),InitPops(i)
            enddo
            WRITE(6,*) ""
            WRITE(15,*) ""

!Can deallocate the eigenvector, but we want to keep the ExcitStore in order to compare the populations at a later date
            DEALLOCATE(Eigenvector)
            CALL LogMemDealloc(this_routine,EigenvectorTag)

        ELSE

            DEALLOCATE(Eigenvector)
            CALL LogMemDealloc(this_routine,EigenvectorTag)
            DEALLOCATE(ExcitStore)
            CALL LogMemDealloc(this_routine,ExcitStoreTag)

        ENDIF

        RETURN

    END SUBROUTINE StartWavevector


!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalc()
        USE DetCalc , only : NDet
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot
        TYPE(HElement) :: rh
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMC'


        IF(TStartMP1) THEN
!Start the initial distribution off at the distribution of the MP1 eigenvector

            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."
            CALL StartWavevector(1)

        ELSEIF(TReadPops) THEN

!Set the maximum number of walkers allowed
            MaxWalkers=MemoryFac*InitWalkers

!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

            CurrentDets=>WalkVecDets
            CurrentSign=>WalkVecSign
            NewDets=>WalkVec2Dets
            NewSign=>WalkVec2Sign

            IF((ABS(ScaleWalkers-1.D0)).lt.1.D-8) THEN
!Read in walker positions
                do i=1,InitWalkers
                    READ(17,*) CurrentDets(:,i),CurrentSign(i)
                enddo
            ELSE
!Read in walker positions - we will scale these later...
                do i=1,InitWalkers
                    READ(17,*) NewDets(:,i),NewSign(i)
                enddo
                WRITE(6,*) "Scaling number of walkers by: ",ScaleWalkers
                ReadWalkers=InitWalkers
                InitWalkers=0
!First, count the total number of initial walkers on each determinant - sort into list
                CALL SortDets(ReadWalkers,NewDets(:,1:ReadWalkers),NEl,NewSign(1:ReadWalkers),1)

                j=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                do k=1,NEl
!DetCurr is the current determinant
                    DetCurr(k)=NewDets(k,j)
                enddo

                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((DetLT(NewDets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(NewSign(j)) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
!Now update the current determinant
                    do k=1,NEl
                        DetCurr(k)=NewDets(k,j)
                    enddo
!Count total number of initial walkers
                    InitWalkers=InitWalkers+abs(nint((TotWalkersDet+0.D0)*ScaleWalkers))
                enddo
                WRITE(6,*) "Total number of walkers is now: ",InitWalkers
!Set the new maximum number of walkers allowed
                MaxWalkers=MemoryFac*InitWalkers

!Deallocate old memory block for WalkVec
                DEALLOCATE(WalkVecDets)
                CALL LogMemDealloc(this_routine,WalkVecDetsTag)
                DEALLOCATE(WalkVecSign)
                CALL LogMemDealloc(this_routine,WalkVecSignTag)

!Allocate memory to hold new maximum number of walkers
                ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
                ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
                
                CurrentDets=>WalkVecDets
                CurrentSign=>WalkVecSign

!Now multiply them up...
                j=1
                VecSlot=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                do k=1,NEl
                    DetCurr(k)=NewDets(k,j)
                enddo
!DetCurr is the current determinant
                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((DetLT(NewDets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(NewSign(j)) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
!Now multiply up the number of walkers, and insert into CurrentDets
                    TotWalkersDet=nint((TotWalkersDet+0.D0)*ScaleWalkers)
                    IF(TotWalkersDet.gt.0) THEN
                        do l=1,abs(TotWalkersDet)
                            do k=1,NEl
                                CurrentDets(k,VecSlot)=DetCurr(k)
                            enddo
                            CurrentSign(VecSlot)=.true.
                            VecSlot=VecSlot+1
                        enddo
                    ELSE
                        do l=1,abs(TotWalkersDet)
                            do k=1,NEl
                                CurrentDets(k,VecSlot)=DetCurr(k)
                            enddo
                            CurrentSign(VecSlot)=.false.
                            VecSlot=VecSlot+1
                        enddo
                    ENDIF
                    do k=1,NEl
                        DetCurr(k)=NewDets(k,j)
                    enddo
                enddo
                IF((VecSlot-1).ne.InitWalkers) THEN
                    WRITE(6,*) "Problem scaling up walker number - exiting..."
                    STOP 'Problem scaling up walker number - exiting...'
                ENDIF

!Now deallocate and reallocate WalkVec2 with correct number of total walkers
                DEALLOCATE(WalkVec2Dets)
                CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
                DEALLOCATE(WalkVec2Sign)
                CALL LogMemDealloc(this_routine,WalkVec2SignTag)
                ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
                ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

                NewDets=>WalkVec2Dets
                NewSign=>WalkVec2Sign

            ENDIF

!End of reading in POPSFILE
            CLOSE(17)

        ELSE
!If not reading in from POPSFILE, then we need to initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkers=MemoryFac*InitWalkers

!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

            CurrentDets=>WalkVecDets
            CurrentSign=>WalkVecSign
            NewDets=>WalkVec2Dets
            NewSign=>WalkVec2Sign

            do j=1,InitWalkers
                do k=1,NEl
                    CurrentDets(k,j)=FDet(k)
                enddo
                CurrentSign(j)=.true.
            enddo

            IF(TNoBirth) THEN
                
                WRITE(6,*) ""
                WRITE(6,*) "    Step   Components   WalkersOnDet   Expected"
                WRITE(15,*) "#    Step   Components   WalkersOnDet   Expected"

                ALLOCATE(InitPops(1),stat=ierr)
                CALL LogMemAlloc('InitPops',1,4,this_routine,InitPopsTag)
                ALLOCATE(ExcitStore(NEl,1),stat=ierr)
                CALL LogMemAlloc('ExcitStore',NEl,4,this_routine,ExcitStoreTag)

                NoComps=1
                do j=1,NEl
                    ExcitStore(j,1)=FDet(j)
                enddo
                InitPops(1)=InitWalkers
                WRITE(6,"(2I9,2G16.7)") 0,1,InitWalkers,InitWalkers
                WRITE(15,"(2I9,2G16.7)") 0,1,InitWalkers,InitWalkers

            ELSEIF(TDetPops) THEN
                
                SizeofSpace=NDET
                WRITE(6,*) "Size of space is: ", SizeOfSpace

                ALLOCATE(PopsVec(SizeofSpace),stat=ierr)
                CALL LogMemAlloc('PopsVec',SizeofSpace,8,this_routine,PopsVecTag)
                CALL AZZERO(PopsVec,SizeofSpace)
                ALLOCATE(TransMat(SizeofSpace,SizeofSpace),stat=ierr)
                CALL LogMemAlloc('TransMat',SizeofSpace**2,8,this_routine,TransMatTag)
                CALL AZZERO(TransMat,SizeofSpace**2)

                IF(DetLT(NMRKS(:,1),FDet,NEl).ne.0) THEN
                    WRITE(6,*) "Problem with NMRKS"
                    STOP "Problem with NMRKS"
                ENDIF
                
                do j=1,SizeOfSpace
                    rh=GetHElement2(FDet,NMRKS(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
                    WRITE(67,"(F10.5,$)") rh%v
                enddo
                WRITE(67,*) ""
                WRITE(67,*) ""

            ENDIF

        ENDIF

!TotWalkers contains the number of current walkers at each step
        TotWalkers=InitWalkers
        TotWalkersOld=InitWalkers

        CALL IAZZERO(CullInfo,30)
        NoCulls=0

        IF(TNodalCutoff) CALL CalcNodalSurface()

        RETURN

    END SUBROUTINE InitFCIMCCalc

!This routine calculates the normalisation for the MP1 wavefunction. This is needed if a nodal structure is being applied, and can also calculate the number of determinants which are being constrained.
    SUBROUTINE CalcNodalSurface()
        USE Calc , only : i_P
        USE System , only : Beta
        USE Integrals , only : nTay
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,Doubs,FixedSign
        REAL*8 :: EigenComp
        CHARACTER(len=*), PARAMETER :: this_routine='CalcNodalSurface'
        INTEGER :: nStore(6),nExcitMemLen,nJ(NEl),iMaxExcit,nExcitTag=0,iExcit
        INTEGER , ALLOCATABLE :: nExcit(:)
        TYPE(HElement) :: Fj,Hamij!,rhij,rhjj

        WRITE(6,"(A,F19.9)") "Calculating the nodal structure of the MP1 wavefunction with a normalised cutoff of ",NodalCutoff

!        IF(nTay(2).ne.3) THEN
!This is no longer needed since the MP1 components are calculated exactly
!Fock-partition-lowdiag is not set - it must be in order to use the given formulation for the MP1 wavefunction
!            WRITE(6,*) "FOCK-PARTITION-LOWDIAG is not specified. It must be to use NODALCUTOFF."
!            WRITE(6,*) "Resetting all rho integrals to use FOCK-PARTITION-LOWDIAG"
!            nTay(2)=3
!        ENDIF

!First, generate all excitations, and store their determianants, and rho matrix elements 
        nStore(1)=0
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,2)
        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
        nExcit(1)=0
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,2)

!Calculate F0
        CALL GetH0Element(FDet,NEl,Arr,nBasis,ECore,FZero)
!        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)

        Doubs=0
!The HF determinant has a component 1         
        MPNorm=1.D0
        do while (.true.)
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,2)
            IF(nJ(1).eq.0) EXIT
            Doubs=Doubs+1
            Hamij=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fj)
!            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
!            CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)

!We want the value of rho_jj/rho_ii
!            rhjj=rhjj/rhii
!EigenComp is now the component of the MP1 wavefunction for nJ, but unnormalised
!            EigenComp=(rhij%v)/((rhjj%v)-1.D0)
            EigenComp=(Hamij%v)/(Fj%v-FZero%v)
            MPNorm=MPNorm+(EigenComp**2)
        enddo

        WRITE(6,"(A,I15)") "Number of double excitations found in MP Wavevector: ",Doubs

!Find the normalisation for the MP1 wavevector
        MPNorm=SQRT(MPNorm)

!Reset the excitation generator, so that we can run through the excitations again and calculate the number of components of the MP1 wavefunction which will be sign-constrained
        CALL ResetExit2(FDet,NEl,G1,nBasis,nBasisMax,nExcit,0)

        FixedSign=0     !FixedSign is the number of determinants which are fixed in sign with the given value of the NodalCutoff

        IF((1.D0/MPNorm).gt.abs(NodalCutoff)) THEN
!This indicates that the HF is included in the fixed sign approximation, and so will always be constrained to have net positive particles.
            FixedSign=FixedSign+1
            WRITE(6,*) "Hartree-Fock determinant constrained to always have net positive sign"
        ELSE
            WRITE(6,"(A,F17.9)") "Hartree-Fock determinant NOT constrained to have net positive sign, since normalised component is only: ", 1.D0/MPNorm
        ENDIF

        do while(.true.)
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,2)
            IF(nJ(1).eq.0) EXIT
            Hamij=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fj)
!            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
!            CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)

!We want the value of rho_jj/rho_ii
!            rhjj=rhjj/rhii
!EigenComp is now the component of the MP1 wavefunction for nJ, but now normalised
!            EigenComp=abs(((rhij%v)/((rhjj%v)-1.D0))/MPNorm)
            EigenComp=abs(((Hamij%v)/(Fj%v-FZero%v))/MPNorm)
            IF(EigenComp.gt.abs(NodalCutoff)) THEN
!Increase the counter of number of determinants with fixed sign
                FixedSign=FixedSign+1
            ENDIF

        enddo

        WRITE(6,*) "Total number of determinants constrained by fixed sign: ", FixedSign

        DEALLOCATE(nExcit)
        CALL LogMemDealloc(this_routine,nExcitTag)

        RETURN

    END SUBROUTINE CalcNodalSurface

!This function tells us whether we want to diffuse from DetCurr to nJ
    SUBROUTINE AttemptDiffuse(DetCurr,nJ,Prob,IC,WSign,KeepOrig,CreateAtI,CreateAtJ)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,i,CreateAtI,CreateAtJ
        REAL*8 :: Prob,Ran2,rat
        LOGICAL :: WSign,KeepOrig
        TYPE(HElement) :: rh

!Calculate off-diagonal hamiltonian matrix element between determinants
        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        
!rat is the probability of diffusing to nJ
        rat=Tau*Lambda*abs(rh%v)/Prob

        IF(rat.gt.1.D0) CALL STOPGM("AttemptDiffuse","*** Probability > 1 to diffuse.")

        IF(rat.gt.Ran2(Seed)) THEN
            IF(TExtraPartDiff) THEN
!We want to perform the non-total number conserving diffusion matrix - anti-diffusion creates 2 new particles
                IF(real(rh%v).gt.0.D0) THEN
!Perform anti-diffusion - particle number will increase
                    IF(WSign) THEN
!We have a positive walker
                        KeepOrig=.true.
                        CreateAtI=1
                        CreateAtJ=-1
                    ELSE
!We have a negative walker
                        KeepOrig=.true.
                        CreateAtI=-1
                        CreateAtJ=1
                    ENDIF
                ELSE
!Perform diffusion - particle number will remain constant
                    IF(WSign) THEN
                        KeepOrig=.false.    !Particle is annihilated by newly created opposing signed particle
                        CreateAtI=0
                        CreateAtJ=1
                    ELSE
                        KeepOrig=.false.
                        CreateAtI=0
                        CreateAtJ=-1
                    ENDIF
                ENDIF
            ELSE
!We are performing a number-conserving diffusion process - this will have a different diagonal birth/death unbiasing probability
                IF(real(rh%v).gt.0.D0) THEN
!Perform anti-diffusion, but conserve total particle number
                    IF(WSign) THEN
                        KeepOrig=.false.
                        CreateAtI=0
                        CreateAtJ=-1
                    ELSE
                        KeepOrig=.false.
                        CreateAtI=0
                        CreateAtJ=1
                    ENDIF
                ELSE
!Perform diffusion - this should conserve particle number, and be the same as the other diffusion process
                    IF(WSign) THEN
                        KeepOrig=.false.
                        CreateAtI=0
                        CreateAtJ=1
                    ELSE
                        KeepOrig=.false.
                        CreateAtI=0
                        CreateAtJ=-1
                    ENDIF
                ENDIF
            ENDIF
        ELSE
!No diffusion will occur...
            KeepOrig=.true.        !We don't want to copy it accross, because it still can die - wait to see if it dies before copying accross
            CreateAtI=0
            CreateAtJ=0
        ENDIF

        RETURN

    END SUBROUTINE AttemptDiffuse

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
!Kik returns the a value which can be used to unbias the diffusion at the birth/death stage
    INTEGER FUNCTION AttemptCreate(DetCurr,WSign,nJ,Prob,IC,Kik)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate
        LOGICAL :: WSign
        REAL*8 :: Prob,Ran2,rat,Kik
        TYPE(HElement) :: rh

!Calculate off diagonal hamiltonian matrix element between determinants
        IF(TNoBirth) THEN
            rh=HElement(0.D0)
        ELSE
            rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        ENDIF

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
        IF(TDiffuse) THEN
            IF(TExtraPartDiff) THEN
                Kik=Tau*rh%v/Prob
            ELSE
                Kik=Tau*abs(rh%v)/Prob
            ENDIF
            rat=(1.D0-Lambda)*abs(Kik)
        ELSE
            rat=Tau*abs(rh%v)/Prob
        ENDIF

!        IF(rat.lt.0.D0) THEN
!            CALL STOPGM("AttemptCreate","*** Probability < 0 to create child.")
!        ENDIF

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

        IF(Tau.lt.0.D0) THEN
!If tau is negative, we are going back in time, and so will actually create antiparticles - flip sign again...
            AttemptCreate=-AttemptCreate
        ENDIF

        IF(TDetPops) THEN
!Here, we want to record the details of every spawning connection

            StoreNumTo=0
            StoreNumFrom=0
            do i=1,SizeOfSpace
                IF((DetLT(nJ,NMRKS(:,i),NEl).eq.0).and.(StoreNumTo.eq.0)) THEN
!Found position of determinant
                    StoreNumTo=i
                ENDIF
                IF((DetLT(DetCurr,NMRKS(:,i),NEl).eq.0).and.(StoreNumFrom.eq.0)) THEN
                    StoreNumFrom=i
                ENDIF
                IF((StoreNumTo.ne.0).and.(StoreNumFrom.ne.0)) EXIT
            enddo

            IF(AttemptCreate.lt.0) THEN
!If creating a negative particle, reduce the connection
                TransMat(StoreNumFrom,StoreNumTo)=TransMat(StoreNumFrom,StoreNumTo)-(AttemptCreate*0.00001)
                TransMat(StoreNumTo,StoreNumFrom)=TransMat(StoreNumTo,StoreNumFrom)-(AttemptCreate*0.00001)
            ELSE
!If creating a positive particle, increase connection
                TransMat(StoreNumFrom,StoreNumTo)=TransMat(StoreNumFrom,StoreNumTo)+(AttemptCreate*0.00001)
                TransMat(StoreNumTo,StoreNumFrom)=TransMat(StoreNumTo,StoreNumFrom)+(AttemptCreate*0.00001)
            ENDIF

        ENDIF

        RETURN

    END FUNCTION AttemptCreate

!This function tells us whether we should kill the particle at determinant DetCurr
!If also diffusing, then we need to know the probability with which we have spawned. This will reduce the death probability.
!The function allows multiple births(if +ve shift) or deaths from the same particle.
!The returned number is the number of deaths if positive, and the number of births if negative.
    INTEGER FUNCTION AttemptDie(DetCurr,Kik,nExcitMemLen,nExcit,nStore)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),DetLT,iKill,nExcitMemLen,nStore(6),nExcit(nExcitMemLen)
        INTEGER :: nJ(NEl),IC,iMaxExcit,ierr,nExcitTag=0
        CHARACTER(len=*) , PARAMETER :: this_routine='AttemptDie'
        TYPE(HElement) :: rh,rhij
        REAL*8 :: Ran2,rat,Kik

!Test if determinant is FDet - in a strongly single-configuration problem, this will save time
        IF(DetLT(DetCurr,FDet,NEl).eq.0) THEN
            rh=0.D0
        ELSE
!Calculate the diagonal hamiltonian matrix element for the determinant
            rh=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!Subtract from the diagonal the value of the lowest hamiltonian matrix element
            rh=rh-Hii
        ENDIF

        IF(TDiffuse) THEN
!If also diffusing, then the probability of dying must be modified, since the diagonal elements have been altered
            IF(TExtraPartDiff) THEN

                IF(TFullUnbias) THEN


!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
                    
                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
                    
                    Kik=0.D0    !If we are fully unbiasing, then the unbiasing factor is reset and recalculated from all excitations of DetCurr
                    do while(.true.)
                        CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
                        IF(nJ(1).eq.0) EXIT
                        rhij=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
                        Kik=Kik+(rhij%v)
                    enddo

                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)

                    Kik=Kik*Tau    !Unbias with the tau*sum of connected elements
                ENDIF

                rat=(Tau*((rh%v)-DiagSft))+(Lambda*Kik)     !This is now the probability with the correct unbiasing

            ELSE
                IF(TFullUnbias) THEN

                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)
                    
                    Kik=0.D0    !If we are fully unbiasing, then the unbiasing factor is reset and recalculated from all excitations of DetCurr
                    do while(.true.)
                        CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
!                        CALL WRITEDET(6,nJ,NEl,.true.)
                        IF(nJ(1).eq.0) EXIT
                        rhij=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
                        Kik=Kik+abs(rhij%v)
                    enddo

                    CALL ResetExIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,0)

                    Kik=Kik*Tau    !Unbias with the tau*sum of connected elements
                ENDIF

                rat=(Tau*((rh%v)-DiagSft))-(Lambda*Kik)     !This is now the probability with the correct unbiasing

            ENDIF

        ELSE
!Subtract the current value of the shift and multiply by tau
            rat=Tau*((rh%v)-DiagSft)
        ENDIF

!        IF(rat.gt.1.D0) THEN
!If probs of dying is greater than one, reduce tau
!            CALL STOPGM("AttemptDie","*** Death probability > 1. *** Tau too large")
!        ENDIF

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
!        IF(AttemptDie.le.-1) WRITE(6,*) Iter,AttemptDie,rat

        RETURN

    END FUNCTION AttemptDie

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
