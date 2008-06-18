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
    USE Determinants , only : FDet,GetHElement2
    USE DetCalc , only : NMRKS
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE Logging , only : TPopsFile,TCalcWavevector,WavevectorPrint,TDetPops
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    IMPLICIT NONE
    SAVE

    INTEGER , POINTER :: WalkVecDets(:,:),WalkVec2Dets(:,:)
    LOGICAL , POINTER :: WalkVecSign(:),WalkVec2Sign(:)
    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0

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

    REAL*8 :: GrowRate,DieRat

    TYPE(HElement) :: Hii

    contains

    SUBROUTINE FciMC(Weight,Energyxw)
        Use MCDets, only : MCDetsCalc
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j,iSub,WalkOnDet,DetLT,DetCurr(NEl),ExpectedDets
        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'
        TYPE(HElement) :: Hamii

        if(NMCyc.lt.0.or.TMCDets) then

           CALL MCDetsCalc(FDet,G_VMC_Seed,-NMCyc,Tau,SftDamp,10000*InitWalkers,InitWalkers,StepsSft,DiagSft)
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
            WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers"
            WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers"
!TotWalkersOld is the number of walkers last time the shift was changed
            IF(TReadPops) THEN
                WRITE(6,"(I9,G16.7,I9,G16.7,I9)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                WRITE(15,"(I9,G16.7,I9,G16.7,I9)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
            ELSE
                WRITE(6,"(I9,G16.7,I9,G16.7,I9)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                WRITE(15,"(I9,G16.7,I9,G16.7,I9)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
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
                            IF(DetLT(DetCurr,WalkVecDets(:,j),NEl).eq.0) THEN
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
                        WRITE(15,"(I9,G16.7,I9,G16.7,I9)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                        WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                    ELSE
                        IF(Tau.gt.0.D0) THEN
                            WRITE(15,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                        ELSE
                            WRITE(15,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
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
                WRITE(16,*) WalkVecDets(:,i),WalkVecSign(i)
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
        INTEGER :: nJ(NEl),ierr,nExcitTag=0,IC,Child,iSubCyc,TotWalkersNew,iCount!,NoDie
        REAL*8 :: Prob,rat,Kik
        INTEGER , ALLOCATABLE :: nExcit(:)
        LOGICAL :: TDiffused        !Indicates whether a particle has diffused or not
        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
        LOGICAL :: WSign
        CHARACTER(len=*), PARAMETER :: this_routine='PerformFCIMCyc'
        
        CALL TISET('MCyc',iSubCyc)
        
!VecSlot indicates the next free position in WalkVec2
        VecSlot=1
!        NoDie=0

        do j=1,TotWalkers
!j runs through all current walkers
            do k=1,NEl
                DetCurr(k)=WalkVecDets(k,j)
            enddo
            TDiffused=.false.

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

                    Child=AttemptCreate(DetCurr,WalkVecSign(j),nJ,Prob,IC,Kik)
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
                            WalkVec2Dets(k,VecSlot)=nJ(k)
                        enddo
                        WalkVec2Sign(VecSlot)=WSign
                        VecSlot=VecSlot+1
                    enddo

                enddo

            ELSE
!Run through all possible excitations of each walker

                do while(.true.)
                    CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
                    IF(nJ(1).eq.0) EXIT

                    Child=AttemptCreate(DetCurr,WalkVecSign(j),nJ,1.D0,IC,Kik)
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
                            WalkVec2Dets(k,VecSlot)=nJ(k)
                        enddo
                        WalkVec2Sign(VecSlot)=WSign
                        VecSlot=VecSlot+1
                    enddo

                enddo

            ENDIF

            IF(TDiffuse) THEN
!Next look at possibility of diffusion to another determinant
                CALL GenRandSymExcitIt3(DetCurr,nExcit,nJ,Seed,IC,0,Prob,iCount)
                IF(AttemptDiffuse(DetCurr,nJ,Prob,IC)) THEN
!Walker wants to diffuse to nJ from DetCurr, and keep the same sign - copy across nJ, rather than the original walker
                    do k=1,NEl
                        WalkVec2Dets(k,VecSlot)=nJ(k)
                    enddo
                    WalkVec2Sign(VecSlot)=WalkVecSign(j)
                    VecSlot=VecSlot+1
                    TDiffused=.true.
!                    NoDie=NoDie+1
                ENDIF
            ENDIF

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
            iDie=AttemptDie(DetCurr,Kik)
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
            IF((iDie.le.0).AND.(.NOT.TDiffused)) THEN
!This indicates that the particle is spared and we are not diffusing...copy him across to WalkVec2
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle
!NoDie actually counts number of particles spared - later it is transformed into the number that actually die
!                NoDie=NoDie+1+abs(iDie)
                do l=1,(abs(iDie)+1)
                    do k=1,NEl
                        WalkVec2Dets(k,VecSlot)=DetCurr(k)
                    enddo
                    WalkVec2Sign(VecSlot)=WalkVecSign(j)
                    VecSlot=VecSlot+1
                enddo
            ELSEIF((iDie.gt.1).AND.(.NOT.TDiffused)) THEN
!This indicates that particles on DetCurr want to be killed. The first kill will simply be performed by not copying accross the original particle.
!Therefore, if iDie = 1, then we can simply ignore it.
!However, after that anti-particles will need to be created on the same determinant.
!                NoDie=NoDie-(iDie-1)
                do l=1,(iDie-1)
                    do k=1,NEl
                        WalkVec2Dets(k,VecSlot)=DetCurr(k)
                    enddo
                    IF(WalkVecSign(j)) THEN
!Copy accross new anti-particles
                        WalkVec2Sign(VecSlot)=.FALSE.
                    ELSE
                        WalkVec2Sign(VecSlot)=.TRUE.
                    ENDIF
                    VecSlot=VecSlot+1
                enddo
            ELSEIF((iDie.gt.0).AND.TDiffused) THEN
!This means that the particle is destroyed, as well as possibly more, but the actual particle has already diffused to a different determinant (nJ).
!We need to destroy on the original determinant by creating an opposing signed particle on the original determinant.
                do l=1,iDie
                    do k=1,NEl
                        WalkVec2Dets(k,VecSlot)=DetCurr(k)
                    enddo
                    IF(WalkVecSign(j)) THEN
                        WalkVec2Sign(VecSlot)=.FALSE.
                    ELSE
                        WalkVec2Sign(VecSlot)=.TRUE.
                    ENDIF
                    VecSlot=VecSlot+1
                
!Since we have already indicated that the particle is spared since it is copied accross, just in a different determinant, we have to subtract to ge the correct number of particles spared.
!                    NoDie=NoDie-1
                enddo
            ELSEIF((iDie.lt.0).AND.TDiffused) THEN
!Here, the particle is is recreated. However, the original particle has moved to a new determinant, and so additional particles need to be created on the original determinant, which are identical to the original.
                do l=1,(abs(iDie))
                    do k=1,NEl
                        WalkVec2Dets(k,VecSlot)=DetCurr(k)
                    enddo
                    WalkVec2Sign(VecSlot)=WalkVecSign(j)
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

!        NoDie=TotWalkers-NoDie
!        DieRat=(NoDie+0.D0)/(TotWalkers+0.D0)
                
!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
        TotWalkersNew=VecSlot-1
        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
        ENDIF

!This routine now cancels down the particles with opposing sign on each determinant
!This routine does not necessarily need to be called every Iter, but it does at the moment, since it is the only way to 
!transfer the residual particles back onto WalkVec

        CALL AnnihilatePairs(TotWalkersNew)
!        WRITE(6,*) "Number of annihilated particles= ",TotWalkersNew-TotWalkers
        
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
                    WalkVecDets(i,Chosen)=WalkVecDets(i,TotWalkers)
                enddo
                WalkVecSign(Chosen)=WalkVecSign(TotWalkers)

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
                    WalkVecDets(j,VecSlot)=WalkVecDets(j,i)
                enddo
                WalkVecSign(VecSlot)=WalkVecSign(i)

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
        CALL SortDets(TotWalkersNew,WalkVec2Dets(:,1:TotWalkersNew),NEl,WalkVec2Sign(1:TotWalkersNew),1)

!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to WalkVec
        j=1
!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
        do k=1,NEl
!DetCurr is the current determinant
            DetCurr(k)=WalkVec2Dets(k,j)
        enddo
        VecSlot=1

        do while(j.le.TotWalkersNew)
!Loop over all walkers
            TotWalkersDet=0
            do while ((DetLT(WalkVec2Dets(:,j),DetCurr,NEl).eq.0).and.(j.le.TotWalkersNew))
!Loop over all walkers on DetCurr and count residual number after cancelling
                IF(WalkVec2Sign(j)) THEN
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
                        WalkVecDets(k,VecSlot)=DetCurr(k)
                    enddo
                    WalkVecSign(VecSlot)=.true.
                    VecSlot=VecSlot+1
                enddo
            ELSE
!Negative sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    do k=1,NEl
                        WalkVecDets(k,VecSlot)=DetCurr(k)
                    enddo
                    WalkVecSign(VecSlot)=.false.
                    VecSlot=VecSlot+1
                enddo
            ENDIF
!Now update the current determinant
            do k=1,NEl
                DetCurr(k)=WalkVec2Dets(k,j)
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
        TYPE(HElement) :: rhii,rhij,rhjj

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
!Need to ensure that Low-Diag is being used to recover MP1 wavefunction, where c_j=rh_ij/(rh_jj - 1)
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
    
        VecSlot=1
!Cycle over components
        do i=1,NoComps
!Cycle over number of initial walkers wanted in each component
            IF(Eigenvector(i).gt.0.D0) THEN
!Walkers in this determinant are positive
                WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
                do j=1,WalkersOnDet
                    do k=1,NEl
                        WalkVecDets(k,VecSlot)=ExcitStore(k,i)
                    enddo
                    WalkVecSign(VecSlot)=.true.
                    VecSlot=VecSlot+1
                enddo

            ELSE
                WalkersOnDet=abs(nint(Eigenvector(i)*GrowFactor))
                do j=1,WalkersOnDet
                    do k=1,NEl
                        WalkVecDets(k,VecSlot)=ExcitStore(k,i)
                    enddo
                    WalkVecSign(VecSlot)=.false.
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

            IF((ABS(ScaleWalkers-1.D0)).lt.1.D-8) THEN
!Read in walker positions
                do i=1,InitWalkers
                    READ(17,*) WalkVecDets(:,i),WalkVecSign(i)
                enddo
            ELSE
!Read in walker positions - we will scale these later...
                do i=1,InitWalkers
                    READ(17,*) WalkVec2Dets(:,i),WalkVec2Sign(i)
                enddo
                WRITE(6,*) "Scaling number of walkers by: ",ScaleWalkers
                ReadWalkers=InitWalkers
                InitWalkers=0
!First, count the total number of initial walkers on each determinant - sort into list
                CALL SortDets(ReadWalkers,WalkVec2Dets(:,1:ReadWalkers),NEl,WalkVec2Sign(1:ReadWalkers),1)

                j=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                do k=1,NEl
!DetCurr is the current determinant
                    DetCurr(k)=WalkVec2Dets(k,j)
                enddo

                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((DetLT(WalkVec2Dets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(WalkVec2Sign(j)) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
!Now update the current determinant
                    do k=1,NEl
                        DetCurr(k)=WalkVec2Dets(k,j)
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

!Now multiply them up...
                j=1
                VecSlot=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                do k=1,NEl
                    DetCurr(k)=WalkVec2Dets(k,j)
                enddo
!DetCurr is the current determinant
                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((DetLT(WalkVec2Dets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(WalkVec2Sign(j)) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
!Now multiply up the number of walkers, and insert into WalkVec
                    TotWalkersDet=nint((TotWalkersDet+0.D0)*ScaleWalkers)
                    IF(TotWalkersDet.gt.0) THEN
                        do l=1,abs(TotWalkersDet)
                            do k=1,NEl
                                WalkVecDets(k,VecSlot)=DetCurr(k)
                            enddo
                            WalkVecSign(VecSlot)=.true.
                            VecSlot=VecSlot+1
                        enddo
                    ELSE
                        do l=1,abs(TotWalkersDet)
                            do k=1,NEl
                                WalkVecDets(k,VecSlot)=DetCurr(k)
                            enddo
                            WalkVecSign(VecSlot)=.false.
                            VecSlot=VecSlot+1
                        enddo
                    ENDIF
                    do k=1,NEl
                        DetCurr(k)=WalkVec2Dets(k,j)
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

            do j=1,InitWalkers
                do k=1,NEl
                    WalkVecDets(k,j)=FDet(k)
                enddo
                WalkVecSign(j)=.true.
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


        RETURN

    END SUBROUTINE InitFCIMCCalc

!This function tells us whether we want to diffuse from DetCurr to nJ
    LOGICAL FUNCTION AttemptDiffuse(DetCurr,nJ,Prob,IC)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,i
        REAL*8 Prob,Ran2,rat
        TYPE(HElement) :: rh

!Calculate off-diagonal hamiltonian matrix element between determinants
        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        
!rat is the probability of diffusing to nJ
        rat=Tau*Lambda*abs(rh%v)/Prob

        IF(rat.gt.1.D0) CALL STOPGM("AttemptDiffuse","*** Probability > 1 to diffuse.")

        IF(rat.gt.Ran2(Seed)) THEN
            AttemptDiffuse=.true.
        ELSE
            AttemptDiffuse=.false.
        ENDIF

        RETURN

    END FUNCTION AttemptDiffuse

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
!Kik returns the probability of creating the child without any diffusion term.
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
        Kik=abs(Tau)*abs(rh%v)/Prob

        IF(TDiffuse) THEN
            rat=(1.D0-Lambda)*Kik
        ELSE
            rat=Kik
        ENDIF
!        IF(rat.lt.0.D0) THEN
!            CALL STOPGM("AttemptCreate","*** Probability < 0 to create child.")
!        ENDIF

!If probability is > 1, then we can just create multiple children at the chosen determinant
        ExtraCreate=INT(rat)
        rat=rat-DREAL(ExtraCreate)


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
    INTEGER FUNCTION AttemptDie(DetCurr,Kik)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),DetLT,iKill
        TYPE(HElement) :: rh
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
!If also diffusing, then the probability of dying must be reduced, since the probability of annihilation has been increased by it.
            rat=(Tau*((rh%v)-DiagSft))-(Lambda*Kik)
        ELSE
!Subtract the current value of the shift and multiply by tau
            rat=Tau*((rh%v)-DiagSft)
        ENDIF

!        IF(rat.gt.1.D0) THEN
!If probs of dying is greater than one, reduce tau
!            CALL STOPGM("AttemptDie","*** Death probability > 1. *** Tau too large")
!        ENDIF

        iKill=INT(rat)
        rat=rat-DREAL(iKill)

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
