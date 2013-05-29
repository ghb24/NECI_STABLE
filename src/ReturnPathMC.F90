! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
MODULE ReturnPathMCMod
!    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr,nMsh,Beta
!    use CalcData , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
!    use CalcData , only : GrowMaxFactor,CullFactor,i_P
!    use CalcData , only : PRet  !This is the probability of generating the return determinant to spawn to
!    use CalcData , only : CLMax   !This is the maximum allowed chain length
!    use CalcData , only : TRhoElems     !This tells us to use rho elements rather than H-elements
!    use CalcData , only : NEquilSteps
!    USE Determinants , only : FDet,GetHElement2
!    USE DetCalc , only : NMRKS
!    use IntegralsData , only : fck,NMax,UMat,nTay
!    USE global_utilities
!    use constants, only: dp
    use constants, only: dp
!    IMPLICIT NONE
!    SAVE
!
!    INTEGER , ALLOCATABLE :: HFDet(:)           !This is the HF determinant - do not want to use FDet as this is pre-frozen
!    INTEGER :: HFDetTag=0
!
!    
!    TYPE ExcitGenerator                         !Derived type for an excitation generator
!        INTEGER , ALLOCATABLE :: ExcitData(:)   !This is the actual excitation generator data
!        INTEGER :: nStore(6)                    !This holds the nStore for the generator
!        INTEGER :: ExcitLen                     !This is the length of the excitation generator
!        LOGICAL :: ForCurrentDet                !This is a logical which tells us whether the 
                                        !excitation generator stored is for use with the current det.
!    END TYPE
!
!    TYPE Part                   !This is the type for a spawning particle
!        INTEGER , ALLOCATABLE :: Det(:)     !This is the current determinant of the particle - (NEl)
!        INTEGER :: ChainLength  !This is the current length of the chain back to HF (If at HF, ChainLength = 0)
!        INTEGER , ALLOCATABLE :: IC0(:)   !This is the excitation level of the particle w.r.t. HF determinant for 
!now (IC0(ChainLength)) and all its histories - (CLMax)
!        real(dp) , ALLOCATABLE :: Kii(:)    !This is the diagonal K-matrix element: Hii-H00 for now (Kii(CLMax)) 
!and all its histories - (CLMax)
!        INTEGER , ALLOCATABLE :: HistExcit(:)             !HistExcit(i) is the number of excitations between the 
!determinant at History(:,i) and History(:,i-1) - size(CLMax)
!!The number of excitations between the current determinant the one it was spawned from is given in HistExcit(ChainLength)
!        real(dp) :: Hi0           !This is the off-diagonal hamiltonian matrix element between the determinant and HF
!!IC0(1) contains information about the double excitation
!        LOGICAL :: WSign        !This is the sign of the particle
!        INTEGER , ALLOCATABLE :: History(:,:)         !This is the history of the particle. HF is not 
!included in the history - this is the start for all particles
!!The first double excitation from HF is stored in History(:,1)
!!The previous excitation (i.e. the one the particle was spawned from) is stored in History(:,ChainLength-1)
!!Is of size History(NEl,CLMax-1)
!        TYPE(ExcitGenerator) :: ExGen           !This is the excitation generator for the determinant 
!the particle is currently at
!
!!Might also be worth saving the excitation generators for all previous determinants
!    END TYPE
!    
!    INTEGER :: PartSize
!
!
!    TYPE(Part) , ALLOCATABLE , TARGET :: WalkVec(:),WalkVec2(:)      !These are the two arrays to hold the particle data
!    INTEGER :: WalkVecTag=0,WalkVec2Tag=0
!
!    TYPE(Part) , POINTER :: ActiveVec(:),NewVec(:)  !These are pointers to the two arrays
!
!    TYPE(ExcitGenerator) :: FDetExGen   !This is the excitation generator for the HF determinant
!
!    !MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
!    INTEGER :: MemoryFac=300
!
!    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,Iter
!    INTEGER :: exFlag=3
!
!    !This is information needed by the thermostating, so that the correct change in walker number can 
!be calculated, and hence the correct shift change.
!    !NoCulls is the number of culls in a given shift update cycle
!    INTEGER :: NoCulls=0
!    
!    !CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is 
!the previous number of steps before this cull...
!    !Only 15 culls/growth increases are allowed in a given shift cycle
!    INTEGER :: CullInfo(15,3)
!
!    real(dp) :: GrowRate
!    real(dp) :: ProjectionE,SumENum
!    integer(int64) :: SumNoatHF
!    INTEGER :: MinExit,MaxExit,SumWalkersCyc
!    real(dp) :: MeanExit
!
!    real(dp) :: Hii,rhii
!    HElement_t :: HFDiag

    contains

    SUBROUTINE ReturnPathMC(Weight,Energyxw)
        real(dp) :: Weight,Energyxw

        ! Avoid warnings
        weight = weight; energyxw = energyxw
!        INTEGER :: j
!        CHARACTER , PARAMETER :: this_routine='ReturnPathMC'
!        HElement_t :: HiiHEl,rhiiHEl

        CALL Stop_All("ReturnPathMC","This code has been commented out.")

!        OPEN(15,file='ReturnPathMCStats',status='unknown')
!
!        WRITE(6,*) ""
!        WRITE(6,*) "Performing particle spawning algorithm on returning paths..."
!        WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
!        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
!        WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
!
!        CALL InitRetPathMC()    !Initialise variables and allocate memory for calculation
!        
!        WRITE(6,*) ""
!        WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers        Proj.E   MeanExit   MinExit   MaxExit"
!        WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers         Proj.E   MeanExit   MinExit   MaxExit"
!
!        WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,G16.7,2I6)") 0,DiagSft,TotWalkers-TotWalkersOld,
                !GrowRate,TotWalkers,ProjectionE,MeanExit,0,MaxExit
!        WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,G16.7,2I6)") 0,DiagSft,TotWalkers-TotWalkersOld,
                !GrowRate,TotWalkers,ProjectionE,MeanExit,0,MaxExit
!
!!Start MC simulation
!        do Iter=1,NMCyc
!
!            CALL DoNMCyc()
!
!            IF(mod(Iter,StepsSft).eq.0) THEN
!
!                CALL UpdateSft()        !Update shift
!
!!Write out info
!                ProjectionE=SumENum/REAL(SumNoatHF,dp)
!                MeanExit=MeanExit/real(SumWalkersCyc,dp)
!                WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,G16.7,2I6)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,
                            !TotWalkers,ProjectionE,MeanExit,MinExit,MaxExit
!                WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,G16.7,2I6)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,
                                !TotWalkers,ProjectionE,MeanExit,MinExit,MaxExit
!
!                CALL neci_flush(15)
!                CALL neci_flush(6)       !Probably remove neci_flushes for big systems
!
!                TotWalkersOld=TotWalkers    !Reset 'old' number of walkers for next shift update cycle
!                MeanExit=0.0_dp
!                MinExit=NEl+10
!                MaxExit=0
!                SumWalkersCyc=0
!                
!
!            ENDIF
!
!        enddo   !End MC Cycle
!
!        Weight=(0.0_dp)
!        Energyxw=(SumENum/REAL(SumNoatHF,dp))
!
!!Deallocate Memory
!        do j=1,MaxWalkers
!            DEALLOCATE(WalkVec(j)%Det)
!            DEALLOCATE(WalkVec(j)%IC0)
!            DEALLOCATE(WalkVec(j)%Kii)
!            DEALLOCATE(WalkVec(j)%HistExcit)
!            DEALLOCATE(WalkVec(j)%History)
!            DEALLOCATE(WalkVec2(j)%Det)
!            DEALLOCATE(WalkVec2(j)%IC0)
!            DEALLOCATE(WalkVec2(j)%Kii)
!            DEALLOCATE(WalkVec2(j)%HistExcit)
!            DEALLOCATE(WalkVec2(j)%History)
!            IF(Allocated(WalkVec(j)%ExGen%ExcitData)) DEALLOCATE(WalkVec(j)%ExGen%ExcitData)
!            IF(Allocated(WalkVec2(j)%ExGen%ExcitData)) DEALLOCATE(WalkVec2(j)%ExGen%ExcitData)
!        enddo
!        DEALLOCATE(WalkVec)
!        DEALLOCATE(WalkVec2)
!        DEALLOCATE(HFDet)
!        CALL LogMemDealloc(this_routine,HFDetTag)
!
!        CLOSE(15)
!
!        RETURN

    END SUBROUTINE ReturnPathMC

!!This subroutine will perform the actual MC Cycles
!    SUBROUTINE DoNMCyc()
!        real(dp) :: Preturn,Ran2,Hij,hHi0,DiagElem,rat
!        INTEGER :: VecSlot,j,k,ToSpawn,IC,ExcitLevel,nJ(NEl),iDie
!        INTEGER :: iGetExcitLevel
!
!        VecSlot=1   !This indicates the next free slot in NewVec
!
!        do j=1,TotWalkers
!
!!First, sum in the energy contribution from the walker.
!            CALL SumInE(ActiveVec(j))
!
!!Now attempt to spawn from the determinant we are on - this can involve attempting to spawn back, with prob PRet, or to an excit
!            Preturn=FindPRet(ActiveVec(j))       !The PReturn may change depending on various parameters - calculate it here
!
!            IF((Preturn.eq.1.0_dp).or.(Preturn.gt.Ran2(Seed))) THEN
!!We are only allowed to generate the return determinant
!
!                ToSpawn=SpawnReturn(ActiveVec(j),Preturn)   !This tells us how many particles to spawn back, and their sign
!
!                do k=1,abs(ToSpawn)
!!We have decided that we want to spawn at the return determinant - copy through the correct number of new particles
!
!                    IF(ToSpawn.gt.0) THEN
!!New particles to copy through are positive
!                        CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength-1,.true.,VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!                    ELSE
!!New particles are negative
!                        CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength-1,.false.,VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!                    ENDIF
!
!                    VecSlot=VecSlot+1   !Increment the VecSlot to put particles into NewVec
!
!                enddo
!
!            ELSE
!!We are now going to attempt to spawn to all neighbours, rather than attempting to return to a previous determinant in the chain
!            
!                CALL SetupExitGen(ActiveVec(j)%Det,ActiveVec(j)%ExGen)      !Setup the excitation generator for this particle 
!if it isn't already there
!
!                do while(.true.)
!                    CALL GenSymExcitIt2(ActiveVec(j)%Det,NEl,G1,nBasis,nBasisMax,.FALSE.,  
                                !ActiveVec(j)%ExGen%ExcitData,nJ,IC,0,ActiveVec(j)%ExGen%nStore,exFlag)
!                    IF(nJ(1).eq.0) EXIT
!                    ExcitLevel=iGetExcitLevel(HFDet,nJ,NEl)
!                    IF((ActiveVec(j)%ChainLength.eq.0).or.(ExcitLevel.gt.ActiveVec(j)%IC0(ActiveVec(j)%ChainLength))) THEN
!!Excitation generated is one which is deeper into excitation space w.r.t. HF than the parent particle...see if we can spawn there.
!                        
!                        ToSpawn=SpawnForward(ActiveVec(j),PReturn,nJ,IC,Hij)    !Calculate number to spawn deeper into excit space
!
!                        IF(ToSpawn.ne.0) THEN
!!We are creating at least one new particle - find things out about it...
!                            DiagElem=GetConnection(nJ,nJ,0)   !Find diagonal element
!                            IF(ExcitLevel.eq.2) THEN
!!Particle will have connection to HF, therefore it must have been spawned from HFDet, and so hHi0=Hij
!                                hHi0=Hij
!                            ENDIF
!
!                            do k=1,abs(ToSpawn)
!!We have succesfully spawned deeper into excit space - copy these accross
!                                IF(ToSpawn.gt.0) THEN
!                                    CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength+1,.true.,
                                                !VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!                                ELSE
!                                    CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength+1,.false.,i
                                            !VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!                                ENDIF
!                    
!                                VecSlot=VecSlot+1
!
!                            enddo
!
!                        ENDIF   !End if we are spawning
!
!                    ENDIF   !End if excit is deeper into the space
!
!                enddo
! !Reset exgen so we can run through it again
!                CALL ResetExIt2(ActiveVec(j)%Det,NEl,G1,nBasis,nBasisMax,ActiveVec(j)%ExGen%ExcitData,0)   
!
!            ENDIF   !End of attempted spawning
!
!            iDie=AttemptDestruct(ActiveVec(j))  !Find how many particles are going to die...
!
!            IF(iDie.le.0) THEN
!!Particle survives (or possibly increases in number), and wants to be copied across to NewVec
!
!                do k=1,abs(iDie)+1
!                    CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength,ActiveVec(j)%WSign,VecSlot,
                                                !nJ,DiagElem,ExcitLevel,IC,hHi0)
!                    VecSlot=VecSlot+1
!                enddo
!            
!            ELSEIF(iDie.gt.1) THEN
!!Extra particles want to be killed
!                
!                do k=1,(iDie-1)
!                    IF(ActiveVec(j)%WSign) THEN
!                        CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength,.false.,VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!                        VecSlot=VecSlot+1
!                    ELSE
!                        CALL CreateParticle(ActiveVec(j),ActiveVec(j)%ChainLength,.true.,VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!                        VecSlot=VecSlot+1
!                    ENDIF
!                enddo
!            
!            ENDIF
!!If iDie=1, we simply do not copy the particle across and it is gone
!
!        enddo   !End cycling over walkers
!
!        SumWalkersCyc=SumWalkersCyc+TotWalkers
!        
!!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
!        TotWalkers=VecSlot-1
!        rat=(TotWalkers+0.0_dp)/(MaxWalkers+0.0_dp)
!        IF(rat.gt.0.9) THEN
!            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
!        ENDIF
!
!!Switch around the arrays that we are working on
!        IF(associated(ActiveVec,target=WalkVec)) THEN
!            ActiveVec=>WalkVec2
!            NewVec=>WalkVec
!        ELSE
!            ActiveVec=>WalkVec
!            NewVec=>WalkVec2
!        ENDIF
!
!        IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!!Particle number is too large - kill them randomly
!
!!Log the fact that we have made a cull
!            NoCulls=NoCulls+1
!            IF(NoCulls.gt.15) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!            CullInfo(NoCulls,1)=TotWalkers
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!
!            WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number..."
!            WRITE(6,*) "Killing randomly selected particles in order to reduce total number..."
!            WRITE(6,"(A,F8.2)") "Population will reduce by a factor of ",CullFactor
!            CALL ThermostatParticles(.true.)
!
!        ELSEIF(TotWalkers.lt.(InitWalkers/2)) THEN
!!Particle number is too small - double every particle in its current position
!
!!Log the fact that we have made a cull
!            NoCulls=NoCulls+1
!            IF(NoCulls.gt.15) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!            CullInfo(NoCulls,1)=TotWalkers
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!
!            WRITE(6,*) "Doubling particle population to increase total number..."
!            CALL ThermostatParticles(.false.)
!
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE DoNMCyc
!   
!!Sum in the energy contributions and excitation level statistics for a particle
!    SUBROUTINE SumInE(Particle)
!        TYPE(Part) :: Particle
!        INTEGER :: ExLevel
!
!        IF(Particle%ChainLength.eq.1) THEN
!!We are at a double excitation - sum the energy into SumENum
!            IF(Particle%WSign) THEN
!                IF(Iter.gt.NEquilSteps) SumENum=SumENum+Particle%Hi0
!            ELSE
!                IF(Iter.gt.NEquilSteps) SumENum=SumENum-Particle%Hi0
!            ENDIF
!            ExLevel=Particle%IC0(1)
!            MeanExit=MeanExit+ExLevel
!            IF(ExLevel.gt.MaxExit) MaxExit=ExLevel
!            IF(ExLevel.lt.MinExit) MinExit=ExLevel
!        ELSEIF(Particle%ChainLength.eq.0) THEN
!!We are at HF - sum in energy and correction to SumNoatHF
!!                IF(Particle%Hi0.ne.Hii) CALL Stop_All("DoNMCyc","Problem with particles at HF")
!            IF(Particle%WSign) THEN
!                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+1
!            ELSE
!                CALL Stop_All("DoNMCyc","Should not have negative particles at HF")
!                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF-1
!            ENDIF
!            ExLevel=0
!            IF(ExLevel.lt.MinExit) MinExit=ExLevel
!        ELSE
!!Sum into mean excitation level
!            ExLevel=Particle%IC0(Particle%ChainLength)
!            MeanExit=MeanExit+ExLevel
!            IF(ExLevel.gt.MaxExit) MaxExit=ExLevel
!            IF(ExLevel.lt.MinExit) MinExit=ExLevel
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE SumInE
!
!!This routine acts as a thermostat for the simulation - killing random particles if the population becomes too large, or 
!!Doubling them if it gets too low...
!    SUBROUTINE ThermostatParticles(HighLow)
!        IMPLICIT NONE
!        LOGICAL :: HighLow
!        INTEGER :: VecSlot,i,j,ToCull,Culled,OrigWalkers,Chosen
!        real(dp) :: Ran2
!
!        IF(HighLow) THEN
!!The population is too large - cull TotWalkers/CullFactor randomly selected particles
!
!            OrigWalkers=TotWalkers
!            ToCull=TotWalkers-nint((TotWalkers+0.0_dp)/CullFactor)
!            Culled=0
!
!            do while (Culled.lt.ToCull)
!
!!Pick a random walker between 1 and TotWalkers
!                Chosen=int((Ran2(Seed)*TotWalkers)+1.0_dp)
!
!!Move the particle at the end of the list to the position of the walker we have chosen to destroy
!                CALL CopyParticle(ActiveVec(TotWalkers),ActiveVec(Chosen))
!
!!Then deallocate the excitation generator for the particle at the end of the list
!                IF(Allocated(ActiveVec(TotWalkers)%ExGen%ExcitData)) DEALLOCATE(ActiveVec(TotWalkers)%ExGen%ExcitData)
!                ActiveVec(TotWalkers)%ExGen%ForCurrentDet=.false.
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
!                CALL CopyParticle(ActiveVec(i),ActiveVec(VecSlot))
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
!!This subroutine will create a particle on the NewVec array, given an initial template from ActiveVec, i
!which it may modify if it has moved up/down the chain
!!The routine will create a particle at the same place in a chain, the next link, or the previous history determinant.
!!If we want to create at a new determinant, further into excitation space than the Particle, then we have to 
!provide the extra information in the last 5 arguments
!    SUBROUTINE CreateParticle(Particle,NewChainLength,WSign,VecSlot,nJ,DiagElem,ExcitLevel,IC,hHi0)
!        TYPE(Part) :: Particle
!        INTEGER :: NewChainLength,VecSlot,nJ(NEl),ExcitLevel,IC
!        real(dp) :: DiagElem,hHi0
!        HElement_t :: ConntoHF
!        LOGICAL :: WSign
!
!        IF(Particle%ChainLength.eq.NewChainLength) THEN
!!We want to create an identical copy of the particle in NewVec (perhaps change the sign, but we shouldn't ever be doing this)
!            
!            IF(NewChainLength.eq.0) THEN
!                CALL SetupHFParticle(NewVec(VecSlot),WSign)     !Simply creating a particle at HF with no history
!            ELSE
!
!!Simply copy the data from Particle into an identical copy in NewVec
!                NewVec(VecSlot)%Kii(1:NewChainLength)=Particle%Kii(1:NewChainLength)
!                NewVec(VecSlot)%Hi0=Particle%Hi0
!                NewVec(VecSlot)%ChainLength=NewChainLength    !Though this should obviously be the same as the original particle
!                NewVec(VecSlot)%IC0(1:NewChainLength)=Particle%IC0(1:NewChainLength)
!                NewVec(VecSlot)%WSign=WSign
!                IF(NewChainLength.eq.1) THEN
!!Do not need to zero here, but may be useful for debugging purposes - only creating at a double excit, so no history
!                    NewVec(VecSlot)%History(:,:)=0.0_dp
!                ELSE
!                    NewVec(VecSlot)%History(:,1:(NewChainLength-1))=Particle%History(:,1:(NewChainLength-1))
!                ENDIF
!                NewVec(VecSlot)%HistExcit(1:NewChainLength)=Particle%HistExcit(1:NewChainLength)
!                NewVec(VecSlot)%Det(:)=Particle%Det(:)
!
!!Need to also copy accross the excitation generator for Det.
!                CALL CopyExGen(Particle%ExGen,NewVec(VecSlot)%ExGen) !This will copy the information exactly - IF it is correct
!
!            ENDIF
!
!        ELSEIF(Particle%ChainLength.eq.(NewChainLength+1)) THEN
!!We want to create a new particle which is one link closer to HF than Particle
!    
!            IF(NewChainLength.eq.0) THEN
!                CALL SetupHFParticle(NewVec(VecSlot),WSign)   !We are simply creating a particle at HF with no history
!            ELSE
!
!                NewVec(VecSlot)%Det=Particle%History(:,NewChainLength)  !Det is the last history
!
!                NewVec(VecSlot)%Kii(1:NewChainLength)=Particle%Kii(1:NewChainLength)
!                NewVec(VecSlot)%IC0(1:NewChainLength)=Particle%IC0(1:NewChainLength)
!                
!                IF(Particle%IC0(NewChainLength).gt.2) THEN
!!Excitation is more than a double, so connection to HF = 0
!                    NewVec(VecSlot)%Hi0=0.0_dp
!                ELSE
!!Need calculate the connection to HF - we are at a double excitation
!                    ConntoHF=GetHElement2(NewVec(VecSlot)%Det,FDet,NEl,nBasisMax,G1,nBasis,Brr,
                    !nMsh,fck,NMax,ALat,UMat,Particle%IC0(NewChainLength),ECore)
!                    NewVec(VecSlot)%Hi0=REAL(ConntoHF,dp)
!                ENDIF
!                
!                NewVec(VecSlot)%ChainLength=NewChainLength
!                NewVec(VecSlot)%WSign=WSign
!                IF(NewChainLength.eq.1) THEN
!!Again, do not have to zero this...
!                    NewVec(VecSlot)%History(:,:)=0.0_dp
!                ELSE
!                    NewVec(VecSlot)%History(:,1:(NewChainLength-1))=Particle%History(:,1:(NewChainLength-1))
!                ENDIF
!                NewVec(VecSlot)%HistExcit(1:NewChainLength)=Particle%HistExcit(1:NewChainLength)
!
!                IF(Allocated(NewVec(VecSlot)%ExGen%ExcitData)) DEALLOCATE(NewVec(VecSlot)%ExGen%ExcitData)
!                NewVec(VecSlot)%ExGen%ForCurrentDet=.false.
!
!            ENDIF
!
!
!        ELSEIF(Particle%ChainLength.eq.(NewChainLength-1)) THEN
!!We want to create a new particle which is one link further away from HF than Particle. We i
!need the additional information for this
!
!            NewVec(VecSlot)%Det=nJ(:)   !This is the new determinant
!
!            IF(NewChainLength.gt.1) THEN
!!Copy through the old stuff first, before we update with next link
!                NewVec(VecSlot)%Kii(1:(NewChainLength-1))=Particle%Kii(1:(NewChainLength-1))
!                NewVec(VecSlot)%IC0(1:(NewChainLength-1))=Particle%IC0(1:(NewChainLength-1))
!                NewVec(VecSlot)%HistExcit(1:(NewChainLength-1))=Particle%HistExcit(1:(NewChainLength-1))
!            ENDIF
!
!!Now update the new info
!            NewVec(VecSlot)%Kii(NewChainLength)=DiagElem        
!            NewVec(VecSlot)%IC0(NewChainLength)=ExcitLevel
!            NewVec(VecSlot)%ChainLength=NewChainLength
!            NewVec(VecSlot)%HistExcit(NewChainLength)=IC
!            NewVec(VecSlot)%Hi0=hHi0
!            NewVec(VecSlot)%WSign=WSign
!
!!Add det of previous particle to the history
!            IF(NewChainLength.gt.2) THEN
!                NewVec(VecSlot)%History(:,1:NewChainLength-2)=Particle%History(:,1:NewChainLength-2)
!                NewVec(VecSlot)%History(:,NewChainLength-1)=Particle%Det(:)
!            ELSEIF(NewChainLength.eq.2) THEN
!                NewVec(VecSlot)%History(:,NewChainLength-1)=Particle%Det(:)
!            ELSEIF(NewChainLength.eq.1) THEN
!!No history, since we are at double excit
!                IF(ExcitLevel.ne.2) CALL Stop_All("CreateParticle","First link in chain is not a double excit")
!            ELSE
!                CALL Stop_All("CreateParticle","Wanting to create more distant link in chain, but is trying to create HF")
!            ENDIF
!
!!The excitgen has not been created
!            IF(Allocated(NewVec(VecSlot)%ExGen%ExcitData)) DEALLOCATE(NewVec(VecSlot)%ExGen%ExcitData)
!            NewVec(VecSlot)%ExGen%ForCurrentDet=.false.
!
!        ELSE
!!Not sure what we are trying to do here - we are creating a particle which is > 1 link away 
!from the original - should not be possible
!            CALL Stop_All("CreateParticle","Trying to create a particle which is > 1 link away - error")
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE CreateParticle
!
!
!!This routine sets up a HF particle - assuming that the FDetExGen has been proparly initialised
!    SUBROUTINE SetupHFParticle(Particle,WSign)
!        TYPE(Part) :: Particle
!        LOGICAL :: WSign
!        INTEGER :: ierr,j,k
!        
!!Since we are at HF, there is no Kii or IC0 contribution since they start at the first link in the chain
!        Particle%Det(:)=HFDet(:)     !All walkers start at HF
!        Particle%Hi0=Hii            !Connection of HF is simply HF energy
!        Particle%ChainLength=0      !Initial Chain Length = 0 as at HF
!        Particle%WSign=WSign
!!Can remove the subsequent zeroing if we are sure that the code is bug-free. Otherwise, it may help for debugging        
!        Particle%History(:,:)=0.0_dp
!        Particle%HistExcit(:)=0
!        Particle%IC0(:)=0
!        Particle%Kii(:)=0.0_dp
!        CALL CopyExGen(FDetExGen,Particle%ExGen)
!
!    END SUBROUTINE SetupHFParticle
!
!!Copy a particle
!    SUBROUTINE CopyParticle(OrigPart,NewPart)
!        TYPE(Part) :: OrigPart,NewPart
!
!        NewPart%Det(:)=OrigPart%Det(:)
!        NewPart%ChainLength=OrigPart%ChainLength
!        NewPart%IC0(:)=OrigPart%IC0(:)
!        NewPart%Kii(:)=OrigPart%Kii(:)
!        NewPart%HistExcit(:)=OrigPart%HistExcit(:)
!        NewPart%Hi0=OrigPart%Hi0
!        NewPart%WSign=OrigPart%WSign
!        NewPart%History(:,:)=OrigPart%History(:,:)
!        IF(Allocated(NewPart%ExGen%ExcitData)) DEALLOCATE(NewPart%ExGen%ExcitData)
!        NewPart%ExGen%ForCurrentDet=.false.
!        CALL CopyExGen(OrigPart%ExGen,NewPart%ExGen)    !Attempt to copy across exGen
!
!        RETURN
!
!    END SUBROUTINE CopyParticle
!
!!This routine will copy an exgen from OrigExGen to NewExGen, providing that OrigExGen is said i
!to be allocated for the 'correct' determinant
!    SUBROUTINE CopyExGen(OrigExGen,NewExGen)
!        TYPE(ExcitGenerator) :: OrigExGen,NewExGen
!        INTEGER :: ierr,i
!
!        IF(.not.OrigExGen%ForCurrentDet) THEN
!!The particle we are trying to create does not have the correct corresponding excitation generator
!            NewExGen%ForCurrentDet=.false.
!!            CALL Stop_All("CopyExGen","Trying to copy accross an exGen, but it is not for the desired determinant")
!            RETURN
!        ENDIF
!        IF(Allocated(NewExGen%ExcitData)) THEN
!            DEALLOCATE(NewExGen%ExcitData)
!        ENDIF
!        ALLOCATE(NewExGen%ExcitData(OrigExGen%ExcitLen),stat=ierr)
!        IF(ierr.ne.0) CALL Stop_All("CopyExGen","Problem in allocation of ExGen")
!        do i=1,OrigExGen%ExcitLen
!            NewExGen%ExcitData(i)=OrigExGen%ExcitData(i)
!        enddo
!        NewExGen%nStore(:)=OrigExGen%nStore(:)
!        NewExGen%ExcitLen=OrigExGen%ExcitLen
!        NewExGen%ForCurrentDet=.true.
!        RETURN
!    END SUBROUTINE CopyExGen
!
!!Function which returns the number of particles to kill (or birth if -ve - possible if shift is positive)    
!!Zero indicates keep particle. One is destroy particle. -1 is create extra particle. >1 is destroy i
!particle, and create anti-particles
!    INTEGER FUNCTION AttemptDestruct(Particle)
!        TYPE(Part) :: Particle
!        real(dp) :: rat,Ran2
!
!        IF(TRhoElems) THEN
!            IF(Particle%ChainLength.eq.0) THEN
!!We are at HF, so Kii is zero
!                rat=EXP(Tau*DiagSft)
!            ELSE
!                rat=GetSpawnRhoEl(Particle%Det,Particle%Det,.true.,Particle%Kii(Particle%ChainLength))
!            ENDIF
!!GetSpawnRhoEl will recover the spawning rate. The death rate will be 1 - this.
!            rat=1.0_dp-rat
!
!        ELSE
!            IF(Particle%ChainLength.eq.0) THEN
!!We are at HF, so Kii is zero
!                rat=-Tau*DiagSft
!            ELSE
!                rat=Tau*((Particle%Kii(Particle%ChainLength))-DiagSft)  !Prob of death
!            ENDIF
!
!        ENDIF
!
!        AttemptDestruct=INT(rat)
!        rat=rat-REAL(AttemptDestruct,dp)
!
!        IF(abs(rat).gt.Ran2(Seed)) THEN
!            IF(rat.ge.0.0_dp) THEN
!                AttemptDestruct=AttemptDestruct+1
!            ELSE
!                AttemptDestruct=AttemptDestruct-1
!            ENDIF
!        ENDIF
!
!        RETURN
!
!    END FUNCTION AttemptDestruct
!
!!This is redone since we want to change the diagonal elements due to the shift + we want tau to be as defined in the module
!    real(dp) FUNCTION GetSpawnRhoEl(nI,nJ,LSame,Conn)
!        INTEGER :: nI(NEl),nJ(NEl),iGetExcitLevel
!        LOGICAL :: LSame
!        real(dp) :: Conn      !For off-diagonal, this = Hij. For Diagonal, this equals Ei-E0
!        real(dp) :: UExp
!        HElement_t :: EDiag,RH,EDiag2
!
!        IF(LSame) THEN
!
!!Calculate a diagonal rho matrix element
!!            IF(NTAY(2).eq.3) THEN
!!Partition with Trotter with H(0) having just the Fock Operators
!!Fock-Partition-Lowdiag
!                RH=EXP(-Tau*(Conn-DiagSft))
!!            ELSEIF(NTAY(2).eq.2) THEN
!!Partition with Trotter with H(0) having just the Fock Operators
!!Fock-Partition
!!                CALL Stop_All("GetSpawnRhoElement","This is not functional yet")
!!            ENDIF
!
!        ELSE
!
!!We want a diagonal element - this is exactly the same as in calcrho2
!            call GetH0Element(nI,nEl,Arr,nBasis,ECore,EDiag2)
!            call GetH0Element(nJ,nEl,Arr,nBasis,ECore,EDiag)
!            EDiag=(EDiag2+EDiag)/(2.0_dp)
!            UExp=-Tau*Conn
!            RH=EXP(-Tau*REAL(EDiag,dp))*UExp
!
!        ENDIF
!
!        GetSpawnRhoEl=REAL(RH,dp)
!
!    END FUNCTION GetSpawnRhoEl
!
!
!!Call this function when we have decided to spawn to a determinant further down the chain and have 
!decided that it is deeper into excitation space w.r.t HF
!!Will return the number and sign of the particles spawned there (nJ). IC is the no. of 
!excitations between particle and nJ.
!    INTEGER FUNCTION SpawnForward(Particle,Preturn,nJ,IC,Hij)
!        TYPE(Part) :: Particle
!        real(dp) :: Preturn,Hij,rat,Ran2,rhoel
!        INTEGER :: nJ(:),IC
!
!        IF(Preturn.eq.1.0_dp) CALL Stop_All("SpawnForward","Preturn=1, but trying to spawn forward")
!
!!First calculate connection to parent determinant - this is the hamiltonian matrix element
!        Hij=GetConnection(Particle%Det,nJ,IC)
!
!        IF(TRhoElems) THEN
!            IF(IC.eq.0) CALL Stop_All("SpawnForward","IC should not be zero")
!            rhoel=GetSpawnRhoEl(Particle%Det,nJ,.false.,Hij)
!            rat=abs(rhoel)/(1.0_dp-Preturn)
!        ELSE
!            rat=Tau*abs(Hij)/(1.0_dp-Preturn)
!        ENDIF
!        SpawnForward=INT(rat)
!        rat=rat-REAL(SpawnForward,dp)
!        IF(rat.gt.Ran2(Seed)) SpawnForward=SpawnForward+1   !Stochastic step successful at creating another particle
!        
!        IF(SpawnForward.gt.0) THEN
!!Attempt to spawn at return determinant is successful - determine sign of new particles
!            IF(TRhoElems) THEN
!                IF(.not.Particle%WSign) SpawnForward=-SpawnForward
!                IF(rhoel.lt.0.0_dp) SpawnForward=-SpawnForward
!            ELSE
!                IF((Particle%WSign).and.(Hij.gt.0.0_dp)) THEN             !Positive particle & connection
!                    SpawnForward=-SpawnForward                          !New particles negative
!                ELSEIF((.not.Particle%WSign).and.(Hij.lt.0.0_dp)) THEN    !Negative particle & connection
!                    SpawnForward=-SpawnForward                          !New particles negative
!                ENDIF
!            ENDIF
!        ENDIF
!
!        RETURN
!
!    END FUNCTION SpawnForward
!
!                
!!Call this function after we have decided to attempt to spawn at the determinant that we ourselves were spawned from.
!!SpawnReturn will then tell us if the attempt is successful by indicating the number of particles 
!spawned and sign of the resultant particle(s)
!    INTEGER FUNCTION SpawnReturn(Particle,Preturn)
!        TYPE(Part) :: Particle
!        real(dp) :: Hij,Ran2,Preturn,rat,rhoel
!        INTEGER :: nJ(NEl)
!
!!First, need to find return determinant, and calculate connection to it
!        IF(Particle%ChainLength.eq.1) THEN
!!Particle has a chainlength of 1, (i.e. is at a double excitation) - therefore it wants to try to return to HF
!            IF((Particle%HistExcit(1).ne.2).or.(Particle%IC0(1).ne.2)) THEN
!!Quick check that we are actually at a double excitation
!                CALL Stop_All("SpawnReturn","Chainlength is one, but we are not at a double excitation")
!            ENDIF
!
!            Hij=GetConnection(HFDet,Particle%Det,2)
!            IF(TRhoElems) THEN
!                nJ(:)=HFDet(:)
!            ENDIF
!
!        ELSEIF(Particle%ChainLength.gt.1) THEN
!!Particle is some way along a chain. Its previous determinant can be ascertained from its history.
!!The previous determinant is stored in History(:,ChainLength-1), and the number of excitations 
!between them in HistExcit(ChainLength)
!
!            Hij=GetConnection(Particle%Det,Particle%History(:,(Particle%ChainLength)-1),Particle%HistExcit(Particle%ChainLength))
!            IF(TRhoElems) THEN
!                nJ(:)=Particle%History(:,(Particle%ChainLength)-1)
!            ENDIF
!
!        ELSE
!!Particle has a chainlength of 0 (or less than zero!) - it cannot spawn to a previous determinant, as there are none - error here
!            CALL Stop_All("SpawnReturn","Cannot return if at HF")
!        ENDIF
!
!!Prob of accepting to spawn to a previous determinant given by tau*abs(Hij)/Preturn (rho matrix elements already has tau included)
!        IF(TRhoElems) THEN
!            rhoel=GetSpawnRhoEl(Particle%Det,nJ,.false.,Hij)
!            rat=abs(rhoel)/Preturn
!        ELSE
!            rat=Tau*abs(Hij)/Preturn
!        ENDIF
!        SpawnReturn=INT(rat)              !Number of particles we are definatly creating 
!(can create multiple new particles if prob > 1)
!        rat=rat-REAL(SpawnReturn,dp)
!
!        IF(rat.gt.Ran2(Seed)) SpawnReturn=SpawnReturn+1     !Create an extra particle from stochastic step
!
!        IF(SpawnReturn.gt.0) THEN
!!Attempt to spawn at return determinant is successful - determine sign of new particles
!            IF(TRhoElems) THEN
!                IF(.not.Particle%WSign) SpawnReturn=-SpawnReturn
!                IF(rhoel.lt.0.0_dp) SpawnReturn=-SpawnReturn
!            ELSE
!                IF((Particle%WSign).and.(Hij.gt.0.0_dp)) THEN             !Positive particle & connection
!                    SpawnReturn=-SpawnReturn                            !New particles negative
!                ELSEIF((.not.Particle%WSign).and.(Hij.lt.0.0_dp)) THEN    !Negative particle & connection
!                    SpawnReturn=-SpawnReturn                            !New particles negative
!                ENDIF
!            ENDIF
!        ENDIF
!
!        RETURN
!
!    END FUNCTION SpawnReturn
!
!!This function simply gets the connection between two determinants - i.e. the Hij element
!!For diagonal elements (IC.eq.0) then the Hii element is found and the EHF subtracted from it, to give the Kii element
!!However, if TRhoElems is on, then it returns Ei-E0 for diagonal elements.
!    real(dp) FUNCTION GetConnection(nI,nJ,IC)
!        INTEGER :: nI(NEl),nJ(NEl),IC
!        HElement_t :: rhiiHEl,HiiHEl
!
!        IF(TRhoElems.and.(IC.eq.0)) THEN
!!We want to calculate rho transition matrix elements
!
!!            CALL CalcRho2(nI,nJ,Beta,i_P,NEl,G1,nBasis,nMsh,fck,Arr,ALat,UMat,rhiiHEl,nTay,IC,ECore)
!!            GetConnection=REAL(rhiiHEl,dp)
!
!!            IF(IC.eq.0) THEN
!!We are after Hii-H00 so subtract the reference energy
!!                GetConnection=GetConnection-rhii
!!            ENDIF
!
!            call GetH0Element(nI,NEl,Arr,nBasis,ECore,rhiiHEl)
!            rhiiHEl=rhiiHEl-HFDiag
!            GetConnection=Real(rhiiHEl,dp)
!
!        ELSE
!!We want to calculate Hamiltonian transition matrix elements
!                
!            HiiHEl=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,IC,ECore)
!            GetConnection=REAL(HiiHEl,dp)
!            
!            IF(IC.eq.0) THEN
!!We are after Hii-H00 so subtract the reference energy
!                GetConnection=GetConnection-Hii
!            ENDIF
!
!        ENDIF
!
!        RETURN
!
!    END FUNCTION GetConnection
!
!!This function finds the probability of returning the way it came
!    real(dp) FUNCTION FindPRet(Particle)
!        TYPE(Part) :: Particle
!
!        IF(Particle%ChainLength.eq.0) THEN
!!We are at the HF - there is no possibility to 'return' - set Pret to zero
!            FindPRet=0.0_dp
!            RETURN
!        ELSEIF(Particle%ChainLength.eq.CLMax) THEN
!!Length of chain is the maximum allowed chain length - force the particle to attempt to return - set PRet to one
!            FindPRet=1.0_dp
!            RETURN
!        ELSEIF(Particle%IC0(Particle%ChainLength).eq.NEl) THEN
!!Particle is already at the highest excitation level possible - it has to return as it can no longer increase its 
!excitation level
!            FindPRet=1.0_dp
!            RETURN
!        ELSEIF(Particle%IC0(Particle%ChainLength).eq.0) THEN
!!Particle says it is at the HF - however, it should only be allowed to go deeper into excit space, and so 
!chainlength should be 0 - error here
!            CALL Stop_All("FindPRet","Particle should only be allowed to go deeper into excit space")
!        ELSE
!!Particle is free to attempt to return, or proceed further into excit space
!            FindPRet=PRet
!            RETURN
!        ENDIF
!        RETURN
!
!    END FUNCTION FindPRet
!
!!This routine is used to initialise the calculation for the return pathMC - allocate memory and put all initial particles on HF
!    SUBROUTINE InitRetPathMC()
!        INTEGER :: ierr,j
!        CHARACTER , PARAMETER :: this_routine='InitRetPath'
!        HElement_t :: HiiHel,rhiiHel
!
!        SumWalkersCyc=0
!        MeanExit=0.0_dp
!        MaxExit=0
!        MinExit=NEl+10
!
!        ALLOCATE(HFDet(NEl),stat=ierr)
!        CALL LogMemAlloc('HFDet',NEl,4,this_routine,HFDetTag)
!        do j=1,NEl
!            HFDet(j)=FDet(j)
!        enddo
!
!!Provide various tests that variables are within allowed ranges
!        IF((PRet.gt.1.0_dp).or.(PRet.lt.0.0_dp)) CALL Stop_All("ReturnPathMC","PRet must be a normalised probability")
!        IF(HElement_t_size.gt.1) CALL Stop_All("ReturnPathMC","ReturnPathMC cannot function with complex orbitals.")
!
!        MaxWalkers=InitWalkers*MemoryFac    !Set maximum number of allowed walkers
!
!!Calculate Hii and rhii
!        HiiHEl=GetHElement2(HFDet,HFDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
!        Hii=REAL(HiiHEl,dp)
!        CALL CalcRho2(HFDet,HFDet,Beta,i_P,NEl,G1,nBasis,nMsh,fck,Arr,ALat,UMat,rhiiHEl,nTay,0,ECore)
!        rhii=REAL(rhiiHEl,dp)
!        call GetH0Element(HFDet,NEl,Arr,nBasis,ECore,HFDiag)
!
!        IF(TRhoElems) THEN
!            WRITE(6,"(A,F20.10)") "Root diagonal element is: ",rhii
!        ELSE
!            WRITE(6,"(A,F20.10)") "Root diagonal element is: ",Hii
!        ENDIF
!
!        PartSize=((CLMax*(NEl+4))+4)*4      !This is the size of a particle (without exgen stuff)
!
!        WRITE(6,"(A)",advance='no') "Attempting to allocate initial memory....."
!        CALL neci_flush(6)
!
!!Allocate memory to hold all walkers
!        ALLOCATE(WalkVec(MaxWalkers),stat=ierr)
!        CALL LogMemAlloc('WalkVec',MaxWalkers,PartSize,this_routine,WalkVecTag,ierr)
!        ALLOCATE(WalkVec2(MaxWalkers),stat=ierr)
!        CALL LogMemAlloc('WalkVec2',MaxWalkers,PartSize,this_routine,WalkVec2Tag,ierr)
!
!        do j=1,MaxWalkers
!!Run through all the arrays, and allocate memory for the walkers in the various arrays
!            ALLOCATE(WalkVec(j)%Det(NEl),stat=ierr)
!            ALLOCATE(WalkVec(j)%IC0(CLMax),stat=ierr)
!            ALLOCATE(WalkVec(j)%Kii(CLMax),stat=ierr)
!            ALLOCATE(WalkVec(j)%HistExcit(CLMax),stat=ierr)
!            ALLOCATE(WalkVec(j)%History(NEl,CLMax-1),stat=ierr)
!            ALLOCATE(WalkVec2(j)%Det(NEl),stat=ierr)
!            ALLOCATE(WalkVec2(j)%IC0(CLMax),stat=ierr)
!            ALLOCATE(WalkVec2(j)%Kii(CLMax),stat=ierr)
!            ALLOCATE(WalkVec2(j)%HistExcit(CLMax),stat=ierr)
!            ALLOCATE(WalkVec2(j)%History(NEl,CLMax-1),stat=ierr)
!            IF(ierr.ne.0) CALL Stop_All("InitRetPath","Problem allocating initial memory")
!        enddo
!
!        ActiveVec=>WalkVec      !Point ActiveVec to WalkVec initially
!        NewVec=>WalkVec2
!
!        WRITE(6,*) "DONE"
!        CALL neci_flush(6)
!
!        FDetExGen%ForCurrentDet=.false.
!        CALL SetupExitGen(HFDet,FDetExGen)
!
!        do j=1,InitWalkers
!            CALL SetupHFParticle(ActiveVec(j),.true.)
!        enddo
!
!        TotWalkers=InitWalkers      !Setup the variable for the number of walkers
!        TotWalkersOld=InitWalkers
!        
!!Initialise other variables needed
!        Seed=G_VMC_Seed     !Initialise random number seed
!        GrowRate=0.0_dp
!        CullInfo=0
!        NoCulls=0
!!Initialise variables for calculation of the running average
!        ProjectionE=0.0_dp
!        SumENum=0.0_dp
!        SumNoatHF=0
!
!        RETURN
!    END SUBROUTINE InitRetPathMC
!
!!This routine sets up an excitation generator for nI in ExcitGen, of length nExcitMemLen
!!ExcitGen%ForCurrentDet is a logical which tells us whether the excitation generator is already created for the determinant
!    SUBROUTINE SetupExitGen(nI,ExcitGen)
!        INTEGER :: ierr,nI(NEl),iMaxExcit,nJ(NEl)
!        TYPE(ExcitGenerator) :: ExcitGen
!
!        IF(ExcitGen%ForCurrentDet) THEN
!
!            IF(.not.Allocated(ExcitGen%ExcitData)) THEN
!                CALL Stop_All("SetupExitGen","ExGen supposedly already created, but not present")
!            ENDIF
!
!            RETURN
!
!        ELSE
!!Need to create excitation generator
!            IF(Allocated(ExcitGen%ExcitData)) THEN
!                DEALLOCATE(ExcitGen%ExcitData)      !Deallocate old data if it is there
!            ENDIF
!
!!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
!            iMaxExcit=0
!            ExcitGen%nStore(1:6)=0
!            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGen%ExcitLen,nJ,iMaxExcit,0,ExcitGen%nStore,exFlag)
!            ALLOCATE(ExcitGen%ExcitData(ExcitGen%ExcitLen),stat=ierr)
!            IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
!            ExcitGen%ExcitData(1)=0
!            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGen%ExcitData,nJ,iMaxExcit,0,ExcitGen%nStore,exFlag)
!            ExcitGen%ForCurrentDet=.true.
!
!        ENDIF
!
!        RETURN
!    END SUBROUTINE SetupExitGen
!
!!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!!value of the diagonal shift in the hamiltonian in order to compensate for this
!    SUBROUTINE UpdateSft()
!        IMPLICIT NONE
!        INTEGER :: j,k,GrowthSteps
!
!        IF(NoCulls.eq.0) THEN
!            GrowRate=(TotWalkers+0.0_dp)/(TotWalkersOld+0.0_dp)
!        ELSEIF(NoCulls.eq.1) THEN
!!GrowRate is the sum of the individual grow rates for each uninterrupted growth sequence, multiplied by the 
!fraction of the cycle which was spent on it
!            GrowRate=((CullInfo(1,3)+0.0_dp)/(StepsSft+0.0_dp))*((CullInfo(1,1)+0.0_dp)/(TotWalkersOld+0.0_dp))
!            GrowRate=GrowRate+(((StepsSft-CullInfo(1,3))+0.0_dp)/(StepsSft+0.0_dp))*((TotWalkers+0.0_dp)/(CullInfo(1,2)+0.0_dp))
!
!            NoCulls=0
!            CullInfo=0
!        ELSE
!            GrowRate=((CullInfo(1,3)+0.0_dp)/(StepsSft+0.0_dp))*((CullInfo(1,1)+0.0_dp)/(TotWalkersOld+0.0_dp))
!            do j=2,NoCulls
!
!!This is needed since the steps between culls are stored cumulatively
!                GrowthSteps=CullInfo(j,3)-CullInfo(j-1,3)
!                GrowRate=GrowRate+((GrowthSteps+0.0_dp)/(StepsSft+0.0_dp))*((CullInfo(j,1)+0.0_dp)/(CullInfo(j-1,2)+0.0_dp))
!
!            enddo
!
!            GrowthSteps=StepsSft-CullInfo(NoCulls,3)
!            GrowRate=GrowRate+((GrowthSteps+0.0_dp)/(StepsSft+0.0_dp))*((TotWalkers+0.0_dp)/(CullInfo(NoCulls,2)+0.0_dp))
!
!            NoCulls=0
!            CullInfo=0
!
!        ENDIF
!        DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.0_dp))
!!        IF((DiagSft).gt.0.0_dp) THEN
!!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!!            STOP
!!        ENDIF
!
!    END SUBROUTINE UpdateSft
!
END MODULE ReturnPathMCMod
