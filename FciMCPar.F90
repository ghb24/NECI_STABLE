#ifdef PARALLEL
!This is a parallel MPI version of the FciMC code. It picks random excitations to work with, and can fully diagonalise the 
!resultant 2v graph, or apply it many times.
!Excitation generators are now stored along with the particles (but in a separate array)
!All variables refer to values per processor

MODULE FciMCParMod
    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr
    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    USE Calc , only : TStartMP1
    USE Calc , only : GrowMaxFactor,CullFactor
    USE Calc , only : RhoApp,TResumFCIMC
    USE Determinants , only : FDet,GetHElement2
    USE DetCalc , only : NMRKS
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    USE Parallel
    IMPLICIT NONE
    SAVE

    INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
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
    INTEGER :: MemoryFac=3000

    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps
    INTEGER :: exFlag=3

!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle for each variable
    INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
    INTEGER :: CullInfo(10,3)

!The following variables are calculated as per processor, but at the end of each update cycle, are combined to the root processor
    REAL*8 :: GrowRate,DieRat,ProjectionE,SumENum
    INTEGER*8 :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
    REAL*8 :: PosFrac           !This is the fraction of positive particles on each node
    INTEGER :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
    REAL*8 :: MeanExcitLevel    
    INTEGER :: MinExcitLevel
    INTEGER :: MaxExcitLevel

!These are the global variables, calculated on the root processor, from the values above
    REAL*8 :: AllGrowRate,AllMeanExcitLevel
    INTEGER :: AllMinExcitLevel,AllMaxExcitLevel,AllTotWalkers,AllTotWalkersOld,AllSumWalkersCyc
    REAL*8 :: AllSumNoatHF,AllSumENum,AllPosFrac

    REAL*8 :: MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector

    TYPE(HElement) :: rhii,FZero
    REAL*8 :: Hii

    REAL*8 , ALLOCATABLE :: GraphRhoMat(:,:)    !This stores the rho matrix for the graphs in resumFCIMC
    INTEGER :: GraphRhoMatTag=0

    REAL*8 , ALLOCATABLE :: GraphVec(:)         !This stores the final components for the propagated graph in ResumFCIMC
    INTEGER :: GraphVecTag=0

    INTEGER , ALLOCATABLE :: DetsinGraph(:,:)   !This stores the determinants in the graph created for ResumFCIMC
    INTEGER :: DetsinGraphTag=0

    REAL*8 :: RootExcitProb     !This is the probability of generating an excitation from the current root in the ResumFCIMC current graph.

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j
        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'
        TYPE(HElement) :: Hamii

        CALL InitFCIMCCalcPar()

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(6,*) ""
            WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers        Proj.E      Net+veWalk     MeanExcitLevel   MinExcitLevel   MaxExcitLevel"
            WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers         Proj.E      Net+veWalk     MeanExcitLevel   MinExcitLevel   MaxExcitLevel"
!TotWalkersOld is the number of walkers last time the shift was changed
            WRITE(6,"(I9,G16.7,I9,G16.7,I9,2G16.7,G16.7,F16.7,2I5)") 0,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,ProjectionE,AllPosFrac,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
            WRITE(15,"(I9,G16.7,I9,G16.7,I9,2G16.7,G16.7,F16.7,2I5)") 0,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,ProjectionE,AllPosFrac,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
        ENDIF
        

!Start MC simulation...
        do Iter=1,NMCyc
            
            CALL PerformFCIMCyc()

            IF(mod(Iter,StepsSft).eq.0) THEN
!This will communicate between all nodes, find the new shift (and other parameters) and broadcast them to the other nodes.
                CALL CalcNewShift()
            ENDIF

!End of MC cycle
        enddo

        Weight=HDElement(0.D0)
        Energyxw=HDElement(DiagSft)

!Deallocate memory
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

        IF(iProcIndex.eq.Root) CLOSE(15)

        CALL MPIEnd(.false.)    !Finalize MPI

        RETURN

    END SUBROUTINE FciMCPar

!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCycPar()
        INTEGER :: VecSlot,i,j,k,l
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount
        REAL*8 :: Prob,rat
        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,iGetExcitLevel
        LOGICAL :: WSign
        TYPE(HElement) :: HDiag,HOffDiag
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1

        do j=1,TotWalkers
!j runs through all current walkers

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(CurrentDets(:,j),CurrentIC(j),CurrentH(2,j),CurrentSign(j))

!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
            CALL SetupExitgenPar(CurrentDets(:,j),CurrentExcits(j))

            CALL GenRandSymExcitIt3(CurrentDets(:,j),CurrentExcits(j)%ExcitData,nJ,Seed,IC,0,Prob,iCount)

            IF(TResumFCIMC) THEN
                WRITE(6,*) "Resum facility not yet operational"
                CALL FLUSH(6)
                CALL MPIStopAll(1)

                CALL ResumFciMCPar(CurrentDets(:,j),CurrentSign(j),CurrentH(1,j),nJ,IC,Prob,VecSlot,CurrentExcits(j),j)
            ELSE
!Calculate number of children to spawn
 
                Child=AttemptCreatePar(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC)
                
                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store

                    IF(Child.gt.0) THEN
!We have successfully created at least one positive child at nJ
                        WSign=.true.
                    ELSE
!We have successfully created at least one negative child at nJ
                        WSign=.false.
                    ENDIF
!Calculate excitation level, connection to HF and diagonal ham element
                    ExcitLevel=iGetExcitLevel(HFDet,nJ,NEl)
                    IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
                        HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                    ENDIF
                    IF(ExcitLevel.eq.0) THEN
!We know we are at HF - HDiag=0
                        HDiag=HElement(0.D0)
                    ELSE
                        HDiag=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                    ENDIF

            
                    do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them
                        NewDets(:,VecSlot)=nJ(:)
                        NewSign(VecSlot)=WSign
                        NewExcits(VecSlot)%ExitGenForDet=.false.
                        NewIC(VecSlot)=ExcitLevel
                        NewH(1,VecSlot)=REAL(HDiag%v,r2)-Hii      !Diagonal H-element-Hii
                        NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                        VecSlot=VecSlot+1
                    enddo
                
                ENDIF   !End if child created

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
                iDie=AttemptDiePar(CurrentDets(:,j),CurrentH(1,j))
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births

                IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle
        
                    do l=1,abs(iDie)+1    !We need to copy accross one more, since we need to include the original spared particle
                        NewDets(:,VecSlot)=CurrentDets(:,j)
                        NewSign(VecSlot)=CurrentSign(j)
!Copy excitation generator accross
                        CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot))
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
                        CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot))
                        NewIC(VecSlot)=CurrentIC(j)
                        NewH(:,VecSlot)=CurrentH(:,j)
                        VecSlot=VecSlot+1
                    enddo
            
                ENDIF   !To kill if

            ENDIF   !Resum if

!Finish cycling over walkers
        enddo
        
!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+TotWalkers

!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
        TotWalkers=VecSlot-1
        rat=(TotWalkers+0.D0)/(MaxWalkers+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
        ENDIF
        
!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
        IF(associated(CurrentDets,target=WalkVecDets)) THEN
            CurrentDets=>WalkVec2Dets
            CurrentSign=>WalkVec2Sign
            CurrentExcits=>WalkVec2Excits
            NewDets=>WalkVecDets
            NewSign=>WalkVecSign
            NewExcits=>WalkVecExcits
        ELSE
            CurrentDets=>WalkVecDets
            CurrentSign=>WalkVecSign
            CurrentExcits=>WalkVecExcits
            NewDets=>WalkVec2Dets
            NewSign=>WalkVec2Sign
            NewExcits=>WalkVec2Excits
        ENDIF

        IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!Particle number is too large - kill them randomly

!Log the fact that we have made a cull
            NoCulls=NoCulls+1
            IF(NoCulls.gt.10) THEN
                WRITE(6,*) "Too Many Culls"
                CALL FLUSH(6)
                CALL MPIStopAll(2)
            ENDIF
!CullInfo(:,1) is walkers before cull
            CullInfo(NoCulls,1)=TotWalkers
!CullInfo(:,3) is MC Steps into shift cycle before cull
            CullInfo(NoCulls,3)=mod(Iter,StepsSft)

!            WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number..."
!            WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number..."
!            WRITE(6,"(A,F8.2)") "Population will reduce by a factor of ",CullFactor
            CALL ThermostatParticlesPar(.true.)

        ELSEIF(TotWalkers.lt.(InitWalkers/2)) THEN
!Particle number is too small - double every particle in its current position

!Log the fact that we have made a cull
            NoCulls=NoCulls+1
            IF(NoCulls.gt.10) CALL STOPGM("PerformFCIMCyc","Too Many Culls")
!CullInfo(:,1) is walkers before cull
            CullInfo(NoCulls,1)=TotWalkers
!CullInfo(:,3) is MC Steps into shift cycle before cull
            CullInfo(NoCulls,3)=mod(Iter,StepsSft)
            
!            WRITE(6,*) "Doubling particle population to increase total number..."
            CALL ThermostatParticlesPar(.false.)

        ENDIF

        RETURN

    END SUBROUTINE PerformFCIMCycPar

!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,Hij0,WSign)
        INTEGER :: DetCurr(NEl),ExcitLevel
        LOGICAL :: WSign
        REAL*8 :: Hij0      !This is the hamiltonian matrix element between DetCurr and HF

        MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,r2)
        IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
        IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
        IF(ExcitLevel.eq.0) THEN
            IF(WSign) THEN
                SumNoatHF=SumNoatHF+1
                PosFrac=PosFrac+1.D0
            ELSE
                SumNoatHF=SumNoatHF-1
            ENDIF
        ELSEIF(ExcitLevel.eq.2) THEN
!At double excit - sum in energy
            IF(WSign) THEN
                SumENum=SumENum+Hij0
                PosFrac=PosFrac+1.D0
            ELSE
                SumENum=SumENum-Hij0
            ENDIF
        ELSE
            IF(WSign) THEN
                PosFrac=PosFrac+1.D0
            ENDIF
        ENDIF

        RETURN

    END SUBROUTINE SumEContrib
    

!This performs a resummed FCIMC calculation, where small graphs are created from each walker at the space, and the true matrix propagated around
    SUBROUTINE ResumFciMCPar(nI,WSign,Kii,nJ,IC,Prob,VecSlot,nIExcitGen,VecInd)
        INTEGER :: IC,nI(NEl),nJ(NEl),VecSlot,VecInd
        TYPE(ExcitGenerator) :: nIExcitGen
        TYPE(HElement) :: Hij,Hjj
        LOGICAL :: WSign
        REAL*8 :: rat,Prob,Kii,Kjj,RhoMat(2,2),Vector(2)

!First, we need to explicitly create the matrix we are using
        Hij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
        Kjj=real(Hjj%v,r2)-Hii

        RhoMat(1,1)=1.D0-Tau*(Kii-DiagSft)
        RhoMat(1,2)=-Tau*real(Hij%v,r2)
        RhoMat(2,1)=RhoMat(1,2)
        RhoMat(2,2)=1.D0-Tau*(Kjj-DiagSft)

        CALL ApplyRhoMatPar(RhoMat,Vector)  !Successivly apply the rho matrix to the particle RhoApp times

        Vector(2)=Vector(2)/Prob         !Divide the final weight of the excitation by the probability of creating it

!Create new particles according to the components of GraphVec, and put them into NewVec
        CALL CreateNewPartsPar(Vector,nI,nJ,WSign,Kii,Kjj,nIExcitGen,VecSlot,VecInd)

        RETURN
    END SUBROUTINE ResumFciMCPar

!This routine creates new particles from the vector which results from the 
    SUBROUTINE CreateNewPartsPar(Vector,nI,nJ,WSign,Kii,Kjj,nIExcitGen,VecSlot,VecInd)
        IMPLICIT NONE
        LOGICAL :: WSign,TempSign
        TYPE(ExcitGenerator) :: nIExcitGen
        INTEGER :: nI(NEl),nJ(NEl),VecInd
        INTEGER :: i,j,VecSlot,Create,ExcitLevel,iGetExcitLevel
        TYPE(HElement) :: HOffDiag
        REAL*8 :: Ran2,rat,Vector(2),Kii,Kjj

!First deal with root
        Create=INT(abs(Vector(1)))

        rat=abs(GraphVec(1))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
        IF(rat.gt.Ran2(Seed)) Create=Create+1
        
        IF((abs(Create)).gt.0) THEN
        
            IF(.not.WSign) Create=-Create
            IF(Vector(1).lt.0.D0) Create=-Create

!Test since the root should not change sign - comment out later
            IF(WSign.and.(Create.lt.0)) THEN
                WRITE(6,*) "Root determinant should not change sign"
                CALL MPIStopAll(3)
            ELSEIF((.not.WSign).and.(Create.gt.0)) THEN
                WRITE(6,*) "Root determinant should not change sign"
                CALL MPIStopAll(3)
            ENDIF
            
            IF(Create.lt.0) THEN
                TempSign=.false.
            ELSE
                TempSign=.true.
            ENDIF

    !Now actually create the particles in NewDets and NewSign
            do j=1,abs(Create)
                NewDets(:,VecSlot)=nI(:)
                NewSign(VecSlot)=TempSign
    !Copy excitation generator accross
                CALL CopyExitgenPar(nIExcitGen,NewExcits(VecSlot))
                NewIC(VecSlot)=CurrentIC(VecInd)
                NewH(:,VecSlot)=CurrentH(:,VecInd)
                VecSlot=VecSlot+1
            enddo

        ENDIF

!Now do the same for the excitation...
        Create=INT(abs(Vector(2)))

        rat=abs(GraphVec(2))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
        IF(rat.gt.Ran2(Seed)) Create=Create+1
        
        IF((abs(Create)).gt.0) THEN

            IF(.not.WSign) Create=-Create
            IF(Vector(2).lt.0.D0) Create=-Create
        
            IF(Create.lt.0) THEN
                TempSign=.false.
            ELSE
                TempSign=.true.
            ENDIF

!Calculate excitation level, connection to HF and diagonal ham element
            ExcitLevel=iGetExcitLevel(HFDet,nJ,NEl)
            IF(ExcitLevel.eq.2) THEN
!Only need off-diag conn for double excitations, since these are the only ones which contribute to energy
                HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
            ENDIF
            
            do j=1,abs(Create)
!Copy across children - cannot copy excitation generators, as do not know them
                NewDets(:,VecSlot)=nJ(:)
                NewSign(VecSlot)=TempSign
                NewExcits(VecSlot)%ExitGenForDet=.false.
                NewIC(VecSlot)=ExcitLevel
                NewH(1,VecSlot)=Kjj                     !Diagonal H-element-Hii
                NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                VecSlot=VecSlot+1
            enddo

        ENDIF

        RETURN

    END SUBROUTINE CreateNewPartsPar

!This applies the rho matrix successive times to a root determinant. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
    SUBROUTINE ApplyRhoMatPar(RhoMat,FinalVec)
        IMPLICIT NONE
        REAL*8 :: RhoMat(2,2),FinalVec(2),TempVec(2)
        INTEGER :: i,j,k
        
        FinalVec(1)=1.D0    !Set the initial vector to be 1 at the root (i.e. for one walker initially)
        FinalVec(2)=0.D0

        do i=1,RhoApp   !Cycle over the number of times we want to apply the rho matrix

!            CALL DGEMV('n',Components,Components,1.D0,GraphRhoMat,Components,GraphVec,1,0.D0,TempVec,1)
!            CALL DCOPY(Components,TempVec,1,GraphVec,1)
!            CALL AZZERO(TempVec,Components)

            TempVec(1)=0.D0
            TempVec(2)=0.D0
            
            TempVec(1)=TempVec(1)+GraphRhoMat(1,1)*FinalVec(1)
            TempVec(1)=TempVec(1)+GraphRhoMat(1,2)*FinalVec(2)
            
            TempVec(2)=TempVec(2)+GraphRhoMat(2,1)*FinalVec(1)
            TempVec(2)=TempVec(2)+GraphRhoMat(2,2)*FinalVec(2)
            
            FinalVec(1)=TempVec(1)
            FinalVec(2)=TempVec(2)
            
        enddo

        RETURN
    END SUBROUTINE ApplyRhoMatPar


!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
    SUBROUTINE CalcNewShift()
        INTEGER :: error,rc
        INTEGER :: inpair(2),outpair(2)
        REAL*8 :: TempSumNoatHF

!This first call will calculate the GrowRate for each processor, taking culling into account
        CALL UpdateDiagSftPar()

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!We need to collate the information from the different processors
        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)

        CALL MPI_Reduce(SumWalkersCyc,AllSumWalkersCyc,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        
!We want to calculate the mean growth rate over the update cycle, weighted by the total number of walkers
        GrowRate=GrowRate*SumWalkersCyc                    
        CALL MPIDSumRoot(GrowRate,1,AllGrowRate,Root)   
        IF(iProcIndex.eq.Root) THEN
            AllGrowRate=AllGrowRate/(real(AllSumWalkersCyc,r2))
        ENDIF

!Do the same for the mean excitation level of all walkers, and the total positive particles
!MeanExcitLevel here is just the sum of all the excitation levels - it needs to be divided by the total walkers in the update cycle first.
        MeanExcitLevel=(MeanExcitLevel/(real(SumWalkersCyc,r2)))
        CALL MPIDSumRoot(MeanExcitLevel,1,AllMeanExcitLevel,Root)
        IF(iProcIndex.eq.Root) THEN
            AllMeanExcitLevel=AllMeanExcitLevel/real(nProcessors,r2)
        ENDIF

        PosFrac=real(PosFrac,r2)/real(SumWalkersCyc,r2)
        CALL MPI_Reduce(PosFrac,AllPosFrac,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        IF(iProcIndex.eq.Root) THEN
            AllPosFrac=AllPosFrac/real(nProcessors,r2)
        ENDIF

!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,r2)
        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)

!To find minimum and maximum excitation levels, search for them using MPI_Reduce
        inpair(1)=MaxExcitLevel
        inpair(2)=iProcIndex

        CALL MPI_Reduce(inpair,outpair,1,MPI_2INTEGER,MPI_MAXLOC,Root,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in finding max excitation level"
            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
        ENDIF
!Max Excit Level is found on processor outpair(2) and is outpair(1)
        IF(iProcIndex.eq.Root) THEN
            AllMaxExcitLevel=outpair(1)
        ENDIF

        inpair(1)=MinExcitLevel
        inpair(2)=iProcIndex
        CALL MPI_Reduce(inpair,outpair,1,MPI_2INTEGER,MPI_MINLOC,Root,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in finding min excitation level"
            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
        ENDIF
        IF(iProcIndex.eq.Root) THEN
            AllMinExcitLevel=outpair(1)
        ENDIF

!We now want to find how the shift should change for the entire ensemble of processors
        IF(iProcIndex.eq.Root) THEN
            DiagSft=DiagSft-(log(AllGrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
            ProjectionE=AllSumENum/AllSumNoatHF
        ENDIF
!We wan to now broadcast this new shift to all processors
        CALL MPI_Bcast(DiagSft,1,MPI_DOUBLE_PRECISION,Root,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in broadcasting new shift"
            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
        ENDIF

        IF(iProcIndex.eq.Root) THEN
!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
            WRITE(15,"(I12,G16.7,I9,G16.7,I12,2G16.7,G16.7,2I6)") Iter,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,ProjectionE,AllPosFrac,AllMeanExcitLevel,AllMaxExcitLevel,AllMinExcitLevel
            WRITE(6,"(I12,G16.7,I9,G16.7,I12,2G16.7,G16.7,2I6)") Iter,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,ProjectionE,AllPosFrac,AllMeanExcitLevel,AllMaxExcitLevel,AllMinExcitLevel
            CALL FLUSH(15)
            CALL FLUSH(6)
        ENDIF

!Now need to reinitialise all variables on all processers
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        MeanExcitLevel=0.D0
        SumWalkersCyc=0
        PosFrac=0.D0
!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllSumNoatHF=0.D0
        AllTotWalkersOld=AllTotWalkers
        AllTotWalkers=0
        AllGrowRate=0.D0
        AllMeanExcitLevel=0.D0
        AllPosFrac=0.D0
        AllSumWalkersCyc=0

        RETURN
    END SUBROUTINE CalcNewShift


!This routine acts as a thermostat for the simulation - killing random particles if the population becomes too large, or 
!Doubling them if it gets too low...
    SUBROUTINE ThermostatParticlesPar(HighLow)
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

    END SUBROUTINE ThermostatParticlesPar


!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!value of the diagonal shift in the hamiltonian in order to compensate for this
    SUBROUTINE UpdateDiagSftPar()
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
    
!This is needed since the steps between culling is stored cumulatively
                GrowthSteps=CullInfo(j,3)-CullInfo(j-1,3)
                GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((CullInfo(j,1)+0.D0)/(CullInfo(j-1,2)+0.D0))

            enddo

            GrowthSteps=StepsSft-CullInfo(NoCulls,3)
            GrowRate=GrowRate+((GrowthSteps+0.D0)/(StepsSft+0.D0))*((TotWalkers+0.D0)/(CullInfo(NoCulls,2)+0.D0))

            NoCulls=0
            CALL IAZZERO(CullInfo,30)

        ENDIF
        
!        DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
!        IF((DiagSft).gt.0.D0) THEN
!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!            STOP
!        ENDIF

    END SUBROUTINE UpdateDiagSftPar


!This routine calculates the MP1 eigenvector, and uses it as a guide for setting the initial walker configuration
!    SUBROUTINE StartWavevectorPar(WaveType)
!        USE Calc , only : i_P
!        USE System , only : Beta
!        USE Integrals , only : nTay
!        IMPLICIT NONE
!        INTEGER :: ierr,i,j,WaveType,EigenvectorTag=0,k,VecSlot,NoDoublesWalk
!        CHARACTER(len=*), PARAMETER :: this_routine='StartWavevectorPar'
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
!        DEALLOCATE(Eigenvector)
!        CALL LogMemDealloc(this_routine,EigenvectorTag)
!        DEALLOCATE(ExcitStore)
!        CALL LogMemDealloc(this_routine,ExcitStoreTag)
!
!        RETURN
!
!    END SUBROUTINE StartWavevectorPar


!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        USE DetCalc , only : NDet
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot
        TYPE(HElement) :: rh,TempHii
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'

        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors

        IF(HElementSize.gt.1) THEN
            CALL STOPGM("FCIMCPar","FciMCPar cannot function with complex orbitals.")
        ENDIF
        
        IF(iProcIndex.eq.Root) THEN
            OPEN(15,file='FCIMCStats',status='unknown')
        ENDIF

!Store information specifically for the HF determinant
        ALLOCATE(HFDet(NEl),stat=ierr)
        CALL LogMemAlloc('HFDet',NEl,4,this_routine,HFDetTag)
        do i=1,NEl
            HFDet(i)=FDet(i)
        enddo

!Setup excitation generator for the HF determinant
        CALL SetupExcitgenPar(HFDet,HFExcit)

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=G_VMC_Seed-iProcIndex
        WRITE(8,*) iProcIndex

!Calculate Hii
        TempHii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
        Hii=REAL(TempHii%v,r2)

!Initialise variables for calculation on each node
        ProjectionE=0.D0
        PosFrac=0.D0
        SumENum=0.D0
        SumNoatHF=0
        MeanExcitLevel=0.D0
        MinExcitLevel=NEl+10
        MaxExcitLevel=0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllSumNoatHF=0.D0
        AllGrowRate=0.D0
        AllMeanExcitLevel=0.D0
        AllSumWalkersCyc=0
        AllPosFrac=0.D0

!Initialise global variables for calculation on the root node
        IF(iProcIndex.eq.root) THEN
            AllTotWalkers=InitWalkers*nProcessors
            AllTotWalkersOld=InitWalkers*nProcessors
        ENDIF

        IF(TResumFciMC) THEN
            IF(iProcIndex) THEN
                WRITE(6,*) "Resumming in multiple transitions to/from each excitation - only ONE excitation supported in each graph"
            ENDIF
        ENDIF
        WRITE(6,*) ""
        WRITE(6,*) "Performing FCIMC...."
        WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
        WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
        
        IF(TStartMP1) THEN
!Start the initial distribution off at the distribution of the MP1 eigenvector

            WRITE(6,*) "This feature not ready yet"
            CALL FLUSH(6)
            CALL MPIStopAll(1)
            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."
            CALL StartWavevectorPar(1)

        ELSE
!initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkers=MemoryFac*InitWalkers

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

            do j=1,InitWalkers
                CurrentDets(:,j)=HFDet(:)
                CurrentSign(j)=.true.
                CurrentIC(j)=0
            enddo

            ALLOCATE(WalkVecExcits(MaxWalkers),stat=ierr)
            ALLOCATE(WalkVec2Excits(MaxWalkers),stat=ierr)
            IF(ierr.ne.0) CALL STOPGM("InitFCIMMCCalcPar","Error in allocating walker excitation generators")

!Allocate pointers to the correct excitation arrays
            CurrentExcits=>WalkVecExcits
            NewExcits=>WalkVec2Excits

            do j=1,InitWalkers
!Copy the HF excitation generator accross to each initial particle
                CALL CopyExitGenPar(HFExcit,CurrentExcits(j))
            enddo

        ENDIF

!TotWalkers contains the number of current walkers at each step
        TotWalkers=InitWalkers
        TotWalkersOld=InitWalkers

        CALL IAZZERO(CullInfo,30)
        NoCulls=0

        RETURN

    END SUBROUTINE InitFCIMCCalcPar

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,WSign,nJ,Prob,IC)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate
        LOGICAL :: WSign
        REAL*8 :: Prob,Ran2,rat
        TYPE(HElement) :: rh

!Calculate off diagonal hamiltonian matrix element between determinants
        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)

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
                    AttemptCreatePar=-1     !-ve walker created
                ELSE
                    AttemptCreatePar=1      !+ve walker created
                ENDIF

            ELSE
!Parent particle is negative
                IF(real(rh%v).gt.0.D0) THEN
                    AttemptCreatePar=1      !+ve walker created
                ELSE
                    AttemptCreatePar=-1     !-ve walker created
                ENDIF
            ENDIF

        ELSE
!No child particle created
            AttemptCreatePar=0
        ENDIF

        IF(ExtraCreate.ne.0) THEN
!Need to include the definitely create additional particles from a initial probability > 1

            IF(AttemptCreatePar.lt.0) THEN
!In this case particles are negative
                AttemptCreatePar=AttemptCreatePar-ExtraCreate
            ELSEIF(AttemptCreatePar.gt.0) THEN
!Include extra positive particles
                AttemptCreatePar=AttemptCreatePar+ExtraCreate
            ELSEIF(AttemptCreatePar.eq.0) THEN
!No particles were stochastically created, but some particles are still definatly created - we need to determinant their sign...
                IF(WSign) THEN
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreatePar=-1*ExtraCreate    !Additional particles are negative
                    ELSE
                        AttemptCreatePar=ExtraCreate       !Additional particles are positive
                    ENDIF
                ELSE
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreatePar=ExtraCreate
                    ELSE
                        AttemptCreatePar=-1*ExtraCreate
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        RETURN

    END FUNCTION AttemptCreatePar

!This function tells us whether we should kill the particle at determinant DetCurr
!If also diffusing, then we need to know the probability with which we have spawned. This will reduce the death probability.
!The function allows multiple births(if +ve shift) or deaths from the same particle.
!The returned number is the number of deaths if positive, and the number of births if negative.
    INTEGER FUNCTION AttemptDiePar(DetCurr,Kii)
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

        AttemptDiePar=iKill

        RETURN

    END FUNCTION AttemptDiePar

!This routine copies an excitation generator from origExcit to NewExit, if the original claims that it is for the correct determinant
    SUBROUTINE CopyExitgenPar(OrigExit,NewExit)
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
            IF(ierr.ne.0) CALL STOPGM("CopyExitgenPar","Problem with allocating memory for new excitation generator")
            NewExit%ExcitData(:)=OrigExit%ExcitData(:)
            NewExit%ExitGenForDet=.true.
        ENDIF

        RETURN

    END SUBROUTINE CopyExitgenPar

    SUBROUTINE SetupExitgenPar(nI,ExcitGen)
        TYPE(ExcitGenerator) :: ExcitGen
        INTEGER :: ierr,iMaxExcit,nExcitMemLen,nJ(NEl)
        INTEGER :: nI(NEl),nStore(6)

        IF(ExcitGen%ExitGenForDet) THEN
!The excitation generator is already allocated for the determinant in question - no need to recreate it
            IF(.not.Allocated(ExcitGen%ExcitData)) THEN
                CALL STOPGM("SetupExitgenPar","Excitation generator meant to already be set up")
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
            IF(ierr.ne.0) CALL STOPGM("SetupExcitGen","Problem allocating excitation generator")
            ExcitGen%ExcitData(1)=0
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGen%ExcitData,nJ,iMaxExcit,0,nStore,3)

!Indicate that the excitation generator is now correctly allocated.
            ExcitGen%ExitGenForDet=.true.
        
        ENDIF

    END SUBROUTINE SetupExitgenPar

END MODULE FciMCParMod

#else

MODULE FciMCParMod
!Dummy module so we can use it in serial - contains all global variables
    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr
    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    USE Calc , only : TStartMP1
    USE Calc , only : GrowMaxFactor,CullFactor
    USE Calc , only : RhoApp,TResumFCIMC
    USE Determinants , only : FDet,GetHElement2
    USE DetCalc , only : NMRKS
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    USE Parallel
    IMPLICIT NONE
    SAVE

    INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
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
    INTEGER :: MemoryFac=3000

    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps
    INTEGER :: exFlag=3

!This is information needed by the thermostating, so that the correct change in walker number can be calculated, and hence the correct shift change.
!NoCulls is the number of culls in a given shift update cycle for each variable
    INTEGER :: NoCulls=0
!CullInfo is the number of walkers before and after the cull (elements 1&2), and the third element is the previous number of steps before this cull...
!Only 10 culls/growth increases are allowed in a given shift cycle
    INTEGER :: CullInfo(10,3)

!The following variables are calculated as per processor, but at the end of each update cycle, are combined to the root processor
    REAL*8 :: GrowRate,DieRat,ProjectionE,SumENum
    INTEGER*8 :: SumNoatHF      !This is the sum over all previous cycles of the number of particles at the HF determinant
    REAL*8 :: PosFrac           !This is the fraction of positive particles on each node
    INTEGER :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
    REAL*8 :: MeanExcitLevel    
    INTEGER :: MinExcitLevel
    INTEGER :: MaxExcitLevel

!These are the global variables, calculated on the root processor, from the values above
    REAL*8 :: AllGrowRate,AllMeanExcitLevel
    INTEGER :: AllMinExcitLevel,AllMaxExcitLevel,AllTotWalkers,AllTotWalkersOld,AllSumWalkersCyc
    REAL*8 :: AllSumNoatHF,AllSumENum,AllPosFrac

    REAL*8 :: MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector

    TYPE(HElement) :: rhii,FZero
    REAL*8 :: Hii

    REAL*8 , ALLOCATABLE :: GraphRhoMat(:,:)    !This stores the rho matrix for the graphs in resumFCIMC
    INTEGER :: GraphRhoMatTag=0

    REAL*8 , ALLOCATABLE :: GraphVec(:)         !This stores the final components for the propagated graph in ResumFCIMC
    INTEGER :: GraphVecTag=0

    INTEGER , ALLOCATABLE :: DetsinGraph(:,:)   !This stores the determinants in the graph created for ResumFCIMC
    INTEGER :: DetsinGraphTag=0

    REAL*8 :: RootExcitProb     !This is the probability of generating an excitation from the current root in the ResumFCIMC current graph.

END MODULE FciMCParMod
    
#endif
