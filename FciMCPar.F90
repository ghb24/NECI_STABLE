#ifdef PARALLEL
!This is a parallel MPI version of the FciMC code. It picks random excitations to work with, and can fully diagonalise the 
!resultant 2v graph, or apply it many times.
!Excitation generators are now stored along with the particles (but in a separate array)
!All variables refer to values per processor

!!!!! TO DO !!!!
! Package up variables into one array to communicated to avoid latency overheads
! Test different hashes
! Option for calculating non-essential information
! Other optimisation (Annihilation - do we want to be allocating?)


MODULE FciMCParMod
    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr
    use CalcData , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    use CalcData , only : TStartMP1,NEquilSteps,TReadPops
    use CalcData , only : GrowMaxFactor,CullFactor,TStartSinglePart
    use CalcData , only : NDets,RhoApp,TResumFCIMC,TNoAnnihil
    USE Determinants , only : FDet,GetHElement2
    USE DetCalc , only : NMRKS,ICILevel
    use IntegralsData , only : fck,NMax,UMat
    USE global_utilities
    USE HElem
    USE Parallel
    IMPLICIT NONE
    SAVE

    INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
    INTEGER , PARAMETER :: r2=kind(0.d0)
    INTEGER , PARAMETER :: i2=SELECTED_INT_KIND(18)
    
    TYPE ExcitGenerator
        INTEGER , ALLOCATABLE :: ExcitData(:)      !This stores the excitation generator
        INTEGER :: nExcitMemLen                    !This is the length of the excitation generator
        LOGICAL :: ExitGenForDet=.false.           !This is true when the excitation generator stored corresponds to the determinant
    END TYPE
    
    TYPE(ExcitGenerator) , ALLOCATABLE , TARGET :: WalkVecExcits(:),WalkVec2Excits(:)   !This will store the excitation generators for the particles on each node
    INTEGER , ALLOCATABLE , TARGET :: WalkVecDets(:,:),WalkVec2Dets(:,:)                !Contains determinant list
    LOGICAL , ALLOCATABLE , TARGET :: WalkVecSign(:),WalkVec2Sign(:)                    !Contains sign list
    INTEGER , ALLOCATABLE , TARGET :: WalkVecIC(:),WalkVec2IC(:)                        !Contains excit level list
    REAL(KIND=r2) , ALLOCATABLE , TARGET :: WalkVecH(:),WalkVec2H(:)                    !Diagonal hamiltonian element
    INTEGER , ALLOCATABLE :: IndexTable(:),Index2Table(:)                               !Indexing for the annihilation
    INTEGER , ALLOCATABLE :: ProcessVec(:),Process2Vec(:)                               !Index for process rank of original walker
    INTEGER(KIND=i2) , ALLOCATABLE :: HashArray(:),Hash2Array(:)                         !Hashes for the walkers when annihilating
    LOGICAL , ALLOCATABLE :: TempSign(:)                                                         !Temp array to hold sign of walkers when annihilating
    INTEGER(KIND=i2) , ALLOCATABLE :: TempHash(:)
    
    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0
    INTEGER :: WalkVecICTag=0,WalkVec2ICTag=0,WalkVecHTag=0,WalkVec2HTag=0
    INTEGER :: HashArrayTag=0,Hash2ArrayTag=0,IndexTableTag=0,Index2TableTag=0,ProcessVecTag=0,Process2VecTag=0
    INTEGER :: TempSignTag=0,TempHashTag=0

!Pointers to point at the correct arrays for use
    INTEGER , POINTER :: CurrentDets(:,:), NewDets(:,:)
    LOGICAL , POINTER :: CurrentSign(:), NewSign(:)
    INTEGER , POINTER :: CurrentIC(:), NewIC(:)
    REAL*8 , POINTER :: CurrentH(:), NewH(:)
    TYPE(ExcitGenerator) , POINTER :: CurrentExcits(:), NewExcits(:)
    
    INTEGER , ALLOCATABLE :: HFDet(:)       !This will store the HF determinant
    INTEGER :: HFDetTag=0
    TYPE(ExcitGenerator) :: HFExcit         !This is the excitation generator for the HF determinant
    INTEGER(KIND=i2) :: HFHash               !This is the hash for the HF determinant

!MemoryFac is the factor by which space will be made available for extra walkers compared to InitWalkers
    INTEGER :: MemoryFac=30

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
    INTEGER :: Annihilated      !This is the number annihilated on one processor
    INTEGER :: NoatHF           !This is the instantaneous number of particles at the HF determinant
    INTEGER :: NoatDoubs
    INTEGER :: Acceptances      !This is the number of accepted spawns - this is only calculated per node.
    REAL*8 :: AccRat            !Acceptance ratio for each node over the update cycle

!These are the global variables, calculated on the root processor, from the values above
    REAL*8 :: AllGrowRate,AllMeanExcitLevel
    INTEGER :: AllMinExcitLevel,AllMaxExcitLevel,AllTotWalkers,AllTotWalkersOld,AllSumWalkersCyc,AllTotSign,AllTotSignOld
    INTEGER :: AllAnnihilated,AllNoatHF,AllNoatDoubs
    REAL*8 :: AllSumNoatHF,AllSumENum,AllPosFrac

    REAL*8 :: MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector

    TYPE(HElement) :: rhii,FZero
    REAL*8 :: Hii

    REAL*8 , ALLOCATABLE :: GraphRhoMat(:,:)    !This stores the rho matrix for the graphs in resumFCIMC
    INTEGER :: GraphRhoMatTag=0

    REAL*8 , ALLOCATABLE :: GraphVec(:)         !This stores the final components for the propagated graph in ResumFCIMC
    INTEGER :: GraphVecTag=0

    REAL*8 , ALLOCATABLE :: GraphKii(:)         !This stores the diagonal Kii matrix elements for the determinants in the graph
    INTEGER :: GraphKiiTag=0

    INTEGER , ALLOCATABLE :: DetsinGraph(:,:)   !This stores the determinants in the graph created for ResumFCIMC
    INTEGER :: DetsinGraphTag=0

    LOGICAL :: TSinglePartPhase                 !This is true if TStartSinglePart is true, and we are still in the phase where the shift is fixed and particle numbers are growing

    INTEGER :: mpilongintegertype               !This is used to create an MPI derived type to cope with 8 byte integers

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j,error
        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'
        TYPE(HElement) :: Hamii

        CALL InitFCIMCCalcPar()

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(6,*) ""
            WRITE(6,*) "       Step     Shift      WalkerCng    GrowRate      TotWalkers   Annihil    Proj.E        NoatHF NoatDoubs  +veWalkFrac       AccRat   MeanEx     MinEx MaxEx"
            WRITE(15,*) "#       Step     Shift      WalkerCng    GrowRate      TotWalkers  Annihil    Proj.E        NoatHF NoatDoubs   +veWalkFrac       AccRat   MeanEx     MinEx MaxEx"
!TotWalkersOld is the number of walkers last time the shift was changed
            WRITE(15,"(I12,G16.7,I9,G16.7,I12,I8,G16.7,2I10,F13.5,2G13.5,2I6)") Iter,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,AllAnnihilated,ProjectionE,AllNoatHF,AllNoatDoubs,1.D0,AccRat,AllMeanExcitLevel,AllMaxExcitLevel,AllMinExcitLevel
            WRITE(6,"(I12,G16.7,I9,G16.7,I12,I8,G16.7,2I10,F13.5,2G13.5,2I6)") Iter,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,AllAnnihilated,ProjectionE,AllNoatHF,AllNoatDoubs,1.D0,AccRat,AllMeanExcitLevel,AllMaxExcitLevel,AllMinExcitLevel
        ENDIF
        

!Start MC simulation...
        do Iter=1,NMCyc
            
            CALL PerformFCIMCycPar()

            IF(mod(Iter,StepsSft).eq.0) THEN
!This will communicate between all nodes, find the new shift (and other parameters) and broadcast them to the other nodes.
                CALL CalcNewShift()
            ENDIF
            
!End of MC cycle
        enddo

        Weight=HDElement(0.D0)
        Energyxw=HDElement(ProjectionE)

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
        
        IF(TResumFCIMC.and.(NDets.gt.2)) THEN
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

        IF(.not.TNoAnnihil) THEN
!Free the mpi derived type that we have created for the hashes.
            CALL MPI_Type_free(mpilongintegertype,error)
        ENDIF

        IF(iProcIndex.eq.Root) CLOSE(15)

!        CALL MPIEnd(.false.)    !Finalize MPI

        RETURN

    END SUBROUTINE FciMCPar

!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCycPar()
        INTEGER :: VecSlot,i,j,k,l
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount
        REAL*8 :: Prob,rat,HDiag
        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,iGetExcitLevel_2,TotWalkersNew
        LOGICAL :: WSign
        INTEGER(KIND=i2) :: HashTemp
        TYPE(HElement) :: HDiagTemp,HOffDiag
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        
        do j=1,TotWalkers
!j runs through all current walkers

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(CurrentDets(:,j),CurrentIC(j),CurrentSign(j))

!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
            CALL SetupExitgenPar(CurrentDets(:,j),CurrentExcits(j))
            
            IF(TResumFCIMC) THEN

                IF(NDets.eq.2) THEN

                    CALL GenRandSymExcitIt3(CurrentDets(:,j),CurrentExcits(j)%ExcitData,nJ,Seed,IC,0,Prob,iCount)
                    CALL ResumFciMCPar(CurrentDets(:,j),CurrentSign(j),CurrentH(j),nJ,IC,Prob,VecSlot,CurrentExcits(j),j)

                ELSE

                    CALL ResumGraphPar(CurrentDets(:,j),CurrentSign(j),VecSlot,CurrentExcits(j),j)

                ENDIF

            ELSE

                CALL GenRandSymExcitIt3(CurrentDets(:,j),CurrentExcits(j)%ExcitData,nJ,Seed,IC,0,Prob,iCount)
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
                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
!                    IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
!                        HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
!                    ENDIF
                    IF(ExcitLevel.eq.0) THEN
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                    ELSE
                        HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                        HDiag=(REAL(HDiagTemp%v,r2))-Hii
                    ENDIF

                    IF(.not.TNoAnnihil) THEN
                        HashTemp=CreateHash(nJ)
                    ENDIF

                    do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them
                        NewDets(:,VecSlot)=nJ(:)
                        NewSign(VecSlot)=WSign
                        NewExcits(VecSlot)%ExitGenForDet=.false.
                        NewIC(VecSlot)=ExcitLevel
                        NewH(VecSlot)=HDiag                     !Diagonal H-element-Hii
!                        NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                        IF(.not.TNoAnnihil) THEN
                            Hash2Array(VecSlot)=HashTemp
                        ENDIF
                        VecSlot=VecSlot+1
                    enddo

                    Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
                
                ENDIF   !End if child created

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
                iDie=AttemptDiePar(CurrentDets(:,j),CurrentH(j))
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
                        NewH(VecSlot)=CurrentH(j)
                        IF(.not.TNoAnnihil) THEN
                            Hash2Array(VecSlot)=HashArray(j)
                        ENDIF
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
                        NewH(VecSlot)=CurrentH(j)
                        IF(.not.TNoAnnihil) THEN
                            Hash2Array(VecSlot)=HashArray(j)
                        ENDIF
                        VecSlot=VecSlot+1
                    enddo
            
                ENDIF   !To kill if

            ENDIF   !Resum if

!Finish cycling over walkers
        enddo
        
!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+TotWalkers

!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1
        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
            CALL FLUSH(6)
        ENDIF
        IF(TotWalkersNew.eq.0) THEN
            CALL Stop_All("PerformFciMCycPar","All walkers have died on a node")
        ENDIF
        
        IF(TNoAnnihil) THEN
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

            TotWalkers=TotWalkersNew

        ELSE
!This routine now cancels down the particles with opposing sign on each determinant

            CALL AnnihilatePartPar(TotWalkersNew)
            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)

        ENDIF
        
!Find the total number of particles at HF (x sign) across all nodes. If this is negative, flip the sign of all particles.
        AllNoatHF=0
!Find sum of noathf, and then use an AllReduce to broadcast it to all nodes
        CALL MPI_AllReduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        IF(AllNoatHF.lt.0) THEN
!Flip the sign if we're beginning to get a negative population on the HF
            WRITE(6,*) "No. at HF < 0 - flipping sign of entire ensemble of particles..."
            CALL FlipSign()
        ENDIF

        IF(TSinglePartPhase) THEN
!Do not allow culling if we are still in the single particle phase.
            IF(iProcIndex.eq.root) THEN     !Only exit phase if particle number is sufficient on head node.
                IF(TotWalkers.gt.InitWalkers) THEN
                    WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
                    TSinglePartPhase=.false.
                ENDIF
            ENDIF
!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration
            CALL MPI_Bcast(TSinglePartPhase,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
            IF(ierr.ne.MPI_SUCCESS) THEN
                WRITE(6,*) "Problem in broadcasting new TSinglePartPhase"
                CALL Stop_All("PerformFciMCycPar","Problem in broadcasting new TSinglePartPhase")
            ENDIF
        ELSE

            IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!Particle number is too large - kill them randomly

!Log the fact that we have made a cull
                NoCulls=NoCulls+1
                IF(NoCulls.gt.10) THEN
                    WRITE(6,*) "Too Many Culls"
                    CALL FLUSH(6)
                    call Stop_All("PerformFCIMCyc","Too Many Culls")
                ENDIF
!CullInfo(:,1) is walkers before cull
                CullInfo(NoCulls,1)=TotWalkers
!CullInfo(:,3) is MC Steps into shift cycle before cull
                CullInfo(NoCulls,3)=mod(Iter,StepsSft)

                WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number on this node..."
                WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number on this node..."
                WRITE(6,"(A,F8.2)") "Population on this node will reduce by a factor of ",CullFactor
                CALL ThermostatParticlesPar(.true.)

            ELSEIF(TotWalkers.lt.(InitWalkers/2)) THEN
!Particle number is too small - double every particle in its current position

!Log the fact that we have made a cull
                NoCulls=NoCulls+1
                IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!CullInfo(:,1) is walkers before cull
                CullInfo(NoCulls,1)=TotWalkers
!CullInfo(:,3) is MC Steps into shift cycle before cull
                CullInfo(NoCulls,3)=mod(Iter,StepsSft)
                
                WRITE(6,*) "Doubling particle population on this node to increase total number..."
                CALL ThermostatParticlesPar(.false.)

            ENDIF
        
        ENDIF

        RETURN

    END SUBROUTINE PerformFCIMCycPar

        
    SUBROUTINE AnnihilatePartPar(TotWalkersNew)
        INTEGER :: i,j,k,ToAnnihilateIndex,TotWalkersNew,ierr,error,sendcounts(nProcessors)
        INTEGER :: TotWalkersDet,InitialBlockIndex,FinalBlockIndex,ToAnnihilateOnProc,VecSlot
        INTEGER :: disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors),MaxIndex
        INTEGER :: Minsendcounts,Maxsendcounts,DebugIter
        INTEGER(KIND=i2) :: HashCurr!MinBin,RangeofBins,NextBinBound
        CHARACTER(len=*), PARAMETER :: this_routine='AnnihilatePartPar'

        DebugIter=0
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "Printing out annihilation debug info for Iteration: ",Iter,DebugIter
        ENDIF

!First, allocate memory to hold the signs and the hashes while we annihilate
        ALLOCATE(TempSign(TotWalkersNew),stat=ierr)
!Comment out the memallocs later
        CALL LogMemAlloc('TempSign',TotWalkersNew,4,this_routine,TempSignTag,ierr)
        ALLOCATE(TempHash(TotWalkersNew),stat=ierr)
        CALL LogMemAlloc('TempHash',TotWalkersNew,8,this_routine,TempHashTag,ierr)
        
!Temporary arrays, storing the signs and Hashes need ot be kept, as both these arrays are going to be mixed
        do i=1,TotWalkersNew
!I can probably use icopy/dcopy here instead
            TempSign(i)=NewSign(i)
            TempHash(i)=Hash2Array(i)
        enddo
    
!Create the arrays for index and process
        do i=1,TotWalkersNew
            IndexTable(i)=i
            ProcessVec(i)=iProcIndex
        enddo

        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) TotWalkersNew
            do i=1,TotWalkersNew
                WRITE(6,*) i,Hash2Array(i),IndexTable(i),ProcessVec(i),NewSign(i)
            enddo
        ENDIF

!Next, order the hash array, taking the index, CPU and sign with it...
!Order the array by mod(Hash,nProcessors). This will result in a more load-balanced system
        CALL SortMod3I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew),nProcessors)
!        CALL Sort3I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew))
        
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "***************"
            WRITE(6,*) TotWalkersNew
            do i=1,TotWalkersNew
                WRITE(6,*) Hash2Array(i),mod(Hash2Array(i),nProcessors),IndexTable(i),ProcessVec(i),NewSign(i)
            enddo
        ENDIF
        
!        IF(nProcessors.eq.1) THEN
!            CALL STOP_All("AnnihilatePartPar","One processor annihilation not available yet...")
!        ENDIF
 
!Create the send counts and disps for the AlltoAllv. Work out equal ranges of bins for the hashes
!        Rangeofbins=HUGE(Rangeofbins)/(nProcessors/2)
!        MinBin=HUGE(MinBin)*-1
!        NextBinBound=MinBin+Rangeofbins

!Send counts is the size of each block of ordered dets which are going to each processor. This could be binary searched for extra speed
        j=1
        do i=0,nProcessors-1    !Search through all possible values of mod(Hash,nProcessors)

            do while((mod(Hash2Array(j),nProcessors).eq.i).and.(j.le.TotWalkersNew))
                j=j+1 
            enddo
            sendcounts(i+1)=j-1
!            IF(i.eq.nProcessors) THEN
!                NextBinBound=HUGE(NextBinBound)
!            ELSE
!                NextBinBound=NextBinBound+Rangeofbins
!            ENDIF

        enddo
        
        IF(sendcounts(nProcessors).ne.TotWalkersNew) THEN
            CALL Stop_All("AnnihilatePartPar","Incorrect calculation of sendcounts")
        ENDIF

!Oops, we have calculated them cumulativly - undo this
        maxsendcounts=sendcounts(1)
        minsendcounts=sendcounts(1)     !Find max & min sendcounts, so that load-balancing can be checked
!        WRITE(6,*) maxsendcounts,minsendcounts
        do i=2,nProcessors
            do j=1,i-1
                sendcounts(i)=sendcounts(i)-sendcounts(j)
            enddo
            IF(sendcounts(i).gt.maxsendcounts) THEN
                maxsendcounts=sendcounts(i)
            ELSEIF(sendcounts(i).lt.minsendcounts) THEN
                minsendcounts=sendcounts(i)
            ENDIF
        enddo
        
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!The disps however do want to be cumulative - this is the array indexing the start of the data block
        disps(1)=0      !Starting element is always the first element
        do i=2,nProcessors
            disps(i)=disps(i-1)+sendcounts(i-1)
        enddo
        
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "SENDCOUNTS: "
            WRITE(6,*) sendcounts(:)
            WRITE(6,*) "DISPS: "
            WRITE(6,*) disps(:)
            CALL FLUSH(6)
        ENDIF

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        CALL IAZZERO(recvcounts,nProcessors)
!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.95*MaxWalkers)) THEN
            CALL Stop_All("AnnihilatePartPar","Maximum index is close to maximum length of array")
        ENDIF

        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "RECVCOUNTS: "
            WRITE(6,*) recvcounts(:)
            WRITE(6,*) "RECVDISPS: "
            WRITE(6,*) recvdisps(:),MaxIndex
            CALL FLUSH(6)
        ENDIF

!Insert a load-balance check here...maybe find the s.d. of the sendcounts array - maybe just check the range first.
        IF(TotWalkersNew.gt.200) THEN
            IF((Maxsendcounts-Minsendcounts).gt.(TotWalkersNew/3)) THEN
                WRITE(6,"(A,I12)") "**WARNING** Parallel annihilation not optimally balanced on this node, for iter = ",Iter
                WRITE(6,*) "Sendcounts is: ",sendcounts(:)
                CALL FLUSH(6)
            ENDIF
        ENDIF
        
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Now send the chunks of hashes to the corresponding processors
        CALL MPI_AlltoAllv(Hash2Array(1:TotWalkersNew),sendcounts,disps,MPI_DOUBLE_PRECISION,HashArray(1:MaxIndex),recvcounts,recvdisps,mpilongintegertype,MPI_COMM_WORLD,error)        
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in dividing up hashes in annihilation"
            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
        ENDIF

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPI_AlltoAllv(NewSign(1:TotWalkersNew),sendcounts,disps,MPI_LOGICAL,CurrentSign,recvcounts,recvdisps,MPI_LOGICAL,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in dividing up hashes in annihilation"
            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
        ENDIF
        CALL MPI_AlltoAllv(IndexTable(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in dividing up hashes in annihilation"
            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
        ENDIF
        CALL MPI_AlltoAllv(ProcessVec(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Process2Vec,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in dividing up hashes in annihilation"
            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
        ENDIF
        
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "AFTER DIVISION:   - No. on processor is: ",MaxIndex
            do i=1,MaxIndex
                WRITE(6,*) HashArray(i),mod(HashArray(i),nProcessors),Index2Table(i),Process2Vec(i),CurrentSign(i)
            enddo
            CALL FLUSH(6)
        ENDIF

!The hashes now need to be sorted again - this time by their number
        CALL Sort3I1LLong(MaxIndex,HashArray(1:MaxIndex),Index2Table(1:MaxIndex),Process2Vec(1:MaxIndex),CurrentSign(1:MaxIndex))
        
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "AFTER DIVISION & ORDERING:   - No. on processor is: ",MaxIndex
            do i=1,MaxIndex
                WRITE(6,*) HashArray(i),mod(HashArray(i),nProcessors),Index2Table(i),Process2Vec(i),CurrentSign(i)
            enddo
            CALL FLUSH(6)
        ENDIF

!Work out the index of the particles which want to be annihilated
        j=1
        ToAnnihilateIndex=1
        do while(j.le.MaxIndex)
            TotWalkersDet=0
            InitialBlockIndex=j
            FinalBlockIndex=j-1         !Start at j-1 since we are increasing FinalBlockIndex even with the first det in the next loop
            HashCurr=HashArray(j)
            do while((HashArray(j).eq.HashCurr).and.(j.le.MaxIndex))
!First loop counts walkers in the block
                IF(CurrentSign(j)) THEN
                    TotWalkersDet=TotWalkersDet+1
                ELSE
                    TotWalkersDet=TotWalkersDet-1
                ENDIF
                FinalBlockIndex=FinalBlockIndex+1
                j=j+1
            enddo

            IF(Iter.eq.DebugIter) THEN
                WRITE(6,*) "Common block of dets found from ",InitialBlockIndex," ==> ",FinalBlockIndex
                WRITE(6,*) "Sum of signs in block is: ",TotWalkersDet
                CALL FLUSH(6)
            ENDIF

            do k=InitialBlockIndex,FinalBlockIndex
!Second run through the block of same determinants marks walkers for annihilation
                IF(TotWalkersDet.eq.0) THEN
!All walkers in block want to be annihilated
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 1",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
                ELSEIF((TotWalkersDet.lt.0).and.(CurrentSign(k))) THEN
!Annihilate if block has a net negative walker count, and current walker is positive
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 2",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
                ELSEIF((TotWalkersDet.gt.0).and.(.not.CurrentSign(k))) THEN
!Annihilate if block has a net positive walker count, and current walker is negative
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 3",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
                ELSE
!If net walkers is positive, and we have a positive walkers, then remove one from the net positive walkers and continue through the block
                    IF(CurrentSign(k)) THEN
                        TotWalkersDet=TotWalkersDet-1
                    ELSE
                        TotWalkersDet=TotWalkersDet+1
                    ENDIF
                ENDIF
            enddo
            
!            j=j+1   !Increment counter

        enddo

        ToAnnihilateIndex=ToAnnihilateIndex-1   !ToAnnihilateIndex now tells us the total number of particles to annihilate from the list on this processor
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "Number of particles to annihilate from hashes on this processor: ",ToAnnihilateIndex
            CALL FLUSH(6)
        ENDIF

!The annihilation is complete - particles to be annihilated are stored in IndexTable and need to be sent back to their original processor
!To know which processor that is, we need to order the particles to be annihilated in terms of their CPU, i.e. ProcessVec(1:ToAnnihilateIndex)
!Is the list already ordered according to CPU? Is this further sort even necessary?

        IF(ToAnnihilateIndex.gt.1) THEN
            CALL Sort2IILongL(ToAnnihilateIndex,ProcessVec(1:ToAnnihilateIndex),IndexTable(1:ToAnnihilateIndex),Hash2Array(1:ToAnnihilateIndex),NewSign(1:ToAnnihilateIndex))
        ENDIF

!We now need to regenerate sendcounts and disps
        CALL IAZZERO(sendcounts,nProcessors)
        do i=1,ToAnnihilateIndex
            IF(ProcessVec(i).gt.(nProcessors-1)) THEN
                CALL Stop_All("AnnihilatePartPar","Annihilation error")
            ENDIF
            sendcounts(ProcessVec(i)+1)=sendcounts(ProcessVec(i)+1)+1
        enddo
!The disps however do want to be cumulative
        disps(1)=0      !Starting element is always the first element
        do i=2,nProcessors
            disps(i)=disps(i-1)+sendcounts(i-1)
        enddo

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        CALL IAZZERO(recvcounts,nProcessors)
!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        ToAnnihilateonProc=recvdisps(nProcessors)+recvcounts(nProcessors)
        
        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "FOR RETURN OF ANNIHILATED PARTICLES, SENDCOUNTS: ",sendcounts(:)
            WRITE(6,*) "DISPS: ",disps(:)
            WRITE(6,*) "RECVCOUNTS: ",recvcounts(:)
            WRITE(6,*) "RECVDISPS: ",recvdisps(:)
            WRITE(6,*) "ToAnnihilateOnProc: ",ToAnnihilateonProc
            CALL FLUSH(6)
        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Perform another matrix transpose of the annihilation data using MPI_AlltoAllv, to send the data back to its correct Processor
        CALL MPI_AlltoAllv(Hash2Array(1:TotWalkersNew),sendcounts,disps,MPI_DOUBLE_PRECISION,HashArray,recvcounts,recvdisps,mpilongintegertype,MPI_COMM_WORLD,error)        
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in sending back annihilated particles"
            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
        ENDIF

!The signs of the hashes, index and CPU also need to be taken with them. (CPU does not need to be taken - every element of CPU should be equal to the rank of the processor+1)
!Hash also does not need to be taken, but will be taken as a precaution
        CALL MPI_AlltoAllv(NewSign(1:TotWalkersNew),sendcounts,disps,MPI_LOGICAL,CurrentSign,recvcounts,recvdisps,MPI_LOGICAL,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in sending back annihilated particles"
            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
        ENDIF
        CALL MPI_AlltoAllv(IndexTable(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in sending back annihilated particles"
            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
        ENDIF
        CALL MPI_AlltoAllv(ProcessVec(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Process2Vec,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in sending back annihilated particles"
            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
        ENDIF

!TEST
        do i=1,ToAnnihilateonProc
            IF(Process2Vec(i).ne.(iProcIndex)) THEN
                CALL Stop_All("AnnihilatePartPar","AlltoAllv performed incorrectly")
            ENDIF
        enddo

!Index2Table now is a list, of length "ToAnnihilateonProc", of walkers which should NOT be transferred to the next array. 
!Order the list according to this index (Hash and sign does not need to be sorted, but will for debugging purposes)

        CALL SORTIILongL(ToAnnihilateonProc,Index2Table(1:ToAnnihilateonProc),HashArray(1:ToAnnihilateonProc),CurrentSign(1:ToAnnihilateonProc))

        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "Number of hashes originally on processor which need to be removed=",ToAnnihilateonProc
            WRITE(6,*) "To annihilate from processor: "
            do i=1,ToAnnihilateonProc
                WRITE(6,*) Index2Table(i),HashArray(i),CurrentSign(i)
            enddo
        ENDIF

!TEST - do the hashes and signs match the ones that are returned?
        do i=1,ToAnnihilateonProc
            IF(TempHash(Index2Table(i)).ne.(HashArray(i))) THEN
                CALL Stop_All("AnnihilatePartPar","Incorrect Hash returned")
            ENDIF
            IF(TempSign(Index2Table(i))) THEN
                IF(.not.CurrentSign(i)) THEN
                    CALL Stop_All("AnnihilatePartPar","Incorrect Sign returned")
                ENDIF
            ELSE
                IF(CurrentSign(i)) THEN
                    CALL Stop_All("AnnihilatePartPar","Incorrect Sign returned")
                ENDIF
            ENDIF
        enddo
        

        IF(ToAnnihilateonProc.ne.0) THEN
!Copy across the data, apart from ones which have an index given by the indicies in Index2Table(1:ToAnnihilateonProc)
            VecSlot=1       !VecSlot is the index in the final array of TotWalkers
            i=1             !i is the index in the original array of TotWalkersNew
            do j=1,ToAnnihilateonProc
!Loop over all particles to be annihilated
                IF(Iter.eq.DebugIter) WRITE(6,*) Index2Table(j)
                do while(i.lt.Index2Table(j))
!Copy accross all particles less than this number
                    CALL ICOPY(NEl,NewDets(:,i),1,CurrentDets(:,VecSlot),1)
                    CurrentIC(VecSlot)=NewIC(i)
                    CurrentH(VecSlot)=NewH(i)
                    CALL CopyExitgenPar(NewExcits(i),CurrentExcits(VecSlot))
                    HashArray(VecSlot)=TempHash(i)
                    CurrentSign(VecSlot)=TempSign(i)
                    i=i+1
                    VecSlot=VecSlot+1
                enddo
                i=i+1
            enddo

!Now need to copy accross the residual - from Index2Table(ToAnnihilateonProc) to TotWalkersNew
            do i=Index2Table(ToAnnihilateonProc)+1,TotWalkersNew
                CALL ICOPY(NEl,NewDets(:,i),1,CurrentDets(:,VecSlot),1)
                CurrentIC(VecSlot)=NewIC(i)
                CurrentH(VecSlot)=NewH(i)
                CALL CopyExitgenPar(NewExcits(i),CurrentExcits(VecSlot))
                HashArray(VecSlot)=TempHash(i)
                CurrentSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo

        ELSE
!No particles annihilated
            VecSlot=1
            do i=1,TotWalkersNew
                CALL ICOPY(NEl,NewDets(:,i),1,CurrentDets(:,VecSlot),1)
                CurrentIC(VecSlot)=NewIC(i)
                CurrentH(VecSlot)=NewH(i)
                CALL CopyExitgenPar(NewExcits(i),CurrentExcits(VecSlot))
                HashArray(VecSlot)=TempHash(i)
                CurrentSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo
        ENDIF
                
        TotWalkers=VecSlot-1

        IF(Iter.eq.DebugIter) THEN
            WRITE(6,*) "FINAL CONFIGURATION: "
            do i=1,TotWalkers
                WRITE(6,*) i,HashArray(i),CurrentSign(i)
            enddo
        ENDIF

        IF((TotWalkersNew-TotWalkers).ne.ToAnnihilateonProc) THEN
            WRITE(6,*) TotWalkers,TotWalkersNew,ToAnnihilateonProc,Iter
            CALL FLUSH(6)
            CALL Stop_All("AnnihilatePartPar","Problem with numbers when annihilating")
        ENDIF

        DEALLOCATE(TempSign)
        CALL LogMemDealloc(this_routine,TempSignTag)
        DEALLOCATE(TempHash)
        CALL LogMemDealloc(this_routine,TempHashTag)
        
        RETURN

    END SUBROUTINE AnnihilatePartPar


!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,WSign)
        INTEGER :: DetCurr(NEl),ExcitLevel
        LOGICAL :: WSign
        TYPE(HElement) :: HOffDiag

        MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,r2)
        IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
        IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
        IF(ExcitLevel.eq.0) THEN
            IF(WSign) THEN
                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+1
                NoatHF=NoatHF+1
                PosFrac=PosFrac+1.D0
            ELSE
                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF-1
                NoatHF=NoatHF-1
            ENDIF
        ELSEIF(ExcitLevel.eq.2) THEN
            NoatDoubs=NoatDoubs+1
!At double excit - find and sum in energy
            HOffDiag=GetHElement2(HFDet,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
            IF(WSign) THEN
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+REAL(HOffDiag%v,r2)
                PosFrac=PosFrac+1.D0
            ELSE
                IF(Iter.gt.NEquilSteps) SumENum=SumENum-REAL(HOffDiag%v,r2)
            ENDIF
        ELSE
            IF(WSign) THEN
                PosFrac=PosFrac+1.D0
            ENDIF
        ENDIF

        RETURN

    END SUBROUTINE SumEContrib
    
    SUBROUTINE ResumGraphPar(nI,WSign,VecSlot,nIExcitGen,VecInd)
        INTEGER :: nI(NEl),VecSlot,VecInd,ExcitLevel,iGetExcitLevel,Create,i,j
        TYPE(ExcitGenerator) :: nIExcitGen
        TYPE(HElement) :: HOffDiag
        LOGICAL :: WSign,ChildSign
        REAL*8 :: Prob,rat,Ran2

        CALL CreateGraphPar(nI,nIExcitGen,Prob,VecInd)      !Create graph with NDets distinct determinants
        CALL ApplyRhoMatPar()   !Apply the rho matrix successive times

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
            NewH(VecSlot)=CurrentH(VecInd)
            CALL CopyExitgenPar(CurrentExcits(VecInd),NewExcits(VecSlot))
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
                
                ExcitLevel=iGetExcitLevel(HFDet,DetsinGraph(:,i),NEl)
!                IF(ExcitLevel.eq.2) THEN
!Only need it for double excitations, since these are the only ones which contribute to energy
!                    HOffDiag=GetHElement2(HFDet,DetsinGraph(:,i),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                IF(ExcitLevel.eq.0) THEN
                    IF(ABS(GraphKii(i)).gt.1.D-07) THEN
                        CALL Stop_All("ResumGraphPar","Diagonal K-mat element should be zero for HF particles")
                    ENDIF
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
                    NewH(VecSlot)=GraphKii(i)       !Diagonal H El previously stored
!                    NewH(2,VecSlot)=REAL(HOffDiag%v,r2)
                    NewExcits(VecSlot)%ExitGenForDet=.false.
                    VecSlot=VecSlot+1
                enddo

            ENDIF

        enddo

        RETURN

    END SUBROUTINE ResumGraphPar

    SUBROUTINE CreateGraphPar(nI,nIExcitGen,Prob,VecInd)
        INTEGER :: nI(NEl),VecInd,nJ(NEl),iCount,IC,i,j,Attempts
        TYPE(ExcitGenerator) :: nIExcitGen
        REAL*8 :: Prob,Kii,ExcitProb
        LOGICAL :: SameDet,CompiPath
        TYPE(HElement) :: Hamij,Hamii

        CALL AZZERO(GraphRhoMat,NDets*NDets)

!Do not need to put the root determinant in the first column of DetsinGraph -
!just assume its there.
        
        Kii=CurrentH(VecInd)      !This is now the Kii element of the root
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

    END SUBROUTINE CreateGraphPar


!This performs a resummed FCIMC calculation, where two-vertex graphs are created from each walker at the space, and the true matrix propagated around
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

!        WRITE(6,*) "Rhomat is : ",RhoMat

        CALL ApplyRhoMatTwoPar(RhoMat,Vector)  !Successivly apply the rho matrix to the particle RhoApp times

!        WRITE(6,*) "Rhomat applied to get : ",Vector

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
        INTEGER :: i,j,VecSlot,Create,ExcitLevel,iGetExcitLevel_2
        TYPE(HElement) :: HOffDiag
        REAL*8 :: Ran2,rat,Vector(2),Kii,Kjj

!First deal with root
        Create=INT(abs(Vector(1)))

        rat=abs(Vector(1))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
        IF(rat.gt.Ran2(Seed)) Create=Create+1
        
        IF((abs(Create)).gt.0) THEN
        
            IF(.not.WSign) Create=-Create
            IF(Vector(1).lt.0.D0) Create=-Create

!Test since the root should not change sign - comment out later
            IF(WSign.and.(Create.lt.0)) THEN
                call Stop_All("CreateNewPartsPar","Root determinant should not change sign")
            ELSEIF((.not.WSign).and.(Create.gt.0)) THEN
                call Stop_All("CreateNewPartsPar","Root determinant should not change sign")
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
                NewH(VecSlot)=CurrentH(VecInd)
                VecSlot=VecSlot+1
            enddo

        ENDIF

!Now do the same for the excitation...
        Create=INT(abs(Vector(2)))

        rat=abs(Vector(2))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
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
            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
!            IF(ExcitLevel.eq.2) THEN
!Only need off-diag conn for double excitations, since these are the only ones which contribute to energy
!                HOffDiag=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
!            ENDIF
            
            do j=1,abs(Create)
!Copy across children - cannot copy excitation generators, as do not know them
                NewDets(:,VecSlot)=nJ(:)
                NewSign(VecSlot)=TempSign
                NewExcits(VecSlot)%ExitGenForDet=.false.
                NewIC(VecSlot)=ExcitLevel
                NewH(VecSlot)=Kjj                     !Diagonal H-element-Hii
!                NewH(2,VecSlot)=REAL(HOffDiag%v,r2)       !Off-diagonal H-element
                VecSlot=VecSlot+1
            enddo

        ENDIF

        RETURN

    END SUBROUTINE CreateNewPartsPar

!This applies the rho matrix successive times to a root determinant. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
    SUBROUTINE ApplyRhoMatPar()
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

    END SUBROUTINE ApplyRhoMatPar


!This applies the rho matrix successive times to a root determinant for two vertex graphs. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
    SUBROUTINE ApplyRhoMatTwoPar(RhoMat,FinalVec)
        REAL*8 :: RhoMat(2,2),FinalVec(2),TempVec(2)
        INTEGER :: i,j,k
        
        FinalVec(1)=1.D0    !Set the initial vector to be 1 at the root (i.e. for one walker initially)
        FinalVec(2)=0.D0

        do i=1,RhoApp   !Cycle over the number of times we want to apply the rho matrix

!            CALL DGEMV('n',Components,Components,1.D0,RhoMat,Components,GraphVec,1,0.D0,TempVec,1)
!            CALL DCOPY(Components,TempVec,1,GraphVec,1)
!            CALL AZZERO(TempVec,Components)

            TempVec(1)=(RhoMat(1,1)*FinalVec(1))+(RhoMat(1,2)*FinalVec(2))
            TempVec(2)=(RhoMat(2,1)*FinalVec(1))+(RhoMat(2,2)*FinalVec(2))
            
            FinalVec(1)=TempVec(1)
            FinalVec(2)=TempVec(2)

!            WRITE(6,*) "TempVec for i=",i," is: ",TempVec
            
        enddo

        RETURN
    END SUBROUTINE ApplyRhoMatTwoPar


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

!Find total Annihilated,Total at HF and Total at doubles
        CALL MPI_Reduce(Annihilated,AllAnnihilated,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(NoatDoubs,AllNoatDoubs,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)

!Do the same for the mean excitation level of all walkers, and the total positive particles
!MeanExcitLevel here is just the sum of all the excitation levels - it needs to be divided by the total walkers in the update cycle first.
        MeanExcitLevel=(MeanExcitLevel/(real(SumWalkersCyc,r2)))
        CALL MPIDSumRoot(MeanExcitLevel,1,AllMeanExcitLevel,Root)
        IF(iProcIndex.eq.Root) THEN
            AllMeanExcitLevel=AllMeanExcitLevel/real(nProcessors,r2)
        ENDIF

        PosFrac=PosFrac/real(SumWalkersCyc,r2)
        CALL MPI_Reduce(PosFrac,AllPosFrac,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
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
            IF(.not.TSinglePartPhase) DiagSft=DiagSft-(log(AllGrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
            ProjectionE=AllSumENum/AllSumNoatHF
        ENDIF
!We wan to now broadcast this new shift to all processors
        CALL MPI_Bcast(DiagSft,1,MPI_DOUBLE_PRECISION,Root,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            WRITE(6,*) "Error in broadcasting new shift"
            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
        ENDIF

        AccRat=(REAL(Acceptances,r2))/(REAL(SumWalkersCyc,r2))      !The acceptance ratio which is printed is only for the current node - not summed over all nodes

        IF(iProcIndex.eq.Root) THEN
!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
            WRITE(15,"(I12,G16.7,I9,G16.7,I12,I8,G16.7,2I10,F13.5,2G13.5,2I6)") Iter,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,AllAnnihilated,     &
                    ProjectionE,AllNoatHF,AllNoatDoubs,AllPosFrac,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
            WRITE(6,"(I12,G16.7,I9,G16.7,I12,I8,G16.7,2I10,F13.5,2G13.5,2I6)") Iter,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,AllTotWalkers,AllAnnihilated,      &
                    ProjectionE,AllNoatHF,AllNoatDoubs,AllPosFrac,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
            CALL FLUSH(15)
            CALL FLUSH(6)
        ENDIF

!Now need to reinitialise all variables on all processers
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        MeanExcitLevel=0.D0
        SumWalkersCyc=0
        PosFrac=0.D0
        Annihilated=0
        Acceptances=0
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
        AllAnnihilated=0
        AllNoatHF=0
        AllNoatDoubs=0

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
                CurrentDets(:,Chosen)=CurrentDets(:,TotWalkers)
                CurrentSign(Chosen)=CurrentSign(TotWalkers)
                CurrentH(Chosen)=CurrentH(TotWalkers)
                CurrentIC(Chosen)=CurrentIC(TotWalkers)
                CALL CopyExitgenPar(CurrentExcits(TotWalkers),CurrentExcits(Chosen))
                CurrentExcits(TotWalkers)%ExitGenForDet=.false.

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
                CurrentDets(:,VecSlot)=CurrentDets(:,i)
                CurrentSign(VecSlot)=CurrentSign(i)
                CurrentH(VecSlot)=CurrentH(i)
                CurrentIC(VecSlot)=CurrentIC(i)
                CALL CopyExitgenPar(CurrentExcits(i),CurrentExcits(VecSlot))

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
!        use CalcData , only : i_P
!        use SystemData , only : Beta
!        use IntegralsData , only : nTay
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
        use CalcData, only : EXCITFUNCS
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot,error,HFConn,MemoryAlloc
        TYPE(HElement) :: rh,TempHii
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors

        IF(HElementSize.gt.1) THEN
            CALL Stop_All("FCIMCPar","FciMCPar cannot function with complex orbitals.")
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
        HFHash=CreateHash(HFDet)

!Setup excitation generator for the HF determinant
        CALL SetupExitgenPar(HFDet,HFExcit)
        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=G_VMC_Seed-iProcIndex

!Calculate Hii
        TempHii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
        Hii=REAL(TempHii%v,r2)

!Initialise variables for calculation on each node
        ProjectionE=0.D0
        PosFrac=0.D0
        SumENum=0.D0
        SumNoatHF=0
        NoatHF=0
        NoatDoubs=0
        MeanExcitLevel=0.D0
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        Annihilated=0
        Acceptances=0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllNoatHF=0
        AllNoatDoubs=0
        AllSumNoatHF=0.D0
        AllGrowRate=0.D0
        AllMeanExcitLevel=0.D0
        AllSumWalkersCyc=0
        AllPosFrac=0.D0

        IF(TNoAnnihil) THEN
            WRITE(6,*) "No Annihilation to occur. Results are likely not to converge on right value. Proceed with caution. "
        ENDIF
        IF(ICILEVEL.ne.0) THEN
            CALL Stop_All("InitFCIMCCalcPar","Truncated FCIMC not yet available in parallel.")
        ENDIF
        IF(TReadPops) THEN
            CALL Stop_All("InitFCIMCCalcPar","POPFILE read facility not yet available in parallel.")
        ENDIF

        IF(TResumFciMC) THEN
            IF(NDets.gt.2) THEN
                IF(.not.EXCITFUNCS(10)) THEN
                    WRITE(6,*) "Cannot have an excitation bias with multiple determinant graphs...exiting."
                    CALL Stop_All("InitFCIMCCalcPar","Cannot have biasing with Graphsizes > 2")
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
                CALL Stop_All("InitFCIMCCalcPar","Graphs cannot be smaller than two vertices")
            ENDIF
            IF(iProcIndex.eq.root) THEN
                WRITE(6,*) "Resumming in multiple transitions to/from each excitation"
                WRITE(6,"(A,I5,A)") "Graphs to resum will consist of ",NDets," determinants."
            ENDIF
            IF(.not.TNoAnnihil) THEN
                CALL Stop_All("InitFCIMCCalcPar","Annihilation is not currently compatable with ResumFCIMC in parallel")
            ENDIF
        ENDIF
        WRITE(6,*) ""
        WRITE(6,*) "Performing FCIMC...."
        IF(TStartSinglePart) THEN
            WRITE(6,"(A,F9.3,A,I9)") "Initial number of particles set to 1, and shift will be held at ",DiagSft," until particle number on root node gets to ",InitWalkers
        ELSE
            WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
        ENDIF
        WRITE(6,*) "Maximum connectivity of HF determinant is: ",HFConn
        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
        WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
        
        IF(TStartMP1) THEN
!Start the initial distribution off at the distribution of the MP1 eigenvector

            call Stop_All("InitFCIMCCalcPar","This feature not ready yet")
            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."
!            CALL StartWavevectorPar(1)

        ELSE
!initialise the particle positions - start at HF with positive sign
            IF(TStartSinglePart) THEN
                TSinglePartPhase=.true.
                IF(TReadPops) THEN
                    CALL Stop_All("InitFciMCCalcPar","Cannot read in POPSFILE as well as starting with a single particle")
                ENDIF
            ELSE
                TSinglePartPhase=.false.
            ENDIF

!Set the maximum number of walkers allowed
            MaxWalkers=MemoryFac*InitWalkers
            WRITE(6,*) "Memory allocated for a maximum particle number per node of: ",MaxWalkers

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("InitFciMCCalc","MPI error code ne 0")
            ENDIF

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
            ALLOCATE(WalkVecH(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkers,8,this_routine,WalkVecHTag,ierr)
            CALL AZZERO(WalkVecH,MaxWalkers)
            ALLOCATE(WalkVec2H(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2H',MaxWalkers,8,this_routine,WalkVec2HTag,ierr)
            CALL AZZERO(WalkVec2H,MaxWalkers)
            
            MemoryAlloc=((2*NEl+8)*MaxWalkers)*4    !Memory Allocated in bytes

            IF(.not.TNoAnnihil) THEN
                ALLOCATE(HashArray(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('HashArray',MaxWalkers,8,this_routine,HashArrayTag,ierr)
                HashArray(:)=0
                ALLOCATE(Hash2Array(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('Hash2Array',MaxWalkers,8,this_routine,Hash2ArrayTag,ierr)
                Hash2Array(:)=0
                ALLOCATE(IndexTable(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('IndexTable',MaxWalkers,4,this_routine,IndexTableTag,ierr)
                CALL IAZZERO(IndexTable,MaxWalkers)
                ALLOCATE(Index2Table(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('Index2Table',MaxWalkers,4,this_routine,Index2TableTag,ierr)
                CALL IAZZERO(Index2Table,MaxWalkers)
                ALLOCATE(ProcessVec(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('ProcessVec',MaxWalkers,4,this_routine,ProcessVecTag,ierr)
                CALL IAZZERO(ProcessVec,MaxWalkers)
                ALLOCATE(Process2Vec(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('Process2Vec',MaxWalkers,4,this_routine,Process2VecTag,ierr)
                CALL IAZZERO(Process2Vec,MaxWalkers)

                MemoryAlloc=MemoryAlloc+32*MaxWalkers
            ENDIF

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
                CurrentH(1)=0.D0
                IF(.not.TNoAnnihil) THEN
                    HashArray(1)=HFHash
                ENDIF
            ELSE

                do j=1,InitWalkers
                    CurrentDets(:,j)=HFDet(:)
                    CurrentSign(j)=.true.
                    CurrentIC(j)=0
                    CurrentH(j)=0.D0
                    IF(.not.TNoAnnihil) THEN
                        HashArray(j)=HFHash
                    ENDIF
                enddo
            ENDIF

            WRITE(6,"(A,F14.6,A)") "Initial memory (without excitgens + temp arrays) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial memory allocation sucessful..."
            CALL FLUSH(6)

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("InitFciMCCalc","MPI error code ne 0")
            ENDIF

            ALLOCATE(WalkVecExcits(MaxWalkers),stat=ierr)
            ALLOCATE(WalkVec2Excits(MaxWalkers),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("InitFCIMMCCalcPar","Error in allocating walker excitation generators")

!Allocate pointers to the correct excitation arrays
            CurrentExcits=>WalkVecExcits
            NewExcits=>WalkVec2Excits

            IF(TStartSinglePart) THEN
                CALL CopyExitGenPar(HFExcit,CurrentExcits(1))
            ELSE
                do j=1,InitWalkers
!Copy the HF excitation generator accross to each initial particle
                    CALL CopyExitGenPar(HFExcit,CurrentExcits(j))
                enddo
            ENDIF
            MemoryAlloc=((HFExcit%nExcitMemLen)+2)*4*MaxWalkers

            WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial allocation of excitation generators successful..."
            WRITE(6,*) "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkers*12,r2)/1048576.D0," Mb/Processor"
            CALL FLUSH(6)

        ENDIF

        IF(TStartSinglePart) THEN
            TotWalkers=1
            TotWalkersOld=1
!Initialise global variables for calculation on the root node
            IF(iProcIndex.eq.root) THEN
                AllTotWalkers=nProcessors
                AllTotWalkersOld=nProcessors
            ENDIF
        ELSE
!TotWalkers contains the number of current walkers at each step
            TotWalkers=InitWalkers
            TotWalkersOld=InitWalkers
!Initialise global variables for calculation on the root node
            IF(iProcIndex.eq.root) THEN
                AllTotWalkers=InitWalkers*nProcessors
                AllTotWalkersOld=InitWalkers*nProcessors
            ENDIF

        ENDIF

        CALL IAZZERO(CullInfo,30)
        NoCulls=0

        IF(.not.TNoAnnihil) THEN
!Need to declare a new MPI type to deal with the long integers we use in the hashing
    
            CALL MPI_Type_create_f90_integer(18,mpilongintegertype,error)
            CALL MPI_Type_commit(mpilongintegertype,error)
        
        ENDIF

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        RETURN

    END SUBROUTINE InitFCIMCCalcPar

!This routine flips the sign of all particles on the node
    SUBROUTINE FlipSign()
        INTEGER :: i

        do i=1,TotWalkers
            CurrentSign(i)=.not.CurrentSign(i)
        enddo
        RETURN
    
    END SUBROUTINE FlipSign

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
    
    FUNCTION CreateHash(DetCurr)
        INTEGER :: DetCurr(NEl),i
        INTEGER(KIND=i2) :: CreateHash

        CreateHash=0
        do i=1,NEl
            CreateHash=13*CreateHash+i*DetCurr(i)
!            CreateHash=(1099511628211*CreateHash)!+i*DetCurr(i)
            
!            CreateHash=mod(1099511628211*CreateHash,2**64)
!            CreateHash=XOR(CreateHash,DetCurr(i))
        enddo
!        WRITE(6,*) CreateHash
        RETURN

    END FUNCTION CreateHash

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
!            IF(OrigExit%nExcitMemLen.eq.0) THEN
!                CALL Stop_All("CopyExitgenPar","Problem allocating memory for new excit")
!            ENDIF
            IF(ierr.ne.0) CALL Stop_All("CopyExitgenPar","Problem with allocating memory for new excitation generator")
            NewExit%ExcitData(:)=OrigExit%ExcitData(:)
            NewExit%nExcitMemLen=OrigExit%nExcitMemLen
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
                CALL Stop_All("SetupExitgenPar","Excitation generator meant to already be set up")
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

!Indicate that the excitation generator is now correctly allocated.
            ExcitGen%ExitGenForDet=.true.
        
        ENDIF

    END SUBROUTINE SetupExitgenPar

END MODULE FciMCParMod

#else

MODULE FciMCParMod
!Dummy module so we can use it in serial - contains all global variables
    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,Arr,nMsh
    use CalcData , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    use CalcData , only : TStartMP1
    use CalcData , only : GrowMaxFactor,CullFactor
    use CalcData , only : RhoApp,TResumFCIMC
    USE Determinants , only : FDet,GetHElement2
    USE DetCalc , only : NMRKS
    use IntegralsData , only : fck,NMax,UMat
    USE global_utilities
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
    REAL*8 , ALLOCATABLE , TARGET :: WalkVecH(:),WalkVec2H(:)       !First element is diagonal hamiltonian element - second is the connection to HF determinant
    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0
    INTEGER :: WalkVecICTag=0,WalkVec2ICTag=0,WalkVecHTag=0,WalkVec2HTag=0

!Pointers to point at the correct arrays for use
    INTEGER , POINTER :: CurrentDets(:,:), NewDets(:,:)
    LOGICAL , POINTER :: CurrentSign(:), NewSign(:)
    INTEGER , POINTER :: CurrentIC(:), NewIC(:)
    REAL*8 , POINTER :: CurrentH(:), NewH(:)
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

    contains

    SUBROUTINE FciMCPar(Weight,Energyxw)
    TYPE(HDElement) :: Weight,Energyxw

        CALL Stop_All("FciMCPar","Entering the wrong FCIMCPar parallel routine")

    END SUBROUTINE FciMCPar

END MODULE FciMCParMod
    
#endif
